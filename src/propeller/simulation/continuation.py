"""Continue an archived serial solve in a new, independently recorded attempt."""
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import time

from ..config import validate_config
from ..experiments import Store, atomic_json, canonical, evaluation_key
from .parser import summarize
from .runner import run_command, runtime_info
from .field_io import normalize_ascii_fields


def archived_inputs(parent):
    parent = Path(parent).resolve()
    result = json.loads((parent / "result.json").read_text())
    manifest = json.loads((parent / "manifest.json").read_text())
    original = manifest["config"]
    if result["status"] not in ("SUCCESS", "NOT_CONVERGED"):
        raise ValueError("Continuation requires a completed solver evaluation")
    if original["cfd"].get("mpi_ranks", 1) != 1:
        raise ValueError("Continuation currently supports serial cases")
    case = parent / "case"
    if json.loads((case / "configuration.json").read_text()) != original:
        raise ValueError("Parent case configuration differs from its manifest")
    times = [(float(p.name), p) for p in case.iterdir()
             if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)]
    latest, fields = max(times)
    if latest <= 0:
        raise ValueError("No saved solution fields are available")
    for field in ("U", "p", "k", "omega", "nut"):
        if not (fields / field).is_file():
            raise ValueError("Missing restart field: " + field)
    if not result.get("mesh", {}).get("valid"):
        raise ValueError("Parent mesh has no passing quality record")
    rotor = case / "constant/triSurface/rotor.stl"
    if hashlib.sha256(rotor.read_bytes()).hexdigest() != result["geometry_sha256"]:
        raise ValueError("Parent geometry checksum differs")
    inputs = []
    for folder in ("constant", "system", "0", fields.name):
        inputs.extend(p for p in (case / folder).rglob("*") if p.is_file())
    hashes = {str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(inputs)}
    return parent, result, manifest, fields.name, hashes


def restart_inputs(parent, config):
    inputs = archived_inputs(parent)
    original = inputs[2]["config"]
    validate_config(config)
    expected = copy.deepcopy(original)
    for key in ("iterations", "timeout_seconds"):
        expected["cfd"][key] = config["cfd"][key]
    if config['cfd']['force_window']<original['cfd']['force_window']:
        raise ValueError('Continuation cannot shorten the averaging window')
    expected['cfd']['force_window']=config['cfd']['force_window']
    if expected != config:
        raise ValueError("Continuation may change only iterations, timeout_seconds and a larger force_window")
    if config["cfd"]["iterations"] - float(inputs[3]) < config["cfd"]["force_window"]:
        raise ValueError("Continuation must add at least one complete force window")
    return inputs


def continue_evaluation(parent, config, root="runs/evaluations"):
    parent, previous, original, latest, hashes = restart_inputs(parent, config)
    runtime = runtime_info(config)
    if not runtime["available"]:
        raise ValueError("Native OpenFOAM runtime is unavailable")
    if runtime.get("binary_sha256") != original["runtime"].get("binary_sha256"):
        raise ValueError("Continuation runtime differs from the parent calculation")
    key, manifest = evaluation_key(config)
    provenance = {"parent": str(parent), "restart_time": latest,
                  "parent_result_sha256": hashlib.sha256((parent / "result.json").read_bytes()).hexdigest(),
                  "input_sha256": hashes}
    manifest.update(runtime=runtime, continuation=provenance)
    key = hashlib.sha256((key + canonical(runtime) + canonical(provenance)).encode()).hexdigest()
    store = Store(root)
    identity, directory = store.start(key, manifest)
    started = time.monotonic()
    commands = []
    try:
        case = directory / "case"
        case.mkdir()
        for folder in ("constant", "system", "0", latest):
            shutil.copytree(parent / "case" / folder, case / folder)
        provenance['field_normalization'] = normalize_ascii_fields(case / latest)
        for file in ("geometry.json", "case.foam", "log.checkMesh", "log.snappyHexMesh"):
            shutil.copyfile(parent / "case" / file, case / file)
        for file in ("validation.json", "preview.html", "mesh.json"):
            if (parent / file).exists():
                shutil.copyfile(parent / file, directory / file)
        control = case / "system/controlDict"
        text = control.read_text()
        text, nstart = re.subn(r"\bstartFrom\s+\w+\s*;", "startFrom latestTime;", text)
        text, nend = re.subn(r"\bendTime\s+[^;]+;", "endTime %s;" % config["cfd"]["iterations"], text)
        if nstart != 1 or nend != 1:
            raise ValueError("Ambiguous restart control dictionary")
        # Keep restart fields readable and ensure the final requested state is
        # written even when the old write interval does not divide the new end.
        end = config['cfd']['iterations']
        interval = math.gcd(end, max(50, end//5))
        if interval < 50:interval = end
        text, count = re.subn(r'\bwriteFormat\s+\w+\s*;', 'writeFormat binary;', text)
        if not count:text += '\nwriteFormat binary;\n'
        text, count = re.subn(r'\bwriteInterval\s+\d+\s*;', 'writeInterval %s;' % interval, text, count=1)
        if not count:text += '\nwriteInterval %s;\n' % interval
        control.write_text(text)
        atomic_json(case / "configuration.json", config)
        atomic_json(directory / "restart.json", provenance)
        command = run_command(config, case, ["simpleFoam"])
        commands.append(command)
        atomic_json(directory / "commands.json", commands)
        if command["status"] != "SUCCESS":
            result = {"status": "TIMEOUT" if command["status"] == "TIMEOUT" else "SOLVER_FAILED", "reasons": [command]}
        else:
            result = summarize(case, config, previous["geometry"]["diameter"])
        result.update(geometry=previous["geometry"], mesh=previous["mesh"],
                      geometry_sha256=previous["geometry_sha256"])
    except KeyboardInterrupt:
        try:
            store.finish(identity, directory, {"status": "INTERRUPTED", "reasons": ["User interruption"], "qualified": False})
        finally:
            store.close()
        raise
    except Exception as exc:
        result = {"status": "INTERNAL_ERROR", "reasons": [type(exc).__name__ + ": " + str(exc)]}
    finally:
        if "result" not in locals():
            store.close()
    result.update(qualified=False, commands=commands, seconds=time.monotonic()-started,
                  cache_hit=False, cache_key=key, continuation=provenance)
    try:
        return store.finish(identity, directory, result)
    finally:
        store.close()

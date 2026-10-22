"""Map a converged solution onto a sensitivity mesh as an initial guess."""
import hashlib
import json
from pathlib import Path
import shutil

from .continuation import archived_inputs
from .case_builder import header
from .field_io import normalize_ascii_fields
from ..experiments import atomic_json


def initialization_inputs(parent, config, runtime):
    parent, result, manifest, latest, hashes = archived_inputs(parent)
    if result["status"] != "SUCCESS":
        raise ValueError("Mesh sensitivity initialization requires a converged source")
    original = manifest["config"]
    for group in ("design", "parameters", "geometry"):
        if original[group] != config[group]:
            raise ValueError("Sensitivity initialization requires identical " + group)
    def physics(c):
        return {k: v for k, v in c["operating"].items()
                if k not in ("target_thrust", "reference_power")}
    if physics(original) != physics(config):
        raise ValueError("Sensitivity initialization requires identical operating conditions")
    if runtime.get("binary_sha256") != manifest["runtime"].get("binary_sha256"):
        raise ValueError("Sensitivity initialization runtime differs")
    if config["cfd"].get("mpi_ranks", 1) != 1:
        raise ValueError("Sensitivity initialization currently supports serial cases")
    return {"parent": str(parent), "source_time": latest,
            "parent_result_sha256": hashlib.sha256((parent / "result.json").read_bytes()).hexdigest(),
            "input_sha256": hashes,
            "consistent_domain": original["cfd"].get("domain") == config["cfd"].get("domain"),
            "method": "interpolate", "purpose": "initial guess only; target must independently converge"}


def mapping_command(case, initialization):
    # Different domains preserve target farfield conditions and interpolate
    # overlapping interior cells. Uncovered cells keep their configured guess.
    source = Path(initialization["parent"]) / "case"
    for relative, expected in initialization.get("input_sha256", {}).items():
        if hashlib.sha256((source / relative).read_bytes()).hexdigest() != expected:
            raise ValueError("Initialization source changed: " + relative)
    staged = case.parent / "initialization-source"
    record = staged / "normalization.json"
    if not staged.exists():
        staged.mkdir()
        for folder in ("constant", "system", "0", initialization["source_time"]):
            shutil.copytree(source / folder, staged / folder)
        normalization = normalize_ascii_fields(staged / initialization["source_time"])
        atomic_json(record,normalization)
    else:
        normalization = json.loads(record.read_text())
    initialization['staged_source'] = str(staged.resolve())
    initialization['field_normalization'] = normalization
    (case / "system/mapFieldsDict").write_text(
        header("mapFieldsDict") + "patchMap (rotor rotor); cuttingPatches (farfield);\n")
    command = ["mapFields", str(staged.resolve()),
               "-sourceTime", initialization["source_time"], "-mapMethod", "interpolate"]
    if initialization["consistent_domain"]:
        command.append("-consistent")
    return command


def check_geometry(geometry, initialization):
    """Guard against geometry code changes hidden behind identical parameters."""
    import trimesh
    from scipy.spatial import cKDTree
    source = trimesh.load(Path(initialization["parent"]) / "case/constant/triSurface/rotor.stl")
    candidate = geometry.mesh
    distance = float(max(cKDTree(source.vertices).query(candidate.vertices)[0].max(),
                         cKDTree(candidate.vertices).query(source.vertices)[0].max()))
    volume_error = abs(source.volume-candidate.volume)/max(abs(source.volume),1e-30)
    area_error = abs(source.area-candidate.area)/max(abs(source.area),1e-30)
    # Exported STL vertices are float32. This tolerance is below mesh resolution
    # and protects against changed geometry while allowing storage roundoff.
    tolerance = 1e-7 * geometry.measurements['diameter']
    if distance > tolerance or max(volume_error,area_error)>1e-6:
        raise ValueError("Mapped sensitivity geometry differs from the source rotor")
    return {"maximum_nearest_vertex_distance_m":distance,"tolerance_m":tolerance,
            "relative_volume_error":float(volume_error),"relative_area_error":float(area_error)}

"""Geometry-first CFD evaluator; unavailable and failed runs never have a score."""
import hashlib
from pathlib import Path
import time
from ..experiments import Store,evaluation_key,atomic_json
from ..geometry import generate_rotor,validate_geometry
from ..geometry.preview import write_preview
from .case_builder import build_case
from .runner import runtime_info,run_command
from .parser import summarize,mesh_diagnostics


def evaluate_cfd(config,root="runs/evaluations",use_cache=True):
    store=Store(root); key,manifest=evaluation_key(config)
    # Verify actual runtime before retrieving a stored calculation.
    runtime=runtime_info(config);manifest["runtime"]=runtime
    if not runtime["available"]:
        store.close()
        return {"status":"RUNTIME_UNAVAILABLE","reasons":[runtime],"qualified":False}
    key=hashlib.sha256((key+str(runtime)).encode()).hexdigest()
    if use_cache:
        cached=store.cached(key)
        if cached: store.close();return dict(cached,cache_hit=True)
    identity,directory=store.start(key,manifest);commands=[];start=time.monotonic();mesh_info=None
    def finish(result):
        result.update(commands=commands,seconds=time.monotonic()-start,cache_hit=False,cache_key=key)
        try: return store.finish(identity,directory,result)
        finally: store.close()
    try:
        try: geometry=generate_rotor(config)
        except Exception as exc: return finish({"status":"GEOMETRY_INVALID","reasons":[str(exc)]})
        validation=validate_geometry(geometry,config,deep=True)
        atomic_json(directory/"validation.json",validation.to_dict())
        write_preview(geometry,validation,config,directory/"preview.html")
        if not validation.valid: return finish({"status":"GEOMETRY_INVALID","reasons":validation.violations})
        try: case=build_case(geometry,config,directory/"case")
        except Exception as exc: return finish({"status":"EXPORT_FAILED","reasons":[str(exc)]})
        pipeline=[["blockMesh"],["surfaceFeatureExtract"],["snappyHexMesh","-overwrite"],["checkMesh","-allTopology","-meshQuality"],["topoSet"]]
        ranks=config["cfd"].get("mpi_ranks",1)
        if ranks>1:pipeline.extend([["decomposePar"],["mpirun","-np",str(ranks),"simpleFoam","-parallel"]])
        else:pipeline.append(["simpleFoam"])
        for command in pipeline:
            attempt=run_command(config,case,command);commands.append(attempt)
            atomic_json(directory/"commands.json",commands)
            if attempt["status"]!="SUCCESS":
                return finish({"status":"TIMEOUT" if attempt["status"]=="TIMEOUT" else ("SOLVER_FAILED" if command[0] in ("simpleFoam","mpirun") else "MESH_FAILED"),"reasons":[attempt]})
            if command[0]=="checkMesh":
                mesh_info=mesh_diagnostics(case,config)
                atomic_json(directory/"mesh.json",mesh_info)
                if not mesh_info["valid"]:return finish({"status":"MESH_FAILED","reasons":mesh_info["errors"]})
        try:result=summarize(case,config,geometry.measurements["diameter"])
        except ValueError as exc:result={"status":"NOT_CONVERGED","reasons":["Unusable force/convergence evidence: "+str(exc)],"qualified":False}
        result["geometry"]=geometry.measurements
        result["mesh"]=mesh_info
        result["geometry_sha256"]=hashlib.sha256((case/"constant/triSurface/rotor.stl").read_bytes()).hexdigest()
        # Numerical convergence does not imply experimental or mesh qualification.
        result["qualified"]=False
        return finish(result)
    except KeyboardInterrupt:
        finish({"status":"INTERRUPTED","reasons":["User interruption"]});raise
    except Exception as exc:
        return finish({"status":"INTERNAL_ERROR","reasons":[type(exc).__name__+": "+str(exc)]})

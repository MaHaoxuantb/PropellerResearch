"""Conservative comparison reports. No confidence claim from a lone CFD solve."""
import json
from pathlib import Path


def compare_results(paths):
    rows=[]
    for path in paths:
        result=json.loads(Path(path).read_text())
        rows.append({"file":str(path),"status":result.get("status"),"metrics":result.get("metrics"),"mesh":result.get("mesh"),"qualified":result.get("qualified",False)})
    if len(rows)<2:raise ValueError("At least two evaluations are required")
    valid=[r for r in rows if r["status"]=="SUCCESS" and r["metrics"]]
    changes=[]
    for before,after in zip(valid,valid[1:]):
        changes.append({"from":before["file"],"to":after["file"],"relative_change":{
            key:(after["metrics"][key]-before["metrics"][key])/abs(before["metrics"][key]) if before["metrics"][key] else None for key in ("thrust","shaft_power")}})
    return {"evaluations":rows,"changes":changes,"interpretation":"Descriptive differences only. Establish matching operating conditions, mesh/domain/time sensitivity and uncertainty before ranking designs."}

"""Parse OpenCFD force histories, solver residuals and convergence evidence."""
import math
from pathlib import Path
import re
import numpy as np

NUMBER=r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?"


def residual_history(log):
    records={};iteration=None
    for line in log.splitlines():
        match=re.match(r'Time = ('+NUMBER+')',line)
        if match:iteration=float(match.group(1))
        match=re.search(r'Solving for (\w+), Initial residual = ('+NUMBER+')',line)
        if match and iteration is not None:
            fields=records.setdefault(iteration,{})
            fields[match.group(1)]=max(fields.get(match.group(1),0),float(match.group(2)))
    return records


def read_vectors(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"): continue
        values=[float(v) for v in re.findall(NUMBER,line)]
        if len(values)!=10: raise ValueError("Expected time + total/pressure/viscous vectors in "+str(path))
        rows.append(values)
    a=np.asarray(rows,dtype=float)
    if a.ndim!=2 or not np.all(np.isfinite(a)): raise ValueError("Empty/non-finite force history")
    if np.any(np.diff(a[:,0])<=0): raise ValueError("Force times must be strictly increasing")
    return a


def summarize(case,config,diameter):
    case=Path(case); folders=list((case/"postProcessing/forces").glob("*"))
    folders=[p for p in folders if (p/"force.dat").exists() and (p/"moment.dat").exists()]
    if not folders: raise ValueError("Missing paired force.dat / moment.dat")
    directory=max(folders,key=lambda p:float(p.name))
    force,moment=read_vectors(directory/"force.dat"),read_vectors(directory/"moment.dat")
    if force.shape!=moment.shape or not np.array_equal(force[:,0],moment[:,0]): raise ValueError("Force/moment timestamps differ")
    c,o=config["cfd"],config["operating"];window=c["force_window"]
    omega=o["rotation_sign"]*o["rpm"]*math.pi/30
    thrust=force[:,3];power=-moment[:,3]*omega
    tailT,tailP=thrust[-window:],power[-window:]
    T,P=float(tailT.mean()),float(tailP.mean())
    def stability(values):
        scale=max(abs(float(values.mean())),1e-12)
        mid=len(values)//2
        return {"relative_std":float(np.std(values)/scale),"relative_drift":float(abs(values[mid:].mean()-values[:mid].mean())/scale) if mid else float("inf")}
    tstats,pstats=stability(tailT),stability(tailP)
    log=(case/"log.simpleFoam").read_text(errors="replace")
    times=set(force[-window:,0]);history=residual_history(log)
    history={t:fields for t,fields in history.items() if t in times}
    residuals={field:max(fields[field] for fields in history.values() if field in fields)
               for field in ('p','Ux','Uy','Uz','k','omega') if any(field in fields for fields in history.values())}
    continuity=re.findall(r"sum local = ("+NUMBER+r")",log)
    continuity=float(continuity[-1]) if continuity else None
    yplus=re.findall(r"patch rotor y\+ : min = ("+NUMBER+r"), max = ("+NUMBER+r"), average = ("+NUMBER+r")",log)
    wall=dict(zip(("min","max","average"),map(float,yplus[-1]))) if yplus else None
    reasons=[]
    if len(force)<window: reasons.append("Insufficient force samples")
    if not log.rstrip().endswith("End"): reasons.append("Solver did not finish cleanly")
    if len(history)!=len(times) or any(not {"p","Ux","Uy","Uz","k","omega"}.issubset(fields) for fields in history.values()): reasons.append("Missing residual evidence in averaging window")
    if any(v>c["residual_tolerance"] for v in residuals.values()): reasons.append("Residual threshold exceeded")
    if continuity is None or abs(continuity)>c["continuity_tolerance"]: reasons.append("Continuity threshold missing/exceeded")
    if wall is None:reasons.append("Missing rotor y+ evidence")
    if max(tstats.values())>c["force_relative_tolerance"] or max(pstats.values())>c["force_relative_tolerance"]: reasons.append("Thrust/power not stable")
    n=o["rpm"]/60;rho=o["density"];area=math.pi*(diameter/2)**2
    metrics={"thrust":T,"shaft_power":P,"fluid_torque_z":float(moment[-window:,3].mean()),
             "CT":T/(rho*n*n*diameter**4),"CP":P/(rho*n**3*diameter**5),"J":o["inflow_speed"]/(n*diameter),
             "efficiency":T*o["inflow_speed"]/P if P>0 and o["inflow_speed"]>0 and T>0 else None,
             "figure_of_merit":T**1.5/(P*math.sqrt(2*rho*area)) if T>0 and P>0 and o["inflow_speed"]==0 else None}
    return {"status":"NOT_CONVERGED" if reasons else "SUCCESS","metrics":metrics,"reasons":reasons,
            "diagnostics":{"samples":len(force),"last_time":float(force[-1,0]),"window":window,"thrust_stability":tstats,"power_stability":pstats,"residuals":residuals,"residual_window_samples":len(history),"local_continuity":continuity,"y_plus":wall},
            "qualified":False,"provenance":"native-openfoam-"+c["version"]}


def mesh_diagnostics(case,config):
    case=Path(case)
    log=(case/"log.checkMesh").read_text(errors="replace")
    meshlog=(case/"log.snappyHexMesh").read_text(errors="replace")
    def last(pattern):
        values=re.findall(pattern,log)
        return float(values[-1]) if values else None
    cells=last(r"\bcells:\s*(\d+)")
    nonorth=last(r"Mesh non-orthogonality Max:\s*("+NUMBER+r")")
    skew=last(r"Max skewness\s*=\s*("+NUMBER+r")")
    coverage=re.findall(r"Extruding \d+ out of \d+ faces \(("+NUMBER+r")%\)",meshlog)
    errors=[]
    if "Mesh OK." not in log:errors.append("checkMesh did not pass")
    if cells is None or cells>config['cfd']['max_cells']:errors.append("Cell count missing or over resource limit")
    if nonorth is None or nonorth>config['cfd']['max_non_orthogonality']:errors.append("Non-orthogonality limit missing/exceeded")
    if skew is None or skew>config['cfd']['max_skewness']:errors.append("Skewness limit missing/exceeded")
    return {"valid":not errors,"errors":errors,"cells":cells,"max_non_orthogonality":nonorth,"max_skewness":skew,
            "layer_extrusion_percent":float(coverage[-1]) if coverage else None,
            "wall_treatment":"SST with Spalding nut wall function; actual y+ must be inspected"}

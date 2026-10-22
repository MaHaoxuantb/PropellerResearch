"""Conservative comparison reports. No confidence claim from a lone CFD solve."""
import json
import math
from pathlib import Path
import numpy as np
from .simulation.parser import read_vectors,residual_history


def assess_baseline(path,window=500):
    """Check a longer force window before promoting a numerically accepted seed."""
    path=Path(path).resolve();result=json.loads(path.read_text())
    config=json.loads((path.parent/'manifest.json').read_text())['config']
    if not isinstance(window,int) or window<4:raise ValueError('Assessment window must be an integer >=4')
    case=path.parent/'case'
    folders=[p for p in (case/'postProcessing/forces').glob('*') if (p/'force.dat').exists() and (p/'moment.dat').exists()]
    if not folders:raise ValueError('Missing force histories')
    folder=max(folders,key=lambda p:float(p.name))
    force=read_vectors(folder/'force.dat');moment=read_vectors(folder/'moment.dat')
    if force.shape!=moment.shape or not np.array_equal(force[:,0],moment[:,0]):raise ValueError('Force/moment histories differ')
    omega=config['operating']['rotation_sign']*config['operating']['rpm']*math.pi/30
    reasons=[];statistics={}
    if result.get('status')!='SUCCESS':reasons.append('Numerical convergence has not been accepted')
    if len(force)<window:reasons.append('Insufficient samples for the longer assessment window')
    for name,values in [('thrust',force[-window:,3]),('shaft_power',-moment[-window:,3]*omega)]:
        mean=float(values.mean());scale=max(abs(mean),1e-12);middle=len(values)//2
        statistics[name]={'mean':mean,'relative_std':float(np.std(values)/scale),
                          'relative_drift':float(abs(values[middle:].mean()-values[:middle].mean())/scale) if middle else None}
        if mean<=0:reasons.append(name+' must be positive')
        if max(v for k,v in statistics[name].items() if k.startswith('relative_') and v is not None)>config['cfd']['force_relative_tolerance']:
            reasons.append(name+' fails the longer-window stability check')
    times=set(force[-window:,0])
    residuals={t:fields for t,fields in residual_history((case/'log.simpleFoam').read_text(errors='replace')).items() if t in times}
    fields_required={'p','Ux','Uy','Uz','k','omega'}
    if len(residuals)!=len(times) or any(not fields_required.issubset(fields) for fields in residuals.values()):
        reasons.append('Incomplete residual evidence in the longer assessment window')
    maxima={field:max(fields[field] for fields in residuals.values() if field in fields)
            for field in fields_required if any(field in fields for fields in residuals.values())}
    if any(v>config['cfd']['residual_tolerance'] for v in maxima.values()):
        reasons.append('Residuals have not stayed below threshold throughout the assessment window')
    return {'status':'BASELINE_REJECTED' if reasons else 'READY_FOR_SENSITIVITY','source':str(path),
            'window':window,'available_samples':len(force),'statistics':statistics,'reasons':reasons,
            'residual_window_maximum':maxima,'residual_window_samples':len(residuals),
            'qualified':False,'interpretation':'A longer-window numerical settling check; mesh/domain and physical qualification remain separate.'}


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

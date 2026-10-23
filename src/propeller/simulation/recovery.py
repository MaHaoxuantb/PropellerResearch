"""Audit completed orphaned cases without fabricating process exit codes."""
import hashlib
import json
from pathlib import Path
import re
import time
from .runner import process_alive
from .parser import summarize, mesh_diagnostics
from ..experiments import Store


def live_case_process(case):
    case=Path(case)
    marker=case/'active-process.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        if record.get('status')=='RUNNING' and process_alive(record.get('pid')):
            return True
    # Legacy cases did not record subprocess ownership. A live log PID is
    # conservatively treated as busy, including permission-denied checks.
    for log in case.glob('log.*'):
        head=log.read_bytes()[:2500].decode(errors='replace')
        match=re.search(r'^PID\s*:\s*(\d+)',head,re.M)
        if match and process_alive(int(match.group(1))):return True
    return False


def recover_completed(directory):
    directory=Path(directory).resolve()
    if (directory/'result.json').exists():return json.loads((directory/'result.json').read_text())
    state=json.loads((directory/'state.json').read_text())
    if process_alive(state.get('pid')) or live_case_process(directory/'case'):
        raise ValueError('Evaluation or native process is still alive; cannot recover')
    manifest=json.loads((directory/'manifest.json').read_text());config=manifest['config'];case=directory/'case'
    if json.loads((case/'configuration.json').read_text())!=config:raise ValueError('Archived configuration mismatch')
    validation=json.loads((directory/'validation.json').read_text())
    if not validation['valid'] or validation['checks'].get('self_intersection_scan')!='completed':
        raise ValueError('Missing full passing geometry validation')
    mesh=mesh_diagnostics(case,config)
    if not mesh['valid']:raise ValueError('Mesh is not accepted')
    log=(case/'log.simpleFoam').read_text(errors='replace')
    times=re.findall(r'^Time = (\S+)',log,re.M)
    if not log.rstrip().endswith('End') or 'FOAM FATAL' in log or not times or float(times[-1])!=config['cfd']['iterations']:
        raise ValueError('Solver has no clean completion at its requested end')
    geometry=json.loads((case/'geometry.json').read_text())
    result=summarize(case,config,geometry['diameter'])
    if result['diagnostics']['last_time']!=config['cfd']['iterations']:raise ValueError('Final force sample missing')
    evidence={str(p.relative_to(directory)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [directory/'manifest.json',directory/'validation.json',case/'log.simpleFoam',case/'log.checkMesh',case/'constant/triSurface/rotor.stl']}
    result.update(geometry=geometry,mesh=mesh,geometry_sha256=evidence['case/constant/triSurface/rotor.stl'],
                  recovered=True,cache_hit=False,qualified=False,
                  recovery={'time':time.time(),'evidence_sha256':evidence,'observed_exit_code':None,
                            'basis':'Clean complete native logs and re-parsed geometry/mesh/convergence evidence; original supervisor exited.'})
    if 'continuation' in manifest:result['continuation']=manifest['continuation']
    if 'initialization' in manifest:result['initialization']=manifest['initialization']
    commands=directory/'commands.json'
    result['commands']=json.loads(commands.read_text()) if commands.exists() else []
    result['commands'].append({'status':'RECOVERED_COMPLETION','returncode':None,'log':str(case/'log.simpleFoam')})
    store=Store(directory.parent.parent)
    try:return store.finish(directory.name,directory,result)
    finally:store.close()

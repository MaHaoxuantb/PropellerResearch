"""Attach an objective to accepted CFD data without rerunning its baseline."""
import copy
import hashlib
import json
import math
from pathlib import Path

from .analysis import assess_baseline
from .experiments import atomic_json,canonical
from .geometry import generate_rotor
from .simulation.initialization import check_geometry
from .simulation.parser import summarize
from .simulation.runner import runtime_info


def simulation_contract(config):
    value={name:copy.deepcopy(config[name]) for name in ('design','parameters','geometry','operating','cfd')}
    for key in ('target_thrust','reference_power'):value['operating'].pop(key)
    for key in ('timeout_seconds','startup_timeout_seconds','force_window'):value['cfd'].pop(key,None)
    return value


def _digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze_baseline(config,result_path,output,study_path=None,calibrate_timeout=False):
    result_path=Path(result_path).resolve();output=Path(output).resolve()
    result=json.loads(result_path.read_text());manifest=json.loads((result_path.parent/'manifest.json').read_text())
    if result.get('status')!='SUCCESS':raise ValueError('Baseline must have SUCCESS convergence status')
    if config!=manifest['config']:raise ValueError('Selected configuration does not match the evaluated baseline')
    assessment=assess_baseline(result_path,config['cfd']['force_window'])
    if assessment['status']!='READY_FOR_SENSITIVITY':raise ValueError('Baseline fails sustained convergence: '+str(assessment['reasons']))
    parsed=summarize(result_path.parent/'case',config,result['geometry']['diameter'])
    T,P=parsed['metrics']['thrust'],parsed['metrics']['shaft_power']
    if parsed['status']!='SUCCESS' or T<=0 or P<=0:raise ValueError('Baseline must have accepted positive thrust and shaft power')
    record={'schema_version':1,'source':str(result_path),'source_result_sha256':_digest(result_path),
            'simulation_contract_sha256':hashlib.sha256(canonical(simulation_contract(config)).encode()).hexdigest(),
            'target_thrust':T,'reference_power':P,'window':config['cfd']['force_window'],'qualified':False}
    if study_path:
        study_path=Path(study_path).resolve();study=json.loads(study_path.read_text())
        if study['status']!='READY_FOR_EXPLORATORY_SEARCH' or study.get('recommended_search_source')!=str(result_path):
            raise ValueError('Study does not select this result as its pilot search fidelity')
        record['study']={'path':str(study_path),'sha256':_digest(study_path)}
    reference=Path(str(output)+'.baseline.json')
    if output.exists() or reference.exists():raise FileExistsError('Baseline configuration/reference already exists')
    frozen=copy.deepcopy(config);frozen['operating'].update(target_thrust=T,reference_power=P)
    solver=[c for c in result.get('commands',[]) if Path(c.get('log','')).name=='log.simpleFoam']
    samples=result.get('diagnostics',{}).get('samples',0)
    if solver and samples>0:
        estimate=solver[-1]['seconds']/samples*config['cfd']['iterations']
        timeout=max(config['cfd']['timeout_seconds'],math.ceil(1.5*estimate))
        if calibrate_timeout:frozen['cfd']['timeout_seconds']=timeout
        record['runtime_calibration']={'observed_solver_seconds':solver[-1]['seconds'],'observed_iterations':samples,
                                       'estimated_full_solve_seconds':estimate,'recommended_timeout_seconds':timeout,
                                       'applied_timeout_seconds':frozen['cfd']['timeout_seconds'],'safety_factor':1.5}
    atomic_json(reference,record);atomic_json(output,frozen)
    return {'configuration':str(output),'reference':str(reference),'target_thrust':T,'reference_power':P,
            'timeout_seconds':frozen['cfd']['timeout_seconds'],'runtime_calibration':record.get('runtime_calibration'),'qualified':False}


def load_baseline_seed(config,reference):
    record=json.loads(Path(reference).read_text())
    if record.get('schema_version')!=1:raise ValueError('Legacy baseline reference; re-freeze it or use --fresh-baseline')
    path=Path(record['source']).resolve()
    if _digest(path)!=record['source_result_sha256']:raise ValueError('Baseline result checksum differs')
    result=json.loads(path.read_text());manifest=json.loads((path.parent/'manifest.json').read_text())
    original=simulation_contract(manifest['config'])
    if hashlib.sha256(canonical(original).encode()).hexdigest()!=record['simulation_contract_sha256']:
        raise ValueError('Baseline manifest contract checksum differs')
    if simulation_contract(config)!=original:
        raise ValueError('Baseline simulation contract differs from the requested search')
    if config['cfd']['force_window']<record['window']:raise ValueError('Cannot weaken baseline averaging-window acceptance')
    if 'study' in record and _digest(record['study']['path'])!=record['study']['sha256']:
        raise ValueError('Baseline study checksum differs')
    runtime=runtime_info(config)
    if not runtime['available'] or runtime.get('binary_sha256')!=manifest['runtime'].get('binary_sha256'):
        raise ValueError('Baseline native runtime differs or is unavailable')
    # Continuations inherit their mesh; use its actual generation ancestor.
    origin=path.parent;seen=set()
    while True:
        if origin in seen or len(seen)>=100:raise ValueError('Invalid baseline mesh history')
        seen.add(origin);ancestor=json.loads((origin/'manifest.json').read_text())
        parent=ancestor.get('continuation',{}).get('parent')
        if not parent:break
        origin=Path(parent).resolve()
    template=origin/'source/propeller/simulation/case_builder.py'
    if not template.exists() or _digest(template)!=_digest(Path(__file__).parent/'simulation/case_builder.py'):
        raise ValueError('Baseline mesh-generation template differs or is missing; use --fresh-baseline')
    validation=json.loads((path.parent/'validation.json').read_text())
    if not validation.get('valid') or validation.get('checks',{}).get('self_intersection_scan')!='completed' or not result.get('mesh',{}).get('valid'):
        raise ValueError('Baseline geometry/mesh acceptance evidence is incomplete')
    equivalence=check_geometry(generate_rotor(config),{'parent':str(path.parent)})
    assessed=assess_baseline(path,config['cfd']['force_window'])
    parsed=summarize(path.parent/'case',config,result['geometry']['diameter'])
    if parsed['status']!='SUCCESS' or assessed['status']!='READY_FOR_SENSITIVITY':
        raise ValueError('Baseline does not pass current sustained convergence policy')
    return dict(parsed,cache_hit=True,source_kind='accepted-baseline-reuse',path=str(path.parent),
                attempt_id=result.get('attempt_id'),baseline_source=str(path),baseline_source_sha256=record['source_result_sha256'],
                geometry_equivalence=equivalence,mesh_origin=str(origin))

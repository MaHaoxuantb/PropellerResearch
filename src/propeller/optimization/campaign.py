"""Durable, bounded pilot campaign. Native evaluations run in owned workers."""
import argparse
import copy
import fcntl
import hashlib
import html
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import shutil
import numpy as np
from ..config import load_config
from ..experiments import atomic_json,canonical,source_fingerprint
from ..geometry.screening import free_coordinates
from ..simulation.runner import process_alive,runtime_info
from .bayesian import accepted_metrics,candidate_pool,propose,pareto_indices,hypervolume


def read(path):return json.loads(Path(path).read_text())


def outputs(observations,baseline):
    t,p=baseline
    return [[o['metrics'][0]/t,(o['metrics'][0]/o['metrics'][1])/(t/p)] for o in observations]


def write_status(root,state):
    state['updated']=time.time();atomic_json(root/'state.json',state)
    obs=state['observations'];y=outputs(obs,state['baseline_metrics']);front=pareto_indices(y)
    atomic_json(root/'pareto.json',{'qualified':False,'objectives':['thrust_N','thrust_per_watt_N_W'],
        'radius':state['settings']['radius'],'frontier':[obs[i] for i in front],
        'all_accepted':obs,'status':state['status'],'stop_reason':state.get('stop_reason')})
    rows=[]
    for i,o in enumerate(obs):
        t,p=o['metrics'];rows.append('<tr><td>%s</td><td>%.6f</td><td>%.6f</td><td>%.6f</td><td>%s</td></tr>'%(html.escape(o['id']),t,p,t/p,'frontier' if i in front else 'dominated'))
    jobrows=''.join('<tr><td>%s</td><td>%s</td><td>%s</td></tr>'%(j['id'],html.escape(j.get('method','')),html.escape(j['status'])) for j in state['jobs'])
    # Plot only measured accepted points; no invented interpolated frontier.
    xmax=max(o['metrics'][0] for o in obs)*1.15;ymax=max(o['metrics'][0]/o['metrics'][1] for o in obs)*1.15
    dots=''.join('<circle cx="%.2f" cy="%.2f" r="6" fill="%s"><title>%s: %.5f N, %.5f N/W</title></circle>'%(60+o['metrics'][0]/xmax*600,320-(o['metrics'][0]/o['metrics'][1])/ymax*270,'#007f70' if i in front else '#999',html.escape(o['id']),o['metrics'][0],o['metrics'][0]/o['metrics'][1]) for i,o in enumerate(obs) if o['metrics'][0]>0)
    content='''<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="30"><title>Propeller optimization pilot</title>
<style>body{font:16px system-ui;max-width:960px;margin:40px auto;padding:0 20px;color:#20303a}td,th{padding:8px 16px;border-bottom:1px solid #ddd;text-align:left}svg{width:100%%;max-width:720px;background:#f6f8fa}code{overflow-wrap:anywhere}</style>
<h1>Fixed-radius propeller optimization</h1><p>Status: <b>%s</b>. %s</p><p>Native OpenFOAM + Bayesian optimization. Numerically accepted points only; mesh/domain and physical qualification remain pending.</p>
<p>The plot shows positive-thrust points; the table also retains converged reverse-thrust outcomes for learning.</p>
<p>Radius %.3f mm · %d new CFD evaluations completed · %d geometry rejections · maximum %.1f hours</p>
<svg viewBox="0 0 720 380" role="img" aria-label="Thrust per watt versus thrust"><path d="M60 45V320H680" fill="none" stroke="#333"/>%s<text x="260" y="365">Thrust (N), axis 0–%.3f</text><text x="65" y="25">Thrust per watt (N/W), axis 0–%.4f</text></svg>
<table><tr><th>Design</th><th>Thrust N</th><th>Power W</th><th>N/W</th><th>Comparison</th></tr>%s</table>
<h2>Jobs</h2><table><tr><th>ID</th><th>Proposal</th><th>Status</th></tr>%s</table>
<p>Stop gracefully by creating <code>STOP</code> in this campaign directory, or use the campaign stop command.</p>'''%(html.escape(state['status']),html.escape(state.get('stop_reason','')),state['settings']['radius']*1000,state['evaluations'],state['geometry_rejections'],state['settings']['max_hours'],dots,xmax,ymax,''.join(rows),jobrows)
    (root/'report.html').write_text(content)


def baseline_record(base,path,radius):
    from ..simulation.recovery import recover_completed
    from ..simulation.parser import summarize,mesh_diagnostics
    from ..simulation.initialization import check_geometry
    from ..geometry import generate_rotor,validate_geometry
    path=Path(path).resolve();r=recover_completed(path)
    source=read(path/'manifest.json');expected=copy.deepcopy(source['config'])
    # The optimizer settings and runtime timeout do not alter CFD physics.
    expected['optimization']=base['optimization']
    expected['cfd']['timeout_seconds']=base['cfd']['timeout_seconds']
    if expected!=base:raise ValueError('Baseline and search physics/mesh configurations differ')
    if runtime_info(base).get('binary_sha256')!=source['runtime'].get('binary_sha256'):
        raise ValueError('Baseline native runtime differs')
    g=generate_rotor(base);check_geometry(g,{'parent':str(path)})
    validation=validate_geometry(g,base,deep=True)
    if not validation.valid:raise ValueError('Baseline fails current full geometry validation')
    if abs(g.measurements['diameter']/2-radius)>5e-5:raise ValueError('Baseline radius differs')
    if not mesh_diagnostics(path/'case',base)['valid']:raise ValueError('Baseline mesh rejected')
    parsed=summarize(path/'case',base,g.measurements['diameter']);metrics=accepted_metrics(parsed)
    if metrics is None or metrics[0]<=0:raise ValueError('Baseline is not an accepted positive-thrust result')
    return {'id':'baseline','q':free_coordinates(base),'parameters':base['parameters'],'metrics':metrics,'path':str(path),
            'result_sha256':hashlib.sha256((path/'result.json').read_bytes()).hexdigest()}


def validate_settings(s):
    required={'base_config','baseline','radius','seed','workers','memory_budget_gb','job_memory_gb','max_hours',
              'candidate_hours','max_evaluations','max_attempts','initial_evaluations','candidate_pool',
              'max_consecutive_failures','stagnation_patience','min_accepted','improvement_tolerance','min_free_gb','max_growth_gb'}
    if set(s)!=required:raise ValueError('Unexpected/missing campaign settings')
    if not 0<s['max_hours']<=168 or not 0<s['candidate_hours']<=s['max_hours']:raise ValueError('Invalid wall-clock limits')
    for key in ('workers','max_evaluations','max_attempts','initial_evaluations','candidate_pool','max_consecutive_failures','stagnation_patience','min_accepted'):
        if not isinstance(s[key],int) or s[key]<1:raise ValueError('Invalid '+key)
    if s['workers']>4 or s['workers']*s['job_memory_gb']>s['memory_budget_gb'] or s['memory_budget_gb']>12:
        raise ValueError('Pilot memory/worker limit exceeded')
    if s['initial_evaluations']>=s['max_evaluations']:raise ValueError('Leave evaluations for Bayesian proposals')
    if s['max_attempts']<s['max_evaluations']:raise ValueError('Attempt budget too small')
    for key in ('radius','memory_budget_gb','job_memory_gb','improvement_tolerance','min_free_gb','max_growth_gb'):
        if not isinstance(s[key],(int,float)) or not np.isfinite(s[key]) or s[key]<=0:raise ValueError('Invalid '+key)


def prepare(settings_file,root):
    s=read(settings_file);validate_settings(s);root=Path(root).resolve()
    if root.exists() and any(root.iterdir()):raise FileExistsError('Campaign directory must be empty')
    root.mkdir(parents=True,exist_ok=True)
    started=time.time();base=load_config(s['base_config']);baseline=baseline_record(base,s['baseline'],s['radius'])
    atomic_json(root/'base.json',base);atomic_json(root/'settings.json',s)
    state={'schema_version':1,'status':'PREPARED','started':started,'deadline':started+s['max_hours']*3600,
           'settings':s,'source':source_fingerprint(),'python':sys.executable,'base_sha256':hashlib.sha256(canonical(base).encode()).hexdigest(),
           'observations':[baseline],'baseline_metrics':baseline['metrics'],'jobs':[],'evaluations':0,
           'geometry_rejections':0,'failure_streak':0,'stagnation':0,'best_hypervolume':1.,
           'pool_rejections':{},'rng':np.random.default_rng(s['seed']).bit_generator.state}
    (root/'jobs').mkdir();write_status(root,state)
    return state


def owned_stop(process):
    if process.poll() is not None:return
    process.terminate()
    try:process.wait(timeout=8)
    except subprocess.TimeoutExpired:process.kill();process.wait(timeout=5)


def worker(jobfile):
    from ..simulation.evaluator import evaluate_cfd
    job=read(jobfile)
    def terminate(signum,frame):raise KeyboardInterrupt()
    signal.signal(signal.SIGTERM,terminate)
    os.environ['PROPELLER_DEADLINE_EPOCH']=str(job['deadline'])
    os.environ['PROPELLER_SUPERVISOR_PID']=str(job['supervisor_pid'])
    try:
        result=evaluate_cfd(job['config'],job['evaluations_root'],radius_constraint={'radius':job['radius'],'tolerance':5e-5})
    except KeyboardInterrupt:
        result={'status':'INTERRUPTED','reasons':['Worker stopped by supervisor/deadline']}
    except Exception as exc:result={'status':'WORKER_FAILED','reasons':[type(exc).__name__+': '+str(exc)]}
    atomic_json(Path(jobfile).with_suffix('.result.json'),result)


def disk_usage(root):
    size=0
    for path in root.rglob('*'):
        try:
            if path.is_file():size+=path.stat().st_size
        except FileNotFoundError:pass  # Atomic state-file replacement is normal.
    return size/1024**3


def ingest(state,j,result):
    j['status']=result['status'];j['result']=result.get('path');j['reasons']=result.get('reasons',[])
    j['seconds']=result.get('seconds')
    j['peak_child_gb']=max((c.get('cumulative_child_max_rss_bytes',0) for c in result.get('commands',[])),default=0)/1024**3
    if result['status']=='GEOMETRY_INVALID':state['geometry_rejections']+=1
    else:state['evaluations']+=1
    metrics=accepted_metrics(result)
    if metrics is not None:
        state['observations'].append({'id':j['id'],'q':j['q'],'parameters':j['parameters'],'metrics':metrics,'path':result['path']})
        state['failure_streak']=0
        hv=hypervolume(outputs(state['observations'],state['baseline_metrics']))
        gain=(hv-state['best_hypervolume'])/max(state['best_hypervolume'],1e-12)
        state['stagnation']=0 if gain>state['settings']['improvement_tolerance'] else state['stagnation']+1
        state['best_hypervolume']=max(hv,state['best_hypervolume'])
    elif result['status']=='SUCCESS':state['failure_streak']=0
    elif result['status']!='GEOMETRY_INVALID':state['failure_streak']+=1


def worker_limit(state):
    settings=state['settings']
    measured=max((j.get('peak_child_gb',0) for j in state['jobs']),default=0)
    estimate=max(settings['job_memory_gb'],measured*1.5)
    capacity=int(settings['memory_budget_gb']//estimate)
    # Establish resource use on two jobs before allowing greater concurrency.
    warmup=2 if len(state['observations'])<3 else settings['workers']
    return min(settings['workers'],capacity,warmup)


def run(root):
    root=Path(root).resolve();lock=(root/'supervisor.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=read(root/'state.json');s=state['settings'];base=read(root/'base.json')
    if state['source']!=source_fingerprint():raise ValueError('Campaign source changed; preserve snapshot before resuming')
    if hashlib.sha256(canonical(base).encode()).hexdigest()!=state['base_sha256']:raise ValueError('Campaign config changed')
    # Explicit restart is allowed only after all old workers have exited. Never
    # launch a duplicate while an orphaned worker may still own a native solver.
    for j in state['jobs']:
        if j['status']=='RUNNING' and process_alive(j.get('pid')):raise ValueError('Prior worker still alive')
    rng=np.random.default_rng();rng.bit_generator.state=state['rng']
    running={};state.update(status='RUNNING',pid=os.getpid());write_status(root,state)
    def stop_signal(signum,frame):raise KeyboardInterrupt()
    signal.signal(signal.SIGTERM,stop_signal)
    def alarm(signum,frame):raise TimeoutError('Proposal time budget exceeded')
    signal.signal(signal.SIGALRM,alarm)
    reason=None
    try:
        while True:
            # Ingest completed results, including acknowledged work after restart.
            for j in state['jobs']:
                if j['status']!='RUNNING':continue
                resultfile=root/'jobs'/(j['id']+'.result.json')
                proc=running.get(j['id'])
                if proc is not None and proc.poll() is None:continue
                if resultfile.exists():r=read(resultfile)
                elif proc is None or proc.poll() is not None:r={'status':'WORKER_LOST','reasons':['No worker acknowledgement']}
                else:continue
                running.pop(j['id'],None);ingest(state,j,r)
                write_status(root,state)
            now=time.time()
            if (root/'STOP').exists():reason='user stop requested'
            elif now>=state['deadline']:reason='wall-clock ceiling reached'
            elif state['failure_streak']>=s['max_consecutive_failures']:reason='repeated numerical failures; inspect before spending more'
            elif state.get('proposal_failures',0)>=2:reason='repeated Bayesian proposal failures; inspect implementation'
            elif worker_limit(state)<1:reason='measured per-job memory exceeds campaign allowance'
            elif len(state['observations'])>=s['min_accepted'] and state['stagnation']>=s['stagnation_patience']:reason='accepted frontier improvement stalled'
            elif state['evaluations']>=s['max_evaluations']:reason='pilot evaluation budget reached'
            elif len(state['jobs'])>=s['max_attempts']:reason='attempt budget reached'
            elif shutil.disk_usage(root).free/1024**3<s['min_free_gb'] or disk_usage(root)>s['max_growth_gb']:reason='disk reserve reached'
            # Budget/no-improvement stops drain active work; explicit stop/deadline
            # stops owned workers immediately. No new work is launched meanwhile.
            if reason and (not running or reason in ('user stop requested','wall-clock ceiling reached','disk reserve reached')):break
            if reason:
                state['status']='DRAINING';state['stop_reason']=reason;write_status(root,state);time.sleep(2);continue
            state['worker_limit']=worker_limit(state)
            if len(running)>=state['worker_limit'] or state['evaluations']+len(running)>=s['max_evaluations']:
                time.sleep(2);continue
            remaining=state['deadline']-now
            if remaining<s['candidate_hours']*3600:
                reason='insufficient time reserved for another complete candidate'
                if not running:break
                state['status']='DRAINING';state['stop_reason']=reason;write_status(root,state);time.sleep(2);continue
            existing=[o['q'] for o in state['observations']]+[j['q'] for j in state['jobs']]
            pool,configs,reject=candidate_pool(base,s['radius'],rng,existing,s['candidate_pool'])
            for code,count in reject.items():state['pool_rejections'][code]=state['pool_rejections'].get(code,0)+count
            pending=[j['q'] for j in state['jobs'] if j['status']=='RUNNING']
            use_bo=state['evaluations']>=s['initial_evaluations']
            proposal_seed=int(rng.integers(2**31-1));warning=None
            try:
                signal.setitimer(signal.ITIMER_REAL,min(180,max(1,state['deadline']-time.time())))
                index,method=propose(pool,state['observations'],state['baseline_metrics'],pending,proposal_seed,use_bo)
                state['proposal_failures']=0
            except (Exception,) as exc:
                warning=type(exc).__name__+': '+str(exc)
                index,method=propose(pool,state['observations'],state['baseline_metrics'],pending,proposal_seed,False)
                method='space-filling-fallback'
                state['proposal_failures']=state.get('proposal_failures',0)+1
            finally:signal.setitimer(signal.ITIMER_REAL,0)
            if time.time()+s['candidate_hours']*3600>state['deadline']:continue
            identity='design-%03d'%(len(state['jobs'])+1)
            candidate=configs[index]
            j={'id':identity,'q':pool[index].tolist(),'parameters':candidate['parameters'],'status':'RUNNING','method':method,'proposal_warning':warning,'started':time.time()}
            job={'config':candidate,'radius':s['radius'],'supervisor_pid':os.getpid(),
                 'deadline':min(state['deadline'],time.time()+s['candidate_hours']*3600),'evaluations_root':str(root/'evaluations')}
            jobfile=root/'jobs'/(identity+'.json');atomic_json(jobfile,job)
            with (root/'jobs'/(identity+'.log')).open('w') as out:
                process=subprocess.Popen([sys.executable,'-m','propeller.optimization.campaign','worker',str(jobfile)],stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
            j['pid']=process.pid;running[identity]=process;state['jobs'].append(j);state['rng']=rng.bit_generator.state;write_status(root,state)
        state['stop_reason']=reason or 'finished'
    except KeyboardInterrupt:state['stop_reason']='supervisor interrupted'
    except Exception as exc:state['stop_reason']=type(exc).__name__+': '+str(exc);state['status']='FAILED'
    finally:
        for proc in running.values():owned_stop(proc)
        for j in state['jobs']:
            if j['status']=='RUNNING':
                resultfile=root/'jobs'/(j['id']+'.result.json')
                ingest(state,j,read(resultfile) if resultfile.exists() else {'status':'INTERRUPTED'})
        if state['status']!='FAILED':state['status']='STOPPED'
        state['finished']=time.time();write_status(root,state);lock.close()


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['start','run','worker','status','stop']);p.add_argument('path');p.add_argument('--settings')
    a=p.parse_args();root=Path(a.path).resolve()
    if a.action=='worker':worker(root);return
    if a.action=='run':run(root);return
    if a.action=='status':print(json.dumps(read(root/'state.json'),indent=2));return
    if a.action=='stop':(root/'STOP').touch();print('Stop requested');return
    if not a.settings:raise ValueError('start requires --settings')
    prepare(a.settings,root)
    env=os.environ.copy();env['PYTHONPATH']=str(Path(__file__).resolve().parents[2]);env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1';env['VECLIB_MAXIMUM_THREADS']='1'
    with (root/'supervisor.log').open('a') as out:
        proc=subprocess.Popen([sys.executable,'-m','propeller.optimization.campaign','run',str(root)],stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True,env=env,cwd=str(Path.cwd()))
    atomic_json(root/'launch.json',{'pid':proc.pid,'started':time.time(),'command':'detached pilot supervisor'})
    print(json.dumps({'campaign':str(root),'pid':proc.pid,'report':str(root/'report.html')}))


if __name__=='__main__':main()

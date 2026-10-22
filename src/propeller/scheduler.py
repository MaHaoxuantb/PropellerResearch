"""Bounded independent candidate execution; a single SA chain stays serial."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing
import os
import math
from pathlib import Path
from .config import load_config
from .experiments import atomic_json


def check_budget(configs,workers,memory_budget_gb=None,job_memory_gb=None,cpu_count=None):
    if not isinstance(workers,int) or workers<1:raise ValueError('workers must be a positive integer')
    cores=cpu_count or os.cpu_count() or 1
    ranks=max(c['cfd'].get('mpi_ranks',1) for c in configs)
    if workers*ranks>cores:raise ValueError('Requested simultaneous MPI ranks exceed available logical CPUs')
    if workers>1:
        if memory_budget_gb is None or job_memory_gb is None or not all(math.isfinite(v) and v>0 for v in (memory_budget_gb,job_memory_gb)):
            raise ValueError('Parallel candidate runs require explicit memory budget and measured per-job memory estimate')
        if workers*job_memory_gb>memory_budget_gb:raise ValueError('Concurrent jobs exceed the memory budget')
    return {'workers':workers,'ranks_per_job_upper_bound':ranks,'logical_cpus':cores,'memory_budget_gb':memory_budget_gb,'job_memory_gb':job_memory_gb}


def _evaluate(job):
    from .simulation.evaluator import evaluate_cfd
    config,root=job
    return evaluate_cfd(config,root)


def evaluate_batch(paths,root,workers=1,memory_budget_gb=None,job_memory_gb=None):
    if not paths:raise ValueError('No batch configurations supplied')
    configs=[load_config(path) for path in paths]
    budget=check_budget(configs,workers,memory_budget_gb,job_memory_gb)
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    manifest={'budget':budget,'jobs':[{'configuration':str(Path(path).resolve()),'status':'QUEUED'} for path in paths]}
    import uuid
    record=root/('batch-'+uuid.uuid4().hex+'.json')
    atomic_json(record,manifest)
    if workers==1:
        for i,c in enumerate(configs):
            manifest['jobs'][i]['status']='RUNNING';atomic_json(record,manifest)
            result=_evaluate((c,root));manifest['jobs'][i].update(status=result['status'],result=result);atomic_json(record,manifest)
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            tasks={pool.submit(_evaluate,(c,root)):i for i,c in enumerate(configs)}
            for job in manifest['jobs']:job['status']='SUBMITTED'
            atomic_json(record,manifest)
            for future in as_completed(tasks):
                i=tasks[future]
                try:result=future.result()
                except Exception as exc:result={'status':'INTERNAL_ERROR','reasons':[str(exc)]}
                manifest['jobs'][i].update(status=result['status'],result=result);atomic_json(record,manifest)
    return {'record':str(record.resolve()),**manifest}

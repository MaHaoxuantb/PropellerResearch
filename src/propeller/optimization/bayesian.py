"""Multi-objective Bayesian proposals over a geometrically screened pool."""
import numpy as np
from ..geometry.screening import fixed_radius_config


def pareto_indices(values):
    y=np.asarray(values,dtype=float)
    if not len(y):return []
    return [i for i,row in enumerate(y) if not np.any(np.all(y>=row,axis=1)&np.any(y>row,axis=1))]


def hypervolume(values):
    """Two positive maximization outputs, with a fixed origin reference."""
    points=np.asarray(values,dtype=float)
    if not len(points):return 0.
    points=points[np.all(np.isfinite(points)&(points>0),axis=1)]
    points=points[pareto_indices(points)]
    area=0.;previous=0.
    for x,y in sorted(points.tolist()):
        area+=(x-previous)*y;previous=x
    return float(area)


def accepted_metrics(result):
    if result.get('status')!='SUCCESS':return None
    m=result.get('metrics',{});t=m.get('thrust');p=m.get('shaft_power')
    # Negative thrust with positive consumed power is a valid observation of a
    # poor propeller, not a failed solve. The acquisition's T>0 constraint and
    # positive baseline keep it off the useful propulsion frontier.
    if t is None or p is None or not np.isfinite([t,p]).all() or p<=0:return None
    return [float(t),float(p)]


def candidate_pool(base,radius,rng,existing,count=256):
    from ..geometry.screening import free_coordinates
    center=np.array(free_coordinates(base));pool=[];configs=[];rejections={}
    for _ in range(count*30):
        # Include broad exploration while retaining a useful local pilot around
        # the accepted baseline. Reflection avoids point masses at the bounds.
        q=rng.random(8) if rng.random()<.25 else center+rng.normal(0,.16,8)
        q=np.mod(q,2);q=np.where(q>1,2-q,q)
        prior=np.asarray(list(existing)+pool)
        if len(prior) and np.min(np.linalg.norm(prior-q,axis=1))<1e-5:continue
        try:c=fixed_radius_config(base,q,radius)
        except ValueError as exc:
            key=str(exc);rejections[key]=rejections.get(key,0)+1;continue
        pool.append(q.tolist());configs.append(c)
        if len(pool)>=count:break
    if not pool:raise ValueError('No section-feasible fixed-radius candidates')
    return np.asarray(pool),configs,rejections


def propose(pool,observations,baseline_metrics,pending,seed,bayesian=True):
    """Return pool index and algorithm label; actual CFD is always authoritative."""
    previous=np.asarray([o['q'] for o in observations]+list(pending))
    if not bayesian or len(observations)<4:
        distances=np.linalg.norm(pool[:,None,:]-previous[None,:,:],axis=2).min(axis=1)
        if len(observations)<3:
            local=np.linalg.norm(pool-np.asarray(observations[0]['q']),axis=1)<=.5
            if np.any(local):
                distances=np.where(local,distances,-np.inf)
                return int(np.argmax(distances)),'space-filling-local'
        return int(np.argmax(distances)),'space-filling'
    import torch
    from botorch.models import SingleTaskGP,ModelListGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.fit import fit_gpytorch_mll
    from gpytorch.mlls import SumMarginalLogLikelihood
    from botorch.acquisition.multi_objective.logei import qLogNoisyExpectedHypervolumeImprovement
    from botorch.acquisition.multi_objective.objective import GenericMCMultiOutputObjective
    from botorch.sampling.normal import SobolQMCNormalSampler
    torch.set_num_threads(1);torch.manual_seed(seed)
    tensor=lambda x:torch.as_tensor(x,dtype=torch.double,device='cpu')
    x=tensor([o['q'] for o in observations]);raw=np.array([o['metrics'] for o in observations])
    t0,p0=baseline_metrics
    y=tensor(np.column_stack((raw[:,0]/t0,np.log(raw[:,1]/p0))))
    models=[SingleTaskGP(x,y[:,i:i+1],outcome_transform=Standardize(m=1)) for i in range(2)]
    model=ModelListGP(*models)
    fit_gpytorch_mll(SumMarginalLogLikelihood(model.likelihood,model),optimizer_kwargs={'options':{'maxiter':75}})
    # Shared sampled T appears in both outputs. P is strictly positive because
    # the second GP models log(P/P0). No unconverged force values enter training.
    objective=GenericMCMultiOutputObjective(lambda samples,X=None:torch.stack((samples[...,0],samples[...,0]*torch.exp(-samples[...,1])),dim=-1))
    acq=qLogNoisyExpectedHypervolumeImprovement(model=model,ref_point=[0.,0.],X_baseline=x,
        objective=objective,constraints=[lambda samples:-samples[...,0]],
        sampler=SobolQMCNormalSampler(sample_shape=torch.Size([64]),seed=seed),
        X_pending=tensor(pending) if len(pending) else None,prune_baseline=True)
    values=[]
    with torch.no_grad():
        for batch in torch.split(tensor(pool),32):values.append(acq(batch.unsqueeze(-2)))
    scores=torch.cat(values)
    if not torch.isfinite(scores).any():raise ValueError('Acquisition has no finite scores')
    return int(torch.argmax(scores).item()),'qLogNEHVI'

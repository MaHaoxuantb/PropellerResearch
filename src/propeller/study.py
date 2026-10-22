"""Matched baseline/mesh/domain evidence, with explicit exploratory gates."""
import hashlib
import html
import json
import math
import os
from pathlib import Path

from .analysis import assess_baseline
from .experiments import atomic_json
from .simulation.parser import read_vectors, NUMBER


def _plot_history(result_path):
    """Join same-case continuations; mapped initial guesses start a new history."""
    import re
    import numpy as np
    chain=[];seen=set();current=Path(result_path).resolve()
    while current.exists():
        if current in seen or len(chain)>=100:raise ValueError('Invalid continuation history chain')
        seen.add(current);chain.append(current)
        result=json.loads(current.read_text());parent=result.get('continuation',{}).get('parent')
        if not parent:break
        current=Path(parent)/'result.json'
    forces={};moments={};residuals={}
    for path in reversed(chain):
        case=path.parent/'case'
        folders=[p for p in (case/'postProcessing/forces').glob('*') if (p/'force.dat').exists() and (p/'moment.dat').exists()]
        if not folders:continue
        folder=max(folders,key=lambda p:float(p.name));force=read_vectors(folder/'force.dat');moment=read_vectors(folder/'moment.dat')
        if force.shape!=moment.shape or not np.array_equal(force[:,0],moment[:,0]):raise ValueError('Study force/moment timestamps differ')
        # A restart may replay from a saved field before its parent's last log.
        for table in (forces,moments,residuals):
            for time in [t for t in table if t>=force[0,0]]:del table[time]
        forces.update((row[0],row) for row in force);moments.update((row[0],row) for row in moment)
        iteration=None
        for line in (case/'log.simpleFoam').read_text(errors='replace').splitlines():
            match=re.match(r'Time = ('+NUMBER+')',line)
            if match:iteration=float(match.group(1))
            match=re.search(r'Solving for (\w+), Initial residual = ('+NUMBER+')',line)
            if match and iteration is not None:residuals[iteration]=max(residuals.get(iteration,0),float(match.group(2)))
    times=sorted(forces)
    return np.array([forces[t] for t in times]),np.array([moments[t] for t in times]),residuals,[str(p) for p in reversed(chain)]


def study_summary(paths, window=500, relative_limit=.05):
    if len(paths)<3:
        raise ValueError("Provide a baseline, a finer-mesh result and a larger-domain result")
    if not 0<relative_limit<1:
        raise ValueError("relative_limit must lie in (0,1)")
    if not isinstance(window,int) or window<4:raise ValueError('Assessment window must be an integer >=4')
    if len(set(str(Path(p).resolve()) for p in paths))!=len(paths):raise ValueError('Study requires distinct evaluation records')
    rows=[];configs=[];reference=None;reference_dictionaries=None;reference_runtime=None
    for path in paths:
        path=Path(path).resolve();result=json.loads(path.read_text())
        manifest=json.loads((path.parent/'manifest.json').read_text());config=manifest['config']
        configs.append(config)
        physics={k:config[k] for k in ('design','parameters','geometry')}
        physics['operating']={k:v for k,v in config['operating'].items() if k not in ('target_thrust','reference_power')}
        physics['acceptance']={k:config['cfd'][k] for k in ('residual_tolerance','continuity_tolerance','force_relative_tolerance','max_non_orthogonality','max_skewness')}
        dictionaries={name:hashlib.sha256((path.parent/'case'/name).read_bytes()).hexdigest()
                      for name in ('system/fvSchemes','system/fvSolution','constant/transportProperties','constant/turbulenceProperties','constant/MRFProperties')}
        runtime=manifest['runtime']['binary_sha256']
        if reference is None:
            reference=physics;reference_dictionaries=dictionaries;reference_runtime=runtime
        elif physics!=reference or dictionaries!=reference_dictionaries or runtime!=reference_runtime:
            raise ValueError('Study physics, solver dictionaries or native runtime differ: '+str(path))
        assessment=None
        if result.get('metrics'):
            assessment=assess_baseline(path,window)
        row={'source':str(path),'fidelity':config['cfd']['fidelity'],'status':result['status'],
             'mesh':result.get('mesh'),'assessment':assessment,'numerical_metrics':result.get('metrics'),
             'reasons':result.get('reasons',[]),'relative_change':None,'role':'initial baseline' if not rows else 'other',
             'continuation_parent':manifest.get('continuation',{}).get('parent')}
        rows.append(row)
    def dimensions(c):
        d=c['cfd'].get('domain',{'radial_diameters':3,'upstream_diameters':3,'downstream_diameters':8})
        return d,[2*d['radial_diameters']/c['cfd']['background_cells'][0],2*d['radial_diameters']/c['cfd']['background_cells'][1],(d['upstream_diameters']+d['downstream_diameters'])/c['cfd']['background_cells'][2]]
    base_domain,base_spacing=dimensions(configs[0]);finer=False;larger=False
    base_h=math.prod(base_spacing)**(1/3)
    for row,config in zip(rows,configs):
        _,spacing=dimensions(config)
        row['nominal_spacing_m']=[v*config['design']['diameter_max'] for v in spacing]
        row['resolution_factor']=base_h/(math.prod(spacing)**(1/3))
        row['label']=row['role']
    for index,config in enumerate(configs[1:],1):
        domain,spacing=dimensions(config)
        is_finer=domain==base_domain and all(a<b for a,b in zip(spacing,base_spacing))
        is_larger=all(domain[k]>=base_domain[k] for k in domain) and any(domain[k]>base_domain[k] for k in domain) and all(math.isclose(a,b,rel_tol=1e-8) for a,b in zip(spacing,base_spacing))
        finer |= is_finer;larger |= is_larger
        if is_finer:
            ratios=[a/b for a,b in zip(spacing,base_spacing)]
            if not all(math.isclose(v,ratios[0],rel_tol=1e-8) for v in ratios):raise ValueError('Use a uniform background refinement family')
            rows[index]['role']='finer mesh'
        elif is_larger:rows[index]['role']='larger domain'
        elif {k:v for k,v in config['cfd'].items() if k not in ('iterations','timeout_seconds')}=={k:v for k,v in configs[0]['cfd'].items() if k not in ('iterations','timeout_seconds')}:
            rows[index]['role']='same-mesh control'
        rows[index]['label']='mesh x'+format(rows[index]['resolution_factor'],'.3g') if is_finer else rows[index]['role']
    if not finer or not larger:raise ValueError('Study requires a finer background mesh and a larger domain at matched spacing')
    controls=[i for i,r in enumerate(rows) if r['role']=='same-mesh control']
    if len(controls)>1:raise ValueError('Select one same-mesh control reference')
    reference_index=controls[0] if controls else 0
    base=rows[reference_index]['assessment']
    essential=[r for i,r in enumerate(rows) if i==reference_index or r['role'] in ('finer mesh','larger domain')]
    incomplete=any(not r['assessment'] or r['assessment']['status']!='READY_FOR_SENSITIVITY' for r in essential)
    if base and base['status']=='READY_FOR_SENSITIVITY':
        for index,row in enumerate(rows):
            if index==reference_index:continue
            if row['assessment'] and row['assessment']['status']=='READY_FOR_SENSITIVITY':
                row['relative_change']={key:(row['assessment']['statistics'][key]['mean']/base['statistics'][key]['mean']-1)
                                        for key in ('thrust','shaft_power') if base['statistics'][key]['mean']>0}
    domain_ok=all(r['relative_change'] and all(abs(v)<=relative_limit for v in r['relative_change'].values()) for r in rows if r['role']=='larger domain')
    grid_indices=[reference_index]+[i for i,r in enumerate(rows) if r['role']=='finer mesh']
    grid_indices.sort(key=lambda i:rows[i]['resolution_factor'])
    factors=[round(rows[i]['resolution_factor'],10) for i in grid_indices]
    if len(set(factors))!=len(factors):raise ValueError('Select one latest record per grid resolution')
    for key in ('surface_refinement','rotor_refinement','wake_refinement','layers'):
        if any(configs[i]['cfd'][key]!=configs[0]['cfd'][key] for i in grid_indices):raise ValueError('Uniform refinement study requires matched '+key)
    def difference(coarser,finer):
        a,b=rows[coarser]['assessment'],rows[finer]['assessment']
        if not a or not b or any(v['status']!='READY_FOR_SENSITIVITY' for v in (a,b)):return None
        return {key:b['statistics'][key]['mean']/a['statistics'][key]['mean']-1 for key in ('thrust','shaft_power')}
    pairs=[]
    for a,b in zip(grid_indices,grid_indices[1:]):
        change=difference(a,b)
        cells_a=(rows[a]['mesh'] or {}).get('cells');cells_b=(rows[b]['mesh'] or {}).get('cells')
        if cells_a and cells_b and cells_b<=cells_a:raise ValueError('Refined case has not increased cell count')
        pairs.append({'coarser':rows[a]['source'],'finer':rows[b]['source'],'relative_change':change,
                      'within_pilot_limit':change is not None and all(abs(v)<=relative_limit for v in change.values())})
    selected=None
    if not incomplete and domain_ok:
        for position,index in enumerate(grid_indices[:-1]):
            change=difference(index,grid_indices[-1])
            if all(p['within_pilot_limit'] for p in pairs[position:]) and change and all(abs(v)<=relative_limit for v in change.values()):
                selected=index;break
    settling=None
    if controls and rows[0]['assessment'] and base:
        settling={key:base['statistics'][key]['mean']/rows[0]['assessment']['statistics'][key]['mean']-1
                  for key in ('thrust','shaft_power') if rows[0]['assessment']['statistics'][key]['mean']>0}
    return {'status':'VALIDATION_INCOMPLETE' if incomplete else 'SENSITIVITY_TOO_LARGE' if selected is None else 'READY_FOR_EXPLORATORY_SEARCH',
            'window':window,'relative_limit':relative_limit,'cases':rows,'qualified':False,
            'comparison_reference':rows[reference_index]['source'],'same_mesh_settling_change':settling,
            'grid_pairs':pairs,'domain_within_pilot_limit':bool(domain_ok),
            'recommended_search_source':rows[selected]['source'] if selected is not None else None,
            'confirmation_source':rows[grid_indices[-1]]['source'],
            'interpretation':'Matched operating point and solver equations. Uniform grid refinements and a domain enlargement are a pilot sensitivity screen, not a formal uncertainty estimate or physical qualification. Search fidelity is selected only from a stable finest-grid sequence and must agree with the finest tested grid. Domain sensitivity is checked at the original control resolution. Force settling remains part of the comparison.'}


def write_study(paths, output, window=500, relative_limit=.05):
    summary=study_summary(paths,window,relative_limit)
    output=Path(output).resolve()
    if output.suffix.lower()!='.html':raise ValueError('Study report output must end in .html')
    output.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(output.with_suffix('.json'),summary)
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/propeller-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True,layout='constrained')
    colors=('#147d92','#a65d1b','#7b4daf')
    baseline_end=None
    baseline_config=json.loads((Path(summary['cases'][0]['source']).parent/'manifest.json').read_text())['config']
    for index,row in enumerate(summary['cases']):
        if not row['assessment']:continue
        case=Path(row['source']).parent/'case'
        force,moment,residuals,segments=_plot_history(row['source'])
        if not len(force):continue
        row['history_segments']=segments
        config=json.loads((case/'configuration.json').read_text())
        if force.shape!=moment.shape or not (force[:,0]==moment[:,0]).all():raise ValueError('Study force/moment timestamps differ')
        if index==0:baseline_end=force[-1,0]
        is_handoff_control=row['role']=='same-mesh control' and row['continuation_parent']==str(Path(summary['cases'][0]['source']).parent)
        shift=baseline_end if index==0 or is_handoff_control else 0
        # Baseline's final window lies left of zero; new meshes start at zero.
        keep=force[:,0]>=baseline_end-window+1 if index==0 else force[:,0]>baseline_end if is_handoff_control else force[:,0]>=0
        x=force[keep,0]-shift
        omega=config['operating']['rotation_sign']*config['operating']['rpm']*math.pi/30
        color='#147d92' if row['role']=='same-mesh control' else '#568d37' if row['role']=='larger domain' else colors[index%len(colors)];label=row['label']
        if row['assessment']['status']!='READY_FOR_SENSITIVITY':label+=' [unaccepted]'
        style='--' if row['role']=='same-mesh control' else '-'
        axes[0].plot(x,force[keep,3],label=label,color=color,linewidth=1.2,linestyle=style)
        axes[1].plot(x,-moment[keep,3]*omega,color=color,linewidth=1.2,linestyle=style)
        samples=[(t-shift,r) for t,r in sorted(residuals.items()) if (index!=0 or t>=baseline_end-window+1) and (not is_handoff_control or t>baseline_end)]
        if samples:axes[2].semilogy([s[0] for s in samples],[max(s[1],1e-16) for s in samples],color=color,linewidth=.8)
    axes[0].set_ylabel('Thrust [N]');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Shaft power [W]')
    axes[2].set_ylabel('Largest initial residual');axes[2].set_xlabel('Steady iteration relative to baseline handoff (not physical time)')
    axes[2].axhline(baseline_config['cfd']['residual_tolerance'],color='#9f3d32',linestyle='--',linewidth=.9,label='Residual acceptance limit')
    axes[2].legend(fontsize=8)
    for axis in axes:axis.grid(alpha=.2);axis.axvline(0,color='#697586',linewidth=.7,linestyle='--')
    fig.suptitle('Baseline sensitivity evidence · '+summary['status'],fontsize=13)
    fig.savefig(output.with_suffix('.svg'));fig.savefig(output.with_suffix('.png'),dpi=150);plt.close(fig)
    atomic_json(output.with_suffix('.json'),summary)
    rows=[]
    for i,row in enumerate(summary['cases']):
        stats=row['assessment']['statistics'] if row['assessment'] else {}
        metrics='; '.join(k+': '+format(v['mean'],'.5g') for k,v in stats.items()) or 'Unavailable'
        if row['assessment'] and row['assessment']['status']!='READY_FOR_SENSITIVITY':metrics+=' (unaccepted)'
        delta='; '.join(k+': '+format(100*v,'+.2f')+'%' for k,v in (row['relative_change'] or {}).items()) or ('Reference' if row['source']==summary['comparison_reference'] else 'Not accepted')
        rows.append('<tr><td>'+html.escape(row['label'])+'</td><td>'+html.escape(row['status'])+'</td><td>'+html.escape(metrics)+'</td><td>'+html.escape(delta)+'</td></tr>')
    reference_note=' The same-mesh control is the comparison reference, separating further settling from mesh/domain changes.' if summary['same_mesh_settling_change'] is not None else ''
    labels={r['source']:r['label'] for r in summary['cases']}
    grid_rows=[]
    for pair in summary['grid_pairs']:
        delta='; '.join(k+': '+format(100*v,'+.2f')+'%' for k,v in (pair['relative_change'] or {}).items()) or 'Not accepted'
        grid_rows.append('<tr><td>'+html.escape(labels[pair['coarser']]+' → '+labels[pair['finer']])+'</td><td>'+html.escape(delta)+'</td><td>'+('Pass' if pair['within_pilot_limit'] else 'Not passed')+'</td></tr>')
    selection='Selected pilot search fidelity: '+labels[summary['recommended_search_source']]+'.' if summary['recommended_search_source'] else 'No search fidelity selected yet.'
    output.write_text('<!doctype html><meta charset="utf-8"><title>Baseline sensitivity study</title><style>body{max-width:1100px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#20344a}img{width:100%}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #d8dfe7;padding:10px;text-align:left}pre{white-space:pre-wrap;font-size:12px}</style><h1>Baseline sensitivity study</h1><p>Status: <strong>'+summary['status']+'</strong></p><p>'+html.escape(summary['interpretation']+reference_note)+'</p><p>Means and changes use the final '+str(window)+' samples in every case. Pilot change limit: '+format(100*relative_limit,'.1f')+'%.</p><p><strong>'+html.escape(selection)+'</strong></p><table><tr><th>Case</th><th>Numerical status</th><th>Window means [N, W]</th><th>Change from control</th></tr>'+''.join(rows)+'</table><table><tr><th>Grid refinement</th><th>Change</th><th>Pilot limit</th></tr>'+''.join(grid_rows)+'</table><img alt="Thrust, power and residual sensitivity histories" src="'+output.with_suffix('.svg').name+'"><details><summary>Full assessment and provenance</summary><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre></details>')
    return output,summary

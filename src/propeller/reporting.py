"""Standalone numerical evidence reports with standard Matplotlib figures."""
import csv
import html
import json
import os
from pathlib import Path
import re
from .simulation.parser import read_vectors,NUMBER


def write_report(result_path,output):
    result_path=Path(result_path).resolve();output=Path(output).resolve()
    result=json.loads(result_path.read_text());output.parent.mkdir(parents=True,exist_ok=True)
    case=result_path.parent/'case'
    configuration=case/'configuration.json'
    plot=None
    if configuration.exists():
        config=json.loads(configuration.read_text())
        folders=[p for p in (case/'postProcessing/forces').glob('*') if (p/'force.dat').exists() and (p/'moment.dat').exists()]
        if folders:
            folder=max(folders,key=lambda p:float(p.name))
            force,moment=read_vectors(folder/'force.dat'),read_vectors(folder/'moment.dat')
            if force.shape!=moment.shape or not (force[:,0]==moment[:,0]).all():raise ValueError('Force/moment histories differ')
            import numpy as np
            os.environ.setdefault('MPLCONFIGDIR','/private/tmp/propeller-matplotlib')
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
            power=-moment[:,3]*config['operating']['rotation_sign']*config['operating']['rpm']*np.pi/30
            fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True,layout='constrained')
            axes[0].plot(force[:,0],force[:,3],color='#147d92',linewidth=1.2)
            axes[0].set_ylabel('Thrust [N]')
            axes[1].plot(force[:,0],power,color='#9355b5',linewidth=1.2)
            axes[1].set_ylabel('Shaft power [W]')
            log=(case/'log.simpleFoam').read_text(errors='replace')
            samples={};iteration=None
            for line in log.splitlines():
                match=re.match(r'Time = ('+NUMBER+')',line)
                if match:iteration=float(match.group(1))
                match=re.search(r'Solving for (\w+), Initial residual = ('+NUMBER+')',line)
                if match and iteration is not None:
                    key=(iteration,match.group(1));samples[key]=max(samples.get(key,0),float(match.group(2)))
            for field in ('Ux','Uy','Uz','p','k','omega'):
                rows=sorted((t,r) for (t,f),r in samples.items() if f==field)
                if rows:axes[2].semilogy([r[0] for r in rows],[max(r[1],1e-16) for r in rows],label=field,linewidth=.8)
            axes[2].axhline(config['cfd']['residual_tolerance'],color='#b04c35',linestyle='--',label='acceptance threshold')
            axes[2].set_ylabel('Initial residual');axes[2].set_xlabel('Steady solver iteration');axes[2].legend(ncol=4,fontsize=8)
            window=config['cfd']['force_window'];start=force[max(0,len(force)-window),0]
            for axis in axes:
                axis.grid(alpha=.2);axis.axvspan(start,force[-1,0],alpha=.08,color='#354b62')
            fig.suptitle('Native OpenFOAM evidence · '+result.get('status','UNKNOWN'),fontsize=14)
            plot=output.with_suffix('.svg');fig.savefig(plot);fig.savefig(output.with_suffix('.png'),dpi=150);plt.close(fig)
            with output.with_suffix('.csv').open('w',newline='') as stream:
                writer=csv.writer(stream);writer.writerow(['iteration','thrust_N','fluid_torque_z_Nm','shaft_power_W'])
                writer.writerows(zip(force[:,0],force[:,3],moment[:,3],power))
    status=html.escape(result.get('status','UNKNOWN'))
    parts=['<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Rotor evaluation report</title>',
           '<style>body{max-width:1000px;margin:40px auto;padding:0 24px;color:#20344a;font:16px/1.6 system-ui;background:#f6f8fa}h1{line-height:1.2}section{background:white;border:1px solid #dde4eb;padding:24px;border-radius:10px;margin:20px 0}pre{white-space:pre-wrap;font-size:12px;overflow-wrap:anywhere}img{width:100%}.notice{border-left:4px solid #c78d3f;padding:12px 20px;background:#fff3df}a{color:#147d92}</style>',
           '<h1>Rotor evaluation</h1><p>Numerical status: <strong>'+status+'</strong></p>',
           '<p class="notice">This is a pilot CFD evaluation. Numerical convergence, if achieved, does not establish grid independence, experimental accuracy, structural adequacy or an optimized design.</p>']
    if result.get('reasons'):parts.append('<section><h2>Acceptance findings</h2><pre>'+html.escape(json.dumps(result['reasons'],indent=2))+'</pre></section>')
    if plot:parts.append('<section><h2>Convergence evidence</h2><p>The shaded region is the final force-averaging window.</p><img alt="Thrust, power and solver residual histories" src="'+html.escape(plot.name)+'"></section>')
    for key in ('metrics','geometry','mesh','diagnostics'):
        if key in result:parts.append('<section><h2>'+key.title()+'</h2><pre>'+html.escape(json.dumps(result[key],indent=2))+'</pre></section>')
    parts.append('<p>Source record: <a href="'+html.escape(os.path.relpath(result_path,output.parent))+'">'+html.escape(result_path.name)+'</a></p>')
    output.write_text('\n'.join(parts))
    return output

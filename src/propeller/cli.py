"""Command-line research workflow; every artifact is local and inspectable."""
import argparse
import copy
import json
import os
from pathlib import Path
import platform
import sys
import time
import numpy as np
from .config import load_config,physical,with_parameters
from .experiments import atomic_json,versions
from .geometry import generate_rotor,validate_geometry
from .geometry.preview import write_preview
from .simulation.runner import runtime_info
from .simulation.case_builder import build_case
from .simulation.evaluator import evaluate_cfd
from .optimization.annealing import optimize


def emit(value):print(json.dumps(value,indent=2,allow_nan=False))


def main(argv=None):
    parser=argparse.ArgumentParser(description="Toroidal rotor geometry, native OpenFOAM and reproducible optimization")
    parser.add_argument("--config",default="config/demo.json")
    sub=parser.add_subparsers(dest="command",required=True)
    sub.add_parser("doctor",help="Inspect Python and native Homebrew OpenFOAM")
    geo=sub.add_parser("geometry");geos=geo.add_subparsers(dest="action",required=True)
    for action in ("preview","validate","export"):
        p=geos.add_parser(action);p.add_argument("--output",default="runs/baseline");p.add_argument("--fast",action="store_true",help="Skip expensive self-intersection scan; never valid for CFD")
    p=geos.add_parser("sample");p.add_argument("--count",type=int,default=200);p.add_argument("--seed",type=int,default=42);p.add_argument("--output",default="runs/geometry-sample.json")
    p=sub.add_parser("case",help="Generate a native OpenCFD case without running it");p.add_argument("--output",required=True)
    p=sub.add_parser("evaluate");p.add_argument("--output",default="runs/evaluations");p.add_argument("--no-cache",action="store_true");p.add_argument("--initialize-from",help="Map a converged identical rotor/operating point onto a new sensitivity mesh")
    p=sub.add_parser("continue-evaluation",help="Continue an archived serial solve in a new attempt");p.add_argument("parent");p.add_argument("--output",default="runs/evaluations")
    for name in ("optimize","resume"):
        p=sub.add_parser(name);p.add_argument("--checkpoint",default="runs/optimization/checkpoint.json");p.add_argument("--output",default="runs/evaluations");p.add_argument("--exploratory",action="store_true",help="Explicitly run an unqualified pilot; results are not confirmed improvements")
        p.add_argument('--baseline-reference',help='Accepted CFD reference written by freeze-baseline');p.add_argument('--fresh-baseline',action='store_true',help='Re-evaluate the initial design instead of using its accepted reference')
    p=sub.add_parser("report");p.add_argument("result");p.add_argument("--output",default="runs/report.html")
    p=sub.add_parser("compare");p.add_argument("results",nargs="+");p.add_argument("--output",default="runs/comparison.json")
    p=sub.add_parser("assess-baseline",help="Check longer-window settling before baseline promotion");p.add_argument("result");p.add_argument("--window",type=int,default=500);p.add_argument("--output",default="runs/baseline-assessment.json")
    p=sub.add_parser("study-report",help="Compare matched baseline/mesh/domain evidence");p.add_argument("results",nargs="+");p.add_argument("--window",type=int,default=500);p.add_argument("--relative-limit",type=float,default=.05);p.add_argument("--output",default="runs/reports/baseline-study.html")
    p=sub.add_parser("surface-audit",help="Compare the native CFD rotor wall with the generated solid");p.add_argument("attempt");p.add_argument("--output",required=True)
    p=sub.add_parser("freeze-baseline",help="Freeze target thrust and reference power from a successful baseline");p.add_argument("result");p.add_argument("--output",required=True);p.add_argument('--study',help='Require this result to be the selected pilot fidelity in a study JSON');p.add_argument('--calibrate-timeout',action='store_true',help='Increase the full-run timeout from measured iteration cost with a 50% buffer')
    p=sub.add_parser("recover",help="Record orphaned attempts as interrupted without deleting artifacts");p.add_argument("--output",default="runs/evaluations")
    p=sub.add_parser("status",help="Inspect recorded attempts and current solver iteration");p.add_argument("--output",default="runs/evaluations")
    p=sub.add_parser("batch",help="Evaluate independent configurations within explicit resource limits");p.add_argument("configs",nargs="+");p.add_argument("--output",default="runs/evaluations");p.add_argument("--workers",type=int,default=1);p.add_argument("--memory-budget-gb",type=float);p.add_argument("--job-memory-gb",type=float)
    args=parser.parse_args(argv)
    try:
        c=load_config(args.config)
        if args.command=="doctor":
            info={"python":sys.version,"platform":platform.platform(),"cpu_count":os.cpu_count(),"dependencies":versions(),"openfoam":runtime_info(c)}
            emit(info);return 0 if info["openfoam"]["available"] else 2
        if args.command=="recover":
            from .experiments import Store
            store=Store(args.output)
            try:emit({"recovered":store.recover_interrupted()})
            finally:store.close()
            return 0
        if args.command=="status":
            from .progress import run_status
            emit(run_status(args.output));return 0
        if args.command=="batch":
            from .scheduler import evaluate_batch
            result=evaluate_batch(args.configs,args.output,args.workers,args.memory_budget_gb,args.job_memory_gb)
            emit(result);return 0 if all(j['status']=='SUCCESS' for j in result['jobs']) else 2
        if args.command=="geometry" and args.action=="sample":
            if not 1<=args.count<=10000:raise ValueError("count must lie in [1,10000]")
            rng=np.random.default_rng(args.seed);rows=[]
            # Sample section geometry first; avoid expensive Booleans for impossible envelopes.
            from .geometry.core import sections
            for i in range(args.count):
                q=rng.random(9);candidate=with_parameters(c,physical(c,q))
                try:
                    s=sections(candidate);points=s["points"]
                    diameter=2*float(np.linalg.norm(points[:,:,:2],axis=2).max())
                    codes=[]
                    if diameter>c["design"]["diameter_max"]:codes.append("DIAMETER")
                    for root in (points[0],points[-1]):
                        if np.linalg.norm(root[:,:2],axis=1).max()>=c["design"]["hub_radius"] or np.abs(root[:,2]).max()>=c["design"]["hub_length"]/2:codes.append("ROOT_CAP_EXPOSED")
                    rows.append({"parameters":candidate["parameters"],"screen_pass":not codes,"codes":sorted(set(codes)),"diameter":diameter})
                except ValueError as exc:rows.append({"parameters":candidate["parameters"],"screen_pass":False,"codes":[str(exc)]})
            summary={"kind":"geometry-envelope-screen-only","seed":args.seed,"count":args.count,"passed":sum(r["screen_pass"] for r in rows),"note":"Passing this cheap screen does not establish valid solids or aerodynamic performance.","rows":rows}
            atomic_json(args.output,summary);emit({k:v for k,v in summary.items() if k!="rows"});return 0
        if args.command in ("geometry","case"):
            g=generate_rotor(c);report=validate_geometry(g,c,deep=not getattr(args,"fast",False))
            if args.command=="case":
                if not report.valid:emit(report.to_dict());return 2
                build_case(g,c,args.output);emit({"case":str(Path(args.output).resolve()),"validated_geometry":True});return 0
            out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
            atomic_json(out/"configuration.json",c);atomic_json(out/"validation.json",report.to_dict());atomic_json(out/"geometry.json",g.measurements)
            g.mesh.export(out/"rotor.stl");g.mesh.export(out/"rotor.ply")
            np.savez_compressed(out/"sections.npz",**g.samples)
            path=write_preview(g,report,c,out/"preview.html")
            emit({"preview":str(path.resolve()),"geometry":g.measurements,"validation":report.to_dict()});return 0 if report.valid else 2
        if args.command=="evaluate":
            result=evaluate_cfd(c,args.output,not args.no_cache,args.initialize_from);emit(result);return 0 if result["status"]=="SUCCESS" else 2
        if args.command=="continue-evaluation":
            from .simulation.continuation import continue_evaluation
            result=continue_evaluation(args.parent,c,args.output);emit(result);return 0 if result["status"]=="SUCCESS" else 2
        if args.command=="compare":
            from .analysis import compare_results
            result=compare_results(args.results);atomic_json(args.output,result);emit(result);return 0
        if args.command=="assess-baseline":
            from .analysis import assess_baseline
            result=assess_baseline(args.result,args.window);atomic_json(args.output,result);emit(result);return 0 if result['status']=='READY_FOR_SENSITIVITY' else 2
        if args.command=="study-report":
            from .study import write_study
            output,result=write_study(args.results,args.output,args.window,args.relative_limit);emit({'report':str(output),'status':result['status'],'qualified':False});return 0 if result['status']=='READY_FOR_EXPLORATORY_SEARCH' else 2
        if args.command=="surface-audit":
            from .simulation.surface_audit import surface_audit
            emit(surface_audit(args.attempt,args.output));return 0
        if args.command=="freeze-baseline":
            from .baseline import freeze_baseline
            emit(freeze_baseline(c,args.result,args.output,args.study,args.calibrate_timeout));return 0
        if args.command in ("optimize","resume"):
            if not args.exploratory:
                raise ValueError("CFD qualification has not been established. Use --exploratory only for an explicitly unconfirmed pilot after inspecting baseline convergence.")
            initial=None
            if args.command=='optimize' and not args.fresh_baseline:
                reference=Path(args.baseline_reference or (args.config+'.baseline.json'))
                if reference.exists():
                    from .baseline import load_baseline_seed
                    initial=load_baseline_seed(c,reference)
                elif args.baseline_reference:raise FileNotFoundError(reference)
            result=optimize(c,lambda candidate:evaluate_cfd(candidate,args.output),args.checkpoint,resume=args.command=="resume",initial_result=initial)
            emit({"best":result["best"],"proposals":result["proposals"],"evaluations":result["evaluations"],"qualified":False,"checkpoint":args.checkpoint});return 0
        if args.command=="report":
            if Path(args.output).suffix.lower()==".html":
                from .reporting import write_report
                emit({"report":str(write_report(args.result,args.output))});return 0
            result=json.loads(Path(args.result).read_text())
            text="# Evaluation report\n\nStatus: **"+result.get("status","OPTIMIZATION_CHECKPOINT")+"**\n\n"
            text+="This is a simulation research result. It is not a mesh-converged or experimentally validated performance claim.\n\n"
            for name in ("metrics","geometry","diagnostics","reasons","best"):
                if name in result:text+="## "+name.title()+"\n\n```json\n"+json.dumps(result[name],indent=2)+"\n```\n\n"
            path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text);emit({"report":str(path.resolve())});return 0
    except (ValueError,FileNotFoundError,FileExistsError,KeyError) as exc:
        parser.exit(2,"error: "+str(exc)+"\n")


if __name__=="__main__":sys.exit(main())

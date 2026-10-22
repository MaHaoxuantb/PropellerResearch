import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from propeller.config import load_config
from propeller.experiments import evaluation_key,Store,atomic_json
from propeller.optimization.annealing import optimize,reflect,score
from propeller.simulation.parser import read_vectors,summarize
from propeller.simulation.runner import run_command
from propeller.simulation.case_builder import build_case
from propeller.scheduler import check_budget
from propeller.simulation.continuation import restart_inputs,continue_evaluation
from propeller.simulation.initialization import initialization_inputs,mapping_command,check_geometry
from propeller.analysis import assess_baseline
from propeller.simulation.field_io import normalize_ascii_fields
from propeller.study import study_summary,_plot_history
from propeller.baseline import freeze_baseline,load_baseline_seed
from propeller.geometry import generate_rotor
import shutil
import hashlib
import sys

ROOT=Path(__file__).resolve().parents[1]


class WorkflowTests(unittest.TestCase):
    def setUp(self):self.c=load_config(ROOT/'config/demo.json')

    def test_meshing_seed_is_inside_a_farfield_cell_for_even_and_odd_grids(self):
        import re
        from types import SimpleNamespace
        import trimesh
        geometry=SimpleNamespace(mesh=trimesh.creation.box(extents=[.1,.1,.02]),measurements={'diameter':.15})
        for counts,domain in [([20,20,32],(3,3,8)),([25,25,40],(3,3,8)),([30,30,48],(4.5,3,13.5)),([4,4,4],(1.01,1.01,1.01))]:
            with self.subTest(counts=counts,domain=domain),tempfile.TemporaryDirectory() as directory:
                c=copy.deepcopy(self.c)
                c['cfd'].update(background_cells=counts,domain=dict(zip(('radial_diameters','upstream_diameters','downstream_diameters'),domain)))
                case=build_case(geometry,c,directory)
                text=(case/'system/snappyHexMeshDict').read_text()
                point=np.fromstring(re.search(r'locationInMesh\s*\(([^)]+)\)',text).group(1),sep=' ')
                D=c['design']['diameter_max'];radial,upstream,downstream=domain
                lower=np.array([-radial*D,-radial*D,-downstream*D])
                upper=np.array([radial*D,radial*D,upstream*D])
                indices=(point-lower)/(upper-lower)*counts
                self.assertTrue(np.all(point>lower) and np.all(point<upper))
                np.testing.assert_allclose(indices%1,.5,atol=1e-10)
                self.assertGreater(np.linalg.norm(point[:2]),D/2)

    def test_cache_identity_includes_physics_and_fidelity(self):
        key=evaluation_key(self.c)[0]
        for section,name,value in [('operating','rpm',5001),('cfd','fidelity','fine'),('geometry','span_stations',100),('parameters','H',.0762000000001)]:
            c=copy.deepcopy(self.c);c[section][name]=value
            self.assertNotEqual(key,evaluation_key(c)[0])

    def test_failures_are_not_success_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            s=Store(directory);identity,path=s.start('key',{})
            s.finish(identity,path,{'status':'MESH_FAILED'})
            self.assertIsNone(s.cached('key'))
            identity,path=s.start('key',{});s.finish(identity,path,{'status':'SUCCESS','metrics':{'thrust':1,'shaft_power':10}})
            self.assertEqual(s.cached('key')['metrics']['thrust'],1);s.close()

    def test_native_command_timeout_retains_log(self):
        with tempfile.TemporaryDirectory() as d:
            c=copy.deepcopy(self.c)
            c['cfd'].update(launcher=sys.executable,timeout_seconds=.15)
            result=run_command(c,d,['-c','import time; time.sleep(30)'])
            self.assertEqual(result['status'],'TIMEOUT')
            self.assertTrue(Path(result['log']).exists())
            self.assertLess(result['seconds'],5)

    def test_reflection_handles_large_jumps(self):
        np.testing.assert_allclose(reflect(np.array([-2.2,-.2,0,.5,1,1.2,4.2])),[.2,.2,0,.5,1,.8,.2])

    def test_parallel_budget_is_enforced(self):
        with self.assertRaises(ValueError):check_budget([self.c],2,cpu_count=8)
        with self.assertRaises(ValueError):check_budget([self.c],4,4,2,cpu_count=8)
        with self.assertRaises(ValueError):check_budget([self.c],4,16,2,cpu_count=2)
        self.assertEqual(check_budget([self.c],2,8,2,cpu_count=8)['workers'],2)

    def test_recover_only_dead_attempts(self):
        with tempfile.TemporaryDirectory() as directory:
            store=Store(directory);identity,path=store.start('key',{})
            with patch('propeller.experiments.os.kill',side_effect=PermissionError):self.assertEqual(store.recover_interrupted(),[])
            with patch('propeller.experiments.os.kill',side_effect=ProcessLookupError):self.assertEqual(store.recover_interrupted(),[identity])
            self.assertEqual(json.loads((path/'result.json').read_text())['status'],'INTERRUPTED')
            self.assertIsNone(store.cached('key'));store.close()

    def test_invalid_results_have_no_objective(self):
        self.c['operating'].update(target_thrust=1,reference_power=10)
        self.assertIsNone(score({'status':'MESH_FAILED'},self.c))
        self.assertIsNone(score({'status':'SUCCESS','metrics':{'thrust':2,'shaft_power':-1}},self.c))
        s=score({'status':'SUCCESS','metrics':{'thrust':.5,'shaft_power':1}},self.c)
        self.assertFalse(s['feasible'])
        self.assertGreater(s['score'],.1)

    def test_checkpoint_resume_reproduces_uninterrupted_chain(self):
        self.c['operating'].update(target_thrust=1,reference_power=10)
        self.c['optimization'].update(max_proposals=20,max_evaluations=30)
        def synthetic(c):
            # Test fixture only; not a propeller model or CFD claim.
            p=c['parameters'];return {'status':'SUCCESS','metrics':{'thrust':1.2,'shaft_power':8+100*(p['H']-.06)**2+4*p['alpha0']**2}}
        with tempfile.TemporaryDirectory() as directory:
            full=optimize(self.c,synthetic,Path(directory)/'full.json')
            optimize(self.c,synthetic,Path(directory)/'split.json',proposal_limit=7)
            split=optimize(self.c,synthetic,Path(directory)/'split.json',resume=True)
            for key in ('history','q','best','rng','temperature','evaluations'):self.assertEqual(full[key],split[key])

    def test_parser_torque_sign_and_stability(self):
        self.c['cfd'].update(force_window=4,iterations=10)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);folder=root/'postProcessing/forces/0';folder.mkdir(parents=True)
            for name,z in [('force',2),('moment',-.1)]:
                (folder/(name+'.dat')).write_text('# time total pressure viscous\n'+''.join(f'{i} (0 0 {z}) (0 0 {z}) (0 0 0)\n' for i in range(1,11)))
            (root/'log.simpleFoam').write_text(''.join(f'Time = {i}\n'+'\n'.join('Solving for '+name+', Initial residual = 1e-8, Final residual = 1e-9' for name in ('p','Ux','Uy','Uz','k','omega'))+'\n' for i in range(1,11))+'sum local = 1e-9\npatch rotor y+ : min = 1, max = 40, average = 20\nEnd\n')
            result=summarize(root,self.c,.2)
            self.assertEqual(result['status'],'SUCCESS',result)
            self.assertAlmostEqual(result['metrics']['shaft_power'],.1*5000*np.pi/30)
            self.assertIsNone(result['metrics']['efficiency'])
            log=root/'log.simpleFoam';log.write_text(log.read_text().replace('Time = 8\nSolving for p, Initial residual = 1e-8','Time = 8\nSolving for p, Initial residual = 1e-3'))
            self.assertEqual(summarize(root,self.c,.2)['status'],'NOT_CONVERGED')
            (root/'log.simpleFoam').write_text('FOAM FATAL ERROR')
            self.assertEqual(summarize(root,self.c,.2)['status'],'NOT_CONVERGED')

    def test_longer_baseline_window_rejects_slow_settling(self):
        # Synthetic history: stable final 100 samples can conceal earlier drift.
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);folder=root/'case/postProcessing/forces/0';folder.mkdir(parents=True)
            atomic_json(root/'manifest.json',{'config':self.c});atomic_json(root/'result.json',{'status':'SUCCESS'})
            log=root/'case/log.simpleFoam'
            log.write_text(''.join(f'Time = {i}\n'+''.join('Solving for '+field+', Initial residual = 1e-8\n' for field in ('p','Ux','Uy','Uz','k','omega')) for i in range(1,501)))
            (folder/'moment.dat').write_text(''.join(f'{i} (0 0 -.1) (0 0 -.1) (0 0 0)\n' for i in range(1,501)))
            def thrust(drift):
                (folder/'force.dat').write_text(''.join(f'{i} (0 0 {1 if drift and i<400 else 1.5}) (0 0 1) (0 0 0)\n' for i in range(1,501)))
            thrust(True)
            self.assertEqual(assess_baseline(root/'result.json')['status'],'BASELINE_REJECTED')
            thrust(False)
            self.assertEqual(assess_baseline(root/'result.json')['status'],'READY_FOR_SENSITIVITY')
            log.write_text(log.read_text().replace('Initial residual = 1e-8','Initial residual = 1e-3',1))
            self.assertEqual(assess_baseline(root/'result.json')['status'],'BASELINE_REJECTED')
            self.assertEqual(assess_baseline(root/'result.json',100)['status'],'READY_FOR_SENSITIVITY')
            for name in ('force','moment'):
                p=folder/(name+'.dat');p.write_text(p.read_text().splitlines()[0]+'\n')
            self.assertEqual(assess_baseline(root/'result.json')['status'],'BASELINE_REJECTED')

    def test_sensitivity_initialization_checks_actual_geometry(self):
        import trimesh
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);surface=root/'case/constant/triSurface';surface.mkdir(parents=True)
            mesh=trimesh.creation.box(extents=[.1,.1,.1]);mesh.export(surface/'rotor.stl')
            geometry=SimpleNamespace(mesh=mesh,measurements={'diameter':.15})
            self.assertLess(check_geometry(geometry,{'parent':directory})['relative_volume_error'],1e-6)
            geometry.mesh=mesh.copy();geometry.mesh.apply_translation([.001,0,0])
            with self.assertRaisesRegex(ValueError,'differs'):check_geometry(geometry,{'parent':directory})

    def test_ascii_subnormal_repair_is_limited_and_audited(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);p=root/'p'
            p.write_text('FoamFile { format ascii; class volScalarField; }\n(-1e-318 1e-309 1e-300 2e-7 0 0e0)')
            binary=root/'U';binary.write_bytes(b'FoamFile { format binary; class volVectorField; }\n\xff\x00')
            before=binary.read_bytes();report=normalize_ascii_fields(root)
            self.assertEqual(report['changes'][0]['count'],2)
            self.assertIn('(0 0 1e-300 2e-7 0 0e0)',p.read_text())
            self.assertEqual(before,binary.read_bytes())
            self.assertEqual(normalize_ascii_fields(root)['changes'],[])

    def test_study_gate_requires_matched_converged_distinct_cases(self):
        # Synthetic study fixtures exercise the gate, not a CFD prediction.
        with tempfile.TemporaryDirectory() as directory:
            paths=[]
            configs=[copy.deepcopy(self.c) for _ in range(5)]
            configs[1]['cfd']['background_cells']=[25,25,40]
            configs[2]['cfd'].update(background_cells=[30,30,48],domain={'radial_diameters':4.5,'upstream_diameters':3,'downstream_diameters':13.5})
            configs[3]['cfd']['iterations']=2000
            configs[4]['cfd']['background_cells']=[30,30,48]
            for index,c in enumerate(configs):
                root=Path(directory)/str(index);folder=root/'case/postProcessing/forces/0';folder.mkdir(parents=True)
                atomic_json(root/'manifest.json',{'config':c,'runtime':{'binary_sha256':{'solver':'fixture'}}})
                for name in ('system/fvSchemes','system/fvSolution','constant/transportProperties','constant/turbulenceProperties','constant/MRFProperties'):
                    target=root/'case'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_text('same fixture equations')
                for name,value in [('force',1),('moment',-.1)]:
                    (folder/(name+'.dat')).write_text(''.join(f'{i} (0 0 {value}) (0 0 {value}) (0 0 0)\n' for i in range(1,501)))
                (root/'case/log.simpleFoam').write_text(''.join(f'Time = {i}\n'+''.join('Solving for '+field+', Initial residual = 1e-8\n' for field in ('p','Ux','Uy','Uz','k','omega')) for i in range(1,501)))
                path=root/'result.json';atomic_json(path,{'status':'SUCCESS','metrics':{'thrust':1,'shaft_power':50}});paths.append(path)
            study=study_summary(paths)
            self.assertEqual(study['status'],'READY_FOR_EXPLORATORY_SEARCH')
            self.assertEqual(study['comparison_reference'],str(paths[3].resolve()))
            with self.assertRaisesRegex(ValueError,'distinct'):study_summary([paths[0]]*3)
            atomic_json(paths[1],{'status':'NOT_CONVERGED','metrics':{'thrust':1,'shaft_power':50}})
            incomplete=study_summary(paths)
            self.assertEqual(incomplete['status'],'VALIDATION_INCOMPLETE')
            self.assertIsNone(incomplete['cases'][1]['relative_change'])
            atomic_json(paths[1],{'status':'SUCCESS','metrics':{'thrust':1.2,'shaft_power':50}})
            (paths[1].parent/'case/postProcessing/forces/0/force.dat').write_text(''.join(f'{i} (0 0 1.2) (0 0 1.2) (0 0 0)\n' for i in range(1,501)))
            self.assertEqual(study_summary(paths)['status'],'SENSITIVITY_TOO_LARGE')
            (paths[4].parent/'case/postProcessing/forces/0/force.dat').write_text(''.join(f'{i} (0 0 1.21) (0 0 1.21) (0 0 0)\n' for i in range(1,501)))
            selected=study_summary(paths)
            self.assertEqual(selected['status'],'READY_FOR_EXPLORATORY_SEARCH')
            self.assertEqual(selected['recommended_search_source'],str(paths[1].resolve()))
            manifest=paths[2].parent/'manifest.json';data=json.loads(manifest.read_text());data['config']['operating']['rpm']+=1;atomic_json(manifest,data)
            with self.assertRaisesRegex(ValueError,'differ'):study_summary(paths)

    def test_plot_history_joins_restarts_without_duplicate_iterations(self):
        with tempfile.TemporaryDirectory() as directory:
            parent=Path(directory)/'parent';child=Path(directory)/'child'
            for root,times,thrust in [(parent,range(1,5),1),(child,range(4,7),2)]:
                folder=root/'case/postProcessing/forces/0';folder.mkdir(parents=True)
                for name,value in [('force',thrust),('moment',-.1)]:
                    (folder/(name+'.dat')).write_text(''.join(f'{i} (0 0 {value}) (0 0 {value}) (0 0 0)\n' for i in times))
                (root/'case/log.simpleFoam').write_text(''.join(f'Time = {i}\nSolving for Ux, Initial residual = 1e-8\n' for i in times))
            atomic_json(parent/'result.json',{'initialization':{'parent':'not-a-continuation'}})
            atomic_json(child/'result.json',{'continuation':{'parent':str(parent)}})
            force,moment,residuals,segments=_plot_history(child/'result.json')
            self.assertEqual(force[:,0].tolist(),[1,2,3,4,5,6])
            self.assertEqual(force[:,3].tolist(),[1,1,1,2,2,2])
            self.assertEqual(len(segments),2)
            self.assertEqual(sorted(residuals),[1,2,3,4,5,6])

    def test_frozen_baseline_is_reused_without_a_solver_call(self):
        # Synthetic CFD files verify reference integrity and optimizer reuse.
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)/'baseline';case=root/'case';folder=case/'postProcessing/forces/0';folder.mkdir(parents=True)
            runtime={'available':True,'binary_sha256':{'solver':'fixture'}}
            atomic_json(root/'manifest.json',{'config':self.c,'runtime':runtime})
            atomic_json(root/'validation.json',{'valid':True,'checks':{'self_intersection_scan':'completed'}})
            geometry=generate_rotor(self.c);surface=case/'constant/triSurface';surface.mkdir(parents=True);geometry.mesh.export(surface/'rotor.stl')
            snapshot=root/'source/propeller/simulation';snapshot.mkdir(parents=True)
            shutil.copyfile(ROOT/'src/propeller/simulation/case_builder.py',snapshot/'case_builder.py')
            for name,value in [('force',2),('moment',-.1)]:
                (folder/(name+'.dat')).write_text(''.join(f'{i} (0 0 {value}) (0 0 {value}) (0 0 0)\n' for i in range(1,501)))
            (case/'log.simpleFoam').write_text(''.join(f'Time = {i}\n'+''.join('Solving for '+field+', Initial residual = 1e-8\n' for field in ('p','Ux','Uy','Uz','k','omega')) for i in range(1,501))+'sum local = 1e-9\npatch rotor y+ : min = 1, max = 40, average = 20\nEnd\n')
            parsed=summarize(case,self.c,geometry.measurements['diameter']);parsed.update(geometry=geometry.measurements,mesh={'valid':True},attempt_id='synthetic-fixture')
            atomic_json(root/'result.json',parsed)
            config_path=Path(directory)/'study.json';frozen=freeze_baseline(self.c,root/'result.json',config_path)
            config=load_config(config_path);config['optimization'].update(max_proposals=1,max_evaluations=1)
            config['cfd']['timeout_seconds']*=2
            with patch('propeller.baseline.runtime_info',return_value=runtime):seed=load_baseline_seed(config,frozen['reference'])
            self.assertEqual(seed['source_kind'],'accepted-baseline-reuse')
            state=optimize(config,lambda _:self.fail('Baseline was evaluated again'),Path(directory)/'checkpoint.json',proposal_limit=0,initial_result=seed)
            self.assertEqual(state['evaluations'],0)
            self.assertTrue(state['best']['feasible'])
            bad=copy.deepcopy(config);bad['operating']['rpm']+=1
            with self.assertRaisesRegex(ValueError,'contract differs'):load_baseline_seed(bad,frozen['reference'])

    def test_parser_rejects_missing_and_nonmonotone_rows(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'force.dat';path.write_text('0 1 2\n')
            with self.assertRaises(ValueError):read_vectors(path)
            path.write_text('1 0 0 1 0 0 1 0 0 0\n1 0 0 1 0 0 1 0 0 0\n')
            with self.assertRaises(ValueError):read_vectors(path)

    def test_continuation_preserves_parent_and_rejects_physics_change(self):
        # Synthetic file fixture; no numerical CFD result is asserted here.
        with tempfile.TemporaryDirectory() as directory:
            parent=Path(directory)/'parent';case=parent/'case'
            for folder in ('constant/triSurface','system','0','1000'):(case/folder).mkdir(parents=True,exist_ok=True)
            for field in ('U','p','k','omega','nut'):(case/'1000'/field).write_text('restart field fixture')
            (case/'constant/triSurface/rotor.stl').write_text('surface fixture')
            (case/'system/controlDict').write_text('startFrom startTime; endTime 1000;')
            for name in ('geometry.json','case.foam','log.checkMesh','log.snappyHexMesh'):(case/name).write_text('fixture')
            atomic_json(case/'configuration.json',self.c)
            runtime={'available':True,'binary_sha256':{'solver':'fixture'}}
            atomic_json(parent/'manifest.json',{'config':self.c,'runtime':runtime})
            atomic_json(parent/'result.json',{'status':'NOT_CONVERGED','geometry':{'diameter':.2},'mesh':{'valid':True},'geometry_sha256':hashlib.sha256(b'surface fixture').hexdigest()})
            c=copy.deepcopy(self.c);c['cfd']['iterations']=2000
            strict=copy.deepcopy(c);strict['cfd']['force_window']=200
            self.assertEqual(restart_inputs(parent,strict)[3],'1000')
            strict['cfd']['force_window']=50
            with self.assertRaisesRegex(ValueError,'shorten'):restart_inputs(parent,strict)
            bad=copy.deepcopy(c);bad['operating']['rpm']+=1
            with self.assertRaisesRegex(ValueError,'only iterations'):restart_inputs(parent,bad)
            with patch('propeller.simulation.continuation.runtime_info',return_value=runtime),patch('propeller.simulation.continuation.run_command',return_value={'status':'SUCCESS'}),patch('propeller.simulation.continuation.summarize',return_value={'status':'NOT_CONVERGED'}):
                result=continue_evaluation(parent,c,Path(directory)/'runs')
            self.assertEqual(result['status'],'NOT_CONVERGED')
            self.assertEqual((case/'system/controlDict').read_text(),'startFrom startTime; endTime 1000;')
            new=Path(result['path'])/'case'
            self.assertIn('startFrom latestTime;', (new/'system/controlDict').read_text())
            self.assertEqual(result['continuation']['restart_time'],'1000')
            self.assertTrue((new/'1000/U').exists())
            with self.assertRaisesRegex(ValueError,'converged source'):initialization_inputs(parent,c,runtime)
            source=json.loads((parent/'result.json').read_text());source['status']='SUCCESS';atomic_json(parent/'result.json',source)
            c['cfd']['background_cells']=[25,25,40]
            initialization=initialization_inputs(parent,c,runtime)
            self.assertIn('-consistent',mapping_command(new,initialization))
            bad=copy.deepcopy(c);bad['parameters']['H']+=.001
            with self.assertRaisesRegex(ValueError,'identical parameters'):initialization_inputs(parent,bad,runtime)
            c['cfd']['domain']={'radial_diameters':4.5,'upstream_diameters':3,'downstream_diameters':13.5}
            initialization=initialization_inputs(parent,c,runtime)
            self.assertNotIn('-consistent',mapping_command(new,initialization))


if __name__=='__main__':unittest.main()

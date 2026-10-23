"""Synthetic fixtures stay in temporary directories, never the CFD archive."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np
from propeller.config import load_config
from propeller.geometry import generate_rotor
from propeller.geometry.screening import fixed_radius_config,free_coordinates,section_screen
from propeller.optimization.bayesian import pareto_indices,hypervolume,accepted_metrics,propose
from propeller.optimization.campaign import ingest,run,write_status,worker_limit
from propeller.simulation.evaluator import evaluate_cfd
from propeller.simulation.runner import run_command
from propeller.experiments import source_fingerprint,canonical,atomic_json
import hashlib

ROOT=Path(__file__).resolve().parents[1]


class BayesianTests(unittest.TestCase):
    def setUp(self):self.c=load_config(ROOT/'config/flow-oriented-fine-4500.json')

    def test_radius_mapping_preserves_finished_solid_and_hub(self):
        radius=.11175999790430069
        q=np.asarray(free_coordinates(self.c));q[2]+=.02
        c=fixed_radius_config(self.c,q,radius)
        g=generate_rotor(c)
        self.assertLess(abs(g.measurements['diameter']/2-radius),5e-5)
        self.assertEqual(c['design'],self.c['design'])
        self.assertNotEqual(c['parameters']['H'],self.c['parameters']['H'])
        with self.assertRaises(ValueError):fixed_radius_config(self.c,q,.021)

    def test_bad_geometry_never_launches_native_mesh_or_solver(self):
        self.c['design']['diameter_max']=.15
        self.assertFalse(section_screen(self.c)['valid'])
        with tempfile.TemporaryDirectory() as tmp,patch('propeller.simulation.evaluator.runtime_info',return_value={'available':True}),patch('propeller.simulation.evaluator.run_command') as native:
            result=evaluate_cfd(self.c,tmp)
            self.assertEqual(result['status'],'GEOMETRY_INVALID');native.assert_not_called()
            self.assertTrue((Path(result['path'])/'screening.json').exists())

    def test_pareto_frontier_and_area_exclude_dominated_points(self):
        points=[[1,1],[2,.5],[.5,.5]]
        self.assertEqual(pareto_indices(points),[0,1]);self.assertAlmostEqual(hypervolume(points),1.5)
        self.assertIsNone(accepted_metrics({'status':'NOT_CONVERGED','metrics':{'thrust':2,'shaft_power':1}}))
        self.assertIsNone(accepted_metrics({'status':'SUCCESS','metrics':{'thrust':2,'shaft_power':0}}))
        self.assertEqual(accepted_metrics({'status':'SUCCESS','metrics':{'thrust':-1,'shaft_power':2}}),[-1,2])
        self.assertEqual(pareto_indices([[1,1],[-1,-.5]]),[0])

    def test_initial_trials_start_near_baseline(self):
        obs=[{'q':[.5]*8,'metrics':[1,2]}]
        pool=np.array([[.9]*8,[.6]*8,[.55]*8])
        index,method=propose(pool,obs,[1,2],[],7,False)
        self.assertEqual(index,1);self.assertEqual(method,'space-filling-local')

    def test_failed_case_is_not_a_fake_objective(self):
        state={'observations':[{'metrics':[1,1]}],'baseline_metrics':[1,1],'evaluations':0,'geometry_rejections':0,
               'failure_streak':0,'stagnation':0,'best_hypervolume':1,'settings':{'improvement_tolerance':.005}}
        job={'id':'test','q':[.5]*8,'parameters':{}}
        ingest(state,job,{'status':'NOT_CONVERGED','metrics':{'thrust':100,'shaft_power':.001}})
        self.assertEqual(len(state['observations']),1);self.assertEqual(state['failure_streak'],1)
        ingest(state,job,{'status':'GEOMETRY_INVALID'})
        self.assertEqual(state['evaluations'],1);self.assertEqual(state['geometry_rejections'],1)
        ingest(state,job,{'status':'SUCCESS','metrics':{'thrust':2,'shaft_power':1},'path':'fixture'})
        self.assertEqual(len(state['observations']),2);self.assertEqual(state['failure_streak'],0)
        self.assertEqual(state['best_hypervolume'],4)

    def test_worker_ramp_respects_observed_memory(self):
        state={'settings':{'workers':4,'job_memory_gb':3,'memory_budget_gb':12},'jobs':[],'observations':[{}]}
        self.assertEqual(worker_limit(state),2)
        state['observations']=[{},{},{}];self.assertEqual(worker_limit(state),4)
        state['jobs']=[{'peak_child_gb':4}];self.assertEqual(worker_limit(state),2)
        state['jobs']=[{'peak_child_gb':10}];self.assertEqual(worker_limit(state),0)

    def test_expired_deadline_launches_no_process(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict('os.environ',{'PROPELLER_DEADLINE_EPOCH':str(time.time()-1)}),patch('subprocess.Popen') as spawn:
            r=run_command(self.c,tmp,['simpleFoam']);self.assertEqual(r['status'],'TIMEOUT');spawn.assert_not_called()

    def test_deadline_terminates_owned_command(self):
        c=copy.deepcopy(self.c);c['cfd'].update(launcher=sys.executable,timeout_seconds=20)
        with tempfile.TemporaryDirectory() as tmp,patch.dict('os.environ',{'PROPELLER_DEADLINE_EPOCH':str(time.time()+.3),'PROPELLER_SUPERVISOR_PID':'0'}):
            r=run_command(c,tmp,['-c','import time; time.sleep(20)'])
            self.assertEqual(r['status'],'TIMEOUT');self.assertLess(r['seconds'],3)

    def test_completed_orphan_recovery_refuses_active_child(self):
        from propeller.simulation.recovery import recover_completed
        with tempfile.TemporaryDirectory() as tmp:
            atomic_json(Path(tmp)/'state.json',{'status':'RUNNING','pid':123})
            with patch('propeller.simulation.recovery.process_alive',return_value=False),patch('propeller.simulation.recovery.live_case_process',return_value=True):
                with self.assertRaisesRegex(ValueError,'still alive'):recover_completed(tmp)

    def test_campaign_stop_is_durable_and_launches_nothing(self):
        settings=json.loads((ROOT/'config/mobo-pilot.json').read_text())
        obs={'id':'baseline','q':[.5]*8,'parameters':{},'metrics':[1,1],'path':'fixture'}
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'jobs').mkdir();atomic_json(root/'base.json',self.c)
            state={'status':'PREPARED','source':source_fingerprint(),'base_sha256':hashlib.sha256(canonical(self.c).encode()).hexdigest(),
                   'settings':settings,'observations':[obs],'baseline_metrics':[1,1],'jobs':[],'evaluations':0,'geometry_rejections':0,
                   'failure_streak':0,'stagnation':0,'best_hypervolume':1.,'deadline':time.time()+100,'rng':np.random.default_rng(42).bit_generator.state}
            write_status(root,state);(root/'STOP').touch()
            with patch('subprocess.Popen') as spawn:run(root);spawn.assert_not_called()
            final=json.loads((root/'state.json').read_text());self.assertEqual(final['status'],'STOPPED')
            self.assertEqual(final['stop_reason'],'user stop requested')

    def test_real_bayesian_acquisition_smoke(self):
        try:import botorch
        except ImportError:self.skipTest('Run in .venv-bo for the real acquisition test')
        rng=np.random.default_rng(7)
        observations=[{'q':q.tolist(),'metrics':[float(1+q[0]),float(2+q[1])]} for q in rng.random((6,8))]
        pool=rng.random((12,8))
        index,method=propose(pool,observations,[1,2],[[.5]*8],seed=7)
        self.assertEqual(method,'qLogNEHVI');self.assertIn(index,range(len(pool)))


if __name__=='__main__':unittest.main()

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
from propeller.scheduler import check_budget
import sys

ROOT=Path(__file__).resolve().parents[1]


class WorkflowTests(unittest.TestCase):
    def setUp(self):self.c=load_config(ROOT/'config/demo.json')

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
            (root/'log.simpleFoam').write_text('\n'.join('Solving for '+name+', Initial residual = 1e-8, Final residual = 1e-9' for name in ('p','Ux','Uy','Uz','k','omega'))+'\nsum local = 1e-9\npatch rotor y+ : min = 1, max = 40, average = 20\nEnd\n')
            result=summarize(root,self.c,.2)
            self.assertEqual(result['status'],'SUCCESS',result)
            self.assertAlmostEqual(result['metrics']['shaft_power'],.1*5000*np.pi/30)
            self.assertIsNone(result['metrics']['efficiency'])
            (root/'log.simpleFoam').write_text('FOAM FATAL ERROR')
            self.assertEqual(summarize(root,self.c,.2)['status'],'NOT_CONVERGED')

    def test_parser_rejects_missing_and_nonmonotone_rows(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'force.dat';path.write_text('0 1 2\n')
            with self.assertRaises(ValueError):read_vectors(path)
            path.write_text('1 0 0 1 0 0 1 0 0 0\n1 0 0 1 0 0 1 0 0 0\n')
            with self.assertRaises(ValueError):read_vectors(path)


if __name__=='__main__':unittest.main()

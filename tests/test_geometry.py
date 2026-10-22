import copy
from pathlib import Path
import unittest
import numpy as np
import trimesh
from propeller.config import load_config,normalized,physical,validate_config
from propeller.geometry.core import curve,sections,airfoil,generate_rotor,oriented_airfoils
from propeller.geometry.validation import validate_geometry,self_intersections

ROOT=Path(__file__).resolve().parents[1]


class GeometryMathTests(unittest.TestCase):
    def setUp(self):self.c=load_config(ROOT/'config/demo.json')

    def test_endpoints_apex(self):
        p=self.c['parameters'];rh=self.c['design']['root_radius']
        C,D,DD=curve([0,.5,1],p,rh)
        np.testing.assert_allclose(np.linalg.norm(C[:,:2],axis=1),[rh,rh+p['H'],rh])
        np.testing.assert_allclose(C[:,2],[0,p['Z'],0],atol=1e-14)
        np.testing.assert_allclose(C[0,:2],[rh*np.cos(p['beta']),-rh*np.sin(p['beta'])])

    def test_derivatives(self):
        u=np.linspace(.02,.98,33);h=1e-5;p=self.c['parameters'];rh=self.c['design']['root_radius']
        C,D,DD=curve(u,p,rh);before=curve(u-h,p,rh)[0];after=curve(u+h,p,rh)[0]
        np.testing.assert_allclose(D,(after-before)/(2*h),atol=1e-8)
        np.testing.assert_allclose(DD,(after-2*C+before)/h**2,atol=1e-6)

    def test_frame_and_section_shape(self):
        s=sections(self.c)
        frame=np.stack([s['T'],s['e1'],s['e2']],axis=1)
        np.testing.assert_allclose(frame@frame.transpose(0,2,1),np.broadcast_to(np.eye(3),frame.shape),atol=1e-12)
        np.testing.assert_allclose(np.linalg.det(frame),1,atol=1e-12)
        delta=(s['points']-s['C'][:,None,:])/s['width'][:,None,None]
        x=np.einsum('ijk,ik->ij',delta,s['d'])+s['root_shift'][:,None]
        y=np.einsum('ijk,ik->ij',delta,s['n'])
        np.testing.assert_allclose(x,np.broadcast_to(s['airfoil'][:,0],x.shape),atol=1e-12)
        np.testing.assert_allclose(y,np.broadcast_to(s['airfoil'][:,1],y.shape),atol=1e-12)
        self.assertTrue(np.all(np.diff(s['s'])>0))
        self.assertGreater(np.min(np.diff(s['u'])),1e-8)

    def test_bernstein_convex_hull(self):
        self.c['parameters'].update(b0=.01,b1=.025,b2=.015)
        s=sections(self.c)
        self.assertGreaterEqual(s['width'].min(),.01)
        self.assertLessEqual(s['width'].max(),.025)
        self.assertAlmostEqual(s['width'][0],.01)
        self.assertAlmostEqual(s['width'][-1],.015)

    def test_shape_exponent_singularity_rejected(self):
        self.c['parameters']['shape_exponent']=1.5
        with self.assertRaises(ValueError):sections(self.c)

    def test_normalization_and_nonfinite(self):
        q=normalized(self.c);p=physical(self.c,q)
        for k,v in p.items():self.assertAlmostEqual(v,self.c['parameters'][k])
        self.c['parameters']['H']=float('nan')
        with self.assertRaises(ValueError):validate_config(self.c)

    def test_airfoil_thickness_and_trailing_edge(self):
        a=airfoil(201,.12)
        self.assertAlmostEqual(2*a[:,1].max(),.12,delta=.0001)
        self.assertGreater(a[0,1]-a[-1,1],0)

    def test_opposed_profiles_keep_correspondence_and_thickness(self):
        u=np.linspace(0,1,201);profiles,flip=oriented_airfoils(u,25,.12)
        np.testing.assert_allclose(profiles[0,:25,1],profiles[-1,:25,1][::-1],atol=1e-14)
        np.testing.assert_allclose(np.ptp(profiles[:,:,1],axis=1),np.ptp(profiles[0,:,1]),atol=1e-14)
        np.testing.assert_allclose(profiles[:,:,0].min(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(profiles[:,:,0].max(axis=1),1,atol=1e-14)
        self.assertTrue(np.all(np.linalg.norm(np.roll(profiles,-1,axis=1)-profiles,axis=2)>0))
        self.assertTrue(np.all(np.diff(flip)>=0))
        self.assertEqual(flip[0],0);self.assertEqual(flip[-1],1)

    def test_opposed_section_coordinates_follow_the_profile(self):
        c=load_config(ROOT/'config/flow-oriented.json');s=sections(c)
        delta=(s['points']-s['C'][:,None,:])/s['width'][:,None,None]
        x=np.einsum('ijk,ik->ij',delta,s['d'])+s['root_shift'][:,None]
        y=np.einsum('ijk,ik->ij',delta,s['n'])
        np.testing.assert_allclose(x,s['airfoil_profiles'][:,:,0],atol=1e-12)
        np.testing.assert_allclose(y,s['airfoil_profiles'][:,:,1],atol=1e-12)
        active=(s['u']>.15)&(s['u']<.85)&(s['section_speed']>.1)
        self.assertTrue(np.all(s['chord_flow'][active]>=0))

    def test_crossing_triangles_detected(self):
        vertices=[[-1,-1,0],[1,-1,0],[0,1,0],[0,0,-1],[0,0,1],[.5,.5,0]]
        m=trimesh.Trimesh(vertices=vertices,faces=[[0,1,2],[3,4,5]],process=False)
        self.assertTrue(self_intersections(m))
        m.vertices[3:]+=10
        self.assertFalse(self_intersections(m))

    def test_shared_vertex_does_not_hide_a_fold(self):
        m=trimesh.Trimesh(vertices=[[0,0,0],[2,0,0],[0,2,0],[1,1,-1],[1,1,1]],faces=[[0,1,2],[0,3,4]],process=False)
        self.assertTrue(self_intersections(m))

    def test_shared_edge_coplanar_overlap_is_invalid(self):
        m=trimesh.Trimesh(vertices=[[0,0,0],[2,0,0],[0,2,0],[1,1,0]],faces=[[0,1,2],[0,1,3]],process=False)
        self.assertTrue(self_intersections(m))


class SolidTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=load_config(ROOT/'config/demo.json');cls.g=generate_rotor(cls.c)

    def test_connected_watertight_positive_solid(self):
        r=validate_geometry(self.g,self.c,deep=False)
        self.assertTrue(r.valid,r.violations)
        self.assertEqual(len(self.g.mesh.split()),1)
        self.assertGreater(self.g.mesh.volume,0)

    def test_opposed_profile_forms_one_valid_solid(self):
        c=load_config(ROOT/'config/flow-oriented.json');g=generate_rotor(c);r=validate_geometry(g,c,deep=False)
        self.assertTrue(r.valid,r.violations)
        self.assertEqual(g.measurements['reverse_flow_station_fraction'],0)
        self.assertEqual(g.measurements['geometric_reverse_flow_station_fraction'],.5)

    def test_rotation_and_determinism(self):
        a,b=self.g.blades
        np.testing.assert_allclose(a.vertices[:,:2],-b.vertices[:,:2],atol=1e-12)
        np.testing.assert_allclose(a.vertices[:,2],b.vertices[:,2],atol=1e-12)
        g2=generate_rotor(self.c)
        np.testing.assert_array_equal(self.g.mesh.vertices,g2.mesh.vertices)
        np.testing.assert_array_equal(self.g.mesh.faces,g2.mesh.faces)

    def test_diameter_violation(self):
        c=copy.deepcopy(self.c);c['design']['diameter_max']=.1
        r=validate_geometry(self.g,c,deep=False)
        self.assertIn('DIAMETER',[v['code'] for v in r.violations])

    def test_open_surface_rejected(self):
        g=copy.deepcopy(self.g);g.mesh.update_faces(np.arange(len(g.mesh.faces)-1))
        r=validate_geometry(g,self.c,deep=False)
        self.assertIn('OPEN_SURFACE',[v['code'] for v in r.violations])

    def test_colliding_blades_rejected(self):
        g=copy.deepcopy(self.g);g.blades[1]=g.blades[0].copy()
        report=validate_geometry(g,self.c,deep=False)
        self.assertIn('BLADE_COLLISION',[v['code'] for v in report.violations])

    def test_clearance_constraint(self):
        c=copy.deepcopy(self.c);c['geometry']['minimum_clearance']=.1
        report=validate_geometry(self.g,c,deep=False)
        self.assertIn('BLADE_CLEARANCE',[v['code'] for v in report.violations])


if __name__=='__main__':unittest.main()

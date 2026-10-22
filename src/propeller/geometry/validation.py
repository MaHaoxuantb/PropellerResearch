"""Staged solid validation with explicit triangle-intersection checks."""
from dataclasses import dataclass, field, asdict
import numpy as np
from shapely.geometry import Polygon
import trimesh
import manifold3d


@dataclass
class ValidationReport:
    valid: bool = True
    violations: list = field(default_factory=list)
    checks: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def reject(self, code, message):
        self.valid=False
        self.violations.append({"code":code,"message":message})

    def to_dict(self): return asdict(self)


def _segments_hit_triangle(segments, triangle, eps):
    edge1,edge2=triangle[1]-triangle[0],triangle[2]-triangle[0]
    ray=segments[:,1]-segments[:,0]
    h=np.cross(ray,edge2)
    det=h@edge1
    usable=np.abs(det)>1e-12*np.linalg.norm(ray,axis=1)*np.linalg.norm(edge1)*np.linalg.norm(edge2)
    if not np.any(usable): return False
    h,ray,det,segments=h[usable],ray[usable],det[usable],segments[usable]
    s=segments[:,0]-triangle[0]
    u=np.einsum("ij,ij->i",s,h)/det
    q=np.cross(s,edge1)
    v=np.einsum("ij,ij->i",ray,q)/det
    t=q@edge2/det
    return bool(np.any((u>=-1e-8)&(v>=-1e-8)&(u+v<=1+1e-8)&(t>=-1e-8)&(t<=1+1e-8)))


def triangles_intersect(a,b,eps):
    na=np.cross(a[1]-a[0],a[2]-a[0]); nb=np.cross(b[1]-b[0],b[2]-b[0])
    la,lb=np.linalg.norm(na),np.linalg.norm(nb)
    if min(la,lb)<eps**2: return False
    if np.linalg.norm(np.cross(na/la,nb/lb))<1e-8 and np.max(np.abs((b-a[0])@(na/la)))<eps:
        axes=[i for i in range(3) if i!=int(np.argmax(np.abs(na)))]
        return not Polygon(a[:,axes]).intersection(Polygon(b[:,axes])).is_empty
    return _segments_hit_triangle(np.stack((a,np.roll(a,-1,axis=0)),axis=1),b,eps) or _segments_hit_triangle(np.stack((b,np.roll(b,-1,axis=0)),axis=1),a,eps)


def self_intersections(mesh,eps=1e-9,limit=8):
    """AABB sweep, allowing shared topological edges but not interior overlap."""
    triangles=mesh.triangles
    lower=triangles.min(axis=1); upper=triangles.max(axis=1)
    axis=int(np.argmax(np.ptp(mesh.vertices,axis=0)))
    order=np.argsort(lower[:,axis]); found=[]
    for position,i in enumerate(order):
        end=np.searchsorted(lower[order,axis],upper[i,axis]+eps,side="right")
        js=order[position+1:end]
        js=js[np.all(upper[js]>=lower[i]-eps,axis=1)&np.all(lower[js]<=upper[i]+eps,axis=1)]
        for j in js:
            shared=np.intersect1d(mesh.faces[i],mesh.faces[j]).size
            a,b=triangles[i],triangles[j]
            if shared:
                na,nb=mesh.face_normals[i],mesh.face_normals[j]
                if shared==2 and abs(float(na@nb))<1-1e-8:continue
                # Remove the legal shared boundary, then test the interiors.
                a=a.mean(axis=0)+(a-a.mean(axis=0))*(1-1e-6)
                b=b.mean(axis=0)+(b-b.mean(axis=0))*(1-1e-6)
            if triangles_intersect(a,b,eps):
                found.append([int(i),int(j)])
                if len(found)>=limit: return found
    return found


def validate_geometry(geometry,config,deep=True):
    r=ValidationReport(warnings=list(geometry.warnings)); mesh=geometry.mesh
    d,g=config["design"],config["geometry"]
    if not np.all(np.isfinite(mesh.vertices)): r.reject("NONFINITE","Non-finite vertices")
    if not mesh.is_watertight: r.reject("OPEN_SURFACE","Surface is not watertight")
    if not mesh.is_winding_consistent: r.reject("WINDING","Inconsistent winding")
    if mesh.volume<=0: r.reject("VOLUME","Non-positive volume")
    components=len(mesh.split(only_watertight=False))
    if components!=1: r.reject("DISCONNECTED","Expected one connected rotor, found %d"%components)
    if geometry.measurements["diameter"]>d["diameter_max"]+1e-7:
        r.reject("DIAMETER","Finished solid exceeds maximum diameter")
    if np.any(mesh.area_faces<1e-16): r.reject("DEGENERATE_FACES","Degenerate triangles")
    r.checks.update(watertight=bool(mesh.is_watertight),components=components,diameter=geometry.measurements["diameter"])
    s=geometry.samples
    for chord_fraction in (.01,.25,.5,.75,.99):
        mid=s["C"]+s["width"][:,None]*(chord_fraction-s["root_shift"])[:,None]*s["d"]
        jacobian=np.sum(np.gradient(mid,s["u"],axis=0)*s["T"],axis=1)
        active=np.linalg.norm(mid[:,:2],axis=1)>d["hub_radius"]
        if np.any(jacobian[active]<=0):r.reject("LOCAL_FOLD","Local swept section folds outside the hub")
    minimum_body=float(s["width"].min()*g["thickness_ratio"])
    minimum_te=float(s["width"].min()*(s["airfoil"][0,1]-s["airfoil"][-1,1]))
    if minimum_body<g["minimum_body_thickness"]: r.reject("BODY_THICKNESS","Nominal section body thickness too small")
    if minimum_te<g["minimum_trailing_edge_thickness"]: r.reject("TRAILING_EDGE","Trailing edge below minimum thickness")
    for endpoint in (0,-1):
        root=s["points"][endpoint]
        if np.max(np.linalg.norm(root[:,:2],axis=1))>=d["hub_radius"] or np.max(np.abs(root[:,2]))>=d["hub_length"]/2:
            r.reject("ROOT_CAP_EXPOSED","Root cap is not contained in the hub")
        root_area=Polygon(s["airfoil"]).area*s["width"][endpoint]**2
        if root_area<g["minimum_root_area"]: r.reject("ROOT_AREA","Root section area too small")
    # Boolean intersection checks distinguish the intended hub overlap.
    outside=[]
    for blade in geometry.blades:
        connected=trimesh.boolean.intersection([blade,geometry.hub],engine="manifold")
        if connected.volume<1e-12: r.reject("ROOT_DETACHED","Blade has no finite overlap with hub")
        outside.append(trimesh.boolean.difference([blade,geometry.hub],engine="manifold"))
    for i in range(len(outside)):
        for j in range(i+1,len(outside)):
            common=trimesh.boolean.intersection([outside[i],outside[j]],engine="manifold")
            if len(common.faces) and common.volume>1e-12: r.reject("BLADE_COLLISION","Blades overlap outside hub")
            if g["minimum_clearance"]>0 and len(outside[i].faces) and len(outside[j].faces):
                def as_manifold(m):
                    return manifold3d.Manifold(manifold3d.Mesh(vert_properties=np.asarray(m.vertices,dtype=np.float32),tri_verts=np.asarray(m.faces,dtype=np.uint32)))
                gap=as_manifold(outside[i]).min_gap(as_manifold(outside[j]),g["minimum_clearance"])
                if gap<g["minimum_clearance"]*(1-1e-5):r.reject("BLADE_CLEARANCE","Blades do not meet configured clearance outside hub")
    # The solid-level check is authoritative; always report whether it ran.
    if deep:
        crossings=self_intersections(mesh,g["intersection_tolerance"])
        if crossings: r.reject("SELF_INTERSECTION","Triangle crossings: "+str(crossings))
        r.checks["self_intersection_scan"]="completed"
    else:
        r.checks["self_intersection_scan"]="not_run"
        r.warnings.append("Fast validation is not sufficient for CFD acceptance.")
    r.checks.update(minimum_nominal_body_thickness=minimum_body,minimum_trailing_edge_thickness=minimum_te)
    return r

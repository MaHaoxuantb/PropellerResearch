"""Analytic reference curve, explicit section loft and exact mesh Booleans.

The only root-specific deformation is a smooth chordwise retraction into the hub.
No voxelization, smoothing or automatic topology repair is performed.
"""
from dataclasses import dataclass, field
import math
import numpy as np
import trimesh


@dataclass
class GeometryResult:
    mesh: trimesh.Trimesh
    blades: list
    hub: trimesh.Trimesh
    samples: dict
    measurements: dict
    warnings: list = field(default_factory=list)


def unit(v):
    length = np.linalg.norm(v, axis=-1, keepdims=True)
    if np.any(length < 1e-14):
        raise ValueError("Degenerate tangent/frame")
    return v / length


def curve(u, p, root_radius):
    u = np.asarray(u, dtype=float)
    sn = np.maximum(np.sin(np.pi*u), 0)
    cs = np.cos(np.pi*u)
    exponent = p["shape_exponent"]
    if exponent != 1 and exponent < 2:
        raise ValueError("Exponent in (1,2) gives singular root curvature")
    r = root_radius + p["H"]*sn**exponent
    rp = p["H"]*exponent*np.pi*sn**(exponent-1)*cs
    if exponent == 1:
        rpp = -p["H"]*np.pi**2*sn
    else:
        rpp = p["H"]*exponent*np.pi**2*((exponent-1)*sn**(exponent-2)*cs**2-sn**exponent)
    theta = p["beta"]*(2*u-1)
    tp = 2*p["beta"]
    co, si = np.cos(theta), np.sin(theta)
    C = np.column_stack((r*co, r*si, p["Z"]*sn**2))
    D = np.column_stack((rp*co-r*si*tp, rp*si+r*co*tp, p["Z"]*np.pi*np.sin(2*np.pi*u)))
    DD = np.column_stack((rpp*co-2*rp*si*tp-r*co*tp**2,
                          rpp*si+2*rp*co*tp-r*si*tp**2,
                          2*p["Z"]*np.pi**2*np.cos(2*np.pi*u)))
    return C, D, DD


def airfoil(points=25, thickness=.12):
    """NACA symmetric 4-digit section, standard finite trailing edge.

    Returned contour is upper TE -> LE -> lower TE, then closure at TE.
    Formula uses half thickness, so there is no second factor of 1/2.
    """
    x = .5*(1-np.cos(np.linspace(0, np.pi, points)))
    yt = 5*thickness*(.2969*np.sqrt(x)-.1260*x-.3516*x**2+.2843*x**3-.1015*x**4)
    return np.column_stack((np.r_[x[::-1], x[1:]], np.r_[yt[::-1], -yt[1:]]))


def oriented_airfoils(u,points,thickness,transition=.015):
    """Mirror the NACA thickness distribution smoothly on the return leg.

    Equal finite end caps retain vertex correspondence as the round nose and
    thin trailing region exchange ends. The nose truncation is the standard
    finite TE half-gap. Only the apex connector morphs between profiles.
    """
    x=.5*(1-np.cos(np.linspace(0,np.pi,points)))
    half=5*thickness*(.2969*np.sqrt(x)-.1260*x-.3516*x*x+.2843*x**3-.1015*x**4)
    half=np.maximum(half,half[-1])
    t=np.clip((u-(.5-transition))/(2*transition),0,1)
    flip=t*t*t*(10+t*(-15+6*t))
    y=(1-flip[:,None])*half[None,:]+flip[:,None]*half[::-1][None,:]
    y*=half.max()/y.max(axis=1)[:,None]
    profiles=np.empty((len(u),2*points,2))
    profiles[:,:,0]=np.r_[x[::-1],x]
    profiles[:,:,1]=np.concatenate((y[:,::-1],-y),axis=1)
    return profiles,flip


def sections(config):
    p, g, d = (config[k] for k in ("parameters", "geometry", "design"))
    # Curvature-weighted sampling retains endpoints/apex and resolves turning roots.
    dense = np.linspace(0, 1, max(2001, 10*g["span_stations"]))
    C0, D0, DD0 = curve(dense, p, d["root_radius"])
    ds = np.linalg.norm(np.diff(C0[:,:2], axis=0), axis=1)
    s0 = np.r_[0, np.cumsum(ds)]; s0 /= s0[-1]
    turn = np.linalg.norm(np.diff(unit(D0), axis=0), axis=1)
    weight = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C0, axis=0),axis=1)/d["diameter_max"] + .08*turn)]
    weight /= weight[-1]
    u = np.r_[np.interp(np.linspace(0,1,g["span_stations"]), weight, dense), .5]
    mapping=g.get('airfoil_mapping','fixed_naca')
    if mapping=='opposed_naca':
        band=g.get('airfoil_transition_fraction',.015)
        u=np.r_[u,np.linspace(.5-band,.5+band,9)]
    u = np.sort(u)
    u = u[np.r_[True,np.diff(u)>1e-10]]
    u[np.argmin(np.abs(u-.5))]=.5
    C, D, DD = curve(u, p, d["root_radius"])
    s = np.interp(u, dense, s0)
    T = unit(D)
    e1 = unit(np.column_stack((D[:,1], -D[:,0], np.zeros(len(u)))))
    e2 = np.cross(T, e1)
    alpha = p["alpha0"] + p["alpha1"]*(2*u-1)
    direction = np.cos(alpha)[:,None]*e1 + np.sin(alpha)[:,None]*e2
    normal = np.cross(T, direction)
    w = p["b0"]*(1-s)**2 + 2*p["b1"]*s*(1-s)+p["b2"]*s*s
    section = airfoil(g["chord_points"], g["thickness_ratio"])
    if mapping=='opposed_naca':
        profiles,flip=oriented_airfoils(u,g['chord_points'],g['thickness_ratio'],g.get('airfoil_transition_fraction',.015))
        section=profiles[0]
    else:
        profiles=np.broadcast_to(section,(len(u),len(section),2));flip=np.zeros(len(u))
    # Retract the root chord into the hub with C2-smooth quintic blending.
    fraction = g["root_transition_fraction"]
    if not 0 < fraction < .25: raise ValueError("root_transition_fraction must lie in (0,.25)")
    t = np.clip(np.minimum(u,1-u)/fraction, 0, 1)
    smooth = t*t*t*(10+t*(-15+6*t))
    shift = .95*(1-smooth)
    points = C[:,None,:] + w[:,None,None]*((profiles[:,:,0]-shift[:,None])[:,:,None]*direction[:,None,:]+profiles[:,:,1,None]*normal[:,None,:])
    omega = config["operating"]["rotation_sign"]*config["operating"]["rpm"]*2*np.pi/60
    flow = -np.cross(np.array([0,0,omega]), C)
    flow[:,2] -= config["operating"]["inflow_speed"]
    section_flow = flow - np.sum(flow*T,axis=1)[:,None]*T
    speed = np.linalg.norm(section_flow,axis=1)
    chord_flow = np.sum(section_flow*direction,axis=1)
    geometric_chord_flow=chord_flow.copy()
    if mapping=='opposed_naca':chord_flow=chord_flow*np.where(flip<.5,1,-1)
    return dict(u=u,s=s,C=C,T=T,e1=e1,e2=e2,d=direction,n=normal,width=w,
                points=points,curvature=np.linalg.norm(np.cross(D,DD),axis=1)/np.linalg.norm(D,axis=1)**3,
                section_flow=section_flow,chord_flow=chord_flow,section_speed=speed,
                root_shift=shift,airfoil=section,airfoil_profiles=profiles,profile_flip_weight=flip,
                geometric_chord_flow=geometric_chord_flow)


def loft(points):
    n, m, _ = points.shape
    vertices = list(points.reshape(-1,3))
    faces=[]
    for i in range(n-1):
        for j in range(m):
            a=i*m+j; b=i*m+(j+1)%m; c=(i+1)*m+(j+1)%m; d=(i+1)*m+j
            faces.extend(((a,b,c),(a,c,d)))
    # Airfoil contour is star-shaped around this interior center.
    for i in (0,n-1):
        center=len(vertices); vertices.append(points[i].mean(axis=0))
        for j in range(m): faces.append((center,i*m+(j+1)%m,i*m+j))
    mesh=trimesh.Trimesh(vertices=np.array(vertices),faces=np.array(faces),process=False)
    trimesh.repair.fix_normals(mesh, multibody=True)
    return mesh


def generate_rotor(config):
    samples=sections(config)
    blade=loft(samples["points"])
    d,g=config["design"],config["geometry"]
    hub=trimesh.creation.cylinder(radius=d["hub_radius"],height=d["hub_length"],sections=g["hub_sections"])
    if d["shaft_radius"] > 0:
        hole=trimesh.creation.cylinder(radius=d["shaft_radius"],height=2*d["hub_length"],sections=g["hub_sections"])
        hub=trimesh.boolean.difference([hub,hole],engine="manifold")
    blades=[]
    for k in range(d["blade_count"]):
        rotated=blade.copy()
        rotated.apply_transform(trimesh.transformations.rotation_matrix(2*np.pi*k/d["blade_count"],[0,0,1]))
        blades.append(rotated)
    mesh=trimesh.boolean.union([hub]+blades,engine="manifold",check_volume=True)
    if mesh is None or len(mesh.faces)==0: raise ValueError("Empty rotor Boolean union")
    diameter=2*float(np.max(np.linalg.norm(mesh.vertices[:,:2],axis=1)))
    flow=samples["chord_flow"]
    active=(samples["section_speed"]>.1)&(samples["u"]>.15)&(samples["u"]<.85)
    reverse=float(np.mean(flow[active]<0)) if np.any(active) else 0.
    warnings=[]
    if reverse>.05:
        warnings.append("Nominated profile LE sees reverse chordwise flow at %.1f%% of active stations; induced flow is omitted."%(100*reverse))
    if g.get('airfoil_mapping')=='opposed_naca':
        warnings.append('NACA-derived fore-aft mapping uses finite nose caps and an apex profile morph; aerodynamic improvement is unverified.')
    warnings.append("Root retraction is a geometric transition, not a structural certification or a tangent fillet.")
    measurements={"diameter":diameter,"volume":float(mesh.volume),"surface_area":float(mesh.area),
                  "triangles":len(mesh.faces),"blade_count":d["blade_count"],
                  "maximum_reference_curvature":float(samples["curvature"].max()),
                  "reverse_flow_station_fraction":reverse,"backend":"explicit-section-loft+manifold",
                  "root_transition":"quintic-chord-retraction-v1"}
    measurements['airfoil_mapping']=g.get('airfoil_mapping','fixed_naca')
    measurements['geometric_reverse_flow_station_fraction']=float(np.mean(samples['geometric_chord_flow'][active]<0)) if np.any(active) else 0.
    return GeometryResult(mesh,blades,hub,samples,measurements,warnings)

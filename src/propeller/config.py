"""Strict, explicit SI configuration. JSON avoids an extra parser dependency."""
import copy
import json
import math
from pathlib import Path

PARAMETERS = ("H", "beta", "shape_exponent", "b0", "b1", "b2", "Z", "alpha0", "alpha1")


def load_config(path):
    config = json.loads(Path(path).read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    validate_config(config)
    return config


def validate_config(c):
    required = {"schema_version", "design", "parameters", "bounds", "operating", "geometry", "cfd", "optimization"}
    if set(c) != required or c["schema_version"] != 1:
        raise ValueError("Configuration requires schema_version=1 and exactly: " + ", ".join(sorted(required)))
    def finite_tree(value):
        if isinstance(value, dict):
            for v in value.values(): finite_tree(v)
        elif isinstance(value, list):
            for v in value: finite_tree(v)
        elif isinstance(value, float) and not math.isfinite(value):
            raise ValueError("Non-finite configuration value")
    finite_tree(c)
    d, p, g, o = (c[k] for k in ("design", "parameters", "geometry", "operating"))
    def keys(name,required,optional=()):
        value=c[name]
        if not isinstance(value,dict):raise ValueError(name+' must be an object')
        missing=set(required)-set(value);unknown=set(value)-set(required)-set(optional)
        if missing or unknown:raise ValueError(name+': missing='+str(sorted(missing))+', unknown='+str(sorted(unknown)))
    keys('design',('diameter_max','blade_count','hub_radius','hub_length','root_radius','shaft_radius'))
    keys('operating',('rpm','rotation_sign','inflow_speed','density','viscosity','sound_speed','target_thrust','reference_power'))
    keys('geometry',('span_stations','chord_points','thickness_ratio','hub_sections','minimum_clearance','section_anchor','root_transition_fraction','minimum_root_area','minimum_body_thickness','minimum_trailing_edge_thickness','intersection_tolerance'))
    keys('cfd',('backend','launcher','version','fidelity','background_cells','surface_refinement','wake_refinement','rotor_refinement','layers','max_cells','iterations','timeout_seconds','force_window','force_relative_tolerance','residual_tolerance','continuity_tolerance','max_non_orthogonality','max_skewness','qualified'),('mpi_ranks','startup_timeout_seconds','domain'))
    keys('optimization',('seed','temperature','cooling','proposal_scale','thrust_penalty','max_proposals','max_evaluations','max_seconds'))
    for key in ("diameter_max", "hub_radius", "hub_length", "root_radius"):
        if d[key] <= 0: raise ValueError(key + " must be positive")
    if not isinstance(d["blade_count"], int) or d["blade_count"] < 1:
        raise ValueError("blade_count must be a positive integer")
    if not 0 <= d["shaft_radius"] < d["root_radius"] < d["hub_radius"] < d["diameter_max"]/2:
        raise ValueError("Require shaft_radius < root_radius < hub_radius < diameter_max/2")
    if set(p) != set(PARAMETERS) or set(c["bounds"]) != set(PARAMETERS):
        raise ValueError("Exactly nine named parameters and bounds are required")
    for key in PARAMETERS:
        lo, hi = c["bounds"][key]
        if not lo < hi: raise ValueError("Invalid bounds for " + key)
        if not lo <= p[key] <= hi: raise ValueError(key + " is outside configured bounds")
    if p["shape_exponent"] != 1 and p["shape_exponent"] < 2:
        raise ValueError("Root curvature requires exponent=1 or exponent>=2")
    if p["H"] <= 0 or not 0 < p["beta"] < math.pi or min(p[k] for k in ("b0","b1","b2")) <= 0:
        raise ValueError("Invalid loop dimensions")
    if o["rpm"] <= 0 or o["density"] <= 0 or o["viscosity"] <= 0 or o["sound_speed"] <= 0:
        raise ValueError("Operating rpm and fluid properties must be positive")
    if o["rotation_sign"] not in (-1, 1) or o["inflow_speed"] < 0:
        raise ValueError("Invalid rotation sign/inflow speed")
    for key in ("target_thrust", "reference_power"):
        if o[key] is not None and o[key] <= 0: raise ValueError(key + " must be positive or null")
    if g["span_stations"] < 17 or g["chord_points"] < 9 or not 0 < g["thickness_ratio"] < .5:
        raise ValueError("Insufficient geometry resolution or invalid thickness")
    if c["cfd"]["backend"] != "native" or c["cfd"]["version"] != "v2512":
        raise ValueError("This implementation supports native OpenCFD v2512")
    if c["optimization"]["cooling"] <= 0 or c["optimization"]["cooling"] > 1:
        raise ValueError("cooling must lie in (0,1]")
    f=c["cfd"];opt=c["optimization"]
    if f['qualified']:raise ValueError('Qualification must be established from evidence; it cannot be asserted in the configuration')
    if not isinstance(f.get("mpi_ranks",1),int) or not 1<=f.get("mpi_ranks",1)<=64:raise ValueError("mpi_ranks must be an integer in [1,64]")
    if f["force_window"]<4 or f["iterations"]<f["force_window"]:
        raise ValueError("iterations must cover force_window >= 4")
    if len(f["background_cells"])!=3 or any(not isinstance(n,int) or n<4 for n in f["background_cells"]):
        raise ValueError("background_cells requires three integers >=4")
    if len(f["surface_refinement"])!=2 or not 0<=f["surface_refinement"][0]<=f["surface_refinement"][1]<=10:
        raise ValueError("Invalid surface refinement levels")
    for value in f['surface_refinement']+[f['rotor_refinement'],f['wake_refinement'],f['layers']]:
        if not isinstance(value,int) or value<0:raise ValueError('Refinement/layer counts must be non-negative integers')
    if 'domain' in f:
        if set(f['domain'])!={'radial_diameters','upstream_diameters','downstream_diameters'} or any(v<=1 for v in f['domain'].values()):
            raise ValueError('Domain extents require radial/upstream/downstream diameters > 1')
    for key in ("timeout_seconds","max_cells","force_relative_tolerance","residual_tolerance","continuity_tolerance"):
        if f[key]<=0:raise ValueError(key+" must be positive")
    for key in ("temperature","proposal_scale","thrust_penalty","max_proposals","max_evaluations","max_seconds"):
        if opt[key]<=0:raise ValueError(key+" must be positive")
    if g["section_anchor"]!="inner_leading_edge":raise ValueError("Unsupported section anchoring convention")
    if not 0<g["root_transition_fraction"]<.25:raise ValueError("Invalid root transition fraction")
    if any(g[k]<0 for k in ("minimum_clearance","minimum_root_area","minimum_body_thickness","minimum_trailing_edge_thickness")):
        raise ValueError("Feature constraints cannot be negative")
    return c


def with_parameters(config, parameters):
    result = copy.deepcopy(config)
    result["parameters"] = dict(parameters)
    return validate_config(result)


def normalized(config):
    return [(config["parameters"][k]-config["bounds"][k][0])/(config["bounds"][k][1]-config["bounds"][k][0]) for k in PARAMETERS]


def physical(config, q):
    return {k: float(config["bounds"][k][0]+q[i]*(config["bounds"][k][1]-config["bounds"][k][0])) for i,k in enumerate(PARAMETERS)}

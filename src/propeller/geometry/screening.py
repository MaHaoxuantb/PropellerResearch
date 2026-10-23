"""Cheap section checks and a fixed-radius map; no flow predictions."""
import copy
import math
import numpy as np
from scipy.optimize import brentq
from ..config import PARAMETERS, validate_config
from .core import sections

FREE_PARAMETERS = tuple(k for k in PARAMETERS if k != 'H')


def section_screen(config, samples=None):
    s = sections(config) if samples is None else samples
    d, g = config['design'], config['geometry']
    points = s['points']; codes = []
    if not np.isfinite(points).all():
        return {'valid': False, 'codes': ['NONFINITE'], 'kind': 'section-screen'}
    radius = max(float(np.linalg.norm(points[..., :2], axis=-1).max()), d['hub_radius'])
    if 2*radius > d['diameter_max']+1e-7: codes.append('DIAMETER')
    for root in (points[0], points[-1]):
        if np.linalg.norm(root[:, :2], axis=1).max() >= d['hub_radius'] or np.abs(root[:, 2]).max() >= d['hub_length']/2:
            codes.append('ROOT_CAP_EXPOSED')
    if s['width'].min()*g['thickness_ratio'] < g['minimum_body_thickness']:
        codes.append('BODY_THICKNESS')
    for fraction in (.01, .25, .5, .75, .99):
        mid = s['C']+s['width'][:, None]*(fraction-s['root_shift'])[:, None]*s['d']
        active = np.linalg.norm(mid[:, :2], axis=1) > d['hub_radius']
        jacobian = np.sum(np.gradient(mid, s['u'], axis=0)*s['T'], axis=1)
        if np.any(jacobian[active] <= 0): codes.append('LOCAL_FOLD')
    o = config['operating']
    mach = math.hypot(o['rpm']*math.pi/30*radius, o['inflow_speed'])/o['sound_speed']
    if mach > .3: codes.append('MACH_UNSUPPORTED')
    return {'valid': not codes, 'codes': sorted(set(codes)), 'radius': radius,
            'kind': 'section-screen', 'note': 'Passing requires full solid validation before CFD.'}


def fixed_radius_config(base, q, radius, tolerance=5e-5):
    q = np.asarray(q, dtype=float)
    if q.shape != (len(FREE_PARAMETERS),) or not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError('Expected eight normalized parameters in [0,1]')
    if not base['design']['hub_radius'] < radius <= base['design']['diameter_max']/2:
        raise ValueError('Fixed radius must lie outside hub and inside diameter bound')
    c = copy.deepcopy(base)
    for k, v in zip(FREE_PARAMETERS, q):
        lo, hi = c['bounds'][k]; c['parameters'][k] = float(lo+v*(hi-lo))
    def error(h):
        c['parameters']['H'] = float(h)
        p = sections(c)['points']
        return float(np.linalg.norm(p[..., :2], axis=-1).max())-radius
    lo, hi = c['bounds']['H']
    a, b = error(lo), error(hi)
    if a*b > 0: raise ValueError('FIXED_RADIUS_UNREACHABLE')
    h = brentq(error, lo, hi, xtol=1e-10)
    c['parameters']['H'] = float(h)
    report = section_screen(c)
    if not report['valid']: raise ValueError(','.join(report['codes']))
    if abs(report['radius']-radius) > tolerance: raise ValueError('FIXED_RADIUS_MISMATCH')
    return validate_config(c)


def free_coordinates(config):
    return [(config['parameters'][k]-config['bounds'][k][0])/(config['bounds'][k][1]-config['bounds'][k][0]) for k in FREE_PARAMETERS]

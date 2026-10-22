"""Read-only audit of the CFD wall's geometric approximation to the rotor."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import trimesh

from ..experiments import atomic_json


def surface_audit(attempt, output):
    attempt=Path(attempt).resolve();output=Path(output).resolve();case=attempt/'case'
    manifest=json.loads((attempt/'manifest.json').read_text())
    if not json.loads((attempt/'mesh.json').read_text()).get('valid'):
        raise ValueError('Surface audit requires an accepted static mesh')
    if output.exists() and any(output.iterdir()):raise FileExistsError('Audit directory must be empty')
    output.mkdir(parents=True,exist_ok=True)
    wall=output/'rotor-wall.stl'
    command=[manifest['config']['cfd']['launcher'],'surfaceMeshExtract',str(wall),'-patches','rotor','-time','0','-case',str(case)]
    with (output/'extraction.log').open('w') as log:
        process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    if process.returncode or not wall.exists():raise ValueError('Native wall extraction failed; inspect '+str(output/'extraction.log'))
    source=trimesh.load(case/'constant/triSurface/rotor.stl',force='mesh')
    target=trimesh.load(wall,force='mesh')
    def distances(points,mesh):
        values=[]
        for start in range(0,len(points),256):
            _,distance,_=trimesh.proximity.closest_point(mesh,points[start:start+256])
            values.append(distance)
        distance=np.concatenate(values)
        return {'samples':len(distance),'maximum_m':float(distance.max()),'p95_m':float(np.percentile(distance,95)),
                'median_m':float(np.median(distance)),'rms_m':float(np.sqrt(np.mean(distance**2)))}
    # Vertices expose lost edges; face centroids expose faceting between them.
    edge_counts=np.bincount(target.edges_unique_inverse)
    result={'kind':'sampled-wall-geometry-audit','source_attempt':str(attempt),'command':command,
            'cad_to_wall':distances(np.vstack((source.vertices,source.triangles_center)),target),
            'wall_to_cad':distances(np.vstack((target.vertices,target.triangles_center)),source),
            'cad_volume_m3':float(abs(source.volume)),'oriented_wall_volume_estimate_m3':float(abs(target.volume)),
            'relative_oriented_volume_difference':float(abs(target.volume)/abs(source.volume)-1),
            'relative_area_error':float(target.area/source.area-1),'wall_watertight':bool(target.is_watertight),
            'export_boundary_edges':int((edge_counts==1).sum()),'export_nonmanifold_edges':int((edge_counts>2).sum()),
            'wall_triangles':len(target.faces),'wall_sha256':hashlib.sha256(wall.read_bytes()).hexdigest(),
            'qualified':False,'interpretation':'Distances at all vertices and face centroids, not an exact Hausdorff bound or a flow-discretization uncertainty estimate. STL coordinate merging/triangulation can affect export topology; export watertightness is not a verdict on native polyMesh closure. Volume is an oriented-surface estimate. The archived case is read only.'}
    atomic_json(output/'audit.json',result)
    return result

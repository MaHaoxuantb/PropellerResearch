"""Generate a self-contained OpenCFD v2512 incompressible SST/MRF case."""
import json
import math
from pathlib import Path


def header(name, cls="dictionary"):
    return 'FoamFile { version 2.0; format ascii; class %s; object %s; }\n'%(cls,name)


def build_case(geometry, config, directory):
    root=Path(directory)
    if root.exists() and any(root.iterdir()): raise FileExistsError("Case directory must be empty")
    for part in ("0","constant/triSurface","system"): (root/part).mkdir(parents=True,exist_ok=True)
    c,o=config["cfd"],config["operating"]
    D=config["design"]["diameter_max"]; rho=o["density"]; nu=o["viscosity"]/rho
    omega=o["rotation_sign"]*o["rpm"]*math.pi/30
    mach=math.hypot(abs(omega)*geometry.measurements["diameter"]/2,o["inflow_speed"])/o["sound_speed"]
    if mach>.3: raise ValueError("Estimated Mach >0.3: compressible backend is required and not implemented")
    geometry.mesh.export(root/"constant/triSurface/rotor.stl")
    def put(path,text,cls="dictionary"):
        (root/path).write_text(header(Path(path).name,cls)+text)
    put("constant/transportProperties",f"transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] {nu:.12g};\n")
    put("constant/turbulenceProperties","simulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }\n")
    put("constant/MRFProperties",f"rotatingFrame {{ cellZone rotorZone; active yes; nonRotatingPatches (farfield); origin (0 0 0); axis (0 0 1); omega {omega:.12g}; }}\n")
    put("system/controlDict",f'''application simpleFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {c['iterations']}; deltaT 1;
writeControl timeStep; writeInterval {max(50,c['iterations']//5)}; purgeWrite 2;
writeFormat binary; writePrecision 10; writeCompression off; timeFormat general; timePrecision 8; runTimeModifiable false;
functions
{{
 forces {{ type forces; libs (forces); patches (rotor); rho rhoInf; rhoInf {rho}; CofR (0 0 0);
    writeControl timeStep; writeInterval 1; log true; }}
 yPlus {{ type yPlus; libs (fieldFunctionObjects); writeControl writeTime; log true; }}
}}
''')
    # +z thrust corresponds to a -z wake; leave a long downstream volume.
    domain=c.get('domain',{'radial_diameters':3,'upstream_diameters':3,'downstream_diameters':8})
    x=domain['radial_diameters']*D; z0=-domain['downstream_diameters']*D; z1=domain['upstream_diameters']*D
    vertices=[(-x,-x,z0),(x,-x,z0),(x,x,z0),(-x,x,z0),(-x,-x,z1),(x,-x,z1),(x,x,z1),(-x,x,z1)]
    vertices_text="\n".join("(%g %g %g)"%v for v in vertices)
    cells=" ".join(map(str,c["background_cells"]))
    # Use the centre of a far-field background cell, not an axis or a grid
    # plane. Native point location can reject points on coincident cell faces.
    # The positive-x/y corner is outside the diameter-constrained rotor even
    # for the smallest supported domain and background resolution.
    nx,ny,nz=c['background_cells']
    fluid_seed=(x-x/nx, x-x/ny, z1-(z1-z0)/(2*nz))
    fluid_seed_text=' '.join(format(value,'.15g') for value in fluid_seed)
    put("system/blockMeshDict",f'''scale 1; vertices ({vertices_text});
blocks (hex (0 1 2 3 4 5 6 7) ({cells}) simpleGrading (1 1 1)); edges ();
boundary (farfield {{ type patch; faces ((0 3 2 1) (4 5 6 7) (0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7)); }});
mergePatchPairs ();
''')
    put("system/surfaceFeatureExtractDict",'''rotor.stl { extractionMethod extractFromSurface; extractFromSurfaceCoeffs { includedAngle 150; } writeObj no; }''')
    mesh_max=max(abs(geometry.mesh.bounds[:,2]).max(),.1*D)
    zone_z=float(mesh_max+.12*D); zone_r=.58*D
    put("system/topoSetDict",f'''actions
(
 {{ name rotatingCells; type cellSet; action new; source cylinderToCell; p1 (0 0 {-zone_z}); p2 (0 0 {zone_z}); radius {zone_r}; }}
 {{ name rotorZone; type cellZoneSet; action new; source setToCellZone; set rotatingCells; }}
);
''')
    lo,hi=c["surface_refinement"]
    quality=f'''maxNonOrtho {c['max_non_orthogonality']}; maxBoundarySkewness {c['max_skewness']}; maxInternalSkewness {c['max_skewness']}; maxConcave 80;
minVol 1e-16; minTetQuality 1e-15; minArea -1; minTwist .02; minDeterminant .005;
minFaceWeight .05; minVolRatio .01; minTriangleTwist -1; nSmoothScale 4; errorReduction .75;'''
    put("system/meshQualityDict",quality)
    put("system/snappyHexMeshDict",f'''castellatedMesh true; snap true; addLayers {'true' if c['layers'] else 'false'};
geometry
{{
 rotor.stl {{ type triSurfaceMesh; name rotor; }}
 rotatingRegion {{ type searchableCylinder; point1 (0 0 {-zone_z}); point2 (0 0 {zone_z}); radius {zone_r}; }}
 wake {{ type searchableCylinder; point1 (0 0 {-4*D}); point2 (0 0 {zone_z}); radius {zone_r}; }}
}}
castellatedMeshControls
{{
 maxLocalCells {c['max_cells']}; maxGlobalCells {c['max_cells']}; minRefinementCells 0;
 maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
 features ({{ file "rotor.eMesh"; level {hi}; }});
 refinementSurfaces {{ rotor {{ level ({lo} {hi}); patchInfo {{ type wall; }} }} }}
 resolveFeatureAngle 30;
 refinementRegions {{ rotatingRegion {{ mode inside; levels ((1e15 {c['rotor_refinement']})); }} wake {{ mode inside; levels ((1e15 {c['wake_refinement']})); }} }}
 locationInMesh ({fluid_seed_text}); allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2; nSolveIter 30; nRelaxIter 5; nFeatureSnapIter 10; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false; }}
addLayersControls
{{
 relativeSizes true; layers {{ rotor {{ nSurfaceLayers {c['layers']}; }} }}
 expansionRatio 1.2; finalLayerThickness .4; minThickness .08;
 nGrow 0; featureAngle 60; nRelaxIter 5; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
 maxFaceThicknessRatio .5; maxThicknessToMedialRatio .3; minMedialAxisAngle 90;
 nBufferCellsNoExtrude 0; nLayerIter 50; nRelaxedIter 20;
}}
meshQualityControls
{{
 {quality}
}}
mergeTolerance 1e-7;
''')
    put("system/fvSchemes",'''ddtSchemes { default steadyState; }
gradSchemes { default cellLimited Gauss linear 1; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }
wallDist { method meshWave; }
''')
    put("system/fvSolution",'''solvers
{
 p { solver GAMG; tolerance 1e-8; relTol .05; smoother GaussSeidel; }
 "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol .1; }
}
SIMPLE { nNonOrthogonalCorrectors 1; consistent yes; }
relaxationFactors { fields { p .25; } equations { U .5; k .5; omega .5; } }
''')
    put("system/decomposeParDict",f"numberOfSubdomains {c.get('mpi_ranks',1)}; method scotch;\n")
    speed=max(o["inflow_speed"],1.0)
    k=1.5*(.01*speed)**2; turbulent_omega=math.sqrt(k)/(.09**.25*.07*D)
    if o["inflow_speed"]>0:
        farU=f'type freestreamVelocity; freestreamValue uniform (0 0 {-o["inflow_speed"]});'
        farP='type freestreamPressure; freestreamValue uniform 0;'
    else:
        farU='type pressureInletOutletVelocity; value uniform (0 0 0);'
        farP='type totalPressure; p0 uniform 0; value uniform 0;'
    fields={
      "U":("volVectorField","0 1 -1 0 0 0 0",f'(0 0 {-o["inflow_speed"]})',farU,"type noSlip;"),
      "p":("volScalarField","0 2 -2 0 0 0 0","0",farP,"type zeroGradient;"),
      "k":("volScalarField","0 2 -2 0 0 0 0",str(k),f"type inletOutlet; inletValue uniform {k}; value uniform {k};","type kqRWallFunction; value uniform "+str(k)+";"),
      "omega":("volScalarField","0 0 -1 0 0 0 0",str(turbulent_omega),f"type inletOutlet; inletValue uniform {turbulent_omega}; value uniform {turbulent_omega};","type omegaWallFunction; value uniform "+str(turbulent_omega)+";"),
      "nut":("volScalarField","0 2 -1 0 0 0 0","0","type calculated; value uniform 0;","type nutUSpaldingWallFunction; value uniform 0;")}
    for name,(cls,dim,initial,far,wall) in fields.items():
        put("0/"+name,f"dimensions [{dim}];\ninternalField uniform {initial};\nboundaryField {{ farfield {{ {far} }} rotor {{ {wall} }} }}\n",cls)
    (root/"case.foam").touch()
    (root/"configuration.json").write_text(json.dumps(config,indent=2,allow_nan=False))
    (root/"geometry.json").write_text(json.dumps(geometry.measurements,indent=2))
    return root

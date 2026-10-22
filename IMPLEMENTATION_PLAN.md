# Toroidal propeller implementation plan

Planning baseline: 22 October 2026. Based on the supplied specification v1.0 and the user's follow-up decisions. This document plans implementation; no geometry or CFD results are claimed.

Implementation update: the user subsequently required the Homebrew-installed native OpenFOAM. The implemented backend is the discovered OpenCFD v2512 macOS application, not Docker/v2606. See README.md and docs/IMPLEMENTATION_NOTES.md for implemented scope and verification evidence; the original phase gates below remain research milestones.

## 1. Scope and decisions

Build an automated pipeline for a continuous U-shaped blade, attached to the hub at both ends, with rotationally repeated identical units. Preserve the nine design variables, quadratic Bernstein chord, analytic axial rise, simulated annealing (SA), and staged CFD validation in the proposal.

The reference curve is a **geometric edge**, not an assumed aerodynamic leading edge. Section orientation must be established before accepting a baseline. Use explicit airfoil sections in local frames and loft them into a solid, rather than surface-normal thickening as the production method.

The implementation should find designs, without requiring the user to specify an optimized geometry. Operating requirements belong to reusable experiment configurations. The initial objective remains minimum shaft power subject to target thrust, fixed RPM, maximum diameter, and geometric constraints. Do not add the constant-blade-area constraint suggested before the full proposal was supplied: it is not part of this specification. Record area and volume as outputs.

No universal optimum is implied. A geometry optimized for one scale, Reynolds number, RPM, or inflow can perform differently elsewhere. Later experiments can vary these conditions or use multiple operating points. If RPM is eventually optimized, it needs its own bounded variable and motor/operating constraints; do not silently add it to the nine-variable search.

### Execution environment

Use Python orchestration on macOS and the native Homebrew-installed OpenCFD v2512 application at `/Applications/OpenFOAM-v2512.app`. Invoke its bundled launcher to mount and activate the ARM64 runtime. This supersedes the original Docker/v2606 plan at the user's explicit request. No container runtime is required.

Pin the discovered distribution/version in each configuration and fingerprint executable/library files in evaluation records. Serial native meshing and solving have been verified. The optional MPI path is not qualified on this host because runtime startup stalled; the supplied configurations use one rank.

The existing Python 3.9.6 environment now contains the dependencies pinned in `requirements.lock`. The original `Model.py` and `shape.stl` remain historical prototypes; `main.py` is the new package entry point.

## 2. Mathematical decisions required before solid generation

### Root regularity

For r(u) = r_h + H sin(pi u)^p, near u=0:

    r(u) - r_h ~ H pi^p u^p
    r''(u) ~ H pi^p p(p-1) u^(p-2)

Thus 1 < p < 2 has unbounded second derivative at the endpoint. Since beta and r_h are positive, the projected tangent is nonzero there, and this produces unbounded reference-curve curvature. The specification's claim of smoothness over the full proposed interval needs qualification. p=1 is a special regular case; the asymptotic coefficient above vanishes.

Implementation decision: begin continuous optimization with **2 <= p <= 3**, using p=2 for the first baseline. Preserve p=1 as a separately testable geometry. Supporting the entire interval [1,3] requires explicit finite-length root regularization with matched position/tangent/curvature and a versioned transition rule. This is a deliberate restriction of the provisional bounds, not an unnoticed change.

### Frame orientation and section topology

Compute the specified projected-normal frame analytically, validate orthonormality, and make its sign a documented convention. At the apex u=1/2, for positive beta, the specified e1 points radially inward. It must be negated if the intended chord there extends outward; then recompute e2 = T cross e1. Never flip signs independently at individual stations.

At the same apex, r'=z'=0 and T is azimuthal. The rotor's rotational velocity is therefore parallel to the local span tangent. This is an inherent strongly three-dimensional part of the shape. A conventional two-dimensional section interpretation cannot be assumed throughout the loop.

For a diagnostic flow estimate, compute:

    U_rel = U_infinity - Omega_vector cross C
    U_section = U_rel - (U_rel dot T) T

Plot its direction against the chord at both legs and the apex. This neglects induced flow and is only an orientation diagnostic, not an aerodynamic prediction. Later repeat the inspection using CFD velocities.

Use continuous section correspondence: the same section vertex indices connect between neighboring stations. Do not reverse leading/trailing edges station by station to follow a changing flow direction, since this can create a discontinuity. If a consistent NACA section mapping produces substantial reverse-flow operation on one leg, revisit the twist law or section topology before optimization. NACA 0012 is symmetric above/below its chord, not fore/aft.

### Airfoil mapping

Use the proposal's explicit-section alternative. With b = T cross d, map dimensionless 2D section coordinates into the local section plane:

    X(u,j) = C(u) + w(s(u)) [x_j d(u) + y_j b(u)]

Here the definition of x_j includes which section edge is anchored to C. Make that convention explicit in the geometry model. A surface normal computed from S_u cross S_v generally varies along the chord; displacing along it does not preserve an exact planar NACA section. Keep that construction only as an optional future geometry backend.

Choose and record the NACA formula, thickness convention, trailing-edge closure, and sampling method. The supplied f_t convention is ambiguous unless its normalization is specified. Validate the resulting thickness/chord ratio directly. Use cosine-spaced chord samples and adaptive span sampling, with consistent correspondence.

### Hub attachment and feature limits

The root reference radius and hub outer radius are different quantities. Two capped solids touching at an edge or plane do not guarantee a robust union. Add finite root overlap and deterministic transition geometry into the fixed hub. Limit the overlap exemption to named attachment regions; blade collisions elsewhere remain invalid.

A sharp airfoil edge has locally zero thickness. Therefore a blanket minimum-thickness test would reject every sharp-edged airfoil. Define body thickness, trailing-edge thickness/radius, leading-edge radius, and root ligament dimensions separately. Use a fixed, recorded finite trailing-edge treatment for meshable solids, with sensitivity checks before claiming small gains.

Require a valid CAD solid, valid exported surface mesh, and explicit intersection tests. `is_watertight` alone is insufficient. CAD tessellation tolerances must enter the manifest and be small enough that diameter, clearance, and root measurements are stable under refinement.

## 3. Configuration and demonstration experiment

Separate four immutable configuration groups: geometry/bounds, operating conditions, numerical fidelity, and optimization policy. Use SI internally, radians in the model, and explicit degree conversion at user input. Give the shape exponent a clear code name such as `shape_exponent`, separate from pressure and the parameter vector.

Use dimensionless geometric inputs relative to D_max for exploration, but convert to physical lengths before CFD. A dimensionless geometry alone does not determine Reynolds number or physical performance.

Provide an **illustrative software demonstration**, not a user-approved engineering requirement:

| Setting | Initial demonstration value |
|---|---|
| D_max | 0.254 m |
| RPM | 5,000 |
| Inflow | Static air |
| Air properties | rho=1.225 kg/m^3, mu=1.81e-5 Pa s, sound speed=343 m/s |
| Rotation / thrust convention | Positive rotation about +z; desired force on rotor +z |
| Loop count | 2 initially; 3 and 4 after pipeline qualification |
| Geometry seed | H/D=.30, beta=.5 pi/N_b, p=2, b_i/D=(.06,.08,.06), Z/D=.06 |
| Twist seed | alpha_0=15 degrees, alpha_1=0; provisional until orientation inspection |
| Hub seed | Outer radius .08D, axial length .08D; root anchor .07D plus explicit attachment blend |

These values only seed the geometry and runtime tests. They are not assertions of aerodynamic or mechanical adequacy. The approximate tip-speed Mach upper bound at D_max is 0.194, before induced/local velocity effects; check actual flow results before accepting incompressibility.

For the demonstration, first establish a positive-thrust, positive-consumed-power, converged baseline, then set a fixed demonstration target thrust from that baseline and record it before launching optimization. Never move the target or P_ref during the run. If the seed fails, revise the baseline openly and record the revision; do not fabricate a thrust target from failed output. A physical application can later replace these example requirements without changing code.

## 4. Implementation milestones and acceptance gates

| Milestone | Implementation and deliverables | Required evidence to advance |
|---|---|---|
| M0: Project contracts | `pyproject.toml`, locked environment, configuration schemas, typed records, CLI, version manifest, source specification snapshot | Invalid/missing configurations rejected clearly; units and parameter normalization round-trip; environment report distinguishes installed vs usable tools |
| M1: Geometry prototype | Analytic curve/derivatives, planar arc-length mapping, Bernstein chord, frames, twist, replication; annotated 3D preview and parameter sweeps | Endpoints, apex, width bounds, frame handedness, continuity and rotational symmetry verified; both loop legs' edge/flow orientation reviewed |
| M2: One robust solid | Airfoil section mapping, loft, finite trailing edge, root transitions, hub union and replication; STEP/STL plus region metadata | Positive-volume connected rotor with intended holes, no unintended intersections, oriented manifold surface, valid root connections and stable dimensions |
| M3: Geometry screening | Staged validator, structured violations, adaptive tessellation, randomized/boundary samples, failure classification | Known bad fixtures rejected; hundreds of sampled candidates characterized; every accepted candidate passes all mandatory checks; feasible neighborhood exists around baseline |
| M4: Reproducible CFD case | Pinned runtime, one full-rotor steady MRF template, automatic mesh, force extraction and convergence policy | One baseline completes repeatedly; forces/torques/signs/density treatment verified; mesh/layer quality and domain sensitivity assessed |
| M5: Evaluation service | Candidate IDs, immutable run folders, cache, timeout/retry handling, resource limits, durable state and resume | Interrupted runs recover; failures never become performance values; cache keys distinguish fidelity/physics/version changes |
| M6: Pilot and SA | Pilot runtime/noise measurements, sensitivities, normalized proposals, temperature calibration, checkpointing and independent restarts | SA logic passes analytic/mock tests; separate best-feasible archive; reproducible resumed chain; affordable measured CFD budget |
| M7: Confirmation and report | Coarse/medium/fine studies, transient finalists, uncertainty, robustness, matched comparisons | Claimed gain survives refinement/model checks; target thrust verified with uncertainty considered; negative findings reported honestly |
| M8: Extensions | Separate blade-count studies, asymmetric geometry, splines, CST, surrogate and alternative optimizers; structural/manufacturing workflow | Each extension uses controlled budgets and demonstrates value relative to the validated baseline |

The critical path is M0 → M1 → M2 → M3 → M4 → M5 → M6 → M7. M5's storage/contracts can be implemented while M2–M4 are developed. Physics validation remains a gate before costly optimization. Keep M1 as the first reviewable deliverable, not an entire optimizer before the blade shape can be inspected.

### CAD backend decision

Use a short CadQuery/OpenCascade feasibility spike for lofting, root blends and union, with Trimesh for mesh inspection/export checks. CadQuery provides the relevant operations, but their availability does not establish robustness for these shapes; see its [class reference](https://cadquery.readthedocs.io/en/latest/classreference.html).

Test one baseline and several severe twist/width/curvature examples before committing to the CAD backend. Prefer deterministic constructed root transitions where generic fillets fail. Any alternative union/repair backend must be versioned and its geometric changes measured. Do not silently voxelize or smooth the candidate into a different design to make CFD succeed.

### Meaningful geometry verification

- Analytic endpoint/apex results and derivative comparisons away from singular cases.
- Arc-length mapping is monotone and stable as sampling is refined.
- Width stays in the coefficient convex hull; b1 is a control coefficient, not generally w(.5).
- Local frames are orthonormal and continuous; replicated vertices preserve radii and z.
- Section reconstruction recovers the intended profile and thickness to tolerance.
- Fixtures exercise self-intersection, adjacent-blade collision, detached roots, open seams, inverted faces and near-limit diameter.
- Refinement establishes stability of volume, diameter, clearance and curvature estimates.

Report rejection rates by reason. No valid-candidate percentage is assumed in advance. Tighten proposal bounds or improve parameterization based on observed failures, with a new experiment/configuration version.

## 5. CFD details that must be explicit

Start with the full rotor; periodic-sector savings require separate validation. The flow regime determines incompressible versus compressible configuration. Maintain separate backend templates rather than attempting to toggle density/energy behavior with a few ambiguous keys.

Build the initial case from a working v2512 rotating-flow tutorial/configuration. Generate the background mesh, extract features, run snappyHexMesh, inspect layers and mesh quality, then run the solver. Capture actual layer coverage and y+ distributions: a mesher exit code is not evidence of adequate near-wall resolution. [OpenFOAM's meshing guide](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.4-mesh-generation-with-the-snappyhexmesh-utility) documents the staged process.

Record the rotating cell zone, axis/origin, angular speed and rotating-wall treatment. Include the complete intended wetted rotor surface in the force integration, with moments about the rotor center. For kinematic-pressure solvers, configure density correctly; OpenFOAM distinguishes this from physical pressure in its [forces documentation](https://doc.openfoam.com/2306/tools/post-processing/function-objects/forces/forces/).

Record signed fluid moment M_z and signed Omega_z. Consumed shaft power is P = -M_z Omega_z when M is fluid torque on the rotor under this convention. Do not apply an unconditional absolute value that could hide a windmilling case or a sign error. Verify the function object's convention on a known case before use.

Every fidelity configuration specifies residual thresholds, force/torque drift and variability windows, conservation limits, maximum iterations/time, mesh limits and wall treatment. Calibrate thresholds in the pilot and freeze them before optimization. Persistent unsteady loads mean the steady solution is not adequate; do not accept an arbitrary last value.

Verify baseline mesh/domain sensitivity and compare a reproducible conventional-propeller reference against measured data where matching geometry is available. This checks the workflow, but does not experimentally validate a loop rotor. The [UIUC propeller database](https://m-selig.ae.illinois.edu/props/propDB.html) is a candidate reference source. Use a three-grid study and time-step refinement for transient finalists; a small fine/medium difference alone is not proof of convergence. See [NASA's verification and validation tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html).

For static operation report thrust, shaft power and figure of merit; mark conventional forward-flight efficiency as not useful for static comparisons. Report actual diameter in coefficients, while preserving D_max as the constraint. Retain pressure, viscous and total force contributions for diagnostics.

## 6. Evaluation records and optimization behavior

Use immutable parameter/configuration records and an evaluation state machine. Preserve the specified terminal states: GEOMETRY_INVALID, EXPORT_FAILED, MESH_FAILED, SOLVER_FAILED, NOT_CONVERGED and SUCCESS. Add orchestration states such as QUEUED, RUNNING, TIMEOUT, CANCELLED and INTERRUPTED without confusing them with aerodynamic infeasibility.

Separate four concepts: candidate geometry, evaluation at one fidelity, an execution attempt, and an optimizer proposal/acceptance event. One candidate may have multiple fidelities or retries, and a cache hit is still an optimizer event.

Use canonical finite-valued serialized inputs, generator/source version, operating conditions, airfoil/root conventions, geometry tolerances, complete numerical configuration, runtime version/digest and fidelity in the cache key. Hash exported geometry as an additional check. Never round nearby parameter vectors together. Failed attempts retain records and must not masquerade as cached successes.

Start with SQLite metadata and per-evaluation artifact directories, using a single writer on the host and atomic result publication. Retain configuration, parameters, generated geometry, commands, logs, histories, quality metrics and result records for every evaluation. Control large field/checkpoint output explicitly; retain finalist fields and enough restart/reproduction material. No automatic `clean()` after evaluation as in the current starter script.

For SA:

1. Normalize variables, use symmetric reflected proposals, and start from a validated baseline.
2. Use the proposed fixed normalized power plus thrust-deficit penalty. Calibrate P_ref, penalty and initial temperature before the production chain.
3. Keep the current search state separate from the best feasible design. A low penalized score does not establish feasibility.
4. Treat invalid geometry and failed CFD as rejected/unavailable candidates, never as zero-power results.
5. Use one frozen CFD fidelity within a chain. High-fidelity confirmation is a separate evaluation/ranking stage; do not mix incomparable scores in the acceptance rule.
6. Freeze proposal adaptation after warm-up, or explicitly document the adaptive algorithm. Record RNG state, temperature, current/best states, proposal scales and counters in each checkpoint.
7. Define cooling on successfully scored proposal events, including cache hits; invalid/failed proposals do not cool the chain. Independently cap total proposals, solver launches, failures and wall time to prevent an infinite failure loop.
8. Parallelize independent restarts or blade-count studies first. A single chain remains sequential. Limit concurrent jobs from measured CPU/RAM use, including meshing peaks and Linux runtime memory limits.

Use mock objectives and parser fixtures for unit/integration tests, clearly labelled as synthetic. Never put their values in the CFD result archive as physical data.

## 7. Pilot budget, selection and expansion

Begin with one loop count and one baseline. Run approximately 200–500 cheap geometry samples, concentrating additional samples near boundaries and the feasible baseline. Then use a small CFD pilot of roughly 10–20 geometrically accepted cases, adjusting downward if runtime is high. These are planning quantities, not a promise to launch jobs now.

Measure wall time, core-hours, peak RAM, disk use, failure rate, force uncertainty and coarse/fine ranking consistency. Estimate the campaign using measured case cost, not assumed throughput. Start with one CFD job at a time, then choose either more MPI ranks or more independent jobs from benchmarks. Do not divide elapsed time by core count and assume linear speedup.

After the pilot, choose restarts and per-chain evaluation budget. Compare runs by solver launches/core-hours as well as scored candidates, because cache hits and failures have different costs. Compare alternate optimizers only after SA and the evaluator work reliably, at equivalent budgets.

For finalist selection, require positive useful thrust, positive consumed power, all hard geometry constraints and adequate numerical evidence. Where a defensible thrust uncertainty bound U_T is available, require T_estimated - U_T >= T_target. Do not count an uncertain threshold crossing as demonstrated feasibility. A power improvement must exceed the relevant comparison uncertainty; re-evaluate the baseline and candidate at matched high fidelity.

If no gain survives confirmation, report that no improvement was established within the explored family and budget. That is a valid research outcome.

## 8. Planned package layout and interfaces

    pyproject.toml
    config/{design,constraints,operating_point,cfd,optimization}.yaml
    src/propeller/{config,cli,models}.py
    src/propeller/geometry/{loop,arclength,width,frames,twist,airfoil,loft,hub,assembly,validation}.py
    src/propeller/simulation/{case_builder,mesh,runner,parser,convergence}.py
    src/propeller/optimization/{objective,annealing,proposals}.py
    src/propeller/experiments/{database,scheduler,analysis}.py
    simulation/case_builder.py  # native OpenCFD v2512 dictionaries
    tests/{unit,integration,fixtures}/
    docs/
    runs/

Preserve the proposal's interfaces: generate_rotor, validate_geometry, evaluate_cfd and optimize. Define structured GeometryResult, ValidationReport and CFDResult records, including status, units, diagnostics and artifact paths. Introduce an evaluator interface so mock tests, coarse CFD and final CFD share orchestration without confusing provenance. Avoid empty surrogate/CST/benchmark frameworks until those phases begin.

Planned commands: `doctor`, `geometry preview`, `geometry validate`, `evaluate`, `optimize`, `resume`, and `report`. The first milestone should deliver a readable configuration, a dimensionless parameter preview, geometry diagnostics and reproducible files. CFD and optimization commands become usable only as their acceptance gates pass.

## 9. Remaining decisions and scope limits

The user delegated geometry choices and specified macOS with a current OpenFOAM. No further answer is necessary to begin M0–M1 under the decisions above.

Discover runtime resources and installed OpenFOAM details during environment setup. Finalize section orientation from the M1 diagnostic preview, and establish the baseline's usable parameter neighborhood during M2–M3. Determine the production compute budget from M4–M6 measurements rather than asking the user to guess it now.

Physical material, shaft mounting, manufacturing tolerances and allowable stresses are not specified. V1 is therefore a simulation research tool with configurable geometry constraints, not a structurally approved printable rotor. Do not claim mechanical acceptance from minimum root width alone. Before manufacture, obtain material/process/mounting inputs and assess centrifugal plus aerodynamic loads, deformation, fatigue and balance, followed by controlled experimental validation.

Noise optimization, airfoil optimization, aeroelasticity, multi-machine scheduling and surrogate-guided searches remain later phases as requested. Their omission from the first implementation is staged delivery, not removal from the research roadmap.

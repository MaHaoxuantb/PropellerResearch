# Implementation notes

## Decisions and departures from the proposal

- Native Homebrew OpenCFD v2512 replaces the planned Docker/v2606 backend at the user's request. No Homebrew tap trust settings were changed.
- Configuration uses strict JSON rather than YAML to avoid an additional runtime parser; SI units and explicit schemas remain.
- Explicit section triangulation plus Manifold union replaces the proposed CAD feasibility route for this first implementation. This backend was already installed, preserves the generated vertices except for Boolean intersections, and does not apply voxel reconstruction or smoothing. STEP export and CAD fillets are not claimed.
- Root chords retract into the hub through a quintic transition over the configured end interval. Root caps must be inside the hub. This is an explicit local departure from using C(u) as the section edge everywhere; away from root transitions the specified mapping is used exactly.
- The first airfoil is conventional symmetric NACA 0012 with its standard finite trailing-edge coefficient. The section-flow inspector explicitly exposes the opposite chordwise flow on the two loop legs. Reversing a section abruptly was rejected because it would introduce a discontinuity.
- Continuous shape exponents use [2,3]. The p=1 curve remains available with a separately edited configuration. No intermediate exponent with singular root curvature is accepted.
- Pitch bounds allow positive and negative alpha0. That enables both aerodynamic handedness possibilities without changing the prescribed rotation convention.
- CFD acceptance uses explicit mesh-quality limits with `checkMesh -allTopology -meshQuality`. Extended `-allGeometry` diagnostics are useful separately: they include convexity tests on cut/polyhedral cells that are not identical to the configured finite-volume quality criteria. No extended report is treated as equivalent to a passing configured check.

## Working features

Geometry generation, explicit sections, root assembly, rotational replication, mesh export, staged validation, offline interactive inspection, envelope sampling, native runtime discovery, complete steady MRF case generation, native meshing and solver execution, force/moment parsing, convergence classification, immutable records/cache, baseline target freezing, simulated annealing, checkpoint resume, bounded independent candidate scheduling, live run status, interrupted-attempt recovery and plotted evidence reports.

The inspector adds a useful diagnostic beyond the proposal: it links position on the loop to the section shape, chordwise flow sign and curvature. It deliberately displays geometry quality and aerodynamic uncertainty separately.

## Evidence collected during development

- Native macOS ARM64 OpenCFD v2512 binaries were discovered and executed through the Homebrew cask app's launcher.
- The initial two-loop baseline is approximately 0.22232 m in diameter and 3.6206e-5 m³ in solid volume, represented by 14,436 triangles. The complete rotor is connected, oriented and watertight; the triangle scan found no nonadjacent crossings.
- A seeded screen of 250 parameter vectors found 140 passing the cheap diameter/root-envelope checks. These are not 140 certified solids or CFD candidates.
- Initial native snappyHexMesh generation produced roughly 192,000 cells in about 27 seconds. Extended checks exposed weak cells/weights, so the generated quality settings were tightened before the integrated evaluation.
- The integrated case passed native `checkMesh -allTopology -meshQuality` and reached a real `simpleFoam` run with force and moment histories. Exact runtime, convergence and metrics are retained in the run directory, not inferred from a screenshot or solver exit code.
- A four-iteration native forward-inflow run verified that the generated freestream boundary conditions execute; it is a boundary-condition smoke test, not performance validation.
- Serial native execution works. MPI startup stalled even outside the sandbox and was terminated; optional MPI execution is not qualified on this host.
- All 26 unit/integration tests pass and exercise geometry invariants, deliberate invalid cases, exact resumed SA histories, cache identity and signed power extraction.

## Explicit remaining research/development scope

The implementation has not established grid/domain independence, transition-model accuracy or experimental validation. The initial no-induced-flow diagnostic exposes reverse chordwise flow on one leg; a geometry-valid baseline must not be described as an aerodynamically optimized blade. Strict convergence checks can therefore return NOT_CONVERGED, and this is a legitimate result rather than a reason to weaken acceptance thresholds.

Automatic transient/sliding-mesh CFD, compressible CFD, STEP/CAD fillets, production structural analysis, manufacturing constraints tied to a real material/process, distributed scheduling, uncertainty-qualified finalist selection, surrogates, CST, spline extensions and alternative optimizer benchmarking remain later milestones. The software refuses to call exploratory optimization qualified. No physical rotor has been manufactured or certified.

The retained original scripts are historical references. `main.py` now invokes the package CLI and no longer references a missing `trials` case or deletes results after running.

## Completed native CFD trials

The original constant-pitch seed completed 1,000 iterations but failed residual and force-stability acceptance. Its report is `runs/reports/original-baseline.html`.

The antisymmetric pitch seed completed 1,500 iterations on 193,726 cells in about 982 seconds including meshing. The final 100-sample thrust drift was 0.656% and power drift 0.192%, both below the 2% limits. Its axial-velocity initial residual was 1.1948e-4, exceeding the unchanged 1e-4 threshold. The result is **NOT_CONVERGED**, not an accepted baseline. Rotor y+ ranged from 1.30 to 194.09 (average 21.10); near-wall sensitivity remains unresolved. The report and plotted histories are `runs/reports/antisymmetric.html`.

Neither trial supplies an optimization objective or a frozen target. The next research step is a separately recorded longer run of the antisymmetric hypothesis, then mesh/domain and wall-resolution sensitivity if it converges. Preserve the failed trials and existing acceptance thresholds. A full optimization campaign has not been launched.

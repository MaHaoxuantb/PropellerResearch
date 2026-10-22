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

## Baseline continuation and sensitivity protocol

`continue-evaluation` creates a new attempt from an archived serial run. It checks the parent manifest/configuration, passing mesh record, rotor checksum, restart fields and native runtime fingerprints. It copies rather than edits the parent. The manifest records the parent result checksum and hashes of inherited dictionaries, mesh and fields. Only iteration count and timeout may change; all acceptance thresholds remain fixed. This path is covered by a synthetic file-level regression test in addition to the native continuation run.

The next baseline uses `config/antisymmetric-extended.json`: 3,000 total iterations, restarted from the saved 1,500 field. Its target thrust and reference power remain unset pending acceptance.

Prepared sensitivity configurations are:

- `config/antisymmetric-fine.json`: background resolution 25 × 25 × 40 instead of 20 × 20 × 32. At unchanged refinement levels the nominal cell edge is reduced by 20%, including the relative layer scale. This is one practical refinement step, not a formal three-grid uncertainty estimate.
- `config/antisymmetric-large-domain.json`: radial extent 4.5 diameters and downstream extent 13.5 diameters, with the upstream extent held at 3. Background resolution 30 × 30 × 48 preserves the original cell spacing and grid alignment around the rotor. Rotor/wake refinement geometry is unchanged.

Run each sensitivity case to numerical acceptance before comparing thrust and power. A relative difference above 5% in either metric triggers further investigation before an optimization campaign; differences below 5% support an exploratory search only, not physical qualification. All cases retain identical operating conditions and shape parameters. Record geometry equivalence explicitly because STL byte order can change across processes without changing the solid. A development comparison found maximum nearest-vertex distance below 8e-17 m between regenerated and archived baseline geometry, with equal volume.

The sensitivity configurations allow 1,000 new iterations after mapping a converged baseline using native `mapFields`. Mapping is an initial guess only; it does not transfer convergence status. Both consistent and explicit-patch mapping commands completed native compatibility tests in `runs/mapping-smoke`. Those utility tests used an unconverged archived field and are not sensitivity or performance results. All 32 regression tests pass.

`assess-baseline` adds a separate 500-sample settling check before promotion: both thrust and signed shaft power must be positive and satisfy the original relative drift/deviation tolerance over that longer window. This does not change the solver's original numerical acceptance policy. A synthetic regression fixture verifies that a stable last 100 samples cannot hide a longer settling trend.

Sensitivity initialization additionally compares the regenerated rotor with the archived STL using nearest-vertex distance, surface area and volume. The configured parameter match alone is insufficient if geometry code has changed. A regression test rejects a translated surface despite otherwise identical initialization metadata.

## Mesh-quality consistency correction

The first 1.25-resolution sensitivity mesh had 331,390 cells and passed the mesher's internal error report, but native `checkMesh` found two faces with global skewness up to 4.1489, exceeding the existing limit of 4. It was recorded as `MESH_FAILED` and never solved. The builder previously allowed boundary skewness up to 20 while the evaluator enforced a global maximum of 4. The boundary construction limit now uses the same configured value as the internal/global limit. This tightens generation rather than relaxing acceptance. The failed attempt is retained; the rebuilt mesh is a separate attempt.

## Native field round-trip correction

At iteration 3,000, the native solver wrote ASCII pressure values down to approximately 1e-323 in the far field. Its own scalar reader rejected subnormal tokens when `mapFields` reopened that solution. The archived run was left unchanged. Initialization now stages a complete source copy and replaces only nonzero ASCII values below the IEEE-754 minimum normal magnitude (2.2250738585072014e-308) with zero. The record retains original/staged hashes, replacement counts and maximum replaced magnitudes. The same audited copy-only repair applies to continuation fields. This representational correction is negligible relative to any physical/numerical tolerance; ordinary small values are preserved. Binary fields are never text-processed.

New cases write binary solution fields to avoid ASCII underflow on later mapping/restarts. `mapFields` failures are now classified as `INITIALIZATION_FAILED`, rather than mesh failures. Early development attempts retain their original recorded classification. The rebuilt fine mesh itself passed with maximum skewness 3.64964; its initial field-transfer attempt failed on the native ASCII read. The larger-domain attempt also reached field transfer before the same reader issue. These attempts remain in the archive.

The subnormal staging repair was verified natively on the accepted finer mesh: mapping succeeded after 390 pressure tokens were replaced, with maximum absolute replacement 3.732770915e-309. A four-iteration solve and a further four-iteration restart from binary output both completed. This is an I/O compatibility test, not a sensitivity performance result. Its files are retained in `runs/mapping-subnormal-smoke`.

The matched `study-report` checks identical operating/shape definitions, solver-equation dictionaries, native runtime fingerprints and acceptance limits. It requires distinct finer-background and larger-domain records, applies the same 500-sample force window, and withholds accepted percentage changes for failed or unconverged cases. Its output includes a decision record plus static PNG/SVG histories; the x axis is explicitly steady iteration, not physical time. Regression fixtures reject mismatched physics, duplicate records, incomplete convergence and changes above the pilot limit.

## Accepted continuation and matched settling control

The antisymmetric seed reached numerical `SUCCESS` at 3,000 iterations. Its final 100-sample means were 0.690865 N and 12.078022 W; all initial residuals were below 1e-4. The separate 500-sample check passed with thrust drift 1.934% and power drift 0.629%.

The larger-domain case then reached `SUCCESS` after 1,000 new iterations from that source. Comparing it directly to the old baseline confounds domain size with further settling. A same-mesh continuation from 3,000 to 4,000 was therefore run with identical equations, mesh and acceptance limits. Both passed the 500-sample check. Their 500-sample means were:

| Case | Thrust [N] | Shaft power [W] |
| --- | ---: | ---: |
| Original baseline at 3,000 | 0.701925794 | 12.138188701 |
| Same-mesh control at 4,000 | 0.645767380 | 11.879161531 |
| Larger domain, 1,000 new iterations | 0.656364839 | 11.915255598 |

Against the matched same-mesh control, domain differences are about +1.64% thrust and +0.30% power. This supports the pilot domain check, but also exposes slow settling in the original baseline; it is not physical qualification. No target or reference power has yet been frozen for optimization.

The finer mesh's first 1,000 iterations were rejected for residuals and force drift, and a continuation to 2,000 was started. Its apparent lower power was not treated as an improvement. The final study report uses accepted comparisons only and includes the settling control as the reference when supplied.

## CFD wall geometry audit

`surface-audit` extracts the native rotor patch to a new analysis directory and measures nearest-surface distances in both directions using every vertex and face centroid. It never edits the archived case. The maximum sampled CAD-to-wall distance was about 0.596 mm for the baseline and 0.328 mm for the finer mesh; 95th percentiles were approximately 0.141 and 0.129 mm. These are sampled geometry deviations, not Hausdorff bounds or flow-uncertainty estimates. The generated CAD STL remains one oriented, watertight solid.

The extracted STL surfaces have nonmanifold edges after coordinate merging/triangulation despite passing native polyMesh closure checks. Their oriented volume estimates are reported as estimates; export topology alone does not establish native fluid-mesh leakage. Near-edge/root geometric fidelity warrants further review before optimization. Audits and exported walls are retained under `runs/surface-audits/*-final`.

## Sustained convergence and continuation histories

At 2,500 iterations, the finer case passed the original final-point acceptance policy, with a maximum final residual of about 7.9e-5. A separate 500-sample check found earlier residuals up to 1.78e-4 in that same window. The saved `SUCCESS` record remains unchanged, but it was not promoted for optimization. A further 500-step continuation was launched with `force_window=500`.

The evaluator now applies residual acceptance throughout its configured force-averaging window, requiring complete evidence for all six fields at every sample. This prevents a transient dip from producing an objective. Source fingerprints separate this policy from historical final-point results. The promotion/study assessment uses a uniform 500-sample residual and force check across all cases, so legacy 100-sample records are rechecked consistently without rewriting them. The original baseline, same-mesh control and larger-domain case all pass this sustained 500-sample residual check (maximum initial residuals approximately 9.45e-5, 7.21e-5 and 7.18e-5).

Continuation may lengthen its averaging window but cannot shorten it or relax thresholds. Copied fields are normalized only for unreadable ASCII subnormals, and new output is binary. The write interval is derived to retain the final requested state; the native control run confirmed a saved 4,000-iteration checkpoint. The study plots join same-case continuation histories, replacing replayed iterations with the newer record and keeping mapped initialization histories separate. Regression fixtures cover both sustained residual rejection and correct restart-history joining. All 32 tests pass.


## Flow-oriented profile repair

The flow-oriented (`opposed_naca`) variant passed full solid validation, but its original coarse attempt failed at fluid-region point location and its refined 2,500-iteration attempt failed the sustained residual threshold. Neither supplies an accepted optimization objective. The coarse builder now chooses a far-field background-cell centre instead of an axis point on grid planes. A fresh coarse evaluation and an unchanged-solver refined continuation to 4,500 iterations were launched. All 37 tests pass, including generated seed checks on even/odd grids and custom domains. The complete acceptance and follow-up sequence is recorded in [CFD_FIX_PLAN.md](CFD_FIX_PLAN.md); final native outcomes remain in each immutable attempt record.

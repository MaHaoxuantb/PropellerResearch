# Current technical pathways

Last reviewed: 2026-10-23 (Japan time). Scope: the implementation in this checkout, its local configuration, and recorded research evidence. This is the maintained technical map; the specification and implementation plan describe intent and may contain superseded choices.

Update this file whenever geometry, CFD, acceptance, optimization, execution budgets, provenance, or qualification status changes. Keep implemented behavior separate from proposed work. Record completed outcomes from durable evidence, never from a launch message. Updating documentation alone does not require launching CFD or restarting campaigns.

## 1. End-to-end map

```text
Strict JSON configuration
  -> analytic loop and section construction
  -> cheap section screen
  -> triangulated lofts + hub Boolean union
  -> full solid validation
  -> native OpenCFD case construction
  -> mesh generation + native mesh acceptance
  -> steady SST/MRF solution
  -> sustained residual / force / conservation acceptance
  -> immutable attempt result + evidence
       -> assessment / matched sensitivity / surface audit
       -> baseline freezing -> simulated annealing
       -> fixed-radius campaign -> initial sampling -> Bayesian acquisition
       -> reports / comparisons / recovery
```

Continuation reuses an archived mesh and saved fields in a new attempt. Mapped initialization builds an independent mesh and transfers a converged source flow guess. Both return to numerical acceptance; neither inherits acceptance merely by copying fields.

All performance remains research evidence with `qualified: false`. Numerical `SUCCESS`, optimization feasibility, mesh/domain sensitivity, and experimental qualification are separate gates.

## 2. Entry points and configuration

- `main.py`, `src/propeller/__main__.py`, and the installed `propeller` command dispatch to [cli.py](../src/propeller/cli.py).
- [config.py](../src/propeller/config.py) enforces schema version 1, exact keys, finite values, bounds and supported settings. Configuration is JSON in SI units; angles are radians.
- Required groups: `design`, `parameters`, `bounds`, `operating`, `geometry`, `cfd`, `optimization`, plus `schema_version`.
- Design fixes diameter envelope, loop count, hub/root/shaft dimensions. Operating conditions fix RPM, handedness, inflow and fluid properties. Target thrust/reference power are initially nullable and must be positive before annealing.
- The native backend supports OpenCFD v2512. Configuration cannot assert qualification. The estimated tip Mach limit is 0.3; there is no compressible pathway.
- Base environment: `.venv`, pinned by [requirements.lock](../requirements.lock). Bayesian environment: `.venv-bo`, pinned by [requirements-bo.lock](../requirements-bo.lock); its Torch/BoTorch dependencies are separate from the core package.
- `Model.py` is a historical prototype, not the production entry point.

Commands below run from the repository root. Path placeholders denote existing records or new output locations, not literal filenames.

## 3. Geometry and inspection

Implementation: [core.py](../src/propeller/geometry/core.py), [screening.py](../src/propeller/geometry/screening.py), [validation.py](../src/propeller/geometry/validation.py), [preview.py](../src/propeller/geometry/preview.py).

The nine controls are `H`, `beta`, `shape_exponent`, `b0`, `b1`, `b2`, `Z`, `alpha0`, `alpha1`. For loop coordinate u in [0,1]:

```text
r(u) = root_radius + H sin(pi u)^shape_exponent
theta(u) = beta (2u - 1)
C(u) = (r cos(theta), r sin(theta), Z sin(pi u)^2)
alpha(u) = alpha0 + alpha1 (2u - 1)
w(s) = b0 (1-s)^2 + 2 b1 s(1-s) + b2 s^2
```

Here s is normalized projected arc length. Analytic derivatives supply tangents and curvature; curvature-weighted sampling resolves roots and the apex. Exponents in (1,2) are rejected because of singular root curvature; exponent 1 is separately supported.

Sections use a symmetric NACA thickness law with finite trailing edge. `fixed_naca` retains the same profile orientation; `opposed_naca` smoothly exchanges profile ends around the apex to address opposite chordwise flow on the two legs. Pitch uses u, while chord width uses s. Root chords retract through a quintic transition into the hub. This is a geometric connection, not a stress-qualified fillet.

Explicit triangulated lofts are rotationally replicated and united with the hub using Manifold. STL/PLY and sampled section arrays are available; STEP/CAD export is not implemented.

Validation stages:

1. Cheap section screen checks finite points, diameter/root enclosure, nominal body thickness, local folds and Mach. `geometry sample` has its own narrower diameter/root envelope screen.
2. Full validation checks closed, consistently oriented, positive-volume, connected solids; degenerate faces; diameter, thickness and root area; hub overlap; inter-loop collision/clearance; local folds; and triangle crossings.
3. `--fast` omits the triangle-crossing scan and marks that omission. CFD always requests full validation.

```sh
.venv/bin/python main.py geometry preview --output runs/baseline
.venv/bin/python main.py geometry validate --output runs/validation
.venv/bin/python main.py geometry export --output runs/exports
.venv/bin/python main.py geometry sample --count 250
```

Preview/validate/export share the artifact-writing branch: configuration, measurements, validation, STL, PLY, compressed sections and standalone interactive HTML. The inspector shows section geometry, estimated local flow direction and curvature; it does not predict aerodynamic performance.

## 4. Native CFD pathway

Implementation: [case_builder.py](../src/propeller/simulation/case_builder.py), [runner.py](../src/propeller/simulation/runner.py), [evaluator.py](../src/propeller/simulation/evaluator.py).

The configured Homebrew launcher is `/Applications/OpenFOAM-v2512.app/Contents/Resources/etc/openfoam`. Runtime discovery verifies distribution/version and solver binary fingerprints before cache lookup. The pathway runs natively on macOS, without Docker.

The generated model is a full-rotor, steady incompressible SST/MRF case with an explicit rotation zone, static/open-flow boundaries or configured inflow, Spalding viscosity wall treatment, refinements and layers. Native command order:

```text
blockMesh -> surfaceFeatureExtract -> snappyHexMesh -overwrite
 -> checkMesh -allTopology -meshQuality -> topoSet
 -> [mapFields if initialized] -> simpleFoam
```

Optional MPI inserts `decomposePar` and `mpirun ... simpleFoam -parallel`. Serial execution is verified; the installed MPI startup stalled in development, so the practical configuration remains `mpi_ranks: 1`.

`checkMesh` must report `Mesh OK.` and satisfy configured cell-count, nonorthogonality and skewness limits. Layer coverage and later y+ are retained for review. Passing these limits does not establish grid independence or adequate near-wall physics.

The fluid-region seed is a far-field background-cell centre, avoiding the former seed on grid planes. New solution fields are binary. Native logs and command timing/failure evidence remain in each case. Per-command timeouts, campaign candidate deadlines and supervisor-liveness checks can independently stop execution; current budget conflicts are described below.

```sh
.venv/bin/python main.py doctor
.venv/bin/python main.py --config config/flow-oriented-fine-4500.json case --output runs/new-case
.venv/bin/python main.py --config config/flow-oriented-fine-4500.json evaluate --output runs/evaluations
```

## 5. Metrics, numerical acceptance and failures

Implementation: [parser.py](../src/propeller/simulation/parser.py), [analysis.py](../src/propeller/analysis.py).

Paired force/moment histories require finite values and matching, strictly increasing timestamps. Thrust is total axial force; signed shaft power is `-fluid_torque_z * signed_angular_velocity`. Metrics include CT, CP, advance ratio, forward-flight efficiency where defined, and static figure of merit where defined.

Acceptance requires a clean solver ending, enough force samples, complete initial residual evidence for p/Ux/Uy/Uz/k/omega throughout the averaging window, residuals below the configured threshold, final local continuity below its limit, rotor y+ evidence, and both relative standard deviation and half-window drift within force/power tolerances. The current refined baseline uses 500 samples, residual tolerance 1e-4 and force/power tolerance 2%. These are configured numerical limits, not physical error bounds. Steady iterations are not physical time.

Numerical `SUCCESS` alone does not require positive thrust. Baseline promotion requires positive thrust and power; optimization applies its own scoring/observation rules.

| Outcome | Meaning / downstream behavior |
| --- | --- |
| `RUNTIME_UNAVAILABLE` | Native backend unavailable; evaluator returns before starting an attempt. |
| `GEOMETRY_INVALID` | Section, solid or fixed-radius checks reject the candidate; no accepted flow objective. |
| `EXPORT_FAILED`, `MESH_FAILED` | Case construction or native mesh stage failed. |
| `INITIALIZATION_FAILED` | Field transfer failed or its interpolation evidence is incomplete. |
| `SOLVER_FAILED`, `TIMEOUT` | Native command failure or enforced time limit. |
| `NOT_CONVERGED` | Solver evidence fails acceptance, including incomplete/unusable histories. |
| `INTERRUPTED`, `INTERNAL_ERROR` | Interrupted work or unexpected evaluator exception. |
| `SUCCESS` | Numerical criteria pass; may enter cache and applicable optimizer checks. |

Campaign acknowledgements can also report `WORKER_FAILED` or `WORKER_LOST`. Failed/unconverged force estimates may be examined diagnostically but never become fabricated performance observations.

## 6. Restart and mapped-initialization branches

Implementation: [continuation.py](../src/propeller/simulation/continuation.py), [initialization.py](../src/propeller/simulation/initialization.py), [field_io.py](../src/propeller/simulation/field_io.py).

**Continuation:** a completed serial `SUCCESS` or `NOT_CONVERGED` attempt must have accepted mesh evidence, matching configuration/geometry, and saved U/p/k/omega/nut fields. The native runtime must match. A new attempt copies inherited inputs and records hashes. Only end iterations, timeout and a non-shortened force window may change; enough new iterations must cover a full window. Changing physics, geometry, mesh or thresholds requires a new evaluation. Timeout results are not directly accepted by this continuation entry point.

**Mapped initialization:** a numerically converged source with identical shape and operating conditions can initialize a newly generated mesh/domain study. The source geometry/runtime and input evidence are checked. `mapFields` must record interpolation of all five fields and a clean ending. The new solution must independently converge.

Older ASCII fields containing nonzero subnormal values may be unreadable by the native parser. Only staging copies are normalized below the minimum normal double, retaining original/staged hashes, counts and magnitudes. Binary fields are not text-processed.

```sh
.venv/bin/python main.py --config config/antisymmetric-extended.json continue-evaluation runs/evaluations/KEY/ATTEMPT
.venv/bin/python main.py --config config/antisymmetric-fine.json evaluate --initialize-from runs/evaluations/KEY/ATTEMPT
```

Configuration examples must match the actual parent contract; the commands do not authorize arbitrary parent/configuration combinations.

## 7. Qualification and baseline promotion

Implementation: [study.py](../src/propeller/study.py), [baseline.py](../src/propeller/baseline.py), [surface_audit.py](../src/propeller/simulation/surface_audit.py).

1. `assess-baseline` rechecks longer-window settling/residual evidence and positive propulsion metrics; readiness is `READY_FOR_SENSITIVITY`.
2. `study-report` checks matched shape/physics, dictionaries, runtime and acceptance limits. It requires distinct finer-mesh and larger-domain records. A same-mesh settling control can replace the original reference to separate continuing settling from sensitivity. Accepted comparisons within the pilot limit (default 5%) can yield `READY_FOR_EXPLORATORY_SEARCH`; unconverged records do not provide accepted percentage changes.
3. `surface-audit` extracts the native rotor wall into a separate analysis directory and samples bidirectional nearest-surface deviations using vertices and face centroids. These are sampled distances, not Hausdorff bounds. Extracted STL topology alone does not prove native fluid-mesh leakage.
4. `freeze-baseline` requires exact evaluated configuration, sustained accepted positive thrust/power, and writes a new configuration plus `.baseline.json` provenance. `--study` additionally enforces that the study selects this result; that argument is optional in code, so freezing alone does not prove sensitivity qualification. Optional timeout calibration uses measured iteration cost with a 50% buffer.
5. Reusing a frozen annealing seed checks result/study/contract hashes, current runtime, geometry, mesh-generation ancestor template and current acceptance. Otherwise use an explicitly fresh evaluation.

```sh
.venv/bin/python main.py assess-baseline PATH/result.json --window 500
.venv/bin/python main.py study-report BASE/result.json FINE/result.json DOMAIN/result.json CONTROL/result.json
.venv/bin/python main.py surface-audit PATH/ATTEMPT --output runs/surface-audits/new-audit
.venv/bin/python main.py --config PATH/case/configuration.json freeze-baseline PATH/result.json --output config/study.json --study PATH/study.json
```

Mesh/domain/transition-model sensitivity and experimental validation remain necessary before physical performance claims. The campaign baseline checker validates numerical/geometry/mesh compatibility; it does not require a completed sensitivity-study record.

## 8. Optimization pathway A: simulated annealing

Implementation: [annealing.py](../src/propeller/optimization/annealing.py).

All nine normalized coordinates are free inside bounds. Proposals reflect at bounds; 75% perturb one coordinate and the remainder perturb all coordinates. Accepted finite positive-power `SUCCESS` results receive:

```text
score = P/reference_power + thrust_penalty * max(0, 1 - T/target_thrust)^2
feasible = T >= target_thrust
```

Metropolis acceptance permits uphill moves according to temperature. The separate best-feasible archive selects minimum power among thrust-feasible results. Failed results receive no score. Temperature cools after scored proposals. Proposal/evaluation/time limits bound the chain.

Pending proposals are saved before expensive work. Checkpoints retain current state, best point, RNG, temperature, history and elapsed time. Resume rejects changed source/configuration. Distinct restarts use separate seeds/checkpoints. The CLI requires `--exploratory`.

```sh
.venv/bin/python main.py --config config/study.json optimize --exploratory --checkpoint runs/study/checkpoint.json
.venv/bin/python main.py --config config/study.json resume --exploratory --checkpoint runs/study/checkpoint.json
```

## 9. Optimization pathway B: fixed-radius Pareto pilot

Implementation: [bayesian.py](../src/propeller/optimization/bayesian.py), [campaign.py](../src/propeller/optimization/campaign.py), [mobo-pilot.json](../config/mobo-pilot.json).

Eight normalized coordinates are proposed; Brent root solving determines H so sampled section radius matches the fixed radius. Full finished-solid radius is independently checked to 5e-5 m. The configured radius is approximately 0.111760 m. Each candidate passes cheap screening, full solid checks, native meshing and CFD.

The candidate pool mixes 25% uniform proposals with Gaussian proposals around the baseline (normalized standard deviation 0.16), reflecting at bounds. Exact/near duplicate coordinates are excluded. Initial selection maximizes distance from accepted/pending observations; with fewer than three observations it prefers points within normalized distance 0.5 of baseline. That bound permits large simultaneous changes.

Actual qLogNEHVI starts only after at least six counted CFD evaluations and four accepted observations including baseline. Independent standardized GPs model normalized thrust and log normalized power. Sampled outputs maximize thrust and thrust per watt, normalized to baseline, using a zero-reference frontier and a positive-thrust acquisition constraint. Pending candidates enter acquisition. Accepted finite positive-power negative-thrust results remain valid poor observations; useful positive propulsion/hypervolume excludes them. Failed solves never train the performance model. Proposal errors fall back to space filling; two consecutive proposal failures stop new work.

Current configured budgets:

| Control | Value / behavior |
| --- | --- |
| Whole campaign / candidate | 12 hours / 3 hours; native per-command timeout still applies. |
| Evaluations / selected attempts | 12 new counted evaluations / 36 jobs including invalid solids. |
| Candidate pool | 128 section-feasible points requested per proposal. |
| Workers / memory plan | Starts with at most 2; may grow to 4 after 2 new accepted observations; 12 GB total, 3 GB/job floor, observed child memory ×1.5 can reduce capacity. This is planning, not an OS quota. |
| Failure stop | 4 consecutive numerical failures; active jobs ordinarily drain. |
| Stagnation stop | 4 accepted results without >0.5% area gain, after at least 8 accepted observations. |
| Disk | 8 GB free reserve and 8 GB campaign-growth cap. |

The supervisor has an exclusive lock, durable RNG/state and owned workers. Resume requires unchanged source/base and exited prior workers; the original deadline persists. STOP/deadline/disk-reserve stops terminate owned work; ordinary failure/budget/stagnation stops drain it. A full candidate allowance must remain before another launch.

```sh
PYTHONPATH=src .venv-bo/bin/python -m propeller.optimization.campaign start runs/NEW-CAMPAIGN --settings config/mobo-pilot.json
PYTHONPATH=src .venv-bo/bin/python -m propeller.optimization.campaign status runs/mobo-pilot-02
PYTHONPATH=src .venv-bo/bin/python -m propeller.optimization.campaign stop runs/EXISTING-CAMPAIGN
PYTHONPATH=src .venv-bo/bin/python -m propeller.optimization.campaign run runs/EXISTING-CAMPAIGN
```

`start` requires an empty destination and launches a detached supervisor. `run` is explicit resume, not a way to bypass fingerprints or expired deadlines. State, per-job input/acknowledgement/logs, `pareto.json` and refreshing `report.html` persist.

## 10. Provenance, scheduling and reporting

Implementation: [experiments.py](../src/propeller/experiments.py), [scheduler.py](../src/propeller/scheduler.py), [progress.py](../src/propeller/progress.py), [recovery.py](../src/propeller/simulation/recovery.py), [reporting.py](../src/propeller/reporting.py).

Evaluation identity includes exact configuration, package Python source fingerprint, recorded dependency versions, actual runtime, and optional initialization/radius constraints. Nearby float vectors are not merged. Only successful records are cached. Attempts live under `runs/evaluations/KEY/ATTEMPT` (or a campaign-local evaluations root), indexed by SQLite with WAL.

New attempts preserve manifest, source snapshot, process state and available screening/validation/mesh/geometry evidence, case dictionaries, commands, native logs, histories and final results. Atomic JSON replacement protects state writes. Failed early stages naturally lack downstream artifacts. Early development records may contain fingerprints without source snapshots. Solution write/purge policy limits retained fields; cases are not deleted after execution.

`batch` runs independent configurations within CPU and explicit concurrency-memory budgets; one annealing chain remains serial. Concurrent throughput is not equivalent to physical qualification.

`recover` checks that the recorded owner is absent and that no live native case child owns the solve. Complete orphaned native evidence can publish an audited result; otherwise the attempt becomes interrupted. Recovery retains evidence and does not invent a known original process exit code or overwrite an existing published result.

```sh
.venv/bin/python main.py status --output runs/evaluations
.venv/bin/python main.py recover --output runs/evaluations
.venv/bin/python main.py report PATH/result.json --output runs/report.html
.venv/bin/python main.py compare COARSE/result.json FINE/result.json
.venv/bin/python main.py batch config/demo.json config/antisymmetric.json --workers 1
```

HTML evaluation reports produce force/power/residual plots, PNG/SVG and CSV evidence. Non-HTML `report` output is Markdown. General comparisons are descriptive; matched study reports provide the stronger compatibility/acceptance gate. Study histories join same-case continuations without treating mapped-source iterations as a continuation of the new mesh.

## 11. Current durable evidence and unresolved pathway work

Reviewed against local records on 2026-10-23; this section is a snapshot, not live process monitoring.

- Flow-oriented refined baseline: [result](../runs/evaluations/6388c6b553280fec/ec0cf993e6f34d22b482e80e207476a5/result.json) records `SUCCESS`, approximately 0.683206 N and 9.361454 W after the 4,500-iteration continuation. It remains physically unqualified.
- [Pilot 02 state](../runs/mobo-pilot-02/state.json) records `STOPPED` for repeated numerical failures, with 8 counted CFD attempts. The [investigation](research/pilot-diagnostics/INVESTIGATION.md) records 1 accepted new CFD candidate, 5 unconverged, 2 timed out, and 1 additional invalid-solid job. All selected proposals were initial sampling; zero Bayesian acquisition proposals ran. The accepted new candidate was worse than baseline in both objectives. No optimized design was established.
- Earlier antisymmetric mesh/domain/settling evidence and surface deviations are retained in [implementation notes](IMPLEMENTATION_NOTES.md). They do not automatically qualify the different flow-oriented geometry.
- [CFD repair plan](CFD_FIX_PLAN.md) contains historical launched/pending statements; durable results and the newer investigation take precedence for completed outcomes.

Proposed next work, **not implemented by this document**:

1. Reconcile the 3-hour candidate allowance with inherited 5,400-second native command timeouts, and report the actual limiting budget.
2. Bootstrap with smaller parameter-scaled local changes; contract exploration after repeated failures. Current sampling lacks a failure-feasibility model or trust-region contraction.
3. Use bounded standalone continuation experiments on a borderline residual case and a still-drifting case, preserving original attempts and thresholds. The campaign currently starts candidates cold and does not orchestrate continuation.
4. Build enough accepted observations to reach real acquisition; report proposal method and distinct outcome counts explicitly.
5. Complete matched sensitivity and physical validation before promoting any claimed improvement.

Detailed numerical evidence, original result hashes and recommended repair order: [pilot investigation](research/pilot-diagnostics/INVESTIGATION.md) and [evidence.json](research/pilot-diagnostics/evidence.json). These diagnostic files existed before this technical-map update.

## 12. Verification and maintenance checklist

```sh
PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -v
PYTHONPATH=src .venv-bo/bin/python -m unittest discover -s tests -v
```

The second environment supports the real Bayesian acquisition smoke test. Tests cover geometry derivatives/orientation/constraints, exports and crossings, case generation, signs/convergence/failure handling, cache/restart contracts, deterministic annealing resume and Bayesian/campaign behavior. Synthetic optimization evidence is separate from CFD records. Historical README reports 48 passing tests in the Bayesian environment; that is historical evidence, not a fresh test run for this documentation change.

For future updates: inspect changed modules and configuration; verify affected commands and gates; inspect durable results for status claims; update this map and relevant focused notes together. Preserve unresolved limitations until evidence closes them. Do not silently reinterpret old acceptance policies or rewrite archived outcomes.

Related references: [README](../README.md), [specification](SPECIFICATION.md), [implementation plan](../IMPLEMENTATION_PLAN.md), [pilot implementation](research/PILOT_IMPLEMENTATION.md), [one-week research plan](research/ONE_WEEK_OPTIMIZATION.md).

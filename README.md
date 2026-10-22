# Toroidal propeller research

A nine-parameter loop-rotor generator, offline 3D inspector, native Homebrew OpenFOAM evaluator, and checkpointable simulated annealing engine.

The current system is a research implementation. Geometry and optimizer tests are separate from CFD qualification; it does not claim to have found a validated optimal propeller.

## Run from this checkout

The existing `.venv` has the dependencies installed. Commands run from the repository root:

```sh
.venv/bin/python main.py doctor
.venv/bin/python main.py geometry preview --output runs/baseline
.venv/bin/python main.py geometry sample --count 250
.venv/bin/python main.py case --output runs/my-case
.venv/bin/python main.py evaluate --output runs/evaluations
```

To extend a completed serial solve without rebuilding the mesh or overwriting its record, increase `cfd.iterations` (and optionally its timeout) in a new configuration. The averaging window may also be lengthened, never shortened:

```sh
.venv/bin/python main.py --config config/antisymmetric-extended.json continue-evaluation runs/evaluations/KEY/ATTEMPT
.venv/bin/python main.py assess-baseline runs/evaluations/KEY/NEW_ATTEMPT/result.json --window 500
```

The continuation copies the latest saved fields into a new attempt, verifies the geometry and native runtime checksums, and records hashes of all inherited case inputs. Changed physics, geometry, mesh settings or convergence thresholds require a new evaluation.

For a mesh or domain study, an accepted baseline can supply the initial flow guess:

```sh
.venv/bin/python main.py --config config/antisymmetric-fine.json evaluate --initialize-from runs/evaluations/KEY/ATTEMPT
```

This requires a numerically converged source with identical shape and operating conditions. The new mesh is independently generated and checked, and its solution must independently converge. Field interpolation and its source hashes are recorded; a mapped guess is not an accepted performance result. New solution fields use binary output. When an older ASCII source contains unreadable subnormal values, a staging copy zeros only values below the minimum normal double; the original files and the exact repair audit are retained.

Open `runs/baseline/preview.html` in a browser. It is a standalone file with no CDN or server dependency. Drag to orbit, scroll to zoom, choose top/side views, and move the section slider to inspect chord width, curvature and the estimated local flow direction. STL, PLY, configuration, sampled sections and validation results accompany the viewer.

`geometry sample` is a cheap envelope/root screen, **not** certification that sampled solids are valid. `geometry preview --fast` skips the triangle-crossing scan and is explicitly marked incomplete. CFD always invokes the full validator.

To create a new environment, use Python 3.9–3.12 with the tested pinned dependencies:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
```

The checkout entry point requires no editable installation. A conventional package and `propeller` entry point are also defined in `pyproject.toml`.

## Native Homebrew OpenFOAM

The configured launcher is:

```text
/Applications/OpenFOAM-v2512.app/Contents/Resources/etc/openfoam
```

This is the already-installed Homebrew cask application. Its launcher mounts the bundled read-only runtime and activates the native environment. `doctor` verifies the actual distribution/version and hashes the solver binaries. No Docker dependency is used. Homebrew itself currently reports the third-party cask tap as untrusted; this tool does not alter Homebrew's trust settings or reinstall the application.

Change `cfd.launcher` if this installation moves. Case dictionaries currently target **OpenCFD v2512**; another distribution/version is rejected rather than silently treated as compatible.

The generated case uses a full-rotor, steady incompressible SST/MRF model, explicit rotation zone, open static-flow boundaries, Spalding viscosity wall treatment, automatic surface/region refinement and layers. It runs `blockMesh`, `surfaceFeatureExtract`, `snappyHexMesh`, `checkMesh -allTopology -meshQuality`, `topoSet`, and `simpleFoam`. Native logs remain beside the case. Estimated tip Mach above 0.3 stops case generation; a compressible backend is not yet implemented.

Serial native execution is verified. The optional `mpi_ranks` path is implemented, but the installed MPI runtime stalled during startup on this machine; keep `mpi_ranks: 1`. Startup timeouts preserve the failure logs.

An evaluation can complete numerically without being physically qualified. Inspect mesh quality, y+, force drift, residuals and conservation. Grid/domain/transition-model sensitivity and experimental validation remain required before performance claims.

## Configure a design

`config/demo.json` is an illustrative 254 mm maximum-diameter, 5,000 RPM static-air experiment. It does not encode manufacturing requirements. Dimensions are meters and angles radians. The nine controls are:

```text
H, beta, shape_exponent, b0, b1, b2, Z, alpha0, alpha1
```

The optimizer changes those variables within the configured bounds. The operating point and constraints stay fixed for each experiment. The continuous exponent range starts at 2 to avoid singular root curvature. Both signs of pitch are allowed to explore handedness. Initial airfoil geometry is fixed NACA 0012 with a finite trailing edge; the reference edge is geometrically defined, and its leading-edge interpretation is checked diagnostically rather than assumed correct for both legs.

`config/antisymmetric.json` is a second test hypothesis: pitch changes sign along the loop (alpha0 = 0, alpha1 = -15 degrees). It uses the same nine controls and is not an established aerodynamic improvement.

The root chord retracts smoothly into the hub over a fixed end interval. This creates finite overlap without smoothing away the candidate geometry. It is not a stress-qualified fillet. Manifold unions explicit triangulated section lofts; STEP/CAD output is not implemented.

## Optimize and resume

After a baseline returns `SUCCESS` with positive thrust/power, passes `assess-baseline`, and the matched pilot sensitivity study is acceptable, freeze its target and power normalization into a new configuration:

```sh
.venv/bin/python main.py --config runs/evaluations/KEY/ATTEMPT/case/configuration.json freeze-baseline runs/evaluations/KEY/ATTEMPT/result.json --output config/study.json
.venv/bin/python main.py --config config/study.json optimize --exploratory --checkpoint runs/study/checkpoint.json
.venv/bin/python main.py --config config/study.json resume --exploratory --checkpoint runs/study/checkpoint.json
```

Replace `KEY/ATTEMPT` with the actual returned evaluation path. `freeze-baseline` requires the configuration to match the evaluated manifest. To set application-specific requirements, edit a separate configuration with positive `target_thrust` and `reference_power`.

`--exploratory` is required because the numerical model is not yet qualified. Best candidates are **unconfirmed** until higher-fidelity checks establish feasibility and improvement. Failed or unconverged evaluations never get an objective. A candidate with low power but inadequate thrust cannot enter the best-feasible archive.

CFD acceptance requires residuals to remain below the configured threshold throughout the force-averaging window. The current fine-study template uses 500 samples; a single low final residual is insufficient.

The checkpoint stores the current state, best feasible point, RNG state, temperature, proposals and acceptance history. It checkpoints pending proposals before expensive work. Resume rejects changed configuration or source code. For independent restarts, use different configured seeds and separate checkpoints. A single annealing chain runs serially. The `batch` command supports independent candidates with CPU limits and explicit memory budgets for concurrent jobs. Concurrent CFD throughput has not yet been qualified on this machine.

## Records and reports

Every attempted CFD evaluation has its own directory and SQLite index. New attempts also retain a package source snapshot; early development attempts retain source fingerprints only. Inputs, source/dependency/runtime fingerprints, geometry, validation, commands, force histories, mesh checks, diagnostics and failures are retained. Only successful evaluations are cacheable. Nearby float vectors are not merged. Large solution fields use the explicit OpenFOAM write/purge policy; the case is not deleted after execution.

```sh
.venv/bin/python main.py report runs/evaluations/KEY/ATTEMPT/result.json --output runs/report.html
.venv/bin/python main.py compare path/to/coarse/result.json path/to/fine/result.json
.venv/bin/python main.py study-report path/to/baseline/result.json path/to/fine/result.json path/to/larger-domain/result.json path/to/same-mesh-control/result.json
.venv/bin/python main.py surface-audit runs/evaluations/KEY/ATTEMPT --output runs/surface-audits/my-audit
.venv/bin/python main.py status --output runs/evaluations
.venv/bin/python main.py recover --output runs/evaluations
.venv/bin/python main.py batch config/demo.json config/antisymmetric.json --workers 1
```

HTML reports include thrust, signed shaft power and residual plots, with PNG/SVG figures and CSV histories alongside them. `recover` marks attempts interrupted only when their recorded owner process no longer exists, preserving all files.

`study-report` checks matched physics, solver dictionaries, runtime and acceptance limits, then assesses the same longer force window in each case. A same-mesh continuation can serve as the comparison reference to separate further settling from mesh/domain changes. It requires distinct finer-mesh and larger-domain configurations, and does not compute accepted percentage changes from unconverged results. A 5% pilot sensitivity screen can permit an exploratory search; it does not qualify the physical model. General `compare` output remains descriptive.

## Current verification result

All 37 tests pass. The legacy antisymmetric baseline continuation reached `SUCCESS` at 3,000 iterations and passed a separate 500-iteration settling check. The newer flow-oriented profile remains under evaluation: its coarse seed-placement failure has been repaired and a fresh evaluation launched; its refined solution is being continued from 2,500 to 4,500 iterations with unchanged acceptance thresholds. See the [CFD repair plan](docs/CFD_FIX_PLAN.md) for the active record paths and qualification sequence. No optimized or physically qualified design is claimed. Reports are in `runs/reports/`.

## Verify the code

```sh
PYTHONPATH=src .venv/bin/python -m unittest discover -s tests -v
```

Tests cover analytic derivatives, section mapping, frame orientation, constraints, triangle crossings, connected solids, deterministic exports, cache identity, failure handling, torque signs and deterministic optimization resume. Synthetic optimizer tests are labelled and never written to the CFD results archive.

See [the source specification](docs/SPECIFICATION.md), [implementation plan](IMPLEMENTATION_PLAN.md), and [implementation notes](docs/IMPLEMENTATION_NOTES.md). The original `Model.py` remains as a historical prototype.

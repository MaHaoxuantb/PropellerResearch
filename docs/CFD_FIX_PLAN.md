# Flow-oriented rotor CFD repair

The immediate goal is to obtain independently accepted coarse and refined
evaluations of the same flow-oriented rotor. Optimization remains gated on
convergence and mesh/domain sensitivity; a solver exit code alone is insufficient.

## 1. Repair fluid-region seed placement

The original coarse attempt stopped in `snappyHexMesh` before refinement:
`Point (0.6096 0 0) is not inside the mesh or on a face or edge.`
That point lies on background-grid planes for the 20 × 20 × 32 grid.

Replace the fixed axis point with the centre of the background cell at the
positive-x/y/z far-field corner. The centre stays strictly inside the box,
off its grid planes, and outside the diameter-constrained rotor for all
supported domain sizes and cell counts. Preserve the physical domain, mesh
resolution, CAD and mesh-quality limits. Check the generated dictionary for
even, odd, expanded-domain and minimum-resolution backgrounds, then require
native `snappyHexMesh` and `checkMesh` to pass on a new coarse attempt.

## 2. Test whether the refined flow needs more settling

The original refined mesh passed quality checks with 332,191 cells. Its
2,500-iteration solution failed sustained residual acceptance, although its
last 500 force samples passed the 2% stability check. Maximum axial-velocity
initial residuals fell from 1.82e-3 over iterations 501–1000 to 3.72e-4 over
2001–2500. This supports testing additional settling; it does not prove
eventual convergence.

Continue from the binary 2,500 checkpoint to 4,500 using
`config/flow-oriented-fine-4500.json`. Copy the archived mesh, dictionaries and
fields into a new recorded attempt, verify checksums/runtime, and retain the
original result. Change only the end iteration; keep the original 5,400-second
timeout and every acceptance threshold. Run the repaired coarse case from
its own initial fields at the original 2,500-iteration budget.

## 3. Evaluate results before further solver changes

Each run must finish cleanly, pass native mesh quality, provide all six
residual histories below 1e-4 throughout its final 500 samples, satisfy the
existing continuity check, record rotor y+, and pass the unchanged 2% force
and power stability criteria. Positive thrust and shaft power are required
before baseline promotion. Acceptance remains numerical, not physical validation.

If either case remains rejected, inspect residual and force trends by
500-iteration windows. Extend only where continuing decay supports it. A
plateau, persistent oscillation or growing residuals requires investigating
mesh/wall resolution and steady-MRF applicability rather than relaxing
thresholds or repeatedly extending the same calculation. Any change to
discretization or turbulence treatment must be a separate recorded experiment
and applied consistently before mesh comparisons.

## 4. Return to the optimization qualification sequence

Compare accepted coarse/refined thrust and power at identical operating
conditions. Differences above the existing 5% pilot limit require another
grid level. Check a larger domain and same-mesh settling control, then freeze
a baseline and use accepted results to seed an exploratory optimization.
No unconverged result supplies an optimization score.

## Verification and active records

All 37 unit tests passed after the seed repair, including the new generated-case
regression. Native reruns were launched using the installed Homebrew OpenFOAM
v2512 launcher:

- Coarse: `runs/evaluations/1f042f76336086a7/11dacc538c2c4e8c98713f2fc9f7b5bd`.
- Refined continuation: `runs/evaluations/6388c6b553280fec/ec0cf993e6f34d22b482e80e207476a5`.

Inspect current durable state and latest solver iteration with
`.venv/bin/python main.py status`. The record paths above describe launched
work; final acceptance is determined by their `result.json` files.

The repaired coarse case subsequently completed native mesh generation and
passed `checkMesh` (`Mesh OK.`). Both cases were observed advancing in
`simpleFoam`; the coarse meshing failure is resolved. Refined convergence and
the final coarse solution remain pending.

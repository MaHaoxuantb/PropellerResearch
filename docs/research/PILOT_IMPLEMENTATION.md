# Implemented pilot

The research recommendation is now implemented as a bounded pilot rather than a
mandatory seven-day campaign. Configuration: `config/mobo-pilot.json`.

- Radius: 111.759997904 mm; generator uses nine parameters with eight independent
  coordinates and a bounded root solve for H. The completed solid is rechecked.
- Objectives: maximize thrust and thrust per watt; baseline is a guide rather
  than a numerical minimum-thrust constraint.
- Operating condition: existing 5,000 RPM static-air, two-loop setup.
- Proposal methods: initial feasible space-filling selection followed by BoTorch
  qLogNEHVI on thrust and log-power Gaussian processes, transformed to the two
  objectives using shared posterior samples. Pending jobs influence acquisition.
- Cheap rejection: impossible radius, exposed root caps, local folds, insufficient
  nominal body thickness, non-finite coordinates and unsupported Mach regime.
- Authoritative pre-CFD checks: watertightness, connectivity, orientation, positive
  volume, finished radius/diameter, degenerate triangles, minimum thickness/root
  area, finite hub overlap, outside-hub clearance/collisions and triangle crossings.
- Scope limits: no structural solver or trusted pre-CFD aerodynamic model. The
  nominal flow indicators are geometric diagnostics, not an airfoil stall model.
- Pilot limits: 12 new evaluations, 36 total attempts, 12 hours, three hours per
  candidate, six initial CFD attempts before Bayesian selection, minimum four
  accepted observations to fit the Bayesian model.
- Stops: four consecutive numerical failures; four accepted observations without
  meaningful frontier-area improvement once eight accepted observations exist;
  two proposal failures; disk reserve/growth limits; elapsed time or evaluation
  cap; user STOP. Geometry rejection does not invent an aerodynamic score.
- Execution: detached supervisor and individually acknowledged workers; two jobs
  initially, up to four after observed successful runs under memory estimates.
  Whole-campaign and per-candidate deadlines reach the native process loop.
- Artifacts: source/config/runtime provenance for each native attempt, campaign
  state, proposal method, screened rejection counts, measured Pareto points and
  an HTML progress report. Restart retains the original absolute deadline.

Tests include a real qLogNEHVI proposal on a synthetic fixture; the fixture is
never added to the CFD archive. Integration checks establish that rejected
geometry cannot invoke native meshing/solving, deadline expiration terminates an
owned command, fixed-radius mapping preserves the final solid and hub, and a
stop request launches no work. Acceptance still requires the full sustained
residual and force checks used by the existing CFD parser.

## Live pilot audit (2026-10-22)

The active campaign is `runs/mobo-pilot-02`. The initial trial's two CFD attempts
and elapsed time are deducted: the continuation permits ten more CFD attempts
within the original twelve-hour ceiling. Both first candidates passed the full
geometry and native mesh checks and entered simpleFoam. No new accepted result
was available at startup verification.

Initial proposals stay near the baseline while fewer than three accepted
observations exist. Converged negative-thrust results with positive power remain
valid model observations, but do not enter the positive-thrust Pareto frontier.

A launch-time race left an unlaunched controller edit in the workspace. The
working source was restored to the exact live source fingerprint, with an
explicit `source-reconciliation.json` audit; original verification is retained.
The deferred initial-phase cap is preserved in `deferred/initial-phase-cap.patch`.
The live controller switches methods after six completed CFD attempts, but
asynchronous initial proposals may exceed six. It still enforces the total
evaluation ceiling; it does not strictly reserve a fixed number of adaptive
slots. Apply the deferred refinement only to a new campaign or an audited restart.

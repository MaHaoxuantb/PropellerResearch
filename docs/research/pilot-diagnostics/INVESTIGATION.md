# Pilot investigation — 2026-10-23

## Finding

The run exhausted its numerical-failure allowance during initial sampling. It
never made a Bayesian acquisition proposal and did not establish an optimized
design. The stop was intentional, after approximately 5 h 45 min, rather than
the 12-hour ceiling. There were eight native CFD attempts: one accepted, five
not converged, and two timed out. A ninth job failed solid geometry validation.
436 additional pool proposals were rejected by the cheap geometry screen.

The attached PDF faithfully reflects the archived table, but its generic
“Native OpenFOAM + Bayesian optimization” label does not disclose that zero
Bayesian proposals actually ran. “8 evaluations completed” also includes two
terminated solves. Both descriptions should be made explicit in future reports.

## Confirmed causes

1. **Conflicting runtime budgets.** Campaign `candidate_hours=3` sets an outer
   deadline, but each native command also inherits `cfd.timeout_seconds=5400`.
   The runner takes whichever limit arrives first. Designs 005 and 008 were
   terminated at 5,400 solver seconds, at iterations 4,253 and 4,283 of 4,500.
   This is a configuration/orchestration defect, not a three-hour CPU limit.
   Both partial histories also fail residual and force checks; finishing the
   missing iterations would not guarantee acceptance.
2. **The initial search was too aggressive.** The sampler maximizes distance
   from accepted and pending points inside a normalized eight-dimensional ball
   of radius 0.5. All nine selected points have baseline distances 0.490–0.500.
   This is exploration along the edge of that ball, not small perturbations.
   Design 003 reduces b0 and b2 by about 46%, changes alpha0 from 0 to 14 degrees,
   and alpha1 from -15 to -7.66 degrees simultaneously. The configuration uses
   these as real geometry changes: alpha(u)=alpha0+alpha1*(2u-1).
3. **No learning from failed regions.** Failed jobs are excluded from exact
   duplicate proposals but not included in the maximin distance calculation,
   which uses accepted observations and pending points. There is no failure
   feasibility model or trust-region contraction. The sampler kept selecting
   similarly distant changes after failures.
4. **The Bayesian bootstrap never completed.** Acquisition requires six CFD
   attempts AND at least four accepted observations. Only the baseline and
   design 003 were available. Even after attempt six, `propose` returned
   `space-filling-local`. The previous asynchronous phase-cap issue was not
   the primary blocker: there were insufficient accepted observations anyway.
5. **The solver policy has no campaign-level continuation path.** Each candidate
   starts cold, runs to a fixed 4,500 iterations, then passes or fails a strict
   final-500-sample test. Standalone continuation support exists, but this
   campaign does not use it. It does not distinguish borderline residuals,
   ongoing drift and time-budget termination when counting failure streaks.

## Per-design numerical evidence

Thresholds: maximum initial residual over the final 500 samples <= 1e-4;
relative force/power standard deviation and half-window drift <= 2%.
All eight CFD jobs passed geometry and native mesh checks. All mesh cell counts
were between 291,501 and 358,098. No accepted status below has been changed.

| Design | Outcome | Evidence |
|---|---|---|
| 001 | Not converged | Thrust standard deviation 2.22%, narrowly above 2%. Reverse thrust about -0.148 N; power about 52.5 W. Residuals pass. |
| 002 | Not converged | Stable forces (thrust std 0.26%); peak residual 1.81e-4. Residual threshold alone rejects it. |
| 003 | Accepted | 0.232674 N, 26.708742 W, 0.008712 N/W. Worse in both objectives than baseline. |
| 004 | Not converged | Near-zero thrust (~0.00783 N), 12.29% within-window drift (~0.000962 N). Relative normalization amplifies this near-zero case, but successive 500-sample means are also drifting: 0.02075, 0.01208, 0.00783 N. |
| 005 | Timeout | Stopped at iteration 4,253; thrust std 4.67%, peak residual 2.72e-4. |
| 006 | Not converged | Thrust drift 2.91%, peak residual 3.41e-4. Last three 500-sample means rise 0.20088 -> 0.21597 -> 0.23455 N. Still settling. |
| 007 | Geometry invalid | Triangle self-intersections. No native mesh/flow solver launched. |
| 008 | Timeout | Stopped at iteration 4,283; thrust std 10.27%, peak residual 3.16e-4. |
| 009 | Not converged | Stable forces (thrust std 0.34%); peak residual 7.65e-4. Residual threshold alone rejects it. |

Numbers for failed/time-limited runs are diagnostic estimates, not optimization
observations or qualified performance. The accepted candidate has 34.1% of
baseline thrust, 2.853 times its power, and 11.9% of its thrust per watt.
The precise aerodynamic cause of its loss has not been isolated: simultaneous
changes prevent attribution to twist, width, sweep, or their interactions.
Stall, physical unsteadiness, or hardware insufficiency are not established by
these logs alone.

The fourth consecutive failure after design 003 triggered draining. Another
already-running candidate then failed, explaining a final streak of five.
Concurrency stayed at two because expansion to four required two new accepted
results. Thus the campaign did not reach its configured four-worker maximum.

## Recommended repair order

1. Unify stage/solver and whole-candidate deadlines and report the actual stop
   cause. Add a regression check that a three-hour candidate allowance cannot
   silently inherit an unintended 90-minute solver ceiling.
2. Bootstrap with small, parameter-scaled changes around the converged baseline,
   initially varying one parameter or a small group at a time. Contract the
   proposal region after repeated numerical failures. Keep fixed-radius solving
   and both existing geometry-validation stages.
3. Diagnose convergence with a bounded continuation of a borderline case
   (002) and a still-drifting case (006), preserving thresholds and original
   evidence in new attempts. This distinguishes insufficient iterations from
   a residual plateau. Do not rerun the entire pool or grant every poor design
   more time. Neither case currently demonstrates an improvement over baseline.
4. Establish enough reliable local observations before resuming Bayesian
   acquisition. Reserve adaptive slots, and distinguish geometry infeasibility,
   poor but converged performance, numerical nonconvergence and time exhaustion.
   Failed solves must not be fabricated into thrust/power observations.
5. Report the active proposal method and counts of accepted, nonconverged,
   timed-out and geometry-rejected jobs separately. Continue to label mesh/domain
   and experimental qualification as pending.

Do not simply relax thresholds to fill the plot. OpenFOAM residuals are
normalized equation measures whose normalization depends on the solution;
force stability and residual history provide different evidence. See the
[OpenFOAM residual documentation](https://doc.openfoam.com/2306/tools/processing/numerics/solvers/residuals/)
and [solution/relaxation controls](https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control).
These references support the interpretation of diagnostics, not a claim that a
specific alternate threshold or solver setting is already validated here.

## Audit trail

- Archived run: `runs/mobo-pilot-02/state.json`, per-job acknowledgements, native
  logs, `result.json`, `manifest.json`, mesh records and force/moment histories.
- All nine saved candidate parameter dictionaries match their native attempt
  manifests. No job-to-result association error was found.
- `evidence.json` beside this note contains read-only parser results, exact paths,
  source result hashes, distances and mesh diagnostics for independent review.
- Relevant implementation: `optimization/bayesian.py` (`candidate_pool`,
  `propose`), `optimization/campaign.py` (`ingest`, `worker_limit`, `run`),
  `simulation/runner.py` (`run_command`), `simulation/parser.py` (`summarize`).
- No native solves were launched, archived results changed, or convergence
  criteria relaxed during this investigation. Only diagnostic files were added.

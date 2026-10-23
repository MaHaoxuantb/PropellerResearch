# One-week shape optimization recommendation

Research date: 2026-10-22. Status: recommendation and implementation plan; the Bayesian optimizer and campaign have not been launched.

## Decision and objective

Use parallel multi-objective Bayesian optimization (MOBO), with BoTorch's qLogNoisyExpectedHypervolumeImprovement acquisition function. The user's clarified deliverable is a curve of thrust per watt versus thrust at a fixed propeller radius, with continuous use of the Mac's available compute for no more than seven days.

Choose an actual radius of 111.76 mm to retain the existing baseline. Keep 5,000 RPM, static air, two identical loop units, hub, fluid properties and profile mapping fixed in the first campaign. These are the existing study conditions, not a universal operating optimum.

Maximize the two outputs simultaneously:

- thrust T in N;
- power loading T/P in N/W, with positive shaft power P.

Display T on the horizontal axis and T/P on the vertical axis. Retain nondominated designs: no other accepted design has both at least as much thrust and at least as much thrust per watt, with one strict improvement. This is an estimated Pareto frontier across different shapes at the same operating point, not an RPM sweep of one propeller. It may be a sparse set of points rather than a smooth curve.

No numerical minimum thrust was supplied. Mark the baseline thrust as a practical reference, report the best ratio at or above that reference, a compromise region, and the highest-thrust confirmed design. Do not invent a unique optimum or silently impose a thrust requirement. A power/motor limit can later be added if one is specified. Static propulsive efficiency T*V/P is zero at V=0 and is not the search objective.

## Variables

The generator retains H, beta, shape_exponent, b0, b1, b2, Z, alpha0 and alpha1. Search eight independent normalized variables and derive H to satisfy the fixed finished-solid radius, initially to within 0.05 mm. Use a bounded root solve and reject candidates with no feasible H or ambiguous discontinuous geometry. Verify the final Boolean solid radius; do not assume the reference-curve radius equals the rotor radius and do not rescale the fixed hub.

Use the existing bounds initially, with cheap envelope/root screening followed by full geometry validation before CFD. Screen a large candidate pool cheaply and select a space-filling valid initial set. Filtering destroys exact Sobol balance, so do not claim an unfiltered low-discrepancy guarantee for the surviving designs. Keep blade count, airfoil family and additional spline variables outside this first campaign.

## Why this method

MOBO learns a probabilistic approximation from previous expensive evaluations and chooses new shapes expected to expand the trade-off frontier. NEHVI was developed for noisy, parallel multi-objective searches [1]; BoTorch recommends the log implementation for numerical robustness [2]. This is a reasoned choice for this problem, not a measured optimizer comparison on this rotor.

SA is implemented and useful as a fallback, but its current scalar power-plus-thrust-penalty objective does not produce the requested frontier. It would require multiple scalarizations or a multi-objective adaptation. Its local acceptance mechanism also does not use a learned response model to allocate expensive CFD calls.

CMA-ES is a strong derivative-free method for continuous shape problems [3]. A population-based search or multi-objective evolutionary method remains viable, but I would prioritize evaluation efficiency with a surrogate for this budget. Ordinary finite-difference optimization is less attractive with meshing failures and remeshing discontinuities; adjoint infrastructure does not currently exist here.

Use ordinary GP-based MOBO first at eight independent variables and roughly a hundred observations. MORBO's multiple local models target much larger dimensions [4]; add that complexity only if global surrogate diagnostics or search stagnation justify it. SCBO is suitable for the earlier single-objective constrained formulation, not the selected two-objective frontier [5].

## CFD evidence and prerequisites

Read-only reanalysis of the refined continuation's completed native log found clean completion at iteration 4,500, all existing acceptance checks passed, T=0.683205696 N and P=9.361453662 W. Thus T/P=0.072980727 N/W (P/T=13.702248 W/N). Maximum initial residual across the final 500 samples was about 7.03e-6, below 1e-4. Thrust drift was about 0.0754%, below 2%.

Source: runs/evaluations/6388c6b553280fec/ec0cf993e6f34d22b482e80e207476a5. Its parent supplied the first 2,500 iterations. The combined solver time is about 80 minutes, or about 82 minutes including original setup/geometry/meshing. This is one measured geometry, not a guaranteed candidate runtime.

The repaired coarse mesh passed, but the coarse solution at 2,500 iterations remained unconverged. Do not interpret its difference from the refined solution as an accepted mesh error. Neither mesh/domain independence nor physical validation has been established for this profile family.

Both launch supervisors exited before publishing final result.json files while native child solvers continued. Before autonomous optimization, implement supervised durable jobs and audited recovery of completed native outputs. Preserve original logs, distinguish recovered completion from a witnessed process exit code, and never mark an active native child interrupted just because its Python owner is gone.

Bound the remaining preflight: recover results, verify the fixed-radius parameter map, test a few perturbations, and check baseline sensitivity. Use one consistent refined search fidelity; reserve finer/domain checks for baseline and frontier finalists. If baseline sensitivity is unresolved after the preflight budget, label the search exploratory at the selected fidelity and retain that limitation. If even stable numerical evaluation fails, stop and report it rather than spend the week producing unrankable cases.

## Search loop

1. Recover the accepted numerical baseline; establish a compatible frozen search case and schema.
2. Generate approximately 24–32 diverse valid initial designs including baseline-neighbour designs. Start with 32 as the planning target, reducing count if measured cost demands it.
3. Fit normalized GP models to accepted CFD outputs. One useful implementation models T and log(P), then computes T/P from posterior samples; this keeps the shared thrust term consistent across objectives and prevents negative predicted power in the ratio. Require positive thrust for the propulsion frontier. Treat numerical failures separately from aerodynamic infeasibility; never fabricate zero power or a large objective as if it were an observation.
4. Choose batches of 2–4 using qLogNEHVI, accounting for pending jobs. Refit as observations complete. Freeze the baseline normalization and hypervolume reference after initialization. The reference point is an acquisition setting, not a hidden thrust requirement.
5. Apply full geometry checks, native meshing and sustained convergence checks to every proposed CFD evaluation. Continue only promising settling histories within a bounded per-candidate allowance. Keep numerically rejected samples out of the accepted frontier.
6. Re-evaluate 3–5 representative frontier designs and the baseline at matched confirmation settings, with independent initialization where useful. Compare ratios at adequate thrust and verify that apparent gains survive resolution/domain changes.
7. Deliver accepted points, uncertainty/sensitivity evidence, geometry previews, parameters and exports. Show surrogate predictions separately; do not draw an interpolating curve that implies untested designs are confirmed.

GP uncertainty is not a CFD discretization-error estimate. Steady iteration samples are correlated settling histories, not independent experimental repetitions. Their standard deviation must not be divided by sqrt(N) and presented as physical confidence.

## Mac mini budget and execution

Measured hardware: 16 GiB RAM, eight physical/logical CPU cores. Free disk at inspection: about 25 GiB; existing evaluation archive about 2.1 GiB.

Benchmark two versus four simultaneous native serial cases and choose the highest measured throughput without swapping. Four workers is a candidate setting, not an assumed fourfold speedup. Current refined peak child RSS was about 1.47 GB; plan conservatively around 3 GB per search worker and a total campaign memory ceiling near 12 GB, reducing workers for larger meshes. The native MPI path previously failed; parallel independent designs avoids making MPI repair a prerequisite. The M2 GPU is not assumed to accelerate this OpenFOAM installation.

Working estimate: 80–120 distinct search designs if measured parallel throughput supports it, plus a small set of confirmation evaluations. A slow/failure-heavy campaign may produce substantially fewer. One candidate taking 1.5 hours with two effective workers yields about 1.33 candidates/hour; 100 search hours then permits approximately 133 attempts before overhead and failures. Use live measured costs and remaining time, not this extrapolation, to schedule work.

| Elapsed allocation | Purpose |
| --- | --- |
| First 12 hours | Reliable supervision/recovery, fixed-radius adapter, baseline checks and throughput benchmark |
| Next 36 hours | Initial diverse CFD designs; begin fitting models as results arrive |
| Next 72 hours | Adaptive MOBO; reduce or expand candidate count to fit actual throughput |
| Next 36 hours | Confirm frontier representatives and baseline |
| Final 12 hours | Contingency, plots, exports and report |

Total: 168 hours including preparation. This is a ceiling, not an instruction to consume all seven days. If preparation takes longer, reduce search time; never extend the deadline. Reserve confirmation time before scheduling more exploration. Limit acquisition/model-fitting time as well as solver time.

Implement an absolute campaign deadline that survives restarts, checks remaining time before launches and continuations, and terminates owned process groups when their bounded allocation expires. The current SA elapsed-time check runs only between evaluations and is insufficient for this hard cap. Existing per-command timeouts also do not enforce a whole-campaign deadline.

Monitor disk before each launch with an initial campaign growth cap around 12 GiB and a free-space reserve. Preserve existing archives. Prefer fewer evaluations over silently deleting evidence; any compact archive design must retain provenance and permit explicitly documented rehydration. Stop early when insufficient time/disk remains for a complete evaluation and confirmation.

## Implementation order

1. Durable worker supervisor, recovery and absolute wall-clock/resource accounting.
2. Fixed-radius adapter and revised output/objective schema; retain SA as an optional backend.
3. BoTorch MOBO backend, initial valid design selection, pending-job awareness and reproducible checkpoints. Check and pin a compatible Python/PyTorch/BoTorch environment separately from the existing CFD environment before use.
4. Campaign manifest, live Pareto report and finalist confirmation queue.
5. Short end-to-end smoke campaign, then the bounded week campaign.

This research turn does not install dependencies or start the week-long search.

## Sources

[1] Daulton, Balandat & Bakshy, Parallel Bayesian Optimization of Multiple Noisy Objectives with Expected Hypervolume Improvement (2021): https://arxiv.org/abs/2105.08195

[2] BoTorch multi-objective tutorial and recommendation to replace qNEHVI with qLogNEHVI: https://botorch.org/docs/v0.17.2/tutorials/multi_objective_bo

[3] Hansen's CMA-ES overview: https://cma-es.github.io/

[4] Daulton et al., Multi-objective Bayesian optimization over high-dimensional search spaces (MORBO, 2022): https://proceedings.mlr.press/v180/daulton22a.html

[5] Eriksson & Poloczek, Scalable Constrained Bayesian Optimization (2021): https://proceedings.mlr.press/v130/eriksson21a.html

[6] BoTorch multi-fidelity KG tutorial: https://botorch.org/docs/v0.15.1/tutorials/multi_fidelity_bo . Multi-fidelity modelling is a future option only after paired, converged observations establish useful cross-fidelity correlation. It is not a justification for mixing unconverged coarse outputs with accepted refined values.

# Filter Utility Extraction Report

## Purpose and scope

This report identifies code in `examples/tetherFilterDemo.m` and
`examples/tetherCovarianceDemo.m` that is likely to be useful in future
simulation, estimation, and analysis scripts. It also proposes an end-to-end
flow for a refactored filter demonstration that can run independent trials in
parallel.

The recommendations are intentionally broad, but they distinguish reusable
filter and analysis utilities from scenario-specific setup and presentation
code. This is an architectural recommendation, not a claim that the utilities
have already been extracted.

## Findings from the current demos

The filter demo already contains a useful separation of concerns inside one
file: initial conditions and truth propagation, filter propagation and
measurement updates, Monte Carlo collection, statistical analyses, and
visualization. The most portable calculations are currently local functions,
so other scripts cannot call them directly.

There are two important behavioral constraints to address before using the
demo as a parallel Monte Carlo driver:

1. `runMonteCarlo` is a **state-only comparison run**, not an independent
   full-filter trial. It reuses the reference run's covariance and Kalman
   gains, while generating new truth and measurements. That is useful for
   comparing state-error samples under fixed gains, but it does not produce
   trial-specific covariance histories for a full-filter NEES comparison.
2. The demo seeds MATLAB's global random stream and consumes random numbers
   during truth perturbation and measurement generation. Parallel trials need
   independent, reproducible streams assigned by trial index rather than a
   shared global stream.

The covariance demo provides further evidence for what is reusable and what
is domain-specific. It has its own tether initial-state and single-body
propagation helpers, while the filter demo has overlapping implementations.
Those are candidates for shared simulation support if both demos should use
identical behavior. The covariance demo also includes useful plot-specific
calculations, such as extracting a covariance block's standard deviation and
decomposing relative velocity uncertainty along and across the tether.

## Recommended utility inventory

### Filter execution and numerical safeguards

| Candidate utility | Existing implementation | Recommendation |
|---|---|---|
| Continuous state-and-covariance propagation between measurement times | `runFilter` in `tetherFilterDemo.m`; dynamics wrapper in `dynamicsFunction` | Extract a reusable propagation/update driver that accepts dynamics, initial state/covariance, measurement times, observations, and model callbacks. Keep model-specific parameters outside the driver. |
| EKF measurement update | `measurementUpdate` in `tetherFilterDemo.m` | Separate the generic update from synthetic observation generation. A reusable update should accept an observation, predicted measurement, Jacobian, and noise covariance; truth-based noisy measurement generation belongs in simulation support. Preserve the Joseph covariance update. |
| Packed covariance conversion and symmetry cleanup | `unpackCovariance` in `tetherFilterDemo.m` | Extract, with explicit state dimension or a documented fixed 18-state contract. Consider a matching pack helper to keep the vectorization convention explicit. |
| Covariance health validation | Currently spread across covariance use and propagation | Add reusable checks for finite values, symmetry error, and positive definiteness/Cholesky success. Return explicit health status and diagnostic values rather than silently repairing an invalid covariance. |
| Angle innovation wrapping | `wrapAngleInnovation` in `tetherFilterDemo.m` | Extract as a measurement-model-aware helper. The current row-count logic assumes the implemented measurement-group order; make the angular row indices explicit or derive them from shared measurement metadata. |
| Unit scaling and state summaries | `stateNorm` in `tetherFilterDemo.m` | Extract conversion/scaling utilities for position and velocity. Prefer explicit scale vectors or named conversion functions over a generic string selector. |

### Simulation and initial-condition support

| Candidate utility | Existing implementation | Recommendation |
|---|---|---|
| Tethered-body initial state | `tetherInitialState` in both demos | Consolidate if both demos are expected to share the same initialization convention. Document the input frame, units, mass-fraction convention, and output ordering. |
| Single-body propagation | `propagateSingleBody` in both demos | Consolidate as a simulation helper with solver tolerances supplied explicitly or through a small options structure. |
| Scenario truth propagation | `propagateScenario` and `perturbInitialState` in `tetherFilterDemo.m` | Keep as simulation utilities, separate from the production filter API. Pass a per-trial random stream to perturbation and observation generation. |
| Tether/rendezvous scenario construction | `makeInitialEstimate` and `rendezvousInitialState` in `tetherFilterDemo.m`, plus equivalent setup in the covariance demo | Keep scenario-specific construction in an example or simulation-support layer. It depends on this tethered CR3BP scenario rather than on EKF mechanics. |
| Dynamics callback construction | `dynamicsFunction` in `tetherFilterDemo.m` | Extract only if more than one caller needs the same callback assembly. Keep the underlying dynamics in `toolbox/+coa/+orbit/`. |

### Diagnostics and statistical analysis

| Candidate utility | Existing implementation | Recommendation |
|---|---|---|
| NEES per trial and state block | `computeNees` and `neesBlocks` in `tetherFilterDemo.m` | Extract as truth-based consistency metrics. Accept arbitrary block definitions so orbit, body, position/velocity, or other subsystem partitions can be analyzed. |
| ANEES chi-square bounds | `neesBounds` in `tetherFilterDemo.m` | Extract alongside NEES; parameterize run count, degrees of freedom, and confidence level. |
| NIS and innovation diagnostics | Innovations and innovation covariances are stored by `measurementUpdate`, but NIS is not currently computed | Add per-update NIS using a linear solve rather than an explicit matrix inverse. Include degrees of freedom and a configurable statistical test/window in the result. Add innovation mean/whiteness diagnostics as offline or windowed analyses. |
| Error and covariance histories | `emptyHistory` and `appendPredictionHistory` in `tetherFilterDemo.m` | Extract an optional history collector. Make sample retention configurable so large parallel batches can return compact summaries when full histories are unnecessary. |
| Relative-state covariance mapping | `relativeStateDistribution` and `tetherFrame` in `tetherFilterDemo.m` | Keep the generic linear map \(y=Mx\), \(P_y=MP_xM^\mathsf{T}\) reusable; keep the tether-frame construction in domain-specific analysis support. Validate degenerate tether geometry before normalizing vectors. |
| Covariance block summaries | `diagBlocks` in `tetherFilterDemo.m`; `blockSigma` in `tetherCovarianceDemo.m` | Consolidate as small analysis helpers for diagonal and block covariance summaries. Define whether each output is componentwise sigma, trace-based sigma, or another statistic. |
| Measurement group metadata | `activeMeasurementGroups` in `tetherFilterDemo.m` | Centralize group names, row counts, and row ordering with the measurement model or its documentation. Reuse it for angle wrapping, diagnostic attribution, and plots to avoid duplicated row assumptions. |

### Reporting and visualization

The following routines are helpful but should not be required by a worker
running an individual trial:

- ANEES report and plot: `reportAnees`, `plotAnees`.
- RMSE versus predicted standard deviation: `plotRmse`.
- Innovation and post-fit residual plots: `plotInnovations`, `plotResiduals`,
  and `plotMeasurementHistory`.
- State-error/covariance corner plots: `plotCorner` and `drawCornerPlot`.
- Relative-state reporting and plots: `reportRelativeState` and
  `plotRelativeCorner`.
- Signed component error and uncertainty plots:
  `plotPredictionHistory`, `symmetricPlotLimit`, and `plotSignedComponent`.
- Human-readable state labels and body mapping: `stateLabels` and `bodyOf`.
- Analysis registration and dispatch: `monteCarloAnalyses` and
  `runMonteCarloAnalyses`.

These are appropriate for an offline analysis layer. Workers should return
data, not create figures, write shared output files, or interleave console
reports.

## Suggested organization

Follow the repository's existing `toolbox/+coa/` package convention. One
possible organization is:

```text
toolbox/+coa/+filter/
    ekfUpdate.m
    propagate.m
    wrapAngleInnovation.m
    covarianceHealth.m
    nis.m
    nees.m
    consistencyBounds.m
    mapGaussian.m

examples/support/
    tetherInitialState.m
    propagateSingleBody.m
    makeTetherScenario.m
    runFilterTrial.m
    runMonteCarloTrials.m
    collectFilterHistory.m
    analyzeFilterRuns.m
    plotFilterDiagnostics.m
```

This is a suggested boundary, not a requirement to create every file. Put
generally reusable estimation mathematics under `toolbox/+coa/+filter/`.
Scenario generation, trial orchestration, demo-specific block definitions,
and graphics can remain under examples/support. Avoid moving tether physics
or measurement models out of their existing domain packages.

## Proposed flow of the final script

The top-level example should become a short orchestration script. Its flow
would be:

1. **Load configuration and validate it.** Define physical parameters,
   solver tolerances, measurement schedule, noise models, Monte Carlo count,
   diagnostic thresholds, and which outputs to retain. Validate dimensions,
   finite values, and time ordering before starting workers.
2. **Construct deterministic scenario inputs.** Create the nominal initial
   estimate, initial covariance, dynamics/model callbacks, and any shared
   scenario data once. These shared inputs should be read-only during trials.
3. **Assign trial identities and random streams.** Derive one reproducible
   stream/substream per trial index. The same trial index should reproduce
   the same truth perturbations and measurement noise independent of worker
   scheduling.
4. **Run independent full-filter trials in parallel.** Each worker receives
   the trial index, configuration, nominal inputs, and its own random stream.
   It independently generates truth and observations, propagates its own
   state and covariance, and performs its own measurement updates. It returns
   a result struct without graphics or shared-file writes.
5. **Compute per-update online diagnostics within each trial.** At minimum,
   return NIS and its degrees of freedom, plus covariance-health results.
   Optionally maintain a sequential/windowed consistency monitor. Record
   alarm state, trigger time/update, and the metric values that caused it;
   distinguish a numerical failure from a statistical inconsistency.
6. **Collect trial results in deterministic order.** Aggregate by trial index,
   not completion order. Retain the full trajectory only when needed; otherwise
   return update-level metrics and selected endpoint summaries to limit memory.
7. **Compute ensemble diagnostics after collection.** When truth and the
   trial-specific covariance are available, compute NEES/ANEES and confidence
   bounds using each trial's own covariance. Also compute RMSE, innovation
   mean/whiteness, alarm rates, false-alarm rates on nominal cases, and
   detection delay on injected-fault cases as appropriate.
8. **Report and visualize centrally.** Produce concise tables and plots after
   all workers finish. Keep plotting and file output out of the parallel loop.
9. **Save versioned results and configuration.** Persist enough information
   to reproduce a run: configuration, seed/stream policy, trial identifiers,
   software revision, diagnostic settings, and summarized outputs.

In pseudocode:

```matlab
config = makeConfig();
validateConfig(config);
inputs = makeScenarioInputs(config);
trial_ids = 1:config.trial_count;

parfor trial_id = trial_ids
    stream = makeTrialStream(config.seed, trial_id);
    trials(trial_id) = runFilterTrial(inputs, config, stream);
end

summary = analyzeFilterRuns(trials, config);
reportFilterRuns(summary, config);
plotFilterDiagnostics(summary, config);
saveFilterResults(trials, summary, config);
```

The real `runFilterTrial` should be an independent filter run, not the current
state-only Monte Carlo comparison. A separate state-only/fixed-gain experiment
can remain available under an explicitly descriptive name if it is still
useful.

## Practical extraction order

1. Extract and test the generic EKF update, covariance conversion, angle
   wrapping, and state/covariance health checks.
2. Extract NIS and the truth-based NEES/ANEES calculations and bounds. Keep
   statistical tests separate from their display.
3. Refactor the filter demo so synthetic observation generation and filter
   updates are separate, then implement one independently seeded full-filter
   trial.
4. Add deterministic per-trial random streams and parallel orchestration.
   First verify that serial and parallel execution produce the same result
   for each trial ID.
5. Extract repeated tether initialization and free-body propagation only if
   both demos need the same shared implementation.
6. Move reusable history aggregation and relative-state mapping out of the
   demo; keep specialized plot layouts and labels in analysis support.

## Validation criteria

Before relying on the parallel version, validate both numerical behavior and
the automation contract:

- A fixed trial ID and seed reproduce the same result in serial and parallel
  execution, regardless of worker scheduling.
- Every trial runs its own covariance propagation and gain/update sequence.
- NIS agrees with a direct calculation on small test cases and uses the
  correct measurement degrees of freedom.
- NEES uses the covariance from the same trial and update as its state error.
- Covariance checks detect nonfinite, asymmetric, and non-positive-definite
  cases without silently converting them into healthy results.
- Simulated nominal runs establish the false-alarm rate; injected biases,
  outliers, and dynamics mismatch measure detection rate and delay.
- Results can be aggregated and reported without relying on worker execution
  order or worker-side graphics/file side effects.

The resulting boundary keeps the numerical filter reusable, makes each
parallel trial self-contained and reproducible, and reserves simulation and
visualization details for the example/analysis layer.

# Changelog

## [Unreleased]
- `coa.orbit.tetherAccel`, `tetherSpringDamperAccel`, and
  `tetheredThreeBodyEkfPropagation` now use acceleration callbacks with
  signature `accel_fn(time, state)`. `tetherAccel` now takes time as its
  first argument; covariance propagation still requires state Jacobians.
- `derivations/ekf_equations.tex`: the sensor table now follows
  `coa.meas.measurementModels`. Sensors are grouped and ordered by mask bit,
  with each group's bit, mask value and row count, and each sensor's stacked
  rows and nominal noise from `coa.meas.measurementCovariances`. Crosslink
  rows are renumbered to the code's order (13-14 range, 15-16 Doppler) and use
  its vectors, `r2 - r3` and `r2 - r1`. The unimplemented accelerometer
  candidates (old rows 20-22) are removed. Code comments and README updated to
  match.
- Initial repository structure, build file, and local Git hooks.
- Added `coa.filter.ekfUpdate`, a reusable observation-driven Joseph-form EKF
  measurement update, extracted from `examples/tetherFilterDemo.m` (which now
  generates its synthetic observations separately and calls it).
- `coa.filter.ekfUpdate` accepts an optional `max_iterations` for an iterated
  EKF update (default 1, the ordinary update); `tetherFilterDemo` uses 5.
- `examples/tetherFilterDemo.m`: each plot can now be set to `"off"`, `"save"`
  (export to file and close) or `"show"` via `config.plot_modes`; saved files go
  to `config.plot_output_folder`.
- Added `coa.orbit.tetherSpringDamperAccel`, a states-only spring-damper tether
  acceleration that drives the m1-m2 separation to `L(t) = L_0 + Ldot * t`, and
  `examples/tetherSpringDamperAccelDemo.m`, which compares a critically damped
  tether (`c = 2*sqrt(k)`) against one with no spring or damper over a one-day
  CR3BP integration starting from a broken length constraint, for several Ldot.
  Figures cover length and rate error against the exact theory, tension and its
  feedback split, trajectory/spin divergence, geometry, and solver cost.
  `tetherSpringDamperAccel` also returns Jacobians (`[a, da_dp, da_dv, tension]`),
  and `coa.orbit.tetheredThreeBodyEkfPropagation` accepts an optional
  `spring_damper` struct (`initial_length`, `stiffness`, `damping`) to use it in
  place of the rigid `tetherAccel`. `examples/tetherHaloDemo.m` now uses the
  critically damped spring-damper by default (`use_spring_damper = false`
  restores the rigid constraint). `examples/tetherCovarianceDemo.m` now
  uses the spring-damper tether with a COM/relative initial covariance
  (`coa.filter.tetherInitialCovariance`), checks it against the rigid tether,
  and compares the linear covariance (including spin phase) with a Monte Carlo
  run, for zero and nonzero (reel in/out) length rates.

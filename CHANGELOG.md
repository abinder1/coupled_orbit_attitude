# Changelog

## [Unreleased]
- Initial repository structure, build file, and CI workflow.
- Added `coa.filter.ekfUpdate`, a reusable observation-driven Joseph-form EKF
  measurement update, extracted from `examples/tetherFilterDemo.m` (which now
  generates its synthetic observations separately and calls it).
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

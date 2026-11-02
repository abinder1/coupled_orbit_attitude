# Changelog

## [Unreleased]
- Initial repository structure, build file, and CI workflow.
- Added `coa.filter.ekfUpdate`, a reusable observation-driven Joseph-form EKF
  measurement update, extracted from `examples/tetherFilterDemo.m` (which now
  generates its synthetic observations separately and calls it).
- `examples/tetherFilterDemo.m`: each plot can now be set to `"off"`, `"save"`
  (export to file and close) or `"show"` via `config.plot_modes`; saved files go
  to `config.plot_output_folder`.
- Added `coa.orbit.tetherSpringDamperDynamics` (state-only spring-damper tether
  with optional centripetal/gravity-gradient feedforward) and the comparison and
  stiffness/damping sweep `examples/tetherSpringDamperDemo.m`.

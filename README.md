# coupled_orbit_attitude

Simulation code and derivations for coupled orbit-attitude dynamics.

## Layout

| Folder | Contents |
|---|---|
| `toolbox/+coa/` | Production code, in the `coa` namespace (`coa.orbit.*`, `coa.meas.*`, `coa.filter.*`, `coa.util.*`) |
| `tests/unit/` | Fast tests of derivation results (limits, symmetries, identities) |
| `tests/regression/` | Short fixed-seed runs compared against references in `tests/data/` |
| `derivations/` | LaTeX derivations; equation labels (`eq:...`) are referenced in code comments |
| `examples/` | Short runnable demos |
| `scratch/` | Exploratory work; never called from `toolbox/` |
| `__archive/` | Legacy scripts kept for reference |

## Usage

Requires MATLAB R2024a or later. From the repository root:

```matlab
buildtool          % add toolbox to path, run Code Analyzer, run all tests
buildtool clean    % delete the XML test and coverage results
```

Test results and coverage are written to `test_results/`
(`test-results.xml` and `coverage.xml`).

## Test and coverage report

An interactive HTML report lets you browse test results and line coverage
instead of reading XML. It needs only Python 3 (standard library) on `PATH`; no
packages to install. `--run` also needs `matlab` on `PATH`.

### Quick start (recommended: live mode)

1. From the repository root, start the report in live mode:

   ```powershell
   python test_results/test_report.py --watch --open
   ```

   This opens `test_results/report.html` in your default browser and keeps
   running in the terminal. Leave it running (Ctrl+C stops it).

2. In a second terminal (or MATLAB), run the tests:

   ```powershell
   matlab -batch "buildtool clean test"
   ```

3. Watch the browser. When the run finishes the page updates itself within
   a few seconds; no refresh needed. The header shows when the results were
   written, followed by "live: updates automatically".

In VS Code you can do the same from Terminal > Run Task: **Run tests and open
report** (starts live mode and runs the tests), or **Test report (live)** and
**Run tests** separately. The tasks are defined in `.vscode/tasks.json`.

Notes:

- Plain `buildtool test` skips the run as "up to date" when nothing changed, so
  no new results are written and the page correctly stays the same. Use
  `buildtool clean test` to force a fresh run.
- In live mode an update keeps your current tab and open source file; text typed
  into the test filter box is cleared.

### Other modes

```powershell
python test_results/test_report.py              # write report.html once
python test_results/test_report.py --open       # ...and open it
python test_results/test_report.py --run --open # run buildtool test, then open
```

Without `--watch`, the page is a snapshot of the results at the time you ran the
script. Re-run the script and refresh the browser to see newer results. The
script works from any working directory. Options: `--results DIR` (folder
containing the XML files) and `--output FILE` (HTML path).

### Reading the report

- **Unit tests tab:** summary cards (tests, passed, failed/errors, skipped,
  suites, total time); one expandable section per test suite; a text filter and
  a status filter (all, passed, failed, skipped); failing suites open
  automatically with their failure messages; a "Slowest tests" table.
- **Code coverage tab:** overall line coverage and counts of covered and missed
  lines; a per-file table with coverage bars (green at 90% or more, orange at
  75% or more, red below), sortable by clicking a column header. Click a file to
  see its source with per-line hit counts: green lines were executed, red lines
  were never executed. Use "< All files" to go back.

Coverage is line-based only; MATLAB's Cobertura output has no branch data, and
only code under `toolbox/` is measured.

### Files

`test_results/test_report.py` is the only tracked file in `test_results/`.
`report.html`, `report_data.js` and the XML results are regenerated and
git-ignored. `buildtool clean` deletes only the two XML result files; the
report files remain until the next run of the script overwrites them.

## Tether initial state

`coa.orbit.tetherInitialState(com_state, length_nd, mass_2_fraction, ...
spin_rate_nd, length_dot, direction)` constructs the rotating-frame state
`[p1; p2; v1; v2]` from a center-of-mass state. The direction is a unit vector
from the center of mass toward body 1 and is used without normalization. The
spin is inertial about +z; the nondimensional frame rotates about +z at rate 1.
See the function help for the input conventions.

`coa.orbit.tetherFrame(relative_position, relative_velocity)` constructs a
right-handed tether-aligned orientation using relative vectors in the same
coordinate frame. Its rows are the frame axes; it does not transform velocities
or covariances or add rotating-frame transport terms.

## Measurement models

The `coa.meas` package models measurements for the 18-element Cartesian state
`[r1; r2; r3; v1; v2; v3]`, in the same frame and nondimensional units used
by the dynamics. `coa.meas.measurementModels(state, mask)` returns the selected
predicted measurements; optional second and third outputs are
the state Jacobian and measurement-noise covariance, respectively. The
Jacobian is evaluated only when requested. Noise covariances are needed and
validated only when the third output is requested.

Each bit in the integer `mask` selects one group (group 1 is bit 0):

| Group | Mask bit | Measurement rows |
|---|---:|---|
| 1 | `2^0` | Ground angles for bodies 1-3 (6) |
| 2 | `2^1` | Relative camera and star-tracker angles (4) |
| 3 | `2^2` | Tether angular size (1) |
| 4 | `2^3` | DSN ranges for bodies 1-3 (3) |
| 5 | `2^4` | DSN range rates for bodies 1-3 (3) |
| 6 | `2^5` | Crosslink ranges, `r2-r3` then `r2-r1` (2) |
| 7 | `2^6` | Crosslink range rates, `r2-r3` then `r2-r1` (2) |
| 8 | `2^7` | Lunar PN&T/GNSS position for bodies 1-3 (9) |

Selected groups are stacked in ascending group order; rows within each group
follow the order above and the corresponding function documentation. A zero
mask returns no measurement rows. The mask must be a nonnegative integer.
The third output is block diagonal across
selected rows/groups; covariance fields are required only for selected groups
and must be finite, real, symmetric, positive definite matrices of the
documented size.

`[groups, angle_rows] = coa.meas.measurementGroups(mask)` returns the selected
group names, group indices, and stacked row indices, plus the angular rows to
wrap during an EKF update. It uses the same group order as `measurementModels`.

`coa.meas.measurementCovariances(length_scale, speed_scale)` supplies nominal
measurement-noise covariances using characteristic scales in km and km/s.
Both scales must be finite and positive. These covariances describe sample
noise, not sensor bias or drift. Candidate
accelerometers in sensor-table rows 20-22 of `derivations/ekf_equations.tex`
are not implemented in the current measurement model.

## Filter update

`coa.filter.ekfUpdate(prior_state, prior_covariance, observation, measurement_model, angle_rows, max_iterations)`
applies a Joseph-form EKF measurement update to a prior estimate. The optional
`max_iterations` (default 1) relinearizes the model about the latest iterate for
an iterated EKF, with early exit on a small Mahalanobis step. It is
deterministic and takes an *observation*, not truth, so the same function serves
simulated and real data; simulated observations (and their random numbers) are
generated by the caller. `measurement_model` is a handle returning
`[h, H, R]` for a state, for example
`@(x) coa.meas.measurementModels(x, mask, covs)`. `angle_rows` lists the
angular measurement rows to wrap to `[-pi, pi)`.

Besides the posterior state and covariance, an `info` struct returns the
innovation, innovation covariance, gain, and post-fit residual and its
covariance, which are the inputs to normalized-innovation-squared (NIS) health
checks. See `help coa.filter.ekfUpdate` for the full contract, and
`derivations/ekf_health_diagnostics.tex` for the diagnostics.

`coa.filter.computeNees(errors, covariances, blocks)` computes normalized
estimation error squared for each trial, update, and caller-defined state
block. Errors have shape `[run, update, state]`; covariances have shape
`[state, state, update, run]` and use each trial's own covariance. Each block
is a struct with a `rows` vector of state indices. Selected covariance blocks
must be nonsingular.

`coa.filter.meanChiSquareBounds(run_count, dof, confidence)` returns two-sided
acceptance bounds for the mean of independent chi-square statistics, such as
ANEES or ANIS. Degrees of freedom may be supplied as a scalar or matrix.

Run `examples/tetherFilterDemo.m` in MATLAB for an example that propagates the
tethered three-body estimate and covariance for 24 hours, applying six noisy
measurement updates, with the last update 30 minutes before the end. The
iterated EKF uses Joseph-form covariance updates and all eight implemented
measurement groups. The script holds the configuration and high-level
workflow; its scenario, filter, Monte Carlo, analysis, and plotting helpers
are grouped in `examples/+tetherFilterDemoSupport`.

The demo displays the truth trajectories of the three bodies, alongside the
payload-to-catch range and inbound speed. The prediction history creates one
figure per body, with absolute componentwise position and velocity errors and
shaded regions below the positive 3-sigma bounds on logarithmic y-axes. Exact
zero values cannot appear on a logarithmic axis. The range plot marks
measurement times, the RMSE plot uses logarithmic y-axes, and relative-state
corner plots show velocity uncertainty in m/s.

Each plot is controlled by `config.plot_modes` in `demoConfig`
(`prediction_history`, `trajectory`, `anees`, `rmse`, `innovations`, `residuals`, `corner`,
`relative_corner`). Set a plot to `"off"` to skip it (including its data
preparation), `"save"` to export it to `config.plot_output_folder` as
`config.plot_file_format` at `config.plot_resolution` dpi and close it, or
`"show"` to display it. For batch or parallel runs, set every mode to `"off"`
or `"save"`.

## Developer setup (Windows)

Quick steps to get the repository ready for development and automated checks.

1. Install local git hooks (PowerShell):

```powershell
pwsh scripts/install_git_hooks.ps1
```

2. Install `pre-commit` for formatter + linter hooks:

```powershell
py -m pip install --user pre-commit
py -m pre_commit run --all-files
```

The MATLAB pre-commit pipeline is defined in `.pre-commit-config.yaml` and uses:
- `mh_style --fix` (auto-format/style fixes)
- `mh_lint` (lint/static analysis)

Project defaults for MISS_HIT are in `miss_hit.cfg`.

The test/coverage report (`test_results/test_report.py`) needs Python 3 on
`PATH` (already required above for `pre-commit`) and no additional packages.
Its `--run` option also needs `matlab` on `PATH`.

3. Run MATLAB build checks locally (same commands used by hook tasks):

```powershell
pwsh -Command "matlab -batch \"buildtool('check')\""
```

4. Run full tests (same as pre-push):

```powershell
pwsh -Command "matlab -batch \"buildtool('test')\""
```

5. See repository developer docs and checklist:

- `docs/HOOKS.md` — hook behavior and installer notes (Windows specifics)
- `docs/COMMIT_CHECKLIST.md` — concise commit/push checklist (Windows-focused)

CI strategy (local-only)

This repository uses local hooks (`.githooks/pre-commit` and `.githooks/pre-push`) for automated checks on your machine. This avoids hosted runner costs and MATLAB licensing issues.

If in future you want remote PR gating or cross-release testing, use a self-hosted runner with an activated MATLAB license or re-add a workflow and provide a license secret.

If you need to bypass hooks for a one-off commit or push, use `--no-verify` with `git commit` or `git push`.

## Conventions

- Every implemented equation carries a `\label{eq:...}` in `derivations/` and an
  `% IMPL:` comment pointing to its code; the code cites the label back.
  `grep -r "eq:name"` finds the equation, its implementation, and its tests.
- Symbols and their code names are defined in `derivations/nomenclature.tex`.
- Promote work from `scratch/` into `toolbox/` only with tests.

# coupled_orbit_attitude

Simulation code and derivations for coupled orbit-attitude dynamics.

## Layout

| Folder | Contents |
|---|---|
| `toolbox/+coa/` | Production code, in the `coa` namespace (`coa.orbit.*`, `coa.attitude.*`, `coa.util.*`) |
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
buildtool clean    % remove results/
```

Test results and coverage are written to `results/`.

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

`coa.meas.measurementCovariances(length_scale, speed_scale)` supplies nominal
measurement-noise covariances using characteristic scales in km and km/s.
Both scales must be finite and positive. These covariances describe sample
noise, not sensor bias or drift. Candidate
accelerometers in sensor-table rows 20-22 of `derivations/ekf_equations.tex`
are not implemented in the current measurement model.

Run `examples/tetherFilterDemo.m` in MATLAB for an example that propagates the
tethered three-body estimate and covariance for 24 hours, applying 72 noisy
measurement updates at 20-minute intervals with Joseph-form EKF updates. The
demo selects all eight implemented measurement groups. It creates one figure
per body, arranged as a 3-by-2 grid with x, y, and z component rows and
position and velocity columns. Each subplot shows the signed error with a
shaded componentwise +/-3-sigma bound.

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

3. Run MATLAB build checks locally (same commands used by hook tasks):

```powershell
pwsh -Command "matlab -batch \"buildtool('check')\""
```

4. Run full tests (same as pre-push):

```powershell
pwsh -Command "matlab -batch \"buildtool('test')\""
```

5. See repository developer docs and checklist:

- `docs/HOOKS.md` ??? hook behavior and installer notes (Windows specifics)
- `docs/COMMIT_CHECKLIST.md` ??? concise commit/push checklist (Windows-focused)

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

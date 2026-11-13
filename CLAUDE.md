# CLAUDE.md — coupled_orbit_attitude

MATLAB simulation code and LaTeX derivations for coupled orbit-attitude
dynamics: a tethered three-body system in the CR3BP, its measurement models,
and an (iterated) EKF with consistency diagnostics. README.md is the
authoritative user-facing reference; read it before changing public behavior.

## Layout

| Path | Purpose |
|---|---|
| `toolbox/+coa/` | Production code. Namespaces `coa.orbit`, `coa.meas`, `coa.filter`, `coa.util`. |
| `tests/unit/<package>/` | `matlab.unittest` class tests mirroring the toolbox packages (plus `examples/`). |
| `tests/regression/`, `tests/data/` | Fixed-seed short runs and their references (currently empty). |
| `derivations/` | LaTeX derivations (`main.tex` + topic files, `.bib`s). `nomenclature.tex` is the symbol ↔ code-name ↔ units source of truth. |
| `examples/` | Runnable demo scripts; `examples/+tetherFilterDemoSupport/` holds helpers for `tetherFilterDemo.m`. |
| `scratch/` | Exploratory work. Never called from `toolbox/`. |
| `test_results/` | Only `test_report.py` is tracked; XML/HTML outputs are generated and ignored. |
| `CHANGELOG.md`, `LOG.md` | Changelog (`[Unreleased]`) and dated lab log. |

## Build and test

Requires MATLAB R2024a+. From the repo root:

```matlab
buildtool                 % setup -> check (Code Analyzer on toolbox/) -> test
buildtool clean test      % force a fresh run (plain `buildtool test` may be skipped as up to date)
runtests('tests/unit', 'IncludeSubfolders', true)   % unit tests directly
```

- `setupTask` adds the whole repo (`genpath(root)`) to the path.
- Results go to `test_results/test-results.xml` and `coverage.xml` (coverage
  measured on `toolbox/` only). HTML report: `python test_results/test_report.py --open`
  (`--watch` for live updates, `--run` to run MATLAB first).
- MATLAB may not be available in every environment (e.g. cloud sessions). If you
  cannot run tests, say so explicitly rather than implying they passed.

## Hooks and linting

- `pre-commit` (via `.githooks/`, installed with `scripts/install_git_hooks.{sh,ps1}`):
  end-of-file / trailing-whitespace / YAML checks, then MISS_HIT `mh_style --fix`
  and `mh_lint` on staged `.m` files. `mh_style --fix` rewrites files, so re-stage
  after it runs.
- `pre-push`: `buildtool('test')`.
- No hosted CI (MATLAB licensing); hooks are the only automated gate.
- Run MISS_HIT locally: `pip install miss_hit==0.9.44`, then `mh_style --fix <files>`
  and `mh_lint <files>`.

## Code conventions

- `miss_hit.cfg`: 88-char lines, 4-space indent. Functions/scripts/methods
  `lowerCamelCase`, classes `PascalCase`, **parameters and attributes `snake_case`**
  (e.g. `length_dot`, `accel_fn`).
- Every source and test file starts with `% Copyright (c) 2026 Andrew Binder`.
- Function help: `NAME One-line summary.` followed by usage lines written with the
  fully qualified name (`coa.orbit.tetherAccel(...)`), inputs/outputs with sizes
  and units, and a pointer to the derivation section implemented.
- Validate inputs with `arguments` blocks (or `validateattributes`).
- Compute Jacobians/extra outputs only when requested (`if nargout > 1`).
- State conventions: nondimensional CR3BP rotating frame (rotates about +z at
  rate 1); measurement models use the 18-element state `[r1; r2; r3; v1; v2; v3]`;
  tether helpers use `[p1; p2; v1; v2]` or per-body `[p; v]` — check the
  function help before assuming an ordering.
- Filter code is observation-driven and deterministic: random draws for
  simulated measurements belong to the caller (demos, Monte Carlo), not to
  `coa.filter.*`.

## Derivation ↔ code traceability

- Each implemented equation has `\label{eq:...}` in `derivations/` and an
  `% IMPL: toolbox/...` comment beneath it; the code and its tests cite the label
  back (e.g. `% Tests for coa.util.skew (eq:skew)`). `grep -r "eq:name"` should
  find the equation, implementation, and tests.
- When adding or changing an implemented equation, update all three, and add
  any new symbols to `derivations/nomenclature.tex`.

## Workflow expectations

- New toolbox functionality needs tests in the matching `tests/unit/<package>/`
  folder. Promote code from `scratch/` into `toolbox/` only with tests.
- When behavior changes, update README.md, the relevant example, and the
  `[Unreleased]` section of `CHANGELOG.md` (see `docs/COMMIT_CHECKLIST.md`).
  Record notable experiments in `LOG.md` (date, what was tried, commit, output location).
- `examples/tetherFilterDemo.m` is configured through `config.plot_modes`
  (`"off"`/`"save"`/`"show"` per plot); use `"off"`/`"save"` for batch runs.
- Keep commits focused; commit messages are short imperative summaries
  (e.g. "Extract tether-aligned frame construction").
- Ask before consequential modeling, filter-design, or scientific choices.

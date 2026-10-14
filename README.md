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

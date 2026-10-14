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

Quick steps to get the repository ready for development and to enable automated checks/hooks.

1. Install local git hooks (PowerShell):

```powershell
pwsh scripts/install_git_hooks.ps1
```

2. Run the quick pre-commit checks locally (same as hooks):

```powershell
pwsh -Command "matlab -batch \"buildtool('check')\""
```

3. Run full tests (same as pre-push):

```powershell
pwsh -Command "matlab -batch \"buildtool('test')\""
```

4. See repository developer docs and checklist:

- `docs/HOOKS.md` — hook behavior and installer notes (Windows specifics)
- `docs/COMMIT_CHECKLIST.md` — concise commit/push checklist (Windows-focused)

CI notes

- The repository CI uses `matlab-actions/run-build@v3` to run `buildtool` on GitHub Actions and upload test artifacts.
- If a CI job fails with a MATLAB license error on GitHub-hosted runners, either add a base64-encoded license secret and write it to `$HOME/.matlab/<RELEASE>_licenses/license.dat` in the workflow, or run the job on a self-hosted runner with an activated MATLAB license. See `docs/COMMIT_CHECKLIST.md` for short troubleshooting steps.

If you need to bypass hooks for a one-off commit or push, use `--no-verify` with `git commit` or `git push`.

## Conventions

- Every implemented equation carries a `\label{eq:...}` in `derivations/` and an
  `% IMPL:` comment pointing to its code; the code cites the label back.
  `grep -r "eq:name"` finds the equation, its implementation, and its tests.
- Symbols and their code names are defined in `derivations/nomenclature.tex`.
- Promote work from `scratch/` into `toolbox/` only with tests.

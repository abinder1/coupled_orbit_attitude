# Commit & Push Checklist (Windows-focused)

This checklist is tuned for a single Windows developer workflow. Items marked **(automated)** are handled by the repository hooks or CI when hooks are installed.

## Pre-commit (runs automatically if hooks are installed)

- [x] Static analysis / code issues — **automated** via `buildtool('check')` (pre-commit hook)
- [ ] Update documentation and examples (manual)
- [ ] Update `CHANGELOG.md` / version metadata if applicable (manual)
- [ ] Add/update unit tests for new functionality (manual)
- [ ] Quick smoke-test main example(s) (recommended manual check)

## Pre-push (runs automatically if hooks are installed)

- [x] Run unit tests — **automated** via `buildtool('test')` (pre-push hook)
- [x] CI will run full build/check/test across configured MATLAB releases on push/PR — **automated**
- [ ] Validate release artifacts / packaging (manual)
- [ ] Confirm that examples and documentation build without errors (manual)

## CI / GitHub Actions

- The workflow at `.github/workflows/ci.yml` runs `buildtool` and uploads `test_results/` as artifacts.
- Resolve CI failures before merging PRs.

## Run checks locally (PowerShell — recommended on Windows)

Quick static checks (same as pre-commit):

```powershell
pwsh -Command "matlab -batch \"buildtool('check')\""
```

Run full tests (same as pre-push):

```powershell
pwsh -Command "matlab -batch \"buildtool('test')\""
```

## Install hooks (once per clone)

Run the PowerShell installer to copy the hooks into `.git/hooks`:

```powershell
pwsh scripts/install_git_hooks.ps1
```

If you use Git Bash or a Unix-like shell on Windows and need the hooks to be executable, run:

```bash
chmod +x .git/hooks/pre-commit .git/hooks/pre-push
```

## Bypass hooks (if necessary)

To skip hooks for a commit or push:

```powershell
git commit -m "..." --no-verify
git push --no-verify
```

## Notes

- The PowerShell installer copies hook scripts but does not attempt to set Git-tracked executable bits (that operation can fail on Windows). The hooks will run when executable in your shell environment (Git Bash) or when invoked by Git on Windows.
- Keep the checklist items up-to-date when adding new CI steps, tests, or documentation requirements.

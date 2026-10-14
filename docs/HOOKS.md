# Git hooks and CI — install and usage

This repository includes local git hook scripts under `.githooks/` and a GitHub Actions workflow at `.github/workflows/ci.yml`.

Install the local hooks (once per clone):

PowerShell (Windows):

```powershell
pwsh scripts/install_git_hooks.ps1
```

Git Bash / macOS / Linux:

```bash
bash scripts/install_git_hooks.sh
```

What the hooks do
- `pre-commit`: runs `buildtool('check')` for quick static checks.
- `pre-push`: runs `buildtool('test')` (full test run; can be slow).

If you need faster checks, edit `.githooks/pre-commit` to run a lighter command.

CI
- The GitHub Actions workflow runs `buildtool` on push and PRs and uploads `test_results/` as artifacts.

Windows notes
- The PowerShell installer (`scripts/install_git_hooks.ps1`) copies hook scripts into the repository `.git/hooks` directory but does not try to set the executable bit via `git update-index` (this step can fail on Windows and for files inside `.git/hooks`).
- If you use Git Bash or another Unix-like shell on Windows and need the hooks to be executable, run:

```bash
chmod +x .git/hooks/pre-commit .git/hooks/pre-push
```

- It's normal for the installer to skip `git update-index` for `.git/hooks`; the hooks will still run when they are executable in your shell environment. If you prefer not to run hooks locally, commit with `--no-verify`.

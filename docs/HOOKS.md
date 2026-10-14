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

CI (local-only)

This repository no longer uses GitHub Actions for CI. Instead, local git hooks implement automated checks at commit and push time.

- `pre-commit` runs lightweight static checks (`buildtool('check')`) on staged MATLAB files.
- `pre-push` runs the full test task (`buildtool('test')`) before pushing (can be adjusted if you want faster checks).

If you prefer remote CI in the future, re-add a workflow and provide either a license secret for hosted runners or a self-hosted runner with an activated MATLAB license.

Windows notes
- The PowerShell installer (`scripts/install_git_hooks.ps1`) copies hook scripts into the repository `.git/hooks` directory but does not try to set the executable bit via `git update-index` (this step can fail on Windows and for files inside `.git/hooks`).
- If you use Git Bash or another Unix-like shell on Windows and need the hooks to be executable, run:

```bash
chmod +x .git/hooks/pre-commit .git/hooks/pre-push
```

- It's normal for the installer to skip `git update-index` for `.git/hooks`; the hooks will still run when they are executable in your shell environment. If you prefer not to run hooks locally, commit with `--no-verify`.

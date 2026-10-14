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

CI notes
- This repository uses the `matlab-actions/run-build@v3` workflow to run the project's `buildtool` tasks on GitHub Actions. That action performs MATLAB setup and runs the requested tasks (for example, `tasks: test`).
- Common CI failure: "Unable to find a license for MATLAB" — hosted GitHub runners may have MATLAB installed but require a license to run jobs. Two ways to resolve this:
	- Provide a license via a repository secret: store a base64-encoded `license.dat` in a secret (e.g. `MATLAB_LICENSE_B64`) and write it to `$HOME/.matlab/<RELEASE>_licenses/license.dat` in a workflow step before `run-build`.
	- Use a self-hosted runner that already has MATLAB installed and activated under your license (set `runs-on: self-hosted` or a custom label).
- If CI fails with license errors, check the Actions run log to identify the failing step (Setup vs Run) and follow the above fixes. See `docs/COMMIT_CHECKLIST.md` for concise troubleshooting steps.

Windows notes
- The PowerShell installer (`scripts/install_git_hooks.ps1`) copies hook scripts into the repository `.git/hooks` directory but does not try to set the executable bit via `git update-index` (this step can fail on Windows and for files inside `.git/hooks`).
- If you use Git Bash or another Unix-like shell on Windows and need the hooks to be executable, run:

```bash
chmod +x .git/hooks/pre-commit .git/hooks/pre-push
```

- It's normal for the installer to skip `git update-index` for `.git/hooks`; the hooks will still run when they are executable in your shell environment. If you prefer not to run hooks locally, commit with `--no-verify`.

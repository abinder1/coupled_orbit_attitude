# Git hooks (quick reference)

This repo uses local git hooks for automated checks.

## One-time setup

Windows (PowerShell):

```powershell
pwsh scripts/install_git_hooks.ps1
py -m pip install --user pre-commit
py -m pre_commit run --all-files
```

macOS/Linux (Bash):

```bash
bash scripts/install_git_hooks.sh
python3 -m pip install --user pre-commit
pre-commit run --all-files
```

## Hook behavior

- `pre-commit`: runs `.pre-commit-config.yaml` hooks
  - `mh_style --fix` on staged `.m` files
  - `mh_lint` on staged `.m` files
- `pre-push`: runs `buildtool('test')`

## Notes

- MISS_HIT rules are configured in `miss_hit.cfg`.
- Bypass hooks when necessary with `--no-verify`.

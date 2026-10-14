# Commit & Push Checklist (Windows)

Use this as a quick pre-flight list before committing.

## Automated (when hooks are installed)

- `pre-commit`: `mh_style --fix` and `mh_lint` on staged `.m` files
- `pre-push`: `buildtool('test')`

## Manual checks

- [ ] Update docs/examples if behavior changed
- [ ] Add or update tests for new functionality
- [ ] Update `CHANGELOG.md` / version metadata when applicable
- [ ] Smoke-test key example(s) when behavior changed

## One-time setup (per clone)

```powershell
pwsh scripts/install_git_hooks.ps1
py -m pip install --user pre-commit
py -m pre_commit run --all-files
```

## Run checks manually (optional)

```powershell
pwsh -Command "matlab -batch \"buildtool('check')\""
pwsh -Command "matlab -batch \"buildtool('test')\""
```

## Emergency bypass

```powershell
git commit -m "..." --no-verify
git push --no-verify
```

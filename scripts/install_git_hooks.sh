#!/usr/bin/env bash
# Install git hooks from .githooks into .git/hooks (Unix / Git Bash)
set -e
repo_root="$(git rev-parse --show-toplevel)"
src="$repo_root/.githooks"
dest="$repo_root/.git/hooks"
mkdir -p "$dest"
for f in "$src"/*; do
  cp -f "$f" "$dest/$(basename "$f")"
  chmod +x "$dest/$(basename "$f")"
done
echo "Git hooks installed."

# Install git hooks from .githooks into .git/hooks (Windows / PowerShell)
param()
try {
    $repo = git rev-parse --show-toplevel
} catch {
    Write-Error "Not a git repository or git not found."
    exit 1
}
Push-Location $repo
$src = Join-Path $repo '.githooks'
$dest = Join-Path $repo '.git\hooks'
if (-not (Test-Path $src)) {
    Write-Error ".githooks not found in $repo"
    Pop-Location
    exit 1
}
Get-ChildItem -Path $src -File | ForEach-Object {
    Copy-Item -Path $_.FullName -Destination (Join-Path $dest $_.Name) -Force
    # Do not call `git update-index` on files inside .git/hooks: these files are not tracked
    # and attempting to set the executable bit via git on Windows can produce errors.
    # If running on a Unix-like shell, make the hook executable as a courtesy.
    try {
        if ($IsLinux -or $IsMacOS) {
            & chmod +x (Join-Path $dest $_.Name) 2>$null
        }
    } catch {
        # ignore chmod failures
    }
}
Pop-Location
Write-Output "Git hooks installed."

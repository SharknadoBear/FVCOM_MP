param(
    [string]$VenvName = ".venv"
)

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$venvPath = Join-Path $scriptDir $VenvName

python -m venv --system-site-packages $venvPath

Write-Host "Created virtual environment at $venvPath"
Write-Host "Activate with:"
Write-Host "  $venvPath\\Scripts\\Activate.ps1"

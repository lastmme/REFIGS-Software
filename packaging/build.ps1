param([string]$Version = '1.0.0')
$ErrorActionPreference = 'Stop'
Set-Location (Split-Path $PSScriptRoot -Parent)
$python = Join-Path (Get-Location) '.venv\Scripts\python.exe'
& $python -m PyInstaller --noconfirm --windowed --onedir --name RE-FIGS --icon "$PWD\REFIGS.ico" --distpath out/dist --workpath out/build --specpath out --collect-data qfluentwidgets --collect-all pyteomics --collect-all cvxopt --hidden-import packaging_smoke --exclude-module tkinter main.py
if ($LASTEXITCODE -ne 0) { throw 'Application build failed' }
$report = Join-Path $PWD ('out/smoke-' + [guid]::NewGuid().ToString('N') + '.json')
$application = Join-Path $PWD 'out/dist/RE-FIGS/RE-FIGS.exe'
$smoke = Start-Process -FilePath $application -ArgumentList @('--smoke-test', ('"' + $report + '"')) -WindowStyle Hidden -PassThru
$deadline = (Get-Date).AddSeconds(180)
while (-not $smoke.WaitForExit(1000)) {
    if ((Get-Date) -gt $deadline) { $smoke.Kill(); throw 'Application smoke test timed out' }
}
if ($smoke.ExitCode -ne 0 -or -not (Test-Path $report)) { throw 'Application smoke test failed; inspect %LOCALAPPDATA%/RE-FIGS/startup.log' }
Get-Content $report
$compiler = @("${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe", "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $compiler) { throw 'Install Inno Setup 6 first' }
& $compiler "/DAppVersion=$Version" packaging/installer.iss
if ($LASTEXITCODE -ne 0) { throw 'Installer build failed' }
Compress-Archive -Path out/dist/RE-FIGS -DestinationPath "out/RE-FIGS-$Version-windows-x64-portable.zip" -Force
Get-FileHash "out/RE-FIGS-$Version-windows-x64-setup.exe", "out/RE-FIGS-$Version-windows-x64-portable.zip" -Algorithm SHA256 | Format-Table -AutoSize
$lines = Get-FileHash "out/RE-FIGS-$Version-windows-x64-setup.exe", "out/RE-FIGS-$Version-windows-x64-portable.zip" -Algorithm SHA256 | ForEach-Object { $_.Hash.ToLowerInvariant() + '  ' + (Split-Path $_.Path -Leaf) }
[IO.File]::WriteAllLines((Join-Path $PWD 'out/SHA256SUMS.txt'), $lines, (New-Object Text.UTF8Encoding($false)))

# Copies the three examples into a CHAI3D tree and registers them in examples\GLFW\CMakeLists.txt.
# Usage (PowerShell):  .\scripts\install.ps1 -Chai3dRoot "E:\chai3d-master\chai3d-master"
# Then: cd <chai3d>\build ; cmake .. ; cmake --build . --config Release --target 44-soft-organ-scalpel
param([string]$Chai3dRoot = "E:\chai3d-master\chai3d-master")

$ex   = Join-Path $Chai3dRoot "examples\GLFW"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = Split-Path -Parent $here
if (-not (Test-Path $ex)) { Write-Error "Not found: $ex"; exit 1 }

$names = @('42-bone-probe', '43-soft-organ-probe', '44-soft-organ-scalpel')
foreach ($n in $names) {
    $dst = Join-Path $ex $n
    New-Item -ItemType Directory -Path $dst -Force | Out-Null
    Get-ChildItem (Join-Path $repo "examples\$n") -File |
        Where-Object { $_.Extension -in '.h', '.cpp', '.bin' } |
        Copy-Item -Destination $dst -Force
    Write-Host "copied $n"
}

$cm = Join-Path $ex 'CMakeLists.txt'
$t  = Get-Content $cm -Raw
$missing = $names | Where-Object { $t -notmatch [regex]::Escape($_) }
if ($missing.Count -eq 0) {
    Write-Host "CMakeLists.txt already lists all three examples."
} elseif ($t -match '40-surgical-sim') {
    Copy-Item $cm "$cm.bak" -Force
    $re = [regex]::new('40-surgical-sim')
    $t  = $re.Replace($t, ('40-surgical-sim ' + ($missing -join ' ')), 1)
    Set-Content $cm $t -NoNewline
    Write-Host "added to CMakeLists.txt: $($missing -join ', ')  (backup: CMakeLists.txt.bak)"
} else {
    Write-Warning "Could not find 40-surgical-sim in CMakeLists.txt. Add these names to the example list by hand: $($missing -join ' ')"
}
Write-Host "Next: cd $Chai3dRoot\build ; cmake .. ; cmake --build . --config Release --target <example name>"

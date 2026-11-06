# CRISPRme+ one-line installer (Windows).
#   irm https://pinellolab.github.io/crisprme-plus/install.ps1 | iex
#
# Installs a clickable "CRISPRme" app (3 big buttons: Start / Update / Stop). The
# first Start downloads the reference + variant data automatically, then opens the
# web interface — so after this one command the user never needs the terminal
# again. Docker Desktop required.
$ErrorActionPreference = 'Stop'
$Image   = 'pinellolab/crisprme:latest'
$DataDir = Join-Path $env:USERPROFILE 'CRISPRme'
$AppDir  = Join-Path $env:LOCALAPPDATA 'CRISPRme'

Write-Host "=== CRISPRme+ installer ===" -ForegroundColor Cyan

# ---- 1. Docker check -------------------------------------------------------
if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
  Write-Host "Docker is not installed. CRISPRme runs inside Docker." -ForegroundColor Yellow
  Write-Host "Please install Docker Desktop, then re-run this installer:"
  Write-Host "  https://www.docker.com/products/docker-desktop/" -ForegroundColor Blue
  Start-Process "https://www.docker.com/products/docker-desktop/"
  return
}
Write-Host "Docker found." -ForegroundColor Green

# ---- 1b. Pre-pull the engine image (best-effort, visible progress) --------
# Front-load the engine image now (a ~5 GB download: Docker Desktop on Windows pulls the
# linux/amd64 one, ~7 GB+ once unpacked) with progress here, so the app's first Start
# doesn't appear to hang while it pulls. Skipped if the daemon isn't up yet.
try { docker info *> $null; if ($LASTEXITCODE -eq 0) { Write-Host "Fetching the CRISPRme engine image (~5 GB download, one time)..." -ForegroundColor Cyan; docker pull $Image } } catch {}

# ---- 2. Config dir --------------------------------------------------------
# The app stores its data (a ~70 GB download, ~150 GB once installed) in a folder the user
# picks on first run and remembers in $DataDir\.data_location (so a low-space home disk can
# point the data at another drive). $DataDir itself only holds that tiny pointer file.
New-Item -ItemType Directory -Force -Path $DataDir | Out-Null
Write-Host "Config folder ready: $DataDir  (data location chosen on first Start)" -ForegroundColor Green

# ---- 3. Write the clickable app (3 big buttons) ---------------------------
New-Item -ItemType Directory -Force -Path $AppDir | Out-Null
$AppScript = Join-Path $AppDir 'CRISPRme.ps1'
@'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$ConfigDir = Join-Path $env:USERPROFILE 'CRISPRme'
$LocFile   = Join-Path $ConfigDir '.data_location'
$Image     = 'pinellolab/crisprme:latest'

# Data (a ~70 GB download, ~150 GB once installed) lives in a user-chosen folder, remembered
# across launches in .data_location (so a low-space home disk can point the data at another drive).
function Get-DataDir { if (Test-Path $LocFile) { (Get-Content -Raw $LocFile).Trim() } else { Join-Path $ConfigDir 'crisprme-data' } }
function Set-DataDir($p) { New-Item -ItemType Directory -Force -Path $ConfigDir | Out-Null; Set-Content -NoNewline -Encoding ascii -Path $LocFile -Value $p }
function Has-Stored { Test-Path $LocFile }
# first run: let the user choose where the data lives (e.g. a big external drive)
function Pick-DataDir {
  $def = Join-Path $ConfigDir 'crisprme-data'
  $r = [System.Windows.Forms.MessageBox]::Show("Where should CRISPRme store its data (a ~70 GB download, ~150 GB once installed: pick a drive with ~200 GB free)?`n`nDefault: $def`n`nClick Yes to choose another drive/folder (if your home disk is low on space), or No to use the default.","CRISPRme+ - data location",'YesNo')
  if ($r -eq 'Yes') {
    $fb = New-Object System.Windows.Forms.FolderBrowserDialog
    $fb.Description = 'Select a folder for the CRISPRme data (~150 GB once installed; keep ~200 GB free)'
    if ($fb.ShowDialog() -eq 'OK' -and $fb.SelectedPath) { $p = Join-Path $fb.SelectedPath 'crisprme-data' } else { $p = $def }
  } else { $p = $def }
  New-Item -ItemType Directory -Force -Path $p | Out-Null
  Set-DataDir $p
  return $p
}
$Data = Get-DataDir

# The first-time setup marks the data folder until both of its downloads have ended well, so
# a Start after one failed sets up again (download picks up what is there) instead of opening
# the app on part of the data: Genomes/ is there from the first minutes of a download.
function Setup-Mark { Join-Path $Data '.setup_unfinished' }
function Has-Data { (Test-Path (Join-Path $Data 'Genomes')) -and -not (Test-Path (Setup-Mark)) }
function Docker-Running { try { docker info *> $null; return ($LASTEXITCODE -eq 0) } catch { return $false } }
# a bind mount actually works (freshly-started Docker can report info OK while its
# file-sharing is still initializing -> the download's makedirs fails with ENOENT)
function Mount-Ready {
  try { docker run --rm -v "${Data}:/DATA" $Image bash -c 'touch /DATA/.probe && rm -f /DATA/.probe' *> $null; return ($LASTEXITCODE -eq 0) } catch { return $false }
}
function Ensure-Docker {
  if (-not (Docker-Running)) {
    $dd = Join-Path $env:ProgramFiles 'Docker\Docker\Docker Desktop.exe'
    if (Test-Path $dd) { Start-Process $dd } else { [System.Windows.Forms.MessageBox]::Show('Docker Desktop does not seem to be installed. Install it from https://www.docker.com/products/docker-desktop/ then try again.'); return $false }
    $up = $false
    for ($i=0; $i -lt 30; $i++) { Start-Sleep 3; if (Docker-Running) { $up = $true; break } }
    if (-not $up) { [System.Windows.Forms.MessageBox]::Show('Docker Desktop is still starting. Give it a few more seconds, then click Start again.'); return $false }
  }
  for ($i=0; $i -lt 20; $i++) { if (Mount-Ready) { return $true }; Start-Sleep 2 }
  return $true
}
function Run-InConsole($title, $cmd) {
  Start-Process powershell -ArgumentList @('-NoExit','-Command',
    "Write-Host '== $title =='; $cmd; Write-Host ''; Write-Host 'DONE - you can close this window.'")
}

$form = New-Object System.Windows.Forms.Form
$form.Text = 'CRISPRme+'
$form.Size = New-Object System.Drawing.Size(440,320)
$form.StartPosition = 'CenterScreen'; $form.FormBorderStyle = 'FixedSingle'; $form.MaximizeBox = $false
$status = New-Object System.Windows.Forms.Label
$status.SetBounds(20,12,400,24); $form.Controls.Add($status)
$big = New-Object System.Drawing.Font('Segoe UI',15,[System.Drawing.FontStyle]::Bold)

function Mk($text,$y,$action) {
  $b = New-Object System.Windows.Forms.Button
  $b.SetBounds(30,$y,380,64); $b.Text = $text; $b.Font = $big; $b.Add_Click($action)
  $form.Controls.Add($b)
}
Mk 'Start' 44 {
  # first run with no data AND no remembered location: let the user pick the folder
  if ((-not (Has-Data)) -and (-not (Has-Stored))) { $script:Data = Pick-DataDir }
  if (-not (Ensure-Docker)) { return }
  # Auto-use an NVIDIA GPU when Docker can expose one (no-op on CPU-only hosts): probe
  # '--gpus all' and, if it works, run the scorer on the GPU (cuda backend). The
  # CRISPR-Bulge scorer still CPU-falls-back internally if the GPU turns out unusable.
  $gpu = ''
  try { docker run --rm --gpus all $Image true *> $null; if ($LASTEXITCODE -eq 0) { $gpu = '--gpus all -e CRISPRME_COMPUTE_BACKEND=cuda ' } } catch {}
  $run = "docker rm -f crisprme 2>`$null; docker run -d --name crisprme ${gpu}-v `"${Data}:/DATA`" -w /DATA -p 8080:8080 $Image crisprme.py web-interface"
  if (Has-Data) {
    Invoke-Expression $run; Start-Sleep 3; Start-Process 'http://localhost:8080'
  } else {
    $first = if (Test-Path (Setup-Mark)) { 'The first-time setup did not finish: CRISPRme will download what is missing of' } else { 'First run: CRISPRme will download' }
    $r = [System.Windows.Forms.MessageBox]::Show("$first the reference genome + HGDP + 1000 Genomes variant data (~70 GB) into:`n$Data`n`nOnce only. A window shows progress, then your browser opens automatically.","CRISPRme+ - first-time setup",'OKCancel')
    if ($r -eq 'OK') {
      $mnt = "-v `"${Data}:/DATA`" -w /DATA $Image"
      $mark = Setup-Mark
      New-Item -ItemType File -Force -Path $mark | Out-Null
      Run-InConsole 'First-time setup - downloading ~70 GB, then starting CRISPRme' "docker run --rm $mnt crisprme.py download --what all --path /DATA; if (`$LASTEXITCODE -eq 0) { docker run --rm $mnt crisprme.py download --what index --index-name NRR_3_hg38+hg38_HGDP1kGP --path /DATA }; if (`$LASTEXITCODE -eq 0) { Remove-Item -Force -LiteralPath `"$mark`"; $run; Start-Sleep 3; Start-Process 'http://localhost:8080' } else { Write-Host 'The download did not finish: see above, then click Start again to finish setting up.' -ForegroundColor Red }"
    }
  }
}
Mk 'Update' 120 { if (Ensure-Docker) { Run-InConsole 'Updating CRISPRme' "docker pull $Image" } }
Mk 'Stop' 196 { try { docker rm -f crisprme } catch {}; [System.Windows.Forms.MessageBox]::Show('CRISPRme stopped.') }

if (Has-Data) { $status.Text = 'Ready. Click Start to open CRISPRme.' } elseif (Test-Path (Setup-Mark)) { $status.Text = 'Setup did not finish: Start finishes it.' } else { $status.Text = 'First run: Start will set everything up (a ~70 GB download, once).' }
[void]$form.ShowDialog()
'@ | Set-Content -Encoding UTF8 $AppScript

# ---- 4. Shortcuts (Desktop + Start Menu) ----------------------------------
$ws = New-Object -ComObject WScript.Shell
foreach ($dir in @([Environment]::GetFolderPath('Desktop'),
                   (Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs'))) {
  $lnk = $ws.CreateShortcut((Join-Path $dir 'CRISPRme.lnk'))
  $lnk.TargetPath = 'powershell.exe'
  $lnk.Arguments  = "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File `"$AppScript`""
  $lnk.WorkingDirectory = $DataDir
  $lnk.IconLocation = 'powershell.exe,0'
  $lnk.Save()
}

Write-Host ""
Write-Host "Done! Open 'CRISPRme' from your Desktop or Start Menu and click Start." -ForegroundColor Green
Write-Host "The first Start downloads the data (once) and opens the web app."

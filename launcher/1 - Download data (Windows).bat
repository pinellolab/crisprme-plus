@echo off
REM Double-click this ONCE to download the CRISPRme+ reference data + index.
REM Requires Docker Desktop installed and running.
cd /d "%~dp0"
set IMG=pinellolab/crisprme:v2.4.0

echo ===================================================================
echo   CRISPRme+  -  one-time data download
echo ===================================================================
docker info >nul 2>&1
if errorlevel 1 (
  echo.
  echo   Docker Desktop does not appear to be running.
  echo   Open Docker Desktop, wait for it to start, then run this again.
  echo.
  pause
  exit /b 1
)

if not exist crisprme-data mkdir crisprme-data
echo.
echo Downloading the reference genome, annotations and reference index
echo (~25 GB, one time - this can take a while)...
echo.
docker run --rm -v "%cd%/crisprme-data:/DATA" -w /DATA %IMG% crisprme.py download --what all --path /DATA

echo.
echo -------------------------------------------------------------------
echo   Reference data is ready - reference-only searches will work now.
echo.
echo   For VARIANT-AWARE search (1000G + HGDP; needs ~64 GB RAM in Docker
echo   Desktop and ~85 GB free disk), also double-click:
echo     "1b - Download variant index (Windows).bat"
echo -------------------------------------------------------------------
echo.
echo   Next: double-click  "2 - Start CRISPRme (Windows).bat"
echo.
pause

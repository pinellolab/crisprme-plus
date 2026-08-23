@echo off
REM OPTIONAL, one time: download the variant-aware (1000G + HGDP) index.
REM Needs ~64 GB RAM in Docker Desktop and ~85 GB free disk. Run "1 - Download data" first.
cd /d "%~dp0"
set IMG=pinellolab/crisprme:v2.4.0

echo ===================================================================
echo   CRISPRme+  -  variant index download (1000G + HGDP)
echo ===================================================================
docker info >nul 2>&1
if errorlevel 1 (
  echo.
  echo   Docker Desktop is not running - open it and try again.
  echo.
  pause
  exit /b 1
)
if not exist crisprme-data (
  echo.
  echo   Please run "1 - Download data" first.
  echo.
  pause
  exit /b 1
)

echo.
echo Downloading the dict-less 1000G + HGDP variant index (~16 GB download,
echo expands to ~60 GB). One time - this can take a while...
echo.
docker run --rm -v "%cd%/crisprme-data:/DATA" -w /DATA %IMG% crisprme.py download --what index --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path /DATA

echo.
echo   Variant index ready. In Docker Desktop, set Settings -^> Resources -^>
echo   Memory to 64 GB before running a genome-wide variant search.
echo.
echo   Next: double-click  "2 - Start CRISPRme (Windows).bat"
echo.
pause

@echo off
REM Double-click to START CRISPRme+. Your browser opens at http://localhost:8080.
REM Keep the window that appears OPEN while you use CRISPRme; close it to stop.
REM Requires Docker Desktop running and "1 - Download data" done.
cd /d "%~dp0"

echo ===================================================================
echo   CRISPRme+  -  starting the web interface
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
  echo   No data found. Please run "1 - Download data" first.
  echo.
  pause
  exit /b 1
)

echo.
echo   Starting... your browser will open at  http://localhost:8080
echo   KEEP THIS WINDOW OPEN while you use CRISPRme.
echo   To stop CRISPRme, close this window.
echo.

REM open the browser shortly after the server starts
start "" cmd /c "timeout /t 8 >nul & start http://localhost:8080"

docker compose version >nul 2>&1
if errorlevel 1 (
  docker-compose up
) else (
  docker compose up
)

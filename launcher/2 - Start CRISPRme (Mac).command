#!/bin/bash
# Double-click to START CRISPRme+. Your browser opens at http://localhost:8080.
# Keep the window that appears OPEN while you use CRISPRme; close it (or press
# Ctrl-C) to stop. Requires Docker Desktop running and '1 - Download data' done.
cd "$(dirname "$0")" || exit 1

echo "==================================================================="
echo "  CRISPRme+  —  starting the web interface"
echo "==================================================================="
if ! docker info >/dev/null 2>&1; then
  echo; echo "  Docker Desktop is not running — open it and try again."; echo
  read -r -p "Press Enter to close."; exit 1
fi
if [ ! -d crisprme-data ]; then
  echo; echo "  No data found. Please run '1 - Download data' first."; echo
  read -r -p "Press Enter to close."; exit 1
fi

echo
echo "  Starting... your browser will open at  http://localhost:8080"
echo "  KEEP THIS WINDOW OPEN while you use CRISPRme."
echo "  To stop CRISPRme, close this window (or press Ctrl-C)."
echo

# open the browser a few seconds after the server comes up
( sleep 8; open "http://localhost:8080" >/dev/null 2>&1 ) &

# 'docker compose' (v2) or the legacy 'docker-compose'
if docker compose version >/dev/null 2>&1; then
  docker compose up
else
  docker-compose up
fi

#!/bin/bash
# Double-click this ONCE to download the CRISPRme+ reference data + index.
# Requires Docker Desktop to be installed and running.
cd "$(dirname "$0")" || exit 1
IMG="pinellolab/crisprme:v2.4.0"

echo "==================================================================="
echo "  CRISPRme+  —  one-time data download"
echo "==================================================================="
if ! docker info >/dev/null 2>&1; then
  echo
  echo "  Docker Desktop does not appear to be running."
  echo "  Please open Docker Desktop, wait for it to start, then run this again."
  echo
  read -r -p "Press Enter to close."
  exit 1
fi

mkdir -p crisprme-data
echo
echo "Downloading the reference genome, annotations and the reference index"
echo "(~25 GB, one time — this can take a while on the first run)..."
echo
docker run --rm -v "$(pwd)/crisprme-data:/DATA" -w /DATA "$IMG" \
  crisprme.py download --what all --path /DATA

echo
echo "-------------------------------------------------------------------"
echo "  Reference data is ready — reference-only searches will work now."
echo
echo "  For VARIANT-AWARE search (1000G + HGDP; needs ~64 GB RAM in Docker"
echo "  Desktop and ~85 GB free disk), also download the variant index by"
echo "  double-clicking:  '1b - Download variant index (Mac).command'"
echo "-------------------------------------------------------------------"
echo
echo "  Next: double-click  '2 - Start CRISPRme (Mac).command'"
echo
read -r -p "Press Enter to close."

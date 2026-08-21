#!/usr/bin/env bash
#
# retrieve_AoU.sh
# ===============
# Build the All of Us (AoU) aggregate variant dataset for CRISPRme+ by scraping
# the AoU public Data Browser API and turning the results into MAF-filtered VCFs.
# Unlike 1000G / HGDP / gnomAD, AoU offers no bulk VCF download; the aggregate
# allele frequencies are only reachable through the Data Browser's public API, so
# this "retrieval" is a scrape-and-build pipeline. It runs four phases per
# chromosome:
#
#   1. SCRAPE   the AoU variant API window-by-window into raw JSON pages
#              (scrape_aou.py; resumable, rate-limited).
#   2. CSV      combine the JSON pages into one de-duplicated CSV
#              (json_to_csv.py).
#   3. VCF      convert the CSV into an aggregate VCF with a single `AllOfUs`
#              pseudo-sample and AF/AC/AN/HOM INFO (csv_to_vcf.py), then
#              `bcftools sort` + `tabix` -> full/.
#   4. FILTER   keep variants with INFO/AF > MAF (default 0.001) -> more_than_001/
#              (the same threshold used for gnomAD / TOPMed).
#
# Dataset .......... All of Us Research Program, aggregate allele frequencies
# Source ........... AoU public Data Browser API
#                    https://public.api.researchallofus.org/v1/genomics/search-variants
# Reference build .. GRCh38 / hg38 (AoU Data Browser coordinates)
# Data type ........ Aggregate, sites-only — single pseudo-sample, no genotypes
#
# Outputs (under $OUTDIR)
# -----------------------
#   json/<chrom>/file_*.json                          (1: raw scraped pages)
#   csv/<chrom>.csv                                   (2: combined CSV)
#   full/<chrom>.allofus.unphased.sorted.vcf.gz       (3: full VCF + .tbi)
#   more_than_001/<chrom>.allofus.unphased.sorted.vcf.gz  (4: MAF-filtered + .tbi)
#
# The pipeline is resumable and idempotent: the scrape checkpoints its progress,
# and each phase is rebuilt only when its input is newer than its output (or
# FORCE=1). NOTE: scraping the whole genome is SLOW (the API is paginated and
# politely rate-limited); run inside tmux/screen and expect many hours.
#
# Dependencies
# ------------
#   * python3 with `requests`, `pandas`, `pysam` (the three stage scripts).
#   * bcftools + tabix (sort / filter / index).
#   * a GRCh38/hg38 reference FASTA (REFERENCE), for the VCF contig set.
#
# Usage
# -----
#   REFERENCE=/path/to/hg38.fa ./retrieve_AoU.sh [chrom ...]
#
#   With no chromosome arguments all 24 are built (chr1-22, chrX, chrY). Pass one
#   or more labels to build only those, e.g.:
#       REFERENCE=hg38.fa ./retrieve_AoU.sh chr22
#
# Configuration (environment variables)
# -------------------------------------
#   REFERENCE      REQUIRED — path to the GRCh38/hg38 reference FASTA.
#   OUTDIR         destination directory (default: ./hg38_AoU).
#   MAF            allele-frequency threshold; keep INFO/AF > MAF (default: 0.001).
#   FILTER_SUBDIR  subfolder for the MAF-filtered VCFs (default: more_than_001).
#   WINDOW         scrape window size in bp (default: 1000000).
#   ROWCOUNT       rows per API page (default: 50000).
#   FORCE          set to 1 to rebuild all phases regardless of freshness.
#   PYTHON         python interpreter (default: python3).
#   BCFTOOLS       bcftools command (default: bcftools).
#   TABIX          tabix command (default: tabix).
#
# Exit status: 0 on success; non-zero if a phase failed for any chromosome.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
REFERENCE="${REFERENCE:-}"
OUTDIR="${OUTDIR:-hg38_AoU}"
MAF="${MAF:-0.001}"
FILTER_SUBDIR="${FILTER_SUBDIR:-more_than_001}"
WINDOW="${WINDOW:-1000000}"
ROWCOUNT="${ROWCOUNT:-50000}"
FORCE="${FORCE:-0}"
PYTHON="${PYTHON:-python3}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"

ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# ------------------------------------------------------------------------------
# GRCh38 / hg38 primary-assembly chromosome lengths (bp).
# Only used to bound the scrape windows; an exact value is not critical (the API
# simply returns nothing past the real chromosome end).
# ------------------------------------------------------------------------------
chrom_length() {
    case "$1" in
        chr1)  echo 248956422 ;;
        chr2)  echo 242193529 ;;
        chr3)  echo 198295559 ;;
        chr4)  echo 190214555 ;;
        chr5)  echo 181538259 ;;
        chr6)  echo 170805979 ;;
        chr7)  echo 159345973 ;;
        chr8)  echo 145138636 ;;
        chr9)  echo 138394717 ;;
        chr10) echo 133797422 ;;
        chr11) echo 135086622 ;;
        chr12) echo 133275309 ;;
        chr13) echo 114364328 ;;
        chr14) echo 107043718 ;;
        chr15) echo 101991189 ;;
        chr16) echo 90338345 ;;
        chr17) echo 83257441 ;;
        chr18) echo 80373285 ;;
        chr19) echo 58617616 ;;
        chr20) echo 64444167 ;;
        chr21) echo 46709983 ;;
        chr22) echo 50818468 ;;
        chrX)  echo 156040895 ;;
        chrY)  echo 57227415 ;;
        *)     return 1 ;;
    esac
}

# ------------------------------------------------------------------------------
# Prerequisite checks
# ------------------------------------------------------------------------------
check_prereqs() {
    local missing=0
    command -v "$PYTHON" >/dev/null 2>&1 || { echo "ERROR: python ('$PYTHON') not found" >&2; missing=1; }
    "$PYTHON" -c "import requests, pandas, pysam" 2>/dev/null || {
        echo "ERROR: python packages requests/pandas/pysam are required" >&2; missing=1; }
    command -v "$BCFTOOLS" >/dev/null 2>&1 || { echo "ERROR: bcftools ('$BCFTOOLS') not found" >&2; missing=1; }
    command -v "$TABIX" >/dev/null 2>&1 || { echo "ERROR: tabix ('$TABIX') not found" >&2; missing=1; }
    if [ -z "$REFERENCE" ]; then
        echo "ERROR: REFERENCE is required — set REFERENCE=/path/to/hg38.fa" >&2; missing=1
    elif [ ! -f "$REFERENCE" ]; then
        echo "ERROR: reference not found: $REFERENCE" >&2; missing=1
    fi
    [ "$missing" -eq 0 ]
}

# ------------------------------------------------------------------------------
# Per-phase helpers (each rebuilds only when stale, unless FORCE=1)
# ------------------------------------------------------------------------------
stale() {
    # stale <output> <input>: true if FORCE, or output missing, or input newer.
    [ "$FORCE" = "1" ] && return 0
    [ -f "$1" ] || return 0
    [ "$2" -nt "$1" ]
}

build_chrom() {
    local chrom="$1" length="$2"
    local jdir="$OUTDIR/json/$chrom"
    local csv="$OUTDIR/csv/$chrom.csv"
    local full="$OUTDIR/full/$chrom.allofus.unphased.sorted.vcf.gz"
    local filt="$OUTDIR/$FILTER_SUBDIR/$chrom.allofus.unphased.sorted.vcf.gz"
    local complete="$jdir/.complete"

    echo "===== $chrom (length $length) ====="

    # Phase 1: scrape (resumable; scrape_aou.py skips a completed chromosome).
    "$PYTHON" "$SCRIPT_DIR/scrape_aou.py" \
        --chrom "$chrom" --length "$length" --outdir "$jdir" \
        --window "$WINDOW" --rowcount "$ROWCOUNT"

    # Phase 2: combine JSON -> CSV (rebuild if the scrape finished more recently).
    mkdir -p "$OUTDIR/csv"
    if [ ! -s "$csv" ] || [ "$FORCE" = "1" ] || { [ -f "$complete" ] && [ "$complete" -nt "$csv" ]; }; then
        "$PYTHON" "$SCRIPT_DIR/json_to_csv.py" --json-dir "$jdir" --out-csv "$csv"
    else
        echo "[skip] CSV up to date: $csv"
    fi

    # Phase 3: CSV -> VCF -> sort -> index.
    mkdir -p "$OUTDIR/full"
    if stale "$full" "$csv" || [ ! -f "$full.tbi" ]; then
        local tmp="$OUTDIR/full/$chrom.unsorted.vcf.gz"
        "$PYTHON" "$SCRIPT_DIR/csv_to_vcf.py" --csv "$csv" --out-vcf "$tmp" --reference "$REFERENCE"
        "$BCFTOOLS" sort "$tmp" -Oz -o "$full"
        "$TABIX" -p vcf "$full"
        rm -f "$tmp"
        echo "[ok] full VCF: $full"
    else
        echo "[skip] full VCF up to date: $full"
    fi

    # Phase 4: MAF filter -> more_than_001/
    mkdir -p "$OUTDIR/$FILTER_SUBDIR"
    if stale "$filt" "$full" || [ ! -f "$filt.tbi" ]; then
        "$BCFTOOLS" view -i "INFO/AF > $MAF" "$full" -Oz -o "$filt"
        "$TABIX" -p vcf "$filt"
        echo "[ok] filtered (INFO/AF > $MAF): $filt"
    else
        echo "[skip] filtered VCF up to date: $filt"
    fi
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
case "${1:-}" in
    -h|--help)
        sed -n '2,67p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

if [ "$#" -gt 0 ]; then
    CHROMS="$*"
else
    CHROMS="$ALL_CHROMS"
fi

check_prereqs || exit 1
mkdir -p "$OUTDIR"

n_ok=0
n_fail=0
failed=""
for chrom in $CHROMS; do
    length="$(chrom_length "$chrom" || true)"
    if [ -z "$length" ]; then
        echo "[SKIP] $chrom — not a GRCh38 primary chromosome (chr1-22, chrX, chrY)" >&2
        n_fail=$((n_fail + 1)); failed="$failed $chrom"; continue
    fi
    if build_chrom "$chrom" "$length"; then
        n_ok=$((n_ok + 1))
    else
        echo "[FAIL] $chrom" >&2
        n_fail=$((n_fail + 1)); failed="$failed $chrom"
    fi
done

echo
echo "Summary: $n_ok built, $n_fail failed."
if [ "$n_fail" -gt 0 ]; then
    echo "Failed:$failed" >&2
    exit 1
fi
echo "All requested AoU chromosomes built and MAF-filtered (INFO/AF > $MAF) in $OUTDIR."

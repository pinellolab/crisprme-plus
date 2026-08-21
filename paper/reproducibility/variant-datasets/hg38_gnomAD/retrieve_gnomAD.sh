#!/usr/bin/env bash
#
# retrieve_gnomAD.sh
# ==================
# Retrieve and prepare the gnomAD v4.1 "joint" (exome + genome) SITES-ONLY data
# for CRISPRme+. The script runs up to three phases:
#
#   A. DOWNLOAD  the raw gnomAD v4.1 joint sites VCFs (GRCh38), MD5-verified.
#   B. CONVERT   them to a CRISPRme-compatible format with the built-in
#                gnomAD-converter (raw gnomAD sites-only VCFs are not directly
#                usable by CRISPRme).
#   C. FILTER    the converted variants by allele frequency, keeping only those
#                with MAF > 0.001 (the "more_than_001" subset).
#
# Phases B and C run automatically after a successful download when the gnomAD
# samplesID file is present — by default the samplesIDs.gnomad.v41.txt bundled
# next to this script. If that file is absent (and SAMPLESID is not pointed at
# another path), the script stops after phase A, so a plain download still works.
#
# Dataset .......... Genome Aggregation Database (gnomAD), v4.1 joint sites-only
# Release .......... v4.1 (April 2024; files dated 2024-05-01)
# Reference build .. GRCh38 / hg38
# Files ............ gnomad.joint.v4.1.sites.<chrom>.vcf.bgz
# Chromosomes ...... chr1-22, chrX, chrY  (24 files)
# Source ........... gnomAD public bucket on Google Cloud Storage
#                    https://storage.googleapis.com/gcp-public-data--gnomad/
#                      release/4.1/vcf/joint/
#
# gnomAD is an AGGREGATE, SITES-ONLY resource: the VCFs carry population allele
# counts/frequencies in INFO and contain NO per-sample genotypes. The download
# phase fetches the RAW, unfiltered v4.1 joint sites files; the convert phase adds
# the per-ancestry pseudo-sample columns CRISPRme expects (reading INFO/AF_joint).
#
# The per-file MD5 checksums embedded below are the objects' own MD5s as reported
# by Google Cloud Storage (the `x-goog-hash: md5=...` metadata, base64-decoded to
# hex), so the download is self-contained and independently verifiable.
#
# SIZE WARNING: these are very large. The 24 files total ~817 GB (chr1 alone is
# ~72 GB; the smallest, chrY, is ~0.8 GB). Ensure ample free disk, and run inside
# tmux/screen.
#
# Outputs
# -------
#   $OUTDIR/gnomad.joint.v4.1.sites.<chrom>.vcf.bgz            (A: raw download)
#   $OUTDIR/gnomad.joint.v4.1.sites.<chrom>.biallelic.vcf.gz  (B: converted)
#   $OUTDIR/more_than_001/…biallelic.vcf.gz (+ .tbi)          (C: MAF-filtered)
#
# Download guarantees: integrity (per-file MD5 vs the embedded manifest;
# already-valid files are skipped), resumable (wget -c / curl -C -, never
# restarted from byte 0), and retrying (bounded per-attempt timeout, capped
# attempts). The raw .tbi indexes are not downloaded (regenerate with
# `tabix -p vcf <file>` if needed); the MAF-filtered outputs are indexed here.
#
# Dependencies
# ------------
#   * wget (preferred) OR curl — resumable, retrying HTTPS downloads (phase A).
#   * one of: md5sum | md5 | openssl — MD5 verification (phase A).
#   * the CRISPRme environment: crisprme.py + pysam + bcftools — convert (phase B).
#   * bcftools — the MAF filter (phase C; also used internally by the converter).
#
# Usage
# -----
#   ./retrieve_gnomAD.sh [chrom ...]
#
#   Full pipeline — download + convert + MAF-filter (uses the bundled samplesID):
#       ./retrieve_gnomAD.sh
#   A subset of chromosomes (strongly recommended given the size):
#       ./retrieve_gnomAD.sh chr22 chrY
#
# Configuration (environment variables, all optional)
# ---------------------------------------------------
#   OUTDIR         destination directory (default: ./hg38_gnomAD).
#   BASE_URL       override the download base URL.
#   DL_TRIES       max attempts per file           (default: 10).
#   DL_TIMEOUT     per-attempt connect/read timeout, seconds (default: 60).
#   DL_WAITRETRY   base wait between retries, seconds (default: 10).
#   SAMPLESID      path to the gnomAD samplesID file (default: the bundled
#                  samplesIDs.gnomad.v41.txt next to this script). Phases B and C
#                  run when this file exists; if it is missing, they are skipped.
#   CRISPRME       CRISPRme entry point for the converter (default: crisprme.py).
#   BCFTOOLS       bcftools command (default: bcftools).
#   THREADS        threads for the gnomAD-converter (default: 4).
#   MAF            allele-frequency threshold; keep INFO/AF > MAF (default: 0.001).
#   FILTER_SUBDIR  subfolder for the MAF-filtered VCFs (default: more_than_001).
#
# Exit status: 0 on success; non-zero if any download failed verification or a
# processing phase failed.
#
set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
BASE_URL="${BASE_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint}"
OUTDIR="${OUTDIR:-hg38_gnomAD}"
FILE_TEMPLATE="gnomad.joint.v4.1.sites.{C}.vcf.bgz"

DL_TRIES="${DL_TRIES:-10}"
DL_TIMEOUT="${DL_TIMEOUT:-60}"
DL_WAITRETRY="${DL_WAITRETRY:-10}"

# Processing steps (phases B and C) — run after download when the samplesID exists
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # this script's directory
SAMPLESID="${SAMPLESID:-$SCRIPT_DIR/samplesIDs.gnomad.v41.txt}"  # default: bundled samplesID
CRISPRME="${CRISPRME:-crisprme.py}"              # CRISPRme entry point (gnomAD-converter)
BCFTOOLS="${BCFTOOLS:-bcftools}"                 # bcftools (MAF filter; also used by converter)
THREADS="${THREADS:-4}"                          # threads for the converter
MAF="${MAF:-0.001}"                              # keep INFO/AF > MAF
FILTER_SUBDIR="${FILTER_SUBDIR:-more_than_001}"  # output subfolder for filtered VCFs

ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# ------------------------------------------------------------------------------
# MD5 manifest (from the GCS object metadata: x-goog-hash md5, base64 -> hex)
# Keyed by chromosome; the expected digest is for the .vcf.bgz of that chrom.
# ------------------------------------------------------------------------------
expected_md5_for_chrom() {
    case "$1" in
        chr1)  echo "11c62331b0a654fce6a9cd43838de648" ;;
        chr2)  echo "563a8fe6f148621169b0215ac9f19602" ;;
        chr3)  echo "c44f1661bafc15685f1eee593b4886ea" ;;
        chr4)  echo "e6d438b84539c6adad5bc67d6febb33b" ;;
        chr5)  echo "d226db0e0055b5e87e72f4b63158664a" ;;
        chr6)  echo "2820f13a2439ebbdf55066a0c320cdb5" ;;
        chr7)  echo "d1d50e4fa082246a5787eee76c036189" ;;
        chr8)  echo "195f2825e94c9b5e43b34bb2b1ab5c7b" ;;
        chr9)  echo "1c739fb01fd9de816e3fddc958668627" ;;
        chr10) echo "e2f174f150b5d709d5d7349ac241c438" ;;
        chr11) echo "b8651f2e5a0aafa23d7fc3406b35bc69" ;;
        chr12) echo "62795bafd326eae566ef49781a26bc91" ;;
        chr13) echo "92244327bee6d45973f6077a8134ccd9" ;;
        chr14) echo "f7a8344b03a4162cb71cdb628dd1e15b" ;;
        chr15) echo "40c8ab829f973688d2ef891ce5acabb7" ;;
        chr16) echo "58a5f920fc191b2069126c41278cc077" ;;
        chr17) echo "aa48657f45d8db7711c69fbf71a25cdc" ;;
        chr18) echo "80dd729bc61be464c964d3a3bfb0f41a" ;;
        chr19) echo "1853ca4993ceb25bd6f3a4554173f7cf" ;;
        chr20) echo "09263d3c29b822760c61607c6398f5c4" ;;
        chr21) echo "2ec2d9876d61fc9c5b1e84ab6841b62e" ;;
        chr22) echo "df15a5ea8ae2e3090eae112f548c74ef" ;;
        chrX)  echo "a5288ced0c2fe893fcfae4d2022b9cd9" ;;
        chrY)  echo "7b882f00919d582139acbc116a7a559f" ;;
        *)     return 1 ;;
    esac
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Print the MD5 hex digest of a file, using whatever tool is available.
# Portable across Linux (md5sum), macOS/BSD (md5), and openssl
md5_of() {
    local file="$1"
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$file" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then
        md5 -q "$file"
    elif command -v openssl >/dev/null 2>&1; then
        openssl md5 "$file" | awk '{print $NF}'
    else
        echo "ERROR: no MD5 tool found (need md5sum, md5, or openssl)" >&2
        return 2
    fi
}

# Return 0 if $file exists and its MD5 equals $expected, non-zero otherwise
verify_md5() {
    local file="$1" expected="$2" actual
    [ -f "$file" ] || return 1
    actual="$(md5_of "$file")" || return 2
    [ "$actual" = "$expected" ]
}

# Resumable, retrying download of $url into directory $dir (keeping the remote
# filename). Prefers wget; falls back to curl. Never restarts a partial file
download_resumable() {
    local url="$1" dir="$2" dest="$3"
    if command -v wget >/dev/null 2>&1; then
        # -c continue partial; --tries cap attempts; --timeout bounds each
        # stalled connect/read; --waitretry backs off between attempts
        wget -c \
             --tries="$DL_TRIES" \
             --timeout="$DL_TIMEOUT" \
             --waitretry="$DL_WAITRETRY" \
             --retry-connrefused \
             -P "$dir" "$url"
    elif command -v curl >/dev/null 2>&1; then
        # -C - resume; --retry cap attempts; --connect-timeout bounds connect
        curl -fL \
             -C - \
             --retry "$DL_TRIES" \
             --retry-delay "$DL_WAITRETRY" \
             --connect-timeout "$DL_TIMEOUT" \
             -o "$dest" "$url"
    else
        echo "ERROR: neither wget nor curl is available" >&2
        return 3
    fi
}

# ------------------------------------------------------------------------------
# Phase A — download the raw v4.1 joint sites VCFs (MD5-verified, resumable)
# Returns non-zero if any requested file is missing or fails verification.
# ------------------------------------------------------------------------------
phase_download() {
    echo "===== PHASE A: download raw gnomAD v4.1 joint sites ====="
    mkdir -p "$OUTDIR"
    echo "Destination : $OUTDIR"
    echo "Source      : $BASE_URL"
    echo

    local n_ok=0 n_skip=0 n_fail=0 failed_files="" chrom expected fname dest url
    for chrom in $CHROMS; do
        expected="$(expected_md5_for_chrom "$chrom" || true)"
        if [ -z "$expected" ]; then
            echo "[SKIP] $chrom — not part of this release (chr1-22, chrX, chrY only)"
            n_fail=$((n_fail + 1)); failed_files="$failed_files $chrom"; continue
        fi
        fname="${FILE_TEMPLATE/\{C\}/$chrom}"
        dest="$OUTDIR/$fname"
        url="$BASE_URL/$fname"

        if verify_md5 "$dest" "$expected"; then
            echo "[SKIP] $fname — already present, MD5 OK"
            n_skip=$((n_skip + 1)); continue
        fi
        echo "[GET ] $fname"
        if ! download_resumable "$url" "$OUTDIR" "$dest"; then
            echo "[FAIL] $fname — download did not complete after $DL_TRIES attempts" >&2
            n_fail=$((n_fail + 1)); failed_files="$failed_files $fname"; continue
        fi
        if verify_md5 "$dest" "$expected"; then
            echo "[ OK ] $fname — MD5 verified ($expected)"
            n_ok=$((n_ok + 1))
        else
            echo "[FAIL] $fname — MD5 mismatch (expected $expected, got $(md5_of "$dest" 2>/dev/null || echo '?'))" >&2
            echo "       The local copy is corrupt or a stale partial. Delete it and re-run:" >&2
            echo "         rm -f '$dest'" >&2
            n_fail=$((n_fail + 1)); failed_files="$failed_files $fname"
        fi
    done

    echo
    echo "Phase A summary: $n_ok downloaded, $n_skip already-valid, $n_fail failed."
    if [ "$n_fail" -gt 0 ]; then
        echo "Failed:$failed_files" >&2
        return 1
    fi
    echo "All requested raw gnomAD v4.1 joint sites VCFs are present and MD5-verified."
}

# ------------------------------------------------------------------------------
# Phase B — convert to CRISPRme format with the built-in gnomAD-converter.
# Produces $OUTDIR/gnomad.joint.v4.1.sites.<chrom>.biallelic.vcf.gz (biallelic,
# PASS-only by default), reading INFO/AF_joint (the --joint mode).
# ------------------------------------------------------------------------------
phase_convert() {
    echo "===== PHASE B: convert to CRISPRme format (gnomAD-converter) ====="
    command -v "$CRISPRME" >/dev/null 2>&1 || {
        echo "ERROR: '$CRISPRME' not found — needed to convert. Activate the CRISPRme" >&2
        echo "       environment or set CRISPRME=/path/to/crisprme.py." >&2
        return 1
    }
    [ -f "$SAMPLESID" ] || {
        echo "ERROR: samplesID file not found: '${SAMPLESID:-<unset>}'" >&2
        echo "       Set SAMPLESID to CRISPRme's samplesIDs.gnomad.v41.txt." >&2
        return 1
    }
    # v4.1 files carry INFO/AF_joint -> use --joint. Defaults keep biallelic
    # output and discard non-PASS variants.
    "$CRISPRME" gnomAD-converter \
        --gnomAD_VCFdir "$OUTDIR" \
        --samplesID "$SAMPLESID" \
        --joint \
        --threads "$THREADS"
    echo "Phase B done: converted VCFs written as *.biallelic.vcf.gz in $OUTDIR"
}

# ------------------------------------------------------------------------------
# Phase C — MAF filter: keep converted variants with INFO/AF > $MAF.
# Writes $OUTDIR/$FILTER_SUBDIR/<same name> (bgzipped + tabix-indexed).
# ------------------------------------------------------------------------------
phase_filter() {
    echo "===== PHASE C: MAF filter (keep INFO/AF > $MAF) ====="
    command -v "$BCFTOOLS" >/dev/null 2>&1 || {
        echo "ERROR: '$BCFTOOLS' not found — needed for the MAF filter." >&2
        return 1
    }
    local fdir="$OUTDIR/$FILTER_SUBDIR"
    mkdir -p "$fdir"
    local converted base out n_done=0 n_skip=0
    for converted in "$OUTDIR"/gnomad.joint.v4.1.sites.*.biallelic.vcf.gz; do
        [ -e "$converted" ] || continue  # no matches -> glob stays literal
        base="$(basename "$converted")"
        out="$fdir/$base"
        if [ -f "$out" ]; then
            echo "[SKIP] $base — already filtered"
            n_skip=$((n_skip + 1)); continue
        fi
        echo "[FILT] $base"
        "$BCFTOOLS" view -i "INFO/AF > $MAF" -Oz -o "$out" "$converted"
        "$BCFTOOLS" index -t "$out"
        n_done=$((n_done + 1))
    done
    if [ "$n_done" -eq 0 ] && [ "$n_skip" -eq 0 ]; then
        echo "ERROR: no converted *.biallelic.vcf.gz found in $OUTDIR — run phase B first." >&2
        return 1
    fi
    echo "Phase C done: $n_done filtered, $n_skip already-filtered -> $fdir (INFO/AF > $MAF)"
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
case "${1:-}" in
    -h|--help)
        sed -n '2,87p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

# Chromosomes: command-line arguments if given, else the full release set
if [ "$#" -gt 0 ]; then
    CHROMS="$*"
else
    CHROMS="$ALL_CHROMS"
fi

phase_download || exit 1

if [ -f "$SAMPLESID" ]; then
    echo "samplesID    : $SAMPLESID"
    phase_convert
    phase_filter
    echo
    echo "Done: downloaded, converted, and MAF-filtered (INFO/AF > $MAF) in $OUTDIR."
else
    echo
    echo "Note: samplesID file not found ($SAMPLESID) — stopping after download (phase A)."
    echo "      Provide it (or set SAMPLESID=/path/to/samplesIDs.gnomad.v41.txt) to"
    echo "      also convert and MAF-filter."
fi

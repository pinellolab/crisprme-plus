#!/usr/bin/env bash
#
# retrieve_1000G_2019.sh
# ======================
# Download the 1000 Genomes Project (1KGP) Phase 3 GRCh38 2019 reanalysis
# the March 2019 *biallelic SNV + INDEL* phased call set for use with
# CRISPRme+.
#
# Dataset .......... 1000 Genomes Project, GRCh38 2019 biallelic reanalysis
# Release ID ....... 20190312_biallelic_SNV_and_INDEL
# Reference build .. GRCh38 / hg38
# Files ............ ALL.<chrom>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
# Chromosomes ...... chr1-22, chrX  (23 files; no chrY / chrM in this release)
# Source ........... International Genome Sample Resource (IGSR) @ EBI
#                    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
#                      1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
#                    (an ftp:// mirror of the same tree also exists on that host)
#
# This is the exact release CRISPRme's built-in `crisprme.py setup` downloads
# for its default 1000G data. The per-file MD5 checksums embedded below are the
# same known-good digests shipped in the tool (PostProcess/utils.py: MD51000G),
# reproduced here so this script is self-contained and independently verifiable.
#
# What the script guarantees
# --------------------------
#   * Integrity   — every downloaded .vcf.gz is MD5-verified against the embedded
#                   manifest; a file that already exists and matches is skipped.
#   * Resumable   — an interrupted transfer is CONTINUED from where it stopped
#                   (wget -c / curl -C -), never restarted from byte 0.
#   * Retrying    — each transfer retries on transient network errors, with a
#                   bounded per-attempt timeout and a capped number of attempts.
#
# What the script does NOT do (by design)
# ----------------------------------------
#   * It does not download the .tbi tabix indexes: the release ships them, but
#     they carry no published MD5 to verify, and CRISPRme (via tabix) rebuilds
#     them deterministically from the .vcf.gz. Index locally afterwards if you
#     want them, e.g.:  for f in "$OUTDIR"/*.vcf.gz; do tabix -p vcf "$f"; done
#   * It does not build the CRISPRme search index (that is a separate step).
#
# Dependencies
# ------------
#   * wget (preferred) OR curl — for resumable, retrying HTTPS downloads.
#   * one of: md5sum (Linux) | md5 (macOS/BSD) | openssl — for MD5 verification.
#
# Usage
# -----
#   ./retrieve_1000G_2019.sh [chrom ...]
#
#   With no arguments all 23 chromosomes are fetched. Pass one or more
#   chromosome labels to fetch only those (handy for a quick test), e.g.:
#       ./retrieve_1000G_2019.sh chr22
#       ./retrieve_1000G_2019.sh chr21 chr22 chrX
#
# Configuration (environment variables, all optional)
# ---------------------------------------------------
#   OUTDIR         destination directory (default: ./hg38_1000G_2019).
#                  For a CRISPRme run, point it at <workdir>/VCFs/hg38_1000G_2019.
#   BASE_URL       override the source base URL (e.g. to use an ftp:// mirror).
#   DL_TRIES       max attempts per file           (default: 10).
#   DL_TIMEOUT     per-attempt connect/read timeout, seconds (default: 60).
#   DL_WAITRETRY   base wait between retries, seconds (default: 10).
#
# Exit status: 0 if every requested file is present and MD5-verified, non-zero
# if any file failed to download or failed verification.
#
set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
BASE_URL="${BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL}"
OUTDIR="${OUTDIR:-hg38_1000G_2019}"
FILE_TEMPLATE="ALL.{C}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

DL_TRIES="${DL_TRIES:-10}"
DL_TIMEOUT="${DL_TIMEOUT:-60}"
DL_WAITRETRY="${DL_WAITRETRY:-10}"

ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

# ------------------------------------------------------------------------------
# MD5 manifest (mirrors PostProcess/utils.py: MD51000G)
# Keyed by chromosome; the expected digest is for the .vcf.gz of that chrom.
# ------------------------------------------------------------------------------
expected_md5_for_chrom() {
    case "$1" in
        chr1)  echo "77f154e53c2b7c36b04d03bab3af8b74" ;;
        chr2)  echo "f9d29c4935e591b2b269eed7cd7e35d8" ;;
        chr3)  echo "6e59d00235de71562b4199e09b7e5934" ;;
        chr4)  echo "70a2c1ede97eceb7baeea06c8e46cf3c" ;;
        chr5)  echo "74d5486c0fd29b0e6add24d3740fc3b4" ;;
        chr6)  echo "8c5d83c1a9253058120368af39baf0c8" ;;
        chr7)  echo "dfaa282712fc1292146173dd2ffeb1d9" ;;
        chr8)  echo "ddf7b370fcee63462037c237f12b4444" ;;
        chr9)  echo "5ade69521dc50d88ad7c91bf4ec6fcd8" ;;
        chr10) echo "1c409a674426eda2fd29b49078137c5d" ;;
        chr11) echo "65339bffc61bc97f2130832fe9f84d7c" ;;
        chr12) echo "9a1bda389121140d30c768ef6a1b1370" ;;
        chr13) echo "47b0463541be137a8bbfe40f6aade864" ;;
        chr14) echo "241aedf0792c45d5345d421105c782af" ;;
        chr15) echo "b48e7c64e35b727d34786faa76467f94" ;;
        chr16) echo "1ce7d66799cab6718852d78dd2aab765" ;;
        chr17) echo "ecc22783fd1ee7a1c66b053491873192" ;;
        chr18) echo "fdf3e460e91cd955a9e8cebf01b5d815" ;;
        chr19) echo "a2f17e4ec552fc07cbd05c1eac0cf7ec" ;;
        chr20) echo "155c3b440d7990630132e4756f7fcc85" ;;
        chr21) echo "52882490028507e5d4e606b0905072b1" ;;
        chr22) echo "57a1722e6ed7d9df08cb3c0e42b62d53" ;;
        chrX)  echo "e6a3d41811faee60de177061edcd6fe6" ;;
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
# Main
# ------------------------------------------------------------------------------
case "${1:-}" in
    -h|--help)
        sed -n '2,80p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

# Chromosomes: command-line arguments if given, else the full release set
if [ "$#" -gt 0 ]; then
    CHROMS="$*"
else
    CHROMS="$ALL_CHROMS"
fi

mkdir -p "$OUTDIR"
echo "Destination : $OUTDIR"
echo "Source      : $BASE_URL"
echo

n_ok=0
n_skip=0
n_fail=0
failed_files=""

for chrom in $CHROMS; do
    expected="$(expected_md5_for_chrom "$chrom" || true)"
    if [ -z "$expected" ]; then
        echo "[SKIP] $chrom — not part of this release (chr1-22, chrX only)"
        n_fail=$((n_fail + 1))
        failed_files="$failed_files $chrom"
        continue
    fi

    fname="${FILE_TEMPLATE/\{C\}/$chrom}"
    dest="$OUTDIR/$fname"
    url="$BASE_URL/$fname"

    # Already present and correct -> nothing to do
    if verify_md5 "$dest" "$expected"; then
        echo "[SKIP] $fname — already present, MD5 OK"
        n_skip=$((n_skip + 1))
        continue
    fi

    # Download (resumes a partial file; retries transient failures)
    echo "[GET ] $fname"
    if ! download_resumable "$url" "$OUTDIR" "$dest"; then
        echo "[FAIL] $fname — download did not complete after $DL_TRIES attempts" >&2
        n_fail=$((n_fail + 1))
        failed_files="$failed_files $fname"
        continue
    fi

    # Verify integrity of the freshly downloaded file
    if verify_md5 "$dest" "$expected"; then
        echo "[ OK ] $fname — MD5 verified ($expected)"
        n_ok=$((n_ok + 1))
    else
        echo "[FAIL] $fname — MD5 mismatch (expected $expected, got $(md5_of "$dest" 2>/dev/null || echo '?'))" >&2
        echo "       The local copy is corrupt or a stale partial. Delete it and re-run:" >&2
        echo "         rm -f '$dest'" >&2
        n_fail=$((n_fail + 1))
        failed_files="$failed_files $fname"
    fi
done

echo
echo "Summary: $n_ok downloaded, $n_skip already-valid, $n_fail failed."
if [ "$n_fail" -gt 0 ]; then
    echo "Failed:$failed_files" >&2
    exit 1
fi
echo "All requested 1000G 2019 VCFs are present and MD5-verified in $OUTDIR."

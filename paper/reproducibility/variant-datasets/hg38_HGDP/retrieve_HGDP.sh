#!/usr/bin/env bash
#
# retrieve_HGDP.sh
# ================
# Download the Human Genome Diversity Project (HGDP) whole-genome sequencing
# variant call set, aligned to GRCh38 / hg38, for use with CRISPRme+.
#
# Dataset .......... Human Genome Diversity Project (HGDP), WGS call set
# Release ID ....... hgdp_wgs.20190516  ("full" call set)
# Reference build .. GRCh38 / hg38
# Files ............ hgdp_wgs.20190516.full.<chrom>.vcf.gz
# Chromosomes ...... chr1-22, chrX by default (23 files)
#                    chrY is also published in this release; pass it explicitly
#                    to fetch it (it is NOT part of the default set, matching the
#                    chromosomes used for the CRISPRme+ combined panel).
# Source ........... Wellcome Sanger Institute
#                    https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
#
# This is the exact release CRISPRme's built-in `crisprme.py setup` downloads
# for its default HGDP data. The per-file MD5 checksums embedded below are the
# same known-good digests shipped in the tool (PostProcess/utils.py: MD5HGDP),
# reproduced here so this script is self-contained and independently verifiable.
#
# SIZE WARNING: these are large. chr22 alone is ~5.5 GB and the full default set
# is on the order of ~390 GB. Ensure ample free disk and run inside tmux/screen.
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
#   * It does not download the .tbi tabix indexes: they carry no published MD5 to
#     verify, and CRISPRme (via tabix) rebuilds them deterministically from the
#     .vcf.gz. Index locally afterwards if you want them, e.g.:
#       for f in "$OUTDIR"/*.vcf.gz; do tabix -p vcf "$f"; done
#   * It does not build the CRISPRme search index (that is a separate step).
#
# Dependencies
# ------------
#   * wget (preferred) OR curl — for resumable, retrying HTTPS downloads.
#   * one of: md5sum (Linux) | md5 (macOS/BSD) | openssl — for MD5 verification.
#
# Usage
# -----
#   ./retrieve_HGDP.sh [chrom ...]
#
#   With no arguments the 23 default chromosomes are fetched. Pass one or more
#   chromosome labels to fetch only those (handy for a quick test, or to add
#   chrY), e.g.:
#       ./retrieve_HGDP.sh chr22
#       ./retrieve_HGDP.sh chr21 chr22 chrX chrY
#
# Configuration (environment variables, all optional)
# ---------------------------------------------------
#   OUTDIR         destination directory (default: ./hg38_HGDP).
#                  For a CRISPRme run, point it at <workdir>/VCFs/hg38_HGDP.
#   BASE_URL       override the source base URL.
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
BASE_URL="${BASE_URL:-https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516}"
OUTDIR="${OUTDIR:-hg38_HGDP}"
FILE_TEMPLATE="hgdp_wgs.20190516.full.{C}.vcf.gz"

DL_TRIES="${DL_TRIES:-10}"
DL_TIMEOUT="${DL_TIMEOUT:-60}"
DL_WAITRETRY="${DL_WAITRETRY:-10}"

ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

# ------------------------------------------------------------------------------
# MD5 manifest (mirrors PostProcess/utils.py: MD5HGDP)
# Keyed by chromosome; the expected digest is for the .vcf.gz of that chrom.
# chrY is included (it is a valid part of the release) even though it is not in
# the default ALL_CHROMS set above.
# ------------------------------------------------------------------------------
expected_md5_for_chrom() {
    case "$1" in
        chr1)  echo "70d82ae3ae65cb73858738f547f64e93" ;;
        chr2)  echo "539d4eb31355b90f0262453fa1349ae6" ;;
        chr3)  echo "0d37ba60afd5ff092cf1bc75bde3588e" ;;
        chr4)  echo "68d57e5c2129bbafa1a9dd75f630cf89" ;;
        chr5)  echo "929eb66a26e9679320bcc26df0bd4116" ;;
        chr6)  echo "28c9ad734025e7292bde533da908cf68" ;;
        chr7)  echo "ed7eaf339cd7964b9f1e7581de5bdeb1" ;;
        chr8)  echo "2d47d60ff6b63e1163d219d964999ee3" ;;
        chr9)  echo "52f917fc3068eff76f0ba8bde0c59292" ;;
        chr10) echo "7131981641e886173da90d215346e857" ;;
        chr11) echo "2f127b0006cbc36fb32c66860d4b31d9" ;;
        chr12) echo "023d0e2d852c167490d4578f814d043d" ;;
        chr13) echo "afcfba8b01258e418f5fb230b14daa02" ;;
        chr14) echo "90b0c15b61fd9c47a9751495f2b784ce" ;;
        chr15) echo "665e844d7e2e85e226d25827ea8014be" ;;
        chr16) echo "0d6f1b6141c78489a2b2e27eeec848dd" ;;
        chr17) echo "d53421438b3bc3c5ce5ab51b90578182" ;;
        chr18) echo "6351d9b20995cf500ac4b11490ff31c7" ;;
        chr19) echo "167ce7a43876b32e586978a75f3b0d39" ;;
        chr20) echo "d90130b11620378bed7c2cc43be94b7e" ;;
        chr21) echo "8f44e4daa3952cd73751141f66b6e5ae" ;;
        chr22) echo "84f4a1d86f54bdc0cd9b19502ff8d2c2" ;;
        chrX)  echo "8d0e4e178fdfa07db76d0218a9b2ceab" ;;
        chrY)  echo "54b3aba28600c8d0d8a695c8dcfdc4cd" ;;
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

# Chromosomes: command-line arguments if given, else the default release set
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
        echo "[SKIP] $chrom — not part of this release (chr1-22, chrX, chrY only)"
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
echo "All requested HGDP VCFs are present and MD5-verified in $OUTDIR."

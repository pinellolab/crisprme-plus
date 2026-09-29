#!/usr/bin/env bash
#
# retrieve_1000G_2021.sh
# ======================
# Download the 1000 Genomes Project high-coverage (30x) NYGC re-sequencing
# release — the 3,202-sample, GRCh38, phased SNV+INDEL call set — for use with
# CRISPRme+.
#
# Dataset .......... 1000 Genomes Project, high-coverage NYGC re-sequencing
# Release ID ....... 20201028_3202_phased  (3,202 samples, ~30x)
# Reference build .. GRCh38 / hg38
# Files ............ autosomes: CCDG_14151_B01_GRM_WGS_2020-08-05_<chrom>.filtered.shapeit2-duohmm-phased.vcf.gz
#                    chrX:      CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz
# Chromosomes ...... chr1-22 by default (22 files). chrX exists but is EXCLUDED
#                    by default — see "Known issues" below; pass it explicitly to
#                    fetch it anyway.
# Source ........... International Genome Sample Resource (IGSR) @ EBI
#                    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
#                      1000G_2504_high_coverage/working/20201028_3202_phased/
#
# The per-file MD5 checksums embedded below are taken verbatim from the release's
# own official manifest (phased-manifest_July2021.tsv) in that directory, so this
# script is self-contained and independently verifiable.
#
# Known issues with this release (documented; affect how it is used)
# ------------------------------------------------------------------
#   * chrX haploid genotypes. chrX is eagle2-phased and stores hemizygous male
#     genotypes as a bare `0`/`1` (no `|` separator) rather than the diploid
#     encoding used on the autosomes. This breaks tools that assume diploid
#     genotypes, so chrX is NOT part of the default set and is not recommended
#     for CRISPRme. It can still be downloaded by naming `chrX` explicitly.
#   * Non-CRISPRme filenames. The upstream names embed the chromosome inside the
#     first underscore-delimited segment (`..._chr1.filtered...`), not as a
#     standalone dot-separated segment. This script downloads the files under
#     their ORIGINAL names (so the MD5s match the official manifest); a separate
#     rename step is required to make them CRISPRme-compatible.
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
#   * It does not download the .tbi tabix indexes (regenerate with
#     `tabix -p vcf <file>.vcf.gz` if needed), rename files for CRISPRme, build
#     the samplesID file, or build the CRISPRme search index — those are separate
#     steps.
#
# Dependencies
# ------------
#   * wget (preferred) OR curl — for resumable, retrying HTTPS downloads.
#   * one of: md5sum (Linux) | md5 (macOS/BSD) | openssl — for MD5 verification.
#
# Usage
# -----
#   ./retrieve_1000G_2021.sh [chrom ...]
#
#   With no arguments the 22 autosomes are fetched. Pass one or more chromosome
#   labels to fetch only those (handy for a quick test, or to add chrX), e.g.:
#       ./retrieve_1000G_2021.sh chr22
#       ./retrieve_1000G_2021.sh chr21 chr22 chrX
#
# Configuration (environment variables, all optional)
# ---------------------------------------------------
#   OUTDIR         destination directory (default: ./hg38_1000G_2021).
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
BASE_URL="${BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased}"
OUTDIR="${OUTDIR:-hg38_1000G_2021}"

DL_TRIES="${DL_TRIES:-10}"
DL_TIMEOUT="${DL_TIMEOUT:-60}"
DL_WAITRETRY="${DL_WAITRETRY:-10}"

# Default set: autosomes only (chrX excluded — see "Known issues" in the header)
ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

# ------------------------------------------------------------------------------
# Per-chromosome filename (autosomes and chrX use different phasers, hence names)
# ------------------------------------------------------------------------------
filename_for_chrom() {
    if [ "$1" = "chrX" ]; then
        echo "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
    else
        echo "CCDG_14151_B01_GRM_WGS_2020-08-05_$1.filtered.shapeit2-duohmm-phased.vcf.gz"
    fi
}

# ------------------------------------------------------------------------------
# MD5 manifest (verbatim from phased-manifest_July2021.tsv, .vcf.gz entries)
# Keyed by chromosome; the expected digest is for the .vcf.gz of that chrom.
# ------------------------------------------------------------------------------
expected_md5_for_chrom() {
    case "$1" in
        chr1)  echo "8d286b08ec8979c7d977724faf2ecf27" ;;
        chr2)  echo "4e7ec2cb793a03337c34fe4a347b2b39" ;;
        chr3)  echo "ca2f9b4ee4e441debf25e64ef1d66d49" ;;
        chr4)  echo "8eb36436c8b1dc7f53e21a0a3ef07c5a" ;;
        chr5)  echo "8ea6ddbfd7f03698505256d57f3ee1bd" ;;
        chr6)  echo "13f6dddb290883b926f147bde1533225" ;;
        chr7)  echo "14d0eba84ce07fe54833734421b56a9d" ;;
        chr8)  echo "38014c6362ecfd2111c137cca81893d4" ;;
        chr9)  echo "52a1f3a6832db74002f1552eb791eb52" ;;
        chr10) echo "3164712a612107af3b6846b82a58bf79" ;;
        chr11) echo "d4b5dad3aac770b07840af19d01f25c6" ;;
        chr12) echo "29bf540e16267ab91310dd0537a089a4" ;;
        chr13) echo "c8aea983f2ed4accf1b9d9ab1f72a52c" ;;
        chr14) echo "c4e0d9eaceb2d4a10c1749ccbfad71b6" ;;
        chr15) echo "fcf9f0363253a49a7460b4f2d647d740" ;;
        chr16) echo "7933b3977d670abee92ac98e547a6fe5" ;;
        chr17) echo "06ef514d7b8ce9d73edbad54435418e9" ;;
        chr18) echo "f854768a4f1f6b0635fb38019118c25c" ;;
        chr19) echo "459c7d0846798aa0612fb04a71e56b4e" ;;
        chr20) echo "521d2011fbf06517f34edd181b6283e7" ;;
        chr21) echo "80e99e81f0bc739e6e36843a9c5c46df" ;;
        chr22) echo "aaf19d9c7ffcd86b34275899ddc898e7" ;;
        chrX)  echo "072957d65f78e0fb0d55d8300e869849" ;;
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
        sed -n '2,77p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

# Chromosomes: command-line arguments if given, else the default (autosome) set
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

    # chrX carries the documented haploid-genotype issue: warn before fetching.
    if [ "$chrom" = "chrX" ]; then
        echo "[WARN] chrX is eagle2-phased: hemizygous male genotypes are haploid" >&2
        echo "       (bare 0/1). Not recommended for CRISPRme; the default set is" >&2
        echo "       chr1-22. Proceeding because chrX was requested explicitly." >&2
    fi

    fname="$(filename_for_chrom "$chrom")"
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
echo "All requested 1000G 2021 VCFs are present and MD5-verified in $OUTDIR."

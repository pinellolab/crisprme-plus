#!/usr/bin/env bash
#
# build_combined_panel.sh — end-to-end: merge all chromosomes into one panel,
# assemble its samplesID, then enrich hg38 and build the (pamless) index.
#
# Resumable: per-chromosome merges that already produced <panel>/merged.<chr>.vcf.gz.tbi
# are skipped. Aborts before enrichment if the merge is incomplete (never builds a
# variant-less index by accident).
#
# See merge_vcf_panels.sh for the merge/provenance logic and dependencies (notably
# CRISPRitz PR #36 for enrichment). Configure SOURCES in that script.
#
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

CWD=${1:?usage: build_combined_panel.sh <crisprme_working_dir> [panel_name] [pam] [bdna] [brna] [parallel]}
PANEL=${2:-hg38_1000G_HGDP}
PAM=${3:-20bp-NNN-NO-PAM}          # pamless NNN by default (one index serves all PAMs)
BDNA=${4:-2}; BRNA=${5:-2}
PAR=${6:-4}
CHROMS=${CHROMS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"}
NCHR=$(wc -w <<<"$CHROMS")
VCFDIR="$CWD/VCFs/$PANEL"
BCFTOOLS=${BCFTOOLS:-bcftools}     # honor an explicit bcftools (e.g. a conda env's)

echo "===== PHASE A: per-chromosome merge ($PAR-way parallel) [$(date +%T)] ====="
printf '%s\n' $CHROMS | xargs -P "$PAR" -I{} bash "$HERE/merge_vcf_panels.sh" "$VCFDIR" {}
NMERGED=$(ls "$VCFDIR"/merged.*.vcf.gz.tbi 2>/dev/null | wc -l)
echo "merged: $NMERGED / $NCHR"
[ "$NMERGED" -eq "$NCHR" ] || { echo "ABORT: merge incomplete — not enriching/indexing"; exit 1; }

echo "===== PHASE B: panel samplesID (from the merged VCF's ACTUAL samples) [$(date +%T)] ====="
# A source samplesID can OVER-LIST its VCF (e.g. the 1000G metadata lists ~3500
# samples but the phased VCF contains 2548), so a blind union produces a
# samplesID with samples that have no genotypes in the merged panel — wrong
# population denominators (phantom hom-ref individuals inflate AN). Instead take the
# ACTUAL samples in the merged VCF and attach each one's population metadata from the
# source samplesIDs.
SID="$CWD/samplesIDs/$PANEL.samplesID.txt"
DATA=${CRISPRME_SID_DATA:-/srv/local/crisprme/data/samplesIDs}
_anychr=$(ls "$VCFDIR"/merged.*.vcf.gz 2>/dev/null | head -1)
# The ONE authoritative set of genotyped samples: exactly what the merged VCF holds.
GENOTYPED="$CWD/samplesIDs/.$PANEL.genotyped.txt"
"$BCFTOOLS" query -l "$_anychr" | sort -u > "$GENOTYPED"
_NGENO=$(wc -l < "$GENOTYPED")

# combined samplesID (unchanged): merged-VCF samples with their source metadata.
head -1 "$DATA/hg38_1000G.samplesID.txt" > "$SID"
{ tail -n +2 "$DATA/hg38_1000G.samplesID.txt"; tail -n +2 "$DATA/hg38_HGDP.samplesID.txt"; } \
  | awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($1 in keep)' "$GENOTYPED" - >> "$SID"
echo "combined samplesID rows: $(($(wc -l < "$SID") - 1))  (from $_NGENO merged-VCF samples)"

# PER-DB VCF-FILTERED samplesID lists + a LISTING naming them, so tier emission
# (crisprme.py _build_db_to_samplesid -> tier0_compile.build_sample_meta) builds the
# panel over GENOTYPED samples ONLY, keeping each dataset's provenance SEPARATE (never
# conflate 1000G vs HGDP). Each per-db file is the source samplesID intersected with
# the SAME genotyped set used for the combined samplesID above -> AN = 2 * (VCF-
# filtered 1000G + HGDP), not the over-listed source count.
SIDDIR="$CWD/samplesIDs"
LISTING="$SIDDIR/$PANEL.dblist.txt"
: > "$LISTING"
for _src in hg38_1000G hg38_HGDP; do
  _srcfile="$DATA/$_src.samplesID.txt"
  [ -f "$_srcfile" ] || { echo "WARN: missing source samplesID $_srcfile — skipping"; continue; }
  # Clean, panel-scoped filename so _db_name_from_samplesid yields a clean dataset
  # label (e.g. hg38_1000G.<panel>.samplesID.txt -> "1000G.<panel>"); the two files
  # keep DISTINCT labels so 1000G / HGDP provenance is never conflated.
  _out="$SIDDIR/$_src.$PANEL.samplesID.txt"
  head -1 "$_srcfile" > "$_out"
  tail -n +2 "$_srcfile" \
    | awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($1 in keep)' "$GENOTYPED" - >> "$_out"
  _rows=$(($(wc -l < "$_out") - 1))
  _srcrows=$(($(wc -l < "$_srcfile") - 1))
  echo "  $_src: $_rows genotyped / $_srcrows listed"
  # Build-time WARNING: a per-db list that still OVER-LISTS the VCF would re-inflate
  # AN with phantom hom-ref samples. This can only happen if the source is somehow
  # NOT a superset of the genotyped set; the intersection above prevents it, but we
  # assert it loudly so a wiring regression is visible.
  if [ "$_rows" -gt "$_srcrows" ]; then
    echo "  WARN [$_src]: filtered rows ($_rows) exceed source rows ($_srcrows) — panel denominator may be inflated"
  fi
  # LISTING entries are bare filenames resolved under samplesIDs/ (the format
  # _build_db_to_samplesid expects), one per line, distinct dataset labels.
  basename "$_out" >> "$LISTING"
done
echo "per-db samplesID listing -> $LISTING:"; sed 's/^/    /' "$LISTING"

echo "===== PHASE C: enrich + build index [$(date +%T)] ====="
cd "$CWD"
# Pass the PER-DB VCF-FILTERED listing as --samplesID so build-index-only emits the
# dictless tiers (registry + genotypes) with the CORRECT genotyped-only panel
# denominator. The combined $SID above stays the search-side samplesID; the LISTING
# is the tier-emission input.
crisprme.py build-index-only --genome "$CWD/Genomes/hg38" --pam "$CWD/PAMs/$PAM.txt" \
  --bDNA "$BDNA" --bRNA "$BRNA" --vcf "$VCFDIR" --path "$CWD" --samplesID "$LISTING"
echo "build-index-only exit=$? [$(date +%T)]"
echo "indexes: $(ls "$CWD/genome_library/" | grep -iE "$PANEL" | tr '\n' ' ')"
echo "BUILD_COMBINED_DONE [$(date +%T)]"

#!/bin/bash
# 1000G-2021 single-source, PHASED variant-index build (as-run, resumable) — a clean
# "real-haplotypes" default index NRG_3_hg38+hg38_1000G2021.
#
# WHAT THIS IS: a standard CRISPRme variant-aware index built by ENRICHING hg38 with the
# 1000 Genomes 2021 high-coverage (30x) NYGC/IGSR callset (CCDG_14151_B01_GRM_WGS_2020-08-05,
# filtered, shapeit2-duohmm PHASED; chrX is eagle2-phased). 3,202 samples, GRCh38,
# chr-prefixed, chr1..22 + chrX (the callset has NO chrY -> chrY is reference-only in the
# index, as in every 1000G index). Because the genotypes are PHASED, the index supports
# CONFIRMED cis co-occurrence + named per-sample carriers (like 1000G-2021+HGDP), and unlike
# the sites-only mega every off-target here corresponds to an OBSERVED haplotype in a real
# panel. This is the "simple/clean default" index; the mega remains the worst-case sites-only
# fallback (a haplotype it lists may or may not exist in any individual).
#
# WHY SINGLE-SOURCE: the 2021 callset is already genotyped+phased, so there is no cross-source
# merge (skip merge_panels Mode 1/2) -- just per-chrom normalization + AF recompute + a
# VCF-derived samplesID, then the standard build-index-only. (Same shape as hprc_build.sh.)
#
# CORRECTNESS: the samplesID is DERIVED FROM THE VCF HEADER (exactly the 3,202 genotyped
# samples) with population metadata joined from the existing 3,501-row 1000G-2021 samplesID.
# Using the raw 3,501-row file would over-list the panel and DEFLATE every reported AF
# (the known samplesID-over-listing bug); AN must be 2x the genotyped panel = 6,404.
#
# OUTPUT: genome_library/NRG_3_hg38+hg38_1000G2021 (+ _INDELS) + Dictionaries/{registry,
# genotypes,indel_genotypes,log_indels}_hg38_1000G2021 + samplesIDs/hg38_1000G2021.samplesID.txt.
# Registry is RAW (CRISPRME_REGISTRY_COMPRESS=0), SNP+indel co-occurrence ON (CRISPRME_INDEL_SNP)
# -- consistent with the other production indexes. Built with the current dev code overlaid on
# the v2.5.5 SIF so the reg_<chrom>.idx carry the data_type/phased manifest (genotyped-phased).
#
# PATHS are as-run cluster paths; override via env. Requires the crisprme dev overlay ($OV)
# with the 4 scoring pkls copied in, and an apptainer SIF ($SIF).
set -u
SIF=${SIF:-/srv/local/lp698/crisprme_v255.sif}                 # apptainer image w/ crisprme env
OV=${OV:-/srv/local/lp698/dev_build_src}                       # dev crisprme overlay (data_type stamp)
B=${B:-/srv/local/lp698/onekg2021_build}                       # build working dir (--path)
SRC=${SRC:-/srv/local/lp698/mode1_2021_gw/DATA/SRC_1000G_2021} # CCDG_14151_..._<chr>.filtered.*phased*.vcf.gz
GENOME=${GENOME:-/srv/local/lp698/mode1_2021_gw/DATA/Genomes/hg38}  # per-chrom hg38 fastas (+ .fai)
PAM=${PAM:-/srv/local/lp698/mode1_2021_gw/DATA/PAMs/20bp-NRG-SpCas9.txt}  # NRG (NAG+NGG) SpCas9 PAM
SID_META=${SID_META:-/srv/local/lp698/mode1_2021_gw/DATA/samplesIDs/hg38_1000G2021.samplesID.txt}  # pop-metadata source (over-listed; used only for join)
THREADS=${THREADS:-64}
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

mkdir -p "$B/VCFs/hg38_1000G2021" "$B/samplesIDs" "$B/tmp"

# dev overlay: bind crisprme.py + the crisprme opt tree (data_type stamp, RAW default, INDEL_SNP).
# Bind a REAL on-disk /tmp ($B/tmp on the 6.7 TB volume) instead of --writable-tmpfs: the tiny
# in-memory tmpfs overflows when bcftools sort spills the large 3202-sample chromosomes
# ("bgzf_flush File write failed"). bcftools sort also gets an explicit -T on the same volume.
run() { apptainer exec -B /srv/local/lp698 -B /data/pinello -B "$B/tmp:/tmp" \
  -B "$OV/crisprme.py:/opt/conda/bin/crisprme.py" \
  -B "$OV:/opt/conda/opt/crisprme" \
  "$SIF" bash -c "export PATH=/opt/conda/bin:\$PATH TMPDIR=/tmp; $*"; }

# 1) NORMALIZE each chrom: split multiallelics + left-align (norm -m -any), recompute AF/AC/AN
#    FROM THE GENOTYPES (+fill-tags; so AN excludes any missing alleles and AF is exact for the
#    genotyped panel), SORT (left-alignment can reorder positions), bgzip + tabix. Resumable
#    (skips a chrom whose .tbi already exists). The source filename suffix differs for chrX.
echo "=== NORM $(date +%H:%M:%S) ==="
for c in $CHRS; do
  out="$B/VCFs/hg38_1000G2021/onekg.$c.norm.vcf.gz"
  [ -s "$out.tbi" ] && { echo "[$c] exists -> skip"; continue; }
  src=$(ls "$SRC"/*_"$c".filtered.*phased*.vcf.gz 2>/dev/null | head -1)
  [ -z "$src" ] && { echo "[$c] NO SOURCE VCF -> skip"; continue; }
  run "mkdir -p /tmp/sort_$c; bcftools norm -m -any -f $GENOME/$c.fa $src -Ou \
        | bcftools +fill-tags -Ou -- -t AF,AC,AN \
        | bcftools sort -T /tmp/sort_$c -Oz -o $out - && tabix -f -p vcf $out; rm -rf /tmp/sort_$c" \
    && echo "[$c] norm+sort done ($(basename $src))"
done

# 2) SAMPLES ID (exactly the genotyped panel): 4-col tab (#SAMPLE_ID POP SUPERPOP SEX). Sample
#    order/set from the VCF header (chr22), population metadata joined from the existing
#    1000G-2021 samplesID; default unknowns if a sample is missing there.
echo "=== SAMPLESID $(date +%H:%M:%S) ==="
run "zcat $B/VCFs/hg38_1000G2021/onekg.chr22.norm.vcf.gz | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n'" > "$B/samplesIDs/.onekg_order.txt"
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{ if(FNR>1){p[$1]=$2;s[$1]=$3;x[$1]=$4} next }
  { if($1 in p) print $1,p[$1],s[$1],x[$1]; else print $1,"1000G","unknown","unknown" }' \
  "$SID_META" "$B/samplesIDs/.onekg_order.txt" > "$B/samplesIDs/.body.txt"
{ printf '#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n'; cat "$B/samplesIDs/.body.txt"; } \
  > "$B/samplesIDs/hg38_1000G2021.samplesID.txt"
rm -f "$B/samplesIDs/.onekg_order.txt" "$B/samplesIDs/.body.txt"
echo "hg38_1000G2021.samplesID.txt" > "$B/onekg.samplesID.config.txt"
echo "samplesID rows (genotyped panel): $(($(wc -l < "$B/samplesIDs/hg38_1000G2021.samplesID.txt") - 1))  [expect 3202 -> AN 6404]"

# 3) BUILD the variant index. RAW registry + SNP+indel ON. Emits SNP index + _INDELS + registry/
#    genotypes/indel_genotypes/log_indels tiers + variant_count.json. Resumable.
echo "=== BUILD-INDEX-ONLY (RAW, NRG, SNP+indel, dev overlay) $(date +%H:%M:%S) ==="
run "export CRISPRME_REGISTRY_COMPRESS=0 CRISPRME_INDEL_SNP=1; cd $B; \
     crisprme.py build-index-only --genome $GENOME --pam $PAM --bDNA 2 --bRNA 2 --thread $THREADS \
       --vcf $B/VCFs/hg38_1000G2021 --samplesID $B/onekg.samplesID.config.txt --path $B"
echo "BUILD_EXIT=$? $(date +%H:%M:%S)"
echo "index: $B/genome_library/NRG_3_hg38+hg38_1000G2021 (+ _INDELS)"

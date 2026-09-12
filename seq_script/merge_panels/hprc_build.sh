#!/bin/bash
# HPRC pangenome variant-index build (as-run, resumable) — the genotyped, PHASED
# third production index NRG_3_hg38+hg38_HPRC.
#
# WHAT THIS IS: a standard CRISPRme variant-aware index built by ENRICHING hg38 with the
# HPRC Minigraph-Cactus pangenome, decomposed to a per-chromosome VCF with `vg deconstruct`
# -> `vcfbub` -> `vcfwave` (INFO carries LV/PS/CONFLICT snarl fields; FILTER is '.'; GT is
# diploid and PHASED, "a|b", with graph-coverage half-missingness like "1|."). 232 samples
# = 231 population individuals + CHM13 (the T2T graph backbone). GRCh38, chr-prefixed, all
# chr1..22,X,Y. This is NOT Ann's `assembly-search` (which searches per-individual assembly
# FASTAs and takes no VCF) — it is orthogonal to that feature.
#
# WHY IT IS ITS OWN SCRIPT (vs merge_panels Mode 1/2): the HPRC vcfwave VCF is a SINGLE
# already-genotyped source, so there is no cross-source merge — just per-chrom normalization
# + AF fill + a samplesID, then the standard `build-index-only`. Because the genotypes are
# phased, the index supports CONFIRMED cis co-occurrence + per-sample carriers (like the
# 1000G-2021+HGDP index), unlike the sites-only mega.
#
# OUTPUT: genome_library/NRG_3_hg38+hg38_HPRC (+ _INDELS) + Dictionaries/{registry,genotypes,
# indel_genotypes,log_indels}_hg38_HPRC + samplesIDs/hg38_HPRC.samplesID.txt. Registry is RAW
# (uncompressed, CRISPRME_REGISTRY_COMPRESS=0) and SNP+indel co-occurrence is ON
# (CRISPRME_INDEL_SNP default) — consistent with the other two production indexes.
#
# VALIDATED (v2.5.4, ml007): build clean on all 24 contigs; a --per-sample search produced
# 68,057 CONFIRMED phased-cis off-targets, 855 SNP+indel and 140 SNP+SNP co-occurring rows,
# with exact joint AF and named cis carriers ground-truth-verified against the source GTs
# (e.g. HG03579 carries chr10:12896898 TG>T + chr10:12896900 T>A in cis, joint AF 1/464).
#
# PATHS are the as-run cluster paths — override via env. Requires the CRISPRme code in an
# apptainer SIF ($SIF) with bcftools/tabix/python3 + the crisprme conda env.
set -u
SIF=${SIF:-/srv/local/lp698/crisprme_v254.sif}                 # apptainer image w/ crisprme env
B=${B:-/srv/local/lp698/hprc_build}                            # build working dir (--path)
SRC=${SRC:-/data/pinello/SHARED_DATA/CRISPRme_data/variants_datasets/archive/hg38_HPRC_vcfwave}  # ALL.<chr>.HPRC.vcfwave.vcf.gz
GENOME=${GENOME:-$B/Genomes/hg38}                              # per-chrom hg38 fastas (+ .fai)
PAM=${PAM:-$B/PAMs/20bp-NRG-SpCas9.txt}                        # NRG (NAG+NGG) SpCas9 PAM
P1KG=${P1KG:-/srv/local/lp698/mode1_2021_gw/DATA/samplesIDs/hg38_1000G2021.samplesID.txt}  # 1000G-2021 panel for pop-metadata join
THREADS=${THREADS:-32}
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

mkdir -p "$B/VCFs/hg38_HPRC" "$B/samplesIDs"

run() { apptainer exec --writable-tmpfs -B /srv/local/lp698 -B /data/pinello "$SIF" bash -c "export PATH=/opt/conda/bin:\$PATH; $*"; }

# 1) NORMALIZE each chrom: split multiallelics + left-align (norm -m -any), then recompute
#    AF/AC/AN FROM THE GENOTYPES (+fill-tags; the vcfwave VCF ships no AF and AN must exclude
#    the half-missing "." alleles), then SORT (left-alignment can reorder positions -> not
#    tabix-able otherwise; add-variants tolerates unsorted input but sorted is safer + indexable),
#    then bgzip + tabix. Resumable (skips a chrom whose .tbi already exists).
echo "=== NORM $(date +%H:%M:%S) ==="
for c in $CHRS; do
  out="$B/VCFs/hg38_HPRC/hprc.$c.norm.vcf.gz"
  [ -s "$out.tbi" ] && { echo "[$c] exists -> skip"; continue; }
  run "bcftools norm -m -any -f $GENOME/$c.fa $SRC/ALL.$c.HPRC.vcfwave.vcf.gz -Ou \
        | bcftools +fill-tags -Ou -- -t AF,AC,AN \
        | bcftools sort -Oz -o $out - && tabix -f -p vcf $out" && echo "[$c] norm+sort done"
done

# 2) SAMPLES ID (232 samples): 4-col tab format (#SAMPLE_ID POPULATION_ID SUPERPOPULATION_ID SEX),
#    sample order from the VCF header, population metadata joined from the 1000G-2021 panel where
#    the sample ID matches (most HPRC individuals are 1000G HG*/NA* samples); default POP=HPRC /
#    SUPERPOP=unknown / SEX=unknown for the rest (incl. HG002/HG005 and CHM13).
echo "=== SAMPLESID $(date +%H:%M:%S) ==="
run "zcat $SRC/ALL.chr22.HPRC.vcfwave.vcf.gz | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n'" > "$B/samplesIDs/.hprc_order.txt"
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{ if(FNR>1){p[$1]=$2;s[$1]=$3;x[$1]=$4} next }
  { if($1 in p) print $1,p[$1],s[$1],x[$1]; else print $1,"HPRC","unknown","unknown" }' \
  "$P1KG" "$B/samplesIDs/.hprc_order.txt" > "$B/samplesIDs/.body.txt"
{ printf '#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n'; cat "$B/samplesIDs/.body.txt"; } \
  > "$B/samplesIDs/hg38_HPRC.samplesID.txt"
rm -f "$B/samplesIDs/.hprc_order.txt" "$B/samplesIDs/.body.txt"
echo "hg38_HPRC.samplesID.txt" > "$B/hprc.samplesID.config.txt"
echo "samplesID rows: $(($(wc -l < "$B/samplesIDs/hg38_HPRC.samplesID.txt") - 1))"

# 3) BUILD the variant index. RAW registry (CRISPRME_REGISTRY_COMPRESS=0) + SNP+indel ON
#    (CRISPRME_INDEL_SNP default). Emits SNP index + _INDELS + registry/genotypes/indel_genotypes/
#    log_indels tiers + variant_count.json. Reuses any prior enriched genome/index (resumable).
echo "=== BUILD-INDEX-ONLY (RAW, NRG, SNP+indel) $(date +%H:%M:%S) ==="
run "export CRISPRME_REGISTRY_COMPRESS=0 CRISPRME_INDEL_SNP=1; cd $B; \
     crisprme.py build-index-only --genome $GENOME --pam $PAM --bDNA 2 --bRNA 2 --thread $THREADS \
       --vcf $B/VCFs/hg38_HPRC --samplesID $B/hprc.samplesID.config.txt --path $B"
echo "BUILD_EXIT=$? $(date +%H:%M:%S)"
echo "index: $B/genome_library/NRG_3_hg38+hg38_HPRC (+ _INDELS)"

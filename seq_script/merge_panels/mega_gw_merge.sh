#!/bin/bash
# Genome-wide 5-source "mega" sites-only merge (Mode 2), all 24 contigs, resumable.
# This is the as-run orchestrator that produced the shipped NRG_3_hg38+hg38_mega source
# VCFs (mega.<chr>.afmax.vcf.gz). It is the whole-genome driver around the per-chromosome
# logic in merge_mega_sites.sh; see merge_panels/README.md (Mode 2) + docs/DESIGN_mega_index.md.
#
# INPUTS (override for your environment — the defaults below are the Pinello-lab cluster
# paths the shipped index was built from; the 5 sources are private/controlled-access —
# gnomAD v4.1, TOPMed freeze, All-of-Us tier require registered dbGaP/Terra access, 1000G-2021
# + HGDP are public):
#   $VCF_DATA : dir with one subdir per source, each holding per-chrom VCFs
#   $REFDIR   : per-chromosome hg38 reference FASTAs (chr<N>.fa)
#   $OUT      : output dir for mega.<chr>.afmax.vcf.gz + <chr>.done markers
# Requires bcftools + tabix + bgzip on PATH.
set -o pipefail
VCF_DATA=${VCF_DATA:-/data/pinello/SHARED_DATA/CRISPRme_data/variants_datasets_20260827}
REFDIR=${REFDIR:-./Genomes/hg38}
OUT=${OUT:-./mega_gw}
mkdir -p "$OUT"
MAF="AF<=0.001 || AF>=0.999"     # keep MAF>0.1% variants (min-AF screening bound)
# tag:subdir for each source (per-dataset AF is renamed AF_<tag>; AF_max added after merge)
ORDER="1000G2021:hg38_1000G_2021 HGDP:hg38_HGDP gnomAD:hg38_gnomAD TOPMed:hg38_TOPMed AoU:hg38_AoU"

find_vcf() { local dir=$1 c=$2; ls "$VCF_DATA/$dir"/*.vcf.gz 2>/dev/null | grep -E "chr${c}([^0-9]|$)" | head -1; }

do_chrom() {
  local cN=$1; local c=chr$cN; local W=$OUT/work_$c
  if [ -f "$OUT/$c.done" ]; then echo "[$c] done, skip"; return 0; fi
  rm -rf "$W"; mkdir -p "$W"
  local REF=$REFDIR/$c.fa
  [ -f "$REF" ] || { echo "[$c] MISSING REF $REF"; return 1; }
  local inputs=""
  for pair in $ORDER; do
    local tag=${pair%%:*} dir=${pair#*:}
    local vcf; vcf=$(find_vcf "$dir" "$cN")
    [ -z "$vcf" ] && { echo "[$c] no $tag"; continue; }
    echo "INFO/AF AF_${tag}" > "$W/ren_$tag.txt"
    # normalize (split multiallelics, left-align), drop genotypes (-G, sites-only),
    # keep only INFO/AF, rename AF -> AF_<tag>
    bcftools norm -m -any -f "$REF" "$vcf" -Ou 2>>"$W/err_$tag.log" \
      | bcftools view -e "$MAF" -G -Ou 2>>"$W/err_$tag.log" \
      | bcftools annotate -x "^INFO/AF" -Ou 2>>"$W/err_$tag.log" \
      | bcftools annotate --rename-annots "$W/ren_$tag.txt" -Oz -o "$W/$tag.$c.vcf.gz" 2>>"$W/err_$tag.log"
    if [ $? -ne 0 ] || [ ! -s "$W/$tag.$c.vcf.gz" ]; then echo "[$c] FAIL proc $tag"; return 1; fi
    tabix -f -p vcf "$W/$tag.$c.vcf.gz" 2>>"$W/err_$tag.log"
    inputs="$inputs $W/$tag.$c.vcf.gz"
  done
  [ -n "$inputs" ] || { echo "[$c] no datasets"; return 1; }
  bcftools merge -m none $inputs -Oz -o "$W/mega.$c.vcf.gz" 2>>"$W/err_merge.log" || { echo "[$c] FAIL merge"; return 1; }
  tabix -f -p vcf "$W/mega.$c.vcf.gz" 2>>"$W/err_merge.log"
  # AF_max = max per-dataset AF across the 5 sources
  bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_1000G2021\t%INFO/AF_HGDP\t%INFO/AF_gnomAD\t%INFO/AF_TOPMed\t%INFO/AF_AoU\n" "$W/mega.$c.vcf.gz" 2>/dev/null \
    | awk -F"\t" 'BEGIN{OFS="\t"}{m=0;for(i=5;i<=9;i++){if($i!="."){v=$i+0;if(v>m)m=v}};print $1,$2,$3,$4,m}' | bgzip > "$W/afmax.tsv.gz"
  tabix -s1 -b2 -e2 -f "$W/afmax.tsv.gz"
  echo '##INFO=<ID=AF_max,Number=A,Type=Float,Description="Max per-dataset AF across the 5 merged datasets">' > "$W/afmax.hdr"
  bcftools annotate -a "$W/afmax.tsv.gz" -h "$W/afmax.hdr" -c CHROM,POS,REF,ALT,INFO/AF_max "$W/mega.$c.vcf.gz" -Oz -o "$OUT/mega.$c.afmax.vcf.gz" 2>>"$W/err_merge.log" || { echo "[$c] FAIL afmax"; return 1; }
  tabix -f -p vcf "$OUT/mega.$c.afmax.vcf.gz"
  local n; n=$(bcftools index -n "$OUT/mega.$c.afmax.vcf.gz" 2>/dev/null)
  echo "[$c] DONE sites=$n" | tee "$OUT/$c.done"
  rm -rf "$W"
}
export -f do_chrom find_vcf; export VCF_DATA REFDIR OUT MAF ORDER
echo "MEGA GW MERGE START $(date -u)"
printf "%s\n" 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y | xargs -P "${PARALLEL:-6}" -I{} bash -c "do_chrom {}"
echo "MEGA GW MERGE END $(date -u)"
tot=0; for c in $(seq 1 22) X Y; do [ -f "$OUT/chr$c.done" ] && { cat "$OUT/chr$c.done"; n=$(grep -oE "sites=[0-9]+" "$OUT/chr$c.done"|cut -d= -f2); tot=$((tot+n)); }; done
echo "TOTAL mega sites = $tot"

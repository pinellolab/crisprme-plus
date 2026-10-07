#!/bin/bash
# gnomAD HGDP+1kGP RE-PHASED single-source, PHASED variant-index build (as-run, resumable) —
# the candidate NEW DEFAULT "real-haplotypes" index NRG_3_hg38+hg38_HGDP1kGP.
#
# WHAT THIS IS: a standard CRISPRme variant-aware index built by ENRICHING hg38 with the gnomAD
# v3.1.2 HGDP+1kGP jointly-called callset, RE-PHASED by us with SHAPEIT5 (phase_common scaffold
# + phase_rare) to KEEP SINGLETONS (the public phased_haplotypes_v2 release drops them at MAC>=2).
# 4,091 samples (929 HGDP + 3,162 1kGP post-QC), GRCh38, chr-prefixed, chr1..22 + chrX (NO chrY/chrM
# -> chrY is reference-only in the index, as in every 1000G/HGDP index). Because BOTH cohorts are
# co-phased, the whole panel is genotyped-phased: HGDP off-targets now flip PUTATIVE -> CONFIRMED
# cis with named per-sample carriers + exact joint AF (the old NRG_3_hg38+hg38_1000G2021_HGDP mixed
# a phased 1kGP with an UNphased HGDP -> HGDP was PUTATIVE-only). 144,858,810 records total.
#
# UPSTREAM DEPENDENCY: this recipe starts from the 23 RE-PHASED BCFs
# ($RP/hgdp1kg.chr{1..22,X}.rephased.bcf). The SHAPEIT5 re-phase itself (dense callset download,
# per-chrom phase_common+phase_rare, the chrX 3-way PAR split + --haploids ploidy fix, the full
# validation vs the public release) is documented separately in
# docs/HGDP_1KGP_PHASED_INDEX_METHODS.md Section 6. Re-phase verified CLEAN: 23/23 read-OK, 4,091
# identical samples, 0 unphased genome-wide (full-file census), singletons kept.
#
# CORRECTNESS NOTES:
#  * samplesID is DERIVED FROM the atgu/gnomAD-PRESCRIBED metadata file gnomad_meta_updated.tsv
#    (release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/ -- the exact file the HGDP+1kGP
#    tutorial notebooks nb1-nb5 load), in VCF-header (sample) order: POP = hgdp_tgp_meta.Population,
#    SUPERPOP = hgdp_tgp_meta.Genetic.region (the 7 geographic regions AFR/AMR/CSA/EAS/EUR/MID/OCE),
#    SEX = sex_imputation.is_female (true->female, false->male). We deliberately AVOID the bare
#    convenience columns `population` (relabels 98 samples) and `sex` (NA for all 1kGP samples).
#    SEX must be lowercase male/female -- tier0_registry.make_chr_ploidy lowercases+matches "male".
#  * chrX non-PAR MALES are hemizygous: the re-phased BCF encodes them HOM-DIPLOID (x|x, SHAPEIT5
#    --haploids convention, AN=8182) but CRISPRme's chrX ploidy model (ploidy_of_for_chrom) treats
#    SEX=male as HAPLOID. We therefore `bcftools +fixploidy` non-PAR males x|x -> x BEFORE norm so
#    AC/AN are consistent + biologically correct: non-PAR AN=5990 (2192 males*1 + 1899 females*2),
#    PAR AN=8182 (all diploid). Verified: every gnomAD-male is hom on non-PAR (0 hets) so the
#    collapse is unambiguous; the 3 het-phased-male/metadata-female edge samples stay diploid.
#
# OUTPUT: genome_library/NRG_3_hg38+hg38_HGDP1kGP (+ _INDELS) + Dictionaries/{registry,genotypes,
# indel_genotypes,log_indels}_hg38_HGDP1kGP + variant_count.json + samplesIDs/hg38_HGDP1kGP.samplesID.txt.
# Registry RAW (CRISPRME_REGISTRY_COMPRESS=0), SNP+indel co-occurrence ON (CRISPRME_INDEL_SNP) --
# consistent with the other production indexes. The v2.6.2 SIF already stamps the data_type/phased
# manifest (genotyped-phased), so NO dev overlay is needed (unlike onekg2021_build.sh).
#
# PATHS are as-run cluster paths; override via env.
set -u
SIF=${SIF:-/srv/local/lp698_crisprme_v262_cleanroom/crisprme.sif}       # apptainer image w/ crisprme v2.6.2 env
# PAR-aware chrX ploidy overlay: 4 PostProcess files carrying the position-routed
# ploidy (PAR1/PAR2 diploid for all, non-PAR haploid-male) so chrX AF uses the
# correct denominator (PAR AN=8182, non-PAR AN=5990). Bound over the SIF's copies
# until a CRISPRme+ release bakes the fix into the image. Minimal delta vs v2.6.2
# (autosomes/chrY/mega byte-identical). See docs/HGDP_1KGP_PHASED_INDEX_METHODS.md Sec 9.
OV=${OV:-/srv/local/lp698_hgdp1kgp_rephase/par_overlay/PostProcess}      # PAR-aware PostProcess
B=${B:-/srv/local/lp698_hgdp1kgp_rephase}                               # build working dir (--path)
RP=${RP:-$B}                                                            # re-phased BCFs: $RP/hgdp1kg.<chr>.rephased.bcf
GENOME=${GENOME:-/srv/local/lp698/mode1_2021_gw/DATA/Genomes/hg38}      # per-chrom hg38 fastas (+ .fai)
PAM=${PAM:-/srv/local/lp698_crisprme_v262_cleanroom/PAMs/20bp-NRG-SpCas9.txt}  # NRG (NAG+NGG) SpCas9 PAM
META=${META:-/srv/local/lp698_hgdp1kgp_verify/gnomad_meta_updated.tsv}  # PRESCRIBED gnomAD metadata (v2/updated)
THREADS=${THREADS:-64}
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
# GRCh38 chrX non-PAR region (PAR1 ends 2781479, PAR2 starts 155701383); non-PAR males -> haploid.
NONPAR_FROM=2781480; NONPAR_TO=155701382

mkdir -p "$B/VCFs/hg38_HGDP1kGP" "$B/samplesIDs" "$B/tmp" "$B/sorttmp"

# Known-good SIF invocation: plain `bash -c` (NOT -l, which would source the bind-mounted host
# ~/.bashrc and clobber PATH back to the host), explicit /opt/conda/bin on PATH, and a REAL on-disk
# /tmp bind ($B/tmp) -- the in-memory tmpfs overflows when bcftools sort spills the big chromosomes.
run() { apptainer exec -B /srv/local -B "$B/tmp:/tmp" \
  "$SIF" bash -c "export PATH=/opt/conda/bin:\$PATH TMPDIR=/tmp; $*"; }

# 0) samplesID derivation (prescribed metadata, VCF-header order) + chrX ploidy/sex files.
echo "=== SAMPLESID + chrX ploidy/sex $(date +%H:%M:%S) ==="
run "bcftools query -l $RP/hgdp1kg.chr22.rephased.bcf" > "$B/samplesIDs/.order.txt"
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR{ if(FNR>1){ pop[$1]=($176==""?"unknown":$176); sup[$1]=($177==""?"unknown":$177);
                       f=tolower($80); sex[$1]=(f=="true"||f=="1"?"female":(f=="false"||f=="0"?"male":"unknown")) } next }
  { if($1 in pop) print $1,pop[$1],sup[$1],sex[$1]; else print $1,"unknown","unknown","unknown" }' \
  "$META" "$B/samplesIDs/.order.txt" > "$B/samplesIDs/.body.txt"
{ printf '#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n'; cat "$B/samplesIDs/.body.txt"; } \
  > "$B/samplesIDs/hg38_HGDP1kGP.samplesID.txt"
rm -f "$B/samplesIDs/.order.txt" "$B/samplesIDs/.body.txt"
echo "hg38_HGDP1kGP.samplesID.txt" > "$B/hgdp1kgp.samplesID.config.txt"
echo "samplesID rows: $(($(wc -l < "$B/samplesIDs/hg38_HGDP1kGP.samplesID.txt") - 1))  [expect 4091]"
# chrX +fixploidy inputs: non-PAR males haploid; PAR males + all females diploid (fall-through)
printf 'chrX\t%s\t%s\tM\t1\n*\t*\t*\tM\t2\n*\t*\t*\tF\t2\n' "$NONPAR_FROM" "$NONPAR_TO" > "$B/chrX.ploidy.txt"
awk -F'\t' 'NR>1{print $1"\t"(($4=="male")?"M":"F")}' "$B/samplesIDs/hg38_HGDP1kGP.samplesID.txt" > "$B/chrX.sex.txt"

# 1) VCF-PREP each re-phased BCF -> build-ready VCF. Strip the mis-declared PP FORMAT field, split
#    multiallelics + left-align against the reference (norm -m -any -f), recompute AF/AC/AN from the
#    genotypes (+fill-tags), SORT (left-alignment can reorder), bgzip + tabix. chrX additionally runs
#    +fixploidy FIRST (non-PAR males x|x -> x). Resumable (skips a chrom whose .tbi exists).
echo "=== VCF-PREP $(date +%H:%M:%S) ==="
for c in $CHRS; do
  out="$B/VCFs/hg38_HGDP1kGP/hgdp1kgp.$c.norm.vcf.gz"
  [ -s "$out.tbi" ] && { echo "[$c] exists -> skip"; continue; }
  in="$RP/hgdp1kg.$c.rephased.bcf"
  [ -s "$in" ] || { echo "[$c] NO re-phased BCF -> skip"; continue; }
  if [ "$c" = "chrX" ]; then
    run "mkdir -p /tmp/sort_$c; bcftools view $in -Ou \
          | bcftools +fixploidy -Ou -- -p $B/chrX.ploidy.txt -s $B/chrX.sex.txt \
          | bcftools annotate -x FORMAT/PP -Ou \
          | bcftools norm -m -any -f $GENOME/$c.fa -Ou \
          | bcftools +fill-tags -Ou -- -t AN,AC,AF \
          | bcftools sort -T /tmp/sort_$c -Oz -o $out - && tabix -f -p vcf $out; rm -rf /tmp/sort_$c"
  else
    run "mkdir -p /tmp/sort_$c; bcftools view $in -Ou \
          | bcftools annotate -x FORMAT/PP -Ou \
          | bcftools norm -m -any -f $GENOME/$c.fa -Ou \
          | bcftools +fill-tags -Ou -- -t AN,AC,AF \
          | bcftools sort -T /tmp/sort_$c -Oz -o $out - && tabix -f -p vcf $out; rm -rf /tmp/sort_$c"
  fi
  echo "[$c] prep done"
done

# 2) BUILD the variant index. RAW registry + SNP+indel ON. Emits SNP index + _INDELS + registry/
#    genotypes/indel_genotypes/log_indels tiers + variant_count.json. Resumable.
echo "=== BUILD-INDEX-ONLY (RAW, NRG, bDNA2/bRNA2, SNP+indel) $(date +%H:%M:%S) ==="
run "export CRISPRME_REGISTRY_COMPRESS=0 CRISPRME_INDEL_SNP=1; cd $B; \
     crisprme.py build-index-only --genome $GENOME --pam $PAM --bDNA 2 --bRNA 2 --thread $THREADS \
       --vcf $B/VCFs/hg38_HGDP1kGP --samplesID $B/hgdp1kgp.samplesID.config.txt --path $B"
echo "BUILD_EXIT=$? $(date +%H:%M:%S)"
echo "index: $B/genome_library/NRG_3_hg38+hg38_HGDP1kGP (+ _INDELS)"

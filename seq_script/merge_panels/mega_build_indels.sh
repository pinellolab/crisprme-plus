#!/bin/bash
# Mega sites-only INDELS build (all 24 contigs, resumable) — the caller of
# PostProcess/synth_sites_gt.py that makes the sites-only mega index's indels SEARCHABLE.
#
# WHY THIS EXISTS: the CRISPRitz enricher only materializes a fake-indel contig when a VCF
# SAMPLE column carries the indel. A sites-only (genotype-stripped) mega VCF has no sample
# columns, so the _INDELS genome would be EMPTY and indels unsearchable. Per chromosome this
# script: (1) synth_sites_gt.py adds ONE synthetic 0/1 "MEGA" pseudo-sample at every indel
# site; (2) build-index-only builds the fake-indel genome + INDELS index + log_indels via the
# classic pipeline; (3) STRIP the MEGA pseudo-sample from log_indels (col2 -> empty) so the
# shipped store is clean sites-only; (4) drop the fake-indel genome + INDELS bins + stripped
# log_indels into the near-complete mega index (which already has the SNP index + Tier-0
# registry + indel_af sidecar). Resumable (skips a chr whose INDELS bins + real log_indels
# already exist). See merge_panels/README.md (Mode 2) + docs/DESIGN_mega_index.md.
#
# PATHS below are the as-run cluster paths — override for your environment. Requires the
# CRISPRme code checkout ($CK), an apptainer SIF with the CRISPRme conda env ($SIF),
# tabix/python3 in the SIF.
set -u
CK=${CK:-/srv/local/lp698/crisprme}                # CRISPRme code checkout (crisprme.py + PostProcess/)
SIF=${SIF:-/srv/local/lp698/crisprme.sif}          # apptainer image with the crisprme conda env
MC=${MC:-/srv/local/lp698/mega_cleanroom}          # near-complete sites-only mega index (destination)
SRC=${SRC:-/srv/local/lp698/mega_gw}               # per-chr merged sites VCFs mega.chrN.afmax.vcf.gz (from mega_gw_merge.sh)
SCR=${SCR:-/srv/local/lp698/mega_gw_indel_scratch} # per-chr scratch builds
PAM=${PAM:-$(find /srv/local/lp698 -name "20bp-NRG-SpCas9.txt" 2>/dev/null | head -1)}
mkdir -p "$SCR" "$MC/Genomes/hg38+hg38_mega_INDELS" "$MC/genome_library/NRG_3_hg38+hg38_mega_INDELS"
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
strip_and_place() {  # $1=chr  $2=build_root
  local c=$1 W=$2
  local L; L=$(find "$W/Dictionaries/log_indels_hg38_mega" -name "log${c}.txt.gz" 2>/dev/null | head -1)
  [ -z "$L" ] && { echo "  [$c] NO log_indels in build -> SKIP"; return 1; }
  zcat "$L" | awk -F'\t' 'NR==1{print;next}{$2="";print}' OFS='\t' | gzip > "$MC/Dictionaries/log_indels_hg38_mega/log${c}.txt.gz"
  cp -f "$W"/Genomes/hg38+hg38_mega_INDELS/*fakechr${c#chr}* "$MC/Genomes/hg38+hg38_mega_INDELS/" 2>/dev/null
  cp -rf "$W"/genome_library/NRG_3_hg38+hg38_mega_INDELS/*fakechr${c#chr}* "$MC/genome_library/NRG_3_hg38+hg38_mega_INDELS/" 2>/dev/null
  local nb; nb=$(find "$MC/genome_library/NRG_3_hg38+hg38_mega_INDELS" -path "*fakechr${c#chr}*" -name "*.bin" 2>/dev/null | wc -l)
  echo "  [$c] placed: $(zcat "$MC/Dictionaries/log_indels_hg38_mega/log${c}.txt.gz"|wc -l) log recs (MEGA stripped), $nb INDELS bin(s)"
}
for c in $CHRS; do
  n=${c#chr}
  DONE=$(find "$MC/genome_library/NRG_3_hg38+hg38_mega_INDELS" -path "*fakechr${n}*" -name "*.bin" 2>/dev/null | wc -l)
  RLOG=$(zcat "$MC/Dictionaries/log_indels_hg38_mega/log${c}.txt.gz" 2>/dev/null | wc -l)
  if [ "$DONE" -ge 1 ] && [ "$RLOG" -ge 2 ]; then echo "[$c] already done ($DONE bins, $RLOG log recs) -> resume-skip"; continue; fi
  VCF="$SRC/mega.${c}.afmax.vcf.gz"
  [ -f "$VCF" ] || { echo "[$c] NO source VCF $VCF -> skip"; continue; }
  echo "[$c] BUILD start $(date +%H:%M:%S)"
  W="$SCR/$c"; rm -rf "$W"; mkdir -p "$W/Genomes/hg38" "$W/VCFs/hg38_mega" "$W/PAMs"
  cp "$MC/Genomes/hg38/${c}.fa" "$W/Genomes/hg38/" 2>/dev/null
  [ -f "$MC/Genomes/hg38/${c}.fa.fai" ] && cp "$MC/Genomes/hg38/${c}.fa.fai" "$W/Genomes/hg38/"
  cp "$PAM" "$W/PAMs/"
  apptainer exec -B /srv/local/lp698 "$SIF" /opt/conda/bin/python3 "$CK/PostProcess/synth_sites_gt.py" "$VCF" "$W/VCFs/hg38_mega/${c}.vcf.gz" MEGA >/dev/null 2>&1
  apptainer exec -B /srv/local/lp698 "$SIF" /opt/conda/bin/tabix -f -p vcf "$W/VCFs/hg38_mega/${c}.vcf.gz"
  apptainer exec --writable-tmpfs -B "$CK/crisprme.py":/opt/conda/bin/crisprme.py -B "$CK/PostProcess":/opt/conda/opt/crisprme/PostProcess -B /srv/local/lp698 "$SIF" \
    bash -c "export PATH=/opt/conda/bin:\$PATH; cd $W; crisprme.py build-index-only --genome $W/Genomes/hg38 --pam $W/PAMs/20bp-NRG-SpCas9.txt --bDNA 2 --bRNA 2 --thread 16 --vcf $W/VCFs/hg38_mega --path $W" > "$W/build.log" 2>&1
  echo "[$c] build exit=$? $(date +%H:%M:%S)"
  strip_and_place "$c" "$W" && rm -rf "$W/Genomes/run_"* "$W/genome_library/NRG_3_hg38" 2>/dev/null
done
echo "GW MEGA INDEL BUILD COMPLETE $(date +%H:%M:%S)"
echo "final INDELS bins: $(find "$MC/genome_library/NRG_3_hg38+hg38_mega_INDELS" -name '*.bin' 2>/dev/null | wc -l)"

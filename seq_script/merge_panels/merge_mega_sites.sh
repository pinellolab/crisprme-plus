#!/usr/bin/env bash
#
# merge_mega_sites.sh — merge N population resources into ONE per-chromosome
# SITES-ONLY panel (the "mega" all-source index), keeping each source's allele
# frequency separately plus a global AF_max. NO genotypes are kept.
#
# WHY A SEPARATE SCRIPT FROM merge_vcf_panels.sh:
#   merge_vcf_panels.sh is the GENOTYPED cis-capable merge (1000G-2021 + HGDP):
#   it keeps per-sample GT, recomputes a POOLED INFO/AF over the union panel with
#   `+fill-tags`, and feeds the genotype-based Tier-0 registry. That only works
#   when EVERY merged source is genotyped.
#   The mega merges heterogeneous AGGREGATE sources (gnomAD is AF-only; TOPMed
#   ships AN=0; All-of-Us is a single pseudo-sample) where a pooled AC/AN is NOT
#   computable and there are no cross-source haplotypes to preserve. So the mega
#   is SITES-ONLY: it keeps each source's ORIGINAL AF as AF_<src>, reports the
#   per-site max as AF_max, and the registry is built from those AFs directly
#   (PostProcess/build_mega_registry.py -> tier0_registry.compile_registry_from_info_af),
#   NOT from genotypes.
#
# WHAT THIS PRODUCES (per variant, provenance recoverable):
#   INFO/AF_<src>   each source's ORIGINAL allele frequency (missing where that
#                   source does not contain the variant). Per-alt (Number=A),
#                   because records are split to biallelic first.
#   INFO/AF_max     max over the present AF_<src> at the site (the global AF the
#                   mega registry's GLOBAL group carries; a pooled AC/AN is not
#                   computable across these sources).
#
# CORRECTNESS NOTES (why each step is where it is):
#   * `norm -m -any -f REF` FIRST: split multiallelics to biallelic AND left-align
#     against the reference, so the SAME indel from two sources gets an IDENTICAL
#     representation and `merge` actually combines them (not left as duplicate
#     records). AF is Number=A, so norm carries the correct per-alt AF — no
#     genotype read needed. (merge_vcf_panels.sh omits this; it relies on
#     pre-normalized genotyped sources.)
#   * MAF filter AFTER norm: `-e "AF<=0.001 || AF>=0.999"` keeps 0.001 < AF < 0.999
#     (MAF > 0.001, both rare-alt and rare-ref extremes dropped). Valid per-alt
#     only because norm already made AF scalar.
#   * `view -G`: drop all genotype columns -> sites-only.
#   * `annotate -x "^INFO/AF"`: keep ONLY INFO/AF (the "^" = complement). NB: this
#     MUST be "^INFO/AF" and NOT "INFO,^INFO/AF" — the bare "INFO" clause would
#     strip everything including AF (a bug fixed during development).
#   * `merge -m none`: union of sites; each record keeps its distinct AF_<src>
#     fields (no collision — they were renamed apart before merge).
#   * AF_max via query -> awk max (missing "." skipped) -> annotate: bcftools has
#     no native max-across-INFO, so compute it explicitly.
#
# USAGE:
#   merge_mega_sites.sh <out_vcf_dir> <chr> [ref_fasta_dir]
#   Configure SOURCES below (name | dir | filename glob with {C}). Produces
#   <out_vcf_dir>/mega.<chr>.afmax.vcf.gz (+ .tbi). Resumable (skips if .tbi
#   exists). Genome-wide: `printf '%s\n' chr1 .. chrY | xargs -P N -I{} \
#   merge_mega_sites.sh <out_vcf_dir> {}`.
#
# DEPENDENCIES: bcftools (norm/view/annotate/merge/query), bgzip, tabix. A
#   per-chrom reference FASTA <ref_dir>/<chr>.fa (+ .fai) for `norm -f`.
#
set -euo pipefail

OUTDIR=${1:?usage: merge_mega_sites.sh <out_vcf_dir> <chr> [ref_dir]}
CHR=${2:?usage: merge_mega_sites.sh <out_vcf_dir> <chr> [ref_dir]}
REFDIR=${3:-${CRISPRME_REF_DIR:?set ref_dir arg or CRISPRME_REF_DIR (per-chrom <chr>.fa)}}
BCFTOOLS=${BCFTOOLS:-bcftools}
DATA=${CRISPRME_VCF_DATA:-/data/pinello/SHARED_DATA/CRISPRme_data/variants_datasets_20260827}
MAF=${MAF_EXPR:-"AF<=0.001 || AF>=0.999"}

# name | directory | per-chr filename glob ({C} = chr token, matched so chr2 does
# not also catch chr20/22). Order fixes the AF_max query column order below.
SOURCES=(
  "1000G2021|$DATA/hg38_1000G_2021|*{C}*.vcf.gz"
  "HGDP|$DATA/hg38_HGDP|*{C}*.vcf.gz"
  "gnomAD|$DATA/hg38_gnomAD|*{C}*.vcf.gz"
  "TOPMed|$DATA/hg38_TOPMed|*{C}*.vcf.gz"
  "AoU|$DATA/hg38_AoU|*{C}*.vcf.gz"
)

mkdir -p "$OUTDIR"
OUT="$OUTDIR/mega.$CHR.afmax.vcf.gz"
[ -f "$OUT.tbi" ] && { echo "[skip] $CHR already merged"; exit 0; }
REF="$REFDIR/$CHR.fa"
[ -f "$REF" ] || { echo "[error] missing reference $REF"; exit 1; }
WORK=$(mktemp -d -p "${TMPDIR:-$OUTDIR}")
trap 'rm -rf "$WORK"' EXIT

prepped=(); names=()
for entry in "${SOURCES[@]}"; do
  IFS='|' read -r name dir glob <<<"$entry"
  # Match the chrom token exactly, ANCHORED TO THE BASENAME: "${CHR}([^0-9]|$)"
  # keeps chr2 from matching chr20/21/22; the trailing "[^/]*$" requires the token
  # be after the last slash, so a chr-token in a PARENT DIR name cannot false-match.
  # Ambiguity (>1 match, e.g. a sharded source or a stray sidecar) is FATAL, never
  # a silent head -1 truncation.
  hits=$(ls ${dir}/${glob//\{C\}/$CHR} 2>/dev/null | grep -E "${CHR}([^0-9]|\$)[^/]*\$")
  nhit=$(printf '%s' "$hits" | grep -c .)
  if [ "$nhit" -gt 1 ]; then
    echo "[error] $name: $nhit VCFs match $CHR (ambiguous) — refusing to pick one:" >&2
    printf '  %s\n' $hits >&2; exit 1
  fi
  vcf=$(printf '%s\n' "$hits" | head -1)
  [ -n "$vcf" ] || { echo "[warn] $name has no VCF for $CHR — skipping this source"; continue; }
  echo "[$name] $(basename "$vcf")"
  printf "INFO/AF AF_%s\n" "$name" > "$WORK/ren_$name.txt"
  pre="$WORK/$name.$CHR.vcf.gz"
  "$BCFTOOLS" norm -m -any -f "$REF" "$vcf" -Ou \
    | "$BCFTOOLS" view -e "$MAF" -G -Ou \
    | "$BCFTOOLS" annotate -x "^INFO/AF" -Ou \
    | "$BCFTOOLS" annotate --rename-annots "$WORK/ren_$name.txt" -Oz -o "$pre"
  "$BCFTOOLS" index -t "$pre"
  prepped+=("$pre"); names+=("$name")
done
[ "${#prepped[@]}" -ge 1 ] || { echo "[error] no sources for $CHR"; exit 1; }

# union merge -> per-dataset AF_<src> co-exist at each site
"$BCFTOOLS" merge -m none "${prepped[@]}" -Oz -o "$WORK/mega.noafmax.vcf.gz"
"$BCFTOOLS" index -t "$WORK/mega.noafmax.vcf.gz"

# AF_max = max over the present AF_<src> (missing "." skipped)
qfmt='%CHROM\t%POS\t%REF\t%ALT'; for n in "${names[@]}"; do qfmt="$qfmt\t%INFO/AF_$n"; done
# emit the max TWICE (col 5 -> INFO/AF_max, col 6 -> plain INFO/AF); annotate maps
# columns to targets positionally, so one value cannot fan out to two fields.
"$BCFTOOLS" query -f "$qfmt\n" "$WORK/mega.noafmax.vcf.gz" \
  | awk -F'\t' 'BEGIN{OFS="\t"}{m=""; for(i=5;i<=NF;i++){ if($i!="."){ v=$i+0; if(m==""||v>m)m=v } } if(m=="")m="."; print $1,$2,$3,$4,m,m}' \
  | bgzip > "$WORK/afmax.tab.gz"
tabix -s1 -b2 -e2 "$WORK/afmax.tab.gz"
# Annotate BOTH AF_max AND a plain INFO/AF (= AF_max). The plain AF is the single
# global frequency CRISPRitz `add-variants` embeds into the fake-indel contig
# metadata, so INDEL off-targets carry the global AF_max through the existing path
# with no post-analysis change (per-dataset indel AF comes from the separate
# build_mega_indel_af.py sidecar). For SNPs the dict-less registry is authoritative,
# so this plain AF only matters for indels / dict-based use.
printf '##INFO=<ID=AF_max,Number=A,Type=Float,Description="Max per-dataset AF across the merged sources">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Global allele frequency (= AF_max) for the enricher/fake-indel path">\n' > "$WORK/afmax.hdr"
"$BCFTOOLS" annotate -a "$WORK/afmax.tab.gz" -h "$WORK/afmax.hdr" \
  -c CHROM,POS,REF,ALT,INFO/AF_max,INFO/AF "$WORK/mega.noafmax.vcf.gz" -Oz -o "$OUT"
"$BCFTOOLS" index -t "$OUT"
echo "[done] $CHR: $("$BCFTOOLS" index -n "$OUT") sites, sources=${names[*]} -> $OUT"

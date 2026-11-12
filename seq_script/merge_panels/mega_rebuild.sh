#!/usr/bin/env bash
#
# mega_rebuild.sh — rebuild the sites-only "mega" index (CRISPRme 2.10.0, CRISPRitz 2.10.0)
# from the SAME merged VCFs the shipped mega was built from (mega.<chr>.afmax.vcf.gz, as
# mega_gw_merge.sh wrote them and vcfs/hg38_mega on HuggingFace holds them). Every part of
# the index then reads one set of records, the ones CRISPRitz's enricher reads:
#
#   1. the build input (PostProcess/prepare_mega_build_vcf.py): only the records whose
#      FILTER is PASS or '.', INFO/AF = AF_max for SNVs and indels alike (the global
#      frequency the enricher writes into the indel log; the as-run merge wrote no AF, so
#      every indel row of the shipped mega showed none), and one record of each REF>ALT,
#      the more frequent (the enricher merges two SNP records of one allele, but gives each
#      indel record a window and a log row of its own; the as-run merge wrote 8 indels on
#      two lines genome-wide, none on chr22 or chr14);
#   2. the presence genotype (PostProcess/synth_sites_gt.py): one pseudo-sample, MEGA,
#      carrying every ALT of every record, so the enricher windows every indel (it windows
#      an indel a sample carries; a sites-only VCF has none);
#   3. ONE build-index-only of the genome with those VCFs, without --samplesID: the
#      enriched genome (every PASS SNV allele coded, Q at four-allele positions), the SNP
#      index, the fake-indel genome and its INDELS index (index format 2), the indel log.
#      The SNP index and the enriched genome do not depend on the genotypes, so they are
#      those of a build of the sites-only VCFs themselves;
#   4. MEGA taken out again: the indel log's SAMPLES column emptied, the per-sample
#      dictionaries (MEGA in every entry) removed; build-index-only builds no genotype
#      store or phased indel store for a sites-only dataset;
#   5. the registry (PostProcess/build_mega_registry.py) and the per-source indel
#      frequencies (PostProcess/build_mega_indel_af.py) from the build input, both reading
#      it as the enricher does: FILTER PASS or '.', the SNP ALTs it codes, every allele of a
#      line, an allele once; no indel within the windows' left flank of the chromosome's
#      start, the flank read from the indel log just built;
#   6. samplesIDs/<dataset>.samplesID.txt, the file a search names (--samplesID): the registry's
#      sources as pseudo-samples, as the shipped install has it (one already there is kept);
#   7. the build input VCFs left in VCFs/<dataset>/ (sites-only), the checks of
#      PostProcess/verify_mega_chrom.py run on every chromosome (the source VCF read
#      independently against the registry, the indel frequencies, the indel log and the
#      enriched genome).
#
# USAGE
#   mega_rebuild.sh <install_dir> <merged_vcf_dir> [chrom ...]
#     <install_dir>     a CRISPRme working directory holding Genomes/$GENOME/<chrom>.fa (one
#                       FASTA per chromosome) and the PAM file; the index is built there
#     <merged_vcf_dir>  the folder of mega.<chrom>.afmax.vcf.gz (the source VCFs)
#     chrom ...         the chromosomes to build (default: every mega.*.afmax.vcf.gz there).
#                       Genomes/$GENOME must hold exactly these chromosomes' FASTAs: the
#                       index is that of the genome folder.
#   Environment (defaults):
#     PAM=PAMs/20bp-NRR-SpCas9.txt   (relative to <install_dir>; one PAM per run. For another
#                                     PAM, run build-index-only --pam <it> on the finished install,
#                                     which reuses its enrichment, as for any published dataset:
#                                     docs/PRECOMPUTED_INDEXES.md, "Building a new-PAM index")
#     BDNA=2 BRNA=2 THREADS=16 GENOME=hg38 DATASET=hg38_mega PRESENCE=MEGA
#     PP=<this repository>/PostProcess   CRISPRME=crisprme.py   PYTHON=python3
#     BGZIP=bgzip TABIX=tabix
#     RUN=''   a prefix for every command, e.g. 'apptainer exec --bind /data crisprme.sif'
#
#   Requires CRISPRitz 2.10.0 (format-2 indexes, every allele coded) and the CRISPRme of
#   this repository. Measured on chr22 and chr14 (METHODS section 2; CHANGELOG, 2.10.0).
#   Do not run it genome-wide before the release's readiness audit (docs/RELEASE_CHECKLIST).
set -euo pipefail

INSTALL=$(cd "${1:?usage: mega_rebuild.sh <install_dir> <merged_vcf_dir> [chrom ...]}" && pwd)
SRC=$(cd "${2:?usage: mega_rebuild.sh <install_dir> <merged_vcf_dir> [chrom ...]}" && pwd)
shift 2
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PP=${PP:-$(cd "$HERE/../../PostProcess" && pwd)}
PAM=${PAM:-PAMs/20bp-NRR-SpCas9.txt}
BDNA=${BDNA:-2}; BRNA=${BRNA:-2}; THREADS=${THREADS:-16}
GENOME=${GENOME:-hg38}; DATASET=${DATASET:-hg38_mega}; PRESENCE=${PRESENCE:-MEGA}
CRISPRME=${CRISPRME:-crisprme.py}; PYTHON=${PYTHON:-python3}
BGZIP=${BGZIP:-bgzip}; TABIX=${TABIX:-tabix}
RUN=${RUN:-}
say() { echo "[mega_rebuild $(date -u +%FT%TZ)] $*"; }
run() { $RUN "$@"; }

[ -f "$INSTALL/$PAM" ] || { echo "[error] no PAM file $INSTALL/$PAM" >&2; exit 1; }
[ -d "$INSTALL/Genomes/$GENOME" ] || { echo "[error] no genome folder $INSTALL/Genomes/$GENOME" >&2; exit 1; }
if [ "$#" -gt 0 ]; then
  CHROMS=("$@")
else
  CHROMS=()
  for f in "$SRC"/mega.*.afmax.vcf.gz; do
    [ -e "$f" ] || continue
    c=${f##*/mega.}; CHROMS+=("${c%.afmax.vcf.gz}")
  done
fi
[ "${#CHROMS[@]}" -ge 1 ] || { echo "[error] no mega.<chrom>.afmax.vcf.gz in $SRC" >&2; exit 1; }
for c in "${CHROMS[@]}"; do
  case "$c" in *[!A-Za-z0-9_.]*|'') echo "[error] bad chromosome name '$c'" >&2; exit 1;; esac
  [ -f "$SRC/mega.$c.afmax.vcf.gz" ] || { echo "[error] no $SRC/mega.$c.afmax.vcf.gz" >&2; exit 1; }
  [ -f "$INSTALL/Genomes/$GENOME/$c.fa" ] || { echo "[error] no $INSTALL/Genomes/$GENOME/$c.fa" >&2; exit 1; }
done
for f in "$INSTALL/Genomes/$GENOME"/*.fa; do
  c=$(basename "$f" .fa)
  printf '%s\n' "${CHROMS[@]}" | grep -qxF "$c" \
    || { echo "[error] $f is not one of the chromosomes built: the index is that of the genome folder" >&2; exit 1; }
done
VDIR="$INSTALL/VCFs/$DATASET"
STAGE="$INSTALL/.mega_rebuild_$DATASET"
D="$INSTALL/Dictionaries"
mkdir -p "$VDIR" "$STAGE/prepared" "$D"
if ls "$VDIR"/*.vcf* >/dev/null 2>&1; then
  echo "[error] $VDIR already holds VCFs: this recipe writes the dataset's VCFs itself" >&2; exit 1
fi

# 1-2. the build input, and its presence-genotype copy the enricher reads
for c in "${CHROMS[@]}"; do
  B="$STAGE/prepared/mega.$c.build.vcf.gz"
  run "$PYTHON" "$PP/prepare_mega_build_vcf.py" "$SRC/mega.$c.afmax.vcf.gz" "$B"
  run "$TABIX" -f -p vcf "$B"
  run "$PYTHON" "$PP/synth_sites_gt.py" "$B" "$VDIR/$c.vcf.gz" "$PRESENCE"
  run "$TABIX" -f -p vcf "$VDIR/$c.vcf.gz"
done
say "build input of ${#CHROMS[@]} chromosome(s) in $STAGE/prepared; presence genotype in $VDIR"

# 3. one build: enriched genome, SNP index, fake-indel genome, INDELS index, indel log
( cd "$INSTALL" && run "$CRISPRME" build-index-only --genome "Genomes/$GENOME" --pam "$PAM" \
    --bDNA "$BDNA" --bRNA "$BRNA" --thread "$THREADS" --vcf "VCFs/$DATASET" --path "$INSTALL" )
say "build-index-only done"

# 4. the presence genotype out again
for g in "genotypes_$DATASET" "indel_genotypes_$DATASET"; do
  if [ -e "$D/$g" ]; then
    echo "[error] $D/$g exists: build-index-only must build no genotype store for a sites-only dataset" >&2
    exit 1
  fi
done
rm -rf "${D:?}/dictionaries_${DATASET:?}"
L="$D/log_indels_$DATASET"
# (each inline python runs isolated, -I: it imports no module of the folder it runs in, which
# others may write; docs/TROUBLESHOOTING.md section 7, test_no_code_from_files)
for c in "${CHROMS[@]}"; do
  for lf in "$L/log$c.txt.gz" "$L/log$c.txt"; do
    [ -f "$lf" ] || continue
    if [ "${lf%.gz}" != "$lf" ]; then
      run "$PYTHON" -I - "$lf" <<'PY'
import gzip, os, sys
p = sys.argv[1]
with gzip.open(p, "rt") as fh, gzip.open(p + ".tmp", "wt") as out:
    out.write(fh.readline())
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > 1:
            f[1] = ""
        out.write("\t".join(f) + "\n")
os.replace(p + ".tmp", p)
PY
    else
      run "$PYTHON" -I - "$lf" <<'PY'
import os, sys
p = sys.argv[1]
with open(p) as fh, open(p + ".tmp", "w") as out:
    out.write(fh.readline())
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > 1:
            f[1] = ""
        out.write("\t".join(f) + "\n")
os.replace(p + ".tmp", p)
PY
    fi
  done
done
rm -f "$VDIR"/*.vcf.gz "$VDIR"/*.vcf.gz.tbi
say "presence genotype $PRESENCE taken out of the indel log; per-sample dictionaries removed"

# 5. registry and indel frequencies, from the build input; it becomes the dataset's VCFs
for c in "${CHROMS[@]}"; do
  B="$STAGE/prepared/mega.$c.build.vcf.gz"
  run "$PYTHON" "$PP/build_mega_registry.py" --vcf "$B" --chrom "$c" --out-dir "$D/registry_$DATASET"
  run "$PYTHON" "$PP/build_mega_indel_af.py" --vcf "$B" --chrom "$c" --out-dir "$D/indel_af_$DATASET" \
    --log-dir "$L"
  mv "$B" "$B.tbi" "$VDIR/"
done
rmdir "$STAGE/prepared" "$STAGE"
say "registry_$DATASET and indel_af_$DATASET written; VCFs/$DATASET holds the build input"

# 6. the samplesID file a search names (--samplesID) and its post-analysis reads: the
# registry's sources as pseudo-samples, as the shipped install has it; one already there is kept
SID="$INSTALL/samplesIDs/$DATASET.samplesID.txt"
if [ ! -e "$SID" ]; then
  mkdir -p "$INSTALL/samplesIDs"
  run "$PYTHON" -I - "$PP" "$SID" <<'PY'
import os, sys
sys.path.insert(0, sys.argv[1])
import build_mega_registry as b
with open(sys.argv[2] + ".tmp", "w") as out:
    out.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
    for ds in b.DEFAULT_DATASET_META:
        out.write("%s\t%s\t%s\tNA\n" % (ds, ds, ds))
os.replace(sys.argv[2] + ".tmp", sys.argv[2])
PY
  say "samplesIDs/$DATASET.samplesID.txt written (the registry's sources)"
fi

# 7. checks, the source VCFs read independently
run "$PYTHON" "$PP/verify_mega_chrom.py" "$SRC" "$D/registry_$DATASET" "$D/indel_af_$DATASET" \
  "${CHROMS[@]}" --log "$L" --enriched "$INSTALL/Genomes/$GENOME+$DATASET" \
  --reference "$INSTALL/Genomes/$GENOME" --build-vcf-dir "$VDIR"
say "done: $(ls -d "$INSTALL"/genome_library/*"$GENOME+$DATASET"* | xargs -n1 basename | tr '\n' ' ')"

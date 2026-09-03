#!/bin/bash
# 2.5.1 fast-mode benchmark: FAST vs SLOW new_simple_analysis.py on the REAL dense
# 2021+HGDP observed-haplotype path (chr22). Isolates the post-analysis (the thing
# --fast changes) by re-running it on the same cluster both ways. Runs INSIDE the SIF.
set -uo pipefail

CODE=/opt/conda/opt/crisprme/PostProcess
BENCH=/srv/local/lp698/fastbench
mkdir -p "$BENCH"
cd "$CODE" || exit 2

E2E=/srv/local/lp698/e2e2021/Results/TGCTTGGTCGGCACTGATAG_2021_feature_e2e/crispritz_targets
REF="$E2E/hg38_20bp-NRG-SpCas9.txt_guides.txt_5_2_2.targets.txt"
ALT="$E2E/hg38+hg38_1000G2021_HGDP_20bp-NRG-SpCas9.txt_guides.txt_5_2_2.targets.txt"
DICT=/srv/local/lp698/mode1_2021_gw/DATA/Dictionaries/dictionaries_hg38_1000G2021_HGDP/my_dict_chr22.json
GENOME=/srv/local/lp698/mode1_2021_gw/DATA/Genomes/hg38/chr22.fa
GUIDES=/srv/local/lp698/e2e2021/Results/TGCTTGGTCGGCACTGATAG_2021_feature_e2e/guides.txt
PAM=/srv/local/lp698/mode1_2021_gw/DATA/PAMs/20bp-NRG-SpCas9.txt
CHROM=chr22

echo "==== FASTBENCH start $(date -u) ===="
echo "code HEAD: $(cd /srv/local/lp698/fastbench_ck && git rev-parse --short HEAD 2>/dev/null) + overlay"
python3 -c "import sys; sys.path.insert(0,'$CODE'); import twopass_emit; print('twopass_emit import OK')" || exit 3

# --- Phase 0: produce the chr22 observed cluster (once) --------------------------
CLUSTER="$BENCH/$CHROM.total.cluster.txt"
if [ ! -s "$CLUSTER" ]; then
  echo "[phase0] grep $CHROM targets $(date -u)"
  LC_ALL=C grep -F -w "$CHROM" "$REF" > "$BENCH/refc.txt"
  LC_ALL=C grep -F -w "$CHROM" "$ALT" > "$BENCH/altc.txt"
  echo "[phase0] ref=$(wc -l < "$BENCH/refc.txt") alt=$(wc -l < "$BENCH/altc.txt") lines"
  echo "[phase0] extraction.sh $(date -u)"
  ./extraction.sh "$BENCH/refc.txt" "$BENCH/altc.txt" "$BENCH/$CHROM"
  rm -f "$BENCH/$CHROM.common_targets.txt"
  awk '{rg=$2; gsub("-","",rg); print $0"\tn\tn\tn\tn\t"rg"\tn\tn\tn"}' "$BENCH/$CHROM.semi_common_targets.txt" > "$BENCH/$CHROM.scm.txt"
  awk '{rg=$2; gsub("-","",rg); print $0"\tn\tn\tn\tn\t"rg"\tn\tn\tn"}' "$BENCH/$CHROM.unique_targets.txt"       > "$BENCH/$CHROM.uq.txt"
  cat "$BENCH/$CHROM.uq.txt" "$BENCH/$CHROM.scm.txt" > "$BENCH/$CHROM.total.txt"
  echo "[phase0] cluster.dict.py $(date -u)"
  ./cluster.dict.py "$BENCH/$CHROM.total.txt" 'no' 'True' 'True' "$GUIDES" 'total' 'orderChr'
  # cluster.dict.py writes <input w/ .total.txt -> .total.cluster.txt>
  [ -s "$BENCH/$CHROM.total.cluster.txt" ] || { echo "CLUSTER NOT PRODUCED"; ls -la "$BENCH"; exit 4; }
  rm -f "$BENCH/$CHROM.total.txt" "$BENCH/$CHROM.scm.txt" "$BENCH/$CHROM.uq.txt" \
        "$BENCH/$CHROM.semi_common_targets.txt" "$BENCH/$CHROM.unique_targets.txt" \
        "$BENCH/refc.txt" "$BENCH/altc.txt"
fi
echo "[phase0] cluster lines=$(wc -l < "$CLUSTER") $(date -u)"

# --- Phase 1: FAST (quick; validates the path) -----------------------------------
echo "[FAST] start $(date -u)"; F0=$(date +%s)
CRISPRME_FAST_MODE=1 ./new_simple_analysis.py "$GENOME" "$DICT" "$CLUSTER" "$PAM" "$BENCH/out_fast" 5 > "$BENCH/fast.log" 2>&1
FRC=$?; F1=$(date +%s)
echo "[FAST] end $(date -u) rc=$FRC elapsed=$((F1-F0))s rows=$(wc -l < "$BENCH/out_fast.bestmmblg.txt" 2>/dev/null)"
tail -2 "$BENCH/fast.log"

# --- Phase 2: SLOW (the enumeration; the thing --fast eliminates) -----------------
echo "[SLOW] start $(date -u)"; S0=$(date +%s)
./new_simple_analysis.py "$GENOME" "$DICT" "$CLUSTER" "$PAM" "$BENCH/out_slow" 5 > "$BENCH/slow.log" 2>&1
SRC=$?; S1=$(date +%s)
echo "[SLOW] end $(date -u) rc=$SRC elapsed=$((S1-S0))s rows=$(wc -l < "$BENCH/out_slow.bestmmblg.txt" 2>/dev/null)"
tail -2 "$BENCH/slow.log"

echo "==== TIMING: FAST=$((F1-F0))s  SLOW=$((S1-S0))s  speedup=$(python3 -c "print('%.1fx' % (($S1-$S0)/max(1,($F1-$F0))))") ===="

# --- Phase 3: no-miss + worst-possible comparison --------------------------------
./fastbench_compare.py "$BENCH/out_slow.bestmmblg.txt" "$BENCH/out_fast.bestmmblg.txt"
echo "==== FASTBENCH done $(date -u) ===="

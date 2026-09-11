#!/bin/bash
# Bounded, sequential (no-contention) fast-vs-slow slice benchmark: run the SNP
# post-analysis both ways on the SAME dense chr22 slice and compare row-count collapse,
# no-miss, worst-possible, and wall time. Runs INSIDE the SIF.
set -uo pipefail
CODE=/opt/conda/opt/crisprme/PostProcess
BENCH=/srv/local/lp698/fastbench
SLICE=/srv/local/lp698/fastbench_slice
N="${1:-50000}"
mkdir -p "$SLICE"
cd "$CODE" || exit 2
DICT=/srv/local/lp698/mode1_2021_gw/DATA/Dictionaries/dictionaries_hg38_1000G2021_HGDP/my_dict_chr22.json
GENOME=/srv/local/lp698/mode1_2021_gw/DATA/Genomes/hg38/chr22.fa
PAM=/srv/local/lp698/mode1_2021_gw/DATA/PAMs/20bp-NRG-SpCas9.txt

echo "==== SLICE BENCH start $(date -u)  N=$N ===="
head -n "$N" "$BENCH/chr22.total.cluster.txt" > "$SLICE/slice.cluster.txt"
echo "slice lines=$(wc -l < "$SLICE/slice.cluster.txt")"

echo "[SLICE FAST] start $(date -u)"; F0=$(date +%s)
CRISPRME_FAST_MODE=1 ./new_simple_analysis.py "$GENOME" "$DICT" "$SLICE/slice.cluster.txt" "$PAM" "$SLICE/out_fast" 5 > "$SLICE/fast.log" 2>&1
FRC=$?; F1=$(date +%s)
echo "[SLICE FAST] end $(date -u) rc=$FRC elapsed=$((F1-F0))s rows=$(wc -l < "$SLICE/out_fast.bestmmblg.txt" 2>/dev/null)"

echo "[SLICE SLOW] start $(date -u)"; S0=$(date +%s)
./new_simple_analysis.py "$GENOME" "$DICT" "$SLICE/slice.cluster.txt" "$PAM" "$SLICE/out_slow" 5 > "$SLICE/slow.log" 2>&1
SRC=$?; S1=$(date +%s)
echo "[SLICE SLOW] end $(date -u) rc=$SRC elapsed=$((S1-S0))s rows=$(wc -l < "$SLICE/out_slow.bestmmblg.txt" 2>/dev/null)"

echo "==== SLICE TIMING: FAST=$((F1-F0))s  SLOW=$((S1-S0))s  speedup=$(python3 -c "print('%.2fx' % (($S1-$S0)/max(1,($F1-$F0))))") ===="
./fastbench_compare.py "$SLICE/out_slow.bestmmblg.txt" "$SLICE/out_fast.bestmmblg.txt"
echo "==== SLICE BENCH done $(date -u) ===="

#!/bin/bash
# The SNP+indel oracle through the real CLI, inside a CRISPRme image:
#   run_oracle.sh <empty work dir> [<oracle dir>] [cas12a]
# generates the synthetic phased genome + VCF (oracle.py generate), builds its
# indexes and dictionaries with build-index-only, runs complete-search per-sample and
# population-level, then checks every indel row, and the carriers of every per-sample
# SNP-only row, against the per-haplotype truth
# (oracle.py check) and the ALT-PAM and no-N release gate (scripts/check_alt_pam.py).
# The default profile plants SpCas9 sites (3' PAM NRG); "cas12a" plants Cas12a sites
# (5' PAM TTTV) whose first base is one an insertion adds.
# Exits non-zero on the first failing step. A few minutes on 2 threads.
set -uo pipefail
wd=$1
here=${2:-$(cd "$(dirname "$0")" && pwd)}
profile=${3:-cas9}
pp=$(dirname "$(command -v crisprme.py)")/../opt/crisprme/PostProcess
[ -d "$pp" ] || pp=/opt/conda/opt/crisprme/PostProcess
gate=$here/../../scripts/check_alt_pam.py
mkdir -p "$wd" && cd "$wd" || exit 1
case "$profile" in
	cas9) python3 -I "$here/oracle.py" generate "$wd" || exit 1 ;;
	cas12a) python3 -I "$here/oracle.py" generate "$wd" --cas12a || exit 1 ;;
	*) echo "unknown profile $profile (cas9, cas12a)"; exit 2 ;;
esac
pam=PAMs/$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["pam_file"])' spec.json)
bgzip -f VCFs/synth_vcf/synth.chr21.vcf || exit 1
echo "=== build-index-only $(date -u)"
crisprme.py build-index-only --genome Genomes/synth --pam "$pam" \
	--bDNA 1 --bRNA 1 --vcf VCFs/synth_vcf --samplesID list_samplesID.txt --path "$wd" \
	--thread 2 || { echo "FAIL build-index-only"; exit 1; }
rc=0
for mode in per-sample population-level; do
	out=orc_${mode:0:3}
	echo "=== complete-search --$mode $(date -u)"
	crisprme.py complete-search --genome Genomes/synth --vcf list_vcf.txt \
		--samplesID list_samplesID.txt --guide guides.txt --pam "$pam" \
		--mm 4 --bDNA 1 --bRNA 1 --output "$out" --"$mode" --thread 2 \
		|| { echo "FAIL complete-search --$mode"; rc=1; continue; }
	err=Results/$out/log_error.txt
	[ -e "$err" ] || err=Results/$out/log_error_no_check.txt
	if [ ! -e "$err" ] || [ -s "$err" ]; then
		echo "FAIL $err is missing or not empty:"; cat "$err" 2>/dev/null; rc=1
	fi
	grep -h "indel-snp:" Results/"$out"/log_verbose.txt || true
	python3 -I "$gate" Results/"$out" || rc=1
done
python3 -I "$here/oracle.py" check "$wd" "$pp" Results/orc_per Results/orc_pop || rc=1
# The same rows through the independent real-data checker (scripts/check_indel_rows.py),
# which imports and copies no code of oracle.py or CRISPRme: its own VCF reader,
# haplotype builder and exhaustive aligner, driven only by what a user has. Soundness in
# both modes, and completeness in both: per sample every (sample, target) must be
# reported, while the population-level run (one row per cluster, no alternative
# alignments) is checked per cluster -- the checker reads that from the run's .alt_mode.
# On the oracle the planted two-indel sites O and Q (and the second records of W and V)
# come out as limit-two-indels and site P's phasing artefact as limit-overlapping-records,
# non-gating; every gating finding fails the oracle, as oracle.py's do, except the ones
# known_check_indel_rows.tsv lists (a CRISPRme defect reported upstream: they print as
# known-defect, and as known-not-seen once fixed).
chk=$here/../../scripts/check_indel_rows.py
if [ -f "$chk" ]; then
	for mode in soundness completeness; do
		for out in orc_per orc_pop; do
			echo "=== check_indel_rows --$mode $out $(date -u)"
			known=()
			# the listed defect is in the per-sample run's alternative alignments
			[ "$out/$mode" = orc_per/soundness ] && known=(--known "$here/known_check_indel_rows.tsv")
			python3 -I "$chk" --results "Results/$out" --reference Genomes/synth \
				--vcf VCFs/synth_vcf/synth.chr21.vcf.gz --pam PAMs/20bp-NRG-SpCas9.txt \
				--guide guides.txt --mm 4 --bDNA 1 --bRNA 1 --max-total-edits 6 \
				--chrom chr21 --mode "$mode" "${known[@]}" \
				--out "check_${out}_${mode}.tsv" || rc=1
		done
	done
else
	echo "FAIL: $chk not found; the independent real-data check cannot run"; rc=1
fi
[ "$rc" = 0 ] && echo "ORACLE PASS" || echo "ORACLE FAIL"
exit "$rc"

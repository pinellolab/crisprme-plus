#!/bin/bash
# The SNP+indel oracle through the real CLI, inside a CRISPRme image:
#   run_oracle.sh <empty work dir> [<oracle dir>]
# generates the synthetic phased genome + VCF (oracle.py generate), builds its
# indexes and dictionaries with build-index-only, runs complete-search per-sample and
# population-level, then checks every indel row against the per-haplotype truth
# (oracle.py check) and the ALT-PAM and no-N release gate (scripts/check_alt_pam.py).
# Exits non-zero on the first failing step. A few minutes on 2 threads.
set -uo pipefail
wd=$1
here=${2:-$(cd "$(dirname "$0")" && pwd)}
pp=$(dirname "$(command -v crisprme.py)")/../opt/crisprme/PostProcess
[ -d "$pp" ] || pp=/opt/conda/opt/crisprme/PostProcess
gate=$here/../../scripts/check_alt_pam.py
mkdir -p "$wd" && cd "$wd" || exit 1
python3 -I "$here/oracle.py" generate "$wd" || exit 1
bgzip -f VCFs/synth_vcf/synth.chr21.vcf || exit 1
echo "=== build-index-only $(date -u)"
crisprme.py build-index-only --genome Genomes/synth --pam PAMs/20bp-NRG-SpCas9.txt \
	--bDNA 1 --bRNA 1 --vcf VCFs/synth_vcf --samplesID list_samplesID.txt --path "$wd" \
	--thread 2 || { echo "FAIL build-index-only"; exit 1; }
rc=0
for mode in per-sample population-level; do
	out=orc_${mode:0:3}
	echo "=== complete-search --$mode $(date -u)"
	crisprme.py complete-search --genome Genomes/synth --vcf list_vcf.txt \
		--samplesID list_samplesID.txt --guide guides.txt --pam PAMs/20bp-NRG-SpCas9.txt \
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
[ "$rc" = 0 ] && echo "ORACLE PASS" || echo "ORACLE FAIL"
exit "$rc"

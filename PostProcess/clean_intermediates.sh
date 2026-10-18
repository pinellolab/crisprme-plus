#!/bin/bash
# Remove the bulky pre-integration intermediate files at the end of a
# complete-search run, once the finalized *_integrated_results.tsv and the
# report have been produced. Invoked by
# submit_job_automated_new_multiple_vcfs.sh. Factored out as its own script so
# the destructive behaviour can be unit-tested in isolation (test_clean_intermediates.py).
#
# Usage: clean_intermediates.sh <output_folder> <cicd_test>
#
# No-op (keeps everything) when CRISPRME_KEEP_INTERMEDIATES=True
# (crisprme.py complete-search --keep-intermediates) OR when <cicd_test> == "True"
# (complete-test / validate-benchmarks keep every intermediate for the brute-force gate).
#
# REMOVES only: <job>.altMerge.txt (the raw non-best alignments working table; not
# read after the report is built), the raw
# *_all_results_with_alternative_alignments.tsv (the .zip is kept), and the
# per-chromosome *.population_summary.tsv companions.
#
# KEPT regardless: <job>.bestMerge.txt -- the web results page reads it LIVE for the
# per-position cluster drill-down (pages/results_page.py cluster_page greps
# <job>.bestMerge.txt each time a cluster is opened), and a web run has no way to opt
# out of cleanup, so this file must persist for the job to stay browsable;
# *_integrated_results.tsv and .zip, *_report.zip,
# *_all_results_with_alternative_alignments.zip, the per-chrom co-occurrence /
# high-variant-density sidecars (snp_snp_cooc / indel_snp_cooc /
# high_variant_density_regions / indel_af), the per-guide/per-sample summaries,
# imgs/, and the sqlite .<job>.db web results-page index.
output_folder="$1"
cicd_test="$2"

if [ "${CRISPRME_KEEP_INTERMEDIATES:-False}" = "True" ] || [ "$cicd_test" = "True" ]; then
	exit 0
fi

_job=$(basename "${output_folder}")
echo -e "Cleaning intermediate files (set CRISPRME_KEEP_INTERMEDIATES=True to keep them)"
# NOTE: <job>.bestMerge.txt is deliberately NOT removed -- the web cluster drill-down
# greps it live (see header). Only altMerge (non-best, not read post-report) goes.
rm -f "${output_folder}/${_job}.altMerge.txt"
rm -f "${output_folder}"/*_all_results_with_alternative_alignments.tsv  # the .zip is kept
rm -f "${output_folder}"/*.population_summary.tsv

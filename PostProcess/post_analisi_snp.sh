#!/bin/bash

output_folder=$1
ref_folder=$2
ref_name=$(basename $2)
vcf_name=$3
guide_file=$4
guide_name=$(basename $4)
mm=$5
bDNA=$6
bRNA=$7
annotation_file=$8
annotation_name=$(basename $8)
pam_file=$9
pam_name=$(basename $9)
dict_folder=${10}
final_res=${11}
final_res_alt=${12}
key=${13}
# 2.10.0, option P: all (default), cluster or finish (scriptAnalisiNNN_v3.sh)
stage=${14:-all}

# reference and alternative crispritz targets files
targets_tsv_ref="$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_alt="$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"

if [[ ! -f "$targets_tsv_alt" ]]; then 
    touch "$targets_tsv_alt"  # keep consistency through following steps
fi

# start chromosome data processing
chrom=$key
targets_chrom_prefix="$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_${key}"
if [ "$stage" = "finish" ]; then
    # the steps after the analysis, once the pool has merged its pieces
    ./scriptAnalisiNNN_v3.sh "" "" "$targets_chrom_prefix" "$annotation_file" \
        "${dict_folder}/my_dict_${key}.json" "${ref_folder}/${key}.fa" "$mm" "$bDNA" "$bRNA" \
        "$guide_file" "$pam_file" "$output_folder" finish
    exit $?
fi
echo "Processing SNPs for ${chrom}"

# extract chrom-specific reference and alternative targets (the name as text: unquoted,
# the shell took '*', '?' and '[' in it for a pattern of the files of the folder this runs
# in, which grep then read, and grep took a first '-' for an option)
targets_tsv_ref_chrom="${targets_tsv_ref}.${chrom}"
targets_tsv_alt_chrom="${targets_tsv_alt}.${chrom}"
LC_ALL=C grep -F -w -e "$chrom" "$targets_tsv_ref" > "$targets_tsv_ref_chrom"
LC_ALL=C grep -F -w -e "$chrom" "$targets_tsv_alt" > "$targets_tsv_alt_chrom"

# perform targets analysis by chromosome (scores, annotation, etc.)
./scriptAnalisiNNN_v3.sh \
    "$targets_tsv_ref_chrom" \
    "$targets_tsv_alt_chrom" \
    "$targets_chrom_prefix" \
    "$annotation_file" \
    "${dict_folder}/my_dict_${key}.json" \
    "${ref_folder}/${key}.fa" \
    "$mm" "$bDNA" "$bRNA" \
    "$guide_file" \
    "$pam_file" \
    "$output_folder" \
    "$stage"

# remove chrom-specific targets tsv files
rm $targets_tsv_ref_chrom $targets_tsv_alt_chrom  

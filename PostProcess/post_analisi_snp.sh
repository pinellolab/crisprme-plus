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

# reference and alternative crispritz targets files
targets_tsv_ref="$output_folder/crispritz_targets/${ref_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"
targets_tsv_alt="$output_folder/crispritz_targets/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}.targets.txt"

if [[ ! -f "$targets_tsv_alt" ]]; then 
    touch "$targets_tsv_alt"  # keep consistency through following steps
fi

# start chromosome data processing
chrom=$key
echo "Processing SNPs for ${chrom}"

# extract chrom-specific reference and alternative targets (the name as text: unquoted,
# the shell took '*', '?' and '[' in it for a pattern of the files of the folder this runs
# in, which grep then read, and grep took a first '-' for an option)
targets_tsv_ref_chrom="${targets_tsv_ref}.${chrom}"
targets_tsv_alt_chrom="${targets_tsv_alt}.${chrom}"
LC_ALL=C grep -F -w -e "$chrom" "$targets_tsv_ref" > "$targets_tsv_ref_chrom"
LC_ALL=C grep -F -w -e "$chrom" "$targets_tsv_alt" > "$targets_tsv_alt_chrom"

# the sites-only density guard (sites_only_guard.py, CRISPRme 2.10.0): of a sites-only
# dataset (the mega), the candidates whose window spells 2^bits sequences or more
# (CRISPRME_SITES_ONLY_GUARD, 20 by default) leave the targets before the post-analysis
# reads them, and their regions are listed in the BED; a genotyped dataset's are untouched
python sites_only_guard.py snp "$targets_tsv_alt_chrom" "$(dirname "${dict_folder%/}")" \
    "$vcf_name" "$chrom" \
    "$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${mm}_${bDNA}_${bRNA}_${chrom}.snp.sites_only_density_guard.bed"

# perform targets analysis by chromosome (scores, annotation, etc.)
targets_chrom_prefix="$output_folder/${ref_name}+${vcf_name}_${pam_name}_${guide_name}_${annotation_name}_${mm}_${bDNA}_${bRNA}_${key}"
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
    "$output_folder"

# remove chrom-specific targets tsv files
rm $targets_tsv_ref_chrom $targets_tsv_alt_chrom  

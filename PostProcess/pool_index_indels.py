#!/usr/bin/env python

import os
import sys
from multiprocessing import Pool
from datetime import datetime

indels_folder = sys.argv[1]
pam_file = sys.argv[2]
# true_pam / bMax here are the RESOLVED variant index's PAM and N (its <PAM>_<N>
# prefix), used for the per-chrom leaf name (<PAM>_<N>_fake<chr>) that the search
# path expects. IMPORTANT: crispritz index-genome AUTOMATICALLY prepends
# "<PAM>_<bMax>_" to the output-name's outer dir (that's how the SNP index gets its
# prefix from a name passed WITHOUT one). So the outer name must be passed WITHOUT a
# manual prefix — passing "<PAM>_<N>_<ref>+<vcf>_INDELS" here would make crispritz
# produce a DOUBLE-prefixed "<PAM>_<N>_<PAM>_<N>_..._INDELS" (a regression that broke
# the on-demand indel build / complete-test). Outer = "<ref>+<vcf>_INDELS" ->
# crispritz -> "<PAM>_<N>_<ref>+<vcf>_INDELS", matching detection + search.
true_pam = sys.argv[3]
ref_name = sys.argv[4]
vcf_name = sys.argv[5]
bMax = sys.argv[6]
ncpus = sys.argv[7]


def index_indels(chrom):
    print("Indexing INDELs in", chrom)
    os.system(
        f"crispritz.py index-genome {ref_name}+{vcf_name}_INDELS/{true_pam}_{bMax}_fake{chrom} {indels_folder}/fake_{vcf_name}_{chrom} {pam_file} -bMax {bMax} -th 1  >/dev/null"
    )  # {indels_folder}/fake_{vcf_name}_{chrom}
    print("Indexing ended for INDELs in", chrom)


chrs = []
for f in os.listdir(indels_folder):
    if "chr" in f:
        chrs.append(f.split("_")[-1])

with Pool(processes=int(ncpus)) as pool:
    pool.map(index_indels, chrs)


# os.system('echo "Indexing INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')

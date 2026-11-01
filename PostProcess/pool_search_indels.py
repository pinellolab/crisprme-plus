#!/usr/bin/env python

import gzip
import os
import sys
from multiprocessing import Pool
from datetime import datetime


def _chrom_from_vcf(vcf_path: str) -> str:
    """Return the chromosome name from the first data line of a VCF (plain or gzipped)."""
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t")[0]
    raise ValueError(f"No data lines found in VCF: {vcf_path}")


def _normalize_chrom(chrom: str) -> str:
    """Ensure chromosome name has chr prefix (e.g. '22' -> 'chr22').
    Some VCF datasets (e.g. 1000G GRCh38) store chromosomes without the
    prefix in the CHROM field, while the genome indices use 'chr'-prefixed names.
    """
    return chrom if chrom.startswith("chr") else "chr" + chrom


# CRISPRitz's targets header, written on its own when every chromosome was skipped,
# so the stages downstream still read a well-formed (empty) targets file.
_TARGETS_HEADER = (
    "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
    "Direction\tMismatches\tBulge Size\tTotal"
)


def _logged_no_indels(log_dir: str, chrom: str) -> bool:
    """True only when the enricher's indel log for ``chrom`` exists and is empty.

    ``Dictionaries/log_indels_<vcf>/log<chrom>.txt[.gz]`` holds one row per indel
    window the enricher wrote into that chromosome's fake-indel genome, under a
    ``CHR`` header, and ships with every variant index. Header-only means that
    genome is empty -- a chromosome whose VCF carries no indel in any sample.
    buildTST cannot index an empty genome, so its sub-index holds no ``.bin``, and
    CRISPRitz 2.9.1 refuses to search an index directory without one: rightly for
    a broken index, wrongly here, where no indel off-target is the answer. Before
    2.9.1 the same state silently produced an empty result.

    A missing log returns False -- unknown, so the search runs, and an index that
    is empty although its log lists indels still fails loudly.
    """
    for name in (f"log{chrom}.txt.gz", f"log{chrom}.txt"):
        path = os.path.join(log_dir, name)
        if not os.path.isfile(path):
            continue
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for i, line in enumerate(fh):
                if not line.strip():
                    continue
                if i == 0 and line.split("\t", 1)[0] == "CHR":
                    continue
                return False
        return True
    return False


def _dataset_chroms(vcf_dir: str, log_indels_dir: str):
    """Chromosome list for a variant dataset.

    Prefer the raw VCFs (chrom = first data record of each ``<chr>.vcf.gz``). When
    the raw VCFs are absent -- e.g. a batteries-included install that shipped only
    the precomputed index + dictionaries, not the multi-GB source VCFs -- fall back
    to the per-chromosome indel logs bundled with the index
    (``Dictionaries/log_indels_<vcf>/log<chrom>.txt``). Both are produced from the
    same VCFs at index-build time, so the chromosome set is identical.
    """
    if vcf_dir and os.path.isdir(vcf_dir):
        vcfs = sorted(f for f in os.listdir(vcf_dir) if f.endswith("vcf.gz"))
        if vcfs:
            return [_normalize_chrom(_chrom_from_vcf(os.path.join(vcf_dir, f))) for f in vcfs]
    if log_indels_dir and os.path.isdir(log_indels_dir):
        chroms = []
        for f in sorted(os.listdir(log_indels_dir)):
            if not f.startswith("log"):
                continue
            # accept gzipped logs (log<chrom>.txt.gz) as well as plain (log<chrom>.txt)
            stem = f[:-3] if f.endswith(".gz") else f
            if not stem.endswith(".txt"):
                continue
            chroms.append(_normalize_chrom(stem[len("log"):-len(".txt")]))
        if chroms:
            return chroms
    raise ValueError(
        f"No VCFs found in '{vcf_dir}' and no indel logs under '{log_indels_dir}'; "
        "cannot determine chromosomes for the indel search."
    )


ref_folder = sys.argv[1]
ref_name = os.path.basename(sys.argv[1])
vcf_dir = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
guide_name = os.path.basename(sys.argv[4])
pam_file = sys.argv[5]
pam_name = os.path.basename(sys.argv[5])
# bMax / true_pam here are the RESOLVED variant index's N and PAM (its <PAM>_<N>
# prefix), so the INDELS index path below matches the installed/precomputed index
# even when it differs from the requested budget (e.g. a batteries index with N>bMax,
# or a pamless NNN index). bDNA/bRNA below remain the user's requested bulges.
bMax = sys.argv[6]
mm = sys.argv[7]
bDNA = sys.argv[8]
bRNA = sys.argv[9]
output_folder = sys.argv[10]
true_pam = sys.argv[11]
current_working_directory = sys.argv[12]
threads = int(sys.argv[13])
# optional total-edits cap (mm+bulges), forwarded to the in-search TST prune (#107)
max_edits = sys.argv[14] if len(sys.argv) > 14 else "-1"
# optional --prune-pam pass-through ("1"=on): filter the PAM-less INDELS index down to
# pam_file at CRISPRitz leaf-emit, so the indel OTs match the SNP OTs' pruned PAM. OFF
# -> no flag -> byte-identical indel search. Needs a CRISPRitz build with --prune-pam.
prune_pam = sys.argv[15] if len(sys.argv) > 15 else "0"
prune_pam_flag = "--prune-pam" if prune_pam == "1" else ""


def search_indels(chrom):
    """Run the per-chromosome indel off-target search for one fake-indel chromosome.

    The indel search ALWAYS routes through the SHIPPED indexed indel genome
    (``genome_library/<PAM>_<N>_<ref>+<vcf>_INDELS/<PAM>_<N>_fake<chrom>/``) via
    ``crispritz.py search -index`` -- including the 0-bulge case. ``crispritz -index``
    fully supports ``-bDNA 0 -bRNA 0`` (verified: exit 0, empty stderr, well-formed
    targets), so the same indexed command serves any bulge depth. This is the correct,
    non-hack fix for issue: on a batteries / dict-less install ONLY the indexed
    ``_INDELS`` genome is shipped -- the raw ``Genomes/<ref>+<vcf>_INDELS/`` fake genome
    is absent -- so the old 0-bulge branch searched a missing directory, produced no
    ``fake<chrom>...targets.txt``, and post_analisi_indel.sh then choked on
    tail/head/rm/sed of the missing file.

    The raw-genome brute-force command is kept ONLY as a fallback for a non-batteries
    install that has the raw fake-indel genome but no indexed one (the indexed genome is
    preferred whenever present, which it always is on batteries/dict-less installs; the
    calling shell guarantees the indexed ``_INDELS`` dir exists or errors out before this
    script runs).

    Returns True when the chromosome was skipped because it has no indels (see
    ``_logged_no_indels``); it then writes no part file, and the merge knows why.
    """
    if _logged_no_indels(log_indels_dir, chrom):
        print(f"No indels on {chrom}: skipping its indel search (no indel off-target possible)", flush=True)
        return True
    print("Searching for INDELs in", chrom)
    indexed_fakechrom = os.path.join(
        current_working_directory,
        "genome_library",
        f"{true_pam}_{bMax}_{ref_name}+{vcf_name}_INDELS",
        f"{true_pam}_{bMax}_fake{chrom}",
    )
    raw_fakechrom = os.path.join(
        current_working_directory,
        "Genomes",
        f"{ref_name}+{vcf_name}_INDELS",
        f"fake_{vcf_name}_{chrom}",
    )
    out_name = f"fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}"
    if os.path.isdir(indexed_fakechrom) or not os.path.isdir(raw_fakechrom):
        # Preferred / batteries path: indexed search (valid for 0 bulges too).
        os.system(
            f"crispritz.py search {indexed_fakechrom}/ {pam_file} {guide_file} "
            f"{out_name} -index -mm {mm} -bDNA {bDNA} -bRNA {bRNA} -t -th 1 "
            f"--max-edits {max_edits} {prune_pam_flag} >/dev/null"
        )
    else:
        # Fallback: only a raw fake-indel genome is installed (no indexed one). This
        # path cannot serve bulges, so it is used strictly for a 0-bulge search.
        os.system(
            f"crispritz.py search {raw_fakechrom}/ {pam_file} {guide_file} "
            f"{out_name} -mm {mm} -t -th 1 >/dev/null"
        )
    print("Search ended for INDELs in", chrom)
    return False


log_indels_dir = os.path.join(current_working_directory, "Dictionaries", "log_indels_" + vcf_name)
chrs = _dataset_chroms(vcf_dir, log_indels_dir)

# cpus = len(os.sched_getaffinity(0))
# if cpus - 3 < 10:
#     if cpus - 3 < 0:
#         t = 1
#     else:
#         t = cpus - 3
# else:
#     t = 10

os.chdir(output_folder)
# with Pool(processes=t) as pool:
with Pool(processes=threads) as pool:
    skipped = {c for c, was_skipped in zip(chrs, pool.map(search_indels, chrs)) if was_skipped}


header = ""
merged = f"{output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
for chrom in chrs:
    # Only a chromosome skipped for having no indels may lack its part file. Any
    # other missing part file still reaches tail/head/rm below, whose stderr fails
    # the stage -- a search that produced nothing must not pass as an empty one.
    if chrom in skipped:
        continue
    os.system(
        f"tail -n +2 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt >> {output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    )
    header = os.popen(
        f"head -1 {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    ).read()
    os.system(
        f"rm {output_folder}/fake{chrom}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
    )

if len(skipped) == len(chrs):
    # no chromosome has an indel: there is nothing to merge, and no part file to
    # take the header from
    with open(merged, "w") as fh:
        fh.write(_TARGETS_HEADER + "\n")
else:
    os.system(f'sed -i 1i"{header}" {merged}')


# os.system('echo "Search INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')

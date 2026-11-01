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

# The record of a chromosome's sub-index that buildTST built with no target, which
# pool_search_indels.py reads. Its fake-indel genome holds indels, but no window of
# one holds a site of the PAM (a rare PAM on a chromosome with few indels, such as
# Cas12a's TTTV): buildTST writes no .bin, rightly, and CRISPRitz 2.9.1 refuses to
# search an index directory without one, as it does a broken or partial build. This
# file tells the two apart. index-genome removes the directory it rebuilds, record
# and all, and crispritz.py search passes over a dot file.
NO_TARGETS = ".no_targets"


def _fake_genome_is_empty(fake_dir):
    """True when every FASTA in this chromosome's fake-indel genome is header-only.

    The enricher writes one fake contig per chromosome from the windows of the
    indels some sample carries, so a chromosome with none gets a FASTA with no
    sequence. There is nothing to index: buildTST writes no .bin for it (before
    CRISPRitz 2.9.1 it segfaulted, which the exit-status check 2.9.1 adds would
    report as a failed build). The search skips the same
    chromosome on its indel log (pool_search_indels._logged_no_indels). No FASTA at
    all is not "empty": it is left to fail as before.
    """
    try:
        fastas = [f for f in os.listdir(fake_dir) if f.endswith((".fa", ".fasta"))]
    except OSError:
        return False
    if not fastas:
        return False
    for f in fastas:
        with open(os.path.join(fake_dir, f)) as fh:
            for line in fh:
                if line.strip() and not line.startswith(">"):
                    return False
    return True


def _index_pam(pam_file):
    """The PAM index-genome names an index by, read from ``pam_file`` as crispritz.py
    reads it: the first |n| bases of its sequence for a 5' PAM (n < 0), the last n
    for a 3' one."""
    with open(pam_file) as fh:
        seq, size = fh.read().split()[:2]
    size = int(size)
    return seq[:-size] if size < 0 else seq[-size:]


def _record_no_targets(chrom, fake_dir):
    """After index-genome has built ``chrom``'s sub-index: when it holds no .bin and
    no contig of the fake-indel genome holds a target, by buildTST's own rule (the
    crispritz_bins check build-index-only makes), write NO_TARGETS into it. Nothing is
    written for a sub-index with a .bin, nor for one whose genome holds a target (the
    engine left out its .bin, as with the disk's inode table full: build-index-only
    fails on it) or is read otherwise than buildTST reads it: a search refuses such a
    sub-index. Returns the status of this step: 1 when the record could not be
    written."""
    try:
        import crispritz_bins

        # where index-genome writes it: its outer directory takes the PAM of the PAM
        # file, the sub-index the name given
        sub = os.path.join(
            "genome_library",
            f"{_index_pam(pam_file)}_{bMax}_{ref_name}+{vcf_name}_INDELS",
            f"{true_pam}_{bMax}_fake{chrom}",
        )
        if any(f.endswith(".bin") for f in os.listdir(sub)):
            return 0
        # (index-genome indexes each .fa there, and none is not "no target")
        if not any(f.strip().split(".")[-1] == "fa" for f in os.listdir(fake_dir)):
            return 0
        missing, unread = crispritz_bins.missing_bins(sub, fake_dir, pam_file, int(bMax))
    except (ImportError, OSError, ValueError, IndexError):
        return 0
    if missing or unread:
        return 0
    record = os.path.join(sub, NO_TARGETS)
    try:
        with open(record, "w") as fh:
            fh.write(
                f"built, 0 targets: no window of {os.path.basename(os.path.normpath(fake_dir))} "
                f"holds a target of {os.path.basename(pam_file)} with -bMax {bMax}, so "
                f"buildTST wrote no .bin\n"
            )
    except OSError as exc:
        print(f"Indexing INDELs in {chrom}: could not write {record} ({exc})", flush=True)
        return 1
    print(f"No {true_pam} target in the indels of {chrom}: its index holds no .bin "
          f"({record} records it)", flush=True)
    return 0


def index_indels(chrom):
    fake_dir = f"{indels_folder}/fake_{vcf_name}_{chrom}"
    if _fake_genome_is_empty(fake_dir):
        print(f"No indels on {chrom}: nothing to index", flush=True)
        return
    print("Indexing INDELs in", chrom)
    status = os.system(
        f"crispritz.py index-genome {ref_name}+{vcf_name}_INDELS/{true_pam}_{bMax}_fake{chrom} {indels_folder}/fake_{vcf_name}_{chrom} {pam_file} -bMax {bMax} -th 1  >/dev/null"
    )  # {indels_folder}/fake_{vcf_name}_{chrom}
    if not status:
        status = _record_no_targets(chrom, fake_dir)
    print("Indexing ended for INDELs in", chrom)
    return status


chrs = []
for f in os.listdir(indels_folder):
    if "chr" in f:
        chrs.append(f.split("_")[-1])

with Pool(processes=int(ncpus)) as pool:
    statuses = pool.map(index_indels, chrs)

# A chromosome whose index-genome failed (a buildTST error, an OOM kill) leaves its
# part of the _INDELS index absent or partial: say so by the exit status, which
# build-index-only checks before it puts the index in place. (The search goes by
# the engine's error on stderr instead; this line goes to stdout.)
failed = sorted(c for c, s in zip(chrs, statuses) if s)
if failed:
    print("Indexing failed for INDELs in", ", ".join(failed), flush=True)
    sys.exit(1)


# os.system('echo "Indexing INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')

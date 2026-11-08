#!/usr/bin/env python
"""For a search that enriches a dataset or builds an INDELS index itself: what
build-index-only does for SNP+indel targets.

The indel search runs on a fake-indel genome that add-variants cuts from the plain
reference. Only the SNP overlay (overlay_indel_snps) lets it find a target that needs both
an indel and a SNP beside it, and the phased indel genotype store (build_indel_genotypes)
tells whose haplotype holds both. build-index-only has done both since 2.5.0. A search
that enriched the dataset itself (a download plus its VCFs asked for another PAM or bulge
budget, the web's first search of a dataset that has only its VCFs) did neither: it lost
those targets with no message, and every later run kept the INDELS index it had built.

Usage (from submit_job_automated_new_multiple_vcfs.sh; stderr fails the stage, so progress
goes to stdout, and only a failure writes to stderr):

    search_indel_snp.py overlay <fake_indel_genome_dir> <enriched_genome_dir> <log_indels_dir>
    search_indel_snp.py store <vcf_dir> <store_dir> <scratch_dir> <threads>
"""

import glob
import gzip
import multiprocessing
import os
import shutil
import sys

import build_indel_genotypes
import overlay_indel_snps


def overlay(indels_dir, enriched_dir, log_dir):
    """Put the SNP codes on a fake-indel genome. A FASTA none of whose bases changes, as
    one overlaid before, is left as it is. Run before every INDELS build, so a search of
    another PAM building from the same genome meanwhile needs no lock: each codes what is
    not coded yet, each FASTA replaced whole, and builds only once all of them are."""
    overlay_indel_snps.main(["overlay_indel_snps", indels_dir, enriched_dir, log_dir],
                            before_write=lambda: None)


def _chrom_of(vcf):
    """The chromosome of a VCF, from its first record, with its chr prefix; None for a VCF
    with none."""
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                chrom = line.split("\t", 1)[0]
                return chrom if chrom.startswith("chr") else "chr" + chrom
    return None


def store(vcf_dir, store_dir, scratch_dir, threads):
    """Write into ``store_dir`` the phased indel genotype store of each VCF of ``vcf_dir``
    that lacks a complete one written under the FILTER rule of 2.10.0
    (build_indel_genotypes.store_is_current), as build-index-only does. Written in
    ``scratch_dir`` and moved in once whole, so a run stopped meanwhile leaves no part of
    it. Every sample of
    the VCF is kept, as the indel log lists them. Guarded as build-index-only's is: a
    chromosome that fails is skipped with a WARNING, as the store is optional (its
    carriers are then PUTATIVE), and only a full disk or quota, an I/O error or no memory
    left ends the run."""
    jobs = []
    for vcf in sorted(glob.glob(os.path.join(vcf_dir, "*.vcf.gz"))
                      + glob.glob(os.path.join(vcf_dir, "*.vcf"))):
        chrom = _chrom_of(vcf)
        if chrom is None:
            continue
        name = f"gt_indel_{chrom}.tsv.gz"
        if not build_indel_genotypes.store_is_current(os.path.join(store_dir, name)):
            jobs.append((vcf, os.path.join(scratch_dir, name), None))
    if not jobs:
        return
    os.makedirs(scratch_dir, exist_ok=True)
    try:
        workers = int(os.environ.get("CRISPRME_INDELGT_WORKERS", "8"))
    except ValueError:
        workers = 8
    workers = max(1, min(threads, workers, len(jobs)))
    if workers > 1:
        with multiprocessing.get_context("fork").Pool(workers) as pool:
            pool.starmap(build_indel_genotypes.compile_indel_genotypes_safe, jobs)
    else:
        for job in jobs:
            build_indel_genotypes.compile_indel_genotypes_safe(*job)
    os.makedirs(store_dir, exist_ok=True)
    for name in sorted(os.listdir(scratch_dir)):  # (the marker of each after its store)
        shutil.move(os.path.join(scratch_dir, name), os.path.join(store_dir, name))
    print(f"[indel-snp] phased indel genotype store: {len(jobs)} chromosome(s) -> "
          f"{os.path.basename(store_dir)}", flush=True)


def main(argv):
    if len(argv) == 5 and argv[1] == "overlay":
        try:
            overlay(*argv[2:5])
        except OSError as err:  # a read-only genome folder, a full disk: one line, no traceback
            sys.stderr.write(
                f"ERROR: the SNP codes cannot be put on the fake-indel genome in {argv[2]} ({err}), "
                "and an INDELS index built from it would miss every target that needs a SNP beside "
                "an indel. Make that folder writable, or install an index that holds its INDELS "
                "index\n")
            return 1
    elif len(argv) == 6 and argv[1] == "store":
        store(argv[2], argv[3], argv[4], int(argv[5]))
    else:
        sys.stderr.write(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

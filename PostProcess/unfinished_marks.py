"""The parts of a variant dataset that a search refuses while build-index-only marks them
unfinished.

build-index-only writes some of a dataset's outputs in place, chromosome by chromosome or
FASTA by FASTA, and marks each meanwhile with a hidden ``.building`` directory that holds
its owner record (crisprme.py _mark_unfinished): a run stopped there leaves a part that
looks whole, which a search would use as if it were. complete-search refuses each part
while it is marked (crisprme.py _refuse_unfinished), and so does complete-test's
one-chromosome smoke, which runs complete-search from a view of the working directory
that the marks are not in, so it looks at them in the working directory first. Both take
the parts from ``dataset_parts``, so that they refuse the same ones. (The search script
the web runs has its own copy of these checks, in bash: STEP 1 and STEP 2.3 of
submit_job_automated_new_multiple_vcfs.sh.)
"""

import os

import index_fit


def indel_snp_on() -> bool:
    """Whether the SNP+indel co-occurrence is on (CRISPRME_INDEL_SNP, on unless set to
    anything but 1, true, True or yes): the indel genotype store is read for it alone."""
    return os.environ.get("CRISPRME_INDEL_SNP", "1") in ("1", "true", "True", "yes")


def is_marked(directory: str) -> bool:
    """Whether ``directory`` holds build-index-only's unfinished mark."""
    return os.path.isdir(os.path.join(directory, ".building"))


def _pam_serves(pam: str, motif: str, prune: bool, three_prime: bool) -> bool:
    """Whether an index built for PAM ``pam`` serves a search for PAM ``motif``: the
    motif itself, or an all-N PAM of its width, or -- under --prune-pam -- a dedicated
    PAM that holds it, or, for a 3' PAM only, a wider all-N one (index_fit.pam_serves:
    the search script's _resolve_index_dir, and crisprme.py's _installed_index_bulge_cap)."""
    return index_fit.pam_serves(pam, motif, prune, three_prime) is not None


def indels_index_in_place(library: str, dataset: str, motif: str, n_want: int,
                          prune: bool = False, three_prime: bool = True) -> bool:
    """Whether ``library`` (a genome_library) holds a variant index of ``dataset``
    (``<genome>+<vcf>``) that serves a search for PAM ``motif`` with index bulge budget
    ``n_want`` (max(bDNA, bRNA) + 1), with its INDELS index beside it, neither of them
    locked by a run that builds it, as the search script takes them (STEP 2.2, 2.3).
    Such a search does not read the dataset's fake-indel genome, which the script reads
    only to build an INDELS index it does not find (STEP 2.3). With several indexes
    that serve, the script searches one of them (_resolve_index_dir), whose INDELS index
    may be missing: the script then refuses a marked fake-indel genome at STEP 2.3."""
    suffix = "_" + dataset
    try:
        names = os.listdir(library)
    except OSError:
        return False
    for name in names:
        if not name.endswith(suffix) or name.startswith("."):
            continue
        pam, sep, n_str = name[: -len(suffix)].partition("_")
        if not sep or not n_str.isdigit() or int(n_str) < n_want:
            continue
        if not _pam_serves(pam, motif, prune, three_prime):
            continue
        index = os.path.join(library, name)
        indels = index + "_INDELS"
        if not (os.path.isdir(index) and os.path.isdir(indels)) or any(
                os.path.isdir(d + ".lock") for d in (index, indels)):
            continue
        try:
            if any(f.endswith(".bin") for f in os.listdir(index)):
                return True
        except OSError:
            continue
    return False


def _real(path: str) -> str:
    """``path``, or what it links to when it is a symlink: complete-test's view links
    each Dictionaries entry, as a user may link one, and the path to remove is the
    real one."""
    return os.path.realpath(path) if os.path.islink(path) else path


def dataset_parts(workdir: str, genome: str, vcf: str, indels_in_place: bool) -> list:
    """The parts of variant dataset ``vcf`` of genome ``genome``, in working directory
    ``workdir``, that a search refuses while build-index-only marks them unfinished, as
    ``(directory, what, rerun)`` in the order complete-search checks them: what
    build-index-only writes there, and what to do once that run has stopped. The indel
    genotype store only with the SNP+indel co-occurrence on (indel_snp_on), which alone
    reads it; the fake-indel genome only for a search without its INDELS index in place
    (``indels_in_place``: indels_index_in_place), which alone reads it."""
    dicts = os.path.join(workdir, "Dictionaries")
    reg = _real(os.path.join(dicts, f"registry_{vcf}"))
    gt = _real(os.path.join(dicts, f"genotypes_{vcf}"))
    rerun = f"Run build-index-only again with {vcf}"
    parts = [(reg, f"the dict-less tiers of {vcf}",
              f"{rerun} and its --samplesID, which finishes them, or remove them to "
              f"search without them (rm -rf '{reg}' '{gt}')")]
    if indel_snp_on():
        igt = _real(os.path.join(dicts, f"indel_genotypes_{vcf}"))
        parts.append((igt, f"the phased indel genotype store of {vcf}",
                      f"{rerun}, which finishes it, or remove it to search without it "
                      f"(rm -rf '{igt}')"))
    enriched = os.path.join(workdir, "Genomes", f"{genome}+{vcf}")
    parts.append((_real(enriched), f"the enrichment of {vcf}", f"{rerun}, which enriches it again"))
    if not indels_in_place:
        parts.append((_real(enriched + "_INDELS"), f"the SNP codes of the {vcf} fake-indel genome",
                      f"{rerun}, which finishes them"))
    return parts

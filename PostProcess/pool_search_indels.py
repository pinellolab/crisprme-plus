#!/usr/bin/env python

import gzip
import os
import sys
import tempfile
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


def _checked_chrom(chrom: str) -> str:
    """``chrom`` back, unless its name holds a character a shell reads. Each name here
    comes from a dataset VCF's CHROM column (or a bundled indel log named after it) and
    is taken, unquoted, into the search and merge commands this script runs through a
    shell: a broken or hostile dataset whose CHROM holds such a character would run it. A
    real chromosome name holds none; refused as a genome's sequence id is (genome_layout).
    See the twin in ``pool_post_analisi_indel.py``."""
    from genome_layout import PATTERN_CHARS, UNUSABLE_ID_START, UNUSABLE_IN_ID

    bad = next((c for c in chrom if c in UNUSABLE_IN_ID + PATTERN_CHARS or c.isspace()), None)
    if bad is not None or chrom[:1] in UNUSABLE_ID_START:
        raise ValueError(
            f"chromosome name {chrom!r} holds {(bad if bad else chrom[:1])!r}, which "
            "CRISPRme cannot take into a shell command: fix the dataset's CHROM column"
        )
    return chrom


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
    is empty although its log lists indels still fails loudly, unless it is shown to
    hold no target (``_empty_index``), as does a bulge search of one that is missing
    (``_missing_index``).
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
            return [_checked_chrom(_normalize_chrom(_chrom_from_vcf(os.path.join(vcf_dir, f)))) for f in vcfs]
    if log_indels_dir and os.path.isdir(log_indels_dir):
        chroms = []
        for f in sorted(os.listdir(log_indels_dir)):
            if not f.startswith("log"):
                continue
            # accept gzipped logs (log<chrom>.txt.gz) as well as plain (log<chrom>.txt)
            stem = f[:-3] if f.endswith(".gz") else f
            if not stem.endswith(".txt"):
                continue
            chroms.append(_checked_chrom(_normalize_chrom(stem[len("log"):-len(".txt")])))
        if chroms:
            return chroms
    raise ValueError(
        f"No VCFs found in '{vcf_dir}' and no indel logs under '{log_indels_dir}'; "
        "cannot determine chromosomes for the indel search."
    )


# What search_indels did with a chromosome.
SEARCHED, SKIPPED, REFUSED = "searched", "skipped", "refused"

# The record pool_index_indels.py writes into a chromosome's sub-index that buildTST
# built with no target: its indels are there, but no window of one holds a site of
# the PAM, so there is no .bin.
NO_TARGETS = ".no_targets"

# The prefix of the folder complete-test --chrom runs its search from
# (complete_test.CHROMVIEWPREFIX): see _outside_view.
CHROMVIEWPREFIX = ".crisprme-test-view-"


def _holds_bin(index_dir: str) -> bool:
    """Whether ``index_dir`` holds a tree, a regular file ``<name>.bin``: what
    CRISPRitz 2.9.1's searchTST loads, which refuses a directory with none. True
    where it cannot be listed: the engine then says what is wrong, as before."""
    try:
        names = os.listdir(index_dir)
    except OSError:
        return True
    return any(
        len(f) > 4 and f.endswith(".bin") and os.path.isfile(os.path.join(index_dir, f))
        for f in names
    )


def _pam_line(pam_path: str):
    """The sequence and the number of a PAM file, as crispritz.py reads them."""
    with open(pam_path) as fh:
        seq, size = fh.read().split()[:2]
    return seq, int(size)


def _fake_genome(fake_dir: str):
    """What a chromosome's fake-indel genome holds, by pool_index_indels.py's test of it
    (_fake_genome_is_empty): "sequence", "empty" (every FASTA header-only: no indel,
    and no sub-index built for it), or None (no FASTA there, or one not read)."""
    try:
        fastas = [f for f in os.listdir(fake_dir) if f.endswith((".fa", ".fasta"))]
        for f in fastas:
            with open(os.path.join(fake_dir, f)) as fh:
                if any(line.strip() and not line.startswith(">") for line in fh):
                    return "sequence"
    except (OSError, UnicodeDecodeError):
        return None
    return "empty" if fastas else None


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
# optional: the library that holds the variant index the search resolved, and its INDELS
# index beside it -- a library given with --index-path, or the working directory's
index_library = (sys.argv[16] if len(sys.argv) > 16
                 else os.path.join(current_working_directory, "genome_library"))
# optional: that variant index's name, <true_pam>_<bMax>_<ref>+<vcf> or, since 2.9.0, the
# form that carries the guide length and the side (<guide>-<PAM>_... for a 3' PAM,
# <PAM>-<guide>_... for a 5' one: index_fit.parse_index_name); its INDELS index is
# <name>_INDELS. The sub-indexes in it keep the engine's names, <true_pam>_<bMax>_fake<chrom>.
variant_index = (sys.argv[17] if len(sys.argv) > 17
                 else f"{true_pam}_{bMax}_{ref_name}+{vcf_name}")


def _in_chrom_view():
    """Whether this search runs from complete-test --chrom's view of the working directory
    (complete_test._chrom_view), a folder in it that the run removes when it ends."""
    folder = os.path.basename(os.path.normpath(current_working_directory))
    return folder.startswith(CHROMVIEWPREFIX)


def _outside_view(path):
    """``path`` as the user knows it: where it lies in complete-test --chrom's view, the
    same path in the working directory that view mirrors. The search reads the view,
    which holds this chromosome alone. But what a message names must outlast the run,
    and what tells of an index is read in the working directory too: the view leaves out
    an index's .pam_build, and need not hold the fake-indel genome."""
    if not _in_chrom_view():
        return path
    view = os.path.normpath(current_working_directory)
    inside = os.path.relpath(path, view)
    if inside == os.pardir or inside.startswith(os.pardir + os.sep):
        return path
    return os.path.normpath(os.path.join(os.path.dirname(view), inside))


def _own_pam():
    """``(line, name)`` when the PAM of the search's PAM file is not ``true_pam``, the PAM
    of the index: the PAM line (``"<sequence> <n>"``) of the PAM file the index was built
    with, and that file's name (None where the index does not record it). Else None.

    The index is then a PAM-less one standing in for the search's PAM (NNN_3 serving
    NRG, or NNNNNNNN_3 under --prune-pam, as the search resolves one), whose INDELS index
    was built with its own PAM (pool_index_indels._own_pam): build-index-only records the
    PAM file of a build in the variant index's .pam_build (``<name> <length> <n>``);
    without that record, the line is the search PAM file's guide around ``true_pam``."""
    try:
        seq, size = _pam_line(pam_file)
    except (OSError, ValueError):
        return None
    if (seq[:-size] if size < 0 else seq[-size:]) == true_pam:
        return None
    k = len(true_pam)
    guide, n, name = len(seq) - abs(size), (k if size > 0 else -k), None
    built = os.path.join(_outside_view(index_library), variant_index, ".pam_build")
    try:
        with open(built) as fh:
            got = fh.read().rsplit(None, 2)
        if len(got) == 3 and abs(int(got[2])) == k and int(got[1]) > k:
            name, guide, n = got[0], int(got[1]) - k, int(got[2])
    except (OSError, ValueError):
        pass
    return ("N" * guide + true_pam if n > 0 else true_pam + "N" * guide) + f" {n}", name


def _pam_dirs():
    """Where a PAM file is looked for: the working directory's PAMs/, and the folder of
    the search's."""
    dirs = []
    for d in (os.path.join(_outside_view(current_working_directory), "PAMs"),
              os.path.dirname(os.path.abspath(pam_file))):
        d = os.path.normpath(d)
        if d not in dirs:
            dirs.append(d)
    return dirs


def _pam_file_of(line, name):
    """A PAM file whose first line is ``line``, in _pam_dirs: the one named ``name``
    first, then any other there; None when there is none."""
    want = line.split()
    for named in (True, False):
        for d in _pam_dirs():
            try:
                names = sorted(f for f in os.listdir(d) if f.endswith(".txt"))
            except OSError:
                continue
            for f in names:
                if (f == f"{name}.txt") != named:
                    continue
                path = os.path.join(d, f)
                try:
                    with open(path) as fh:
                        if fh.readline().split() == want:
                            return path
                except (OSError, UnicodeDecodeError):
                    continue
    return None


def _holds_vcfs(vcf_folder):
    """Whether ``vcf_folder`` holds VCFs, as add-variants and _dataset_chroms read one."""
    try:
        return any(f.endswith("vcf.gz") for f in os.listdir(vcf_folder))
    except OSError:
        return False


def _rebuild(raw_fakechrom, fasta):
    """How to rebuild the INDELS index this search reads, the end of "Rebuild it", by what
    is here. With the dataset's VCFs, the build-index-only command, with the PAM file the
    index was built with, that does it, with the enrichment that makes the fake-indel
    genome too when that genome is not here (``fasta`` false), or else the download, for
    an index that came that way. Without those VCFs (an index downloaded, or copied from
    where it was built) that command cannot run, and what it removes first it could not
    make again: with the fake-indel genome here, the search builds a missing INDELS index
    from it with the index's own PAM, so removing it is enough; without it, there is the
    download. So too where no such PAM file is to be found (a PAM-less index searched for
    another PAM, built with one that is not in _pam_dirs): the steps are said in words,
    never as a command with a blank to fill in. complete-test --chrom's view may not hold
    that genome, and an INDELS index built there would go with the view, so that run sends
    the user to a run that does build it."""
    cwd = os.path.normpath(_outside_view(current_working_directory))
    indels_index = os.path.join(os.path.normpath(_outside_view(index_library)),
                                f"{variant_index}_INDELS")
    # (a downloaded index comes without the fake-indel genome. download skips an index
    # whose INDELS index is in place, and else swaps in both from one archive: a variant
    # index left there survives a download that fails, as for an index that is not on
    # the hub)
    download = (f"rm -rf '{indels_index}' && crisprme.py download --what index --index-name "
                f"{os.path.basename(indels_index)[: -len('_INDELS')]} --path '{cwd}'")
    # (where the search's --vcf list spells the dataset as VCFs/<name>, the search looks in
    # VCFs/VCFs/<name>, and takes the chromosomes from the indel logs)
    named = os.path.join(cwd, "VCFs", vcf_name)
    vcfs = next((d for d in (_outside_view(vcf_dir), named) if _holds_vcfs(d)), None)
    own = _own_pam()
    pam = _pam_file_of(*own) if own else pam_file
    if fasta and not (vcfs and pam):
        search = ("complete-test without --chrom, or a complete-search of this index,"
                  if _in_chrom_view() else "this search again,")
        return (f": remove the INDELS index with rm -rf '{indels_index}' and run {search} which "
                f"builds it anew from the fake-indel genome in {os.path.dirname(raw_fakechrom)}")
    if not vcfs:
        return (f": for an index you downloaded, download it again: {download}; for one you "
                f"built, put its VCFs back in '{named}' and run this search again, which then "
                f"gives the command that rebuilds it")
    remove = [indels_index]
    again = ""
    if not fasta:
        # build-index-only makes that genome as it enriches, and with an enriched
        # genome in place it takes the dataset for one without indels
        remove.append(os.path.join(cwd, "Genomes", f"{ref_name}+{vcf_name}"))
        again = f"; or, for an index you downloaded, download it again: {download}"
    # the bulges build-index-only builds the index's N from
    n = int(bMax) - 1 if bMax.isdigit() and int(bMax) > 0 else bMax
    rest = f"--bDNA {n} --bRNA {n} --vcf '{vcfs}' --path '{cwd}'"
    if pam:
        removed = " ".join("'%s'" % p for p in remove)
        return (f": rm -rf {removed} && crisprme.py build-index-only --genome "
                f"'{_outside_view(ref_folder)}' --pam '{pam}' {rest}{again}")
    named = f"{own[1]}.txt, " if own[1] else ""
    return (f" with build-index-only and a PAM file of the PAM it was built with ({named}the "
            f"line '{own[0]}'), which is not in {' or '.join(_pam_dirs())}: remove "
            + " and ".join("'%s'" % p for p in remove)
            + f", write that PAM file, and run crisprme.py build-index-only with it as --pam, "
              f"and --genome '{_outside_view(ref_folder)}' {rest}{again}")


def _empty_index(chrom, indexed_fakechrom, raw_fakechrom):
    """``(True, why)`` when ``chrom``'s sub-index, which holds no .bin, holds no target
    either, so no indel off-target is the answer there; else ``(False, error)``, an
    ERROR that says how to rebuild it.

    pool_index_indels.py records a build that found no target (NO_TARGETS). An index
    built before it did has no record: there the fake-indel genome it was built from
    tells, read by buildTST's own rule (crispritz_bins, as build-index-only checks a
    build): when no contig of it holds a target of the index's own PAM (the search's,
    or the one a PAM-less index standing in for it was built with: _own_pam), a build
    holds no .bin. Anything else may be a broken or partial build, and the search must
    not pass it off as an empty result.
    """
    if os.path.isfile(os.path.join(indexed_fakechrom, NO_TARGETS)):
        return True, (f"No {true_pam} target in the indels of {chrom}: skipping its indel "
                      f"search (its index was built with none)")
    fasta = False
    try:
        import crispritz_bins

        # (the FASTAs index-genome would index there)
        fasta = os.path.isdir(raw_fakechrom) and any(
            f.strip().split(".")[-1] == "fa" for f in os.listdir(raw_fakechrom))
        if not fasta:
            why = (f"the fake-indel genome it was built from ({raw_fakechrom}) is not here to "
                   f"check it against")
        else:
            own = _own_pam()
            with tempfile.TemporaryDirectory(prefix="crisprme_indels_pam_") as tmp:
                built_with = pam_file
                if own:
                    built_with = os.path.join(tmp, f"{true_pam}.txt")
                    with open(built_with, "w") as fh:
                        fh.write(own[0] + "\n")
                missing, unread = crispritz_bins.missing_bins(
                    indexed_fakechrom, raw_fakechrom, built_with, int(bMax))
            if not missing and not unread:
                return True, (f"No {true_pam} target in the indels of {chrom}: skipping its "
                              f"indel search (its index holds no .bin, and {raw_fakechrom} "
                              f"no target)")
            why = (f"{os.path.join(raw_fakechrom, missing[0][0])} holds a target of "
                   f"{true_pam}: the index is incomplete" if missing else
                   f"{os.path.join(raw_fakechrom, unread[0])} cannot be read as buildTST reads it")
    except (ImportError, OSError, ValueError, IndexError) as exc:
        why = f"it could not be checked ({exc})"
    return False, (f"ERROR: the indel index of {chrom} ({_outside_view(indexed_fakechrom)}) "
                   f"holds no .bin and no record that its build found no target there (an "
                   f"index built before CRISPRme 2.8.1 has none), and {why}. Rebuild it"
                   f"{_rebuild(raw_fakechrom, fasta)}")


def _bulged():
    """Whether this search allows a bulge."""
    try:
        return int(bDNA) > 0 or int(bRNA) > 0
    except ValueError:
        return True


def _missing_index(chrom, indexed_fakechrom, raw_fakechrom):
    """For a search with bulges: the ERROR for ``chrom``, which has indels (its log lists
    them, or its fake-indel genome holds them) but no sub-index in the INDELS index,
    saying how to rebuild the index; None where none is due.

    Without its sub-index, a chromosome's search falls back to its fake-indel genome,
    which it searches with no bulge: the bulged indel off-targets there were lost with
    no error, after an INDELS index was copied in part, or after a search had built a
    missing one for a PAM-less index under the search's PAM (as pool_index_indels.py
    did before _own_pam). A sub-index is rightly absent where that genome is
    header-only (pool_index_indels.py indexes none), and for a contig with an underscore
    in its name, which no build gives one (pool_index_indels.py says so as it builds),
    so that no rebuild would help: its search goes on, as on every release before this
    check. Where nothing shows that the chromosome has indels, the search goes on as
    before (the engine then says what is wrong)."""
    if "_" in chrom:
        return None
    held = _fake_genome(raw_fakechrom)
    if held == "empty":
        return None
    logged = any(os.path.isfile(os.path.join(log_indels_dir, f"log{chrom}.txt{gz}"))
                 for gz in (".gz", ""))
    if held is None and not logged:
        return None
    try:  # (the FASTAs index-genome would index there)
        fasta = any(f.strip().split(".")[-1] == "fa" for f in os.listdir(raw_fakechrom))
    except OSError:
        fasta = False
    shown = "its indel log lists them" if logged else f"{raw_fakechrom} holds them"
    return (f"ERROR: the indel index of {chrom} ({_outside_view(indexed_fakechrom)}) is not "
            f"there, though {chrom} has indels ({shown}): the INDELS index is incomplete, and "
            f"a search with bulges cannot do without it. Rebuild it"
            f"{_rebuild(raw_fakechrom, fasta)}")


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

    The raw-genome brute-force command is kept ONLY as a fallback for a 0-bulge search
    on a non-batteries install that has the raw fake-indel genome but no indexed one
    (the indexed genome is preferred whenever present, which it always is on
    batteries/dict-less installs; the calling shell guarantees the indexed ``_INDELS``
    dir exists or errors out before this script runs). It cannot search with bulges, so
    a search with bulges whose chromosome has indels but no sub-index stops there
    instead (``_missing_index``).

    Returns SKIPPED when the chromosome has no indels (see ``_logged_no_indels``),
    or its sub-index no .bin and no target (``_empty_index``): it then writes no part
    file, and the merge knows why. REFUSED when that sub-index may be broken, or is
    missing: it has said why on stderr, and nothing is merged.
    """
    if _logged_no_indels(log_indels_dir, chrom):
        print(f"No indels on {chrom}: skipping its indel search (no indel off-target possible)", flush=True)
        return SKIPPED
    indexed_fakechrom = os.path.join(
        index_library,
        f"{variant_index}_INDELS",
        f"{true_pam}_{bMax}_fake{chrom}",
    )
    raw_fakechrom = os.path.join(
        current_working_directory,
        "Genomes",
        f"{ref_name}+{vcf_name}_INDELS",
        f"fake_{vcf_name}_{chrom}",
    )
    if os.path.isdir(indexed_fakechrom) and not _holds_bin(indexed_fakechrom):
        empty, said = _empty_index(chrom, indexed_fakechrom, _outside_view(raw_fakechrom))
        if empty:
            print(said, flush=True)
            return SKIPPED
        sys.stderr.write(said + "\n")
        sys.stderr.flush()
        return REFUSED
    if not os.path.isdir(indexed_fakechrom) and _bulged():
        said = _missing_index(chrom, indexed_fakechrom, _outside_view(raw_fakechrom))
        if said:
            sys.stderr.write(said + "\n")
            sys.stderr.flush()
            return REFUSED
    # (flushed, as every line of a worker: the pool's exit terminates the workers, and
    # what one still holds in its buffer is lost to log_verbose.txt)
    print("Searching for INDELs in", chrom, flush=True)
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
    print("Search ended for INDELs in", chrom, flush=True)
    return SEARCHED


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
    outcomes = pool.map(search_indels, chrs)
if REFUSED in outcomes:
    # (each said why on stderr, which fails the stage): no merged file, which would
    # pass for a search of every chromosome
    sys.exit(1)
skipped = {c for c, outcome in zip(chrs, outcomes) if outcome == SKIPPED}


header = ""
merged = f"{output_folder}/indels_{ref_name}+{vcf_name}_{pam_name}_{guide_name}_{mm}_{bDNA}_{bRNA}.targets.txt"
for chrom in chrs:
    # Only a chromosome skipped for having no indels, or no target, may lack its part
    # file. Any other missing part file still reaches tail/head/rm below, whose stderr
    # fails the stage -- a search that produced nothing must not pass as an empty one.
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

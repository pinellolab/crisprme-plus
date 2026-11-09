#!/usr/bin/env python

import os
import shlex
import shutil
import sys
import tempfile
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
# optional: the variant index this INDELS index goes with, whose .pam_build _own_pam reads:
# its folder, or its name in genome_library/ (<true_pam>_<bMax>_<ref>+<vcf> or, since 2.10.0,
# the form that carries the guide length and the side: index_fit.parse_index_name). A search
# passes its folder: it runs this in a folder of the INDELS index's build lock, and puts what
# it writes there in place (submit_job_automated_new_multiple_vcfs.sh _build_in_lock).
variant_index = (sys.argv[8] if len(sys.argv) > 8
                 else f"{true_pam}_{bMax}_{ref_name}+{vcf_name}")
# the dataset's indel logs, which record each fake-indel window (Dictionaries/log_indels_<vcf>
# of the working directory that holds Genomes/, where both the search and build-index-only
# keep them)
LOG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(indels_folder))),
                       "Dictionaries", f"log_indels_{vcf_name}")

# The record of a chromosome's sub-index that buildTST built with no target, which
# pool_search_indels.py reads. Its fake-indel genome holds indels, but no window of
# one holds a site of the PAM (a rare PAM on a chromosome with few indels, such as
# Cas12a's TTTV): buildTST writes no .bin, rightly, and CRISPRitz 2.9.1 refuses to
# search an index directory without one, as it does a broken or partial build. This
# file tells the two apart. index-genome removes the directory it rebuilds, record
# and all, and crispritz.py search passes over a dot file. (Each sub-index also records
# what it holds, its index format among it: index_fit.write_index_meta, _record_index.)
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


def _pam_line(pam_file):
    """The sequence and the number of a PAM file, as crispritz.py reads them."""
    with open(pam_file) as fh:
        seq, size = fh.read().split()[:2]
    return seq, int(size)


def _index_pam(pam_file):
    """The PAM index-genome names an index by, read from ``pam_file`` as crispritz.py
    reads it: the first |n| bases of its sequence for a 5' PAM (n < 0), the last n
    for a 3' one."""
    seq, size = _pam_line(pam_file)
    return seq[:-size] if size < 0 else seq[-size:]


def _own_pam():
    """``(line, name)`` when the PAM of ``pam_file`` is not ``true_pam``, the PAM of the
    variant index this INDELS index goes with: the PAM line (``"<sequence> <n>"``) of the
    PAM file that index was built with, and that file's name (None where the index does
    not record it). Else None.

    That variant index is then a PAM-less one standing in for the search's PAM (NNN_3
    serving NRG, or NNNNNNNN_3 under --prune-pam, as complete-search resolves one), or
    under --prune-pam a dedicated one of a PAM that holds it (NRR_3 serving NGG), and
    the search builds its INDELS index, when it is missing, with the PAM file of the
    search. index-genome names an index by the PAM of its PAM file, so that build went
    into <search PAM>_<N>_<ref>+<vcf>_INDELS, where the search does not look: finding no
    sub-index in <true_pam>_<N>_<ref>+<vcf>_INDELS, it searched each chromosome's fake-indel
    genome without bulges, and lost the bulged indel off-targets with no error.
    build-index-only records the PAM file of a build in the variant index's .pam_build
    (``<name> <length> <n>``); without that record, the line is the search PAM file's
    guide around ``true_pam``."""
    try:
        seq, size = _pam_line(pam_file)
    except (OSError, ValueError):
        return None  # (index-genome says what is wrong with it)
    if (seq[:-size] if size < 0 else seq[-size:]) == true_pam:
        return None
    k = len(true_pam)
    guide, n, name = len(seq) - abs(size), (k if size > 0 else -k), None
    try:
        with open(os.path.join("genome_library", variant_index, ".pam_build")) as fh:
            built = fh.read().rsplit(None, 2)
        if len(built) == 3 and abs(int(built[2])) == k and int(built[1]) > k:
            name, guide, n = built[0], int(built[1]) - k, int(built[2])
    except (OSError, ValueError):
        pass
    return ("N" * guide + true_pam if n > 0 else true_pam + "N" * guide) + f" {n}", name


def _sub_index(chrom):
    """Where index-genome writes ``chrom``'s sub-index: its outer directory takes the PAM
    of the PAM file, the sub-index the name given."""
    return os.path.join(
        "genome_library",
        f"{_index_pam(build_pam)}_{bMax}_{ref_name}+{vcf_name}_INDELS",
        f"{true_pam}_{bMax}_fake{chrom}",
    )


def _record_no_targets(chrom, fake_dir):
    """After index-genome has built ``chrom``'s sub-index: when it holds no .bin and
    no contig of the fake-indel genome holds a target, by buildTST's own rule (the
    crispritz_bins check build-index-only makes, by the rule of the index format the
    installed CRISPRitz writes, the one that has just built it), write NO_TARGETS into
    it. Nothing is written for a sub-index with a .bin, nor for one whose genome holds a
    target (the engine left out its .bin, as with the disk's inode table full:
    build-index-only fails on it) or is read otherwise than buildTST reads it: a search
    refuses such a sub-index. Returns the status of this step: 1 when the record could
    not be written."""
    try:
        import crispritz_bins
        import index_fit

        sub = _sub_index(chrom)
        if any(f.endswith(".bin") for f in os.listdir(sub)):
            return 0
        # (index-genome indexes each .fa there, and none is not "no target")
        if not any(f.strip().split(".")[-1] == "fa" for f in os.listdir(fake_dir)):
            return 0
        missing, unread = crispritz_bins.missing_bins(sub, fake_dir, build_pam, int(bMax),
                                                      index_fit.engine_index_format())
    except (ImportError, OSError, ValueError, IndexError):
        return 0
    if missing or unread:
        return 0
    record = os.path.join(sub, NO_TARGETS)
    try:
        with open(record, "w") as fh:
            fh.write(
                f"built, 0 targets: no window of {os.path.basename(os.path.normpath(fake_dir))} "
                f"holds a target of {os.path.basename(build_pam)} with -bMax {bMax}, so "
                f"buildTST wrote no .bin\n"
            )
    except OSError as exc:
        print(f"Indexing INDELs in {chrom}: could not write {record} ({exc})", flush=True)
        return 1
    print(f"No {true_pam} target in the indels of {chrom}: its index holds no .bin "
          f"({record} records it)", flush=True)
    return 0


def _contigs_left_out():
    """The contigs with indels that get no sub-index: those with an underscore in their
    name. The chromosome the loop below takes from such a directory is what follows its
    last underscore, a fake-indel genome that is not there: index-genome says so, on the
    stdout it sends away, and exits 0, on every release. build-index-only expects no
    sub-index of them, and a search with bulges requires none. A header-only genome
    (no indels) is not left out of anything."""
    prefix = f"fake_{vcf_name}_"
    contigs = []
    for f in sorted(os.listdir(indels_folder)):
        if f.startswith(prefix) and "_" in f[len(prefix):]:
            try:
                if _fake_genome_is_empty(os.path.join(indels_folder, f)):
                    continue
            except (OSError, UnicodeDecodeError):
                pass
            contigs.append(f[len(prefix):])
    return contigs


def index_indels(chrom):
    fake_dir = f"{indels_folder}/fake_{vcf_name}_{chrom}"
    if _fake_genome_is_empty(fake_dir):
        print(f"No indels on {chrom}: nothing to index", flush=True)
        return
    # (flushed, as every line of a worker: the pool's exit terminates the workers, and
    # what one still holds in its buffer is lost)
    print("Indexing INDELs in", chrom, flush=True)
    status = os.system(
        f"crispritz.py index-genome {ref_name}+{vcf_name}_INDELS/{true_pam}_{bMax}_fake{chrom} {indels_folder}/fake_{vcf_name}_{chrom} {build_pam_arg} -bMax {bMax} -th 1  >/dev/null"
    )  # {indels_folder}/fake_{vcf_name}_{chrom}
    if not status:
        status = _record_no_targets(chrom, fake_dir)
    if not status:
        _record_index(chrom)
    print("Indexing ended for INDELs in", chrom, flush=True)
    return status


def _log_window(chrom):
    """The fake-indel windows (left, right) of ``chrom``'s indel log, which records each
    window's start and end (indel_window): the windows its fake-indel genome holds, which
    the sub-index records. None where the log is not there, holds no row, or holds two
    widths (it is then no log of this genome)."""
    import indel_window

    path = indel_window.log_path(LOG_DIR, chrom)
    try:
        found = indel_window.log_windows(path) if path else set()
    except (OSError, ValueError, EOFError):
        return None
    return next(iter(found)) if len(found) == 1 else None


def _record_index(chrom):
    """Record in ``chrom``'s sub-index, which index-genome has just built, what it holds
    (index_fit.write_index_meta): its index format among it, which a sub-index with no .bin
    records nowhere else, and the fake-indel windows it was built from, which its indel log
    tells (_log_window): a search refuses the sub-index beside a log of other windows. A
    record that cannot be written is said on stdout."""
    sub = _sub_index(chrom)
    if not os.path.isdir(sub):
        return
    try:
        import index_fit

        why = index_fit.write_index_meta(sub, build_pam, built=True,
                                         indel_window=_log_window(chrom))
    except (ImportError, OSError) as exc:
        why = str(exc)
    if why:
        print(f"Indexing INDELs in {chrom}: {why}", flush=True)


chrs = []
for f in os.listdir(indels_folder):
    if "chr" in f:
        chrs.append(f.split("_")[-1])

left_out = _contigs_left_out()
if left_out:
    print(f"No INDELS sub-index for {', '.join(left_out)}: a contig with an underscore in its "
          f"name gets none, so a search does not find the indel off-targets there that need a "
          f"bulge", flush=True)

# The PAM file index-genome builds with: the one given, but for a PAM-less variant
# index that stands in for its PAM, a file of the PAM that index was built with
# (_own_pam), so that this INDELS index is the one build-index-only builds with it,
# where the search looks for it.
build_pam = build_pam_arg = pam_file
own = _own_pam()
own_dir = None
if own:
    try:
        own_dir = tempfile.mkdtemp(prefix="crisprme_indels_pam_")
        build_pam = os.path.join(own_dir, f"{true_pam}.txt")
        with open(build_pam, "w") as fh:
            fh.write(own[0] + "\n")
    except OSError as exc:
        # (stderr: the search fails this step on it, as on an engine error)
        sys.stderr.write(f"ERROR: could not write the PAM file of {true_pam} to index the "
                         f"INDELs with ({exc})\n")
        if own_dir:
            shutil.rmtree(own_dir, ignore_errors=True)
        sys.exit(1)
    build_pam_arg = shlex.quote(build_pam)
    # (a PAM-less one, or a dedicated one of a PAM that holds the search's: NRR for NGG)
    _what = "PAM-less" if set(true_pam) == {"N"} else f"of the PAM {true_pam}"
    print(f"The variant index {os.path.basename(variant_index)} is {_what} and stands in "
          f"for {_index_pam(pam_file)}: indexing the INDELs with its own PAM ({own[0]}"
          + (f", as {own[1]}" if own[1] else "") + ")", flush=True)

try:
    with Pool(processes=int(ncpus)) as pool:
        statuses = pool.map(index_indels, chrs)
    if left_out:
        # (where those contigs alone have indels nothing has made the index's directory,
        # and build-index-only would take it for a failed build)
        os.makedirs(os.path.join("genome_library", f"{_index_pam(build_pam)}_{bMax}_{ref_name}"
                                 f"+{vcf_name}_INDELS"), exist_ok=True)
finally:
    if own_dir:
        shutil.rmtree(own_dir, ignore_errors=True)

# A chromosome whose index-genome failed (a buildTST error, an OOM kill) leaves its
# part of the _INDELS index absent or partial: say so by the exit status, which
# build-index-only checks before it puts the index in place. (The search goes by
# the engine's error on stderr instead; this line goes to stdout.)
failed = sorted(c for c, s in zip(chrs, statuses) if s)
if failed:
    print("Indexing failed for INDELs in", ", ".join(failed), flush=True)
    sys.exit(1)


# os.system('echo "Indexing INDELs End: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'" >> '+output_folder+'/../log.txt')

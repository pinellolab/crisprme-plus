"""The fake-indel window, read from the data: how many reference bases CRISPRitz's
add-variants writes on either side of each indel (and MNV) record.

CRISPRitz's enricher cuts each indel allele's contig, its fake-indel window, from the
reference: LEFT bases before the record's REF, the REF (from POS on) with the ALT in its
place, then RIGHT bases after the REF (fewer where the chromosome ends), one N between two
contigs. CRISPRitz 2.9.1 and earlier write LEFT = 25 and RIGHT = 27 (OLD_WINDOW); 2.10.0
writes 30 and 30, which its crispritz.py states as INDEL_LEFT and INDEL_RIGHT
(engine_window). The indel log, ``Dictionaries/log_indels_<vcf>/log<chrom>.txt[.gz]``,
records every window: its CHR field ``<chrom>_<start>-<end>_<id>`` (start and end 0-based,
end exclusive and not clipped at the chromosome's end) and its indel field
``<chrom>_<POS>_<REF>_<ALT>``. So a reader takes each window's left flank as
``POS - 1 - start`` and its right flank as ``end - (POS - 1 + len(REF))`` (row_flanks), and
never assumes a width: the windows of either release read right, and a reader that took
2.10.0's 30 bases for 2.9.1's 25 would put every flank base, SNP code and row of a window
five bases off.

A target of L bases (guide + PAM) with b DNA bulges holds a record's change on its window
only where L + b <= LEFT + 2 and L + b <= RIGHT + 1, or L + b <= LEFT + 1 for a record whose
REF and ALT differ at their first base, an MNV or complex record (reach): an insertion's
first added base, or the base after a deletion's anchor, is at offset LEFT + 1, an MNV's
first replaced base at LEFT, and the last base the ALT adds or replaces has RIGHT bases after
it. A target of more is held by no index, whatever the engine. CRISPRitz 2.10.0's window is
LEFT = RIGHT = Lmax - 1 + bmax, Lmax = 29 (the general K=8 index's flank and a 21-nt guide)
and bmax = 2 (the most DNA bulges the web offers): L + b <= 31.

Stdlib only.
"""

import gzip
import os

# The window of CRISPRitz 2.9.1 and earlier (enricher.py: genome[POS-26 : POS+26+len(REF)]),
# of every index published up to CRISPRme 2.8.1: 25 reference bases before the record's REF,
# 27 after it.
OLD_WINDOW = (25, 27)
# The window of CRISPRitz 2.10.0 (its crispritz.py INDEL_LEFT, INDEL_RIGHT), taken where no
# crispritz.py on PATH says (engine_window): the CRISPRitz CRISPRme 2.10.0 installs with.
WINDOW_2100 = (30, 30)
# The rows of a log read for its window: the left flank is the same for every window of one
# enrichment (one whose start would fall before base 0 is not written), the right one too
# (its end is not clipped at the chromosome's end), so the first rows say it.
ROWS = 200
# The record of an INDELS sub-index (index_fit.write_index_meta) that says the window of
# the fake-indel genome it was built from, [left, right].
META_KEY = "indel_window"


def window_start(chr_field):
    """The 0-based start of a window (``<chrom>_<start>-<end>_<id>``), read from the right,
    so a chromosome name holding '_' or '-' reads right too."""
    return int(chr_field.rsplit("_", 1)[0].rsplit("-", 1)[0].rsplit("_", 1)[1])


def window_end(chr_field):
    """The 0-based end (exclusive, unclipped) of a window, as window_start reads its start."""
    return int(chr_field.rsplit("_", 1)[0].rsplit("-", 1)[1])


def left_flank(pos1, start):
    """The reference bases a window holds before its record's REF: ``POS - 1 - start``. A
    ValueError for a row whose start is past its POS, which no enricher writes."""
    left = int(pos1) - 1 - int(start)
    if left < 0:
        raise ValueError("a window that starts at %d, after its record's POS %d" % (start, pos1))
    return left


def row_flanks(chr_field, indel_field):
    """(left, right): the reference bases a log row's window holds before its record's REF
    and after it (row_flanks of the CHR and indel fields)."""
    pos1, ref = indel_field.rsplit("_", 3)[1:3]
    pos1 = int(pos1)
    left = left_flank(pos1, window_start(chr_field))
    return left, window_end(chr_field) - (pos1 - 1 + len(ref))


def log_path(log_dir, chrom):
    """The indel log of ``chrom`` in ``log_dir`` (log<chrom>.txt.gz, else .txt), or None."""
    for name in ("log%s.txt.gz" % chrom, "log%s.txt" % chrom):
        p = os.path.join(log_dir, name)
        if os.path.isfile(p):
            return p
    return None


def log_windows(path, rows=ROWS):
    """The set of (left, right) of the windows of the indel log ``path``, from its first
    ``rows`` rows (every row with None): empty for a log with no row (a chromosome without
    indels). A row that is not of the log's layout raises ValueError."""
    out = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        n = 0
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0 and f[0] == "CHR":
                continue
            if len(f) < 5 or not f[0]:
                continue
            try:
                out.add(row_flanks(f[0], f[4]))
            except (IndexError, ValueError) as exc:
                raise ValueError("%s: a row that is no indel log row (%s)" % (path, exc))
            n += 1
            if rows is not None and n >= rows:
                break
    return out


def log_chroms(log_dir):
    """The chromosomes ``log_dir`` holds a log of (log<chrom>.txt[.gz]), sorted."""
    out = set()
    try:
        names = os.listdir(log_dir)
    except OSError:
        return []
    for name in names:
        stem = name[:-3] if name.endswith(".gz") else name
        if stem.startswith("log") and stem.endswith(".txt") and len(stem) > len("log.txt"):
            out.add(stem[len("log"):-len(".txt")])
    return sorted(out)


def dataset_windows(log_dir, chroms=None, rows=ROWS):
    """{chromosome: set of (left, right)} of the logs of ``log_dir`` (or of ``chroms``), from
    each one's first ``rows`` rows; a chromosome with no log, or no row, maps to an empty
    set."""
    out = {}
    for ch in (chroms if chroms is not None else log_chroms(log_dir)):
        p = log_path(log_dir, ch)
        out[ch] = log_windows(p, rows) if p else set()
    return out


def dataset_window(log_dir, chroms=None, rows=ROWS):
    """The one (left, right) of the windows of ``log_dir``'s logs, None where none holds a
    row; a ValueError where they hold two widths (logs of two enrichments)."""
    found = set()
    for ws in dataset_windows(log_dir, chroms, rows).values():
        found |= ws
    if len(found) > 1:
        raise ValueError("the indel logs of %s hold windows of %s: they come from more than one "
                         "enrichment" % (log_dir, " and ".join(describe(w) for w in sorted(found))))
    return next(iter(found)) if found else None


def reach(window, mnv=False):
    """The longest target (L + b: its bases, guide and PAM, and its DNA bulges) that holds a
    record's change on a window of ``window`` (left, right) reference bases: min(LEFT + 2,
    RIGHT + 1), min(LEFT + 1, RIGHT + 1) for an MNV or complex record (REF and ALT differ at
    their first base)."""
    left, right = window
    return min(left + (1 if mnv else 2), right + 1)


def describe(window):
    """'25/27' and the like: a window's left and right flanks."""
    return "%d/%d" % tuple(window)


def _index_fit():
    """index_fit, imported as a sibling module (crisprme.py, the search's scripts) or, where
    this module was imported as PostProcess.indel_window, from that package."""
    try:
        import index_fit
    except ImportError:
        from PostProcess import index_fit
    return index_fit


def engine_window():
    """The window the installed CRISPRitz's add-variants writes: its crispritz.py's
    INDEL_LEFT and INDEL_RIGHT (CRISPRitz 2.10.0 and later state them); OLD_WINDOW for a
    crispritz.py that states none (2.9.1 and earlier); WINDOW_2100 where no crispritz.py is
    on PATH to tell."""
    index_fit = _index_fit()

    left, right = index_fit._engine_value("INDEL_LEFT"), index_fit._engine_value("INDEL_RIGHT")
    if left is None or right is None:
        return WINDOW_2100
    if left.isdigit() and right.isdigit():
        return int(left), int(right)
    return OLD_WINDOW


def recorded_window(index_dir):
    """The window an INDELS sub-index records it was built from (META_KEY of its
    index_fit record), or None."""
    index_fit = _index_fit()

    w = (index_fit.read_index_meta(index_dir) or {}).get(META_KEY)
    if isinstance(w, list) and len(w) == 2 and all(isinstance(v, int) and not isinstance(v, bool)
                                                  for v in w):
        return tuple(w)
    return None


def index_window(index_dir):
    """The window the fake-indel genome of an INDELS sub-index held when it was built: the
    one it records (recorded_window, CRISPRme 2.10.0 and later), else OLD_WINDOW for one
    whose trees are of index format 1, which only CRISPRitz 2.9.1 and earlier build, whose
    enricher wrote no other; None where neither tells (no .bin, or format 2 with no
    record)."""
    index_fit = _index_fit()

    w = recorded_window(index_dir)
    if w is not None:
        return w
    try:
        formats = index_fit.bin_formats(index_dir)
    except OSError:
        return None
    return OLD_WINDOW if set(formats) == {1} else None


def mismatched_indexes(library, dataset, log_dir):
    """[(INDELS index, its window, the logs' windows)] of the INDELS indexes of ``dataset``
    in ``library`` (<...>+<dataset>_INDELS) whose sub-index of a chromosome was built from
    fake-indel windows (index_window) of another width than that chromosome's indel log in
    ``log_dir`` records: the index and the logs come from two enrichments, and every search
    refuses the index (pool_search_indels' _window_mismatch). One entry an index, of its
    first such sub-index."""
    out = []
    try:
        names = sorted(os.listdir(library))
    except OSError:
        return out
    for name in names:
        idx = os.path.join(library, name)
        if not name.endswith("+%s_INDELS" % dataset) or not os.path.isdir(idx):
            continue
        for sub in sorted(os.listdir(idx)):
            if "_fake" not in sub or not os.path.isdir(os.path.join(idx, sub)):
                continue
            w = index_window(os.path.join(idx, sub))
            path = log_path(log_dir, sub.rsplit("_fake", 1)[1])
            try:
                logged = log_windows(path) if path else set()
            except (OSError, ValueError, EOFError):
                continue
            if w is not None and logged and logged != {w}:
                out.append((idx, w, sorted(logged)))
                break
    return out


def target_note(window, target_len, bdna, who, dataset):
    """What a search of targets ``target_len`` long (guide + PAM) with up to ``bdna`` DNA
    bulges says of ``dataset``, whose fake-indel windows are ``window`` (left, right): None
    where every such target fits (reach); else a WARNING (``who``: the search's name for
    itself) that the indel off-targets of more bases than the window holds are not
    reported."""
    if window is None:
        return None
    want = target_len + bdna
    full, mnv = reach(window), reach(window, mnv=True)
    if want <= mnv:
        return None
    holds = ("up to %d bases" % full if full == mnv else
             "up to %d bases (%d at a record whose REF and ALT differ at their first base, an "
             "MNV or complex record)" % (full, mnv))
    how = "a %d-nt guide and PAM" % target_len + (
        " with %d DNA bulge%s" % (bdna, "" if bdna == 1 else "s") if bdna else "")
    advice = ("" if tuple(window) == WINDOW_2100 else
              " CRISPRitz 2.10.0 writes windows of 30 bases on either side, which hold up to "
              "31: enrich %s again with it (build-index-only) to search them." % dataset)
    return ("WARNING [%s]: the fake-indel windows of %s hold %d reference bases before each "
            "indel and %d after it, so a target that holds an indel there holds %s, guide, PAM "
            "and DNA bulges counted; this search's targets need up to %d (%s): such indel "
            "off-targets are held by no index and not reported.%s"
            % (who, dataset, window[0], window[1], holds, want, how, advice))

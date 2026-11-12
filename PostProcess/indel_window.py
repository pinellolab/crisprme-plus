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
never assumes a width: the windows of either release read right. A reader that took 2.10.0's
30 bases for 2.9.1's 25 would map a window's bases from its 26th through the ALT's last to
the wrong reference position, or to none where they have one (the bases before and after them
map alike under either width): the five flank bases before the REF would lose their code or
take a neighbour's, and the ALT's last five bases (all of a shorter ALT) would take codes,
where none belongs. Two exceptions: at an MNV, whose REF and ALT have one length, the flank
bases past the ALT's length map to their own position by chance (the 28th to the 30th at a
2-base MNV), and an ALT longer than five bases has its bases before its last five map to none
under either width, as they should. So the SNP codes on and around the REF and ALT, and the
SNP+indel rows next to the indel, would be wrong.

A target of L bases (guide + PAM) with b DNA bulges holds a record's change on its window
only where L + b <= LEFT + 2 and L + b <= RIGHT + 1, or L + b <= LEFT + 1 for a record whose
REF and ALT differ at their first base, an MNV or complex record (reach): an insertion's
first added base, or the base after a deletion's anchor, is at offset LEFT + 1, an MNV's
first replaced base at LEFT, and the last base the ALT adds or replaces has RIGHT bases after
it. A target of more is held by no index, whatever the engine. CRISPRitz 2.10.0's window is
LEFT = RIGHT = Lmax - 1 + bmax, Lmax = 29 (the general K=8 index's flank and a 21-nt guide)
and bmax = 2 (the most DNA bulges the web offers): L + b <= 31.

The windows also tell an enrichment of CRISPRitz 2.9.1 or earlier (25/27) from one of 2.10.0,
which differs in its SNPs too (2.9.1 enriched only the last ALT of a site written on several
records). The dictionaries, indel logs, enriched genome and fake-indel genome of a dataset, its
dict-less tiers and every index of it must come from one enrichment (enrichment_stores): a
dataset enriched anew beside the stores of an earlier enrichment of other windows is refused,
by build-index-only and by a search that enriches it itself (reenrichment_refusal; the search
script calls ``indel_window.py reenrichment``).

A sites-only dataset (sites_only: its VCFs name no sample; the mega) is the exception to
enriching anew: CRISPRitz writes a fake-indel window only for an indel a sample of the VCF
carries, so an enrichment of its VCFs holds none, and the windows it has come from an
enriched-genome kit whose indels a synthetic sample carries
(seq_script/merge_panels/mega_build_indels.sh). It is not enriched anew beside windows of
such a kit, whatever their width (reenrichment_refusal), and every advice on it points to its
archive and kit of CRISPRme 2.10.0, or that script, not to build-index-only (sites_only_way).
Its registry, built from its VCFs' INFO AF fields, is no store of an enrichment
(registry_data_type).

Stdlib only.
"""

import gzip
import json
import os
import sys

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


def target_note(window, target_len, bdna, who, dataset, sites=False):
    """What a search of targets ``target_len`` long (guide + PAM) with up to ``bdna`` DNA
    bulges says of ``dataset``, whose fake-indel windows are ``window`` (left, right): None
    where every such target fits (reach); else a WARNING (``who``: the search's name for
    itself) that the indel off-targets of more bases than the window holds are not
    reported, and, for windows narrower than CRISPRitz 2.10.0's, how to have its windows: by
    enriching the dataset again with build-index-only, or, for a sites-only one (``sites``:
    sites_only), whose windows no enrichment of its VCFs writes, as sites_only_way says."""
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
    if tuple(window) == WINDOW_2100:
        advice = ""
    elif sites:
        advice = (" CRISPRitz 2.10.0 writes windows of 30 bases on either side, which hold up "
                  "to 31. To search them: %s." % sites_only_way(dataset))
    else:
        advice = (" CRISPRitz 2.10.0 writes windows of 30 bases on either side, which hold up "
                  "to 31: enrich %s again with it (build-index-only) to search them." % dataset)
    return ("WARNING [%s]: the fake-indel windows of %s hold %d reference bases before each "
            "indel and %d after it, so a target that holds an indel there holds %s, guide, PAM "
            "and DNA bulges counted; this search's targets need up to %d (%s): such indel "
            "off-targets are held by no index and not reported.%s"
            % (who, dataset, window[0], window[1], holds, want, how, advice))


# ---- a dataset enriched anew: what an earlier enrichment of it left -----------------------

def quote(path):
    """``path`` single-quoted for a shell command line, a quote in it escaped."""
    return "'%s'" % str(path).replace("'", "'\\''")


def dataset_indexes(library, genome, dataset):
    """The variant and INDELS indexes of variant dataset ``dataset`` of genome ``genome`` in
    ``library`` (<PAM>_<N>_<genome>+<dataset>[_INDELS], either form of the name:
    index_fit.parse_index_name), as paths, sorted: not a hidden folder, nor a build lock."""
    index_fit = _index_fit()
    try:
        names = sorted(os.listdir(library))
    except OSError:
        return []
    rests = ("%s+%s" % (genome, dataset), "%s+%s_INDELS" % (genome, dataset))
    return [os.path.join(library, name) for name in names
            if not name.startswith(".") and (index_fit.parse_index_name(name) or "xxx")[2] in rests
            and os.path.isdir(os.path.join(library, name))]


def registry_data_type(registry):
    """The data type the registry folder ``registry`` gives in a chromosome's manifest
    (reg_<chrom>.idx; dataset-wide, so the first that tells is enough): its ``data_type``,
    else "sites-only" where its ``aggregation`` is "info_af" -- a registry built from VCFs'
    INFO AF fields (PostProcess/build_mega_registry.py) before the manifest had a data_type,
    as every published mega's -- as crisprme.py's _registry_data_type and
    tier0_registry._classify_data_type read it. None where no manifest tells."""
    try:
        names = sorted(n for n in os.listdir(registry) if n.startswith("reg_")
                       and n.endswith(".idx"))
    except OSError:
        return None
    for name in names:
        try:
            with open(os.path.join(registry, name)) as fh:
                manifest = json.load(fh)
        except (OSError, ValueError):
            continue
        if not isinstance(manifest, dict):
            continue
        kind = manifest.get("data_type") or (
            "sites-only" if manifest.get("aggregation") == "info_af" else None)
        if kind:
            return kind
    return None


def _sites_only_registry(registry):
    """Whether the registry ``registry`` is of a sites-only dataset (registry_data_type):
    built from its VCFs' INFO AF fields, not from an enrichment's dictionaries."""
    return registry_data_type(registry) == "sites-only"


def _vcf_names_a_sample(path):
    """Whether the bgzipped VCF ``path`` names a sample on its #CHROM line (no for one that
    cannot be read, as the input validator's _vcf_header_samples has it)."""
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    return len(line.rstrip("\n").split("\t")) > 9
                if not line.startswith("#"):
                    break
    except (OSError, EOFError, ValueError):
        pass
    return False


def sites_only(workdir, dataset, vcf_dir=None):
    """Whether variant dataset ``dataset`` of working directory ``workdir`` is sites-only:
    its VCFs name no sample, so an enrichment of them writes no fake-indel window (CRISPRitz's
    add-variants writes one for an indel a sample carries), and the windows it has come from
    an enriched-genome kit whose indels a synthetic sample carries
    (seq_script/merge_panels/mega_build_indels.sh, as the published mega's). Told by its
    VCFs, those of ``vcf_dir`` (as crisprme.py's _sites_only_vcfs: .vcf.gz files, none naming
    a sample), where there are any to read, else by its registry (registry_data_type)."""
    if vcf_dir:
        try:
            vcfs = sorted(os.path.join(vcf_dir, n) for n in os.listdir(vcf_dir)
                          if n.endswith(".vcf.gz"))
        except OSError:
            vcfs = []
        if vcfs:
            return not any(_vcf_names_a_sample(v) for v in vcfs)
    return _sites_only_registry(os.path.join(workdir, "Dictionaries", "registry_%s" % dataset))


def sites_only_why(dataset):
    """Why sites-only dataset ``dataset`` (sites_only) is not given fake-indel windows by
    enriching it anew from its VCFs."""
    return ("%s is sites-only: its VCFs name no sample, and CRISPRitz writes a fake-indel "
            "window only for an indel a sample carries, so an enrichment of them holds none"
            % dataset)


def sites_only_way(dataset):
    """How sites-only dataset ``dataset`` (sites_only) gets the fake-indel windows of
    CRISPRitz 2.10.0, not from build-index-only (sites_only_why): from its archive and
    enriched-genome kit of CRISPRme 2.10.0, or its kit built again with a synthetic sample
    (seq_script/merge_panels/mega_build_indels.sh)."""
    return ("%s (so build-index-only would leave it with no indel off-target); its windows "
            "come from its enriched-genome kit, whose indels a synthetic sample carries: "
            "install its archive and enriched-genome kit of CRISPRme 2.10.0 "
            "(docs/PRECOMPUTED_INDEXES.md), or build its kit again that way "
            "(seq_script/merge_panels/mega_build_indels.sh)" % sites_only_why(dataset))


def _holds_a_file(directory):
    try:
        return any(not n.startswith(".") for n in os.listdir(directory))
    except OSError:
        return False


def _indexes(workdir, genome, dataset, libraries):
    """The indexes of the dataset (dataset_indexes) in the working directory's genome_library
    and in each of ``libraries``, each once."""
    out, seen = [], set()
    for lib in [os.path.join(workdir, "genome_library")] + list(libraries):
        for idx in dataset_indexes(lib, genome, dataset):
            if os.path.realpath(idx) not in seen:
                seen.add(os.path.realpath(idx))
                out.append(idx)
    return out


def _tiers(workdir, dataset):
    """The dict-less tiers build-index-only emits from the dataset's dictionaries:
    Dictionaries/registry_<dataset> (not a sites-only one, built from the VCFs' INFO AF
    fields: _sites_only_registry) and genotypes_<dataset>."""
    dicts = os.path.join(workdir, "Dictionaries")
    registry = os.path.join(dicts, "registry_%s" % dataset)
    return ([] if _sites_only_registry(registry) else [registry]) + [
        os.path.join(dicts, "genotypes_%s" % dataset)]


def enrichment_stores(workdir, genome, dataset, libraries=(), present=True):
    """What of variant dataset ``dataset`` of genome ``genome``, in working directory
    ``workdir``, comes from its enrichment, as paths, in the order a removal names them:
    its enriched genome and fake-indel genome (Genomes/<genome>+<dataset>[_INDELS]), its
    dictionaries and indel logs (Dictionaries/dictionaries_<dataset>, log_indels_<dataset>),
    its dict-less tiers, which build-index-only emits from those dictionaries
    (Dictionaries/registry_<dataset>, genotypes_<dataset>), and its indexes, built from those
    genomes (dataset_indexes, of the working directory's genome_library and of each of
    ``libraries``). With ``present``, those that are there (a tier where it holds a file).
    Not what comes from its VCFs alone, the same from any enrichment of them: its phased
    indel genotype store, and a sites-only registry (the mega's, which
    PostProcess/build_mega_registry.py writes from its VCFs' INFO AF fields, as its manifest
    says: its data_type, or, in every published mega's, its aggregation, registry_data_type)
    with its indel frequencies."""
    own = [os.path.join(workdir, "Genomes", "%s+%s%s" % (genome, dataset, s))
           for s in ("", "_INDELS")] + [
        os.path.join(workdir, "Dictionaries", "%s_%s" % (k, dataset))
        for k in ("dictionaries", "log_indels")]
    tiers = _tiers(workdir, dataset)
    if present:
        own = [d for d in own if os.path.lexists(d)]
        tiers = [d for d in tiers if _holds_a_file(d)]
    return own + tiers + _indexes(workdir, genome, dataset, libraries)


def earlier_windows(workdir, genome, dataset, libraries=(), variant=True):
    """The fake-indel windows that the records of the enrichment of ``dataset`` of genome
    ``genome`` in ``workdir`` tell -- of an earlier one, where a new one is to be made: its
    indel logs' (Dictionaries/log_indels_<dataset>: dataset_windows), each INDELS sub-index's
    of the dataset (index_window: the one it records, or 25/27 for one of index format 1),
    and, with ``variant``, 25/27 for a variant index of the dataset of format 1 where the
    installed CRISPRitz writes format 2 (only CRISPRitz 2.9.1 and earlier build one, from an
    enrichment of theirs, whose windows were 25/27; without, only what holds windows tells,
    a log row or an INDELS sub-index). A set: empty where nothing tells (no log row, no such
    index)."""
    index_fit = _index_fit()
    found = set()
    try:
        for ws in dataset_windows(os.path.join(workdir, "Dictionaries",
                                               "log_indels_%s" % dataset)).values():
            found |= ws
    except (OSError, ValueError, EOFError):
        pass
    newer = (index_fit.engine_index_format() or 2) >= 2
    for idx in _indexes(workdir, genome, dataset, libraries):
        if idx.endswith("_INDELS"):
            try:
                subs = sorted(os.listdir(idx))
            except OSError:
                subs = []
            for sub in subs:
                if "_fake" in sub and os.path.isdir(os.path.join(idx, sub)):
                    w = index_window(os.path.join(idx, sub))
                    if w is not None:
                        found.add(w)
        elif variant and newer and set(index_fit.bin_formats(idx)) == {1}:
            found.add(OLD_WINDOW)
    return found


def reenrichment_refusal(workdir, genome, dataset, libraries=(), who="build-index-only",
                         vcf_dir=None):
    """Why ``who`` (build-index-only, or "this search": the search script's STEP 1, where a
    search enriches a dataset itself) must not enrich variant dataset ``dataset`` of genome
    ``genome`` anew in ``workdir`` from its VCFs (``vcf_dir``), its finished enrichment not in
    place; None where nothing stands in the way. Stores built from an earlier enrichment of
    it are in place, its tiers or its indexes (enrichment_stores; ``libraries``: those a
    search reads besides the working directory's), and that enrichment was another
    CRISPRitz's, as its records tell by fake-indel windows other than those the installed
    CRISPRitz writes (earlier_windows, engine_window): one of 2.9.1 or earlier, whose
    enrichment differs in its SNPs too (it enriched only the last ALT of a site written on
    several records). Enriched anew, the dataset would hold stores of two enrichments, which a
    search reads together: the alleles its indexes match and the carriers its tiers give
    would disagree. An earlier enrichment of the installed CRISPRitz's windows gives the same
    enrichment of the same VCFs: what was built from it is kept, as before. The removal named
    is of everything that comes from the earlier enrichment (enrichment_stores), with its
    dictionaries and indel logs, which a new enrichment writes again (a search's own moves
    its files in beside them).

    But not that of a sites-only dataset (sites_only, the mega) whose indel logs or INDELS
    indexes hold windows, whatever their width: an enrichment of its VCFs holds none, they
    came from an enriched-genome kit whose indels a synthetic sample carries, and enriched
    anew it would lose every indel off-target, its new indel logs holding no row (CRISPRme
    2.8.1 did so, silently). That kit is to be put back in place, or the dataset's archive
    and kit of CRISPRme 2.10.0 installed (sites_only_way)."""
    enriched = os.path.join(workdir, "Genomes", "%s+%s" % (genome, dataset))
    engine = engine_window()
    if sites_only(workdir, dataset, vcf_dir):
        kit = sorted(earlier_windows(workdir, genome, dataset, libraries, variant=False))
        if kit:
            return ("no finished enrichment of %s is in place (%s), so %s would enrich it anew "
                    "from its VCFs; but %s, while its indel logs or its INDELS indexes record "
                    "windows of %s, which came from an enriched-genome kit whose indels a "
                    "synthetic sample carries (seq_script/merge_panels/mega_build_indels.sh, "
                    "as the published mega's). Enriched anew, it would lose every indel "
                    "off-target: the new enrichment's indel logs hold no row, and its INDELS "
                    "indexes no window. Put that kit back in place (%s and %s_INDELS: "
                    "of a published dataset, its enriched-genome kit, "
                    "docs/PRECOMPUTED_INDEXES.md), from which a build or a search builds the "
                    "index it needs%s."
                    % (dataset, enriched, who, sites_only_why(dataset),
                       " and ".join(describe(w) for w in kit), enriched, enriched,
                       "" if kit == [engine] else
                       "; or, for the windows of the installed CRISPRitz (%s), install the "
                       "dataset's archive and enriched-genome kit of CRISPRme 2.10.0 in place "
                       "of all that came from that kit (rm -rf %s), or build its kit again that "
                       "way" % (describe(engine), " ".join(quote(d) for d in enrichment_stores(
                           workdir, genome, dataset, libraries)))))
    built = ([d for d in _tiers(workdir, dataset) if _holds_a_file(d)]
             + _indexes(workdir, genome, dataset, libraries))
    if not built:
        return None
    other = sorted(earlier_windows(workdir, genome, dataset, libraries) - {engine})
    if not other:
        return None
    then = ("then run build-index-only again, which enriches the dataset anew and builds its "
            "index, and its tiers with --samplesID, from that enrichment (a build or a search "
            "builds each other index of it again where it needs it)"
            if who == "build-index-only" else
            "then run build-index-only with the dataset's VCFs and its --samplesID, which "
            "enriches it anew and builds its tiers and an index from that enrichment, and "
            "run this search again")
    return ("no finished enrichment of %s is in place (%s), so %s would enrich it anew with "
            "the installed CRISPRitz, whose fake-indel windows are %s; but an earlier "
            "enrichment of it, whose records tell windows of %s (CRISPRitz 2.9.1 and earlier "
            "wrote 25/27), left what was built from it: %s. An enrichment of another CRISPRitz "
            "differs in its SNPs too (2.9.1 enriched only the last ALT of a site written on "
            "several records, 2.10.0 enriches every one), and a search would read the stores "
            "of the two enrichments together, with carriers or rows wrong or lost. Remove "
            "them, with the dataset's dictionaries and indel logs, which the new enrichment "
            "writes again (rm -rf %s), %s (of a published dataset, its archive and "
            "enriched-genome kit of CRISPRme 2.10.0 may be installed instead: "
            "docs/PRECOMPUTED_INDEXES.md); or put back the enrichment they were built from "
            "(%s and %s_INDELS), on which they are searched as they are."
            % (dataset, enriched, who, describe(engine), " and ".join(describe(w) for w in other),
               ", ".join(built), " ".join(quote(d) for d in enrichment_stores(
                   workdir, genome, dataset, libraries)), then, enriched, enriched))


# The bases an IUPAC code of the enriched genome stands for (the SNP overlay writes a flank
# base's code over it), Q (CRISPRitz 2.10.0's four bases) and N among them.
_STANDS_FOR = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
               "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
               "N": "ACGTN", "Q": "ACGT"}


def _fake_sequence(fake_fasta, upto):
    """The first ``upto`` + 1 bases of the fake-indel FASTA ``fake_fasta`` (its record's
    lines joined), or all of it where it is shorter."""
    out, n = [], 0
    with open(fake_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            line = line.strip()
            out.append(line)
            n += len(line)
            if n > upto:
                break
    return "".join(out)


def fake_window_disagreement(fake_fasta, log_path, rows=ROWS):
    """How the fake-indel genome of a chromosome, ``fake_fasta``, was not cut by the windows
    its indel log ``log_path`` records -- the two come from two enrichments, as where a
    download installs an archive's indel logs beside an enrichment made here by another
    CRISPRitz, whose windows were another width -- or None: each of the first ``rows`` rows'
    window (every row's with None) lies at the offsets of the fake sequence its FAKEPOS field
    gives, an N (the separator) or the sequence's end after it, and its left flank (POS - 1 -
    start) holds the reference bases its refseq field records, a SNP code standing for its own
    base (where the field holds them all). None too where there is nothing to compare, or
    nothing that tells: no log, no row, no fake FASTA, a row of no indel log's layout."""
    if not log_path or not os.path.isfile(fake_fasta):
        return None
    windows = []
    opener = gzip.open if log_path.endswith(".gz") else open
    with opener(log_path, "rt") as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0 and f[0] == "CHR":
                continue
            if len(f) < 7 or not f[0]:
                continue
            try:
                a, b = (int(x) for x in f[5].split(","))
                left = row_flanks(f[0], f[4])[0]
            except (IndexError, ValueError):
                return None   # (no log of CRISPRitz's layout: nothing to tell by)
            windows.append((f[4], a, b, left, f[6]))
            if rows is not None and len(windows) >= rows:
                break
    if not windows:
        return None
    seq = _fake_sequence(fake_fasta, max(b for _, _, b, _, _ in windows))
    for indel, a, b, left, refseq in windows:
        where = "%d-%d" % (a, b)
        if b > len(seq) or (b < len(seq) and seq[b] not in "Nn"):
            return ("the window of %s, which the log places at %s of the fake sequence, ends "
                    "elsewhere there (no separator after it)" % (indel, where))
        flank = seq[a:a + left]
        if len(refseq) >= left and (len(flank) < left or any(
                r.upper() not in _STANDS_FOR.get(c.upper(), c.upper())
                for c, r in zip(flank, refseq[:left]))):
            return ("the window of %s, which the log places at %s of the fake sequence, does "
                    "not hold there the %d reference bases before the record that the log "
                    "records" % (indel, where, left))
    return None


def fake_genome_disagreement(fake_dir, log_dir, dataset):
    """(chromosome, how) of the first chromosome of fake-indel genome ``fake_dir``
    (Genomes/<genome>+<dataset>_INDELS: fake_<dataset>_<chrom>/fake<chrom>.fa) whose
    windows are not those its indel log in ``log_dir`` records (fake_window_disagreement),
    or None."""
    prefix = "fake_%s_" % dataset
    try:
        names = sorted(n for n in os.listdir(fake_dir) if n.startswith(prefix))
    except OSError:
        return None
    for name in names:
        chrom = name[len(prefix):]
        how = fake_window_disagreement(os.path.join(fake_dir, name, "fake%s.fa" % chrom),
                                       log_path(log_dir, chrom))
        if how:
            return chrom, how
    return None


def fake_genome_refusal(workdir, genome, dataset, who="build-index-only", vcf_dir=None):
    """Why ``who`` (build-index-only, or "this search": the search script's STEP 1 with the
    enrichment in place and no variant index of the dataset that serves it, which it builds
    from that enrichment, and its STEP 2 before it puts the SNP codes on the fake-indel
    genome in place and builds the INDELS index it lacks from it) must not build an index of
    ``dataset`` from the enrichment in place in ``workdir``: its fake-indel genome and the
    dataset's indel logs come from two enrichments (fake_genome_disagreement), and so then do
    its enriched genome and the dataset's dictionaries and tiers, which a download installs
    with the logs; an index built from the one would be read with the other, every indel row
    misplaced (its windows are another width) and the SNP side another enrichment's, and the
    SNP overlay, which places each code by the logs, would rewrite the fake-indel genome with
    codes at the wrong bases. None where they agree. The second way out it names, enriching
    the dataset anew, is not that of a sites-only one (sites_only: ``vcf_dir`` its VCFs),
    whose windows no enrichment of its VCFs writes: building its kit again is
    (sites_only_way)."""
    enriched = os.path.join(workdir, "Genomes", "%s+%s" % (genome, dataset))
    logs = os.path.join(workdir, "Dictionaries", "log_indels_%s" % dataset)
    found = fake_genome_disagreement(enriched + "_INDELS", logs, dataset)
    if not found:
        return None
    every = " ".join(quote(d) for d in enrichment_stores(workdir, genome, dataset))
    anew = ("or build that kit again with a synthetic sample "
            "(seq_script/merge_panels/mega_build_indels.sh), in place of all of the dataset "
            "that comes from an enrichment (rm -rf %s): %s, so build-index-only would enrich "
            "it anew with no indel off-target" % (every, sites_only_why(dataset))
            if sites_only(workdir, dataset, vcf_dir) else
            "or remove all of the dataset that comes from an enrichment (rm -rf %s) and run "
            "build-index-only again, which enriches it anew" % every)
    return ("the fake-indel genome %s_INDELS and the indel logs %s of %s come from two "
            "enrichments: on %s, %s. With the fake-indel genome goes the enriched genome %s, "
            "and with the logs the dataset's dictionaries and dict-less tiers (a download "
            "installs an index's dictionaries and logs, not the enrichment they come from), so "
            "an index %s built from the enrichment in place would be read with another "
            "enrichment's logs, dictionaries and tiers: every indel row misplaced, and the SNP "
            "alleles and carriers of two enrichments mixed. Put the enrichment those come with "
            "in place of this one (rm -rf %s %s; of a published index, its enriched-genome "
            "kit, docs/PRECOMPUTED_INDEXES.md), %s."
            % (enriched, logs, dataset, found[0], found[1], enriched,
               "that build-index-only" if who == "build-index-only" else "that this search",
               quote(enriched), quote(enriched + "_INDELS"), anew))


def main(argv):
    """The search script's checks of a dataset, each the refusal on stdout, exit 3, or
    nothing, exit 0: ``indel_window.py reenrichment [--vcf <dir>] <workdir> <genome>
    <dataset> [<library>...]`` before it enriches the dataset itself (STEP 1:
    reenrichment_refusal); ``indel_window.py fakecheck [--vcf <dir>] <workdir> <genome>
    <dataset>`` before it builds an index of it from the enrichment in place (STEP 1, and
    STEP 2 before the SNP overlay of the fake-indel genome and its INDELS index:
    fake_genome_refusal). ``--vcf``: the dataset's VCF folder, which tells a sites-only one
    (sites_only)."""
    vcf_dir = None
    if len(argv) > 3 and argv[2] == "--vcf":
        vcf_dir, argv = argv[3], argv[:2] + argv[4:]
    if len(argv) < 5 or argv[1] not in ("reenrichment", "fakecheck"):
        sys.stderr.write("usage: indel_window.py reenrichment [--vcf <dir>] <workdir> <genome> "
                         "<dataset> [<library>...] | fakecheck [--vcf <dir>] <workdir> "
                         "<genome> <dataset>\n")
        return 2
    if argv[1] == "reenrichment":
        why = reenrichment_refusal(argv[2], argv[3], argv[4], argv[5:], who="this search",
                                   vcf_dir=vcf_dir)
    else:
        why = fake_genome_refusal(argv[2], argv[3], argv[4], who="this search",
                                  vcf_dir=vcf_dir)
    if why:
        print(why)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

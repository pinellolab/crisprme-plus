#!/usr/bin/env python
"""Phased indel genotype store (feature/indel-snp, gated).

The ``log_indels`` file lists each indel's carrier SAMPLES as a FLAT, UNPHASED
list, so it cannot prove that an indel and a nearby SNP sit on the SAME haplotype
(cis). CONFIRMED-cis SNP+indel reporting needs the indel's PHASED genotype per
carrier. The SNP side already has this (Tier-1 genotypes); this module produces
the indel-side equivalent.

It streams the build-time VCF (the same VCFs ``add-variants`` consumed), and for
every INDEL alt (len(REF) != len(alt)) of a record whose FILTER is PASS or '.' emits,
keyed by ``<pos>_<REF>_<ALT>``, the per-carrier phased genotype NORMALIZED to a
biallelic "this-indel vs not" call (so it is directly comparable to a SNP's biallelic
phased GT when checking cis). A record of any other FILTER is skipped, as CRISPRitz's
enricher skips it: no fake-indel window holds its indel, so a carrier is never taken
off a row for it (from 2.10.0; a store written before kept every indel of the VCF,
and build-index-only writes such a store again: store_is_current).

From 2.10.1 the store holds the MNVs too (_is_mnv: REF and ALT of one length above 1,
the equal-length records CRISPRitz's enricher cuts a fake-indel window for, as for an
indel), in a file of their own beside the indel file, in the same format and with the
same .done marker: ``gt_indel_<chrom>.tsv.gz`` and ``gt_mnv_<chrom>.tsv.gz``
(store_files). Up to 2.10.0 no MNV was stored, so nothing a search read said which
haplotype holds one. A reader of 2.10.0 opens the indel file by name, so the MNV file
leaves its readings as they were, and its presence tells a reader of 2.10.1 that the
MNVs were stored (IndelGenotypeReader.has_mnv). compile_store writes both files in one
pass over the VCF, each only where it is asked for: a store of 2.10.0 gains its MNV file
and keeps its indel file byte for byte (build-index-only's STEP 1c and
search_indel_snp.store ask for the files that are not current: files_to_write). A
reader takes an MNV file it cannot read to its end for none, the indel file as it is
(read_mnv_lines).

STDLIB ONLY (gzip). Runs at build time; the reader loads one chromosome's indels
into a dict (indels are ~10x fewer than SNPs, so this is cheap). Gated behind
CRISPRME_INDEL_SNP; classic build untouched.
"""

import contextlib
import errno
import os
import shlex
import sys
import gzip

import site_alleles

# The FILTER values whose records the store holds: those CRISPRitz's enricher reads
# (Enrichment/enricher.py and enricher.cpp, `if line[6] not in ('PASS', '.')`, the line
# skipped; validate_inputs.ENRICHER_PASS_VALUES), so that every indel of the store is one
# a fake-indel window holds.
PASS_FILTERS = ("PASS", ".")
# The second line of a store's .done marker (after its record count) from 2.10.0: the
# store holds only the records of PASS_FILTERS. A store written before has no such line.
FILTER_MARK = "FILTER\tPASS,."


def _is_indel(ref, alt):
    """True iff (ref, alt) is an indel (length differs) and not a symbolic allele."""
    return len(ref) != len(alt) and "<" not in alt and "*" not in alt


def _is_mnv(ref, alt):
    """True iff (ref, alt) is an MNV: REF and ALT of one length above 1 (AA>GG), not a
    symbolic allele. These are the equal-length records CRISPRitz's enricher cuts a
    fake-indel window for, as for an indel (an SNV, REF and ALT of one base, it codes in the
    SNP genome instead), and that the phased store holds in its MNV file (from 2.10.1)."""
    return len(ref) == len(alt) > 1 and "<" not in alt and "*" not in alt


def store_files(store_dir, chrom):
    """(indel file, MNV file) of the phased store of ``chrom`` in its folder ``store_dir``
    (``Dictionaries/indel_genotypes_<vcf>/``): ``gt_indel_<chrom>.tsv.gz`` and, from
    2.10.1, ``gt_mnv_<chrom>.tsv.gz``."""
    return (os.path.join(store_dir, f"gt_indel_{chrom}.tsv.gz"),
            os.path.join(store_dir, f"gt_mnv_{chrom}.tsv.gz"))


def mnv_store_path(indel_path):
    """The MNV file beside a store's indel file (store_files: ``gt_mnv_<chrom>.tsv.gz`` for
    ``gt_indel_<chrom>.tsv.gz``), or None for a file not so named."""
    folder, name = os.path.split(indel_path)
    if not name.startswith("gt_indel_"):
        return None
    return os.path.join(folder, "gt_mnv_" + name[len("gt_indel_"):])


def normalize_gt_for_alt(gt, alt_idx):
    """Normalize a phased/unphased VCF genotype string to a biallelic
    "this-indel-alt vs everything-else" call, preserving the phase separator.

    ``gt`` is the GT subfield (e.g. "1|0", "0|1", "1|2", "2|2", "1/2", "0|0").
    ``alt_idx`` is the 1-based ALT index of the indel of interest. Each allele
    slot becomes "1" iff it equals ``alt_idx``, else "0". Returns the normalized
    2-slot string, or ``None`` if the sample carries the indel on NO slot (so it
    is not a carrier of THIS indel).

    Examples (alt_idx=1): "1|0"->"1|0", "0|1"->"0|1", "1|1"->"1|1", "0|0"->None.
    Examples (alt_idx=2): "1|2"->"0|1", "2|2"->"1|1", "1|0"->None, "1/2"->"0/1".
    A '/' (unphased) separator is preserved so the downstream cis check forces
    PUTATIVE (never CONFIRMED) for unphased samples.
    """
    gt = gt.split(":", 1)[0]  # GT is always the first FORMAT subfield
    if "|" in gt:
        sep, slots = "|", gt.split("|")
    elif "/" in gt:
        sep, slots = "/", gt.split("/")
    else:
        slots = [gt]  # haploid (e.g. chrY / male chrX)
        sep = ""
    tgt = str(alt_idx)
    norm = ["1" if s == tgt else "0" for s in slots]
    if "1" not in norm:
        return None  # not a carrier of this alt
    return sep.join(norm) if sep else norm[0]


def iter_store_genotypes(vcf_path, keep_samples=None, kinds=("indel", "mnv")):
    """Stream a (b)gzipped or plain VCF and yield, per indel or MNV alt (``kinds``: those
    of "indel" and "mnv" to yield): ``(kind, pos:int, ref:str, alt:str, {sample_id:
    normalized_gt})``, in the VCF's order.

    Multiallelic sites emit one tuple per such ALT (SNP alts are skipped), and a record
    whose FILTER is not PASS_FILTERS emits none. ``keep_samples`` (a set) restricts the
    emitted carriers to that panel; None keeps all. Only carriers (normalize_gt_for_alt !=
    None) are stored: a genotype is read as "this ALT or not", its separator kept, a
    missing allele read as not this ALT (1|. is 1|0, 1/. is 1/0) and a haploid call as
    one slot.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        sample_cols = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_cols = line.rstrip("\n").split("\t")[9:]
                continue
            if sample_cols is None:
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10 or c[6] not in PASS_FILTERS:
                continue
            pos, ref, alt_field = c[1], c[3], c[4]
            gts = c[9:]
            for a_i, alt in enumerate(alt_field.split(","), start=1):
                kind = ("indel" if _is_indel(ref, alt) else
                        "mnv" if _is_mnv(ref, alt) else None)
                if kind is None or kind not in kinds:
                    continue
                carriers = {}
                for col, gt in enumerate(gts):
                    sid = sample_cols[col]
                    if keep_samples is not None and sid not in keep_samples:
                        continue
                    ng = normalize_gt_for_alt(gt, a_i)
                    if ng is not None:
                        carriers[sid] = ng
                if carriers:
                    yield kind, int(pos), ref, alt, carriers


def iter_indel_genotypes(vcf_path, keep_samples=None):
    """Stream a (b)gzipped or plain VCF and yield, per indel alt:
    ``(pos:int, ref:str, alt:str, {sample_id: normalized_gt})``.

    Multiallelic sites emit one tuple per indel ALT (SNP and MNV alts are skipped), and a
    record whose FILTER is not PASS_FILTERS emits none. ``keep_samples`` (a set) restricts
    the emitted carriers to that panel; None keeps all. Only carriers (normalize_gt_for_alt
    != None) are stored.
    """
    for _kind, pos, ref, alt, carriers in iter_store_genotypes(vcf_path, keep_samples,
                                                               ("indel",)):
        yield pos, ref, alt, carriers


def iter_mnv_genotypes(vcf_path, keep_samples=None):
    """iter_indel_genotypes for the MNV alts (_is_mnv)."""
    for _kind, pos, ref, alt, carriers in iter_store_genotypes(vcf_path, keep_samples,
                                                               ("mnv",)):
        yield pos, ref, alt, carriers


def _done_marker(out_path):
    """Sidecar path written ONLY after ``out_path`` is fully in place (record count)."""
    return out_path + ".done"


def store_is_complete(out_path):
    """True iff the indel store at ``out_path`` is a complete, non-truncated file.

    Fast path: a ``.done`` marker (written only after the atomic rename) exists and is
    at least as new as the store. Otherwise fall back to a streaming gzip integrity
    check (equivalent to ``gzip -t``) so a legacy store written by the old non-atomic
    code, or one truncated by a hard kill/OOM mid-write, is detected and re-built on
    resume instead of being silently skipped as "done" (the gt_indel_chr1 truncation
    bug). A size>0 file that fails to fully decompress returns False."""
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        return False
    marker = _done_marker(out_path)
    try:
        if os.path.isfile(marker) and os.path.getmtime(marker) >= os.path.getmtime(out_path):
            return True
    except OSError:
        pass
    # No trusted marker -> verify the gzip stream end-to-end (catches truncation).
    try:
        with gzip.open(out_path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def store_reads_filter(out_path):
    """True iff the store at ``out_path`` was written under the FILTER rule of 2.10.0 (its
    records of PASS_FILTERS only): its .done marker, as new as the store, says so
    (FILTER_MARK). A store written before 2.10.0, or adopted without a marker, does not."""
    marker = _done_marker(out_path)
    try:
        if os.path.getmtime(marker) < os.path.getmtime(out_path):
            return False
        with open(marker) as m:
            return FILTER_MARK in m.read().splitlines()[1:]
    except OSError:
        return False


def store_is_current(out_path):
    """True iff the store at ``out_path`` is complete (store_is_complete) and written under
    the FILTER rule (store_reads_filter): one a build need not write again."""
    return store_is_complete(out_path) and store_reads_filter(out_path)


def compile_indel_genotypes(vcf_path, out_path, keep_samples=None):
    """Write the phased indel genotype store for one chromosome's VCF.

    Format (gzipped TSV, one line per indel alt):
        <pos>_<REF>_<ALT> \t sample:gt,sample:gt,...
    Returns the number of indel records written. The .done marker holds that count, then
    FILTER_MARK.

    The write is ATOMIC: rows stream into ``<out_path>.tmp`` and the temp file is
    ``os.replace``d onto the final path only after a clean close, then a ``.done``
    marker is dropped. A SIGKILL/OOM/interrupt mid-write leaves only the temp file
    (removed here), never a truncated final file that a resume would accept as done.
    """
    tmp_path = out_path + ".tmp"
    marker = _done_marker(out_path)
    try:
        os.remove(marker)  # a stale marker must never outlive a rebuild
    except OSError:
        pass
    n = 0
    try:
        with gzip.open(tmp_path, "wt") as out:
            for pos, ref, alt, carriers in iter_indel_genotypes(vcf_path, keep_samples):
                toks = ",".join(f"{s}:{g}" for s, g in carriers.items())
                out.write(f"{pos}_{ref}_{alt}\t{toks}\n")
                n += 1
        os.replace(tmp_path, out_path)  # atomic: final path is complete-or-absent
    except BaseException:  # noqa: BLE001 - incl. KeyboardInterrupt/SystemExit/kill
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass
        raise
    _mark_done(out_path, n)
    return n


def _mark_done(out_path, n):
    """Drop the .done marker of a store file written whole at ``out_path`` (its record
    count ``n``, then FILTER_MARK)."""
    marker = _done_marker(out_path)
    try:
        with open(marker, "w") as m:
            m.write(f"{n}\n{FILTER_MARK}\n")
    except OSError as err:
        # a marker the disk could not take is no marker; on a full disk or an I/O
        # error (ends_the_build) its store goes too, as the build ends on it, for the
        # next run to write both again
        for _p in (marker, out_path) if ends_the_build(err) else (marker,):
            try:
                os.remove(_p)
            except OSError:
                pass
        if ends_the_build(err):
            raise


def files_to_write(store_dir, chrom, write_dir=None):
    """(indel file, MNV file) for compile_store to write of ``chrom``'s store in
    ``store_dir`` (store_files): each one missing or not current there (store_is_current),
    at its name in ``write_dir`` (None: ``store_dir`` itself, written in place); None for
    one that is current. A store of 2.10.0 asks for its MNV file alone, so its indel file
    stays byte for byte; a store of 2.10.1 whose files are current asks for none."""
    out = []
    for path in store_files(store_dir, chrom):
        if store_is_current(path):
            out.append(None)
        else:
            out.append(path if write_dir is None
                       else os.path.join(write_dir, os.path.basename(path)))
    return tuple(out)


def compile_store(vcf_path, indel_out=None, mnv_out=None, keep_samples=None):
    """Write the files of one chromosome's phased store that are asked for, in ONE pass
    over its VCF: the indel file ``indel_out`` (as compile_indel_genotypes writes it) and
    the MNV file ``mnv_out`` (the same format and marker, for its MNV alts: _is_mnv); None
    asks for no such file. A caller asks only for the files that are missing or not current
    (files_to_write), so that a store of 2.10.0 gains its MNV file and keeps its indel
    file byte for byte. A VCF without an MNV gets an MNV file of no record, which says that
    the MNVs were stored (none). Returns (indel records, MNV records), None for a file not
    asked for.

    Each file is written as compile_indel_genotypes writes its own: into ``<file>.tmp``,
    moved onto its path once both are written whole, then marked done; a stale marker is
    removed first, and a failure removes the temporary files and leaves no part of a file
    at its path."""
    outs = [(kind, path) for kind, path in (("indel", indel_out), ("mnv", mnv_out)) if path]
    if not outs:
        return None, None
    for _kind, path in outs:
        try:
            os.remove(_done_marker(path))  # a stale marker must never outlive a rebuild
        except OSError:
            pass
    tmps = {kind: path + ".tmp" for kind, path in outs}
    counts = {kind: 0 for kind, _path in outs}
    try:
        with contextlib.ExitStack() as stack:
            fhs = {kind: stack.enter_context(gzip.open(tmps[kind], "wt"))
                   for kind, _path in outs}
            for kind, pos, ref, alt, carriers in iter_store_genotypes(vcf_path, keep_samples,
                                                                      tuple(fhs)):
                toks = ",".join(f"{s}:{g}" for s, g in carriers.items())
                fhs[kind].write(f"{pos}_{ref}_{alt}\t{toks}\n")
                counts[kind] += 1
        for kind, path in outs:
            os.replace(tmps[kind], path)  # atomic: each final path is complete-or-absent
    except BaseException:  # noqa: BLE001 - incl. KeyboardInterrupt/SystemExit/kill
        for tmp in tmps.values():
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass
        raise
    for kind, path in outs:
        _mark_done(path, counts[kind])
    return counts.get("indel"), counts.get("mnv")


def ends_the_build(exc):
    """True for a failure that every chromosome after this one would hit too: a full
    disk or quota (ENOSPC, EDQUOT), an I/O error (EIO), no memory left. The guarded
    writer raises it again rather than skip the chromosome, and build-index-only
    ends with an ERROR on it, the store left marked unfinished."""
    return isinstance(exc, MemoryError) or (
        isinstance(exc, OSError) and exc.errno in (errno.ENOSPC, errno.EDQUOT, errno.EIO))


def compile_indel_genotypes_safe(vcf_path, out_path, keep_samples=None):
    """Guarded ``compile_indel_genotypes`` for parallel fan-out: returns the record
    count, or -1 on any other failure so one bad chromosome cannot sink the whole
    store; it raises only a failure that ends the build (``ends_the_build``). The
    atomic writer already cleans up its temp file; here we also drop any stale
    final/marker so a later resume re-attempts cleanly."""
    try:
        return compile_indel_genotypes(vcf_path, out_path, keep_samples=keep_samples)
    except Exception as err:  # noqa: BLE001 - store is optional (-> PUTATIVE)
        for _p in (out_path + ".tmp", _done_marker(out_path)):
            try:
                if os.path.exists(_p):
                    os.remove(_p)
            except OSError:
                pass
        if ends_the_build(err):
            raise
        print(
            f"WARNING [indel-snp]: compile_indel_genotypes failed for "
            f"{os.path.basename(vcf_path)} ({err})",
            flush=True,
        )
        return -1


def compile_store_safe(vcf_path, indel_out=None, mnv_out=None, keep_samples=None):
    """Guarded ``compile_store``, as compile_indel_genotypes_safe guards its writer: its
    (indel records, MNV records) tuple, each None for a file not asked for, or None on any
    other failure, the temporary files and the markers of the files asked for removed; it
    raises only a failure that ends the build (``ends_the_build``). So a caller tells a
    chromosome written (a tuple) from one skipped with a WARNING (None) by ``is not None``
    (build-index-only's STEP 1c). compile_indel_genotypes_safe keeps its own contract (a
    count, or -1)."""
    try:
        return compile_store(vcf_path, indel_out, mnv_out, keep_samples=keep_samples)
    except Exception as err:  # noqa: BLE001 - store is optional (-> PUTATIVE)
        for out in (indel_out, mnv_out):
            if not out:
                continue
            for _p in (out + ".tmp", _done_marker(out)):
                try:
                    if os.path.exists(_p):
                        os.remove(_p)
                except OSError:
                    pass
        if ends_the_build(err):
            raise
        print(
            f"WARNING [indel-snp]: compile_store failed for "
            f"{os.path.basename(vcf_path)} ({err})",
            flush=True,
        )
        return None


def merge_store_carriers(earlier, later):
    """The carriers ("sample:gt,...") of one indel that two lines of the store give, the
    earlier line's first: two VCF records of one REF>ALT at a position are one allele
    (decision B2, as CRISPRitz 2.10.0's enricher merges two listings of one SNV), so every
    carrier of either is kept, and a carrier of both gets its two calls merged as the
    enricher merges them (site_alleles.merge_gt: a slot holds the indel where either call
    puts it there when they are phased alike; else the call with more copies, the earlier
    on a tie). The store writes one line per record, and up to 2.10.0's integration every
    reader kept the later line alone, so a carrier only the earlier line listed lost its
    phased genotype: CRISPRme listed it on no row of the indel, or on some alignments only.
    The readers key both lines in upper case (store_key), so two lines in two cases are
    merged too. HPRC chr22 writes 179 indels on two lines: 39 in one case (1,400 listings on
    the earlier only, which every reader lost) and 140 in two (read as two indels before)."""
    if not earlier:
        return later
    if not later:
        return earlier
    a, b = {}, {}
    for into, toks in ((a, earlier), (b, later)):
        for t in toks.split(","):
            sid, sep, gt = t.partition(":")
            if sid and sep:
                into[sid] = site_alleles.merge_gt(into[sid], gt) if sid in into else gt
    return ",".join("%s:%s" % (s, g) for s, g in site_alleles.merge_carriers(a, b).items())


def store_key(pos, ref, alt):
    """The key a reader of the store gives an indel: ``<pos>_<REF>_<ALT>`` with REF and
    ALT in upper case. The store writes each record's key as its VCF line does, and VCF
    4.x reads REF and ALT whatever their case, as the fake-indel window does (CRISPRitz
    writes the ALT in place of the REF, the search reads the window in upper case), so
    two lines of one indel in two cases (tg>t beside TG>T) are one allele here as two
    lines written alike are (merge_store_carriers); the rows name it in upper case
    (indel_snp_cis.record_name)."""
    return f"{pos}_{ref}_{alt}".upper()


def store_lines(path):
    """Every line of a file of the store, read to its end before any is taken: a file cut
    short (a copy stopped partway, a disk that lost a block) or corrupt raises here
    (EOFError, OSError, zlib.error, UnicodeDecodeError...), so that a reader takes none of
    its lines rather than a part of them."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.readlines()


def marked_count(path):
    """The record count that the .done marker of the store file at ``path`` gives (its first
    line), where that marker is as new as the file (the marker store_is_complete trusts);
    None where there is no such marker, or its first line is not a count."""
    marker = _done_marker(path)
    try:
        if os.path.getmtime(marker) < os.path.getmtime(path):
            return None
        with open(marker) as m:
            return int(m.readline())
    except (OSError, ValueError):
        return None


def read_mnv_lines(mnv_path):
    """(the lines of a chromosome's MNV file, None), read whole (store_lines), or (None, why)
    where it cannot be read to its end: corrupt or cut short. Each reader of the store (the
    indel post-analysis' IndelGenotypeReader, the SNP rows' new_simple_analysis._SnpRowIndels,
    the carrier gate's haplotype_prioritizer.IndelIndex) then reads that chromosome as a
    store of 2.10.0 is read, with no MNV file, its indel file as it is, and says so on
    stdout (mnv_unreadable_note); a corrupt indel file fails the store as before 2.10.1.
    Before, a corrupt MNV file failed the whole store: the indel post-analysis and the SNP
    rows read the chromosome with no store at all, and the gate's worker raised out of its
    pool.

    Two files that gzip reads to their end are not read either. A file of 0 bytes (a copy or
    a write stopped at once): gzip reads it as a file of no line, but compile_store writes
    even a file of no MNV as a gzip stream, and store_is_complete already takes an empty
    file for a missing one, so build-index-only writes it again and the download step
    fetches it. Read, it said that the chromosome's MNVs were read by their phase while none
    was. And a file whose .done marker, as new as the file (marked_count), counts other
    records than the file holds lines: it is not the file the marker was written for."""
    try:
        if os.path.getsize(mnv_path) == 0:
            return None, "0 bytes: even a file of no MNV is a gzip stream"
        lines = store_lines(mnv_path)
    except MemoryError:
        raise
    except Exception as err:  # noqa: BLE001 - not read to its end: no MNV file
        return None, (str(err) or type(err).__name__)
    n = marked_count(mnv_path)
    if n is not None and n != len(lines):
        return None, f"its .done marker counts {n} record(s), the file holds {len(lines)} line(s)"
    return lines, None


def mnv_unreadable_note(vcf, chrom, mnv_path, why, reading):
    """The NOTE (stdout: the run's own record of it, never stderr) of a reader that read
    the MNV file ``mnv_path`` of ``vcf``'s store on ``chrom`` as absent (read_mnv_lines:
    ``why``); ``reading``: how that reader then reads the chromosome's MNVs. It opens as
    the NOTE of a store with no MNV file does ("the phased indel genotype store of <vcf> on
    <chrom> holds no MNV"), the reading being the same, and says how to get the file
    whole: a corrupt file can keep a .done marker, so build-index-only takes it for
    current until it is removed."""
    q = shlex.quote
    return (f"NOTE [indel-snp]: the phased indel genotype store of {vcf} on {chrom} holds no "
            f"MNV that can be read: its MNV file {os.path.basename(mnv_path)} is corrupt or "
            f"cut short ({why}), so {chrom} is read as with a store written before CRISPRme "
            f"2.10.1, which has no MNV file: {reading}. To read them by their phase, remove "
            f"it (rm -f {q(mnv_path)} {q(_done_marker(mnv_path))}) and run build-index-only "
            f"with the dataset's VCFs, which writes it again, or, for an index downloaded "
            f"from Hugging Face, run crisprme.py download --what indel-genotypes --dataset "
            f"{vcf}, which installs the store whole again.")


class IndelGenotypeReader(object):
    """Read a phased indel genotype store (one chromosome) into a dict keyed by
    ``<pos>_<REF>_<ALT>``, REF and ALT in upper case (store_key). Small enough to hold
    in memory (indels << SNPs). ``reads_filter``: whether the store holds only the
    records of PASS_FILTERS (store_reads_filter); one written before 2.10.0 may hold the
    others too. A key the store holds twice (an indel the VCF writes on two lines, in
    one case or in two) is one allele, its carriers merged (merge_store_carriers).

    ``store_path`` is the store's indel file; its MNV file (``mnv_path``; None: the one
    beside it, mnv_store_path) is read too where it is there (from 2.10.1), its MNVs
    records of the store as its indels are (carriers, carriers_dict, and indels(), which
    gives both kinds). ``has_mnv``: whether an MNV file was read, so whether the store's
    MNVs can be told from no MNV at all (a store of 2.10.0 has no MNV file: its readers
    read the MNVs, from the indel log, as 2.10.0 did).

    An MNV file that cannot be read to its end (corrupt, cut short, 0 bytes, or not the file
    its .done marker counts: read_mnv_lines) is read as no MNV file, the indel file as it
    is: ``has_mnv`` False, and ``mnv_error`` says why (None where the file was read, or is
    not there), for the reader's NOTE (mnv_unreadable_note). An indel file that cannot be
    read raises, as before 2.10.1."""

    __slots__ = ("_by_key", "reads_filter", "has_mnv", "mnv_error", "mnv_path")

    def __init__(self, store_path, mnv_path=None):
        self._by_key = {}
        self.reads_filter = store_reads_filter(store_path)
        if mnv_path is None:
            mnv_path = mnv_store_path(store_path)
        self.mnv_path = mnv_path or None
        self.has_mnv = False
        self.mnv_error = None
        opener = gzip.open if store_path.endswith(".gz") else open
        with opener(store_path, "rt") as fh:
            self._take(fh)
        if mnv_path and os.path.isfile(mnv_path):
            lines, self.mnv_error = read_mnv_lines(mnv_path)
            if lines is not None:
                self._take(lines)
                self.has_mnv = True

    def _take(self, lines):
        for line in lines:
            key, _, toks = line.rstrip("\n").partition("\t")
            if not key:
                continue
            key = key.upper()  # store_key: the position is digits
            prev = self._by_key.get(key)
            self._by_key[key] = toks if prev is None else merge_store_carriers(prev, toks)

    def carriers(self, pos, ref, alt):
        """Return "sample:gt" tokens (comma-joined) for the indel, or "" if absent."""
        return self._by_key.get(store_key(pos, ref, alt), "")

    def carriers_dict(self, pos, ref, alt):
        """Return {sample_id: normalized_gt} for the indel (empty if absent)."""
        toks = self.carriers(pos, ref, alt)
        out = {}
        for t in toks.split(","):
            if not t:
                continue
            s, _, g = t.partition(":")
            out[s] = g
        return out

    def indels(self):
        """Every (pos, REF, ALT) the store holds, sorted by position: its indels, and its
        MNVs where it has an MNV file (has_mnv)."""
        out = []
        for key in self._by_key:
            pos, _, rest = key.partition("_")
            ref, _, alt = rest.partition("_")
            try:
                out.append((int(pos), ref, alt))
            except ValueError:
                continue
        out.sort()
        return out


def _load_panel(samplesid_path):
    """Load the set of panel sample IDs from a samplesID file (skip header)."""
    keep = set()
    opener = gzip.open if samplesid_path.endswith(".gz") else open
    with opener(samplesid_path, "rt") as fh:
        for i, line in enumerate(fh):
            s = line.strip()
            if not s or (i == 0 and (s.lower().startswith("sample") or "\t" in s)):
                # tolerate a header row; otherwise treat the first token as an ID
                pass
            keep.add(s.split("\t")[0])
    keep.discard("SAMPLE_ID")
    return keep


def vcf_chrom(vcf_path):
    """The chromosome of a VCF as build-index-only names its store's files: its first
    record's CHROM, with a chr prefix; None for a VCF with no record."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                chrom = line.split("\t", 1)[0].strip()
                return chrom if chrom.startswith("chr") else "chr" + chrom
    return None


USAGE = ("usage: build_indel_genotypes.py <vcf> <out.tsv.gz> [samplesID]\n"
         "       build_indel_genotypes.py --mnv-only <vcf> <store folder> [samplesID]\n")


def main(argv):
    """CLI: build_indel_genotypes.py <vcf> <out_store.tsv.gz> [<samplesID>]
    writes the indel file of a chromosome's store at <out_store.tsv.gz>, as up to 2.10.0.

    build_indel_genotypes.py --mnv-only <vcf> <store folder> [<samplesID>]
    writes the MNV file (gt_mnv_<chrom>.tsv.gz, store_files) of the VCF's chromosome
    (vcf_chrom) in the store's folder (Dictionaries/indel_genotypes_<vcf>/), and nothing
    else: the folder's other files stay as they are, so a store of 2.10.0 gains its MNVs
    (HPRC's rebuild of 2.10.1). The samplesID must be the build's (build-index-only keeps
    the samples of its --samplesID files): a samplesID file named but not found is an
    error here, where the first form reads every sample of the VCF.

    Enabled by default (opt-out): no-op only if CRISPRME_INDEL_SNP=0.
    """
    if os.environ.get("CRISPRME_INDEL_SNP", "1") not in ("1", "true", "True", "yes"):
        return 0
    if len(argv) > 1 and argv[1] == "--mnv-only":
        if len(argv) not in (4, 5):
            sys.stderr.write(USAGE)
            return 2
        vcf_path, store_dir = argv[2], argv[3]
        if not os.path.isdir(store_dir):
            sys.stderr.write(f"build_indel_genotypes: no store folder {store_dir}\n")
            return 2
        keep = None
        if len(argv) == 5:
            if not os.path.isfile(argv[4]):
                sys.stderr.write(f"build_indel_genotypes: no samplesID file {argv[4]}\n")
                return 2
            keep = _load_panel(argv[4])
        chrom = vcf_chrom(vcf_path)
        if chrom is None:
            print(f"build_indel_genotypes: {vcf_path} holds no record: no MNV file written",
                  flush=True)
            return 0
        out_path = store_files(store_dir, chrom)[1]
        _, n = compile_store(vcf_path, mnv_out=out_path, keep_samples=keep)
        print(f"build_indel_genotypes: {n} phased MNV records -> {out_path}", flush=True)
        return 0
    if len(argv) < 3:
        sys.stderr.write(USAGE)
        return 2
    vcf_path, out_path = argv[1], argv[2]
    keep = _load_panel(argv[3]) if len(argv) > 3 and os.path.isfile(argv[3]) else None
    n = compile_indel_genotypes(vcf_path, out_path, keep_samples=keep)
    print(f"build_indel_genotypes: {n} phased indel records -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

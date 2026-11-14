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

STDLIB ONLY (gzip). Runs at build time; the reader loads one chromosome's indels
into a dict (indels are ~10x fewer than SNPs, so this is cheap). Gated behind
CRISPRME_INDEL_SNP; classic build untouched.
"""

import errno
import os
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


def iter_indel_genotypes(vcf_path, keep_samples=None):
    """Stream a (b)gzipped or plain VCF and yield, per indel alt:
    ``(pos:int, ref:str, alt:str, {sample_id: normalized_gt})``.

    Multiallelic sites emit one tuple per indel ALT (SNP alts are skipped), and a record
    whose FILTER is not PASS_FILTERS emits none. ``keep_samples`` (a set) restricts the
    emitted carriers to that panel; None keeps all. Only carriers (normalize_gt_for_alt !=
    None) are stored.
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
                if not _is_indel(ref, alt):
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
                    yield int(pos), ref, alt, carriers


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
    return n


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


class IndelGenotypeReader(object):
    """Read a phased indel genotype store (one chromosome) into a dict keyed by
    ``<pos>_<REF>_<ALT>``, REF and ALT in upper case (store_key). Small enough to hold
    in memory (indels << SNPs). ``reads_filter``: whether the store holds only the
    records of PASS_FILTERS (store_reads_filter); one written before 2.10.0 may hold the
    others too. A key the store holds twice (an indel the VCF writes on two lines, in
    one case or in two) is one allele, its carriers merged (merge_store_carriers)."""

    __slots__ = ("_by_key", "reads_filter")

    def __init__(self, store_path):
        self._by_key = {}
        self.reads_filter = store_reads_filter(store_path)
        opener = gzip.open if store_path.endswith(".gz") else open
        with opener(store_path, "rt") as fh:
            for line in fh:
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
        """Every (pos, REF, ALT) the store holds, sorted by position."""
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


def main(argv):
    """CLI: build_indel_genotypes.py <vcf> <out_store.tsv.gz> [<samplesID>]

    Enabled by default (opt-out): no-op only if CRISPRME_INDEL_SNP=0.
    """
    if os.environ.get("CRISPRME_INDEL_SNP", "1") not in ("1", "true", "True", "yes"):
        return 0
    if len(argv) < 3:
        sys.stderr.write("usage: build_indel_genotypes.py <vcf> <out.tsv.gz> [samplesID]\n")
        return 2
    vcf_path, out_path = argv[1], argv[2]
    keep = _load_panel(argv[3]) if len(argv) > 3 and os.path.isfile(argv[3]) else None
    n = compile_indel_genotypes(vcf_path, out_path, keep_samples=keep)
    print(f"build_indel_genotypes: {n} phased indel records -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

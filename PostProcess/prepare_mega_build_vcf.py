#!/usr/bin/env python3
"""The build input of the sites-only "mega" index, from its merged VCF.

``mega_gw_merge.sh`` (or ``merge_mega_sites.sh``) writes one ``mega.<chrom>.afmax.vcf.gz``
per chromosome: the union of the five sources' sites, each source's frequency as
``AF_<source>`` and their maximum as ``AF_max``, the FILTER of the one source that sets
one (HGDP-2019, ``LOW_VQSLOD`` / ``ExcHet``) carried onto the merged record. CRISPRme 2.10.0
rebuilds the mega from those same VCFs, and every part of the index reads the same records:

  * only the records whose FILTER is ``PASS`` or ``.``, the ones CRISPRitz's enricher
    reads (it skips every other line): the genome, the indel log and fake-indel windows,
    the registry (build_mega_registry.py) and the indel frequencies
    (build_mega_indel_af.py) then hold one set of records;
  * ``INFO/AF`` = ``AF_max``, for SNVs and indels alike: the one global frequency the
    enricher writes into the indel log and the dictionaries (it reads ``INFO/AF`` by its
    key since CRISPRitz 2.10.0), so an indel row shows its frequency, as a SNP row shows
    the registry's GLOBAL group (``AF_max`` too). The as-run merge wrote no ``INFO/AF``,
    so every indel row of the shipped mega showed none. An ``AF`` the input already
    has is replaced. A record without ``AF_max`` gets the one its sources give, the largest
    of its ``AF_<source>`` values per ALT, as ``AF_max`` and ``AF`` both: the as-run merge
    wrote none on chrY (its 177,968 PASS records), and declares it. One with no source
    frequency either gets no ``AF``.

  * one record of each REF>ALT at a position: of an allele on two PASS/'.' records, the
    more frequent record's (its INFO/AF, AF_max; on a tie, the earlier), as the enricher
    merges two SNP records of one allele in its dictionaries. The enricher merges no
    indel records: each gets a fake-indel window and an indel-log row of its own, while the
    indel frequencies (build_mega_indel_af.py), the registry and the checks
    (verify_mega_chrom.py) hold an allele once. So the build input holds one record each,
    and every part reads one. The as-run merge wrote 8 such indels genome-wide (long
    deletions on chr1, chr4, chr7, chr10 (3), chr13 and chr21, each on two lines with one
    AF_max) and no SNV: with both lines, chr21's check failed ("1 allele(s) on more than
    one row"). A duplicated ALT of a line listing several leaves that line, with its
    Number=A and Number=R INFO values (a line with sample columns stops the run: its
    genotypes would need recoding).

Every other byte of a kept record is the input's: its ALTs (every allele of a line that
lists several), ID, QUAL and the other INFO keys, in their order. The header gains the
``AF`` INFO line (where it has none) and a ``##CRISPRme_mega_build=`` line saying what was
done. Counts go to stderr.

Usage:  prepare_mega_build_vcf.py IN.vcf[.gz] OUT.vcf[.gz]
A ``.gz`` output is BGZF (pysam), so ``tabix`` indexes it, as synth_sites_gt.py writes.
STDLIB only, but for the BGZF writer.
"""
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_mega_registry import af_rank  # noqa: E402  (the enricher's rank of two records)

PASS_FILTERS = ("PASS", ".")
AF_HEADER = ('##INFO=<ID=AF,Number=A,Type=Float,Description="Global allele frequency: '
             'the maximum per-source frequency (AF_max) of the merged sources">')
BUILD_MARK = ("##CRISPRme_mega_build=records whose FILTER is PASS or '.' only; "
              "INFO/AF = AF_max; one record of each REF>ALT, the more frequent "
              "(prepare_mega_build_vcf.py)")
_INFO_NUMBER = re.compile(r'^##INFO=<ID=([^,>]+),Number=([^,>]+)')


def _open_read(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _open_write(path):
    """(writer, is_bgzf): a ``.gz`` path -> pysam BGZFile (tabix-indexable); otherwise a
    plain text file ('-' is stdout)."""
    if path == "-":
        return sys.stdout, False
    if path.endswith(".gz"):
        try:
            from pysam import BGZFile
        except ImportError:  # pragma: no cover
            sys.exit("prepare_mega_build_vcf: writing a BGZF .vcf.gz needs pysam "
                     "(install pysam, or give a plain .vcf output path)")
        return BGZFile(path, "wb"), True
    return open(path, "w"), False


def _number(value):
    """A frequency of the VCF as a float, or None ('.', empty, not a number, not positive)."""
    try:
        v = float(value)
    except ValueError:
        return None
    return v if 0.0 < v <= 1.0 else None


def derived_af_max(info, n_alts):
    """The ``AF_max`` a record without one has: per ALT, the largest of its sources'
    frequencies (the ``AF_<source>`` keys, Number=A), as the merge computes it, the value as
    the source wrote it; '.' for an ALT no source gives one; None when no ALT has one. The
    as-run merge declares AF_max on every chromosome and writes it on none of chrY's records
    (its 177,968 PASS records; every other chromosome has it on all)."""
    best = [None] * n_alts
    for kv in (info.split(";") if info and info != "." else []):
        key, _s, value = kv.partition("=")
        if not key.startswith("AF_") or key == "AF_max":
            continue
        vals = value.split(",")
        if len(vals) != n_alts:
            continue
        for k, text in enumerate(vals):
            v = _number(text)
            if v is not None and (best[k] is None or v > best[k][0]):
                best[k] = (v, text)
    if all(b is None for b in best):
        return None
    return ",".join("." if b is None else b[1] for b in best)


def with_af_max(info, n_alts=1):
    """``info`` with its ``AF`` set to its ``AF_max`` value (an ``AF`` it has is dropped,
    the new one appended last), and whether it got one. A record with no ``AF_max`` (or one
    of '.') gets the one its sources give (derived_af_max), written as ``AF_max`` too, so
    every part of the index reads it; one with no source frequency either gets no ``AF``.
    Returns (info, got an AF, derived the AF_max)."""
    keys = info.split(";") if info and info != "." else []
    afmax = None
    for kv in keys:
        if kv.startswith("AF_max="):
            afmax = kv[len("AF_max="):]
    kept = [kv for kv in keys if kv.partition("=")[0] != "AF"]
    derived = False
    if afmax in (None, "", "."):
        afmax = derived_af_max(info, n_alts)
        if afmax is None:
            return (";".join(kept) or "."), False, False
        kept = [kv for kv in kept if kv.partition("=")[0] != "AF_max"] + ["AF_max=" + afmax]
        derived = True
    return ";".join(kept + ["AF=" + afmax]), True, derived


def _without_alts(f, drop, number_of):
    """Record ``f`` without its ALTs of the indices ``drop``, and their values of each INFO
    key the header gives one value an ALT (Number=A) or an allele (Number=R)."""
    if len(f) > 8:
        sys.exit("prepare_mega_build_vcf: %s:%s lists a REF>ALT another record has, on a line "
                 "with sample columns: its genotypes would need recoding" % (f[0], f[1]))
    alts = f[4].split(",")
    keep = [k for k in range(len(alts)) if k not in drop]
    info = []
    for kv in (f[7].split(";") if f[7] not in ("", ".") else []):
        key, sep, value = kv.partition("=")
        values = value.split(",")
        number = number_of.get(key)
        if sep and number == "A" and len(values) == len(alts):
            kv = key + "=" + ",".join(values[k] for k in keep)
        elif sep and number == "R" and len(values) == len(alts) + 1:
            kv = key + "=" + ",".join([values[0]] + [values[k + 1] for k in keep])
        info.append(kv)
    return f[:4] + [",".join(alts[k] for k in keep)] + f[5:7] + [";".join(info) or "."]


def collapse(records, number_of, stats):
    """The records of one position (FILTER PASS or '.', INFO/AF set), one record of each
    REF>ALT: of an allele on several, the more frequent record's (af_rank, the enricher's
    rank of INFO/AF; on a tie, the earlier). Another record of it loses the allele: the
    whole line when it is its only ALT (``stats["collapsed"]``), else that ALT
    (``stats["alts_collapsed"]``)."""
    best = {}
    for i, f in enumerate(records):
        alts = f[4].split(",")
        for k, alt in enumerate(alts):
            rank = af_rank(f[7], k, len(alts))
            have = best.get((f[3], alt))
            if have is None or rank > have[0]:
                best[(f[3], alt)] = (rank, i, k)
    out = []
    for i, f in enumerate(records):
        alts = f[4].split(",")
        drop = {k for k, alt in enumerate(alts) if best[(f[3], alt)][1:] != (i, k)}
        if not drop:
            out.append(f)
        elif len(drop) == len(alts):
            stats["collapsed"] += 1
        else:
            stats["alts_collapsed"] += len(drop)
            out.append(_without_alts(f, drop, number_of))
    return out


def prepare(inp, out):
    """Write the build input of ``inp`` to ``out``; return its counts: ``records`` (in),
    ``kept``, ``filtered`` (another FILTER), ``collapsed`` (a PASS/'.' record whose every
    ALT another record of its position has, more frequent: left out), ``alts_collapsed``
    (such ALTs taken off a line that keeps others), ``af`` (kept records given an AF),
    ``no_af_max`` (kept records with no AF_max and no source frequency, so no AF),
    ``af_max_derived`` (PASS/'.' records whose AF_max is the one their sources give, as the
    VCF wrote none)."""
    stats = {"records": 0, "kept": 0, "filtered": 0, "collapsed": 0, "alts_collapsed": 0,
             "af": 0, "no_af_max": 0, "af_max_derived": 0}
    w, is_bgzf = _open_write(out)

    def emit(s):
        w.write(s.encode() if is_bgzf else s)

    def flush(block):
        for f in collapse(block, number_of, stats):
            stats["kept"] += 1
            stats["af" if "AF_max=" in f[7] and any(
                kv.startswith("AF=") for kv in f[7].split(";")) else "no_af_max"] += 1
            emit("\t".join(f) + "\n")

    has_af_header = False
    number_of = {"AF": "A"}
    block = []   # the PASS/'.' records of one position
    try:
        with _open_read(inp) as fh:
            for line in fh:
                if line.startswith("##"):
                    if line.startswith("##INFO=<ID=AF,"):
                        has_af_header = True
                    m = _INFO_NUMBER.match(line)
                    if m:
                        number_of[m.group(1)] = m.group(2)
                    emit(line)
                    continue
                if line.startswith("#CHROM"):
                    if not has_af_header:
                        emit(AF_HEADER + "\n")
                    emit(BUILD_MARK + "\n")
                    emit(line)
                    continue
                if not line.strip():
                    continue
                stats["records"] += 1
                f = line.rstrip("\n").split("\t")
                if len(f) < 8 or f[6] not in PASS_FILTERS:
                    stats["filtered"] += 1
                    continue
                f[7], _got, derived = with_af_max(f[7], len(f[4].split(",")))
                stats["af_max_derived"] += derived
                if block and (f[0], f[1]) != (block[0][0], block[0][1]):
                    flush(block)
                    block = []
                block.append(f)
            flush(block)
    finally:
        if w is not sys.stdout:
            w.close()
    return stats


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        sys.exit("usage: prepare_mega_build_vcf.py IN.vcf[.gz] OUT.vcf[.gz]")
    st = prepare(argv[0], argv[1])
    sys.stderr.write(
        "prepare_mega_build_vcf: %d records -> %d kept (FILTER PASS or '.'), %d of another "
        "FILTER left out, %d left out as a second record of one REF>ALT (%d ALT(s) taken off "
        "a line that keeps others); INFO/AF = AF_max on %d (%d with the AF_max their sources "
        "give, the VCF writing none), %d with no AF_max -> %s\n"
        % (st["records"], st["kept"], st["filtered"], st["collapsed"], st["alts_collapsed"],
           st["af"], st["af_max_derived"], st["no_af_max"], argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

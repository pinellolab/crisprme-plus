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
    has is replaced. A record without ``AF_max`` gets no ``AF``.

Every other byte of a kept record is the input's: its ALTs (every allele of a line that
lists several), ID, QUAL and the other INFO keys, in their order. The header gains the
``AF`` INFO line (where it has none) and a ``##CRISPRme_mega_build=`` line saying what was
done. Counts go to stderr.

Usage:  prepare_mega_build_vcf.py IN.vcf[.gz] OUT.vcf[.gz]
A ``.gz`` output is BGZF (pysam), so ``tabix`` indexes it, as synth_sites_gt.py writes.
STDLIB only, but for the BGZF writer.
"""
import gzip
import sys

PASS_FILTERS = ("PASS", ".")
AF_HEADER = ('##INFO=<ID=AF,Number=A,Type=Float,Description="Global allele frequency: '
             'the maximum per-source frequency (AF_max) of the merged sources">')
BUILD_MARK = ("##CRISPRme_mega_build=records whose FILTER is PASS or '.' only; "
              "INFO/AF = AF_max (prepare_mega_build_vcf.py)")


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


def with_af_max(info):
    """``info`` with its ``AF`` set to its ``AF_max`` value (an ``AF`` it has is dropped,
    the new one appended last), and whether it got one. No ``AF_max``, or one of '.', gives
    no ``AF``."""
    keys = info.split(";") if info and info != "." else []
    afmax = None
    for kv in keys:
        if kv.startswith("AF_max="):
            afmax = kv[len("AF_max="):]
    kept = [kv for kv in keys if kv.partition("=")[0] != "AF"]
    if afmax in (None, "", "."):
        return (";".join(kept) or "."), False
    return ";".join(kept + ["AF=" + afmax]), True


def prepare(inp, out):
    """Write the build input of ``inp`` to ``out``; return its counts: ``records`` (in),
    ``kept``, ``filtered`` (another FILTER), ``af`` (kept records given an AF),
    ``no_af_max`` (kept records with no AF_max, so no AF)."""
    stats = {"records": 0, "kept": 0, "filtered": 0, "af": 0, "no_af_max": 0}
    w, is_bgzf = _open_write(out)

    def emit(s):
        w.write(s.encode() if is_bgzf else s)

    has_af_header = False
    try:
        with _open_read(inp) as fh:
            for line in fh:
                if line.startswith("##"):
                    if line.startswith("##INFO=<ID=AF,"):
                        has_af_header = True
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
                f[7], got = with_af_max(f[7])
                stats["kept"] += 1
                stats["af" if got else "no_af_max"] += 1
                emit("\t".join(f) + "\n")
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
        "FILTER left out; INFO/AF = AF_max on %d, %d with no AF_max -> %s\n"
        % (st["records"], st["kept"], st["filtered"], st["af"], st["no_af_max"], argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

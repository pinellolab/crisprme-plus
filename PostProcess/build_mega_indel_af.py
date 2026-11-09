#!/usr/bin/env python3
"""Build (and read) the sites-only INDEL allele-frequency sidecar for the "mega"
all-source index.

The Tier-0 registry is SNP-only (it packs single-base ref/alt), so indel off-targets
get no per-dataset AF from it. In the genotyped pipeline an indel's frequency is
embedded in the fake-indel log at `add-variants` time (from the VCF's `INFO/AF`); the
mega's build input carries `INFO/AF` = `AF_max` (seq_script/merge_panels/mega_rebuild.sh),
which gives the log the global frequency, while this sidecar keeps every source's: for
every indel-path allele of the merged VCF it records each source's `AF_<dataset>` plus
`AF_max`, keyed by (pos, ref, alt), so the indel post-analysis and the carrier gate can
annotate an indel off-target with per-dataset AF exactly like the SNP registry does for
SNPs.

The VCF is read exactly as CRISPRitz 2.10.0's enricher reads it
(build_mega_registry.read_as_enricher): only the records whose FILTER is `PASS` or `.`;
the alleles it puts in a fake-indel window (a REF or ALT longer than one base, an MNV
included, no symbolic allele), each ALT of a line listing several with its own AF; one
REF>ALT on two records is one entry, with the frequencies of the more frequent record
(the higher `INFO/AF`; on a tie, the earlier), the rule by which the enricher merges an
allele's entries; and no indel within 26 bases of the chromosome's start, for which the
enricher writes no window. Up to CRISPRme 2.8.1 it read every record whatever its
FILTER, and merged a repeated key by the highest frequency of each source.

Format: one gzipped TSV per chromosome, `indel_af_<chrom>.tsv.gz`, columns
`pos<TAB>ref<TAB>alt<TAB>AF_<ds1>..AF_<dsN><TAB>AF_max` with a `#`-comment header
naming the dataset columns. Missing per-dataset AF is `.`. The reader loads it into a
dict keyed by (pos, ref, alt) for O(1) lookup during post-analysis.

Usage:
  build_mega_indel_af.py --vcf mega.chrN.build.vcf.gz --chrom chrN --out-dir DIR
      [--datasets 1000G2021,HGDP,gnomAD,TOPMed,AoU]
Writes <out-dir>/indel_af_<chrom>.tsv.gz.
"""
import argparse
import gzip
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from build_mega_registry import (  # noqa: E402
    DEFAULT_DATASET_META, _parse_info_afs, indel_alleles, read_as_enricher,
)

__all__ = ["DEFAULT_DATASET_META", "_parse_info_afs", "IndelAfReader",
           "build_indel_af_store", "iter_indel_af_records"]

# The enricher writes no fake-indel window for an indel whose window would start before
# the chromosome (indel_to_fasta: POS - 26 < 0), so no row can name it.
WINDOW_LEFT = 26


def _is_indel(ref, alt):
    """True for an allele CRISPRitz's enricher puts in a fake-indel window: a REF or ALT
    longer than one base, of another length or both longer than one base (an MNV), not
    symbolic (build_mega_registry.indel_alleles). A SNP is not, nor a one-base ALT the
    enricher does not code (a '*' over a one-base REF)."""
    return bool(indel_alleles(ref, [alt]))


def _parse_af_max(info, k=0, n_alts=1):
    """Read INFO/AF_max (already computed by the merge) of ALT ``k`` of ``n_alts`` (a
    record listing several ALTs gives one value per ALT). '.'/missing/garbage, or a
    list of another length -> None."""
    for kv in info.split(";"):
        if kv.startswith("AF_max="):
            values = kv[len("AF_max="):].split(",")
            if len(values) != n_alts:
                return None
            try:
                v = float(values[k])
            except ValueError:
                return None
            return v if (0.0 < v <= 1.0 and v == v) else None
    return None


def iter_indel_af_records(vcf, datasets, chrom=None, stats=None):
    """Yield (pos:str, ref, alt, af_by_dataset:dict, af_max:float|None) for each
    indel-path allele of a sites-only merged VCF, read as CRISPRitz's enricher reads it
    (build_mega_registry.read_as_enricher, kind "indel"): FILTER PASS or '.'; each ALT
    of a line listing several with its own AFs (the k-th of each list); one REF>ALT on
    two records once, with the more frequent record's frequencies; none within
    WINDOW_LEFT bases of the chromosome's start (``stats["no_window"]``). ``chrom``
    (optional) restricts to one contig (defends against a multi-contig VCF whose same-POS
    records on different chroms would otherwise collide by key)."""
    if stats is None:
        stats = {}
    stats.setdefault("no_window", 0)
    for pos, ref, alt, _rsid, afs, _rank, info, k, n in read_as_enricher(
            vcf, "indel", datasets, chrom=chrom, stats=stats):
        if pos - WINDOW_LEFT < 0:
            stats["no_window"] += 1
            continue
        yield str(pos), ref, alt, afs, _parse_af_max(info, k, n)


def build_indel_af_store(vcf, chrom, out_dir, dataset_meta=None, stats=None):
    """Write <out_dir>/indel_af_<chrom>.tsv.gz. Returns (record count, path).

    Atomic: writes a .tmp then os.replace, so an interrupted build never leaves a
    truncated store that a resume would trust (the gt_indel_chr1 lesson). Keys are
    unique (one REF>ALT read once, iter_indel_af_records), so IndelAfReader's dict
    cannot silently drop a record. ``stats`` (a dict): read_as_enricher's counts."""
    dataset_meta = dataset_meta or DEFAULT_DATASET_META
    labels = list(dataset_meta)
    datasets = set(labels)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "indel_af_%s.tsv.gz" % chrom)
    tmp = out + ".tmp"
    if stats is None:
        stats = {}
    # (pos:int, REF, ALT) upper-cased -> ([af|None per label], af_max|None, rank). The
    # reader looks an indel up upper-cased: two records of one indel written in other
    # cases (which the enricher windows apart) are one entry, the more frequent one's.
    acc = {}
    stats.setdefault("no_window", 0)
    for pos, ref, alt, _rsid, afs, rank, info, k, n in read_as_enricher(
            vcf, "indel", datasets, chrom=chrom, stats=stats):
        if pos - WINDOW_LEFT < 0:
            stats["no_window"] += 1
            continue
        key = (int(pos), ref.upper(), alt.upper())
        if key in acc and rank <= acc[key][2]:
            continue
        acc[key] = ([afs.get(d) for d in labels], _parse_af_max(info, k, n), rank)
    try:
        with gzip.open(tmp, "wt") as fh:
            fh.write("#pos\tref\talt\t%s\tAF_max\n"
                     % "\t".join("AF_" + d for d in labels))
            for pos, ref, alt in sorted(acc):
                vals, af_max, _rank = acc[(pos, ref, alt)]
                fh.write("%d\t%s\t%s\t%s\t%s\n" % (
                    pos, ref, alt,
                    "\t".join("." if v is None else repr(v) for v in vals),
                    "." if af_max is None else repr(af_max)))
        os.replace(tmp, out)
    except BaseException:
        if os.path.isfile(tmp):
            os.remove(tmp)
        raise
    if stats.get("merged"):
        sys.stderr.write("  [indel-af] %s: %d second record(s) of one REF>ALT merged into "
                         "the more frequent one\n" % (chrom, stats["merged"]))
    return len(acc), out


class IndelAfReader(object):
    """Loads indel_af_<chrom>.tsv.gz into a dict keyed by (pos:int, ref, alt) for
    O(1) lookup. lookup() returns {dataset: af, 'AF_max': af} (only present keys),
    or None if the indel is absent. ``wanted``: the positions (POS) to keep, so a
    chromosome's whole store is not held for a few indels; None keeps all."""

    def __init__(self, store_path, wanted=None):
        self._d = {}
        self._labels = []
        opener = gzip.open if store_path.endswith(".gz") else open
        with opener(store_path, "rt") as fh:
            header = fh.readline().rstrip("\n")
            # "#pos ref alt AF_<ds>.. AF_max" -> the dataset labels between alt and AF_max.
            cols = header.lstrip("#").split("\t")
            if len(cols) < 4 or cols[0] != "pos" or cols[-1] != "AF_max":
                raise ValueError("indel-AF store: unexpected header %r" % header)
            self._labels = [c[len("AF_"):] for c in cols[3:-1]]
            nfields = len(cols)  # pos, ref, alt, <labels...>, AF_max
            afm_col = 3 + len(self._labels)  # AF_max column by POSITION (not f[-1])
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) != nfields:
                    continue  # skip a malformed row rather than mis-parse AF_max
                pos, ref, alt = int(f[0]), f[1], f[2]
                if wanted is not None and pos not in wanted:
                    continue
                rec = {}
                for i, lab in enumerate(self._labels):
                    v = f[3 + i]
                    if v != ".":
                        rec[lab] = float(v)
                afm = f[afm_col]
                if afm != ".":
                    rec["AF_max"] = float(afm)
                self._d[(pos, ref.upper(), alt.upper())] = rec

    def lookup(self, pos, ref, alt):
        return self._d.get((int(pos), ref.upper(), alt.upper()))

    def __len__(self):
        return len(self._d)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--datasets", default=",".join(DEFAULT_DATASET_META))
    args = ap.parse_args(argv)
    labels = [d.strip() for d in args.datasets.split(",") if d.strip()]
    meta = {d: DEFAULT_DATASET_META[d] for d in labels if d in DEFAULT_DATASET_META}
    stats = {}
    n, out = build_indel_af_store(args.vcf, args.chrom, args.out_dir,
                                  dataset_meta=meta, stats=stats)
    sys.stderr.write("mega indel-AF store %s: %d indels -> %s (%d alleles of records "
                     "outside FILTER PASS/. skipped; %d with no AF; %d within %d bases of "
                     "the start, which no window holds)\n"
                     % (args.chrom, n, out, stats.get("filtered", 0), stats.get("no_af", 0),
                        stats.get("no_window", 0), WINDOW_LEFT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

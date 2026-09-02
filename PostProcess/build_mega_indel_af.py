#!/usr/bin/env python3
"""Build (and read) the sites-only INDEL allele-frequency sidecar for the "mega"
all-source index.

The Tier-0 registry is SNP-only (it packs single-base ref/alt), so indel off-targets
get no per-dataset AF from it. In the genotyped pipeline an indel's frequency is
embedded in the fake-indel log at `add-variants` time (from the VCF's `INFO/AF` +
genotypes); the sites-only mega VCF has neither a plain `AF` nor genotypes, so that
path yields nothing. This sidecar closes the gap: for every INDEL in the merged
`mega.<chrom>.afmax.vcf.gz` it records each source's `AF_<dataset>` plus `AF_max`,
keyed by (pos, ref, alt), so the indel post-analysis can annotate an indel off-target
with per-dataset AF exactly like the SNP registry does for SNPs.

Format: one gzipped TSV per chromosome, `indel_af_<chrom>.tsv.gz`, columns
`pos<TAB>ref<TAB>alt<TAB>AF_<ds1>..AF_<dsN><TAB>AF_max` with a `#`-comment header
naming the dataset columns. Missing per-dataset AF is `.`. The reader loads it into a
dict keyed by (pos, ref, alt) for O(1) lookup during post-analysis.

Usage:
  build_mega_indel_af.py --vcf mega.chrN.afmax.vcf.gz --chrom chrN --out-dir DIR
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
    DEFAULT_DATASET_META, _parse_info_afs,
)


def _is_indel(ref, alt):
    """True for a normalized biallelic INDEL/MNV (anything not a single-base SNP).
    Mirrors build_mega_registry's SNP gate, complemented."""
    ru, au = ref.upper(), alt.upper()
    snp = len(ru) == 1 and len(au) == 1 and ru in "ACGT" and au in "ACGT"
    return not snp


def _parse_af_max(info):
    """Read INFO/AF_max (already computed by the merge). '.'/missing/garbage -> None."""
    for kv in info.split(";"):
        if kv.startswith("AF_max="):
            try:
                v = float(kv[len("AF_max="):])
            except ValueError:
                return None
            return v if (0.0 < v <= 1.0 and v == v) else None
    return None


def iter_indel_af_records(vcf, datasets):
    """Yield (pos:str, ref, alt, af_by_dataset:dict, af_max:float|None) for each
    biallelic INDEL record of a sites-only merged VCF. Comma-ALT rows are skipped."""
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            _c, pos, _vid, ref, alt, _q, _fl, info = f[:8]
            if "," in alt or not _is_indel(ref, alt):
                continue
            afs = _parse_info_afs(info, datasets)
            if not afs:
                continue
            yield (pos, ref, alt, afs, _parse_af_max(info))


def build_indel_af_store(vcf, chrom, out_dir, dataset_meta=None):
    """Write <out_dir>/indel_af_<chrom>.tsv.gz. Returns the record count.

    Atomic: writes a .tmp then os.replace, so an interrupted build never leaves a
    truncated store that a resume would trust (the gt_indel_chr1 lesson)."""
    dataset_meta = dataset_meta or DEFAULT_DATASET_META
    labels = list(dataset_meta)
    datasets = set(labels)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "indel_af_%s.tsv.gz" % chrom)
    tmp = out + ".tmp"
    n = 0
    try:
        with gzip.open(tmp, "wt") as fh:
            fh.write("#pos\tref\talt\t%s\tAF_max\n"
                     % "\t".join("AF_" + d for d in labels))
            for pos, ref, alt, afs, af_max in iter_indel_af_records(vcf, datasets):
                cols = [afs.get(d) for d in labels]
                fh.write("%s\t%s\t%s\t%s\t%s\n" % (
                    pos, ref.upper(), alt.upper(),
                    "\t".join("." if c is None else repr(c) for c in cols),
                    "." if af_max is None else repr(af_max)))
                n += 1
        os.replace(tmp, out)
    except BaseException:
        if os.path.isfile(tmp):
            os.remove(tmp)
        raise
    return n, out


class IndelAfReader(object):
    """Loads indel_af_<chrom>.tsv.gz into a dict keyed by (pos:int, ref, alt) for
    O(1) lookup. lookup() returns {dataset: af, 'AF_max': af} (only present keys),
    or None if the indel is absent."""

    def __init__(self, store_path):
        self._d = {}
        self._labels = []
        opener = gzip.open if store_path.endswith(".gz") else open
        with opener(store_path, "rt") as fh:
            header = fh.readline().rstrip("\n")
            # "#pos ref alt AF_<ds>.. AF_max" -> the dataset labels between alt and AF_max
            cols = header.lstrip("#").split("\t")
            self._labels = [c[len("AF_"):] for c in cols[3:-1]]
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                pos, ref, alt = int(f[0]), f[1], f[2]
                rec = {}
                for i, lab in enumerate(self._labels):
                    v = f[3 + i]
                    if v != ".":
                        rec[lab] = float(v)
                afm = f[-1]
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
    n, out = build_indel_af_store(args.vcf, args.chrom, args.out_dir,
                                  dataset_meta=meta)
    sys.stderr.write("mega indel-AF store %s: %d indels -> %s\n"
                     % (args.chrom, n, out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build the Tier-0 registry for the sites-only "mega" all-source merged index.

The mega index is a merged, genotype-stripped VCF whose INFO carries per-dataset
allele frequencies (``AF_<dataset>``) plus ``AF_max``. This driver reads those AFs
directly and compiles the mmap registry via
``tier0_registry.compile_registry_from_info_af`` -- NO genotypes required (unlike
the genotype-counting ``tier0_compile`` path used by the per-sample panels).

The VCF is parsed with a dependency-free streaming gzip reader (only the 8 fixed
columns + INFO are needed), so this module imports and unit-tests without pysam.

Usage:
  build_mega_registry.py --vcf mega.chrN.afmax.vcf.gz --chrom chrN --out-dir DIR
      [--datasets 1000G2021,HGDP,gnomAD,TOPMed,AoU] [--compress]

Writes ``<out-dir>/reg_<chrom>.bin`` + ``<out-dir>/reg_<chrom>.idx`` (the per-chrom
registry layout the dict-less search reads).
"""
import argparse
import gzip
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from tier0_registry import compile_registry_from_info_af  # noqa: E402

# The five mega datasets and their nominal allele numbers (2 x documented cohort
# N). AN_nominal is the denominator SHOWN for that dataset and the precision of
# the reconstructed AF (AC = round(AF x AN)); using each cohort's true 2N keeps
# both meaningful. Override cohort sizes here if a dataset's N is re-snapshotted.
DEFAULT_DATASET_META = {
    "1000G2021": {"sample_count": 3202, "an_nominal": 6404},
    "HGDP": {"sample_count": 929, "an_nominal": 1858},
    "gnomAD": {"sample_count": 807162, "an_nominal": 1614324},
    "TOPMed": {"sample_count": 53831, "an_nominal": 107662},
    "AoU": {"sample_count": 535662, "an_nominal": 1071324},
}


def _parse_info_afs(info, datasets):
    """Extract {dataset -> AF_float} from a VCF INFO string. AF_max is skipped
    (the registry recomputes the global from the per-dataset AFs)."""
    afs = {}
    for kv in info.split(";"):
        eq = kv.find("=")
        if eq < 0:
            continue
        key = kv[:eq]
        if not key.startswith("AF_") or key == "AF_max":
            continue
        ds = key[3:]
        if ds not in datasets:
            continue
        try:
            val = float(kv[eq + 1:])
        except ValueError:
            continue
        # Accept only a real frequency in (0, 1]. This rejects a comma-list
        # ("0.01,0.02" -> ValueError above), "." (missing), 0.0 (monomorphic), and
        # -- critically -- inf / a 1e400-style overflow literal (float("inf") > 0.0
        # is True, and would later OverflowError in round(inf*AN), aborting the whole
        # chromosome build). AF > 1 (garbage) is dropped too, so AC can never exceed AN.
        if 0.0 < val <= 1.0 and val == val:  # (val==val excludes NaN defensively)
            afs[ds] = val
    return afs


def iter_vcf_af_records(path, datasets, chrom=None, stats=None):
    """Yield ``(pos, ref, alt, rsid, af_by_dataset)`` for each biallelic record of
    a sites-only merged VCF. Records must be biallelic (the merge splits with
    ``bcftools norm -m -any``); a multiallelic ALT is skipped defensively.

    ``datasets`` is the set of dataset labels to read (the INFO ``AF_<label>``
    fields). ``chrom`` (optional) restricts to one contig. ``stats`` (optional dict)
    accumulates drop counts so nothing is silently discarded: ``multiallelic`` (an
    un-split comma ALT) and ``no_af`` (no parseable AF for any known dataset).
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            c, pos, vid, ref, alt, _qual, _filt, info = f[:8]
            if chrom is not None and c != chrom:
                continue
            if "," in alt:  # not split to biallelic -- skip (should not happen)
                if stats is not None:
                    stats["multiallelic"] = stats.get("multiallelic", 0) + 1
                continue
            afs = _parse_info_afs(info, datasets)
            if not afs:  # no parseable AF for any known dataset at this site
                if stats is not None:
                    stats["no_af"] = stats.get("no_af", 0) + 1
                continue
            rsid = vid if (vid and vid != ".") else ""
            yield (int(pos), ref, alt, rsid, afs)


def build(vcf, chrom, out_dir, dataset_meta=None, compress=False):
    """Build ``reg_<chrom>.bin``/``.idx`` in ``out_dir`` from the merged VCF.

    The Tier-0 registry is SNP-only (single-char ref/alt); INDEL records go
    through CRISPRme's separate fake-indel path, not this registry. So SNPs are
    compiled here and indels are COUNTED and skipped (never silently dropped --
    the count is returned in ``stats`` and logged by ``main``). Aggregate-indel
    AF is a follow-up path (the mega ships SNP AF first).

    Returns ``(manifest, bin_path, idx_path, stats)`` where ``stats`` accounts for
    EVERY input record: ``snps`` (compiled), ``indels`` (routed to fake-indel path),
    ``multiallelic`` (un-split comma ALT), ``no_af`` (no parseable AF).
    """
    dataset_meta = dataset_meta or DEFAULT_DATASET_META
    datasets = set(dataset_meta)
    os.makedirs(out_dir, exist_ok=True)
    binp = os.path.join(out_dir, "reg_%s.bin" % chrom)
    idxp = os.path.join(out_dir, "reg_%s.idx" % chrom)
    stats = {"snps": 0, "indels": 0, "multiallelic": 0, "no_af": 0}

    def _snp_records():
        for (pos, ref, alt, rsid, afs) in iter_vcf_af_records(vcf, datasets,
                                                              stats=stats):
            ref_u, alt_u = ref.upper(), alt.upper()
            if len(ref_u) == 1 and len(alt_u) == 1 and ref_u in "ACGT" and alt_u in "ACGT":
                stats["snps"] += 1
                yield (pos, ref_u, alt_u, rsid, afs)
            else:
                stats["indels"] += 1

    manifest = compile_registry_from_info_af(
        _snp_records(), dataset_meta, binp, idxp, compress=compress)
    return manifest, binp, idxp, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, help="merged sites-only VCF (AF_<ds>)")
    ap.add_argument("--chrom", required=True, help="contig label, e.g. chr22")
    ap.add_argument("--out-dir", required=True, help="registry output directory")
    ap.add_argument("--datasets", default=",".join(DEFAULT_DATASET_META),
                    help="comma-joined dataset labels (default: all five). NOTE: the "
                    "GLOBAL group's AF is recomputed as the max over the SELECTED "
                    "datasets, so a subset here yields a GLOBAL that can differ from "
                    "the VCF's INFO/AF_max (which is over all datasets present).")
    ap.add_argument("--compress", action="store_true",
                    help="write the v3 block-compressed registry")
    args = ap.parse_args(argv)
    labels = [d.strip() for d in args.datasets.split(",") if d.strip()]
    meta = {d: DEFAULT_DATASET_META[d] for d in labels if d in DEFAULT_DATASET_META}
    missing = [d for d in labels if d not in DEFAULT_DATASET_META]
    if missing:
        sys.stderr.write("WARNING: unknown dataset label(s) ignored: %s\n"
                         % ",".join(missing))
    manifest, binp, idxp, stats = build(args.vcf, args.chrom, args.out_dir,
                                        dataset_meta=meta, compress=args.compress)
    sys.stderr.write(
        "mega registry %s: %d SNP records compiled | %d indels -> fake-indel path"
        " | %d multiallelic skipped | %d no-AF skipped | %d datasets | "
        "aggregation=%s -> %s\n"
        % (args.chrom, stats["snps"], stats["indels"], stats["multiallelic"],
           stats["no_af"], len(meta), manifest["aggregation"], binp))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

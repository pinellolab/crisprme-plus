#!/usr/bin/env python3
"""Build the Tier-0 registry for the sites-only "mega" all-source merged index.

The mega index is a merged, genotype-stripped VCF whose INFO carries per-dataset
allele frequencies (``AF_<dataset>``) plus ``AF_max``. This driver reads those AFs
directly and compiles the mmap registry via
``tier0_registry.compile_registry_from_info_af`` -- NO genotypes required (unlike
the genotype-counting ``tier0_compile`` path used by the per-sample panels).

It reads the VCF's records exactly as CRISPRitz 2.10.0's enricher (``add-variants``)
reads them, so that the registry holds the SNP alleles the enriched genome codes, and
nothing else (``read_as_enricher``):

  * a record whose FILTER is neither ``PASS`` nor ``.`` is skipped whole (the enricher
    skips the line), its ALTs counted (``filtered``);
  * a SNP allele is an ALT that is one upper-case A, C, G or T of a record whose REF is
    one base: the bases the enricher codes (its ``_BASE_BIT``). Every such ALT of a
    record listing several is read, each with its own AF (``AF_<dataset>`` is
    Number=A) and ID. A one-character ALT it does not code -- the ``*`` of a spanning
    deletion, an ``N``, a lower-case base, a symbolic allele -- is in no store
    (``not_coded``);
  * an indel-path ALT (a REF or ALT longer than one base: an indel, or an MNV, which the
    enricher puts in a fake-indel window) goes to the fake-indel path and to the indel
    frequencies of ``build_mega_indel_af.py`` (``indels``);
  * one REF>ALT on two records is one allele, as the enricher merges it into one
    dictionary entry: the entry takes the frequencies (and ID) of the more frequent
    record -- the higher ``INFO/AF``, found by its key (exactly ``AF``), a missing or
    non-numeric one counting 0; on a tie, the earlier record (``merged``).

Up to CRISPRme 2.8.1 the registry read every record of the VCF whatever its FILTER,
and kept a one-base lower-case ALT (upper-cased): the shipped mega's registry lists
2,442,989 SNVs outside PASS that its genome never coded, and its SNP rows named some.

The VCF is parsed with a dependency-free streaming gzip reader (only the 8 fixed
columns + INFO are needed), so this module imports and unit-tests without pysam.

Usage:
  build_mega_registry.py --vcf mega.chrN.build.vcf.gz --chrom chrN --out-dir DIR
      [--datasets 1000G2021,HGDP,gnomAD,TOPMed,AoU] [--compress]

Writes ``<out-dir>/reg_<chrom>.bin`` + ``<out-dir>/reg_<chrom>.idx`` (the per-chrom
registry layout the dict-less search reads).
"""
import argparse
import gzip
import os
import re
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


def _parse_info_afs(info, datasets, k=0, n_alts=1):
    """Extract {dataset -> AF_float} of ALT ``k`` (of ``n_alts``) from a VCF INFO
    string. AF_max is skipped (the registry recomputes the global from the
    per-dataset AFs). A record listing several ALTs gives one AF per ALT, in their
    order (Number=A); a list of another length, or a list on a record of one ALT,
    gives that dataset no AF."""
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
        values = kv[eq + 1:].split(",")
        if len(values) != n_alts:
            continue
        try:
            val = float(values[k])
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


def _alt_ids(vid, n_alts):
    """The ID of each ALT of a record: one each when the ID column lists one per ALT
    (';'-separated, in the ALTs' order), else the whole column for every ALT (one ID
    covering the site, or IDs that cannot be told apart); '' for none."""
    if not vid or vid == ".":
        return [""] * n_alts
    ids = vid.split(";")
    if n_alts > 1 and len(ids) == n_alts:
        return ["" if i == "." else i for i in ids]
    return [vid] * n_alts


# --------------------------------------------------------------------------- #
# The records as CRISPRitz 2.10.0's enricher reads them (Enrichment/enricher.py and
# enricher.cpp, both ports alike; its docs/MULTIALLELIC_ENRICHMENT.md sections 2.1-2.5)
# --------------------------------------------------------------------------- #
# The FILTER values of the lines the enricher reads: any other line is skipped whole
# (`if line[6] not in ('PASS', '.')`), so nothing of it is coded, listed or windowed.
PASS_FILTERS = ("PASS", ".")
# The SNP ALTs the enricher codes into the genome: one upper-case base (its _BASE_BIT).
# A one-character ALT that is not one of them (the '*' of a spanning deletion, an 'N',
# a lower-case base) adds nothing to the code.
CODED_BASES = ("A", "C", "G", "T")
# A plain decimal number, the only AF string the enricher ranks (its _AF_NUMBER): 'nan',
# 'inf', '1_000' or a list are not one and rank 0.
_AF_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def snp_alleles(ref, alts):
    """[(k, ALT)] of the SNP alleles the enricher codes on a record of REF ``ref`` and
    ALTs ``alts`` (k: the ALT's index on the line, from 0): a REF of one base, an ALT of
    one upper-case A, C, G or T (SNPsProcess: a longer REF codes nothing)."""
    if len(ref) != 1:
        return []
    return [(k, a) for k, a in enumerate(alts) if a in CODED_BASES]


def indel_alleles(ref, alts):
    """[(k, ALT)] of the alleles the enricher puts in a fake-indel window (indel_to_fasta):
    on a record whose REF or ALT is longer than one base and whose REF holds no '<', an
    ALT of another length than REF, or of the same length when both are longer than one
    base (an MNV), that is not symbolic (holds no '<')."""
    if not (len(ref) > 1 or any(len(a) > 1 for a in alts)) or "<" in ref:
        return []
    return [(k, a) for k, a in enumerate(alts)
            if (len(a) != len(ref) or (len(a) > 1 and len(ref) > 1)) and "<" not in a]


def af_rank(info, k=0, n_alts=1):
    """The rank of ALT ``k``'s INFO/AF, as the enricher tells the more frequent of two
    records of one allele apart (_af_rank over _af_value): the value of the key exactly
    ``AF`` (an ``AF_max`` or ``AF_afr`` is not it), the k-th of a line listing several
    ALTs, the whole value on a line of one; 0 when it is missing or not a plain decimal
    number."""
    for kv in info.split(";"):
        key, _sep, value = kv.partition("=")
        if key == "AF":
            if n_alts > 1:
                values = value.split(",")
                value = values[k] if k < len(values) else ""
            return float(value) if _AF_NUMBER.match(value) else 0.0
    return 0.0


def iter_records(path, chrom=None):
    """Yield ``(pos, id, ref, [alts], filter, info)`` for each data line of a VCF
    (plain or gzipped; the 8 fixed columns are all that is read). ``chrom`` (optional)
    restricts to one contig."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            c, pos, vid, ref, alt, _qual, filt, info = f[:8]
            if chrom is not None and c != chrom:
                continue
            yield int(pos), vid, ref, alt.split(","), filt, info


def read_as_enricher(path, kind, datasets, chrom=None, stats=None):
    """Yield ``(pos, ref, alt, rsid, af_by_dataset, af_rank, info, k, n)`` of each allele
    of ``kind``
    ("snp": the SNP alleles the enricher codes; "indel": those it puts in a fake-indel
    window) of a sites-only merged VCF, ONE per REF>ALT at a position: of an allele on
    two records, the more frequent record's (af_rank; on a tie, the earlier), as the
    enricher merges one allele's entries; ``info``, ``k`` and ``n`` are that record's
    INFO, the allele's index on its line and the line's ALT count. In file order of each
    allele's first record.

    ``stats`` (a dict) counts, by ALT, every allele of the VCF, so nothing is silently
    dropped: ``filtered`` (a record of another FILTER than PASS_FILTERS), ``snps`` /
    ``indels`` (the alleles of each kind read), ``not_coded`` (a one-character ALT the
    enricher does not code, or a symbolic one), ``no_af`` (an allele of ``kind`` with
    no parseable AF of a dataset of ``datasets``: in no store), ``merged`` (a second
    record of one REF>ALT, merged into the first), and ``multiallelic`` (the records of
    several ALTs, each ALT read)."""
    if stats is None:
        stats = {}
    for key in ("filtered", "snps", "indels", "not_coded", "no_af", "merged",
                "multiallelic"):
        stats.setdefault(key, 0)
    pick = snp_alleles if kind == "snp" else indel_alleles
    kept = {}
    for pos, vid, ref, alts, filt, info in iter_records(path, chrom):
        n = len(alts)
        if n > 1:
            stats["multiallelic"] += 1
        if filt not in PASS_FILTERS:
            stats["filtered"] += n
            continue
        snps = dict(snp_alleles(ref, alts))
        indels = dict(indel_alleles(ref, alts))
        stats["snps"] += len(snps)
        stats["indels"] += len(indels)
        stats["not_coded"] += n - len(snps) - len(indels)
        ids = _alt_ids(vid, n)
        for k, alt in pick(ref, alts):
            afs = _parse_info_afs(info, datasets, k, n)
            if not afs:  # no parseable AF for any known dataset at this allele
                stats["no_af"] += 1
                continue
            key = (pos, ref, alt)
            rank = af_rank(info, k, n)
            have = kept.get(key)
            if have is not None:
                stats["merged"] += 1
                if rank <= have[5]:
                    continue
            kept[key] = (pos, ref, alt, ids[k], afs, rank, info, k, n)
    for rec in kept.values():
        yield rec


def iter_vcf_af_records(path, datasets, chrom=None, stats=None):
    """Yield ``(pos, ref, alt, rsid, af_by_dataset)`` for each SNP allele of a sites-only
    merged VCF that CRISPRitz's enricher codes, read as it reads them
    (read_as_enricher): FILTER PASS or '.', an upper-case A, C, G or T over a REF of
    one base, one REF>ALT on two records merged into the more frequent one. Each ALT of
    a record listing several is read with its own AF (the k-th of each ``AF_<label>``
    list) and ID (``_alt_ids``). ``stats``: read_as_enricher's counts."""
    for rec in read_as_enricher(path, "snp", datasets, chrom=chrom, stats=stats):
        yield rec[:5]


def build(vcf, chrom, out_dir, dataset_meta=None, compress=False):
    """Build ``reg_<chrom>.bin``/``.idx`` in ``out_dir`` from the merged VCF.

    The Tier-0 registry is SNP-only (single-char ref/alt): it holds the SNP alleles
    CRISPRitz 2.10.0's enricher codes, read as it reads them (read_as_enricher); indel
    records go through CRISPRme's separate fake-indel path and the indel frequencies of
    build_mega_indel_af.py.

    Returns ``(manifest, bin_path, idx_path, stats)`` where ``stats`` accounts for
    EVERY input ALT (read_as_enricher's counts): ``snps`` (SNP alleles read; those with
    an AF are compiled, a second record of one allele merged into the first: ``merged``),
    ``indels`` (to the fake-indel path), ``not_coded``, ``filtered``, ``no_af``; and
    ``written``, the registry's records."""
    dataset_meta = dataset_meta or DEFAULT_DATASET_META
    datasets = set(dataset_meta)
    os.makedirs(out_dir, exist_ok=True)
    binp = os.path.join(out_dir, "reg_%s.bin" % chrom)
    idxp = os.path.join(out_dir, "reg_%s.idx" % chrom)
    stats = {}
    # chrom= filters to the target contig: the Tier-0 registry is keyed by (pos, alt)
    # with NO chromosome, so a multi-contig VCF would collide same-POS records across
    # chroms. The per-chrom mega VCFs carry one contig, so this is a safety guard.
    records = list(iter_vcf_af_records(vcf, datasets, chrom=chrom, stats=stats))
    stats["written"] = len(records)
    manifest = compile_registry_from_info_af(
        records, dataset_meta, binp, idxp, compress=compress,
        manifest_extra={"vcf_records": VCF_RECORDS_READ})
    return manifest, binp, idxp, stats


# The manifest's account of which VCF records the registry holds (a registry of CRISPRme
# 2.8.1 or earlier, which read every record, has no such key).
VCF_RECORDS_READ = ("FILTER PASS or '.'; the SNP ALTs CRISPRitz 2.10.0 add-variants codes "
                    "(an upper-case A, C, G or T over a REF of one base); one REF>ALT on two "
                    "records merged into the more frequent one (INFO/AF)")


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
        "mega registry %s: %d SNP records compiled (%d SNP alleles read, %d of them a "
        "second record of one REF>ALT merged into the more frequent, %d with no AF) | "
        "%d indel-path alleles -> fake-indel path | %d alleles of records outside FILTER "
        "PASS/. skipped | %d one-base ALTs CRISPRitz does not code skipped | %d records of "
        "several ALTs read per ALT | %d datasets | aggregation=%s -> %s\n"
        % (args.chrom, stats["written"], stats["snps"], stats["merged"], stats["no_af"],
           stats["indels"], stats["filtered"], stats["not_coded"], stats["multiallelic"],
           len(meta), manifest["aggregation"], binp))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build a CRISPRme-compatible VCF from an All of Us per-chromosome CSV.

Converts the combined AoU CSV (from ``json_to_csv.py``) into a bgzipped
``VCFv4.2`` with a single aggregate pseudo-sample ``AllOfUs`` and four INFO
fields (``AF``/``AC``/``AN``/``HOM``). Columns are selected **by name**, so the
converter is robust to changes in column order. Records are written in input
order; sort the result downstream with ``bcftools sort``.

Usage
-----
    csv_to_vcf.py --csv csv/chr22.csv --out-vcf chr22.unsorted.vcf.gz --reference hg38.fa
"""
from __future__ import annotations

import argparse
import csv as csvmod
import sys
from typing import Optional, Tuple

import pysam

REQUIRED_COLUMNS = (
    "variantId",
    "alleleCount",
    "alleleNumber",
    "alleleFrequency",
    "homozygoteCount",
)


def parse_variant_id(variant_id: str) -> Optional[Tuple[str, int, str, str]]:
    """Parse an AoU ``variantId`` of the form ``chrom-pos-ref-alt``.

    Parameters
    ----------
    variant_id : str
        Variant identifier, e.g. ``"22-10510077-C-T"`` (chromosome unprefixed).

    Returns
    -------
    tuple of (str, int, str, str) or None
        ``(chrom, pos, ref, alt)`` with a ``chr``-prefixed chromosome, or
        ``None`` if the identifier is malformed.
    """
    parts = variant_id.split("-")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    try:
        return chrom, int(pos), ref, alt
    except ValueError:
        return None


def safe_int(value: str, default: int = 0) -> int:
    """Return *value* as an int, or *default* if empty/invalid.

    Parameters
    ----------
    value : str
        Raw string from the CSV.
    default : int, optional
        Fallback value (default 0).

    Returns
    -------
    int
        Parsed integer or the default.
    """
    try:
        return int(value) if str(value).strip() else default
    except (ValueError, TypeError):
        return default


def safe_float(value: str, default: float = 0.0) -> float:
    """Return *value* as a float, or *default* if empty/invalid.

    Parameters
    ----------
    value : str
        Raw string from the CSV.
    default : float, optional
        Fallback value (default 0.0).

    Returns
    -------
    float
        Parsed float or the default.
    """
    try:
        return float(value) if str(value).strip() else default
    except (ValueError, TypeError):
        return default


def build_header(fasta: "pysam.FastaFile") -> "pysam.VariantHeader":
    """Build the minimal aggregate VCF header (AF/AC/AN/HOM, one pseudo-sample).

    Parameters
    ----------
    fasta : pysam.FastaFile
        Open reference; its contigs (with lengths) are copied into the header.

    Returns
    -------
    pysam.VariantHeader
        Header declaring the four INFO fields, ``PASS`` filter, ``GT`` format,
        every reference contig, and a single sample ``AllOfUs``.
    """
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.2")
    header.add_meta("INFO", items=[("ID", "AF"), ("Number", "A"), ("Type", "Float"), ("Description", "Allele frequency")])
    header.add_meta("INFO", items=[("ID", "AC"), ("Number", "A"), ("Type", "Integer"), ("Description", "Allele count")])
    header.add_meta("INFO", items=[("ID", "AN"), ("Number", "1"), ("Type", "Integer"), ("Description", "Allele number")])
    header.add_meta("INFO", items=[("ID", "HOM"), ("Number", "1"), ("Type", "Integer"), ("Description", "Homozygote count")])
    header.add_meta("FILTER", items=[("ID", "PASS"), ("Description", "All filters passed")])
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")])
    for contig in fasta.references:
        header.contigs.add(contig, length=fasta.get_reference_length(contig))
    header.add_sample("AllOfUs")
    return header


def convert(csv_path: str, out_vcf: str, reference_fasta: str) -> int:
    """Convert an AoU CSV into a bgzipped VCF and return the record count.

    Parameters
    ----------
    csv_path : str
        Input CSV with the AoU field names (see :data:`REQUIRED_COLUMNS`).
    out_vcf : str
        Output path; written bgzip-compressed (``.vcf.gz``), unsorted.
    reference_fasta : str
        Reference FASTA providing the contig set/lengths.

    Returns
    -------
    int
        Number of variant records written.

    Raises
    ------
    KeyError
        If the CSV is missing a required column.
    """
    fasta = pysam.FastaFile(reference_fasta)
    contigs = set(fasta.references)
    header = build_header(fasta)

    written = 0
    skipped = 0
    with pysam.VariantFile(out_vcf, "wz", header=header) as out, open(csv_path, newline="") as handle:
        reader = csvmod.DictReader(handle)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise KeyError(f"CSV {csv_path} is missing columns: {missing}")

        for row in reader:
            parsed = parse_variant_id(row["variantId"])
            if parsed is None:
                skipped += 1
                continue
            chrom, pos, ref, alt = parsed
            if chrom not in contigs:
                skipped += 1
                continue

            ac = safe_int(row["alleleCount"])
            an = safe_int(row["alleleNumber"])
            af = safe_float(row["alleleFrequency"])
            hom = safe_int(row["homozygoteCount"])

            record = out.new_record(
                contig=chrom,
                start=pos - 1,  # pysam is 0-based
                stop=pos - 1 + len(ref),
                alleles=(ref, alt),
                info={"AF": af, "AC": ac, "AN": an, "HOM": hom},
            )
            record.filter.add("PASS")
            # Synthesised genotype for the single aggregate pseudo-sample:
            # homozygous-alt only when every allele is ALT and homozygotes exist.
            record.samples["AllOfUs"]["GT"] = (1, 1) if hom > 0 and ac == an else (0, 1)
            out.write(record)
            written += 1

    if skipped:
        sys.stderr.write(f"{csv_path}: skipped {skipped} malformed/off-contig rows\n")
    return written


def main() -> None:
    """Parse arguments and run the CSV-to-VCF conversion."""
    parser = argparse.ArgumentParser(
        description="Build an AoU aggregate VCF (AllOfUs pseudo-sample) from a CSV."
    )
    parser.add_argument("--csv", required=True, help="Input per-chromosome CSV")
    parser.add_argument("--out-vcf", required=True, help="Output bgzipped VCF (unsorted)")
    parser.add_argument("--reference", required=True, help="Reference FASTA (hg38)")
    args = parser.parse_args()

    n = convert(args.csv, args.out_vcf, args.reference)
    print(f"wrote {n} records -> {args.out_vcf}")


if __name__ == "__main__":
    main()

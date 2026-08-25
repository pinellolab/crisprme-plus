#!/usr/bin/env python
"""Tests for build_indel_genotypes: phase normalization + VCF parse + reader.

Run: python PostProcess/test_indel_genotypes.py   (or via pytest)
"""

import os
import gzip
import tempfile

import build_indel_genotypes as big


def test_normalize_biallelic_phased():
    assert big.normalize_gt_for_alt("1|0", 1) == "1|0"
    assert big.normalize_gt_for_alt("0|1", 1) == "0|1"
    assert big.normalize_gt_for_alt("1|1", 1) == "1|1"
    assert big.normalize_gt_for_alt("0|0", 1) is None  # not a carrier


def test_normalize_multiallelic():
    # indel is ALT #2 (e.g. ALT="A,AT"); slot carrying '2' -> '1'
    assert big.normalize_gt_for_alt("1|2", 2) == "0|1"
    assert big.normalize_gt_for_alt("2|1", 2) == "1|0"
    assert big.normalize_gt_for_alt("2|2", 2) == "1|1"
    assert big.normalize_gt_for_alt("1|0", 2) is None  # carries alt1, not alt2


def test_normalize_unphased_preserved():
    # a '/' separator must survive so the cis check downgrades to PUTATIVE
    assert big.normalize_gt_for_alt("0/1", 1) == "0/1"
    assert big.normalize_gt_for_alt("1/2", 2) == "0/1"
    assert big.normalize_gt_for_alt("0/0", 1) is None


def test_normalize_strips_format_subfields_and_haploid():
    assert big.normalize_gt_for_alt("1|0:0.99:35", 1) == "1|0"  # GT is first subfield
    assert big.normalize_gt_for_alt("1", 1) == "1"              # haploid carrier
    assert big.normalize_gt_for_alt("0", 1) is None             # haploid non-carrier


_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr22>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3
chr22\t100\trs1\tAT\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0
chr22\t200\trs2\tG\tGATAA\t.\tPASS\t.\tGT:DP\t0|1:9\t1|1:12\t0|0:8
chr22\t300\trsSNP\tC\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t1|1
chr22\t400\trs3\tG\tA,GT\t.\tPASS\t.\tGT\t1|2\t0|0\t2|2
"""


def _write_vcf(gz=False):
    suffix = ".vcf.gz" if gz else ".vcf"
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(_VCF)
    else:
        with open(path, "w") as fh:
            fh.write(_VCF)
    return path


def test_iter_indels_skips_snp_and_parses_carriers():
    vcf = _write_vcf(gz=True)
    try:
        recs = list(big.iter_indel_genotypes(vcf))
    finally:
        os.remove(vcf)
    # pos 300 is a SNP (C>T) -> skipped; pos 400 has a SNP alt (A) skipped + an indel alt (GT)
    by_key = {(p, r, a): car for p, r, a, car in recs}
    assert (100, "AT", "A") in by_key
    assert by_key[(100, "AT", "A")] == {"HG1": "1|0", "HG2": "0|1"}  # HG3 0|0 dropped
    assert (200, "G", "GATAA") in by_key
    assert by_key[(200, "G", "GATAA")] == {"HG1": "0|1", "HG2": "1|1"}
    # the SNP at 300 must NOT appear
    assert not any(p == 300 for p, r, a, car in recs)
    # pos 400: ALT="A,GT"; 'A' is a SNP (skip), 'GT' (alt #2) is an insertion
    assert (400, "G", "GT") in by_key
    assert by_key[(400, "G", "GT")] == {"HG1": "0|1", "HG3": "1|1"}  # HG1 1|2 -> slot1, HG3 2|2 -> 1|1


def test_compile_and_reader_roundtrip():
    vcf = _write_vcf(gz=True)
    fd, store = tempfile.mkstemp(suffix=".tsv.gz")
    os.close(fd)
    try:
        n = big.compile_indel_genotypes(vcf, store)
        assert n == 3  # three indel alts (100, 200, 400)
        r = big.IndelGenotypeReader(store)
        assert r.carriers_dict(100, "AT", "A") == {"HG1": "1|0", "HG2": "0|1"}
        assert r.carriers_dict(400, "G", "GT") == {"HG1": "0|1", "HG3": "1|1"}
        assert r.carriers(300, "C", "T") == ""   # SNP absent from store
    finally:
        os.remove(vcf)
        os.remove(store)


def test_keep_samples_panel_filter():
    vcf = _write_vcf(gz=True)
    try:
        recs = {(p, r, a): car for p, r, a, car in
                big.iter_indel_genotypes(vcf, keep_samples={"HG1"})}
    finally:
        os.remove(vcf)
    assert recs[(100, "AT", "A")] == {"HG1": "1|0"}  # HG2 filtered out
    assert (200, "G", "GATAA") in recs and recs[(200, "G", "GATAA")] == {"HG1": "0|1"}


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

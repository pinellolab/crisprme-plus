#!/usr/bin/env python
"""Guard: --samplesID over-listing the VCF-genotyped panel is WARNed (#46).

Phantom samplesID entries (listed but genotyped in no VCF) are counted as hom-ref
when the Tier-0 panel AN is built, inflating AN and silently DEFLATING every
reported allele frequency (Samples/AF/MAF) -- shared by the classic SNP path and
the indel+SNP joint-AF. validate_inputs.check_samplesid_overlisting surfaces it at
pre-flight, and is correct for multi-dataset runs (union of all VCF headers).

Run: python PostProcess/test_samplesid_overlisting.py
"""

import gzip
import os
import tempfile

import validate_inputs as vi


def _mkvcf(d, name, samples):
    p = os.path.join(d, name)
    with gzip.open(p, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")
    return p


def test_overlisting_warns_and_quantifies():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    iss = vi.check_samplesid_overlisting([v], ["S1", "S2", "S3", "PH1", "PH2"])
    assert len(iss) == 1 and iss[0].severity == vi.WARN, iss
    assert "2 sample" in iss[0].message and "registry_fix_an" in iss[0].message


def test_matched_panel_is_silent():
    d = tempfile.mkdtemp()
    v = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    assert vi.check_samplesid_overlisting([v], ["S1", "S2", "S3"]) == []


def test_multidataset_union_no_false_positive():
    d = tempfile.mkdtemp()
    v1 = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2", "S3"])
    v2 = _mkvcf(d, "chr2.vcf.gz", ["S4", "S5"])
    # every sample is genotyped in SOME vcf -> no phantom, even though each is
    # absent from the other dataset's header
    assert vi.check_samplesid_overlisting([v1, v2], ["S1", "S2", "S3", "S4", "S5"]) == []


def test_multidataset_real_phantom_flagged():
    d = tempfile.mkdtemp()
    v1 = _mkvcf(d, "chr1.vcf.gz", ["S1", "S2"])
    v2 = _mkvcf(d, "chr2.vcf.gz", ["S3", "S4"])
    iss = vi.check_samplesid_overlisting([v1, v2], ["S1", "S2", "S3", "S4", "GHOST"])
    assert len(iss) == 1 and "GHOST" in iss[0].message


def test_empty_inputs_are_silent():
    assert vi.check_samplesid_overlisting([], ["S1"]) == []
    assert vi.check_samplesid_overlisting(["x.vcf.gz"], []) == []


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

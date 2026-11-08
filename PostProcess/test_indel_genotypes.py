#!/usr/bin/env python
"""Tests for build_indel_genotypes: phase normalization + VCF parse + reader.

Run: python PostProcess/test_indel_genotypes.py   (or via pytest)
"""

import unittest
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


def _mkstore():
    fd, store = tempfile.mkstemp(suffix=".tsv.gz")
    os.close(fd)
    os.remove(store)  # start absent; compile creates it
    return store


def _cleanup(*paths):
    for p in paths:
        for q in (p, p + ".tmp", big._done_marker(p)):
            try:
                os.remove(q)
            except OSError:
                pass


def test_compile_is_atomic_and_marks_done():
    # A clean compile leaves the final store + a .done marker (record count) and NO
    # leftover .tmp; store_is_complete must accept it via the fast marker path.
    vcf = _write_vcf(gz=True)
    store = _mkstore()
    try:
        n = big.compile_indel_genotypes(vcf, store)
        assert n == 3
        assert os.path.isfile(store)
        assert not os.path.exists(store + ".tmp")           # atomic: no temp residue
        marker = big._done_marker(store)
        assert os.path.isfile(marker)
        with open(marker) as m:                             # the count, then the FILTER rule
            assert m.read().splitlines() == ["3", big.FILTER_MARK]
        assert big.store_is_complete(store) is True
        assert big.store_is_current(store) is True
    finally:
        _cleanup(store)
        os.remove(vcf)


def test_store_is_complete_detects_truncation():
    # THE gt_indel_chr1 bug: a size>0 but truncated gzip must NOT be accepted as done,
    # so a resume rebuilds it instead of silently shipping a short store.
    vcf = _write_vcf(gz=True)
    store = _mkstore()
    try:
        big.compile_indel_genotypes(vcf, store)
        full = os.path.getsize(store)
        assert full > 20
        with open(store, "r+b") as fh:                      # lop off the gzip tail
            fh.truncate(full // 2)
        os.remove(big._done_marker(store))                  # marker gone -> integrity path
        assert os.path.getsize(store) > 0                   # size check WOULD have passed
        assert big.store_is_complete(store) is False        # integrity check catches it
    finally:
        _cleanup(store)
        os.remove(vcf)


def test_store_is_complete_adopts_valid_legacy_without_marker():
    # A store written by the old (pre-atomic) code has a valid gzip but no marker;
    # the integrity fallback must accept it so we don't needlessly rebuild good chroms.
    vcf = _write_vcf(gz=True)
    store = _mkstore()
    try:
        big.compile_indel_genotypes(vcf, store)
        os.remove(big._done_marker(store))
        assert big.store_is_complete(store) is True         # valid gzip -> adopted
    finally:
        _cleanup(store)
        os.remove(vcf)


def test_store_is_complete_rejects_absent_and_empty():
    store = _mkstore()
    try:
        assert big.store_is_complete(store) is False        # absent
        open(store, "wb").close()                           # zero-byte
        assert big.store_is_complete(store) is False
    finally:
        _cleanup(store)


_VCF_FILTER = """\
##fileformat=VCFv4.2
##contig=<ID=chr22>
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3
chr22\t100\trs1\tAT\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0
chr22\t150\trs4\tCAG\tC\t.\tLowQual\t.\tGT\t1|0\t0|0\t0|1
chr22\t200\trs2\tG\tGATAA\t.\t.\t.\tGT\t0|1\t1|1\t0|0
chr22\t250\trs5\tT\tTA,TAA\t.\tLowQual;ExcHet\t.\tGT\t1|2\t0|0\t2|0
chr22\t260\trs6\tGC\tG\t.\tpass\t.\tGT\t1|1\t0|0\t0|0
"""


def test_a_record_whose_filter_is_not_pass_or_dot_is_not_stored():
    # CRISPRitz's enricher skips a line whose FILTER is neither 'PASS' nor '.' (exact
    # strings): no fake-indel window holds its indel, so the store must not hold it either,
    # or a carrier is taken off a row for an indel no window holds (2.10.0)
    fd, vcf = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    with gzip.open(vcf, "wt") as fh:
        fh.write(_VCF_FILTER)
    store = _mkstore()
    try:
        keys = {(p, r, a) for p, r, a, _car in big.iter_indel_genotypes(vcf)}
        assert keys == {(100, "AT", "A"), (200, "G", "GATAA")}, keys
        assert big.compile_indel_genotypes(vcf, store) == 2
        r = big.IndelGenotypeReader(store)
        assert r.indels() == [(100, "AT", "A"), (200, "G", "GATAA")]
        assert r.carriers(150, "CAG", "C") == ""            # LowQual
        assert r.carriers(250, "T", "TAA") == ""            # LowQual;ExcHet, a multi-ALT line
        assert r.carriers(260, "GC", "G") == ""             # 'pass' is not 'PASS'
        assert r.reads_filter is True
    finally:
        _cleanup(store)
        os.remove(vcf)


def test_a_store_written_before_the_filter_rule_is_complete_but_not_current():
    # a store of 2.8.x: its .done marker holds the record count alone, and it may hold the
    # indels of records of any FILTER; a build writes it again (store_is_current), and a
    # reader says which it read (reads_filter)
    vcf = _write_vcf(gz=True)
    store = _mkstore()
    try:
        big.compile_indel_genotypes(vcf, store)
        with open(big._done_marker(store), "w") as m:
            m.write("3\n")
        assert big.store_is_complete(store) is True
        assert big.store_reads_filter(store) is False
        assert big.store_is_current(store) is False
        assert big.IndelGenotypeReader(store).reads_filter is False
        os.remove(big._done_marker(store))                  # adopted without a marker
        assert big.store_is_complete(store) is True
        assert big.store_is_current(store) is False
        big.compile_indel_genotypes(vcf, store)             # written again: current
        assert big.store_is_current(store) is True
        # a marker older than its store marks no rule (as it marks no completeness)
        os.utime(big._done_marker(store), (1, 1))
        assert big.store_reads_filter(store) is False
    finally:
        _cleanup(store)
        os.remove(vcf)


def test_compile_leaves_no_partial_final_on_kill(monkeypatch=None):
    # Simulate an OOM/SIGKILL mid-write: the row generator dies after one record.
    # compile must NOT leave a partial FINAL file (only a temp, which it cleans up),
    # so the final path stays absent and a resume re-attempts from scratch.
    vcf = _write_vcf(gz=True)
    store = _mkstore()

    def _boom(_vcf, _keep=None):
        yield 100, "AT", "A", {"HG1": "1|0"}
        raise RuntimeError("simulated OOM kill mid-write")

    orig = big.iter_indel_genotypes
    big.iter_indel_genotypes = _boom
    try:
        raised = False
        try:
            big.compile_indel_genotypes(vcf, store)
        except RuntimeError:
            raised = True
        assert raised
        assert not os.path.exists(store)                    # no truncated final file
        assert not os.path.exists(store + ".tmp")           # temp cleaned up
        assert big.store_is_complete(store) is False        # resume will rebuild
    finally:
        big.iter_indel_genotypes = orig
        _cleanup(store)
        os.remove(vcf)



def load_tests(loader, tests, pattern):
    """Under ``python -m unittest`` too, each test_ function above, which it found none of,
    and said OK (suite P4)."""
    return unittest.TestSuite(unittest.FunctionTestCase(fn) for name, fn in
                              sorted(globals().items()) if name.startswith("test_")
                              and callable(fn))

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

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


def test_an_indel_written_on_two_lines_is_one_allele():
    """R3B-N1: two VCF lines of one indel REF>ALT (HPRC chr22 writes 39 indels so) are one
    allele (decision B2): the store keeps a line per record, and its reader merges them,
    every carrier of either kept, a carrier of both with its two calls merged as CRISPRitz
    2.10.0's enricher merges two listings of one SNV (site_alleles.merge_gt). The reader
    kept the later line alone: HG2 and HG3, on the earlier only, lost their phased
    genotype, and HG1 its slot 0 (the review's world WD: S02 and S03 on no row)."""
    vcf_text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr22>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3\tHG4\n"
        "chr22\t100\t.\tAT\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t1/0\t0/1\n"
        "chr22\t100\t.\tAT\tA\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t1|0\n")
    fd, vcf = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    with gzip.open(vcf, "wt") as fh:
        fh.write(vcf_text)
    fd, store = tempfile.mkstemp(suffix=".tsv.gz")
    os.close(fd)
    try:
        assert big.compile_indel_genotypes(vcf, store) == 2   # one line per record
        r = big.IndelGenotypeReader(store)
        # HG4: 0/1 then 1|0, not phased alike, one copy each: the earlier, as it is
        assert r.carriers_dict(100, "AT", "A") == {"HG1": "1|1", "HG2": "0|1", "HG3": "1/0",
                                                   "HG4": "0/1"}
        assert r.indels() == [(100, "AT", "A")]
    finally:
        os.remove(vcf)
        _cleanup(store)
    assert big.merge_store_carriers("A:1|0,B:0|1", "A:0|1,C:1") == "A:1|1,B:0|1,C:1"
    assert big.merge_store_carriers("A:1|0", "A:0/1") == "A:1|0"
    assert big.merge_store_carriers("", "A:0|1") == "A:0|1"
    assert big.merge_store_carriers("A:0|1", "") == "A:0|1"


def test_an_indel_written_in_two_cases_is_one_allele():
    """C-R1-1: VCF 4.x reads REF and ALT whatever their case, and so does the fake-indel
    window; HPRC's vcfwave writes 140 indels of chr22 on two lines in two cases
    (chr22:23541792 TG>T beside tg>t). The store keeps each line as the VCF writes it,
    and its reader keyed it so: the two lines were two indels, and the lower-case one was
    found by no upper-case lookup. Now one allele, keyed in upper case, its carriers
    merged as for two lines written alike: HG1 holds the deletion on both slots, HG2 (the
    lower-case line's alone) keeps its call."""
    vcf_text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr22>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3\n"
        "chr22\t100\t.\tTG\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|1\n"
        "chr22\t100\t.\ttg\tt\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\n"
        "chr22\t300\t.\tAc\ta\t.\tPASS\t.\tGT\t0|0\t0|1\t0|0\n")
    fd, vcf = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    with gzip.open(vcf, "wt") as fh:
        fh.write(vcf_text)
    fd, store = tempfile.mkstemp(suffix=".tsv.gz")
    os.close(fd)
    try:
        assert big.compile_indel_genotypes(vcf, store) == 3   # one line per record
        r = big.IndelGenotypeReader(store)
        want = {"HG1": "1|1", "HG3": "0|1", "HG2": "1|0"}
        assert r.carriers_dict(100, "TG", "T") == want
        assert r.carriers_dict(100, "tg", "t") == want
        assert r.carriers_dict(300, "AC", "A") == {"HG2": "0|1"}
        assert r.indels() == [(100, "TG", "T"), (300, "AC", "A")]
    finally:
        os.remove(vcf)
        _cleanup(store)
    assert big.store_key(100, "tg", "t") == "100_TG_T"


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



# --- the store's MNV file (2.10.1, the plan of the MNVs' phase on HPRC) --------------------
# Up to 2.10.0 the store kept the records whose REF and ALT lengths differ, so nothing a
# search read said which haplotype holds an MNV. From 2.10.1 each chromosome's store has an
# MNV file beside its indel file, in the same format and with the same marker. Each test
# below fails at c64e62e1 (2.10.0 and 2.10.1's quick wins).

_VCF_MNV = """\
##fileformat=VCFv4.2
##contig=<ID=chr22>
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3
chr22\t100\trs1\tAT\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0
chr22\t120\trsM1\tAA\tGG\t.\tPASS\t.\tGT\t1|.\t1/.\t1
chr22\t130\trsM2\tTCA\tGAT\t.\t.\t.\tGT\t0|1\t0|0\t.|1
chr22\t140\trsM3\tCG\tTA\t.\tLowQual\t.\tGT\t1|1\t0|0\t0|0
chr22\t150\trsMX\tAC\tGT,A,AG\t.\tPASS\t.\tGT\t1|2\t3|0\t0|0
chr22\t160\trsS\tC\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t1|1
chr22\t170\trsSV\tA\t<DEL>\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0
chr22\t180\trsSY\tACGTA\t<INV>\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0
"""


def _vcf_file(text):
    fd, vcf = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    with gzip.open(vcf, "wt") as fh:
        fh.write(text)
    return vcf


def _lines(path):
    with gzip.open(path, "rt") as fh:
        return fh.read().splitlines()


def test_an_mnv_is_ref_and_alt_of_one_length_above_one():
    assert big._is_mnv("AA", "GG") and big._is_mnv("TCA", "TGA")
    assert not big._is_mnv("A", "G")            # an SNV: the SNP genome codes it
    assert not big._is_mnv("AT", "A") and not big._is_mnv("A", "AT")   # indels
    assert not big._is_mnv("ACGTA", "<INV>")    # a symbolic allele of five characters
    assert not big._is_mnv("AC", "*C")
    assert not big._is_indel("AA", "GG")        # the indel file keeps none


def test_the_mnv_file_its_marker_and_filter():
    """compile_store writes the indel file and the MNV file of one chromosome in one pass:
    each MNV ALT of a record whose FILTER is PASS or '.', read as "this ALT or not" with its
    separator (1|. is 1|0, 1/. is 1/0, a haploid 1 is 1), carriers only; a line of several
    ALTs gives each MNV ALT its own line (its genotype read against its own index; AC>AG,
    one base changed but REF and ALT of one length, is an MNV by the rule, as CRISPRitz's
    enricher windows it) and its indel ALT to the indel file; SNVs and symbolic alleles go
    to neither. Each file atomic,
    with the marker of 2.10.0's (its count, then FILTER_MARK), current; the reader reads
    both (has_mnv), indels() giving both kinds."""
    vcf = _vcf_file(_VCF_MNV)
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    try:
        assert big.compile_store(vcf, indel_out, mnv_out) == (2, 4)
        assert _lines(mnv_out) == ["120_AA_GG\tHG1:1|0,HG2:1/0,HG3:1",
                                   "130_TCA_GAT\tHG1:0|1,HG3:0|1",
                                   "150_AC_GT\tHG1:1|0", "150_AC_AG\tHG2:1|0"]
        assert _lines(indel_out) == ["100_AT_A\tHG1:1|0,HG2:0|1", "150_AC_A\tHG1:0|1"]
        for out, n in ((indel_out, 2), (mnv_out, 4)):
            with open(big._done_marker(out)) as m:
                assert m.read().splitlines() == [str(n), big.FILTER_MARK]
            assert big.store_is_current(out)
            assert not os.path.exists(out + ".tmp")
        r = big.IndelGenotypeReader(indel_out)
        assert r.has_mnv is True and r.reads_filter is True
        assert r.carriers_dict(120, "AA", "GG") == {"HG1": "1|0", "HG2": "1/0", "HG3": "1"}
        assert r.carriers_dict(130, "tca", "gat") == {"HG1": "0|1", "HG3": "0|1"}
        assert r.carriers(140, "CG", "TA") == ""    # LowQual
        assert r.indels() == [(100, "AT", "A"), (120, "AA", "GG"), (130, "TCA", "GAT"),
                              (150, "AC", "A"), (150, "AC", "AG"), (150, "AC", "GT")]
        # the indel file alone, as a reader of 2.10.0 opens it: no MNV
        r0 = big.IndelGenotypeReader(indel_out, mnv_path="")
        assert r0.has_mnv is False and r0.carriers(120, "AA", "GG") == ""
        assert big.mnv_store_path(indel_out) == mnv_out
        assert big.mnv_store_path(os.path.join(d, "store.tsv.gz")) is None
        assert ([(k, p, r, a) for k, p, r, a, _c in big.iter_store_genotypes(vcf)]
                == [("indel", 100, "AT", "A"), ("mnv", 120, "AA", "GG"),
                    ("mnv", 130, "TCA", "GAT"), ("mnv", 150, "AC", "GT"),
                    ("indel", 150, "AC", "A"), ("mnv", 150, "AC", "AG")])
        assert [p for p, _r, _a, _c in big.iter_mnv_genotypes(vcf)] == [120, 130, 150, 150]
        assert [p for p, _r, _a, _c in big.iter_indel_genotypes(vcf)] == [100, 150]
    finally:
        os.remove(vcf)
        _cleanup(indel_out, mnv_out)
        os.rmdir(d)


def test_one_mnv_on_two_lines_or_in_two_cases_is_one_allele():
    """An MNV the VCF writes on two lines, in one case or in two, is one allele, as an
    indel is (decision B2): the MNV file keeps a line per record, and the reader merges
    them (merge_store_carriers)."""
    vcf = _vcf_file(
        "##fileformat=VCFv4.2\n##contig=<ID=chr22>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3\n"
        "chr22\t120\t.\tAA\tGG\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\n"
        "chr22\t120\t.\tAA\tGG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n"
        "chr22\t120\t.\taa\tgg\t.\tPASS\t.\tGT\t0|0\t0|0\t1|0\n")
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    try:
        assert big.compile_store(vcf, indel_out, mnv_out) == (0, 3)
        r = big.IndelGenotypeReader(indel_out)
        assert r.carriers_dict(120, "AA", "GG") == {"HG1": "1|1", "HG2": "0|1", "HG3": "1|0"}
        assert r.indels() == [(120, "AA", "GG")]
    finally:
        os.remove(vcf)
        _cleanup(indel_out, mnv_out)
        os.rmdir(d)


def test_a_store_of_2100_gains_its_mnv_file_and_keeps_its_indel_file():
    """A store of 2.10.0 (an indel file, current) gains its MNV file: compile_store asked
    for the files that are not current (files_to_write: the MNV file alone, as
    search_indel_snp.store asks) leaves the indel file and its marker byte for byte, mtimes
    included, and then asks for none. Read before, the store had no MNV (has_mnv False)."""
    vcf = _vcf_file(_VCF_MNV)
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    try:
        assert big.compile_indel_genotypes(vcf, indel_out) == 2
        os.utime(indel_out, (1000000, 1000000))
        os.utime(big._done_marker(indel_out), (1000001, 1000001))
        before = {p: (open(p, "rb").read(), os.path.getmtime(p))
                  for p in (indel_out, big._done_marker(indel_out))}
        assert big.IndelGenotypeReader(indel_out).has_mnv is False
        todo = big.files_to_write(d, "chr22")
        assert todo == (None, mnv_out)
        assert big.files_to_write(d, "chr22", "/scratch") == (None, "/scratch/gt_mnv_chr22.tsv.gz")
        assert big.compile_store(vcf, *todo) == (None, 4)
        assert big.files_to_write(d, "chr22") == (None, None)
        assert {p: (open(p, "rb").read(), os.path.getmtime(p)) for p in before} == before
        assert big.store_is_current(mnv_out)
        r = big.IndelGenotypeReader(indel_out)
        assert r.has_mnv and r.carriers_dict(130, "TCA", "GAT") == {"HG1": "0|1",
                                                                     "HG3": "0|1"}
        assert big.compile_store(vcf) == (None, None)  # nothing asked, nothing read
    finally:
        os.remove(vcf)
        _cleanup(indel_out, mnv_out)
        os.rmdir(d)


def test_a_vcf_with_no_mnv_gets_an_mnv_file_of_no_record():
    """Its MNV file says the MNVs were stored (none): current, read (has_mnv), and the
    store's records its indels alone, as from its indel file."""
    vcf = _write_vcf(gz=True)
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    try:
        assert big.compile_store(vcf, indel_out, mnv_out) == (3, 0)
        assert _lines(mnv_out) == []
        assert big.store_is_current(mnv_out)
        r = big.IndelGenotypeReader(indel_out)
        assert r.has_mnv is True
        assert r.indels() == big.IndelGenotypeReader(indel_out, mnv_path="").indels()
        # the indel file is the one compile_indel_genotypes writes, byte for byte
        alone = os.path.join(d, "alone.tsv.gz")
        big.compile_indel_genotypes(vcf, alone)
        assert _lines(alone) == _lines(indel_out)
        _cleanup(alone)
    finally:
        os.remove(vcf)
        _cleanup(indel_out, mnv_out)
        os.rmdir(d)


def test_compile_store_leaves_no_part_of_a_file_on_a_kill():
    """A pass that dies (an OOM kill, an interrupt) leaves neither file nor a temporary
    one, and no marker; compile_store_safe returns None on any other failure (a tuple for a
    chromosome written, so build-index-only's STEP 1c counts the written ones by `is not
    None`; it returned -1, which STEP 1c's `>= 0` could not compare with a tuple), the
    temporary files and markers gone, and raises a failure that ends the build (a full
    disk). compile_indel_genotypes_safe keeps its contract: a count, or -1."""
    import errno
    vcf = _vcf_file(_VCF_MNV)
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    orig = big.iter_store_genotypes

    def boom(_vcf, _keep=None, _kinds=("indel", "mnv")):
        yield "mnv", 120, "AA", "GG", {"HG1": "1|0"}
        raise RuntimeError("simulated OOM kill mid-write")

    def full(_vcf, _keep=None, _kinds=("indel", "mnv")):
        yield "mnv", 120, "AA", "GG", {"HG1": "1|0"}
        raise OSError(errno.ENOSPC, "No space left on device")

    try:
        for out in (indel_out, mnv_out):
            with open(big._done_marker(out), "w") as m:
                m.write("1\n%s\n" % big.FILTER_MARK)    # a stale marker
        big.iter_store_genotypes = boom
        try:
            big.compile_store(vcf, indel_out, mnv_out)
            raise AssertionError("no failure")
        except RuntimeError:
            pass
        assert sorted(os.listdir(d)) == []
        assert big.compile_store_safe(vcf, indel_out, mnv_out) is None
        assert sorted(os.listdir(d)) == []
        assert big.compile_indel_genotypes_safe(vcf, indel_out) == -1
        assert sorted(os.listdir(d)) == []
        big.iter_store_genotypes = full
        try:
            big.compile_store_safe(vcf, indel_out, mnv_out)
            raise AssertionError("a full disk did not end the build")
        except OSError as err:
            assert big.ends_the_build(err)
        assert sorted(os.listdir(d)) == []
    finally:
        big.iter_store_genotypes = orig
        os.remove(vcf)
        _cleanup(indel_out, mnv_out)
        os.rmdir(d)


def test_the_mnv_only_command_line():
    """build_indel_genotypes.py --mnv-only <vcf> <store folder> [<samplesID>] (HPRC's
    rebuild): the MNV file of the VCF's chromosome (named as build-index-only names it, a
    chr prefix added) in the folder, with the build's samples; every other file of the
    folder as it was. A folder or a samplesID file named but not found is an error."""
    import contextlib
    import io
    vcf = _vcf_file(_VCF_MNV.replace("chr22\t", "22\t"))
    d = tempfile.mkdtemp()
    indel_out, mnv_out = big.store_files(d, "chr22")
    panel = os.path.join(d, "samples.txt")
    with open(panel, "w") as fh:
        fh.write("SAMPLE_ID\tPOP\nHG1\tX\nHG2\tX\n")
    try:
        with open(indel_out, "wb") as fh:
            fh.write(b"a 2.10.0 store's indel file")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            assert big.main(["build_indel_genotypes.py", "--mnv-only", vcf, d, panel]) == 0
        assert "4 phased MNV records -> %s" % mnv_out in out.getvalue()
        assert _lines(mnv_out) == ["120_AA_GG\tHG1:1|0,HG2:1/0", "130_TCA_GAT\tHG1:0|1",
                                   "150_AC_GT\tHG1:1|0", "150_AC_AG\tHG2:1|0"]
        assert big.store_is_current(mnv_out)
        with open(indel_out, "rb") as fh:
            assert fh.read() == b"a 2.10.0 store's indel file"
        assert sorted(os.listdir(d)) == sorted([os.path.basename(indel_out), "samples.txt",
                                                os.path.basename(mnv_out),
                                                os.path.basename(mnv_out) + ".done"])
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            assert big.main(["x", "--mnv-only", vcf, d, panel + ".missing"]) == 2
            assert big.main(["x", "--mnv-only", vcf, os.path.join(d, "nofolder")]) == 2
            assert big.main(["x", "--mnv-only", vcf]) == 2
        assert "no samplesID file" in err.getvalue() and "no store folder" in err.getvalue()
        assert big.vcf_chrom(vcf) == "chr22"
    finally:
        os.remove(vcf)
        _cleanup(indel_out, mnv_out, panel)
        os.rmdir(d)


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

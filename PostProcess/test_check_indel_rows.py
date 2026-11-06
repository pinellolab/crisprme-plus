#!/usr/bin/env python3
"""Unit tests for ``scripts/check_indel_rows.py``, the independent real-data check of
CRISPRme's indel (and MNV) off-target rows.

The checker must CATCH a wrong row and PASS a right one, so the suite is built both
ways: hand-built true rows that must raise nothing, and mutation tests that corrupt each
of those rows one field at a time -- the start shifted by 1 and by an insertion's
offset, one target base changed, a carrier swapped for a non-carrier, a named variant
dropped, a cross-haplotype SNP added, and a phantom of the shape Manuel Tognon reported
(an insertion's only carrier listed under SNP alleles only other individuals carry) --
and assert that each one is reported, with the class it belongs to.

Everything is built here: a small reference, a phased VCF written and bgzipped with
Python's own gzip (BGZF is a gzip stream, which the checker's reader accepts), and rows
in the column layout of ``integrated_results.tsv``. No CRISPRme module is imported, as
the checker imports none.
"""

import gzip
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                                "scripts"))
import check_indel_rows as C  # noqa: E402

GUIDE = "ACTGAAATCTGTAAGCAGGC"
SAMPLES = ["S%02d" % i for i in range(1, 9)]


def write_fasta(path, seq, name="chr21"):
    with open(path, "w") as fh:
        fh.write(">%s\n" % name)
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def write_vcf(path, chrom, records, samples=SAMPLES):
    """``records``: (pos1, ref, alt, {sample: gt}). Written bgzip-less: a plain gzip
    stream is a legal single BGZF-free member and the checker's reader streams it."""
    body = ["##fileformat=VCFv4.2", "##contig=<ID=%s>" % chrom,
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for pos1, ref, alt, gts in records:
        row = [chrom, str(pos1), ".", ref, alt, ".", "PASS", "AC=1", "GT"]
        row += [gts.get(s, "0|0") for s in samples]
        body.append("\t".join(row))
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(body) + "\n")


HEAD = (["Spacer+PAM", "Chromosome"]
        + sum(([
            "Start_coordinate_(%s)" % b, "Strand_(%s)" % b, "Aligned_spacer+PAM_(%s)" % b,
            "Aligned_protospacer+PAM_REF_(%s)" % b, "Aligned_protospacer+PAM_ALT_(%s)" % b,
            "PAM_(%s)" % b, "Mismatches_(%s)" % b, "Bulges_(%s)" % b,
            "Mismatches+bulges_(%s)" % b, "Bulge_type_(%s)" % b,
            "REF/ALT_origin_(%s)" % b, "Variant_info_genome_(%s)" % b,
            "Variant_samples_(%s)" % b] for b in C.BLOCKS), []))


def row_fields(chrom, start, strand, crrna, ref, alt, pam, mm, bul, btype, var, samples):
    out = [GUIDE + "NNN", chrom]
    for _b in C.BLOCKS:
        out += [str(start), strand, crrna, ref, alt, pam, str(mm), str(bul),
                str(mm + bul), btype, "alt", var, samples]
    return out


def write_results(folder, rows):
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "G+NRG_hg38+v_4+2_integrated_results.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(HEAD) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path


class Case(object):
    """One built world: a reference, a VCF, a results folder and the argv to check it."""

    def __init__(self, tmp, seq, records, rows, chrom="chr21"):
        self.tmp = tmp
        self.chrom = chrom
        os.makedirs(os.path.join(tmp, "Genomes"), exist_ok=True)
        write_fasta(os.path.join(tmp, "Genomes", chrom + ".fa"), seq, chrom)
        self.vcf = os.path.join(tmp, "v.vcf.gz")
        write_vcf(self.vcf, chrom, records)
        self.pam = os.path.join(tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")
        self.res = os.path.join(tmp, "Results")
        write_results(self.res, rows)

    def run(self, mode="soundness", extra=()):
        out = os.path.join(self.tmp, "fail.tsv")
        argv = ["--results", self.res, "--reference", os.path.join(self.tmp, "Genomes"),
                "--vcf", self.vcf, "--pam", self.pam, "--mm", "4", "--bDNA", "1",
                "--bRNA", "1", "--max-total-edits", "6", "--mode", mode,
                "--chrom", self.chrom, "--guide", GUIDE, "--out", out] + list(extra)
        rep_classes = []
        rc = C.main(argv)
        if os.path.exists(out):
            with open(out) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = dict(zip(head, line.rstrip("\n").split("\t")))
                    rep_classes.append((f["gating"], f["class"], f["detail"]))
        return rc, rep_classes


# --------------------------------------------------------------------------- #
# the worlds
# --------------------------------------------------------------------------- #
def world_deletion(tmp):
    """A deletion restores a perfect target: ...G[c]UIDE... with one extra base deleted.
    S01 and S02 carry it, S03 does not."""
    pre = "GGCATTACGGTTACCATGCAAT"
    seq = pre + GUIDE[:10] + "C" + GUIDE[10:] + "TGG" + "ACGTTGCAATTGGCCAATTGCA"
    pos1 = len(pre) + 10            # 1-based POS of the anchor base GUIDE[9]
    rec = (pos1, GUIDE[9] + "C", GUIDE[9], {"S01": "1|0", "S02": "0|1"})
    start = len(pre)                # the target's first base, 0-based
    var = "chr21_%d_%s_%s" % (pos1, GUIDE[9] + "C", GUIDE[9])
    ref_text = seq[start:start + 23]
    rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S01,S02")]
    return Case(tmp, seq, [rec], rows), dict(
        seq=seq, pos1=pos1, start=start, var=var, ref_text=ref_text, rows=rows)


def world_insertion(tmp):
    """An insertion holds the target's first bases: the row must start at the anchor, the
    last reference base the inserted bases follow. S04 carries it."""
    pre = "GGCATTACGGTTACCATGCAAT"
    tail = GUIDE[6:] + "TGG"
    seq = pre + tail + "ACGTTGCAATTGGCCAATTGCA"
    anchor0 = len(pre) - 1          # the base the insertion follows, 0-based
    pos1 = anchor0 + 1
    added = "TTCAT" + GUIDE[:6]
    rec = (pos1, seq[anchor0], seq[anchor0] + added, {"S04": "1|0"})
    var = "chr21_%d_%s_%s" % (pos1, seq[anchor0], seq[anchor0] + added)
    # the target's first base is the 6th base the insertion adds; the row starts at the
    # anchor and its REF text is read from there
    ref_text = seq[anchor0:anchor0 + 23]
    rows = [row_fields("chr21", anchor0, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S04")]
    return Case(tmp, seq, [rec], rows), dict(
        seq=seq, pos1=pos1, anchor0=anchor0, var=var, ref_text=ref_text, rows=rows,
        added=added)


def world_indel_snp(tmp):
    """A deletion and a spacer SNP the guide needs, in cis on S05's haplotype and in
    trans on S06's. The row names both variants and lists S05 only."""
    pre = "GGCATTACGGTTACCATGCAAT"
    seq = pre + GUIDE[:5] + "C" + GUIDE[6:10] + "C" + GUIDE[10:] + "TGG" + "ACGTTGCAATT"
    dpos = len(pre) + 10
    spos = len(pre) + 6             # 1-based POS of the SNP (GUIDE[5] -> C in the ref)
    recs = [(spos, "C", GUIDE[5], {"S05": "1|0", "S06": "0|1"}),
            (dpos, GUIDE[9] + "C", GUIDE[9], {"S05": "1|0", "S06": "1|0"})]
    start = len(pre)
    var = "chr21_%d_%s_%s,chr21_%d_C_%s" % (dpos, GUIDE[9] + "C", GUIDE[9], spos,
                                            GUIDE[5])
    ref_text = seq[start:start + 23]
    rows = [row_fields("chr21", start, "+", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                       "TGG", 0, 0, "X", var, "S05")]
    return Case(tmp, seq, recs, rows), dict(
        seq=seq, dpos=dpos, spos=spos, start=start, var=var, ref_text=ref_text,
        rows=rows)


class TrueRows(unittest.TestCase):
    """A right row raises nothing."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="chkidl_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_deletion_row_passes(self):
        case, _ = world_deletion(self.tmp)
        rc, found = case.run()
        self.assertEqual((rc, found), (0, []), found)

    def test_insertion_row_at_the_anchor_passes(self):
        case, _ = world_insertion(self.tmp)
        rc, found = case.run()
        self.assertEqual((rc, found), (0, []), found)

    def test_indel_plus_cis_snp_row_passes(self):
        case, _ = world_indel_snp(self.tmp)
        rc, found = case.run()
        self.assertEqual((rc, found), (0, []), found)

    def test_minus_strand_row_passes(self):
        pre = "GGCATTACGGTTACCATGCAAT"
        fwd = C.rc(GUIDE + "TGG")
        seq = pre + fwd[:9] + "A" + fwd[9:] + "ACGTTGCAATTGGCC"
        pos1 = len(pre) + 9         # anchor of the 1-bp deletion inside the target
        rec = (pos1, fwd[8] + "A", fwd[8], {"S07": "1|0"})
        start = len(pre)
        var = "chr21_%d_%s_%s" % (pos1, fwd[8] + "A", fwd[8])
        ref_text = C.rc(seq[start:start + 23])
        rows = [row_fields("chr21", start, "-", GUIDE + "NNN", ref_text, GUIDE + "TGG",
                           "TGG", 0, 0, "X", var, "S07")]
        case = Case(self.tmp, seq, [rec], rows)
        rc, found = case.run()
        self.assertEqual((rc, found), (0, []), found)


class Mutations(unittest.TestCase):
    """Each corruption of a right row is caught, with its class."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="chkidl_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _mutate(self, build, **kw):
        """Rebuild a world with one field of its single row replaced."""
        case, info = build(self.tmp)
        row = list(info["rows"][0])
        per = (len(row) - 2) // len(C.BLOCKS)
        for key, value in kw.items():
            off = {"start": 0, "strand": 1, "crrna": 2, "ref": 3, "alt": 4, "pam": 5,
                   "mm": 6, "bul": 7, "total": 8, "btype": 9, "origin": 10, "var": 11,
                   "samples": 12}[key]
            for b in range(len(C.BLOCKS)):
                row[2 + b * per + off] = str(value)
        write_results(case.res, [row])
        return case.run()

    # --- the start -------------------------------------------------------- #
    def test_start_shifted_by_one_is_caught(self):
        case, info = world_deletion(self.tmp)
        rc, found = self._mutate(world_deletion, start=info["start"] + 1)
        self.assertEqual(rc, 1)
        self.assertIn("start-not-on-carrier-haplotype", [f[1] for f in found], found)

    def test_start_shifted_by_the_insertion_offset_is_caught(self):
        """The 2.8.1 defect itself: the row placed at the first base's place on the
        haplotype instead of the insertion's anchor."""
        case, info = world_insertion(self.tmp)
        shifted = info["anchor0"] + 6   # the 6 bases the target's first base follows
        rc, found = self._mutate(world_insertion, start=shifted)
        self.assertEqual(rc, 1)
        self.assertIn("start-not-on-carrier-haplotype", [f[1] for f in found], found)

    def test_start_shifted_to_the_base_after_the_insertion_is_caught(self):
        case, info = world_insertion(self.tmp)
        rc, found = self._mutate(world_insertion, start=info["anchor0"] + 1)
        self.assertEqual(rc, 1)
        self.assertIn("start-not-on-carrier-haplotype", [f[1] for f in found], found)

    # --- the sequence ----------------------------------------------------- #
    def test_one_target_base_changed_is_caught(self):
        case, info = world_deletion(self.tmp)
        bad = "A" if GUIDE[4] != "A" else "C"
        alt = GUIDE[:4] + bad + GUIDE[5:] + "TGG"
        rc, found = self._mutate(world_deletion, alt=alt, mm=1)
        self.assertEqual(rc, 1)
        classes = [f[1] for f in found]
        self.assertTrue(
            {"target-not-on-carrier-haplotype", "start-not-on-carrier-haplotype"}
            & set(classes), found)

    def test_counts_not_recounting_is_caught(self):
        rc, found = self._mutate(world_deletion, mm=2)
        self.assertEqual(rc, 1)
        self.assertIn("counts-mismatch", [f[1] for f in found], found)

    def test_a_pam_the_search_never_matched_is_caught(self):
        case, info = world_deletion(self.tmp)
        alt = GUIDE + "TGA"
        seq = info["seq"]
        rc, found = self._mutate(world_deletion, alt=alt, pam="TGA")
        self.assertEqual(rc, 1)
        self.assertIn("pam-mismatch", [f[1] for f in found], found)

    def test_ref_text_read_at_the_wrong_place_is_caught(self):
        case, info = world_deletion(self.tmp)
        wrong = info["seq"][info["start"] + 3:info["start"] + 26]
        rc, found = self._mutate(world_deletion, ref=wrong)
        self.assertEqual(rc, 1)
        self.assertIn("ref-text-mismatch", [f[1] for f in found], found)

    # --- the carriers ----------------------------------------------------- #
    def test_a_carrier_swapped_for_a_non_carrier_is_caught(self):
        """The row's only carrier does not have the target: nobody listed does, so the
        row is reported as unplaced rather than once per carrier."""
        rc, found = self._mutate(world_deletion, samples="S03")
        self.assertEqual(rc, 1)
        self.assertIn("target-not-on-carrier-haplotype", [f[1] for f in found], found)

    def test_a_carrier_added_beside_the_real_ones_is_caught(self):
        rc, found = self._mutate(world_deletion, samples="S01,S02,S03")
        self.assertEqual(rc, 1)
        self.assertIn("carrier-without-target", [f[1] for f in found], found)

    def test_a_missing_carrier_is_caught(self):
        rc, found = self._mutate(world_deletion, samples="S01")
        self.assertEqual(rc, 1)
        self.assertIn("carrier-missing", [f[1] for f in found], found)

    # --- the variants ----------------------------------------------------- #
    def test_a_named_variant_dropped_is_caught(self):
        case, info = world_indel_snp(self.tmp)
        only_indel = info["var"].split(",")[0]
        rc, found = self._mutate(world_indel_snp, var=only_indel)
        self.assertEqual(rc, 1)
        self.assertIn("variant-missing", [f[1] for f in found], found)

    def test_a_variant_that_is_not_a_record_is_caught(self):
        case, info = world_deletion(self.tmp)
        rc, found = self._mutate(world_deletion,
                                 var=info["var"] + ",chr21_9_A_AAAA")
        self.assertEqual(rc, 1)
        self.assertIn("variant-not-in-vcf", [f[1] for f in found], found)

    def test_a_cross_haplotype_snp_added_to_the_row_is_caught(self):
        """S05 carries the SNP in cis, S06 in trans. Naming the SNP on a row that lists
        S06 is the cross-person pairing of the 2.8.0 defect, in miniature."""
        case, info = world_indel_snp(self.tmp)
        rc, found = self._mutate(world_indel_snp, samples="S05,S06")
        self.assertEqual(rc, 1)
        classes = [f[1] for f in found]
        self.assertTrue({"carrier-without-target", "carrier-without-variants"}
                        & set(classes), found)

    def test_the_phantom_manuel_reported_is_caught(self):
        """An insertion's only carrier listed for a sequence that pairs the insertion
        with SNP alleles other individuals carry: the 2.8.0 row at chr22:22322768 in
        miniature. S04 carries the insertion; S05 and S06 carry the two SNPs, each on
        its own, and nobody carries any two of the three."""
        pre = "GGCATTACGGTTACCATGCAAT"
        tail = GUIDE[6:] + "TGG"
        seq = pre + tail + "ACGTTGCAATTGGCCAATTGCA"
        anchor0 = len(pre) - 1
        pos1 = anchor0 + 1
        # the insertion alone gives the guide's first 6 bases with two bases wrong
        wrong = list(GUIDE[:6])
        wrong[2] = "A" if wrong[2] != "A" else "C"
        added = "TTCAT" + "".join(wrong)
        recs = [(pos1, seq[anchor0], seq[anchor0] + "".join(added), {"S04": "1|0"})]
        # two SNPs downstream of the insertion that other people carry, which together
        # with it would spell the guide
        p1 = len(pre) + 2
        recs.append((p1, seq[p1 - 1], "A" if seq[p1 - 1] != "A" else "C",
                     {"S05": "1|0"}))
        p2 = len(pre) + 5
        recs.append((p2, seq[p2 - 1], "A" if seq[p2 - 1] != "A" else "C",
                     {"S06": "0|1"}))
        phantom = GUIDE + "TGG"
        var = "chr21_%d_%s_%s" % (pos1, seq[anchor0], seq[anchor0] + added)
        ref_text = seq[anchor0:anchor0 + 23]
        rows = [row_fields("chr21", anchor0, "+", GUIDE + "NNN", ref_text, phantom,
                           "TGG", 0, 0, "X", var, "S04")]
        case = Case(self.tmp, seq, recs, rows)
        rc, found = case.run()
        self.assertEqual(rc, 1)
        classes = [f[1] for f in found]
        self.assertTrue({"target-not-on-carrier-haplotype",
                         "start-not-on-carrier-haplotype"} & set(classes), found)


class Completeness(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="chkidl_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_a_reported_target_satisfies_completeness(self):
        case, _ = world_deletion(self.tmp)
        rc, found = case.run(mode="completeness")
        self.assertEqual([f for f in found if f[0] == "yes"], [], found)

    def test_a_dropped_row_is_caught(self):
        case, info = world_deletion(self.tmp)
        write_results(case.res, [])
        rc, found = case.run(mode="completeness")
        self.assertEqual(rc, 1)
        self.assertIn("missing-row", [f[1] for f in found], found)

    def test_a_dropped_carrier_is_caught(self):
        case, info = world_deletion(self.tmp)
        row = list(info["rows"][0])
        per = (len(row) - 2) // len(C.BLOCKS)
        for b in range(len(C.BLOCKS)):
            row[2 + b * per + 12] = "S01"
        write_results(case.res, [row])
        rc, found = case.run(mode="completeness")
        self.assertEqual(rc, 1)
        self.assertTrue(any("S02" in f[2] for f in found if f[1] == "missing-row"), found)

    def test_a_flank_target_is_not_demanded(self):
        """A target in a deletion's downstream flank holds no base of the indel: the
        reference and SNP searches report it, and the checker must not ask the indel
        path for it."""
        pre = "GGCATTACGGTTACCATGCAAT"
        seq = pre + "GT" + "CA" + GUIDE + "TGG" + "ACGTTGCAATTGGCC"
        pos1 = len(pre) + 1
        rec = (pos1, "GT", "G", {"S01": "1|0"})
        case = Case(self.tmp, seq, [rec], [])
        rc, found = case.run(mode="completeness")
        self.assertEqual([f for f in found if f[0] == "yes"], [], found)


class Pieces(unittest.TestCase):
    """The checker's own primitives."""

    def test_the_junction_rule_for_a_deletion(self):
        # GT>G at POS 100 (0-based 99): the deleted base is 100, the junction between 99
        # and 101. A target must hold a base on each side.
        self.assertTrue(C.holds_indel(95, 105, 100, "GT", "G"))
        self.assertTrue(C.holds_indel(99, 101, 100, "GT", "G"))
        self.assertFalse(C.holds_indel(101, 120, 100, "GT", "G"))   # wholly downstream
        self.assertFalse(C.holds_indel(80, 99, 100, "GT", "G"))     # wholly upstream

    def test_the_junction_rule_for_an_insertion(self):
        # G>GAAA at POS 100: the added bases sit between 99 and 100
        self.assertTrue(C.holds_indel(95, 105, 100, "G", "GAAA"))
        self.assertFalse(C.holds_indel(100, 120, 100, "G", "GAAA"))
        self.assertFalse(C.holds_indel(80, 98, 100, "G", "GAAA"))

    def test_the_junction_rule_read_on_the_target_bases(self):
        """holds_record is exact where a span cannot be: a target holding the anchor and
        the bases an insertion adds has the span (anchor, anchor)."""
        # G>GAAA at POS 100: coordinates 97, 98, 99, then three added bases, then 100
        self.assertTrue(C.holds_record([97, 98, 99, None, None, None], 100, "G", "GAAA"))
        self.assertTrue(C.holds_record([None, None, 100, 101], 100, "G", "GAAA"))
        self.assertFalse(C.holds_record([97, 98, 99], 100, "G", "GAAA"))
        self.assertFalse(C.holds_record([100, 101, 102], 100, "G", "GAAA"))
        # GT>G at POS 100: the junction is the jump from 99 to 101
        self.assertTrue(C.holds_record([98, 99, 101, 102], 100, "GT", "G"))
        self.assertFalse(C.holds_record([101, 102, 103], 100, "GT", "G"))
        self.assertFalse(C.holds_record([96, 97, 98, 99], 100, "GT", "G"))
        # an MNV covers a base its ALT replaces
        self.assertTrue(C.holds_record([99, 100, 101], 100, "AC", "TG"))
        self.assertFalse(C.holds_record([101, 102], 100, "AC", "TG"))

    def test_the_junction_rule_for_an_mnv(self):
        self.assertTrue(C.holds_indel(99, 100, 100, "AC", "TG"))
        self.assertFalse(C.holds_indel(101, 120, 100, "AC", "TG"))

    def test_changes_bases(self):
        self.assertTrue(C.changes_bases(100, "GT", "G", 99, 101))
        self.assertFalse(C.changes_bases(100, "GT", "G", 102, 120))
        self.assertTrue(C.changes_bases(100, "G", "GAAA", 99, 100))
        self.assertFalse(C.changes_bases(100, "G", "GAAA", 100, 120))

    def test_the_aligner_finds_a_perfect_target_and_its_bulges(self):
        full = GUIDE + "NRG"
        seq = "TTTT" + GUIDE + "TGG" + "TTTT"
        hits = list(C.scan_sequence(full, seq, 4, 1, 1, 6, 3, False))
        perfect = [h for h in hits if h[4:] == (0, 0, 0)]
        self.assertEqual(len(perfect), 1)
        self.assertEqual((perfect[0][0], perfect[0][1]), ("+", 4))
        self.assertEqual(perfect[0][3], GUIDE + "TGG")

    def test_the_aligner_puts_no_gap_in_the_pam(self):
        full = GUIDE + "NRG"
        seq = "TTTT" + GUIDE + "TGG" + "TTTT"
        for _st, _off, sa, ta, _mm, _dn, _rn in C.scan_sequence(
                full, seq, 4, 1, 1, 6, 3, False):
            self.assertNotIn("-", sa[-3:])
            self.assertNotIn("-", ta[-3:])

    def test_the_aligner_reports_a_dna_bulge_right_before_the_pam(self):
        """CRISPRme does (``CTAACAGTTGCTTTTATCAC-NNN``), so the checker must too."""
        full = GUIDE + "NRG"
        seq = "TTTT" + GUIDE + "A" + "TGG" + "TTTT"
        hits = [h for h in C.scan_sequence(full, seq, 4, 1, 1, 6, 3, False)
                if h[2].endswith("-NRG")]
        self.assertTrue(hits)
        self.assertEqual(hits[0][5], 1)  # one DNA bulge

    def test_recount_matches_a_crisprme_pair(self):
        mm, dn, rn, bad = C.recount("CTAACAGTTGCTTTTATCAC-NNN",
                                    "-TAACtGTTtCTgTTATCAgTAAG", 3, False)
        self.assertIsNone(bad)
        self.assertEqual((mm, dn, rn), (4, 1, 1))

    def test_first_coord_takes_an_inserted_base_to_its_anchor(self):
        coords = [10, 11, None, None, 12, 13]
        self.assertEqual(C.first_coord(coords, 0), 10)
        self.assertEqual(C.first_coord(coords, 2), 11)  # an inserted base -> the anchor
        self.assertEqual(C.first_coord(coords, 3), 11)
        self.assertEqual(C.first_coord(coords, 4), 12)

    def test_the_pam_file_is_read_both_ways(self):
        tmp = tempfile.mkdtemp(prefix="chkidl_")
        try:
            p3 = os.path.join(tmp, "nrg.txt")
            with open(p3, "w") as fh:
                fh.write("N" * 20 + "NRG 3\n")
            self.assertEqual(C.read_pam(p3), ("NRG", 3, False))
            p5 = os.path.join(tmp, "ttt.txt")
            with open(p5, "w") as fh:
                fh.write("TTTV" + "N" * 20 + " -4\n")
            self.assertEqual(C.read_pam(p5), ("TTTV", 4, True))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class HaplotypeBuild(unittest.TestCase):
    """The haplotype builder and its coordinate map, the checker's own truth."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="chkidl_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _window(self, seq, records):
        os.makedirs(os.path.join(self.tmp, "G"), exist_ok=True)
        write_fasta(os.path.join(self.tmp, "G", "chr21.fa"), seq)
        vcf = os.path.join(self.tmp, "v.vcf.gz")
        write_vcf(vcf, "chr21", records)
        return C.Window(C.Reference(os.path.join(self.tmp, "G")), "chr21", 0, len(seq),
                        [C.Vcf(vcf)], pad=0)

    def test_an_insertion_gives_its_added_bases_no_coordinate(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "A", "ATTT", {"S01": "1|0"})])
        hap, coords = win.haplotype("S01", 0, 0, 20)
        self.assertEqual(hap, seq[:5] + "TTT" + seq[5:20])
        self.assertEqual(coords[4], 4)
        self.assertEqual(coords[5:8], [None, None, None])
        self.assertEqual(coords[8], 5)

    def test_a_deletion_jumps_the_deleted_bases(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "ACG", "A", {"S01": "1|0"})])
        hap, coords = win.haplotype("S01", 0, 0, 20)
        self.assertEqual(hap, seq[:5] + seq[7:20])
        self.assertEqual(coords[4:6], [4, 7])

    def test_an_mnv_keeps_its_coordinates(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "AC", "TG", {"S01": "1|0"})])
        hap, coords = win.haplotype("S01", 0, 0, 20)
        self.assertEqual(hap[4:6], "TG")
        self.assertEqual(coords[4:6], [4, 5])

    def test_the_other_haplotype_is_the_reference(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "A", "ATTT", {"S01": "1|0"})])
        hap, coords = win.haplotype("S01", 1, 0, 20)
        self.assertEqual(hap, seq[:20])
        self.assertEqual(coords, list(range(20)))

    def test_an_unphased_sample_has_no_haplotype(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "A", "ATTT", {"S01": "0/1"})])
        self.assertEqual(win.haplotype("S01", 0, 0, 20), (None, None))
        self.assertFalse(win.phased("S01"))
        self.assertTrue(win.carries_allele("S01", 5, "A", "ATTT"))

    def test_a_vcf_line_with_two_alts_keeps_each_genotype_on_its_own_alt(self):
        seq = "ACGT" * 20
        win = self._window(seq, [(5, "A", "C,G", {"S01": "0|1", "S02": "2|1",
                                                  "S03": "2|2"})])
        self.assertEqual(set(win.carry[(5, "A", "C")]), {"S01", "S02"})
        self.assertEqual(set(win.carry[(5, "A", "G")]), {"S02", "S03"})
        self.assertEqual(win.slot_alleles("S02", 0), [(5, "A", "G")])
        self.assertEqual(win.slot_alleles("S02", 1), [(5, "A", "C")])


if __name__ == "__main__":
    unittest.main(verbosity=2)

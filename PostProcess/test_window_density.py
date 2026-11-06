#!/usr/bin/env python3
"""The ambiguity budget of a variant search (2.9.0): how many sequences a window of the
variant genome spells, and how many of its ambiguous columns a reported row rests on.

Coding every allele of a split multi-allelic site, and a four-allele position as Q, makes
windows of the immunoglobulin loci spell far more sequences than CRISPRME_IUPAC_CAP's count
of positions could see: a three-allele column was one position, as a two-allele one is.
  * CRISPRME_IUPAC_CAP counts bits (log2 of the sequences a window spells): a window of
    seven three-allele columns (2^11.09 sequences) is over the default 10, though it holds
    seven positions;
  * a window over the cap is dense, as since 2.8.1; one spelling 2^HYPERDENSE_BITS
    sequences or more is hyperdense, and is listed in the high-variant-density BED whatever
    the cap, and the BED gives each window's bits and four-allele positions;
  * the four-allele positions of the reported windows go to <out>.four_allele_positions.bed,
    and resultIntegrator flags a row whose window holds one in High_variant_density_region,
    which the report shows in High_complexity_region, with a hyperdense window's sequences;
  * complete-search's --max-degeneracy K: a worst-case alignment is cut back to rest on at
    most K ambiguous columns, and any other row resting on more is withheld, unless its
    carriers hold its variants together on one phased haplotype.

The real new_simple_analysis.py, resultIntegrator.py and crisprme.py run on hand-built inputs
(the tier path's observed haplotypes through test_phased_haplotype's harness), and the
report's functions on hand-built rows.

Run with:
    cd PostProcess && python3 -m unittest test_window_density -v
"""

import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

import generate_report as gr
import test_build_index_only_atomic as bio
import twopass_emit
import window_density as wd
from test_phased_haplotype import _load_pure_functions

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

HERE = os.path.dirname(os.path.abspath(__file__))
GUIDE = "CTAACAGTTGCTTTTATCAC"
CRRNA = GUIDE + "NNN"
FLANK = "ACGTGCATTC" * 10
PAM_NRG = "NNNNNNNNNNNNNNNNNNNNNRG 3\n"


def _degeneracy_of_row(row):
    """D of a finalized 3' NRG row of the post-analysis output: its DNA (column 2) against its
    Reference (column 18)."""
    return wd.degeneracy(row[2], row[18], row[1][:-3], "NRG", 0, -3, -3, None)


class TheSequencesAWindowSpells(unittest.TestCase):
    def test_each_letter_counts_the_bases_it_stands_for(self):
        self.assertEqual(wd.sequences("ACGT"), 1)
        self.assertEqual(wd.sequences("R"), 2)
        self.assertEqual(wd.sequences("y"), 2)
        self.assertEqual(wd.sequences("B"), 3)
        self.assertEqual(wd.sequences("Q"), 4)
        self.assertEqual(wd.sequences("RBQ"), 24)
        # an assembly gap matches nothing, and a bulge gap is no base
        self.assertEqual(wd.sequences("NNN--"), 1)
        self.assertEqual(wd.bits(wd.sequences("RYSW")), 4.0)

    def test_three_allele_columns_weigh_more_than_one_position(self):
        seven = wd.sequences("BDHVBDH")
        self.assertEqual(seven, 3 ** 7)
        self.assertAlmostEqual(wd.bits(seven), 11.09, places=2)
        self.assertTrue(wd.over_cap(seven, 10))
        # ten two-allele columns are not over a cap of 10 (2^10 is not more than 2^10),
        # eleven are: as the count of positions was
        self.assertFalse(wd.over_cap(wd.sequences("R" * 10), 10))
        self.assertTrue(wd.over_cap(wd.sequences("R" * 11), 10))

    def test_a_negative_cap_is_none(self):
        self.assertFalse(wd.over_cap(wd.sequences("Q" * 23), -1))

    def test_hyperdense_from_two_to_the_twenty(self):
        self.assertEqual(wd.HYPERDENSE_BITS, 20)
        self.assertTrue(wd.hyperdense(wd.sequences("R" * 20)))
        self.assertFalse(wd.hyperdense(wd.sequences("R" * 19)))
        self.assertTrue(wd.hyperdense(wd.sequences("Q" * 10)))

    def test_four_allele_offsets(self):
        self.assertEqual(wd.four_allele_offsets("ACQTqG-R"), [2, 4])
        self.assertEqual(wd.four_allele_offsets("ACGT"), [])


class TheColumnsARowRestsOn(unittest.TestCase):
    # 3' NRG: the protospacer is [0:-3], the PAM [-3:]
    GEOM = (0, -3, -3, None)

    def d(self, row, ref, guide=GUIDE, pam="NRG", geom=GEOM):
        return wd.degeneracy(row, ref, guide, pam, *geom)

    def test_a_mismatch_the_row_matches_counts(self):
        ref = "CTgACAGTTcCTTTTATCACTGG"   # two mismatches (lower case, as the finalizer has them)
        self.assertEqual(self.d(GUIDE + "TGG", ref), 2)
        self.assertEqual(self.d("CTAACAGTTcCTTTTATCACTGG", ref), 1)
        self.assertEqual(self.d(ref, ref), 0)

    def test_a_match_the_row_breaks_or_a_mismatch_it_changes_does_not(self):
        ref = "CTgACAGTTGCTTTTATCACTGG"
        # column 2: G -> C, a mismatch still; column 5: A -> t, breaking a match
        self.assertEqual(self.d("CTcACtGTTGCTTTTATCACTGG", ref), 0)

    def test_a_pam_the_row_mends_counts(self):
        # NRG: the reference's T at R does not fit, the row's A does
        self.assertEqual(self.d(GUIDE + "TAG", GUIDE + "TTG"), 1)
        self.assertEqual(self.d(GUIDE + "TGG", GUIDE + "TAG"), 0)

    def test_a_bulge_column_counts_nothing(self):
        # an RNA bulge (a gap in the target): the column after it counts as without one
        self.assertEqual(self.d("CT-ACAGTTGCTTTTATCACTGG", "CT-gCAGTTGCTTTTATCACTGG"), 1)
        # a DNA bulge (a gap in the guide): its base has no guide base to fit
        self.assertEqual(
            wd.degeneracy("CTAAACAGTTGCTTTTATCACTGG", "CTAAgCAGTTGCTTTTATCACTGG",
                          "CTAA-CAGTTGCTTTTATCAC", "NRG", 0, -3, -3, None), 0)

    def test_a_five_prime_pam(self):
        # TTTV then the protospacer: [4:] and [0:4]
        guide = "ACGTACGTACGTACGTACGTACG"
        self.assertEqual(
            wd.degeneracy("TTTA" + guide, "TTTT" + "AgGTACGTACGTACGTACGTACG", guide, "TTTV",
                          4, None, 0, 4), 2)


class TheCutToTheBudget(unittest.TestCase):
    """twopass_emit.cut_to_degeneracy on the greedy worst case of an all-A guide."""

    GUIDE_NO_PAM = "A" * 20

    def setUp(self):
        # C at columns 3, 9 and 15 (mismatches), each with an alt A that matches
        ref = list("A" * 20 + "AGG")
        for c in (3, 9, 15):
            ref[c] = "C"
        self.ref = "".join(ref)
        self.columns = [
            {"pos_c": c, "candidates": [{"alt": "A", "carriers": {"S%d" % c},
                                         "info": ["rs%d" % c, "0.1", "chrT_%d_C_A" % c]}]}
            for c in (3, 9, 15)]

    def mm(self, seq):
        return twopass_emit.aligned_mismatches(seq, "x" * 23, self.GUIDE_NO_PAM, False, 0, -3,
                                               lambda s: s)

    def deg(self, seq):
        return wd.degeneracy("".join(seq), self.ref, self.GUIDE_NO_PAM, "NGG", 0, -3, -3,
                             None)

    def greedy(self):
        return twopass_emit.greedy_worst_case(self.columns, self.ref, "x" * 23,
                                              self.GUIDE_NO_PAM, False, 0, -3, lambda s: s)

    def test_the_greedy_rests_on_all_three(self):
        rep = self.greedy()
        self.assertEqual(self.mm(rep["seq"]), 0)
        self.assertEqual(self.deg(rep["seq"]), 3)
        self.assertEqual(sorted(rep["chosen"]), [3, 9, 15])

    def test_cut_to_one_column(self):
        rep = self.greedy()
        cut = twopass_emit.cut_to_degeneracy(rep, self.columns, self.ref, 1, self.deg,
                                             lambda s: self.mm(s))
        self.assertEqual(self.deg(cut["seq"]), 1)
        self.assertEqual(self.mm(cut["seq"]), 2)
        # the lowest column goes first on a tie, so the last stays: its carriers and info
        self.assertEqual(cut["carriers"], {"S15"})
        self.assertEqual(cut["info"], [["rs15", "0.1", "chrT_15_C_A"]])
        self.assertEqual({c: (v or {}).get("alt") for c, v in cut["chosen"].items()},
                         {3: None, 9: None, 15: "A"})

    def test_the_rank_chooses(self):
        rep = self.greedy()
        # rank the sequence keeping column 3 first
        cut = twopass_emit.cut_to_degeneracy(rep, self.columns, self.ref, 1, self.deg,
                                             lambda s: (s[3] != "A",))
        self.assertEqual(cut["carriers"], {"S3"})

    def test_within_the_budget_it_is_the_rep_itself(self):
        rep = self.greedy()
        self.assertIs(twopass_emit.cut_to_degeneracy(rep, self.columns, self.ref, 3,
                                                     self.deg, self.mm), rep)
        cut = twopass_emit.cut_to_degeneracy(rep, self.columns, self.ref, 0, self.deg,
                                             self.mm)
        self.assertEqual(cut["seq"], list(self.ref))
        self.assertEqual((cut["carriers"], cut["info"]), (set(), []))


def _site_window(site, codes):
    """``site`` with the letters ``codes`` ({column: letter}) at its columns."""
    window = list(site)
    for col, letter in codes.items():
        window[col] = letter
    return "".join(window)


class _PostAnalysis(unittest.TestCase):
    """A genome of sites, the variant search's rows for them, a classic per-sample
    dictionary, and the real new_simple_analysis.py on them."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="window_density_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.pam = os.path.join(self.tmp, "pam.txt")
        with open(self.pam, "w") as fh:
            fh.write(PAM_NRG)

    def write_inputs(self, sites, rows, entries):
        """``sites``: the sites' reference sequences, laid out between flanks (their 0-based
        starts are returned); ``rows``: (site index, window, strand); ``entries``: the
        dictionary's values by (site index, column)."""
        genome, starts = FLANK, []
        for site in sites:
            starts.append(len(genome))
            genome += site + FLANK
        self.fasta = os.path.join(self.tmp, "chrT.fa")
        with open(self.fasta, "w") as fh:
            fh.write(">chrT\n" + genome + "\n")
        self.dict = os.path.join(self.tmp, "my_dict_chrT.json")
        with open(self.dict, "w") as fh:
            json.dump({"chrT,%d" % (starts[i] + col + 1): value
                       for (i, col), value in entries.items()}, fh)
        self.targets = os.path.join(self.tmp, "chrT.targets.txt")
        with open(self.targets, "w") as fh:
            fh.write("#header\n")
            for i, window, strand in rows:
                fh.write("\t".join(
                    ["X", CRRNA, window, "chrT", str(starts[i]), str(starts[i]), strand,
                     "0", "0", "0", "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"]) + "\n")
        return starts

    def analyse(self, allowed_mms="4", **environ):
        """The rows of .bestCFD.txt, stdout, and the lines of the high-variant-density and
        four-allele BEDs ([] for a BED not written)."""
        out = os.path.join(self.tmp, "out")
        for old in os.listdir(self.tmp):
            if old.startswith("out."):
                os.remove(os.path.join(self.tmp, old))
        env = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        env.update(CRISPRME_SKIP_SCORER_ENV="1", **environ)
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), self.fasta,
             self.dict, self.targets, self.pam, out, allowed_mms],
            cwd=HERE, env=env, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        with open(out + ".bestCFD.txt") as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh.readlines()[1:]]
        beds = []
        for name in ("high_variant_density_regions", "four_allele_positions"):
            path = "%s.%s.bed" % (out, name)
            lines = []
            if os.path.exists(path):
                with open(path) as fh:
                    lines = [line.rstrip("\n").split("\t") for line in fh]
            beds.append(lines)
        return rows, done.stdout, beds[0], beds[1]


# The reference of a site with mismatches at the guide's columns 2, 8 and 14
MISMATCHED = _site_window(GUIDE + "TGG", {2: "G", 8: "C", 14: "A"})


class TheCapCountsBits(_PostAnalysis):
    """Seven three-allele columns: 3^7 = 2^11.09 sequences, over a cap of 10 bits; seven
    positions, under a cap of 10 positions."""

    # (column, reference base, two ALTs): the first ALT matches the guide at 2, 8 and 14
    COLUMNS = ((2, "G", "A", "C"), (8, "C", "T", "G"), (14, "A", "T", "G"),
               (4, "C", "A", "G"), (10, "C", "T", "G"), (16, "T", "A", "C"),
               (18, "A", "C", "T"))
    CODE = {frozenset("ACG"): "V", frozenset("CGT"): "B", frozenset("AGT"): "D",
            frozenset("ACT"): "H"}

    def setUp(self):
        super().setUp()
        window = _site_window(MISMATCHED, {
            col: self.CODE[frozenset((ref, a, b))] for col, ref, a, b in self.COLUMNS})
        entries = {}
        for col, ref, a, b in self.COLUMNS:
            entries[(0, col)] = "S1:1|0;%s,%s;rs%da;0.1$S2:0|1;%s,%s;rs%db;0.1" % (
                ref, a, col, ref, b, col)
        self.start = self.write_inputs([MISMATCHED], [(0, window, "+")], entries)[0]
        self.window = window

    def test_seven_three_allele_columns_are_over_the_cap(self):
        rows, said, hvdr, _q = self.analyse(CRISPRME_FAST_MODE="0")
        self.assertEqual(hvdr[0], ["#chrom", "start", "end", "guide", "n_variants",
                                   "samples_with_alt", "iupac_protospacer", "bits",
                                   "four_allele"])
        self.assertEqual(
            hvdr[1:],
            [["chrT", str(self.start), str(self.start + 23), CRRNA, "7", "S1,S2",
              self.window, "11.09", "0"]])
        self.assertIn("High-variant-density windows of chrT: 1 spell more than 2^10 "
                      "sequences (dense: over CRISPRME_IUPAC_CAP), 0 spell 2^20 or more "
                      "(hyperdense)", said)
        self.assertTrue(rows)
        # a cap of 12 bits holds them
        _rows, said, hvdr, _q = self.analyse(CRISPRME_FAST_MODE="0", CRISPRME_IUPAC_CAP="12")
        self.assertEqual(hvdr[1:], [])
        self.assertNotIn("High-variant-density windows", said)


class AHyperdenseWindow(_PostAnalysis):
    """Twenty two-allele columns and a four-allele PAM column: 2^22 sequences, hyperdense."""

    def setUp(self):
        super().setUp()
        site = GUIDE + "TGG"
        codes = {col: "R" if base in "AG" else "Y" for col, base in enumerate(GUIDE)}
        codes[20] = "Q"
        self.window = _site_window(site, codes)
        # nobody carries the ALTs (the transition of each base)
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}
        entries = {(0, col): ";%s,%s;rs%d;0.1" % (base, alt[base], col)
                   for col, base in enumerate(GUIDE)}
        self.start = self.write_inputs([site], [(0, self.window, "+")], entries)[0]

    def test_it_is_listed_with_its_bits_and_four_allele_position(self):
        _rows, said, hvdr, _q = self.analyse(CRISPRME_FAST_MODE="1")
        self.assertEqual(hvdr[1:], [["chrT", str(self.start), str(self.start + 23), CRRNA,
                                     "20", ".", self.window, "22.00", "1"]])
        self.assertIn("1 spell more than 2^10 sequences (dense: over CRISPRME_IUPAC_CAP), 1 "
                      "spell 2^20 or more (hyperdense)", said)

    def test_whatever_the_cap(self):
        for cap in ("-1", "30"):
            with self.subTest(cap=cap):
                _rows, said, hvdr, _q = self.analyse(CRISPRME_FAST_MODE="1",
                                                     CRISPRME_IUPAC_CAP=cap)
                self.assertEqual([line[7] for line in hvdr[1:]], ["22.00"])
                self.assertIn("1 spell 2^20 or more (hyperdense)", said)


class TheFourAllelePositions(_PostAnalysis):
    def setUp(self):
        super().setUp()
        site = GUIDE + "TGG"
        # a Q at the PAM's N beside an R, and a minus-strand window holding a Q alone (its
        # fourth letter, read from the PAM's end: the 20th base of the site from its start)
        with_r = _site_window(site, {20: "Q", 12: "Y"})
        minus = _site_window(site, {3: "Q"})
        entries = {(0, 12): "S1:1|0;T,C;rs12;0.1"}
        self.starts = self.write_inputs([site, site], [(0, with_r, "+"), (1, minus, "-")],
                                        entries)

    def test_they_are_listed_0_based_on_the_forward_strand(self):
        _rows, said, _hvdr, q = self.analyse(CRISPRME_FAST_MODE="0")
        self.assertEqual(q, [["#chrom", "start", "end"],
                             ["chrT", str(self.starts[0] + 20), str(self.starts[0] + 21)],
                             ["chrT", str(self.starts[1] + 19), str(self.starts[1] + 20)]])
        self.assertIn("Four-allele positions (Q: every base is an observed allele there) in "
                      "the windows of chrT: 2", said)
        # with none, no file
        with open(self.targets) as fh:
            text = fh.read().replace("Q", "A")
        with open(self.targets, "w") as fh:
            fh.write(text)
        _rows, said, _hvdr, q = self.analyse(CRISPRME_FAST_MODE="0")
        self.assertEqual(q, [])
        self.assertNotIn("Four-allele positions", said)


class MaxDegeneracy(_PostAnalysis):
    """Three mismatches, each a two-allele column whose ALT matches the guide: rows resting
    on 1, 2 or 3 of them."""

    def write(self, gt):
        window = _site_window(MISMATCHED, {2: "R", 8: "Y", 14: "W"})
        entries = {(0, 2): "S1:%s;G,A;rs2;0.1" % gt, (0, 8): "S1:%s;C,T;rs8;0.1" % gt,
                   (0, 14): "S1:%s;A,T;rs14;0.1" % gt}
        self.write_inputs([MISMATCHED], [(0, window, "+")], entries)

    def degeneracies(self, rows):
        return sorted(_degeneracy_of_row(r) for r in rows if r[18] != "n")

    def test_the_worst_case_is_cut_back(self):
        self.write("0/1")
        # off by default: the worst case rests on all three
        rows, said, _hvdr, _q = self.analyse(CRISPRME_FAST_MODE="1")
        self.assertIn(3, self.degeneracies(rows))
        self.assertNotIn("--max-degeneracy", said)
        rows, said, _hvdr, _q = self.analyse(CRISPRME_FAST_MODE="1", CRISPRME_MAX_DEGENERACY="1")
        self.assertTrue(rows)
        self.assertEqual(set(self.degeneracies(rows)), {1})
        # the min-mismatch one: two mismatches left
        self.assertIn("2", [r[7] for r in rows])
        self.assertIn("--max-degeneracy 1 on chrT: 1 worst-case alignment(s) cut back to rest "
                      "on at most 1 ambiguous column(s); 0 row(s) resting on more, which no "
                      "phased carrier confirms, withheld", said)

    def test_an_unphased_combination_resting_on_more_is_withheld(self):
        self.write("0/1")
        rows, _said, _hvdr, _q = self.analyse(CRISPRME_FAST_MODE="0")
        self.assertEqual(self.degeneracies(rows), [1, 1, 1, 2, 2, 2, 3])
        rows, said, _hvdr, _q = self.analyse(CRISPRME_FAST_MODE="0", CRISPRME_MAX_DEGENERACY="1")
        self.assertEqual(self.degeneracies(rows), [1, 1, 1])
        self.assertIn("0 worst-case alignment(s) cut back to rest on at most 1 ambiguous "
                      "column(s); 4 row(s) resting on more", said)

    def test_a_phased_haplotype_keeps_its_row(self):
        self.write("1|0")
        rows, said, _hvdr, _q = self.analyse(CRISPRME_FAST_MODE="0", CRISPRME_MAX_DEGENERACY="1")
        self.assertEqual([(r[7], r[12], _degeneracy_of_row(r)) for r in rows], [("0", "S1", 3)])
        self.assertIn("0 row(s) resting on more", said)


class _GtSentinel(object):
    """Selects the tier path's observed haplotypes (mygt present)."""

    def __len__(self):
        return 1


class MaxDegeneracyOnTheTierPath(unittest.TestCase):
    """The observed haplotypes of the dict-less tiers: the guide is all A, the reference has a
    C at columns 0 to 2, and each ALT A matches."""

    GUIDE = "A" * 20 + "NGG"

    def run_path(self, gt, fast=False, max_degeneracy=1):
        ref = "CCC" + "A" * 17 + "AGG"
        split = ["X", self.GUIDE, "WWW" + ref[3:], "chrT", "0", "0", "+", "3", "0", "3",
                 "NGG", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        ns = _load_pure_functions(dict(
            genomeStr=ref, current_chr="chrT", myreg=None, mygt=_GtSentinel(),
            mydict={"chrT,%d" % (i + 1): "S:%s;C,A;rs%d;0.10" % (gt, i) for i in range(3)},
            haplotype_check=True, IUPAC_CAP=10, hvdr_bed=io.StringIO(), pam="NGG",
            pos_beg=0, pos_end=-3, pam_begin=-3, pam_end=None, allowed_mms=4,
            _FAST_MODE=fast, MAX_DEGENERACY=max_degeneracy))
        ns["_degeneracy_counts"].update(cut=0, withheld=0)
        cluster = []
        ns["iupac_decomposition"](split, self.GUIDE, "A" * 20, cluster)
        variant = [r for r in cluster if r[-2] != 55]
        return ([wd.degeneracy(r[2], r[-3], "A" * 20, "NGG", 0, -3, -3, None)
                 for r in variant], ns["_degeneracy_counts"])

    def test_a_confirmed_haplotype_keeps_its_row(self):
        self.assertEqual(self.run_path("1|0"), ([3], {"cut": 0, "withheld": 0}))

    def test_a_putative_combination_resting_on_more_is_withheld(self):
        d, counts = self.run_path("0/1")
        self.assertEqual(sorted(d), [1, 1, 1])
        self.assertEqual(counts, {"cut": 0, "withheld": 4})
        d, counts = self.run_path("0/1", max_degeneracy=-1)
        self.assertEqual(sorted(d), [1, 1, 1, 2, 2, 2, 3])

    def test_the_population_level_worst_case_is_cut_back(self):
        d, counts = self.run_path("1|0", fast=True)
        self.assertEqual(d, [1])
        self.assertEqual(counts["cut"], 1)
        d, _counts = self.run_path("1|0", fast=True, max_degeneracy=-1)
        self.assertEqual(d, [3])


# intervaltree is not among the light CI deps: a stand-in that answers as it does
INTERVALTREE_STUB = '''import collections

Interval = collections.namedtuple("Interval", "begin end data")


class IntervalTree:
    """Stand-in for the test: the intervals a point lies in, or a range overlaps."""

    def __init__(self):
        self.items = []

    def __setitem__(self, key, value):
        self.items.append(Interval(key.start, key.stop, value))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return {i for i in self.items if i.begin < key.stop and key.start < i.end}
        return {i for i in self.items if i.begin <= key < i.end}

    def __len__(self):
        return len(self.items)
'''


def _gnu_cut():
    """resultIntegrator drops its unused columns with GNU ``cut --complement``."""
    try:
        res = subprocess.run(["cut", "--complement", "-f1"], input=b"a\tb\n",
                             capture_output=True)
    except OSError:
        return False
    return res.returncode == 0


@unittest.skipUnless(_gnu_cut(), "resultIntegrator.py needs GNU cut --complement")
class TheIntegratedRow(unittest.TestCase):
    """The real resultIntegrator.py on bestMerge rows at chrT 1000, 2000, 3000 and 4000,
    beside the post-analysis's BEDs: a dense window over 1000 (2.9.0's BED), one over 4000
    (2.8.1's), and four-allele positions in the windows at 1000 and 2000, and just past the
    one at 3000."""

    ALN = GUIDE + "TGG"

    @classmethod
    def setUpClass(cls):
        from test_seed_counts_dna_bulge import _block
        cls.d = tempfile.mkdtemp(prefix="window_density_ri_")
        cls.addClassCleanup(shutil.rmtree, cls.d, ignore_errors=True)
        with open(os.path.join(cls.d, "intervaltree.py"), "w") as fh:
            fh.write(INTERVALTREE_STUB)
        best = os.path.join(cls.d, "job.bestMerge.txt")
        starts = (1000, 2000, 3000, 4000)
        with open(best, "w") as fh:
            for start in starts:
                block = _block("X", CRRNA, CRRNA, cls.ALN, cls.ALN.lower(), 0, 0, "S1",
                               "chrT_%d_A_G" % (start + 3))
                block[5] = block[6] = str(start)
                fh.write("\t".join(block * 3) + "\n")
        with open(os.path.join(cls.d, "job_chrT.high_variant_density_regions.bed"), "w") as fh:
            fh.write("#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\t"
                     "iupac_protospacer\tbits\tfour_allele\n"
                     "chrT\t1000\t1023\t%s\t24\tS1\tRRRRRRRRRRRRRRRRRRRRQRR\t24.00\t1\n" % CRRNA)
        with open(os.path.join(cls.d, "old_chrT.high_variant_density_regions.bed"), "w") as fh:
            fh.write("#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\t"
                     "iupac_protospacer\n"
                     "chrT\t4000\t4023\t%s\t12\tS1\tRRRRRRRRRRRRAAAAAAAATGG\n" % CRRNA)
        with open(os.path.join(cls.d, "job_chrT.four_allele_positions.bed"), "w") as fh:
            fh.write("#chrom\tstart\tend\n"
                     "chrT\t1020\t1021\nchrT\t2010\t2011\nchrT\t3023\t3024\n")
        found = os.path.join(cls.d, "job.found.bed")
        with open(found, "w") as fh:
            fh.write("NA\n" * len(starts))
        empirical = os.path.join(cls.d, "dummy.txt")
        open(empirical, "w").close()
        guides = os.path.join(cls.d, "guides.txt")
        with open(guides, "w") as fh:
            fh.write(CRRNA + "\n")
        env = dict(os.environ, PYTHONPATH=cls.d)
        res = subprocess.run(
            [sys.executable, os.path.join(HERE, "resultIntegrator.py"), best, empirical,
             found, guides, cls.d + "/", "true", "hg38", "vuota", "1", "0", "none"],
            cwd=HERE, env=env, capture_output=True, text=True)
        if res.returncode or res.stderr:
            raise AssertionError("resultIntegrator.py failed (stderr is fatal in the "
                                 "pipeline): rc=%s\n%s" % (res.returncode, res.stderr))
        with open(best + ".integrated_results.tsv", newline="") as fh:
            cls.notes = {int(r["Start_coordinate_(highest_CFD)"]): r["High_variant_density_region"]
                         for r in csv.DictReader(fh, delimiter="\t")}

    def test_a_dense_window_with_its_sequences_and_four_allele_position(self):
        self.assertEqual(
            self.notes[1000],
            "high_variant_density (24 variants, 2^24.00 sequences): a greedy worst-case "
            "alignment is reported here; additional haplotype alignments may exist; "
            "full_IUPAC=RRRRRRRRRRRRRRRRRRRRQRR; four_allele_positions=chrT:1021")
        # a BED of 2.8.1, which gave no bits
        self.assertEqual(
            self.notes[4000],
            "high_variant_density (12 variants): a greedy worst-case alignment is reported "
            "here; additional haplotype alignments may exist; "
            "full_IUPAC=RRRRRRRRRRRRAAAAAAAATGG")

    def test_a_four_allele_position_in_the_window_and_not_past_it(self):
        self.assertEqual(self.notes[2000], "four_allele_positions=chrT:2011")
        self.assertEqual(self.notes[3000], "NA")


class TheReport(unittest.TestCase):
    NOTE = ("high_variant_density (%s variants%s): a greedy worst-case alignment is reported "
            "here; additional haplotype alignments may exist; full_IUPAC=ACGT")
    COLS = {"complex_region": "High_variant_density_region",
            "var_genome": "Variant_info_genome_(highest_CFD)"}

    def cell(self, note):
        return gr._curated_cell("complex_region", {"High_variant_density_region": note},
                                self.COLS)

    def test_the_dense_region_column(self):
        self.assertEqual(self.cell(self.NOTE % ("12", "")), "12 in window")
        self.assertEqual(self.cell(self.NOTE % ("12", ", 2^19.99 sequences")), "12 in window")
        self.assertEqual(self.cell(self.NOTE % ("24", ", 2^24.00 sequences")),
                         "24 in window (2^24.0 sequences)")
        self.assertEqual(
            self.cell(self.NOTE % ("24", ", 2^24.00 sequences")
                      + "; four_allele_positions=chrT:1021"),
            "24 in window (2^24.0 sequences); 1 four-allele position")
        self.assertEqual(self.cell("four_allele_positions=chrT:2011,chrT:2015"),
                         "2 four-allele positions")
        self.assertEqual(self.cell("NA"), gr.CURATED_MISSING)

    def test_a_four_allele_position_alone_is_no_dense_window(self):
        """A --per-sample report of untiered datasets reads a dense window's row of two or
        more SNPs as its greedy worst-case alignment (_greedy_union): not a row whose
        window only holds a four-allele position."""
        saved = gr._PER_SAMPLE, gr._UNTIERED
        self.addCleanup(setattr, gr, "_UNTIERED", saved[1])
        self.addCleanup(setattr, gr, "_PER_SAMPLE", saved[0])
        gr._PER_SAMPLE = gr._UNTIERED = True
        row = {"Variant_info_genome_(highest_CFD)": "chrT_1006_A_G,chrT_1010_C_T",
               "High_variant_density_region": "four_allele_positions=chrT:1006"}
        self.assertFalse(gr._greedy_union(row, self.COLS))
        row["High_variant_density_region"] = self.NOTE % ("12", "")
        self.assertTrue(gr._greedy_union(row, self.COLS))

    def test_the_inputs_card_counts_the_rows_and_says_the_cap(self):
        hyperdense = self.NOTE % ("24", ", 2^24.00 sequences")
        df = pd.DataFrame({"High_variant_density_region": [
            hyperdense, hyperdense + "; four_allele_positions=chrT:1021",
            "four_allele_positions=chrT:5",
            self.NOTE % ("12", ""), "NA"]})
        self.assertEqual(gr._density_counts(df, self.COLS), (2, 2))
        base = {"datasets": "HGDP1kGP", "mm": "4", "bdna": 1, "brna": 1, "max_edits": "6",
                "search_mode": "population-level"}
        html = gr.render_inputs_criteria({**base, "n_hyperdense_rows": 2,
                                          "n_four_allele_rows": 1, "max_degeneracy": "2"})
        self.assertIn("Hyperdense and four-allele windows", html)
        self.assertIn("2 reported off-targets lie in a hyperdense window, one spelling "
                      "2<sup>20</sup> sequences or more", html)
        self.assertIn("1 reported off-target lies in a window holding a four-allele "
                      "position", html)
        self.assertIn("A SNP variant row rests on at most 2 ambiguous columns "
                      "(<code>--max-degeneracy 2</code>) unless its carriers hold its variants "
                      "together on one phased haplotype.", html)
        self.assertIn("seldom exists in any one genome", html)
        plain = gr.render_inputs_criteria(base)
        self.assertNotIn("Hyperdense and four-allele windows", plain)
        self.assertNotIn("--max-degeneracy", plain)
        # a per-sample row there names the individuals who carry it: no worst case to doubt
        per_sample = gr.render_inputs_criteria({**base, "search_mode": "per-sample",
                                                "n_hyperdense_rows": 1})
        self.assertIn("1 reported off-target lies in a hyperdense window, one spelling "
                      "2<sup>20</sup> sequences or more", per_sample)
        self.assertNotIn("seldom exists in any one genome", per_sample)

    def test_the_cap_is_stated_with_its_exception(self):
        """--max-degeneracy keeps a row whose carriers hold its variants together on one
        phased haplotype whatever the columns it rests on (a per-sample report's CONFIRMED
        row: MaxDegeneracy.test_a_phased_haplotype_keeps_its_row keeps one resting on 3 under
        a cap of 1). The search criteria state the cap with that exception in the same
        sentence, never as a bound every SNP variant row keeps."""
        base = {"datasets": "HGDP1kGP", "mm": "4", "bdna": 1, "brna": 1, "max_edits": "6",
                "max_degeneracy": "1"}
        for mode in ("per-sample", "population-level"):
            with self.subTest(mode=mode):
                html = gr.render_inputs_criteria({**base, "search_mode": mode})
                at = html.index("rests on at most 1 ambiguous column")
                sentence = html[html.rfind(". ", 0, at) + 2:html.index(". ", at) + 1]
                self.assertEqual(
                    sentence,
                    "A SNP variant row rests on at most 1 ambiguous column (<code>"
                    "--max-degeneracy 1</code>) unless its carriers hold its variants "
                    "together on one phased haplotype.")
                self.assertNotIn("Each SNP variant row", html)

    def test_the_legend_says_both(self):
        """One word each: a window is dense over CRISPRME_IUPAC_CAP, as since 2.8.1, and
        hyperdense at 2^HYPERDENSE_BITS sequences or more, flagged whatever the cap."""
        legend = gr.build_annotation_legend_html()
        self.assertIn("A window is dense when it spells more sequences than "
                      "<code>CRISPRME_IUPAC_CAP</code> allows (2<sup>10</sup> by default",
                      legend)
        self.assertIn("2<sup>20</sup> or more is hyperdense, and is flagged whatever that cap "
                      "is", legend)
        self.assertIn("<code>1 four-allele position</code>", legend)


class TheCommandLine(bio._Fixture):
    """complete-search --max-degeneracy, as a user runs it, against
    test_build_index_only_atomic's stub engine, up to the search it launches."""

    def params(self):
        with open(self.path("Results", "search", "Params.txt")) as fh:
            return dict(line.rstrip("\n").split("\t", 1) for line in fh if "\t" in line)

    def test_it_is_written_to_the_run_s_params(self):
        self.run_bio()
        rc, out, err = self.run_search("--max-degeneracy", "2")
        self.assertIn("Launching job ", out, out + err)
        self.assertEqual(self.params().get("Max_degeneracy"), "2")
        rc, out, err = self.run_search()
        self.assertIn("Launching job ", out, out + err)
        self.assertNotIn("Max_degeneracy", self.params())

    def test_it_takes_a_count(self):
        for value, said in (("-1", "--max-degeneracy must be a non-negative integer"),
                            ("two", "Please provide a non-negative integer for "
                                    "--max-degeneracy")):
            with self.subTest(value=value):
                rc, out, err = self.run_search("--max-degeneracy", value)
                self.assertNotEqual(rc, 0)
                self.assertIn(said, err)

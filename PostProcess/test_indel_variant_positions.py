#!/usr/bin/env python
"""Variant_info_spacer+PAM on indel rows (resultIntegrator.py, 2.8.2), through the REAL
integrator on bestMerge rows.

A SNP row's variants are numbered by their distance from the row's start along the
reference: POS - start on the + strand, len(aligned spacer) + 1 - (POS - start) on the -
strand. An indel row's target is bases of the haplotype carrying its indel, which holds
an insertion's bases and not a deletion's, and its variants are numbered the same way
along that haplotype (indel_snp_cis.variant_places). Up to 2.8.1 they were numbered along
the reference: a SNP past an insertion was len(ALT) - len(REF) bases short (long on the -
strand), one past a deletion as many long (short), and on a row starting on a base an
insertion adds (which 2.8.2 places at the insertion's anchor) shorter again by the bases
of the insertion the target does not hold. SNP rows, and an indel row's own indel, keep
their numbers.

Needs GNU cut (resultIntegrator drops its unused columns with ``cut --complement``);
intervaltree is stubbed. Run: python -m unittest test_indel_variant_positions -v
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_seed_counts_dna_bulge import INTERVALTREE_STUB, _gnu_cut

HERE = os.path.dirname(os.path.abspath(__file__))
G = "ACTGAAATCTGTAAGCAGGC"


def _block(start, strand, crrna, dna, reference, variants, bulge="X", mm=0, bulges=0):
    """One 24-column block of a bestMerge row of alt origin."""
    n_var = len(variants.split(","))
    return [bulge, crrna, dna, reference, "chrT", str(start), str(start), strand,
            str(mm), str(bulges), str(mm + bulges), "NA", "y", "S1", "CDS_gencode",
            crrna.replace("-", ""), ",".join(["rs"] * n_var), ",".join(["0.1"] * n_var),
            variants, "0", "0.5", "0.4", "0.1", "0.1"]


def _rc(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


# name -> (block, expected Variant_info_spacer+PAM, what 2.8.1 wrote)
def _cases():
    cases = {}
    # Z3 of the SNP+indel oracle (+): AAG + the guide's first six bases inserted after the
    # anchor T at 0-based 99; the target (guide + TGG) starts at the 4th added base, the
    # row at the anchor; a SNP past the insertion holds the guide's 13th base. It is the
    # target's 13th base (2.8.1, at the haplotype place 103: 4; the start alone at 99: 8).
    ins = "AAG" + G[:6]
    ref = "N" * 99 + "T" + G[6:12] + "C" + G[13:] + "TGG" + "GGGGGG"
    cases["inside_insertion_plus"] = (
        _block(99, "+", G + "NNN", G + "TGG", ref[99:122].lower(),
               "chrT_100_T_T%s,chrT_107_C_%s" % (ins, G[12])),
        "1T>T%s,13C>%s" % (ins, G[12]), "1T>T%s,4C>%s" % (ins, G[12]))
    # chr22:21052639 (-) of the HGDP + 1000 Genomes README search: the SNP is the 22nd
    # column on the - strand (2.8.1 at 21052642: clamped to the last, 24; at 21052639: 23)
    cases["inside_insertion_minus"] = (
        _block(21052639, "-", "CTAACAGT-TGCTTTTATCACNNN", "C-AACAGTGaGCaTgTgTCACTGG",
               "C-cAacagTgagcaTgtgtcaCCA", "chrT_21052640_T_TCACC,chrT_21052641_G_A",
               bulge="RNA,DNA", mm=4, bulges=2),
        "24A>GGTGA,22C>T", "24A>GGTGA,24C>T")
    # a deletion (+): the guide on the reference with a C after its 10th base, which the
    # carrier deletes, and a SNP at the PAM's last base past it (the oracle's site D): the
    # 23rd base, not the 24th (a number past the target's end)
    d_ref = "N" * 1000 + G[:10] + "C" + G[10:] + "TGA" + "TTT"
    cases["snp_past_a_deletion_plus"] = (
        _block(1000, "+", G + "NNN", G + "TGG", d_ref[1000:1023].lower(),
               "chrT_1010_%sC_%s,chrT_1024_A_G" % (G[9], G[9])),
        "10%sC>%s,23A>G" % (G[9], G[9]), "10%sC>%s,24A>G" % (G[9], G[9]))
    # the same on the - strand (the oracle's site R): the SNP, forward right of the
    # deletion, is the target's 6th base (2.8.1: 5)
    site = G[:5] + "C" + G[6:10] + "C" + G[10:] + "TGG"
    fwd = _rc(site)
    r_ref = "N" * 2000 + fwd + "TTTT"
    cases["snp_past_a_deletion_minus"] = (
        _block(2000, "-", G + "NNN", G + "TGG", _rc(r_ref[2000:2023]).lower(),
               "chrT_2013_%s_%s,chrT_2019_G_T" % (fwd[12:14], fwd[12])),
        "10%s>%s,6C>A" % (_rc(fwd[12:14]), _rc(fwd[12])),
        "10%s>%s,5C>A" % (_rc(fwd[12:14]), _rc(fwd[12])))
    # an insertion in the protospacer (+, row starting before it) and a SNP past it: 2 on
    p_ref = "N" * 3000 + G[:10] + G[13:15] + "T" + G[16:] + "TGGAAA"
    cases["snp_past_an_insertion_plus"] = (
        _block(3000, "+", G + "NNN", G + "TGG", p_ref[3000:3023].lower(),
               "chrT_3010_%s_%s,chrT_3013_T_%s" % (G[9], G[9:13], G[15])),
        "10%s>%s,16T>%s" % (G[9], G[9:13], G[15]), "10%s>%s,13T>%s" % (G[9], G[9:13], G[15]))
    # SNP rows keep their numbers: + strand, and - strand with an RNA bulge
    s_ref = "N" * 4000 + G + "AGG"
    cases["snp_row_plus"] = (
        _block(4000, "+", G + "NNN", G[:5] + "t" + G[6:] + "AGG", s_ref[4000:4023].lower(),
               "chrT_4006_%s_T" % G[5], mm=1),
        "6%s>T" % G[5], "6%s>T" % G[5])
    cases["snp_row_minus_rna_bulge"] = (
        _block(5000, "-", G + "NNN", G[:8] + "-" + G[9:15] + "a" + G[16:] + "TGG",
               G[:8] + "-" + G[9:] + "TGG", "chrT_5003_G_T", bulge="RNA", mm=1, bulges=1),
        "21C>A", "21C>A")
    return cases


@unittest.skipUnless(_gnu_cut(), "resultIntegrator.py needs GNU cut --complement")
class IndelVariantPositions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="indel_variant_positions_")
        cls.addClassCleanup(shutil.rmtree, cls.d, ignore_errors=True)
        with open(os.path.join(cls.d, "intervaltree.py"), "w") as fh:
            fh.write(INTERVALTREE_STUB)
        cls.cases = _cases()
        cls.names = sorted(cls.cases)
        best = os.path.join(cls.d, "job.bestMerge.txt")
        with open(best, "w") as fh:
            for name in cls.names:
                fh.write("\t".join(cls.cases[name][0] * 3) + "\n")
        with open(os.path.join(cls.d, "job.found.bed"), "w") as fh:
            fh.write("NA\n" * len(cls.names))
        open(os.path.join(cls.d, "dummy.txt"), "w").close()
        with open(os.path.join(cls.d, "guides.txt"), "w") as fh:
            fh.write(G + "NNN\n")
        res = subprocess.run(
            [sys.executable, os.path.join(HERE, "resultIntegrator.py"), best,
             os.path.join(cls.d, "dummy.txt"), os.path.join(cls.d, "job.found.bed"),
             os.path.join(cls.d, "guides.txt"), cls.d + "/", "true", "hg38", "vuota",
             "1", "0", "none"],
            cwd=HERE, env=dict(os.environ, PYTHONPATH=cls.d), capture_output=True,
            text=True)
        assert res.returncode == 0 and res.stderr == "", res.stdout[-2000:] + res.stderr
        with open(best + ".integrated_results.tsv", newline="") as fh:
            cls.rows = list(csv.DictReader(fh, delimiter="\t"))

    def test_every_block_counts_along_the_target(self):
        self.assertEqual(len(self.rows), len(self.names))
        for name, row in zip(self.names, self.rows):
            want, old = self.cases[name][1:]
            for crit in ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE"):
                self.assertEqual(row["Variant_info_spacer+PAM_(%s)" % crit], want,
                                 "%s %s (2.8.1 wrote %s)" % (name, crit, old))

    def test_an_insertion_keeps_its_place_when_the_row_starts_inside_it(self):
        """+ strand: position 1; - strand: the last, len(aligned spacer)."""
        by = dict(zip(self.names, self.rows))
        self.assertTrue(by["inside_insertion_plus"]["Variant_info_spacer+PAM_(highest_CFD)"]
                        .startswith("1T>T"))
        self.assertTrue(by["inside_insertion_minus"]["Variant_info_spacer+PAM_(highest_CFD)"]
                        .startswith("24A>"))

    def test_the_cases_change_where_2_8_1_was_wrong_only(self):
        changed = sorted(n for n in self.names if self.cases[n][1] != self.cases[n][2])
        self.assertEqual(changed, ["inside_insertion_minus", "inside_insertion_plus",
                                   "snp_past_a_deletion_minus", "snp_past_a_deletion_plus",
                                   "snp_past_an_insertion_plus"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

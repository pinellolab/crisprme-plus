"""Deterministic representative TIE-BREAK regression (issue #139).

CONTEXT
-------
``remove_contiguous_samples.get_best_targets`` picks each off-target CLUSTER's
representative row with ~30 ``sorted(..., key=lambda x: (-float(x[cfd]),
int(x[total...]), int(x[pos])))`` calls (score branch) plus the ~30 mirror calls in
the mm+bulges branch. The FINAL pre-existing key element ``int(x[pos])`` is the
Cluster_Position column, which is CONSTANT for every row in a cluster. So when two
alignments tie EXACTLY on CFD and on mm+bulges, the winner fell to Python's
stable-sort INPUT order -- i.e. the alt/enumeration order. The dict path (VCF order)
and the dict-less path (sorted) enumerate a locus's alts in DIFFERENT orders, so for
the ~12 exactly-tied clusters they picked a DIFFERENT representative alignment even
though the cluster + carriers + CFD were identical.

THE FIX (under test)
--------------------
Every sort key now APPENDS an intrinsic final element ``*_rep_tiebreak(x)`` =
``(most carriers, aligned-target-DNA at col2, true genomic Position at col pos-1,
has-indel)``. Because it is consulted only when all earlier key elements tie exactly,
it CANNOT reorder any cluster whose best row is already unique (output byte-identical
for those); for exactly-tied clusters it makes both paths converge on the
lexicographically-smallest aligned-DNA representative.

2.8.1: a row of SNPs and a row with an indel can give the SAME aligned DNA at the same
position, with the same number of carriers. Those two fell through to input order, the
whole-line sort of the merge step, where the REF alignment's letter case and the sample
IDs decide. The has-indel element settles it: the row of SNPs is the representative.

WHAT THIS TEST DOES
-------------------
Drives the REAL ``remove_contiguous_samples.py`` (pure stdlib) via subprocess on an
on-disk POST-ADJUST bestCFD file (the exact layout the merge consumes, 22 cols, DNA at
col2, Cluster_Position at col5, SNP at col18, CFD at col20), with the same 1-based
indices ``merge_close_targets_cfd.sh`` passes (chrom=5 pos=7 total=11 true_guide=16
snp_info=19 cfd=21, sort_order=score, criteria=mm+bulges).

ASSERTS:
  * EXACT tie (same CFD, same total, same Cluster_Position, DIFFERENT aligned DNA):
    feeding the two rows in BOTH input orders yields the SAME representative in the
    best file -- and it is the lexicographically-SMALLER aligned-DNA row. Pre-fix, the
    representative flipped with input order.
  * NON-tie control (different CFD): the higher-CFD row is the representative
    regardless of input order -- the new key never overrides a real CFD difference.
  * SNP/indel tie (same aligned DNA, position, CFD, total, number of carriers): the row
    of SNPs is the representative in BOTH input orders and for all the variant forms
    (deletion, insertion, indel + SNP, a contig name with '_'), and the row with the
    indel is listed in the alternative alignments, on the score and on the total
    branch. Controls: more carriers, a smaller aligned DNA, a smaller Position or a
    higher CFD still win over it.

Run with:
    cd PostProcess && python3 -m unittest test_representative_tiebreak -v
"""

import os
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))

# POST-ADJUST layout the merge consumes (see adjust_cols.py / merge_close_targets_cfd.sh):
#   col2  = DNA (aligned target sequence; only a row of SNPs and a row with an indel can
#           share it at the same position within a cluster)
#   col4  = Chromosome
#   col5  = Cluster_Position (CONSTANT within a cluster -> the pre-fix final key element)
#   col10 = Total (mm+bulges)
#   col15 = Real_Guide (true guide; cluster is keyed on this + chrom + pos-window)
#   col18 = SNP  ("n" == reference row)
#   col20 = CFD
# 1-based indices merge_close_targets_cfd.sh passes to remove_contiguous_samples.py.
POST_ADJUST_HEADER = [
    "#Bulge_type", "crRNA", "DNA", "Reference", "Chromosome", "Position",
    "Cluster_Position", "Direction", "Mismatches", "Bulge_Size", "Total",
    "PAM_gen", "Var_uniq", "Samples", "Annotation_Type", "Real_Guide", "rsID",
    "AF", "SNP", "#Seq_in_cluster", "CFD", "CFD_ref",
]
DNA_COL = 2
CLUSTER_POS_COL = 6
TOTAL_COL = 10
SNP_COL = 18
CFD_COL = 20

GUIDE = "GAGTCCGAGCAGAAGAAGAANNN"


def _row(dna, cfd, total=4, cluster_pos=1000, position=1000, snp="n",
         chrom="chr1", guide=GUIDE, samples="n"):
    """Build one 22-col POST-ADJUST reference row (SNP=='n')."""
    r = ["X"] * 22
    r[0] = "X"                       # Bulge_type
    r[1] = guide                     # crRNA
    r[DNA_COL] = dna                 # aligned target DNA (the discriminator)
    r[3] = "n"                       # Reference marker (ref row)
    r[4] = chrom                     # Chromosome
    r[5] = str(position)             # true genomic Position (col pos-1 == 5)
    r[CLUSTER_POS_COL] = str(cluster_pos)  # Cluster_Position (CONSTANT in a cluster)
    r[7] = "+"                       # Direction
    r[8] = "0"                       # Mismatches
    r[9] = "0"                       # Bulge_Size
    r[TOTAL_COL] = str(total)        # Total (mm+bulges)
    r[11] = "NGG"                    # PAM_gen
    r[12] = "n"                      # Var_uniq
    r[13] = samples                  # Samples roster (default "n" == ref)
    r[14] = "NA"                     # Annotation_Type
    r[15] = guide.replace("-", "")   # Real_Guide (true guide)
    r[16] = "n"                      # rsID
    r[17] = "n"                      # AF
    r[SNP_COL] = snp                 # SNP ("n" -> reference row)
    r[19] = "0"                      # #Seq_in_cluster
    r[CFD_COL] = str(cfd)            # CFD
    r[21] = str(cfd)                 # CFD_ref
    return r


class RepresentativeTiebreak(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="rep_tiebreak_")

    def _run_merge(self, rows, sort_order="score"):
        """Write rows to a POST-ADJUST bestCFD file, run the REAL
        remove_contiguous_samples.py, return the list of data rows in the BEST file."""
        infile = os.path.join(self.d, "in.bestCFD.txt")
        outfile = os.path.join(self.d, "out.bestCFD.txt")
        with open(infile, "w") as fh:
            fh.write("\t".join(POST_ADJUST_HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        cmd = [sys.executable, os.path.join(PP, "remove_contiguous_samples.py"),
               infile, outfile, "3", "5", "7", "11", "16", "19", "21",
               sort_order, "mm+bulges"]
        res = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, "merge failed: %s" % res.stderr)
        with open(outfile) as fh:
            lines = [l.rstrip("\n") for l in fh]
        return [l.split("\t") for l in lines[1:] if l.strip()]

    def _representative(self, rows, sort_order="score"):
        """Single-cluster input -> the one representative row in the best file."""
        best = self._run_merge(rows, sort_order)
        self.assertEqual(len(best), 1,
                         "expected exactly ONE representative for the single cluster, "
                         "got %d" % len(best))
        return best[0]

    def _alternatives(self):
        """The data rows the last merge listed as alternative alignments."""
        with open(os.path.join(self.d, "out.bestCFD.txt.discarded_samples")) as fh:
            lines = [l.rstrip("\n") for l in fh]
        return [l.split("\t") for l in lines[1:] if l.strip()]

    def _assert_representative_both_orders(self, winner, loser, why):
        """``winner`` is the representative whichever row comes first, on the score
        branch (highest_CFD, highest_CRISPR_BULGE) and on the total branch (fewest_mm+b),
        and ``loser`` is the one alternative alignment."""
        for sort_order in ("score", "total"):
            for rows in ((winner, loser), (loser, winner)):
                rep = self._representative([list(r) for r in rows], sort_order)
                self.assertEqual(
                    rep[SNP_COL], winner[SNP_COL],
                    "%s: on the %s branch with %s first the representative is %s"
                    % (why, sort_order, rows[0][SNP_COL], rep[SNP_COL]))
                alts = self._alternatives()
                self.assertEqual([a[SNP_COL] for a in alts], [loser[SNP_COL]])

    def test_exact_tie_is_input_order_independent(self):
        """Two rows tie EXACTLY on CFD, on Total (mm+bulges) and on Cluster_Position but
        have DIFFERENT aligned DNA. The representative must be the SAME (and the
        smaller-DNA one) regardless of input order. Pre-fix it flipped with the order."""
        dna_small = "AAAATCCGAGCAGAAGAAGAAGG"
        dna_large = "TTTTTCCGAGCAGAAGAAGAAGG"
        self.assertLess(dna_small, dna_large)
        row_small = _row(dna_small, cfd="0.900", total=4, cluster_pos=1000, position=1001)
        row_large = _row(dna_large, cfd="0.900", total=4, cluster_pos=1000, position=1002)

        rep_ab = self._representative([row_small, row_large])
        rep_ba = self._representative([row_large, row_small])

        self.assertEqual(
            rep_ab[DNA_COL], rep_ba[DNA_COL],
            "representative flipped with input order: %r vs %r"
            % (rep_ab[DNA_COL], rep_ba[DNA_COL]))
        self.assertEqual(
            rep_ab[DNA_COL], dna_small,
            "the representative must be the lexicographically-smallest aligned DNA "
            "(%r), got %r" % (dna_small, rep_ab[DNA_COL]))

    def test_non_tie_best_cfd_still_wins(self):
        """Control: with DIFFERENT CFD the new key must NOT override the real CFD
        difference -- the higher-CFD row wins regardless of input order, even when its
        aligned DNA is lexicographically LARGER (so a naive DNA-only sort would have
        picked the wrong row)."""
        dna_hi = "TTTTTCCGAGCAGAAGAAGAAGG"   # higher CFD, larger DNA
        dna_lo = "AAAATCCGAGCAGAAGAAGAAGG"   # lower  CFD, smaller DNA
        self.assertGreater(dna_hi, dna_lo)
        row_hi = _row(dna_hi, cfd="0.950", total=4, cluster_pos=1000, position=1001)
        row_lo = _row(dna_lo, cfd="0.500", total=4, cluster_pos=1000, position=1002)

        rep_ab = self._representative([row_hi, row_lo])
        rep_ba = self._representative([row_lo, row_hi])

        self.assertEqual(rep_ab[CFD_COL], "0.950")
        self.assertEqual(rep_ba[CFD_COL], "0.950")
        self.assertEqual(rep_ab[DNA_COL], dna_hi,
                         "best CFD must win over the smaller-DNA tie-break")
        self.assertEqual(rep_ba[DNA_COL], dna_hi,
                         "best CFD must win regardless of input order")

    def test_no_samples_sentinel_counts_as_zero_carriers(self):
        """A "NO_SAMPLES" roster (the empty-roster sentinel used across PostProcess,
        e.g. process_summaries.py) must count as ZERO carriers in the 2.7.10
        most-carriers tie-break -- not as one carrier from a naive comma split. The
        no-carrier row is given the lexicographically-SMALLER aligned DNA, so if
        "NO_SAMPLES" were miscounted as 1 carrier the two rows would tie on carriers and
        the smaller-DNA (no-carrier) row would win. With the sentinel correctly at 0,
        the genuine single-carrier row wins despite its larger DNA, in both input orders."""
        dna_nosamp = "AAAATCCGAGCAGAAGAAGAAGG"   # no carriers, SMALLER DNA
        dna_carrier = "TTTTTCCGAGCAGAAGAAGAAGG"  # one real carrier, LARGER DNA
        self.assertLess(dna_nosamp, dna_carrier)
        row_nosamp = _row(dna_nosamp, cfd="0.900", total=4, cluster_pos=1000,
                          position=1001, samples="NO_SAMPLES")
        row_carrier = _row(dna_carrier, cfd="0.900", total=4, cluster_pos=1000,
                           position=1002, samples="HG00001")

        rep_ab = self._representative([row_nosamp, row_carrier])
        rep_ba = self._representative([row_carrier, row_nosamp])

        self.assertEqual(rep_ab[DNA_COL], dna_carrier,
                         "the genuine 1-carrier row must win; NO_SAMPLES must count as 0")
        self.assertEqual(rep_ba[DNA_COL], dna_carrier,
                         "NO_SAMPLES-as-0 must hold regardless of input order")

    # 2.8.1: a row of SNPs and a row with an indel with the same alignment.

    DNA_TIE = "CTAACAGTTGCTTTTATCACTGG"

    def _variant_row(self, snp, samples, dna=None, cfd="0.500", position=1001):
        return _row(dna or self.DNA_TIE, cfd=cfd, total=4, cluster_pos=1000,
                    position=position, snp=snp, samples=samples)

    def test_snp_row_is_the_representative_at_a_tie_with_an_indel_row(self):
        """Same aligned DNA, position, CFD, total and number of carriers (other
        carriers): the row of SNPs in both input orders. Before this element the
        representative was the first row of the input, so one of the two orders failed
        for every form below."""
        forms = [
            ("deletion", "chr1_1010_G_C", "chr1_1015_ACT_A"),
            ("insertion", "chr1_1010_G_C", "chr1_1015_A_ACT"),
            ("SNP that sorts after the indel", "chr1_1090_G_C", "chr1_1015_ACT_A"),
            ("two SNPs, an indel with a SNP", "chr1_1010_G_C,chr1_1012_A_T",
             "chr1_1010_G_C,chr1_1015_ACT_A"),
            ("contig name with '_'", "chr1_KI270706v1_random_1010_G_C",
             "chr1_KI270706v1_random_1015_ACT_A"),
        ]
        for name, snps, indel in forms:
            with self.subTest(form=name):
                self._assert_representative_both_orders(
                    self._variant_row(snps, "HG00002,HG00003"),
                    self._variant_row(indel, "HG00001,HG00004"),
                    name)

    def test_has_indel_is_only_the_last_word(self):
        """A row with an indel still wins on everything that comes before: more carriers,
        a smaller aligned DNA, a smaller Position, a higher CFD."""
        snp, indel = "chr1_1010_G_C", "chr1_1015_ACT_A"
        controls = [
            ("more carriers",
             self._variant_row(indel, "HG00001,HG00002"),
             self._variant_row(snp, "HG00003")),
            ("smaller aligned DNA",
             self._variant_row(indel, "HG00001", dna="ATAACAGTTGCTTTTATCACTGG"),
             self._variant_row(snp, "HG00003", dna="GTAACAGTTGCTTTTATCACTGG")),
            ("smaller Position",
             self._variant_row(indel, "HG00001", position=1001),
             self._variant_row(snp, "HG00003", position=1002)),
        ]
        for name, winner, loser in controls:
            with self.subTest(control=name):
                self._assert_representative_both_orders(winner, loser, name)
        # a higher CFD decides on the score branch only: the total branch does not read it
        winner = self._variant_row(indel, "HG00001", cfd="0.900")
        loser = self._variant_row(snp, "HG00003", cfd="0.500")
        for rows in ((winner, loser), (loser, winner)):
            rep = self._representative([list(r) for r in rows], "score")
            self.assertEqual(rep[SNP_COL], indel)

    def test_reference_row_is_not_a_variant_row(self):
        """The marker of a reference row ("n") is not a variant token, so it is not an
        indel; a reference row and a variant row go through the same key without error."""
        ref_row = _row(self.DNA_TIE, cfd="0.700", total=4, cluster_pos=1000, position=1001)
        var_row = self._variant_row("chr1_1015_ACT_A", "HG00001", cfd="0.500")
        for rows in ((ref_row, var_row), (var_row, ref_row)):
            rep = self._representative([list(r) for r in rows], "score")
            self.assertEqual(rep[SNP_COL], "n")


if __name__ == "__main__":
    unittest.main()

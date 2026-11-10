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
IDs decide. The has-indel element settles it among rows with carriers: the row of SNPs
is the representative.

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
    higher CFD still win over it; rows without carriers (a sites-only index) keep the
    order of the merge step's sort.
  * Alignments that start inside one insertion (2.10.0: all at its anchor): two targets
    of the insertion with one carrier set both stay, the site's row and an alternative
    alignment, on both strands; two copies of one target are one row; a bulged alignment
    of one target, before the same PAM, still folds into the best (as up to 2.8.1).
  * Two alignments of one indel target string (2.10.0: the merge keys an indel row by its
    aligned spacer too): a carrier resolved PUTATIVE stays on the alignment it was resolved
    on (the SNP+indel oracle's site U, S19); one alignment from two VCF datasets is still
    one entry; alignments of one target with one carrier set are one entry, the best; carrier
    sets are compared whole (an alignment whose carriers are a subset of a better one's keeps
    its entry), the same entries in every order.
  * Cost (2.10.0): one target's alignments inside a 54 kb insertion, each its own aligned
    spacer, have their place on the haplotype read once per target, not once per alignment
    (indel_snp_cis.row_first_base, counted in a fresh interpreter).

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
_COMP = str.maketrans("ACGTacgt-", "TGCAtgca-")


def _rc(seq):
    return seq.translate(_COMP)[::-1]


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

    def test_rows_without_carriers_keep_their_order(self):
        """The rows of a sites-only index (mega) have no carriers, and their REF alignment
        was never upper case: their ties stay in the order of the merge step's sort, so
        these searches give what they gave before the SNP-first rule. The two rows are
        told apart here by the order alone, whichever comes first is reported."""
        snp, indel = "chr1_1010_G_C", "chr1_1015_ACT_A"
        for samples in ("", "NO_SAMPLES"):
            first_snp = self._variant_row(snp, samples)
            first_indel = self._variant_row(indel, samples)
            for rows in ((first_snp, first_indel), (first_indel, first_snp)):
                for sort_order in ("score", "total"):
                    with self.subTest(samples=samples, first=rows[0][SNP_COL],
                                      sort_order=sort_order):
                        rep = self._representative([list(r) for r in rows], sort_order)
                        self.assertEqual(rep[SNP_COL], rows[0][SNP_COL])

    def test_reference_row_is_not_a_variant_row(self):
        """The marker of a reference row ("n") is not a variant token, so it is not an
        indel; a reference row and a variant row go through the same key without error."""
        ref_row = _row(self.DNA_TIE, cfd="0.700", total=4, cluster_pos=1000, position=1001)
        var_row = self._variant_row("chr1_1015_ACT_A", "HG00001", cfd="0.500")
        for rows in ((ref_row, var_row), (var_row, ref_row)):
            rep = self._representative([list(r) for r in rows], "score")
            self.assertEqual(rep[SNP_COL], "n")


class IndelRowsKeepTheirOwnCarriers(unittest.TestCase):
    """The merge step and SNP+indel rows (2.8.1).

    Rows of one cluster naming the same variants at one Cluster_Position merge into the
    first (their samples unioned), whatever their sequence. On an indel row the carriers
    belong to the sequence (each carrier haplotype is resolved on its own), and two
    alignments of a site can name the same variants with different carriers: an
    RNA-bulge alignment does not reach the base a SNP sits on, so it lists that SNP's
    carriers under the indel alone; the union then put them on the full alignment, a
    sequence they do not have."""

    INDEL = "chr1_1010_TC_T"
    X = "ACTGAAATCTGTAAGCAGGCTGG"
    RNA = "-CTGAAATCTGTAAGCAGGCTGG"

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="merge_indel_rows_")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.d, ignore_errors=True)

    def _row(self, dna, cfd, total, samples, snp=None, rsid="rs1", af="0.1"):
        r = _row(dna, cfd=cfd, total=total, samples=samples, snp=snp or self.INDEL)
        r[3] = dna.upper()  # the Reference column holds the reference alignment
        r[16], r[17] = rsid, af
        return r

    def _merge(self, rows, sort_order="score", criteria="mm+bulges"):
        infile = os.path.join(self.d, "in.txt")
        outfile = os.path.join(self.d, "out.txt")
        with open(infile, "w") as fh:
            fh.write("\t".join(POST_ADJUST_HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        cmd = [sys.executable, os.path.join(PP, "remove_contiguous_samples.py"),
               infile, outfile, "3", "5", "7", "11", "16", "19", "21",
               sort_order, criteria]
        res = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        out = []
        for path in (outfile, outfile + ".discarded_samples"):
            with open(path) as fh:
                out.append([l.rstrip("\n").split("\t") for l in fh.readlines()[1:]])
        return out

    def test_a_carrier_stays_off_a_sequence_it_does_not_have(self):
        best, alt = self._merge([self._row(self.X, "1.000", 0, "S23"),
                                 self._row(self.RNA, "0.900", 1, "S23,S24")])
        self.assertEqual((best[0][DNA_COL], best[0][13]), (self.X, "S23"))
        self.assertEqual([(a[DNA_COL], a[13]) for a in alt], [(self.RNA, "S23,S24")])
        self.assertEqual(best[0][19], "1")  # the cluster's other rows

    def test_alignments_with_one_carrier_set_stay_one_entry(self):
        """The rule as before: one row per variant list at a site, the best kept."""
        best, alt = self._merge([self._row(self.X, "1.000", 0, "S23,S24"),
                                 self._row(self.RNA, "0.900", 1, "S24,S23")])
        self.assertEqual((best[0][DNA_COL], best[0][13]), (self.X, "S23,S24"))
        self.assertEqual(alt, [])

    def test_one_sequence_from_two_datasets_is_still_merged(self):
        """What the merge is for: the same target from two VCF datasets."""
        best, alt = self._merge([
            self._row(self.X, "1.000", 0, "S23", rsid="rs1", af="0.1"),
            self._row(self.X, "1.000", 0, "T01", rsid=".", af="0.2")])
        self.assertEqual(best[0][13], "S23,T01")
        # each rsID beside its AF ("." is T01's, 0.2): (rsID, AF) pairs, by variant
        self.assertEqual(best[0][16], ".,rs1")
        self.assertEqual(best[0][17], "0.2,0.1")
        self.assertEqual(alt, [])

    def test_folded_alignments_still_add_their_rsid_and_af(self):
        """Another sequence with the carriers the first has: one entry, as before."""
        best, alt = self._merge([
            self._row(self.X, "1.000", 0, "S23", rsid="rs1", af="0.1"),
            self._row(self.RNA, "0.900", 1, "S23", rsid=".", af="0.2")])
        self.assertEqual((best[0][DNA_COL], best[0][16], best[0][17]),
                         (self.X, ".,rs1", "0.2,0.1"))
        self.assertEqual(alt, [])

    def test_snp_only_rows_merge_as_before(self):
        """A row with no indel keeps the old key, whatever its sequence."""
        snp = "chr1_1006_C_A"
        best, alt = self._merge([self._row(self.X, "1.000", 0, "S23", snp=snp),
                                 self._row(self.RNA, "0.900", 1, "S23,S24", snp=snp)])
        self.assertEqual((best[0][DNA_COL], best[0][13]), (self.X, "S23,S24"))
        self.assertEqual(alt, [])

    def test_a_snp_named_after_the_indel_counts_as_an_indel_row(self):
        snps = self.INDEL + ",chr1_1006_C_A"
        best, alt = self._merge([self._row(self.X, "1.000", 0, "S23", snp=snps),
                                 self._row(self.RNA, "0.900", 1, "S24", snp=snps)])
        self.assertEqual(best[0][13], "S23")
        self.assertEqual([a[13] for a in alt], ["S24"])

    def test_an_mnv_row_keeps_its_own_carriers(self):
        """A multi-base record of equal length (AC>TG) goes to a fake-indel window like
        an indel, and its rows are resolved per carrier: S24's slot holds a SNP the full
        alignment reaches and the RNA-bulge alignment does not."""
        for mnv in ("chr1_1010_AC_TG", "chr1_1006_C_A,chr1_1010_AC_TG",
                    "chrUn_KI270742v1_1010_AC_TG"):
            best, alt = self._merge([self._row(self.X, "1.000", 0, "S23", snp=mnv),
                                     self._row(self.RNA, "0.900", 1, "S23,S24", snp=mnv)])
            self.assertEqual((best[0][DNA_COL], best[0][13]), (self.X, "S23"), mnv)
            self.assertEqual([(a[DNA_COL], a[13]) for a in alt], [(self.RNA, "S23,S24")], mnv)

    def test_snp_rows_on_a_contig_with_underscores_merge_as_before(self):
        snp = "chrUn_KI270742v1_1006_C_A"
        best, alt = self._merge([self._row(self.X, "1.000", 0, "S23", snp=snp),
                                 self._row(self.RNA, "0.900", 1, "S23,S24", snp=snp)])
        self.assertEqual((best[0][DNA_COL], best[0][13]), (self.X, "S23,S24"))
        self.assertEqual(alt, [])

    def _tied(self, bulge_type, dna, mm, bulges, rsid="rs1", af="0.1"):
        r = self._row(dna, "0.001", mm + bulges, "S1,S2,S3", rsid=rsid, af=af)
        r[0], r[8], r[9] = bulge_type, str(mm), str(bulges)
        return r

    def test_the_folded_entry_is_the_best_by_the_criteria(self):
        """Two alignments of one site with the same variants and carriers, tied on CFD
        and on mismatches + bulges. The input is in the merge step's file order, where
        'RNA<tab>' sorts before 'RNA,DNA'. The entry kept is the one the file's sort
        ranks first, in the score files and in the mismatches + bulges file alike, and
        it takes the other's rsID and AF."""
        rna = self._tied("RNA", "tCTGAAActTGTAAtCAGG-AAG", 4, 1, rsid=".", af="0.2")
        rna_dna = self._tied("RNA,DNA", "ATCTGAAActTGTAAtCAGG-AAG", 3, 2)
        for criteria, want in (("mm+bulges,mm", rna_dna), ("mm+bulges,bulges", rna),
                               ("mm", rna_dna), ("bulges", rna)):
            for sort_order in ("score", "total"):
                best, alt = self._merge([list(rna), list(rna_dna)], sort_order, criteria)
                self.assertEqual((best[0][0], best[0][DNA_COL], best[0][8], best[0][9]),
                                 (want[0], want[DNA_COL], want[8], want[9]),
                                 (criteria, sort_order))
                self.assertEqual((best[0][13], best[0][16], best[0][17]),
                                 ("S1,S2,S3", ".,rs1", "0.2,0.1"), (criteria, sort_order))
                self.assertEqual(alt, [], (criteria, sort_order))
                self.assertEqual(best[0][19], "0")

    def test_the_row_count_does_not_depend_on_the_sort(self):
        """The three files of a search hold the same rows, each sorted by its own
        criterion, and their alternative alignments are joined line by line (pr -m), so
        a cluster must keep one row count in all three."""
        rows = [self._row(self.X, "1.000", 0, "S1"),
                self._row(self.RNA, "0.900", 1, "S1,S2"),
                self._row("ACTGAAATCTGTAAGCAGG-CAG", "0.300", 1, "S1"),
                self._row("A-TGAAATCTGTAAGCAGGCTGG", "0.200", 1, "S1,S2"),
                self._row("ACTGAA-TCTGTAAGCAGGCTGG", "0.100", 1, "S2")]
        counts = set()
        for order in (rows, rows[::-1], rows[2:] + rows[:2], rows[1::2] + rows[::2]):
            best, alt = self._merge([list(r) for r in order])
            self.assertEqual(len(best), 1)
            counts.add(len(alt))
        self.assertEqual(counts, {2})  # the carrier sets {S1}, {S1,S2}, {S2}


class OneTargetTwoAlignments(unittest.TestCase):
    """2.10.0: an indel row merges only with rows of its own alignment, its target and its
    aligned spacer the same.

    A carrier resolved PUTATIVE (phase unknown) is listed with the strongest target its own
    alleles allow on each alignment (indel_snp_cis.worst_case_bases): the guide's base at a
    spacer column, the CFD's choice in a DNA bulge. So two alignments of one target string
    can list different carriers. The rows are the SNP+indel oracle's site U as the indel
    post-analysis writes them: the deletion chr21:8010 TC>T, the SNP chr21:8006 C>A whose A
    the guide needs, S19 0/1 for both, S20 0/1 for the deletion alone. The target
    AC-TGA?ATCTGTAAGCAGGCTGG (? the SNP's base) has two alignments: A-CTGAAATCTG..., where ?
    faces the guide's A (S19 with A, CFD 0.93; S20 with C, 0.423), and ACTGAA-ATCTG..., where
    ? is a DNA bulge and the CFD takes C (r-:dG,7 0.059, r-:dT,7 0), so S19 and S20 with C
    (0.014). Up to 8d3f278 the merge keyed an indel row by its target alone, and the two C
    rows were one entry, the first in the file's order: the CFD 0.423 alignment, which then
    listed S19 with C, where its alleles give that alignment A (0.93)."""

    DEL, SNP = "chr21_8010_TC_T", "chr21_8006_C_A"
    G = "ACTGAAATCTGTAAGCAGGCNNN"
    BULGED_A = "A-CTGAAATCTGTAAGCAGGCNNN"  # the SNP's base faces the guide's A
    BULGED_B = "ACTGAA-ATCTGTAAGCAGGCNNN"  # the SNP's base is a DNA bulge
    REF = "ACTGACATCTCGTAAGCAGGCTGG"  # the reference from 8000 (the deleted C at 8010)

    setUp = IndelRowsKeepTheirOwnCarriers.setUp
    tearDown = IndelRowsKeepTheirOwnCarriers.tearDown
    _merge = IndelRowsKeepTheirOwnCarriers._merge

    def _u(self, crrna, dna, cfd, mm, samples, snp=False, rsid=".", af="0.0128205"):
        bulges = crrna.count("-") + dna.count("-")
        r = _row(dna, cfd=cfd, total=mm + bulges, cluster_pos=8000, position=8000,
                 snp=(self.SNP + "," + self.DEL) if snp else self.DEL, chrom="chr21",
                 guide=crrna, samples=samples)
        ref = self.REF[:len(dna.replace("-", ""))]
        for i, c in enumerate(dna):
            if c == "-":
                ref = ref[:i] + "-" + ref[i:]
        r[0] = {0: "X", 1: "DNA" if "-" in crrna else "RNA"}.get(bulges, "RNA,DNA")
        r[3], r[8], r[9] = ref, str(mm), str(bulges)
        r[16], r[17] = (".,rs1" if snp else rsid), ("0.00641026," + af if snp else af)
        return r

    def _site_u(self):
        return [self._u(self.G, "ACTGAAATCTGTAAGCAGGCTGG", "1.000", 0, "S19", snp=True),
                self._u(self.BULGED_A, "AC-TGAAATCTGTAAGCAGGCTGG", "0.930", 0, "S19", snp=True),
                self._u(self.BULGED_A, "AC-TGAcATCTGTAAGCAGGCTGG", "0.423", 1, "S20"),
                self._u(self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG", "0.014", 2, "S19,S20")]

    def _all(self, rows, sort_order):
        best, alt = self._merge([list(r) for r in rows], sort_order)
        return best + alt

    def test_a_putative_carrier_stays_on_the_alignment_it_was_resolved_on(self):
        """In the merge step's order of each file (CFD then mismatches + bulges, or
        mismatches + bulges), and reversed: S19 is listed under the deletion alone only on
        the alignment that makes the SNP's base a DNA bulge, and the other alignment of that
        target lists S20 alone. One entry per carrier set and variant list: S19 under both
        variants (the perfect alignment, A-CTGAAATCTG...'s A folded into it, as up to 2.8.1),
        S20 and S19 + S20 under the deletion alone; the same three whatever the sort, so the
        three files keep one row count."""
        rows = self._site_u()
        by_cfd = sorted(rows, key=lambda r: (-float(r[CFD_COL]), int(r[TOTAL_COL])))
        by_total = sorted(rows, key=lambda r: int(r[TOTAL_COL]))
        for sort_order, order in (("score", by_cfd), ("score", by_cfd[::-1]),
                                  ("total", by_total), ("total", by_total[::-1])):
            with self.subTest(sort_order=sort_order, first=order[0][1]):
                out = self._all(order, sort_order)
                got = {(r[1], r[DNA_COL]): (r[SNP_COL], r[13]) for r in out}
                self.assertEqual(got, {
                    (self.G, "ACTGAAATCTGTAAGCAGGCTGG"): (self.SNP + "," + self.DEL, "S19"),
                    (self.BULGED_A, "AC-TGAcATCTGTAAGCAGGCTGG"): (self.DEL, "S20"),
                    (self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG"): (self.DEL, "S19,S20")})

    def test_one_alignment_from_two_datasets_is_still_merged(self):
        """The same alignment from two VCF datasets (a carrier T01 of the second): one entry
        with the carriers, rsIDs and AFs of both."""
        one = self._u(self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG", "0.014", 2, "S19,S20")
        two = self._u(self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG", "0.014", 2, "T01",
                      rsid="rs7", af="0.2")
        for sort_order in ("score", "total"):
            out = self._all([one, two], sort_order)
            self.assertEqual(len(out), 1)
            self.assertEqual((out[0][13], out[0][16], out[0][17]),
                             ("S19,S20,T01", ".,rs7", "0.0128205,0.2"))

    def test_alignments_of_one_target_with_one_carrier_set_are_one_entry(self):
        """Two alignments of one target string (ACTG-AcATCTG...) listing one carrier set
        (S20, its slot's C facing the guide's A in both) are one entry, the best by the order
        that picks each block's representative (the CFD on the score branch, whichever comes
        first), with the other's rsID and AF; and the alignment of that target where the C
        is a DNA bulge, which lists S19 too, stays an entry of its own."""
        strong = self._u("ACT-GAAATCTGTAAGCAGGCNNN", "ACTG-AcATCTGTAAGCAGGCTGG", "0.042", 1,
                         "S20", rsid="rs9", af="0.5")
        weak = self._u("ACTGA-AATCTGTAAGCAGGCNNN", "ACTG-AcATCTGTAAGCAGGCTGG", "0.012", 1,
                       "S20")
        bulge = self._u(self.BULGED_B, "ACTG-ACATCTGTAAGCAGGCTGG", "0.019", 0, "S19,S20")
        for sort_order in ("score", "total"):
            for order in ((strong, weak, bulge), (weak, strong, bulge), (bulge, weak, strong)):
                with self.subTest(sort_order=sort_order, first=order[0][1]):
                    out = self._all(order, sort_order)
                    self.assertEqual(len(out), 2)
                    s20 = [r for r in out if r[13] == "S20"]
                    self.assertEqual(len(s20), 1)
                    self.assertEqual((s20[0][16], s20[0][17]), (".,rs9", "0.0128205,0.5"))
                    if sort_order == "score":
                        self.assertEqual(s20[0][1], strong[1])
                    self.assertEqual([r[1] for r in out if r[13] == "S19,S20"], [bulge[1]])

    def test_carrier_sets_are_compared_whole(self):
        """The fold compares whole carrier sets, so that every block keeps one row count
        whatever its order. S04 (homozygous for the SNP) and S06 (phased, the SNP in cis)
        carry the A wherever the target is; the unphased S01 and S03 take it where the
        alignment pairs it with the guide's A (the perfect alignment) and the C where it is
        a DNA bulge (ACTGAA-ATCTG..., whose CFD the A makes 0). So ACTGAA-ATCTG... with the
        A lists S04 and S06 alone and keeps an entry of its own, though both are listed on
        the perfect alignment of that target and variants (METHODS: an extra alternative
        alignment, never a wrong carrier or score). The same three entries in every order."""
        rows = [self._u(self.G, "ACTGAAATCTGTAAGCAGGCTGG", "1.000", 0, "S01,S03,S04,S06",
                        snp=True),
                self._u(self.BULGED_B, "AC-tgAAATCTGTAAGCAGGCTGG", "0.000", 2, "S04,S06",
                        snp=True),
                self._u(self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG", "0.014", 2, "S01,S03")]
        both = self.SNP + "," + self.DEL
        for sort_order in ("score", "total"):
            for order in (rows, rows[::-1], rows[1:] + rows[:1]):
                with self.subTest(sort_order=sort_order, first=order[0][DNA_COL]):
                    out = self._all(order, sort_order)
                    got = {(r[1], r[DNA_COL]): (r[SNP_COL], r[13]) for r in out}
                    self.assertEqual(got, {
                        (self.G, "ACTGAAATCTGTAAGCAGGCTGG"): (both, "S01,S03,S04,S06"),
                        (self.BULGED_B, "AC-tgAAATCTGTAAGCAGGCTGG"): (both, "S04,S06"),
                        (self.BULGED_B, "AC-tgACATCTGTAAGCAGGCTGG"): (self.DEL, "S01,S03")})


class AlignmentsStartingInsideOneInsertion(unittest.TestCase):
    """2.10.0: a row whose first base is one an insertion adds starts at the insertion's
    anchor (indel_snp_cis.reference_start), so alignments that start on different bases
    of one insertion have one Position and one Cluster_Position, and are one site.

    The merge folds a site's alignments naming the same variants with the same carriers
    into one entry, and compares them where they are on the haplotype carrying the indel
    (remove_contiguous_samples._haplotype_cluster_position): there, two targets that start
    on different inserted bases are as many bases apart as 2.8.1 placed them. So they stay,
    the site's row and its alternative alignments (on the Cluster_Position alone they fold
    into one: on the HGDP + 1000 Genomes index, an alignment of the TRAC guide's left the
    tables), and one target's alignments with different bulges before one PAM still fold,
    as up to 2.8.1. Two copies of one target inside an insertion (a repeat it holds) are
    the same row at the anchor: one entry.

    The rows are those the indel post-analysis writes: the insertion after the anchor at
    0-based 100 (VCF POS 101), the target's bases on the haplotype carrying it, the REF
    text the reference's from the anchor, as many bases as the target holds."""

    T1 = "ACTGAAATCTGTAAGCAGGCTGG"  # the guide + TGG
    T2 = "ACTGAAATCTGTAAGCAGGCAGG"  # the guide + AGG
    REFSEQ = "GATTACAGCT" * 10 + "G" + "CCTAGGTCAA" * 10  # the anchor at 100: G

    setUp = IndelRowsKeepTheirOwnCarriers.setUp
    tearDown = IndelRowsKeepTheirOwnCarriers.tearDown
    _row = IndelRowsKeepTheirOwnCarriers._row
    _merge = IndelRowsKeepTheirOwnCarriers._merge
    INDEL = IndelRowsKeepTheirOwnCarriers.INDEL

    def _site(self, ins, anchor=100):
        """The insertion's token, and the haplotype carrying it (the insertion after the
        0-based ``anchor``)."""
        base = self.REFSEQ[anchor]
        return ("chr1_%d_%s_%s" % (anchor + 1, base, base + ins),
                self.REFSEQ[:anchor + 1] + ins + self.REFSEQ[anchor + 1:])

    def _inside(self, ins, k, n, strand="+", cfd="1.000", total=0, dna_bulges=0,
                samples="S1", rsid="rs1", af="0.1", rna_bulges=0, anchor=100):
        """A row of the alignment whose first base (the leftmost) is the k-th base the
        insertion adds, its target the n bases from there; a DNA bulge one base of them,
        an RNA bulge one base of the guide facing none (a gap in the target)."""
        token, hap = self._site(ins, anchor)
        bases, ref = hap[anchor + k:anchor + k + n], self.REFSEQ[anchor:anchor + n]
        if strand == "-":
            bases, ref = _rc(bases), _rc(ref)
        if rna_bulges:
            bases = bases[:5] + "-" + bases[5:]
        r = self._row(bases, cfd, total, samples, snp=token, rsid=rsid, af=af)
        r[3], r[7] = ref, strand
        r[5] = str(anchor)  # the anchor
        # the Cluster_Position: the start, plus the DNA bulges and less the RNA bulges on
        # the + strand, 3' PAM
        r[CLUSTER_POS_COL] = str(anchor + ((dna_bulges - rna_bulges) if strand == "+" else 0))
        if dna_bulges:
            r[0], r[1], r[9] = "DNA", "-" + GUIDE, str(dna_bulges)
        if rna_bulges:
            r[0], r[9] = "RNA", str(rna_bulges)
        return r

    def test_other_targets_of_one_insertion_stay_as_alternatives(self):
        """Two targets inside one insertion, 25 bases apart on its haplotype (2.8.1: two
        sites), the same variant and carrier: one site at the anchor, the second target
        one of its alternative alignments, on both strands and both branches."""
        for strand, ins in (("+", "AT" + self.T1 + "CA" + self.T2 + "TC"),
                            ("-", "AT" + _rc(self.T1) + "CA" + _rc(self.T2) + "TC")):
            rows = [self._inside(ins, 3, 23, strand), self._inside(ins, 28, 23, strand)]
            want = {r[DNA_COL] for r in rows}
            for sort_order in ("score", "total"):
                for order in (rows, rows[::-1]):
                    with self.subTest(strand=strand, sort_order=sort_order):
                        best, alt = self._merge([list(r) for r in order], sort_order)
                        self.assertEqual((len(best), len(alt)), (1, 1))
                        self.assertEqual({best[0][DNA_COL], alt[0][DNA_COL]}, want)
                        self.assertEqual(best[0][19], "1")  # the site's other row
                        self.assertEqual((best[0][13], alt[0][13]), ("S1", "S1"))

    def test_two_copies_of_one_target_are_one_row(self):
        """The insertion holds the target twice (a repeat): the two alignments are the
        same row at the anchor, its bases, REF text, variant and carriers, so one entry
        (2.8.1 reported each copy at its own start)."""
        ins = "AT" + self.T1 + "CA" + self.T1 + "TC"
        rows = [self._inside(ins, 3, 23), self._inside(ins, 28, 23)]
        self.assertEqual(rows[0], rows[1])
        for sort_order in ("score", "total"):
            best, alt = self._merge([list(r) for r in rows], sort_order)
            self.assertEqual((len(best), len(alt)), (1, 0))
            self.assertEqual((best[0][DNA_COL], best[0][13]), (self.T1, "S1"))

    def test_two_copies_are_created_when_either_copys_pam_is(self):
        """Two copies of one target inside one insertion are one row (above), and the indel
        post-analysis judges each copy's PAM on its own: the copy that a placement of the
        insertion moves onto the reference's own PAM is not created ("n"), the copy the
        insertion adds can be (its ALT PAM). The row is created when either copy is,
        whichever the merge meets first (2.10.0; up to 50803b6 the first one met
        was kept, and only the pre-merge sort's byte order, an upper-case PAM before "n",
        put the created copy first: TRAC, HGDP + 1000 Genomes, chr22:50193928 (-))."""
        ins = "AT" + self.T1 + "CA" + self.T1 + "TC"
        rows = [self._inside(ins, 3, 23), self._inside(ins, 28, 23)]
        rows[0][11], rows[1][11] = "n", "TGG"
        for sort_order in ("score", "total"):
            for order in (rows, rows[::-1]):
                with self.subTest(sort_order=sort_order, first=order[0][11]):
                    best, alt = self._merge([list(r) for r in order], sort_order)
                    self.assertEqual((len(best), len(alt)), (1, 0))
                    self.assertEqual((best[0][DNA_COL], best[0][11]), (self.T1, "TGG"))
        # neither created: "n"
        rows[1][11] = "n"
        best, alt = self._merge([list(r) for r in rows])
        self.assertEqual(best[0][11], "n")

    def test_bulged_alignments_of_one_target_still_fold(self):
        """A DNA-bulge alignment of the same target, one base longer and before the same
        PAM, starts one inserted base earlier: on the haplotype the same site (2.8.1: one
        Cluster_Position), here Cluster_Positions 100 and 101. One entry, the better one,
        with the other's rsID and AF, each pair in its place (2.10.0: the rsID and AF lists
        follow the variants, so the "." of the bulged row goes with its 0.2)."""
        ins = "AT" + self.T1 + "CA"
        perfect = self._inside(ins, 3, 23)
        bulged = self._inside(ins, 2, 24, cfd="0.500", total=1, dna_bulges=1, rsid=".",
                              af="0.2")
        self.assertEqual(bulged[DNA_COL], "T" + self.T1)
        for sort_order in ("score", "total"):
            for order in ((perfect, bulged), (bulged, perfect)):
                best, alt = self._merge([list(r) for r in order], sort_order)
                self.assertEqual((best[0][DNA_COL], best[0][16], best[0][17]),
                                 (self.T1, ".,rs1", "0.2,0.1"))
                self.assertEqual(alt, [])

    def test_a_target_on_the_anchor_and_one_inside_stay_apart(self):
        """The target of a row starting on the anchor base itself (its first base a
        reference one) and one starting inside the insertion, the same carrier: two
        alignments, as up to 2.8.1."""
        ins = "AT" + self.T2 + "CA"
        token, hap = self._site(ins)
        on_anchor = self._row(hap[100:123], "0.300", 3, "S1", snp=token)
        on_anchor[3], on_anchor[5] = self.REFSEQ[100:123], "100"
        on_anchor[CLUSTER_POS_COL] = "100"
        rows = [self._inside(ins, 3, 23), on_anchor]
        best, alt = self._merge([list(r) for r in rows])
        self.assertEqual((best[0][DNA_COL], [a[DNA_COL] for a in alt]),
                         (self.T2, [hap[100:123]]))


    def test_an_insertion_row_at_its_anchor_chains_no_two_sites(self):
        """A row starting inside an insertion stands at the anchor, k bases before where
        2.8.1 placed it, and chained there it joined two sites 2.8.1 kept apart: on
        hg38_mega chr22, population-level, chr22:36,838,374 + of VEGFA site 3 (a deletion
        row, 2 bases before such a row's anchor) left the results for the cluster of the
        reference row 3 bases after it, its place on the haplotype. It joins one site --
        the one nearer its place on the haplotype -- and chains none: two clusters, two
        rows, on both branches. Without the row starting inside, the same two."""
        ins = "AT" + self.T2 + "CA"
        inside = self._inside(ins, 3, 23)            # Cluster_Position 100, haplotype 103
        before = _row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.200", 3, cluster_pos=98,
                      position=98)
        after = _row(self.T1, "0.400", 2, cluster_pos=103, position=103)
        for rows in ([before, inside, after], [before, after]):
            for sort_order in ("score", "total"):
                with self.subTest(n=len(rows), sort_order=sort_order):
                    best, alt = self._merge([list(r) for r in rows], sort_order)
                    # the site before the anchor keeps its row
                    self.assertEqual(len(best), 2)
                    self.assertIn(before[DNA_COL], [r[DNA_COL] for r in best])
                    if len(rows) == 3:
                        # the inside row is in the cluster of the site after its anchor,
                        # whose row it is (the better), that site's alignment its
                        # alternative
                        self.assertIn(inside[DNA_COL], [r[DNA_COL] for r in best])
                        self.assertEqual([r[DNA_COL] for r in alt], [after[DNA_COL]])
                    else:
                        self.assertEqual(alt, [])

    def test_an_insertion_row_with_no_site_near_is_a_site_of_its_own(self):
        """The rows that start inside an insertion with no other row within --merge bp
        are a cluster of their own, the insertion's targets one site, and they chain no
        site to another insertion's: here a second insertion 3 bases on (its anchor at
        103) has a site 3 bases past its anchor (106), and joins it; the first one's
        rows, 6 bases from that site, stay a site of their own. Chained like any other
        row (up to the merge of this release), the first insertion's rows reached the
        site through the second's: one cluster, one row."""
        ins = "AT" + self.T1 + "CA" + self.T2 + "TC"
        rows = [self._inside(ins, 3, 23), self._inside(ins, 28, 23, cfd="0.900"),
                self._inside("GA" + self.T2 + "TT", 3, 23, cfd="0.300", total=2, anchor=103),
                _row(self.T1, "0.400", 2, cluster_pos=106, position=106)]
        for sort_order in ("score", "total"):
            with self.subTest(sort_order=sort_order):
                best, alt = self._merge([list(r) for r in rows], sort_order)
                self.assertEqual(sorted(r[CLUSTER_POS_COL] for r in best), ["100", "106"])
                self.assertEqual(sorted(r[CLUSTER_POS_COL] for r in alt), ["100", "103"])
                # the first insertion's site: its two targets, the best the row
                self.assertIn(self.T1, [r[DNA_COL] for r in best])
                self.assertEqual([r[19] for r in best if r[CLUSTER_POS_COL] == "100"], ["1"])

    def test_one_targets_bulged_alignments_between_two_sites_stay_one_row(self):
        """One target inside an insertion (on the haplotype at 103, 3 bases past the
        anchor) aligned without a bulge (Cluster_Position 100) and with a DNA bulge
        (Cluster_Position 101), between a site 3 bases before the anchor (97) and one 3
        bases after the bulged alignment (104): the rows of one insertion are placed
        together, with the site nearer their place on the haplotype, and the bulged
        alignment folds into the other. Placed each by its own Cluster_Position, the plain
        one joined the site before and the bulged one the site after: two rows of one
        off-target (chr22:45,890,432 (+) of VEGFA site 3 on HGDP + 1000 Genomes, CFD 0.179
        and 0.172, one insertion and one carrier; review r1 of lane M, M-R1-1)."""
        ins = "AT" + self.T2 + "CA"
        plain = self._inside(ins, 3, 23, cfd="0.400", total=2)
        bulged = self._inside(ins, 2, 24, cfd="0.300", total=3, dna_bulges=1)
        before = _row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.100", 4, cluster_pos=97,
                      position=97)
        after = _row(self.T1, "0.100", 4, cluster_pos=104, position=104)
        self.assertEqual((plain[CLUSTER_POS_COL], bulged[CLUSTER_POS_COL]), ("100", "101"))
        for sort_order in ("score", "total"):
            with self.subTest(sort_order=sort_order):
                best, alt = self._merge([list(r) for r in (before, plain, bulged, after)],
                                        sort_order)
                self.assertEqual(sorted(r[DNA_COL] for r in best),
                                 sorted([before[DNA_COL], plain[DNA_COL]]))
                self.assertEqual([r[DNA_COL] for r in alt], [after[DNA_COL]])

    def test_an_rna_bulge_alignment_of_one_target_stays_with_the_plain_one(self):
        """One target inside an insertion aligned without a bulge (Cluster_Position 100)
        and with an RNA bulge (99: one base fewer, from the next inserted base), a site 3
        bases before the RNA-bulge alignment (96) and one 4 bases after the plain one
        (104): the insertion's rows are within --merge bp of the first site and join it,
        both; the RNA-bulge alignment folds into the plain one, and the site after keeps
        its row. Placed each by its own Cluster_Position, the RNA-bulge alignment joined
        the site before and the plain one, 4 bases from either, was a site of its own: two
        rows of one off-target."""
        ins = "AT" + self.T2 + "CA"
        plain = self._inside(ins, 3, 23, cfd="0.400", total=2)
        rna = self._inside(ins, 4, 22, cfd="0.300", total=3, rna_bulges=1)
        before = _row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.100", 4, cluster_pos=96,
                      position=96)
        after = _row(self.T1, "0.200", 4, cluster_pos=104, position=104)
        self.assertEqual((rna[CLUSTER_POS_COL], plain[CLUSTER_POS_COL]), ("99", "100"))
        for sort_order in ("score", "total"):
            with self.subTest(sort_order=sort_order):
                best, alt = self._merge([list(r) for r in (before, rna, plain, after)],
                                        sort_order)
                self.assertEqual(sorted(r[DNA_COL] for r in best),
                                 sorted([plain[DNA_COL], after[DNA_COL]]))
                self.assertEqual([r[DNA_COL] for r in alt], [before[DNA_COL]])

    def test_the_clusters_do_not_depend_on_the_order_of_rows_of_one_place(self):
        """(A guard of what the merge kept, not a change.) The three merges sort the rows
        of one Cluster_Position differently (the CFD, the edits, the CRISPR-Bulge score),
        and their alternative alignments are pasted line by line, so they must make the
        same clusters: rows of one Cluster_Position in any order, the same rows written."""
        ins = "AT" + self.T2 + "CA"
        unit = [self._inside(ins, 3, 23, cfd="0.400", total=2),
                self._inside(ins, 2, 24, cfd="0.300", total=3, dna_bulges=1),
                self._inside(ins, 4, 22, cfd="0.200", total=3, rna_bulges=1)]
        chain = [_row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.100", 4, cluster_pos=97,
                      position=97),
                 _row(self.T1, "0.100", 4, cluster_pos=101, position=101),
                 _row(self.T1[:-1] + "A", "0.050", 4, cluster_pos=101, position=101),
                 _row(self.T1, "0.150", 4, cluster_pos=104, position=104)]
        rows = sorted(unit + chain, key=lambda r: int(r[CLUSTER_POS_COL]))
        seen = set()
        for flip in (False, True):
            ordered = sorted(rows, key=lambda r: (int(r[CLUSTER_POS_COL]),
                                                  r[DNA_COL] if not flip else
                                                  "".join(chr(255 - ord(c)) for c in r[DNA_COL])))
            best, alt = self._merge([list(r) for r in ordered])
            seen.add((tuple(sorted(r[DNA_COL] for r in best)),
                      tuple(sorted(r[DNA_COL] for r in alt))))
        self.assertEqual(len(seen), 1, seen)

    def _on_reference(self, ins, start, n, cfd, total, dna_bulges=0, samples="S1"):
        """A row naming the insertion whose first base is the reference base ``start`` (the
        anchor, 100, or one before it), its target the n bases of the haplotype from there:
        a row of the chain like any other, its Cluster_Position the start plus its DNA
        bulges (+ strand, 3' PAM)."""
        token, hap = self._site(ins)
        r = self._row(hap[start:start + n], cfd, total, samples, snp=token)
        r[3], r[5] = self.REFSEQ[start:start + n], str(start)
        r[CLUSTER_POS_COL] = str(start + dna_bulges)
        if dna_bulges:
            r[0], r[1], r[9] = "DNA", "-" * dna_bulges + GUIDE, str(dna_bulges)
        return r

    def _target_written(self, rows, target):
        """The rows of ``target`` (its aligned DNA) the best file holds, on each branch."""
        out = {}
        for sort_order in ("score", "total"):
            best, _alt = self._merge([list(r) for r in rows], sort_order)
            out[sort_order] = sorted(r[DNA_COL] for r in best if r[DNA_COL] in target)
        return out

    def test_a_targets_alignment_from_the_anchor_base_goes_with_its_alignments_inside(self):
        """2 DNA + 2 RNA bulges (review r2 of lane M, M-R2-1). Target T, its first base the
        insertion's 2nd (on the haplotype at 102), aligned plain (Cluster_Position 100) and
        with a DNA bulge from the 1st inserted base (101), both inside the insertion, and
        with 2 DNA bulges from the anchor base itself: a row of the chain like any other,
        at 102. Another target's 2-RNA-bulge alignment from the 1st inserted base (98, on
        the haplotype at 99) and a reference row at 98. The rows inside the insertion were
        one unit, which went to the site before, nearer its first place (99), while T's
        alignment from the anchor base was a site of its own: T written twice, CFD 0.500
        and 0.300 (the chain of every row, and 2.8.1, wrote it once). A target with an
        alignment on a reference base is placed by it: its alignments inside go with it,
        and T is one row, the best of its alignments."""
        ins = "AT" + self.T2 + "CA"
        plain = self._inside(ins, 2, 23, cfd="0.500", total=2)
        dna1 = self._inside(ins, 1, 24, cfd="0.400", total=3, dna_bulges=1)
        dna2 = self._on_reference(ins, 100, 25, "0.300", 4, dna_bulges=2)
        other = self._inside(ins, 1, 21, cfd="0.100", total=4, rna_bulges=2)
        ref = _row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.200", 4, cluster_pos=98, position=98)
        self.assertEqual([r[CLUSTER_POS_COL] for r in (ref, other, plain, dna1, dna2)],
                         ["98", "98", "100", "101", "102"])
        target = {plain[DNA_COL], dna1[DNA_COL], dna2[DNA_COL]}
        self.assertEqual(len(target), 3)
        rows = [ref, other, plain, dna1, dna2]
        self.assertEqual(self._target_written(rows, target),
                         {"score": [plain[DNA_COL]], "total": [plain[DNA_COL]]})
        # the reference row keeps its site, the other target's alignment its alternative
        best, alt = self._merge([list(r) for r in rows])
        self.assertEqual(sorted(r[DNA_COL] for r in best), sorted([ref[DNA_COL], plain[DNA_COL]]))
        self.assertEqual([r[DNA_COL] for r in alt], [other[DNA_COL]])

    def test_a_targets_alignment_from_a_base_before_the_anchor_goes_with_them_too(self):
        """The same with 3 DNA bulges: T's alignment from the base before the anchor (99,
        Cluster_Position 102) is the one on a reference base, and the alignments inside go
        with it all the same (a rule for the anchor base alone would leave it a site of its
        own): T is one row."""
        ins = "AT" + self.T2 + "CA"
        plain = self._inside(ins, 2, 23, cfd="0.500", total=2)
        dna1 = self._inside(ins, 1, 24, cfd="0.400", total=3, dna_bulges=1)
        dna3 = self._on_reference(ins, 99, 26, "0.300", 5, dna_bulges=3)
        other = self._inside(ins, 1, 21, cfd="0.100", total=4, rna_bulges=2)
        ref = _row("GGACTGAAATCTGTAAGCAGGCAGG"[:23], "0.200", 4, cluster_pos=98, position=98)
        self.assertEqual(dna3[CLUSTER_POS_COL], "102")
        target = {plain[DNA_COL], dna1[DNA_COL], dna3[DNA_COL]}
        self.assertEqual(self._target_written([ref, other, plain, dna1, dna3], target),
                         {"score": [plain[DNA_COL]], "total": [plain[DNA_COL]]})


# Runs the real merge script in a fresh interpreter with indel_snp_cis.row_first_base
# counting its calls, and prints their number.
_COUNTING_MERGE = """
import runpy, sys
import indel_snp_cis as isc
calls = []
_read = isc.row_first_base
def _counted(*args, **kwargs):
    calls.append(args[1])
    return _read(*args, **kwargs)
isc.row_first_base = _counted
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
sys.stderr.write("ROW_FIRST_BASE_CALLS %d\\n" % len(calls))
"""


class OneTargetsAlignmentsInsideALongInsertion(unittest.TestCase):
    """2.10.0, cost: the merge reads a row's place on the haplotype carrying its indel
    (_haplotype_cluster_position) once per target, REF text, place and variants, not once
    per aligned spacer, which the merge also keys an indel row by.

    The rows start inside an insertion of 54 kb, as at chr22:10690626 of HPRC, where 402
    alleles of 53.9 to 54.4 kb hold 16,058 rows of the README guide's chr22 search: 20
    DNA-bulge alignments of one target (the gap in a different place of the guide, the
    mismatches in a different case), its perfect alignment and a second target. Read per
    aligned spacer, 22 times; per target, 3. The fold is the same either way: the bulged
    alignments fold into the perfect one, and the second target is its alternative."""

    T1 = AlignmentsStartingInsideOneInsertion.T1
    T2 = AlignmentsStartingInsideOneInsertion.T2
    INS = "C" * 30000 + "AT" + T1 + "CA" + T2 + "TC" + "G" * 23950
    setUp = IndelRowsKeepTheirOwnCarriers.setUp
    tearDown = IndelRowsKeepTheirOwnCarriers.tearDown
    _site = AlignmentsStartingInsideOneInsertion._site
    _inside = AlignmentsStartingInsideOneInsertion._inside
    _row = IndelRowsKeepTheirOwnCarriers._row
    REFSEQ = AlignmentsStartingInsideOneInsertion.REFSEQ
    INDEL = IndelRowsKeepTheirOwnCarriers.INDEL

    def _merge_counting(self, rows):
        infile = os.path.join(self.d, "in.txt")
        outfile = os.path.join(self.d, "out.txt")
        with open(infile, "w") as fh:
            fh.write("\t".join(POST_ADJUST_HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        cmd = [sys.executable, "-c", _COUNTING_MERGE,
               os.path.join(PP, "remove_contiguous_samples.py"), infile, outfile,
               "3", "5", "7", "11", "16", "19", "21", "score", "mm+bulges"]
        res = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        calls = [int(l.split()[1]) for l in res.stderr.splitlines()
                 if l.startswith("ROW_FIRST_BASE_CALLS ")]
        out = []
        for path in (outfile, outfile + ".discarded_samples"):
            with open(path) as fh:
                out.append([l.rstrip("\n").split("\t") for l in fh.readlines()[1:]])
        return calls, out

    def test_the_place_is_read_once_per_target(self):
        k1 = 30000 + 3  # T1's first base: the k-th base the insertion adds
        perfect = self._inside(self.INS, k1, 23)
        bulged = []
        for i in range(20):
            r = self._inside(self.INS, k1 - 1, 24, cfd="0.%03d" % (100 + i), total=1,
                             dna_bulges=1, rsid=".", af="0.2")
            r[1] = GUIDE[:i] + "-" + GUIDE[i:]
            r[DNA_COL] = r[DNA_COL][:i + 1] + r[DNA_COL][i + 1].lower() + r[DNA_COL][i + 2:]
            bulged.append(r)
        self.assertEqual(len({(r[1], r[DNA_COL]) for r in bulged}), 20)
        self.assertEqual(len({r[DNA_COL].upper() for r in bulged}), 1)
        second = self._inside(self.INS, k1 + 25, 23, cfd="0.800", total=1)
        calls, (best, alt) = self._merge_counting([perfect] + bulged + [second])
        self.assertEqual(calls, [3])
        self.assertEqual((len(best), len(alt)), (1, 1))
        # (the rsID and AF of each folded row together, as the lists follow the variants)
        self.assertEqual((best[0][DNA_COL], best[0][16], best[0][17]), (self.T1, ".,rs1", "0.2,0.1"))
        self.assertEqual(alt[0][DNA_COL], self.T2)



def _merge_module():
    """remove_contiguous_samples.py's definitions, in this process, with the columns
    merge_close_targets_cfd.sh passes (0-based); its run as a script is cut off."""
    import types
    path = os.path.join(PP, "remove_contiguous_samples.py")
    with open(path) as fh:
        src = fh.read()
    cut = src.find('\nif __name__ == "__main__":')
    if cut < 0:
        cut = src.index("\ntau = int(sys.argv[3])")
    if PP not in sys.path:
        sys.path.insert(0, PP)
    mod = types.ModuleType("remove_contiguous_samples_defs")
    mod.__file__ = path
    exec(compile(src[:cut], path, "exec"), mod.__dict__)
    mod.chrom, mod.pos, mod.total, mod.true_guide, mod.snp_info, mod.cfd = 4, 6, 10, 15, 18, 20
    return mod


class TheClustersHoldWholeTargets(unittest.TestCase):
    """(2.10.0, review r2 of lane M.) The merge's clusters (_clusters) of rows drawn from a
    model of targets: a random reference, random insertions of 1 to 12 bases after an
    anchor (now and then two at one anchor), and targets on them, each a PAM's placement on
    a haplotype (the reference, or one carrying one insertion). Each alignment of a target
    with nd DNA and nr RNA bulges within the budget (some left out at random) is a row as
    the indel post-analysis writes it: its bases, REF text, Position and Cluster_Position,
    so the merge reads each row's place on the haplotype itself (_inside_place,
    row_first_base). A search has one PAM side: with a 3' PAM a + strand alignment's first
    base is s - nd + nr and its Cluster_Position its start plus nd less nr, a - strand
    one's first base is s; with a 5' PAM the other way round. A row whose first base is an
    inserted one starts at the anchor; one past the insertion holds none of its bases and
    is a reference row (another target, on the reference).

    Over 1+1 to 3+3 bulges and --merge 0 to 6: every target's rows are in one cluster (the
    first merge of this release wrote a target twice at 2 DNA + 2 RNA bulges and more,
    M-R2-1); the clusters do not depend on the order of the rows of one Cluster_Position;
    with no indel row they are the plain chain of every row; restricted to the rows of
    the targets that have a row on a reference base they are the plain chain of their
    targets' places (a target with a row on a reference base is placed by it, and the
    rows of the targets wholly inside an insertion join one site and chain none); and the
    rows of the targets wholly inside the insertions at one anchor are in one cluster
    (decision (d))."""

    N = 23  # a target's bases on its haplotype aligned with no bulge (20 + the PAM)
    BUDGETS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 3), (0, 2), (2, 0), (3, 1), (1, 3), (0, 0)]
    TAUS = [0, 1, 2, 3, 3, 3, 4, 5, 6]

    @classmethod
    def setUpClass(cls):
        cls.mod = _merge_module()

    def _trial(self, rng):
        """(rows, {row id: (target, its place, 'in' | 'chain' | 'ref', anchor)}, budget)."""
        bd, br = rng.choice(self.BUDGETS)
        kinds = ["L+", "L+", "F-"] if rng.random() < 0.7 else ["F+", "L-", "L-"]
        ref = "".join(rng.choice("ACGT") for _ in range(420))
        rows, model = [], {}

        def target(ins, s):
            kind = rng.choice(kinds)
            hap = ins[3] if ins else ref
            combos = [(nd, nr) for nd in range(bd + 1) for nr in range(br + 1)]
            for nd, nr in [c for c in combos if rng.random() < 0.6] or [rng.choice(combos)]:
                if kind[0] == "L":  # the left end moves with the bulges
                    f, e, adj = s - nd + nr, s + self.N, nd - nr
                else:
                    f, e, adj = s, s + self.N + nd - nr, 0
                if f < 0 or e > len(hap):
                    continue
                bases, m = hap[f:e], e - f
                if ins is not None and f <= ins[0] + ins[1] and e - 1 >= ins[0] + 1:
                    a, n_ins, token, _h = ins
                    start = a if f > a else f
                    reftext, snp = ref[start:start + m], token
                    tag, place = "T|%s|%s|%d" % (token, kind, s), s
                    where = "in" if f > a else "chain"
                else:
                    start = f - ins[1] if ins is not None and f > ins[0] + ins[1] else f
                    reftext, snp = "n", "n"
                    tag, place, where = "R|%s|%d" % (kind, start + adj), start + adj, "ref"
                strand = kind[1]
                dna = bases if strand == "+" else _rc(bases)
                if reftext != "n" and strand == "-":
                    reftext = _rc(reftext)
                for _ in range(nr):  # an RNA bulge: a gap in the target
                    i = rng.randint(1, len(dna) - 1)
                    dna = dna[:i] + "-" + dna[i:]
                r = ["X"] * 22
                r[2], r[3], r[4], r[5] = dna, reftext, "chr1", str(start)
                r[6], r[7], r[10], r[13] = str(start + adj), strand, str(nd + nr), str(len(rows))
                r[15], r[18], r[20] = "GUIDE", snp, "0.5"
                rows.append(r)
                model[r[13]] = (tag, place, where, ins[0] if ins else None)

        x = 40
        while x < 330:
            x += rng.randint(1, 9)
            u = rng.random()
            if u < 0.6:
                for allele in range(1 if rng.random() < 0.85 else 2):
                    n_ins = rng.randint(1, 12)
                    ins = "".join(rng.choice("ACGT") for _ in range(n_ins))
                    token = "chr1_%d_%s_%s" % (x + 1, ref[x], ref[x] + ins)
                    site = (x, n_ins, token, ref[:x + 1] + ins + ref[x + 1:])
                    for _t in range(rng.randint(1, 4) if allele == 0 else 1):
                        target(site, x + rng.randint(-4, n_ins + 3))
            elif u < 0.85:
                target(None, x)
        return rows, model, (bd, br)

    def _clusters(self, rows, tau):
        mod = self.mod
        mod.tau = tau
        mod._FIRST_BASE_MEMO.clear()
        shift = max(abs(int(r[6]) - int(r[5])) for r in rows)
        return list(mod._clusters(["\t".join(r) for r in rows], shift))

    @staticmethod
    def _chain(items, tau):
        out, cur, last = [], [], None
        for place, rid in sorted(items):
            if cur and place - last > tau:
                out.append(sorted(cur))
                cur = []
            cur.append(rid)
            last = place
        if cur:
            out.append(sorted(cur))
        return sorted(out)

    def test_the_clusters_of_rows_of_random_targets(self):
        import random
        rng = random.Random(20261110)
        seen = {"targets split": 0, "inside rows read elsewhere": 0}
        trials = 0
        while trials < 600:
            rows, model, budget = self._trial(rng)
            if not rows:
                continue
            memo = {}
            if any(self.mod._haplotype_cluster_position(r, memo) != model[r[13]][1]
                   for r in rows if model[r[13]][2] != "ref"):
                # a repeat the model does not follow: the row's bases are at an earlier place
                seen["inside rows read elsewhere"] += 1
                continue
            trials += 1
            tau = rng.choice(self.TAUS)
            got = self._clusters(sorted(rows, key=lambda r: (int(r[6]), rng.random())), tau)
            part = sorted(sorted(r[13] for r in c) for c in got)
            cluster_of = {r[13]: i for i, c in enumerate(got) for r in c}
            with self.subTest(trial=trials, budget=budget, tau=tau):
                by_target = {}
                for rid, (tag, _p, _w, _a) in model.items():
                    by_target.setdefault(tag, set()).add(cluster_of[rid])
                split = sorted(t for t, cs in by_target.items() if len(cs) > 1)
                seen["targets split"] += bool(split)
                self.assertEqual(split, [], "a target's rows in two clusters")
                again = self._clusters(sorted(rows, key=lambda r: (int(r[6]), rng.random())), tau)
                self.assertEqual(sorted(sorted(r[13] for r in c) for c in again), part,
                                 "the clusters follow the order of rows of one place")
                if all(w == "ref" for _t, _p, w, _a in model.values()):
                    self.assertEqual(part, self._chain([(int(r[6]), r[13]) for r in rows], tau))
                kinds = {}
                for tag, _p, w, _a in model.values():
                    kinds.setdefault(tag, set()).add(w)
                wholly_inside = {t for t, k in kinds.items() if k == {"in"}}
                placed = [(p, rid) for rid, (t, p, _w, _a) in model.items()
                          if t not in wholly_inside]
                ids = {rid for _p, rid in placed}
                self.assertEqual(sorted(c for c in (sorted(i for i in c if i in ids)
                                                    for c in part) if c),
                                 self._chain(placed, tau),
                                 "an insertion's unit chained two sites, or a target with a "
                                 "row on a reference base was not placed by it")
                units = {}
                for rid, (t, _p, _w, a) in model.items():
                    if t in wholly_inside:
                        units.setdefault(a, set()).add(cluster_of[rid])
                self.assertEqual({a: cs for a, cs in units.items() if len(cs) > 1}, {},
                                 "the targets wholly inside one anchor's insertions in two sites")
        # the trials are not idle: most hold rows inside an insertion and bulged alignments
        self.assertLess(seen["inside rows read elsewhere"], 60)


if __name__ == "__main__":
    unittest.main()

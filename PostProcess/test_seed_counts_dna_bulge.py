"""Seed / non-seed counts of the fewest_mm+b block include the DNA bulge (2.8.1).

BUG (pre-existing, every search with --bDNA >= 1, dictionary and tier path alike)
---------------------------------------------------------------------------------
``resultIntegrator.py`` fills Seed_mismatches+bulges_{REF,ALT}_(<criterion>) and
Non_seed_mismatches+bulges_{REF,ALT}_(<criterion>) from the aligned strings:
``seed_processing`` counts the lower case and '-' of the protospacer, which leaves out a
DNA bulge (a '-' in the aligned spacer, not in the protospacer), and each block counts
those columns apart in ``dna_bulge_count_seed`` / ``dna_bulge_count_non_seed``. The
highest_CFD and highest_CRISPR_BULGE blocks add them. The fewest_mm+b block computed
them and wrote ``seed_list`` alone, so on a row with a DNA bulge its seed + non-seed was
the row's mismatches + bulges minus the DNA bulge, and the same alignment read
differently in the fewest_mm+b block and in the other two.

WHAT THIS TEST ASSERTS
----------------------
Drives the REAL ``resultIntegrator.py`` through a subprocess on hand-built rows of the
bestMerge layout (three 24-column blocks holding the same alignment) -- DNA bulge in the
seed half, in the non-seed half, together with an RNA bulge, on a reference row, with a
5' PAM, and an RNA-bulge row and a no-bulge row as controls:
  * each of the four fewest_mm+b columns holds the count expected by hand;
  * the fewest_mm+b block equals the highest_CFD and the highest_CRISPR_BULGE blocks;
  * in all three blocks seed + non-seed adds up to Mismatches+bulges (of the ALT
    alignment on a row of alt origin, of the REF alignment on a row of ref origin);
  * a row of ref origin still reads 0 / 0 in the ALT columns.

No annotation, empirical or VCF inputs are read (``vuota`` skips the VCF directory, the
annotation file holds "NA" per row). ``intervaltree`` is not among the light CI deps and
its trees stay empty here, so a stand-in goes on PYTHONPATH.

Run with:
    cd PostProcess && python3 -m unittest test_seed_counts_dna_bulge -v
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))

CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
SPACER_3P = "CTAACAGTTGCTTTTATCACNNN"
SPACER_5P = "NNNNACGTTGCAGTCAGTACGTCA"

INTERVALTREE_STUB = '''class IntervalTree:
    """Stand-in for the test: no empirical data, no high-variant-density bed."""
    def __setitem__(self, key, value):
        pass
    def __getitem__(self, key):
        return set()
'''


def _block(bulge_type, spacer, aligned_guide, dna, reference, mismatches, bulges,
           samples, snp):
    """One 24-column block of a bestMerge row (column names as in its header)."""
    return [
        bulge_type, aligned_guide, dna, reference, "chrT", "1000", "1000", "+",
        str(mismatches), str(bulges), str(mismatches + bulges),
        "NA", "NA", samples, "CDS_gencode", spacer,
        "rs1" if snp != "NA" else "NA", "0.25" if snp != "NA" else "NA", snp, "0",
        "0.5", "0.4", "0.1", "0.1",
    ]


# name: spacer, bulge type, aligned guide, ALT alignment, REF alignment (None on a row of
# ref origin), mismatches, bulges, then the expected
# (REF seed, REF non-seed, ALT seed, ALT non-seed).
# 3' PAM: the protospacer has 21 columns here, the seed is columns 10 to 20; the DNA-bulge
# column keeps its base in upper case and is not counted by seed_processing.
CASES = {
    "dna_bulge_in_seed_half": (
        SPACER_3P, "DNA", "CTAACAGTTGCTTTT-ATCACNNN",
        "CTAACtGTTGCcTTTAATCACTGG",
        "CTAACtGTTGCcTTTAATCgCTGG",
        2, 1, (3, 1, 2, 1)),
    "dna_bulge_in_non_seed_half": (
        SPACER_3P, "DNA", "CTAAC-AGTTGCTTTTATCACNNN",
        "CTAACAAGTTGCTTcTATCACTGG",
        "CTaACAAGTTGCTTcTATCACTGG",
        1, 1, (1, 2, 1, 1)),
    "rna_and_dna_bulge": (
        SPACER_3P, "RNA,DNA", "CT-AACAGTTGCTTTTATCACNNN",
        "-TTAAC-GTTGCcTTTATCACTGG",
        "-TTAAC-GTTGCcTTTATCgCTGG",
        1, 3, (2, 3, 1, 3)),
    "dna_bulge_on_reference_row": (
        SPACER_3P, "DNA", "CTAACAGTTGCTTTT-ATCACNNN",
        "CTAACtGTTGCcTTTAATCACTGG",
        None,
        2, 1, (2, 1, 0, 0)),
    "five_prime_pam_dna_bulge_in_seed_half": (
        SPACER_5P, "DNA", "NNNNACG-TTGCAGTCAGTACGTCA",
        "TTTAACGATTGCAGTCAGTAgGTCA",
        "TTTAACGATTGCAGTCAGTAgGcCA",
        1, 1, (1, 2, 1, 1)),
    "five_prime_pam_dna_bulge_in_non_seed_half": (
        SPACER_5P, "DNA", "NNNNACGTTGCAGTCAGT-ACGTCA",
        "TTTAACcTTGCAGTCAGTAACGTCA",
        "TTTAtCcTTGCAGTCAGTAACGTCA",
        1, 1, (2, 1, 1, 1)),
    "control_rna_bulge_only": (
        SPACER_3P, "RNA", "CTAACAGTTGCTTTTATCACNNN",
        "CT-ACtGTTGCTTTTATCACTGG",
        "CT-ACtGTTGCTTTTATCgCTGG",
        1, 1, (1, 2, 0, 2)),
    "control_no_bulge": (
        SPACER_3P, "X", "CTAACAGTTGCTTTTATCACNNN",
        "CTAACtGTTGCTTTTATCACTGG",
        "CTAACtGTTGCTTTTATCgCTGG",
        1, 0, (1, 1, 0, 1)),
}


def _row(case):
    spacer, bulge_type, guide, dna, ref, mismatches, bulges, _ = CASES[case]
    if ref is None:  # a reference row, after remove_n_and_dots: Samples, SNP etc. are NA
        block = _block(bulge_type, spacer, guide, dna, "NA", mismatches, bulges, "NA", "NA")
    else:
        block = _block(bulge_type, spacer, guide, dna, ref, mismatches, bulges,
                       "HG00001", "chrT_1005_A_T")
    return block + block + block


def _gnu_cut():
    """resultIntegrator drops its unused columns with GNU ``cut --complement``."""
    try:
        res = subprocess.run(["cut", "--complement", "-f1"], input=b"a\tb\n",
                             capture_output=True)
    except OSError:
        return False
    return res.returncode == 0


@unittest.skipUnless(_gnu_cut(), "resultIntegrator.py needs GNU cut --complement")
class SeedCountsIncludeTheDnaBulge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="seed_counts_")
        cls.addClassCleanup(shutil.rmtree, cls.d, ignore_errors=True)
        with open(os.path.join(cls.d, "intervaltree.py"), "w") as fh:
            fh.write(INTERVALTREE_STUB)
        names = sorted(CASES)
        best = os.path.join(cls.d, "job.bestMerge.txt")
        with open(best, "w") as fh:
            for name in names:
                fh.write("\t".join(_row(name)) + "\n")
        found = os.path.join(cls.d, "job.found.bed")
        with open(found, "w") as fh:
            fh.write("NA\n" * len(names))
        empirical = os.path.join(cls.d, "dummy.txt")
        open(empirical, "w").close()
        guides = os.path.join(cls.d, "guides.txt")
        with open(guides, "w") as fh:
            fh.write(SPACER_3P + "\n")
        env = dict(os.environ, PYTHONPATH=cls.d)
        res = subprocess.run(
            [sys.executable, os.path.join(PP, "resultIntegrator.py"), best, empirical,
             found, guides, cls.d + "/", "true", "hg38", "vuota", "1", "0", "none"],
            cwd=PP, env=env, capture_output=True, text=True)
        if res.returncode or res.stderr:
            raise AssertionError("resultIntegrator.py failed (stderr is fatal in the "
                                 "pipeline): rc=%s\n%s" % (res.returncode, res.stderr))
        with open(best + ".integrated_results.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        cls.rows = dict(zip(names, rows))

    @staticmethod
    def _counts(row, criterion):
        return (int(row["Seed_mismatches+bulges_REF_(%s)" % criterion]),
                int(row["Non_seed_mismatches+bulges_REF_(%s)" % criterion]),
                int(row["Seed_mismatches+bulges_ALT_(%s)" % criterion]),
                int(row["Non_seed_mismatches+bulges_ALT_(%s)" % criterion]))

    def test_every_case_was_integrated(self):
        self.assertEqual(sorted(self.rows), sorted(CASES))

    def test_fewest_block_counts(self):
        for name, case in CASES.items():
            with self.subTest(case=name):
                self.assertEqual(self._counts(self.rows[name], "fewest_mm+b"), case[-1])

    def test_all_three_blocks_agree(self):
        """Same alignment in the three blocks, same counts."""
        for name, case in CASES.items():
            for criterion in CRITERIA:
                with self.subTest(case=name, criterion=criterion):
                    self.assertEqual(self._counts(self.rows[name], criterion), case[-1])

    def test_seed_plus_non_seed_is_mismatches_plus_bulges(self):
        """On the row's own alignment: ALT for a row of alt origin, REF for ref origin."""
        for name in CASES:
            row = self.rows[name]
            for criterion in CRITERIA:
                with self.subTest(case=name, criterion=criterion):
                    ref_seed, ref_non, alt_seed, alt_non = self._counts(row, criterion)
                    origin = row["REF/ALT_origin_(%s)" % criterion]
                    own = alt_seed + alt_non if origin == "alt" else ref_seed + ref_non
                    self.assertEqual(own, int(row["Mismatches+bulges_(%s)" % criterion]))

    def test_reference_row_alt_columns_stay_zero(self):
        row = self.rows["dna_bulge_on_reference_row"]
        for criterion in CRITERIA:
            self.assertEqual(row["REF/ALT_origin_(%s)" % criterion], "ref")
            self.assertEqual(self._counts(row, criterion)[2:], (0, 0))


if __name__ == "__main__":
    unittest.main()

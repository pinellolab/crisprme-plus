"""Tests for the position-aware CFD bulge penalty (bulge_cfd_penalty.py)."""
import unittest

import bulge_cfd_penalty as bp
import cfd_score

# Clean-room CFD tables (committed JSON), loaded once.
_MM, _PAM = cfd_score.cfd_tables()


def _c(cr, pr):
    return bp.bulge_cfd_correction(cr, pr, _MM)


def _s(cr, pr):
    return bp.cfd_score_bulge_aware(cr, pr, _MM, _PAM)


class TestBulgeCfdPenalty(unittest.TestCase):
    def test_mismatch_only_unchanged(self):
        # no gaps -> correction 1.0 (CFD mismatch scoring stays byte-identical / CRISPOR-faithful)
        self.assertEqual(_c("ACTGAAATCTGTAAGCAGGCNGG", "ACTGAAATCTGTAAGCAGGCAGG"), 1.0)

    def test_distal_dna_bulge_uses_doench_and_stays_lenient(self):
        # real DNA bulge at PAM-distal position 2 -> Doench r-:dA,2 (lenient, correct)
        c = _c("A-CTGAAATCTGTAAGCAGGCNNN", "ATCTGAAATCTGTAAGCAGGCGGG")
        self.assertAlmostEqual(c, _MM["r-:dA,2"], places=6)
        self.assertGreater(c, 0.5)  # distal bulges remain tolerated

    def test_seed_dna_bulge_heavily_penalized(self):
        # DNA bulge in the seed (position 15) -> near-zero Doench factor
        c = _c("ACTGAAATCTGTAA-GCAGGCNNN", "ACTGAAATCTGTAATGCAGGCGGG")
        self.assertLess(c, 0.1)

    def test_pam_proximal_bulge_abolished(self):
        # chr19:42345264 pattern: PAM-proximal DNA bulge (col 21) -> 0 (CFD 0.872 -> ~0),
        # matching CRISPR-Bulge = 0. The RNA bulge at col 2 is already handled by the
        # clean-room base scorer (in-range) so it contributes no extra correction here.
        c = _c("ACTGAAATCTGTAAGCAGGC-NNN", "t-TaAAATCTaTAAGaAGGCCAGG")
        self.assertEqual(c, 0.0)

    def test_in_range_rna_bulge_not_double_counted(self):
        # RNA bulge at an in-range position (gap in protospacer, col 5): the clean-room base
        # scorer already applied its rX:d-,5 factor, so the correction must be 1.0 (no double penalty).
        c = _c("ACTGAAATCTGTAAGCAGGCNGG", "ACTG-AATCTGTAAGCAGGCAGG")
        self.assertEqual(c, 1.0)

    def test_combined_dna_then_rna_bulge_guide_frame_no_double_count(self):
        # Regression for the frame-inconsistency bug on combined (DNA+RNA) bulge rows that
        # crispritz DOES emit under bDNA>=1,bRNA>=1 (e.g. the real chr19/chr13 6,1,1 rows).
        # DNA bulge at display col2 + RNA bulge at display col21: deleting the DNA-bulge
        # column shifts the RNA bulge to its TRUE guide position 20 (in-range), which the
        # clean-room base scorer already scores. The correction must NOT also treat it as
        # terminal (0.0); it must contribute only the distal DNA-bulge factor.
        corr = _c("A-CTGAAATCTGTAAGCAGGCNGG", "AACTGAAATCTGTAAGCAGG-AGG")
        self.assertGreaterEqual(corr, 0.5)   # a distal DNA-bulge factor, NOT 0.0
        self.assertNotEqual(corr, 0.0)

    def test_combined_chr19_dna_bulge_at_pam_junction_still_zero(self):
        # chr19 has the DNA bulge at display col21 with NO preceding crRNA gap, so its
        # guide-frame position is a genuine 21 (PAM junction) -> abolished. The guide-frame
        # fix must NOT dilute this to position 20.
        self.assertEqual(_c("ACTGAAATCTGTAAGCAGGC-NNN", "t-TaAAATCTaTAAGaAGGCCAGG"), 0.0)

    def test_terminal_rna_bulge_distal_lenient_proximal_abolished(self):
        # RNA bulge free-passed by the clean-room base scorer at the terminal columns:
        # distal (1) lenient, proximal/PAM (21) abolished.
        self.assertEqual(bp._terminal_factor(1), 0.9)
        self.assertEqual(bp._terminal_factor(21), 0.0)
        self.assertIsNone(bp._terminal_factor(10))  # in-range -> use the matrix


class TestCfdScoreBulgeAware(unittest.TestCase):
    """End-to-end: cfd_score_bulge_aware on real display alignments (ml007 6,1,1 rows)."""

    def test_mismatch_only_byte_identical_to_clean_room_cfd(self):
        cr = "ACTGAAATCTGTAAGCAGGCNGG"
        pr = "ACTGAAATCTGTAAGCAGGCAGG"
        self.assertEqual(_s(cr, pr), cfd_score.cfd(cr, pr[:-3], pr[-2:], (_MM, _PAM)))

    def test_chr19_pam_proximal_dna_bulge_goes_to_zero(self):
        # real chr19:42345264 combined row: DNA bulge at PAM-proximal col21 + RNA bulge col2.
        # Reported CFD under the old free-pass was 0.872; the position-aware score must be ~0.
        s = _s("ACTGAAATCTGTAAGCAGGC-NNN", "t-TaAAATCTaTAAGaAGGCCAGG")
        self.assertEqual(s, 0.0)

    def test_distal_dna_bulge_stays_lenient(self):
        # real chr6:122734951 pure DNA bulge at distal col2: old 0.528 -> *0.96 (r-:dG,2).
        s = _s("A-CTGAAATCTGTAAGCAGGCNNN", "ACCTGgAAgCTGTAAGaAGGCTGG")
        self.assertAlmostEqual(s, 0.528 * _MM["r-:dG,2"], places=3)
        self.assertGreater(s, 0.4)  # distal bulge remains actionable

    def test_terminal_rna_bulge_pos1_scaled_0p9(self):
        # real pos-1 RNA bulge: in-range mismatches unchanged, whole score * 0.9.
        cr = "ACTGAAATCTGTAAGCAGGCNNN"
        pr = "-CTGAAAcaaaTAAaCAGGCGGG"
        base = cfd_score.cfd(cr, pr.upper()[:-3], pr.upper()[-2:], (_MM, _PAM))
        self.assertAlmostEqual(_s(cr, pr), base * 0.9, places=6)

    def test_untranslatable_iupac_short_circuits(self):
        # the clean-room base scorer returns 0 for an untranslatable base; the wrapper
        # must not multiply it.
        self.assertEqual(_s("ACTGAAATCTGTAAGCAGGCNGG", "ANTGAAATCTGTAAGCAGGCAGG"), 0)


if __name__ == "__main__":
    unittest.main()

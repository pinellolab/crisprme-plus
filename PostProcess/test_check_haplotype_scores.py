"""``scripts/check_haplotype_scores.py`` -- the standing gate.

A checker nobody checks is a rubber stamp. Every invariant it claims to enforce is
exercised here BOTH ways: a file that violates it must fail, and a correct file must
pass. The second half matters as much as the first -- on its first run against real data
this script failed a CORRECT file, because it compared rows to haplotypes where the
companion's counters are per haplotype. A false alarm in a release gate is not a safe
error: it teaches whoever sees it to ignore the gate.
"""

import io
import os
import shutil
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.normpath(os.path.join(_HERE, "..", "scripts"))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import check_haplotype_scores as chk  # noqa: E402


HEADER_COMMENTS = [
    "# CRISPRme+ per-haplotype scores.",
    "# THRESHOLD. above_threshold is the CRISPR-Bulge sensitive operating point,",
    "#   0.00019 at 95% recall, taken BY NAME from cb_operating_point.py.",
    "#     scored                                   -- at least one axis produced a number",
    "#     out-of-domain:model-refused(n/length/multi-gap) -- reached and declined",
    "#     needs-realignment:<site>                 -- the frame cannot hold it",
    "#     yes                                      -- measured, at or above",
    "#     no                                       -- measured, below",
    "#     unknown                                  -- not measured",
    "#     unknown(model-refused)                   -- reached and declined",
    "#     n/a                                      -- the table says ref",
]

COLS = ["site_key", "variant_set", "cfd", "crispr_bulge", "score_status",
        "above_threshold", "row_reported_alignment_found", "n_haplotypes_at_site",
        "n_above_threshold", "n_not_scored", "site_max_crispr_bulge"]


def _row(**kw):
    """A consistent row by default.

    site_max_crispr_bulge used to be hard-coded "NA" here while callers set
    crispr_bulge="1", so the gate's own suite BLESSED the shape the gate was failing to
    catch -- a per-site roll-up saying nothing was measured on a site whose rows carry
    twenty times the operating point. A fixture that encodes the defect is why the hole
    survived: it now defaults to the row's own score, and a test that wants them
    inconsistent must say so explicitly.
    """
    base = {
        "site_key": "chr1:100:+:GUIDE", "variant_set": "chr1_100_A_C",
        "cfd": "0.5", "crispr_bulge": "NA", "score_status": "scored",
        "above_threshold": "unknown", "row_reported_alignment_found": "unknown",
        "n_haplotypes_at_site": "1", "n_above_threshold": "0", "n_not_scored": "0",
    }
    base.update(kw)
    base.setdefault("site_max_crispr_bulge", base["crispr_bulge"])
    return base


class GateFixture(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="hapgate_")
        self.addCleanup(shutil.rmtree, self.dir, True)

    def _write(self, rows, comments=None):
        path = os.path.join(self.dir, "c.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            for c in (HEADER_COMMENTS if comments is None else comments):
                fh.write(c + "\n")
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(r.get(c, "") for c in COLS) + "\n")
        return path

    def _names(self, rows, comments=None):
        return set(v.name for v in chk.check(self._write(rows, comments)))


class TestACorrectFilePasses(GateFixture):
    """The half that stops the gate becoming a rubber stamp."""

    def test_a_plain_correct_file_has_no_violations(self):
        self.assertEqual(self._names([_row()]), set())

    def test_one_haplotype_on_several_frames_is_not_a_counter_error(self):
        """The real-data false alarm. n_above_threshold counts HAPLOTYPES; the same
        variant set appears once per frame that named it, and twice on one frame under
        different alignments, so rows legitimately outnumber the counter."""
        rows = [_row(crispr_bulge="1", above_threshold="yes", n_above_threshold="1",
                     site_max_crispr_bulge="1"),
                _row(crispr_bulge="0.9", above_threshold="yes", n_above_threshold="1",
                     site_max_crispr_bulge="1"),
                _row(crispr_bulge="0.5", above_threshold="yes", n_above_threshold="1",
                     site_max_crispr_bulge="1")]
        self.assertEqual(self._names(rows), set())

    def test_two_distinct_haplotypes_above_with_a_counter_of_two(self):
        rows = [_row(variant_set="chr1_100_A_C", crispr_bulge="1",
                     above_threshold="yes", n_haplotypes_at_site="2",
                     n_above_threshold="2", site_max_crispr_bulge="1"),
                _row(variant_set="chr1_101_G_T", crispr_bulge="0.9",
                     above_threshold="yes", n_haplotypes_at_site="2",
                     n_above_threshold="2", site_max_crispr_bulge="1")]
        self.assertEqual(self._names(rows), set())


class TestEachInvariantActuallyFires(GateFixture):

    def test_a_blank_cell(self):
        self.assertIn("blank-cell", self._names([_row(cfd="")]))

    def test_the_sentinel_printed(self):
        self.assertIn("sentinel-printed", self._names([_row(crispr_bulge="-1.0")]))

    def test_above_threshold_on_an_unscored_row(self):
        got = self._names([_row(crispr_bulge="1", above_threshold="yes",
                                score_status="needs-realignment:chr1:100",
                                n_above_threshold="1")])
        self.assertIn("above-but-unscored", got)

    def test_a_verdict_with_no_number_beside_it(self):
        self.assertIn("verdict-without-number",
                      self._names([_row(crispr_bulge="NA", above_threshold="yes",
                                        n_above_threshold="1")]))

    def test_a_verdict_that_contradicts_its_number(self):
        got = self._names([_row(crispr_bulge="1", above_threshold="no")])
        self.assertIn("verdict-disagrees", got)

    def test_a_value_the_header_does_not_explain(self):
        self.assertIn("undeclared-verdict",
                      self._names([_row(above_threshold="probably")]))

    def test_a_total_below_what_was_written(self):
        rows = [_row(variant_set="chr1_100_A_C", n_haplotypes_at_site="1"),
                _row(variant_set="chr1_101_G_T", n_haplotypes_at_site="1")]
        self.assertIn("total-below-shown", self._names(rows))

    def test_a_counter_exceeding_the_total(self):
        self.assertIn("counter-exceeds-total",
                      self._names([_row(n_haplotypes_at_site="1",
                                        n_not_scored="5")]))

    def test_more_haplotypes_above_than_the_counter_claims(self):
        rows = [_row(variant_set="chr1_100_A_C", crispr_bulge="1",
                     above_threshold="yes", n_haplotypes_at_site="2",
                     n_above_threshold="1"),
                _row(variant_set="chr1_101_G_T", crispr_bulge="1",
                     above_threshold="yes", n_haplotypes_at_site="2",
                     n_above_threshold="1")]
        self.assertIn("above-exceeds-counter", self._names(rows))

    def test_rows_of_one_site_disagreeing_about_that_site(self):
        rows = [_row(variant_set="chr1_100_A_C", n_haplotypes_at_site="2"),
                _row(variant_set="chr1_101_G_T", n_haplotypes_at_site="3")]
        self.assertIn("site-disagrees-with-itself", self._names(rows))

    def test_an_empty_file_is_not_a_clean_file(self):
        self.assertIn("no-rows", self._names([]))

    def test_a_file_with_no_header(self):
        path = os.path.join(self.dir, "empty.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        self.assertIn("no-header", set(v.name for v in chk.check(path)))


class TestTheExitCodeMatchesTheVerdict(GateFixture):

    def test_a_clean_file_exits_zero(self):
        self.assertEqual(chk.main(["--tsv", self._write([_row()]), "--quiet"]), 0)

    def test_a_violating_file_exits_one(self):
        path = self._write([_row(cfd="")])
        self.assertEqual(chk.main(["--tsv", path, "--quiet"]), 1)

    def test_a_missing_file_exits_two_not_zero(self):
        self.assertEqual(
            chk.main(["--tsv", os.path.join(self.dir, "nope.tsv"), "--quiet"]), 2)


class TestTheGateSaysWhenItCannotCheck(GateFixture):
    """A check that silently skips is worse than one that fails: it reports success for
    work it did not do. The threshold parser matched only within a single line, while
    the real header puts the phrase on one line and the number on the next -- so on
    every real file the above_threshold-vs-score check was inert, and the gate passed
    460,714 rows partly vacuously."""

    def test_a_header_with_no_threshold_is_a_violation_not_a_pass(self):
        got = self._names([_row()], comments=["# a header that states no threshold"])
        self.assertIn("threshold-not-stated", got)

    def test_the_real_two_line_header_shape_is_read(self):
        got = self._names([_row(crispr_bulge="1", above_threshold="no")])
        self.assertNotIn("threshold-not-stated", got)
        self.assertIn("verdict-disagrees", got)

    def test_an_explicit_threshold_overrides_a_silent_header(self):
        path = self._write([_row(crispr_bulge="1", above_threshold="no")],
                           comments=["# no threshold here"])
        names = set(v.name for v in chk.check(path, threshold=0.00019))
        self.assertNotIn("threshold-not-stated", names)
        self.assertIn("verdict-disagrees", names)


class TestTheRollUpMustMatchTheRows(GateFixture):
    """A cap only ever REMOVES rows, so the per-site CRISPR-Bulge roll-up can never say
    less than the rows under it. This was checked nowhere: on a site whose rows carry
    0.004 -- twenty times the operating point -- setting site_max_crispr_bulge to NA
    passed with "OK: every invariant holds"."""

    def test_a_roll_up_of_na_over_scored_rows_fails(self):
        got = self._names([_row(crispr_bulge="0.004", above_threshold="yes",
                                n_above_threshold="1", site_max_crispr_bulge="NA")])
        self.assertIn("site-max-missing", got)

    def test_a_roll_up_below_the_rows_fails(self):
        got = self._names([_row(crispr_bulge="0.004", above_threshold="yes",
                                n_above_threshold="1",
                                site_max_crispr_bulge="0.000001")])
        self.assertIn("site-max-below-rows", got)

    def test_a_roll_up_above_the_rows_is_fine(self):
        """A cap legitimately hides the higher-scoring haplotype."""
        got = self._names([_row(crispr_bulge="0.004", above_threshold="yes",
                                n_haplotypes_at_site="2", n_above_threshold="2",
                                site_max_crispr_bulge="0.9")])
        self.assertEqual(got, set())

    def test_a_site_with_no_numbers_is_not_flagged(self):
        self.assertEqual(self._names([_row(crispr_bulge="NA",
                                           site_max_crispr_bulge="NA")]), set())


class TestTheSentinelCheckFiresOnEverySpelling(GateFixture):
    """It tested the literal "-1.0" while %.6g renders -1.0 as "-1" and the pipeline
    writes "-1.000" -- so it could not fire on anything anything actually produces."""

    def test_each_spelling_is_caught(self):
        for spelling in ("-1", "-1.0", "-1.000", "-1.0000"):
            got = self._names([_row(crispr_bulge=spelling,
                                    site_max_crispr_bulge=spelling)])
            self.assertIn("sentinel-printed", got, "missed %r" % spelling)

    def test_a_real_low_score_is_not_mistaken_for_it(self):
        got = self._names([_row(crispr_bulge="0.000001",
                                site_max_crispr_bulge="0.000001")])
        self.assertNotIn("sentinel-printed", got)

    def test_a_non_numeric_score_is_flagged(self):
        got = self._names([_row(crispr_bulge="probably",
                                site_max_crispr_bulge="probably")])
        self.assertIn("score-not-a-number", got)


if __name__ == "__main__":
    unittest.main(verbosity=2)

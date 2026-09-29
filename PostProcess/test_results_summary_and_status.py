"""Tests for two web-page behaviours found while running a reference-only,
non-human (pig) search:

* the Result Summary matrix showed bulge rows 0-5 for a 1 DNA + 1 RNA search
  (``results_page._drop_empty_variant_block``);
* a search that finds nothing left every later step as "To do" forever, with no
  way to tell it from a hang (``load_page.refresh_search`` "no hits" message).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_results_summary_and_status.py' -v
"""

import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    raise unittest.SkipTest("dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

import pages.load_page as lp  # noqa: E402
import pages.results_page as rp  # noqa: E402


def _counts(rows):
    return pd.DataFrame(rows, columns=["0", "1", "2", "3", "4"])


REF_BLOCK = [[1, 0, 1, 2, 32], [0, 1, 4, 74, 2099], [0, 0, 20, 400, 16067]]
ZERO_BLOCK = [[0] * 5] * 3


class TestDropEmptyVariantBlock(unittest.TestCase):
    def test_reference_only_keeps_just_the_reference_block(self):
        out = rp._drop_empty_variant_block(_counts(REF_BLOCK + ZERO_BLOCK), "ref")
        self.assertEqual(len(out), 3)  # bulges 0..2, not 0..5
        self.assertEqual(out.values.tolist(), REF_BLOCK)
        self.assertEqual(list(out.index), [0, 1, 2])

    def test_search_with_variants_is_untouched(self):
        df = _counts(REF_BLOCK + [[0, 0, 0, 4, 11], [0, 2, 37, 843, 34478], [0, 0, 10, 122, 3414]])
        self.assertEqual(len(rp._drop_empty_variant_block(df, "both")), 6)

    def test_real_data_in_the_second_half_is_never_discarded(self):
        df = _counts(REF_BLOCK + [[0, 0, 0, 0, 1]] + ZERO_BLOCK[:2])
        self.assertEqual(len(rp._drop_empty_variant_block(df, "ref")), 6)

    def test_other_genome_types_and_odd_shapes_untouched(self):
        df = _counts(REF_BLOCK + ZERO_BLOCK)
        self.assertEqual(len(rp._drop_empty_variant_block(df, "var")), 6)
        odd = _counts(REF_BLOCK)  # 3 rows: cannot be two equal blocks
        self.assertEqual(len(rp._drop_empty_variant_block(odd, "ref")), 3)
        single = _counts(REF_BLOCK[:1])
        self.assertEqual(len(rp._drop_empty_variant_block(single, "ref")), 1)

    def test_larger_bulge_budget(self):
        block = [[i, 0, 0, 0, 0] for i in range(5)]  # e.g. 2 DNA + 2 RNA -> 5 rows
        df = _counts(block + [[0] * 5] * 5)
        self.assertEqual(len(rp._drop_empty_variant_block(df, "ref")), 5)


class TestNoHitsStatus(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        os.makedirs(os.path.join(self.tmp, lp.RESULTS_DIR))

    def _job(self, name, log_text):
        d = os.path.join(self.tmp, lp.RESULTS_DIR, name)
        os.makedirs(d)
        with open(os.path.join(d, ".Params.txt"), "w") as f:
            f.write("Genome_selected\tsusScr11\nRef_comp\tFalse\n")
        with open(os.path.join(d, ".guides.txt"), "w") as f:
            f.write("ACACCCCCCAGGTTTTTGTGNNN\n")
        with open(os.path.join(d, lp.LOG_FILE), "w") as f:
            f.write(log_text)
        return d

    def _status(self, name):
        with patch.object(lp, "current_working_directory", self.tmp + os.sep):
            return lp.refresh_search(1, f"?job={name}")

    LOG_EMPTY = (
        "Job\tStart\tFri\nIndex-genome Reference\tEnd\tFri\nOff-targets search\tStart\tFri\n"
        "Off-targets search Reference\tEnd\tFri\nOff-targets search\tEnd\tFri\n"
        "Post-analysis\tStart\tFri\nPost-analysis\tEnd\tFri\n"
        "No off-targets found\nJob\tEnd\tFri\n"
    )

    def test_genuine_zero_hit_job_says_so_and_offers_no_results_link(self):
        self._job("nohits", self.LOG_EMPTY)
        out = self._status("nohits")
        self.assertEqual(out[0], {"visibility": "hidden"})  # no VIEW RESULTS button
        self.assertEqual(out[8], "")  # no results link
        alert = out[9]
        self.assertIn("no off-targets were found", str(alert.children))
        self.assertEqual(alert.color, "info")

    def test_unfinished_job_is_not_reported_as_no_hits(self):
        self._job("running", "Job\tStart\tFri\nOff-targets search\tStart\tFri\n")
        out = self._status("running")
        self.assertNotIn("no off-targets", str(out[9]).lower())

    def test_normal_finished_job_still_shows_results(self):
        log = self.LOG_EMPTY.replace("No off-targets found\n", "").replace(
            "Job\tEnd\tFri\n", "Creating database\tEnd\tFri\nJob\tDone\tFri\n"
        )
        self._job("finished", log)
        out = self._status("finished")
        self.assertEqual(out[0], {"visibility": "visible"})
        self.assertIn("result?job=finished", out[8])


if __name__ == "__main__":
    unittest.main()

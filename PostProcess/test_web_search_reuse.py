"""The web form reuses an earlier search only when it is the same search.

change_url sends a submission to an earlier job when .Params.txt and the guides
match. Those do not record the analysis mode, the full-report option or which
annotations are enabled (.Params.txt names only the assembled ".active" bed), so a
Population-level request went to an earlier per-sample job (and the reverse), a
full-report request to a slim one, a search after an annotation change to the old
annotations, and a new job name was dropped -- all without a word. Each job now
records those (and the CRISPRme version) in .reuse_key.json, a job is reused only
when they match too, and its status page is opened with "#reused", which says so.
Identical resubmissions are still reused.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_search_reuse.py' -v
"""

import json
import os
import shutil
import sys
import tempfile
import unittest
import warnings
from unittest.mock import MagicMock, patch

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    raise unittest.SkipTest("dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

import pages.load_page as lp  # noqa: E402
import pages.main_page as mp  # noqa: E402
import pages.pages_utils as pu  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCAC"


class TestSearchReuse(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)  # files the page code leaves open
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        for d in ("Results", "PAMs", "Annotations", "Genomes/hg38"):
            os.makedirs(os.path.join(self.tmp, d))
        with open(os.path.join(self.tmp, "PAMs", "20bp-NRG-SpCas9.txt"), "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        for name, rows in (("a.hg38.bed", "chr1\t10\t500\tpromoter\n"),
                           ("b.hg38.bed", "chr1\t900\t1500\tenhancer\n")):
            with open(os.path.join(self.tmp, "Annotations", name), "w") as f:
                f.write(rows)
        for mod in (mp, pu):
            p = patch.object(mod, "current_working_directory", self.tmp + os.sep)
            p.start()
            self.addCleanup(p.stop)
        self.pool = MagicMock()
        p = patch.object(mp, "pool_executor", self.pool)
        p.start()
        self.addCleanup(p.stop)
        pu.write_enabled_annotations("hg38", ["a.hg38.bed"])

    def submit(self, *, mode="per-sample", full=False, name=None):
        """One click on Submit, as the form sends it (Advanced options open)."""
        return mp.change_url(
            1, "http://127.0.0.1:8080/index", "hg38", "ref", "20bp-NRG-SpCas9", "IP",
            GUIDE, 4, 1, 1, "N", None, None, None, [], None, name, 6, True,
            mode, ["full"] if full else [],
        )

    def finish(self, url):
        """The job of a submission runs to its end, cleanly."""
        job = url[1].split("job=")[-1]
        with open(os.path.join(self.tmp, "Results", job, "log.txt"), "w") as f:
            f.write("Job\tStart\tFri\nCreating database\tEnd\tFri\nJob\tDone\tFri\n")
        return job

    def jobs(self):
        return sorted(os.listdir(os.path.join(self.tmp, "Results")))

    def test_identical_search_is_reused_and_says_so(self):
        first = self.finish(self.submit())
        self.assertEqual(self.pool.submit.call_count, 1)
        again = self.submit()
        self.assertEqual(again, ("/load", f"?job={first}", lp.REUSED_HASH))
        self.assertEqual(self.pool.submit.call_count, 1)  # not run again
        self.assertEqual(self.jobs(), [first])  # the new job directory was removed
        note = lp.show_reused_search_note(lp.REUSED_HASH)
        self.assertIn("not run again", note.children)
        self.assertEqual(lp.show_reused_search_note(""), "")

    def test_each_difference_runs_a_new_search(self):
        first = self.finish(self.submit())
        for label, kwargs in (
            ("population-level", dict(mode="population-level")),
            ("full report", dict(full=True)),
            ("job name", dict(name="other")),
        ):
            with self.subTest(differs=label):
                url = self.submit(**kwargs)
                self.assertNotEqual(url[1], f"?job={first}")
                self.assertEqual(url[2], "")
                job = self.finish(url)
                with open(os.path.join(self.tmp, "Results", job, ".search_mode")) as f:
                    self.assertEqual(f.read(), kwargs.get("mode", "per-sample"))
                # ... and an identical one of THAT is reused
                self.assertEqual(self.submit(**kwargs), ("/load", f"?job={job}", lp.REUSED_HASH))
        self.assertEqual(self.pool.submit.call_count, 4)

    def test_annotation_change_runs_a_new_search(self):
        first = self.finish(self.submit())
        pu.write_enabled_annotations("hg38", ["a.hg38.bed", "b.hg38.bed"])
        url = self.submit()
        self.assertNotEqual(url[1], f"?job={first}")
        with open(os.path.join(self.tmp, "Results", self.finish(url), mp.REUSE_KEY_FILE)) as f:
            self.assertEqual(json.load(f)["annotations"], ["a.hg38.bed", "b.hg38.bed"])

    def test_job_without_a_reuse_key_is_not_reused(self):
        # a job from before 2.8.1: what it ran with is not known
        first = self.finish(self.submit())
        os.remove(os.path.join(self.tmp, "Results", first, mp.REUSE_KEY_FILE))
        self.assertNotEqual(self.submit()[1], f"?job={first}")

    def test_failed_job_is_not_reused(self):
        first = self.finish(self.submit())
        with open(os.path.join(self.tmp, "Results", first, "log_error.txt"), "w") as f:
            f.write("ERROR: something\n")
        self.assertNotEqual(self.submit()[1], f"?job={first}")

    def test_queued_identical_search_is_joined_without_a_stray_job(self):
        first = self.submit()[1].split("job=")[-1]  # submitted, not started (no log)
        self.assertEqual(self.submit(), ("/load", f"?job={first}", lp.REUSED_HASH))
        self.assertEqual(self.jobs(), [first])  # no job directory that never runs

    def test_reuse_key_records_what_params_do_not(self):
        job = self.finish(self.submit(mode="population-level", full=True, name="my job"))
        with open(os.path.join(self.tmp, "Results", job, mp.REUSE_KEY_FILE)) as f:
            key = json.load(f)
        self.assertEqual(key["search_mode"], "population-level")
        self.assertTrue(key["report_full"])
        self.assertEqual(key["job_name"], "my_job")
        self.assertEqual(key["annotations"], ["a.hg38.bed"])
        self.assertEqual(key["crisprme_version"], mp._crisprme_version())
        self.assertNotEqual(key["crisprme_version"], "unknown")


if __name__ == "__main__":
    unittest.main()

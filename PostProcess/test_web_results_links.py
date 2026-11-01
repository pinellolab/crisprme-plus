"""Results-page links, found by the 2.8.1 web walkthrough:

* a results link to a job with no results (no such job, failed, refused for an
  unfinished build, found nothing, not finished) failed on the missing results file
  with a server error and a blank page (``results_page._result_not_ready``);
* the job status page held the "view-results" link twice (one id, two elements).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_results_links.py' -v
"""

import os
import shutil
import sqlite3
import sys
import tempfile
import time
import unittest
import warnings
from unittest.mock import patch

import pandas as pd

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
    from dash.development.base_component import Component
except ImportError:  # pragma: no cover
    raise unittest.SkipTest("dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

import pages.load_page as lp  # noqa: E402
import pages.results_page as rp  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCACNNN"
CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _walk(component):
    """Every component in a Dash layout."""
    yield component
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for child in children:
            yield from _walk(child)
    elif isinstance(children, Component):
        yield from _walk(children)


def _by_id(layout, wanted):
    items = layout if isinstance(layout, (list, tuple)) else [layout]
    for item in items:
        for c in _walk(item):
            if getattr(c, "id", None) == wanted:
                return c
    return None


def _text(layout):
    items = layout if isinstance(layout, (list, tuple)) else [layout]
    out = []
    for item in items:
        for c in _walk(item):
            if isinstance(c, str):
                out.append(c)
            elif isinstance(getattr(c, "children", None), str):
                out.append(c.children)
    return " ".join(out)


class _Tmp(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)  # files the page code leaves open
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        os.makedirs(os.path.join(self.tmp, rp.RESULTS_DIR))
        cwd = patch.object(rp, "current_working_directory", self.tmp + os.sep)
        cwd.start()
        self.addCleanup(cwd.stop)
        for mod in ("pages.pages_utils",):
            p = patch.object(sys.modules[mod], "current_working_directory", self.tmp + os.sep)
            p.start()
            self.addCleanup(p.stop)

    def _job(self, job, *, log=None, log_error=None, queued=False, finished=True):
        d = os.path.join(self.tmp, rp.RESULTS_DIR, job)
        os.makedirs(d)
        with open(os.path.join(d, ".Params.txt"), "w") as f:
            f.write("Genome_selected\thg38\nGenome_ref\thg38\nGenome_idx\tNRG_2_hg38+hg38_x\n"
                    "Pam\tNRG\nMax_bulges\t2\nMismatches\t4\nDNA\t1\nRNA\t1\n"
                    "Max_total_edits\t6\nThreshold_mode\tadvanced\nNuclease\tSpCas9\n"
                    "Ref_comp\tTrue\n")
        with open(os.path.join(d, ".guides.txt"), "w") as f:
            f.write(GUIDE + "\n")
        if log is not None:
            with open(os.path.join(d, "log.txt"), "w") as f:
                f.write(log)
        if log_error is not None:
            with open(os.path.join(d, "log_error.txt"), "w") as f:
                f.write(log_error)
        if queued:
            open(os.path.join(d, "queue.txt"), "w").close()
        if finished:
            self._finished_files(job, d)
        return d

    def _finished_files(self, job, d):
        """What the results pages read for a finished variant job (a minimal copy)."""
        rows = []
        for i, (mmb, samples) in enumerate(((0, "S1,S2"), (0, "S1"), (3, "S2"), (4, "NA"))):
            row = {"Spacer+PAM": GUIDE, "Chromosome": "chr1"}
            for crit in CRITERIA:
                row[f"Start_coordinate_({crit})"] = 1000 * (i + 1)
                row[f"Mismatches_({crit})"] = mmb
                row[f"Bulges_({crit})"] = 0
                row[f"Bulge_type_({crit})"] = "X"
                row[f"Mismatches+bulges_({crit})"] = mmb
                row[f"Variant_samples_({crit})"] = samples
            row["CFD_score_(highest_CFD)"] = 1.0 - 0.2 * i
            row["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"] = 1.0 - 0.2 * i
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(d, f"{GUIDE}+NRG_hg38+hg38_x_4+2_integrated_results.tsv"),
                  sep="\t", index=False)
        with sqlite3.connect(os.path.join(d, f".{job}.db")) as conn:
            df.to_sql("final_table", conn, index=False)
        with open(os.path.join(d, ".sampleID.txt"), "w") as f:
            f.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                    "S1\tPOPA\tSUPX\tmale\nS2\tPOPB\tSUPX\tfemale\n")
        for crit in ("CFD", "CRISPR-Bulge", "fewest"):
            with open(os.path.join(d, f"{job}.summary_by_samples.{GUIDE}_{crit}.txt"), "w") as f:
                f.write(f"{GUIDE}\nSample\tSex\tPopulation\tSuper Population\tTargets in Sample\t"
                        "Targets in Population\tTargets in Super Population\tPAM Creation\n"
                        "S1\tmale\tPOPA\tSUPX\t2\t2\t3\t0\nS2\tfemale\tPOPB\tSUPX\t2\t2\t3\t1\n")
        with open(os.path.join(d, ".dropdown.json"), "w") as f:
            f.write("CFD")


class TestResultLinkWithoutResults(_Tmp):
    REFUSED = ("/data/Genomes/hg38+x is unfinished: a build-index-only run is writing it, "
               "or stopped before it was done. Let that run finish, or run build-index-only "
               "again to its end, then run this search again")

    def _pages(self, job):
        """Every results route for a job: the page and its three sub-pages."""
        return {
            "result": rp.result_page(job),
            "guide": rp.guidePagev3(job, f"{GUIDE}newX00"),
            "sample": rp.sample_page(job, f"{GUIDE}-Sample-S1"),
            "region": rp.cluster_page(job, f"{GUIDE}-Pos-chr1-1000"),
        }

    def test_unknown_job(self):
        for route, page in self._pages("NOSUCHJOB").items():
            with self.subTest(route=route):
                self.assertIn("The selected result does not exist", _text(page))

    def test_job_without_results_says_why_and_links_its_status_page(self):
        cases = {
            "failed": (dict(log="Job\tStart\tFri\n", log_error="ERROR: search failed\n"),
                       "The selected result encountered some errors"),
            "refused": (dict(log="Job\tStart\tFri\n", log_error=f"ERROR: {self.REFUSED}\n"),
                        f"The search was not run: {self.REFUSED}."),
            "nohits": (dict(log="Job\tStart\tFri\nNo off-targets found\nJob\tEnd\tFri\n"),
                       "no off-targets were found"),
            "queued": (dict(queued=True), "has not finished yet"),
            "running": (dict(log="Job\tStart\tFri\nOff-targets search\tStart\tFri\n"),
                        "has not finished yet"),
        }
        for job, (kwargs, says) in cases.items():
            self._job(job, finished=False, **kwargs)
            for route, page in self._pages(job).items():
                with self.subTest(job=job, route=route):
                    self.assertIn(says, _text(page))
                    links = [c.href for c in _walk(page) if getattr(c, "href", None)]
                    self.assertIn(f"/load?job={job}", links)

    def test_finished_job_still_gets_its_results_page(self):
        self._job("done", log="Job\tDone\tFri\n")
        with open(os.path.join(self.tmp, rp.RESULTS_DIR, "done", ".done.acfd_CFD.txt"), "w") as f:
            f.write(f"{GUIDE}\t99\t1\n")
        self.assertIsNone(rp._result_not_ready("done"))
        page = rp.guidePagev3("done", f"{GUIDE}newX00")
        self.assertIsNotNone(_by_id(page, "table-subset-target"))


class TestStatusPageIds(unittest.TestCase):
    def test_no_id_twice(self):
        for name, layout in (("load_page", lp.load_page("x")),
                             ("load_page_assembly", lp.load_page_assembly("x"))):
            with self.subTest(page=name):
                ids = [c.id for c in _walk(layout) if getattr(c, "id", None)]
                self.assertEqual(sorted(i for i in set(ids) if ids.count(i) > 1), [])


if __name__ == "__main__":
    unittest.main()

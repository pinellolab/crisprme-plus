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

import gzip
import json
import os
import random
import shutil
import sys
import tempfile
import unittest
import warnings
from unittest.mock import MagicMock, patch

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
except ImportError:  # pragma: no cover
    dash = None
# Without Dash each class skips, not the module: a SkipTest raised at import is a traceback
# under `python -m unittest <module>` (skips D7).
_NEEDS_DASH = unittest.skipUnless(dash, "dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

if dash:  # (the pages import Dash)
    import pages.load_page as lp  # noqa: E402
    import pages.main_page as mp  # noqa: E402
    import pages.pages_utils as pu  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCAC"


@_NEEDS_DASH
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

    def submit(self, *, mode="per-sample", full=False, name=None, variant="ref",
               gate=False):
        """One click on Submit, as the form sends it (Advanced options open)."""
        return mp.change_url(
            1, "http://127.0.0.1:8080/index", "hg38", variant, "20bp-NRG-SpCas9", "IP",
            GUIDE, 4, 1, 1, "N", None, None, None, [], None, name, 6, True,
            mode, ["full"] if full else [], ["on"] if gate else [],
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

    def test_the_reused_job_s_log_keeps_its_lines_and_takes_today_s_date(self):
        # a reuse wrote the log back without its last line, the report's (the search
        # writes it after "Job Done"), and with a second "Job Done" with no line end
        first = self.finish(self.submit())
        log = os.path.join(self.tmp, "Results", first, "log.txt")
        lines = ["Job\tStart\tFri", "Creating database\tEnd\tFri", "Job\tDone\tFri",
                 "Building shareable report\tStart\tFri", "Building shareable report\tEnd\tFri"]
        with open(log, "w") as f:
            f.write("\n".join(lines) + "\n")
        self.assertEqual(self.submit()[1], f"?job={first}")
        with open(log) as f:
            again = f.read()
        self.assertTrue(again.endswith("\n"))
        again = again.split("\n")[:-1]
        self.assertEqual(len(again), len(lines))
        self.assertEqual([l for i, l in enumerate(again) if i != 2], lines[:2] + lines[3:])
        self.assertTrue(again[2].startswith("Job\tDone\t"))
        self.assertNotEqual(again[2], "Job\tDone\tFri")  # (today's date)

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

    def test_the_carrier_gate_check_box_reaches_the_search_and_is_a_difference(self):
        """The check box that asks the carrier gate to check a sites-only panel's rows
        against the phased HGDP + 1000 Genomes genotypes (complete-search
        --prioritize-observed) goes to the search as CRISPRME_PRIORITIZE_OBSERVED, and a job
        run with it is not reused for one without it (nor the reverse)."""
        first = self.finish(self.submit())
        self.assertIn(" CRISPRME_PRIORITIZE_OBSERVED=auto ", self.pool.submit.call_args[0][1])
        url = self.submit(gate=True)
        self.assertNotEqual(url[1], f"?job={first}")
        self.assertIn(" CRISPRME_PRIORITIZE_OBSERVED=on ", self.pool.submit.call_args[0][1])
        job = self.finish(url)
        with open(os.path.join(self.tmp, "Results", job, mp.REUSE_KEY_FILE)) as f:
            self.assertTrue(json.load(f)["prioritize_observed"])
        self.assertEqual(self.submit(gate=True), ("/load", f"?job={job}", lp.REUSED_HASH))
        self.assertEqual(self.submit(), ("/load", f"?job={first}", lp.REUSED_HASH))

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


@_NEEDS_DASH
class TestReuseOnceTheSearchHasRun(unittest.TestCase):
    """An identical request is sent to an earlier job that the search has run to its end,
    as the job left its folder. 2.10.0's first cut had the search record the indexes it took
    in the job's .Params.txt, which change_url compares byte for byte with a new request's,
    so every identical resubmission started a new search: the record is a file of its own.
    The search the image runs, its script and its engine, on a small genome."""

    def setUp(self):
        TestSearchReuse.setUp(self)
        # the guide's site and two of one mismatch, on a chromosome of 30 kb
        rng = random.Random(7)
        seq = [rng.choice("ACGT") for _ in range(30000)]
        for i, at in enumerate((1000, 9000, 20000)):
            site = list(GUIDE + "AGG")
            if i:
                site[3 + i] = "A" if site[3 + i] != "A" else "C"
            seq[at:at + len(site)] = site
        with open(os.path.join(self.tmp, "Genomes", "hg38", "chr1.fa"), "w") as f:
            f.write(">chr1\n" + "".join(seq) + "\n")
        os.makedirs(os.path.join(self.tmp, "samplesIDs"))
        os.makedirs(os.path.join(self.tmp, "VCFs"))
        # the job runs here, to its end, as the web's pool runs it
        self.pool.submit.side_effect = lambda fn, *args, **kwargs: fn(*args, **kwargs)

    def test_an_identical_request_is_sent_to_it(self):
        url = TestSearchReuse.submit(self)
        job = url[1].split("job=")[-1]
        folder = os.path.join(self.tmp, "Results", job)
        with open(os.path.join(folder, "log.txt")) as f:
            self.assertIn("Job\tDone", f.read())
        # (a finished search moves its stderr log aside, as log_error_no_check.txt)
        self.assertFalse(os.path.exists(os.path.join(folder, "log_error.txt")))
        self.assertEqual(TestSearchReuse.submit(self), ("/load", f"?job={job}", lp.REUSED_HASH))
        self.assertEqual(self.pool.submit.call_count, 1)
        self.assertEqual(TestSearchReuse.jobs(self), [job])
        # what the search recorded of the index it took (it built it, under the name that
        # says its guide length and side), beside the request's parameters
        with open(os.path.join(folder, ".indexes_searched.txt")) as f:
            self.assertIn("Reference_index\t20-NRG_2_hg38 (in ", f.read())
        with open(os.path.join(folder, mp.PARAMS_FILE)) as f:
            self.assertNotIn("_index", f.read())


@_NEEDS_DASH
class TestDatasetsWithoutTheTiers(unittest.TestCase):
    """A search of a variant dataset with neither Dictionaries/registry_<vcf>/ nor
    genotypes_<vcf>/ (an index built without --samplesID, such as the one the benchmark
    complete-test builds) reports the allele frequencies of its classic dictionaries (an
    indel's, of its indel log), which CRISPRitz add-variants copied from its VCF, and the form
    runs it population-level (it offers per-sample only with a genotype store). complete-search
    records such datasets (.untiered_datasets, beside .variant_datasets), so that the report
    says so, and what --per-sample gives there; the web job, which launches the search
    script itself, now records them too."""

    setUp, submit, finish = TestSearchReuse.setUp, TestSearchReuse.submit, TestSearchReuse.finish
    jobs = TestSearchReuse.jobs

    def records(self, url):
        """The job's .untiered_datasets and .variant_datasets (None: not written)."""
        job = os.path.join(self.tmp, "Results", self.finish(url))
        out = []
        for name in (".untiered_datasets", ".variant_datasets"):
            try:
                with open(os.path.join(job, name)) as f:
                    out.append(f.read())
            except FileNotFoundError:
                out.append(None)
        return out

    def test_a_dataset_without_the_tiers_is_recorded_as_complete_search_records_it(self):
        os.makedirs(os.path.join(self.tmp, "VCFs", "hg38_ds"))
        os.makedirs(os.path.join(self.tmp, "samplesIDs"))
        with open(os.path.join(self.tmp, "samplesIDs", "hg38_ds.samplesID.txt"), "w") as f:
            f.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\nS1\tGBR\tEUR\tmale\n")
        self.assertEqual(self.records(self.submit(mode="population-level", variant="ds")),
                         ["hg38_ds\n", "hg38_ds\n"])
        # with the tiers, or with no dataset, nothing
        os.makedirs(os.path.join(self.tmp, "Dictionaries", "registry_hg38_ds"))
        self.assertEqual(self.records(self.submit(mode="population-level", variant="ds",
                                                  name="tiered")), [None, None])
        self.assertEqual(self.records(self.submit(mode="population-level")), [None, None])

    def test_a_sites_only_dataset_without_the_tiers_is_not_recorded(self):
        # its VCFs have no sample columns: with no allele-frequency registry its search
        # reports none of its variants, so the report says nothing of it. The form offers one
        # dataset at a time, and Submit runs only a dataset it offers, so a request for two
        # ("gt+so") makes no job: the case beside a dataset whose VCFs have samples is
        # complete-search's (test_build_index_only_atomic)
        os.makedirs(os.path.join(self.tmp, "samplesIDs"))
        head = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
        for name, columns, genotypes, listed in (
                ("so", "", "", ""), ("gt", "\tFORMAT\tS1", "\tGT\t0|1", "S1\tGBR\tEUR\tmale\n")):
            os.makedirs(os.path.join(self.tmp, "VCFs", f"hg38_{name}"))
            with gzip.open(os.path.join(self.tmp, "VCFs", f"hg38_{name}",
                                        f"{name}.chr1.vcf.gz"), "wt") as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                        f"{columns}\nchr1\t11\trs1\tA\tG\t.\tPASS\tAF=0.25{genotypes}\n")
            with open(os.path.join(self.tmp, "samplesIDs", f"hg38_{name}.samplesID.txt"),
                      "w") as f:
                f.write(head + listed)
        self.assertEqual(self.records(self.submit(mode="population-level", variant="so")),
                         [None, None])
        self.assertEqual(self.records(self.submit(mode="population-level", variant="gt")),
                         ["hg38_gt\n", "hg38_gt\n"])
        from dash.exceptions import PreventUpdate
        before = self.jobs()
        with self.assertRaises(PreventUpdate):
            self.submit(mode="population-level", variant="gt+so")
        self.assertEqual(self.jobs(), before)


if __name__ == "__main__":
    unittest.main()

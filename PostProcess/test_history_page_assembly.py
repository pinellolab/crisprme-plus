"""Regression test for assembly-search jobs in the History table
(pages/history_page.py).

Before this, the History page crashed for *any* assembly-search job --
including ordinary, web-submitted ones already sitting in a real Results/
directory, not just CLI-only jobs: `read_job_info()` expects a
`"Job\\tStart"` line in `log.txt` that only complete-search's structured
per-stage log ever contains, and `read_params()`'s whitespace `.split()`
mis-parses assembly's 3-column `.Params.txt` before `construct_history_
summary()` KeyErrors on `params["Genome_selected"]` -- a key that doesn't
exist in assembly's schema.

This test extracts the real functions from source (so it needs neither Dash
nor a configured app), the same technique test_sanitize_job_name.py already
uses for pages/main_page.py.
"""

import ast
import os
import re
import tempfile
import unittest
from datetime import datetime
from typing import Dict, Optional

import pandas as pd

import safe_files  # (PostProcess/: how the History page reads a job's files)

_HISTORY_PAGE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pages", "history_page.py"
)
_PAGES_UTILS = os.path.join(os.path.dirname(_HISTORY_PAGE), "pages_utils.py")

_NAMES = [
    "count_guides",
    "read_job_info",
    "_is_assembly_params",
    "_read_assembly_params_file",
    "_assembly_genome_display",
    "_assembly_max_total_edits",
    "_int_or_dash",
    "_assembly_history_row",
]


def _load_functions():
    with open(_HISTORY_PAGE) as fh:  # (left open, it printed a ResourceWarning)
        src = fh.read()
    tree = ast.parse(src)
    # PARAMS_FILE/GUIDES_FILE are module-level constants history_page.py
    # imports from pages_utils.py, which itself needs a live Dash `app` --
    # inlined here with their known, stable values instead of importing it.
    ns = {
        "os": os, "datetime": datetime, "Dict": Dict, "Optional": Optional,
        "PARAMS_FILE": ".Params.txt", "GUIDES_FILE": ".guides.txt", "safe_files": safe_files,
    }
    found = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in _NAMES:
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
            found.add(node.name)
    missing = set(_NAMES) - found
    if missing:
        raise AssertionError(f"Not found in pages/history_page.py: {sorted(missing)}")
    return ns


_NS = _load_functions()


class TestIsAssemblyParams(unittest.TestCase):
    is_assembly_params = staticmethod(_NS["_is_assembly_params"])

    def test_true_for_assembly_params_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, ".Params.txt")
            with open(path, "w") as f:
                f.write("1\tGenome_type\tassembly\n2\tPam\tNGG\n")
            self.assertTrue(self.is_assembly_params(path))

    def test_false_for_complete_search_params_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, ".Params.txt")
            with open(path, "w") as f:
                f.write("Genome_selected\tHG38\nGenome_idx\tNone\n")
            self.assertFalse(self.is_assembly_params(path))

    def test_false_when_file_missing(self):
        self.assertFalse(self.is_assembly_params("/no/such/file"))


class TestReadAssemblyParamsFile(unittest.TestCase):
    read_params = staticmethod(_NS["_read_assembly_params_file"])

    def test_parses_three_column_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, ".Params.txt")
            with open(path, "w") as f:
                f.write("1\tGenome_type\tassembly\n2\tPam\tNGG\n3\tMismatches\t4\n")
            params = self.read_params(path)
            self.assertEqual(
                params, {"Genome_type": "assembly", "Pam": "NGG", "Mismatches": "4"}
            )

    def test_two_column_lines_are_skipped_not_misparsed(self):
        # a complete-search .Params.txt fed through this parser must not
        # silently produce a garbled index->key mapping
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, ".Params.txt")
            with open(path, "w") as f:
                f.write("Genome_selected\tHG38\n")
            self.assertEqual(self.read_params(path), {})


class TestAssemblyGenomeDisplay(unittest.TestCase):
    display = staticmethod(_NS["_assembly_genome_display"])

    def test_legacy_layout_same_individual_collapses(self):
        self.assertEqual(
            self.display("HG01255_paternal", "HG01255_maternal"), "HG01255 (assembly)"
        )

    def test_bundle_layout_same_individual_collapses(self):
        self.assertEqual(
            self.display("HG01255 (paternal)", "HG01255 (maternal)"), "HG01255 (assembly)"
        )

    def test_different_individuals_shown_plainly(self):
        self.assertEqual(
            self.display("HG01255_paternal", "HG002_maternal"),
            "HG01255_paternal / HG002_maternal",
        )

    def test_unrecognized_label_shape_falls_back_plainly(self):
        self.assertEqual(self.display("weird_label", "other_label"), "weird_label / other_label")


class TestAssemblyMaxTotalEdits(unittest.TestCase):
    max_total_edits = staticmethod(_NS["_assembly_max_total_edits"])

    def test_reads_from_paternal_haplotype_params(self):
        with tempfile.TemporaryDirectory() as results_dir:
            pat_dir = os.path.join(results_dir, "job_paternal")
            os.makedirs(pat_dir)
            with open(os.path.join(pat_dir, ".Params.txt"), "w") as f:
                f.write("Mismatches\t4\nMax_total_edits\t5\n")
            params = {"Paternal_dir": "job_paternal", "Maternal_dir": "job_maternal"}
            self.assertEqual(self.max_total_edits(params, results_dir), "5")

    def test_falls_back_to_maternal_when_paternal_missing(self):
        with tempfile.TemporaryDirectory() as results_dir:
            mat_dir = os.path.join(results_dir, "job_maternal")
            os.makedirs(mat_dir)
            with open(os.path.join(mat_dir, ".Params.txt"), "w") as f:
                f.write("Max_total_edits\t6\n")
            params = {"Paternal_dir": "job_paternal", "Maternal_dir": "job_maternal"}
            self.assertEqual(self.max_total_edits(params, results_dir), "6")

    def test_dash_when_nowhere_to_find_it(self):
        with tempfile.TemporaryDirectory() as results_dir:
            self.assertEqual(self.max_total_edits({}, results_dir), "-")


class TestIntOrDash(unittest.TestCase):
    int_or_dash = staticmethod(_NS["_int_or_dash"])

    def test_valid_int_string(self):
        self.assertEqual(self.int_or_dash("4"), 4)

    def test_none_and_garbage_become_dash(self):
        self.assertEqual(self.int_or_dash(None), "-")
        self.assertEqual(self.int_or_dash("not a number"), "-")


class TestAssemblyHistoryRow(unittest.TestCase):
    row_fn = staticmethod(_NS["_assembly_history_row"])

    def _make_job(self, results_dir, jobid, params_lines, guides=("GAAACAGTCGATTTTATCACNNN",)):
        job_dir = os.path.join(results_dir, jobid)
        os.makedirs(job_dir)
        with open(os.path.join(job_dir, ".guides.txt"), "w") as f:
            f.write("\n".join(guides) + "\n")
        return job_dir

    def test_full_row_with_job_start_present(self):
        with tempfile.TemporaryDirectory() as results_dir:
            jobid = "myjob_combined"
            self._make_job(results_dir, jobid, None)
            params = {
                "Genome_paternal": "HG01255_paternal",
                "Genome_maternal": "HG01255_maternal",
                "Mismatches": "4", "DNA": "1", "RNA": "1", "Pam": "NGG",
                "Job_start": "2026-09-22 03:29:00",
                "Paternal_dir": "myjob_paternal", "Maternal_dir": "myjob_maternal",
            }
            row = self.row_fn(jobid, params, results_dir)
            self.assertEqual(row["Job"], jobid)
            self.assertEqual(row["Genome"], "HG01255 (assembly)")
            self.assertEqual(row["Variants"], "-")
            self.assertEqual(row["Mismatches"], 4)
            self.assertEqual(row["DNA bulge"], 1)
            self.assertEqual(row["RNA bulge"], 1)
            self.assertEqual(row["PAM"], "NGG")
            self.assertEqual(row["Number of Guides"], 1)
            self.assertEqual(row["Start"], "2026-09-22 03:29:00")
            self.assertEqual(row["Max edits"], "-")  # no haplotype dirs on disk here

    def test_missing_job_start_falls_back_to_params_file_mtime(self):
        # real case: a real, already-existing web-submitted assembly job
        # created before this fix, whose .Params.txt has no Job_start field
        with tempfile.TemporaryDirectory() as results_dir:
            jobid = "oldjob_combined"
            job_dir = self._make_job(results_dir, jobid, None)
            with open(os.path.join(job_dir, ".Params.txt"), "w") as f:
                f.write("1\tGenome_type\tassembly\n")
            params = {
                "Genome_paternal": "HG01255_paternal", "Genome_maternal": "HG01255_maternal",
                "Mismatches": "4", "DNA": "1", "RNA": "1", "Pam": "NGG",
            }
            row = self.row_fn(jobid, params, results_dir)
            # must be a real, parseable timestamp string, not a crash or "-"
            datetime.strptime(row["Start"], "%Y-%m-%d %H:%M:%S")

    def test_missing_guides_file_does_not_crash(self):
        with tempfile.TemporaryDirectory() as results_dir:
            jobid = "nogiudes_combined"
            os.makedirs(os.path.join(results_dir, jobid))
            row = self.row_fn(jobid, {"Genome_paternal": "?", "Genome_maternal": "?"}, results_dir)
            self.assertEqual(row["Number of Guides"], 0)


class TestMixedStartTimeColumnSortable(unittest.TestCase):
    """Once assembly-search jobs coexist with complete-search jobs in the
    same Results/ dir, construct_history_summary()'s "Start" column mixes
    complete-search's ctime-style string ("Fri 28 Aug 2026 04:09:43 PM UTC")
    with assembly rows' "%Y-%m-%d %H:%M:%S" -- confirmed against real
    production job data to crash pd.to_datetime()'s default format
    inference (locks onto the first row's shape) and, even with
    format="mixed", crash sort_values() on tz-aware vs. tz-naive datetimes
    in the same column. This pins the exact fix (format="mixed", utc=True)
    against the exact real strings that triggered it."""

    def test_mixed_formats_parse_and_sort_without_crashing(self):
        mixed = [
            "Fri 28 Aug 2026 04:09:43 PM UTC",  # complete-search (tz-aware once parsed)
            "2026-09-22 05:17:46",  # assembly-search Job_start (naive)
            "Tue 04 Aug 2026 12:10:05 AM UTC",
        ]
        parsed = pd.to_datetime(mixed, format="mixed", utc=True)
        df = pd.DataFrame({"Start": mixed}).assign(Start=parsed)
        sorted_df = df.sort_values(["Start"], ascending=False)
        self.assertEqual(
            sorted_df["Start"].tolist()[0].strftime("%Y-%m-%d"), "2026-09-22"
        )

    def test_default_format_inference_crashes_on_this_exact_mix(self):
        # documents the failure this test class's fix addresses -- guards
        # against someone reverting to plain pd.to_datetime(col) later.
        # Order matters: pandas infers a format from the first element, so
        # an assembly row (this format) sorted/inserted before a
        # complete-search row is what actually reproduces the crash.
        mixed = ["2026-09-22 05:17:46", "Fri 28 Aug 2026 04:09:43 PM UTC"]
        with self.assertRaises(ValueError):
            pd.to_datetime(mixed)


if __name__ == "__main__":
    unittest.main()


class TestReadJobInfoQueued(unittest.TestCase):
    """A QUEUED job (one submitted while every search slot is busy: 2 by default,
    configurable 1-16) has a log.txt but no 'Job\\tStart' line yet. read_job_info()
    must return None (not raise StopIteration/IOError) so construct_history_summary
    falls back to the Params mtime and the whole History table keeps rendering --
    the bug Manuel reported."""

    read_job_info = staticmethod(_NS["read_job_info"])

    def test_started_job_returns_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "log.txt")
            with open(log, "w") as f:
                f.write("Job\tStart\tFri 28 Aug 2026 04:09:43 PM UTC\n")
            self.assertEqual(
                self.read_job_info(log, "job1"), "Fri 28 Aug 2026 04:09:43 PM UTC"
            )

    def test_queued_job_log_without_start_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "log.txt")
            with open(log, "w") as f:
                f.write("Queued at\tFri 28 Aug 2026 04:00:00 PM UTC\n")  # queue note, no Start
            self.assertIsNone(self.read_job_info(log, "queued_job"))

    def test_missing_log_returns_none(self):
        self.assertIsNone(self.read_job_info("/no/such/log.txt", "new_job"))


class _Tag:
    """Stands in for a dash html component: keeps its children and keywords."""

    def __init__(self, children=None, **kwargs):
        self.children = children
        self.kwargs = kwargs


class _Html:
    A = Td = Th = Tr = Thead = Tbody = Table = Div = Br = Small = _Tag


def _load_history_table():
    """construct_history_summary and history_table_body_, with what they call, from
    pages/history_page.py, the html components replaced by _Tag (no Dash here), and the
    form of a job's id, from pages/pages_utils.py."""
    with open(_HISTORY_PAGE) as fh:  # (left open, it printed a ResourceWarning)
        src = fh.read()
    ns = {
        "os": os, "sys": __import__("sys"), "pd": pd, "datetime": datetime, "Dict": Dict,
        "List": __import__("typing").List, "Optional": Optional,
        "Tuple": __import__("typing").Tuple, "html": _Html, "URL": "",
        "PARAMS_FILE": ".Params.txt", "GUIDES_FILE": ".guides.txt", "LOG_FILE": "log.txt",
        "RESULTS_DIR": "Results", "re": re, "safe_files": safe_files,
    }
    names = set(_NAMES) | {"read_params", "process_genome", "construct_history_summary",
                           "history_table_body_", "_job_link"}
    found = set()
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102 - trusted own source
            found.add(node.name)
        elif isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "SUMMARYTABCOLS" for t in node.targets):
            exec(ast.get_source_segment(src, node), ns)  # noqa: S102
            found.add("SUMMARYTABCOLS")
    with open(_PAGES_UTILS) as fh:
        utils = fh.read()
    for node in ast.parse(utils).body:
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "_JOB_ID_RE" for t in node.targets):
            exec(ast.get_source_segment(utils, node), ns)  # noqa: S102 - trusted own source
            found.add("_JOB_ID_RE")
    # (_job_link is the History cell of a job; the code before it had none)
    missing = (names | {"SUMMARYTABCOLS", "_JOB_ID_RE"}) - found - {"_job_link"}
    if missing:
        raise AssertionError(f"Not found in pages/: {sorted(missing)}")
    return ns


class TestHistoryTableNewestFirst(unittest.TestCase):
    """The History table lists the jobs newest first. construct_history_summary sorted
    them by Start but kept the labels the rows had before the sort, and
    history_table_body_ reads row i by its label: the table showed the jobs in the order
    of the Results/ folder listing (since upstream ab6cae8, 2025; v2.8.0 too)."""

    def test_rows_render_newest_first_with_their_own_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            ns = _load_history_table()
            ns["current_working_directory"] = tmp
            starts = {  # the folder listing order below is neither this one nor its reverse
                "job_mid": "Tue Nov 3 09:52:45 UTC 2026",
                "job_new": "Tue Nov 3 09:53:59 UTC 2026",
                "job_old": "Tue Nov 3 09:49:01 UTC 2026",
                "job_newer": "Tue Nov 3 09:53:10 UTC 2026",
            }
            for jobid, start in starts.items():
                job_dir = os.path.join(tmp, "Results", jobid)
                os.makedirs(job_dir)
                with open(os.path.join(job_dir, ".Params.txt"), "w") as f:
                    f.write("Genome_selected\thg38\nGenome_idx\tNRG_2_hg38\nPam\tNRG\n"
                            "Mismatches\t4\nDNA\t1\nRNA\t1\nMax_total_edits\t6\n")
                with open(os.path.join(job_dir, "log.txt"), "w") as f:
                    f.write(f"Job\tStart\t{start}\n")
                with open(os.path.join(job_dir, ".guides.txt"), "w") as f:
                    f.write("CTAACAGTTGCTTTTATCACNNN\n")
            summary = ns["construct_history_summary"](list(starts))
            rows = ns["history_table_body_"](summary, 1, 1000000, summary.shape[0])
            shown = []
            for row in rows:
                link = row.children[0].children  # the Job cell's link
                start = row.children[-1].children  # the Start cell
                self.assertEqual(link.kwargs["href"], f"load?job={link.children}")
                shown.append((link.children, start))
            self.assertEqual([job for job, _ in shown],
                             ["job_new", "job_newer", "job_mid", "job_old"])
            # each row's Start is its own job's
            for job, start in shown:
                self.assertEqual(start, pd.to_datetime(starts[job], utc=True))


class TestHistoryJobsItsPagesDoNotOpen(unittest.TestCase):
    """A job's pages open only a folder of Results/ whose id has the form the search form
    gives its jobs (pages_utils.results_job_dir: they take the id of their address into
    paths, and the region page into shell commands). A command-line search names its
    folder after its --output, which may hold other characters: History listed such a job
    with a link to pages that said there was no such job. It says why they do not open
    now, and what to do."""

    def test_a_job_named_otherwise_is_listed_with_why_its_pages_do_not_open(self):
        with tempfile.TemporaryDirectory() as tmp:
            ns = _load_history_table()
            ns["current_working_directory"] = tmp
            # the search form's ids, and an assembly search's, and three --output names
            jobs = ["web_G05B8KHU0H", "8CYCH5VMZ5_combined", "run.2-b", "run+2", "_test",
                    "my run"]
            for n, jobid in enumerate(jobs):
                job_dir = os.path.join(tmp, "Results", jobid)
                os.makedirs(job_dir)
                with open(os.path.join(job_dir, ".Params.txt"), "w") as f:
                    f.write("Genome_selected\thg38\nGenome_idx\tNRG_2_hg38\nPam\tNRG\n"
                            "Mismatches\t4\nDNA\t1\nRNA\t1\nMax_total_edits\t6\n")
                with open(os.path.join(job_dir, "log.txt"), "w") as f:
                    f.write(f"Job\tStart\tTue Nov 3 09:5{n}:00 UTC 2026\n")
                with open(os.path.join(job_dir, ".guides.txt"), "w") as f:
                    f.write("CTAACAGTTGCTTTTATCACNNN\n")
            summary = ns["construct_history_summary"](jobs)
            rows = ns["history_table_body_"](summary, 1, 1000000, summary.shape[0])
            linked, notes = {}, {}
            for row in rows:
                cell = row.children[0].children  # the Job cell
                if isinstance(cell, list):  # the name, a line break, why and what to do
                    name, _br, note = cell
                    notes[name] = note.children
                else:
                    linked[cell.children] = cell.kwargs["href"]
            self.assertEqual(linked, {j: f"load?job={j}" for j in jobs[:3]})
            self.assertEqual(sorted(notes), sorted(jobs[3:]))
            for name, note in notes.items():
                self.assertIn("letters, digits, '.', '_' and '-' only", note)
                self.assertIn(f"Results/{name}/", note)
                self.assertIn("--output", note)


if __name__ == "__main__":
    unittest.main()

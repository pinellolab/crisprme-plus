"""Guard tests for clean_intermediates.sh (CRISPRme+ 2.7.10).

The intermediate-cleanup step runs by DEFAULT at the end of every complete-search
and deletes ~tens of GB of working files on a genome-wide run. Because it is a
default-on destructive `rm`, the exact set it removes -- and, more importantly,
the set it must NEVER touch (the finalized integrated results, the report, the
alt-alignments .zip, the co-occurrence / high-variant-density sidecars, the
per-guide/per-sample summaries, and the sqlite .db web index) -- is pinned here.

Each test populates a throwaway output folder with the real file names produced
by submit_job_automated_new_multiple_vcfs.sh, invokes the actual shell script,
and asserts the survivors/removals. The two short-circuit guards
(CRISPRME_KEEP_INTERMEDIATES=True and the cicd_test argument) are exercised too,
since a regression in either silently turns cleanup on during the brute-force
gate or off for real users.
"""

import os
import subprocess
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "clean_intermediates.sh")

# A representative job. Folder basename == job name; the bestMerge/altMerge files
# are named "<job>.bestMerge.txt" / "<job>.altMerge.txt" (submit_job line ~1086).
JOB = "GUIDEA+NGG_hg38+1000G_6+2"
STEM = "sgRNA1+NGG_hg38+1000G_6+2"  # the per-guide finalized-output stem

# Files the cleanup MUST remove (the documented patterns).
REMOVED = [
    "{job}.altMerge.txt",
    "{stem}_all_results_with_alternative_alignments.tsv",
    "{stem}.population_summary.tsv",
    "sgRNA2+NGG_hg38+1000G_6+2.population_summary.tsv",  # prove the glob is plural
]

# The web per-position cluster drill-down (pages/results_page.py cluster_page) greps
# <job>.bestMerge.txt LIVE, and a web run cannot opt out of cleanup, so this file must
# survive -- deleting it leaves the drill-down table header-only. Pinned separately so a
# regression that re-adds it to the removal set fails loudly.
WEB_DRILLDOWN_SOURCE = "{job}.bestMerge.txt"

# Files the cleanup MUST keep. Note the adversarial pair: the alt-alignments .zip
# shares the exact stem of the removed .tsv and must survive the .tsv glob.
KEPT = [
    WEB_DRILLDOWN_SOURCE,
    "{stem}_integrated_results.tsv",
    "{stem}_integrated_results.zip",
    "{stem}_report.zip",
    "{stem}_all_results_with_alternative_alignments.zip",   # same stem as a REMOVED .tsv
    "{stem}.summary_by_guide.txt",
    "{stem}.summary_by_samples.txt",
    "{stem}.snp_snp_cooc.tsv",
    "{stem}.indel_snp_cooc.tsv",
    "{stem}.high_variant_density_regions.bed",
    "{stem}.indel_af.tsv",
    ".{job}.db",   # sqlite web results-page index (leading-dot hidden file)
    "{job}.db",
]


def _populate(folder):
    """Create every REMOVED + KEPT file (plus an imgs/ dir) in `folder`."""
    names = [n.format(job=JOB, stem=STEM) for n in (REMOVED + KEPT)]
    for name in names:
        with open(os.path.join(folder, name), "w") as fh:
            fh.write("x")
    imgs = os.path.join(folder, "imgs")
    os.makedirs(imgs, exist_ok=True)
    with open(os.path.join(imgs, "plot.png"), "w") as fh:
        fh.write("x")
    return names


def _run(folder, cicd_test="", keep_env=None):
    env = dict(os.environ)
    env.pop("CRISPRME_KEEP_INTERMEDIATES", None)
    if keep_env is not None:
        env["CRISPRME_KEEP_INTERMEDIATES"] = keep_env
    proc = subprocess.run(
        ["bash", SCRIPT, folder, cicd_test],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return proc


def _present(folder, name):
    return os.path.exists(os.path.join(folder, name.format(job=JOB, stem=STEM)))


class TestCleanIntermediates(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.folder = os.path.join(self._tmp.name, JOB)
        os.makedirs(self.folder)
        _populate(self.folder)

    def tearDown(self):
        self._tmp.cleanup()

    def test_default_removes_only_intermediates(self):
        """Default run (no env, no cicd flag): remove the four patterns, keep the rest."""
        proc = _run(self.folder)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in REMOVED:
            self.assertFalse(
                _present(self.folder, name),
                f"intermediate {name} should have been removed",
            )
        for name in KEPT:
            self.assertTrue(
                _present(self.folder, name),
                f"output {name} must NOT be removed by cleanup",
            )
        # imgs/ and its contents are untouched.
        self.assertTrue(os.path.isfile(os.path.join(self.folder, "imgs", "plot.png")))

    def test_web_cluster_drilldown_source_survives(self):
        """The file the web per-position drill-down greps live (<job>.bestMerge.txt)
        must survive a default cleanup -- a web search cannot opt out, so deleting it
        would leave the cluster table header-only. Guards against re-adding it to the
        removal set (the 2.7.10 critical regression)."""
        proc = _run(self.folder)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(
            _present(self.folder, WEB_DRILLDOWN_SOURCE),
            f"{WEB_DRILLDOWN_SOURCE} is the web cluster drill-down source and must "
            "NOT be removed by default cleanup",
        )

    def test_keep_env_short_circuits(self):
        """CRISPRME_KEEP_INTERMEDIATES=True keeps absolutely everything."""
        proc = _run(self.folder, keep_env="True")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in REMOVED + KEPT:
            self.assertTrue(
                _present(self.folder, name),
                f"{name} must survive when CRISPRME_KEEP_INTERMEDIATES=True",
            )

    def test_cicd_guard_short_circuits(self):
        """cicd_test=True (brute-force gate) keeps absolutely everything."""
        proc = _run(self.folder, cicd_test="True")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in REMOVED + KEPT:
            self.assertTrue(
                _present(self.folder, name),
                f"{name} must survive when cicd_test=True",
            )

    def test_keep_env_only_exact_true(self):
        """Only the literal 'True' keeps; any other value (e.g. 'False') still cleans."""
        proc = _run(self.folder, keep_env="False")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertFalse(_present(self.folder, "{job}.altMerge.txt"))
        self.assertTrue(_present(self.folder, "{stem}_integrated_results.tsv"))

    def test_empty_folder_is_safe(self):
        """A run over a folder with none of the targets must not error."""
        empty = os.path.join(self._tmp.name, "empty_job")
        os.makedirs(empty)
        proc = _run(empty)
        self.assertEqual(proc.returncode, 0, proc.stderr)


if __name__ == "__main__":
    unittest.main()

"""The web form's personal-assembly submit (pages/main_page.py
submit_assembly_search_job), called the way Dash calls it.

* With no individual picked, the warning lists one item, "select an
  individual", not the six hidden genome/chain/chromAlias fields the picker
  fills (they are not on the form).
* The job's .Params.txt names each haplotype by the bundle's assembly name --
  the results header and History showed "Paternal genome: genome | Maternal
  genome: genome", the name of the folder every bundle haplotype is kept in.
* The "Max edits" slider is passed on as --max-total-edits under the reference
  tab's rule (Advanced closed: the slider is the cap; open: the per-type caps
  govern). It was shown and ignored. A submission left at the defaults runs the
  same assembly-search command as before.

Needs a real (if minimal) Dash `app` import, same as the rest of
pages/main_page.py, so it runs with the web image's tests (web-e2e.yml), with
current_working_directory patched to a tmpdir and the background job launcher
replaced by a recorder.
"""

import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "PostProcess"))

import pages.main_page as mp  # noqa: E402

PAM = "20bp-NRG-SpCas9"
GUIDE = "CTAACAGTTGCTTTTATCAC"


def _make_bundle(cwd, individual):
    root = os.path.join(cwd, "Assemblies", individual)
    haps = {}
    for hap in ("paternal", "maternal"):
        gdir = os.path.join(root, hap, "genome")
        os.makedirs(gdir)
        with open(os.path.join(gdir, "chr1.fa"), "w") as fh:
            fh.write(">chr1\nACGT\n")
        for fname in (f"{individual}_{hap}_vs_GRCh38.chain", f"{individual}_{hap}.chromAlias.txt"):
            with open(os.path.join(root, hap, fname), "w") as fh:
                fh.write("x\n")
        haps[hap] = {"assembly_name": f"{individual}_{hap}_CM0{len(haps)}.1",
                     "genome_dir": f"{hap}/genome",
                     "chain": f"{hap}/{individual}_{hap}_vs_GRCh38.chain",
                     "chromalias": f"{hap}/{individual}_{hap}.chromAlias.txt"}
    with open(os.path.join(root, "metadata.json"), "w") as fh:
        json.dump({"format": "crisprme-personal-assembly", "format_version": 1,
                   "individual": individual, "haplotypes": haps}, fh)


class _SubmitHarness(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)
        os.makedirs(os.path.join(self.cwd, "Results"))
        os.makedirs(os.path.join(self.cwd, "PAMs"))
        shutil.copy(os.path.join(_REPO_ROOT, "PAMs", f"{PAM}.txt"), os.path.join(self.cwd, "PAMs"))
        _make_bundle(self.cwd, "WALK01")
        self.jobs = []
        for target, value in (
            ("current_working_directory", self.cwd + "/"),
            ("build_active_annotation", lambda genome: ("vuoto.txt", "vuoto.txt")),
            ("pool_executor", type("Pool", (), {"submit": lambda _s, fn, *a: self.jobs.append(a)})()),
        ):
            patcher = patch.object(mp, target, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def fields(self, individual):
        """The six hidden fields as populate_assembly_fields_from_individual fills
        them from the picker: cwd-relative paths, or all None with none picked."""
        if not individual:
            return [None] * 6
        rel = lambda hap, name: os.path.join("Assemblies", individual, hap, name)  # noqa: E731
        return [rel("paternal", "genome"), rel("maternal", "genome"),
                rel("paternal", f"{individual}_paternal_vs_GRCh38.chain"),
                rel("maternal", f"{individual}_maternal_vs_GRCh38.chain"),
                rel("paternal", f"{individual}_paternal.chromAlias.txt"),
                rel("maternal", f"{individual}_maternal.chromAlias.txt")]

    def submit(self, individual="WALK01", guides=GUIDE, mms=4, dna=1, rna=1,
               max_edits=6, advanced_open=False):
        # the form's defaults: 4 mismatches, 1 DNA + 1 RNA bulge, slider 6, Advanced closed
        return mp.submit_assembly_search_job(
            1, *self.fields(individual), PAM, guides, mms, dna, rna, None, individual,
            max_edits, advanced_open)

    def command(self, **kwargs):
        """The assembly-search command line the submit hands to the job runner."""
        self.assertIsNone(self.warnings(self.submit(**kwargs)))
        cmd, combined_dir, note = self.jobs[-1]
        return cmd.split()

    @staticmethod
    def warnings(result):
        """The items of the 'Missing or wrong input' list, or None when the job
        was submitted."""
        pathname, search, is_open, children = result
        if not is_open:
            return None
        return [li.children for li in children.children]


class NoIndividualPicked(_SubmitHarness):
    def test_one_item_says_what_to_do(self):
        self.assertEqual(self.warnings(self.submit(individual=None)),
                         ["select an individual on the Personal assembly tab"])
        self.assertEqual(self.jobs, [])

    def test_other_missing_input_is_still_listed(self):
        self.assertEqual(self.warnings(self.submit(individual=None, guides="")),
                         ["select an individual on the Personal assembly tab",
                          "at least one guide sequence"])

    def test_a_picked_individual_submits(self):
        self.assertIsNone(self.warnings(self.submit()))
        self.assertEqual(len(self.jobs), 1)


class HaplotypeLabels(_SubmitHarness):
    def test_params_name_each_haplotype_by_its_assembly_name(self):
        pathname, search, _, _ = self.submit()
        job = search.split("=", 1)[1]
        params = {}
        with open(os.path.join(self.cwd, "Results", job, ".Params.txt")) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                params[fields[1]] = fields[2]
        self.assertEqual(params["Genome_paternal"], "WALK01_paternal_CM00.1")
        self.assertEqual(params["Genome_maternal"], "WALK01_maternal_CM01.1")


class MaxEditsSlider(_SubmitHarness):
    def cap(self, **kwargs):
        argv = self.command(**kwargs)
        return argv[argv.index("--max-total-edits") + 1] if "--max-total-edits" in argv else None

    def test_the_slider_caps_the_total_edits(self):
        self.assertEqual(self.cap(max_edits=3), "3")

    def test_the_defaults_pass_no_cap(self):
        # assembly-search's own default is mm + bDNA + bRNA = 6: the same search,
        # by the same command, as before the slider was read
        self.assertIsNone(self.cap())
        argv = self.command()
        self.assertEqual(argv[argv.index("--mm"):argv.index("--mm") + 8],
                         ["--mm", "4", "--bDNA", "1", "--bRNA", "1", "--merge", "3"])

    def test_with_advanced_open_the_per_type_caps_govern(self):
        self.assertIsNone(self.cap(max_edits=3, advanced_open=True))

    def test_a_slider_at_or_above_the_per_type_total_adds_nothing(self):
        self.assertIsNone(self.cap(max_edits=6, mms=4, dna=1, rna=1))
        self.assertIsNone(self.cap(max_edits=8, mms=4, dna=1, rna=1))
        self.assertEqual(self.cap(max_edits=6, mms=6, dna=1, rna=1), "6")


if __name__ == "__main__":
    unittest.main()

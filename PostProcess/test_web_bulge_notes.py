"""The bulge note and defaults of the search form match what is installed and what a
search does.

What a web search does (submit_job_automated_new_multiple_vcfs.sh, launched without
--index-path): with no installed reference index deep enough for the bulges asked
for, it builds one from the installed genome first, once. A variant panel's index
is built the same way, at the bulges asked for, but the form caps a panel's bulges
at its installed index: with none yet only a 0-bulge search is offered, and that
search first prepares the variant genome and a 0-bulge index. A deeper one is built
in Settings.

The note said "the first bulge search will build the reference index" whenever the
reference depth relied on the raw genome -- also with an installed index that serves
the bulges selected -- and "no raw genome to build one" for a panel whose genome and
VCF files are installed. The cap set DNA/RNA to 0 (and so Max edits to 4), and lifting
it did not give them back.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_bulge_notes.py' -v
"""

import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
    from dash import no_update
except ImportError:  # pragma: no cover
    raise unittest.SkipTest("dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

import pages.main_page as mp  # noqa: E402
import pages.pages_utils as pu  # noqa: E402

PAM = "20bp-NRG-SpCas9"


class _Install(unittest.TestCase):
    """A data folder: genome hg38 (raw FASTA), panels "walk" (VCFs + a 1-bulge
    index), "fresh" (VCFs, no index), "zero" (VCFs + a 0-bulge index) and "shipped"
    (a 0-bulge index, no VCF files)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        for d in ("PAMs", "Genomes/hg38", "genome_library", "VCFs/hg38_walk",
                  "VCFs/hg38_fresh", "VCFs/hg38_zero"):
            os.makedirs(os.path.join(self.tmp, d))
        with open(os.path.join(self.tmp, "PAMs", f"{PAM}.txt"), "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        open(os.path.join(self.tmp, "Genomes", "hg38", "chr1.fa"), "w").close()
        for panel in ("walk", "fresh", "zero"):
            open(os.path.join(self.tmp, "VCFs", f"hg38_{panel}", "x.chr1.vcf.gz"), "w").close()
        for d in ("NRG_2_hg38+hg38_walk", "NRG_1_hg38+hg38_zero", "NRG_1_hg38+hg38_shipped"):
            self.index(d)
        for mod in (mp, pu):
            p = patch.object(mod, "current_working_directory", self.tmp + os.sep)
            p.start()
            self.addCleanup(p.stop)

    def index(self, name):
        os.makedirs(os.path.join(self.tmp, "genome_library", name))


class TestBulgeNote(_Install):
    def note(self, variant, dna, rna):
        return mp.update_bulge_note("hg38", PAM, variant, dna, rna)

    def test_installed_reference_index_serves_the_bulges_selected(self):
        self.index("NRG_2_hg38")
        for variant in ("ref", "walk"):
            with self.subTest(variant=variant):
                note = self.note(variant, 1, 1)
                self.assertNotIn("build", note)
                self.assertTrue(note.startswith("Up to "))

    def test_more_bulges_than_the_installed_reference_index(self):
        self.index("NRG_2_hg38")
        note = self.note("ref", 2, 1)
        self.assertIn("Up to 2 DNA/RNA bulge(s) available.", note)
        self.assertIn("The installed reference index supports 1 bulge(s), so a search with 2 "
                      "first builds a deeper one from the installed genome", note)

    def test_no_reference_index_at_all(self):
        for dna in (0, 1):  # a 0-bulge search builds one too
            with self.subTest(dna=dna):
                self.assertIn("No reference index for this PAM is installed, so the search "
                              "first builds one", self.note("ref", dna, 0))

    def test_panel_without_an_index_but_with_its_vcf_files(self):
        self.index("NRG_2_hg38")
        note = self.note("fresh", 0, 0)
        self.assertIn("No fresh variant index for this PAM is installed yet, so only a "
                      "0-bulge search is offered, and the first one prepares", note)
        self.assertIn("build its index in Settings (Build an index locally) first", note)
        self.assertNotIn("raw genome", note)

    def test_panel_with_a_0_bulge_index(self):
        self.index("NRG_2_hg38")
        note = self.note("zero", 0, 0)
        self.assertIn("The installed zero variant index for this PAM supports no bulges", note)
        self.assertIn("build its index in Settings", note)
        note = self.note("shipped", 0, 0)  # no VCF files to build one from
        self.assertIn("download a variant index for it in Settings", note)

    def test_no_index_and_no_genome_files(self):
        os.remove(os.path.join(self.tmp, "Genomes", "hg38", "chr1.fa"))
        self.assertIn("no FASTA files", self.note("ref", 0, 0))


class TestBulgeCap(_Install):
    def limit(self, variant, dna, rna, edits, held, mms=4):
        return mp.limit_bulges_to_index("hg38", PAM, variant, dna, rna, mms, edits, held)

    def test_values_the_cap_lowered_come_back_when_it_is_lifted(self):
        self.index("NRG_2_hg38")
        # a panel with no index caps the bulges at 0 (and sync_max_edits_bound then
        # pulls Max edits from 6 down to 4 + 0 + 0)
        opts, _, dna, rna, edits, held = self.limit("fresh", 1, 1, 6, {})
        self.assertEqual((dna, rna, [o["value"] for o in opts]), (0, 0, [0]))
        self.assertIs(edits, no_update)
        new_max, _marks, slider = mp.sync_max_edits_bound(4, dna, rna, 6)
        self.assertEqual((new_max, slider), (4, 4))
        # back to "Reference only": 1/1 and 6 again
        opts, _, dna, rna, edits, held = self.limit("ref", dna, rna, slider, held)
        self.assertEqual((dna, rna, edits, held), (1, 1, 6, {}))
        self.assertEqual([o["value"] for o in opts], [0, 1, 2])
        self.assertEqual(mp.sync_max_edits_bound(4, dna, rna, edits)[2], 6)

    def test_a_value_the_user_changed_since_stays(self):
        self.index("NRG_2_hg38")
        *_, dna, rna, edits, held = self.limit("fresh", 1, 1, 6, {})
        # (on the capped panel the user cannot pick another bulge count, but can
        # move the slider:) Max edits moved to 3 -- it stays at 3
        *_, dna2, rna2, edits2, _ = self.limit("ref", 0, 0, 3, held)
        self.assertEqual((dna2, rna2, edits2), (1, 1, no_update))
        # DNA set to 0 by the user on a 1-bulge panel: it is not raised again
        *_, dna, rna, edits, held = self.limit("walk", 2, 2, 6, {})
        self.assertEqual((dna, rna), (1, 1))
        *_, dna, rna, edits, held = self.limit("ref", 0, 1, 6, held)
        self.assertEqual((dna, rna, held), (0, 1, {}))

    def test_partly_lifted_cap(self):
        self.index("NRG_3_hg38")
        *_, dna, rna, _e, held = self.limit("fresh", 2, 2, 8, {})
        self.assertEqual((dna, rna), (0, 0))
        *_, dna, rna, _e, held = self.limit("walk", dna, rna, 4, held)  # walk: 1 bulge
        self.assertEqual((dna, rna), (1, 1))
        self.assertEqual(held["dna"], 2)  # still owed
        slider = mp.sync_max_edits_bound(4, dna, rna, 4)[2]  # 4 (it was pulled down to 4)
        *_, dna, rna, edits, held = self.limit("ref", dna, rna, slider, held)
        self.assertEqual((dna, rna, edits, held), (2, 2, 8, {}))

    def test_nothing_lowered_nothing_changed(self):
        self.index("NRG_2_hg38")
        opts, _, dna, rna, edits, held = self.limit("walk", 1, 0, 6, {})
        self.assertEqual((dna, rna, edits, held), (1, 0, no_update, {}))


if __name__ == "__main__":
    unittest.main()

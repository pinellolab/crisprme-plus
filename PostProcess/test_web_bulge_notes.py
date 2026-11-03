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
it did not give them back. They come back also when the browser runs the cap twice for
one change, whatever the order of the two runs and their answers (_Form), and Max edits
comes back as far as the bulges allow when a cap is lifted only partly (a panel with a
1-bulge index after one with none): it stayed at the 4 of the 0-bulge cap.

A sites-only dataset was offered a search that cannot run (TestRefusedDatasets): Submit
refuses every search of it, with or without its index, and the note offered a 0-bulge
search of it and a build in Settings, or the bulges of its index. Now it says what
Submit says.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_bulge_notes.py' -v
"""

import gzip
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

    def test_panel_with_no_index_for_this_pam_and_nothing_to_build_one(self):
        # a wider PAM-less index and no VCF files: the search script refuses it, so the note
        # does not promise that a first search prepares anything
        self.index("NRG_2_hg38")
        self.index("NNNNNNNN_3_hg38+hg38_k8")
        note = self.note("k8", 0, 0)
        self.assertNotIn("prepares", note)
        self.assertIn("neither the panel's VCF files nor its enriched genome are here", note)
        self.assertIn("download a variant index for it in Settings", note)

    def test_panel_with_its_enriched_genome_can_still_build_one(self):
        self.index("NRG_2_hg38")
        self.index("NNNNNNNN_3_hg38+hg38_k8")
        os.makedirs(os.path.join(self.tmp, "Genomes", "hg38+hg38_k8"))
        open(os.path.join(self.tmp, "Genomes", "hg38+hg38_k8", "chr1.enriched.fa"), "w").close()
        self.assertIn("the first one prepares", self.note("k8", 0, 0))

    def test_no_index_and_no_genome_files(self):
        os.remove(os.path.join(self.tmp, "Genomes", "hg38", "chr1.fa"))
        self.assertIn("no FASTA files", self.note("ref", 0, 0))


def _vcf(path, samples=("S1", "S2")):
    """A VCF file with a real header, gzipped as bgzip writes it: with sample columns,
    or with none (a sites-only file, as the hg38_mega VCFs)."""
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        cols += "\tFORMAT\t" + "\t".join(samples)
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n" + cols + "\n")


class TestRefusedDatasets(_Install):
    """A dataset whose search Submit refuses whatever the bulges is not offered as
    searchable: the note says what Submit says."""

    def setUp(self):
        super().setUp()
        self.index("NRG_2_hg38")

    def note(self, variant, dna=0, rna=0):
        return mp.update_bulge_note("hg38", PAM, variant, dna, rna)

    def submit(self, variant, dna=0, rna=0):
        """check_input, as the Submit button runs it: (runs, messages)."""
        out = mp.check_input(1, None, "hg38", PAM, "IP", "CTAACAGTTGCTTTTATCAC", 4, dna, rna,
                             variant, False)
        return out[0] == 1, [li.children for li in out[-1].children[1].children]

    def dataset(self, name, files, samples=("S1", "S2"), sample_list=True):
        """A dataset as Settings adds it: its VCF files and, when they name samples, the
        sample list written from their header."""
        folder = os.path.join(self.tmp, "VCFs", f"hg38_{name}")
        os.makedirs(folder)
        for f in files:
            _vcf(os.path.join(folder, f), samples)
        if sample_list:
            os.makedirs(os.path.join(self.tmp, "samplesIDs"), exist_ok=True)
            with open(os.path.join(self.tmp, "samplesIDs", f"hg38_{name}.samplesID.txt"),
                      "w") as fh:
                fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                         + "".join(f"{s}\thg38_{name}\tunknown\tunknown\n" for s in samples))
        return folder

    def test_sites_only_dataset_with_or_without_its_index(self):
        # Submit refuses a search of a sites-only dataset, which has no sample list, and
        # building its index does not change that; the note offered a 0-bulge search and a
        # build in Settings, then the bulges of the index built
        self.dataset("sites", ["sites.chr1.vcf.gz"], samples=(), sample_list=False)
        for built in (False, True):
            if built:
                self.index("NRG_2_hg38+hg38_sites")
            with self.subTest(built=built):
                note = self.note("sites")
                self.assertEqual(note, pu.samplesid_problem("hg38_sites", "sites"))
                self.assertTrue(note.startswith("The variant dataset 'sites' is sites-only"))
                self.assertNotIn("Build an index locally", note)
                self.assertEqual(self.submit("sites"), (False, [note]))

    def test_datasets_submit_does_not_refuse_keep_their_note(self):
        # a sites-only dataset with a sample list, as the published hg38_mega index
        # installs one: searchable, so the note is the bulge note
        self.dataset("mega", ["mega.chr1.vcf.gz"], samples=(), sample_list=False)
        os.makedirs(os.path.join(self.tmp, "samplesIDs"), exist_ok=True)
        with open(os.path.join(self.tmp, "samplesIDs", "hg38_mega.samplesID.txt"), "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                     "gnomAD\tgnomAD\tgnomAD\tunknown\n")
        self.index("NRG_2_hg38+hg38_mega")
        self.assertEqual(self.note("mega", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        # a genotyped dataset with no sample list: Submit refuses it, saying how to get
        # one (a build in Settings writes it); the note is the bulge note
        self.dataset("nolist", ["nolist.chr1.vcf.gz"], sample_list=False)
        self.assertIn("the first one prepares", self.note("nolist"))
        ok, msgs = self.submit("nolist")
        self.assertFalse(ok)
        self.assertEqual(len(msgs), 1)
        self.assertTrue(msgs[0].startswith("No sample list for the variant dataset 'nolist'"))


class _Form:
    """The form values the browser holds. A callback's answer is applied as Dash does
    it (no_update leaves a value as it is), and a request carries the values held
    when it is sent, which need not be the latest: the browser runs
    limit_bulges_to_index twice for one change of the variant dataset (for the
    dataset, then for the PAM that update_pam_options writes back), and the second
    request can go out before or after the first answer, and before or after
    sync_max_edits_bound has answered."""

    def __init__(self, dna=1, rna=1, mms=4, edits=6):
        self.v = {"dna": dna, "rna": rna, "mms": mms, "edits": edits, "held": {}}

    def state(self, **override):
        return dict(self.v, **override)

    @staticmethod
    def limit(variant, s):
        return mp.limit_bulges_to_index(
            "hg38", PAM, variant, s["dna"], s["rna"], s["mms"], s["edits"], s["held"])

    def apply_limit(self, answer):
        _, _, dna, rna, edits, held = answer
        for key, value in (("dna", dna), ("rna", rna), ("edits", edits), ("held", held)):
            if value is not no_update:
                self.v[key] = value

    def sync(self):
        s = self.v
        self.v["edits"] = mp.sync_max_edits_bound(s["mms"], s["dna"], s["rna"], s["edits"])[2]

    def change(self, variant, order):
        """One change of the variant dataset, with its two runs of the cap in `order`:
        "before" (the second request goes out before the first answer lands),
        "unlowered" (after it, with Max edits still where it was before the change,
        as when only sync_max_edits_bound lowered it), "after" (after it) or "synced"
        (after sync_max_edits_bound answered too)."""
        first = self.limit(variant, self.state())
        if order == "before":
            second = self.limit(variant, self.state())
            self.apply_limit(first)
            self.apply_limit(second)
        else:
            edits_before = self.v["edits"]
            self.apply_limit(first)
            if order == "synced":
                self.sync()
            override = {"edits": edits_before} if order == "unlowered" else {}
            self.apply_limit(self.limit(variant, self.state(**override)))
        self.sync()
        return self.v["dna"], self.v["rna"], self.v["edits"]


ORDERS = ("before", "unlowered", "after", "synced")


class TestBulgeCap(_Install):
    def limit(self, variant, dna, rna, edits, held, mms=4):
        return mp.limit_bulges_to_index("hg38", PAM, variant, dna, rna, mms, edits, held)

    def test_values_the_cap_lowered_come_back_when_it_is_lifted(self):
        self.index("NRG_2_hg38")
        # a panel with no index caps the bulges at 0, and Max edits from 6 down to
        # 4 + 0 + 0 (sync_max_edits_bound then sets it to the same)
        opts, _, dna, rna, edits, held = self.limit("fresh", 1, 1, 6, {})
        self.assertEqual((dna, rna, [o["value"] for o in opts]), (0, 0, [0]))
        self.assertEqual(edits, 4)
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
        *_, dna, rna, edits, held = self.limit("fresh", 2, 2, 8, {})
        self.assertEqual((dna, rna, edits), (0, 0, 4))
        *_, dna, rna, edits, held = self.limit("walk", dna, rna, edits, held)  # walk: 1 bulge
        # 1/1, and Max edits up to 4 + 1 + 1 = 6 again, not left at the 4 the first
        # cap set; DNA/RNA 2 and Max edits 8 are still owed
        self.assertEqual((dna, rna, edits), (1, 1, 6))
        self.assertEqual((held["dna"], held["rna"], held["edits"]), (2, 2, 8))
        slider = mp.sync_max_edits_bound(4, dna, rna, edits)[2]
        self.assertEqual(slider, 6)
        *_, dna, rna, edits, held = self.limit("ref", dna, rna, slider, held)
        self.assertEqual((dna, rna, edits, held), (2, 2, 8, {}))

    def test_nothing_lowered_nothing_changed(self):
        self.index("NRG_2_hg38")
        opts, _, dna, rna, edits, held = self.limit("walk", 1, 0, 6, {})
        self.assertEqual((dna, rna, edits, held), (1, 0, no_update, {}))

    def test_two_runs_for_one_change_keep_what_the_cap_owes(self):
        # The second run of the cap for one change, sent after the first answer but
        # with Max edits not yet lowered, took Max edits 6 for a value the user had
        # set, and forgot it: Reference only then gave back 1/1 but left Max edits at
        # 4, and a search with Advanced closed ran with at most 4 mismatches + bulges.
        self.index("NRG_2_hg38")
        for capping in ORDERS:
            for lifting in ORDERS:
                with self.subTest(capping=capping, lifting=lifting):
                    form = _Form()
                    self.assertEqual(form.change("fresh", capping), (0, 0, 4))
                    self.assertEqual(form.change("ref", lifting), (1, 1, 6))

    def test_two_runs_for_each_of_several_changes(self):
        # Reference only with DNA 2, then walk (1 bulge), then panels with no
        # bulge index (0), then Reference only: 2/1 and Max edits 6 again
        self.index("NRG_2_hg38")
        for order in ORDERS:
            with self.subTest(order=order):
                form = _Form(dna=2)
                self.assertEqual(form.change("walk", order), (1, 1, 6))
                self.assertEqual(form.change("fresh", order), (0, 0, 4))
                self.assertEqual(form.change("zero", order), (0, 0, 4))
                self.assertEqual(form.change("ref", order), (2, 1, 6))

    def test_max_edits_moved_by_the_user_stays_with_two_runs(self):
        self.index("NRG_2_hg38")
        for order in ORDERS:
            with self.subTest(order=order):
                form = _Form()
                form.change("fresh", order)
                form.v["edits"] = 3  # moved while the bulges are capped
                form.change("zero", order)  # another panel capped at 0
                self.assertEqual(form.change("ref", order), (1, 1, 3))

    def test_a_partly_lifted_cap_gives_back_the_max_edits_it_allows(self):
        # Reference only with DNA 2, then walk2 (no index: 0 bulges), then walk (a
        # 1-bulge index): the cap gave back 1/1 but left Max edits at the 4 the 0-bulge
        # cap had set, while walk straight from Reference only gave 6, so the same form
        # ran a search with at most 4 or 6 mismatches + bulges (Advanced closed),
        # depending on the way there. Both ways give 6 now, with the two runs of the
        # cap for each change in any order.
        self.index("NRG_2_hg38")
        for capping in ORDERS:
            for lifting in ORDERS:
                with self.subTest(capping=capping, lifting=lifting):
                    direct = _Form(dna=2)
                    self.assertEqual(direct.change("walk", lifting), (1, 1, 6))
                    form = _Form(dna=2)
                    self.assertEqual(form.change("fresh", capping), (0, 0, 4))
                    self.assertEqual(form.change("walk", lifting), (1, 1, 6))
                    self.assertEqual(form.change("ref", lifting), (2, 1, 6))
                    form = _Form(dna=2)  # walk, a panel with no index, walk again
                    self.assertEqual(form.change("walk", capping), (1, 1, 6))
                    self.assertEqual(form.change("fresh", capping), (0, 0, 4))
                    self.assertEqual(form.change("walk", lifting), (1, 1, 6))
                    self.assertEqual(form.change("ref", lifting), (2, 1, 6))

    def test_a_partly_lifted_cap_keeps_a_lower_max_edits_and_one_moved_by_hand(self):
        self.index("NRG_2_hg38")
        for order in ORDERS:
            with self.subTest(order=order):
                # Max edits 5, below the 4 + 2 + 1 the bulges allow: 4 on the 0-bulge
                # panel, 5 again on walk (4 + 1 + 1 = 6 allows it) and on Reference only
                form = _Form(dna=2, edits=5)
                self.assertEqual(form.change("fresh", order), (0, 0, 4))
                self.assertEqual(form.change("walk", order), (1, 1, 5))
                self.assertEqual(form.change("ref", order), (2, 1, 5))
                # moved by hand on the 0-bulge panel: it stays where the user put it
                form = _Form(dna=2)
                form.change("fresh", order)
                form.v["edits"] = 3
                self.assertEqual(form.change("walk", order), (1, 1, 3))
                self.assertEqual(form.change("ref", order), (2, 1, 3))


if __name__ == "__main__":
    unittest.main()

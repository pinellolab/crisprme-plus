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

Two kinds of dataset were offered a search that cannot run (TestRefusedDatasets): one
added as a single multi-chromosome VCF file, with no variant index yet, whose search
fails as it prepares the variant genome (only a build of its index splits the file per
chromosome), and a sites-only one, whose search Submit refuses with or without its
index. The note offered a 0-bulge search of either and a build in Settings, or the
bulges of its index; now it says what Submit says, and Submit refuses the first too.
So it does of a dataset with no sample list, which Submit refused while the note offered
its search: one whose VCFs have no list yet, and one whose index stays after its VCF
dataset was removed in Settings, which takes the list with it. For that one Submit said
to download its index again, which a locally built index cannot be, and which Settings
does not do for an index that is installed: it now says to remove a downloaded one first.

Submit itself failed with a server error, and the page did nothing, when the genome or
PAM selected was not installed (TestSubmitWithoutAGenomeOrPAM), as after Load Example on
the empty data folder the web app starts with: it said what is missing only once both
were installed. Now it says so, and where to add them.

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_bulge_notes.py' -v
"""

import gzip
import json
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
import pages.settings_page as sp  # noqa: E402
import index  # noqa: E402  (check_directories: the data folder the web app starts with)

PAM = "20bp-NRG-SpCas9"


class _Install(unittest.TestCase):
    """A data folder: genome hg38 (raw FASTA), panels "walk" (VCFs + a 1-bulge
    index), "fresh" (VCFs, no index), "zero" (VCFs + a 0-bulge index) and "shipped"
    (a 0-bulge index, no VCF files), each with its sample list."""

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
        for panel in ("walk", "fresh", "zero", "shipped"):
            self.sample_list(panel)
        for mod in (mp, pu, sp):
            p = patch.object(mod, "current_working_directory", self.tmp + os.sep)
            p.start()
            self.addCleanup(p.stop)

    def index(self, name):
        os.makedirs(os.path.join(self.tmp, "genome_library", name))

    def sample_list(self, panel, samples=("S1", "S2")):
        """The sample list a search of panel reads, as a download or Settings installs it:
        without it Submit refuses the search, and the note says so."""
        os.makedirs(os.path.join(self.tmp, "samplesIDs"), exist_ok=True)
        with open(os.path.join(self.tmp, "samplesIDs", f"hg38_{panel}.samplesID.txt"), "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
                     + "".join(f"{s}\thg38_{panel}\tunknown\tunknown\n" for s in samples))


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
        self.sample_list("k8")
        note = self.note("k8", 0, 0)
        self.assertNotIn("prepares", note)
        self.assertIn("neither the panel's VCF files nor its enriched genome are here", note)
        self.assertIn("download a variant index for it in Settings", note)

    def test_panel_with_its_enriched_genome_can_still_build_one(self):
        self.index("NRG_2_hg38")
        self.index("NNNNNNNN_3_hg38+hg38_k8")
        self.sample_list("k8")
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

    SINGLE = ("VCF is a single multi-chromosome file, which a search cannot prepare the "
              "variant genome from: the file is split per chromosome only when its index is "
              "built. Build its index in Settings (Build an index locally) first, or remove "
              "the dataset in Settings and add it again as one VCF file per chromosome "
              "(<name>.chrN.vcf.gz).")

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

    def test_single_multi_chromosome_vcf_without_its_index(self):
        # the search fails as it begins (CRISPRitz add-variants takes the chromosome from a
        # '.chrN.' token in the file name: "cannot open Genomes/hg38/.fa"), and the note
        # said that the first 0-bulge search prepares the variant genome
        folder = self.dataset("one", ["walk.vcf.gz"])
        note = self.note("one")
        self.assertEqual(note, "No one variant index for this PAM is installed yet, and the "
                               "dataset's " + self.SINGLE)
        self.assertNotIn("prepares the variant genome and", note)
        self.assertEqual(self.submit("one"), (False, [note]))
        # with Reference only, or with the dataset selected with another one
        self.assertEqual(self.submit("ref"), (True, []))
        self.dataset("per", ["walk.chr1.vcf.gz", "walk.chr2.vcf.gz"])
        self.assertEqual(self.note("per+one"), note)
        self.assertEqual(self.submit("per+one"), (False, [note]))
        # one file per chromosome: offered as before, and Submit runs it
        self.assertIn("No per variant index for this PAM is installed yet, so only a 0-bulge "
                      "search is offered, and the first one prepares", self.note("per"))
        self.assertEqual(self.submit("per"), (True, []))
        # its variant genome already prepared: a search does not read its VCF
        os.makedirs(os.path.join(self.tmp, "Genomes", "hg38+hg38_one"))
        self.assertIn("the first one prepares", self.note("one"))
        self.assertEqual(self.submit("one"), (True, []))
        os.rmdir(os.path.join(self.tmp, "Genomes", "hg38+hg38_one"))
        # built in Settings: the file split per chromosome, and its index there
        os.rename(os.path.join(folder, "walk.vcf.gz"), os.path.join(folder, "walk.chr1.vcf.gz"))
        _vcf(os.path.join(folder, "walk.chr2.vcf.gz"))
        self.index("NRG_2_hg38+hg38_one")
        self.assertEqual(self.note("one", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        self.assertEqual(self.submit("one", 1, 1), (True, []))
        # a single-file dataset with a variant index of its own for this PAM (as a
        # downloaded one, whose dictionaries a search reads instead of the VCFs)
        self.dataset("dl", ["walk.vcf.gz"])
        self.index("NRG_1_hg38+hg38_dl")
        self.assertIn("The installed dl variant index for this PAM supports no bulges",
                      self.note("dl"))
        self.assertEqual(self.submit("dl"), (True, []))

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
        # a single multi-chromosome sites-only file: the sites-only message alone, as a
        # build does not make it searchable
        self.dataset("agg", ["agg.vcf.gz"], samples=(), sample_list=False)
        note = self.note("agg")
        self.assertTrue(note.startswith("The variant dataset 'agg' is sites-only"))
        self.assertEqual(self.submit("agg"), (False, [note]))

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
        # a merged panel, as the combined 1000G + HGDP index, with the lists of its datasets
        # and none of its own: Submit makes it from theirs, so the note does not refuse it
        self.sample_list("1000G")
        self.sample_list("HGDP")
        self.index("NRG_2_hg38+hg38_1000G_HGDP")
        self.assertEqual(self.note("1000G_HGDP", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        self.assertEqual(self.submit("1000G_HGDP", 1, 1), (True, []))

    def test_a_dataset_without_its_sample_list(self):
        # a genotyped dataset with no sample list: Submit refuses it, saying how to get one
        # (a build in Settings writes it), and the note offered its search
        self.dataset("nolist", ["nolist.chr1.vcf.gz"], sample_list=False)
        note = self.note("nolist")
        self.assertEqual(note, pu.samplesid_problem("hg38_nolist", "nolist"))
        self.assertTrue(note.startswith("No sample list for the variant dataset 'nolist'"))
        self.assertIn("Build its variant index in Settings, which writes the list", note)
        self.assertEqual(self.submit("nolist"), (False, [note]))
        # added as a single multi-chromosome file too: the note says both things Submit says
        self.dataset("nolistone", ["nolistone.vcf.gz"], sample_list=False)
        ok, msgs = self.submit("nolistone")
        self.assertEqual((ok, len(msgs)), (False, 2))
        self.assertEqual(self.note("nolistone"), " ".join(msgs))

    def test_a_dataset_removed_in_settings_whose_index_stays(self):
        # Settings: per-chromosome VCFs uploaded as hg38_walkint, its variant index built,
        # then 'VCF dataset: hg38_walkint' removed, which takes the sample list the page
        # wrote and leaves the index. The search form offered the dataset as before, the
        # note said nothing, and Submit said to download its index again, which a local
        # build cannot be
        def upload():
            part = os.path.join(self.tmp, "walkint.chr1.vcf.gz.part")
            _vcf(part)
            sp._finalize_upload("vcf", part, "walkint.chr1.vcf.gz", "hg38_walkint", "hg38")

        upload()
        sid = os.path.join(self.tmp, "samplesIDs", "hg38_walkint.samplesID.txt")
        self.assertTrue(os.path.isfile(sid))
        # its build: the index, and the tiers that say what it holds
        self.index("NRG_2_hg38+hg38_walkint")
        os.makedirs(os.path.join(self.tmp, "Dictionaries", "registry_hg38_walkint"))
        with open(os.path.join(self.tmp, "Dictionaries", "registry_hg38_walkint",
                               "reg_chr1.idx"), "w") as fh:
            json.dump({"data_type": "genotyped-phased"}, fh)

        def label():
            return {o["value"]: o["label"] for o in pu.get_variant_dataset_options("hg38")}

        built = label()["walkint"]
        self.assertTrue(built.endswith("[genotyped, phased]"), built)
        self.assertEqual(self.note("walkint", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        self.assertEqual(self.submit("walkint", 1, 1), (True, []))
        self.assertIn("searches using them will break until you re-add it",
                      sp._delete_summary("vcf:hg38_walkint"))
        sp.do_delete(1, "vcf:hg38_walkint")
        self.assertFalse(os.path.exists(sid))
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, "genome_library",
                                                   "NRG_2_hg38+hg38_walkint")))
        # offered as before, with the same label; the note says what Submit says, and what
        # to do fits an index built here as well as a downloaded one
        self.assertEqual(label()["walkint"], built)
        note = self.note("walkint", 1, 1)
        self.assertEqual(note, pu.samplesid_problem("hg38_walkint", "walkint"))
        self.assertEqual(note, (
            "No sample list for the variant dataset 'walkint': a search with it reads "
            "samplesIDs/hg38_walkint.samplesID.txt (one row per sample its VCFs genotype) and "
            "fails at its end without it. If its index came from HuggingFace, remove the "
            "index in Settings and download it again there, which installs the list; add the "
            "VCF dataset hg38_walkint again in Settings, which writes the list from the VCF "
            "header; or add the file yourself, then submit again."))
        self.assertEqual(self.submit("walkint", 1, 1), (False, [note]))
        # added again in Settings: its list is written, and its index serves the search
        upload()
        self.assertTrue(os.path.isfile(sid))
        self.assertEqual(self.note("walkint", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        self.assertEqual(self.submit("walkint", 1, 1), (True, []))

    def test_a_dataset_named_otherwise_removed_in_settings(self):
        # the same with the dataset named walkint, not hg38_walkint: a Settings build names
        # its index after the folder (NRG_2_hg38+walkint), and the search form reads a
        # dataset off an index name only as a download names it (+hg38_<dataset>), so it
        # offers the dataset no more once its folder is gone; added again under the same
        # name, the dataset is offered and its index serves the search
        # (docs/SETTINGS_DATA_MANAGER.md says both)
        def upload():
            part = os.path.join(self.tmp, "walkint.chr1.vcf.gz.part")
            _vcf(part)
            sp._finalize_upload("vcf", part, "walkint.chr1.vcf.gz", "walkint", "hg38")

        def offered():
            return [o["value"] for o in pu.get_variant_dataset_options("hg38")]

        upload()
        self.index("NRG_2_hg38+walkint")
        self.assertIn("walkint", offered())
        sp.do_delete(1, "vcf:walkint")
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, "genome_library",
                                                   "NRG_2_hg38+walkint")))
        self.assertNotIn("walkint", offered())
        upload()
        self.assertIn("walkint", offered())
        self.assertEqual(self.note("walkint", 1, 1), "Up to 1 DNA/RNA bulge(s) available.")
        self.assertEqual(self.submit("walkint", 1, 1), (True, []))

    def test_a_downloaded_index_without_its_sample_list(self):
        # a downloaded index whose list is gone, with no VCF folder. The message said to
        # download the index again, which installs the list (test_samplesid_emit), but
        # Settings downloads no index that is installed, and `download --what index`
        # skips a complete install whatever its list (test_crisprme_hf): it says to
        # remove the index in Settings first, after which its download starts
        name = "NRG_1_hg38+hg38_shipped"
        os.remove(os.path.join(self.tmp, "samplesIDs", "hg38_shipped.samplesID.txt"))
        note = self.note("shipped")
        self.assertEqual(note, pu.samplesid_problem("hg38_shipped", "shipped"))
        self.assertIn("If its index came from HuggingFace, remove the index in Settings and "
                      "download it again there, which installs the list", note)
        self.assertEqual(self.submit("shipped"), (False, [note]))
        self.assertEqual(sp.add_index_hf(1, name)[2],
                         f"Index '{name}' is already installed — nothing to download.")
        sp.do_delete(1, f"index:{name}")
        with patch.object(sp, "_catalog_size", return_value=0), \
                patch.object(sp, "_preflight_disk", return_value=None), \
                patch.object(sp, "launch_settings_job", return_value="job") as launch:
            self.assertEqual(sp.add_index_hf(1, name)[0], "job")
        launch.assert_called_once_with(
            ["download", "--what", "index", "--index-name", name], "Download index")


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


class TestSubmitWithoutAGenomeOrPAM(unittest.TestCase):
    """Submit with a genome or PAM selected that is not installed, or none installed: on the
    data folder the web app starts with (index.check_directories: its folders, all empty),
    Load Example selects hg38 and 20bp-NRG-SpCas9, and check_input opened
    PAMs/20bp-NRG-SpCas9.txt (the default one, too, with no PAM selected), which is not
    there: a server error (HTTP 500), and the page did nothing. It read the genome folder of
    BED coordinates the same way. Now the missing-input dialog says what is missing, and
    where to add it."""

    NO_GENOME = ("Genome: no reference genome is installed. Add one in Settings (Add a "
                 "reference genome).")
    NO_PAM = "PAM: no PAM is installed. Add one in Settings (Define a nuclease / PAM)."

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        index.check_directories(self.tmp)
        for mod in (mp, pu):
            p = patch.object(mod, "current_working_directory", self.tmp + os.sep)
            p.start()
            self.addCleanup(p.stop)

    def submit(self, genome, pam, guides="CTAACAGTTGCTTTTATCAC", guide_type="IP",
               variant="ref", mms=4, dna=1, rna=1):
        """check_input, as the Submit button runs it: its outputs, and the dialog's lines."""
        out = mp.check_input(1, None, genome, pam, guide_type, guides, mms, dna, rna, variant,
                             False)
        return out, [li.children for li in out[-1].children[1].children]

    def install(self, genome, pam, motif):
        os.makedirs(os.path.join(self.tmp, "Genomes", genome))
        with open(os.path.join(self.tmp, "Genomes", genome, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "T" * 10 + "CTAACAGTTGCTTTTATCACAGG" + "T" * 10 + "\n")
        with open(os.path.join(self.tmp, "PAMs", f"{pam}.txt"), "w") as fh:
            fh.write("N" * 20 + motif + " 3\n")

    def test_load_example_then_submit_on_a_fresh_data_folder(self):
        guide, pam, genome, variant, mms, dna, rna = mp.load_example_data(1)[:7]
        self.assertEqual((genome, pam), ("hg38", PAM))
        # the PAM dropdown, which offers nothing, can lose the value Load Example set
        for chosen in (pam, None):
            with self.subTest(pam=chosen):
                out, msgs = self.submit(genome, chosen, guide, "IP", variant, mms, dna, rna)
                # not run: the dialog opens, the genome and PAM boxes in red
                self.assertEqual(out[:4], (None, True, "missing-input", "missing-input"))
                self.assertEqual(msgs, [self.NO_GENOME, self.NO_PAM])
        # nothing selected, and a genomic sequence given by BED coordinates
        out, msgs = self.submit(None, None, ">s\nchr1 0 43", "GS")
        self.assertEqual((out[0], msgs), (None, [self.NO_GENOME, self.NO_PAM]))

    def test_a_genome_or_pam_selected_that_is_not_installed(self):
        self.install("susScr11", "20bp-NGG-SpCas9", "NGG")
        no_hg38 = ("Genome: hg38 is not installed. Select an installed reference genome, or "
                   "add it in Settings (Add a reference genome).")
        out, msgs = self.submit("hg38", PAM)
        self.assertEqual(out[:4], (None, True, "missing-input", "missing-input"))
        self.assertEqual(msgs, [no_hg38, "PAM: 20bp-NRG-SpCas9 is not installed. Select an "
                                         "installed PAM, or add it in Settings (Define a "
                                         "nuclease / PAM)."])
        # BED coordinates in a genome not installed are not read
        out, msgs = self.submit("hg38", "20bp-NGG-SpCas9", ">s\nchr1 0 43", "GS")
        self.assertEqual((out[0], msgs), (None, [no_hg38]))
        # none selected, with some installed: as before
        self.assertEqual(self.submit(None, None)[1], ["Genome", "PAM"])
        # both installed: Submit runs the search
        self.assertEqual(self.submit("susScr11", "20bp-NGG-SpCas9")[0][0], 1)

    def test_bed_coordinates_the_genome_folder_holds_no_sequence_at(self):
        # extract_seq reads the genome from the working directory: the web app's is its
        # data folder
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.tmp)
        self.install("hg38", PAM, "NRG")
        out, msgs = self.submit("hg38", PAM, ">s\nchr1 0 43", "GS")
        self.assertEqual((out[0], msgs), (1, []))  # read, and its guide searched
        for bed in ("chrZ 0 43", "chr1 100 200"):  # no FASTA of chrZ; past chr1's end
            with self.subTest(bed=bed):
                chrom, start, stop = bed.split()
                out, msgs = self.submit("hg38", PAM, f">s\n{bed}", "GS")
                self.assertIsNone(out[0])
                self.assertIn(f"Cannot read {chrom}:{start}-{stop} from genome hg38: its "
                              f"folder holds no FASTA of {chrom}, or the coordinates are past "
                              f"the end of {chrom}", msgs)


if __name__ == "__main__":
    unittest.main()

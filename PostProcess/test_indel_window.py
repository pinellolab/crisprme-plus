#!/usr/bin/env python
"""Tests for indel_window: the fake-indel window read from the data.

CRISPRitz 2.10.0 writes 30 reference bases on either side of an indel, 2.9.1 and earlier 25
before and 27 after; the indel log's CHR field (<chrom>_<start>-<end>_<id>) holds each
window's start and end, so a reader takes each window's flanks from its own row and never
assumes a width. A target of L bases (guide and PAM) with b DNA bulges holds a record's
change only where L + b <= LEFT + 2 and L + b <= RIGHT + 1 (LEFT + 1 at an MNV or complex
record).

Run: python -m unittest test_indel_window -v
"""

import gzip
import os
import shutil
import tempfile
import unittest
from unittest import mock

import index_fit
import indel_window as iw


def row(chrom, pos, ref, alt, left, right, k=1):
    """An indel log row as CRISPRitz's add-variants writes it, of a window of left/right."""
    return "%s_%d-%d_%d\tS1\trs%d\t0.1\t%s_%d_%s_%s\t0,%d\tN" % (
        chrom, pos - 1 - left, pos - 1 + len(ref) + right, k, k, chrom, pos, ref, alt,
        left + len(alt) + right)


class TheRowsWindow(unittest.TestCase):
    def test_flanks_of_either_release(self):
        # chr22_10524814_TTAG_T, a row of CRISPRitz 2.9.1 (25/27) and one of 2.10.0 (30/30)
        self.assertEqual(iw.row_flanks("chr22_10524788-10524844_1", "chr22_10524814_TTAG_T"),
                         (25, 27))
        self.assertEqual(iw.row_flanks("chr22_10524783-10524847_1", "chr22_10524814_TTAG_T"),
                         (30, 30))

    def test_a_chromosome_name_holding_underscores_and_dashes(self):
        chr_field = "chrUn_KI2707-42v1_100-160_5"
        self.assertEqual((iw.window_start(chr_field), iw.window_end(chr_field)), (100, 160))
        self.assertEqual(iw.row_flanks(chr_field, "chrUn_KI2707-42v1_131_A_AT"), (30, 29))

    def test_a_start_past_the_pos_is_no_window(self):
        with self.assertRaises(ValueError):
            iw.left_flank(10, 20)


class ALogsWindows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def write(self, name, rows, gz=True):
        path = os.path.join(self.tmp, name)
        opener = gzip.open if gz else open
        with opener(path, "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.write("".join(r + "\n" for r in rows))
        return path

    def test_a_log_of_one_width(self):
        path = self.write("logchr1.txt.gz", [row("chr1", 100, "A", "AT", 30, 30, 1),
                                             row("chr1", 200, "ACG", "A", 30, 30, 2),
                                             row("chr1", 300, "AC", "TG", 30, 30, 3)])
        self.assertEqual(iw.log_windows(path), {(30, 30)})
        self.assertEqual(iw.log_path(self.tmp, "chr1"), path)

    def test_a_plain_log_a_header_alone_and_no_log(self):
        self.write("logchr2.txt", [row("chr2", 100, "A", "AT", 25, 27)], gz=False)
        self.write("logchr3.txt.gz", [])
        self.assertEqual(iw.dataset_windows(self.tmp),
                         {"chr2": {(25, 27)}, "chr3": set()})
        self.assertEqual(iw.dataset_windows(self.tmp, ["chr2", "chr9"]),
                         {"chr2": {(25, 27)}, "chr9": set()})
        self.assertEqual(iw.dataset_window(self.tmp), (25, 27))
        self.assertEqual(iw.log_chroms(self.tmp), ["chr2", "chr3"])

    def test_logs_of_two_enrichments(self):
        self.write("logchr1.txt.gz", [row("chr1", 100, "A", "AT", 30, 30)])
        self.write("logchr2.txt.gz", [row("chr2", 100, "A", "AT", 25, 27)])
        with self.assertRaises(ValueError):
            iw.dataset_window(self.tmp)

    def test_the_first_rows_say_the_width(self):
        rows = [row("chr1", 100 + 10 * k, "A", "AT", 30, 30, k) for k in range(5)]
        path = self.write("logchr1.txt.gz",
                          rows + ["chr1_x-y_9\tS1\t.\t0.1\tchr1_900_A_AT\t0,1\tN"])
        self.assertEqual(iw.log_windows(path, rows=5), {(30, 30)})
        with self.assertRaises(ValueError):
            iw.log_windows(path, rows=None)


class WhatAWindowHolds(unittest.TestCase):
    def test_reach(self):
        self.assertEqual((iw.reach((25, 27)), iw.reach((25, 27), mnv=True)), (27, 26))
        self.assertEqual((iw.reach((30, 30)), iw.reach((30, 30), mnv=True)), (31, 31))

    def test_a_search_that_fits_says_nothing(self):
        # every published index's motif at up to 2 DNA bulges: 29 + 2 on 30/30
        self.assertIsNone(iw.target_note((30, 30), 29, 2, "complete-search", "ds"))
        self.assertIsNone(iw.target_note((25, 27), 23, 3, "complete-search", "ds"))
        self.assertIsNone(iw.target_note(None, 40, 5, "complete-search", "ds"))

    def test_a_search_past_the_window_says_so(self):
        note = iw.target_note((30, 30), 29, 3, "complete-search", "ds")
        self.assertTrue(note.startswith("WARNING [complete-search]: the fake-indel windows of ds "
                                        "hold 30 reference bases before each indel and 30 "
                                        "after it"), note)
        self.assertIn("holds up to 31 bases, guide, PAM and DNA bulges counted; this search's "
                      "targets need up to 32 (a 29-nt guide and PAM with 3 DNA bulges)", note)
        self.assertNotIn("CRISPRitz 2.10.0 writes", note)
        # an archive of 2.9.1's windows: an MNV's limit too, and how to search them
        note = iw.target_note((25, 27), 27, 0, "indel search", "ds")
        self.assertIn("holds up to 27 bases (26 at a record whose REF and ALT differ at their "
                      "first base, an MNV or complex record)", note)
        self.assertIn("this search's targets need up to 27 (a 27-nt guide and PAM)", note)
        self.assertIn("enrich ds again with it (build-index-only)", note)


class _EngineCase(unittest.TestCase):
    """A temporary folder, and engine(): a crispritz.py of the lines given on PATH."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        index_fit._engine_cache.clear()
        self.addCleanup(index_fit._engine_cache.clear)

    def engine(self, lines):
        bindir = os.path.join(self.tmp, "bin")
        os.makedirs(bindir, exist_ok=True)
        exe = os.path.join(bindir, "crispritz.py")
        with open(exe, "w") as fh:
            fh.write("".join(x + "\n" for x in lines))
        os.chmod(exe, 0o755)
        index_fit._engine_cache.clear()
        return mock.patch.dict(os.environ, {"PATH": bindir})


class TheEngineAndTheIndex(_EngineCase):
    def test_the_engine_states_its_window(self):
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertEqual(iw.engine_window(), (30, 30))
        with self.engine(['VERSION = "2.9.1"']):
            self.assertEqual(iw.engine_window(), iw.OLD_WINDOW)
        with mock.patch.dict(os.environ, {"PATH": os.path.join(self.tmp, "none")}):
            index_fit._engine_cache.clear()
            self.assertEqual(iw.engine_window(), iw.WINDOW_2100)

    def test_an_indels_sub_index_says_its_window(self):
        sub = os.path.join(self.tmp, "NRG_2_fakechr1")
        os.makedirs(sub)
        self.assertIsNone(iw.index_window(sub))       # nothing to tell by
        with open(os.path.join(sub, "NRG_fakechr1_1.bin"), "wb") as fh:
            fh.write(b"\x01\x00\x00\x00\x14\x00\x00\x00")   # format 1: 2.9.1's, 25/27
        self.assertEqual(iw.index_window(sub), iw.OLD_WINDOW)
        with open(os.path.join(sub, index_fit.INDEX_META), "w") as fh:
            fh.write('{"index_format": 1, "indel_window": [30, 30]}')
        self.assertEqual((iw.recorded_window(sub), iw.index_window(sub)), ((30, 30), (30, 30)))


class ADatasetEnrichedAnew(_EngineCase):
    """W-R2-F1: what an enrichment left, which a new one of another CRISPRitz must not be
    mixed with: its genomes, dictionaries, indel logs, dict-less tiers and indexes, by
    absolute path (enrichment_stores); the windows its records tell (earlier_windows); and
    the refusal of build-index-only and of a search that would enrich the dataset anew beside
    the stores of an earlier enrichment of other windows (reenrichment_refusal)."""

    def put(self, *rel, text="x", data=None):
        path = os.path.join(self.tmp, "work", *rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb" if data is not None else "w") as fh:
            fh.write(data if data is not None else text)
        return path

    def dataset(self, window=(25, 27), data_type="genotyped-phased", samples=True):
        """A dataset ds of genome hg38 as build-index-only leaves it, its indel log of
        ``window``, its variant index of format 1 (2.9.1's), beside what is not its own; its
        registry's manifest of ``data_type`` (None: none, and the aggregation of a registry
        built from INFO AF fields, as every published mega's), its VCF naming a sample or,
        without ``samples``, none."""
        w = os.path.join(self.tmp, "work")
        shutil.rmtree(w, ignore_errors=True)
        self.put("Genomes", "hg38", "chr1.fa")
        self.put("Genomes", "hg38+ds", "chr1.enriched.fa")
        self.put("Genomes", "hg38+ds_INDELS", "fake_ds_chr1", "fakechr1.fa")
        self.put("Dictionaries", "dictionaries_ds", "my_dict_chr1.json.gz")
        with gzip.open(self.put("Dictionaries", "log_indels_ds", "logchr1.txt.gz", data=b""),
                       "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
                     + row("chr1", 100, "A", "AT", *window) + "\n")
        self.put("Dictionaries", "registry_ds", "reg_chr1.idx",
                 text='{"n_records": 2, "data_type": "%s"}' % data_type if data_type else
                 '{"n_records": 2, "aggregation": "info_af", "databases": {}}')
        with gzip.open(self.put("VCFs", "ds", "ds.chr1.vcf.gz", data=b""), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                     + ("\tFORMAT\tS1" if samples else "") + "\n")
        self.put("Dictionaries", "registry_ds", "reg_chr1.bin")
        self.put("Dictionaries", "genotypes_ds", "gt_chr1.idx")
        self.put("Dictionaries", "indel_genotypes_ds", "gt_indel_chr1.tsv.gz")
        fmt1 = b"\x01\x00\x00\x00\x14\x00\x00\x00"
        for name in ("NRG_3_hg38+ds", "20-NRG_2_hg38+ds", "20-NRG_2_hg38",
                     "20-NRG_2_x_hg38+ds", "20-NRG_2_hg38+ds2"):
            self.put("genome_library", name, "NRG_chr1_1.bin", data=fmt1)
        self.put("genome_library", "20-NRG_2_hg38+ds_INDELS", "NRG_2_fakechr1",
                 "NRG_fakechr1_1.bin", data=fmt1)
        os.makedirs(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds.lock"))
        os.makedirs(os.path.join(w, "genome_library", ".build-index-only-x", "NRG_3_hg38+ds"))
        return w

    def test_what_comes_from_an_enrichment(self):
        w = self.dataset()
        g, d, lib = (os.path.join(w, x) for x in ("Genomes", "Dictionaries", "genome_library"))
        indexes = [os.path.join(lib, n) for n in ("20-NRG_2_hg38+ds", "20-NRG_2_hg38+ds_INDELS",
                                                 "NRG_3_hg38+ds")]
        self.assertEqual(iw.enrichment_stores(w, "hg38", "ds"), [
            os.path.join(g, "hg38+ds"), os.path.join(g, "hg38+ds_INDELS"),
            os.path.join(d, "dictionaries_ds"), os.path.join(d, "log_indels_ds"),
            os.path.join(d, "registry_ds"), os.path.join(d, "genotypes_ds")] + indexes)
        # another library a search reads, each index once; the tiers there or not
        other = os.path.join(self.tmp, "staged")
        os.makedirs(os.path.join(other, "NRR_3_hg38+ds"))
        os.symlink(os.path.join(lib, "NRG_3_hg38+ds"), os.path.join(other, "NRG_3_hg38+ds"))
        shutil.rmtree(os.path.join(d, "genotypes_ds"))
        os.makedirs(os.path.join(d, "genotypes_ds", ".building"))
        self.assertEqual(iw.enrichment_stores(w, "hg38", "ds", [other])[4:],
                         [os.path.join(d, "registry_ds")] + indexes
                         + [os.path.join(other, "NRR_3_hg38+ds")])
        self.assertIn(os.path.join(d, "genotypes_ds"),
                      iw.enrichment_stores(w, "hg38", "ds", present=False))
        # a sites-only registry is its VCFs' (the mega's), not the enrichment's
        w = self.dataset(data_type="sites-only")
        self.assertNotIn(os.path.join(d, "registry_ds"), iw.enrichment_stores(w, "hg38", "ds"))
        # W-R3-F1: nor one whose manifest has no data_type, of INFO AF fields, as every
        # published mega's registry
        w = self.dataset(data_type=None)
        self.assertEqual(iw.registry_data_type(os.path.join(d, "registry_ds")), "sites-only")
        self.assertNotIn(os.path.join(d, "registry_ds"), iw.enrichment_stores(w, "hg38", "ds"))
        self.assertNotIn(os.path.join(d, "registry_ds"),
                         iw.enrichment_stores(w, "hg38", "ds", present=False))

    def test_the_windows_an_earlier_enrichment_tells(self):
        w = self.dataset((25, 27))
        lib = os.path.join(w, "genome_library")
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            # each told alone: its log, its INDELS index (format 1: 2.9.1's), its variant index
            os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            shutil.rmtree(os.path.join(lib, "20-NRG_2_hg38+ds_INDELS"))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            for n in ("NRG_3_hg38+ds", "20-NRG_2_hg38+ds"):
                shutil.rmtree(os.path.join(lib, n))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), set())
        # a variant index of format 1 tells nothing where the installed CRISPRitz writes it
        w = self.dataset((30, 30))
        os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
        shutil.rmtree(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds_INDELS"))
        with self.engine(['VERSION = "2.9.1"']):
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), set())

    def test_a_dataset_enriched_anew_beside_another_crispritzs_stores_is_refused(self):
        w = self.dataset((25, 27))
        g, d, lib = (os.path.join(w, x) for x in ("Genomes", "Dictionaries", "genome_library"))
        shutil.rmtree(os.path.join(g, "hg38+ds"))
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            why = iw.reenrichment_refusal(w, "hg38", "ds")
        built = [os.path.join(d, "registry_ds"), os.path.join(d, "genotypes_ds")] + [
            os.path.join(lib, n) for n in ("20-NRG_2_hg38+ds", "20-NRG_2_hg38+ds_INDELS",
                                           "NRG_3_hg38+ds")]
        self.assertTrue(why.startswith(
            "no finished enrichment of ds is in place (%s), so build-index-only would enrich "
            "it anew with the installed CRISPRitz, whose fake-indel windows are 30/30; but an "
            "earlier enrichment of it, whose records tell windows of 25/27 (CRISPRitz 2.9.1 "
            "and earlier wrote 25/27), left what was built from it: %s." % (
                os.path.join(g, "hg38+ds"), ", ".join(built))), why)
        self.assertIn("(rm -rf %s), then run build-index-only again" % " ".join(
            "'%s'" % p for p in [os.path.join(g, "hg38+ds_INDELS"),
                                 os.path.join(d, "dictionaries_ds"),
                                 os.path.join(d, "log_indels_ds")] + built), why)
        self.assertNotIn("indel_genotypes_ds", why)
        # of the installed CRISPRitz's windows: the same enrichment, nothing refused
        with self.engine(['VERSION = "2.9.1"']):
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))
        # nothing built from it: nothing to refuse
        for p in built:
            shutil.rmtree(p)
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))

    def test_the_search_scripts_check(self):
        """indel_window.py reenrichment, as the search script runs it before it enriches a
        dataset itself (STEP 1): the refusal on stdout, exit 3; nothing, exit 0."""
        import subprocess
        import sys

        w = self.dataset((25, 27))
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "indel_window.py")
        for lines, rc in (([ 'VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                            "INDEL_RIGHT = 30"], 3), (['VERSION = "2.9.1"'], 0)):
            with self.engine(lines):
                done = subprocess.run([sys.executable, script, "reenrichment", w, "hg38", "ds",
                                       os.path.join(self.tmp, "staged")],
                                      capture_output=True, text=True, env=dict(os.environ))
            self.assertEqual((done.returncode, done.stderr), (rc, ""), done.stdout)
            self.assertEqual(done.stdout.startswith(
                "no finished enrichment of ds is in place (%s), so this search would enrich "
                "it anew" % os.path.join(w, "Genomes", "hg38+ds")), rc == 3, done.stdout)
        # --vcf: the dataset's VCF folder, which tells a sites-only one (W-R3-F1), beside the
        # windows of the installed CRISPRitz too
        w = self.dataset((30, 30), samples=False)   # (its registry a genotyped dataset's)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        for n in os.listdir(os.path.join(w, "genome_library")):   # (2.9.1's indexes)
            shutil.rmtree(os.path.join(w, "genome_library", n))
        for vcf, rc in ((["--vcf", os.path.join(w, "VCFs", "ds")], 3), ([], 0)):
            with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                              "INDEL_RIGHT = 30"]):
                done = subprocess.run([sys.executable, script, "reenrichment", *vcf, w, "hg38",
                                       "ds"], capture_output=True, text=True,
                                      env=dict(os.environ))
            self.assertEqual((done.returncode, done.stderr), (rc, ""), done.stdout)
            self.assertEqual("so this search would enrich it anew from its VCFs; but ds is "
                             "sites-only" in done.stdout, rc == 3, done.stdout)


    def test_a_registry_says_its_data_type(self):
        """registry_data_type reads a manifest as crisprme.py's _registry_data_type does: its
        data_type, else sites-only for an aggregation of INFO AF fields."""
        reg = os.path.join(self.tmp, "registry_x")
        os.makedirs(reg)
        self.assertIsNone(iw.registry_data_type(reg))
        self.assertIsNone(iw.registry_data_type(os.path.join(self.tmp, "none")))
        for text, kind in (('{"aggregation": "info_af"}', "sites-only"),
                           ('{"aggregation": "carriers"}', None),
                           ('{"aggregation": "info_af", "data_type": "sites-only"}', "sites-only"),
                           ('{"aggregation": "panel", "data_type": "genotyped-phased"}',
                            "genotyped-phased"),
                           ("not json", None), ("[1, 2]", None)):
            with open(os.path.join(reg, "reg_chr1.idx"), "w") as fh:
                fh.write(text)
            self.assertEqual(iw.registry_data_type(reg), kind, text)
        # the first manifest that tells
        with open(os.path.join(reg, "reg_chr2.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        self.assertEqual(iw.registry_data_type(reg), "sites-only")

    def test_a_sites_only_dataset(self):
        """sites_only: its VCFs name no sample, where there are any to read; else its
        registry says so. A VCF of the synthetic MEGA sample (synth_sites_gt.py) names one."""
        w = self.dataset(data_type=None, samples=False)
        vcfs = os.path.join(w, "VCFs", "ds")
        self.assertTrue(iw.sites_only(w, "ds", vcfs))
        self.assertTrue(iw.sites_only(w, "ds"))            # (by its registry)
        w = self.dataset(data_type=None, samples=True)
        self.assertFalse(iw.sites_only(w, "ds", vcfs))     # (the VCFs tell first)
        self.assertTrue(iw.sites_only(w, "ds", os.path.join(w, "VCFs", "none")))
        w = self.dataset(data_type="genotyped-phased", samples=False)
        self.assertTrue(iw.sites_only(w, "ds", vcfs))
        self.assertFalse(iw.sites_only(w, "ds"))

    def test_a_sites_only_dataset_is_not_enriched_anew_beside_its_kits_windows(self):
        """W-R3-F1: the mega's windows come from its kit (a synthetic sample carries its
        indels), and an enrichment of its VCFs, which name no sample, holds none. Enriched
        anew beside the windows of its logs or INDELS indexes, whatever their width, it
        would lose every indel off-target: refused, pointing to the kit, its registry (no
        data_type, as the published mega's) not named, nor build-index-only."""
        engine30 = ['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                    "INDEL_RIGHT = 30"]
        for window in ((25, 27), (30, 30)):
            w = self.dataset(window, data_type=None, samples=False)
            g, d = os.path.join(w, "Genomes"), os.path.join(w, "Dictionaries")
            vcfs = os.path.join(w, "VCFs", "ds")
            for p in (os.path.join(g, "hg38+ds"), os.path.join(g, "hg38+ds_INDELS"),
                      os.path.join(d, "genotypes_ds"), os.path.join(d, "dictionaries_ds")):
                shutil.rmtree(p)
            for n in os.listdir(os.path.join(w, "genome_library")):   # (MG1: logs alone)
                shutil.rmtree(os.path.join(w, "genome_library", n))
            with self.engine(engine30):
                why = iw.reenrichment_refusal(w, "hg38", "ds", vcf_dir=vcfs)
            self.assertTrue(why.startswith(
                "no finished enrichment of ds is in place (%s), so build-index-only would "
                "enrich it anew from its VCFs; but ds is sites-only: its VCFs name no sample, "
                "and CRISPRitz writes a fake-indel window only for an indel a sample carries, "
                "so an enrichment of them holds none, while its indel logs or its INDELS "
                "indexes record windows of %s" % (os.path.join(g, "hg38+ds"),
                                                  iw.describe(window))), why)
            self.assertIn("Put that kit back in place (%s and %s_INDELS: of a published "
                          "dataset, its enriched-genome kit" % ((os.path.join(g, "hg38+ds"),) * 2),
                          why)
            self.assertNotIn("registry_ds", why)
            self.assertNotIn("build-index-only again", why)
            self.assertEqual("rm -rf" in why, window != (30, 30), why)
            if window != (30, 30):
                self.assertIn("install the dataset's archive and enriched-genome kit of "
                              "CRISPRme 2.10.0 in place of all that came from that kit (rm -rf "
                              "'%s')" % os.path.join(d, "log_indels_ds"), why)
            # by its registry alone, its VCFs not there
            with self.engine(engine30):
                self.assertTrue(iw.reenrichment_refusal(w, "hg38", "ds").startswith(
                    "no finished enrichment of ds is in place"))
        # its INDELS index alone tells it had windows
        w = self.dataset((30, 30), data_type=None, samples=False)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
        with self.engine(engine30):
            self.assertIn("is sites-only", iw.reenrichment_refusal(w, "hg38", "ds"))
            # nothing that holds a window, a variant index of format 1 alone (2.9.1's): no
            # window to lose, refused as any dataset beside another CRISPRitz's stores is,
            # its registry not named
            shutil.rmtree(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds_INDELS"))
            why = iw.reenrichment_refusal(w, "hg38", "ds")
            self.assertTrue(why.startswith(
                "no finished enrichment of ds is in place (%s), so build-index-only would "
                "enrich it anew with the installed CRISPRitz" % os.path.join(
                    w, "Genomes", "hg38+ds")), why)
            self.assertNotIn("registry_ds", why)
            for n in os.listdir(os.path.join(w, "genome_library")):
                shutil.rmtree(os.path.join(w, "genome_library", n))
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))
        # VCFs that name a sample (the synthetic MEGA sample) write windows: enriched anew as
        # any dataset beside the stores of the installed CRISPRitz's windows; by its registry
        # alone (its VCFs not given), refused
        w = self.dataset((30, 30), data_type=None, samples=True)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        for n in os.listdir(os.path.join(w, "genome_library")):
            shutil.rmtree(os.path.join(w, "genome_library", n))
        with self.engine(engine30):
            self.assertIsNone(iw.reenrichment_refusal(
                w, "hg38", "ds", vcf_dir=os.path.join(w, "VCFs", "ds")))
            self.assertIn("is sites-only", iw.reenrichment_refusal(w, "hg38", "ds"))

    def test_the_advice_on_a_sites_only_dataset(self):
        """W-R3-F1: a search's window note and the fake-genome refusal point a sites-only
        dataset to its archive and kit, or its kit built again, not to build-index-only."""
        note = iw.target_note((25, 27), 27, 1, "complete-search", "hg38_mega", sites=True)
        self.assertIn("CRISPRitz 2.10.0 writes windows of 30 bases on either side, which hold "
                      "up to 31. To search them: hg38_mega is sites-only: its VCFs name no "
                      "sample, and CRISPRitz writes a fake-indel window only for an indel a "
                      "sample carries, so an enrichment of them holds none (so build-index-only "
                      "would leave it with no indel off-target); its windows come from its "
                      "enriched-genome kit, whose indels a synthetic sample carries: install its "
                      "archive and enriched-genome kit of CRISPRme 2.10.0 "
                      "(docs/PRECOMPUTED_INDEXES.md), or build its kit again that way "
                      "(seq_script/merge_panels/mega_rebuild.sh).", note)
        self.assertNotIn("enrich hg38_mega again", note)
        self.assertIn("enrich ds again with it (build-index-only) to search them.",
                      iw.target_note((25, 27), 27, 1, "complete-search", "ds"))


class AFakeGenomeAndItsLogs(unittest.TestCase):
    """The fake-indel genome of a chromosome and its indel log must come from one enrichment
    (fake_window_disagreement): each row's window at the offsets its FAKEPOS field gives,
    an N after it, its left flank the reference its refseq field records (a SNP code of the
    overlay standing for its base, Q for any). A download installs an archive's logs beside an
    enrichment made by another CRISPRitz: a 25/27 fake-indel genome read with 30/30 logs, or
    the other way round, disagrees at its first window."""

    REF = (lambda r: "".join(r.choice("ACGT") for _ in range(400)))(
        __import__("random").Random(7))

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_fake_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def cut(self, left, right, name, overlay=False):
        """The fake-indel FASTA and indel log an enricher of windows left/right writes of an
        insertion at 100, a deletion at 200 and an MNV at 300 (overlay: every flank base
        coded, as the SNP overlay codes it, the first as Q)."""
        recs = [(100, self.REF[99], self.REF[99] + "TT"), (200, self.REF[199:202], self.REF[199]),
                (300, self.REF[299:301], "".join({"A": "C"}.get(b, "A") for b in self.REF[299:301]))]
        code = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
        flat, rows = "", []
        for k, (pos, ref, alt) in enumerate(recs):
            start = pos - 1 - left
            lf = self.REF[start:pos - 1]
            if overlay:
                lf = "Q" + "".join(code[b] for b in lf[1:])
            win = lf + alt + self.REF[pos - 1 + len(ref):pos - 1 + len(ref) + right]
            rows.append("chrT_%d-%d_%d\tS1\trs%d\t0.1\tchrT_%d_%s_%s\t%d,%d\t%s" % (
                start, pos - 1 + len(ref) + right, k + 1, k, pos, ref, alt, len(flat),
                len(flat) + len(win), self.REF[start:start + len(win)]))
            flat += win + "N"
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        fa = os.path.join(d, "fakechrT.fa")
        with open(fa, "w") as fh:
            fh.write(">fakechrT\n" + "".join(w + "\nN\n" for w in flat.split("N")[:-1]))
        log = os.path.join(d, "logchrT.txt.gz")
        with gzip.open(log, "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n" + "\n".join(rows) + "\n")
        return fa, log

    def test_one_enrichment_agrees(self):
        for left, right in ((25, 27), (30, 30)):
            for overlay in (False, True):
                fa, log = self.cut(left, right, "w%d_%d_%s" % (left, right, overlay), overlay)
                self.assertIsNone(iw.fake_window_disagreement(fa, log), (left, overlay))
                self.assertIsNone(iw.fake_window_disagreement(fa, log, rows=None))

    def test_two_enrichments_disagree(self):
        fa27, log27 = self.cut(25, 27, "w25", overlay=True)
        fa30, log30 = self.cut(30, 30, "w30", overlay=True)
        how = iw.fake_window_disagreement(fa27, log30)
        self.assertEqual(how, "the window of chrT_100_%s_%sTT, which the log places at 0-63 of "
                              "the fake sequence, ends elsewhere there (no separator after it)"
                         % (self.REF[99], self.REF[99]))
        self.assertIn("ends elsewhere there", iw.fake_window_disagreement(fa30, log27))
        # the right length, but another flank: a window of 30 cut 3 bases later
        with open(fa30) as fh:
            text = fh.read()
        lines = text.split("\n")
        lines[1] = self.REF[99 - 27:99 - 27 + len(lines[1])]
        with open(fa30, "w") as fh:
            fh.write("\n".join(lines))
        self.assertIn("does not hold there the 30 reference bases before the record",
                      iw.fake_window_disagreement(fa30, log30))
        # nothing to compare: no log, a log of no row, no FASTA
        self.assertIsNone(iw.fake_window_disagreement(fa27, None))
        empty = os.path.join(self.tmp, "logempty.txt")
        with open(empty, "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
        self.assertIsNone(iw.fake_window_disagreement(fa27, empty))
        self.assertIsNone(iw.fake_window_disagreement(os.path.join(self.tmp, "none.fa"), log27))

    def test_a_dataset_of_two_enrichments_is_refused(self):
        """fake_genome_refusal, as build-index-only and the search script take it."""
        w = os.path.join(self.tmp, "work")
        fa, log = self.cut(25, 27, "w25", overlay=True)
        _, log30 = self.cut(30, 30, "w30")
        fake = os.path.join(w, "Genomes", "hg38+ds_INDELS", "fake_ds_chrT")
        os.makedirs(fake)
        shutil.copy(fa, fake)
        logs = os.path.join(w, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        shutil.copy(log, logs)
        self.assertIsNone(iw.fake_genome_refusal(w, "hg38", "ds"))
        shutil.copy(log30, logs)
        why = iw.fake_genome_refusal(w, "hg38", "ds", who="this search")
        self.assertTrue(why.startswith(
            "the fake-indel genome %s_INDELS and the indel logs %s of ds come from two "
            "enrichments: on chrT, the window of chrT_100_" % (
                os.path.join(w, "Genomes", "hg38+ds"), logs)), why)
        self.assertIn("so an index that this search built from the enrichment in place", why)
        self.assertIn("(rm -rf '%s' '%s_INDELS'; of a published index" % (
            os.path.join(w, "Genomes", "hg38+ds"), os.path.join(w, "Genomes", "hg38+ds")), why)
        self.assertIn("and run build-index-only again, which enriches it anew.", why)
        # W-R3-F1: a sites-only dataset (its registry of INFO AF fields, no data_type) is not
        # sent to build-index-only, which would enrich it with no window
        reg = os.path.join(w, "Dictionaries", "registry_ds")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chrT.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        why = iw.fake_genome_refusal(w, "hg38", "ds")
        self.assertNotIn("build-index-only again", why)
        self.assertNotIn(reg, why)
        self.assertTrue(why.endswith(
            "or build that kit again with a synthetic sample "
            "(seq_script/merge_panels/mega_rebuild.sh), in place of all of the dataset "
            "that comes from an enrichment (rm -rf '%s' '%s'): ds is sites-only: its VCFs name "
            "no sample, and CRISPRitz writes a fake-indel window only for an indel a sample "
            "carries, so an enrichment of them holds none, so build-index-only would enrich it "
            "anew with no indel off-target." % (os.path.join(w, "Genomes", "hg38+ds_INDELS"),
                                                logs)), why)

    def test_the_recipe_the_advice_names_builds_the_whole_kit(self):
        """The advice on a sites-only dataset says to build its kit again with a synthetic
        sample, in place of all of it that comes from an enrichment (the enriched genome and
        its indexes among them). The recipe it names must build the whole kit of 2.10.0, the
        one the index set builds the mega with (seq_script/rebuild/index_set_2.10.0.tsv):
        lane M's mega_rebuild.sh, one build-index-only of the enriched genome, the SNP index,
        the fake-indel genome and its INDELS index. It named mega_build_indels.sh, up to
        2.8.1's recipe of the fake-indel half alone, after which the next search found no
        enriched genome and was refused, the dataset sites-only and its new logs holding
        windows (the integration of stage 2b)."""
        import csv
        import re
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        named = set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh",
                               iw.sites_only_way("hg38_mega")))
        # the fake-genome refusal of a sites-only dataset (test_a_dataset_of_two_enrichments_...)
        w = os.path.join(self.tmp, "work")
        fa, log = self.cut(25, 27, "w25", overlay=True)
        _, log30 = self.cut(30, 30, "w30")
        fake = os.path.join(w, "Genomes", "hg38+ds_INDELS", "fake_ds_chrT")
        os.makedirs(fake)
        shutil.copy(fa, fake)
        logs = os.path.join(w, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        shutil.copy(log30, logs)
        reg = os.path.join(w, "Dictionaries", "registry_ds")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chrT.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        why = iw.fake_genome_refusal(w, "hg38", "ds")
        self.assertIn("or build that kit again with a synthetic sample", why)
        named |= set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh", why))
        manifest = os.path.join(repo, "seq_script", "rebuild", "index_set_2.10.0.tsv")
        with open(manifest) as fh:
            lines = [x for x in fh if x.strip() and not x.startswith("#")]
        mega = [r for r in csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)
                if r["panel"] == "hg38_mega" and r["kind"] == "variant"]
        self.assertTrue(mega)
        builds = set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh",
                                " ".join(r["build"] for r in mega)))
        self.assertEqual(named, builds)
        for script in named:
            self.assertTrue(os.path.isfile(os.path.join(repo, script)), script)

    def test_the_overlay_refuses_a_genome_its_log_does_not_describe(self):
        """W-R3-F2: the SNP overlay places each code by the log's windows. Given the logs of
        another enrichment (30/30 beside a 25/27 fake-indel genome) it wrote 641,930 codes
        at the wrong bases of a real chr22 genome; it raises NotItsLog now, the FASTA left as
        it was, and the search's overlay step says so on stderr, exit 1. With the genome's
        own log it overlays as before."""
        import io
        from contextlib import redirect_stderr, redirect_stdout

        import overlay_indel_snps
        import search_indel_snp

        fa27, log27 = self.cut(25, 27, "w25")
        _, log30 = self.cut(30, 30, "w30")
        enriched = os.path.join(self.tmp, "chrT.enriched.fa")
        code = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
        with open(enriched, "w") as fh:
            fh.write(">chrT\n" + "".join(code[b] for b in self.REF) + "\n")
        with open(fa27) as fh:
            before = fh.read()
        with self.assertRaises(overlay_indel_snps.NotItsLog) as said:
            overlay_indel_snps.overlay_fake_chromosome(fa27, enriched, log30)
        self.assertIn("come from two enrichments: the window of chrT_100_", str(said.exception))
        with open(fa27) as fh:
            self.assertEqual(fh.read(), before)
        # the search's step (search_indel_snp.py overlay <genome> <enriched> <logs>)
        top = os.path.join(self.tmp, "search")
        for d, f, src in (("genome/fake_ds_chrT", "fakechrT.fa", fa27),
                          ("enriched", "chrT.enriched.fa", enriched),
                          ("logs", "logchrT.txt.gz", log30)):
            os.makedirs(os.path.join(top, d))
            shutil.copy(src, os.path.join(top, d, f))
        err = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(err):
            rc = search_indel_snp.main(["search_indel_snp.py", "overlay",
                                        os.path.join(top, "genome"), os.path.join(top, "enriched"),
                                        os.path.join(top, "logs")])
        self.assertEqual(rc, 1)
        self.assertTrue(err.getvalue().startswith("ERROR: the fake-indel genome "), err.getvalue())
        self.assertIn("come from two enrichments", err.getvalue())
        with open(os.path.join(top, "genome", "fake_ds_chrT", "fakechrT.fa")) as fh:
            self.assertEqual(fh.read(), before)
        # its own log: overlaid
        self.assertGreater(overlay_indel_snps.overlay_fake_chromosome(fa27, enriched, log27), 0)
        self.assertIsNone(iw.fake_window_disagreement(fa27, log27))


if __name__ == "__main__":
    unittest.main(verbosity=2)

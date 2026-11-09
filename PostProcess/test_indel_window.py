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


class TheEngineAndTheIndex(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main(verbosity=2)

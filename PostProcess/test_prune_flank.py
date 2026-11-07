#!/usr/bin/env python3
"""An all-N index built before CRISPRitz 2.9.1 is told by its trees, and a search under
--prune-pam that reads an index whose flank is wider than its PAM says what that costs.

CRISPRitz before 2.9.1 (CRISPRme 2.8.0) dropped every target whose flank held an N, so such an
index searched under --prune-pam for a narrower PAM misses the sites whose flank holds an N
(idx audit IDX-5: 3 of 36 rows on the synthetic genome, 2.6% of the indel rows on chr22). The
index records nothing of its engine; its trees do, as an N nibble (0) in a flank, which no
earlier build stores (prune_flank). A flank wider than the PAM also cannot index a PAM at a
contig's end (IDX-4): complete-search says both on stdout, and into the run's Params.txt; the
search script, of a search it prunes by itself (one the web launched), in its log
(test_defaults_and_index_flex).

The trees are fabricated as buildTST writes them (measured on real ones: for an odd flank the
last byte's high nibble is padding, its low nibble the last base).
"""

import contextlib
import io
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import prune_flank  # noqa: E402

CODE = {"A": 1, "C": 2, "G": 4, "T": 8, "N": 0}


def tree(path, k, flanks, guide=20, twins=()):
    """A .bin of one target a flank (a string of k bases), packed as buildTST packs them;
    the targets of ``twins`` carry the '_' flag and its int."""
    width = (k + 1) // 2
    out = struct.pack("=ii", len(flanks), guide)
    for i, flank in enumerate(flanks):
        codes = [CODE[b] for b in flank]
        if k % 2:
            codes = codes[:-1] + [0, codes[-1]]  # (the padding before the last base)
        out += struct.pack("=i", i) + bytes(codes[j] << 4 | codes[j + 1]
                                            for j in range(0, 2 * width, 2))
        out += b"_" + struct.pack("=i", -1) if i in twins else b"0"
    out += struct.pack("=i", 1) + b"\x10"  # (a tree after the targets, not read)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(out)
    return path


class _Dir(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="prune_flank_")
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.root, *parts)


class AFlank(_Dir):
    def test_an_even_flank_with_an_n_or_none(self):
        self.assertTrue(prune_flank.flank_holds_n(
            tree(self.path("a.bin"), 8, ["ACGTACGT", "ACGTNCGT"]), 8))
        self.assertFalse(prune_flank.flank_holds_n(
            tree(self.path("b.bin"), 8, ["ACGTACGT", "TTTTAAAA"]), 8))

    def test_an_odd_flank_is_not_read_into_its_padding(self):
        self.assertFalse(prune_flank.flank_holds_n(
            tree(self.path("a.bin"), 3, ["ACG", "TTT"]), 3))
        self.assertTrue(prune_flank.flank_holds_n(tree(self.path("b.bin"), 3, ["ACN"]), 3))
        self.assertTrue(prune_flank.flank_holds_n(tree(self.path("c.bin"), 3, ["NCG"]), 3))

    def test_a_target_of_a_sequence_seen_before_is_stepped_over(self):
        self.assertTrue(prune_flank.flank_holds_n(
            tree(self.path("a.bin"), 8, ["ACGTACGT", "ACGTACGT", "NNNNACGT"], twins=(1,)), 8))
        self.assertFalse(prune_flank.flank_holds_n(
            tree(self.path("b.bin"), 8, ["ACGTACGT", "ACGTACGT"], twins=(1,)), 8))

    def test_a_file_that_is_no_tree_tells_nothing(self):
        with open(self.path("x.bin"), "wb") as fh:
            fh.write(b"\x00\x00")
        self.assertIsNone(prune_flank.flank_holds_n(self.path("x.bin"), 8))
        self.assertIsNone(prune_flank.flank_holds_n(self.path("none.bin"), 8))


class AReferenceIndex(_Dir):
    def setUp(self):
        super().setUp()
        self.genome = self.path("Genomes", "g")
        os.makedirs(self.genome)
        for contig, seq in (("chrA", "ACGT" * 50), ("chrN", "ACGT" * 20 + "N" * 30 + "ACGT" * 20),
                            ("chrB", "N" * 500 + "ACGT" * 50)):
            with open(os.path.join(self.genome, f"{contig}.fa"), "w") as fh:
                fh.write(f">{contig}\n{seq}\n")
        self.index = self.path("genome_library", "NNNNNNNN_3_g")

    def test_a_fixed_build_holds_an_n_in_a_flank_of_a_contig_with_one(self):
        tree(os.path.join(self.index, "NNNNNNNN_chrA_1.bin"), 8, ["ACGTACGT"])
        tree(os.path.join(self.index, "NNNNNNNN_chrN_1.bin"), 8, ["ACGTACGT", "ACGNNNNN"])
        self.assertFalse(prune_flank.built_before_n_fix(self.index, self.genome))

    def test_one_built_before_holds_none(self):
        tree(os.path.join(self.index, "NNNNNNNN_chrA_1.bin"), 8, ["ACGTACGT"])
        tree(os.path.join(self.index, "NNNNNNNN_chrN_1.bin"), 8, ["ACGTACGT", "TTTTACGT"])
        self.assertTrue(prune_flank.built_before_n_fix(self.index, self.genome))

    def test_a_contig_of_several_trees_or_with_no_n_tells_nothing(self):
        tree(os.path.join(self.index, "NNNNNNNN_chrA_1.bin"), 8, ["ACGTACGT"])
        tree(os.path.join(self.index, "NNNNNNNN_chrN_1.bin"), 8, ["ACGTACGT"])
        tree(os.path.join(self.index, "NNNNNNNN_chrN_2.bin"), 8, ["ACGTACGT"])
        self.assertIsNone(prune_flank.built_before_n_fix(self.index, self.genome))


class AnIndelsIndex(_Dir):
    def test_a_tree_of_several_windows_tells_and_a_small_one_does_not(self):
        old = prune_flank.ONE_WINDOW
        self.addCleanup(setattr, prune_flank, "ONE_WINDOW", old)
        prune_flank.ONE_WINDOW = 2
        indels = self.path("genome_library", "NNNNNNNN_3_g+ds_INDELS")
        tree(os.path.join(indels, "NNNNNNNN_3_fakechr1", "NNNNNNNN_fakechr1_1.bin"), 8,
             ["ACGTACGT", "ACGTACGT"])  # one window: no N to tell by
        tree(os.path.join(indels, "NNNNNNNN_3_fakechr2", "NNNNNNNN_fakechr2_1.bin"), 8,
             ["ACGTACGT", "ACGTACGT", "TTTTTTTT"])
        self.assertTrue(prune_flank.indels_built_before_n_fix(indels))
        tree(os.path.join(indels, "NNNNNNNN_3_fakechr2", "NNNNNNNN_fakechr2_1.bin"), 8,
             ["ACGTACGT", "ACGTNNNN", "TTTTTTTT"])
        self.assertFalse(prune_flank.indels_built_before_n_fix(indels))
        prune_flank.ONE_WINDOW = 10
        self.assertIsNone(prune_flank.indels_built_before_n_fix(indels))


def _load_crisprme():
    """crisprme.py's definitions, as test_bulge_cap loads them."""
    import test_bulge_cap
    return test_bulge_cap.crisprme


class TheSearch(_Dir):
    """complete-search under --prune-pam says it once a search reads an all-N index wider
    than its PAM, and warns of one built before the fix (crisprme._wider_flank_notes), of
    the index _resolve_index_dir takes (crisprme._index_the_search_takes)."""

    def setUp(self):
        super().setUp()
        self.crisprme = _load_crisprme()
        self.lib = self.path("genome_library")
        self.genome = self.path("Genomes", "g")
        os.makedirs(self.genome)
        with open(os.path.join(self.genome, "chrN.fa"), "w") as fh:
            fh.write(">chrN\n" + "ACGT" * 20 + "N" * 30 + "ACGT" * 20 + "\n")

    def index(self, name, flanks=("ACGTACGT", "ACGNNNNN")):
        k = len(name.split("_")[0])
        tree(os.path.join(self.lib, name, f"{name.split('_')[0]}_chrN_1.bin"), k,
             [f[:k].ljust(k, "A") for f in flanks])
        return os.path.join(self.lib, name)

    def takes(self, prune=True, three_prime=True, motif="NRG"):
        return self.crisprme._index_the_search_takes(
            [self.lib], motif, "g", 2, prune, three_prime, 20, 1, 1, False)

    def test_the_index_taken_is_the_search_scripts(self):
        k8 = self.index("NNNNNNNN_3_g")
        self.assertEqual(self.takes(), k8)
        self.assertIsNone(self.takes(prune=False))
        self.assertIsNone(self.takes(three_prime=False, motif="TTTV"))
        k3 = self.index("NNN_3_g")
        self.assertEqual(self.takes(), k3)  # (the narrowest flank)
        nrg = self.index("NRG_3_g")
        self.assertEqual(self.takes(), nrg)  # (the search's own PAM first)

    def test_a_wider_flank_is_said_once_and_a_build_before_the_fix_warned_of(self):
        k8 = self.index("NNNNNNNN_3_g")
        notes = self.crisprme._wider_flank_notes([k8, k8], "NRG", self.genome)
        self.assertEqual(len(notes), 1)
        self.assertTrue(notes[0].startswith(
            "NOTE [complete-search]: --prune-pam serves the 3-nt PAM NRG from NNNNNNNN_3_g, "
            "whose stored flank is 8 nt: a + strand PAM with 1 to 5 nt after it to the end of a "
            "contig, and a - strand PAM in the first 5 nt of a contig, are not indexed"), notes)
        old = self.index("NNNNNNNN_4_g", flanks=("ACGTACGT", "TTTTACGT"))
        notes = self.crisprme._wider_flank_notes([old], "NRG", self.genome)
        self.assertEqual(len(notes), 2)
        self.assertTrue(notes[1].startswith(
            "WARNING [complete-search]: NNNNNNNN_4_g stores no N in the flank of any target "
            "where its genome has one: it was built before CRISPRitz 2.9.1"), notes)

    def test_the_warning_names_a_route_that_serves_the_search_whole(self):
        """Its advice was to "search without --prune-pam and --index-path, which builds a
        dedicated index": since 2.9.0 such a search prunes by itself from that very index
        (_prune_on_its_own) and builds none, so the advice sent the user round the same
        warning. It names the two routes that serve the search whole: the index built
        again, or a dedicated index of the PAM beside it, which the search takes first."""
        from unittest import mock
        old = self.index("NNNNNNNN_4_g", flanks=("ACGTACGT", "TTTTACGT"))
        warning = self.crisprme._wider_flank_notes([old], "NRG", self.genome)[1]
        self.assertNotIn("without --prune-pam", warning)
        self.assertTrue(warning.endswith(
            "Remove it and build it again with build-index-only (which keeps an index it "
            "finds in place), or put a dedicated NRG index beside it, which a search takes "
            "before it (build one with build-index-only, or download one if it is "
            "published)."), warning)
        # the search the older advice named takes that same index, pruning by itself (the
        # library is the working directory's), and builds none; with a dedicated index
        # beside it, that one answers it, unpruned
        with mock.patch.object(self.crisprme, "current_working_directory", self.root), \
                mock.patch.object(self.crisprme, "_engine_prunes", return_value=True):
            self.assertEqual(
                self.crisprme._prune_on_its_own("_", "g", [], "NRG", 2, True, 20, 1, 1), old)
            dedicated = self.index("NRG_3_g", flanks=("TTT",))
            self.assertIsNone(
                self.crisprme._prune_on_its_own("_", "g", [], "NRG", 2, True, 20, 1, 1))
        self.assertEqual(self.takes(prune=False), dedicated)
        self.assertEqual(self.crisprme._wider_flank_notes([dedicated], "NRG", self.genome), [])

    def test_an_index_of_the_pams_width_says_nothing(self):
        for name in ("NNN_3_g", "NRG_3_g"):
            self.assertEqual(self.crisprme._wider_flank_notes(
                [self.index(name, flanks=("TTT",))], "NRG", self.genome), [])


class TheSearchScriptsCommand(_Dir):
    """The search script runs prune_flank.py for a search it prunes by itself (one the web
    launched, which never passes --prune-pam): the notes complete-search says on its console,
    on stdout, with the web app's route in the warning's advice, and nothing on stderr, which
    fails the stage, whatever it is handed (an index or a genome that is not there, none)."""

    def run_it(self, *args):
        return subprocess.run([sys.executable, os.path.join(HERE, "prune_flank.py"), *args],
                              capture_output=True, text=True)

    def test_its_notes_on_stdout_and_nothing_on_stderr(self):
        genome = self.path("Genomes", "g")
        os.makedirs(genome)
        with open(os.path.join(genome, "chrN.fa"), "w") as fh:
            fh.write(">chrN\n" + "ACGT" * 20 + "N" * 30 + "ACGT" * 20 + "\n")
        old = self.path("genome_library", "NNNNNNNN_3_g")
        tree(os.path.join(old, "NNNNNNNN_chrN_1.bin"), 8, ["ACGTACGT", "TTTTACGT"])
        run = self.run_it("NGG", genome, old, "", old)  # (once each; "": none taken)
        self.assertEqual((run.returncode, run.stderr), (0, ""))
        lines = run.stdout.splitlines()
        self.assertEqual(len(lines), 2, run.stdout)
        self.assertTrue(lines[0].startswith(
            "NOTE: --prune-pam serves the 3-nt PAM NGG from NNNNNNNN_3_g, whose stored flank is "
            "8 nt"), lines[0])
        self.assertTrue(lines[1].startswith(
            "WARNING: NNNNNNNN_3_g stores no N in the flank of any target where its genome has "
            "one"), lines[1])
        self.assertIn("(build one with crisprme.py build-index-only or in the web app's "
                      "Settings, or download one if it is published)", lines[1])
        # an index and a genome that are not there: the flank, from the index's name, alone
        run = self.run_it("NGG", self.path("none"), self.path("genome_library", "NNNNNNNN_3_x"),
                          self.path("genome_library", "NNNNNNNN_3_x+ds"))
        self.assertEqual((run.returncode, run.stderr), (0, ""))
        self.assertEqual(len(run.stdout.splitlines()), 1, run.stdout)
        # none handed, or none wider than the PAM: nothing
        for args in (("NGG", genome), ("NGG", genome, "", self.path("genome_library", "NNN_3_g"))):
            with self.subTest(args=args[2:]):
                run = self.run_it(*args)
                self.assertEqual((run.returncode, run.stdout, run.stderr), (0, "", ""))


if __name__ == "__main__":
    unittest.main(verbosity=2)

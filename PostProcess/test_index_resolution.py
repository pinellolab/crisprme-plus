#!/usr/bin/env python3
"""Which installed index answers a search must not depend on the environment.

An index directory is named ``<PAM>_<N>_<suffix>``. Several can be sufficient for
one search -- a dedicated PAM, a PAM-less stand-in, the same PAM at a larger bulge
budget -- so the pipeline picks one. It used to pick whichever the shell glob
produced first, and glob order follows ``LC_COLLATE``: ``NNNNNNNN_3_hg38`` sorts
before ``NNN_3_hg38`` under ``LC_ALL=C`` (as inside the container) and after it
under ``en_US.UTF-8``. Publishing a second PAM-less index of a different width
would therefore have made the choice depend on the locale.

These tests run the *shipped* ``_resolve_index_dir`` out of
``submit_job_automated_new_multiple_vcfs.sh`` -- extracted, not copied -- against
real directory trees, so they fail if the shell function and this file drift apart.
``EveryScanInTheShippedScript`` goes one step further and runs the whole script,
with a stub engine on PATH, up to its variant search, so the three places that call
the resolver are tested as called, not only the function they share.

The resolver skips an index whose build lock is still there, so the locks are tested
here too: ``ADeadRunsLockDoesNotWedgeTheNext`` and ``ARunThatEndsMidBuildLeavesItsLock``
run the shipped lock helpers, and ``BuildLocksInTheShippedScript`` kills, interrupts
and races whole runs. They read the process table from /proc, as the script does, for
the image has no ps (``TheLockTestsNeedNoPs``). ``TheThreadGuardSizesTheIndexSearched``
covers the RAM guard, which must size a search by the index the resolver picked.
``AFinishedRunRestoresOnlyRealAnnotations`` runs the script's last annotation step,
which must leave the empty placeholder alone, and
``TheSamplesPlaceholderIsNeverWritten`` runs the script with the empty samples
placeholder, which it must leave alone too, and its samples-list blocks with a list of
the user's own that has the placeholder's name, which they must still read whole.
"""

import gzip
import os
import re
import shlex
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

_SCRIPT = Path(__file__).resolve().parent / "submit_job_automated_new_multiple_vcfs.sh"
_PAMS = _SCRIPT.parent.parent / "PAMs"


def _extract_resolver():
    """Pull `_resolve_index_dir` verbatim out of the shipped script, with the
    `_index_side` it calls, which the script defines just before it."""
    text = _SCRIPT.read_text()
    start = text.index("\t_index_side() {")
    resolver = text.index("\t_resolve_index_dir() {", start)
    # the body is tab-indented one level; its closing brace is the first "\t}"
    end = text.index("\n\t}\n", resolver) + len("\n\t}\n")
    body = text[start:end]
    # strip the one leading tab so it is a top-level function in the harness
    return "\n".join(line[1:] if line.startswith("\t") else line
                     for line in body.splitlines())


RESOLVER = _extract_resolver()


def resolve(libdir, pam_want, suffix, bmax, locale="C", prune=0, pos=1, any_side="",
            guide=0, bdna=0, brna=0):
    """Run the shipped resolver; return the chosen basename, or "" for none.

    prune/pos mirror the script's own `prune_pam` and `pos` (the PAM-file offset:
    > 0 is a 3' PAM). The resolver reads them to decide whether a flank WIDER than
    the PAM may be used. any_side is its optional fifth argument. guide, bdna and brna
    are the search's guide length (its PAM file then has that spacer) and bulges, the
    script's own `fullseqpam`, `bDNA` and `bRNA`, which it fits an index to.
    """
    fullseq = "N" * (guide + abs(pos)) if guide else ""
    script = (f"prune_pam={prune}\npos={pos}\nfullseqpam={fullseq}\nbDNA={bdna}\nbRNA={brna}\n"
              f"{RESOLVER}\n"
              f"_resolve_index_dir \"$1\" \"$2\" \"$3\" \"$4\" \"$5\"\n")
    env = dict(os.environ, LC_ALL=locale, LANG=locale)
    out = subprocess.run(
        ["bash", "-c", script, "bash", str(libdir), pam_want, suffix, str(bmax), any_side],
        capture_output=True, text=True, env=env,
    )
    if out.returncode not in (0, 1):
        raise AssertionError(f"resolver errored rc={out.returncode}: {out.stderr}")
    chosen = out.stdout.strip()
    return Path(chosen).name if chosen else ""


def resolver_stderr(libdir, pam_want, suffix, bmax, prune=0, pos=1):
    """What the shipped resolver writes to stderr, which fails the enclosing stage."""
    script = (f"prune_pam={prune}\npos={pos}\n{RESOLVER}\n"
              f"_resolve_index_dir \"$1\" \"$2\" \"$3\" \"$4\"\n")
    return subprocess.run(
        ["bash", "-c", script, "bash", str(libdir), pam_want, suffix, str(bmax)],
        capture_output=True, text=True, env=dict(os.environ, LC_ALL="C", LANG="C"),
    ).stderr


def built_for_guide(index_dir, guide):
    """The trees of <index_dir> record the guide length `guide`, as buildTST writes it:
    the second int32 of a .bin."""
    Path(index_dir, "chr1.bin").write_bytes(struct.pack("=ii", 7, guide) + b"\0" * 16)


def built_for(index_dir, pos):
    """Record in <index_dir>/.pam_build, as build-index-only does, that the index was
    built from an all-N PAM file of width |pos| after a 20-nt guide: a 3' PAM for
    pos > 0, a 5' one for pos < 0."""
    width = abs(pos)
    Path(index_dir, ".pam_build").write_text(
        f"{20 + width}bp-{'N' * width}-NO-PAM {20 + width} {pos}\n")


class ResolverHarness(unittest.TestCase):
    def test_the_function_was_actually_found(self):
        self.assertIn("_resolve_index_dir() {", RESOLVER)
        self.assertIn("_pam_want", RESOLVER)
        # the whole point: no reliance on glob order
        self.assertIn("_best_exact", RESOLVER)


class _Tree(unittest.TestCase):
    """Builds a genome_library with the given index directory names."""

    def lib(self, *names):
        """Each name becomes a FINISHED index: a directory that holds a .bin.

        The .bin matters -- the resolver treats a directory with none as absent,
        because that is exactly what an index looks like before its build writes
        the first chromosome.
        """
        d = Path(self.tmp.name)
        for n in names:
            (d / n).mkdir(parents=True, exist_ok=True)
            (d / n / "chr1.bin").write_bytes(b"")
        return d

    def lib_unfinished(self, *names):
        """Each name becomes a directory with NO .bin -- a build not yet started."""
        d = Path(self.tmp.name)
        for n in names:
            (d / n).mkdir(parents=True, exist_ok=True)
        return d

    def lock(self, libdir, name):
        """Take the build lock on one index, as `index-genome` holds it while writing."""
        (Path(libdir) / f"{name}.lock").mkdir()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)


class PicksASufficientIndex(_Tree):
    def test_exact_pam(self):
        lib = self.lib("NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_nothing_when_budget_too_small(self):
        lib = self.lib("NRG_2_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_smallest_sufficient_budget_wins(self):
        lib = self.lib("NRG_3_hg38", "NRG_5_hg38", "NRG_4_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_smallest_sufficient_budget_wins_against_glob_order(self):
        # NRG_10_hg38 globs before NRG_3_hg38 ('1' < '3'), so a resolver that kept
        # its first hit would take the larger budget; with single-digit budgets the
        # smallest one always globs first and the rule goes untested
        lib = self.lib("NRG_10_hg38", "NRG_3_hg38")
        for loc in ChoiceIsLocaleIndependent.LOCALES:
            self.assertEqual(resolve(lib, "NRG", "hg38", 3, locale=loc), "NRG_3_hg38", loc)

    def test_larger_budget_used_when_it_is_the_only_one(self):
        lib = self.lib("NRG_5_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_5_hg38")

    def test_indels_companion_is_never_selected(self):
        lib = self.lib("NRG_3_hg38_INDELS")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_non_numeric_budget_is_skipped(self):
        lib = self.lib("NRG_x_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_variant_suffix(self):
        lib = self.lib("NRG_3_hg38", "NRG_3_hg38+hg38_HGDP1kGP")
        self.assertEqual(resolve(lib, "NRG", "hg38+hg38_HGDP1kGP", 3),
                         "NRG_3_hg38+hg38_HGDP1kGP")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")


class PamLessStandIn(_Tree):
    def test_pamless_serves_a_dedicated_search(self):
        lib = self.lib("NNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NNN_3_hg38")

    def test_a_different_dedicated_pam_does_not(self):
        lib = self.lib("NGG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_nnnrrt_is_a_dedicated_pam_not_a_pamless_index(self):
        # PAMs/21bp-NNNRRT-SaCas9KKH.txt ships, so NNNRRT_3_* can exist. A glob of
        # "NNN*" matches it, but it constrains positions 4-6 and cannot stand in
        # for another PAM.
        lib = self.lib("NNNRRT_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")
        # ... while it still serves its own search
        self.assertEqual(resolve(lib, "NNNRRT", "hg38", 3), "NNNRRT_3_hg38")

    def test_exact_pam_beats_a_pamless_stand_in(self):
        lib = self.lib("NNN_3_hg38", "NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_exact_pam_wins_even_at_a_larger_budget(self):
        # no padding and no reliance on the post-search PAM filter beats a tighter
        # bulge budget
        lib = self.lib("NNN_3_hg38", "NRG_4_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_4_hg38")

    def test_narrowest_pamless_flank_wins(self):
        # least N-padding for --prune-pam
        lib = self.lib("NNNNNNNN_3_hg38", "NNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=1), "NNN_3_hg38")
        # ... and with pruning off only the equal-width one is eligible at all
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=0), "NNN_3_hg38")

    def test_wide_pamless_used_only_under_prune(self):
        """The engine reads as many stored PAM nibbles per leaf as the search PAM is
        long, and only --prune-pam reconciles that width. Picking a wider index for a
        non-pruned search hands the engine an index it must refuse -- and before
        --prune-pam the scan could not see these at all, so a run that used to build a
        dedicated index would start failing."""
        lib = self.lib("NNNNNNNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=0), "")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")

    def test_wide_pamless_refused_for_a_5_prime_pam(self):
        """The engine refuses a shorter 5' PAM off a wider flank: the emitted sequence
        is trimmed at the end only, so a front pad lands inside the PAM region."""
        lib = self.lib("NNNNNNNN_3_hg38")
        built_for(lib / "NNNNNNNN_3_hg38", -8)  # a 5' build: only its width refuses it
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=1, pos=-4), "")

    def test_equal_width_pamless_serves_its_own_side_without_prune(self):
        lib = self.lib("NNN_3_hg38", "NNNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=0), "NNN_3_hg38")
        built_for(lib / "NNNN_3_hg38", -4)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=0, pos=-4), "NNNN_3_hg38")


class APamLessIndexServesOnlyTheSideItWasBuiltFor(_Tree):
    """An all-N index stores its flank on the side of the guide where the PAM file it
    was built from puts the PAM, and the engine reads the flank on the side the search
    PAM file names. Searched with a PAM on the other side it returns no target, and no
    error: build-index-only with a 3' all-N PAM file of width 4 gives NNNN_2, on which
    a TTTV (5') search reported 0 rows and exited 0, where a dedicated TTTV_2 index
    finds 3 sites. The resolver took any all-N index of the PAM's width for either
    side. It now reads the side from the index's .pam_build, and takes an index
    without one for a 3' build, as every published one is."""

    def test_a_3_prime_build_does_not_serve_a_5_prime_pam(self):
        lib = self.lib("NNNN_3_hg38")
        built_for(lib / "NNNN_3_hg38", 4)
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=prune, pos=-4), "")

    def test_a_5_prime_build_serves_a_5_prime_pam_of_its_width(self):
        lib = self.lib("NNNN_3_hg38")
        built_for(lib / "NNNN_3_hg38", -4)
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=prune, pos=-4),
                             "NNNN_3_hg38")

    def test_an_index_without_a_record_is_a_3_prime_build(self):
        # at the widths of the published ones, 3 and 8 (the K=8 archive has no .pam_build)
        lib = self.lib("NNN_3_hg38", "NNNNNNNN_3_hg38")
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "TTT", "hg38", 3, prune=prune, pos=-3), "")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=1), "NNN_3_hg38")
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")

    def test_one_of_another_width_without_a_record_is_taken_for_neither_side(self):
        # A search before 2.8.1 built NNNN_2 from a 5' PAM file (23bp-NNNN-NO-PAM) with no
        # record; a 3' NGCG search then took it, and reported a site that is in no
        # chromosome, where 2.8.0 built NGCG_2 and found 2 (flows OF-8, new in 2.8.1).
        lib = self.lib("NNNN_3_hg38", "NNNNNN_3_hg38", "NN_3_hg38")
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "NGCG", "hg38", 3, prune=prune), "")
            self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=prune), "")
            self.assertEqual(resolve(lib, "NG", "hg38", 3, prune=prune), "")
            self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=prune, pos=-4), "")
        # ... but under the search PAM's own name, as 2.8.0 took it
        self.assertEqual(resolve(lib, "NNNN", "hg38", 3), "NNNN_3_hg38")
        # and with its record, on the side it records
        built_for(lib / "NNNN_3_hg38", 4)
        self.assertEqual(resolve(lib, "NGCG", "hg38", 3), "NNNN_3_hg38")
        built_for(lib / "NNNN_3_hg38", -4)
        self.assertEqual(resolve(lib, "NGCG", "hg38", 3), "")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4), "NNNN_3_hg38")

    def test_a_5_prime_build_does_not_serve_a_3_prime_pam(self):
        # at the PAM's width, and wider under --prune-pam: neither stores a 3' flank
        lib = self.lib("NNN_3_hg38", "NNNNNNNN_3_hg38")
        built_for(lib / "NNN_3_hg38", -3)
        built_for(lib / "NNNNNNNN_3_hg38", -8)
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=prune), "")

    def test_a_3_prime_build_still_serves_a_3_prime_pam(self):
        lib = self.lib("NNN_3_hg38", "NNNNNNNN_3_hg38")
        built_for(lib / "NNN_3_hg38", 3)
        built_for(lib / "NNNNNNNN_3_hg38", 8)
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NNN_3_hg38")
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")

    def test_an_all_n_search_pam_needs_its_side_too(self):
        # an all-N search PAM names its index exactly, and is still read on its side
        lib = self.lib("NNN_3_hg38")
        built_for(lib / "NNN_3_hg38", 3)
        self.assertEqual(resolve(lib, "NNN", "hg38", 3, pos=3), "NNN_3_hg38")
        self.assertEqual(resolve(lib, "NNN", "hg38", 3, pos=-3), "")

    def test_any_side_takes_one_built_for_the_other_side(self):
        # only to name it in an error: the search never uses what it returns
        lib = self.lib("NNNN_3_hg38", "TTTV_2_hg38")
        built_for(lib / "NNNN_3_hg38", 4)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4), "")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4, any_side="any"),
                         "NNNN_3_hg38")
        # and nothing else is relaxed: the budget, the width and a dedicated PAM hold
        self.assertEqual(resolve(lib, "TTTV", "hg38", 4, pos=-4, any_side="any"), "")
        self.assertEqual(resolve(lib, "NGG", "hg38", 3, any_side="any"), "")
        self.assertEqual(resolve(lib, "NNNGRRT", "hg38", 3, prune=1, any_side="any"), "")

    def test_a_dedicated_index_needs_no_record(self):
        # a search writes no .pam_build when it builds a dedicated index: the name of
        # a dedicated index says what it serves
        lib = self.lib("TTTV_3_hg38", "NRG_3_hg38")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4), "TTTV_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_a_bad_record_is_a_3_prime_build_and_writes_nothing_to_stderr(self):
        # the resolver's stderr fails the stage that runs it
        lib = self.lib("NNNN_3_hg38")
        rec = lib / "NNNN_3_hg38" / ".pam_build"
        for text, served in (("", False), ("\n", False), ("garbage\n", False),
                             ("27bp-NNNN-5p 27 four\n", False),
                             ("27bp-NNNN-5p 27 -4", True),
                             ("27bp-NNNN-5p 27 -4\r\n", True)):
            with self.subTest(record=text):
                rec.write_bytes(text.encode())
                self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4),
                                 "NNNN_3_hg38" if served else "")
                self.assertEqual(resolver_stderr(lib, "TTTV", "hg38", 3, pos=-4), "")
        rec.unlink()
        rec.mkdir()  # not a file
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4), "")
        self.assertEqual(resolver_stderr(lib, "TTTV", "hg38", 3, pos=-4), "")
        rec.rmdir()
        if os.geteuid() != 0:  # root reads it anyway
            built_for(lib / "NNNN_3_hg38", -4)
            rec.chmod(0)
            self.addCleanup(rec.chmod, 0o644)
            self.assertEqual(resolve(lib, "TTTV", "hg38", 3, pos=-4), "")
            self.assertEqual(resolver_stderr(lib, "TTTV", "hg38", 3, pos=-4), "")


class ADedicatedIndexServesOnlyTheSideItRecords(_Tree):
    """A dedicated index is named after its PAM alone, and a PAM file of either side names
    it alike: 'NNNNNNNNNNNNNNNNNNNNTTTV 4' builds the TTTV_<N>_<genome> the 5' Cas12a file
    does. The side check covered all-N indexes only, so a 5' TTTV search took a 3' build
    and reported one site that is in no chromosome instead of 20 (cli5 D2). It now reads
    the side any index records; one with no record has no side to contradict."""

    def record(self, lib, name, pos):
        Path(lib, name, ".pam_build").write_text(f"x {20 + abs(pos)} {pos}\n")

    def test_a_dedicated_index_recorded_for_the_other_side_does_not_serve(self):
        lib = self.lib("TTTV_2_hg38")
        self.record(lib, "TTTV_2_hg38", 4)
        for prune in (0, 1):
            self.assertEqual(resolve(lib, "TTTV", "hg38", 2, prune=prune, pos=-4), "")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=4), "TTTV_2_hg38")

    def test_the_one_of_its_side_is_taken_beside_it(self):
        # (the smaller budget of the wrong side used to win)
        lib = self.lib("TTTV_2_hg38", "TTTV_3_hg38")
        self.record(lib, "TTTV_2_hg38", 4)
        self.record(lib, "TTTV_3_hg38", -4)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4), "TTTV_3_hg38")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=4), "TTTV_2_hg38")

    def test_one_with_no_record_serves_as_before(self):
        lib = self.lib("TTTV_2_hg38", "NRG_3_hg38")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4), "TTTV_2_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_any_side_still_finds_it_to_say_why(self):
        lib = self.lib("TTTV_2_hg38")
        self.record(lib, "TTTV_2_hg38", 4)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4, any_side="any"),
                         "TTTV_2_hg38")


class FlankMustBeWideEnough(_Tree):
    """A stored flank of width K serves a PAM of length L only when L <= K.

    The index never stored the positions a longer PAM would constrain, so picking a
    too-narrow PAM-less index just defers the failure to the engine.
    """

    def test_narrow_pamless_does_not_serve_a_longer_pam(self):
        lib = self.lib("NNN_3_hg38")
        # a 6-nt PAM (e.g. NNGRRT, SaCas9) needs a flank of at least 6
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=1), "")

    def test_wide_pamless_serves_a_shorter_pam(self):
        lib = self.lib("NNNNNNNN_3_hg38")
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")

    def test_equal_width_serves(self):
        lib = self.lib("NNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NNN_3_hg38")

    def test_narrowest_sufficient_flank_wins(self):
        lib = self.lib("NNN_3_hg38", "NNNNNNNN_3_hg38")
        # NRG (3) -> the width-3 one; NNGRRT (6) -> only the width-8 one qualifies
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NNN_3_hg38")
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 3, prune=1), "NNNNNNNN_3_hg38")


def _glob_order(locale):
    """The order bash globs two PAM-less names in under ``locale``."""
    with tempfile.TemporaryDirectory() as d:
        for name in ("NNNNNNNN_3_hg38", "NNN_3_hg38"):
            os.mkdir(os.path.join(d, name))
        return subprocess.run(
            ["bash", "-c", 'cd "$1" && printf "%s\\n" *', "bash", d],
            capture_output=True, text=True,
            env=dict(os.environ, LC_ALL=locale, LANG=locale),
        ).stdout.split()


class AGuideMustFitTheIndex(_Tree):
    """An index name says nothing of the guide it was built for, and the engine searches
    any guide: past what the tree stores it gives wrong alignments or loses targets, at
    exit 0. The resolver reads the length from the .bin, and takes an index only where
    built + N >= guide + bDNA + (bRNA > 0), the rule measured through the real CLI."""

    def built(self, name, guide):
        lib = self.lib(name)
        built_for_guide(lib / name, guide)
        return lib

    def test_the_rule_at_its_edges(self):
        # (built, N, guide, bDNA, bRNA, served): each is a cell measured with the real CLI
        for built, n, guide, bdna, brna, served in (
                (20, 3, 21, 1, 1, True),   # on any N=3 index, byte-identical to a dedicated one
                (20, 3, 21, 2, 2, False),
                (20, 3, 21, 2, 0, True),
                (20, 3, 21, 0, 2, True),
                (20, 3, 21, 2, 1, False),
                (20, 3, 22, 1, 0, True),
                (20, 3, 22, 0, 1, True),
                (20, 3, 22, 1, 1, False),
                (20, 3, 23, 1, 0, False),
                (20, 3, 23, 0, 1, False),
                (20, 2, 21, 1, 1, False),
                (20, 2, 21, 1, 0, True),
                (20, 2, 21, 0, 1, True),
                (20, 2, 23, 1, 1, False),
                (20, 2, 20, 1, 1, True),
                (23, 2, 20, 1, 1, True),   # a longer build holds a shorter guide
                (18, 2, 20, 1, 1, False)):
            with self.subTest(built=built, n=n, guide=guide, bdna=bdna, brna=brna):
                with tempfile.TemporaryDirectory() as tmp:
                    d = Path(tmp) / f"NGG_{n}_hg38"
                    d.mkdir()
                    built_for_guide(d, built)
                    self.assertEqual(
                        resolve(tmp, "NGG", "hg38", max(bdna, brna) + 1, pos=3, guide=guide,
                                bdna=bdna, brna=brna),
                        d.name if served else "")

    def test_it_holds_for_a_5_prime_pam_too(self):
        lib = self.built("TTTV_3_hg38", 20)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4, guide=21, bdna=1, brna=1),
                         "TTTV_3_hg38")
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4, guide=23, bdna=1, brna=1), "")

    def test_it_holds_for_a_pam_less_index_under_prune(self):
        lib = self.built("NNNNNNNN_3_hg38", 20)
        self.assertEqual(resolve(lib, "NNGRRT", "hg38", 2, prune=1, pos=6, guide=21, bdna=1,
                                 brna=1), "NNNNNNNN_3_hg38")  # the shipped 21-nt SaCas9 PAM
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, prune=1, pos=3, guide=23, bdna=1,
                                 brna=1), "")

    def test_a_reference_search_with_no_bulges_opens_no_index_so_any_fits(self):
        lib = self.built("NGG_2_hg38", 20)
        self.assertEqual(resolve(lib, "NGG", "hg38", 1, pos=3, guide=24), "NGG_2_hg38")

    def test_a_variant_search_with_no_bulges_reads_the_index_of_its_indels(self):
        lib = self.built("NGG_2_hg38+ds", 20)
        self.assertEqual(resolve(lib, "NGG", "hg38+ds", 1, pos=3, guide=22), "NGG_2_hg38+ds")
        self.assertEqual(resolve(lib, "NGG", "hg38+ds", 1, pos=3, guide=23), "")

    def test_an_index_whose_length_cannot_be_read_fits_as_before(self):
        lib = self.lib("NGG_2_hg38")  # an empty .bin
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, pos=3, guide=23, bdna=1, brna=1),
                         "NGG_2_hg38")

    def test_a_pam_file_with_no_spacer_fits_as_before(self):
        lib = self.built("NGG_2_hg38", 20)
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, pos=3, bdna=1, brna=1), "NGG_2_hg38")

    def test_the_length_is_read_through_a_symlinked_index(self):
        lib = self.built("NGG_2_hg38", 20)
        shared = Path(self.tmp.name, "shared")
        shared.mkdir()
        (shared / "NGG_2_hg38").symlink_to(lib / "NGG_2_hg38")
        self.assertEqual(resolve(shared, "NGG", "hg38", 2, pos=3, guide=23, bdna=1, brna=1), "")

    def test_a_correct_index_is_not_passed_over_for_a_shallow_one(self):
        # the smallest N used to win: a 20-nt TTTV_2 answered a 23-nt guide beside a 23-nt TTTV_3
        lib = self.built("TTTV_2_hg38", 20)
        (lib / "TTTV_3_hg38").mkdir()
        built_for_guide(lib / "TTTV_3_hg38", 23)
        self.assertEqual(resolve(lib, "TTTV", "hg38", 2, pos=-4, guide=23, bdna=1, brna=1),
                         "TTTV_3_hg38")

    def test_an_index_of_this_very_guide_is_preferred_to_a_longer_one(self):
        # a 23-nt NGG_2 beside a 20-nt NGG_3, searched for 20 nt: both fit, and only one
        # was built for it
        lib = self.built("NGG_2_hg38", 23)
        (lib / "NGG_3_hg38").mkdir()
        built_for_guide(lib / "NGG_3_hg38", 20)
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, pos=3, guide=20, bdna=1, brna=1),
                         "NGG_3_hg38")

    def test_the_requested_pam_still_comes_first(self):
        lib = self.built("NGG_3_hg38", 20)
        (lib / "NNN_3_hg38").mkdir()
        built_for_guide(lib / "NNN_3_hg38", 21)
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, pos=3, guide=21, bdna=1, brna=1),
                         "NGG_3_hg38")

    def test_any_side_finds_one_that_does_not_fit_to_say_why(self):
        lib = self.built("NGG_2_hg38", 20)
        self.assertEqual(resolve(lib, "NGG", "hg38", 2, pos=3, guide=23, bdna=1, brna=1,
                                 any_side="any"), "NGG_2_hg38")

    def test_stderr_stays_silent(self):
        lib = self.built("NGG_2_hg38", 20)
        script = (f"prune_pam=0\npos=3\nfullseqpam={'N' * 26}\nbDNA=1\nbRNA=1\n{RESOLVER}\n"
                  f"_resolve_index_dir \"$1\" NGG hg38 2\n")
        out = subprocess.run(["bash", "-c", script, "bash", str(lib)], capture_output=True,
                             text=True)
        self.assertEqual(out.stderr, "")


class ChoiceIsLocaleIndependent(_Tree):
    """The regression this file exists for.

    Under C (and POSIX) the glob puts the wrong candidate first for both tie rules
    -- NNNNNNNN_3 before NNN_3, NNN_3 before NRG_3 -- so those arms are the demanding
    ones. The en_US.UTF-8 arm adds the opposite order only where that locale is
    installed and really collates differently: glibc silently falls back to C for a
    missing locale, and macOS collates UTF-8 locales by code point.
    """

    LOCALES = ("C", "POSIX", "en_US.UTF-8")

    def test_the_en_us_arm_really_reorders_the_glob(self):
        if _glob_order("en_US.UTF-8") == _glob_order("C"):
            self.skipTest("en_US.UTF-8 is not installed here, or collates like C: the "
                          "locale tests ran under C order only")
        self.assertEqual(_glob_order("C"), ["NNNNNNNN_3_hg38", "NNN_3_hg38"])

    def test_two_pamless_widths_agree_across_locales(self):
        lib = self.lib("NNNNNNNN_3_hg38", "NNN_3_hg38")
        picks = {loc: resolve(lib, "NRG", "hg38", 3, locale=loc) for loc in self.LOCALES}
        self.assertEqual(set(picks.values()), {"NNN_3_hg38"}, picks)

    def test_dedicated_vs_pamless_agrees_across_locales(self):
        lib = self.lib("NNNNNNNN_3_hg38", "NNN_3_hg38", "NRG_3_hg38")
        picks = {loc: resolve(lib, "NRG", "hg38", 3, locale=loc) for loc in self.LOCALES}
        self.assertEqual(set(picks.values()), {"NRG_3_hg38"}, picks)

    def test_full_library_agrees_across_locales(self):
        lib = self.lib(
            "NRG_3_hg38", "NRG_3_hg38_INDELS", "NGG_3_hg38", "NNN_3_hg38",
            "NNNNNNNN_3_hg38", "NNNRRT_3_hg38", "NRG_5_hg38",
        )
        for want, expected in (("NRG", "NRG_3_hg38"),
                               ("NGG", "NGG_3_hg38"),
                               ("NNNRRT", "NNNRRT_3_hg38"),
                               ("NNN", "NNN_3_hg38")):
            picks = {loc: resolve(lib, want, "hg38", 3, locale=loc) for loc in self.LOCALES}
            self.assertEqual(set(picks.values()), {expected}, f"{want}: {picks}")


class AnIndexStillBeingBuiltIsNotSelected(_Tree):
    """A half-written index is indistinguishable from a finished one by NAME.

    `index-genome` creates the directory, takes "<dir>.lock", and then writes one
    .bin per chromosome. A run that starts in that window used to see the name,
    print "Reference index present", and search a PARTIAL index -- returning fewer
    off-targets with no error and exit 0. That is the same failure shape as a
    crashed search reported as a clean negative, so the resolver must not pick it.
    """

    def test_a_locked_index_is_skipped(self):
        lib = self.lib("NRG_3_hg38")
        self.lock(lib, "NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_a_finished_index_is_used_while_another_is_building(self):
        lib = self.lib("NRG_3_hg38", "NRG_4_hg38")
        self.lock(lib, "NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_4_hg38")

    def test_an_unlocked_index_is_still_used(self):
        """The guard must not make a normal library look empty."""
        lib = self.lib("NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_a_directory_with_no_bin_is_treated_as_absent(self):
        lib = self.lib_unfinished("NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_a_bin_less_directory_does_not_shadow_a_real_one(self):
        lib = self.lib("NRG_4_hg38")
        self.lib_unfinished("NRG_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_4_hg38")

    def test_a_pamless_stand_in_needs_a_bin_too(self):
        lib = self.lib_unfinished("NNN_3_hg38")
        self.assertEqual(resolve(lib, "NGG", "hg38", 3, prune=1), "")

    def test_a_lock_without_its_index_is_harmless(self):
        """Locks are named "<dir>.lock"; the glob must not pick one up as an index."""
        lib = self.lib("NRG_3_hg38")
        self.lock(lib, "NRG_9_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")


class ASymlinkedIndexIsStillAnIndex(_Tree):
    """genome_library/<name> may be a symlink to an index stored elsewhere.

    Sharing one copy of a 100 GB index across working directories that way is
    ordinary, and it is how this release's own validation layouts were built. The
    .bin check must look INSIDE the link: GNU find does not follow a symlink given
    as its starting point, so a plain `find "$_d"` saw no .bin and the resolver
    dropped every symlinked index as unbuilt.
    """

    def store(self, name, finished=True):
        d = Path(self.tmp.name) / "store" / name
        d.mkdir(parents=True)
        if finished:
            (d / "chr1.bin").write_bytes(b"")
        return d

    def libdir(self):
        d = Path(self.tmp.name) / "lib"
        d.mkdir(exist_ok=True)
        return d

    def test_a_symlinked_finished_index_is_selected(self):
        lib = self.libdir()
        (lib / "NRG_3_hg38").symlink_to(self.store("NRG_3_hg38"))
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "NRG_3_hg38")

    def test_a_symlinked_pamless_stand_in_is_selected(self):
        lib = self.libdir()
        (lib / "NNN_3_hg38").symlink_to(self.store("NNN_3_hg38"))
        self.assertEqual(resolve(lib, "NGG", "hg38", 3, prune=1), "NNN_3_hg38")

    def test_a_symlink_to_an_unbuilt_index_is_still_absent(self):
        lib = self.libdir()
        (lib / "NRG_3_hg38").symlink_to(self.store("NRG_3_hg38", finished=False))
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")

    def test_a_dangling_symlink_is_absent(self):
        lib = self.libdir()
        (lib / "NRG_3_hg38").symlink_to(Path(self.tmp.name) / "nowhere")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3), "")


class ADedicatedIndexIsTheFastPath(_Tree):
    """An index built for the requested PAM is the one that must answer.

    Measured genome-wide (HGDP+1kGP, 4 mm + 1 DNA + 1 RNA bulge, per-sample) the
    search stage took 44 min on the dedicated NRG_3 index, against 2 h 59 min on
    NNN_3 and 2 h 35 min on NNNNNNNN_3 under --prune-pam, for byte-identical
    reports. So neither installing a PAM-less index beside a dedicated one nor
    turning --prune-pam on may move a search off the dedicated index.

    Every case runs for both suffixes the script resolves: the reference scan
    (<ref>) and the variant scan (<ref>+<vcf>), which is also what the
    precomputed-variant detection looks for.
    """

    SUFFIXES = ("hg38", "hg38+hg38_HGDP1kGP")
    EVERY_KIND = ("NRG_3", "NNN_3", "NNNNNNNN_3")

    def library(self, *prefixes):
        """A fresh library holding <PAM>_<N> for each prefix, under both suffixes."""
        d = Path(tempfile.mkdtemp(dir=self.tmp.name))
        for p in prefixes:
            for s in self.SUFFIXES:
                (d / f"{p}_{s}").mkdir()
                (d / f"{p}_{s}" / "chr1.bin").write_bytes(b"")
        return d

    def assertPicks(self, lib, pam, expected_prefix, prune, pos=3, bmax=(2, 3)):
        for sfx in self.SUFFIXES:
            for want in bmax:
                with self.subTest(suffix=sfx, bmax=want, prune=prune):
                    expected = f"{expected_prefix}_{sfx}" if expected_prefix else ""
                    self.assertEqual(
                        resolve(lib, pam, sfx, want, prune=prune, pos=pos), expected)

    def test_every_kind_installed_without_prune(self):
        # the NNNNNNNN flank is not even eligible here; NNN_3 is, and loses
        self.assertPicks(self.library(*self.EVERY_KIND), "NRG", "NRG_3", prune=0)

    def test_every_kind_installed_with_prune(self):
        # pruning makes the wide flank eligible too; both stand-ins must lose
        self.assertPicks(self.library(*self.EVERY_KIND), "NRG", "NRG_3", prune=1)

    def test_every_kind_installed_with_prune_agrees_across_locales(self):
        lib = self.library(*self.EVERY_KIND)
        for sfx in self.SUFFIXES:
            picks = {loc: resolve(lib, "NRG", sfx, 2, locale=loc, prune=1)
                     for loc in ChoiceIsLocaleIndependent.LOCALES}
            self.assertEqual(set(picks.values()), {f"NRG_3_{sfx}"}, picks)

    def test_dedicated_alone_serves_a_pruned_search(self):
        # --prune-pam over a dedicated index is a no-op filter, not a reason
        # to look for a PAM-less one
        self.assertPicks(self.library("NRG_3"), "NRG", "NRG_3", prune=1)

    def test_pamless_pair_without_a_dedicated_index(self):
        # narrowest flank under pruning; without it only the equal width is
        # eligible -- the same pick either way, in every locale
        lib = self.library("NNN_3", "NNNNNNNN_3")
        for prune in (0, 1):
            self.assertPicks(lib, "NRG", "NNN_3", prune=prune)
            for sfx in self.SUFFIXES:
                picks = {loc: resolve(lib, "NRG", sfx, 2, locale=loc, prune=prune)
                         for loc in ChoiceIsLocaleIndependent.LOCALES}
                self.assertEqual(set(picks.values()), {f"NNN_3_{sfx}"}, picks)

    def test_an_insufficient_dedicated_index_is_never_chosen(self):
        # 2 bulges of each type needs N >= 3; NRG_2 supports one
        for others, prune, expected in (
            (("NNN_3",), 0, "NNN_3"),
            (("NNN_3",), 1, "NNN_3"),
            (("NNNNNNNN_3",), 0, ""),
            (("NNNNNNNN_3",), 1, "NNNNNNNN_3"),
            (("NNN_3", "NNNNNNNN_3"), 1, "NNN_3"),
            ((), 0, ""),
            ((), 1, ""),
        ):
            lib = self.library("NRG_2", *others)
            with self.subTest(installed=("NRG_2",) + others):
                self.assertPicks(lib, "NRG", expected, prune=prune, bmax=(3,))

    def test_a_dedicated_index_returns_once_its_budget_suffices(self):
        # same library as above, one bulge of each type: NRG_2 is sufficient again
        lib = self.library("NRG_2", "NNN_3", "NNNNNNNN_3")
        for prune in (0, 1):
            self.assertPicks(lib, "NRG", "NRG_2", prune=prune, bmax=(2,))

    def test_a_different_dedicated_pam_never_stands_in(self):
        # NGG cannot serve NRG: it never stored the NAG sites
        for others, expected in (((), ""), (("NNN_3",), "NNN_3"), (("NRG_3",), "NRG_3")):
            lib = self.library("NGG_3", *others)
            for prune in (0, 1):
                with self.subTest(installed=("NGG_3",) + others):
                    self.assertPicks(lib, "NRG", expected, prune=prune)

    def test_a_5_prime_pam_only_takes_an_equal_width_flank(self):
        # TTTV (Cas12a, PAM-file offset -4): a wider flank is refused by the
        # engine, a narrower one never stored the PAM. Every all-N index here is a
        # 5' build, so its width alone decides.
        for installed, expected in (
            (("NNNNNNNN_3",), ""),
            (("NNN_3",), ""),
            (("NNNN_3", "NNNNNNNN_3"), "NNNN_3"),
            (("TTTV_3", "NNNN_3", "NNNNNNNN_3"), "TTTV_3"),
        ):
            lib = self.library(*installed)
            for d in lib.iterdir():
                flank = d.name.split("_")[0]
                if set(flank) == {"N"}:
                    built_for(d, -len(flank))
            for prune in (0, 1):
                with self.subTest(installed=installed):
                    self.assertPicks(lib, "TTTV", expected, prune=prune, pos=-4)


# Stand-in for crispritz.py. It records every call and does nothing else, so the
# script's own scans make every decision.
_STUB_ENGINE = r"""#!/bin/bash
# Accepts --prune-pam: the script refuses to prune with an engine whose source
# never mentions the flag, so this line is load-bearing.
printf '%s\n' "$*" >> "$STUB_LOG"
# STUB_HOLD_BUILD=1 keeps an index build running, as a real one does for hours,
# until the test signals the run.
if [ "$1" = index-genome ] && [ "$STUB_HOLD_BUILD" = 1 ]; then
	exec sleep 300
fi
# STUB_BUILD_UNTIL / STUB_SEARCH_UNTIL=<file> hold an index build / a search until the
# test creates <file>, which ends it as usual.
if [ "$1" = index-genome ] && [ -n "$STUB_BUILD_UNTIL" ]; then
	while ! [ -e "$STUB_BUILD_UNTIL" ]; do sleep 0.05; done
fi
# STUB_BUILD_FAILS=1: the index build then fails, saying so on stderr as crispritz does.
if [ "$1" = index-genome ] && [ "$STUB_BUILD_FAILS" = 1 ]; then
	echo "ERROR: buildTST (index build) failed (exit code 1)." >&2
	exit 1
fi
# STUB_BUILD_WRITES=1: index-genome <name> <genome> <PAM file> -bMax <N> leaves what
# crispritz's does in the current directory, genome_library/<PAM>_<N>_<name>, made
# anew, with a .bin in it.
if [ "$1" = index-genome ] && [ "$STUB_BUILD_WRITES" = 1 ]; then
	_seq=$(cut -f1 -d' ' "$4") _pos=$(cut -f2 -d' ' "$4")
	if [ "$_pos" -gt 0 ]; then _pam=${_seq:${#_seq}-$_pos}; else _pam=${_seq:0:-$_pos}; fi
	_out="genome_library/${_pam}_${6}_$2"
	rm -rf "$_out" && mkdir -p "$_out" && : >"$_out/chr1.bin"
fi
if [ "$1" = search ] && [ -n "$STUB_SEARCH_UNTIL" ]; then
	while ! [ -e "$STUB_SEARCH_UNTIL" ]; do sleep 0.05; done
fi
# STUB_ENRICH_UNTIL=<file> holds add-variants the same way; STUB_ENRICH_FAILS=1 then
# fails it, on stderr, before it writes anything.
if [ "$1" = add-variants ] && [ -n "$STUB_ENRICH_UNTIL" ]; then
	while ! [ -e "$STUB_ENRICH_UNTIL" ]; do sleep 0.05; done
fi
if [ "$1" = add-variants ] && [ "$STUB_ENRICH_FAILS" = 1 ]; then
	echo "ERROR: add-variants failed" >&2
	exit 1
fi
# STUB_ENRICH=1: add-variants leaves in the current directory what crispritz's does for
# VCF folder $2 and reference folder $3 -- the enriched genome, its dictionaries and
# indel logs, one fake-indel genome per chromosome. Otherwise it writes nothing, as
# with no source VCFs installed.
# STUB_ENRICH=snp: the same, for chromosome chr1 of 120 bases, as the enricher writes it: a
# SNP as its IUPAC code at offset 10 of the enriched genome, and an insertion after offset
# 25, cut out with its flanks as a fake contig, whose log row maps it back (what the SNP
# overlay reads).
if [ "$1" = add-variants ] && [ "$STUB_ENRICH" = snp ]; then
	_ref=${3%/}; _ref=${_ref##*/}; _vcf=${2%/}; _vcf=${_vcf##*/}
	_sn=variants_genome/SNPs_genome
	mkdir -p "$_sn/${_ref}_enriched" "variants_genome/fake_${_vcf}_chr1"
	_seq=ACGTTGCAAGACGTTGCAAGACGTTGCAAGACGTTGCAAGACGTTGCAAGACGTTGCAAG
	printf '>chr1\n%sR%s\n' "${_seq:0:10}" "${_seq:11}" >"$_sn/${_ref}_enriched/chr1.enriched.fa"
	printf '{}\n' >"$_sn/my_dict_chr1.json"
	printf 'CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\nchr1_0-60_1\tS2\trs1\t0.25\tchr1_26_%s_%sT\t0,61\t%s\n' \
		"${_seq:25:1}" "${_seq:25:1}" "$_seq" >"$_sn/logchr1.txt"
	printf '>fakechr1\n%sT%s\n' "${_seq:0:26}" "${_seq:26}" >"variants_genome/fake_${_vcf}_chr1/fakechr1.fa"
fi
if [ "$1" = add-variants ] && [ "$STUB_ENRICH" = 1 ]; then
	_ref=${3%/}; _ref=${_ref##*/}; _vcf=${2%/}; _vcf=${_vcf##*/}
	mkdir -p "variants_genome/SNPs_genome/${_ref}_enriched" "variants_genome/fake_${_vcf}_chr1"
	printf '>chr1\nACGT\n' >"variants_genome/SNPs_genome/${_ref}_enriched/chr1.fa"
	printf '{}\n' >"variants_genome/SNPs_genome/my_dict_chr1.json"
	: >"variants_genome/SNPs_genome/logchr1.txt"
	printf '>fakechr1\nACGTACGTACGTACGTACGTACGTACGT\n' >"variants_genome/fake_${_vcf}_chr1/fakechr1.fa"
fi
if [ "$1" = search ]; then
	# STUB_ENGINE_STOPS=reference|variant: that search says what stops it on stdout and
	# exits 0, with no targets file, as crispritz.py does for a folder that is not there.
	_which=reference
	case "$2" in *+*) _which=variant ;; esac  # (an index of <ref>+<vcf>, or that enriched genome)
	if [ "$STUB_ENGINE_STOPS" = "$_which" ]; then
		echo "ERROR! $2 is not a directory."
		exit 0
	fi
	# STUB_FINDS_NOTHING=1: every search writes the targets file of one that found
	# nothing, the header alone, and the run goes on past them.
	if [ "$STUB_FINDS_NOTHING" = 1 ]; then
		printf '#Bulge type\tcrRNA\n' >"$5.targets.txt"
		exit 0
	fi
	# Stop at the last search: the variant one, or the reference one in a
	# reference-only run. Every index decision has been made by then.
	[ "$_which" = variant ] && exit 1
	[ "$STUB_LAST" = reference ] && exit 1
	printf '#Bulge type\tcrRNA\n' >"$5.targets.txt"
fi
exit 0
"""


class _Run:
    """What one pipeline run decided, read off its stdout and the stub's log."""

    def __init__(self, stdout, calls, stderr):
        self.stdout = stdout
        self.stderr = stderr
        self.calls = calls
        self.precomputed = "Using precomputed variant index" in stdout
        self.reference = self._present("Reference index present")
        self.variant = self._present("Variant index present")
        self.indels = self._present("Indels index present")
        self.enriched = any(c[0] == "add-variants" for c in calls)
        self.builds = [c for c in calls if c[0] == "index-genome"]

    def _present(self, what):
        m = re.search(re.escape(what) + r" \((.+?)\)", self.stdout)
        return m.group(1) if m else None

    def searched(self, variant):
        """(index basename, --prune-pam passed?) for one search, or None if it never ran.

        The script may hand the engine a .bin-only view of the index; the view keeps
        the index's basename.
        """
        for c in self.calls:
            if c[0] == "search" and ("-var" in c) == variant:
                return Path(c[1]).name, "--prune-pam" in c
        return None


def _gnu_realpath():
    """True when `realpath` resolves a path that does not exist, as GNU's does.

    The script runs `realpath` on the "_" placeholder that complete-search passes for
    a reference-only search; BSD realpath (macOS) fails on it and writes to stderr,
    which the script then reads as a failed stage. The script targets the Linux
    image, so these tests run where it does.
    """
    try:
        probe = subprocess.run(["realpath", "_no_such_path_"], capture_output=True)
    except OSError:
        return False
    return probe.returncode == 0


class _Launched:
    """A run of the shipped script started by `_ShippedScript.launch`."""

    def __init__(self, proc, work):
        self.proc = proc
        self.work = work

    def stdout(self):
        return (self.work / "stdout.txt").read_text()

    def stderr(self):
        return (self.work / "out" / "log_error.txt").read_text()

    def calls(self):
        log = self.work / "engine_calls.txt"
        return [line.split() for line in log.read_text().splitlines()] if log.exists() else []


def _reap(proc):
    """Kill whatever is left of a run a test started (it runs in its own session)."""
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


@unittest.skipUnless(shutil.which("python") and _gnu_realpath(),
                     "the shipped script needs python and GNU realpath on PATH")
class _ShippedScript(unittest.TestCase):
    """Runs the shipped script with a stub engine on PATH. The layout is a
    downloaded install: indexes and Dictionaries/ tiers, no VCFs."""

    VCF = "hg38_HGDP1kGP"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.home = self.root / "home"  # the working directory crisprme runs from
        self.library = self.home / "genome_library"
        self.library.mkdir(parents=True)
        genome = self.home / "Genomes" / "hg38"
        genome.mkdir(parents=True)
        (genome / "chr1.fa").write_text(">chr1\nACGT\n")
        dicts = self.home / "Dictionaries"
        (dicts / f"registry_{self.VCF}").mkdir(parents=True)
        (dicts / f"log_indels_{self.VCF}").mkdir()
        (dicts / f"log_indels_{self.VCF}" / "logchr1.txt").write_text("")
        stub_dir = self.root / "bin"
        stub_dir.mkdir()
        (stub_dir / "crispritz.py").write_text(_STUB_ENGINE)
        (stub_dir / "crispritz.py").chmod(0o755)
        self.path = f"{stub_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        self.runs = 0

    def install(self, *prefixes, reference=True, variant=True, library=None):
        """Per <PAM>_<N>: a finished reference index and/or a finished variant index
        with its INDELS companion."""
        library = library or self.library
        library.mkdir(parents=True, exist_ok=True)
        for p in prefixes:
            names = ([f"{p}_hg38"] if reference else []) + (
                [f"{p}_hg38+{self.VCF}"] if variant else [])
            for name in names:
                (library / name).mkdir()
                (library / name / "chr1.bin").write_bytes(b"")
            if variant:
                (library / f"{p}_hg38+{self.VCF}_INDELS").mkdir()

    def run_pipeline(self, pam="20bp-NRG-SpCas9", bulges=1, prune=0, variant=True,
                     index_path=None, threads=1, env=None, timeout=300, samples=None,
                     postprocess=None, pam_text=None):
        """Run the script the way complete-search launches it, stopping at its last
        search."""
        return self.finish(self.launch(pam, bulges, prune, variant, index_path,
                                       threads, env, samples=samples,
                                       postprocess=postprocess,
                                       pam_text=pam_text), timeout)

    def finish(self, run, timeout=300):
        """Wait for a launched run; a run still going after `timeout` is a failure."""
        try:
            run.proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            _reap(run.proc)
            self.fail(f"the run had not finished after {timeout} s -- still waiting "
                      f"on a lock?\n{run.stdout()}")
        return _Run(run.stdout(), run.calls(), run.stderr())

    def launch(self, pam="20bp-NRG-SpCas9", bulges=1, prune=0, variant=True,
               index_path=None, threads=1, env=None, wrap=(), samples=None,
               postprocess=None, pam_text=None):
        """Start the script the way complete-search launches it (argument order of
        crisprme.py's crisprme_run), in a session of its own, so that a test can
        signal the whole run as Ctrl-C or a container stop does. `env` adds to its
        environment; the RAM thread guard is off unless `env` turns it on. `wrap` is a
        command that runs the script's bash as its child. `samples` is the samplefile
        to hand it in place of the run's own. `postprocess` is the PostProcess folder
        to run it from, which complete-search also passes as its folder (the helpers
        the script runs, and the samples placeholder); the shipped one by default.
        `pam_text`: the PAM file <pam>.txt holds it, for a PAM file none of the shipped
        ones is."""
        postprocess = Path(postprocess or _SCRIPT.parent)
        if variant and (env or {}).get("STUB_ENRICH"):
            # (add-variants has VCFs to enrich from: the script stops first without any)
            folder = self.home / "VCFs" / self.VCF
            folder.mkdir(parents=True, exist_ok=True)
            if not any(folder.iterdir()):
                (folder / "x.chr1.vcf.gz").write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n"))
        self.runs += 1
        work = self.root / f"run{self.runs}"
        out = work / "out"
        out.mkdir(parents=True)
        pam_file = work / f"{pam}.txt"
        if pam_text is None:
            shutil.copy(_PAMS / f"{pam}.txt", pam_file)
        else:
            pam_file.write_text(pam_text)
        guides = work / "guides.txt"
        guides.write_text("N" * len(pam_file.read_text().split()[0]) + "\n")
        empty = work / "empty.txt"
        empty.write_text("")
        given = samples
        vcfs = samples = "_"
        if variant:
            vcfs, samples = work / "vcfs.txt", work / "samples.txt"
            vcfs.write_text(self.VCF + "\n")
            samples.write_text("")
        samples = given or samples
        # complete-search refuses --prune-pam without --index-path; the usual value
        # is the working directory's own genome_library
        if index_path is None:
            index_path = self.library if prune else "_"
        args = [self.home / "Genomes" / "hg38", vcfs, guides, pam_file, empty, samples,
                bulges + 1, 4, bulges, bulges, 3, out, postprocess, threads, self.home,
                empty, "_", 1, 0, "_", "mm+bulges", "mm+bulges,mm", "False", "PASS,.",
                index_path, 6, prune]
        run_env = dict(os.environ, PATH=self.path, STUB_LOG=str(work / "engine_calls.txt"),
                       STUB_LAST="variant" if variant else "reference",
                       LC_ALL="C", LANG="C", CRISPRME_NO_THREAD_GUARD="1")
        run_env.update(env or {})
        # stderr goes where complete-search sends it: the script fails a stage
        # whenever log_error.txt is non-empty
        with open(out / "log_error.txt", "w") as err, open(work / "stdout.txt", "w") as so:
            proc = subprocess.Popen([*wrap, "bash", str(postprocess / _SCRIPT.name),
                                     *map(str, args)], cwd=work,
                                    env=run_env, stdout=so, stderr=err, text=True,
                                    start_new_session=True)
        self.addCleanup(_reap, proc)
        return _Launched(proc, work)


class EveryScanInTheShippedScript(_ShippedScript):
    """The three scans, run by the shipped script itself with a stub engine.

    The resolver tests above prove the rule. These prove that each call site --
    the precomputed-variant detection, the reference scan and the variant scan --
    applies it with the right library, suffix, bulge budget and pruning state; that
    the reference and variant searches are handed the index their scan picked; that
    the script announces the picked variant index's own INDELS companion; and that
    --prune-pam reaches those two searches exactly when asked. The stub stops the
    run at the variant search, so the indel search itself (pool_search_indels.py,
    which is handed the companion and the flag) is not reached here.
    """

    def assertServedBy(self, run, prefix, prune):
        """Every scan picked <prefix>, and both searches were handed it."""
        ref, var = f"{prefix}_hg38", f"{prefix}_hg38+{self.VCF}"
        self.assertTrue(run.precomputed, run.stdout)
        self.assertFalse(run.enriched)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.reference, ref)
        self.assertEqual(run.variant, var)
        self.assertEqual(run.indels, f"{var}_INDELS")
        self.assertEqual(run.searched(variant=False), (ref, bool(prune)))
        self.assertEqual(run.searched(variant=True), (var, bool(prune)))
        # nothing failed before the variant search, where the stub stops the run
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed",
                         run.stderr)

    def assertNoVariantIndex(self, run):
        """No installed variant index qualified, and nothing here can build one: the run
        says so before it does anything, rather than fail in an enrichment with no VCFs
        ('mv: cannot stat')."""
        self.assertFalse(run.precomputed, run.stdout)
        self.assertFalse(run.enriched)
        self.assertIsNone(run.searched(variant=False))
        self.assertIsNone(run.searched(variant=True))
        self.assertIn("ERROR: no variant index here serves this search", run.stderr)
        self.assertNotIn("Genome enrichment failed", run.stderr)
        self.assertNotIn("cannot stat", run.stderr)

    def test_every_kind_installed_without_prune(self):
        self.install("NRG_3", "NNN_3", "NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(prune=0), "NRG_3", prune=0)

    def test_every_kind_installed_with_prune(self):
        self.install("NRG_3", "NNN_3", "NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(prune=1), "NRG_3", prune=1)

    def test_dedicated_alone_with_prune_is_searched_not_refused(self):
        self.install("NRG_3")
        self.assertServedBy(self.run_pipeline(prune=1), "NRG_3", prune=1)

    def test_pamless_pair_without_a_dedicated_index(self):
        self.install("NNN_3", "NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(prune=1), "NNN_3", prune=1)
        self.assertServedBy(self.run_pipeline(prune=0), "NNN_3", prune=0)

    def test_a_wide_flank_alone_is_searched_pruned(self):
        # without --prune-pam the search prunes by itself (2.9.0): it stopped for want of
        # a variant index, and said to pass the flag
        self.install("NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(prune=1), "NNNNNNNN_3", prune=1)
        run = self.run_pipeline(prune=0)
        self.assertServedBy(run, "NNNNNNNN_3", prune=1)
        self.assertIn("searching with --prune-pam, as if it was given", run.stdout)

    def test_a_dedicated_index_at_a_larger_budget_beats_a_tighter_stand_in(self):
        self.install("NRG_4", "NNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(prune=prune), "NRG_4", prune)

    def test_precomputed_detection_looks_at_the_variant_index(self):
        # a variant index with no reference index beside it: the reference one is
        # built from the raw genome, the variant one is used as installed
        self.install("NRG_3", reference=False)
        run = self.run_pipeline(prune=0)
        self.assertTrue(run.precomputed, run.stdout)
        self.assertFalse(run.enriched)
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))  # (built)
        self.assertEqual(run.searched(variant=True), (f"NRG_3_hg38+{self.VCF}", False))

    def test_a_reference_index_alone_is_not_a_precomputed_variant_index(self):
        self.install("NRG_3", variant=False)
        self.assertNoVariantIndex(self.run_pipeline(prune=0))

    def test_index_path_supplies_the_reference_index_only(self):
        # --index-path is the REFERENCE library; the variant index is always read
        # from the working directory's genome_library
        shared = self.root / "shared_library"
        self.install("NRG_3", "NNN_3", variant=False, library=shared)
        self.install("NRG_3", "NNN_3", reference=False)
        for prune in (0, 1):
            run = self.run_pipeline(prune=prune, index_path=shared)
            self.assertTrue(run.precomputed, run.stdout)
            self.assertEqual(run.builds, [])
            self.assertEqual(run.searched(variant=False), ("NRG_3_hg38", bool(prune)))
            self.assertEqual(run.searched(variant=True),
                             (f"NRG_3_hg38+{self.VCF}", bool(prune)))

    def test_a_variant_index_at_index_path_is_searched_there(self):
        # It used to be looked for in the working directory alone: without VCFs the search
        # stopped, and with them it enriched the genome and built its own (OF-5).
        shared = self.root / "shared_library"
        self.install("NRG_3", library=shared)
        for prune in (0, 1):
            run = self.run_pipeline(prune=prune, index_path=shared)
            self.assertServedBy(run, "NRG_3", prune)
            var = [c for c in run.calls if c[0] == "search" and "-var" in c][0]
            self.assertEqual(var[1], str(shared / f"NRG_3_hg38+{self.VCF}"))
            self.assertIn(f"Variant index library: {shared}", run.stdout)
        # without --index-path, the working directory's library alone, as before
        self.assertNoVariantIndex(self.run_pipeline())

    def test_the_index_path_library_comes_first(self):
        shared = self.root / "shared_library"
        self.install("NRG_3", library=shared)
        self.install("NRG_3", reference=False)
        run = self.run_pipeline(index_path=shared)
        var = [c for c in run.calls if c[0] == "search" and "-var" in c][0]
        self.assertEqual(var[1], str(shared / f"NRG_3_hg38+{self.VCF}"))

    def test_an_indels_index_missing_beside_a_variant_index_at_index_path_is_refused(self):
        shared = self.root / "shared_library"
        self.install("NRG_3", library=shared)
        shutil.rmtree(shared / f"NRG_3_hg38+{self.VCF}_INDELS")
        (self.home / "Genomes" / f"hg38+{self.VCF}_INDELS").mkdir()  # (a search could build it)
        run = self.run_pipeline(index_path=shared)
        self.assertEqual(run.builds, [])
        self.assertIsNone(run.searched(variant=False))
        self.assertIn(f"ERROR: INDELS index 'NRG_3_hg38+{self.VCF}_INDELS' is not beside its "
                      f"variant index in '{shared}' (--index-path), and a search builds no index "
                      "there", run.stderr)
        self.assertFalse((self.library / f"NRG_3_hg38+{self.VCF}_INDELS").exists())

    def test_an_insufficient_dedicated_index_is_never_chosen(self):
        self.install("NRG_2", "NNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(bulges=2, prune=prune), "NNN_3", prune)
        # ... and it is the one used again for a search it can serve
        self.assertServedBy(self.run_pipeline(bulges=1, prune=0), "NRG_2", prune=0)

    def test_an_insufficient_dedicated_index_alone_is_not_used(self):
        self.install("NRG_2")
        self.assertNoVariantIndex(self.run_pipeline(bulges=2, prune=0))
        ref_only = self.run_pipeline(bulges=2, prune=1, variant=False)
        self.assertIsNone(ref_only.searched(variant=False))
        self.assertIn("no matching reference index", ref_only.stderr)

    def test_a_different_dedicated_pam_never_stands_in(self):
        self.install("NGG_3")
        self.assertNoVariantIndex(self.run_pipeline(prune=1))
        # reference-only: under --index-path the missing index is an error ...
        pruned = self.run_pipeline(prune=1, variant=False)
        self.assertIsNone(pruned.searched(variant=False))
        self.assertIn("no matching reference index", pruned.stderr)
        # ... otherwise a dedicated NRG index is built and searched, never NGG_3
        plain = self.run_pipeline(prune=0, variant=False)
        self.assertEqual([c[1] for c in plain.builds], ["hg38"])
        self.assertEqual(plain.searched(variant=False), ("20-NRG_2_hg38", False))

    def test_a_different_dedicated_pam_beside_a_stand_in(self):
        self.install("NGG_3", "NNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(prune=prune), "NNN_3", prune)

    def built_for(self, prefix, pos):
        """Record the side <prefix>'s reference and variant indexes were built for."""
        for d in self.library.glob(f"{prefix}_hg38*"):
            if not d.name.endswith("_INDELS"):
                built_for(d, pos)

    def test_a_5_prime_pam_with_a_wider_flank_only(self):
        self.install("NNNNNNNN_3")
        self.built_for("NNNNNNNN_3", -8)  # a 5' build: only its width refuses it
        self.assertNoVariantIndex(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1))
        ref_only = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1, variant=False)
        self.assertIsNone(ref_only.searched(variant=False))
        self.assertIn("no matching reference index", ref_only.stderr)

    def test_a_5_prime_pam_takes_the_equal_width_flank(self):
        self.install("NNNN_3", "NNNNNNNN_3")
        self.built_for("NNNN_3", -4)
        self.built_for("NNNNNNNN_3", -8)
        self.assertServedBy(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1),
                            "NNNN_3", prune=1)

    def test_a_5_prime_pam_never_takes_a_3_prime_build(self):
        # It used to, and the engine then found nothing, with exit 0.
        self.install("NNNN_3")
        self.built_for("NNNN_3", 4)
        self.assertNoVariantIndex(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1))
        ref_only = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1, variant=False)
        self.assertIsNone(ref_only.searched(variant=False))
        self.assertIn("no matching reference index", ref_only.stderr)
        # without --index-path the search builds a dedicated index and searches it
        plain = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=0, variant=False)
        self.assertEqual([c[1] for c in plain.builds], ["hg38"])
        self.assertEqual(plain.searched(variant=False), ("TTTV-23_2_hg38", False))

    def test_a_5_prime_pam_prefers_its_dedicated_index(self):
        self.install("TTTV_3", "NNNN_3", "NNNNNNNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=prune),
                                "TTTV_3", prune)


# All-N PAM files of a 20-nt guide, as a user writes one: none ships for a 5' PAM.
_PAM_5P_4 = ("27bp-NNNN-5p", "N" * 27 + " -4\n")
_PAM_5P_3 = ("23bp-NNN-5p", "N" * 23 + " -3\n")


class APamLessIndexKeepsTheSideItWasBuiltFor(_ShippedScript):
    """The side of the guide an all-N index was built for is known from its name since
    2.9.0, and else only from its .pam_build (_index_side). A search that built one wrote
    none, so a 5' PAM-less index a search had built was taken for a 3' build and refused by
    every later 5' search under --index-path, which --prune-pam needs: TTTV with
    --prune-pam on the NNNN_2 a 27bp-NNNN-5p search had built stopped with "no matching
    reference index", where the resolver before the side check took it and found the
    dedicated index's 3 rows. And an all-N search PAM built its index under the name of an
    all-N index, so one in place under that name, built for the other side, was searched as
    if just built: no target, exit 0; 2.8.1 refused it as in the way. A search now records
    the side of an index it builds, builds it under a name that says its guide length and
    its side (_build_name), beside an index of the older name built for the other side,
    which it neither searches nor touches, and says why an all-N index under --index-path
    was not used."""

    def pam_build(self, name):
        return (self.library / name / ".pam_build").read_text()

    def test_a_search_records_the_side_of_an_all_n_index_it_builds(self):
        pam, text = _PAM_5P_4
        run = self.run_pipeline(pam=pam, pam_text=text, variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("NNNN-23_2_hg38", False))
        self.assertEqual(self.pam_build("NNNN-23_2_hg38"), "27bp-NNNN-5p 27 -4\n")
        # and a 5' search under --index-path takes it, as it took one before
        for prune in (0, 1):
            run = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=prune, variant=False,
                                    index_path=self.library)
            self.assertEqual(run.builds, [])
            self.assertEqual(run.searched(variant=False), ("NNNN-23_2_hg38", bool(prune)))
            self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed")

    def test_a_3_prime_all_n_build_records_its_side_and_a_dedicated_one_too(self):
        # (a dedicated index was named alike by PAM files of both sides: cli5 D2)
        run = self.run_pipeline(pam="20bp-NRG-SpCas9", variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))
        self.assertEqual(self.pam_build("20-NRG_2_hg38"), "20bp-NRG-SpCas9 23 3\n")
        run = self.run_pipeline(pam="20bp-NNN-NO-PAM", variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("20-NNN_2_hg38", False))
        self.assertEqual(self.pam_build("20-NNN_2_hg38"), "20bp-NNN-NO-PAM 23 3\n")

    def test_a_variant_build_records_its_side_too(self):
        pam, text = _PAM_5P_4
        run = self.run_pipeline(pam=pam, pam_text=text,
                                env={"STUB_BUILD_WRITES": "1", "STUB_ENRICH": "1"})
        self.assertTrue(run.enriched)
        self.assertEqual([c[1] for c in run.builds][:2], ["hg38", f"hg38+{self.VCF}"])
        self.assertEqual(run.searched(variant=True), (f"NNNN-23_2_hg38+{self.VCF}", False))
        self.assertEqual(self.pam_build(f"NNNN-23_2_hg38+{self.VCF}"), "27bp-NNNN-5p 27 -4\n")
        # its INDELS index under its name, the engine's renamed: nothing under the older one
        self.assertTrue((self.library / f"NNNN-23_2_hg38+{self.VCF}_INDELS").is_dir())
        self.assertEqual(sorted(p.name for p in self.library.iterdir()),
                         ["NNNN-23_2_hg38", f"NNNN-23_2_hg38+{self.VCF}",
                          f"NNNN-23_2_hg38+{self.VCF}_INDELS"])

    def assertBuiltBeside(self, run, index, built, pam_file, variant=False, record=None):
        """The run built and searched its own index `built`, under the name of 2.9.0,
        beside `index`, an index of the older name for the other side, which it neither
        searched nor changed (its .pam_build stays `record`): 2.8.1 stopped on it as in the
        way of the index it built under that name."""
        self.assertIn(f"hg38+{self.VCF}" if variant else "hg38", [c[1] for c in run.builds])
        self.assertEqual(run.searched(variant=variant), (built, False))
        self.assertNotIn("is in the way", run.stderr)
        self.assertEqual(self.pam_build(built).split()[0], pam_file)
        if record is not None:
            self.assertEqual(self.pam_build(index), record)

    def test_an_all_n_search_builds_its_own_beside_one_built_for_the_other_side(self):
        # It once searched it: the build path found the name in place. 2.8.1 refused it.
        self.install("NNN_2", variant=False)
        built_for(self.library / "NNN_2_hg38", -3)
        run = self.run_pipeline(pam="20bp-NNN-NO-PAM", variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertBuiltBeside(run, "NNN_2_hg38", "20-NNN_2_hg38", "20bp-NNN-NO-PAM",
                               record="23bp-NNN-NO-PAM 23 -3\n")
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed")

    def test_and_a_5_prime_one_beside_a_3_prime_build_without_a_record(self):
        self.install("NNN_2", variant=False)
        pam, text = _PAM_5P_3
        run = self.run_pipeline(pam=pam, pam_text=text, variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertBuiltBeside(run, "NNN_2_hg38", "NNN-20_2_hg38", "23bp-NNN-5p")
        self.assertFalse((self.library / "NNN_2_hg38" / ".pam_build").exists())
        # one in place recorded for its side is searched, and one of the 2.9.0 name ahead of it
        built_for(self.library / "NNN_2_hg38", -3)
        shutil.rmtree(self.library / "NNN-20_2_hg38")
        run = self.run_pipeline(pam=pam, pam_text=text, variant=False)
        self.assertEqual(run.reference, "NNN_2_hg38")
        self.assertEqual(run.searched(variant=False), ("NNN_2_hg38", False))

    def test_the_variant_build_path_builds_its_own_beside_one_too(self):
        # the reference index serves; the variant one in place was built for the 5' side
        self.install("NNN_2")
        built_for(self.library / "NNN_2_hg38", 3)
        built_for(self.library / f"NNN_2_hg38+{self.VCF}", -3)
        run = self.run_pipeline(pam="20bp-NNN-NO-PAM",
                                env={"STUB_ENRICH": "1", "STUB_BUILD_WRITES": "1"})
        self.assertEqual(run.reference, "NNN_2_hg38")
        self.assertBuiltBeside(run, f"NNN_2_hg38+{self.VCF}", f"20-NNN_2_hg38+{self.VCF}",
                               "20bp-NNN-NO-PAM", variant=True, record="23bp-NNN-NO-PAM 23 -3\n")
        self.assertEqual(run.indels, None)  # (built, not found in place)
        self.assertTrue((self.library / f"20-NNN_2_hg38+{self.VCF}_INDELS").is_dir())
        self.assertTrue((self.library / f"NNN_2_hg38+{self.VCF}_INDELS").is_dir())  # its own

    def test_an_index_this_run_may_use_is_not_refused(self):
        # the same side, and a dedicated PAM, whose index is never all-N
        self.install("NNN_2", variant=False)
        built_for(self.library / "NNN_2_hg38", 3)
        run = self.run_pipeline(pam="20bp-NNN-NO-PAM", variant=False)
        self.assertEqual(run.reference, "NNN_2_hg38")
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed")
        run = self.run_pipeline(pam="23bp-TTTV-Cas12a", variant=False)
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("TTTV-23_2_hg38", False))

    def assertNoMatchSaysWhy(self, run, note, pam, n=2, guide=20, side=3):
        """Under --index-path: the error, why the all-N index there was not used
        (`note`, if any), and what would serve: the published name (the indexes published
        are of N = 3), and the one a build gives it since 2.9.0, of a `guide`-nt guide on
        the `side` of the guide."""
        lib = os.path.realpath(self.library)
        named = f"{guide}-{pam}" if side == 3 else f"{pam}-{guide}"
        published = f"{pam}_{max(n, 3)}_hg38"
        self.assertIsNone(run.searched(variant=False))
        self.assertEqual(run.stderr,
                         f"ERROR: no matching reference index under --index-path '{lib}'\n"
                         + note +
                         f"       expected an index '{pam}_N_hg38' with N >= {n} "
                         f"(e.g. {published}, or {named}_{n}_hg38 as a build names it "
                         "since 2.9.0)\n"
                         f"       re-run 'crisprme.py download --what index --index-name "
                         f"{published}' if it is published (it fetches the shared reference "
                         "index too), or build it with 'crisprme.py build-index-only'.\n")

    def test_index_path_says_why_an_all_n_index_there_was_not_used(self):
        self.install("NNN_2", variant=False)
        built_for(self.library / "NNN_2_hg38", -3)
        for prune in (0, 1):
            self.assertNoMatchSaysWhy(
                self.run_pipeline(pam="20bp-NNN-NO-PAM", prune=prune, variant=False,
                                  index_path=self.library),
                "       NNN_2_hg38 is a PAM-less index for a PAM on the 5' side of the guide "
                "(as its .pam_build records), and 20bp-NNN-NO-PAM.txt puts the PAM on the "
                "3' side, so a search of it would find nothing\n", "NNN")

    def test_and_an_index_without_a_record_is_named_with_how_to_record_its_side(self):
        # a 5' PAM-less index a search built before this release: named for its side, a
        # name records it, as for a 3' one (2.8.1 said to record it with build-index-only,
        # which still does)
        self.install("NNNN_2", variant=False)
        lib = os.path.realpath(self.library)
        self.assertNoMatchSaysWhy(
            self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1, variant=False),
            "       NNNN_2_hg38 is a PAM-less index for a PAM on the 3' side of the guide "
            "(it records no side in a .pam_build, and an index that does not is taken for "
            "a 3' build), and 23bp-TTTV-Cas12a.txt puts the PAM on the 5' side, so a search "
            "of it would find nothing\n"
            f"       if it was built from a 5' PAM file, name it so: mv '{lib}/NNNN_2_hg38' "
            f"'{lib}/NNNN-<guide length>_2_hg38', which records that\n", "TTTV",
            guide=23, side=5)

    def test_a_dedicated_index_of_the_other_side_has_one_built_beside_it(self):
        # a 3' build of TTTV under the name a 5' TTTV search built: it was searched, and
        # reported a site that is in no chromosome, at exit 0 (cli5 D2); 2.8.1 refused it
        self.install("TTTV_2", variant=False)
        Path(self.library, "TTTV_2_hg38", ".pam_build").write_text("20bp-TTTV-3p 24 4\n")
        run = self.run_pipeline(pam="23bp-TTTV-Cas12a", variant=False,
                                env={"STUB_BUILD_WRITES": "1"})
        self.assertBuiltBeside(run, "TTTV_2_hg38", "TTTV-23_2_hg38", "23bp-TTTV-Cas12a",
                               record="20bp-TTTV-3p 24 4\n")

    def test_an_unrecorded_all_n_index_of_another_width_is_not_taken_for_a_3_prime_pam(self):
        # flows OF-8: the NNNN_2 a 5' search of 2.8.0 left, then a 3' NGCG search
        self.install("NNNN_2", variant=False)
        plain = self.run_pipeline(pam="20bp-NGCG-SpCas9VRER", variant=False)
        self.assertEqual([c[1] for c in plain.builds], ["hg38"])
        self.assertEqual(plain.searched(variant=False), ("20-NGCG_2_hg38", False))
        lib = os.path.realpath(self.library)
        # (and, since 2.9.0, how to use one built for a 3' PAM: a name that says its side,
        # here for the guide length its trees say, which these stand-ins do not)
        self.assertNoMatchSaysWhy(
            self.run_pipeline(pam="20bp-NGCG-SpCas9VRER", prune=1, variant=False),
            "       NNNN_2_hg38 is a PAM-less index that records no side of the guide in a "
            ".pam_build, and a search before 2.8.1 built one of that width from a 5' PAM file as "
            "well as from a 3' one, so it is not taken for a 3' PAM\n"
            f"       if it was built from a 3' PAM file, name it so: mv '{lib}/NNNN_2_hg38' "
            f"'{lib}/<guide length>-NNNN_2_hg38', which records that\n", "NGCG")

    def test_a_plain_miss_says_what_it_said(self):
        # (but the index to download, of N = 3 for up to 2 bulges: it named the one of the
        # search's N, bulges + 1, and no index of N = 2 is published)
        self.install("NGG_3", variant=False)
        self.assertNoMatchSaysWhy(self.run_pipeline(prune=1, variant=False), "", "NRG")
        self.assertNoMatchSaysWhy(self.run_pipeline(prune=1, variant=False, bulges=3), "",
                                  "NRG", n=4)


class AnUnfinishedBuildIsNotSearched(_ShippedScript):
    """build-index-only marks each part of a dataset that it writes in place -- the
    dict-less tiers, the phased indel genotype store, the enriched genome as it goes
    into place, the SNP codes of the fake-indel genome -- with a hidden .building
    until that part is whole (crisprme.py _mark_unfinished). A run stopped there
    leaves a part that looks whole. complete-search refuses it in its pre-flight, but
    a search launched from the web runs this script straight, and used it as if it
    were whole, with exit status 0: so the script refuses it too, before it uses or
    builds anything. The fake-indel genome is read only to build an INDELS index the
    search does not find, and is refused there."""

    MARKED = ("Dictionaries/registry_{}", "Dictionaries/indel_genotypes_{}",
              "Genomes/hg38+{}")

    def run_marked(self, marked, env=None):
        """A run with the part `marked` marked unfinished by a run in another
        container; returns (exit status, the run, the marked directory)."""
        where = self.home / marked.format(self.VCF)
        made = not where.exists()
        (where / ".building").mkdir(parents=True)
        (where / ".building" / "owner").write_text(
            "host=otherhost\npidns=pid:[1]\npid=4242\nstarted=2026-10-31T00:00:00Z\n")
        launched = self.launch(env=env)
        run = self.finish(launched)
        shutil.rmtree(where if made else where / ".building")
        return launched.proc.returncode, run, where

    def test_each_part_marked_unfinished_is_refused_before_anything_is_used(self):
        self.install("NRG_3")
        # (and how to search without the tiers, which only a run given --samplesID
        # writes, or without the indel genotype store)
        then = {self.MARKED[0]: " with {} and its --samplesID (or build it again in the web "
                                "Settings page, which passes the dataset's sample list), or "
                                "remove it and {} to search without them,",
                self.MARKED[1]: ", or remove it to search without it,"}
        for marked in self.MARKED:
            with self.subTest(marked=marked):
                rc, run, where = self.run_marked(marked)
                self.assertEqual(rc, 1, run.stdout)
                do = then.get(marked, ",").format(
                    self.VCF, self.home / "Dictionaries" / f"genotypes_{self.VCF}")
                self.assertEqual(run.stderr, f"ERROR: {where} is unfinished: a build-index-only "
                                             "run is writing it, or stopped before it was done. "
                                             "Let that run finish, or run build-index-only "
                                             f"again to its end{do} then run this search again\n")
                self.assertEqual(run.calls, [])
                self.assertFalse(run.precomputed, run.stdout)
                self.assertNotIn("Variants already added", run.stdout)

    def test_the_fake_indel_genome_is_refused_only_where_an_indels_index_is_built_from_it(self):
        """With its INDELS index in place the search goes on, as while build-index-only
        builds a second PAM or bulge budget on the dataset (it used to be refused here
        too, though it never reads the fake-indel genome); without it, it is refused at
        STEP 2.3, before anything is built from it."""
        self.install("NRG_3")
        marked = "Genomes/hg38+{}_INDELS"
        rc, run, where = self.run_marked(marked)
        self.assertTrue(run.precomputed, run.stdout)
        self.assertEqual(run.indels, f"NRG_3_hg38+{self.VCF}_INDELS")
        self.assertEqual(run.searched(variant=True), (f"NRG_3_hg38+{self.VCF}", False))
        shutil.rmtree(self.library / f"NRG_3_hg38+{self.VCF}_INDELS")
        rc, run, where = self.run_marked(marked)
        self.assertEqual(rc, 1, run.stdout)
        self.assertEqual(run.stderr, f"ERROR: {where} is unfinished: a build-index-only run is "
                                     "writing it, or stopped before it was done. Let that run "
                                     "finish, or run build-index-only again to its end, then run "
                                     "this search again\n")
        self.assertEqual(run.calls, [])

    def test_an_indels_index_built_while_its_fake_indel_genome_changed_is_not_kept(self):
        """build-index-only marks the fake-indel genome only just before the first FASTA it
        changes as it puts the SNP codes on it in place: a search past its check of the mark
        before then built its INDELS index from the FASTAs not coded yet, blind to their
        SNPs, and every later run kept it (R2-F1). An INDELS index during whose build a
        FASTA of the genome changed, or after which the mark is there, is not kept."""
        self.install("NRG_3")
        indels = self.library / f"NRG_3_hg38+{self.VCF}_INDELS"
        shutil.rmtree(indels)
        genome = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS"
        fake = genome / f"fake_{self.VCF}_chr1" / "fakechr1.fa"
        fake.parent.mkdir(parents=True)
        fake.write_text(">fakechr1\nACGTACGTACGTACGTACGTACGTACGT\n")
        for overlay in ("coded a FASTA", "is at it"):
            with self.subTest(overlay=overlay):
                built = self.root / f"built {overlay}"
                launched = self.launch(env={"STUB_BUILD_UNTIL": str(built),
                                            "STUB_BUILD_WRITES": "1"})
                _wait_for(lambda: any(c[0] == "index-genome" for c in launched.calls()), 60,
                          "the INDELS build")
                if overlay == "coded a FASTA":  # as overlay_fake_chromosome writes one
                    staged = fake.parent / ".overlay-x"
                    staged.mkdir()
                    (staged / "fakechr1.fa.tmp").write_text(
                        ">fakechr1\nACGTACGTACGTRCGTACGTACGTACGT\n")
                    os.replace(staged / "fakechr1.fa.tmp", fake)
                    staged.rmdir()
                else:
                    (genome / ".building").mkdir()
                built.touch()
                run = self.finish(launched)
                self.assertEqual(launched.proc.returncode, 1)
                self.assertEqual(run.stderr, (
                    f"ERROR: {genome} changed while this search built NRG_3_hg38+{self.VCF}"
                    "_INDELS from it: a build-index-only run puts the SNP codes on it, and the "
                    "index would be blind to those not on yet, so it is not kept. Let that run "
                    "finish, then run this search again\n"))
                self.assertFalse(indels.exists())
                self.assertFalse(Path(f"{indels}.lock").exists())
                self.assertIsNone(run.searched(variant=True))
                shutil.rmtree(genome / ".building", ignore_errors=True)
        # with nothing changed meanwhile, it is kept and searched, as before
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"})
        self.assertTrue(indels.is_dir())
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed")

    def test_the_indel_genotype_store_is_not_refused_when_it_is_not_read(self):
        # CRISPRME_INDEL_SNP=0: no SNP+indel co-occurrence, so the store is not read
        self.install("NRG_3")
        rc, run, where = self.run_marked(self.MARKED[1], env={"CRISPRME_INDEL_SNP": "0"})
        self.assertTrue(run.precomputed, run.stdout)
        self.assertEqual(run.searched(variant=True), (f"NRG_3_hg38+{self.VCF}", False))

    def test_a_finished_build_and_a_reference_search_are_as_before(self):
        self.install("NRG_3")
        run = self.run_pipeline()
        self.assertTrue(run.precomputed, run.stdout)
        self.assertEqual(run.searched(variant=True), (f"NRG_3_hg38+{self.VCF}", False))
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed")
        # a reference search reads none of the dataset's parts
        (self.home / "Dictionaries" / f"registry_{self.VCF}" / ".building").mkdir()
        run = self.run_pipeline(variant=False)
        self.assertEqual(run.searched(variant=False), ("NRG_3_hg38", False))


def _lock_helpers():
    """The shipped build-lock block, verbatim: from its header comment to its last trap,
    the TERM one.

    Read when a test needs it, not at import, so a script without the block fails
    those tests instead of the whole module."""
    text = _SCRIPT.read_text()
    start = text.index("# Build locks.")
    last = re.compile(r"^trap '[^\n]*' TERM$", re.M).search(text, start)
    return text[start:last.end()] + "\n"


def _this_host():
    """What the script records for this host: `uname -n` and its PID namespace."""
    try:
        pidns = os.readlink("/proc/self/ns/pid")
    except OSError:  # no /proc (macOS): the script records it empty too
        pidns = ""
    return os.uname().nodename, pidns


def owner_record(pid, host=None, pidns=None, started="2026-10-30T00:00:00Z", pstart=None):
    """An owner record as the script writes one into a lock; for this host by default.
    `pstart`: the process's start time, as the script records its own (none: an owner
    this run can tell by its PID alone)."""
    this_host, this_pidns = _this_host()
    return (f"host={this_host if host is None else host}\n"
            f"pidns={this_pidns if pidns is None else pidns}\n"
            f"pid={pid}\nstarted={started}\n"
            + ("" if pstart is None else f"pstart={pstart}\n"))


def _pstart(pid):
    """`pid`'s start time: field 22 of /proc/<pid>/stat, after its name in parentheses."""
    stat = Path(f"/proc/{pid}/stat").read_text()
    return stat[stat.rindex(")") + 2:].split()[19]


def dead_pid():
    """The PID of a process that has exited and been reaped."""
    proc = subprocess.Popen(["true"])
    proc.wait()
    return proc.pid


def _wait_for(predicate, timeout, what):
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError(f"waited {timeout} s for {what}")
        time.sleep(0.05)


def _q(path):
    return shlex.quote(str(path))


def _age(path, minutes):
    """Set `path`'s mtime `minutes` back: a record its holder last refreshed then."""
    then = time.time() - minutes * 60
    os.utime(path, (then, then))


# The process table is read from /proc, as the script reads it: the image has no ps
# (no procps), and every test that asked ps failed there, each leaving its lock behind
# for the next. ps is asked only where there is no /proc (macOS).
_PROC = os.path.exists("/proc/self/stat")


def _proc_table():
    """(pid, ppid, pgid, state, command line) of every process, from /proc; the fields
    follow the last ") " of /proc/<pid>/stat, as the name before it may hold either."""
    table = []
    for pid in filter(str.isdigit, os.listdir("/proc")):
        try:
            stat = Path(f"/proc/{pid}/stat").read_text(errors="replace")
            args = Path(f"/proc/{pid}/cmdline").read_bytes()
        except OSError:  # gone meanwhile
            continue
        state, ppid, group = stat[stat.rindex(")") + 2:].split()[:3]
        table.append((int(pid), int(ppid), int(group), state,
                      args.rstrip(b"\0").replace(b"\0", b" ").decode(errors="replace")))
    return table


def _ps_table():
    """The same table as ps gives it."""
    out = subprocess.run(["ps", "-ww", "-e", "-o", "pid=,ppid=,pgid=,stat=,args="],
                         capture_output=True, text=True).stdout
    table = []
    for line in out.splitlines():
        pid, ppid, group, stat, *args = line.split(None, 4)
        table.append((int(pid), int(ppid), int(group), stat[:1], args[0] if args else ""))
    return table


def _processes(pgid):
    """(pid, ppid, command line) of every process still running -- zombies left out --
    in process group `pgid`."""
    return [(pid, ppid, args) for pid, ppid, group, state, args in
            (_proc_table() if _PROC else _ps_table()) if group == pgid and state != "Z"]


def _state(pid):
    """`pid`'s state ("Z": exited but not reaped), or "" once it is gone."""
    if not _PROC:
        out = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], capture_output=True,
                             text=True).stdout
        return out.strip()[:1]
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(errors="replace")
    except OSError:
        return ""
    return stat[stat.rindex(")") + 2:][:1]


def _subshells(pid):
    """The background subshells of bash process `pid`: its children that run its own
    command line -- a lock's refresh loop is one. Seen in two looks 0.2 s apart, as a
    child caught between its fork and its exec shows its parent's command line too."""
    def look():
        procs = _processes(os.getpgid(pid))
        me = [args for p, _, args in procs if p == pid]
        return {p for p, ppid, args in procs if ppid == pid and me and args == me[0]}
    first = look()
    time.sleep(0.2)
    return sorted(first & look())


def _wait_until_gone(pgid, timeout=10):
    """Wait for every process of group `pgid` to end; what is still running then."""
    deadline = time.monotonic() + timeout
    while _processes(pgid) and time.monotonic() < deadline:
        time.sleep(0.05)
    return _processes(pgid)


@unittest.skipUnless(_PROC, "the helpers read /proc where there is one")
class TheLockTestsNeedNoPs(unittest.TestCase):
    """In the shipped image, which has no ps, the lock tests below errored: 43 of them on
    ps itself and the rest on the locks those left behind. The helpers that list and
    check processes read /proc now, and are run here with no ps on PATH, as there."""

    def sleeper(self):
        """A `sleep 60` that leads a process group of its own, left unreaped once dead."""
        proc = subprocess.Popen([shutil.which("sleep"), "60"], start_new_session=True)
        self.addCleanup(_reap, proc)
        return proc

    def test_a_process_is_followed_until_it_is_gone_with_no_ps_on_path(self):
        proc = self.sleeper()
        row = (proc.pid, os.getpid(), f"{shutil.which('sleep')} 60")
        with tempfile.TemporaryDirectory() as empty, \
                mock.patch.dict(os.environ, {"PATH": empty}):
            _wait_for(lambda: _state(proc.pid) == "S", 10, "it to sleep")
            self.assertEqual(_processes(proc.pid), [row])
            os.kill(proc.pid, signal.SIGSTOP)
            _wait_for(lambda: _state(proc.pid) == "T", 10, "it to stop")
            self.assertEqual(_processes(proc.pid), [row])
            os.kill(proc.pid, signal.SIGKILL)  # and no wait(): left unreaped
            _wait_for(lambda: _state(proc.pid) == "Z", 10, "it to die")
            self.assertEqual(_processes(proc.pid), [])  # a zombie is left out
            proc.wait()
            self.assertEqual(_state(proc.pid), "")

    @unittest.skipUnless(shutil.which("ps"), "no ps here to compare with, as in the image")
    def test_where_there_is_ps_it_gives_the_same_answers(self):
        proc = self.sleeper()

        def both():
            """What the helpers say through /proc, and what they say through ps."""
            with mock.patch.object(sys.modules[__name__], "_PROC", False):
                via_ps = (_processes(proc.pid), _state(proc.pid))
            return (_processes(proc.pid), _state(proc.pid)), via_ps

        for signum, state in ((None, "S"), (signal.SIGSTOP, "T"), (signal.SIGKILL, "Z")):
            if signum is not None:
                os.kill(proc.pid, signum)
            _wait_for(lambda: _state(proc.pid) == state, 10, f"state {state}")
            via_proc, via_ps = both()
            self.assertEqual(via_proc, via_ps)
        proc.wait()
        self.assertEqual(both(), (([], ""), ([], "")))


class _Locks(unittest.TestCase):
    """Runs the shipped lock helpers in bash against real lock directories, polling
    every 0.2 s and repeating the waiting notice every 1 s (the script: 5 s, 5 min)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.lock = self.dir / "NRG_2_hg38.lock"
        self.out = self.dir / "NRG_2_hg38"
        self.started = 0

    def start(self, body, as_nobody=False):
        """Start `body` in bash after the shipped helpers; its output goes to files.
        `as_nobody` runs it as uid/gid 65534 (only root can), with the lock directory
        open to that user."""
        script = f"_lock_poll_s=0.2\n_lock_notice_s=1\n{_lock_helpers()}{body}\n"
        self.started += 1
        out, err = self.dir / f"stdout{self.started}", self.dir / f"stderr{self.started}"
        drop = None
        if as_nobody:
            self.dir.chmod(0o755)
            drop = lambda: (os.setgid(65534), os.setuid(65534))  # noqa: E731
        with open(out, "w") as so, open(err, "w") as se:
            proc = subprocess.Popen(["bash", "-c", script], stdout=so, stderr=se,
                                    text=True, start_new_session=True, preexec_fn=drop)
        self.addCleanup(_reap, proc)
        proc.out_path, proc.err_path = out, err
        return proc

    def result(self, proc, stop_after=None):
        """(returncode, stdout, stderr) once `proc` ends. With `stop_after`, a run still
        waiting then is stopped with SIGKILL to its process group: what INT and TERM do
        is tested on a run that holds a lock, not on one waiting for it."""
        try:
            proc.wait(timeout=stop_after or 60)
        except subprocess.TimeoutExpired:
            if stop_after is None:
                raise
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=30)
        return proc.returncode, proc.out_path.read_text(), proc.err_path.read_text()

    def take_body(self, then=""):
        """_lock_take on the lock, then `then`; exits with what _lock_take returned."""
        return (f'_lock_take {_q(self.lock)} "indexing {self.out}" {_q(self.out)}\n'
                f"rc=$?\n{then}\nexit $rc")

    def take(self, then="", stop_after=None):
        return self.result(self.start(self.take_body(then)), stop_after)

    def hold(self, record):
        """Leave the lock as another run holds it, `record` its owner file -- None for
        no owner file, as an older CRISPRme left every lock."""
        self.lock.mkdir()
        if record is not None:
            (self.lock / "owner").write_text(record)

    def sleeper(self, leader=False):
        """A live process, to stand for a run that holds a lock; with `leader`, one that
        leads a process group of its own."""
        proc = subprocess.Popen(["sleep", "60"], start_new_session=leader)
        self.addCleanup(lambda: (proc.kill(), proc.wait()))
        return proc


class ADeadRunsLockDoesNotWedgeTheNext(_Locks):
    """A run that died holding a build lock -- SIGKILL, an OOM kill, a container stop,
    a reboot -- used to leave the lock for good, and every later run that needed the
    same index waited on it forever: no error, no message, a search slot held. The
    lock now names its owner. A dead owner on this host is taken over, and so is a
    record that has gone 15 minutes without a refresh, here or elsewhere; one this run
    cannot check is waited on, never broken, and the wait is explained."""

    def notices(self, out):
        return [line for line in out.splitlines() if line.startswith("Still waiting for ")]

    def test_the_lock_helpers_were_actually_found(self):
        block = _lock_helpers()
        for part in ("_lock_take() {", "_lock_stale() {", "_lock_steal() {",
                     "_lock_release() {", "trap '_lock_abandon_all' EXIT"):
            self.assertIn(part, block)

    def test_a_free_lock_is_taken_with_this_runs_owner_record(self):
        rc, out, err = self.take(then=f"cat {_q(self.lock)}/owner; echo self=$$")
        self.assertEqual(rc, 0, out + err)
        host, pidns = _this_host()
        me = re.search(r"^self=(\d+)$", out, re.M).group(1)
        self.assertRegex(out, rf"host={re.escape(host)}\npidns={re.escape(pidns)}\n"
                              rf"pid={me}\nstarted=\S+\n")
        self.assertEqual(err, "")

    def test_a_dead_owner_on_this_host_is_taken_over(self):
        self.hold(owner_record(dead_pid()))
        rc, out, err = self.take(then=f"cat {_q(self.lock)}/owner; echo self=$$")
        self.assertEqual(rc, 2, out + err)  # 2: taken over, so the caller redoes the build
        self.assertIn(f"Stale lock {self.lock}: the run that held it", out)
        me = re.search(r"^self=(\d+)$", out, re.M).group(1)
        self.assertIn(f"\npid={me}\n", out)  # the record names this run now
        self.assertNotIn("waiting", out)
        self.assertEqual(err, "")

    def test_a_live_owner_is_waited_for_and_taken_over_once_it_dies(self):
        owner = self.sleeper()
        record = owner_record(owner.pid)
        self.hold(record)
        waiter = self.start(self.take_body())
        time.sleep(1.5)  # several polls
        self.assertIsNone(waiter.poll(), "a live owner's lock was not waited for")
        self.assertEqual((self.lock / "owner").read_text(), record)
        owner.kill()
        owner.wait()
        rc, out, err = self.result(waiter)
        self.assertEqual(rc, 2, out + err)
        self.assertIn(f"Another run is indexing {self.out}", out)
        self.assertIn(f"pid {owner.pid} on {_this_host()[0]}", self.notices(out)[0])
        self.assertIn("which is still running", self.notices(out)[0])
        # a run known to be alive: wait for it or stop it -- never "clear its lock",
        # which would leave its build writing into a directory another run rebuilds
        self.assertTrue(self.notices(out)[0].endswith(
            "which is still running. Wait for that run to finish, or stop it."), out)
        self.assertNotIn("rm -rf", out)
        self.assertIn(f"Stale lock {self.lock}", out)
        self.assertEqual(err, "")

    def assertWaitsAndExplains(self, record, why, aged=0, live=False):
        """The lock is waited on and never broken, and the wait is explained -- at the
        notice interval, not at every poll: with how to clear it when this run cannot
        check the owner, with "wait for it" when the owner is known to be alive (`live`).
        `aged`: minutes since the record (or, with none, the lock) was last written."""
        self.hold(record)
        _age(self.lock / "owner" if record is not None else self.lock, aged)
        rc, out, err = self.take(stop_after=3.2)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertNotIn("Stale lock", out)
        notices = self.notices(out)
        self.assertGreaterEqual(len(notices), 2, out)  # repeated ...
        self.assertLessEqual(len(notices), 6, out)  # ... but not at each of ~16 polls
        for notice in notices:
            self.assertIn(f"Still waiting for {self.lock}: ", notice)
            self.assertIn(why, notice)
            if live:
                self.assertTrue(notice.endswith(
                    ". Wait for that run to finish, or stop it."), notice)
                self.assertNotIn("rm -rf", notice)
            else:
                self.assertTrue(notice.endswith(
                    f"clear the lock and what it left with: rm -rf '{self.lock}' "
                    f"'{self.out}'"), notice)
        # untouched, and not released by the trap: the lock was never this run's
        self.assertEqual(sorted(p.name for p in self.lock.iterdir()),
                         [] if record is None else ["owner"])
        if record is not None:
            self.assertEqual((self.lock / "owner").read_text(), record)
        self.assertEqual(err, "")

    def test_a_lock_with_no_owner_record_is_never_broken(self):
        # as an older CRISPRme leaves every lock: nothing says whose it is, nor that
        # its run is gone -- however long ago it was taken
        self.assertWaitsAndExplains(None, "it has no owner record", aged=60)

    # what a waiting run says of a record that another host or container refreshed lately
    ELSEWHERE = ("on another host or in another container, which refreshed it in the "
                 "last 15 minutes; this run takes it over once it goes 15 minutes "
                 "without a refresh")

    def test_a_lock_held_on_another_host_is_never_broken(self):
        # a PID from another host -- or another Docker container -- means nothing here,
        # but its run refreshes the record while it lives
        self.assertWaitsAndExplains(owner_record(dead_pid(), host="another-host"),
                                    self.ELSEWHERE, live=True)

    def test_a_lock_held_in_another_pid_namespace_is_never_broken(self):
        # a container that shares the hostname (docker run --network host) still has
        # PIDs of its own
        self.assertWaitsAndExplains(owner_record(dead_pid(), pidns="pid:[1]"),
                                    self.ELSEWHERE, live=True)

    def test_a_lock_held_elsewhere_and_refreshed_lately_is_never_broken(self):
        # refreshed 14 minutes ago: its run may be alive, just slow to refresh
        for record in (owner_record(dead_pid(), host="another-host"),
                       owner_record(dead_pid(), pidns="pid:[1]")):
            with self.subTest(record=record):
                self.assertWaitsAndExplains(record, self.ELSEWHERE, aged=14, live=True)
                shutil.rmtree(self.lock)

    def test_a_lock_held_elsewhere_and_not_refreshed_for_15_minutes_is_taken_over(self):
        # The Docker case: every `docker run` has a hostname and PID namespace of its
        # own, so the next run cannot check the PID; but a killed container's run never
        # refreshes its record again. Taken over as a dead run's.
        for host, pidns in (("another-host", None), (None, "pid:[1]")):
            with self.subTest(host=host, pidns=pidns):
                self.hold(owner_record(4242, host=host, pidns=pidns))
                _age(self.lock / "owner", 16)
                rc, out, err = self.take(then=f"cat {_q(self.lock)}/owner; echo self=$$\n"
                                              f"_lock_release {_q(self.lock)}",
                                         stop_after=10)
                self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
                self.assertIn(f"Stale lock {self.lock}: the run that held it (pid 4242 on "
                              f"{host or _this_host()[0]}, started 2026-10-30T00:00:00Z) "
                              "has not refreshed its lock for 15 minutes. Taking it over",
                              out)
                me = re.search(r"^self=(\d+)$", out, re.M).group(1)
                self.assertIn(f"\npid={me}\n", out)  # the record names this run now
                self.assertNotIn("waiting", out)
                self.assertFalse(self.lock.exists())  # and released, its build done
                self.assertEqual(err, "")

    def test_a_live_owner_this_run_may_not_signal_is_never_broken(self):
        # kill -0 fails on a live process of another user as well as on a dead one;
        # only "No such process" means gone. PID 1 is another user's to a run that is
        # not root -- and root runs the waiter as nobody here to make it so. With its
        # record refreshed lately it is a live run: waited for, never offered to clear.
        self.hold(owner_record(1))
        rc, out, err = self.result(self.start(self.take_body(), as_nobody=os.geteuid() == 0),
                                   stop_after=1.5)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertNotIn("Stale lock", out)
        notice = self.notices(out)[0]
        self.assertIn("a process this run may not signal (another user's), which refreshed "
                      "it in the last 15 minutes", notice)
        self.assertTrue(notice.endswith(". Wait for that run to finish, or stop it."), notice)
        self.assertNotIn("rm -rf", out)
        self.assertEqual((self.lock / "owner").read_text(), owner_record(1))
        self.assertEqual(err, "")

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a stopped run is told from a live one through /proc")
    def test_a_stopped_owner_on_this_host_is_waited_for_however_old_its_record(self):
        # Stopped (Ctrl-Z, kill -STOP, a SLURM suspend), a run's refresh loop is stopped
        # with it, so its record ages. A run on this host took the lock over after 15
        # minutes and rebuilt under the live build; resumed, that build wrote on, and its
        # run searched the index being rebuilt -- or removed it.
        # How to go on: kill -CONT of the run's PID resumes bash alone, its build step and
        # refresh loop left stopped, so a run that leads its own process group is resumed
        # or ended as a group. One that does not -- in the group of crisprme.py's job, or
        # of the web server -- is offered no signal to that group, nor one to bash alone.
        for leader, tip in (
                (True, "Resume it with fg or kill -CONT -- -{pid}, or end it with kill -KILL "
                       "-- -{pid}."),
                (False, "It must be resumed together with its stopped child processes, as fg "
                        "does: kill -CONT {pid} resumes the script alone.")):
            with self.subTest(leads_its_group=leader):
                shutil.rmtree(self.lock, ignore_errors=True)
                owner = self.sleeper(leader)
                os.kill(owner.pid, signal.SIGSTOP)
                _wait_for(lambda: _state(owner.pid) == "T", 10, "the owner to stop")
                record = owner_record(owner.pid)
                self.hold(record)
                _age(self.lock / "owner", 16)
                rc, out, err = self.take(stop_after=3.2)
                self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
                self.assertNotIn("Stale lock", out)
                notices = self.notices(out)
                self.assertGreaterEqual(len(notices), 2, out)
                for notice in notices:  # stopped, not finishing: how to go on
                    self.assertTrue(notice.endswith(
                        f"it is held by pid {owner.pid} on {_this_host()[0]} (a run started "
                        f"2026-10-30T00:00:00Z), which is stopped. {tip.format(pid=owner.pid)}"),
                        notice)
                self.assertNotIn("rm -rf", out)
                self.assertEqual((self.lock / "owner").read_text(), record)
                self.assertEqual(err, "")

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a dead run that was never reaped is told from a live one through /proc")
    def test_a_zombie_owner_on_this_host_is_taken_over_at_once(self):
        # Dead but never reaped, a run still answers kill -0: a zombie. Its record, just
        # refreshed, is a dead run's all the same.
        zombie = subprocess.Popen(["true"])
        self.addCleanup(zombie.wait)
        _wait_for(lambda: _state(zombie.pid) == "Z", 10, "the owner to die")
        self.hold(owner_record(zombie.pid))
        rc, out, err = self.take(stop_after=10)
        self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {zombie.pid} on "
                      f"{_this_host()[0]}, started 2026-10-30T00:00:00Z) is gone. Taking it "
                      "over", out)
        self.assertNotIn("waiting", out)
        self.assertEqual(err, "")

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a process's state is read through /proc")
    def test_an_owner_whose_name_holds_a_parenthesis_is_read_rightly(self):
        # /proc/<pid>/stat is "<pid> (<name>) <state> ...", and the name may hold spaces and
        # parentheses: a live owner named "x) Z y" is no zombie.
        link = self.dir / "x) Z y"
        link.symlink_to(shutil.which("sleep"))
        owner = subprocess.Popen([str(link), "60"])
        self.addCleanup(lambda: (owner.kill(), owner.wait()))
        _wait_for(lambda: Path(f"/proc/{owner.pid}/stat").read_text().startswith(
            f"{owner.pid} (x) Z y) "), 10, "the owner's name")
        self.hold(owner_record(owner.pid))
        rc, out, err = self.take(stop_after=1.5)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertIn(", which is still running. Wait for that run to finish", out)
        self.assertEqual(err, "")

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a process's start time is read through /proc")
    def test_the_very_run_is_waited_for_however_old_its_record(self):
        # Its record names its start time: this host's process of that PID with that start
        # time is that run. Resumed after 15 minutes stopped (Ctrl-Z, then fg), it was taken
        # over at once, its record not refreshed yet -- not until its refresh loop's next
        # beat, a minute later at most -- and both runs wrote one index.
        owner = self.sleeper()
        record = owner_record(owner.pid, pstart=_pstart(owner.pid))
        self.hold(record)
        _age(self.lock / "owner", 16)
        rc, out, err = self.take(stop_after=1.5)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertNotIn("Stale lock", out)
        self.assertTrue(self.notices(out)[0].endswith(
            "which is still running. Wait for that run to finish, or stop it."), out)
        self.assertNotIn("rm -rf", out)
        self.assertEqual((self.lock / "owner").read_text(), record)
        self.assertEqual(err, "")

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a process's start time is read through /proc")
    def test_a_pid_given_since_to_another_process_is_taken_over_at_once(self):
        # The run that held the lock is gone and its PID reused: the process there now
        # started at another time. Its record, refreshed a moment ago, is a dead run's.
        owner = self.sleeper()
        self.hold(owner_record(owner.pid, pstart=int(_pstart(owner.pid)) - 1))
        rc, out, err = self.take(stop_after=10)
        self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {owner.pid} on "
                      f"{_this_host()[0]}, started 2026-10-30T00:00:00Z) is gone. Taking it "
                      "over", out)
        self.assertEqual(err, "")

    def test_a_lock_with_this_runs_own_record_is_not_waited_on(self):
        # Left by a release whose rm -rf failed, say: its record is this very run's, which
        # must not wait for itself.
        rc, out, err = self.take(
            then=f'_lock_take {_q(self.lock)} "indexing {self.out}" {_q(self.out)}; rc=$?',
            stop_after=10)
        self.assertEqual(rc, 2, out + err)
        self.assertIn(f"Stale lock {self.lock}: the run that held it", out)
        self.assertNotIn("waiting", out)
        self.assertEqual(err, "")

    def test_a_record_not_refreshed_for_15_minutes_is_abandoned_whatever_its_pid_says(self):
        # A live holder refreshes its record every minute, on this host as elsewhere, so
        # one 16 minutes old is a dead run's even though its PID answers kill -0 -- a
        # zombie's does, and so does a reused PID: a process of this user's, or of
        # another's (EPERM, as nobody under root).
        live = self.sleeper()
        for pid, as_nobody in ((live.pid, False), (1, os.geteuid() == 0)):
            with self.subTest(pid=pid):
                shutil.rmtree(self.lock, ignore_errors=True)
                self.hold(owner_record(pid))
                _age(self.lock / "owner", 16)
                if as_nobody:
                    self.lock.chmod(0o777)  # where the takeover writes its own record
                rc, out, err = self.result(self.start(self.take_body(), as_nobody=as_nobody),
                                           stop_after=10)
                self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
                self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {pid} on "
                              f"{_this_host()[0]}, started 2026-10-30T00:00:00Z) has not "
                              "refreshed its lock for 15 minutes. Taking it over", out)
                self.assertNotIn("waiting", out)
                self.assertEqual(err, "")

    def test_a_takeover_already_under_way_is_left_to_it(self):
        # two runs that find the same dead owner must not both take its lock over
        record = owner_record(dead_pid())
        self.hold(record)
        (self.lock / ".takeover").mkdir()  # the other run, between its check and its write
        rc, out, err = self.take(stop_after=1.5)
        self.assertNotIn("Stale lock", out)
        self.assertEqual((self.lock / "owner").read_text(), record)
        self.assertEqual(err, "")

    def test_a_takeover_more_than_a_minute_old_is_abandoned(self):
        # A takeover takes milliseconds. A run that ended inside one (SIGKILL, say) left
        # ".takeover" in the lock, and every later takeover failed on it: no run ever took
        # that lock over, and each waited for good.
        self.hold(owner_record(dead_pid()))
        (self.lock / ".takeover").mkdir()
        _age(self.lock / ".takeover", 60)
        rc, out, err = self.take(then=f"ls -A {_q(self.lock)}", stop_after=10)
        self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
        self.assertIn(f"Stale lock {self.lock}: the run that held it", out)
        self.assertEqual(out.splitlines()[-1], "owner")  # this run's record, nothing else
        self.assertEqual(err, "")

    def test_a_run_that_ends_inside_its_takeover_leaves_no_takeover_behind(self):
        # Once its mkdir has made ".takeover", a run counts it as its own: ended inside its
        # takeover from then on -- Ctrl-C as it signs the lock, say -- it removes it, and
        # the next run takes the lock over at once. Stood for by the signing's mv, hung.
        record = owner_record(dead_pid())
        self.hold(record)
        made = self.dir / "made"
        proc = self.start(f'mv() {{ : >{_q(made)}; sleep 30; command mv "$@"; }}\n'
                          + self.take_body())
        _wait_for(made.exists, 10, "the takeover")
        self.assertTrue((self.lock / ".takeover").is_dir())
        os.killpg(proc.pid, signal.SIGINT)
        rc, out, err = self.result(proc)
        self.assertEqual(rc, -signal.SIGINT, out + err)
        self.assertNotIn(".takeover", os.listdir(self.lock))
        self.assertEqual((self.lock / "owner").read_text(), record)
        self.assertEqual(_wait_until_gone(proc.pid), [])
        self.assertEqual(err, "")
        rc, out, err = self.take()  # and the next run takes it over at once
        self.assertEqual(rc, 2, out + err)

    def test_a_run_that_ends_as_its_takeover_begins_leaves_it_to_the_minute_rule(self):
        # Ctrl-C lands while bash waits for the takeover's mkdir, which a takeover spends
        # most of its time in, and bash runs the trap as soon as the mkdir returns, before
        # the run can count ".takeover" as its own: it is left. The next run waits, as for a
        # takeover under way, and takes the lock over once that ".takeover" is a minute old.
        # Stood for by a mkdir that hangs once it has made ".takeover".
        self.hold(owner_record(dead_pid()))
        made = self.dir / "made"
        proc = self.start(f'mkdir() {{ command mkdir "$@" || return; '
                          f'[ "${{1##*/}}" = .takeover ] || return 0; : >{_q(made)}; sleep 30; }}\n'
                          + self.take_body())
        _wait_for(made.exists, 10, "the takeover")
        os.killpg(proc.pid, signal.SIGINT)
        rc, out, err = self.result(proc)
        self.assertEqual(rc, -signal.SIGINT, out + err)
        self.assertEqual(sorted(os.listdir(self.lock)), [".takeover", "owner"])
        self.assertEqual(_wait_until_gone(proc.pid), [])
        self.assertEqual(err, "")
        rc, out, err = self.take(stop_after=1.5)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertNotIn("Stale lock", out)
        _age(self.lock / ".takeover", 2)
        rc, out, err = self.take(stop_after=10)
        self.assertEqual(rc, 2, out + err)
        self.assertIn(f"Stale lock {self.lock}: the run that held it", out)
        self.assertEqual(err, "")

    def test_a_failed_takeover_never_clears_another_runs_takeover(self):
        # Another run is inside its own takeover of the same dead run's lock, so this run's
        # mkdir of ".takeover" fails -- and Ctrl-C lands meanwhile, its trap run as soon as
        # the mkdir returns. This run removed the other's ".takeover", and a third run could
        # then take the lock over beside it: two builds of one index. Stood for by a mkdir
        # that sends the run INT when it fails on ".takeover".
        self.hold(owner_record(dead_pid()))
        (self.lock / ".takeover").mkdir()  # the other run's
        proc = self.start('mkdir() { command mkdir "$@"; local rc=$?\n'
                          '\t[ "${1##*/}" = .takeover ] && [ $rc != 0 ] && kill -INT $$\n'
                          '\treturn $rc; }\n' + self.take_body())
        rc, out, err = self.result(proc, stop_after=10)
        self.assertEqual(rc, -signal.SIGINT, out + err)
        self.assertEqual(sorted(os.listdir(self.lock)), [".takeover", "owner"])
        self.assertEqual(err, "")

    def test_a_lock_released_with_its_output_in_place_is_not_built_again(self):
        owner = self.sleeper()
        self.hold(owner_record(owner.pid))
        waiter = self.start(self.take_body())
        _wait_for(lambda: "Another run is" in waiter.out_path.read_text(), 10, "the wait")
        self.out.mkdir()  # the holder finished ...
        shutil.rmtree(self.lock)  # ... and released its lock
        rc, out, err = self.result(waiter)
        self.assertEqual(rc, 1, out + err)  # 1: the caller builds nothing
        self.assertEqual(err, "")

    def test_a_free_lock_whose_output_is_in_place_is_let_go(self):
        # The caller found the lock held, and its holder finished and released it before
        # this run's mkdir: the lock is free, but the output is final. Rebuilding it
        # would rewrite (or nest into) what the holder and its waiters are reading.
        self.out.mkdir()
        rc, out, err = self.take(then=f"[ -d {_q(self.lock)} ] && echo HELD || echo FREE")
        self.assertEqual(rc, 1, out + err)  # 1: the caller builds nothing
        self.assertIn("FREE", out)  # and the lock is let go at once
        self.assertNotIn("waiting", out)
        self.assertEqual(err, "")

    def test_a_dead_runs_lock_is_still_taken_over_with_its_output_in_place(self):
        # its output is that run's half-written work, not a finished one
        self.hold(owner_record(dead_pid()))
        self.out.mkdir()
        rc, out, err = self.take()
        self.assertEqual(rc, 2, out + err)
        self.assertEqual(err, "")

    def test_an_index_is_in_place_only_once_it_holds_a_bin(self):
        # -bin: index-genome creates the directory before it writes any .bin, and a
        # directory without one is no index (_resolve_index_dir skips it): build it
        take = (f'_lock_take -bin {_q(self.lock)} "indexing {self.out}" {_q(self.out)}\n'
                f"rc=$?\n_lock_release {_q(self.lock)}\nexit $rc")
        self.out.mkdir()
        rc, out, err = self.result(self.start(take))
        self.assertEqual(rc, 0, out + err)  # a free lock, the index to build
        self.assertFalse(self.lock.exists())  # (released, as once that build is done)
        (self.out / "chr1.bin").write_bytes(b"")
        rc, out, err = self.result(self.start(take))
        self.assertEqual(rc, 1, out + err)  # in place: nothing to build
        self.assertEqual(err, "")

    def test_a_lock_released_beside_a_bin_less_index_is_built_here(self):
        # cleared by hand, say, while its index directory had no .bin yet
        owner = self.sleeper()
        self.hold(owner_record(owner.pid))
        self.out.mkdir()
        waiter = self.start(f'_lock_take -bin {_q(self.lock)} "indexing {self.out}" '
                            f'{_q(self.out)}\nexit $?')
        _wait_for(lambda: "Another run is" in waiter.out_path.read_text(), 10, "the wait")
        shutil.rmtree(self.lock)
        rc, out, err = self.result(waiter)
        self.assertEqual(rc, 0, out + err)
        self.assertIn(f"{self.out} is not there", out)
        self.assertEqual(err, "")

    def test_a_lock_released_without_its_output_is_built_here(self):
        # The holder failed (it removes what it wrote), or the lock was cleared by hand.
        # A search on an index that is not there comes back EMPTY, not as an error, so
        # the waiter must not go on as if it had been built.
        owner = self.sleeper()
        self.hold(owner_record(owner.pid))
        waiter = self.start(self.take_body())
        _wait_for(lambda: "Another run is" in waiter.out_path.read_text(), 10, "the wait")
        shutil.rmtree(self.lock)
        rc, out, err = self.result(waiter)
        self.assertEqual(rc, 0, out + err)
        self.assertIn(f"{self.out} is not there", out)
        self.assertEqual(err, "")


class ARunThatEndsMidBuildLeavesItsLock(_Locks):
    """A run lets go of a lock once its build is done, or has failed and removed what it
    wrote (_lock_release). One that ends any other way -- an exit, Ctrl-C, SIGTERM, a
    hangup -- leaves that build half-written, so it leaves the lock and its record as
    they are, for the next run to take over and build again, and stops only the lock's
    refresh. Letting go would hand the half-written index on as a finished one. It
    never touches a lock it does not hold: another run's, or one taken over since."""

    def test_exit_leaves_this_runs_locks_and_keeps_its_status(self):
        theirs = self.dir / "theirs.lock"
        theirs.mkdir()
        (theirs / "owner").write_text(owner_record(1))
        lost = self.dir / "lost.lock"
        body = (f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                f"_lock_take {_q(lost)} 'indexing b' {_q(self.dir / 'b')} || exit 9\n"
                # another run has taken this one over since: it is no longer this run's
                f"printf 'host=elsewhere\\npid=1\\n' > {_q(lost)}/owner\n"
                f"cat {_q(self.lock)}/owner\nexit 3")
        proc = self.start(body)
        rc, out, err = self.result(proc)
        self.assertEqual(rc, 3, out + err)  # the run's own exit status survives the trap
        self.assertEqual((self.lock / "owner").read_text(), out,
                         "a lock whose build was not done was let go")
        self.assertEqual(_wait_until_gone(proc.pid), [])  # its refresh loops stopped
        self.assertEqual((theirs / "owner").read_text(), owner_record(1))
        self.assertEqual((lost / "owner").read_text(), "host=elsewhere\npid=1\n")
        self.assertEqual(err, "")
        # the next run takes it over, as a dead run's, and does the build
        rc, out, err = self.take()
        self.assertEqual(rc, 2, out + err)
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {proc.pid} ", out)

    def test_a_signal_leaves_this_runs_lock_and_still_ends_the_run(self):
        # Ctrl-C and a terminal's hangup reach the whole process group; `kill` reaches
        # bash alone
        for sig, to_group in ((signal.SIGINT, True), (signal.SIGTERM, False),
                              (signal.SIGHUP, True), (signal.SIGHUP, False)):
            with self.subTest(signal=sig.name, to_group=to_group):
                body = (f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                        "echo taken\nwhile :; do sleep 0.2; done")
                proc = self.start(body)
                _wait_for(lambda: "taken" in proc.out_path.read_text(), 10, "the lock")
                record = (self.lock / "owner").read_text()
                if to_group:
                    os.killpg(proc.pid, sig)
                else:
                    proc.send_signal(sig)
                rc, out, err = self.result(proc)
                self.assertEqual(rc, -sig, out + err)  # ended by the signal, as without a trap
                self.assertEqual((self.lock / "owner").read_text(), record,
                                 "an interrupted run let go of its lock")
                self.assertEqual(_wait_until_gone(proc.pid), [])
                # bash reports the sleep a hangup to the group killed ("Hangup"); nothing
                # else may reach stderr
                self.assertEqual([line for line in err.splitlines() if "Hangup" not in line],
                                 [])
                shutil.rmtree(self.lock)

    def test_a_lock_taken_over_meanwhile_is_neither_used_nor_let_go(self):
        # Another run took the lock over while this one was stopped, or could not refresh
        # it, and builds again under it. After a build step (_lock_lost) and when letting
        # go (_lock_release), this run says so on stderr, as a run that is to end in
        # failure, and leaves the lock as it is. A lock with no record is still this
        # run's: its own record failed to write.
        theirs = owner_record(4242, host="3f9c2a1b7d10")
        taken = (f"ERROR: the lock on {self.out} was taken over by another run (pid 4242 on "
                 "3f9c2a1b7d10, started 2026-10-30T00:00:00Z) while this run held it; that "
                 "run builds it again, and this run uses none of what it built. Re-run once "
                 "that run is done.")
        gone = (f"ERROR: the lock on {self.out} was taken over or cleared while this run "
                "held it; this run uses none of what it built. Re-run once no other run is "
                "building it.")
        for lock, change, says, released, left in (
                ("taken over", f"printf %s {_q(theirs)} >{_q(self.lock)}/owner", taken, 1,
                 ["owner"]),
                ("gone", f"rm -rf {_q(self.lock)}", gone, 1, None),
                ("with no record", f"rm {_q(self.lock)}/owner", None, 0, None)):
            with self.subTest(lock=lock):
                proc = self.start(f"_lock_beat_s=0.2\n"
                                  f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} "
                                  f"|| exit 9\n{change}\n"
                                  f"_lock_lost {_q(self.lock)}; echo lost=$?\n"
                                  f"_lock_release {_q(self.lock)}; echo released=$?")
                rc, out, err = self.result(proc)
                self.assertEqual(out, f"lost={0 if says else 1}\nreleased={released}\n", err)
                self.assertEqual(err, f"{says}\n" * 2 if says else "")
                if left is None:
                    self.assertFalse(self.lock.exists())
                else:
                    self.assertEqual(sorted(os.listdir(self.lock)), left)
                    self.assertEqual((self.lock / "owner").read_text(), theirs)
                    shutil.rmtree(self.lock)
                self.assertEqual(_wait_until_gone(proc.pid), [])  # its refresh stopped

    def test_a_waiting_run_expands_no_command_substitution(self):
        """bash 5.2 can fail to parse a trap that fires while it expands a $(...) --
        "trap: line 2: unexpected EOF while looking for matching `)'" -- and the signal
        is lost: a run waiting on a lock ignored SIGTERM, or abandoned the script's
        main loop, ran the code after it and exited 0. A waiting run spends its time in
        these helpers, so they run none: no $(...), no backticks. ($((...)) is
        arithmetic. _lock_me is set once, before any trap is.)"""
        code, lines = [], iter(_lock_helpers().splitlines())
        for line in lines:
            if line.startswith("_lock_me="):
                while line.endswith("\\"):
                    line = next(lines)
                continue
            code.append(re.sub(r"(^|\s)#.*", "", line))  # drop comments
        code = "\n".join(code)
        self.assertEqual(re.findall(r".*\$\((?!\().*", code), [])
        self.assertNotIn("`", code)

    def test_each_trap_is_one_function_call(self):
        traps = re.findall(r"^trap '([^']*)' (\w+)$", _lock_helpers(), re.M)
        self.assertEqual(sorted(sig for _, sig in traps), ["EXIT", "HUP", "INT", "TERM"])
        for body, sig in traps:
            self.assertRegex(body, r"^_lock_\w+( [A-Z]+)?$", sig)

    def test_a_child_of_the_run_never_acts_for_it(self):
        # bash 5.2 can run the run's traps in a child it has just forked, before the child
        # resets them: a refresh loop stopped at once ran the EXIT trap's abandon. A child
        # must neither stop the run's refresh loops (nor wait for what is not its own,
        # which bash reports on stderr) nor end the run: its INT, TERM or HUP, re-raised
        # on $$, reached the run itself.
        proc = self.start(f"_lock_beat_s=0.2\n"
                          f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                          '( _lock_abandon_all ) & wait "$!"\n'
                          '( _lock_on_signal TERM ) & wait "$!"\n'
                          'echo "child: $?"\nwhile :; do sleep 0.2; done')
        _wait_for(lambda: "child:" in proc.out_path.read_text(), 10, "the children")
        # ended by its own TERM, not the run: 143 -- or 0 on bash 5.2.21 (Ubuntu 24.04's,
        # as in CI), where an async subshell that is sent TERM exits 0 when an EXIT trap
        # was set before the TERM trap, as here
        self.assertRegex(proc.out_path.read_text(), r"child: (0|143)\n")
        owner = self.lock / "owner"
        record = owner.read_text()
        _age(owner, 60)
        _wait_for(lambda: time.time() - owner.stat().st_mtime < 30, 10, "a refresh")
        self.assertEqual(len(_subshells(proc.pid)), 1, _processes(proc.pid))  # the loop
        os.killpg(proc.pid, signal.SIGINT)
        rc, out, err = self.result(proc)
        self.assertEqual(rc, -signal.SIGINT, out + err)
        self.assertEqual(owner.read_text(), record)
        self.assertEqual(_wait_until_gone(proc.pid), [])
        self.assertEqual(err, "")


class ALockIsKeptFreshByTheRunThatHoldsIt(_Locks):
    """A run that cannot check the holder's PID -- in another container, on another host
    -- takes a lock whose record has gone 15 minutes without a refresh for a dead run's.
    So the holder refreshes the record every minute (here every 0.2 s) for as long as
    it holds the lock, from a background loop that stops when the lock is released or
    the holder is gone, and leaves nothing running behind it."""

    def hold_then(self, after_take):
        """A run that takes the lock with a 0.2 s refresh, says "taken", waits for the
        test's go and then runs `after_take`."""
        self.go = self.dir / "go"
        return self.start(f"_lock_beat_s=0.2\n"
                          f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                          f"echo taken\nwhile ! [ -e {_q(self.go)} ]; do sleep 0.05; done\n"
                          f"{after_take}")

    def assertRefreshed(self, proc):
        """The record's mtime is moved on, and nothing else about it changes."""
        _wait_for(lambda: "taken" in proc.out_path.read_text(), 10, "the lock")
        owner = self.lock / "owner"
        record = owner.read_text()
        _age(owner, 60)
        _wait_for(lambda: time.time() - owner.stat().st_mtime < 30, 10, "a refresh")
        self.assertEqual(owner.read_text(), record)
        self.assertEqual(len(_subshells(proc.pid)), 1, _processes(proc.pid))  # the loop

    def test_the_record_is_refreshed_until_the_lock_is_released(self):
        jobs = self.dir / "jobs"
        proc = self.hold_then(f"_lock_release {_q(self.lock)}\njobs -p >{_q(jobs)}\n"
                              "echo released\nsleep 1")
        self.assertRefreshed(proc)
        self.go.touch()
        _wait_for(lambda: "released" in proc.out_path.read_text(), 10, "the release")
        self.assertFalse(self.lock.exists())
        self.assertEqual(jobs.read_text(), "")  # the release stopped and reaped the loop
        self.assertEqual(_subshells(proc.pid), [])
        rc, out, err = self.result(proc)
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(_wait_until_gone(proc.pid), [])
        self.assertEqual(err, "")

    def test_the_refresh_stops_when_the_run_ends(self):
        # an exit, or INT, TERM or HUP, before the build is done: the lock is left for
        # the next run to take over, and no longer refreshed
        for how, status in (("exit 3", 3), ("kill -s TERM $$", -signal.SIGTERM),
                            ("kill -s INT $$", -signal.SIGINT),
                            ("kill -s HUP $$", -signal.SIGHUP)):
            with self.subTest(how=how):
                proc = self.hold_then(how)
                self.assertRefreshed(proc)
                self.go.touch()
                rc, out, err = self.result(proc)
                self.assertEqual(rc, status, out + err)
                self.assertEqual(_wait_until_gone(proc.pid), [])
                _age(self.lock / "owner", 60)
                time.sleep(1)  # five refresh intervals
                self.assertGreater(time.time() - (self.lock / "owner").stat().st_mtime, 3000)
                self.assertEqual(err, "")
                shutil.rmtree(self.lock)
                self.go.unlink()

    def test_the_refresh_stops_when_the_holder_is_killed_outright(self):
        # SIGKILL to bash alone runs no trap: the loop sees its run gone, and stops
        proc = self.hold_then("")
        self.assertRefreshed(proc)
        proc.kill()
        proc.wait()
        self.assertEqual(_wait_until_gone(proc.pid), [])
        stopped = (self.lock / "owner").stat().st_mtime
        time.sleep(1)  # five refresh intervals
        self.assertEqual((self.lock / "owner").stat().st_mtime, stopped)

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a dead run that was never reaped is told from a live one through /proc")
    def test_the_refresh_stops_when_the_holder_dies_unreaped(self):
        # Killed outright and never reaped -- its parent was killed first, say, and the
        # container's PID 1 does not reap -- the run is a zombie, which still answers
        # kill -0. Its lock must not stay fresh for good, or no run anywhere could ever
        # take it over.
        proc = self.hold_then("")
        self.assertRefreshed(proc)
        os.kill(proc.pid, signal.SIGKILL)  # bash alone, and no wait(): left unreaped
        _wait_for(lambda: _state(proc.pid) == "Z", 10, "the holder to die")
        self.assertEqual(_wait_until_gone(proc.pid), [])  # its refresh loop stopped
        _age(self.lock / "owner", 16)
        time.sleep(1)  # five refresh intervals
        self.assertGreater(time.time() - (self.lock / "owner").stat().st_mtime, 15 * 60)
        # and a run on this host, where the zombie's PID still answers, takes it over
        self.assertEqual(_state(proc.pid), "Z")
        rc, out, err = self.take(stop_after=10)
        self.assertEqual(rc, 2, out + err)
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {proc.pid} ", out)
        self.assertIn(") is gone. Taking it over", out)  # a zombie: dead, whatever its age
        self.assertEqual(err, "")

    @unittest.skipUnless(os.path.isdir("/proc/self/fd"),
                         "what a process has open is read through /proc")
    def test_the_loop_never_holds_the_runs_stderr(self):
        # A loop stopped as soon as it is forked -- its lock let go at once (a free lock
        # whose output is in place) or abandoned at once -- can run the run's EXIT trap
        # before it reaches a redirection of its own, and bash 5.2 then wrote "wait_for: No
        # record of process ..." on the stderr it still shared with the run: log_error.txt,
        # which fails the stage although nothing failed. So the loop has /dev/null for
        # its output from its fork on: looked at right after it.
        body = (f"_lock_beat_s=0.2\ni=0\nwhile [ $i -lt 100 ]; do i=$((i+1))\n"
                f"\tmkdir {_q(self.lock)}; _lock_sign {_q(self.lock)}; _lock_hold {_q(self.lock)}\n"
                '\t[ "/proc/$!/fd/2" -ef /dev/null ] || echo "loop $i holds the stderr"\n'
                f"\t_lock_release {_q(self.lock)}\ndone\necho done")
        rc, out, err = self.result(self.start(body))
        self.assertEqual((rc, out), (0, "done\n"), out + err)
        self.assertEqual(err, "")

    def test_a_loop_stopped_before_it_notes_its_sleep_still_reaps_it(self):
        # Stopped between forking its sleep and noting the sleep's PID -- a lock let go the
        # moment it is taken (a free lock whose output is in place), about 5% of the time
        # -- the loop left that sleep to PID 1, which in the web container never reaps it.
        # It now reaps what $! names then -- but only what it forked: before its first
        # sleep, $! names the run's own last background process (here a sleep 30). A DEBUG
        # trap, inherited by the loop under set -T, lands the TERM as the loop is about to
        # note its sleep, or to take its first look.
        for command, forked in (("_s=$!", True), ("_lock_run_alive", False)):
            with self.subTest(term_before=command):
                shutil.rmtree(self.lock, ignore_errors=True)
                pid_file = self.dir / f"pid_{forked}"
                trigger = command.replace("$", "\\$")
                proc = self.start(
                    f"_sp={_q(pid_file)}\nsleep 30 &\nset -T\n"
                    f"trap '[ -n \"$_hit\" ] || [ \"$BASH_COMMAND\" != \"{trigger}\" ] || "
                    '{ _hit=1; echo "$!" >"$_sp"; kill -s TERM "$BASHPID"; :; }\' DEBUG\n'
                    f"mkdir {_q(self.lock)}; _lock_sign {_q(self.lock)}\n"
                    f"_lock_hold {_q(self.lock)}\n"
                    'wait "${_lock_beats[0]}"; echo "loop: $?"')
                rc, out, err = self.result(proc)
                self.addCleanup(lambda g=proc.pid: [os.kill(p, signal.SIGKILL)
                                                    for p, _, _ in _processes(g)])
                self.assertTrue(pid_file.exists(), f"no TERM landed before {command}")
                pid = int(pid_file.read_text())
                self.assertEqual((rc, out), (0, "loop: 0\n"), out + err)
                self.assertEqual(_state(pid), "" if forked else "S",
                                 "the loop left its sleep running" if forked else
                                 "the loop killed a process it did not fork")
                self.assertEqual(err, "")

    def test_a_loop_stopped_the_moment_it_is_forked_still_stops(self):
        # A TERM that reached a refresh loop before it had set its own trap was lost: the
        # loop had inherited the run's TERM trap, and bash drops a signal caught that way.
        # The lock stayed the run's, so the loop went on refreshing it, and the run waited
        # for that loop for good -- a run never ended on a Ctrl-C right after a take, and
        # every later run waited on a lock kept fresh. A lock taken and abandoned at once,
        # as a signal right after the take does, over and over in six runs at a time: a
        # rare race, which this hit in every round before the fix.
        procs, deadline = [], time.monotonic() + 120
        for k in range(6):
            lock = self.dir / f"x{k}.lock"
            procs.append(self.start(
                f"i=0\nwhile [ $i -lt 1000 ]; do i=$((i+1))\n"
                f"\t_lock_take {_q(lock)} 'indexing x' {_q(self.dir / f'x{k}')} || exit 9\n"
                f"\t_lock_abandon_all\n\trm -rf {_q(lock)}\ndone\necho done"))
        for proc in procs:
            rc, out, err = self.result(proc, stop_after=max(deadline - time.monotonic(), 1))
            self.assertEqual((rc, out), (0, "done\n"),
                             "a run waited for good on a loop that missed its TERM")
            self.assertEqual(err, "")

    def test_the_refresh_stops_when_the_lock_is_no_longer_this_runs(self):
        # taken over meanwhile -- the holder was frozen for 15 minutes, say: the new
        # holder's record is not this run's to refresh
        proc = self.hold_then("sleep 30")
        self.assertRefreshed(proc)
        theirs = owner_record(4242, host="another-host")
        (self.lock / "owner").write_text(theirs)
        _wait_for(lambda: _subshells(proc.pid) == [], 10, "the loop to stop")
        _age(self.lock / "owner", 60)
        time.sleep(1)
        self.assertGreater(time.time() - (self.lock / "owner").stat().st_mtime, 3000)
        self.assertEqual((self.lock / "owner").read_text(), theirs)


# A container's PID 1 that reaps only its own child, as the web server does (its
# subprocess waits for the run it started): as a subreaper, it is handed every process
# that a run leaves behind, and never waits for one. Runs argv[1:], then writes each
# process it was handed to $ORPHANS ("<pid> <state> <name>"; Z for a zombie) and clears
# them up.
_NON_REAPING_PID1 = r"""
import ctypes, os, signal, subprocess, sys, time
if ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
    sys.exit("prctl(PR_SET_CHILD_SUBREAPER) failed")
rc = subprocess.call(sys.argv[1:])
time.sleep(0.5)  # for what the run stopped last to end
left = []
for pid in filter(str.isdigit, os.listdir("/proc")):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            stat = fh.read()
    except OSError:
        continue
    state, ppid = stat[stat.rindex(")") + 2:].split()[:2]
    if int(ppid) == os.getpid():
        left.append((int(pid), state, stat[stat.index("(") + 1:stat.rindex(")")]))
with open(os.environ["ORPHANS"], "w") as fh:
    fh.writelines(f"{pid} {state} {name}\n" for pid, state, name in left)
for pid, _, _ in left:
    try:
        os.kill(pid, signal.SIGKILL)
        os.waitpid(pid, 0)
    except OSError:
        pass
sys.exit(rc if rc >= 0 else 128 - rc)
"""


class BuildLocksInTheShippedScript(_ShippedScript):
    """The same locks, taken by whole runs of the shipped script. A reference-only
    search with no index installed builds 20-NRG_2_hg38 (since 2.9.0 a name that says the
    guide length and the side) under genome_library/20-NRG_2_hg38.lock, in a folder of that
    lock (_build_in_lock); STUB_HOLD_BUILD keeps that build running until the test kills or
    interrupts the run."""

    def setUp(self):
        super().setUp()
        self.index = self.library / "20-NRG_2_hg38"
        self.lock = self.library / "20-NRG_2_hg38.lock"

    def start_build(self):
        """A run that is building 20-NRG_2_hg38, and so holds its lock."""
        run = self.launch(variant=False, env={"STUB_HOLD_BUILD": "1"})
        _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                  "the build to start")
        return run

    def assertBuiltAndSearched(self, run):
        """`run` rebuilt 20-NRG_2_hg38 itself, searched it, and released its lock."""
        self.assertIsNone(run.reference, run.stdout)  # never taken for a finished index
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))
        self.assertFalse(self.lock.exists(), "the run left its lock behind")
        # nothing failed before the reference search, where the stub stops the run
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed",
                         run.stderr)

    def test_a_lock_left_by_a_killed_build_is_taken_over(self):
        killed = self.start_build()
        # It had written part of the index when it was killed outright, as by the OOM
        # killer or a container stop: nothing of it ran to release the lock.
        self.index.mkdir()
        (self.index / "chr1.bin").write_bytes(b"")
        os.killpg(killed.proc.pid, signal.SIGKILL)
        killed.proc.wait()
        self.assertTrue(self.lock.is_dir())
        run = self.run_pipeline(variant=False, timeout=60)
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {killed.proc.pid}",
                      run.stdout)
        self.assertBuiltAndSearched(run)

    def test_a_lock_left_in_another_container_is_taken_over_once_stale(self):
        # Each `docker run` has a hostname and PID namespace of its own, so the next
        # container cannot check this PID -- but the run that held the lock stopped
        # refreshing it when its container went away, 16 minutes ago.
        self.lock.mkdir()
        (self.lock / "owner").write_text(owner_record(4242, host="3f9c2a1b7d10",
                                                      pidns="pid:[4026532123]"))
        _age(self.lock / "owner", 16)
        self.index.mkdir()  # its half-written index
        (self.index / "chr1.bin").write_bytes(b"")
        run = self.run_pipeline(variant=False, timeout=60)
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid 4242 on "
                      "3f9c2a1b7d10, started 2026-10-30T00:00:00Z) has not refreshed its "
                      "lock for 15 minutes.", run.stdout)
        self.assertBuiltAndSearched(run)

    def test_the_lock_is_kept_fresh_while_building_and_let_go_after(self):
        built, searching = self.root / "built", self.root / "searched"
        run = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built),
                                              "STUB_SEARCH_UNTIL": str(searching),
                                              "_lock_beat_s": "0.2"})
        _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                  "the build to start")
        owner = self.lock / "owner"
        self.assertEqual(len(_subshells(run.proc.pid)), 1, _processes(run.proc.pid))
        _age(owner, 60)
        _wait_for(lambda: time.time() - owner.stat().st_mtime < 30, 10, "a refresh")
        built.touch()  # the build ends, and the run releases its lock and searches
        _wait_for(lambda: any(c[0] == "search" for c in run.calls()), 60, "the search")
        self.assertFalse(self.lock.exists())
        self.assertEqual(_subshells(run.proc.pid), [])  # the refresh stopped with it
        searching.touch()
        result = self.finish(run, timeout=60)
        self.assertEqual(_wait_until_gone(run.proc.pid), [])  # nothing left running
        self.assertEqual(result.stderr.strip(), "ERROR: Off-targets search Reference failed",
                         result.stderr)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "a process handed a run's orphans is a Linux subreaper")
    def test_a_run_leaves_no_zombie_under_a_pid_1_that_does_not_reap(self):
        # In the web container the server is PID 1, and it reaps only the runs it starts.
        # A refresh loop that was stopped left its sleep unreaped, so every lock a run
        # took and let go left one zombie sleep there for the life of the container.
        # With nothing installed, the run takes and lets go of all four locks.
        orphans = self.root / "orphans.txt"
        run = self.finish(self.launch(env={"STUB_ENRICH": "1", "ORPHANS": str(orphans)},
                                      wrap=(sys.executable, "-c", _NON_REAPING_PID1)),
                          timeout=120)
        self.assertTrue(run.enriched)
        self.assertEqual([c[1] for c in run.builds],
                         ["hg38", f"hg38+{self.VCF}", f"hg38+{self.VCF}_INDELS/NRG_2_fakechr1"])
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed",
                         run.stderr)
        self.assertEqual(orphans.read_text(), "")

    def test_a_live_build_is_waited_for_and_taken_over_when_it_dies(self):
        holder = self.start_build()
        waiter = self.launch(variant=False, env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60,
                  "the second run to wait")
        time.sleep(1)
        # the live holder's lock is honoured: no second build, its record untouched
        self.assertEqual([c for c in waiter.calls() if c[0] == "index-genome"], [])
        self.assertIn(f"\npid={holder.proc.pid}\n", (self.lock / "owner").read_text())
        # and the wait names it, without the how-to-clear that would break its lock
        self.assertRegex(waiter.stdout(), rf"Still waiting for {re.escape(str(self.lock))}: "
                         rf"it is held by pid {holder.proc.pid} on .+, which is still "
                         r"running\. Wait for that run to finish, or stop it\.\n")
        self.assertNotIn("rm -rf", waiter.stdout())
        os.killpg(holder.proc.pid, signal.SIGKILL)
        holder.proc.wait()
        run = self.finish(waiter, timeout=60)
        self.assertIn(f"Stale lock {self.lock}", run.stdout)
        self.assertBuiltAndSearched(run)

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a stopped run is told from a live one through /proc")
    def test_a_stopped_build_on_this_host_is_waited_for_however_long_it_is_stopped(self):
        # Ctrl-Z, kill -STOP, a SLURM suspend: the run is stopped, its refresh loop with it,
        # and its record ages. A run on this host took the lock over after 15 minutes and
        # rebuilt under the stopped build, which went on writing once resumed.
        built = self.root / "built"
        holder = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built)})
        _wait_for(lambda: any(c[0] == "index-genome" for c in holder.calls()), 60,
                  "the build to start")
        os.killpg(holder.proc.pid, signal.SIGSTOP)
        _wait_for(lambda: _state(holder.proc.pid) == "T", 10, "the run to stop")
        record = (self.lock / "owner").read_text()
        _age(self.lock / "owner", 16)  # stopped for 16 minutes: never refreshed meanwhile
        waiter = self.launch(variant=False, env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the second run to wait")
        time.sleep(1)  # several polls
        self.assertEqual([c for c in waiter.calls() if c[0] == "index-genome"], [])
        self.assertEqual((self.lock / "owner").read_text(), record)
        # bash leads its process group here, so the run is resumed, or ended, as a group
        self.assertRegex(waiter.stdout(), rf"Still waiting for {re.escape(str(self.lock))}: "
                         rf"it is held by pid {holder.proc.pid} on .+, which is stopped\. Resume "
                         rf"it with fg or kill -CONT -- -{holder.proc.pid}, or end it with kill "
                         rf"-KILL -- -{holder.proc.pid}\.\n")
        self.assertNotIn("rm -rf", waiter.stdout())
        self.assertEqual(waiter.stderr(), "")
        # Resumed, it is the same run, however old its record -- whose refresh waits for the
        # loop's next beat -- whether bash alone is resumed (kill -CONT <pid>: its build step
        # and refresh loop stay stopped) or the whole run (fg, or kill -CONT -- -<pgid>).
        for resume in (lambda: os.kill(holder.proc.pid, signal.SIGCONT),
                       lambda: os.killpg(holder.proc.pid, signal.SIGCONT)):
            resume()
            _wait_for(lambda: _state(holder.proc.pid) != "T", 10, "the run to resume")
            time.sleep(1)  # several polls
            self.assertEqual([c for c in waiter.calls() if c[0] == "index-genome"], [])
            self.assertEqual((self.lock / "owner").read_text(), record)
        # its build ends, it lets go, and the waiter searches what it built
        self.index.mkdir()
        (self.index / "chr1.bin").write_bytes(b"")
        built.touch()
        run = self.finish(waiter, timeout=60)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed",
                         run.stderr)
        self.assertEqual(self.finish(holder, timeout=60).searched(variant=False),
                         ("20-NRG_2_hg38", False))

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a stopped run is told from a live one through /proc")
    def test_a_stopped_build_in_a_group_it_does_not_lead_is_offered_no_group_signal(self):
        # Launched by crisprme.py from a shell, or by the web server, the run's bash is in a
        # process group it does not lead -- the job's, the server's -- which a signal to the
        # group would stop or kill as well. And resumed alone (kill -CONT <pid>, as the
        # notice said), bash waits on for its build step and refresh loop, still stopped:
        # the run never finished. Stood for by a launcher that leads the group.
        built = self.root / "built"
        launcher = "import subprocess, sys; sys.exit(subprocess.call(sys.argv[1:]))"
        holder = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built)},
                             wrap=(sys.executable, "-c", launcher))
        _wait_for(lambda: any(c[0] == "index-genome" for c in holder.calls()), 60,
                  "the build to start")
        bash = int(re.search(r"^pid=(\d+)$", (self.lock / "owner").read_text(), re.M).group(1))
        self.assertEqual(os.getpgid(bash), holder.proc.pid)  # its launcher's group
        os.killpg(holder.proc.pid, signal.SIGSTOP)  # the whole job: Ctrl-Z, a SLURM suspend
        _wait_for(lambda: _state(bash) == "T", 10, "the run to stop")
        waiter = self.launch(variant=False, env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the second run to wait")
        self.assertRegex(waiter.stdout(), rf"Still waiting for {re.escape(str(self.lock))}: "
                         rf"it is held by pid {bash} on .+, which is stopped\. It must be resumed "
                         r"together with its stopped child processes, as fg does: kill -CONT "
                         rf"{bash} resumes the script alone\.\n")
        for advice in ("-- -", "kill -KILL", "rm -rf"):
            self.assertNotIn(advice, waiter.stdout())
        self.assertEqual(waiter.stderr(), "")
        # kill -CONT <pid> resumes the script alone: its build step stays stopped
        os.kill(bash, signal.SIGCONT)
        _wait_for(lambda: _state(bash) != "T", 10, "bash to resume")
        # bash's own child only: a child the group stop caught between its fork and its
        # exec shows the build step's command line too
        build = [p for p, ppid, args in _processes(holder.proc.pid)
                 if ppid == bash and "crispritz.py index-genome" in args]
        self.assertEqual([_state(p) for p in build], ["T"])
        # resumed with it, the run finishes, and the waiter searches what it built
        os.killpg(holder.proc.pid, signal.SIGCONT)
        self.index.mkdir()
        (self.index / "chr1.bin").write_bytes(b"")
        built.touch()
        run = self.finish(waiter, timeout=60)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))
        self.assertEqual(self.finish(holder, timeout=60).searched(variant=False),
                         ("20-NRG_2_hg38", False))

    @unittest.skipUnless(os.path.exists("/proc/self/stat"),
                         "a dead run that was never reaped is told from a live one through /proc")
    def test_a_build_whose_run_died_unreaped_is_taken_over_at_once(self):
        # Killed outright and never reaped -- its parent killed first, and the container's
        # PID 1 does not reap -- the run is a zombie, which answers kill -0. Its record,
        # refreshed a moment ago, is a dead run's all the same.
        killed = self.start_build()
        self.index.mkdir()
        (self.index / "chr1.bin").write_bytes(b"")
        os.killpg(killed.proc.pid, signal.SIGKILL)  # and no wait(): left unreaped
        _wait_for(lambda: _state(killed.proc.pid) == "Z", 10, "the run to die")
        run = self.run_pipeline(variant=False, timeout=30)
        self.assertRegex(run.stdout, rf"Stale lock {re.escape(str(self.lock))}: the run that "
                         rf"held it \(pid {killed.proc.pid} on \S+, started \S+\) is gone\. ")
        self.assertBuiltAndSearched(run)

    def test_a_run_whose_lock_was_taken_over_meanwhile_leaves_it_to_that_run(self):
        # Stopped for 15 minutes where no other run can see it stopped -- in another
        # container, on another host -- a run has its lock taken over, and the taker builds
        # again. Resumed, its own build step ended: with success, it searched the index
        # being rebuilt; with a failure (its directory removed under it), it removed the
        # taker's build. Now it ends with an ERROR and touches neither.
        theirs = owner_record(4242, host="3f9c2a1b7d10", pidns="pid:[4026532123]")
        for fails in (False, True):
            with self.subTest(build_fails=fails):
                for left in (self.lock, self.index):
                    shutil.rmtree(left, ignore_errors=True)
                built = self.root / f"built{int(fails)}"
                run = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built),
                                                      "STUB_BUILD_FAILS": str(int(fails))})
                _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                          "the build to start")
                os.killpg(run.proc.pid, signal.SIGSTOP)
                _wait_for(lambda: _state(run.proc.pid) == "T", 10, "the run to stop")
                # what the run in another container leaves, once it has taken over
                (self.lock / "owner").write_text(theirs)
                self.index.mkdir()
                (self.index / "chr1.bin").write_bytes(b"theirs")
                os.killpg(run.proc.pid, signal.SIGCONT)
                built.touch()
                result = self.finish(run, timeout=60)
                self.assertEqual(run.proc.returncode, 1)
                self.assertEqual(result.stderr.splitlines()[-1],
                                 f"ERROR: the lock on {self.index} was taken over by another "
                                 "run (pid 4242 on 3f9c2a1b7d10, started 2026-10-30T00:00:00Z) "
                                 "while this run held it; that run builds it again, and this run "
                                 "uses none of what it built. Re-run once that run is done.")
                self.assertEqual(len(result.stderr.splitlines()), 1 + fails, result.stderr)
                self.assertIsNone(result.searched(variant=False))
                self.assertEqual((self.lock / "owner").read_text(), theirs)
                self.assertEqual(os.listdir(self.index), ["chr1.bin"])
                self.assertEqual((self.index / "chr1.bin").read_bytes(), b"theirs")
                self.assertEqual(_wait_until_gone(run.proc.pid), [])

    # the record of a run in another container that has taken a lock over
    THEIRS = owner_record(4242, host="3f9c2a1b7d10", pidns="pid:[4026532123]")

    def take_over_while_stopped(self, run, lock, theirs, release):
        """Stop `run` in its build step, and do meanwhile what a run in another container
        does once it has taken `lock` over: its record, and `theirs()` for what it has
        built so far. Then resume `run` and end its build step, by creating `release`."""
        os.killpg(run.proc.pid, signal.SIGSTOP)
        _wait_for(lambda: _state(run.proc.pid) == "T", 10, "the run to stop")
        (lock / "owner").write_text(self.THEIRS)
        theirs()
        os.killpg(run.proc.pid, signal.SIGCONT)
        release.touch()

    def assertLeftToTheTaker(self, run, output, lock, build_errors):
        """`run` ended in failure, the ERROR of a lost lock on `output` its last line on
        stderr after `build_errors` of its build step's, searched nothing, left the taker's
        record as it was, and left nothing running."""
        result = self.finish(run, timeout=60)
        self.assertEqual(run.proc.returncode, 1, result.stdout[-2000:] + result.stderr)
        self.assertEqual(result.stderr.splitlines()[-1],
                         f"ERROR: the lock on {output} was taken over by another run (pid 4242 "
                         "on 3f9c2a1b7d10, started 2026-10-30T00:00:00Z) while this run held "
                         "it; that run builds it again, and this run uses none of what it "
                         "built. Re-run once that run is done.")
        self.assertEqual(len(result.stderr.splitlines()), 1 + build_errors, result.stderr)
        self.assertEqual([c for c in result.calls if c[0] == "search"], [])
        self.assertEqual((lock / "owner").read_text(), self.THEIRS)
        self.assertEqual(_wait_until_gone(run.proc.pid), [])

    def test_a_run_whose_enrichment_lock_was_taken_over_meanwhile_leaves_it_to_that_run(self):
        # The enrichment site, as the reference index one above: resumed with its lock lost,
        # the run moved its enriched genome and fake-indel genomes into the taker's.
        enriched = self.home / "Genomes" / f"hg38+{self.VCF}"
        indels = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS"
        lock = Path(f"{enriched}.lock")
        for fails in (False, True):
            with self.subTest(enrichment_fails=fails):
                for left in (enriched, indels, lock):
                    shutil.rmtree(left, ignore_errors=True)
                done = self.root / f"enriched{int(fails)}"
                run = self.launch(env={"STUB_ENRICH": "1", "STUB_ENRICH_UNTIL": str(done),
                                       "STUB_ENRICH_FAILS": str(int(fails))})
                _wait_for(lambda: any(c[0] == "add-variants" for c in run.calls()), 60,
                          "the enrichment to start")

                def theirs():
                    enriched.mkdir()
                    (enriched / "chr1.fa").write_text("theirs")
                    indels.mkdir()
                    (indels / "theirs").write_text("")
                self.take_over_while_stopped(run, lock, theirs, done)
                self.assertLeftToTheTaker(run, enriched, lock, fails)
                self.assertEqual(os.listdir(enriched), ["chr1.fa"])
                self.assertEqual((enriched / "chr1.fa").read_text(), "theirs")
                self.assertEqual(os.listdir(indels), ["theirs"])
                self.assertEqual(list((self.home / "Genomes").glob("run_*")), [])  # its scratch
                self.assertFalse(any(c[0] == "index-genome" for c in run.calls()))

    def test_a_run_whose_variant_index_lock_was_taken_over_meanwhile_leaves_it_to_that_run(self):
        # resumed with its lock lost, a failed build removed the taker's variant index
        self.install("NRG_2", variant=False)  # a reference index; the variant one is built
        var = self.library / f"20-NRG_2_hg38+{self.VCF}"
        lock = Path(f"{var}.lock")
        for fails in (False, True):
            with self.subTest(build_fails=fails):
                for left in (var, lock, self.home / "Genomes" / f"hg38+{self.VCF}",
                             self.home / "Genomes" / f"hg38+{self.VCF}_INDELS"):
                    shutil.rmtree(left, ignore_errors=True)
                built = self.root / f"built{int(fails)}"
                run = self.launch(env={"STUB_ENRICH": "1", "STUB_BUILD_UNTIL": str(built),
                                       "STUB_BUILD_FAILS": str(int(fails))})
                _wait_for(lambda: ["index-genome", f"hg38+{self.VCF}"] in
                          [c[:2] for c in run.calls()], 60, "the variant index build")

                def theirs():
                    var.mkdir(exist_ok=True)
                    (var / "chr1.bin").write_bytes(b"theirs")
                self.take_over_while_stopped(run, lock, theirs, built)
                self.assertLeftToTheTaker(run, var, lock, fails)
                self.assertEqual(os.listdir(var), ["chr1.bin"])
                self.assertEqual((var / "chr1.bin").read_bytes(), b"theirs")

    def test_a_run_whose_indels_lock_was_taken_over_meanwhile_leaves_it_to_that_run(self):
        # resumed with its lock lost, a failed build removed the taker's INDELS index
        self.install("NRG_3")  # a downloaded variant index; its INDELS companion built here
        indels = self.library / f"NRG_3_hg38+{self.VCF}_INDELS"
        lock = Path(f"{indels}.lock")
        source = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS" / f"fake_{self.VCF}_chr1"
        source.mkdir(parents=True)
        (source / "fakechr1.fa").write_text(">fakechr1\nACGTACGTACGTACGTACGTACGTACGT\n")
        for fails in (False, True):
            with self.subTest(build_fails=fails):
                for left in (indels, lock):
                    shutil.rmtree(left, ignore_errors=True)
                built = self.root / f"built{int(fails)}"
                run = self.launch(env={"STUB_BUILD_UNTIL": str(built),
                                       "STUB_BUILD_FAILS": str(int(fails))})
                _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                          "the INDELS build")

                def theirs():
                    (indels / "NRG_3_fakechr1").mkdir(parents=True, exist_ok=True)
                    (indels / "NRG_3_fakechr1" / "fakechr1.bin").write_bytes(b"theirs")
                self.take_over_while_stopped(run, lock, theirs, built)
                self.assertLeftToTheTaker(run, indels, lock, fails)
                self.assertEqual(os.listdir(indels), ["NRG_3_fakechr1"])
                self.assertEqual((indels / "NRG_3_fakechr1" / "fakechr1.bin").read_bytes(),
                                 b"theirs")

    def test_index_path_names_a_dead_builds_partial_index(self):
        # --index-path never builds, and the scan skips a locked index. A run that died
        # building one there used to leave every later --index-path run failing with
        # "re-download it" and no word of the lock -- which nothing ever cleared.
        self.index.mkdir()  # as complete as it looks
        (self.index / "chr1.bin").write_bytes(b"")
        for record, aged, gone in (
                (owner_record(dead_pid()), 0, "is gone"),
                (owner_record(4242, host="3f9c2a1b7d10"), 16,
                 "has not refreshed its lock for 15 minutes")):
            with self.subTest(gone=gone):
                self.lock.mkdir()
                (self.lock / "owner").write_text(record)
                _age(self.lock / "owner", aged)
                run = self.run_pipeline(variant=False, index_path=self.library, timeout=60)
                err = run.stderr.splitlines()
                self.assertRegex(err[0], rf"^ERROR: {re.escape(str(self.index))} is a partial "
                                         rf"index: the run that was building it \(pid \d+ on "
                                         rf"\S+, started \S+\) {gone}\. Clear it, and its "
                                         rf"lock, with: rm -rf '{re.escape(str(self.lock))}' "
                                         rf"'{re.escape(str(self.index))}'$")
                self.assertEqual(err[1], "ERROR: no matching reference index under "
                                         f"--index-path '{self.library}'")
                self.assertEqual(len(err), 4, run.stderr)  # then how to get the index again
                self.assertEqual(run.builds, [])
                self.assertEqual((self.lock / "owner").read_text(), record)  # left as it was
                shutil.rmtree(self.lock)

    def test_index_path_waits_for_a_live_build_and_then_uses_it(self):
        # A run building in the working library -- the usual --index-path -- is waited
        # for, as by a run that builds, instead of failing at once.
        built = self.root / "built"
        holder = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built)})
        _wait_for(lambda: any(c[0] == "index-genome" for c in holder.calls()), 60,
                  "the build to start")
        waiter = self.launch(variant=False, index_path=self.library,
                             env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the second run to wait")
        self.assertIn(f"Another run is indexing {self.index} — waiting...", waiter.stdout())
        self.assertIsNone(waiter.proc.poll())
        self.assertEqual(waiter.stderr(), "")
        self.index.mkdir()  # the build finishes ...
        (self.index / "chr1.bin").write_bytes(b"")
        built.touch()  # ... and its run lets go of the lock
        run = self.finish(waiter, timeout=60)
        self.assertEqual(run.reference, "20-NRG_2_hg38")
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))
        self.assertEqual(run.builds, [])
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed",
                         run.stderr)
        self.finish(holder, timeout=60)

    def test_index_path_waits_on_a_lock_it_cannot_check(self):
        # no owner record: as when building, waited on and explained, never broken
        self.index.mkdir()
        (self.index / "chr1.bin").write_bytes(b"")
        self.lock.mkdir()
        waiter = self.launch(variant=False, index_path=self.library,
                             env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the notice")
        self.assertIn(f"Still waiting for {self.lock}: it has no owner record", waiter.stdout())
        self.assertIn(f"rm -rf '{self.lock}' '{self.index}'", waiter.stdout())
        self.assertEqual(waiter.stderr(), "")
        os.killpg(waiter.proc.pid, signal.SIGKILL)
        waiter.proc.wait(timeout=30)
        self.assertEqual(list(self.lock.iterdir()), [])

    def test_index_path_ignores_the_lock_of_an_index_that_could_not_serve(self):
        # a dead run's NGG build in the same library is no part of an NRG search
        other = self.library / "NGG_2_hg38.lock"
        other.mkdir()
        (other / "owner").write_text(owner_record(dead_pid()))
        run = self.run_pipeline(variant=False, index_path=self.library, timeout=60)
        self.assertTrue(run.stderr.startswith("ERROR: no matching reference index under "
                                              "--index-path"), run.stderr)
        self.assertNotIn("NGG_2_hg38", run.stderr)
        self.assertNotIn("waiting", run.stdout)

    def test_index_path_names_a_dead_lock_exactly_where_its_index_could_serve(self):
        # The resolver's own width rule, both ways: an all-N index serves at the PAM's
        # width, and a wider one only under --prune-pam with a 3' PAM. A dead run's lock
        # on one that could serve guards a partial index, which the error names; one on
        # an index that could not is no part of the search. (The side it is built for
        # is in its .pam_build, written once its build is over, so a lock is judged by
        # width alone: the finished copy here is built for the search PAM's side.)
        for pam, prefix, prune, serves in (
                ("20bp-NRG-SpCas9", "NNN_2", 0, True),
                ("20bp-NRG-SpCas9", "NNNNNNNN_2", 0, False),
                ("20bp-NRG-SpCas9", "NNNNNNNN_2", 1, True),
                ("20bp-NRG-SpCas9", "NN_2", 1, False),
                ("23bp-TTTV-Cas12a", "NNNN_2", 0, True),
                ("23bp-TTTV-Cas12a", "NNNNNNNN_2", 1, False)):
            with self.subTest(pam=pam, index=prefix, prune=prune):
                for left in self.library.glob("*.lock"):
                    shutil.rmtree(left)
                index = self.library / f"{prefix}_hg38"
                # the resolver, given that index finished, serves this search from it or not
                finished = Path(tempfile.mkdtemp(dir=self.root))
                (finished / index.name).mkdir()
                (finished / index.name / "chr1.bin").write_bytes(b"")
                seq, pos = (_PAMS / f"{pam}.txt").read_text().split()
                pos = int(pos)
                want = seq[-pos:] if pos > 0 else seq[:-pos]
                flank = prefix.split("_")[0]
                built_for(finished / index.name, len(flank) if pos > 0 else -len(flank))
                self.assertEqual(resolve(finished, want, "hg38", 2, prune=prune, pos=pos),
                                 index.name if serves else "")
                lock = Path(f"{index}.lock")
                lock.mkdir()
                (lock / "owner").write_text(owner_record(dead_pid()))
                run = self.run_pipeline(pam=pam, prune=prune, variant=False,
                                        index_path=self.library, timeout=60)
                self.assertEqual(f"ERROR: {index} is a partial index" in run.stderr, serves,
                                 run.stderr)
                self.assertIn("ERROR: no matching reference index under --index-path",
                              run.stderr)
                self.assertNotIn("waiting", run.stdout)

    def test_index_path_waits_only_for_a_build_that_could_serve(self):
        live = subprocess.Popen(["sleep", "120"])
        self.addCleanup(lambda: (live.kill(), live.wait()))
        lock = self.library / "NNNNNNNN_2_hg38.lock"
        lock.mkdir()
        (lock / "owner").write_text(owner_record(live.pid))
        # without --prune-pam an 8-N flank never serves an NRG search: no wait for that
        # build, which would end in the same error
        run = self.run_pipeline(variant=False, index_path=self.library, timeout=20,
                                env={"_lock_poll_s": "0.2"})
        self.assertNotIn("waiting", run.stdout)
        self.assertTrue(run.stderr.startswith("ERROR: no matching reference index under "
                                              "--index-path"), run.stderr)
        # under --prune-pam it does: the run waits for it
        waiter = self.launch(variant=False, prune=1, index_path=self.library,
                             env={"_lock_poll_s": "0.2"})
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the wait")
        self.assertIn(f"Another run is indexing {self.library / 'NNNNNNNN_2_hg38'} — "
                      "waiting...", waiter.stdout())
        self.assertEqual(waiter.stderr(), "")
        os.killpg(waiter.proc.pid, signal.SIGKILL)
        waiter.proc.wait(timeout=30)
        self.assertEqual((lock / "owner").read_text(), owner_record(live.pid))

    def test_a_locked_indels_index_that_cannot_be_rebuilt_says_whose_lock_it_is(self):
        # A precomputed variant index whose INDELS companion is locked, with no source
        # genome to rebuild it from: not "missing", but a dead run's partial build, or a
        # live run's build in progress.
        self.install("NRG_3")
        indels = self.library / f"NRG_3_hg38+{self.VCF}_INDELS"
        lock = Path(f"{indels}.lock")
        source = f"hg38+{self.VCF}_INDELS"
        live = subprocess.Popen(["sleep", "60"])
        self.addCleanup(lambda: (live.kill(), live.wait()))
        for record, says in (
                (owner_record(dead_pid()),
                 rf"^ERROR: INDELS index '{re.escape(str(indels))}' is a partial build: the "
                 rf"run that was building it \(pid \d+ on \S+, started \S+\) is gone, and it "
                 rf"cannot be rebuilt here \(source genome '{re.escape(source)}' is not "
                 rf"installed\)\. Clear it, and its lock, with rm -rf "
                 rf"'{re.escape(str(lock))}' '{re.escape(str(indels))}', then re-download "
                 rf"it: crisprme\.py download --what index --index-name "
                 rf"NRG_3_hg38\+{self.VCF}$"),
                (owner_record(live.pid),
                 rf"^ERROR: INDELS index '{re.escape(str(indels))}' is locked, and it "
                 rf"cannot be built here \(source genome '{re.escape(source)}' is not "
                 rf"installed\): it is held by pid {live.pid} on \S+ \(a run started "
                 r"\S+\), which is still running\. Wait for that run to finish, or stop "
                 r"it\.$")):
            with self.subTest(record=record):
                lock.mkdir()
                (lock / "owner").write_text(record)
                run = self.run_pipeline(timeout=60)
                self.assertRegex(run.stderr.strip(), says)
                self.assertEqual((lock / "owner").read_text(), record)
                shutil.rmtree(lock)

    def test_a_dead_runs_enrichment_is_redone_from_scratch(self):
        # The enrichment site: a run died holding Genomes/<ref>+<vcf>.lock after it had
        # moved part of its output into place -- the enriched genome, a fake-indel
        # genome. The takeover clears both before it enriches again: a `mv` into an
        # existing directory nests the new genome inside the old one, or fails on stderr.
        self.install("NRG_2", variant=False)  # a reference index; no variant index
        enriched = self.home / "Genomes" / f"hg38+{self.VCF}"
        indels = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS"
        fake = indels / f"fake_{self.VCF}_chr1"
        lock = Path(f"{enriched}.lock")
        lock.mkdir()
        (lock / "owner").write_text(owner_record(dead_pid()))
        enriched.mkdir()
        (enriched / "chr1.fa").write_text(">chr1\nPARTIAL\n")
        fake.mkdir(parents=True)
        (fake / "fakechr1.fa").write_text(">fakechr1\nPARTIAL\n")
        run = self.run_pipeline(timeout=60, env={"STUB_ENRICH": "1"})
        self.assertIn(f"Stale lock {lock}: the run that held it (pid ", run.stdout)
        self.assertEqual([c[0] for c in run.calls].count("add-variants"), 1)
        # this run's output, each in its place: nothing nested, nothing of the dead run's
        self.assertEqual([p.name for p in enriched.iterdir()], ["chr1.fa"])
        self.assertEqual((enriched / "chr1.fa").read_text(), ">chr1\nACGT\n")
        self.assertEqual([p.name for p in indels.iterdir()], [fake.name])
        self.assertEqual([p.name for p in fake.iterdir()], ["fakechr1.fa"])
        self.assertNotIn("PARTIAL", (fake / "fakechr1.fa").read_text())
        dicts = self.home / "Dictionaries" / f"dictionaries_{self.VCF}"
        self.assertEqual([p.name for p in dicts.iterdir()], ["my_dict_chr1.json"])
        self.assertFalse(lock.exists())
        self.assertEqual(list((self.home / "Genomes").glob("run_*")), [])  # its scratch
        # and the new genome went on to be indexed, its indels too, and searched
        self.assertIn(["index-genome", f"hg38+{self.VCF}"], [c[:2] for c in run.calls])
        self.assertIn(["index-genome", f"hg38+{self.VCF}_INDELS/NRG_2_fakechr1"],
                      [c[:2] for c in run.calls])
        self.assertEqual(run.searched(variant=True), (f"20-NRG_2_hg38+{self.VCF}", False))
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed",
                         run.stderr)

    def test_an_index_directory_with_no_bin_and_no_lock_is_built(self):
        # The scan takes a directory with no .bin for no index (a download or a build
        # that never wrote one), and the build replaces it. A lock found free with its
        # "output" in place must not count such a directory as built.
        self.index.mkdir()
        run = self.run_pipeline(variant=False, timeout=60, env={"STUB_BUILD_WRITES": "1"})
        self.assertBuiltAndSearched(run)
        self.assertEqual(sorted(p.name for p in self.index.iterdir()), [".pam_build", "chr1.bin"])

    def test_a_build_leaves_nothing_under_its_name_and_touches_no_index_of_the_older_one(self):
        # The engine writes an index under the older name, <PAM>_<N>_<genome>, and removes
        # a folder of that name first. It runs in a folder of the build lock, and its index
        # is renamed into place, under the name of 2.9.0, once whole: an index of the older
        # name in place (here one of the other side, which 2.8.1 stopped on, as in the way)
        # is never touched, and a build killed partway leaves its part in the lock, which
        # the run that takes the lock over clears.
        older = self.library / "NRG_2_hg38"
        older.mkdir()
        (older / "chr1.bin").write_bytes(b"theirs")
        (older / ".pam_build").write_text("23bp-NRG-5p 23 -3\n")
        built = self.root / "built"
        run = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built)})
        _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                  "the build to start")
        (scratch,) = self.lock.glob("build.*")
        part = scratch / "genome_library" / "NRG_2_hg38"  # what the engine writes there
        part.mkdir(parents=True)
        (part / "chr1.bin").write_bytes(b"half")
        os.killpg(run.proc.pid, signal.SIGKILL)  # an OOM kill, a container stop
        run.proc.wait()
        self.assertFalse(self.index.exists())
        self.assertEqual((older / "chr1.bin").read_bytes(), b"theirs")
        after = self.run_pipeline(variant=False, timeout=60, env={"STUB_BUILD_WRITES": "1"})
        self.assertIn(f"Stale lock {self.lock}: the run that held it (pid {run.proc.pid} ",
                      after.stdout)
        self.assertBuiltAndSearched(after)
        self.assertEqual((self.index / "chr1.bin").read_bytes(), b"")  # the stub's, not "half"
        self.assertEqual((older / "chr1.bin").read_bytes(), b"theirs")
        self.assertEqual((older / ".pam_build").read_text(), "23bp-NRG-5p 23 -3\n")
        self.assertEqual(sorted(p.name for p in self.library.iterdir()),
                         ["20-NRG_2_hg38", "NRG_2_hg38"])

    def test_a_lock_with_no_owner_record_is_waited_on_never_broken(self):
        # an older CRISPRme's lock: whether its run is still alive cannot be told
        self.lock.mkdir()
        waiter = self.launch(variant=False)
        _wait_for(lambda: "Still waiting for" in waiter.stdout(), 60, "the notice")
        out = waiter.stdout()
        self.assertIn(f"Still waiting for {self.lock}: it has no owner record", out)
        self.assertIn(f"rm -rf '{self.lock}' '{self.index}'", out)
        self.assertEqual([c for c in waiter.calls() if c[0] == "index-genome"], [])
        self.assertEqual(waiter.stderr(), "")
        os.killpg(waiter.proc.pid, signal.SIGKILL)  # stop it while it waits
        waiter.proc.wait(timeout=30)
        self.assertEqual(list(self.lock.iterdir()), [])  # not its lock: left as it was

    # Ctrl-C reaches the whole run; so do a stop sent to all of its processes (SLURM's
    # scancel or time limit) and a terminal's hangup
    INTERRUPTS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)

    def interrupt(self, run, sig, partial, lock, quiet=True, again=False):
        """Signal the whole of `run` once the build it holds `lock` for has written
        `partial`, one .bin of its index, and check that the run left the lock as it
        held it. `quiet`: nothing reaches stderr but bash's word for a step the signal
        killed ("Terminated", "Hangup"; none for Ctrl-C). `again`: signal once more if
        the first leaves the step hung, as a user presses Ctrl-C again."""
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.write_bytes(b"")
        record = (lock / "owner").read_text()
        os.killpg(run.proc.pid, sig)
        try:
            run.proc.wait(timeout=10 if again else 30)
        except subprocess.TimeoutExpired:
            if not again:
                raise
            os.killpg(run.proc.pid, sig)
            run.proc.wait(timeout=30)
        self.assertEqual(run.proc.returncode, -sig)  # still ended by the signal
        self.assertEqual((lock / "owner").read_text(), record,
                         "an interrupted build let go of its lock")
        self.assertEqual(_wait_until_gone(run.proc.pid), [])  # its refresh stopped too
        if quiet:
            self.assertEqual([line for line in run.stderr().splitlines()
                              if not re.search(r"\b(Terminated|Hangup)\b", line)], [])

    def test_an_interrupted_build_leaves_its_lock_and_the_next_run_redoes_it(self):
        # Letting go of the lock would leave the half-written index looking finished, and
        # the next run would search it: fewer off-targets, no error.
        for sig in self.INTERRUPTS:
            with self.subTest(signal=sig.name):
                run = self.start_build()
                self.interrupt(run, sig, self.index / "chr1.bin", self.lock)
                after = self.run_pipeline(variant=False, timeout=60)
                self.assertIn(f"Stale lock {self.lock}: the run that held it (pid "
                              f"{run.proc.pid} ", after.stdout)
                self.assertBuiltAndSearched(after)  # never "Reference index present"
                shutil.rmtree(self.index)

    def test_a_hangup_or_term_to_bash_alone_waits_for_the_build_and_leaves_the_lock(self):
        # `kill -HUP <pid>` or `kill <pid>`, or a hangup that reaches bash alone.
        # Untrapped, either ran the EXIT trap at once and ended bash mid-build, its build
        # step left writing as an orphan into what the next run takes over and rebuilds.
        # Trapped, the run waits for that step, then leaves the lock for the next run to
        # redo the build. TERM is trapped again after each refresh loop's fork.
        for sig in (signal.SIGHUP, signal.SIGTERM):
            with self.subTest(signal=sig.name):
                built = self.root / "built"
                run = self.launch(variant=False, env={"STUB_BUILD_UNTIL": str(built)})
                _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                          "the build to start")
                record = (self.lock / "owner").read_text()
                os.kill(run.proc.pid, sig)  # bash alone
                time.sleep(1)
                self.assertIsNone(run.proc.poll(), "bash did not wait for its build step")
                self.assertEqual((self.lock / "owner").read_text(), record)
                built.touch()  # the build step ends
                run.proc.wait(timeout=30)
                self.assertEqual(run.proc.returncode, -sig)  # ended by the signal
                self.assertEqual((self.lock / "owner").read_text(), record)
                self.assertEqual(_wait_until_gone(run.proc.pid), [])
                self.assertEqual(run.stderr(), "")
                after = self.run_pipeline(variant=False, timeout=60)
                self.assertIn(f"Stale lock {self.lock}: the run that held it (pid "
                              f"{run.proc.pid} ", after.stdout)
                self.assertBuiltAndSearched(after)
                built.unlink()  # the stub builds nothing, so this is all there is to undo

    def test_an_interrupted_variant_index_build_leaves_its_lock_and_the_next_run_redoes_it(self):
        self.install("NRG_2", variant=False)  # a reference index; the variant one is built
        var = self.library / f"20-NRG_2_hg38+{self.VCF}"
        lock = Path(f"{var}.lock")
        for sig in self.INTERRUPTS:
            with self.subTest(signal=sig.name):
                run = self.launch(env={"STUB_ENRICH": "1", "STUB_HOLD_BUILD": "1"})
                _wait_for(lambda: ["index-genome", f"hg38+{self.VCF}"] in
                          [c[:2] for c in run.calls()], 60, "the variant index build")
                self.interrupt(run, sig, var / "chr1.bin", lock)
                after = self.run_pipeline(timeout=60)
                self.assertIn(f"Stale lock {lock}: the run that held it (pid {run.proc.pid} ",
                              after.stdout)
                # never "Variant index present", nor taken for a precomputed one
                self.assertIsNone(after.variant, after.stdout)
                self.assertFalse(after.precomputed, after.stdout)
                self.assertIn(["index-genome", f"hg38+{self.VCF}"], [c[:2] for c in after.calls])
                self.assertEqual(after.searched(variant=True), (var.name, False))
                self.assertFalse(lock.exists())
                self.assertEqual(after.stderr.strip(), "ERROR: Off-targets search Variant failed",
                                 after.stderr)
                shutil.rmtree(var)

    def test_an_interrupted_indels_build_leaves_its_lock_and_the_next_run_redoes_it(self):
        # a downloaded variant index whose INDELS companion is built here, from its source
        self.install("NRG_3")
        indels = self.library / f"NRG_3_hg38+{self.VCF}_INDELS"
        shutil.rmtree(indels)
        lock = Path(f"{indels}.lock")
        source = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS" / f"fake_{self.VCF}_chr1"
        source.mkdir(parents=True)
        (source / "fakechr1.fa").write_text(">fakechr1\nACGTACGTACGTACGTACGTACGTACGT\n")
        for sig in self.INTERRUPTS:
            with self.subTest(signal=sig.name):
                run = self.launch(env={"STUB_HOLD_BUILD": "1"})
                _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                          "the INDELS build")
                # pool_index_indels.py may say on stderr how it was stopped. On Ctrl-C it
                # can also hang, as at b5d6cc1: its pool's terminate() can kill a worker
                # that holds the result queue's lock while handing back its result.
                # Another Ctrl-C ends it.
                self.interrupt(run, sig, indels / "NRG_3_fakechr1" / "fakechr1.bin", lock,
                               quiet=False, again=sig == signal.SIGINT)
                after = self.run_pipeline(timeout=60)
                self.assertIn(f"Stale lock {lock}: the run that held it (pid {run.proc.pid} ",
                              after.stdout)
                self.assertIsNone(after.indels, after.stdout)  # never "Indels index present"
                self.assertIn(["index-genome", f"hg38+{self.VCF}_INDELS/NRG_3_fakechr1"],
                              [c[:2] for c in after.calls])
                self.assertFalse(lock.exists())
                self.assertEqual(after.stderr.strip(), "ERROR: Off-targets search Variant failed",
                                 after.stderr)
                shutil.rmtree(indels)

    def test_a_run_interrupted_as_it_takes_a_lock_over_leaves_it_to_the_next(self):
        # Ctrl-C as soon as a takeover had begun (".takeover" appeared) left ".takeover" in
        # the lock, 40 times out of 40, and no run ever took that lock over again. One that
        # lands as the takeover's mkdir returns may still leave it, before the run counts it
        # as its own; the next run clears it once it is a minute old (aged here for that).
        takeover, hits = self.lock / ".takeover", 0
        for _ in range(20):
            for left in (self.lock, self.index):
                shutil.rmtree(left, ignore_errors=True)
            self.lock.mkdir()
            (self.lock / "owner").write_text(owner_record(dead_pid()))
            self.index.mkdir()
            (self.index / "chr1.bin").write_bytes(b"")
            run = self.launch(variant=False, env={"STUB_HOLD_BUILD": "1",
                                                  "_lock_beat_s": "0.2"})
            calls, deadline = run.work / "engine_calls.txt", time.monotonic() + 30
            while not (takeover.is_dir() or calls.exists() or run.proc.poll() is not None
                       or time.monotonic() > deadline):
                pass
            if not takeover.is_dir():  # missed: the takeover was over before it was seen
                _reap(run.proc)
                continue
            os.killpg(run.proc.pid, signal.SIGINT)
            try:
                run.proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                # Late (this test descheduled, under load): the takeover was over when
                # the Ctrl-C came ("Stale lock" is out), and bash 5.2 can lose a trapped
                # signal as it expands a $(...) -- the line after the takeover expands
                # $(date), and the script's note on its traps says they stay best-effort
                # there -- so the run went on to its build. Missed, as above.
                if f"Stale lock {self.lock}: " not in (run.work / "stdout.txt").read_text():
                    raise
                _reap(run.proc)
                continue
            self.assertEqual(run.proc.returncode, -signal.SIGINT)
            self.assertTrue((self.lock / "owner").is_file())  # the lock is left, not let go
            self.assertEqual(_wait_until_gone(run.proc.pid), [])
            hits += 1
            if hits == 5:
                break
        self.assertGreaterEqual(hits, 3, "the takeover was too quick to be interrupted")
        if takeover.is_dir():
            _age(takeover, 2)
        after = self.run_pipeline(variant=False, timeout=60)
        self.assertIn(f"Stale lock {self.lock}: the run that held it", after.stdout)
        self.assertBuiltAndSearched(after)

    def test_a_locked_indels_index_is_not_taken_as_present(self):
        # An INDELS build that died left its lock beside its half-written directory. The
        # presence check looked at the directory alone, so every later run searched the
        # partial index; with the source genome installed it is now rebuilt.
        self.install("NRG_3")
        indels = self.library / f"NRG_3_hg38+{self.VCF}_INDELS"
        lock = Path(f"{indels}.lock")
        lock.mkdir()
        (lock / "owner").write_text(owner_record(dead_pid()))
        source = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS" / f"fake_{self.VCF}_chr1"
        source.mkdir(parents=True)
        (source / "fakechr1.fa").write_text(">fakechr1\nACGTACGTACGTACGTACGTACGTACGT\n")
        run = self.run_pipeline(timeout=60)
        self.assertIsNone(run.indels, run.stdout)
        self.assertIn(f"Stale lock {lock}", run.stdout)
        self.assertIn(["index-genome", f"hg38+{self.VCF}_INDELS/NRG_3_fakechr1"],
                      [c[:2] for c in run.calls])
        self.assertFalse(lock.exists())
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed",
                         run.stderr)


def _mem_available_mb():
    """MemAvailable in MB, as the guard reads it; 0 where there is none."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return 0


@unittest.skipUnless(_mem_available_mb() > 0,
                     "the RAM thread guard sizes a search from MemAvailable in /proc/meminfo")
class AVariantSearchWithNoBulges(_ShippedScript):
    """With no bulges the variant search scans the enriched genome. A download ships no
    enriched genome (Genomes/<ref>+<vcf>): crispritz.py said so on stdout and exited 0, and
    every variant-only target was lost, with no error. The installed index answers it now,
    as a search with bulges does; a whole enriched genome is still scanned."""

    def assertSearchedIndex(self, run, prefix="NRG_3", prune=False):
        call = next(c for c in run.calls if c[0] == "search" and "-var" in c)
        self.assertEqual(Path(call[1]).name, f"{prefix}_hg38+{self.VCF}")
        self.assertIn("-index", call)
        self.assertNotIn("-r", call)
        self.assertEqual(call[call.index("-bDNA") + 1], "0")
        self.assertEqual(call[call.index("-bRNA") + 1], "0")
        self.assertEqual(call[call.index("--max-edits") + 1], "6")
        self.assertEqual("--prune-pam" in call, prune)

    def enriched(self, *fastas):
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir(parents=True, exist_ok=True)
        for name in fastas:
            (folder / name).write_text(">chr\nACGT\n")
        return folder

    def test_with_no_enriched_genome_it_searches_the_installed_index(self):
        self.install("NRG_3")
        self.assertSearchedIndex(self.run_pipeline(bulges=0))

    def test_a_pam_less_stand_in_answers_it_too(self):
        self.install("NNN_3")
        self.assertSearchedIndex(self.run_pipeline(bulges=0), "NNN_3")

    def test_a_wider_flank_answers_it_under_prune(self):
        self.install("NNNNNNNN_3")
        self.assertSearchedIndex(self.run_pipeline(bulges=0, prune=1), "NNNNNNNN_3", prune=True)

    def test_an_empty_enriched_genome_is_no_genome_to_scan(self):
        self.install("NRG_3")
        self.enriched()
        self.assertSearchedIndex(self.run_pipeline(bulges=0))

    def test_an_enriched_genome_lacking_a_chromosome_is_not_scanned(self):
        (self.home / "Genomes" / "hg38" / "chr2.fa").write_text(">chr2\nACGT\n")
        self.install("NRG_3")
        self.enriched("chr1.enriched.fa")
        self.assertSearchedIndex(self.run_pipeline(bulges=0))

    def test_a_whole_enriched_genome_is_still_scanned(self):
        self.install("NRG_3")
        folder = self.enriched("chr1.enriched.fa")
        run = self.run_pipeline(bulges=0)
        calls = [c for c in run.calls if c[0] == "search"]
        self.assertEqual(calls[1][1].rstrip("/"), str(folder))
        self.assertIn("-r", calls[1])
        self.assertNotIn("-index", calls[1])

    def test_a_search_with_bulges_goes_on_taking_the_index(self):
        self.install("NRG_3")
        self.enriched("chr1.enriched.fa")
        run = self.run_pipeline(bulges=1)
        call = next(c for c in run.calls if c[0] == "search" and "-var" in c)
        self.assertEqual(Path(call[1]).name, f"NRG_3_hg38+{self.VCF}")
        self.assertIn("-index", call)


class AnEngineThatStopsBeforeItSearchesFailsTheRun(_ShippedScript):
    """crispritz.py says what stops it before it searches on stdout and exits 0, so a
    search that never ran looked like one that found nothing: the run went on with the
    other search's rows alone, or ended 'No off-targets found'. It writes a targets file
    with its header alone when it does find nothing."""

    def ran_past_the_searches(self, run):
        return any("Search INDELs Start" == line or line == "Start post-analysis"
                   for line in run.stdout.splitlines())

    def test_the_reference_search_stopped(self):
        self.install("NRG_3", variant=False)
        run = self.run_pipeline(variant=False, env={"STUB_ENGINE_STOPS": "reference"})
        self.assertIn("ERROR: off-targets search Reference stopped before it searched",
                      run.stderr)
        self.assertIn("log_verbose.txt", run.stderr)
        self.assertFalse(self.ran_past_the_searches(run), run.stdout)

    def test_the_variant_search_stopped(self):
        self.install("NRG_3")
        run = self.run_pipeline(env={"STUB_ENGINE_STOPS": "variant"})
        self.assertIn("ERROR: off-targets search Variant stopped before it searched",
                      run.stderr)
        self.assertFalse(self.ran_past_the_searches(run), run.stdout)

    def test_a_search_that_finds_nothing_goes_on(self):
        self.install("NRG_3", variant=False)
        run = self.run_pipeline(variant=False, env={"STUB_FINDS_NOTHING": "1"})
        self.assertNotIn("stopped before it searched", run.stderr)
        self.assertIn("Search Reference completed", run.stdout)
        self.assertIn("Start post-analysis", run.stdout)


class ALibraryTheRunCannotWriteTo(_ShippedScript):
    """A library that takes no lock -- a read-only mount, no write permission, no inode
    left -- and lacks the index a search needs. The failed mkdir read as "another run
    built it", so the search went on with a folder that is not there: the engine says
    that on stdout and exits 0, and the run ended "No off-targets found", or without its
    variant rows. A stub mkdir stands in for such a mount, as these tests also run as
    root, for whom a folder's permissions change nothing."""

    def setUp(self):
        super().setUp()
        # (fails each mkdir of a lock; STUB_MKDIR_ONCE=<file>: only the first)
        stub = self.root / "bin" / "mkdir"
        stub.write_text(
            '#!/bin/bash\n'
            'real=$(PATH=/usr/bin:/bin command -v mkdir)\n'
            'for a in "$@"; do\n'
            '  case "$a" in *.lock)\n'
            '    if [ -n "$STUB_MKDIR_ONCE" ]; then\n'
            '      [ -e "$STUB_MKDIR_ONCE" ] && exec "$real" "$@"\n'
            '      : >"$STUB_MKDIR_ONCE"\n'
            '    fi\n'
            '    echo "mkdir: cannot create directory: Read-only file system" >&2\n'
            '    exit 1 ;;\n'
            '  esac\n'
            'done\n'
            'exec "$real" "$@"\n')
        stub.chmod(0o755)

    def assertRefused(self, run, index):
        self.assertRegex(run.stderr, rf"^ERROR: index '{re.escape(index)}' is not installed, "
                                     r"and this run cannot build it: nothing can be created in ")
        self.assertIn(str(self.library), run.stderr)
        self.assertEqual(len(run.stderr.splitlines()), 1, run.stderr)  # the web shows one line
        self.assertEqual(run.builds, [])
        self.assertFalse([c for c in run.calls if c[0] == "search"], run.calls)

    def test_a_reference_search_with_bulges_stops(self):
        self.assertRefused(self.run_pipeline(variant=False), "20-NRG_2_hg38")
        self.assertIn("--index-path", self.run_pipeline(variant=False).stderr)

    def test_a_variant_search_stops_at_the_missing_variant_index(self):
        # the worse shape: the reference index is there, so only the variant rows went
        self.install("NRG_2", variant=False)
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir()
        (folder / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        run = self.run_pipeline()
        self.assertRefused(run, f"20-NRG_2_hg38+{self.VCF}")
        # --index-path does not serve a variant index
        self.assertNotIn("--index-path", run.stderr)

    def test_a_variant_search_stops_at_the_missing_indels_index(self):
        # (it fell back to a bulge-less scan of the fake-indel genome)
        self.install("NRG_2")
        shutil.rmtree(self.library / f"NRG_2_hg38+{self.VCF}_INDELS")
        fakes = self.home / "Genomes" / f"hg38+{self.VCF}_INDELS"
        (fakes / f"fake_{self.VCF}_chr1").mkdir(parents=True)
        self.assertRefused(self.run_pipeline(), f"NRG_2_hg38+{self.VCF}_INDELS")

    def test_a_search_of_an_index_in_place_needs_no_lock(self):
        self.install("NRG_2", variant=False)
        run = self.run_pipeline(variant=False)
        self.assertEqual(run.reference, "NRG_2_hg38")
        self.assertEqual(run.searched(variant=False), ("NRG_2_hg38", False))
        self.assertNotIn("is not installed", run.stderr)

    def test_a_0_bulge_search_opens_no_index_and_goes_on(self):
        run = self.run_pipeline(bulges=0, variant=False)
        self.assertNotIn("is not installed", run.stderr)
        calls = [c for c in run.calls if c[0] == "search"]
        self.assertEqual(len(calls), 1, run.calls)
        self.assertEqual(Path(calls[0][1]).name, "hg38")

    def test_a_0_bulge_variant_search_of_the_enriched_genome_goes_on(self):
        self.install("NRG_2", variant=False)
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir()
        (folder / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        (self.home / "Genomes" / f"hg38+{self.VCF}_INDELS").mkdir()
        run = self.run_pipeline(bulges=0)
        self.assertNotIn("is not installed", run.stderr)
        self.assertEqual([Path(c[1]).name for c in run.calls if c[0] == "search"],
                         ["hg38", f"hg38+{self.VCF}"])

    def test_a_lock_a_holder_lets_go_of_as_this_run_looks_is_taken_on_the_second_try(self):
        # mkdir fails on a lock that is gone by the next look: its holder failed and let
        # go, and the build is this run's. Not a library that takes no lock.
        run = self.run_pipeline(variant=False, env={
            "STUB_MKDIR_ONCE": str(self.root / "once"), "STUB_BUILD_WRITES": "1"})
        self.assertNotIn("is not installed", run.stderr)
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("20-NRG_2_hg38", False))


class NoVariantIndexAndNothingToBuildOneFrom(_ShippedScript):
    """A download holds its variant index, and neither VCFs nor an enriched genome: asked
    for another PAM or a deeper budget the search tried to enrich, and failed in 'mv: cannot
    stat', with nothing that named the PAM, the indexes there or what to do; on the web too,
    which launches this script straight."""

    def refused(self, **kw):
        launched = self.launch(**kw)
        run = self.finish(launched)
        self.assertEqual(launched.proc.returncode, 1, run.stdout + run.stderr)
        self.assertEqual(run.calls, [])
        return run

    def engine_cannot_prune(self):
        """The stub engine as a CRISPRitz without --prune-pam: its source never names it."""
        (self.root / "bin" / "crispritz.py").write_text(
            _STUB_ENGINE.replace("--prune-pam", "the flag"))

    def assertNothingMade(self):
        self.assertEqual(sorted(p.name for p in (self.home / "Genomes").iterdir()), ["hg38"])
        self.assertEqual(sorted(p.name for p in (self.home / "Dictionaries").iterdir()),
                         [f"log_indels_{self.VCF}", f"registry_{self.VCF}"])

    def test_it_names_what_was_asked_what_is_installed_and_what_to_do(self):
        # (an NGG index, which does not hold NRG: an NRG one holds NGG, and serves it
        # pruned, ADedicatedIndexServesAPamItHolds in test_defaults_and_index_flex)
        self.install("NGG_3")
        run = self.refused()
        self.assertEqual(
            run.stderr.splitlines()[0],
            "ERROR: no variant index here serves this search (genome hg38, dataset "
            f"{self.VCF}, PAM NRG, 1 DNA and 1 RNA bulge(s), a 20-nt guide), and none can be "
            f"built: neither its VCFs (VCFs/{self.VCF}) nor its enriched genome (Genomes/hg38+"
            f"{self.VCF}) is installed. Variant indexes installed for {self.VCF}: NGG_3_hg38+"
            f"{self.VCF} (up to 2 bulge(s) of each type). To install one: crisprme.py download "
            f"--what index --index-name NRG_3_hg38+{self.VCF} if it is published, or crisprme.py "
            "build-index-only with the dataset's VCFs (docs/PRECOMPUTED_INDEXES.md)")
        self.assertEqual(len(run.stderr.splitlines()), 1, run.stderr)
        self.assertNothingMade()

    def test_a_deeper_budget_is_told_the_bulges_the_index_holds(self):
        self.install("NRG_2")
        run = self.refused(bulges=2)
        self.assertIn(f"NRG_2_hg38+{self.VCF} holds this PAM with up to 1 bulge(s) of each type: "
                      "search with no more.", run.stderr)

    def test_a_wider_pam_less_index_is_named_where_the_engine_cannot_prune(self):
        # (where it can, the search prunes by itself: EveryScanInTheShippedScript)
        self.engine_cannot_prune()
        self.install("NNNNNNNN_3")
        run = self.refused()
        self.assertIn(f"NNNNNNNN_3_hg38+{self.VCF} serves it pruned to this PAM at leaf time "
                      "(--prune-pam), which the installed CRISPRitz cannot do: install a "
                      "CRISPRitz build that supports --prune-pam (CRISPRme pins one).", run.stderr)

    def test_not_one_that_would_not_serve_it(self):
        # too shallow for the budget, or built for the other side
        self.install("NNNNNNNN_2")
        self.assertNotIn("--prune-pam", self.refused(bulges=2).stderr)

    def test_a_pam_less_index_of_the_other_side_says_so(self):
        self.install("NNN_3")
        built_for(self.library / f"NNN_3_hg38+{self.VCF}", -3)
        run = self.refused()
        self.assertIn("is a PAM-less index for a PAM on the 5' side of the guide", run.stderr)

    def test_a_guide_no_index_fits_is_told_their_guide_lengths(self):
        self.install("NGG_3")
        built_for_guide(self.library / f"NGG_3_hg38+{self.VCF}", 20)
        run = self.refused(pam="24bp-NGG-X", pam_text="N" * 24 + "NGG 3\n")
        self.assertIn("a 24-nt guide", run.stderr)
        self.assertIn(f"NGG_3_hg38+{self.VCF} (up to 2 bulge(s) of each type, built for a 20-nt "
                      "guide)", run.stderr)
        self.assertNotIn("search with no more", run.stderr)  # (fewer bulges would not serve it)

    def test_under_index_path_it_names_the_variant_indexes_there(self):
        # (a variant index at --index-path serves the search: EveryScanInTheShippedScript)
        self.install("NGG_3", library=self.root / "shared")
        run = self.refused(index_path=self.root / "shared")
        self.assertIn(f"Variant indexes installed for {self.VCF}: NGG_3_hg38+{self.VCF} ",
                      run.stderr)
        self.assertNotIn("supplies the reference index only", run.stderr)

    def test_an_index_that_serves_with_its_dictionaries_gone_says_so(self):
        # and how to download it again: a download of an installed index brings what of it
        # is missing, the dictionaries too since 2.9.0 (FR :163), and Settings downloads
        # none that is installed
        self.install("NRG_3")
        shutil.rmtree(self.home / "Dictionaries" / f"registry_{self.VCF}")
        run = self.refused()
        self.assertIn(
            f"Its dictionaries are not installed (Dictionaries/registry_{self.VCF} or "
            f"dictionaries_{self.VCF}, and log_indels_{self.VCF}): download the index again, "
            f"which installs them (crisprme.py download --what index --index-name NRG_3_hg38+"
            f"{self.VCF}; in the web app, remove the index in Settings, then download it "
            "there). ", run.stderr)
        self.assertNotIn("--prune-pam", run.stderr)  # (it needs no flag)

    def test_the_enriched_genome_alone_can_build_it(self):
        # the build-inputs kit: no VCFs, and an index for another PAM builds from it
        self.install("NGG_3")
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir()
        (folder / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"})
        self.assertNotIn("no variant index here", run.stderr)
        self.assertIn(f"hg38+{self.VCF}", [c[1] for c in run.builds])

    def test_vcfs_alone_can_build_it(self):
        self.install("NGG_3")
        run = self.run_pipeline(env={"STUB_ENRICH": "1", "STUB_BUILD_WRITES": "1"})
        self.assertNotIn("no variant index here", run.stderr)
        self.assertTrue(run.enriched)

    def test_an_empty_vcf_folder_is_no_vcfs(self):
        self.install("NGG_3")
        (self.home / "VCFs" / self.VCF).mkdir(parents=True)
        self.assertIn("no variant index here serves this search", self.refused().stderr)

    def test_a_reference_search_is_not_held_to_it(self):
        run = self.run_pipeline(variant=False, env={"STUB_BUILD_WRITES": "1"})
        self.assertNotIn("no variant index here", run.stderr)

    def test_an_index_that_serves_searches(self):
        self.install("NRG_3")
        run = self.run_pipeline()
        self.assertNotIn("no variant index here", run.stderr)
        self.assertEqual(run.variant, f"NRG_3_hg38+{self.VCF}")


class TheSnpCodesOfASearchsOwnEnrichment(_ShippedScript):
    """build-index-only puts the SNP codes on the fake-indel genome, and writes the phased
    indel genotype store; a search that enriched a dataset itself, or built an INDELS index
    from a fake-indel genome in place, did neither: it lost every target that needs an
    indel and a SNP beside it, with no message, and every later run kept that index. The
    overlay and the store are the real modules here, on what the stub writes."""

    SEQ = "ACGTTGCAAG" * 6
    PLAIN = f">fakechr1\n{SEQ[:26]}T{SEQ[26:]}\n"
    CODED = f">fakechr1\n{SEQ[:10]}R{SEQ[11:26]}T{SEQ[26:]}\n"
    LOG = ("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
           f"chr1_0-60_1\tS2\trs1\t0.25\tchr1_26_{SEQ[25]}_{SEQ[25]}T\t0,61\t{SEQ}\n")

    def setUp(self):
        super().setUp()
        self.fake = (self.home / "Genomes" / f"hg38+{self.VCF}_INDELS" / f"fake_{self.VCF}_chr1"
                     / "fakechr1.fa")
        self.store = self.home / "Dictionaries" / f"indel_genotypes_{self.VCF}"

    def vcfs(self):
        folder = self.home / "VCFs" / self.VCF
        folder.mkdir(parents=True)
        with gzip.open(folder / "x.chr1.vcf.gz", "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
                     f"FORMAT\tS1\tS2\nchr1\t26\trs1\t{self.SEQ[25]}\t{self.SEQ[25]}T\t.\tPASS\t"
                     "AF=0.25\tGT\t0|0\t1|0\n")

    def enrichment_in_place(self, text):
        """An enrichment that came another way: its fake-indel genome holds ``text``."""
        enriched = self.home / "Genomes" / f"hg38+{self.VCF}"
        enriched.mkdir(parents=True)
        (enriched / "chr1.enriched.fa").write_text(f">chr1\n{self.SEQ[:10]}R{self.SEQ[11:]}\n")
        (self.home / "Dictionaries" / f"dictionaries_{self.VCF}").mkdir()
        (self.home / "Dictionaries" / f"log_indels_{self.VCF}" / "logchr1.txt").write_text(self.LOG)
        self.fake.parent.mkdir(parents=True)
        self.fake.write_text(text)

    def indels_built(self, run):
        return [c for c in run.builds if "_INDELS/" in c[1]]

    def test_a_search_that_enriches_codes_the_fake_genome_before_it_goes_into_place(self):
        self.install("NRG_2", variant=False)
        self.vcfs()
        shutil.rmtree(self.home / "Dictionaries" / f"log_indels_{self.VCF}")
        run = self.run_pipeline(env={"STUB_ENRICH": "snp", "STUB_BUILD_WRITES": "1"})
        self.assertTrue(run.enriched)
        self.assertEqual(self.fake.read_text().split("\n")[1], self.CODED.split("\n")[1],
                         run.stdout + run.stderr)
        self.assertIn("overlay_indel_snps: chr1 -> 1 SNP overlays", run.stdout)
        self.assertEqual(len(self.indels_built(run)), 1, run.calls)

    def test_and_writes_the_phased_indel_genotype_store(self):
        self.install("NRG_2", variant=False)
        self.vcfs()
        self.assertFalse(self.store.exists())
        run = self.run_pipeline(env={"STUB_ENRICH": "snp", "STUB_BUILD_WRITES": "1"})
        with gzip.open(self.store / "gt_indel_chr1.tsv.gz", "rt") as fh:
            self.assertIn("S2", fh.read(), run.stdout + run.stderr)
        self.assertTrue((self.store / "gt_indel_chr1.tsv.gz.done").exists())
        self.assertEqual(sorted(p.name for p in (self.home / "Genomes").iterdir()),
                         sorted(["hg38", f"hg38+{self.VCF}", f"hg38+{self.VCF}_INDELS"]))

    def test_a_store_in_place_is_kept(self):
        self.install("NRG_2", variant=False)
        self.vcfs()
        self.store.mkdir()
        (self.store / "gt_indel_chr1.tsv.gz").write_bytes(gzip.compress(b"kept\n"))
        self.run_pipeline(env={"STUB_ENRICH": "snp", "STUB_BUILD_WRITES": "1"})
        with gzip.open(self.store / "gt_indel_chr1.tsv.gz", "rt") as fh:
            self.assertEqual(fh.read(), "kept\n")

    def test_a_fake_genome_in_place_is_coded_before_its_index_is_built(self):
        # one a search before 2.8.1 enriched, or that came another way
        self.install("NRG_2")
        shutil.rmtree(self.library / f"NRG_2_hg38+{self.VCF}_INDELS")
        self.enrichment_in_place(self.PLAIN)
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"})
        self.assertEqual(self.fake.read_text().split("\n")[1], self.CODED.split("\n")[1],
                         run.stdout + run.stderr)
        self.assertEqual(len(self.indels_built(run)), 1, run.calls)

    def test_a_failed_overlay_builds_no_index_from_that_genome(self):
        # an index built after it would be blind where the overlay did not get to
        self.install("NRG_2")
        shutil.rmtree(self.library / f"NRG_2_hg38+{self.VCF}_INDELS")
        self.enrichment_in_place(self.PLAIN)
        tool = self.root / "postprocess"
        shutil.copytree(_SCRIPT.parent, tool, ignore=shutil.ignore_patterns("__pycache__", "test_*"))
        (tool / "search_indel_snp.py").write_text("import sys\nsys.exit('the overlay failed')\n")
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"}, postprocess=tool)
        self.assertIn("the overlay failed", run.stderr)
        self.assertIn(f"ERROR: indels indexing failed on {self.VCF}", run.stderr)
        self.assertEqual(self.indels_built(run), [])

    def test_one_coded_already_is_left_as_it_is(self):
        self.install("NRG_2")
        shutil.rmtree(self.library / f"NRG_2_hg38+{self.VCF}_INDELS")
        self.enrichment_in_place(self.CODED)
        before = self.fake.stat().st_mtime_ns
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"})
        self.assertEqual(self.fake.stat().st_mtime_ns, before, run.stdout + run.stderr)

    def test_an_index_in_place_is_used_and_the_genome_left_alone(self):
        self.install("NRG_2")
        self.enrichment_in_place(self.PLAIN)
        run = self.run_pipeline()
        self.assertEqual(run.indels, f"NRG_2_hg38+{self.VCF}_INDELS")
        self.assertEqual(self.fake.read_text(), self.PLAIN)

    def test_the_opt_out_leaves_all_of_it_as_it_was(self):
        off = {"CRISPRME_INDEL_SNP": "0", "STUB_ENRICH": "snp", "STUB_BUILD_WRITES": "1"}
        self.install("NRG_2", variant=False)
        self.vcfs()
        shutil.rmtree(self.home / "Dictionaries" / f"log_indels_{self.VCF}")
        self.run_pipeline(env=off)
        self.assertEqual(self.fake.read_text().split("\n")[1], self.PLAIN.split("\n")[1])
        self.assertFalse(self.store.exists())

    def test_nothing_goes_to_stderr(self):
        self.install("NRG_2", variant=False)
        self.vcfs()
        run = self.run_pipeline(env={"STUB_ENRICH": "snp", "STUB_BUILD_WRITES": "1"})
        # (the stub stops the run at its variant search, which it reaches)
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed")


class AnIndexBuiltForAnotherGuideInTheScript(_ShippedScript):
    """The script's own use of the fit: a same-named index it would build over, the
    message under --index-path, and a search it does serve."""

    def setUp(self):
        super().setUp()
        self.p23 = {"pam": "23bp-NGG-X", "pam_text": "N" * 23 + "NGG 3\n"}
        self.p21 = {"pam": "21bp-NGG-X", "pam_text": "N" * 21 + "NGG 3\n"}

    def install_built(self, prefix, guide, variant=False):
        self.install(prefix, variant=variant)
        built_for_guide(self.library / f"{prefix}_hg38", guide)
        if variant:
            built_for_guide(self.library / f"{prefix}_hg38+{self.VCF}", guide)

    def test_an_index_of_the_older_name_for_a_shorter_guide_has_its_own_built_beside_it(self):
        # the plain two-search sequence: a 20-nt search of 2.8.1 built NGG_2, the next one
        # is 23-nt. It stopped on it: "This search builds its index under that name"
        self.install_built("NGG_2", 20)
        run = self.run_pipeline(variant=False, env={"STUB_BUILD_WRITES": "1"}, **self.p23)
        self.assertNotIn("too few bases", run.stderr)
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("23-NGG_2_hg38", False))
        self.assertTrue((self.library / "NGG_2_hg38" / "chr1.bin").exists())

    def test_so_does_a_variant_index_of_the_older_name(self):
        self.install_built("NGG_2", 23)  # the reference fits
        (self.library / f"NGG_2_hg38+{self.VCF}").mkdir()
        built_for_guide(self.library / f"NGG_2_hg38+{self.VCF}", 20)
        (self.library / f"NGG_2_hg38+{self.VCF}_INDELS").mkdir()
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir()
        (folder / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        (self.home / "Genomes" / f"hg38+{self.VCF}_INDELS").mkdir()
        run = self.run_pipeline(env={"STUB_BUILD_WRITES": "1"}, **self.p23)
        self.assertNotIn("too few bases", run.stderr)
        self.assertEqual(run.reference, "NGG_2_hg38")
        self.assertEqual([c[1] for c in run.builds], [f"hg38+{self.VCF}"])
        self.assertEqual(run.searched(variant=True), (f"23-NGG_2_hg38+{self.VCF}", False))
        self.assertTrue((self.library / f"NGG_2_hg38+{self.VCF}" / "chr1.bin").exists())

    def test_an_index_of_the_2_9_0_name_built_for_another_guide_stops_the_search(self):
        # one named for a 23-nt guide whose trees were built for a 20-nt one (a rename by
        # hand) is in the way of the index this search builds under that name: refused, with
        # the name of what it is
        self.install_built("NGG_2", 20)
        (self.library / "NGG_2_hg38").rename(self.library / "23-NGG_2_hg38")
        run = self.run_pipeline(variant=False, **self.p23)
        self.assertEqual(run.stderr.splitlines()[0],
                         "ERROR: 23-NGG_2_hg38 was built for a 20-nt guide, and this search has a "
                         "23-nt guide with 1 DNA and 1 RNA bulge(s): it stores too few bases per "
                         "site, so a search of it would report wrong alignments and lose "
                         "off-targets. This search builds its index under that name: name it for "
                         f"the guide it was built for (mv '{self.library}/23-NGG_2_hg38' "
                         f"'{self.library}/20-NGG_2_hg38'), move it out of {self.library}, or run "
                         "this search from another working directory")
        self.assertEqual(run.builds, [])
        self.assertIsNone(run.searched(variant=False))
        self.assertTrue((self.library / "23-NGG_2_hg38" / "chr1.bin").exists())

    def test_another_budget_of_a_shorter_guide_is_passed_over_and_the_search_builds_its_own(self):
        self.install_built("NGG_3", 20)
        run = self.run_pipeline(variant=False, env={"STUB_BUILD_WRITES": "1"}, **self.p23)
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("23-NGG_2_hg38", False))

    def test_a_guide_the_index_serves_is_not_rebuilt_or_refused(self):
        self.install_built("NGG_3", 20)
        run = self.run_pipeline(variant=False, **self.p21)
        self.assertEqual(run.builds, [])
        self.assertEqual(run.reference, "NGG_3_hg38")
        self.assertEqual(run.searched(variant=False), ("NGG_3_hg38", False))

    def test_index_path_says_why_the_index_there_does_not_serve(self):
        self.install_built("NGG_2", 20)
        run = self.run_pipeline(variant=False, index_path=self.library, **self.p23)
        self.assertIn("no matching reference index under --index-path", run.stderr)
        self.assertIn("NGG_2_hg38 was built for a 20-nt guide, and this search has a 23-nt guide",
                      run.stderr)
        self.assertIsNone(run.searched(variant=False))

    def test_a_reference_search_with_no_bulges_goes_on(self):
        self.install_built("NGG_1", 20)
        run = self.run_pipeline(variant=False, bulges=0, **self.p23)
        self.assertNotIn("too few bases", run.stderr)
        self.assertEqual(Path(next(c for c in run.calls if c[0] == "search")[1]).name, "hg38")

    def test_a_search_of_this_very_guide_changes_nothing(self):
        self.install_built("NGG_2", 20)
        run = self.run_pipeline(pam="20bp-NGG-SpCas9", variant=False)
        self.assertEqual(run.searched(variant=False), ("NGG_2_hg38", False))
        self.assertNotIn("too few bases", run.stderr)


class TheThreadGuardSizesTheIndexSearched(_ShippedScript):
    """The RAM guard caps the search threads to what MemAvailable can hold: ~12 GB per
    loaded index + 2.2 GB a thread for a pamless search, ~4 GB + 0.8 GB for a
    constrained one. A search on an all-N index loads and walks every position,
    whatever PAM was asked for, but the guard keyed its model on the requested PAM:
    an NRG search on NNN_3 -- or, under --prune-pam, on NNNNNNNN_3 -- was sized as
    constrained and could run at enough threads to be OOM-killed. Only the thread
    count changes; what is searched does not."""

    # the guard on, and no user ceiling under it
    GUARD_ON = {"CRISPRME_NO_THREAD_GUARD": "0", "CRISPRME_MAX_SEARCH_THREADS": ""}

    def capped(self, **kw):
        """The guard's (threads, memory model) for a 100000-thread request, checked
        against what each search was handed and the guard's own arithmetic."""
        before = _mem_available_mb()
        run = self.run_pipeline(threads=100000, env=self.GUARD_ON, **kw)
        after = _mem_available_mb()
        m = re.search(r"NOTE: capping search threads 100000 -> (\d+) to fit ~\d+GB free "
                      r"RAM \((.+?)\)\. Override", run.stdout)
        self.assertIsNotNone(m, run.stdout)
        threads, kind = int(m.group(1)), m.group(2)
        handed = {c[c.index("-th") + 1] for c in run.calls if c[0] == "search"}
        self.assertEqual(handed, {str(threads)}, run.calls)
        self.assertRegex(run.stderr.strip(), r"^ERROR: Off-targets search \w+ failed$")
        self.avail = (before, after)
        return threads, kind

    def assertSized(self, threads, per_thread_mb, base_mb, searches):
        """`threads` is what the guard's formula gives for MemAvailable as read around
        the run, give or take 1 GB for the memory moving meanwhile."""
        def cap(avail):
            budget = avail - 4000 - base_mb * searches
            return 1 if budget < per_thread_mb else budget // per_thread_mb
        lo, hi = cap(min(self.avail) - 1024), cap(max(self.avail) + 1024)
        self.assertTrue(lo <= threads <= hi, (threads, lo, hi, self.avail))

    def test_an_nrg_search_on_a_pamless_stand_in_is_sized_as_pamless(self):
        self.install("NNN_3")
        for prune in (0, 1):
            with self.subTest(prune=prune):
                threads, kind = self.capped(prune=prune)
                self.assertTrue(kind.startswith("pamless (memory-heavy) index NNN_3_hg38"),
                                kind)
                self.assertSized(threads, 2200, 12000, searches=2)

    def test_a_wider_pamless_flank_under_prune_is_sized_as_pamless(self):
        self.install("NNNNNNNN_3")
        threads, kind = self.capped(prune=1)
        self.assertTrue(kind.startswith("pamless (memory-heavy) index NNNNNNNN_3_hg38"), kind)
        self.assertSized(threads, 2200, 12000, searches=2)

    def test_a_dedicated_index_is_sized_as_constrained(self):
        # NRG_3 answers even with a stand-in beside it (ADedicatedIndexIsTheFastPath)
        self.install("NRG_3", "NNN_3")
        threads, kind = self.capped(prune=0)
        self.assertEqual(kind, "standard PAM")
        self.assertSized(threads, 800, 4000, searches=2)

    def test_a_no_bulge_search_does_not_open_the_index(self):
        # it scans the genome FASTA for NRG itself, so the stand-in index changes nothing
        self.install("NNN_3", variant=False)
        threads, kind = self.capped(bulges=0, variant=False)
        self.assertEqual(kind, "standard PAM")
        self.assertSized(threads, 800, 4000, searches=1)

    def test_a_no_bulge_variant_search_with_no_enriched_genome_opens_the_index(self):
        # it has no FASTA of the variant genome to scan, so it walks the stand-in
        self.install("NNN_3")
        threads, kind = self.capped(bulges=0)
        self.assertTrue(kind.startswith(f"pamless (memory-heavy) index NNN_3_hg38+{self.VCF}"),
                        kind)
        self.assertSized(threads, 2200, 12000, searches=2)

    def test_a_no_bulge_variant_search_of_a_whole_enriched_genome_does_not(self):
        self.install("NNN_3")
        folder = self.home / "Genomes" / f"hg38+{self.VCF}"
        folder.mkdir()
        (folder / "chr1.enriched.fa").write_text(">chr1\nACGT\n")
        threads, kind = self.capped(bulges=0)
        self.assertEqual(kind, "standard PAM")
        self.assertSized(threads, 800, 4000, searches=2)


def _annotation_restore_block():
    """The shipped end-of-run block that gunzips the annotation files back, verbatim:
    from its comment to the log_error rename that follows it."""
    text = _SCRIPT.read_text()
    start = text.index("# restore an unzipped copy of the annotation files")
    return text[start:text.index("# keep log_error but no block visualization", start)]


class AFinishedRunRestoresOnlyRealAnnotations(unittest.TestCase):
    """A run with no --annotation / --gene_annotation is handed the empty placeholder
    vuoto.txt as both files, as a path (complete-search: <install>/PostProcess/vuoto.txt,
    the web: <job>/vuoto.txt). The block compared that path with the bare name, so it
    gunzipped the placeholder too, and every such run's error log got two
    "gzip: .../vuoto.txt: unknown suffix -- ignored" lines. The block runs at the very
    end of a run, which the stub-engine runs above never reach, so it is run on its
    own here, with the real gunzip."""

    def run_block(self, annotation_file, gene_proximity):
        script = (f'annotation_file={shlex.quote(str(annotation_file))}\n'
                  f'gene_proximity={shlex.quote(str(gene_proximity))}\n'
                  f'{_annotation_restore_block()}')
        return subprocess.run(["bash", "-c", script], capture_output=True, text=True)

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_the_placeholder_is_left_alone(self):
        vuoto = self.tmp / "vuoto.txt"
        vuoto.write_text("")
        out = self.run_block(vuoto, vuoto)
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertEqual(sorted(p.name for p in self.tmp.iterdir()), ["vuoto.txt"])

    def test_real_annotations_are_still_restored_beside_their_gz(self):
        import gzip
        beds = {"annotation.sorted.bed": "chr1\t0\t10\texon\n",
                "genes.sorted.bed": "chr1\t5\t20\tgene\n"}
        for name, body in beds.items():
            with gzip.open(self.tmp / f"{name}.gz", "wt") as fh:
                fh.write(body)
        out = self.run_block(self.tmp / "annotation.sorted.bed.gz",
                             self.tmp / "genes.sorted.bed.gz")
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        for name, body in beds.items():
            self.assertEqual((self.tmp / name).read_text(), body)
            self.assertTrue((self.tmp / f"{name}.gz").is_file())  # -k: the .gz stays


def _samplefile_blocks():
    """The shipped blocks that handle the samplefile ($6), verbatim: the one that copies it
    into the output folder with its trailing newline as a run starts (from its comment to
    its fi), and the loop that reads the datasets it lists. (Up to 2.8.1 the newline went
    to $6 itself, and a third block stripped it as a finished run ended.)"""
    text = _SCRIPT.read_text()
    blocks = []
    for first, last in (("# $6 is the samplefile", "\nfi\n"),
                        ("while read samples; do", "\ndone <$sampleID\n")):
        start = text.index(first)
        blocks.append(text[start:text.index(last, start) + len(last)])
    return blocks


class TheSamplesPlaceholderIsNeverWritten(_ShippedScript):
    """A search with no samples -- reference-only, or variant-aware without --samplesID
    -- is handed the install's shared empty placeholder PostProcess/vuoto.txt as its
    samplefile. The script added a newline to it as the run started and stripped it only
    at the run's very end, so every run that ended before -- no hits, an error, a Ctrl-C
    -- left one more newline in it for good. The stub ends these runs at their last
    search, before that end. The script tells the placeholder apart as that very file,
    the vuoto.txt of the PostProcess folder it is run from, so each test runs it from a
    copy of that folder: a samples list of the user's own may be called vuoto.txt too,
    and must still be read whole."""

    def setUp(self):
        super().setUp()
        self.postprocess = self.root / "PostProcess"
        shutil.copytree(_SCRIPT.parent, self.postprocess,
                        ignore=shutil.ignore_patterns("__pycache__"))
        self.vuoto = self.postprocess / "vuoto.txt"
        self.vuoto.write_bytes(b"")

    def run_blocks(self, samplefile, *blocks):
        """Run shipped blocks in one bash with `samplefile` as $6 and the variables they
        read set as the script sets them, for a run from this test's PostProcess copy."""
        out = self.root / "blocks"
        out.mkdir(exist_ok=True)
        script = (f"starting_dir={shlex.quote(str(self.postprocess))}\n"
                  'sampleID=$(realpath "$6")\n'
                  f"current_working_directory={shlex.quote(str(self.home))}\n"
                  f"output_folder={shlex.quote(str(out))}\n"
                  "logerror=$output_folder/log_error.txt\n" + "".join(blocks))
        return subprocess.run(["bash", "-c", script, "bash", *["_"] * 5, str(samplefile)],
                              capture_output=True, text=True)

    def test_a_reference_only_search_leaves_it_empty(self):
        self.install("NRG_3", variant=False)
        run = self.run_pipeline(variant=False, samples=self.vuoto,
                                postprocess=self.postprocess)
        self.assertEqual(run.searched(variant=False), ("NRG_3_hg38", False))
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Reference failed",
                         run.stderr)
        self.assertEqual(self.vuoto.read_bytes(), b"")

    def test_a_variant_search_without_samples_leaves_it_empty(self):
        self.install("NRG_3")
        run = self.run_pipeline(samples=self.vuoto, postprocess=self.postprocess)
        self.assertEqual(run.searched(variant=True), (f"NRG_3_hg38+{self.VCF}", False))
        self.assertEqual(run.stderr.strip(), "ERROR: Off-targets search Variant failed",
                         run.stderr)
        self.assertEqual(self.vuoto.read_bytes(), b"")

    def test_a_real_samplefile_is_read_from_a_copy_with_its_newline(self):
        # the newline that `while read` needs for the last line goes to a copy in the
        # output folder, and the list itself is left as it is (RV4-1), as is the
        # placeholder, even one an older run left a newline in
        add, _ = _samplefile_blocks()
        samples = self.root / "samples.txt"
        samples.write_text(f"{self.VCF}.samplesID.txt")
        self.vuoto.write_bytes(b"\n")
        for samplefile in (samples, self.vuoto):
            out = self.run_blocks(samplefile, add)
            self.assertEqual((out.returncode, out.stderr), (0, ""))
        self.assertEqual(samples.read_text(), f"{self.VCF}.samplesID.txt")
        self.assertEqual((self.root / "blocks" / ".samplesID_list_searched.txt").read_text(),
                         f"{self.VCF}.samplesID.txt\n")
        self.assertEqual(self.vuoto.read_bytes(), b"\n")

    def test_a_samples_list_of_ones_own_called_vuoto_txt_is_read_whole(self):
        # its one line has no newline, so `while read` reads it only with the one the
        # script adds; the list is left as it was, and the placeholder alone. Told apart
        # by name, the list got no newline, its dataset was never read, and the run
        # failed only at its end, on the samples its results name.
        add, read = _samplefile_blocks()
        rows = "S1\tPOP1\tSUP1\tF\nS2\tPOP1\tSUP1\tM\n"
        (self.home / "samplesIDs").mkdir()
        (self.home / "samplesIDs" / f"{self.VCF}.samplesID.txt").write_text(
            "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n" + rows)
        mine = self.root / "lists" / "vuoto.txt"
        mine.parent.mkdir()
        mine.write_text(f"{self.VCF}.samplesID.txt")
        out = self.run_blocks(mine, add, read)
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        read_in = self.root / "blocks" / ".sampleID.txt"
        self.assertEqual(read_in.read_text() if read_in.exists() else "", rows)
        self.assertEqual(mine.read_text(), f"{self.VCF}.samplesID.txt")
        self.assertEqual(self.vuoto.read_bytes(), b"")


class ShippedScriptCarriesTheFix(unittest.TestCase):
    def test_every_lock_site_goes_through_the_shared_helpers(self):
        """A lock holds its owner record now, so a bare `rmdir` of it fails -- on stderr,
        which fails the stage -- and a bare wait loop is the one that waited forever."""
        text = _SCRIPT.read_text().replace(_lock_helpers(), "")
        for bare in ('rmdir "$_lock"', 'mkdir "$_lock"', 'while [ -d "$_lock" ]'):
            self.assertNotIn(bare, text)
        # the enrichment, reference index, variant index and INDELS index builds
        self.assertEqual(len(re.findall(r'^\s*_lock_take (-bin )?"\$_lock" ', text, re.M)), 4)
        # the two genome indexes count as built only once they hold a .bin
        self.assertEqual(len(re.findall(r'^\s*_lock_take -bin "\$_lock" "indexing '
                                        r'\$idx_folder1" ', text, re.M)), 2)
        self.assertNotIn("trap ", text)  # the helpers' traps are the only ones

    def test_every_build_site_stops_once_its_lock_is_lost(self):
        """A build site checks, as soon as its build step returns and before its failure
        check (whose rm -r would remove the other run's build), that its lock is still
        its own; and lets go of it on success only while it is."""
        lines = [line.strip() for line in
                 _SCRIPT.read_text().replace(_lock_helpers(), "").splitlines()]
        checks = [i for i, line in enumerate(lines) if '_lock_lost "$_lock"' in line]
        self.assertEqual(len(checks), 4, "enrichment, reference, variant and INDELS builds")
        for i in checks:
            j = lines.index("fi", i) + 1 if lines[i].startswith("if ") else i + 1
            while not lines[j] or lines[j].startswith("#"):
                j += 1
            self.assertRegex(lines[j], r'^if \[ -s "?\$logerror"? \]; then$', lines[i])
        self.assertEqual(lines.count('_lock_release "$_lock" || exit 1'), 4)

    def test_no_bare_nnn_glob_remains(self):
        # `NNN*_*` as a glob is what matched NNNRRT and made width ordering
        # locale-dependent; every scan must go through the resolver instead.
        text = _SCRIPT.read_text()
        self.assertNotIn('NNN*_*_"${ref_name}"', text)
        self.assertNotIn('NNN*_*_"${ref_name}+${vcf_name}"', text)

    def test_all_three_scan_sites_use_the_resolver(self):
        text = _SCRIPT.read_text()
        calls = len(re.findall(r'\$\(_resolve_index_dir ', text))
        self.assertGreaterEqual(calls, 3, "reference, variant and precomputed scans")

    def test_script_is_syntactically_valid(self):
        out = subprocess.run(["bash", "-n", str(_SCRIPT)], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)

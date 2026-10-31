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
here too: ``ADeadRunsLockDoesNotWedgeTheNext`` and ``ALockIsReleasedByTheRunThatHoldsIt``
run the shipped lock helpers, and ``BuildLocksInTheShippedScript`` kills, interrupts
and races whole runs. ``TheThreadGuardSizesTheIndexSearched`` covers the RAM guard,
which must size a search by the index the resolver picked.
"""

import os
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parent / "submit_job_automated_new_multiple_vcfs.sh"
_PAMS = _SCRIPT.parent.parent / "PAMs"


def _extract_resolver():
    """Pull `_resolve_index_dir` verbatim out of the shipped script."""
    text = _SCRIPT.read_text()
    start = text.index("\t_resolve_index_dir() {")
    # the body is tab-indented one level; its closing brace is the first "\t}"
    end = text.index("\n\t}\n", start) + len("\n\t}\n")
    body = text[start:end]
    # strip the one leading tab so it is a top-level function in the harness
    return "\n".join(line[1:] if line.startswith("\t") else line
                     for line in body.splitlines())


RESOLVER = _extract_resolver()


def resolve(libdir, pam_want, suffix, bmax, locale="C", prune=0, pos=1):
    """Run the shipped resolver; return the chosen basename, or "" for none.

    prune/pos mirror the script's own `prune_pam` and `pos` (the PAM-file offset:
    > 0 is a 3' PAM). The resolver reads them to decide whether a flank WIDER than
    the PAM may be used.
    """
    script = (f"prune_pam={prune}\npos={pos}\n{RESOLVER}\n"
              f"_resolve_index_dir \"$1\" \"$2\" \"$3\" \"$4\"\n")
    env = dict(os.environ, LC_ALL=locale, LANG=locale)
    out = subprocess.run(
        ["bash", "-c", script, "bash", str(libdir), pam_want, suffix, str(bmax)],
        capture_output=True, text=True, env=env,
    )
    if out.returncode not in (0, 1):
        raise AssertionError(f"resolver errored rc={out.returncode}: {out.stderr}")
    chosen = out.stdout.strip()
    return Path(chosen).name if chosen else ""


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
        self.assertEqual(resolve(lib, "TTTV", "hg38", 3, prune=1, pos=-4), "")

    def test_equal_width_pamless_serves_either_orientation_without_prune(self):
        lib = self.lib("NNN_3_hg38")
        self.assertEqual(resolve(lib, "NRG", "hg38", 3, prune=0), "NNN_3_hg38")
        lib2 = self.lib("NNNN_3_hg38")
        self.assertEqual(resolve(lib2, "TTTV", "hg38", 3, prune=0, pos=-4), "NNNN_3_hg38")


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
        # engine, a narrower one never stored the PAM
        for installed, expected in (
            (("NNNNNNNN_3",), ""),
            (("NNN_3",), ""),
            (("NNNN_3", "NNNNNNNN_3"), "NNNN_3"),
            (("TTTV_3", "NNNN_3", "NNNNNNNN_3"), "TTTV_3"),
        ):
            lib = self.library(*installed)
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
if [ "$1" = search ] && [ -n "$STUB_SEARCH_UNTIL" ]; then
	while ! [ -e "$STUB_SEARCH_UNTIL" ]; do sleep 0.05; done
fi
if [ "$1" = search ]; then
	# Stop at the last search: the variant one, or the reference one in a
	# reference-only run. Every index decision has been made by then.
	for _a in "$@"; do [ "$_a" = -var ] && exit 1; done
	[ "$STUB_LAST" = reference ] && exit 1
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
                     index_path=None, threads=1, env=None, timeout=300):
        """Run the script the way complete-search launches it, stopping at its last
        search."""
        return self.finish(self.launch(pam, bulges, prune, variant, index_path,
                                       threads, env), timeout)

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
               index_path=None, threads=1, env=None):
        """Start the script the way complete-search launches it (argument order of
        crisprme.py's crisprme_run), in a session of its own, so that a test can
        signal the whole run as Ctrl-C or a container stop does. `env` adds to its
        environment; the RAM thread guard is off unless `env` turns it on."""
        self.runs += 1
        work = self.root / f"run{self.runs}"
        out = work / "out"
        out.mkdir(parents=True)
        pam_file = work / f"{pam}.txt"
        shutil.copy(_PAMS / f"{pam}.txt", pam_file)
        guides = work / "guides.txt"
        guides.write_text("N" * len(pam_file.read_text().split()[0]) + "\n")
        empty = work / "empty.txt"
        empty.write_text("")
        vcfs = samples = "_"
        if variant:
            vcfs, samples = work / "vcfs.txt", work / "samples.txt"
            vcfs.write_text(self.VCF + "\n")
            samples.write_text("")
        # complete-search refuses --prune-pam without --index-path; the usual value
        # is the working directory's own genome_library
        if index_path is None:
            index_path = self.library if prune else "_"
        args = [self.home / "Genomes" / "hg38", vcfs, guides, pam_file, empty, samples,
                bulges + 1, 4, bulges, bulges, 3, out, _SCRIPT.parent, threads, self.home,
                empty, "_", 1, 0, "_", "mm+bulges", "mm+bulges,mm", "False", "PASS,.",
                index_path, 6, prune]
        run_env = dict(os.environ, PATH=self.path, STUB_LOG=str(work / "engine_calls.txt"),
                       STUB_LAST="variant" if variant else "reference",
                       LC_ALL="C", LANG="C", CRISPRME_NO_THREAD_GUARD="1")
        run_env.update(env or {})
        # stderr goes where complete-search sends it: the script fails a stage
        # whenever log_error.txt is non-empty
        with open(out / "log_error.txt", "w") as err, open(work / "stdout.txt", "w") as so:
            proc = subprocess.Popen(["bash", str(_SCRIPT), *map(str, args)], cwd=work,
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
        """No installed variant index qualified: the run falls back to building one
        from source VCFs a downloaded install does not have, and searches nothing."""
        self.assertFalse(run.precomputed, run.stdout)
        self.assertTrue(run.enriched)
        self.assertIsNone(run.searched(variant=False))
        self.assertIsNone(run.searched(variant=True))
        self.assertIn("Genome enrichment failed", run.stderr)

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

    def test_a_wide_flank_alone_is_used_only_under_prune(self):
        self.install("NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(prune=1), "NNNNNNNN_3", prune=1)
        self.assertNoVariantIndex(self.run_pipeline(prune=0))

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
        self.assertEqual(run.searched(variant=False), ("NRG_2_hg38", False))
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
        self.assertEqual(plain.searched(variant=False), ("NRG_2_hg38", False))

    def test_a_different_dedicated_pam_beside_a_stand_in(self):
        self.install("NGG_3", "NNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(prune=prune), "NNN_3", prune)

    def test_a_5_prime_pam_with_a_wider_flank_only(self):
        self.install("NNNNNNNN_3")
        self.assertNoVariantIndex(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1))
        ref_only = self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1, variant=False)
        self.assertIsNone(ref_only.searched(variant=False))
        self.assertIn("no matching reference index", ref_only.stderr)

    def test_a_5_prime_pam_takes_the_equal_width_flank(self):
        self.install("NNNN_3", "NNNNNNNN_3")
        self.assertServedBy(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=1),
                            "NNNN_3", prune=1)

    def test_a_5_prime_pam_prefers_its_dedicated_index(self):
        self.install("TTTV_3", "NNNN_3", "NNNNNNNN_3")
        for prune in (0, 1):
            self.assertServedBy(self.run_pipeline(pam="23bp-TTTV-Cas12a", prune=prune),
                                "TTTV_3", prune)


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


def owner_record(pid, host=None, pidns=None, started="2026-10-30T00:00:00Z"):
    """An owner record as the script writes one into a lock; for this host by default."""
    this_host, this_pidns = _this_host()
    return (f"host={this_host if host is None else host}\n"
            f"pidns={this_pidns if pidns is None else pidns}\n"
            f"pid={pid}\nstarted={started}\n")


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


def _processes(pgid):
    """(pid, ppid, command line) of every process still running -- zombies left out --
    in process group `pgid`."""
    out = subprocess.run(["ps", "-ww", "-e", "-o", "pid=,ppid=,pgid=,stat=,args="],
                         capture_output=True, text=True).stdout
    procs = []
    for line in out.splitlines():
        pid, ppid, group, stat, *args = line.split(None, 4)
        if int(group) == pgid and not stat.startswith("Z"):
            procs.append((int(pid), int(ppid), args[0] if args else ""))
    return procs


def _subshells(pid):
    """The background subshells of bash process `pid`: forks that run its command
    line. A lock's heartbeat is one."""
    procs = _processes(os.getpgid(pid))
    me = [args for p, _, args in procs if p == pid]
    return [p for p, _, args in procs if p != pid and me and args == me[0]]


def _wait_until_gone(pgid, timeout=10):
    """Wait for every process of group `pgid` to end; what is still running then."""
    deadline = time.monotonic() + timeout
    while _processes(pgid) and time.monotonic() < deadline:
        time.sleep(0.05)
    return _processes(pgid)


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

    def sleeper(self):
        """A live process, to stand for a run that holds a lock."""
        proc = subprocess.Popen(["sleep", "60"])
        self.addCleanup(lambda: (proc.kill(), proc.wait()))
        return proc


class ADeadRunsLockDoesNotWedgeTheNext(_Locks):
    """A run that died holding a build lock -- SIGKILL, an OOM kill, a container stop,
    a reboot -- used to leave the lock for good, and every later run that needed the
    same index waited on it forever: no error, no message, a search slot held. The
    lock now names its owner. A dead owner on this host is taken over; one this run
    cannot check is waited on, never broken, and the wait is explained."""

    def notices(self, out):
        return [line for line in out.splitlines() if line.startswith("Still waiting for ")]

    def test_the_lock_helpers_were_actually_found(self):
        block = _lock_helpers()
        for part in ("_lock_take() {", "_lock_stale() {", "_lock_steal() {",
                     "_lock_release() {", "trap '_lock_release_all' EXIT"):
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
                rc, out, err = self.take(then=f"cat {_q(self.lock)}/owner; echo self=$$",
                                         stop_after=10)
                self.assertEqual(rc, 2, out + err)  # 2: taken over, the build redone
                self.assertIn(f"Stale lock {self.lock}: the run that held it (pid 4242 on "
                              f"{host or _this_host()[0]}, started 2026-10-30T00:00:00Z) "
                              "has not refreshed its lock for 15 minutes. Taking it over",
                              out)
                me = re.search(r"^self=(\d+)$", out, re.M).group(1)
                self.assertIn(f"\npid={me}\n", out)  # the record names this run now
                self.assertNotIn("waiting", out)
                self.assertFalse(self.lock.exists())  # and it released it at its exit
                self.assertEqual(err, "")

    def test_a_live_owner_this_run_may_not_signal_is_never_broken(self):
        # kill -0 fails on a live process of another user as well as on a dead one;
        # only "No such process" means gone. PID 1 is another user's to a run that is
        # not root -- and root runs the waiter as nobody here to make it so.
        self.hold(owner_record(1))
        rc, out, err = self.result(self.start(self.take_body(), as_nobody=os.geteuid() == 0),
                                   stop_after=1.5)
        self.assertEqual(rc, -signal.SIGKILL, out + err)  # still waiting when stopped
        self.assertNotIn("Stale lock", out)
        self.assertIn("a process this run may not signal (another user's)",
                      self.notices(out)[0])
        self.assertEqual((self.lock / "owner").read_text(), owner_record(1))
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
        take = f'_lock_take -bin {_q(self.lock)} "indexing {self.out}" {_q(self.out)}\nexit $?'
        self.out.mkdir()
        rc, out, err = self.result(self.start(take))
        self.assertEqual(rc, 0, out + err)  # a free lock, the index to build
        self.assertFalse(self.lock.exists())  # (released when that run exited)
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


class ALockIsReleasedByTheRunThatHoldsIt(_Locks):
    """A run releases every lock it holds however it ends -- an exit, Ctrl-C, SIGTERM
    -- and no lock it does not hold: another run's, or one taken over from it since."""

    def test_exit_releases_this_runs_locks_and_keeps_its_status(self):
        theirs = self.dir / "theirs.lock"
        theirs.mkdir()
        (theirs / "owner").write_text(owner_record(1))
        lost = self.dir / "lost.lock"
        body = (f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                f"_lock_take {_q(lost)} 'indexing b' {_q(self.dir / 'b')} || exit 9\n"
                # another run has taken this one over since: it is no longer this run's
                f"printf 'host=elsewhere\\npid=1\\n' > {_q(lost)}/owner\n"
                "exit 3")
        rc, out, err = self.result(self.start(body))
        self.assertEqual(rc, 3, out + err)  # the run's own exit status survives the trap
        self.assertFalse(self.lock.exists(), "a lock this run held outlived it")
        self.assertEqual((theirs / "owner").read_text(), owner_record(1))
        self.assertEqual((lost / "owner").read_text(), "host=elsewhere\npid=1\n")
        self.assertEqual(err, "")

    def test_a_signal_releases_this_runs_locks_and_still_ends_the_run(self):
        # Ctrl-C reaches the whole process group; `kill` reaches bash alone
        for sig, to_group in ((signal.SIGINT, True), (signal.SIGTERM, False)):
            with self.subTest(signal=sig.name):
                body = (f"_lock_take {_q(self.lock)} 'indexing a' {_q(self.out)} || exit 9\n"
                        "echo taken\nwhile :; do sleep 0.2; done")
                proc = self.start(body)
                _wait_for(lambda: "taken" in proc.out_path.read_text(), 10, "the lock")
                if to_group:
                    os.killpg(proc.pid, sig)
                else:
                    proc.send_signal(sig)
                rc, out, err = self.result(proc)
                self.assertEqual(rc, -sig, out + err)  # ended by the signal, as without a trap
                self.assertFalse(self.lock.exists(), "an interrupted run left its lock")
                self.assertEqual(err, "")

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
        self.assertEqual(sorted(sig for _, sig in traps), ["EXIT", "INT", "TERM"])
        for body, sig in traps:
            self.assertRegex(body, r"^_lock_\w+( [A-Z]+)?$", sig)


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
        # an exit, or INT/TERM, which reach the loop through the release
        for how, status in (("exit 3", 3), ("kill -s TERM $$", -signal.SIGTERM),
                            ("kill -s INT $$", -signal.SIGINT)):
            with self.subTest(how=how):
                proc = self.hold_then(how)
                self.assertRefreshed(proc)
                self.go.touch()
                rc, out, err = self.result(proc)
                self.assertEqual(rc, status, out + err)
                self.assertFalse(self.lock.exists())
                self.assertEqual(_wait_until_gone(proc.pid), [])
                self.assertEqual(err, "")
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


class BuildLocksInTheShippedScript(_ShippedScript):
    """The same locks, taken by whole runs of the shipped script. A reference-only
    search with no index installed builds NRG_2_hg38 under
    genome_library/NRG_2_hg38.lock; STUB_HOLD_BUILD keeps that build running until
    the test kills or interrupts the run."""

    def setUp(self):
        super().setUp()
        self.index = self.library / "NRG_2_hg38"
        self.lock = self.library / "NRG_2_hg38.lock"

    def start_build(self):
        """A run that is building NRG_2_hg38, and so holds its lock."""
        run = self.launch(variant=False, env={"STUB_HOLD_BUILD": "1"})
        _wait_for(lambda: any(c[0] == "index-genome" for c in run.calls()), 60,
                  "the build to start")
        return run

    def assertBuiltAndSearched(self, run):
        """`run` rebuilt NRG_2_hg38 itself, searched it, and released its lock."""
        self.assertIsNone(run.reference, run.stdout)  # never taken for a finished index
        self.assertEqual([c[1] for c in run.builds], ["hg38"])
        self.assertEqual(run.searched(variant=False), ("NRG_2_hg38", False))
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

    def test_an_index_directory_with_no_bin_and_no_lock_is_built(self):
        # The scan takes a directory with no .bin for no index (a download or a build
        # that never wrote one), and the build replaces it. A lock found free with its
        # "output" in place must not count such a directory as built.
        self.index.mkdir()
        run = self.run_pipeline(variant=False, timeout=60)
        self.assertBuiltAndSearched(run)

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

    def test_an_interrupted_build_releases_its_lock(self):
        # Ctrl-C reaches the whole run; so does a stop sent to its process group
        for sig in (signal.SIGINT, signal.SIGTERM):
            with self.subTest(signal=sig.name):
                run = self.start_build()
                self.assertTrue(self.lock.is_dir())
                os.killpg(run.proc.pid, sig)
                run.proc.wait(timeout=30)
                self.assertEqual(run.proc.returncode, -sig)  # still ended by the signal
                self.assertFalse(self.lock.exists(), "an interrupted build left its lock")
                # bash reports a build killed by SIGTERM ("... Terminated ..."), never one
                # killed by Ctrl-C; nothing else may reach stderr
                self.assertEqual([line for line in run.stderr().splitlines()
                                  if "Terminated" not in line], [])

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

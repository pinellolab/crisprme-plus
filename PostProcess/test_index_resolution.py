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
"""

import os
import re
import shutil
import subprocess
import tempfile
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

    def __init__(self, proc, calls, stderr):
        self.stdout = proc.stdout
        self.stderr = stderr
        self.calls = calls
        self.precomputed = "Using precomputed variant index" in proc.stdout
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


@unittest.skipUnless(shutil.which("python") and _gnu_realpath(),
                     "the shipped script needs python and GNU realpath on PATH")
class EveryScanInTheShippedScript(unittest.TestCase):
    """The three scans, run by the shipped script itself with a stub engine.

    The resolver tests above prove the rule. These prove that each call site --
    the precomputed-variant detection, the reference scan and the variant scan --
    applies it with the right library, suffix, bulge budget and pruning state; that
    the reference and variant searches are handed the index their scan picked; that
    the script announces the picked variant index's own INDELS companion; and that
    --prune-pam reaches those two searches exactly when asked. The stub stops the
    run at the variant search, so the indel search itself (pool_search_indels.py,
    which is handed the companion and the flag) is not reached here.
    The layout is a downloaded install: indexes and Dictionaries/ tiers, no VCFs.
    """

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
                     index_path=None):
        """Run the script the way complete-search launches it (argument order of
        crisprme.py's crisprme_run), stopping at its last search."""
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
                bulges + 1, 4, bulges, bulges, 3, out, _SCRIPT.parent, 1, self.home,
                empty, "_", 1, 0, "_", "mm+bulges", "mm+bulges,mm", "False", "PASS,.",
                index_path, 6, prune]
        log = work / "engine_calls.txt"
        env = dict(os.environ, PATH=self.path, STUB_LOG=str(log),
                   STUB_LAST="variant" if variant else "reference",
                   LC_ALL="C", LANG="C", CRISPRME_NO_THREAD_GUARD="1")
        # stderr goes where complete-search sends it: the script fails a stage
        # whenever log_error.txt is non-empty
        with open(out / "log_error.txt", "w") as err:
            proc = subprocess.run(["bash", str(_SCRIPT), *map(str, args)], cwd=work,
                                  env=env, stdout=subprocess.PIPE, stderr=err,
                                  text=True, timeout=300)
        calls = ([line.split() for line in log.read_text().splitlines()]
                 if log.exists() else [])
        return _Run(proc, calls, (out / "log_error.txt").read_text())

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


class ShippedScriptCarriesTheFix(unittest.TestCase):
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

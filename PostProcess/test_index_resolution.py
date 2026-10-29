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
"""

import os
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parent / "submit_job_automated_new_multiple_vcfs.sh"


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


class ChoiceIsLocaleIndependent(_Tree):
    """The regression this file exists for."""

    LOCALES = ("C", "POSIX", "en_US.UTF-8")

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

#!/usr/bin/env python3
"""A wide-flank PAM-less index must not land in the repo old clients read.

A PAM-less index with a 3-nt flank stores exactly what a dedicated 3-nt PAM index
stores, so any CRISPRme release can read it and ``NNN_3`` ships in the default
dataset repo. A WIDER flank cannot be read correctly by a release that predates
``--prune-pam``, so those indexes go in a separate repo that older defaults never
point at. Routing is derived from the flank width rather than kept as a list of
names, so it cannot be forgotten at publish time.

Hermetic: these are pure functions over names, no network.
"""

import os
import pathlib
import unittest

import crisprme_hf as hf


class FlankWidth(unittest.TestCase):
    def test_pamless_widths(self):
        self.assertEqual(hf.index_flank_width("NNN_3_hg38"), 3)
        self.assertEqual(hf.index_flank_width("NNNNNNNN_3_hg38"), 8)
        self.assertEqual(hf.index_flank_width("NNNN_3_hg38"), 4)

    def test_variant_index_names(self):
        self.assertEqual(hf.index_flank_width("NNN_3_hg38+hg38_HGDP1kGP"), 3)
        self.assertEqual(hf.index_flank_width("NNNNNNNN_3_hg38+hg38_HGDP1kGP"), 8)

    def test_dedicated_pams_are_not_pamless(self):
        for name in ("NRG_3_hg38", "NGG_3_hg38", "NRR_3_hg38", "TTTV_4_hg38"):
            self.assertIsNone(hf.index_flank_width(name), name)

    def test_a_pam_that_merely_starts_with_n_is_not_pamless(self):
        # PAMs/21bp-NNNRRT-SaCas9KKH.txt ships, so this name can exist
        self.assertIsNone(hf.index_flank_width("NNNRRT_3_hg38"))

    def test_publish_marker_is_stripped_first(self):
        # canonical_index_name drops the -dictless marker from the ref segment
        self.assertEqual(
            hf.index_flank_width("NNNNNNNN_3_hg38-dictless+hg38_HGDP1kGP"), 8
        )

    def test_empty_and_malformed(self):
        self.assertIsNone(hf.index_flank_width(""))
        self.assertIsNone(hf.index_flank_width("nonsense"))


class Routing(unittest.TestCase):
    def setUp(self):
        self._saved = os.environ.pop("CRISPRME_HF_REPO", None)
        self.addCleanup(self._restore)

    def _restore(self):
        os.environ.pop("CRISPRME_HF_REPO", None)
        if self._saved is not None:
            os.environ["CRISPRME_HF_REPO"] = self._saved

    def test_compat_flank_stays_in_the_default_repo(self):
        for name in ("NNN_3_hg38", "NNN_3_hg38+hg38_HGDP1kGP"):
            self.assertEqual(hf.resolve_index_repo(name), hf.DEFAULT_HF_REPO, name)

    def test_dedicated_pam_stays_in_the_default_repo(self):
        for name in ("NRG_3_hg38", "NRG_3_hg38+hg38_HGDP1kGP", "NNNRRT_3_hg38"):
            self.assertEqual(hf.resolve_index_repo(name), hf.DEFAULT_HF_REPO, name)

    def test_wide_flank_routes_to_the_general_index_repo(self):
        for name in ("NNNN_3_hg38", "NNNNNNNN_3_hg38",
                     "NNNNNNNN_3_hg38+hg38_HGDP1kGP"):
            self.assertEqual(
                hf.resolve_index_repo(name), hf.GENERAL_INDEX_HF_REPO, name
            )

    def test_the_boundary_is_exactly_three(self):
        self.assertEqual(hf._COMPAT_FLANK_WIDTH, 3)
        self.assertEqual(hf.resolve_index_repo("NNN_3_hg38"), hf.DEFAULT_HF_REPO)
        self.assertEqual(
            hf.resolve_index_repo("NNNN_3_hg38"), hf.GENERAL_INDEX_HF_REPO
        )

    def test_explicit_flag_overrides_routing(self):
        self.assertEqual(
            hf.resolve_index_repo("NNNNNNNN_3_hg38", "me/my-repo"), "me/my-repo"
        )
        self.assertEqual(hf.resolve_index_repo("NNN_3_hg38", "me/my-repo"), "me/my-repo")

    def test_env_var_overrides_routing(self):
        os.environ["CRISPRME_HF_REPO"] = "me/env-repo"
        self.assertEqual(hf.resolve_index_repo("NNNNNNNN_3_hg38"), "me/env-repo")
        self.assertTrue(hf.repo_overridden())

    def test_not_overridden_by_default(self):
        self.assertFalse(hf.repo_overridden())
        self.assertFalse(hf.repo_overridden(None))

    def test_the_two_repos_are_distinct(self):
        self.assertNotEqual(hf.DEFAULT_HF_REPO, hf.GENERAL_INDEX_HF_REPO)


class ExistenceProbeFailsSoft(unittest.TestCase):
    def test_unreadable_repo_reads_as_absent(self):
        """"Cannot tell" must behave like "not there", so download falls back."""
        self.assertFalse(
            hf._index_archive_exists("definitely/not-a-real-repo", "NNN_3_hg38")
        )


class GenotypeCompanionNoteIsAttachedToTheRightResult(unittest.TestCase):
    """The "no genotype companion" note must follow the companion fetch, not the
    reference-index fetch.

    It had been written as the `else:` of the shared-reference-index try/except, so
    it fired whenever THAT download succeeded -- telling the user per-sample Samples
    were degraded moments after the companion was installed -- and, worse, a
    genuinely missing companion was reported by nothing at all whenever the reference
    index was already present and its block skipped. Observed on a real
    `download --what index --index-name NNN_3_hg38+hg38_HGDP1kGP`.
    """

    NOTE = "has no genotypes_"

    def setUp(self):
        import ast
        self.src = pathlib.Path(hf.__file__).read_text()
        self.tree = ast.parse(self.src)

    def _note_parents(self):
        """Every statement that directly contains the note, by kind."""
        import ast
        hits = []
        for node in ast.walk(self.tree):
            if not isinstance(node, (ast.If, ast.Try)):
                continue
            for field, branch in (("body", getattr(node, "body", [])),
                                  ("orelse", getattr(node, "orelse", []))):
                for stmt in branch:
                    text = ast.get_source_segment(self.src, stmt) or ""
                    if self.NOTE in text:
                        hits.append((type(node).__name__, field,
                                     ast.get_source_segment(self.src, node) or ""))
        return hits

    def test_the_note_exists_exactly_once(self):
        self.assertEqual(self.src.count(self.NOTE), 1)

    def test_the_note_is_not_in_a_try_else(self):
        """A Try/orelse is what mis-attached it to the reference-index download."""
        for kind, field, _ in self._note_parents():
            self.assertNotEqual(
                (kind, field), ("Try", "orelse"),
                "the companion note is attached to a try/except else again",
            )

    def test_the_note_is_guarded_on_the_companion_result(self):
        """It must sit under a branch testing `got`, not the reference index."""
        parents = self._note_parents()
        self.assertTrue(parents, "note not found inside any if/try")
        self.assertTrue(
            any("got is False" in src for _, _, src in parents),
            f"note is not guarded on the companion result: {[k for k, _, _ in parents]}",
        )

    def test_an_errored_fetch_is_not_reported_as_absent(self):
        """got=None marks "errored and already reported"; only False means absent."""
        self.assertIn("got = None", self.src)
        self.assertIn("elif got is False:", self.src)

    def test_reference_index_success_message_is_still_there(self):
        self.assertIn("Also installed the shared reference index", self.src)


if __name__ == "__main__":
    unittest.main(verbosity=2)

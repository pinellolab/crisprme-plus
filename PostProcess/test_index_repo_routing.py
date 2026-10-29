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


if __name__ == "__main__":
    unittest.main(verbosity=2)

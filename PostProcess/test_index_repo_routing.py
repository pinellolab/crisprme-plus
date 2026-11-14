#!/usr/bin/env python3
"""A wide-flank PAM-less index must not land in the repo old clients read.

A PAM-less index with a 3-nt flank stores exactly what a dedicated 3-nt PAM index
stores, so any CRISPRme release can read it and ``NNN_3`` ships in the default
dataset repo. A WIDER flank cannot be read correctly by a release that predates
``--prune-pam``, so those indexes go in a separate repo that older defaults never
point at. Routing is derived from the flank width rather than kept as a list of
names, so it cannot be forgotten at publish time.

Since 2.10.0 an index of index format 2, or one named for its guide length and side, goes to
and comes from a repo of its own, lucapinello/crisprme_plus_indices, whose id this file pins
(TheRepoOfThisReleasesIndexesIsTheOneOnTheHub).

Hermetic: these are pure functions over names, no network.
"""

import os
import pathlib
import types
import unittest
from unittest import mock

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


class TheRepoOfThisReleasesIndexesIsTheOneOnTheHub(unittest.TestCase):
    """The ids of the dataset repos, as they are on the Hub. The indexes of 2.10.0 (index format
    2) are published in, and fetched from, lucapinello/crisprme_plus_indices: its owner made it
    under that name on 2026-11-14, and the name the plan had given it (PLANNED) is no repo's. The
    archives of 2.8.1 stay in lucapinello/crisprme-data and lucapinello/crisprme-data-general-index.
    The tests of the download, the catalog and the publish (test_panel_companions,
    test_legacy_index_access) read crisprme_hf.FORMAT2_HF_REPO, whatever it is; this pins it, and
    that no file of the repository names the planned repo, where a download or a publish of 2.10.0
    would find no repo."""

    NEW = "lucapinello/crisprme_plus_indices"
    PLANNED = "crisprme-indexes-" + "format2"  # (in two parts, so that this file does not name it)

    def setUp(self):
        saved = os.environ.pop("CRISPRME_HF_REPO", None)  # (a repo named would win)
        if saved is not None:
            self.addCleanup(os.environ.__setitem__, "CRISPRME_HF_REPO", saved)

    def test_the_ids(self):
        self.assertEqual(hf.FORMAT2_HF_REPO, self.NEW)
        self.assertEqual(hf.DEFAULT_HF_REPO, "lucapinello/crisprme-data")
        self.assertEqual(hf.GENERAL_INDEX_HF_REPO, "lucapinello/crisprme-data-general-index")

    def test_an_index_of_2_10_0_goes_there_and_one_of_2_8_1_to_its_own_repo(self):
        for name in ("20-NRR_3_hg38", "20-NRR_3_hg38+hg38_HGDP1kGP", "20-NRR_3_hg38+hg38_HPRC",
                     "20-NRR_3_hg38+hg38_mega", "TTTN-23_3_hg38", "22-NNGRRN_3_hg38+hg38_HGDP1kGP",
                     "20-NNNNNNNN_3_hg38+hg38_HGDP1kGP", "20-NGG_3_hg38-dictless+hg38_HGDP1kGP"):
            self.assertEqual(hf.resolve_index_repo(name), self.NEW, name)
        # (publish-index reads the format of the .bin files: one of format 2, whatever its name)
        self.assertEqual(hf.resolve_index_repo("NRR_3_hg38", index_format=2), self.NEW)
        for name in ("NRR_3_hg38", "NRR_3_hg38+hg38_HGDP1kGP", "NRG_3_hg38+hg38_HPRC",
                     "NGG_3_hg38+hg38_HGDP1kGP", "NRG_3_hg38+hg38_mega", "NNN_3_hg38"):
            self.assertEqual(hf.resolve_index_repo(name), "lucapinello/crisprme-data", name)
        self.assertEqual(hf.resolve_index_repo("NNNNNNNN_3_hg38+hg38_HGDP1kGP"),
                         "lucapinello/crisprme-data-general-index")

    def test_the_catalog_lists_it(self):
        """list_available_downloads('index'), the Settings dropdown's list, reads it alone."""
        asked = []

        class HfApi:
            def list_repo_tree(self, repo, **kw):
                asked.append((repo, kw.get("repo_type"), kw.get("path_in_repo")))
                return iter(())

        with mock.patch.object(hf, "_require_hf", return_value=types.SimpleNamespace(HfApi=HfApi)):
            self.assertEqual(hf.list_available_downloads("index"), [])
        self.assertEqual(asked, [(self.NEW, "dataset", "indexes")])

    def test_no_file_of_the_repository_names_the_planned_repo(self):
        root = pathlib.Path(__file__).resolve().parent.parent
        files = [root / n for n in ("README.md", "METHODS.md", "CHANGELOG.md", "crisprme.py",
                                    "index.py", "Dockerfile", "environment.yml")]
        for top in ("docs", "seq_script", "scripts", "PostProcess", "pages", ".github", "changelog"):
            if (root / top).is_dir():
                files += sorted(p for p in (root / top).rglob("*") if p.is_file() and p.suffix in (
                    ".md", ".py", ".tsv", ".sh", ".html", ".yml", ".txt", ".json"))
        files = [p for p in files if p.is_file() and p.resolve() != pathlib.Path(__file__).resolve()]
        self.assertGreater(len(files), 300)  # (the scan reads the docs and the code)
        named = [str(p.relative_to(root)) for p in files
                 if self.PLANNED in p.read_text(encoding="utf-8", errors="replace")]
        self.assertEqual(named, [], f"files that name the planned repo, which is no repo: {named}")
        self.assertTrue((root / "docs" / "hf" / "crisprme_plus_indices_README.md").is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)

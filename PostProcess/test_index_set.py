#!/usr/bin/env python3
"""The index set of CRISPRme 2.10.0 (the rebuild manifest, seq_script/rebuild/index_set_2.10.0.tsv)
holds against the code that builds, names, publishes and downloads its archives
(seq_script/rebuild/check_index_set.py, which the D7 readiness audit runs before any build),
and the check finds what it is for: a name build-index-only would not give, an N that is not
the build's bulges + 1, a panel without its bundle, a genotyped panel without its store, an
archive of 2.8.1 nothing replaces (NNN_3's pointer to the general index), an index the web's
2 + 2 bulges do not fit.

Hermetic: reads the repository's files only.
"""

import copy
import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_spec = importlib.util.spec_from_file_location(
    "check_index_set", os.path.join(ROOT, "seq_script", "rebuild", "check_index_set.py"))
cis = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cis)


def _row(rows, archive):
    return next(r for r in rows if r["archive"] == archive)


class TheSetHolds(unittest.TestCase):
    def setUp(self):
        self.rows = cis.read_manifest()

    def test_no_problem(self):
        self.assertEqual(cis.check(self.rows), [])

    def test_it_is_the_approved_set(self):
        indexes = sorted(r["archive"] for r in self.rows if r["kind"] in cis.INDEX_KINDS)
        hg = "+hg38_HGDP1kGP"
        self.assertEqual(indexes, sorted(
            [f"20-NRR_3_hg38{p}" for p in ("", hg, "+hg38_HPRC", "+hg38_mega")]
            + [f"{name}{p}" for p in ("", hg) for name in (
                "20-NGG_3_hg38", "20-NRG_3_hg38", "20-NNNNNNNN_3_hg38", "22-NNGRRN_3_hg38",
                "TTTN-20_3_hg38", "TTTN-23_3_hg38")]))
        self.assertFalse([a for a in indexes if a.startswith(("NNN_", "20-NNN_"))],
                         "NNN_3 is not built again for 2.10.0")
        self.assertEqual({r["bmax"] for r in self.rows if r["kind"] in cis.INDEX_KINDS}, {"3"})

    def test_each_variant_panel_ships_its_dictionaries_once(self):
        for panel, companions in (("hg38_HGDP1kGP", ("bundle", "genotypes", "indel-genotypes")),
                                  ("hg38_HPRC", ("bundle", "genotypes", "indel-genotypes")),
                                  ("hg38_mega", ("bundle",))):
            got = sorted(r["kind"] for r in self.rows if r["panel"] == panel
                         and r["kind"] not in cis.INDEX_KINDS + ("build-inputs",))
            self.assertEqual(got, sorted(companions), panel)

    def test_nnn_3_points_to_the_general_index(self):
        for rest in ("hg38", "hg38+hg38_HGDP1kGP"):
            self.assertIn(f"NNN_3_{rest}",
                          _row(self.rows, f"20-NNNNNNNN_3_{rest}")["replaces"].split(","))

    def test_the_commands_print(self):
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(cis.main(["--commands"]), 0)
        text = out.getvalue()
        self.assertLess(text.index("# build"), text.index("# publish"))
        self.assertIn("crisprme.py publish-index --index genome_library/20-NRR_3_hg38+hg38_HGDP1kGP",
                      text)


class TheCheckFindsWhatItIsFor(unittest.TestCase):
    """Each mutation of the set is caught (and the set itself passes: TheSetHolds)."""

    def setUp(self):
        self.rows = cis.read_manifest()

    def mutated(self, which, **changes):
        rows = copy.deepcopy(self.rows)
        _row(rows, which).update(changes)
        return cis.check(rows)

    def test_a_name_build_index_only_would_not_give(self):
        problems = self.mutated("TTTN-23_3_hg38", archive="23-TTTN_3_hg38")
        self.assertTrue([p for p in problems if "build-index-only names this build TTTN-23_3_hg38"
                         in p], problems)

    def test_an_n_that_is_not_the_bulges_plus_one(self):
        problems = self.mutated("20-NGG_3_hg38", bdna="1", brna="1")
        self.assertTrue([p for p in problems if "is not max(bDNA 1, bRNA 1) + 1" in p], problems)

    def test_the_web_s_two_bulges_must_fit(self):
        # a 23-nt guide on an index of N = 2 is refused at 2 + 2 (and its name says N = 2)
        problems = self.mutated("TTTN-23_3_hg38", bmax="2", bdna="1", brna="1",
                                archive="TTTN-23_2_hg38")
        self.assertTrue([p for p in problems if "does not serve its own 23-nt guide" in p],
                        problems)

    def test_a_pam_file_that_says_otherwise(self):
        problems = self.mutated("22-NNGRRN_3_hg38", spacer="21")
        self.assertTrue([p for p in problems if "holds NNGRRN 3' 22 nt" in p], problems)

    def test_a_panel_without_its_bundle(self):
        rows = [r for r in self.rows if r["archive"] != "dictionaries-bundle_hg38_HPRC"]
        self.assertIn("panel hg38_HPRC: 0 bundle archive(s), not 1", cis.check(rows))

    def test_a_genotyped_panel_without_its_store(self):
        rows = [r for r in self.rows if r["archive"] != "genotypes_hg38_HGDP1kGP"]
        self.assertIn("panel hg38_HGDP1kGP: 0 genotypes archive(s), not 1", cis.check(rows))

    def test_an_archive_of_2_8_1_nothing_replaces(self):
        problems = self.mutated("20-NNNNNNNN_3_hg38+hg38_HGDP1kGP",
                                replaces="NNNNNNNN_3_hg38+hg38_HGDP1kGP")
        self.assertTrue([p for p in problems if p.startswith("NNN_3_hg38+hg38_HGDP1kGP (2.8.1): "
                                                             "replaced by 0")], problems)

    def test_a_repo_of_2_8_1(self):
        problems = self.mutated("20-NRR_3_hg38", repo="lucapinello/crisprme-data")
        self.assertTrue([p for p in problems if "repo lucapinello/crisprme-data" in p], problems)

    def test_a_variant_index_without_its_reference_index(self):
        rows = [r for r in self.rows if r["archive"] != "22-NNGRRN_3_hg38"]
        self.assertTrue([p for p in cis.check(rows)
                         if "its reference index 22-NNGRRN_3_hg38" in p])


if __name__ == "__main__":
    unittest.main()

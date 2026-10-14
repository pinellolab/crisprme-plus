#!/usr/bin/env python3
"""Unit tests for targets_scoring (the score-a-list CLI backend).

Covers the CasOFFinder + generic-TSV parsers (both schemas, bulge-from-alignment
inference, 0-vs-1-based positions), the CFD-per-site slicing (vs calc_cfd direct,
incl. the DNA-bulge front-trim convention), the CRISPR-Bulge null rules (stubbed
scorer so the test is hermetic and needs no cbulge env), and the variant carrier
overlay against a synthetic in-memory registry built with compile_registry.

STDLIB only (+ the repo's tier0_registry for the overlay fixture). Hermetic: no
installed PAMs/index, no cbulge env.
"""
import os
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import targets_scoring as ts  # noqa: E402

G = "GAGTCCGAGCAGAAGAAGAA"  # 20-nt spacer used across tests
PAMGG = "AGG"


def _co_file(rows, header=True):
    """Write a CasOFFinder-style TSV and return its path."""
    fd, path = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        if header:
            fh.write("#Id\tBulge Type\tcrRNA\tDNA\tChromosome\tLocation\t"
                     "Direction\tMismatches\tBulge Size\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


class TestParsers(unittest.TestCase):
    def test_casoffinder_modern_schema(self):
        p = _co_file([
            ["g", "X", G + PAMGG, G + "AGG", "chr8", 1000, "+", 0, 0],
            ["g", "X", G + PAMGG, G + "TGG", "8", 2000, "-", 1, 0],
        ])
        sites = ts.parse_casoffinder(p)
        self.assertEqual(len(sites), 2)
        self.assertEqual(sites[0].chrom, "chr8")
        self.assertEqual(sites[1].chrom, "chr8")  # "8" normalized to "chr8"
        self.assertEqual(sites[0].pos, 1001)       # 0-based 1000 -> 1-based 1001
        self.assertEqual(sites[0].bulge_type, "X0")
        self.assertEqual(sites[0].strand, "+")
        os.unlink(p)

    def test_casoffinder_legacy_6col(self):
        fd, p = tempfile.mkstemp(suffix=".txt")
        with os.fdopen(fd, "w") as fh:
            # legacy: crRNA Chromosome Location DNA Direction Mismatches (no header)
            fh.write("\t".join([G + PAMGG, "chr3", "500", G + "AGG", "+", "0"]) + "\n")
        sites = ts.parse_casoffinder(p)
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].chrom, "chr3")
        self.assertEqual(sites[0].pos, 501)
        os.unlink(p)

    def test_bulge_inferred_from_alignment(self):
        # DNA bulge: '-' in crRNA (guide gap), extra base in DNA
        dna_bulge_guide = "G-" + G[1:] + PAMGG      # gap in guide
        dna_bulge_off = "GA" + G[1:] + "AGG"
        p = _co_file([["g", "", dna_bulge_guide, dna_bulge_off, "chr1", 10, "+", 0, ""]])
        s = ts.parse_casoffinder(p)[0]
        self.assertEqual(s.bulge_type, "DNA")
        self.assertEqual(s.bulge_size, 1)
        os.unlink(p)
        # RNA bulge: '-' in DNA
        rna_guide = G + PAMGG
        rna_off = G[:5] + "-" + G[5:] + "GG"
        p = _co_file([["g", "", rna_guide, rna_off, "chr1", 10, "+", 0, ""]])
        s = ts.parse_casoffinder(p)[0]
        self.assertEqual(s.bulge_type, "RNA")
        self.assertEqual(s.bulge_size, 1)
        os.unlink(p)

    def test_position_base_1(self):
        p = _co_file([["g", "X", G + PAMGG, G + "AGG", "chr8", 1000, "+", 0, 0]])
        s = ts.parse_casoffinder(p, position_base=1)[0]
        self.assertEqual(s.pos, 1000)  # already 1-based -> unchanged
        os.unlink(p)

    def test_generic_tsv(self):
        fd, p = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(fd, "w") as fh:
            fh.write("guide\tofftarget\tchr\tpos\tstrand\n")
            fh.write("\t".join([G + PAMGG, G + "AGG", "chr5", "123", "+"]) + "\n")
        sites = ts.parse_tsv(p, "guide", "offtarget", chrom_col="chr",
                             pos_col="pos", strand_col="strand", position_base=1)
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0].chrom, "chr5")
        self.assertEqual(sites[0].pos, 123)
        os.unlink(p)


class TestCFD(unittest.TestCase):
    def setUp(self):
        self.mm, self.pam = ts.get_mm_pam_scores()

    def test_perfect_match_equals_pam_term(self):
        s = ts.Site(G + "GGG", G + "GGG", bulge_type="X0")
        self.assertAlmostEqual(ts.cfd_for_site(s, self.mm, self.pam), self.pam["GG"])

    def test_noncanonical_pam_scores_low(self):
        # perfect spacer, NGA PAM -> CFD == pam['GA'] (the low non-canonical term)
        s = ts.Site(G + "GGG", G + "TGA", bulge_type="X0")
        self.assertAlmostEqual(ts.cfd_for_site(s, self.mm, self.pam), self.pam["GA"])
        self.assertLess(ts.cfd_for_site(s, self.mm, self.pam), 0.1)

    def test_matches_calc_cfd_direct_mismatch(self):
        off = G[:18] + "T" + G[19:]           # one spacer mismatch
        s = ts.Site(G + "GGG", off + "AGG", bulge_type="X0")
        expect = ts.calc_cfd(G + "GGG", (off + "AGG")[:-3], (off + "AGG")[-2:],
                             self.mm, self.pam, True)
        self.assertEqual(ts.cfd_for_site(s, self.mm, self.pam), expect)

    def test_dna_bulge_front_trim(self):
        # DNA-bulge slicing must match calc_cfd(crRNA[bs:], DNA[bs:-3], DNA[-2:])
        crrna = "-" + G + "GGG"
        dna = "A" + G + "AGG"
        s = ts.Site(crrna, dna, bulge_type="DNA", bulge_size=1)
        expect = ts.calc_cfd(crrna[1:], dna[1:-3], dna[-2:], self.mm, self.pam, True)
        self.assertEqual(ts.cfd_for_site(s, self.mm, self.pam), expect)


class TestCrisprBulgeNullRules(unittest.TestCase):
    """CRISPR-Bulge via a STUBBED scorer_runner (hermetic -- no cbulge env)."""

    def setUp(self):
        import types
        self._stub = types.ModuleType("scorer_runner")
        # echo a deterministic score == len(sg) sentinel so we can see which rows
        # were scored vs nulled; the real scorer returns [0,1] floats.
        self._stub.CRISPR_BULGE_predict_list = (
            lambda sg, off, device=None: [0.5] * len(sg)
        )
        sys.modules["scorer_runner"] = self._stub

    def tearDown(self):
        sys.modules.pop("scorer_runner", None)

    def test_scoreable_vs_nulled(self):
        sites = [
            ts.Site(G + "GGG", G + "AGG", bulge_type="X0"),                 # ok
            ts.Site(G + "GGG", G[:19] + "N" + "AGG", bulge_type="X0"),      # N -> null
            ts.Site(G + "GGG", "G-" + G[1:] + "A-G", bulge_type="DNA"),     # 2 gaps -> null
        ]
        scores = ts.crispr_bulge_for_sites(sites)
        self.assertEqual(scores[0], 0.5)
        self.assertEqual(scores[1], -1.0)
        self.assertEqual(scores[2], -1.0)

    def test_env_absent_graceful(self):
        sys.modules.pop("scorer_runner", None)
        # force import to fail by shadowing with a bad finder is overkill; instead
        # point at a name that cannot import, emulated by removing it and relying on
        # the real absence in CI (no cbulge). If a real scorer_runner imports, this
        # still returns floats; assert only that it never raises and len matches.
        sites = [ts.Site(G + "GGG", G + "AGG", bulge_type="X0")]
        out = ts.crispr_bulge_for_sites(sites)
        self.assertEqual(len(out), 1)


class TestCarrierOverlay(unittest.TestCase):
    """Overlay against a synthetic single-chromosome registry (compile_registry)."""

    def setUp(self):
        import tier0_registry as t0
        self.t0 = t0
        self.meta = {
            "S1": ("1000G", "EUR", "female"),
            "S2": ("1000G", "EUR", "male"),
            "H1": ("HGDP", "EAS", "female"),
        }
        # a common variant at pos 1005 (inside a 1000..1022 window), one rare at 1010
        self.records = [
            (1005, "A", "G", "rs1005", {"S1": "0|1", "S2": "1|1", "H1": "0|1"}),
            (1010, "C", "T", "rs1010", {"S1": "0|1"}),
        ]
        self.d = tempfile.mkdtemp()
        binp = os.path.join(self.d, "reg_chr8.bin")
        idxp = os.path.join(self.d, "reg_chr8.idx")
        t0.compile_registry(self.records, self.meta, None,
                            t0.autosomal_ploidy, binp, idxp)

    def test_overlay_counts(self):
        ov = ts._RegistryOverlay(self.d)
        # window covers 1001..1023 (23-nt aligned off, no gaps)
        s = ts.Site(G + "GGG", G + "AGG", chrom="chr8", pos=1001, strand="+")
        res = ov.overlay(s)
        ov.close()
        self.assertEqual(res["n_overlapping_variants"], 2)
        # both 1005 and 1010 fall in the window; 1005 is the common one
        self.assertIn("rs1005", res["variants"])
        self.assertIn("rs1010", res["variants"])
        self.assertGreater(float(res["max_global_AF"]), 0.0)
        self.assertGreaterEqual(int(res["max_carrier_count"]), 3)  # S1,S2,H1 @1005
        self.assertIn("1000G=", res["per_dataset_AF"])
        self.assertIn("HGDP=", res["per_dataset_AF"])

    def test_overlay_no_variants_outside_window(self):
        ov = ts._RegistryOverlay(self.d)
        s = ts.Site(G + "GGG", G + "AGG", chrom="chr8", pos=2000, strand="+")
        res = ov.overlay(s)
        ov.close()
        self.assertEqual(res["n_overlapping_variants"], 0)
        self.assertEqual(res["variants"], "")

    def test_overlay_missing_chrom_is_blank(self):
        ov = ts._RegistryOverlay(self.d)
        s = ts.Site(G + "GGG", G + "AGG", chrom="chrZ", pos=1001, strand="+")
        res = ov.overlay(s)
        ov.close()
        self.assertEqual(res["n_overlapping_variants"], 0)


class TestEndToEnd(unittest.TestCase):
    def test_run_tsv_output(self):
        import types
        stub = types.ModuleType("scorer_runner")
        stub.CRISPR_BULGE_predict_list = lambda sg, off, device=None: [0.9] * len(sg)
        sys.modules["scorer_runner"] = stub
        try:
            p = _co_file([
                ["g", "X", G + "GGG", G + "GGG", "chr8", 1000, "+", 0, 0],
                ["g", "X", G + "GGG", G + "TGA", "chr8", 2000, "+", 0, 0],
            ])
            out = tempfile.mktemp(suffix=".tsv")
            n = ts.run(p, out, casoffinder=True)
            self.assertEqual(n, 2)
            lines = open(out).read().splitlines()
            self.assertEqual(len(lines), 3)  # header + 2
            header = lines[0].split("\t")
            self.assertIn("CFD", header)
            self.assertIn("CRISPR_BULGE", header)
            # sorted by CFD desc: perfect NGG site first
            first = dict(zip(header, lines[1].split("\t")))
            self.assertEqual(first["PAM"], "GGG")
            self.assertGreater(float(first["CFD"]), float(
                dict(zip(header, lines[2].split("\t")))["CFD"]))
            os.unlink(p); os.unlink(out)
        finally:
            sys.modules.pop("scorer_runner", None)


if __name__ == "__main__":
    unittest.main(verbosity=2)

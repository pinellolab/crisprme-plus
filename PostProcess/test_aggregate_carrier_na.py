"""Aggregate (AF-only "mega") registry -> carrier/hom fields render 'NA'.

A sites-only aggregate panel (compile_registry_from_info_af, manifest
aggregation=="info_af") carries real per-dataset allele frequencies but NO
genotypes, so its carrier / homozygote counts are 0 by construction. The
population-summary companion must surface those as "NA" (unknown), never a
fabricated 0, while leaving the allele-frequency fields intact. STDLIB only.
"""
import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)
import population_summary_companion as psc  # noqa: E402


class _Reader:
    def __init__(self, manifest):
        self.manifest = manifest


class _GS:
    def __init__(self, af, defined, cf, n):
        self.allele_freq = af
        self.allele_freq_defined = defined
        self.carrier_freq = cf
        self.n_carrier = n


class _Summary:
    def __init__(self, groups):
        self.groups = groups


class TestAggregateDetector(unittest.TestCase):
    def test_info_af_is_aggregate(self):
        self.assertTrue(psc._is_aggregate_af_only(_Reader({"aggregation": "info_af"})))

    def test_genotyped_is_not_aggregate(self):
        self.assertFalse(psc._is_aggregate_af_only(_Reader({"aggregation": "carriers"})))
        self.assertFalse(psc._is_aggregate_af_only(_Reader({"aggregation": "panel"})))

    def test_no_manifest_is_not_aggregate(self):
        self.assertFalse(psc._is_aggregate_af_only(object()))
        self.assertFalse(psc._is_aggregate_af_only(None))


class TestGroupBreakdownNA(unittest.TestCase):
    def _summary(self):
        return _Summary({
            "global": _GS(0.20, True, 0.30, 100),
            "1000G2021": _GS(0.17, True, 0.0, 0),   # aggregate: carrier freq/n are artifacts
            "gnomAD": _GS(0.22, True, 0.0, 0),
        })

    def test_aggregate_renders_carrier_na_af_intact(self):
        af_by, cf_by, n_by = psc._encode_group_breakdowns(
            self._summary(), "global", "::", aggregate_af_only=True)
        # allele-freq breakdown is REAL (unchanged)
        self.assertIn("1000G2021=0.17", af_by)
        self.assertIn("gnomAD=0.22", af_by)
        # carrier freq + count are NA per group (never a fabricated 0)
        self.assertEqual(cf_by, "1000G2021=NA;gnomAD=NA")
        self.assertEqual(n_by, "1000G2021=NA;gnomAD=NA")

    def test_genotyped_renders_real_carrier(self):
        s = _Summary({"global": _GS(0.2, True, 0.3, 100),
                      "EUR": _GS(0.1, True, 0.15, 42)})
        af_by, cf_by, n_by = psc._encode_group_breakdowns(
            s, "global", "::", aggregate_af_only=False)
        self.assertEqual(cf_by, "EUR=0.15")
        self.assertEqual(n_by, "EUR=42")


if __name__ == "__main__":
    unittest.main(verbosity=2)

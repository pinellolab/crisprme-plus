#!/usr/bin/env python3
"""Tests for compile_registry_from_info_af -- the sites-only "mega" INFO-AF registry.

The mega all-source index is a merged, genotype-stripped VCF: each site carries a
per-dataset allele frequency (INFO/AF_<dataset>) and AF_max. This builder turns
those AFs directly into the mmap registry (no genotypes to count). Guarantees:

  (1) AF round-trips: each per-dataset db group's AF == the input AF within
      0.5/AN_nominal (exact for practical purposes at real cohort sizes).
  (2) The GLOBAL group's AF == AF_max (the max per-dataset AF at the site),
      within 0.5/pooled_AN.
  (3) Every returned group is a well-formed Counts: n_hom <= n_carrier <= n_called
      and AC <= AN.
  (4) Sparsity + floor: a dataset absent at a site yields no group; a positive AF
      that would round AC below 1 is floored to a present AC=1 (a MAF-passing site
      is never silently dropped).
  (5) It round-trips through the SHIPPED RegistryReader -- the binary format is
      byte-for-byte the compile_registry format (this reuses _write_registry).
  (6) _af_only_counts unit properties (AF exactness; carrier/hom are 0, never an
      estimate).
  (7) No shipped page says those counts are Hardy-Weinberg estimates.

STDLIB only.
"""
import os
import re
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from tier0_registry import (  # noqa: E402
    compile_registry_from_info_af,
    _af_only_counts,
    RegistryReader,
    Counts,
    GLOBAL_GROUP_ID,
    db_group_id,
)

# The five mega datasets with their nominal allele numbers (2 x documented N).
DATASET_META = {
    "1000G2021": {"sample_count": 3202, "an_nominal": 6404},
    "HGDP": {"sample_count": 929, "an_nominal": 1858},
    "gnomAD": {"sample_count": 807162, "an_nominal": 1614324},
    "TOPMed": {"sample_count": 53831, "an_nominal": 107662},
    "AoU": {"sample_count": 535662, "an_nominal": 1071324},
}
_AN = {ds: m["an_nominal"] for ds, m in DATASET_META.items()}


def _well_formed(cnt):
    return (0 <= cnt.n_hom_indiv <= cnt.n_carrier_indiv <= cnt.n_called_indiv
            and 0 <= cnt.AC <= cnt.AN)


class TestAfOnlyCountsUnit(unittest.TestCase):
    def test_af_exact_within_half_allele(self):
        # In-range af (>= the MAF floor 0.001, so >= 1/(2*AN) for every real AN) is
        # reproduced within 0.5/AN -- the floor-to-1 branch is NOT exercised here.
        for an in (6404, 1858, 1614324, 107662, 1071324):
            for af in (0.001, 0.0137, 0.05, 0.2764, 0.5, 0.91, 0.999):
                c = _af_only_counts(af, an, an // 2)
                self.assertLessEqual(abs(c.allele_freq() - af), 0.5 / an + 1e-12,
                                     f"an={an} af={af}: got {c.allele_freq()}")

    def test_carrier_hom_are_zero_af_only(self):
        # AF-only: per-individual counts are NOT modelled (never HWE) -> always 0.
        for af in (0.001, 0.05, 0.5, 0.999):
            c = _af_only_counts(af, 6404, 3202)
            self.assertEqual((c.n_carrier_indiv, c.n_hom_indiv), (0, 0))
            self.assertEqual(c.n_called_indiv, 3202)
            self.assertTrue(_well_formed(c))

    def test_boundary_af(self):
        c0 = _af_only_counts(0.0, 6404, 3202)
        self.assertEqual(c0.as_tuple(), (0, 6404, 0, 0, 3202))
        c1 = _af_only_counts(1.0, 6404, 3202)
        self.assertEqual((c1.AC, c1.AN), (6404, 6404))
        self.assertEqual((c1.n_carrier_indiv, c1.n_hom_indiv), (0, 0))

    def test_tiny_af_is_floored_present(self):
        # af that rounds AC below 1 (only reachable OUTSIDE the recipe) -> AC=1 so
        # the site is never dropped; carrier/hom stay 0 (AF-only).
        c = _af_only_counts(1e-7, 6404, 3202)
        self.assertEqual(c.AC, 1)
        self.assertEqual((c.n_carrier_indiv, c.n_hom_indiv), (0, 0))
        self.assertTrue(_well_formed(c))

    def test_af_gt_one_is_capped(self):
        # a caller-supplied garbage af>1 must not make AC exceed AN.
        c = _af_only_counts(1.5, 6404, 3202)
        self.assertLessEqual(c.AC, c.AN)
        self.assertTrue(_well_formed(c))

    def test_n_called_zero_stays_bounded(self):
        # a misconfigured dataset (sample_count 0) must not yield n_carrier>n_called.
        c = _af_only_counts(0.05, 0, 0)
        self.assertTrue(_well_formed(c))
        self.assertEqual((c.n_carrier_indiv, c.n_hom_indiv, c.n_called_indiv),
                         (0, 0, 0))


class TestCompileFromInfoAf(unittest.TestCase):
    def _build(self, records):
        d = tempfile.mkdtemp()
        binp = os.path.join(d, "reg.bin")
        idxp = os.path.join(d, "reg.idx")
        manifest = compile_registry_from_info_af(records, DATASET_META, binp, idxp)
        return RegistryReader(binp, idxp), manifest

    def test_roundtrip_af_and_global_is_afmax(self):
        # One site in all five datasets with distinct AFs.
        afs = {"1000G2021": 0.0377889, "HGDP": 0.030541, "gnomAD": 0.00139027,
               "TOPMed": 0.034036, "AoU": 0.00322}
        reader, manifest = self._build([(10516173, "A", "G", "rs1", dict(afs))])
        groups = reader.lookup(10516173, "G")
        # every present dataset -> one db group + the GLOBAL group, nothing else.
        self.assertEqual(set(groups),
                         {db_group_id(ds) for ds in afs} | {GLOBAL_GROUP_ID})
        for ds, af in afs.items():
            g = groups[db_group_id(ds)]
            self.assertTrue(_well_formed(g))
            self.assertLessEqual(abs(g.allele_freq() - af), 0.5 / _AN[ds] + 1e-12)
        # GLOBAL AF == AF_max over the pooled nominal denominator.
        af_max = max(afs.values())
        pooled_an = sum(_AN.values())
        gg = groups[GLOBAL_GROUP_ID]
        self.assertLessEqual(abs(gg.allele_freq() - af_max), 0.5 / pooled_an + 1e-12)
        self.assertEqual(manifest["aggregation"], "info_af")
        self.assertTrue(manifest["databases"]["gnomAD"]["aggregate_af_only"])

    def test_partial_dataset_presence_is_sparse(self):
        # Only two datasets carry this site -> only their groups (+ GLOBAL) exist.
        reader, _ = self._build(
            [(200, "C", "T", "rs2", {"HGDP": 0.276596, "TOPMed": 0.077978})])
        groups = reader.lookup(200, "T")
        self.assertEqual(set(groups),
                         {db_group_id("HGDP"), db_group_id("TOPMed"),
                          GLOBAL_GROUP_ID})
        # GLOBAL AF == the larger of the two (HGDP), pooled over the two ANs only.
        pooled = _AN["HGDP"] + _AN["TOPMed"]
        self.assertLessEqual(
            abs(groups[GLOBAL_GROUP_ID].allele_freq() - 0.276596),
            0.5 / pooled + 1e-12)

    def test_multiallelic_distinct_alts(self):
        recs = [
            (500, "C", "T", "rs3", {"1000G2021": 0.10, "gnomAD": 0.02}),
            (500, "C", "A", "rs4", {"AoU": 0.004}),
        ]
        reader, _ = self._build(recs)
        gt = reader.lookup(500, "T")
        ga = reader.lookup(500, "A")
        self.assertEqual(set(gt),
                         {db_group_id("1000G2021"), db_group_id("gnomAD"),
                          GLOBAL_GROUP_ID})
        self.assertEqual(set(ga), {db_group_id("AoU"), GLOBAL_GROUP_ID})
        self.assertLessEqual(
            abs(gt[db_group_id("1000G2021")].allele_freq() - 0.10),
            0.5 / _AN["1000G2021"] + 1e-12)

    def test_tiny_af_site_survives(self):
        # A rare-but-MAF-passing AF that rounds AC<1 at small AN is still present.
        reader, _ = self._build([(999, "G", "A", "rs5", {"HGDP": 5e-5})])
        groups = reader.lookup(999, "A")
        self.assertIn(db_group_id("HGDP"), groups)
        self.assertEqual(groups[db_group_id("HGDP")].AC, 1)

    def test_all_groups_well_formed_over_many_records(self):
        import random
        random.seed(7)
        recs = []
        for i in range(500):
            present = random.sample(list(DATASET_META), random.randint(1, 5))
            afs = {ds: round(random.uniform(0.001, 0.999), 6) for ds in present}
            recs.append((1000 + i, "A", "C", f"rs{i}", afs))
        reader, _ = self._build(recs)
        for i in range(500):
            for gid, cnt in reader.lookup(1000 + i, "C").items():
                self.assertTrue(_well_formed(cnt), f"rec {i} group {gid}: {cnt}")


class TheDocsSayWhatTheRegistryHolds(unittest.TestCase):
    """METHODS section 2 said that a sites-only panel's carrier and homozygote counts "are
    reported as Hardy-Weinberg expectations from the allele frequency (the frequency itself is
    exact; carrier/hom are flagged as estimates)", and docs/DESIGN_mega_index.md and
    seq_script/merge_panels/README.md that they are HWE expectations, from a helper
    (_hwe_counts_from_af) that tier0_registry does not have. compile_registry_from_info_af
    estimates none: it stores 0 for both (test_carrier_hom_are_zero_af_only above), and the
    population-summary companion gives them as NA (test_aggregate_carrier_na)."""

    CLAIM = re.compile(
        r"(?:Hardy[-\u2013]Weinberg|HWE) (?:expectations|estimates)"
        r"|carrier/hom(?:ozygote)? (?:counts )?are (?:flagged as )?(?:HWE )?"
        r"(?:expectations|estimates)"
        r"|_hwe_counts_from_af|HWE bounds", re.I)

    @staticmethod
    def _prose(text):
        """Markdown emphasis dropped and whitespace runs folded, so a claim split over two
        lines or set in bold reads as one."""
        return " ".join(text.replace("**", "").split())

    def test_the_scan_finds_what_the_pages_said(self):
        for said in ("so they are reported as\nHardy\u2013Weinberg expectations from the "
                     "allele frequency (the frequency itself is exact;\ncarrier/hom are "
                     "flagged as estimates)",
                     "counts do not exist for aggregate data \u2192 **Hardy-Weinberg "
                     "expectations** from AF",
                     "bounded + self-consistent; **AF is\nexact, carrier/hom are estimates**",
                     "- **Carrier/hom counts are HWE expectations** from AF",
                     "helper\n`_hwe_counts_from_af`) reuses `_write_registry`"):
            self.assertTrue(self.CLAIM.search(self._prose(said)), said)
        for not_said in ("the registry (never a fabricated Hardy-Weinberg estimate)",
                         "view -i'HWE>=1e-30 && F_MISSING<=0.1 && ExcHet>=0.5",
                         "**Carrier/hom counts are not estimated**: the registry stores 0"):
            self.assertIsNone(self.CLAIM.search(self._prose(not_said)), not_said)

    def test_no_shipped_page_says_they_are_estimated(self):
        import tier0_registry
        self.assertFalse(hasattr(tier0_registry, "_hwe_counts_from_af"))
        repo = os.path.dirname(_HERE)
        pages = [os.path.join(repo, "README.md"), os.path.join(repo, "METHODS.md")]
        for top in ("docs", "seq_script"):
            for root, _dirs, files in os.walk(os.path.join(repo, top)):
                pages += [os.path.join(root, f) for f in sorted(files) if f.endswith(".md")]
        self.assertIn(os.path.join(repo, "docs", "DESIGN_mega_index.md"), pages)
        self.assertIn(os.path.join(repo, "seq_script", "merge_panels", "README.md"), pages)
        found = []
        for page in pages:
            with open(page, encoding="utf-8") as fh:
                for m in self.CLAIM.finditer(self._prose(fh.read())):
                    found.append("%s: %s" % (os.path.relpath(page, repo), m.group(0)))
        self.assertEqual(found, [])
        # and METHODS section 2 says what it holds
        with open(os.path.join(repo, "METHODS.md"), encoding="utf-8") as fh:
            methods = self._prose(fh.read())
        self.assertIn("Per-individual carrier and homozygote counts do not exist for aggregate "
                      "data, and none is estimated: the registry stores 0 for both", methods)


if __name__ == "__main__":
    unittest.main(verbosity=2)

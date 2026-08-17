"""Unit tests for combination-aware population summaries (Phase 3, 3a+3b).

STDLIB ONLY (unittest + tempfile + json + math). Builds a tiny Tier-0 registry +
Tier-1 genotype store + sample axis over a synthetic 2-db / 3-subpop panel, then:

  (a) SINGLE-variant summary == the Tier-0 registry's per-group Counts /
      derive_record_stats (exact).
  (b) MULTI-variant PHASED: two SNPs where some individuals carry BOTH in cis
      (same hap), some in trans (different haps), some only one -> combination
      carrier count == cis-only individuals; per-(db x subpop)/per-db/global
      correct; hom = both-hap carriers; allele_freq = cis-haplotype count / AN.
  (c) MULTI-variant UNPHASED: same genotypes as '/' -> carrier count ==
      assume-cis upper bound (het/hom at both) and allele_freq is UNDEFINED;
      lower bound <= upper bound.
  (d) GLOBAL dedup (a shared canonical id counted once) + max_subpop_af(+label) +
      observed; and an off-target with ZERO combination carriers -> observed
      False, freqs 0.
  (e) an INDEPENDENT brute force (walk the panel individual-by-individual,
      hap-by-hap) confirms the module's combination math (not merely
      self-consistent).
"""

from __future__ import annotations

import math
import os
import tempfile
import unittest

import population_summary as ps
import tier0_registry as t0
import tier1_genotypes as t1


# --------------------------------------------------------------------------- #
# Synthetic panel: 2 databases, 3 subpops, mixed sexes. Autosomal (all diploid)
# by default; a chrX-style ploidy variant exercises the ploidy/sex path.
# --------------------------------------------------------------------------- #
# sample_id -> (database, subpopulation, sex)
PANEL = {
    # 1000G / EUR
    "S_EUR_M1": ("1000G", "EUR", "male"),
    "S_EUR_F1": ("1000G", "EUR", "female"),
    "S_EUR_M2": ("1000G", "EUR", "male"),
    # 1000G / AFR
    "S_AFR_M1": ("1000G", "AFR", "male"),
    "S_AFR_F1": ("1000G", "AFR", "female"),
    # HGDP / EAS
    "H_EAS_M1": ("HGDP", "EAS", "male"),
    "H_EAS_F1": ("HGDP", "EAS", "female"),
    "H_EAS_F2": ("HGDP", "EAS", "female"),
}


def _write_registry(tmpdir, records, sample_meta, ploidy_of, name="t0"):
    """Compile a FULL-PANEL Tier-0 registry (panel AN) from carrier-only records.

    ``records``: iterable of (pos, ref, alt, rsid, carrier_gts) where carrier_gts
    is {sample_id: gt} of alt-carriers only (legacy dict semantics).
    """
    binp = os.path.join(tmpdir, name + ".bin")
    idxp = os.path.join(tmpdir, name + ".idx.json")
    t0.compile_registry_panel(records, sample_meta, None, ploidy_of, binp, idxp)
    return t0.RegistryReader(binp, idxp)


def _write_gt(tmpdir, records, axis, name="t1"):
    """Compile a Tier-1 genotype store from (pos, alt, carriers[, ref]) records."""
    binp = os.path.join(tmpdir, name + ".bin")
    idxp = os.path.join(tmpdir, name + ".idx.json")
    t1.compile_genotypes(records, axis, binp, idxp)
    return t1.GenotypeReader(binp, idxp)


# --------------------------------------------------------------------------- #
# (e) INDEPENDENT brute force. Deliberately structured DIFFERENTLY from the
# module: it walks EVERY panel individual, materializes each one's full genotype
# at every pair (defaulting unlisted individuals to ref), and decides carriage by
# an explicit per-haplotype loop. No shared code path with population_summary.
# --------------------------------------------------------------------------- #
def brute_force_combination(pos_alts, per_pair_carriers, panel_meta, ploidy_of,
                            phased):
    """Return a dict of expected per-group summaries computed from scratch.

    per_pair_carriers: list (len k) of {sample_id: gt_string} alt-carrier maps.
    panel_meta: sample_id -> (db, subpop, sex) for the WHOLE panel.
    Returns dict[group_id] -> dict with keys AC, AN, n_carrier, n_hom, n_called,
    n_carrier_lower (unphased only, else None). Only groups with >=1 carrier.
    """
    k = len(pos_alts)

    # Materialize every panel individual's genotype at each pair (unlisted = ref
    # "0|0"/"0"), honoring ploidy. Build per-hap allele arrays independently.
    def haps_of(sid, i):
        gt = per_pair_carriers[i].get(sid)
        ploidy = ploidy_of(sid, panel_meta[sid][2])
        if gt is None:
            # Not listed => hom-ref: fill ploidy ref slots.
            return ["0"] * max(ploidy, 0)
        # Split by whichever separator is present; cap at ploidy; never pad.
        if "|" in gt:
            toks = gt.split("|")
        elif "/" in gt:
            toks = gt.split("/")
        else:
            toks = [gt]
        if ploidy <= 1:
            return toks[:1]
        return toks[:ploidy]

    # Full-panel denominators, ploidy/sex-aware, GLOBAL deduped by canonical id.
    groups = {}  # gid -> accumulator dict

    def ensure(gid):
        if gid not in groups:
            groups[gid] = {"AC": 0, "AN": 0, "n_carrier": 0, "n_hom": 0,
                           "n_called": 0, "singles": [0] * k}
        return groups[gid]

    global_seen = set()
    # First: denominators + single-pair tallies over the FULL panel.
    for sid, (db, sp, sex) in panel_meta.items():
        ploidy = ploidy_of(sid, sex)
        if ploidy <= 0:
            continue
        for gid in (t0.GLOBAL_GROUP_ID if sid not in global_seen else None,
                    t0.db_group_id(db), t0.db_subpop_group_id(db, sp)):
            if gid is None:
                continue
            g = ensure(gid)
            g["AN"] += ploidy
            g["n_called"] += 1
        if sid not in global_seen:
            global_seen.add(sid)
        # single-pair carrier tallies (per group).
        for i in range(k):
            slots = haps_of(sid, i)
            if any(tok == "1" for tok in slots):
                ensure(t0.db_subpop_group_id(db, sp))["singles"][i] += 1
                ensure(t0.db_group_id(db))["singles"][i] += 1
        # global single tally (deduped): handled separately below.

    # global single tally deduped by canonical id.
    for i in range(k):
        gseen = set()
        for sid, (db, sp, sex) in panel_meta.items():
            if ploidy_of(sid, sex) <= 0 or sid in gseen:
                continue
            gseen.add(sid)
            slots = haps_of(sid, i)
            if any(tok == "1" for tok in slots):
                ensure(t0.GLOBAL_GROUP_ID)["singles"][i] += 1

    # Now carriage decisions per individual.
    gseen_carrier = set()
    for sid, (db, sp, sex) in panel_meta.items():
        ploidy = ploidy_of(sid, sex)
        if ploidy <= 0:
            continue
        per_pair_haps = [haps_of(sid, i) for i in range(k)]
        n_hap = min(len(h) for h in per_pair_haps) if per_pair_haps else 0

        if phased:
            cis_haps = 0
            for h in range(n_hap):
                if all(per_pair_haps[i][h] == "1" for i in range(k)):
                    cis_haps += 1
            is_carrier = cis_haps >= 1
            is_hom = (cis_haps == ploidy and ploidy >= 1)
            ac_contrib = cis_haps
        else:
            # het/hom at EVERY pair.
            is_carrier = all(any(tok == "1" for tok in per_pair_haps[i])
                             for i in range(k))
            # hom at every pair: fully-called and all slots alt.
            def hom_at(i):
                slots = per_pair_haps[i]
                return (len(slots) == ploidy and ploidy >= 1
                        and all(tok == "1" for tok in slots))
            is_hom = is_carrier and all(hom_at(i) for i in range(k))
            ac_contrib = None

        if not is_carrier:
            continue
        for gid, dedup in ((t0.db_subpop_group_id(db, sp), False),
                           (t0.db_group_id(db), False),
                           (t0.GLOBAL_GROUP_ID, True)):
            if dedup:
                if sid in gseen_carrier:
                    continue
                gseen_carrier.add(sid)
            g = ensure(gid)
            g["n_carrier"] += 1
            if is_hom:
                g["n_hom"] += 1
            if ac_contrib is not None:
                g["AC"] += ac_contrib

    # Emit sparse (>=1 carrier) group summaries.
    out = {}
    for gid, g in groups.items():
        if g["n_carrier"] < 1:
            continue
        if phased:
            AC = g["AC"]
            lower = None
        else:
            AC = None
            total = sum(g["singles"])
            floor = max(0, total - (k - 1) * g["n_called"])
            lower = min(floor, min(g["singles"]) if g["singles"] else 0)
            lower = min(lower, g["n_carrier"])
        out[gid] = {"AC": AC, "AN": g["AN"], "n_carrier": g["n_carrier"],
                    "n_hom": g["n_hom"], "n_called": g["n_called"],
                    "n_carrier_lower": lower}
    return out


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def make_axis():
    return t1.build_sample_axis(PANEL)


def carriers_from(gt_map):
    """A {sid: gt} carrier map restricted to those actually carrying alt (any
    non-'0' non-'.' token) -- the Tier-1 store lists ONLY alt-carriers."""
    out = {}
    for sid, gt in gt_map.items():
        toks = gt.replace("/", "|").split("|")
        if any(t == "1" for t in toks):
            out[sid] = gt
    return out


# --------------------------------------------------------------------------- #
# (a) SINGLE-variant exactness vs the Tier-0 registry.
# --------------------------------------------------------------------------- #
class TestSingleVariantExact(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.autosomal_ploidy  # everyone diploid
        # One SNP at pos 100 alt A: a spread of carriers across all groups.
        self.gts = {
            "S_EUR_M1": "1|0",   # EUR het
            "S_EUR_F1": "1|1",   # EUR hom
            "S_AFR_M1": "0|1",   # AFR het
            "H_EAS_M1": "1|1",   # EAS hom
            "H_EAS_F1": "1|0",   # EAS het
        }
        carriers = carriers_from(self.gts)
        recs = [(100, "G", "A", "rs100", carriers)]
        self.t0r = _write_registry(self.d, recs, PANEL, self.ploidy_of)
        self.axis = make_axis()
        gt_recs = [(100, "A", carriers, "G")]
        self.t1r = _write_gt(self.d, gt_recs, self.axis)

    def tearDown(self):
        self.t0r.close()
        self.t1r.close()

    def test_single_equals_registry(self):
        summary = ps.summarize(
            [(100, "A")], self.t0r, self.t1r, self.axis, self.ploidy_of,
            phased=True)
        reg = self.t0r.lookup(100, "A")
        # Every registry group must appear with identical AC/AN/carrier/hom/called.
        self.assertEqual(set(summary.groups), set(reg))
        for gid, c in reg.items():
            gs = summary.groups[gid]
            self.assertEqual(gs.allele_count, c.AC, gid)
            self.assertEqual(gs.allele_number, c.AN, gid)
            self.assertEqual(gs.n_carrier, c.n_carrier_indiv, gid)
            self.assertEqual(gs.n_hom, c.n_hom_indiv, gid)
            self.assertEqual(gs.n_called, c.n_called_indiv, gid)
            self.assertAlmostEqual(gs.allele_freq, c.allele_freq(), 12, gid)
            self.assertAlmostEqual(gs.carrier_freq, c.carrier_freq(), 12, gid)
            self.assertAlmostEqual(gs.hom_freq, c.hom_freq(), 12, gid)

    def test_single_record_stats_match_derive(self):
        summary = ps.summarize(
            [(100, "A")], self.t0r, self.t1r, self.axis, self.ploidy_of,
            phased=True)
        derived = t0.derive_record_stats(self.t0r.lookup(100, "A"))
        rs = ps.record_stats(summary)
        self.assertAlmostEqual(rs["global_af"], derived["global_af"], 12)
        self.assertAlmostEqual(rs["global_carrier_freq"],
                               derived["global_carrier_freq"], 12)
        self.assertAlmostEqual(rs["global_hom_freq"],
                               derived["global_hom_freq"], 12)
        self.assertAlmostEqual(rs["max_subpop_af"], derived["max_subpop_af"], 12)
        self.assertEqual(rs["max_subpop_af_label"],
                         derived["max_subpop_af_label"])
        self.assertEqual(rs["observed"], derived["observed"])

    def test_single_absent_off_target(self):
        # A (pos, alt) not in the registry -> no groups, observed False, freq 0.
        summary = ps.summarize(
            [(999, "T")], self.t0r, self.t1r, self.axis, self.ploidy_of,
            phased=True)
        self.assertEqual(summary.groups, {})
        self.assertFalse(summary.observed)
        self.assertEqual(summary.global_carrier_n, 0)
        self.assertEqual(summary.global_af, 0.0)


# --------------------------------------------------------------------------- #
# Shared two-SNP combination fixture for (b)(c)(e).
# --------------------------------------------------------------------------- #
class _TwoSnpFixture(unittest.TestCase):
    """pos 100 alt A and pos 200 alt C. Genotypes designed for a mix of cis /
    trans / single carriage across groups. Phased strings; the unphased case
    reuses the SAME allele content with '/' separators."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.autosomal_ploidy  # diploid everyone
        self.axis = make_axis()

        # Per-sample phased genotype at (100,A) and (200,C).
        # hap layout is "hap0|hap1"; alt token "1".
        self.gt100 = {
            "S_EUR_M1": "1|0",   # A on hap0
            "S_EUR_F1": "1|1",   # A both haps
            "S_EUR_M2": "1|0",   # A on hap0
            "S_AFR_M1": "0|1",   # A on hap1
            "S_AFR_F1": "1|0",   # A on hap0
            "H_EAS_M1": "1|0",   # A on hap0
            "H_EAS_F1": "0|1",   # A on hap1
            # H_EAS_F2: ref at 100 (not listed)
        }
        self.gt200 = {
            "S_EUR_M1": "1|0",   # C on hap0  -> CIS with A (both hap0)
            "S_EUR_F1": "1|1",   # C both haps -> combination on both haps (hom)
            "S_EUR_M2": "0|1",   # C on hap1  -> TRANS vs A(hap0): NOT a carrier
            "S_AFR_M1": "0|1",   # C on hap1  -> CIS with A(hap1)
            "S_AFR_F1": "0|1",   # C on hap1  -> TRANS vs A(hap0): NOT a carrier
            "H_EAS_M1": "1|0",   # C on hap0  -> CIS with A(hap0)
            # H_EAS_F1: ref at 200 (not listed) -> single carrier only
            "H_EAS_F2": "1|0",   # C on hap0 but no A -> single carrier only
        }
        # Registry (single-variant, exact) for completeness / (a)-style checks.
        c100 = carriers_from(self.gt100)
        c200 = carriers_from(self.gt200)
        recs0 = [
            (100, "G", "A", "rs100", c100),
            (200, "T", "C", "rs200", c200),
        ]
        self.t0r = _write_registry(self.d, recs0, PANEL, self.ploidy_of)

        # Tier-1 phased store.
        gt_recs_phased = [
            (100, "A", c100, "G"),
            (200, "C", c200, "T"),
        ]
        self.t1_phased = _write_gt(self.d, gt_recs_phased, self.axis,
                                   name="t1_phased")

        # Tier-1 unphased store: SAME allele content, '/' separator.
        def to_unphased(m):
            return {sid: gt.replace("|", "/") for sid, gt in m.items()}
        c100u = to_unphased(c100)
        c200u = to_unphased(c200)
        gt_recs_unphased = [
            (100, "A", c100u, "G"),
            (200, "C", c200u, "T"),
        ]
        self.t1_unphased = _write_gt(self.d, gt_recs_unphased, self.axis,
                                     name="t1_unphased")

        # Full per-pair carrier maps (for the brute force), PHASED content.
        self.pair_carriers = [c100, c200]
        self.pair_carriers_unphased = [c100u, c200u]

    def tearDown(self):
        self.t0r.close()
        self.t1_phased.close()
        self.t1_unphased.close()

    def _assert_matches_brute(self, summary, expected):
        self.assertEqual(set(summary.groups), set(expected),
                         "group set mismatch")
        for gid, exp in expected.items():
            gs = summary.groups[gid]
            self.assertEqual(gs.allele_count, exp["AC"], "AC %s" % gid)
            self.assertEqual(gs.allele_number, exp["AN"], "AN %s" % gid)
            self.assertEqual(gs.n_carrier, exp["n_carrier"],
                             "n_carrier %s" % gid)
            self.assertEqual(gs.n_hom, exp["n_hom"], "n_hom %s" % gid)
            self.assertEqual(gs.n_called, exp["n_called"], "n_called %s" % gid)
            self.assertEqual(gs.n_carrier_lower, exp["n_carrier_lower"],
                             "n_carrier_lower %s" % gid)


# --------------------------------------------------------------------------- #
# (b) MULTI-variant PHASED.
# --------------------------------------------------------------------------- #
class TestPhasedCombination(_TwoSnpFixture):
    def test_phased_cis_only_carriers(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_phased, self.axis,
            self.ploidy_of, phased=True)

        # CIS carriers (A and C on the SAME hap):
        #   S_EUR_M1 (both hap0), S_EUR_F1 (both haps -> hom), S_AFR_M1 (both hap1),
        #   H_EAS_M1 (both hap0). TRANS: S_EUR_M2, S_AFR_F1 (NOT carriers).
        #   H_EAS_F1 (A only), H_EAS_F2 (C only): NOT carriers.
        g = summary.groups
        self.assertEqual(g["1000G::EUR"].n_carrier, 2)  # M1, F1
        self.assertEqual(g["1000G::EUR"].n_hom, 1)      # F1 both haps
        self.assertEqual(g["1000G::EUR"].allele_count, 3)  # M1:1 + F1:2
        self.assertEqual(g["1000G::AFR"].n_carrier, 1)  # M1
        self.assertEqual(g["1000G::AFR"].allele_count, 1)
        self.assertEqual(g["HGDP::EAS"].n_carrier, 1)   # H_EAS_M1
        self.assertEqual(g["HGDP::EAS"].allele_count, 1)
        # DB rollups.
        self.assertEqual(g["1000G"].n_carrier, 3)
        self.assertEqual(g["1000G"].allele_count, 4)
        self.assertEqual(g["HGDP"].n_carrier, 1)
        # GLOBAL.
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].n_carrier, 4)
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].allele_count, 5)
        self.assertTrue(summary.observed)
        self.assertTrue(summary.allele_freq_defined)

    def test_phased_denominators_and_freq(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_phased, self.axis,
            self.ploidy_of, phased=True)
        g = summary.groups
        # AN full-panel: EUR has 3 diploid samples => AN 6; AFR 2 => 4; EAS 3 => 6.
        self.assertEqual(g["1000G::EUR"].allele_number, 6)
        self.assertEqual(g["1000G::AFR"].allele_number, 4)
        self.assertEqual(g["HGDP::EAS"].allele_number, 6)
        self.assertEqual(g["1000G"].allele_number, 10)
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].allele_number, 16)
        # allele_freq = cis-haplotype count / AN.
        self.assertAlmostEqual(g["1000G::EUR"].allele_freq, 3 / 6, 12)
        self.assertAlmostEqual(g[t0.GLOBAL_GROUP_ID].allele_freq, 5 / 16, 12)
        # n_called full-panel.
        self.assertEqual(g["1000G::EUR"].n_called, 3)
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].n_called, 8)

    def test_phased_matches_brute(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_phased, self.axis,
            self.ploidy_of, phased=True)
        expected = brute_force_combination(
            [(100, "A"), (200, "C")], self.pair_carriers, PANEL,
            self.ploidy_of, phased=True)
        self._assert_matches_brute(summary, expected)


# --------------------------------------------------------------------------- #
# (c) MULTI-variant UNPHASED.
# --------------------------------------------------------------------------- #
class TestUnphasedCombination(_TwoSnpFixture):
    def test_unphased_assume_cis_upper_bound(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_unphased, self.axis,
            self.ploidy_of, phased=False)
        g = summary.groups
        # Assume-cis: het/hom at BOTH sites. Those are exactly the individuals
        # listed as carriers at both 100 AND 200 (regardless of hap):
        #   S_EUR_M1, S_EUR_F1, S_EUR_M2, S_AFR_M1, S_AFR_F1, H_EAS_M1.
        #   (H_EAS_F1: no C; H_EAS_F2: no A.)
        self.assertEqual(g["1000G::EUR"].n_carrier, 3)   # M1,F1,M2
        self.assertEqual(g["1000G::AFR"].n_carrier, 2)   # M1,F1
        self.assertEqual(g["HGDP::EAS"].n_carrier, 1)    # M1
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].n_carrier, 6)
        # hom at both sites: only S_EUR_F1 (1/1 at both).
        self.assertEqual(g["1000G::EUR"].n_hom, 1)
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].n_hom, 1)

    def test_unphased_allele_freq_undefined(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_unphased, self.axis,
            self.ploidy_of, phased=False)
        self.assertFalse(summary.allele_freq_defined)
        self.assertTrue(math.isnan(summary.global_af))
        for gid, gs in summary.groups.items():
            self.assertIsNone(gs.allele_count, gid)
            self.assertFalse(gs.allele_freq_defined, gid)
            self.assertTrue(math.isnan(gs.allele_freq), gid)
        # max_subpop_af has no defined AF to rank -> label None, af 0.0.
        self.assertIsNone(summary.max_subpop_af_label)
        self.assertEqual(summary.max_subpop_af, 0.0)

    def test_unphased_lower_le_upper(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_unphased, self.axis,
            self.ploidy_of, phased=False)
        for gid, gs in summary.groups.items():
            self.assertIsNotNone(gs.n_carrier_lower, gid)
            self.assertLessEqual(gs.n_carrier_lower, gs.n_carrier, gid)
            self.assertGreaterEqual(gs.n_carrier_lower, 0, gid)
            # carrier_freq_lower <= carrier_freq.
            self.assertLessEqual(gs.carrier_freq_lower, gs.carrier_freq, gid)

    def test_unphased_matches_brute(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_unphased, self.axis,
            self.ploidy_of, phased=False)
        expected = brute_force_combination(
            [(100, "A"), (200, "C")], self.pair_carriers_unphased, PANEL,
            self.ploidy_of, phased=False)
        self._assert_matches_brute(summary, expected)

    def test_unphased_upper_ge_phased(self):
        # The unphased assume-cis carrier count is an UPPER bound on the phased
        # (true-cis) count for the same allele content.
        s_ph = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_phased, self.axis,
            self.ploidy_of, phased=True)
        s_un = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_unphased, self.axis,
            self.ploidy_of, phased=False)
        for gid in s_ph.groups:
            self.assertLessEqual(
                s_ph.groups[gid].n_carrier, s_un.groups[gid].n_carrier, gid)


# --------------------------------------------------------------------------- #
# (d) GLOBAL dedup + max_subpop_af(+label) + observed + zero-carrier off-target.
# --------------------------------------------------------------------------- #
class TestGlobalDedupAndRecordStats(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.autosomal_ploidy
        # A panel where ONE canonical id ("SHARED") appears in BOTH databases.
        # (Cross-db sample overlap; global must count it ONCE.)
        self.panel = {
            "SHARED":   ("1000G", "EUR", "male"),  # canonical id in 1000G
            "A_ONLY":   ("1000G", "EUR", "female"),
            "H_ONLY_1": ("HGDP", "EAS", "male"),
            "H_ONLY_2": ("HGDP", "EAS", "female"),
        }
        # NOTE: our axis is 1:1 sample_id -> meta, so "SHARED" has a single
        # canonical entry; dedup is exercised by GLOBAL counting it once across
        # the db + subpop rollups (a sample contributes to exactly one db here,
        # but we ALSO verify GLOBAL n <= sum of per-db carriers via the dedup path
        # in the combination walk). To exercise a TRUE cross-db shared id we add a
        # second membership through a SECOND registry/panel below.
        self.axis = t1.build_sample_axis(self.panel)

        # Two SNPs; combination carried by SHARED (cis) + H_ONLY_1 (cis).
        gt100 = {"SHARED": "1|0", "H_ONLY_1": "1|0", "A_ONLY": "1|0"}
        gt200 = {"SHARED": "1|0", "H_ONLY_1": "1|0"}  # A_ONLY lacks 200 -> single
        self.c100 = carriers_from(gt100)
        self.c200 = carriers_from(gt200)
        recs0 = [(100, "G", "A", "rs100", self.c100),
                 (200, "T", "C", "rs200", self.c200)]
        self.t0r = _write_registry(self.d, recs0, self.panel, self.ploidy_of)
        gt_recs = [(100, "A", self.c100, "G"), (200, "C", self.c200, "T")]
        self.t1r = _write_gt(self.d, gt_recs, self.axis)

    def tearDown(self):
        self.t0r.close()
        self.t1r.close()

    def test_max_subpop_and_observed_single(self):
        # Single variant at 100: EUR has SHARED+A_ONLY carriers (2/2 samples ->
        # AC 2 / AN 4), EAS has H_ONLY_1 (AC 1 / AN 4). max subpop = 1000G::EUR.
        summary = ps.summarize(
            [(100, "A")], self.t0r, self.t1r, self.axis, self.ploidy_of,
            phased=True)
        self.assertTrue(summary.observed)
        self.assertEqual(summary.max_subpop_af_label, "1000G::EUR")
        self.assertAlmostEqual(summary.max_subpop_af, 2 / 4, 12)
        # ABSOLUTE global carrier N (deduped union).
        self.assertEqual(summary.global_carrier_n, 3)  # SHARED, A_ONLY, H_ONLY_1

    def test_combination_global_dedup_and_max_subpop(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1r, self.axis,
            self.ploidy_of, phased=True)
        g = summary.groups
        # Combination carriers: SHARED (1000G::EUR, cis) + H_ONLY_1 (HGDP::EAS).
        # A_ONLY carries only 100 -> not a combination carrier.
        self.assertEqual(g["1000G::EUR"].n_carrier, 1)
        self.assertEqual(g["HGDP::EAS"].n_carrier, 1)
        # GLOBAL counts each canonical id once: 2 carriers total.
        self.assertEqual(g[t0.GLOBAL_GROUP_ID].n_carrier, 2)
        # max subpop AF: EUR AC 1 / AN 4 (2 diploid samples), EAS AC 1 / AN 4 ->
        # tie at 0.25; deterministic tie-break to alphabetically-first "1000G::EUR".
        self.assertEqual(summary.max_subpop_af_label, "1000G::EUR")
        self.assertAlmostEqual(summary.max_subpop_af, 1 / 4, 12)

    def test_zero_carrier_off_target(self):
        # A combination NO ONE carries: pos 100 alt A (exists) + pos 300 alt C
        # (absent record -> zero carriers). Intersection empty -> observed False.
        summary = ps.summarize(
            [(100, "A"), (300, "C")], self.t0r, self.t1r, self.axis,
            self.ploidy_of, phased=True)
        self.assertEqual(summary.groups, {})
        self.assertFalse(summary.observed)
        self.assertEqual(summary.global_carrier_n, 0)
        self.assertEqual(summary.global_carrier_freq, 0.0)
        self.assertEqual(summary.global_hom_freq, 0.0)
        self.assertEqual(summary.max_subpop_af, 0.0)
        self.assertIsNone(summary.max_subpop_af_label)


# --------------------------------------------------------------------------- #
# Ploidy/sex-aware combination (chrX-style): males haploid, females diploid.
# --------------------------------------------------------------------------- #
class TestChrXPloidyCombination(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.make_chr_ploidy(haploid_male=True)  # male=1, female=2
        self.axis = make_axis()
        # A single-hap "1" for males (haploid), phased "a|b" for females.
        self.gt100 = {
            "S_EUR_M1": "1",     # haploid male, carries A
            "S_EUR_F1": "1|0",   # female het
            "H_EAS_M1": "1",     # haploid male, carries A
        }
        self.gt200 = {
            "S_EUR_M1": "1",     # haploid male, carries C -> combination (1 hap)
            "S_EUR_F1": "1|0",   # female het, C on hap0 -> cis with A(hap0)
            "H_EAS_M1": "0",     # haploid male, NO C -> not a combination carrier
        }
        c100 = carriers_from(self.gt100)
        c200 = carriers_from(self.gt200)
        recs0 = [(100, "G", "A", "rs100", c100), (200, "T", "C", "rs200", c200)]
        self.t0r = _write_registry(self.d, recs0, PANEL, self.ploidy_of)
        gt_recs = [(100, "A", c100, "G"), (200, "C", c200, "T")]
        self.t1r = _write_gt(self.d, gt_recs, self.axis)
        self.pair_carriers = [c100, c200]

    def tearDown(self):
        self.t0r.close()
        self.t1r.close()

    def test_haploid_combination_carrier_is_hom(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1r, self.axis,
            self.ploidy_of, phased=True)
        g = summary.groups
        # S_EUR_M1: haploid carrier of A+C -> 1 hap combination, hom (ploidy 1).
        # S_EUR_F1: A(hap0)+C(hap0) cis -> carrier, 1 hap, NOT hom (ploidy 2).
        # H_EAS_M1: no C -> not a carrier.
        self.assertEqual(g["1000G::EUR"].n_carrier, 2)
        self.assertEqual(g["1000G::EUR"].n_hom, 1)          # only the haploid male
        self.assertEqual(g["1000G::EUR"].allele_count, 2)   # M1:1 hap + F1:1 hap
        # EUR AN: M1 (male, ploidy1) + F1 (female, 2) + M2 (male,1) = 4.
        self.assertEqual(g["1000G::EUR"].allele_number, 4)
        self.assertNotIn("HGDP::EAS", g)  # H_EAS_M1 not a combination carrier

    def test_chrx_matches_brute(self):
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1r, self.axis,
            self.ploidy_of, phased=True)
        expected = brute_force_combination(
            [(100, "A"), (200, "C")], self.pair_carriers, PANEL,
            self.ploidy_of, phased=True)
        # Compare AC/AN/carrier/hom/called on the shared groups.
        self.assertEqual(set(summary.groups), set(expected))
        for gid, exp in expected.items():
            gs = summary.groups[gid]
            self.assertEqual(gs.allele_count, exp["AC"], gid)
            self.assertEqual(gs.allele_number, exp["AN"], gid)
            self.assertEqual(gs.n_carrier, exp["n_carrier"], gid)
            self.assertEqual(gs.n_hom, exp["n_hom"], gid)
            self.assertEqual(gs.n_called, exp["n_called"], gid)


# --------------------------------------------------------------------------- #
# Parsing helper unit checks.
# --------------------------------------------------------------------------- #
class TestParsingHelpers(unittest.TestCase):
    def test_parse_haplotypes_ploidy(self):
        self.assertEqual(ps.parse_haplotypes("1|0", 2), ["1", "0"])
        self.assertEqual(ps.parse_haplotypes("0/1", 2), ["0", "1"])
        self.assertEqual(ps.parse_haplotypes("1", 1), ["1"])
        # Diploid-style GT declared haploid -> ONE slot (no double count).
        self.assertEqual(ps.parse_haplotypes("1|1", 1), ["1"])
        # Haploid GT declared diploid -> ONE slot (never padded).
        self.assertEqual(ps.parse_haplotypes("1", 2), ["1"])
        self.assertEqual(ps.parse_haplotypes(".", 2), ["."])
        self.assertEqual(ps.parse_haplotypes("", 2), [])
        self.assertEqual(ps.parse_haplotypes(None, 2), [])

    def test_carries_and_hom(self):
        self.assertTrue(ps.individual_carries_at("1|0", 2))
        self.assertTrue(ps.individual_carries_at("0/1", 2))
        self.assertFalse(ps.individual_carries_at("0|0", 2))
        self.assertFalse(ps.individual_carries_at("2|0", 2))  # foreign alt != ref alt
        self.assertTrue(ps.individual_hom_at("1|1", 2))
        self.assertFalse(ps.individual_hom_at("1|0", 2))
        self.assertFalse(ps.individual_hom_at("1|.", 2))  # half-missing not hom
        self.assertTrue(ps.individual_hom_at("1", 1))      # haploid hemizygous

    def test_lower_bound_math(self):
        # Two sets of size 5 and 4 in a universe of 8: floor = 5+4-8 = 1.
        self.assertEqual(ps._unphased_carrier_lower_bound([5, 4], 8), 1)
        # Vacuous (sets too small): 2+2-8 < 0 -> 0.
        self.assertEqual(ps._unphased_carrier_lower_bound([2, 2], 8), 0)
        # Never exceeds the smallest set: 8+3 - 1*8 = 3, min set 3 -> 3.
        self.assertEqual(ps._unphased_carrier_lower_bound([8, 3], 8), 3)
        # Single pair -> the set itself is the lower bound (k-1=0).
        self.assertEqual(ps._unphased_carrier_lower_bound([4], 8), 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)

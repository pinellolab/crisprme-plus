"""Tests for the ADDITIVE companion population-summary WIRING (Phase 3c).

STDLIB ONLY (unittest + tempfile + os). Exercises the PURE, importable companion
writer (``population_summary_companion``) end-to-end against a REAL compiled Tier-0
registry + Tier-1 genotype store over a synthetic 2-db / 3-subpop panel, so the
frequencies are a true round-trip through ``population_summary.summarize`` (the
writer delegates all frequency math there -- these tests assert the wiring, not the
freq derivation, which ``test_population_summary`` already covers).

Asserts (per the Phase-3c spec):
  (a) a SINGLE-variant off-target row's global / per-db / per-subpop freqs in the
      companion row match ``population_summary.summarize`` exactly.
  (b) a MULTI-variant off-target row uses the COMBINATION (cis) frequency, not the
      per-SNP product.
  (c) a REFERENCE / PAM-only off-target (no SNP token) is SKIPPED (no companion row).
  (d) tier0_reader=None -> ``write_companion`` returns False and writes NO file
      (byte-identical legacy behavior).
Plus wiring guards:
  - the SNP-field parser recovers (pos, alt) pairs (single + multi + junk-skip);
  - the phased flag resolves (explicit > manifest > detect > conservative default);
  - the companion header's first columns are the bestMerge identity columns;
  - a per-row error is isolated (on_error called, file still written);
  - the unphased combination row carries a labeled carrier LOWER BOUND + NA AF.
"""

from __future__ import annotations

import os
import tempfile
import unittest

import population_summary as ps
import population_summary_companion as psc
import tier0_registry as t0
import tier1_genotypes as t1


# --------------------------------------------------------------------------- #
# Synthetic panel (same shape as test_population_summary): 2 db, 3 subpop.
# --------------------------------------------------------------------------- #
PANEL = {
    "S_EUR_M1": ("1000G", "EUR", "male"),
    "S_EUR_F1": ("1000G", "EUR", "female"),
    "S_EUR_M2": ("1000G", "EUR", "male"),
    "S_AFR_M1": ("1000G", "AFR", "male"),
    "S_AFR_F1": ("1000G", "AFR", "female"),
    "H_EAS_M1": ("HGDP", "EAS", "male"),
    "H_EAS_F1": ("HGDP", "EAS", "female"),
    "H_EAS_F2": ("HGDP", "EAS", "female"),
}
CHROM = "chr22"


def carriers_from(gt_map):
    out = {}
    for sid, gt in gt_map.items():
        toks = gt.replace("/", "|").split("|")
        if any(t == "1" for t in toks):
            out[sid] = gt
    return out


def _write_registry(tmpdir, records, name="t0"):
    binp = os.path.join(tmpdir, name + ".bin")
    idxp = os.path.join(tmpdir, name + ".idx.json")
    t0.compile_registry_panel(records, PANEL, None, t0.autosomal_ploidy, binp, idxp)
    return t0.RegistryReader(binp, idxp)


def _write_gt(tmpdir, records, axis, name="t1"):
    binp = os.path.join(tmpdir, name + ".bin")
    idxp = os.path.join(tmpdir, name + ".idx.json")
    t1.compile_genotypes(records, axis, binp, idxp)
    return t1.GenotypeReader(binp, idxp)


def _snp_token(pos, ref, alt, chrom=CHROM):
    return "%s_%d_%s_%s" % (chrom, pos, ref, alt)


def _off_target(chrom, position, direction, crrna, dna, snp):
    """A finalized-target-shaped record (dict) as the companion writer consumes."""
    return {
        "Chromosome": chrom,
        "Position": str(position),
        "Direction": direction,
        "crRNA": crrna,
        "DNA": dna,
        "SNP": snp,
    }


def _parse_group_col(cell):
    """Decode a "g1=v1;g2=v2" breakdown cell into {group: value_str}."""
    out = {}
    if not cell:
        return out
    for part in cell.split(";"):
        if "=" in part:
            g, v = part.split("=", 1)
            out[g] = v
    return out


def _read_companion(path):
    """Read the companion TSV into (header_list, list_of_row_dicts)."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    header = lines[0].lstrip("#").split("\t")
    rows = []
    for ln in lines[1:]:
        if not ln:
            continue
        rows.append(dict(zip(header, ln.split("\t"))))
    return header, rows


# --------------------------------------------------------------------------- #
# (a) SINGLE-variant off-target row matches summarize() exactly.
# --------------------------------------------------------------------------- #
class TestSingleVariantCompanion(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.autosomal_ploidy
        self.axis = t1.build_sample_axis(PANEL)
        # SNP at pos 100 alt A (phased 1000G-style '|').
        gts = {
            "S_EUR_M1": "1|0", "S_EUR_F1": "1|1", "S_AFR_M1": "0|1",
            "H_EAS_M1": "1|1", "H_EAS_F1": "1|0",
        }
        c = carriers_from(gts)
        self.t0r = _write_registry(self.d, [(100, "G", "A", "rs100", c)])
        self.t1r = _write_gt(self.d, [(100, "A", c, "G")], self.axis)

    def tearDown(self):
        self.t0r.close()
        self.t1r.close()

    def test_single_row_matches_summarize(self):
        snp = _snp_token(100, "G", "A")
        ot = _off_target(CHROM, 100, "+", "GUIDEAAAAAAAAAAAAAAA", "GUIDEAAAAAAAAAAAAAAt", snp)
        out = os.path.join(self.d, "out.population_summary.tsv")
        wrote = psc.write_companion(
            out, [ot], self.t0r, self.t1r, self.axis, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=True,
            global_group_id=t0.GLOBAL_GROUP_ID)
        self.assertTrue(wrote)
        header, rows = _read_companion(out)
        self.assertEqual(len(rows), 1)
        row = rows[0]

        # Identity columns are the bestMerge join key, verbatim.
        self.assertEqual(row["Chromosome"], CHROM)
        self.assertEqual(row["Position"], "100")
        self.assertEqual(row["Direction"], "+")
        self.assertEqual(row["crRNA"], "GUIDEAAAAAAAAAAAAAAA")
        self.assertEqual(row["DNA"], "GUIDEAAAAAAAAAAAAAAt")
        self.assertEqual(row["SNP"], snp)
        self.assertEqual(row["n_variants"], "1")

        # Ground truth straight from summarize().
        summary = ps.summarize([(100, "A")], self.t0r, self.t1r, self.axis,
                               self.ploidy_of, phased=True)
        self.assertEqual(row["global_carrier_n"], str(summary.global_carrier_n))
        self.assertAlmostEqual(float(row["global_allele_freq"]),
                               summary.global_af, places=6)
        self.assertAlmostEqual(float(row["global_carrier_freq"]),
                               summary.global_carrier_freq, places=6)
        self.assertAlmostEqual(float(row["global_hom_freq"]),
                               summary.global_hom_freq, places=6)
        self.assertEqual(row["observed"], "1" if summary.observed else "0")
        self.assertEqual(row["allele_freq_defined"], "1")
        self.assertEqual(row["max_subpop_label"],
                         summary.max_subpop_af_label or "")
        self.assertAlmostEqual(float(row["max_subpop_af"]),
                               summary.max_subpop_af, places=6)

        # Per-(db x subpop) + per-db breakdowns match each group's summary.
        af_by = _parse_group_col(row["allele_freq_by_group"])
        cf_by = _parse_group_col(row["carrier_freq_by_group"])
        n_by = _parse_group_col(row["carrier_n_by_group"])
        for gid, gs in summary.groups.items():
            if gid == t0.GLOBAL_GROUP_ID:
                continue
            self.assertIn(gid, af_by)
            self.assertAlmostEqual(float(af_by[gid]), gs.allele_freq, places=6)
            self.assertAlmostEqual(float(cf_by[gid]), gs.carrier_freq, places=6)
            self.assertEqual(n_by[gid], str(gs.n_carrier))
        # No global group leaked into the per-group breakdown columns.
        self.assertNotIn(t0.GLOBAL_GROUP_ID, af_by)


# --------------------------------------------------------------------------- #
# Two-SNP combination fixture (cis vs trans), phased + unphased.
# --------------------------------------------------------------------------- #
class _TwoSnp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ploidy_of = t0.autosomal_ploidy
        self.axis = t1.build_sample_axis(PANEL)
        self.gt100 = {
            "S_EUR_M1": "1|0", "S_EUR_F1": "1|1", "S_EUR_M2": "1|0",
            "S_AFR_M1": "0|1", "S_AFR_F1": "1|0", "H_EAS_M1": "1|0",
            "H_EAS_F1": "0|1",
        }
        self.gt200 = {
            "S_EUR_M1": "1|0", "S_EUR_F1": "1|1", "S_EUR_M2": "0|1",
            "S_AFR_M1": "0|1", "S_AFR_F1": "0|1", "H_EAS_M1": "1|0",
            "H_EAS_F2": "1|0",
        }
        c100 = carriers_from(self.gt100)
        c200 = carriers_from(self.gt200)
        self.t0r = _write_registry(
            self.d, [(100, "G", "A", "rs100", c100), (200, "T", "C", "rs200", c200)])
        self.t1_phased = _write_gt(
            self.d, [(100, "A", c100, "G"), (200, "C", c200, "T")], self.axis,
            name="t1_phased")

        def unph(m):
            return {sid: gt.replace("|", "/") for sid, gt in m.items()}
        self.t1_unphased = _write_gt(
            self.d, [(100, "A", unph(c100), "G"), (200, "C", unph(c200), "T")],
            self.axis, name="t1_unphased")

    def tearDown(self):
        self.t0r.close()
        self.t1_phased.close()
        self.t1_unphased.close()


# --------------------------------------------------------------------------- #
# (b) MULTI-variant off-target row uses the combination (cis) frequency.
# --------------------------------------------------------------------------- #
class TestMultiVariantCompanion(_TwoSnp):
    def test_combination_row_uses_cis_freq(self):
        snp = "%s,%s" % (_snp_token(100, "G", "A"), _snp_token(200, "T", "C"))
        ot = _off_target(CHROM, 100, "+", "GUIDE", "guide", snp)
        out = os.path.join(self.d, "combo.population_summary.tsv")
        psc.write_companion(
            out, [ot], self.t0r, self.t1_phased, self.axis, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=True, global_group_id=t0.GLOBAL_GROUP_ID)
        _header, rows = _read_companion(out)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["n_variants"], "2")

        # The combination (cis) summary is ground truth.
        summary = ps.summarize(
            [(100, "A"), (200, "C")], self.t0r, self.t1_phased, self.axis,
            self.ploidy_of, phased=True)
        self.assertEqual(row["global_carrier_n"], str(summary.global_carrier_n))
        self.assertAlmostEqual(float(row["global_allele_freq"]),
                               summary.global_af, places=6)

        # It is the COMBINATION freq, NOT the per-SNP product: only 4 cis carriers
        # globally (see test_population_summary::TestPhasedCombination), 5 cis haps.
        self.assertEqual(row["global_carrier_n"], "4")
        self.assertAlmostEqual(float(row["global_allele_freq"]), 5 / 16, places=6)

        # Per-group cis carrier counts (M1+F1 in EUR, M1 in AFR, M1 in EAS).
        n_by = _parse_group_col(row["carrier_n_by_group"])
        self.assertEqual(n_by["1000G::EUR"], "2")
        self.assertEqual(n_by["1000G::AFR"], "1")
        self.assertEqual(n_by["HGDP::EAS"], "1")

    def test_unphased_combination_row_marks_undefined_af_and_lower_bound(self):
        snp = "%s,%s" % (_snp_token(100, "G", "A"), _snp_token(200, "T", "C"))
        ot = _off_target(CHROM, 100, "+", "GUIDE", "guide", snp)
        out = os.path.join(self.d, "combo_unph.population_summary.tsv")
        psc.write_companion(
            out, [ot], self.t0r, self.t1_unphased, self.axis, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=False, global_group_id=t0.GLOBAL_GROUP_ID)
        _header, rows = _read_companion(out)
        row = rows[0]
        # allele freq UNDEFINED without phase -> NA + flag 0.
        self.assertEqual(row["global_allele_freq"], "NA")
        self.assertEqual(row["allele_freq_defined"], "0")
        self.assertEqual(row["phased"], "0")
        # per-group AF all NA, but carrier counts / freqs present.
        af_by = _parse_group_col(row["allele_freq_by_group"])
        for v in af_by.values():
            self.assertEqual(v, "NA")
        # carrier lower bound present + <= carrier freq (labeled column).
        self.assertIn("global_carrier_freq_lower_bound", row)
        self.assertNotEqual(row["global_carrier_freq_lower_bound"], "")
        self.assertLessEqual(float(row["global_carrier_freq_lower_bound"]),
                             float(row["global_carrier_freq"]))


# --------------------------------------------------------------------------- #
# (c) REFERENCE / PAM-only off-target is SKIPPED. + (d) tier0=None -> no file.
# --------------------------------------------------------------------------- #
class TestSkipAndGate(_TwoSnp):
    def test_reference_only_rows_skipped(self):
        variant = _off_target(CHROM, 100, "+", "G", "g",
                              _snp_token(100, "G", "A"))
        # Reference / PAM-only rows: empty SNP, "n" placeholder, "." placeholder.
        ref_empty = _off_target(CHROM, 5, "+", "G", "g", "")
        ref_n = _off_target(CHROM, 6, "+", "G", "g", "n")
        ref_dot = _off_target(CHROM, 7, "+", "G", "g", ".")
        out = os.path.join(self.d, "skip.population_summary.tsv")
        psc.write_companion(
            out, [ref_empty, variant, ref_n, ref_dot], self.t0r, self.t1_phased,
            self.axis, self.ploidy_of, ps, panel_cls=ps.Panel, phased=True,
            global_group_id=t0.GLOBAL_GROUP_ID)
        _header, rows = _read_companion(out)
        # Only the single VARIANT row survives; the 3 reference rows are skipped.
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Position"], "100")
        self.assertEqual(rows[0]["SNP"], _snp_token(100, "G", "A"))

    def test_no_registry_writes_nothing(self):
        ot = _off_target(CHROM, 100, "+", "G", "g", _snp_token(100, "G", "A"))
        out = os.path.join(self.d, "none.population_summary.tsv")
        wrote = psc.write_companion(
            out, [ot], None, self.t1_phased, self.axis, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=True, global_group_id=t0.GLOBAL_GROUP_ID)
        # GATE: tier0_reader None -> returns False and creates NO file.
        self.assertFalse(wrote)
        self.assertFalse(os.path.exists(out))

    def test_empty_off_targets_writes_header_only(self):
        out = os.path.join(self.d, "empty.population_summary.tsv")
        wrote = psc.write_companion(
            out, [], self.t0r, self.t1_phased, self.axis, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=True, global_group_id=t0.GLOBAL_GROUP_ID)
        self.assertTrue(wrote)  # registry present -> header-only file
        header, rows = _read_companion(out)
        self.assertEqual(rows, [])
        self.assertEqual(header[:5],
                         ["Chromosome", "Position", "Direction", "crRNA", "DNA"])


# --------------------------------------------------------------------------- #
# Per-row error isolation: a bad off-target must NOT sink the whole file.
# --------------------------------------------------------------------------- #
class TestPerRowErrorIsolation(_TwoSnp):
    def test_multi_variant_row_degrades_without_gt_tier(self):
        # A multi-variant (k>=2) off-target on a registry-only build (gt_reader=None,
        # axis=None) cannot compute a haplotype frequency, but must DEGRADE
        # gracefully (row EMITTED with undefined freqs) rather than be dropped via
        # on_error -- this INVERTS the old behavior where it raised and was skipped.
        good = _off_target(CHROM, 100, "+", "G", "g", _snp_token(100, "G", "A"))
        multi = _off_target(CHROM, 100, "+", "G", "g",
                            "%s,%s" % (_snp_token(100, "G", "A"),
                                       _snp_token(200, "T", "C")))
        seen = []
        out = os.path.join(self.d, "degrade.population_summary.tsv")
        wrote = psc.write_companion(
            out, [multi, good], self.t0r, None, None, self.ploidy_of, ps,
            panel_cls=ps.Panel, phased=True, global_group_id=t0.GLOBAL_GROUP_ID,
            on_error=lambda ot, err: seen.append((ot, err)))
        self.assertTrue(wrote)
        _header, rows = _read_companion(out)
        # BOTH rows emitted; on_error NOT called (graceful degrade, not an error).
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(seen), 0)
        _snps = {r["SNP"] for r in rows}
        self.assertIn(_snp_token(100, "G", "A"), _snps)             # single
        self.assertIn("%s,%s" % (_snp_token(100, "G", "A"),
                                 _snp_token(200, "T", "C")), _snps)  # combination

    def test_genuine_row_error_isolated(self):
        # A row that GENUINELY raises during summarize is still isolated via
        # on_error (the file is written; only the bad row is skipped) -- preserves
        # the isolation contract independent of the graceful-degradation path.
        ot = _off_target(CHROM, 100, "+", "G", "g", _snp_token(100, "G", "A"))
        seen = []

        def _boom(*a, **k):
            raise RuntimeError("synthetic summarize failure")

        _orig = ps.summarize
        ps.summarize = _boom
        out = os.path.join(self.d, "isolate.population_summary.tsv")
        try:
            wrote = psc.write_companion(
                out, [ot], self.t0r, self.t1_phased, self.axis, self.ploidy_of, ps,
                panel_cls=ps.Panel, phased=True, global_group_id=t0.GLOBAL_GROUP_ID,
                on_error=lambda o, err: seen.append((o, err)))
        finally:
            ps.summarize = _orig
        self.assertTrue(wrote)
        _header, rows = _read_companion(out)
        self.assertEqual(len(rows), 0)   # the only row raised -> isolated
        self.assertEqual(len(seen), 1)


# --------------------------------------------------------------------------- #
# SNP-field parser + phased resolution unit checks.
# --------------------------------------------------------------------------- #
class TestSnpFieldParser(unittest.TestCase):
    def test_single(self):
        self.assertEqual(psc.parse_snp_field("chr1_100_G_A"), [(100, "A")])

    def test_multi_combination_order_preserved(self):
        got = psc.parse_snp_field("chr1_100_G_A,chr1_200_T_C")
        self.assertEqual(got, [(100, "A"), (200, "C")])

    def test_chrom_with_underscores(self):
        # A contig name containing '_' must not corrupt pos/ref/alt (last 3 fields).
        got = psc.parse_snp_field("chr1_KI270711v1_random_100_G_A")
        self.assertEqual(got, [(100, "A")])

    def test_reference_placeholders_skipped(self):
        self.assertEqual(psc.parse_snp_field(""), [])
        self.assertEqual(psc.parse_snp_field("n"), [])
        self.assertEqual(psc.parse_snp_field("."), [])
        self.assertEqual(psc.parse_snp_field(None), [])

    def test_junk_token_skipped(self):
        got = psc.parse_snp_field("chr1_100_G_A,garbage,chr1_xx_T_C")
        self.assertEqual(got, [(100, "A")])


class TestPhasedResolution(unittest.TestCase):
    def test_explicit_wins(self):
        self.assertTrue(psc.resolve_phased(True, None))
        self.assertFalse(psc.resolve_phased(False, None))

    def test_detect_from_gt_strings(self):
        self.assertTrue(psc.resolve_phased(None, None, ["1|0", "0|1"]))
        self.assertFalse(psc.resolve_phased(None, None, ["1/0", "0/1"]))

    def test_conservative_default_when_unknown(self):
        # Nothing observed, no reader -> conservative default: unphased (False).
        self.assertFalse(psc.resolve_phased(None, None, []))

    def test_manifest_flag_used(self):
        class _FakeGT(object):
            manifest = {"phased": True}

            def gt_vocab(self):
                return ["0/1"]  # would say unphased, but the manifest wins
        self.assertTrue(psc.resolve_phased(None, _FakeGT()))

    def test_sniff_from_gt_vocab(self):
        class _FakeGT(object):
            manifest = {}

            def gt_vocab(self):
                return ["1|0", "0|1"]
        self.assertTrue(psc.resolve_phased(None, _FakeGT()))


class TestCompanionHeader(unittest.TestCase):
    def test_identity_columns_first(self):
        self.assertEqual(psc.COMPANION_HEADER[:5],
                         ["Chromosome", "Position", "Direction", "crRNA", "DNA"])

    def test_has_labeled_lower_bound_and_breakdowns(self):
        self.assertIn("global_carrier_freq_lower_bound", psc.COMPANION_HEADER)
        self.assertIn("allele_freq_by_group", psc.COMPANION_HEADER)
        self.assertIn("carrier_freq_by_group", psc.COMPANION_HEADER)
        self.assertIn("carrier_n_by_group", psc.COMPANION_HEADER)


if __name__ == "__main__":
    unittest.main(verbosity=2)

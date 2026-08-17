"""Tests for the ADDITIVE Tier-1 GENOTYPE-TIER wiring in new_simple_analysis
(CRISPRme+ dictless redesign, PHASE 2 step 2c).

This step wires the per-sample genotype tier (``tier1_genotypes.GenotypeReader``)
in as the ``sample_list`` provider so a TRUE dictless install (Tier-0 registry +
genotype tier, NO per-sample dict) reconstructs the legacy "Samples" column.

We exercise the PURE selection-matrix core
(``simple_analysis_registry.retrieve_5tuple`` with the NEW ``gtreader=`` argument),
which ``new_simple_analysis.retrieveFromDict`` delegates to (passing the
module-level ``mygt``). Testing the core directly avoids new_simple_analysis.py's
heavy module-level side effects WITHOUT changing runtime behavior.

New selection matrix asserted (WIRING spec):
  (a) registry + gtreader + NO dict            -> DICTLESS FULL: sample_list per alt
      EXACTLY equals the genotype tier's carrier_tokens, aligned to reader.alts_at.
  (b) registry + gtreader + NO dict, multiallelic -> per-alt carriers aligned, no
      cross-alt leak.
  (c) registry-only + NO gtreader               -> still samples=[] (degraded),
      UNCHANGED.
  (d) no registry                               -> legacy path unchanged (reuse the
      independent legacy re-derivation).
Plus:
  - gtreader present but the tier has NO record at a (pos, alt) the registry knows
    -> [] for that alt (carrier_tokens None), still aligned.
  - the dict-PRESENT augment path IGNORES gtreader (Samples come from the dict).
  - _resolve_genotype_paths mirrors _resolve_registry_paths (folder/file swap) and
    is detection-by-existence only.

STDLIB ONLY (unittest + tempfile). Uses a REAL compiled Tier-0 registry AND a REAL
compiled Tier-1 genotype store over the SAME synthetic (pos, alt, carriers) so both
mmap bisect readers are exercised end-to-end and the carrier tokens are a true
round-trip (no hand-mocked reader).
"""

from __future__ import annotations

import ast
import os
import tempfile
import unittest

import tier0_registry as t0
from tier0_registry import RegistryReader, GLOBAL_GROUP_ID, autosomal_ploidy

import tier1_genotypes as t1
import simple_analysis_registry as sar


# --------------------------------------------------------------------------- #
# Load the path resolvers WITHOUT importing new_simple_analysis (it opens
# sys.argv files at import time -- see the registry test's module docstring). We
# parse the source and exec ONLY the two resolver function defs into a namespace
# with ``os`` available, so we test the SHIPPING function bodies (not a copy) with
# zero module-level side effects.
# --------------------------------------------------------------------------- #
def _load_resolvers():
    src_path = os.path.join(os.path.dirname(__file__), "new_simple_analysis.py")
    with open(src_path, "r") as fh:
        tree = ast.parse(fh.read(), filename=src_path)
    wanted = {"_resolve_registry_paths", "_resolve_genotype_paths"}
    defs = [n for n in tree.body
            if isinstance(n, ast.FunctionDef) and n.name in wanted]
    module = ast.Module(body=defs, type_ignores=[])
    ns = {"os": os}
    exec(compile(module, src_path, "exec"), ns)
    return ns


_RESOLVERS = _load_resolvers()
_resolve_genotype_paths = _RESOLVERS["_resolve_genotype_paths"]
_resolve_registry_paths = _RESOLVERS["_resolve_registry_paths"]


# --------------------------------------------------------------------------- #
# Fixtures: a tiny panel shared by the registry AND the genotype tier.
# --------------------------------------------------------------------------- #
SAMPLE_META = {
    "S1": ("1000G", "EUR", "male"),
    "S2": ("1000G", "EUR", "female"),
    "S3": ("1000G", "AFR", "male"),
    "S4": ("HGDP", "EAS", "female"),
}
CHROM = "chr22"


def _build_registry(records, tmpdir, name="reg_chr22"):
    """Compile ``records`` (pos, ref, alt, rsid, alt_genotypes) into a panel-aware
    registry, return an open reader."""
    out_bin = os.path.join(tmpdir, name + ".bin")
    out_idx = os.path.join(tmpdir, name + ".idx")
    t0.compile_registry_panel(
        records, SAMPLE_META, None, autosomal_ploidy, out_bin, out_idx,
    )
    return RegistryReader(out_bin, out_idx)


def _build_gtstore(records, tmpdir, name="gt_chr22"):
    """Compile the SAME synthetic (pos, alt, carriers) into a Tier-1 genotype store,
    return an open GenotypeReader.

    ``records``: iterable of (pos:int, ref, alt, rsid, alt_genotypes:dict). We reuse
    the registry-shaped records and project them to the genotype tier's
    (pos, alt, carriers, ref) 4-tuple so the two stores describe the identical panel
    at the identical (pos, alt) keys.
    """
    axis = t1.build_sample_axis(SAMPLE_META)
    gt_records = [(pos, alt, dict(gts), ref)
                  for (pos, ref, alt, _rsid, gts) in records]
    out_bin = os.path.join(tmpdir, name + ".bin")
    out_idx = os.path.join(tmpdir, name + ".idx.json")
    t1.compile_genotypes(gt_records, axis, out_bin, out_idx)
    return t1.GenotypeReader(out_bin, out_idx)


def _expected_tokens(carrier_dict, axis):
    """The genotype tier's deterministic carrier_tokens for a carrier dict:
    "sid:gt" in ascending sample-index order."""
    ordered = sorted(carrier_dict.items(), key=lambda kv: axis.index_of(kv[0]))
    return ["%s:%s" % (sid, gt) for sid, gt in ordered]


# --------------------------------------------------------------------------- #
# Independent re-derivation of the ORIGINAL retrieveFromDict body (pre-wiring),
# to prove the no-registry legacy path is byte-identical (constraint (d)).
# --------------------------------------------------------------------------- #
def _legacy_reference(entry, current_chr, chr_pos):
    if entry is None:
        return (["C"], [[]], ["."], ["0"],
                [current_chr + "_" + str(chr_pos + 1) + "_" + "C" + "_" + "G"])
    snp_list, sample_list, AF_list, rsID_list, snp_info_list = [], [], [], [], []
    for e in entry.split("$"):
        split_entry = e.split(";")
        samples = split_entry[0].strip().split(",")
        if samples[0] == "":
            samples = []
        sample_list.append(samples)
        snp_list.append(split_entry[1].strip().split(",")[1])
        rsID_list.append(split_entry[2].strip())
        AF_list.append(split_entry[3].strip())
        snp_info_list.append(
            current_chr + "_" + str(chr_pos + 1) + "_"
            + split_entry[1].split(",")[0] + "_" + split_entry[1].split(",")[1]
        )
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


def _dict_entry(*per_alt):
    """Legacy '$'-joined dict value; each per_alt=(samples, ref, alt, rsid, af)."""
    return "$".join(
        "{};{},{};{};{}".format(samples, ref, alt, rsid, af)
        for (samples, ref, alt, rsid, af) in per_alt
    )


class TestGenotypeTierWiring(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self._closeables = []
        self.axis = t1.build_sample_axis(SAMPLE_META)

    def tearDown(self):
        for r in self._closeables:
            try:
                r.close()
            except Exception:
                pass
        self._tmp.cleanup()

    def _reg(self, records):
        r = _build_registry(records, self.tmpdir)
        self._closeables.append(r)
        return r

    def _gt(self, records):
        r = _build_gtstore(records, self.tmpdir)
        self._closeables.append(r)
        return r

    # ---- (a) registry + gtreader + NO dict: DICTLESS FULL ----------------- #
    def test_registry_gtreader_no_dict_biallelic(self):
        pos1 = 1000
        chr_pos = pos1 - 1
        carriers = {"S1": "1|0", "S3": "0|1"}
        records = [(pos1, "G", "A", "rs777", carriers)]
        reg = self._reg(records)
        gt = self._gt(records)

        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID, gtreader=gt
        )

        self.assertEqual(snp, ["A"])
        self.assertEqual(rsid, ["rs777"])
        # AF still from the registry (AC/AN over the full panel): 2/8 = 0.25.
        self.assertAlmostEqual(float(af[0]), 0.25)
        self.assertEqual(info, [CHROM + "_1000_G_A"])
        # Samples reconstructed from the GENOTYPE TIER, EXACTLY.
        expect = _expected_tokens(carriers, self.axis)
        self.assertEqual(samples, [expect])
        # And identical to what the tier itself returns for this (pos, alt).
        self.assertEqual(samples[0], gt.carrier_tokens(pos1, "A"))

    # ---- (b) registry + gtreader + NO dict, MULTIALLELIC ------------------ #
    def test_registry_gtreader_no_dict_multiallelic_no_leak(self):
        pos1 = 5000
        chr_pos = pos1 - 1
        cA = {"S1": "1|0"}
        cT = {"S3": "1|0", "S4": "0|1"}
        records = [
            (pos1, "G", "A", "rsA", cA),
            (pos1, "G", "T", "rsT", cT),
        ]
        reg = self._reg(records)
        gt = self._gt(records)

        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID, gtreader=gt
        )

        # Dictless order is reader.alts_at (sorted): A then T.
        self.assertEqual(snp, ["A", "T"])
        self.assertEqual(rsid, ["rsA", "rsT"])
        self.assertAlmostEqual(float(af[0]), 1.0 / 8.0)   # A: S1 het
        self.assertAlmostEqual(float(af[1]), 2.0 / 8.0)   # T: two het alleles
        self.assertEqual(info, [CHROM + "_5000_G_A", CHROM + "_5000_G_T"])
        # Per-alt carriers aligned; NO cross-alt leak.
        self.assertEqual(samples[0], _expected_tokens(cA, self.axis))
        self.assertEqual(samples[1], _expected_tokens(cT, self.axis))
        self.assertEqual(samples[0], gt.carrier_tokens(pos1, "A"))
        self.assertEqual(samples[1], gt.carrier_tokens(pos1, "T"))
        # Sanity: the two carrier SETS are disjoint (proving no leak).
        self.assertEqual(set(samples[0]) & set(samples[1]), set())

    def test_registry_gtreader_no_dict_alt_missing_in_tier(self):
        # The registry knows an alt the genotype tier has NO record for (e.g. a
        # carrier-less alt the tier skipped): carrier_tokens -> None -> [] for that
        # alt, still aligned to the other alts. We simulate by compiling the tier
        # with ONLY the "A" record while the registry carries both alts.
        pos1 = 5500
        chr_pos = pos1 - 1
        cA = {"S2": "1|1"}
        reg = self._reg([
            (pos1, "C", "A", "rsA", cA),
            (pos1, "C", "G", "rsG", {"S4": "0|1"}),
        ])
        gt = self._gt([(pos1, "C", "A", "rsA", cA)])  # tier lacks the "G" record

        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID, gtreader=gt
        )
        self.assertEqual(snp, ["A", "G"])
        self.assertEqual(samples[0], _expected_tokens(cA, self.axis))
        self.assertIsNone(gt.carrier_tokens(pos1, "G"))
        self.assertEqual(samples[1], [])  # None -> [] fallback, still aligned

    # ---- (c) registry-only + NO gtreader: still degraded [] (unchanged) --- #
    def test_registry_only_no_gtreader_degraded(self):
        pos1 = 2000
        chr_pos = pos1 - 1
        reg = self._reg([(pos1, "C", "T", "rs42", {"S2": "0|1"})])
        # gtreader defaults to None -> UNCHANGED registry-only behavior.
        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID
        )
        self.assertEqual(snp, ["T"])
        self.assertAlmostEqual(float(af[0]), 0.125)
        self.assertEqual(rsid, ["rs42"])
        self.assertEqual(samples, [[]])   # degraded, unchanged
        self.assertEqual(info, [CHROM + "_2000_C_T"])
        # Passing gtreader=None explicitly is identical.
        self.assertEqual(
            sar.retrieve_5tuple(reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID,
                                gtreader=None),
            (snp, samples, rsid, af, info),
        )

    # ---- (d) no registry: legacy path unchanged (even WITH a gtreader) ---- #
    def test_no_registry_legacy_unchanged_hit(self):
        pos1 = 3000
        chr_pos = pos1 - 1
        entry = _dict_entry(("S1:1|0,S4:0|1", "G", "A", "rs9", "0.03"))
        gt = self._gt([(pos1, "G", "A", "rs9", {"S1": "1|0", "S4": "0|1"})])
        exp = _legacy_reference(entry, CHROM, chr_pos)
        # reader None -> legacy branch; a stray gtreader must NOT be consulted.
        got = sar.retrieve_5tuple(None, entry, CHROM, chr_pos, GLOBAL_GROUP_ID,
                                  gtreader=gt)
        self.assertEqual(got, exp)

    def test_no_registry_legacy_unchanged_miss(self):
        pos1 = 4000
        chr_pos = pos1 - 1
        exp = _legacy_reference(None, CHROM, chr_pos)
        got = sar.retrieve_5tuple(None, None, CHROM, chr_pos, GLOBAL_GROUP_ID,
                                  gtreader=None)
        self.assertEqual(got, exp)

    # ---- dict PRESENT + gtreader present: dict wins (augment path) -------- #
    def test_dict_present_ignores_gtreader(self):
        pos1 = 6000
        chr_pos = pos1 - 1
        carriers = {"S1": "1|0", "S3": "0|1"}
        records = [(pos1, "G", "A", "rs777", carriers)]
        reg = self._reg(records)
        # Genotype tier with DIFFERENT carriers, to prove the dict path ignores it.
        gt = self._gt([(pos1, "G", "A", "rs777", {"S2": "1|1"})])
        entry = _dict_entry(("S1:1|0,S3:0|1", "G", "A", "rs777", ""))
        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, entry, CHROM, chr_pos, GLOBAL_GROUP_ID, gtreader=gt
        )
        # Samples come from the DICT, NOT the genotype tier.
        self.assertEqual(samples, [["S1:1|0", "S3:0|1"]])
        # AF still overridden from the registry.
        self.assertAlmostEqual(float(af[0]), 0.25)


# --------------------------------------------------------------------------- #
# _resolve_genotype_paths: mirrors _resolve_registry_paths (folder/file swap),
# detection-by-existence only.
# --------------------------------------------------------------------------- #
class TestResolveGenotypePaths(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def _layout(self, vcf="1000G", chrom="chr22"):
        dicts_dir = os.path.join(self.root, "Dictionaries",
                                 "dictionaries_" + vcf)
        gt_dir = os.path.join(self.root, "Dictionaries", "genotypes_" + vcf)
        os.makedirs(dicts_dir)
        os.makedirs(gt_dir)
        dict_path = os.path.join(dicts_dir, "my_dict_" + chrom + ".json")
        return dict_path, gt_dir

    def test_absent_returns_none(self):
        dict_path, _gt_dir = self._layout()
        self.assertIsNone(_resolve_genotype_paths(dict_path, "chr22"))

    def test_present_returns_bin_idx(self):
        dict_path, gt_dir = self._layout()
        binp = os.path.join(gt_dir, "gt_chr22.bin")
        idxp = os.path.join(gt_dir, "gt_chr22.idx")
        open(binp, "w").close()
        open(idxp, "w").close()
        got = _resolve_genotype_paths(dict_path, "chr22")
        self.assertEqual(got, (binp, idxp))

    def test_only_bin_present_returns_none(self):
        # Both files must exist (detection-by-existence, like the registry).
        dict_path, gt_dir = self._layout()
        open(os.path.join(gt_dir, "gt_chr22.bin"), "w").close()
        self.assertIsNone(_resolve_genotype_paths(dict_path, "chr22"))

    def test_mirrors_registry_folder_swap(self):
        # The genotype resolver's folder/file swap is the registry resolver's, with
        # genotypes_/gt_ in place of registry_/reg_. Compile the SAME dict layout
        # and confirm the genotype paths track the registry paths byte-for-byte
        # (folder prefix + file stem swapped) at the same chromosome.
        dict_path, gt_dir = self._layout(vcf="hg38_1000G_HGDP")
        dict_path_chr1 = os.path.join(
            os.path.dirname(dict_path), "my_dict_chr1.json")

        # Genotype store present (gt_dir already created by _layout).
        open(os.path.join(gt_dir, "gt_chr1.bin"), "w").close()
        open(os.path.join(gt_dir, "gt_chr1.idx"), "w").close()
        got_gt = _resolve_genotype_paths(dict_path_chr1, "chr1")
        self.assertEqual(
            got_gt,
            (os.path.join(gt_dir, "gt_chr1.bin"),
             os.path.join(gt_dir, "gt_chr1.idx")),
        )

        # Registry store present in the mirror location; confirm the two resolvers
        # differ ONLY by the genotypes_/gt_ vs registry_/reg_ swap.
        reg_dir = os.path.join(self.root, "Dictionaries",
                               "registry_hg38_1000G_HGDP")
        os.makedirs(reg_dir)
        open(os.path.join(reg_dir, "reg_chr1.bin"), "w").close()
        open(os.path.join(reg_dir, "reg_chr1.idx"), "w").close()
        got_reg = _resolve_registry_paths(dict_path_chr1, "chr1")
        self.assertEqual(
            [p.replace("genotypes_", "registry_").replace("/gt_", "/reg_")
             for p in got_gt],
            list(got_reg),
        )


if __name__ == "__main__":
    unittest.main()

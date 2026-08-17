"""Tests for the ADDITIVE Tier-0 registry wiring in new_simple_analysis
(CRISPRme+ dictless redesign, Phase 1 FINAL step).

We exercise the PURE, importable retrieveFromDict machinery
(``simple_analysis_registry.retrieve_5tuple`` and its verbatim-legacy helpers),
which new_simple_analysis.retrieveFromDict delegates to. Testing the helper
directly avoids new_simple_analysis.py's heavy module-level side effects (it opens
sys.argv files at import) WITHOUT changing runtime behavior -- retrieveFromDict is
now a thin wrapper that supplies (myreg, dict-entry, current_chr).

Selection matrix asserted (WIRING spec):
  (a) registry+dict   -> AF from registry (AC/AN), samples from dict, ordering aligned
  (b) registry-only   -> AF from registry, samples=[] per alt
  (c) no registry     -> legacy path unchanged (mydict entry parsed exactly as before)
  (d) multiallelic    -> one 5-tuple element per alt, registry vs dict aligned

STDLIB ONLY (unittest + tempfile). Uses a REAL compiled registry (panel-aware)
so AF is the true AC/AN and the mmap bisect reader is exercised end-to-end.
"""

import os
import tempfile
import unittest

import tier0_registry as t0
from tier0_registry import RegistryReader, GLOBAL_GROUP_ID, autosomal_ploidy

import simple_analysis_registry as sar


# --------------------------------------------------------------------------- #
# Fixtures: a tiny panel-aware registry + matching legacy dict entries.
# --------------------------------------------------------------------------- #
# Panel: 4 diploid samples across two (db x subpop) cells.
SAMPLE_META = {
    "S1": ("1000G", "EUR", "male"),
    "S2": ("1000G", "EUR", "female"),
    "S3": ("1000G", "AFR", "male"),
    "S4": ("HGDP", "EAS", "female"),
}
CHROM = "chr22"


def _build_registry(records, tmpdir, name="reg_chr22"):
    """Compile ``records`` into a panel-aware registry, return an open reader.

    records: iterable of (pos:int, ref, alt, rsid, alt_genotypes:dict).
    """
    out_bin = os.path.join(tmpdir, name + ".bin")
    out_idx = os.path.join(tmpdir, name + ".idx")
    t0.compile_registry_panel(
        records, SAMPLE_META, None, autosomal_ploidy, out_bin, out_idx,
    )
    return RegistryReader(out_bin, out_idx)


def _dict_entry(*per_alt):
    """Build a legacy '$'-joined dict value from per-alt tuples.

    each per_alt: (samples_field, ref, alt, rsid, af)
      samples_field = "sampleID:genotype,sampleID:genotype" (only carriers) or "".
    Format per entry: "<samples>;<ref,alt>;<rsID>;<AF>".
    """
    return "$".join(
        "{};{},{};{};{}".format(samples, ref, alt, rsid, af)
        for (samples, ref, alt, rsid, af) in per_alt
    )


class TestSelectionMatrix(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self._readers = []

    def tearDown(self):
        for r in self._readers:
            try:
                r.close()
            except Exception:
                pass
        self._tmp.cleanup()

    def _reg(self, records):
        r = _build_registry(records, self.tmpdir)
        self._readers.append(r)
        return r

    # ---- (a) registry + dict: AF from registry, samples from dict --------- #
    def test_registry_plus_dict(self):
        pos1 = 1000  # 1-based genomic position
        chr_pos = pos1 - 1  # retrieveFromDict takes 0-based
        # Registry: at pos 1000, alt "A" carried het by S1 and S3.
        reg = self._reg([
            (pos1, "G", "A", "rs777", {"S1": "1|0", "S3": "0|1"}),
        ])
        # Legacy dict entry (carriers-only); its AF field is the OLD (empty/wrong)
        # one -- we must IGNORE it and take the registry AC/AN.
        entry = _dict_entry(("S1:1|0,S3:0|1", "G", "A", "rs777", ""))

        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, entry, CHROM, chr_pos, GLOBAL_GROUP_ID
        )

        self.assertEqual(snp, ["A"])
        # AF from registry: 2 alt alleles / (4 diploid samples * 2) = 2/8 = 0.25.
        groups = reg.lookup(pos1, "A")
        self.assertAlmostEqual(float(af[0]), 0.25)
        self.assertAlmostEqual(groups[GLOBAL_GROUP_ID].allele_freq(), 0.25)
        self.assertEqual(rsid, ["rs777"])
        # Samples come from the DICT (carrier tokens, per-alt aligned).
        self.assertEqual(samples, [["S1:1|0", "S3:0|1"]])
        self.assertEqual(info, [CHROM + "_1000_N_A"])

    # ---- (b) registry only: AF from registry, samples=[] per alt ---------- #
    def test_registry_only(self):
        pos1 = 2000
        chr_pos = pos1 - 1
        reg = self._reg([
            (pos1, "C", "T", "rs42", {"S2": "0|1"}),
        ])
        # No dict entry for this position.
        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID
        )
        self.assertEqual(snp, ["T"])
        # 1 alt allele / 8 = 0.125.
        self.assertAlmostEqual(float(af[0]), 0.125)
        self.assertEqual(rsid, ["rs42"])
        # Degraded Samples: [] per alt (genotype tier is a later phase).
        self.assertEqual(samples, [[]])
        self.assertEqual(info, [CHROM + "_2000_N_T"])

    # ---- (c) no registry: legacy path unchanged --------------------------- #
    def test_no_registry_legacy_hit(self):
        pos1 = 3000
        chr_pos = pos1 - 1
        entry = _dict_entry(("S1:1|0,S4:0|1", "G", "A", "rs9", "0.03"))
        # Independent legacy re-derivation (the ORIGINAL retrieveFromDict body).
        exp = _legacy_reference(entry, CHROM, chr_pos)
        got = sar.retrieve_5tuple(None, entry, CHROM, chr_pos, GLOBAL_GROUP_ID)
        self.assertEqual(got, exp)
        # Concretely: alt "A", samples verbatim, rsID "rs9", AF is the DICT's
        # "0.03" (NOT recomputed -- there is no registry), info ref/alt from dict.
        snp, samples, rsid, af, info = got
        self.assertEqual(snp, ["A"])
        self.assertEqual(samples, [["S1:1|0", "S4:0|1"]])
        self.assertEqual(rsid, ["rs9"])
        self.assertEqual(af, ["0.03"])
        self.assertEqual(info, [CHROM + "_3000_G_A"])

    def test_no_registry_legacy_miss(self):
        pos1 = 4000
        chr_pos = pos1 - 1
        # No registry AND no dict entry -> the legacy fake-C>G fallback, verbatim.
        exp = _legacy_reference(None, CHROM, chr_pos)
        got = sar.retrieve_5tuple(None, None, CHROM, chr_pos, GLOBAL_GROUP_ID)
        self.assertEqual(got, exp)
        snp, samples, rsid, af, info = got
        self.assertEqual(snp, ["C"])
        self.assertEqual(samples, [[]])
        self.assertEqual(rsid, ["."])
        self.assertEqual(af, ["0"])
        self.assertEqual(info, [CHROM + "_4000_C_G"])

    # ---- (d) multiallelic: one element per alt, registry vs dict aligned -- #
    def test_multiallelic_registry_plus_dict(self):
        pos1 = 5000
        chr_pos = pos1 - 1
        # Two alts at the same pos: "A" (carried by S1) and "T" (carried by S3,S4).
        reg = self._reg([
            (pos1, "G", "A", "rsA", {"S1": "1|0"}),
            (pos1, "G", "T", "rsT", {"S3": "1|0", "S4": "0|1"}),
        ])
        # Dict lists the two alts in a DIFFERENT order than the registry's sorted
        # (A before T) order, to prove we align by alt, not by position in the dict.
        entry = _dict_entry(
            ("S3:1|0,S4:0|1", "G", "T", "rsT", ""),
            ("S1:1|0", "G", "A", "rsA", ""),
        )
        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, entry, CHROM, chr_pos, GLOBAL_GROUP_ID
        )
        # Registry emits alts in sorted on-disk order: A, then T.
        self.assertEqual(snp, ["A", "T"])
        self.assertEqual(rsid, ["rsA", "rsT"])
        # AF: A = 1/8 = 0.125 ; T = 2/8 = 0.25.
        self.assertAlmostEqual(float(af[0]), 0.125)
        self.assertAlmostEqual(float(af[1]), 0.25)
        # Samples aligned to the registry alt order (NOT the dict's entry order).
        self.assertEqual(samples, [["S1:1|0"], ["S3:1|0", "S4:0|1"]])
        self.assertEqual(info, [CHROM + "_5000_N_A", CHROM + "_5000_N_T"])

    def test_multiallelic_registry_only(self):
        pos1 = 6000
        chr_pos = pos1 - 1
        reg = self._reg([
            (pos1, "C", "A", "rsA", {"S1": "1|0"}),
            (pos1, "C", "G", "rsG", {"S2": "1|1"}),
        ])
        snp, samples, rsid, af, info = sar.retrieve_5tuple(
            reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID
        )
        self.assertEqual(snp, ["A", "G"])
        self.assertEqual(rsid, ["rsA", "rsG"])
        self.assertAlmostEqual(float(af[0]), 1.0 / 8.0)     # S1 het
        self.assertAlmostEqual(float(af[1]), 2.0 / 8.0)     # S2 hom
        self.assertEqual(samples, [[], []])                 # degraded Samples

    # ---- registry present but NO record at pos -> legacy fallback --------- #
    def test_registry_present_position_absent_falls_back_to_dict(self):
        pos1 = 7000
        chr_pos = pos1 - 1
        # Registry has a record elsewhere, but nothing at pos 7000.
        reg = self._reg([
            (9999, "G", "A", "rsX", {"S1": "1|0"}),
        ])
        entry = _dict_entry(("S2:0|1", "G", "A", "rsLegacy", "0.5"))
        # Because the registry has no record at 7000, we must NOT lose the dict's
        # off-target: fall back to the exact legacy decode.
        exp = _legacy_reference(entry, CHROM, chr_pos)
        got = sar.retrieve_5tuple(reg, entry, CHROM, chr_pos, GLOBAL_GROUP_ID)
        self.assertEqual(got, exp)

    def test_registry_present_position_absent_no_dict_fake_snp(self):
        pos1 = 8000
        chr_pos = pos1 - 1
        reg = self._reg([
            (9999, "G", "A", "rsX", {"S1": "1|0"}),
        ])
        exp = _legacy_reference(None, CHROM, chr_pos)
        got = sar.retrieve_5tuple(reg, None, CHROM, chr_pos, GLOBAL_GROUP_ID)
        self.assertEqual(got, exp)


# --------------------------------------------------------------------------- #
# Independent re-derivation of the ORIGINAL retrieveFromDict body (pre-wiring).
# Kept here, written a different way, so the legacy path is proven byte-identical
# rather than merely "equal to my own helper".
# --------------------------------------------------------------------------- #
def _legacy_reference(entry, current_chr, chr_pos):
    if entry is None:
        snp_list = []
        sample_list = []
        AF_list = []
        rsID_list = []
        snp_info_list = []
        sample_list.append([])
        snp_list.append("C")
        rsID_list.append(".")
        AF_list.append("0")
        snp_info_list.append(
            current_chr + "_" + str(chr_pos + 1) + "_" + "C" + "_" + "G"
        )
        return snp_list, sample_list, rsID_list, AF_list, snp_info_list
    multi_entry = entry.split("$")
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for e in multi_entry:
        split_entry = e.split(";")
        samples = split_entry[0].strip().split(",")
        if samples[0] == "":
            samples = []
        sample_list.append(samples)
        snp_list.append(split_entry[1].strip().split(",")[1])
        rsID_list.append(split_entry[2].strip())
        AF_list.append(split_entry[3].strip())
        snp_info_list.append(
            current_chr
            + "_"
            + str(chr_pos + 1)
            + "_"
            + split_entry[1].split(",")[0]
            + "_"
            + split_entry[1].split(",")[1]
        )
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


class TestLegacyHelperVerbatim(unittest.TestCase):
    """The legacy helpers in simple_analysis_registry are a verbatim lift; assert
    they reproduce the independent reference across hit / miss / multiallelic."""

    def test_hit(self):
        entry = "S1:1|0,S2:0|1;G,A;rs1;0.1$;C,T;.;"
        self.assertEqual(
            sar.legacy_5tuple_from_entry(entry, "chr1", 41),
            _legacy_reference(entry, "chr1", 41),
        )

    def test_miss(self):
        self.assertEqual(
            sar.no_entry_5tuple("chr1", 41),
            _legacy_reference(None, "chr1", 41),
        )

    def test_empty_carrier_field(self):
        # An entry whose carrier field is empty -> samples == [] (samples[0] == "").
        entry = ";G,A;rs1;0.2"
        got = sar.legacy_5tuple_from_entry(entry, "chrX", 0)
        self.assertEqual(got, _legacy_reference(entry, "chrX", 0))
        self.assertEqual(got[1], [[]])


if __name__ == "__main__":
    unittest.main()

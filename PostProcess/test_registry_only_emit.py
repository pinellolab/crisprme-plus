"""Integration test for the REGISTRY-ONLY degraded-emit path of
``new_simple_analysis.iupac_decomposition`` (#136).

A registry-only install has a Tier-0 registry (for AF / rsID / metadata) but NO
per-sample dict AND NO Tier-1 genotype tier (e.g. ``download --no-genotypes``).
Historically the variant finalizer's ``len(samples) > 0`` guard DROPPED such
variant-created off-target sites entirely (the carrier list is unavoidably empty).
The fix emits the site anyway with a degraded Samples column ("NA"), still carrying
the creating-variant's registry rsID / AF, gated on ``registry_only_mode`` so every
other install mode stays byte-identical.

These drive the REAL ``iupac_decomposition`` through the same AST-load harness the
observed-haplotype integration test uses (exec only the pre-``inFasta`` prologue,
inject runtime globals), but with ``myreg`` = a REAL compiled registry, ``mygt`` =
None, ``mydict`` = {} -- i.e. mode 3.

STDLIB ONLY (+ the AST harness's numpy/pandas/CRISTA stubs) + tier0_registry.

Run with:
    cd PostProcess && python3 -m unittest test_registry_only_emit -v
"""

import io
import os
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN

import tier0_registry as t0
from tier0_registry import RegistryReader, GLOBAL_GROUP_ID, autosomal_ploidy

# 4 diploid samples across two (db x subpop) cells (mirrors the registry unit test).
SAMPLE_META = {
    "S1": ("1000G", "EUR", "male"),
    "S2": ("1000G", "EUR", "female"),
    "S3": ("1000G", "AFR", "male"),
    "S4": ("HGDP", "EAS", "female"),
}


def _build_registry(records, tmpdir):
    out_bin = os.path.join(tmpdir, "reg_chrT.bin")
    out_idx = os.path.join(tmpdir, "reg_chrT.idx")
    t0.compile_registry_panel(
        records, SAMPLE_META, None, autosomal_ploidy, out_bin, out_idx,
    )
    return RegistryReader(out_bin, out_idx)


def _build_mode3(reg, registry_only_mode):
    """A single variant-created off-target: ref 'C' at column 0 (mismatch vs the
    all-'A' guide), alt 'A' (matches -> 0-mismatch after substitution), NGG PAM.
    Registry supplies the alt at pos1=1 with a real carrier (for AF), but with NO
    genotype tier the Samples column cannot be resolved."""
    ref = ["A"] * TARGET_LEN
    ref[0] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"
    genome = "".join(ref)

    dna = list(genome)
    dna[0] = "W"  # IUPAC ambiguity -> variant-decomposition branch
    dna = "".join(dna)

    split = [
        "X", GUIDE, dna, "chrT", "0", "0", "+",
        "1", "0", "1", "NGG", "y", "NA",
        "NA", "NA", "NA", "NA", "NA",
    ]
    overrides = dict(
        genomeStr=genome,
        current_chr="chrT",
        mydict={},                    # mode 3: NO per-sample dict
        myreg=reg,                    # Tier-0 registry present
        mygt=None,                    # mode 3: NO genotype tier
        haplotype_check=False,        # no dict -> unphased/empty
        dict_tier_present=False,
        registry_only_mode=registry_only_mode,
        IUPAC_CAP=10,
        hvdr_bed=io.StringIO(),
        pam="NGG",
        pos_beg=0,
        pos_end=-3,
        pam_begin=-3,
        pam_end=None,
        allowed_mms=6,
    )
    return split, overrides


class TestRegistryOnlyEmit(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        # registry: pos1=1, ref 'C' -> alt 'A', one het carrier S2 (AF = 1/8).
        self.reg = _build_registry(
            [(1, "C", "A", "rs42", {"S2": "0|1"})], self._tmp.name
        )

    def tearDown(self):
        try:
            self.reg.close()
        except Exception:
            pass
        self._tmp.cleanup()

    def _run(self, registry_only_mode):
        split, overrides = _build_mode3(self.reg, registry_only_mode)
        ns = _load_pure_functions(overrides)
        ns["_phase_confirmation_rows"].clear()
        ns["_phase_confirmation_keys"].clear()
        cluster = []
        ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
        return cluster

    def test_registry_only_emits_degraded_variant_row(self):
        """mode 3 with registry_only_mode ON: the site is EMITTED, Samples='NA',
        and the registry rsID/AF ride along."""
        cluster = self._run(registry_only_mode=True)
        self.assertEqual(
            len(cluster), 1, "registry-only should emit exactly one degraded row"
        )
        row = cluster[0]
        self.assertEqual(row[12], "NA", "degraded Samples sentinel")
        self.assertEqual(row[15], "rs42", "registry rsID surfaced")
        self.assertAlmostEqual(float(row[16]), 0.125, msg="registry AC/AN AF")
        self.assertEqual(row[17], "chrT_1_C_A", "creating-variant identity")
        # ref-score sentinel tail intact (scoring not corrupted).
        self.assertEqual(row[-2], 33)

    def test_registry_only_mode_off_drops_variant_row(self):
        """REGRESSION: with the flag OFF, the empty-carrier variant row is DROPPED
        exactly as before the fix -- so no other install mode changes."""
        cluster = self._run(registry_only_mode=False)
        self.assertEqual(
            len(cluster), 0,
            "flag off must preserve the legacy drop (byte-identical behaviour)",
        )


if __name__ == "__main__":
    unittest.main()

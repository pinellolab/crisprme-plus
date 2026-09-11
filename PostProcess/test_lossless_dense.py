"""P1 (2.5.2): CRISPRME_LOSSLESS_DENSE emits the maximal co-located PUTATIVE haplotype in
a CAPPED registry-only dense window, so a multi-variant off-target the min-mismatch greedy
drops is not MISSED. Default OFF -> byte-identical (the extra row is not emitted).

Fixture: all-'A' guide; cols 0,1,2 are ref 'C' -> alt 'A' (the alt LOWERS mismatch, so the
greedy takes them), col 3 is ref 'A' -> alt 'T' (the alt RAISES mismatch, so the greedy
leaves it at reference). The greedy rep is therefore all-'A' (0 mm) and the genuine
4-variant haplotype (col3='T', 1 mm, still within budget) is DROPPED -- unless lossless-dense
emits the co-located union. IUPAC_CAP=2 < 4 variants forces the capped path."""
import io
import os
import sys
import unittest
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN
from test_registry_only_emit import _build_registry


def _dense_fixture(reg, lossless):
    ref = ["A"] * TARGET_LEN
    ref[0] = ref[1] = ref[2] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"   # NGG PAM
    genome = "".join(ref)
    dna = list(genome)
    dna[0] = dna[1] = dna[2] = "M"              # {C,A}
    dna[3] = "W"                                # {A,T}
    split = ["X", GUIDE, "".join(dna), "chrT", "0", "0", "+",
             "1", "0", "1", "NGG", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
    overrides = dict(
        genomeStr=genome, current_chr="chrT", mydict={}, myreg=reg, mygt=None,
        haplotype_check=False, dict_tier_present=False, registry_only_mode=True,
        IUPAC_CAP=2, hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
        pam_begin=-3, pam_end=None, allowed_mms=6, _LOSSLESS_DENSE=lossless,
    )
    return split, overrides


class TestLosslessDense(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.reg = _build_registry([
            (1, "C", "A", "rs1", {"S1": "0|1"}),
            (2, "C", "A", "rs2", {"S1": "0|1"}),
            (3, "C", "A", "rs3", {"S1": "0|1"}),
            (4, "A", "T", "rs4", {"S2": "0|1"}),
        ], self._tmp.name)

    def tearDown(self):
        try:
            self.reg.close()
        except Exception:
            pass
        self._tmp.cleanup()

    def _run(self, lossless):
        split, overrides = _dense_fixture(self.reg, lossless)
        ns = _load_pure_functions(overrides)
        ns["_phase_confirmation_rows"].clear()
        ns["_phase_confirmation_keys"].clear()
        ns["_snp_snp_cooc_rows"].clear()
        ns["_snp_snp_cooc_keys"].clear()
        cluster = []
        ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
        self._ns = ns
        return cluster

    def _has_col3_alt(self, cluster):
        # the union haplotype carries the mm-raising alt 'T' at column 3, which the
        # finalizer lowercases to 't' as a mismatch.
        return any(len(row[2]) > 3 and row[2][3].lower() == "t" for row in cluster)

    def test_flag_off_drops_union(self):
        off = self._run(lossless=False)
        self.assertFalse(self._has_col3_alt(off),
                         "baseline (flag off) must NOT emit the col3-alt union haplotype")

    def test_flag_on_emits_union(self):
        off = self._run(lossless=False)
        on = self._run(lossless=True)
        self.assertGreater(len(on), len(off),
                           "lossless-dense should add the dropped union haplotype")
        self.assertTrue(self._has_col3_alt(on),
                        "flag on must emit the multi-variant union haplotype (col3 alt)")

    def test_snp_snp_cooc_recorded_putative(self):
        # a capped registry-only multi-SNP off-target (the greedy applies alts at cols
        # 0,1,2) must be surfaced in the SNP+SNP co-occurrence companion, PUTATIVE, with
        # >=2 SNP positions and a conservative min-AF bound.
        self._run(lossless=True)
        rows = self._ns["_snp_snp_cooc_rows"]
        self.assertTrue(rows, "a multi-SNP registry-only off-target should record a row")
        self.assertTrue(all(r["Phase"] == "PUTATIVE" for r in rows),
                        "registry-only path must be PUTATIVE")
        multi = [r for r in rows if len(set(
            p for p in str(r["SNP_positions"]).split(",") if p not in ("", "NA", "."))) >= 2]
        self.assertTrue(multi, "at least one row must carry >=2 distinct SNP positions")
        self.assertTrue(all(r["N_carriers"] == "NA" for r in rows),
                        "registry-only path has no countable carriers")
        # min-AF bound is a real float on the recorded rows
        for r in multi:
            self.assertNotEqual(r["MinAF_bound"], ".")
            float(r["MinAF_bound"])


if __name__ == "__main__":
    unittest.main()

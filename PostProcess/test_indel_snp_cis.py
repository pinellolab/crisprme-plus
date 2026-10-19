#!/usr/bin/env python
"""Tests for indel_snp_cis.cis_cooccurrence (the CONFIRMED-cis join).

Run: python PostProcess/test_indel_snp_cis.py   (or via pytest)
"""

import indel_snp_cis as isc


def test_phased_cis_single_snp_confirmed():
    # indel and SNP both on slot 0 -> in cis, CONFIRMED, 1 copy
    cis, phase, ac = isc.cis_cooccurrence({"S": "1|0"}, [{"S": "1|0"}])
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 1


def test_phased_trans_not_cis():
    # indel on slot 0, SNP on slot 1 -> NOT cis; no carrier, still CONFIRMED (no putative)
    cis, phase, ac = isc.cis_cooccurrence({"S": "1|0"}, [{"S": "0|1"}])
    assert cis == set() and ac == 0 and phase == isc.CONFIRMED


def test_homozygous_two_copies():
    cis, phase, ac = isc.cis_cooccurrence({"S": "1|1"}, [{"S": "1|1"}])
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 2


def test_unphased_is_putative():
    cis, phase, ac = isc.cis_cooccurrence({"S": "0/1"}, [{"S": "0/1"}])
    assert cis == {"S"} and phase == isc.PUTATIVE and ac == 1


def test_mixed_phased_and_unphased_forces_putative():
    indel = {"A": "1|0", "B": "0/1"}
    snp = {"A": "1|0", "B": "0/1"}
    cis, phase, ac = isc.cis_cooccurrence(indel, [snp])
    assert cis == {"A", "B"} and phase == isc.PUTATIVE and ac == 2  # A:1 (phased) + B:1 (putative)


def test_carries_indel_but_not_snp_excluded():
    # S carries the indel but is absent from the SNP carriers -> not co-occurring
    cis, phase, ac = isc.cis_cooccurrence({"S": "1|0"}, [{"OTHER": "1|0"}])
    assert cis == set() and ac == 0


def test_multi_snp_all_same_slot_cis():
    indel = {"S": "1|0"}
    snps = [{"S": "1|0"}, {"S": "1|0"}]
    cis, phase, ac = isc.cis_cooccurrence(indel, snps)
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 1


def test_multi_snp_one_trans_breaks_cis():
    indel = {"S": "1|0"}
    snps = [{"S": "1|0"}, {"S": "0|1"}]  # second SNP on the other slot
    cis, phase, ac = isc.cis_cooccurrence(indel, snps)
    assert cis == set() and ac == 0


# --- HIGH-MULTIPLICITY (>=4-variant) cis: the indel-path analog of the SNP path's
# test_N4/test_trans_N4 (test_phased_haplotype.py). The SNP lattice once dropped the
# maximal k>=4 cis haplotype (in-loop-subtraction STARVATION, fixed by the deferred
# peel). This path is structurally immune -- cis_cooccurrence does a DIRECT per-sample
# same-slot check over the indel + ALL used SNP alts (no incremental lattice, nothing
# to starve) -- these tests pin that it scales past 4 variants with no under-report. ---

def test_hi_mult_indel_plus_4_snps_all_cis_confirmed():
    # indel + 4 SNPs (5 variants) all on slot 0 -> CONFIRMED cis, 1 copy. The maximal
    # combo is recovered directly; no k>=4 starvation under-report.
    indel = {"S": "1|0"}
    snps = [{"S": "1|0"}] * 4
    cis, phase, ac = isc.cis_cooccurrence(indel, snps)
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 1


def test_hi_mult_one_of_four_snps_trans_breaks_cis():
    # indel + 4 SNPs, the 3rd on the other slot -> no single slot carries all -> excluded.
    indel = {"S": "1|0"}
    snps = [{"S": "1|0"}, {"S": "1|0"}, {"S": "0|1"}, {"S": "1|0"}]
    cis, phase, ac = isc.cis_cooccurrence(indel, snps)
    assert cis == set() and ac == 0


def test_hi_mult_homozygous_two_cis_copies():
    # indel + 4 SNPs all 1|1 -> both slots carry the full set -> 2 cis copies.
    indel = {"S": "1|1"}
    snps = [{"S": "1|1"}] * 4
    cis, phase, ac = isc.cis_cooccurrence(indel, snps)
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 2


def test_hi_mult_population_mix_cis_trans_unphased():
    # 3 samples x (indel + 3 SNPs = 4 variants):
    #   A: all slot-0 phased          -> CONFIRMED cis (1 copy)
    #   B: one SNP on the other slot  -> not all-cis -> excluded
    #   C: carries all but unphased   -> PUTATIVE cis (1 copy), forces overall PUTATIVE
    indel = {"A": "1|0", "B": "1|0", "C": "0/1"}
    s1 = {"A": "1|0", "B": "1|0", "C": "0/1"}
    s2 = {"A": "1|0", "B": "1|0", "C": "0/1"}
    s3 = {"A": "1|0", "B": "0|1", "C": "0/1"}
    cis, phase, ac = isc.cis_cooccurrence(indel, [s1, s2, s3])
    assert cis == {"A", "C"} and phase == isc.PUTATIVE and ac == 2


def test_no_snps_indel_only_trivially_cis():
    cis, phase, ac = isc.cis_cooccurrence({"S": "1|0", "T": "1|1"}, [])
    assert cis == {"S", "T"} and phase == isc.CONFIRMED and ac == 3  # S:1 + T:2


def test_population_mix_af():
    # 3 samples: A cis(1), B trans(0), C homozygous-cis(2) -> ac=3
    indel = {"A": "1|0", "B": "1|0", "C": "1|1"}
    snp = {"A": "1|0", "B": "0|1", "C": "1|1"}
    cis, phase, ac = isc.cis_cooccurrence(indel, [snp])
    assert cis == {"A", "C"} and phase == isc.CONFIRMED and ac == 3
    assert abs(isc.joint_af(ac, 6954) - 3 / 6954) < 1e-12


def test_putative_joint_af_min_bound():
    # no genotypes: joint cis AF <= the rarest participating allele's AF -> report min
    assert abs(isc.putative_joint_af([0.5, 0.1, 0.3]) - 0.1) < 1e-12
    assert isc.putative_joint_af([0.2]) == 0.2
    assert isc.putative_joint_af([]) == 0.0          # empty -> 0
    assert isc.putative_joint_af([0.5, None]) == 0.0  # any marginal missing -> 0


def _snp_at(table):
    return lambda pos: table.get(pos)


def test_used_snp_guide_matches_alt():
    # guide[1]='C' == alt 'C' (ref 'A') at an IUPAC 'M' target column -> USED
    guide, target = "ACGT", "AMGT"
    snps = {1: ("A", "C", "rs1", {"S": "1|0"})}
    used = isc.used_snps_for_target(guide, target, lambda j: j, _snp_at(snps))
    assert used == [(1, "rs1", {"S": "1|0"})]


def test_ref_satisfied_ambiguity_skipped():
    # guide[1]='A' == ref 'A' -> the reference works, SNP is incidental -> NOT used
    guide, target = "AAGT", "AMGT"
    snps = {1: ("A", "C", "rs1", {"S": "1|0"})}
    used = isc.used_snps_for_target(guide, target, lambda j: j, _snp_at(snps))
    assert used == []


def test_bulge_column_and_missing_snp_skipped():
    guide, target = "ACGT", "AMMT"
    snps = {2: ("A", "G", "rs2", {"S": "1|0"})}  # only col 2 has a SNP
    # col1 IUPAC maps to None (bulge); col2 maps to a real pos with no matching alt
    def o2r(j):
        return None if j == 1 else j
    used = isc.used_snps_for_target(guide, target, o2r, _snp_at(snps))
    # col1 skipped (None); col2 guide 'G' == alt 'G' -> used
    assert used == [(2, "rs2", {"S": "1|0"})]


def test_decode_then_cis_end_to_end():
    # off-target uses one SNP alt; join with the indel genotype -> CONFIRMED cis
    guide, target = "ACGT", "AMGT"
    snps = {1: ("A", "C", "rs1", {"S": "1|0", "T": "0|1"})}
    used = isc.used_snps_for_target(guide, target, lambda j: j, _snp_at(snps))
    snp_gts = [gt for _, _, gt in used]
    indel_gt = {"S": "1|0", "T": "1|0"}
    cis, phase, ac = isc.cis_cooccurrence(indel_gt, snp_gts)
    # S: indel 1|0 + SNP 1|0 -> cis; T: indel 1|0 + SNP 0|1 -> trans (excluded)
    assert cis == {"S"} and phase == isc.CONFIRMED and ac == 1


def test_offset_to_real_plus_strand_left_flank():
    # '+' strand, no bulge, target starts in the left flank (fst<25) -> 1:1 map
    # indel at real start 100, ref/alt len 2/1 (delta +1). target fake_start = indel
    # fake_start + 3 (so fst=3, left flank).
    real = isc.build_offset_to_real("ACGT", fake_start=13, strand="+",
                                    indel_fake_start=10, indel_start=100,
                                    ref="AT", alt="A")
    # fst=3 -> real 103, then 104,105,106 (all left flank, 1:1)
    assert real == [103, 104, 105, 106]


def test_offset_to_real_minus_strand_reverses():
    # '-' strand: leftmost column is the HIGHEST forward fake position
    real = isc.build_offset_to_real("ACGT", fake_start=13, strand="-",
                                    indel_fake_start=10, indel_start=100,
                                    ref="AT", alt="A")
    # 4 bases at forward fake 13,14,15,16 (fst 3,4,5,6 -> real 103..106); '-' reverses
    assert real == [106, 105, 104, 103]


def test_offset_to_real_rna_bulge_column_is_none():
    # a '-' in the DNA (RNA bulge) consumes no fake position
    real = isc.build_offset_to_real("AC-GT", fake_start=13, strand="+",
                                    indel_fake_start=10, indel_start=100,
                                    ref="AT", alt="A")
    assert real[2] is None
    assert real[0] == 103 and real[1] == 104 and real[3] == 105 and real[4] == 106


def test_complement():
    assert isc.complement("A") == "T" and isc.complement("g") == "c"
    assert isc.complement("M") == "M"  # IUPAC/other unchanged


# --- forward_snp_tokens + the CONFIRMED-cis ROW ENRICHMENT (the fix) -------------
# forward_snp_tokens re-reads REF/ALT/rsID/AF from the registry (forward strand),
# NOT from the decode's strand-complemented snp_at. These tests pin the token /
# rsID / MAF format (== the SNP path) and the 1:1 alignment, plus the row-level
# semantics the analyzer applies in the CONFIRMED branch (indel FIRST, carriers
# narrowed to the exact cis set) and that the PUTATIVE/indel-only row is untouched.

class _Cnt:
    def __init__(self, af):
        self._af = af

    def allele_freq(self):
        return self._af


class _FakeReg:
    """Minimal duck-typed Tier-0 registry: {pos1: (alt, ref, rsid_or_None, af)}."""

    def __init__(self, table):
        self._t = table  # pos1 -> (alt, ref, rsid, af)

    def alts_at(self, pos1):
        rec = self._t.get(pos1)
        return [rec[0]] if rec else []

    def ref(self, pos1, alt):
        rec = self._t.get(pos1)
        return rec[1] if rec else None

    def rsid(self, pos1, alt):
        rec = self._t.get(pos1)
        return rec[2] if rec else None

    def lookup(self, pos1, alt):
        rec = self._t.get(pos1)
        return {"global": _Cnt(rec[3])} if rec else None


def test_forward_snp_tokens_format_and_alignment():
    # Two cis SNPs at 0-based real positions 199 (pos1=200) and 299 (pos1=300).
    reg = _FakeReg({200: ("A", "T", "rs9", 0.0123),
                    300: ("T", "C", None, 0.5)})
    used = [(199, "x", {}), (299, "y", {})]  # decode's rsid/gt ignored on purpose
    tok, rs, maf = isc.forward_snp_tokens(reg, "chr18", used)
    assert tok == ["chr18_200_T_A", "chr18_300_C_T"]   # forward-strand chrom_pos_ref_alt
    assert rs == ["rs9", "."]                            # missing rsID -> "."
    assert maf == ["%.6g" % 0.0123, "%.6g" % 0.5]        # SNP-path "%.6g" formatting
    assert len(tok) == len(rs) == len(maf)               # 1:1 aligned


def test_forward_snp_tokens_ignores_target_strand_complement():
    # The token MUST be the forward-strand ref/alt from the registry, regardless of
    # the target strand -- i.e. the decode's complemented alleles must NOT leak in.
    reg = _FakeReg({501: ("G", "C", "rs1", 0.1)})
    used = [(500, "rsZZZ", {})]  # a wrong rsid in `used` must be overridden by the reg
    tok, rs, maf = isc.forward_snp_tokens(reg, "chr1", used)
    assert tok == ["chr1_501_C_G"] and rs == ["rs1"]


def test_forward_snp_tokens_absent_record_skipped_keeps_alignment():
    reg = _FakeReg({200: ("A", "T", "rs9", 0.01)})  # pos1=300 absent
    used = [(199, "x", {}), (299, "y", {})]
    tok, rs, maf = isc.forward_snp_tokens(reg, "chr2", used)
    # the absent SNP drops from ALL three lists together (no misalignment)
    assert tok == ["chr2_200_T_A"] and rs == ["rs9"] and len(maf) == 1


def _enrich_confirmed_row(final_result, indel_data, reg, chrom, used, cis, phase):
    """Replicate the analyzer's CONFIRMED-cis enrichment (the exact mutation the fix
    applies) so the row-level contract is testable without driving the full loop.
    Mirrors analisi_indels_NNN.py: gated on phase == CONFIRMED; indel token FIRST."""
    if phase != isc.CONFIRMED:
        return  # PUTATIVE / non-confirmed -> row stays indel-only (sidecar-only)
    snp_tok, snp_rs, snp_maf = isc.forward_snp_tokens(reg, chrom, used)
    if snp_tok:
        final_result[17] = indel_data[4] + "," + ",".join(snp_tok)
        final_result[15] = ",".join([indel_data[2]] + snp_rs)
        final_result[16] = ",".join([indel_data[3]] + snp_maf)
        final_result[12] = ",".join(sorted(cis))


def test_confirmed_row_enrichment_full_haplotype():
    # indel_data = [fakechr, samples, rsID, MAF, INDELinfo, fakestart]
    indel_data = ["chr18_x", "A,B,C,HGzz", "rs_indel", "0.002",
                  "chr18_12386282_TTTTC_T", 0]
    final_result = [""] * 18
    final_result[12] = indel_data[1]   # indel-all-carriers (to be REPLACED)
    final_result[15] = indel_data[2]
    final_result[16] = indel_data[3]
    final_result[17] = indel_data[4]
    reg = _FakeReg({12386297: ("A", "T", "rs_snp1", 0.3),
                    12386300: ("T", "TTTTC", "rs_snp2", 0.1)})
    used = [(12386296, "x", {}), (12386299, "y", {})]  # 0-based real -> pos1 +1
    cis = {"B", "A"}  # exact cis carriers (narrower than indel-all-carriers)
    _enrich_confirmed_row(final_result, indel_data, reg, "chr18", used, cis,
                          isc.CONFIRMED)
    # col-17: indel token FIRST, then forward-strand SNP tokens
    assert final_result[17] == (
        "chr18_12386282_TTTTC_T,chr18_12386297_T_A,chr18_12386300_TTTTC_T")
    # carriers narrowed to the exact cis set (sorted), replacing indel-all-carriers
    assert final_result[12] == "A,B"
    # rsID + MAF stay 1:1 aligned with col-17 (indel first)
    assert final_result[15] == "rs_indel,rs_snp1,rs_snp2"
    assert final_result[16] == "0.002," + ("%.6g" % 0.3) + "," + ("%.6g" % 0.1)
    assert (len(final_result[17].split(",")) == len(final_result[15].split(","))
            == len(final_result[16].split(",")))


def test_putative_row_stays_indel_only():
    # a PUTATIVE (unphased / no-genotype) call must NOT mutate the row columns
    indel_data = ["chr1_x", "A,B", "rs_i", "0.01", "chr1_100_AT_A", 0]
    final_result = [""] * 18
    final_result[12], final_result[15] = indel_data[1], indel_data[2]
    final_result[16], final_result[17] = indel_data[3], indel_data[4]
    reg = _FakeReg({201: ("C", "A", "rs1", 0.2)})
    used = [(200, "x", {})]
    _enrich_confirmed_row(final_result, indel_data, reg, "chr1", used,
                          {"A", "B"}, isc.PUTATIVE)
    assert final_result[17] == "chr1_100_AT_A"   # indel-only, unchanged
    assert final_result[12] == "A,B"             # carriers unchanged
    assert final_result[15] == "rs_i" and final_result[16] == "0.01"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

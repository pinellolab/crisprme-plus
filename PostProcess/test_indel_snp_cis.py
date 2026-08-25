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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

#!/usr/bin/env python
"""Tests for overlay_indel_snps: the fake->real piecewise map + IUPAC overlay.

Run: python PostProcess/test_indel_snp_overlay.py   (or via pytest)
"""

import overlay_indel_snps as ov


def test_map_deletion_AT_to_A():
    # Real chr22 row: indel chr22_10510355_AT_A, CHR start = 10510329 (= POS-26).
    sp = 10510329
    ref, alt = "AT", "A"
    rl, al = len(ref), len(alt)  # 2, 1 -> delta = +1
    # left flank 1:1
    assert ov.map_fake_offset_to_real(0, sp, rl, al) == sp
    assert ov.map_fake_offset_to_real(24, sp, rl, al) == sp + 24
    # ALT anchor base -> no reference position
    assert ov.map_fake_offset_to_real(25, sp, rl, al) is None
    # downstream shifted by +delta (=+1): fake offset 26 -> real sp+27
    assert ov.map_fake_offset_to_real(26, sp, rl, al) == sp + 27
    assert ov.map_fake_offset_to_real(52, sp, rl, al) == sp + 53


def test_map_insertion_G_to_GATAA():
    # Real chr22 row: indel chr22_10511378_GATAA_G is a DELETION in that row, but
    # exercise the symmetric INSERTION G>GATAA here (ref_len=1, alt_len=5 -> delta=-4).
    sp = 500
    ref, alt = "G", "GATAA"
    rl, al = len(ref), len(alt)
    assert ov.map_fake_offset_to_real(0, sp, rl, al) == sp
    assert ov.map_fake_offset_to_real(24, sp, rl, al) == sp + 24
    # the 5 ALT bases (anchor G + inserted ATAA) have no reference position
    for k in range(25, 30):
        assert ov.map_fake_offset_to_real(k, sp, rl, al) is None
    # downstream: fake offset 30 -> real sp+30+(1-5) = sp+26
    assert ov.map_fake_offset_to_real(30, sp, rl, al) == sp + 26
    assert ov.map_fake_offset_to_real(31, sp, rl, al) == sp + 27


def test_map_snp_like_ref1_alt1():
    # Degenerate case len(ref)==len(alt)==1 (delta 0): fully 1:1 except the single
    # ALT base at offset 25.
    sp = 0
    assert ov.map_fake_offset_to_real(25, sp, 1, 1) is None
    assert ov.map_fake_offset_to_real(26, sp, 1, 1) == 26  # delta 0
    assert ov.map_fake_offset_to_real(10, sp, 1, 1) == 10


def test_overlay_substitutes_flank_iupac_only():
    sp = 100
    ref, alt = "AT", "A"  # delta +1
    # fake sub: 25 'A' left flank + anchor 'A' (offset 25) + 10 'C' downstream
    sub = "A" * 25 + "A" + "C" * 10  # length 36
    enriched = ["A"] * 300
    enriched[105] = "M"   # left-flank k=5  (real sp+5 = 105)
    enriched[127] = "R"   # downstream k=26 (real sp+26+1 = 127)
    enriched[100 + 25] = "Y"  # real position of the ALT anchor -> must NOT overlay (k=25 is None)
    out = ov.overlay_sub_fasta(sub, sp, ref, alt, "".join(enriched))
    assert len(out) == len(sub)                 # length preserved
    assert out[5] == "M"                        # flank SNP overlaid
    assert out[26] == "R"                       # downstream SNP overlaid
    assert out[25] == "A"                       # ALT base untouched (k=25 -> None)
    # everything else unchanged
    for i, ch in enumerate(out):
        if i not in (5, 26):
            assert ch == sub[i]


def test_overlay_ignores_non_iupac_and_N():
    sp = 0
    sub = "N" + "A" * 24 + "A" + "G" * 10   # a leading N pad in the flank
    enriched = ["A"] * 100
    enriched[0] = "M"   # real 0 maps to fake k=0, but fake base is 'N' -> skip
    enriched[3] = "S"   # real 3, fake k=3 is 'A' -> overlay
    out = ov.overlay_sub_fasta(sub, sp, "A", "A", "".join(enriched))
    assert out[0] == "N"    # N pad never overlaid
    assert out[3] == "S"    # plain flank base overlaid
    assert len(out) == len(sub)


def test_overlay_out_of_bounds_real_position():
    # downstream real position past the enriched chromosome end must be skipped, not crash
    sp = 90
    sub = "A" * 25 + "A" + "C" * 10
    enriched = "A" * 95  # shorter than sp+downstream real positions
    out = ov.overlay_sub_fasta(sub, sp, "AT", "A", enriched)
    assert out == sub  # nothing overlaid (all downstream reals >= 95), no exception


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

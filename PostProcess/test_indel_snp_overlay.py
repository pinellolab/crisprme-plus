#!/usr/bin/env python
"""Tests for overlay_indel_snps: the fake->real piecewise map + IUPAC overlay.

The map and the overlay take each window's left flank from its own log row (POS - 1 -
start, indel_window.left_flank): 25 bases for an enrichment of CRISPRitz 2.9.1 or earlier
(the windows of every index published up to CRISPRme 2.8.1, which must still read right),
30 for one of 2.10.0.

Run: python PostProcess/test_indel_snp_overlay.py   (or via pytest)
"""

import unittest

import overlay_indel_snps as ov


def test_map_deletion_AT_to_A():
    # Real chr22 row of CRISPRitz 2.9.1: indel chr22_10510355_AT_A, CHR start = 10510329
    # (= POS-26: a left flank of 25).
    sp = 10510329
    ref, alt = "AT", "A"
    rl, al = len(ref), len(alt)  # 2, 1 -> delta = +1
    # left flank 1:1
    assert ov.map_fake_offset_to_real(0, sp, rl, al, 25) == sp
    assert ov.map_fake_offset_to_real(24, sp, rl, al, 25) == sp + 24
    # ALT anchor base -> no reference position
    assert ov.map_fake_offset_to_real(25, sp, rl, al, 25) is None
    # downstream shifted by +delta (=+1): fake offset 26 -> real sp+27
    assert ov.map_fake_offset_to_real(26, sp, rl, al, 25) == sp + 27
    assert ov.map_fake_offset_to_real(52, sp, rl, al, 25) == sp + 53


def test_map_insertion_G_to_GATAA():
    # Real chr22 row: indel chr22_10511378_GATAA_G is a DELETION in that row, but
    # exercise the symmetric INSERTION G>GATAA here (ref_len=1, alt_len=5 -> delta=-4).
    sp = 500
    ref, alt = "G", "GATAA"
    rl, al = len(ref), len(alt)
    assert ov.map_fake_offset_to_real(0, sp, rl, al, 25) == sp
    assert ov.map_fake_offset_to_real(24, sp, rl, al, 25) == sp + 24
    # the 5 ALT bases (anchor G + inserted ATAA) have no reference position
    for k in range(25, 30):
        assert ov.map_fake_offset_to_real(k, sp, rl, al, 25) is None
    # downstream: fake offset 30 -> real sp+30+(1-5) = sp+26
    assert ov.map_fake_offset_to_real(30, sp, rl, al, 25) == sp + 26
    assert ov.map_fake_offset_to_real(31, sp, rl, al, 25) == sp + 27


def test_map_snp_like_ref1_alt1():
    # Degenerate case len(ref)==len(alt)==1 (delta 0): fully 1:1 except the single
    # ALT base at offset 25.
    sp = 0
    assert ov.map_fake_offset_to_real(25, sp, 1, 1, 25) is None
    assert ov.map_fake_offset_to_real(26, sp, 1, 1, 25) == 26  # delta 0
    assert ov.map_fake_offset_to_real(10, sp, 1, 1, 25) == 10


def test_overlay_substitutes_flank_iupac_only():
    sp = 100
    ref, alt = "AT", "A"  # delta +1
    # fake sub: 25 'A' left flank + anchor 'A' (offset 25) + 10 'C' downstream
    sub = "A" * 25 + "A" + "C" * 10  # length 36
    enriched = ["A"] * 300
    enriched[105] = "M"   # left-flank k=5  (real sp+5 = 105)
    enriched[127] = "R"   # downstream k=26 (real sp+26+1 = 127)
    enriched[100 + 25] = "Y"  # real position of the ALT anchor -> must NOT overlay (k=25 is None)
    out = ov.overlay_sub_fasta(sub, sp, ref, alt, 25, "".join(enriched))
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
    out = ov.overlay_sub_fasta(sub, sp, "A", "A", 25, "".join(enriched))
    assert out[0] == "N"    # N pad never overlaid
    assert out[3] == "S"    # plain flank base overlaid
    assert len(out) == len(sub)


def test_overlay_out_of_bounds_real_position():
    # downstream real position past the enriched chromosome end must be skipped, not crash
    sp = 90
    sub = "A" * 25 + "A" + "C" * 10
    enriched = "A" * 95  # shorter than sp+downstream real positions
    out = ov.overlay_sub_fasta(sub, sp, "AT", "A", 25, enriched)
    assert out == sub  # nothing overlaid (all downstream reals >= 95), no exception


def test_crispr_bulge_never_sees_an_iupac_code():
    """The target CRISPR-Bulge scores is the per-haplotype resolution of the overlaid
    one (indel_snp_cis.resolve_target), and the model's agct2numerals KeyErrors on an
    IUPAC code: no code may survive into a resolved row. A target with no code (the
    classic, non-overlaid fake genome) is never resolved at all."""
    import indel_snp_cis as isc

    assert isc.ambiguous_columns("ACGTACGT-ACGT") == []  # classic: untouched
    # a realistic overlaid minus-strand window: K at a spacer column, R in the PAM
    tgt = "-TAAKCAGTTGCTATTTTAAAGRG"
    gd = "CTAA-CAGTTGCTTTTATCACNNN"
    amb = isc.ambiguous_columns(tgt)
    assert amb == [4, 22]
    cols = [isc.SnpColumn(4, 100, "A", {"C": {"S1": "1|0"}}, "K", "-"),
            isc.SnpColumn(22, 83, "C", {"T": {"S2": "0|1"}}, "R", "-")]
    rows = isc.resolve_target(
        gd, tgt, cols, {"S1": "1|0", "S2": "1|1", "S3": "0|1"},
        isc.pam_columns(len(tgt), "NRG"), range(len(tgt) - 3), 6, -1, 2)
    assert rows
    for r in rows:
        assert set(r["seq"]) <= set("ACGTacgt-"), r["seq"]
    # and so for the worst case of a sites-only panel (no genotypes)
    cols = [isc.SnpColumn(c.col, c.real, c.ref_fwd, None, tgt[c.col], "-",
                          known=False) for c in cols]
    rows = isc.resolve_target(gd, tgt, cols, None, isc.pam_columns(len(tgt), "NRG"),
                              range(len(tgt) - 3), 6, -1, 2)
    assert len(rows) == 1 and set(rows[0]["seq"]) <= set("ACGTacgt-"), rows


# --- the windows of CRISPRitz 2.10.0: 30 reference bases on either side, read from the
# log row, as are 2.9.1's 25 / 27 ----------------------------------------------------------

def _fake_genome(tmp, ref, records, left, right):
    """A fake-indel genome and its log, cut from ``ref`` as CRISPRitz's enricher cuts them
    with a window of ``left`` / ``right`` reference bases: [(pos1, REF, ALT)] -> (fake
    FASTA, log, [(fake start, fake end)])."""
    import os

    flat, rows, spans = "", [], []
    for k, (pos1, r, a) in enumerate(records):
        start = pos1 - 1 - left
        assert ref[pos1 - 1:pos1 - 1 + len(r)] == r
        win = ref[start:pos1 - 1] + a + ref[pos1 - 1 + len(r):pos1 - 1 + len(r) + right]
        spans.append((len(flat), len(flat) + len(win)))
        rows.append("chrT_%d-%d_%d\tS1\trs%d\t0.1\tchrT_%d_%s_%s\t%d,%d\t%s" % (
            start, pos1 - 1 + len(r) + right, k + 1, k, pos1, r, a, len(flat),
            len(flat) + len(win), ref[start:start + len(win)]))
        flat += win + "N"
    fa = os.path.join(tmp, "fakechrT.fa")
    with open(fa, "w") as fh:
        fh.write(">fakechrT\n" + "".join(flat[a:b] + "\nN\n" for a, b in spans))
    log = os.path.join(tmp, "logchrT.txt")
    with open(log, "w") as fh:
        fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n" + "\n".join(rows) + "\n")
    return fa, log, spans


def test_map_a_window_of_30():
    # chr22_10510355_AT_A in a window of CRISPRitz 2.10.0: CHR start 10510324 (= POS-31)
    sp = 10510324
    assert ov.map_fake_offset_to_real(0, sp, 2, 1, 30) == sp
    assert ov.map_fake_offset_to_real(29, sp, 2, 1, 30) == sp + 29   # the last flank base
    assert ov.map_fake_offset_to_real(30, sp, 2, 1, 30) is None      # the ALT's anchor
    assert ov.map_fake_offset_to_real(31, sp, 2, 1, 30) == sp + 32   # past the deletion
    # the same row read with 2.9.1's 25-base flank puts flank bases on the ALT and off by
    # the deletion: what a fixed constant did to a window of 30
    assert ov.map_fake_offset_to_real(27, sp, 2, 1, 25) == sp + 28 != sp + 27


def test_parse_log_indels_takes_each_rows_left_flank():
    import tempfile

    ref = "".join("ACGT"[(i * 7 + i // 5) % 4] for i in range(400))
    with tempfile.TemporaryDirectory() as tmp:
        for left, right in ((25, 27), (30, 30)):
            fa, log, spans = _fake_genome(tmp, ref, [(100, ref[99], ref[99] + "TT"),
                                                     (200, ref[199:202], ref[199])],
                                          left, right)
            got = list(ov.parse_log_indels(log))
            assert [g[5] for g in got] == [left, left], (left, got)
            assert [g[0] for g in got] == [100 - 1 - left, 200 - 1 - left]
            assert [(g[3], g[4]) for g in got] == spans


def test_overlay_of_a_fake_genome_of_either_window():
    """overlay_fake_chromosome codes every flank base whose reference position carries a
    SNP code in the enriched genome, and no ALT base, on the windows of CRISPRitz 2.9.1
    (25 / 27) and 2.10.0 (30 / 30) alike: each window read with its own left flank. With
    the 25-base flank of 2.9.1 taken for every window, a window of 30 gets codes on its
    inserted bases and none on the five flank bases before its REF."""
    import os
    import tempfile

    ref = "".join("ACGT"[(i * 7 + i // 5) % 4] for i in range(400))
    code = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
    for left, right in ((25, 27), (30, 30)):
        with tempfile.TemporaryDirectory() as tmp:
            records = [(100, ref[99], ref[99] + "TTTTT"), (200, ref[199:202], ref[199])]
            fa, log, spans = _fake_genome(tmp, ref, records, left, right)
            # a SNP code at every reference position: each flank base must get its own
            enriched = "".join(code[b] for b in ref)
            en = os.path.join(tmp, "chrT.enriched.fa")
            with open(en, "w") as fh:
                fh.write(">chrT\n" + enriched + "\n")
            ov.overlay_fake_chromosome(fa, en, log)
            with open(fa) as fh:
                fh.readline()
                flat = fh.read().replace("\n", "")
            for (pos1, r, a), (f0, f1) in zip(records, spans):
                win = flat[f0:f1]
                start = pos1 - 1 - left
                want = (enriched[start:pos1 - 1] + a
                        + enriched[pos1 - 1 + len(r):pos1 - 1 + len(r) + right])
                assert win == want, (left, pos1, win, want)



def load_tests(loader, tests, pattern):
    """Under ``python -m unittest`` too, each test_ function above, which it found none of,
    and said OK (suite P4)."""
    return unittest.TestSuite(unittest.FunctionTestCase(fn) for name, fn in
                              sorted(globals().items()) if name.startswith("test_")
                              and callable(fn))

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

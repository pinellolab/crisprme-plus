#!/usr/bin/env python
"""Tests for indel_snp_cis.cis_cooccurrence (the CONFIRMED-cis join).

Run: python PostProcess/test_indel_snp_cis.py   (or via pytest)
"""

import unittest

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


def load_tests(loader, tests, pattern):
    """Under ``python -m unittest`` too, each test_ function above, which it found none of,
    and said OK (suite P4)."""
    return unittest.TestSuite(unittest.FunctionTestCase(fn) for name, fn in
                              sorted(globals().items()) if name.startswith("test_")
                              and callable(fn))


# --- where an indel row starts (reference_start). A window starts at W (0-based, the
# log's real start); its REF is at offset 25, the indel's POS (1-based) at W + 26.
W = 1000


def _start(offset, ref, alt):
    return isc.reference_start(offset, W + 26, W, ref, alt)


def test_reference_start_flank_and_anchor_are_their_own_place():
    for ref, alt in (("G", "GCAATCATTAT"), ("AGCAGCT", "A"), ("AC", "TG")):
        assert _start(3, ref, alt) == W + 3      # left flank
        assert _start(25, ref, alt) == W + 25    # the anchor (an MNV's first base)


def test_reference_start_on_an_added_base_is_the_insertions_anchor():
    # chr22:22322762 G>GCAATCATTAT (HGDP1kGP): a target starting at the 7th added base
    # was placed 6 bases past the POS, at its haplotype place, before 2.10.0
    for k in range(1, 11):
        assert _start(25 + k, "G", "GCAATCATTAT") == W + 25
    assert _start(25 + 11, "G", "GCAATCATTAT") == W + 26     # past it: the reference


def test_reference_start_deletion_and_mnv_unchanged():
    assert _start(26, "AGCAGCT", "A") == W + 26 + 6           # past a deletion
    assert _start(26, "AC", "TG") == W + 26                   # an MNV, base for base
    assert _start(26, "ATTT", "AG") == W + 26                 # a base standing for REF's


def test_reference_start_complex_and_unnormalized_insertions():
    # AT>ACCC: the first C stands in T's place, the others follow it
    assert [_start(25 + k, "AT", "ACCC") for k in range(5)] == [
        W + 25, W + 26, W + 26, W + 26, W + 27]
    # ATG>ATCCG (a shared tail G): the CC follow T, the tail is the reference's G
    assert [_start(25 + k, "ATG", "ATCCG") for k in range(6)] == [
        W + 25, W + 26, W + 26, W + 26, W + 27, W + 28]


def test_reference_start_agrees_with_offset_to_real_on_reference_bases():
    # wherever the window's base has a reference position, the row starts there
    for ref, alt in (("G", "GCAATCATTAT"), ("AGCAGCT", "A"), ("AC", "TG")):
        n = 25 + len(alt) + 27
        for off in range(n):
            real = isc.map_fake_offset_to_real(off, W, len(ref), len(alt))
            if real is not None:
                assert _start(off, ref, alt) == real, (ref, alt, off)


# --- counting along an indel row's target (variant_places, 2.10.0): a variant's place on
# the haplotype carrying the indel, and the row's first base there.


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_haplotype_place_before_and_past_the_indel():
    # an insertion of 4 at POS 101 (anchor 0-based 100): a base past it is 4 on
    assert isc.haplotype_place(90, 101, "T", "TCACC") == 90
    assert isc.haplotype_place(101, 101, "T", "TCACC") == 101      # its own anchor
    assert isc.haplotype_place(102, 101, "T", "TCACC") == 106
    # a deletion of 3: a base past it is 3 back, a deleted one is on no haplotype
    assert isc.haplotype_place(105, 101, "TCAC", "T") == 102
    assert isc.haplotype_place(102, 101, "TCAC", "T") is None
    # an MNV: base for base; an unnormalized record's shared tail moves with the change
    assert isc.haplotype_place(102, 101, "AC", "TG") is None
    assert isc.haplotype_place(103, 101, "AC", "TG") == 103
    assert isc.haplotype_place(104, 101, "ATG", "ATCCG") == 106
    assert isc.haplotype_place(103, 101, "ATG", "ATCCG") == 105    # the tail G


def test_haplotype_base_kind():
    # T>TCACC at POS 101: the anchor (100) is a reference base, 101..104 added, then the
    # reference from 101 on
    kinds = [isc.haplotype_base_kind(h, 101, "T", "TCACC") for h in range(99, 107)]
    assert kinds == [("ref", 99), ("ref", 100), ("added", None), ("added", None),
                     ("added", None), ("added", None), ("ref", 101), ("ref", 102)]
    # AT>ACCC: the first C stands in T's place, the other two are added
    assert [isc.haplotype_base_kind(h, 101, "AT", "ACCC")[0] for h in range(100, 105)] == [
        "ref", "alt", "added", "added", "ref"]
    # a deletion and an MNV add nothing
    assert isc.haplotype_base_kind(101, 101, "TCAC", "T") == ("ref", 104)
    assert isc.haplotype_base_kind(101, 101, "AC", "TG") == ("alt", 101)


def test_row_first_base_on_an_inserted_base():
    """The oracle's Z3 geometry: TAAG + the guide's first six bases inserted after the
    anchor T (0-based 99); the target (the guide + TGG) starts at the 4th added base, so
    the row starts at the anchor and its first base is 4 on, at 103."""
    g = "ACTGAAATCTGTAAGCAGGC"
    ins = "AAG" + g[:6]
    ref = "N" * 99 + "T" + g[6:12] + "C" + g[13:] + "TGG" + "GGGGGG"
    snp = (100 + 7, "C", g[12])  # the reference C at 0-based 106, the guide's base in cis
    assert ref[106] == "C" and g[12] != "C"
    alt = g + "TGG"
    ref_aligned = ref[99:99 + 23].lower()
    toks = [(100, "T", "T" + ins), snp]
    first, places = isc.variant_places("+", 99, alt, ref_aligned, toks)
    assert first == 103 and places == [100, 107 + 9]
    assert [max(1, p - first) for p in places] == [1, 13]   # the SNP is the 13th base
    # with no SNP named the first base does not matter: the reference's
    assert isc.variant_places("+", 99, alt, ref_aligned, toks[:1]) == (99, [100])


def test_row_first_base_minus_strand_chr22():
    """chr22:21052639 (-), HGDP + 1000 Genomes, README guide: T>TCACC at 21052640 and
    G>A at 21052641 (in LD); the target's leftmost base is the 3rd added base."""
    toks = [(21052640, "T", "TCACC"), (21052641, "G", "A")]
    first, places = isc.variant_places(
        "-", 21052639, "C-AACAGTGaGCaTgTgTCACTGG", "C-cAacagTgagcaTgtgtcaCCA", toks)
    assert first == 21052642 and places == [21052640, 21052645]
    # resultIntegrator's - strand number: len(aligned spacer) + 1 - distance
    assert abs(places[1] - first - 1 - len("CTAACAGT-TGCTTTTATCACNNN")) == 22


def test_row_first_base_reference_start_and_ambiguity():
    # a row starting on the anchor itself: the anchor's own place
    toks = [(101, "T", "TCA"), (103, "G", "C")]
    ref = "N" * 100 + "T" + "AGGTTACCAGTACCAGTTAGGACCTA"
    hap = "T" + "CA" + "AC" + ref[103:]     # the anchor, CA added, A, the SNP's C ...
    alt = hap[:23]
    assert isc.row_first_base("+", 100, alt, ref[100:123], toks[0], toks[1:]) == 100
    # the 2nd added base: one on
    assert isc.row_first_base("+", 100, hap[2:25], ref[100:123], toks[0], toks[1:]) == 102
    # a run the target lies wholly in reads the same bases from several places: the
    # smallest is taken (the anchor's own base here)
    run = "N" * 100 + "A" * 40
    assert isc.row_first_base("+", 100, "A" * 23, run[100:123], (101, "A", "AAA"),
                              [(110, "A", "A")]) == 100
    # bases the row's target is not: its start
    assert isc.row_first_base("+", 100, "G" * 23, ref[100:123], toks[0], toks[1:]) == 100
    # rows that start before the anchor, deletions, MNVs: the start
    assert isc.row_first_base("+", 90, alt, ref[90:113], toks[0], toks[1:]) == 90
    assert isc.row_first_base("+", 100, alt, ref[100:123], (101, "TAG", "T"), ()) == 100
    assert isc.row_first_base("+", 100, alt, ref[100:123], (101, "TA", "GC"), ()) == 100


def test_variant_places_snps_past_a_deletion_and_snp_rows():
    # a SNP past a 3-base deletion is 3 back; one before it keeps its place
    toks = [(101, "TCAC", "T"), (95, "A", "G"), (110, "C", "T")]
    assert isc.variant_places("+", 90, "X", "X", toks) == (90, [101, 95, 107])
    # SNP rows and MNV rows: the reference's
    assert isc.variant_places("+", 90, "X", "X", [(95, "A", "G"), (99, "C", "T")]) == (
        90, [95, 99])
    assert isc.variant_places("-", 90, "X", "X", [(101, "AC", "TG"), (110, "C", "T")]) == (
        90, [101, 110])


def _row_first_base_base_by_base(strand, start, alt_aligned, ref_aligned, indel, snps=()):
    """row_first_base as c3a75fe wrote it: the haplotype built one base at a time
    (haplotype_base_kind), then a slice compared at each k. The model the sliced one must
    give back."""
    ipos, iref, ialt = indel
    if len(ialt) <= len(iref) or start != ipos - 2 + len(iref) - isc.shared_tail(iref, ialt):
        return start
    target = alt_aligned.replace("-", "").upper()
    refb = ref_aligned.replace("-", "").upper()
    if not target or len(refb) != len(target):
        return start
    if strand == "-":
        target, refb = _rc(target), _rc(refb)
    n, more = len(target), len(ialt) - len(iref)
    hap = []
    for h in range(start, start + n + more):
        kind, p = isc.haplotype_base_kind(h, ipos, iref, ialt)
        if kind == "ref":
            hap.append(refb[p - start] if 0 <= p - start < n else "?")
        else:
            hap.append(ialt[h - ipos + 1].upper())
    for pos1, _r, alt in snps:
        at = isc.haplotype_place(pos1, ipos, iref, ialt)
        if at is not None and 0 <= at - 1 - start < len(hap):
            hap[at - 1 - start] = alt.upper()
    hap = "".join(hap)
    for k in range(more + 1):
        if hap[k:k + n] == target:
            return start + k
    return start


def test_row_first_base_equals_the_base_by_base_reading():
    """2.10.0: row_first_base cuts the haplotype in slices and finds the target with
    str.find. On 4,000 random rows (plain, unnormalized, complex and anchorless insertions,
    deletions and MNVs; rows at the anchor, on it and near it; both strands; gaps and lower
    case; SNPs in and out of the window, a multi-base one now and then; low-complexity
    sequence, where several k hold the target; targets that are nowhere) it gives what the
    base-by-base reading gives."""
    import random
    rng = random.Random(282)

    def rnd(n, alpha):
        return "".join(rng.choice(alpha) for _ in range(n))

    placed = 0
    for _ in range(4000):
        alpha = rng.choice(["ACGT", "AC", "A", "ACGT"])
        ref = rnd(260, alpha)
        ipos = rng.randint(100, 140)
        a = ipos - 1
        shape = rng.randrange(6)
        if shape == 0:    # plain insertion
            iref = ref[a]; ialt = iref + rnd(rng.randint(1, 60), alpha)
        elif shape == 1:  # unnormalized, a shared tail
            iref = ref[a:a + 3]; ialt = iref[:2] + rnd(rng.randint(1, 30), alpha) + iref[2:]
        elif shape == 2:  # complex: a replaced base, then added ones
            iref = ref[a:a + 2]; ialt = iref[0] + rnd(rng.randint(2, 40), alpha)
        elif shape == 3:  # no first base shared
            iref = ref[a]; ialt = rnd(rng.randint(1, 20), alpha) + iref
        elif shape == 4:  # several first bases shared
            iref = ref[a:a + 4]; ialt = iref + rnd(rng.randint(1, 50), alpha)
        else:             # a deletion or an MNV: the start
            iref = ref[a:a + 3]; ialt = rng.choice([iref[0], rnd(3, "ACGT")])
        hap = ref[:a] + ialt + ref[a + len(iref):]
        anchor = ipos - 2 + len(iref) - isc.shared_tail(iref, ialt)
        start = anchor if rng.random() < 0.85 else anchor + rng.choice([-3, -1, 1, 2])
        n = rng.choice([20, 23, 24])
        k = rng.randint(0, max(0, len(ialt) - len(iref) + 2))
        snps = []
        for _s in range(rng.choice([0, 0, 1, 2, 3])):
            p1 = rng.randint(a - 5, a + len(iref) + n + 10)
            snps.append((p1, ref[p1 - 1], rng.choice("ACGT") if rng.random() < 0.9
                         else rnd(2, "ACGT")))
        target = hap[start + k:start + k + n] if rng.random() < 0.9 else rnd(n, "ACGT")
        ref_text = ref[start:start + n]
        if len(target) != n:
            continue
        if rng.random() < 0.3:  # an RNA bulge: a gap at one place in both
            g = rng.randint(1, n - 2)
            target, ref_text = target[:g] + "-" + target[g:], ref_text[:g] + "-" + ref_text[g:]
        target = "".join(c.lower() if rng.random() < 0.2 else c for c in target)
        strand = rng.choice("+-")
        if strand == "-":
            target, ref_text = _rc(target.upper()), _rc(ref_text)
        args = (strand, start, target, ref_text, (ipos, iref, ialt), snps)
        want = _row_first_base_base_by_base(*args)
        assert isc.row_first_base(*args) == want, args
        placed += want != start
    assert placed > 800  # about a third of the rows start on an added base (1,241)


def test_row_first_base_cost_does_not_grow_with_the_insertion():
    """2.10.0: a row at the anchor of a long insertion is placed without a Python step per
    inserted base. HPRC holds 402 insertions of 53.9 to 54.4 kb at chr22:10690626, and built
    base by base (as up to c3a75fe) a row there took 83 ms: the merge's 16,058 rows of the
    README guide inside them took most of its 27 minutes. The lines of indel_snp_cis run for
    a row 20 bases from the end of a 1 kb and of a 20 kb insertion are as many, and a 1 Mb
    one is placed in well under a second (process time; base by base, about one second)."""
    import sys
    import time

    def run_counting_lines(*args):
        count = [0]

        def local(frame, event, arg):
            if event == "line":
                count[0] += 1
            return local

        def tracer(frame, event, arg):
            return local if frame.f_globals is vars(isc) else None

        before = sys.gettrace()
        sys.settrace(tracer)
        try:
            out = isc.row_first_base(*args)
        finally:
            sys.settrace(before)
        return out, count[0]

    ref = ("GATTACAGCT" * 30)[:300]
    a = 100  # the anchor (0-based), VCF POS 101

    def row(length, strand):
        ins = "C" * (length - 30) + "ATGCATGCAATTGGCCAAGTTCAGTACGTA"
        alt = ref[a] + ins
        hap = ref[:a + 1] + ins + ref[a + 1:]
        k = length - 20  # its first base: the k-th added base
        target, ref_text = hap[a + k:a + k + 23], ref[a:a + 23]
        if strand == "-":
            target, ref_text = _rc(target), _rc(ref_text)
        return (strand, a, target, ref_text, (a + 1, ref[a], alt), [(a + 8, ref[a + 7], "T")]), k

    lines = {}
    for length in (1000, 20000):
        for strand in "+-":
            args, k = row(length, strand)
            first, lines[length, strand] = run_counting_lines(*args)
            assert first == a + k, (length, strand, first)
    assert lines[1000, "+"] == lines[20000, "+"] and lines[1000, "-"] == lines[20000, "-"], lines
    t0 = time.process_time()
    for strand in "+-":
        args, k = row(1000000, strand)
        assert isc.row_first_base(*args) == a + k
    assert time.process_time() - t0 < 0.5


# --- where an indel row's PAM is judged on the reference (pam_reference_place, 2.10.0):
# the reference bases its PAM sits on in a placement of the indel, through its repeat, that
# leaves it on the reference; else, a PAM made wholly of inserted bases in every placement,
# created (D4); else the reference PAM beside the protospacer's base that touches the PAM,
# placed on the reference; else created. ------------------------------------------------

def _trims(ref, alt):
    c = 0
    while c < min(len(ref), len(alt)) and ref[c].upper() == alt[c].upper():
        c += 1
    t = 0
    while t < min(len(ref), len(alt)) - c and ref[-1 - t].upper() == alt[-1 - t].upper():
        t += 1
    return c, t


def _pam_place_every_placement(genome, indel, snps, first, n, k, pam_right):
    """pam_reference_place's rule read the long way: every place of the indel that gives
    the same haplotype (the reference with the row's indel and SNP alts) is built whole,
    each of its bases mapped to the reference, and the rule's four steps are read off
    those maps: the VCF's placement first, then the ones right of it, then left. A complex
    record's first min(len(REF), len(ALT)) changed bases are written into the reference
    as the SNP alts are, in the places of the REF bases they replace (whose reference
    bases are read), and its other ALT bases added after them, or its other REF bases
    deleted; the haplotype so built is the record's (checked)."""
    ipos, iref, ialt = indel
    c, t = _trims(iref, ialt)
    rm, am = iref[c:len(iref) - t], ialt[c:len(ialt) - t].upper()
    if len(rm) == len(am):
        return None
    rs = list(genome.upper())
    for p, _r, a in snps:
        rs[p - 1] = a.upper()
    whole = "".join(rs[:ipos - 1]) + ialt.upper() + "".join(rs[ipos - 1 + len(iref):])
    r = min(len(rm), len(am))
    for i in range(r):  # a complex record's bases in REF bases' places
        rs[ipos - 1 + c + i] = am[i]
    rs = "".join(rs)
    rm, am = rm[r:], am[r:]
    a0, d, m = ipos - 1 + c + r, len(rm), len(am)

    def build(at, ins):
        if m:
            return (rs[:at] + ins + rs[at:],
                    list(range(at)) + [None] * m + list(range(at, len(rs))))
        return rs[:at] + rs[at + d:], list(range(at)) + list(range(at + d, len(rs)))

    h0 = build(a0, am)[0]
    assert h0 == whole, (indel, snps)
    maps = []
    span = m + d + n + k + 8
    for at in sorted(range(a0 - span, a0 + span + 1),
                     key=lambda x: (x != a0, x < a0, abs(x - a0))):
        if at < 0 or at + d > len(rs):
            continue
        hj, pos = build(at, h0[at:at + m])
        if hj == h0:
            maps.append(pos)
    p0 = first + n - k if pam_right else first
    tb = p0 - 1 if pam_right else p0 + k
    free = set()
    for pos in maps:
        refs = pos[p0:p0 + k]
        if None not in refs and refs == list(range(refs[0], refs[0] + k)):
            free.add(genome[refs[0]:refs[0] + k].upper())
    if free:
        return "ref", frozenset(free)
    # D4 (2.10.0): a PAM every placement makes wholly of the inserted bases is created,
    # whatever PAM the reference holds beside the protospacer
    if m and all(all(x is None for x in pos[p0:p0 + k]) for pos in maps):
        return "created", frozenset()
    for pos in maps:
        if pos[tb] is not None:
            lo = pos[tb] + 1 if pam_right else pos[tb] - k
            return "ref", frozenset([genome[lo:lo + k].upper()])
    return "created", frozenset()


def _bases(genome, k, place):
    """pam_reference_place's answer with each place read: the reference bases there."""
    if place is None:
        return None
    return place[0], frozenset(genome[lo:lo + k].upper() for lo in place[1])


def _random_pam_rows(seed, count, complex_share=0.0):
    """Random indel rows: a reference holding a tandem repeat (a unit of 1 to 6 bases, 1 to
    12 times) or none, an insertion (of whole units, a rotated or a broken unit, or random
    bases) or a deletion in or beside it, written left-aligned as bcftools norm writes it
    or at another place of its repeat, with its anchor, without one, or with a shared tail;
    0 to 2 SNP alts near it; a target of 20 to 24 bases holding the indel's changed bases
    or its junction, with a PAM of 3 or 4 bases at either end. With ``complex_share``, that
    share of the rows have a complex record instead: the 1 or 2 REF bases from the anchor
    on replaced (the first by another base) and the same bases added, or deleted. Yields
    (genome, indel, snps, first, n, k, pam_right)."""
    import random
    rng = random.Random(seed)

    def rnd(size, alpha="ACGT"):
        return "".join(rng.choice(alpha) for _ in range(size))

    made = 0
    while made < count:
        unit = rnd(rng.randint(1, 6))
        reps = rng.randint(0, 12)
        left = rnd(rng.randint(60, 90))
        genome = left + unit * reps + rnd(rng.randint(60, 90))
        lo = len(left)
        a = rng.randint(lo - 3, lo + len(unit) * reps + 2)  # where the change starts
        shape = rng.randrange(7)
        if shape <= 1:    # whole units, rotated to where it is written
            r = (a - lo) % len(unit)
            ins = (unit[r:] + unit[:r]) * rng.randint(1, 3)
        elif shape == 2:  # a unit with one base changed
            ins = list(unit * rng.randint(1, 2))
            ins[rng.randrange(len(ins))] = rng.choice("ACGT")
            ins = "".join(ins)
        elif shape == 3:  # random bases
            ins = rnd(rng.randint(1, 30))
        if shape <= 3:
            kind = rng.randrange(3)
            if kind == 0:    # with its anchor, the VCF's way
                ipos, iref, ialt = a, genome[a - 1], genome[a - 1] + ins
            elif kind == 1:  # no first base shared
                ipos, iref, ialt = a + 1, genome[a], ins + genome[a]
            else:            # unnormalized: a shared tail
                ipos, iref, ialt = a, genome[a - 1:a + 1], genome[a - 1] + ins + genome[a]
        else:             # a deletion: whole units or random bases
            dl = len(unit) * rng.randint(1, 2) if shape <= 5 else rng.randint(1, 8)
            ipos, iref, ialt = a, genome[a - 1:a + dl], genome[a - 1]
        if complex_share and rng.random() < complex_share:
            # a complex record: the j REF bases from the anchor on replaced
            j = rng.choice([1, 1, 2])
            put = rng.choice([b for b in "ACGT" if b != genome[a - 1]]) + rnd(j - 1)
            if shape <= 3:
                ipos, iref, ialt = a, genome[a - 1:a - 1 + j], put + ins
            else:
                ipos, iref, ialt = a, genome[a - 1:a - 1 + j + dl], put
        c, t = _trims(iref, ialt)
        gone, added = len(iref) - c - t, len(ialt) - c - t
        if gone == added or (not complex_share and (gone > 0) == (added > 0)):
            continue
        if rng.random() < 0.7:  # left-aligned, as bcftools norm writes it
            while ipos > 2 and len(iref) == 1 + gone and len(ialt) == 1 + added and \
                    iref[-1].upper() == ialt[-1].upper() and genome[ipos - 2] == iref[-1]:
                ipos -= 1
                iref, ialt = genome[ipos - 1] + iref[:-1], genome[ipos - 1] + ialt[:-1]
        e = ipos - 1 + c
        n, k = rng.randint(20, 24), rng.choice([3, 4])
        # the target holds an added base, or both sides of the deleted ones
        first = rng.randint(e - n + 1, e + max(added, 1) - 1) if added \
            else rng.randint(e - n + 1, e - 1)
        if first < 0 or first + n + 30 > len(genome) - gone + added:
            continue
        snps = []
        for _s in range(rng.choice([0, 0, 0, 1, 2])):
            q = rng.randint(max(1, first - 10), first + n + 10)  # 1-based
            if ipos - 1 + c < q <= ipos - 1 + c + gone or ipos <= q <= ipos - 1 + len(iref):
                continue  # a base the indel deletes or holds
            snps.append((q, genome[q - 1], rng.choice("ACGT")))
        made += 1
        yield genome, (ipos, iref, ialt), snps, first, n, k, rng.random() < 0.5


def test_pam_reference_place_equals_every_placement_read():
    """2.10.0: pam_reference_place reads only the shifts its tests need, each one slice; on
    3,000 random rows (insertions and deletions of repeat units, broken units and random
    bases, left-aligned or not, with and without an anchor, with a shared tail; SNP alts
    near the PAM; PAMs of 3 and 4 bases at either end) it gives what building every
    placement of the indel and reading its rule off them gives."""
    shifted = created = 0
    for row in _random_pam_rows(2100, 3000):
        want = _pam_place_every_placement(*row)
        got = _bases(row[0], row[5], isc.pam_reference_place(*row))
        assert got == want, (row, isc.pam_reference_place(*row), want)
        created += got == ("created", frozenset())
        genome, (ipos, iref, ialt), snps, first, n, k, right = row
        c, t = _trims(iref, ialt)
        e, added = ipos - 1 + c, len(ialt) - len(iref)
        p0 = first + n - k if right else first
        # the VCF's placement holds an added base in the PAM, or the junction inside it
        if (added > 0 and e < p0 + k and e + added > p0) or (added < 0 and p0 < e < p0 + k):
            shifted += 1
    assert shifted > 600 and created > 150, (shifted, created)


def test_pam_reference_place_does_not_depend_on_where_the_vcf_writes_the_indel():
    """2.10.0: the same haplotype written with its indel at every place of its repeat (the
    VCF left-aligns it; a record left unnormalized, or a caller's other choice, writes it
    further right) gives one answer: the same reference bases decide, or none."""
    def haplotype(genome, p, gone, ins, snps):  # the reference with the SNP alts, then
        rs = list(genome)                         # the record at POS p written in
        for q, _r, a in snps:
            rs[q - 1] = a
        rs = "".join(rs)
        return rs[:p] + ins + rs[p + gone:]

    checked = 0
    for genome, (ipos, iref, ialt), snps, first, n, k, right in _random_pam_rows(2101, 1500):
        c, t = _trims(iref, ialt)
        if c != 1 or t:
            continue
        gone, added = len(iref) - 1, len(ialt) - 1
        hap = haplotype(genome, ipos, gone, ialt[1:], snps)
        answers = set()
        for p in range(ipos - 40, ipos + 41):  # every record that writes this haplotype,
            if p < 1 or p + gone > len(genome) - 1:  # its SNP alts where they are
                continue
            if any(p < q <= p + gone for q, _r, _a in snps):
                continue
            r = genome[p - 1:p + gone]
            a = genome[p - 1] + (hap[p:p + added] if added else "")
            if haplotype(genome, p, gone, a[1:], snps) != hap:
                continue
            answers.add(_bases(genome, k, isc.pam_reference_place(genome, (p, r, a), snps,
                                                                  first, n, k, right)))
        assert len(answers) == 1, (genome, ipos, iref, ialt, snps, first, n, k, right, answers)
        checked += 1
    assert checked > 500, checked


def test_pam_reference_place_geometries():
    """2.10.0, by hand: + strand, 3' PAM (NGG at the target's right end) unless said."""
    g = "ACTGAAATCTGTAAGCAGGC"
    # a tandem duplication of the PAM's bases: the VCF writes the copy at the protospacer's
    # end, a placement 3 bases right puts the PAM on the reference's own TGG
    genome = "T" * 30 + g + "TGG" + "CATTACG" + "T" * 30
    assert isc.pam_reference_place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == \
        ("ref", (50,))
    # TGG before the reference's AGG, no repeat: the PAM wholly inserted wherever the
    # insertion sits, an additional PAM beside a protospacer that had a valid one, AGG:
    # created (D4; before it, the reference PAM beside the protospacer decided, AGG: NA)
    genome = "T" * 30 + g + "AGG" + "CATTACG" + "T" * 30
    assert isc.pam_reference_place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == \
        ("created", ())
    # the insertion holds the protospacer's last two bases and the PAM's first wherever it
    # is placed (it can move one base right, no further): created
    genome = "T" * 30 + g[:18] + "GG" + "AGG" + "T" * 30
    assert isc.pam_reference_place(genome, (48, "G", "GGCT"), [], 30, 23, 3, True) == \
        ("created", ())
    # - strand, 3' PAM (CCN forward at the target's left end): C>CCA after the anchor C,
    # no repeat: the protospacer's touching base is a reference one, the PAM beside it CTC
    # (GAG on -)
    genome = "T" * 27 + "CTC" + _rc(g) + "TTAA" + "T" * 30
    assert isc.pam_reference_place(genome, (30, "C", "CCA"), [], 29, 23, 3, False) == \
        ("ref", (27,))
    # the same where the insertion CTGC repeats the four bases after the anchor: placed 4
    # bases right, the PAM is reference bases, CCT (AGG on -)
    genome = "T" * 28 + "GAC" + _rc(g)[2:] + "GGAA" + "T" * 30
    assert isc.pam_reference_place(genome, (31, "C", "CCTGC"), [], 30, 23, 3, False) == \
        ("ref", (30,))
    # a deletion whose junction is inside the PAM: on the - strand the PAM beside the
    # protospacer is the reference bases before the protospacer's first, past the deleted
    # ones (the REF alignment, read from the row's start, reads the deleted ones)
    genome = "T" * 30 + "C" + "GTAC" + "CA" + _rc(g) + "T" * 30
    assert isc.pam_reference_place(genome, (31, "CGTAC", "C"), [], 30, 23, 3, False) == \
        ("ref", (34,))
    # an MNV: not this rule's (the caller reads the reference under its PAM)
    assert isc.pam_reference_place(genome, (31, "CG", "TA"), [], 30, 23, 3, False) is None
    # a complex record is: CG>TAA is CG replaced by TA and an A added after them; the PAM
    # (TAA forward) holds the added A, and the protospacer's touching base is the reference
    # one after it: the PAM beside it is the reference's TCG, the REF bases CG included
    assert isc.pam_reference_place(genome, (31, "CG", "TAA"), [], 30, 23, 3, False) == \
        ("ref", (29,))
    # Cas12a (TTTV at the target's left end, + strand): T>TTT, the reference TTTA beside
    # the spacer
    genome = "C" * 28 + "TTTA" + "GATCCAGTCACGTTAGCATG" + "C" * 30
    assert isc.pam_reference_place(genome, (31, "T", "TTT"), [], 30, 24, 4, False) == \
        ("ref", (28,))


def test_pam_reference_place_an_additional_inserted_pam_is_created():
    """2.10.0, decision D4: an additional PAM, made wholly of the bases an insertion adds
    wherever the insertion sits, beside a reference protospacer is created, even where the
    reference already holds a valid PAM beside that protospacer (the haplotype has a PAM the
    reference does not have there). Before D4 the reference PAM beside the protospacer
    decided, and a valid one there read NA. Unchanged by D4: a placement that puts the PAM
    on the reference's own bases (a tandem duplication of the PAM: NA), a PAM that only some
    placements make wholly of inserted bases or that holds a reference base (the reference
    PAM beside the protospacer decides), and a deletion's junction."""
    g = "ACTGAAATCTGTAAGCAGGC"
    place = isc.pam_reference_place
    # + strand, 3' PAM: TGG inserted after the protospacer, beside the reference's AGG (a
    # valid NGG) and beside CAT (none): created both ways
    for after in ("AGG", "CAT"):
        genome = "T" * 30 + g + after + "CATTACG" + "T" * 30
        assert place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == ("created", ()), after
    # beside TCA the inserted TGG can sit one base right (GGT after the reference's T), where
    # the PAM holds that reference T: not every placement makes it wholly of inserted bases,
    # and the reference PAM beside the protospacer decides, TCA (no NGG: created by the
    # caller)
    genome = "T" * 30 + g + "TCA" + "CATTACG" + "T" * 30
    assert place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == ("ref", (50,))
    # the tandem duplication of the reference's own TGG: a placement 3 bases right puts the
    # PAM on it (NA, the exception D4 keeps)
    genome = "T" * 30 + g + "TGG" + "CATTACG" + "T" * 30
    assert place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == ("ref", (50,))
    # TGG inserted beside the reference's TGA: placed one or two bases right (GGT, GTG) the
    # PAM holds a reference base, so not every placement makes it wholly of inserted bases:
    # the reference PAM beside the protospacer decides, TGA (no NGG: created by the caller)
    genome = "T" * 30 + g + "TGA" + "CATTACG" + "T" * 30
    assert place(genome, (50, "C", "CTGG"), [], 30, 23, 3, True) == ("ref", (50,))
    # a PAM of two inserted bases and a reference one (TG, then the reference's G): the
    # reference PAM beside the protospacer decides, GTA (no NGG: created by the caller)
    genome = "T" * 30 + g + "GTA" + "CATTACG" + "T" * 30
    assert place(genome, (50, "C", "CTG"), [], 30, 23, 3, True) == ("ref", (50,))
    # - strand, 3' PAM (forward CCN at the target's left end): TTCCA inserted after the
    # anchor G, a target from its 3rd added base, the PAM CCA (TGG on -) the last three
    # added bases, beside the reference's CCG (CGG on -, a valid NGG): created (the oracle's
    # QW, TRAC chr17:77419671 (-))
    rcg = _rc(g)
    genome = "T" * 30 + "CCG" + rcg + "T" * 30
    assert place(genome, (33, "G", "GTTCCA"), [], 35, 23, 3, False) == ("created", ())
    # Cas12a, 5' PAM, - strand (the PAM at the target's right end, forward): TAAA (TTTA on
    # -) inserted after the spacer's last forward base, beside the reference's GAAA (TTTC
    # on -, a valid TTTV): created; beside its own TAAA, a tandem duplication: NA
    sp = "GATCCAGTCACGTTAGCATGCAA"
    f = _rc(sp)
    for after, want in (("GAAA", ("created", ())), ("TAAA", ("ref", (53,)))):
        genome = "C" * 30 + f + after + "CCGT" + "C" * 30
        got = place(genome, (53, f[-1], f[-1] + "TAAA"), [], 30, 27, 4, True)
        assert got == want, (after, got)
    # + strand, 5' PAM (at the target's left end): TTTA inserted before the spacer, beside
    # the reference's TTTC (valid): created
    genome = "G" * 26 + "CA" + "TTTC" + sp + "G" * 30
    assert place(genome, (32, "C", "CTTTA"), [], 32, 27, 4, False) == ("created", ())
    # a complex record's added bases are an insertion's: A>CTGG on the protospacer's last
    # base, its added TGG beside the reference's AGG: created
    genome = "T" * 40 + g[:19] + "A" + "AGG" + "CATTAC" + "T" * 30
    assert place(genome, (60, "A", g[19] + "TGG"), [], 40, 23, 3, True) == ("created", ())
    # a deletion is not an additional PAM: CA>C, its junction inside the PAM CCG (CGG on
    # -), keeps the reference PAM beside the protospacer, CAG (CTG on -)
    genome = "T" * 30 + "CCAG" + rcg + "T" * 30
    assert place(genome, (32, "CA", "C"), [], 30, 23, 3, False) == ("ref", (31,))


def test_pam_reference_place_reads_a_complex_record_as_replaced_and_added_bases():
    """2.10.0: a complex record (REF bases replaced and bases added or lost, as HPRC writes
    many insertions: A>T plus inserted bases) is its replaced bases, in the places of the
    REF bases they replace and read as those (a SNP alt's rule), and its other bases, added
    after them (or deleted), placed as an insertion's (a deletion's). On 3,000 random rows,
    half of them complex records, pam_reference_place gives what building every placement
    of those bases, and reading the rule off them, gives. On 50803b6 a complex
    record was not this rule's: a PAM holding its added bases was read on the REF
    alignment's window, at the anchor on a row that starts inside them (NA where that
    window happens to hold a valid PAM, though the PAM has no place on the reference)."""
    complex_rows = created = beside = 0
    for row in _random_pam_rows(2102, 3000, complex_share=0.5):
        want = _pam_place_every_placement(*row)
        got = _bases(row[0], row[5], isc.pam_reference_place(*row))
        assert got == want, (row, isc.pam_reference_place(*row), want)
        genome, (ipos, iref, ialt), snps, first, n, k, right = row
        c, t = _trims(iref, ialt)
        if len(iref) - c - t and len(ialt) - c - t:
            complex_rows += 1
            created += got == ("created", frozenset())
            beside += got[0] == "ref" and len(iref) < len(ialt)
    assert complex_rows > 900 and created > 80 and beside > 300, \
        (complex_rows, created, beside)


def test_pam_reference_place_complex_records_by_hand():
    """2.10.0, complex records by hand (3' PAM; pamc-review-r1's F1 and its residual)."""
    g = "ACTGAAATCTGTAAGCAGGC"
    # + strand: A>T and 33 bases added after it, which hold a whole target: its PAM and the
    # protospacer's end are added bases wherever they sit: created, as for the insertion of
    # the same bases after the A (A>A+I); the REF alignment's window, read from the anchor,
    # holds AGG
    ins = "CCGTA" + g + "TGG" + "CATCC"
    genome = "C" * 40 + "A" + "GTCAGTCAGTCAGTCAGTC" + "AGG" + "C" * 40
    for alt in ("T" + ins, "A" + ins):
        assert isc.pam_reference_place(genome, (41, "A", alt), [], 46, 23, 3, True) == \
            ("created", ())
    # - strand (forward CCA at the target's left end): C>T and the mirror bases added; the
    # reference after the anchor C holds C (the window, CCN, is NGG on -)
    ins = "GTCAG" + _rc(g + "TGG") + "GATCA"
    genome = "T" * 40 + "C" + "C" + "GTAGTAGTAGTAGTAGTAGTAGTA" + "T" * 40
    for alt in ("T" + ins, "C" + ins):
        assert isc.pam_reference_place(genome, (41, "C", alt), [], 46, 23, 3, False) == \
            ("created", ())
    # - strand, G>CCA: the row starts on the replaced base, the PAM CCA (TGG on -) is it
    # and the two added bases, and the protospacer's touching base is a reference one: the
    # reference PAM beside it, CTG (CAG on -), holds the motif (the window, read from the
    # row's start, GGC, is GCC on -)
    genome = "T" * 38 + "CT" + "G" + _rc(g) + "TTAA" + "T" * 30
    assert isc.pam_reference_place(genome, (41, "G", "CCA"), [], 40, 23, 3, False) == \
        ("ref", (38,))
    # - strand, TG>C: T replaced by C and G deleted, the junction inside the PAM CCA (TGG
    # on -): the reference PAM beside the protospacer is TGA (TCA on -), no motif (the
    # window, read from the row's start, CTG, is CAG on -, with the deleted G)
    genome = "T" * 38 + "CTGA" + _rc(g) + "T" * 30
    assert isc.pam_reference_place(genome, (40, "TG", "C"), [], 38, 23, 3, False) == \
        ("ref", (39,))
    # + strand, A>GTTT on the PAM's last base: the PAM is reference bases and the replaced
    # one, read as the reference's A (TGA: created by the record); the added TTT follow it
    genome = "T" * 40 + g + "TGA" + "CATTAC" + "T" * 30
    assert isc.pam_reference_place(genome, (63, "A", "GTTT"), [], 40, 23, 3, True) == \
        ("ref", (60,))
    # + strand, A>CTGG on the protospacer's last base: the PAM is the added TGG wherever it
    # sits, beside the protospacer whose touching base is the replaced one, in its REF base's
    # place, where the reference holds AGG: created, as for the insertion of TGG (D4; before
    # it, the reference PAM beside it decided, AGG: NA)
    genome = "T" * 40 + g[:19] + "A" + "AGG" + "CATTAC" + "T" * 30
    assert isc.pam_reference_place(genome, (60, "A", g[19] + "TGG"), [], 40, 23, 3, True) \
        == ("created", ())


def test_pam_reference_place_reads_a_snp_alt_in_the_pam_as_its_reference_base():
    """2.10.0: the placements are read on the reference with the row's SNP alts, and a SNP
    alt in a PAM that a placement puts on the reference is read as the reference base it
    replaces (the SNP path's rule). + strand: TGG inserted after the protospacer, where the
    reference holds TGA; the row's A>G on that A makes the next bases TGG too, so the
    insertion can be placed past the PAM, whose place on the reference holds TGA: created
    either way (without the SNP the insertion stays, and the PAM beside the protospacer is
    that TGA)."""
    g = "ACTGAAATCTGTAAGCAGGC"
    genome = "T" * 30 + g + "TGA" + "CATTACG" + "T" * 30
    for snps in ([(53, "A", "G")], []):
        assert isc.pam_reference_place(genome, (50, "C", "CTGG"), snps, 30, 23, 3, True) \
            == ("ref", (50,))
    # a row whose PAM the VCF's placement leaves on the reference (the insertion before it,
    # in the protospacer), with a SNP alt in it: that place; and where the insertion can
    # also be placed past the PAM (a G run), both places, which can differ there
    genome = "T" * 30 + "C" + "G" * 25 + "ATTACG" + "T" * 30
    assert isc.pam_reference_place(genome, (31, "C", "CG"), [], 35, 23, 3, True) == \
        ("ref", (54,))
    snps = [(51, "G", "A")]  # stops the run: the insertion cannot pass it
    assert isc.pam_reference_place(genome, (31, "C", "CG"), snps, 35, 23, 3, True) == \
        ("ref", (54,))


def test_pam_reference_place_cost_does_not_grow_with_the_insertion():
    """2.10.0: a target inside a long insertion is judged with a few slice comparisons, not
    one Python step per inserted base: the lines of indel_snp_cis run for a target 30 bases
    from the end of a 1 kb and a 20 kb insertion are as many, in a repeat (a tandem
    duplication, which the placements walk through) and out of one, and a 1 Mb one takes
    well under a second (process time); so for a complex record whose 1 kb, 20 kb or 1 Mb
    of replaced bases come before its added ones."""
    import sys
    import time

    def count_lines(*args):
        count = [0]

        def local(frame, event, arg):
            if event == "line":
                count[0] += 1
            return local

        def tracer(frame, event, arg):
            return local if frame.f_globals is vars(isc) else None

        before = sys.gettrace()
        sys.settrace(tracer)
        try:
            out = isc.pam_reference_place(*args)
        finally:
            sys.settrace(before)
        return out, count[0]

    def row(length, dup):
        unit = ("ACGTTGCAAC" * (length // 10 + 1))[:length]
        tail = unit if dup else "T" + unit[1:]  # a copy of the inserted bases, or not
        genome = "G" * 200 + "C" + tail + "G" * 200
        first = 200 + 1 + length - 30
        return (genome, (201, "C", "C" + unit), [], first, 23, 3, True)

    def complex_row(length, dup):  # ``length`` REF bases replaced, then 100 bases added
        unit = ("ACGTTGCAAC" * 10)[:100]
        tail = unit if dup else "T" + unit[1:]
        genome = "G" * 200 + "A" * length + tail + "G" * 200
        first = 200 + length + 100 - 30
        return (genome, (201, genome[200:200 + length], "C" * length + unit), [], first, 23,
                3, True)

    lines = {}
    for dup in (True, False):
        for length in (1000, 20000):
            got, lines[dup, length] = count_lines(*row(length, dup))
            assert got == (("ref", (200 + 1 + length - 30 + 20,)) if dup else ("created", ()))
            got, lines["cx", dup, length] = count_lines(*complex_row(length, dup))
            assert got == (("ref", (200 + length + 90,)) if dup else ("created", ()))
        assert lines[dup, 1000] == lines[dup, 20000], lines
        assert lines["cx", dup, 1000] == lines["cx", dup, 20000], lines
    t0 = time.process_time()
    for dup in (True, False):
        isc.pam_reference_place(*row(1000000, dup))
        isc.pam_reference_place(*complex_row(1000000, dup))
    assert time.process_time() - t0 < 0.5


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")

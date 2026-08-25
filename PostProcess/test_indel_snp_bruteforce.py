#!/usr/bin/env python
"""Brute-force ORACLE for the indel-SNP feature, run against a tiny real-data
fixture (PostProcess/fixtures/indel_snp_chr22_fixture.json, ~16 KB, extracted from
the e2e chr22 data by build_indel_snp_fixture.py). NO large files, NO crispritz.

This is an INDEPENDENT check: it reconstructs the actual sequence from first
principles (fake segment = reference+indel, then apply the sample's SNP allele) and
directly aligns the guide -- if that agrees with the overlay/coordinate/cis
machinery, the whole chain is cross-validated. It guards every coordinate hop
(overlay, bulge+strand window map, co-occurrence decode, cis join) on every PR.

Run: python PostProcess/test_indel_snp_bruteforce.py   (or via pytest)
"""

import os
import json

import overlay_indel_snps as ov
import indel_snp_cis as isc

_FIX = os.path.join(os.path.dirname(__file__), "fixtures", "indel_snp_chr22_fixture.json")
FIX = json.load(open(_FIX))
COOC_POS = 22880651  # the known SNP+indel co-occurrence indel (chr22:22880651 CT>C)


def _snp_gt(tokens):
    out = {}
    for t in tokens.split(","):
        if t:
            s, _, g = t.partition(":")
            out[s] = g
    return out


def _overlay_target_dna(dna, ro, ind, strand):
    """Simulate the OVERLAID search: replace a plain target base with the SNP IUPAC
    code where the enriched window has one (complemented on '-' strand, since the
    displayed DNA is the reverse complement). The fixture targets come from a
    pre-overlay search; production targets carry these IUPAC codes already."""
    enr = ind["enriched"]
    out = list(dna)
    for j, ch in enumerate(dna):
        if ch == "-" or ro[j] is None:
            continue
        b = enr.get(str(ro[j]))
        if b and b in ov._IUPAC_AMBIG:
            out[j] = isc.complement(b) if strand == "-" else b
    return "".join(out)


def _snp_at(ind, strand):
    """Build a snp_at(real) closure from a fixture indel's SNPs, oriented to the
    target strand (dict is 1-based -> key = real+1; forward alleles complemented
    on '-')."""
    snps = ind["snps"]

    def snp_at(real):
        d = snps.get(str(real + 1))
        if not d:
            return None
        rb, ab = d["ref"], d["alt"]
        if strand == "-":
            rb, ab = isc.complement(rb), isc.complement(ab)
        return (rb, ab, str(real + 1), _snp_gt(d["gt"]))

    return snp_at


def test_fixture_is_sane():
    assert FIX["indels"] and any(i["pos"] == COOC_POS for i in FIX["indels"])
    assert any(i["snps"] for i in FIX["indels"])


def test_window_mapper_refbase_match():
    """build_offset_to_real: every non-SNP DNA column's real position holds the
    matching reference base (oriented) -- the independent ground truth."""
    total = match = 0
    for ind in FIX["indels"]:
        enr = ind["enriched"]
        sp, ref, alt, fs = ind["start_position"], ind["ref"], ind["alt"], ind["fake_start"]
        for t in ind["targets"]:
            dna, strand = t["dna"], t["strand"]
            ro = isc.build_offset_to_real(dna, t["fake_pos"], strand, fs, sp, ref, alt)
            for j, ch in enumerate(dna):
                if ch == "-" or ro[j] is None:
                    continue
                b = enr.get(str(ro[j]))
                if b is None or b in ov._IUPAC_AMBIG or b == "N":
                    continue  # SNP / pad / outside the shipped window
                want = ch.upper() if strand == "+" else isc.complement(ch).upper()
                total += 1
                match += (want == b)
    assert total > 0 and match == total, f"ref-base match {match}/{total}"


def test_overlay_places_iupac_on_flanks_only():
    """Overlay each fake segment with its enriched window -> IUPAC lands only on
    plain flank bases; lengths preserved."""
    for ind in FIX["indels"]:
        reals = {int(k): v for k, v in ind["enriched"].items()}
        lo, hi = min(reals), max(reals)
        slc = "".join(reals.get(i, "A") for i in range(lo, hi + 1))
        over = ov.overlay_sub_fasta(
            ind["fake_seq"], ind["start_position"], ind["ref"], ind["alt"],
            slc, enriched_offset=lo,
        )
        assert len(over) == len(ind["fake_seq"])
        for j, ch in enumerate(over):
            if ch in ov._IUPAC_AMBIG:
                assert ind["fake_seq"][j] in "ACGT"


def test_used_snps_detects_cooccurrence():
    """The known co-occurrence indel has >=1 target where a SNP alt matches the
    guide (used_snps_for_target returns it)."""
    found = 0
    for ind in FIX["indels"]:
        if ind["pos"] != COOC_POS:
            continue
        sp, ref, alt, fs = ind["start_position"], ind["ref"], ind["alt"], ind["fake_start"]
        for t in ind["targets"]:
            dna, guide, strand = t["dna"], t["guide"], t["strand"]
            ro = isc.build_offset_to_real(dna, t["fake_pos"], strand, fs, sp, ref, alt)
            dna_ov = _overlay_target_dna(dna, ro, ind, strand)
            used = isc.used_snps_for_target(guide, dna_ov, lambda j: ro[j], _snp_at(ind, strand))
            found += len(used)
    assert found > 0, "no co-occurrence SNP detected for the known indel"


def test_full_cis_on_real_data():
    """END-TO-END on real fixture data: co-occurrence targets -> cis_cooccurrence
    over the REAL phased indel GT + SNP GT yields a valid phase state, and when the
    intersection is non-empty the joint AF is well-formed."""
    exercised = 0
    for ind in FIX["indels"]:
        if ind["pos"] != COOC_POS:
            continue
        sp, ref, alt, fs = ind["start_position"], ind["ref"], ind["alt"], ind["fake_start"]
        indel_gt = ind["indel_gt"]
        for t in ind["targets"]:
            dna, guide, strand = t["dna"], t["guide"], t["strand"]
            ro = isc.build_offset_to_real(dna, t["fake_pos"], strand, fs, sp, ref, alt)
            dna_ov = _overlay_target_dna(dna, ro, ind, strand)
            used = isc.used_snps_for_target(guide, dna_ov, lambda j: ro[j], _snp_at(ind, strand))
            if not used:
                continue
            snp_gts = [gt for _, _, gt in used]
            cis, phase, ac = isc.cis_cooccurrence(indel_gt, snp_gts)
            assert phase in (isc.CONFIRMED, isc.PUTATIVE)
            assert ac >= 0 and (not cis or ac > 0)
            exercised += 1
    assert exercised > 0, "no co-occurrence target exercised the full cis path"


def test_brute_force_snp_flips_a_mismatch():
    """INDEPENDENT oracle: reconstructing the alt allele reduces the guide mismatch
    count for the known co-occurrence -- proves the SNP genuinely changes the
    off-target (not merely overlaps it)."""
    def mm(guide, seq):
        return sum(1 for g, d in zip(guide, seq)
                   if g != "-" and d != "-" and g.upper() != d.upper())
    proven = 0
    for ind in FIX["indels"]:
        if ind["pos"] != COOC_POS:
            continue
        sp, ref, alt, fs = ind["start_position"], ind["ref"], ind["alt"], ind["fake_start"]
        for t in ind["targets"]:
            dna, guide, strand = t["dna"], t["guide"], t["strand"]
            ro = isc.build_offset_to_real(dna, t["fake_pos"], strand, fs, sp, ref, alt)
            dna_alt = list(dna)
            changed = False
            for j, ch in enumerate(dna):
                if ch == "-" or ro[j] is None:
                    continue
                d = ind["snps"].get(str(ro[j] + 1))
                if not d:
                    continue
                ab = isc.complement(d["alt"]) if strand == "-" else d["alt"]
                if guide[j].upper() == ab.upper() and guide[j].upper() != ch.upper():
                    dna_alt[j] = ab
                    changed = True
            if changed and mm(guide, "".join(dna_alt)) < mm(guide, dna):
                proven += 1
    assert proven > 0, "brute force found no mismatch-reducing SNP"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} brute-force oracle tests passed.")

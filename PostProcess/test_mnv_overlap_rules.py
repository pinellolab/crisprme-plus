#!/usr/bin/env python
"""The 2.10.1 readings of the record shapes HPRC v2.1's decomposition leaves behind.

Three shapes put two records of the phased store, or a record and a SNP, on one phased
slot in a way no real haplotype has. They are artefacts of decomposing the pangenome
graph: `vcfwave` cuts each bubble into pieces, HPRC's merge joins pieces of different
bubbles at one POS, and our own preparation (`bcftools norm -m -any -f hg38`) left-aligns
each indel against hg38 alone, ignoring the haplotype's other variants, so in a repeat it
slides across them. The reading of each was measured against HPRC's pre-vcfwave original
records (`hprc-v2.1-mc-grch38.raw.vcf.gz`, read with vcfbub's rule), 921 slots in all,
16 of them confirmed character for character against the carriers' release-2 assemblies.

  1. A complex record whose REF lies wholly inside the bases a deletion deletes (31
     pairs, 249 slots). The slot leaves both records' rows and is listed nowhere, as up
     to 2.10.0. Reading the complex record as changing nothing, so that the slot stays on
     the deletion's rows, would list 4,864 23-mers of which 692 are false (48 slots);
     outside the chr1:977,063 VNTR that reading is 676 real against 692 false.
  2. An SNV and an MNV on one slot (28 pairs, 150 slots): the slot is read base by base.
     The MNV's bases stand where it replaces them; on a base it keeps, the SNV's base
     wins (the originals decide 16 of the 20 such slots, and name the SNV's base in all
     16: chr6 x14, chr9, chr19, two of them confirmed against assemblies). An SNV whose
     base is neither the MNV's base there nor the reference base it keeps is a true
     third-base conflict: one slot of HPRC v2.1 (chrY:21,604,661 AAC>CAG with
     21,604,663 C>T, one chrY copy, the records' CONFLICT field naming another sample,
     the assembly reading CAGATACACC), and that slot leaves both records' rows.
  3. An indel at a complex record's POS, one of the two keeping no anchor base (50
     pairs, 522 slots): the slot is read with the complex record alone, which is 84.8 %
     real (98.0 % for the insertions) and catches 90.4 % of the real targets. Reading
     each record without the other is 48.5 % real with 12,009 false 23-mers, and the
     indel alone is 23.4 % real.

Every "real" above means the 23-mer occurs somewhere in the carrier haplotype's true
sequence over the window, so in a tandem repeat a wrong reading can still score no false
listing: each false count is a lower bound and each "% real" an upper bound. Of the false
listings left by the reading of shape 3, 79 % come from one locus, chr1:108,452,980,
where no reading of our records is right.

The fixtures are those shapes: HPRC v2.1 records at their own coordinates, on the
reference bases of hg38 around them, with the carriers written as plain sample ids. Each
test runs against a mutant of the shipped code -- its reading before the rule, or the
reading the rule turned down -- and shows the mutant gives the other answer.

Hermetic, stdlib only. Run: python3 test_mnv_overlap_rules.py
"""

import ast
import bisect
import os
import types
import unittest

import indel_snp_cis as isc

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = "S1"
GT = {SAMPLE: "1|1"}


# --- the shipped code, and mutants of it ------------------------------------------

def _source(module, subs=()):
    """A module's source with each (text, replacement) of ``subs`` applied once: a
    mutant of the shipped code, to show that a test catches it."""
    path = os.path.join(HERE, module + ".py")
    with open(path) as fh:
        src = fh.read()
    for old, new in subs:
        if old not in src:
            raise AssertionError("%s no longer holds %r" % (path, old))
        src = src.replace(old, new, 1)
    return path, src


def mutant_module(module, subs):
    """``module`` loaded from a mutated source, under a name of its own."""
    path, src = _source(module, subs)
    mod = types.ModuleType(module + "_mutant")
    mod.__file__ = path
    exec(compile(src, path, "exec"), mod.__dict__)
    return mod


def functions(module, names, namespace, subs=()):
    """The named top-level functions of a module's source, compiled into ``namespace``
    (``subs``: a mutant of the shipped code). For the modules a test cannot import:
    analisi_indels_NNN.py runs a whole post-analysis at import time."""
    path, src = _source(module, subs)
    body = [n for n in ast.parse(src).body
            if isinstance(n, ast.FunctionDef) and n.name in names]
    if sorted(n.name for n in body) != sorted(names):
        raise AssertionError("%s: %s" % (path, sorted(n.name for n in body)))
    exec(compile(ast.Module(body=body, type_ignores=[]), path, "exec"), namespace)
    return namespace


class Reference(str):
    """The chromosome's bases, indexed by their own 0-based positions (genomeStr)."""

    def __new__(cls, at, seq):
        self = str.__new__(cls, seq)
        self.at = at
        return self

    def __getitem__(self, k):
        if isinstance(k, slice):
            return str.__getitem__(self, slice(k.start - self.at, k.stop - self.at))
        return str.__getitem__(self, k - self.at)

    def __len__(self):
        return self.at + str.__len__(self)

    def at_1based(self, pos1, n=1):
        """The ``n`` bases from the 1-based position ``pos1``."""
        return str.__getitem__(self, slice(pos1 - 1 - self.at, pos1 - 1 - self.at + n))

    def seq_at(self, a, b):
        """The bases a..b-1 (0-based), as the carrier gate reads them."""
        return "".join(self[q] for q in range(max(a, self.at), min(b, len(self))))


class Store(object):
    """A phased indel genotype store: {(POS, REF, ALT): {sample: genotype}}."""

    def __init__(self, records, has_mnv=True):
        self.records = dict(records)
        self.has_mnv = has_mnv

    def indels(self):
        return sorted(self.records)

    def carriers_dict(self, pos, ref, alt):
        return dict(self.records.get((pos, ref, alt), {}))


ANALYZER_NAMES = ("_store_indels", "_other_indels", "_on_slot", "_slot_allele", "_frame",
                  "_breaking_slots", "_cache_put", "_ref_snps", "_ref_snp_positions")


def analyzer(records, reference, snps=None, has_mnv=True, subs=(), cis=isc):
    """analisi_indels_NNN's record readers, run as the analyzer defines them."""
    snps = snps or {}
    ns = {"bisect": bisect, "_isc": cis, "_jgt_cache": {}, "_cache_held": {},
          "_CACHE_CARRIERS": 1000000, "_indel_index": None, "_indel_starts": [],
          "_indel_maxref": 0, "_indelgt": Store(records, has_mnv),
          "genomeStr": reference, "_snp_record": lambda pos1: snps.get(pos1, {})}
    return functions("analisi_indels_NNN", ANALYZER_NAMES, ns, subs)


def gate(subs=(), cis=isc):
    """The carrier gate's reading of a record against a row's SNP alleles."""
    return functions("haplotype_prioritizer", ("breakers_of", "holds_target"),
                     {"isc": cis}, subs)


# --- the shapes (HPRC v2.1, hg38) --------------------------------------------------

# chr9:67,744,557-67,744,573. The MNV 67,744,565 TTTC>ATAA keeps its second base (T) and
# the SNV 67,744,566 T>C stands on it; the carrier's own assembly holds the SNV's C
# there, and the pre-vcfwave original of that base (>202574658>202574661) names C.
CHR9B = Reference(67744556, "TGTTCCTCTTTCCTTTT")
MNV_KEEPS = (67744565, "TTTC", "ATAA")
SNV_ON_A_KEPT_BASE = (67744566, "T", "C")

# chr13:29,389,294-29,389,377. The MNV 29,389,333 ATATGTATG>GTATGTATA with 29,389,333
# A>G at its own POS and 29,389,341 G>A at its last base: each names the base the MNV
# puts there, so the two readings agree (21 of the 28 pairs are of this kind).
CHR13 = Reference(29389293, "TGTGTATATGTGTATATATGTATACACATATGTGTGTATATATGTATGCACGTGTGTG"
                            "TATATATGTATGCACGTGTGTGTATA")
MNV_SAME = (29389333, "ATATGTATG", "GTATGTATA")
SNV_SAME_AT_POS = (29389333, "A", "G")
SNV_SAME_LATER = (29389341, "G", "A")

# chrY:21,604,661-21,604,670. The one true third-base conflict: the MNV replaces the C at
# 21,604,663 with G and the SNV names T. The slot is one chrY copy (GT 1|.), the records'
# CONFLICT field names another sample, and the assembly reads CAGATACACC, neither
# reading: a decomposition artefact, so the slot leaves both records' rows.
CHRY = Reference(21604660, "AACACACACA")
MNV_CONFLICT = (21604661, "AAC", "CAG")
SNV_CONFLICT = (21604663, "C", "T")

# chr9:18,267,773-18,267,794, an A run. The complex record 18,267,781 C>AAAA and the
# insertion C>CA at its POS: 94 slots, the largest of the 35 same-POS insertion pairs.
CHR9A = Reference(18267772, "AACAAAAACAAAAACAAAAAAA")
CPLX_INS = (18267781, "C", "AAAA")
INS_AT_ITS_POS = (18267781, "C", "CA")

# chr12:117,194,726-117,194,795, a CTTT repeat. Three records of one locus:
#   the deletion 117,194,730 TTTTCTTTC>T, which keeps its anchor;
#   the complex record at its POS, 117,194,730 TTTTCTTTCTTT>GCTTG, which keeps none
#     (33 slots, one of the 15 same-POS deletion pairs);
#   the complex record 117,194,731 TTT>G, whose REF lies wholly inside the bases the
#     deletion deletes (20 slots, one of the 31 pairs of that shape).
CHR12 = Reference(117194725, "CTTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT"
                             "CTTTCCTTCT")
DEL_WITH_ANCHOR = (117194730, "TTTTCTTTC", "T")
CPLX_DEL = (117194730, "TTTTCTTTCTTT", "GCTTG")
CPLX_IN_DELETED = (117194731, "TTT", "G")


class TheSlotIsReadBaseByBase(unittest.TestCase):
    """Rule 2. An SNV and an MNV on one phased slot: the MNV's bases stand where it
    replaces them, the SNV's base wins on a base the MNV keeps, and a true third-base
    conflict leaves the slot off both records' rows (METHODS section 4)."""

    # --- the bases of the slot (indel_snp_cis.local_haplotype)

    def test_a_base_the_mnv_keeps_stays_a_reference_place(self):
        """The slot's own allele there is its base, so the SNV's C is read at
        67,744,566; the bases the MNV replaces are its own."""
        hap = dict(isc.local_haplotype(67744563, 67744568, None, [MNV_KEEPS]))
        self.assertEqual([hap[(q, 0)] for q in range(67744564, 67744568)],
                         ["A", 67744565, "A", "A"])
        # with the slot's alleles: the SNV's base on the base the MNV keeps. The frame
        # is the reference C, the MNV's ATAA with the SNV's C in its second place, and
        # the reference C past it (hg38 there: CTTTCC)
        slot = {67744565: "C"}  # 0-based 67,744,565 is the 1-based 67,744,566
        self.assertEqual(CHR9B.at_1based(67744564, 6), "CTTTCC")
        bases = "".join(s if isinstance(s, str) else slot.get(s, CHR9B[s])
                        for _k, s in isc.local_haplotype(67744563, 67744568, None,
                                                         [MNV_KEEPS]))
        self.assertEqual(bases, "CACAAC")

    def test_the_mutant_writes_the_mnvs_copy_of_the_kept_base(self):
        """2.10.0's reading, which hid the SNV: every base of the MNV a literal."""
        mut = mutant_module("indel_snp_cis", [(
            "b.upper() if b.upper() != cut[k].upper() else start + k", "b.upper()")])
        hap = dict(mut.local_haplotype(67744563, 67744568, None, [MNV_KEEPS]))
        self.assertEqual([hap[(q, 0)] for q in range(67744564, 67744568)],
                         ["A", "T", "A", "A"])

    # --- which SNP alleles leave the slot with no sequence (allele_breaks_record)

    def _breaks(self, record, snv, cis=isc):
        pos1, _ref, alt = snv
        ref = CHR9B if record is MNV_KEEPS else (
            CHRY if record is MNV_CONFLICT else CHR13)
        return cis.allele_breaks_record(record[0], record[1], record[2], pos1 - 1, alt,
                                        ref[pos1 - 1])

    def test_only_a_true_third_base_conflict_breaks_an_mnv(self):
        self.assertFalse(self._breaks(MNV_KEEPS, SNV_ON_A_KEPT_BASE))
        self.assertFalse(self._breaks(MNV_SAME, SNV_SAME_AT_POS))
        self.assertFalse(self._breaks(MNV_SAME, SNV_SAME_LATER))
        self.assertTrue(self._breaks(MNV_CONFLICT, SNV_CONFLICT))
        # an allele the genotypes cannot put on the slot breaks it, as before
        self.assertTrue(isc.allele_breaks_record(*MNV_KEEPS, 67744565, None, "T"))

    def test_an_indel_is_read_as_before(self):
        """Any allele but the reference base on a base of an indel's or a complex
        record's REF: a VCF writes a haplotype holding both as one record."""
        for rec in ((67744565, "TTTC", "T"), (67744565, "T", "TAA"),
                    (67744565, "TTTC", "G")):
            self.assertTrue(isc.allele_breaks_record(*rec, 67744565, "C", "T"))
            self.assertFalse(isc.allele_breaks_record(*rec, 67744565, "T", "T"))

    def test_the_mutant_breaks_an_mnv_on_any_snp_alt(self):
        """2.10.0's reading: a SNP alt anywhere on the REF leaves no one sequence."""
        mut = mutant_module("indel_snp_cis", [(
            """    if not len(ref) == len(alt) > 1:  # an indel, a complex record, a record of one base
        return base.upper() != ref_base.upper()""",
            "    if True:\n        return base.upper() != ref_base.upper()")])
        self.assertTrue(self._breaks(MNV_KEEPS, SNV_ON_A_KEPT_BASE, cis=mut))
        self.assertTrue(self._breaks(MNV_SAME, SNV_SAME_LATER, cis=mut))

    # --- the MNV's own rows (analisi_indels_NNN._ref_snps)

    def _ref_snps(self, window, reference, snps, subs=()):
        """The carriers an MNV row loses to the SNP alts ``snps``, the target covering
        the whole record and four bases on either side."""
        ns = analyzer({window: GT}, reference, snps=snps, subs=subs)
        return ns["_ref_snps"](window[0], window[1], window[2], -4, len(window[1]) + 3)

    def test_an_mnv_row_loses_the_carrier_of_an_snv_on_a_base_it_keeps(self):
        self.assertEqual(self._ref_snps(MNV_KEEPS, CHR9B, {67744566: {"C": GT}}), [GT])

    def test_an_mnv_row_keeps_the_carrier_of_an_snv_that_names_its_own_base(self):
        """At the MNV's own POS and at a later base: 21 of the 28 pairs."""
        self.assertEqual(self._ref_snps(MNV_SAME, CHR13, {29389333: {"G": GT}}), [])
        self.assertEqual(self._ref_snps(MNV_SAME, CHR13, {29389341: {"A": GT}}), [])
        self.assertEqual(self._ref_snps(MNV_CONFLICT, CHRY, {21604663: {"G": GT}}), [])

    def test_an_mnv_row_loses_the_carrier_of_a_true_third_base_conflict(self):
        self.assertEqual(self._ref_snps(MNV_CONFLICT, CHRY, {21604663: {"T": GT}}), [GT])

    def test_an_mnv_row_keeps_a_carrier_whose_snv_lies_outside_the_target(self):
        ns = analyzer({MNV_KEEPS: GT}, CHR9B, snps={67744566: {"C": GT}})
        self.assertEqual(ns["_ref_snps"](*MNV_KEEPS, 2, 6), [])  # the target past it

    def test_an_indel_row_reads_its_own_pos_as_before(self):
        """A SNP at an indel's POS is another allele of that base, a phasing artefact,
        and one on a base it deletes is not in the window's sequence."""
        ins = (67744565, "T", "TAA")
        self.assertEqual(self._ref_snps(ins, CHR9B, {67744565: {"C": GT}}), [])
        dele = (67744565, "TTTC", "T")
        self.assertEqual(self._ref_snps(dele, CHR9B, {67744566: {"C": GT}}), [])

    def test_the_mutant_reads_an_mnv_window_as_2100_did(self):
        """2.10.0 read only the bases an indel keeps past its POS, of which an MNV has
        none, so an MNV row kept every carrier whatever its own alleles."""
        subs = [("""    if len(iref) != len(ialt):
        return range(ipos + 1, ipos + _isc.kept_bases(iref, ialt))""",
                 "    if True:\n        return range(ipos + 1, ipos + _isc.kept_bases("
                 "iref, ialt))")]
        self.assertEqual(self._ref_snps(MNV_KEEPS, CHR9B, {67744566: {"C": GT}},
                                        subs=subs), [])
        self.assertEqual(self._ref_snps(MNV_CONFLICT, CHRY, {21604663: {"T": GT}},
                                        subs=subs), [])

    # --- the SNP row (indel_snp_cis.snp_target_held, new_simple_analysis)

    def _snp_row_holds(self, reference, snv, record, lo, hi, cis=isc):
        """Whether a SNP row's target lo..hi (0-based) is still on the haplotype that
        holds the row's SNV and the store's record."""
        alleles = {snv[0] - 1: snv[2]}

        def value(q):
            return alleles.get(q, reference[q]) if 0 <= q < len(reference) else None

        return cis.snp_target_held(lo, hi, sorted(alleles), [[record], [record]], value,
                                   reference)

    def test_a_snp_row_keeps_a_carrier_whose_snv_the_mnv_keeps_out_of_its_target(self):
        """The MNV keeps the base the row's SNV stands on, so the SNV's base is the
        slot's; with the MNV's other bases outside the target the row's target is
        still there."""
        self.assertTrue(self._snp_row_holds(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                            67744565, 67744565))

    def test_a_snp_row_loses_the_carrier_where_the_mnv_changes_a_base_of_the_target(self):
        self.assertFalse(self._snp_row_holds(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                             67744564, 67744567))

    def test_a_snp_row_loses_the_carrier_of_a_true_third_base_conflict(self):
        """The slot has no one sequence there, so it leaves this row too: with the
        conflict the target cannot be read even where the MNV changes nothing of it."""
        self.assertFalse(self._snp_row_holds(CHRY, SNV_CONFLICT, MNV_CONFLICT,
                                             21604662, 21604662))

    def test_the_mutant_drops_the_carrier_of_an_snv_on_a_base_the_mnv_keeps(self):
        """2.10.0's reading: a SNP alt on any base of the record's REF leaves the slot
        with no sequence, so the row loses the carrier."""
        mut = mutant_module("indel_snp_cis", [(
            """    if not len(ref) == len(alt) > 1:  # an indel, a complex record, a record of one base
        return base.upper() != ref_base.upper()""",
            "    if True:\n        return base.upper() != ref_base.upper()")])
        self.assertFalse(self._snp_row_holds(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                             67744565, 67744565, cis=mut))

    # --- the carrier gate (haplotype_prioritizer.breakers_of)

    def _gate_breaks(self, reference, snv, record, lo, hi, subs=(), cis=isc):
        """Whether the gate reads ``record`` as taking the target's bases lo..hi off a
        haplotype that holds the row's SNV."""
        ns = gate(subs, cis)
        changing = [(record[0], record[1], record[2], GT)]
        return bool(ns["breakers_of"](changing, {snv[0]: snv[2]}, lo, hi,
                                      reference.seq_at))

    def test_the_gate_agrees_on_a_base_the_mnv_keeps(self):
        """Outside the target the MNV takes no haplotype off; over it, it does."""
        self.assertFalse(self._gate_breaks(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                           67744565, 67744565))
        self.assertTrue(self._gate_breaks(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                          67744564, 67744567))

    def test_the_gate_agrees_at_the_mnvs_own_pos(self):
        """A row SNP at the MNV's POS naming the base it puts there: no artefact, and
        the MNV is read against the target like any other record."""
        self.assertFalse(self._gate_breaks(CHR13, SNV_SAME_AT_POS, MNV_SAME,
                                           29389320, 29389325))
        self.assertTrue(self._gate_breaks(CHR13, SNV_SAME_AT_POS, MNV_SAME,
                                          29389332, 29389340))

    def test_the_gate_agrees_on_a_true_third_base_conflict(self):
        self.assertTrue(self._gate_breaks(CHRY, SNV_CONFLICT, MNV_CONFLICT,
                                          21604661, 21604661))

    def test_the_gate_reads_an_indel_as_before(self):
        """A row SNP at an indel's POS is another allele of that base (no breaker), one
        on another base of its REF takes the haplotype off whatever the reference."""
        self.assertFalse(self._gate_breaks(CHR9B, (67744565, "T", "C"),
                                           (67744565, "TTTC", "T"), 67744560, 67744570))
        self.assertTrue(self._gate_breaks(CHR9B, (67744566, "T", "C"),
                                          (67744565, "TTTC", "T"), 67744560, 67744570))

    def test_the_gate_mutant_splits_the_pos_from_the_bases_after_it(self):
        """2.10.0's reading: a row SNP at the record's POS took it off the gate's
        reading, and one anywhere else on its REF took the haplotype off the row."""
        subs = [("""        mnv = len(r) == len(a) > 1""", "        mnv = False")]
        self.assertTrue(self._gate_breaks(CHR9B, SNV_ON_A_KEPT_BASE, MNV_KEEPS,
                                          67744565, 67744565, subs=subs))
        self.assertFalse(self._gate_breaks(CHR13, SNV_SAME_AT_POS, MNV_SAME,
                                           29389332, 29389340, subs=subs))


class AnIndelAtAComplexRecordsPos(unittest.TestCase):
    """Rule 3. One of the two keeps no anchor base, so the two describe the same bases:
    the slot is read with the complex record alone, on whose rows the indel is not read,
    while the indel's rows read the complex record and the slot leaves them."""

    def lost(self, window, others, reference, span, subs=()):
        """The carriers the window's rows lose to ``others`` over the target's span
        (_window_breakers: the records _other_indels gives, with their carriers)."""
        ns = analyzer({window: GT, **{o: GT for o in others}}, reference, subs=subs)
        return [key for key, gts, _overlapping
                in ns["_other_indels"](span[0], span[1], *window) if gts]

    def test_an_insertion_at_a_complex_records_pos(self):
        """chr9:18,267,781 C>AAAA with C>CA at its POS. The complex record's rows keep
        the slot, the insertion's rows lose it. Before the rule the slot left both."""
        span = (18267778, 18267790)
        self.assertEqual(self.lost(CPLX_INS, [INS_AT_ITS_POS], CHR9A, span), [])
        self.assertEqual(self.lost(INS_AT_ITS_POS, [CPLX_INS], CHR9A, span), [CPLX_INS])

    def test_a_deletion_at_a_complex_records_pos_is_unchanged(self):
        """chr12:117,194,730 TTTTCTTTCTTT>GCTTG with TTTTCTTTC>T at its POS: the
        deletion deletes bases the complex record replaces, so the complex record's rows
        already kept the slot, and its own rows already lost it."""
        span = (117194729, 117194745)
        self.assertEqual(self.lost(CPLX_DEL, [DEL_WITH_ANCHOR], CHR12, span), [])
        self.assertEqual(self.lost(DEL_WITH_ANCHOR, [CPLX_DEL], CHR12, span), [CPLX_DEL])

    def test_two_records_that_both_keep_their_anchor_are_unchanged(self):
        """Two alleles of one site (an insertion and a deletion at one POS): neither is
        read on the other's rows, as since 2.8.1."""
        span = (18267778, 18267790)
        deletion = (18267781, "CA", "C")
        self.assertEqual(self.lost(INS_AT_ITS_POS, [deletion], CHR9A, span), [])
        self.assertEqual(self.lost(deletion, [INS_AT_ITS_POS], CHR9A, span), [])

    def test_two_records_that_both_keep_none_are_unchanged(self):
        """chrX:35,966,203 has the only pairs of two complex records at one POS (3
        slots): neither is a displaced copy of the other, so each is read on the other's
        rows and the slot leaves both."""
        span = (18267778, 18267790)
        other = (18267781, "C", "TTTT")
        self.assertEqual(self.lost(CPLX_INS, [other], CHR9A, span), [other])
        self.assertEqual(self.lost(other, [CPLX_INS], CHR9A, span), [CPLX_INS])

    def test_the_mutant_reads_the_insertion_on_the_complex_records_rows(self):
        """The skip before the rule, which needed both records to keep an anchor: the
        insertion is read, the slot leaves the complex record's rows too, and the pair's
        targets go unreported (the reading measured at 48.5 % real)."""
        subs = [("""        if p == ipos and _isc.kept_bases(r, a) and (kept or len(iref) != len(ialt)):""",
                 "        if p == ipos and kept and _isc.kept_bases(r, a):")]
        span = (18267778, 18267790)
        self.assertEqual(self.lost(CPLX_INS, [INS_AT_ITS_POS], CHR9A, span, subs=subs),
                         [INS_AT_ITS_POS])
        # and the deletion, which the rule does not move, is read the same way
        self.assertEqual(self.lost(CPLX_DEL, [DEL_WITH_ANCHOR], CHR12,
                                   (117194729, 117194745), subs=subs), [])


class AComplexRecordInsideDeletedBases(unittest.TestCase):
    """Rule 1 (item 23). A complex record whose REF lies wholly inside the bases a
    deletion deletes, on one phased slot: the slot leaves both records' rows and is
    listed nowhere, as up to 2.10.0. No code change -- this pins the reading so that a
    later refactor cannot move it, and runs the reading that was turned down."""

    SPAN = (117194729, 117194745)

    def lost(self, window, others, subs=()):
        ns = analyzer({window: GT, **{o: GT for o in others}}, CHR12, subs=subs)
        return [key for key, gts, _overlapping
                in ns["_other_indels"](self.SPAN[0], self.SPAN[1], *window) if gts]

    def test_the_slot_is_listed_nowhere(self):
        """chr12:117,194,730 TTTTCTTTC>T with 117,194,731 TTT>G inside its deleted
        bases. The complex record puts its own base among the bases the deletion
        deletes, as an insertion there does, so the slot leaves the deletion's rows;
        the deletion deletes bases the complex record keeps, so it leaves its rows
        too."""
        self.assertEqual(self.lost(DEL_WITH_ANCHOR, [CPLX_IN_DELETED]),
                         [CPLX_IN_DELETED])
        self.assertEqual(self.lost(CPLX_IN_DELETED, [DEL_WITH_ANCHOR]),
                         [DEL_WITH_ANCHOR])

    def test_the_reading_that_was_turned_down(self):
        """"It changes nothing", so the slot would stay on the deletion's rows: an ALT
        base put on bases the haplotype has lost already read as changing nothing.
        Measured over the 31 pairs (249 slots) against HPRC's pre-vcfwave originals,
        that reading lists 4,864 23-mers of which 692 are false, in 48 slots; outside
        the chr1:977,063 VNTR it is 676 real against 692 false. Today's reading lists
        none of them, and no false one."""
        subs = [("""    if len(alt) - t > c and len(alt) != len(ref):  # bases inserted between a+c-1 and a+c
        return lo <= a + c - 1 and a + c <= hi""",
                 """    if (len(alt) - t > c and len(alt) != len(ref)
            and (gone is None or not gone[0] <= a + c <= gone[1] + 1)):
        return lo <= a + c - 1 and a + c <= hi""")]
        mut = mutant_module("indel_snp_cis", subs)
        self.assertTrue(isc.changes_span(*CPLX_IN_DELETED, *self.SPAN,
                                         (117194730, 117194737)))
        self.assertFalse(mut.changes_span(*CPLX_IN_DELETED, *self.SPAN,
                                          (117194730, 117194737)))
        ns = analyzer({DEL_WITH_ANCHOR: GT, CPLX_IN_DELETED: GT}, CHR12, cis=mut)
        self.assertEqual([key for key, gts, _o in ns["_other_indels"](
            self.SPAN[0], self.SPAN[1], *DEL_WITH_ANCHOR) if gts], [])
        # and the reading of a variant on bases the deletion deletes is untouched: a
        # second deletion of those bases, and an MNV on them, still change nothing
        for rec in ((117194731, "TTTT", "T"), (117194731, "TTT", "GGG")):
            self.assertFalse(isc.changes_span(*rec, *self.SPAN, (117194730, 117194737)))


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""A SNP row's carrier whose haplotype holding the row's SNPs also holds an indel that
changes the target's bases (2.10.0).

The enriched genome holds SNPs only, so a SNP row's carriers come from SNP genotypes. Up
to 2.8.1 a carrier whose haplotype also held an indel inside the target was listed,
though its own sequence there is the indel's (chr22 of the HGDP + 1000 Genomes panel:
21052641 G>A is on the same 1,448 haplotypes as 21052640 T>TCACC, and a SNP-only row at
21052638- listed their 1,292 carriers). In a per-sample search such a carrier now leaves
the row, read on the phased indel genotype store as the 2.8.1 rule reads an indel row's
second indel:

  * phased: decided by the haplotype (cis: off the row; trans: stays);
  * homozygous for an indel that changes the target's bases: on every haplotype, off
    the row;
  * an unphased heterozygous indel (or SNP) call: stays, PUTATIVE;
  * an indel leaving the target's bases as they are (a repeat unit at the target's
    edge) or whose deleted bases all lie outside it (its junction at the edge), or an
    insertion before the first base: stays;
  * an indel at the position of a SNP the row names: another allele of that base, a
    phasing artefact, stays;
  * a row left with no carrier is not reported;
  * population-level, and CRISPRME_INDEL_SNP=0 (the opt-out of the SNP+indel
    co-occurrence): the store is not read, unchanged.

Covered here on both per-sample paths (the observed-haplotype path of the dict-less
tiers and the dictionary path), through the AST harness of test_phased_haplotype, plus
the pure helpers (indel_snp_cis), the enumerator's slot filter (observed_haplotypes) and
the store reader.

Run with:
    cd PostProcess && python3 -m unittest test_snp_row_indel_carriers -v
"""

import ast
import contextlib
import gzip
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

import indel_snp_cis as isc  # noqa: E402
import observed_haplotypes as oh  # noqa: E402
from test_phased_haplotype import _load_pure_functions  # noqa: E402

G = "ACTGAAATCTGTAAGCAGGC"  # the guide
T0 = 40  # the target's first base (0-based)
LEFT = ("CATG" * 10)[:T0]
SITE = G[:5] + "C" + G[6:] + "TGG"  # a C>A SNP at the 6th base makes it perfect
RIGHT = "G" + "ACATCAGTCA" * 4  # a G after the PAM: TGGG
GENOME = LEFT + SITE + RIGHT
SNP_POS1 = T0 + 6  # the SNP's VCF POS
LO, HI = T0, T0 + len(SITE) - 1


def _ins_inside():
    """A 2-bp insertion after the target's 13th base (inside the protospacer)."""
    p0 = T0 + 12
    return (p0 + 1, GENOME[p0], GENOME[p0] + "TT")


def _gene_run_del():
    """TG>T at the PAM's T: one G of the TGGG run goes, the run puts a G back."""
    p0 = T0 + 20
    assert GENOME[p0:p0 + 4] == "TGGG"
    return (p0 + 1, "TG", "T")


def _spacer_del():
    """A 1-bp deletion of the target's 4th base."""
    p0 = T0 + 2
    return (p0 + 1, GENOME[p0:p0 + 2], GENOME[p0])


def _edge_del():
    """A 3-bp deletion whose last deleted base is the one before the target."""
    p0 = T0 - 4
    return (p0 + 1, GENOME[p0:T0], GENOME[p0])


def _write_store(path, indels):
    with gzip.open(path, "wt") as fh:
        for (p, r, a), gts in indels:
            fh.write("%d_%s_%s\t%s\n" % (p, r, a, ",".join("%s:%s" % kv for kv in gts.items())))


COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _mnv_inside():
    """An MNV on the target's 10th and 11th bases (every base changed: TG>CA)."""
    p0 = T0 + 9
    ref = GENOME[p0:p0 + 2]
    return (p0 + 1, ref, ref.translate(str.maketrans("ACGT", "TGCA")))


def _mnv_first():
    """An MNV on the target's first two bases."""
    ref = GENOME[T0:T0 + 2]
    return (T0 + 1, ref, ref.translate(str.maketrans("ACGT", "TGCA")))


def _mnv_at_the_rows_snp():
    """An MNV at the position of the SNP the row names, which puts the row's own alt
    there (as 21 of HPRC v2.1's 28 SNV-MNV pairs do) and changes the base after it,
    inside the target (chr13:29,389,333 ATATGTATG>GTATGTATA beside 29,389,333 A>G)."""
    p0 = SNP_POS1 - 1
    return (SNP_POS1, GENOME[p0:p0 + 2], "A" + COMP[GENOME[p0 + 1]])


def _mnv_keeping_the_rows_snp():
    """An MNV that spans the whole target (HPRC v2.1 has MNVs up to 75,266 bp), keeping
    every base of it -- the one the row's SNP stands on included -- and changing two
    bases on either side. Read base by base (2.10.1) it changes nothing of the target
    and the row's SNP is the slot's base there, as the pre-vcfwave originals read the
    bases an MNV keeps (chr9:67,744,565 TTTC>ATAA beside 67,744,566 T>C)."""
    p0, p1 = T0 - 2, HI + 3
    ref = GENOME[p0:p1]
    alt = "".join(COMP[b] if q < T0 or q > HI else b
                  for q, b in zip(range(p0, p1), ref))
    return (p0 + 1, ref, alt)


class _Fixture(object):
    """A SNP target at T0 (+ strand) whose SNP C>A makes it perfect, a phased indel
    store, and the harness globals of one per-sample (or population-level) run."""

    def __init__(self, snp_gts, indels, observed=True, phased=True, fast=False, mnvs=None):
        self.tmp = tempfile.TemporaryDirectory()
        store = os.path.join(self.tmp.name, "gt_indel_chrT.tsv.gz")
        _write_store(store, indels)
        if mnvs is not None:  # the store's MNV file, beside its indel file (2.10.1)
            _write_store(os.path.join(self.tmp.name, "gt_mnv_chrT.tsv.gz"), mnvs)
        dna = G[:5] + "M" + G[6:] + "TGG"
        self.split = ["X", G + "NGG", dna, "chrT", str(T0), str(T0), "+", "1", "0", "1",
                      "NGG", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        mydict = {"chrT,%d" % SNP_POS1: ",".join("%s:%s" % kv for kv in snp_gts.items())
                  + ";C,A;rs1;0.10"}
        keys = {"chrT,%d" % (p + 1) for p in range(len(GENOME))}
        overrides = dict(
            genomeStr=GENOME, current_chr="chrT", mydict=mydict, myreg=None,
            mygt=_GtSentinel() if observed else None, haplotype_check=phased,
            IUPAC_CAP=10, hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
            pam_begin=-3, pam_end=None, allowed_mms=4, _dict_keys_read=keys,
            _FAST_MODE=fast,
        )
        self.ns = _load_pure_functions(overrides)
        self.ns["_snp_indels"] = self.ns["_SnpRowIndels"](store, [(LO, HI)], 500)
        self.ns["_phase_confirmation_rows"].clear()
        self.ns["_phase_confirmation_keys"].clear()

    def run(self):
        cluster = []
        self.ns["iupac_decomposition"](self.split, G + "NGG", G, cluster)
        self.tmp.cleanup()
        return cluster

    def snp_carriers(self, rows):
        """Carriers of the perfect SNP row (0 mismatches), or None when there is none."""
        perfect = [r for r in rows if r[12] not in ("n", "NA") and int(r[7]) == 0]
        if not perfect:
            return None
        return set().union(*[set(r[12].split(",")) for r in perfect])

    def phase(self):
        return {r["Phase_Confirmation"] for r in self.ns["_phase_confirmation_rows"]}


class _GtSentinel(object):
    """A truthy stand-in for a Tier-1 reader: selects the observed path, whose carriers
    come from the dictionary entries here."""

    def __len__(self):
        return 1


class PureHelpers(unittest.TestCase):
    def value(self, src):
        return src if isinstance(src, str) else GENOME[src]

    def test_window_none_is_the_reference(self):
        hap = isc.local_haplotype(10, 14, None)
        self.assertEqual(hap, [((p, 0), p) for p in range(10, 15)])

    def test_window_none_with_indels(self):
        p, r, a = _ins_inside()
        hap = isc.local_haplotype(T0, T0 + 15, None, [(p, r, a)])
        self.assertEqual("".join(s if isinstance(s, str) else GENOME[s] for _k, s in hap),
                         GENOME[T0:T0 + 13] + "TT" + GENOME[T0 + 13:T0 + 16])

    def test_changes_span_junction_rule(self):
        self.assertTrue(isc.changes_span(*_ins_inside(), LO, HI))
        self.assertTrue(isc.changes_span(*_spacer_del(), LO, HI))
        self.assertFalse(isc.changes_span(*_edge_del(), LO, HI))  # junction at the edge
        p0 = T0 - 1  # an insertion before the first base
        self.assertFalse(isc.changes_span(p0 + 1, GENOME[p0], GENOME[p0] + "AA", LO, HI))
        p0 = HI  # an insertion after the last base
        self.assertFalse(isc.changes_span(p0 + 1, GENOME[p0], GENOME[p0] + "AA", LO, HI))
        p0 = T0 - 4  # a deletion whose last deleted base is the target's first
        self.assertTrue(isc.changes_span(p0 + 1, GENOME[p0:T0 + 1], GENOME[p0], LO, HI))

    def test_held_on_bases(self):
        snps = [SNP_POS1 - 1]
        for rec, want in ((_gene_run_del(), True), (_spacer_del(), False),
                          (_ins_inside(), False)):
            self.assertEqual(isc.snp_target_held(LO, HI, snps, ([rec], [rec]),
                                                 self.value, GENOME), want, rec)

    def test_unreadable_way(self):
        rec = _gene_run_del()
        self.assertFalse(isc.snp_target_held(LO, HI, [SNP_POS1 - 1], (None, None),
                                             self.value, GENOME))
        # an unknown allele where the reading compares bases: not held
        unknown = lambda s: None if s == HI + 1 else self.value(s)  # noqa: E731
        self.assertFalse(isc.snp_target_held(LO, HI, [SNP_POS1 - 1], ([rec], [rec]),
                                             unknown, GENOME))

    def test_snp_pushed_out_of_the_stretch(self):
        # the target TGGGG (10..14) needs the row's SNP A>G at its last base; a G
        # inserted after the T reads TGGGG from the first base on with that SNP's base
        # pushed out of the stretch: the haplotype has those bases thanks to the
        # insertion, not the SNP, so it is not the SNP row's carrier
        ref = "AC" * 5 + "TGGGA" + "CA" * 5
        ins = (11, "T", "TG")
        row = {14: "G"}
        value = lambda s: s if isinstance(s, str) else row.get(s, ref[s])  # noqa: E731
        self.assertFalse(isc.snp_target_held(10, 14, [14], ([ins], [ins]), value, ref))
        # the row's SNP inside the stretch instead (A>G at 12): the bases from the
        # first base on are the target's, through the SNP
        ref = "AC" * 5 + "TGAGG" + "CA" * 5
        row = {12: "G"}
        self.assertTrue(isc.snp_target_held(10, 14, [12], ([ins], [ins]), value, ref))

    def test_two_records_of_one_base(self):
        rec = _gene_run_del()
        alt_on_ref = lambda s: "A" if s == rec[0] - 1 else self.value(s)  # noqa: E731
        self.assertFalse(isc.snp_target_held(LO, HI, [], ([rec], [rec]), alt_on_ref,
                                             GENOME))

    def test_target_pairs_through(self):
        hap = isc.local_haplotype(0, 30, None)
        ways = isc.target_pairs(hap, hap, (5, 0), (10, 0), through=[(7, 0)])
        self.assertEqual(ways, [[], []])
        ways = isc.target_pairs(hap, hap, (5, 0), (10, 0), through=[(12, 0)])
        self.assertEqual(ways, [None, None])


def _two_records_at(p0):
    """Two records of one base at 0-based ``p0`` (one VCF POS): a deletion of the next
    base and an insertion after it. A slot phased with both has no one sequence there."""
    return [(p0 + 1, GENOME[p0:p0 + 2], GENOME[p0]), (p0 + 1, GENOME[p0], GENOME[p0] + "A")]


class AVariantPastTheStretch(unittest.TestCase):
    """2.10.0's release gates (B10: VEGFA site 3, per sample, HGDP + 1000 Genomes) found
    an indel row's carrier left off for a variant its slot cannot be written with (two
    records of one base, a SNP alt on an indel's REF, an unphased call) 2 to 41 bases
    past the target, which the bases its reading compares never reach. A SNP row's rule
    read the slot the same way: such a variant anywhere in the frame on a way's side
    stopped that way. The slot is now written without it, and only a way that reaches a
    base of its REF -- a base of its stretch, or of the target's span -- cannot be read.

    The world: _gene_run_del (TG>T at the PAM's T) keeps the target's bases from its
    first base on, through one G of the TGGG run the run puts back; that reading's
    stretch ends one base past the target (HI + 1). Up to the last base it does not hold
    them. Each case puts a variant the slot cannot be written with on the reading from
    the first base on, at HI + 5 (past the stretch: kept) or at HI + 1 (the stretch's
    last base: not). On 2.10.0's code before the fix, every such slot was dropped."""

    def value(self, src):
        return src if isinstance(src, str) else GENOME[src]

    def held(self, held, value=None, unread=None):
        return isc.snp_target_held(LO, HI, [SNP_POS1 - 1], held, value or self.value,
                                   GENOME, *(() if unread is None else (unread,)))

    def test_two_records_of_one_base_past_the_stretch(self):
        run = _gene_run_del()
        self.assertTrue(self.held(([run] + _two_records_at(HI + 5), [run])))
        # on the stretch's last base: that way cannot be read
        self.assertFalse(self.held(([run] + _two_records_at(HI + 1), [run])))

    def test_an_unphased_indel_past_the_stretch(self):
        run = _gene_run_del()
        far = (HI + 6, GENOME[HI + 5], GENOME[HI + 5] + "TT")
        self.assertTrue(self.held(([run], [run]), unread=([far], [])))
        near = (HI + 2, GENOME[HI + 1], GENOME[HI + 1] + "TT")
        self.assertFalse(self.held(([run], [run]), unread=([near], [])))

    def test_a_snp_alt_on_an_indel_ref_past_the_stretch(self):
        run = _gene_run_del()
        for at, want in ((HI + 5, True), (HI + 1, False)):
            dele = (at + 1, GENOME[at:at + 2], GENOME[at])
            snp = "C" if GENOME[at + 1] != "C" else "A"  # an alt on the deleted base
            alt_on_ref = lambda s: snp if s == at + 1 else self.value(s)  # noqa: E731
            self.assertEqual(self.held(([run, dele], [run]), alt_on_ref), want, at)

    def test_without_them_the_same(self):
        """A slot with no such variant reads as it did: the run deletion kept, a deletion
        inside the protospacer dropped."""
        self.assertTrue(self.held(([_gene_run_del()], [_gene_run_del()])))
        self.assertFalse(self.held(([_spacer_del()], [_spacer_del()])))
        self.assertFalse(self.held(([_spacer_del()] + _two_records_at(HI + 5),
                                    [_spacer_del()])))


class EnumeratorSlotFilter(unittest.TestCase):
    """observed_haplotypes: the slot filter decides per haplotype."""

    @staticmethod
    def _positions(gts):
        return [{"pos_c": 5, "alts": ["A"], "carrier_gts": [gts], "alt_index": ["1"],
                 "info": [["rs1", "0.1", "chrT_46_C_A"]], "ps": None}]

    def _enum(self, gts, slot_filter):
        from population_summary import parse_haplotypes
        return oh.enumerate_observed_haplotypes(
            self._positions(gts), "C" * 23, parse_haplotypes, lambda s, x: 2, "chrT",
            slot_filter=slot_filter)

    def test_none_is_unchanged(self):
        gts = {"S1": "1|0", "S2": "0/1"}
        a = self._enum(gts, None)
        b = self._enum(gts, lambda s, h, v: oh.KEEP)
        self.assertEqual([(sorted(h.variant_set), sorted(h.carriers), h.phase_state)
                          for h in a],
                         [(sorted(h.variant_set), sorted(h.carriers), h.phase_state)
                          for h in b])

    def test_drop_and_putative(self):
        gts = {"S1": "1|0", "S2": "1|0", "S3": "0/1"}
        verdicts = {"S1": oh.DROP, "S2": oh.PUTATIVE, "S3": oh.DROP}
        haps = self._enum(gts, lambda s, h, v: verdicts[s])
        self.assertEqual(len(haps), 1)
        self.assertEqual(haps[0].carriers, {"S2"})
        self.assertEqual(haps[0].phase_state, oh.PUTATIVE)

    def test_homozygous_set_kept_by_one_slot(self):
        gts = {"S1": "1|1"}
        haps = self._enum(gts, lambda s, h, v: oh.DROP if h == 0 else oh.KEEP)
        self.assertEqual([(h.carriers, h.phase_state) for h in haps],
                         [({"S1"}, oh.CONFIRMED)])
        haps = self._enum(gts, lambda s, h, v: oh.DROP)
        self.assertEqual(haps, [])

    def test_unphased_candidates_asked_without_slot(self):
        asked = []
        self._enum({"S1": "0/1"}, lambda s, h, v: asked.append(h) or oh.KEEP)
        self.assertEqual(asked, [None])


class StoreReader(unittest.TestCase):
    def test_reads_only_near_the_targets(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "s.tsv.gz")
            near = _ins_inside()
            far = (2000, "A", "AT")
            _write_store(path, [(near, {"S1": "1|0"}), (far, {"S2": "0|1"})])
            ns = _load_pure_functions({})
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual(st.indels, [near])
            self.assertEqual(st.carriers(near), {"S1": "1|0"})
            self.assertEqual(st.touching(LO, HI), [near])
            self.assertTrue(st.covers(LO - 50, HI + 50))
            self.assertFalse(st.covers(LO - 51, HI))

    def test_an_indel_the_store_holds_twice_is_one(self):
        """R3B-N1: an indel written on two VCF lines has two lines in the store; the SNP rows'
        reader takes it once, its carriers merged (build_indel_genotypes.
        merge_store_carriers). It listed the indel twice and kept the later line's carriers:
        S1, on the earlier line only, kept a SNP row its slot's indel changes."""
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "s.tsv.gz")
            near = _ins_inside()
            _write_store(path, [(near, {"S1": "1|0", "S2": "0|1"}), (near, {"S2": "1|0"})])
            ns = _load_pure_functions({})
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual(st.indels, [near])
            self.assertEqual(st.carriers(near), {"S1": "1|0", "S2": "1|1"})

    def test_an_indel_the_store_holds_in_two_cases_is_one(self):
        """C-R1-1: two lines of one indel in two cases (HPRC chr22: 140, tg>t beside TG>T)
        are one allele, read in upper case (build_indel_genotypes.store_key), their
        carriers merged. Keyed as written, the SNP rows' reader held two indels at one
        position, a slot with both of them one with two records of one base."""
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "s.tsv.gz")
            near = _ins_inside()
            low = (near[0], near[1].lower(), near[2].lower())
            _write_store(path, [(near, {"S1": "1|0", "S2": "0|1"}), (low, {"S2": "1|0"})])
            ns = _load_pure_functions({})
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual(st.indels, [near])
            self.assertEqual(st.carriers(near), {"S1": "1|0", "S2": "1|1"})

    def test_reads_the_mnv_file_beside_the_indel_file(self):
        """From 2.10.1 the store's MNVs are in a file beside its indel file
        (gt_mnv_<chrom>.tsv.gz), read with it: an MNV inside a target touches it. A store
        of 2.10.0 has none. At c64e62e1 the reader read the indel file alone."""
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "gt_indel_chrT.tsv.gz")
            _write_store(path, [(_edge_del(), {"S1": "1|0"})])
            ns = _load_pure_functions({})
            self.assertEqual(ns["_SnpRowIndels"](path, [(LO, HI)], 50).indels, [_edge_del()])
            _write_store(os.path.join(d, "gt_mnv_chrT.tsv.gz"),
                         [(_mnv_inside(), {"S2": "0|1"}), ((2000, "AC", "GT"), {"S3": "1|0"})])
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual(st.indels, sorted([_edge_del(), _mnv_inside()]))
            self.assertEqual(st.carriers(_mnv_inside()), {"S2": "0|1"})
            self.assertEqual(st.touching(LO, HI), [_mnv_inside()])

    def test_an_mnv_file_that_cannot_be_read_is_none(self):
        """The review's item 4: an MNV file cut short, not gzip at all, or of 0 bytes, is
        read as no MNV file (a store of 2.10.0: no SNP row loses a carrier for an MNV), the
        indel file as it is, mnv_error saying why; at 516d42a0 the reader raised and the run
        read the chromosome with no store at all, and at 6b5ab4d3 a file of 0 bytes was read
        as one of no MNV, with no NOTE. A corrupt indel file still raises."""
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "gt_indel_chrT.tsv.gz")
            mnv = os.path.join(d, "gt_mnv_chrT.tsv.gz")
            _write_store(path, [(_edge_del(), {"S1": "1|0"})])
            _write_store(mnv, [(_mnv_inside(), {"S2": "0|1"})])
            with open(mnv, "rb") as fh:
                whole = fh.read()
            ns = _load_pure_functions({})
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual((st.indels, st.mnv_error), (sorted([_edge_del(), _mnv_inside()]),
                                                         None))
            for data in (whole[:len(whole) // 2], b"plain text\n", b""):
                with open(mnv, "wb") as fh:
                    fh.write(data)
                st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
                self.assertEqual(st.indels, [_edge_del()])
                self.assertEqual(st.carriers(_edge_del()), {"S1": "1|0"})
                self.assertEqual(st.touching(LO, HI), [])
                self.assertTrue(st.mnv_error)
                self.assertEqual(st.mnv_path, mnv)
            with open(path, "wb") as fh:
                fh.write(b"plain text\n")
            with self.assertRaises(OSError):
                ns["_SnpRowIndels"](path, [(LO, HI)], 50)

    def test_touching_junction_rule(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "s.tsv.gz")
            _write_store(path, [(_edge_del(), {"S1": "1|0"}),
                                (_spacer_del(), {"S1": "1|0"})])
            ns = _load_pure_functions({})
            st = ns["_SnpRowIndels"](path, [(LO, HI)], 50)
            self.assertEqual(st.touching(LO, HI), [_spacer_del()])


class ObservedPath(unittest.TestCase):
    """The per-sample observed-haplotype path (Tier-1 genotypes)."""

    observed = True

    def _carriers(self, snp_gts, indels, **kw):
        fx = _Fixture(snp_gts, indels, observed=self.observed, **kw)
        rows = fx.run()
        return fx.snp_carriers(rows), fx

    def test_cis_trans_homozygous(self):
        got, fx = self._carriers(
            {"S1": "1|0", "S2": "1|0", "S3": "1|0", "S6": "1|0"},
            [(_ins_inside(), {"S1": "1|0", "S2": "0|1", "S3": "1|1"})])
        self.assertEqual(got, {"S2", "S6"})
        st = fx.ns["_snp_indel_stats"]
        self.assertEqual((st["targets"], st["listings"], st["gone"]), (1, 2, 0))

    def test_unphased_stays_putative(self):
        got, fx = self._carriers(
            {"S4": "0/1", "S5": "1|0"},
            [(_ins_inside(), {"S4": "0/1", "S5": "0/1"})])
        self.assertEqual(got, {"S4", "S5"})
        if self.observed:
            self.assertEqual(fx.phase(), {"PUTATIVE"})

    def test_unphased_homozygous_indel_dropped(self):
        got, _fx = self._carriers({"S4": "0/1", "S6": "1|0"},
                                  [(_ins_inside(), {"S4": "1/1"})])
        self.assertEqual(got, {"S6"})

    def test_row_left_with_no_carrier_not_reported(self):
        got, fx = self._carriers({"S1": "1|0"}, [(_ins_inside(), {"S1": "1|0"})])
        self.assertIsNone(got)
        self.assertEqual(fx.ns["_snp_indel_stats"]["gone"], 1)

    def test_bases_left_as_they_are(self):
        got, fx = self._carriers({"S7": "1|0", "S8": "1|0"},
                                 [(_gene_run_del(), {"S7": "1|0"}),
                                  (_spacer_del(), {"S8": "1|0"})])
        self.assertEqual(got, {"S7"})
        self.assertEqual(fx.ns["_snp_indel_stats"]["kept"], 1)

    def test_reading_depends_on_the_haplotypes_own_alleles(self):
        # two carriers of the G-run deletion; S8 also has G>T right after the run, on the
        # same slot: the run no longer puts a G back, so only S7 has the target (the
        # cached reading is keyed by each haplotype's alleles around the target)
        fx = _Fixture({"S7": "1|0", "S8": "1|0"}, [(_gene_run_del(), {"S7": "1|0", "S8": "1|0"})],
                      observed=self.observed)
        fx.ns["mydict"]["chrT,%d" % (HI + 2)] = "S8:1|0;G,T;rs2;0.10"
        self.assertEqual(fx.snp_carriers(fx.run()), {"S7"})

    def test_a_variant_past_the_stretch_read_does_not_stop_it(self):
        """2.10.0's release gates (AVariantPastTheStretch): S7's and S9's slot of the SNP
        holds the G-run deletion, whose bases read from the target's first base on are
        the target's, and past that stretch two records of one base (S7) or an unphased
        insertion (S9): kept. S8's two records of one base stand on the stretch's last
        base: that reading cannot be read, and S8 leaves the row. Up to the gates all
        three left it."""
        run = _gene_run_del()
        far = (HI + 9, GENOME[HI + 8], GENOME[HI + 8] + "TT")
        indels = [(run, {"S7": "1|0", "S8": "1|0", "S9": "1|0"}), (far, {"S9": "0/1"})]
        indels += [(rec, {"S7": "1|0"}) for rec in _two_records_at(HI + 5)]
        indels += [(rec, {"S8": "1|0"}) for rec in _two_records_at(HI + 1)]
        indels.sort(key=lambda x: x[0])
        got, fx = self._carriers({"S7": "1|0", "S8": "1|0", "S9": "1|0", "S6": "1|0"},
                                 indels)
        self.assertEqual(got, {"S6", "S7", "S9"})
        self.assertEqual(fx.ns["_snp_indel_stats"]["kept"], 2)

    def test_junction_outside_and_artefact(self):
        art = (SNP_POS1, "C", "CAA")  # an insertion at the SNP's own position
        got, _fx = self._carriers({"S9": "1|0", "S10": "1|0"},
                                  [(_edge_del(), {"S9": "1|0"}), (art, {"S10": "1|0"})])
        self.assertEqual(got, {"S9", "S10"})

    def test_population_level_unchanged(self):
        got, _fx = self._carriers({"S1": "1|0", "S2": "1|0"},
                                  [(_ins_inside(), {"S1": "1|0", "S2": "1|1"})], fast=True)
        self.assertEqual(got, {"S1", "S2"})

    def test_an_mnv_in_the_target_takes_the_carrier_off(self):
        """2.10.1, the store's MNV file: a carrier whose haplotype holding the row's SNP
        also holds an MNV inside the target (the 10th and 11th bases; and on the first two,
        where an indel row's slot read "from the first base on" started past the MNV up to
        2.10.0) does not have the row's sequence: S1 (cis) and S3 (homozygous) leave the
        row, S2 (the MNV on its other slot) and S6 (no MNV) stay. A population-level
        search reads no store: every carrier stays. At c64e62e1 no MNV was read: all
        stayed (the plan's 3,334 SNP-row listings on HPRC)."""
        snps = {"S1": "1|0", "S2": "1|0", "S3": "1|0", "S6": "1|0"}
        for mnv in (_mnv_inside(), _mnv_first()):
            with self.subTest(mnv=mnv):
                mnvs = [(mnv, {"S1": "1|0", "S2": "0|1", "S3": "1|1"})]
                got, fx = self._carriers(snps, [], mnvs=mnvs)
                self.assertEqual(got, {"S2", "S6"})
                st = fx.ns["_snp_indel_stats"]
                self.assertEqual((st["targets"], st["listings"], st["gone"]), (1, 2, 0))
                got, _fx = self._carriers(snps, [], mnvs=mnvs, fast=True)
                self.assertEqual(got, {"S1", "S2", "S3", "S6"})

    def test_an_mnv_is_read_base_by_base_against_the_rows_snp(self):
        """2.10.1, rule 2. An MNV at the position of the SNP the row names is not
        another allele of that base, as an indel there is: it replaces that base and as
        many after it, so the slot is read base by base. Here it names the row's own alt
        at that base and changes the next one, inside the target: S1 leaves the row, S2
        (the MNV on its other slot) stays. Up to 2.10.0 such a record was skipped as a
        phasing artefact and S1 stayed."""
        got, _fx = self._carriers({"S1": "1|0", "S2": "1|0"}, [],
                                  mnvs=[(_mnv_at_the_rows_snp(),
                                         {"S1": "1|0", "S2": "0|1"})])
        self.assertEqual(got, {"S2"})

    def test_an_mnv_that_keeps_the_base_the_rows_snp_stands_on(self):
        """2.10.1, rule 2. On a base an MNV keeps, the SNV's base wins, so an MNV that
        keeps every base of the target changes none of them and its carrier stays. Up to
        2.10.0 a SNP alt anywhere on a record's REF left the slot with no one sequence
        (a VCF would write a haplotype holding both as one record) and the carrier was
        dropped; that rule still holds for an indel, and for a true third-base
        conflict."""
        got, _fx = self._carriers({"S1": "1|0", "S6": "1|0"}, [],
                                  mnvs=[(_mnv_keeping_the_rows_snp(), {"S1": "1|0"})])
        self.assertEqual(got, {"S1", "S6"})

    def test_an_mnv_past_the_target_or_unphased(self):
        """An MNV past the target changes none of its bases: its carrier (S1) stays; an
        unphased heterozygous MNV call (S4) stays, PUTATIVE, as an unphased indel's; S5's
        phased one in the target takes it off."""
        past = (HI + 3, GENOME[HI + 2:HI + 4],
                GENOME[HI + 2:HI + 4].translate(str.maketrans("ACGT", "TGCA")))
        got, _fx = self._carriers({"S1": "1|0", "S4": "1|0", "S5": "1|0"}, [],
                                  mnvs=[(past, {"S1": "1|0"}),
                                        (_mnv_inside(), {"S4": "0/1", "S5": "1|0"})])
        self.assertEqual(got, {"S1", "S4"})

    def test_no_store_unchanged(self):
        fx = _Fixture({"S1": "1|0"}, [(_ins_inside(), {"S1": "1|0"})],
                      observed=self.observed)
        fx.ns["_snp_indels"] = None
        self.assertEqual(fx.snp_carriers(fx.run()), {"S1"})


class DictionaryPath(ObservedPath):
    """The per-sample dictionary path (no Tier-1 store: the self-test's index)."""

    observed = False

    def test_unphased_dataset(self):
        got, _fx = self._carriers(
            {"S1": "0/1", "S2": "0/1", "S3": "0/1"},
            [(_ins_inside(), {"S1": "1/1", "S2": "0/1", "S3": "1|0"})], phased=False)
        self.assertEqual(got, {"S2", "S3"})


class DenseWindowBed(unittest.TestCase):
    """The high-variant-density BED lists a dense window's samples with an alt, also one
    whose haplotype an indel took off the window's rows: it is unchanged by the filter.
    (Since 2.10.0 the run lists the window from its raw target before it is decomposed,
    _record_dense_window, and writes the BED at the end, one line a window: hvdr_line.)"""

    def test_bed_keeps_the_dropped_sample(self):
        fx = _Fixture({"S1": "1|0", "S2": "1|0"},
                      [(_ins_inside(), {"S1": "1|0", "S2": "0|1"})])
        fx.ns["IUPAC_CAP"] = 0  # every window is logged
        fx.ns["_record_dense_window"](fx.split)
        rows = fx.run()
        self.assertEqual(fx.snp_carriers(rows), {"S2"})
        lines = [fx.ns["hvdr_line"]("chrT", key, rec)
                 for key, rec in sorted(fx.ns["_hvdr_records"].items())]
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].rstrip("\n").split("\t")[5], "S1,S2")


class RunTimeWiring(unittest.TestCase):
    """The helpers the run uses to find and read the store, and the dictionary keys."""

    def test_store_path_mirrors_the_tiers(self):
        with tempfile.TemporaryDirectory() as d:
            dicts = os.path.join(d, "Dictionaries", "dictionaries_hg38_X")
            store = os.path.join(d, "Dictionaries", "indel_genotypes_hg38_X")
            os.makedirs(dicts)
            os.makedirs(store)
            ns = _load_pure_functions({})
            resolve = ns["_resolve_indel_genotype_path"]
            dict_path = os.path.join(dicts, "my_dict_chr7.json")
            self.assertIsNone(resolve(dict_path, "chr7"))
            path = os.path.join(store, "gt_indel_chr7.tsv.gz")
            _write_store(path, [])
            self.assertEqual(resolve(dict_path, "chr7"), path)

    def test_spans_and_frame_keys(self):
        with tempfile.TemporaryDirectory() as d:
            cluster = os.path.join(d, "c.txt")
            dna = G[:5] + "M" + G[6:] + "TGG"
            with open(cluster, "w") as fh:
                fh.write("header\n")
                fh.write("\t".join(["X", G + "NGG", dna, "chrT", str(T0)]) + "\n")
                fh.write("\t".join(["X", G + "NGG", G + "TGG", "chrT", "200"]) + "\n")
            ns = _load_pure_functions({})
            spans = []
            keys = ns["_collect_needed_dict_keys"](cluster, "chrT", spans)
            self.assertEqual(spans, [(LO, HI)])  # the IUPAC target only
            self.assertEqual(len(keys), HI - LO + 1)
            store = os.path.join(d, "s.tsv.gz")
            _write_store(store, [(_ins_inside(), {"S1": "1|0"})])
            st = ns["_SnpRowIndels"](store, spans, 500)
            extra = ns["_snp_indel_frame_keys"](st, spans, "chrT")
            first, last = isc.snp_target_frame(LO, HI, [_ins_inside()])
            self.assertTrue({"chrT,%d" % (p + 1) for p in range(first, last + 1)} <= extra)
            _write_store(store, [(_edge_del(), {"S1": "1|0"})])  # changes nothing
            st = ns["_SnpRowIndels"](store, spans, 500)
            self.assertEqual(ns["_snp_indel_frame_keys"](st, spans, "chrT"), set())

    def test_a_span_many_alignments_share_is_read_once(self):
        """A dense window has a hundred raw targets or more, one per alignment, over the same
        span: the store is read around each span once (2.10.0; reading it once per raw target
        took 227 s on chr22 and an hour on IGH), with the same keys."""
        with tempfile.TemporaryDirectory() as d:
            cluster = os.path.join(d, "c.txt")
            dna = G[:5] + "M" + G[6:] + "TGG"
            low = dna[:2] + dna[2].lower() + dna[3:]
            with open(cluster, "w") as fh:
                fh.write("header\n")
                for _ in range(3):
                    for window in (dna, low):
                        fh.write("\t".join(["X", G + "NGG", window, "chrT", str(T0)]) + "\n")
            ns = _load_pure_functions({})
            spans = []
            keys = ns["_collect_needed_dict_keys"](cluster, "chrT", spans)
            self.assertEqual(spans, [(LO, HI)])
            self.assertEqual(len(keys), HI - LO + 1)


class OptOutGate(unittest.TestCase):
    """_read_snp_row_indels, which the run calls before the dictionary: it reads the store
    in a per-sample search only. CRISPRME_INDEL_SNP=0, the opt-out of the SNP+indel
    co-occurrence, turns this off too (None: no SNP row loses a carrier, as in 2.8.1), as
    does a population-level search. The function sits past the harness's cut (it runs in
    the module body), so it is taken from the source by name."""

    def _reader(self, fast):
        src = os.path.join(PP, "new_simple_analysis.py")
        with open(src) as fh:
            tree = ast.parse(fh.read())
        node = [n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "_read_snp_row_indels"]
        self.assertEqual(len(node), 1)
        ns = _load_pure_functions({})
        exec(compile(ast.Module(body=node, type_ignores=[]), src, "exec"), ns)
        ns.update(current_chr="chrT", _FAST_MODE=fast)
        return ns["_read_snp_row_indels"]

    def _read(self, env, fast=False):
        with tempfile.TemporaryDirectory() as d:
            dicts = os.path.join(d, "Dictionaries", "dictionaries_hg38_X")
            store = os.path.join(d, "Dictionaries", "indel_genotypes_hg38_X")
            os.makedirs(dicts)
            os.makedirs(store)
            _write_store(os.path.join(store, "gt_indel_chrT.tsv.gz"),
                         [(_ins_inside(), {"S1": "1|0"})])
            argv = ["new_simple_analysis.py", "chrT.fa",
                    os.path.join(dicts, "my_dict_chrT.json")]
            keys = set()
            with mock.patch.dict(os.environ, env, clear=False), \
                    mock.patch.object(sys, "argv", argv), \
                    contextlib.redirect_stdout(io.StringIO()):
                if "CRISPRME_INDEL_SNP" not in env:
                    os.environ.pop("CRISPRME_INDEL_SNP", None)
                got = self._reader(fast)([(LO, HI)], keys)
            return got, keys

    def test_read_by_default(self):
        got, keys = self._read({})
        self.assertIsNotNone(got)
        self.assertEqual(len(got), 1)
        self.assertTrue(keys)  # the frame around the target the insertion changes

    def test_on_values(self):
        for value in ("1", "true", "True", "yes"):
            with self.subTest(CRISPRME_INDEL_SNP=value):
                self.assertIsNotNone(self._read({"CRISPRME_INDEL_SNP": value})[0])

    def test_opt_out_turns_it_off(self):
        for value in ("0", "false", "no", ""):
            with self.subTest(CRISPRME_INDEL_SNP=value):
                got, keys = self._read({"CRISPRME_INDEL_SNP": value})
                self.assertIsNone(got)
                self.assertEqual(keys, set())

    def test_population_level_does_not_read_it(self):
        got, keys = self._read({}, fast=True)
        self.assertIsNone(got)
        self.assertEqual(keys, set())

    def test_an_mnv_file_that_cannot_be_read_is_said_on_stdout(self):
        """The review's item 4: a per-sample run whose store's MNV file is cut short reads
        the store's indels (the frame around the target the insertion changes is read) and
        writes one NOTE on stdout, never stderr, that opens as a store of 2.10.0's does."""
        with tempfile.TemporaryDirectory() as d:
            dicts = os.path.join(d, "Dictionaries", "dictionaries_hg38_X")
            store = os.path.join(d, "Dictionaries", "indel_genotypes_hg38_X")
            os.makedirs(dicts)
            os.makedirs(store)
            _write_store(os.path.join(store, "gt_indel_chrT.tsv.gz"),
                         [(_ins_inside(), {"S1": "1|0"})])
            with open(os.path.join(store, "gt_mnv_chrT.tsv.gz"), "wb") as fh:
                fh.write(gzip.compress(b"41_AC_GT\tS2:1|0\n" * 50)[:30])
            argv = ["new_simple_analysis.py", "chrT.fa", os.path.join(dicts, "my_dict_chrT.json")]
            keys, out, err = set(), io.StringIO(), io.StringIO()
            with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(out), \
                    contextlib.redirect_stderr(err):
                os.environ.pop("CRISPRME_INDEL_SNP", None)
                got = self._reader(False)([(LO, HI)], keys)
        self.assertEqual(len(got), 1)
        self.assertTrue(keys)
        self.assertEqual(err.getvalue(), "")
        self.assertEqual(out.getvalue().count(
            "NOTE [indel-snp]: the phased indel genotype store of hg38_X on chrT holds no MNV "
            "that can be read: its MNV file gt_mnv_chrT.tsv.gz is corrupt or cut short ("), 1,
            out.getvalue())
        self.assertIn("no per-sample SNP row loses a carrier for an MNV inside its target",
                      out.getvalue())


class Constants(unittest.TestCase):
    def test_verdicts_match_the_enumerator(self):
        ns = _load_pure_functions({})
        self.assertEqual((ns["_SLOT_KEEP"], ns["_SLOT_DROP"], ns["_SLOT_PUTATIVE"]),
                         (oh.KEEP, oh.DROP, oh.PUTATIVE))


if __name__ == "__main__":
    unittest.main()

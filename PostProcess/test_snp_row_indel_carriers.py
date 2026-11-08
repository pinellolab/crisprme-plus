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


class _Fixture(object):
    """A SNP target at T0 (+ strand) whose SNP C>A makes it perfect, a phased indel
    store, and the harness globals of one per-sample (or population-level) run."""

    def __init__(self, snp_gts, indels, observed=True, phased=True, fast=False):
        self.tmp = tempfile.TemporaryDirectory()
        store = os.path.join(self.tmp.name, "gt_indel_chrT.tsv.gz")
        _write_store(store, indels)
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

    def test_junction_outside_and_artefact(self):
        art = (SNP_POS1, "C", "CAA")  # an insertion at the SNP's own position
        got, _fx = self._carriers({"S9": "1|0", "S10": "1|0"},
                                  [(_edge_del(), {"S9": "1|0"}), (art, {"S10": "1|0"})])
        self.assertEqual(got, {"S9", "S10"})

    def test_population_level_unchanged(self):
        got, _fx = self._carriers({"S1": "1|0", "S2": "1|0"},
                                  [(_ins_inside(), {"S1": "1|0", "S2": "1|1"})], fast=True)
        self.assertEqual(got, {"S1", "S2"})

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
    whose haplotype an indel took off the window's rows: it is unchanged by the filter."""

    def test_bed_keeps_the_dropped_sample(self):
        fx = _Fixture({"S1": "1|0", "S2": "1|0"},
                      [(_ins_inside(), {"S1": "1|0", "S2": "0|1"})])
        fx.ns["IUPAC_CAP"] = 0  # every window is logged
        bed = fx.ns["hvdr_bed"]
        rows = fx.run()
        self.assertEqual(fx.snp_carriers(rows), {"S2"})
        self.assertEqual(bed.getvalue().rstrip("\n").split("\t")[5], "S1,S2")


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


class Constants(unittest.TestCase):
    def test_verdicts_match_the_enumerator(self):
        ns = _load_pure_functions({})
        self.assertEqual((ns["_SLOT_KEEP"], ns["_SLOT_DROP"], ns["_SLOT_PUTATIVE"]),
                         (oh.KEEP, oh.DROP, oh.PUTATIVE))


if __name__ == "__main__":
    unittest.main()

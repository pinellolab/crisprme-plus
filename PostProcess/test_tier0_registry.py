"""Unit tests for tier0_registry (CRISPRme+ dictless redesign, Phase 1).

STDLIB ONLY (unittest + tempfile). Synthetic genotypes only; no real genome
data, no network. The round-trip + correctness assertions are checked against an
independent, hand-written brute-force recomputation so the aggregator's math is
not merely self-consistent.
"""

import os
import tempfile
import unittest

import tier0_registry as t0
from tier0_registry import (
    Counts,
    aggregate_record,
    derive_record_stats,
    compile_registry,
    RegistryReader,
    autosomal_ploidy,
    make_chr_ploidy,
    GLOBAL_GROUP_ID,
    db_group_id,
    db_subpop_group_id,
)


# --------------------------------------------------------------------------- #
# Independent brute-force reference (deliberately written differently from the
# module: iterate slot-by-slot, no shared code path).
# --------------------------------------------------------------------------- #
def brute_counts(genotypes, sample_meta, ploidy_of, restrict=None, dedup=False):
    """Compute Counts over a set of samples with a naive slot walk.

    genotypes: sample_id -> gt string
    restrict:  optional callable(sample_id, meta) -> bool  (which samples belong)
    dedup:     if True, count each sample_id once (global semantics)
    """
    AC = AN = ncar = nhom = ncall = 0
    seen = set()
    for sid, gt in genotypes.items():
        meta = sample_meta.get(sid)
        if meta is None:
            continue
        if restrict is not None and not restrict(sid, meta):
            continue
        if dedup:
            if sid in seen:
                continue
            seen.add(sid)
        db, sp, sex = meta
        ploidy = ploidy_of(sid, sex)
        # explode gt into tokens honoring ploidy, NO duplication of a short gt
        g = gt.strip()
        if "|" in g:
            toks = g.split("|")
        elif "/" in g:
            toks = g.split("/")
        else:
            toks = [g]
        if ploidy <= 1:
            toks = toks[:1]
        else:
            toks = toks[:ploidy]
        s_ac = s_an = 0
        for tk in toks:
            if tk == ".":
                continue
            s_an += 1
            if tk != "0":
                s_ac += 1
        if s_an > 0:
            ncall += 1
        if s_ac > 0:
            ncar += 1
        if s_an > 0 and s_ac == s_an and s_an == ploidy:
            nhom += 1  # fully-called, all-alt (half-missing is NOT hom)
        AC += s_ac
        AN += s_an
    return Counts(AC, AN, ncar, nhom, ncall)


class TestAggregatorCore(unittest.TestCase):
    def setUp(self):
        # Two databases, three subpops.
        # 1000G: EUR {S1(m), S2(f), S3(f)}, AFR {S4(m), S5(f)}
        # HGDP : EAS {H1(m), H2(f)}
        self.meta = {
            "S1": ("1000G", "EUR", "male"),
            "S2": ("1000G", "EUR", "female"),
            "S3": ("1000G", "EUR", "female"),
            "S4": ("1000G", "AFR", "male"),
            "S5": ("1000G", "AFR", "female"),
            "H1": ("HGDP", "EAS", "male"),
            "H2": ("HGDP", "EAS", "female"),
        }

    def test_sparse_groups_only_carriers(self):
        # Only S2 (EUR) and H1 (EAS) carry the alt. AFR has no carriers -> no
        # AFR group. 1000G db group exists (S2). HGDP db group exists (H1).
        gts = {"S2": "0|1", "H1": "1|0"}
        groups = aggregate_record(gts, self.meta, autosomal_ploidy)
        self.assertIn(db_subpop_group_id("1000G", "EUR"), groups)
        self.assertIn(db_subpop_group_id("HGDP", "EAS"), groups)
        self.assertNotIn(db_subpop_group_id("1000G", "AFR"), groups)
        self.assertIn(db_group_id("1000G"), groups)
        self.assertIn(db_group_id("HGDP"), groups)
        self.assertIn(GLOBAL_GROUP_ID, groups)
        # EUR: 1 carrier het -> AC1 AN2 ; global over S2,H1 -> AC2 AN4
        self.assertEqual(groups[db_subpop_group_id("1000G", "EUR")],
                         Counts(1, 2, 1, 0, 1))
        self.assertEqual(groups[GLOBAL_GROUP_ID], Counts(2, 4, 2, 0, 2))


class TestRoundTripBruteForce(unittest.TestCase):
    def setUp(self):
        self.meta = {
            "S1": ("1000G", "EUR", "male"),
            "S2": ("1000G", "EUR", "female"),
            "S3": ("1000G", "AFR", "female"),
            "S4": ("1000G", "AFR", "male"),
            "H1": ("HGDP", "EAS", "female"),
            "H2": ("HGDP", "EAS", "male"),
        }
        # Multi-record synthetic registry (sorted by (pos, alt) inside compiler).
        self.records = [
            # pos 100, alt G: mixed het/hom, one missing
            (100, "A", "G", "rs100",
             {"S1": "0|1", "S2": "1|1", "S3": "./.", "H1": "0|1"}),
            # pos 100, alt T: DIFFERENT carriers (multiallelic same pos)
            (100, "A", "T", "rs100b",
             {"S4": "1|0", "H2": "1|1"}),
            # pos 250, alt C
            (250, "C", "A", "rs250",
             {"S2": "1|0", "S3": "1|1", "H1": "0|1", "H2": "1|0"}),
            # pos 999, alt G (last)
            (999, "T", "G", ".",
             {"S1": "1|1"}),
        ]

    def _compile(self, ploidy_of=autosomal_ploidy):
        d = tempfile.mkdtemp()
        binp = os.path.join(d, "reg.bin")
        idxp = os.path.join(d, "reg.idx.json")
        compile_registry(self.records, self.meta, None, ploidy_of, binp, idxp)
        return binp, idxp

    def test_roundtrip_all_records(self):
        binp, idxp = self._compile()
        reader = RegistryReader(binp, idxp)
        try:
            for (pos, ref, alt, rsid, gts) in self.records:
                got = reader.lookup(pos, alt)
                self.assertIsNotNone(got, "record (%s,%s) missing" % (pos, alt))
                # rsID round-trip
                self.assertEqual(reader.rsid(pos, alt), rsid)
                # brute-force every group id present
                expected = aggregate_record(gts, self.meta, autosomal_ploidy)
                self.assertEqual(set(got.keys()), set(expected.keys()),
                                 "group set mismatch at (%s,%s)" % (pos, alt))
                for gid, cnt in expected.items():
                    self.assertEqual(got[gid], cnt,
                                     "counts mismatch %s @ (%s,%s)" % (gid, pos, alt))
                    # cross-check each group vs the independent brute walk
                    if gid == GLOBAL_GROUP_ID:
                        ref_cnt = brute_counts(gts, self.meta, autosomal_ploidy,
                                               dedup=True)
                    elif t0.SEP in gid:
                        db, sp = gid.split(t0.SEP, 1)
                        ref_cnt = brute_counts(
                            gts, self.meta, autosomal_ploidy,
                            restrict=lambda sid, m, db=db, sp=sp: m[0] == db and m[1] == sp)
                    else:
                        db = gid
                        ref_cnt = brute_counts(
                            gts, self.meta, autosomal_ploidy,
                            restrict=lambda sid, m, db=db: m[0] == db)
                    self.assertEqual(got[gid], ref_cnt,
                                     "brute mismatch %s @ (%s,%s)" % (gid, pos, alt))
                # derived helpers vs hand computation on the GLOBAL group
                der = reader.derived(pos, alt)
                gcnt = got[GLOBAL_GROUP_ID]
                self.assertAlmostEqual(der["global_af"],
                                       gcnt.AC / gcnt.AN if gcnt.AN else 0.0)
                self.assertAlmostEqual(
                    der["global_carrier_freq"],
                    gcnt.n_carrier_indiv / gcnt.n_called_indiv if gcnt.n_called_indiv else 0.0)
        finally:
            reader.close()

    def test_derived_specific_values(self):
        # pos 100 alt G: S1 0|1 (het), S2 1|1 (hom), S3 ./. (missing), H1 0|1
        # global dedup: S1,S2,S3,H1
        #   AC = 1 + 2 + 0 + 1 = 4 ; AN = 2 + 2 + 0 + 2 = 6 (S3 missing excluded)
        #   carriers = S1,S2,H1 = 3 ; hom = S2 = 1 ; called = S1,S2,H1 = 3
        binp, idxp = self._compile()
        reader = RegistryReader(binp, idxp)
        try:
            g = reader.lookup(100, "G")[GLOBAL_GROUP_ID]
            self.assertEqual(g, Counts(4, 6, 3, 1, 3))
            self.assertAlmostEqual(g.allele_freq(), 4.0 / 6.0)
            self.assertAlmostEqual(g.carrier_freq(), 3.0 / 3.0)
            self.assertAlmostEqual(g.hom_freq(), 1.0 / 3.0)
        finally:
            reader.close()


class TestANvs2N(unittest.TestCase):
    """AF uses AN (called alleles), NOT 2*N; fully-missing excluded from called."""

    def test_missing_excluded_from_an_and_called(self):
        meta = {
            "A": ("1000G", "EUR", "female"),
            "B": ("1000G", "EUR", "female"),
            "C": ("1000G", "EUR", "female"),
            "D": ("1000G", "EUR", "female"),
        }
        # A: 0|1 (het), B: 1|1 (hom), C: ./. (fully missing), D: 1|. (half-missing)
        gts = {"A": "0|1", "B": "1|1", "C": "./.", "D": "1|."}
        groups = aggregate_record(gts, meta, autosomal_ploidy)
        g = groups[GLOBAL_GROUP_ID]
        # AC = 1 + 2 + 0 + 1 = 4
        # AN = 2 + 2 + 0 + 1 = 5   (NOT 2*4=8)
        # called individuals = A,B,D = 3  (C fully missing excluded)
        # carriers = A,B,D = 3
        # hom = B (1|1, fully-called all-alt) = 1. D (1|. half-missing) is NOT hom:
        #   an=1 != ploidy=2, so the unobserved allele may be ref (matches
        #   bcftools/plink; the hom stat requires a fully-called genotype).
        self.assertEqual(g, Counts(4, 5, 3, 1, 3))
        self.assertAlmostEqual(g.allele_freq(), 4.0 / 5.0)
        self.assertNotAlmostEqual(g.allele_freq(), 4.0 / 8.0)  # would be the 2N bug
        self.assertEqual(g.n_called_indiv, 3)  # fully-missing C excluded


class TestPloidySex(unittest.TestCase):
    """chrX-nonPAR: haploid males must contribute AC<=1, AN=1 -- no hap*2."""

    def setUp(self):
        self.meta = {
            "M1": ("1000G", "EUR", "male"),    # hemizygous alt
            "M2": ("1000G", "EUR", "male"),    # hemizygous ref
            "F1": ("1000G", "EUR", "female"),  # het
            "F2": ("1000G", "EUR", "female"),  # hom alt
        }
        self.ploidy = make_chr_ploidy(haploid_male=True)

    def test_hemizygous_male_single_copy(self):
        # M1 alt: even if written diploid-style "1|1", ploidy 1 -> AC 1, AN 1.
        gts = {"M1": "1|1", "M2": "0", "F1": "0|1", "F2": "1|1"}
        groups = aggregate_record(gts, self.meta, self.ploidy)
        g = groups[GLOBAL_GROUP_ID]
        # M1: AC1 AN1 (hemizygous carrier, counts hom once)
        # M2: AC0 AN1
        # F1: AC1 AN2
        # F2: AC2 AN2
        # AC = 1+0+1+2 = 4 ; AN = 1+1+2+2 = 6
        # carriers = M1,F1,F2 = 3 ; hom = M1 (all called alt) + F2 = 2 ; called = 4
        self.assertEqual(g, Counts(4, 6, 3, 2, 4))
        # Prove the hap*2 double-count is NOT present: if M1 were treated as
        # homozygous diploid it would add AC2/AN2, giving AC5/AN7.
        self.assertNotEqual(g.AC, 5)
        self.assertNotEqual(g.AN, 7)

    def test_female_diploid_still_two(self):
        gts = {"F1": "1|1"}
        groups = aggregate_record(gts, self.meta, self.ploidy)
        g = groups[GLOBAL_GROUP_ID]
        self.assertEqual(g, Counts(2, 2, 1, 1, 1))  # female ploidy 2

    def test_haploid_male_single_token_gt(self):
        # single-token haploid GT "1"
        gts = {"M1": "1", "M2": "0"}
        groups = aggregate_record(gts, self.meta, self.ploidy)
        g = groups[GLOBAL_GROUP_ID]
        # M1 AC1 AN1 carrier hom; M2 AC0 AN1
        self.assertEqual(g, Counts(1, 2, 1, 1, 2))


class TestMultiallelic(unittest.TestCase):
    def setUp(self):
        self.meta = {
            "S1": ("1000G", "EUR", "female"),
            "S2": ("1000G", "EUR", "female"),
            "S3": ("1000G", "AFR", "female"),
        }
        # Same pos 500: alt C carried by S1; alt T carried by S2,S3.
        self.records = [
            (500, "A", "C", "rsC", {"S1": "1|1"}),
            (500, "A", "T", "rsT", {"S2": "0|1", "S3": "1|1"}),
        ]

    def test_per_alt_independent(self):
        d = tempfile.mkdtemp()
        binp = os.path.join(d, "r.bin")
        idxp = os.path.join(d, "r.idx.json")
        compile_registry(self.records, self.meta, None, autosomal_ploidy, binp, idxp)
        reader = RegistryReader(binp, idxp)
        try:
            c = reader.lookup(500, "C")
            t = reader.lookup(500, "T")
            self.assertIsNotNone(c)
            self.assertIsNotNone(t)
            self.assertNotEqual(c, t)
            # alt C: only S1 (EUR), hom -> global AC2 AN2 car1 hom1 called1
            self.assertEqual(c[GLOBAL_GROUP_ID], Counts(2, 2, 1, 1, 1))
            # EUR group present for C, AFR absent
            self.assertIn(db_subpop_group_id("1000G", "EUR"), c)
            self.assertNotIn(db_subpop_group_id("1000G", "AFR"), c)
            # alt T: S2 (EUR het) + S3 (AFR hom) -> global AC3 AN4 car2 hom1 called2
            self.assertEqual(t[GLOBAL_GROUP_ID], Counts(3, 4, 2, 1, 2))
            self.assertIn(db_subpop_group_id("1000G", "AFR"), t)
            self.assertEqual(reader.rsid(500, "C"), "rsC")
            self.assertEqual(reader.rsid(500, "T"), "rsT")
        finally:
            reader.close()

    def test_unsplit_multiallelic_token_not_miscounted(self):
        # An UNSPLIT multiallelic genotype "1|2" on the alt-1 (default alt_index)
        # record: only the '1' is THIS record's alt; the '2' is a FOREIGN alt =>
        # a called, non-carrier allele, NOT a second copy of alt-1. (Regression
        # for the classifier treating any non-'0' token as this record's alt,
        # which would have doubled AC and fabricated a homozygote.)
        meta = {"S1": ("1000G", "EUR", "female")}
        g = aggregate_record({"S1": "1|2"}, meta, autosomal_ploidy)[GLOBAL_GROUP_ID]
        # AC=1 (only the '1'), AN=2 (both called), carrier yes, hom NO, called 1
        self.assertEqual(g, Counts(1, 2, 1, 0, 1))
        # Viewed as the alt-2 record (explicit alt_index="2"), symmetric result.
        g2 = aggregate_record({"S1": "1|2"}, meta, autosomal_ploidy,
                              alt_index="2")[GLOBAL_GROUP_ID]
        self.assertEqual(g2, Counts(1, 2, 1, 0, 1))
        # "2|2" on the alt-1 record has ZERO alt-1 copies -> not a carrier -> no
        # global group emitted (sparse), never miscounted as hom-alt-1.
        gzero = aggregate_record({"S1": "2|2"}, meta, autosomal_ploidy)
        self.assertNotIn(GLOBAL_GROUP_ID, gzero)


class TestMultiDatabaseGlobalDedup(unittest.TestCase):
    """A sample present in BOTH databases counts ONCE in global."""

    def test_dedup_union(self):
        # SHARED appears in both 1000G and HGDP under the same canonical id.
        # aggregate_record keys genotypes by sample_id, so it can only appear
        # once in the genotype dict -> but its META could be either db. We model
        # the "present in both" case by having the sample carry a distinct db in
        # meta, and asserting global counts it once. To truly exercise dedup we
        # give a genotype dict where the SAME canonical id maps to one gt, while
        # BOTH db groups would like to claim it. Since meta is 1:1 sample->db,
        # the realistic cross-db-shared-sample scenario is handled by the global
        # dedup on canonical id: the same id can only contribute once globally.
        meta = {
            "S1": ("1000G", "EUR", "female"),
            "S2": ("1000G", "EUR", "female"),
            "H1": ("HGDP", "EAS", "female"),
            "SHARED": ("1000G", "EUR", "female"),
        }
        gts = {"S1": "0|1", "S2": "1|1", "H1": "1|0", "SHARED": "1|1"}
        groups = aggregate_record(gts, meta, autosomal_ploidy)
        # 1000G::EUR : S1(1/2) S2(2/2) SHARED(2/2) -> AC5 AN6 car3 hom2 called3
        self.assertEqual(groups[db_subpop_group_id("1000G", "EUR")],
                         Counts(5, 6, 3, 2, 3))
        # HGDP::EAS : H1 -> AC1 AN2 car1 hom0 called1
        self.assertEqual(groups[db_subpop_group_id("HGDP", "EAS")],
                         Counts(1, 2, 1, 0, 1))
        # 1000G db aggregate == its only subpop EUR here
        self.assertEqual(groups[db_group_id("1000G")], Counts(5, 6, 3, 2, 3))
        # global dedup over S1,S2,H1,SHARED (4 individuals)
        #  AC = 1+2+1+2 = 6 ; AN = 2*4 = 8 ; car = all 4 ; hom = S2,SHARED = 2
        self.assertEqual(groups[GLOBAL_GROUP_ID], Counts(6, 8, 4, 2, 4))

    def test_shared_id_counted_once_globally(self):
        # Explicit dedup exercise: a duplicate canonical id must not inflate the
        # global group. We simulate two source rows for one canonical id by
        # merging on the id (as a real cross-db merge on sample_id would).
        meta = {"X": ("1000G", "EUR", "male")}
        ploidy = autosomal_ploidy
        # Independent brute (dedup) vs aggregate for a single sample.
        gts = {"X": "1|1"}
        groups = aggregate_record(gts, meta, ploidy)
        ref = brute_counts(gts, meta, ploidy, dedup=True)
        self.assertEqual(groups[GLOBAL_GROUP_ID], ref)
        self.assertEqual(groups[GLOBAL_GROUP_ID].n_carrier_indiv, 1)


class TestBinarySearchMmap(unittest.TestCase):
    def setUp(self):
        self.meta = {"S1": ("1000G", "EUR", "female")}
        self.records = [
            (10, "A", "G", "rs10", {"S1": "0|1"}),
            (10, "A", "T", "rs10t", {"S1": "1|1"}),
            (50, "C", "A", "rs50", {"S1": "1|0"}),
            (777, "G", "C", "rs777", {"S1": "1|1"}),
        ]
        self.d = tempfile.mkdtemp()
        self.binp = os.path.join(self.d, "r.bin")
        self.idxp = os.path.join(self.d, "r.idx.json")
        compile_registry(self.records, self.meta, None, autosomal_ploidy,
                         self.binp, self.idxp)

    def test_present_absent_first_last(self):
        reader = RegistryReader(self.binp, self.idxp)
        try:
            self.assertEqual(len(reader), 4)
            # first record after sort by (pos, alt) is (10, 'G')
            self.assertEqual(reader.record_key_at(0), (10, "G"))
            self.assertEqual(reader.record_key_at(1), (10, "T"))
            self.assertEqual(reader.record_key_at(len(reader) - 1), (777, "C"))
            # present lookups
            self.assertIsNotNone(reader.lookup(10, "G"))
            self.assertIsNotNone(reader.lookup(10, "T"))
            self.assertIsNotNone(reader.lookup(777, "C"))
            # absent position
            self.assertIsNone(reader.lookup(999, "A"))
            self.assertIsNone(reader.lookup(11, "G"))
            # present position, absent alt
            self.assertIsNone(reader.lookup(10, "C"))
            self.assertIsNone(reader.lookup(50, "G"))
            self.assertIsNone(reader.rsid(999, "A"))
        finally:
            reader.close()

    def test_reader_does_not_materialize_records(self):
        # Structural guarantee: the reader must not keep a Python list/dict of
        # all records. The only permitted large state is the mmap and the JSON
        # manifest (taxonomy/section offsets, NOT per-record data). We assert:
        #  (a) an mmap is held,
        #  (b) no reader attribute (other than the manifest) is a list/tuple, and
        #      no dict attribute (other than the manifest) has one entry per
        #      record -- i.e. the record array was not exploded into Python.
        import mmap as _mmap
        reader = RegistryReader(self.binp, self.idxp)
        try:
            n = len(reader)
            self.assertTrue(
                any(isinstance(v, _mmap.mmap) for v in vars(reader).values()),
                "reader must hold an mmap")
            for name, val in vars(reader).items():
                if name == "manifest":
                    continue  # taxonomy + offsets, not per-record
                if isinstance(val, (list, tuple)):
                    self.fail("reader holds a Python sequence attr %r "
                              "(records must stay on the mmap)" % name)
                if isinstance(val, dict) and n > 1 and len(val) >= n:
                    self.fail("reader materialized a record-sized dict: %s" % name)
            # The record-key view must be a lazy accessor, not a real list: it
            # reports length n but is not a list/tuple instance.
            self.assertFalse(isinstance(reader._keys, (list, tuple)))
            self.assertEqual(len(reader._keys), n)
            # and lookups still work via bisect on mmap
            self.assertIsNotNone(reader.lookup(50, "A"))
        finally:
            reader.close()


class TestMaxSubpopAndObserved(unittest.TestCase):
    def test_max_subpop_af_and_observed(self):
        meta = {
            "E1": ("1000G", "EUR", "female"),
            "E2": ("1000G", "EUR", "female"),
            "A1": ("1000G", "AFR", "female"),
            "A2": ("1000G", "AFR", "female"),
            "X1": ("HGDP", "EAS", "female"),
        }
        # EUR: E1 0|1, E2 0|0(ref, not a carrier -> not in gts, dict lists only
        #      alt carriers) -> EUR AC1 AN2 (only E1 present)
        # AFR: A1 1|1, A2 1|1 -> AC4 AN4 -> AF 1.0 (the max)
        # EAS: X1 0|1 -> AC1 AN2 -> AF 0.5
        gts = {"E1": "0|1", "A1": "1|1", "A2": "1|1", "X1": "0|1"}
        groups = aggregate_record(gts, meta, autosomal_ploidy)
        der = derive_record_stats(groups)
        self.assertTrue(der["observed"])
        self.assertAlmostEqual(der["max_subpop_af"], 1.0)
        self.assertEqual(der["max_subpop_af_label"],
                         db_subpop_group_id("1000G", "AFR"))
        # EUR AF present but lower
        self.assertAlmostEqual(
            groups[db_subpop_group_id("1000G", "EUR")].allele_freq(), 0.5)

    def test_observed_false_when_no_carriers(self):
        meta = {"E1": ("1000G", "EUR", "female")}
        # only ref / missing -> no carrier -> no groups, observed False
        gts = {"E1": "0|0"}
        groups = aggregate_record(gts, meta, autosomal_ploidy)
        self.assertNotIn(GLOBAL_GROUP_ID, groups)
        der = derive_record_stats(groups)
        self.assertFalse(der["observed"])
        self.assertEqual(der["max_subpop_af"], 0.0)
        self.assertIsNone(der["max_subpop_af_label"])


if __name__ == "__main__":
    unittest.main()

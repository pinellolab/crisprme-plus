"""Golden tests for tier0_compile (CRISPRme+ dictless redesign, Phase 1 step 2).

STDLIB ONLY (unittest + tempfile + json). Hand-written synthetic fixtures: a small
``my_dict_<chrom>.json`` (carriers-only, incl. a multiallelic '$' site) and two
per-db samplesID files (1000G with EUR/AFR both sexes, HGDP with EAS). We compile
into a panel-aware registry, open it with RegistryReader, and assert every
per-(db x subpop) / per-db / global slot -- and the derived AF / carrier_freq /
max_subpop -- EXACTLY match an INDEPENDENT brute force.

The brute force is written a DIFFERENT way from the module: it reconstructs the
FULL PANEL (carriers from the fixture dict + EVERY OTHER samplesID sample as
hom-ref) and walks it slot-by-slot. This proves the full-panel AN correctness: a
1-het-carrier record in a 10-sample panel gives AF = 1/20, NOT 1/2.
"""

import os
import tempfile
import unittest

import tier0_registry as t0
from tier0_registry import (
    Counts,
    RegistryReader,
    autosomal_ploidy,
    GLOBAL_GROUP_ID,
    db_group_id,
    db_subpop_group_id,
    SEP,
)
import tier0_compile as tc


# --------------------------------------------------------------------------- #
# Independent brute-force FULL-PANEL reference.
# --------------------------------------------------------------------------- #
def _explode(gt, ploidy):
    """Tokens of a genotype honoring ploidy, no duplication of a short GT."""
    g = (gt or "").strip()
    if "|" in g:
        toks = g.split("|")
    elif "/" in g:
        toks = g.split("/")
    else:
        toks = [g] if g else []
    if ploidy <= 1:
        return toks[:1]
    return toks[:ploidy]


def brute_full_panel(carrier_gts, sample_meta, ploidy_of, restrict=None,
                     dedup=False, alt_index="1"):
    """Slot-by-slot Counts over the FULL PANEL for one (pos, alt) record.

    EVERY sample in ``sample_meta`` (filtered by ``restrict``) is walked. A sample
    listed in ``carrier_gts`` uses its listed genotype; ANY OTHER panel sample is
    reconstructed as HOM-REF (all-ref, fully called, ploidy slots of "0"). This is
    the independent re-derivation of the full-panel AN that tier0_compile must
    match.

    restrict: optional callable(sid, meta) -> bool (group membership filter).
    dedup:    if True, count each canonical sample_id once (global semantics).
    """
    AC = AN = ncar = nhom = ncall = 0
    seen = set()
    for sid, meta in sample_meta.items():
        if restrict is not None and not restrict(sid, meta):
            continue
        if dedup:
            if sid in seen:
                continue
            seen.add(sid)
        db, sp, sex = meta
        ploidy = ploidy_of(sid, sex)
        if sid in carrier_gts:
            toks = _explode(carrier_gts[sid], ploidy)
        else:
            # Unlisted -> hom-ref: ploidy slots of ref "0", all called.
            toks = ["0"] * ploidy
        s_ac = s_an = 0
        for tk in toks:
            if tk == ".":
                continue
            s_an += 1
            if tk == alt_index:
                s_ac += 1
        if s_an > 0:
            ncall += 1
        if s_ac > 0:
            ncar += 1
        if s_an > 0 and s_ac == s_an and s_an == ploidy:
            nhom += 1
        AC += s_ac
        AN += s_an
    return Counts(AC, AN, ncar, nhom, ncall)


# --------------------------------------------------------------------------- #
# Fixture writers.
# --------------------------------------------------------------------------- #
def write_samplesid(path, rows):
    """rows: list of (sample_id, population, superpopulation, sex)."""
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for sid, pop, spop, sex in rows:
            fh.write("\t".join([sid, pop, spop, sex]) + "\n")


def write_dict(path, mapping):
    """mapping: {"<chrom>,<pos>": "<entry>$<entry>..."} -> json file."""
    import json
    with open(path, "w") as fh:
        json.dump(mapping, fh)


# 1000G panel: EUR (both sexes) + AFR (both sexes); HGDP panel: EAS (both sexes).
G1000_ROWS = [
    ("S_EUR_M1", "GBR", "EUR", "male"),
    ("S_EUR_M2", "GBR", "EUR", "male"),
    ("S_EUR_F1", "GBR", "EUR", "female"),
    ("S_EUR_F2", "GBR", "EUR", "female"),
    ("S_AFR_M1", "YRI", "AFR", "male"),
    ("S_AFR_F1", "YRI", "AFR", "female"),
    ("S_AFR_F2", "YRI", "AFR", "female"),
]
HGDP_ROWS = [
    ("H_EAS_M1", "Han", "EAS", "male"),
    ("H_EAS_F1", "Han", "EAS", "female"),
    ("H_EAS_F2", "Han", "EAS", "female"),
]
# Total panel: 7 (1000G) + 3 (HGDP) = 10 individuals.


class _CompileFixture(unittest.TestCase):
    """Shared setUp: writes samplesID files + a dict, builds sample_meta."""

    CHROM = "chr1"

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.g1000_sid = os.path.join(self.d, "hg38_1000G.samplesID.txt")
        self.hgdp_sid = os.path.join(self.d, "hg38_HGDP.samplesID.txt")
        write_samplesid(self.g1000_sid, G1000_ROWS)
        write_samplesid(self.hgdp_sid, HGDP_ROWS)
        self.db_to_sid = {"1000G": self.g1000_sid, "HGDP": self.hgdp_sid}
        # Independent sample_meta reconstruction for the brute force (do NOT reuse
        # tier0_compile.build_sample_meta here -- build it by hand).
        self.meta = {}
        for sid, pop, spop, sex in G1000_ROWS:
            self.meta[sid] = ("1000G", spop, sex)
        for sid, pop, spop, sex in HGDP_ROWS:
            self.meta[sid] = ("HGDP", spop, sex)

    def _dict_path(self):
        return os.path.join(self.d, "my_dict_%s.json" % self.CHROM)

    def _compile(self, dict_mapping):
        dpath = self._dict_path()
        write_dict(dpath, dict_mapping)
        binp = os.path.join(self.d, "reg.bin")
        idxp = os.path.join(self.d, "reg.idx.json")
        stats = tc.compile_from_dict(dpath, self.db_to_sid, self.CHROM, binp, idxp)
        return binp, idxp, stats

    def _ploidy_of(self):
        return tc.ploidy_of_for_chrom(self.CHROM)

    def _assert_all_groups(self, reader, pos, alt, carrier_gts):
        """Assert the registry's groups for (pos,alt) match the FULL-PANEL brute
        force slot-by-slot, and that groups are sparse-on-carrier."""
        ploidy_of = self._ploidy_of()
        got = reader.lookup(pos, alt)
        self.assertIsNotNone(got, "record (%s,%s) missing" % (pos, alt))

        # Expected group ids: those (db x subpop) / db / global groups that have
        # >= 1 carrier in carrier_gts. Compute independently.
        expected_gids = set()
        # collect all (db, subpop) and db present in the panel
        subpops = set((m[0], m[1]) for m in self.meta.values())
        dbs = set(m[0] for m in self.meta.values())
        for (db, sp) in subpops:
            ref = brute_full_panel(
                carrier_gts, self.meta, ploidy_of,
                restrict=lambda sid, m, db=db, sp=sp: m[0] == db and m[1] == sp)
            if ref.n_carrier_indiv >= 1:
                expected_gids.add(db_subpop_group_id(db, sp))
        for db in dbs:
            ref = brute_full_panel(
                carrier_gts, self.meta, ploidy_of,
                restrict=lambda sid, m, db=db: m[0] == db)
            if ref.n_carrier_indiv >= 1:
                expected_gids.add(db_group_id(db))
        gref = brute_full_panel(carrier_gts, self.meta, ploidy_of, dedup=True)
        if gref.n_carrier_indiv >= 1:
            expected_gids.add(GLOBAL_GROUP_ID)

        self.assertEqual(set(got.keys()), expected_gids,
                         "group set mismatch at (%s,%s)" % (pos, alt))

        # Slot-by-slot equality for every present group.
        for gid, cnt in got.items():
            if gid == GLOBAL_GROUP_ID:
                ref = brute_full_panel(carrier_gts, self.meta, ploidy_of,
                                       dedup=True)
            elif SEP in gid:
                db, sp = gid.split(SEP, 1)
                ref = brute_full_panel(
                    carrier_gts, self.meta, ploidy_of,
                    restrict=lambda sid, m, db=db, sp=sp: m[0] == db and m[1] == sp)
            else:
                db = gid
                ref = brute_full_panel(
                    carrier_gts, self.meta, ploidy_of,
                    restrict=lambda sid, m, db=db: m[0] == db)
            self.assertEqual(cnt, ref,
                             "counts mismatch %s @ (%s,%s): got %r want %r"
                             % (gid, pos, alt, cnt, ref))
        return got


# --------------------------------------------------------------------------- #
# (1) FULL-PANEL AN: AF denominator is the whole panel, not carriers-only.
# --------------------------------------------------------------------------- #
class TestFullPanelAN(_CompileFixture):
    def test_one_het_carrier_af_is_one_over_2N(self):
        # A single het carrier (S_EUR_F1) at pos 100. Panel is 10 individuals.
        # Global AN = 2 * 10 = 20 (all diploid autosome). AC = 1. AF = 1/20.
        # NOT 1/2 (which is the carriers-only bug).
        carrier_gts = {"S_EUR_F1": "0|1"}
        mapping = {"%s,100" % self.CHROM: "S_EUR_F1:0|1;A,G;rs100;0.05"}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            got = self._assert_all_groups(reader, 100, "G", carrier_gts)
            g = got[GLOBAL_GROUP_ID]
            self.assertEqual(g.AC, 1)
            self.assertEqual(g.AN, 20)             # full panel, NOT 2
            self.assertAlmostEqual(g.allele_freq(), 1.0 / 20.0)
            self.assertNotAlmostEqual(g.allele_freq(), 1.0 / 2.0)  # carriers-only bug
            # n_called = full panel (10), n_carrier = 1
            self.assertEqual(g.n_called_indiv, 10)
            self.assertEqual(g.n_carrier_indiv, 1)
            self.assertAlmostEqual(g.carrier_freq(), 1.0 / 10.0)
            # 1000G::EUR: 4 individuals -> AN 8, AC 1 -> AF 1/8
            eur = got[db_subpop_group_id("1000G", "EUR")]
            self.assertEqual((eur.AC, eur.AN, eur.n_called_indiv), (1, 8, 4))
            self.assertAlmostEqual(eur.allele_freq(), 1.0 / 8.0)
            # sparse: AFR / EAS / HGDP have no carriers -> absent
            self.assertNotIn(db_subpop_group_id("1000G", "AFR"), got)
            self.assertNotIn(db_subpop_group_id("HGDP", "EAS"), got)
            self.assertNotIn(db_group_id("HGDP"), got)
            # 1000G db group present (its EUR subpop carries)
            self.assertIn(db_group_id("1000G"), got)
            g1 = got[db_group_id("1000G")]
            self.assertEqual((g1.AC, g1.AN, g1.n_called_indiv), (1, 14, 7))
        finally:
            reader.close()

    def test_stats_counts(self):
        mapping = {"%s,100" % self.CHROM: "S_EUR_F1:0|1;A,G;rs100;0.05"}
        _, _, stats = self._compile(mapping)
        self.assertEqual(stats["n_written"], 1)
        self.assertEqual(stats["n_skipped_indel"], 0)
        self.assertEqual(stats["n_skipped_empty"], 0)
        self.assertEqual(stats["manifest"]["aggregation"], "panel")
        self.assertEqual(stats["overlaps"], [])


# --------------------------------------------------------------------------- #
# (2) Ploidy / sex: chrX hemizygous male -> AC1/AN1; panel AN sums males at 1.
# --------------------------------------------------------------------------- #
class TestChrXPloidy(_CompileFixture):
    CHROM = "chrX"

    def test_hemizygous_male_and_panel_an(self):
        # A hemizygous male carrier (S_EUR_M1) on chrX. Panel: 10 individuals.
        # Males (S_EUR_M1, S_EUR_M2, S_AFR_M1, H_EAS_M1) are haploid (AN 1 each);
        # females diploid (AN 2 each). Panel AN = 4*1 (males) + 6*2 (females) = 16.
        # Male carrier gives AC 1 / (its own AN) 1, and IS hom (fully-called all-alt
        # at ploidy 1).
        carrier_gts = {"S_EUR_M1": "1"}
        mapping = {"chrX,500": "S_EUR_M1:1;A,T;rsX;0.02"}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            got = self._assert_all_groups(reader, 500, "T", carrier_gts)
            g = got[GLOBAL_GROUP_ID]
            self.assertEqual(g.AC, 1)
            self.assertEqual(g.AN, 16)     # 4 males*1 + 6 females*2
            self.assertEqual(g.n_carrier_indiv, 1)
            self.assertEqual(g.n_hom_indiv, 1)  # hemizygous male fully-called alt
            self.assertEqual(g.n_called_indiv, 10)
            self.assertAlmostEqual(g.allele_freq(), 1.0 / 16.0)
            # 1000G::EUR: M1,M2 (haploid) + F1,F2 (diploid) -> AN = 1+1+2+2 = 6
            eur = got[db_subpop_group_id("1000G", "EUR")]
            self.assertEqual((eur.AC, eur.AN, eur.n_called_indiv), (1, 6, 4))
            self.assertAlmostEqual(eur.allele_freq(), 1.0 / 6.0)
        finally:
            reader.close()

    def test_diploid_style_male_gt_not_double_counted(self):
        # Even if the dict wrote the male GT diploid-style "1|1", ploidy 1 means
        # AC 1 / AN 1, NOT AC 2 / AN 2 (anti-double-count invariant).
        carrier_gts = {"S_EUR_M1": "1|1"}
        mapping = {"chrX,600": "S_EUR_M1:1|1;A,T;rsX2;0.02"}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            got = self._assert_all_groups(reader, 600, "T", carrier_gts)
            g = got[GLOBAL_GROUP_ID]
            self.assertEqual(g.AC, 1)   # NOT 2
            self.assertEqual(g.AN, 16)  # panel unchanged; male still ploidy 1
        finally:
            reader.close()


# --------------------------------------------------------------------------- #
# (2b) chrY: females carry NO chrY -> ABSENT (ploidy 0), excluded from the panel
#      so they add no phantom alleles to the chrY AN denominator.
# --------------------------------------------------------------------------- #
class TestChrYPloidy(_CompileFixture):
    CHROM = "chrY"

    def test_females_excluded_from_chrY_panel(self):
        # Panel has 4 males + 6 females. On chrY, males are haploid (AN 1) and
        # females are ABSENT (ploidy 0). So the chrY panel = 4 males, global AN = 4
        # (NOT 10, and NOT 4+6=10 if females were wrongly treated as haploid).
        # A single male carrier -> AC 1, AF = 1/4 (biologically correct), not 1/10.
        carrier_gts = {"S_EUR_M1": "1"}
        mapping = {"chrY,900": "S_EUR_M1:1;A,T;rsY;0.01"}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            got = self._assert_all_groups(reader, 900, "T", carrier_gts)
            g = got[GLOBAL_GROUP_ID]
            self.assertEqual(g.AC, 1)
            self.assertEqual(g.AN, 4)              # 4 males*1; 6 females EXCLUDED
            self.assertNotEqual(g.AN, 10)          # the phantom-female-AN bug
            self.assertEqual(g.n_called_indiv, 4)  # females not called on chrY
            self.assertEqual(g.n_carrier_indiv, 1)
            self.assertAlmostEqual(g.allele_freq(), 1.0 / 4.0)
            # 1000G::EUR chrY panel = M1, M2 only (F1, F2 excluded) -> AN 2
            eur = got[db_subpop_group_id("1000G", "EUR")]
            self.assertEqual((eur.AC, eur.AN, eur.n_called_indiv), (1, 2, 2))
        finally:
            reader.close()


# --------------------------------------------------------------------------- #
# (3) Multiallelic '$' entry: per-(pos, alt) independent.
# --------------------------------------------------------------------------- #
class TestMultiallelic(_CompileFixture):
    def test_dollar_split_per_alt(self):
        # pos 300: alt G carried by S_EUR_F1 (het) + S_AFR_F1 (hom);
        #          alt T carried by H_EAS_M1 (hemi? no, autosome so diploid het).
        gts_G = {"S_EUR_F1": "0|1", "S_AFR_F1": "1|1"}
        gts_T = {"H_EAS_M1": "0|1"}
        entry_G = "S_EUR_F1:0|1,S_AFR_F1:1|1;A,G;rsG;0.1"
        entry_T = "H_EAS_M1:0|1;A,T;rsT;0.05"
        mapping = {"%s,300" % self.CHROM: entry_G + "$" + entry_T}
        binp, idxp, stats = self._compile(mapping)
        self.assertEqual(stats["n_written"], 2)
        reader = RegistryReader(binp, idxp)
        try:
            gotG = self._assert_all_groups(reader, 300, "G", gts_G)
            gotT = self._assert_all_groups(reader, 300, "T", gts_T)
            self.assertNotEqual(gotG, gotT)
            # alt G global: het F1 (AC1) + hom F1 AFR (AC2) -> AC3, AN 20
            gG = gotG[GLOBAL_GROUP_ID]
            self.assertEqual((gG.AC, gG.AN, gG.n_carrier_indiv, gG.n_hom_indiv), (3, 20, 2, 1))
            # alt T global: one het EAS male on autosome -> AC1 AN20
            gT = gotT[GLOBAL_GROUP_ID]
            self.assertEqual((gT.AC, gT.AN, gT.n_carrier_indiv), (1, 20, 1))
            # rsID round-trips per alt
            self.assertEqual(reader.rsid(300, "G"), "rsG")
            self.assertEqual(reader.rsid(300, "T"), "rsT")
            # alt G present in 1000G groups (EUR+AFR carriers), NOT HGDP;
            # alt T present in HGDP::EAS, NOT 1000G.
            self.assertIn(db_subpop_group_id("1000G", "EUR"), gotG)
            self.assertIn(db_subpop_group_id("1000G", "AFR"), gotG)
            self.assertNotIn(db_group_id("HGDP"), gotG)
            self.assertIn(db_subpop_group_id("HGDP", "EAS"), gotT)
            self.assertNotIn(db_group_id("1000G"), gotT)
        finally:
            reader.close()


# --------------------------------------------------------------------------- #
# (4) Two databases: per-db + global (deduped) correct; global = union of panels.
# --------------------------------------------------------------------------- #
class TestTwoDatabases(_CompileFixture):
    def test_per_db_and_global_union(self):
        # Carriers in BOTH databases at the same pos, same alt.
        # 1000G: S_EUR_F1 het, S_AFR_M1 hemi? autosome -> diploid het.
        # HGDP:  H_EAS_F1 hom.
        carrier_gts = {
            "S_EUR_F1": "0|1",
            "S_AFR_M1": "0|1",
            "H_EAS_F1": "1|1",
        }
        entry = ("S_EUR_F1:0|1,S_AFR_M1:0|1,H_EAS_F1:1|1;C,A;rsAB;0.2")
        mapping = {"%s,777" % self.CHROM: entry}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            got = self._assert_all_groups(reader, 777, "A", carrier_gts)
            # 1000G db: panel of 7 -> AN 14. AC = het(1) + het(1) = 2.
            d1 = got[db_group_id("1000G")]
            self.assertEqual((d1.AC, d1.AN, d1.n_called_indiv, d1.n_carrier_indiv),
                             (2, 14, 7, 2))
            # HGDP db: panel of 3 -> AN 6. AC = hom(2). carrier 1, hom 1.
            d2 = got[db_group_id("HGDP")]
            self.assertEqual((d2.AC, d2.AN, d2.n_called_indiv, d2.n_carrier_indiv,
                              d2.n_hom_indiv), (2, 6, 3, 1, 1))
            # global = union of BOTH panels (10 individuals, AN 20). AC = 2 + 2 = 4.
            g = got[GLOBAL_GROUP_ID]
            self.assertEqual((g.AC, g.AN, g.n_called_indiv, g.n_carrier_indiv,
                              g.n_hom_indiv), (4, 20, 10, 3, 1))
            self.assertAlmostEqual(g.allele_freq(), 4.0 / 20.0)
            # Cross-check: db AN's sum to global AN (disjoint panels).
            self.assertEqual(d1.AN + d2.AN, g.AN)
            self.assertEqual(d1.n_called_indiv + d2.n_called_indiv, g.n_called_indiv)
        finally:
            reader.close()

    def test_derived_max_subpop_full_panel(self):
        # max_subpop_af must use full-panel AN denominators per subpop.
        # HGDP::EAS: H_EAS_F1 hom (AC2) over a 3-individual EAS panel (AN 6) -> 2/6.
        # 1000G::EUR: S_EUR_F1 het (AC1) over 4-individual EUR panel (AN 8) -> 1/8.
        # 1000G::AFR: S_AFR_M1 het (AC1) over 3-individual AFR panel (AN 6) -> 1/6.
        # max = EAS 2/6 = 1/3.
        carrier_gts = {
            "S_EUR_F1": "0|1",
            "S_AFR_M1": "0|1",
            "H_EAS_F1": "1|1",
        }
        entry = ("S_EUR_F1:0|1,S_AFR_M1:0|1,H_EAS_F1:1|1;C,A;rsAB;0.2")
        mapping = {"%s,777" % self.CHROM: entry}
        binp, idxp, stats = self._compile(mapping)
        reader = RegistryReader(binp, idxp)
        try:
            der = reader.derived(777, "A")
            self.assertTrue(der["observed"])
            self.assertAlmostEqual(der["max_subpop_af"], 2.0 / 6.0)
            self.assertEqual(der["max_subpop_af_label"],
                             db_subpop_group_id("HGDP", "EAS"))
            self.assertAlmostEqual(der["global_af"], 4.0 / 20.0)
        finally:
            reader.close()


# --------------------------------------------------------------------------- #
# Indel skipping + empty-carrier skipping + absent-record semantics.
# --------------------------------------------------------------------------- #
class TestSkipAndAbsent(_CompileFixture):
    def test_indel_skipped_not_crash(self):
        # A multi-char alt (insertion) and a multi-char ref (deletion) must be
        # SKIPPED, counted, and not crash. A valid SNP alongside is still written.
        entry_snp = "S_EUR_F1:0|1;A,G;rsS;0.05"
        entry_ins = "S_AFR_F1:0|1;A,GT;rsI;0.01"   # multi-char alt -> indel
        entry_del = "S_AFR_F2:0|1;AT,A;rsD;0.01"   # multi-char ref -> indel
        mapping = {"%s,100" % self.CHROM: entry_snp,
                   "%s,200" % self.CHROM: entry_ins + "$" + entry_del}
        binp, idxp, stats = self._compile(mapping)
        self.assertEqual(stats["n_written"], 1)
        self.assertEqual(stats["n_skipped_indel"], 2)
        reader = RegistryReader(binp, idxp)
        try:
            self.assertIsNotNone(reader.lookup(100, "G"))
            # indel positions produced no SNP records
            self.assertIsNone(reader.lookup(200, "G"))
            self.assertIsNone(reader.lookup(200, "T"))
        finally:
            reader.close()

    def test_star_spanning_deletion_alt_skipped(self):
        # A '*' (spanning/overlapping-deletion placeholder) is single-char but NOT
        # a real substitution -> must be SKIPPED, not written as a SNP record.
        entry_snp = "S_EUR_F1:0|1;A,G;rsS;0.05"
        entry_star = "S_AFR_F1:0|1;A,*;rsStar;0.01"
        mapping = {"%s,100" % self.CHROM: entry_snp,
                   "%s,150" % self.CHROM: entry_star}
        binp, idxp, stats = self._compile(mapping)
        self.assertEqual(stats["n_written"], 1)
        self.assertGreaterEqual(stats["n_skipped_indel"], 1)  # '*' counted as non-SNP
        reader = RegistryReader(binp, idxp)
        try:
            self.assertIsNotNone(reader.lookup(100, "G"))
            self.assertIsNone(reader.lookup(150, "*"))
        finally:
            reader.close()

    def test_empty_carriers_skipped(self):
        # An entry with NO listed carriers (empty samples field) carries no signal.
        entry = ";A,G;rsE;0.0"
        mapping = {"%s,100" % self.CHROM: entry}
        binp, idxp, stats = self._compile(mapping)
        self.assertEqual(stats["n_written"], 0)
        self.assertEqual(stats["n_skipped_empty"], 1)
        reader = RegistryReader(binp, idxp)
        try:
            self.assertEqual(len(reader), 0)
            self.assertIsNone(reader.lookup(100, "G"))
        finally:
            reader.close()


# --------------------------------------------------------------------------- #
# samplesID reader sanity (columns + subpop_field).
# --------------------------------------------------------------------------- #
class TestSamplesIDReader(unittest.TestCase):
    def test_superpopulation_default(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "s.txt")
        write_samplesid(p, G1000_ROWS)
        meta = tc.read_samplesid(p, "1000G")
        self.assertEqual(meta["S_EUR_M1"], ("1000G", "EUR", "male"))
        self.assertEqual(meta["S_AFR_F1"], ("1000G", "AFR", "female"))
        self.assertEqual(len(meta), len(G1000_ROWS))

    def test_population_field_option(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "s.txt")
        write_samplesid(p, G1000_ROWS)
        meta = tc.read_samplesid(p, "1000G", subpop_field="population")
        self.assertEqual(meta["S_EUR_M1"], ("1000G", "GBR", "male"))
        self.assertEqual(meta["S_AFR_F1"], ("1000G", "YRI", "female"))

    def test_build_sample_meta_union(self):
        d = tempfile.mkdtemp()
        g = os.path.join(d, "g.txt")
        h = os.path.join(d, "h.txt")
        write_samplesid(g, G1000_ROWS)
        write_samplesid(h, HGDP_ROWS)
        meta, overlaps = tc.build_sample_meta({"1000G": g, "HGDP": h})
        self.assertEqual(len(meta), len(G1000_ROWS) + len(HGDP_ROWS))
        self.assertEqual(overlaps, [])
        self.assertEqual(meta["H_EAS_M1"], ("HGDP", "EAS", "male"))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Every allele of a variant position, read the same way by every reader (CRISPRme 2.9.0).

CRISPRitz 2.10.0 codes every ALT of a split multi-allelic site (up to 2.9.1 only the last
record's), lists each in the per-sample dictionary with its own record's carriers, and
writes Q at a position whose reference base and SNP ALTs hold all four bases (N stays a
gap). The readers of CRISPRme had to change with it:

  * site_alleles: one position's genotypes, each read against its own ALT's allele index
    (a record listing several ALTs on one VCF line keeps its raw 1|2 under each ALT), and
    two records of one REF>ALT merged into one allele;
  * the Tier-0 / Tier-1 compiles and the dictionary decode read through it: a 1|2
    carrier's second ALT counted as reference, and two records of one REF>ALT stopped
    the Tier-1 compile of the whole chromosome;
  * observed_haplotypes: a phased slot holding two ALTs of one position is no haplotype
    (PUTATIVE), and no combination takes two ALTs of one position; above
    CRISPRME_IUPAC_CAP variants a sample gets at most 16 unions, past that one, holding at
    each column the ALT nearest the guide;
  * population_summary decides a carrier with two ALTs on one slot, or an unphased
    genotype, without phase;
  * the SNP+SNP companion compares positions, not variant tokens;
  * the dictionary path keeps only the ALTs the target's own column holds, as the
    observed path does (#139), and a Q column admits its position's ALTs only where they
    and the reference hold all four bases, and of them only the bases (a '*' is none);
  * Q in every reader of variant codes (pam_filter, the SNP overlay of the fake-indel
    genome, the indel path's resolver, the site rule of crispritz_bins, twopass_find,
    the cluster sort);
  * the merge of a cluster's rows keeps each rsID and AF in the place of its variant;
  * the sites-only (mega) registry reads a line listing several ALTs one ALT at a time.

The SNP post-analysis is the real new_simple_analysis.py, run as the pipeline runs it,
on tiers and dictionaries the repository's own compilers build from dictionary values
written as CRISPRitz writes them; the merge is the real remove_contiguous_samples.py.

Run with:
    cd PostProcess && python3 -m unittest test_site_alleles -v
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_dictless_tiers  # noqa: E402
import build_mega_indel_af  # noqa: E402
import build_mega_registry  # noqa: E402
import crispritz_bins  # noqa: E402
import indel_snp_cis as isc  # noqa: E402
import overlay_indel_snps  # noqa: E402
import pam_filter  # noqa: E402
import population_summary  # noqa: E402
import simple_analysis_registry as sar  # noqa: E402
import site_alleles  # noqa: E402
import tier0_compile  # noqa: E402
import tier0_registry  # noqa: E402
import tier1_genotypes  # noqa: E402
import twopass_find  # noqa: E402

SPACER = "CTAACAGTTGCTTTTATCAC"
CRRNA = SPACER + "NNN"
SITE = SPACER + "TGG"
FLANK = "ACGTGCATTC" * 10
AT = len(FLANK)  # 0-based start of the site on chrT
SAMPLES = ("O", "P", "S", "U", "V", "X", "Y", "Z")


def _pos1(col):
    """1-based chrT position of the site's column ``col``."""
    return AT + col + 1


# --------------------------------------------------------------------------- #
# site_alleles
# --------------------------------------------------------------------------- #
class OnePositionsAlleles(unittest.TestCase):
    def test_split_and_biallelic_positions_come_back_as_they_are(self):
        for entries in ([("A", {"S": "1|0"})],
                        [("A", {"S": "1|0"}), ("T", {"S": "0|1", "X": "1|1"})],
                        [("A", {"S": "1/0", "H": "1", "M": ".|1"})], [],
                        [("A", {"S": "1|0", "E": ""})]):  # a carrier with no genotype
            self.assertIs(site_alleles.per_alt_entries(entries), entries)

    def test_a_line_of_several_alts_is_read_per_alt(self):
        got = dict(site_alleles.per_alt_entries([
            ("C", {"S": "1|2", "X": "1|0"}), ("G", {"S": "1|2", "Y": "2/2"})]))
        self.assertEqual(got["C"], {"S": "1|0", "X": "1|0"})
        self.assertEqual(got["G"], {"S": "0|1", "Y": "1/1"})

    def test_an_indel_alt_before_the_snp_one(self):
        # A>AT,C: C is allele 2 and its own entry; the indel has none in the dict
        got = dict(site_alleles.per_alt_entries([("C", {"S": "0|2", "X": "2|1"})]))
        self.assertEqual(got["C"], {"S": "0|1", "X": "1|0"})

    def test_genotypes_that_cannot_tell_lose_their_phase(self):
        got = dict(site_alleles.per_alt_entries([("C", {"S": "1|2"}), ("G", {"S": "2|1"})]))
        self.assertEqual(got, {"C": {"S": "1/0"}, "G": {"S": "1/0"}})

    def test_merge_gt(self):
        cases = [("1|0", "0|1", "1|1"), ("1|0", "1|0", "1|0"), ("0|0", "0|1", "0|1"),
                 ("1|.", ".|0", "1|0"), ("1", "1", "1"), ("1|0", "0/1", "1/0"),
                 ("1/1", "0|1", "1/1"), ("1", "0|1", "1/0")]
        for a, b, want in cases:
            with self.subTest(a=a, b=b):
                self.assertEqual(site_alleles.merge_gt(a, b), want)

    def test_two_records_of_one_ref_alt_are_one_allele(self):
        recs = [("C", "T", "rs1", "0.1", {"S": "1|0", "X": "1|1"}),
                ("C", "G", ".", "0.2", {"Y": "0|1"}),
                ("C", "T", "rs9", "0.3", {"S": "0|1", "Z": "0|1"})]
        got = site_alleles.read_records(recs)
        self.assertEqual([(r[0], r[1]) for r in got], [("C", "T"), ("C", "G")])
        ref, alt, rsid, af, carriers = got[0]
        self.assertEqual((rsid, af), ("rs9", "0.3"))  # the more frequent record's
        self.assertEqual(carriers, {"S": "1|1", "X": "1|1", "Z": "0|1"})
        two = recs[:2]
        self.assertIs(site_alleles.read_records(two), two)

    def test_two_alts_on_a_slot(self):
        site = {"G": {"S": "0|1", "O": "1|0"}, "T": {"S": "0|1", "O": "0|1", "U": "1/0"}}
        self.assertEqual(site_alleles.two_alts_on_a_slot(site), {"S": {1}})


# --------------------------------------------------------------------------- #
# The tiers and the dictionary decode
# --------------------------------------------------------------------------- #
def _write_samplesid(path):
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for s in SAMPLES:
            fh.write(f"{s}\tPOP\tSUP\tfemale\n")


def _dict_value(entries):
    """A dictionary value as CRISPRitz writes it: ``(ref, alt, rsid, af, {sample: gt})``
    entries joined by '$', carriers sorted."""
    out = []
    for ref, alt, rsid, af, carriers in entries:
        samples = ",".join(sorted(f"{s}:{g}" for s, g in carriers.items()))
        out.append(f"{samples};{ref},{alt};{rsid};{af}")
    return "$".join(out)


def _build_tiers(root, values, keep_dict):
    """Dictionaries/ of dataset ``ds`` under ``root``: my_dict_chrT.json holding
    ``values`` ({pos1: [entries]}), its Tier-0 registry and Tier-1 store compiled by
    build_dictless_tiers, the dict removed unless ``keep_dict``. Returns the dict path."""
    ddir = os.path.join(root, "Dictionaries", "dictionaries_ds")
    os.makedirs(ddir)
    dict_path = os.path.join(ddir, "my_dict_chrT.json")
    with open(dict_path, "w") as fh:
        json.dump({f"chrT,{p}": _dict_value(e) for p, e in values.items()}, fh)
    sid = os.path.join(root, "ds.samplesID.txt")
    _write_samplesid(sid)
    build_dictless_tiers.emit_dictless_tiers(dict_path, {"ds": sid}, "chrT")
    if not keep_dict:
        os.remove(dict_path)
    return dict_path


def _store_of_2_8_1(root, values):
    """Rewrite the Tier-1 store ``_build_tiers`` made under ``root`` as CRISPRme 2.8.1
    compiled it: each dictionary entry's carriers with their genotypes as written (a line
    of several ALTs keeps its raw 1|2 under each of its ALTs), not read per ALT."""
    d = os.path.join(root, "Dictionaries", "genotypes_ds")
    paths = (os.path.join(d, "gt_chrT.bin"), os.path.join(d, "gt_chrT.idx"))
    reader = tier1_genotypes.GenotypeReader(*paths)
    axis = reader._axis
    reader.close()
    tier1_genotypes.compile_genotypes(
        [(p, alt, dict(carriers), ref) for p, entries in values.items()
         for ref, alt, _rsid, _af, carriers in entries if carriers], axis, *paths)


class TheTiersOfALineOfSeveralAlts(unittest.TestCase):
    """A VCF line C>A,T: S is 1|2, X 1|0, Y 2|2. The registry counted every allele of a
    carrier against "1", so T's carriers, written 2, were reference (registry AC 2
    against a raw AC 2977 at chr22:10570455 in the 2019 HGDP data)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="sa_tiers_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.pos = _pos1(5)
        _build_tiers(self.tmp, {self.pos: [
            ("C", "A", "rs1", "0.25", {"S": "1|2", "X": "1|0"}),
            ("C", "T", "rs1", "0.25", {"S": "1|2", "Y": "2|2"})]}, keep_dict=False)
        d = os.path.join(self.tmp, "Dictionaries")
        self.reg = tier0_registry.RegistryReader(
            os.path.join(d, "registry_ds", "reg_chrT.bin"),
            os.path.join(d, "registry_ds", "reg_chrT.idx"))
        self.gt = tier1_genotypes.GenotypeReader(
            os.path.join(d, "genotypes_ds", "gt_chrT.bin"),
            os.path.join(d, "genotypes_ds", "gt_chrT.idx"))
        self.addCleanup(self.reg.close)
        self.addCleanup(self.gt.close)

    def test_each_alt_counts_its_own_copies(self):
        g = tier0_registry.GLOBAL_GROUP_ID
        a = self.reg.lookup(self.pos, "A")[g]
        t = self.reg.lookup(self.pos, "T")[g]
        self.assertEqual((a.AC, a.AN, a.n_carrier_indiv, a.n_hom_indiv), (2, 16, 2, 0))
        self.assertEqual((t.AC, t.AN, t.n_carrier_indiv, t.n_hom_indiv), (3, 16, 2, 1))

    def test_the_store_holds_each_alt_s_slots(self):
        self.assertEqual(dict(self.gt.carriers(self.pos, "A")), {"S": "1|0", "X": "1|0"})
        self.assertEqual(dict(self.gt.carriers(self.pos, "T")), {"S": "0|1", "Y": "1|1"})


class TwoRecordsOfOneRefAlt(unittest.TestCase):
    """Decision B2: the carriers of two records of one REF>ALT are merged, a slot holding
    the ALT if either record puts it there. Two entries of one (pos, alt) stopped the
    Tier-1 compile of the chromosome ("duplicate (pos, alt) record"), and the registry
    kept a group per record."""

    def test_the_tiers_hold_one_allele_with_both_records_carriers(self):
        tmp = tempfile.mkdtemp(prefix="sa_dup_")
        self.addCleanup(shutil.rmtree, tmp, True)
        pos = _pos1(5)
        _build_tiers(tmp, {pos: [
            ("C", "T", "rs1", "0.1", {"S": "1|0", "X": "1|1"}),
            ("C", "G", "rs2", "0.1", {"O": "0|1"}),
            ("C", "T", "rs1", "0.1", {"S": "0|1", "Z": "0|1"})]}, keep_dict=False)
        d = os.path.join(tmp, "Dictionaries")
        self.assertTrue(os.path.exists(os.path.join(d, "genotypes_ds", "gt_chrT.idx")),
                        "the Tier-1 store was not written")
        gt = tier1_genotypes.GenotypeReader(os.path.join(d, "genotypes_ds", "gt_chrT.bin"),
                                            os.path.join(d, "genotypes_ds", "gt_chrT.idx"))
        self.addCleanup(gt.close)
        self.assertEqual(gt.alts_at(pos), ["G", "T"])
        self.assertEqual(dict(gt.carriers(pos, "T")), {"S": "1|1", "X": "1|1", "Z": "0|1"})
        reg = tier0_registry.RegistryReader(os.path.join(d, "registry_ds", "reg_chrT.bin"),
                                            os.path.join(d, "registry_ds", "reg_chrT.idx"))
        self.addCleanup(reg.close)
        t = reg.lookup(pos, "T")[tier0_registry.GLOBAL_GROUP_ID]
        self.assertEqual((t.AC, t.n_carrier_indiv), (5, 3))


class TheRegistryListsTheFourBasePositions(unittest.TestCase):
    """CRISPRitz 2.10.0 codes every PASS SNP ALT, one no sample carries too, and writes Q
    where they and the reference base hold all four bases. The registry holds no record
    of an ALT nobody carries (it has no group), so ``alts_at`` gave such a position three
    bases, and a reader of its Q column could not tell it from a neighbour's."""

    def test_an_alt_no_sample_carries_counts(self):
        tmp = tempfile.mkdtemp(prefix="sa_four_")
        self.addCleanup(shutil.rmtree, tmp, True)
        p5, p9 = _pos1(5), _pos1(9)
        _build_tiers(tmp, {
            p5: [("C", "A", "rsA", "0.2", {"S": "1|0"}),
                 ("C", "G", "rsG", "0", {}),
                 ("C", "T", "rsT", "0.1", {"V": "1|1"})],
            p9: [("A", "C", "rsC", "0.1", {"S": "1|0"}),
                 ("A", "G", "rsG2", "0", {})]}, keep_dict=False)
        d = os.path.join(tmp, "Dictionaries", "registry_ds")
        reg = tier0_registry.RegistryReader(os.path.join(d, "reg_chrT.bin"),
                                            os.path.join(d, "reg_chrT.idx"))
        self.addCleanup(reg.close)
        self.assertEqual(reg.alts_at(p5), ["A", "T"])  # no record of G, as before
        self.assertTrue(reg.four_base(p5))
        self.assertFalse(reg.four_base(p9))  # A, C and G
        self.assertFalse(reg.four_base(p5 + 1))
        self.assertEqual(reg.manifest["four_base_positions"], [p5])


class TheDictionaryDecode(unittest.TestCase):
    def test_a_line_of_several_alts_and_a_repeated_allele(self):
        entry = _dict_value([("C", "A", "rs1", "0.25", {"S": "1|2", "X": "1|0"}),
                             ("C", "T", "rs1,rs7", "0.25", {"S": "1|2", "Y": "2|2"}),
                             ("C", "A", "rs1", "0.3", {"Z": "0|1"})])
        snps, samples, rsids, afs, infos = sar.legacy_5tuple_from_entry(entry, "chrT", 104)
        self.assertEqual(snps, ["A", "T"])
        self.assertEqual(samples, [["S:1|0", "X:1|0", "Z:0|1"], ["S:0|1", "Y:1|1"]])
        self.assertEqual(rsids, ["rs1", "rs1;rs7"])
        self.assertEqual(afs, ["0.3", "0.25"])
        self.assertEqual(infos, ["chrT_105_C_A", "chrT_105_C_T"])

    def test_split_records_decode_as_before(self):
        entry = _dict_value([("C", "A", "rs1", "0.25", {"S": "1|0"}),
                             ("C", "T", "rs2", "0.5", {"S": "0|1", "Y": "1|1"})])
        self.assertEqual(sar.legacy_5tuple_from_entry(entry, "chrT", 104),
                         (["A", "T"], [["S:1|0"], ["S:0|1", "Y:1|1"]], ["rs1", "rs2"],
                          ["0.25", "0.5"], ["chrT_105_C_A", "chrT_105_C_T"]))

    def test_the_indel_path_reads_the_dictionary_the_same_way(self):
        value = _dict_value([("C", "A", "rs1", "0.25", {"S": "1|2", "X": "1|0"}),
                             ("C", "T", "rs1", "0.25", {"S": "1|2", "Y": "2|2"}),
                             ("C", "A", "rs1", "0.3", {"Z": "0|1"})])
        got = isc.parse_dict_value(value)
        self.assertEqual([(r[1], r[4]) for r in got],
                         [("A", {"S": "1|0", "X": "1|0", "Z": "0|1"}),
                          ("T", {"S": "0|1", "Y": "1|1"})])


class ThePopulationSummaryOfTwoAltsOnASlot(unittest.TestCase):
    """A two-variant row G at column 5 + C at column 9: S holds G and C on slot 1, and T
    of column 5 on slot 1 too (split records phased one by one). That slot is no
    haplotype: S is a carrier of the row, decided without phase, with no haplotype
    count. It was counted as one cis haplotype."""

    def test_the_slot_is_not_counted_as_a_haplotype(self):
        tmp = tempfile.mkdtemp(prefix="sa_pop_")
        self.addCleanup(shutil.rmtree, tmp, True)
        p5, p9 = _pos1(5), _pos1(9)
        _build_tiers(tmp, {
            p5: [("C", "G", "rs1", "0.1", {"S": "0|1", "O": "0|1"}),
                 ("C", "T", "rs2", "0.1", {"S": "0|1"})],
            p9: [("T", "C", "rs3", "0.1", {"S": "0|1", "O": "0|1"})]}, keep_dict=False)
        d = os.path.join(tmp, "Dictionaries")
        reg = tier0_registry.RegistryReader(os.path.join(d, "registry_ds", "reg_chrT.bin"),
                                            os.path.join(d, "registry_ds", "reg_chrT.idx"))
        gt = tier1_genotypes.GenotypeReader(os.path.join(d, "genotypes_ds", "gt_chrT.bin"),
                                            os.path.join(d, "genotypes_ds", "gt_chrT.idx"))
        self.addCleanup(reg.close)
        self.addCleanup(gt.close)
        summary = population_summary.summarize(
            [(p5, "G"), (p9, "C")], reg, gt, gt._axis,
            tier0_compile.ploidy_of_for_chrom("chrT"), phased=True)
        g = summary.groups[tier0_registry.GLOBAL_GROUP_ID]
        self.assertEqual((g.n_carrier, g.allele_count), (2, 1))

    def test_a_store_of_an_earlier_release(self):
        """A line C>A,T at column 5 (S 1|2: T on slot 1) and A>G at column 9 (S 0|1), in a
        Tier-1 store compiled by 2.8.1, which keeps the raw 1|2 under A and T: the row T + G
        is S's haplotype 1, one cis copy, as in the store this release compiles."""
        p5, p9 = _pos1(5), _pos1(9)
        values = {p5: [("C", "A", "rs1", "0.1", {"S": "1|2", "X": "1|0"}),
                       ("C", "T", "rs1", "0.1", {"S": "1|2", "Y": "2|2"})],
                  p9: [("A", "G", "rs9", "0.1", {"S": "0|1"})]}
        for old_store in (False, True):
            with self.subTest(store_of_2_8_1=old_store):
                tmp = tempfile.mkdtemp(prefix="sa_pop281_")
                self.addCleanup(shutil.rmtree, tmp, True)
                _build_tiers(tmp, values, keep_dict=False)
                if old_store:
                    _store_of_2_8_1(tmp, values)
                d = os.path.join(tmp, "Dictionaries")
                reg = tier0_registry.RegistryReader(
                    os.path.join(d, "registry_ds", "reg_chrT.bin"),
                    os.path.join(d, "registry_ds", "reg_chrT.idx"))
                gt = tier1_genotypes.GenotypeReader(
                    os.path.join(d, "genotypes_ds", "gt_chrT.bin"),
                    os.path.join(d, "genotypes_ds", "gt_chrT.idx"))
                self.addCleanup(reg.close)
                self.addCleanup(gt.close)
                summary = population_summary.summarize(
                    [(p5, "T"), (p9, "G")], reg, gt, gt._axis,
                    tier0_compile.ploidy_of_for_chrom("chrT"), phased=True)
                g = summary.groups[tier0_registry.GLOBAL_GROUP_ID]
                self.assertEqual((g.n_carrier, g.allele_count), (1, 1))

    def test_a_carrier_of_an_unphased_genotype(self):
        """A carrier whose genotype at a variant of the row is unphased ('/') is decided
        without phase: its slot order is no phase. A phased store holds such a genotype
        for a line C>A,T whose only carrier S is 1|2 (which of S's slots holds A, the
        genotypes cannot tell: 1/0 under A and under T) and for the merge of a phased and
        an unphased record of one REF>ALT. With G at column 9 on S's slot 0, A + G and
        T + G were each one cis copy, though A and T cannot both be on slot 0; and S was
        no carrier of the merged G + C, slot 0 against slot 1."""
        p5, p9 = _pos1(5), _pos1(9)
        cases = (
            ({p5: [("C", "A", "rs1", "0.1", {"S": "1|2"}),
                   ("C", "T", "rs1", "0.1", {"S": "1|2"})],
              p9: [("A", "G", "rs9", "0.1", {"S": "1|0"})]},
             ([(p5, "A"), (p9, "G")], [(p5, "T"), (p9, "G")])),
            ({p5: [("C", "G", "rs1", "0.1", {"S": "1|0"}),
                   ("C", "G", "rs2", "0.1", {"S": "0/1"})],
              p9: [("A", "C", "rs9", "0.1", {"S": "0|1"})]},
             ([(p5, "G"), (p9, "C")],)),
        )
        for values, combinations in cases:
            tmp = tempfile.mkdtemp(prefix="sa_popun_")
            self.addCleanup(shutil.rmtree, tmp, True)
            _build_tiers(tmp, values, keep_dict=False)
            d = os.path.join(tmp, "Dictionaries")
            reg = tier0_registry.RegistryReader(
                os.path.join(d, "registry_ds", "reg_chrT.bin"),
                os.path.join(d, "registry_ds", "reg_chrT.idx"))
            gt = tier1_genotypes.GenotypeReader(
                os.path.join(d, "genotypes_ds", "gt_chrT.bin"),
                os.path.join(d, "genotypes_ds", "gt_chrT.idx"))
            self.addCleanup(reg.close)
            self.addCleanup(gt.close)
            for pos_alts in combinations:
                with self.subTest(pos_alts=pos_alts):
                    summary = population_summary.summarize(
                        pos_alts, reg, gt, gt._axis,
                        tier0_compile.ploidy_of_for_chrom("chrT"), phased=True)
                    g = summary.groups.get(tier0_registry.GLOBAL_GROUP_ID)
                    # a carrier, with no haplotype count
                    self.assertEqual((g.n_carrier, g.allele_count) if g else None, (1, 0))


# --------------------------------------------------------------------------- #
# Q in the readers of variant codes
# --------------------------------------------------------------------------- #
class TheFourAlleleCode(unittest.TestCase):
    def test_pam_filter_reads_q_as_every_base_and_n_as_none(self):
        self.assertTrue(pam_filter._pam_can_match("QGG", "NGG"))
        self.assertTrue(pam_filter._pam_can_match("AQG", "NRG"))
        self.assertFalse(pam_filter._pam_can_match("NGG", "NGG"))

    def test_the_overlay_copies_q_onto_the_fake_indel_flanks(self):
        enriched = "ACGTQACGTACGTACGTACGTACGTACGTACGT"
        sub = "ACGTAACGTACGTACGTACGTACGTAACGTACG"
        got = overlay_indel_snps.overlay_sub_fasta(sub, 0, "T", "TA", enriched)
        self.assertEqual(got[4], "Q")

    def test_the_indel_resolver_reads_a_q_column(self):
        self.assertIn(5, isc.ambiguous_columns(SPACER[:5] + "Q" + SPACER[6:] + "TGG"))
        col = isc.SnpColumn(5, 104, "C", {"A": {"S": "1|0"}, "G": {"U": "0|1"},
                                          "T": {"V": "1|1"}}, "Q", "+")
        self.assertEqual(col.code, frozenset("ACGT"))

    def test_crispritz_bins_reads_a_q_contig_and_a_q_pam_base(self):
        tmp = tempfile.mkdtemp(prefix="sa_bins_")
        self.addCleanup(shutil.rmtree, tmp, True)
        fa = os.path.join(tmp, "chrT.fa")
        seq = "ACGT" * 10 + SPACER + "TQG" + "ACGT" * 10
        with open(fa, "w") as fh:
            fh.write(">chrT\n" + seq + "\n")
        self.assertEqual(crispritz_bins._read_contig(fa), ("chrT", seq))
        self.assertTrue(crispritz_bins.holds_target(seq, "NGG", 23, False, 0))

    def test_twopass_find_reads_q_as_every_base(self):
        sets = twopass_find.iupac_window_to_sets("AQN")
        self.assertEqual(sets[1], frozenset("ACGT"))

    def test_the_brute_force_truth_reads_q_as_every_base(self):
        """test/benchmark/generate_brute_force.py, the benchmarks' truth (it imports
        Biopython, which this job does not install: its table is read from the file)."""
        import ast
        path = os.path.join(os.path.dirname(HERE), "test", "benchmark", "generate_brute_force.py")
        with open(path) as fh:
            tree = ast.parse(fh.read())
        table = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                     and getattr(n.targets[0], "id", "") == "iupac_dict")
        self.assertEqual(table.get("Q"), {"A", "C", "G", "T"})

    def test_the_cluster_sort_takes_q(self):
        with open(os.path.join(HERE, "cluster.dict.py")) as fh:
            line = next(l for l in fh if l.startswith("alphabet = "))
        alphabet = eval(line.split("=", 1)[1].split("#")[0])
        self.assertIn("Q", alphabet)
        self.assertLess(alphabet.index("Q"), alphabet.index("A"))


# --------------------------------------------------------------------------- #
# The SNP post-analysis: the real new_simple_analysis.py
# --------------------------------------------------------------------------- #
# 0-based columns of a row new_simple_analysis.py writes (bestCFD.txt)
DNA, SAMPLES_COL, RSID, AF, SNP = 2, 12, 15, 16, 17


def _site(**bases):
    """The site's reference: SITE with ``c<col>=<base>`` replaced."""
    s = list(SITE)
    for k, b in bases.items():
        s[int(k[1:])] = b
    return "".join(s)


def _window(ref, codes):
    """The DNA column of the variant search's hit: ``ref`` with each {col: code}."""
    w = list(ref)
    for col, code in codes.items():
        w[col] = code
    return "".join(w)


class _SnpPostAnalysis(unittest.TestCase):
    """Runs the real new_simple_analysis.py as post_analisi_snp.sh runs it: the
    reference chromosome, the dataset's dictionary path (its tiers beside it), the
    variant search's targets, the PAM, the output stem and the mismatch budget."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="sa_nsa_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def analyse(self, ref, windows, values, tiers=True, keep_dict=False, mode="0",
                mms=4, sites_only_vcf=None, genotypes=True, store_of_2_8_1=False):
        fasta = os.path.join(self.tmp, "chrT.fa")
        with open(fasta, "w") as fh:
            fh.write(">chrT\n" + FLANK + ref + FLANK + "\n")
        if sites_only_vcf is not None:
            ddir = os.path.join(self.tmp, "Dictionaries")
            build_mega_registry.build(sites_only_vcf, "chrT",
                                      os.path.join(ddir, "registry_ds"))
            dict_path = os.path.join(ddir, "dictionaries_ds", "my_dict_chrT.json")
        elif tiers:
            dict_path = _build_tiers(self.tmp, values, keep_dict)
            if store_of_2_8_1:
                _store_of_2_8_1(self.tmp, values)
            if not genotypes:  # the registry alone, as `download --no-genotypes` leaves it
                shutil.rmtree(os.path.join(self.tmp, "Dictionaries", "genotypes_ds"))
        else:
            ddir = os.path.join(self.tmp, "Dictionaries", "dictionaries_ds")
            os.makedirs(ddir)
            dict_path = os.path.join(ddir, "my_dict_chrT.json")
            with open(dict_path, "w") as fh:
                json.dump({f"chrT,{p}": _dict_value(e) for p, e in values.items()}, fh)
        pam = os.path.join(self.tmp, "pam.txt")
        with open(pam, "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")
        targets = os.path.join(self.tmp, "chrT.targets.txt")
        with open(targets, "w") as fh:
            fh.write("#header\n")
            for w in windows:
                fh.write("\t".join(["X", CRRNA, w, "chrT", str(AT), str(AT), "+", "0", "0",
                                    "0", "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"])
                         + "\n")
        out = os.path.join(self.tmp, "out")
        env = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        env.update(CRISPRME_FAST_MODE=mode, CRISPRME_SKIP_SCORER_ENV="1")
        done = subprocess.run(
            [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), fasta, dict_path,
             targets, pam, out, str(mms)],
            cwd=HERE, env=env, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(done.stderr, "", "a pipeline stage must write nothing on stderr")
        with open(out + ".bestCFD.txt") as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh][1:]
        companions = {}
        for name in ("phase_confirmation", "snp_snp_cooc"):
            path = "%s.%s.tsv" % (out, name)
            if os.path.exists(path):
                with open(path) as fh:
                    lines = [ln.rstrip("\n") for ln in fh]
                # the column names are the last '#' line
                head = [ln for ln in lines if ln.startswith("#")][-1].lstrip("#").split("\t")
                companions[name] = [dict(zip(head, ln.split("\t")))
                                    for ln in lines if not ln.startswith("#")]
        return rows, companions

    def variant_rows(self, rows):
        """{tuple of variant tokens: row} of the rows that name variants."""
        out = {}
        for r in rows:
            if r[SNP] not in ("n", "NA"):
                out[tuple(r[SNP].split(","))] = r
        return out

    def assert_no_row_names_two_alts_of_a_position(self, rows):
        for r in rows:
            if r[SNP] in ("n", "NA"):
                continue
            positions = [t.rsplit("_", 3)[1] for t in r[SNP].split(",")]
            self.assertEqual(len(positions), len(set(positions)), r[SNP])


class TwoAltsOnOneSlotThroughThePostAnalysis(_SnpPostAnalysis):
    """MULTIALLELIC section 6.6's case, made small: column 5 of the site is a reference C
    that the guide reads A; split records C>A and C>T, phased one by one, put both on
    slot 1 of S (0|1 and 0|1), as NA18869's at chr22:10,661,629. The row naming both
    alleles was reported CONFIRMED, for S alone, and the SNP+SNP companion took it for
    a co-occurrence."""

    REF = _site(c5="C")
    VALUES = {_pos1(5): [("C", "A", "rsA", "0.2", {"O": "1|0", "S": "0|1"}),
                         ("C", "T", "rsT", "0.1", {"S": "0|1", "P": "0|1"})]}

    def test_per_sample_on_the_tiers(self):
        rows, comp = self.analyse(self.REF, [_window(self.REF, {5: "H"})], self.VALUES)
        self.assert_no_row_names_two_alts_of_a_position(rows)
        got = self.variant_rows(rows)
        a, t = ("chrT_%d_C_A" % _pos1(5),), ("chrT_%d_C_T" % _pos1(5),)
        self.assertEqual(set(got), {a, t})
        self.assertEqual(set(got[a][SAMPLES_COL].split(",")), {"O", "S"})
        self.assertEqual(set(got[t][SAMPLES_COL].split(",")), {"P", "S"})
        phases = {row["SNP"]: row["Phase_Confirmation"] for row in comp["phase_confirmation"]}
        self.assertEqual(phases, {a[0]: "PUTATIVE", t[0]: "PUTATIVE"})
        self.assertNotIn("snp_snp_cooc", comp)


class ADenseWindowOfTwoAltsAColumn(_SnpPostAnalysis):
    """An unphased sample S carries two ALTs (1/2, split into 1/0 and 0/1) at each of 12
    columns of the site's protospacer and at the last base of its PAM; one of the two is
    the guide's base (the PAM's G there), the alphabetically later of the two at the
    five T columns and in the PAM. Above CRISPRME_IUPAC_CAP (10) variants a PUTATIVE
    sample's sets are its union, which was one per choice of ALT: 8,192 here, 794 of
    them within 4 mismatches and the PAM, each a row (NA19042, PUTATIVE in an IGL window
    of chr22 for a slot holding two ALTs, carries two or three at 16 of its 26
    positions: 8.5 million sets and 90 GB before any row). At most 16 now; past that one
    union, holding at each such column the ALT nearest the guide: one the PAM allows,
    then the fewest mismatches (each column's first ALT would put five mismatches in it,
    and an A in its PAM)."""

    COLS = list(range(2, 14)) + [22]
    # the guide's base -> (the reference base, the other ALT, the column's code)
    OF = {"A": ("C", "G", "V"), "C": ("A", "T", "H"), "G": ("A", "T", "D"),
          "T": ("A", "C", "H")}

    def _spec(self, c):
        """(the reference base, the ALT nearest the guide, the other ALT, the code)."""
        if c == 22:  # the PAM's last G: the reference C breaks the PAM, A too, G not
            return "C", "G", "A", "V"
        r, other, code = self.OF[SPACER[c]]
        return r, SPACER[c], other, code

    def test_one_union_nearest_the_guide(self):
        ref = _site(**{"c%d" % c: self._spec(c)[0] for c in self.COLS})
        values = {}
        for c in self.COLS:
            r, near, other, _code = self._spec(c)
            values[_pos1(c)] = [(r, near, "rs%dn" % c, "0.1", {"S": "1/0"}),
                                (r, other, "rs%do" % c, "0.1", {"S": "0/1"})]
        window = _window(ref, {c: self._spec(c)[3] for c in self.COLS})
        rows, _comp = self.analyse(ref, [window], values)
        self.assert_no_row_names_two_alts_of_a_position(rows)
        got = self.variant_rows(rows)
        nearest = tuple("chrT_%d_%s_%s" % ((_pos1(c),) + self._spec(c)[:2])
                        for c in self.COLS)
        self.assertEqual(list(got), [nearest])
        self.assertEqual((got[nearest][DNA], got[nearest][SAMPLES_COL]), (SITE, "S"))


class TheSnpSnpCompanionComparesPositions(unittest.TestCase):
    """_record_snp_snp_cooc counted distinct variant tokens: a row of two alts of one
    position was a co-occurrence."""

    def test_two_alts_of_one_position_are_no_cooccurrence(self):
        from test_phased_haplotype import _load_pure_functions
        ns = _load_pure_functions({"myreg": object()})
        line = ["X", CRRNA, SITE, "chrT", "100", "100", "+", "0", "0", "0", "NA", "y",
                "S", "NA", "NA", "rsA,rsT", "0.2,0.1", "chrT_106_C_A,chrT_106_C_T"]
        ns["_record_snp_snp_cooc"](line, "PUTATIVE")
        self.assertEqual(ns["_snp_snp_cooc_rows"], [])
        line[17] = "chrT_106_C_A,chrT_110_A_G"
        ns["_record_snp_snp_cooc"](line, "PUTATIVE")
        self.assertEqual(len(ns["_snp_snp_cooc_rows"]), 1)
        line[17] = "chrUn_KI270302v1_106_C_A,chrUn_KI270302v1_106_C_T"
        ns["_record_snp_snp_cooc"](line[:3] + ["chrUn"] + line[4:], "PUTATIVE")
        self.assertEqual(len(ns["_snp_snp_cooc_rows"]), 1)


class TheFourAlleleCodeThroughThePostAnalysis(_SnpPostAnalysis):
    """Column 5 holds C (reference), A, G and T: CRISPRitz 2.10.0 writes Q there. A
    target whose only code was Q was taken for a reference row: no carrier, Q in its
    DNA column."""

    REF = _site(c5="C")
    VALUES = {_pos1(5): [("C", "A", "rsA", "0.2", {"S": "1|0"}),
                         ("C", "G", "rsG", "0.1", {"U": "0|1"}),
                         ("C", "T", "rsT", "0.1", {"V": "1|1"})]}

    def _check(self, rows):
        self.assertFalse([r for r in rows if "Q" in r[DNA].upper()], rows)
        got = self.variant_rows(rows)
        p = _pos1(5)
        self.assertEqual({k: got[k][SAMPLES_COL] for k in got},
                         {("chrT_%d_C_A" % p,): "S", ("chrT_%d_C_G" % p,): "U",
                          ("chrT_%d_C_T" % p,): "V"})
        self.assertEqual(got[("chrT_%d_C_A" % p,)][DNA], SITE)

    def test_on_the_tiers(self):
        rows, _comp = self.analyse(self.REF, [_window(self.REF, {5: "Q"})], self.VALUES)
        self._check(rows)

    def test_on_the_dictionary(self):
        rows, _comp = self.analyse(self.REF, [_window(self.REF, {5: "Q"})], self.VALUES,
                                   tiers=False)
        self._check(rows)

    def test_population_level(self):
        rows, _comp = self.analyse(self.REF, [_window(self.REF, {5: "Q"})], self.VALUES,
                                   mode="1")
        got = self.variant_rows(rows)
        self.assertIn(("chrT_%d_C_A" % _pos1(5),), got)

    def test_an_alt_no_sample_carries(self):
        """Column 5 holds C (reference), A (S), G (a PASS record nobody carries) and T (V):
        CRISPRitz codes G all the same, so the window holds Q. The Tier-0 registry holds no
        record of G, and on the tiers alone (the published layout) the position read three
        bases: the Q column admitted no alt, and S and V were listed for column 9's G on a
        haplotype holding the reference C at column 5 (HPRC chr22:17,877,081 T>A, 8
        copies, T>C none, T>G 46: no row of A or G there, and G's carriers listed on
        haplotypes holding T). Every layout and mode now reads it as the dictionary does."""
        ref = _site(c5="C", c9="A")
        values = {_pos1(5): [("C", "A", "rsA", "0.2", {"S": "1|0"}),
                             ("C", "G", "rsG", "0", {}),
                             ("C", "T", "rsT", "0.1", {"V": "1|1"})],
                  _pos1(9): [("A", "G", "rs9", "0.1", {"S": "1|0", "V": "1|0"})]}
        a5, t5 = "chrT_%d_C_A" % _pos1(5), "chrT_%d_C_T" % _pos1(5)
        g9 = "chrT_%d_A_G" % _pos1(9)
        for mode in ("0", "1"):
            for tiers, keep_dict in ((True, False), (True, True), (False, False)):
                with self.subTest(mode=mode, tiers=tiers, keep_dict=keep_dict):
                    self.tmp = tempfile.mkdtemp(prefix="sa_q0_")
                    self.addCleanup(shutil.rmtree, self.tmp, True)
                    rows, _comp = self.analyse(ref, [_window(ref, {5: "Q", 9: "R"})],
                                               values, tiers=tiers, keep_dict=keep_dict,
                                               mode=mode)
                    got = self.variant_rows(rows)
                    self.assertNotIn((g9,), got)  # no carrier holds C at column 5 with G
                    if mode == "0":
                        self.assertEqual({k: got[k][SAMPLES_COL] for k in got},
                                         {(a5, g9): "S", (t5, g9): "V", (t5,): "V"})
                    else:
                        self.assertEqual(set(got[(a5, g9)][SAMPLES_COL].split(",")),
                                         {"S", "V"})

    def test_an_alt_that_is_no_base(self):
        """The dictionary of CRISPRitz 2.10.0 lists a one-character ALT that is no base,
        such as the '*' of a spanning deletion (a GATK joint call's). The Q column
        admitted it with the position's ALTs: a row of X, '*' in its DNA column, CFD 0
        (every other code's column drops it). The tiers hold no such ALT."""
        values = {_pos1(5): self.VALUES[_pos1(5)] + [("C", "*", "rsS", "0.1", {"X": "0|1"})]}
        for mode in ("0", "1"):
            for tiers, keep_dict in ((True, True), (False, False)):
                with self.subTest(mode=mode, tiers=tiers, keep_dict=keep_dict):
                    self.tmp = tempfile.mkdtemp(prefix="sa_qstar_")
                    self.addCleanup(shutil.rmtree, self.tmp, True)
                    rows, _comp = self.analyse(self.REF, [_window(self.REF, {5: "Q"})],
                                               values, tiers=tiers, keep_dict=keep_dict,
                                               mode=mode)
                    self.assertFalse([r for r in rows if "*" in r[DNA] + r[SNP]], rows)
                    if mode == "0":
                        self._check(rows)

    def test_a_q_column_at_a_position_of_three_bases_admits_none(self):
        """#139 at a Q column: Q permits every base, so it cannot tell its position's
        alts from a neighbour's; it admits them only where they and the reference hold
        all four bases. Here column 5's position holds C, A and T (the Q of the hit is
        another position's, as at a bulge boundary): no row names an alt there."""
        values = {_pos1(5): [("C", "A", "rsA", "0.2", {"S": "1|0"}),
                             ("C", "T", "rsT", "0.1", {"V": "1|1"})],
                  _pos1(9): [("A", "G", "rsG", "0.1", {"S": "1|0"})]}
        ref = _site(c5="C", c9="A")
        rows, _comp = self.analyse(ref, [_window(ref, {5: "Q", 9: "R"})], values)
        for r in rows:
            self.assertNotIn("_%d_" % _pos1(5), r[SNP], r)
        self.assertIn(("chrT_%d_A_G" % _pos1(9),), self.variant_rows(rows))

    def test_a_window_whose_only_column_admits_none(self):
        """The same Q, the window's only variant column, on the registry alone (no
        dictionary, no Tier-1 store), where a row needs no carrier: population-level, the
        post-analysis stopped (IndexError) on the entry naming no variant."""
        values = {_pos1(5): [("C", "A", "rsA", "0.2", {"S": "1|0"}),
                             ("C", "T", "rsT", "0.1", {"V": "1|1"})]}
        ref = _site(c5="C")
        for mode in ("0", "1"):
            with self.subTest(mode=mode):
                self.tmp = tempfile.mkdtemp(prefix="sa_q3_")
                self.addCleanup(shutil.rmtree, self.tmp, True)
                rows, _comp = self.analyse(ref, [_window(ref, {5: "Q"})], values,
                                           mode=mode, genotypes=False)
                self.assertEqual(self.variant_rows(rows), {})


class TheDictionaryPathReadsTheColumn(_SnpPostAnalysis):
    """#139 on the dictionary path too: it substituted every alt the dictionary lists at
    a position, also one the target's column does not hold, where the observed path
    keeps only those the column holds; the two paths now agree."""

    def test_an_alt_the_column_does_not_hold_is_no_allele_of_the_target(self):
        ref = _site(c5="C")
        values = {_pos1(5): [("C", "A", "rsA", "0.2", {"O": "1|0"}),
                             ("C", "T", "rsT", "0.1", {"P": "0|1"})]}
        for tiers in (False, True):
            with self.subTest(tiers=tiers):
                self.tmp = tempfile.mkdtemp(prefix="sa_col_")
                self.addCleanup(shutil.rmtree, self.tmp, True)
                rows, _comp = self.analyse(ref, [_window(ref, {5: "M"})], values,
                                           tiers=tiers)
                self.assertEqual(set(self.variant_rows(rows)),
                                 {("chrT_%d_C_A" % _pos1(5),)})


class ALineOfSeveralAltsThroughThePostAnalysis(_SnpPostAnalysis):
    """A VCF line C>A,T at column 5 with S 1|2 (A on slot 0, T on slot 1) and Y 2|2, and
    A>G at column 9 with S 0|1: S's slot 1 holds T and G. Read against "1", S's 1|2 put
    A on both slots: the dictionary path combined A with G for S. The tiers counted Y's
    two copies of T as reference (T's AF 1/16, not 3/16), and the observed path's guess
    of A's allele index from A's own carriers (S alone, 1|2: a tie, taken as 2) put A
    on T's slot, a row naming both."""

    REF = _site(c5="C", c9="A")
    VALUES = {_pos1(5): [("C", "A", "rs1", "0.25", {"S": "1|2"}),
                         ("C", "T", "rs1", "0.25", {"S": "1|2", "Y": "2|2"})],
              _pos1(9): [("A", "G", "rs9", "0.1", {"S": "0|1"})]}

    def _check(self, rows):
        self.assert_no_row_names_two_alts_of_a_position(rows)
        got = self.variant_rows(rows)
        p5, p9 = _pos1(5), _pos1(9)
        a, t, g = "chrT_%d_C_A" % p5, "chrT_%d_C_T" % p5, "chrT_%d_A_G" % p9
        self.assertEqual(got[(t, g)][SAMPLES_COL], "S")
        self.assertNotIn("S", got.get((a, g), [""] * 18)[SAMPLES_COL].split(","))
        self.assertEqual(got[(a,)][SAMPLES_COL], "S")
        self.assertIn("Y", got[(t,)][SAMPLES_COL].split(","))
        return got

    def test_on_the_dictionary(self):
        self._check(self.analyse(self.REF, [_window(self.REF, {5: "H", 9: "R"})],
                                 self.VALUES, tiers=False)[0])

    def test_on_the_tiers(self):
        rows, comp = self.analyse(self.REF, [_window(self.REF, {5: "H", 9: "R"})],
                                  self.VALUES)
        got = self._check(rows)
        t = "chrT_%d_C_T" % _pos1(5)
        g = "chrT_%d_A_G" % _pos1(9)
        # T's frequency from the registry: S's copy and Y's two, of 16
        self.assertEqual(got[(t, g)][AF].split(",")[0], "%.6g" % (3 / 16.0))
        phases = {row["SNP"]: row["Phase_Confirmation"] for row in comp["phase_confirmation"]}
        self.assertEqual(phases[t + "," + g], "CONFIRMED")

    def test_on_a_store_of_an_earlier_release(self):
        """A Tier-1 store compiled by 2.8.1 holds the raw 1|2 under each ALT of the line:
        the observed path reads it per ALT, as the store this release compiles."""
        rows, comp = self.analyse(self.REF, [_window(self.REF, {5: "H", 9: "R"})],
                                  self.VALUES, store_of_2_8_1=True)
        self._check(rows)
        phases = {row["SNP"]: row["Phase_Confirmation"] for row in comp["phase_confirmation"]}
        self.assertEqual(phases["chrT_%d_C_T,chrT_%d_A_G" % (_pos1(5), _pos1(9))],
                         "CONFIRMED")


def _mega_vcf(path, records):
    """A sites-only VCF as the mega's: ``records`` of (1-based position, REF, ALT, ID,
    FILTER), each with an ``AF_HGDP``."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n"
                 '##INFO=<ID=AF_HGDP,Number=A,Type=Float,Description="AF">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for p, ref, alt, rs, flt in records:
            fh.write(f"chrT\t{p}\t{rs}\t{ref}\t{alt}\t.\t{flt}\tAF_HGDP=0.1\n")
    return path


class TheSitesOnlyUnionTakesOneAltAColumn(_SnpPostAnalysis):
    """A sites-only dataset's window: column 5 holds C>A and C>T (split records), column
    9 A>C and A>G. The co-located union applied every alt of a column (the sequence
    holding the later, the row naming both): one union per choice of alt now, each
    choice its own row, though the greedy representative gives one of them (A and G,
    the guide's bases) and the others only their own unions."""

    def test_one_union_per_choice(self):
        vcf = _mega_vcf(os.path.join(self.tmp, "mega.chrT.vcf"), (
            (_pos1(5), "C", "A", "rsA", "PASS"), (_pos1(5), "C", "T", "rsT", "PASS"),
            (_pos1(9), "A", "C", "rsC", "PASS"), (_pos1(9), "A", "G", "rsG", "PASS")))
        ref = _site(c5="C", c9="A")
        rows, _comp = self.analyse(ref, [_window(ref, {5: "H", 9: "V"})], None,
                                   sites_only_vcf=vcf, mode="0", mms=2)
        self.assert_no_row_names_two_alts_of_a_position(rows)
        p5, p9 = _pos1(5), _pos1(9)
        self.assertEqual(set(self.variant_rows(rows)),
                         {("chrT_%d_C_%s" % (p5, a5), "chrT_%d_A_%s" % (p9, a9))
                          for a5 in "AT" for a9 in "CG"})


class TheSitesOnlyPathTakesEveryAltItsRegistryLists(_SnpPostAnalysis):
    """A sites-only (mega) registry is read from its VCF, every record of it: it lists
    the ALTs an index built by CRISPRitz 2.9.1 or earlier never coded (the earlier
    records of a split site; chr22:10,752,820 C>A, whose later record C>T is the one
    coded) and non-PASS records. 2.8.1 substituted each at its position's column, and so
    does this release; with the dictionary path's column filter, a search of the
    published mega index lost them (on chr22, for five guides, the highest CFD of 114
    sites fell: TRAC's chr22:22,900,049 from 0.524 to 0.127)."""

    def test_an_alt_the_enrichment_did_not_code(self):
        # column 5: C>A then C>T, coded Y (the later record's alone); column 9: A>C, coded
        # M, and A>G, a non-PASS record CRISPRitz never codes
        vcf = _mega_vcf(os.path.join(self.tmp, "mega.chrT.vcf"), (
            (_pos1(5), "C", "A", "rsA", "PASS"), (_pos1(5), "C", "T", "rsT", "PASS"),
            (_pos1(9), "A", "C", "rsC", "PASS"), (_pos1(9), "A", "G", "rsG", "LowQual")))
        ref = _site(c5="C", c9="A")
        rows, _comp = self.analyse(ref, [_window(ref, {5: "Y", 9: "M"})], None,
                                   sites_only_vcf=vcf, mms=2)
        p5, p9 = _pos1(5), _pos1(9)
        got = self.variant_rows(rows)
        # the guide's own bases: the worst case of the site
        self.assertIn(("chrT_%d_C_A" % p5, "chrT_%d_A_G" % p9), got)
        self.assertEqual(set(got), {("chrT_%d_C_%s" % (p5, a5), "chrT_%d_A_%s" % (p9, a9))
                                    for a5 in "AT" for a9 in "CG"})

    def test_a_q_column(self):
        """A Q column admits every ALT the registry lists at its position, whether or not
        they and the reference hold four bases (the registry drops an ALT with no
        frequency; a Q at a bulge boundary can be a neighbouring position's): each is
        the panel's own allele at its own position."""
        vcf = _mega_vcf(os.path.join(self.tmp, "mega.chrT.vcf"), (
            (_pos1(5), "C", "A", "rsA", "PASS"), (_pos1(5), "C", "T", "rsT", "PASS")))
        ref = _site(c5="C")
        rows, _comp = self.analyse(ref, [_window(ref, {5: "Q"})], None, sites_only_vcf=vcf)
        self.assertEqual(set(self.variant_rows(rows)),
                         {("chrT_%d_C_A" % _pos1(5),), ("chrT_%d_C_T" % _pos1(5),)})


# --------------------------------------------------------------------------- #
# The merge of a cluster's rows: the real remove_contiguous_samples.py
# --------------------------------------------------------------------------- #
class TheRsidAndAfListsFollowTheVariants(unittest.TestCase):
    """The merge sorted the variant, rsID and AF lists each on its own, as text, and
    kept each value once (the README guide's chr22:16858220 listed AF 0.000981162 for
    the variant of 0.805338 and the reverse; ACTGA had 953 of 1,967 multi-SNP rows out
    of order)."""

    HEAD = ("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\t"
            "Mismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\t"
            "Real_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster").split("\t")

    def _row(self, samples, rsids, afs, snps, cfd="0.9"):
        v = {"#Bulge_type": "X", "crRNA": CRRNA, "DNA": "CTAACAGTTGCaTTTATCACTGG",
             "Chromosome": "chr22", "Position": "16858220", "Cluster_Position": "16858220",
             "Direction": "+", "Mismatches": "1", "Bulge_Size": "0", "Total": "1",
             "PAM_gen": "n", "Var_uniq": "y", "Samples": samples, "Annotation_Type": "n",
             "Real_Guide": CRRNA, "rsID": rsids, "AF": afs, "SNP": snps,
             "Reference": SITE, "CFD_ref": "0.1", "CFD": cfd, "#Seq_in_cluster": "0"}
        return "\t".join(v[c] for c in self.HEAD)

    def _merge(self, rows):
        d = tempfile.mkdtemp(prefix="sa_merge_")
        self.addCleanup(shutil.rmtree, d, True)
        path = os.path.join(d, "x.bestCFD.txt")
        with open(path, "w") as fh:
            fh.write("\t".join(self.HEAD) + "\n" + "".join(r + "\n" for r in rows))
        subprocess.run([sys.executable, "adjust_cols.py", path], cwd=HERE, check=True,
                       capture_output=True, text=True)
        done = subprocess.run(
            [sys.executable, "remove_contiguous_samples.py", path, path + ".trimmed", "3",
             "5", "7", "11", "16", "19", "21", "score", "mm+bulges"],
            cwd=HERE, capture_output=True, text=True)
        self.assertEqual((done.returncode, done.stderr), (0, ""))
        with open(path + ".trimmed") as fh:
            best = [ln.rstrip("\n").split("\t") for ln in fh][1:]
        self.assertEqual(len(best), 1)
        # after adjust_cols.py: Samples 13, rsID 16, AF 17, SNP 18
        return best[0][16], best[0][17], best[0][18]

    def test_one_row(self):
        snps = "chr22_16858239_C_G,chr22_16858242_G_T"
        self.assertEqual(self._merge([self._row("S1", "rs9,rs1", "0.805338,0.000981162",
                                                snps)]),
                         ("rs9,rs1", "0.805338,0.000981162", snps))

    def test_two_variants_of_one_frequency_keep_it_each(self):
        snps = "chr22_16858239_C_G,chr22_16858242_G_T"
        self.assertEqual(self._merge([self._row("S1", "rs9,.", "0.5,0.5", snps)]),
                         ("rs9,.", "0.5,0.5", snps))

    def test_the_row_two_datasets_give(self):
        """Each dataset's rsID and AF stay on the row, by variant: (rsID, AF) pairs."""
        snps = "chr22_16858239_C_G,chr22_16858242_G_T"
        rsids, afs, got = self._merge([
            self._row("S1", "rs9,rs1", "0.805338,0.000981162", snps),
            self._row("S2", "rs9,.", "0.7,0.002", snps)])
        self.assertEqual((rsids, afs, got),
                         ("rs9,rs9,.,rs1", "0.7,0.805338,0.002,0.000981162", snps))
        self.assertEqual(self._merge([self._row("S1", "rs9,rs1", "0.8,0.001", snps),
                                      self._row("S2", "rs9,rs1", "0.8,0.001", snps)]),
                         ("rs9,rs1", "0.8,0.001", snps))

    def test_an_rsid_holding_two_ids(self):
        snps = "chr22_16858239_C_G,chr22_16858242_G_T"
        self.assertEqual(self._merge([self._row("S1", "rs9;rs8,rs1", "0.8,0.001", snps)]),
                         ("rs9;rs8,rs1", "0.8,0.001", snps))


# --------------------------------------------------------------------------- #
# The sites-only builders read a line of several ALTs one ALT at a time
# --------------------------------------------------------------------------- #
class TheMegaIndelAfOfALineOfSeveralAlts(unittest.TestCase):
    def test_each_indel_alt_its_own_af(self):
        d = tempfile.mkdtemp(prefix="sa_megaindel_")
        self.addCleanup(shutil.rmtree, d, True)
        vcf = os.path.join(d, "mega.chrT.vcf.gz")
        with gzip.open(vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("chrT\t100\t.\tCA\tC,CAA,G\t.\tPASS\t"
                     "AF_HGDP=0.1,0.2,0.3;AF_max=0.1,0.2,0.3\n")
        got = list(build_mega_indel_af.iter_indel_af_records(vcf, {"HGDP"}, "chrT"))
        self.assertEqual(got, [("100", "CA", "C", {"HGDP": 0.1}, 0.1),
                               ("100", "CA", "CAA", {"HGDP": 0.2}, 0.2),
                               ("100", "CA", "G", {"HGDP": 0.3}, 0.3)])


if __name__ == "__main__":
    unittest.main()

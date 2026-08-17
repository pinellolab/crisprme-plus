"""Tests for build_dictless_tiers (CRISPRme+ dictless redesign, Phase-3d BUILD).

These tests pin the #1 correctness property of the build-side tier emission: the
directory + file names it writes MUST EXACTLY match what the SEARCH-side resolvers
(``new_simple_analysis._resolve_registry_paths`` / ``_resolve_genotype_paths``)
read, or a freshly-built install silently ignores the tiers.

We do NOT import ``new_simple_analysis`` (it executes module-level code reading
sys.argv and imports numpy/pandas/CRISTA_score -- not available in the light CI).
Instead we:
  1. Re-derive the resolver-expected paths with an INDEPENDENT re-implementation of
     the resolver algorithm (``_resolver_registry_paths`` / ``_resolver_genotype_
     paths`` below, hand-copied from the resolver source), and assert the helper's
     emitted paths equal them byte-for-byte.
  2. Additionally GREP the resolver SOURCE for the load-bearing string literals
     ("registry_", "genotypes_", "reg_", "gt_", ".bin", ".idx", the
     "dictionaries_"-prefix swap) so a future edit to the resolvers that changes
     the naming trips this test rather than silently drifting from the emitter.

STDLIB ONLY (unittest, tempfile, os, json, gzip).
"""

import gzip
import json
import os
import tempfile
import unittest

import build_dictless_tiers as bdt
import tier0_registry as t0reg
import tier1_genotypes as t1gt


# --------------------------------------------------------------------------- #
# Independent re-implementation of the SEARCH-side resolvers (hand-copied from
# new_simple_analysis._resolve_registry_paths / _resolve_genotype_paths). If the
# emitter drifts from these, the emitted files land where the search side does NOT
# look and the tiers are silently ignored -- so these ARE the contract.
# --------------------------------------------------------------------------- #
def _resolver_registry_paths(dict_path, chrom):
    dict_dir = os.path.dirname(dict_path)
    parent = os.path.dirname(dict_dir)
    dict_folder_name = os.path.basename(dict_dir)
    if dict_folder_name.startswith("dictionaries_"):
        reg_folder_name = "registry_" + dict_folder_name[len("dictionaries_"):]
    else:
        reg_folder_name = "registry_" + dict_folder_name
    reg_dir = os.path.join(parent, reg_folder_name)
    return (os.path.join(reg_dir, "reg_" + str(chrom) + ".bin"),
            os.path.join(reg_dir, "reg_" + str(chrom) + ".idx"))


def _resolver_genotype_paths(dict_path, chrom):
    dict_dir = os.path.dirname(dict_path)
    parent = os.path.dirname(dict_dir)
    dict_folder_name = os.path.basename(dict_dir)
    if dict_folder_name.startswith("dictionaries_"):
        gt_folder_name = "genotypes_" + dict_folder_name[len("dictionaries_"):]
    else:
        gt_folder_name = "genotypes_" + dict_folder_name
    gt_dir = os.path.join(parent, gt_folder_name)
    return (os.path.join(gt_dir, "gt_" + str(chrom) + ".bin"),
            os.path.join(gt_dir, "gt_" + str(chrom) + ".idx"))


# --------------------------------------------------------------------------- #
# Fixture writers (mirror test_tier0_compile).
# --------------------------------------------------------------------------- #
G1000_ROWS = [
    ("S_EUR_M1", "GBR", "EUR", "male"),
    ("S_EUR_F1", "GBR", "EUR", "female"),
    ("S_AFR_M1", "YRI", "AFR", "male"),
    ("S_AFR_F1", "YRI", "AFR", "female"),
]
HGDP_ROWS = [
    ("H_EAS_M1", "Han", "EAS", "male"),
    ("H_EAS_F1", "Han", "EAS", "female"),
]
CHROM = "chr1"
VCF = "hg38_1000G_HGDP"   # combined-panel convention


def _write_samplesid(path, rows):
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
        for sid, pop, spop, sex in rows:
            fh.write("\t".join([sid, pop, spop, sex]) + "\n")


def _write_dict(path, mapping, *, gz=False):
    if gz:
        with gzip.open(path, "wt") as fh:
            json.dump(mapping, fh)
    else:
        with open(path, "w") as fh:
            json.dump(mapping, fh)


class _EmitFixture(unittest.TestCase):
    """Builds a tmp Dictionaries/dictionaries_<vcf>/ layout + samplesID files."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        # <root>/Dictionaries/dictionaries_<vcf>/my_dict_<chrom>.json
        self.dictionaries = os.path.join(self.root, "Dictionaries")
        self.dict_folder = os.path.join(self.dictionaries, "dictionaries_" + VCF)
        os.makedirs(self.dict_folder)
        # samplesID files (two databases -> combined panel).
        self.sids = os.path.join(self.root, "samplesIDs")
        os.makedirs(self.sids)
        self.g1000_sid = os.path.join(self.sids, "hg38_1000G.samplesID.txt")
        self.hgdp_sid = os.path.join(self.sids, "hg38_HGDP.samplesID.txt")
        _write_samplesid(self.g1000_sid, G1000_ROWS)
        _write_samplesid(self.hgdp_sid, HGDP_ROWS)
        self.db_to_sid = {"1000G": self.g1000_sid, "HGDP": self.hgdp_sid}

        # A small dict with a multiallelic '$' SNP site + a plain SNP site.
        #   pos 100 : alt G carried het by S_EUR_F1 (1000G::EUR)
        #   pos 300 : alt G (S_EUR_F1 het) $ alt T (H_EAS_M1 het, HGDP::EAS)
        self.mapping = {
            "%s,100" % CHROM: "S_EUR_F1:0|1;A,G;rs100;0.05",
            "%s,300" % CHROM: (
                "S_EUR_F1:0|1;C,G;rsG300;0.1"
                "$"
                "H_EAS_M1:0|1;C,T;rsT300;0.02"
            ),
        }
        self.dict_path = os.path.join(self.dict_folder, "my_dict_%s.json" % CHROM)

    def _write_plain(self):
        _write_dict(self.dict_path, self.mapping, gz=False)
        return self.dict_path

    def _write_gz(self):
        gzp = self.dict_path + ".gz"
        _write_dict(gzp, self.mapping, gz=True)
        return gzp


# --------------------------------------------------------------------------- #
# (1) EMITTED PATHS == RESOLVER-EXPECTED PATHS (the #1 correctness property).
# --------------------------------------------------------------------------- #
class TestEmittedPathsMatchResolvers(_EmitFixture):
    def test_helper_paths_equal_resolver_paths(self):
        # Both for a plain dict path AND the .gz path -- the resolver keys off the
        # FOLDER name, not the file suffix, so both must map to the same tier dirs.
        for dpath in (self.dict_path, self.dict_path + ".gz"):
            got_reg = bdt.registry_paths_for(dpath, CHROM)
            got_gt = bdt.genotype_paths_for(dpath, CHROM)
            want_reg = _resolver_registry_paths(dpath, CHROM)
            want_gt = _resolver_genotype_paths(dpath, CHROM)
            self.assertEqual(got_reg, want_reg,
                             "registry path drift for %s" % dpath)
            self.assertEqual(got_gt, want_gt,
                             "genotype path drift for %s" % dpath)
        # And spell out the EXACT expected on-disk names for documentation.
        reg_bin, reg_idx = bdt.registry_paths_for(self.dict_path, CHROM)
        gt_bin, gt_idx = bdt.genotype_paths_for(self.dict_path, CHROM)
        self.assertEqual(
            os.path.relpath(reg_bin, self.root),
            os.path.join("Dictionaries", "registry_" + VCF, "reg_%s.bin" % CHROM))
        self.assertEqual(
            os.path.relpath(reg_idx, self.root),
            os.path.join("Dictionaries", "registry_" + VCF, "reg_%s.idx" % CHROM))
        self.assertEqual(
            os.path.relpath(gt_bin, self.root),
            os.path.join("Dictionaries", "genotypes_" + VCF, "gt_%s.bin" % CHROM))
        self.assertEqual(
            os.path.relpath(gt_idx, self.root),
            os.path.join("Dictionaries", "genotypes_" + VCF, "gt_%s.idx" % CHROM))

    def test_resolver_source_still_uses_the_pinned_names(self):
        # Guard against a future edit to the resolvers that renames the dirs/files:
        # grep the resolver SOURCE for the load-bearing literals. If someone changes
        # "registry_" -> "reg_dir_" in new_simple_analysis, this fails LOUDLY here
        # (and the emitter must be updated in lock-step).
        here = os.path.dirname(os.path.abspath(__file__))
        src_path = os.path.join(here, "new_simple_analysis.py")
        if not os.path.isfile(src_path):
            self.skipTest("new_simple_analysis.py not present in this checkout")
        with open(src_path) as fh:
            src = fh.read()
        # both resolvers must exist and swap the dictionaries_ prefix
        self.assertIn("def _resolve_registry_paths", src)
        self.assertIn("def _resolve_genotype_paths", src)
        self.assertIn('"registry_" + dict_folder_name[len("dictionaries_"):]', src)
        self.assertIn('"genotypes_" + dict_folder_name[len("dictionaries_"):]', src)
        # per-chromosome file naming
        self.assertIn('"reg_" + str(chrom) + ".bin"', src)
        self.assertIn('"reg_" + str(chrom) + ".idx"', src)
        self.assertIn('"gt_" + str(chrom) + ".bin"', src)
        self.assertIn('"gt_" + str(chrom) + ".idx"', src)


# --------------------------------------------------------------------------- #
# (2) EMISSION lands the files + they OPEN + a known (pos,alt) resolves.
# --------------------------------------------------------------------------- #
class TestEmitAndRead(_EmitFixture):
    def _run_and_check(self, dpath):
        res = bdt.emit_dictless_tiers(dpath, self.db_to_sid, CHROM,
                                      self.dict_folder)
        # emitted paths are exactly the resolver-expected paths
        self.assertEqual((res["registry_bin"], res["registry_idx"]),
                         _resolver_registry_paths(dpath, CHROM))
        self.assertEqual((res["genotype_bin"], res["genotype_idx"]),
                         _resolver_genotype_paths(dpath, CHROM))
        # all four files exist on disk in the sibling dirs
        for p in (res["registry_bin"], res["registry_idx"],
                  res["genotype_bin"], res["genotype_idx"]):
            self.assertTrue(os.path.exists(p), "missing emitted file %s" % p)

        # OPEN the Tier-0 registry and assert the known (pos, alt) global AF.
        reg = t0reg.RegistryReader(res["registry_bin"], res["registry_idx"])
        try:
            # pos 100 alt G: one het carrier (S_EUR_F1) over a 6-individual panel
            # (4 1000G + 2 HGDP), all diploid autosome -> AN 12, AC 1 -> AF 1/12.
            der = reg.derived(100, "G")
            self.assertIsNotNone(der, "registry has no (100,G) record")
            self.assertTrue(der["observed"])
            self.assertAlmostEqual(der["global_af"], 1.0 / 12.0)
            # multiallelic site: both alts present & independent
            self.assertIsNotNone(reg.lookup(300, "G"))
            self.assertIsNotNone(reg.lookup(300, "T"))
            self.assertEqual(reg.rsid(300, "T"), "rsT300")
        finally:
            reg.close()

        # OPEN the Tier-1 genotype store and assert carrier_tokens round-trip.
        gt = t1gt.GenotypeReader(res["genotype_bin"], res["genotype_idx"])
        try:
            toks = gt.carrier_tokens(100, "G")
            self.assertIsNotNone(toks, "genotype store has no (100,G) record")
            self.assertTrue(toks, "carrier_tokens empty for (100,G)")
            self.assertIn("S_EUR_F1:0|1", toks)
            # multiallelic alt T at 300 carried by the HGDP EAS male
            toksT = gt.carrier_tokens(300, "T")
            self.assertIsNotNone(toksT)
            self.assertIn("H_EAS_M1:0|1", toksT)
        finally:
            gt.close()

    def test_emit_from_plain_dict(self):
        self._run_and_check(self._write_plain())

    def test_emit_from_gzipped_dict(self):
        # gz DECISION check: the emitter reads a .json.gz dict transparently (the
        # compilers are gz-aware), so a build that already gzipped its dicts still
        # emits usable tiers. Same assertions as the plain case.
        self._run_and_check(self._write_gz())

    def test_emit_finds_gz_when_asked_for_plain(self):
        # The build gzips in place; a caller that hands the PLAIN name after gzip
        # must still find the .gz and emit. (Emitter's _resolve_dict_path fallback.)
        self._write_gz()  # only the .gz exists
        self.assertFalse(os.path.exists(self.dict_path))
        res = bdt.emit_dictless_tiers(self.dict_path, self.db_to_sid, CHROM,
                                      self.dict_folder)
        self.assertTrue(os.path.exists(res["registry_bin"]))
        self.assertTrue(os.path.exists(res["genotype_bin"]))


# --------------------------------------------------------------------------- #
# (3) GUARD: a bad/missing dict raises in the plain helper, but the GUARDED
#     wrapper logs + returns None WITHOUT crashing the caller.
# --------------------------------------------------------------------------- #
class TestGuard(_EmitFixture):
    def test_missing_dict_raises_in_plain_helper(self):
        missing = os.path.join(self.dict_folder, "my_dict_chrNOPE.json")
        with self.assertRaises(FileNotFoundError):
            bdt.emit_dictless_tiers(missing, self.db_to_sid, "chrNOPE",
                                    self.dict_folder)

    def test_guarded_wrapper_swallows_missing_dict(self):
        missing = os.path.join(self.dict_folder, "my_dict_chrNOPE.json")
        # Must NOT raise; returns None (the caller's build continues).
        res = bdt.emit_dictless_tiers_guarded(missing, self.db_to_sid, "chrNOPE",
                                              self.dict_folder)
        self.assertIsNone(res)
        # and it did NOT create any registry/genotype dir/files for chrNOPE
        reg_bin, _ = bdt.registry_paths_for(missing, "chrNOPE")
        self.assertFalse(os.path.exists(reg_bin))

    def test_guarded_wrapper_swallows_bad_dict_content(self):
        # A corrupt (non-JSON) dict body -> compile raises -> guarded returns None.
        badp = os.path.join(self.dict_folder, "my_dict_chrBAD.json")
        with open(badp, "w") as fh:
            fh.write("this is not json {{{")
        res = bdt.emit_dictless_tiers_guarded(badp, self.db_to_sid, "chrBAD",
                                              self.dict_folder)
        self.assertIsNone(res)


if __name__ == "__main__":
    unittest.main()

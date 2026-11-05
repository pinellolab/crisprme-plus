"""Tier-path SNP rows carry the REF alignment the dictionary path gives them.

BUG (2.8.1, pre-existing since the dict-less tiers)
---------------------------------------------------
A search that reads the dict-less tiers (a Hugging Face tiered index, a
``build-index-only --samplesID`` build, a Settings variant build) finalizes its SNP
variant rows in ``new_simple_analysis._iupac_decomposition_observed``. Each such row
carries the candidate's REFERENCE alignment in its Reference column (target[-3],
``refSeq_with_bulges``). The dictionary-path block of ``iupac_decomposition``
lower-cases every protospacer base of that alignment that differs from the guide;
the per-sample (observed-haplotype) block and the population-level (fast) block of
the tier path did not. ``resultIntegrator.seed_processing`` counts lower case (and
'-') to fill Seed_/Non_seed_mismatches+bulges_REF, and CRISPRme_plots places the
REF point of the fewest-mismatches plot from the same lower case, so on every tier-
path SNP row Aligned_protospacer+PAM_REF had no lower case and the REF seed /
non-seed counts left out the mismatches (they counted only the bulges; the
fewest_mm+b block, only the RNA bulges, see test_seed_counts_dna_bulge).

FIX
---
Both tier-path blocks build the Reference column with ``_ref_with_bulges``, which
does what the dictionary-path block does: re-insert the RNA-bulge gaps, lower-case
each protospacer base that differs from the guide, leave a DNA-bulge column (guide
'-') and the PAM bases as they are.

WHAT THIS TEST ASSERTS
----------------------
For six hand-built candidates -- plus strand with mismatches, an RNA bulge, a DNA
bulge, the minus strand, a PAM-creating SNP, and a 5' PAM -- driving the REAL
``iupac_decomposition`` (the AST harness of ``test_phased_haplotype``) on the
dictionary path (``mygt`` None) and on the tier path (``mygt`` present), each in
per-sample mode and in population-level mode (``_FAST_MODE``):
  * every variant row of the tier path carries the expected REF alignment (lower
    case at the mismatches only; PAM and DNA-bulge columns upper case) -- this fails
    on aa9c54b for both tier-path blocks, which wrote the alignment all upper case;
  * the dictionary path gives the same string (it is the reference behaviour, and
    the test pins it so the two paths cannot drift apart again);
  * the tier path's REF alignment equals the dictionary path's, row for row;
  * the PAM-creation check that reads this string still reports the created PAM;
  * resultIntegrator's own ``seed_processing`` (AST-extracted, so no annotation or
    VCF inputs are needed) counts the expected REF seed / non-seed mismatches from
    the tier-path string, and the same counts as from the dictionary-path string.

Run with:
    cd PostProcess && python3 -m unittest test_tier_ref_columns -v
"""

import ast
import io
import os
import sys
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions

SPACER = "CTAACAGTTGCTTTTATCAC"
GUIDE_3P = SPACER + "NNN"  # 20 nt + NRG, as CRISPRitz writes the crRNA column

# NRG at the 3' end: the runtime globals new_simple_analysis.py derives from the PAM file
GEOM_3P = dict(pam="NRG", pos_beg=0, pos_end=-3, pam_begin=-3, pam_end=None)
# TTTV at the 5' end (Cas12a-like): pos_beg/pam_end move to the front
GEOM_5P = dict(pam="TTTV", pos_beg=4, pos_end=None, pam_begin=0, pam_end=4)

IUPAC_OF = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
            frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}


class _GtSentinel(object):
    """Truthy stand-in for the Tier-1 GenotypeReader: the tier path is chosen on
    ``mygt is not None``; carriers come from the dict entries via retrieveFromDict."""

    def __len__(self):
        return 1


def _load_seed_processing():
    """resultIntegrator.seed_processing, AST-extracted (the module itself reads its
    inputs from sys.argv at import time)."""
    src_path = os.path.join(PP, "resultIntegrator.py")
    with open(src_path) as fh:
        tree = ast.parse(fh.read())
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "seed_processing")
    ns = {}
    exec(compile(ast.Module(body=[node], type_ignores=[]), src_path, "exec"), ns)
    return ns["seed_processing"]


SEED_PROCESSING = _load_seed_processing()


def _candidate(name):
    """(split, genome, mydict, geometry, expected REF alignment, extra checks).

    ``window`` is the forward-strand reference; ``snps`` are (offset in window, ref,
    alt, {sample: GT}). The DNA column is the window with each SNP as its IUPAC code
    (reverse-complemented for the minus strand) and RNA-bulge gaps inserted, exactly
    as CRISPRitz reports a variant-genome hit."""
    rc = None
    gap_at = None          # RNA bulge: '-' in the DNA column at this index
    crrna = GUIDE_3P
    bulge_type, bulges = "X", "0"
    strand = "+"
    geom = GEOM_3P
    extra = {}
    if name == "plus_mismatches":
        window = "CTTACTGTTGCCTTTACCACTGG"
        snps = [(5, "T", "A", {"S1": "1|0", "S2": "0|1"}),
                (11, "C", "T", {"S1": "1|0", "S3": "1|1"})]
        expected = "CTtACtGTTGCcTTTAcCACTGG"
    elif name == "plus_rna_bulge":
        window = "TTACTGTTGCCTTTACCACTGG"
        snps = [(4, "T", "A", {"S1": "1|0"}), (10, "C", "T", {"S1": "1|0"})]
        gap_at = 2
        bulge_type, bulges = "RNA", "1"
        expected = "tT-ACtGTTGCcTTTAcCACTGG"
    elif name == "plus_dna_bulge":
        # extra genomic A at window[11] (the guide has '-' there); window[5] T is a
        # reference mismatch that the SNP T>A removes
        window = "CTAACTGTTGCATTTTATCACTGG"
        snps = [(5, "T", "A", {"S1": "1|1"})]
        crrna = "CTAACAGTTGC-TTTTATCACNNN"
        bulge_type, bulges = "DNA", "1"
        expected = "CTAACtGTTGCATTTTATCACTGG"  # bulge column A stays upper case
    elif name == "minus_strand":
        window = "CCCGTGATAAAAGCAAATGTTAG"  # target strand: CTAACATTTGCTTTTATCACGGG
        snps = [(10, "A", "G", {"S1": "0|1", "S3": "1|1"}),
                (16, "A", "C", {"S2": "1|0"})]
        strand = "-"
        rc = True
        expected = "CTAACAtTTGCTTTTATCACGGG"
    elif name == "pam_creation":
        window = "CTAACAGTTGCTTTTATCTCTGC"  # PAM TGC is not NRG; C>G makes TGG
        snps = [(22, "C", "G", {"S4": "0|1", "S5": "1|1"})]
        expected = "CTAACAGTTGCTTTTATCtCTGC"
        extra["pam_gen"] = "TGG"
    elif name == "five_prime_pam":
        spacer = "ACGTTGCAGTCAGTACGTCA"
        crrna = "NNNN" + spacer
        window = "TTTA" + "ACCTTGCAGTCAGGACGTCA"  # mismatches at spacer 2 and 13
        snps = [(17, "G", "T", {"S1": "1|0"})]    # removes the spacer-13 mismatch
        geom = GEOM_5P
        expected = "TTTAACcTTGCAGTCAGgACGTCA"
    else:
        raise KeyError(name)

    dna = list(window)
    mydict = {}
    for off, ref, alt, gts in snps:
        assert window[off] == ref, (name, off)
        dna[off] = IUPAC_OF[frozenset(ref + alt)]
        tokens = ",".join("%s:%s" % (s, g) for s, g in sorted(gts.items()))
        mydict["chrT,%d" % (off + 1)] = "%s;%s,%s;rs%s%d;0.25" % (tokens, ref, alt, name, off)
    dna = "".join(dna)
    if rc:
        tab = str.maketrans("ACGTRYSWKM", "TGCAYRSWMK")
        dna = dna.translate(tab)[::-1]
    if gap_at is not None:
        dna = dna[:gap_at] + "-" + dna[gap_at:]
    split = [bulge_type, crrna, dna, "chrT", "0", "0", strand, "0", bulges, bulges,
             "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
    return split, window, mydict, geom, expected, extra


CANDIDATES = ["plus_mismatches", "plus_rna_bulge", "plus_dna_bulge", "minus_strand",
              "pam_creation", "five_prime_pam"]


def _run(name, tier, population_level):
    split, genome, mydict, geom, expected, extra = _candidate(name)
    overrides = dict(
        genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
        mygt=_GtSentinel() if tier else None, haplotype_check=True, IUPAC_CAP=10,
        hvdr_bed=io.StringIO(), allowed_mms=6, **geom)
    ns = _load_pure_functions(overrides)
    ns["_FAST_MODE"] = population_level
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    guide = split[1]
    guide_no_pam = guide[geom["pos_beg"]:geom["pos_end"]]
    cluster = []
    ns["iupac_decomposition"](list(split), guide.replace("-", ""), guide_no_pam, cluster)
    # variant rows carry the 33 sentinel and the Reference alignment just before it
    return [r for r in cluster if r[-2] == 33], expected, extra, geom


def _seed_counts(ref, geom):
    """REF seed / non-seed counts, sliced as resultIntegrator slices them for a 3' PAM
    (or a 5' PAM), then counted by its own seed_processing."""
    n = geom["pam_end"] if geom["pam_begin"] == 0 else -geom["pam_begin"]
    if geom["pam_begin"] == 0:          # PAM at the start
        real = ref[n:]
        seed, non_seed = real[: int(len(real) / 2)], real[int(len(real) / 2):]
    else:
        real = ref[: len(ref) - n]
        seed, non_seed = real[int(len(real) / 2):], real[: int(len(real) / 2)]
    out = SEED_PROCESSING(seed, "", non_seed, "")
    return out[0], out[2]


EXPECTED_SEED = {  # (REF seed, REF non-seed) from the lower case + '-' of `expected`
    "plus_mismatches": (2, 2), "plus_rna_bulge": (2, 3), "plus_dna_bulge": (0, 1),
    "minus_strand": (0, 1), "pam_creation": (1, 0), "five_prime_pam": (1, 1),
}


class TierRefColumnMatchesDictionaryPath(unittest.TestCase):

    def _check(self, name, population_level):
        mode = "population-level" if population_level else "per-sample"
        legacy, expected, extra, geom = _run(name, tier=False, population_level=population_level)
        tier, _, _, _ = _run(name, tier=True, population_level=population_level)
        self.assertTrue(legacy, "%s/%s: the dictionary path emitted no variant row" % (name, mode))
        self.assertTrue(tier, "%s/%s: the tier path emitted no variant row" % (name, mode))
        for r in legacy:
            self.assertEqual(r[-3], expected,
                             "%s/%s: dictionary-path REF alignment changed" % (name, mode))
        for r in tier:
            self.assertEqual(
                r[-3], expected,
                "%s/%s: tier-path REF alignment %r lacks the mismatch lower case of the "
                "dictionary path (%r)" % (name, mode, r[-3], expected))
        # tier == dictionary, row for row (same haplotype DNA -> same REF alignment)
        by_dna = {r[2]: r[-3] for r in legacy}
        for r in tier:
            if r[2] in by_dna:
                self.assertEqual(r[-3], by_dna[r[2]], "%s/%s: %s" % (name, mode, r[2]))
        self.assertTrue(set(by_dna) & {r[2] for r in tier},
                        "%s/%s: no haplotype in common to compare" % (name, mode))
        # PAM positions stay upper case (the PAM-creation lookup reads them)
        p = slice(geom["pam_begin"], geom["pam_end"])
        for r in tier:
            self.assertEqual(r[-3][p], r[-3][p].upper())
        if "pam_gen" in extra:
            for r in legacy + tier:
                self.assertEqual(r[10], extra["pam_gen"],
                                 "%s/%s: the created PAM is no longer reported" % (name, mode))
        # resultIntegrator's own counting on the REF alignment
        want = EXPECTED_SEED[name]
        for r in tier:
            self.assertEqual(_seed_counts(r[-3], geom), want,
                             "%s/%s: REF seed/non-seed counts" % (name, mode))
        for r in legacy:
            self.assertEqual(_seed_counts(r[-3], geom), want)

    def test_per_sample_observed_block(self):
        for name in CANDIDATES:
            with self.subTest(candidate=name):
                self._check(name, population_level=False)

    def test_population_level_fast_block(self):
        for name in CANDIDATES:
            with self.subTest(candidate=name):
                self._check(name, population_level=True)

    def test_reference_row_unchanged(self):
        # the candidate's own REFERENCE row (sentinel 55, Reference 'n') is built by
        # _finalize_reference_entry, which already lower-cased its DNA: unchanged.
        split, genome, mydict, geom, expected, extra = _candidate("plus_mismatches")
        overrides = dict(
            genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
            mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
            hvdr_bed=io.StringIO(), allowed_mms=6, **geom)
        ns = _load_pure_functions(overrides)
        cluster = []
        ns["iupac_decomposition"](list(split), split[1], split[1][:-3], cluster)
        refs = [r for r in cluster if r[-2] == 55]
        self.assertEqual(len(refs), 1)
        self.assertEqual(refs[0][2], expected)
        self.assertEqual(refs[0][-3], "n")


class AMinusStrandCandidateGetsNoReferenceRow(unittest.TestCase):
    """The reference search holds its alignment. _finalize_reference_entry judged the
    plus-strand letters (a strand-corrected string reverse-complemented again), and on a
    near-palindromic guide they pass its gates: a row with strand '-' and the '+' letters,
    in the alternative alignments and one more in Other_motifs (TRAC guide, genome-wide,
    5 sites: round 5 review, RV3-A). The candidate is chr1:147073620 of that search."""

    WINDOW = "CTCTCCTCTGGTGCAGCTGAGAG"  # the forward strand there
    CRRNA = "CTCTCAGCTGGTACA-CGGCANNN"
    # a SNP G>A on the last forward base: the first base of the minus-strand target
    DNA = "YTCTCAGCTGcacCAGaGG-AGAG"

    def cluster(self, population_level, strand="-", dna=None):
        ns = _load_pure_functions(dict(
            genomeStr=self.WINDOW, current_chr="chrT", mydict={"chrT,23": "S1:1|0;G,A;rsX;0.25"},
            myreg=None, mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
            hvdr_bed=io.StringIO(), allowed_mms=4, **GEOM_3P))
        ns["_FAST_MODE"] = population_level
        split = ["RNA,DNA", self.CRRNA, dna or self.DNA, "chrT", "0", "0", strand, "4", "2", "6",
                 "NA", "y", "NA", "NA", "NA", "NA", "NA", "NA"]
        out = []
        ns["iupac_decomposition"](list(split), self.CRRNA.replace("-", ""), self.CRRNA[:-3], out)
        return out

    def test_in_either_mode(self):
        for population_level in (False, True):
            with self.subTest(population_level=population_level):
                rows = self.cluster(population_level)
                # (the plus-strand letters read CTCTCctCTGGTgCAGCtG-AGAG, within budget)
                self.assertEqual([r[2] for r in rows if r[-2] == 55], [])
                self.assertEqual({r[6] for r in rows} - {"-"}, set())


if __name__ == "__main__":
    unittest.main()

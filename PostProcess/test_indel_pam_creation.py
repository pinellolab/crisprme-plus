"""PAM-creation on the INDEL path (analisi_indels_NNN.py).

The block under test is inline in the analyzer's main loop, not a function, and
the module's top level reads sys.argv and opens the fake genome -- so, following
test_indel_snp_overlay.py, the source is extracted and exec'd against a synthetic
namespace instead of imported.

Why this file exists: the first version of the fix wrote the created PAM to
final_result[11] instead of [10]. The predicate itself was right, so every
standalone check of "which rows should be flagged" passed, and only a full
end-to-end run exposed it. These tests pin the *destination index* as well as the
predicate, so that class of error is caught in milliseconds.

Written as a unittest.TestCase so the project's `python -m unittest` CI harness
collects it (a prior pytest-style module-level `def test_*` layout collected zero
tests under unittest).
"""

import os
import unittest

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "analisi_indels_NNN.py")

# The real indel cluster row: 18 fields, NO Chromosome_fake column.
# 0 Bulge_type, 1 crRNA, 2 DNA(ALT target), 3 Chromosome, 4 Position,
# 5 Cluster Position, 6 Direction, 7 Mismatches, 8 Bulge_Size, 9 Total,
# 10 PAM_gen, 11 Var_uniq, 12 Samples, 13 Annotation Type, 14 Real Guide,
# 15 rsID, 16 AF, 17 indel info
PAM_GEN = 10
VAR_UNIQ = 11

IUPAC = {
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "A": {"A"}, "T": {"T"}, "C": {"C"}, "G": {"G"},
    "a": {"a"}, "t": {"t"}, "c": {"c"}, "g": {"g"},
    "N": {"A", "T", "G", "C"},
}


def _block():
    """Extract just the PAM-creation block (ref_pam assignment through the loop)."""
    src = open(SRC).read()
    start = src.index("    ref_pam = ref_aligned[pam_begin:pam_end]")
    end = src.index("    final_result.append(ref_aligned)", start)
    body = src[start:end]
    # dedent one level: the block lives inside a `for line in inResult:` loop
    return "\n".join(l[4:] if l.startswith("    ") else l
                     for l in body.splitlines())


def _row(alt_target):
    r = ["RNA,DNA", "GUIDE", alt_target, "chr1", "100", "100", "+",
         "2", "0", "2", "n", "n", "SAMPLE1", "n", "GUIDE", "rs1", "0.01",
         "chr1_100_A_AT"]
    assert len(r) == 18
    return r


def _run(alt_target, ref_aligned, pam="NRG", pam_begin=-3, pam_end=None):
    ns = {"iupac_code_set": IUPAC, "final_result": _row(alt_target),
          "ref_aligned": ref_aligned, "pam": pam,
          "pam_begin": pam_begin, "pam_end": pam_end}
    exec(compile(_block(), "<pam_creation>", "exec"), ns)
    return ns["final_result"]


class IndelPamCreation(unittest.TestCase):

    def test_fires_into_pam_gen_and_leaves_var_uniq_alone(self):
        # reference window GAT is not NRG (position 2: R={A,G} vs T) -> creation
        r = _run("ACGTACGTACGTACGTACGTTGG", "acgtacgtacgtacgtacgtGAT")
        self.assertEqual(r[PAM_GEN], "TGG")
        # the original bug wrote here instead; Var_uniq must be untouched
        self.assertEqual(r[VAR_UNIQ], "n")

    def test_noop_when_reference_pam_already_valid(self):
        # CAG IS a valid NRG (N=C, R=A, G=G) -> not a creation event
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtCAG")
        self.assertEqual(r[PAM_GEN], "n")
        self.assertEqual(r[VAR_UNIQ], "n")

    def test_second_real_window_from_live_run(self):
        # GCA -> AGG, observed in the instrumented complete-test run
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtGCA")
        self.assertEqual(r[PAM_GEN], "AGG")

    def test_dash_in_pam_window_is_skipped_not_raised(self):
        # a bulge dash inside the PAM window: '-' is not in the IUPAC table
        r = _run("ACGTACGTACGTACGTACGTT-G", "acgtacgtacgtacgtacgtG-T")
        self.assertEqual(r[PAM_GEN], "n")

    def test_short_window_is_skipped(self):
        r = _run("ACGTTGG", "AC")  # ref_aligned shorter than the PAM
        self.assertEqual(r[PAM_GEN], "n")

    def test_lowercase_reference_window_is_not_a_false_positive(self):
        # lowercase IUPAC keys map to lowercase sets, so without .upper() this
        # would wrongly score as a creation event
        r = _run("ACGTACGTACGTACGTACGTAGG", "acgtacgtacgtacgtacgtcag")
        self.assertEqual(r[PAM_GEN], "n")

    def test_five_prime_pam_window(self):
        # Cas12a-style leading PAM: TTTN at [0:4]
        r = _run("TTTAACGTACGTACGTACGTACGT", "GGGAacgtacgtacgtacgtacgt",
                 pam="TTTN", pam_begin=0, pam_end=4)
        self.assertEqual(r[PAM_GEN], "TTTA")

    def test_source_targets_index_10_not_11(self):
        """Static guard: the assignment must name index 10.

        The sed header in analisi_indels_NNN.sh lists a Chromosome_fake column at 10,
        which invites re-deriving this as 11. The rows the analyzer reads have no such
        column (18 fields; verified by instrumenting a live run), and this function's
        own [12]/[15]/[16]/[17] assignments confirm the SNP layout.
        """
        src = open(SRC).read()
        self.assertIn("final_result[10] = final_result[2][pam_begin:pam_end].upper()", src)
        self.assertNotIn("final_result[11] = final_result[2][pam_begin:pam_end].upper()", src)


if __name__ == "__main__":
    unittest.main()

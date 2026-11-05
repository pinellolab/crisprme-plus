#!/usr/bin/env python3
"""The genome-string IUPAC collapse must mirror the code table, and be the only one.

IUPAC ambiguity bases reach the indel post-analysis (``analisi_indels_NNN.py``) from
the SNP-overlaid fake-indel genome (``CRISPRME_INDEL_SNP``, default on), in the DNA
column of the engine's targets. Up to 2.8.0 that column was collapsed column by column
too, by ``_resolve_overlaid_iupac`` (the guide's base when it was an allele, else
``alleles[0]``), without any genotype: it reported one sample's indel with SNP alleles
only other samples carry, and read an R at NRG's last PAM base as A. Since 2.8.1 the
target is resolved per carrier haplotype (``indel_snp_cis.resolve_target``), so what
is left here is ``_IUPAC_FIRST_ALLELE_TAB``: a ``str.maketrans`` table applied once to
``genomeStr``, the PLAIN reference chromosome, where it is a no-op kept as a guard.

The table is a hand-written mirror of every ``iupac_code`` entry's first element, so
adding a code to the dict, or reordering one entry's alleles, would silently leave a
wrong base in a report column (reference alignment, CRISPR-Bulge context, motif/radar
counts) should a code ever reach ``genomeStr``. These tests derive the expected table
FROM ``iupac_code`` and compare, so the mirror cannot rot, and pin that no target
column is collapsed by first allele again. Hermetic: the module is parsed, never
imported -- it reads ``sys.argv`` at import time and would exit.
"""

import ast
import unittest
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "analisi_indels_NNN.py"


def _literal(name):
    """Return the module-level literal assigned to `name`, without importing."""
    tree = ast.parse(_SRC.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == name:
                    return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found as a module-level literal in {_SRC.name}")


def _maketrans_args(name):
    """Return the (from, to) string args of a `name = str.maketrans("..", "..")`."""
    tree = ast.parse(_SRC.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            continue
        call = node.value
        assert isinstance(call, ast.Call), f"{name} is not a call"
        args = [ast.literal_eval(a) for a in call.args]
        assert len(args) == 2, f"{name}: expected str.maketrans(from, to)"
        return args[0], args[1]
    raise AssertionError(f"{name} not found")


class TheTableMirrorsTheDict(unittest.TestCase):
    def setUp(self):
        self.code = _literal("iupac_code")
        self.frm, self.to = _maketrans_args("_IUPAC_FIRST_ALLELE_TAB")

    def test_table_is_well_formed(self):
        self.assertEqual(len(self.frm), len(self.to), "maketrans args differ in length")
        self.assertEqual(len(set(self.frm)), len(self.frm), "a source char repeats")

    def test_every_table_entry_matches_the_dicts_first_allele(self):
        """iupac_code lists alleles in UPPERCASE even under a lowercase key, while
        the table must keep a lowercase code lowercase. Compare against the dict's
        first allele put back into the code's own case."""
        for src, dst in zip(self.frm, self.to):
            self.assertIn(src, self.code, f"{src!r} is in the table but not iupac_code")
            first = self.code[src][0]
            expected = first.upper() if src.isupper() else first.lower()
            self.assertEqual(
                dst, expected,
                f"table maps {src!r}->{dst!r} but iupac_code's first allele, "
                f"in the code's own case, is {expected!r}",
            )

    def test_every_ambiguity_code_is_in_the_table(self):
        """A code in the dict but missing from the table survives the fast path as
        IUPAC and then reaches an ACGT-only consumer (CRISPR-Bulge's agct2numerals
        does a bare dict lookup -> KeyError, aborting scoring)."""
        for code, alleles in self.code.items():
            if len(alleles) == 1:
                continue        # plain A/C/G/T, nothing to collapse
            if code in ("N", "n"):
                continue        # deliberately absent: preserved for the N-guard
            self.assertIn(
                code, self.frm,
                f"{code!r} ({alleles}) is an ambiguity code with no table entry",
            )

    def test_N_is_deliberately_absent(self):
        """'N' must survive the collapse so the existing N-guard still nulls the
        window. Putting it in the table would silently turn unknown sequence into
        a concrete base."""
        self.assertNotIn("N", self.frm)
        self.assertNotIn("n", self.frm)
        self.assertIn("N", self.code)

    def test_case_is_preserved(self):
        """Lowercase input must collapse to lowercase: the genome string is
        upper-cased before this runs, but the per-sequence function preserves case
        and the two must agree."""
        for src, dst in zip(self.frm, self.to):
            if src.islower():
                self.assertTrue(dst.islower(), f"{src!r}->{dst!r} changed case")
            else:
                self.assertTrue(dst.isupper(), f"{src!r}->{dst!r} changed case")

    def test_both_cases_of_each_code_are_covered(self):
        for src in self.frm:
            self.assertIn(
                src.swapcase(), self.frm,
                f"{src!r} is in the table but {src.swapcase()!r} is not",
            )

    def test_the_table_comment_claim_is_true(self):
        """The table's comment says it maps each code to its first-listed allele."""
        src = _SRC.read_text()
        self.assertIn("first-listed-allele", src)

    def test_no_target_column_is_collapsed_by_first_allele(self):
        """A target's DNA column is resolved per carrier haplotype, never by a
        per-column pick that ignores the genotypes (the 2.5.0-2.8.0 resolver)."""
        src = _SRC.read_text()
        self.assertNotIn("_resolve_overlaid_iupac", src)
        self.assertNotIn("pick = alleles[0]", src)
        self.assertIn("_isc.resolve_target(", src)
        self.assertNotIn("final_result[2].translate(", src)

    def test_the_collapse_is_a_no_op_on_plain_acgt(self):
        for base in "ACGTacgt":
            self.assertNotIn(base, self.frm, f"{base!r} must pass through untouched")


class TheOverlayIsDefaultOn(unittest.TestCase):
    """If the SNP+indel overlay is default-on, the collapse is on the DEFAULT path,
    not an opt-in corner -- which is why the mirror above has to be pinned."""

    def test_indel_snp_defaults_to_enabled(self):
        for f in ("analisi_indels_NNN.py", "overlay_indel_snps.py"):
            src = (_SRC.parent / f).read_text()
            self.assertIn(
                'os.environ.get("CRISPRME_INDEL_SNP", "1")', src,
                f"{f}: the overlay default is no longer on",
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)

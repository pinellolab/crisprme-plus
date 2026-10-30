#!/usr/bin/env python3
"""The two IUPAC-collapse implementations must not drift apart.

``analisi_indels_NNN.py`` collapses IUPAC ambiguity bases -- which reach the indel
post-analysis only from the SNP-overlaid fake-indel genome (``CRISPRME_INDEL_SNP``,
default on) -- in two places, by two different mechanisms:

* ``_resolve_overlaid_iupac(seq, guide_aligned)`` walks ``iupac_code`` and picks the
  guide-matched allele when it can, else ``alleles[0]``;
* ``_IUPAC_FIRST_ALLELE_TAB``, a ``str.maketrans`` table applied to the whole
  genome string in one pass, which can only do the ``alleles[0]`` half.

The table is a hand-written mirror of every ``iupac_code`` entry's first element.
Nothing tied the two together, so adding a code to the dict, or reordering one
entry's alleles, would silently leave the fast path disagreeing with the slow one
-- and the disagreement would show up as a wrong base in a report column
(reference alignment, CRISPR-Bulge context, motif/radar counts), not as an error.

These tests derive the expected table FROM ``iupac_code`` and compare, so the
mirror cannot rot. Hermetic: the module is parsed, never imported -- it reads
``sys.argv`` at import time and would exit.
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
        """iupac_code lists alleles in UPPERCASE even under a lowercase key, and
        _resolve_overlaid_iupac re-cases its pick (`pick if c.isupper() else
        pick.lower()`). The table has to bake that in, so compare against the
        dict's first allele put back into the code's own case."""
        for src, dst in zip(self.frm, self.to):
            self.assertIn(src, self.code, f"{src!r} is in the table but not iupac_code")
            first = self.code[src][0]
            expected = first.upper() if src.isupper() else first.lower()
            self.assertEqual(
                dst, expected,
                f"table maps {src!r}->{dst!r} but iupac_code's first allele, "
                f"re-cased as _resolve_overlaid_iupac does, is {expected!r}",
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

    def test_the_docstring_claim_is_true(self):
        """_resolve_overlaid_iupac's docstring states the fallback is alleles[0];
        the table exists only to do that same thing faster."""
        src = _SRC.read_text()
        self.assertIn("pick = alleles[0]", src)
        self.assertIn("first-listed-allele", src)

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

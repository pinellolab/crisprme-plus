#!/usr/bin/env python
"""``scripts/check_alt_pam.py``: the release gate must fail a run whose ALT row shows a PAM
the search never matched, or that reports a target holding an N, and pass one that shows
neither.

Up to 2.8.0 the indel post-analysis read an R at NRG's last PAM base as A (a genome-wide
HGDP + 1000 Genomes run reported 1,088 such rows), and on a classic dictionary install
left the IUPAC codes themselves in the ALT column. Both must fail the gate; a valid NRG
or 5' TTTV PAM, a REF row, and an NA ALT must pass. An N of the genome is an assembly gap:
an N (or an n) in the REF or ALT target of any block, a REF row's too, must fail it.
Hermetic, stdlib only.

Run: python -m unittest test_check_alt_pam -v
"""

import importlib.util
import io
import os
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout

_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts",
                       "check_alt_pam.py")
_spec = importlib.util.spec_from_file_location("check_alt_pam", _SCRIPT)
gate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gate)

BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _table(rows):
    """An integrated-results-shaped table: rows of (spacer, [(origin, alt[, ref]) x 3]);
    the REF target is the ALT one, with no N, unless given."""
    head = ["Spacer+PAM", "Chromosome"]
    for b in BLOCKS:
        head += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                 "Aligned_spacer+PAM_(%s)" % b, "Aligned_protospacer+PAM_REF_(%s)" % b,
                 "Aligned_protospacer+PAM_ALT_(%s)" % b, "REF/ALT_origin_(%s)" % b]
    lines = ["\t".join(head)]
    for spacer, blocks in rows:
        cells = [spacer, "chr1"]
        for block in blocks:
            origin, alt = block[:2]
            ref = block[2] if len(block) > 2 else (alt if alt != "NA" else spacer[:-3] + "TGG")
            cells += ["100", "+", spacer, ref, alt, origin]
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


class AltPamGate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="check_alt_pam_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _run(self, rows, motif=None, params=None):
        path = os.path.join(self.tmp, "x_integrated_results.tsv")
        with open(path, "w") as fh:
            fh.write(_table(rows))
        if params:
            with open(os.path.join(self.tmp, ".Params.txt"), "w") as fh:
                fh.write("Genome_ref\thg38\nPam\t%s\nMismatches\t4\n" % params)
        argv = ["check_alt_pam.py"] + (["--motif", motif] if motif else []) + [path]
        out = io.StringIO()
        with redirect_stdout(out):
            rc = gate.main(argv)
        return rc, out.getvalue()

    G = "ACTGAAATCTGTAAGCAGGC"

    def test_a_misread_r_fails(self):
        rc, out = self._run([(self.G + "NNN", [("alt", self.G + "TGA"), ("alt", self.G + "TGG"),
                                               ("ref", self.G + "TGA")])], motif="NRG")
        self.assertEqual(rc, 1, out)
        self.assertIn("1 with a PAM that violates NRG", out)
        self.assertIn("highest_CFD chr1:100(+)", out)

    def test_a_left_iupac_code_fails(self):
        rc, out = self._run([(self.G + "NNN", [("alt", self.G + "TGR")] * 3)], motif="NRG")
        self.assertEqual(rc, 1, out)
        self.assertIn("3 with a PAM", out)

    def test_valid_rows_pass(self):
        rows = [(self.G + "NNN", [("alt", self.G + "AGG"), ("alt", "A-" + self.G[1:] + "CGG"),
                                  ("alt", "NA")]),
                (self.G + "NNN", [("ref", self.G + "TTT"), ("alt", self.G.lower() + "TGG"),
                                  ("alt", self.G + "GAG")])]
        rc, out = self._run(rows, motif="NRG")
        self.assertEqual(rc, 0, out)
        self.assertIn("4 ALT blocks, 0 with", out)

    def test_the_motif_from_params_and_a_5prime_pam(self):
        g = "NNNN" + "CCTTGTCAAGGCTATTGGTC"
        ok = "TTTA" + "CCTTGTCAAGGCTATTGGTC"
        bad = "TTTT" + "CCTTGTCAAGGCTATTGGTC"  # V allows no T
        rc, out = self._run([(g, [("alt", ok), ("alt", bad), ("alt", ok)])], params="TTTV")
        self.assertEqual(rc, 1, out)
        self.assertIn("1 with a PAM that violates TTTV", out)

    def test_valid_rows_hold_no_n(self):
        rows = [(self.G + "NNN", [("alt", self.G + "AGG"), ("ref", "NA", self.G + "TGG"),
                                  ("alt", self.G.lower() + "TGG")])]
        rc, out = self._run(rows, motif="NRG")
        self.assertEqual(rc, 0, out)
        self.assertIn("5 reported targets (REF and ALT, every block), 0 with an N", out)

    def test_an_n_in_a_reported_target_fails(self):
        gap = self.G[:10] + "N" + self.G[11:]
        for name, blocks, where in (
                ("a REF row's target", [("ref", "NA", gap + "TGG"), ("ref", "NA", self.G + "TGG"),
                                        ("ref", "NA", self.G + "TGG")], "highest_CFD REF"),
                ("an ALT target, as an n", [("alt", self.G + "TGG"),
                                            ("alt", self.G[:-1] + "n" + "TGG", self.G + "TGG"),
                                            ("alt", self.G + "TGG")], "fewest_mm+b ALT"),
                ("the REF of an ALT row, at its PAM", [("alt", self.G + "TGG"),
                                                       ("alt", self.G + "TGG"),
                                                       ("alt", self.G + "AGG", self.G + "ANG")],
                 "highest_CRISPR_BULGE REF")):
            with self.subTest(name):
                rc, out = self._run([(self.G + "NNN", blocks)], motif="NRG")
                self.assertEqual(rc, 1, out)
                self.assertIn("1 with an N", out)
                self.assertIn("0 with a PAM that violates NRG", out)
                self.assertIn("    %s chr1:100(+)" % where, out)

    def test_no_motif_is_a_usage_error(self):
        rc, out = self._run([(self.G + "NNN", [("alt", self.G + "TGG")] * 3)])
        self.assertEqual(rc, 2, out)

    def test_the_docstring_names_its_blind_spot(self):
        self.assertIn("NRR", gate.__doc__)
        self.assertIn("BLIND SPOT", gate.__doc__)


if __name__ == "__main__":
    unittest.main(verbosity=2)

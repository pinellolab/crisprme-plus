#!/usr/bin/env python
"""``scripts/check_alt_pam.py``: the release gate must fail a run whose ALT row shows a PAM
the search never matched, or that reports a target holding an N, and pass one that shows
neither.

Up to 2.8.0 the indel post-analysis read an R at NRG's last PAM base as A (a genome-wide
HGDP + 1000 Genomes run reported 1,088 such rows), and on a classic dictionary install
left the IUPAC codes themselves in the ALT column. Both must fail the gate; a valid NRG
or 5' TTTV PAM, a REF row, and an NA ALT must pass. An N of the genome is an assembly gap:
an N (or an n) in the REF or ALT target of any block, a REF row's too, must fail it, but
for the REF text of an indel row where the reference genome itself holds an N (2.10.0: the
reference at the row's coordinates can run into a gap next to the indel; HPRC's chr22,
chr22:11978143 (-) and 11980225 (+)). Hermetic, stdlib only.

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

    # -- an indel row's REF text over the reference's own N (2.10.0) ---------------------

    def _gap_run(self, cells, genome_arg=True, params_genome=False):
        """A run's table of indel-row blocks (start, strand, REF, ALT, variants[, origin:
        "alt" unless given]) over chr1, whose reference holds N at 100-129 and bases
        elsewhere; --genome given, or found as the run's Genomes/<Genome_ref> beside its
        Results folder."""
        data = os.path.join(self.tmp, "DATA")
        gdir = os.path.join(data, "Genomes", "hgT")
        os.makedirs(gdir)
        with open(os.path.join(gdir, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 25 + "N" * 30 + "ACGT" * 25 + "\n")
        res = os.path.join(data, "Results", "job")
        os.makedirs(res)
        if params_genome:
            with open(os.path.join(res, ".Params.txt"), "w") as fh:
                fh.write("Genome_ref\thgT\nPam\tNRG\n")
        head = ["Spacer+PAM", "Chromosome"]
        for b in BLOCKS:
            head += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                     "Aligned_spacer+PAM_(%s)" % b, "Aligned_protospacer+PAM_REF_(%s)" % b,
                     "Aligned_protospacer+PAM_ALT_(%s)" % b, "REF/ALT_origin_(%s)" % b,
                     "Variant_info_genome_(%s)" % b]
        lines = ["\t".join(head)]
        g = self.G + "NNN"
        for start, strand, ref, alt, var, *origin in cells:
            row = [g, "chr1"]
            for _b in BLOCKS:
                row += [str(start), strand, g, ref, alt, (origin or ["alt"])[0], var]
            lines.append("\t".join(row))
        path = os.path.join(res, "x_integrated_results.tsv")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        argv = ["check_alt_pam.py", "--motif", "NRG"] + (
            ["--genome", gdir] if genome_arg else []) + [path]
        out = io.StringIO()
        with redirect_stdout(out):
            rc = gate.main(argv)
        return rc, out.getvalue()

    INDEL = "chr1_99_A_ACTTAG"
    ALT = "ACTGAAATCTGTAAGCAGGCTGG"

    def test_an_indel_rows_ref_text_over_the_references_n_passes(self):
        """The REF text of an indel row is the reference at its start, as many bases as its
        target holds (gapped as it, its mismatches in lower case): where the reference holds
        N there, so may it: + strand from 105, a gap and lower-case n's; - strand from 100,
        reverse-complemented, ending on a reference base past the gap; the genome given, or
        the run's own (Genomes/<Genome_ref> beside its Results folder)."""
        plus = "nnnn-nnnnnnnnnnnnNnnnNNN"           # 23 bases: 105-127, all N
        minus = "nnnnnnnnnnnnnnnnnnnnNNN"           # 100-122 reverse-complemented: all N
        edge = "TnnnnnnnnnnnnnnnnnnnNNN"            # - strand from 108: 130 is A (T on -)
        for genome_arg, params in ((True, False), (False, True)):
            with self.subTest(genome_arg=genome_arg):
                shutil.rmtree(os.path.join(self.tmp, "DATA"), ignore_errors=True)
                rc, out = self._gap_run([(105, "+", plus, self.ALT[:4] + "-" + self.ALT[4:],
                                          self.INDEL), (100, "-", minus, self.ALT, self.INDEL),
                                         (108, "-", edge, self.ALT, self.INDEL)],
                                        genome_arg, params)
                self.assertEqual(rc, 0, out)
                self.assertIn("0 with an N; 9 indel rows' REF texts hold an N only where "
                              "the reference genome does", out)

    def test_an_n_where_the_reference_holds_a_base_fails(self):
        """An N of an indel row's REF text at a position where the reference holds a base
        fails: - strand from 109, the text's first base is 131 (C on +, G shown), an N here."""
        rc, out = self._gap_run([(100, "-", "N" * 23, self.ALT, self.INDEL),
                                 (109, "-", "N" * 23, self.ALT, self.INDEL)])
        self.assertEqual(rc, 1, out)
        self.assertIn("3 with an N; 3 indel rows' REF texts hold an N only where the "
                      "reference genome does", out)
        self.assertIn("    highest_CFD REF chr1:109(-) " + "N" * 23, out)

    def test_an_n_in_an_alt_text_fails_even_over_the_references_n(self):
        rc, out = self._gap_run([(100, "+", "N" * 23, "N" * 23, self.INDEL)])
        self.assertEqual(rc, 1, out)
        self.assertIn("3 with an N; 3 indel rows' REF texts hold an N only where the "
                      "reference genome does", out)
        self.assertIn("    fewest_mm+b ALT chr1:100(+) " + "N" * 23, out)

    def test_a_snp_or_mnv_rows_ref_text_over_the_references_n_fails(self):
        """Only an indel row's REF text may hold the reference's N: a SNP row's or an MNV
        row's target is the reference at its own coordinates, which holds no gap."""
        rc, out = self._gap_run([(100, "+", "N" * 23, self.ALT, self.INDEL),
                                 (100, "+", "N" * 23, self.ALT, "chr1_110_A_G"),
                                 (100, "+", "N" * 23, self.ALT, "chr1_110_AC_GT")])
        self.assertEqual(rc, 1, out)
        self.assertIn("6 with an N; 3 indel rows' REF texts", out)

    def test_a_ref_origin_blocks_ref_text_over_the_references_n_fails(self):
        """The exception is an indel row's: an ALT-origin block, whose REF text is the
        reference at its start and not its target. A REF-origin block's REF text is the
        reported target itself, so an N there fails whatever its Variant_info_genome names
        (no table holds such a block, whose Variant_info_genome is NA; pamc-review-r1's G2)."""
        rc, out = self._gap_run([(100, "+", "N" * 23, self.ALT, self.INDEL),
                                 (100, "+", "N" * 23, "NA", self.INDEL, "ref")])
        self.assertEqual(rc, 1, out)
        self.assertIn("3 with an N; 3 indel rows' REF texts", out)
        self.assertIn("    highest_CFD REF chr1:100(+) " + "N" * 23, out)

    def test_an_indel_rows_ref_n_without_a_genome_fails(self):
        rc, out = self._gap_run([(100, "+", "N" * 23, self.ALT, self.INDEL)],
                                genome_arg=False)
        self.assertEqual(rc, 1, out)
        self.assertIn("3 with an N", out)
        self.assertIn("no reference genome to read an indel row's REF text against", out)

    def test_no_motif_is_a_usage_error(self):
        rc, out = self._run([(self.G + "NNN", [("alt", self.G + "TGG")] * 3)])
        self.assertEqual(rc, 2, out)

    def test_the_docstring_names_its_blind_spot(self):
        self.assertIn("NRR", gate.__doc__)
        self.assertIn("BLIND SPOT", gate.__doc__)


if __name__ == "__main__":
    unittest.main(verbosity=2)

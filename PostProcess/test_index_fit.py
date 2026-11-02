#!/usr/bin/env python3
"""Which guide lengths an index answers (index_fit.py), and that the search script, which
has its own copy of the rule in bash, says the same in every measured cell."""

import os
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path

import index_fit

_SCRIPT = Path(__file__).resolve().parent / "submit_job_automated_new_multiple_vcfs.sh"


def _functions(*names):
    """Those functions of the shipped script, verbatim, as top-level functions."""
    text = _SCRIPT.read_text()
    out = []
    for name in names:
        start = text.index(f"\t{name}() {{")
        end = text.index("\n\t}\n", start) + len("\n\t}\n")
        out.append("\n".join(line[1:] if line.startswith("\t") else line
                             for line in text[start:end].splitlines()))
    return "\n".join(out)


def _bin(folder, guide):
    with open(os.path.join(folder, "NGG_chr1_1.bin"), "wb") as fh:
        fh.write(struct.pack("=ii", 7, guide) + b"\0" * 8)


class TheHeaderOfAnIndex(unittest.TestCase):
    def test_the_second_int32_of_a_bin_is_the_built_guide_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            for guide in (18, 20, 21, 23):
                _bin(tmp, guide)
                self.assertEqual(index_fit.built_guide_len(tmp), guide)

    def test_no_bin_or_one_too_short_tells_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(index_fit.built_guide_len(tmp))
            Path(tmp, "x.bin").write_bytes(b"")
            self.assertIsNone(index_fit.built_guide_len(tmp))
            Path(tmp, "x.bin").write_bytes(struct.pack("=i", 7))
            self.assertIsNone(index_fit.built_guide_len(tmp))
        self.assertIsNone(index_fit.built_guide_len("/no/such/folder"))

    def test_a_sidecar_is_not_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, ".pam_build").write_text("23bp-NGG-X 26 3\n")
            _bin(tmp, 20)
            self.assertEqual(index_fit.built_guide_len(tmp), 20)


class ThePamFile(unittest.TestCase):
    def test_its_spacer(self):
        with tempfile.TemporaryDirectory() as tmp:
            for line, spacer in (("NNNNNNNNNNNNNNNNNNNNNGG 3\n", 20),
                                 ("TTTVNNNNNNNNNNNNNNNNNNNNNNN -4\n", 23),
                                 ("NGG 3\n", 0)):
                p = os.path.join(tmp, "p.txt")
                Path(p).write_text(line)
                self.assertEqual(index_fit.pam_file_spacer(p), spacer)
            Path(p).write_text("garbage\n")
            self.assertIsNone(index_fit.pam_file_spacer(p))
        self.assertIsNone(index_fit.pam_file_spacer("/no/such/file"))


class TheRule(unittest.TestCase):
    def test_the_bulge_cap_is_the_most_of_each_type_that_fits(self):
        for built in (18, 20, 23):
            for n in (1, 2, 3, 4):
                for guide in range(18, 26):
                    cap = index_fit.bulge_cap(built, n, guide)
                    self.assertLessEqual(cap, n - 1)
                    if cap:
                        self.assertTrue(index_fit.fits(built, n, guide, cap, cap))
                    if cap < n - 1:
                        self.assertFalse(index_fit.fits(built, n, guide, cap + 1, cap + 1))

    def test_what_is_not_known_fits(self):
        self.assertTrue(index_fit.fits(None, 2, 23, 1, 1))
        self.assertTrue(index_fit.fits(20, 2, None, 1, 1))
        self.assertTrue(index_fit.fits(20, 2, 0, 1, 1))
        self.assertEqual(index_fit.bulge_cap(None, 3, 23), 2)
        self.assertEqual(index_fit.bulge_cap(20, 3, None), 2)


class TheScriptSaysTheSame(unittest.TestCase):
    """The bash copy of the rule, on an index folder, over a grid wider than the cells
    measured."""

    HELPERS = _functions("_index_guide_len", "_index_fits_guide")

    def script_says(self, folder, n, guide, bdna, brna):
        script = (f"pos=3\nfullseqpam={'N' * (guide + 3)}\nbDNA={bdna}\nbRNA={brna}\n"
                  f"{self.HELPERS}\n_index_fits_guide \"$1\" {n}\n")
        out = subprocess.run(["bash", "-c", script, "bash", folder], capture_output=True, text=True)
        self.assertEqual(out.stderr, "")
        return out.returncode == 0

    def test_in_every_cell(self):
        with tempfile.TemporaryDirectory() as tmp:
            for built in (20, 22):
                folder = os.path.join(tmp, f"NGG_x_hg38+ds{built}")
                os.makedirs(folder)
                _bin(folder, built)
                for n in (2, 3, 4):
                    for guide in (20, 21, 22, 23, 24):
                        for bdna, brna in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2),
                                           (2, 1), (1, 2), (2, 2), (3, 0), (0, 3)):
                            with self.subTest(built=built, n=n, guide=guide, b=(bdna, brna)):
                                self.assertEqual(self.script_says(folder, n, guide, bdna, brna),
                                                 index_fit.fits(built, n, guide, bdna, brna))


if __name__ == "__main__":
    unittest.main(verbosity=2)

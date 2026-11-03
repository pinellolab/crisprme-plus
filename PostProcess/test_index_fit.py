#!/usr/bin/env python3
"""Which guide lengths an index answers (index_fit.py), and that the search script, which
has its own copy of the rule in bash, says the same in every measured cell. An index a
newer CRISPRitz wrote (index format 2, CRISPRitz 2.9.2 and later: the magic "CZI" 0xFE,
then the format version) holds no guide length this release reads: neither reads its
version as one (a format-2 index read as built for a 2-nt guide, and refused as too short
for any search), and the search script stops on it with one ERROR."""

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


def _bin2(folder, guide=20, version=2, name="NGG_chr1_1.bin"):
    """A .bin in index format 2 (CRISPRitz's docs/INDEX_FORMAT.md, section 4): the magic,
    the version, the header length, the flags, the guide length, the PAM length, the max
    bulges and the PAM, then the body (leaf count and guide length) and the magic again."""
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, name), "wb") as fh:
        fh.write(b"CZI\xfe" + struct.pack("=6i", version, 31, 0, guide, 3, 2) + b"NGG"
                 + struct.pack("=ii", 7, guide) + b"\0" * 8 + b"CZI\xfe")


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

    def test_a_newer_index_format_is_said_not_misread(self):
        """bytes 4..7 of a format-2 .bin are its version, 2: never a guide length."""
        with tempfile.TemporaryDirectory() as tmp:
            index = os.path.join(tmp, "NGG_2_hg38")
            _bin2(index)
            self.assertIsNone(index_fit.built_guide_len(index))
            self.assertTrue(index_fit.fits(index_fit.built_guide_len(index), 2, 20, 1, 1))
            self.assertEqual(index_fit.bin_format(os.path.join(index, "NGG_chr1_1.bin")), 2)
            self.assertEqual(index_fit.newer_bin(index),
                             (os.path.join(index, "NGG_chr1_1.bin"), 2))
            msg = index_fit.newer_format_error(index)
            self.assertEqual(msg, (
                f"{index} was built by a newer CRISPRitz (NGG_chr1_1.bin is in index format 2, "
                f"which CRISPRitz 2.9.2 and later write): this CRISPRme reads index format 1 "
                f"only (CRISPRitz up to 2.9.1), and needs a newer release to search it. Use a "
                f"CRISPRme that reads index format 2, or install CRISPRitz 2.9.1 and build the "
                f"index again, with {index} moved out of the way first"))
            # an INDELS index holds one folder a chromosome
            indels = os.path.join(tmp, "NGG_2_hg38+ds_INDELS")
            os.makedirs(os.path.join(indels, "NGG_2_fakechr1"))
            _bin(os.path.join(indels, "NGG_2_fakechr1"), 20)
            _bin2(os.path.join(indels, "NGG_2_fakechr2"), name="NGG_fakechr2_1.bin")
            self.assertEqual(index_fit.newer_bin(indels)[1], 2)
            # a format-1 index, or none, is not refused
            plain = os.path.join(tmp, "NGG_3_hg38")
            os.makedirs(plain)
            _bin(plain, 20)
            self.assertEqual(index_fit.bin_format(os.path.join(plain, "NGG_chr1_1.bin")), 1)
            self.assertIsNone(index_fit.newer_format_error(plain))
            self.assertIsNone(index_fit.newer_format_error(os.path.join(tmp, "none")))


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

    HELPERS = _functions("_index_guide_len", "_index_fits_guide", "_index_format",
                         "_refuse_newer_index")

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


    def run_helpers(self, folder, call):
        script = f"pos=3\nfullseqpam={'N' * 23}\nbDNA=1\nbRNA=1\n{self.HELPERS}\n{call}\n"
        return subprocess.run(["bash", "-c", script, "bash", folder], capture_output=True,
                              text=True)

    def test_a_newer_index_format_is_not_read_as_a_guide_length(self):
        """Up to this fix the script read a format-2 index as built for a 2-nt guide: it
        stopped a search with 'NGG_2_hg38 was built for a 2-nt guide, and this search has a
        20-nt guide with 1 DNA and 1 RNA bulge(s)...'."""
        with tempfile.TemporaryDirectory() as tmp:
            index = os.path.join(tmp, "NGG_2_hg38")
            _bin2(index)
            out = self.run_helpers(index, '_index_guide_len "$1"')
            self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "", ""))
            out = self.run_helpers(index, '_index_fits_guide "$1" 2')
            self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "", ""))
            out = self.run_helpers(index, '_refuse_newer_index "$1"; echo "went on"')
            self.assertEqual(out.returncode, 1)
            self.assertEqual(out.stdout, "")
            self.assertEqual(out.stderr, index_fit.newer_format_error(index).replace(
                f"{index} was", f"ERROR: {index} was", 1) + "\n")

    def test_an_indels_index_with_a_newer_chromosome_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            indels = os.path.join(tmp, "NGG_2_hg38+ds_INDELS")
            os.makedirs(os.path.join(indels, "NGG_2_fakechr1"))
            _bin(os.path.join(indels, "NGG_2_fakechr1"), 20)
            _bin2(os.path.join(indels, "NGG_2_fakechr2"), name="NGG_fakechr2_1.bin")
            out = self.run_helpers(indels, '_refuse_newer_index "$1"')
            self.assertEqual(out.returncode, 1)
            self.assertIn("NGG_2_fakechr2/NGG_fakechr2_1.bin is in index format 2", out.stderr)
            self.assertEqual(out.stderr.count("ERROR: "), 1)

    def test_a_format_1_index_goes_on(self):
        with tempfile.TemporaryDirectory() as tmp:
            index = os.path.join(tmp, "NGG_2_hg38")
            os.makedirs(index)
            _bin(index, 20)
            out = self.run_helpers(index, '_refuse_newer_index "$1"; echo "went on"')
            self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "went on\n", ""))


if __name__ == "__main__":
    unittest.main(verbosity=2)

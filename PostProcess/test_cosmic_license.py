"""Tests for the COSMIC licence gate (PostProcess/cosmic_license.py)."""

import gzip
import os
import tempfile
import unittest

import cosmic_license as cl


_BED_ROWS = [
    "chr1\t100\t200\tENCODE_cCRE_1\n",
    "chr1\t150\t250\tTP53_COSMIC\n",
    "chr1\t300\t400\tDHS_peak_3\n",
    "chr2\t500\t600\tBRCA1_COSMIC\n",
    "chr2\t700\t800\tgene_X_gencode\n",
]


class TestCosmicFlag(unittest.TestCase):
    def test_default_disabled_and_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            # missing flag file -> disabled (licence-safe default)
            self.assertFalse(cl.cosmic_enabled(d))
            cl.set_cosmic_license(d, True, "attested via test")
            self.assertTrue(cl.cosmic_enabled(d))
            self.assertEqual(cl.get_cosmic_license(d).get("attestation"), "attested via test")
            cl.set_cosmic_license(d, False)
            self.assertFalse(cl.cosmic_enabled(d))

    def test_corrupt_flag_file_is_disabled(self):
        with tempfile.TemporaryDirectory() as d:
            with open(cl.license_path(d), "w") as fh:
                fh.write("not json{")
            self.assertFalse(cl.cosmic_enabled(d))


class TestCosmicRowDetection(unittest.TestCase):
    def test_is_cosmic_row(self):
        self.assertTrue(cl.is_cosmic_row("chr1\t150\t250\tTP53_COSMIC\n"))
        self.assertTrue(cl.is_cosmic_row("chr2\t500\t600\tBRCA1_COSMIC"))
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tENCODE_cCRE_1\n"))
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tgene_gencode\n"))
        # the tag only counts in the name column, not coordinates/other cols
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tDHS\tnote_COSMIC_ish\n"))


class TestStripCosmic(unittest.TestCase):
    def _write(self, path, rows, gz=False):
        op = gzip.open if gz else open
        with op(path, "wt") as fh:
            fh.writelines(rows)

    def test_strip_plain_bed(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed")
            out = os.path.join(d, "ann.nocosmic.bed")
            self._write(src, _BED_ROWS)
            self.assertTrue(cl.bed_has_cosmic(src))
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual(dropped, 2)   # TP53_COSMIC + BRCA1_COSMIC
            self.assertEqual(kept, 3)
            with open(out) as fh:
                body = fh.read()
            self.assertNotIn("_COSMIC", body)
            self.assertIn("ENCODE_cCRE_1", body)
            self.assertIn("gene_X_gencode", body)
            self.assertFalse(cl.bed_has_cosmic(out))

    def test_strip_gzip_bed(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed.gz")
            out = os.path.join(d, "ann.nocosmic.bed.gz")
            self._write(src, _BED_ROWS, gz=True)
            self.assertTrue(cl.bed_has_cosmic(src))
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual((kept, dropped), (3, 2))
            with gzip.open(out, "rt") as fh:
                body = fh.read()
            self.assertNotIn("_COSMIC", body)

    def test_strip_keeps_comments(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed")
            out = os.path.join(d, "out.bed")
            self._write(src, ["#track cosmic annotations\n"] + _BED_ROWS)
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual(dropped, 2)
            with open(out) as fh:
                self.assertIn("#track cosmic annotations", fh.read())


if __name__ == "__main__":
    unittest.main()

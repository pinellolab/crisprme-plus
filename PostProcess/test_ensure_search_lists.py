"""Tests for _ensure_search_lists: makes a downloaded variant index searchable
from the CLI (list_vcf.txt / list_samplesID.txt) out of the box, matching what
the web form builds per-search.
"""

import os
import tempfile
import unittest

import crisprme_hf


def _mk_install(tmp, sid=True, vcf_name="hg38_1000G_HGDP"):
    os.makedirs(os.path.join(tmp, "samplesIDs"), exist_ok=True)
    if sid:
        open(
            os.path.join(tmp, "samplesIDs", f"{vcf_name}.samplesID.txt"), "w"
        ).close()


def _read(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


class TestEnsureSearchLists(unittest.TestCase):
    def test_creates_both_lists(self):
        with tempfile.TemporaryDirectory() as d:
            _mk_install(d)
            crisprme_hf._ensure_search_lists(d, "hg38_1000G_HGDP")
            self.assertEqual(_read(os.path.join(d, "list_vcf.txt")), ["hg38_1000G_HGDP"])
            self.assertEqual(
                _read(os.path.join(d, "list_samplesID.txt")),
                ["hg38_1000G_HGDP.samplesID.txt"],
            )

    def test_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            _mk_install(d)
            for _ in range(3):
                crisprme_hf._ensure_search_lists(d, "hg38_1000G_HGDP")
            self.assertEqual(_read(os.path.join(d, "list_vcf.txt")), ["hg38_1000G_HGDP"])
            self.assertEqual(
                _read(os.path.join(d, "list_samplesID.txt")),
                ["hg38_1000G_HGDP.samplesID.txt"],
            )

    def test_second_dataset_appends(self):
        with tempfile.TemporaryDirectory() as d:
            _mk_install(d, vcf_name="hg38_1000G_HGDP")
            _mk_install(d, vcf_name="hg38_gnomAD")
            crisprme_hf._ensure_search_lists(d, "hg38_1000G_HGDP")
            crisprme_hf._ensure_search_lists(d, "hg38_gnomAD")
            self.assertEqual(
                _read(os.path.join(d, "list_vcf.txt")),
                ["hg38_1000G_HGDP", "hg38_gnomAD"],
            )
            self.assertEqual(
                _read(os.path.join(d, "list_samplesID.txt")),
                ["hg38_1000G_HGDP.samplesID.txt", "hg38_gnomAD.samplesID.txt"],
            )

    def test_missing_samplesid_still_writes_vcf_list(self):
        with tempfile.TemporaryDirectory() as d:
            _mk_install(d, sid=False)
            crisprme_hf._ensure_search_lists(d, "hg38_1000G_HGDP")
            self.assertEqual(_read(os.path.join(d, "list_vcf.txt")), ["hg38_1000G_HGDP"])
            self.assertFalse(os.path.exists(os.path.join(d, "list_samplesID.txt")))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for gzip-compressed variant dictionaries.

Batteries variant indexes now ship the per-sample dicts GZIPPED
(``my_dict_<chrom>.json.gz``, ``log<chrom>.txt.gz``) so they stay ~3.5x smaller on
disk (~40-50GB not ~152GB for 1000G+HGDP) and are read on the fly. These tests guard:
  (a) chromosome discovery from gzipped indel logs (the ``.txt`` filter/slice must
      strip an optional ``.gz``), and
  (b) that publish/download route gzipped dict files exactly like the plain ones
      (the packaging is extension-agnostic).
Lightweight: no numpy/pandas, no real HuggingFace (the snapshot fetch is monkeypatched).
"""
import gzip
import os
import shutil
import tarfile
import tempfile
import unittest

import crisprme_hf as hf


class TestBundleGzippedDicts(unittest.TestCase):
    """publish (_make_index_tarball) + download must route gzipped dict files exactly
    like plain ones (dir-level, extension-agnostic)."""

    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="gzbundle_")
        self.name = "NGG_2_hg38+hg38_1000G"
        vcf = self.name.partition("+")[2]
        idx = os.path.join(self.root, "genome_library", self.name)
        ind = idx + "_INDELS"
        d_snp = os.path.join(self.root, "Dictionaries", f"dictionaries_{vcf}")
        d_ind = os.path.join(self.root, "Dictionaries", f"log_indels_{vcf}")
        for d in (idx, ind, d_snp, d_ind):
            os.makedirs(d)
        with open(os.path.join(idx, "TSTgenome.NGG.bin"), "w") as fh:
            fh.write("x")
        with open(os.path.join(ind, "NGG_2_fakechr22"), "w") as fh:
            fh.write("x")
        with gzip.open(os.path.join(d_snp, "my_dict_chr22.json.gz"), "wt") as fh:
            fh.write("{}")
        with gzip.open(os.path.join(d_ind, "logchr22.txt.gz"), "wt") as fh:
            fh.write("chr22\n")
        self.idx, self.ind, self.dicts = idx, ind, [d_snp, d_ind]

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def test_download_routes_gzipped_dicts(self):
        tarball = os.path.join(self.root, f"{self.name}.tar.gz")
        hf._make_index_tarball(
            tarball, self.idx, self.name, self.ind, {"name": self.name}, self.dicts
        )
        names = tarfile.open(tarball).getnames()
        self.assertTrue(
            any(n.endswith("my_dict_chr22.json.gz") for n in names),
            f"gzipped SNP dict not bundled: {names}",
        )

        workdir = tempfile.mkdtemp(prefix="gzdst_")
        self.addCleanup(shutil.rmtree, workdir, True)

        def fake_snapshot(repo, allow_patterns, local_dir, token=None):
            dest = os.path.join(local_dir, "indexes")
            os.makedirs(dest, exist_ok=True)
            shutil.copy(tarball, os.path.join(dest, self.name + ".tar.gz"))
            return local_dir

        orig = hf._hf_snapshot
        hf._hf_snapshot = fake_snapshot
        try:
            hf.download_component(
                "index", workdir, repo="local/test", index_name=self.name, token="x"
            )
        finally:
            hf._hf_snapshot = orig

        self.assertTrue(
            os.path.isfile(os.path.join(
                workdir, "Dictionaries", "dictionaries_hg38_1000G", "my_dict_chr22.json.gz")),
            "gzipped SNP dict not installed into Dictionaries/",
        )
        self.assertTrue(
            os.path.isfile(os.path.join(
                workdir, "Dictionaries", "log_indels_hg38_1000G", "logchr22.txt.gz")),
            "gzipped indel log not installed into Dictionaries/",
        )


if __name__ == "__main__":
    unittest.main()

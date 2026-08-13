#!/usr/bin/env python3
"""Regression tests for the batteries variant-index bundle: a published variant index
must carry its per-sample dictionaries so a `download --what index` install can run a
variant search WITHOUT the multi-GB source VCFs.

Guards the alpha bug where a variant search on a batteries install failed at the genome-
enrichment step (`ERROR: Genome enrichment failed`) because the index-only download
shipped neither the raw VCFs nor the dictionaries the search needs. These tests cover
the packaging half (publish bundles the dicts; download routes them to Dictionaries/);
the pipeline half (submit_job skips enrichment + the indel pipeline derives chroms from
the bundled logs) is covered by the end-to-end search test.

Runnable in the lightweight unit-tests CI: crisprme_hf imports huggingface_hub lazily,
and these tests never touch HF (the snapshot fetch is monkeypatched to a local copy).
"""
import os
import shutil
import tarfile
import tempfile
import unittest

import crisprme_hf as hf


def _mklayout(root, index_name):
    """Create a minimal built variant index + its dictionaries under `root`."""
    vcf = index_name.partition("+")[2]
    idx = os.path.join(root, "genome_library", index_name)
    ind = idx + "_INDELS"
    d_snp = os.path.join(root, "Dictionaries", f"dictionaries_{vcf}")
    d_ind = os.path.join(root, "Dictionaries", f"log_indels_{vcf}")
    for d in (idx, ind, d_snp, d_ind):
        os.makedirs(d)
    for path, content in (
        (os.path.join(idx, "TSTgenome.NGG.bin"), "x"),
        (os.path.join(ind, "NGG_2_fakechr22"), "x"),
        (os.path.join(d_snp, "my_dict_chr22.json"), "{}"),
        (os.path.join(d_ind, "logchr22.txt"), "chr22\n"),
    ):
        with open(path, "w") as fh:
            fh.write(content)
    return idx, ind, [d_snp, d_ind]


class TestIndexBundle(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="bundle_src_")
        self.name = "NGG_2_hg38+hg38_1000G"
        self.idx, self.ind, self.dicts = _mklayout(self.root, self.name)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def test_tarball_bundles_dictionaries(self):
        """_make_index_tarball must include the index, _INDELS, and the dicts under a
        Dictionaries/ prefix (so the extractor can route them separately)."""
        tarball = os.path.join(self.root, f"{self.name}.tar.gz")
        hf._make_index_tarball(tarball, self.idx, self.name, self.ind, {"name": self.name}, self.dicts)
        names = tarfile.open(tarball).getnames()
        self.assertTrue(any(n.startswith(f"{self.name}/") for n in names), "index dir missing")
        self.assertTrue(any(n.startswith(f"{self.name}_INDELS/") for n in names), "_INDELS missing")
        self.assertTrue(any(n.startswith("Dictionaries/dictionaries_hg38_1000G") for n in names),
                        f"SNP dicts not bundled: {names}")
        self.assertTrue(any(n.startswith("Dictionaries/log_indels_hg38_1000G") for n in names),
                        f"indel logs not bundled: {names}")

    def test_download_routes_dicts_to_dictionaries(self):
        """download_component('index') must install the index + _INDELS into
        genome_library/ and the bundled dicts into <workdir>/Dictionaries/."""
        tarball = os.path.join(self.root, f"{self.name}.tar.gz")
        hf._make_index_tarball(tarball, self.idx, self.name, self.ind, {"name": self.name}, self.dicts)

        workdir = tempfile.mkdtemp(prefix="bundle_dst_")
        try:
            def fake_snapshot(repo, allow_patterns, local_dir, token=None):
                dest = os.path.join(local_dir, "indexes")
                os.makedirs(dest, exist_ok=True)
                shutil.copy(tarball, os.path.join(dest, self.name + ".tar.gz"))
                return local_dir
            orig = hf._hf_snapshot
            hf._hf_snapshot = fake_snapshot
            try:
                hf.download_component("index", workdir, repo="local/test",
                                      index_name=self.name, token="x")
            finally:
                hf._hf_snapshot = orig

            self.assertTrue(os.path.isdir(os.path.join(workdir, "genome_library", self.name)),
                            "index not installed into genome_library/")
            self.assertTrue(os.path.isdir(os.path.join(workdir, "genome_library", self.name + "_INDELS")),
                            "_INDELS not installed")
            self.assertTrue(os.path.isfile(os.path.join(
                workdir, "Dictionaries", "dictionaries_hg38_1000G", "my_dict_chr22.json")),
                "SNP dicts not installed into Dictionaries/")
            self.assertTrue(os.path.isfile(os.path.join(
                workdir, "Dictionaries", "log_indels_hg38_1000G", "logchr22.txt")),
                "indel logs not installed into Dictionaries/")
            # the dicts must NOT be left under genome_library/
            self.assertFalse(os.path.isdir(os.path.join(workdir, "genome_library", "Dictionaries")),
                             "Dictionaries leaked into genome_library/")
        finally:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

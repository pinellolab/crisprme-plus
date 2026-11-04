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
import contextlib
import io
import os
import shutil
import stat
import tarfile
import tempfile
import unittest
import warnings
from unittest import mock

import crisprme_hf as hf
import safe_tar


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
            def fake_snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
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


class TestADownloadedArchiveGoesInWithSafeModes(unittest.TestCase):
    """download_component unpacks a downloaded index, its genotype store and a build-inputs
    kit into the data folder, which the web app's Settings download too: tarfile gave each
    file the mode its member has, setuid and setgid too, and, run as root (as the Docker
    image runs, on a data folder mounted from the host), the owner it names. Each goes in
    with safe modes, owned by the user unpacking it."""

    NAME = "NGG_2_hg38+hg38_1000G"

    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="modes_")
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.workdir = os.path.join(self.root, "data")
        os.makedirs(self.workdir)
        self.served = {}  # path in the repo -> the archive it is

    def serve(self, remote, entries):
        """The repo holds at ``remote`` a .tar.gz of ``entries``, (name, mode, payload)
        tuples (payload None for a folder), each member owned by another user (1234)."""
        path = os.path.join(self.root, os.path.basename(remote))
        with tarfile.open(path, "w:gz") as tar:
            for name, mode, payload in entries:
                info = tarfile.TarInfo(name)
                info.mode, info.uid, info.gid = mode, 1234, 1234
                if payload is None:
                    info.type = tarfile.DIRTYPE
                    tar.addfile(info)
                else:
                    info.size = len(payload)
                    tar.addfile(info, io.BytesIO(payload))
        self.served[remote] = path

    @contextlib.contextmanager
    def as_root(self, fallback):
        """download_component's fetches served from self.served (nothing else is reached),
        as root, the owners tarfile sets recorded instead, as (path, uid, gid), in the list
        yielded; ``fallback``: as on Python before 3.11.4, with no 'data' filter."""
        def snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            for pattern in allow_patterns:
                if pattern in self.served:
                    os.makedirs(os.path.dirname(os.path.join(local_dir, pattern)), exist_ok=True)
                    shutil.copy(self.served[pattern], os.path.join(local_dir, pattern))
            return local_dir

        def no_hub():
            raise ImportError("no huggingface_hub here")

        owners = []
        record = lambda path, uid, gid, *a, **kw: owners.append((path, uid, gid))  # noqa: E731
        with mock.patch.object(hf, "_hf_snapshot", snapshot), \
                mock.patch.object(hf, "_require_hf", no_hub), \
                mock.patch("os.geteuid", return_value=0), \
                mock.patch("os.chown", record), mock.patch("os.lchown", record), \
                mock.patch.object(safe_tar, "_DATA_FILTER",
                                  None if fallback else safe_tar._DATA_FILTER), \
                warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)  # (3.12's, for no filter)
            yield owners

    def modes(self, folder):
        """{path under folder: its mode bits} of each file and folder there."""
        return {os.path.relpath(os.path.join(top, name), folder):
                stat.S_IMODE(os.stat(os.path.join(top, name)).st_mode)
                for top, dirs, files in os.walk(folder) for name in dirs + files}

    def test_an_index_and_its_genotype_store(self):
        self.serve(f"indexes/{self.NAME}.tar.gz", [
            (self.NAME, 0o3777, None),
            (f"{self.NAME}/TSTgenome.NGG.bin", 0o4755, b"x"),
            (f"{self.NAME}_INDELS/NGG_2_fakechr22", 0o2666, b"x"),
            ("Dictionaries/dictionaries_hg38_1000G/my_dict_chr22.json", 0o6644, b"{}")])
        self.serve("indexes/genotypes_hg38_1000G.tar.gz", [
            ("Dictionaries/genotypes_hg38_1000G/chr22.gt", 0o4755, b"gt")])
        for fallback in (False, True):
            with self.subTest(fallback=fallback):
                shutil.rmtree(self.workdir)
                os.makedirs(self.workdir)
                with self.as_root(fallback) as owners:
                    hf.download_component("index", self.workdir, repo="local/test",
                                          index_name=self.NAME, token="x")
                library = os.path.join(self.workdir, "genome_library")
                dicts = os.path.join(self.workdir, "Dictionaries")
                self.assertEqual(
                    (self.modes(os.path.join(library, self.NAME)),
                     self.modes(os.path.join(library, f"{self.NAME}_INDELS")),
                     self.modes(os.path.join(dicts, "dictionaries_hg38_1000G")),
                     self.modes(os.path.join(dicts, "genotypes_hg38_1000G"))),
                    ({"TSTgenome.NGG.bin": 0o755}, {"NGG_2_fakechr22": 0o644},
                     {"my_dict_chr22.json": 0o644}, {"chr22.gt": 0o755}))
                self.assertFalse(stat.S_IMODE(os.stat(os.path.join(library, self.NAME)).st_mode)
                                 & 0o7022)
                self.assertEqual([o for o in owners if o[1:] != (-1, -1)], [])

    def test_a_build_inputs_kit(self):
        self.serve("build-inputs/X.tar.gz", [
            ("hg38+hg38_X", 0o3777, None), ("hg38+hg38_X/chr22.enriched.fa", 0o4755, b">c\nA\n")])
        for fallback in (False, True):
            with self.subTest(fallback=fallback):
                shutil.rmtree(self.workdir)
                os.makedirs(self.workdir)
                with self.as_root(fallback) as owners:
                    hf.download_component("build-inputs", self.workdir, repo="local/test",
                                          dataset="X", token="x")
                kit = os.path.join(self.workdir, "Genomes", "hg38+hg38_X")
                self.assertEqual(self.modes(kit), {"chr22.enriched.fa": 0o755})
                self.assertFalse(stat.S_IMODE(os.stat(kit).st_mode) & 0o7022)
                self.assertEqual([o for o in owners if o[1:] != (-1, -1)], [])


if __name__ == "__main__":
    unittest.main()

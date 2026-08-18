#!/usr/bin/env python3
"""Regression tests for the ADDITIVE dictless publish/download of the Tier-0
registry + Tier-1 genotype tiers.

The classic per-sample SNP dictionaries (dictionaries_<vcf>/, ~152GB) are being
replaced by two additive siblings under Dictionaries/:
    registry_<vcf>/reg_<chrom>.{bin,idx}   (Tier-0, small; out-of-the-box detection)
    genotypes_<vcf>/gt_<chrom>.{bin,idx}   (Tier-1, ~22GB; per-sample Samples)
These tests guard the packaging/transport half added to crisprme_hf:

  * DICTLESS publish -> main tarball CONTAINS registry_<vcf>/ + log_indels_<vcf>/
    and EXCLUDES dictionaries_<vcf>/, and a SEPARATE genotypes_<vcf>.tar.gz is
    produced + uploaded (a second upload_file call).
  * DEFAULT publish -> main-tarball MEMBER NAMES are byte-for-byte the classic set
    (dicts present) plus registry when it exists; genotypes never rides inside the
    main tarball (golden member-list assertion).
  * download extracts registry_ from the main tarball AND fetches + places
    genotypes_<vcf>/ into Dictionaries/; --no-genotypes / a missing gt tarball
    warns + continues with genotypes_ absent.
  * a dictless index (no dictionaries_) downloads without error.

Fully offline: the HF client (HfApi/upload_file/snapshot_download) is stubbed, and
crisprme_hf imports huggingface_hub lazily, so importing this test needs no HF dep.
"""
import os
import shutil
import tarfile
import tempfile
import unittest

import crisprme_hf as hf


# --- helpers ----------------------------------------------------------------

def _mkfile(path, content="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(content)


def _mklayout(root, index_name, registry=True, genotypes=True, dicts=True):
    """Build a variant index + a chosen subset of the Dictionaries/ siblings.

    Returns the index dir path. Mirrors the installed layout that
    new_simple_analysis._resolve_registry_paths/_resolve_genotype_paths read.
    """
    vcf = index_name.partition("+")[2]
    idx = os.path.join(root, "genome_library", index_name)
    ind = idx + "_INDELS"
    _mkfile(os.path.join(idx, "TSTgenome.NGG.bin"))
    _mkfile(os.path.join(ind, "NGG_2_fakechr22"))
    droot = os.path.join(root, "Dictionaries")
    _mkfile(os.path.join(droot, f"log_indels_{vcf}", "logchr22.txt"), "chr22\n")
    if dicts:
        _mkfile(os.path.join(droot, f"dictionaries_{vcf}", "my_dict_chr22.json"), "{}")
    if registry:
        _mkfile(os.path.join(droot, f"registry_{vcf}", "reg_chr22.bin"))
        _mkfile(os.path.join(droot, f"registry_{vcf}", "reg_chr22.idx"))
    if genotypes:
        _mkfile(os.path.join(droot, f"genotypes_{vcf}", "gt_chr22.bin"))
        _mkfile(os.path.join(droot, f"genotypes_{vcf}", "gt_chr22.idx"))
    return idx


class _FakeApi:
    """Records upload_file(path_or_fileobj, path_in_repo) calls, copying each
    uploaded tarball aside so a later download test can consume it."""

    def __init__(self, sink_dir):
        self.sink_dir = sink_dir
        self.uploads = []  # list of (local_copy_path, remote_path)

    def upload_file(self, path_or_fileobj, path_in_repo, repo_id, repo_type, token):
        dst = os.path.join(self.sink_dir, os.path.basename(path_in_repo))
        shutil.copy(path_or_fileobj, dst)
        self.uploads.append((dst, path_in_repo))


class _FakeHF:
    def __init__(self, api):
        self._api = api

    def HfApi(self):
        return self._api


def _publish(index_dir, sink_dir, dictless=False):
    """Run publish_index fully offline; return the _FakeApi (with .uploads)."""
    api = _FakeApi(sink_dir)
    orig_require, orig_token = hf._require_hf, hf.resolve_token
    hf._require_hf = lambda: _FakeHF(api)
    hf.resolve_token = lambda t=None: "faketoken"
    try:
        hf.publish_index(index_dir, repo="local/test", token="x", dictless=dictless)
    finally:
        hf._require_hf, hf.resolve_token = orig_require, orig_token
    return api


# --- publish ----------------------------------------------------------------

class TestDictlessPublish(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="dictless_src_")
        self.sink = tempfile.mkdtemp(prefix="dictless_hf_")
        self.name = "NGG_2_hg38+hg38_1000G"
        self.vcf = self.name.partition("+")[2]

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        shutil.rmtree(self.sink, ignore_errors=True)

    def test_dictless_publish_main_tarball_and_gt_companion(self):
        """(a) DICTLESS: main tarball has registry_ + log_indels_, NO
        dictionaries_; a separate genotypes_<vcf>.tar.gz is uploaded."""
        idx = _mklayout(self.root, self.name)
        api = _publish(idx, self.sink, dictless=True)

        remotes = [r for _, r in api.uploads]
        self.assertIn("indexes/NGG_2_hg38+hg38_1000G.tar.gz", remotes,
                      f"main tarball not uploaded: {remotes}")
        self.assertIn(f"indexes/genotypes_{self.vcf}.tar.gz", remotes,
                      f"genotype companion not uploaded separately: {remotes}")

        main = next(p for p, r in api.uploads if r.endswith(f"{self.name}.tar.gz")
                    and "genotypes_" not in r)
        members = tarfile.open(main).getnames()
        self.assertTrue(
            any(m.startswith(f"Dictionaries/registry_{self.vcf}") for m in members),
            f"registry_ missing from dictless main tarball: {members}")
        self.assertTrue(
            any(m.startswith(f"Dictionaries/log_indels_{self.vcf}") for m in members),
            f"log_indels_ missing from dictless main tarball: {members}")
        self.assertFalse(
            any(m.startswith(f"Dictionaries/dictionaries_{self.vcf}") for m in members),
            f"SNP dicts must be EXCLUDED in dictless mode: {members}")
        # the genotype store must NOT ride inside the main tarball
        self.assertFalse(
            any(m.startswith(f"Dictionaries/genotypes_{self.vcf}") for m in members),
            f"genotypes_ must NOT be inside the main tarball: {members}")

        # the separate companion carries the genotype store, routed under Dictionaries/
        gt = next(p for p, r in api.uploads if r.endswith(f"genotypes_{self.vcf}.tar.gz"))
        gt_members = tarfile.open(gt).getnames()
        self.assertTrue(
            any(m.startswith(f"Dictionaries/genotypes_{self.vcf}") for m in gt_members),
            f"genotype companion mislaid: {gt_members}")

    def test_default_publish_member_names_unchanged(self):
        """(b) DEFAULT (classic): main-tarball MEMBER NAMES are exactly the
        pre-change set (index + _INDELS + dicts + indel logs + manifest) plus the
        additive registry when present; genotypes never inside the main tarball."""
        # classic layout WITHOUT registry/genotypes: must match the historical set
        idx = _mklayout(self.root, self.name, registry=False, genotypes=False)
        api = _publish(idx, self.sink, dictless=False)
        self.assertEqual(len(api.uploads), 1,
                         "no genotype dir => only the main tarball is uploaded")
        main = api.uploads[0][0]
        got = set(tarfile.open(main).getnames())
        golden = {
            self.name,
            f"{self.name}/TSTgenome.NGG.bin",
            f"{self.name}_INDELS",
            f"{self.name}_INDELS/NGG_2_fakechr22",
            f"Dictionaries/dictionaries_{self.vcf}",
            f"Dictionaries/dictionaries_{self.vcf}/my_dict_chr22.json",
            f"Dictionaries/log_indels_{self.vcf}",
            f"Dictionaries/log_indels_{self.vcf}/logchr22.txt",
            "manifest.json",
        }
        self.assertEqual(got, golden,
                         f"default main tarball member names changed:\n got={sorted(got)}\n"
                         f" golden={sorted(golden)}")

    def test_default_publish_adds_registry_and_gt_companion(self):
        """DEFAULT with registry + genotypes present: dicts still bundled (classic),
        registry_ ADDED to the main tarball, genotypes_ uploaded separately."""
        idx = _mklayout(self.root, self.name, registry=True, genotypes=True)
        api = _publish(idx, self.sink, dictless=False)
        remotes = [r for _, r in api.uploads]
        self.assertIn(f"indexes/genotypes_{self.vcf}.tar.gz", remotes)
        main = next(p for p, r in api.uploads
                    if r.endswith(f"{self.name}.tar.gz") and "genotypes_" not in r)
        members = tarfile.open(main).getnames()
        self.assertTrue(any(m.startswith(f"Dictionaries/dictionaries_{self.vcf}") for m in members),
                        "default mode must still bundle the SNP dicts")
        self.assertTrue(any(m.startswith(f"Dictionaries/registry_{self.vcf}") for m in members),
                        "registry_ must be ADDED to the default main tarball when present")

    def test_publish_no_genotype_dir_uploads_only_main(self):
        """No genotypes_<vcf>/ present => no companion produced (registry-only /
        classic publish), and a missing registry is simply not added."""
        idx = _mklayout(self.root, self.name, registry=False, genotypes=False)
        api = _publish(idx, self.sink, dictless=False)
        self.assertEqual([r for _, r in api.uploads],
                         [f"indexes/{self.name}.tar.gz"])


# --- download ---------------------------------------------------------------

class TestDictlessDownload(unittest.TestCase):
    """Drives download_component('index') against locally-built tarballs, stubbing
    _hf_snapshot to serve them (mirrors test_index_bundle)."""

    def setUp(self):
        self.src = tempfile.mkdtemp(prefix="dictless_dlsrc_")
        self.sink = tempfile.mkdtemp(prefix="dictless_dlhf_")
        self.workdir = tempfile.mkdtemp(prefix="dictless_dldst_")
        self.name = "NGG_2_hg38+hg38_1000G"
        self.vcf = self.name.partition("+")[2]
        self._orig_snapshot = hf._hf_snapshot

    def tearDown(self):
        hf._hf_snapshot = self._orig_snapshot
        for d in (self.src, self.sink, self.workdir):
            shutil.rmtree(d, ignore_errors=True)

    def _publish_to_sink(self, _dictless=False, **layout_kw):
        idx = _mklayout(self.src, self.name, **layout_kw)
        return _publish(idx, self.sink, dictless=_dictless)

    def _serve_from_sink(self):
        """snapshot stub: copy any requested indexes/* file out of the sink."""
        def fake_snapshot(repo, allow_patterns, local_dir, token=None):
            dest = os.path.join(local_dir, "indexes")
            os.makedirs(dest, exist_ok=True)
            for pat in allow_patterns:
                base = os.path.basename(pat)  # e.g. NGG_..tar.gz or genotypes_...tar.gz
                srcf = os.path.join(self.sink, base)
                if os.path.isfile(srcf):
                    shutil.copy(srcf, os.path.join(dest, base))
            return local_dir
        hf._hf_snapshot = fake_snapshot

    def test_download_installs_registry_and_genotypes(self):
        """(c) download extracts registry_ from the main tarball AND fetches +
        places genotypes_<vcf>/ into Dictionaries/ (sibling of dictionaries_)."""
        self._publish_to_sink(registry=True, genotypes=True)
        self._serve_from_sink()
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, token="x", genotypes=True)
        dd = os.path.join(self.workdir, "Dictionaries")
        self.assertTrue(os.path.isfile(os.path.join(dd, f"registry_{self.vcf}", "reg_chr22.bin")),
                        "registry_ not extracted from the main tarball")
        self.assertTrue(os.path.isfile(os.path.join(dd, f"genotypes_{self.vcf}", "gt_chr22.bin")),
                        "genotypes_ not fetched + placed as a Dictionaries/ sibling")
        self.assertTrue(os.path.isfile(os.path.join(dd, f"genotypes_{self.vcf}", "gt_chr22.idx")))
        # the index + _INDELS still land in genome_library/
        self.assertTrue(os.path.isdir(os.path.join(self.workdir, "genome_library", self.name)))

    def test_download_no_genotypes_flag_skips_store(self):
        """(c) with genotypes=False (--no-genotypes): warn + continue; genotypes_
        absent, but registry_ + index still installed."""
        self._publish_to_sink(registry=True, genotypes=True)
        self._serve_from_sink()
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, token="x", genotypes=False)
        dd = os.path.join(self.workdir, "Dictionaries")
        self.assertFalse(os.path.isdir(os.path.join(dd, f"genotypes_{self.vcf}")),
                         "--no-genotypes must NOT install the genotype store")
        self.assertTrue(os.path.isfile(os.path.join(dd, f"registry_{self.vcf}", "reg_chr22.bin")),
                        "registry_ should still install with --no-genotypes")

    def test_download_missing_gt_tarball_warns_and_continues(self):
        """(c) a registry-only publish (no genotypes_ companion): default download
        warns + continues; genotypes_ absent, no error."""
        self._publish_to_sink(registry=True, genotypes=False)
        self._serve_from_sink()
        # must not raise even though no genotypes_<vcf>.tar.gz exists in the sink
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, token="x", genotypes=True)
        dd = os.path.join(self.workdir, "Dictionaries")
        self.assertFalse(os.path.isdir(os.path.join(dd, f"genotypes_{self.vcf}")))
        self.assertTrue(os.path.isfile(os.path.join(dd, f"registry_{self.vcf}", "reg_chr22.bin")))

    def test_download_dictless_index_no_dicts_ok(self):
        """(d) a dictless index (no dictionaries_ in the main tarball) downloads
        without error: registry + genotype tiers cover it."""
        api = self._publish_to_sink(registry=True, genotypes=True, _dictless=True)
        # sanity: the published main tarball really is dictless
        main = next(p for p, r in api.uploads
                    if r.endswith(f"{self.name}.tar.gz") and "genotypes_" not in r)
        self.assertFalse(
            any(m.startswith(f"Dictionaries/dictionaries_{self.vcf}")
                for m in tarfile.open(main).getnames()))
        self._serve_from_sink()
        hf.download_component("index", self.workdir, repo="local/test",
                              index_name=self.name, token="x", genotypes=True)
        dd = os.path.join(self.workdir, "Dictionaries")
        self.assertFalse(os.path.isdir(os.path.join(dd, f"dictionaries_{self.vcf}")),
                         "dictless index must not ship the per-sample SNP dicts")
        self.assertTrue(os.path.isfile(os.path.join(dd, f"registry_{self.vcf}", "reg_chr22.bin")))
        self.assertTrue(os.path.isfile(os.path.join(dd, f"genotypes_{self.vcf}", "gt_chr22.bin")))
        self.assertTrue(os.path.isfile(os.path.join(dd, f"log_indels_{self.vcf}", "logchr22.txt")),
                        "dictless mode must KEEP the indel logs")


if __name__ == "__main__":
    unittest.main()

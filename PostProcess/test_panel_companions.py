#!/usr/bin/env python3
"""CRISPRme 2.10.0's indexes on Hugging Face: a location of their own (decision C2), and each
panel's Dictionaries bundle and genotype store published and fetched once per panel.

* Routing: an index of index format 2, or one whose name carries its guide length and side
  (20-NRR_3_hg38, TTTN-23_3_hg38), goes to and comes from crisprme_hf.FORMAT2_HF_REPO, which no
  release before 2.10.0 reads (2.8.1 refuses an index of format 2, 2.8.0 aborts on one); an
  archive of 2.8.1 (index format 1, the older name) is still fetched by its own name from the
  repos of 2.8.1, and the download names the index of this release that replaces it (the index
  set, seq_script/rebuild/index_set_2.10.0.tsv). An index of format 2 is never published in
  the repos 2.8.1 reads.
* The catalog (the web's Settings list) offers the indexes of this release alone, and no
  companion as an index.
* publish-index of a variant index of format 2 leaves the panel's Dictionaries out of the
  index's archive and uploads them, once per panel, in dictionaries-bundle_<vcf>.tar.gz, as it
  does the genotype store (genotypes_<vcf>) and the phased indel genotypes
  (indel_genotypes_<vcf>), each with a record of what it holds (<name>.json: a fingerprint of
  its files); a companion the repo holds already, the same, is not uploaded again, and one the
  repo holds in another build is refused (the panel's indexes there read it) unless
  replace_companions. Before 2.10.0 every variant archive carried its own copy of the
  Dictionaries (16.5 GB per HGDP + 1000 Genomes archive), each install replaced them, and
  publish-index re-uploaded the genotype store with every index.
* download fetches the bundle and the genotype store that an index's manifest names once per
  panel: kept, not fetched again, where the install holds them with the same fingerprint;
  checked against the manifest's fingerprint; refused over the Dictionaries of another build
  where an archive of 2.8.1 would replace a bundle.

Hermetic: huggingface_hub is replaced through crisprme_hf._require_hf by a fake hub that keeps
its repos in a temporary folder and records every call. No network, no token.
"""

import gzip
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

import crisprme_hf as hf

# (read so, the tests run, and each fails, on a crisprme_hf that predates the location)
FMT2 = getattr(hf, "FORMAT2_HF_REPO", "lucapinello/crisprme_plus_indices")
DEFAULT = hf.DEFAULT_HF_REPO
GENERAL = hf.GENERAL_INDEX_HF_REPO
VCF = "hg38_HGDP1kGP"
NRR = f"20-NRR_3_hg38+{VCF}"
NGG = f"20-NGG_3_hg38+{VCF}"


class EntryNotFoundError(Exception):
    """What the hub raises for a file a repo does not hold (its class name is what
    crisprme_hf._hf_absent reads)."""


class _Entry:
    def __init__(self, path, size):
        self.path, self.size = path, size


class FakeHub:
    """Repos as folders under ``root``; records (call, repo, detail) for every hub call."""

    def __init__(self, root):
        self.root = root
        self.calls = []
        hub = self

        class HfApi:
            def __init__(self, *a, **k):
                pass

            def upload_file(self, path_or_fileobj, path_in_repo, repo_id, repo_type=None,
                            token=None, **k):
                hub.calls.append(("upload_file", repo_id, path_in_repo))
                dst = hub._path(repo_id, path_in_repo)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copyfile(path_or_fileobj, dst)

            def list_repo_tree(self, repo_id, repo_type=None, path_in_repo=None,
                               recursive=False, token=None, **k):
                hub.calls.append(("list_repo_tree", repo_id, path_in_repo))
                base = hub._repo(repo_id)
                folder = os.path.join(base, path_in_repo.strip("/"))
                if not os.path.isdir(folder):
                    return []
                return [_Entry(f"{path_in_repo.strip('/')}/{n}",
                               os.path.getsize(os.path.join(folder, n)))
                        for n in sorted(os.listdir(folder))
                        if os.path.isfile(os.path.join(folder, n))]

            def file_exists(self, repo_id, filename, repo_type=None, token=None, **k):
                hub.calls.append(("file_exists", repo_id, filename))
                return os.path.isfile(hub._path(repo_id, filename))

        self.HfApi = HfApi

    def _repo(self, repo_id):
        base = os.path.join(self.root, repo_id.replace("/", "__"))
        if not os.path.isdir(base):
            raise RuntimeError(f"Repository Not Found: {repo_id}")
        return base

    def _path(self, repo_id, path):
        return os.path.join(self.root, repo_id.replace("/", "__"), path)

    def make_repo(self, repo_id):
        os.makedirs(os.path.join(self.root, repo_id.replace("/", "__")), exist_ok=True)

    def snapshot_download(self, repo_id, repo_type=None, allow_patterns=None, local_dir=None,
                          token=None, **k):
        self.calls.append(("snapshot_download", repo_id, tuple(allow_patterns)))
        base = self._repo(repo_id)
        for path in allow_patterns:
            src = os.path.join(base, path)
            if os.path.isfile(src):
                dst = os.path.join(local_dir, path)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copyfile(src, dst)
        return local_dir

    def get_paths_info(self, repo_id, paths, repo_type=None, token=None, **k):
        base = self._repo(repo_id)
        return [_Entry(p, os.path.getsize(os.path.join(base, p))) for p in paths
                if os.path.isfile(os.path.join(base, p))]

    def hf_hub_download(self, repo_id, filename, repo_type=None, token=None, **k):
        self.calls.append(("hf_hub_download", repo_id, filename))
        src = os.path.join(self._repo(repo_id), filename)
        if not os.path.isfile(src):
            raise EntryNotFoundError(f"404 Client Error. Entry Not Found for url: {filename}")
        dst = os.path.join(tempfile.mkdtemp(dir=self.root), os.path.basename(filename))
        shutil.copyfile(src, dst)
        return dst

    def uploads(self, since=0):
        return [c[2] for c in self.calls[since:] if c[0] == "upload_file"]

    def fetched(self, since=0):
        return [p for c in self.calls[since:] if c[0] == "snapshot_download" for p in c[2]]

    def repos_contacted(self):
        return {c[1] for c in self.calls}


def _write(path, data=b"x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)


def _tree_bin(path, fmt=2):
    """A file index_fit.bin_format reads as an index tree of format ``fmt``."""
    head = (b"CZI\xfe" + (2).to_bytes(4, sys.byteorder, signed=True)) if fmt == 2 \
        else (5).to_bytes(4, sys.byteorder, signed=True) + (20).to_bytes(4, sys.byteorder)
    _write(path, head + b"\0" * 24)


def _gz(path, text, mtime):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=mtime) as fh:
            fh.write(text.encode())


def make_panel(root, names, vcf=VCF, registry=b"reg-1", fmt=2, genotypes=True, mtime=1):
    """A working directory holding the variant indexes ``names`` of panel ``vcf`` (their
    _INDELS indexes and the shared reference index) and the panel's Dictionaries and
    samplesID list, as build-index-only leaves them."""
    for name in names:
        pam = hf.index_fit.index_pam(name)
        _tree_bin(os.path.join(root, "genome_library", name, f"{pam}_chr22_1.bin"), fmt)
        _tree_bin(os.path.join(root, "genome_library", name + "_INDELS",
                               f"{pam}_fakechr22_1.bin"), fmt)
        base = name.partition("+")[0]
        _tree_bin(os.path.join(root, "genome_library", base, f"{pam}_chr22_1.bin"), fmt)
    d = os.path.join(root, "Dictionaries")
    _gz(os.path.join(d, f"log_indels_{vcf}", "logchr22.txt.gz"), "chr22\t1\n", mtime)
    _write(os.path.join(d, f"registry_{vcf}", "reg_chr22.bin"), registry)
    _write(os.path.join(d, f"registry_{vcf}", "reg_chr22.idx"), b"idx")
    _gz(os.path.join(d, f"indel_genotypes_{vcf}", "gt_indel_chr22.tsv.gz"), "S1\t1|0\n", mtime)
    if genotypes:
        _write(os.path.join(d, f"genotypes_{vcf}", "gt_chr22.bin"), b"gt")
    _write(os.path.join(root, "samplesIDs", f"{vcf}.samplesID.txt"), b"#SAMPLE_ID\nS1\n")
    return root


class _Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="panelcomp_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        saved = os.environ.pop("CRISPRME_HF_REPO", None)
        self.addCleanup(lambda: saved is None or os.environ.__setitem__("CRISPRME_HF_REPO", saved))
        self.hub = FakeHub(os.path.join(self.tmp, "hub"))
        for repo in (FMT2, DEFAULT, GENERAL):
            self.hub.make_repo(repo)

    def publish(self, index_dir, **kw):
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()):
            hf.publish_index(index_dir, token="t", **kw)
        return out.getvalue()

    def download(self, workdir, name, **kw):
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()):
            hf.download_component("index", workdir, index_name=name, token=None, **kw)
        return out.getvalue()

    def remote(self, repo, path):
        return self.hub._path(repo, path)

    def members(self, repo, path):
        with tarfile.open(self.remote(repo, path)) as tf:
            return sorted(tf.getnames())

    def manifest(self, repo, name):
        with tarfile.open(self.remote(repo, f"indexes/{name}.tar.gz")) as tf:
            return json.load(tf.extractfile("manifest.json"))


class Routing(_Base):
    def test_a_name_of_2_10_0_goes_to_the_format_2_repo(self):
        for name in ("20-NRR_3_hg38", NRR, "TTTN-23_3_hg38+" + VCF, "22-NNGRRN_3_hg38",
                     "20-NNNNNNNN_3_hg38+" + VCF):
            self.assertEqual(hf.resolve_index_repo(name), FMT2, name)
            self.assertTrue(hf.is_named_for_guide(name), name)

    def test_an_index_of_format_2_goes_there_whatever_its_name(self):
        self.assertEqual(hf.resolve_index_repo("NRR_3_hg38", index_format=2), FMT2)
        self.assertEqual(hf.resolve_index_repo("NNNNNNNN_3_hg38", index_format=2), FMT2)

    def test_the_archives_of_2_8_1_keep_their_repos(self):
        self.assertEqual(hf.resolve_index_repo("NRR_3_hg38+" + VCF), DEFAULT)
        self.assertEqual(hf.resolve_index_repo("NRR_3_hg38", index_format=1), DEFAULT)
        self.assertEqual(hf.resolve_index_repo("NNNNNNNN_3_hg38"), GENERAL)
        self.assertFalse(hf.is_named_for_guide("NRR_3_hg38+" + VCF))

    def test_a_repo_named_still_wins(self):
        self.assertEqual(hf.resolve_index_repo(NRR, "me/mine"), "me/mine")

    def test_a_kit_comes_from_the_format_2_repo(self):
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                self.assertRaises(ValueError) as caught:
            hf.download_component("build-inputs", self.tmp, dataset=VCF)
        self.assertEqual({c[1] for c in self.hub.calls if c[0] == "snapshot_download"}, {FMT2})
        self.assertIn(FMT2, str(caught.exception))


class TheCatalog(_Base):
    def _list(self, repo=None):
        with mock.patch.object(hf, "_require_hf", return_value=self.hub):
            return hf.list_available_downloads("index", repo)

    def test_it_offers_this_release_s_indexes_and_no_companion(self):
        for path in (f"indexes/{NRR}.tar.gz", "indexes/20-NRR_3_hg38.tar.gz",
                     f"indexes/dictionaries-bundle_{VCF}.tar.gz", f"indexes/genotypes_{VCF}.tar.gz",
                     f"indexes/indel_genotypes_{VCF}.tar.gz",
                     f"indexes/dictionaries-bundle_{VCF}.json"):
            _write(self.remote(FMT2, path))
        _write(self.remote(DEFAULT, f"indexes/NRR_3_hg38+{VCF}.tar.gz"))
        rows = {r["name"]: r for r in self._list()}
        self.assertEqual(sorted(rows), ["20-NRR_3_hg38", NRR])
        self.assertTrue(rows[NRR]["has_bundle"] and rows[NRR]["has_genotypes"])
        self.assertFalse(rows["20-NRR_3_hg38"]["has_bundle"])
        self.assertEqual(self.hub.repos_contacted(), {FMT2})

    def test_an_unreadable_repo_is_the_offline_signal(self):
        # (not a catalog of 2.8.1's archives, which this release does not offer)
        _write(self.remote(DEFAULT, f"indexes/NRR_3_hg38+{VCF}.tar.gz"))
        shutil.rmtree(self.hub._repo(FMT2))
        self.assertEqual(self._list(), [])


class PublishOncePerPanel(_Base):
    def setUp(self):
        super().setUp()
        self.work = make_panel(os.path.join(self.tmp, "work"), [NRR, NGG])

    def test_the_first_index_ships_the_companions_and_carries_no_dictionaries(self):
        said = self.publish(os.path.join(self.work, "genome_library", NRR), dictless=True)
        self.assertEqual(self.hub.uploads(), [
            f"indexes/dictionaries-bundle_{VCF}.tar.gz", f"indexes/dictionaries-bundle_{VCF}.json",
            f"indexes/genotypes_{VCF}.tar.gz", f"indexes/genotypes_{VCF}.json",
            f"indexes/indel_genotypes_{VCF}.tar.gz", f"indexes/indel_genotypes_{VCF}.json",
            f"indexes/{NRR}.tar.gz"])
        self.assertEqual({c[1] for c in self.hub.calls if c[0] == "upload_file"}, {FMT2})
        self.assertIn(f"Published index to {FMT2}:indexes/{NRR}.tar.gz", said)
        names = self.members(FMT2, f"indexes/{NRR}.tar.gz")
        self.assertFalse([n for n in names if n.startswith(("Dictionaries", "samplesIDs"))], names)
        self.assertIn("manifest.json", names)
        bundle = self.members(FMT2, f"indexes/dictionaries-bundle_{VCF}.tar.gz")
        for want in (f"Dictionaries/log_indels_{VCF}/logchr22.txt.gz",
                     f"Dictionaries/registry_{VCF}/reg_chr22.bin",
                     f"Dictionaries/indel_genotypes_{VCF}/gt_indel_chr22.tsv.gz",
                     f"samplesIDs/{VCF}.samplesID.txt", "bundle.json"):
            self.assertIn(want, bundle)
        self.assertFalse([n for n in bundle if "/genotypes_" in n], bundle)
        manifest = self.manifest(FMT2, NRR)
        record = json.load(open(self.remote(FMT2, f"indexes/dictionaries-bundle_{VCF}.json")))
        self.assertEqual(manifest["dictionaries_bundle"],
                         {"name": f"dictionaries-bundle_{VCF}", "fingerprint": record["fingerprint"]})
        gt = json.load(open(self.remote(FMT2, f"indexes/genotypes_{VCF}.json")))
        self.assertEqual(manifest["genotypes"]["fingerprint"], gt["fingerprint"])
        self.assertEqual(manifest["index_format"], 2)

    def test_the_next_index_of_the_panel_ships_its_own_archive_alone(self):
        self.publish(os.path.join(self.work, "genome_library", NRR), dictless=True)
        mark = len(self.hub.calls)
        said = self.publish(os.path.join(self.work, "genome_library", NGG), dictless=True)
        self.assertEqual(self.hub.uploads(mark), [f"indexes/{NGG}.tar.gz"])
        self.assertIn(f"dictionaries-bundle_{VCF} is in {FMT2} already, the same", said)
        self.assertEqual(self.manifest(FMT2, NGG)["dictionaries_bundle"],
                         self.manifest(FMT2, NRR)["dictionaries_bundle"])

    def test_another_build_of_the_panel_is_refused_unless_replaced(self):
        self.publish(os.path.join(self.work, "genome_library", NRR), dictless=True)
        other = make_panel(os.path.join(self.tmp, "other"), [NGG], registry=b"reg-2")
        mark = len(self.hub.calls)
        with self.assertRaises(ValueError) as caught:
            self.publish(os.path.join(other, "genome_library", NGG), dictless=True)
        self.assertIn(f"dictionaries-bundle_{VCF}", str(caught.exception))
        self.assertIn("Nothing was uploaded", str(caught.exception))
        self.assertEqual(self.hub.uploads(mark), [])
        self.publish(os.path.join(other, "genome_library", NGG), dictless=True,
                     replace_companions=True)
        self.assertIn(f"indexes/dictionaries-bundle_{VCF}.tar.gz", self.hub.uploads(mark))

    def test_a_gzip_timestamp_is_no_other_build(self):
        a = make_panel(os.path.join(self.tmp, "a"), [NRR], mtime=1)
        b = make_panel(os.path.join(self.tmp, "b"), [NRR], mtime=99)
        members = lambda root: [(f"Dictionaries/{d}", os.path.join(root, "Dictionaries", d))
                                for d in (f"log_indels_{VCF}", f"registry_{VCF}")]
        fa, fb = hf.companion_fingerprint(members(a)), hf.companion_fingerprint(members(b))
        self.assertEqual(fa["fingerprint"], fb["fingerprint"])
        _write(os.path.join(b, "Dictionaries", f"registry_{VCF}", "reg_chr22.bin"), b"reg-3")
        self.assertNotEqual(fa["fingerprint"], hf.companion_fingerprint(members(b))["fingerprint"])

    def test_an_index_of_format_2_never_goes_to_a_repo_of_2_8_1(self):
        for repo in (DEFAULT, GENERAL):
            with self.assertRaises(ValueError) as caught:
                self.publish(os.path.join(self.work, "genome_library", NRR), repo=repo)
            self.assertIn("decision C2", str(caught.exception))
        self.assertEqual(self.hub.uploads(), [])

    def test_an_index_of_format_1_keeps_its_dictionaries_and_its_repo(self):
        old = make_panel(os.path.join(self.tmp, "old"), [f"NRR_3_hg38+{VCF}"], fmt=1)
        self.publish(os.path.join(old, "genome_library", f"NRR_3_hg38+{VCF}"), dictless=True)
        self.assertEqual({c[1] for c in self.hub.calls if c[0] == "upload_file"}, {DEFAULT})
        self.assertIn(f"indexes/genotypes_{VCF}.tar.gz", self.hub.uploads())
        names = self.members(DEFAULT, f"indexes/NRR_3_hg38+{VCF}.tar.gz")
        self.assertIn(f"Dictionaries/registry_{VCF}/reg_chr22.bin", names)
        self.assertNotIn("dictionaries_bundle", self.manifest(DEFAULT, f"NRR_3_hg38+{VCF}"))


class DownloadOncePerPanel(_Base):
    def setUp(self):
        super().setUp()
        work = make_panel(os.path.join(self.tmp, "work"), [NRR, NGG])
        for name in (NRR, NGG):
            self.publish(os.path.join(work, "genome_library", name), dictless=True)
        for base in ("20-NRR_3_hg38", "20-NGG_3_hg38"):
            self.publish(os.path.join(work, "genome_library", base))
        self.dest = os.path.join(self.tmp, "dest")
        os.makedirs(self.dest)

    def marker(self, name):
        return json.load(open(hf.companion_marker(self.dest, name)))

    def test_the_first_index_installs_the_panel_s_bundle_and_store(self):
        self.download(self.dest, NRR)
        d = os.path.join(self.dest, "Dictionaries")
        for sub in (f"registry_{VCF}/reg_chr22.bin", f"log_indels_{VCF}/logchr22.txt.gz",
                    f"indel_genotypes_{VCF}/gt_indel_chr22.tsv.gz", f"genotypes_{VCF}/gt_chr22.bin"):
            self.assertTrue(os.path.isfile(os.path.join(d, sub)), sub)
        self.assertTrue(os.path.isfile(os.path.join(self.dest, "samplesIDs", f"{VCF}.samplesID.txt")))
        self.assertTrue(os.path.isdir(os.path.join(self.dest, "genome_library", NRR + "_INDELS")))
        self.assertTrue(os.path.isdir(os.path.join(self.dest, "genome_library", "20-NRR_3_hg38")))
        self.assertEqual(self.marker(f"dictionaries-bundle_{VCF}")["indexes"], [NRR])
        self.assertEqual(self.marker(f"genotypes_{VCF}")["indexes"], [NRR])
        self.assertEqual({c[1] for c in self.hub.calls if c[0] == "snapshot_download"}, {FMT2})

    def test_the_next_index_of_the_panel_fetches_its_archive_alone(self):
        self.download(self.dest, NRR)
        mark = len(self.hub.calls)
        said = self.download(self.dest, NGG)
        fetched = self.hub.fetched(mark)
        self.assertIn(f"indexes/{NGG}.tar.gz", fetched)
        self.assertFalse([p for p in fetched if "dictionaries-bundle_" in p or "genotypes_" in p],
                         fetched)
        self.assertIn(f"dictionaries-bundle_{VCF} is installed already", said)
        self.assertEqual(self.marker(f"dictionaries-bundle_{VCF}")["indexes"], [NGG, NRR])

    def test_no_genotypes_leaves_the_store_out(self):
        self.download(self.dest, NRR, genotypes=False)
        self.assertFalse(os.path.isdir(os.path.join(self.dest, "Dictionaries", f"genotypes_{VCF}")))
        self.assertFalse([p for p in self.hub.fetched() if p.endswith(f"/genotypes_{VCF}.tar.gz")])
        # the bundle is installed all the same, from the format-2 repo
        self.assertEqual(self.marker(f"dictionaries-bundle_{VCF}")["indexes"], [NRR])
        self.assertIn(f"indexes/dictionaries-bundle_{VCF}.tar.gz", self.hub.fetched())

    def test_a_complete_install_costs_no_hub_call(self):
        self.download(self.dest, NRR)
        self.assertIsNotNone(hf.installed_companion(self.dest, f"dictionaries-bundle_{VCF}"))
        mark = len(self.hub.calls)
        self.download(self.dest, NRR)
        self.assertEqual(self.hub.calls[mark:], [])

    def test_a_bundle_that_is_not_the_index_s_is_refused(self):
        record = self.remote(FMT2, f"indexes/dictionaries-bundle_{VCF}.json")
        rec = json.load(open(record))
        rec["fingerprint"] = "0" * 64
        json.dump(rec, open(record, "w"))
        with self.assertRaises(hf.DownloadIntegrityError):
            self.download(self.dest, NRR)
        self.assertFalse(os.path.isdir(os.path.join(self.dest, "Dictionaries", f"registry_{VCF}")))

    def test_an_archive_of_2_8_1_is_refused_over_a_bundle(self):
        self.download(self.dest, NRR)
        old = make_panel(os.path.join(self.tmp, "old"), [f"NRR_3_hg38+{VCF}"], fmt=1,
                         registry=b"reg-281")
        self.publish(os.path.join(old, "genome_library", f"NRR_3_hg38+{VCF}"), dictless=True)
        with self.assertRaises(ValueError) as caught:
            self.download(self.dest, f"NRR_3_hg38+{VCF}")
        self.assertIn("another working directory", str(caught.exception))
        self.assertFalse(os.path.isdir(os.path.join(self.dest, "genome_library",
                                                    f"NRR_3_hg38+{VCF}")))
        with open(os.path.join(self.dest, "Dictionaries", f"registry_{VCF}", "reg_chr22.bin"),
                  "rb") as fh:
            self.assertEqual(fh.read(), b"reg-1")

    def test_a_bundle_over_another_build_s_dictionaries_names_its_indexes(self):
        old = make_panel(os.path.join(self.tmp, "old"), [f"NRR_3_hg38+{VCF}"], fmt=1,
                         registry=b"reg-281")
        self.publish(os.path.join(old, "genome_library", f"NRR_3_hg38+{VCF}"), dictless=True)
        self.publish(os.path.join(old, "genome_library", "NRR_3_hg38"))
        self.download(self.dest, f"NRR_3_hg38+{VCF}")
        said = self.download(self.dest, NRR)
        self.assertIn(f"WARNING: NRR_3_hg38+{VCF}", said)
        self.assertIn("another build", said)
        # and that index is no longer a complete install: downloaded again, it is refused
        with self.assertRaises(ValueError):
            self.download(self.dest, f"NRR_3_hg38+{VCF}")


class AnArchiveOf281NamesItsSuccessor(_Base):
    def _old(self, name):
        old = make_panel(os.path.join(self.tmp, "old"), [name],
                         vcf=name.partition("+")[2] or VCF, fmt=1)
        index = os.path.join(old, "genome_library", name)
        self.publish(index, dictless=True)
        if "+" in name:
            self.publish(os.path.join(old, "genome_library", name.partition("+")[0]))
        return name

    def test_its_build_of_2_10_0(self):
        said = self.download(self.tmp, self._old(f"NRR_3_hg38+{VCF}"))
        self.assertIn(f"This release's build of it is {NRR}, in {FMT2}", said)

    def test_nnn_3_is_served_by_the_general_index(self):
        said = self.download(self.tmp, self._old(f"NNN_3_hg38+{VCF}"))
        self.assertIn(f"It is not built again for 2.10.0: 20-NNNNNNNN_3_hg38+{VCF}", said)
        # (a PAM-less index costs the same whatever its flank: every position is a site)
        self.assertIn("serves every search it serves, pruned to the search's PAM, at the same "
                      "speed", said)
        self.assertNotIn("or faster", said)
        # (the note is read from the index set, offline: the format-2 repo is not asked)
        self.assertNotIn(FMT2, {c[1] for c in self.hub.calls if c[0] != "upload_file"})

    def test_nrg_on_hprc_and_mega_is_served_by_nrr_more_slowly(self):
        """NRG is not built again on HPRC and mega: the panel's NRR index serves its searches,
        pruned to NRG. The note said it did so "at the same speed or faster", which holds of
        the general index for NNN_3 alone: an NRR index holds the sites of NAA, NAG, NGA and
        NGG and prunes at leaf time (the engine's search of chr22's reference sequence, three
        NRG guides at 4 mismatches + 1 DNA + 1 RNA bulge, the same 85,563 rows: 107 s wall and
        380 s CPU, against 57 s and 110 s for an NRG index; the stage-2b review of lane D,
        D-R1-2). The note said "on chr22 about twice the wall time", which reads as a search
        on an HPRC or mega index, whose post-analysis was not measured: it names the setting
        of the measurement (D-R2-4)."""
        for vcf in ("hg38_HPRC", "hg38_mega"):
            name = self._old(f"NRG_3_hg38+{vcf}")
            said = self.download(os.path.join(self.tmp, vcf), name)
            self.assertIn(f"It is not built again for 2.10.0: 20-NRR_3_hg38+{vcf}, in {FMT2}, "
                          "serves every search it serves, pruned to the search's PAM: the same "
                          "rows, more slowly than a dedicated NRG index", said)
            self.assertIn("about twice the wall time in the engine's search of chr22's "
                          "reference sequence", said)
            self.assertNotIn("same speed", said)
            self.assertNotIn("faster", said)
            shutil.rmtree(os.path.join(self.tmp, "old"))

    def test_it_says_so_before_the_archive_is_fetched(self):
        """The note came after the install: after the whole archive of 2.8.1 was fetched
        (tens of GB) for an index the user may rather replace. It now comes first, and a
        download stopped during the fetch has said it."""
        name = self._old(f"NNN_3_hg38+{VCF}")
        out = io.StringIO()

        def stopped(*a, **k):
            raise KeyboardInterrupt

        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch.object(hf, "_hf_snapshot", side_effect=stopped), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()), \
                self.assertRaises(KeyboardInterrupt):
            hf.download_component("index", self.tmp, index_name=name, token=None)
        self.assertIn(f"It is not built again for 2.10.0: 20-NNNNNNNN_3_hg38+{VCF}",
                      out.getvalue())
        # (and once, when the download goes on to the end)
        self.assertEqual(self.download(self.tmp, name).count(f"NOTE: {name} is an index of "
                                                               "CRISPRme 2.8.1"), 1)

    def test_a_name_of_2_10_0_not_yet_published_says_where_it_was_looked_for(self):
        """Asked for an index of this release that its repo does not hold (before the
        rebuild is published), the error described the repo of index format 2 as the one
        "where indexes with a flank wider than 3 nt are published", 2.8.1's general index's."""
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", io.StringIO()), mock.patch("sys.stderr", io.StringIO()), \
                self.assertRaises(ValueError) as caught:
            hf.download_component("index", self.tmp, index_name="20-NRR_3_hg38", token=None)
        said = str(caught.exception)
        self.assertIn(f"No index '20-NRR_3_hg38' was found in HuggingFace repo '{FMT2}' (where "
                      "the indexes of CRISPRme 2.10.0 and later, in index format 2, are "
                      f"published), nor in the fallback '{DEFAULT}'", said)
        self.assertNotIn("flank wider than 3 nt", said)
        # (an older wide-flank name: 2.8.1's general index's repo, as before)
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", io.StringIO()), mock.patch("sys.stderr", io.StringIO()), \
                self.assertRaises(ValueError) as caught:
            hf.download_component("index", self.tmp, index_name="NNNNNNNN_3_hg38", token=None)
        self.assertIn(f"'{GENERAL}' (where indexes with a flank wider than 3 nt are published)",
                      str(caught.exception))

    def test_the_index_set_replaces_every_archive_of_2_8_1(self):
        self.assertEqual(hf.successor_index("NRG_3_hg38+hg38_HPRC"), "20-NRR_3_hg38+hg38_HPRC")
        self.assertEqual(hf.successor_index("NNNNNNNN_3_hg38"), "20-NNNNNNNN_3_hg38")
        self.assertIsNone(hf.successor_index(NRR))


def add_mnv_files(root, vcf=VCF, chroms=("chr22",), mtime=1):
    """The panel's phased indel genotype store gains its MNV files (2.10.1: gt_mnv_<chrom>
    and its marker, as build-index-only's STEP 1c and build_indel_genotypes --mnv-only write
    them), its indel files left as they are."""
    store = os.path.join(root, "Dictionaries", f"indel_genotypes_{vcf}")
    for chrom in chroms:
        _gz(os.path.join(store, f"gt_mnv_{chrom}.tsv.gz"), "120_AA_GG\tS1:1|0\n", mtime)
        _write(os.path.join(store, f"gt_mnv_{chrom}.tsv.gz.done"), b"1\nFILTER\tPASS,.\n")


def _json(path):
    with open(path) as fh:
        return json.load(fh)


def hub_files(hub, repo):
    """{path: bytes} of every file of a repo of the fake hub."""
    base = hub._repo(repo)
    out = {}
    for d, _dirs, names in os.walk(base):
        for n in names:
            with open(os.path.join(d, n), "rb") as fh:
                out[os.path.relpath(os.path.join(d, n), base)] = fh.read()
    return out


class _NoNetwork(object):
    """Any attempt to reach the network fails the test: a socket's connect, a name lookup,
    urllib's urlopen. The fake hub answers in process, so with it in place a staged
    publish and the download step make no network call."""

    def __enter__(self):
        import socket
        import urllib.request

        def refuse(*a, **k):
            raise AssertionError("a network call was made: %r" % (a[:2],))

        self._patches = [mock.patch.object(socket.socket, "connect", refuse),
                         mock.patch.object(socket, "create_connection", refuse),
                         mock.patch.object(socket, "getaddrinfo", refuse),
                         mock.patch.object(urllib.request, "urlopen", refuse)]
        for p in self._patches:
            p.start()
        return self

    def __exit__(self, *exc):
        for p in reversed(self._patches):
            p.stop()
        return False


class StageOnlyPublish(_Base):
    """publish-index --stage-only (2.10.1, the one-commit publish of the plan of 2.10.1's
    MNVs, section 4): everything the publish would upload is packed into a folder, the
    index's archive, each companion it would upload and its record, with the SHA-256 lines
    of the archives, and NOTHING is uploaded: no write token is asked for, no HfApi is made,
    no network is reached; the only calls to the hub are anonymous reads of the companions'
    records (token False, whatever HF_TOKEN says), to keep a companion the repo holds the
    same and refuse one it holds in another build. It prints the commands that publish the
    folder in one commit with indexes/SHA256SUMS, which it does not run. Here the panel's
    store gains its MNV files after 2.10.0's publish, as HPRC's does: the bundle and the
    phased indel genotypes are staged with them, the genotype store (the same) is not."""

    def setUp(self):
        super().setUp()
        self.work = make_panel(os.path.join(self.tmp, "work"), [NRR])
        self.index = os.path.join(self.work, "genome_library", NRR)
        self.publish(self.index, dictless=True)  # (2.10.0's publish)
        self.before = hub_files(self.hub, FMT2)
        add_mnv_files(self.work)
        self.stage = os.path.join(self.tmp, "stage")

    def stage_it(self, hub=None, **kw):
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=hub or self.hub), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()), \
                _NoNetwork():
            hf.publish_index(self.index, stage_dir=self.stage, **kw)
        return out.getvalue()

    def staged(self):
        return sorted(os.listdir(os.path.join(self.stage, "indexes")))

    def test_it_packs_what_the_publish_uploads_and_uploads_nothing(self):
        mark = len(self.hub.calls)
        said = self.stage_it(dictless=True, replace_companions=True)
        calls = self.hub.calls[mark:]
        self.assertEqual(sorted(c[:3] for c in calls), sorted(
            ("hf_hub_download", FMT2, f"indexes/{n}.json")
            for n in (f"dictionaries-bundle_{VCF}", f"genotypes_{VCF}", f"indel_genotypes_{VCF}")))
        self.assertEqual(self.hub.uploads(mark), [])
        self.assertEqual(hub_files(self.hub, FMT2), self.before)
        self.assertEqual(self.staged(), sorted([
            f"{NRR}.tar.gz", f"dictionaries-bundle_{VCF}.json",
            f"dictionaries-bundle_{VCF}.tar.gz", f"indel_genotypes_{VCF}.json",
            f"indel_genotypes_{VCF}.tar.gz"]))
        at = lambda n: os.path.join(self.stage, "indexes", n)  # noqa: E731
        # the SHA-256 lines of the archives staged, sha256sum's format, by name
        with open(os.path.join(self.stage, "staged.sha256")) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines, [f"{hf._sha256_file(at(n))}  {n}" for n in sorted(
            n for n in self.staged() if n.endswith(".tar.gz"))])
        # the companions with their MNV files, their records what they hold; the index's
        # manifest names the new bundle, and the genotype store the hub's (kept)
        d = os.path.join(self.work, "Dictionaries")
        ig = _json(at(f"indel_genotypes_{VCF}.json"))
        self.assertEqual(ig["fingerprint"], hf.companion_fingerprint(
            [(f"Dictionaries/indel_genotypes_{VCF}", os.path.join(d, f"indel_genotypes_{VCF}"))]
        )["fingerprint"])
        with tarfile.open(at(f"indel_genotypes_{VCF}.tar.gz")) as tf:
            self.assertIn(f"Dictionaries/indel_genotypes_{VCF}/gt_mnv_chr22.tsv.gz", tf.getnames())
        with tarfile.open(at(f"{NRR}.tar.gz")) as tf:
            manifest = json.load(tf.extractfile("manifest.json"))
        bundle = _json(at(f"dictionaries-bundle_{VCF}.json"))
        self.assertEqual(manifest["dictionaries_bundle"]["fingerprint"], bundle["fingerprint"])
        self.assertNotEqual(bundle["fingerprint"], json.loads(
            self.before[f"indexes/dictionaries-bundle_{VCF}.json"])["fingerprint"])
        self.assertEqual(manifest["genotypes"]["fingerprint"], json.loads(
            self.before[f"indexes/genotypes_{VCF}.json"])["fingerprint"])
        # the commands, printed and not run
        title = (f"Publish dictionaries-bundle_{VCF}, indel_genotypes_{VCF}, {NRR} with "
                 "indexes/SHA256SUMS, in one commit (publish-index --stage-only)")
        cmds = hf.stage_commands(self.stage, FMT2, title)
        self.assertEqual(cmds[1], f"hf upload {FMT2} {self.stage}/indexes indexes --repo-type "
                                  f"dataset --commit-message '{title}'")
        self.assertTrue(cmds[0].startswith("python3 ") and cmds[0].endswith(
            f" {FMT2} --staged {self.stage}/staged.sha256 > {self.stage}/indexes/SHA256SUMS"))
        self.assertTrue(cmds[2].endswith(f" {FMT2} | diff - {self.stage}/indexes/SHA256SUMS"))
        for cmd in cmds:
            self.assertIn(f"  {cmd}\n", said)
        self.assertIn(f"Staged in {self.stage}; NOTHING was uploaded to {FMT2}", said)
        self.assertNotIn("Published", said)

    def test_no_token_no_hfapi_and_anonymous_reads(self):
        """With HF_TOKEN set, the reads are anonymous all the same (token False), and no
        HfApi, the hub's writer, is made; the socket layer is never reached (_NoNetwork)."""
        tokens = []

        class Hub(FakeHub):
            def hf_hub_download(hub, repo_id, filename, repo_type=None, token=None, **k):
                tokens.append(token)
                return FakeHub.hf_hub_download(hub, repo_id, filename, repo_type, token, **k)

        hub = Hub(self.hub.root)

        class NoWriter(object):
            def __init__(self, *a, **k):
                raise AssertionError("an HfApi was made: the publish could upload")

        hub.HfApi = NoWriter
        with mock.patch.dict(os.environ, {"HF_TOKEN": "hf_not_to_be_sent"}):
            self.stage_it(hub=hub, dictless=True, replace_companions=True)
        self.assertEqual(tokens, [False, False, False])
        self.assertEqual(hub_files(self.hub, FMT2), self.before)

    def test_the_folder_in_one_commit_is_what_a_download_installs(self):
        """The staged folder copied into the repo whole (what `hf upload` of it does, in one
        commit): a fresh download of the index installs the bundle with the MNV files, and
        an install of 2.10.0 gains them by the download step."""
        old = os.path.join(self.tmp, "old")
        os.makedirs(old)
        self.download(old, NRR)  # an install of 2.10.0's publish
        self.stage_it(dictless=True, replace_companions=True)
        for n in self.staged():
            shutil.copyfile(os.path.join(self.stage, "indexes", n),
                            self.remote(FMT2, f"indexes/{n}"))
        new = os.path.join(self.tmp, "new")
        os.makedirs(new)
        self.download(new, NRR)
        mnv = f"Dictionaries/indel_genotypes_{VCF}/gt_mnv_chr22.tsv.gz"
        self.assertTrue(os.path.isfile(os.path.join(new, mnv)))
        self.assertFalse(os.path.isfile(os.path.join(old, mnv)))
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", out), _NoNetwork():
            hf.download_indel_genotypes(old, VCF)
        with open(os.path.join(old, mnv), "rb") as a, open(os.path.join(new, mnv), "rb") as b:
            self.assertEqual(a.read(), b.read())

    def test_a_folder_that_holds_files_is_refused(self):
        os.makedirs(os.path.join(self.stage, "indexes"))
        _write(os.path.join(self.stage, "indexes", "leftover.tar.gz"))
        with self.assertRaises(ValueError) as caught:
            self.stage_it(dictless=True, replace_companions=True)
        self.assertIn("stage into a new, empty folder", str(caught.exception))
        self.assertEqual(os.listdir(os.path.join(self.stage, "indexes")), ["leftover.tar.gz"])

    def test_another_build_is_refused_unless_replaced_and_nothing_staged(self):
        with self.assertRaises(ValueError) as caught:
            self.stage_it(dictless=True)
        self.assertIn("Nothing was uploaded", str(caught.exception))
        self.assertEqual(self.staged(), [])
        self.assertFalse(os.path.exists(os.path.join(self.stage, "staged.sha256")))

    def test_the_command_line(self):
        """crisprme.py publish-index --stage-only <folder> passes the folder (made absolute),
        and asks for one."""
        import test_bulge_cap as tbc
        crisprme = tbc.crisprme
        self.addCleanup(setattr, crisprme, "input_args", crisprme.input_args)
        got = []
        with mock.patch.object(crisprme, "publish_index",
                               side_effect=lambda *a, **k: got.append(k["stage_dir"])):
            crisprme.input_args = ["crisprme.py", "publish-index", "--index", self.index,
                                   "--dictless", "--stage-only", "rel/stage"]
            crisprme.publish_index_cmd()
            self.assertEqual(got, [os.path.abspath("rel/stage")])
            crisprme.input_args = ["crisprme.py", "publish-index", "--index", self.index,
                                   "--stage-only"]
            with mock.patch("sys.stderr", io.StringIO()) as err, self.assertRaises(SystemExit):
                crisprme.publish_index_cmd()
            self.assertIn("--stage-only needs the folder", err.getvalue())
        self.assertEqual(got, [os.path.abspath("rel/stage")])


class TheMnvDownloadStep(_Base):
    """download --what indel-genotypes --dataset <vcf> (2.10.1, the plan of 2.10.1's MNVs,
    question 1, option A): an install of 2.10.0, whose phased indel genotype store has no MNV
    file, fetches the panel's indel_genotypes_<vcf> companion and installs it whole, where it
    is the documented upgrade: the store installed, file for file, with an MNV file for each
    chromosome. Any other build is refused with DownloadIntegrityError (nothing installed),
    as a download refuses a companion the index's manifest does not name; a repo whose store
    has no MNV yet is said; a store that has them costs no hub call. Hermetic: the fake hub,
    no network (_NoNetwork)."""

    def setUp(self):
        super().setUp()
        self.work = make_panel(os.path.join(self.tmp, "work"), [NRR])
        self.publish(os.path.join(self.work, "genome_library", NRR), dictless=True)
        self.publish(os.path.join(self.work, "genome_library", "20-NRR_3_hg38"))
        self.dest = os.path.join(self.tmp, "dest")
        os.makedirs(self.dest)
        self.download(self.dest, NRR)  # an install of 2.10.0
        self.store = os.path.join(self.dest, "Dictionaries", f"indel_genotypes_{VCF}")
        self.indel = os.path.join(self.store, "gt_indel_chr22.tsv.gz")
        with open(self.indel, "rb") as fh:
            self.indel_bytes = fh.read()

    def republish(self, work=None):
        """The panel's republish of 2.10.1 (--replace-companions): its store with the MNV
        files."""
        work = work or self.work
        add_mnv_files(work)
        self.publish(os.path.join(work, "genome_library", NRR), dictless=True,
                     replace_companions=True)

    def step(self, **kw):
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()), \
                _NoNetwork():
            hf.download_indel_genotypes(self.dest, VCF, **kw)
        return out.getvalue()

    def assert_untouched(self):
        self.assertEqual(sorted(os.listdir(self.store)),
                         ["gt_indel_chr22.tsv.gz"])
        with open(self.indel, "rb") as fh:
            self.assertEqual(fh.read(), self.indel_bytes)
        self.assertFalse([n for n in os.listdir(os.path.join(self.dest, "Dictionaries"))
                          if n.startswith((".extract_", f".indel_genotypes_{VCF}"))])
        self.assertFalse([n for n in os.listdir(self.dest) if n.startswith(".hf_stage_")])

    def test_it_installs_the_store_with_its_mnv_files(self):
        self.republish()
        mark = len(self.hub.calls)
        said = self.step()
        self.assertEqual(self.hub.fetched(mark), [f"indexes/indel_genotypes_{VCF}.tar.gz",
                                                  f"indexes/indel_genotypes_{VCF}.json"])
        self.assertEqual(sorted(os.listdir(self.store)), [
            "gt_indel_chr22.tsv.gz", "gt_mnv_chr22.tsv.gz", "gt_mnv_chr22.tsv.gz.done"])
        with open(self.indel, "rb") as fh:
            self.assertEqual(fh.read(), self.indel_bytes)
        import build_indel_genotypes as big
        reader = big.IndelGenotypeReader(self.indel)
        self.assertTrue(reader.has_mnv)
        self.assertEqual(reader.carriers_dict(120, "AA", "GG"), {"S1": "1|0"})
        marker = _json(hf.companion_marker(self.dest, f"indel_genotypes_{VCF}"))
        record = _json(self.remote(FMT2, f"indexes/indel_genotypes_{VCF}.json"))
        self.assertEqual((marker["fingerprint"], marker["repo"], marker["indexes"]),
                         (record["fingerprint"], FMT2, [NRR]))
        self.assertIn(f"Installed indel_genotypes_{VCF} from {FMT2}", said)
        # a store that has them: no hub call
        mark = len(self.hub.calls)
        said = self.step()
        self.assertEqual(self.hub.calls[mark:], [])
        self.assertIn("nothing fetched", said)
        # and the index is a complete install still: no hub call either
        mark = len(self.hub.calls)
        self.download(self.dest, NRR)
        self.assertEqual(self.hub.calls[mark:], [])

    def test_a_repo_whose_store_has_no_mnv_says_so(self):
        with self.assertRaises(ValueError) as caught:
            self.step()
        self.assertIn("holds no MNV file", str(caught.exception))
        self.assertIn("Nothing was installed", str(caught.exception))
        self.assert_untouched()

    def test_another_build_of_the_store_is_refused(self):
        other = make_panel(os.path.join(self.tmp, "other"), [NRR])
        _gz(os.path.join(other, "Dictionaries", f"indel_genotypes_{VCF}", "gt_indel_chr22.tsv.gz"),
            "S2\t1|0\n", 1)
        self.republish(other)
        with self.assertRaises(hf.DownloadIntegrityError) as caught:
            self.step()
        self.assertIn("another build of the phased indel genotypes", str(caught.exception))
        self.assertIn("gt_indel_chr22.tsv.gz", str(caught.exception))
        self.assertIn("Nothing of it was installed", str(caught.exception))
        self.assert_untouched()

    def test_an_archive_that_is_not_what_its_record_says_is_refused(self):
        self.republish()
        record = self.remote(FMT2, f"indexes/indel_genotypes_{VCF}.json")
        rec = _json(record)
        rec["fingerprint"] = "0" * 64
        with open(record, "w") as fh:
            json.dump(rec, fh)
        with self.assertRaises(hf.DownloadIntegrityError) as caught:
            self.step()
        self.assertIn("other files than its record says", str(caught.exception))
        self.assert_untouched()

    def test_a_chromosome_left_without_its_mnv_file_is_refused(self):
        self.republish()
        # the install has a chromosome the repo's store gives no MNV file
        _gz(os.path.join(self.store, "gt_indel_chr21.tsv.gz"), "S1\t1|0\n", 1)
        _gz(os.path.join(self.work, "Dictionaries", f"indel_genotypes_{VCF}",
                         "gt_indel_chr21.tsv.gz"), "S1\t1|0\n", 1)
        self.publish(os.path.join(self.work, "genome_library", NRR), dictless=True,
                     replace_companions=True)
        with self.assertRaises(hf.DownloadIntegrityError) as caught:
            self.step()
        self.assertIn("no whole MNV file for 1 of the 2 chromosome(s)", str(caught.exception))
        self.assertFalse(os.path.exists(os.path.join(self.store, "gt_mnv_chr22.tsv.gz")))

    def test_an_mnv_file_that_cannot_be_read_is_fetched_whole_again(self):
        """A store whose MNV file is corrupt (its .done marker left as it was) is not taken
        for one that has its MNVs: the step installs the repo's whole again. So too one of 0
        bytes, and a whole gzip file of no line where its marker, as new as the file, counts
        one record (build_indel_genotypes.read_mnv_lines): at 6b5ab4d3 the step took the
        latter for whole, and fetched nothing."""
        self.republish()
        self.step()
        mnv = os.path.join(self.store, "gt_mnv_chr22.tsv.gz")
        with open(mnv, "rb") as fh:
            whole = fh.read()
        for what, data in (("cut short", whole[:10]), ("0 bytes", b""),
                           ("not the file its marker counts", gzip.compress(b""))):
            with open(mnv, "wb") as fh:
                fh.write(data)
            os.utime(mnv, (1, 1))
            said = self.step()
            self.assertIn(f"Installed indel_genotypes_{VCF}", said, what)
            with open(mnv, "rb") as fh:
                self.assertEqual(fh.read(), whole, what)

    def test_no_store_or_an_unfinished_one_is_an_error(self):
        os.makedirs(os.path.join(self.store, ".building"))
        with self.assertRaises(ValueError) as caught:
            self.step()
        self.assertIn("is unfinished", str(caught.exception))
        shutil.rmtree(self.store)
        with self.assertRaises(ValueError) as caught:
            self.step()
        self.assertIn(f"No phased indel genotype store of {VCF} is installed", str(caught.exception))
        self.assertEqual(self.hub.fetched(len(self.hub.calls)), [])

    def test_the_command_line_and_its_help_name_it(self):
        """crisprme.py download --what indel-genotypes --dataset <vcf> runs the step (its
        refusals an Error, exit 1), --dataset required; download --help names it."""
        import test_bulge_cap as tbc
        crisprme = tbc.crisprme
        self.addCleanup(setattr, crisprme, "input_args", crisprme.input_args)
        self.republish()
        crisprme.input_args = ["crisprme.py", "download", "--what", "indel-genotypes",
                               "--dataset", VCF, "--path", self.dest]
        out = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", out), _NoNetwork():
            crisprme.download_data()
        self.assertIn(f"Phased indel genotype store: {self.store}", out.getvalue())
        self.assertTrue(os.path.isfile(os.path.join(self.store, "gt_mnv_chr22.tsv.gz")))
        crisprme.input_args = ["crisprme.py", "download", "--what", "indel-genotypes",
                               "--path", self.dest]
        with mock.patch("sys.stderr", io.StringIO()) as err, self.assertRaises(SystemExit):
            crisprme.download_data()
        self.assertIn("--what indel-genotypes requires --dataset", err.getvalue())
        # a refusal (another build) is the command's Error, as a download's
        shutil.rmtree(self.store)
        os.makedirs(self.store)
        _gz(os.path.join(self.store, "gt_indel_chr22.tsv.gz"), "S9\t1|0\n", 1)
        crisprme.input_args = ["crisprme.py", "download", "--what", "indel-genotypes",
                               "--dataset", VCF, "--path", self.dest]
        with mock.patch.object(hf, "_require_hf", return_value=self.hub), \
                mock.patch("sys.stdout", io.StringIO()), \
                mock.patch("sys.stderr", io.StringIO()) as err, \
                self.assertRaises(SystemExit) as code, _NoNetwork():
            crisprme.download_data()
        self.assertEqual(code.exception.code, 1)
        self.assertIn("another build of the phased indel genotypes", err.getvalue())
        with mock.patch("sys.stderr", io.StringIO()) as err, self.assertRaises(SystemExit):
            crisprme.print_help_download()
        self.assertIn("crisprme.py download --what indel-genotypes --dataset hg38_HPRC",
                      " ".join(err.getvalue().split()))


class TheDefaultReferenceIndex(unittest.TestCase):
    """download --what all fetches 20-NRR_3_hg38, and 2.8.1's NRR_3_hg38 where it cannot
    (before the rebuild is published). That download then said "This release's build of it is
    20-NRR_3_hg38 ...: download that one instead", the index that had just failed: the
    fallback names no successor (its own note says why it is fetched); a download of
    NRR_3_hg38 asked for by name still does."""

    def fetch(self, fails):
        import test_bulge_cap as tbc
        crisprme = tbc.crisprme
        calls = []

        def fake(component, workdir, **kw):
            calls.append((kw["index_name"], kw.get("_note_successor", True)))
            if kw["index_name"] in fails:
                raise ValueError(f"No index '{kw['index_name']}' was found")
            return os.path.join(workdir, "genome_library", kw["index_name"])

        out = io.StringIO()
        with mock.patch.object(crisprme, "download_component", side_effect=fake), \
                mock.patch("sys.stdout", out), mock.patch("sys.stderr", io.StringIO()):
            crisprme._download_default_reference_index("/w", None, "hg38")
        return calls, out.getvalue()

    def test_the_fallback_names_no_successor(self):
        calls, said = self.fetch({"20-NRR_3_hg38"})
        self.assertEqual(calls, [("20-NRR_3_hg38", True), ("NRR_3_hg38", False)])
        self.assertIn("fetching NRR_3_hg38, CRISPRme 2.8.1's, in its place", said)
        self.assertIn("Downloaded reference index NRR_3_hg38", said)

    def test_the_default_alone_where_it_is_published(self):
        calls, said = self.fetch(set())
        self.assertEqual(calls, [("20-NRR_3_hg38", True)])
        self.assertIn("Downloaded reference index 20-NRR_3_hg38", said)


if __name__ == "__main__":
    unittest.main()

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
FMT2 = getattr(hf, "FORMAT2_HF_REPO", "lucapinello/crisprme-indexes-format2")
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
        old = make_panel(os.path.join(self.tmp, "old"), [name], fmt=1)
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
        # (the note is read from the index set, offline: the format-2 repo is not asked)
        self.assertNotIn(FMT2, {c[1] for c in self.hub.calls if c[0] != "upload_file"})

    def test_the_index_set_replaces_every_archive_of_2_8_1(self):
        self.assertEqual(hf.successor_index("NRG_3_hg38+hg38_HPRC"), "20-NRR_3_hg38+hg38_HPRC")
        self.assertEqual(hf.successor_index("NNNNNNNN_3_hg38"), "20-NNNNNNNN_3_hg38")
        self.assertIsNone(hf.successor_index(NRR))


if __name__ == "__main__":
    unittest.main()

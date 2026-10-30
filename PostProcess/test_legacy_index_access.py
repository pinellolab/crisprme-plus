#!/usr/bin/env python3
"""Dedicated-PAM ("legacy") indexes stay downloadable exactly as before the repo split.

2.8.1 sends PAM-less indexes with a flank wider than 3 nt to a second dataset repo.
The dedicated-PAM indexes (NRG_3 / NGG_3 / NRR_3 and their variant panels) and the
3-nt PAM-less NNN_3 search several times faster than a wide PAM-less index pruned at
search time, so they must keep being listed from, and downloaded from, the default
repo: the archive, its genotype companion and the shared reference index alike, with
no lookup in the second repo, even when that repo holds an archive of the same name.

test_index_repo_routing.py pins the name -> repo function. This file pins what the
download and the catalog actually DO with it, by recording every hub call.

Hermetic: huggingface_hub is replaced through crisprme_hf._require_hf (every hub
access in the module goes through it) by a fake that serves small real archives per
repo and records each call. No network, no token.
"""

import ast
import io
import json
import os
import pathlib
import shutil
import tarfile
import tempfile
import unittest
from typing import List, Optional
from unittest import mock

import crisprme_hf as hf

DEFAULT = hf.DEFAULT_HF_REPO
GENERAL = hf.GENERAL_INDEX_HF_REPO

# The index archives in lucapinello/crisprme-data/indexes/ as listed by
# `scripts/hf_sha256sums.py lucapinello/crisprme-data` on 2026-10-30.
LEGACY = (
    "NGG_3_hg38", "NGG_3_hg38+hg38_HGDP1kGP",
    "NNN_3_hg38", "NNN_3_hg38+hg38_HGDP1kGP",
    "NRG_3_hg38", "NRG_3_hg38+hg38_HGDP1kGP", "NRG_3_hg38+hg38_HPRC", "NRG_3_hg38+hg38_mega",
    "NRR_3_hg38", "NRR_3_hg38+hg38_HGDP1kGP", "NRR_3_hg38+hg38_HPRC", "NRR_3_hg38+hg38_mega",
)
GENOTYPED_PANELS = ("hg38_HGDP1kGP", "hg38_HPRC")  # mega is sites-only: no companion
WIDE = ("NNNNNNNN_3_hg38", "NNNNNNNN_3_hg38+hg38_HGDP1kGP")


class _NotFound(Exception):
    """What the hub raises for a repo that does not exist."""


class _Entry:
    def __init__(self, path, size):
        self.path, self.size = path, size


class _FakeHub:
    """Serves ``{repo: {path_in_repo: local_file}}`` and records every call."""

    def __init__(self, repos):
        self.repos = repos
        self.calls = []  # (function, repo, detail, token)
        hub = self

        class HfApi:
            def __init__(self, *a, **k):
                pass

            def list_repo_files(self, repo_id, repo_type=None, token=None, **k):
                hub.calls.append(("list_repo_files", repo_id, None, token))
                return sorted(hub._repo(repo_id))

            def file_exists(self, repo_id, filename, repo_type=None, token=None, **k):
                hub.calls.append(("file_exists", repo_id, filename, token))
                return filename in hub.repos.get(repo_id, {})

            def list_repo_tree(self, repo_id, repo_type=None, path_in_repo=None,
                               recursive=False, token=None, **k):
                hub.calls.append(("list_repo_tree", repo_id, path_in_repo, token))
                prefix = path_in_repo.strip("/") + "/"
                return [_Entry(p, os.path.getsize(f))
                        for p, f in sorted(hub._repo(repo_id).items())
                        if p.startswith(prefix) and "/" not in p[len(prefix):]]

        self.HfApi = HfApi

    def _repo(self, repo_id):
        if repo_id not in self.repos:
            raise _NotFound(f"Repository Not Found: {repo_id}")
        return self.repos[repo_id]

    def snapshot_download(self, repo_id, repo_type=None, allow_patterns=None,
                          local_dir=None, token=None, **k):
        self.calls.append(("snapshot_download", repo_id, tuple(allow_patterns), token))
        for path, src in self._repo(repo_id).items():
            if path in allow_patterns:  # the download path only asks for exact paths
                dst = os.path.join(local_dir, path)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copyfile(src, dst)
        return local_dir

    def get_paths_info(self, repo_id, paths, repo_type=None, token=None, **k):
        self.calls.append(("get_paths_info", repo_id, None, token))
        files = self._repo(repo_id)
        return [_Entry(p, os.path.getsize(files[p])) for p in paths if p in files]

    def repos_contacted(self):
        return {c[1] for c in self.calls}


def _add(tf, name, data=b""):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tf.addfile(info, io.BytesIO(data))


class _Fixture(unittest.TestCase):
    """Builds per-repo archives stamped with the repo that serves them."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        saved = os.environ.pop("CRISPRME_HF_REPO", None)
        self.addCleanup(self._restore_env, saved)

    @staticmethod
    def _restore_env(saved):
        os.environ.pop("CRISPRME_HF_REPO", None)
        if saved is not None:
            os.environ["CRISPRME_HF_REPO"] = saved

    def _archive(self, repo, name, pad=0):
        path = os.path.join(self.tmp, "remote", repo.replace("/", "__"), f"{name}.tar.gz")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        stamp = repo.encode()
        with tarfile.open(path, "w:gz") as tf:
            _add(tf, f"{name}/served_by", stamp)
            _add(tf, f"{name}/TSTgenome.bin", b"\0" * (8 + pad))
            vcf = name.partition("+")[2]
            if vcf:
                _add(tf, f"{name}_INDELS/served_by", stamp)
                _add(tf, f"Dictionaries/registry_{vcf}/served_by", stamp)
            _add(tf, "manifest.json", json.dumps({"created_at": "test"}).encode())
        return path

    def _genotypes(self, repo, vcf):
        path = os.path.join(self.tmp, "remote", repo.replace("/", "__"), f"genotypes_{vcf}.tar.gz")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with tarfile.open(path, "w:gz") as tf:
            _add(tf, f"Dictionaries/genotypes_{vcf}/served_by", repo.encode())
        return path

    def _repo(self, repo, names, panels, pad=0):
        files = {f"indexes/{n}.tar.gz": self._archive(repo, n, pad) for n in names}
        files.update({f"indexes/genotypes_{v}.tar.gz": self._genotypes(repo, v) for v in panels})
        return files

    def _download(self, hub, name, **kw):
        workdir = tempfile.mkdtemp(dir=self.tmp)
        with mock.patch.object(hf, "_require_hf", return_value=hub), \
             mock.patch("sys.stdout", new_callable=io.StringIO), \
             mock.patch("sys.stderr", new_callable=io.StringIO):
            returned = hf.download_component("index", workdir, index_name=name, **kw)
        return workdir, returned

    @staticmethod
    def _stamp(workdir, rel):
        with open(os.path.join(workdir, rel, "served_by")) as fh:
            return fh.read()


class EveryPublishedLegacyNameRoutesToTheDefaultRepo(unittest.TestCase):
    def setUp(self):
        saved = os.environ.pop("CRISPRME_HF_REPO", None)
        self.addCleanup(_Fixture._restore_env, saved)

    def test_each_published_archive(self):
        for name in LEGACY:
            self.assertEqual(hf.resolve_index_repo(name), DEFAULT, name)

    def test_the_dictless_spelling_of_each(self):
        for name in LEGACY:
            base, plus, vcf = name.partition("+")
            spelled = base + "-dictless" + plus + vcf
            self.assertEqual(hf.resolve_index_repo(spelled), DEFAULT, spelled)

    def test_published_names_are_already_canonical(self):
        """The install dir is genome_library/<canonical name>; for these it is the name."""
        for name in LEGACY:
            self.assertEqual(hf.canonical_index_name(name), name)


class LegacyDownloadNeverTouchesTheGeneralIndexRepo(_Fixture):
    """Archive, genotype companion and reference index all come from the default repo,
    with no routing probe, even when the general-index repo carries the same names."""

    def test_every_legacy_index(self):
        repos = {
            DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS),
            # decoys: every legacy name and both companions, stamped GENERAL
            GENERAL: self._repo(GENERAL, LEGACY + WIDE, GENOTYPED_PANELS, pad=4096),
        }
        for name in LEGACY:
            with self.subTest(name=name):
                hub = _FakeHub(repos)
                workdir, returned = self._download(hub, name)
                self.assertEqual(hub.repos_contacted(), {DEFAULT}, hub.calls)
                self.assertFalse([c for c in hub.calls if c[0] == "list_repo_files"],
                                 "a legacy name must not trigger the routing probe")
                self.assertEqual({c[3] for c in hub.calls}, {None}, "anonymous download")
                self.assertEqual(returned, os.path.join(workdir, "genome_library", name))
                self.assertEqual(self._stamp(workdir, f"genome_library/{name}"), DEFAULT)
                base, _, vcf = name.partition("+")
                if not vcf:
                    continue
                self.assertEqual(self._stamp(workdir, f"genome_library/{name}_INDELS"), DEFAULT)
                self.assertEqual(self._stamp(workdir, f"genome_library/{base}"), DEFAULT,
                                 "shared reference index from the wrong repo")
                fetched = [c[2] for c in hub.calls if c[0] == "snapshot_download"]
                self.assertIn((f"indexes/genotypes_{vcf}.tar.gz",), fetched)
                if vcf in GENOTYPED_PANELS:
                    self.assertEqual(
                        self._stamp(workdir, f"Dictionaries/genotypes_{vcf}"), DEFAULT)

    def test_the_general_index_repo_being_absent_changes_nothing(self):
        hub = _FakeHub({DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS)})
        workdir, _ = self._download(hub, "NRG_3_hg38+hg38_HGDP1kGP")
        self.assertEqual(hub.repos_contacted(), {DEFAULT})
        self.assertEqual(self._stamp(workdir, "genome_library/NRG_3_hg38"), DEFAULT)


class WideFlankDownloadKeepsItsPiecesTogether(_Fixture):
    """For contrast: a K=8 index takes all three pieces from its own repo, and falls
    back to the default repo only when its own repo lacks the archive."""

    def test_all_pieces_from_the_general_index_repo(self):
        hub = _FakeHub({
            DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS),
            GENERAL: self._repo(GENERAL, WIDE, ("hg38_HGDP1kGP",)),
        })
        workdir, _ = self._download(hub, "NNNNNNNN_3_hg38+hg38_HGDP1kGP")
        fetched = {c[1] for c in hub.calls if c[0] == "snapshot_download"}
        self.assertEqual(fetched, {GENERAL})
        self.assertEqual(self._stamp(workdir, "genome_library/NNNNNNNN_3_hg38"), GENERAL)
        self.assertEqual(self._stamp(workdir, "Dictionaries/genotypes_hg38_HGDP1kGP"), GENERAL)

    def test_falls_back_to_the_default_repo_when_its_own_lacks_it(self):
        hub = _FakeHub({DEFAULT: self._repo(DEFAULT, LEGACY + WIDE, GENOTYPED_PANELS)})
        workdir, _ = self._download(hub, "NNNNNNNN_3_hg38")
        self.assertEqual({c[1] for c in hub.calls if c[0] == "snapshot_download"}, {DEFAULT})
        self.assertEqual(self._stamp(workdir, "genome_library/NNNNNNNN_3_hg38"), DEFAULT)

    def test_absent_from_both_repos_names_both(self):
        hub = _FakeHub({DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS)})
        with self.assertRaises(ValueError) as caught:
            self._download(hub, "NNNNNNNN_3_hg38")
        message = str(caught.exception)
        self.assertIn(GENERAL, message)
        self.assertIn(DEFAULT, message)
        self.assertIn(f"--hf-repo {GENERAL}", message)

    def test_an_installed_wide_flank_index_costs_no_hub_call(self):
        hub = _FakeHub({GENERAL: self._repo(GENERAL, WIDE, ("hg38_HGDP1kGP",))})
        workdir = tempfile.mkdtemp(dir=self.tmp)
        os.makedirs(os.path.join(workdir, "genome_library", "NNNNNNNN_3_hg38"))
        with mock.patch.object(hf, "_require_hf", return_value=hub), \
             mock.patch("sys.stdout", new_callable=io.StringIO):
            hf.download_component("index", workdir, index_name="NNNNNNNN_3_hg38")
        self.assertEqual(hub.calls, [])


class TheFallbackProbeLooksForWhatTheDownloadFetches(_Fixture):
    """The probe that redirects a routed index to the default repo must look up the
    exact archive(s) the download then asks for, with the same token.

    It used to list the whole repo untokened, for the archive under the CANONICAL
    name -- so a wide-flank index the general-index repo holds only under its
    `-dictless` name (the spelling the merged catalog offers), or a repo needing a
    token, read as absent and was redirected to a default repo without it."""

    REQUESTED = "NNNNNNNN_3_hg38+hg38_HGDP1kGP"
    DICTLESS = "NNNNNNNN_3_hg38-dictless+hg38_HGDP1kGP"
    TOKEN = "fake-token-for-tests"

    def _general(self, *names):
        return _FakeHub({
            DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS),
            GENERAL: self._repo(GENERAL, names, ("hg38_HGDP1kGP",)),
        })

    def test_an_archive_only_under_its_dictless_name_is_not_redirected(self):
        hub = self._general(self.DICTLESS, "NNNNNNNN_3_hg38")
        workdir, returned = self._download(hub, self.REQUESTED)
        self.assertEqual({c[1] for c in hub.calls if c[0] == "snapshot_download"},
                         {GENERAL}, hub.calls)
        # -dictless is a publish marker: it installs under the canonical name
        self.assertEqual(returned, os.path.join(workdir, "genome_library", self.REQUESTED))
        self.assertEqual(self._stamp(workdir, f"genome_library/{self.REQUESTED}"), GENERAL)

    def test_the_probe_uses_the_download_token_and_exact_paths(self):
        hub = self._general(self.REQUESTED, "NNNNNNNN_3_hg38")
        self._download(hub, self.REQUESTED, token=self.TOKEN)
        self.assertFalse([c for c in hub.calls if c[0] == "list_repo_files"],
                         "the probe still lists the whole repo")
        self.assertEqual({c[3] for c in hub.calls}, {self.TOKEN}, hub.calls)
        probed = [c[2] for c in hub.calls if c[0] == "file_exists"]
        self.assertIn(f"indexes/{self.REQUESTED}.tar.gz", probed)
        self.assertTrue(all(p.startswith("indexes/") and p.endswith(".tar.gz")
                            for p in probed), probed)

    def test_an_archive_in_neither_repo_still_falls_back_and_names_both(self):
        hub = self._general("NNNNNNNN_3_hg38")
        with self.assertRaises(ValueError) as caught:
            self._download(hub, self.REQUESTED)
        message = str(caught.exception)
        for expected in (GENERAL, DEFAULT, f"--hf-repo {GENERAL}"):
            self.assertIn(expected, message)


class AnExplicitRepoStillWins(_Fixture):
    def _hub(self, repo):
        return _FakeHub({
            DEFAULT: self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS),
            repo: self._repo(repo, LEGACY, GENOTYPED_PANELS),
        })

    def test_hf_repo_flag(self):
        hub = self._hub("me/mirror")
        workdir, _ = self._download(hub, "NRG_3_hg38+hg38_HGDP1kGP", repo="me/mirror")
        self.assertEqual(hub.repos_contacted(), {"me/mirror"})
        self.assertEqual(self._stamp(workdir, "genome_library/NRG_3_hg38"), "me/mirror")

    def test_env_var(self):
        os.environ["CRISPRME_HF_REPO"] = "me/env-repo"
        hub = self._hub("me/env-repo")
        workdir, _ = self._download(hub, "NRG_3_hg38+hg38_HGDP1kGP")
        self.assertEqual(hub.repos_contacted(), {"me/env-repo"})
        self.assertEqual(self._stamp(workdir, "Dictionaries/genotypes_hg38_HGDP1kGP"),
                         "me/env-repo")


class CatalogKeepsEveryLegacyIndex(_Fixture):
    """list_available_downloads('index') is what the web Settings dropdown shows."""

    def _list(self, hub, repo=None):
        with mock.patch.object(hf, "_require_hf", return_value=hub):
            return hf.list_available_downloads("index", repo)

    def setUp(self):
        super().setUp()
        self.default_files = self._repo(DEFAULT, LEGACY, GENOTYPED_PANELS)

    def _default_only(self):
        return _FakeHub({DEFAULT: self.default_files})

    def test_both_repos_listed_and_every_legacy_row_kept_verbatim(self):
        expected = {r["name"]: r for r in self._list(self._default_only(), DEFAULT)}
        self.assertEqual(sorted(expected), sorted(LEGACY))
        hub = _FakeHub({
            DEFAULT: self.default_files,
            GENERAL: self._repo(GENERAL, LEGACY + WIDE, GENOTYPED_PANELS, pad=4096),
        })
        rows = {r["name"]: r for r in self._list(hub)}
        self.assertEqual({c[1] for c in hub.calls}, {DEFAULT, GENERAL})
        self.assertEqual(sorted(rows), sorted(LEGACY + WIDE))
        for name in LEGACY:  # a same-named general-index row never shadows the default's
            self.assertEqual(rows[name], expected[name], name)

    def test_companions_are_not_offered_as_indexes(self):
        rows = self._list(self._default_only())
        self.assertFalse([r for r in rows if r["name"].startswith("genotypes_")])
        flags = {r["name"]: r["has_genotypes"] for r in rows}
        self.assertTrue(flags["NRG_3_hg38+hg38_HGDP1kGP"])
        self.assertTrue(flags["NRR_3_hg38+hg38_HPRC"])
        self.assertFalse(flags["NRG_3_hg38+hg38_mega"])
        self.assertFalse(flags["NRG_3_hg38"])

    def test_an_unreachable_general_index_repo_drops_no_legacy_row(self):
        alone = self._list(self._default_only(), DEFAULT)
        self.assertEqual(self._list(self._default_only()), alone)

    def test_an_explicit_repo_lists_only_that_repo(self):
        hub = _FakeHub({"me/mirror": self._repo("me/mirror", LEGACY, GENOTYPED_PANELS)})
        rows = self._list(hub, "me/mirror")
        self.assertEqual(hub.repos_contacted(), {"me/mirror"})
        self.assertEqual(sorted(r["name"] for r in rows), sorted(LEGACY))


class CatalogRowsComeFromTheCopyDownloadFetches(CatalogKeepsEveryLegacyIndex):
    """Each row describes the archive `download` would fetch for that name, so the
    size the disk preflight uses and has_genotypes are that copy's."""

    def _rows(self, repos):
        return {r["name"]: r for r in self._list(_FakeHub(repos))}

    def test_a_wide_flank_row_comes_from_the_general_index_repo(self):
        general_files = self._repo(GENERAL, WIDE, ("hg38_HGDP1kGP",), pad=4096)
        own = {r["name"]: r for r in self._list(_FakeHub({GENERAL: general_files}), GENERAL)}
        stale = self._repo(DEFAULT, LEGACY + WIDE, GENOTYPED_PANELS)
        stale_rows = {r["name"]: r for r in self._list(_FakeHub({DEFAULT: stale}), DEFAULT)}
        for name in WIDE:  # precondition: the two copies are distinguishable
            self.assertNotEqual(own[name]["size"], stale_rows[name]["size"], name)
        rows = self._rows({DEFAULT: stale, GENERAL: general_files})
        for name in WIDE:
            self.assertEqual(rows[name], own[name], name)

    def test_a_wide_flank_row_only_in_the_default_repo_is_offered(self):
        # download falls back to the default repo, so the catalog does too
        rows = self._rows({DEFAULT: self._repo(DEFAULT, LEGACY + WIDE, GENOTYPED_PANELS)})
        self.assertEqual(sorted(rows), sorted(LEGACY + WIDE))

    def test_a_legacy_name_only_in_the_general_index_repo_is_not_offered(self):
        # download looks for a legacy name in the default repo only
        rows = self._rows({
            DEFAULT: self._repo(DEFAULT, ("NRG_3_hg38",), ()),
            GENERAL: self._repo(GENERAL, ("NGG_3_hg38",) + WIDE, ()),
        })
        self.assertEqual(sorted(rows), sorted(("NRG_3_hg38",) + WIDE))

    def test_an_unreadable_default_repo_offers_nothing(self):
        # the offline signal 2.8.0 gave, not a catalog of only the slower indexes
        hub = _FakeHub({GENERAL: self._repo(GENERAL, WIDE, ("hg38_HGDP1kGP",))})
        self.assertEqual(self._list(hub), [])
        self.assertIn(DEFAULT, hub.repos_contacted())


class WebSettingsOffersAndRoutesLegacyIndexes(unittest.TestCase):
    """The Settings page lists the merged catalog unfiltered and downloads by name.

    The callbacks need Dash; their bodies are compiled from the source with the Dash
    pieces stubbed, the way test_index_repo_routing reads crisprme_hf with ast.
    """

    PAGE = pathlib.Path(__file__).resolve().parent.parent / "pages" / "settings_page.py"
    WANT = ("_fmt_size", "_catalog_size", "_hf_options", "add_index_hf",
            "launch_settings_job", "_shlex_quote")

    def setUp(self):
        self.src = self.PAGE.read_text()
        self.tree = ast.parse(self.src)
        self.defs = {n.name: n for n in self.tree.body if isinstance(n, ast.FunctionDef)}
        missing = [n for n in self.WANT + ("_hf_catalog",) if n not in self.defs]
        self.assertFalse(missing, f"settings_page.py no longer defines {missing}")

    def _load(self, catalog, workdir):
        launched = []
        nodes = []
        for name in self.WANT:
            node = self.defs[name]
            node.decorator_list = []  # the @app.callback registration is not under test
            nodes.append(node)
        env = {
            "os": os, "List": List, "Optional": Optional, "ONLINE": False,
            "PreventUpdate": RuntimeError, "no_update": "<no_update>",
            "current_working_directory": workdir,
            "_hf_catalog": lambda component: catalog if component == "index" else [],
            "_preflight_disk": lambda need: None,
            "_launch_job_cmd": lambda cmd, stage: launched.append(cmd) or "job",
            "_start": lambda job: ("started", job),
        }
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(self.PAGE), "exec"), env)
        return env, launched

    def test_every_catalog_row_is_offered(self):
        catalog = [{"name": n, "size": 1} for n in LEGACY + WIDE]
        env, _ = self._load(catalog, tempfile.gettempdir())
        values = [o["value"] for o in env["_hf_options"]("index", [])]
        self.assertEqual(values, list(LEGACY + WIDE))

    def test_download_is_by_name_with_no_repo_flag(self):
        workdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, workdir, True)
        env, launched = self._load([{"name": n, "size": 1} for n in LEGACY], workdir)
        for name in LEGACY:
            launched.clear()
            env["add_index_hf"](1, name)
            self.assertEqual(len(launched), 1, name)
            self.assertIn(f"download --what index --index-name {env['_shlex_quote'](name)}",
                          launched[0])
            self.assertNotIn("--hf-repo", launched[0])

    def test_the_catalog_is_the_merged_one(self):
        """_hf_catalog must ask for list_available_downloads(component) with no repo."""
        calls = [n for n in ast.walk(self.defs["_hf_catalog"])
                 if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "list_available_downloads"]
        self.assertEqual(len(calls), 1)
        self.assertEqual((len(calls[0].args), calls[0].keywords), (1, []))


if __name__ == "__main__":
    unittest.main(verbosity=2)

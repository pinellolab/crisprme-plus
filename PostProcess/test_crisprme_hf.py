"""Unit tests for the pure (network-free) logic of crisprme_hf.

These exercise repo/token resolution, gzip decompression, the index tar
round-trip, and argument validation without touching HuggingFace. The
network-facing calls (snapshot_download / upload_file) are covered separately by
an end-to-end test against a populated repo once one exists.
"""

import contextlib
import gzip
import importlib.util
import os
import tarfile
import tempfile
import unittest
from unittest import mock

import crisprme_hf as hf

# The network-free CI deliberately does NOT install huggingface_hub (see
# unit-tests.yml). Patching "huggingface_hub.get_token" is only possible when the
# module is importable; when it is absent, resolve_token already yields None via
# its own ImportError fallback, so no patch is needed.
_HAS_HF = importlib.util.find_spec("huggingface_hub") is not None


@contextlib.contextmanager
def _no_cached_hf_token():
    """Neutralize a cached ``huggingface-cli login`` token for the no-token path.

    A no-op when huggingface_hub is absent (resolve_token returns None on its own);
    patches get_token to None when it is present so a logged-in dev machine can't
    leak a real token into these tests.
    """
    if _HAS_HF:
        with mock.patch("huggingface_hub.get_token", return_value=None):
            yield
    else:
        yield


class TestResolveRepo(unittest.TestCase):
    def setUp(self):
        self._saved = {k: os.environ.get(k) for k in ("CRISPRME_HF_REPO", "HF_TOKEN",
                                                       "HUGGING_FACE_HUB_TOKEN")}
        for k in self._saved:
            os.environ.pop(k, None)

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_cli_wins(self):
        os.environ["CRISPRME_HF_REPO"] = "env/repo"
        self.assertEqual(hf.resolve_repo("cli/repo"), "cli/repo")

    def test_env_over_default(self):
        os.environ["CRISPRME_HF_REPO"] = "env/repo"
        self.assertEqual(hf.resolve_repo(None), "env/repo")

    def test_default(self):
        self.assertEqual(hf.resolve_repo(None), hf.DEFAULT_HF_REPO)

    def test_token_precedence(self):
        os.environ["HF_TOKEN"] = "tok_env"
        self.assertEqual(hf.resolve_token("tok_cli"), "tok_cli")
        self.assertEqual(hf.resolve_token(None), "tok_env")
        os.environ.pop("HF_TOKEN")
        os.environ["HUGGING_FACE_HUB_TOKEN"] = "tok_alt"
        self.assertEqual(hf.resolve_token(None), "tok_alt")

    def test_token_none(self):
        # neutralize a cached `huggingface-cli login` token so the no-token path
        # is exercised even on a machine that is logged into HuggingFace
        with _no_cached_hf_token():
            self.assertIsNone(hf.resolve_token(None))


class TestDecompress(unittest.TestCase):
    def test_gunzip(self):
        import tempfile

        d = tempfile.mkdtemp()
        raw = os.path.join(d, "chr1.fa")
        with open(raw, "w") as f:
            f.write(">chr1\nACGT\n")
        gzp = raw + ".gz"
        with open(raw, "rb") as src, gzip.open(gzp, "wb") as dst:
            dst.write(src.read())
        os.remove(raw)
        out = hf.decompress_gz(gzp)
        self.assertEqual(out, raw)
        self.assertTrue(os.path.exists(raw))
        self.assertFalse(os.path.exists(gzp))  # source removed by default
        with open(raw) as f:
            self.assertIn("ACGT", f.read())

    def test_noop_on_plain(self):
        self.assertEqual(hf.decompress_gz("/x/y/chr1.fa"), "/x/y/chr1.fa")


class TestComponentValidation(unittest.TestCase):
    def test_unknown_component(self):
        with self.assertRaises(ValueError):
            hf.download_component("nonsense", "/tmp")

    def test_vcf_requires_dataset(self):
        with self.assertRaises(ValueError):
            hf.download_component("vcf", "/tmp", dataset=None)

    def test_index_requires_name(self):
        with self.assertRaises(ValueError):
            hf.download_component("index", "/tmp", index_name=None)

    def test_component_maps_are_consistent(self):
        self.assertEqual(
            set(hf._COMPONENT_PREFIXES), set(hf._COMPONENT_LOCALDIR)
        )


class TestPublishValidation(unittest.TestCase):
    def test_missing_dir(self):
        with self.assertRaises(ValueError):
            hf.publish_index("/no/such/index/dir", token="x")

    def test_missing_token(self):
        import tempfile

        d = tempfile.mkdtemp()
        idx = os.path.join(d, "NGG_2_hg38")
        os.makedirs(idx)
        open(os.path.join(idx, "chr1.bin"), "wb").close()
        # no token in env, none passed, and no cached login -> must raise before
        # any network call (patch get_token so a logged-in machine can't cause a
        # real upload during the test)
        for k in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN"):
            os.environ.pop(k, None)
        with _no_cached_hf_token():
            with self.assertRaises(ValueError):
                hf.publish_index(idx, token=None)


class TestIndexTarRoundTrip(unittest.TestCase):
    """The tar layout must unpack to genome_library/<name>/ (arcname=name)."""

    def test_arcname_layout(self):
        import tempfile

        d = tempfile.mkdtemp()
        idx = os.path.join(d, "genome_library", "NGG_2_hg38")
        os.makedirs(idx)
        with open(os.path.join(idx, "chr1.bin"), "wb") as f:
            f.write(b"\x00\x01")
        tarball = os.path.join(d, "NGG_2_hg38.tar.gz")
        with tarfile.open(tarball, "w:gz") as tf:
            tf.add(idx, arcname="NGG_2_hg38")
        out = os.path.join(d, "unpack")
        os.makedirs(out)
        with tarfile.open(tarball) as tf:
            names = tf.getnames()
            tf.extractall(out)
        self.assertIn("NGG_2_hg38", names[0])
        self.assertTrue(os.path.exists(os.path.join(out, "NGG_2_hg38", "chr1.bin")))


class TestDownloadIntegrity(unittest.TestCase):
    """Post-download verification: truncation/corruption must be caught (the
    class of bug where a partial download installed silently and only failed
    later at search time with a cryptic htslib error)."""

    def _tmp(self):
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(d, ignore_errors=True))
        return d

    def test_valid_bgzf_passes(self):
        d = self._tmp()
        p = os.path.join(d, "x.bed.gz")
        with open(p, "wb") as fh:
            fh.write(b"X" * 200 + hf._BGZF_EOF)  # ends with the BGZF EOF marker
        self.assertIsNone(hf.verify_download_file(p))

    def test_truncated_bgzf_detected(self):
        d = self._tmp()
        p = os.path.join(d, "x.bed.gz")
        with open(p, "wb") as fh:
            fh.write(b"X" * 200)  # no BGZF EOF marker -> truncated
        reason = hf.verify_download_file(p)
        self.assertIsNotNone(reason)
        self.assertIn("EOF", reason)

    def test_empty_file_detected(self):
        d = self._tmp()
        p = os.path.join(d, "x.bed.gz")
        open(p, "wb").close()
        self.assertIn("empty", hf.verify_download_file(p) or "")

    def test_size_mismatch_detected(self):
        d = self._tmp()
        p = os.path.join(d, "y.fa.gz")
        with gzip.open(p, "wb") as fh:
            fh.write(b">c\nACGT\n")
        actual = os.path.getsize(p)
        self.assertIsNone(hf.verify_download_file(p, expected_size=actual))
        self.assertIn("truncated", hf.verify_download_file(p, expected_size=actual + 999) or "")

    def test_bad_gzip_magic_detected(self):
        d = self._tmp()
        p = os.path.join(d, "z.gz")
        with open(p, "wb") as fh:
            fh.write(b"not a gzip stream")
        self.assertIn("gzip", hf.verify_download_file(p) or "")

    def test_valid_targz_passes(self):
        d = self._tmp()
        p = os.path.join(d, "a.tar.gz")
        with tarfile.open(p, "w:gz"):
            pass
        self.assertIsNone(hf.verify_download_file(p))

    def test_deep_scan_catches_interior_corruption(self):
        d = self._tmp()
        p = os.path.join(d, "c.gz")
        with gzip.open(p, "wb") as fh:
            fh.write(b"A" * 4096)
        self.assertIsNone(hf.verify_download_file(p, deep=True))  # intact
        raw = bytearray(open(p, "rb").read())
        raw[len(raw) // 2] ^= 0xFF  # flip an interior byte
        with open(p, "wb") as fh:
            fh.write(raw)
        self.assertIsNotNone(hf.verify_download_file(p, deep=True))


class _FakeResp:
    def __init__(self, status, retry_after=None):
        self.status_code = status
        self.headers = {} if retry_after is None else {"Retry-After": retry_after}


class _HTTPErr(Exception):
    def __init__(self, msg, status=None, retry_after=None):
        super().__init__(msg)
        if status is not None:
            self.response = _FakeResp(status, retry_after)


class TestHfRetry(unittest.TestCase):
    """The retry-on-429/5xx helper (network-free): classification, Retry-After
    honoring, backoff, exhaustion, and the CRISPRME_HF_MAX_RETRIES override."""

    def test_retryable_classification(self):
        self.assertTrue(hf._hf_retryable(_HTTPErr("rate", status=429)))
        self.assertTrue(hf._hf_retryable(_HTTPErr("boom", status=503)))
        self.assertTrue(hf._hf_retryable(_HTTPErr("429 Too Many Requests")))
        self.assertTrue(hf._hf_retryable(_HTTPErr("503 Service Unavailable")))
        # 404 / auth / generic errors are NOT retried
        self.assertFalse(hf._hf_retryable(_HTTPErr("not found", status=404)))
        self.assertFalse(hf._hf_retryable(_HTTPErr("401 Unauthorized")))
        self.assertFalse(hf._hf_retryable(ValueError("bad arg")))

    def test_retry_after_parsed_and_capped(self):
        self.assertEqual(hf._hf_retry_after(_HTTPErr("x", 429, "12")), 12.0)
        self.assertIsNone(hf._hf_retry_after(_HTTPErr("x", 429, "99999")))  # capped
        self.assertIsNone(hf._hf_retry_after(_HTTPErr("x", 429, "nope")))  # unparsable
        self.assertIsNone(hf._hf_retry_after(_HTTPErr("x", 429)))  # no header

    def test_succeeds_after_transient_429(self):
        calls = {"n": 0}

        def flaky(**kw):
            calls["n"] += 1
            if calls["n"] < 3:
                raise _HTTPErr("429 Too Many Requests", status=429)
            return "ok"

        with mock.patch.object(hf.time, "sleep") as slept:
            out = hf._hf_call_with_retry(flaky, _what="test")
        self.assertEqual(out, "ok")
        self.assertEqual(calls["n"], 3)
        self.assertEqual(slept.call_count, 2)  # slept before each retry

    def test_honors_retry_after_header(self):
        seq = [_HTTPErr("429", 429, "7"), "ok"]

        def flaky(**kw):
            v = seq.pop(0)
            if isinstance(v, Exception):
                raise v
            return v

        with mock.patch.object(hf.time, "sleep") as slept:
            hf._hf_call_with_retry(flaky)
        slept.assert_called_once_with(7.0)

    def test_non_retryable_reraised_immediately(self):
        def boom(**kw):
            raise _HTTPErr("404 not found", status=404)

        with mock.patch.object(hf.time, "sleep") as slept:
            with self.assertRaises(_HTTPErr):
                hf._hf_call_with_retry(boom)
        slept.assert_not_called()

    def test_exhaustion_reraises_last(self):
        def always429(**kw):
            raise _HTTPErr("429", status=429)

        with mock.patch.dict(os.environ, {"CRISPRME_HF_MAX_RETRIES": "3"}):
            with mock.patch.object(hf.time, "sleep") as slept:
                with self.assertRaises(_HTTPErr):
                    hf._hf_call_with_retry(always429)
        self.assertEqual(slept.call_count, 2)  # 3 attempts => 2 sleeps

    def test_disable_via_env(self):
        calls = {"n": 0}

        def always429(**kw):
            calls["n"] += 1
            raise _HTTPErr("429", status=429)

        with mock.patch.dict(os.environ, {"CRISPRME_HF_MAX_RETRIES": "1"}):
            with mock.patch.object(hf.time, "sleep"):
                with self.assertRaises(_HTTPErr):
                    hf._hf_call_with_retry(always429)
        self.assertEqual(calls["n"], 1)  # no retry


class TestIndexAlreadyInstalled(unittest.TestCase):
    """`--what index` must NO-OP when the index is already installed.

    Without this, re-running a documented download command re-fetches and
    re-extracts the whole tarball (8.86 GB for NRG_3_hg38). The guard runs before
    any network call, so these tests make one by raising from both network entry
    points: reaching either means the guard did not short-circuit.
    """

    class _Reached(Exception):
        pass

    def _run(self, workdir, name, genotypes=True):
        def boom(*_a, **_k):
            raise self._Reached()

        with mock.patch.object(hf, "_hf_snapshot", boom), \
             mock.patch.object(hf, "_require_hf", boom):
            try:
                return hf.download_component(
                    "index", workdir, index_name=name, genotypes=genotypes
                )
            except self._Reached:
                return None  # proceeded to the network

    def _library(self, *names):
        root = tempfile.mkdtemp()
        for n in names:
            os.makedirs(os.path.join(root, "genome_library", n))
        return root

    def _genotype_store(self, root, vcf):
        # the Tier-1 genotype store lands under <workdir>/Dictionaries/genotypes_<vcf>/
        os.makedirs(os.path.join(root, "Dictionaries", f"genotypes_{vcf}"))

    # A COMPLETE variant install = main dir + _INDELS + shared reference index
    # (NRG_3_hg38) + genotype store (genotypes_hg38_X). Reference-only has none.

    def test_reference_only_installed_is_noop(self):
        wd = self._library("NRG_3_hg38")
        self.assertEqual(
            self._run(wd, "NRG_3_hg38"),
            os.path.join(wd, "genome_library", "NRG_3_hg38"),
        )

    def test_complete_variant_install_is_noop(self):
        wd = self._library(
            "NRG_3_hg38+hg38_X", "NRG_3_hg38+hg38_X_INDELS", "NRG_3_hg38"
        )
        self._genotype_store(wd, "hg38_X")
        self.assertEqual(
            self._run(wd, "NRG_3_hg38+hg38_X"),
            os.path.join(wd, "genome_library", "NRG_3_hg38+hg38_X"),
        )

    def test_dictless_name_resolves_to_canonical_dir(self):
        # a "-dictless" publish marker installs to the canonical folder, so a
        # request by either name must see the same (complete) install.
        wd = self._library(
            "NRG_3_hg38+hg38_X", "NRG_3_hg38+hg38_X_INDELS", "NRG_3_hg38"
        )
        self._genotype_store(wd, "hg38_X")
        self.assertEqual(
            self._run(wd, "NRG_3_hg38-dictless+hg38_X"),
            os.path.join(wd, "genome_library", "NRG_3_hg38+hg38_X"),
        )

    def test_variant_missing_indels_companion_still_downloads(self):
        # the two dirs are renamed into place by SEPARATE os.rename calls, so a
        # main dir without its companion is an interrupted install, not a
        # complete one -- it must NOT be skipped.
        wd = self._library("NRG_3_hg38+hg38_X", "NRG_3_hg38")
        self._genotype_store(wd, "hg38_X")
        self.assertIsNone(self._run(wd, "NRG_3_hg38+hg38_X"))

    def test_variant_missing_reference_index_still_downloads(self):
        # the shared reference index travels separately from the variant tarball;
        # absent -> re-run must heal it, not no-op.
        wd = self._library("NRG_3_hg38+hg38_X", "NRG_3_hg38+hg38_X_INDELS")
        self._genotype_store(wd, "hg38_X")
        self.assertIsNone(self._run(wd, "NRG_3_hg38+hg38_X"))

    def test_variant_missing_genotype_store_still_downloads(self):
        # genotypes default ON: a complete-but-for-the-store install must re-run to
        # heal the store (the documented "Samples degraded" recovery), not no-op.
        # This is the regression the coarse guard introduced.
        wd = self._library(
            "NRG_3_hg38+hg38_X", "NRG_3_hg38+hg38_X_INDELS", "NRG_3_hg38"
        )
        self.assertIsNone(self._run(wd, "NRG_3_hg38+hg38_X"))

    def test_no_genotypes_does_not_require_store(self):
        # with --no-genotypes the genotype store is not part of a complete install,
        # so an otherwise-complete install still no-ops.
        wd = self._library(
            "NRG_3_hg38+hg38_X", "NRG_3_hg38+hg38_X_INDELS", "NRG_3_hg38"
        )
        self.assertEqual(
            self._run(wd, "NRG_3_hg38+hg38_X", genotypes=False),
            os.path.join(wd, "genome_library", "NRG_3_hg38+hg38_X"),
        )

    def test_absent_index_still_downloads(self):
        wd = self._library("NRG_3_hg38")
        self.assertIsNone(self._run(wd, "NRG_3_hg38+hg38_ABSENT"))


class TestProgressReporter(unittest.TestCase):
    """The env-gated download-progress JSON file (2.8.1 web progress bar)."""

    def _tmp(self):
        import shutil
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(d, ignore_errors=True))
        return d

    def test_disabled_is_noop(self):
        # no CRISPRME_PROGRESS_FILE and no explicit path -> enabled False, writes
        # nothing (the current happy path must be byte-for-byte unchanged).
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("CRISPRME_PROGRESS_FILE", None)
            r = hf._ProgressReporter()
            self.assertFalse(r.enabled)
            r.set_phase("x")
            r.add_bytes(100)
            r.write(force=True)  # must not raise, must not create a file

    def test_env_gated_path(self):
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        with mock.patch.dict(os.environ, {"CRISPRME_PROGRESS_FILE": p}):
            r = hf._ProgressReporter()
            self.assertTrue(r.enabled)
            self.assertEqual(r.path, p)

    def test_bytes_done_pct_written(self):
        import json
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        r = hf._ProgressReporter(path=p, min_interval=0.0)
        r.set_total(1000)
        r.set_phase("Downloading index", "NRG_3_hg38")
        r.add_bytes(250)
        r.add_bytes(250)
        with open(p) as fh:
            data = json.load(fh)
        self.assertEqual(data["bytes_done"], 500)
        self.assertEqual(data["bytes_total"], 1000)
        self.assertEqual(data["pct"], 50.0)
        self.assertEqual(data["phase"], "Downloading index")
        self.assertEqual(data["detail"], "NRG_3_hg38")
        self.assertIn("updated", data)

    def test_throttle_and_force(self):
        import json
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        r = hf._ProgressReporter(path=p, min_interval=1000.0)  # effectively never
        r.set_total(100)
        r.set_phase("p")      # force=True -> writes
        r.add_bytes(10)       # throttled -> no write
        with open(p) as fh:
            self.assertEqual(json.load(fh)["bytes_done"], 0)
        r.write(force=True)   # explicit force -> flush latest
        with open(p) as fh:
            self.assertEqual(json.load(fh)["bytes_done"], 10)

    def test_eta_is_null_without_total(self):
        import json
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        r = hf._ProgressReporter(path=p, min_interval=0.0)
        r.set_phase("decompress")
        r.add_bytes(123)  # no total set -> pct 0, eta null, still records bytes
        with open(p) as fh:
            data = json.load(fh)
        self.assertEqual(data["bytes_total"], 0)
        self.assertEqual(data["pct"], 0.0)
        self.assertIsNone(data["eta_seconds"])
        self.assertEqual(data["bytes_done"], 123)

    def _eta_after(self, steps):
        """Drive a reporter on a fake clock; return the last eta_seconds written.

        ``steps`` are ``(advance_seconds, method, args)``."""
        import json
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        clock = [1000.0]
        with mock.patch.object(hf.time, "time", lambda: clock[0]):
            r = hf._ProgressReporter(path=p, min_interval=0.0)
            for advance, method, args in steps:
                clock[0] += advance
                getattr(r, method)(*args)
        with open(p) as fh:
            return json.load(fh)["eta_seconds"]

    def test_eta_is_this_stage_bytes_over_this_stage_time(self):
        # 5 GB in 10 min, then one minute into a 20 GB stage at the same pace: about
        # 9 min remain, not the ~99 min the time since the download began gives.
        eta = self._eta_after([
            (0, "set_total", (5 * 10**9,)), (0, "set_phase", ("Downloading index 1/3",)),
            (600, "add_bytes", (5 * 10**9,)),
            (0, "set_total", (20 * 10**9,)), (0, "set_phase", ("Downloading genotypes 2/3",)),
            (60, "add_bytes", (2 * 10**9,)),
        ])
        self.assertAlmostEqual(eta, 540, delta=1)

    def test_carried_bytes_keep_the_time_they_took(self):
        # carry_done keeps the running byte total, so the clock must not restart:
        # 2 GB in 200 s is 10 MB/s, leaving 800 s for the other 8 GB.
        eta = self._eta_after([
            (0, "set_total", (10 * 10**9,)), (0, "set_phase", ("part 1",)),
            (100, "add_bytes", (10**9,)),
            (0, "set_phase", ("part 2", "", True)),
            (100, "add_bytes", (10**9,)),
        ])
        self.assertAlmostEqual(eta, 800, delta=1)

    def test_tqdm_class_accumulates_bytes(self):
        # the custom tqdm_class' update(n) must feed the reporter so
        # snapshot_download byte progress lands in progress.json.
        import json
        d = self._tmp()
        p = os.path.join(d, "progress.json")
        r = hf._ProgressReporter(path=p, min_interval=0.0)
        r.set_total(300)
        r.set_phase("dl")

        # stub tqdm so the test needs neither tqdm nor network. The modules
        # themselves are faked: patching the attribute "tqdm.auto.tqdm" would
        # import the real tqdm.auto, which the light CI image does not have.
        import sys
        import types

        class _FakeTqdm:
            def __init__(self, *a, **k):
                pass
            def update(self, n=1):
                return None
        fake_auto = types.ModuleType("tqdm.auto")
        fake_auto.tqdm = _FakeTqdm
        fake_pkg = types.ModuleType("tqdm")
        fake_pkg.auto = fake_auto
        with mock.patch.dict(sys.modules, {"tqdm": fake_pkg, "tqdm.auto": fake_auto}):
            cls = hf._make_progress_tqdm(r)
            self.assertIsNotNone(cls)
            # unit="B" exactly as hf builds its byte bars -- that is the bar whose
            # updates are bytes, and the only kind that feeds the reporter
            bar = cls(total=300, unit="B", unit_scale=True)
            bar.update(120)
            bar.update(80)
        with open(p) as fh:
            self.assertEqual(json.load(fh)["bytes_done"], 200)


class TestProgressHookNeverBreaksDownload(unittest.TestCase):
    """Progress reporting is best-effort: the bar must never fail the download it
    reports on.

    huggingface_hub builds the snapshot bars from the caller's ``tqdm_class`` with
    its OWN keyword arguments -- ``name="huggingface_hub.snapshot_download"`` on
    versions 1.1.0 through 1.9.2, which plain tqdm rejects with TqdmKeyError unless
    the bar happens to be auto-disabled. huggingface_hub is unpinned
    (environment.yml, Dockerfile), so the hook has to survive whatever it is given.
    """

    class _TqdmKeyError(KeyError):
        """Stands in for tqdm.TqdmKeyError (also a KeyError)."""

    def _fake_tqdm_modules(self, update_raises=None):
        """Fake ``tqdm`` / ``tqdm.auto`` whose tqdm rejects unknown kwargs the way
        the real one does (the light CI image has no tqdm). ``update_raises``, if
        given, is raised by every ``update()`` of the base bar.

        Shaped like the real thing: ``tqdm.auto.tqdm`` is tqdm_asyncio, whose own
        ``__init__`` is ``(self, iterable=None, *args, **kwargs)`` and forwards to
        tqdm.std.tqdm, which is where total/unit/... are declared. A fake that
        declared them on the class itself would hide a keyword filter that reads
        only the immediate signature and so drops every one of them.

        And, as in tqdm >= 4.66, tqdm.std.tqdm's ``__init__`` is a
        ``functools.partialmethod`` (its ``@envwrap`` decorator returns one), which
        ``inspect.signature`` cannot read straight out of the class ``__dict__``.
        """
        import functools
        import sys
        import types

        unknown = self._TqdmKeyError

        def _std_init(self, iterable=None, desc=None, total=None, disable=False,
                      unit=None, unit_scale=False, initial=0, bar_format=None,
                      **kwargs):
            if kwargs:
                raise unknown("Unknown argument(s): " + str(kwargs))
            self.given = {"desc": desc, "total": total, "unit": unit,
                          "initial": initial}
            self.total = total
            self.n = initial

        class _FakeStdTqdm:
            __init__ = functools.partialmethod(_std_init)  # what @envwrap leaves

            def update(self, n=1):
                if update_raises is not None:
                    raise update_raises
                self.n += n or 0

        class _FakeAutoTqdm(_FakeStdTqdm):  # stands in for tqdm_asyncio
            def __init__(self, iterable=None, *args, **kwargs):
                super().__init__(iterable, *args, **kwargs)

        auto = types.ModuleType("tqdm.auto")
        auto.tqdm = _FakeAutoTqdm
        pkg = types.ModuleType("tqdm")
        pkg.auto = auto
        return mock.patch.dict(sys.modules, {"tqdm": pkg, "tqdm.auto": auto})

    def _tmp(self):
        import shutil
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(d, ignore_errors=True))
        return d

    def _reporter(self):
        return hf._ProgressReporter(path=os.path.join(self._tmp(), "p.json"),
                                    min_interval=0.0)

    def test_the_bar_drops_keywords_plain_tqdm_rejects(self):
        r = self._reporter()
        with self._fake_tqdm_modules():
            cls = hf._make_progress_tqdm(r)
            # exactly how hub 1.9.2 builds its aggregated bytes bar; the fake checks
            # keywords as a LIVE bar does (a disabled real one returns before that)
            bar = cls(desc="Downloading (incomplete total...)", disable=None, total=0,
                      initial=0, unit="B", unit_scale=True,
                      name="huggingface_hub.snapshot_download")
            bar.update(7)
        self.assertEqual(bar.given["unit"], "B")  # known kwargs still get through
        # hf does `bar.total += total` on this bar, so dropping total (not just the
        # keyword it does not know) would break the download a different way
        self.assertEqual(bar.given["total"], 0)
        self.assertEqual(bar.total, 0)
        self.assertEqual(bar.n, 7)                # and the base bar still advances
        self.assertEqual(r.bytes_done, 7)         # and the reporter still sees bytes

    def test_only_one_byte_bar_feeds_the_reporter(self):
        # hf drives SEVERAL bars from this one class: since hub 1.23.0 a network
        # "Downloading bytes" bar AND a "Reconstructing" bar, each fed the full byte
        # count (measured 2.00x the real size on 1.23.0-2.4.0), plus a "Fetching N
        # files" bar counting files. Summing them all put the web bar at 100% when
        # half the data had moved.
        r = self._reporter()
        with self._fake_tqdm_modules():
            cls = hf._make_progress_tqdm(r)
            transfer = cls(desc="Downloading bytes", total=0, unit="B", unit_scale=True)
            reconstruct = cls(desc="Reconstructing", total=0, unit="B", unit_scale=True)
            files = cls(desc="Fetching 2 files", total=2)  # no unit -> not bytes
            for bar in (transfer, reconstruct):
                bar.update(600)
            files.update(1)
        self.assertEqual(r.bytes_done, 600)
        self.assertEqual((transfer.n, reconstruct.n, files.n), (600, 600, 1))

    # --- the download goes ahead even when the hook cannot be used -------------

    class _FakeHub:
        """Records each snapshot_download. Given a ``tqdm_class`` it builds the bytes
        bar from it exactly as hub 1.9.2 does; ``fail_with`` is then raised by every
        call, as a download that really fails would."""

        def __init__(self, fail_with=None):
            self.calls = []
            self._fail_with = fail_with

        def snapshot_download(self, **kwargs):
            self.calls.append(kwargs)
            if "tqdm_class" in kwargs:
                kwargs["tqdm_class"](desc="Downloading (incomplete total...)",
                                     disable=None, total=0, initial=0, unit="B",
                                     unit_scale=True,
                                     name="huggingface_hub.snapshot_download")
            if self._fail_with is not None:
                raise self._fail_with
            return kwargs["local_dir"]

    class _BarBuildingHub:
        """Builds and feeds its byte bars from ``tqdm_class`` on EVERY call, as hub
        2.3.0 does: a network-transfer bar, then a reconstruction bar, both fed each
        chunk; resumed bytes go to the reconstruction bar only. The first call is
        rate-limited (429) after ``cut`` of ``size`` bytes."""

        def __init__(self, size=1000, cut=300, chunk=100):
            self.size, self.cut, self.chunk = size, cut, chunk
            self.calls = []

        def snapshot_download(self, **kwargs):
            self.calls.append(kwargs)
            cls = kwargs["tqdm_class"]
            transfer = cls(desc="Downloading bytes", total=0, initial=0, unit="B",
                           unit_scale=True)
            reconstruct = cls(desc="Reconstructing", total=0, initial=0, unit="B",
                              unit_scale=True)
            done = 0 if len(self.calls) == 1 else self.cut
            reconstruct.update(done)  # what was already on disk (resume)
            while done < self.size:
                if len(self.calls) == 1 and done == self.cut:
                    raise _HTTPErr("429 Too Many Requests", status=429)
                for bar in (transfer, reconstruct):
                    bar.update(self.chunk)
                done += self.chunk
            return kwargs["local_dir"]

    def _snapshot(self, hub, reporter):
        """(stdout, stderr) of one _hf_snapshot; also kept on self.printed for a
        call that raises."""
        import io
        self.printed = io.StringIO()
        with mock.patch.object(hf, "_require_hf", return_value=hub), \
             mock.patch.dict(os.environ, {"CRISPRME_SKIP_VERIFY": "1"}), \
             mock.patch.object(hf.time, "sleep"), \
             mock.patch("sys.stdout", self.printed), \
             mock.patch("sys.stderr", new_callable=io.StringIO) as err, \
             self._fake_tqdm_modules():
            hf._hf_snapshot("me/repo", ["indexes/x.tar.gz"],
                            os.path.join(self._tmp(), "dl"), None, reporter=reporter)
        return self.printed.getvalue(), err.getvalue()

    def test_a_retried_download_keeps_feeding_the_bar(self):
        # _hf_call_with_retry exists for this: a 429 mid-download, then the rest.
        # hf builds new bars on every attempt, so the bar must keep counting through
        # them; it used to stay at the first attempt's 30% until extraction began.
        hub = self._BarBuildingHub(size=1000, cut=300)
        r = self._reporter()
        r.set_total(1000)
        printed, _ = self._snapshot(hub, r)
        self.assertEqual(len(hub.calls), 2)
        self.assertTrue(all("tqdm_class" in call for call in hub.calls))
        self.assertEqual(r.bytes_done, 1000)
        self.assertEqual(printed, "")

    def test_a_rejected_bar_retries_the_download_without_it(self):
        # a tqdm whose signature cannot be read gets every keyword passed on; when
        # it then rejects one, the BAR has failed, so the data is fetched plain
        hub = self._FakeHub()
        with mock.patch("inspect.signature", side_effect=ValueError("no signature")):
            printed, errored = self._snapshot(hub, self._reporter())
        self.assertEqual(len(hub.calls), 2)
        self.assertIn("tqdm_class", hub.calls[0])
        self.assertNotIn("tqdm_class", hub.calls[1])  # degraded progress, real download
        self.assertEqual(hub.calls[1]["allow_patterns"], ["indexes/x.tar.gz"])
        self.assertIn("progress bar", printed)
        self.assertIn("Unknown argument", printed)  # names what the bar raised
        self.assertEqual(errored, "")  # diagnostics on stdout, like the rest here

    def test_a_working_bar_is_kept_and_nothing_is_retried(self):
        hub = self._FakeHub()
        printed, _ = self._snapshot(hub, self._reporter())
        self.assertEqual(len(hub.calls), 1)
        self.assertIn("tqdm_class", hub.calls[0])
        self.assertEqual(printed, "")

    def test_a_download_that_really_fails_still_fails(self):
        hub = self._FakeHub(fail_with=RuntimeError("404 Client Error: Not Found"))
        with self.assertRaises(RuntimeError):
            self._snapshot(hub, self._reporter())

    # --- ...and only when it is the BAR that failed ----------------------------
    # Anything else propagates exactly as without a bar: the same attempts as on
    # b5d6cc1, and no NOTE blaming the bar for a failure it did not cause.

    def test_what_the_bar_itself_raises_is_recorded(self):
        failures = []
        with self._fake_tqdm_modules(), \
             mock.patch("inspect.signature", side_effect=ValueError("no signature")):
            cls = hf._make_progress_tqdm(self._reporter(), failures)
            with self.assertRaises(self._TqdmKeyError) as built:
                cls(unit="B", name="huggingface_hub.snapshot_download")
        boom = RuntimeError("the bar could not advance")
        with self._fake_tqdm_modules(update_raises=boom):
            bar = hf._make_progress_tqdm(self._reporter(), failures)(unit="B")
            with self.assertRaises(RuntimeError):
                bar.update(5)
        self.assertEqual(failures, [built.exception, boom])

    BUDGET = 4  # CRISPRME_HF_MAX_RETRIES in the tests below

    def _fails_every_attempt(self, exc):
        hub = self._FakeHub(fail_with=exc)
        with mock.patch.dict(os.environ, {"CRISPRME_HF_MAX_RETRIES": str(self.BUDGET)}), \
             self.assertRaises(type(exc)) as caught:
            self._snapshot(hub, self._reporter())
        self.assertIs(caught.exception, exc)  # re-raised unchanged
        return hub

    def test_a_persistent_503_uses_the_retry_budget_once(self):
        hub = self._fails_every_attempt(_HTTPErr("503 Service Unavailable", status=503))
        self.assertEqual(len(hub.calls), self.BUDGET)
        self.assertEqual(self.printed.getvalue(), "")

    def test_a_401_fails_at_once(self):
        hub = self._fails_every_attempt(_HTTPErr("401 Client Error: Unauthorized",
                                                 status=401))
        self.assertEqual(len(hub.calls), 1)
        self.assertEqual(self.printed.getvalue(), "")

    def test_with_progress_off_the_call_is_unchanged(self):
        hub = self._FakeHub()
        r = hf._ProgressReporter(path=None)
        self.assertFalse(r.enabled)
        self._snapshot(hub, r)
        self.assertEqual(len(hub.calls), 1)
        self.assertNotIn("tqdm_class", hub.calls[0])


def _mountinfo_line(n, root, mount_point):
    """One /proc/self/mountinfo line: field 4 is the mounted folder's path inside
    its own filesystem, field 5 the mount point, both escaped as the kernel writes
    them (space, tab, newline and backslash in octal, e.g. \\040)."""
    def escape(path):
        return "".join(f"\\{ord(c):03o}" if c in " \t\n\\" else c for c in path)
    return (f"{30 + n} 1 0:52 {escape(root)} {escape(mount_point)} rw,relatime "
            f"shared:1 - ext4 /dev/vda1 rw\n")


class TestAtomicInstallDir(unittest.TestCase):
    """_atomic_install_dir: a clean build swaps in; a failed build leaves NO
    partial in the final location (the key 2.8.1 robustness fix)."""

    def _tmp(self):
        import shutil
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(d, ignore_errors=True))
        return d

    def test_clean_build_installs(self):
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")

        def build(stage):
            with open(os.path.join(stage, "chr1.fa"), "w") as fh:
                fh.write(">chr1\nACGT\n")

        out = hf._atomic_install_dir(dest, build)
        self.assertEqual(out, dest)
        self.assertTrue(os.path.isfile(os.path.join(dest, "chr1.fa")))
        # no staging / backup leftovers
        self.assertFalse(os.path.exists(os.path.join(d, "Genomes", ".hg38.staging")))
        self.assertFalse(os.path.exists(os.path.join(d, "Genomes", ".hg38.replaced")))

    def test_failed_build_leaves_no_partial(self):
        d = self._tmp()
        dest = os.path.join(d, "VCFs", "mydata")

        def build(stage):
            with open(os.path.join(stage, "chr1.vcf.gz"), "wb") as fh:
                fh.write(b"\x00" * 10)
            raise RuntimeError("simulated mid-download interruption")

        with self.assertRaises(RuntimeError):
            hf._atomic_install_dir(dest, build)
        # the final location must NOT exist (no half-installed dataset)
        self.assertFalse(os.path.exists(dest))
        # staging must be cleaned up
        self.assertFalse(os.path.exists(os.path.join(d, "VCFs", ".mydata.staging")))

    def test_failed_build_preserves_existing_dest(self):
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")
        os.makedirs(dest)
        with open(os.path.join(dest, "old.fa"), "w") as fh:
            fh.write("existing")

        def build(stage):
            with open(os.path.join(stage, "new.fa"), "w") as fh:
                fh.write("new")
            raise RuntimeError("boom")

        with self.assertRaises(RuntimeError):
            hf._atomic_install_dir(dest, build)
        # the pre-existing good copy must be untouched
        self.assertTrue(os.path.isfile(os.path.join(dest, "old.fa")))
        self.assertFalse(os.path.exists(os.path.join(dest, "new.fa")))

    def test_verify_rejects_before_swap(self):
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")

        def build(stage):
            open(os.path.join(stage, "truncated.fa.gz"), "wb").close()

        def verify(stage):
            raise hf.DownloadIntegrityError("truncated")

        with self.assertRaises(hf.DownloadIntegrityError):
            hf._atomic_install_dir(dest, build, verify=verify)
        self.assertFalse(os.path.exists(dest))

    # A dest that is a symlink to a shared copy (Genomes/hg38 -> /shared/hg38): the
    # download replaces the LINK with a real directory and never touches the target.

    def _linked_dest(self):
        d = self._tmp()
        shared = os.path.join(d, "shared", "hg38")
        os.makedirs(shared)
        with open(os.path.join(shared, "old.fa"), "w") as fh:
            fh.write("shared")
        dest = os.path.join(d, "Genomes", "hg38")
        os.makedirs(os.path.dirname(dest))
        os.symlink(shared, dest)
        return d, shared, dest

    @staticmethod
    def _build_new(stage):
        with open(os.path.join(stage, "new.fa"), "w") as fh:
            fh.write("new")

    def _install(self, dest, build):
        import io
        with mock.patch("sys.stdout", new_callable=io.StringIO) as out:
            hf._atomic_install_dir(dest, build)
        return out.getvalue()

    def test_symlinked_dest_replaces_the_link_not_its_target(self):
        d, shared, dest = self._linked_dest()
        printed = self._install(dest, self._build_new)
        self.assertFalse(os.path.islink(dest))
        self.assertEqual(sorted(os.listdir(dest)), ["new.fa"])
        self.assertEqual(sorted(os.listdir(shared)), ["old.fa"])  # target untouched
        genomes = os.path.join(d, "Genomes")
        self.assertEqual(os.listdir(genomes), ["hg38"])  # no .replaced link left behind
        self.assertEqual(printed.count("\n"), 1, printed)
        self.assertIn(shared, printed)
        self.assertIn("untouched", printed)

    def test_reinstalling_over_a_replaced_link_succeeds(self):
        # the second download used to hit ENOTDIR renaming the real dir onto the
        # leftover .hg38.replaced link, after the whole build, and leak the stage
        d, shared, dest = self._linked_dest()
        self._install(dest, self._build_new)
        self._install(dest, self._build_new)
        self.assertEqual(sorted(os.listdir(dest)), ["new.fa"])
        self.assertEqual(sorted(os.listdir(shared)), ["old.fa"])
        self.assertEqual(os.listdir(os.path.join(d, "Genomes")), ["hg38"])

    def test_a_leftover_backup_link_is_unlinked_not_followed(self):
        # what an earlier run left on disk: a private dest plus a backup that is a
        # link to the shared copy
        d, shared, dest = self._linked_dest()
        os.unlink(dest)
        os.makedirs(dest)
        backup = os.path.join(d, "Genomes", ".hg38.replaced")
        os.symlink(shared, backup)
        self._install(dest, self._build_new)
        self.assertFalse(os.path.lexists(backup))
        self.assertEqual(sorted(os.listdir(dest)), ["new.fa"])
        self.assertEqual(sorted(os.listdir(shared)), ["old.fa"])

    def _failing_rename(self, failing_src, err=16):
        real_rename = os.rename

        def rename(src, dst):
            if src == failing_src:
                raise OSError(err, os.strerror(err))
            return real_rename(src, dst)
        return mock.patch.object(hf.os, "rename", rename)

    def test_a_failed_move_aside_discards_the_stage(self):
        import errno
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")
        os.makedirs(dest)
        with open(os.path.join(dest, "old.fa"), "w") as fh:
            fh.write("existing")
        with self._failing_rename(dest, errno.EACCES), \
             self.assertRaises(OSError) as caught:
            hf._atomic_install_dir(dest, self._build_new)
        self.assertEqual(caught.exception.errno, errno.EACCES)  # passed on as it was
        self.assertEqual(os.listdir(os.path.join(d, "Genomes")), ["hg38"])  # no stage
        self.assertEqual(sorted(os.listdir(dest)), ["old.fa"])

    def test_a_mount_point_dest_is_refused_before_anything_is_touched(self):
        # Measured on a real macOS mount: the kernel renames a mount point (the
        # mount moves with it), and clearing the moved-aside copy then deleted the
        # mounted folder's own files. So it is refused before any build or rename,
        # here a bind mount from the same filesystem that only the mount table shows.
        import errno
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")
        os.makedirs(dest)
        with open(os.path.join(dest, "old.fa"), "w") as fh:
            fh.write("shared")
        mountinfo = os.path.join(d, "mountinfo")
        with open(mountinfo, "w") as fh:
            fh.write(_mountinfo_line(0, "/", "/"))
            fh.write(_mountinfo_line(1, "/srv/hg38", os.path.realpath(dest)))
        built = []
        with mock.patch.object(hf, "_MOUNTINFO", mountinfo), \
             self.assertRaises(OSError) as caught:
            hf._atomic_install_dir(dest, built.append)
        self.assertEqual(caught.exception.errno, errno.EBUSY)
        self.assertIn(f"{dest} is a mount point", str(caught.exception))
        self.assertEqual(built, [])
        self.assertEqual(sorted(os.listdir(dest)), ["old.fa"])
        self.assertEqual(os.listdir(os.path.join(d, "Genomes")), ["hg38"])

    def test_a_mount_point_dest_says_why_it_cannot_be_replaced(self):
        # a folder mounted on its own cannot be renamed aside (EBUSY). Without a
        # readable mount table, os.path.ismount misses a bind mount from the same
        # filesystem, so this rename is where it shows: say what to do, and leave
        # everything as it was
        import errno
        d = self._tmp()
        dest = os.path.join(d, "Genomes", "hg38")
        os.makedirs(dest)
        with open(os.path.join(dest, "old.fa"), "w") as fh:
            fh.write("existing")
        with mock.patch.object(hf, "_MOUNTINFO", os.path.join(d, "absent")), \
             self._failing_rename(dest, errno.EBUSY), \
             self.assertRaises(OSError) as caught:
            hf._atomic_install_dir(dest, self._build_new)
        self.assertEqual(caught.exception.errno, errno.EBUSY)
        self.assertIn(f"{dest} is a mount point", str(caught.exception))
        self.assertIn("mount its parent folder instead", str(caught.exception))
        self.assertEqual(os.listdir(os.path.join(d, "Genomes")), ["hg38"])  # no stage
        self.assertEqual(sorted(os.listdir(dest)), ["old.fa"])

    def test_a_failed_swap_puts_the_link_back(self):
        d, shared, dest = self._linked_dest()
        stage = os.path.join(d, "Genomes", ".hg38.staging")
        with self._failing_rename(stage), self.assertRaises(OSError):
            hf._atomic_install_dir(dest, self._build_new)
        self.assertTrue(os.path.islink(dest))
        self.assertEqual(os.readlink(dest), shared)
        self.assertEqual(os.listdir(os.path.join(d, "Genomes")), ["hg38"])
        self.assertEqual(sorted(os.listdir(shared)), ["old.fa"])


class TestMountPointInstallFolder(unittest.TestCase):
    """A Genomes/<ref> or VCFs/<dataset> folder mounted on its own (e.g. docker run
    -v shared:/data/Genomes/hg38) cannot be replaced: the atomic install renames it
    aside and the kernel refuses (EBUSY). That used to surface as a bare errno
    after the whole multi-GB download; it is now refused before anything is
    fetched, with what to do instead. A bind mount whose source is on the same
    filesystem as its parent (Docker Desktop host binds, Apptainer --bind, mount
    --bind) has its parent's device, so os.path.ismount misses it; on Linux the
    process's mount table (/proc/self/mountinfo, faked here) shows it."""

    class _Reached(Exception):
        """The download was started."""

    def setUp(self):
        import shutil
        real = tempfile.mkdtemp(prefix="crisprme mount ")  # the table writes \040
        self.addCleanup(shutil.rmtree, real, True)
        # the work dir is reached through a link; the table lists real paths
        self.workdir = real + " link"
        os.symlink(real, self.workdir)
        self.addCleanup(os.unlink, self.workdir)
        self.mountinfo = os.path.join(real, "mountinfo")  # unreadable until written

    def _folder(self, component, name):
        """The install folder of ``name``, holding one file already."""
        dest = os.path.join(self.workdir, hf._COMPONENT_LOCALDIR[component], name)
        os.makedirs(dest)
        with open(os.path.join(dest, "old"), "w") as fh:
            fh.write("old")
        return dest

    def _mount_table(self, *mounts):
        """Write the fake mount table: the root mount plus ``mounts``, each a
        (folder inside its filesystem, mount point) pair."""
        with open(self.mountinfo, "w") as fh:
            for n, (root, point) in enumerate((("/", "/"),) + mounts):
                fh.write(_mountinfo_line(n, root, point))

    def _download(self, component, ismount=False, **kw):
        """download_component under the fake mount table, with os.path.ismount
        answering ``ismount`` and the network cut off at its first call."""
        def reached(*_a, **_k):
            raise self._Reached()

        with mock.patch.object(hf, "_MOUNTINFO", self.mountinfo), \
             mock.patch.object(hf.os.path, "ismount", return_value=ismount), \
             mock.patch.object(hf, "_hf_snapshot", reached), \
             mock.patch.object(hf, "_require_hf", reached):
            hf.download_component(component, self.workdir, **kw)

    def _refused(self, component, **kw):
        with self.assertRaises(ValueError) as caught:
            self._download(component, **kw)
        return str(caught.exception)

    def test_a_mounted_genome_folder_is_refused_before_the_download(self):
        # a bind mount from the same filesystem: os.path.ismount says no, and the
        # mount table lists the folder's real path, its spaces written as \040
        dest = self._folder("genome", "hg38")
        self._mount_table(("/srv/hg38", os.path.realpath(dest)))
        message = self._refused("genome", ref="hg38")
        self.assertIn(f"{dest} is a mount point, so the download cannot replace it",
                      message)
        self.assertIn("mount its parent folder instead", message)
        self.assertEqual(os.listdir(dest), ["old"])

    def test_a_mounted_vcf_folder_is_refused_before_the_download(self):
        dest = self._folder("vcf", "hg38_1000G")
        self._mount_table(("/srv/hg38_1000G", os.path.realpath(dest)))
        message = self._refused("vcf", dataset="hg38_1000G")
        self.assertIn(f"{dest} is a mount point", message)
        self.assertIn("mount its parent folder instead", message)
        self.assertEqual(os.listdir(dest), ["old"])

    def test_a_mount_point_with_a_carriage_return_is_still_seen(self):
        # the kernel escapes newline but not \r, so the table holds a \r raw; a
        # parser that also split lines on \r cut that line and missed the mount
        self.workdir = os.path.join(os.path.realpath(self.workdir), "work\rdir")
        os.makedirs(self.workdir)
        dest = self._folder("genome", "hg38")
        self._mount_table(("/srv/hg38", os.path.realpath(dest)))
        message = self._refused("genome", ref="hg38")
        self.assertIn(f"{dest} is a mount point", message)
        self.assertEqual(os.listdir(dest), ["old"])

    def test_an_ordinary_folder_still_downloads(self):
        # the table lists the work dir, and the folder only as what another mount
        # point shows (field 4); neither makes the folder a mount point
        for component, kw in (("genome", {"ref": "hg38"}),
                              ("vcf", {"dataset": "hg38_1000G"})):
            with self.subTest(component=component):
                dest = self._folder(component, kw.get("ref") or kw.get("dataset"))
                self._mount_table(("/srv/work", os.path.realpath(self.workdir)),
                                  (os.path.realpath(dest), "/elsewhere"))
                with self.assertRaises(self._Reached):
                    self._download(component, **kw)

    def test_without_a_readable_mount_table_os_path_ismount_decides(self):
        # as on macOS, which has no /proc/self/mountinfo
        self._folder("genome", "hg38")
        self.assertFalse(os.path.exists(self.mountinfo))
        message = self._refused("genome", ref="hg38", ismount=True)
        self.assertIn("/Genomes/hg38 is a mount point", message)
        with self.assertRaises(self._Reached):
            self._download("genome", ref="hg38", ismount=False)

    def test_a_link_to_a_mounted_folder_still_downloads(self):
        # the install replaces the link itself and never touches its target, so a
        # link is never a mount point, as for os.path.ismount
        shared = os.path.join(os.path.realpath(self.workdir), "shared")
        os.makedirs(shared)
        os.makedirs(os.path.join(self.workdir, "Genomes"))
        os.symlink(shared, os.path.join(self.workdir, "Genomes", "hg38"))
        self._mount_table(("/srv/hg38", shared))
        with self.assertRaises(self._Reached):
            self._download("genome", ref="hg38")


class TestStructuralCheckWhenSizesEmpty(unittest.TestCase):
    """(2.8.1 piece C) when the HF size API fails (_remote_sizes -> {}), the
    download must STILL run the BGZF/gzip structural checks — only the size
    comparison is skipped."""

    def _tmp(self):
        import shutil
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(d, ignore_errors=True))
        return d

    def test_truncated_bgzf_caught_with_no_sizes(self):
        d = self._tmp()
        # expected_size=None (the sizes.get(rel) value when _remote_sizes=={}) must
        # NOT bypass the structural check
        p = os.path.join(d, "x.vcf.gz")
        with open(p, "wb") as fh:
            fh.write(b"X" * 100)  # no BGZF EOF -> truncated
        reason = hf.verify_download_file(p, expected_size=None)
        self.assertIsNotNone(reason)
        self.assertIn("EOF", reason)

    def test_remote_sizes_empty_on_api_failure(self):
        # _remote_sizes swallows any HF error and returns {} (structural-only mode)
        class _Boom:
            def get_paths_info(self, *a, **k):
                raise RuntimeError("size API down")
        self.assertEqual(
            hf._remote_sizes(_Boom(), "repo", ["a", "b"], None), {}
        )


if __name__ == "__main__":
    unittest.main()

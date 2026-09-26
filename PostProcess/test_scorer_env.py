"""Unit tests for scorer_env (env-manager lifecycle + state, no conda needed)."""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scorer_env as se  # noqa: E402


class TestSpec(unittest.TestCase):
    def test_cbulge_pins_present(self):
        spec = se.SCORER_ENVS["cbulge"]
        cpu = spec["cpu_packages"]
        self.assertIn("numpy=1.23.5", cpu)          # CRISPR-Bulge requirement
        self.assertIn("setuptools<81", cpu)          # pkg_resources restored
        self.assertIn("xgboost", cpu)                # pulled by encoder import chain
        self.assertIn("catboost", cpu)
        self.assertTrue(any(p.startswith("tensorflow-cpu") for p in cpu))
        self.assertTrue(any("cuda" in p for p in spec["gpu_packages"]))
        self.assertEqual(spec["python"], "3.10")


class TestDetect(unittest.TestCase):
    def test_override_wins(self):
        with mock.patch.dict(os.environ, {"CRISPRME_CONDA_EXE": "/opt/x/micromamba"}, clear=True), \
             mock.patch("os.path.exists", return_value=True), \
             mock.patch("os.access", return_value=True):
            exe, kind = se.detect_env_manager()
            self.assertEqual(kind, "micromamba")
            self.assertTrue(exe.endswith("micromamba"))

    def test_none_when_absent(self):
        with mock.patch.dict(os.environ, {}, clear=True), \
             mock.patch("shutil.which", return_value=None), \
             mock.patch("os.path.exists", return_value=False):
            self.assertIsNone(se.detect_env_manager())

    def test_kind_of(self):
        self.assertEqual(se._kind_of("/a/micromamba"), "micromamba")
        self.assertEqual(se._kind_of("/a/mamba"), "mamba")
        self.assertEqual(se._kind_of("/a/conda"), "conda")


class TestCreateCommand(unittest.TestCase):
    def test_cpu_command_shape(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/micromamba", "micromamba")):
            cmd = se.build_create_command("cbulge", gpu=False)
        self.assertEqual(cmd[:5], ["/m/micromamba", "create", "-y", "-n", "cbulge"])
        self.assertIn("python=3.10", cmd)
        self.assertIn("-c", cmd)
        self.assertIn("conda-forge", cmd)
        self.assertIn("numpy=1.23.5", cmd)
        self.assertIn("setuptools<81", cmd)
        self.assertFalse(any("cuda" in c for c in cmd))

    def test_gpu_command_uses_cuda(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/mamba", "mamba")):
            cmd = se.build_create_command("cbulge", gpu=True)
        self.assertTrue(any("cuda" in c for c in cmd))

    def test_no_manager_returns_none(self):
        with mock.patch.object(se, "detect_env_manager", return_value=None):
            self.assertIsNone(se.build_create_command("cbulge"))


class TestProvision(unittest.TestCase):
    def test_default_repo_precedence(self):
        with mock.patch.dict(os.environ, {"CBULGE_REPO": "/explicit"}, clear=False):
            self.assertEqual(se.default_cbulge_repo(), "/explicit")
        with mock.patch.dict(os.environ, {"CRISPRME_CBULGE_HOME": "/home2"}, clear=False):
            os.environ.pop("CBULGE_REPO", None)
            self.assertEqual(se.default_cbulge_repo(), "/home2")
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("CBULGE_REPO", None)
            os.environ.pop("CRISPRME_CBULGE_HOME", None)
            # derived: <prefix>/opt/CRISPR-Bulge, ends with CRISPR-Bulge
            self.assertTrue(se.default_cbulge_repo().endswith("CRISPR-Bulge"))

    def test_pin_is_full_sha(self):
        self.assertEqual(len(se.CBULGE_PIN), 40)
        self.assertTrue(se.CBULGE_URL.endswith("CRISPR-Bulge.git"))

    def test_provision_noop_when_present(self):
        with mock.patch.object(se, "source_present", return_value=True):
            ok, msg = se.provision_source("/some/repo")
            self.assertTrue(ok)
            self.assertIn("present", msg)

    def test_provision_fails_without_git(self):
        with mock.patch.object(se, "source_present", return_value=False), \
             mock.patch("shutil.which", return_value=None):
            ok, msg = se.provision_source("/some/repo")
            self.assertFalse(ok)
            self.assertIn("git", msg)


class TestState(unittest.TestCase):
    def test_roundtrip_and_version_stamp(self):
        with tempfile.TemporaryDirectory() as d:
            se.set_scorer_env_state(d, {"status": "ok", "env": "cbulge"})
            rec = se.get_scorer_env_state(d)
            self.assertEqual(rec["status"], "ok")
            self.assertEqual(rec["version"], 1)

    def test_missing_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(se.get_scorer_env_state(d), {})

    def test_corrupt_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            with open(se.state_path(d), "w") as fh:
                fh.write("{not json")
            self.assertEqual(se.get_scorer_env_state(d), {})

    def test_atomic_write_leaves_no_tmp(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertTrue(se.set_scorer_env_state(d, {"status": "warn"}))
            leftovers = [f for f in os.listdir(d) if f.endswith(".tmp")]
            self.assertEqual(leftovers, [])

    def test_write_to_readonly_is_best_effort(self):
        # a write failure must return False, not raise (diagnostic cache)
        ok = se.set_scorer_env_state("/proc/nonexistent/cannot/create", {"status": "ok"})
        self.assertFalse(ok)

    def test_health_check_unknown_env(self):
        rec = se.health_check("does-not-exist")
        self.assertEqual(rec["status"], se.ERROR)
        self.assertTrue(any("unknown scorer env" in m for _, m in rec["issues"]))


class TestHealthCheck(unittest.TestCase):
    def _probe_result(self, bad=None, vers=None):
        cp = mock.Mock()
        cp.returncode = 0
        cp.stdout = json.dumps({"py": "3.10.13", "vers": vers or {"tensorflow": "2.13.1"},
                                "bad": bad or []})
        cp.stderr = ""
        return cp

    def test_no_manager_is_error(self):
        with mock.patch.object(se, "detect_env_manager", return_value=None):
            rec = se.health_check("cbulge")
        self.assertEqual(rec["status"], se.ERROR)

    def test_missing_env_is_error(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/mamba", "mamba")), \
             mock.patch.object(se, "env_python", return_value=None):
            rec = se.health_check("cbulge")
        self.assertEqual(rec["status"], se.ERROR)
        self.assertTrue(any("does not exist" in m for _, m in rec["issues"]))

    def test_bad_import_is_error(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/mamba", "mamba")), \
             mock.patch.object(se, "env_python", return_value="/fake/py"), \
             mock.patch.object(se, "_run", return_value=self._probe_result(bad=[["tensorflow", "No module"]])):
            rec = se.health_check("cbulge")
        self.assertEqual(rec["status"], se.ERROR)

    def test_ok_when_deps_and_weights_present(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/mamba", "mamba")), \
             mock.patch.object(se, "env_python", return_value="/fake/py"), \
             mock.patch.object(se, "_run", return_value=self._probe_result()), \
             mock.patch.dict(os.environ, {"CBULGE_REPO": "/repo"}, clear=False), \
             mock.patch("os.path.isdir", return_value=True):
            rec = se.health_check("cbulge")
        self.assertEqual(rec["status"], se.OK)
        self.assertTrue(rec["weights_present"])

    def test_warn_when_weights_absent(self):
        with mock.patch.object(se, "detect_env_manager", return_value=("/m/mamba", "mamba")), \
             mock.patch.object(se, "env_python", return_value="/fake/py"), \
             mock.patch.object(se, "_run", return_value=self._probe_result()), \
             mock.patch.dict(os.environ, {"CBULGE_REPO": ""}, clear=False):
            rec = se.health_check("cbulge")
        self.assertEqual(rec["status"], se.WARN)

    def test_render_health_smoke(self):
        rec = {"env": "cbulge", "status": "warn", "issues": [("warn", "no weights")],
               "manager": {"kind": "mamba", "exe": "/m/mamba"}, "versions": {"numpy": "1.23.5"},
               "env_python": "/p", "python_version": "3.10.13"}
        text = se.render_health(rec)
        self.assertIn("WARN", text)
        self.assertIn("numpy=1.23.5", text)


if __name__ == "__main__":
    unittest.main()

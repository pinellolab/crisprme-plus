"""Unit tests for compute_backend: the cpu|cuda|metal selection + fallback matrix.

Pure stdlib + mocks (no tensorflow/cupy needed) — validates that a requested backend
degrades to CPU when the hardware isn't present, and that the pre-import env is set
correctly. The tensorflow-metal GRU miscompute + numerical self-test fallback is
validated separately (that lives in crispr_bulge_score.load_models, needs the model)."""
import os
import unittest
from unittest import mock

import compute_backend as cb


class TestResolveBackend(unittest.TestCase):
    def _resolve(self, requested, apple=False, cuda=False):
        with mock.patch.object(cb, "_is_apple_silicon", return_value=apple), \
             mock.patch.object(cb, "_has_cuda", return_value=cuda):
            return cb.resolve_backend(requested)

    def test_cpu_is_default_and_explicit(self):
        self.assertEqual(self._resolve("cpu"), cb.CPU)
        self.assertEqual(self._resolve(None), cb.CPU)          # unset env default
        self.assertEqual(self._resolve(""), cb.CPU)

    def test_cuda_present_and_absent(self):
        self.assertEqual(self._resolve("cuda", cuda=True), cb.CUDA)
        self.assertEqual(self._resolve("cuda", cuda=False), cb.CPU)   # graceful fallback

    def test_metal_present_and_absent(self):
        self.assertEqual(self._resolve("metal", apple=True), cb.METAL)
        self.assertEqual(self._resolve("metal", apple=False), cb.CPU)  # not a Mac -> cpu

    def test_auto_prefers_cuda_then_metal_then_cpu(self):
        self.assertEqual(self._resolve("gpu", cuda=True, apple=True), cb.CUDA)
        self.assertEqual(self._resolve("auto", cuda=False, apple=True), cb.METAL)
        self.assertEqual(self._resolve("gpu", cuda=False, apple=False), cb.CPU)

    def test_case_insensitive(self):
        self.assertEqual(self._resolve("METAL", apple=True), cb.METAL)
        self.assertEqual(self._resolve("Cuda", cuda=True), cb.CUDA)

    def test_unknown_value_falls_back_to_cpu(self):
        self.assertEqual(self._resolve("gpu-turbo"), cb.CPU)

    def test_env_var_is_read_when_no_arg(self):
        with mock.patch.dict(os.environ, {"CRISPRME_COMPUTE_BACKEND": "metal"}, clear=False), \
             mock.patch.object(cb, "_is_apple_silicon", return_value=True), \
             mock.patch.object(cb, "_has_cuda", return_value=False):
            self.assertEqual(cb.resolve_backend(), cb.METAL)


class TestPreImportEnv(unittest.TestCase):
    def test_cpu_hides_cuda_devices(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("CUDA_VISIBLE_DEVICES", None)
            cb.pre_import_env(cb.CPU)
            self.assertEqual(os.environ.get("CUDA_VISIBLE_DEVICES"), "")

    def test_accelerators_leave_cuda_env_untouched(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("CUDA_VISIBLE_DEVICES", None)
            cb.pre_import_env(cb.CUDA)
            self.assertNotIn("CUDA_VISIBLE_DEVICES", os.environ)
            cb.pre_import_env(cb.METAL)
            self.assertNotIn("CUDA_VISIBLE_DEVICES", os.environ)


class TestCpuThreadCap(unittest.TestCase):
    def test_explicit_override_wins(self):
        with mock.patch.dict(os.environ, {"CRISPRME_SCORER_THREADS": "8"}, clear=False):
            self.assertEqual(cb._cpu_thread_cap(), 8)

    def test_omp_num_threads_honored(self):
        env = {k: v for k, v in os.environ.items() if k not in ("CRISPRME_SCORER_THREADS",)}
        env["OMP_NUM_THREADS"] = "4"
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(cb._cpu_thread_cap(), 4)

    def test_default_caps_at_16(self):
        env = {k: v for k, v in os.environ.items()
               if k not in ("CRISPRME_SCORER_THREADS", "OMP_NUM_THREADS")}
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch("os.cpu_count", return_value=256):
            self.assertEqual(cb._cpu_thread_cap(), 16)

    def test_default_uses_cpu_count_when_small(self):
        env = {k: v for k, v in os.environ.items()
               if k not in ("CRISPRME_SCORER_THREADS", "OMP_NUM_THREADS")}
        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch("os.cpu_count", return_value=8):
            self.assertEqual(cb._cpu_thread_cap(), 8)

    def test_configure_tf_sets_thread_caps(self):
        cfg = TestConfigureTf._FakeConfig(gpus=[])
        cfg.experimental = mock.Mock()
        cfg.threading = mock.Mock()
        tf = mock.Mock(config=cfg)
        cb.configure_tf(tf, cb.CPU)
        cfg.threading.set_intra_op_parallelism_threads.assert_called()
        cfg.threading.set_inter_op_parallelism_threads.assert_called_with(2)


class TestArrayModule(unittest.TestCase):
    def test_metal_and_cpu_use_numpy(self):
        import numpy as np
        self.assertIs(cb.array_module(cb.METAL), np)
        self.assertIs(cb.array_module(cb.CPU), np)

    def test_cuda_without_cupy_falls_back_to_numpy(self):
        import numpy as np
        # cupy isn't installed on the dev Mac; the import fails -> numpy fallback
        self.assertIs(cb.array_module(cb.CUDA), np)


class TestConfigureTf(unittest.TestCase):
    class _FakeConfig:
        def __init__(self, gpus):
            self._gpus = gpus
            self.hidden = None
            self.grown = []

        def list_physical_devices(self, kind):
            return self._gpus if kind == "GPU" else []

        def set_visible_devices(self, devs, kind):
            self.hidden = (devs, kind)

        class experimental:  # noqa: N801
            pass

    def _tf(self, gpus):
        cfg = self._FakeConfig(gpus)
        cfg.experimental = mock.Mock()
        return mock.Mock(config=cfg), cfg

    def test_cpu_hides_gpus_when_present(self):
        tf, cfg = self._tf(gpus=["gpu0"])
        self.assertEqual(cb.configure_tf(tf, cb.CPU), cb.CPU)
        self.assertEqual(cfg.hidden, ([], "GPU"))

    def test_accelerator_with_no_gpu_degrades_to_cpu(self):
        tf, _ = self._tf(gpus=[])
        self.assertEqual(cb.configure_tf(tf, cb.METAL), cb.CPU)
        self.assertEqual(cb.configure_tf(tf, cb.CUDA), cb.CPU)

    def test_accelerator_with_gpu_keeps_backend(self):
        tf, cfg = self._tf(gpus=["gpu0"])
        self.assertEqual(cb.configure_tf(tf, cb.METAL), cb.METAL)
        cfg.experimental.set_memory_growth.assert_called()


if __name__ == "__main__":
    unittest.main()

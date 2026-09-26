"""Unit tests for scorer_runner — exercises the real subprocess protocol against a
FAKE worker (stdlib only, no TensorFlow), plus graceful-degradation paths."""
import os
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scorer_env  # noqa: E402
import scorer_runner  # noqa: E402

# A protocol-faithful worker with NO ML deps: ready handshake, then 1.0 if the pair
# matches else 0.25, honoring the quit command.
_GOOD_WORKER = textwrap.dedent(
    """
    import json, sys
    sys.stdout.write(json.dumps({"ready": True, "load_s": 0.0}) + "\\n"); sys.stdout.flush()
    while True:
        line = sys.stdin.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        if req.get("cmd") == "quit":
            break
        pairs = req.get("pairs", [])
        scores = [1.0 if p[0] == p[1] else 0.25 for p in pairs]
        sys.stdout.write(json.dumps({"scores": scores}) + "\\n"); sys.stdout.flush()
    """
)

# Returns the WRONG number of scores -> runner must disable + sentinel.
_BADCOUNT_WORKER = textwrap.dedent(
    """
    import json, sys
    sys.stdout.write(json.dumps({"ready": True, "load_s": 0.0}) + "\\n"); sys.stdout.flush()
    line = sys.stdin.readline()
    sys.stdout.write(json.dumps({"scores": [0.5]}) + "\\n"); sys.stdout.flush()
    """
)

# Never sends the ready handshake (exits immediately) -> spawn/handshake failure.
_NOREADY_WORKER = "import sys; sys.exit(0)\n"

# Reports ready:false (model-load failure fail-fast path).
_NOTREADY_WORKER = (
    'import json,sys\n'
    'sys.stdout.write(json.dumps({"ready": False, "error": "no weights"})+"\\n"); sys.stdout.flush()\n'
)

# Right count but non-numeric values -> must degrade (never raise).
_NONNUMERIC_WORKER = textwrap.dedent(
    """
    import json, sys
    sys.stdout.write(json.dumps({"ready": True, "load_s": 0.0}) + "\\n"); sys.stdout.flush()
    sys.stdin.readline()
    sys.stdout.write(json.dumps({"scores": ["x", "y"]}) + "\\n"); sys.stdout.flush()
    """
)

# NaN / out-of-range values -> must degrade.
_NAN_WORKER = textwrap.dedent(
    """
    import json, sys
    sys.stdout.write(json.dumps({"ready": True, "load_s": 0.0}) + "\\n"); sys.stdout.flush()
    sys.stdin.readline()
    sys.stdout.write(json.dumps({"scores": [float("nan"), 2.5]}) + "\\n"); sys.stdout.flush()
    """
)

# Per-row -1.0 sentinels are VALID (worker's own failure marker) -> preserved.
_SENTINEL_WORKER = textwrap.dedent(
    """
    import json, sys
    sys.stdout.write(json.dumps({"ready": True, "load_s": 0.0}) + "\\n"); sys.stdout.flush()
    while True:
        line = sys.stdin.readline()
        if not line: break
        req = json.loads(line)
        if req.get("cmd") == "quit": break
        n = len(req.get("pairs", []))
        sys.stdout.write(json.dumps({"scores": [-1.0] * n}) + "\\n"); sys.stdout.flush()
    """
)


def _write(tmpdir, name, body):
    p = os.path.join(tmpdir, name)
    with open(p, "w") as fh:
        fh.write(body)
    return p


class TestRunnerProtocol(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _runner(self, body, env_python=sys.executable):
        worker = _write(self.tmp, "fake_worker.py", body)
        r = scorer_runner.ScorerRunner(scorer="cbulge", device="cpu", worker=worker)
        self._patch = mock.patch.object(scorer_env, "env_python", return_value=env_python)
        self._patch.start()
        self.addCleanup(self._patch.stop)
        self.addCleanup(r.close)
        return r

    def test_happy_path(self):
        r = self._runner(_GOOD_WORKER)
        scores = r.predict(["AAAA", "GGGG"], ["AAAA", "CCCC"])
        self.assertEqual(scores, [1.0, 0.25])
        self.assertFalse(r.disabled)
        self.assertIsNotNone(r.load_s)

    def test_reuses_one_worker_across_batches(self):
        r = self._runner(_GOOD_WORKER)
        r.predict(["AA"], ["AA"])
        proc1 = r.proc
        r.predict(["GG"], ["CC"])
        self.assertIs(r.proc, proc1)  # persistent, not respawned

    def test_empty_batch(self):
        r = self._runner(_GOOD_WORKER)
        self.assertEqual(r.predict([], []), [])

    def test_length_mismatch_raises(self):
        r = self._runner(_GOOD_WORKER)
        with self.assertRaises(ValueError):
            r.predict(["AA", "GG"], ["AA"])

    def test_missing_env_degrades(self):
        r = self._runner(_GOOD_WORKER, env_python=None)
        scores = r.predict(["AA", "GG"], ["AA", "CC"])
        self.assertEqual(scores, [-1.0, -1.0])
        self.assertTrue(r.disabled)

    def test_bad_count_degrades(self):
        r = self._runner(_BADCOUNT_WORKER)
        scores = r.predict(["AA", "GG"], ["AA", "CC"])
        self.assertEqual(scores, [-1.0, -1.0])
        self.assertTrue(r.disabled)

    def test_no_ready_handshake_degrades(self):
        r = self._runner(_NOREADY_WORKER)
        scores = r.predict(["AA"], ["AA"])
        self.assertEqual(scores, [-1.0])
        self.assertTrue(r.disabled)

    def test_disabled_runner_is_sticky(self):
        r = self._runner(_NOREADY_WORKER)
        r.predict(["AA"], ["AA"])
        # second call must not attempt a respawn; still sentinel
        self.assertEqual(r.predict(["GG"], ["GG"]), [-1.0])

    def test_notready_handshake_degrades(self):
        r = self._runner(_NOTREADY_WORKER)
        self.assertEqual(r.predict(["AA"], ["AA"]), [-1.0])
        self.assertTrue(r.disabled)

    def test_nonnumeric_scores_degrade(self):
        r = self._runner(_NONNUMERIC_WORKER)
        self.assertEqual(r.predict(["AA", "GG"], ["AA", "CC"]), [-1.0, -1.0])
        self.assertTrue(r.disabled)

    def test_nan_or_out_of_range_degrades(self):
        r = self._runner(_NAN_WORKER)
        self.assertEqual(r.predict(["AA", "GG"], ["AA", "CC"]), [-1.0, -1.0])
        self.assertTrue(r.disabled)

    def test_per_row_sentinel_preserved(self):
        # -1.0 from the worker is a VALID per-row failure marker, not corruption
        r = self._runner(_SENTINEL_WORKER)
        self.assertEqual(r.predict(["AA", "GG"], ["AA", "CC"]), [-1.0, -1.0])
        self.assertFalse(r.disabled)

    def test_large_batch_is_chunked(self):
        r = self._runner(_GOOD_WORKER)
        with mock.patch.object(scorer_runner, "_MAX_BATCH", 2):
            sg = ["AA", "GG", "TT", "CC", "AA"]
            off = ["AA", "CC", "TT", "GG", "AA"]
            scores = r.predict(sg, off)
        self.assertEqual(scores, [1.0, 0.25, 1.0, 0.25, 1.0])


class TestModuleSingleton(unittest.TestCase):
    def test_predict_list_graceful_without_env(self):
        with mock.patch.object(scorer_env, "env_python", return_value=None):
            scorer_runner.close_runner()
            scores = scorer_runner.CRISPR_BULGE_predict_list(["AA", "GG"], ["AA", "CC"])
            self.assertEqual(scores, [-1.0, -1.0])
        scorer_runner.close_runner()


if __name__ == "__main__":
    unittest.main()

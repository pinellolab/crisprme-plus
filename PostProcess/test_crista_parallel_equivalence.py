"""CRISPRME_CRISTA_PARALLEL: the opt-in parallel CRISTA feature build must produce a matrix
element-for-element identical to the serial path (order + values), so the downstream RF
predict is byte-identical. Default OFF => serial, no executor created. The parallel path only
touches feature building (get_features); the 276MB RF model is never loaded here, so this test
needs no model file and runs anywhere the module imports."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import CRISTA_score as C  # noqa: E402


def _make_inputs(n):
    """Deterministic varied (guide, aligned_off, 29nt-window) triples (clean ACGT)."""
    import random
    rng = random.Random(0)  # seeds INPUT generation only (not any scoring RNG)
    b = "ACGT"
    return [("".join(rng.choice(b) for _ in range(20)),
             "".join(rng.choice(b) for _ in range(20)),
             "".join(rng.choice(b) for _ in range(29))) for _ in range(n)]


def _feat(inputs, workers, minbatch):
    os.environ["CRISPRME_CRISTA_PARALLEL"] = str(workers)
    os.environ["CRISPRME_CRISTA_PARALLEL_MINBATCH"] = str(minbatch)
    return C._build_crista_features(inputs)


def _key(feats):
    # stringify each element (feature rows mix floats + string numerals) for exact compare
    return [tuple(str(x) for x in row) for row in feats]


class TestCristaParallelEquivalence(unittest.TestCase):
    def tearDown(self):
        C._shutdown_crista_executor()
        os.environ.pop("CRISPRME_CRISTA_PARALLEL", None)
        os.environ.pop("CRISPRME_CRISTA_PARALLEL_MINBATCH", None)

    def test_serial_vs_parallel_identical(self):
        inp = _make_inputs(53)  # not divisible by the worker count
        serial = _key(_feat(inp, 1, 1))
        parallel = _key(_feat(inp, 4, 1))
        self.assertEqual(len(serial), 53)
        self.assertEqual(serial, parallel,
                         "parallel feature matrix must be element-wise + order identical")

    def test_chunk_boundary_uneven(self):
        inp = _make_inputs(37)  # 37 % 5 != 0 -> ragged final chunk
        self.assertEqual(_key(_feat(inp, 1, 1)), _key(_feat(inp, 5, 1)))

    def test_workers_gt_inputs(self):
        inp = _make_inputs(3)  # fewer inputs than workers -> nchunks clamps to n
        self.assertEqual(_key(_feat(inp, 8, 1)), _key(_feat(inp, 1, 1)))

    def test_minbatch_gate_stays_serial(self):
        # below minbatch -> serial path even with workers>1: no executor is created
        C._CRISTA_EXECUTOR = None
        _feat(_make_inputs(10), 4, 1000)  # n=10 < minbatch=1000
        self.assertIsNone(C._CRISTA_EXECUTOR, "small batch must stay serial (no executor)")

    def test_default_is_serial(self):
        os.environ.pop("CRISPRME_CRISTA_PARALLEL", None)
        os.environ.pop("CRISPRME_CRISTA_PARALLEL_MINBATCH", None)
        self.assertEqual(C._crista_parallel_cfg(), (1, 20000))  # default OFF
        C._CRISTA_EXECUTOR = None
        C._build_crista_features(_make_inputs(20))
        self.assertIsNone(C._CRISTA_EXECUTOR, "default must be serial (no executor)")


if __name__ == "__main__":
    unittest.main()

"""Fix-1 regression: the 276 MB CRISTA RandomForest model must be loaded ONCE per
process (module-level cache in predict_crista_score), and the reported CRISTA scores
must be byte-identical to the pre-cache behaviour.

Skips automatically when the CRISTA model isn't present (e.g. plain CI checkout with
no unpacked CRISTA_predictors.pkl/.zip), so it only runs where the model exists."""
import os
import pickle


def test_crista_model_loaded_once_and_scores_stable():
    here = os.path.dirname(os.path.abspath(__file__))
    if not (
        os.path.exists(os.path.join(here, "CRISTA_predictors.pkl"))
        or os.path.exists(os.path.join(here, "CRISTA_predictors.zip"))
    ):
        import pytest

        pytest.skip("CRISTA_predictors model not present in this checkout")

    import sys

    if here not in sys.path:
        sys.path.insert(0, here)

    _orig = pickle.load
    _crista_loads = [0]

    def _counted(f, *a, **k):
        try:
            if "CRISTA_predictors" in str(getattr(f, "name", "")):
                _crista_loads[0] += 1
        except Exception:
            pass
        return _orig(f, *a, **k)

    import CRISTA_score as C

    # start from a cold cache so this test controls the load count deterministically
    C._CRISTA_PREDICTORS = None
    pickle.load = _counted
    try:
        g = "TGCTTGGTCGGCACTGATAGNGG"
        off = "TGCTTGGTCGGCACTGATAGAGG"
        win = "AAAAATGCTTGGTCGGCACTGATAGGGGG"
        assert len(win) == 29
        s1 = ["{:.3f}".format(float(x)) for x in C.CRISTA_predict_list([g], [off], [win])]
        s2 = ["{:.3f}".format(float(x)) for x in C.CRISTA_predict_list([g], [off], [win])]
        s3 = ["{:.3f}".format(float(x)) for x in C.CRISTA_predict_list([g], [off], [win])]
    finally:
        pickle.load = _orig

    # the model is unpickled exactly once across three scoring calls
    assert _crista_loads[0] == 1, "CRISTA model loaded %d times, expected 1 (cache)" % _crista_loads[0]
    assert C._CRISTA_PREDICTORS is not None, "cache not populated"
    # reported (3-decimal) scores are stable across calls
    assert s1 == s2 == s3, "reported CRISTA scores unstable: %s / %s / %s" % (s1, s2, s3)

#!/usr/bin/env python
"""
CRISPR-Bulge off-target scorer — standalone prototype (Phase 1).

Mirrors the CRISPR-Bulge batch contract (`CRISPR_BULGE_score.CRISPR_BULGE_predict_list`) so it can
drop into `calculate_scores` behind the scorer-runner, EXCEPT it needs no 29-nt
genomic context (CRISPR-Bulge computes DNAshape over flanks; CRISPR-Bulge scores the
aligned pair alone):

    CRISPR_BULGE_predict_list(sgseq_aligned_list, offseq_aligned_list) -> list[float]

  * order preserved; each score is the 5-model ensemble mean probability in [0, 1]
  * a row that fails to encode/score yields -1.0 (never raises for one bad row)
  * models load ONCE (module-global cache); safe to call repeatedly per run/chrom

Runs inside the dedicated `cbulge` conda env (TensorFlow 2.x). The CRISPR-Bulge
source + weights are located via CBULGE_REPO (Phase 0 layout); Phase 4 will vendor
the 53 MB weights + fetch from Hugging Face on first use.

Device: CPU by default; set CRISPRME_COMPUTE_BACKEND=gpu|cuda|metal (or pass device=)
to use a visible GPU. Falls back to CPU with a warning if no GPU is present. A load-time
numerical self-test additionally guards against accelerators that MISCOMPUTE this model:
tensorflow-metal (Apple Silicon) silently returns garbage for the GRU kernel (every
score collapses to ~1.0), so if the selected accelerator fails to discriminate two
reference pairs, scoring transparently falls back to a CPU device context -- correct
scores, no crash, one warning (see load_models / _selftest_margin).

Reference: Yaish & Orenstein, NAR 2024 (doi:10.1093/nar/gkae428). MIT license.
"""
import os
import sys


def _repo():
    """Resolve the CRISPR-Bulge source+weights dir: $CBULGE_REPO, else the install
    location scorer_env provisions (pinned clone), else a dev fallback."""
    r = os.environ.get("CBULGE_REPO")
    if r:
        return r
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import scorer_env
        return scorer_env.default_cbulge_repo()
    except Exception:
        return "/srv/local/lp698/cbulge_bench/CRISPR-Bulge"


# 5-model classification ensemble (c_2, aligned, GUIDE-seq finetuned) — the config
# upstream's main_predict.py uses for the Refined_TrueOT ensemble prediction.
_ENSEMBLE_REL = (
    "files/bulges/1_folds/5_revision_ensemble_{}_exclude_RHAMPseq_continue_from_change_seq/"
    "read_ts_0/cleavage_models/aligned/FullGUIDEseq/classification/c_2/"
    "ln_x_plus_one_trans/model_fold_0"
)

_MODELS = None       # cached list of 5 loaded model instances
_PREDICT = None      # cached (predict fn, Padding_type, Encoding_type)
_BACKEND = None       # effective compute backend after load ('cpu'|'cuda'|'metal')
_TF = None            # cached tensorflow module (for the device-context fallback)
_DEVICE_CTX = None    # None = default placement; '/CPU:0' = forced-CPU fallback

# Load-time numerical self-test. tensorflow-metal (and, defensively, any future
# accelerator build) can silently miscompute this model's GRU kernel and collapse
# every score to ~1.0 -- catastrophic for an off-target tool (every site would look
# like a perfect cut). So after loading on an accelerator we score two reference
# pairs whose CORRECT scores are far apart; if the backend fails to discriminate them
# (or emits NaN), we transparently force scoring onto the CPU device (the model
# weights are fine -- only op PLACEMENT changes -- so no reload is needed).
_SELFTEST_SG = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTTGCTCCACAGG"]
_SELFTEST_OFF = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTT-CTCCACAGG"]
_SELFTEST_MIN_MARGIN = 0.5   # perfect-match minus bulge-disrupted; CPU truth ~1.0


def _validated_repo():
    """Resolve+validate CBULGE_REPO before it goes on sys.path.

    CBULGE_REPO is trusted config, but we still refuse to inject a path that isn't a
    real CRISPR-Bulge checkout: this both hardens against a bogus/hostile value and
    turns a misconfiguration into a clear error instead of a confusing ImportError.
    """
    src = _repo()
    repo = os.path.realpath(os.path.abspath(src))
    if not os.path.isdir(repo) or not os.path.isdir(os.path.join(repo, "OT_deep_score_src")):
        raise RuntimeError(
            f"CRISPR-Bulge source not found (missing OT_deep_score_src): {src!r} "
            f"-- provision with: crisprme.py scorer-env create"
        )
    return repo


def _ensemble_paths():
    return [os.path.join(_validated_repo(), _ENSEMBLE_REL.format(i)) for i in range(5)]


def _score_df(df):
    """Run the 5-model ensemble over a prepared DataFrame under the active device
    context (_DEVICE_CTX). Returns a numpy array of ensemble-mean probabilities."""
    import contextlib
    import numpy as np
    predict, Padding_type, Encoding_type = _PREDICT
    cm = _TF.device(_DEVICE_CTX) if (_DEVICE_CTX and _TF is not None) else contextlib.nullcontext()
    preds = []
    with cm:
        for m in _MODELS:
            y = predict(
                df.copy(), model=m,
                include_distance_feature=False, include_sequence_features=True,
                include_gmt_score=False, include_nuclea_seq_score=False,
                padding_type=Padding_type.GAP, aligned=True, bulges=True,
                encoding_type=Encoding_type.ONE_HOT, flat_encoding=False,
            )
            preds.append(np.asarray(y).ravel())
    return np.mean(preds, axis=0)


def _selftest_margin():
    """Score the reference pairs under the CURRENT device context and return
    (margin, has_nan) where margin = perfect-match score - bulge-disrupted score.
    A correct backend gives margin ~1.0; a miscomputing one collapses it toward 0."""
    import numpy as np
    import pandas as pd
    from OT_deep_score_src.general_utilities import OFF_TARGET, SG_RNA_SEQ
    df = pd.DataFrame({SG_RNA_SEQ: _SELFTEST_SG, OFF_TARGET: _SELFTEST_OFF})
    s = _score_df(df)
    has_nan = bool(np.any(~np.isfinite(s)))
    return (float(s[0] - s[1]), has_nan)


def load_models(device=None):
    """Load + cache the ensemble on the resolved compute backend (cuda|metal|cpu).
    Idempotent. Returns the cached model list. The effective backend is recorded in
    the module global _BACKEND; if an accelerator fails the numerical self-test we set
    _DEVICE_CTX='/CPU:0' so scoring runs correctly on the CPU (see module notes)."""
    global _MODELS, _PREDICT, _BACKEND, _TF, _DEVICE_CTX
    if _MODELS is not None:
        return _MODELS
    # compute_backend lives beside this file (PostProcess); ensure it's importable
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import compute_backend as cb

    backend = cb.resolve_backend(device)          # cuda | metal | cpu
    cb.pre_import_env(backend)                     # set CUDA_VISIBLE_DEVICES before TF (cpu/cuda)
    repo = _validated_repo()
    if repo not in sys.path:
        sys.path.insert(0, repo)
    import tensorflow as tf
    _TF = tf
    _BACKEND = cb.configure_tf(tf, backend)        # device placement; effective backend
    from OT_deep_score_src.general_utilities import Encoding_type, Padding_type
    from OT_deep_score_src.models_inter import Model
    from OT_deep_score_src.predict_utilities import predict

    _MODELS = [Model.load_model_instance(p) for p in _ensemble_paths()]
    _PREDICT = (predict, Padding_type, Encoding_type)

    # Numerical self-test: only distrust actual accelerators (cpu is the oracle).
    if _BACKEND != cb.CPU:
        try:
            margin, has_nan = _selftest_margin()
            if has_nan or margin < _SELFTEST_MIN_MARGIN:
                _DEVICE_CTX = "/CPU:0"
                margin2, has_nan2 = _selftest_margin()  # confirm CPU fallback is sane
                ok = (not has_nan2) and margin2 >= _SELFTEST_MIN_MARGIN
                sys.stderr.write(
                    f"[crispr_bulge] WARNING: {_BACKEND} miscomputes this model "
                    f"(self-test margin={margin:.4g}, nan={has_nan}); forcing CPU "
                    f"device for scoring (CPU margin={margin2:.4g}). "
                    f"GPU acceleration is unavailable for the scorer on this stack.\n"
                )
                _BACKEND = f"{_BACKEND}->cpu" if ok else _BACKEND
        except Exception as e:
            sys.stderr.write(f"[crispr_bulge] self-test skipped ({e})\n")

    sys.stderr.write(f"[crispr_bulge] backend={_BACKEND} (requested {backend})\n")
    return _MODELS


def CRISPR_BULGE_predict_list(sgseq_aligned_list, offseq_aligned_list, device=None):
    """Score aligned (sgRNA, off-target) pairs. See module docstring for the contract."""
    n = len(sgseq_aligned_list)
    if n != len(offseq_aligned_list):
        raise ValueError("sgseq and offseq lists differ in length")
    if n == 0:
        return []

    try:
        import pandas as pd

        # load-once (raises on missing env/weights/bad CBULGE_REPO) -> degrade to -1.
        # MUST precede the OT_deep_score_src import: load_models() puts the repo on
        # sys.path, so importing its symbols before this would fail with ModuleNotFound.
        load_models(device)
        from OT_deep_score_src.general_utilities import OFF_TARGET, SG_RNA_SEQ

        df = pd.DataFrame(
            {SG_RNA_SEQ: [s.upper() for s in sgseq_aligned_list],
             OFF_TARGET: [o.upper() for o in offseq_aligned_list]}
        )
        return _score_df(df).astype(float).tolist()
    except Exception as e:  # never let one bad batch kill scoring
        sys.stderr.write(f"[crispr_bulge] batch scoring failed: {e}\n")
        return [-1.0] * n


if __name__ == "__main__":
    # tiny self-check
    sg = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTTCTCCACAGG"]
    off = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTTCTCCTCAGG"]
    print(CRISPR_BULGE_predict_list(sg, off))

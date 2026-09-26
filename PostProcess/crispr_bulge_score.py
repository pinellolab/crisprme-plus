#!/usr/bin/env python
"""
CRISPR-Bulge off-target scorer — standalone prototype (Phase 1).

Mirrors the CRISTA batch contract (`CRISTA_score.CRISTA_predict_list`) so it can
drop into `calculate_scores` behind the scorer-runner, EXCEPT it needs no 29-nt
genomic context (CRISTA computes DNAshape over flanks; CRISPR-Bulge scores the
aligned pair alone):

    CRISPR_BULGE_predict_list(sgseq_aligned_list, offseq_aligned_list) -> list[float]

  * order preserved; each score is the 5-model ensemble mean probability in [0, 1]
  * a row that fails to encode/score yields -1.0 (never raises for one bad row)
  * models load ONCE (module-global cache); safe to call repeatedly per run/chrom

Runs inside the dedicated `cbulge` conda env (TensorFlow 2.x). The CRISPR-Bulge
source + weights are located via CBULGE_REPO (Phase 0 layout); Phase 4 will vendor
the 53 MB weights + fetch from Hugging Face on first use.

Device: CPU by default; set CRISPRME_COMPUTE_BACKEND=gpu (or pass device="gpu") to
use a visible GPU. Falls back to CPU with a warning if no GPU is present.

Reference: Yaish & Orenstein, NAR 2024 (doi:10.1093/nar/gkae428). MIT license.
"""
import os
import sys

REPO = os.environ.get("CBULGE_REPO", "/srv/local/lp698/cbulge_bench/CRISPR-Bulge")

# 5-model classification ensemble (c_2, aligned, GUIDE-seq finetuned) — the config
# upstream's main_predict.py uses for the Refined_TrueOT ensemble prediction.
_ENSEMBLE_REL = (
    "files/bulges/1_folds/5_revision_ensemble_{}_exclude_RHAMPseq_continue_from_change_seq/"
    "read_ts_0/cleavage_models/aligned/FullGUIDEseq/classification/c_2/"
    "ln_x_plus_one_trans/model_fold_0"
)

_MODELS = None       # cached list of 5 loaded model instances
_PREDICT = None      # cached (predict fn, Padding_type, Encoding_type)


def _ensemble_paths():
    return [os.path.join(REPO, _ENSEMBLE_REL.format(i)) for i in range(5)]


def _select_device(device):
    """Resolve requested device -> set CUDA_VISIBLE_DEVICES BEFORE importing TF."""
    dev = (device or os.environ.get("CRISPRME_COMPUTE_BACKEND", "cpu")).lower()
    if dev == "gpu":
        # leave CUDA_VISIBLE_DEVICES as-is (let TF see the GPU)
        return "gpu"
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    return "cpu"


def load_models(device=None):
    """Load + cache the ensemble. Idempotent. Returns the cached model list."""
    global _MODELS, _PREDICT
    if _MODELS is not None:
        return _MODELS
    requested = _select_device(device)
    if REPO not in sys.path:
        sys.path.insert(0, REPO)
    import tensorflow as tf  # noqa: F401  (import after CUDA_VISIBLE_DEVICES is set)
    from OT_deep_score_src.general_utilities import Encoding_type, Padding_type
    from OT_deep_score_src.models_inter import Model
    from OT_deep_score_src.predict_utilities import predict

    if requested == "gpu" and not tf.config.list_physical_devices("GPU"):
        sys.stderr.write("[crispr_bulge] GPU requested but none visible; using CPU.\n")

    _MODELS = [Model.load_model_instance(p) for p in _ensemble_paths()]
    _PREDICT = (predict, Padding_type, Encoding_type)
    return _MODELS


def CRISPR_BULGE_predict_list(sgseq_aligned_list, offseq_aligned_list, device=None):
    """Score aligned (sgRNA, off-target) pairs. See module docstring for the contract."""
    n = len(sgseq_aligned_list)
    if n != len(offseq_aligned_list):
        raise ValueError("sgseq and offseq lists differ in length")
    if n == 0:
        return []

    import numpy as np
    import pandas as pd
    from OT_deep_score_src.general_utilities import OFF_TARGET, SG_RNA_SEQ

    models = load_models(device)
    predict, Padding_type, Encoding_type = _PREDICT

    df = pd.DataFrame(
        {SG_RNA_SEQ: [s.upper() for s in sgseq_aligned_list],
         OFF_TARGET: [o.upper() for o in offseq_aligned_list]}
    )
    try:
        # exact upstream ensemble_predict kwargs -> guaranteed fidelity
        preds = []
        for m in models:
            y = predict(
                df.copy(), model=m,
                include_distance_feature=False, include_sequence_features=True,
                include_gmt_score=False, include_nuclea_seq_score=False,
                padding_type=Padding_type.GAP, aligned=True, bulges=True,
                encoding_type=Encoding_type.ONE_HOT, flat_encoding=False,
            )
            preds.append(np.asarray(y).ravel())
        return np.mean(preds, axis=0).astype(float).tolist()
    except Exception as e:  # never let one bad batch kill scoring
        sys.stderr.write(f"[crispr_bulge] batch scoring failed: {e}\n")
        return [-1.0] * n


if __name__ == "__main__":
    # tiny self-check
    sg = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTTCTCCACAGG"]
    off = ["GTAACGGCAGACTTCTCCACAGG", "GTAACGGCAGACTTCTCCTCAGG"]
    print(CRISPR_BULGE_predict_list(sg, off))

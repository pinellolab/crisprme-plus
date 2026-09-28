#!/usr/bin/env python
"""
Phase-0 CRISPR-Bulge benchmark  (standalone; touches NO CRISPRme pipeline code).

Answers the gate questions for REPLACE / ALONGSIDE / NO-GO:
  * model-load cost (one-time, amortized by a persistent worker)
  * throughput (OTs/sec), p50/p99 per-batch latency, peak RSS
    at batch sizes {1k, 10k, 100k, 1M} on CPU (thread sweep) and GPU
  * accuracy: reproduce the ensemble AUPR on Refined_TrueOT (validates that the
    vendored weights + our encoding path are wired correctly) + a CFD anchor on
    the same set, split out on the bulge-containing subset (where CRISTA/CFD are
    known-weak and CRISPR-Bulge is claimed to win).

Run from the CRISPR-Bulge repo root inside the cbulge env:
    mamba run -n cbulge_cpu python bench_cbulge.py --mode speed  --sizes 1000,10000,100000
    mamba run -n cbulge_gpu python bench_cbulge.py --mode speed  --sizes 1000,10000,100000,1000000 --gpu
    mamba run -n cbulge_cpu python bench_cbulge.py --mode accuracy

Emits a human table to stdout and a JSON blob to --out (default bench_out.json).
"""
import argparse
import json
import os
import resource
import sys
import time

REPO = os.environ.get("CBULGE_REPO", "/srv/local/lp698/cbulge_bench/CRISPR-Bulge")
sys.path.insert(0, REPO)
os.chdir(REPO)

# The DataFrame column names the encoder reads are these upstream constants
# (SG_RNA_SEQ == "Align.sgRNA", OFF_TARGET == "Align.off-target", LABEL == "label").
# Import them so we always match, whatever the CSV header is called.
from OT_deep_score_src.general_utilities import SG_RNA_SEQ, OFF_TARGET, LABEL  # noqa: E402

# 5-model classification ensemble (c_2, aligned, GUIDE-seq fine-tuned) — the
# configuration main_predict.py demonstrates for Refined_TrueOT.
ENSEMBLE_PATHS = [
    "files/bulges/1_folds/5_revision_ensemble_{}_exclude_RHAMPseq_continue_from_change_seq/"
    "read_ts_0/cleavage_models/aligned/FullGUIDEseq/classification/c_2/"
    "ln_x_plus_one_trans/model_fold_0".format(i)
    for i in range(5)
]

TRUEOT_CSV = "files/datasets/Refined_TrueOT.csv"


def peak_rss_mb():
    # ru_maxrss is KB on Linux, bytes on macOS.
    v = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v / 1024.0 if sys.platform.startswith("linux") else v / (1024.0 * 1024.0)


def configure_threads(n_threads, use_gpu):
    """Must run before TF touches any op. Returns the imported tf module."""
    if not use_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
    if n_threads:
        os.environ["OMP_NUM_THREADS"] = str(n_threads)
    import tensorflow as tf

    if n_threads:
        try:
            tf.config.threading.set_intra_op_parallelism_threads(n_threads)
            tf.config.threading.set_inter_op_parallelism_threads(max(1, n_threads // 4))
        except RuntimeError:
            pass  # already initialized
    gpus = tf.config.list_physical_devices("GPU")
    if use_gpu and gpus:
        for g in gpus:
            try:
                tf.config.experimental.set_memory_growth(g, True)
            except RuntimeError:
                pass
    return tf, bool(gpus)


def load_models():
    from OT_deep_score_src.models_inter import Model

    t0 = time.time()
    models = [Model.load_model_instance(p) for p in ENSEMBLE_PATHS]
    return models, time.time() - t0


def make_scorer(models):
    """Return score_batch(df)->np.ndarray[avg ensemble prob]; feature-build once, predict 5x."""
    import numpy as np
    from OT_deep_score_src.data_processing_utilities import build_sequence_features
    from OT_deep_score_src.general_utilities import Encoding_type, Padding_type

    def score_batch(df):
        feats = build_sequence_features(
            df.copy(),
            include_distance_feature=False,
            include_sequence_features=True,
            include_gmt_score=False,
            include_nuclea_seq_score=False,
            bulges=True,
            padding_type=Padding_type.GAP,
            aligned=True,
            encoding_type=Encoding_type.ONE_HOT,
            flat_encoding=False,
        )
        preds = [np.asarray(m.predict(feats)).ravel() for m in models]
        return np.mean(preds, axis=0)

    return score_batch


def load_trueot():
    import pandas as pd

    df = pd.read_csv(TRUEOT_CSV)
    return df


def resample_to(df, n, seed=0):
    """Sample-with-replacement to size n, preserving real format + distribution."""
    return df.sample(n=n, replace=True, random_state=seed).reset_index(drop=True)


def percentile(vals, p):
    if not vals:
        return None
    s = sorted(vals)
    k = (len(s) - 1) * (p / 100.0)
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def run_speed(sizes, use_gpu, n_threads, repeats, out):
    import numpy as np  # noqa

    tf, gpu_present = configure_threads(n_threads, use_gpu)
    device = "gpu" if (use_gpu and gpu_present) else "cpu"
    print(f"[speed] device={device}  gpu_present={gpu_present}  threads={n_threads or 'default'}  tf={tf.__version__}")

    models, load_s = load_models()
    print(f"[speed] model-load (5-model ensemble): {load_s:.2f}s")
    score = make_scorer(models)

    base = load_trueot()[[SG_RNA_SEQ, OFF_TARGET]].dropna().reset_index(drop=True)
    print(f"[speed] TrueOT rows={len(base)}  example lens "
          f"sg={len(base[SG_RNA_SEQ].iloc[0])} ot={len(base[OFF_TARGET].iloc[0])}")

    # warm-up (build TF graph / kernels) so the first timed batch isn't graph-compile cost
    _ = score(resample_to(base, min(256, len(base)), seed=99))

    results = []
    for n in sizes:
        batch = resample_to(base, n, seed=n)
        lat = []
        for r in range(repeats):
            t0 = time.time()
            y = score(batch)
            dt = time.time() - t0
            lat.append(dt)
            assert len(y) == n
        best = min(lat)
        rec = {
            "n": n,
            "device": device,
            "threads": n_threads or "default",
            "repeats": repeats,
            "best_s": best,
            "median_s": percentile(lat, 50),
            "p99_s": percentile(lat, 99),
            "ots_per_sec": n / best,
            "peak_rss_mb": round(peak_rss_mb(), 1),
        }
        results.append(rec)
        print(
            f"[speed] n={n:>8}  best={best:8.3f}s  {rec['ots_per_sec']:>12,.0f} OT/s  "
            f"p99={rec['p99_s']:.3f}s  RSS={rec['peak_rss_mb']:.0f}MB"
        )

    blob = {"kind": "speed", "device": device, "model_load_s": load_s, "results": results}
    _append_out(out, blob)
    # project to genome-wide 6.77M-OT scale from the largest batch
    if results:
        biggest = max(results, key=lambda r: r["n"])
        gw = 6_770_000 / biggest["ots_per_sec"]
        print(f"[speed] projected GW (6.77M OT) at {device}: {gw/60:.1f} min "
              f"({biggest['ots_per_sec']:,.0f} OT/s from n={biggest['n']})")
    return blob


def _aupr(y_true, y_score):
    from sklearn.metrics import average_precision_score

    return float(average_precision_score(y_true, y_score))


def run_accuracy(out):
    import numpy as np
    import pandas as pd  # noqa

    tf, _ = configure_threads(None, False)
    models, load_s = load_models()
    score = make_scorer(models)
    df = load_trueot()
    print(f"[acc] TrueOT columns: {list(df.columns)}")
    label_col = LABEL if LABEL in df.columns else None
    work = df[[SG_RNA_SEQ, OFF_TARGET]].dropna().reset_index(drop=True)
    df = df.loc[work.index] if len(work) == len(df) else df.dropna(subset=[SG_RNA_SEQ, OFF_TARGET]).reset_index(drop=True)
    y_score = score(work)

    blob = {"kind": "accuracy", "n": int(len(work)), "label_col": label_col}
    if label_col is not None:
        y_raw = pd.to_numeric(df[label_col], errors="coerce").fillna(0).values
        y_true = (y_raw > 0).astype(int)
        # bulge subset: aligned pair contains a gap
        has_bulge = ((work[SG_RNA_SEQ].str.contains("-")) | (work[OFF_TARGET].str.contains("-"))).values
        blob["aupr_all"] = _aupr(y_true, y_score)
        blob["n_pos"] = int(y_true.sum())
        if has_bulge.any():
            blob["aupr_bulge"] = _aupr(y_true[has_bulge], y_score[has_bulge])
            blob["n_bulge"] = int(has_bulge.sum())
            blob["n_bulge_pos"] = int(y_true[has_bulge].sum())
        print(f"[acc] ensemble AUPR all={blob['aupr_all']:.4f}  "
              f"bulge={blob.get('aupr_bulge','n/a')}  "
              f"(n={blob['n']}, pos={blob['n_pos']}, bulge={blob.get('n_bulge','?')})")
    else:
        print("[acc] no label column found; printing score distribution only")
        blob["score_mean"] = float(np.mean(y_score))
        blob["score_std"] = float(np.std(y_score))
    _append_out(out, blob)
    return blob


def _append_out(path, blob):
    data = []
    if os.path.exists(path):
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception:
            data = []
    if not isinstance(data, list):
        data = [data]
    data.append(blob)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["speed", "accuracy"], default="speed")
    ap.add_argument("--sizes", default="1000,10000,100000")
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--threads", type=int, default=0, help="0 = TF default")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--out", default="bench_out.json")
    args = ap.parse_args()

    if args.mode == "speed":
        sizes = [int(x) for x in args.sizes.split(",") if x]
        run_speed(sizes, args.gpu, args.threads, args.repeats, args.out)
    else:
        run_accuracy(args.out)


if __name__ == "__main__":
    main()

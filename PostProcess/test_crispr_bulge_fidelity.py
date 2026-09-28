#!/usr/bin/env python
"""
Phase-1 fidelity gate: our CRISPR_BULGE_predict_list must reproduce the upstream
ensemble prediction (`pred_averege_ensemble`) bit-for-bit (within float tol).

Runs in the cbulge env from the CRISPR-Bulge repo root:
    CBULGE_REPO=/srv/local/lp698/cbulge_bench/CRISPR-Bulge \
      mamba run -n cbulge_cpu python test_crispr_bulge_fidelity.py
Exit 0 = match. Non-zero = drift.
"""
import os
import sys

REPO = os.environ.get("CBULGE_REPO", "/srv/local/lp698/cbulge_bench/CRISPR-Bulge")
sys.path.insert(0, REPO)
os.chdir(REPO)
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

N = int(os.environ.get("FIDELITY_N", "500"))


def main():
    import numpy as np
    import pandas as pd
    from OT_deep_score_src.general_utilities import OFF_TARGET, SG_RNA_SEQ
    from train_and_predict_scripts.utilities import ensemble_predict

    import crispr_bulge_score as cb

    df = pd.read_csv("files/datasets/Refined_TrueOT.csv").dropna(
        subset=[SG_RNA_SEQ, OFF_TARGET]
    ).head(N).reset_index(drop=True)
    # include a healthy number of bulge rows in the sample
    bulge = df[df[OFF_TARGET].str.contains("-") | df[SG_RNA_SEQ].str.contains("-")]
    sample = pd.concat([df.head(N // 2), bulge.head(N // 2)]).drop_duplicates().reset_index(drop=True)
    print(f"[fidelity] sample n={len(sample)} (bulge rows in sample: "
          f"{int((sample[OFF_TARGET].str.contains('-') | sample[SG_RNA_SEQ].str.contains('-')).sum())})")

    # upstream reference (writes predictions.csv as a side effect; that's fine)
    ref = ensemble_predict(cb._ensemble_paths(), sample.copy())
    ref_scores = ref["pred_averege_ensemble"].to_numpy()

    # our module
    ours = np.asarray(cb.CRISPR_BULGE_predict_list(
        sample[SG_RNA_SEQ].tolist(), sample[OFF_TARGET].tolist()))

    max_abs = float(np.max(np.abs(ref_scores - ours)))
    ok = np.allclose(ref_scores, ours, atol=1e-6, rtol=0)
    print(f"[fidelity] max|Δ| = {max_abs:.3e}  ->  {'MATCH' if ok else 'DRIFT'}")
    print(f"[fidelity] examples ref/ours: "
          f"{list(zip(np.round(ref_scores[:4],5), np.round(ours[:4],5)))}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

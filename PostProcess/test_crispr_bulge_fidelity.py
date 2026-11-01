#!/usr/bin/env python
"""
Phase-1 fidelity gate: our CRISPR_BULGE_predict_list must reproduce the upstream
ensemble prediction (`pred_averege_ensemble`) bit-for-bit (within float tol).

Runs in the cbulge env, against a CRISPR-Bulge checkout that holds the Refined_TrueOT
dataset (CBULGE_REPO, an absolute path; default: the copy `crisprme.py scorer-env
create` provisions, <prefix>/opt/CRISPR-Bulge):
    CBULGE_REPO=/path/to/CRISPR-Bulge \
      mamba run -n cbulge python test_crispr_bulge_fidelity.py
(micromamba run in the Docker image). Exit 0 = match. 1 = drift. 2 = it cannot run
here (no TensorFlow, or no dataset).

The provisioned copy holds no datasets, and no git metadata to fetch them with: it is
cloned with Git LFS smudge off, so its files/datasets.zip is only an LFS pointer, and
its .git is then deleted. Fetch the one file this reads in a clone of your own at the
same commit, with git-lfs installed (a 524 MB download):
    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/OrensteinLab/CRISPR-Bulge.git cb
    git -C cb checkout 3eddcd5bfcaff00b2bdf29425116ec9756ace870
    git -C cb lfs install --local
    git -C cb lfs pull --include files/datasets.zip
    unzip cb/files/datasets.zip 'datasets/Refined_TrueOT.csv' -d cb/files
then run with CBULGE_REPO="$PWD/cb". Or unzip into the provisioned copy instead
(-d <prefix>/opt/CRISPR-Bulge/files) and run with no CBULGE_REPO. Where `git lfs
install` was never run, `git lfs pull` without the `lfs install` step above leaves the
pointer in place and still exits 0.

Under `python -m unittest` (CI, the main env, the image) the gate is one test, skipped
with the reason wherever TensorFlow or the dataset is absent; a second checks that the
steps above clone the repository and commit scorer-env provisions. Nothing runs at
import.
"""
import os
import sys
import unittest

N = int(os.environ.get("FIDELITY_N", "500"))
DATASET = os.path.join("files", "datasets", "Refined_TrueOT.csv")
HERE = os.path.dirname(os.path.abspath(__file__))


def cbulge_repo():
    """The CRISPR-Bulge checkout: $CBULGE_REPO, else where scorer_env provisions it."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import scorer_env
    return scorer_env.default_cbulge_repo()


def cannot_run(repo):
    """Why the gate cannot run here, or None."""
    import importlib.util
    if importlib.util.find_spec("tensorflow") is None:
        return "no TensorFlow here: run it in the cbulge env"
    if not os.path.isfile(os.path.join(repo, DATASET)):
        return f"no {DATASET} in {repo} (this module's docstring says how to fetch it)"
    return None


def fidelity(repo):
    """Score a sample of Refined_TrueOT with upstream and with ours; True on a match.
    Runs from the repo root, as upstream's code expects."""
    sys.path.insert(0, repo)
    os.chdir(repo)
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

    import numpy as np
    import pandas as pd
    from OT_deep_score_src.general_utilities import OFF_TARGET, SG_RNA_SEQ
    from train_and_predict_scripts.utilities import ensemble_predict

    import crispr_bulge_score as cb

    df = pd.read_csv(DATASET).dropna(
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
    return bool(ok)


class CrisprBulgeFidelity(unittest.TestCase):
    def test_our_scorer_reproduces_the_upstream_ensemble(self):
        repo = cbulge_repo()
        why = cannot_run(repo)
        if why:
            self.skipTest(why)
        self.addCleanup(os.chdir, os.getcwd())
        self.assertTrue(fidelity(repo), "CRISPR_BULGE_predict_list drifted from upstream")


@unittest.skipIf(os.environ.get("CRISPRME_CBULGE_URL") or os.environ.get("CRISPRME_CBULGE_PIN"),
                 "the CRISPR-Bulge source scorer-env provisions is overridden here")
class TheDocumentedFetch(unittest.TestCase):
    """The docstring's steps fetch the dataset into a clone of the repository and commit
    that scorer-env provisions, so a bump of the pin cannot leave them on the old one."""

    def test_clones_the_provisioned_repository_at_the_provisioned_commit(self):
        if HERE not in sys.path:
            sys.path.insert(0, HERE)
        import scorer_env
        doc = sys.modules[__name__].__doc__
        self.assertIn(f"git clone {scorer_env.CBULGE_URL} cb\n", doc)
        self.assertIn(f"git -C cb checkout {scorer_env.CBULGE_PIN}\n", doc)


def main():
    repo = cbulge_repo()
    why = cannot_run(repo)
    if why:
        print(f"[fidelity] cannot run: {why}", file=sys.stderr)
        sys.exit(2)
    sys.exit(0 if fidelity(repo) else 1)


if __name__ == "__main__":
    main()

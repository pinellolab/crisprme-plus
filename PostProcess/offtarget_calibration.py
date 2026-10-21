#!/usr/bin/env python3
"""Off-target probability calibration — runtime application.

Loads a pre-fit calibration artifact (``<scorer>_calibration.npz`` from
``scripts/fit_calibration.py``) and turns a search's raw per-site scores into,
per score band, the EXPECTED number of true off-targets plus a coverage-valid
90% prediction interval.

Design (locked after an adversarial statistician review + a bake-off):
  * point estimate = isotonic score->P(true off-target), monotone, leave-guide-out
    validated; fair smooth challengers under-count the high-risk decision band, so
    isotonic is kept.
  * interval = two-level guide-cluster bootstrap over the bundled ensemble of
    guide-resampled isotonic maps, PLUS a negative-binomial / Gamma between-guide
    overdispersion multiplier (``phi_bands``), PLUS a Poisson count draw. Validated
    per-held-out-guide (= one novel guide) to ~90% coverage.
  * EACH reported band -- including aggregates such as >=0.2 -- is computed on its
    OWN pooled sites with its OWN overdispersion ``phi_bands[label]``; aggregates are
    NEVER formed by summing independent sub-band draws (a guide's promiscuity is
    correlated across bands ~+0.66, so independent summation wrongly inflates the
    aggregate's lower bound).
  * both scorers are sequence-based / variant-agnostic, so the calibration applies to
    any <=1-gap alignment (reference or variant/ALT allele). Only >=2-gap sites and the
    per-sample carrier-incidence dimension are out of scope (see METHODS).
  * counts are at the benchmark prevalence; callers present a prevalence-sensitivity
    (``prevalence_scaled``) as a separate deployment axis.

numpy only at runtime.
"""
from __future__ import annotations
import json
import os
from typing import Dict, List, Sequence

import numpy as np

_ARTIFACT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration")
_CACHE: Dict[str, dict] = {}


def set_artifact_dir(path: str) -> None:
    """Override the artifact directory (e.g. for out-of-tree testing)."""
    global _ARTIFACT_DIR
    _ARTIFACT_DIR = path
    _CACHE.clear()


def artifact_path(scorer: str) -> str:
    return os.path.join(_ARTIFACT_DIR, f"{scorer}_calibration.npz")


def load_artifact(scorer: str) -> dict:
    """Load (and cache) the calibration artifact for ``scorer``."""
    if scorer in _CACHE:
        return _CACHE[scorer]
    path = artifact_path(scorer)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"calibration artifact not found: {path}")
    npz = np.load(path)
    art = {
        "grid": npz["grid"],
        "p_point": npz["p_point"],
        "p_ens": npz["p_ens"],
        "meta": json.loads(str(npz["meta"])),
    }
    _CACHE[scorer] = art
    return art


def _in_range(scores: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return (scores >= lo) & (scores < hi) if hi < 1.0 else (scores >= lo)


def expected_counts(
    scores: Sequence[float],
    scorer: str,
    n_boot: int = 4000,
    seed: int = 0,
) -> List[dict]:
    """Per-report-band expected true-off-target count + coverage-valid 90% interval.

    ``scores`` are the raw per-site scores of ONE search, restricted by the caller
    to the calibration support (<=1-gap off-targets; on-targets excluded). Returns one
    dict per report band in ``meta['report_bands']``:
      ``band`` (label), ``lo``/``hi`` (score range), ``n_sites``, ``point`` (expected
      count), ``lo90``/``hi90`` (90% prediction interval for the realized count).
    """
    art = load_artifact(scorer)
    grid, p_point, ens = art["grid"], art["p_point"], art["p_ens"]
    meta = art["meta"]
    report_bands = meta["report_bands"]
    phi_bands = meta.get("phi_bands", {})
    s = np.asarray(scores, dtype=float)
    rng = np.random.default_rng(seed)
    # sample n_boot ensemble maps WITH replacement (do NOT cap at the ensemble size -- that
    # pins the interval to the number of bundled maps and leaves the heavy-tailed upper
    # bound under-resolved; each draw adds a fresh Gamma + Poisson realization)
    ens_idx = rng.integers(0, ens.shape[0], n_boot)

    out = []
    for label, (lo, hi) in report_bands.items():
        m = _in_range(s, float(lo), float(hi))
        sb = s[m]
        point = float(np.interp(sb, grid, p_point).sum()) if sb.size else 0.0
        phi = float(phi_bands.get(label, 0.0))
        if sb.size == 0:
            out.append({"band": label, "lo": lo, "hi": hi, "n_sites": 0,
                        "point": 0.0, "lo90": 0.0, "hi90": 0.0})
            continue
        draws = np.empty(len(ens_idx))
        for j, ei in enumerate(ens_idx):
            lam = float(np.interp(sb, grid, ens[ei]).sum())
            if phi > 1e-9 and lam > 0:                     # between-guide overdispersion
                lam = rng.gamma(1.0 / phi, lam * phi)
            draws[j] = rng.poisson(max(lam, 0.0))
        lo90, hi90 = np.percentile(draws, [5, 95])
        out.append({"band": label, "lo": lo, "hi": hi, "n_sites": int(m.sum()),
                    "point": point, "lo90": float(lo90), "hi90": float(hi90)})
    return out


def prevalence_scaled(point: float, factors=(1 / 3.0, 3.0)) -> Dict[str, float]:
    """Prevalence-sensitivity: the absolute count scales ~linearly with the guide's
    own off-target propensity relative to the benchmark average (a separate deployment
    axis, not a coverage interval)."""
    return {f"x{f:.2g}": point * f for f in factors}


def artifact_summary(scorer: str) -> dict:
    """Expose the bundled IND quantities (AUPR + CI, per-held-out-guide coverage,
    per-band guide/positive counts, prevalence, provenance) for the report/modal."""
    m = load_artifact(scorer)["meta"]
    return {k: m.get(k) for k in (
        "scorer", "n_guides", "n_positives", "prevalence", "aupr",
        "loo_total_count_coverage_nominal90", "band_counts", "phi_bands", "provenance")}

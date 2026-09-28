# Phase 0 — CRISPR-Bulge vs CRISTA: benchmark + feasibility (verdict: **REPLACE**)

**Date:** 2026-09-26 · **Host:** ml007 (2× A100-PCIE-40GB) · **Scratch:** `/srv/local/lp698/cbulge_bench`
**Scope:** throwaway envs + standalone drivers only — **no CRISPRme pipeline code touched.**
**Drivers (in `PostProcess/`):** `bench_cbulge.py`, `bench_crista_speed.py`, `bench_ipc.py`, `scorer_worker.py`.

CRISPR-Bulge = OrensteinLab, *NAR* 2024 (doi:10.1093/nar/gkae428). TensorFlow GRU-embedding
5-model ensemble; **MIT license**. Reference: Yaish & Orenstein, Fig 4, Refined_TrueOT.

---

## 1. Dependency solve — modular env is CLEAN and decoupled

Isolated env solves with no conflict, carrying **modern** pandas/sklearn alongside the pinned
numpy/tf — proving the scorer's deps never have to co-solve with the main image's stack:

```
cbulge_cpu:  python=3.10 tensorflow-cpu=2.13 numpy=1.23.5 pandas=2.3.3 scikit-learn=1.7.2
             xgboost=3.2 catboost=1.2.10 setuptools<81       (279 MB download)
cbulge_gpu:  same + tensorflow=2.12=cuda* (CUDA build)       (455 MB download, sees A100)
```

**⇒ This is the key structural win:** with CRISTA (and later azimuth) moved out of the main env,
the `scikit-learn=1.1.3 / numpy=1.24.4 / scipy=1.10.1` pin can finally be dropped.

Gotchas (all handled; noted for the Dockerfile + prod worker):
- The encoder import chain (`base_models`) pulls **xgboost + catboost + pkg_resources** even for the
  NN-only path → include them, or trim via lazy-import for a lean worker.
- **setuptools ≥ 81 removes `pkg_resources`** → pin `setuptools<81`.

## 2. Model weights — small, committed, fast

5-model c_2 GUIDE-seq-finetuned classification ensemble = **53 MB**, committed in-repo (not LFS).
Load ~**1.0 s** in-process. Input = aligned `Align.sgRNA`/`Align.off-target` (23-mer, `-` gaps),
`seq_len=24` fixed for bulges with left-pad. **No 29-nt genomic flank needed** (CRISTA needs a 29-nt
DNAshape context; CRISPR-Bulge is simpler to feed from CRISPRme's aligned pairs).

## 3. Accuracy — reproduced exactly; decisively beats CRISTA on bulges

Refined_TrueOT (354,352 aligned pairs, 100 positives, 121,039 bulge rows), AUPR:

| Metric | CRISPR-Bulge (measured) | Paper | CRISTA (paper) | Advantage |
|---|---|---|---|---|
| Full-set AUPR | **0.4976** | 0.498 ✅ | ≤ 0.237 | **2.1×** |
| Bulge-only AUPR | **0.2451** | ≥ 0.213 ✅ | ≤ 0.045 | **5.4×** |

Reproducing the published number confirms the vendored weights + our encoding path are correct.
The bulge regime is exactly CRISPRme+'s differentiator (SNP+indel co-occurrence, bulge off-targets).

## 4. Speed — CPU beats CRISTA with proper threads; GPU is a large bonus

Throughput at 100k batch (OTs/sec) and projected genome-wide time (6.77M OTs):

| Scorer / config | OT/s @100k | GW projection | peak RSS |
|---|---|---|---|
| CRISTA (v2.5.5 SIF, 1 proc) | 9,218 | 12.2 min | 1.7 GB |
| CRISPR-Bulge CPU (TF default threads) | 8,975 | 12.6 min | 2.7 GB |
| **CRISPR-Bulge CPU (intra-op=16)** | **14,451** | **7.8 min** | 1.3 GB |
| **CRISPR-Bulge A100 GPU** | **44,604** (51,184 @1M) | **2.2 min** | 4.5 GB |

CPU thread sweep @100k: t1=4,675 · t4=8,984 · t16=14,451. The model is small enough that ~4–16
threads saturate, so per-worker CPU footprint is modest → parallel workers per chromosome remain
viable (same pattern as CRISTA's current parallelization).

## 5. IPC / modular-env boundary — negligible

Persistent worker in its own env (invoked by **absolute env-python path**, fed newline-JSON batches):

| batch | worker round-trip | in-process | serialize |
|---|---|---|---|
| 1k | 1,873 OT/s | 1,878 | 0.6 ms |
| 10k | 6,563 OT/s | 6,380 | 5.3 ms |
| 100k | 9,200 OT/s | 8,975 | 53.7 ms (~0.5%) |

Round-trip ≈ in-process; one-time spawn+load ~4.5 s per run/chromosome (amortized over ~68 batches
at GW scale). **The dedicated-env-over-subprocess design costs essentially nothing at batch scale.**

Worker plumbing (fixed in `scorer_worker.py`; required for any prod runner):
- invoke the env's python by **absolute path** — `mamba run` captures/buffers stdio and breaks streaming;
- **dup fd1 to a private protocol channel + redirect fd1→stderr** — the encoder/Keras print to stdout;
- use `sys.stdin.readline()` loop, not `for line in sys.stdin` (read-ahead deadlock).

---

## Verdict: **REPLACE**

All REPLACE gates cleared:

| Gate | Result |
|---|---|
| Accuracy ≥ CRISTA on bulges | ✅ **5.4×** (0.245 vs ≤0.045); 2.1× overall |
| CPU throughput ≥ CRISTA | ✅ **1.57×** faster @16 threads (7.8 vs 12.2 min GW); ≈parity at TF default |
| IPC overhead small | ✅ ~0.5% at 100k + one-time 4.5 s spawn |
| Per-worker RSS acceptable | ✅ 1.3–2.7 GB (≈ CRISTA's 1.7 GB) |
| Clean dependency path | ✅ isolated env frees the main-env sklearn/numpy pin |
| License | ✅ MIT |
| GPU (optional) upside | ✅ 3.5–5.7× (2.2 min GW), one shared compute-backend switch |

**Recommendation:** retire CRISTA → CRISPR-Bulge as the off-target ML score, running in its own
`cbulge` conda env behind the persistent-worker boundary; keep CFD as the fast in-process default;
make GPU optional (also accelerates CFD). Then isolate/retire azimuth and modernize the main env.

**Next (on GO):** P1 standalone prototype scorer (`CRISPR_BULGE_predict_list`) → P2 scorer-runner →
P3 compute-backend → P4 worker + HF-hosted weights → P5 GPU-CFD → P6 retire CRISTA →
P7 isolate azimuth + modernize main env → P8 report/flags/deps. Lands on `dev` for a later version;
independent of the held v2.5.6 release.

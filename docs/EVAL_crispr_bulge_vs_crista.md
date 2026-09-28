# Evaluation: replace CRISTA with OrensteinLab/CRISPR-Bulge?

> **OUTCOME (2026-09) — DONE: CRISTA is REPLACED by CRISPR-Bulge.** The go/no-go spike below
> passed, so the recommendation that follows ("do not replace yet") is **superseded** and kept
> only as the historical rationale. Final implementation + results:
> - CRISPR-Bulge (TF GRU ensemble, MIT) is now the sole ML off-target score; **CFD remains the
>   primary score**. It runs in a dedicated `cbulge` conda env over a persistent-worker batch
>   protocol, so the main environment no longer carries CRISTA's scikit-learn/numpy pin (those
>   pins were retired and the numerical stack modernized).
> - Fidelity to the upstream model reproduced bit-for-bit on CPU; accuracy ~2.1× overall / ~5.4×
>   on bulge off-targets by AUPR. Report/web thresholds re-derived on its own 0–1 scale
>   (0.5/0.2/0.1); the physical `CRISTA_score` column name is kept only as a stable identifier.
> - Compute backend `cpu|cuda|metal` (opt-in GPU, graceful CPU fallback). Validated on an NVIDIA
>   A100 (numerically correct, ~14.7k off-targets/s); tensorflow-metal is guarded by a load-time
>   numerical self-test (it miscomputes the GRU) that falls back to CPU.
> - **Multi-bulge:** the model scores ≤1-bulge alignments; ≥2-bulge off-targets are nulled to −1
>   and carried by CFD + edit distance (an adversarial review rejected the earlier "collapse"
>   reduction as unsafe on a non-monotonic model). Validated on real data: no ≥2-bulge alignment
>   is ever ML-scored.
>
> See `CHANGELOG.md` for the shipped summary. The original evaluation follows.

Grounded eval (workflow `wf_37c7f6d7`, 9 agents: web research + repo grounding +
adversarial verify). **Motivation is accuracy, not speed** — CRISTA is ~4% of runtime.

## Bottom line (ORIGINAL — superseded, see OUTCOME above)
**Do NOT replace CRISTA now. Do NOT commit it to a release yet.** The right next step is a
**dependency go/no-go spike**. If that passes, **ADD** CRISPR-Bulge as a *bulge-specialized*
score in a **future** release (not 2.5.3, not a patch); defer any **replace** to a later
release gated on a repo-local benchmark.

## The strong case FOR it
- **Fills a real gap.** CFD *cannot* score bulges — the gap `-` triggers a KeyError swallowed
  by a bare `pass`, so a bulge position silently contributes factor 1.0
  (`new_simple_analysis.py:236`, `analisi_indels_NNN.py:135`); METHODS §8 admits the CFD-on-bulge
  value is "extrapolation beyond the training domain." CRISTA isn't bulge-specialized either.
- **Bulge-first-class + published accuracy.** Gap-aware one-hot + GRU-Emb; NAR 2024
  (PMC11229338) bulge-only AUPR **~0.21–0.29 vs CRISTA ≤0.045** (~5–6× on CRISPRme's exact weak
  spot); full-dataset ~2× CRISTA. (Paper does NOT benchmark CFD.)
- **Smaller + permissive.** ~14 MB Keras ensemble vs CRISTA's **276 MB** pkl; **MIT** licensed.
- **Clean seam.** `CRISTA_predict_list(sg, off, 29nt)` (`CRISTA_score.py:684`), 3 importers —
  the eventual swap is mechanically easy.

## Decisive blockers / caveats
1. **DEPENDENCY COLLISION (primary go/no-go, verified).** CRISPR-Bulge needs **TF 2.12 →
   numpy<1.24**, but CRISPRme+ **hard-pins numpy=1.24.4** because the azimuth + CRISTA pickles
   are version-locked to it (`environment.yml:17-19`; `CRISTA_score.py:473` disclaims
   bit-identical output across versions). In an ADD scenario you keep CRISTA, so you **cannot**
   relax the pin. Realistic path = **process/env isolation** of TF inference (separate conda env
   / subprocess so TF's numpy never shares the pickle stack) OR a newer TF that loads the `.h5`
   on numpy≥1.24 — unproven across linux-64 / aarch64 / macOS-arm64 + the crispritz build. If
   none works cleanly, **don't force it**.
2. **No repo-local accuracy validation exists** (METHODS §8): the AUPR is on external assays
   (CHANGE-seq/GUIDE-seq, SpCas9), untested on CRISPRme's NRG/NGG default PAM + our own targets.
   A retrospective head-to-head is the gate before any REPLACE.
3. **Scope caveat** — actionable (high-CFD/CRISTA) off-targets are structurally *low-edit* with
   few/no bulges; bulge-heavy sites cluster in the ~0-CFD weak tail, so the *practical*
   decision-making lift is likely smaller than the headline AUPR ratio.
4. **Calibration** — must use the **classification head (0-1)**, whose scale is NOT calibrated to
   CRISTA's; report thresholds (`CRISTA_THRESHOLDS=(0.6,0.4,0.2)`) must be **re-derived**
   (gates `crista_computed()` / the tier system — a correctness blocker, not polish).
5. **Scoring ≠ search.** It only re-scores bulge sites the engine already emits; it does NOT fix
   the search-side ≥2-cis-indel candidate-generation gap (METHODS §5).
6. **Unmaintained** research artifact (16 commits, 2 stars, custom `Model.load_model_instance`
   wrapper) → we'd own any future TF/keras load-format break.

## Recommended path
1. **Dependency spike FIRST** (standalone go/no-go): can CPU-only TF inference be isolated
   from the numpy=1.24.4 pickle stack across all 3 platforms + crispritz? If no → stop, keep CRISTA.
2. If yes → **ADD** as a 4th `.best<Model>` projection with **new named columns** (avoid the
   brittle positional integrator indices `resultIntegrator.py:66/68/69/70`), re-derived
   thresholds, weights-only bundle (~14 MB, MIT), behind the graceful-absent report path.
3. **Retrospective benchmark** on our own targets (e.g. the TRAC rhAMP-Seq ~150-site panel +
   a public GUIDE-seq/CIRCLE-seq set): CFD vs CRISTA vs CRISPR-Bulge AUPR, overall + bulge-only.
4. **REPLACE only later**, gated on (3) confirming CRISPR-Bulge ≥ CRISTA on our bulge targets —
   at which point removal also drops 276 MB→14 MB + the sklearn-unpickle-shim debt.

**Not a fit for 2.5.3** (adds a TF runtime + needs the benchmark first). 2.5.3 ships as-is.

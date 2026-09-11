# Design: fast-mode CRISTA-skip / CFD-only fallback

> **⚠️ SUPERSEDED by measurement (do not implement as a perf fix).** A direct CRISTA
> throughput bench (on-SIF, 40k targets) measured **~7,500 targets/s serial** and only
> **1.27× at 8 workers**. Cross-referenced against the version-matrix (179,843 bestCRISTA
> rows in the 100-min-timeout run), **CRISTA is ~1–3% of the post-analysis tail** — the
> ~58-min SNP step is row production (enumeration + Tier-0/1 lookups + CFD + emit) and the
> indel post-analysis is single-threaded. A CRISTA-skip would save ~1–2 min of a 100-min
> timeout. The real lever is **#174** (parallelize per-contig row production + the indel
> path). This doc is retained only as the (correct-but-moot) design that *would* apply if
> CRISTA were ever the bottleneck; the `--cfd-only` flag survives at most as a minor
> convenience, not a performance mechanism. See `RELEASE_REPORT_2.5.2.md` §4.


Grounded design (workflow `wf_f83df9e3`, 13 agents) for letting fast mode drop the
expensive CRISTA best-effort screen so a search always **finishes** (CFD stays the
exact worst case). Motivated by the version-matrix: a very-dense guide
(`TGCTTGGTCGGCACTGATAG`, mm5/b2/b2, chr22 genotyped) timed out at 100 min in **both**
full and `--fast` — the tail is candidate-volume CRISTA scoring + the single-threaded
indel post-analysis.

## Key grounding facts (code-verified)
- **A built-in CRISTA-bypass primitive already exists:** `preprocess_CRISTA_score`'s
  `if do_scores: … else:` branch (`new_simple_analysis.py:1518-1527`) writes the
  `-1.000` sentinel into the CRISTA column (index `-2` + appended dup) and returns
  **without loading the 276 MB model**. The indel twin has the mirror at
  `analisi_indels_NNN.py:582-591`. But `do_scores` also gates CFD, so a CRISTA-only
  skip needs a **separate** flag.
- **CFD is independent of CRISTA** (`preprocess_CFD_score`) and stays byte-exact.
- **The report already degrades gracefully to CFD-only:** `crista_computed()`
  (`generate_report.py:371`) masks `-1.000`→NaN and returns False when *all* values are
  out of range → 2 scatter panels instead of 4, CRISTA column dropped from tables, the
  Section-4 panel falls back to CFD ranking. Tested (`test_generate_report.py:713`).
- **Invocation is per-100k-batch inside a per-contig process** (`new_simple_analysis.py:2228`),
  CRISTA runs alt+ref = ~2× the batch rows.

## Why a candidate-count cap does NOT guarantee termination (the correction)
A per-contig count cap (my first instinct) fails as a guarantee: it doesn't bound the
CFD/IUPAC-decomposition cost (that's `CRISPRME_IUPAC_CAP`'s job), it trips only at a
100k batch boundary *after* crossing (the first over-cap batch still runs full CRISTA),
and it under-counts (alt+ref ≈ 2× the counted rows). It's a heuristic, not a guarantee.

## Recommended design (hybrid)
1. **Deterministic explicit lever (the primary, reproducible answer):** `--cfd-only`
   (alias `--no-crista`) → `CRISPRME_SKIP_CRISTA=1` → routes ALL CRISTA (SNP + indel)
   through the existing null-emit branch, emitting `-1.000`; CFD byte-exact; never
   touches `--full`. Fully deterministic.
2. **Optional wall-clock backstop ("always finishes"):** `CRISPRME_CRISTA_BUDGET_SEC`
   (0 = off). When it trips mid-contig, an **all-or-nothing** rewrite forces that
   contig's CRISTA rows all to `-1.000` — so `crista_computed()` sees an *honest*
   all-absent state (CFD-only report), never the non-deterministic **partial** column
   that a naive time-budget would leave (the killer flaw the judges rejected in the
   partial-budget design). Non-deterministic *whether* it trips (host speed) → CI /
   paper runs set `CRISPRME_CRISTA_BUDGET_SEC=0`.
3. **Rejected:** "skip CRISTA when parallelism is ineffective" — inversely correlated
   with risk (skips on the few-core boxes where CRISTA is cheapest; attempts it on the
   many-core boxes where dense GW actually hangs) + non-reproducible.

## Never silent
`--cfd-only` flag; launch banner line; `.search_mode` gains a `cfd-only` token +
`<out>.crista_skipped.txt` audit companion; report "Analysis inputs" note: *"CRISTA
scoring was skipped (CFD-only) to guarantee completion; CFD is the exact worst case."*
Web fast/full toggle gets a matching `--cfd-only` checkbox.

## Must-dos
- **Mirror the change in the indel twin** (`analisi_indels_NNN.py`) — it's a co-equal
  bottleneck; a SNP-only fix is incomplete. Keep SNP+indel guard logic in a tiny shared
  helper to avoid drift.
- **Guard `--full --cfd-only`** (contradictory) → error/warn; `--full` forces budget 0.
- **Tests:** build-report-from-all-`-1.000`-fixture e2e; assert `CRISPRME_SKIP_CRISTA=1`
  never calls `CRISTA_predict_list` while CFD is byte-identical; CI sets budget=0.
- **Validate e2e on the motivating guide** before shipping (don't assume finish from the
  SNP measurement — the indel path must be exercised).

## Open questions
- Backstop default posture: OFF (reproducible default, user opts in) vs ON at ~30 min
  (always-finishes default, host-dependent when tripped). **Product call.**
- Budget value calibration (30 min guess) against the motivating guide + a normal dense guide.
- Web UX: `--cfd-only` visible vs behind an "advanced" expander.

## Scope
Small surface (≈2-line guard per file reusing the verified null branch + the backstop
latch/rewrite). Synth recommends fold-into-2.5.3; the counter-argument is 2.5.3 is
already validated + staged and this needs its own e2e validation + budget calibration.
**Scope + backstop-posture are the two decisions to make before implementing.**

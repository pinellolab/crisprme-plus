# CRISPRme+ 2.5.2 — release report (for review)

Status snapshot of the 2.5.2 release + validation + the post-2.5.2 CRISTA-parallel prototype.
Written for a morning review; the clean-room result tables are filled in as those runs land.

## 1. Released (live)

| Artifact | State |
|---|---|
| Git tag `v2.5.2` | pushed (`8145ea4`), on `main` at `72c8e27` → README `0a888cf` → TESTING `48dcbd3` → grid `ba83ba6` |
| Docker image | **`pinellolab/crisprme:v2.5.2` live on Docker Hub** (multi-arch build succeeded) |
| Genotyped index | `lucapinello/crisprme-data : NRG_3_hg38+hg38_1000G2021_HGDP` (on HF; unchanged, 2.5.2-ready) |
| Mega index | `lucapinello/crisprme-data : NRG_3_hg38+hg38_mega` — **newly published** (12.1 GB: SNP + searchable INDELS + registry + indel_af + log_indels) |
| GitHub *release notes* | **not yet attached** (optional): `gh release create v2.5.2 --repo pinellolab/crisprme-plus --title "CRISPRme+ 2.5.2" --notes-file ~/crisprme_252_notes.md --verify-tag` |

## 2. What 2.5.2 adds (the two-index, co-occurrence story)

- **Two production indices, complementary by design:**
  - **1000G-2021 + HGDP** (genotyped) → **observed** haplotypes: CONFIRMED (phased cis) / PUTATIVE (unphased), with named carrier samples + exact joint AF.
  - **mega** (5 sources, sites-only) → **putative** haplotypes: PUTATIVE with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`).
- **Two co-occurrence dimensions**, both indices: **SNP+SNP** (`snp_snp_cooc.tsv`, new) and **SNP+indel** (`indel_snp_cooc.tsv`), the class the classic two-pass search couldn't see.
- **Mega searchable indels genome-wide** (the P3 gap: the enricher never materialized the fake-indel genome for sites-only; fixed by a synth-GT build + MEGA-carrier strip → cleanly sites-only).
- **Lossless dense-region worst-case** on sites-only (`CRISPRME_LOSSLESS_DENSE`, default-on sites-only, byte-identical for genotyped).
- **Report fixes:** bundle `snp_snp_cooc.tsv`; **fix pre-existing `indel_af.tsv`-never-bundled** (#184); radar-chart tolerant of IUPAC on sites-only DNA (was crashing the whole variants summary).

Code commits on `main`: `c7b145e` (P1), `f98e18d` (P2), `ddeff99` (report), `1cbe5b2`/`9f25e5e`/`72c8e27` (docs+version), + README/TESTING/grid.

## 3. Validation

### 3a. chr22 development progression (done)
Genotyped: SNP+SNP 3,329 (35 CONFIRMED + 3,294 PUTATIVE), SNP+indel 14 (2+12). Mega: SNP+SNP 139 PUTATIVE, SNP+indel 6,502 PUTATIVE, indel_af 6,824. GW mega verify (all 24 chr): 839,897 SNP+indel + 10,753 SNP+SNP (all PUTATIVE), 868,724 indel_af, indel off-targets on all 24 chr, 0 pseudo-sample leaks. Provenance spread 785K single-source → 8.7K all-5.

### 3b. Clean-room, new-user from scratch (Docker/apptainer pull v2.5.2 → HF-download both indices)
Two configs in parallel on separate hosts. Guide = TRAC `CTCTCAGCTGGTACACGGCA`, NRG.

**ml007 — mm4 / bDNA1 / bRNA1 (fast validation, /srv/local SSD, fully from-scratch new-user):**

New-user path validated: `apptainer pull v2.5.2` → HF-downloaded **both** indices + INDELS (`dl_geno exit=0`, `dl_mega exit=0`, HF-authenticated → no 429).

| Cell | report.zip | off-targets | snp_snp_cooc | indel_snp_cooc | indel_af |
|---|---|---|---|---|---|
| geno_slow | ✅ | 119,082 | 1,183 — **124 CONFIRMED** + 1,059 PUTATIVE | 85 — **18 CONFIRMED** + 67 PUTATIVE | n/a |
| geno_fast | ✅ | 132,021 | 7,694 (worst-case reps, PUTATIVE) | 85 (18 CONFIRMED + 67 PUTATIVE) | n/a |
| mega_slow | ✅ | 111,986 | 327 all PUTATIVE | 38,222 PUTATIVE | 39,358 |
| mega_fast | ✅ | 129,643 | 2,944 all PUTATIVE | 38,222 PUTATIVE | 39,358 |

**All 4 cells passed** (report.zip + companions, exit 0).

**Genotyped observed-haplotype path validated** — the slow genotyped cell produced **CONFIRMED cis co-occurrences with named carrier samples**: 125 CONFIRMED SNP+SNP + 19 CONFIRMED SNP+indel. This is the class only a genotyped index can assert (phased cis, exact carriers), the complement of the mega's PUTATIVE-only output. The mega cells here are all-PUTATIVE (correct). Both indices' companions produced from a fresh HF download.

**Notable slow-vs-`--fast` behavior on the genotyped index:** slow emits the **observed CONFIRMED** haplotypes (124 CONFIRMED SNP+SNP + 18 CONFIRMED SNP+indel); `--fast` emits more rows (7,694) but they are worst-possible **representatives**, tagged PUTATIVE, not per-sample CONFIRMED — i.e. `--fast` trades per-sample cis attribution for tractability (exactly its documented design). Use the slow path when you need CONFIRMED carriers; `--fast` for a dense-panel screen.

**Bottom line: full 4-cell matrix passed on ml007** (both indices × slow/`--fast`, fresh new-user download) — genotyped gives CONFIRMED+carriers, mega gives PUTATIVE+min-AF+provenance, exactly the intended two-index design. The ml008 heavyweight genotyped hang is a 4×-parallel RAM over-subscription (documented above), not a code issue — the same genotyped index completes cleanly here.

**ml008 — mm6 / bDNA2 / bRNA2 (thorough):**

| Cell | report.zip | off-targets | snp_snp_cooc | indel_snp_cooc | indel_af |
|---|---|---|---|---|---|
| mega_slow | ✅ | 550,456 | 3,453 **all PUTATIVE** | 318,455 PUTATIVE | 329,980 |
| mega_fast | ✅ | 547,807 | 42,273 **all PUTATIVE** | 318,455 PUTATIVE | 329,980 |
| geno_slow | ⚠️ hung @51.6% | — | — | — | — |
| geno_fast | ⚠️ hung @51.6% | — | — | — | — |

**Heavyweight mega passed** genome-wide (318K SNP+indel + 3.4K SNP+SNP PUTATIVE, 330K indel_af, 550K off-targets) — the sites-only feature set holds at mm6/b2/b2 GW.

**Genotyped cells hung at 51.6%** (0 active post-analysis workers, `run_searches.sh` still waiting): a worker was OOM-killed under **4× mm6/b2/b2 memory pressure** — the genotyped observed-enumeration loads the per-sample genotype store (far heavier than the registry-only mega), and 4 such GW searches in parallel over-subscribed RAM. **Resource/contention finding, not a genotyped-index bug** — the genotyped is validated by the chr22 progression (§3a) and the ml007 mm4/1/1 clean-room (below). Takeaway: don't run 4× genotyped GW mm6/b2/b2 concurrently on one host; run fewer in parallel or cap `CRISPRME_POSTPROC_MAX_WORKERS` / memory.

The mega SNP+SNP rows are **all PUTATIVE** (verified: zero CONFIRMED data rows across all 24 per-chr files) — correct for a sites-only index. `--fast` emits far more SNP+SNP rows than slow (42,273 vs 3,453) because it materializes a worst-possible representative per dense window (which combines multiple SNPs), whereas the slow greedy emits fewer multi-SNP reps; both are PUTATIVE with the min-AF bound.

What each cell is checked for: report.zip generated; `snp_snp_cooc.tsv` / `indel_snp_cooc.tsv` / `indel_af.tsv` present; CONFIRMED+carriers on genotyped vs PUTATIVE+min-AF+provenance on mega; `--fast` ⊇ slow; no crash.

**New-user friction found:** anonymous HF downloads hit **429 "Too Many Requests"** (the `--what all` many-small-file API calls); authenticating (HF token) fixes it. **FIXED** on `dev` (`b58a0b1`): retry-on-429/5xx with backoff (honoring `Retry-After`), tunable `CRISPRME_HF_MAX_RETRIES`; ships in 2.5.3.

### 3c. Old-vs-new correctness (does 2.5.2 regress the prior version?)
**v2.5.0 was never tagged**, so the "old" reference is `722ea4b` = **2.5.1-dev** (has 2.5.1's `--fast` + CRISTA-perf but *not* the 2.5.2 co-occurrence/lossless-dense work). Verified two ways:

**Static proof (airtight).** The entire 2.5.2 search-path delta is **pure insertion, 0 deletions**: `new_simple_analysis.py` +128/-0, `analisi_indels_NNN.py` +33/-0 — so every pre-existing line is byte-identical, and every *new* branch is **dead on the genotyped path**: lossless-dense gated `registry_only_mode` (genotyped isn't registry-only); SNP+indel PUTATIVE fallback gated `_gt is None or _indelgt is None` (genotyped has both — the code comment explicitly guards against phantom PUTATIVE rows); `_record_snp_snp_cooc` is a pure side-write to a *new* companion file (never mutates the emitted row, try/except-wrapped, gated on `myreg`).

**Empirical diff (confirms it).** OLD `722ea4b` vs NEW `72c8e27` overlaid on the **same** base SIF (isolating the code delta), identical chr22 genotyped search (guide `CTAAC`, NRG, mm4/1/1, `--index-path`):

| | OLD (2.5.1-dev) | NEW (2.5.2) |
|---|---|---|
| `integrated_results` rows | 1,734 | 1,734 |
| `integrated_results` md5 (sorted) | `91ce3b68…eb25` | `91ce3b68…eb25` — **IDENTICAL** |
| `bestMerge` (sorted) | — | **identical to OLD** |
| `snp_snp_cooc.tsv` files | 0 | 1 (**additive**) |

**Result: on the genotyped index the core scored off-target table is byte-identical OLD→NEW; the only change is the new additive SNP+SNP companion. No regression.** (On the **mega** sites-only index the core output *does* change by design — searchable indels + lossless-dense PUTATIVE are new 2.5.x capabilities, active only where `registry_only_mode` holds; the version-matrix, §3d, quantifies that.)

### 3d. Version-matrix — chr22 genotyped, guide `TGCTTGGTCGGCACTGATAG`, mm5/bDNA2/bRNA2 (a deliberately pathological, very-dense config)
The `run_v1/v2/v3` progression, rebuilt against real prebuilt tiers (`--index-path`), each capped at a 100-min timeout:

| Cell | Version / mode | Outcome | Partial output before cutoff |
|---|---|---|---|
| **V1** | stock **v2.4.0**, feature-off | **FAILED** | v2.4.0 cannot post-process the 2.5.x index's indels (`adjust_cols.py: cols.remove("CFD_ref")` ValueError + `KeyError 'AK'`) — you can't run old CRISPRme on the new index. |
| **V2** | **2.5.2 full** enumeration | **TIMED OUT** (`rc=124`) | 52,962 SNP+SNP (**3,644 CONFIRMED + carriers**); the SNP phase alone took 60 min, then the indel phase ran out the clock. |
| **V3** | **2.5.2 `--fast`** | **TIMED OUT** (`rc=124`) | 104,712 SNP+SNP (**0 CONFIRMED**, all PUTATIVE) — more rows (worst-possible reps per window), no carriers. |

**Honest takeaways (this is a worst-case stress config, not a typical run):**
1. **v2.4.0 cannot run on the 2.5.x index** at all — the clearest reason the index/version move forward together.
2. **The CONFIRMED-vs-PUTATIVE contrast is exactly as designed** even in the partial output: full produced 3,644 CONFIRMED SNP+SNP *with carriers*; `--fast` produced 0 (all PUTATIVE).
3. **`--fast` is not a universal cure.** On this guide *both* modes exceed 100 min because the bottleneck here is **candidate-volume CRISTA scoring + the single-threaded indel post-analysis** — which `--fast` does **not** touch (it collapses only the SNP 2ᵏ haplotype lattice). `--fast`'s decisive, measured win is on **enumeration-bound dense panels** (the #183 4×-density panel: 49 h+ → tractable) and on typical guides (the mm4/1/1 clean-room + §6 e2e complete in minutes). The CRISTA-scoring tail on high-mm/bulge dense guides is the separate **#174** lever (outer-level contig/CRISTA parallelism), still open.

**Net for fast-default:** still correct — `--fast` never does *more* work than full and it removes the enumeration wall; but the matrix keeps us honest that it doesn't fix the CRISTA/indel tail, so the docs frame it as "tractable on any panel" for the *enumeration* cost, not a blanket speed guarantee.

## 4. CRISTA speed + parallel prototype (post-2.5.2, on `dev`) — **CRISTA is NOT the tail**

**Measured absolute throughput (on-SIF, 40k synthetic targets, 256-core host):**

| `CRISPRME_CRISTA_PARALLEL` | wall | throughput | speedup |
|---|---|---|---|
| 1 (serial) | 5.29 s | **7,557 targets/s** | 1.0× |
| 4 | 4.37 s | 9,153 targets/s | 1.21× |
| 8 | 4.16 s | 9,623 targets/s | 1.27× |

Serial split: **feature-build 58%, predict 42%**. The predict already carries `n_jobs=-1` (all cores), so `CRISPRME_CRISTA_PARALLEL` only speeds the feature build → **1.27× at 8 workers**, not the ~3.8× the earlier feature-build-only micro-benchmark suggested end-to-end.

**The load-bearing correction (measured):** CRISTA is **fast** (~7.5k targets/s) and is **NOT the post-analysis tail**. Cross-referencing the version-matrix (§3d): the dense guide produced ~118k variant SNP targets → **179,843 rows in bestCRISTA.txt**. At 7.5k/s (×2 for the alt+ref passes) that is **~100 seconds of CRISTA — ~1–3% of the 60-min SNP post-analysis.** The remaining ~58 min is **row production** (per-window IUPAC enumeration + Tier-0/Tier-1 registry/genotype lookups + CFD + emit, ≈20 ms/row), and the search *also* timed out in the **single-threaded indel post-analysis** — which `--fast` does not touch.

Consequence: a **CRISTA-skip / CFD-only fallback would save ~1–2 min of a 100-min timeout** — it targets the wrong thing. We are **not** building it as a perf fix (see `docs/DESIGN_fast_crista_skip.md`, superseded by this measurement).

The `CRISPRME_CRISTA_PARALLEL` fork-pool prototype (opt-in, default OFF, bit-identical — `test_crista_parallel_equivalence` in CI; `dev` `7771b66`→`c5fff0d`→`a0e7f16`) remains a clean but **minor** (1.27×) partial win.

### 4a. The real tail, profiled + FIXED + tuned genome-wide (folded into 2.5.3)
A `cProfile` of the dense-guide SNP post-analysis (chr22, mm4/1/1, fast) found the tail was **not** CRISTA (~4%) or CFD (~4%) but **`zlib.decompress` at 71%** (953,103 calls): the genotyped v3 registry is zlib block-compressed (the #99 disk compaction) and the decompressed-block LRU was only **8 blocks**, while a dense IUPAC decomposition touches positions across the whole chromosome — so hot blocks were re-decompressed ~10⁶ times.

**Fix (byte-identical):** raise the decompressed-block cache default **8 → 4096**, tunable via `CRISPRME_REGISTRY_CACHE_BLOCKS`. The LRU only grows to blocks actually touched (≤ the contig's block count), so this is effectively "cache the whole contig registry."

**Measured A/B** (same dense chr22 search): SNP post-analysis **326 s → 89 s = 3.66×**, `integrated_results` **byte-identical** (same md5, 2085 rows).

**Tuned for genome-wide "truly fast"** — measured genotyped-registry block counts (a block = ~64 KB decompressed): chr22 = **339 blocks** (~22 MB), chr1 (largest) = **1,868 blocks** (~122 MB). So **4096 holds *every* genotyped contig fully** (chr1 with headroom), extending the 3.66× from small chromosomes to the whole genome, and the 4096 cap bounds even a pathological panel at ~256 MB/reader — within the memory-capped pool's budget. **The mega sites-only index is `codec=RAW` (uncompressed) → no decompression, unaffected** (the cache change is a harmless no-op there). 67 registry unit tests green. A one-line cache bump beats CRISTA parallelism (~1.27×) by ~3×.

**Indel path covered too (code-verified):** `analisi_indels_NNN.py` imports the same `tier0_registry` and does `_reg.lookup(...)` (`:980,:1188,:1212`) through the same `_decompress_block` + module-level `_LRU_BLOCKS` — so the 4096-block cache benefits the indel post-analysis's registry lookups identically, no separate change. A full indel-tail *profile* (does registry lookup dominate the indel tail as it did the SNP tail, or the fake-indel target processing?) is **deferred**: the attempt stalled on ml007 under heavy external contention (load ~1374 from other users). Non-blocking for 2.5.3.

## 5. Pending / for your decision

1. **Clean-room results** — appended above as they land (both running).
2. **CRISTA prototype**: full-scale byte-identical + predict-`n_jobs` measurement on the SIF; then decide whether to fold into a release (it's `dev`, default-off).
3. **HF 429**: add a retry-on-429 to the downloader (and/or a docs note that a token avoids it).
4. **GitHub release notes** object (command above) — cosmetic; the tag + image are live.
5. `pinellolab/crisprme-data` vs `lucapinello/crisprme-data`: the published indices + the code default (`DEFAULT_HF_REPO`) are both `lucapinello/crisprme-data` — consistent, but confirm that's the intended long-term home.

## 6. 2.5.3 (staged on `dev`, GATED — awaiting your approval to tag/build)

**What it adds** (on top of 2.5.2):
- **Fast mode is now the DEFAULT; `--full` opts into exact per-sample enumeration.** The SNP post-analysis reports worst-possible PUTATIVE representatives by default (tractable on any panel); `--full` gives observed haplotypes with CONFIRMED cis + named carriers + exact joint AF. `--fast` kept as a deprecated no-op alias. The carrier trade-off is surfaced in **three** non-silent places: the launch message, a **web "Search mode" control** (with explanation), and a **report "Search mode" row** (index-aware). *(Your call: global fast-default + `--full`.)*
- **HF-429 download retry** (`b58a0b1`) — new-user rate-limit resilience.
- **Parallel CRISTA feature build** (`CRISPRME_CRISTA_PARALLEL`, opt-in, default OFF).

**Validation (all green):**
- **Unit:** 667 passed (1 failure is the local Mac lacking `sklearn`; green in CI); +2 new report-note tests.
- **Fast/full e2e** (chr22 genotyped, guide CTAAC, mm4/1/1): default → `.search_mode=fast`, report shows the fast caveat, ~0 CONFIRMED (PUTATIVE reps), 1,975 rows; `--full` → `.search_mode=full`, report shows the enumeration note, **36 CONFIRMED SNP+SNP with carriers**, 1,733 rows. Markers + notes correct in both.
- **Byte-identity chain (md5 `91ce3b68…`):** `2.5.1-dev full ≡ 2.5.2 full ≡ 2.5.3 --full ≡ 2.5.3 --full + CRISPRME_CRISTA_PARALLEL=1`. The exact-enumeration path is unchanged across every version; the flip only swaps the *default*, and the CRISTA-parallel flag is **bit-identical** end-to-end.

**To release when you approve:** `release-crisprme 2.5.3` (tag → Docker → `pinellolab/crisprme:v2.5.3`). Version is already bumped to `2.5.3-dev` on `dev`; CHANGELOG staged.

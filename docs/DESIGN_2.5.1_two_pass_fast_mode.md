# CRISPRme+ 2.5.1 — Two-Pass Fast-Mode Variant Search (design + implementation plan)

**Status:** design-complete, NOT yet implemented. Target: 2.5.1 (after 2.5.0 ships).
Lynchpin adversarially verified + empirically bounded (2026-09-01). This document is
the buildable spec; nothing here is in the shipped 2.5.0 code.

Companion artifacts already in-repo:
- `PostProcess/test_twopass_lynchpin_counterexamples.py` — the two verified lower-bound
  counterexamples + the fixes, as an executable regression fixture (commit 113b9b7).
- Memory: `two-pass-dense-search-design` (the full verdict + measured numbers).

---

## 1. Problem — the dense-panel wall (measured)

On the 4× 1000G-2021+HGDP panel (guide `TGCTTGGTCGGCACTGATAG`, NRG, mm5/bDNA2/bRNA2/max-edits7):

| phase | time |
|---|---|
| index SEARCH (genome-wide) | ~16 h 47 m — cheap, fine |
| **SNP haplotype post-analysis** (`new_simple_analysis.py`) | **49 h+ and did NOT finish** (stuck on chr1) |
| indel cooc post-analysis (`analisi_indels_NNN.py`) | ~12.5 h for chr1 alone |

Root cause: post-analysis enumerates the **2^k IUPAC haplotype lattice** per dense window
(the observed-haplotype path + the greedy `CRISPRME_IUPAC_CAP` fallback). Search is cheap;
**post-analysis is intractable at 4× density.** The interim stopgap is a MAF≥0.1% index; the
real fix is this two-pass, score-agnostic design.

## 2. Design overview

Two passes, both **score-agnostic** and **enumeration-free**:

- **Pass 1 — score-free region find.** Compute the *collapsed-IUPAC min edit* `D` for every
  candidate window (each variant position = its allele SET, each fake-indel contig materialized
  once). `D = min over alignments of [#positions where guide base ∉ allele-set] + [#bulges]`.
  No haplotype enumeration. `D` is a **lower bound** on the true edit distance of any real
  haplotype (proof §5), so thresholding `D ≤ k` is **lossless for DETECTION** (for the cases the
  bound covers — see §6). Output: region + `D` + per-position ambiguity map.

- **Pass 2 — per-score worst-POSSIBLE annotation** (only on Pass-1 hits, and only what's asked):
  - **min-edit**: `D` directly. Parity tiebreak = greedy over PAM-region alleles (prefer
    PAM-creating) + bulge placement (prefer distal-from-seed) — score-free biology.
  - **CFD**: EXACT worst-case via factorization (§7), `O(guide)`, no enumeration.
  - **CRISTA**: non-factorizable (RandomForest) → evaluate on a low-edit candidate SET; flag
    CRISTA-approximate where the set is large. NOT exhaustive — state this in the report.
  - Scope is **worst-POSSIBLE, not worst-observed**: unphased (HGDP) / pop-only (gnomAD, TOPMed,
    AoU) have no/partial haplotypes; the phased store (SNP tiers + indel_genotypes) is only an
    ANNOTATION ("is this worst-case actually carried?"), never a dependency.

Counterintuitive but load-bearing: a fully-ambiguous `NNNN…N` MHC window is the EASY case —
a perfect protospacer is formable ⇒ min-edit 0 ⇒ max for every score ⇒ just report it.

## 3. Pass 1 — score-free region find (implementation)

- Reuse the existing **variant-enriched reference** (SNP → IUPAC) and the **SNP-overlaid
  fake-indel genome** (already built for the 2.5.0 indel+SNP feature — no new genome build).
- Compute `D` with a single banded DP over `guide × window` with **set-membership** at each cell
  (`match if guide_i ∈ S_i`), allowing ≤ bulge-cap gaps. `O(guide · window)`, NOT `O(2^k)`.
- `D = min` over `{ enriched-reference } ∪ { each single-indel fake-contig, SNP-overlaid }`.
- Candidate generation MUST be additive/union (never drop the reference-allele k-mer — see the
  `variant-enrichment-masks-reference-offtarget` risk) and the edit cap must not prune below the
  biological `k`.

## 4. Pass 2 — exact/worst-case scoring (see §7 for CFD)

## 5. Lower-bound correctness — VERIFIED (adversarial panel, 15 agents, 2026-09-01)

**Holds (brute-forced, 0 failures):**
- **SNPs**: per-position set-containment `S_i ⊇ {H_i}` ⇒ `[guide_i ∉ S_i] ≤ [guide_i ≠ H_i]`;
  summed over a shared alignment frame ⇒ `D ≤ d(guide,H)`. (~945K triples, 0 violations.)
- **Bulges**: reuse H's own optimal alignment under set-evaluation ⇒ `collapsed-cost(A*) ≤
  concrete-cost(A*) = d`. (~356K different-placement cases, 0 violations.) Bulges do NOT break it.

**Load-bearing assumptions** (each is a build-time invariant to keep):
1. Membership/union-completeness for SNPs (enricher unions ref+all single-base alts).
2. Bulge invariance (unit cost, content-independent).
3. Shared coordinate frame (SNPs don't change length). ← INDELS violate this.
4. Indel coverage: every carriable indel haplotype must map to ≥1 enumerated fake contig.
5. No lossy IUPAC cap that drops a real allele from `S_i` before seeding the index.
6. Additive (never subtractive) enrichment + adequate edit cap (candidate generation).

## 6. Where it BREAKS (verifier-confirmed) + the fixes

**Break 1 — ≥2 cis indels in one protospacer (the real gap).** Each fake contig materializes ONE
indel, so a haplotype carrying two lives on NO single contig; `D` over-counts by up to the total
bulge cost of the un-modeled indels → a false negative, even at k=0 (perfect match).
Verified: `+1ins/−1del` length-cancelling pair → `D=1` vs real `d=0`; two 2 bp cis deletions →
`D=2` vs `d=0`. (Both pinned in the counterexample fixture.)

**Break 2 — equal-length MNV / block substitution.** The enricher DROPS `len(REF)==len(ALT)>1`
and `_is_indel` excludes it → alt bases in no `S_i`, no contig → false negative (verified `D=2`
vs `d=0`). NOT triggered on the as-shipped (already-normalized) 1000G VCF; risk is UN-normalized
merges. **Fix: `bcftools norm -a` atomization pre-pass** in the enrichment path.

**Fixes for Break 1** (fixture-validated; count-only `n−1` slack is INSUFFICIENT for multi-bp
indels — must be size-aware):
- **(a) flag-all-≥2-cis-indel-windows** — report any such window regardless of `D`. PROVABLY
  lossless; cheapest; over-reports only in dense-indel windows (≈ all repeat loci, already flagged).
- **(b) size-aware slack** — `D' = max(0, D − (Σ bulges − min bulge))`; restores the lower bound.
- **(c) bounded combined-contig enumeration** — materialize the FEW multi-indel contigs in the
  (post-repeat-mask, few) dense windows.

**Empirical justification the gap is negligible (why any of (a)/(b)/(c) is enough):**
- The gap is a **pre-existing** limitation — 2.5.0's single-indel search already can't find these.
- Multi-indel-cis co-occurrence is **rare outside repeats**: raw ~15% → phased-confirmed ~19% →
  non-overlapping ~13% → **RepeatMasker soft-mask excluded ~3–5%** → dedup ~1–2%. Unphased
  (loose) 33–47% is an overestimate (trans carriers can't be cis) — use PHASED only.
- **Genuinely-MISSED off-targets** (H matches guide ≤ budget AND no single-indel rep does):
  chr9/21/22 = 32/6/13 = **~0.1–0.2% of indel off-targets**, and **100% sit at the edit-budget
  ceiling** (all `dH=7, dsingle=8`) — structurally: a miss needs `d(H) > budget − δ`, so for
  1–2 bp indels only the weakest tier (≈0 CFD). Conservative (no PAM, Levenshtein proxy) → true
  number lower. ⇒ single-indel + SNP captures every biologically-relevant off-target.
- **Always RepeatMasker soft-mask repeats** (off-target analysis should regardless) — that alone
  removes ~75% of the apparent multi-indel signal (STR/VNTR length variants stored as many records).

## 7. CFD — exact via factorization (bit-verified) + C++ port

`max_CFD = (∏_i max_{a∈S_i} f(guide_i, a, i)) × max_{p∈PamSet} pam_score(p)` — **bit-exact (0 ULP)**
vs the shipped `calc_cfd` across 30K + ~945K cases. Valid because `calc_cfd` is a pure left-to-right
product of per-position factors (each ∈ [0,1], position-independent) × one PAM factor.

Caveats the port MUST honor:
- Factorize **per bulge-alignment**, then max over alignments — never factorize one collapsed string.
- Only in the `do_scores` regime (guide_len==20, len_pam==3, 3′-PAM); else emit the `-1` sentinel.
- **N-handling DIVERGES between the two `calc_cfd` copies** (`analisi_indels_NNN.py` treats N-in-DNA
  and N-in-PAM as free match=1.0; `new_simple_analysis.py` gives 0.0) — replicate PER FILE/PATH.
- Position-1 gap is free (silent 1.0); T→U on both sides before compare+key; single revcom on the
  DNA base only; key `= 'r'+guide+':d'+revcom(DNA)+','+(i+1)`, 1-based; bulges ARE scored (penalty keys).

**C++ equivalence gate (must pass BEFORE the C++ path is authoritative):**
1. exhaustive per-position table (pos 1..20 × guide{ACGT} × DNA{ACGT,-}) == Python factor;
2. exhaustive PAM (16 dinuc + N-containing + non-canonical) in BOTH modes;
3. ≥1e6 random (guide, concrete DNA, PAM) triples → **raw-double bit-equality** (memcmp, 0 ULP);
4. bulge differential (RNA/DNA gaps, pos-1 free, re-index);
5. factorization oracle: random IUPAC allele-SETS == Python brute-force `itertools.product` max;
6. `do_scores==False` → both emit `-1.000`;
7. **`CRISPRME_CFD_SHADOW`** real-run flag: compute both ways per row, assert bit-equality on the
   chr22 1000G+HGDP indel-snp fixture AND a dense-panel region — ZERO divergences; keep as a CI
   regression against future pickle changes.
Requirements: load the shipped pickles as binary (don't retype float literals); accumulate strictly
pos1→pos20 then ×PAM in `double`, NO `-ffast-math`/FMA/reassociation; assert all mm/pam values ∈ [0,1]
at load (required for the factorization max to be valid). **CRISTA is OUT OF SCOPE for the C++ path.**

## 8. Implementation plan (phased)

1. **Enrichment pre-pass**: add `bcftools norm -a` (atomize MNVs) to the enrichment path (fixes Break 2).
2. **RepeatMasker soft-mask** consumption: flag/annotate soft-masked (lowercase) loci; downgrade
   confidence there (dense-indel windows collapse here).
3. **Pass-1 module** (`twopass_find.py`): banded set-membership min-edit DP over enriched ref +
   single-indel contigs; emit region + `D` + ambiguity map. Multi-indel handling = fix (a) initially
   (flag-all-≥2-cis-indel-windows), (c) as an optimization later.
4. **Pass-2 min-edit + CFD** (`twopass_score.py`): factorized exact worst-case CFD (Python first,
   validated against `calc_cfd` via the oracle test), CRISTA candidate-set eval with the approximate flag.
5. **C++ CFD** (`cfd_exact.cpp`) gated by §7's equivalence suite + `CRISPRME_CFD_SHADOW`.
6. **Fast-mode CLI**: a `--fast` mode that runs Pass-1 only (region + min-edit, no ML score) as a
   screen; full mode adds Pass-2. Report states CRISTA is shortlist-only / not exhaustive.
7. **Regression**: extend `test_twopass_lynchpin_counterexamples.py`; add the CFD oracle + shadow tests.

## 9. Risks / residuals

- Break-1 fix (a) over-reports in dense-indel windows — bounded by repeat-masking; measure the
  over-report rate on a real run before choosing (a) vs (c).
- Candidate-generation masking (enriched index dropping the ref k-mer) is an EXTERNAL no-miss
  dependency, not a property of `D` — audit the enricher/index for union-only behavior.
- CFD C++ IEEE-754 order + per-file N divergence are the equivalence-breakers — the shadow gate catches them.
- CRISTA remains approximate in the fast path — acceptable if flagged; exact CRISTA only via local enum in flagged windows if ever contractually required.

## 10. Interim (pre-2.5.1) stopgap
MAF≥0.1% index for routine runs (collapses the dense-region blowup); full 95.6M-SNP panel overnight
for confirmatory/IND. This is the bridge until the two-pass lands.

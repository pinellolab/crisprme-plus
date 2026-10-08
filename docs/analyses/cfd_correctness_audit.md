# CFD correctness audit

Internal note / paper-supporting. Audits CRISPRme+'s CFD (Doench 2016) implementation for bugs, and
explains why heavily-edited off-targets (e.g. **4 mismatches + 1 bulge**) can still score **CFD > 0.5**.

**Trigger.** A collaborator reviewing the SBDSP1/SDS **6,1,1** search (guide `ACTGAAATCTGTAAGCAGGC`,
SpCas9/NRG, hg38, 5,540,022 off-targets) was surprised that so many top-ranked sites carry CFD > 0.5,
and asked whether CRISPRme+'s CFD has a bug. Short answer: **no bug in the CFD computation** — the
numbers are faithful to the Doench matrix. The high counts come from two intrinsic properties of CFD
(one a known leniency, one an inherited bulge-scoring artifact), both of which argue for ranking
bulge/multi-mismatch off-targets by the empirically-trained **CRISPR-Bulge** score, not CFD.

---

## Method

CFD was checked three independent ways:

1. **Reference equivalence.** A from-scratch reimplementation of the Doench-2016 / CRISPOR CFD
   algorithm (loaded from the same shipped `mismatch_score.pkl` / `PAM_scores.pkl`) was run against
   `new_simple_analysis.calc_cfd` over **25,240 mismatch-only cases** (an exhaustive single-mismatch
   sweep across all 20 positions × all substitutions, 20,000 random 1–6-mismatch targets, and all PAM
   keys). Result: **zero differences, max |Δ| = 0.0** — bit-identical.
2. **Faithful-value recompute.** The reported CFD of real high-CFD bulge sites from the 6,1,1 run was
   recomputed from the stored aligned spacer/protospacer strings, replicating the call-site
   preprocessing (`.upper()` + `T→U` on both strands, PAM = last 2 bases). Every site reproduced the
   reported value to the third decimal (see table below) — the printed CFD is exactly what `calc_cfd`
   computes.
3. **Twin consistency.** CRISPRme+ carries two copies of the scorer — the SNP/variant path
   (`new_simple_analysis.calc_cfd`) and the INDEL path (`analisi_indels_NNN.calc_cfd`). These were
   diffed and pinned together (see "Twin reconciliation").

---

## Finding 1 — the CFD core is correct

Mismatch-only CFD is **bit-identical to the canonical reference** (check 1). The published
product-of-per-position-factors algorithm is implemented correctly; there is no off-by-one, no wrong
key format, no rounding divergence.

## Finding 2 — the reported CFD values are faithful, including the surprising ones

The "4 mismatches + 1 bulge, CFD > 0.5" sites are **not** a recompute error. Decomposition of the
representative CFD-0.918 site (`proto g-TGAAATCTacAAGaAGGCAGG`):

| position | event | Doench matrix factor |
|---|---|---|
| 1  | mismatch `rA:dC` | **1.000** (fully tolerated) |
| 11 | mismatch `rG:dT` | **1.000** (fully tolerated) |
| 16 | mismatch `rC:dT` | **1.000** (fully tolerated) |
| 12 | mismatch `rU:dG` | 0.947 |
| 2  | RNA bulge `rC:d-` | 0.969 |
| —  | PAM `GG` | 1.000 |

→ 1.000 × 1.000 × 1.000 × 0.947 × 0.969 = **0.918** (reported: 0.918). Five such 4mm+1bulge sites were
decomposed; all matched the report exactly (0.918, 0.703, 0.688, 0.678, 0.711).

## Finding 3 — *why* heavily-edited sites score high: dominated by tolerated-mismatch cells

CFD is a **product** of per-position factors, so any factor of 1.0 contributes nothing. The dominant
driver is (A); (B) is a minor, second-order contributor.

**(A) Tolerated-mismatch cells — a pure property of the Doench matrix (the dominant driver).** Of the
240 mismatch cells, **14 (6%) are exactly 1.0** (fully tolerated — Doench measured zero activity loss
for that substitution at that position) and **33 (14%) are ≥0.9**. If a site's mismatches land on these
cells, "4 mismatches" barely dents the score. In the limit, four mismatches all on 1.0 cells leave
CFD ≈ 1.0 *before* the bulge even applies. This is **correct CFD** — but it means **mismatch count is a
poor proxy for CFD**; CFD weights edits by empirical tolerance, and many edits carry ~zero weight. (This
is CFD being lenient — matching Doench's and the collaborator's own 2020 caveat that CFD "may
overestimate.") This mechanism is present in all 335 of the CFD > 0.5 sites below.

**(B) PAM-distal terminal-bulge free pass — a minor, inherited effect.** The Doench bulge-penalty
matrix has gap keys (`r-:dX`, `rX:d-`) for **positions 2–20 only** — none at the PAM-distal terminus
(aligned column 1). A bulge that lands at column 1 has no key, so the lookup falls through
`except: pass` to **factor 1.0 — no penalty**. This is the *correct, intended* behavior: it is exactly
the "solution 3" John Doench endorsed in 2020 ("just use the CFD values for bulges as-is, Fig. 5c,d")
— for a DNA bulge the code trims the elongated DNA back to 20 nt by removing the PAM-distal
nucleotide, so a PAM-distal bulge is *removed*, not mis-scored. Crucially, **interior and PAM-proximal
bulges ARE penalized** (gap keys exist for columns 2–20; e.g. `r-:dA,20` = 0.6, `rC:d-,2` = 0.969) — so
this is not a "free pass for all bulges," only for the most PAM-distal (and biologically most tolerant)
position. In the SBDSP1 run it affects only **51 of the 329** bulge-bearing CFD > 0.5 sites; the other
**278 carry a matrix-penalized bulge** and are high because of (A). The DNA-bulge scoring path was
adversarially verified to implement solution 3 correctly — bit-identical to an independent
implementation over 50,000 synthetic DNA-bulge alignments (0 divergences) and 3,000/3,000 real sites.
It is **inherited from the Doench/CRISPOR convention**, not a CRISPRme+ regression, and not a bug.

### The counts, in one table (6,1,1 SBDSP1 run, 5,540,022 off-targets)

| | count |
|---|--:|
| off-targets total | 5,540,022 |
| **CFD > 0.5** | **335** (0.006%) |
| &nbsp;&nbsp;↳ containing ≥1 bulge | 329 (98%) |
| &nbsp;&nbsp;&nbsp;&nbsp;↳ PAM-distal column-1 bulge (genuinely free, mechanism B) | 51 |
| &nbsp;&nbsp;&nbsp;&nbsp;↳ interior / PAM-proximal bulge (matrix-penalized) | 278 |
| &nbsp;&nbsp;↳ mismatch-only (0 bulge) | 6 |
| mismatch-count among the 335 | 0mm:2 · 2mm:11 · 3mm:20 · **4mm:77 · 5mm:100 · 6mm:125** |
| bulge-count among the 335 | 0:6 · 1:160 · 2:169 |

Read-out: the CFD > 0.5 set is **not** near-perfect sites — **302 of 335 carry ≥4 mismatches** and
**329 of 335 carry a bulge** (and **278 of those 329 carry a matrix-penalized bulge**, so the high
score comes from the tolerated mismatches, not a free bulge). They are heavily-edited sites that CFD
rates high chiefly via mechanism A (tolerated cells, all 335); mechanism B (PAM-distal free bulge)
touches only 51. CRISPR-Bulge, trained on
CHANGE-seq/GUIDE-seq editing data including bulges, correctly collapses the vast majority of these
(the same run reports only **4** sites at CRISPR-Bulge ≥ 0.5). This is exactly why CRISPRme+ ships both
scores and why **bulge off-targets should be prioritized by CRISPR-Bulge, not CFD**.

## Finding 4 — one real (minor) bug found and fixed: INDEL-path N-in-PAM divergence

The INDEL-path twin had drifted from the SNP path in two ways:

- `if "N" == sl: score *= 1` in the mismatch loop — **dead code** (`revcom("N")` is `None`, so the key
  build raises `TypeError` and the score is zeroed+broken before this branch runs; both paths already
  returned 0 for an N off-target base), but confusing.
- `if "N" in pam: score *= 1 else: ...get(pam, 0.0)` for the PAM factor — a **real divergence**: an N in
  the scored 2 bp PAM kept a non-zero CFD on the INDEL path, whereas the SNP path zeroes a
  non-canonical PAM via `pam_scores.get(pam, 0.0)` (the issue-#94 guard). Rare on a clean hg38 with an
  NRG PAM (the N wobble base is excluded by `[-2:]`), but a true cross-path inconsistency.

**Fix.** `analisi_indels_NNN.calc_cfd` was reconciled to be byte-for-byte equivalent to
`new_simple_analysis.calc_cfd` (an N in the genome is an unknown base → zeroing is both consistent and
conservative). A regression test (`test_calc_cfd_twins_agree.py`) exec's `revcom`+`calc_cfd` from both
modules and asserts identical raw-double output over 50,000 random + edge inputs (N-in-PAM,
N/IUPAC-in-off-target, bulge gaps, non-canonical PAM), so the twins cannot silently diverge again.

> Not fixed (deliberately): mechanism B (terminal-bulge free pass) is **documented, not changed** —
> altering it would silently deviate CRISPRme+'s CFD from the field-standard CRISPOR convention in an
> IND-relevant tool. The correct mitigation is to rank bulge off-targets by CRISPR-Bulge.

---

## Verdict

- **CFD core: no bug.** Bit-identical to the canonical reference; reported values faithful.
- **"4mm+1bulge, CFD>0.5" is true, not a bug** — driven chiefly by (A) Doench's fully-tolerated
  mismatch cells (6% of cells = 1.0; present in all 335 sites), with (B) the PAM-distal terminal-bulge
  free pass a minor second-order effect (51/329 sites). Both are intrinsic to CFD, not CRISPRme+
  miscomputations. The DNA-bulge scoring path was adversarially verified (50k synthetic + 3k real
  sites, 0 divergences) to implement Doench's "solution 3" correctly.
- **One minor bug fixed:** the INDEL-path N-in-PAM divergence, now reconciled + regression-tested.
- **Actionable takeaway:** CFD over-nominates heavily-edited/bulge sites (conservative but noisy);
  CRISPR-Bulge collapses them correctly. Prioritize bulge off-targets by CRISPR-Bulge. CFD remains the
  primary score for the mismatch-only, low-edit regime where it is empirically calibrated.

## Reproduce

```
# per-position decomposition of the 4mm+1bulge CFD>0.5 sites + the count table
PostProcess/ : mismatch_score.pkl, PAM_scores.pkl
scripts used in-session: /tmp/cfd_decomp.py (decomposition), /tmp/cfd_stats.py (counts)
# twin equivalence:
python -m pytest PostProcess/test_calc_cfd_twins_agree.py -q
```

# CRISPRme+ 2.8.0 Design Plan — Sites-Only Chimera-Aware Prioritization

> Status: **PLANNING ONLY** (no code). Approved concept by Luca 2026-10-22; targeted for 2.8.0.
> Produced by a 5-agent adversarial exploration (statistics / red-team / implementation / framing →
> synthesis), grounded in the live codebase. This is the sites-only counterpart to the 2.8.0
> context-aware default (per-sample on genotyped indexes).

## The idea (Luca)

When running a **sites-only** database (the "mega" index = union of 1000G-2021 + HGDP + gnomAD v4.1 +
TOPMed + AoU; allele frequencies only, NO genotypes), add a post-processing prioritization step that:
1. **Puts real haplotypes on top** — multi-variant combinations observed as actual *cis* haplotypes in
   the genotyped/phased 1000G+HGDP panel (4,091 samples / 8,182 haplotypes).
2. **Trusts single-variant sites** — which "for sure are not chimeras."

Rationale: get everything (max sensitivity, incl. rare variants absent from 1000G+HGDP), then use the
phased panel to recover real haplotypes; and a single rare variant with AF>0, even if unresolvable to a
named individual, is known to exist in ≥1 person (so real, not chimeric). A "chimera" = a fabricated
*cis* combination of ≥2 variants that co-occur in no real individual.

## 1. Verdict — **tight-with-caveats**

Luca's two load-bearing claims are correct; the implicit chimera-detection half rests on an
ancestry-skewed panel.

- **(a) "single variant ⇒ never a chimera"** — holds **unconditionally** (definitional; a single ALT has
  no co-occurrence to fabricate). Strongest link in the chain.
- **(b) "AF>0 ⇒ observed in ≥1 individual"** — a **near-identity**, not a probabilistic leap, given a
  MAF/quality floor. The shipped mega's **MAF>0.001 filter** supplies exactly that floor (≥1,600-allele
  margin in gnomAD) → removes the singleton worry. Honest statement: *existence in the source cohort,
  not nameability, modulo the source's QC*.
- **(c) "observed-cis-in-panel ⇒ real"** — a **sufficient (not necessary)** condition, enforced by the
  existing phased-same-slot CONFIRMED gate.

**Three caveats (in priority order):**
1. **Ancestry bias (most severe).** Phased panel superpopulations: AFR 986 / AMR 549 / CSA 782 / EAS 816
   / EUR 771 / **MID 157 / OCE 30**. A real *cis*-haplotype private to an under-represented group looks
   "not confirmable" and would be deprioritized *identically to a true chimera* — a false-negative that
   harms under-represented patients. The aggregate sources carry no phase, so they can't rescue a
   thin-ancestry *cis* case. **FATAL for an IND if implicit; MITIGABLE** with a named-count disclaimer +
   a never-penalize (deprioritize-never-drop) rule.
2. **Single-variant artifact floor.** The anti-chimera guarantee says nothing about *call validity*
   (reference error, normalization/strand, per-cohort noise). MAF>0.001 makes a pure artifact
   implausible but not impossible. Scope the "trust single variants" claim to *"not a fabricated
   cis-combination"* only; call quality is an orthogonal axis.
3. **Informative vs uninformative zero.** A zero co-occurrence count is chimera evidence **only** when
   *every* participating variant is itself in-panel (tier C1). If any is out-of-panel (C2), the zero is
   ascertainment missingness, not repulsion. Even a C1 zero is informative only when expected count
   `E = 8182·∏pᵢ ≫ 1`; when `E ≪ 1` a zero is expected and carries no information. Rare *cis* calls carry
   ~6-7% phasing-switch-error risk at AC=2-3 (~35% for singletons).

## 2. Model — discrete tiers A/B/C1/C2 (NOT a continuous headline score)

No calibration corpus exists for "true vs chimeric haplotype," so a continuous headline score would be
uncalibrated (violates the tri-signed risk-card precedent). The truth is categorical. Tiers also make
"NOT dropped" visually obvious.

| Tier | Definition | Meaning | Rank |
|------|-----------|---------|------|
| **A** | k≥2, exact *cis* combo observed in ≥1 of 8,182 phased haplotypes | CONFIRMED cis — exists; named panel carriers + joint AF=count/8,182 | top |
| **B** | k=1, AF>0 above MAF floor | single variant — real in ≥1 individual, non-chimeric by construction | next |
| **C1** | k≥2, **all** constituents in-panel, cis count **=0** | LIKELY CHIMERA — looked where confirmable, absent | deprioritized |
| **C2** | k≥2, **≥1** constituent **not** in panel | UNCONFIRMABLE — absence of evidence ≠ evidence of absence; ancestry-safe | mildly deprioritized |

A/B/C1/C2 **extend, never fork** the shipped CONFIRMED/PUTATIVE vocabulary (A = CONFIRMED; B/C1/C2 refine
PUTATIVE). Rigor requirements: A has a phase sub-state (phased-cis → A-confirmed; unphased co-occurrence
→ A-putative-in-panel; count==1 → "single-haplotype support"); C1≠C2 never collapsed; concordance from
the **observed joint haplotype count**, never marginal allele presence (catches repulsion). Optional C1
diagnostic `E = 8182·∏pᵢ` with the **E≥3 gate** (below it: "panel underpowered"), captioned
"independence assumption; LD not modeled." No hard calibrated chimera probability.

## 3. Algorithm — one mega search → decompose → classify → cis-query → tier → re-rank

**Detection is unchanged; only ranking changes.** Every off-target stays in the set and in a
downloadable tier.
1. Search (unchanged) — single sites-only mega search, worst-case PUTATIVE representative as today.
2. Decompose each off-target into its `(pos, alt)` variant tuples.
3. Classify by k: k=1 → **B** (no panel query); k≥2 → continue.
4. Cis-query the external phased store (the one new step — a cross-index JOIN): all constituents
   in-panel? → look up observed exact-cis phased haplotype count → count>0 **A** (+phase sub-state);
   count==0 **C1** (+diagnostic E, E≥3 gate). Any constituent out-of-panel → **C2** (leave at rung-3
   min-AF; do not interpret the zero).
5. Tier/score: append `Variant_Combination_Tier` + `Tier_Basis` evidence string.
6. Re-rank — **stable, never a filter** (byte-identical off-target set pre/post; only order changes).

## 4. Implementation — reuse ~80%, one new pure module

- **New: `PostProcess/haplotype_prioritizer.py`** (pure, styled like `observed_haplotypes.py`):
  `decompose_off_target`, `classify_combo`, `prioritize_off_targets`.
- **Reuses:** `tier1_genotypes.py:460-600` (`GenotypeReader` bisect JOIN oracle), `observed_haplotypes.py`
  (cis-resolution), `indel_snp_cis.py:43-99` (phased-same-slot CONFIRMED gate) + `:163-174`
  (`putative_joint_af` for C2), variant provenance `new_simple_analysis.py:500-525` →
  `resultIntegrator.py:620-750`, sort `generate_report.py:3000-3032`.
- **Deliver the tier as a joinable companion TSV** (keyed like `phase_confirmation_companion.py`), **NOT**
  a trailing bestMerge column — `change_headers_bestMerge.py` silently drops unknown columns and scorers
  address the tail by negative index (hard constraint).
- **Data input (HYBRID two-source requirement — Luca 2026-10-22):** chimera-aware mode is inherently a
  hybrid op — **detection** uses the sites-only **mega** (genome_library + registry); **cis-confirmation**
  needs the **1000G+HGDP phased genotype store** as the JOIN oracle. It does **NOT** need the full
  HGDP1kGP *search* index (the ~19 GB enriched genome) — only the **genotype store + registry**
  (`Dictionaries/genotypes_hg38_HGDP1kGP` ~28 GB + `registry_hg38_HGDP1kGP` ~8 GB + indel genotypes if
  indel combos). Download shapes:
  1. **Lean "genotype-oracle" artifact (RECOMMENDED):** ship just the phased store + registry (~36 GB) via
     `download --what genotype-oracle --panel HGDP1kGP` — cheapest for a mega-only user (~64 GB total vs
     ~78 GB for two full indexes).
  2. **Reuse full HGDP1kGP index if already present** — default users already have it; chimera-aware mode
     reads its store with no extra download.
  3. **Download-both UX:** `download` for mega chimera-aware mode auto-requires + fetches the oracle so a
     user can't end up half-equipped.
  Resolve the oracle path via manifest sidecar / `CRISPRME_GENOTYPE_TIER` env / conventional path.
  **Fallback when the oracle is absent: no-op WITH a clear message** ("install the HGDP1kGP genotype
  oracle to enable chimera-aware prioritization") — the report is still correct, just un-tiered
  (today's worst-case PUTATIVE). Ties to the existing registry `data_type=hybrid` manifest stamp.
- **Gating:** default OFF; `--prioritize-observed` / web checkbox; **sites-only indexes only** (on a
  genotyped index the 2.8.0 `--per-sample` default already suppresses the chimera at search time).
- **Precondition:** assert mega `bcftools norm -m -any` normalization parity with the panel before the
  join.

## 5. Report / messaging

- Columns (companion): `Variant_Combination_Tier` ∈ {A,B,C1,C2} + `Tier_Basis` evidence string; for A,
  joint AF=count/8,182 + phase sub-state + panel carrier ids; for C1, diagnostic E.
- Stable re-rank with a visible tier chip; **C1/C2 greyed, not hidden** (deprioritized-not-dropped).
- **Honest labels (non-negotiable):** rename user-facing "CONFIRMED" → **"cis-observed in HGDP+1kGP panel
  (n=4,091)"**; three-state visual (`cis-observed (panel)` / `single-variant (not a cis-combination)` /
  `cis-unconfirmable — not in this 4,091-sample panel; MAY be real in an under-represented ancestry OR a
  chimera`). Carriers are **panel individuals, not the patient**. Frequency and chimera-status are
  orthogonal axes.
- **Orthogonal to the calibrated CFD/CRISPR-Bulge expected-count card** — a provenance/plausibility
  badge, never blended into one number. Symmetric story: genotyped index → `--per-sample` (chimera never
  emitted); sites-only → worst-case emitted then A/B/C1/C2 re-ranks via the panel oracle.

## 6. Risks + mitigations
See caveats §1. Summary: ancestry bias → deprioritize-never-drop + named-count disclaimer + C2 only
mildly deprioritized; artifact floor → scope the B claim + optional QC/segdup-CNV axis; over-trust →
rename + confidence-by-AC; LD/E → observed count + E≥3 gate; C1/C2 collapse → two-sided gate; buried =
suppressed → byte-identical-set test; normalization mismatch → parity assertion.

## 7. Work breakdown (rough)
- P0 Decisions + methods pre-write (S) — resolve §8, write definitions + ancestry rationale, tri-sign.
- P1 `haplotype_prioritizer.py` module (M).
- P2 Cross-index JOIN + data wiring + no-op fallback + normalization-parity assertion (M).
- P3 Companion emitter (S).
- P4 Report integration + stable re-rank + honest-label rename + byte-identical-set test (M).
- P5 Gate flag + sites-only gating (S).
- P6 QC/artifact + segdup/CNV flag axis (M→L, optional; the only part on single-variant call validity).
- P7 Tests + real-data validation on mega (M).
- P8 Release-surface consistency + ship (S).
Total ≈ 5 S + 4 M + 1 optional L.

## 8b. DECISIONS (Luca, 2026-10-22)

- **Q1 — sort precedence: DECIDED.** Keep the calibrated CFD/CRISPR-Bulge score as the PRIMARY ranking.
  A possible chimera must NOT sit on top. Validation leads with high-score **AND** observed-real sites.
  Refinement agreed: **C1 (likely chimera: all variants in panel, combo unseen) → HARD demote to bottom
  (demote-not-drop, as shipped for ≥2-bulge); C2 (unconfirmable: out-of-panel variant) → stays in score
  order / only a light nudge, clearly labeled "couldn't confirm"** — because C2 may be a real haplotype in
  an under-represented ancestry (OCE n=30 / MID n=157) and hard-demoting it would recreate the ancestry
  bias. A (confirmed cis) + B (single real variant) stay in calibrated score order on top.
- **Q2 — entry point: DECIDED = Option B (demote re-rank) for 2.8.0.** The companion TSV is the under-the-
  hood plumbing (bestMerge format can't take new columns without breaking scorers), NOT a user choice.
  Visible re-rank (C1 hard-demote, C2 soft) is required to achieve "chimeras off the top." Byte-identical
  off-target SET guard still mandatory (nothing dropped; only order).
- **Q6 — ancestry disclaimer: DECIDED = static warning only.** Name the thin-coverage groups (OCE n=30,
  MID n=157); no automatic ancestry-aware flagging. "Caveat; little we can do about panel coverage."

## 8. Open questions still for Luca
1. **Sort precedence:** tier **primary** ("put real on top": A>B>C1>C2 then detectability within tier),
   or tier as a **tertiary tiebreak** under the calibrated CFD/CRISPR-Bulge sort? (Materially different
   top-1000 orderings. Luca's verbatim "PUT ON TOP" reads tier-primary; the analyses lean badge +
   detectability-primary to preserve the calibrated ranking.)
2. **Entry point:** in-report re-rank (higher impact, touches generate_report.py) vs companion-only
   enrichment first (no re-rank), deferring re-rank to 2.8.x?
3. **E≥3 chimera-flag gate** — bless the threshold (parameter, not magic constant).
4. **Does C1 get a visible "likely chimera" word, or only deprioritization?** (How assertive in an
   IND-facing report.)
5. **QC/artifact axis (P6) in 2.8.0 or deferred?** (Only part addressing single-variant call validity.)
6. **Ancestry disclaimer:** automatic known-patient-ancestry → panel-coverage flag (needs user input) vs
   static methods-doc disclaimer naming OCE n=30 / MID n=157?
7. **Phase-unconfirmed A / count==1 A** — annotated sub-states rather than full-weight A? (Confirm.)

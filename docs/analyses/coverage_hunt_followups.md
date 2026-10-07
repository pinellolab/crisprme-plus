# Coverage hunt — "declared-but-dead" bug class (follow-ups)

An adversarial coverage workflow (2026-10-07, prompted by the IntOGen-column bug: a declared report
column with no populating code that silently rendered `-` in every row) hunted the same class across
report / scorer / registry-AF / CLI-env / web paths. It found **3 confirmed bugs, 3 likely, 14
coverage gaps**. This file records what was fixed and what remains for review.

## Fixed in v2.7.1

1. **IntOGen_cancer_driver `-` in every row** — missing `intogen` branch in
   `generate_report._curated_cell`; + a guard that every `CURATED_COLUMNS` kind has a branch. (#58)
2. **`generate-report --no-maf` silently ignored** by the top-level CLI wrapper — now threaded into
   `build_report(drop_maf=...)` + help + test.
3. **Per-cluster web view/download dropped all 24 CRISPR-Bulge columns** —
   `change_headers_bestMerge.py`'s pre-CRISPR-Bulge `new_order` silently discarded the block; now
   preserved + `test_change_headers_bestmerge.py`.
4. **Reset button never cleared the threshold dropdown** — `Output("thresh_drop", "value ")`
   trailing-space typo; fixed + `test_dash_output_props.py` guards all `pages/` props.

## Deferred — need careful review / real-data validation (NOT rushed into v2.7.1)

### A. Tier-1 genotype axis not filtered to VCF-genotyped samples (#46 AF-deflation on the co-occurrence path)
`build_dictless_tiers.py:237` / `tier1_genotypes.py:410` — the #46 phantom-sample filter
(`genotyped_samples`) is wired on the **Tier-0** path but the **Tier-1** genotype axis
(`build_sample_meta`) is built unfiltered, so an over-listing samplesID can still deflate the joint
AF denominator on the per-sample co-occurrence path. **Correctness-sensitive; needs a fixture-based
AN assertion + real-data check before touching.** Proposed test: `test_tier1_axis_drops_phantom_samplesid`.

### B. Annotation columns can render all-`-` when an annotation source is all-blank
`generate_report.py:3880` — the column-**drop** logic keys on column *name* (always emitted by
resultIntegrator), while the sibling legend logic is value-aware; a fully-blank annotation (e.g.
COSMIC-off default) is kept as an all-`-` column instead of being dropped. Needs a product decision
(drop vs keep-with-note) before changing. Proposed test: `test_annotation_column_dropped_when_all_values_blank`.

## Test-coverage debt (works today, but untested value paths — add tests)

The highest-value gaps, all in the output/scoring path, where a wrong **value** (not structure) would
pass current tests:

- **`resultIntegrator.py` positional `target[N]` projection** (≈45–94 columns, incl. the CRISPR-Bulge
  block at :698) — hard-coded indices with no by-name guard; `intermediate_schema.py` (the by-name
  resolver built to fix this) is itself **unwired**. A shifted upstream column would silently
  mis-populate. Add a fixture asserting a known row projects to the expected named values; wire
  `intermediate_schema`.
- **`add_risk_score.py`** — computes `Highest_CFD_Risk_Score` from fixed `fields[20]/[21]` even when
  scoring the bestCRISPR_BULGE file; no test. Add a value test; confirm the CRISPR-Bulge risk column
  is correct (or intentionally CFD-based).
- **CRISPR-Bulge `-1` sentinel vs genuine all-row degradation** (`generate_report.py:561`,
  `crispr_bulge_score.py:157`) — no test distinguishes the expected 2-bulge/N/over-length NA from a
  silently-degraded all-`-1` scorer; the compute self-test is skipped on the default CPU backend and
  never run on real search output. Add a test asserting 0/1-bulge rows are non-`-1` on a known input.
- **`resolve_phased`** (`tier0_registry.py:1158`) never consults the registry manifest's
  authoritative `phased`/`data_type` field. Add a test.
- **`pam_filter.py:30`** (PAM-satisfiability filter for degenerate/NNN indexes, in the production
  search path) — untested. Add satisfiable/unsatisfiable cases.

> Source of truth: adversarial workflow run `wf_aa3f0bd9-4a4` (2026-10-07). The lesson generalized:
> tests that assert *structure/column-names* are blind to a dead *value*; prefer value-asserting
> tests on representative data for every declared output field/flag.

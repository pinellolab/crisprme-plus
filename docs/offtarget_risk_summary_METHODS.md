# CRISPRme+ Off-target Risk Summary: Methods

*This document specifies the OFF-TARGET RISK card **shipped** at the top of the
variant-search report. It **replaced** the legacy `100/(100+Σ)` aggregate specificity
score there, which is no longer computed or displayed in that report (the two are never
presented together). The separate personal-assembly report is a different output path and
is out of scope here. The card has been empirically bench-tested and adversarially
reviewed. Numbers below are from the real calibration benchmark and a real search; nothing
is illustrative unless labelled.*

*Independent three-way sign-off (2026-10-20): a **statistician**, a **clinician**, and a **genome editor** each reviewed this document and returned **SIGN** after their required fixes were applied (all fixes incorporated below). Reviewers verified the shipped `.npz` artifacts (φ, guide counts, prevalence) match the text.*

*Revised 2026-10-21, after the sign-off above: the deployment prior-shift caveat was made
**qualitative** — the illustrative ×⅓…×3 multipliers (and `offtarget_calibration.prevalence_scaled`)
were removed, because an individual deployed guide's own off-target propensity is unknown and a
specific ratio would imply knowledge we do not have (§5, §6, §7). The three reviewers
**(statistician, clinician, genome-editor) RE-CONFIRMED their SIGN on this revision (2026-10-21)**
after a second-round fix (the §6 magnitude reframed as a cross-guide population envelope, distinct
from a per-guide multiplier; the §7 interval corrected to [2, 81]).*

*Revised 2026-11-22 (2.10.1), after an audit of the risk-card modal and an adversarial
re-verification of that audit. **The 2026-10-20 sign-off line above says the reviewers
verified that the shipped `.npz` artifacts match this text. For the CFD artifact that
stopped being true when its score column was refitted and this document was not.** Every
CFD φ, the decision-band φ, the ≥0.2 aggregate φ, the per-held-out-guide coverage table,
the §4 global LOGO residual, the §4 CFD [0.2,0.5) band error and the §7 dominant-guide φ
move were the superseded pre-2.7.11 artifact's; all are now the shipped ones, each named
as corrected where it appears, and
`PostProcess/test_methods_calibration_numbers.py` re-reads them from the `.npz` in CI so
the drift cannot recur. Also corrected in this revision: the benchmark's positives include
each guide's own on-target while the report's ranked list excludes them (§3, now disclosed
with the deployment-matched values rather than silently restated); the two CFD panel arms
are measured rather than asserted (§3, §3.5); the 95%→99% list ratio (§3.5); "held-out"
applied to the CFD column, which has no training set (§3); 287/359 printed as sites rather
than positives (§3); the "1–3 guides" decision-band support, which matched no artifact
(§5, §8); and the "≈ +0.66" cross-band correlation (§7). What was **not** recomputed in
this pass, and is backlog item 45: the three challenger columns of §4's bake-off table, the
§7 ACTGA/NRG worked-example counts, and a refit of the page on the deployment-matched
population.*

---

## 1. What the summary reports (and what it does not)

For a given guide + PAM + search parameters, the summary reports (**per scorer, never
blended into one number**) the count of candidate off-target sites in each score band
and a calibrated **expected number of true off-targets** with an honest uncertainty
range. The two scorers answer different questions and are shown side by side:

- **CFD** (Doench et al., *Nat. Biotechnol.* 2016), a broad screen: it flags nearly
  every real off-target but most flags are false positives.
- **CRISPR-Bulge** (Yaish & Orenstein, *NAR* 2024), a prioritizer: when it scores a
  site high the site is usually real, but it confirms only about half of **assay-detected**
  off-targets (see §2 for what "real" means) and is defined only on ≤1-gap alignments.

**The expected count is an estimate of how many of the enumerated candidates are *true*
off-targets. It is NOT an upper bound, a ceiling, or a shortlist.** Every enumerated
in-band site remains an assay candidate regardless of the expected count: a row that reads
"217 sites | expected true ≈6 | [1, 14]" does **not** mean "assay ~6–14 of the 217"; it
means that among the 217 enumerated sites an estimated ≈6 (plausibly 1–14) are true
off-targets, and all 217 are still candidates to be assayed. The calibrator re-ranks and
quantifies; it never licenses dropping enumerated sites.

The calibrated numbers apply to **reference-genome, ≤1-gap off-targets** (the domain
of the validation data). In CRISPRme terminology **"gap" means *bulge*, not total edit
distance**: a ≤1-gap site may carry up to three mismatches *and is calibrated*. The
calibration therefore covers **0–1 bulge with any number of in-domain mismatches**; only
**multi-bulge (≥2-gap)** sites are excluded. Both scorers are **sequence-based and
variant-agnostic** (CFD is a function of the aligned target sequence, CRISPR-Bulge of the
guide + target sequence), so the score→probability calibration applies to **any ≤1-gap
alignment, whether the target carries a reference or a variant (ALT) allele**: a
variant-allele off-target with CFD 0.6 is biophysically equivalent to a reference one with
CFD 0.6, and the detection probability is a property of the sequence, not of how the
sequence arose. Variant-/ALT-allele ≤1-gap sites are therefore **included** in the
calibrated expected-count. What is *not* calibrated is (a) **≥2-gap** sites (outside the
scorer domain; the large majority of the raw enumeration, rarely real as each extra bulge
sharply reduces cleavage; zero benchmark support), and (b) the **per-sample
carrier-incidence** dimension: the expected count is the expected number of true
off-target *sequences*, and for a variant site that presumes the individual carries the
variant; which of the panel's haplotypes carry each one, and how many, is a separate,
genotype-dependent axis the score does not address. (Score-conditional transfer presumes a
variant-allele site's score lies within the benchmark's score-feature support; the §6
out-of-support flag is applied **per site (to ALT alleles as to any guide)**, so an ALT
allele that pushes a site into an unusual score regime is flagged. "Score-conditional" is
not "support-free.") This is stated so the clean calibrated
numbers are not read as covering the carrier dimension, and so "≤1-gap calibrated" is not
misread as "only near-perfect matches calibrated."

## 2. Validation benchmark

Calibration and all performance numbers use a benchmark of experimentally measured
off-targets assembled from four published assays (CRISPR-Bulge repository,
OrensteinLab; Zenodo 10.5281/zenodo.10902445):

| source | confirmed positives |
|---|---|
| Tsai et al. 2015 (GUIDE-seq) | 146 |
| Refined TrueOT | 99 |
| Chen et al. 2017 | 66 |
| Listgarten et al. 2018 | 48 |
| **total** | **359** across **46 guides**, screened over **2,940,588** candidate sites |

Prevalence ≈ 1.22×10⁻⁴.

**Assay-sensitivity floor (what "real off-target" means).** The benchmark positives are
**assay-*detected*** off-targets (GUIDE-seq/Tsai, Chen, Listgarten, Refined TrueOT). No
assay recovers all true double-strand breaks (GUIDE-seq alone recovers on the order of
~30–50% of true DSBs), so throughout this document **"real off-target" means "detected by
at least one of these assays,"** and the "negatives" conflate true negatives with assay
misses. Undetected-but-real sites therefore inflate the apparent miss rate: this biases
CRISPR-Bulge recall (the "confirms only about half," §1/§3) *optimistically* with respect
to true recall and CFD precision *pessimistically*. **A recall of 0.50–0.60 must not be
read as calibrated against ground truth**; it is recall against the detected set.

Every confirmed positive is a reference-genome, ≤1-gap site (297 with 0 bulges, 62 with 1
bulge; zero ≥2-gap). As noted in §1, "≤1-gap" is a bulge count, not an edit distance: the
62 single-bulge positives carry 1–3 mismatches in addition to their single bulge and are
fully in the calibrated domain.

The positive set is **guide-clustered**: **one guide contributes 79 of the 359 positives
(79/359 = 22.0%)**, and 16 guides contribute exactly one. This single dominant cluster is
the most influential element of the calibration; its leave-that-guide-out effect is
reported as a full result in §5 and §7, not as a parenthetical. CRISPR-Bulge is evaluated
only on guides held out of its training and on ≤1-gap sites (287 positives / 40 guides,
bulge split 240 + 47; CFD's support is all 359, bulge split 297 + 62), so its table is on
a different support than CFD's and is footnoted as such.

## 3. Discrimination (precision–recall)

AUPR is primary (imbalance-appropriate). On the external benchmark, with
guide-cluster-bootstrap 90% CIs (300 guide-resamples): **CFD AUPR 0.168 (90% CI
0.13–0.24)**, **CRISPR-Bulge AUPR 0.613 (90% CI 0.56–0.70)**. CRISPR-Bulge's is measured on
the 40 guides held out of its training; **CFD has no training set, so nothing is held out
for it** — its support is all 46 benchmark guides, the 6 CRISPR-Bulge trained on included.
"External" is the word that is true of both columns; "held-out" is true only of
CRISPR-Bulge's, and the card, the figures and this section now say which.

**The benchmark counts each guide's own on-target as a positive; the report's ranked list
does not.** 46 of CFD's 359 positives and 40 of CRISPR-Bulge's 287 are perfect (0-mismatch,
0-gap) matches — one per guide, CFD exactly 1.000 on every one of them, CRISPR-Bulge
0.998–1.000 — while the card's score arrays are built as "off-targets (exclude perfect)"
(`generate_report.py`). Every discrimination number in this section is therefore measured
on a population the deployed ranking does not contain. Re-running the same estimators on
the same benchmark with those rows dropped gives:

| quantity | as published | deployment-matched (on-targets dropped) |
|---|---|---|
| CFD AUPR (trapezoid, the artifact's rule) | 0.168 | **0.025** |
| CRISPR-Bulge AUPR (tie-grouped, full precision — the figure's 0.62) | 0.622 | **0.521** |
| CRISPR-Bulge AUPR (trapezoid, 3-decimal column — the table's 0.613) | 0.613 | **0.511** |
| CFD bulge=0 AP | 0.272 | 0.076 |
| CFD t≥0.5 precision / recall | 0.085 / 0.276 | **0.048** / 0.169 |
| CFD t≥0.2 precision / recall | 0.0073 / 0.721 | 0.0060 / 0.681 |
| CRISPR-Bulge t≥0.5 precision / recall | 0.728 / 0.495 | 0.658 / 0.413 |
| CRISPR-Bulge t≥1.9×10⁻⁴ precision / recall | 0.0449 / 0.951 | 0.0386 / 0.943 |
| guides carrying a decision-band (≥0.5) positive | 46 / 40 | 19 / 19 |

CFD moves furthest because its 46 on-targets tie at exactly 1.0 at the head of the ranking:
that single tie group (51 rows, 46 of them positive — precision 0.902 at recall 0.128)
carries **0.116 of the 0.168, 69% of the headline**, before CFD ranks one off-target. The
effect on the card's *expected counts* is small: refit the isotonic map without those rows
and the expected total CFD would print moves ×1.08 in the ≥0.5 band and ×1.02 over ≥0.2
(×1.00 in 0.2–0.5), while **no CRISPR-Bulge band moves at all** — on-targets populate only
the top bins. This is a discrimination-claim problem, not an expected-count problem. **The published values are kept in this release** so they stay
comparable with the benchmark literature, and the deployment-matched values are printed
beside them here and in the modal; re-fitting and re-rendering the page on the
deployment-matched population is backlog item 45, not a release fix.

**Caveats on these headline numbers (required reading alongside them).** Both AUPRs are
point estimates from clustered positives (359 / 46 guides for CFD, 287 / 40 for CB) in
which a single guide supplies 79 positives — 22% of CFD's, and 27.5% of CB's held-out 287,
the same guide in both; AUPR is high-variance and
cluster-dependent, so each number carries a wide **guide-cluster-bootstrap 90% CI**
(reported above; resampling guides, not sites), CFD's spans 0.13–0.24, and the point
value must always be read with it. The **no-skill AUPR baseline equals the prevalence,
≈1.22×10⁻⁴**: CFD's 0.168 and CB's 0.613 are to be read against that floor, not against 0.5. **The two AUPRs are not
comparable across scorers**: they are computed on different support (CFD 359 positives /
46 guides vs CB 287 / 40), so the 0.168-vs-0.613 gap partly reflects different positive
sets and held-out populations, not purely discrimination.

**The same caveat binds harder on the bulge-stratified APs, and hardest on bulge=1.** The
figure's CRISPR-Bulge bulge=1 curve (AP 0.285) rests on **47 positives from 12 of the 40
held-out guides, 25 of them (53%) from one guide**; the CFD bulge=1 curve (AP ≈0.01) on 62
positives from 17 of 46. Restricting the CRISPR-Bulge panels to genuinely held-out guides
made that stratum *more* guide-concentrated, not less (it was 62 positives / 17 guides,
40% from the same guide, while the 6 training guides were still in it). **No bootstrap CI
is computed for the stratified APs**, so their three printed decimals are precision, not
accuracy: 0.285 supports "CRISPR-Bulge ranks single-bulge sites far above CFD's prevalence
floor, on a dozen guides", and must not be read as a point estimate of the same standing
as the headline AUPRs above, which do carry CIs. A single-guide-driven bulge AP is also
the weakest link in §3.5's panel logic — and the union's CFD single-bulge arm does not
repair it: **measured on the 40 held-out guides, CFD ≥0.5 on single-bulge sites adds 614
sites to the sensitive tier and recovers none of the 10 bulge-carrying validated
off-targets that tier misses** (see §3.5 for the full arm-by-arm measurement). The arm is
retained as insurance against the benchmark's own blind spots, not on evidence that it
rescues anything.

AUROC (≈0.97/0.98) is reported only as a footnote: it is prevalence-invariant and does
not reflect operating-point precision at this imbalance.

Operating points (external benchmark; CRISPR-Bulge on its 40 held-out guides, CFD on all
46): CFD t≥0.2 → recall 0.72 / precision 0.007; t≥0.5 → recall
0.28 / precision 0.085. CRISPR-Bulge t≥0.26 → recall 0.60 / precision 0.50; t≥0.5 →
recall 0.50 / precision 0.72. **These CRISPR-Bulge recalls (0.50–0.60) are recall against
the assay-detected set (§2), not against ground truth, and a low CRISPR-Bulge count is not
clearance (§7).** These curves (overall and stratified by bulge count) are shown in the modal.
The modal's curve is titled **AUPR 0.62, not 0.613**, on the same 287 held-out positives
(2,703,869 sites): the
figure re-scores CRISPR-Bulge at full precision and groups tied scores, while the
calibration artifact the headline number and the card's table come from is fitted on the
report's 3-decimal CRISPR-Bulge column (0.613 trapezoid / 0.615 tie-grouped on that
column; 0.618 / 0.622 at full precision). Same sites, same model, different score
precision and AP rule — the gap is far inside the 0.56–0.70 CI, and the 3-decimal fit is
the one the deployed probability map is built on. The CRISPR-Bulge operating points just
above are measured on that same 3-decimal column; at full precision the t≥0.5 point reads
precision 0.728 rather than 0.724 (recalls are identical either way).

**CFD position-aware bulge penalty (2.8.0).** CFD now applies Doench's position-weighted
DNA/RNA-bulge factors at each bulge's true spacer position (previously a DNA bulge was realigned
by a 5' trim that discarded the gap, so every DNA bulge — and any bulge at the two undefined
terminal columns — was scored as if absent; a PAM-proximal bulge scored like a perfect match).
This is a **per-site correctness fix, not a bulge-ranking improvement**: mismatch-only and
in-range RNA-bulge scores are unchanged (bitwise), overall AUPR is unchanged (0.168; the
benchmark is mismatch-dominated, 297/359 positives bulge-free), and CFD's bulge discrimination
stays at the floor (bulge=1 AP ≈ 0.01). What it buys is **precision**: removing the
PAM-proximal over-call drops ~9,000 spurious high-CFD bulge sites at t≥0.2 (precision
0.006→0.007) and roughly doubles precision at t≥0.5 (0.05→0.085). It therefore *reinforces*
the §3.5 panel logic (CRISPR-Bulge is the bulge ranker; CFD rescues a single-bulge site only at
high confidence) rather than changing it. A ≥2-bp bulge (opt-in, outside the single-bulge
calibration domain) applies the position factor per looped-out base — conservative (lower CFD),
and flagged here because it is an extrapolation beyond the published matrix.

## 3.5. Recommended default validation panel

CRISPR-Bulge assigns a *non-zero* score to almost every candidate site. On the held-out
benchmark, scored at **full precision**, **no** validated off-target scores exactly 0 — but
neither do the negatives: **99.9% of them sit in a tiny noise floor below ≈1.9×10⁻⁴**. So a
literal `>0` cutoff would list essentially the whole genome, and the earlier report (which
rounded CRISPR-Bulge to 3 decimals) was in effect thresholding at the rounding floor, 5×10⁻⁴.
The sensitive operating point is instead a **principled, documented threshold**: the
CRISPR-Bulge score that recovers **~95% of the assay-detected off-targets** on the 40-guide
held-out benchmark —

> **CRISPR-Bulge ≥ 1.9×10⁻⁴ (the 95%-recall operating point)**

(point 1.9×10⁻⁴; guide-cluster bootstrap 95% CI [5.8×10⁻⁶, 3.1×10⁻³]; selectable to 90% or 99%
via `CRISPRME_CB_RECALL_TARGET`, single source of truth `PostProcess/cb_operating_point.py`). The
threshold is applied to **full-precision** CRISPR-Bulge scores — the report now writes `{:.6g}`,
not 3 decimals, so the sub-5×10⁻⁴ operating point is representable (3-decimal rounding would have
floored it at 5×10⁻⁴ and dropped ~7% of validated off-targets from the sensitive tier). The
practical default panel is a **union**, not a single cutoff:

> **Default panel = (CRISPR-Bulge ≥ 1.9×10⁻⁴, 95%-recall) ∪ (CFD ≥ 0.2 on 0-bulge sites) ∪ (CFD ≥ 0.5 on single-bulge sites) ∪ must-includes**

where *must-includes* are all perfect / near-perfect matches, any site in a functional element
or cancer-driver gene, and — on a genotyped/phased panel — observed *cis* haplotypes with named
carriers.

- **CRISPR-Bulge ≥ the 95%-recall threshold** carries the bulge sites and the ML-ranked core —
  CFD cannot rank bulges (bulge=1 AP ≈ 0.01 vs CRISPR-Bulge 0.285 on the 40 held-out guides —
  47 positives from 12 of them, 53% from one; read it with the support caveat in §3).
  CFD's 2.8.0 position-aware bulge penalty fixes the *per-site* over-call (a PAM-proximal bulge
  no longer scores like a perfect match) but does not make CFD a bulge ranker.
- **CFD ≥ 0.2 on 0-bulge** adds the high-CFD no-bulge sites CRISPR-Bulge ranks low (CFD is
  reliable off-bulge).
- **CFD ≥ 0.5 on a single bulge** adds the *high-CFD* 1-bulge sites: CFD now position-penalizes
  the bulge but is still not calibrated to rank bulges, so only high-confidence calls are added
  (CRISPR-Bulge remains the primary bulge ranker; this only *adds* sites it scored low).
- **≥2-bulge sites** are excluded from scoring-driven selection (uncalibrated for both scores)
  but kept-and-demoted in the report, never dropped.

**What the two CFD arms are worth, measured.** On the 40 held-out guides, with CRISPR-Bulge
at full precision:

| arm | sites added | validated off-targets recovered |
|---|---|---|
| CRISPR-Bulge ≥ 1.9×10⁻⁴ (the base tier) | 6,076 | 273 of 287 (recall 0.9512) |
| + CFD ≥ 0.2 on 0-bulge sites | +9,660 | **+1** |
| + CFD ≥ 0.5 on single-bulge sites | +614 | **0** |
| union | 16,350 | 274 (recall 0.9547) |

The two arms **2.7× the panel for +0.35 pp of recall**. Of the 14 positives the base tier
misses, **10 carry a bulge and the CFD single-bulge arm recovers none of them**; the one
recovery is a no-bulge site. The same answer on all 46 guides (+12,188/+1 and +704/0) and
with the on-target rows excluded. Earlier drafts of this section and of the modal stated
the rescue as a fact; it is not one on this benchmark. **The arms are kept anyway**, for a
reason the table cannot show: a site no assay ever detected is scored here as a false
positive (§2), so 1 in 10,274 is an *upper bound* on what the arms find on the detected
set, not an estimate of what they would find in a new experiment. They are cheap insurance
against the benchmark's blind spots, and the panel text must say so rather than claim a
measured benefit.

**Recall vs. validation size.** Recall above ~95% is reachable but the list grows **steeply**,
because the last few percent of true off-targets sit in the noise floor shared with the
negatives. On the benchmark the trade is ≈**9 sites per true off-target at 90% recall, ≈21 at
95%, ≈354 at 99%**, ≈790 at 100%. In list sizes on the 40 held-out guides, the 99%-recall list
is **≈16× the 95% one — 6,076 sites against 99,477** at the shipped thresholds (16.4× at the
exact recall quantiles, 6,048 → 99,406). Earlier drafts said "≈6×" here; that figure is the
ACTGA/NRG worked example's ratio, not the benchmark's, and the "≈21 → ≈354" in the same
sentence already implies ≈17×. (`PostProcess/cb_operating_point.py`'s docstring still carries a
third version of the same claim, "~6x … 206 → 1,286 loci"; that file belongs to another lane
and is listed for correction there.) For the ACTGA/NRG worked example (189,957
candidate loci) the 95% threshold yields **162 off-target validation loci** (90% → 58; 99% → 1,158; the shipped crispr_bulge_gt_0.tsv tier; +2 perfect-match candidate on-targets), measured
directly on the regenerated full-precision report (not the per-true-OT ratios above, which are
benchmark-wide averages).
Panel depth is a **risk/validation-size optimization** read off the expected-true counts (§5),
**not a recall guarantee** — detection is lossless (every site is listed); only score-ranked
selection carries this trade.

## 4. Point calibration: isotonic, and why

The raw score is a ranking number, not a probability. We map score → P(true off-target)
with **isotonic regression** (monotone, pool-adjacent-violators). **The map that ships is
fit on all benchmark guides at once** — it is a deployment map, and it closes on its own
training data by construction. **Leave-one-guide-out (LOGO) is the validation, not the
fit**: the band-error table below, the overdispersion φ of §5 and the per-held-out-guide
coverage check are each computed from refits in which a guide never informs its own
numbers. Saying the shipped map is "fit leave-one-guide-out" (as earlier drafts and the
modal did) presents an in-sample closure as an out-of-sample guarantee; the out-of-sample
statement is the LOGO total-count residual at the end of this section.
We chose isotonic after a bake-off
against fairly-specified smooth monotone challengers (logistic-on-logit, monotone spline,
beta calibration), judged on **per-band signed expected-count error**.

**Evaluation protocol for the table below (stated explicitly).** The per-band signed
errors are computed **on this single 359-positive benchmark**. *Selection among the four
families is therefore in-sample model selection on one dataset.* Isotonic is the most
flexible (step-function) family and will tend to fit band-level counts best on the same
data used to judge it; if the band errors below are in-sample rather than out-of-fold, the
comparison is **biased toward isotonic**, and that qualification must travel with the
table. [The band-error figures should be reported as LOGO/out-of-fold wherever the fit
pipeline supports it; the global LOGO residuals in the last paragraph *are* out-of-fold.]

| decision band | isotonic | logit-logistic | monotone spline | beta |
|---|---|---|---|---|
| CFD [0.5,1] | **+1.7%** | +10.6% | +10.6% | −19.5% |
| CFD [0.2,0.5) | **+14.3%** | −16.6% | −16.6% | +27.6% |
| CRISPR-Bulge [0.5,1] | **+0.0%** | +14.1% | +14.5% | −6.9% |

*Fit provenance of this table (read before comparing columns).* The **isotonic** column is
the shipped fit, re-derived from `{cfd,crispr_bulge}_calibration.npz`: CFD [0.5,1] +1.68%,
CFD [0.2,0.5) **+14.32%**, CRISPR-Bulge [0.5,1] +0.03%. The CFD [0.2,0.5) cell previously
read +3.6%, which is this band's error on the **superseded pre-2.7.11 CFD score column**
(+3.57% there) — one of a family of pre-2.7.11 values this document carried; the rest are
corrected in §5 and §7. **The three challenger columns have not been recomputed on the
current CFD column** and are retained as the historical bake-off; recomputing them is part
of backlog item 45. The decision-band conclusion below does not depend on them: it rests on
isotonic's own [0.5,1] errors, which are shipped-fit values.

The smooth maps buy a better *global* figure by **under-counting the high-risk decision
band (the unsafe direction)**; beta is rejected specifically because it under-counts the
decision band (CFD −19.5%, CB −6.9%), not for "compensating errors." **On this benchmark
isotonic showed the smallest decision-band count error (+1.7% CFD, +0.0% CB); its sign is
not guaranteed.** A monotone step function is not conservative "by design": the direction
of its top-band bias is data-dependent and happened to be +1.7% on this sample. We
therefore describe isotonic as the best of the fair candidates for *this* data, not as an
absolute ceiling and not as a method with a guaranteed conservative direction.

**Global vs decision-band bias (do not conflate).** Isotonic's *global/total-count* LOGO
residual on the shipped fits is an over-estimate of **+25.0% for CFD and +4.9% for
CRISPR-Bulge** (448.7 expected against 359 observed; 301.0 against 287). Earlier drafts
gave CFD as +3.8%, which is the superseded pre-2.7.11 column's residual (+3.80% there); on
the current column it is six times that. **Neither figure describes the decision band.**
For CFD the +25.0% is almost entirely the screen bands — the 0.1–0.2 band over-predicts by
**+114%** (130.7 expected against 61 observed) and 0–0.1 under-predicts by −11.8% — while
in the safety-relevant decision band the signed band error is **+1.7% (CFD) and +0.0%
(CB)**, i.e. effectively unbiased, and across the whole ≥0.2 aggregate **+9.5% (CFD) and
+1.8% (CB)**. **An IND reviewer must not read the global residual either as a safety margin
or as a defect in the band where a clinician acts: the decision-band bias there is ≈0–1.7%
and the ≥0.2 bias ≈+2–9%. CFD's screen bands below 0.2 are not calibrated to a standard the
decision band meets, and the card does not display them.** The deployment map (fit on all
guides) is in-sample unbiased: summed through the runtime `np.interp` path over every
benchmark row, Σp 358.9 vs 359 for CFD and 284.3 vs 287 for CRISPR-Bulge (exactly 359.000
and 287.000 at the bin level the isotonic fit closes on — the small shortfall is the linear
interpolation between bin centres that `offtarget_calibration` applies at report time).

## 5. Uncertainty: prediction interval (what it covers and what it does not)

The reported range is a **90% prediction interval for the realized number of true
off-targets**, built as a **two-level guide-cluster bootstrap**: resample the 46 (CB: 40)
guides with replacement, **refit** the isotonic map each replicate (map + guide-selection
uncertainty), apply it to the search's sites, add a **negative-binomial / Gamma
between-guide overdispersion** multiplier (φ estimated per band by method of moments),
then draw the realized count (Poisson).

**What this interval does and does not cover.** It captures **within-population
statistical uncertainty only**: calibration-map uncertainty, guide-selection uncertainty,
and between-guide overdispersion, *among guides resembling the 46-guide benchmark
support*. It **does not** cover, and **cannot be validated to cover, a novel clinical
guide**, because the per-guide multiplicative prior shift (§6) is outside the bootstrap.
**The deployment-guide prior shift (§6), not this bootstrap interval, is what addresses a
novel guide's uncertainty; the summary flags it qualitatively rather than with a specific
multiplier, because an individual guide's own off-target propensity is unknown.**

**Overdispersion φ is weakly estimated in the decision band, and the band's positives come
from few guides once the on-targets are counted out.** φ is estimated per band by method of
moments. The shipped values, read from the artifacts, are **CFD φ ≈ 1.47 / 0.98 / 0.89**
for 0–0.1 / 0.1–0.2 / 0.2–0.5 and **0.28 in the decision band**; CRISPR-Bulge's decision
band is **0.01**. (Earlier drafts gave 1.83 / 1.92 / 0.99 and a decision-band 0.19: those
are the superseded pre-2.7.11 CFD artifact's values throughout.)

**Correction to the support claim.** Earlier drafts said the decision band "rests on only
1–3 guides carrying positives". **No artifact, current or superseded, contains 1–3.** The
shipped `band_counts` say the ≥0.5 band carries **99 positives from all 46 guides** (CFD)
and **142 positives from all 40 guides** (CRISPR-Bulge). The band whose support really is thin is
**CRISPR-Bulge 0.2–0.5: 41 positives from 12 of 40 guides**, and the decision band's
support does thin sharply once each guide's own on-target is excluded (§3) — **53 positives
from 19 guides** (CFD ≥0.5) and **102 from 19** (CRISPR-Bulge ≥0.5), i.e. roughly half the
guides, not 1–3. The underlying argument stands on the corrected numbers: a
method-of-moments φ from 19 effective clusters, in a band whose per-guide counts are small
and heavily tied to one dominant guide (§2), has few degrees of freedom and large sampling
error, so the **decision-band interval width remains assumption-driven rather than
tightly estimated**, and the decision-band coverage figures below must not be presented as
strongly data-validated. *The card now states, per band, the number of positives and the
number of guides carrying one (`band_counts`, displayed in the modal); absent a pooled /
hierarchical / shrinkage estimate of φ with a stated prior, the decision-band interval
width stays a modelling assumption.* (Direction of the resulting values: φ is large in the
low/action bands (interval honestly wide where a single guide's rate is least
predictable) and small at the top.)

**Validation: per-held-out-guide coverage (a guide new to the fit, drawn from the
benchmark's own population).** We validate on the nearest thing to the deployment case the
data supports: leave **one** guide out, build *its* per-band 90% prediction
interval from the other 45 (CB: 39) guides, and check whether that held-out guide's **own**
realized true-positive count falls inside. This is the single-guide estimand, not
the aggregate total-count, which would average out exactly the between-guide variation §6
identifies as dominant. **It is not the same as covering a clinically novel guide**: the
held-out guide is still one of the benchmark's 46, so this check bounds map and
guide-selection uncertainty within that population, and the paragraph above — the interval
cannot be validated to cover a novel clinical guide — continues to hold. The modal states
it in those two parts. Averaged over the 46 (CB: 40) folds, realized coverage of the
nominal 90% interval is:

| coverage of nominal 90% | CFD 0.2–0.5 | CFD ≥0.5 | **CFD ≥0.2** | CB ≥0.2 |
|---|---|---|---|---|
| naive Poisson-binomial (1.3–6.3× too narrow) † | ~72% | 83% | n/a | n/a |
| two-level, no overdispersion † | 76% | 87% | n/a | n/a |
| **two-level + overdispersion (shipped)** | **89.1%** | **95.7%** | **95.7%** | **100%** |

† The two retired-estimator rows were measured on the superseded pre-2.7.11 CFD column and
have not been recomputed; they are kept only as the historical comparison. The shipped row
is read directly from `loo_total_count_coverage_nominal90` in the two artifacts. Earlier
drafts printed that row as 95 / 93 / 95 / 100, which is the pre-2.7.11 CFD artifact's
coverage (95.7 / 93.5 / 95.7), not what ships.

**One shipped band falls below nominal: CFD 0.2–0.5, at 89.1% (41 of 46 folds).** The other
bands are **conservative**, 95.7–100%. The earlier claim that "all shipped bands reach ≥90%"
was true of the superseded artifact, not of this one. How to read these numbers:
- They are empirical coverages from only **46 (CB: 40) folds**, binomial SE ≈4.4pp (n=46) /
  ≈4.7pp (n=40) at p = 0.90, so **89.1% is one fold below nominal and is not statistically
  distinguishable from 90%** — two independent reruns of the same check with different
  Monte-Carlo streams put it at 87.0% (40/46) and 91.3% (42/46). Read it as "at nominal,
  within sampling error, and the only band not visibly conservative," not as validated
  ≥90%. The over-90% values indicate a mildly **conservative (slightly too wide)** interval,
  acceptable for an IND but labelled conservative rather than "nominal." **"100%"
  at n = 40** implies an exact one-sided 90% Clopper–Pearson lower bound of **94.4%**
  (95.1% for 46/46), not literal full coverage.
- Coverage is also only **partly** out-of-sample: the fold refits the isotonic map without
  guide *q*, but the φ that sets the interval's width is estimated on all guides including
  *q*. With 40–46 clusters the leakage is small, but "validated per held-out guide" is a
  statement about the map, not about the width.
- Overdispersion φ is well-supported for the **≥0.2 aggregate** (**0.826**, estimated on all
  46 guides) but the **decision-band** φ (**≈0.28 CFD, ≈0.01 CB**) comes from a band whose
  positives, excluding on-targets, sit on 19 guides (§3, above), so the decision-band
  interval width is better read as an **assumption than as a data estimate**, and the
  per-band positive and guide counts are now displayed in the modal.

The naive Poisson-binomial interval (treating sites as independent, 1.3–6.3× too narrow)
under-covers to ~72% and is retired.

**Dominant-guide (22%) sensitivity, reported as a full result.** The single guide
contributing 79/359 = 22% of positives (§2) drives the calibration. Removing it
(leave-that-guide-out) shifts the CFD action-band interval ~20%; the full recomputed
point estimate and interval for the worked example under this guide's removal are given in
§7. The number of benchmark guides carrying **any** decision-band positive (the count on
which both the decision-band calibration and φ hinge) is **46 (CFD ≥0.5) and 40 (CB ≥0.5)
as the artifacts record it, and 19 and 19 once each guide's own on-target is excluded**
(§3) — not the "1–3" earlier drafts gave. It is stated on the face of the decision-band
rows, in the modal's per-band support table.

## 6. Prevalence / prior-shift

All counts are at the benchmark's average off-target prevalence. A deployed guide at its
own locus can be systematically more or less off-target-prone; that multiplicative prior
shift is **the dominant source of absolute-count uncertainty**, larger than the entire
within-population prediction interval of §5. The summary therefore does **not** headline a
bare count: it shows the within-population statistical range (§5) **and** flags, as a
**separate axis**, that this prior shift can move the count substantially in either
direction and is usually the larger uncertainty. As a population-level scale, the plausible
envelope of that shift across guides spans roughly an order of magnitude or more in either
direction — but **we do not assign a specific per-guide multiplier** (an earlier draft used
an illustrative ×⅓…×3): an individual deployed guide's own propensity is unknown, so a point
ratio would imply per-guide knowledge we lack. The caveat is therefore qualitative — it
tells the reader the estimate is not a ceiling and the true count for their guide can lie
well outside the §5 range — and it, not the §5 bootstrap, is the axis that speaks to a novel
deployed guide. The summary also flags a guide whose score distribution is far from the
46-guide support as out-of-support.

## 7. Worked example (real search: ACTGAAATCTGTAAGCAGGC, NRG, HGDP+1kGP, per-sample)

**Scope and site breakdown (read before the table).** This per-sample HGDP+1kGP search
examined **6,151,531 sites**, which split as:

| | sites | calibrated? |
|---|---|---|
| perfect (0-mismatch, 0-gap) on-target matches | 2 | n/a (on-target) |
| ≤1-gap off-targets | 728,505 | **yes** (412,433 reference-origin + 316,072 variant/ALT-allele) |
| ≥2-gap off-targets | 5,423,024 | no (outside scorer domain) |

Because the scorers are sequence-based (§1), **both the reference and the variant-allele
≤1-gap sites are calibrated**; only the ≥2-gap sites (88% of the enumeration) are not. The
rows below therefore cover the full ≤1-gap off-target set (5,671 at CFD ≥0.2 = 217 at ≥0.5 +
5,454 at 0.2–0.5), **reference and variant-allele together**. What they do *not* resolve
is per-sample **carrier incidence**: which of the 4,091 haplotypes carries each variant
site (a separate genotype axis, §1). The 2 perfect matches mean **no single unambiguous
on-target**.

| scorer / band | sites (of which variant-allele) | expected true | 90% interval |
|---|---|---|---|
| CFD ≥0.5 | 217 (109 variant) | ≈6 | [1, 14] |
| CFD 0.2–0.5 | 5,454 (2,548 variant) | ≈21 | [1, 64] |
| **CFD ≥0.2** | 5,671 (2,657 variant) | **≈28** | **[2, 81]** |
| CRISPR-Bulge ≥0.2 | 6 | ≈2 | [0, 6] |

(Intervals are the converged 90% prediction interval, bootstrap n=4000, seed-stable; a
preliminary low-n draw read the ≥0.2 upper bound as ~69; the stable value is ~81, the
heavy between-guide overdispersion tail.)

**~47% of the ≥0.2 sites (2,657/5,671) are variant/ALT-allele**, present only in carriers.
The ≈28 is therefore an expected count of true off-target *sequences* across the **union of
the 4,091 panel haplotypes, not a single-patient quantity.** For a per-product /
per-patient validation panel, the reference-allele sites (≈3,014 at ≥0.2) are present in
every genome; each variant-allele site is present only in individuals carrying that variant
(carrier incidence, §1). The split is shown per band so the carrier-contingent fraction of
the headline is visible, not just stated abstractly.

**Reading notes on this table (each forecloses a misread an IND reviewer will otherwise
make):**

- **The expected count is not a ceiling or a shortlist (§1).** "CFD ≥0.5 | 217 | ≈6 |
  [1, 14]" means ≈6 (plausibly 1–14) of the 217 enumerated sites are estimated true
  off-targets; **all 217 remain assay candidates.** Do not read it as "assay ~6–14 of the
  217."
- **The ≥0.2 point (≈28) is the jointly-computed estimate** over the pooled 5,671 sites
  (sub-bands ≈6 + ≈21); the "≈28" is rounding, not an arithmetic error.
- **The ≥0.2 interval [2, 81] is computed on the pooled ≥0.2 band with its own
  overdispersion (shipped φ ≈ **0.826**; earlier drafts said 0.806, the pre-2.7.11 value),
  not by combining the sub-bands.** Naively summing the
  sub-band draws *as if independent* would wrongly raise the lower bound to ~5, because a
  guide's promiscuity is correlated across bands. Earlier drafts put that correlation at
  "≈ +0.66", which reproduces under no definition we can find; measured on the 46 benchmark
  guides it is **+0.91 Pearson on per-guide positive counts (0.2–0.5 vs ≥0.5), +0.73
  Spearman, +0.79 on the LOGO residuals** (CRISPR-Bulge: +0.90 / +0.62 / +0.79). The
  correlation is stronger than the retired figure claimed, so the argument for pooling is
  stronger too; modelling it (pooled band) puts the floor at **2**, where it belongs.
- **CRISPR-Bulge ≥0.2 | 6 | ≈2 | [0, 6]: a low CB count is NOT clearance.** The lower
  bound of 0 on 6 sites must not be read as "nothing here." Per §1–§3, CRISPR-Bulge
  confirms only ~half of assay-detected off-targets and is defined only on ≤1-gap
  alignments, so **real off-targets can sit entirely outside what CRISPR-Bulge counts.** A
  low CB number means CB did not prioritize these sites, not that the locus is clear.
- **Dominant-guide (22%) sensitivity.** Dropping the single 79-positive benchmark guide
  (§2, §5) and refitting on the other 45 moves the CFD ≥0.2 row from ≈28 [2, 81] to
  **≈23** (~−18% point); on the shipped fit the matching overdispersion move is
  **φ_≥0.2 0.826 → 0.737** (earlier drafts said 0.806 → 0.62, both pre-2.7.11 values). The
  headline does not hinge on that one cluster. *The ≈28 / [2, 81] / ≈23 figures themselves
  are from the ACTGA/NRG search and were not recomputed against the shipped artifact in
  this pass; see backlog item 45.*

The deployment prior-shift axis (§6) applies on top of the CFD ≥0.2 estimate: this guide's
own off-target propensity at its locus is unknown and can move the count substantially in
either direction — usually a larger uncertainty than the statistical range above. The
summary states this qualitatively and attaches no specific multiplier. The two scorer
estimates are different views of the same sites and are **not additive**.

## 8. Limitations (carried into the report)

1. **Information ceiling**: 359 positives in 46 guides (one holding 22%). No calibrator
   manufactures information the scores lack; CFD's AUPR 0.168 (baseline prevalence
   1.22×10⁻⁴ — and **0.025 on the population the report actually ranks**, once each guide's
   own on-target is excluded, §3) makes its bands inherently wide. The honest interval
   width is the deliverable, not a defect.
2. **Top-band resolution**: the decision band carries 99 (CFD) / 142 (CB) positives from
   every benchmark guide, but **only 53 and 102 from 19 guides each once the on-targets are
   excluded** (§3, §5; the "1–3 guides" of earlier drafts matched no artifact). Its point
   is accurate but its interval is wide by construction, and its φ (hence its width) is
   assumption-driven more than data-estimated (§5).
3. **Generalization**: single-guide coverage is validated per-held-out-**benchmark**-guide
   (§5; 89.1% in the CFD 0.2–0.5 band, 95.7–100% in the others, against a nominal 90%, so
   one band sits at nominal and the rest are conservative); a chemically/structurally novel
   clinical guide at a
   different locus nonetheless remains an extrapolation, addressed by the prevalence axis (§6).
4. **Prevalence transfer** (§6) dominates the absolute count — across guides the plausible
   envelope spans roughly an order of magnitude or more — and is larger than the §5 interval;
   no specific per-guide multiplier is assigned (the guide's own propensity is unknown).
5. **Support and assay floor**: calibration applies to any ≤1-gap alignment (0–1 bulge,
   any number of in-domain mismatches), **reference or variant/ALT allele**, since the
   scorers are sequence-based (§1); only **≥2-bulge** sites and the **per-sample
   carrier-incidence** dimension are out of scope. "Real off-target" throughout means "assay-detected"
   (§2): undetected-but-real sites inflate the apparent miss rate, so reported
   recall/precision are relative to the detected set, not ground truth.
6. **CFD on gapped sites** is a heuristic extension of a mismatch model and over-calls
   bulges; the calibration is applied to the same CFD variant the report displays.

## 9. Implementation

- Fit (offline, re-run if the benchmark changes): `scripts/fit_calibration.py` →
  `PostProcess/calibration/{cfd,crispr_bulge}_calibration.npz` (isotonic point map +
  1,000-map guide-bootstrap ensemble + per-band φ + provenance/prevalence metadata; ~1.9 MB).
  The `.npz` contains `grid, p_point, p_ens, phi`, and `meta` carrying **`phi_bands`** (per
  report band incl. the ≥0.2 aggregate), **`band_counts`** (per-band positives + contributing
  guides), **`aupr`** (point + 90% CI), and **`loo_total_count_coverage_nominal90`**
  (per-held-out-guide coverage per band): the quantities §3/§5 require.
- Runtime: `PostProcess/offtarget_calibration.py`, `expected_counts(scores, scorer)`
  applies the point map and the two-level + overdispersion + Poisson bootstrap as a
  cheap 1-D lookup; no fitting at report time.
- The summary card replaces the legacy `100/(100+Σ)` rows (and the legacy number must not
  be shipped alongside this component, see header); the "How to read this" modal
  (reusing the observed-haplotypes modal) holds the PR curves with their bootstrap CIs and
  prevalence baseline, the calibration explanation, the expected-count-is-not-a-ceiling and
  low-CB-is-not-clearance notes, the priority-order guidance, the 2-gap / uncalibrated-
  support note, the on-target disclosure of §3 with its deployment-matched values, and a
  per-band support table rendered from `band_counts` (positives, guides carrying one,
  guides contributing any site) so no band shows an expected count and a range without
  the support behind them.


## 10. Data provenance & references

**Calibration benchmark provenance.** All calibration and performance numbers use the four
external-evaluation datasets distributed with the CRISPR-Bulge repository (Yaish & Orenstein
2024; code + data archived at Zenodo doi:10.5281/zenodo.10902445;
github.com/OrensteinLab/CRISPR-Bulge), under `files/datasets/`: `Refined_TrueOT.csv`,
`Tsai_2015_GUIDE_seq.csv`, `Chen_2017_GUIDE_seq.csv`, `Listgarten_2018_GUIDE_seq.csv`. These
are the repository's held-out evaluation sets; the CHANGE-seq training data is excluded.
Binary labels follow the repository rule: read-count GUIDE-seq files → reads ≥ 100 active,
reads = 0 inactive, 0 < reads < 100 dropped as ambiguous, off-targets containing `N` dropped;
TrueOT → editing rate > 0.1 % active. Each dataset traces to an original study (below).

**References.**
1. **CFD score.** Doench JG, Fusi N, Sullender M, et al. Optimized sgRNA design to maximize
   activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol.* 2016;34(2):184–191.
   doi:10.1038/nbt.3437.
2. **CRISPR-Bulge score & benchmark.** Yaish O, Orenstein Y. Generating, modeling and evaluating
   a large-scale set of CRISPR/Cas9 off-target sites with bulges. *Nucleic Acids Res.*
   2024;52(12):6777–6790. doi:10.1093/nar/gkae428.
3. **Tsai 2015 (GUIDE-seq).** Tsai SQ, Zheng Z, Nguyen NT, et al. GUIDE-seq enables genome-wide
   profiling of off-target cleavage by CRISPR-Cas nucleases. *Nat Biotechnol.* 2015;33(2):187–197.
   doi:10.1038/nbt.3117.
4. **Chen 2017 (high-fidelity Cas9 / GUIDE-seq).** Chen JS, Dagdas YS, Kleinstiver BP, et al.
   Enhanced proofreading governs CRISPR-Cas9 targeting accuracy. *Nature.* 2017;550(7676):407–410.
   doi:10.1038/nature24268.
5. **Listgarten 2018 (Elevation).** Listgarten J, Weinstein M, Kleinstiver BP, et al. Prediction
   of off-target activities for the end-to-end design of CRISPR guide RNAs. *Nat Biomed Eng.*
   2018;2(1):38–47. doi:10.1038/s41551-017-0178-6.
6. **TrueOT (refined as `Refined_TrueOT`).** Kota PK, Pan Y, Vu H-A, Cao M, Baraniuk RG, Bao G.
   The Need for Transfer Learning in CRISPR-Cas Off-Target Scoring. bioRxiv 2021.08.28.457846
   (2021); dataset github.com/baolab-rice/CRISPR_OT_scoring.

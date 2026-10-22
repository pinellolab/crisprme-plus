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
only on guides held out of its training and on ≤1-gap sites (287 positives / 40 guides;
bulge split 297+62 = 359 for CFD support), so its table is on a different support than
CFD's and is footnoted as such.

## 3. Discrimination (precision–recall)

AUPR is primary (imbalance-appropriate). Held-out, with guide-cluster-bootstrap 90% CIs
(300 guide-resamples): **CFD AUPR 0.168 (90% CI 0.13–0.24)**, **CRISPR-Bulge AUPR 0.613
(90% CI 0.56–0.70)**.

**Caveats on these headline numbers (required reading alongside them).** Both AUPRs are
point estimates from clustered positives (359 / 46 guides for CFD, 287 / 40 for CB) in
which a single guide supplies 22% of CFD positives; AUPR is high-variance and
cluster-dependent, so each number carries a wide **guide-cluster-bootstrap 90% CI**
(reported above; resampling guides, not sites), CFD's spans 0.13–0.24, and the point
value must always be read with it. The **no-skill AUPR baseline equals the prevalence,
≈1.22×10⁻⁴**: CFD's 0.168 and CB's 0.613 are to be read against that floor, not against 0.5. **The two AUPRs are not
comparable across scorers**: they are computed on different support (CFD 359 positives /
46 guides vs CB 287 / 40), so the 0.168-vs-0.613 gap partly reflects different positive
sets and held-out populations, not purely discrimination.

AUROC (≈0.97/0.98) is reported only as a footnote: it is prevalence-invariant and does
not reflect operating-point precision at this imbalance.

Operating points (held-out): CFD t≥0.2 → recall 0.74 / precision 0.006; t≥0.5 → recall
0.29 / precision 0.05. CRISPR-Bulge t≥0.26 → recall 0.60 / precision 0.50; t≥0.5 →
recall 0.50 / precision 0.72. **These CRISPR-Bulge recalls (0.50–0.60) are recall against
the assay-detected set (§2), not against ground truth, and a low CRISPR-Bulge count is not
clearance (§7).** These curves are shown in the modal.

## 3.5. Recommended default validation panel

CRISPR-Bulge is a high-precision screen with a **hard recall ceiling**: it scores ~8% of the
assay-detected off-targets as *exactly* 0, so no single CRISPR-Bulge threshold is both
high-recall and manageable. A cutoff of `>0` already captures essentially everything the model
can see (~92% recall against the assay-detected set); chasing higher recall collapses the
threshold toward 0 and returns the whole genome. The practical default panel is therefore a
**union, not a cutoff**:

> **Default panel = (CRISPR-Bulge > 0) ∪ (CFD ≥ 0.2 on 0-bulge sites) ∪ must-includes**

where *must-includes* are all perfect / near-perfect matches, any site in a functional element
or cancer-driver gene, and — on a genotyped/phased panel — observed *cis* haplotypes with named
carriers.

- **CRISPR-Bulge > 0** carries the bulge sites and the ML-ranked core — CFD cannot rank bulges
  (bulge=1 AP ≈ 0.02 vs CRISPR-Bulge ≈ 0.25; §3).
- **CFD ≥ 0.2 on 0-bulge** rescues the high-CFD no-bulge sites CRISPR-Bulge misses (CFD is
  defensible only off-bulge).
- **≥2-bulge sites** are excluded from scoring-driven selection (uncalibrated for both scores)
  but kept-and-demoted in the report, never dropped.

On the benchmark this union reaches **~93% recall** against the assay-detected set; the CFD
rescue adds only a few sites, and **~7% (dominated by bulge sites) score below both scores** and
are recovered only by the must-include criteria or by validating deeper. Panel depth is a
**risk/budget optimization** read off the expected-true counts (§5), **not a recall guarantee** —
detection is lossless (every site is listed), but score-ranked selection carries this ceiling.

## 4. Point calibration: isotonic, and why

The raw score is a ranking number, not a probability. We map score → P(true off-target)
with **isotonic regression** (monotone, pool-adjacent-violators), fit **leave-one-guide-out**
(LOGO; a guide never informs its own probabilities). We chose isotonic after a bake-off
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
| CFD [0.2,0.5) | **+3.6%** | −16.6% | −16.6% | +27.6% |
| CRISPR-Bulge [0.5,1] | **+0.0%** | +14.1% | +14.5% | −6.9% |

The smooth maps buy a better *global* figure by **under-counting the high-risk decision
band (the unsafe direction)**; beta is rejected specifically because it under-counts the
decision band (CFD −19.5%, CB −6.9%), not for "compensating errors." **On this benchmark
isotonic showed the smallest decision-band count error (+1.7% CFD, +0.0% CB); its sign is
not guaranteed.** A monotone step function is not conservative "by design": the direction
of its top-band bias is data-dependent and happened to be +1.7% on this sample. We
therefore describe isotonic as the best of the fair candidates for *this* data, not as an
absolute ceiling and not as a method with a guaranteed conservative direction.

**Global vs decision-band bias (do not conflate).** Isotonic's *global/total-count* LOGO
residual is a small over-estimate (**+3.8% CFD, +4.9% CB**), i.e. mildly conservative in
aggregate. **This +4–5% figure is a global, total-count number and does NOT describe the
decision band.** In the safety-relevant decision band the signed band error is **+1.7%
(CFD) and +0.0% (CB)**, i.e. **effectively unbiased, not +4–5% conservative**. The
aggregate +3.8%/+4.9% is dominated by the low/screen bands. **An IND reviewer must not
infer a +4–5% safety margin in the band where a clinician acts; the decision-band bias
there is ≈0–1.7%.** The deployment map (fit on all guides) is in-sample unbiased (Σp
358.1 vs 359; 284.3 vs 287).

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

**Overdispersion φ is weakly estimated and must be read as assumption-driven, not
data-validated, in the decision band.** φ is estimated per band by method of moments. In
the low/action bands it rests on many guides (CFD φ ≈ 1.83 / 1.92 / 0.99) and is
reasonably supported. In the **decision band it rests on only 1–3 guides carrying positives
(§8.2)**: a method-of-moments overdispersion estimate from 1–3 clusters has essentially no
degrees of freedom and enormous sampling error. The reported decision-band values
(**φ ≈ 0.19 CFD, 0.01 CB**) are therefore **not meaningfully estimated**, and the resulting
decision-band interval width is **assumption-driven, not empirically estimated**. The
decision-band coverage figures below must not be presented as data-validated. *The report
should state, per band, the number of guides contributing any positive; absent a pooled /
hierarchical / shrinkage estimate of φ with a stated prior, the decision-band interval
width is a modelling assumption.* (Direction of the resulting values: φ is large in the
low/action bands (interval honestly wide where a single guide's rate is least
predictable) and small at the top.)

**Validation: per-held-out-guide (single novel guide) coverage.** We validate on the
deployment case directly: leave **one** guide out, build *its* per-band 90% prediction
interval from the other 45 (CB: 39) guides, and check whether that held-out guide's **own**
realized true-positive count falls inside. This is the single-deployed-guide estimand, not
the aggregate total-count, which would average out exactly the between-guide variation §6
identifies as dominant. Averaged over the 46 (CB: 40) folds, realized coverage of the
nominal 90% interval is:

| coverage of nominal 90% | CFD 0.2–0.5 | CFD ≥0.5 | **CFD ≥0.2** | CB ≥0.2 |
|---|---|---|---|---|
| naive Poisson-binomial (1.3–6.3× too narrow) | ~72% | 83% | n/a | n/a |
| two-level, no overdispersion | 76% | 87% | n/a | n/a |
| **two-level + overdispersion (shipped)** | **95%** | **93%** | **95%** | **100%** |

**All shipped bands reach ≥90% on this single-guide estimand**, so the interval is
coverage-valid for a single deployed guide, not merely for the aggregate total. How to read
these numbers:
- They are empirical coverages from only **46 (CB: 40) folds**, binomial SE ≈4.4pp (n=46) /
  ≈4.7pp (n=40) at p = 0.90, so **"95% vs 93%" are not statistically distinguishable**, and
  the three-of-four over-90% values indicate a mildly **conservative (slightly too wide)**
  interval, acceptable for an IND but labelled conservative rather than "nominal." **"100%"
  at n = 40** implies an exact one-sided 90% Clopper–Pearson lower bound of **94.4%**
  (95.1% for 46/46), not literal full coverage.
- Overdispersion φ is well-supported for the **≥0.2 aggregate** (0.806, estimated on all 46
  guides) but the **decision-band** φ (≈0.19 CFD, ≈0.01 CB) rests on only **1–3 guides
  (§8.2)** (a method-of-moments estimate with essentially no degrees of freedom), so the
  decision-band interval width is an **assumption, not a data estimate**, and the
  contributing-guide count is stated on the decision-band rows.

The naive Poisson-binomial interval (treating sites as independent, 1.3–6.3× too narrow)
under-covers to ~72% and is retired.

**Dominant-guide (22%) sensitivity, reported as a full result.** The single guide
contributing 79/359 = 22% of positives (§2) drives the calibration. Removing it
(leave-that-guide-out) shifts the CFD action-band interval ~20%; the full recomputed
point estimate and interval for the worked example under this guide's removal are given in
§7. The number of benchmark guides carrying **any** decision-band positive (the count on
which both the decision-band calibration and φ hinge) is **1–3 (§8.2)** and must be stated
on the face of the decision-band rows.

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
rows below therefore cover the full ≤1-gap off-target set (5,671 at CFD ≥0.2 = 217 at ≥0.5
+ 5,454 at 0.2–0.5), **reference and variant-allele together**. What they do *not* resolve
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
  overdispersion (φ ≈ 0.806), not by combining the sub-bands.** Naively summing the
  sub-band draws *as if independent* would wrongly raise the lower bound to ~5, because a
  guide's promiscuity is correlated across bands (empirical cross-band correlation ≈ +0.66);
  modelling that correlation (pooled band) puts the floor at **2**, where it belongs.
- **CRISPR-Bulge ≥0.2 | 6 | ≈2 | [0, 6]: a low CB count is NOT clearance.** The lower
  bound of 0 on 6 sites must not be read as "nothing here." Per §1–§3, CRISPR-Bulge
  confirms only ~half of assay-detected off-targets and is defined only on ≤1-gap
  alignments, so **real off-targets can sit entirely outside what CRISPR-Bulge counts.** A
  low CB number means CB did not prioritize these sites, not that the locus is clear.
- **Dominant-guide (22%) sensitivity.** Dropping the single 79-positive benchmark guide
  (§2, §5) and refitting on the other 45 moves the CFD ≥0.2 row from ≈28 [2, 81] to
  **≈23** (~−18% point; φ_≥0.2 0.806 → 0.62); the headline does not hinge on that one
  cluster.

The deployment prior-shift axis (§6) applies on top of the CFD ≥0.2 estimate: this guide's
own off-target propensity at its locus is unknown and can move the count substantially in
either direction — usually a larger uncertainty than the statistical range above. The
summary states this qualitatively and attaches no specific multiplier. The two scorer
estimates are different views of the same sites and are **not additive**.

## 8. Limitations (carried into the report)

1. **Information ceiling**: 359 positives in 46 guides (one holding 22%). No calibrator
   manufactures information the scores lack; CFD's AUPR 0.168 (baseline prevalence
   1.22×10⁻⁴) makes its bands inherently wide. The honest interval width is the
   deliverable, not a defect.
2. **Top-band resolution**: the decision band rests on **1–3 guides carrying positives**;
   its point is accurate but its interval is wide by construction, and its φ (hence its
   width) is assumption-driven, not data-estimated (§5).
3. **Generalization**: single-guide deployment coverage is validated per-held-out-guide
   (§5; ≥90% in all shipped bands); a chemically/structurally novel clinical guide at a
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
  low-CB-is-not-clearance notes, the priority-order guidance, and the 2-gap / uncalibrated-
  support note.


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

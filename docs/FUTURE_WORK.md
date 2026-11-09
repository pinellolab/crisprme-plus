# CRISPRme+ — Future work / pickup-ready ideas

Deferred, non-blocking enhancements. None are correctness bugs; each is a self-contained
improvement another contributor (or agent) can pick up. For each: **why**, **current behavior**,
**what to build**, **code pointers**, **acceptance**. Keep this list current — move an item to the
CHANGELOG when it ships, and delete it from here.

---

## 1. Mega report: per-database (+ superpopulation) AF breakdown
**Why.** The sites-only `mega` index's whole value proposition is *per-dataset allele-frequency
provenance* (`AF_1000G2021`, `AF_HGDP`, `AF_gnomAD`, `AF_TOPMed`, `AF_AoU`, `AF_max`), and the Tier-0
registry already stores it — but the SNP report surfaces only a single merged global MAF, so a user
can't see which database a variant's frequency came from.

**Current behavior (mega search report).** One global merged **MAF** column, a **REF/ALT origin**
flag, and a *"By dataset (provenance)"* bar chart that shows only per-source off-target **counts**
(not frequencies). Per-dataset AF is already exposed for **indels** via the `*.indel_af.tsv`
companion, but **not for SNPs** and not in the main table/panel.

**What to build.** Surface the per-database AFs (and, where available, a per-db × superpop analog)
for SNP off-targets in the report — e.g. `AF_<db>` columns in the curated table and/or a
per-database frequency panel analogous to the genotyped superpop panel. Data already lives in the
registry db-groups; no new computation at search time.

**Code pointers.** `PostProcess/generate_report.py`: `CURATED_COLUMNS` (~L196–222), `_curated_cell`
(MAF/origin cells ~L663–679), `plot_population` (~L2259–2329, the `use_superpop` branch that already
switches genotyped→superpop vs mega→"by dataset"), the indel per-db companion (~L4255–4291) as the
precedent; registry db-group / `db_subpop_group_id` in the Tier-0 code. (Was tracking task #244.)

**Acceptance.** A mega search report shows per-dataset AF per SNP off-target (matching the registry
values, and matching the indel companion's shape); genotyped indices unchanged.

---

## 2. Dense-panel search: what is left of the post-analysis' parallelism
**Why.** Dense variant panels are the slowest searches. The combinatorial blow-up itself is tamed
(the `--max-edits` bound, the high-variant-density greedy cap `CRISPRME_IUPAC_CAP`, the Tier-0
registry cache, and since 2.10.0 the exact haplotype pre-check), and since 2.10.0 the post-analysis
runs each chromosome in pieces over one pool of `--thread` workers (`PostProcess/postproc_pool.py`:
threads that each run one subprocess at a time, no longer the `multiprocessing` `Pool` whose
`pool.map` deadlocked at `--thread 200` in earlier releases; the files are those of one process a
chromosome, byte for byte: METHODS.md §5). What is left is serial within a chromosome.

**Current behavior.** Per chromosome, the extraction, sorting and clustering of the raw targets
(about 5 min for IGH's 37.6 million indel raw targets), the merge of the pieces, and the
CRISPR-Bulge scoring of all of the chromosome's rows, which the merge does in one process to keep
the batches of one process (13 of the 20 min of a population-level SNP stage of chr22 on a loaded
node), run on one core. Equal line ranges are not equal costs (IGH: one SNP piece 294 s, the
others 77 s or less). Each piece repeats the start-up one process does (about 30% more CPU per
sample on chr22).

**What to build.** Score the CRISPR-Bulge batches in parallel, each in the batch one process forms,
with a test that the floats do not depend on the process that scores them; cut the ranges by an
estimate of their cost (the windows' variant columns) rather than by lines; measure a genome-wide
search on the indexes of the 2.10.0 rebuild at `--thread 16` and `--thread 200`.

**Code pointers.** `PostProcess/postproc_pool.py` (`run_stage`, `n_pieces`, `ranges`),
`PostProcess/merge_pieces.py` (`score_cbulge`), `PostProcess/cbulge_rows.py`; the drivers'
`cluster` and `finish` stages (`scriptAnalisiNNN_v3.sh`, `analisi_indels_NNN.sh`).

**Acceptance.** The same files, byte for byte, as `CRISPRME_POSTPROC_PIECES=0` (RELEASE_CHECKLIST
step 4d), and a measured speed-up on IGH alone and genome-wide.

---

## 3. Two-pass lossless dense-region search (exploration)
**Why.** A principled alternative to the greedy cap for ultra-dense regions: lossless worst-case
scoring without full 2^k-IUPAC haplotype enumeration.

**What to build (design, not yet implemented).** Pass 1 = score-free **minimum-edit** find; Pass 2 =
per-score **worst-possible** reconstruction. The min-edit lower bound is proven for
SNP + bulge + single-indel, but **BREAKS** (verified false-negative even at k=0) for ≥2 cis indels in
a protospacer and for equal-length MNVs — a size-aware slack / flag-all-multi-indel-windows / `norm -a`
fix is needed; multi-indel-in-cis is common in raw counts (~17–19% of real indel off-targets), most
of it in STR/VNTR repeats: ~1–2% of indel loci after repeat-masking, and ~0.1–0.2% of indel
off-targets genuinely missed (METHODS §8). CFD worst-case is
already exact via factorization (`twopass_cfd_exact.py`, 0-ULP verified); CRISPR-Bulge would be
shortlist-only.

**Code pointers.** `PostProcess/twopass_cfd_exact.py`, the fixture at commit `113b9b7`; memory note
`two-pass-dense-search-design` (the lynchpin analysis + caveats). (Was tracking task #175.)

**Acceptance.** A design that provably loses no off-target vs. full enumeration on the hard cases
(≥2-cis-indel, MNV), benchmarked against the greedy cap.

---

## 4. HPRC high-edit-budget validation search
**Why.** The guide `ACTGAAATCTGTAAGCAGGC` has been searched at high edit budget on 1000G+HGDP and on
the `mega` cap-8 worst-case (6.77M OTs), and at 4/1/1 on HGDP1kGP vs mega — but not on the **HPRC**
pangenome index at a high budget. A quick coverage item, not a gap.

**What to build.** Run `ACTGAAATCTGTAAGCAGGC` on `NRG_3_hg38+hg38_HPRC` at a high edit budget
(e.g. 6/2/2 with `--max-total-edits`) and record the result for the validation matrix. (Was tracking
task #242.)

**Acceptance.** HPRC ACTGA high-budget result captured alongside the other indices'.

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

## 2. Dense-panel search: per-chromosome parallelism + benchmark
**Why.** Dense / pamless variant panels are the slowest searches. The combinatorial blow-up itself is
already tamed (the `--max-edits` bound + the high-variant-density greedy cap `CRISPRME_IUPAC_CAP` +
the Tier-0 registry cache / O(1) LRU that shipped in 2.5.3 and remains in place as of 2.7.10). What's left is raw throughput: the
search still processes contigs largely serially.

**What to build.** Run the dense-panel search **per chromosome in parallel** and **measure** the
wall-clock speedup vs. the current path on a dense panel (e.g. the mega or a pamless index). Fold in
the known post-analysis fix while here: `pool_post_analisi_snp.py` spins up `Pool(cores)` (~200 on a
big host) across hundreds of contigs and FD-inheritance **deadlocks** the parent `pool.map` at very
high `--thread`; cap workers (~16–32) + `maxtasksperchild`.

**Code pointers.** `PostProcess/pool_post_analisi_snp.py` (the `Pool` worker count), the
`complete-search` submit chain in `crisprme.py`, `CRISPRME_IUPAC_CAP` + the high-variant-density
region handling. See memory notes `crisprme-postanalysis-pool-deadlock`,
`variant-search-complexity-explosion`. (Was tracking task #174.)

**Acceptance.** Measured speedup on a dense/pamless genome-wide search with byte-identical results,
and no `Pool` deadlock at `--thread 200`.

---

## 3. Two-pass lossless dense-region search (exploration)
**Why.** A principled alternative to the greedy cap for ultra-dense regions: lossless worst-case
scoring without full 2^k-IUPAC haplotype enumeration.

**What to build (design, not yet implemented).** Pass 1 = score-free **minimum-edit** find; Pass 2 =
per-score **worst-possible** reconstruction. The min-edit lower bound is proven for
SNP + bulge + single-indel, but **BREAKS** (verified false-negative even at k=0) for ≥2 cis indels in
a protospacer and for equal-length MNVs — a size-aware slack / flag-all-multi-indel-windows / `norm -a`
fix is needed; multi-indel-in-cis is common (~17–19% of real indel off-targets). CFD worst-case is
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

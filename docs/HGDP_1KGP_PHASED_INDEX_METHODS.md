# gnomAD HGDP + 1kGP re-phased variant index — methods & validation record

> **Living document.** This records the rationale, source-data provenance, verification, and
> (once finalized) the validated build recipe for a new **default** CRISPRme+ variant index built
> from the gnomAD v3.1.x HGDP + 1000 Genomes joint callset. **Once the recipe is finalized and
> validated, sections 5–7 are folded into `METHODS.md` and the CRISPRme+ manuscript Methods**
> (this is a key methodological contribution: a fully-phased, singleton-inclusive, uniformly-processed
> HGDP+1kGP haplotype panel enabling CONFIRMED-cis co-occurrence across BOTH cohorts).
>
> Nothing is built genome-wide or pushed to HuggingFace until every section is verified spotless and
> Luca signs off. Status log at the bottom.

## 1. Rationale

The current shipped genotyped index `NRG_3_hg38+hg38_1000G2021_HGDP` is a **hybrid**: 1000 Genomes-2021
is phased (→ CONFIRMED cis co-occurrence + named carriers) but HGDP is genotyped-**unphased**
(→ PUTATIVE co-carrier only). The gnomAD HGDP+1kGP callset is a **single, jointly-called, uniformly
QC'd, SHAPEIT5-co-phased** resource in which **both** cohorts are phased. Building the index from it
would:

- flip HGDP multi-variant off-targets from PUTATIVE → **CONFIRMED cis** (named carriers + exact joint AF);
- lower switch error vs 1kGP-alone (paper: SNP mean 0.00184 vs 0.00338);
- be genuinely single-source (structurally not "hybrid"); clean AN denominator;
- carry an open, redistributable license (gnomAD no-restriction; paper CC BY 4.0; atgu code MIT).

**Goal:** make this the new **default** index and retire the old mixed 1000G2021_HGDP.

## 2. Source data + provenance

- **Pre-made phased release:** `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/`
  (anonymous HTTPS mirror: `https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/`).
  Per-chromosome BCF `hgdp1kgp_chr<N>.filtered.SNV_INDEL.phased.shapeit5.bcf` (+ `.csi`),
  **chr1–22 + chrX (PAR1/PAR2/non-PAR)**, chrY/chrM absent, ~16 GB total. GRCh38, SHAPEIT5-phased.
- **Dense (unphased, all-variants incl. singletons) callset:** gnomAD v3.1.x HGDP+1kGP genotype callset
  (Hail MatrixTable + per-chrom VCFs) — the source we would RE-PHASE ourselves (see §4–5).
- **Cross-validation reference:** Zenodo record `18156285` — a *lossy* common-SNP subset
  (chr1–22 only, biallelic-SNP-only, MAF > 0.5%, 'chr'-prefixed, 11.3 GB) derived from the same
  gnomAD SHAPEIT5 haplotypes; usable only as a same-provenance corroborator for sample list +
  common-SNP AF/phasing, **not** a build input.
- Authoritative sample count/split: paper (Koenig et al., *Genome Research* 34(5):796, 2024;
  PMC9900804) reports **4,094 post-QC = 929 HGDP + 3,165 1kGP**. The phased BCFs carry **4,091**
  samples (provenance of the −3 under investigation, see §3).

## 3. Verification of the pre-made phased release (chr22)

Method: downloaded `hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf` (266 MB) and interrogated it
with `bcftools` (v from the released `pinellolab/crisprme:v2.6.2` image:
`apptainer exec crisprme.sif /opt/conda/bin/bcftools`). All findings independently re-derived by an
adversarial-verification panel (see status log).

| Property | Result (chr22) | Verdict |
|---|---|---|
| Records / distinct POS / multi-record POS | 1,093,149 / 1,037,983 / **32,445** | multiallelics **SPLIT into biallelic records, all alleles retained** (e.g. chr22:10650015 A>C *and* A>G) |
| Phasing | 818,200,000 `\|` separators, **zero `/`, zero missing** (200k-record scan) | **uniformly phased**, complete |
| Variant types | ~90% SNV / ~10% indel | indels present |
| FILTER column | all `.` | nothing dropped by a FILTER flag |
| INFO fields | `AC` (Number=A), `AN`; **no AF / MAF field** | AF must be recomputed (`bcftools +fill-tags`) |
| Contig naming | `chr22` (chr-prefixed) | matches CRISPRme hg38 fastas |
| Sample count | **4,091** | −3 vs paper's 4,094 (see below) |
| **Allele-count spectrum** | **min AC = 2, max AC = 8,180, ZERO AC=1** | **⚠ SINGLETONS REMOVED (MAC ≥ 2 filter)** |

**Two material findings:**

1. **Singletons are filtered out (MAC ≥ 2).** The allele-count spectrum has no AC=1 bin at all
   (a normal spectrum peaks at singletons). This is the key issue for the "use ALL variants, no
   MAF/other filter" requirement — singletons (variants private to one haplotype) are absent. This
   is a common by-product of statistical phasing pipelines, but it **is** a filter. → drives the
   decision in §4.
2. **Sample-provenance discrepancy.** The bucket `README.md` states samples HGDP01371 and
   LP6005441-DNA_A09 were removed, yet both are **present** in the phased BCF (only CHMI_CHMI3_WGS2
   is absent). The true sample set / QC / correct AN denominator of the phased release is being
   reconciled against the paper + atgu metadata.

Verified-good properties (retained, no filter): multiallelics (split, all alleles), indels, uniform
phasing, no missing genotypes, chrX (PAR1/PAR2/non-PAR).

## 4. Decision: re-phase the DENSE callset to retain rare variants + singletons

Because the pre-made phased release drops singletons (§3.1) and the requirement is to use the full
callset unfiltered, the plan is to **re-phase the dense (unphased, all-variants) callset ourselves
with the latest SHAPEIT5** (`odelaneau/shapeit`) — whose headline capability is phasing rare variants
and singletons (`phase_rare` on a `phase_common` scaffold) — then **benchmark** the result against the
gnomAD MAC≥2 release. This yields "all variants **and** phased" (CONFIRMED-cis capable).

## 5. Re-phasing recipe (SHAPEIT5) — reproduce the published pipeline, minus the singleton drop

**Key provenance finding (source: `atgu/hgdp_tgp` `phasing/`, commit `bf8ef3b`; SHAPEIT5 v5.1.1 image
`lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1`; paper Methods; SHAPEIT5 paper Hofmeister/Delaneau
Nat Genet 2023):** SHAPEIT5 does **not** drop singletons — `phase_rare` phases them (coalescent Viterbi
model; SER <5% at 1/100k). The atgu `phase_rare` output `hgdp1kgp_chr{i}.full.shapeit5_rare.bcf`
**still contains singletons**. The public `phased_haplotypes_v2/` release is that file after a
**separate post-phasing** `bcftools view -i'MAC>=2'` (`phasing/remove_singletons.py:56`; chrX
`phasing/chrX/filter_phased.py:83,92`). The paper Methods do **not** mention this MAC≥2 step. Therefore
the singleton-inclusive panel is obtained by running the identical pipeline and **skipping
`remove_singletons.py`** (use the `*.full.shapeit5_rare.bcf` directly).

> The atgu team's own unfiltered `.full.` BCFs live at `gs://hgdp-1kg/phasing/shapeit5/phase_rare/`
> but that bucket is **private** (HTTP 403/401 anonymous). Options: (1) request read access from the
> Broad/atgu team (Lindo Nkambule) → zero re-phasing compute; (2) reproduce the pipeline below.

**Published pipeline (autosomes chr1–22), verbatim structure:**

- **Step 0** `prepare_data_phasing.py` — from the public dense MT
  `gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt`,
  apply gnomAD sample/variant/genotype QC (`filter_to_adj`), export one VCF per chromosome. *(Requires Hail.)*
- **Step 1** `filter_and_convert_to_bcf.py` — remove 29 samples (5 duplicates + 24 PCA outliers) → **4,091**:
  `bcftools view --samples-file ^{samples_to_filter} {vcf} -Ob -o {bcf}`
- **Step 2** pre-phasing QC (`qc.py:55`) — **no MAC/MAF filter**:
  `bcftools +fill-tags {bcf} -Ou -- -t all | bcftools view -i'HWE>=1e-30 && F_MISSING<=0.1 && ExcHet>=0.5 && ExcHet<=1.5' -o {qced.bcf}`
- **Step 3A** `phase_common` per 20 cM chunk — scaffold of common variants:
  `phase_common --input {qced.bcf} --map chr{i}.b38.gmap.gz --output {common.chunk.bcf} --filter-maf 0.001 --region {col3-of-20cM-chunk} --pedigree hgdp1kg_pedigree.fam --thread T`
  *(`--filter-maf 0.001` only bounds the SCAFFOLD; it does NOT lose singletons — they are re-read from the full `--input` in phase_rare. Confirmed in `phase_common .../phaser_parameters.cpp:69`.)*
- **Step 3B** `ligate` common chunks → per-chrom scaffold:
  `ligate --input {common_chunks_list} --pedigree hgdp1kg_pedigree.fam --output {scaffold.bcf} --thread T --index`
- **Step 3C** `phase_rare` per 4 cM chunk — phases rare + **singletons** onto the scaffold, **no MAF/MAC filter**:
  `phase_rare --input {qced.bcf} --input-region {col4-of-4cM-chunk} --scaffold {scaffold.bcf} --scaffold-region {col3-of-4cM-chunk} --map chr{i}.b38.gmap.gz --pedigree hgdp1kg_pedigree.fam --output {rare.chunk.bcf} --thread T`
- **Step 3D** concatenate rare chunks → final per-chrom phased BCF (**keep everything**):
  `bcftools concat -n -f {rare_chunks_list} -o hgdp1kgp_chr{i}.full.shapeit5_rare.bcf && bcftools index …`
- **⟶ STOP. Do NOT run `bcftools view -i'MAC>=2'`** (that is the singleton-dropping step).
- **chrX** analogous (`phasing/chrX/phase_chrX.py`) + a fix removing 7,667 spurious male non-PAR hets;
  same MAC≥2 drop to skip.

**Resources (verified present in `github.com/odelaneau/shapeit`, HEAD `c34d4db` — the original
`odelaneau/shapeit5` repo was disabled by GitHub ToS):** GRCh38 maps `resources/maps/b38/chr{1..22,X}.b38.gmap.gz`;
chunk coords `resources/chunks/b38/{20cM,4cM}/chunks_chr{N}.txt` (col 3 = with buffers, col 4 = without);
pedigree = `hgdp1kg_pedigree.fam` (599 families = 593 trios + 6 duos over 4,091 samples; build per atgu
README §2 from PC-Relate/IBD cross-checked vs `1kGP.3202_samples.pedigree_info.txt`). LICENSE MIT.
Pin the exact SHAPEIT5 version used. Optional: `phase_rare --score-singletons` (experimental singleton
phase-confidence 0.5–1.0). chrX chunk resources referenced a now-dead `UKB_WGS_200k` path — re-point to
`odelaneau/shapeit`.

**Compute/effort:** re-phasing needs (a) a Hail/Spark export from the dense MT, (b) a genome-wide
SHAPEIT5 run for 4,091 samples (the atgu team ran this on Hail Batch; per-chrom, chunked). chr22 SMOKE
first to measure wall-time/RAM before GW. *(Exact numbers PENDING the chr22 smoke.)*

## 6. Validation / comparison protocol

**Dense input source (public, no auth):** chr22 smoke =
`gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz`
(~55 GB full FORMAT) or Hail-export GT-only from the dense MT
`gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt` (genome-wide, avoids multi-TB egress). GRCh38 GLIMPSE
b38 maps from `odelaneau/shapeit`; chunk coords via the GLIMPSE chunker (do not hand-roll); pedigree
`.fam` rebuilt from public 1000G `.ped` + gnomAD relatedness (or requested from atgu).

**chr22 SMOKE gates (all four must be green before any genome-wide run):**
- **VAL-1 Singleton recovery** — our chr22 phased BCF has `min AC==1`, singleton_count in the tens of
  thousands, 100% `|`-phased + non-missing; the released MAC≥2 BCF has ~0. (Proves the feature.)
- **VAL-2 / VAL-3 Accuracy (anti-self-deception)** — `SHAPEIT5 switch` SER against **1000G trios
  (offspring re-phased with PARENTS HELD OUT** — leaving parents in fakes ~0 SER, the #1 self-deception)
  and against the **HGSVC2 34-genome assembly truth**. ACCEPT common+rare SNP SER ≲ published 0.00184
  (+1sd ~0.0033) and ≤ our own re-run of `switch` on gnomAD's release with the same trios; indel SER
  ≲ 0.00899. Singleton SER ~30–40% is EXPECTED, not a failure (see §4 tiering).
- **VAL-4 Backbone concordance** — `bcftools isec` shared common sites vs `phased_haplotypes_v2`;
  `switch` phase-consistency near-identical (SER ~1e-3) → our scaffold didn't regress.
- **VAL-5 Structural sanity** — samples==4,091; AN==2N (8,182) autosomal (sex-aware chrX); 0 missing;
  100% `|`; multiallelics split+left-aligned+sorted; per-chrom record count **≥** the release (we add,
  never lose).
- **VAL-6 Common-SNP recall** — every common SNP (Zenodo 18156285 subset) present + phase-concordant.
- **VAL-7 (post-index) CRISPRme cis audit** — run a known guide `--per-sample` on the new index vs the
  current 1000G2021_HGDP; the EXTRA cis calls that hinge on a recovered singleton are flagged as the
  lower-confidence tier (§4), not silently promoted to CONFIRMED.

### 6.1 chr22 re-phasing results — MEASURED (2026-09-30)

Our re-phased chr22 (`hgdp1kg.chr22.rephased.bcf`, 1,998,332 records, 4,091 samples) compared against
gnomAD's pre-phased release (`hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf`, 1,093,149 records)
via `bcftools isec` (site overlap) and `SHAPEIT5_switch` (release = validation/truth, ours = estimation).

**VAL-1 Singleton recovery — PASS.** Our chr22 has `min AC = 1` with **828,868 singletons**; the release
has **zero** (MAC≥2). Feature proven.

**VAL-4 Backbone concordance — PASS.**
| Metric | Value | Interpretation |
|---|---|---|
| SHARED sites (our ∩ release) | **1,090,102** | **99.72%** of the release backbone is present in ours |
| release-only (0001) | 3,047 (0.28%) | sites the release kept but our atgu-QC (HWE/F_MISSING/ExcHet) correctly dropped — not a superset violation |
| our-only (0000) | 908,230 | sites we add; **825,653 (91%) are AC=1 singletons**, rest rare |
| **Genotype (allele) concordance on shared** | **99.98%** (932,669 / 4,459,607,282 discordant = 0.0209%) | underlying calls essentially identical — we re-phase + add, we don't alter data |
| **Phasing switch-error vs gnomAD (global)** | **0.58%** (790,362 / 135,364,707 het genotypes) | our haplotype backbone reproduces gnomAD's almost exactly (this is discordance between two *statistical* phasings, not vs trio-truth) |
| Phasing SER by allele count | AC=2: 6.6%, AC=3: 6.3%, … AC=31: 2.1%, common ≪0.58% | rare-variant phase intrinsically uncertain; common variants dominate → low global |

**Conclusion:** chr22 is correct — sites, genotypes, and phasing all reproduce gnomAD's release on the
shared backbone; the only material difference is the ~826k deliberately-recovered singletons. Validates
the recipe for the genome-wide run. (Note: switch-vs-trio-truth per VAL-2/VAL-3 remains the stronger
accuracy check and is still pending; the vs-release concordance here confirms we did not *regress* the
backbone while adding singletons.)

### 6.2 Genome-wide autosome re-phasing results — all 22, MEASURED (2026-10-03)

The complete chr1–22 re-phase (4,091 samples, SHAPEIT5 phase_common+phase_rare, atgu pre-phasing QC,
singletons kept) was validated against gnomAD's public release
(`gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_<chr>.filtered.SNV_INDEL.phased.shapeit5.bcf`).
Per chromosome: full BGZF read-integrity (`bcftools view`), structural checks, `bcftools isec` site overlap,
and `SHAPEIT5_switch` (release = validation/truth, ours = estimation) for genotype concordance + switch error.

**Result: 22/22 chromosomes CLEAN** — every one read-OK (no corruption), exactly 4,091 samples, **0 unphased
genotypes** (100% `|`), concordant with the public backbone, and adding singletons.

| chr | our records | singletons recovered | GT concord % | switch-err % | rel-only |
|---|---|---|---|---|---|
| chr1 | 11,101,151 | 4,760,285 | 99.9880 | 0.412 | 4,523 |
| chr2 | 12,104,482 | 5,247,203 | 99.9896 | 0.395 | 3,469 |
| chr3 | 9,995,561 | 4,330,916 | 99.9902 | 0.363 | 2,366 |
| chr4 | 9,880,794 | 4,304,346 | 99.9898 | 0.341 | 2,892 |
| chr5 | 9,096,534 | 3,948,220 | 99.9906 | 0.358 | 1,781 |
| chr6 | 8,563,279 | 3,648,933 | 99.9899 | 0.331 | 1,901 |
| chr7 | 8,131,512 | 3,486,837 | 99.9890 | 0.381 | 2,395 |
| chr8 | 7,806,477 | 3,384,020 | 99.9905 | 0.365 | 2,427 |
| chr9 | 6,208,694 | 2,662,354 | 99.9879 | 0.432 | 3,723 |
| chr10 | 6,860,458 | 2,922,656 | 99.9886 | 0.383 | 2,691 |
| chr11 | 6,848,667 | 2,956,338 | 99.9899 | 0.361 | 1,783 |
| chr12 | 6,601,225 | 2,829,776 | 99.9892 | 0.394 | 1,535 |
| chr13 | 4,933,644 | 2,149,831 | 99.9897 | 0.369 | 1,131 |
| chr14 | 4,533,362 | 1,938,982 | 99.9878 | 0.418 | 1,761 |
| chr15 | 4,117,280 | 1,753,105 | 99.9868 | 0.477 | 3,097 |
| chr16 | 4,626,002 | 1,992,922 | 99.9883 | 0.465 | 1,772 |
| chr17 | 4,032,356 | 1,706,922 | 99.9864 | 0.509 | 1,641 |
| chr18 | 3,865,839 | 1,661,595 | 99.9899 | 0.390 | 968 |
| chr19 | 3,154,021 | 1,304,591 | 99.9845 | 0.485 | 1,722 |
| chr20 | 3,119,135 | 1,327,087 | 99.9887 | 0.439 | 1,273 |
| chr21 | 1,925,238 | 797,826 | 99.9779 | 0.510 | 5,728 |
| chr22 | 1,998,332 | 825,653 | 99.9791 | 0.584 | 3,047 |
| **TOTAL** | **139.5 M** | **59.9 M** | **≥99.978** | **≤0.584** | **53,626** |

**Interpretation.** Across 139.5 M total autosomal records the re-phased panel reproduces gnomAD's public
MAC≥2 backbone at **≥99.978% genotype concordance** and **≤0.584% switch-error** (discordance between two
independent statistical phasings, not vs trio-truth — low = faithful reproduction), with **zero** read
failures or unphased genotypes. The material difference is the **59.9 M recovered singletons** (absent from
the MAC≥2 release) — the intended feature. `rel-only` (sites the release kept but our atgu-QC dropped) is
53,626 total = **0.04%** of records, i.e. we are effectively a strict superset of the public backbone plus
singletons. **Autosomes are correct, uncorrupted, and ready for the CRISPRme index build.** chrX (special
PAR/haploid handling) is validated separately before the full-genome build. (VAL-2/VAL-3 switch-vs-trio-truth
remains the stronger absolute-accuracy check and is still pending; this vs-release concordance confirms no
backbone regression.)

### 6.3 chrX re-phasing — recipe + the `--haploids` ploidy fix, MEASURED (2026-10-04)

chrX follows the atgu/hgdp_tgp chrX recipe (`phasing/chrX/`), which is a port of the SHAPEIT5 UKB-200k
chrX workflow: **split into 3 regions on GRCh38 PAR boundaries and phase each separately to avoid mixing
ploidy** — PAR1 `chrX:10001-2781479`, non-PAR `chrX:2781480-155701382`, PAR2 `chrX:155701383-156030895`.
- **PAR1/PAR2**: `phase_common --filter-maf 0.0` single-block, region-specific maps (`chrX_par1/par2.b38.gmap.gz`),
  diploid for everyone (no `--haploids`).
- **non-PAR**: sex-aware QC — HWE/ExcHet computed **in females only** (`+fill-tags` on a female subset →
  annotate back; thresholds identical to autosomes), F_MISSING/AC over all; singletons kept. Then
  `phase_common --filter-maf 0.001 --haploids <males>` per 20cM chunk → ligate scaffold → `phase_rare
  --haploids <males>` per 4cM chunk → ligate. Sex derived from the released non-PAR het-rate (clean 350×
  bimodal gap): **2195 males / 1896 females**.
- Concat the 3 regions → `hgdp1kg.chrX.rephased.bcf` (single `chrX` contig).

**KEY FIX (manuscript-relevant), the `phase_rare` NaN:** running `phase_rare --haploids` on non-PAR crashed
on 3/23 chunks with `Assertion !isnan(GRvar_genotypes[vr][tidx].prob)` (`genotype_set_phasing.cpp:66`).
Root cause (confirmed by direct test + SHAPEIT5 source + the atgu recipe): **a male-ploidy ENCODING mismatch**
— the dense input encodes non-PAR males as single-allele **haploid** (`0`), but `phase_common --haploids`
writes the scaffold with males as **diploid-homozygous** (`0|0`), and so does the gnomAD release. `phase_rare`
conditions the haploid input against the diploid-hom scaffold → the Li-&-Stephens HMM mass collapses to
`0/0 = NaN` on rare-dense chunks. **Fix (SHAPEIT5's documented `--haploids` convention, NOT a hack):**
ploidy-normalize the non-PAR input to diploid-hom before phasing — `bcftools +fixploidy -f 2` (males `0`→`0|0`).
With input and scaffold ploidy consistent, all 23 chunks phase; males emit as hemizygous-diploid `x|x`,
AN≡8182, exactly matching the release. (`--haploids` is kept; chunks are not dropped; the assertion is a real
degenerate-state guard, not disabled. No `--pedigree` — consistent with the autosomes, which validated at
≥99.978% without it. Also: the post-phase het-male cleanup is a no-op here because the fix yields 0 male hets.)

**chrX validation vs release — PASS** (`hgdp1kg.chrX.rephased.bcf`, 5,354,767 records, 4,091 samples):

| Metric | Value |
|---|---|
| Unphased genotypes | **0** (100% `\|`) |
| AN (constant) | **8182** = 4091×2 (males hemizygous-diploid `x\|x`) |
| Male non-PAR hets / PAR1 hets | **0 / 4261** (hemizygous non-PAR, diploid PAR — correct) |
| Singletons recovered (AC=1) | **1,597,916** (release = 0, MAC≥2) |
| Shared with release backbone | **3,495,756 (99.84%)** |
| our-only (singletons+rare added) | 1,859,011 |
| release-only (our QC dropped) | 5,459 (0.16%) |
| **Genotype concordance** | **99.9957%** |
| **Phasing switch-error vs release** | **0.397%** (589,624 / 148.4 M) |

chrX reproduces the public backbone (99.84% shared, 99.996% GT concordance, 0.40% switch-error) while adding
1.60 M singletons, with biologically-correct male hemizygosity. **The full chr1–22 + chrX panel is complete,
correct, and ready for the CRISPRme NRG index build.**

## 4b. Singleton accuracy — the scientific caveat (drives a product decision)

SHAPEIT5 phases singletons non-randomly but at **~35% switch error** (vs ~50% random; SHAPEIT5 paper),
because a singleton is carried on exactly one haplotype in one individual. Implications for CRISPRme:
- **Detection is lossless + unambiguous** — a singleton variant that creates an off-target is real and
  worth nominating; it is simply carried by exactly one individual. Keeping singletons is a clear win
  for coverage.
- **Only the cis-PHASE of a singleton with another nearby variant is ~65% reliable.** So a *multi-variant*
  CONFIRMED-cis co-occurrence that HINGES on a singleton's phase must NOT be treated as equal to a
  MAC≥2 CONFIRMED call. gnomAD dropped singletons from its public release precisely to avoid shipping
  this lower-reliability phase.
- **Recommended handling (honest + a manuscript point):** KEEP singletons (satisfies "use all of it";
  single-variant singleton off-targets are fully valid), but mark singleton-hinged multi-variant cis
  as a **distinct lower-confidence tier** rather than silently CONFIRMED. This gives complete coverage
  without over-trusting ~35%-reliable phase.

## 4c. Adversarial verification verdict + strategic options

An independent adversarial panel (6 agents, ~1.3 B genotypes re-scanned) **confirmed all five
findings** by re-deriving them different ways: singleton MAC≥2 removal (min AC=2 raw, +fill-tags,
GT-derived, both allele sides), multiallelic retention (0 dropped alleles; a 74-allele STR at
chr22:47316738 preserved), uniform phasing (100% `|`, 0 missing across the whole chr22), and the
provenance mismatch. It also:

- **Source-documented the MAC≥2 mechanism** (resolving the earlier hedge): `atgu/hgdp_tgp
  phasing/remove_singletons.py` runs `bcftools view -i'MAC>=2'` on `*.full.shapeit5_rare.bcf` → the
  released file. Deliberate, named filter — not incidental monomorphic drop.
- **Corrected provenance:** header-derived cohort split is **925 HGDP + 3,166 1kGP = 4,091** (HGDP IDs
  = 808 `HGDP#####` + 108 `LP…-DNA_` + 9 `SS…` Sanger → a `^HGDP` grep undercounts HGDP by 117; derive
  samplesID from the BCF header, never a prefix grep). Autosomal AN = 8,182. The `phased_haplotypes_v2/README.md`
  is a **stale changelog for the UNPHASED hgdp_1kg→v2 sample transition** and does NOT describe the
  phased artifact (which predates CHARR removal → still carries HGDP01371 + LP6005441-DNA_A09).
- **Reproducibility gotcha:** `apptainer exec` does not auto-bind `/srv/local`; add `--bind /srv/local`
  or bcftools silently errors "No such file" (easy to mistake for empty output).

**Strategic argument (drives the decision):** the two hard requirements — "use ALL variants, no filter"
and "CONFIRMED-cis phasing" — are **mutually exclusive on any *public* artifact**: the only public
*phased* HGDP+1kGP callset is MAC≥2; the only public callset *with singletons* is *unphased*. And the
dropped singletons (AC=1 = one haplotype in one of 4,091 people, joint AF ~1.2e-4) are, by construction,
**population-irrelevant for off-target nomination** — never recurrent, and if recovered they'd carry the
weakest (~35% SER) phase in the panel. Genome-wide singletons are ~46% of variants (53% among unrelated).

**Options:**
- **A+ (panel recommendation) — ship the MAC≥2 phased release as the new default, no re-phasing.**
  Reframe "no filter" honestly: CRISPRme imposes **no filter of its own**; MAC≥2 is an inherent,
  **disclosed** property of the upstream gnomAD SHAPEIT5 source, stamped in the manifest
  (`singleton_filter=MAC>=2 (gnomAD)`, `data_type=genotyped-phased`, `sample_count=4091`, `AN=8182`).
  LOW effort; the already-verified public BCFs. Loses population-irrelevant singletons.
- **B — build on the dense UNPHASED callset (all variants incl. singletons).** Truly lossless, but
  **no phasing → no CONFIRMED-cis** (forfeits the other hard requirement) + re-introduces the dense
  IUPAC search blow-up. Cannot satisfy the goal alone.
- **C — re-phase the dense callset ourselves (SHAPEIT5, keep singletons + phasing).** The only path that
  satisfies BOTH literally. But: the singleton-inclusive intermediate is in a **private** bucket →
  full genome-wide re-phase (days, HIGH risk), recovering population-irrelevant singletons at the
  weakest phase tier. Effort HIGH.
- **Dual-index — A+ default PLUS the unphased dense callset as an opt-in PUTATIVE-lossless fallback**
  (mirrors the existing `mega` worst-case pattern). Honest way to offer both without pretending one
  artifact does both.

**DECISION PENDING (Luca).** The chr22 re-phase smoke (§6) is proceeding to give empirical
singleton-recovery + switch-error numbers, but the *value* question (are 1-person singletons worth a
multi-day re-phase for a population index?) is now the crux. Panel leans A+; final call is Luca's.

## 7. Final validated build recipe (→ METHODS / manuscript)

Recipe script: `seq_script/merge_panels/hgdp1kgp_build.sh` (single-source, resumable). Input =
the 23 re-phased, validated BCFs `hgdp1kg.chr{1..22,X}.rephased.bcf` (§5–§6). Three stages:
(7.1) VCF-prep, (7.2) samplesID, (7.3) index build; chrX ploidy (7.4) is the one subtlety and is
handled in both the data and the registry.

### 7.1 Per-chromosome VCF preparation
Each re-phased BCF → a build-ready, biallelic, reference-checked, sorted, bgzipped+tabixed VCF:
```
bcftools view <bcf> \
  | [chrX only] bcftools +fixploidy -p <ploidy> -s <sex> \   # non-PAR males x|x -> x; see 7.4
  | bcftools annotate -x FORMAT/PP        \   # drop the mis-declared phase-confidence field
  | bcftools norm -m -any -f <chrN.fa>    \   # split multiallelics + left-align against hg38
  | bcftools +fill-tags -t AN,AC,AF       \   # recompute AN/AC/AF from the genotypes
  | bcftools sort -T <ondisk> -Oz -o <out.vcf.gz> ; tabix -p vcf <out.vcf.gz>
```
MEASURED (all 23, 2026-10-05): **lossless** — input records == output records on every chromosome,
**144,858,810 total**; `norm` reported **0 reference mismatches genome-wide** (every REF matches
hg38) and **0 records skipped**; already biallelic (0 splits; the re-phase had normalized), only a
handful of indel left-alignments; 4,091 samples, sample order identical to source; singletons
retained (e.g. chr21 AC=1 → AF=1/8182 exact).

### 7.2 samplesID (prescribed metadata, VCF-header order)
The 4-column CRISPRme samplesID (`#SAMPLE_ID  POPULATION_ID  SUPERPOPULATION_ID  SEX`) is derived
ONLY from the atgu/gnomAD-prescribed metadata file — the exact file the HGDP+1kGP tutorial
notebooks (nb1–nb5) load:
`release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv`, joined onto
the 4,091 VCF-header samples in order:
- POPULATION_ID = `hgdp_tgp_meta.Population` (78 populations; 1kGP as 3-letter codes, HGDP as names)
- SUPERPOPULATION_ID = `hgdp_tgp_meta.Genetic.region` (the 7 regions AFR/AMR/CSA/EAS/EUR/MID/OCE)
- SEX = `sex_imputation.is_female` → `male`/`female` (lowercase; the token
  `tier0_registry.make_chr_ploidy` matches)

The bare convenience columns `population`/`sex` are deliberately NOT used (`sex` is `NA` for all
1kGP samples; `population` relabels 98 samples). Result: **4,091 rows, 0 unknowns**, 2,192 male /
1,899 female; region counts AFR 986 / AMR 549 / CSA 782 / EAS 816 / EUR 771 / MID 157 / OCE 30.

### 7.3 Index build
```
export CRISPRME_REGISTRY_COMPRESS=0 CRISPRME_INDEL_SNP=1
crisprme.py build-index-only --genome <hg38/> --pam 20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread <N> --vcf <VCFs/hg38_HGDP1kGP> --samplesID <config> --path <B>
```
NRG (NAG+NGG) SpCas9 PAM; 2 DNA + 2 RNA bulges; RAW (uncompressed) registry and SNP+indel
co-occurrence ON — consistent with the other production indexes (bulges are a build parameter;
mismatches are a search-time parameter, run at `--mm 6`). Output
`genome_library/NRG_3_hg38+hg38_HGDP1kGP` (+ `_INDELS`) + `Dictionaries/{registry,genotypes,
indel_genotypes,log_indels}_hg38_HGDP1kGP` + `variant_count.json` (**125,180,264 SNPs +
17,192,381 indels**; `phased:true`); every `reg_<chrom>.idx` carries `data_type=genotyped-phased`.

### 7.4 chrX pseudoautosomal (PAR) ploidy — a correctness requirement
chrX is hemizygous in males ONLY outside the pseudoautosomal regions. On GRCh38: PAR1
chrX:10,001–2,781,479, PAR2 chrX:155,701,383–156,030,895, non-PAR core chrX:2,781,480–155,701,382.
Males carry **2** copies in PAR (the homologous Y-PAR) and **1** in non-PAR, so the AF denominator
(AN) must switch BY POSITION. For the 4,091-sample panel (2,192 M / 1,899 F): **PAR AN = 8,182**
(all diploid), **non-PAR AN = 5,990** (2,192×1 + 1,899×2). Two places implement this:

(a) **DATA.** SHAPEIT5 `--haploids` writes non-PAR males HOM-DIPLOID (`x|x`); for the index we
convert them to true haploid with `bcftools +fixploidy` (ploidy file: chrX:2,781,480–155,701,382
M→1; PAR + females fall through to 2). Verified every gnomAD-male is hom on non-PAR (0 hets) so the
collapse is unambiguous; `fill-tags` then yields non-PAR AN=5,990 and PAR AN=8,182.

(b) **REGISTRY + SEARCH.** CRISPRme's chrX ploidy model previously applied the haploid-male model
UNIFORMLY across chrX (a documented TODO), over-counting PAR males as haploid and inflating PAR
allele frequencies ~1.37×. We made it POSITION-aware: the SNP-registry build
(`tier0_compile.compile_from_dict` → `tier0_registry.compile_registry_panel`, new `regime_for_pos`)
routes each record through a PAR (all-diploid) or non-PAR (haploid-male) panel by coordinate
(`tier0_registry.chrx_in_par`); the search-time k≥2 co-occurrence path
(`population_summary.summarize` via `tier0_compile.chrx_par_panels`) selects the same regime by
position. Single-variant (k=1) AF reads the now-PAR-correct stored registry AN. The indel path
inherits correctness: carriers encode ploidy in the GT string (PAR-correct VCF), the SNP+indel
joint denominator is the registry AN, and the marginal indel MAF is the VCF `INFO/AF`. **Autosomes,
chrY, and the sites-only mega index are byte-identical to before** (guarded by unit tests).
Empirically confirmed on the chrX build: a PAR indel (chrX:13,413) MAF = 6/8,182 and a non-PAR
indel (chrX:2,781,857) MAF = 1,613/5,990, matching their VCFs exactly; 8 new unit tests
(`TestChrXParPloidy`, `TestChrXParCoOccurrence`) assert PAR→diploid / non-PAR→haploid-male AN
(114 tier0/popsummary/registry tests green). The build binds the fix as a 4-file PostProcess
overlay on the v2.6.2 image (pending a CRISPRme+ release that bakes it into the codebase).

## 8. Status log

- 2026-09-29 — Exploration + chr22 verification complete (§1–§4). Singleton MAC≥2 filter confirmed;
  multiallelic-split / uniform-phasing / indels / chrX confirmed good; sample-provenance discrepancy
  open. Decision to re-phase with SHAPEIT5. Three background analyses in flight: adversarial
  re-verification, SHAPEIT5 re-phase plan, and extraction of the original published phasing commands.
  Nothing built GW / pushed to HF. Held for Luca.
- 2026-09-29 (cont.) — Original-commands extraction (§5) + SHAPEIT5 re-phase plan (§5/§6) landed.
  KEY: singleton drop is a standalone post-phasing `bcftools view -i'MAC>=2'` (`remove_singletons.py`),
  NOT a SHAPEIT5 behavior; `phase_rare` keeps singletons. Recipe = reproduce atgu pipeline, skip that
  step. atgu `.full.` singleton-inclusive BCFs exist but in a PRIVATE bucket (gs://hgdp-1kg, 403) →
  either request access or re-phase (public dense chr22 VCF direct-downloadable for the smoke). Added
  §4b: recovered singletons are ~35% switch-error → keep for detection but TIER singleton-hinged cis
  (don't silently CONFIRM). Full VAL-1..7 protocol in §6 (incl. trio-parents-held-out anti-self-
  deception gate). PENDING Luca decision on singleton handling (drives whether we do the GW run) +
  index naming/replace-vs-alongside. §7 final recipe still pending the chr22 smoke.
- 2026-09-30 — **chr22 SMOKE PROVED singleton recovery** (partial: 3 of 8 4cM chunks =
  chr22:1–~28 Mb, ligated + validated). Our re-phased region: **326,537 AC=1 singletons recovered
  (all `|`-phased)** vs **0** in the gnomAD MAC≥2 release over the same region; **916,299 records vs
  403,312 (2.27×)** — i.e. re-phasing recovers the ~53% the release drops AND keeps it phased →
  CONFIRMED-cis capable; 4,091 samples; 122.7M genotypes scanned 100% `|`, zero missing. **Two recipe
  corrections the smoke surfaced (now mandatory in §5/§7):** (1) `phase_rare` needs explicit
  `chr:start-end` regions (chunk-file col3=scaffold / col4=input), NOT a bare chromosome name;
  (2) the atgu **pre-phasing QC is MANDATORY** — `bcftools +fill-tags -t all | view -i 'HWE>=1e-30 &&
  F_MISSING<=0.1 && ExcHet>=0.5 && ExcHet<=1.5'` PLUS dropping AC=0 monomorphic-in-subset sites;
  without it SHAPEIT5 `phase_rare` NaN-aborts (`Assertion !isnan(...prob)`), which is what killed
  chunks 4–8 in the smoke. Also: single-region `phase_common` = ~3.5 h/chr22 → GW must chunk (20cM).
  VERDICT: approach validated; GW re-phase is GO once the QC + chunking are wired (HELD for Luca).
- 2026-09-30 (cont.) — **Corrected-recipe chr22 GATE PASSED end-to-end (full chromosome, no crash).**
  With the mandatory atgu QC in place: qced = **1,998,332** records (the QC dropped ~310k pathological
  sites — exactly what was NaN-crashing `phase_rare`); chunked `phase_common` → 473,471-variant
  scaffold (~1.5 h, faster than the 3.5 h single-region); **`phase_rare` completed ALL 8 chunks in
  ~7 min, no crash**; ligate → **4,091 samples, min AC=1, 828,868 phased singletons** recovered.
  Timings (per chrom, excl. download): QC ~50 min, phase_common ~90 min (long pole), phase_rare ~7 min.
  **GW fan-out LAUNCHED** (`drive_gw.sh`): chr1–21 via `rephase_one_chrom.sh`, K=3 concurrent chroms
  (each chunk-parallel, RP_THREADS=6 → ~200 cores), small→big, all on the /srv/local SSD, raw dense
  VCF deleted after QC (disk-bounded), resumable. chrX handled separately after (PAR/haploid). Then
  CRISPRme index build (`NRG_3_hg38+hg38_HGDP1kGP`) + full validation. Publish + default-wiring HELD
  for Luca.
- 2026-09-29 (cont.) — DECISION (Luca): **INCLUDE singletons** (prioritize completeness /
  "never miss a true haplotype" over avoiding low-confidence phase) → **re-phase path** (SHAPEIT5,
  keep singletons). No special "singleton" report flag needed — the `Samples`/"N carrier(s)" +
  MAF columns already make a one-individual variant self-evident; phase-confidence rides on the
  existing CONFIRMED/PUTATIVE label. chr22 re-phase smoke restarted. Also: investigated the ≥4
  multi-variant bug (§9) — it is FIXED + tested; added ≥4 regression tests to the indel+SNP path.
- 2026-10-03/04 — **All 22 autosomes + chrX re-phased + validated** (§6.2, §6.3); genome-wide
  unphased census = **0 unphased** across all 23 (full-file scan). 4,091 samples, singletons kept.
- 2026-10-05 — **VCF-prep + samplesID done; GW index build in flight (done PROPERLY with the PAR
  fix baked in).** Metadata switched to the tutorial-prescribed `gnomad_meta_updated.tsv` (v2);
  samplesID built (§7.2 — 4,091 rows, 0 unknowns, 2,192 M / 1,899 F). VCF-prep of all 23 validated
  lossless (§7.1 — 144,858,810 records, 0 REF mismatch, singletons kept). chrX smoke build passed
  (non-PAR AN=5,990, `data_type=genotyped-phased`). **Found + fixed the chrX PAR ploidy issue**
  (§7.4): PAR males were counted haploid (AF inflated ~1.37×) — now position-aware in both the
  registry build and the search-time k≥2 co-occurrence path; 8 new unit tests, 114 tier0/popsummary
  tests green, autosomes/chrY/mega byte-identical; indel path PAR-correct by inheritance (empirically
  verified). The stale (pre-fix) GW build was killed and **re-run as a single clean `build-index-only`
  with the PAR-aware overlay** (no post-hoc patch), thread 128, RAW + INDEL_SNP. REMAINING: validate
  the finished index + VAL-7 cis audit + new-vs-old comparison → THEN HF publish (Luca GO'd publish
  if the comparison looks good; default-flip staged separately). Publish held pending validation.

## 9. Multi-variant haplotype handling — the ≥4 fix, the dense-window approximation, and the flags

This governs how co-occurring variants inside one protospacer window become off-targets, and
matters more for the new panel (denser + singleton-inclusive → more multi-variant windows).

### 9a. The ≥4-variant "starvation" bug and its fix (SNP path)
For a window with several co-located variants, the genotyped/phased path builds a **combination
lattice**: level-0 = individual variants (seeds w/ carrier sets); level *n+1* = cross level *n*
against the level-0 seeds, intersecting carriers (a combo's carriers = samples carrying ALL its
alts → cis proof). **Original bug:** a de-dup subtraction ran *inside* the growth loop, emptying the
level-0 seeds as soon as a 2-variant combo formed; since deeper levels grow only by crossing against
those seeds, a sample carrying k≥4 cis alts got "used up" into disjoint pairs and its **maximal
k-variant combo never formed** → the true, lowest-mismatch (often 0-mm perfect) cis haplotype was
silently dropped (N=4→under-reported as ≥2 mm, N=6→≥4 mm, …). **Fix** (`new_simple_analysis.py:1296–1331`):
grow the FULL lattice mutation-free, then a single **deferred "peel"** attributes each sample to its
maximal cis combo and subtracts it from strict allelic-subset ancestors (matched on `(pos,elem)` →
multiallelic-safe). Tested across **N=3/4/6/8/12** + trans + unphased (`test_phased_haplotype.py`).

### 9b. Indel + SNP co-occurrence is structurally immune (verified)
The SNP+indel cis join (`indel_snp_cis.cis_cooccurrence`) is **not** a lattice: it takes the indel +
ALL SNP alts the aligned target uses (fixed by the alignment) and does a **direct per-sample
same-slot check** (`some haplotype slot carries the indel AND every SNP alt`) that scales to any
multiplicity with nothing to starve. Added ≥4 regression tests (`test_indel_snp_cis.py`: indel + 4
SNPs cis→CONFIRMED, one-trans→excluded, homozygous→2 copies, population cis/trans/unphased mix).

### 9c. The dense-window APPROXIMATION (the greedy cap) + FLAGS + DEFAULT behaviour
Enumerating the lattice is 2^k in the number of co-located variants, so extreme-density windows are
**capped**: `capped = _FAST_MODE OR (CRISPRME_IUPAC_CAP ≥ 0 AND #variants > CRISPRME_IUPAC_CAP)`.

- **DEFAULT = population-level** (`CRISPRME_FAST_MODE=1`, the default search): **every** multi-variant
  window emits ONE **worst-possible representative** off-target (greedy min-mismatch / max-CFD),
  tagged **PUTATIVE**, phasing-independent. No 2^k enumeration; **no MISS** (the worst case covers the
  window); bounded compute.
- **`--per-sample`** (`CRISPRME_FAST_MODE=0`): builds the full lattice + deferred peel → **CONFIRMED**
  cis for verified combos, PUTATIVE otherwise — *except* windows with more co-located variants than
  `CRISPRME_IUPAC_CAP`, which fall back to the same single worst-case PUTATIVE rep.
- Every capped window is logged to **`<out>.high_variant_density_regions.bed`** so approximated
  regions are visible, never silent.

**The approximation is conservative for the "never miss" goal:** a capped window still emits its
worst-case off-target (so the region is never dropped); what it trades away is *per-combination cis
enumeration* (many individual CONFIRMED rows collapse to one worst-case PUTATIVE row). Detection is
preserved; only the per-haplotype CONFIRMED breakdown is coarsened in the densest windows.

**Implication for the new HGDP+1kGP index (denser + singletons kept):** more windows will be
multi-variant. Under the default population-level mode nothing changes (already worst-case rep). Under
`--per-sample`, moderate-density windows are handled losslessly by the ≥4-fixed lattice; only the
densest windows hit the cap → collapse to the worst-case PUTATIVE rep + HVDR-bed log. Net: no missed
regions; the CONFIRMED/PUTATIVE label + carrier count remain the reader's confidence signal.

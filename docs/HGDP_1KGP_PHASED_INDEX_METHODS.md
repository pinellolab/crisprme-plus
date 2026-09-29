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

## 7. Final validated build recipe (→ METHODS / manuscript)  *(PENDING)*

The finalized, reproduced-and-validated pipeline (re-phasing + CRISPRme index build
`NRG_3_hg38+hg38_HGDP1kGP`) goes here and is copied into `METHODS.md` + the manuscript Methods.

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

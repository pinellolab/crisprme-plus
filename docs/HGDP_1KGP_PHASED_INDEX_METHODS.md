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

## 5. Re-phasing recipe (SHAPEIT5)  *(PENDING — being finalized)*

To be filled from: (a) the exact published commands used to produce `phased_haplotypes_v2`
(from the `atgu/hgdp_tgp` phasing code + paper Methods/supplement — including the specific parameter
that dropped singletons), and (b) the official SHAPEIT5 tutorials. Will include: dense-callset
acquisition, normalization (keep singletons), `phase_common` scaffold, `phase_rare` (singletons),
`ligate`, GRCh38 genetic map, trio/pedigree usage, chunking, and compute/RAM/time.

## 6. Validation / comparison protocol  *(PENDING)*

To be filled: singleton-recovery confirmation (min AC=1 + count recovered), switch-error rate vs
1000G trios (`SHAPEIT5 switch`), common-variant haplotype concordance vs `phased_haplotypes_v2`,
per-chrom variant counts, cross-check vs Zenodo 18156285, and internal source-GT ground-truth for a
handful of CONFIRMED-cis HGDP carriers. Acceptance criteria gate the genome-wide run + index build.

## 7. Final validated build recipe (→ METHODS / manuscript)  *(PENDING)*

The finalized, reproduced-and-validated pipeline (re-phasing + CRISPRme index build
`NRG_3_hg38+hg38_HGDP1kGP`) goes here and is copied into `METHODS.md` + the manuscript Methods.

## 8. Status log

- 2026-09-29 — Exploration + chr22 verification complete (§1–§4). Singleton MAC≥2 filter confirmed;
  multiallelic-split / uniform-phasing / indels / chrX confirmed good; sample-provenance discrepancy
  open. Decision to re-phase with SHAPEIT5. Three background analyses in flight: adversarial
  re-verification, SHAPEIT5 re-phase plan, and extraction of the original published phasing commands.
  Nothing built GW / pushed to HF. Held for Luca.

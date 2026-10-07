# Why re-phasing HGDP mattered: phased HGDP+1kGP vs the unphased hybrid

**Guide** `ACTGAAATCTGTAAGCAGGC` · **PAM** NRG (SpCas9) · **search** mm6 + 1 DNA + 1 RNA
bulge (`6+2`), `--per-sample`, all functional/cancer annotations, both scorers (CFD +
CRISPR-Bulge) · genome hg38. Identical parameters on both indices.

| | OLD — hybrid | NEW — re-phased |
|---|---|---|
| Index | `NRG_3_hg38+hg38_1000G2021_HGDP` | `NRG_3_hg38+hg38_HGDP1kGP` |
| Panel | 1000G-2021 (3,202, **phased**) + HGDP (929, **UNphased**) | gnomAD HGDP+1kGP, 4,091, **jointly SHAPEIT5 re-phased** |
| HGDP co-occurrence | can only be **PUTATIVE** (cis unprovable) | **CONFIRMED** cis + named carriers |

## The motivation (the story)

The previously shipped default (`…1000G2021_HGDP`) is a **hybrid**: the 1000-Genomes
portion is phased, but the **HGDP portion is unphased**. CRISPRme reports a co-occurring
off-target (two nearby variants needed together) as **CONFIRMED** only when it can prove the
variants sit *in cis* on one haplotype in a real individual — which **requires phased
genotypes**. On the hybrid index, every HGDP-involving co-occurrence therefore collapses to
**PUTATIVE** ("possible, but cis unproven"), with no named HGDP carrier. To fix this we
re-phased the gnomAD HGDP+1kGP callset end-to-end with SHAPEIT5 (so **both** cohorts are
phased), rebuilt the index, and re-ran the identical search.

## What changed — headline

**Re-phasing's impact is in co-occurrence interpretation, not top-site reshuffling.**

### 1. The CONFIRMED/PUTATIVE flip (this is why phasing mattered)

| phase_confirmation | CONFIRMED | PUTATIVE | % CONFIRMED |
|---|---|---|---|
| OLD (hybrid) | 8,393,667 | 22,515,383 | 27.2% |
| **NEW (re-phased)** | **36,319,150** | **971,014** | **97.4%** |

PUTATIVE collapsed **23×** (22.5M → 0.97M); CONFIRMED rose **4.3×**. The residual 0.97M
PUTATIVE in NEW is *correct* — genuinely-unobserved cis combinations (a worst-case haplotype
no individual carries), which should stay PUTATIVE. Per chromosome the flip is uniform on the
autosomes (chr1/7/13 all ~99.6–99.9% CONFIRMED); chrX improves 5× (106,204 → 526,278
CONFIRMED) but stays partly PUTATIVE by design (haploid non-PAR worst-case enumeration).

### 2. The top-1000 off-targets are STABLE (the switch is safe)

Comparing the two reports' top-1000 single-site off-targets (ranked by CFD):

| Metric | Value |
|---|---|
| Shared top-1000 sites (chrom:pos:strand) | **793 / 1000 (79.3%)** |
| Shared sites with **identical** CFD (<5e-4) | **767 / 793 (96.7%)** (max ΔCFD 0.22) |
| Observed-category change on shared sites | 520 reference→reference, 254 carrier→carrier, 17 reference→carrier, 2 carrier→reference |

So the **actionable, highest-risk predictions do not move** when switching to the phased
index — reassuring for making it the default. The 207 site differences and the richer
variant-driven fraction (NEW 478 vs OLD 388 variant-driven in the top-1000) are driven by
**panel composition** (gnomAD's QC'd variant set differs from the raw 1000G+HGDP set), **not
by phasing** — phasing is orthogonal to single-site detection and scoring.

### 3. Ground-truth: a CONFIRMED cis call traced to the phased genotypes

Off-target `chr13:100001388(-)` uses two co-occurring SNPs `chr13_100001399_A_T` +
`chr13_100001406_G_A` (joint AF 9/8182 = 0.00110). Carriers include the HGDP sample
`LP6005443-DNA_E02`, whose phased genotypes are **`0|1` / `0|1`** — both ALT alleles on the
**same haplotype (cis)**, confirmed directly in the re-phased VCF. On the unphased hybrid this
same locus could only be PUTATIVE.

### 4. Detection is a superset, with one documented caveat

NEW total off-targets 5,996,092 vs OLD 5,540,022 (**+8.2%**); NEW ≥ OLD in every
mismatch/bulge bucket. Caveat (panel-composition, **not** a phasing or software effect): ~192k
*variant-driven* OLD off-targets (incl. 333 with CFD ≥ 0.2) are **absent** from NEW because
their driving alleles are not in gnomAD's QC'd re-phased callset; reference-genome detection is
99.997% preserved and NEW adds ~4× more high-CFD sites than it drops.

## Take-home for the paper

Re-phasing the HGDP cohort did **not** reshuffle the top off-target hits — those are robust to
both the panel change and phasing. Its value is a **completeness/correctness** gain in the
*interpretation* of co-occurring variants: it converts 22.5M → 0.97M PUTATIVE calls into
**CONFIRMED cis** with named carriers spanning HGDP's seven genetic regions and exact joint
allele frequencies — information the unphased hybrid simply could not provide. The default flip
is therefore safe (stable top sites) and strictly more informative (confirmed cis + carriers).

## Artifacts & reproducibility

- **Reports (full, 1.4 GB each):** `NEW_phased_HGDP1kGP_report.zip`,
  `OLD_hybrid_1000G2021_HGDP_report.zip`
- **Top-1000 tables:** `NEW_phased_HGDP1kGP_top1000.tsv`, `OLD_hybrid_1000G2021_HGDP_top1000.tsv`
- **Phase tally:** `phase_confirmation_tally.txt`
- Durable copy on ml007: `/srv/local/lp698_PAPER_phasing_importance/`
- Index build recipe: `seq_script/merge_panels/hgdp1kgp_build.sh`; re-phasing + validation
  methods: `docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`
- Published index: HuggingFace `lucapinello/crisprme-data` →
  `indexes/NRG_3_hg38+hg38_HGDP1kGP.tar.gz` (+ `genotypes_hg38_HGDP1kGP.tar.gz`)

Replicate the search on either index:

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt --guide guide.txt \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --annotation Annotations/dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz \
  --gene_annotation Annotations/gencode.protein_coding.bed.gz \
  --mm 6 --bDNA 1 --bRNA 1 --per-sample --output cmp --thread 32 \
  --index-path genome_library
```

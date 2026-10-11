# Precomputed CRISPRme indexes on HuggingFace

> **Reproducing the shipped production indexes.** The end-to-end recipes for the
> released variant indexes — the phased default **`NRG_3_hg38+hg38_HGDP1kGP`** (HGDP +
> 1000 Genomes), the single-source phased **`NRG_3_hg38+hg38_1000G2021`**, the HPRC pangenome
> **`NRG_3_hg38+hg38_HPRC`**, and sites-only **`NRG_3_hg38+hg38_mega`** (5 sources) —
> live in [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md) (Mode 1,
> the HPRC section with `hprc_build.sh`, and Mode 2, with the as-run `mega_gw_merge.sh` +
> `mega_build_indels.sh` drivers). This page covers the generic single-dataset
> build/publish/download flow. **Naming note:** `--dictless` is a *publish flag* (it drops the
> ~152 GB per-sample SNP dictionaries; see below), **not** a name marker — the shipped index
> names are simply `<pam>_<N>_<ref>+<vcf>` (`NRG_3_hg38+hg38_HGDP1kGP`,
> `NRG_3_hg38+hg38_1000G2021`, `NRG_3_hg38+hg38_HPRC`, `NRG_3_hg38+hg38_mega`).

Bulge-enabled CRISPRme searches need a CRISPRitz **index** of the reference
genome. Building it is the single most expensive one-time step of a search. That
index depends only on the **genome + PAM + bulge count** (not on the guides or
the variant dataset), so it can be built once and reused — or built once by a
maintainer, published to HuggingFace, and downloaded by everyone else.

A **variant-aware** index additionally bundles everything the variant search
needs at post-analysis time so it works with **no source VCFs** and **no
separate samples download**:

- the additive **dict-less tiers** — a small Tier-0 **registry** (`registry_<vcf>/`,
  out-of-the-box off-target detection with corrected AF/rsID) plus the larger
  Tier-1 **genotype store** (`genotypes_<vcf>/`, per-sample `Samples`), which
  together replace the ~152 GB per-sample SNP dictionaries; and
- the **samplesID** files the search's `--samplesID` listing expects (the
  combined `<vcf>.samplesID.txt` plus the per-db `<ref>_<db>.samplesID.txt`
  lists), so a downloaded variant index is self-complete.

This document covers the download → build → publish workflow. It complements
Section 3.5 of the data-setup guide (`docs/crisprme_data_setup_051826.md`).

## The four production indexes

SpCas9 **NRG** (NAG+NGG) variant indexes ship prebuilt on HuggingFace — pick
by whether you want the broad fully-phased default, a clean single-source panel,
pangenome coverage, or the widest allele-frequency provenance. They are
**complementary**, not alternatives:

| Index name | Sources | Data model | Co-occurrence confidence |
|---|---|---|---|
| `NRG_3_hg38+hg38_HGDP1kGP` *(recommended default)* | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000 Genomes), jointly SHAPEIT5-re-phased | **sample-level genotypes, fully phased** | **CONFIRMED** phased cis throughout with **named carrier samples** + exact joint AF (including HGDP's seven genetic regions, with correct chrX PAR allele frequencies) — the broadest fully-phased panel of real observed haplotypes |
| `NRG_3_hg38+hg38_1000G2021` | 1000 Genomes 2021 (3,202 samples) | **sample-level genotypes, fully phased** | **CONFIRMED** phased cis throughout with **named carrier samples** + exact joint AF (all genotyped-phased, no PUTATIVE-from-phasing) — a clean single-source panel of real observed haplotypes |
| `NRG_3_hg38+hg38_HPRC` | HPRC Release 2 pangenome `hprc-v2.0-mc-grch38` (232 assembly-derived genomes incl. CHM13) | **sample-level genotypes, phased** (graph-derived) | **CONFIRMED** phased cis with **named carrier samples** + exact joint AF; captures **pangenome-specific** variation absent from short-read panels |
| `NRG_3_hg38+hg38_mega` | 5 sources — 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us | **sites-only** (allele frequencies, no shared samples) | **PUTATIVE** with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`); a co-occurring haplotype here **may or may not exist in any real individual** — no carriers |

**Choosing an index.** You only need **one**. For most use cases pick the **recommended
default `NRG_3_hg38+hg38_HGDP1kGP`** — HGDP + 1000 Genomes (4,091 samples), jointly
SHAPEIT5-re-phased, so every reported co-occurrence is **CONFIRMED** cis with **named
per-sample carriers** and exact joint AF throughout, including HGDP's seven genetic regions
and with correct chrX PAR allele frequencies: the broadest fully-phased panel of real
observed haplotypes. Pick **`NRG_3_hg38+hg38_1000G2021`** when you want a clean single-source
1000 Genomes 2021 panel (3,202 samples, fully phased → CONFIRMED cis + named carriers). Pick
**HPRC** for assembly-derived / pangenome variation (phased → CONFIRMED cis +
carriers). **Escalate to `mega`** only for a widest-provenance worst-case screen when you
must not miss a rare allele from any of five databases: it is **sites-only** (union of allele
frequencies, no genotypes), so every multi-variant / co-occurring off-target is **PUTATIVE**
(a worst-case reconstruction that may not exist in any real individual) with a min-AF bound
and **no named carriers**.

All four carry **searchable indels genome-wide** and report **SNP+SNP and SNP+indel**
co-occurring off-targets. Download the reference data once, then any index:

```bash
crisprme.py download --what all  --path .                                                 # reference genome + annotations (once)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP       --path .     # recommended default (phased HGDP + 1000 Genomes, 4,091 samples; carriers)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_1000G2021      --path .     # 1000-Genomes-only alternative (single-source, phased; carriers)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HPRC           --path .     # HPRC pangenome (phased assemblies, carriers)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_mega          --path .      # mega (sites-only, per-dataset AF provenance)
```

A first-time unauthenticated download of these many-small-file indexes can hit
HuggingFace's HTTP-429 rate limit; the downloader retries with backoff, and
authenticating (`hf auth login`) avoids it entirely. See METHODS §2 for
the genotyped-vs-sites-only distinction and §4 for how each is scanned.

## Layout on HuggingFace

Indexes live under `indexes/` in the CRISPRme dataset repo (default
`lucapinello/crisprme-data`, override with `--hf-repo` or `CRISPRME_HF_REPO`):

```
indexes/
  NRG_3_hg38.tar.gz                              # SpCas9 (NRG = NAG+NGG) reference index of hg38, up-to-2-bulge (DEFAULT reference index)
  NRG_3_hg38+hg38_HGDP1kGP.tar.gz           # SpCas9 (NRG) dict-less variant-aware index (HGDP + 1000 Genomes, 4,091 samples, phased), up-to-2-bulge (recommended default)
  genotypes_hg38_HGDP1kGP.tar.gz                     # SEPARATE Tier-1 genotype store companion for the HGDP1kGP index (rides along on download)
  NRG_3_hg38+hg38_1000G2021.tar.gz          # SpCas9 (NRG) dict-less variant-aware index (1000G-2021, phased), up-to-2-bulge (1000-Genomes-only alternative)
  genotypes_hg38_1000G2021.tar.gz                    # SEPARATE Tier-1 genotype store companion for the 1000G2021 index (rides along on download)
  NRG_3_hg38+hg38_HPRC.tar.gz                    # SpCas9 (NRG) HPRC pangenome variant index (232 phased assembly-derived genomes)
  genotypes_hg38_HPRC.tar.gz                          # SEPARATE Tier-1 genotype store companion for the HPRC index
  NRG_3_hg38+hg38_mega.tar.gz                    # SpCas9 (NRG) sites-only mega index (5 sources) — NO genotype companion (sites-only)
```

Each variant index stamps its **`data_type`** (`sites-only` / `genotyped-unphased` /
`genotyped-phased` / `hybrid`) and a `phased` flag in every `registry_<vcf>/reg_<chrom>.idx`
manifest, so tooling (e.g. the web variant-dataset selector) can label an index's type
and phasing without scanning the multi-GB per-sample store — `HGDP1kGP` = **genotyped-phased**,
`1000G2021` = **genotyped-phased**, `HPRC` = **genotyped-phased**, `mega` = **sites-only**
(a merged panel mixing phased and unphased sources would be **hybrid**).

The **NRG** default matches SpCas9's broad recognition (NAG + NGG), so variant-created
NAG off-targets (e.g. the CPS1 off-target from the CRISPRme paper) are found out of the
box.

**Two artifacts per variant index:**

- **Main tarball** (`NRG_3_hg38+hg38_HGDP1kGP.tar.gz`) — the index
  itself (`<name>/`), its `_INDELS` companion, the indel logs
  (`Dictionaries/log_indels_<vcf>/`), the Tier-0 `Dictionaries/registry_<vcf>/`,
  the samplesID lists (`samplesIDs/<vcf>.samplesID.txt` + the per-db
  `samplesIDs/<ref>_<db>.samplesID.txt`), and `manifest.json` at the archive
  root. In a **classic** (non-dict-less) publish the main tarball ALSO carries
  the per-sample SNP dicts (`Dictionaries/dictionaries_<vcf>/`).
- **Genotype companion** (`genotypes_hg38_HGDP1kGP.tar.gz`) — the big Tier-1
  store, uploaded separately under the same `indexes/` prefix. `download --what
  index` fetches it automatically unless `--no-genotypes` is given.

**Naming — `--dictless` is a flag, not a name marker.** `publish-index --dictless`
controls *what goes in the tarball* (it drops the ~152 GB per-sample SNP
dictionaries, keeping the compact Tier-0/Tier-1 stores), but it does **not** decorate
the index name. A variant index is published and installed under the plain convention
`<pam>_<N>_<ref>+<vcf>` — e.g. `NRG_3_hg38+hg38_HGDP1kGP`, extracting to
`genome_library/NRG_3_hg38+hg38_HGDP1kGP/`. The search resolves an index by that
same convention (the `<ref>` segment == the genome-folder basename; the `+<vcf>` segment
is shared verbatim with the `genotypes_<vcf>` companion). A reference-only index (no `+`)
unpacks to a single `genome_library/<name>/` directory plus `manifest.json` and is used
with no extra steps.

## Dict-less flow (variant-aware index)

### Build (maintainers)

Build the variant-aware index with `--vcf` and `--samplesID`:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 \
  --vcf VCFs/hg38_HGDP1kGP --samplesID samplesIDs.config.txt \
  --path "$CRISPRME_DIR"
# -> genome_library/NRG_3_hg38+hg38_HGDP1kGP/          (+ _INDELS companion)
#    Dictionaries/registry_hg38_HGDP1kGP/              (Tier-0 allele-freq registry)
#    Dictionaries/registry_hg38_HGDP1kGP/variant_count.json  (SNP + indel counts, for the report)
#    Dictionaries/genotypes_hg38_HGDP1kGP/             (Tier-1 per-sample genotype store)
#    Dictionaries/log_indels_hg38_HGDP1kGP/            (indel logs, for indel post-analysis)
#    samplesIDs/hg38_HGDP1kGP.samplesID.txt            (combined, emitted by the build)
```

**One command, all supporting files.** `build-index-only` writes *everything* a
variant-aware search + report needs — the index and its `_INDELS` companion, the
Tier-0 registry, the Tier-1 genotype store, the indel logs, the combined samplesID
list, and the `variant_count.json` manifest (`n_records` SNPs + `n_indels`, which
feeds the report's *Variants included* line). Nothing else has to be assembled by
hand, so `publish-index` (below) ships a self-complete index.

> **Your own VCF.** `--vcf` accepts a folder of `.vcf.gz`. It need not be pre-split: a **single
> multi-chromosome VCF is auto-split per chromosome** before enrichment (the enricher requires
> per-chromosome files), so a user-provided VCF — including one registered through the web
> **Settings → Add a VCF dataset** card and built with **Build index** — works without manual
> preparation. A **SNP-only VCF** (no indels) builds a complete SNP variant index; the indels index
> is simply skipped (SNP+indel co-occurrence does not apply to an indel-free dataset).

`--samplesID` is a listing file (one samplesID filename per line, resolved under
`samplesIDs/`); a combined panel lists **both** the 1000G and HGDP files, e.g.:

```
hg38_1000G.samplesID.txt
hg38_HGDP.samplesID.txt
```

> **`--samplesID` is required for the dict-less / self-complete flow.** The
> dicts still build without it, but you silently get a **dicts-only** index: no
> Tier-0 registry, no Tier-1 genotype store (so no fast post-analysis), and — for
> a merged panel — no combined `samplesIDs/<vcf>.samplesID.txt`. Pass
> `--samplesID` to emit the tiers **and** (for a merged panel) write the combined
> samplesID list into the install, so the published index is self-complete.

### Merged multi-dataset panels & phasing (maintainers)

To build an index over **several datasets at once** (e.g. 1000G + HGDP, or adding
a new cohort later), first merge the per-dataset VCFs into one combined panel, then
point `build-index-only --vcf` at the merged VCF folder and list every per-db
samplesID file in `--samplesID`.

```bash
# 1) normalize each source VCF to per-alt records (left-align + split multiallelics)
for db in 1000G HGDP NewCohort; do
  bcftools norm -m -any -f hg38.fa "$db.vcf.gz" -Oz -o "$db.norm.vcf.gz" && bcftools index -t "$db.norm.vcf.gz"
done
# 2) merge into one combined panel (consistent reference build + contig naming across datasets)
#    (hypothetical 3-way merge — pick any name for the combined-panel folder)
bcftools merge 1000G.norm.vcf.gz HGDP.norm.vcf.gz NewCohort.norm.vcf.gz -Oz -o VCFs/hg38_combined_panel/merged.vcf.gz
# 3) build the dict-less index (emits registry + genotypes + variant_count.json + combined samplesID)
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --vcf VCFs/hg38_combined_panel --samplesID samplesIDs.config.txt --path "$CRISPRME_DIR"
```

Two things the build handles for you, which matter as new merged panels are added:

- **Per-dataset provenance is preserved.** Allele frequencies are recomputed on the
  merged multiallelic records (requires the multiallelic-AF fix in CRISPRitz PR #36)
  and reported **per native dataset label** *and* as a combined global frequency over
  the union panel — 1000G / HGDP / gnomAD are never conflated. The panel size (AN)
  is the **genotyped** sample set, not the roster (`samplesID` can over-list samples;
  see METHODS §1). The `databases` block in each `registry_<vcf>/reg_*.idx` carries
  the per-dataset `sample_count`.
- **Phasing is resolved per haplotype from the genotypes themselves.** A multi-variant
  off-target is enumerated as an **observed haplotype**; it is reported *confirmed*
  only when its carriers reached it via **phased** (`|`), same-phase-set genotypes,
  and *putative* when the genotypes are **unphased** (`/`) or span different phase
  sets. So a **mixed** merged panel (some datasets phased, some not) is handled
  correctly with no special flags — phased datasets yield confirmed haplotypes,
  unphased yield putative, each on its own. The dataset-wide population-summary
  phasing flag defaults conservatively to *unphased* (reports bounds, never a false
  confirmation) and can be set explicitly at build time if a dataset's phasing is
  known but not detectable from its GT separators.

### Publish (maintainers)

```bash
export HF_TOKEN=hf_...        # your write token, in the shell only
crisprme.py publish-index --index genome_library/NRG_3_hg38+hg38_HGDP1kGP --dictless
# -> indexes/NRG_3_hg38+hg38_HGDP1kGP.tar.gz     (main: index + _INDELS + registry (+ variant_count.json) + indel logs + samplesIDs + manifest)
# -> indexes/genotypes_hg38_HGDP1kGP.tar.gz      (separate Tier-1 companion)
```

`--dictless` **drops the ~152 GB per-sample SNP dictionaries**
(`dictionaries_<vcf>/`) from the main tarball — the Tier-0 registry + Tier-1
genotype tiers replace them — while **keeping the indel logs** (indel
post-analysis still needs them; the tiers are SNP-only). In BOTH modes the small
`registry_<vcf>/` is added to the main tarball when present, a separate
`genotypes_<vcf>.tar.gz` companion is produced and uploaded when a genotype
store exists, and the samplesID files this index needs are bundled into the main
tarball (so `download --what index` is self-complete). Without `--dictless`,
publishing is byte-for-byte the classic path plus these additive members.

### Download (end users)

```bash
# fetch the variant-aware index (main tarball + genotype companion + bundled samplesIDs)
crisprme.py download --what index \
  --index-name NRG_3_hg38+hg38_HGDP1kGP --path "$CRISPRME_DIR"
# installs under the canonical name genome_library/NRG_3_hg38+hg38_HGDP1kGP/
```

- The index installs under its plain `<pam>_<N>_<ref>+<vcf>` name so the search
  resolves it (there is no name decoration to strip — `--dictless` only affected
  what the maintainer put in the tarball).
- The **combined + per-db samplesID lists** are installed into
  `samplesIDs/` from the main tarball, so no `--what samples` / `--what all` is
  needed. (For a legacy index built before this bundling existed, download falls
  back to fetching the per-db lists from HF and synthesizing the combined file.)
- The **Tier-1 genotype store** is fetched automatically. Add `--no-genotypes`
  to skip the big companion: off-target **detection still works** via the Tier-0
  registry, but per-sample `Samples` are degraded until the store is present.
- The **CLI search-list files** `list_vcf.txt` (the dataset, `hg38_HGDP1kGP`)
  and `list_samplesID.txt` (its combined samplesID) are written at the install
  root, so a CLI search works out of the box — the same lists the web form builds
  per-search. Installing several variant indexes appends each dataset once.

Then search, pointing at that library (or just run from `$CRISPRME_DIR`):

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --guide my_guide.txt --mm 4 --bDNA 2 --bRNA 2 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --index-path "$CRISPRME_DIR/genome_library" \
  --output my_search
```

The search finds the prebuilt index (and its tiers + samplesIDs) and skips the
build/enrichment entirely.

## Building a new-PAM index for a published dataset (HGDP+1kGP, mega)

The prebuilt indexes above are all **SpCas9 NRG**. To screen a **different nuclease / PAM**
(SpCas9 **NGG**, SaCas9 **NNGRRT**, Cas12a **TTTN**, or any IUPAC motif) against the same
re-phased **HGDP+1kGP** or sites-only **mega** variant data, you build that one index once and
reuse it for every later search — you do **not** re-derive the variant data.

### What's reusable vs. what you rebuild

Only the genome index itself is PAM-specific; everything derived from the variants is not:

| Artifact | PAM-specific? | For a new PAM |
|---|---|---|
| `genome_library/<PAM>_<N>_hg38+hg38_<ds>/` (+ `_INDELS`) | **yes** | **rebuilt** |
| Enriched genome `Genomes/hg38+hg38_<ds>/` (+ `_INDELS`) | no | reuse (skip enrichment) |
| `Dictionaries/` (registry, genotypes, indel-GT, logs) | no | reuse (from the normal index download) |
| `samplesIDs/` | no | reuse |

> **Why you need an extra file.** The HuggingFace **index** tarball ships the *NRG* index + its
> `Dictionaries/` + `samplesIDs/` — it does **not** contain the enriched genome, so it alone can't
> build a different PAM. The PAM-independent **enriched-genome kit** is published separately under
> **`build-inputs/<dataset>.tar.gz`** in the same repo (`lucapinello/crisprme-data`). Download it and you
> skip straight to indexing the new PAM. If instead you want to **reproduce the kit** (or build one
> for a **different VCF**), the full source-to-kit workflow is documented below.

### Prerequisites

- CRISPRme+ installed; the **hg38 reference** (`crisprme.py download --what all --path .`).
- `hf` (the HuggingFace CLI, ships with `huggingface_hub`), `bcftools`/`htslib`.
- Disk/RAM: the **fast path** needs a few GB + the index output. The **full path** re-runs
  enrichment — a genome-wide genotyped build writes **>100 GB of temp** and needs ample RAM
  (see the "Scratch space" note in [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md)).

### Step 1 — choose or author the PAM file

A PAM file is one line, two whitespace-separated columns — `<IUPAC_motif> <signed_position>` —
where a **positive** position means a **3′** PAM (the last *N* characters of the motif) and a
**negative** one a **5′** PAM (the first *N*). Full spec: [`docs/INPUT_FORMATS.md`](INPUT_FORMATS.md) §1.
Ready-made files already ship in `PAMs/`:

| Nuclease | `PAMs/` file | Contents |
|---|---|---|
| SpCas9 (NGG) | `20bp-NGG-SpCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNGG 3` |
| SpCas9 (NRG = NAG+NGG) | `20bp-NRG-SpCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNRG 3` |
| SaCas9 (NNGRRT) | `21bp-NNGRRT-SaCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNNNGRRT 6` |
| Cas12a (TTTN, 5′) | `23bp-TTTN-Cas12a.txt` | `TTTNNNNNNNNNNNNNNNNNNNNNNNN -4` |

For any **other** nuclease, write one line yourself — e.g. a hypothetical 20 nt-spacer 3′ `NGA`:

```bash
printf 'NNNNNNNNNNNNNNNNNNNNNGA 3\n' > PAMs/20bp-NGA-custom.txt
```

### Step 2 — download the enriched-genome kit

```bash
# ~1 GB/dataset download (a few GB once extracted): the SNP-enriched genome + its fake-indel companion, as one tarball.
hf download lucapinello/crisprme-data build-inputs/hg38_HGDP1kGP.tar.gz \
  --repo-type dataset --local-dir .
tar xzf build-inputs/hg38_HGDP1kGP.tar.gz
#  -> Genomes/hg38+hg38_HGDP1kGP/         (SNP-enriched, hg38 + IUPAC overlay)
#  -> Genomes/hg38+hg38_HGDP1kGP_INDELS/  (fake-indel genome)
```

Check the tarball against the `SHA256SUMS` shipped alongside it under `build-inputs/`. (For
**mega**, substitute `hg38_mega`.) The PAM-independent `Dictionaries/` + `samplesIDs/` you need at
**search** time come from the normal index download (`download --what index --index-name
NRG_3_hg38+hg38_<ds>`).

### Step 3 — build the new-PAM index (reuse the enriched genome, skip enrichment)

Point `crispritz index-genome` at the downloaded enriched genome with your new PAM. This builds the
PAM-specific index directly from the IUPAC-enriched FASTA — **no re-enrichment** (verified: a new
NGG index off the HGDP+1kGP enriched genome builds in seconds per chromosome):

```bash
# SNP index -> genome_library/<PAM>_<N>_hg38+hg38_<ds>/
crispritz.py index-genome hg38+hg38_<ds> Genomes/hg38+hg38_<ds>/ PAMs/<NEW_PAM>.txt \
  -bMax 3 -th <T>          # -bMax = N = max(bDNA,bRNA)+1 = 3 (the same N as in the index name, e.g. NRG_3)
# indel index -> genome_library/<PAM>_<N>_hg38+hg38_<ds>_INDELS/
python PostProcess/pool_index_indels.py \
  Genomes/hg38+hg38_<ds>_INDELS/ PAMs/<NEW_PAM>.txt <PAM_MOTIF> hg38 hg38_<ds> 3 <T>
```

Here `<PAM_MOTIF>` is the PAM region token itself — the IUPAC motif without the spacer `N`s
(e.g. `NGG`, `NRG`, `NNGRRT`, `TTTN`) — the same prefix the index name uses
(`<PAM_MOTIF>_<N>_hg38+hg38_<ds>`).

> **One-command alternative.** With the enriched genome already in place, `crisprme.py
> build-index-only --genome Genomes/hg38 --pam PAMs/<NEW_PAM>.txt --bDNA 2 --bRNA 2 --vcf
> VCFs/hg38_<ds> …` detects it (prints `Enriched genome already present`) and skips straight to the
> same two indexes. Use it if you also have the dataset's `Dictionaries/`/`VCFs/` staged; otherwise
> the two `crispritz`/`pool_index_indels` calls above need only the enriched kit.

### Step 4 — search with the new-PAM index

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/<NEW_PAM>.txt \
  --guide my_guide.txt --mm 4 --bDNA 2 --bRNA 2 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --index-path "$CRISPRME_DIR/genome_library" --output my_search
```

`complete-search` resolves the new index by its `<PAM>_<N>_hg38+hg38_<ds>` name and reuses the
sibling `Dictionaries/` + `samplesIDs/` (from the normal index download) — no rebuild.

## Reproducing the enriched-genome kit from source (or building one for a new VCF)

The published kit is produced **once** from the source variants; reproduce it (or make a kit for
your **own** VCF) with the full workflow below. Each dataset's authoritative, validated recipe lives
in its own document — summarized here with the key provenance and filters, and linked for the
complete protocol.

### HGDP + 1000 Genomes (re-phased, genotyped)

1. **Data source.** The gnomAD **v3.1.2 HGDP+1kGP** jointly-called callset (GRCh38),
   **4,091 samples** (925 HGDP + 3,166 1000 Genomes, post-QC) from the gnomAD public release.
2. **Re-phasing (why our kit differs from the public release).** We re-phase with **SHAPEIT5**
   (`phase_common` scaffold + `phase_rare`) to **keep singletons** — the public
   `phased_haplotypes_v2` release drops variants at MAC≥2 — so rare alleles survive for detection.
   chr1–22 + chrX; chrX uses a 3-way PAR split with non-PAR **males haploid** (`+fixploidy`;
   non-PAR AN=5990, PAR AN=8182). Full protocol + validation:
   [`docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`](HGDP_1KGP_PHASED_INDEX_METHODS.md) §5–§7.
3. **samplesID** is derived from the gnomAD-prescribed `gnomad_meta_updated.tsv` in VCF-header
   order (POP = Population, SUPERPOP = the 7 Genetic regions, SEX = is_female).
4. **VCF-prep (per chrom).** `bcftools +fixploidy` (chrX) → `annotate -x FORMAT/PP` →
   `norm -m -any -f <chr>.fa` (split + left-align) → `+fill-tags -t AN,AC,AF` (recompute AF) →
   `sort` → `bgzip`+`tabix`.
5. **Enrich + index** (produces the kit):
   ```bash
   export CRISPRME_REGISTRY_COMPRESS=0 CRISPRME_INDEL_SNP=1
   crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
     --bDNA 2 --bRNA 2 --thread <N> \
     --vcf VCFs/hg38_HGDP1kGP --samplesID hgdp1kgp.samplesID.config.txt --path .
   # the Genomes/hg38+hg38_HGDP1kGP/ (+ _INDELS/) it writes IS the enriched kit
   ```
   As-run script: `seq_script/merge_panels/hgdp1kgp_build.sh`.

### mega (5-source, sites-only)

1. **Sources (GRCh38).** 1000G-2021 (3,202, phased), HGDP (929, unphased), gnomAD v4.1 (sites-only
   AF), TOPMed (sites-only, `AN=0`), All-of-Us (single pseudo-sample). TOPMed/AoU are
   controlled-access — only their **per-source AF** survives into the sites-only index, never
   individual data. Provenance table: [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md).
2. **Merge strategy + filters (per source, per chrom).** `bcftools norm -m -any -f <chr>.fa`
   (split + left-align so the same indel from two sources merges) → **MAF filter**
   `-e 'AF<=0.001 || AF>=0.999'` (keep 0.001 < AF < 0.999) → `view -G` (strip genotypes) →
   `annotate -x '^INFO/AF'` (keep **only** AF) → rename `INFO/AF → AF_<src>`. Then
   `bcftools merge -m none` the sources and add **`AF_max`** (max over present `AF_<src>`). There is
   **no pooled AF** (ΣAC/ΣAN isn't meaningful across these sources). Driver:
   `seq_script/merge_panels/mega_gw_merge.sh`.
3. **Registry (sites-only).** `PostProcess/build_mega_registry.py` →
   `tier0_registry.compile_registry_from_info_af`: each source becomes a db group with
   `AC=round(AF·AN_nom)`, `AN_nom=2×cohort-N`; the GLOBAL group carries `AF_max`. SNP-only; indels
   (~24%) are routed to the fake-indel path.
4. **Searchable indels.** `seq_script/merge_panels/mega_build_indels.sh` adds a synthetic `MEGA`
   pseudo-sample so the enricher materializes fake contigs, builds the `_INDELS` genome, then strips
   the pseudo-sample back out.
5. **Enrich + index** (produces the kit): `build-index-only` **without** `--samplesID` (sites-only
   has no genotypes to count); drop in the registry afterward. Design rationale:
   [`docs/DESIGN_mega_index.md`](DESIGN_mega_index.md).

Once you have an enriched genome (`Genomes/hg38+hg38_<ds>/` + `_INDELS/`) — downloaded or
rebuilt — **any** PAM reuses it via Step 3.

> **Time & space.** Downloading the kit + indexing a new PAM is minutes–~1 h (PAM/bulges/threads).
> Rebuilding the kit from source re-runs the full enrichment (and, for HGDP+1kGP, the SHAPEIT5
> re-phase) — budget hours–days + >100 GB scratch (see the "Scratch space" note in the merge
> README). We publish the enriched kit precisely so the common case skips that.

## Reference-only index (no variants)

A reference index has no `+<vcf>` segment, no dicts/tiers, and no samplesID
bundle. Build and publish are the simple case:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 --path "$CRISPRME_DIR"
# -> genome_library/NRG_3_hg38/

crisprme.py publish-index --index genome_library/NRG_3_hg38
# -> indexes/NRG_3_hg38.tar.gz

crisprme.py download --what index --index-name NRG_3_hg38 --path "$CRISPRME_DIR"
```

## manifest.json

Every published index carries a small provenance manifest inside its tarball.
A reference index:

```json
{
  "name": "NRG_3_hg38",
  "created_at": "2026-08-05T12:00:00+00:00",
  "pam": "NRG",
  "index_bmax": "3",
  "genome": "hg38"
}
```

A dict-less variant index adds the tier/companion/self-completeness fields:

```json
{
  "name": "NRG_3_hg38+hg38_HGDP1kGP",
  "created_at": "2026-08-05T12:00:00+00:00",
  "pam": "NRG",
  "index_bmax": "3",
  "genome": "HGDP1kGP",
  "display_label": "SpCas9 NRG — hg38 (HGDP + 1000 Genomes)",
  "has_registry": true,
  "dictless": true,
  "has_genotypes": true,
  "has_samplesids": true
}
```

> **Note:** `pam` / `index_bmax` / `genome` are parsed from the index *name* and are
> only meaningful for a **reference** index. For a **variant/merged** index the name
> is `<pam>_<N>_<ref>+<vcf>`, so this name-splitting makes `genome` the last VCF token
> (here `HGDP1kGP`), not the reference assembly — a harmless artifact. The authoritative
> dict-less descriptors are the boolean flags: `dictless`, `has_registry`,
> `has_genotypes`, `has_samplesids`.

- `has_registry` — the main tarball carries the Tier-0 `registry_<vcf>/`.
- `dictless` — the per-sample SNP dicts were excluded (tiers replace them).
- `has_genotypes` — a separate `genotypes_<vcf>.tar.gz` companion was uploaded.
- `display_label` — human-friendly name shown by the web index list / search form.
- `has_samplesids` — the samplesID lists are bundled in the main tarball, so
  `download --what index` yields a self-complete install (no `--what samples`).

The manifest is surfaced (build timestamp / display name) when the index is
downloaded, and old fields are ignored by older consumers — the index directory
itself is what `complete-search` consumes.

## Notes

- An index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`;
  a different PAM or a higher bulge count needs its own index.
- Variant-enriched (genome + VCF) indexes are dataset-specific. The default
  HGDP + 1000 Genomes enriched index (`NRG_3_hg38+hg38_HGDP1kGP`) is prebuilt and
  hosted alongside the reference index; for any other cohort, build and publish the
  enriched index the same way (pass `--samplesID` so it is self-complete).
- A published variant index is **self-complete**: the samplesID lists are
  bundled, so `download --what index` + search works without a separate
  `--what samples`. The genotype store is the only optional piece (skip it with
  `--no-genotypes` for detection-only).
- If `--index-path` is given but no matching index is found there,
  `complete-search` fails fast with a clear message rather than silently
  rebuilding — so a missing/wrong download is caught immediately.

## SNP+indel co-occurring off-targets (on by default)

CRISPRme+ searches **SNP+indel co-occurring off-targets** — those that require **both**
a nearby SNP **and** an indel in the same protospacer window — **by default** since
2.5.0. The prebuilt indexes are already built with this on, so a normal `download` +
`complete-search` reports these co-occurrences **with no extra flag** (a class the
classic two-independent-passes search could not see: SNPs on the IUPAC-enriched genome,
indels on a fake-indel genome cut from the *plain* reference).

To reproduce the classic two-independent-passes behavior, set `CRISPRME_INDEL_SNP=0`
**before both the build and the search**. With the gate off, builds and searches are
byte-identical to classic CRISPRme.

```bash
# default (co-occurrence ON) — nothing to set:
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
    --bDNA 2 --bRNA 2 --vcf VCFs/hg38_HGDP1kGP --samplesID samplesID.listing.txt --path ./
# to opt out and reproduce classic CRISPRme (byte-identical), set the gate to 0 first:
#   export CRISPRME_INDEL_SNP=0   # before BOTH build and search
```

With co-occurrence on, the build (a) compiles a per-chromosome **phased indel
genotype store** (`Dictionaries/indel_genotypes_<vcf>/`) and (b) **overlays SNP IUPAC
codes** onto the fake-indel genome flanks before indexing the `_INDELS` companion, so
the indel search can match SNP+indel haplotypes. Post-analysis then emits, per
co-occurring off-target, a **CONFIRMED-cis** call (all covered variants phased on one
haplotype) or **PUTATIVE** (unphased / can't prove cis), with the per-sample carriers
and the joint allele frequency (`AC_cis / AN` over the VCF-genotyped panel). On a
sites-only panel it falls back to a PUTATIVE min-AF bound with no carriers (see
METHODS §4/§8).

> **Note.** With `mm` at its maximum, raise `--max-total-edits` (e.g. `--max-total-edits
> 6`) so the extra alignment budget for the indel bulge isn't consumed by mismatches,
> or the co-occurring indel off-targets may be pruned (the indel consumes a bulge slot).


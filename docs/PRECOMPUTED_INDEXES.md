# Precomputed CRISPRme indexes on HuggingFace

> **Reproducing the four shipped production indexes.** The end-to-end recipes for the
> released variant indexes — the single-source phased default **`NRG_3_hg38+hg38_1000G2021`**,
> the genotyped **`NRG_3_hg38+hg38_1000G2021_HGDP`**, the HPRC pangenome
> **`NRG_3_hg38+hg38_HPRC`**, and sites-only **`NRG_3_hg38+hg38_mega`** (5 sources) —
> live in [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md) (Mode 1,
> the HPRC section with `hprc_build.sh`, and Mode 2, with the as-run `mega_gw_merge.sh` +
> `mega_build_indels.sh` drivers). This page covers the generic single-dataset
> build/publish/download flow. **Naming note:** `--dictless` is a *publish flag* (it drops the
> ~152 GB per-sample SNP dictionaries; see below), **not** a name marker — the shipped index
> names are simply `<pam>_<N>_<ref>+<vcf>` (`NRG_3_hg38+hg38_1000G2021`,
> `NRG_3_hg38+hg38_1000G2021_HGDP`, `NRG_3_hg38+hg38_HPRC`, `NRG_3_hg38+hg38_mega`).

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

## The four production indexes (2.5.5)

Four SpCas9 **NRG** (NAG+NGG) variant indexes ship prebuilt on HuggingFace — pick
by whether you want a clean single-source phased default, broader population coverage,
pangenome coverage, or the widest allele-frequency provenance. They are
**complementary**, not alternatives:

| Index name | Sources | Data model | Co-occurrence confidence |
|---|---|---|---|
| `NRG_3_hg38+hg38_1000G2021` *(recommended default)* | 1000 Genomes 2021 (3,202 samples) | **sample-level genotypes, fully phased** | **CONFIRMED** phased cis throughout with **named carrier samples** + exact joint AF (all genotyped-phased, no PUTATIVE-from-phasing) — a clean, simple default of real observed haplotypes |
| `NRG_3_hg38+hg38_1000G2021_HGDP` *(broader coverage)* | 1000 Genomes 2021 + HGDP (adds 929 individuals) | **sample-level genotypes** (hybrid: 1000G phased, HGDP unphased) | 1000G → **CONFIRMED** (phased cis); HGDP → **PUTATIVE** co-carrier (unphased), both with **named carrier samples** + joint AF |
| `NRG_3_hg38+hg38_HPRC` | HPRC Release 2 pangenome `hprc-v2.0-mc-grch38` (232 assembly-derived genomes incl. CHM13) | **sample-level genotypes, phased** (graph-derived) | **CONFIRMED** phased cis with **named carrier samples** + exact joint AF; captures **pangenome-specific** variation absent from short-read panels |
| `NRG_3_hg38+hg38_mega` | 5 sources — 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us | **sites-only** (allele frequencies, no shared samples) | **PUTATIVE** with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`); a co-occurring haplotype here **may or may not exist in any real individual** — no carriers |

**Choosing an index.** You only need **one**. For most use cases pick the **recommended
default `NRG_3_hg38+hg38_1000G2021`** — single-source 1000 Genomes 2021, fully **phased**,
so every reported co-occurrence is **CONFIRMED** cis with **named per-sample carriers** and
exact joint AF (all genotyped-phased, no PUTATIVE-from-phasing): a clean, simple default of
real observed haplotypes. Pick **`NRG_3_hg38+hg38_1000G2021_HGDP`** when you want the extra
HGDP diversity (adds HGDP's 929 individuals for broader population coverage; hybrid
confidence — CONFIRMED on the phased 1000G portion, PUTATIVE co-carrier on the unphased HGDP
portion). Pick **HPRC** for assembly-derived / pangenome variation (phased → CONFIRMED cis +
carriers). **Escalate to `mega`** only for a widest-provenance worst-case screen when you
must not miss a rare allele from any of five databases: it is **sites-only** (union of allele
frequencies, no genotypes), so every multi-variant / co-occurring off-target is **PUTATIVE**
(a worst-case reconstruction that may not exist in any real individual) with a min-AF bound
and **no named carriers**.

All four carry **searchable indels genome-wide** and report **SNP+SNP and SNP+indel**
co-occurring off-targets. Download the reference data once, then any index:

```bash
crisprme.py download --what all  --path .                                                 # reference genome + annotations (once)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_1000G2021      --path .     # recommended default (single-source, phased; carriers)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_1000G2021_HGDP --path .     # broader coverage (adds HGDP; hybrid confidence)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_HPRC           --path .     # HPRC pangenome (phased assemblies, carriers)
crisprme.py download --what index --index-name NRG_3_hg38+hg38_mega          --path .      # mega (sites-only, per-dataset AF provenance)
```

A first-time unauthenticated download of these many-small-file indexes can hit
HuggingFace's HTTP-429 rate limit; the downloader retries with backoff, and
authenticating (`huggingface-cli login`) avoids it entirely. See METHODS §2 for
the genotyped-vs-sites-only distinction and §4 for how each is scanned.

## Layout on HuggingFace

Indexes live under `indexes/` in the CRISPRme dataset repo (default
`lucapinello/crisprme-data`, override with `--hf-repo` or `CRISPRME_HF_REPO`):

```
indexes/
  NRG_3_hg38.tar.gz                              # SpCas9 (NRG = NAG+NGG) reference index of hg38, up-to-2-bulge (DEFAULT reference index)
  NRG_3_hg38+hg38_1000G2021.tar.gz          # SpCas9 (NRG) dict-less variant-aware index (1000G-2021, phased), up-to-2-bulge (recommended default)
  genotypes_hg38_1000G2021.tar.gz                    # SEPARATE Tier-1 genotype store companion for the 1000G2021 index (rides along on download)
  NRG_3_hg38+hg38_1000G2021_HGDP.tar.gz     # SpCas9 (NRG) dict-less variant-aware index (1000G + HGDP), up-to-2-bulge (broader coverage)
  genotypes_hg38_1000G2021_HGDP.tar.gz               # SEPARATE Tier-1 genotype store companion for the 1000G2021+HGDP index (rides along on download)
  NRG_3_hg38+hg38_HPRC.tar.gz                    # SpCas9 (NRG) HPRC pangenome variant index (232 phased assembly-derived genomes)
  genotypes_hg38_HPRC.tar.gz                          # SEPARATE Tier-1 genotype store companion for the HPRC index
  NRG_3_hg38+hg38_mega.tar.gz                    # SpCas9 (NRG) sites-only mega index (5 sources) — NO genotype companion (sites-only)
```

Each variant index stamps its **`data_type`** (`sites-only` / `genotyped-unphased` /
`genotyped-phased` / `hybrid`) and a `phased` flag in every `registry_<vcf>/reg_<chrom>.idx`
manifest, so tooling (e.g. the web variant-dataset selector) can label an index's type
and phasing without scanning the multi-GB per-sample store — `1000G2021` = **genotyped-phased**,
`1000G2021_HGDP` = **hybrid**, `HPRC` = **genotyped-phased**, `mega` = **sites-only**.

The **NRG** default matches SpCas9's broad recognition (NAG + NGG), so variant-created
NAG off-targets (e.g. the CPS1 off-target from the CRISPRme paper) are found out of the
box.

**Two artifacts per variant index:**

- **Main tarball** (`NRG_3_hg38+hg38_1000G2021_HGDP.tar.gz`) — the index
  itself (`<name>/`), its `_INDELS` companion, the indel logs
  (`Dictionaries/log_indels_<vcf>/`), the Tier-0 `Dictionaries/registry_<vcf>/`,
  the samplesID lists (`samplesIDs/<vcf>.samplesID.txt` + the per-db
  `samplesIDs/<ref>_<db>.samplesID.txt`), and `manifest.json` at the archive
  root. In a **classic** (non-dict-less) publish the main tarball ALSO carries
  the per-sample SNP dicts (`Dictionaries/dictionaries_<vcf>/`).
- **Genotype companion** (`genotypes_hg38_1000G2021_HGDP.tar.gz`) — the big Tier-1
  store, uploaded separately under the same `indexes/` prefix. `download --what
  index` fetches it automatically unless `--no-genotypes` is given.

**Naming — `--dictless` is a flag, not a name marker.** `publish-index --dictless`
controls *what goes in the tarball* (it drops the ~152 GB per-sample SNP
dictionaries, keeping the compact Tier-0/Tier-1 stores), but it does **not** decorate
the index name. A variant index is published and installed under the plain convention
`<pam>_<N>_<ref>+<vcf>` — e.g. `NRG_3_hg38+hg38_1000G2021_HGDP`, extracting to
`genome_library/NRG_3_hg38+hg38_1000G2021_HGDP/`. The search resolves an index by that
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
  --vcf VCFs/hg38_1000G2021_HGDP --samplesID samplesIDs.config.txt \
  --path "$CRISPRME_DIR"
# -> genome_library/NRG_3_hg38+hg38_1000G2021_HGDP/          (+ _INDELS companion)
#    Dictionaries/registry_hg38_1000G2021_HGDP/              (Tier-0 allele-freq registry)
#    Dictionaries/registry_hg38_1000G2021_HGDP/variant_count.json  (SNP + indel counts, for the report)
#    Dictionaries/genotypes_hg38_1000G2021_HGDP/             (Tier-1 per-sample genotype store)
#    Dictionaries/log_indels_hg38_1000G2021_HGDP/            (indel logs, for indel post-analysis)
#    samplesIDs/hg38_1000G2021_HGDP.samplesID.txt            (combined, emitted by the build)
```

**One command, all supporting files.** `build-index-only` writes *everything* a
variant-aware search + report needs — the index and its `_INDELS` companion, the
Tier-0 registry, the Tier-1 genotype store, the indel logs, the combined samplesID
list, and the `variant_count.json` manifest (`n_records` SNPs + `n_indels`, which
feeds the report's *Variants included* line). Nothing else has to be assembled by
hand, so `publish-index` (below) ships a self-complete index.

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
bcftools merge 1000G.norm.vcf.gz HGDP.norm.vcf.gz NewCohort.norm.vcf.gz -Oz -o VCFs/hg38_1000G2021_HGDP_NewCohort/merged.vcf.gz
# 3) build the dict-less index (emits registry + genotypes + variant_count.json + combined samplesID)
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --vcf VCFs/hg38_1000G2021_HGDP_NewCohort --samplesID samplesIDs.config.txt --path "$CRISPRME_DIR"
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
crisprme.py publish-index --index genome_library/NRG_3_hg38+hg38_1000G2021_HGDP --dictless
# -> indexes/NRG_3_hg38+hg38_1000G2021_HGDP.tar.gz     (main: index + _INDELS + registry (+ variant_count.json) + indel logs + samplesIDs + manifest)
# -> indexes/genotypes_hg38_1000G2021_HGDP.tar.gz      (separate Tier-1 companion)
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
  --index-name NRG_3_hg38+hg38_1000G2021_HGDP --path "$CRISPRME_DIR"
# installs under the canonical name genome_library/NRG_3_hg38+hg38_1000G2021_HGDP/
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
- The **CLI search-list files** `list_vcf.txt` (the dataset, `hg38_1000G2021_HGDP`)
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
  "name": "NRG_3_hg38+hg38_1000G2021_HGDP",
  "created_at": "2026-08-05T12:00:00+00:00",
  "pam": "NRG",
  "index_bmax": "3",
  "genome": "HGDP",
  "display_label": "SpCas9 NRG — hg38 (1000G + HGDP)",
  "has_registry": true,
  "dictless": true,
  "has_genotypes": true,
  "has_samplesids": true
}
```

> **Note:** `pam` / `index_bmax` / `genome` are parsed from the index *name* and are
> only meaningful for a **reference** index. For a **variant/merged** index the name
> is `<pam>_<N>_<ref>+<vcf>`, so this name-splitting makes `genome` the last VCF token
> (here `HGDP`), not the reference assembly — a harmless artifact. The authoritative
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
- Variant-enriched (genome + VCF) indexes are dataset-specific. The 1000G + HGDP
  enriched index is prebuilt and hosted alongside the reference index; for any
  other cohort, build and publish the enriched index the same way (pass
  `--samplesID` so it is self-complete).
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
    --bDNA 2 --bRNA 2 --vcf VCFs/hg38_1000G2021_HGDP --samplesID samplesID.listing.txt --path ./
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


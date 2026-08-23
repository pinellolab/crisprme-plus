# Precomputed CRISPRme indexes on HuggingFace

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

## Layout on HuggingFace

Indexes live under `indexes/` in the CRISPRme dataset repo (default
`lucapinello/crisprme-data`, override with `--hf-repo` or `CRISPRME_HF_REPO`):

```
indexes/
  NRG_3_hg38.tar.gz                              # SpCas9 (NRG = NAG+NGG) reference index of hg38, up-to-2-bulge (DEFAULT)
  NRG_3_hg38-dictless+hg38_1000G_HGDP.tar.gz     # SpCas9 (NRG) dict-less variant-aware index (1000G + HGDP), up-to-2-bulge (web default)
  genotypes_hg38_1000G_HGDP.tar.gz               # SEPARATE Tier-1 genotype store companion for the variant index (rides along on download)
```

The **NRG** default matches SpCas9's broad recognition (NAG + NGG), so variant-created
NAG off-targets (e.g. the CPS1 off-target from the CRISPRme paper) are found out of the
box.

**Two artifacts per variant index:**

- **Main tarball** (`NRG_3_hg38-dictless+hg38_1000G_HGDP.tar.gz`) — the index
  itself (`<name>/`), its `_INDELS` companion, the indel logs
  (`Dictionaries/log_indels_<vcf>/`), the Tier-0 `Dictionaries/registry_<vcf>/`,
  the samplesID lists (`samplesIDs/<vcf>.samplesID.txt` + the per-db
  `samplesIDs/<ref>_<db>.samplesID.txt`), and `manifest.json` at the archive
  root. In a **classic** (non-dict-less) publish the main tarball ALSO carries
  the per-sample SNP dicts (`Dictionaries/dictionaries_<vcf>/`).
- **Genotype companion** (`genotypes_hg38_1000G_HGDP.tar.gz`) — the big Tier-1
  store, uploaded separately under the same `indexes/` prefix. `download --what
  index` fetches it automatically unless `--no-genotypes` is given.

**The `-dictless` marker.** A dict-less variant index is published with a
`-dictless` marker in the REF segment of its name
(`NRG_3_hg38-dictless+hg38_1000G_HGDP`), so it extracts as
`genome_library/NRG_3_hg38-dictless+hg38_1000G_HGDP/`. The search resolves an
index by the convention `<pam>_<N>_<ref>+<vcf>` (ref segment == genome-folder
basename), so **download strips the marker** and installs under the canonical
name `NRG_3_hg38+hg38_1000G_HGDP` (the `+<vcf>` segment — shared with the
`genotypes_<vcf>` companion — is preserved verbatim). A reference-only index
(no `+`) unpacks to a single `genome_library/<name>/` directory plus
`manifest.json` and is used with no extra steps.

## Dict-less flow (variant-aware index)

### Build (maintainers)

Build the variant-aware index with `--vcf` and `--samplesID`:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 \
  --vcf VCFs/hg38_1000G_HGDP --samplesID samplesIDs.config.txt \
  --path "$CRISPRME_DIR"
# -> genome_library/NRG_3_hg38+hg38_1000G_HGDP/          (+ _INDELS companion)
#    Dictionaries/registry_hg38_1000G_HGDP/              (Tier-0 allele-freq registry)
#    Dictionaries/registry_hg38_1000G_HGDP/variant_count.json  (SNP + indel counts, for the report)
#    Dictionaries/genotypes_hg38_1000G_HGDP/             (Tier-1 per-sample genotype store)
#    Dictionaries/log_indels_hg38_1000G_HGDP/            (indel logs, for indel post-analysis)
#    samplesIDs/hg38_1000G_HGDP.samplesID.txt            (combined, emitted by the build)
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
bcftools merge 1000G.norm.vcf.gz HGDP.norm.vcf.gz NewCohort.norm.vcf.gz -Oz -o VCFs/hg38_1000G_HGDP_NewCohort/merged.vcf.gz
# 3) build the dict-less index (emits registry + genotypes + variant_count.json + combined samplesID)
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --vcf VCFs/hg38_1000G_HGDP_NewCohort --samplesID samplesIDs.config.txt --path "$CRISPRME_DIR"
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
crisprme.py publish-index --index genome_library/NRG_3_hg38+hg38_1000G_HGDP --dictless
# -> indexes/NRG_3_hg38+hg38_1000G_HGDP.tar.gz     (main: index + _INDELS + registry (+ variant_count.json) + indel logs + samplesIDs + manifest)
# -> indexes/genotypes_hg38_1000G_HGDP.tar.gz      (separate Tier-1 companion)
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
  --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path "$CRISPRME_DIR"
# installs under the canonical name genome_library/NRG_3_hg38+hg38_1000G_HGDP/
```

- The index installs under its **canonical** name (the `-dictless` marker is
  stripped from the ref segment) so the search resolves it.
- The **combined + per-db samplesID lists** are installed into
  `samplesIDs/` from the main tarball, so no `--what samples` / `--what all` is
  needed. (For a legacy index built before this bundling existed, download falls
  back to fetching the per-db lists from HF and synthesizing the combined file.)
- The **Tier-1 genotype store** is fetched automatically. Add `--no-genotypes`
  to skip the big companion: off-target **detection still works** via the Tier-0
  registry, but per-sample `Samples` are degraded until the store is present.
- The **CLI search-list files** `list_vcf.txt` (the dataset, `hg38_1000G_HGDP`)
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
  "name": "NRG_3_hg38-dictless+hg38_1000G_HGDP",
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
```


# Precomputed CRISPRme indexes on HuggingFace

Bulge-enabled CRISPRme searches need a CRISPRitz **index** of the reference
genome. Building it is the single most expensive one-time step of a search. That
index depends only on the **genome + PAM + bulge count** (not on the guides or
the variant dataset), so it can be built once and reused — or built once by a
maintainer, published to HuggingFace, and downloaded by everyone else.

This document covers the download → build → publish workflow. It complements
Section 3.5 of the data-setup guide (`docs/crisprme_data_setup_051826.md`).

## Layout on HuggingFace

Indexes live under `indexes/` in the CRISPRme dataset repo (default
`lucapinello/crisprme-data`, override with `--hf-repo` or `CRISPRME_HF_REPO`):

```
indexes/
  NRG_3_hg38.tar.gz                     # SpCas9 (NRG = NAG+NGG) reference index of hg38, up-to-2-bulge (DEFAULT)
  NRG_3_hg38+hg38_1000G_HGDP.tar.gz     # SpCas9 (NRG) variant-aware index (1000G + HGDP), up-to-2-bulge (web default)
  NGG_3_hg38.tar.gz                     # SpCas9 (canonical NGG only) reference index of hg38, up-to-2-bulge
  NGG_3_hg38+hg38_1000G_HGDP.tar.gz     # SpCas9 (NGG only) variant-aware index (1000G + HGDP), up-to-2-bulge
  NNN_3_hg38+hg38_1000G_HGDP.tar.gz     # pamless (any-PAM) variant-aware index (1000G + HGDP), up-to-2-bulge (advanced)
```

The **NRG** default matches SpCas9's broad recognition (NAG + NGG), so variant-created
NAG off-targets (e.g. the CPS1 off-target from the CRISPRme paper) are found out of the
box. Use the `NGG_3_*` indexes if you specifically want canonical-NGG-only results.

Each tarball unpacks to a single `genome_library/<name>/` directory plus a
`manifest.json` describing its provenance (PAM, bulge level, genome, build
timestamp). The folder name is `<PAM>_<bMax+1>_<ref>` — exactly what
`complete-search` looks for — so a downloaded index is used with no extra steps.

## Download (end users)

```bash
# fetch a prebuilt index straight into genome_library/
crisprme.py download --what index --index-name NRG_3_hg38 --path "$CRISPRME_DIR"

# then search, pointing at that library (or just run from $CRISPRME_DIR)
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --guide my_guide.txt --mm 4 --bDNA 2 --bRNA 2 \
  --index-path "$CRISPRME_DIR/genome_library" \
  --output my_search
```

The search finds the prebuilt index and skips the build entirely.

## Build (maintainers)

`build-index-only` builds the index without running a search. Pass the **same**
`--genome`/`--pam`/`--bDNA`/`--bRNA` end users will search with, so the folder
name matches:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 --path "$CRISPRME_DIR"
# -> genome_library/NRG_3_hg38/
```

## Publish (maintainers)

Upload the built index to the dataset repo (needs an HF **write** token — via
`--token` or the `HF_TOKEN` env var; never commit it):

```bash
export HF_TOKEN=hf_...        # your write token, in the shell only
crisprme.py publish-index --index genome_library/NRG_3_hg38
# -> uploaded to indexes/NRG_3_hg38.tar.gz (with a manifest.json inside)
```

## manifest.json

Every published index carries a small provenance manifest inside its tarball:

```json
{
  "name": "NRG_3_hg38",
  "created_at": "2026-08-05T12:00:00+00:00",
  "pam": "NRG",
  "index_bmax": "3",
  "genome": "hg38"
}
```

It is surfaced (build timestamp) when the index is downloaded, and otherwise
ignored — the index directory itself is what `complete-search` consumes.

## Notes

- An index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`;
  a different PAM or a higher bulge count needs its own index.
- Variant-enriched (genome + VCF) indexes are dataset-specific. The 1000G + HGDP
  enriched indexes (`NGG_3_hg38+hg38_1000G_HGDP`, `NNN_3_hg38+hg38_1000G_HGDP`)
  are prebuilt and hosted alongside the reference index; for any other cohort,
  build and publish the enriched index the same way.
- If `--index-path` is given but no matching index is found there,
  `complete-search` fails fast with a clear message rather than silently
  rebuilding — so a missing/wrong download is caught immediately.

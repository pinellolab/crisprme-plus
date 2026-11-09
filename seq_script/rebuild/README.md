# The rebuild of the published indexes for CRISPRme 2.10.0

`index_set_2.10.0.tsv` is the **rebuild manifest**: every archive the rebuild for 2.10.0 builds
and publishes in [`lucapinello/crisprme-indexes-format2`](../../docs/hf/crisprme-indexes-format2_README.md),
the location of index format 2 that no release before 2.10.0 reads (decision C2), one line each,
with its panel, PAM file, PAM, side, spacer (guide) length, bMax (N, the name's), the build's
`--bDNA`/`--bRNA`, the archives of 2.8.1 it replaces, and its build and publish commands.

| Who reads it | For |
|---|---|
| the rebuild | what to build, in which working directory, and how to publish it |
| the D7 readiness audit | `check_index_set.py`, before any build: the set against the code |
| `crisprme.py download` (`crisprme_hf.successor_index`) | which index of 2.10.0 replaces an archive of 2.8.1 a user asks for by name |
| the docs | `docs/PRECOMPUTED_INDEXES.md`'s table and the Hugging Face card list the same set |

## The set (approved 2026-11-08, "the lean set")

- SpCas9 **NRR** 20-nt (the default) on the reference, HGDP + 1000 Genomes, HPRC and mega;
- SpCas9 **NGG** 20-nt, the **general** PAM-less index (an 8-nt 3′ flank, 20-nt guide), SaCas9
  **NNGRRN** 22-nt (Manuel's CEP290 guides), Cas12a **TTTN** 5′ at 20 and 23 nt (every shipped
  Cas12a PAM and Manuel's Cas12a guides), and SpCas9 **NRG** 20-nt (Kris's TRAC reports), each on
  the reference and HGDP + 1000 Genomes;
- **NNN_3 is not built again**: the general index serves its searches at the same speed, and a
  download of `NNN_3_*` says so; 2.8.1's archives stay where they are;
- each variant panel's **Dictionaries bundle** and genotype store once
  (`dictionaries-bundle_<dataset>`, `genotypes_<dataset>`, `indel_genotypes_<dataset>`), and its
  enriched-genome kit (`build-inputs/<dataset>`).

Every index is of N = 3 (`--bDNA 2 --bRNA 2`), the web's ceiling of 2 bulges of each type, which
it serves for its own guide length (`index_fit.fits`).

## Rules the rebuild keeps

1. **One working directory per panel**, every index of the panel built in it (`$HGDP`, `$HPRC`,
   `$MEGA`; `$REF` for the reference indexes). The first build enriches the genome and writes the
   Dictionaries; each later PAM re-indexes the same enriched genome ("Enriched genome already
   present"). So the panel's indexes share one Dictionaries, whose fingerprint every one of their
   manifests records: `publish-index` uploads the bundle and the genotype store with the panel's
   first index, skips them for the others (the repo holds them, the same), and refuses to publish
   an index built against another build of them (`--replace-companions` replaces them, for a
   whole panel built again).
2. **Build with the release image only** (the CRISPRme and CRISPRitz the release ships), after the
   D7 readiness audit: no fix that can change index creation is pending.
3. **Publish each panel as a unit**, then write `indexes/SHA256SUMS` with
   `scripts/hf_sha256sums.py lucapinello/crisprme-indexes-format2 --staged staged.sha256` and
   check it against the Hub after the upload (`docs/hf/README.md`).

## Commands

```bash
python seq_script/rebuild/check_index_set.py              # the check (the D7 audit runs it): exit 1 on a problem
python seq_script/rebuild/check_index_set.py --commands   # every build command, then every publish command, in order
```

The check holds when: every PAM file exists and gives the row's PAM, side and spacer; N =
max(bDNA, bRNA) + 1; every index's name is the one `build-index-only` gives its PAM file
(`index_fit.index_name`, guide length and side); every archive routes to the format-2 repo
(`crisprme_hf.resolve_index_repo`); every index serves its own guide at 2 + 2 bulges; every
variant index's reference index is in the set; every variant panel has one bundle and one kit,
and a genotyped one (HGDP + 1000 Genomes, HPRC) one genotype store and one store of phased indel
genotypes; and every index archive of 2.8.1 is replaced by exactly one archive of the set, the
one a download of it names. `PostProcess/test_index_set.py` runs it in the unit tests.

The mega rows follow lane M's full rebuild of `hg38_mega` (decision D2: PASS/'.' records, AF =
AF_max, its registry and indel allele frequencies consistent), whose recipe
(`seq_script/merge_panels/mega_gw_merge.sh`, `mega_build_indels.sh`) that lane owns.

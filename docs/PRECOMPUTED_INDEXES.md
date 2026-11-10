# Precomputed CRISPRme indexes on HuggingFace

> **Reproducing the shipped production indexes.** The end-to-end recipes for the
> released variant indexes — the phased default **`NRR_3_hg38+hg38_HGDP1kGP`** (HGDP +
> 1000 Genomes; NRR is the default PAM since 2.10.0, and the NRG and NGG indexes are built the
> same way), the HPRC pangenome
> **`NRR_3_hg38+hg38_HPRC`**, and sites-only **`NRR_3_hg38+hg38_mega`** (5 sources) —
> live in [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md) (Mode 1,
> the HPRC section with `hprc_build.sh`, and Mode 2, with the as-run `mega_gw_merge.sh` +
> `mega_build_indels.sh` drivers). This page covers the generic single-dataset
> build/publish/download flow. **Naming note:** `--dictless` is a *publish flag* (it drops the
> ~152 GB per-sample SNP dictionaries; see below), **not** a name marker — the shipped index
> names are simply `<pam>_<N>_<ref>+<vcf>` (`NRR_3_hg38+hg38_HGDP1kGP`,
> `NRG_3_hg38+hg38_HPRC`, `NRG_3_hg38+hg38_mega`). Since 2.10.0 an index is built under a name
> that carries the guide length and the side of the PAM too, in the order of the target on the
> strand (`20-NRR_3_hg38` for a 3′ PAM, `TTTV-23_3_hg38` for a 5′ one), and published under it;
> CRISPRme reads both forms, takes one of the new form only for a search of its side, and of
> one index under both names, the newer one.

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

Throughout this page `<vcf>` is the variant-dataset token resolved from the index name's
`+<vcf>` segment — for the recommended default `NRR_3_hg38+hg38_HGDP1kGP` it is
**`hg38_HGDP1kGP`**, so the folders are literally `registry_hg38_HGDP1kGP/`,
`genotypes_hg38_HGDP1kGP/`, and `hg38_HGDP1kGP.samplesID.txt`. Substitute your own dataset's
token where you see `<vcf>`, and note that the resolver expects these exact names.

This document covers the download → build → publish workflow. It complements
Section 3.5 of the data-setup guide (`docs/crisprme_data_setup_051826.md`).

## The three production indexes

SpCas9 **NRR** (NAA, NAG, NGA, NGG; the default PAM since 2.10.0) variant indexes ship
prebuilt on HuggingFace, as do the same panels for **NRG** (NAG + NGG, the default before) —
pick by whether you want the broad fully-phased default, pangenome coverage, or the widest
allele-frequency provenance. They are **complementary**, not alternatives:

| Index name | Sources | Data model | Co-occurrence confidence |
|---|---|---|---|
| `NRR_3_hg38+hg38_HGDP1kGP` *(recommended default)* | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000 Genomes), jointly SHAPEIT5-re-phased | **sample-level genotypes, fully phased** | **CONFIRMED** phased cis throughout with **named carrier samples** + exact joint AF (including HGDP's seven genetic regions, with correct chrX PAR allele frequencies) — the broadest fully-phased panel of real observed haplotypes |
| `NRR_3_hg38+hg38_HPRC` | HPRC Release 2 pangenome `hprc-v2.0-mc-grch38` (232 assembly-derived genomes incl. CHM13) | **sample-level genotypes, phased** (graph-derived) | **CONFIRMED** phased cis with **named carrier samples** + exact joint AF; captures **pangenome-specific** variation absent from short-read panels |
| `NRR_3_hg38+hg38_mega` | 5 sources — 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us | **sites-only** (allele frequencies, no shared samples) | **PUTATIVE** with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`); a co-occurring haplotype here **may or may not exist in any real individual** — no carriers |

**Choosing an index.** You only need **one**. For most use cases pick the **recommended
default `NRR_3_hg38+hg38_HGDP1kGP`** — HGDP + 1000 Genomes (4,091 samples), jointly
SHAPEIT5-re-phased, so every reported co-occurrence is **CONFIRMED** cis with **named
per-sample carriers** and exact joint AF throughout, including HGDP's seven genetic regions
and with correct chrX PAR allele frequencies: the broadest fully-phased panel of real
observed haplotypes. Pick
**HPRC** for assembly-derived / pangenome variation (phased → CONFIRMED cis +
carriers). **Escalate to `mega`** only for a widest-provenance worst-case screen when you
must not miss a rare allele from any of five databases: it is **sites-only** (union of allele
frequencies, no genotypes), so every multi-variant / co-occurring off-target is **PUTATIVE**
(a worst-case reconstruction that may not exist in any real individual) with a min-AF bound
and **no named carriers**.

All three carry **searchable indels genome-wide** and report **SNP+SNP and SNP+indel**
co-occurring off-targets. Download the reference data once, then any index:

```bash
crisprme.py download --what all  --path .                                                 # reference genome + annotations (once)
crisprme.py download --what index --index-name NRR_3_hg38+hg38_HGDP1kGP       --path .     # recommended default (phased HGDP + 1000 Genomes, 4,091 samples; carriers)
crisprme.py download --what index --index-name NRR_3_hg38+hg38_HPRC           --path .     # HPRC pangenome (phased assemblies, carriers)
crisprme.py download --what index --index-name NRR_3_hg38+hg38_mega          --path .      # mega (sites-only, per-dataset AF provenance)
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
  NRG_3_hg38.tar.gz                              # SpCas9 (NRG = NAG+NGG) reference index of hg38, up-to-2-bulge (the default reference index before 2.10.0)
  NRG_3_hg38+hg38_HGDP1kGP.tar.gz           # SpCas9 (NRG) dict-less variant-aware index (HGDP + 1000 Genomes, 4,091 samples, phased), up-to-2-bulge
  genotypes_hg38_HGDP1kGP.tar.gz                     # SEPARATE Tier-1 genotype store companion for the HGDP1kGP index (rides along on download)
  NRG_3_hg38+hg38_HPRC.tar.gz                    # SpCas9 (NRG) HPRC pangenome variant index (232 phased assembly-derived genomes)
  genotypes_hg38_HPRC.tar.gz                          # SEPARATE Tier-1 genotype store companion for the HPRC index
  NRG_3_hg38+hg38_mega.tar.gz                    # SpCas9 (NRG) sites-only mega index (5 sources) — NO genotype companion (sites-only)
  NNN_3_hg38.tar.gz                              # PAM-LESS reference index of hg38, 3-nt flank — serves ANY 3-nt 3′ PAM with `--prune-pam`
  NNN_3_hg38+hg38_HGDP1kGP.tar.gz           # PAM-LESS variant-aware index (HGDP + 1000 Genomes, same variant data as the NRG one)
  NGG_3_hg38.tar.gz                              # NGG-only reference index (NGG cannot serve NAG or NGA off-targets)
  NGG_3_hg38+hg38_HGDP1kGP.tar.gz           # NGG-only variant index
  NRR_3_hg38.tar.gz                              # SpCas9 (NRR = NAA, NAG, NGA, NGG) reference index (DEFAULT reference index since 2.10.0)
  NRR_3_hg38+hg38_{HGDP1kGP,HPRC,mega}.tar.gz    # NRR variant indexes (HGDP1kGP: the recommended default)
  SHA256SUMS                                     # one SHA-256 per archive: `sha256sum -c --ignore-missing SHA256SUMS` inside indexes/
```

The two 8-nt-flank PAM-less indexes are in a second repo,
`lucapinello/crisprme-data-general-index` (see the flank-width note below), under the same
`indexes/` layout:

```
indexes/
  NNNNNNNN_3_hg38.tar.gz                         # PAM-LESS reference index of hg38, 8-nt flank
  NNNNNNNN_3_hg38+hg38_HGDP1kGP.tar.gz           # PAM-LESS variant-aware index (HGDP + 1000 Genomes)
  genotypes_hg38_HGDP1kGP.tar.gz                 # the same genotype companion as in the main repo (same SHA-256)
  SHA256SUMS
```

Download sizes, as the Hugging Face API reports them (GB = 10⁹ bytes; the whole default repo
is about 432 GB, so fetch single archives rather than the repo):

| Archive | Repo | Download |
|---|---|---|
| `NRG_3_hg38` | default | 8.9 GB |
| `NRG_3_hg38+hg38_HGDP1kGP` + `genotypes_hg38_HGDP1kGP` | default | 27.7 GB + 5.7 GB |
| `NRG_3_hg38+hg38_HPRC` + `genotypes_hg38_HPRC` | default | 15.2 GB + 0.6 GB |
| `NRG_3_hg38+hg38_mega` | default | 12.1 GB |
| `NGG_3_hg38`, `NGG_3_hg38+hg38_HGDP1kGP` | default | 3.7 GB, 21.5 GB |
| `NRR_3_hg38`, `+hg38_HGDP1kGP`, `+hg38_HPRC`, `+hg38_mega` | default | 20.1, 42.9, 30.7, 25.2 GB |
| `NNN_3_hg38`, `NNN_3_hg38+hg38_HGDP1kGP` | default | 71.5 GB, 109.5 GB |
| `NNNNNNNN_3_hg38`, `NNNNNNNN_3_hg38+hg38_HGDP1kGP` | `crisprme-data-general-index` | 82.2 GB, 122.7 GB |

**The PAM-less (`NNN_*`) indexes** store the raw genomic flank instead of a matched PAM,
so one index answers any 3′ PAM no longer than that flank when the search passes
`--prune-pam` (see [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)).
A 5′ PAM (Cas12a) needs a dedicated index: these indexes store the flank on the 3′ side of
the guide, so a 5′ search does not use them.
Every index here was built for a 20-nt guide, with N = 3. That also serves a 21-nt guide at up
to 1 bulge of each type, byte-identically to an index built for it (so the shipped 21-nt SaCas9
PAM files work on the 8-nt-flank indexes with `--prune-pam`); a longer guide needs its own
index ([METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)).
They carry exactly the same variant data and `Dictionaries/` as the dedicated index for the
same panel — those tiers are PAM-independent — and a `--prune-pam` run over one gives the
dedicated index's report for that PAM: byte for byte from a flank of the PAM's width (`NNN_3`
for NRG, verified genome-wide), and from a wider one (K = 8 for a 3-nt PAM) but for what it
cannot index at a contig's very end and start: a + strand PAM with 1 to K − L nt after it to the
end, and a − strand PAM in the first K − L nt. On hg38 at K = 8 that is 145 + 101 NGG and
308 + 241 NRG sites, about 7.5 × 10⁻⁷ of all, none on a primary chromosome (whole-genome, K = 8
lost one row of 27,216 for a guide picked for its perfect match there, and one of 216,978 for
VEGFA site 3); a search that reads a wider flank of format 1 says so on stdout (an index of
format 2 has no such limit: with `--prune-pam` it gives a dedicated index's rows for the same
guide length, contig ends included). CRISPRitz 2.9.1's search also reported a − strand PAM
exactly K − L nt from the start on the + strand at another position; 2.10.0's reads it right in
either format. A wider index must have been built by
CRISPRme 2.8.1 (CRISPRitz 2.9.1), which keeps a target whose flank holds an `N`, as every
archive published here was; the search warns of one built before. They are
correspondingly larger, because nothing is excluded at build time: 71 GB for the reference
index against 9 GB for `NRG_3_hg38`, and 110 GB against 28 GB for the HGDP + 1000 Genomes
variant index. Choose one when you expect to search several PAMs over the same genome.

**Limits of every index of format 1 at contig ends and `N` runs** (dedicated or PAM-less, any
flank; [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)). Every index
published here was built by CRISPRitz 2.9.1 or earlier, in index format 1. CRISPRitz 2.10.0,
which CRISPRme 2.10.0 installs, writes format 2, which holds these sites, and reads a − strand PAM
at a contig's first base right in either format; CRISPRme 2.10.0 reads both formats, and says in
the search's `log_verbose.txt` when it takes an index of format 1 (2.8.1 refuses an index of
format 2, and 2.8.0 aborts on one). An index of format 1 of
budget N cannot report a site whose guide lies within N − 1 nt of a contig start, a contig end or
an `N` run, even a perfect match, which a search without bulges (a scan of the FASTA) does
report; but a variant search without bulges on an install without its whole enriched genome (a
download) searches the variant index, so its SNP variant rows have this limit too. So an index of
a larger budget is not a strict superset of a smaller one. For `NRG_3_hg38` that is 546 + strand
and 504 − strand sites (1,050 NRG sites in all, 535 on primary contigs), none within 4 mismatches
of nine therapeutic guides checked: what matters is a guide whose own on-target sits at such a
boundary, whose report then has no perfect match. Under CRISPRitz 2.9.1 (CRISPRme 2.8.1),
a − strand PAM at a contig's first base is reported on the + strand (11 NGG and 37 NRG hg38
contigs), and a + strand PAM on a contig's last bases is found by no search, indexed or not (8
NGG and 31 NRG contigs): all alt, random or chrUn contigs, none a primary chromosome. Under
CRISPRitz 2.10.0 the first is reported on its own strand from an index of either format, and the
second by a search without bulges (a scan of the FASTA) and from an index of format 2, but not
from one of format 1.

> **Flank width and compatibility.** A 3-nt flank is readable by any CRISPRme version: the
> stored geometry is the same as a dedicated 3-nt PAM index, which is why `NNN_3` ships in
> this repo. (Use the rebuilt `NNN_3` archives, published on 2026-10-30, with 2.8.1 or later:
> 2.8.0 keeps the targets with an `N` in their PAM that the rebuild adds, and a search that
> reaches one fails or, with RNA bulges only on CRISPRitz 2.9.0, reports it as a malformed
> 0-mismatch, 0-bulge row; see [`INDEX_PROVENANCE.md`](INDEX_PROVENANCE.md).) A **wider** flank is not — a release predating `--prune-pam` would mis-read it
> — so those indexes live in a separate dataset repo,
> **`lucapinello/crisprme-data-general-index`**, which older defaults never point at.
> `download` and `publish-index` route by the flank width in the index name, so a wide-flank
> index cannot be published to the default repo by accident, and the Settings catalog merges
> both repos so what is listed is what can be fetched. An explicit `--hf-repo` or
> `CRISPRME_HF_REPO` overrides the routing.

> Per-archive provenance and validity, with the measurement behind each claim, is in
> [`INDEX_PROVENANCE.md`](INDEX_PROVENANCE.md).

Each variant index stamps its **`data_type`** (`sites-only` / `genotyped-unphased` /
`genotyped-phased` / `hybrid`) and a `phased` flag in every `registry_<vcf>/reg_<chrom>.idx`
manifest, so tooling (e.g. the web variant-dataset selector) can label an index's type
and phasing without scanning the multi-GB per-sample store — `HGDP1kGP` = **genotyped-phased**,
`HPRC` = **genotyped-phased**, `mega` = **sites-only**
(a merged panel mixing phased and unphased sources would be **hybrid**).

The **NRR** default (since 2.10.0; NRG before) matches SpCas9's broad recognition (NAA, NAG,
NGA and NGG), so variant-created NAG off-targets (e.g. the CPS1 off-target from the CRISPRme
paper) and NGA ones are found out of the box. An NRR index serves an NGG, NAG, NRG or NGA
search too, pruned to that PAM (`--prune-pam`, which the search passes by itself), and an
NRG one an NGG or NAG search.

**Two artifacts per variant index:**

- **Main tarball** (`NRR_3_hg38+hg38_HGDP1kGP.tar.gz`) — the index
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
the index name. A variant index is published and installed under the name of the folder it
was built in: the plain convention `<pam>_<N>_<ref>+<vcf>` for every index published so far
— e.g. `NRG_3_hg38+hg38_HGDP1kGP`, extracting to `genome_library/NRG_3_hg38+hg38_HGDP1kGP/`
— and, for an index built since 2.10.0, that name with its guide length and side, in the order
of the target on the strand: `<guide>-<pam>_<N>_<ref>+<vcf>` for a 3′ PAM
(`20-NRR_3_hg38+hg38_HGDP1kGP`), `<pam>-<guide>_<N>_<ref>+<vcf>` for a 5′ one. The search
resolves an index by either convention (the `<ref>` segment == the genome-folder basename;
the `<vcf>` token after the `+` — e.g. `hg38_HGDP1kGP` — is shared verbatim with the Tier-1
`genotypes_<vcf>` companion, i.e. `genotypes_hg38_HGDP1kGP`). A reference-only index (no `+`)
unpacks to a single `genome_library/<name>/` directory plus `manifest.json` and is used
with no extra steps.

## Dict-less flow (variant-aware index)

### Build (maintainers)

Build the variant-aware index with `--vcf` and `--samplesID`:

```bash
crisprme.py build-index-only \
  --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 \
  --vcf VCFs/hg38_HGDP1kGP --samplesID samplesIDs.config.txt \
  --path "$CRISPRME_DIR"
# -> genome_library/20-NRR_3_hg38+hg38_HGDP1kGP/       (+ _INDELS companion; the name says the guide length and the side)
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
`samplesIDs/`). For `VCFs/hg38_HGDP1kGP` it names that panel's one samples file (the one
`download --what samples` installs):

```
hg38_HGDP1kGP.samplesID.txt
```

A panel merged from several datasets lists each dataset's file instead (next section).

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
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
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
crisprme.py publish-index --index genome_library/20-NRR_3_hg38+hg38_HGDP1kGP --dictless
# -> indexes/20-NRR_3_hg38+hg38_HGDP1kGP.tar.gz  (main: index + _INDELS + registry (+ variant_count.json) + indel logs + samplesIDs + manifest)
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
  --index-name NRR_3_hg38+hg38_HGDP1kGP --path "$CRISPRME_DIR"
# installs under the canonical name genome_library/NRR_3_hg38+hg38_HGDP1kGP/
```

- The index installs under the name it was published under, which the search resolves:
  here the plain `<pam>_<N>_<ref>+<vcf>` of every index published so far (an index
  published from a 2.10.0 build carries its guide length and side, as
  `20-NRR_3_hg38+hg38_HGDP1kGP`). There is no name decoration to strip — `--dictless` only
  affected what the maintainer put in the tarball.
- The **combined + per-db samplesID lists** are installed into
  `samplesIDs/` from the main tarball, so no `--what samples` / `--what all` is
  needed. (For a legacy index built before this bundling existed, download falls
  back to fetching the per-db lists from HF and synthesizing the combined file.)
- The **Tier-1 genotype store** is fetched automatically. Add `--no-genotypes`
  to skip the big companion: a search then runs the **population-level** analysis
  from the Tier-0 registry (each site's worst case, with no CONFIRMED cis phasing
  and no carriers but those of an indel, which the indel logs of the main tarball
  name; where a site needs a SNP too, its row lists every carrier of the indel)
  until the store is present, and `--per-sample` falls back to it with a warning:
  a per-sample analysis without genotypes would lose the off-targets that need two
  SNPs of one sample together.
- The **CLI search-list files** `list_vcf.txt` (the dataset, `hg38_HGDP1kGP`)
  and `list_samplesID.txt` (its combined samplesID) are written at the install
  root, so a CLI search works out of the box — the same lists the web form builds
  per-search. Installing several variant indexes appends each dataset once.

Then search, pointing at that library (or just run from `$CRISPRME_DIR`):

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
  --guide my_guide.txt --mm 4 --bDNA 2 --bRNA 2 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --index-path "$CRISPRME_DIR/genome_library" \
  --output my_search
```

The search finds the prebuilt index (and its tiers + samplesIDs) and skips the
build/enrichment entirely. The variant index and its `_INDELS` index are read from the
`--index-path` library first, then from the working directory's `genome_library/`; the
tiers (`Dictionaries/`), the sample lists (`samplesIDs/`) and the genome (`Genomes/`) are
read from the working directory, so a library kept elsewhere needs those there (links
are fine). A search builds no index in an `--index-path` library: a variant index found
there without its `_INDELS` index beside it stops the search with a message.

## Building a new-PAM index for a published dataset (HGDP+1kGP, mega)

The prebuilt indexes above cover the **NRR**, **NRG** and **NGG** PAMs (not each on every
panel) plus the PAM-less `NNN_3` / `NNNNNNNN_3` indexes, which serve other PAMs through
`--prune-pam`, as an NRR index serves the PAMs it holds. To screen a **different nuclease / PAM** with a dedicated index (SaCas9
**NNGRRT**, Cas12a **TTTN**, any IUPAC motif, or NGG on mega, the one of these two panels it
is not published for) against the same re-phased **HGDP+1kGP** or sites-only **mega** variant data, you build
that one index once and reuse it for every later search — you do **not** re-derive the
variant data.

### What's reusable vs. what you rebuild

Only the genome index itself is PAM-specific; everything derived from the variants is not:

| Artifact | PAM-specific? | For a new PAM |
|---|---|---|
| `genome_library/<PAM>_<N>_hg38+hg38_<ds>/` (+ `_INDELS`) | **yes** | **rebuilt** |
| Enriched genome `Genomes/hg38+hg38_<ds>/` (+ `_INDELS`) | no | reuse (skip enrichment) |
| `Dictionaries/` (registry, genotypes, indel-GT, logs) | no | reuse (from the normal index download) |
| `samplesIDs/` | no | reuse |

> **Why you need an extra file.** The HuggingFace **index** tarball ships the index of its PAM + its
> `Dictionaries/` + `samplesIDs/` — it does **not** contain the enriched genome, so it alone can't
> build a different PAM. The PAM-independent **enriched-genome kit** is published separately under
> **`build-inputs/<dataset>.tar.gz`** in the same repo (`lucapinello/crisprme-data`). Download it and you
> skip straight to indexing the new PAM. If instead you want to **reproduce the kit** (or build one
> for a **different VCF**), the full source-to-kit workflow is documented below.
>
> A search of the index itself needs no enriched genome, with no bulges either: a 0-bulge variant
> search scans the enriched genome when every chromosome of it is in place, and otherwise searches
> the installed variant index, with the same results but for a site whose guide lies within 2 nt
> of a contig start or end or an `N` run (N − 1 for an index of budget N): only the scan finds
> its variant version (see the limits of every index at contig ends and `N` runs, above).

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
| SpCas9 (NRR = NAA, NAG, NGA, NGG; the default) | `20bp-NRR-SpCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNRR 3` |
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
NRR_3_hg38+hg38_<ds>`).

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
> same two indexes, which it names, since 2.10.0, for the guide length and the side
> (`20-<PAM_MOTIF>_3_hg38+hg38_<ds>` and its `_INDELS` for a 20-nt guide and a 3′ PAM). Use it
> if you also have the dataset's `Dictionaries/`/`VCFs/` staged; otherwise the two
> `crispritz`/`pool_index_indels` calls above need only the enriched kit.

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
  --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
  --bDNA 2 --bRNA 2 --thread 16 --path "$CRISPRME_DIR"
# -> genome_library/20-NRR_3_hg38/  (since 2.10.0 the name says the guide length and the side)

crisprme.py publish-index --index genome_library/20-NRR_3_hg38
# -> indexes/20-NRR_3_hg38.tar.gz

crisprme.py download --what index --index-name 20-NRR_3_hg38 --path "$CRISPRME_DIR"
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
- `indel_window` — (a variant index published by CRISPRme 2.10.0 or later) the reference
  bases each fake-indel window holds before an indel and after it, as the index's indel
  logs record them: `[30, 30]` for one CRISPRitz 2.10.0 built (an archive without the field,
  every one published up to 2.8.1, holds 25 and 27). `download --what index` says when an
  archive's windows are not the installed CRISPRitz's ([METHODS §4](../METHODS.md#snpindel-off-targets-each-carriers-own-haplotype)).

The manifest is surfaced (build timestamp / display name) when the index is
downloaded, and old fields are ignored by older consumers — the index directory
itself is what `complete-search` consumes.

## Notes

- An index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`;
  a different PAM or a higher bulge count needs its own index.
- Variant-enriched (genome + VCF) indexes are dataset-specific. The default
  HGDP + 1000 Genomes enriched index (`NRR_3_hg38+hg38_HGDP1kGP`) is prebuilt and
  hosted alongside the reference index; for any other cohort, build and publish the
  enriched index the same way (pass `--samplesID` so it is self-complete).
- A published variant index is **self-complete**: the samplesID lists are
  bundled, so `download --what index` + search works without a separate
  `--what samples`. The genotype store is the only optional piece (skip it with
  `--no-genotypes` for a population-level install, with no carriers but those of an
  indel).
- If `--index-path` is given but no matching reference index is found there,
  `complete-search` fails fast with a clear message rather than silently
  rebuilding — so a missing/wrong download is caught immediately. A variant index
  is looked for there first, then in the working directory's `genome_library/`.

## SNP+indel co-occurring off-targets (on by default)

CRISPRme+ searches **SNP+indel co-occurring off-targets** — those that require **both**
a nearby SNP **and** an indel in the same protospacer window — **by default** since
2.5.0. The prebuilt indexes are already built with this on, so a normal `download` +
`complete-search` reports these co-occurrences **with no extra flag** (a class the
classic two-independent-passes search could not see: SNPs on the IUPAC-enriched genome,
indels on a fake-indel genome cut from the *plain* reference).

To reproduce the classic two-independent-passes behavior, set `CRISPRME_INDEL_SNP=0`
**before both the build and the search**. With the gate off, builds are byte-identical to
classic CRISPRme, and searches run its two independent passes with the fixes of 2.8.1 that
apply to every search, so their results can differ from an earlier release's. Three of those
fixes change indel rows (METHODS §4): a target that holds no base of its indel gives no row;
an indel row with an RNA bulge starts at its first aligned base and shows as REF the reference
its target covers from there, gapped as its target; and of the alignments of one site that name
the same variants with the same carriers, the merge keeps the best by the order that picks each
block's representative, not the first it reads (the `Other_motifs` counts of their clusters
follow).

```bash
# default (co-occurrence ON) — nothing to set:
crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRG-SpCas9.txt \
    --bDNA 2 --bRNA 2 --vcf VCFs/hg38_HGDP1kGP --samplesID samplesID.listing.txt --path ./
# to opt out and reproduce classic CRISPRme, set the gate to 0 first:
#   export CRISPRME_INDEL_SNP=0   # before BOTH build and search
```

With co-occurrence on, the build (a) compiles a per-chromosome **phased indel
genotype store** (`Dictionaries/indel_genotypes_<vcf>/`) and (b) **overlays SNP IUPAC
codes** onto the fake-indel genome flanks before indexing the `_INDELS` companion, so
the indel search can match SNP+indel haplotypes. Each indel's fake-indel window holds 30
reference bases on either side of it in an index CRISPRitz 2.10.0 built (25 before and 27
after in one of 2.9.1, every archive published up to CRISPRme 2.8.1), so a target of guide,
PAM and DNA bulges up to 31 bases (27) holds its indel there; the width is read from each
archive's own indel logs, which every INDELS sub-index records it was built from, and a
search refuses an INDELS index installed beside indel logs of other windows (a 2.8.1
archive's index beside the logs of a 2.10.0 one, say: two archives of one panel and one
release share their logs, and are searched as usual), and says when its targets need more
than its windows hold (METHODS §4). Post-analysis resolves each indel
target on its carriers' own haplotypes (each gets the SNP alleles it carries, so a row
lists only the individuals that have its sequence, or, where the phase is unknown, those
whose own alleles allow it, PUTATIVE, split multi-allelic SNP sites aside on an index
CRISPRitz 2.9.1 or earlier built (METHODS §8);
on an install without SNP genotypes, as after `download --what index --no-genotypes`,
every carrier of the indel, PUTATIVE; a sequence no individual has is dropped; since
2.8.1), then emits, per co-occurring
off-target, a **CONFIRMED-cis** call (all covered variants phased on one haplotype) or
**PUTATIVE** (unphased / can't prove cis), with the per-sample carriers and the joint
allele frequency (`AC_cis / AN` over the VCF-genotyped panel). On a sites-only panel it
falls back to a PUTATIVE min-AF bound with no carriers (see METHODS §4/§8). The indel
genotype store keeps only records whose REF and ALT differ in length, so an equal-length
multi-base record (an MNV, such as AC>TG), which the indel search covers too, has its
carriers resolved with the phase unknown, PUTATIVE: the HGDP + 1000 Genomes index has
none, HPRC's chr22 has 9,079 among its 195,813 indel-search records. Since 2.10.0 a
per-sample search also reads the store for the SNP rows: a SNP row does not list a carrier
whose haplotype holding its SNPs also holds an indel that changes the target (METHODS §4).
With the gate off it lists them, as up to 2.8.1.

> **Note.** With `mm` at its maximum, raise `--max-total-edits` (e.g. `--max-total-edits
> 6`) so the extra alignment budget for the indel bulge isn't consumed by mismatches,
> or the co-occurring indel off-targets may be pruned (the indel consumes a bulge slot).


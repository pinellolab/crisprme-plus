# Precomputed CRISPRme indexes on HuggingFace

> **The indexes of CRISPRme 2.10.0.** Since 2.10.0 the precomputed indexes are built by
> CRISPRitz 2.10.0 in **index format 2** and published in their own dataset repo,
> **[`lucapinello/crisprme-indexes-format2`](https://huggingface.co/datasets/lucapinello/crisprme-indexes-format2)**,
> which no earlier release reads (CRISPRme 2.8.1 refuses an index of format 2, and 2.8.0 aborts
> on one). Each is named for its guide length and the side of its PAM, in the order
> of the target on the strand (`20-NRR_3_hg38` for a 20-nt guide and a 3′ PAM, `TTTN-23_3_hg38`
> for a 23-nt guide and a 5′ PAM, `+<dataset>` for a variant index), and each variant panel's
> Dictionaries and genotype store are published **once**, beside its indexes, and fetched once
> per panel. The set, with how each archive is built and published, is the **index set**,
> [`seq_script/rebuild/index_set_2.10.0.tsv`](../seq_script/rebuild/index_set_2.10.0.tsv)
> ([`seq_script/rebuild/README.md`](../seq_script/rebuild/README.md)). The archives of 2.8.1
> (index format 1, the older names: `NRR_3_hg38+hg38_HGDP1kGP`, `NNN_3_hg38`, ...) stay in
> `lucapinello/crisprme-data` and `lucapinello/crisprme-data-general-index`, unchanged for
> 2.8.1; 2.10.0 reads both formats and still fetches one by its name, naming the index of
> 2.10.0 that replaces it. The end-to-end recipes of the variant panels — the phased default
> HGDP + 1000 Genomes, the HPRC pangenome and sites-only mega (5 sources) — live in
> [`seq_script/merge_panels/README.md`](../seq_script/merge_panels/README.md); this page covers
> the generic build/publish/download flow. `--dictless` is a *publish flag* (it drops the ~152 GB
> per-sample SNP dictionaries; see below), **not** a name marker.

Bulge-enabled CRISPRme searches need a CRISPRitz **index** of the reference
genome. Building it is the single most expensive one-time step of a search. That
index depends only on the **genome + PAM + guide length + bulge count** (not on the guides
or the variant dataset), so it can be built once and reused — or built once by a
maintainer, published to HuggingFace, and downloaded by everyone else.

A **variant-aware** index additionally needs, at post-analysis time, everything the variant
search reads, so that it works with **no source VCFs** and **no separate samples download**:

- the additive **dict-less tiers** — a small Tier-0 **registry** (`registry_<vcf>/`,
  out-of-the-box off-target detection with corrected AF/rsID) plus the larger
  Tier-1 **genotype store** (`genotypes_<vcf>/`, per-sample `Samples`), which
  together replace the ~152 GB per-sample SNP dictionaries;
- the indel logs and the phased indel genotypes (`log_indels_<vcf>/`, `indel_genotypes_<vcf>/`);
  and
- the **samplesID** files the search's `--samplesID` listing expects (the
  combined `<vcf>.samplesID.txt` plus the per-db `<ref>_<db>.samplesID.txt`
  lists), so a downloaded variant index is self-complete.

These belong to the variant panel, not to the PAM: every index of a panel reads the same ones.
Since 2.10.0 they travel once per panel (the Dictionaries bundle and the genotype store, below).

Throughout this page `<vcf>` is the variant-dataset token resolved from the index name's
`+<vcf>` segment — for the recommended default `20-NRR_3_hg38+hg38_HGDP1kGP` it is
**`hg38_HGDP1kGP`**, so the folders are literally `registry_hg38_HGDP1kGP/`,
`genotypes_hg38_HGDP1kGP/`, and `hg38_HGDP1kGP.samplesID.txt`. Substitute your own dataset's
token where you see `<vcf>`, and note that the resolver expects these exact names.

This document covers the download → build → publish workflow. It complements
Section 3.5 of the data-setup guide (`docs/crisprme_data_setup_051826.md`).

## The indexes of 2.10.0

Every index is of N = 3 (built with `--bDNA 2 --bRNA 2`: up to 2 DNA and 2 RNA bulges, the web's
ceiling):

| Index | Guide | PAM | Reference | HGDP + 1000 Genomes | HPRC | mega |
|---|---|---|---|---|---|---|
| SpCas9 **NRR** (NAA, NAG, NGA, NGG; the default) | 20 nt | 3′ | `20-NRR_3_hg38` | `20-NRR_3_hg38+hg38_HGDP1kGP` | `20-NRR_3_hg38+hg38_HPRC` | `20-NRR_3_hg38+hg38_mega` |
| SpCas9 **NGG** (the fast canonical search) | 20 nt | 3′ | `20-NGG_3_hg38` | `20-NGG_3_hg38+hg38_HGDP1kGP` | | |
| SpCas9 **NRG** (NAG + NGG) | 20 nt | 3′ | `20-NRG_3_hg38` | `20-NRG_3_hg38+hg38_HGDP1kGP` | | |
| **General** (PAM-less, 8-nt flank) | 20 nt (21 at 1 + 1) | 3′ | `20-NNNNNNNN_3_hg38` | `20-NNNNNNNN_3_hg38+hg38_HGDP1kGP` | | |
| SaCas9 **NNGRRN** | 22 nt | 3′ | `22-NNGRRN_3_hg38` | `22-NNGRRN_3_hg38+hg38_HGDP1kGP` | | |
| Cas12a **TTTN** | 20 nt | 5′ | `TTTN-20_3_hg38` | `TTTN-20_3_hg38+hg38_HGDP1kGP` | | |
| Cas12a **TTTN** | 23 nt | 5′ | `TTTN-23_3_hg38` | `TTTN-23_3_hg38+hg38_HGDP1kGP` | | |

An index serves a search of a PAM its own holds, pruned to it (`--prune-pam`, which the search
passes by itself): an NRR index serves NGG, NAG, NGA and NRG; NRG serves NGG and NAG; NNGRRN serves
NNGRRT (SaCas9, the CEP290 guides: `PAMs/22bp-NNGRRT-SaCas9.txt`); TTTN serves TTTV and every
Cas12a PAM CRISPRme ships (`PAMs/20bp-TTTV-Cas12a.txt`, `23bp-TTTV-Cas12a.txt`); the general index
serves any 3′ PAM of up to 8 nt (SpRY, VRER, St1Cas9, Nme, a PAM of your own), its search taking
3.3 times a dedicated NGG index's and 1.7 times the default NRR index's (genome-wide on the
reference, `ACTGAAATCTGTAAGCAGGC` at 4 mismatches + 1 DNA + 1 RNA bulge, the engine's search: 26
min for NGG against 8 min on the NGG index, and 15 min for NRR on the NRR index). NRG on HPRC and
mega is served by the NRR index there, pruned: the same rows, more slowly than a dedicated NRG
index (about twice the wall time in the engine's search of chr22's reference sequence: three NRG
guides at 4 mismatches + 1 DNA + 1 RNA bulge, the same 85,563 rows in 107 against 57 s wall, 380
against 110 s CPU; on an HPRC or mega index the post-analysis, the same for both, adds to each),
and a download of 2.8.1's `NRG_3_hg38+hg38_HPRC` or `+hg38_mega` says so. **`NNN_3` is not built
again for 2.10.0**: the general index serves every search it served, at the same speed (a 3-nt
and an 8-nt PAM-less index cost the same, every position being a site); a download of an
`NNN_3` archive of 2.8.1 says so. A guide of another length than its index's is served within
the rule of [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)
(`built + N >= guide + bDNA + (1 if bRNA > 0, else 0)`: the 20-nt general index serves a 21-nt
guide at up to 1 DNA bulge with up to 2 RNA bulges, or at 2 DNA bulges with none, a 22-nt guide
at 1 DNA bulge with none of RNA, or up to 2 RNA bulges with none of DNA, and a 23-nt guide
without bulges; the web, whose cap is one number for both types (`index_fit.bulge_cap`), at up
to 1 of each for a 21-nt guide and none for a longer one), else it needs an index of its own,
which `build-index-only` builds (below). Its rows are not always those of an index built for
the guide's length (**The PAM-less (general) indexes**, below).

The three variant panels are **complementary**, not alternatives:

| Panel | Sources | Data model | Co-occurrence confidence |
|---|---|---|---|
| `hg38_HGDP1kGP` *(recommended default)* | HGDP + 1000 Genomes (4,091 samples = 925 HGDP + 3,166 1000 Genomes), jointly SHAPEIT5-re-phased | **sample-level genotypes, fully phased** | **CONFIRMED** phased cis throughout with **named carrier samples** + exact joint AF (including HGDP's seven genetic regions, with correct chrX PAR allele frequencies) — the broadest fully-phased panel of real observed haplotypes |
| `hg38_HPRC` | HPRC Release 2 pangenome, Minigraph-Cactus graph v2.1 (`hprc-v2.1-mc-grch38`; 2.8.1's archives: v2.0) (232 assembly-derived genomes incl. CHM13) | **sample-level genotypes, phased** (graph-derived) | **CONFIRMED** phased cis with **named carrier samples** + exact joint AF; captures **pangenome-specific** variation absent from short-read panels |
| `hg38_mega` | 5 sources — 1000G-2021 + HGDP + gnomAD v4.1 + TOPMed + All-of-Us | **sites-only** (allele frequencies, no shared samples) | **PUTATIVE** with a conservative **min-AF** joint bound + **per-dataset AF provenance** (`AF_1000G2021 … AF_AoU` + `AF_max`); a co-occurring haplotype here **may or may not exist in any real individual** — no carriers |

**Choosing an index.** You only need **one**. For most use cases pick the **recommended
default `20-NRR_3_hg38+hg38_HGDP1kGP`** — HGDP + 1000 Genomes (4,091 samples), jointly
SHAPEIT5-re-phased, so every reported co-occurrence is **CONFIRMED** cis with **named
per-sample carriers** and exact joint AF throughout, including HGDP's seven genetic regions
and with correct chrX PAR allele frequencies: the broadest fully-phased panel of real
observed haplotypes. Pick
**HPRC** for assembly-derived / pangenome variation (phased → CONFIRMED cis +
carriers). **Escalate to `mega`** only for a widest-provenance worst-case screen when you
must not miss a rare allele from any of five databases: it is **sites-only** (union of allele
frequencies, no genotypes), so every multi-variant / co-occurring off-target is **PUTATIVE**
(a worst-case reconstruction that may not exist in any real individual) with a min-AF bound
and **no named carriers**. A dedicated index of another PAM on the same panel (NGG for speed;
NNGRRN, TTTN) installs beside it in the same working directory and shares its Dictionaries.

All three carry **searchable indels genome-wide** and report **SNP+SNP and SNP+indel**
co-occurring off-targets. Download the reference data once, then any index (the commands of
this page are bash: in zsh, the default shell of macOS, run `setopt interactive_comments` first,
or zsh reads the `#` of a comment as a word):

```bash
# reference genome + annotations + the default reference index 20-NRR_3_hg38 (once)
crisprme.py download --what all  --path .
# recommended default (phased HGDP + 1000 Genomes, 4,091 samples; carriers)
crisprme.py download --what index --index-name 20-NRR_3_hg38+hg38_HGDP1kGP --path .
# HPRC pangenome (phased assemblies, carriers)
crisprme.py download --what index --index-name 20-NRR_3_hg38+hg38_HPRC     --path .
# mega (sites-only, per-dataset AF provenance)
crisprme.py download --what index --index-name 20-NRR_3_hg38+hg38_mega     --path .
# another index of the default panel: its own archive (+ 20-NGG_3_hg38 if not installed)
crisprme.py download --what index --index-name 20-NGG_3_hg38+hg38_HGDP1kGP --path .
```

A first-time unauthenticated download of these many-small-file indexes can hit
HuggingFace's HTTP-429 rate limit; the downloader retries with backoff, and
authenticating (`hf auth login`) avoids it entirely. See METHODS §2 for
the genotyped-vs-sites-only distinction and §4 for how each is scanned.

## Layout on HuggingFace

The indexes of 2.10.0 are under `indexes/` in **`lucapinello/crisprme-indexes-format2`**
(`download` and `publish-index` route there by themselves, an explicit `--hf-repo` or
`CRISPRME_HF_REPO` aside), and the enriched-genome kits they were built from under
`build-inputs/`:

```
indexes/
  20-NRR_3_hg38.tar.gz                           # a reference index: the index alone (+ manifest.json)
  20-NRR_3_hg38+hg38_HGDP1kGP.tar.gz             # a variant index: it, its _INDELS index and manifest.json
  ...                                            # every index of the table above
  dictionaries-bundle_hg38_HGDP1kGP.tar.gz       # the panel's Dictionaries and samplesID lists, once
  dictionaries-bundle_hg38_HGDP1kGP.json         # what the bundle holds: a fingerprint of its files
  genotypes_hg38_HGDP1kGP.tar.gz (+ .json)       # the panel's Tier-1 genotype store, once
  indel_genotypes_hg38_HGDP1kGP.tar.gz (+ .json) # its phased indel genotypes alone (the genotype oracle)
  ... (the same for hg38_HPRC; hg38_mega, sites-only, has a bundle and no genotype store)
  SHA256SUMS                                     # one SHA-256 per .tar.gz: `sha256sum -c --ignore-missing SHA256SUMS`
build-inputs/
  hg38_HGDP1kGP.tar.gz (+ .SHA256SUMS)           # the enriched-genome kits (below)
  hg38_HGDP1kGP.members.md5                      # the md5 of each file a kit holds: `md5sum -c` once unpacked
  hg38_HPRC.tar.gz, hg38_mega.tar.gz (+ .SHA256SUMS, .members.md5)
  hg38_HPRC/                                     # HPRC's normalised VCFs (v2.1), hprc.<chr>.norm.vcf.gz + .tbi, and their SHA256SUMS
```

**Each panel's Dictionaries once.** A variant archive holds the index, its `_INDELS` index and
`manifest.json`, which names the panel's **Dictionaries bundle** (`dictionaries-bundle_<vcf>`:
`registry_<vcf>/`, `log_indels_<vcf>/`, `indel_genotypes_<vcf>/`, mega's `indel_af_<vcf>/`, a
classic publish's `dictionaries_<vcf>/`, and the samplesID lists) and its **genotype store**
(`genotypes_<vcf>`) by name and fingerprint. `download` fetches each once per panel: a second
index of the panel finds them installed, the same, and fetches its own archive (and its PAM's
reference index where that is not installed: `20-NGG_3_hg38` for `20-NGG_3_hg38+hg38_HGDP1kGP`),
not the bundle and the store again. Up to
2.8.1 every variant archive carried its own copy of the Dictionaries (16.5 GB of each HGDP +
1000 Genomes archive) and every install replaced them (`docs/hf/README.md`). Download sizes are
given here once the rebuild has published them.

### The archives of CRISPRme 2.8.1

`lucapinello/crisprme-data` (`indexes/`) and `lucapinello/crisprme-data-general-index` hold the
archives 2.8.1 published (index format 1, the older names), which 2.8.1 and earlier read; each
carries its panel's Dictionaries itself:

```
indexes/                                         # lucapinello/crisprme-data
  NRG_3_hg38.tar.gz, NRG_3_hg38+hg38_{HGDP1kGP,HPRC,mega}.tar.gz
  NGG_3_hg38.tar.gz, NGG_3_hg38+hg38_HGDP1kGP.tar.gz
  NRR_3_hg38.tar.gz, NRR_3_hg38+hg38_{HGDP1kGP,HPRC,mega}.tar.gz
  NNN_3_hg38.tar.gz, NNN_3_hg38+hg38_HGDP1kGP.tar.gz
  genotypes_hg38_HGDP1kGP.tar.gz, genotypes_hg38_HPRC.tar.gz
  SHA256SUMS
indexes/                                         # lucapinello/crisprme-data-general-index
  NNNNNNNN_3_hg38.tar.gz, NNNNNNNN_3_hg38+hg38_HGDP1kGP.tar.gz
  genotypes_hg38_HGDP1kGP.tar.gz                 # the same genotype companion (same SHA-256)
  SHA256SUMS
```

Each is replaced in 2.10.0 by the index the index set names (`NRR_3_hg38+hg38_HGDP1kGP` by
`20-NRR_3_hg38+hg38_HGDP1kGP`; `NNN_3_*` and `NNNNNNNN_3_*` by the general index `20-NNNNNNNN_3_*`;
NRG on HPRC and mega by the NRR index there, more slowly than NRG's own), which a download of it
names. An install holds one build of a panel's Dictionaries: install an archive of 2.8.1 in another
working directory than an index of 2.10.0 of the same panel (`download` refuses to put an archive
of 2.8.1's Dictionaries over a bundle of 2.10.0, and says which installed indexes read another
build's where a bundle replaces them). Their download sizes, as the Hugging Face API reports them
(GB = 10⁹ bytes):

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

**The PAM-less (general) indexes** store the raw genomic flank instead of a matched PAM,
so one index answers any 3′ PAM no longer than that flank when the search passes
`--prune-pam`, which it does by itself (see [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)).
A 5′ PAM (Cas12a) needs a dedicated index: these indexes store the flank on the 3′ side of
the guide, so a 5′ search does not use them (the TTTN indexes of 2.10.0 serve Cas12a).
The general index of 2.10.0, `20-NNNNNNNN_3_*`, and 2.8.1's `NNN_3_*` and `NNNNNNNN_3_*` were
built for a 20-nt guide, with N = 3. That also serves a 21-nt guide within the rule above (so the
shipped 21-nt SaCas9 PAM files work on the 8-nt-flank indexes with `--prune-pam`), but its rows
are not always those of an index built for 21-nt guides: they can differ next to a contig end or
an `N` run, and, with bulges, at a few other sites
([METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering); CRISPRitz's
`docs/INDEX_FORMAT.md`, section 5.1). The fake-indel genome a variant index's `_INDELS` index
searches has an `N` run about every 65 bases (on chr22: 65 for HGDP + 1000 Genomes, 67 for mega,
with CRISPRitz 2.10.0's windows of 30/30; 57 and 59 with 2.9.1's of 25/27). On chr22's of
hg38_mega, with the windows of 25/27, the 21-nt SaCas9 guide at 4 mismatches + 1 DNA + 1 RNA
bulge gave 456 rows of 149,688 in the engine's search of `20-NNNNNNNN_3` that an index built for
21-nt guides does not give (on chr22's reference sequence, the same 2,664,778 rows). For the
rows of a guide length, search an index built for it, which the search takes first where one is
installed; a longer guide needs its own index. They read the same variant data and `Dictionaries/` as the dedicated index for the
same panel — those tiers are PAM-independent — and a `--prune-pam` run over one gives the
dedicated index's report for that PAM: from an index of format 2 (2.10.0's), at any flank
width, contig ends included; from one of format 1 (2.8.1's), byte for byte from a flank of the
PAM's width (`NNN_3` for NRG, verified genome-wide), and from a wider one (K = 8 for a 3-nt PAM)
but for what it cannot index at a contig's very end and start: a + strand PAM with 1 to K − L nt after it to the
end, and a − strand PAM in the first K − L nt. On hg38 at K = 8 that is 145 + 101 NGG and
308 + 241 NRG sites, about 7.5 × 10⁻⁷ of all, none on a primary chromosome (whole-genome, K = 8
lost one row of 27,216 for a guide picked for its perfect match there, and one of 216,978 for
VEGFA site 3); a search that reads a wider flank of format 1 says so on stdout (an index of
format 2 has no such limit: with `--prune-pam` it gives a dedicated index's rows for the same
guide length, contig ends included). CRISPRitz 2.9.1's search also reported a − strand PAM
exactly K − L nt from the start on the + strand at another position; 2.10.0's reads it right in
either format. A wider index of format 1 must have been built by
CRISPRme 2.8.1 (CRISPRitz 2.9.1), which keeps a target whose flank holds an `N`, as every
archive 2.8.1 published was; the search warns of one built before. They are
correspondingly larger, because nothing is excluded at build time: 71 GB for 2.8.1's `NNN_3`
reference index against 9 GB for `NRG_3_hg38`, and 110 GB against 28 GB for the HGDP + 1000
Genomes variant index; a 3-nt and an 8-nt flank cost the same search time (every position is a
site), which is why 2.10.0 publishes the 8-nt one alone. Choose one when you expect to search a
PAM no dedicated index of the set holds.

**Limits of every index of format 1 at contig ends and `N` runs** (dedicated or PAM-less, any
flank; [METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)). Every index
2.8.1 published was built by CRISPRitz 2.9.1 or earlier, in index format 1. CRISPRitz 2.10.0,
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

> **Flank width and compatibility (the repos of 2.8.1).** A 3-nt flank is readable by any
> CRISPRme version: the stored geometry is the same as a dedicated 3-nt PAM index, which is why
> 2.8.1's `NNN_3` is in `lucapinello/crisprme-data`. (Use the rebuilt `NNN_3` archives, published on 2026-10-30, with 2.8.1 or later:
> 2.8.0 keeps the targets with an `N` in their PAM that the rebuild adds, and a search that
> reaches one fails or, with RNA bulges only on CRISPRitz 2.9.0, reports it as a malformed
> 0-mismatch, 0-bulge row; see [`INDEX_PROVENANCE.md`](INDEX_PROVENANCE.md).) A **wider** flank is not — a release predating `--prune-pam` would mis-read it
> — so those indexes live in a separate dataset repo,
> **`lucapinello/crisprme-data-general-index`**, which older defaults never point at.
> `download` and `publish-index` route an index of format 1 by the flank width in its name, so
> a wide-flank index cannot be published to the default repo by accident. Since 2.10.0 an index
> of format 2, or one named for its guide length and side, routes to
> `lucapinello/crisprme-indexes-format2` whatever its flank (no client that cannot read format 2
> reads that repo), is never published in these two, and the Settings catalog lists that
> repo's indexes. An explicit `--hf-repo` or `CRISPRME_HF_REPO` overrides the routing.

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

**What a variant index is published as** (2.10.0, an index of format 2):

- **Its archive** (`20-NRR_3_hg38+hg38_HGDP1kGP.tar.gz`) — the index itself (`<name>/`), its
  `_INDELS` companion, and `manifest.json` at the archive root, which names the panel's
  Dictionaries bundle and genotype store by name and fingerprint.
- **The panel's Dictionaries bundle** (`dictionaries-bundle_hg38_HGDP1kGP.tar.gz`, once per
  panel) — the indel logs (`Dictionaries/log_indels_<vcf>/`), the Tier-0
  `Dictionaries/registry_<vcf>/`, the phased indel genotypes (`indel_genotypes_<vcf>/`), mega's
  indel allele frequencies (`indel_af_<vcf>/`), in a **classic** (non-dict-less) publish the
  per-sample SNP dicts (`dictionaries_<vcf>/`), and the samplesID lists
  (`samplesIDs/<vcf>.samplesID.txt` + the per-db `samplesIDs/<ref>_<db>.samplesID.txt`), with
  its record `dictionaries-bundle_<vcf>.json`.
- **The genotype companion** (`genotypes_hg38_HGDP1kGP.tar.gz`, once per panel) — the big
  Tier-1 store. `download --what index` fetches it unless `--no-genotypes` is given.

`download --what index` fetches the bundle and the companion its manifest names once per panel
(kept, not fetched again, where the install holds them with that fingerprint). An archive of
format 1 (2.8.1's) has two artifacts: its main tarball carries the index, its `_INDELS`
companion, the Dictionaries above but the genotype store, the samplesID lists and the manifest,
and the genotype companion is separate.

**Naming — `--dictless` is a flag, not a name marker.** `publish-index --dictless`
controls *what goes in the bundle* (it drops the ~152 GB per-sample SNP
dictionaries, keeping the compact Tier-0/Tier-1 stores), but it does **not** decorate
the index name. A variant index is published and installed under the name of the folder it
was built in: since 2.10.0 a name with its guide length and side, in the order of the target
on the strand: `<guide>-<pam>_<N>_<ref>+<vcf>` for a 3′ PAM (`20-NRR_3_hg38+hg38_HGDP1kGP`,
extracting to `genome_library/20-NRR_3_hg38+hg38_HGDP1kGP/`), `<pam>-<guide>_<N>_<ref>+<vcf>`
for a 5′ one (`TTTN-23_3_hg38+hg38_HGDP1kGP`); the plain convention `<pam>_<N>_<ref>+<vcf>` for
every index 2.8.1 published (`NRG_3_hg38+hg38_HGDP1kGP`). The search
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
# your write token, in the shell only
export HF_TOKEN=hf_...
crisprme.py publish-index --index genome_library/20-NRR_3_hg38+hg38_HGDP1kGP --dictless
# -> lucapinello/crisprme-indexes-format2 (an index of format 2 goes there by itself):
#    indexes/dictionaries-bundle_hg38_HGDP1kGP.tar.gz (+ .json)  the Dictionaries of the panel (registry (+ variant_count.json),
#                                                               indel logs, indel genotypes) + samplesIDs, once per panel
#    indexes/genotypes_hg38_HGDP1kGP.tar.gz (+ .json)            the Tier-1 genotype store, once per panel
#    indexes/indel_genotypes_hg38_HGDP1kGP.tar.gz (+ .json)      the phased indel genotypes alone (the genotype oracle)
#    indexes/20-NRR_3_hg38+hg38_HGDP1kGP.tar.gz                  the index + _INDELS + manifest (it names the two above)
crisprme.py publish-index --index genome_library/20-NGG_3_hg38+hg38_HGDP1kGP --dictless
# -> indexes/20-NGG_3_hg38+hg38_HGDP1kGP.tar.gz alone: the repo holds the bundle and store of the panel, the same
```

`--dictless` **drops the ~152 GB per-sample SNP dictionaries**
(`dictionaries_<vcf>/`) from the bundle — the Tier-0 registry + Tier-1
genotype tiers replace them — while **keeping the indel logs** (indel
post-analysis still needs them; the tiers are SNP-only). In BOTH modes the small
`registry_<vcf>/` is in the bundle when present, the genotype companion is produced when a
genotype store exists, and the samplesID files this index needs are bundled (so
`download --what index` is self-complete).

Each companion goes up **once per panel**: `publish-index` reads the repo's record of it
(`<name>.json`, a fingerprint of its files) and uploads it only where the repo holds none.
Where the repo holds another build of it — the index was built against other Dictionaries
than the panel's indexes already published read — it refuses, before anything is uploaded;
`--replace-companions` replaces them, for a panel built again as a whole. So build every index
of a panel in one working directory, and publish them from it (`seq_script/rebuild/README.md`).
An index of format 1 is published as up to 2.8.1 (its Dictionaries and samplesIDs in its own
archive, the genotype companion re-uploaded with each), and never an index of format 2 in a
repo 2.8.1 reads.

### Download (end users)

```bash
# fetch the variant-aware index (its archive, the Dictionaries bundle of the panel with the
# samplesIDs, and the genotype companion: the last two once per panel)
crisprme.py download --what index \
  --index-name 20-NRR_3_hg38+hg38_HGDP1kGP --path "$CRISPRME_DIR"
# installs genome_library/20-NRR_3_hg38+hg38_HGDP1kGP/ (+ _INDELS) and the shared
# reference index genome_library/20-NRR_3_hg38/
```

- The index installs under the name it was published under, which the search resolves:
  `20-NRR_3_hg38+hg38_HGDP1kGP` (an archive of 2.8.1, the plain `<pam>_<N>_<ref>+<vcf>`, the
  same way). There is no name decoration to strip — `--dictless` only affected what the
  maintainer put in the archives.
- The **Dictionaries bundle** and the **genotype store** of the panel are fetched once per
  panel: a second index of it (`20-NGG_3_hg38+hg38_HGDP1kGP`) finds them installed, the same
  build its manifest names, and fetches its own archive (and `20-NGG_3_hg38`, its reference
  index, where that is not installed), not them again. A marker,
  `Dictionaries/.<name>.json`, records what is installed and the indexes installed with it; a
  repo whose companion is not the one the index's manifest names is refused, nothing of it
  installed. An archive of 2.8.1, which carries the Dictionaries of its own build, is refused in
  a working directory where a bundle of 2.10.0 of its panel is installed (install it in
  another one); a bundle that replaces another build's Dictionaries names the installed
  indexes that read them.
- The **combined + per-db samplesID lists** are installed into
  `samplesIDs/` from the bundle (an archive of 2.8.1: its main tarball), so no `--what samples`
  / `--what all` is needed. (For a legacy index built before this bundling existed, download
  falls back to fetching the per-db lists from HF and synthesizing the combined file.)
- The **Tier-1 genotype store** is fetched automatically. Add `--no-genotypes`
  to skip the big companion: a search then runs the **population-level** analysis
  from the Tier-0 registry (each site's worst case, with no CONFIRMED cis phasing
  and no carriers but those of an indel, which the indel logs of the bundle
  name; where a site needs a SNP too, its row lists every carrier of the indel)
  until the store is present, and `--per-sample` falls back to it with a warning:
  a per-sample analysis without genotypes would lose the off-targets that need two
  SNPs of one sample together. In a data folder where 2.8.1 built or installed the panel, its
  `Dictionaries/dictionaries_<panel>`, which no download of 2.10.0 replaces, and its
  `Dictionaries/genotypes_<panel>`, where 2.10.0's store is not installed, stay there, and a
  search reads them: remove them. So does 2.8.1's enrichment of the panel (`Genomes/hg38+<panel>`
  and `_INDELS`, from its `build-index-only` or its kit), which a search with no bulges scans in
  place of the index and refuses as another build's: remove it too, or put 2.10.0's kit in its
  place (`docs/TROUBLESHOOTING.md` §16).
- The **CLI search-list files** `list_vcf.txt` (the dataset, `hg38_HGDP1kGP`)
  and `list_samplesID.txt` (its combined samplesID) are written at the install
  root, so a CLI search works out of the box — the same lists the web form builds
  per-search. Installing several variant indexes appends each dataset once.

Then search, pointing at that library (or just run from `$CRISPRME_DIR`):

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
  --guide my_guide.txt --mm 4 --bDNA 1 --bRNA 1 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --index-path "$CRISPRME_DIR/genome_library" \
  --output my_search
```

(4 mismatches + 1 DNA + 1 RNA bulge is the default budget, which `--bDNA`/`--bRNA` given
neither take too; a deeper one, such as `--bDNA 2 --bRNA 2`, which an index of N = 3 serves,
also needs `--max-total-edits` raised to the sum, or alignments over 6 edits are pruned.)

The search finds the prebuilt index (and its tiers + samplesIDs) and skips the
build/enrichment entirely. The variant index and its `_INDELS` index are read from the
`--index-path` library first, then from the working directory's `genome_library/`; the
tiers (`Dictionaries/`), the sample lists (`samplesIDs/`) and the genome (`Genomes/`) are
read from the working directory, so a library kept elsewhere needs those there (links
are fine). A search builds no index in an `--index-path` library: a variant index found
there without its `_INDELS` index beside it stops the search with a message.

## Building a new-PAM index for a published dataset (HGDP+1kGP, mega)

The prebuilt indexes above cover **NRR**, **NGG**, **NRG**, SaCas9 **NNGRRN** (22 nt) and Cas12a
**TTTN** (20 and 23 nt) (not each on every panel) plus the general PAM-less index, which serves
any other 3′ PAM of up to 8 nt through `--prune-pam`, as an NRR index serves the PAMs it holds.
To screen a **different nuclease / PAM or guide length** with a dedicated index (a 5′ PAM other
than a TTTN one, a guide of another length, any IUPAC motif, or NGG on mega, the one of these two
panels it is not published for) against the same re-phased **HGDP+1kGP** or sites-only **mega**
variant data, you build that one index once and reuse it for every later search — you do
**not** re-derive the variant data.

### What's reusable vs. what you rebuild

Only the genome index itself is PAM-specific; everything derived from the variants is not:

| Artifact | PAM-specific? | For a new PAM |
|---|---|---|
| `genome_library/<guide>-<PAM>_<N>_hg38+hg38_<ds>/` (+ `_INDELS`) | **yes** | **rebuilt** |
| Enriched genome `Genomes/hg38+hg38_<ds>/` (+ `_INDELS`) | no | reuse (skip enrichment) |
| `Dictionaries/` (registry, genotypes, indel-GT, logs) | no | reuse (the panel's bundle, from any index download of it) |
| `samplesIDs/` | no | reuse |

> **Why you need an extra file.** The HuggingFace **index** download brings the index of its
> PAM and the panel's `Dictionaries/` + `samplesIDs/` — it does **not** contain the enriched
> genome, so it alone can't build a different PAM. The PAM-independent **enriched-genome kit** is
> published separately under **`build-inputs/<dataset>.tar.gz`**: 2.10.0's, the genome its
> indexes and Dictionaries bundle were built from (every allele of a site coded, its
> fake-indel windows), in `lucapinello/crisprme-indexes-format2`. (2.8.1's kits, in
> `lucapinello/crisprme-data`, were enriched by CRISPRitz 2.9.1: an index built from one does
> not match 2.10.0's Dictionaries.) Download it and you skip straight to indexing the new PAM.
> If instead you want to **reproduce the kit** (or build one for a **different VCF**), the full
> source-to-kit workflow is documented below.
>
> A search of the index itself needs no enriched genome, with no bulges either: a 0-bulge variant
> search scans the enriched genome when every chromosome of it is in place, and otherwise searches
> the installed variant index, with the same results but for a site whose guide lies within 2 nt
> of a contig start or end or an `N` run (N − 1 for an index of budget N): only the scan finds
> its variant version (see the limits of every index at contig ends and `N` runs, above). An
> enrichment in place of another build than the panel's Dictionaries (2.8.1's kit beside an index
> of 2.10.0, say) is not scanned: the search refuses it, as its fake-indel genome tells, and once
> it is removed searches the index ([Troubleshooting §16](TROUBLESHOOTING.md)).

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
| SaCas9 (NNGRRN, 22 nt; the published index's) | `22bp-NNGRRN-SaCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNNNNGRRN 6` |
| SaCas9 (NNGRRT, 22 nt: CEP290; the NNGRRN index serves it, pruned) | `22bp-NNGRRT-SaCas9.txt` | `NNNNNNNNNNNNNNNNNNNNNNNNGRRT 6` |
| Cas12a (TTTN, 5′, 23 nt) | `23bp-TTTN-Cas12a.txt` | `TTTNNNNNNNNNNNNNNNNNNNNNNNN -4` |
| Cas12a (TTTN, 5′, 20 nt) | `20bp-TTTN-Cas12a.txt` | `TTTNNNNNNNNNNNNNNNNNNNNN -4` |
| Cas12a (TTTV, 5′, 23 and 20 nt; the TTTN indexes serve them, pruned) | `23bp-TTTV-Cas12a.txt`, `20bp-TTTV-Cas12a.txt` | `TTTVNNNNNNNNNNNNNNNNNNNNNNN -4`, `TTTVNNNNNNNNNNNNNNNNNNNN -4` |
| PAM-less (the general index: any 3′ PAM of up to 8 nt) | `20bp-NNNNNNNN-NO-PAM.txt` | `NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8` |

For any **other** nuclease, write one line yourself — e.g. a hypothetical 20 nt-spacer 3′ `NGA`:

```bash
printf 'NNNNNNNNNNNNNNNNNNNNNGA 3\n' > PAMs/20bp-NGA-custom.txt
```

### Step 2 — download the enriched-genome kit

```bash
# ~1 GB/dataset download (a few GB once extracted): the SNP-enriched genome + its fake-indel companion, as one tarball.
# the kit of 2.10.0, its SHA-256 checked
crisprme.py download --what build-inputs --dataset hg38_HGDP1kGP --path .
# or by hand:
hf download lucapinello/crisprme-indexes-format2 build-inputs/hg38_HGDP1kGP.tar.gz \
  --repo-type dataset --local-dir .
tar xzf build-inputs/hg38_HGDP1kGP.tar.gz
#  -> Genomes/hg38+hg38_HGDP1kGP/         (SNP-enriched, hg38 + IUPAC overlay)
#  -> Genomes/hg38+hg38_HGDP1kGP_INDELS/  (fake-indel genome)
```

Check the tarball against the `SHA256SUMS` shipped alongside it under `build-inputs/` (`download
--what build-inputs` does), and what it unpacked against its member list,
`build-inputs/<dataset>.members.md5`: `md5sum -c hg38_HGDP1kGP.members.md5` in the folder that
holds its `Genomes/` (the tarball is not byte-reproducible; the list holds for any packing of the
same files). (For **mega** or **HPRC**, substitute `hg38_mega` or `hg38_HPRC`.) The
PAM-independent `Dictionaries/` + `samplesIDs/` you need at **search** time come from any index
download of the panel, which installs its Dictionaries bundle (`download --what index
--index-name 20-NRR_3_hg38+hg38_<ds>`).

### Step 3 — build the new-PAM index (reuse the enriched genome, skip enrichment)

Point `crispritz index-genome` at the downloaded enriched genome with your new PAM. This builds the
PAM-specific index directly from the IUPAC-enriched FASTA — **no re-enrichment** (verified: a new
NGG index off the HGDP+1kGP enriched genome builds in seconds per chromosome):

```bash
# SNP index -> genome_library/<PAM>_<N>_hg38+hg38_<ds>/ (the engine writes the older name; rename it
# to the name build-index-only gives it, <guide>-<PAM>_<N>_... or <PAM>-<guide>_<N>_...)
# -bMax = N = max(bDNA,bRNA)+1 = 3 (the same N as in the index name, e.g. 20-NRG_3)
crispritz.py index-genome hg38+hg38_<ds> Genomes/hg38+hg38_<ds>/ PAMs/<NEW_PAM>.txt \
  -bMax 3 -th <T>
# indel index -> genome_library/<PAM>_<N>_hg38+hg38_<ds>_INDELS/
python PostProcess/pool_index_indels.py \
  Genomes/hg38+hg38_<ds>_INDELS/ PAMs/<NEW_PAM>.txt <PAM_MOTIF> hg38 hg38_<ds> 3 <T>
```

Here `<PAM_MOTIF>` is the PAM region token itself — the IUPAC motif without the spacer `N`s
(e.g. `NGG`, `NRG`, `NNGRRT`, `TTTN`) — the same prefix the index name uses
(`<PAM_MOTIF>_<N>_hg38+hg38_<ds>`).

An index's targets do not depend on the thread count `<T>`, but its `.bin` files do (`buildTST`
sorts in parallel): to build a published index's files byte for byte, give `index-genome` the
threads its archive was built with, `-th 64` for HGDP + 1000 Genomes, `-th 32` for HPRC and `-th
16` for mega and for the reference indexes (`docs/INDEX_PROVENANCE.md`, 2.10.0, "Threads"; the D7
readiness audit's re-index of the HGDP + 1000 Genomes kit was byte-identical at 64 threads, and 6
of its 7 `.bin` files differed at 16). `pool_index_indels.py` builds each `_INDELS` sub-index at
one thread (`-th 1`), whatever its last argument, the number of chromosomes it indexes at once.

The kit and the panel's Dictionaries must come from one release: `pool_index_indels.py` refuses
a chromosome whose fake-indel genome and indel log come from two enrichments (2.8.1's kit, of
CRISPRitz 2.9.1's windows of 25/27, beside 2.10.0's Dictionaries bundle, of 30/30, say: `the
fake-indel genome of <chrom> ... and its indel log ... come from two enrichments`), since every
indel row of that index would be misplaced. Put back the kit that goes with the installed
Dictionaries (`Genomes/hg38+hg38_<ds>` and `Genomes/hg38+hg38_<ds>_INDELS`, from the same
release's `build-inputs/`), or the Dictionaries that go with the kit; `build-index-only`, run on
the dataset, says what to remove, each path absolute ([`TROUBLESHOOTING.md`](TROUBLESHOOTING.md)
§13).

> **One-command alternative.** With the enriched genome already in place, `crisprme.py
> build-index-only --genome Genomes/hg38 --pam PAMs/<NEW_PAM>.txt --bDNA 2 --bRNA 2 --vcf
> VCFs/hg38_<ds> …` detects it (prints `Enriched genome already present`) and skips straight to the
> same two indexes, which it names, since 2.10.0, for the guide length and the side
> (`20-<PAM_MOTIF>_3_hg38+hg38_<ds>` and its `_INDELS` for a 20-nt guide and a 3′ PAM). Use it
> if you also have the dataset's `Dictionaries/`/`VCFs/` staged; otherwise the two
> `crispritz`/`pool_index_indels` calls above need only the enriched kit. For **mega** the kit is
> not optional: an enrichment of its sites-only VCFs holds no fake-indel window (CRISPRitz windows
> an indel only where a sample carries it), so without the kit `build-index-only` refuses the
> dataset beside its indel logs, rather than enrich it anew with no indel off-target (METHODS §4).

### Step 4 — search with the new-PAM index

```bash
crisprme.py complete-search \
  --genome Genomes/hg38 --pam PAMs/<NEW_PAM>.txt \
  --guide my_guide.txt --mm 4 --bDNA 1 --bRNA 1 \
  --vcf list_vcf.txt --samplesID list_samplesID.txt \
  --index-path "$CRISPRME_DIR/genome_library" --output my_search
```

`complete-search` resolves the new index by its name (`<guide>-<PAM>_<N>_hg38+hg38_<ds>`, or the
older `<PAM>_<N>_hg38+hg38_<ds>`) and reuses the sibling `Dictionaries/` + `samplesIDs/` (from
the panel's index download) — no rebuild.

## Reproducing the enriched-genome kit from source (or building one for a new VCF)

The published kit is produced **once** from the source variants; reproduce it (or make a kit for
your **own** VCF) with the full workflow below. Each dataset's authoritative, validated recipe lives
in its own document — summarized here with the key provenance and filters, and linked for the
complete protocol. 2.10.0's kits and indexes of HGDP + 1000 Genomes are built on all of hg38
(455 contigs, as `download --what genome` installs it, as the reference indexes are, and as
2.8.1's NRR, NGG and NNN indexes of the panel were; its NRG index held 24), those of HPRC and
mega on its chromosomes chr1-22, X and Y alone (24 contigs), as 2.8.1's indexes of those panels
were (`seq_script/rebuild/README.md`, rule 1). So a variant search of HGDP + 1000
Genomes finds every site a reference search finds; an alt contig repeats a region of a primary
chromosome (the HLA region, for one), and no VCF of the panel covers it, so a site there can be
reported twice, the second time as a reference row on the alt contig.

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
   crisprme.py build-index-only --genome Genomes/hg38 --pam PAMs/20bp-NRR-SpCas9.txt \
     --bDNA 2 --bRNA 2 --thread <N> \
     --vcf VCFs/hg38_HGDP1kGP --samplesID hgdp1kgp.samplesID.config.txt --path .
   # the Genomes/hg38+hg38_HGDP1kGP/ (+ _INDELS/) it writes IS the enriched kit
   ```
   As-run script: `seq_script/merge_panels/hgdp1kgp_build.sh` (the PAM of its run: `PAM=`; the
   rebuild for 2.10.0 runs it once per PAM of the panel in one working directory,
   `seq_script/rebuild/index_set_2.10.0.tsv`).

### HPRC pangenome (genotyped, phased)

1. **Data source.** HPRC Release 2's Minigraph-Cactus graph **v2.1** on GRCh38, as HPRC decomposes
   it into a VCF (`vg deconstruct` → `vcfbub` → `vcfwave`): `hprc-v2.1-mc-grch38.wave.vcf.gz`
   (2,291,014,302 bytes, md5 `63cc82deb0e010ab98719a9e31dfa352`) and its `.tbi`, 232 samples (231
   individuals and CHM13), phased. 2.8.1's HPRC archives were built from v2.0
   (`hprc-v2.0-mc-grch38`), split the same way; the CHANGELOG (2.10.0, Changed) gives what differs.
   HPRC's experimental L95 VCF is not used.
2. **Split per chromosome** into the form the recipe reads, `ALL.<chr>.HPRC.vcfwave.vcf.gz` and
   its `.tbi` for chr1-22, X and Y (the rebuild's `$HPRC_SRC`), once the download's md5 checks
   (on a mismatch `md5sum -c` prints `FAILED`, and nothing is split):
   ```bash
   U=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38
   curl -fLO $U/hprc-v2.1-mc-grch38.wave.vcf.gz
   curl -fLO $U/hprc-v2.1-mc-grch38.wave.vcf.gz.tbi
   if echo "63cc82deb0e010ab98719a9e31dfa352  hprc-v2.1-mc-grch38.wave.vcf.gz" | md5sum -c -; then
     mkdir -p src
     for c in chr{1..22} chrX chrY; do
       bcftools view --regions $c -Oz -o src/ALL.$c.HPRC.vcfwave.vcf.gz hprc-v2.1-mc-grch38.wave.vcf.gz
       tabix -p vcf src/ALL.$c.HPRC.vcfwave.vcf.gz
     done
   fi
   ```
   Each file holds the records of its chromosome in the genome-wide file, in their order, under
   its header and the two lines bcftools adds for the split, one with the split's date: a split
   made again holds the same records, not the same bytes, and the rebuild checks its own split by
   sha256. v2.1 also has records on chrM and on 32 random and unplaced contigs, which the recipe
   does not read.
3. **Normalise, enrich and index** with `seq_script/merge_panels/hprc_build.sh`, in the image
   (`SIF=`), on a working directory whose `Genomes/hg38` holds chr1-22, X and Y of hg38. Per
   chromosome: `bcftools norm -m -any -f <chr>.fa` (split + left-align) → `+fill-tags -t AF,AC,AN`
   (from the genotypes, the `.` alleles left out) → `sort` → `bgzip` + `tabix`; the samplesID in
   the VCF's sample order, with the populations of the 1000 Genomes 2021 samplesID (`P1KG=`) where
   a sample is there; then `build-index-only` with the PAM of the run (`PAM=`). The rebuild for
   2.10.0 runs it as `seq_script/rebuild/index_set_2.10.0.tsv` writes it (`SRC=$HPRC_SRC`).

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
3. **The 2.10.0 rebuild, from those VCFs** (`seq_script/merge_panels/mega_rebuild.sh <install>
   <merged_vcf_dir> [chrom ...]`, one PAM a run; `PAM`, `BDNA`, `BRNA`, `THREADS` and `RUN`, a
   container prefix, from the environment). Every part of the index reads the records CRISPRitz
   2.10.0's enricher reads:
   - the build input (`PostProcess/prepare_mega_build_vcf.py`): the records whose FILTER is `PASS`
     or `.` only, `INFO/AF` = `AF_max` for SNVs and indels (the indel log's frequency; a record
     with no `AF_max`, as on chrY, the largest of its sources' frequencies, written as both), one
     record of each REF>ALT, the more frequent (the enricher merges two SNP records of one
     allele but gives each indel record a window and a log row of its own; the as-run merge
     wrote 8 long deletions on two lines genome-wide, none on chr22 or chr14);
   - one `build-index-only` **without** `--samplesID` of its copy with a presence genotype,
     `MEGA`, carrying every ALT (`PostProcess/synth_sites_gt.py`), so the enricher windows the
     indels: the enriched genome (the kit), the SNP index, the fake-indel genome, the `_INDELS`
     index and the indel log; then `MEGA` out of the log and the per-sample dictionaries removed
     (`build-index-only` says, on stdout, that a sites-only dataset gets no genotype store);
   - the registry (`PostProcess/build_mega_registry.py` →
     `tier0_registry.compile_registry_from_info_af`: each source a db group with
     `AC=round(AF·AN_nom)`, `AN_nom=2×cohort-N`, the GLOBAL group `AF_max`) and the per-source
     indel frequencies (`PostProcess/build_mega_indel_af.py`), both reading the build input as
     the enricher does (FILTER, every allele of a line, an allele once, no indel within the
     window's left flank of the chromosome's start, read from the indel log); the registry is
     SNP-only, the indels (~24%) the fake-indel path's;
   - `samplesIDs/hg38_mega.samplesID.txt`, the file a search names with `--samplesID`: the five
     sources as pseudo-samples, as the shipped install has it (one already there is kept);
   - `PostProcess/verify_mega_chrom.py`, which reads the source VCF on its own and checks the
     registry, the indel frequencies, the indel log, the enriched genome and the build input
     against it.

   Measured on chr22 and chr14 (METHODS section 2): every check passes, and the indexes and the
   indel log are byte for byte those of the two-build recipe it replaces (step 4, and
   `build-index-only` without `--samplesID` for the SNP index). Every check passes too on chr21,
   chr18, chr10 and chr5, four of the seven chromosomes whose merged VCF holds a record written in
   lower case (chr3, chr4, chr5, chr10, chr12, chr18 and chr21), which the recipe reads in upper
   case, as CRISPRitz 2.10.0 does (read as written, they failed the checks). Run it genome-wide
   only after the release's readiness audit. Design rationale:
   [`docs/DESIGN_mega_index.md`](DESIGN_mega_index.md).

   The recipe builds the dataset whole, its VCFs, registry and indel frequencies too, in a
   working directory that holds none of it: `Genomes/hg38` with the FASTAs of the chromosomes
   to build and no other, and the PAM file. One that holds the dataset's VCFs, or a part the
   recipe builds (its kit, indel logs, registry, indel frequencies or indexes), is refused
   before anything is made there. To put a kit it built in a data folder in place of the
   dataset's (the advice of a build or a search on a sites-only dataset whose kit is not that
   of the installed CRISPRitz), run it in a new working directory, remove from the data folder
   every part of the dataset it built (the kit in place and all that was built from it, which
   the advice names for removal first, and the dataset's VCFs, registry and indel
   frequencies) and move its parts in: [`docs/TROUBLESHOOTING.md`](TROUBLESHOOTING.md) §12 has
   the commands, which remove each of those parts by name, the dataset's indexes by a `find`
   of their name patterns, once the recipe has built the whole dataset and passed its checks,
   and change nothing in the data folder where it fails. They run in bash, pasted into bash or
   into zsh (the default shell of macOS, which reads a pasted `#` as a word and stops a command
   at a pattern that matches nothing), and so does each message's removal of the kit.
4. **Up to 2.8.1: searchable indels** (the shipped mega's). `seq_script/merge_panels/mega_build_indels.sh`
   adds a synthetic `MEGA` pseudo-sample so the enricher materializes fake contigs, builds the
   `_INDELS` genome, then strips the pseudo-sample back out; the registry and the indel
   frequencies were built from every record of the VCF, whatever its FILTER.
5. **Searching it.** A sites-only dataset's candidates in windows of 2<sup>20</sup> sequences or
   more are withheld from its post-analysis by the sites-only density guard
   (`--no-sites-only-density-guard` turns it off), listed in the dense-region files and noted on
   the rows there, and its rows of two or more variants whose frequencies fail the safeguards'
   bound are ranked after the others by the carrier gate (METHODS section 5). With the guard, a
   search of five guides at 4 mismatches and 1 + 1 bulges took 1.35 to 1.46 times the shipped
   mega's on chr22 and 1.53 to 1.73 times on chr14 (IGH), in several
   measurements side by side on shared nodes during 2.10.0's development; unguarded, 4.4 to 4.9
   times on chr22, and the chr14 search was projected at 12 to 16 h against 46 min. The release's
   gates time it genome-wide (`docs/RELEASE_CHECKLIST.md`).

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
# -> lucapinello/crisprme-indexes-format2: indexes/20-NRR_3_hg38.tar.gz

crisprme.py download --what index --index-name 20-NRR_3_hg38 --path "$CRISPRME_DIR"
```

## manifest.json

Every published index carries a small provenance manifest inside its tarball.
A reference index of 2.10.0:

```json
{
  "name": "20-NRR_3_hg38",
  "created_at": "2026-11-20T12:00:00+00:00",
  "index_format": 2,
  "index_meta": {"index_format": 2, "build_pam": "NRR", "side": "3'", "guide_len": 20,
                 "max_bulges": 3, "flank_len": 3, "pamless": false, "four_allele": false,
                 "built_by": {"crispritz": "2.10.0", "crisprme": "2.10.0"}},
  "pam": "NRR",
  "guide_len": 20,
  "pam_side": "3prime",
  "index_bmax": "3",
  "genome": "hg38"
}
```

A dict-less variant index of 2.10.0 adds the tier/companion/self-completeness fields, and
names its panel's companions by name and fingerprint:

```json
{
  "name": "20-NRR_3_hg38+hg38_HGDP1kGP",
  "index_format": 2,
  "pam": "NRR",
  "guide_len": 20,
  "pam_side": "3prime",
  "index_bmax": "3",
  "genome": "HGDP1kGP",
  "has_registry": true,
  "has_indel_genotypes": true,
  "dictless": true,
  "has_genotypes": true,
  "has_samplesids": true,
  "has_bundle": true,
  "dictionaries_bundle": {"name": "dictionaries-bundle_hg38_HGDP1kGP", "fingerprint": "<sha256>"},
  "genotypes": {"name": "genotypes_hg38_HGDP1kGP", "fingerprint": "<sha256>"}
}
```

(An archive of 2.8.1 has the older name, no `index_format` or companions, and its
`Dictionaries/` inside.)

> **Note:** `pam` / `index_bmax` / `genome` are parsed from the index *name* and are
> only meaningful for a **reference** index. For a **variant/merged** index the name
> is `<pam>_<N>_<ref>+<vcf>`, so this name-splitting makes `genome` the last VCF token
> (here `HGDP1kGP`), not the reference assembly — a harmless artifact. The authoritative
> dict-less descriptors are the boolean flags: `dictless`, `has_registry`,
> `has_genotypes`, `has_samplesids`.

- `has_registry` — the Dictionaries (the bundle; an archive of 2.8.1: its main tarball) carry
  the Tier-0 `registry_<vcf>/`.
- `dictless` — the per-sample SNP dicts were excluded (tiers replace them).
- `has_genotypes` — a separate `genotypes_<vcf>.tar.gz` companion was uploaded.
- `display_label` — human-friendly name shown by the web index list / search form.
- `has_samplesids` — the samplesID lists are bundled (in the Dictionaries bundle; an archive
  of 2.8.1: its main tarball), so `download --what index` yields a self-complete install (no
  `--what samples`).
- `dictionaries_bundle`, `genotypes` (2.10.0) — the panel's companions this index was built
  with, by name and fingerprint (the SHA-256 of one line per file, its path, size and SHA-256,
  a gzip header's timestamp read as zeros): `download` fetches each once per panel, and checks
  it against the repo's record of it, `<name>.json`.
- `indel_window` — (a variant index published by CRISPRme 2.10.0 or later) the reference
  bases each fake-indel window holds before an indel and after it, as the index's indel
  logs record them: `[30, 30]` for one CRISPRitz 2.10.0 built (an archive without the field,
  every one published up to 2.8.1, holds 25 and 27). `download --what index` says when an
  archive's windows are not the installed CRISPRitz's, and warns of an index of the dataset
  already installed from another enrichment (an INDELS index of other windows, or a variant
  index of index format 1 beside an archive of 30/30), which a search would read with the
  archive's dictionaries ([METHODS §4](../METHODS.md#snpindel-off-targets-each-carriers-own-haplotype)).

The manifest is surfaced (build timestamp / display name) when the index is
downloaded, and old fields are ignored by older consumers — the index directory
itself is what `complete-search` consumes.

## Notes

- An index is only valid for a matching `--genome`/`--pam`/`--bDNA`/`--bRNA`;
  a different PAM or a higher bulge count needs its own index.
- Variant-enriched (genome + VCF) indexes are dataset-specific. The default
  HGDP + 1000 Genomes enriched index (`20-NRR_3_hg38+hg38_HGDP1kGP`) is prebuilt and
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
PAM and DNA bulges up to 31 bases (27, 26 at an MNV) holds its indel there; the width is read
from each archive's own indel logs, which every INDELS sub-index records it was built from, and a
search refuses an INDELS index installed beside indel logs of other windows (a 2.8.1
archive's index beside the logs of a 2.10.0 one, say: two archives of one panel and one
release share their logs, and are searched as usual), and says when its targets need more
than its windows hold (METHODS §4). A dataset of 2.9.1's windows that you built yourself is
enriched again with 2.10.0's as `build-index-only` run on it says: by removing its enrichment
and all that was built from it (its genomes, dictionaries, indel logs and dict-less tiers, and
every index of it), which the command names by absolute path, then running it again, since a
2.9.1 enrichment differs from a 2.10.0 one in its SNPs too; a dataset is not enriched anew
beside the stores of an earlier enrichment whose records tell other windows (its indel logs, an
INDELS index of it, or a variant index of it in index format 1; a search or a build refuses it,
naming them). The sites-only mega is not enriched anew beside the windows of its kit that its
indel logs or INDELS indexes record, whatever their width: an enrichment of its VCFs, which name
no sample, holds no window, so its windows come from its enriched-genome kit (built with a
synthetic sample: `mega_rebuild.sh`, step 3 of the mega above; up to 2.8.1 its fake-indel half,
`mega_build_indels.sh`), and the advice on it is its kit, its archive and kit of CRISPRme
2.10.0, or that recipe. Post-analysis resolves each indel
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
none, HPRC's chr22 has 9,079 among its 195,813 indel-search records (in v2.0 of its graph,
2.8.1's input). Since 2.10.0 a
per-sample search also reads the store for the SNP rows: a SNP row does not list a carrier
whose haplotype holding its SNPs also holds an indel that changes the target (METHODS §4).
With the gate off it lists them, as up to 2.8.1.

> **Note.** With `mm` at its maximum, raise `--max-total-edits` (e.g. `--max-total-edits
> 6`) so the extra alignment budget for the indel bulge isn't consumed by mismatches,
> or the co-occurring indel off-targets may be pruned (the indel consumes a bulge slot).


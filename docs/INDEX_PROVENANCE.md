# Published-index provenance and validity

Which published index is valid, for which client, and why — with the measurement
behind each claim. Written for the 2.8.1 release, when CRISPRitz 2.9.1 changed
the index builder and a second, non-backward-compatible index family appeared.

Two separate questions per archive:

1. **Readability** — can an older CRISPRme release read it? Decided by the
   **stored flank width**, which equals the PAM length in the archive's name.
   `loadTST` consumes exactly `ceil(len(search_pam)/2)` stored PAM bytes per leaf,
   so a stored width that differs from the search PAM desynchronises the leaf
   stream at the first leaf.
2. **Content** — does the 2.9.1 builder produce the same index? Decided by whether
   the genomic-`N` flank fix (CRISPRitz `5c612c0`) changes that PAM's build.
3. **Completeness** — does the archive hold what a variant search needs? Checked in
   the archive itself: every variant archive was downloaded (SHA-256 equal to the
   Hub's LFS oid) and its members listed, for `_INDELS` bins per sub-index and the
   bundled `Dictionaries/`.

## Verdict

**Every index in `lucapinello/crisprme-data` stores a 3-nt PAM region, so every
one of them is readable by any CRISPRme release, 2.8.0 included.** Wider PAM-less
flanks are the only exception and live in a separate repo.

| Archive | Size | Uploaded | PAM | Readable by 2.8.0 | 2.9.1 builder content | `_INDELS` in the archive |
|---|---|---|---|---|---|---|
| `NRG_3_hg38` | 8.2 GB | 2026-08-25 | NRG | yes | bit-identical | — |
| `NRG_3_hg38+hg38_HGDP1kGP` | 25.8 GB | 2026-10-07 | NRG | yes | bit-identical | 23 sub-indexes, 34 bins |
| `NRG_3_hg38+hg38_HPRC` | 14.2 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 129 |
| `NRG_3_hg38+hg38_mega` | 11.3 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 32 |
| `NRR_3_hg38` | 18.7 GB | 2026-10-25 | NRR | yes | bit-identical | — |
| `NRR_3_hg38+hg38_HGDP1kGP` | 39.9 GB | 2026-10-25 | NRR | yes | bit-identical | 23, 70 |
| `NRR_3_hg38+hg38_HPRC` | 28.6 GB | 2026-10-27 | NRR | yes | bit-identical | 24, 248 — **registry metadata stale, replaced** |
| `NRR_3_hg38+hg38_mega` | 20.8 GB | 2026-10-27 | NRR | yes | bit-identical | **24, 0 — empty, replaced** |
| `NGG_3_hg38` | 3.5 GB | 2026-10-25 | NGG | yes | bit-identical | — |
| `NGG_3_hg38+hg38_HGDP1kGP` | 20.0 GB | 2026-10-25 | NGG | yes | bit-identical | 23, 23 |
| `NNN_3_hg38` | 66.6 GB | 2026-10-25 | none (all-N, width 3) | yes | **differs, see below** | — |
| `NNN_3_hg38+hg38_HGDP1kGP` | 100.9 GB | 2026-10-25 | none (all-N, width 3) | yes | **differs, see below** | 23, 222 |
| `genotypes_hg38_HGDP1kGP` | 5.4 GB | 2026-10-07 | — | yes | not an index | — |
| `genotypes_hg38_HPRC` | 0.5 GB | 2026-10-27 | — | yes | not an index | — |

**2.8.1 replaces `NRR_3_hg38+hg38_mega`**, which shipped an empty indel index, and
**`NRR_3_hg38+hg38_HPRC`**, whose registry lacks the stamp that selects the per-sample
default (both below). The two `NNN_3` archives carry the documented caveat below.
Archives built on one panel now bundle the same `Dictionaries/` (HGDP1kGP: NGG, NRG,
NRR and NNN already agreed on all 117 members by name and exact size; HPRC and mega
did not, and their replacements take the NRG archive's dictionaries byte for byte).

## The 10 dedicated-PAM archives: bit-identical, measured

The N-gate commit claims a dedicated index built by the patched builder is
bit-identical to one from the pristine builder. Measured rather than quoted —
same chr22, same PAM file, both builders, driven through
`crispritz.py index-genome`:

| PAM | pre-fix builder | post-fix builder | verdict |
|---|---|---|---|
| NRG | 276,820,413 B, 3 bins | 276,820,413 B, 3 bins | bit-identical |
| NGG | 136,238,465 B, 2 bins | 136,238,465 B, 2 bins | bit-identical |
| NRR | 533,448,996 B, 5 bins | 533,448,996 B, 5 bins | bit-identical |

Why: the relaxation applies only when the PAM states **no** constraint at any
position. `NGG`/`NRG`/`NRR` each constrain two positions, so a genomic `N`
(encoded `0000`, intersecting nothing) still fails the PAM match exactly as
before.

## The 2 `NNN_3` archives: content differs, and it is reachable only by an
unconstrained search

Built 2026-10-25, before the fix. An all-N PAM is fully unconstrained, so the
relaxation does apply.

**Content** (chr22, same inputs, both builders):

| | bytes | bins |
|---|---|---|
| pre-fix builder | 1,801,785,226 | 16 |
| post-fix builder | 1,801,789,651 | 16 |

A difference of **4,425 bytes in 1.8 GB (0.00025 %)** — the extra entries are
targets whose 3-nt PAM region holds a genomic `N` while the 20-nt guide region
does not, i.e. sites at the edge of an assembly gap.

**Reachability**, measured on a genome with such a site planted on both strands
(the minus-strand case matters: a minus-strand 3′ PAM sits at the *front* of the
extracted window):

| search PAM | pre-fix index | post-fix index | |
|---|---|---|---|
| `20bp-NNN-NO-PAM` (unconstrained) | 1 row | **3 rows** | the two extra rows are the planted `NNN`-PAM sites |
| `20bp-NRG-SpCas9` (`--prune-pam`) | 1 row | 1 row | identical — genomic `N` satisfies neither `R` nor `G` |

So for a **constrained** search PAM — every intended use of a PAM-less index —
the published archive behaves exactly like a freshly built one. Confirmed at full
scale: the published `NNN_3_hg38+hg38_HGDP1kGP`, downloaded from Hugging Face and
searched with `NRG --prune-pam` genome-wide at 4 mismatches + 1 DNA + 1 RNA bulge
with per-sample resolution, reproduces the dedicated-`NRG_3` report
**byte-for-byte**: md5 `6af7a94e4a5f36c26e8c3a0871638ad0`, 189,950 rows, 13,253
indel rows each justified by its own window, 0 unjustified.

**The caveat**, stated plainly: searching the published `NNN_3` with
`20bp-NNN-NO-PAM` as the *search* PAM — "report off-targets regardless of PAM" —
under-reports sites whose PAM bases are assembly `N`. Those are sites whose PAM
sequence is unknown, so a dedicated index would not report them either.

## `NRR_3_hg38+hg38_mega`: published with an empty indel index

The two questions above say nothing about whether an archive actually *contains*
what an indel search needs, so every published variant archive was also checked in
the archive itself: downloaded, SHA-256 equal to the Hub's LFS oid, members listed.
One failed.

| | published 2026-10-27 | replacement (staged) |
|---|---|---|
| `_INDELS` | 24 sub-indexes, **0** `.bin` | 24 sub-indexes, 67 `.bin` (5.29 GB) |
| `Dictionaries/log_indels_hg38_mega` | 24 header-only logs (72–73 B each) | 419.1 MB, as in NRG mega |
| `Dictionaries/indel_genotypes_hg38_mega` | 48 empty files | absent, as in NRG mega |
| `Dictionaries/registry_hg38_mega` | data identical to NRG mega's; index metadata predates the `data_type`/`phased` stamps | NRG mega's |
| SNP index `NRR_3_hg38+hg38_mega` | 346 `.bin` | unchanged |

**Cause.** `build_nrr_ml007.sh` indexed ml007's `mega_gw_build` fake-indel genome,
built 2026-09-02, before the sites-only fix (#191): the enricher writes an indel
contig only when a sample carries the indel, a sites-only VCF has no samples, so
every chromosome's FASTA there is a bare header. `buildTST` crashed on each empty
genome, `os.system` ignored it, and `publish_nrr_mh.sh` uploaded the result with
that build's equally empty `log_indels`. The repo's own driver for this,
`seq_script/merge_panels/mega_build_indels.sh`, adds a synthetic pseudo-sample at
every indel site precisely to avoid it, and was not used.

**Effect.** Every indel search on NRR mega found nothing: silently on 2.8.0; on 2.8.1
the empty logs make it read as a panel with no indels, so also silently. And an
install replaces each `Dictionaries/` directory, so installing NRR mega next to
`NRG_3_hg38+hg38_mega` replaced NRG mega's 419 MB of indel logs with the empty ones,
and NRG mega lost its indel off-targets in that working directory too. (The other two
differences have no effect found: the registry's `.bin` data is byte-identical and,
without the stamps, the CLI and the web both fall back to population-level/sites-only
for this panel; an empty indel genotype store changes nothing for a panel with no
per-sample store, whose co-occurrence rows are PUTATIVE either way.)

**Repair.** The PAM-independent mega fake-indel genome that `mega_build_indels.sh`
produced on 2026-10-10 (ml007 `mega_kit_stage`) was re-indexed with the 2.9.1
builder. As a provenance check it was indexed for NRG first: the result is
**byte-identical to the published NRG mega `_INDELS`** (32 of 32 files,
2,137,365,874 B), so it is that index's exact input. Indexed for NRR, it gives the
67 `.bin` above. The replacement archive is the real `publish_index` output
(upload intercepted) from the published SNP index, that `_INDELS`, and the
`Dictionaries/` + `samplesIDs/` extracted from the published NRG mega archive.

**Verified.** Genome-wide, sites-only, population-level, 4 mismatches + 1 DNA + 1 RNA
bulge, guide `ACTGAAATCTGTAAGCAGGC`, all three runs exit 0:

| run | indel off-target rows (per-dataset AF) | integrated rows |
|---|---|---|
| published NRR mega, CRISPRme 2.8.0 | **0** | 379,288 |
| replacement NRR mega, 2.8.1 | 2,081,489 | 400,459 |
| published NRG mega, 2.8.1 | 1,131,494 | 188,744 |

And at the engine: on each of the 24 `_INDELS` sub-indexes, the replacement's NRR
targets passed through `pam_filter.py` with NRG **equal** the published NRG mega's
targets, 24 of 24 chromosomes, 1,131,494 targets in all (NRR's `N[AG][AG]` contains
NRG's `N[AG]G`, so nothing else would be consistent). Each engine-level total equals
the pipeline's indel row count above, so the pipeline carried every one through.

## `NRR_3_hg38+hg38_HPRC`: complete, with a registry that sets the wrong default

Its indel index is complete (24 sub-indexes, 248 bins) and its index directories
match the publish source member for member. But it was published from ml007's
`hprc_build`, whose registry predates the `data_type`/`phased` stamps the other
archives carry (`reg_*.idx` 1,584 B against NRG HPRC's 1,625 B). `crisprme.py`
chooses the default analysis mode from that stamp — per-sample when it reads
`genotyped-phased`, population-level otherwise — so an NRR HPRC search silently ran
population-level, without named carriers or CONFIRMED cis phasing, on a phased
panel. And an install replaces `Dictionaries/registry_hg38_HPRC`, so installing NRR
HPRC next to NRG HPRC changed NRG HPRC's default too. The replacement takes its
`Dictionaries/` + `samplesIDs/` from the published NRG HPRC archive and changes
nothing else.

## The wide-flank (K = 8) family: a separate repo

`NNNNNNNN_3` stores 8 PAM bytes. A release predating `--prune-pam` reads 2 bytes
per leaf for a 3-nt PAM, desynchronises, and (before 2.9.1's refusal) crashed —
reported as a clean empty result. These therefore go to
`lucapinello/crisprme-data-general-index`, which no older default points at, with
the split derived from the flank width at publish time rather than from a list of
names.

**Publish the 2026-10-28 build**, not the 2026-10-26 one:

| source | built | status |
|---|---|---|
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38` (1,628 bins) | 2026-10-28 17:01 | **publish this** |
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38+hg38_HGDP1kGP` (1,628 bins) + `_INDELS` | 2026-10-28 19:56 / 20:00 | **publish this** |
| `nnn_fromscratch/genome_library/NNNNNNNN_3_*` | 2026-10-26 11:25 | **do not publish** — pre-N-gate-fix |

The repo also needs its own `genotypes_hg38_HGDP1kGP.tar.gz`. `download` fetches
the genotype companion (and the shared reference index) from the repo the index
came from, with no fallback, so without a copy there a per-sample search on the
K = 8 variant index runs with degraded Samples. Upload the main repo's file
unchanged, so both repos hold the same SHA-256.

A wide PAM-less index is the one case the N-gate fix materially changes, and it
bites hardest on the synthetic indel genome, which is `N`-separated roughly every
56 bases: a width-31 window saw 34,690,315 usable starts against a width-26
window's 41,454,105, so **16.3 % of start positions were invisible** to a
pre-fix wide PAM-less build.

**Provenance of the 2026-10-28 build.** `buildTST` compiles from
`mainParallel.cpp` alone (`Makefile_conda`), which includes `pamGenerator.cpp`.
Both are byte-identical between the ml008 tree that produced the index and the
committed 2.9.1 release:

```
1bb392601d1ecd9fee138451e991540b  pamGenerator.cpp
6854efd77db8cd128a3f124422986097  mainParallel.cpp
```

`searchOnTST.cpp` did change afterwards (the width-mismatch refusal) but is not
compiled into `buildTST`, so it cannot affect index content.

## Gaps to close at publish time

- **No `SHA256SUMS` under `indexes/`** in either repo, although `build-inputs/`
  has shipped them since the kits were uploaded. `scripts/hf_sha256sums.py` now
  writes one from the Hub's own LFS checksums (an LFS object is addressed by its
  SHA-256; the oids match both shipped `build-inputs/*.SHA256SUMS`), so the "is
  what is on Hugging Face what we validated?" question is answered by checksum
  instead of by re-running a gate arm. Stage it with the archives, then diff the
  regenerated file after the upload (`docs/hf/README.md`).
- **Nothing stamps the builder version into an index.** The provenance above was
  reconstructed from upload dates, build mtimes and source hashes. A
  `.pam_build`-style sidecar recording the CRISPRitz version and git SHA would
  make this a lookup. Not done in 2.8.1.
- **An indel genotype store built by an enricher at or below 2.9.0 from a
  multiallelic VCF** can hold a fabricated compound allele under a real variant's
  rsID (CRISPRitz `cfcb6a6`). Latent on every CRISPRme-supplied VCF, which are
  split — 0 of 1,998,332 chr22 records carry a comma in `ALT` — and reachable only
  through `gnomAD-converter --multiallelic`. The fix does not retroactively touch
  an already-built store, and nothing invalidates one.

## How to re-run these measurements

The builder comparison needs two engines from the same tree, one with the three
N-gate files reverted to `5c612c0^`, each driven through its own conda-style
layout (`<root>/bin/crispritz.py` + `<root>/opt/crispritz/`). Two gotchas on a
cluster node: `apptainer exec` does **not** activate the image's conda
environment, so call `/opt/conda/bin/python3` by absolute path (`bash -lc` also
re-reads the profile and drops conda's `PATH`), and `/srv/local` is visible inside
the container only if it is the working directory or is bound with
`-B /srv/local`. Guard every comparison against an empty result: two empty hash
listings compare equal, so a failed build reads as "bit-identical" unless the
script counts the files it hashed.

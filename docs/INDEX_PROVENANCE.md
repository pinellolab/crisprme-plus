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
flanks are the only exception and live in a separate repo. Readable is not enough
for the rebuilt `NNN_3` archives: a 2.8.0 search of them can fail or report malformed
rows (*With 2.8.0*, below), so use them with 2.8.1 or later.

Sizes are download sizes as the Hugging Face API reports them, in GB (10⁹ bytes);
"Uploaded" is the date of the commit that last wrote the archive, so for the four
2.8.1 replacements it is the replacement's upload.

| Archive | Size | Uploaded | PAM | Readable by 2.8.0 | 2.9.1 builder content | `_INDELS` in the archive |
|---|---|---|---|---|---|---|
| `NRG_3_hg38` | 8.9 GB | 2026-08-25 | NRG | yes | bit-identical | — |
| `NRG_3_hg38+hg38_HGDP1kGP` | 27.7 GB | 2026-10-07 | NRG | yes | bit-identical | 23 sub-indexes, 34 bins |
| `NRG_3_hg38+hg38_HPRC` | 15.2 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 129 |
| `NRG_3_hg38+hg38_mega` | 12.1 GB | 2026-09-21 | NRG | yes | bit-identical | 24, 32 |
| `NRR_3_hg38` | 20.1 GB | 2026-10-25 | NRR | yes | bit-identical | — |
| `NRR_3_hg38+hg38_HGDP1kGP` | 42.9 GB | 2026-10-25 | NRR | yes | bit-identical | 23, 70 |
| `NRR_3_hg38+hg38_HPRC` | 30.7 GB | 2026-10-30 | NRR | yes | bit-identical | 24, 248; **replaced** (the 2026-10-27 archive's registry metadata was stale) |
| `NRR_3_hg38+hg38_mega` | 25.2 GB | 2026-10-30 | NRR | yes | bit-identical | 24, 67; **replaced** (the 2026-10-27 archive had 24, 0: empty) |
| `NGG_3_hg38` | 3.7 GB | 2026-10-25 | NGG | yes | bit-identical | — |
| `NGG_3_hg38+hg38_HGDP1kGP` | 21.5 GB | 2026-10-25 | NGG | yes | bit-identical | 23, 23 |
| `NNN_3_hg38` | 71.5 GB | 2026-10-30 | none (all-N, width 3) | yes, but a search can fail or report malformed rows (*With 2.8.0*) | **rebuilt with it, replaced** (the 2026-10-25 archive differs) | — |
| `NNN_3_hg38+hg38_HGDP1kGP` | 109.5 GB | 2026-10-30 | none (all-N, width 3) | yes, but a search can fail or report malformed rows (*With 2.8.0*) | **rebuilt with it, replaced** (the 2026-10-25 archive differs) | 23, 245 (the 2026-10-25 archive: 23, 222) |
| `genotypes_hg38_HGDP1kGP` | 5.7 GB | 2026-10-07 | — | yes | not an index | — |
| `genotypes_hg38_HPRC` | 0.6 GB | 2026-10-27 | — | yes | not an index | — |

**2.8.1 replaced `NRR_3_hg38+hg38_mega`**, which shipped an empty indel index,
**`NRR_3_hg38+hg38_HPRC`**, whose registry lacked the stamp that selects the
per-sample default, and **both `NNN_3` archives**, rebuilt with the 2.9.1 builder (all
below). The replacements were uploaded on 2026-10-30, each at its archive's existing
path.

An install replaces each `Dictionaries/` folder of its panel, so every replacement and
the K = 8 archive were checked member by member (SHA-256) against what the published
NRG archive of the same panel installs. The mega and HPRC replacements are
byte-identical to it. The HGDP1kGP archives (the rebuilt `NNN_3` variant archive and
the K = 8 one) carry the published `NNN_3`'s dictionaries byte for byte, and those
match NRG's in content: `registry_` and `samplesIDs/` are byte-identical, and the 46
compressed `log_indels_`/`indel_genotypes_` files differ only in the gzip header's
4-byte timestamp (identical once it is masked), which no reader uses. The published
`NNN_3` has carried exactly these dictionaries since 2026-10-25, and the dedicated NRG
index reproduces its gate report reading them.

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

## The 2 `NNN_3` archives: rebuilt for 2.8.1

The archives published on 2026-10-25 were built before the fix. An all-N PAM is fully
unconstrained, so the relaxation applies to them, and 2.8.1 replaced both with a 2.9.1
build, uploaded on 2026-10-30.

**What the fix changes** (chr22, same inputs, both builders):

| | bytes | bins |
|---|---|---|
| pre-fix builder | 1,801,785,226 | 16 |
| post-fix builder | 1,801,789,651 | 16 |

A difference of **4,425 bytes in 1.8 GB (0.00025 %)** — the extra entries are
targets whose 3-nt PAM region holds a genomic `N` while the 20-nt guide region
does not, i.e. sites at the edge of an assembly gap.

**Reachability in the engine**, counted in CRISPRitz's own search output, before
CRISPRme's post-search PAM filter, on a genome with such a site planted on both strands
(the minus-strand case matters: a minus-strand 3′ PAM sits at the *front* of the
extracted window):

| search PAM | pre-fix index | post-fix index | |
|---|---|---|---|
| `20bp-NNN-NO-PAM` (unconstrained) | 1 row | **3 rows** | the two extra rows are the planted `NNN`-PAM sites |
| `20bp-NRG-SpCas9` (`--prune-pam`) | 1 row | 1 row | identical — genomic `N` satisfies neither `R` nor `G` |

**What a CRISPRme 2.8.1 report shows.** 2.8.1 never reports a target with an `N`
inside its PAM: the post-search PAM filter drops it, and under `--prune-pam` the
engine already does. Measured through `complete-search` (reference-only, 2
mismatches, 1 DNA and 1 RNA bulge, rows of `integrated_results.tsv`) on a synthetic
genome with one control site and four gap-edge sites (PAM region `GNN` and `AGN` on
the plus strand, `TGN` and `GNN` on the minus strand), with an `NNN_2` index from the
pre-fix builder (CRISPRitz `5c612c0^`) and one from the 2.9.1 builder:

| search PAM | pre-fix index | post-fix index | |
|---|---|---|---|
| `20bp-NNN-NO-PAM`, with or without `--prune-pam` | 3 rows | 3 rows | byte-identical; without `--prune-pam` the filter drops the post-fix index's 639 extra targets |
| `20bp-NRG-SpCas9`, without `--prune-pam` | 2 rows | 2 rows | byte-identical |
| `20bp-NG-SpCas9NG` (`--prune-pam`) | 1 row | **5 rows** | the pre-fix index misses the two sites whose PAM (`AG`, `TG`) is followed by an `N`, and a 1-bulge row at each of the other two; a dedicated `NG_2` index gives the same 5 rows, byte-identical |

So in a 2.8.1 search with **any 3-nt 3′ PAM** the 2026-10-25 archives already behaved
exactly like a fresh build. Confirmed at full scale: the published
`NNN_3_hg38+hg38_HGDP1kGP`, downloaded from Hugging Face and searched with
`NRG --prune-pam` genome-wide at 4 mismatches + 1 DNA + 1 RNA bulge with per-sample
resolution, reproduces the dedicated-`NRG_3` report **byte-for-byte**: md5
`6af7a94e4a5f36c26e8c3a0871638ad0`, 189,950 rows, 13,253 indel rows each justified by
its own window, 0 unjustified. What they miss is a **shorter 3′ PAM** searched with
`--prune-pam`: a site whose base just past that PAM is an assembly `N` is a real
target there, since the engine tests the real PAM positions and not the pad
([METHODS §9](../METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering)), and a
dedicated index of that PAM keeps it.

**The rebuild.** Genome-wide, with the 2.9.1 `buildTST` (md5 `76febff1…`), from the
2026-10-25 build's own inputs: the same hg38 FASTAs, the same enriched genome, the
same fake-indel genome (the md5 of all 23 files checked before and at the build), the
same 27-byte PAM file, and `-th 16`, because the builder's parallel sort is not
stable across thread counts. `Dictionaries/` and `samplesIDs/` are the published
archive's own, byte for byte, so only the index directories change. Where they
changed, against the published archive's member list, bin by bin, by exact size:

| index | bins (published → rebuilt) | bins changed | bytes | contigs changed |
|---|---|---|---|---|
| reference `NNN_3_hg38` | 1,628 → 1,628 | 1,222 | +145,540 | 66: exactly the contigs that contain `N` |
| variant `NNN_3_hg38+hg38_HGDP1kGP` | 1,628 → 1,628 | 1,221 | +106,849 | the same 66 |
| `_INDELS` | 222 → 245 | all | +2,105,770,267 (+9 %) | all 23 fake-indel contigs |

The fix only adds a target whose PAM region overlaps assembly `N`, so it predicts
that no contig without an `N` changes. None did: 0 of the 389 `N`-free contigs
changed, in the reference or the variant index. Every fake-indel contig is cut by an
`N` roughly every 56 bases, so all of them change, and that is where nearly all the
growth is; a 3-nt PAM such as NRG still filters those targets out, while a shorter
3′ PAM under `--prune-pam` reports those whose `N` lies just past it.

**Reproducible.** ml008 and ml007 each built the reference index independently, from
the same binary and the same inputs: 1,628 of 1,628 bins byte-identical.

**Verified.** The rebuilt pair, searched by 2.8.1 genome-wide with `NRG --prune-pam`
(4 mm + 1 DNA + 1 RNA bulge, per-sample, HGDP + 1000 Genomes): 189,950 rows, 13,253
indel rows each justified by its own window, 0 unjustified, and `pam_filter` dropping
nothing on any of the three searches. Every column is byte-identical to the gate's
report except `Annotation_COSMIC`: that run used a fresh `download --what
annotations`, where COSMIC stays off until its licence is attested (its log: "COSMIC
EXCLUDED ... dropped 742 COSMIC feature(s)"), so the column reads `NA` throughout.
With that column removed, both reports have md5 `a02f03c5…`. For a user still on
2.8.0, the released 2.8.0 image, unmodified, searching chr22 on the rebuilt reference
index and on the published one gives **byte-identical** reports, both with NRG (the
NNN index as a stand-in, how 2.8.0 serves a PAM it has no dedicated index for; 1,488
rows) and with the unconstrained NNN PAM (8,108 rows). Both hold only because no
hit there has an `N` in its PAM: one that does makes 2.8.0 fail or misreport.

**With 2.8.0.** Use the rebuilt archives with 2.8.1 or later. 2.8.0 reads them, but it
has no `--prune-pam` (it ignores the flag), and its post-search PAM filter treats a
target's `N` as matching any base, so it keeps the targets with an `N` in their PAM
that the rebuild adds. What a search that reaches one does then depends on the
CRISPRitz under 2.8.0. On 2.9.1 the post-analysis fails on the `N` (`cluster.dict.py`).
On 2.9.0, which 2.8.0's image and its `install_from_source.sh` build, `loadTST` leaves
the character of a stored `N` unset (2.9.1 writes `N`), and the engine writes the
target cut short at its `N`: 21 or 22 nt where the others have 23 or 24. A search with
DNA bulges then fails in post-analysis (`IndexError: string index out of range`, then
`adjust_cols.py`: `ValueError: list.remove(x): x not in list`); one with RNA bulges
only finishes, exit 0, and reports the short target as a 0-mismatch, 0-bulge
off-target. Measured with 2.8.0's own code (`crisprme.py`, `PostProcess/`, `pages/`)
on each engine (2.9.1 as the 2.8.1 image ships it; 2.9.0 compiled from its tag with
the 2.8.0 Dockerfile's `g++` line, `searchTST` md5 `c13cfcf0…`), on the synthetic
gap-edge genome and the two `NNN_2` indexes of *What a CRISPRme 2.8.1 report shows*
above, reference-only:

| search PAM and budget | pre-fix index (either engine) | post-fix index, CRISPRitz 2.9.1 | post-fix index, CRISPRitz 2.9.0 |
|---|---|---|---|
| `20bp-NRG-SpCas9`, 2 mm + 1 DNA + 1 RNA (with or without `--prune-pam`) | 2 rows | exit 1: `ValueError: substring not found`, `ERROR: off-targets post-analysis (reference) failed` | exit 1: `IndexError: string index out of range`, the same `ERROR` line |
| `20bp-NNN-NO-PAM`, 2 mm + 1 DNA + 1 RNA | 3 rows | exit 1, the same error | exit 1, the same error |
| `20bp-NRG-SpCas9`, `--mm 2 --bDNA 0 --bRNA 1` (with or without `--prune-pam`) | 2 rows | exit 1, the same error | **exit 0, 4 rows**: 3 of them short 0-mismatch, 0-bulge rows with CFD 0.259, one at `chrA:984` in place of the real 1-bulge alignment there (`GAGTCCGAGCAGAAGAAGA-AAG`, CFD 0.016) |
| `20bp-NNN-NO-PAM`, `--mm 2 --bDNA 0 --bRNA 1` | 3 rows | exit 1, the same error | **exit 0, 5 rows**, 4 of them short 0-mismatch, 0-bulge rows |
| either PAM, `--mm 2 --bDNA 0 --bRNA 0` | 1 row | 1 row | 1 row: a reference search without bulges scans the genome, not the index |

The `--mm 0` runs at `--bDNA 0 --bRNA 1` give the same counts. Genome-wide, the NRR
search behind 2.8.1's PAM-filter fix (CHANGELOG) meets 1 such reference target, 1
variant target and 93,820 indel targets, so a 2.8.0 variant search through the
rebuilt `_INDELS` is likely to reach one. If you searched the rebuilt `NNN_3` with
2.8.0, run that search again with 2.8.1.

**Staged as searched.** The staged archives hold the build that was searched: rsync
checked every transferred file against a whole-file checksum, every file is in the
archive at its path with its exact size, and the gzip stream's CRC held on unpacking.
After the upload every entry of `indexes/SHA256SUMS` equals the Hub's LFS SHA-256 for
that archive, the post-upload check `docs/hf/README.md` prescribes for proving that
the Hub holds what was staged.

## `NRR_3_hg38+hg38_mega`: published with an empty indel index

The two questions above say nothing about whether an archive actually *contains*
what an indel search needs, so every published variant archive was also checked in
the archive itself: downloaded, SHA-256 equal to the Hub's LFS oid, members listed.
One failed.

| | published 2026-10-27 | replacement (on the Hub since 2026-10-30) |
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

## `NRR_3_hg38+hg38_HPRC`: complete, with a registry that set the wrong default

Its indel index is complete (24 sub-indexes, 248 bins) and its index directories
match the publish source member for member. But it was published from ml007's
`hprc_build`, whose registry predates the `data_type`/`phased` stamps the other
archives carry (`reg_*.idx` 1,584 B against NRG HPRC's 1,625 B). `crisprme.py`
chooses the default analysis mode from that stamp — per-sample when it reads
`genotyped-phased`, population-level otherwise — so an NRR HPRC search silently ran
population-level, without named carriers or CONFIRMED cis phasing, on a phased
panel. And an install replaces `Dictionaries/registry_hg38_HPRC`, so installing NRR
HPRC next to NRG HPRC changed NRG HPRC's default too. The replacement, uploaded on
2026-10-30, takes its `Dictionaries/` + `samplesIDs/` from the published NRG HPRC
archive and changes nothing else.

## The wide-flank (K = 8) family: a separate repo

`NNNNNNNN_3` stores 8 PAM bytes. A release predating `--prune-pam` reads 2 bytes
per leaf for a 3-nt PAM, desynchronises, and (before 2.9.1's refusal) crashed —
reported as a clean empty result. These therefore go to
`lucapinello/crisprme-data-general-index`, which no older default points at, with
the split derived from the flank width at publish time rather than from a list of
names.

In `lucapinello/crisprme-data-general-index` (Size and Uploaded as defined above):

| Archive | Size | Uploaded | Readable by 2.8.0 | 2.9.1 builder content |
|---|---|---|---|---|
| `NNNNNNNN_3_hg38` | 82.2 GB | 2026-10-30 | no (8-nt flank) | built with it (the 2026-10-28 build below) |
| `NNNNNNNN_3_hg38+hg38_HGDP1kGP` | 122.7 GB | 2026-10-30 | no (8-nt flank) | built with it (the 2026-10-28 build below) |
| `genotypes_hg38_HGDP1kGP` | 5.7 GB | 2026-10-30 | yes | not an index; the main repo's file, unchanged |

**The 2026-10-28 build was published**, not the 2026-10-26 one:

| source | built | status |
|---|---|---|
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38` (1,628 bins) | 2026-10-28 17:01 | **published** |
| `n8_fixed_install/genome_library/NNNNNNNN_3_hg38+hg38_HGDP1kGP` (1,628 bins) + `_INDELS` | 2026-10-28 19:56 / 20:00 | **published** |
| `nnn_fromscratch/genome_library/NNNNNNNN_3_*` | 2026-10-26 11:25 | **not published**: pre-N-gate-fix |

The repo carries its own `genotypes_hg38_HGDP1kGP.tar.gz`. `download` fetches
the genotype companion (and the shared reference index) from the repo the index
came from, with no fallback, so without a copy there a per-sample search on the
K = 8 variant index would run with degraded Samples. It is the main repo's file
uploaded unchanged: both repos' `SHA256SUMS` list the same SHA-256 for it.

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

## Gaps

- **Closed at the 2026-10-30 upload: `SHA256SUMS` under `indexes/`.** Until then
  neither repo had one, although `build-inputs/` has shipped them since the kits
  were uploaded. Both repos now carry one, written by `scripts/hf_sha256sums.py`
  from the Hub's own LFS checksums (an LFS object is addressed by its SHA-256; the
  oids match both shipped `build-inputs/*.SHA256SUMS`): 14 archives in
  `lucapinello/crisprme-data` and 3 in `lucapinello/crisprme-data-general-index`,
  each entry equal to that archive's LFS SHA-256. The "is what is on Hugging Face
  what we validated?" question is now answered by checksum instead of by re-running
  a gate arm (`docs/hf/README.md`).
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
- **At a position carrying more than one SNP record, only the last record's ALT
  reaches the enriched genome** (CRISPRitz, both enricher ports, every release; in
  its 2.9.1 Known issues). Normalising a VCF splits a multiallelic site into
  biallelic records, and each record rebuilds the position's IUPAC code from the
  reference base and overwrites the SNP dict entry, so an earlier record's ALT is in
  neither the enriched genome nor the dict, and an off-target that needs it is not
  found, or is reported at a lower score, without its carriers (a VCF line that lists
  several ALTs itself is coded with all of them). It applies to any index enriched
  from a VCF with such sites. Measured for HGDP + 1000 Genomes: 2,486,165 of the
  panel's 127,666,429 SNP records (1.95 %), at 2,441,139 positions, are not coded for
  this reason; per guide (NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per-sample,
  genome-wide) this leaves 3 to 7 sites whose true CFD puts them in the top-1000 out
  of it or understated there, and at most one site reaching CFD ≥ 0.5 in the five
  guides tested (EMX1's chr17:12179615, 11 haplotypes, CFD 0.508), and 52 to 324 rows
  per guide (99 % of them at CFD < 0.1) name an individual who carries such an ALT,
  read as the reference base, for a sequence they do not have. On chr22, 45,837 of
  1,757,775 SNP records (2.6 %), at 43,389 positions, lose their ALT this way,
  counted in the normalized source VCF
  (`hgdp1kgp.chr22.norm.vcf.gz`) as the records not last at their position; 5,104
  of them have AC ≥ 41 (AF ≥ 0.5 % of 8,182 haplotypes) and 1,566 AC ≥ 818 (about 10 %).
  Confirmed directly in the shipped data at six such positions, the three of chr22
  with the highest AC, two at AC 257 and 82, and chr20:12296886: the earlier ALT is
  absent from the enriched genome's IUPAC code, the Tier-0 registry, the Tier-1
  genotype store and the per-sample dictionary, and only the last record's ALT is
  held. Two examples: chr22:19142014 has `rs807755` T>A at AC = 8,179 and T>G at
  AC = 1, and the enriched genome reads `K` (G/T), the stores hold G alone;
  chr20:12296886 has `rs6109350` G>A at AC = 1,451 and G>C at AC = 2, the code is
  `S` (C/G), and the stores hold C alone. A fix changes the enriched genomes, and
  with them every affected index: it will be fixed, with rebuilt indexes, in the next
  release.

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

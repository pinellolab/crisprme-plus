# CRISPRme+ 2.10.0

CRISPRme+ 2.10.0 is released with CRISPRitz 2.10.0, its search engine, and with every published
index built again by it. From this release the two share their major.minor version. This page
says what is new, how to upgrade, what the release does not do yet, and what comes next. The full
list, with the measurements behind each entry, is in
[CHANGELOG.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/CHANGELOG.md) (section
2.10.0); [METHODS.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/METHODS.md)
describes the methods, §10 the known limits and §11 how the indexes were built and checked.

## What is new

- **New indexes, in a repo of their own.** Every index is built again by CRISPRitz 2.10.0 and
  published in [`lucapinello/crisprme_plus_indices`](https://huggingface.co/datasets/lucapinello/crisprme_plus_indices),
  named for its guide length and the side of its PAM (`20-NRR_3_hg38+hg38_HGDP1kGP`): SpCas9 NRR
  on the reference, HGDP + 1000 Genomes, HPRC and the mega; NGG, NRG, the general PAM-less index
  (an 8-nt flank), SaCas9 NNGRRN at 22 nt and Cas12a TTTN at 20 and 23 nt on the reference and
  HGDP + 1000 Genomes. Each panel's Dictionaries and genotype store are published once, beside its
  indexes, and fetched once: after `download --what all`, which installs the reference index
  `20-NRR_3_hg38`, the default `20-NRR_3_hg38+hg38_HGDP1kGP` is a 50.1 GB download with them, and a
  second index of the panel its own archive, with its PAM's reference index where that is not
  installed. `NNN_3` is not built again: the general index serves its searches as fast. The archives
  of 2.8.1 stay in `lucapinello/crisprme-data` (and `lucapinello/crisprme-data-general-index`),
  unchanged, and a download by one of their names says which index of 2.10.0 replaces it.
- **What the new indexes hold.** Every allele of a split multi-allelic site is coded: up to 2.9.1
  CRISPRitz coded only the last record's ALT at a position written on several records, so an
  off-target that needs an earlier ALT was missed, or reported lower, without its carriers
  (2,486,165 of the 127,666,429 SNP records of HGDP + 1000 Genomes, 1.95 %). Index format 2 holds
  the sites at contig ends and `N` runs that format 1 loses. Each indel is searched in a window of
  30 bases on either side (25 and 27 up to 2.9.1), so that the indel off-targets of longer targets
  are found (a 21-nt SaCas9 guide with a DNA bulge, Cas12a at 23 nt). REF and ALT are read in
  either case, as VCF 4.x allows, and every record's `AF` by its key (HPRC's chrY indel log had no
  `AF` for 91,950 of its 102,665 indels). The reference indexes and HGDP + 1000 Genomes hold all
  455 contigs of hg38 (the TTTN indexes 454: one unplaced contig has no TTTN site), HPRC and the
  mega chr1-22, X and Y.
- **HPRC v2.1.** The HPRC panel is built from HPRC Release 2's Minigraph-Cactus graph v2.1, where
  2.8.1's was built from v2.0: the same 232 samples, 0.86 % more records, and alleles gained and
  lost, so an HPRC search of 2.10.0 differs from one of 2.8.1 by its data too.
- **The defaults.** SpCas9 NRR is the default PAM and index (it was NRG), and a search given no
  bulges takes 1 DNA + 1 RNA bulge, `assembly-search` included.
- **Indel rows.** A row that starts inside an insertion starts at the insertion's anchor, and the
  merge takes the rows of one insertion for one site; an indel row's variants are numbered along
  its target; `PAM_creation` is read wherever in its repeat the VCF writes the indel; a per-sample
  SNP row no longer lists a carrier whose haplotype holds an indel that changes the target.
- **Population-level rows and the mega.** The carrier gate checks a population-level row's
  variants against phased haplotypes, and the report ranks a combination no haplotype carries
  after the others; it drops no row. The mega is rebuilt from the records the enricher reads, with
  a density guard and a frequency bound that ranks its least likely combinations after the others.
- **Speed.** Coding every allele makes the immunoglobulin loci dense. An exact haplotype pre-check,
  a carriers-first check on the indel side and a post-analysis in pieces over the search's threads
  give the same rows about ten times faster (the README guide's per-sample search of chr22).
  Genome-wide, on 32 threads, the release gate's per-sample searches took 0.95 to 1.3 times 2.8.1's
  time on HGDP + 1000 Genomes and 0.9 on HPRC (Known limitations: the population-level one). A
  search above 1 DNA + 1 RNA bulge on such an index is warned of.
- **One index, many PAMs.** A dedicated index serves the PAMs its own holds (an NRR index serves
  NGG, NAG and NRG), pruned, and a search prunes by itself where only such an index serves it.
- **An independent check of the indel rows on real data**, `scripts/check_indel_rows.py`, which
  writes each sample's haplotypes from the VCF itself; it is a release gate.
- **Ubuntu 26.04.** A search of a Bioconda or source install of 2.8.1 there stopped at
  "Integrating results" (uutils' `pr`); 2.10.0 joins its lists with `paste`, and its conda
  environment installs GNU coreutils. The Docker image was not affected.
- **Security.** The web app follows no link another user put in its data folder, and its Settings
  Retry runs no command read back from a file; a search loads no Python module from the folder it
  is run from. If you serve the web app to other people, upgrade.

**How the indexes were checked.** Before anything was published, the release's gates checked them as
built: the integrity of every index and its counts against the input VCFs; the independent indel
checker on genome-wide searches of each panel, whose findings were fixed or documented below until
none gated; the equivalence of a pruned index with a dedicated one, row for row; and 13 pairs of
searches timed against 2.8.1 ([METHODS §11](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/METHODS.md#11-index-rebuild-of-2100);
[INDEX_PROVENANCE.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/INDEX_PROVENANCE.md)
records each archive, its inputs and their checksums). What a default install takes, measured
from the Hub: a 71.46 GB download, 153.39 GB on disk and a peak of 159.25 GB while the archives
unpack (keep ~200 GB free); the genome-wide variant search 38 GB of memory on 8 threads, 63 GB on
16 (the README's System Requirements).

## Upgrading from 2.8.1

- 2.10.0 reads the indexes of both formats; CRISPRme 2.8.1 reads format 1 only and refuses the
  indexes of 2.10.0.
- **Keep each release's indexes in a data folder of their own.** A download of an index of 2.10.0
  installs its panel's Dictionaries in place of those an archive of 2.8.1 installed. A search of an
  index of 2.8.1 beside them is refused before it writes anything, with what to do: search the
  index of 2.10.0 that replaces it, or keep the old one in a data folder of its own
  ([TROUBLESHOOTING §16](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/TROUBLESHOOTING.md)).
- **After installing an index of 2.10.0 of a panel that 2.8.1 built or installed in that data
  folder**, remove the panel's `Dictionaries/dictionaries_<panel>` (2.8.1's `build-index-only`
  wrote them; 2.10.0's Dictionaries bundle does not carry them, so they stay), and its
  `Dictionaries/genotypes_<panel>` unless the download installed 2.10.0's genotype store: a search
  reads both where they are, and nothing checks them. Remove too 2.8.1's enrichment of the panel,
  `Genomes/hg38+<panel>` and `Genomes/hg38+<panel>_INDELS` (2.8.1's `build-index-only` and its
  enriched-genome kit leave it), or put 2.10.0's kit in its place: a search with no bulges scans it
  in place of the index where it holds every chromosome of the genome (as a `build-index-only` of
  all of hg38 leaves it; the kit holds chr1-22, X and Y alone), a build of an index from it reads
  it, and each refuses it, by its fake-indel genome, as another build's
  ([TROUBLESHOOTING §16](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/TROUBLESHOOTING.md)).
  A new data folder has none of them.
- Re-run any variant search made with 2.8.1 or earlier where its result matters: the coding of
  split sites, the indel windows, the index format, the default PAM and HPRC's data all change
  rows.

## The version policy

CRISPRme and CRISPRitz share their major.minor version from this release: CRISPRme 2.10.z runs on
CRISPRitz 2.10.z, then 2.11 on 2.11. A breaking change in either package (an index format, the
enricher's output, a command line or a file one package reads of the other) bumps the minor of
both; any other change bumps only the patch of the package it changes. So CRISPRme goes from 2.8.1
to 2.10.0, skipping 2.9.x. The Docker image builds CRISPRitz's tag `v2.10.0`; the Bioconda recipe
pins `crispritz >=2.10.0,<2.11`, and CRISPRitz 2.10.0's recipe asks for `crisprme >=2.10.0`, so
that an environment of an earlier CRISPRme does not take an engine whose index format it cannot
read ([RELEASING.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/RELEASING.md#version-policy)).

## Known limitations

These are the main known limits of 2.10.0. METHODS §10 states each, with its scale where it was
measured, and the smaller ones too; TROUBLESHOOTING §12, §14 and §15 say what a user sees and what
to do. The first two concern a search's cost; the next five the rows (which carriers some list, a
target split or written twice, targets left unreported), each with its size where it was measured;
the last three what the real-data check does not check, one refusal's advice on the mega, and what
an index of 2.8.1 keeps.

- **Searches above 1 DNA + 1 RNA bulge in dense regions.** Who: anyone who raises the bulges on a
  variant index. On an index that codes every allele, the search of the immunoglobulin loci and some
  tandem repeats grows by orders of magnitude: on IGL alone, 0.76 million raw targets at 4 + 1 + 1
  and 33.7 million at 4 + 2 + 2, and IGH needed 225 GB at 4 + 2 + 1 (measured during development).
  Keep the default budget there, or `--max-total-edits` 3 or less, and raise the mismatches first.
- **Population-level searches of HGDP + 1000 Genomes in dense loci.** Who: a population-level search
  of that panel, which is per sample by default. The population-level analysis reads no haplotype,
  so it cannot skip the dense loci's windows as a per-sample one does: the README guide's
  genome-wide population-level search with NRG took 2.2 times 2.8.1's time (2 h 29 min against 1 h 8
  min on 32 threads), where the per-sample searches took 0.95 to 1.3 times.
- **On HPRC, an MNV beside an indel on one haplotype.** Who: HPRC searches. No store the search
  reads holds an MNV's phase, so a carrier whose haplotype holding a row's indel also holds an MNV
  in the target is listed for the target that haplotype would have without the MNV: the README
  guide's genome-wide search lists carriers so 1,163 times at 94 sites per sample on indel rows,
  each row at a CFD of 0.085 or less, and 3,334 times on SNP rows. Such a row can describe a
  sequence no haplotype holds: the strongest, chr5:177,802,035 (+) at a CFD of 0.084, lists the 206
  carriers of its deletion, and each of the 317 HPRC haplotypes with the deletion holds the MNV too.
  And a target that a haplotype holds only through an MNV and another indel or MNV is not reported
  for its carriers, as with two indels (below): on chr22, the gate's per-sample check counts 21
  sites of targets held only through two or more indels or MNVs, and, rebuilt from the VCF, 16 such
  sites need an MNV, each target at a CFD of at most 0.0164. Of 2.10.0's panels only HPRC v2.1 has
  MNVs (2.8.1's HPRC archives, of v2.0, have them too).
- **Two indels on one haplotype.** Who: most of all, guides whose sites lie in tandem repeats. A
  fake-indel window holds one indel, so a carrier whose second indel moves a SNP the row names out
  of its target is listed all the same, and one whose own target, through a second indel in a
  repeat, starts a repeat unit or more away is listed on the row of its bases. How often depends on
  the guide: on the release gate's genome-wide checks of HGDP + 1000 Genomes, 9,059 listings of this
  limit for three repeat-rich guides per sample, almost all at VEGFA site 3 in GT repeats (2,972
  carriers on one alternative alignment), 3 for the README guide and none for the other guides
  checked. Rebuilt from the VCF, 8,693 of the 9,059 are of the first shape, and 8 of the second, by
  6 carriers (five on VEGFA site 3, one on EMX1), their own target starting one to three repeat
  units before the row's through one or two more indels in the repeat; most of the rest hold the
  row's target with every variant it names. And a target that a haplotype holds only through two
  indels is not reported for its carriers: on chr22, the gate's per-sample checks find 7 such sites
  for the README guide (3,425 carrier listings, each target at a CFD below 0.0004; 3 with NRG), 343
  for the three repeat-rich guides (VEGFA site 3's up to a CFD of 0.35), 5 on HPRC (each at a CFD of
  0; 16 more there need an MNV, above), 9 for TRAC, 2 for the ACTGA guide, 9 for SaCas9's guide at 2
  mismatches (two of them perfect matches; 43 at 4 + 1 + 1) and none for the CEP290 and Cas12a
  guides. No fix is planned yet: the haplotype scan writes, as the search does, one indel per
  window.
- **A SNP row's carrier beside a long indel.** Who: per-sample searches. A search reads the phased
  indel store within 500 bases of its SNP targets, so a carrier whose indel in the target adds or
  deletes more bases than that, less the target's length, is left off the SNP row even where its
  bases are the row's, unless other SNP targets near it widened what was read: on the release gate,
  one carrier on two rows of HPRC (NA19185, beside a deletion of 1,365 bases), none on HGDP + 1000
  Genomes.
- **A target wholly inside a tandem repeat** can be split over two rows or written twice (in a
  model of microsatellites, 134 and 15 of 10,000 targets at 1 + 1 bulges; no search of real data
  checked has one). No row is lost or wrong.
- **Unphased calls.** A carrier of an unphased call of two ALTs at one position (`1/2`) is listed,
  PUTATIVE, on a row that needs the reference base there; and on a classic install (an index built
  without `--samplesID`) a carrier with an unphased call in a VCF that is otherwise phased is listed
  only on the row of all its variants. No index of 2.10.0 meets either: none of its panels holds a
  `1/2` call (HGDP + 1000 Genomes has no unphased call, HPRC v2.1 only on chrY, each with a missing
  allele, and the mega no genotype; 2.8.1's archives of HPRC v2.0 were not measured), and every
  published variant index of a genotyped panel is built with `--samplesID`. An index you build
  yourself can meet them; their size there was not measured.
- **What the real-data check does not check.** `scripts/check_indel_rows.py` does not read an
  indel row's `Variant_info_spacer+PAM` (each variant's number in its target), nor whether an
  alternative alignment is one the merge should have folded into a better one; the unit tests and
  a planted-site oracle check both. And it refuses (exit 2), unchecked, a run whose VCF holds a
  split SNP site on a chromosome it checks where nothing in the run tells whether its genome codes
  every allele of such a site (METHODS §10 says what tells it).
- **One refusal's advice on the mega.** Who: a data folder that holds the mega's VCFs but not its
  fake-indel genome (`Genomes/hg38+hg38_mega_INDELS`). Where a search there refuses the mega's
  INDELS index (incomplete, or built from other fake-indel windows than its indel logs record), the
  message ends with a command that removes the index and the mega's enriched genome, then builds
  them again with `build-index-only`, which refuses the mega, a sites-only dataset: the kit's
  enriched genome is lost. Do not run it; put the mega's kit back instead
  ([TROUBLESHOOTING §12](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/TROUBLESHOOTING.md)).
- **An index of format 1**, as every archive of 2.8.1, keeps its limits at contig ends and `N`
  runs until it is built again.

Two smaller limits, kept from earlier releases, are in METHODS §10 too: the carriers of an MNV a
row names are listed PUTATIVE, as no store holds an MNV's phase, and an indel row gives its
record's first ID, as does each ALT of a VCF line that lists several.

## Planned for 2.10.1

CRISPRme+ 2.10.1 is to run on CRISPRitz 2.10 and the indexes of 2.10.0, with no index built again
([docs/ROADMAP.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/ROADMAP.md)):

- **The MNVs' phase on HPRC.** Each MNV stored with its haplotype in the phased store (HPRC's
  Dictionaries built again, its index as it is), and read as a second variant of the slot, which
  takes such carriers off the rows of a sequence they do not have. A target held only through an
  MNV and another indel or MNV stays unreported, as one held through two indels does.
- **The haplotype scan of the dense regions.** In a per-sample search above 1 DNA + 1 RNA bulge,
  the samples' own haplotypes searched in the dense regions in place of the variant genome's
  combinations, so that such searches finish; CRISPRitz 2.10.0 already holds what it needs.
- **A faster population-level analysis of the dense loci**, and the timing gate run again.
- **A SNP row's carriers beside a long indel**, the store read as far as a haplotype's own indels
  need; and **the merge in a microsatellite**, one row for one target.
- **The checker and the web app.** The real-data check is to read `Variant_info_spacer+PAM` and
  redundant alternative alignments, hold less memory genome-wide (one process holds the whole
  alternative-alignment table: 106 to 125 GiB on the gate's largest searches), check a search
  beside the mega on chrY, and see a SNP row the rule takes off whole; the web app's handling of
  shared temporary folders and of files read back from a results folder is to be hardened, and its
  History page made to read a `Results/` folder that is a link.

Two indels on one haplotype have no fix planned yet: a fake-indel window holds one indel, and the
haplotype scan writes the lines the search writes. The unphased-call limits, which no index of
2.10.0 meets, and the refusal's advice on the mega have no release named yet (ROADMAP).

## Install

- Docker: `docker pull pinellolab/crisprme:v2.10.0` (amd64 and arm64).
- Bioconda: `crisprme=2.10.0`, once its recipe update is merged.
- From source: `install_from_source.sh` (Linux).

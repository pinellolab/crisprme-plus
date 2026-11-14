# CRISPRme+ 2.10.0

CRISPRme+ 2.10.0 is released with CRISPRitz 2.10.0, its search engine, and with its published
indexes built again by it. From this release the two share their major.minor version. Most of what
changes in results comes from coding every allele of a variant site, from new indel windows and
from a new index format. This page is short: the full list, with the measurements behind each
entry, is in
[CHANGELOG.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/CHANGELOG.md) (section
2.10.0), and [METHODS.md](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/METHODS.md)
describes the methods.

## What changes in results

- **Every allele of a split multi-allelic site is coded.** Up to 2.9.1, CRISPRitz coded only the
  last record's ALT at a position written on several records, so an off-target that needs an
  earlier ALT was missed, or reported lower, without its carriers: 2,486,165 of the 127,666,429
  SNP records of HGDP + 1000 Genomes (1.95 %). CRISPRitz 2.10.0 codes each of them, `Q` where the
  reference base and its ALTs cover A, C, G and T, and every index is built again with it.
- **Index format 2.** It holds the sites at contig ends and `N` runs that format 1 loses: a guide
  within a few bases of a contig end or an `N` run, a PAM on a contig's last base.
- **Fake-indel windows of 30 bases** on either side of an indel (25 and 27 up to 2.9.1), so that
  indel off-targets of longer targets are found: a 21-nt SaCas9 guide with a DNA bulge, Cas12a at
  23 nt, the general index. Every reader takes a window's width from its indel log, so the indexes
  of either release read right.
- **REF and ALT are read in either case**, as VCF 4.x allows: an SNV whose ALT is written in
  lower case was not coded, so its carriers' targets through it were never searched. And every
  record's `AF` is read by its key: HPRC's chrY indel log had no `AF` for 91,950 of its 102,665
  indels.
- **Indel rows.** A row that starts inside an insertion starts at the insertion's anchor, and the
  merge takes the rows of one insertion for one site; an indel row's variants are numbered along
  its target; `PAM_creation` is read wherever in its repeat the VCF writes the indel; a per-sample
  SNP row no longer lists a carrier whose haplotype holds an indel that changes the target.
- **The defaults.** SpCas9 NRR is the default PAM and index (it was NRG), and a search given no
  bulges takes 1 DNA + 1 RNA bulge, `assembly-search` included.
- **Population-level rows and the mega.** The carrier gate checks a population-level row's
  variants against phased haplotypes, and the report ranks a combination no haplotype carries
  after the others; it drops no row. The mega is rebuilt from the records the enricher reads
  (PASS or `.`, AF = AF_max), with a density guard, and a frequency bound that ranks its least
  likely combinations after the others.
- **HPRC v2.1.** The HPRC panel is built from HPRC Release 2's Minigraph-Cactus graph v2.1, where
  2.8.1's was built from v2.0: the same 232 samples, 0.86 % more records, and alleles gained and
  lost.

## New

- **The indexes of 2.10.0**, in a Hugging Face repo of their own,
  [`lucapinello/crisprme-indexes-format2`](https://huggingface.co/datasets/lucapinello/crisprme-indexes-format2),
  named for their guide length and the side of their PAM (`20-NRR_3_hg38+hg38_HGDP1kGP`): SpCas9
  NRR on the reference, HGDP + 1000 Genomes, HPRC and the mega; NGG, NRG, the general PAM-less
  index (an 8-nt flank), SaCas9 NNGRRN at 22 nt and Cas12a TTTN at 20 and 23 nt on the reference
  and HGDP + 1000 Genomes. Each panel's Dictionaries go up once, as a bundle, and each panel's
  enriched-genome kit with the md5 list of its files. `NNN_3` is not built again: the general index
  serves its searches as fast. The archives of 2.8.1 stay where they are, and a download by one of
  their names says which index of 2.10.0 replaces it.
- **Speed.** Coding every allele makes the immunoglobulin loci dense. An exact haplotype
  pre-check, a carriers-first check on the indel side and a post-analysis in pieces over the
  search's threads give the same rows about ten times faster (the README guide's per-sample search
  of chr22). A search above 1 DNA + 1 RNA bulge on such an index is warned of.
- **A dedicated index serves the PAMs its own holds** (an NRR index serves NGG, NAG and NRG,
  pruned), and a search prunes by itself where only such an index serves it.
- **An independent check of the indel rows on real data**, `scripts/check_indel_rows.py`, which
  writes each sample's haplotypes from the VCF itself; it is a release gate.
- **Security.** The web app follows no link another user put in its data folder, and its Settings
  Retry runs no command read back from a file; a search loads no Python module from the folder it
  is run from. If you serve the web app to other people, upgrade.

## Upgrading from 2.8.1

- 2.10.0 reads the indexes of both formats; CRISPRme 2.8.1 reads format 1 only and refuses the
  indexes of 2.10.0.
- **Keep each release's indexes in a data folder of their own.** A download of an index of 2.10.0
  installs its panel's Dictionaries in place of those an archive of 2.8.1 installed. A search of an
  index of 2.8.1 beside them is refused before it writes anything, with what to do: search the
  index of 2.10.0 that replaces it, or keep the old one in a data folder of its own
  ([TROUBLESHOOTING §16](https://github.com/pinellolab/crisprme-plus/blob/v2.10.0/docs/TROUBLESHOOTING.md)).
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

## Known limits

Each is stated with its scale in METHODS §10, and TROUBLESHOOTING §14 and §15 say what a user
sees. The haplotype scan of the dense regions, planned for 2.10.1 with no index rebuild, is to
lift the first two.

- **Searches above 1 DNA + 1 RNA bulge in dense regions** (the immunoglobulin loci, some tandem
  repeats) grow by orders of magnitude on an index that codes every allele: on IGL alone, 0.76
  million raw targets at 4 + 1 + 1 and 33.7 million at 4 + 2 + 2, and IGH needed 225 GB at
  4 + 2 + 1 (measured during development, on fake-indel windows of 25/27). Keep the default budget
  there, or `--max-total-edits` 3 or less.
- **Two indels on one haplotype.** A fake-indel window holds one indel, so a carrier whose second
  indel moves a SNP the row names out of its target is listed all the same (one carrier in seven
  guides on 1000 Genomes chr22, at 2 + 2 bulges).
- **A target wholly inside a tandem repeat** can be split over two rows or written twice (in a
  model of microsatellites; no search of real data checked has one). A fix is planned for 2.10.1.
- **An unphased call of two ALTs at one position** is listed, PUTATIVE, on a row that needs the
  reference base there; on a classic install (an index built without `--samplesID`) an unphased
  carrier is listed only on the row of all its variants.
- **An index of format 1**, as every archive of 2.8.1, keeps its limits at contig ends and `N`
  runs until it is built again.

## Install

- Docker: `docker pull pinellolab/crisprme:v2.10.0` (amd64 and arm64).
- Bioconda: `crisprme=2.10.0`, once its recipe update is merged.
- From source: `install_from_source.sh` (Linux).

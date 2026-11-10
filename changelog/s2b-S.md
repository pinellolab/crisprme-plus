Notes for the merge (not changelog text): stage 2b lane S (speed). Each block goes at the end of
its subsection of [Unreleased]. Three statements already in [Unreleased] were corrected in place
(the --max-degeneracy entry's "the per-sample search changes no row" on the indexes of 2.10.0,
the CRISPRME_IUPAC_CAP entry's BED line count, and the real-data check entry's "97 min on the
indexes this release built", which these steps bring to 8.4 min, both on 8 threads).

### Added
- **A per-sample search skips the raw targets no sample's genotypes can produce (the exact
  haplotype pre-check).** On the indexes CRISPRitz 2.10.0 builds, which code every allele, the
  immunoglobulin loci hold nearly all of a variant search's raw targets (chr22, the README
  guide, NRG, 4 mismatches + 1 DNA + 1 RNA bulge: 766,427 SNP raw targets, 98.6% of them in
  IGL), and almost none is a sequence anybody carries. Before a raw target is expanded into the
  haplotypes its carriers have, it is now tested against the samples' genotypes: every row the
  per-sample post-analysis writes puts, at each variant column, the reference base or an ALT one
  of its carriers carries there, so a raw target gives a row only if one sample, free to take at
  each column the reference base or any ALT it carries there (its union, with no phase), reaches
  the mismatch budget and the PAM; one no sample reaches is skipped, and its reference row is
  still written. The test is per column (a mismatch an ALT fixes, a PAM base an ALT admits) over
  samples held as bits of integers, so it costs a few integer operations a column. It is a
  bound, never a guess: it never skips a target a row comes from (its unit test compares it with
  every choice of alleles of every sample: 3,000 random windows, both strands, bulged or not,
  and 1,500 with a 5′ PAM), and on chr22 the README, NRR, TRAC and VEGFA site 3 searches, and
  the IGL, IGH and IGK loci alone, give the same files with it on and off (44 data files each,
  byte for byte). It runs per sample on the observed-haplotype path and the dictionary path
  within `CRISPRME_IUPAC_CAP`, never on a registry-only dataset (sites-only, or downloaded
  without its genotypes). The real script's tests, each against the same run with the
  pre-check off, skip and keep raw targets with no bulge, an RNA bulge or a DNA bulge, on both
  paths, as one process and in pieces, read a minus-strand target at its forward positions,
  and leave a registry-only dataset's rows alone. On chr22 it skips 755,137 of the README
  guide's 764,435 raw targets with a variant letter (98.8%).
  `complete-search --no-haplotype-precheck` (or `CRISPRME_HAPLOTYPE_PRECHECK` set to 0, false,
  no or off, in any case, as `CRISPRME_POSTPROC_PIECES` below) turns it off; the run says per
  chromosome how many it skipped, or that it is off (stdout, `log_verbose.txt`). See [`METHODS.md`
  §5](METHODS.md#search-speed-on-indexes-that-code-every-allele-the-haplotype-pre-check-the-k-skip-and-the-post-analysis-in-pieces).
- **The indel post-analysis tests an overlaid target's carriers first.** An overlaid target of
  the fake-indel genome gives the rows its indel's carriers have, each holding at a code column
  the reference base or a base the carrier carries there. A second indel that breaks a slot
  takes that slot's carrier away or, unphased, makes the carrier PUTATIVE, which can add a row
  none of its slots holds, but one still made of the carrier's own alleles: a target no carrier,
  free to take any of those at each column, can have within the budget and the PAM gives no
  row, and is now skipped before the second indels of the phased indel store are read for it
  (the window breakers). The indel post-analysis is the same in both analysis modes, and so is
  this check: on chr22 it tells 2,904,141 of the README guide's 2,904,487 overlaid targets (all
  but 346), per sample and population-level. Its unit test compares it with every choice of
  every carrier, DNA bulges at code columns among them, and the real script gives the same
  files with it on and off, on targets that include one at exactly the mismatch budget and one
  whose code column no genotype covers. The switch above turns it off too, and the indel
  post-analysis then says so in `log_verbose.txt` ("indel-snp: carriers-first pre-check off").
- **The K-skip: a population-level search with `--max-degeneracy K` skips, before expanding
  them, the raw targets every row of which would rest on more than K ambiguous columns.** Every
  row the cap lets through rests on at most K, and every row within the budget and the PAM rests
  on at least the mismatches it must fix plus the PAM columns whose reference base the PAM does
  not admit: a raw target needing more than K gives none under the cap, and only its reference
  row is written. It writes what the cap alone writes (its unit test compares the bound with
  every choice; on chr22 the README guide's population-level search with K = 4 skips 753,990 of
  its 764,435 raw targets with a variant letter and takes 8 to 10 min on 16 threads, 17 to 18
  min with the K-skip off, with the same files), and the switch above turns it off too (the log
  then says the K-skip is off). `--max-degeneracy` stays off by default.
- **The post-analysis runs each chromosome in pieces over the search's threads.** Up to now the
  SNP and the indel post-analysis ran one process a chromosome, so the IGL locus of chr22 ran on
  one core. Each chromosome's clustered raw targets are now cut into up to `--thread` ranges of
  lines (at least 16,000 a range for SNPs, 50,000 for indels; `CRISPRME_POSTPROC_PIECE_LINES`),
  the ranges of every chromosome run on one pool of `--thread` workers (bounded by memory as
  before), and a merge (`PostProcess/merge_pieces.py`) joins their files. Each piece starts up
  as one process does over the whole chromosome and analyses its own lines, so it writes the
  rows one process writes for them; the companions keep the first row of each identity, as one
  process does. A CRISPR-Bulge score's last digit depends on which rows share a batch of the
  scorer (80 of 890,496 population-level chr22 rows differed when scored in other batches), so
  the pieces leave their rows unscored and the merge scores them in exactly the batches one
  process forms: the files are byte for byte those of one process (chr22 README per sample and
  population-level against the same searches as one process a chromosome, every data file).
  `CRISPRME_POSTPROC_PIECES=0` (or false, no or off, in any case) runs one process a
  chromosome. The SNP pieces read the chromosome's dictionary keys from a copy written once a
  chromosome: JSON, which loads as data whoever wrote it (it sits in the run's output folder,
  which its group can write under a umask of 002, and loading a pickle runs what it names),
  made 0600 in a pieces' folder made 0700. The indel pieces read the windows of the indel log
  their lines need. The analysis' arguments reach the pieces through a file of a folder only the
  run can write (`tempfile.mkdtemp`), never a name another user of a shared `/tmp` could make
  first.
  `log_verbose.txt` keeps one piece's output and what the others add, but their own run times,
  the notes one process ends a chromosome with (the pre-check's count, the dense windows),
  which the merge writes from the pieces' counts, and the chromosome's run time ("Post-analysis
  SNP of chr22: done in 101.2 s (16 pieces)").
- **A search warns above 1 DNA + 1 RNA bulge on a variant index that codes every allele.** On
  such an index (format 2, CRISPRitz 2.10.0; or none installed yet where the installed
  CRISPRitz builds one) the variant search of the dense loci explodes above the default budget,
  a cost no step after the search removes: on IGL alone (the variant search alone, 16 threads)
  755,963 SNP raw targets at 4 + 1 + 1, 5,097,095 at 4 + 2 + 1 and 33,657,519 at 4 + 2 + 2
  (27.9 GB of memory) at the default cap of 6, where the indel search writes 204,453,499 (157
  GB); IGH needed 225 GB at 4 + 2 + 1; with the cap at the budget's sum, 7,885,211 at 4 + 2 + 1
  (a cap of 7) and 93,753,487 at 4 + 2 + 2 (8; 75 GB), which the warning quotes too, with the
  cap of each count. A `--max-total-edits` of 3 or less keeps a wider budget near the
  default's cost (6 + 2 + 2 wrote 618,209 SNP raw targets there at a cap of 3, 4.5 million at 4,
  14.7 million at 5), so the warning counts the most edits an alignment can have (the smaller
  of the cap and mm + bDNA + bRNA) and is given only above 3: `complete-test`'s 6 + 2 + 2
  benchmarks at a cap of 3 get none, those at 5 one. `complete-search` prints it on stdout and
  keeps it in the run's `.Params.txt` (`Bulge_budget_note`, in the report's bundle); the search
  script prints it in `log_verbose.txt` before the variant search, a web search's too, which
  runs that script alone (`index_fit.py wide-budget`); the web's help and guide say what 2
  bulges of a type cost there.

### Changed
- **`<out>.high_variant_density_regions.bed` holds one line a dense window, written before the
  pre-check.** It held a line per raw target in a dense window, so a window had one line per
  alignment (756,551 lines, 2.2 GB, for chr22's 7,493 windows in the README guide's per-sample
  search; about 13 GB for IGH alone), and its lines depended on which raw targets the
  post-analysis expanded. It is now read from every raw target of the search before the
  pre-check or the K-skip decide anything, one line a window (start, end, guide; 7,493 lines, 19
  MB, for that search), so it is the same with them on and off, and documents the dense regions.
  Its `samples_with_alt` are the samples whose genotype holds an ALT of the window that its
  letters code (each genotype read against its own ALT's allele index; `.` without genotypes):
  per sample the same samples as before; population-level they were the carriers of the window's
  worst-case alignment only. Of the window's lines the one kept is the one `resultIntegrator`'s
  lookup read first (the smallest by variants, letters, bits and greedy: its letters in the case
  of one alignment), so every row's `High_variant_density_region` note is the one it had, and
  the README guide's searches of chr22 give the same integrated results as before, byte for
  byte. The report's count of dense regions counts windows.
- **Which analysis mode is faster depends on the panel.** `complete-search`'s help and its
  per-sample banner, the web form, README, METHODS and the guides called `--population-level`
  "the fast worst-case screen" for a dense, high-edit or genome-wide search. On the indexes
  CRISPRitz 2.10.0 builds it is not, on a dense panel of many samples: the population-level
  analysis reads no haplotype, so nothing skips its dense windows, and on HGDP + 1000 Genomes
  (4,091 samples) the README guide's population-level search of chr22 takes 23 to 32 min against
  8 to 11 min per sample (16 threads). On a small panel it can still be the faster: on HPRC (232
  samples) 9 to 12 min against 14 to 17 min per sample, the post-analysis costing about the same
  in both modes and the per-sample run writing more (the alternative alignments, the
  observed-haplotypes sidecar); with mega beside HGDP + 1000 Genomes, 50 to 55 min against 96 to
  107 min, 54 to 62 min of which the per-sample report (the whole search of both datasets). They
  now say it is a resolution choice whose speed depends on the panel, and what the default
  search costs (README, "What the default search costs"; METHODS.md §5).
- **The CRISPR-Bulge scoring of a batch looks its nulled rows up at once.** The scorer writes -1
  for the rows its model cannot score (two bulges, an `N`, more than 24 columns), nearly every
  row of a dense locus (a DNA and an RNA bulge). It looked each row of a batch up in the list of
  the rows to null: about 40 s a batch of 100,000 such rows (one core of a laptop). It holds
  them as a set now, with the same scores (the README guide's population-level search of chr22
  scores 890,496 rows).
- **A per-sample SNP post-analysis reads the phased indel store around each target's span
  once.** It read it once per raw target, a hundred or more of them a dense window: 227 s on
  chr22 and about an hour on IGH for 10,000 to 25,000 distinct spans. Each span is now read once
  (about 4 s on chr22, 27 s on IGH).

### Fixed
- **A temporary file is made where only the run's user can reach it, never under a name another
  user could make first.** A review of the post-analysis in pieces found it loading a pickle,
  which runs what it names, from the run's output folder (never released: the pieces' copy is
  JSON); a sweep of CRISPRme for the same class found these in released code:
  - The COSMIC-free copy of an annotation whose folder cannot hold it (a read-only or shared
    install) was written into `$CRISPRME_TMPDIR`, the working directory or `/tmp` under a fixed
    name (`.nocosmic.<annotation>`), by a strip writing beside it under one made of its process
    id: another user of a shared `/tmp` could make those names first, a link to a file of the
    user's among them, which the strip then overwrote with the annotation's lines. The copy is
    now made in a folder of the run's own there (`tempfile.mkdtemp`, 0700), removed with the
    annotation temp folders, and the strip never writes through a name that already stands
    (`O_EXCL`): where one does, in an `Annotations/` folder others can write, the copy is made
    in such a folder of the run's own instead.
  - The macOS app installer (`docs/install.sh`) wrote the launcher's AppleScript, which it then
    compiled into the app, under a name it made by adding `.applescript` to the file
    `mktemp -t` made: it writes it in a folder `mktemp -d` makes now.
  - `PostProcess/validate_indel_snp_coords.py`, a developer's check, wrote its overlaid genome
    under a name `tempfile.mktemp` gave: in a folder of its own now.

  A test of the unit-test job (`test_no_code_from_files`) keeps the class out: over every
  module but the tests, no pickle, marshal, shelve or the like (their loads run what a file
  names), no `eval` or `exec`, no YAML load but `safe_load`, no `tempfile.mktemp`; in the shell
  scripts, no `eval`, no `mktemp -u` and no name added to the one `mktemp` made. A folder
  others can write still lets them change a search's inputs and results: on a shared machine,
  keep the data and output folders writable by you alone (`docs/TROUBLESHOOTING.md` §7).

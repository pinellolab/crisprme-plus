<!-- fmt track: [2.9.0] entries for CHANGELOG.md (CRISPRme reads index format 2; the repin to
CRISPRitz 2.10.0), to be merged with the other tracks' fragments. -->

**For the merge, not for the changelog: four notes.**

1. The Changed entry on the repin to CRISPRitz 2.10.0 comes with the last commit of the track, which
   sets the version. Merge that commit only once CRISPRitz v2.10.0 is tagged, and only with or after
   the alleles track (`n290-alleles`): CRISPRitz 2.10.0 writes the four-allele symbol Q, which the
   alleles track's readers take (`819a6d4`) and this track's do not. On this track alone with
   2.10.0, a variant search that reports a target through a position whose reference base and SNP
   ALTs cover A, C, G and T stops in the SNP post-analysis, per-sample and population-level alike:
   `ValueError: substring not found` at `cluster.dict.py` line 164, exit 1. HGDP + 1000 Genomes
   holds 2,448 such positions on chr22 (CRISPRitz's `docs/MULTIALLELIC_ENRICHMENT.md`) and 45,026
   genome-wide; the 300 samples of this track's chr22 runs below hold 134 (in 938 targets of their
   NGG index), and those runs stopped so: per-sample with `feat/2.10.0` 25ee50c, per-sample and
   population-level with a8e52ef. With the alleles track merged in (one conflict,
   `crispritz_bins._BITS`, where both tracks add Q: take the alleles side, which keeps Q out of
   `_PAM_BITS`), the same build and searches, and a reference search, end with exit 0 and nothing on
   stderr (a8e52ef): its reference and variant indexes are this track's, file for file, so its SNP
   post-analysis read the very targets that stopped this track's, and its reference rows are the
   9,211 above. Until the tag, the image and `install_from_source.sh` build CRISPRitz 2.9.1, which
   reads index format 1 only, and every search of an index of format 2 stops with the ERROR the
   first Added entry quotes.
2. The entries name the CRISPRitz of each real run. 2.9.1 is the tag v2.9.1;
   `feat/engine-boundaries` c286bf8 writes format 2 and no Q (the engine's boundary fixes alone);
   `feat/2.10.0` 25ee50c and a8e52ef are CRISPRitz 2.10.0 before its tag, with Q (a8e52ef's
   `buildTST`, `searchTST`, `searchBruteForce` and `enricher` are 25ee50c's, byte for byte; its
   `add-variants` takes only the `.fa` and `.fasta` files of the genome folder); `feat/2.10.0`
   8bda5bb, 7219b55 and 2d1f57b, the engine lane's builds of the last runs, have a8e52ef's
   binaries, byte for byte (7219b55 changes two comments, no engine source; 2d1f57b its CHANGELOG
   and a doc). None ran on the tag: Gate B step 4b of `docs/RELEASE_CHECKLIST.md` runs
   `check_real_engine.py` and `check_boundary_rows.py` on the release's own engine.
3. [2.8.1]'s entries that say CRISPRme 2.8.1 refuses an index of format 2 ("the CRISPRitz planned
   as 2.9.2", "index format 1 only") stay as they are: they say what 2.8.1 does.
4. The note of an index of format 1 says how to build it again where the install can, and never
   to download it again: `download` fetches `lucapinello/crisprme-data`, whose archives are all of
   format 1, and no track points it yet at the rebuilt archives' own location (decision C2). Of a
   variant or INDELS index the install cannot build again from what it holds (a download holds
   neither the dataset's VCFs nor its enriched genome), it says what building it again takes and
   that the dataset's archives rebuilt in index format 2 hold the sites, and gives no command to
   install them. Once a track points `download` at that location, the note (`_how_built_again` in
   the search script, `test_index_format2`) and its entry under Changed can give that command.

### Added
- **CRISPRme reads index format 2, the format CRISPRitz writes from 2.10.0.** An index of format
  2 (CRISPRitz's `docs/INDEX_FORMAT.md`) opens and ends with a magic, records in a header the
  build PAM, the side of the guide it is on, the guide length and the budget, and holds the sites
  an index of format 1 (CRISPRitz up to 2.9.1, every index published until the rebuild for this
  release) cannot: a site whose guide lies within a few bases of a contig end or an `N` run, a
  PAM on a contig's last base, and, PAM-less, a PAM within *K − L* nt of a contig end under
  `--prune-pam` (METHODS §9). CRISPRme 2.8.1 refused every such index ("this CRISPRme reads index
  format 1 only"), and 2.8.0 aborts on one. A search now takes one where the installed CRISPRitz
  reads it, its guide length and side read from its header; under CRISPRitz 2.9.1 it stops with
  one ERROR that names both: `ERROR: <index> was built by a newer CRISPRitz (<file> is in index
  format 2, which CRISPRitz 2.10.0 and later write) than the one installed here, CRISPRitz
  2.9.1, which reads index format 1 only. Install the CRISPRitz this CRISPRme installs with ...`.
  An index of a format newer than 2 is refused as before, naming the formats this release reads.
  Through the real CLI on a synthetic genome of perfect, mismatched and bulged sites at contig ends
  and beside `N` runs, with SNPs and indels by them (`test/index_format2/check_boundary_rows.py`),
  with CRISPRitz 2.10.0 (`feat/2.10.0` 25ee50c and 7219b55) and with the engine of format 2 before Q
  (`feat/engine-boundaries` c286bf8) alike: every search -- reference, and variant per-sample and
  population-level, NGG and the 5′ TTTV, and NGG under `--prune-pam` off an 8-nt flank -- ends with
  nothing on stderr; every row is its genome's bases at its coordinates, in each of its three
  alignments (highest CFD, fewest mismatches and bulges, highest CRISPR-Bulge score); the results
  hold every row of the same searches with CRISPRitz 2.9.1 (13 against 6 NGG rows, 2 against 0 TTTV
  rows, each row more at a boundary) but the one 2.9.1 reports on the wrong strand (a − strand PAM
  at a contig's first base, IDX-14); and the 8-nt flank gives the dedicated index's 13 rows (with
  2.9.1, 5 of its 6). On hg38 chr22 (five guides -- `ACTGAAATCTGTAAGCAGGC`, TRAC, EMX1, VEGFA site 3
  and one whose perfect match starts right after an `N` run -- NGG, 4 mismatches + 1 DNA + 1 RNA
  bulge, `build-index-only` then per-sample over 300 HGDP + 1000 Genomes samples, and a reference
  search): with CRISPRitz 2.9.1 this release gives 2.8.1's integrated results byte for byte (10,942
  and 9,210 rows); with format 2 (`feat/engine-boundaries` c286bf8), every one of those rows and one
  more, that perfect match at chr22:11,118,987, which no index of format 1 holds, and CRISPRitz
  2.10.0 (25ee50c) gives the same 9,211 reference rows, byte for byte; a reference search of the
  format-1 index by 2.10.0 (25ee50c) gives the 2.9.1 rows and the note below; and one of 2.10.0's
  format-2 index by 2.9.1 stops with that ERROR, exit 1.
- **Every index CRISPRme builds records what it holds, in `.index_meta.json`.** `build-index-only`,
  a search's own build and each sub-index of an INDELS index write the index format, the build
  PAM, its side, the flank and guide lengths, the budget, whether the PAM is all `N` and whether the
  index holds the four-allele symbol Q, and the CRISPRitz and CRISPRme that built it, into the
  index folder. Each value comes from the index's own `.bin` headers where they hold it (format 2
  holds them all); an index found in place gets a record of what it holds itself, and keeps one
  it has. It is a dot file, as `.pam_build` is: `crispritz.py search` refuses an index folder that
  holds any other file but a `.bin`, with an ERROR on stdout, no targets file and exit 0 (CRISPRitz
  2.9.1 and 2.10.0 alike, measured). `publish-index` puts the index format, and the record, on the
  archive's manifest, and `download` says which format an archive is in.

### Changed
- **The check that an index is whole follows the index format it was built in.**
  `build-index-only` and the search check an index before they take it: each `.bin` whole, and no
  contig that holds a site without its trees. That check read every `.bin` of format 2 as "not a
  tree (numLeaves=-28747197)", never ran its full-part test on one, and judged every contig by
  format 1's rule, under which a contig whose only sites lie at its ends holds none: on an index of
  format 2 a lost `.bin` of such a contig passed for a contig with no site. It now reads format
  2's header (as searchTST checks it), targets and end magic, and applies the rule of the format
  each index was built in, read from its own `.bin` files, not from the installed CRISPRitz,
  which need not be the one that built it. An index with no `.bin` -- the sub-index of a chromosome
  whose indels hold no site -- is judged by the format it records (`.index_meta.json`), else by its
  sibling sub-indexes', else by format 2's rule (which errs only by asking for a rebuild). A folder
  whose `.bin` files are of both formats, which no build writes, is refused as put together from
  two builds. Checked against the real engine (`test/index_format2/check_real_engine.py`): with
  `feat/engine-boundaries` c286bf8 (format 2 before Q), 224 of 224 checks over 60 random contigs,
  7 PAMs (NGG, NRG, NNGRRT, the 5′ TTTV, PAM-less K = 3, 4 and 8) and budgets 0 to 3, every count
  of a contig's sites equal to the targets the engine wrote and all 45,304 cuts of its files told;
  with CRISPRitz 2.10.0 (`feat/2.10.0` 25ee50c), on contigs holding Q among their bases, 252 of
  252 (those 224 and, for each index, a file of it damaged in nine ways, each read as
  `crispritz.py index-info` reads it) and all 255,442 cuts told; with 2.9.1, 196 of 196.
- **A search takes an index the installed CRISPRitz reads ahead of one it would stop on.** Of a
  library holding an index of format 2 and one of format 1 that serve a search, the scan took the
  format-2 one when its budget was the smaller, and then stopped on it under CRISPRitz 2.9.1
  (PKG1-scan-takes-format2); it took the format-1 one before 2.8.1, misreading the other. The
  search script's scan and complete-search's now rank an index the installed CRISPRitz reads
  first, and take the other only where none serves.
- **A search of an index of format 1 by CRISPRitz 2.10.0 says what it lacks, and how those sites are
  searched.** On stdout, once for each index it takes (a variant search with bulges takes three: the
  reference, variant and INDELS indexes): `NOTE: <index> is in index format 1 (built by CRISPRitz
  2.9.1 or earlier): it lacks the sites whose guide lies within a few bases of a contig end or an N
  run, which an index built by CRISPRitz 2.10.0 or later holds (METHODS section 9).` It goes on `To
  search them, move <path> out of the way and build it again with crisprme.py build-index-only
  (which keeps an index it finds in place)` only where the install builds the index again from what
  it holds, and a search takes the one built again: a reference index from the genome folder, and a
  variant or INDELS index from the dataset's enriched genome or fake-indel genome in place (`, from
  Genomes/<ref>+<vcf>.`) with its VCFs, which `build-index-only` takes too, in the working
  directory's `genome_library`. Of a variant or INDELS index the install cannot build again it says
  what building it again takes that the install does not hold -- the dataset's VCFs; its enriched
  genome, without which the whole dataset is enriched again from the VCFs first; or its fake-indel
  genome -- and that `Where <vcf> is published, its archives rebuilt in index format 2 hold them`.
  An install of downloaded indexes (`download --what index`, the README's Quickstart) holds neither
  VCFs nor an enriched genome: the note had its variant and INDELS indexes moved out of the way too,
  and `build-index-only` then stopped on the dataset (`Error: VCF folder 'hg38_HGDP1kGP' ... not
  found under .../VCFs/`, exit 1). With CRISPRitz 2.10.0 (`feat/2.10.0` 8bda5bb), genome-wide on a
  mirror of such an install of `NRG_3_hg38+hg38_HGDP1kGP` (one guide, NRG, 4 mismatches + 1 DNA + 1
  RNA bulge, per-sample), the search's notes are those three: the reference index's says to move it
  and build it again, and those of the variant and INDELS indexes that `building it again with
  crisprme.py build-index-only takes the dataset's VCFs (VCFs/hg38_HGDP1kGP), which this install
  does not hold`, as `build-index-only` of the dataset says (exit 1, `Error: VCF folder
  'hg38_HGDP1kGP' from .../list_vcf.txt not found under .../VCFs/`). On hg38 chr22 (the guides and
  samples above), the install CRISPRme 2.8.1 built with CRISPRitz 2.9.1, laid out as a download (no
  VCFs, no enriched genome, no per-sample dictionaries), gives 2.8.1's 10,942 rows byte for byte and
  the three notes; `build-index-only` of the dataset stops as its notes say (exit 1, every `.bin`
  unchanged); with the reference index moved out of the way as its note says, `build-index-only`
  builds it in format 2, and the per-sample search gives the 10,943 rows of format 2 above, byte for
  byte (the row more, the perfect match at chr22:11,118,987, is a reference site), with the other
  two notes, and the reference search the 9,211, with none. The same install as built (its VCFs and
  enriched genomes in place): the notes of its variant and INDELS indexes say to build them again
  from `Genomes/hg38+hg38_hgdp300` and `Genomes/hg38+hg38_hgdp300_INDELS`; `build-index-only` alone
  keeps the three indexes, every `.bin` unchanged; with each moved out of the way as its note says,
  it builds all three in format 2, and the per-sample and reference searches give the same 10,943
  and 9,211 rows, with no note. With its VCFs and no enriched genome, those two notes say that
  building them again takes the enriched genome, or else the VCFs, from which the whole dataset is
  enriched again first. Each of these searches ends with exit 0 and nothing on stderr. A search
  takes an index of an `--index-path` library other than the working directory's `genome_library`
  before one built again there -- its reference index from that library alone, and its INDELS index
  from beside its variant index -- and of one there the note said the same: built again as it said,
  the index was not searched, and the next search stopped (`ERROR: no matching reference index under
  --index-path`; `ERROR: INDELS index ... is not beside its variant index in ... (--index-path), and
  a search builds no index there`). Of such an index it now goes on `A search builds no index in the
  --index-path library: to search them, replace <path> with one built again by crisprme.py
  build-index-only (which builds into the genome_library folder of its --path, and keeps an index it
  finds there in place)` (with `from Genomes/<ref>+<vcf>` after `build-index-only` for a variant or
  INDELS index), and `, or run this search without --index-path` (` and --prune-pam (which needs
  it)` under `--prune-pam`) where the search then builds every index it needs: a reference search,
  or a variant one with the dataset's VCFs, enriched genome and fake-indel genome in place; without
  them it finds no variant index, and stops (`Error: no variant index here serves this search ...`).
  With CRISPRitz 2.10.0 (`feat/2.10.0` 7219b55), of indexes CRISPRitz 2.9.1 built, in a library
  outside the working directory: a reference search (NGG, 2 mismatches + 1 DNA + 1 RNA bulge; two
  contigs, each with a perfect match from an `N` run's edge) with `--index-path` gives the 2 rows of
  format 1 and that note; with the index moved out of the way and built again in the working
  directory, as the note said, it stops (exit 1); with that one put in its place, it gives 4 rows
  and no note, as does the search without `--index-path`, which builds the index itself, and the
  search of a library named `genome_library` in another folder, the index built there with
  `build-index-only --path`. A per-sample search of the dataset of `check_boundary_rows.py`, its
  three indexes in such a library and its VCFs, enriched genome and fake-indel genome in the working
  directory, gives the 7 rows of format 1 and the three notes; with the INDELS index moved out of
  the way and built again in the working directory, it stops (exit 1); with each replaced as they
  say, and without `--index-path`, it gives the 13 rows of the same search on indexes built in
  format 2, byte for byte, and no note. With the working directory laid out as a download (the
  dictionaries, no VCFs and no enriched genome), the reference index's note offers no search without
  `--index-path`, which stops (exit 1, `Error: no variant index here serves this search`), and the
  notes of the variant and INDELS indexes are those above; with the reference index replaced as its
  note says, those two remain. Each of these searches that does not stop ends with exit 0 and
  nothing on stderr. `test_index_format2` does the same through `build-index-only` and the search
  script on stub engines. The note `--prune-pam` gives of a flank wider than its PAM, and the
  warning of a flank built before CRISPRitz kept an `N` in it, are given of an index of format 1
  only.
- **The `--prune-pam` note of a flank wider than its PAM says that a PAM is reported on the wrong
  strand only under CRISPRitz 2.9.1.** Of an all-N index of format 1 whose flank of *K* nt is
  wider than the PAM's *L*, the note, on stdout and in the run's `.Params.txt` (`Prune_pam_note`),
  said that a − strand PAM exactly *K − L* nt from a contig's start `is reported on the + strand
  at another position`. CRISPRitz 2.9.1's search does so: that PAM's window starts at the contig's
  first base, which it reads as a window of the + strand. 2.10.0's reads it as the − strand window
  it is, in an index of either format, as it reads that of a − strand PAM at a contig's first
  base, so under 2.10.0 the note said it of a site the search reports right. It now says it only
  under a CRISPRitz that reads index format 1 alone (2.9.1), and the PAMs such a flank does not
  index under either. With CRISPRitz 2.10.0 (`feat/2.10.0` 2d1f57b), of an 8-nt-flank index that
  CRISPRitz 2.9.1 built (`build-index-only` with a 28-nt all-`N` PAM file: index format 1) of a
  contig holding an NGG site on the − strand whose PAM starts 5 nt from its start and one on the +
  strand within it, a `--prune-pam` NGG search gives both on their own strands and at their own
  positions, the two rows a dedicated NGG index gives under either engine, and the note without
  that clause; with 2.9.1 it gives the − strand site at 2 on the + strand (the dedicated index: at
  5 on the − strand), and the note with it. An 8-nt-flank index that 2.10.0 built (format 2) gives
  the same two rows and no note. Each of these builds and searches ends with exit 0 and nothing on
  stderr. METHODS §9, `docs/PRECOMPUTED_INDEXES.md`, `docs/DYNAMIC_PAM_INDEX_COMPAT.md` and the
  general-index card in `docs/hf/` say that this misread is 2.9.1's, and `complete-search --help`
  gives the limits of a wider flank as those of index format 1.
- **The side of an index of format 2 is the one its header records.** The search script,
  complete-search, complete-test, the coverage check, `build-index-only` and the web read it there
  before `.pam_build`, which an index built by `crispritz.py` itself, or downloaded, may lack.
  `build-index-only` read `.pam_build` alone: over a 5′ PAM-less index of format 2 that had none,
  it took the shipped 3′ `PAMs/20bp-NNN-NO-PAM.txt`, said the index was already present and
  recorded that file on it, whose 3′ geometry the web then read, and a 3′ search refused the index
  as a 5′ build "(as its .pam_build records)", which recorded 3′; and it refused a 3′ index of
  flank 4 that had none, given a 3′ PAM file of that flank, as "a PAM-less index that records no
  side of the guide". It now refuses a PAM file of the other side than the header's, saying the
  header records it, and records one of its side; the web takes the index's geometry from its
  header, and the PAM name of `.pam_build` only where that record's geometry is the header's; and a
  search that refuses an index for its side says which of the two it read. With CRISPRitz 2.10.0
  (`feat/2.10.0` 8bda5bb) on a genome of three contigs: of a 5′ PAM-less index that `crispritz.py
  index-genome` built, `build-index-only` with the shipped 3′ PAM file stops (exit 1: `Error: ...
  is already there, a PAM-less index for a PAM on the 5' side of the guide (as the header of its
  .bin files records), and ... puts the PAM on the 3' side. ...`) and writes nothing, the web reads
  `(None, 23, -3)` of the index, and a 3′ search stops, saying the header records that side; with a
  5′ PAM file, `build-index-only` records that file, the web reads its name and the 5′ geometry,
  and a 5′ search of the index ends with exit 0 and nothing on stderr. Of a 3′ index of flank 4,
  `build-index-only` with a 3′ PAM file of that flank keeps it (exit 0) and records that file.
- **The second-indel check of an indel row reads no base beyond its chromosome.** The haplotype
  window `_breaking_slots` writes around an indel target ran past the chromosome's ends by up to
  the target's length and the indels' (the positions it read from the dictionary were clamped at
  the start only); a target of an index of format 2 may end on a chromosome's last base or start
  on its first. The window now stops at both ends. No row changes: a comparison with a base
  beyond the chromosome failed already. Every other window CRISPRme reads around a row is the
  row's own bases (the scores take the aligned pair alone; the variant lookups read the row's
  span, guarded where the reference comes back short), and an indel row lies wholly in its own
  fake-indel window, the windows being joined by one `N` (checked through the real CLI above).
- **METHODS §9, `docs/PRECOMPUTED_INDEXES.md`, `docs/DYNAMIC_PAM_INDEX_COMPAT.md` and the README**
  say which limits at contig ends and `N` runs are format 1's, what format 2 holds, and that the
  equivalences of a PAM-less index with `--prune-pam` and of budgets hold between indexes built
  for the same guide length (an index built for 20-nt guides holds sites at contig ends and `N`
  runs that one built for 21-nt guides does not).

<!-- fmt track: [2.9.0] entries for CHANGELOG.md (CRISPRme reads index format 2; the repin to
CRISPRitz 2.10.0), to be merged with the other tracks' fragments. -->

**For the merge, not for the changelog: two notes.**

1. The Changed entry on the repin to CRISPRitz 2.10.0 comes with the last commit of the track,
   which sets the version; it is to be merged only once CRISPRitz v2.10.0 is tagged. Until then
   the image and `install_from_source.sh` build CRISPRitz 2.9.1, which reads index format 1 only,
   and every search of an index of format 2 stops with the ERROR the first Added entry quotes.
2. [2.8.1]'s entries that say CRISPRme 2.8.1 refuses an index of format 2 ("the CRISPRitz planned
   as 2.9.2", "index format 1 only") stay as they are: they say what 2.8.1 does.

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
  Through the real CLI with CRISPRitz `feat/engine-boundaries` (format 2) on a synthetic genome
  of perfect, mismatched and bulged sites at contig ends and beside `N` runs, with SNPs and
  indels by them (`test/index_format2/check_boundary_rows.py`): every search -- reference, and
  variant per-sample and population-level, NGG and the 5′ TTTV, and NGG under `--prune-pam` off an
  8-nt flank -- ends with nothing on stderr; every row is its genome's bases at its coordinates;
  the results hold every row of the same searches with CRISPRitz 2.9.1 (13 against 6 NGG rows, 2
  against 0 TTTV rows, each row more at a boundary) but the one 2.9.1 reports on the wrong strand
  (a − strand PAM at a contig's first base, IDX-14); and the 8-nt flank gives the dedicated
  index's 13 rows (with 2.9.1, 5 of its 6). On hg38 chr22 (five guides -- `ACTGAAATCTGTAAGCAGGC`,
  TRAC, EMX1, VEGFA site 3 and one whose perfect match starts right after an `N` run -- NGG, 4
  mismatches + 1 DNA + 1 RNA bulge, `build-index-only` then per-sample over 300 HGDP + 1000
  Genomes samples, and a reference search): with CRISPRitz 2.9.1 this release gives 2.8.1's
  integrated results byte for byte (10,942 and 9,210 rows); with format 2, every one of those rows
  and one more, that perfect match at chr22:11,118,987, which no index of format 1 holds; a
  reference search of the format-1 index by CRISPRitz 2.10.0 gives the 2.9.1 rows and the note
  below; and one of the format-2 index by 2.9.1 stops with that ERROR, exit 1.
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
  CRISPRitz `feat/engine-boundaries`, 224 of 224 checks over 60 random contigs, 7 PAMs (NGG, NRG,
  NNGRRT, the 5′ TTTV, PAM-less K = 3, 4 and 8) and budgets 0 to 3, every count of a contig's sites
  equal to the targets the engine wrote and all 45,304 cuts of its files told; with 2.9.1, 196 of
  196.
- **A search takes an index the installed CRISPRitz reads ahead of one it would stop on.** Of a
  library holding an index of format 2 and one of format 1 that serve a search, the scan took the
  format-2 one when its budget was the smaller, and then stopped on it under CRISPRitz 2.9.1
  (PKG1-scan-takes-format2); it took the format-1 one before 2.8.1, misreading the other. The
  search script's scan and complete-search's now rank an index the installed CRISPRitz reads
  first, and take the other only where none serves.
- **A search of an index of format 1 by CRISPRitz 2.10.0 says what it lacks.** On stdout: `NOTE:
  <index> is in index format 1 (built by CRISPRitz 2.9.1 or earlier): it lacks the sites whose
  guide lies within a few bases of a contig end or an N run, which an index built by CRISPRitz
  2.10.0 or later holds (METHODS section 9). Build it again with crisprme.py build-index-only, or
  download it again, to search them.` The note `--prune-pam` gives of a flank wider than its PAM,
  and the warning of a flank built before CRISPRitz kept an `N` in it, are given of an index of
  format 1 only.
- **The side of an index of format 2 is the one its header records.** The search script,
  complete-search, complete-test, the coverage check and the web read it there before
  `.pam_build`, which an index built by `crispritz.py` itself, or downloaded, may lack.
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

# Changelog

All notable changes to CRISPRme are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). From 2.10.0
the version numbers follow the version policy of `docs/RELEASING.md`, not Semantic Versioning:
CRISPRme shares its major.minor with CRISPRitz, and a breaking change in either package bumps
the minor of both.

When cutting a release, move items out of `[Unreleased]` into a new dated
version section and update the link-reference footer. See `docs/RELEASING.md`
and the `release-crisprme` skill.

## [Unreleased]

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
  bulge, per-sample over 300 HGDP + 1000 Genomes samples, and a reference search), the indexes
  CRISPRme 2.8.1 built with CRISPRitz 2.9.1 (format 1) give this release's 10,927 per-sample and
  9,210 reference rows with CRISPRitz 2.9.1, and the same rows, byte for byte, with CRISPRitz 2.10.0
  (`feat/2.10.0` 881f761's build), which adds the note below; built again by `build-index-only`
  with 2.10.0 (format 2), they give every one of those rows and one more, that perfect match at
  chr22:11,118,987, which no index of format 1 holds (10,928 and 9,211 rows); and a search of them
  by 2.9.1 stops with that ERROR, exit 1. The reference rows are 2.8.1's. Of 2.8.1's 10,942
  per-sample rows, 10,752 are this release's byte for byte, and the other 190 differ as entries
  under Fixed and Changed say: 54 started on a base an insertion adds, 35 of which start at its
  anchor (as 34 rows: two were copies of one target inside one insertion) and 19 are now an
  alternative alignment of a row there; 16 are gone or differ by the carriers whose haplotype
  holds an indel in the target, whom a SNP row no longer names (5 and 11); 31 more, within 25
  bases of an indel row, change with the site the merge now gives that row (4 gone, 27 differ; 14
  rows there are new); 16 number an indel's variants along its target or give `PAM_creation` by
  its placement rule (13 and 3); and 73 differ only in the order of their rsID and AF lists (66), a
  dense window's note (3; Changed) or a CRISPR-Bulge score's last digit (4). Such a score can
  differ in its last digit from one batch of the scorer, or one machine, to another: so do 7
  per-sample rows and 1 reference row of the indexes built again.
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

- **Q, the code of a variant position whose reference base and SNP ALTs hold all four bases, is read
  as every base wherever CRISPRme reads a variant code; N stays a gap that matches nothing.**
  CRISPRitz 2.10.0 writes Q at such a position (decision B1), where 2.9.1 coded the reference base
  and the ALT of the position's last split record alone or, for one VCF line listing three SNP ALTs,
  wrote N and lost every target through the position. The SNP post-analysis decomposes a Q column
  into each ALT's carriers and loads the dictionary keys of a target whose only code is Q (it took
  such a target for a reference row, with Q in its DNA column and no carrier); the PAM filter of a
  PAM-less or wider index keeps a target with Q in its PAM; the SNP overlay copies Q onto the
  fake-indel genome's flanks and the indel path resolves it per carrier; the site rule CRISPRme
  checks an index's trees against (`crispritz_bins`) reads a FASTA holding Q and matches it to every
  PAM base; and the population-level lower bound (`twopass_find`), the cluster sort and the guide
  extraction from a sequence take it. A Q stands for every base, so it cannot say which position it
  codes: at a bulge boundary, where the target's column can be read against a neighbouring position
  (#139), a Q column admits that position's ALTs only where they and its reference base hold all
  four bases, an ALT that no sample carries counted, as CRISPRitz codes it (a sites-only dataset's
  admits every ALT its registry lists, as each of its columns does: see Fixed), and of them only
  the bases: an ALT that is no base, such as the `*` of a spanning deletion (in a GATK joint call),
  which the dictionary of CRISPRitz 2.10.0 lists, is no allele of a Q, as of any other code (with
  the dictionary, a row substituted it, `*` in its DNA column, and the CRISPR-Bulge scorer, which
  cannot read it, failed the batch it was scored in: no alignment of that batch, up to 100,000 the
  SNP post-analysis scores together, had a CRISPR-Bulge score, and no row of the review's synthetic
  search had one). The dictionary lists such an ALT; the Tier-0 registry holds no record of it, and
  lists the positions of four bases in its manifest (`four_base_positions`), so that the tiers
  alone, the layout of every published index, read such a position as the dictionary does. Checked
  through the real CLI on an HPRC kit (hg38 chr22:1-17,950,000 and the normalised HPRC records up to
  17,949,000; TAGGGAAGTTAGGGAGGTTA, NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per sample): of its 17
  positions of four bases, chr22:17,877,081 holds T>A (8 copies), T>C (none) and T>G (46). On the
  tiers alone no row named it, and the table held 432,738 rows against the dictionary's 439,108; it
  is now the dictionary's, byte for byte, T>A on 648 rows (8 carriers) and T>G on 6,220 (44). Each
  of the search's 8,129 Q columns, 8,108 of them in an alignment with a bulge, is read against a
  position of four bases. Covered by `PostProcess/test_site_alleles.py`, with the real
  `new_simple_analysis.py` on tiers and dictionaries the repository's own compilers build.

- **The carrier gate: a population-level row's variant combination is checked against phased
  haplotypes, and the report ranks a combination no haplotype carries after every other row, never
  dropping it.** A population-level row takes, per variant window, the alleles that match the guide
  best and lists every carrier of any of them, so two or more of them can sit on no one haplotype:
  no individual has the sequence the row reports, and, being the best-matching, it ranked first.
  After writing the integrated results the search now reads each row's variants (each of its three
  alignments) in the dataset's own phased genotypes, SNPs and indels
  (`Dictionaries/genotypes_<vcf>/`, `indel_genotypes_<vcf>/`;
  `PostProcess/haplotype_prioritizer.py`), and gives each a tier: **A**, two or more variants one
  phased haplotype carries together, with no indel of the store on that slot that takes the
  target's bases off it, as a per-sample SNP row reads its carriers (Fixed): one whose REF holds a
  variant's base at another position does, one at a variant's own position, a phasing artefact,
  does not (states `one haplotype`, and `phase not known` where only unphased calls, or two ALTs of
  one position on one slot, show it), or that a haplotype reaches
  the alignment's CFD without (states `same score`, `same score, phase not known` and `same score
  as the reference`): an allele the worst case takes can leave the CFD the row reports (its three
  decimals), the alignment's mismatches and its PAM as they are -- one at a PAM base the search's
  PAM leaves free, a mismatch for a mismatch the
  CFD table scores alike -- and then a haplotype carrying the rest, and nothing else inside the
  alignment that changes them, has the row's off-target at the row's score (the companion's
  `Scoring_variants` column names the variants the score needs). The guide's side of an alignment
  writes its PAM as `N`, which every base matches, so the PAM region is read against the search's
  own PAM (the search gives the gate its PAM file): without an allele the PAM needs there is no
  off-target at the site at all. The CFD gives this state to the alignments it scores (the
  highest-CFD one, which the report ranks by, and the fewest-mismatches one), not to the one ranked
  by CRISPR-Bulge score, which the gate does not read (the model is not run again) and which an
  allele the CFD never reads does change; **B**, one variant; **C1**,
  two or more variants, each in the panel, that none of its haplotypes carries together nor
  reaches the row's CFD without, with E,
  the haplotypes expected to carry them together if they were independent (below 3, the state is
  `too rare to tell`), or variants, a single one too, every haplotype carrying which has an indel
  that takes the target's bases off it (`indel in the target`); **C2**, a variant not in the panel,
  which cannot say (the frequency bound, the product of the variants' frequencies >= 1e-06 and the
  rarest >= 0.001, is shown and applied to nothing). The gate changes no row: it writes
  `..._carrier_gate.tsv` beside the integrated results (one line per alignment that uses a variant)
  and `.carrier_gate.json`. The report ranks a C1 row after every other row in the views ranked by
  the score of the alignment that is C1: a row whose highest-CFD alignment is C1 leaves the
  top-1000 by CFD, panels (a) and (b) of Section 2, the validation panel's CFD thresholds and the
  risk card's CFD bands, and is listed in a table of its own after the top-1000
  (`unobserved_worst_case.tsv` holds them all; its caption
  counts the perfect matches among them); a row whose highest-CRISPR-Bulge alignment is C1 leaves
  the top-1000 by that score, panels (c) and (d), and the validation panel's CRISPR-Bulge
  thresholds, as the risk card's
  CRISPR-Bulge bands. The two alignments of a row need not be the same one, so a row of the
  demoted table can keep its place in the views ranked by CRISPR-Bulge score, and the other way
  round. The panel's totals (its off-targets, perfect matches and variant-created
  off-targets) and `run_manifest.json`'s counts are of every row, as the Section-1 matrix's, and a
  row of the panel's table says how many rows its thresholds, the panel and the risk card leave
  out. A, B and C2 rows keep their place, and the `Observed` column says what the gate found (see
  Changed). The inputs card has a **Carrier gate** row: its panels, its counts by tier, and the
  panel's samples by superpopulation (AFR 986, AMR 549, CSA 782, EAS 816, EUR 771, MID 157, OCE 30
  on HGDP + 1000 Genomes), as a haplotype found only in an ancestry with few samples there can be
  missing from it. On by default for a population-level search of a genotyped dataset;
  `--no-prioritize-observed` turns it off, `--gate-min-expected` sets the `too rare to tell`
  threshold (default 3; 0 turns it off). A per-sample search, whose rows name their cis carriers,
  needs none; a dataset with the classic dictionaries alone has no phased store to read, and the
  search and the report say so. On the genome-wide population-level search of the README guide
  (HGDP + 1000 Genomes, NRG, 4 mismatches + 1 DNA + 1 RNA bulge, 198,812 rows) the gate took 29 s in
  the search (26 s and 1.6 GB on 24 threads alone) and gave, by each row's highest-CFD alignment, A
  4,153 (689 `one haplotype`, 1 `phase not known`, 2,574 `same score`, none `same score as the
  reference`), B 55,180, C1 40,231 (246 where 3 or more haplotypes were expected, 39,780 `too rare
  to tell`, 205 `indel in the target`, 173 of them of one variant) and C2 0, its reading of an indel
  over a variant's base (above) moving chr5:4507931 (+), whose 58 haplotypes carrying chr5:4507935
  T>A all hold chr5:4507932 CTTTCA>C, from B to `indel in the target`. The report's top-1000 by CFD
  and validation panel then hold no row whose highest-CFD alignment is C1, and its own table all
  40,231, while the panel's totals are the search's (198,812 sites, 1 perfect match, 91,275
  variant-created off-targets) and 43,264 rows are C1 at their highest-CRISPR-Bulge alignment. The
  integrated results are byte for byte those of the same search with `--no-prioritize-observed`
  (TRAC's search of chr22, measured so). Checked against all 8,182 phased haplotypes by a harness
  that rebuilds every alignment's window on each haplotype from the VCF, on 2.8.1's run of the same
  search (198,826 rows) with this release's gate (A 4,155, B 55,190, C1 40,232, C2 0 by each row's
  highest-CFD alignment, in 27 s on 24 threads, the integrated results unchanged; the gate of the
  integration's first stage, before the fixes of 2.8.2 and that reading of an indel, gave B 55,191
  and C1 40,231): the only alignments left in tier A or B that list a combination no haplotype
  carries are the 2,574 of a `same score` state (2,328 at the fewest-mismatches alignment, none at
  the CRISPR-Bulge one, as on this release's run), which is what that state says; every other A or
  B alignment lists a combination a haplotype carries, and every C1 alignment that a haplotype does
  carry is one with an indel in the target (205 at the highest-CFD alignment, 187 at each other).
  The harness's stricter test (a haplotype that reproduces the alignment within its mismatches with
  a valid PAM) reaches 2,572 of the 2,574 (2,326 of the 2,328), and the per-sample run of the same
  search reports a row within three bases of every one of them, at the CFD the row shows. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **`--prioritize-observed` checks a sites-only dataset's rows against the phased HGDP + 1000
  Genomes genotypes, and `download --what genotype-oracle` installs them.** The `mega` panel has no
  genotypes, so its multi-variant rows could not be told apart. With `--prioritize-observed` (the
  web form's check box, offered for a sites-only panel once they are installed) the gate reads its
  rows in the phased HGDP + 1000 Genomes genotypes and indel genotypes (the genotype oracle), which
  an installed `NRG_3_hg38+hg38_HGDP1kGP` index already has, which
  `download --what genotype-oracle [--panel HGDP1kGP]` installs without the index (4 min and 28 GB
  on disk from the Hub), and which `CRISPRME_GENOTYPE_TIER` can name elsewhere. A variant the panel
  does not hold makes its row C2, kept in place, with its frequency bound (a SNP's frequency from
  the dataset's registry, an indel's from its indel frequencies, `AF_max` both). On chr22 of the
  mega panel (five guides, NRG, 4 + 1 + 1, 20,010 rows): A 138 (32 of them `same score`), B 4,318
  (1,479 of a variant not in the panel), C1 587 (58 where 3 or more haplotypes were expected, 513
  `too rare to tell`, 16 `indel in the target`, 13 of them of one variant), C2 3,297 (1,833 within
  the bound, 1,464 below it), the integrated results byte for byte those of the search without the
  flag. Without the
  oracle the search says so (a `WARNING` line on stdout), and the rows are reported unchecked.
  Until the Hub has the panel's phased indel genotypes (see Changed), the download says so in one
  line and installs the SNP genotypes, and an indel inside a target then takes no haplotype off.
  Covered by `PostProcess/test_haplotype_prioritizer.py` and
  `PostProcess/test_web_search_reuse.py`.

- **A window that spells 2^20 sequences or more is flagged hyperdense, whatever `CRISPRME_IUPAC_CAP`
  is.** A window's sequences are the product of the bases each letter of the variant genome stands
  for: 2 for a two-allele column (R, Y, S, W, K, M), 3 for a three-allele one (B, D, H, V), 4 for a
  four-allele one (`Q`), 1 for a base or an `N` (an assembly gap, which matches nothing). Such a
  window matches nearly any guide somewhere, as uninformative as a gap, yet the search reports it,
  as it reports any window it aligns within the mismatches and bulges asked for. It is listed in
  `<out>.high_variant_density_regions.bed`, read from its own letters (a window whose only variant
  letters are four-allele positions too), and the report's `High_complexity_region` column gives its
  count (`24 in window (2^24.0 sequences)`), as the integrated results'
  `High_variant_density_region` note does (`2^24.00 sequences`). The BED gains three columns, `bits`
  (log2 of the window's sequences), `four_allele` (how many four-allele positions it holds) and
  `greedy` (1 when a greedy worst-case alignment stands for the window, as in any
  population-level search, a sites-only dataset's among them), and its
  `n_variants` counts every variant column of the window, a four-allele one too. Only where one
  stands for the window does its note say that a greedy worst-case alignment is reported there:
  within a cap raised to its bits, or with `CRISPRME_IUPAC_CAP=-1`, a `--per-sample` search
  enumerates a hyperdense window's combinations as any other's, and its rows read as in 2.8.1
  (`N cis carriers`). The report's inputs card says how many rows lie in such a window, and
  `complete-search` says per chromosome how many windows were dense (over `CRISPRME_IUPAC_CAP`) and
  how many hyperdense (stdout, `log_verbose.txt`). Among the indexable NRG windows of chr22, the
  reference genome has none, the variant genome of HGDP + 1000 Genomes published with CRISPRitz
  2.9.1 has 1,763 (the worst 2^23.00), and the one CRISPRitz 2.10.0 writes from the same VCF, every
  allele coded, 5,103 (the worst 2^44.75, at chr22:22,881,189 in the IGL locus). On the published
  index, the README guide's search of chr22 lists 565 dense windows, 381 of them hyperdense (the
  worst 2^24.00; each aligned with 2 to 4 columns that do not admit the guide's base), and 164 of
  the population-level report's 2,197 rows lie in one (1 of the per-sample report's 1,570). The
  threshold is one constant, `window_density.HYPERDENSE_BITS`; see
  [`METHODS.md` §5](METHODS.md#5-search-space-control-for-high-variant-density-regions).
- **A row whose window holds a four-allele position says so.** Where the reference base and three
  ALTs cover A, C, G and T, the variant genome writes `Q`, which matches any base (`N` stays an
  assembly gap). The post-analysis lists the four-allele positions of the windows the search
  reported in `<out>.four_allele_positions.bed` (0-based; written when there is one) and says how
  many on stdout; `High_variant_density_region` names those in a row's window
  (`four_allele_positions=<chrom>:<position>`, 1-based), the report's `High_complexity_region`
  column reads `1 four-allele position`, and its inputs card says how many rows have one. A
  four-allele position alone makes no window dense to the report.
- **`complete-search --max-degeneracy K`, off by default: a reported SNP row rests on at most K
  ambiguous columns unless its carriers hold its variants together on one phased haplotype.** An
  ambiguous column a row rests on is one of the variant genome where the reference base would be a
  mismatch, or break the PAM, and the row's allele fits. Few alignments resting on two or more are
  carried by anyone: on chr22 of HGDP + 1000 Genomes (the published variant genome; five guides,
  NRG, 4 mismatches + 1 DNA + 1 RNA bulge) a real phased haplotype carries 96.8% of the alignments
  resting on one, 5.6% of those resting on two, 1.1% on three and 0.02% on four or more. With the
  cap, a worst-case alignment (population-level, or a dense window's greedy one) is cut back to the
  strongest within it (a valid PAM first, then the fewest mismatches or the highest CFD), and any
  other row resting on more is left out, unless its carriers hold its variants together that way (a
  CONFIRMED observed haplotype, or a phased dictionary's combination). The reference row keeps the
  locus where the reference alignment is within the search's limits; a site that only such
  alignments reach goes. `complete-search` says per chromosome how many alignments were cut back and
  how many rows left out, writes the cap to `Params.txt` (`Max_degeneracy`), and the report's search
  criteria state it. An indel row is resolved per carrier haplotype, as before, and is not cut. On
  chr22 of the published index, the README guide's population-level search with `--max-degeneracy 1`
  cut 20,918 worst-case alignments back and reports 1,572 rows instead of 2,197, every SNP variant
  row then resting on at most one column (of the 1,257 such rows without it, 681 rested on two or
  more): 625 sites go, those only such alignments reach (9 more move by 1 or 2 nt, to where their
  cut-back alignment starts), the rows lying in a dense window fall from 281 to 2, and 47 to 49 of
  the 100 sites with the highest CFD leave the top 100, as the ties at the 100th CFD fall (the
  second and third, CFD 0.299 at chr22:22358032 in IGL and 0.275 at chr22:49464611, rest on 11 and 2
  columns). The per-sample search changes no row: its multi-variant rows are CONFIRMED haplotypes.

- **A dedicated index serves the PAMs its own holds.** An `NRR` index stores every site of `NGG`,
  `NAG`, `NGA`, `NAA` and `NRG`, an `NRG` index every site of `NGG` and `NAG`, as concrete bases in
  their stored flank; under `--prune-pam` the engine writes out only those of the search's PAM, at
  leaf time, as it does an all-`N` index's. A search takes such an index after its PAM's own and
  before a PAM-less one (of two, the narrower PAM: `NRG` before `NRR` for `NGG`), and only for a
  PAM on the side of the guide it was built for. It was refused with `--index-path`, and a search
  without it built an index of its PAM (20 to 35 minutes for hg38; a variant search stopped). One
  rule, `index_fit.pam_serves` and the search script's `_pam_holds`, in the search script,
  `crisprme.py`, `unfinished_marks` and the web form. Through the real engine (the release gate's
  `test/index_flex/check_index_flex.py`): `NGG` from an `NRG` or an `NRR` index, and `NRG` from an
  `NRR` index, give the same results as a dedicated index of the PAM. On chr22 with HGDP + 1000
  Genomes (two guides, 4 mismatches + 1 + 1 bulge, per-sample, 16 threads), an `NGG` search that
  the `NRR` index alone served, pruned by itself, wrote the same 2,023 rows, byte for byte, as one
  of the dedicated `NGG` index, in 397 s against 295 s, at 9.6 GB peak memory against 2.6 GB
  (CRISPRitz 2.9.1); with CRISPRitz 2.10.0 both wrote those rows again.
- **A search prunes by itself.** When no installed index serves a search unpruned and one serves it
  pruned (a dedicated index of a PAM that holds the search's, or an all-`N` one wider than the PAM),
  the search prunes, as if `--prune-pam` was given, and says so first, naming the index that needs
  it: from the command line on the console (`NOTE [complete-search]: no installed variant index of
  hg38+hg38_HGDP1kGP serves PAM NGG unpruned, and NRR_3_hg38+hg38_HGDP1kGP (...) serves it pruned
  to it at leaf time: searching with --prune-pam, as if it was given. A dedicated NGG variant index
  is searched faster: ...`; in a library of the PAM's own variant index and another PAM's
  reference index, `no installed reference index of hg38`), from the web, which has no console, in
  the search's `log_verbose.txt` (`NOTE: no installed ...`). The web form, which never passes the
  flag, gets the same: the search script decides it for the search it runs, and the form's bulge
  cap and PAM list count those indexes. Of an all-`N` index wider than the PAM that it reads so, a
  web search also says in its log what such a flank cannot index at a contig's ends, and warns of
  one built before CRISPRitz 2.9.1, as complete-search says both on its console and in `.Params.txt`
  (one text, `prune_flank.wider_flank_notes`; the web search's warning names Settings beside
  `build-index-only`). Before, a search that only a K = 8 index served stopped for want of a variant
  index (or, for the reference index, built one), and said to pass the flag. Where an index serves
  the search unpruned, as an all-`N` one of the PAM's width does, nothing changes. On chr22 with
  HGDP + 1000 Genomes (CRISPRitz 2.10.0), an `NGG` search whose variant index was `NGG`'s own and
  whose reference index only the `NRR` one, under its older name or as `20-NRR_3_hg38` with no
  `.pam_build`, said `no installed reference index of hg38`, and wrote the same 2,023 rows, byte for
  byte, as the dedicated `NGG` indexes. On chr22 with HGDP + 1000 Genomes and the `NRR` indexes
  alone (CRISPRitz 2.9.1), an `NGG` search with one bulge of each type said that it pruned by itself
  on the console from the command line, and wrote nothing of it in its `log_verbose.txt`; the same
  search sent with the web form's Submit (`change_url`) said it in its `log_verbose.txt`. Both
  recorded `Prune_pam auto` and wrote the same 2,023 rows, byte for byte, with nothing on stderr. On
  a 30 kb genome whose `20-NNNN_2_gn` was built before an N run was put into it, an `NGG` search
  from the web form says in its log what the 4-nt flank does not index and that the index was built
  before CRISPRitz 2.9.1, and from the web form of the code before, only that it pruned. Covered by
  `PostProcess/test_defaults_and_index_flex.py` and `test_prune_flank.py`.
- **The run and its report name the index that served the search.** `.Params.txt`'s `Genome_idx`
  says what the search asked for, which another index may serve (`NGG_2_hg38` served by
  `NNNNNNNN_3_hg38`), or none at all (IDX-21). The search script now records the indexes it took,
  in `.indexes_searched.txt` in the result folder, beside the parameters (which the web compares
  byte for byte to reuse an earlier job, so they stay as the request wrote them): `Reference_index`
  and `Variant_index:<dataset>`, each with its library, the PAM it was built for and its side of
  the guide, the guide length its trees were built for, its bulges and how the search narrows it
  to the search's PAM, the variant one with its INDELS index, or `none` for a search with no
  bulges, which scans the genome (and the variant genome, where it is here whole: the record then
  names neither index, whatever is installed, only the INDELS index the indel search reads); and
  `Prune_pam auto` for a search that pruned by itself, whichever decided it — complete-search,
  which tells the search script that it took the flag (it records `Prune_pam`, `auto` or `given`,
  in the parameters too), or the script, for a web search. The report's
  summary shows them, and its zip bundles the record as `indexes_searched.txt`; a run of an
  earlier release has none. On chr22 with HGDP + 1000 Genomes, an `NGG` search of the `NRR`
  indexes with 1 + 1 bulge and no flag records both index lines, each `pruned to NGG at leaf
  time (--prune-pam)`, and `Prune_pam auto` with them, which only a web search wrote before.
- **The web form's PAM picker for a PAM-less index** (plan 2.7, DYNAMIC_PAM_INDEX_COMPAT §8). When an
  all-`N` index is installed for the genome (or the dataset) selected, the form says what it serves
  (any PAM of its guide length on its side, up to its flank) and offers the PAMs of common nucleases
  it serves (SpCas9 NRR, recommended for IND, then NRG and NGG, SpG, SpRY, VQR, VRER, xCas9, ...),
  those it does not serve disabled with the reason, and a field for a PAM of one's own, checked when
  the field is left. "Use this PAM" writes that PAM's file into `PAMs/` and selects it; the search
  prunes the index to it by itself. Only a PAM the index serves, and a nuclease of the table, is
  written. A server online writes nothing into the data folder, as its Settings page does not:
  there the picker selects the PAM file installed for the PAM, whatever its name, offers only the
  nucleases whose PAM file is installed, and says that it adds none for any other. Covered by
  `PostProcess/test_defaults_and_index_flex.py` and `PostProcess/test_web_bulge_notes.py` (in the
  image).
- **`scripts/check_indel_rows.py`: an independent check of the indel path on real data.**
  `validate-test` compares only the engine's **sites** — (chromosome, start, strand) — on
  the reference and on the **SNV-only** 1000 Genomes chr22 genome, so an indel row's
  coordinate, its aligned spacer/target strings, its REF/ALT text, its variants, its
  carriers and the whole alternative-alignments table were never checked against an
  independent truth on real data, and the synthetic SNP+indel oracle checks them only on a
  planted synthetic chromosome. Neither could catch the two defects behind the 2.8.0 row
  reported at chr22:22322768 (a phantom ALT pairing an insertion's only carrier with SNP
  alleles other individuals carry, fixed in 2.8.1; and a row whose first aligned base is
  one an insertion adds taking that base's place on the haplotype as its coordinate, fixed
  below).
  - It imports only Python's standard library and **no code of CRISPRme's variant or indel
    machinery, and copies none**: a unit test compares its code lines with those of every
    Python file of the repository's tree, in every folder (158 files: PostProcess and its
    `supportFunctions` and `fixtures`, `scripts`, `seq_script`, `pages`, the oracle and
    benchmark scripts of `test/`, the paper's scripts, `crisprme.py`), the unit tests aside,
    and no run of four is the same (a run of imports aside; its rule for a line of several
    ALTs had been `site_alleles`' narrowing loop copied line for line, which a test of
    function names let through, and is now written from METHODS' text). Its inputs are what
    a user has — the reference FASTA, the VCF(s) the search read (its own BGZF + tabix
    reader; without a `.tbi` it indexes each
    VCF once in memory, writing nothing), the samples list (`--samples`, the samplesID files
    the search read; by default the run's `.sampleID.txt`, else every VCF sample), the
    guide (`--guide`, required: the table's own column is never trusted), the PAM file, the
    search budget and the result folder; the CFD tables in `PostProcess/` are read as data.
    It writes each sample's haplotypes itself from the VCF 4.x definition of a record
    (every SNV, MNV and indel of a slot), keeps its own map from haplotype bases to
    reference coordinates and aligns with its own exhaustive dynamic-programming aligner. A
    row's variants are read with REF and ALT in upper case, as every record, as VCF 4.x reads
    them (HPRC writes 1,546 indel records of chr22 in lower case, which a release before this
    one names as written and this one in upper case (Fixed); the check looked their rows up
    by their raw tokens: a world of a lower-case deletion gave 21 gating `missing-row`, every
    carrier of every target).
  - **The truth is the VCF.** A slot holds a row's target when it reads the target's bases
    from the row's start through the row's record: each base the record puts in the target
    at its place, and the base after them. Two alleles of one slot whose REF spans meet give
    no sequence; where every allele that meets the record's REF is one METHODS section 4
    reads the slot without (another allele at the record's own position, an allele only on
    bases it deletes, an MNV) the slot is read without them, and a listed carrier must hold
    the target there (`overlap-reading-agrees`, a note, when it does); any other overlap of
    the record leaves the slot no reading, and two other alleles of the slot that overlap
    each other are tried on every reading that can change what it says of the target. A slot
    with no reading gives no truth, so METHODS' rules are the only model there: a listed
    carrier they take off the target for a documented cause fails (`carrier-off-by-rule`),
    any other is counted. Where a call touching the target has unknown phase (unphased, a
    missing allele, or fewer copies than the record's call), or the record's slot holds two
    SNV ALTs of one position that the search's stores both read (every allele coded: no
    sequence, and CRISPRme resolves such a carrier PUTATIVE, decision B4), what the record's
    known slot says comes first, and otherwise the truth is METHODS section 4's PUTATIVE rule: the
    strongest target the sample's own alleles allow on the row's alignment — the guide's
    base wherever it carries it and the alignment pairs the position with a guide base;
    elsewhere (a DNA bulge included), and at the PAM jointly, the alleles that maximize the
    alignment's CFD (a CFD tie left open); a search that scores no CFD (a 5′ PAM, or a guide
    or PAM not 20 and 3 nt long) the reference base where the alleles and the motif allow
    it, else the first of A, C, G, T they allow —, never an allele it does not carry, its
    copies counted against the record's call. The rule chooses only at the fake-indel
    window's code columns, the target's reference bases off the record's own REF span: the
    window writes the ALT verbatim, the anchor included, so a SNV at the record's position
    is never a choice. The records are those the search read: CRISPRitz's enricher skips
    every VCF line whose FILTER is neither `PASS` nor `.` (no such line is in the enriched
    genome, its dictionaries and stores, the indel log or the fake-indel windows), and so,
    from this release, does the phased indel store the carrier resolution reads (Fixed,
    below); a store written before kept every indel whatever its FILTER (a second or
    overlapping indel took a carrier off a target, and a SNP row dropped a carrier for an
    indel of its slot), and a search that reads one says so in its log. So a genotyped
    dataset's line outside `PASS` and `.` is read as the reference, but for its indels where
    the run read a store written before 2.10.0, which stay on their carriers' haplotypes;
    either way no record of the search holds it (a row naming one fails, `variant-filtered`;
    completeness demands none of their targets, `limit-filtered-record`); a sites-only
    dataset's keeps its SNVs for its SNP rows (the `hg38_mega` index's hold them: its
    registry holds all 78,116 of chr22 outside `PASS`, and its SNP rows name some), a SNP
    row through one being a row CRISPRme may write or not, while its fake-indel windows
    code its `PASS` SNVs only, every allele of a split site on a genome CRISPRitz 2.10.0
    built and the last record's on one of 2.9.1 or earlier (as CRISPRitz's overlay does: a
    carrier-less row naming another fails, `variant-filtered`, as does one naming an SNV of
    another dataset, `variant-not-in-dataset`), and it drops its indels and MNVs
    (`hg38_mega`'s indel log holds exactly its `PASS` indels).
    CRISPRme's documented reading of a carrier (METHODS section 4 and the
    2.8.1 Known limits: one indel per window, resolved with each carrier's own SNP alleles;
    a second indel read from either end of the target through the indel; the phasing
    artefact and the complex allele; MNVs and untied multi-ALT alts resolved PUTATIVE (an
    alt's number the one every carrier of it holds, less those the position's other alts of
    one character of a line of several ALTs settle, a `*` among them, an indel's ALT not, an
    alt whose carriers hold only 0 and 1 settling none), and from this release a slot
    holding two ALTs of one position; on a genome of CRISPRitz 2.9.1 or earlier only the
    last record of a split multi-allelic SNP site searched, an earlier record of the same
    REF>ALT too, whose carriers its dictionary lost, and every line before a last line of a
    `*` alone, whose entry 2.9.1's dictionary kept while its genome kept the earlier lines'
    code; on one of 2.10.0 every allele, two listings of one SNV REF>ALT at a position one
    allele, on a line of several ALTs too, their calls merged by the rule CRISPRitz's own
    documentation gives its enricher) is re-implemented from that text and only classifies
    a disagreement: a documented limit (non-gating, listed with its cause, a PUTATIVE tie
    included) or a failure. Two lines of one indel REF>ALT at a position are one allele on
    every release, in one case or in two, their calls merged by that rule (decision B2, as
    2.10.0 reads them: Fixed): the check read the record's call off the first line that
    listed the sample, so a slot only the later line gave was no slot of the record, and a
    listing on both a second allele at the record's position (a world of a deletion on two
    lines: right rows failed, `carrier-without-variants`).
  - **Releases**: where CRISPRme's and CRISPRitz's releases read the data otherwise, the
    check reads which one made the run from the run's own records and prints the reading
    it took (`read as a run of CRISPRme 2.10.0 or later (<dataset> on <chrom>: every allele,
    store of PASS and '.', INDELS format 2)`). A run of 2.10.0 writes `.indexes_searched.txt`
    (and `.carrier_gate.json`); one without them is an earlier CRISPRme's (CRISPRitz 2.9.1
    or earlier: index format 1, the last record of a split site, a store of every FILTER).
    Of a run of 2.10.0, the variant genome codes every allele where its variant index's
    `.bin` headers flag Q, or its enriched genome in the working directory codes an earlier
    record's ALT at a split site, and the last record's where that index is of format 1, or
    the enriched genome holds the last ALT alone, or an index of format 2 flags no Q while
    the VCF has a position of four bases (a line whose ALT is a `*` is no record of such a
    site); the INDELS index's format is its `.bin` files', else the log's format-1 NOTE of
    it; the store applies the FILTER rule but where the log notes a store written before. A
    run whose records do not tell, where its data makes it matter (a split SNP site, two or
    more lines of an SNV ALT at one position, a line of several ALTs among them, or such a
    line before a last line of a `*` alone; a target of
    the guide, PAM and budget that can lie at a fake-indel window's edge; a
    genotyped VCF with an indel record of another FILTER and no log), is refused (exit 2),
    as is a run of every allele whose VCF holds, near what it checks, two lines that each
    list several ALTs at one position (their allele numbers are not told apart, METHODS
    section 4), and a run of any release whose VCF holds there an SNV whose ALT is written in
    lower case (CRISPRitz codes no lower-case ALT, while its dictionary lists the carriers:
    read in upper case, the check demanded targets the search cannot find), the check
    stopping at the first work unit that meets such a position (it refused a run of two lines
    of several ALTs only once every unit had been checked, the whole check time of a
    genome-wide run),
    and so is a run made with the SNP+indel co-occurrence off
    (`CRISPRME_INDEL_SNP=0`, which a search of this release notes in its `log_verbose.txt`):
    every SNP code of its fake-indel genomes reads as the reference, so its indel rows name
    no SNP and list every carrier of their indel, none of the per-carrier resolution the
    check models (the check of 2.8.2's code gave 166 and 77 gating findings on the
    SNP+indel oracle's world searched so, per sample and population-level). An index keeps
    a target's site, its PAM and the guide's length of bases from it, not its alignment, so
    a target is read by its site: one that runs past its fake-indel window (25 bases before
    the record's first REF base, 27 after its last), or whose site's guide region does (an
    RNA-bulged target at the window's first base), is held by no index, and its carriers are
    counted (`limit-window-reach`), as are, on an INDELS index of format 1, those of a
    target whose site's guide region ends within the budget of the window's edge
    (`limit-window-edge`), which an index of format 2 holds. A cluster reached only
    through a SNP row over the run's `--max-degeneracy` is `limit-max-degeneracy`; the
    carrier gate changes nothing the check reads.
  - **Soundness**: every indel/MNV row block of the integrated results (the three blocks)
    and of the alternative alignments must be labelled ALT, name one indel-path record (and
    SNPs), all records of the VCF the search read; have an aligned spacer that is the guide
    with its PAM, counts and a bulge type that recount, no gap in the PAM, a PAM the search
    matched and a REF
    text read from its start; and have a target that is the reference with exactly its named
    variants, read from its start under the placement convention below, holding a base of
    its record. Each listed carrier must be a sample of the panel and hold the target on a
    haplotype at that start with every named variant, or (phase unknown) be listed for its
    strongest target; each carrier whose haplotype holds it, or whose strongest it is, must
    be listed (the merge never moves a carrier off its sequence's row). An ALT-origin block
    that names no variant is reported. A sites-only dataset's rows list no carrier: a row
    that lists none for a record a sites-only VCF holds is checked but for its carriers
    (`limit-no-genotypes`), each record deciding for itself, so that beside a genotyped
    dataset (the default index with `hg38_mega`) the genotyped dataset's rows are checked
    in full. A row naming an earlier SNP record of a split site is a note
    (`variant-shadowed`).
  - **Completeness** (the real-data brute force): for every indel and MNV record in the
    region, each carrier's known slots around it are written with their known alleles (a
    target touching a call of unknown phase is demanded as the carrier's strongest),
    identical sequences aligned once, and every in-budget alignment on both strands with the
    PAM that holds a base the record changes (for a pure deletion, the bases on both sides
    of it) must be reported for that sample **by rows no worse than its best alignment** by
    the merge's own sort keys (`--sorting-criteria`, then the CFD at three decimals): rows
    that list the sample with exactly that target and variants whose fold groups (one per
    carrier set) are together no worse, or rows of its site and strand naming the same
    variants with the target's own carrier sets that are together no worse. A PUTATIVE
    co-carrier listed on some of the target's alignments only splits them between two
    carrier sets (METHODS section 4), the alignment of the fewest edits in one and that of
    the highest CFD in the other: judged one group at a time, both alignments reported for
    the sample failed (`missing-best-alignment`; the review's world WQ, two carriers).
    The site is the cluster position on the haplotype carrying the record (the Cluster
    Position plus the rank of the row's first base among the bases an insertion adds, as
    the merge compares rows that start inside one insertion). A target that needs a second
    indel-path record is demanded exactly where CRISPRme's reading of the slot keeps the
    carrier, the two-indel limit otherwise; one an overlap zone of the slot changes on some
    reading and not on another is counted, not demanded. A run that writes no alternative
    alignments (`.alt_mode` = `no-alt-alignments`: by default a population-level run, as
    the chr22 self-test's) keeps one row per merge cluster: the target's cluster — the
    integrated row with a block within `--merge` bp of its Cluster Position, or of one the
    rows CRISPRme can write chain to it: its reference and SNP rows (in the mode the run
    ran, its `.search_mode`), the targets some carrier of an indel-path record holds, and a
    sites-only dataset's worst case of each alignment of its window (the record's ALT
    written verbatim, a SNV coded only off its REF span; a record a genotyped and a
    sites-only dataset both hold gives the rows of both) — must hold a row no worse than the
    target in both blocks. Rows CRISPRme may write or not, which no documented rule settles
    (the one-record reading of an overlap METHODS does not cover, a target an overlap zone
    changes, a target of a slot read without its overlap that touches a call of unknown
    phase, a target through an SNV its stores read otherwise, a per-sample run's capped
    windows' greedy rows, a per-sample SNP target of a slot that also holds an indel
    touching it), are followed only where the rows it writes reach no row as good as the
    target: a cluster reached only through them is counted under their cause's class,
    neither passed nor failed. Records only sites-only VCFs hold are counted, one
    `limit-no-genotypes` line per chromosome, and those whose FILTER the search skips, one
    `limit-filtered-record` line.
  - **The SNP-row report** (`--snp-rows report`, per-sample runs, report-only): the same
    checks of a SNP row's strings, target and named variants (each must lie in the target),
    its listed carriers and the phased carriers a SNP row of its site leaves out, each read
    as METHODS section 4 reads a SNP row's carriers: an individual with an SNV call of
    unknown phase in the target, or with an alt of a line of several ALTs whose number the
    genotypes cannot tie, is listed for every non-empty subset of its variants there
    (`limit-snp-row-unphased`; such an alt's carrier was `carrier-without-target`); the
    merge joins the carriers of a site's SNP rows of one
    variant list (its key for a SNP row holds no alignment), so a carrier may be listed for
    another alignment of the site its slot holds; and since 2.10.0 a carrier whose slot of
    the row's SNPs holds an indel of the phased store that changes the target is dropped,
    unless its own bases still read the target, from the first base on or up to the last,
    through the row's SNPs (or that indel's call cannot be tied to the slot, or it sits at a
    named SNP's position), and one whose bases cannot be read leaves the row: listing a
    carrier the rule drops is reported (`snp-row-indel-carrier`), an MNV the store does not
    hold is `limit-snp-row-mnv`. Population-level SNP rows, and a sites-only dataset's, are
    worst-case unions, out of scope.
  - **Not checked**: the `CFD_score_REF`/`_ALT`, `CFD_risk_score` and CRISPR-Bulge columns
    and which row the CRISPR-Bulge block keeps, `Variant_info_spacer+PAM`, `PAM_creation`,
    the seed counts, `Variant_MAF`, `Variant_rsID`, `Other_motifs`, `Not_found_in_REF` and
    the annotations; a row's own `CFD_score` is not compared with its alignment (the merge
    checks recompute it from the row's aligned strings).
  - **Classes**: gating ones name the defect (`start-not-on-variant-sequence`,
    `target-not-on-variant-sequence`, `ref-text-mismatch`, `carrier-without-target`,
    `variant-missing`, `carrier-not-strongest`, `carrier-missing`, `carrier-not-in-panel`,
    `carrier-off-by-rule`, `missing-row`, `missing-putative-row`, `missing-best-alignment`,
    `missing-locus`, `origin-not-alt`, `alt-without-variant`, `row-unreadable`,
    `row-not-checked`, `variant-filtered`, `variant-not-in-dataset`, ...); non-gating ones,
    each listed with its samples, are CRISPRme's documented limits (`limit-last-alt-wins`,
    `limit-two-indels`, `limit-putative-carrier`, `limit-overlapping-records`,
    `limit-capped-window`, `limit-contig-end`, `limit-n-bases`, `limit-no-genotypes`,
    `limit-filtered-record`, `limit-window-reach`, `limit-window-edge`,
    `limit-max-degeneracy`, and in the SNP-row report `limit-snp-row-mnv` and
    `limit-snp-row-unphased`), the notes
    `variant-shadowed`, `overlap-reading-agrees` and
    `guide-not-given`, the SNP-row report's own findings (report-only: its
    `snp-row-indel-carrier` and the classes above), and five coverage gaps of the checker
    itself, not CRISPRme limits (`coverage-overlap-cap`, `coverage-putative-cap`,
    `coverage-chain-limit`, `coverage-snp-row-indel`, `coverage-unindexed-long-ref`).
  - Output: a TSV of every finding with its class, gating flag and evidence, a per-class
    summary, exit 1 on any gating finding, exit 2 on a usage error or a configuration it
    does not model (nothing is checked then): two genotyped datasets that share samples, and
    inputs other than the search's where the run's own record of them says so — its
    `.Params.txt` budget, PAM and variant datasets (`Genome_idx`; a search without bulges
    records `Genome_idx None` whatever it searched, and its datasets are then its `--vcf`
    list's lines, read in its working directory), a guide its `.guides.txt` does not hold, a
    `--samples` other than its `.sampleID.txt`, an `--analysis` other than its
    `.search_mode` —, a run whose own records contradict each other, a search of a
    genotyped dataset without its genotype store that says so (below), a run made with the
    SNP+indel co-occurrence off, one whose records do not say how its release read the data
    where that matters (Releases, above), and one that does
    not record what the check needs (whether it wrote alternative alignments
    and, where it wrote none, its mode; `--analysis` can give it). `--chrom` / `--region`,
    `--jobs N`, `--known` (findings already reported upstream turn non-gating, and a listed
    one that stops showing is reported). It models a search that read each genotyped dataset's
    genotypes (its genotype store, or classic dictionaries): an install without them
    (`download --what index --no-genotypes`) lists every carrier of an indel on a SNP+indel
    row. Such a search records it — its `log_verbose.txt` holds the SNP post-analysis's
    `Registry-only install for <chrom>` line for the dataset (`hg38_mega`'s aggregate
    registry's, marked so, is not counted), and a `--per-sample` one writes `.search_mode`
    population-level beside `--per-sample` — and given a genotyped VCF the check refuses
    the run; without that line (no `log_verbose.txt`, or a dataset its SNP post-analysis
    never read), a run made with `--population-level` or no mode flag records nothing else
    that tells it apart, and is not to be checked.
  - Measured on this release itself (its code over the 2.8.1 image with CRISPRitz 2.10.0,
    every index built by `build-index-only` with it, of format 2, each search from a clean
    environment), chr22 of HGDP + 1000 Genomes (4,091 samples), the README guide, NRG, 4
    mismatches + 1 + 1: soundness and completeness pass in both modes (both, with the
    per-sample run's SNP-row report, about 2.5 min with `--jobs 16`; limits: 8 two-indel
    targets, 4 records with overlapping alleles and the checker's coverage gap for one
    carrier's slot near two records; no split-site target, every allele being coded). So do
    NRR and the TRAC guide per sample (19 and 20 two-indel targets, 8 overlapping records
    each), SaCas9's `CACTTTGGGAGGCCGAGGCGG` with `NNGRRT` per sample at 2 mismatches without
    bulges (the 2,103 targets at a fake-indel window's first base reported) and at 4 + 1 + 1
    (16 records whose carriers hold a target that runs past its window,
    `limit-window-reach`, 2,121 carrier listings), three repeat-rich guides (VEGFA site 3,
    EMX1, VEGFA site 1) in both modes, HPRC chr22 (9,125 MNVs, 4.5% of its calls with a
    missing allele) and a copy of it with every fourth sample unphased in both modes, and
    the default index beside the sites-only `hg38_mega` in both modes (2,254,806 and 666 row
    blocks of `hg38_mega` list no carrier and are checked as rows, per sample 615,803 lines
    of the alternative alignments at one (AT)n repeat of chr22; completeness counts the
    127,735 indel records only `hg38_mega` holds and the 32,178 whose FILTER the search
    skips), each read as a run of 2.10.0 from its records, as is the README per-sample
    search with this release's code on the indexes of CRISPRitz 2.9.1 (format 1, the last
    record of a split site, a store of every FILTER), which passes too (4
    `limit-last-alt-wins` lines). That search takes 6.5 min there and 97 min on the indexes
    this release built, its SNP post-analysis 1 min against 76, in the windows every allele
    coded makes dense (7,493 on chr22, 7,255 hyperdense). The SNP+indel oracle with the
    second deletions of its sites O and Q, the insertion of Y1 and an MNV marked `LowQual`
    passes in both modes, and fails where the run read a store of every FILTER (14 and 6
    gating findings); the same world searched with `CRISPRME_INDEL_SNP=0` is refused
    (2.8.2's checker gave 166 and 77 gating findings there). On the documented chr22
    self-test (`complete-test`, as CI runs it) every benchmark passes soundness and
    `cas9_sg1617` completeness: CI's step, replayed verbatim, passes (about 4 min). The
    SNP+indel oracle passes in both profiles and both modes. The released 2.8.1's README
    searches give exactly the 4 rows per mode whose first aligned base is one an insertion
    adds (12 row blocks) and nothing else, and completeness their 1,306 carrier targets per
    sample (population-level the one such row placed more than `--merge` bp past its
    cluster), as they did to the checker of 2.8.2's code, with which the rest was measured
    (its post-analysis over the released 2.8.1 image, indexes of format 1): `cas12a_hbg`
    completeness on the self-test, removing `cas9_sg1617`'s 5 indel rows caught at each of
    their 5 clusters, CI's step failing on 2.8.1's self-test, the README guide and
    ACTGAAATCTGTAAGCAGGC without bulges passing. Corrupting 8 rows of the README per-sample
    run per class — the start moved by one, a REF base, the bulge count, a carrier added, a
    carrier dropped, the strand flipped, a target base — is caught on all 8, in the
    integrated table and in the alternative alignments. On 2.8.0, 10,093 row blocks per
    sample (4,661 targets that are not the reference with their named variants, the row at
    chr22:22322768 among them, 4,292 REF texts, 1,071 PAMs) and 1,931 population-level, and
    1,425 missing carrier targets. Genome-wide, soundness of the 2.8.1 README per-sample,
    population-level and TRAC runs reports that inside-insertion class only (1,111, 682 and
    1,332 row blocks, 70 to 100 s each with `--jobs 24`); completeness of the README
    per-sample run reports it again (63,651 carrier targets at 348 records, and 4 whose
    target at the anchor lost its highest-CFD block to such a row) and nothing else that
    gates, in about an hour on 40 processes (a background job, not the release gate). The
    SNP-row report finds nothing the documented rules do not explain on this release's
    per-sample runs but beside `hg38_mega` (below); on 2.8.1 and d7a6fff it reports 30 row
    blocks that list carriers whose slot of the SNP holds an indel that changes the target,
    which this release drops.
  - **Found with it**: three CRISPRme defects, each fixed in this release. On the SNP+indel
    oracle, the merge listed the unphased S19 on two alternative alignments weaker than its
    strongest: it keyed an indel row by Cluster Position, variant list and aligned target,
    not the aligned spacer, so the row of the alignment on which S19 was resolved with C
    (its SNP's base in a DNA bulge) was merged into an alignment whose column there faces
    the guide's A, which S19 carries (the merge keys an indel or MNV row by its aligned
    spacer too, below); a synthetic world of partly phased carriers searched through the
    real CLI (a deletion phased, the SNP it needs unphased) failed on the earlier merge key
    and passes on this release's, and `test/indel_snp_oracle/known_check_indel_rows.tsv`
    lists nothing. A search whose targets reach 27 nt (SaCas9: a 21-nt guide and a 6-nt PAM)
    lost, on an INDELS index of format 1 (CRISPRitz 2.9.1), those that start at a fake-indel
    window's first base (on chr22 of HGDP + 1000 Genomes, `CACTTTGGGAGGCCGAGGCGG` with
    `NNGRRT`, 2 mismatches, per sample: 2,103 carrier targets at 6 records, 2,088 of them
    the insertion chr22:17325232 A>AT's, the same on two fresh searches; an SpCas9 target,
    24 nt at most in these searches, never starts there holding its record): an index of
    format 2 holds them, and the same search of indexes CRISPRitz 2.10.0 built reports them
    all (its row at chr22:17325206 (+) lists the insertion's 2,088 carriers). CRISPRme's
    phased indel store kept every indel of a VCF whatever its FILTER, so on a VCF with
    records outside `PASS` and `.` a carrier was taken off an indel row's target, and off a
    SNP row, for an indel no fake-indel window of the search holds (the oracle with its
    sites O, Q and Y1 marked `LowQual`; Fixed, below). Not fixed, reported: `CFD_score` is,
    as documented, the higher of the REF and ALT scores, while the merge sorts on the ALT's
    (on an indel row the REF text is the reference read from its start, not a target: its
    score is the one printed on 99, 41 and 60 indel row blocks of the 2.8.1 genome-wide
    README per-sample, population-level and TRAC runs); the published `hg38_mega` index's
    SNP rows read its VCF otherwise than CRISPRitz's enricher does: they name SNVs whose
    FILTER is not `PASS` (chr22:11574550 C>G, `LOW_VQSLOD`; its registry holds all 78,116
    such SNVs of chr22) and an earlier record at 58 split sites of chr22 against the last at
    131 (57 of them a `PASS` record before a later `PASS` one, as chr22:11337076 G>C before
    G>T), so its SNVs outside `PASS` are rows CRISPRme may write or not, while its indel
    rows follow the enricher (the last record at 279 split sites, and at chr22:25093750 the
    `PASS` C>G before the `LOW_VQSLOD` C>T), as does its indel log (exactly its 222,283
    `PASS` indels of chr22); the variant genome CRISPRitz 2.10.0 builds from that VCF holds
    none of those SNVs, but its registry lists every ALT of a position and a sites-only
    dataset's SNP row substitutes each, so the fresh per-sample search beside it names one
    of the 74,025 SNVs of chr22 outside `PASS` that no passing record of either VCF holds in
    1 integrated row and 3,434 alternative alignments (chr22:12810728 (−),
    chr22_12810729_T_C, `ExcHet`, beside the `PASS` T>A the genome codes), the
    population-level one in none; and beside `hg38_mega` the merge, whose key for a SNP row
    is its Cluster Position and variant list, keeps `hg38_mega`'s worst-case alignment of a
    site under the HGDP + 1000 Genomes carriers of the same variants: 144 SNP row blocks of
    the per-sample run (65 lines of the alternative alignments, at 51 sites of chr22; 60 at
    23 sites on the indexes of CRISPRitz 2.9.1) list carriers for a target that holds more
    variant bases than the SNP they carry (chr22:11333697 (−), chr22_11333716_T_G: HG00271,
    HG00318, HG00369, HG00379, HG00864, NA06991, NA06993, NA18547), where the HGDP + 1000
    Genomes run alone lists them on an alignment that holds it (report-only).
  - Wired in: `PostProcess/test_check_indel_rows.py` (226 unit and mutation tests: a right
    row raises nothing gating, each corruption of it is caught, each documented limit and
    overlapping slot comes out under its own class, the checker's coverage gaps under
    theirs; the refusals; each release's reading of the data, told from a run's records;
    the same verdict under every hash seed; and the checker holds no CRISPRme code) runs in
    `unit-tests.yml`, and `test_ci_workflows.py` pins the CI step;
    `validate-benchmarks.yml` runs the checker on the chr22 self-test output after
    `validate-test` (soundness of the six benchmarks, completeness of `cas9_sg1617` on
    chr22, about 4 min; a gating finding fails the step once every benchmark is checked, so
    the uploaded findings hold them all); `test/indel_snp_oracle/run_oracle.sh` runs it
    beside `oracle.py check`, both profiles, both modes, soundness and completeness, and
    fails on any gating finding; `docs/RELEASE_CHECKLIST.md` step 4c makes a real chr22
    HGDP + 1000 Genomes run (both modes, soundness and completeness) and genome-wide
    soundness of the README guide and TRAC a release gate. `test/benchmark/README.md` says
    what `validate-test` covers and what it does not, and METHODS section 8 states the
    method.

### Changed
- **CRISPRme and CRISPRitz share their major.minor version; this release is 2.10.0, after 2.8.1.**
  CRISPRme X.Y.z requires CRISPRitz X.Y.* (2.10.z, then 2.11, ...), patch numbers are
  independent, and a breaking change in either package -- an index format, the enricher's output
  -- bumps the minor of both, so the versions that go together can be told at a glance (Luca,
  2026-11-07; `docs/RELEASING.md`, Version policy). CRISPRme skips 2.9.x to meet CRISPRitz 2.10.0;
  up to 2.8.1 the two numbered apart (CRISPRme 2.8.1 runs on CRISPRitz 2.9.1).
- **The image and `install_from_source.sh` build CRISPRitz 2.10.0, and the Bioconda recipe pins
  `crispritz >=2.10.0,<2.11`.** CRISPRitz 2.10.0 writes index format 2, codes every allele of a
  split multi-allelic site (the four-allele symbol Q where they cover A, C, G and T) and reads
  every record's `AF` (its CHANGELOG). The pin is the Dockerfile's `ARG crispritz_ref`, which
  `install_from_source.sh`'s `CRISPRITZ_REF` default equals; the recipe's bound, the next minor,
  keeps out a CRISPRitz that may change the index format or the enricher again (the version
  policy above; 2.8.1's was `>=2.9.1,<2.9.2`).
- **The two known limits of CRISPRitz `add-variants` that 2.8.1 stated are given as those of an
  index CRISPRitz 2.9.1 or earlier built, every one published until the rebuild for this
  release.** At a split multi-allelic site only the last record's ALT is searched, and a record
  whose `INFO` has `AF` at another place than the VCF's first record has none. METHODS (sections
  2, 3 and 8), the README, `docs/INDEX_PROVENANCE.md`, `docs/INPUT_FORMATS.md`,
  `docs/PRECOMPUTED_INDEXES.md`, `docs/DESIGN_mega_index.md` and the panel recipes of
  `seq_script/merge_panels/` gave them as this release's limits, most of them as to be fixed in
  the next release, and the input validator's warning of a VCF whose records hold `AF` at
  different places in `INFO` said that CRISPRitz reads it at the first record's place, where a
  later record gets a blank `AF`. Each now says that CRISPRitz 2.10.0 codes every ALT of such a
  site and finds each record's `AF` by its key, so that an index it builds, as every one rebuilt
  for this release, has neither limit. Covered by `PostProcess/test_site_alleles.py`,
  `PostProcess/test_generate_report.py` and `PostProcess/test_validate_inputs.py`.
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
  searched.** In the search's `log_verbose.txt` (the search script's stdout: not on the console of
  `complete-search`), once for each index it takes (a variant search with bulges takes three: the
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
  VCFs, no enriched genome, no per-sample dictionaries), gives the 10,927 rows of the install as
  built (above), byte for byte, and the three notes; `build-index-only` of the dataset stops as its
  notes say (exit 1, every `.bin` unchanged); with the reference index moved out of the way as its
  note says, `build-index-only` builds it in format 2, and the per-sample search gives the 10,928
  rows of format 2 above, byte for byte (the row more, the perfect match at chr22:11,118,987, is a
  reference site), with the other two notes, and the reference search the 9,211, with none. The
  same install as built (its VCFs and enriched genomes in place): the notes of its variant and
  INDELS indexes say to build them again from `Genomes/hg38+hg38_hgdp300` and
  `Genomes/hg38+hg38_hgdp300_INDELS`; `build-index-only` alone keeps the three indexes, every
  `.bin` unchanged; with each moved out of the way as its note says, it builds all three in
  format 2, and the per-sample and reference searches give the same 10,928 and 9,211 rows, with no
  note. With its VCFs and no enriched genome, those two notes say that
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
  INDELS index), and `, or run this search without --index-path` (under `--prune-pam` too, which
  no longer needs `--index-path`: an entry below) where the search then builds every index it
  needs: a reference search, or a variant one with the dataset's VCFs, enriched genome and
  fake-indel genome in place; without them it finds no variant index, and stops (`Error: no
  variant index here serves this search ...`). With CRISPRitz 2.10.0 (`feat/2.10.0` 7219b55), of
  indexes CRISPRitz 2.9.1 built, in a library
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
  note says, those two remain. Under `--prune-pam` (CRISPRitz 2.10.0, `feat/2.10.0` 881f761's
  build; a 3 kb genome holding one site of the README's guide, its `NGG_2` index built by
  CRISPRitz 2.9.1 in such a library; NGG, 2 mismatches + 1 DNA + 1 RNA bulge), the note ends as it
  does without it, where it said `, or run this search without --index-path and --prune-pam
  (which needs it)`, and the search it names, without `--index-path` and with `--prune-pam`,
  builds the index in format 2 and gives the same row, with no note. Each of these searches that
  does not stop ends with exit 0 and nothing on stderr. `test_index_format2` does the same through
  `build-index-only` and the search script on stub engines. The note `--prune-pam` gives of a
  flank wider than its PAM, and the warning of a flank built before CRISPRitz kept an `N` in it,
  are given of an index of format 1 only.
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
  side of the guide". It now records a PAM file on such an index only when it is of the side the
  header records: given one of the other side, it builds that file's index beside it, under the
  name that says its side (the names of this release, under Changed), saying the header records
  the other side, and it refuses one of its PAM alone, which says no guide length and so takes the
  index's own name; the web takes the index's geometry from its header, the PAM name of
  `.pam_build` only where that record's geometry is the header's, and offers the index to a
  search of that side alone; and a search that refuses an index for its side says which of the
  two it read. With CRISPRitz 2.10.0 (`feat/2.10.0` 881f761's build) and this release's code
  (cb550f1) on a genome of three contigs: of a 5′ PAM-less index that `crispritz.py index-genome`
  built (`NNN_2_bx`), the web reads `(None, 23, -3)` and offers it to a 5′ search, not a 3′ one;
  `build-index-only` with the shipped 3′ PAM file says `NNN_2_bx is an index for a PAM on the 5'
  side of the guide (as the header of its .bin files records): it is left as it is, and this index
  is built beside it, as 20-NNN_2_bx.`, records nothing on `NNN_2_bx` and the 3′ file on
  `20-NNN_2_bx` (exit 0, nothing on stderr), and a 3′ search takes `20-NNN_2_bx`; with a 5′ PAM
  file, `build-index-only` keeps `NNN_2_bx` and records that file, the web reads its name and the
  5′ geometry, and a 5′ search takes `NNN_2_bx`, each with exit 0 and nothing on stderr. Of a 3′
  index of flank 4, `build-index-only` with a 3′ PAM file of that flank keeps it (exit 0) and
  records that file.
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
  say which limits at contig ends and `N` runs are format 1's, and which were CRISPRitz 2.9.1's
  alone: a − strand PAM at a contig's first base reported on the + strand, which 2.10.0 reports on
  its own strand from an index of either format, and a + strand PAM on a contig's last bases found
  by no search, which 2.10.0 finds without bulges and from an index of format 2, but not from one
  of format 1 (`docs/PRECOMPUTED_INDEXES.md` said both of any search). They say what format 2
  holds, and that the equivalences of a PAM-less index with `--prune-pam` and of budgets hold
  between indexes built for the same guide length (an index built for 20-nt guides holds sites at
  contig ends and `N` runs that one built for 21-nt guides does not).

- **Two VCF records of one REF>ALT at a position are one allele: their carriers are merged, a
  haplotype slot carrying the ALT if either record puts it there, with the rsID and AF of the more
  frequent record (decision B2).** CRISPRitz 2.10.0 merges them in the dictionary, slot by slot
  where a sample's two calls are phased alike, and otherwise keeping the call with more copies of
  the ALT (the earlier on a tie) as it is; CRISPRme's dictionary readers (the Tier-0 and Tier-1
  compiles, the dictionary path of the SNP post-analysis, the indel path) merge two entries of one
  ALT, or a sample listed twice in one entry, that reach them from another enrichment, slot by slot
  where the calls are phased alike, and otherwise keeping the call with more copies of the ALT
  (the earlier on a tie) as it is, as the enricher does (`site_alleles.merge_gt`: CRISPRitz's
  `_gt_union`'s output on 20 pairs of calls, in its test; up to the integration it read a pair not
  phased alike as unphased, `1|0` beside `0/1` giving `1/0` where the enricher's entry holds
  `1|0`). An indel written on two VCF lines is one allele too, in one case or in two: the readers
  of the phased indel genotype store (the indel path's, the SNP rows' and the carrier gate's) read
  its REF and ALT in upper case and merge its two lines by the same rule (see Fixed). Up to
  2.8.1 two entries of one ALT stopped the Tier-1 compile of the whole chromosome (`duplicate
  (pos, alt) record`), and up to CRISPRitz 2.9.1 the later record's
  entry replaced the earlier's. A sample that the merge leaves with two ALTs of the position on
  one slot is resolved without phase (see Fixed). Covered by `PostProcess/test_site_alleles.py`,
  `test_indel_genotypes.py`, `test_snp_row_indel_carriers.py` and `test_haplotype_prioritizer.py`.
- **The input validator says what CRISPRme does with a VCF line that lists several ALTs, and with
  positions that hold several records, and looks for such lines on every run.** It warned only at 10
  ALTs or more, in the opt-in `--full_input_validate`, of an allele-index substring match the
  enricher stopped making in CRISPRitz 2.8.4 (CRISPRitz #44); it now warns on any line of two or
  more ALTs, in the full scan and in the first 1,000 records of each VCF on every run. A position's
  split records (`bcftools norm -m -any`, one ALT each, the published panels' form) no longer warn
  as dropped duplicates; two records of one REF>ALT warn that their carriers are merged; a line of
  several ALTs warns that each of them is given the line's first ID, and that split, each ALT's
  record holds the whole ID column, which a SNP's rows give whole and an indel's by its first ID
  (see the known limit below); and the full scan counts the positions that hold two or more
  lines of several ALTs (one a base each, FILTER `PASS` or `.`), whose allele numbers the readers
  cannot tell apart (see the B3 entry under Fixed), naming up to five, as the bounded scan of
  every run now does among a VCF's first 1,000 records. The survey of records that repeat a
  REF>ALT reads an indel's or MNV's REF and ALT in upper case, as CRISPRme reads them (HPRC chr22:
  179 indels on two lines, of which 140 in two cases, which it compared as written and counted
  none of), and an SNV's as written, as CRISPRitz compares two listings of one; and the full scan
  warns of an SNV whose ALT is written in lower case, which CRISPRitz does not code (so it is never
  searched, though its dictionary lists the carriers), with an `awk` line that writes REF and ALT
  in upper case (none in the published panels' chr22).
  `gnomAD-converter --multiallelic` says, in its help, on stdout and in
  the README, what it writes: a position's SNP records joined into one line and its indel records
  into another (`bcftools norm -m+`), with the converter's genotypes, which are unphased (`1/2`),
  and how CRISPRme reads them. Covered by `PostProcess/test_validate_inputs.py`.
- **A known limit, stated with its scale: an indel row gives its record's first ID, and so does each
  ALT of a VCF line that lists several ALTs.** CRISPRitz `add-variants` (2.10.0, as 2.9.1) writes
  the first ID alone into the indel log for every indel record, whatever the number of its ALTs, and
  into the per-sample dictionary for each ALT of a line of several; no reader of either can recover
  the others. A SNP of a record of one ALT keeps the record's whole ID column, which CRISPRme gives
  as one item (`rs1;rs2`). So where the ID column lists one ID per ALT (`rs1;rs2`, as
  `bcftools norm -m+` and `gnomAD-converter --multiallelic` join them), the rows of the second and
  later ALTs give the first's (`rs1`) as their `Variant_rsID`, with the tiers, with the dictionary
  or both; and a record of a split VCF (`bcftools norm -m -any`, the published panels' form) keeps
  the line's whole ID column, which its SNP rows give whole and its indel rows by the first ID. The
  mega's VCF lists two IDs for 68,017 of the 222,283 PASS indel records of chr22
  (`22:22652455:TTATG:T;rs1372713682`, whose rows give `22:22652455:TTATG:T`), the HGDP + 1000
  Genomes panel's for 8 of its 240,557; the sites-only (mega) registry builder, which names the
  mega's SNPs, gives each ALT of a line its own ID where the ID column lists one per ALT. Lines of
  several ALTs come from a user's VCF, HF's raw HPRC `vcfwave` VCFs (33,993 on chr22, whose IDs are
  graph IDs, not rsIDs), the 2019 HGDP data `complete-test` uses, and
  `gnomAD-converter --multiallelic`. METHODS section 8 and the input validator's warning on such a
  line say so. Checked through the real CLI with CRISPRitz 2.10.0: a line `rsB1;rsB2` A>G,C gives
  the rows of C `rsB1`; and the review's synthetic set, in all four layouts (tiers, tiers and
  dictionary, dictionary, population-level), gives a record `rsI1;rsI2` AA>A of one ALT `rsI1` on
  its indel rows and a record `rsS1;rsS2` A>C `rsS1;rsS2` on its SNP rows.

- **The `Observed` column of a population-level row of two or more variants says what the
  carrier gate found, in place of `N carriers (>=1 of K vars, poss. chimera)`.** The row's
  carriers are those of any of its variants, so the cell could not say whether anyone carries
  them together, and "chimera" named a guess. Where the gate checked the row, the cell gives the
  tier of its highest-CFD alignment: `together on N of H <panel> haplotypes` or `carried
  together by N <panel> samples, phase not known` (A); `not together on any of H <panel>
  haplotypes (expected E)`, the same with `; too rare to tell`, or `all N <panel> haplotypes
  carrying them together have an indel in the target` (C1); `not checked: K of N variants not
  in <panel>`, with the variants' frequency product and rarest frequency against the bound where
  the dataset gives them (C2); `the same CFD on N of H <panel> haplotypes, which carry <variants>`
  where a haplotype reaches the row's CFD without every variant (A). A sites-only panel's row,
  `putative` up to 2.8.1, reads so with
  `--prioritize-observed`. A row of one variant reads its carriers as before, unless every
  haplotype carrying it has an indel that takes the target's bases off it, when it reads `all N
  <panel> haplotypes carrying it have an indel in the target` (C1, ranked after the others as
  well: 172 rows of the genome-wide run, 13 of mega chr22). A per-sample row, and the rows of a
  search the gate did not check (or of a 2.8.1 search, read again), read as before. The report's
  legend, the inputs card, README and METHODS list the labels. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **`publish-index` uploads a variant index's phased indel genotypes on their own too**
  (`indexes/indel_genotypes_<vcf>.tar.gz`, beside `genotypes_<vcf>.tar.gz`; they stay in the
  index's archive as well), for `download --what genotype-oracle`, which installs them without
  the index. Covered by `PostProcess/test_haplotype_prioritizer.py`.
- **Under `--per-sample` a dense window's `High_variant_density_region` note says that the
  haplotypes its carriers have are reported, and so do the report's legend and its note on the
  dense windows near the top-ranked sites.** In either mode the note said, after the window's
  variant count,
  `a greedy worst-case alignment is reported here; additional haplotype alignments may exist`, the
  legend's `High_complexity_region` entry "Such dense windows are reported greedily", and the
  `--report-full` note "a single greedy worst-case alignment is reported for each". A per-sample
  search has reported a genotyped dataset's dense window by the haplotypes its carriers have on the
  dict-less tiers since 2.7.10, and does on the classic dictionaries since this release (see
  Fixed): its rows now say `the haplotypes its carriers have are reported here`, and so does the
  report. A row that names no carrier, a sites-only dataset's, whose rows are the population-level
  analysis's in either mode, keeps the greedy note, and the report says so of it where such a row
  is flagged; a population-level search reads as before. On chr22 of the HGDP + 1000 Genomes panel
  (five guides, NRG, 4 + 1 + 1, per sample) 50 rows carry the note, as many as with 2.8.1. Covered
  by `PostProcess/test_generate_report.py`.

- **`CRISPRME_IUPAC_CAP` counts bits, log2 of the sequences a window spells, where it counted the
  window's ambiguous columns.** Above it a `--per-sample` search of an index without the dict-less
  tiers reports a window by the haplotypes its carriers have, read from the dictionaries as the
  tiers' path reads them (the carrier gate's entry under Fixed), instead of every combination of its
  alleles (a population-level search reports every window by its greedy worst-case alignment, as
  before), and every path lists the window in `high_variant_density_regions.bed`. The default, 10,
  caps a window of two-allele columns where the count of columns did, so that on an index whose
  windows hold no other variant column the unit changes nothing: on the published HGDP + 1000
  Genomes index, whose chr22 holds no three-allele column, the README guide's searches of chr22 list
  the same 565 windows over the cap as with 2.8.1 (8,155 lines, per sample and population-level),
  the notes giving the windows' sequences (and per sample no greedy worst case: Fixed, below); their
  rows differ from 2.8.1's otherwise only as entries under Fixed say: the order of their rsID and AF
  lists; per sample, the carriers whose haplotype holds an indel in the target, whom a SNP row no
  longer names (3 of the 1,574 rows gone, 4 changed); in both modes the 4 rows that started on a
  base an insertion adds, at its anchor now, one of them (chr22:21052642 (−)) an alternative
  alignment of the row at 21052638 (−) (1,570 rows per sample, 2,197 population-level, against
  2.8.1's 1,574 and 2,198), and one row's `Variant_info_spacer+PAM` and one's `PAM_creation`
  (population-level also two CRISPR-Bulge scores in their last digit, scored in another batch). A
  three-allele column weighs 1.58 bits and a four-allele one 2, where each counted 1: coding every
  allele of a split site makes them common, and the immunoglobulin loci explode (on chr22 of HGDP +
  1000 Genomes, five guides' variant search returns 4,372,248 raw targets with CRISPRitz 2.10.0,
  96.6% of them in IGL, where 2.9.1 returned 206,883), which the count of columns could not see. A
  VCF line that lists several ALTs writes a three-allele column with any release: on a small genome
  whose phased VCF has seven such lines in one window of the README guide (3^7 = 2^11.09 sequences),
  2.8.1 listed the window nowhere; it is listed now, and under `--per-sample` an index built without
  `--samplesID` reports it by the haplotypes its carriers have, as one with the tiers does (the
  carrier gate's entry under Fixed, whose test of the cap counts bits as this one does). `-1` still
  disables the cap; a hyperdense window is listed all the same.

- **SpCas9 NRR (`NAA`, `NAG`, `NGA`, `NGG`) is the default PAM and index; it was NRG (`NAG`,
  `NGG`).** The broadest SpCas9 screen of the three published ones, and the one recommended for an
  IND package: an NRR search also reports the `NGA` and `NAA` sites, and the variant-created
  ones, that an NRG search leaves out, as NRG reports the `NAG` ones NGG leaves out. The PAM file
  `PAMs/20bp-NRR-SpCas9.txt` ships with CRISPRme, and a download of the PAM files (`download --what
  pams`, and so `--what all`) installs it where the data repo has none (the repo's `pams/` holds the
  NRG, NGG and PAM-less files). The web form selects it when it is installed (then NRG, then NGG, as
  before), `download --what all` fetches the reference index `NRR_3_hg38`, and every surface that
  names the recommended index names `NRR_3_hg38+hg38_HGDP1kGP` (and its HPRC and mega panels):
  the README's and the website's commands, the one-click installers (`docs/install.sh`,
  `docs/install.ps1`), the Docker Quickstart, the web user guide, the CLI guide,
  PRECOMPUTED_INDEXES, METHODS, the Hugging Face cards, RELEASE_CHECKLIST's default-index gate, the
  help of `download`, `setup` and `complete-test`, complete-test's advice when no index is
  installed, the advice of a search and of the report to download a tiered index, and the notes
  of the panel build recipes (`seq_script/merge_panels/`) and of the paper's dataset table.
  complete-test's smoke searches a PAM-less 8-nt-flank index for NRR (it searched it for NRG). The
  NRG and NGG indexes stay published and searchable. What it costs: the default install is a ~70 GB
  download (it was ~45 GB) and takes ~150 GB on disk (it took ~98 GB; measured 150,048,109,201
  bytes: indexes 94.2 GB, `Dictionaries/` 52.4 GB, genome 3.3 GB), with a peak of 161 GB while the
  archives unpack, so every surface that sizes it (the installers, the README, the website, the
  Quickstart, TROUBLESHOOTING, DATA_LAYOUT, install/README) asks for ~200 GB free (it asked for
  ~150 GB); a reference-only setup takes ~44 GB (it took ~25 GB). A genome-wide search takes
  about twice as long (2 h 26 min against 1 h 14 min for `ACTGAAATCTGTAAGCAGGC` on HGDP + 1000
  Genomes, per-sample, 4 mismatches + 1 + 1 bulge, 32 threads, CRISPRme 2.8.1), at no higher peak
  memory (66 GB, against 82 GB summed over the NRG run's processes). Covered by
  `PostProcess/test_defaults_and_index_flex.py`.
- **`--prune-pam` no longer needs `--index-path`**, and a search passes it by itself (see Added).
  The refusals that told a user to run again with `--prune-pam --index-path` now say, where the
  engine cannot prune, that an index serves the search pruned and the installed CRISPRitz cannot do
  it; the advice to run without `--index-path` no longer says to drop `--prune-pam` too. The
  warning of an all-`N` index built before CRISPRitz 2.9.1, which holds no target whose flank has
  an `N`, said to "search without `--prune-pam` and `--index-path`, which builds a dedicated
  index": such a search now prunes by itself from that very index and builds none, so the advice
  led back to the same warning. It names the two routes that serve the search whole: the index
  built again, or a dedicated index of the PAM beside it, which a search takes first. And the NOTE
  of a search that prunes by itself offers a published index "if it is published", as the refusals
  do — none is published for a genome of one's own, or for `TTTV`, `NGA` and most other PAMs.
  On a 30 kb genome whose `20-NNNN_2_gn` was built before an N run was put into it, the search
  the older advice named — no `--prune-pam`, no `--index-path` — takes that same index, pruned
  by itself, and builds none (`genome_library` holds `20-NNNN_2_gn` alone after it), in 27 s
  with nothing on stderr.
- **An index is named for its guide length and the side of its PAM (E2).** `<PAM>_<N>_<genome>`
  says neither, so a 3′ and a 5′ build of one motif, or builds for two guide lengths (Cas12a
  `TTTV` at 20 and 23 nt, as the rebuild adds), took one name, and the second build stopped on
  the first ("... is already there ... The index for that PAM file has the same name: build it in
  another working directory (--path)"), as a search did. `build-index-only` and a search now
  build an index under a name that carries both, in the order of the target on the strand:
  `20-NRR_3_hg38` for a 20-nt guide and a 3′ PAM, `TTTV-23_3_hg38` for a 5′ PAM and a 23-nt
  guide, the variant index (`20-NRR_3_hg38+hg38_HGDP1kGP`) and its INDELS index
  (`<name>_INDELS`) too. The engine writes the older name, in a build folder (a search's in the
  index's build lock), and the index is renamed into place: an index of the older name in place
  is never touched (in place, the engine removed a folder of the name it writes before writing
  it). Every reader takes both forms (the search script, `crisprme.py`, the web form, `download`
  and its router, which reads a PAM-less name's flank in either form, `publish-index`'s
  manifest, complete-test, the index coverage, FASTA-header and N-flank checks): an index of the
  new form only for a search of the side its name says, which records the side as a `.pam_build`
  does, ahead of one (such an index needs none, at any flank width), and of one index under both
  names, the newer one. An index of the older name (every one published so far, or built by an
  earlier release) is used as before. `build-index-only` keeps one in place that is the index it
  builds (its `.bin` files built for the PAM file's guide length, its `.pam_build` recording that
  side, or none: a dedicated index, an all-`N` one for a 5′ file, or one at the published widths,
  3 and 8), and builds beside one that records another, saying so (`.../NNNN_2_hg38 is an index
  for a PAM on the 3' side of the guide (as its .pam_build records): it is left as it is, and
  this index is built beside it, as NNNN-20_2_hg38.`). An unrecorded all-`N` one of another
  width, for a 3′ file, is refused as before, the ERROR now naming the `mv` to each side's name,
  which records it; so is a folder of the new name that its trees or its record contradict, with
  the name of what it is. A search's refusal of an unrecorded index for a 5′ PAM names that `mv`
  too (it said to run `build-index-only` with that file, which still records it). An index is
  published under its folder's name, so one built by this release under the new one. An earlier
  release does not read the new form, and builds its own. The documented builds and publishes name
  the folder a build writes: PRECOMPUTED_INDEXES' dict-less flow, which the rebuild follows, said
  that its build wrote `NRR_3_hg38+hg38_HGDP1kGP`, and its `publish-index` of that folder stopped
  on `Index directory ... does not exist`, and the release checklist's build said
  `NNN_3_hg38+hg38_HGDP1kGP`; the panel recipes (`seq_script/merge_panels/`) print the index their
  build wrote, of either form, as the CRISPRme they run names it, and the mega INDELS recipe's
  stop on an unexpected name says which to set (`20-NRG_3` from 2.10.0 on). On chr22 with HGDP +
  1000 Genomes, the documented build wrote `20-NRR_3_hg38`, `20-NRR_3_hg38+hg38_HGDP1kGP` and its
  `_INDELS`, with nothing on stderr; the documented `publish-index` of the variant one went past
  its folder to the token it asks for (none was given), where the older name stops on the missing
  folder; the HGDP + 1000 Genomes recipe's lines after its build printed that folder; and the
  page's search under `--index-path` took both indexes. Covered by
  `PostProcess/test_build_index_only_atomic.py` (both sides and two guide lengths built in one
  library, the variant and INDELS indexes too, an index of the older name kept or built beside,
  each refusal and its `mv`), `test_index_resolution.py` (a search builds its own under the new
  name beside an index of the older name for the other side or another guide length, which it
  stopped on), `test_defaults_and_index_flex.py` (every reader takes both forms alike; every
  documented build and publish names the folder of its build, and each recipe prints the one its
  build wrote) and, through
  the real engine, `check_index_flex.py`: the 3′ and 5′ `NNNN` indexes, and the 20- and 23-nt
  `NGG` ones with their variant and INDELS indexes, built in one library, each searched with the
  results of that index alone, byte for byte; an index under the older name gives the same
  results; a 5′ build beside a 3′ index of the older name, each side's search taking its own; one
  named for a 5′ PAM is not taken for a 3′ search; and a PAM-less index named for its side with no
  `.pam_build` serves a search of that side (a 3′ one of width 4 an `NGG` search, pruned; a 5′ one
  a 5′ PAM of its width).
- `PostProcess/pool_search_indels.py` drops `_index_pam`, a helper nothing called (emptyidx F8);
  the results page's sample tab no longer computes paths it never reads (wsearch N8); and
  `remove_contiguous_samples.py` reads a variant token as an indel's in one place, for the
  representative's tie-break and the carrier resolution alike (RV3-D).
- **Alignments that start on bases of one insertion share its anchor, and the merge takes
  them for alignments of one site** (follows from the first fix below; both modes).
  - Each such row now starts at the insertion's anchor, so alignments that start on different
    bases one insertion adds have one start and one cluster position. The merge takes a guide's
    rows on one chromosome for one site while the cluster position of each (its start, less its
    RNA bulges and plus its DNA bulges on the + strand with a 3′ PAM and on the − strand with a
    5′ one) is within `--merge` bases (3 by default) of the one before, on either strand and
    whatever their bases, keeps one of them as the site's row and lists the others as its
    alternative alignments: these alignments are now one site, and such a row can also join the
    site of a reference or SNP row that starts near the anchor. README guide, HGDP + 1000
    Genomes: chr22:21052642 (−) of 2.8.1, a CFD-0 alignment of T>TCACC at 21052640, 4 bases
    from the reference row at 21052638 (−), is now 21052639, one of that row's alternative
    alignments, and the reference row's `Other_motifs` goes 2 → 3 (per-sample, 2 again with the
    SNP-row fix below, which takes off that site's SNP-only alternative alignment; a
    population-level run, which that fix leaves alone, keeps 3).
  - Within a site, the alignments that name the same variants with the same carriers are one
    entry, the best: one target's alignments with different bulges, as on every row. The merge
    compares indel rows there at their cluster position on the haplotype carrying the indel,
    which on a row that starts inside an insertion is past the anchor by the rank of its first
    base among the inserted ones, read from the row's bases. Alignments whose cluster positions
    on that haplotype differ are other targets of it, as many bases apart as 2.8.1 placed them,
    so each is listed, as the site's row or one of its alternative alignments (compared at the
    anchor, all but the best of those that start inside one insertion would leave both tables);
    one target's alignments with different bulges before one PAM still fold into the best, also
    when they start on different bases of one insertion. Two copies of one target inside an
    insertion (a repeat it holds) make the same row at the anchor, its bases, REF alignment and
    variants the same, and that row is listed once, with the carriers of both, where 2.8.1
    listed each copy at its own start (its `PAM_creation` is created when either copy's PAM
    is, below); a bulged alignment of the second copy, which 2.8.1 folded into that copy's
    row, has no row of its own copy left to fold into and is listed as one of the site's
    alternative alignments (in the SNP+indel oracle's X2, the second copy's DNA-bulge
    alignment). TRAC guide (CTCTCAGCTGGTACACGGCA), HGDP + 1000 Genomes, per-sample:
    two alleles of chr22:50193929 that insert 96 bases, a 48-base unit twice, with one carrier
    each, hold the CFD-0 alignment CTACTCAGgaGGgACA-GGCtGAG (−) in both copies; 2.8.1 listed it
    at 50193950 and 50193998, 2.10.0 once per allele at 50193928, so the insertion's 2 rows and
    7 alternative alignments of 2.8.1 are 1 row and 6 alternative alignments. Genome-wide,
    2.8.1's tables merged as 2.10.0 merges them: of the alignments that start inside an
    insertion, the highest-CFD ones, 5 of 371 (229 rows and 142 alternative alignments) are
    such a second copy in the README search per-sample, 1 of 228 rows population-level and 22
    of 444 (251 and 193) with the TRAC guide per-sample, all at CFD 0.004 or less; per-sample,
    every other one is still listed, as a row or an alternative alignment.
  - Cost: the merge reads a row's place on the haplotype once per target, REF alignment, start
    and variants, and cuts the haplotype from the indel's ALT in slices, so a long insertion
    adds no step per inserted base. HPRC holds 402 insertions of 53.9 to 54.4 kb at
    chr22:10690626, and the README guide's chr22 search puts 16,058 rows inside them, now all at
    their anchor (2.8.1: at most 760 rows at one cluster position). On that index, the released
    2.8.1 image and 2.10.0's post-analysis on it side by side on one host, per-sample: each of the
    three merges takes 4.4 to 4.6 s against 3.1 to 3.5 s, and the merge step, which is mostly
    sorting, 77 to 80 s against 54 to 56 s (its sorts by cluster position now meet those rows at one
    place: the highest-CFD file, shuffled, sorts in 15.3 s against 9.8 s); population-level 4.6 and
    77 s against 3.4 and 53 s; TRAC guide, per-sample, 3.8 to 4.0 s and 53 s against 3.0 to 3.1 s
    and 50 s. On chr22 of HGDP + 1000 Genomes every merge of the searches below takes under a
    second, as in 2.8.1. Covered by `test_indel_snp_cis` (the sliced reading against the
    base-by-base one on 4,000 random rows; the Python lines it runs do not grow with the insertion)
    and `test_representative_tiebreak` (one target's 20 alignments inside a 54 kb insertion: one
    reading).
  - The same search on chr22, released 2.8.1 image against 2.10.0's post-analysis on it (CRISPRitz
    2.9.1, index format 1) without the SNP-row fix (the last entry below, which gives the counts
    with it): per-sample, 1,574 → 1,573 integrated rows and 3,756 → 3,757 alternative alignments
    (21052642 joins 21052638's site; that site's two alternative alignments go 1 → 2 and the moved
    row, now one of them, loses `Not_found_in_REF`, its site having a reference row);
    population-level, 2,198 → 2,197 rows (a population-level run writes no alternative alignments).
    NRR, per-sample: 3,987 → 3,985 and 11,581 → 11,583, 21052642 and 33175772 (+), now 33175769,
    joining the reference rows' sites at 21052638 (−) and 33175766 (+), whose `Other_motifs` go 2 →
    3 and 9 → 10, their alternatives' 1 → 2 and 8 → 9. TRAC guide, per-sample: 2,803 → 2,800 and
    7,508 → 7,509. 37564129 (+) and 36738101 (−), now 37564118 and 36738081, join the reference
    rows' sites at 37564117 (+) and 36738082 (−); 49196655 and 49196663 (−), two alleles of one
    insertion, are one site at its anchor, 49196609; an alternative alignment of 49196601 (−) that
    starts inside another insertion is now the row of a site of its own, at that one's anchor,
    49196564; and 50193998 (−) is a second copy (above).
  - The outputs built from the rows follow these changes and the fixes below: the summaries by
    guide and by sample, the report (its counts and top-1000 tables),
    `observed_haplotypes_raw`, the per-site `observed_haplotypes.tsv` (in the report's zip,
    `data/observed_haplotypes.tsv`) and the plots, and the companions that give a row's start
    (`offtarget_start`): the SNP+indel co-occurrence companion (`*.indel_snp_cooc.tsv`) and, on
    a sites-only (mega) index, the per-dataset indel frequency companion (`*.indel_af.tsv`).
    NRG per-sample, without the SNP-row fix: the 1,292 samples listed on 21052642 count one
    target and one created PAM fewer (the row joined 21052638's site, and its PAM is not
    created, below), and the other 2,799 samples one target fewer in their population or
    super-population; the co-occurrence companions and `observed_haplotypes.tsv` of the chr22
    searches above are unchanged (with the SNP-row fix below, `observed_haplotypes.tsv`, and
    with the TRAC guide the SNP+SNP co-occurrence companion, follow the carriers it takes off).
    In the SNP+indel oracle the 197 co-occurrence lines of Z3, which starts inside an insertion,
    go from 15702, 15703 and 15704 to 15699.

### Fixed
- **An index folder that holds a macOS AppleDouble file is read by its own trees.** macOS `tar` and
  `cp` leave an AppleDouble file `._<name>` beside each file they copy onto a volume without
  extended attributes. Every reader of an index took `._NGG_chr22_1.bin`, which opens with no
  magic, for one of its trees, of index format 1 and built for a 512-nt guide (its second int32).
  So CRISPRme 2.8.1's `build-index-only` kept in place an index built for too short a guide:
  over the NGG_3 index of hg38 chr22 built for 20-nt guides, with a 21-nt PAM file at 2 DNA + 1 RNA
  bulges, it said `Reference index already present` (exit 0) with such files beside the trees, and
  stopped without them (`... is already there, built for a 20-nt guide ...`); complete-search's
  scan and the web took such an index as fitting any guide, and the search script did too where
  the file system listed such a file before the trees. With index format 2 read, more went wrong:
  a search of an index of format 2 beside such files stopped, saying it was put together from two
  builds and to remove it (`ERROR: index 'NGG_3_hg22' (...) holds trees of index formats 1 and 2,
  which no build writes together: it was put together from two builds: ... Remove it (rm -rf
  '...') and run this search again ...`), under CRISPRitz 2.9.1 too, in place of the ERROR that
  names its format; `build-index-only` and `publish-index` refused it the same way;
  `build-index-only` recorded nothing of it, and of an index of format 1 it recorded the PAM `.` of
  a 512-nt guide; and the `--prune-pam` note of a wide flank of format 2 was that of format 1. Every
  reader now takes the index's trees alone, a `.bin` that is no dot file: the search script's
  probes, `index_fit`, the coverage check and the check of a build. The reference and variant
  searches hand the engine a view of those trees alone where such a file is there: CRISPRitz
  2.10.0's `searchTST` stops on a folder that holds one (`ERROR: the .bin files in '...' are named
  for more than one PAM ('.', 'NGG'): ... remove the others.`). `publish-index` leaves `._<name>`
  files out of the archive, and has the `tar` of macOS add none of its own. With CRISPRitz 2.10.0
  (`feat/2.10.0` 2d1f57b), on hg38 chr22 (`build-index-only` of NGG at 2 + 2 bulges, then
  `complete-search` of `ACTGAAATCTGTAAGCAGGC` at 4 mismatches + 1 DNA + 1 RNA bulge), with a `._`
  file beside each file of the index (its two trees, `.pam_build` and `.index_meta.json`): of an
  index of format 2, the 779 rows of the index without them, byte for byte, where the search
  stopped as above; the same of an index of those two trees alone, with no record beside them, which
  the engine, handed the folder as it is, refuses; of an index of format 1, those 779 rows and the
  note of format 1; and under CRISPRitz 2.9.1, the ERROR that names the format of the index of
  format 2, and the 779 rows of the one of format 1, as CRISPRme 2.8.1 gives them. With the 21-nt
  PAM file, `build-index-only` stops on either index as above (on 2.9.1, before: `Reference index
  already present`); over either index found in place without its record, it records what it
  holds (its format, NGG, the 3′ side, a 20-nt guide); the archive `publish-index` makes of the
  index with the image's `tar` and `pigz` holds no `._` file (before: four); and a `--prune-pam`
  NGG search off an 8-nt flank of format 2 so copied gives its row and no note. Each search and
  build that does not stop ends with exit 0 and nothing on stderr. The INDELS search hands each
  chromosome's sub-index to the engine as it is, as 2.8.1 does: CRISPRitz 2.10.0 stops on one that
  holds such a file, with the ERROR above.

- **A VCF line that lists several ALTs is read one ALT at a time, each carrier's genotype against
  that ALT's own allele index where its genotypes tell it (decision B3).** Such a line keeps its
  raw genotype (`1|2`) under every ALT its sample carries, and the readers took "1" for "this
  ALT": the Tier-0 registry
  counted the second ALT's copies as reference (on chr22 of HGDP + 1000 Genomes with its split
  records joined into such lines, 42,952 of its 43,389 positions of several SNP ALTs had an ALT
  with a wrong global AC: chr22:10,570,455 T, AC 2 in the registry against 2,977 of 8,182 in the
  VCF, and Variant_MAF was wrong on 3 to 10 rows per guide); the observed path guessed an ALT's
  index from its own carriers' most frequent token (`_alt_index_for_record`, wrong for 449 of
  45,837 ALTs, 2 or 3 carriers lost per search); and the dictionary path put a `1|2` sample's
  first ALT on both its slots. The allele index of an ALT is
  now the token all its carriers share, narrowed by what the position's other ALTs of a line of
  several ALTs need (every ALT of one character the dictionary lists there, the `*` of a spanning
  deletion among them, in every reader: the indel path's reader of an install without a Tier-1
  store dropped the `*` first, so a `1|2` carrier of `G>T,*` whose `*` another sample carries
  alone was read without phase there and with T on its first slot through the Tier-1 store). An
  ALT whose carriers hold no index but 0 and 1 settles none: CRISPRitz 2.10.0 lists every line of
  a position, and the `1` of a line of one ALT, or of the entry it merges from two lines (written
  as one of a line of one ALT), numbers that line, not the others. Read with it, `C>G` (S05 `1|0`)
  beside `C>A,CT` (S06 `1|2`: A on its first slot, the insertion on its second) tied A to `2`, the
  insertion's number, and every reader put S06's A on the insertion's slot: S06, which holds A
  beside the deletion `TC>T` on its first slot, a CFD-1.0 target 2.8.1 reported, was listed only
  for the deletion's 1-mismatch target (CFD 0.455; the same beside `C>A,*`, and on SNP rows a
  carrier listed for two targets its haplotypes hold neither of: the review's worlds X and Y
  through the real CLI, 28 gating findings of the real-data check); the same through a listing
  merged from two lines (`A>C` beside `A>T,C`). Within one line of several ALTs no such ALT
  settles an index another needs (a carrier of two one-base ALTs is listed under both, with both
  indexes), so the narrowing loses nothing there. Where the genotypes cannot tell (every carrier
  of two ALTs is `1|2`, or of an SNV and an indel or a `*` of its line, which no dictionary entry
  holds, beside a line of one ALT of the position too), the carrier holds one copy on an unknown
  slot and is resolved without phase (PUTATIVE): S06 is listed, PUTATIVE, for its CFD-1.0 target
  again, and worlds X and Y pass the check in both modes. Two lines of several ALTs at one
  position are not told apart (an ALT of one can be tied to an index of the other's: the
  pre-flight check of every run warns of such positions among each VCF's first 1,000 records,
  the full input validation of every one, and the real-data check refuses a run that holds
  one); no published panel holds one on chr22 (HGDP + 1000 Genomes, HPRC and mega: 0 lines of
  several ALTs, 0 `*`), and the 33,993 lines of several ALTs of
  HF's raw HPRC `vcfwave` chr22 VCF are each alone at their position. Writing each ALT of such a
  line in split form in the dictionary, which would tie every one, is CRISPRitz's (left to
  CRISPRme by B3); the population summary of a phased panel
  decides it without phase too (it read the unknown slot as slot 0, a cis copy of the row of each
  ALT of the line). The Tier-1 store holds each ALT's genotypes read so (`1|0` and `0|1` for a
  `1|2`), and a store of an earlier release is read the same way. The sites-only (mega) registry
  builder and its indel AF companion read such a line one ALT at a time, each with its own AF and,
  where the ID column lists one per ALT, its own ID; they skipped it whole (the dictionary gives
  each ALT the line's first ID: see the known limit under Changed). A position of split or biallelic
  records reads as before. These lines come from HF's raw HPRC `vcfwave` VCFs (33,993 on chr22), the
  2019 HGDP data `complete-test` uses, and `gnomAD-converter --multiallelic`. Covered by
  `PostProcess/test_site_alleles.py`, `PostProcess/test_build_mega_registry.py` and the SNP+indel
  oracle, which now fails a per-sample SNP row naming two ALTs of one position: the per-sample run
  of rel/2.8.2's code (fe0585a) lists its site N2's lone `2|1` carrier of an ALT on a line of
  several ALTs on 3 such row blocks, which the oracle counted and did not check, as it did for the
  counts it gave before this release; this release's lists none.
- **A haplotype slot holding two ALTs of one position is no haplotype: the sample is PUTATIVE
  there, a carrier of each ALT's row, and no row names two ALTs of one position (decision B4).**
  With every ALT of a split site in the dictionary, a panel's split records phased one by one put
  two ALTs on one slot (672 sample-positions on chr22 of HGDP + 1000 Genomes): the per-sample
  analysis made it one CONFIRMED haplotype holding both, a row naming two alleles at one position
  (chr22:10,661,629 G and T for NA18869, CONFIRMED in 190 alignments), and listed the sample under
  neither ALT's row. The PUTATIVE branch now takes at most one ALT of a position in each
  combination, the union above `CRISPRME_IUPAC_CAP` variants included: one union per choice of ALT,
  at most 16 a sample, and past that one, holding at each such position the ALT nearest the guide
  (uncapped, NA19042, PUTATIVE in the IGL window chr22:22,881,180-22,881,205 and carrying two or
  three ALTs at 16 of its 26 positions, gave 8.5 million unions there, and 90 GB; now 331 distinct
  haplotypes in all for the window's 790 carriers, in 0.1 s). Through the real CLI, on an index
  CRISPRitz 2.10.0 built from the 10,429 HGDP + 1000 Genomes records of chr22:22,860,001-22,899,000,
  a per-sample search of TGGGGGATAAATATGCTTGC (NGG, 2 mismatches; its protospacer is
  chr22:22,881,187-22,881,206) now takes 58 s in all, its largest process at 0.55 GB; before, its
  SNP post-analysis was stopped after 40 minutes, at 15 GB, having treated 76 of the 384 alignments
  that hold more than 10 variants. The population summary decides such a carrier without phase, as
  it now decides a carrier whose genotype at a variant of the row is unphased; the SNP+SNP
  co-occurrence companion counts positions, not variants, so two ALTs of one position are no
  co-occurrence; and a sites-only dataset's co-located union takes one ALT of a position (one union
  per choice, at most 16 a window). The phase label is per row, so a row such a sample shares reads
  PUTATIVE for all its carriers. An indel row resolves a carrier whose indel's slot holds two ALTs
  of a position PUTATIVE too, among those ALTs: its ALT copies there fill its slots, so the
  reference is no option (METHODS section 4). Covered by `PostProcess/test_observed_haplotypes.py`,
  `PostProcess/test_site_alleles.py` and `PostProcess/test_indel_haplotype_resolution.py`.
- **The dictionary path keeps only the ALTs a target's own column holds, as the observed path does
  (#139); a sites-only dataset substitutes every ALT its registry lists, as up to 2.8.1.** The
  dictionary path, and a genotyped panel's registry read alone (`download --no-genotypes`),
  substituted every ALT listed at a position, an ALT the code the engine matched does not hold
  included, so the two paths could report different rows for one window: those sources list what the
  enrichment coded, so such an ALT is another position's. A sites-only (mega) registry is read from
  its VCF, every record of it, and lists ALTs the enrichment never coded: on an index built by
  CRISPRitz 2.9.1 or earlier, as every mega index published until the rebuild, the earlier records
  of a split site (the enrichment codes the last PASS record's), and on any index non-PASS records.
  Each is the panel's own allele at its own position, so a sites-only dataset substitutes every one
  at a variant column, whatever its code, a Q included, as 2.8.1 did. Applied there, the dictionary
  path's filter lowered the highest CFD of 114 sites in a search of the published
  `NRG_3_hg38+hg38_mega` (chr22; the README's guide, ACTGA, TRAC, EMX1 and VEGFA site 3; NRG, 4
  mismatches + 1 DNA + 1 RNA bulge, population-level; against 2.8.1): 110 for a PASS ALT the 2.9.1
  code left out (chr22:10,752,820 C>A, whose later record C>T is the coded one) and 4 for a non-PASS
  one, 22 of them from 0.1 or more and 2 from 0.5 or more (VEGFA site 3's chr22:22,893,912 from
  0.585 to 0.327, TRAC's chr22:22,900,049 from 0.524 to 0.127); the alleles track's code, before the
  fixes of 2.8.2 were merged, gave 2.8.1's 20,338 sites with 2.8.1's highest CFD at every one, and
  no row naming two ALTs of one position (2.8.1: 4 of the integrated table's), with nothing on
  stderr (CRISPRitz 2.9.1, the 2.8.1 image). Those fixes place a row that starts on a base an
  insertion adds at its anchor and compare an indel row at its cluster position on its haplotype,
  where a population-level search writes one row a site (entries below): this release gives 20,010
  rows there (792 of the 833 rows of 2.8.1 it no longer writes name an insertion), none naming two
  ALTs of one position, with nothing on stderr (CRISPRitz 2.10.0; the integration's first stage,
  without those fixes, gives the 20,338). A target that only such an ALT brings within the search's
  budget is still not found on that index: the mega index rebuilt by CRISPRitz 2.10.0 codes every
  PASS ALT (METHODS section 2). Covered by `PostProcess/test_site_alleles.py`.
- **On a row of several variants, the rsID and AF lists follow its variants: one rsID and one AF a
  variant, in the order of `Variant_info_genome`, and two variants of one frequency give it each.**
  The merge of a cluster's rows (`remove_contiguous_samples.py`) sorted the variant, rsID and AF
  lists each on its own, as text, and kept each value once (the known limit 2.8.1 stated: on ACTGA,
  953 of 1,967 rows of two or more SNPs out of order). It now keeps (variant, rsID, AF) together; a
  variant two datasets of a search give with two rsIDs or AFs lists both pairs, one after the other,
  and an rsID field holding two IDs (a VCF ID `rs1;rs2`, which the dictionary holds as `rs1,rs2`)
  reads `rs1;rs2`, one item. The report's MAF, the lowest AF of the row, is unchanged. Covered by
  `PostProcess/test_site_alleles.py`, `PostProcess/test_generate_report.py` and
  `PostProcess/test_representative_tiebreak.py`.

- **A per-sample SNP row no longer lists a carrier whose haplotype holding its SNPs also holds
  an indel that changes the target** (every release with per-sample genotypes; the known limit
  2.8.1 stated).
  - The enriched genome holds SNPs only, so a SNP row's carriers came from SNP genotypes alone,
    and an individual was listed whose own sequence there is the indel's (on the indel row, when
    in budget). On chr22 of the HGDP + 1000 Genomes panel an alternative alignment at
    chr22:21052638 − listed the 1,292 carriers of 21052641 G>A, which sits on the same 1,448
    haplotypes as 21052640 T>TCACC: none of them has that target.
  - Such a carrier now leaves the row, read on the phased indel genotype store the indel
    post-analysis reads, on the dict-less tiers and on the classic dictionaries alike, by the
    rule 2.8.1 applies to an indel row's second indel. An indel
    changes the target when it deletes or replaces a base of its span (protospacer and PAM) or
    inserts bases between two of them; not when its deleted bases all lie outside the span (its
    junction at or past the target's edge), nor when it inserts before the first base or after
    the last. Whether the bases change is read on the haplotype's own bases: a repeat unit gained
    or lost at the target's edge leaves them, and the carrier stays. A phased carrier is decided
    by its haplotype (an indel on the other haplotype takes nothing off); one homozygous for an
    indel that changes the bases leaves the row; one whose indel call is unphased cannot be tied
    to a haplotype and stays, PUTATIVE. An indel at the position of a SNP the row names is
    another allele of that base, a phasing artefact, and takes no carrier off; one whose REF
    holds the base of a SNP the row names at another position puts two records of that base on
    the carrier's slot, whose bases then cannot be read, and the carrier leaves the row, as the
    carrier gate reads such a slot (Added). (The integration's first stage read that indel as a
    phasing artefact too, and kept such carriers: on the genome-wide search below, 89 listings
    more, among them the whole row chr5:4507931 (+), whose 57 carriers all hold chr5:4507932
    CTTTCA>C on its slot; on chr22, TRAC's 46177780 (−) listed HG02652, HG02653 and HG04056,
    whose slot of its 46177786 G>C and 46177787 C>T holds 46177784 AGGC>A.) A haploid call --
    a male's on chrX outside PAR -- is his one haplotype's (see the next entry), so an indel of
    his that changes the target takes him off the row, whatever else that haplotype holds. A
    row left with no carrier is not reported. Population-level searches are unchanged (their
    rows are worst-case composites), and so is a search with `CRISPRME_INDEL_SNP=0`, the opt-out
    of the SNP+indel co-occurrence: it does not read the store, and its SNP rows list their
    carriers as 2.8.1's do.
  - README guide, HGDP + 1000 Genomes, NRG, 4 mismatches + 1 + 1, per-sample (2.8.1 output,
    checked against haplotypes rebuilt from the phased VCF): genome-wide, 280 of 65,184
    highest-CFD SNP rows of the integrated results listed such carriers (28,700 listings,
    13,525 of them carriers homozygous for the indel), and 208 had no other carrier (none at
    CFD ≥ 0.2, highest 0.112); 210 of 39,752 alternative alignments (24,145 listings), 152 with
    no other carrier. TRAC guide: 157 of 49,841 rows (30,340 listings), 97 with no other carrier
    (highest CFD 0.049). chr22 alone: 3 rows and 7 alternative alignments. After the fix,
    checked the same way genome-wide: none, but for 2 listings whose indel sits at the position
    of the row's SNP (a phasing artefact, kept as at an indel row's own position).
  - The same search on chr22, released 2.8.1 image against 2.10.0's post-analysis on it (CRISPRitz
    2.9.1, index format 1), per-sample: 1,574 → 1,570 integrated rows and 3,756 → 3,753 alternative
    alignments (without this fix, 1,573 and 3,757). This fix takes off the rows at 26806221 (−),
    27144006 (+) and 30569221 (−), and the alternative alignments at 21052638 (−), 27795808 (−),
    42848613 (−) and 44756100 (+), whose `Other_motifs` go down by one; none of their listed
    carriers has the target. At 21052638 it is 2, as in 2.8.1, not the 3 it is without this fix:
    21052639 (the row 2.8.1 placed at 21052642, now at its insertion's anchor) takes the place of
    the SNP-only alignment of 21052641 G>A. One row (24441929 (+)) and 3 alternative alignments list
    fewer carriers. The summaries by guide and by sample, `observed_haplotypes_raw`,
    `observed_haplotypes.tsv` and `phase_confirmation` follow; the dense-window BED and the
    co-occurrence companions are unchanged. NRR, per-sample: 3,987 → 3,982 integrated rows and
    11,581 → 11,577 alternative alignments (without this fix, 3,985 and 11,583): the same rows and
    alignments go, and two more alternative alignments, 27795809 (−), another alignment of
    27795808's SNP, and 47366716 (−), the highest-CRISPR-Bulge alignment of its site's row, whose
    block now shows the site's reference alignment; 5 alignments list fewer carriers (27782673 (+):
    1 of 17). TRAC guide, per-sample: 2,803 → 2,796 and 7,508 → 7,504 (without this fix, 2,800 and
    7,509): the rows at 22629803 (+), 23029103 (−), 29864083 (−) and 42909550 (−) and 5 alternative
    alignments go (CFD at most 0.031), 3 alignments list fewer carriers, and so do the SNP+SNP
    co-occurrence companion's 6 lines of 46177780 (−) (729 → 726 carriers; its `MinAF_bound` is the
    pair's marginal bound and does not change). Population-level searches give the same tables with
    and without this fix (2,198 → 2,197 rows, by the fix above). Where rows go, the CRISPR-Bulge
    score of another alignment scored in the same batch can move in its last printed digits (NRR
    38558318 (−): 3.56555e-16 → 3.56556e-16; 3 alignments in these searches); that changes no
    block's choice of alignment (the one block that changes, NRR 47366716 (−) above, changes because
    its alignment went).
  - Cost (per-sample searches only): the SNP post-analysis reads the indels of the store near
    each chromosome's SNP targets (839,234 of its 17.2 million genome-wide, 3 to 46 s a
    chromosome, in parallel) and reads a carrier's bases only where it holds one that changes a
    target, once per case. Genome-wide, side by side on one host, the SNP post-analysis took
    14 min 21 s against 13 min 42 s without this fix (+5%), less than the spread between runs
    on that host (9 min 11 s with this fix and 9 min 32 s for 2.8.1 at other times); its
    largest worker 2.03 GB (1.83 GB; a median of +0.1 GB a chromosome); chr22: 43 s (37 to
    38 s).
  - Covered by `PostProcess/test_snp_row_indel_carriers.py` and the SNP+indel oracle, which now
    checks every per-sample SNP row's carriers against their haplotypes and plants sites Y1-Y4
    (cis, trans, homozygous, unphased, a G run at the target's edge, an always-cis SNP). Without
    this fix that check fails 62 times (57 row blocks at 11 sites, and the Y sites'
    expectations) and 6 in the Cas12a profile (the SNP-only rows of C3 and C5), on top of the
    indel-row checks above: the released v2.8.1 image fails the two 166 and 53 times (the whole
    oracle, with the check the merge fix below adds, 168 and 53).
- **A per-sample row on chrX no longer names a male whose one haplotype holds an allele inside
  the target that the row does not.** On the dict-less tiers a male's calls on chrX outside PAR,
  which are haploid, were read as an unphased sample's: every subset of his alleles in a window
  was a haplotype he could have, so a row could list him for an alignment that holds some of them,
  though his one haplotype holds them all and its sequence there is another. His calls in a
  window are now that one haplotype, as the classic dictionaries read them already, and he is
  listed for the alignment it has, or for none. On the genome-wide per-sample search of the
  README guide (HGDP + 1000 Genomes, NRG, 4 + 1 + 1), 23 chrX rows whose listed carriers were all
  such males, none of them with an indel that changes the target, are gone (at chrX:144234310(-)
  all 171 also hold chrX:144234328 T, at chrX:19417344(+) all 85 chrX:19417355 A), 5 keep their
  place with another alignment and 13 list 1,460 carriers fewer: in all 1,924 carrier listings
  leave 51 SNP alignments of the rows' three blocks (this release's search of chrX with and
  without the fix). With this release's other fixes, the search returns 170,455 rows against
  2.8.1's 170,660 (chrX 8,221 against 8,249). Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_haplotype_prioritizer.py`.
- **A `--per-sample` search of an index without the dict-less tiers no longer reports a dense
  window by its greedy worst-case alignment.** On the classic dictionaries a window denser than
  `CRISPRME_IUPAC_CAP` was reported by an alignment that took every allele matching the guide and
  listed the carriers of any of them, so that nobody need carry it: on the self-test's
  `NGG_2_hg38_chr22+hg38_1000G` the README guide (`20bp-NGG-SpCas9`, 4 mismatches + 1 DNA + 1 RNA
  bulge) gave 5 such rows (chr22:22410074, 22687072, 22687116, 22698190 and 22888934), naming 102
  individuals as carriers 114 times. Such a window now gives the haplotypes its carriers have, read
  from the dictionaries as from the tiers: the same search returns 560 rows against 565, the five
  gone. Two rows that started on a base an insertion adds start at its anchor (chr22:17434485 (+)
  and 21052639 (−), which 2.8.1 placed at 17434488 and 21052642: an entry below), and nine others
  differ: chr22:44864861 lists one carrier fewer and chr22:42848613 one other motif fewer (an
  indel in the target, two entries above: every one of the 87 haplotypes that carry the SNP of
  the alignment gone, chr22:42848616 A>T, carries the insertion at chr22:42848614 too),
  chr22:50763084's CRISPR-Bulge score its last digit (3.61989e-21, 3.6199e-21: a score of another
  batch), and six the order of their rsID and AF lists alone (an entry above). Neither
  `complete-search` nor the report says any longer that a dense window can be such an alignment
  (but the report of a 2.8.1 search's results, read again), and neither does a per-sample row's
  `High_variant_density_region` note (see Changed). Covered by
  `PostProcess/test_haplotype_prioritizer.py`, `PostProcess/test_phased_haplotype.py` and
  `PostProcess/test_generate_report.py`.

- **A per-sample row on an index with the dict-less tiers is no longer noted as a dense window's
  greedy worst-case alignment.** Since 2.7.10 that path reports a dense window by its carriers' own
  haplotypes and its reference row, never by a greedy representative, yet
  `High_variant_density_region` said of each of its rows that "a greedy worst-case alignment is
  reported here; additional haplotype alignments may exist". It now gives the window's variants,
  sequences and IUPAC protospacer, and says that the haplotypes its carriers have are reported
  there (the carrier gate's wording, under Changed), as the BED's new `greedy` column says (on
  chr22 of the published index, the README guide's per-sample search: 2 rows). It says greedy
  where that column, or the row, says a greedy worst case stands for the window: a
  population-level row, or a sites-only dataset's, which names no carrier; a BED of 2.8.1, which
  has no such column, keeps 2.8.1's note. The report reads such a row as before: it takes a dense
  window's row for its greedy alignment only when no dataset of the search has the tiers, in the
  results of a search before this release. Its legend said in every report that "such dense
  windows are reported greedily (one representative alignment is shown)"; a per-sample report now
  says how the search reported its dense windows (under Changed), and doubts a worst case in a
  hyperdense window (under Added) only where one stands for a dense window.

- **A search with no bulges no longer builds an index.** It scans the genome's FASTA files, and a
  variant search the enriched genome's when it is here whole, so it opens no index; it built the
  reference index `<PAM>_1_<genome>` all the same (hours on hg38), and a variant one, which nothing
  read, and under `--index-path` it stopped on a library that held no reference index. Its INDELS
  index, which the indel search reads at any bulges, is still found or built as before. The web's
  progress page still reads the indexing step as done. A personal-assembly search with no bulges
  no longer builds a reference index per haplotype (R4-4). On chr22 with HGDP + 1000 Genomes, a
  search with no bulges of `NGCG`, which no installed index holds, built only the index of its
  indels, `20-NGCG_1_hg38+hg38_HGDP1kGP_INDELS`, in 50 s with nothing on stderr, where 2.8.1 built
  `NGCG_1_hg38` and `NGCG_1_hg38+hg38_HGDP1kGP` as well (54 s); its record of the indexes searched
  says `Reference_index none` and `Variant_index:hg38_HGDP1kGP none`, each with why, the second
  naming that index of its indels. The record says `none` for the variant index of such a search
  wherever its enriched genome is whole, an installed index of the search's PAM beside it or not:
  it named that index, with its bulges, as the index searched, and the report's summary with it,
  though only the genome was scanned (every local or Settings build keeps its enriched genome, so
  a web search with no bulges, "Max edits" at 0 bulges, is the common case). On chr22 with HGDP +
  1000 Genomes, the `NRR` indexes installed and the enriched genome in place, a search with no
  bulges of `NRR`, and one of `NGG` which only those indexes serve, each record
  `Variant_index:hg38_HGDP1kGP none: a search with no bulges scans the variant genome; its
  indels: NRR_3_hg38+hg38_HGDP1kGP_INDELS`, and their reports' summary rows say it; each named
  `NRR_3_hg38+hg38_HGDP1kGP`, with its two bulges, before this.
- **A search its input check refuses leaves its output folder empty.** A guide without its
  PAM's `N` padding, for one, is refused before the search starts, but the run had written
  `.command_line.txt`, `.version.txt`, `Params.txt` and `guides.txt` into `Results/<output>/`,
  so the corrected command, with the same `--output`, stopped on `Output folder ... not
  empty!`. Found running this release's chr22 searches.
- **Load Example selects the PAM and the genome it computed.** It picked the installed default PAM
  and then returned "20bp-NRG-SpCas9" and "hg38" whatever was installed (FR :125).
- **The web form's bulge cap read an index built for the other side of the guide as serving.** Its
  mirror of the search script's choice had no side rule (F2): it now takes an index only for a PAM on
  the side its `.pam_build` (or its name) records, as the search does.
- **`crisprme.py setup` exits non-zero when its download fails.** It caught the legacy script's
  failure, wrote `Legacy database setup exited with error code 1` with no end of line, and exited 0
  (with no network, a `requests` traceback and status 0).
- **setup's progress lines come as they are printed.** Into a pipe or a log file, the download
  script's progress lines came at its end, after its errors, so a merged log read out of order.
- **setup's test-data URLs no longer hold `//`.** The base URL ended in `/` and each URL added
  another (`.../test/data//samplesIDs/...`); GitHub answered 307 and the download worked.
- **`generate-report` says why a folder has nothing to report.** On a search that found no
  off-target, and on an assembly-search's combined folder, it ended on a `FileNotFoundError`
  traceback (OF-12). It now says, in one line, that the folder holds no integrated results (and
  why), or that assembly-search writes the combined report there itself.
- **A PAM-less search records its nuclease as `NO-PAM` from the command line too.** complete-search
  took the third field of the PAM file's name, `NO`, where the web reads `NO-PAM` (OF-13).
- **environment.yml installs pigz**, which the Docker image has and `publish-index` packs an index
  with on every core (install P9).
- **The one-click apps start nothing on half the data after a failed first download.** The
  Windows app's first Start ran `download --what all`, the variant index's download and the app
  with `;` between them, so a failed download went on to the next and started the app on half the
  data; each step now runs only when the one before it succeeded (`$LASTEXITCODE`: Windows
  PowerShell 5.1 has no `&&`) (backlog 4.4). And both apps, the macOS one too, chose between
  setting up and opening the app by `Genomes/` alone, which the download makes in its first
  minutes, so the next Start after a failed download opened the app on part of the data. The
  setup now marks the data folder (`.setup_unfinished`) until both downloads have ended well: a
  Start after a failure sets up again, picking up what is there, and the failure says to click
  Start again to finish setting up. Checked with the macOS app's AppleScript under `osascript`
  and `install.ps1` under PowerShell 7, each with a stub `docker` that fails the index download,
  then not; an install of an earlier installer, with `Genomes/` and no mark, opens as before.
- **A search no longer edits the user's `--vcf` and `--samplesID` lists.** The search script
  appended a newline to each, for `while read`, and cut it off again at its end: a run that
  ended early (no off-target, an error, Ctrl-C) left the newline in them, and a read-only
  `--vcf` list failed the write, on stderr, which fails the stage (RV4-1). It reads copies in
  the output folder (`.vcf_list_searched.txt`, `.samplesID_list_searched.txt`).
- **A `--genome` folder or a working directory with a space in its path is refused up front**,
  by complete-search and by build-index-only, each in one line that names the command: the
  search script takes both paths into its commands unquoted, and failed later with a path that
  is no folder's (IDX-19). The working directory is checked first; for the genome folder, the
  line says to move it or to link it from the working directory under a name without one (`ln
  -s '<folder>' <working directory>/Genomes/<a name without spaces>`), which the search takes.
- **The Settings genome upload says what to do with an archive of gzipped FASTA files.** It was
  refused with `no .fa or .fasta file in X`; it now names the `.fa.gz` files it holds and says
  to gunzip them, or to upload one `.fa.gz` a sequence (M3).
- **A PAM file whose first line is not the motif and the length, one space apart, is refused as
  INPUT_FORMATS §4 says.** complete-search checked only that the last field of the first line is a
  number: a line of the length alone made `Results/<output>` and failed the input check later, one
  with a third field (`N...NGG 3 4`) searched with that field for the PAM's length, and a blank
  first line, with the motif on the next, was refused as a file of two motif lines (the focused
  review of f8dbe73). Each gets the one-line refusal of a first line that is not `<motif>
  <length>`, before anything is made.
- **CRISPR-Bulge's scoring leaves nothing in the temporary folder.** TensorFlow's autograph writes
  the code it generates to a module in `$TMPDIR` and removes it at exit, but importing it left
  `__pycache__/__autograph_generated_file*.pyc` there, one a run: in the data folder itself under
  `CRISPRME_TMPDIR=/DATA` (the focused review of f8dbe73). The scorer's worker writes no bytecode
  now (`PYTHONDONTWRITEBYTECODE`, as a user's own setting is kept): a read-only container's search
  with `-e CRISPRME_TMPDIR=/DATA` leaves `/DATA` as it found it, its results those of 2.8.1.
- **A Settings Direct URL download is laid out as an upload is.** A `.fasta` file goes in as a
  `.fa` file, and the hidden files of an archive are left out: a `.fasta` file that needed
  nothing else stayed so, which no index or search reads, and the `._chr1.fa` that macOS `tar`
  adds beside `chr1.fa` went in, with which a search stops (plan 4.8, R6B
  direct-url-fasta-and-hidden-files). `chr1.fa` and `chr1.fasta` in one archive are two of one
  name, refused as an upload refuses them. docs/SETTINGS_DATA_MANAGER.md said so as a limit.
- **RELEASE_CHECKLIST's HPC route runs from a clean host environment and a fresh sandbox.**
  Apptainer passes the host's environment into the container (a login shell's `CONDA_*`
  reached it on ml008), and steps 2-4b reused the writable sandbox step 1 had left its caches
  in (list_282.md).
- **The web results page no longer warns of "Multiple perfect matches" for a search of several
  guides with one perfect match each.** It counted the perfect-match sites of every guide
  together, so two guides with an on-target each got the red banner that no on-target can be
  told (RV1-F5). It counts them per guide, as the report's banner does, and names the guide.
- **A 5′ search refused for want of an index names the PAM-less index whose flank is too wide
  for it.** A search skips, without a word, a PAM-less index whose flank is wider than a 5′ PAM
  (the engine reads a 5′ PAM off a flank of its own width only), so the refusal named none of
  the indexes there, or blamed the index's recorded side (cli5 D7). The refusal under
  `--index-path`, the search script's refusal of a variant search, and complete-search's
  pre-flight now say, for example, that `NNNNNNNN_2_hg38 is a PAM-less index whose flank is 8
  nt, and TTTV is a 5' PAM of 4 nt`, and that the search needs `NNNN_2_hg38` or `TTTV_2_hg38`.
- **docs/DYNAMIC_PAM_INDEX_COMPAT.md says what the release does.** Its CLI, web and PAM-picker
  sections called planned what 2.10.0 now does, and described a CLI that took any PAM inside an
  index's own and passed `--prune-pam` by itself, which 2.8.1 did not (IDX-22); they now say
  what ships and how it differs from the design, and that a search no installed index serves
  builds an index of its PAM where it can, as before (it said every one was refused). METHODS
  §9 no longer says that one index answers one PAM, nor that a search takes no all-`N` index of
  another width without the flag.
- **validate-benchmarks' read-only smoke can fail.** It made the code folder read-only with
  `chmod -R a-w` in a container that runs as root, which writes through permission bits with
  `CAP_DAC_OVERRIDE`, so a write into the code folder passed and the step could not fail (CI-1).
  Its container drops the capability (`--cap-drop DAC_OVERRIDE`) and works in a folder of its
  own; RELEASE_CHECKLIST T-38 says so. Run on `pinellolab/crisprme:v2.8.1` with this tree: a
  search script planted to write into its own folder fails the step (`Permission denied`), and
  passes it without the flag.
- **A Settings genome upload of another kind of file is refused on its first chunk, and every
  refusal of an upload is answered as the client's.** A file that is not FASTA, an NCBI
  Datasets `.zip` for one, was taken for a FASTA file and refused once whole, for its first line,
  the archive's bytes, with HTTP 500 `finalize failed: ...`, as was every refusal of what an
  upload holds (WEB-9). An upload whose name is not `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`,
  `.tar.gz` or `.tgz`, or whose first bytes are not what its name says (gzip, or FASTA text), is
  refused before the rest of it streams in, and a refusal at its end (an archive of gzipped FASTA
  files, a truncated gzip, a VCF upload into a registered folder) is answered 400 with its reason;
  so is a personal-assembly bundle that is no archive (`not a gzip file`). A failure of the
  server is still 500.
- **validate-test downloads the brute-force references of its own release.** They came from the
  main branch, and are md5-checked against the registry the release ships, so the day main
  regenerated one, validate-test of every installed release (Bioconda leaves the references out
  of its package) failed (C11/P7). They come from the release's tag (`v<version>`, which
  `crisprme.py validate-test` passes on), or from main for a version with no tag yet.
- **`PostProcess/bench_cbulge.py` does nothing on import.** It changed the working folder and the
  module path of whatever imported it, to a developer's path, or failed to import where that path
  is not (RV4-3); it takes the CRISPR-Bulge repository from `--repo` or `CBULGE_REPO` when run.
- **The report's "Panel: N individuals" counts the samples the search ran with.** It counted every
  list in `samplesIDs/`: the datasets not searched, and the lists Settings writes (M2).
- **build-index-only with no `--samplesID` no longer calls the tiers in place missing.** Run as a
  search's ERROR says to run it, to rebuild an INDELS index, it left the dict-less tiers of the
  build before as they were, and said, in a loud WARNING on stderr, that the index was MISSING
  them, INCOMPLETE (emptyidx F6). It says, on stdout, that they are kept; with none in place, or
  ones a stopped run left unfinished, it warns as before.
- **The web status page says why the search was not run when an INDELS sub-index is empty or
  missing, or the genome folder is not as a search needs it.** Each got the general "remove it and
  try to submit again", though a new submission fails the same way (emptyidx F7, and the genome
  layout refusal of the final 2.8.1 review); the banner now gives the ERROR, and for the genome
  folder what `genome_layout.py` says of it, with the command that fixes it.
- **Settings shows an index's display name.** The Installed indexes list and the Remove card named
  an index by its folder alone, never by the display name given it (`build-index-only --name`, or
  the Display name of a build there), which the search form shows (P15.3): both show it beside the
  folder.
- **Six function-based test modules run under `python -m unittest`.** Each reported `Ran 0 tests`
  and `OK` there, a false pass (suite P4); a `load_tests` hook collects every test function. CI
  and the release checklist still run them as scripts too.
- **Four web test modules skip their classes where Dash is not installed.** They raised SkipTest
  at import, which `python -m unittest <module>` reports as a traceback and an error (skips D7).
- **With no `--bDNA`/`--bRNA`, a search under `--index-path` takes the bulges that library
  serves.** complete-search chose 1 bulge of each type when the genome's FASTA files were there
  to build a reference index from, or the working directory's library held one, though under
  `--index-path` a search reads the reference index from that library alone and builds none
  there: given a library with no index of the PAM, it then stopped on `no matching reference index
  under --index-path` (plan 4.8 RV2-5). It searches with no bulges, which needs no index, and its
  WARNING says to install an index of the PAM in that library, or to run without
  `--index-path`.
- **A search leaves no hidden `.binview/` folder in its results.** The engine reads an index
  folder that holds `.bin` files alone, so a search hands it a folder of links to the trees of an
  index with a sidecar (a download's, and since 2.8.1 every index a search or `build-index-only`
  builds), which stayed in the result folder (RV3-4, RV3-6). It goes as each search ends, whether
  the search found targets or failed.
- **Submit makes no job folder for base-editing values the form does not offer.** With base
  editing on, a window bound that is no whole number, or a nucleotide that is no base, was refused
  (by `int()` and an `assert`) only after the job folder and its files were made, so each such
  request left a folder that never ran, as a Max edits value did before 2.8.1 (secfu RV1-N1). The
  values are checked with the others, first.
- **The Personal assembly tab offers every PAM installed and the app's bulges, whatever the
  Reference tab's genome and variant dataset.** Its PAM list, its bulge options and the note
  under them followed the hidden dropdowns of the Reference tab (OF-11): with a VCF dataset
  selected there that had no index, it offered 0 bulges only, with a note about that dataset's
  index, and with an NRG variant index, NRG alone, though an assembly search indexes each
  haplotype itself; its Submit took the same counts. It offers every installed PAM and 0 to 2
  bulges of each type, its note says what its search builds, and the PAM picker of a PAM-less
  index is not shown there.
- **Settings' Build an index locally pairs a VCF dataset with its genome.** It offered every
  genome for a dataset whose genome Settings recorded when the dataset was added, and a build of
  the dataset on another genome ended in the enricher's crash (`std::out_of_range`, SIGABRT;
  FV :328). Choosing such a dataset selects its genome, the one the form then offers, and a
  build of the dataset on another is refused, naming its genome. The form reads that record for
  a dataset it offers alone: a request can send any value. Covered by
  `PostProcess/test_web_bulge_notes.py` and `test_web_shell_safety.py` (in the image).
- **A job an identical search is sent to keeps its log.** Submit gave the finished job's "Job
  Done" line today's date (its results are kept 3 days from it) by writing the log back
  without its last line, the report's, and with a second "Job Done" and no line end
  (FV :335); it changes the date of that line alone.
- **The search form offers a variant dataset installed as an index alone under its own genome,
  and by its name.** It read the dataset off the index's name without the genome before the
  `+`, so an index of another genome, of a dataset named after this one, was offered too
  (FR :278), and it showed each underscore of the dataset's name as " + ", as for a merged
  panel of two datasets (`my + vcf` for `hg38_my_vcf`; FR :286).
- **A Ctrl-C ends a command in one line.** crisprme.py ended in a `KeyboardInterrupt` traceback
  (install P6); it says `crisprme.py <command> interrupted (Ctrl-C)`, for complete-search with
  the folder the run wrote into, and still ends by the signal (status 130 in a shell). The
  pipeline's own stages, which the same Ctrl-C stops, leave theirs in the run's
  `log_error.txt`, as before.
- **A source or Bioconda install's self-test prints no warning the image does not.** Every
  anonymous download of the public data printed "You are sending unauthenticated requests to
  the HF Hub", a server warning that huggingface_hub logs and the Docker image hides
  (`HF_HUB_VERBOSITY=error`): crisprme.py sets the same level, unless the variable is set
  already. And `install_from_source.sh` printed git's advice about a detached HEAD as it checked
  out the CRISPRitz tag (install P10).
- **`--sorting-criteria-scoring` is used.** The merge of close targets chose the criteria of the
  CFD and CRISPR-Bulge merges by a variable nothing set, so every merge took
  `--sorting-criteria` and `--sorting-criteria-scoring` did nothing (wset N1). The CFD and
  CRISPR-Bulge merges take it now; its default is what they took, `mm+bulges,mm` (it said
  `mm+bulges`), so a search with neither flag, and every web search, merges as before. On chr22
  with HGDP + 1000 Genomes (three guides, NGG, 4 mismatches + 1 + 1 bulge, per-sample), a search
  with neither flag wrote the same 2,023 rows, byte for byte, as the code before the change,
  and one with `--sorting-criteria-scoring bulges` 82 of them otherwise; both with nothing on
  stderr.
- **build-index-only refuses, before it enriches anything, a VCF the enrichment cannot read.**
  `crispritz add-variants` takes each VCF's chromosome from the last part of its file name that
  holds `chr` and reads the genome's FASTA file of that name, so a multi-chromosome VCF whose
  chromosomes are named without `chr` (`1` on a `chr1` genome, `Contig1`, `NC_000913.3`), which
  the build split under those names, a VCF named so by chance (`nochr.vcf.gz`), and one of a
  chromosome the genome has no FASTA of (another genome's dataset) each ended in `ERROR: cannot
  open <genome>/<name>.fa` and "Genome enrichment (add-variants) failed", the split one after
  splitting the dataset in place (wset N2). Each is one ERROR that says why and what to do, and
  the dataset is left as it was.
- **The Settings uploader shows a refusal in red.** An upload refused (`Error: ...`), one that
  failed, and "Choose a file first." were grey, as the progress and the success are, so a
  refusal read as if all were well (final confirmation of 2.8.1); they take the red of the
  Settings cards' messages.
- **The advice for a merged panel's missing sample list names the list that is missing.** Submit
  makes a merged panel's list (`hg38_1000G_HGDP`'s) from the lists of the datasets it merges;
  with one of them missing, and no VCF folder of the panel, the note and Submit said to add the
  VCF dataset `hg38_1000G_HGDP` again in Settings, a dataset never added (FR :98). They name the
  missing list (`samplesIDs/hg38_HGDP.samplesID.txt`), and how to get it.
- **A genome folder holding a FASTA file of one sequence named otherwise is refused before the
  search.** A search takes each file for the sequence it is named after, so it found the
  off-targets of `walkone.fa` holding `>chr1`, dropped them all, and ended saying that none
  survived post-analysis and to check the folder with `genome_layout.py`, whose check, which
  the search runs first, passed the file; complete-search's input check only warned that
  results "may be dropped" (M1). The check, the search that runs it first (a web search too)
  and the input check now refuse such a `.fa` file, saying to name it after its sequence with
  `genome_layout.py --split`, as build-index-only already did; a `.fasta` file, which no
  search reads, is still only named.
- **`crisprme.py web-interface --port <N>` serves on another port.** The server listened on
  8080 alone, so where another program held it, as on a shared node, an Apptainer or source
  install could not start it (Docker maps another with `-p`; C19). Its start message, and the
  link of a notification e-mail, name the port.
- **`docker stop` stops a CRISPRme container at once.** Its first process, crisprme.py, had no
  handler for the SIGTERM a stop sends, so `docker stop` waited its 10 s and killed the
  container (exit 137; C16): the image's stop signal is SIGINT, which crisprme.py takes as a
  Ctrl-C (the web interface stops its server and exits 0). Measured with the 2.8.1 image, run
  with each signal: 10.2 s against 0.5 s.
- **The Docker image is built on a pinned base.** `FROM mambaorg/micromamba` took whatever image
  was latest the day of the build (C20); the Dockerfile names the tag and the digest of the one
  2.8.1's image was built on (micromamba 2.10.0, Debian 13).
- **A search refused for want of a reference index under `--index-path` names an all-`N` index
  there that records no side, and how to use it.** Such an index of a width but 3 and 8 is not
  taken for a 3′ PAM (a search before 2.8.1 built one of any width from a 5′ PAM file), and
  `build-index-only` records no side over it, so one really built for a 3′ PAM could only be
  rebuilt (R6B RV1-6); without `--prune-pam` the refusal did not name it at all. It names it
  either way, and says that a name of the 2.10.0 form records the side: `mv
  .../NNNN_3_hg38 .../20-NNNN_3_hg38`, for the guide length its trees were built for, after
  which the search takes it.
- **complete-search stops at once when no temporary folder can be written.** With none
  (`docker run --read-only` without `--tmpfs /tmp`), a search ended with exit 0, CRISPR-Bulge
  scores of -1, no report zip and a WARNING in its error log alone (final confirmation of
  2.8.1; TROUBLESHOOTING §1 said so). It says, before it makes anything, which folders it
  tried, and to set `CRISPRME_TMPDIR` to a writable one or give the container a `/tmp`. A
  folder set with `CRISPRME_TMPDIR` or `TMPDIR` is taken, the working folder included (`docker
  run --read-only -w /DATA -e CRISPRME_TMPDIR=/DATA` runs to its report, its rows those of
  2.8.1); the working folder is not taken when nothing names it, as the scorer and the report
  do not run in it.
- **complete-test is not stopped by the empty results folder a refused search leaves.**
  complete-search makes its output folder, and a refusal of its pre-flight leaves it empty: the
  next complete-test stopped on it, "Complete-test already run once", until the user removed it
  (R2-F2). An empty folder holds no earlier run's results, and complete-search takes one: only
  a folder with something in it stops complete-test now.
- **The base-editing window of BED lines answers no server error.** The form recomputes the
  window's options at each change of the guides, and for genomic coordinates it raised with no
  genome selected (a fresh data folder after Load Example), on a chromosome the genome has no
  FASTA of, on coordinates past its end and on a line short of its three fields: each an HTTP 500
  and a callback error in the browser's console (FR :89, :189). It offers no window for them, as
  for any line Submit then refuses, saying why.
- **`complete-test --chrom` fetches a chromosome but chr22 from Hugging Face too.** The data repo
  holds every hg38 chromosome gzipped (`genomes/hg38/<chrom>.fa.gz`) and chr22 unpacked as well,
  and complete-test asked for the unpacked file alone, so every other chromosome came from UCSC
  (FR :725). It takes the gzipped one from the Hub when the unpacked one is not there, and UCSC
  stays its source when neither is: chr21 came from the Hub, 47,644,190 bytes, as from UCSC.
- **A FASTA whose last line has no line end is refused, and `genome_layout.py --split` ends it
  with one.** The index build drops a last line with no LF, with no error: a sequence written on
  one line lost all its off-targets, and a 60-wrapped chr1 whose last line held a guide's
  on-target gave no row with 2.8.1, and a search with bulges failed (no `.bin`); `--split`,
  which every genome download and Settings upload runs, deleted lone CRs past the first lines,
  which left one long line with no LF (LP unterminated-last-line-dropped). The check a search
  runs first, complete-search's input check and build-index-only refuse such a file (one byte
  a file), saying to run `--split`, which writes a CR between bases as a line end and ends
  every file it writes with LF; after it, the same search found the on-target, with and
  without bulges.
- **The advice for a dataset with no sample list names a VCF file that cannot be read.** A
  Settings build writes the list from the header of every VCF file of the dataset, and writes
  none when one has no readable header (a file cut short, or no VCF); the note and Submit said
  all the same to build the index in Settings, "which writes the list" (M6). They name the file,
  and say to add the dataset again with it whole, or to add the list.
- **A variant search under `--index-path` whose library holds no reference index it can use
  stops before it enriches the genome.** The search script stops on such a library, but only
  after its first step, the enrichment of the genome with the dataset's variants, hours on
  hg38 (IDX-17 e): complete-search stops first, in one line, when that step would run (the
  dataset's VCFs here, no enriched genome, no variant index that serves the search). A
  reference search, which reaches the script's refusal at once and says what the library holds,
  is left to it, and so is a library where a run is building an index, which the script waits
  for. It names the index to download as the script's refusal now does (below).
- **A search refused for want of a reference index under `--index-path` names one to download
  that is published.** The refusal said to download `<PAM>_<N>_<genome>` with N the search's
  bulges + 1: `NRR_2_hg38` for a search of 1 DNA and 1 RNA bulge, the default, where every index
  published is of N = 3, so the download it advised stopped on `No index 'NRR_2_hg38' found`. It
  names the one of N = 3 for up to 2 bulges (`NRR_3_hg38`), "if it is published", as the refusal
  for want of a variant index does, and its other route, to run without `--index-path`, names
  what such a search does since 2.10.0: it takes an index of the working directory's library,
  pruned to the PAM where only such a one serves it, and builds one from the genome where none
  does (it said that it builds one). On chr22 with HGDP + 1000 Genomes, both refusals of a search of
  1 + 1 bulge under an empty `--index-path` name `NRR_3_hg38`, which the Hub holds, and the Hub
  holds no `NRR_2_hg38` (`download` stops on it), nor an `NGCG_3_hg38`. Covered by
  `PostProcess/test_index_resolution.py` and `test_bulge_cap.py`.
- **A variant index of `--index-path` whose dictionaries are not in the working directory is
  named.** complete-search said "No prebuilt variant index found" for a dataset whose variant
  index the library held (IDX-17 c): it says that the index is there, that a search reads the
  dataset's dictionaries, sample list and VCF folder from the working directory, and to run the
  search from the folder that holds them.
- **`crisprme.py` with no argument makes no folder.** It prints the help, and made the seven
  folders of a data folder (`Annotations/`, `Dictionaries/`, `Genomes/`, `PAMs/`, `Results/`,
  `VCFs/`, `samplesIDs/`) wherever one ran it to read it (C14); README §1.1.2, which said so,
  says that it makes none. Covered by `PostProcess/test_complete_search_flags.py`.
- **The web form names the VCF file a search cannot read in a folder of several.** A search of a
  dataset with no variant index for the PAM prepares the variant genome from its VCF folder, which
  takes each file's chromosome from a `.chrN.` part of its name, and the note under the bulges and
  Submit refused a folder with any file named otherwise, saying that "the dataset's VCF is a single
  multi-chromosome file", also of `walk.chr1.vcf.gz` and `walk.chr2.vcf.gz` beside `merged.vcf.gz`
  (RV1-WEB2-5). They say that the folder holds `merged.vcf.gz`, whose name says no chromosome; a
  dataset of one such file is refused as before. Covered by `PostProcess/test_web_bulge_notes.py`
  (in the image).
- **The docs say what a run writes.** The web user guide defined an ALT row (`REF/ALT_origin`) as
  one "created by a variant allele", which the results page's caption contradicts: it is a variant
  alignment, of a site a variant creates (`Not_found_in_REF` `y`) or of a reference site it
  changes (RV1-WEB2-3). METHODS put targets-scoring's CFD-only banner "in its log", which it
  prints on stdout, where a search's `log_verbose.txt` holds it (D3R3-3). The README and
  complete-search's help said that "the CRISPR-Bulge columns read -1" in a CFD-only run, whose
  CRISPR-Bulge risk score, ALT less REF, reads 0 (D3R3-2). And no text said that a line of the
  SNP+indel companion (`indel_snp_cooc.tsv`) lists the first 50 of its carriers by name, which
  its count takes all of: the README and METHODS say so (LP RV4-cooc-carrier-cap). Covered by
  `PostProcess/test_defaults_and_index_flex.py`.
- **The Graphical Reports selector offers the totals of mismatches and bulges the search has.**
  It offered 0 to the mismatches + `Max_bulges`, an index dimension (max(DNA, RNA) + 1 as the
  web writes it, max(DNA, RNA) as complete-search does), and opened on the largest, where the
  search draws the population plots of 0 to its mismatches + DNA + RNA bulges (WEB-4): a web
  search of 4 mismatches and no bulge opened on 5, which has no plot ("No result found for this
  combination of mismatches and bulges"), and a web search of 4 + 2 + 2 could not select its
  sites of 8, nor a complete-search of 4 + 1 + 1 its sites of 6. It offers 0 to mismatches +
  DNA + RNA, up to the cap on them together (Max edits, `--max-total-edits`), and opens on that; a
  job written before the DNA and RNA lines of its parameters is shown as before. On chr22 with
  HGDP + 1000 Genomes, a complete-search of 4 + 1 + 1 drew the plots of totals 4 to 6, and its
  selector offered 0 to 5 and opened on 5: it offers 0 to 6 and opens on 6. Covered by
  `PostProcess/test_web_results_links.py` (in the image).
- **The report of a search of two or more variant datasets names them all.** A search of a
  `--vcf` list of two names its results file after the last dataset alone, and its report
  named that one in its header, inputs card and run manifest, with that dataset's registry
  panel and counts (RV1-F1). The report reads the datasets from the copy of the list the
  search keeps in its result folder: it names each (`1000G, hgdpx`), credits no one registry
  with the run's counts (its panel is the samples the search ran with, each under its
  dataset), and says that each variant's frequency is that of its own dataset's panel, as
  each dataset is searched with its own. A search of one dataset, or a result folder of an
  earlier release, is reported as before: `generate-report` of a chr22 search of HGDP + 1000
  Genomes gives the same labels and registry counts with its list of one dataset, and names both
  datasets, crediting no registry, with a list of two. Covered by
  `PostProcess/test_generate_report.py`.
- **The report of a search with no variant dataset fits a personal haplotype too.** It said
  `none (reference genome only)` and "Sites in the reference genome", also in the report of
  each haplotype of assembly-search, a search of a personal genome (RV1-F2): it says `none
  (no variant dataset searched)`, and "Sites in the genome searched (hg38)", with the name the
  search gives the genome: `WALK01_paternal_575e5055` in the paternal report of an assembly-search
  of the walkthrough's bundle. Covered by `PostProcess/test_generate_report.py`.
- **The report of a search that computed no score says so in place of its score tables.** A
  Cas12a search, or one of any PAM on the 5′ side, computes neither CFD nor CRISPR-Bulge, and
  its report printed the risk card's CFD bands and Section 4's CFD tiers, each 0, and the
  CRISPR-Bulge and CFD sentences of the recommended panel: read as "no site at risk", where
  the panel was chosen by mismatches and bulges alone (cli5 D6). The card and Section 4 say
  that neither score is computed for the nuclease and that the panel is chosen by mismatches
  and bulges, and the panel's note says how; assembly-search's combined report too. A scored
  search's report is unchanged. A Cas12a (`TTTV`) search of chr22 with no bulges reports so.
  Covered by `PostProcess/test_generate_report.py`.
- **`download --what index` brings again the dictionaries an installed variant index has
  lost.** A search of a variant index whose dictionaries were gone
  (`Dictionaries/registry_<dataset>` or `dictionaries_<dataset>`, and `log_indels_<dataset>`) said
  to download the index again,
  and complete-search's check said to download it, but the download took the index for
  installed, as its folder, its INDELS index, the reference index and the genotype store were
  there, and did nothing (FR :163). It counts the dictionaries a search reads as part of the
  install, and brings them with the rest of the index's archive. The search's advice names the
  command, and says that in the web app, whose Settings page downloads no index that is
  installed, the index is removed there first. With the Hub offline, a download of
  `NRR_3_hg38+hg38_HGDP1kGP` over an install of it without the dictionaries said "already
  installed" before, and goes to fetch it now. Covered by `PostProcess/test_crisprme_hf.py` and
  `test_index_resolution.py`.
- **assembly-search searches a haplotype again when its FASTA files changed since its
  result.** A run of the same `--output` reused a finished haplotype search whose parameters
  and genome folder were the run's, whatever the folder held: after a bundle was imported
  again under the same path with other FASTA files, the result of the old sequence was
  reconciled (R4-6). Each haplotype's result records the genome it was searched from (the
  folder and the fingerprint of its FASTA files, as an index records it), and a run reuses it
  only while the genome is that one; a result of CRISPRme 2.8.1 or earlier, which records
  none, is searched again. An assembly-search of the walkthrough's bundle, run again after a base
  of the maternal on-target was changed, reused the paternal result and searched the maternal
  haplotype again (one perfect match, not two), where the code before reused both. Covered by
  `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search keeps no INDELS index it built while build-index-only put the SNP codes on the
  fake-indel genome.** build-index-only marks that genome unfinished only just before the
  first FASTA it changes, so a search past its check of the mark before then built its INDELS
  index from the FASTAs not coded yet, and every later run kept it, blind to a SNP beside an
  indel on them (R2-F1). An INDELS index during whose build a FASTA of the genome changed, or
  after which the mark is there, is not put in place: the search stops, saying to let that
  run finish and to run the search again, which builds it from the coded genome. Covered by
  `PostProcess/test_index_resolution.py`.
- **The top-1000 plots of a run with no score name no perfect match they do not star.** Where
  the score was not computed (a nuclease of no SpCas9 family, or a run without the scorer),
  its column holds -1, and the stars of the perfect matches were drawn there, below the axes,
  while the legend listed them (wsearch N3). A perfect match is starred, and in the legend,
  only at a score of its own; the plot of the count of edits stars them as before. The
  results page's plots and the Personal Risk Cards' alike: the CFD plot of the Cas12a search above,
  whose one site is its on-target, names no perfect match. Covered by
  `PostProcess/test_top_plots.py`.
- **The web's per-sample option, its word on a sites-only dataset and the report's advice say
  what a web user can do.** The search form disabled Per-sample for a genotyped dataset whose
  index is not built yet (listed as "genotyped, samples not indexed") saying that it "needs a
  genotyped panel"; the note under its bulge options, and Submit's refusal, said of a sites-only
  dataset that only an allele-frequency registry reports its variants, "as the shipped hg38_mega
  index has (docs/DESIGN_mega_index.md)", a file of the repository a web user cannot open; and
  the report of a population-level search said to re-run it with `--per-sample` alone, a flag a
  web search cannot pass (FV :342), as it named commands alone to get the tiers of a dataset
  without them, and to bundle what its slim zip leaves out (the variant-created list and subset,
  the complete raw integrated results, the windows near the top-ranked sites, an
  observed-haplotypes companion over the size gate). The option
  says that the dataset's samples are not indexed yet, and to build its index in Settings; the
  note and the refusal name the hg38_mega index published on HuggingFace, which Settings
  downloads (under Download a prebuilt index); the report names the web form's Per-sample
  analysis beside the flag in both sentences that tell the reader to re-run the search per-sample
  (the Search mode row, and the caption of a population-level search's variant combinations),
  Settings' Build an index locally and Download a prebuilt index
  beside the commands that get the tiers, and the form's Include the full report beside
  `--report-full` in every sentence that names it (the website's example report with them).
  The report of a population-level search of chr22 with HGDP + 1000 Genomes holds both
  per-sample sentences and all four `--report-full` ones, each with the form's option beside
  the flag. Covered by
  `PostProcess/test_web_shell_safety.py`, `test_web_bulge_notes.py` and
  `test/web/test_settings_vcf_datasets.py` (in the image), and `test_generate_report.py`.
- **build-index-only refuses an INDELS index in place that lacks the sub-indexes it builds.** It
  said "Indels index already present", then "Variant index ready" at exit 0, of any INDELS index
  under its name: one copied in part, or one where a search before 2.8.1 had put the sub-indexes
  of a PAM-less index it served (`NNN_2_fakechr1` in `NRG_2_<genome>+<dataset>_INDELS`), each of
  which every search with bulges refuses (the index audit's leftover INDELS folder, plan 4.5).
  It now looks for the sub-index of each chromosome whose fake-indel genome holds indels (none is
  due for a contig with an underscore in its name, or a genome of headers alone), and stops with
  an ERROR that names the ones missing, what the folder holds, and how to build it again. With
  the real engine on a synthetic dataset, a run over that folder of `NNN_2_fakechr1` and
  `NNN_2_fakechr2`, or with `NRG_2_fakechr2` gone, ended 0 before and ends 1 now, and a search
  with bulges refuses the same folder; once it is removed, the build makes it whole. On chr22 of
  HGDP1kGP, the published `NRR_3_hg38+hg38_HGDP1kGP` with its `NRR_3_fakechr22` renamed
  `NNN_3_fakechr22`, or gone, is refused the same way, as a search with bulges refuses it, and is
  present again once named back. Covered by `PostProcess/test_build_index_only_atomic.py`.
- **CI keeps the logs of a job its time limit cancels, and the checks of what the workflows run
  read no comment and see every skip.** validate-benchmarks uploaded its logs on `failure()`
  alone, so a job that up to three tries of every benchmark took past its 120 minutes, which
  cancels it, kept no log of what failed (R6B cigate F3): it uploads them on `cancelled()` too.
  `test_ci_workflows`' checks that every test module is named in a workflow, and that each
  function-based one runs as a script, read the workflows' comments too, where unit-tests.yml's
  paths comment names a dozen modules: with the line that runs `test_dash_output_props` gone, the
  first still passed (R6B cigate F2). Its check that the light job runs the flag-gate scans looked
  at a skip of their classes alone, and passed with one of their tests skipped (R6B cigate F1).
  Both fail on those now. And `test_history_page_assembly` reads the page with its file closed
  (it printed six ResourceWarnings: FFR RV1-WEB2-6). Covered by `PostProcess/test_ci_workflows.py`.
- **The docs say what the self-test's opt-out saves, where an empty `_INDELS` comes from, and how
  common two indels in cis are.** The README, the benchmark README and the registry called
  `CRISPRME_SKIP_HEAVY=1` an opt-out "for very small machines", which reads as one for a machine
  short of memory: the run's peak is the variant enrichment, which every benchmark needs, so it
  saves time on a slow machine, not memory, as they now say (R6B cigate F4). RELEASE_CHECKLIST's
  archive gate (T-40) said that a present-but-empty `_INDELS` is always a broken build; a dataset
  whose every indel is on a contig with an underscore in its name gets one from a build that ends
  well, and the gate names it beside the broken build, neither to be published (R6B emptyidx
  RV1-2, RV2-2, R3-1). FUTURE_WORK §3 called two or more indels in cis common, "~17–19% of real
  indel off-targets": it says that this is the raw count, mostly repeats, and gives METHODS §8's
  ~1–2% of indel loci after repeat-masking (plan 4.1). Covered by
  `PostProcess/test_ci_workflows.py` and `test_defaults_and_index_flex.py`.
- **A carrier is no longer taken off a row for an indel whose FILTER the search skips.** The
  phased indel genotype store (`Dictionaries/indel_genotypes_<vcf>/`, `build_indel_genotypes.py`)
  streamed the VCF without reading FILTER and kept every indel of it, while CRISPRitz's enricher
  skips each record whose FILTER is neither `PASS` nor `.` (those exact strings): no fake-indel
  window, dictionary or tier holds such a record. The indel rows' carrier resolution (a second
  indel inside a target, an overlapping one) and the SNP rows' (a carrier whose haplotype holds an
  indel in the target), and the carrier gate, read that store, so on a VCF with records of another
  FILTER a carrier was taken off a row for an indel no window holds: on the SNP+indel oracle's
  world with the second deletions of its sites O and Q and the insertion of Y1 marked `LowQual`,
  searched through the real CLI with CRISPRitz 2.10.0, a store of every FILTER took S21 and S39
  off the indel rows of O and Q in both modes, and S70 and S72 off Y1's per-sample SNP row (the
  real-data indel check reports the two indel rows' carriers, 14 and 6 gating findings). The
  store now holds only the records the enricher reads (on that world the rows list S21, S39, S70
  and S72, and the check passes in both modes), and its `.done` marker says so after the record
  count (`FILTER	PASS,.`). A store written before has no such line: `build-index-only`, and a
  search that enriches a dataset itself, write it again where the dataset's VCFs are, and a search
  that reads one says so in its `log_verbose.txt`, once per chromosome (`NOTE [indel-snp]: the
  phased indel genotype store of <vcf> on <chrom> was written before CRISPRme 2.10.0 ...`), and
  reads it as before. Every panel published (HGDP + 1000 Genomes and HPRC: FILTER `.`; the
  self-test's 1000 Genomes chr22: `PASS`) holds no such record, so their stores are unchanged.
  Covered by `PostProcess/test_indel_genotypes.py`, `test_search_indel_snp.py`,
  `test_build_index_only_atomic.py` (a store of 2.8.1 written again) and
  `test_indel_overlaid_rows.py` (the note, with the real `analisi_indels_NNN.py`).
- **An indel written on two VCF lines of one REF>ALT is one allele: its carriers on either line
  are its carriers** (decision B2; every release with the phased indel store). The store keeps a
  line per record, and its three readers (the indel path's, the SNP rows' and the carrier gate's)
  kept the later line's carriers alone, so a carrier only the earlier line listed lost its phased
  genotype: CRISPRme listed it on no row of the indel, or on some of its alignments only, and a
  carrier listed on both lost the slots only the earlier line gave. HPRC chr22 writes 39 indels so
  in one case (and 140 more in two cases: next entry), with 1,400 carrier listings on the earlier
  line only: chr22:36235742 AGGGT>A, 144 carriers on its
  first line and 4 on its second, was listed for one sample, and the 14 carriers of chr22:32763933
  GC>G on its first line alone were missing from its 24-nt alignment's row (a fresh per-sample
  search of three guides over such deletions, 4 + 1 + 1: 300 gating findings of the real-data
  check; the same search with this fix gives none of them, its 6 left being an MNV on one
  carrier's slot at chr22:46469988, which CRISPRme does not read beside a record).
  The readers now merge the two lines as CRISPRitz merges two listings of one SNV: every carrier of
  either kept, a carrier of both with its two calls merged slot by slot where they are phased
  alike, else with the call of more copies, the earlier on a tie (`build_indel_genotypes.
  merge_store_carriers`). On the review's world WD (a deletion on two lines, beside the SNV its
  perfect target needs) the rows list S02 and S03 (on the first line alone) and the check passes in
  both modes (17 and 3 gating findings before, with the check's own misreading of the two
  lines). No store need be written again. Covered by
  `PostProcess/test_indel_genotypes.py`, `test_snp_row_indel_carriers.py` and
  `test_haplotype_prioritizer.py`.
- **An indel or MNV record is read with its REF and ALT in upper case, as VCF 4.x reads it: two
  lines of one indel written in two cases are one allele, and the carrier gate finds an indel
  written in lower case** (every release with the fake-indel genome; HPRC). CRISPRitz cuts a
  fake-indel window for each line and its indel log names each line as the VCF writes it, while
  the windows hold one sequence (it writes the ALT in place of the REF whatever their case, and the
  search reads the window in upper case). CRISPRme named each window's rows as the log does, and
  its readers of the phased indel store keyed each line as written. So two lines of one indel in
  two cases (HPRC's `vcfwave` writes 140 indels of chr22 so, chr22:23541792 TG>T beside tg>t) were
  two records: each window's rows listed its own line's carriers, and the carriers of one target
  were split over two rows named in two cases (8,755 carrier listings on one of the two lines
  only, 917 carriers on both), which the real-data check, reading the two
  lines as one allele, failed: on a fresh per-sample HPRC chr22 search of two guides over such
  deletions, 4 + 1 + 1, 141 of its 159 gating findings, on three such pairs (126 `missing-row` and
  9 `carrier-missing` at chr22:23541792, 3 and 3 at 17341322 and 34811282; the review's world WDC,
  WD with one of its two lines in lower case: 16 gating). And the carrier gate reads a row's
  variants in upper case (`haplotype_prioritizer.parse_variants`) and looked them up in a store
  keyed as written, so an indel written in lower case (1,546 records of HPRC chr22) was in no
  panel: a population-level row naming one alone was labelled `single variant (not in the
  panel)` (17 alignment lines of the two-guide search above, run population-level), and one
  naming it beside other variants `not checked` (C2), though the panel holds it. Now the
  indel path names an indel or MNV record with REF and ALT in upper case
  (`indel_snp_cis.record_name`), and the three
  readers of the store key each line so (`build_indel_genotypes.store_key`): the two lines are one
  allele, their carriers merged as for two lines written alike (above), and the merge folds the
  two windows' rows into one, which gives the ID and the frequency of each line, as it does for two
  lines written alike. A row of an indel written in lower case names it in upper case, in
  `Variant_info_genome` and `Variant_info_spacer+PAM` (`chr22_36126395_CCCTCTCCCTCCTTTCCTTCCTT_C`,
  which 2.8.1 wrote `chr22_36126395_ccctctccctcctttcctTCCTT_c`). An SNV is read as written,
  as CRISPRitz's dictionary compares two listings of
  one; one whose ALT is written in lower case is not coded by CRISPRitz and is never searched: the
  input validator's full scan warns of it (above), and the real-data check refuses a run whose VCF
  holds one (none in the published panels' chr22). The review's world WDC (WD with one of its
  deletion's two lines in lower case) now gives WD's rows exactly, and the check 0 gating findings in
  both modes; the fresh two-guide per-sample HPRC search gives 18, none at a two-case pair: all are
  an MNV on a carrier's slot, which CRISPRme does not read beside a record (6 sites, 3 blocks each;
  METHODS section 8), as before (its population-level run: 6, the same class). Its population-level
  rows differ from before only in the case of their `Variant_info` and in `Other_motifs`, which
  counts the two windows' copies of one target once (430 to 428 at chr22:23541773). The README
  guide's HPRC searches (both modes, and on the unphased copy) name no indel written in lower case
  and are byte for byte as before. Covered by
  `PostProcess/test_indel_snp_cis.py`, `test_indel_genotypes.py`, `test_snp_row_indel_carriers.py`,
  `test_haplotype_prioritizer.py`, `test_indel_overlaid_rows.py` (the real `analisi_indels_NNN.py`
  on a deletion written on two lines in two cases) and `test_validate_inputs.py`.
- **An indel row whose first base is one an insertion adds is placed at the insertion**
  (every release with the fake-indel genome; both modes, any PAM).
  - The engine reports where the target's first aligned base (the leftmost, on either strand)
    is on the indel's haplotype. A base an insertion adds has no reference coordinate, and the
    row took the haplotype's: 1 to len(ALT) − len(REF) bases past the insertion's anchor. Its
    REF alignment, `CFD_score_REF` and REF seed counts were read there, from reference bases the
    target does not hold.
  - Such a row now starts at the insertion's anchor (its VCF POS; POS − 1 in the 0-based
    `Start_coordinate`), and its REF is read from there. Its ALT begins mid-insertion (on the
    − strand, shown reverse-complemented, it ends there), and its REF begins with the anchor
    base, which the target does not hold (on the − strand it ends with its complement); past the
    insertion the two are out of register by the number of inserted bases the target holds,
    less one (METHODS, "A row that starts inside an insertion").
  - README guide, HGDP + 1000 Genomes, NRG, 4 mismatches + 1 + 1 (2.8.1 code): genome-wide,
    per-sample, 229 highest-CFD rows (233 rows in any block) of 170,660 started so, a median of
    8 and up to 229 bases past the anchor, 95 of them on a tandem duplication; none at CFD ≥ 0.2
    (highest 0.155); and 151 of 468,715 alternative alignments. Population-level: 232 rows.
    chr22 alone, released 2.8.1 image against 2.10.0's post-analysis on it (CRISPRitz 2.9.1, index
    format 1): NRG, per-sample, 4 rows move to their anchor (17434488 (+) → 17434485, 21052642 (−) →
    21052639, 29300189 (+) → 29300175, 32435001 (−) → 32434998), their REF alignment and REF seed
    and non-seed counts with them, and at 32434998 `CRISPR_BULGE_score_REF` and
    `CRISPR_BULGE_risk_score`; population-level 3, 21052642 joining 21052638's site (above); NRR,
    per-sample, 11 (8 integrated rows and 3 alternative alignments, among them 44115947 (+), which
    now starts at 44115944 or 44115946, the anchors of two long insertions it lies in). Every other
    row is unchanged but for the columns the other fixes and the merge above change, a gene distance
    that moves with a start, and the SNP rows and carriers the SNP-row fix below takes off. The
    self-test (`complete-test --chrom chr22 --vcf_dataset 1000G`, then `validate-test`: 6 passed, 0
    failed, 0 skipped): in `cas9_sg1617` 2 of 632 rows move (21052642 → 21052639, whose SNP goes 24
    → 22, and 17434488 → 17434485); its other 630 rows and the other benchmarks' tables are
    unchanged.
  - Found from a 2.8.0 alternative alignment at chr22:22322768, the 7th base of
    G>GCAATCATTAT at 22322762. That row also paired the insertion's only carrier with four SNP
    alleles other individuals carry, which 2.8.1 already fixed: the row is gone.
- **An indel row's variants are numbered along its target in `Variant_info_spacer+PAM`**
  (every release; both modes, any PAM).
  - A SNP row numbers each variant by its distance from the row's start along the reference. An
    indel row did the same, but its target is bases of the haplotype carrying its indel, which
    holds the bases an insertion adds and not those a deletion removes: a SNP past an insertion
    was numbered len(ALT) − len(REF) bases short on the + strand (long on the − strand), one
    past a deletion as many long (short; the PAM's last base of a 23-base target read 24). So
    was a SNP on a row starting inside an insertion, which 2.8.1 placed at its first base's
    place on the haplotype (the SNP+indel oracle's Z3: 4 for 13), and one whose POS was before
    that place was numbered as a variant before the row's first base, 1 on the + strand and the
    last on the − strand. The indel's own number, at its POS base, was right.
  - Each variant is now numbered along the target on that haplotype, as a SNP row numbers it
    on the reference (METHODS, "`Variant_info_spacer+PAM`"); an insertion a row starts inside
    keeps 1 on the + strand, the last on the − strand. chr22:21052639 (−) above: its SNP
    G>A at 21052641 is numbered 22 (the target's 21st base: on the − strand an RNA
    bulge adds one, as on SNP rows), which 2.8.1 clamped to 24.
  - chr22, README guide, released 2.8.1 image against 2.10.0's post-analysis on it (CRISPRitz 2.9.1,
    index format 1): NRG per-sample 2 rows (21052639 (−) and 33013162 (−), a deletion and two SNPs
    past it: 17, 13 → 18, 14), population-level 1 (33013162), NRR per-sample 6 (also SNPs past a
    deletion at 12298466 (−), 32526400 (+), 25271097 (+) and 34192005 (+), the last one read 26 in a
    24-column target, now 21). Genome-wide, 2.8.1's tables of the NRG searches numbered as 2.10.0
    numbers them: README guide, per-sample, of the 5,951 rows whose highest-CFD alignment holds an
    indel, 239 change (165 with a SNP past a deletion, 57 past an insertion, 17 starting inside
    one), and 81 of 2,710 such alternative alignments; population-level, 231 of 5,923; TRAC guide,
    per-sample, 135 of 2,949 and 34 of 1,885. No SNP or MNV row changes.
- **`PAM_creation` of an indel row whose PAM holds bases an insertion adds, or a deletion's
  junction, is judged wherever in its repeat the VCF writes the indel** (every release; both
  modes, any PAM).
  - Up to 2.8.1 such a PAM was judged on the REF alignment's bases, read from the row's start.
    That is the reference PAM beside the protospacer only where the protospacer runs from the
    start to the PAM. On the − strand with a 3′ PAM (the + strand with a 5′ one) the start is
    the PAM's far end, and the bases read were len(ALT) − len(REF) bases from the PAM's place
    beside the protospacer: across a deletion's junction they held the deleted bases, the PAM
    of another protospacer (README guide, HGDP + 1000 Genomes: chr22:21056178 (−), CT>C at
    21056179, read CTC, GAG on −, where the reference holds TCT, AGA on −, beside its
    protospacer: NA, now AGG); past an insertion's anchor, the anchor and the bases after it
    (chr3:154367344 (−), C>CCT: AGG, where the reference holds GGG beside the protospacer, now
    NA; chr22:21052639 (−), T>TCACC: TGG, the reference holding CAG there, now NA). A row that
    starts inside an insertion was read from its first base's place on the haplotype, past the
    anchor, reference bases its target does not hold (with NRR, chr22:44115947 (+), wholly
    inside a 37- or a 33-base insertion: NA, now TAA at 44115944 and 44115946).
  - The PAM is now judged on the row's haplotype, the reference with the row's own variants,
    over every placement of its indel that keeps that haplotype: an insertion's bases rotated
    through its repeat, a deletion's junction moved through it (the VCF writes either at the
    left end of its repeat). A placement that leaves the PAM on reference bases: the bases it
    sits on there, which are its own (not created) but for a SNP alt of the row, read as the
    reference base it replaces. Else a placement that leaves the protospacer's base touching
    the PAM on a reference base: the reference PAM beside it. Else, an insertion that holds
    that base wherever it sits: created. A complex record, which replaces REF bases and adds
    bases or loses some (A>T and inserted bases after it: HPRC writes many insertions so), is
    its replaced bases, read as the REF bases they replace (as a SNP alt is), and its added or
    lost bases, placed as an insertion's or a deletion's. A row whose REF alignment is its ALT
    is never created; a PAM with a bulge, or holding an MNV's ALT bases, keeps the REF
    alignment's (for an MNV, the reference under it). A row that stands for two copies of one
    target inside one insertion is created when either copy's PAM is (the insertion added
    that copy, with a PAM it made). METHODS, "`PAM_creation`".
  - The answer does not depend on where in its repeat the VCF writes the indel. Read where it
    writes it, at the anchor where a row that starts inside an insertion now starts, a repeat
    unit inserted beside a reference PAM would put inserted bases in a PAM the reference
    holds; placed further right, the PAM is the reference's own: NA, as in 2.8.1 (with NRR,
    chr22:25415044 (−), the same alignment as its site's reference row; with the TRAC guide
    (CTCTCAGCTGGTACACGGCA), chr22:49196564 (−) and the row of 50193928 (−); genome-wide, README
    guide, per-sample, 48 alignments that start inside such an insertion). A tandem
    duplication of a PAM's bases after a reference protospacer is not created
    either (TRAC guide, the alternative alignment of chr22:20586277 (+), on the anchor of an
    insertion of two more copies of the 16-base unit the reference repeats after it, the same
    alignment as its site's reference row); a PAM wholly inserted after a protospacer of
    reference bases is created only where the reference lacks a PAM beside that protospacer
    (an inserted TGG before the reference's AGG: NA).
  - chr22, released 2.8.1 image against 2.10.0's post-analysis on it (CRISPRitz 2.9.1, index format
    1): README guide NRG per-sample 2 alignments, 21052639 (−) (TGG → NA; above, now one of the
    alternative alignments of 21052638's site) and the row of 21056178 (−) (NA → AGG);
    population-level the row of 21056178 (−); NRR per-sample 4, 21052639 (TGG → NA) and the three
    alignments of 44115947 (+), now the row of 44115944 and alternative alignments at 44115944 and
    44115946 (NA → TAA). TRAC guide, per-sample: the row of 28013783 (+), wholly inside an insertion
    (NA → CAG), and at 50193928 (−) the alternative alignments of three 48-base alleles that differ
    from the reference's unit there, so that no placement moves them off the PAM (NA → GAG). The two
    96-base alleles there hold the alignment twice, and 2.8.1 listed each copy (NA at 50193950, GAG
    at 50193998); each allele is now one alternative alignment at the anchor, GAG: one allele's
    copies are both created, and the other's left copy, which a placement further right puts on the
    reference's own PAM, is not, but its right copy is, and the row of two copies is created when
    either is. The row of 50193928 (−) and its alternative alignment of the 48-base allele that
    repeats the unit exactly keep NA, as 49196564 (−) does. In the SNP+indel oracle Z2 (wholly
    inserted) reads TGG, Z4 NA and Z5 TGG (partly inserted), T, T2, PR, PM, QN, QW, QZ and the
    Cas12a profile's C6 and CB1 NA, and QC, QD, QA, QX, QY, QV and CB3 created (below).
  - Genome-wide, 2.8.1's tables with each row's `PAM_creation` computed by 2.10.0's code
    (2.8.1's code gives back every value in them), highest-CFD alignments: README guide,
    per-sample, 97 go NA → created (70 rows and 27 alternative alignments): a PAM across a
    deletion's junction on the − strand 56 (CFD at most 0.083), a PAM holding inserted bases on
    a row that starts inside the insertion 27 and on a row that starts on a reference base 14;
    and 4 created → NA (3 rows and 1 alternative alignment); population-level 70 and 3 rows;
    TRAC guide, per-sample, 109 (67 and 42: deletions 29, CFD at most 0.039; rows inside an
    insertion 58; rows on a reference base 22) and 5 (2 and 3). HGDP + 1000 Genomes holds no
    complex record.
  - HPRC's chr22 holds 8,022 complex records with a net insertion and 1,055 with a net
    deletion (against 106,118 insertions and 80,464 deletions). Every block of a row that names
    one reads as an independent re-derivation from the VCF record reads it (the README guide,
    per-sample, 12 blocks of rows and 1,299 of alternative alignments, population-level 12;
    the TRAC guide, per-sample, 111 and 4,299). With the TRAC guide three CFD-0 alignments
    wholly inside the added bases of three complex records are not created (9 blocks): the row
    of chr22:37509610 (−) and alternative alignments at 42134819 (−) and 42135364 (−), whose
    added bases repeat the reference after them, so that placed further right their PAM is
    the reference's own GGG.
- The SNP+indel oracle (`test/indel_snp_oracle`) checks every indel row's start (strictly: at
  an occurrence on a listed carrier's haplotype that holds the row's variants; up to now a row
  with no carrier occurrence within 40 bp of its start was not checked),
  `Variant_info_spacer+PAM` and `PAM_creation` against its carriers' haplotypes, plants Z, Z2,
  Z3 (rows starting inside an insertion), Z4 and Z5 (a PAM partly inserted, the reference PAM
  beside the protospacer valid and not), X (an insertion holding two targets: both are listed),
  X2 (one holding a target twice: one row), T and T2 (an insertion after a reference
  protospacer that copies the reference's own PAM, or adds one before another reference PAM:
  not created), PR, PM, QC, QN, QD, QW and QA (PAMs judged wherever in its repeat the VCF
  writes the indel: tandem duplications of a unit holding a target, an insertion that holds
  the protospacer's end in both its places, a non-repeat insertion on a row starting on its
  anchor, a deletion across a − strand PAM, a wholly inserted PAM beside a reference
  protospacer holding a valid PAM, an insertion of the protospacer's end and the PAM), QX, QY,
  QZ and QV (complex records: A>T and C>T with a whole target in their added bases, G>CCA on a
  row starting on the replaced base, TG>C across a − strand PAM), and has a Cas12a profile (5′
  PAM TTTV, both strands, rows starting inside an insertion, PAMs wholly and partly inserted,
  a tandem duplication of a PAM (C6), a SNP past a deletion, and CB1 and CB3, the placement
  rule on each strand); its `PAM_creation` expectation builds every placement of the row's
  indel on the haplotype the row's variants make. The released v2.8.1 image fails the oracle
  231 times (SpCas9) and 60 times (Cas12a); 2.10.0 passes both. Unit tests:
  `test_indel_variant_positions` (the real integrator), `test_indel_pam_creation`,
  `test_indel_row_coordinates` (rows starting inside an insertion, rows whose PAM an insertion
  adds after a reference protospacer, PAMs judged wherever in its repeat the VCF writes the
  indel, complex records, 3′ and 5′ PAMs, both strands), `test_indel_snp_cis` (the placement
  rule against every placement built whole, on 3,000 random rows and on 3,000 more, over 1,000
  of them complex records, and the same answer for every record that writes one haplotype) and
  `test_representative_tiebreak` (the real merge: alignments that start inside one
  insertion, and two copies of one target, created when either copy is, in either order).
- **The release gate `scripts/check_alt_pam.py` lets an indel row's REF text hold an N where
  the reference genome does** (the gate's no-N check, since 2.8.1).
  - An indel row's REF text is the reference at the row's coordinates, as many bases as its
    target holds, and next to an assembly gap it can run into the gap while the target, on
    the haplotype carrying the indel, holds none. On HPRC's chr22 the README guide's rows at
    chr22:11978143 (−) and 11980225 (+) (3 blocks each) fail 2.8.1's gate, and 15 such
    targets of the released 2.8.1 image's table do.
  - Such an N now passes where the genome itself holds an N at that position, read from the
    row's start and strand and the text's gaps (the genome given with `--genome`, else the
    run's own: `Genomes/<Genome_ref>` beside its `Results` folder); an indel row is an
    ALT-origin block whose `Variant_info_genome` names an indel. Every other N fails as
    before: an ALT text's, a REF-origin block's REF text (the reported target itself), a SNP
    or MNV row's REF text, an N over a reference base, and one the gate cannot check against
    a genome. HPRC chr22, README guide, per-sample: those 6 targets pass and nothing else
    changes (on the released image's table, its 15); the chr22 searches of HGDP + 1000
    Genomes, the self-test and the SNP+indel oracle hold no N. Covered by
    `test_check_alt_pam`.
- **A carrier resolved PUTATIVE is listed only on the alignment it was resolved on** (2.8.1;
  indel and MNV rows, both modes).
  - A carrier whose genotypes cannot put an allele on the indel's haplotype (an unphased call,
    a missing allele, and every carrier of an MNV) gets the strongest target its own alleles
    allow on each alignment: the guide's base at a position the alignment pairs with a guide
    base, the base with the higher CFD at one it leaves as a DNA bulge (the reference base where
    no CFD is scored: a 5′ PAM, or a guide or PAM not 20 and 3 nt long). So two alignments of
    one target string can list it with different bases. The merge keyed an indel or MNV row by
    its cluster position, variants and upper-cased target, not by its aligned spacer: the two
    alignments' rows were one entry, the first in the file's order, with the carriers of both.
    Such a carrier could then be listed on an alignment where its own alleles give a stronger
    target, and the alignment it was resolved on left both tables. SNP+indel oracle, site U (S19
    0/1 for a deletion and for the SNP whose alt the guide needs): the target
    AC-TGA?ATCTGTAAGCAGGCTGG aligned as A-CTGAAATCTG… gives S19 the guide's A (CFD 0.93), and
    aligned as ACTGAA-ATCTG…, where that base is a DNA bulge, C (CFD 0.014; 0 with A); S19 was
    listed with C on the first, CFD 0.423, the row of S20, who has C (an alternative alignment,
    in the highest-CFD block; in the CRISPR-Bulge block S19 was listed so at 7999 too).
  - The merge now keys such a row by its aligned spacer too. Rows of one alignment still merge
    (one target from two VCF datasets, or two copies of one target inside an insertion), and
    the alignments of a site that name the same variants with the same carriers are still one
    entry, the best by the order that picks each block's representative (METHODS, "SNP+indel
    off-targets: each carrier's own haplotype"). A phased carrier has a target string on every
    alignment of it, so only alignments of one target string that list different carriers,
    through a PUTATIVE one, now keep an entry each. At site U, S19 is listed under the deletion
    alone only on alignments where the SNP's base is a DNA bulge, and A-CTGAAATCTG…, with C
    facing the guide's A, lists S20 alone. Carrier sets are compared whole, so that every block
    keeps one row count whatever its order; so an alignment whose carriers are all listed on a
    better alignment of that target and variants still keeps an entry when the better one also
    lists a PUTATIVE carrier resolved to another base on the weaker one: with a SNP whose base
    the weaker alignment leaves as a DNA bulge, only the carriers homozygous for the SNP, or
    phased with it in cis, are listed on the weaker one with the SNP. That entry is an extra
    alternative alignment, its carriers and score right (none in the searches below, nor on
    chr22 of HPRC and its unphased copy).
  - Cost: the merge sorts and folds one row per alignment where it sorted one per target
    string, and still reads each target's place on the haplotype once, whatever its alignments.
    With the README guide on HPRC's chr22, per-sample, each of its three merges takes 4.4 to
    4.6 s against 2.8.1's 3.1 to 3.5 s (the cost note of the first entry above).
  - chr22 of HGDP + 1000 Genomes, the searches above (README guide NRG per-sample and
    population-level, NRR per-sample, TRAC guide per-sample): every table is byte-identical
    with and without this fix. Of their groups of indel or MNV rows with one cluster position,
    variant list and target string, 105 of 233 (NRG, both modes), 259 of 629 (NRR) and 131 of
    300 (TRAC) hold two or more alignments, each group with one carrier set, and each still
    kept as before. On chr22 of HPRC and of a copy of it with every fourth sample unphased (the
    real-data indel check's builds; README guide per-sample and population-level, TRAC guide
    per-sample), every table is byte-identical too. The real-data indel check
    (test/282-indel-realdata-check), which found the defect, found no carrier listed so on chr22
    of HGDP + 1000 Genomes, HPRC or an unphased copy of HPRC.
  - Covered by `test_representative_tiebreak` (site U's rows through the real merge; without
    this fix S19 is on the CFD-0.423 alignment, or that alignment is gone; and an alignment
    whose carriers are a part of a better one's keeping its entry in every order) and by the
    SNP+indel oracle. It checks again every carrier an indel row lists (a phased one at an
    occurrence of the target holding every variant the row names, an unphased one carrying every
    allele it names) and that every variant it names is a record of the VCF, which the stricter
    start check above had stopped doing, and now also that an unphased carrier's target is the
    strongest its own alleles allow on the row's alignment. Without this fix that check fails
    twice, at site U.

## [2.8.1] - 2026-11-05

### Added
- **`--prune-pam`: one index can serve many PAMs.** A *PAM-less* index stores the raw genomic flank
  instead of a matched PAM, and `--prune-pam` applies the PAM you pass during the search, discarding
  a leaf whose stored flank cannot match it. One index therefore answers any 3′ PAM no longer than
  its flank, instead of one index per PAM. Verified genome-wide over HGDP + 1000 Genomes (4,091
  samples, NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per-sample): a 3-nt-flank PAM-less index, an
  8-nt-flank PAM-less index, and a dedicated NRG index produce **byte-identical**
  reports: identical integrated results, 182,369 rows each (the dedicated arm ran on a fresh install
  without the COSMIC licence, and matches once that empty column is set aside), and for NRR 386,012
  rows each, over a dedicated index, both PAM-less indexes, and the 3-nt one searched without
  `--prune-pam`. That was measured for one guide, which has no site at a contig's end: from a flank
  wider than the PAM the build does not index a PAM at a contig's very end or start (on hg38, 549
  NRG and 246 NGG sites, none on a primary chromosome; see Changed). Pruning changes only which
  leaves are written out, never how a guide is aligned, so it cannot move an alignment or a score.
  The cost is size, since nothing is excluded at build time: 110 GB against 28 GB for the HGDP +
  1000 Genomes variant index. Requires `--index-path`, and is command-line only — the web interface
  never passes the flag. The published PAM-less indexes store the flank on the 3′ side of the guide
  and serve no 5′ PAM, so a Cas12a (`TTTV`) search needs a dedicated index. A shorter PAM off a
  wider flank is served for a 3′ PAM only; a 5′ PAM needs an index built from a 5′ all-`N` PAM file
  at exactly its width, which its `.pam_build` records (an all-`N` index without one is taken for a
  3′ build at the published widths, 3 and 8: see Fixed), and the engine refuses any other width
  pairing rather than mis-parsing the index. The `NNN_3` indexes are published alongside the
  dedicated ones; a flank wider than 3 nt is not readable by earlier releases, so those live in a
  separate repository. Needs a CRISPRitz build that supports `--prune-pam`. See [`METHODS.md`
  §9](METHODS.md#9-pam-less-indexes-and-dynamic-pam-filtering).
- **A second dataset repo for wide-flank PAM-less indexes.** A PAM-less index with a 3-nt flank
  stores a 3-nt PAM region, as a dedicated 3-nt PAM index does, so any release can read it and
  `NNN_3` stays in `lucapinello/crisprme-data` (though a 2.8.0 search of the rebuilt `NNN_3`
  archives can fail or report malformed rows: see Fixed). A wider flank cannot be read by a
  release predating `--prune-pam`, so those are published to
  `lucapinello/crisprme-data-general-index`, which older defaults never point at. `download` and
  `publish-index` route on the flank width parsed from the index name — not on a
  list of names, so it cannot be forgotten at publish time — `download` falls back to the default
  repo when the routed one lacks the archive (and, if both lack it, says it looked in both), and the
  catalog listing merges both repos, taking each name's row from the copy `download` would fetch, so
  what is offered is what can be fetched. If the default repo cannot be read the catalog is empty,
  as before, rather than offering only the slower wide-flank indexes. Dedicated-PAM and `NNN_3`
  indexes are listed and downloaded from the default repo exactly as in 2.8.0, never looked up in
  the new one. `--hf-repo` / `CRISPRME_HF_REPO` still override. `download` takes
  an index's genotype companion and reference index from that index's own repo, so the new repo
  carries its own copy of `genotypes_hg38_HGDP1kGP.tar.gz`, byte-identical to the main repo's.
- **`indexes/SHA256SUMS` in both dataset repos**, written by `scripts/hf_sha256sums.py` from the
  checksums the Hub already records (a Git LFS object is addressed by its SHA-256; the oids match
  the independently computed `build-inputs/*.SHA256SUMS`), so, inside `indexes/`,
  `sha256sum -c --ignore-missing SHA256SUMS` verifies a download without anyone re-hashing 400 GB
  (`--ignore-missing` skips the archives you did not download; without it, a one-archive download
  reports every other archive as `FAILED open or read` and exits 1). `--staged` folds in the
  checksums of archives about to be uploaded; regenerating after the upload and diffing proves the
  Hub holds exactly what was staged. It refuses to write an empty file, which would verify nothing
  and read as a pass.

- **Unrecognized `complete-search` flags are now rejected.** `complete-search` parses `sys.argv` by
  hand, so a flag nothing reads was simply ignored: `--prune-pams` instead of `--prune-pam` ran a
  genome-wide search with pruning quietly off, finished without complaint, and was detectable only by
  diffing the results against a dedicated index. Unknown flags now fail up front, with a suggestion
  for a near miss.
- **Settings: configurable number of concurrent search jobs.** The web app ran at most 2 searches at
  once (a hard-coded legacy limit); this is now adjustable from the Settings page (or the
  `CRISPRME_MAX_CONCURRENT_JOBS` environment variable), clamped to `[1, 16]`. The **default stays a
  conservative 2** — each genome-wide search is memory-heavy (~tens of GB of RAM) — so a laptop / the
  public server is never accidentally oversubscribed; raise it only on a machine with enough RAM and
  cores. Extra jobs beyond the limit simply queue. Applies on the next app restart. When
  `CRISPRME_MAX_CONCURRENT_JOBS` is set it takes precedence: the Settings card says so and its input
  is disabled, and a value that is not a whole number from 1 to 16 is ignored or clamped, with a
  one-line warning at startup.
- **Settings: Hugging Face download progress bar.** Index/annotation downloads now show percent,
  size, ETA, and phase in the Settings page, backed by an env-gated progress file
  (`CRISPRME_PROGRESS_FILE`); the progress path is a no-op when the variable is unset, so the CLI
  happy path is byte-for-byte unchanged.

### Changed
- **CRISPRitz pin moved to v2.9.1**, the tag carrying `--prune-pam` plus eight engine fixes: a
  genomic `N` in a stored PAM/flank region no longer drops the target; a non-`.bin` sidecar in the
  index directory no longer crashes the search; a PAM/flank width mismatch is refused up front
  instead of mis-parsing the leaf stream, and `buildTST`/`searchTST`/`searchBruteForce` now have
  their exit status checked, so a crashed search can no longer be reported as an empty result; a
  multiallelic site no longer compounds its indel alleles into one fake contig; the compiled C++
  `enricher` no longer drops a SNP wherever an indel shares its position; `enricher.py` keeps the
  real SNP at a multiallelic site whose ALTs include `*` (it left the base unchanged; the C++ port
  already coded it); and `add-variants` now fails when an enricher fails, instead of exiting 0 with
  that chromosome's enriched genome missing, which the index build would then have indexed as it
  stood (the Python enricher, moreover, no longer needs a `python` on PATH, and no longer races
  its sibling workers to create their shared output folder, which could lose a chromosome); and
  `buildTST` skips a contig too short to hold one target, such as the header-only FASTA the
  enricher writes for a chromosome with no indel, instead of crashing on it.
  Applies to the Docker build (`crispritz_ref`) and `install_from_source.sh`
  (`CRISPRITZ_REF`). The `*` case needs a multiallelic VCF, and every CRISPRme-supplied VCF is
  split, so it moves no published index.
  **No published index is affected** by the C++ enricher fix: the image built only the three search
  binaries, so enrichment ran through `enricher.py`, confirmed by direct measurement on the
  artifacts (chr22:10684299 carries a SNP and an indel as separate records; every enriched
  `chr22.enriched.fa` on the build host holds the IUPAC `M` there, not a plain `C`). A **source**
  install is the exception -- `install_from_source.sh` has installed the compiled binary for some
  time, so an enriched genome produced that way before CRISPRitz 2.9.1 is missing the SNP at every
  position that also carries an indel, and should be re-enriched.
- **The compiled `enricher` is now installed in the image, making it the default enrichment path.**
  It is built beside the three search binaries and copied to both places `crispritz.py`'s
  `_enricher_command()` checks, matching what `install_from_source.sh` already did and removing a
  silent divergence between the two install routes. Measured **1.82x** faster than `enricher.py`
  through the real `add-variants` CLI on the FULL chr22 of HGDP + 1000 Genomes -- 1,998,332 records,
  4,091 samples, 2.7 GB of output: 1,776 s vs 974 s -- on the stage that dominates a from-scratch
  genome-wide build. Gated on CRISPRitz 2.9.1's byte-identity fixes, verified at that same scale:
  all four enricher outputs match `enricher.py` hash for hash, as do the 17 synthetic cases of
  `verify_enricher.sh`, the SNP+indel combinations among them. The same holds off Linux/amd64:
  natively on macOS arm64 (Apple clang, libc++) the compiled port reproduces those full-chr22
  hashes exactly, and so does `enricher.py` under pandas 3 (the recipe leaves pandas unpinned);
  inside the arm64 image, the binary it ships reproduces them too.
- **When a dedicated index for the PAM is installed beside a PAM-less one, the dedicated index is
  now searched.** The index scans settle every choice by rule: the requested PAM first, then an
  index built for this very guide length, then the smallest sufficient bulge budget, then the
  narrowest PAM-less flank. 2.8.0 kept the smallest budget first, so with `NRG_4` and `NNN_3`
  installed it searched `NNN_3`; 2.8.1 searches `NRG_4`,
  which gives the same report faster (genome-wide, the search stage took 44-51 min on a dedicated
  NRG index against 2 h 59 min on `NNN_3` under `--prune-pam`). No rule depends on glob order, so
  the choice cannot depend on the shell's locale either, which matters now that PAM-less indexes
  of two widths can sit side by side: `NNNNNNNN_3_hg38` globs before `NNN_3_hg38` under
  `LC_ALL=C` and after it under `en_US.UTF-8`. `PostProcess/test_index_resolution.py` checks each
  rule with the wrong candidate globbing first, and runs the whole shipped script against a stub
  engine to check that each scan hands the search the index it picked.

- **`build-index-only` makes a missing `--samplesID` impossible to miss.** Building a genotyped
  variant index without `--samplesID` (or with no usable entries) skips the dict-less Tier-0
  registry, so the index lacks fast post-analysis, the report's variant-count panel and the web's
  data-type labels, and cannot be published. The one-line notice is now a red banner on stderr
  that names the consequence and the fix (re-run with `--samplesID`). Behaviour is unchanged.

- **Install docs: Bioconda and macOS.** The README and website said Bioconda installs the old 2.1.x
  line; it has served CRISPRme+ (Python 3.11, Linux) since 2.7.9. The channel setup now uses
  `conda config` (Mamba 2 rejects `mamba config --add`), and the install names its version
  (`crisprme=2.8.1 python=3.11`), since an unpinned solve can silently fall back to a years-old
  release; `scorer-env create` and `doctor` are part of it. macOS is Docker-only: Bioconda has no
  macOS CRISPRitz, and the source build needs GCC (`buildTST` uses GCC's parallel-mode STL).

- **Docs and help text now match what ships and what is on Hugging Face.** The raw source VCFs are
  hosted (`vcfs/hg38_HGDP1kGP`, `vcfs/hg38_HPRC`, `vcfs/hg38_mega`, fetched with
  `download --what vcf --dataset <name>`), although the README's `download` reference, six pages
  under `docs/` and the error `--what vcf` gives for a missing dataset said they were not.
  METHODS, the HGDP + 1000 Genomes methods page and the troubleshooting page described
  population-level as the default analysis; they now describe the context-aware default the code
  has applied since 2.8.0 (per-sample on a phased genotyped index such as the default one,
  population-level otherwise), so on the default index it is `--population-level`, not leaving out
  `--per-sample`, that gives the lighter screen. The index tables now name the NGG, NRR and
  8-nt-flank PAM-less archives with their download sizes, all in GB (10⁹ bytes), and
  `docs/INDEX_PROVENANCE.md` records both dataset repos as they stand after the 2026-10-30 upload
  (the four replacements, `indexes/SHA256SUMS`). Help and warning texts that named unhosted
  datasets or stale values are corrected too: `build-index-only --help` now uses `hg38_HGDP1kGP`,
  `assembly-search --help` gives complete-search's `--max-total-edits` default as 6 (it said 4),
  and the warning for the retired `--fast` no longer says it used to force per-sample.
  The general-index card's search example now runs: it lacked `--index-path`, which `--prune-pam`
  requires, and `--samplesID`, and named a `--vcf` file `download` does not write. The website and
  the index pages say that a PAM-less index serves a shorter PAM only when it is a 3′ PAM (the
  card listed Cas12a's `TTTV`), and the README, the index pages, METHODS and the
  `complete-search --help` entry for `--prune-pam` that the published ones serve no 5′ PAM at any
  width (they said an equal-width one worked). README §1.1.2 checks an install with
  `crisprme.py --version` alone, and §3.1 shows the help with `crisprme.py --help`, since a bare
  `crisprme.py` exits 1 and creates the CRISPRme folders where it runs.
  `RESUME.md`, a scratch file committed by accident in 2.5.4, is gone from the tree, so the Docker
  image no longer ships it.

- **The web Settings page gives sizes in decimal GB, and `complete-search` prints its
  analysis-mode banner only for a variant search.** Settings divided by 1024 while printing GB, so
  every size it shows in GB now reads about 7% larger, as the docs and Hugging Face give it (GB =
  10⁹ bytes): `NNN_3_hg38` goes from 66.6 to 71.5 GB, and the 50 GiB upload limit now reads
  53.7 GB. A reference-only search printed the population-level banner, although it runs no
  variant analysis; it now prints none. The per-sample banner says whether `--per-sample` or the
  phased index's default chose it, and the population-level banner points at `METHODS.md` instead
  of a design note removed in 2.6.0. Search results are unchanged.

- **`--per-sample` on an index without the dict-less tiers now says what it reports.** An index
  built without `--samplesID` has only the classic per-sample dictionaries, and no
  `Dictionaries/registry_<vcf>/` or `genotypes_<vcf>/`. This covers the index that benchmark-mode
  `complete-test` builds, and a build from the web Settings page before this release (which now
  passes `--samplesID`: the entry on a search of a VCF dataset built there). On such an index,
  `complete-search --per-sample` still promised "CONFIRMED cis phasing, named carrier samples and
  exact joint AF", but wrote no phase-confirmation or SNP+SNP co-occurrence table. The banner now
  says that this dataset's results name the carrier samples, with the `AF` of the classic per-sample
  dictionaries (an indel's, of the index's indel log), which CRISPRitz `add-variants` copied from
  the VCF, but have no CONFIRMED/PUTATIVE cis phasing, no SNP+SNP co-occurrence and no exact joint
  AF. It also says to rebuild with `build-index-only --samplesID` or download a prebuilt index. The
  banner of an index with the tiers, and every result, are unchanged. Covered by
  `PostProcess/test_build_index_only_atomic.py`.
- **`complete-test` now says when it ignores `--vcf_dataset`.** The benchmark mode, which
  `--vcf_dataset` configures and `validate-test` checks, runs only with
  `CRISPRME_ALLOW_ONDEMAND_BUILD=1` set. A Docker user who put the variables in front of `docker
  run` set them on the host only, and `complete-test` silently ran the prebuilt-index smoke
  instead. It now prints one line on stdout saying so, and that in Docker the variable goes with
  `-e`. The run and its exit status are unchanged. Covered by `PostProcess/test_complete_test.py`.

- **A search of a VCF dataset built in the web Settings page reads its dict-less tiers, so its
  tables change.** A variant build in Settings now passes the dataset's sample list as `--samplesID`
  when its VCFs have samples (see Fixed), so it also writes the dict-less tiers
  (`Dictionaries/registry_<dataset>/` and `Dictionaries/genotypes_<dataset>/`). Searches of an index
  built (or built again) in Settings now read them, as searches of a command-line `--samplesID`
  build and of a downloaded index do, where they read the dataset's dictionaries before. On a phased
  dataset a search now runs per-sample by default, as on the default index (*Population-level* on
  the search form is the faster worst-possible screen), and a per-sample search names on each row
  only the samples one of whose haplotypes carries all of the row's variants (of two samples with
  the same two SNPs, one that carries both on one haplotype is named on the row of both SNPs, and
  one that carries one on each haplotype is not, where a population-level search names both).
  Measured on a phased 6-sample panel, with a population-level search of one index on its
  dictionaries and then on its tiers, the tables change as follows: a SNP's allele frequency is
  computed from the genotypes of the samples the list names, to six significant digits
  (`0.333333`), where the dictionaries gave the `AF` that `add-variants` copied from the VCF
  (`0.3333`), and an indel's is still the `AF` of the index's indel log; and `Other_motifs` counts
  the alignments the tiers give at a site: one fewer at a site where the dictionaries also listed
  the reference alignment as a variant one, for the carriers of a SNP that destroys its PAM. The row
  of an indel and a SNP that one haplotype carries is the same on both, with the same carriers,
  since this release's indel fix resolves it from the per-sample SNP dictionary too (the entry "An
  indel off-target near a SNP was reported with SNP alleles its carriers do not have" under Fixed).
  The result folder also holds the SNP + indel and SNP + SNP co-occurrence and phase-confirmation
  tables, and the shareable report gives the dataset's SNP and indel counts and its SNP + indel cis
  co-occurrences. A dataset built in Settings before this release has no tiers until its variant
  index is built again there, and its searches read its dictionaries until then; their SNP + indel
  rows change all the same, with that indel fix. Writing the tiers takes memory that a build without
  them did not need: at panel scale ~150-180 GB a chromosome, by default up to 4 chromosomes at once
  (~700 GB, `CRISPRME_TIER_WORKERS` sets how many). A build that runs out of it fails with a message
  that says what to do (see Fixed).
- **A genome uploaded in the web Settings page under the name of an installed genome replaces it,
  and is refused while indexes built on a genome of that name are installed.** The upload went into
  that genome's folder, beside the files already there (and so into the shared copy, when that
  folder is a link to one). It is now built in a hidden folder and replaces the installed genome
  only once whole; a link is replaced and the copy it points to left alone. It is refused, at its
  first chunk, while indexes or variant genomes built on a genome of that name are installed
  (`genome_library/<PAM>_<N>_<genome>[+<dataset>]`, `Genomes/<genome>+<dataset>`): every search
  and build of the new genome would go on using them, with the sequences and positions of the one
  they were built from. Remove them first (an index under *Remove installed data*, a
  `Genomes/<genome>+<dataset>` folder from your data folder), or upload the genome under another
  name; to add the very genome an installed index was built on, use *Download genome*, which this
  does not refuse.

- **The documented self-test now runs every benchmark, and it ends `6 passed, 0 failed, 0 skipped`.**
  **Affected:** anyone who runs `complete-test`'s benchmark mode and `validate-test`, as README
  §1.2, §1.3, §2.2.2, §2.2.3 and §3.2.1, `docs/DOCKER_QUICKSTART.md`, both commands' `--help`, the
  installer's last line and `validate-benchmarks.yml` give them. All of these set
  `CRISPRME_SKIP_HEAVY=1`, on the claim that the bulge-2 index builds need more than 16 GB of RAM,
  so each left out 2 of its 4 benchmarks. That claim was wrong: the run peaks at 3.66 GB of RAM,
  and it passes in a container held to 16 GiB, 4 CPUs and no swap. **What changes:** the commands
  no longer set it, and CI no longer does, so every benchmark runs, and CI fails if one did not.
  All six took about 25 minutes at 8 threads on an Apple Silicon Mac in testing, where the two
  light ones took about 12. **What to do:** nothing. The variable stays as an opt-out for a very
  small machine: `complete-test` says what it leaves out, and `validate-test` ends `2 passed, 0
  failed, 4 skipped`. A folder that ran the old command has no output for the others:
  `validate-test` there ends the same way, saying `complete-test targets not found`. Covered by
  `PostProcess/test_validate_test.py`, `PostProcess/test_complete_test.py` and
  `PostProcess/test_ci_workflows.py`.
- **Two new benchmarks that one lost site fails, `cas9_sg1617_n5` and `cas12a_hbg_n5`.** The two
  bulge-2 benchmarks could hardly catch a lost site: `cas12a_hbg_simple`'s brute-force reference is
  a header with no rows, so that case compares empty with empty, and `cas9_sg1617_simple`'s holds 3
  rows at 2 sites. The same guides at `--max-total-edits 5`, on the same indexes, have 381 and 92
  sites on chr22, and CRISPRme matches brute force on them, on CRISPRitz 2.9.1 and on the 2.9.0 CI
  overlays. No result was lost: this was a gap in the test's power. Their references are in
  `test/benchmark/max-edits-5/`, which Docker, source and Bioconda installs all hold, so they need no
  network. Covered by `PostProcess/test_validate_test.py`.
- **CI runs what it says it runs** (contributors). `unit-tests.yml`'s job limit goes from 10 to 20
  minutes: the hosted job took 523 s at the first release head of 2.8.1, and the steps have grown by
  a quarter since, which predicted about 640 s, past the limit (GitHub has since run it in 14 min
  50 s and 13 min 55 s on two release candidates). `web-e2e.yml` gets a step that runs, in the image, the tests that skip in the light
  job for want of BEDOPS, minimap2, impg or Dash: `test_assembly_annotate`, which no workflow named,
  the impg truth tests, the flag-gate behaviour tests, the `post_process.sh` test, whose own step in
  the light job tested nothing, `test_index_resolution`, and the web app's start-up and shell-safety
  tests (`test_app_settings`, `test_web_shell_safety`): 527 tests, about 100 s. It fails on any
  skip but the two the image cannot help. `validate-benchmarks.yml`'s smoke step runs the real
  Hugging Face download test, which no job ran. The path filters of `unit-tests` and `web-e2e` now
  start, on a change to any file one of their tests reads, a job that runs that test, as found by
  replaying each step with every file read logged: `unit-tests` lists every file its tests read,
  and `web-e2e` every file read by the tests only it runs (the docs and scripts its in-image step
  reads through `test_complete_search_flags`'s scans start `unit-tests`, which runs the same
  scans). Before, a pull request that changed only `pages/`, `scripts/`, `docs/`, the README,
  `test/benchmark/`, `test/data/`, `assets/`, or a script outside `PostProcess/` such as one in
  `seq_script/`, did not start `unit-tests`, and one that changed only `PostProcess/`,
  `seq_script/`, `PAMs/` or `download_hprc_assembly.py` did not start `web-e2e`. Covered by
  `PostProcess/test_ci_workflows.py`, which `unit-tests.yml` runs.
- **`docs/RELEASE_CHECKLIST.md` Gate B gives the local gate** (maintainers): the exact commands that
  leave nothing skipped, on the release image, amd64 and arm64, with the three skips it allows
  named. The six function-based test modules run as scripts there, since `python -m unittest` alone
  reports `Ran 0 tests` and `OK` for them.

- **`--prune-pam` from a flank wider than the PAM: what it cannot index is said, in the docs and
  by the search.** The README, METHODS section 9, docs/PRECOMPUTED_INDEXES.md,
  docs/DYNAMIC_PAM_INDEX_COMPAT.md, the general-index card in `docs/hf/` and `complete-search
  --help` said a PAM-less index with `--prune-pam` gives a report identical to a dedicated index.
  That was measured genome-wide for one guide, which has no site at a contig's end. The index
  build scans for a PAM only where a whole flank of K bases fits, so with K > L a + strand PAM with
  1 to K − L bases after it to the end of a contig, and a − strand PAM in the contig's first K − L
  bases, are not indexed, and a − strand PAM exactly K − L bases from the start is reported on
  the + strand at another position; K = L is unaffected. **Affected:** searches of the K = 8 index
  for a 3-nt PAM (on hg38, 145 + 101 NGG and 308 + 241 NRG sites, about 7.5 × 10⁻⁷ of all, none on
  a primary chromosome for those PAMs), new with K = 8 in 2.8.1. **What to do:** nothing for a
  guide whose sites lie on primary chromosomes; for one near the end of an alt, random or unplaced
  contig, use an index of the PAM's width or a dedicated one. The docs now say so with these
  counts, and a search that reads a flank wider than its PAM prints one `NOTE` on stdout, also
  written to the run's `Params.txt`. The engine is unchanged.
- **The limits every index has at contig ends and `N` runs are documented.** METHODS section 9
  and docs/PRECOMPUTED_INDEXES.md now say that an index of budget N cannot report a site whose
  guide lies within N − 1 nt of a contig start, a contig end or an `N` run, which a search
  without bulges reports from the FASTA (a variant search without bulges of a download reads the
  variant index, so its SNP variant rows have the limit too), so that an index of a larger budget
  is not a strict superset of a smaller one (for `NRG_3_hg38`, 1,050 NRG sites, 535 on primary
  contigs; the exposure that matters is a guide whose own on-target sits at such a boundary);
  that a − strand PAM at a contig's first base is reported on the + strand (11 NGG and 37 NRG
  hg38 contigs); and that a + strand PAM on a contig's last bases is found by no search (8 NGG
  and 31 NRG contigs), these last two only on alt, random or unplaced contigs. These are engine
  limits, in 2.8.0 too; the engine and the published indexes are unchanged.
- **A known limit, stated with its scale: at a split multi-allelic SNP site only the last record's
  ALT is searched.** At a position with several biallelic SNP records only the last one's ALT
  reaches the enriched genome and the stores, so an off-target that needs an earlier one is not
  found, or is reported at a lower score, without its carriers (a VCF
  line that lists several ALTs itself is coded with all of them). On the HGDP + 1000 Genomes
  panel, 2,486,165 of its 127,666,429 SNP records (1.95 %), at 2,441,139 positions, are not coded
  for this reason. Per guide (NRG, 4 mismatches + 1 DNA + 1 RNA bulge, per-sample, genome-wide)
  this leaves 3 to 7 sites whose true CFD puts them in the top-1000 out of it or understated
  there, and at most one site reaching CFD >= 0.5 in the five guides tested (the README's, ACTGA,
  TRAC, EMX1 and VEGFA site 3; the one is EMX1's chr17:12179615, 11 haplotypes, CFD 0.508), and 52
  to 324 rows per guide (99 % of them at CFD < 0.1) name an individual who carries such an ALT,
  read as the reference base, for a sequence they do not have. docs/INDEX_PROVENANCE.md and
  METHODS section 2 give what was measured and how: on chr22, 45,837 of 1,757,775 SNP records
  (2.6 %) are not the last at their position, 5,104 with AC ≥ 41 and 1,566 with AC ≥ 818, counted
  in the normalized source VCF; and, read in the shipped enriched genome, registry, genotype store
  and dictionary at six such positions, the earlier ALT is absent from each, up to AC = 8,179
  (`rs807755` T>A at chr22:19142014, whose position the stores hold as T>G at AC = 1). This is not
  new: 2.8.0 and every earlier release behave the same (the SNP dictionaries have always kept only
  the last record of a position, and the enriched genome since CRISPRitz 2.6.0, 2021), and 2.8.1
  does not change it. It will be fixed, with rebuilt indexes, in the next release.

- **Release docs (maintainers): the CRISPRitz pin of the Bioconda recipe is bounded, and the
  version sweep compares like with like.** `docs/RELEASING.md` and the `release-crisprme` skill
  ask for `crispritz >=X,<Y`, below the next CRISPRitz that changes the index format or the
  enricher: `>=2.9.1,<2.9.2` for 2.8.1, which reads index format 1, where the CRISPRitz planned as
  2.9.2 writes format 2 and enriches multi-allelic sites another way. With a bare `>=` the solver
  would install that engine beside 2.8.1, which cannot read the indexes it builds:
  `build-index-only` fails, and a bulge search refuses such an index when it is reused. The
  checklist's version sweep asked the Dockerfile's `crispritz_ref`
  (v2.9.1) to read CRISPRme's version (2.8.1), which no correct release could pass; it now asks it
  to equal `install_from_source.sh`'s `CRISPRITZ_REF` default, which
  `PostProcess/test_install_from_source.py` checks. The release docs also said the Dockerfile
  carries no CRISPRme version, while its header comment names the release, and
  `scripts/prepare_release.py` bumps `crisprme.py` alone: a release that followed them would have
  left the comment at the previous version, as it read 2.4.0 at every tag from v2.5.0 to v2.7.8.
  The sweep now asks for it, `docs/RELEASING.md` and the skill bump it by hand with a command the
  same test runs, and the test fails while the comment names another release than `crisprme.py`.
- **`docs/RELEASE_CHECKLIST.md` Gate B runs the engine checks on the release's own CRISPRitz**
  (maintainers). CI runs `test/index_flex/check_index_flex.py`, the SNP + indel oracle
  (`test/indel_snp_oracle/run_oracle.sh`) and `scripts/check_alt_pam.py` on search results only
  in `validate-benchmarks.yml`, with the branch overlaid on the published image: for 2.8.1 that is
  2.8.0's, with CRISPRitz 2.9.0 where 2.8.1 ships 2.9.1. There index_flex passes both K = 8 cells
  by refusing them, as on any engine without `--prune-pam`, and still prints `INDEX FLEX PASS`.
  The local gate's step 4b runs the three on the release image, amd64 and arm64, and fails unless
  both K = 8 cells read "same results as" (a dedicated index, and NRG), the oracle prints
  `ORACLE PASS` and `check_alt_pam.py` passes step 4's tables. Trap T-42 in the checklist.
  `PostProcess/test_ci_workflows.py` pins the step and its greps to what the checks print.

- **In a population-level report an indel (or MNV) row's Observed cell reads `N carrier(s)`**, not
  `N carriers (>=1 of K vars, poss. chimera)` as a multi-variant SNP row does, now that such a row
  can name the SNP beside its indel: its carriers are resolved one by one, CONFIRMED where each is
  placed on its phased haplotype and PUTATIVE where the phase is unknown (an unphased call, an MNV).
  On an install without SNP genotypes (`download --what index --no-genotypes`) a SNP+indel row still
  lists every carrier of the indel, with or without the SNP (PUTATIVE in `indel_snp_cooc.tsv`), as
  in 2.8.0. The README and METHODS section 7 say so.

- **A known limit, stated with its scale: a SNP row lists a carrier whose haplotype holding its SNPs
  also carries an indel inside the target.** A SNP row's carriers are resolved from SNP genotypes
  only, so an individual whose haplotype holding the row's SNPs also carries an indel inside the
  target is still listed, though its own sequence there is the indel's (reported on the indel row
  when it is in budget). On chr22 of the HGDP + 1000 Genomes panel (per sample, NRG, 4 mismatches +
  1 DNA + 1 RNA bulge, five guides) this affects 2.2% of SNP row-blocks (990 of 45,921, 194 in the
  integrated results), and 545 of them have no listed carrier that holds the sequence; with no
  bulges, 12 row-blocks. It only over-lists, no site is lost, and 2.8.0 behaves the same. METHODS
  sections 4 and 8 now say so; it is fixed in the next release.

- **A known limit, stated with its scale: on an index without the dict-less tiers, a
  `--per-sample` search can still report a dense window by its greedy worst-case alignment, whose
  listed carriers need not carry every variant in it.** `--per-sample` reports an off-target only
  when a real individual carries it, never a worst case put together from alleles that different
  individuals carry (2.7.10), on the observed-haplotype path, which reads the dict-less tiers. An
  index built without `--samplesID` (such as the one the benchmark `complete-test` builds, so the
  one the README's self-test leaves) has only the classic per-sample dictionaries, and a
  per-sample search of it takes the classic path, which still reports a window of more than
  `CRISPRME_IUPAC_CAP` variants by its greedy worst-case alignment, flagged in
  `High_variant_density_region`, with the carriers of any allele it takes. On the NGG index the
  self-test builds (`NGG_2_hg38_chr22+hg38_1000G`, chr22 of the benchmark's phased 1000 Genomes
  VCF), the README's guide with the `20bp-NGG-SpCas9` PAM, 4 mismatches, 1 DNA and 1 RNA bulge, at
  the default `CRISPRME_IUPAC_CAP`, gives 5 such rows (chr22:22410074, 22687072, 22687116, 22698190
  and 22888934, all on the - strand, 8 to 13 variants each, 4 mismatches and 2 bulges, CFD at most
  0.001). They name 102 individuals as carriers 114 times, and each time the carrier lacks one of
  its row's variants: 342 of the run's 56,730 SNP-row carrier listings (a row lists its carriers
  once for each of its three alignments). After `build-index-only --samplesID` the same search
  returns exactly those 5 rows fewer, and every carrier it lists has all of its row's variants on
  one phased haplotype. The `complete-search` banner and the report's Search mode row say so, and
  to rebuild the index with `build-index-only --samplesID` or use a prebuilt one, as the README and
  METHODS section 5 do; the report's legend (its Observed entry) and the caption of its top-1000
  tables, the `--per-sample` help, METHODS section 7 and `docs/HGDP_1KGP_PHASED_INDEX_METHODS.md`
  state the limit. When no dataset of the search has the tiers, such a row of two or more SNPs
  reads `N carriers (>=1 of K vars, poss. chimera)` in the report's `Observed` column, as in a
  population-level search; the report cannot tell which dataset a row comes from, so beside a
  dataset with the tiers it still reads `N cis carriers`, and the legend and the caption say why.
  Results are unchanged. Pre-existing: 2.8.0 behaves the same; it is fixed in the next release.
  Covered by `PostProcess/test_build_index_only_atomic.py` (the banner) and
  `PostProcess/test_generate_report.py` (the report).
- **A known limit, stated with its scale: on a row of several variants, the rsID and
  allele-frequency lists are not in the order of its variants.** The post-analysis writes a row's
  variant, rsID and AF lists one to one, in the variants' order, and the merge of a cluster's rows
  (`remove_contiguous_samples.get_best_targets`) sorts each of the three on its own, as text, and
  keeps each value once. So in `integrated_results.tsv` and the alternative alignments the n-th
  item of `Variant_rsID` or `Variant_MAF` need not belong to the n-th variant of
  `Variant_info_genome`, and two variants of one frequency show it once. A `--per-sample` search
  of the README's guide (`20bp-NGG-SpCas9` PAM, 4 mismatches, 1 DNA and 1 RNA bulge) on the
  self-test's index rebuilt with `build-index-only --samplesID` lists chr22:16858220 with the
  variants `chr22_16858239_C_G,chr22_16858242_G_T` and the frequencies `0.000981162,0.805338`,
  where the first variant's is 0.805338 and the second's 0.000981162, and two of its seven
  multi-variant rows give one frequency for two variants. On ACTGA (`ACTGAAATCTGTAAGCAGGC`, NRG,
  4 mismatches + 1 DNA + 1 RNA bulge, per sample, genome-wide, on the HGDP + 1000 Genomes index),
  953 of 1,967 rows of two or more SNPs and 118 of 293 rows of an indel and one or more SNPs list
  their rsIDs in another order than their variants (in the highest-CFD columns). The report's MAF,
  the lowest frequency of the list, does not change; the rsIDs its `Variant` column gives after the
  variants, and the rsIDs and AFs of `observed_haplotypes.tsv` and of a site's haplotype table in
  the report, come in the list's order. METHODS section 8 states the limit. Pre-existing: 2.8.0 and
  every release since early 2022 behave so, the 2.1.x releases this changelog lists among them
  (before 2.4.0 each list came in an arbitrary order); it will be fixed in the next release.
  Covered by `PostProcess/test_generate_report.py`.
- **A known limit, stated with its scale: an indel (and, on an index without the dict-less tiers, a
  SNP) whose record has `AF` at another place in `INFO` than the VCF's first record has no allele
  frequency.** CRISPRitz `add-variants`, which writes an index's indel log and classic per-sample
  dictionaries, looks for the `AF` key in the VCF's first record only and reads every record's `AF`
  at that place in `INFO`, so a record that has it elsewhere gets none: its indel has no frequency
  on any index, nor its SNP on an index without the dict-less tiers (with them a SNP's frequency is
  the registry's, which `add-variants` does not write), and the report's MAF leaves it out. A VCF
  merged with `seq_script/merge_panels/merge_vcf_panels.sh` has such records, since a site of one
  source has fewer keys before `AF` than a site of two: the retired `NRG_3_hg38+hg38_1000G2021_HGDP`
  index, built from such a merge, has no `AF` for 14,567,429 of its 15,583,931 indels (224,981 of
  the 232,384 of chr22). Of the shipped indexes only HPRC's chrY is affected, whose VCF has `AF` at
  five places: its log has none for 91,950 of its 102,665 indels. The HGDP1kGP index and HPRC's
  other chromosomes have an `AF` for each of their 17,192,381 and 12,000,599 indels, as the
  self-test's 1000 Genomes index has for its 64,528 indels and 993,880 SNPs. The mega's VCF has no
  `INFO/AF` at all, only `AF_<dataset>` and `AF_max`, so none of its 13,853,761 indels has one and
  the MAF of a row of an indel alone is blank there; the `indel_af.tsv` companion gives their
  per-dataset frequencies. The report's MAF note (and, for a dataset without the dict-less tiers,
  its Search mode and Allele-frequency basis rows and its legend's MAF entry, as the
  `complete-search` banners), the README's MAF column and `--per-sample` text, METHODS sections 3, 7
  and 8, `docs/DESIGN_mega_index.md`, the gaps of `docs/INDEX_PROVENANCE.md` and the VCF
  requirements of `docs/INPUT_FORMATS.md` (with a check) now say so, and the header of
  `merge_vcf_panels.sh` and `seq_script/merge_panels/README.md`, which said that the enricher reads
  the pooled `AF` by its exact key, and a comment of `merge_mega_sites.sh`, which said that the
  plain `AF` it writes gives indels the global `AF_max`, say what the enricher does. Results are
  unchanged. Pre-existing: 2.8.0 behaves the same; it will be fixed in the next CRISPRitz release
  (HPRC's chrY log needs rebuilding with it). Covered by `PostProcess/test_generate_report.py` and
  `PostProcess/test_build_index_only_atomic.py`.
- **METHODS section 2 says what a sites-only panel's registry holds for carriers and
  homozygotes.** It said that the carrier and homozygote counts of the mega index are reported as
  Hardy–Weinberg expectations from the allele frequency, flagged as estimates, and
  `docs/DESIGN_mega_index.md` and `seq_script/merge_panels/README.md` said so too. None is
  estimated: `compile_registry_from_info_af` (`build_mega_registry.py`) stores 0 for both, beside
  the exact allele frequencies, and marks the registry as aggregate (`aggregation` `info_af`, and
  `aggregate_af_only` for each dataset), from which the population-summary companion gives them
  as NA; the report shows no carrier or homozygote count from the registry, and the carrier
  overlay of `targets-scoring` prints the stored 0. METHODS sections 2 and 7, the two docs and the
  README's `targets-scoring` section now say so (the overlay gives each variant's carrier count
  over the panel and each dataset's highest allele frequency, not a carrier count per dataset, as
  the README said). The code is unchanged: it has stored 0 since the mega index came, in 2.5.2,
  whose METHODS already said otherwise. Covered by `PostProcess/test_registry_from_info_af.py`.

### Fixed
- **One indel off-target could be reported on every chromosome, each time with that chromosome's real
  rsID, allele frequency and carrier samples.** Every CRISPRitz targets file carries a `#` header, and
  `post_analisi_indel.sh` copies the merged indel file's first line as the header of each
  per-chromosome subset — it is the only stage that does so, which is why only the indel leg was
  affected. Two filters read that header as a target and dropped it:
  `pam_filter.py` guarded headers with a column count, but the CRISPRitz header has all 10 columns, so
  its `DNA` column was taken for a PAM region and rejected; and the `--max-total-edits` cap compared
  `$10` (`"Total"`) against the numeric cap as strings. With the header gone, `head -1` returned the
  first *data* row, and that single row was prepended to all ~23 subsets. Nothing rejected it, because
  every chromosome's fake genome starts at coordinate 0 with near-uniform windows, so the indel lookup
  (keyed on the fake position alone) resolved it against whichever chromosome was running and overwrote
  the `Chromosome` column with that one. One raw search row therefore became one correct plus N phantom
  report rows, each indistinguishable from a genuine annotated variant off-target. Measured genome-wide
  on `ACTGAAATCTGTAAGCAGGC` at NRG 4+1+1 over HGDP+1kGP: 8 phantom rows on a dedicated NRG index and 11
  on a pamless index, against an identical justified core of indel rows; because the injected row is
  whatever lands at line 1, the count varied with the index. `pam_filter.py` and the edits cap now keep
  `#` lines, `post_analisi_indel.sh` prepends only an actual header, and the indel post-analysis rejects
  and counts rows carrying another contig's label. The injected row also moved one legitimate row's
  CRISPR-Bulge score by one unit in the last place: scoring runs in fixed-size batches, so an extra row
  changes a batch's contents and with it the order of the reductions inside it. That score now agrees
  across indices. The SNP leg was never affected. Covered by
  `PostProcess/test_indel_header_injection.py`.
- **A chromosome with SNPs but no indels no longer aborts the whole search** (a regression in this
  release's CRISPRitz pin, caught before shipping). The enricher builds each chromosome's fake-indel
  genome from the indels some sample carries, so a chromosome with none gets an empty genome;
  `buildTST` cannot index it and leaves that `_INDELS` sub-index without a `.bin`, and CRISPRitz 2.9.1
  -- rightly, for a broken index -- refuses to search such a directory. One SNP-only chromosome
  therefore failed the entire run, reported as an SNP post-analysis failure, where 2.8.0 had silently
  returned the correct answer: no indel off-target there. Reproduced through the real CLI on chr21
  (23,235 SNPs, no indel) + chr22 (2,891 indels): 2.8.0 exits 0 with 284 rows, 2.8.1 without the fix
  exits 1 with none. `pool_index_indels.py` no longer hands an empty fake genome to `buildTST`, and
  `pool_search_indels.py` skips a chromosome whose indel log records no indel, tolerating the missing
  part file for that chromosome only. With the fix the same run exits 0 with a report byte-identical
  to 2.8.0's. An index that is empty although its log lists indels is still a failure, and still
  reported as one. Covered by `PostProcess/test_indel_empty_chromosome.py`.
- **`NRR_3_hg38+hg38_mega` is republished: the 2026-10-27 archive found no indel off-target at
  all.** Measured in the archive itself (downloaded; SHA-256 equal to the Hub's): its `_INDELS` holds
  24 sub-indexes and **0** `.bin` files, and its `log_indels_hg38_mega` logs are header-only. It was
  indexed from a mega fake-indel genome that predates the sites-only fix (#191) -- every chromosome's
  FASTA a bare header -- by a script that bypassed `seq_script/merge_panels/mega_build_indels.sh`;
  `buildTST` crashed on each empty genome and nothing checked. Every NRR mega search therefore
  reported reference and SNP off-targets only, and on 2.8.0 said nothing. It also damaged
  `NRG_3_hg38+hg38_mega` in a shared working directory: an install replaces each `Dictionaries/`
  directory, so it swapped NRG mega's 419 MB of indel logs for the empty ones, and NRG mega then lost
  its indel off-targets too. The replacement's `_INDELS` is rebuilt with CRISPRitz 2.9.1 from the
  2026-10-10 mega kit genome -- which reproduces the published NRG mega `_INDELS` byte-for-byte, 32 of
  32 files -- giving 67 `.bin` across 24 sub-indexes, and its `Dictionaries/` + `samplesIDs/` are
  taken from the published NRG mega archive, so both mega archives now install identical
  dictionaries. Its SNP index is unchanged. Measured genome-wide (sites-only, 4+1+1): the published
  archive gives **0** indel off-target rows on 2.8.0, the replacement 2,081,489 on 2.8.1; and on
  each of the 24 indel sub-indexes, the replacement's NRR targets filtered to NRG equal the
  published NRG mega's exactly (1,131,494). **If you downloaded `NRR_3_hg38+hg38_mega` before
  2.8.1,** remove `genome_library/NRR_3_hg38+hg38_mega{,_INDELS}` and download it again; that also
  restores NRG mega's indel logs where the two share a working directory. Every other published
  variant archive was checked the same way and is complete. A new release gate (T-40) checks this in
  the archive for every variant index.
- **`NRR_3_hg38+hg38_HPRC` is republished so that it defaults to per-sample analysis, like
  `NRG_3_hg38+hg38_HPRC`.** Its index is complete, but its registry was built before the
  `data_type`/`phased` stamps were added, and `complete-search` picks the default mode from that stamp
  (per-sample only on `genotyped-phased`). An NRR HPRC search therefore ran population-level --
  no named carriers, no CONFIRMED cis phasing -- on a phased panel, without saying so, and installing
  it next to NRG HPRC replaced NRG HPRC's stamped registry and changed its default too. The
  replacement changes only `Dictionaries/`, taken from the published NRG HPRC archive; the registry
  data itself is unchanged. Passing `--per-sample` explicitly was never affected.
- **Both `NNN_3` archives are republished, rebuilt with the CRISPRitz 2.9.1 builder.** The archives
  uploaded on 2026-10-25 predate the genomic-`N` flank fix, so they lack every target whose 3-nt
  flank holds an assembly `N`. 2.8.1 never reports a target with an `N` inside its PAM (the
  post-search PAM filter now drops it, and under `--prune-pam` the engine already does), so a 3-nt
  PAM gives the same report on the old and the rebuilt archives: `NGG`, `NRG`, `NRR` and the
  unconstrained `20bp-NNN-NO-PAM` alike. The 2026-10-25 `NNN_3_hg38+hg38_HGDP1kGP`, downloaded and
  searched with NRG `--prune-pam`, reproduced the dedicated NRG report byte-for-byte genome-wide. A
  shorter 3′ PAM under `--prune-pam`, such as `20bp-NG-SpCas9NG`, does change: the old archives miss
  each site whose base just past the PAM is an `N` (at an assembly-gap edge, and in the indel index,
  whose fake-indel genome has an `N` every ~56 bases), and the rebuilt ones report it, as a
  dedicated index for that PAM does. On a synthetic genome with four gap-edge sites, an NG
  `--prune-pam` search gave 1 row on an index from the old builder and 5 on one from the new
  builder, byte-identical to a dedicated NG index. The rebuild uses the old build's own inputs and
  changes the index only where the fix predicts: 1,222 of the reference index's 1,628 bins
  (+145,540 B) and 1,221 of the variant index's (+106,849 B), on exactly the 66 contigs that contain
  `N` and on none of the other 389; and the whole indel index, whose 23 contigs are cut by `N`
  every ~56 bases by construction (222 bins become 245, +2.1 GB). Built independently on two nodes,
  the reference index is byte-identical. Searched genome-wide by a 2.8.1 release candidate (NRG
  `--prune-pam`, 4 mm + 1 DNA + 1 RNA bulge, per-sample), the rebuilt pair matches the gate's
  report in every column but `Annotation_COSMIC`, which
  differs only because that run used freshly downloaded annotations, where COSMIC stays off until
  its licence is attested. The released 2.8.0, searching chr22, gives byte-identical reports on
  the old and the new reference index, with NRG and with NNN, but only because no hit there has an
  `N` in its PAM: use the rebuilt archives with 2.8.1 or later. 2.8.0 keeps the targets with an
  `N` in their PAM that the rebuild adds (it has no `--prune-pam`, ignores the flag, and its
  post-search PAM filter lets the `N` through), and what a search that reaches one does depends on
  its CRISPRitz. On 2.9.1 it fails in post-analysis (`ValueError: substring not found`). On 2.9.0,
  which 2.8.0's image and its `install_from_source.sh` build, the engine writes such a target cut
  short at its `N`, so the search fails in post-analysis (`IndexError: string index out of range`)
  or, with RNA bulges only, finishes and reports the target as a 0-mismatch, 0-bulge off-target
  with a 21- or 22-nt protospacer+PAM. On the synthetic gap-edge genome, reference-only, an index
  from the new builder made 2.8.0 exit 1 at 2 mismatches + 1 DNA + 1 RNA bulge, with NRG and with
  NNN, on either engine (one from the old builder gave 2 and 3 rows); on 2.9.0 with
  `--bDNA 0 --bRNA 1` it gave 4 NRG rows where the old builder's index gave 2, three of them such
  short rows with CFD 0.259, one in place of a real 1-bulge alignment. Genome-wide, the NRR search
  in the post-search PAM filter entry below meets 1 such reference target, 1 variant target and
  93,820 indel targets. If you searched the rebuilt `NNN_3` with 2.8.0, run that search again with
  2.8.1. The variant archive grows from 108.4 GB to 109.5 GB. Hugging Face has served the rebuilt
  archives since 2026-10-30. With 2.8.1 there is nothing to do unless you downloaded `NNN_3`
  before the replacement on 2026-10-30 and search it with `--prune-pam` and a PAM shorter than 3
  nt; then remove `genome_library/NNN_3_hg38`, `genome_library/NNN_3_hg38+hg38_HGDP1kGP` and its
  `_INDELS`, and download again (`download` skips an index that is already installed).
- **A PAM of any length other than 3 could be handed the 3-nt PAM-less index.** The three scans
  that pick an index (the reference scan, the variant scan and the precomputed-variant detection)
  globbed the requested PAM and then a literal `NNN_*`, whatever the PAM's length, and kept the
  smallest sufficient bulge budget. So with `NNN_3` installed (published since 2026-10-25) and no
  dedicated index winning, a search with SaCas9's 6-nt `NNGRRT`, SpCas9-NG's `NG` or Cas12a's
  `TTTV` ran on an index whose stored flank cannot serve it. Reproduced on the released 2.8.0 with a
  synthetic genome holding three planted SaCas9 sites: with `NNN_2` installed the run reported
  `Reference index present (NNN_2_synth)` and failed with `Error on reading the PAM, input char is
  not an IUPAC nucleotide`, where an empty library builds `NNGRRT_2` and finds all three. 2.8.1,
  in the same layout, builds `NNGRRT_2` and returns the same three rows. The scans now share one
  resolver that takes a PAM-less stand-in only when its name is all `N` at the PAM's width (or,
  under `--prune-pam` and for a 3′ PAM, wider) and it was built for the PAM's side of the guide
  (next entry), so a dedicated PAM that merely begins with `N`s, such as the shipped `NNNRRT`, is
  never taken for one either. Covered by `PostProcess/test_index_resolution.py`.
- **A PAM-less index built for the other side of the guide could answer a search, which then found
  nothing.** An all-`N` index stores its flank where the PAM file it was built from puts the PAM,
  and the engine reads the flank on the side the search's PAM file names. The resolver above took
  an all-`N` index of the PAM's width (or wider, with `--prune-pam` and a 3′ PAM) whichever side
  it was built for, as 2.8.0's `NNN_*` glob did, and the search then reported no off-targets and
  exited 0. On a synthetic genome with three planted sites, reference-only, 2 mismatches + 1 DNA +
  1 RNA bulge, Cas12a's `TTTV` on an `NNNN_2` index built from an all-`N` 3′ PAM file returned
  0 rows, with or without `--prune-pam` and with or without `--index-path`, where a dedicated
  `TTTV_2` index finds the three; an 8-nt 5′ PAM on a 3′-built `NNNNNNNN_2`, and `NGG` on an
  `NNN_2` or `NNNN_2` built from a 5′ all-`N` PAM file, returned 0 rows too. CRISPRme now reads
  the side from the index's `.pam_build`, which `build-index-only` writes and, from this release,
  so does a search that builds an all-`N` index; an index without one (built by a search before
  2.8.1, or by `crispritz.py` itself) is taken for a 3′ build at the widths of the published ones,
  3 and 8, and at any other width the search builds its dedicated index (the entry on an index
  taken only for the side it was built for, below). With no index left for its side, a search
  builds a dedicated one, as for any missing index, and under `--index-path` stops with
  `no matching reference index`, naming the all-`N` index it did not take and why. Without
  `--index-path`, the `TTTV` and `NGG` searches above now build the dedicated index and report its
  3 rows, byte-identical to a search of it. A search whose PAM is itself all-`N` builds its index
  under the name an all-`N` index of that width has, and searched one built for the other side
  under that name as if it had built it: `20bp-NNN-NO-PAM` on an `NNN_2` built from a 5′ all-`N`
  PAM file reported no off-targets and exited 0, where a 3′ build gives 4 rows. It now stops there
  and says to move that index out of the way. A 3′ build serving a 3′ PAM, and a 5′ build serving
  a 5′ PAM of its width, are chosen as before, one a search built included: on a synthetic genome
  with three planted Cas12a sites, the 5′ `NNNN_2` a `27bp-NNNN-5p` search builds serves `TTTV`
  with or without `--prune-pam` under `--index-path`, 3 rows, the dedicated `TTTV_2` index's. A 5′
  all-`N` index that a search built before 2.8.1 records no side, so a 5′ search refuses it and
  says so: run `build-index-only` with the PAM file it was built from, the same `--genome` and its
  bulges, which finds it in place and records the side. `build-index-only` no longer records a PAM
  file over an all-`N` index whose `.pam_build` names the other side, which, now that a search
  reads that record, would make a search on its real side refuse it and one on the other side find
  nothing; it stops and says to build in another working directory. No published index searched
  with a shipped PAM file reached any of this (the shipped 5′ PAMs are 4 nt, and the published
  PAM-less flanks are 3 and 8 nt, all built for the 3′ side). A PAM-less index you built yourself
  could: if you searched one with a PAM on the other side, run that search again. Covered by
  `PostProcess/test_index_resolution.py`, `PostProcess/test_bulge_cap.py` and
  `PostProcess/test_build_index_only_atomic.py`.
- **A PAM-less index searched without `--prune-pam` kept targets with a genomic `N` in their PAM,
  and the search then crashed.** With `NNN_3` serving an NRR (or NRG) search and no `--prune-pam`,
  the post-search PAM filter (`pam_filter.py`, in releases since 2.2.0) treated an `N` in a target's
  PAM region as matching any base, so it kept targets a dedicated index never returns: CRISPRitz
  stores a genomic `N` as a code that matches nothing, so neither a dedicated nor a `--prune-pam`
  search reports one. Genome-wide, NRR, 4 mismatches + 1 + 1 bulge, HGDP+1kGP per-sample, the
  filter kept 1 extra reference target, 1 extra variant target and 93,820 extra indel targets (an
  assembly-gap edge, and the `N`s that separate the windows of the fake-indel genome), exactly the
  targets with an `N` in their PAM region, and the SNP post-analysis then failed on the `N`
  (`cluster.dict.py`: `ValueError: substring not found`, `ERROR: off-targets post-analysis (snps)
  failed`). A reference-only search failed the same way, in its reference post-analysis (`ERROR:
  off-targets post-analysis (reference) failed`). Only an index that stores such targets reaches
  this: a PAM-less one from the CRISPRitz 2.9.1 builder, such as the republished `NNN_3` (above);
  the old builder dropped them. A target's `N` now matches nothing, at every position of the PAM,
  even an unconstrained one, so an unconstrained `20bp-NNN-NO-PAM` search on such an index, which
  kept these targets and failed the same way, now drops them too. The filter is unchanged for
  dedicated indexes and under `--prune-pam`, whose targets never carry an `N` there (it dropped 0
  in every such arm). With `--prune-pam`, `NNN_3` and `NNNNNNNN_3` were already byte-identical to
  a dedicated NRR index on that search: 386,012 rows, the alternative alignments and observed
  haplotypes included. Covered by `PostProcess/test_indel_header_injection.py`.
- **`complete-test` passed a flag that did not exist.** The brute-force benchmark command included
  `--bmax`, which `complete-search` has never had — it recomputes bMax from the per-type budgets — so
  the argument was discarded unread. Harmless until unknown flags started being rejected, which then
  broke every benchmark. The flag allowlist was checked against the help text and against what the
  parser reads, but not against the commands this repo builds, so a dead flag in a caller could only
  surface by running the benchmarks; that gap is now a static test.
- **The documented install self-test tested nothing, and `validate-test` then failed with a
  traceback.** README §1.3, §3.2 and the installer's last line presented `crisprme.py complete-test
  --chrom chr22` followed by `validate-test` as a self-test that downloads chr22 and compares with
  the brute-force ground truth. That benchmark mode runs only with `CRISPRME_ALLOW_ONDEMAND_BUILD=1`
  set (`--vcf_dataset` does not select it, whatever the help said). Without it, `complete-test`
  looked for a prebuilt index, found none in a fresh folder and exited 0 having tested nothing,
  and `validate-test` then stopped with two tracebacks (`FileNotFoundError: ... vcf.config.test.txt`,
  then `OSError`). The README (§1.3, §2.2.2, §2.2.3, §3.2), both commands' `--help` and the
  installer now give the commands CI runs:
  `CRISPRME_ALLOW_ONDEMAND_BUILD=1 crisprme.py complete-test --chrom chr22 --vcf_dataset 1000G`,
  then `crisprme.py validate-test --chrom chr22`. Where that run has left nothing to check,
  `validate-test` says so in one line, naming both commands for the chromosome it was asked to
  check (the chr22 pair without `--chrom`) and saying that in Docker its variable goes with `-e`,
  and exits 1; any failed validation now exits 1 without `crisprme.py`'s `OSError` traceback after
  `validate.py`'s own report. A benchmark whose brute-force reference covers chr22 only (the
  Cas12a one and the single-"n edits" ones) is now skipped on any other chromosome and on a
  whole-genome run, saying `SKIPPED: its brute-force reference covers chr22 only`. It used to be
  compared there, with every CRISPRme site off chr22 counted as extra: after the chr21 run that
  line names, `validate-test --chrom chr21` reported `[cas12a_hbg] FAILED` (223 extra sites) where
  it now ends `1 passed, 0 failed, 5 skipped`. Validation on chr22 is unchanged. In Docker, the
  only route on macOS, the variable must be passed with `-e`: set in front of `docker run` it
  stays on the host, and the run tests nothing. README §1.2 and `docs/DOCKER_QUICKSTART.md` now
  give that form. Also in the README: §1.3 clones the `v2.8.1` tag rather than `main`, and §1.2
  says the image is built from source (CRISPRitz v2.9.1 compiled in it) and gives its download and
  on-disk sizes separately. Covered by `PostProcess/test_validate_test.py`.
- **An index still being built could be selected as if it were finished.** A build takes
  `<dir>.lock`, and `index-genome` then creates the index directory and writes one `.bin` per
  chromosome, so for the whole build a half-written index is indistinguishable from a complete one
  by name. A run starting in that
  window printed `Reference index present` and searched a **partial** index — fewer off-targets, no
  error, exit 0 — and did so even when a complete, sufficient index sat right beside it. The scan now
  skips any candidate whose build lock is held, and any directory holding no `.bin` at all, and the
  run waits for that build instead, under `--index-path` as well. A symlinked index directory, as a
  library shared between working directories uses, still counts: the check follows the link. Index
  choice on a finished library is unchanged.
- **A build lock left by a run that died stopped every later run that needed the same build.** A run
  that builds a shared output (an enriched genome, a genome index, an INDELS index) holds
  `<output>.lock` while it writes it, and every other run that needs that output waits for the lock
  to go. In 2.8.0 that wait had no end: a run killed mid-build (SIGKILL, an OOM kill, a container
  stop, a reboot) left its lock behind, and a later run needing the same output waited on it
  forever, or, for an index, took the partial one for finished (above). Each lock now holds an owner
  record (host, PID namespace, PID, start time), refreshed every minute while the build runs. When
  the owner is gone on this host (it died, is a zombie, or its PID now belongs to another process),
  the next run takes the lock over and builds the output again from scratch; from another host or
  container, where the PID cannot be checked, it does so once the record has gone 15 minutes without
  a refresh. With `--index-path`, which never builds a reference index, the run stops instead with
  an ERROR that names the partial index and the command that clears it. A build stopped on this host
  (Ctrl-Z, `kill -STOP`, a SLURM suspend) is waited for however long. A run ended by Ctrl-C, a
  SIGTERM (`scancel`, a SLURM time limit) or a hangup leaves its lock for that takeover rather than
  releasing it, since a released lock would hand on a half-written output as finished. A build
  stopped for 15 minutes or more where no other run can check it (in another container or on
  another host) can lose its lock; resumed, it ends with an ERROR once its build step returns and
  uses none of what it built, so a library shared across hosts or containers must not have a build
  suspended that long. While it waits, a run notes in `log_verbose.txt` who holds the lock and what
  to do. A lock with no record, as 2.8.0 takes, is never broken by age, since a 2.8.0 run may still
  be building under it: the run waits and says every five minutes how to clear it. Covered by
  `PostProcess/test_index_resolution.py`.
- **The RAM thread guard sized a search served by a PAM-less index as if it were a dedicated one.**
  It chose its memory model from the requested PAM, so an NRG bulge search served by an all-N index
  (`NNN_3`, or a wider flank under `--prune-pam`) was sized as constrained (4 GB base, 0.8 GB per
  thread) instead of PAM-less (12 GB base, 2.2 GB per thread), and could be given more threads than
  the host's memory holds, and be OOM-killed. A bulge search is now sized by the index it actually
  searches. So is a variant search with no bulges that lacks its whole enriched genome (a
  download), which searches its variant index. Any other search with no bulges scans the genome
  FASTA, opens no index, and is sized as before. Only the thread count changes, never the results.

- **A downloaded variant index warned that its genotype store was missing right after installing
  it.** The "no `genotypes_<vcf>.tar.gz` companion ... per-sample Samples degraded" note was written
  as the `else:` of the shared-reference-index `try`/`except`, so it fired whenever *that* download
  succeeded — an unrelated event. It was wrong in both directions: false on a complete download, and
  silent when the companion really was absent but the reference index was already installed, so its
  block never ran. The note now follows the companion fetch's own result, and distinguishes
  "genuinely absent" from "the fetch errored" — the latter already reports itself, and calling it
  absent would be wrong.

- **A search that legitimately found nothing was reported as a failure.** Three places sized a
  CRISPRitz targets file with `-s` or `wc -l` to decide whether the search had found anything. Those
  are equivalent to counting data rows only while the file has no header — which was true only
  because the filter bug above was stripping it. With the header correctly preserved, a zero-hit
  search left a header-only file, the zero-hit short-circuit saw a non-empty file, and the run
  aborted with "the search found 1 raw off-target(s) but none survived post-analysis". The two bugs
  had been cancelling. All three sites now count non-header rows through one helper, which also
  handles a file that has no header rather than blind-skipping its first line. Caught by the Cas12a
  brute-force benchmark, whose chr22 result is genuinely empty.

- **`complete-test --chrom chr22` searched the whole genome, and the source installer started on
  macOS.** The README (§2.2.2, §3.2.1) and `complete-test --help` presented `complete-test --chrom
  chr22` against a downloaded index as a quicker single-chromosome test, but the prebuilt-index
  smoke never received `--chrom`: it searched every chromosome of the index's genome, so with
  `NRG_3_hg38+hg38_HGDP1kGP` it was a genome-wide variant search, with the 64 GB or more of RAM
  that §0 gives for one. On a two-chromosome genome it returned 10 rows on `chr22` and 10 on
  `chrA`. The smoke now searches only the chromosome `--chrom` names. It runs the same search from
  a per-run view of the working directory in which that chromosome is the only one: a hidden
  `.crisprme-test-view-*` folder of symlinks to its FASTA, to its trees in every index the search
  could pick, to its `fake<chrom>` INDELS sub-index, to its fake-indel genome and to its indel
  log, so that the search does not load another chromosome's search or post-analysis data (the
  dataset-wide tier folders are linked whole). For a variant index without an INDELS index (that
  of a dataset without indels, or one whose build stopped before it), the search builds a
  throwaway one for that chromosome inside the folder, or finds no indel there to index, and
  leaves none in `genome_library`. Before it makes the folder, the smoke refuses each part of its
  dataset that `build-index-only` has marked unfinished and `complete-search` would refuse (see
  below), in one ERROR line that names its real path, with nothing searched and no
  `Results/crisprme-test-out` left behind. Nothing is copied, and the folder is removed when the
  run ends, on Ctrl-C, SIGTERM or SIGHUP too: a SIGTERM or SIGHUP sent to `complete_test.py` alone
  is passed on to the search, which is waited for (and killed after 10 s) before the folder goes,
  and a hangup that `nohup` ignores stays ignored, as without `--chrom`. Its rows are the
  whole-genome smoke's rows on that chromosome, identical in every column, as checked on each
  sequence of a synthetic genome with a phased variant index and on `chr22` with
  `NRG_3_hg38+hg38_HGDP1kGP`, both searched per-sample, the default on such an index; the one file
  it leaves out is the header-only `fake<chrom>_*.indel_snp_cooc.tsv` that the whole-genome run
  writes for a chromosome with indels but no indel off-target, which holds no row. There
  `--chrom chr22` took about 1.5 to 4 minutes and 5 GB of RAM at 16 threads, the whole genome 11
  to 17 minutes and 44 to 50 GB. A chromosome holding no off-target of the example guide, such as
  `chrM` with `NRG_3_hg38`, passes with an empty result, as the whole-genome test reports none
  there either (complete-search then writes no `integrated_results.tsv`). A chromosome the genome
  or the index does not have is refused before anything is searched, in one line naming those both
  have, and `crisprme.py` exits 2 without a traceback. Without `--chrom`, or with `--chrom all`,
  the smoke runs the search it always ran, and the benchmark mode is unchanged. The help, the
  no-index hints and the README say so, and the README's "Full Genome Test" runs
  `complete-test --thread 4`, since `complete-test` with no option only prints its help.
  `install_from_source.sh` now stops at once on macOS with one line pointing to Docker (README
  §1.2), instead of failing later on a solver or compiler error, and `environment.yml` and the
  README's source-build entry say Linux only. Covered by `PostProcess/test_complete_test.py` and
  `PostProcess/test_install_from_source.py`.

- **`--prune-pam` could run with pruning silently off.** `crispritz.py` builds its search argv
  explicitly from the options it recognizes, so a CRISPRitz build without `--prune-pam` does not
  reject the flag — it drops it, and the search runs unpruned over the whole PAM-less index. The
  post-search PAM filter still enforces the requested PAM, so the report came out right while the run
  read an index many times larger than it needed. `complete-search` now checks that the installed
  CRISPRitz can prune and fails with an explanation if not, rather than assuming it from the pin,
  which is overridable through `CRISPRITZ_REF` or an engine already on `PATH`.

- **"The error log is empty" checks passed after any job that ran to the end, whatever it had
  written.** A job that runs to the end renames its `log_error.txt` to `log_error_no_check.txt`, so
  a check of `log_error.txt` could no longer fail by then. CI's read-only-container smoke checked
  only `log_error.txt`, and printed PASS on this release's candidate although the reference-only
  search it runs had no report (`report.log`: a TypeError) and its error log held three lines. The
  smoke now checks whichever of the two files exists and requires the job's `<job>_report.zip`; the
  zero-hit check, whose job ends before the rename, checks the same way. `docs/DOCKER_QUICKSTART.md`
  and the web interface guide now name `log_error_no_check.txt` where they tell you to check the
  error log of a job that ran to the end, and say that a job still running, failed or with no
  off-targets keeps `log_error.txt`. The rename itself is unchanged.

- **The report manifest could stamp a hardcoded version.** The run manifest in the shareable report
  fell back to a literal `"2.8.0"` when a run recorded no version. It now resolves the real version,
  and reports `unknown` rather than a plausible-but-wrong number — this manifest is the IND
  traceability record.

- **The report manifest named the wrong search-engine version.** `crispritz_note` carried a hardcoded
  `(CRISPRitz v2.8.3)`, a literal that outlived three CRISPRitz pins, so every shipped manifest named
  an engine the install had never run. The manifest now carries a `crispritz_version` field resolved
  from the `crispritz.py` on `PATH` — the engine that actually ran — and `null` when it cannot be
  determined, with no version in the prose. The test that should have caught this was asserting the
  stale literal; it now asserts the property instead.

- **A search with no variant off-target never got its shareable report.** Every reference-only
  search, and any search in which no row carries a variant haplotype, ended its report step with
  `TypeError: cannot unpack non-iterable NoneType object` in `report.log`: no `<job>_report.zip`,
  `Building shareable report FAILED (non-fatal)` in `log.txt`, the "shareable report generation
  failed" WARNING in the error log, and `generate-report --result-dir` exiting 1 with the same
  traceback, while the web results page kept waiting for a report that never came. The
  observed-haplotypes companion returned a bare `None` when it found no variant haplotype to group,
  where its docstring and every other exit return `(None, [], {})`, the three values the report
  unpacks. It now returns those, and the report is built without the companion, as it already was
  for a variant search with no multi-haplotype site. Search results never changed, and the report
  of a search with variant rows is byte-identical. Present since 2.7.10. Covered by
  `PostProcess/test_generate_report.py`.

- **History page no longer crashes when a job is queued.** With every search slot busy, a queued
  job's `log.txt` exists (a queue note is written to it) but has no `Job\tStart` line yet;
  `read_job_info()` raised `StopIteration`, which aborted the entire History-table build — so running
  and completed jobs could not be opened until the queued job started. `read_job_info()` now returns
  `None` for a not-yet-started or missing log (the caller falls back to the `.Params.txt` mtime, the
  same stand-in the assembly path uses), and the per-job summary loop is guarded so a single not-ready
  or malformed job directory can never take down the whole table. (Reported by M. Tognon.)
- **Index/annotation download install is atomic and always structurally verified.** A partial or
  failed download never touches the final install location, and a failure surfaces a copy-paste log
  plus a Retry action instead of a silent or corrupt state. Two genome/VCF download behaviours change
  from 2.8.0, which wrote straight into `Genomes/<ref>` and `VCFs/<ds>`: a folder that is a symbolic
  link (for example to a shared copy) is replaced by the downloaded copy, with a NOTE, and the shared
  copy is left untouched; a folder that is a mount point (including a bind mount) is refused before
  the download starts, saying what to do instead, so `download --what all` stops at that component.
- **Every annotated search left its sorted annotations behind, about 1.7 GB each.**
  `complete-search` and `assembly-search` sort each `--annotation` / `--gene_annotation` into a new
  `crisprme_annot_*` folder in the working directory (`$CRISPRME_TMPDIR` if set), and
  `--personal_annotation` into a `crisprme_pannot_*` one. Nothing removed them, and the pipeline's
  last step gunzipped the plain `.bed` back into them: with the shipped annotations about 280 MB
  plus 1.4 GB per search, and one benchmark-mode `complete-test` left four of them, 3.4 GB. They
  are now removed once the command that made them is over, whether it succeeded or failed (an
  input error, a failed pipeline), and never another job's; `--keep-intermediates` keeps them, as
  its `complete-search --help` entry and the README now say. They are kept, as before, whenever
  the pipeline may still be running on them: after a Ctrl-C, on which `crisprme.py` kills only
  the `/bin/sh` it started the pipeline with (where `/bin/sh` is dash, as in the Docker and
  Apptainer images, an interrupt sent to `crisprme.py` alone, such as `kill -INT <pid>`, leaves
  the whole pipeline running to the end; where it is bash, as on RHEL-family hosts, that shell is
  the pipeline's own, and the step then running finishes on its own), and after a signal that
  killed the pipeline itself, whose workers may outlive it. Search results are unchanged. Covered
  by `PostProcess/test_annotation_tmpdirs.py`.
- **The default variant index could fail to download through Xet, and CRISPRme then gave up at
  once.** On 2026-10-31, `download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP` failed three
  times in a row over an hour, from two networks and with hf_xet 1.7.0 and 1.8.0, with `File
  reconstruction error: CAS Client Error: Format error: I/O error: error decoding response body`
  after 7 to 13 minutes each (the failure is on Hugging Face's side and intermittent: the same
  file downloaded through Xet on 2026-10-18 and again later that evening). Fetched over plain
  HTTPS instead (Hugging Face serves every Xet file through its xet-bridge as well), it arrived
  whole: 27,703,992,074 bytes, SHA-256 `e7046dcd…`, as published. CRISPRme retried
  only a rate limit (HTTP 429) or a server error (5xx), so the first failure ended the download with
  that message and nothing else. Now an hf_xet transfer failure switches the download to plain
  HTTPS for the rest of the run (slower, but independent of Xet) and retries; other transport
  failures (a reset connection, a timeout, a body cut short in requests or httpx) are retried with
  the same backoff; and once the retries run out the error says it is usually a flaky network and
  to re-run the same command. HTTP answers a retry would only repeat (401, 403, 404) and local
  errors such as a full disk still fail at once. A download that works through Xet is unchanged.
- **Pipeline failure summaries are binary-safe.** `summarize_pipeline_failure()` now reads `log.txt`
  and `log_error.txt` with `errors="replace"`, so a single non-UTF-8 byte in a tool's captured
  stderr no longer raises `UnicodeDecodeError` and masks the real error.
- **`--personal_annotation` together with `--annotation` stopped the search before it began.**
  `complete-search` exited 1 with `cat: .../crisprme_annot_*/annotation.sorted.bed: No such file
  or directory` and `Error: Concatenating files failed`, so a personal annotation could only be
  used without `--annotation`. The personal intervals are added to the sorted `--annotation` as
  plain text, and since 2.4.0 that sorted file is bgzipped in its temp folder, which removes the
  plain copy the merge then looked for. The sorted annotation is now decompressed for the merge.
  The search's tables equal those of the same search given, as its only `--annotation`, the two
  files concatenated by hand with each personal name tagged `_personal`: the personal intervals
  appear in `Annotation_personal`. `--annotation` alone and `--personal_annotation` alone are
  unchanged. Covered by `PostProcess/test_annotation_tmpdirs.py`.
- **A search without annotations wrote two `gzip` errors to its error log.** Every
  `complete-search` run without `--annotation` or `--gene_annotation`, and every web search with no
  annotation enabled, ended with `gzip: .../vuoto.txt: unknown suffix -- ignored` twice in
  `log_error_no_check.txt`, so its error log was never empty although nothing had failed. The
  pipeline's last step gunzips the annotation files back beside their `.gz`; it compared each
  path with the bare name of the empty placeholder, `vuoto.txt`, which a path never equals, and so
  gunzipped the placeholder too. It now compares the basename, as the annotation step already
  does. The noise came after the last stage check, so no run failed because of it, and an
  annotated search still gets its `.bed` back beside the `.gz`. Present since 2.2.0. Covered by
  `PostProcess/test_index_resolution.py`.

- **`scorer-env create` failed whenever conda was the env manager, so a Bioconda or source install
  could not run a search.** The CRISPR-Bulge env was created with
  `conda create -y -n cbulge python=3.10 -c conda-forge tensorflow-cpu=2.13 …`. conda parses the
  package specs as one positional argument, so the specs after `-c` were rejected as "unrecognized
  arguments" (conda 25.3 and 26.9 alike; mamba and micromamba accept either order). conda is the
  manager whenever the env was activated with `conda activate` (it sets `CONDA_EXE`), which covers a
  Bioconda install and `install_from_source.sh`: the scorer env was never created, and every
  `complete-search` then refused to start (or ran CFD-only under `CRISPRME_SKIP_SCORER_ENV=1`; a
  search started from the web interface ran CFD-only). The command now puts
  every option before the package specs. The Docker image builds the env with micromamba and was not
  affected. Found by installing `crisprme` 2.8.0 from Bioconda the way a new user would (fresh
  Miniforge, conda 26.9.4); with the fix, `scorer-env create`, its compute self-test and real
  searches all pass.
- **An interrupted or failed `build-index-only` no longer leaves a partial index or enrichment that
  later runs use.** Every later run, `build-index-only` and `complete-search` alike, takes an index
  directory or an enriched genome that is there for a finished one. `build-index-only` ran the engine
  straight into `genome_library/<PAM>_<N>_<genome>[+<vcf>[_INDELS]]`, and moved the enriched genome
  into place before the rest of the enrichment, so a build that ended early (Ctrl-C, an OOM kill, a
  container stop, an engine error) left a partial one that every later search used: a reference
  build killed after `chr1` was searched without `chr2`, with fewer off-targets and no error, and a
  run killed while it compressed the dictionaries left them half compressed beside a finished-looking
  enriched genome. Each index is now built by the same command in a hidden directory of
  `genome_library` and renamed into place once the command has succeeded. A rename within one
  directory never crosses a mount point, so a `genome_library` mounted on its own (`docker -v`,
  `apptainer --bind`) works as before, and the index takes `genome_library`'s setgid group and
  default ACL, as one built in place does. The enrichment is staged the same way: its dictionaries
  are compressed and its fake-indel genome gets its SNP codes (see below) before anything moves, and
  each set of dictionaries is put together whole in a hidden directory of `Dictionaries` (with
  `Dictionaries` mounted on its own, that copies them); then the fake-indel genomes, the enriched
  genome and the dictionaries are renamed into place, the indel logs last, in moments, a Ctrl-C,
  SIGTERM or SIGHUP meanwhile waiting for the last rename. The enriched genome goes in marked
  unfinished (a hidden `.building` in it) until the rest is in place: `build-index-only` does not
  take a marked one for finished, and a search refuses it and needs both sets of dictionaries, so
  a run that ends before then, even by a SIGKILL between two renames, leaves nothing either uses: a
  search refuses the dataset up front, and the next `build-index-only` enriches again. (The
  dictionaries alone used to let a search start, which enriched the genome itself, failed on the
  fake-indel genomes in its way and left its own dictionaries, uncompressed, which later runs kept.)
  A compression that fails now fails the build and puts nothing in place: `pigz`'s exit status was
  ignored, and a Python gzip (used without `pigz`) that failed ended the run with a traceback, the
  enriched genome already in place. An index or an enrichment that another run finished meanwhile
  is kept, with a line on stdout. Ctrl-C, SIGTERM and SIGHUP still end the run by that signal, once
  it has removed its build directories; a second one meanwhile (Ctrl-C pressed twice) waits for
  that. On success the files, their contents (the `.gz` files included), permissions, sidecars,
  stdout lines and exit status are as before.
- **`build-index-only` stops while a search holds the build lock of what it uses or builds.**
  `complete-search` holds `<output>.lock` while it enriches a genome or builds an index, and never
  uses an output whose lock is there. `build-index-only` looked at no lock: started while a search
  built the variant index, it took the half-written index for present and exited 0; started while a
  search enriched the genome, it moved its own enrichment into the places the search was writing.
  Now, while such a lock sits beside the reference, variant or INDELS index or the enriched genome
  that it would use or build, it neither uses nor builds that output: it ends with an ERROR on stderr
  and exit status 1, which names the lock and its owner and says what to do: wait for that run to
  finish, or, if no run is building it, remove the lock and what it guards and run `build-index-only`
  again (not `complete-search`: an enrichment it makes as it takes a dead run's lock over has no SNP
  codes on its fake-indel genome). In the web Settings page the build then shows as failed, with
  that ERROR. The locks are checked before anything is built, before each step, and just before a
  finished build or enrichment goes into place. With no lock there nothing changes, and
  `build-index-only` takes no lock of its own.
- **A chromosome whose INDELS index fails to build now fails the build.** `pool_index_indels.py`
  ignored each chromosome's `index-genome` exit status, so a chromosome that failed (a `buildTST`
  error, an OOM kill) left a partial `_INDELS` index, which `build-index-only` took for built. It now
  prints `Indexing failed for INDELs in <chromosomes>` on stdout and exits 1, and `build-index-only`
  puts no index in place. The search ignores this status, as it goes by the engine's stderr, and is
  unchanged.
- **`build-index-only` checks `--vcf` before it builds the reference index.** A `--vcf` that names no
  dataset, or a listing that names more than one, failed only once the reference index had been
  built; it now fails first, with the same message. A valid `--vcf` builds what it built before.
- **The SNP overlay can no longer leave a fake-indel FASTA cut short** (since 2.5.0). Before it indexed
  the fake-indel genome, `build-index-only` overlaid the SNP codes onto each
  `Genomes/<genome>+<vcf>_INDELS/fake_<vcf>_<chr>/fake<chr>.fa`, rewriting it in place. A run stopped
  in that write left the file empty or cut short, and the next run indexed it as a whole one (an empty
  one as `No indels on <chr>`), so the INDELS index had no indel target past the cut, with no error;
  a write that failed (a full disk) was a WARNING, and the same run indexed what was left. The
  fake-indel genome of an enrichment the run makes is now overlaid before it goes into place (next
  entry); one already in place (a `build-inputs` kit, a run of an earlier release) is still overlaid
  in place, now as the run begins. Either way each FASTA is written in a hidden directory beside
  the original and renamed over it, with the same contents, mode and owner (one in place none of
  whose bases changes, as on a dataset this release enriched, is left as it is), and an overlay
  that fails fails the build with an ERROR: nothing of an enrichment the run makes goes into
  place, and no INDELS index is built. A FASTA that an earlier release cut short is not detected:
  if a `build-index-only` run was stopped while it overlaid, remove
  `genome_library/<PAM>_<N>_<genome>+<vcf>_INDELS`, `Genomes/<genome>+<vcf>` and
  `Genomes/<genome>+<vcf>_INDELS`, then run `build-index-only` again, which enriches again, and do
  not search with that dataset until it has finished.
- **A `build-index-only --vcf` run stopped after the enrichment no longer leaves the INDELS index blind
  to the SNPs beside an indel** (since 2.5.0). `build-index-only` added the SNP codes to the
  fake-indel genome only just before it indexed it, after the tier emission, the indel genotype store
  and the variant index, most of a long build. A run stopped in between (killed, or its variant index
  failing) left the fake-indel genome in place without them, with the dictionaries that let a
  `complete-search` start. That search built the missing INDELS index from it (a search adds no SNP
  codes), and every later `build-index-only` printed `Indels index already present` and kept that
  index: off-targets that need both an indel and a nearby SNP were missing or scored too low, with no
  error. A search that built the INDELS index while `build-index-only` built the variant index made
  the same one. With the real engine, a toy site that needs an insertion and an SNP 5 bp from it on
  one haplotype came out, after a kill while the variant index was built, with 4 mismatches and 1
  bulge, CFD 0.197 and one of its two rsIDs, instead of 3 and 1, CFD 0.281 and both, and was lost at
  `--mm 3`. The SNP codes now go onto the fake-indel genome in the enrichment's hidden build
  directory, before any of the enrichment is in place, so an INDELS index a search builds from it is
  the one `build-index-only` builds. An enrichment already in place without them (a kit, a run of an
  earlier release, a search's own) got them only at that point too, so a run stopped before it left
  the same: it now gets them as the run begins, whenever its INDELS index is still to be built.
  From just before the first FASTA that changes until the overlay ends, the fake-indel genome is
  marked unfinished, as the tiers are (next entry), and a search that would build an INDELS index
  from it (the variant index it searches has none yet) refuses it until a `build-index-only` run
  has finished the codes. `complete-search` refuses it before it starts, unless the dataset has
  another variant index that serves the search with its INDELS index in place: the search then
  goes on, and if the variant index it picks has no INDELS index, it fails as it comes to build
  one, before it builds anything (`CRISPRme run failed`, with that ERROR in its `log_error.txt`),
  as a search launched from the web does (next entry). A search whose INDELS index is in place
  does not read that genome, and the mark does not hold it up. A fake-indel genome none of whose
  FASTAs changes, as when a second PAM or bulge budget is built on a dataset this release
  enriched, is neither marked nor rewritten. If a search built, from that genome meanwhile, the
  INDELS index of the PAM and bulge budget that the `build-index-only` run builds (one that got
  past its check of the mark before the run set it), that index is removed once the overlay ends,
  and built again; that search's own results come from it, with exit status 0, so run such a
  search again once the `build-index-only` run has finished. One that a search built while the
  overlay changed no FASTA, with no mark there, came from the FASTAs the run indexes, and is kept
  (`Indels index already present`). An INDELS index of another PAM or bulge budget that a search
  built that way is not removed: remove it (`genome_library/<PAM>_<N>_<genome>+<vcf>_INDELS`) once
  that run has finished, and the next search or `build-index-only` run with that PAM and budget
  builds it again. With `CRISPRME_INDEL_SNP=0`, which puts no SNP code on, `build-index-only` does
  not build an INDELS index from a fake-indel genome that a stopped run left marked: it stops with
  an ERROR saying to run it again with the overlay on (the default). On success the files, their
  contents, stdout lines (the overlay's still print just before the INDELS index is built) and
  exit status are as before, but for a re-run on a fake-indel genome that carries its codes
  already: its overlay changes no base and rewrites none of its FASTAs, and where its lines said
  `0 SNP overlays` and `0 flank bases IUPAC-coded`, as if it had found nothing, they now say how
  many bases were coded already.  If, with an earlier release, a `build-index-only --vcf` run
  stopped after its enrichment and a `complete-search` with that dataset ran before
  `build-index-only` was run again to the end, or a search built the dataset's INDELS index while
  `build-index-only` ran, remove `genome_library/<PAM>_<N>_<genome>+<vcf>_INDELS` and run
  `build-index-only` again: it first adds the codes to the fake-indel genome in place, then builds
  the INDELS index anew.
- **A `build-index-only --vcf` run stopped while it writes the dict-less tiers or the indel genotype
  store no longer leaves a part of them that searches use as the whole.** Once its enrichment is in
  place, `build-index-only` writes the dict-less tiers (`Dictionaries/registry_<vcf>/` and
  `genotypes_<vcf>/`) and the phased indel genotype store (`Dictionaries/indel_genotypes_<vcf>/`)
  there chromosome by chromosome, and a search reads them chromosome by chromosome. A run stopped
  in between (the tiers take about 150-180 GB of RAM a chromosome at scale, so an OOM kill is
  likely there) left the chromosomes written so far, and the next `complete-search` ran with them,
  exit status 0 and no warning, falling back to the dictionaries for the other chromosomes: there
  it reported no SNP+indel co-occurrence, and other columns differed. With the real engine, after a
  kill with one chromosome's tiers still to write, a toy site there that needs an insertion and an
  SNP was reported only by an alignment one base over, with 4 mismatches and 2 bulges, CFD 0.101
  and one of its two rsIDs, instead of 3 and 1, CFD 0.281 and both, and another site's best
  CRISPR-Bulge alignment came out from the reference, score 0.399, instead of from its variant,
  score 1.0 with its rsID and carriers. A run stopped once its enrichment was in place, before it
  began on them, left the same enrichment without them, which a search took for one built without
  `--samplesID`. Each of them is now marked unfinished (a hidden `.building` directory in it,
  holding the run's owner record) from before a fresh enrichment goes into place until it is
  written, and every search ends with an ERROR while the mark is there. `complete-search` says, for
  a run that has stopped, to run `build-index-only` again, which resumes them and removes the mark,
  and for a run still writing them, to wait for it. A search launched from the web runs the
  pipeline script itself, past those checks, and still used them as whole: it now fails with an
  ERROR in its `log_error.txt` that names the unfinished part and says to let that run finish or run
  `build-index-only` again to its end (for the tiers, with their `--samplesID`, or built again in
  the web Settings page) or, for the tiers and the indel genotype store, to remove them to search
  without them, and the results page shows it. The same goes for the enriched genome marked
  unfinished, and for the fake-indel genome marked unfinished (above) when the search would build
  its INDELS index from it. `complete-search` names each marked part by its real path, the one a
  symlink there points to. `publish-index` refuses such an index too. On success the files, their
  contents, stdout lines and exit status are as before. A chromosome whose tiers or indel
  genotypes fail to be written is still skipped with a WARNING, but the mark now stays, as a line
  after the step's summary says, since a search would take the part written for the whole: run
  `build-index-only` again once the cause is fixed, which writes the rest, or remove them to
  search without them. That run still builds the indexes, but now fails at its end with an ERROR
  (exit status 1) that says so, as does a run that ends with a part of the dataset that another
  run marked unfinished, such as the tiers that a run with `--samplesID` was stopped writing,
  which a run without it (as a build from the web Settings page was before this release) does not
  write: it exited 0, which the Settings page showed as done, while every search of the dataset
  was refused. A full disk or quota, an I/O error or no memory left there, which every chromosome
  after would hit too, ends the run with an ERROR, the mark kept: it was a WARNING, the step took
  its mark away, and once space was freed a search used the part as whole (with the real engine, a
  run that filled the disk in the indel genotype store failed later, at the variant index, and the
  next `complete-search` exited 0 with that store on one or two of its three chromosomes, a site
  losing one of its two rsIDs). So does a `.done` marker the disk could not take beside a whole
  indel genotype store, which was left empty and taken for done for good. Downloaded indexes carry
  no mark. A part that a run of an earlier release left has no mark and is not detected: if such a
  run was stopped while it emitted the tiers or the indel genotype store, run `build-index-only`
  again to its end before you search with that dataset.
- **The build directories of killed runs are reclaimed.** A SIGKILL (an OOM kill, a container stop)
  leaves the build directory, holding a partial index or enrichment as large as the finished one,
  which nothing read or removed. Each build directory now holds an owner record (host, PID namespace,
  PID and start time, as a search's build lock does), and each `build-index-only` run first removes
  those whose owner is gone on this host and PID namespace (dead, a zombie, or its PID now another
  process's), with a line on stdout for each, whatever it goes on to build: after a kill, the next run
  is often a `complete-search`, which builds the missing index or enriches the genome itself. The SNP
  overlay's hidden directories beside the fake-indel FASTAs of the run's dataset hold an owner record
  too, and are removed the same way. One whose owner is running or stopped, or has no record, is left
  alone, as is any directory `build-index-only` did not make. One made in another PID namespace or on
  another host, whose owner cannot be checked from here (a web container's own, once it restarted),
  is left too, with a line on stdout that says how to remove it. The enrichment's directory,
  `Genomes/run_*` before, which the web genome list showed as a genome, is now hidden and reclaimed
  the same way, as is the one in `Dictionaries` that its dictionaries are put together in; the
  search's own `Genomes/run_*` directories are untouched.
- **A full disk, inode table or file-count quota while `build-index-only` builds an index no longer
  leaves that index in place, cut short or short of contigs or of their parts.** `buildTST`
  (CRISPRitz 2.9.1 and earlier) does not check its writes, so on a full disk it left a `.bin` empty
  or cut short and still exited 0. `build-index-only` put that index in place (`Index built`,
  `Variant index ready`) and, when the rest of the run fitted, exited 0 with empty stderr; every
  later run printed `already present` and kept it, and every search failed on it (`is not a valid
  CRISPRitz index (.bin) file (numLeaves=65535)`, `invalid TST node count (numNodes=0)`, or a crash)
  until its directory was removed by hand. Before an index goes into place, each `.bin` the engine
  wrote is now checked against what it records of itself (its number of targets, of tree nodes, and
  the target its tree ends with), the INDELS index's chromosome by chromosome; one that is not whole
  ends the run with an ERROR that names it and says the disk may be full, nothing goes into place,
  and the next run, once there is room, builds it again. A whole index goes into place as before;
  the check reads each `.bin`'s targets, about half a second for every 5 million (one `.bin`), so a
  few minutes for a genome-wide index. An index that a run of an earlier release put in place cut
  short is not detected: if a search fails on an index with one of the messages above, remove its
  directory and run `build-index-only` again. Nor does `buildTST` check that a `.bin` opened: with
  the disk's inode table or a file-count quota full (as on many HPC file systems), it made no `.bin`
  for what it had left to save, the contigs after the one it was on and that contig's later parts
  (it saves a contig's targets 5 million a `.bin`, each hg38 autosome in several), wrote nothing for
  them and exited 0, every `.bin` there whole, so that index went into place and every search lacked
  the off-targets left out, with no error (with the real engine, a ten-contig genome built with a
  few inodes free went into place with 0 to 9 of its 10 `.bin`, and a search found a site that is on
  every contig on 3, 6 or 9 of them; a genome whose last contig takes two `.bin` went into place
  with the first one, and a search lost the site in the second). Once the engine has ended,
  `build-index-only` now reads each contig left without a `.bin` for a target by `buildTST`'s own
  rule (its PAM on either strand, the window inside the contig, no N in it outside the PAM), checks
  that each other contig has its parts from 1 to the last one's number, and, where that last part
  holds 5 million targets, as each part but the last does, counts the contig's targets by the same
  rule. A contig with no `.bin` that holds a target, one that lacks a part, or one that holds more
  targets than its parts ends the run with an ERROR that names it, nothing goes into place, and the
  next run, once there is room, builds it again. Only the FASTAs of contigs without a `.bin` or with
  a full last part are read whole, up to the first that fails. Where it cannot read such a contig as
  `buildTST` does (a byte that is no IUPAC base, CRLF line ends), it makes and removes a file in
  each directory of the index, and a disk or quota that cannot take it fails the build the same way;
  with room left, such a contig is not told from a whole one. A contig with no target gets no `.bin`
  and its index goes into place, as before, as does a whole index whose build left the inode table
  full. The indexes a search builds itself, on demand, are not checked, and CRISPRitz's own writes
  are still unchecked.
- **A full disk or a kill no longer leaves a samplesID list cut short for good.** A merged-panel
  `build-index-only` (two or more datasets) ends by copying each dataset's samplesID list to
  `samplesIDs/<genome>_<dataset>.samplesID.txt` and writing their union,
  `samplesIDs/<vcf>.samplesID.txt`, which `download --what index` and a search launched from the web
  also write when it is missing. Each was written in place, and every later run keeps one that is
  there: one that a full disk or a kill cut short was kept for good, and every later search of that
  dataset read it as the whole list, with no error, `build-index-only` exiting 0. Each is now
  written to a hidden file beside it and renamed into place, whole or not at all, and a full disk or
  quota, an I/O error or no memory left as `build-index-only` writes them ends the run with an
  ERROR, for the next run to write them; any other failure there is still a WARNING. On success the
  files and their contents are as before. A list that a run of an earlier release cut short is not
  detected: if one stopped or failed as it wrote them, remove the lists it wrote and run
  `build-index-only` again.

- **A source install read SNP dictionaries about 45 times slower than the Docker image.**
  `environment.yml` installed `ijson` without `yajl` and `cffi`, so ijson fell back to its
  pure-Python parser where the image uses its C backend (`yajl2_cffi`). The dictionary-based variant
  post-analysis (complete-test's benchmark mode, a search on an index built from your own VCFs) and
  `build-index-only`'s Tier-0 registry compile stream those dictionaries: one pass over chr22's
  1.55 GB dictionary took 189 s instead of 4 s, and the chr22 benchmark's SNP post-analysis about
  200 s instead of 5 to 17 s. `environment.yml` now lists `yajl` and `cffi` beside `ijson`, as the
  Dockerfile does; a linux-64 solve of it adds exactly `yajl`, `cffi` and `pycparser`. Results were
  identical either way, and a search on a shipped prebuilt index reads no such dictionary. A
  Bioconda install lacks the two packages too until its recipe adds them.

- **A search with no samples no longer writes into the install's empty `vuoto.txt`.**
  `complete-search` without `--samplesID` (a reference-only search, or a variant-aware one
  without it) hands the pipeline the shared empty placeholder `PostProcess/vuoto.txt` as its
  samples file. The pipeline added a newline to it as it started and removed it only at the
  very end, so every such run that ended earlier left one more newline in it for good: a
  search that found no off-targets, one that failed, or one stopped with Ctrl-C (measured
  through the CLI: 0 bytes, then 4 after one failed run and three with no hits). A
  personal-assembly search (`assembly-search`, which the web interface's Personal-assembly
  tab runs) does a reference-only search per haplotype with the placeholder, so it was
  affected the same way: 2 bytes after one that found nothing. Only a search from the web
  search form uses a samples file of its own, and never touched the placeholder. This
  happened on every install the run could write to: Docker (as root), source installs, and
  Bioconda, where the file is a hard link into the conda package cache, so the stray bytes
  also reached any env later created from that cache. A read-only install was never written.
  No result changed. The pipeline skips blank lines in a samples file. Where the same file
  stands in for a missing `--annotation`, it is recognised by its name; for a missing
  `--gene_annotation`, its copy goes through `sort-bed`, which drops blank lines. So a search
  run with stray newlines in the placeholder reports what it reports with the file empty.
  Reference-only and variant-aware searches without `--samplesID` report exactly what they
  reported before (every result file compared; the only differences are those between any
  two runs: the plot PDFs' creation dates, the logs' dates, timings, chromosome order and
  free memory, archive timestamps, and the report's dates and plot element ids). The
  pipeline now leaves the placeholder alone. It tells the placeholder apart as that very
  file, not by its name, so a `--samplesID` list of your own that is also called `vuoto.txt`
  is read as before. An install that already has the extra newlines needs nothing; to empty
  the file anyway, run `: > "$CONDA_PREFIX/opt/crisprme/PostProcess/vuoto.txt"`. Present in
  2.8.0 too. Covered by `PostProcess/test_index_resolution.py`.
- **`scorer-env check` no longer lists `pkg_resources=?`.** `pkg_resources` has no version of
  its own, so a healthy scorer env showed `pkg_resources=?`, which read like a missing
  package. It now shows the version of setuptools, which ships it (for example
  `pkg_resources=80.10.2 (setuptools)`; CRISPR-Bulge needs setuptools below 81 for it), or
  `ok` where no version can be read. The same value goes to `Annotations/.scorer_env.json`.
  The status, the other versions and the compute self-test are unchanged. Covered by
  `PostProcess/test_scorer_env.py`.
- **`validate-test` no longer prints pandas `FutureWarning`s, and runs under pandas 3.** Each
  benchmark it compared printed `FutureWarning: Series.__getitem__ treating keys as positions
  is deprecated` (two in the documented chr22 self-test): the brute-force reference's columns
  were read by position from rows labelled by column name. pandas 3 removes that fallback, and
  there `validate-test` stopped with `KeyError: 0`; the shipped environments keep pandas below
  3, so a released install saw only the warnings. The columns are now read by position
  explicitly. The comparison, its report and its exit code are unchanged. Covered by
  `PostProcess/test_validate_test.py`.
- **A source install no longer leaves an untracked `Annotations/` folder in the checkout.**
  `install_from_source.sh` runs `scorer-env create` and `scorer-env check` in the checkout,
  and they save the scorer env's state to `Annotations/.scorer_env.json` there. The
  checkout's `.gitignore` covered the other data folders a search creates in the folder it
  runs in (`Genomes/`, `Results/`, `Dictionaries/`, `VCFs/`, `PAMs/`, `samplesIDs/`), but not
  `Annotations/`, so `git status` showed it after every install. It is now ignored (the
  top-level one only: `test/data/Annotations/` holds tracked test files). Nothing to do; the
  file is a diagnostic record. Covered by `PostProcess/test_install_from_source.py`.
- **The unit tests pass inside the shipped image.** The build-lock tests in
  `test_index_resolution.py` listed processes with `ps`, which the image does not have (no
  procps): run there, the module ended with 44 errors and 6 failures. They now read `/proc`,
  as the pipeline's own lock code does (`ps` only where there is no `/proc`), and pass in the
  image as on CI. `test_crispr_bulge_fidelity.py` changed to a developer's folder as it was
  imported, so `python -m unittest test_crispr_bulge_fidelity` failed everywhere. It now skips,
  saying why, where TensorFlow or its dataset is absent, and in the cbulge env with the
  dataset it runs the fidelity check as before (exit 0 on a match); its docstring says how to
  fetch the one dataset file it needs. Searches are not affected.

- **The report of a reference-only search described a variant search.** A search with no variant
  dataset (no `--vcf`, or "Reference only" on the web form) got a shareable report that gave the
  reference index as its variant dataset (`Variant dataset(s): NGG_2_hg38`; `reference only` for a
  search with no bulges), printed the population-level analysis note with its advice to re-run with
  `--per-sample`, an allele-frequency basis "over the merged NGG_2_hg38 panel" and the criteria of a
  variant search, and, from the one variant registry a default install has, that panel's size and
  "The database contributes N SNPs and M indels, all searched." The MAF notes under the tables named
  the same "merged NGG_2_hg38 union panel", and `data/run_manifest.json`, the IND traceability
  record, carried the index name in `variant_datasets` and the registry's counts in `database`. Such
  a report now says `none (reference genome only)` in the header and in the inputs card, whose only
  other row is the search criteria against the reference genome; the MAF notes say that every MAF is
  blank because no variant dataset was searched; and the manifest has `variant_datasets` and
  `database` null. A search is taken as reference-only when its `.Params.txt` records
  `Ref_comp False`, which the command line and the web both write exactly then (for a folder
  without that key: when neither the results file name nor `Genome_idx` names a dataset), and never
  when a row is variant-created. Affected: every reference-only search, from the command line, from
  the web, and the per-haplotype searches of `assembly-search`, whose reports the combined report
  bundles. The text dates from 2.4.0; 2.7.10 to 2.8.0 built no report for these searches at all, so
  2.8.1 would have been the first release to show it again. To correct a report built earlier,
  rebuild it with `crisprme.py generate-report --result-dir Results/<job>`. Search results are
  unchanged, and the report of a variant search is byte-identical. Covered by
  `PostProcess/test_generate_report.py`.

- **A multi-guide report garbled its variant dataset, and a report could take another dataset's
  registry.** A search with more than one guide names its results
  `Multiple_spacers+<PAM>_<genome>+<dataset>_<mm>+<bulges>_integrated_results.tsv`, and the report
  split that name at the `_` inside `Multiple_spacers`, so the dataset read `NGG+hg38+hg38+1000G`
  instead of `1000G` in the header, the inputs card, the MAF notes and `run_manifest.json`. That
  label matched no variant registry, so the panel size and variant counts were right only while a
  single registry was installed; with more, the report dropped the SNP and indel counts and the
  manifest's `database`, and counted the panel over every installed samplesID file
  (`3,048 individuals (1000G n=2,548, hgdpx n=500)` for a search of a 2,548-sample 1000G panel).
  Separately, the registry was the first whose name held every dataset token, so a search of
  `hg38_hgdpx` installed beside a merged `hg38_1000G_hgdpx` was given the merged panel's size and
  counts (3,048 individuals and 31,130 SNPs instead of 500 and 4,624). The report now reads the
  multi-guide name as it reads a single-guide one, and takes the search's own registry,
  `Dictionaries/registry_<dataset folder>`, before any match by name. Affected: every multi-guide
  variant search (the label), and searches on an install with more than one registry (the counts).
  Present since 2.4.0. To correct a report built earlier, rebuild it with
  `crisprme.py generate-report --result-dir Results/<job>`. Search results are unchanged, and the
  report of a single-guide search is byte-identical unless it had taken another dataset's registry.
  Covered by `PostProcess/test_generate_report.py`.

- **A failed `complete-test` now ends with one line, which names an error log only where one was
  written.** It ended with `crisprme.py`'s own `OSError` traceback, which pointed to
  `Results/crisprme-test-out/log_error.txt`: a search refused before it starts writes no log there,
  and the benchmark mode writes one for each benchmark, in `Results/crisprme-test-out_<name>`. It
  now ends with one line naming each error log, not empty, that this run's searches wrote
  (`complete-test` refuses to start while a folder of an earlier run of its mode is there, and names
  only that mode's: `Results/crisprme-test-out_<name>` in the benchmark mode,
  `Results/crisprme-test-out` in the smoke), or else saying to see the messages above, and exits 1,
  as before. Covered by `PostProcess/test_complete_test.py`.
- **`targets-scoring --index-path` refuses a dict-less registry that `build-index-only` has not
  finished writing.** The carrier overlay read `Dictionaries/registry_<vcf>` as it found it, also
  while a `build-index-only` run was writing it chromosome by chromosome, or after one had stopped
  partway, and for each chromosome not written yet it reported 0 overlapping variants, with no
  warning. It now refuses a registry marked unfinished (see the dict-less tiers entry), as
  `complete-search` and `publish-index` do: it stops with a one-line `Error:` saying to let that run
  finish, or to run `build-index-only` again to its end with its `--samplesID`. A registry that a
  run of an earlier release left partway carries no mark and is not detected. A missing registry is
  now also reported in one line instead of a traceback. Output for a finished registry is
  unchanged. Covered by `PostProcess/test_targets_scoring.py`.

- **A personal-assembly search with bulges searched the maternal haplotype with the paternal
  haplotype's index, and a second individual with the first one's.** `assembly-search` runs one
  `complete-search` per haplotype, and `complete-search` names a genome, and the index it builds or
  reuses (`genome_library/<PAM>_<N>_<name>`), after the genome folder's name. A bundle assembly
  (fetched from HPRC, or imported in Settings) keeps each haplotype in
  `Assemblies/<individual>/<haplotype>/genome`, so every haplotype of every individual was searched
  as `genome`: the paternal search built `NRG_2_genome`, the maternal search logged `Reference index
  present (NRG_2_genome)` and searched the paternal sequence, and the reconciliation then reported
  each site as found on both haplotypes. Haplotype-private off-targets were missed without a
  warning. Two `--genome-paternal` / `--genome-maternal` folders with the same name collided the
  same way. Reproduced with the real CLI on a synthetic bundle (4 mismatches + 1 DNA + 1 RNA bulge,
  the web default): the two haplotypes' result files were byte-identical, and at chr1:35000 the
  maternal results showed the paternal 2-mismatch protospacer instead of the maternal 4-mismatch
  one; a second individual whose paternal haplotype carries that 4-mismatch sequence also got the
  2-mismatch one. Each haplotype is now searched through a link named after it and its folder,
  `haplotype_genomes/<key>/<individual>_<haplotype>_<hash>` in the working directory, where `<hash>`
  is 8 hex digits of the folder's path, so it builds and reuses its own index
  (`NRG_2_HG01255_paternal_5aeb0fc6`, `NRG_2_HG01255_maternal_067bdb25`); a genome folder that is
  not in a bundle is named after the folder plus the same hash. The hash keeps these names apart
  from every genome in `Genomes/`: a bundle search never takes the index of a
  `Genomes/<individual>_<haplotype>` folder of the same individual (the layout before 2.6.0, or a
  custom upload), a search of that folder never takes the bundle's, and two bundles whose folder
  names differ only in characters a file name cannot keep (`HG.002`, `HG_002`) get an index each.
  The per-haplotype result files carry that name too
  (`<guide>+NRG_HG01255_paternal_5aeb0fc6+__4+2_integrated_results.tsv` instead of
  `<guide>+NRG_genome+__4+2_integrated_results.tsv`). Each index a haplotype search builds records
  the genome folder it was built from and the names, sizes and modification times of that folder's
  FASTA files (`.genome_source`, which CRISPRitz does not read), and a later search of that
  haplotype rebuilds the index when they have changed: a FASTA file added, removed, resized or
  modified at another time, as when an individual is deleted and imported again with other files. A
  FASTA file replaced by other content of the same size and modification time is not noticed, nor is
  a change made after a search that was killed before it ended (its container stopped, say): the
  record is written when the search ends, and an index without one is used as it is, as an index of
  a legacy `Genomes/<name>` folder built by an earlier CRISPRme, or downloaded, is used as before. A
  haplotype result written by the colliding search is not reused by a later run with the same
  `--output`: it is searched again. With the fix the paternal results are unchanged byte for byte
  apart from the genome's name (in file names, `.Params.txt`, `.command_line.txt` and the
  haplotype's own report: its `report.html`, `data/Params.txt` and `data/run_manifest.json`, also
  bundled in the combined report), the maternal haplotype reports its own sites, and a legacy
  `Genomes/<name>` pair (the layout before 2.6.0) runs exactly as before. Each haplotype's index is
  the size of a reference index; Settings > Delete > Personal assembly removes the individual's
  haplotype indexes and links with its folder, and its confirmation names those indexes and counts
  their size. **Affected:** every `assembly-search` with
  bulges on a bundle assembly since 2.6.0, from the web or the CLI. The haplotype whose search built
  the shared index (usually the first paternal haplotype searched) has correct results; every other
  haplotype searched with that index does not. **What to do:** run those searches again, and delete
  the index `genome_library/<PAM>_<N>_genome` they left (e.g. `NRG_2_genome`: Settings > Delete >
  `Index: NRG_2_genome`). `assembly-search` no longer uses it, but a `complete-search` with bulges
  of a folder named `genome` still finds it by that name: such a search of one haplotype of a bundle
  is refused while it is there (see below), and a search of any other folder named `genome` takes
  it. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search without bulges of a genome folder outside `Genomes/` came back empty, or scanned
  another genome.** The 0-bulge search reads the genome's FASTA files directly, and read them from
  `Genomes/<name of the --genome folder>/` instead of the folder it was given. Without such a folder
  the search engine reported `ERROR! /DATA/Genomes/<name>/ is not a directory.` but exited 0, and
  the search ended with "The search completed successfully; the result set is empty". With one, that
  folder was scanned in place of the given one (a `--genome other/hg38` whose chr1 carries a
  2-mismatch site at chr1:35000 reported the 4-mismatch protospacer of `Genomes/hg38` there). On the
  command line this hit `complete-search --genome <folder outside Genomes/>` with `--bDNA 0 --bRNA
  0`, and the same search without bulge flags when the library has no index of that genome, since
  the bulges then default to 0 (`complete-search --genome /data/mygenome --guide ... --pam ... --mm
  4`). It also hit every personal-assembly search without bulges whose haplotypes are kept outside
  `Genomes/`, as every bundle's are (`Assemblies/<individual>/<haplotype>/genome`): both haplotypes
  came back empty, and the combined report then failed on the empty set (`TypeError`, see below), so
  `assembly-search` exited 1; `--bDNA` and `--bRNA` are optional and default to 0, so this was the
  documented minimal command (`assembly-search --assembly-individual <name> --guide ... --pam ...
  --mm ...`). The search now scans the folder it is given: a copy of `Genomes/hg38` outside
  `Genomes/` reports the same 9 sites as `Genomes/hg38`, and on the synthetic bundle each haplotype
  reports its 9 sites. A search of a genome in `Genomes/`, which includes every search from the
  web's Reference genome tab, scans the same folder as before; reference and variant search outputs
  were compared byte for byte, with and without bulges. A search of such a folder with bulges takes
  the index of its name, and is now refused when that index, or a folder of that name in `Genomes/`,
  is another folder's (next entry). **Affected:** 0-bulge `complete-search` runs of a genome folder
  outside `Genomes/` from the command line, and 0-bulge personal-assembly searches of haplotypes
  outside `Genomes/` (every bundle assembly) from the CLI or the web. **What to do:** run them
  again. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search with bulges of a genome folder outside `Genomes/` could search another folder's
  sequence.** `complete-search` takes a genome's index (`genome_library/<PAM>_<N>_<name>`) and its
  variant genomes (`Genomes/<name>+<dataset>`) by the genome folder's name alone, whatever folder
  they were built from. With bulges, `--genome /data/projB/ref` took the index a search of
  `/data/projA/ref` had built and reported projA's sites (at chr1:35000 of the synthetic bundle,
  projA's 2-mismatch protospacer instead of projB's own 4-mismatch one), and `--genome other/hg38`
  took the index of `Genomes/hg38` and reported its sites, as the same command without bulges did by
  scanning `Genomes/hg38` itself (previous entry). Given no `--bDNA`/`--bRNA`, that search takes 1 +
  1 bulges when `Genomes/hg38` is there, so the plain command was hit too. The other way round, a
  search of `Genomes/<name>` took the index that a search of a folder of that name outside
  `Genomes/` had built. `complete-search` now refuses a search with bulges, or with `--vcf`, of a
  folder outside `Genomes/` whose name a folder in `Genomes/` has, and of any genome folder whose
  name an index recorded from another folder carries, with a message naming the folder or index in
  the way. The message says how to run the search again: with its `--bDNA`/`--bRNA` (given none, the
  search of a renamed folder, or of one whose index was removed, may take no bulges and miss every
  site that needs one), and, under `--index-path`, without it, since the pipeline builds no index
  there (and so without `--prune-pam`, which needs `--index-path`). An index that
  `complete-search` builds of a folder outside `Genomes/`, and one that `build-index-only` builds
  (the entry on its record, below), records that folder and its FASTA files (`.genome_source`, as
  a personal-assembly haplotype's index does), and is built again when those files change (not
  under `--index-path`, which builds no index and takes it as it is). The record names a folder of
  the machine that built the index, so `publish-index` leaves it out of the archive and `download`
  drops it: a downloaded index records no folder. A search of a genome in `Genomes/` runs exactly
  as before and its index records nothing; an index that records nothing (built by an earlier
  CRISPRme, or downloaded) is taken as before, except by a search of one haplotype of a personal
  assembly (next entry), and so are the indexes of a library given with `--index-path` other than
  the working directory's `genome_library/`. The reference and variant search outputs of
  `Genomes/hg38` were compared byte for byte, with and without bulges. Two cases are not caught:
  two first searches of two such folders that start together (one of them searches the index the
  other builds), and a search from the web's Reference genome tab, which runs the pipeline without
  this check, of a `Genomes/<name>` added after a folder of that name outside `Genomes/` was
  searched with bulges. **Affected:** `complete-search` runs from the command line, with bulges or
  with `--vcf`, of a folder outside `Genomes/` named like a folder in `Genomes/` or like another
  folder searched with bulges in the same working directory, and searches, from the CLI or the
  web, of that `Genomes/` folder made after them. **What to do:** run them again, with a name of
  its own for each folder and with the bulges they ran with (`--bDNA 1 --bRNA 1` for one with
  bulges that gave none). An index built before 2.8.1 records no folder: where folders of one name
  outside `Genomes/` were searched with bulges, remove the indexes of that name (Settings >
  Delete > `Index: <PAM>_<N>_<name>`, or the folder `genome_library/<PAM>_<N>_<name>`), so that
  each search builds its own. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A search with bulges of one haplotype of a personal assembly could search another haplotype's
  sequence.** Searched on its own, `complete-search --genome
  Assemblies/<individual>/<haplotype>/genome` names the genome, and its index, after that folder,
  `genome` for every haplotype of every bundle (first entry), so with bulges it took any
  `genome_library/<PAM>_<N>_genome` there: the one an `assembly-search` with bulges had left, or one
  that a search like this one of another haplotype had built. Given no `--bDNA`/`--bRNA`, the search
  takes 1 + 1 bulges when such an index is there, so the plain command was hit too. Reproduced on
  the synthetic bundle after an `assembly-search` of it: the maternal haplotype logged `Reference
  index present (NRG_2_genome)` and reported the paternal 2-mismatch protospacer at chr1:35000
  instead of its own 4-mismatch one. Such a search is now refused before it starts while an index of
  that name that it would take records no genome folder, as one built before 2.8.1 does, or while
  one records another folder (as any search is: previous entry). The message names the index and
  says to search the individual with `assembly-search`, which searches each haplotype on an index of
  its own, or to remove the index and run the search again, which then builds its own from the
  haplotype's folder. For either it gives the search's `--bDNA`/`--bRNA`: `assembly-search` takes no
  bulges by default, and, given none, the search of the haplotype may take none once the index is
  removed (on the synthetic bundle it then reported 9 sites instead of 17). A search without bulges
  and without `--vcf` is not refused: it scans the haplotype's folder (second entry).
  `assembly-search` and searches of any other folder are not refused this way. **Affected:**
  `complete-search` runs with bulges, from the command line, of one haplotype folder of a bundle, in
  a data folder where another haplotype had been searched with bulges. **What to do:** run them
  again with `assembly-search`, or delete the index (e.g. `NRG_2_genome`) and run them again, giving
  the bulges they ran with (`--bDNA 1 --bRNA 1` for one that gave none; `assembly-search`'s default
  is 0). Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **A genome in `Genomes/` that is a link to a folder of another name failed in the post-analysis.**
  The pipeline resolved the link and named the genome after the link's target in one step and after
  the link in another, so a search of `Genomes/hg38lnk` -> `../shared/GRCh38_copy` stopped with
  `grep: .../crispritz_targets/GRCh38_copy_..._targets.txt: No such file or directory` and `ERROR:
  off-targets post-analysis (reference) failed` (exit 1), with or without bulges. The web runs the
  same pipeline on `Genomes/<name>`, so its searches of such a genome failed the same way. A link
  named like its target (`Genomes/hg38` -> `/shared/genomes/hg38`) was not affected. The genome now
  keeps the name of the folder the search is given (each personal-assembly haplotype is searched
  through such a link, named after the haplotype). The linked genome's results equal those of the
  same sequence in `Genomes/hg38`, apart from the genome's name in the file names. **Affected:**
  every reference search (no `--vcf`) of such a link, from the CLI or the web. **What to do:** run
  them again. A variant-aware search (`--vcf`) of such a link still fails, as before, at the
  enrichment (`mv: cannot stat '.../SNPs_genome/hg38lnk_enriched/'`, then `ERROR: Genome enrichment
  failed on <dataset>`): for it, name the link like its target folder, or put the folder itself in
  `Genomes/`.
- **A personal-assembly search stopped before its combined report when a haplotype had no
  off-target.** When the guide had no off-target on one haplotype (or none there that lifted to
  hg38), or none on either, `assembly-search` wrote the combined table and then failed with
  `TypeError: cannot concatenate object of type '<class 'numpy.float64'>'` (exit 1); on the web the
  job finished without a combined report. The combined table (`<job>_combined_hg38.tsv`) left out
  the columns of a haplotype without off-targets, and the report reads both haplotypes' columns. The
  table now always has both haplotypes' columns, blank for a haplotype without off-targets as on any
  row found on one haplotype only, and the combined report is built (with no off-target on either
  haplotype it says "No off-targets in this run."). A search with off-targets on both haplotypes
  writes the same table and report as before (compared byte for byte, apart from generation times
  and the random ids in the plots). **Affected:** such searches, from the CLI or the web, e.g.
  `--mm 0` with a guide that matches one haplotype only. **What to do:** run them again for the
  combined report. Covered by `PostProcess/test_assembly_reconcile.py`.
- **A search could run on the index of another genome whose name ends in the searched genome's
  name.** The pipeline picks an installed index by the genome's name,
  `genome_library/<PAM>_<N>_<name>`, but the pattern it matched also took
  `<PAM>_<N>_<prefix>_<name>`, the index of another genome: with `NRG_2_mini_hg38` installed, a
  search of `hg38` with bulges used it instead of building `hg38`'s own, and used it ahead of an
  installed `NRG_3_hg38` too, as the smaller bulge budget won; nothing in the output said so. The
  same happened between two individuals of the `Genomes/<individual>_<haplotype>` layout (before
  2.6.0) whose names end alike: reproduced with `assembly-search` of `X_WALK01` first, after which
  both of `WALK01`'s haplotypes logged `Reference index present (NRG_2_X_WALK01_...)` and reported
  `X_WALK01`'s sites. Only an index named for the genome itself is now taken. A library without two
  such names picks the same index as before, and reference and variant search outputs were compared
  byte for byte. **Affected:** searches with bulges, from the CLI or the web, of a genome `<name>`
  while `genome_library/` held an index of a genome named `<prefix>_<name>` (`NRG_2_mini_hg38`
  beside `hg38`; the legacy layout's `X_WALK01_*` beside `WALK01_*`). **What to do:** run them
  again; for `assembly-search`, with a new `--output`, since a completed haplotype result of the
  same `--output` is reused. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **Two searches that started together on a data folder without `genome_library/` could fail.** The
  pipeline makes `Results/`, `Dictionaries/`, `Genomes/` and `genome_library/` when they are
  missing, after checking for each; two runs that started at the same moment (two web jobs, or two
  `docker run`s) both found `genome_library/` missing, and the second to make it stopped with
  `mkdir: cannot create directory '/DATA/genome_library': File exists` and `ERROR: Genome enrichment
  failed on _` (exit 1). `crisprme.py` makes `Genomes/`, `Results/` and its other folders the same
  way, and could stop with `FileExistsError` there. A folder made meanwhile is now no error.
  **Affected:** searches started together on a data folder without those folders, such as the first
  ones before any index is downloaded or built there, from the CLI or the web. **What to do:** run
  the one that failed again. Covered by `PostProcess/test_assembly_haplotype_genomes.py`.
- **The personal-assembly results header, History and combined report named both haplotypes
  "genome".** The header and History showed `Paternal genome: genome | Maternal genome: genome`, and
  the combined report (`report.html` in `<job>_combined_report.zip`) the same in its Paternal genome
  and Maternal genome rows: the name of the folder every bundle haplotype is kept in, which
  `assembly-search` wrote into the job's `.Params.txt` when it finished (the web form's own label
  was overwritten). All three now show the bundle's assembly name (`assembly_name` in its
  `metadata.json`; for an HPRC fetch, the HPRC assembly name), or `<individual>_<haplotype>` when
  the bundle records none, from the CLI and the web alike; a legacy `Genomes/<name>` folder keeps
  its name. Jobs that already finished keep the label they were written with.
- **Submitting the Personal assembly tab with no individual selected listed six fields the form does
  not show.** The warning listed "paternal genome", "maternal genome", the two chain files and the
  two chromAlias files, which the Individual picker fills. It now says "select an individual on the
  Personal assembly tab" once; any other missing input, such as the guide, is still listed.
- **The Personal assembly tab ignored the "Maximum edits" slider.** The web form never passed it on,
  so an assembly search always allowed mismatches + DNA + RNA bulges in total, whatever the slider
  showed (with the slider at 3 and Advanced options closed, the results header said "Max total
  edits: 6"). The slider now applies as on the Reference genome tab: with Advanced options closed it
  caps the total edits (`assembly-search --max-total-edits`), with them open the per-type caps
  govern. A search left at the defaults (slider 6 = 4 mismatches + 1 + 1 bulges) runs the same
  command as before.

- **The web form sent a search to an earlier job that was not the same search.** A submission
  whose guides and `.Params.txt` matched an earlier job was not run: the earlier job was shown
  instead. `.Params.txt` does not record the analysis mode, the full-report option or which
  annotations are enabled (it names only the assembled `.active.<genome>.bed.gz`), so a
  Population-level request was sent to an earlier per-sample job (and the reverse), a
  full-report request to a slim one, and a search made after the annotations were changed in
  Settings to a job with the old annotations; the job name typed was dropped, and nothing said
  that the search had not run. In the walkthrough, a Population-level request named
  `walk_popmode` opened the earlier job `walk_var_ps_...`, whose `.search_mode` is
  `per-sample`. Each web job now records these, the CRISPRme version and the job name in
  `.reuse_key.json`, and an earlier job is reused only when they match too. The status page of a
  reused job says that the search was not run again, and that submitting it with another job
  name runs it again. An identical resubmission is still reused; a job made before 2.8.1 has no
  record and is never reused. Two smaller faults of the same code are gone: the notification
  email of a reused job linked to the job folder that the reuse deletes, and a submission that
  joined a queued identical search left a job folder that never ran, to which a later identical
  submission could be sent. **Affected:** web searches only; the command line never reuses a
  job. **What to do:** for a web job that may have been answered
  by an earlier one, check its `.search_mode` and report against what you asked for, or submit
  it again with another job name. Covered by `PostProcess/test_web_search_reuse.py`.
- **The "Top 1000 candidate off-targets" plots and the Personal Risk Card "Top 100" plots left out
  every site with 0 or 1 mismatches + bulges**, with nothing on the plot to say so: every
  1-mismatch-or-bulge off-target and every perfect match, the on-target and any variant-created
  perfect match among them, the strongest sites of all. In the walkthrough's variant search the CFD
  plot showed 9 of the 18 sites found (the four perfect matches and the five 1-edit sites, CFD
  0.077 to 1.0, were missing), and sample S1's card plot 1 of its 3 personal targets. They are
  plotted now, the perfect matches marked with a black star and a legend entry, and the caption
  above the plot on the Graphical Reports tab says so. A plot of a single site without a variant is
  drawn too: it failed, and no later plot was drawn, so a search whose only off-target with 2 or
  more mismatches + bulges was a reference one got no top-1000 plots (the error was discarded). The
  card's "Top 100 Private Targets" plot was drawn from the personal targets, so it repeated the
  personal plot; it is drawn from the private ones now (an empty plot for a sample with none).
  **What changes:** the plot images a search writes to its `imgs/` folder
  (`CRISPRme_{CFD,CRISPR_BULGE,fewest}_top_1000_log_for_main_text_<guide>` and
  `CRISPRme_{CFD,CRISPR_BULGE}_top_1000_by_variant_effect_<guide>`, `.png` and `.pdf`), from the
  web and from `complete-search`, as shown on the Graphical Reports tab, on the personal-assembly
  results page and in the combined personal-assembly report; and the Personal Risk Card plots, on
  the web and from `generate-personal-card`. Search results and every table are unchanged; the
  shareable report draws its own plots and tables, and still leaves these sites out of them (next
  entry). Plots made before this release stay as they are, and that caption, shown above them too,
  says that a plot drawn before 2.8.1 leaves these sites out: run the search again for new ones (a
  Risk Card made earlier shows its saved plots). Covered by `PostProcess/test_top_plots.py` and
  `PostProcess/test_web_results_links.py`.
- **The shareable report did not say that its top-1000 plots and tables leave out the sites with
  0 or 1 mismatches + bulges.** Its Section 2 plots and its Section 6 tables, with their
  downloads `top1000.tsv` and `top1000_crispr_bulge.tsv` (and the Section 5 tables of the
  personal-assembly report, with `<name>_combined_top1000.tsv`), leave out the perfect matches
  and the off-targets with 1 mismatch or bulge, as before, with nothing on the page to say so. In
  the walkthrough's per-sample search that is 9 of its 18 sites: the four perfect matches and
  the five 1-edit off-targets, CFD 0.077 to 0.636. Each of these views now has a caption that
  says how many sites it leaves out and where they are: all in the recommended validation panel
  (`panel_top100.tsv`), the perfect matches first, and the 1-edit off-targets also in
  `mmb_le_1.tsv`. The tables and files themselves are unchanged. **Affected:** every report of a
  search that found a perfect match or a 1-edit site, from the web, `complete-search`,
  `generate-report` and `assembly-search`. **What to do:** take these sites from the validation
  panel or `mmb_le_1.tsv`, not from the top-1000 files; for a report made before,
  `crisprme.py generate-report --result-dir Results/<job>` writes one with the captions. Covered
  by `PostProcess/test_top_plots.py`.
- **The search form's bulge note and defaults did not match the installed indexes, and the docs
  said the web interface never builds an index.** With the reference index installed and 1 bulge
  selected, the note under Advanced options said that the first bulge search would build the
  reference index; for a variant panel whose genome and VCF files are installed but that has no
  index yet, it said there was "no raw genome to build one", although Settings can build it. That
  panel set the DNA and RNA bulges to 0, and with them Max edits to 4, and choosing Reference only
  again, or a panel with a 1-bulge index, did not give back 1, 1 and 6: a search submitted then ran
  with at most 4 mismatches + bulges, unless the user noticed (in the walkthrough data, a search on
  the panel `walk`, chosen after `walk2`, which has no index, found 15 of the 18 sites it finds
  with `walk` chosen straight away). The note now compares the bulges selected with the installed
  reference index and says when the search will build a deeper one first (a web search does, from
  the installed genome, once). For a panel with no variant index it says that only a 0-bulge search
  is offered and that the first one prepares the variant genome and a 0-bulge index (hours on a
  full genome), and to build the index in Settings ("Build an index locally") when the genome and
  VCF files are installed, or to download one otherwise. A search prepares the variant genome
  from VCF files one per chromosome, and only a build of the index splits a single
  multi-chromosome file, as the Settings card says: for a panel added as one such file (an upload
  or a server folder) the note now says to build its index in Settings first, or to add the panel
  again as one VCF file per chromosome, and Submit refuses its search with the same message. That
  search was accepted and failed within seconds (`ERROR: cannot open .../Genomes/hg38/.fa`,
  `Genome enrichment failed`), the status page saying only that the result encountered some
  errors and to submit it again, which failed the same way (2.8.0 too). The DNA and RNA bulges
  and the Max edits that the cap lowered come back as far as the cap allows them again, unless they were changed by
  hand in between: 1, 1 and 6 by default, on Reference only and on a panel with a 1-bulge index
  alike, also after several panels in a row. The web user guide, the Docker quickstart, the website
  and the data setup guide now say what a web search builds (a missing reference index; for a panel
  without a variant index, a 0-bulge one, from its VCF files one per chromosome), as does the
  message of the submit check. Covered by
  `PostProcess/test_web_bulge_notes.py` and `test/web/test_settings_vcf_datasets.py`.
- **A results link to a job without results gave a server error and a blank page.**
  `/result?job=<id>` for a job that does not exist, failed, was refused (an unfinished index),
  found no off-target or has not finished failed on the missing results file (HTTP 500 and only the
  navigation bar), and so did its Show Targets, sample and region pages. So did a job that failed
  at its last steps (results integration, plots or database) or was still in them: it already has
  its integrated results file, but not the other files these pages read. A link saved or emailed
  earlier could land there. These pages now say which it is, in the words of the job's status page,
  with a link to that page (a finished job whose results files were removed from its folder says
  so); the pages of a finished job are unchanged. Covered by
  `PostProcess/test_web_results_links.py` and `test/web/test_web_e2e.py`.
- **Results-page downloads that never became ready.** The Summary by Sample tab looked for
  `<job>.summary_by_samples.<guide>.txt`, but a search writes one file per criterion
  (`..._CFD.txt`, `..._CRISPR-Bulge.txt`, `..._fewest.txt`), so after 30 s it said "Download link
  unavailable"; it now links the file of the criterion selected, the one its table shows. The
  Show Targets page of Summary by Mismatches/Bulges looked for
  `Results/<job>/<job>/.<type>.<mm>.<size>.<guide>.targets.zip`, a path nothing writes, and said
  "Generating download link, Please wait..." for as long as it was open; it now links
  `Results/<job>/<job>.<type>.<mm>.<size>.<guide>.targets.zip`. Both Show Targets pages (a
  guide's and a sample's) stop asking after 30 polls (2.5 minutes) and then say that the link is
  unavailable, as Summary by Sample does. A job name containing `tsv` made the Show Targets zip
  go to a folder that does not exist (`tsv` was replaced by `zip` in the whole path); only the
  extension is replaced now. As before, a Show Targets download holds the rows of the table page
  shown. Covered by `PostProcess/test_web_results_links.py`.
- **The job status page held its View Results button twice** (two elements with the id
  `view-results`, as on the personal-assembly status page); once now, where it was shown.
  Covered by `PostProcess/test_web_results_links.py`.

- **A VCF dataset added in the web Settings page can be searched from the browser, per sample too.**
  A web variant search reads the dataset's sample list, `samplesIDs/<dataset>.samplesID.txt`, only
  after the search and post-analysis. Nothing in Settings wrote one for a dataset added there, and
  the Settings build passed no `--samplesID`, so a dataset added and built only in Settings ran its
  whole search and then failed (`ERROR: samples IDs processing failed on dataset ...`), was
  labelled `[sites-only (allele frequencies)]` on the search form, and could not run per-sample.
  Adding a dataset in Settings (an upload, or a registered server folder), or building its variant
  index there, now writes that list from the VCF headers when there is none, in the 4-column layout
  of the shipped lists: exactly the samples the VCFs genotype (a command-line search refuses a list
  that names a sample no VCF genotypes, #46, or that leaves out one they genotype), each with the
  population, superpopulation and sex an installed sample list gives it (a downloaded panel's, or
  yours; of the lists that give it a superpopulation, the last by file name), and with the dataset
  as its population and `unknown` superpopulation and sex when no list gives it a superpopulation
  (the shareable report reads every list for a sample's superpopulation, and a list that gave
  `unknown` to a panel's samples turned them unknown in the reports of the panel's searches).
  The list is written again when a file added to the dataset, or a build of its index, gives another
  one (a file with other samples, a sample list installed or removed since), and it goes when the
  dataset is removed in Settings. A list from anywhere else (one you wrote or edited, or a
  download's) is never changed: the page tells its own by the SHA-256 it keeps of it in the hidden
  `samplesIDs/.<dataset>.samplesID.sha256`. The Settings build passes the list, as a hidden one-line
  listing, `samplesIDs/.<dataset>.samplesID.list`, so it also builds the dict-less tiers (which
  searches of its index then read: see Changed): the dataset reads `[genotyped, phased]`
  (`[genotyped, unphased]` for an unphased VCF) and per-sample analysis is available. It passes a
  list only when it names samples the VCFs genotype: a sites-only dataset (no sample columns, as
  the `hg38_mega` VCFs) is built without one, as before, so the per-dataset allele-frequency
  registry put in place for it is left as it is. Removing the dataset in Settings removes the two
  hidden files too. Checked in the browser on a phased 6-sample dataset uploaded and built only in
  Settings, once per chromosome and once as a single multi-chromosome file named without the
  genome: the per-sample search runs to its end and names the same carriers, with the same allele
  frequencies, as the same panel built from the command line with `--samplesID`. Submit now refuses
  a variant search whose sample list is missing, with a message that says how to get one, rather
  than letting it fail at its end (for a sites-only dataset, which has no samples to list, it says
  so, and that a search reports such a VCF's variants only through a per-dataset allele-frequency
  registry, as the shipped `hg38_mega` index has; the note under the bulge options says the same of
  such a dataset, an upload or the `hg38_mega` VCFs a Settings download installs, with or without
  its index built in Settings, where 2.8.0's offered it a 0-bulge search and a build in Settings, or
  the bulges of the index built; the web user guide says so, and that the published
  `NRG_3_hg38+hg38_mega` index is the way to search the mega variants); and a genotyped dataset with no tiers yet reads
  `[genotyped, samples not indexed]`, not sites-only.
  **If you added a dataset in Settings,** build its variant index again there; to give its
  samples' populations and sex (the tiers count a sample of unknown sex as diploid on chrX and
  chrY), put your own list at `samplesIDs/<dataset>.samplesID.txt` first. A list put there, or
  written there again, once the variant index is built reaches its tiers only after you remove
  `Dictionaries/registry_<dataset>` and `Dictionaries/genotypes_<dataset>` from your data folder
  and build the index again: until then a build leaves the tiers as they are, from the list they
  were built with. A search takes its samples' populations and superpopulations from the list that
  is there, and a SNP's allele frequency from the tiers, so on chrX and chrY that frequency counts
  each sample by the sex the tiers were built with.
- **The Settings VCF upload records the reference genome you pick, and the search form offers the
  dataset for it.** The chunked uploader looked the genome up with selectors Dash 2's dropdown does
  not have and sent none, so `VCFs/.<dataset>.refgenome` was never written; and the search form
  offered only datasets whose name holds the genome's, so one named like the card's own example
  (`my_pig_vcf`) never appeared, even once its variant index was built. The uploader now sends the
  dropdown's value; an upload into a dataset with no recorded genome is refused until one is picked
  (the card already said it is required); and the search form offers a dataset for the genome its
  marker names, with the bulges its index supports. A dataset with no marker is offered by its name,
  as before, and now also for the genome of an enriched genome built from it
  (`Genomes/<genome>+<dataset>`).
- **A genome uploaded in Settings is laid out as every search needs it.** A FASTA of several
  sequences was stored as it was, so every search on it was refused with a command-line-only fix,
  and a UCSC `chromFa.tar.gz` kept its `chroms/` folder and an archive of `.fasta` files their
  names, so their index held no `.bin`; all were listed as installed. An upload is now laid out as
  every search needs it: every FASTA of an archive goes to the top of the genome folder as a `.fa`
  file (`chr1.fasta` as `chr1.fa`: an index is built from the `.fa` files only), a FASTA of several
  sequences is split into one file per sequence, and a FASTA of one sequence is named after it: a
  one-chromosome subset (`hg38_chr22.fa` holding `>chr22`) or a bacterial genome went in as
  `<genome>.fa`, and as a search takes each file for the sequence it is named after, every search of
  it found its off-targets and then dropped them all (`none survived post-analysis`), though its
  index build said done. (A genome downloaded by Direct URL or from UCSC as one FASTA of one
  sequence is named after it the same way, and `PostProcess/genome_layout.py --split` renames such a
  file too.) Each header is written as `>id` (see the entry on those). Hidden files of an archive
  are left out: the `._chr1.fa` that macOS `tar` adds beside `chr1.fa` (its file attributes) was
  installed too, and every search read it as a chromosome and failed, though the build said done. An
  archive with two FASTA files that would go in under one name (`chr1.fa` and `chr1.fasta`, or two
  files of one sequence) is refused. (An upload under the name of an installed genome now replaces
  it: see Changed.)
  **A genome uploaded before** this fix that holds a multi-sequence FASTA, a `chroms/` folder,
  `.fasta` files, `._` files or a FASTA not named after its one sequence: remove it in Settings with
  the indexes built on it (*Remove installed data*), and any `Genomes/<genome>+<dataset>` variant
  genome from your data folder, then upload it again. (An index built from a multi-sequence FASTA
  holds the sites of every sequence under the name of the first: kept, it gave a chr2 off-target as
  `chr1:37386` instead of `chr2:17381`.)
- **`build-index-only` fails on a genome folder it cannot index, and never puts in place an index
  with no `.bin`.** `index-genome` indexes the `.fa` files at the top of the genome folder only and
  exits 0 with an empty index when there are none (a `chroms/` folder, `.fasta` names, gzipped
  FASTAs); `build-index-only` put that index in place, and the Settings page showed the build as
  done, though every search refuses an index with no `.bin`. It now ends with an ERROR (exit status
  1) that says what to do: before it builds anything, for a folder with no `.fa` file at its top,
  with a FASTA of several sequences (which every search refuses too; the message says to remove the
  indexes built from the folder while it held that file), with a FASTA of one sequence not named
  after it (`walkone.fa` holding `>chr1`, whose off-targets every search found and then dropped; the
  message gives the `mv`, or `genome_layout.py --split`, which renames them all) or with a header of
  more than `>id` (see the entry on those); and, putting nothing in place, when the engine writes no
  `.bin` for the reference or the variant genome. A reference index with no `.bin` that an earlier
  run left is no longer taken for present (the run said `Reference index already present` and the
  Settings page showed the build done, though no search uses it): the run says to remove it. A
  genome of one FASTA per sequence, each named after its sequence and with a header of `>id` (see
  the entry on those), builds exactly as before.
- **A single multi-chromosome `.vcf.gz` without its `.tbi` builds a variant index.** The
  documented auto-split read the file by region through its index, so a file uploaded alone (the
  Settings uploader takes one file at a time) failed with `Could not load .tbi index ...
  Auto-splitting the VCF dataset into per-chromosome files failed`. A file with no `.tbi` or `.csi`
  is now indexed first, in a scratch folder (a registered server folder is never written to); a
  plain-gzip file, which tabix cannot index, is first recompressed there in bgzip format. What
  tabix and bcftools print goes to stdout (a warning of theirs on stderr fails a pipeline stage),
  and an unsorted file fails with their reason and the `bcftools sort` to run. The split folder now
  gets the same permissions as any other folder made there, not the `0700` of the temporary folder
  it came from. Per-chromosome datasets are not touched, as before.
- **No false "Samples degraded ... Re-download this index" alert after a Settings variant build.**
  The alert (new in this release) is meant for an index download whose genotype store did not
  install, but it matched any `genotype store` text, so the line every variant build prints
  (`[indel-snp] phased indel genotype store: ...`), and a publish's, raised it: wrong advice, as a
  local build has nothing to re-download. It now matches the download's own notes only, and is
  checked for `Download index` jobs only.
- **The Settings page says when a stopped build left its build directory behind.** A `docker stop`
  during a Settings build kills it before it can clean up, and each later container, with a new host
  name, cannot tell whether that build still runs, so it leaves the hidden directory (which can be
  as large as the finished index or enrichment), with a line on stdout that the page showed at most
  in the full log of a failed build. When a Settings job ends, done or failed, the page now names
  each such directory still there, with its size, and says to delete it once no other CRISPRme
  container is building. The Docker quickstart and `docs/TROUBLESHOOTING.md` give the command.
- **An annotation BED whose file name names no genome can be enabled; a second upload of the same
  bytes is added.** The upload said "Enable it under 'Manage annotations'", but that card offered
  only BEDs whose file name holds the genome's, so `walk_regions.bed` was offered for no genome. A
  BED whose file name names no installed genome is now offered for every genome, to enable for the
  one its coordinates are on; one that names an installed genome is still offered for that genome
  only, and the card says so. The upload widget is emptied after each file, so a next file with the
  same bytes (a copy under another name) is added, where nothing happened.
- **Settings: the "Started — progress is shown in the Activity panel" note is no longer shown as an
  error, and goes when its job ends;** the VCF card no longer lists "a URL" as a source it does not
  offer, and says that a single multi-chromosome file's `.tbi` is optional and where the dataset's
  sample list is; the web interface guide no longer says that custom VCFs and PAMs need the command
  line and a server restart, nor that an annotation upload takes a `.bed.gz`.
- **The shareable report's "By superpopulation" chart keeps a sample's superpopulation when a
  sample list read after it gives none.** The report reads every list of `samplesIDs/` and took a
  sample's superpopulation from the last list (by file name) that names it, so a list that gives
  it as `unknown` or `NA` (one you wrote for your own dataset of some of a panel's samples, or one
  the Settings page wrote before a panel naming them was installed) turned those samples unknown in
  the reports of the panel's searches. A list that gives no superpopulation no longer replaces one
  another list gives; where lists give two, the last one's is kept, as before.
- **`Other_motifs` and the alternative alignments count and list each reference alignment once.**
  The merge step of a search, which keeps each site's best alignment and lists the others as
  alternative alignments, took a site's reference rows as they came, and two steps before it give it
  the same alignment more than once: each variant dataset's post-analysis carries every row of the
  one reference search, and on the dict-less tiers (searches of a `--samplesID` build, of a
  downloaded index such as the default one, and of a variant index built in the web Settings page)
  the post-analysis also gives the reference alignment of a variant window on the plus strand again.
  So a search that read the tiers counted each reference alignment over a SNP twice in
  `Other_motifs` (`444` where the dictionaries gave `222`), a search of two variant datasets counted
  every reference alignment once per dataset (`491` where one dataset gave `245`), and a per-sample
  search's file of alternative alignments listed each of them as many times, a copy of the site's
  best alignment among them. On chr22 of the default index, a per-sample search gave 6152 rows for
  its 4941 alternative alignments, and 280 of its 1940 sites a higher `Other_motifs` (280 of 2356 at
  population level); on a 6-sample panel searched with a second dataset, 3779 rows for 1887. Each
  reference alignment is now kept once. Rows equal in every column but `#Seq_in_cluster` and the two
  score columns are one alignment: each dataset scores the rows itself, and CRISPR-Bulge can give
  one alignment two values (`4.3390600000000004e-23` and `4.33904e-23`), of which the higher stays.
  What changes is `Other_motifs` in the integrated table and on the results page, the
  `#Seq_in_cluster` columns of the `bestMerge` table, the file of alternative alignments and the
  `observed_haplotypes_raw` file made from it. The best alignment of each site, its scores, carriers
  and allele frequencies are as before, `observed_haplotypes.tsv` and the report do not change, and
  a search of one dataset on its dictionaries, or of the reference alone, is unchanged (a search
  with a sites-only dataset changes further: see the entry on those). Run a search again for these
  counts and lists.
- **A search with a sites-only dataset lists each site once and loses none; a VCF with no `AF` keeps
  a variant's best alignment.** The merge step of a search sorts each criterion's rows by guide,
  chromosome, position and score before it makes one row of each site, and its sort counted columns
  by blanks, so an empty column moved the later ones up. A row of a sites-only dataset (`hg38_mega`,
  alone or beside another dataset, since 2.5.2) has an empty `Samples`, so it went by its rsID and
  not its guide: such rows came before or after all the others, and the merge, which starts a new
  site only when the position rises, took unrelated sites for one. A site was listed more than once;
  a variant row was marked `Not_found_in_REF` (a `VARIANT` off-target in the report) though its site
  has a reference row; a site could lose its row altogether; and in a per-sample search a row of
  alternative alignments could pair the three criteria's alignments across sites, so
  `observed_haplotypes.tsv` gave a site another site's variants. The sorts now read the
  tab-separated columns, and a variant both datasets give is one row, with the genotyped dataset's
  carriers and both allele frequencies. On chr22 (guide `CTAACAGTTGCTTTTATCACNNN`, 4 mismatches, 1+1
  bulges), `hg38_mega` alone goes from 1621 rows to 1384, rows within 3 bp of another from 556 to 0,
  `Not_found_in_REF` rows from 659 to 422, and it lists 50 sites it did not (the best at
  chr22:50091459, CFD 0.022); the report's `REFERENCE` and `VARIANT` off-targets go from 766 and 659
  to 962 and 422, and its CFD >= 0.05 net from 53 sites to 46. The default index with `hg38_mega`
  goes at population level from 3015 rows to 2520 (within 3 bp 997 to 4, `Not_found_in_REF` 2053 to
  1558, 21 sites more); per sample, with the next entry's fix, from 2406 rows to 2137 (671 to 2,
  1444 to 1175, 73 sites more), and its rows of alternative alignments with the fewest or the
  CRISPR-Bulge block more than 100 bp from the CFD block go from 28 and 27 to 0, and its haplotypes
  with a variant more than 1 kb from their site from 26 to 0. A row with only its `AF` empty (a VCF
  with no `AF`, which the input check warns about: every variant on the dictionaries, every indel on
  the tiers; a record that has `AF` at another place in `INFO` than the VCF's first record, by the
  known limit on the place of `AF` under Changed: its indel on any index, as are 91,950 of the
  102,665 indels of HPRC's chrY, and its SNP on the dictionaries) stayed at its site but went by its
  REF score, so the merge could keep a lower-scoring alignment of a variant: on a 6-sample panel
  with no `AF`, the perfect match a deletion makes at chr2:10000 read CFD 0.9 with a bulge and now
  reads 1.0 with none, and the report's perfect matches go from 3 to 4. A search whose rows have no
  empty column (the default index alone, a VCF with `AF` at one place in `INFO` in every record, the
  reference alone) sorts as before. **If you searched a sites-only dataset, a VCF with no `AF` or
  one with `AF` at more than one place in `INFO` (such as HPRC's chrY, or a VCF merged with
  `merge_vcf_panels.sh`),** run the search again.
- **A sites-only dataset goes the population-level way in a per-sample search too, as the command
  line said it did, and loses no multi-variant worst case there.** A sites-only dataset
  (`hg38_mega`) has no genotypes for a per-sample search to resolve, and `complete-search` said that
  `the population-level analysis runs instead`. It did not: the dataset's rows went through the
  per-sample path, which joins two variants only through carriers they share, and it has none, so a
  site that takes two variants together lost its row, and others scored lower, with no warning.
  Installing the default index and `hg38_mega` lists both for a search
  (`docs/PRECOMPUTED_INDEXES.md`), and the two then ran that way with no mode flag. Now each dataset
  goes by its registry (its `data_type`, or for a registry that has none its aggregation): a
  sites-only dataset's post-analysis is the population-level one whatever the search's mode.
  `--per-sample` on sites-only datasets alone is the population-level analysis (the same files as
  with no flag), with a WARNING on stdout. Beside a genotyped dataset, whose rows stay per-sample,
  `complete-search` says which goes which way. With no mode flag a phased dataset beside a
  sites-only one now runs per-sample whichever name sorts first; with the sites-only name first it
  ran population-level. On chr22 (guide `CTAACAGTTGCTTTTATCACNNN`, 4 mismatches, 1+1 bulges),
  `hg38_mega --per-sample` goes from 1318 rows to 1384, the table of `hg38_mega` with no flag,
  byte for byte, and lists 63 sites it did not (with the TRAC guide `CTCTCAGCTGGTACACGGCA` beside
  it, at 4 mismatches + 1 DNA + 1 RNA bulge, from 3,605 rows to the 3,801 of the run with no
  flag, 3,802 with this release's indel-row fix as well; the other counts in this entry are
  from before that fix); the default index with `hg38_mega`, per-sample as it runs with no
  flag, goes from 2077 rows to 2137 (`Not_found_in_REF` 1115 to 1175, alternative alignments
  5863 to 6226, `observed_haplotypes.tsv` 1556 rows to 1796). It lists 57 sites it did not and loses none, 14
  sites score higher and none lower (chr22:34270481 from CFD 0.173 to 0.233, and chr22:22881199
  reads 0.274 where chr22:22881195 read 0.262), and the report's CFD >= 0.05 and >= 0.2 tables
  each gain one site. On a denser search of the two (chr19 and chr21, 5 mismatches, 2 DNA and 2
  RNA bulges, guides CTAACAGTTGCTTTTATCAC and CTCTCAGCTGGTACACGGCA), 249 sites of `hg38_mega`
  searched alone were not in its table (the best chr19:2433201, CFD 0.533) and 411 scored lower
  there (chr19:388343, 0.503, read 0.114). Now each site that either dataset gives alone is in the
  table with its CFD or a higher one, but for 29 with a CFD of at most 0.002 that the merge step
  joins, as before, to a site within 3 bp on the other strand that scores at least as high. The
  table goes from 49,193 rows to 49,452, and its sites with a CFD >= 0.05, >= 0.2 and >= 0.5 from
  1,273, 138 and 5 to 1,316, 151 and 7. A sites-only dataset's rows have no carriers and read
  `putative` in the report's Observed column (`observed` for one variant). A search with no
  sites-only dataset, `hg38_mega` with no flag and any search at `--population-level` are
  unchanged, byte for byte. **If you ran a per-sample search with a sites-only dataset,** run it
  again. Covered by `PostProcess/test_sites_only_per_sample.py`.
- **A command-line search of two datasets is no longer refused because only one has its VCFs on
  disk.** The input check of `complete-search` compared the samples of all the `--samplesID` lists
  with the VCFs on disk (#46). With `hg38_walk`'s VCFs there and `hg38_walkB` installed as an index,
  it stopped with `--samplesID lists 3 sample(s) genotyped in no VCF (B4, B5, B6)`, and told you to
  filter a list that was right; so too beside a sites-only dataset, on disk or not, whose VCFs have
  no sample columns and whose list names `1000G2021` and `HGDP`. Each step that leads there is
  documented (`download --what vcf` for one dataset, several indexes in one list, your own cohort
  beside a dataset). Each list is now compared on its own: up to as many lists as the search names
  such datasets, those that name no sample of the VCFs on disk, are left out. A list that names some
  of them and others they lack still stops the search, with those others named; so does a list of no
  dataset when every dataset's VCFs are there. Measured on a 6-sample panel, the search that was
  refused gives the same table as with both datasets' VCFs on disk, and as with neither. A web
  search does not run this check and is not affected. Covered by
  `PostProcess/test_samplesid_overlisting.py`.
- **A variant build that runs out of memory as it writes its tiers ends with a message, and no
  longer hangs.** `build-index-only --samplesID`, and in this release every variant build in the web
  Settings page, writes the dict-less tiers and the phased indel genotype store in forked workers,
  one chromosome each. The pool they ran in did not see a worker die: when the OOM killer took one,
  the run waited for good, at 0% CPU, with nothing on stdout or stderr. In Settings the build read
  `running` for ever, with its buttons locked and every later job queued behind it, until the
  container was restarted; the same build finished in 2.8.0, which wrote no tiers. A worker that
  dies now sends the chromosomes left to one worker at a time, with a WARNING on stdout, and the
  build goes on if they fit. If a worker dies then too, the run ends with exit status 1 and an ERROR
  that says memory is the likely cause, names `CRISPRME_TIER_WORKERS` (`CRISPRME_INDELGT_WORKERS`
  for the store), and says how to go without what it could not write; Settings shows the build as
  failed with that message. A run with one worker, or of a dataset of one chromosome, which did not
  hang but died at exit status 137 with no message, now ends with that ERROR too. What it could not
  write stays marked unfinished, so a search refuses it, and the next run writes the rest; no
  variant index is built by that run. To build with no tiers, remove them as the message says and
  run `build-index-only` from the command line with no `--samplesID`: the Settings page has no way
  to. No worker outlives a run that is killed outright either, one that waits for its next job looks
  for the run every few seconds. Measured on a synthetic panel of 3000 samples and a SNP every 3 bp
  on 3 x 50 kb: with a limit of 1.5 GB the Settings build reads `failed` with the message after 5 s,
  and with 6 GB, where three workers at once are too many, the build finishes one chromosome at a
  time. A build with memory enough writes the same files as before, byte for byte, and so does that
  one.
- **A genome whose FASTA headers carry a description, or whose lines end in CRLF, no longer indexes
  as done and then loses every search's results.** CRISPRitz names a chromosome after the whole
  header line, and the pipeline splits on blanks, so the FASTAs NCBI, GENCODE, Ensembl and GATK ship
  went wrong, each in its own way, with the same sequence that gives 7 sites under a plain `>chr1`.
  An NCBI `>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome` or a GATK
  header: exit status 0, `No off-targets found` (in the browser, `Finished: no off-targets were
  found`). GENCODE `>chr1 1`: 7 sites, but with `chr1 1` for a chromosome and `NA` for every
  annotation. Ensembl `>1 dna:chromosome ...`: the search fails at integration (`Start coordinate is
  too large ... BEDOPS`). CRLF: 451 KB of warnings on stderr at exit status 0, and `No off-targets
  found`. A variant build on `>chr1 1` failed with `Indels indexing failed`. An upload in Settings,
  a Direct URL or UCSC download and `genome_layout.py --split` now write each header as `>id`, and
  every line with no CR, also for a FASTA already named after its sequence; a file that needs none
  of this is not touched. Each of those genomes then gives the sites and table of the same sequence
  under a plain header, through the browser and the command line, and a variant build on one gives
  the same table too. A folder left with such a file is refused before anything is built or
  searched: `build-index-only`, `complete-search` and a web search end with an ERROR that gives the
  `--split` command and says to remove the indexes built from the folder as it was (`--split` lists
  them). A search looks for CR line ends in the first two lines of each file, and
  `build-index-only`, which reads the whole genome anyway, in the whole file: one whose first 100
  lines end in LF and the rest in CRLF built at exit status 0, with 402 KB of warnings on stderr,
  and every search of that index found nothing (a search that built its own failed, on those
  warnings); its build is refused now. So is a CR that only ends a file, on its last line or after
  its last LF: the index build takes a CR on a line that ends in LF for a base, and the variant
  enricher keeps every CR, so with a PAM in the last bases of a sequence the build wrote engine
  warnings on stderr at exit status 0 (`The character () of the PAM is not part of the IUPAC
  nucleotide nomenclature`), and a search of that index failed (`searchTST (index search) failed
  (signal SIGSEGV)`) or reported other rows than the file with no CR (a variant search, from a lone
  CR after the last LF that the reference build does not read: 19 rows for 18); an upload, a
  download and `--split` write such a file with no CR. Such an index is refused as well, where its
  `.bin` names show it: it holds the old names once the FASTA is written as it should be, and
  every search of it went on losing results, and a build took it as present. An index built from a
  FASTA with a clean header line and CRLF after it has `.bin` names that do not show it, so it is
  not caught: remove it as the ERROR says. A sequence id that holds `|`, `;`, `'`, `"`, `$`, `(`,
  `)`, `&`, `<`, `>`, a
  backquote or `#` (`gi|49175990|ref|NC_000913.2|`, a PanSN `HG01255#1#CM086702.1`) is refused too:
  a search of one failed in a shell command (`/bin/sh: Syntax error`), or, for `#`, with `none
  survived post-analysis`. So is an id that starts with `-` or `.`: a search of `-x` or `.x` beside
  `b` ended with exit status 0 and the off-targets of `b` alone, none of the 7 of the other, its
  on-target among them (`grep` took `-x` for an option, and the search left out the hidden file
  `.x.fa`). Such an id is refused by a build and a search, and by an upload, a download or `--split`
  when they write its file (a FASTA of several sequences, one not named after its sequence, or one
  whose header is more than `>id` or whose lines end in CR); a file already named after such an id,
  with a plain header, goes in as it is (an upload of an archive leaves out a hidden file such as
  `.x.fa`), and its build is refused. Measured character by character through the real command line,
  `.`, `:`, `_`, `-`, `=`, `,`, `+`, `@`, `%`, `^`, `~` and `!` work inside an id, and all of them
  but `-` and `.` at its start, where `~` works only if no user of the machine has the name that
  follows it (`~root`, and `~` alone, made a search fail in its post-analysis). So do `*`, `?` and
  brackets in an id that matches, read as a shell pattern, no other name of its folder (see the next
  entry) and no file of CRISPRme's `PostProcess` folder, where the post-analysis runs: `HLA-A*01:01`
  and `KIR2DL1*001` worked, but `a*` matches `adjust_cols.py` there, and a search of it failed. A
  file whose lines end in CR alone is refused by `--split`, which says to convert it, as its first
  line holds the whole file. **If you have such a genome,** run `python PostProcess/genome_layout.py
  --split <genome folder>`, remove the indexes it names, and build again, or remove it in Settings
  and upload it again. A search of a genome of plain headers is unchanged. Covered by
  `PostProcess/test_genome_layout.py` and `PostProcess/test_build_index_only_layout.py`.
- **A genome folder in which a sequence id, read as a shell pattern, matches the name of another
  sequence is refused; a search of it failed, or lost that sequence's off-targets.** The search
  script and its post-analysis take each sequence's name into shell commands as it is, where an id
  that holds `*`, `?` or `[` is a pattern that stands for every file name it matches. In a folder
  where it matches the name of another sequence (`HLA-A*02:10` matches `HLA-A*24:02:10`, `c[12]`
  matches `c1`, `seq?1` matches `seq01`), a search failed, at its start (`basename: extra operand`)
  or in its post-analysis, whose steps for the two sequences removed, read or wrote each other's
  files; or, depending on the names and on the order those steps ran in, it ended with exit status
  0 and none of that sequence's off-targets: with `--thread 1`, all 7 rows of `HLA-X*03:01:01`
  beside `HLA-X*01:03:01:01`, its on-target among them, and so for `c[12]` beside `c1` and `seq?1`
  beside `seq01`. The HLA sequences of the GATK GRCh38 FASTA (`Homo_sapiens_assembly38.fasta`),
  whose headers this release writes as `>id` (see the entry above), are such a folder: 10 of their
  525 ids match others, and, downloaded by Direct URL, they indexed as done and every search of
  them failed at its start. `build-index-only`, `complete-search` and a web search now refuse such
  a folder before anything is built or searched (a build in the web Settings page reads `failed`),
  with an ERROR that names the files whose ids match another's (10 at most) and says to rename
  those sequences, header and file, with none of those characters, or to remove them. With those 10
  renamed (`HLA-A_02:10`), the HLA sequences build and search, and the guide `CCAGGCAGTGACAGTGCCCA`
  finds its perfect match in each of the 135 HLA-A sequences. An id that holds one of those
  characters and matches no other name is taken as before (`HLA-X*03:01:01` beside
  `HLA-Y_01:03:01:01` gives the same table), and a search of a genome with no such id is
  unchanged. **If you have such a genome,** rename or remove the sequences the ERROR names, remove
  the indexes built from the folder, and build again. Covered by `PostProcess/test_genome_layout.py`
  and `PostProcess/test_build_index_only_layout.py`.
- **A genome downloaded by Direct URL from an archive with two FASTA files of one name is refused;
  one of them no longer replaces the other.** Each FASTA file of the archive went to the top of the
  genome folder under its file name alone, so of `a/chr1.fa` and `b/chr1.fa` (the two haplotypes of
  an assembly, say) one replaced the other, with exit status 0 and no message (the web Settings page
  read `done`), and every search was of the one left. The download now stops with a message that
  names the two files, and nothing is installed (the page reads `failed`): a genome installed under
  that name is kept as it was. An upload in Settings refused such an archive already. An archive
  with no two FASTA files of one name goes in as before, file for file. **If you downloaded such an
  archive,** remove the genome it gave in Settings, with the indexes built on it, and install the
  genomes it holds each under a name of its own. Covered by `PostProcess/test_genome_layout.py`.
- **A genome download that is refused leaves the genome installed under its name as it was.** A
  Direct URL or UCSC download removed the genome installed under its name before it laid out its
  own files there (split a FASTA of several sequences, named a FASTA of one after it, wrote each
  header as `>id`), so when that step refused it, no genome was left under the name: `crisprme.py
  download --what genome --source url` of a `.fa.gz` holding `>-x`, or two sequences `>a`, over an
  installed genome ended with exit status 1, its message, and no genome under that name (in 2.8.0,
  the one of two sequences `>a` did the same). The download now lays out its files in its staging
  folder, as an upload in Settings does, and puts them in place of the installed genome only once
  that is done: one that is refused leaves it as it was (the web Settings page reads `failed`, with
  the message). A download with a FASTA whose sequence would go in under the name of another file
  of it (`all.fa` holding `>chr1`, beside `chr1.fa`), which ended in a Python traceback, now ends
  with a message: `all.fa: cannot split -- chr1.fa already exists`. A download that is not refused
  goes in as before, file for file. **If a refused download removed a genome,** install it again.
  Covered by `PostProcess/test_genome_layout.py`.

- **On a search that reads the dict-less tiers, the REF alignment of a SNP off-target marks its
  mismatches again, and its REF seed and non-seed counts include them.** A search reads them on
  every chromosome for which its variant index carries the per-sample genotype tier
  (`Dictionaries/genotypes_<dataset>/`): the indexes on Hugging Face that ship it (the default
  `NRG_3_hg38+hg38_HGDP1kGP` among them), every index built with `build-index-only --samplesID`,
  and, in this release, every index built or built again in the web Settings page. A sites-only
  index such as `hg38_mega`, an install downloaded with `--no-genotypes` and an index built
  without `--samplesID` have no genotype tier and are not affected. On an affected search, an
  off-target that SNPs create or change (a row of `alt` origin whose variants are SNPs) had its
  `Aligned_protospacer+PAM_REF` written all in upper case, so `Seed_mismatches+bulges_REF` and
  `Non_seed_mismatches+bulges_REF` left out the reference's mismatches and counted only the row's
  bulges (in `fewest_mm+b` only its RNA bulges, see the next entry), in all three criteria
  (`highest_CFD`, `fewest_mm+b`, `highest_CRISPR_BULGE`). A search of an index's dictionaries marks
  them in lower case, as CRISPRme always has. The fewest mismatches and bulges plot
  (`imgs/CRISPRme_fewest_top_1000_log_for_main_text_<guide>`, `.png` and `.pdf`, and the same plot
  on a Personal Risk Card, in the web or from `generate-personal-card`) drew the reference point of
  such an off-target at its bulge count, 0 for one without a bulge, so a variant that removes
  mismatches looked like one that adds them. Per-sample searches have had this since 2.4.0,
  population-level ones since 2.5.2. The counts of this entry and the next are of their fixes
  alone, on the code before them: the fixes of the merge step and of indel rows in this release
  change these searches too, the second with another count of indel rows. On the default index, a
  per-sample search (guide ACTGAAATCTGTAAGCAGGC, 4 mismatches, 1 DNA and 1 RNA bulge) had 71,110
  such rows of 189,949 for `highest_CFD` (68,348 for `fewest_mm+b`, 65,294 for
  `highest_CRISPR_BULGE`); for example chr7:66994206 (-) read `ACTGAAATCTGTAAGCAGGTGGG` with REF
  seed / non-seed 0 / 0 for its one mismatch. The tier path now writes these columns as the
  dictionary path does: each protospacer base that differs from the guide in lower case
  (`ACTGAAATCTGTAAGCAGGtGGG`, 1 / 0), except at a DNA-bulge column; the PAM stays in upper case.

  The pipeline breaks ties between rows by their whole line, which includes this column, so two
  more things change on these searches. First, rows that its sort holds as tied (on the CFD, then
  on the fewest mismatches and bulges) can come out in another order among themselves. On the
  per-sample search above, 5,124 of the 189,949 rows are at another place in the results and
  otherwise the same (1,285 of 732,270 in the alternative alignments); 4,598 of 205,415 on a
  population-level search of the same index (guide CTAACAGTTGCTTTTATCAC); and with guide
  CTCTCAGCTGGTACACGGCA, 1,549 of 128,730 per-sample and 2,199 of 150,103 population-level. What
  reads the rows in file order follows. The report's tables can list tied rows in another order
  (`data/cfd_ge_0.05.tsv` on all four searches; `data/cfd_ge_0.20.tsv`, `data/mmb_le_4.tsv`,
  `data/top1000.tsv` and `data/top1000_crispr_bulge.tsv` on some; `data/panel_top100.tsv` and
  `data/crispr_bulge_gt_0.tsv` on a fifth, per-sample with guide GAGTCCGAGCAGAAGAAGAA; no row of them
  changed), and so can `observed_haplotypes.tsv` (61 lines at another place on that fifth search,
  where no site reports another row) and the perfect-match banner. Any top-1000 plot can draw tied
  points in another order:
  `imgs/CRISPRme_CFD_top_1000_by_variant_effect_<guide>` on all four,
  `imgs/CRISPRme_CFD_top_1000_log_for_main_text_<guide>` and
  `imgs/CRISPRme_CRISPR_BULGE_top_1000_*_<guide>` on some, and the report's figures from the same
  rows. Where tied rows straddle the cut of 1000 of a list or plot, a by-variant-effect plot for
  instance, which of them make the cut can change too; on these four searches every list and plot
  kept its sites.

  Second, where a row of SNPs and a row with an indel (alone or with a SNP) give the same alignment,
  with the same start, score, mismatches and bulges and the same number of carriers, at least one,
  the row of SNPs is now always the one reported and the other goes to the alternative alignments
  (which a per-sample search lists). Such a tie used to fall to that whole-line sort: on the tier
  path the upper-case REF of the row of SNPs put it first nearly every time, and with the marks back
  the sample IDs would have decided. At the ties that did go to the row with the indel, the row of
  SNPs now takes its place, with its own variant, rsID, allele frequency and carriers and, where the
  two differ, PAM creation and REF and risk scores. That is 4 sites for `highest_CFD`, 4 for
  `fewest_mm+b` and 3 for `highest_CRISPR_BULGE` on the per-sample search above, and one site for
  all three criteria on each of the other three (with guide CTCTCAGCTGGTACACGGCA, chr13:110666221,
  from a row without an rsID to rs534552989), none of them in a report table. It is also what a
  search of the dictionaries reports at these sites, but for `fewest_mm+b` and
  `highest_CRISPR_BULGE` at chr14:43078564. The counts made from the reported rows' carriers move by
  one accordingly (`*.summary_by_samples.*`, `.<job>.PopulationDistribution_*`, on these searches one
  population distribution plot, and the report's population figure), and at population level
  `observed_haplotypes.tsv` can give other lines for such a site. On the dictionary path such a tie
  is rare: none on chr21 and chr22 of the default index with these three guides, where no row moved.
  The rows of a sites-only index have no carriers and keep the order they had. A search of the
  reference alone gives the same results byte for byte but for the columns of the next entry, and so
  does a search of an index without the tiers, except where a row of SNPs and a row with an indel tie
  as described: in a search of two such datasets, one of which gives the SNPs and the other the
  indels (47 of 1,837 rows on a panel with such ties planted), and, with the fix of indel rows in
  this release, in any search of one.

  The set of alignments reported (best and alternative), the scores and every other column are
  otherwise unchanged. To correct results made before, run the search again (`generate-report` alone
  reads the table as it is). Covered by `PostProcess/test_tier_ref_columns.py` and
  `PostProcess/test_representative_tiebreak.py`.

- **On a search with DNA bulges, the `fewest_mm+b` seed and non-seed counts include them, as those of
  the other two criteria do.** `Seed_mismatches+bulges_REF_(fewest_mm+b)`,
  `Non_seed_mismatches+bulges_REF_(fewest_mm+b)` and their two `_ALT_` columns left out a row's DNA
  bulge, so seed + non-seed came short of `Mismatches+bulges_(fewest_mm+b)` by it, and the same
  alignment read differently in this block and in `highest_CFD` and `highest_CRISPR_BULGE`. This
  affects every search with `--bDNA` of 1 or more, of any index or of the reference alone, and the
  haplotype tables of `assembly-search`, since at least 2.2.0. For example chr6:122734951 (+), spacer
  `A-CTGAAATCTGTAAGCAGGCNNN`, read ALT seed / non-seed 1 / 3 in `highest_CFD` and 1 / 2 in
  `fewest_mm+b`, for 4 mismatches and bulges; both now read 1 / 3. Such rows were 170,355 of 189,949
  on the per-sample search of the entry above, 181,836 of 205,415 on the population-level one, and
  114,865 of 128,730 and 133,747 of 150,103 with guide CTCTCAGCTGGTACACGGCA. Only these four columns
  change, on those rows, in `*_integrated_results`, `*_all_results_with_alternative_alignments`,
  `*_observed_haplotypes_raw`, the results database and what reads it (the web tables that list them,
  such as Query Genomic Region with the fewest criterion selected, and a Personal Risk Card's
  targets file); no plot or report table reads them. To correct results made before, run the search
  again. Covered by `PostProcess/test_seed_counts_dna_bulge.py`.

- **A search that reads the dict-less tiers no longer lists a minus-strand reference alignment with
  the letters of the plus strand.** For a variant window on the minus strand the tier path also
  wrote the window's reference alignment, from the plus-strand letters. They usually fail its
  mismatch and PAM checks, and no row came of it, but on a near-palindromic guide they pass: a row
  went out with strand `-`, those letters and the scores they get, and one more in `Other_motifs`.
  The reference search holds the true alignment of such a window when it is within the search's
  mismatches and bulges, so no row is written for a window on the minus strand now. Measured with
  the fix of the merge step in this release and before that of indel rows, on the default index,
  guide CTCTCAGCTGGTACACGGCA, 4 mismatches, 1 DNA and 1 RNA bulge, genome-wide: per sample, 5 of
  309,218 rows of alternative alignments go (chr1:147073620, chr16:12583323, chr16:85966126,
  chr6:5506134, chrX:150119503, none of them a site's best); `Other_motifs` goes down by one on
  the 5 rows of those sites in the table and the database, of 128,730, and on 76 rows of
  alternative alignments, one of which, at chr16:12583323, reads strand `+` for its CRISPR-Bulge
  block, where it paired the wrong row (`*_observed_haplotypes_raw` follows); nothing else changes
  in the table, in `observed_haplotypes.tsv` or in the tables of the report. At population level
  `Other_motifs` goes down by one on the same 5 rows, of 150,103, and one row's CRISPR-Bulge score
  moves in its last digit (chr16:1434667, `3.58341e-18` to `3.58338e-18`, as the scorer's batches
  shift); nothing else changes. With more mismatches and bulges such a row can be a site's best,
  and a site can have several. On chr19 and chr21 with 5 mismatches, 2 DNA and 2 RNA bulges, this
  guide and CTAACAGTTGCTTTTATCAC, per sample, 6 of 46,829 sites reported one as their best:
  chr19:28563886, with a CFD of 0.007, now reports a variant alignment with 0.006, and
  chr19:44589953, chr19:50205996, chr21:14062304, chr21:36964620 and chr21:45983672, with a CFD of
  0, now report a variant alignment and `Not_found_in_REF` `y`, as no reference alignment is left
  at the site (chr21:45983672 on strand `+`). `Other_motifs` goes down at 20 sites, by 1 to 10; 58
  of 316,189 rows of alternative alignments go, and at the sites that lost one the others pair
  differently across the three criteria; and the CRISPR-Bulge scores of 24 rows move in their last
  digit. At population level, with `hg38_mega` beside the default index, 5 sites change a best
  alignment so (4 for `highest_CFD`), and a sixth reads `y` for `Not_found_in_REF`. A search of
  the dictionaries, or of the reference alone, had no such row and is unchanged. Covered by
  `PostProcess/test_tier_ref_columns.py`.

- **A variant search no longer fails when a chromosome's indels hold no site of the PAM** (a
  regression of this release's CRISPRitz pin, like the SNP-only chromosome above, with another
  trigger). `build-index-only --vcf` builds one indel sub-index per chromosome from the short
  windows around that chromosome's indels. When none of them holds a site of the PAM, `buildTST`
  writes no `.bin` there, and CRISPRitz 2.9.1 refuses to search an index directory without one, as
  it should for a broken one. The search therefore stopped with exit status 1, reported at
  "Post-analysis SNPs" with `ERROR: no .bin index files found in '.../TTTV_2_fakechr3/'`, where
  2.8.0 returned the correct report. **Affected:** variant indexes built for a PAM whose sites are
  rare enough, on a dataset with few indels on some chromosome: the shipped Cas12a (TTTV, TTTN) and
  SpCas9-VRER (NGCG) PAMs on a test genome whose chr3 carries a single deletion; NRG, NGG and
  SaCas9 (NNGRRT) had a site there. The same failure followed a search that rebuilt a missing
  INDELS index itself. `pool_index_indels.py` now writes `.no_targets` into such a sub-index, after
  checking by `buildTST`'s own rule (the check build-index-only already makes) that no window holds
  a target, and the search skips that chromosome as it skips one without indels. The record travels
  with the index, so a copied or downloaded index needs nothing else. An index built before this
  fix has no record: the search then checks the fake-indel genome the index was built from
  (`Genomes/<genome>+<dataset>_INDELS/`) by the same rule, for the PAM the index was built with
  (for a PAM-less index searched for another PAM, such as `NNN_2` serving NRG, the one its
  `.pam_build` records), and skips the chromosome when no window holds a target, under
  `complete-test --chrom` too. Otherwise (a window holds one, so the index is incomplete, or that
  genome is not there) the search stops, at its "Off-targets search" stage, with an ERROR alone,
  where it used to go on through the SNP post-analysis and fail at "Post-analysis SNPs" with the
  engine's message. The ERROR says which, and how to rebuild the index, by what is there, with the
  paths of the working directory, under `complete-test --chrom` too. With the dataset's VCFs in
  `VCFs/<dataset>`: the command `rm -rf ... && crisprme.py build-index-only ...` with the PAM file
  the index was built with (for a PAM-less index, a file of its PAM in `PAMs/` or beside the
  search's PAM file) and, when the fake-indel genome is not there either, the download command for
  a downloaded index. Where that command cannot run, the steps in words: with the fake-indel genome
  there (no such VCFs, or no such PAM file), remove the INDELS index and run the search again,
  which rebuilds it from that genome (after `complete-test --chrom`, run `complete-test` without
  `--chrom` or a `complete-search`, as a run of one chromosome does not rebuild the index for later
  runs); with neither those VCFs nor that genome, the download command for a downloaded index,
  which removes only the INDELS index, so that a download that fails leaves the rest in place, and
  for one you built, putting its VCFs back; with those VCFs but no such PAM file or fake-indel
  genome, the PAM line to write into a PAM file for build-index-only. Measured through the real
  CLI: the TTTV, TTTN and NGCG searches now exit 0, both on an index built with the fix and on one
  built before it, and every result file is the same as that of the same code run on CRISPRitz
  2.9.0's search engine, which has no such refusal (but for timestamps, and the engine version
  `run_manifest.json` records), with bulges, and for TTTV without bulges; for TTTV the result
  tables are also the same as 2.8.0's. A search with no off-target at all (NGCG without bulges
  there) ends as on that engine, with `No off-targets found`, but its folder lacks the all-zero
  CRISPRitz profile tables of the skipped chromosome (`crispritz_prof/fake<chrom>_*`). NRG, NGG and
  NNGRRT searches, NNN searches and `--prune-pam` searches over 3- and 8-nt-flank PAM-less indexes
  give byte-identical results to before, and the indexes build-index-only writes are the same file
  for file, apart from the new record where a sub-index holds no target. `complete-test` on an
  intact NRG index, with and without `--chrom`, gives the same files as before, but for what
  differs between any two of its runs: under `--chrom`, the name of its per-run folder in
  `.command_line.txt`, and at `--thread 2`, the row order of CRISPRitz's raw target files. Each
  route the ERROR gives was followed as given, and led to the same result files as an intact index:
  for an NRG index, and for NRG served by the 3- and 8-nt-flank PAM-less indexes, with a sub-index
  missing or holding no `.bin`, with and without the fake-indel genome, and with and without the
  VCF folder; and from `complete-test --chrom` for an NRG index, once its results folder is
  removed, as it asks, but for those two differences. The download command and putting the VCFs
  back were followed for an NRG index only, the download with its archive served from a local
  folder, not from the hub. A run that succeeds does not reach the new stop. **No published index
  is affected:** in the copies of the published indexes on the build machine, every indel sub-index
  holds a `.bin`, except in the superseded `NRR_3_hg38+hg38_mega` archive (see its entry above),
  whose own indel logs are empty, so a search skips those chromosomes as before. **What to do:**
  nothing for an index built with 2.8.1, or for one built before whose fake-indel genome is still
  in place; otherwise follow the ERROR. Covered by `PostProcess/test_indel_empty_chromosome.py`.
- **A variant search with bulges no longer loses the indel off-targets of a chromosome whose indel
  sub-index is missing** (in 2.8.0 too). Such a chromosome was searched in its fake-indel genome
  without bulges, a fallback meant for searches without bulges: the run exited 0 with no message,
  and the bulged indel off-targets of that chromosome were missing from the report. It happened in
  two ways. (1) A search that rebuilt a missing INDELS index for a PAM-less index serving another
  PAM (`NNN_<N>` serving NRG, or `NNNNNNNN_<N>` under `--prune-pam`) built it with the search's PAM
  file, so it went into `<that PAM>_<N>_<genome>+<dataset>_INDELS` (for example
  `NRG_2_..._INDELS/NNN_2_fake<chr>`), where the search does not look; every later such search did
  the same. (2) An INDELS index copied in part. Measured through the real CLI on a test genome
  (NRG, 4 mismatches, 1 DNA and 1 RNA bulge): way 1 lost 5 of the 2446 rows with alternative
  alignments and 190 of the 191 indel+SNP co-occurrence rows of chr3 (2447 -> 2442 rows on 2.8.0),
  way 2 lost 3 rows with alternative alignments (2446 -> 2443; 2446 -> 2444 on 2.8.0), and the
  integrated results changed with both. Now (1) that rebuild uses the PAM the index was built with
  (from its `.pam_build`, else of the search's guide length), and is byte-identical to the INDELS
  index build-index-only builds, so the search's result files are those of an intact index; and (2)
  a search with bulges stops at a chromosome whose indel log, or fake-indel genome, shows indels but
  which has no sub-index, with an ERROR that says how to rebuild the INDELS index (as above), but
  for a contig with an underscore in its name, which no build gives one (next entry). A search
  without bulges still searches the fake-indel genome, and a chromosome without indels is skipped as
  before. Normal runs are unchanged: a build gives every other chromosome with indels its sub-index,
  as it does in each of the 50 INDELS indexes on the build machine that have their indel logs beside
  them. **Affected:** searches with bulges through a PAM-less index serving another PAM whose INDELS
  index was missing, and searches with an INDELS index copied in part. **What to do:** a job of an
  earlier release was affected by (1) when it allowed bulges, its `log_verbose.txt` shows `Variant
  index present (NNN_...)` (or `NNNNNNNN_...`) for another PAM, and its `log.txt` shows `Indexing
  Indels` `Start` (its `log_verbose.txt` then has no `Indels index present` line; the `Indexing
  INDELs in ...` lines tell nothing, as they may be missing from it): run it again with 2.8.1, and
  remove what it left under the other PAM's name, the `NNN_<N>_fake<chr>` (or
  `NNNNNNNN_<N>_fake<chr>`) directories in `genome_library/<that PAM>_<N>_<genome>+<dataset>_INDELS`,
  with that directory when nothing else is in it. For (2), a search now says which index is
  incomplete; once it is rebuilt, run again the searches made with it. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
- **A contig with an underscore in its name (`chrUn_KI270302v1`, `chr1_KI270706v1_random`) gets no
  INDELS sub-index, so the indel off-targets there that need a bulge are not searched; a build now
  says so.** This is so on every release, and `docs/INPUT_FORMATS.md` (section 6) does not list such
  contigs as supported: `pool_index_indels.py` takes what follows the last underscore for the
  chromosome, finds no fake-indel genome of that name, and indexes nothing. `build-index-only`, and
  a search that rebuilds a missing INDELS index (in its `log_verbose.txt`), now print one line that
  names the contigs with indels: `No INDELS sub-index for chrUn_KI270302v1: a contig with an
  underscore in its name gets none, so a search does not find the indel off-targets there that need
  a bulge`. Neither fails for want of that sub-index: `build-index-only` exits 0 (a dataset whose
  only indels are on such contigs gets an INDELS index that holds no sub-index) and a search goes on,
  as on 2.8.0. What a search then does there, measured through the real CLI on a test genome with a
  carried deletion on the contig (NRG, 4 mismatches, 1 DNA and 1 RNA bulge), the same on 2.8.0: it
  searches the contig's fake-indel genome without bulges, so a site that aligns with the guide only
  through a bulge on the carried allele is missing from the report, with exit 0 and no message. An
  indel off-target found there without a bulge stops the search at "Post-analysis INDELs", with
  `ValueError: invalid literal for int() with base 10: 'KI270302v1'`. Where that fake-indel genome is
  not there (an index copied without `Genomes/<genome>+<dataset>_INDELS/`), the search stops at
  "Post-analysis SNPs", with `tail: cannot open ...`. Off-targets on the contig that carry none of
  its variants are reported as usual, but one that carries one of its SNPs stops the search at
  "Integrating results", with the same `ValueError`. A contig with SNPs only gets no such line, and a
  search stops there only for such an off-target. So did 2.8.0, but when the contig's VCF file came
  last of the dataset's by name (as on the test genome, which has no chrX), it stopped every search
  of the dataset at "Post-analysis SNPs", with `sed: -e expression #1, char 2: expected \ after ...`.
  **Affected:** variant datasets with indels on such contigs, or with SNPs there that an off-target
  carries. No published index is: the indel logs of the 50 INDELS indexes on the build machine that
  have them beside them name no such contig, nor do the dictionaries in their working directories.
  **What to do:** to have such a contig's indels searched, and the off-targets that carry its SNPs
  reported, rename it without the underscore in its FASTA (file and header), in the VCF of each
  dataset with records on it (file and records) and in any annotation BED that names it. What was
  built under the old name is then out of date, and build-index-only keeps what it finds: remove each
  reference index of that genome (`genome_library/<PAM>_<N>_<genome>`) and what was built of each
  such dataset (`rm -rf genome_library/<PAM>_<N>_<genome>+<dataset>
  genome_library/<PAM>_<N>_<genome>+<dataset>_INDELS Genomes/<genome>+<dataset>
  Genomes/<genome>+<dataset>_INDELS Dictionaries/dictionaries_<dataset>
  Dictionaries/log_indels_<dataset> Dictionaries/registry_<dataset> Dictionaries/genotypes_<dataset>
  Dictionaries/indel_genotypes_<dataset>`), then build them again with build-index-only, with
  `--vcf` for each such dataset. A dataset without records on the contig needs nothing more.
  Measured on the test genome, the contig renamed `chrUnKI270302v1`: it got its sub-index, the
  search reported the bulged site and the off-target that carries its SNP, and the searches with that
  dataset, with another without records on the contig and without `--vcf` gave the same result files
  as on a build that named the contig so from the start (the other dataset's index keeps the old
  name, which shows only in the name of one link under `.binview/` in its search folder). With
  nothing removed, build-index-only kept everything (`... already present`) and the search stopped at
  "Post-analysis SNPs", with `tail: cannot open ...`. With the reference index kept, the searches
  without `--vcf` or with the other dataset missed the contig's off-targets, with exit 0; with the
  dictionaries kept, the report counted the contig's SNPs twice; with the annotation BED unchanged,
  the contig's off-targets got no annotation from it. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.

- **Tests that failed where nothing was wrong** (contributors; no search or report changes).
  `test_build_index_only_atomic`'s setgid test errored in any unprivileged Apptainer run, where
  `os.getgroups()` gives 65534 for every group the user namespace does not map: it now skips there,
  saying so. `test_generate_report.TestAssemblyExtraColumns` failed 3 of its tests when a combined
  report had been built before it in the same process, as under `python -m unittest discover`: it
  no longer depends on that. `test_complete_search_flags` and `test_post_process_tempdir` failed on
  a read-only code directory such as a SIF: they now pass there, and the first no longer leaves the
  data folders beside the code on a writable install.

- **A 0-bulge variant search on a downloaded index no longer loses every SNP-path variant row.**
  With `--bDNA 0 --bRNA 0` the variant search scanned `Genomes/<genome>+<dataset>/`, which a
  download does not ship. The engine printed `ERROR! ... is not a directory.` on stdout and exited
  0, so the run went on with the reference and indel rows alone, wrote a full report, and exited 0
  with an empty error log. **Affected:** every such search of a downloaded index without the
  enriched genome, on the CLI and in the web form's Advanced panel, in 2.8.0 and back to at least
  2.3.3 (searches with bulges, and installs with the enriched genome, were right). **What to do:**
  run those searches again. Measured on hg38 with the default index (TRAC guide, 4 mismatches, NRG,
  per-sample): 483 rows before, 722 after, which for this guide is byte-identical to a scan of the
  enriched genome. The 240 new sites are all variant sites, and one site listed before is now
  reported 1 nt over, at a variant alignment that scores higher; sites found only in the variant
  genome go from 80 to 319; sites at CFD >= 0.5 from 7 to 11, >= 0.2 from 33 to 50, >= 0.05 from
  156 to 233; and 207 of the sites already listed change too (a variant alignment that scores
  higher, or another motif counted). The search now takes the installed variant index when any
  chromosome's enriched FASTA is missing, with the same `-var`, `--max-edits` and `--prune-pam` a
  search with bulges passes; with the whole enriched genome in place it scans it, as before. The
  index has the limit every index has at contig ends and `N` runs (below): for a site whose guide
  lies within N − 1 nt of one (2 nt for the published indexes), the scan reports the variant
  version and its carriers, and a search through the index only the reference version, from the
  scan of the reference genome, when that is within the mismatches asked.
- **A search the engine never ran now fails the run.** `crispritz.py` says what stops it before it
  searches on stdout and exits 0, so a search that never ran passed for one that found nothing.
  After each reference and variant search, the run now stops with an ERROR when the targets file is
  missing. A search that finds nothing writes a file with its header alone, and still ends in an
  empty result with exit 0.
- **A search with bulges no longer ends empty when its library cannot be written to and lacks the
  index it needs.** On a read-only mount, without write permission or with no inode left, the failed
  lock read as "another run built it", so the search went on with a folder that is not there: the
  engine said so on stdout and exited 0, and the run ended `No off-targets found ... completed
  successfully` (web: "Finished: no off-targets were found"), or, with the reference index in place
  and the variant index short, reported the reference and indel rows alone. **Affected:** every
  such search, in 2.8.0 too. It now stops with one ERROR that names the index and the library and
  says what to do, shown on the web job's banner. A library that holds the right index still serves
  the search, byte-identically, and a 0-bulge search, which opens no index there, still runs.
- **A guide that does not fit an index's stored window is no longer searched on it.** An index name
  carries no guide length, and the engine searches any guide: past what the index stores, it
  reports wrong alignments or loses sites at exit 0, or, on a site planted at a contig start, stops
  with `Non-numeric start coordinate`. A 20-nt search at default flags builds `<PAM>_2`, and the next
  23-nt search took it: on real chr22, 181 rows against 1003, the on-target no longer a perfect
  match. On the web, the same, with no message. `build-index-only` with a 23-nt PAM file over a 20-nt
  index said "already present" and rewrote `.pam_build` with the wrong geometry. **Affected:** a
  search of a guide longer than its index was built for, in 2.8.0 too; every published index is
  20-nt, so a 20-nt search is not (nothing changes for it), a 21-nt one at up to 1 bulge of each
  type on an N=3 index is not, and any longer one is. **What to do:** run those searches again;
  they build their own index, or say how to clear the way. The length comes from the index's own
  `.bin` header, so no index needs rebuilding. An index a newer CRISPRitz built, in the format
  CRISPRitz 2.9.2 and later write (index format 2, whose files open with a magic and the format
  version), holds no length this release reads: a search that takes one, found or built with such a
  CRISPRitz installed, stops with one ERROR that says what it is and what to do (use a CRISPRme that
  reads index format 2, or install CRISPRitz 2.9.1 and build the index again), and so does
  `build-index-only` with one in place that it would keep (the reference, the variant or the INDELS
  index), or one the installed CRISPRitz writes. Measured through the real CLI with that CRISPRitz
  built into the image, they had read its format version as the guide length ("built for a 2-nt
  guide") and its magic as a leaf count ("not a tree ..., as when the disk fills up"); with
  CRISPRitz 2.9.1 the same runs give the same tables. No 2.8.1 install pulls such a CRISPRitz:
  the image and `install_from_source.sh` pin 2.9.1, and the Bioconda recipe `>=2.9.1,<2.9.2`. A search takes an index only where built + N >= guide +
  bDNA + (1 if bRNA > 0), measured through the real CLI on 246 cells: every cell inside is
  byte-identical to a dedicated index, every cell outside differs. When an index of the same name
  is in the way it stops with one ERROR naming the index, both lengths and what to do; the index is
  never removed. An index built for this very guide is now preferred to one that merely fits, so a
  20-nt search beside a 23-nt `NGG_2` and a 20-nt `NGG_3` takes `NGG_3`. `build-index-only` stops
  the same way, and leaves the record of an index built for another guide. The web form's bulge
  caps and PAM list follow the same rule, and two PAM files that read the same show their guide
  length ("Cas12a · TTTV (23 nt)").
- **A guide longer than its PAM file's spacer is refused.** The search kept the spacer's length of
  it and dropped the rest with no message: a 23-nt guide with `20bp-NGG` was searched without its
  last 3 nt (2 rows against 35, and none at 0 bulges, at exit 0), under the full guide's name. The
  pre-flight now says which bases would go and what to do. Padding with extra `N`s drops nothing
  and stays accepted. CLI `--guide` only: the web form and `--sequence` never pass such a guide.
- **A search that enriches a dataset itself no longer loses targets that need an indel and a SNP
  beside it, and no longer keeps a blind INDELS index.** Only `build-index-only` put the SNP codes
  on the fake-indel genome. A search that found the dataset's VCFs and no index for its PAM or
  bulges (a download plus `--what vcf`, or a renamed genome folder), and the web's first 0-bulge
  search of a dataset that has only its VCFs, enriched the genome without them, so a cis SNP+indel
  site read 2 mismatches, CFD 0.096, listing a carrier that does not hold both, against 1 mismatch,
  CFD 0.353 on a built index; the INDELS index stayed in the library, so every later run gave the
  same, and `build-index-only` said "already present". Exit 0, no message. **Affected:** every such
  search since 2.5.0. Now the search script overlays the codes before the fake-indel genome goes
  into place, and again before it builds an INDELS index from one that is there, and writes the
  phased indel genotype store where it enriches. Measured on the walk genome: the on-demand search
  is byte-identical to a search of a built index, in all 41 result files. The first web search gets
  the same row, carrier and scores (it needs the indel track's fix, in this release, to resolve the
  codes without tiers); it has no registry, so its report says no co-occurrence table is written
  rather than "unaffected", and the CLI warning says the same (a report made again from a copy of
  the results away from the install, where that cannot be told, keeps "unaffected" for a
  population-level run, as before). **What to do:** an INDELS index built by such a search stays
  blind until it is rebuilt. Run `build-index-only` with the dataset's VCFs and
  your PAM: it now codes a fake-indel genome it finds without codes, removes the INDELS indexes of
  every PAM and budget of that dataset built before the codes, says which, and builds its own (a
  later search rebuilds the others). A genotyped dataset with no `Dictionaries/indel_genotypes_<vcf>/`
  (an index built before 2.5.0) now gets a `WARNING` that says so and what to run.
- **A downloaded variant index asked for a PAM or bulge budget it does not hold now says so.**
  The search tried to enrich a dataset with no VCFs, and failed in `mv: cannot stat ...` four times
  and `Genome enrichment failed`, after a pre-flight that printed `[OK] ... precomputed variant index
  ... present`, leaving empty folders; the web form promised that a first search would prepare the
  genome, and its job ended "please remove it and try to submit again". This is likely the first
  slip on the default index (searching `NRG_3_hg38+hg38_HGDP1kGP` with NGG), and on a K=8 index
  without `--prune-pam`. **Affected:** a search with no index for its PAM or bulges and neither the
  dataset's VCFs nor its enriched genome installed, in 2.8.0 too. It now stops before it does
  anything, on the CLI before the pre-flight's `[OK]` and on the web job's banner, with one ERROR that
  names the genome, dataset, PAM, bulges and guide length asked for, the variant indexes installed,
  and what to do: search with the bulges an index of this PAM holds, `--prune-pam --index-path
  genome_library` where a wider all-N index serves it, or `download` or `build-index-only` another
  index. Web serving of a wider all-N index is not new here: the web still needs an all-N index of
  exactly the PAM's width. An install with the enriched genome (the build-inputs kit) or the VCFs
  still builds the index and is byte-identical to a dedicated one. The README no longer says a
  PAM-less index must be searched from the CLI: only a wider one must.
- **With no `--bDNA`/`--bRNA`, a search no longer becomes a 0-bulge search with no message.** The
  default of 1 bulge of each type is capped by the indexes a search can reach, but that cap looked
  only at the working directory's `genome_library/` and `Genomes/<name>`. A library at
  `--index-path` (variant, K=8 with `--prune-pam`), or a `--genome` folder elsewhere (the
  shared-mount recipe), gave a cap of 0, and the run was byte-identical to a 0-bulge search: 33 rows
  against 38 with bulges, with a site that needs a bulge reported at a worse alignment or not at all
  (2.8.0 too). **What to do:** run those searches again, or pass the bulges. The cap now reads
  `--index-path` and the `--genome` folder given, so these searches take their 1 bulge, and where no
  usable index holds bulges (K=8 without `--prune-pam`, an index of another PAM) a `WARNING` says
  the search runs with none, and what to pass.
- **`complete-test` works on any index you installed, and no longer breaks when the K=8 index sits
  beside the default.** Its smoke wrote the example guide and a 3' PAM from the index's name, so on
  a PAM wider than 3 nt it stopped with `guide ... is 23 nt but the PAM expects 28 ... Fix: append 8
  'N'`, about a file it rewrites each run; it took the first name in sorted order, so
  `NNNNNNNN_3_hg38` shadowed `NRG_3_hg38`, on the install `DOCKER_QUICKSTART` describes; and on a 5'
  index it searched the wrong side, which under `--chrom` read `PASSED ... holds none of its
  off-targets` with sites there. Each failed run also left a PAM file such as
  `PAMs/20bp-TTTV-SpCas9.txt`, which the web then offered as a nuclease. **Affected:** the
  documented self-test on a K=8, SaCas9, Cas12a or other non-SpCas9-class index, in 2.8.0 too; a
  leftover `PAMs/20bp-<motif>-SpCas9.txt` of such a motif can be removed. It now takes a 3' index of
  a 3-nt PAM first, and an all-N index wider than 3 nt last; searches a wider 3' all-N index for NRG
  with `--prune-pam`, byte-identical to a dedicated NRG index (measured on the synthetic genome,
  and on the published K=8 indexes on chr22); and searches any other for its own PAM, on the side
  its `.pam_build` records (3' without one; for a dedicated index, what the PAM files of its motif
  in `PAMs/` say), with the example guide padded for it. It writes nothing into `PAMs/`: it uses the
  file there that holds the line it needs, or a hidden one for the run. With another nuclease's PAM
  an empty result on the whole genome passes too, as it does on one chromosome. An index whose side
  cannot be told stops in one line. Results on NGG, NRG, NNN and a 2-nt PAM are unchanged.
- **A variant search on a genotyped index without `--samplesID` no longer ends in a `KeyError` after
  the whole search.** The help marks the flag optional and the run warned "using an empty sample
  set", a case meant for sites-only panels; on a genotyped index the rows name carriers, and
  `process_summaries.py` looked each up in an empty table at "Creating images", then deleted the
  results: 65 to 78 minutes on hg38 for a `KeyError`. An empty or a partial list on a folder with no
  VCFs did the same, and so did `complete-test` with a sample list under another name (2.8.0 too).
  **What to do:** run those searches again. Without the flag the search now uses
  `samplesIDs/<dataset>.samplesID.txt` for each dataset that has one, which an index download
  installs, and says so: every result file but its own list `.list_samplesID.txt` is byte-identical
  to passing that list, in both analysis modes. With no such file the search still ends with the
  same `integrated_results`, with empty by-sample and population tables, and the warning says so.
  A carrier a partial list does not name has no row in those tables.
- **`build-index-only` records the folder an index of a genome outside `Genomes/` holds, and a
  Direct URL download no longer replaces a genome something is built on.** A search takes a
  genome's index by the folder's name alone. `complete-search` records, in an index it builds of a
  folder outside `Genomes/`, which folder it holds (`.genome_source`); `build-index-only` did not.
  So after `build-index-only --genome projA/ref`, a search of another folder named `ref` reported
  the results of `projA/ref`, byte for byte, at exit 0; `build-index-only` of a folder named like
  one in `Genomes/` built that one's index from it, which every search of the `Genomes/` folder
  then took; a recorded folder whose FASTA files had changed was "already present"; and the search
  of a personal-assembly haplotype by its own path (`Assemblies/<individual>/<haplotype>/genome`)
  was refused on the index `build-index-only` had just built for it, which it had said would be
  reused. `download --what genome --source url`, the Settings page's Direct URL, put another
  genome's FASTA files under an installed name and kept its indexes, which every search then used,
  with the old genome's sequences and positions. **Affected:** `build-index-only` of a folder
  outside `Genomes/` and a Direct URL download over a genome with an index, in 2.8.0 too; the
  refused haplotype search is new in 2.8.1. **What to do:** a search run on such an index, or on a
  genome downloaded over another, searched the other folder's sequence: remove that index and run
  it again. `build-index-only` now records the folder, refuses one named like a folder in
  `Genomes/` or like another folder an index records, and builds an index again when its folder's
  FASTA files have changed (`Removing index ...`), as `complete-search` does, with the same
  messages. Of a folder in `Genomes/` it records nothing, as before, so a downloaded index, and a
  genome renamed together with its indexes, are served as before. A Direct URL download under the
  name of a genome with an index or a variant genome installed is refused before anything is
  fetched, by the rule a Settings upload has (remove them, or use another name), and the Settings
  page says so at once. `--source ucsc` and the HuggingFace download give the same assembly by
  name, and are not refused.
- **A variant search reads the variant index at `--index-path`.** `--index-path` served the
  reference index alone: the variant index and its INDELS index were looked for in the working
  directory's `genome_library/` only, so a variant search of a complete library kept elsewhere
  (a shared, read-only install) stopped with `mv: cannot stat` and `Genome enrichment failed`
  when the dataset's VCFs were absent, and with them it enriched the genome and built a second
  variant index in the working directory, which every later run there took, also when the one at
  `--index-path` was newer. **Affected:** `complete-search --vcf --index-path <library>` of a
  library other than the working directory's, in 2.8.0 too; the web passes no `--index-path`.
  **What to do:** a variant index such a run built in the working directory can be removed. The
  variant index, with its INDELS index beside it, is now read from the `--index-path` library first
  and then from the working directory's, and the run says which library it read
  (`Variant index library: ...` in `log_verbose.txt`). Measured on the synthetic genome: a library
  elsewhere with no VCFs here, the same with the VCFs here, a K=8 library with `--prune-pam`, and
  a stale index of the same name in the working directory are each byte-identical to the search
  with that library in place, and build nothing. A search builds no index in an `--index-path`
  library: a variant index there without its INDELS index stops the run with a message.
  Dictionaries, sample lists and the genome are still read from the working directory, and a
  search with the library in the working directory, with or without `--index-path`, is unchanged.
- **A search refuses an index that lacks a chromosome's trees, instead of losing that
  chromosome's sites.** A search took an index holding one `.bin` for a whole one: a partial copy,
  a build that an earlier release left stopped between chromosomes, or an index of a genome folder
  that gained a FASTA since, gave no row on the missing chromosome, at exit 0 and with no message,
  on the CLI and the web, and `build-index-only` said "already present". Measured on the synthetic
  genome: a reference index without chr2's tree gave 6 rows instead of 7, a variant index without
  chr1's another result, and a FASTA added after the build was not searched. **Affected:** searches
  of such an index, in 2.8.0 too; a stock install is not (every index on the lab's hg38 install, 32
  reference and variant indexes of NRG, NGG, NRR, NNN and K=8, holds every tree). **What to do:**
  nothing, unless a search now stops on an index: remove it and build it again, or download it
  again, as the message says. Before it opens an index (with bulges, or for a variant search
  without bulges and without the whole enriched genome), the search now checks that it holds the
  trees of every FASTA of the genome folder, and a variant index those of every chromosome of its
  dataset, and stops naming the FASTAs that lack them; `build-index-only` does the same before it
  says "already present". A FASTA with no site of the index's PAM (all N, shorter than a window, no
  PAM) has no tree in a healthy index and is not counted, by the engine's own rule for a site, and
  neither is a variant chromosome the genome folder has no FASTA of. A gap in a chromosome's
  numbered trees, and a full last tree of a chromosome that holds more sites, count as missing
  trees too; a chromosome's last trees lost after a tree that is not full cannot be told from the
  files. The check reads the first line of each FASTA, and a FASTA in full only where it has no
  tree or a full last tree: under a second on hg38.
- **An install without the genotype store runs the population-level analysis, rather than a
  per-sample one that loses off-targets.** On a phased dataset the CLI chose the per-sample
  analysis from the registry's data type alone, so after `download --what index --no-genotypes`,
  or a download whose genotype store did not come, it ran per-sample with no genotypes to read:
  it joined no two SNPs of one site, dropped the off-targets that only two SNPs of a sample
  together reach, and still printed "CONFIRMED cis phasing + named carrier samples". `--per-sample`
  did the same. Measured on the synthetic genome with a site whose PAM only two cis SNPs make (1
  mismatch, CFD 0.6, carrier S4 on a full install): the no-store default reported it at a 2-mismatch
  bulged alignment, CFD 0.095. **Affected:** CLI searches of such an install, in 2.8.0 too; the web
  already offers population-level alone there. **What to do:** run them again, or install the store
  (`download --what index` without `--no-genotypes`) for carriers. A dataset with a registry and
  neither a genotype store with anything in it nor the classic per-sample dictionaries now runs
  population-level, by default and under `--per-sample`, with one `WARNING` on stdout that names
  the dataset and the download. On the synthetic genome the no-store default is then byte-identical
  to its `--population-level` run and reports that site at 1 mismatch, CFD 0.6, with no carriers; it
  lists the same sites, at the same CFD and with the same variants, as a full install's
  `--population-level` run, and a row of SNPs names no carrier. Where an indel and a SNP make a site
  together, its row names both and lists every carrier of the indel, S4 and S5, as no genotype says
  which of them has the SNP too (PUTATIVE in the SNP+indel companion), where a full install lists S4
  alone, whose haplotype holds both (CONFIRMED). A full install, and one with the classic
  dictionaries, are unchanged. The download help, the README, docs/PRECOMPUTED_INDEXES.md and
  docs/TROUBLESHOOTING.md say what such an install runs. A sites-only dataset, such as `hg38_mega`,
  has no genotypes at all, and so no store to install: it is not one of these, and leaves the mode
  of the datasets beside it as it is (a phased one stays per-sample); its own rows are the
  population-level analysis's in any mode (the entry on sites-only datasets).
- **An index is taken only for a search on the side of the guide it was built for, a dedicated
  one too, and an all-N index that records no side is taken for a 3' PAM only at the published
  widths.** 2.8.1 takes an all-N index of a 3' PAM's width, and one with no `.pam_build` for a 3'
  build. A search before 2.8.1 records no side, so after a 5' all-N search (`23bp-NNNN-NO-PAM`,
  which builds `NNNN_2_<genome>`) a 3' search of that width (`NGCG`) read the 5' index: on the
  synthetic genome it found no site against 2 on a dedicated index, which 2.8.0 built. This is new
  in 2.8.1. The side check also covered all-N indexes alone, while a dedicated index is named
  after its PAM by PAM files of both sides (`N...NTTTV 4` and the 5' `TTTV-Cas12a` both build
  `TTTV_<N>_<genome>`): a 5' TTTV search of a 3' build reported nothing, or a site that is in no
  chromosome, at exit 0, as in 2.8.0, and `build-index-only` with the other side's file recorded
  that side over it. **Affected:** a search beside an all-N index of a width other than 3 and 8
  that a search before 2.8.1 built, and a search beside a dedicated index built from the other
  side's PAM file; no published index and no shipped PAM file. **What to do:** run those searches
  again. Now every index's recorded side must be the search's side, and a search records the side
  of every index it builds, dedicated ones too (build-index-only always did). An all-N index with no
  record serves a 3' PAM at widths 3 and 8 (the published `NNN_3` and the K=8 archive, which has
  no `.pam_build`) or under the search PAM's own name, as before; at any other width the search
  builds its dedicated index, as 2.8.0 did, and under `--index-path` it stops and says why. A
  dedicated index with no record is used as before. Where the index of the build name is in the
  way the search stops with a message, and `build-index-only` refuses another side's PAM file over
  a recorded index, and a 3' file over an unrecorded all-N index of another width (a 5' file is
  recorded, as the messages say to do). Measured on the synthetic genome: the NGCG search after the
  5' NNNN search, and the 5' TTTV search beside a 3' `TTTV_2` and a 5' `TTTV_3`, now equal a
  dedicated index; NRG on an unrecorded `NNN_3`, and on an unrecorded K=8 with `--prune-pam`, is
  unchanged. A dedicated index built for the other side's PAM by a search before 2.8.1 records
  nothing and cannot be told apart: remove it if you built one from such a PAM file.
- **A search under `--prune-pam` warns of an all-N index built by CRISPRme 2.8.0 or earlier.**
  CRISPRitz before 2.9.1 dropped every target whose stored PAM/flank region held an `N`. Searched
  for a PAM narrower than its flank, which 2.8.1 newly allows, such an index misses the sites
  whose flank holds an `N` past the bases the PAM covers, at exit 0 (the review measured 3 of 36
  rows on its synthetic genome, and 2.6% of the indel rows on chr22, none above CFD 0.2; on 2.8.0
  such an index could serve only an all-`N` search of its own width, where nothing is lost). The
  index records nothing of the engine that built it, but its trees do: an index built since stores
  an `N` in the flank of a target beside one, and no earlier build does. Measured on the lab's
  hg38 indexes: the `NNN_3` variant and INDELS indexes of the 2.8.0 install hold none, and the
  rebuilt `NNN_3`, the published K=8 and the final K=8 builds do. **Affected:** a K=8 or `NNN_3`
  index you built or downloaded before 2.8.1 (the published replacements are not), searched with
  `--prune-pam` for a narrower PAM. **What to do:** remove such an index and build it again, or
  download it again; `build-index-only` keeps an index it finds in place. The search now reads the
  trees of the index it takes when its flank is wider than the PAM (the smallest tree of a contig
  whose sequence holds an `N`, and of the INDELS index the smallest one of several indel windows)
  and prints a `WARNING` on stdout, also written to the run's `Params.txt`, when it finds no
  `N`; it does not stop the search. An index of the PAM's width is not read, and nothing else
  changes.
- **A working folder the run cannot write to stops it in one line, not a traceback.** Run from a
  read-only folder (a `:ro` mount, no write permission), `complete-search`, `assembly-search` and
  `build-index-only` ended in a raw `OSError` traceback from creating `Genomes/`, the output
  folder or `genome_library/`, and `crisprme.py` with no arguments printed one instead of its help.
  **Affected:** those runs only, in 2.8.0 too; nothing was written or lost. They now stop with
  `Error: cannot create <folder>: <reason>. Run CRISPRme from a working folder you can write to.`,
  and the help prints. A read-only folder that already holds every CRISPRme folder, searched with
  an `--output` elsewhere, runs as before.

- **An indel off-target near a SNP was reported with SNP alleles its carriers do not have, and with
  PAMs the search never matched.** Since 2.5.0 the build writes the IUPAC code of every SNP onto the
  fake-indel genome the indel search runs on, and the search engine matches a code when any of its
  alleles fits, so one raw indel target stands for every combination of the SNP alleles around the
  indel. The indel post-analysis then turned each code into one base, column by column, without
  looking at any genotype: the guide's own base when it was one of the code's alleles, otherwise the
  code's first allele. A row therefore paired one individual's indel with SNP alleles that only
  other individuals carry, scored as if real, under the indel's carriers: on the default HGDP + 1000
  Genomes index, guide `ACTGAAATCTGTAAGCAGGC` (NRG, 4 mismatches + 1 DNA + 1 RNA bulge) reported
  chr18:12386282 at CFD 0.708, fifth in the report's top-1000 table and 18th in the recommended
  validation panel, as "1 cis carrier": NA19072, who is homozygous reference at the SNP the
  displayed sequence uses (its only carrier, NA20541, does not carry the deletion); on NA19072's
  real haplotype the site is 5 mismatches + 1 bulge, past the search's limit. Its mismatch counts
  were understated the same way. And a code at a PAM position always read as its first allele: an R
  (A/G) at NRG's last base was shown as A, a PAM the search never matched, with a CFD of about 0
  (1,088 rows of that search's highest-CFD block). On 9,255 of its 13,622 highest-CFD indel rows no
  listed carrier's real haplotype, built from the phased VCF, holds the reported sequence, mostly
  for this reason (the next two entries and METHODS §8 cover the others). Each carrier haplotype of
  the indel now gets its own allele at every code: the SNP alt it carries, else the reference base,
  several alts on one VCF line each read against its own allele index. Haplotypes left with the same
  sequence form one row that lists only their individuals and names the SNP alts they hold next to
  the indel (`Variant_info`, rsID, MAF), with its mismatches recounted, its CFD, CRISPR-Bulge and
  risk scores computed on that sequence, and its PAM creation flagged from it; a sequence whose PAM
  does not satisfy the searched motif, or that is past `--mm` or `--max-total-edits` (now passed to
  the indel post-analysis), is dropped, so a raw target no individual has gives no row. A row's
  carriers stay on its sequence through the merge of a cluster's rows, which would otherwise add
  those of another alignment of the site naming the same variants; so do those of a row of an
  equal-length multi-base record (an MNV, such as AC>TG), which the search treats as an indel. An
  individual whose genotypes cannot put an allele on the indel's haplotype (an unphased call, a
  missing allele, no phased indel genotype, a haploid call against a diploid one; every carrier of
  an MNV, which the phased indel genotype store leaves out) gets the strongest target its own
  alleles allow, and is never listed for an allele it does not carry; with no genotypes at all (a
  sites-only panel, `hg38_mega`), the row holds the strongest target the codes allow (its PAM
  always one the motif allows) under the carriers the
  indel log lists, as before. Both are PUTATIVE in the SNP+indel co-occurrence companion
  (`*.indel_snp_cooc.tsv`), which now lists the SNP alts the indel post-analysis's rows require, a
  SNP that creates the PAM included, with those rows' carriers (CONFIRMED when each was resolved on
  its phased haplotype), and no longer a SNP set no carrier's haplotype holds. The companion is
  written before the merge picks each cluster's alignments, one line per engine alignment, required
  SNP set and phase (an alignment that needs no SNP writes no line), so a population-level run,
  which writes no alternative alignments, can have lines with no row in its tables. A classic
  dictionary install (an
  index built without `--samplesID`, as the web Settings page built them before this release, so
  with no Tier-1 genotype store) reads the genotypes from the per-sample SNP dictionary: there,
  the codes used to stay in the rows unresolved, scoring CFD 0, so the real SNP+indel targets were
  lost. With `CRISPRME_INDEL_SNP=0` at search time each code reads as the reference base, as on a
  fake-indel genome built without the overlay. A target with no code (every target of a fake-indel
  genome built without the overlay) that holds its indel, as the entry "An indel row could hold no
  base of its indel" defines it, is not changed by this fix. **Affected:** every search since
  2.5.0 on an index that carries the overlay, which every published variant index and every index
  `build-index-only` builds do: from the web or the command line, per-sample and population-level
  (the indel post-analysis does not depend on the mode), with any PAM. A search that enriched its
  dataset on demand (`CRISPRME_ALLOW_ONDEMAND_BUILD=1`, as the benchmark mode of `complete-test`
  does) built its fake-indel genome without the overlay, so no code reached its post-analysis
  (since 2.8.1 it has the overlay: the entry on a search that enriches a dataset itself). Under a
  PAM whose motif allows the base the misreading produced (NRR allows A at its last position) the
  PAM stayed legal while the CFD collapsed, so nothing showed. Reference and SNP-only off-targets
  are not affected, though a cluster whose indel candidate goes away shows its next best in the
  results and its `Other_motifs` count changes. Measured on the search above, rerun genome-wide
  with the fixes of this entry and the next two: its highest-CFD indel rows went from 13,622 to
  5,047, of which 5,037 hold their sequence on a listed carrier's real haplotype, read literally
  from the phased VCF, and every one of those spans its indel; the other 10 hold a SNP at a split
  multi-allelic site, which the enriched genome and its dictionaries read as reference (METHODS
  §8), as does one of the 103 carriers of one more row; no ALT row shows a PAM NRG does not allow
  (3,581 over the table's three blocks before); chr18:12386282 is now the SNP row of its real
  carrier, NA20541, at CFD 0.14, out of the top-1000 table and the panel; that table lost 66 indel
  rows and gained 2 (chr5:170503088 and chr10:103001704), from 88 to 24, the recommended panel
  went from 259 to 234 sites, the CFD >= 0.5 tier from 67 to 57, the CRISPR-Bulge sensitive tier
  from 162 to 156 and the <= 2 mismatches + bulges tier from 9 to 7; the SNP+indel companion went
  from 590 to 336 lines, all CONFIRMED and the same in both analysis modes; all 755,576 reference
  and SNP-only candidates are unchanged, and so are the two perfect matches and the
  multiple-perfect-match banner. HPRC is the only published index with MNVs (chr22: 9,079 of its
  195,813 indel-search records); the search above on HPRC's chr22 lists 8 MNV row blocks per
  block, and the MNV part of the merge fix changes none of their carriers. The search takes as
  long as before and peaks at the same memory, its search stage's; an indel post-analysis worker
  holds about 2% more. **What to do:** run such searches again; results of earlier releases are
  corrected only by running the search again. Covered by
  `PostProcess/test_indel_haplotype_resolution.py`, `PostProcess/test_indel_overlaid_rows.py`,
  `PostProcess/test_representative_tiebreak.py`, the synthetic SNP+indel oracle through the real
  CLI (`test/indel_snp_oracle`, run by `validate-benchmarks.yml`) and a release gate,
  `scripts/check_alt_pam.py` (no ALT row may show a PAM the searched motif does not allow, and no
  REF or ALT target of any block may hold an N, an assembly gap of the genome), run on the chr22
  benchmark outputs.
- **An individual whose haplotype carries a second indel that changes an indel off-target was
  listed for the sequence of the first indel.** The fake-indel genome holds one indel per
  window, so such a haplotype does not have the window's sequence; on the synthetic oracle a
  sample carrying two deletions on one haplotype was listed on both windows' rows, for
  sequences neither of its haplotypes has. Whether a carrier whose haplotype holds a second
  indel inside the target has the row's sequence is read on its bases: the reference around the
  target is written with every indel of the phased indel store the haplotype carries there and
  its own SNP alleles, and the carrier is listed when the target's bases are on it at the
  target's place, from the target's first base on or up to its last, through the indel (not in
  a stretch that another of its variants lengthened, pushing the indel out of it). In a tandem
  repeat a unit gained or lost at the target's edge leaves the bases as they are (the repeat
  puts the same bases back), so such a carrier stays listed; one whose second indel changes
  them does not, nor one whose bases cannot be read either way: each way meets an unphased or
  missing SNP call, an unphased indel that does not itself change the target, or two of its
  variants that overlap, two records of one base: two other indels
  whose REFs share a base (any two at one position do, even a duplication and the deletion of
  the same bases, which together would leave the reference), or a SNP on a base of the second
  indel's REF. A reading from the first base on meets no variant upstream of it, one up to the
  last base none downstream. A carrier whose second indel inside the target has an unphased
  call cannot be tied to a haplotype and stays listed, PUTATIVE, as every carrier on an
  unphased panel is, even when another, phased second indel on its indel haplotype changes
  the target. Nor is one listed whose second indel, on the same phased slot,
  overlaps the first at another position (one complex allele that the VCF splits into
  overlapping records). A variant only on bases the indel deletes changes nothing. Another
  allele at the indel's own position, an indel or a SNP, is an alternative to it (a haplotype
  holding both would be one record), so a slot phased with both is a phasing artefact: on chr22
  of the HGDP + 1000 Genomes panel, two indels at one position share a slot for 40,679 carrier
  pairs against 2,463,324 on opposite slots, a SNP and an indel for 2,404 against 167,951. And
  an unphased carrier of two overlapping records is a multi-allelic heterozygote, one allele
  per haplotype. All of these stay listed. This needs the phased indel genotype store, which
  the build writes with the SNP overlay; a row left with no carrier is dropped. Its own target,
  which needs both indels, is still not searched (METHODS §8). On the search above, 693 targets
  with no SNP code lost such carriers. On chr22 of that panel (per sample, five guides),
  checked against haplotypes built from the phased VCF: of the 7,046 carrier listings the rule
  takes off the integrated results' highest-CFD block (96 rows are left with none, and
  dropped), 7,043 are of a haplotype that lacks the row's sequence within 40 bases of it, 2 of
  one that has it a few bases away without the indel (a repeat's copy, a target of its own) and
  1 of a slot with two overlapping records; no carrier left on an indel row lacks the row's
  sequence because of a second indel.
- **An indel row could hold no base of its indel.** The engine reports every target on a fake-indel
  window, also one wholly in a flank, whose sequence is plain reference (or SNP) sequence that needs
  no indel. The post-analysis computed a row's Start with no shift for the indel and kept the row when
  that range met the indel's, so a target in the downstream flank of a deletion passed: a row with
  the indel's name and carriers, at a coordinate the indel shifts, that could take the real row's
  place in its cluster. On the search above, the top-1000 table showed at rank 628 chr1:69462019 as a
  2-carrier off-target of a deletion and a SNP, hiding the SNP row chr1:69462021, chr1_69462039_A_G,
  with its 3,459 carriers, so a common off-target read as a rare one (it is now at rank 588), and
  at rank 983 chr8:28788278, the sequence of the reference row chr8:28788272 reported again 6 bp off.
  A raw target now has to hold a base the ALT adds or replaces, or the two bases on either side of the
  deleted ones; one wholly inside a long insertion does, and stays. **Affected:** every search with a
  variant dataset that has indels, in every release, on any index, with or without the overlay and
  with `CRISPRME_INDEL_SNP=0` too. No off-target goes with these rows: on the search above, each of the
  1,566 distinct rows this removes has a reference or SNP row at its site that holds its sequence
  (1,560) or names the same SNP alts (6). **What to do:** run such searches again. Covered by the
  same tests.
- **An indel off-target with an RNA bulge was reported one base early, and the reference sequence
  shown for one with both bulge types was a base off.** For an indel row whose bulges are all RNA
  bulges, the indel post-analysis took the search engine's cluster position for the row's start on
  the + strand (on the - strand with a 5' PAM such as Cas12a's), so its `Start_coordinate` was one
  base per RNA bulge before its first aligned base, where the reference and SNP rows of the same
  alignment start, and its variant's position in `Variant_info_spacer+PAM` was off by as much. The
  reference sequence of an indel row (`Aligned_protospacer+PAM_REF`) was read with one base too many
  per RNA bulge and was given its gaps only on a row whose bulges are all RNA bulges: a row with an
  RNA and a DNA bulge (`RNA,DNA`) showed it without the gap and one base too long, out of register
  with the guide from the RNA bulge on, and with a 5' PAM an RNA-bulge row's lost its first PAM
  base. The columns computed from it followed: the REF seed and non-seed mismatches + bulges,
  `CFD_score_REF`, the CFD risk score, `PAM_creation`, and `CFD_score`, which is the higher of the
  REF and ALT scores. Each indel row now starts at its first aligned base and shows as REF the
  reference bases its target covers from there, gapped as its target. Its cluster position, which
  the merge of a cluster's rows sorts and groups on, moves with its start: an RNA-bulge indel row on
  the + strand no longer comes first among rows that tie exactly, and, as for SNP rows, the merge
  lists once the alignments of a site that name the same variants with the same carriers (the best
  of them: the next entry). Measured
  on the search above, rerun genome-wide with this fix alone: before it, 238 of the 5,047 indel rows
  of the highest-CFD block started one base early and 4,467 showed a REF sequence that is not the
  reference at their start; now each of the 5,045 shows the reference at its start. On those rows
  `CFD_score_REF` changed on 785 (`PAM_creation` no longer reads this REF: the entry "An indel
  row was flagged as creating a PAM the reference already has"), and `CFD_score` went down on 218
  and up on 16 (chrX:42118270, 0.149 -> 0.017: the old value was the misaligned REF's, so the row
  left the report's CFD >= 0.05 table, as did 5 others). Four clusters now report the row of SNPs
  that has the same sequence, CFD and carriers as their indel row, and two indel off-targets
  that sat in a reference off-target's cluster because of their early start are rows of their own
  (chr3:55057499, CFD 0.087, 20 cis carriers; chr16:12067029, CFD 0.069, 1). The alternative
  alignments went from 722,564 to 722,347 rows, and 39 rows of the integrated results changed only
  their `Other_motifs` count. The report's tiers keep their sites, some one base later, except for
  those CFD >= 0.05 rows, the CFD >= 0.05 table gaining the two new rows and chr2:33157279 (whose
  REF now scores 0.1), and one tie at the top-1000 table's cut-off (CFD 0.188), now broken the other
  way. The SNP+indel companion keeps its 1,618 lines, one per alignment, the off-target start one
  base later on those rows; counted once each, as in the first entry and in the report, it went
  from 336 to 339 lines, all CONFIRMED, as some of those rows no longer share their line with
  another alignment of their site. The same search of the README's example guide,
  `CTAACAGTTGCTTTTATCAC` (the README search), keeps its 1,672 lines, 351 counted once each
  before and 363 now.
  Apart from `Other_motifs` counts and those four rows of SNPs, no reference or SNP-only row
  changes, nor do the perfect matches and the banner. The search takes as long as before, with the
  same peak memory. **Affected:** every search that allows RNA bulges (`--bRNA` 1 or more; the
  recommended default allows one) with a variant dataset that has indels, in every release, with or
  without the overlay and with `CRISPRME_INDEL_SNP=0` too; with a 5' PAM, every indel row with an
  RNA bulge. **What to do:** run such searches again. Covered by
  `PostProcess/test_indel_row_coordinates.py` (one and two RNA bulges, both bulge types, a DNA bulge
  and none, on both strands, with an NRG and a TTTV PAM, with `CRISPRME_INDEL_SNP=0` too) and the
  SNP+indel oracle, which now checks every indel row's REF sequence against the reference at its
  start.
- **Of the alignments of one site that name the same variants with the same carriers, the merge
  reported the first it read, not the best.** The merge of a cluster's rows lists such alignments of
  an indel row once and drops the others from both tables. It kept the first of them in the order
  its input is sorted, where alignments tied on score and on mismatches + bulges come in text order
  (an `RNA` alignment before an `RNA,DNA` one), so with the default `--sorting-criteria`
  (mismatches + bulges, then mismatches) the alignment reported could have one mismatch more and one
  bulge fewer than a dropped one, and in the fewest mismatches + bulges block another CFD. It now
  keeps the best of them by the order that picks each block's representative; among alignments
  that tie exactly there, the one that order's last tie-break puts first (the aligned sequence first
  in text order), whatever order they came in. Measured on the search above, this fix against the
  code before it: of its 182,367 integrated rows, 243 show another alignment of the same site,
  variants and carriers. In the highest-CFD block that is 56 rows, 12 with one mismatch fewer and
  one bulge more at the same CFD and 44 an alignment of equal rank; in the fewest mismatches +
  bulges block 232, 94 with one mismatch fewer and one bulge more (47 at another CFD) and 138 of
  equal rank (88 at another CFD, as that block's order does not look at CFD), and in 4 clusters the
  better indel alignment now comes before the cluster's SNP row, which moves unchanged to the
  alternative alignments; in the CRISPR-Bulge block 135 of equal rank. The alternative alignments
  keep their 722,347 rows, 191 of them changed in the same ways. Against the run before the
  previous entry's fix, no reported alignment has more mismatches: the 4 rows of the highest-CFD
  block and 18 of the fewest block that fix had moved to the worse alignment show the one reported
  before it again. The README search changes 323 of its 170,657 integrated rows the same way (of
  the 3 + 15 rows that fix had made worse, 16 show the alignment of before and 2 one of equal rank
  to it), and the search above at population level 241 of 211,830 (with no alternative alignments
  there, the 5 SNP rows that give the fewest block to a better indel alignment leave it). The
  report's tables keep their sites (one row of the <= 4 mismatches + bulges table now starts
  elsewhere); no carrier, rsID, AF, row count, reference or SNP-only row changes, the SNP+indel
  companion is the same, and the summaries and the raw observed-haplotypes sidecar follow the
  tables. **Affected:** every search with a variant dataset that has indels, in every release, with
  or without the overlay and with `CRISPRME_INDEL_SNP=0` too (until the previous entry's fix a +
  strand indel row with RNA bulges had a cluster position of its own and did not meet the site's
  other alignments). **What to do:** run such searches again.
  Covered by `PostProcess/test_representative_tiebreak.py`.

- **`--personal_annotation` reads a bgzip- or gzip-compressed BED, as its `--help` asked.** The
  help said the file "must be compressed using bgzip", but since 2.2.0 the personal BED has been
  read as plain text, so a bgzipped one reached the sort garbled and every search given one
  stopped before it began, with or without `--annotation` (`No tabs/spaces found at line 1 ...
  Error: Sorting BED file failed`, exit 1; the message never named compression). A compressed
  file, told by its content and not its name, is now decompressed first, and a damaged one stops
  the search with one line that says so. `--help` and `docs/INPUT_FORMATS.md` now say that a
  plain, bgzipped or gzipped file works. Measured through the real CLI on a toy genome: a
  bgzipped or gzipped personal BED, alone or with a plain or bgzipped `--annotation`, now gives
  the plain BED's results byte for byte, where it stopped at exit 1; with a plain BED the results
  are unchanged, byte for byte. Covered by `PostProcess/test_annotation_tmpdirs.py`.
- **`complete-search --help` and the docs said a search without the CRISPR-Bulge scorer env
  completes on CFD only.** Since 2.7.3 `complete-search` refuses to start while the `cbulge` env
  or its weights are missing or broken, naming the fix (`crisprme.py scorer-env create`;
  `crisprme.py scorer-env doctor` says what is wrong), and runs CFD-only, with the CRISPR-Bulge
  columns at `-1`, only under `CRISPRME_SKIP_SCORER_ENV=1`, as README §1.1.2 says. The
  `--compute-backend` help, its README entry, the README's scorer paragraph,
  `docs/TROUBLESHOOTING.md` (which added "rather than failing") and METHODS now say so, and that
  `targets-scoring` and a search started from the web interface do not check first: with the env
  absent they complete CFD-only, with the CRISPR-Bulge column at `-1` and a banner in the log.
  Nothing else changes. Covered by `PostProcess/test_complete_search_flags.py` (the help) and
  `PostProcess/test_annotation_tmpdirs.py` (the refusal).
- **The installers and the docs give the measured image and install sizes.** The one-line
  installers said the engine image is ~800 MB: `docker pull pinellolab/crisprme:latest`
  downloads 4.84 GB on amd64 (Windows, Intel Macs, most Linux hosts) and 1.81 GB on Apple
  Silicon, summed from the registry manifests, and unpacks to ~7 GB or more. They now say a
  ~5 GB download, ~1.8 GB on Apple Silicon. They also asked for a data folder with ~45 GB free,
  which is the download; the README, the website, `docs/DATA_LAYOUT.md` and the Quickstart said
  ~85 GB on disk and ~100 GB free, and `docs/TROUBLESHOOTING.md` ~3 + 40 + 50 GB. A default
  variant-aware install (`download --what all`, then
  `download --what index --index-name NRG_3_hg38+hg38_HGDP1kGP`) measures 97.9 GB, and peaks
  around 103 GB while the genotype store unpacks beside its archive, so ~100 GB free was full
  before the first result. Every one of them now says ~98 GB on disk and ~150 GB free, the
  Quickstart's breakdown lists the indel genotype store and the indel logs (14 GB it left out),
  and `install/README.md` names the figures to keep in sync. Covered by
  `PostProcess/test_documented_sizes.py`, a new step of `unit-tests.yml`.
- **The `lucapinello/crisprme-data` card said every index there works with older clients.** Its
  overview said so, while the same card, `docs/INDEX_PROVENANCE.md` and this changelog say a
  2.8.0 search of the rebuilt `NNN_3` archives fails or reports malformed rows. It now says any
  release can read them and to use the rebuilt `NNN_3` archives with 2.8.1 or later. The live
  card is corrected when the release uploads `docs/hf/crisprme-data_README.md`.
- **`docs/DYNAMIC_PAM_INDEX_COMPAT.md` read as a description of 2.8.1.** The page, new in 2.8.1,
  is the design written before `--prune-pam` was built, and described what does not ship: a
  `complete-search` that adds `--prune-pam` itself when the PAM is narrower than the index, a web
  PAM picker, an `index_meta.json` sidecar on every PAM-less index, and a `pinellolab` repository
  for PAM-less indexes. A status note under its title now says what 2.8.1 does: the flag is
  explicit, needs `--index-path` and is command-line only, and only `complete-test` passes it
  itself, to search a 3′ all-N index wider than 3 nt such as the K = 8 one; an index serves only
  its own PAM or, all-N, a narrower one; no `index_meta.json` is written or read (the PAM file
  is recorded in a `.pam_build` by `build-index-only` and by a search that builds an index);
  `NNN_3` is in `lucapinello/crisprme-data` and the 8-nt-flank indexes in
  `lucapinello/crisprme-data-general-index`. The sections that describe the rest are marked
  planned, and its example of the engine's subset rule no longer has `NRG` serve `NGA`.

- **History lists the jobs newest first again.** It sorted them by start time but read each row by
  the place the job had in the `Results/` folder listing, so the table showed them in that order,
  with the newest job anywhere in it. Each link opened the job named in its row. Since an upstream
  change of 2025 (2.8.0 too). Covered by `PostProcess/test_history_page_assembly.py`.
- **The caption under the results-page matrix said what VARIANT counts wrongly.** It said VARIANT
  counts the sites that exist only because a variant creates them (`Not_found_in_REF`), and
  REFERENCE the sites in the reference genome, even where a variant alters them. The matrix counts
  as VARIANT every site whose alignment reported for the selected filter criterion is a variant one
  (`REF/ALT_origin` alt): one a variant creates, and one a variant changes, for example to fewer
  mismatches. The shareable report's matrix is the one that counts only the first, so its VARIANT
  total can be smaller: on the walkthrough's per-sample search the page shows REFERENCE 12 /
  VARIANT 6 (3 of the 4 perfect matches under VARIANT) and the report 17 / 1; on a genome-wide
  HGDP + 1000 Genomes search (4 mismatches, 1 DNA and 1 RNA bulge), 105,217 / 84,732 against
  116,092 / 73,857. The caption now says so. The counts, the Download General Table file and the
  report are unchanged. Since 2.4.0. Covered by `PostProcess/test_web_results_links.py`.
- **The shareable report's Search mode row named how the mode was chosen, and was often wrong.** It
  read `Population-level (default)` in every population-level report and
  `Per-sample (--per-sample)` in every per-sample one. Per-sample is the default on a genotyped +
  phased index (such as the default HGDP + 1000 Genomes one), where population-level is the
  `--population-level` opt-out, and a web search passes no flag. A search records only the mode
  that ran, so the row now names that mode, `Population-level` or `Per-sample`, followed by the
  same text as before. The same labels in 2.8.0. Search results are unchanged; to correct a report
  built earlier, rebuild it with `crisprme.py generate-report --result-dir Results/<job>`. Covered
  by `PostProcess/test_generate_report.py`.

- **An indel row was flagged as creating a PAM the reference already has.** The reference bases the
  `PAM_creation` test read for an indel row were its REF alignment's, read from the row's start; on
  the + strand with a 3' PAM (on the - strand with a 5' one, such as Cas12a's) a deletion or an
  insertion in the protospacer moves the PAM `len(REF) - len(ALT)` bases from there, so the test
  read other bases, and a PAM the reference already has was flagged as created: on the WALK panel
  the rows chr2:10000 (`rsW_J_del`) and chr3:20000 (`rsW_N_del` + `rsW_N14`) showed `PAM_creation`
  TGG over a reference TGG, counted in the web's Summary by Mismatches/Bulges and Summary by Sample
  (S2, S4, S6) and on the Personal Risk Card (S4 "PAM creation 2"). The test now reads the reference
  bases the PAM sits on, through the indel; a PAM that holds an ALT base or spans a deletion's
  junction keeps the REF alignment's bases, as before, and so does a PAM on the other side, which no
  indel inside the target moves. **Affected:** indel rows of every variant search since
  `PAM_creation` was first set on them (2.7.10, #94), in 2.8.0 too; only that column and the counts
  built from it change (the report's `PAM_creation`, the summaries, the Personal Risk Card, the
  web), not the detection, scores, carriers or ranking. Measured genome-wide on the default HGDP +
  1000 Genomes index (per sample, 4 mismatches + 1 DNA + 1 RNA bulge), on the highest-CFD block: for
  guide `ACTGAAATCTGTAAGCAGGC` (NRG) the flag is no longer set on 1,954 indel rows and is set on 1
  more; rebuilt from the reference through each row's variants, every one of the 1,834 of them that
  could be checked lies over a reference PAM and the one more over none; of the rows still
  flagged that could be checked, 11 hold a PAM the reference lacks and 120 a PAM that holds a
  base of the indel, flagged as before, and the other 55 (186 in all) have their PAM across the
  indel's junction or inside a long insertion, also flagged as before. For the README's guide, `CTAACAGTTGCTTTTATCAC`, 2,655 no longer set (2,561
  checked, all over a reference PAM; of the rows left that could be checked, 19 + 150), which on its report changes that column, and no
  other, on 15 rows of the top-1000 table and 2 of the recommended panel; for TRAC,
  `CTCTCAGCTGGTACACGGCA`, 1,124 no longer set and 1 more (1,047 checked; of the rows left that could be checked, 10 + 89); for
  `ACTGAAATCTGTAAGCAGGC` on the NRR index, 2,761 and 2 more (2,585 checked; of the rows left that could be checked, 15 + 196). **What
  to do:** run such searches again. Covered by `PostProcess/test_indel_pam_creation.py` (the WALK J
  and N rows among them) and `PostProcess/test_indel_row_coordinates.py`.

- **A web search refused for an index a newer CRISPRitz built now says why on its status and
  results pages.** The search stops, as it should, with one ERROR (`<index> was built by a newer
  CRISPRitz (... is in index format 2, ...): this CRISPRme reads index format 1 only ...`, and what
  to do), but the job status page and its results link (`/load` and `/result`) knew only the other
  refusals' ERRORs, and said "The selected result encountered some errors, please remove it and try
  to submit again.": submitted again, the search was refused the same way, and the reason was shown
  only under Show full log. They now say "The search was not run: " and the ERROR, as they do for
  the other refusals. The search form gives no warning of such an index before Submit: which index
  a search reads is decided by the search script, and the form would need its own copy of that
  choice; after Submit, those pages say why. New in 2.8.1. Covered by
  `PostProcess/test_results_summary_and_status.py` and `PostProcess/test_web_results_links.py`.
- **The search form says before Submit that a variant dataset has no sample list, and how to get
  one, whether its variant index was built in Settings or downloaded.** Removing a VCF dataset in
  Settings takes the sample list the page wrote for it (`samplesIDs/<dataset>.samplesID.txt`) and
  leaves its variant index, and the search form still offers the dataset if its name begins with
  its genome's and an underscore (as `hg38_<name>`; one named otherwise is not offered until it is
  added again). The note under the bulge options then offered its search, and
  Submit refused it, saying to download its index again, which installs the list, or to add the
  file yourself: an index built in Settings cannot be downloaded, and one downloaded is not
  downloaded again while it is installed (Settings answers "already installed — nothing to
  download", and `crisprme.py download --what index` skips a complete install whatever its list).
  The note offered the search of any dataset Submit refuses for want of its list (one whose VCFs
  name samples but which has no list yet, too). It now says what Submit says, as it does of a
  sites-only dataset or of one added as a single multi-chromosome VCF file, and for a dataset
  whose VCF folder is not there the message says, if its index came from HuggingFace, to remove
  the index in Settings and download it again there, which installs the list; to add the VCF
  dataset again in Settings under its name, which writes the list; or to add the file yourself. A
  merged panel, whose list Submit makes from those of its datasets, is offered as before. New in
  2.8.1. Covered by `PostProcess/test_web_bulge_notes.py`.
- **Submit with a genome or PAM selected that is not installed says what is missing and where to
  add it, instead of failing.** On a fresh data folder, which holds no genome or PAM until some
  are added in Settings / Data Manager, Load Example selects `hg38` and `20bp-NRG-SpCas9`, and
  Submit's input check opened `PAMs/20bp-NRG-SpCas9.txt` (the default one, too, with no PAM
  selected), which is not there: the check failed with a server error (HTTP 500), no "Missing or
  wrong input" dialog opened, and the page did nothing. It read the genome folder of a genomic
  sequence given by BED coordinates the same way, so a genome not installed failed too, and so did
  coordinates on a sequence whose FASTA the genome folder does not hold, or past the sequence's
  end. The dialog now lists the genome and the PAM that are not installed, with where to add one
  in Settings (Add a reference genome, Define a nuclease / PAM), and the BED coordinates it cannot
  read; the bulge check, which needs both, waits for them. Pre-existing: 2.8.0 did the same.
  Covered by `PostProcess/test_web_bulge_notes.py`.
- **The web interface no longer prints a `BiopythonWarning` as it starts on a Bioconda install.**
  A stray `site-packages/pyproject.toml`, which conda-forge's `dash-bootstrap-components` ships,
  makes Biopython warn at import that it is imported from inside its source tree. `crisprme.py`
  has filtered that warning since 2.7.10, but `web-interface` starts `index.py` in a new
  interpreter, which inherits no filter and whose pages import Biopython. `index.py` now sets the
  same filter before its first import. The Docker image and `install_from_source.sh` remove the
  stray file, so they never showed it. Pre-existing. Covered by `PostProcess/test_app_settings.py`:
  a check of `index.py`'s code, which runs everywhere, and its imports run as `web-interface` runs
  them, with Biopython beside a stray `pyproject.toml`, which runs where the web stack and
  Biopython are installed: in the image, where `web-e2e.yml` runs it with the other tests that
  skip in the light job.

- **A failed `complete-search` ends on its summary of the failure, with exit status 1 and no
  Python traceback.** It printed where and why the run failed (`CRISPRme run failed.`, the stage,
  the last lines of `log_error.txt`, the log's path), then raised an `OSError`, so every failed
  search, the refusal of an index a newer CRISPRitz built among them, went on with `Traceback (most
  recent call last): ...` and an `OSError` line that only repeated the log's path.
  `assembly-search` added a traceback of its own when a haplotype's search failed; it now ends in
  one line, `Error: Haplotype search failed for --output <name>! See <log> for details`, which the
  web's assembly job reads as the job's failure, as it read the traceback. Exit statuses are
  unchanged, and a failed search still removes its annotation temp folders. Pre-existing, in 2.8.0
  too. Covered by `PostProcess/test_build_index_only_atomic.py`,
  `PostProcess/test_annotation_tmpdirs.py` and `PostProcess/test_assembly_haplotype_genomes.py`.
- **The documented self-test no longer reports a HuggingFace failure for files the data repo does
  not host.** `CRISPRME_ALLOW_ONDEMAND_BUILD=1 crisprme.py complete-test --chrom chr22
  --vcf_dataset 1000G` always printed `HuggingFace fetch failed for
  ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz (404 Client Error ...
  Entry Not Found for url ...); falling back to original source`, then downloaded the VCF from EBI
  and passed. `complete-test` tries the Hugging Face data repo first for the chromosome's FASTA,
  the VCFs and the sample lists it downloads (the annotations come from GitHub), and
  `lucapinello/crisprme-data` holds no raw 1000 Genomes or HGDP VCF (its `vcfs/` holds the
  source VCFs of the three shipped panels), nor `samplesIDs.HGDP.txt`, nor the unpacked FASTA of
  any chromosome but chr22, so `--vcf_dataset HGDP` and `--chrom chr21` printed it too. A file the
  Hub answers it does not hold is now said in one line on stdout, `<file> is not on HuggingFace
  (lucapinello/crisprme-data); downloading it from its original source` (`from UCSC` for a
  FASTA). Any other failure to fetch a file, such as no connection or a checksum that does not
  match, is still a warning on stderr, word for word as before. What is downloaded, and from
  where, is unchanged. Pre-existing, in 2.8.0 too. Covered by `PostProcess/test_complete_test.py`.
- **`setup` no longer reports a HuggingFace failure for files the data repo does not host
  either.** `crisprme.py setup`, which downloads the legacy database of the CRISPRme web server,
  tries the Hugging Face data repo first for the chromosome's FASTA, the 1000 Genomes and HGDP
  VCFs and the two sample lists, and printed `HuggingFace fetch failed for <file> (404 Client
  Error ... Entry Not Found for url ...); falling back to original source` (`to UCSC` for a FASTA)
  on stderr for each one the repo does not hold: every VCF, `samplesIDs.HGDP.txt` and the
  unpacked FASTA of any chromosome but chr22, so `setup --chrom chr22` printed it three times,
  and a whole-genome `setup` 47 times. It now says each in one line on stdout, as
  `complete-test` does: `<file> is not on HuggingFace (lucapinello/crisprme-data); downloading it
  from its original source` (`from UCSC` for a FASTA). Any other failure to fetch a file is still
  a warning on stderr, word for word as before. What is downloaded, and from where, is unchanged.
  The setup guide, `docs/crisprme_data_setup_051826.md`, said `setup` pulls the genome and the
  variants from their original sources, a single chromosome's FASTA from UCSC, and the website and
  `setup --help` that it downloads everything from them; all three now say `setup` takes two files
  from the data repo, `chr22.fa` for `--chrom chr22` and `samplesIDs.1000G.txt`, and everything
  else from UCSC, EBI, Sanger or GitHub. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_setup_legacy_database.py`, a new step of `unit-tests.yml`.
- **The `lucapinello/crisprme-data` card says what its `genomes/hg38/` holds, and lists its
  `guides/`.** Its layout, kept from the live card, gave the folder as a per-chromosome FASTA
  mirror, `chr1.fa … chrX.fa`. The repo holds one gzipped FASTA per hg38 sequence (`chr1.fa.gz` …
  `chrY.fa.gz`, `chrM.fa.gz` and the alt, random and unplaced contigs: 455 files, which
  `download --what genome` unpacks) and one unpacked FASTA, `chr22.fa`, which `complete-test` and
  `setup` fetch as is. The layout also left out the repo's `guides/` folder, which holds
  `sg1617.txt`, the example guide `CTAACAGTTGCTTTTATCACNNN` (no `download` fetches it), and the
  card itself, `README.md`. The layout now lists all of these, as the repo's file list gave them on
  2026-11-04. The live card is corrected when the release uploads
  `docs/hf/crisprme-data_README.md`.
- **`complete-test` run again over an earlier run's results names every folder in its way, and
  exits 2.** It named only the first such folder in sorted order, so after a benchmark-mode run it
  named `Results/crisprme-test-out_cas12a_hbg`, one of six, and a user who deleted it was stopped
  five more times; it stopped the prebuilt-index smoke on those six folders too, though the smoke
  writes only `crisprme-test-out`; and it exited 0, which a script or a CI step reads as a pass
  though nothing ran. It now names, in one line, every folder there of the run's own mode (the
  benchmark mode's `crisprme-test-out_<benchmark>`, the smoke's `crisprme-test-out`) and exits
  with status 2 before it downloads or searches anything; `crisprme.py` passes the 2 on, as for a
  refused `--chrom`. A folder of the other mode no longer stops a run, and the `CRISPRme complete
  test encountered an error` line of a failed run names only its own mode's error logs. README
  §2.2.2 says how to run a mode again. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_complete_test.py`.
- **A search's `log_verbose.txt` keeps every `Searching for INDELs in` and `Search ended for
  INDELs in` line.** `pool_search_indels.py` printed them from its pool's workers without flushing
  them, and leaving the pool terminates each worker it finds idle, with the lines it had not yet
  written: a search's log missed 3 to 23 of them, a different number each run.
  `pool_index_indels.py`, which builds a missing INDELS index for a search and for
  `build-index-only`, lost its `Indexing INDELs in` lines the same way. Both flush them now; no
  result changes. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_indel_empty_chromosome.py`.
- **A source install gets the web stack the Docker image and the Bioconda package run.**
  `environment.yml`, which `install_from_source.sh` installs into and which says it mirrors the
  Dockerfile's validated dependency set, left `flask`, `flask-caching`, `flask-compress` and
  `werkzeug` unpinned, and `cachelib` out. The Dockerfile and the Bioconda recipe pin them (`flask`
  3.1.3, `flask-caching` 2.4.1, `flask-compress` 1.24, `werkzeug` 3.1.8, `cachelib` 0.17.0), as
  an unpinned `flask-caching` once drifted to a release that crashed the web app at import; solved
  for linux-64 on 2026-11-03, the file took `flask-caching` 2.5.1, `flask-compress` 1.26 and
  `werkzeug` 3.1.9. It now pins the five as the Dockerfile does, and a dry-run solve of it for
  linux-64 gives the very builds of the release image. An environment created from the earlier
  file keeps the versions it got: create it again to get these. Pre-existing, in 2.8.0 too.
  Covered by `PostProcess/test_install_from_source.py`.

- **The `--help` of `complete-search` and `assembly-search` no longer says an annotation BED must be
  bgzip-compressed.** `--annotation` and `--gene_annotation` read a plain BED as well (both go
  through `_sort_annotation`, which has always taken either); the help said bgzip was
  required. `docs/DYNAMIC_PAM_INDEX_COMPAT.md` also no longer says that CRISPRme 2.8.0 can search
  the published `NNN_3` with no update: it could search the first archive (2026-10-25), and a 2.8.0
  search of the rebuilt one can fail or report malformed rows (see the `NNN_3` entry above).
  Text only.

- **Settings answers a personal-assembly bundle it refuses with the reason, as a refusal (HTTP
  400).** The upload took each refusal of a bundle for a failure of the server: a `metadata.json`
  naming a path outside the bundle and a FASTA file named with what a shell reads (both refused
  since 2.8.1), and an archive packed otherwise, with more than one top folder, holding a link or
  a device, or named after an assembly already installed, were answered HTTP 500 with
  `finalize failed: ...`, which the browser logs as a failed request, while the upload answers its
  other refusals 400 with the reason. The import now raises its own error for a bundle it refuses
  (`personal_assembly.BundleError`, a `ValueError`), which the upload answers 400 with the reason
  alone, shown on the page as before; any other failure is still answered 500. Pre-existing, since
  2.6.0. Covered by `PostProcess/test_settings_progress.py`.
- **`complete-search` refuses a PAM file it cannot read, and `--index-path` or
  `--vcf-filter-pass-values` with no value, in one line, before it does anything.** A PAM file
  named against the `<length>-<motif>-<CasName>.txt` convention (the nuclease's name is read from
  it), one holding two motif lines or a first line that is not `<motif> <length>`, and
  `--index-path` or `--vcf-filter-pass-values` given no value ended the search in a Python
  traceback (a `ValueError` nothing caught), once `Results/<output>` was made, which stayed,
  empty, and once the annotation was sorted (about 30 s with the shipped one). Each is now one
  `Error: ...` line and exit 1, as the other input errors are, checked before the annotation is
  read or any folder is made; an `--index-path` that names no folder, refused in one line before,
  no longer leaves the output folder either. `docs/INPUT_FORMATS.md`, which said that a misnamed
  PAM file fails with an `IndexError` and that a second line is silently ignored, says so.
  Pre-existing, since 2.2.0. Covered by `PostProcess/test_complete_search_flags.py`.
- **`complete-test` with nothing to search exits 2.** Its prebuilt-index smoke, which it runs
  without `CRISPRME_ALLOW_ONDEMAND_BUILD` (as in Docker when the variable is written in front of
  `docker run`, on the host), said how to download what was missing in a folder with no prebuilt
  index, or with one whose genome is not there, and exited 0, which a script or a CI step reads
  as a pass though nothing ran. It now exits with status 2, as when an earlier run's results are
  in its way, having said the same, and `crisprme.py` passes the 2 on. README §1.2.2 and §2.2.2
  and `docs/DOCKER_QUICKSTART.md` say so. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_complete_test.py`.
- **The report says what its ZIP and the result folder hold.** The default (slim) report said the
  full list of variant-contributed sites "is written to `variant_created.tsv` in the result
  directory (re-run with `--report-full` to bundle it into the archive)", but no run writes that
  file: the tier goes into the ZIP only, and only under `--report-full`. Its Downloads note said
  the bundled per-threshold subsets included the variant-created one, the methods note of its
  panel that the complete raw integrated results were among the bundled downloads, and the legend
  that the variants of a dense window were in the bundled `high_variant_density_regions.bed`: the
  slim ZIP holds none of these, and the full one bundles the windows near the top-ranked sites, as
  `high_variant_density_regions.bed.gz`. Each sentence now says what the slim ZIP holds, where the
  rest is in the result folder (the variant-created sites are the integrated results' rows whose
  `Not_found_in_REF` is `y`, a window's variants are in its `*high_variant_density_regions.bed`
  files) and what `--report-full` adds; the full report names what it bundles; a search with no
  variant-created site says it has none, and a wider CFD ≥ 0.05 net of no site is not called a
  bundled download. On a slim and a full report of a population-level and of a per-sample search,
  only these sentences change, and every other file of the ZIP is the same.
  `generate-report --help`, which said `integrated_results.tsv.gz` is bundled by default, and
  README §2.2.9 say what the slim ZIP holds too, and the website's example report
  (`docs/example-report.html`, a slim 2.8.0 report, which said all four) gets the four sentences
  as 2.8.1 renders them there, and nothing else. Pre-existing, since 2.7.10. Covered by
  `PostProcess/test_generate_report.py`.
- **A home that cannot be written in no longer makes every command print matplotlib's and
  fontconfig's complaints.** With a home CRISPRme cannot write in (`apptainer run --no-home`, or
  `--home` on a read-only folder; `docker run --read-only`), each command printed on stderr that
  matplotlib had made a temporary folder (`... is not a writable directory`,
  `Matplotlib created a temporary cache directory ...`), as `crisprme.py` imports it when it
  starts, and `Fontconfig error: No writable cache directories`, from the `fc-list` that
  matplotlib's font-cache build runs: 815 bytes with the release SIF on an HPC node with a
  read-only home, 635 with no home, from `complete-search` and `generate-report` alike. A search's
  steps write their stderr to its error log, and a step whose log is not empty fails, but they
  printed nothing: matplotlib had passed its folder on to them (`MPLCONFIGDIR`), so a per-sample
  and a population-level variant search of the synthetic genotyped dataset of the SNP+indel oracle
  ran to the end there, their error logs empty (`docker run --read-only` makes `/tmp` read-only
  too, and needs `--tmpfs /tmp` or `CRISPRME_TMPDIR`, as `docs/TROUBLESHOOTING.md` §1 now says).
  `crisprme.py` now gives matplotlib a temporary folder of its own before any command runs (in
  `CRISPRME_TMPDIR` when set, which it now applies before its imports too), removed when the command
  exits (one killed by a signal leaves it); when the home's `~/.cache` cannot be written, that
  folder is `XDG_CACHE_HOME` too, for fontconfig's cache and those of the other programs a command
  runs, but `HF_HOME` stays `~/.cache/huggingface`, where `hf auth login` saves a login, so
  `publish-index` still takes a login saved there and a download on a conda or source install still
  sends it. With a writable home, or `MPLCONFIGDIR` set, nothing changes. Run the same way with this
  change, the same searches and `generate-report` print nothing, with either home, and the searches'
  error logs stay empty. `docs/TROUBLESHOOTING.md` §1 (what goes to the temporary folder, and what
  stays) and the README's Apptainer note say so. Pre-existing, in 2.8.0 too. Covered by
  `PostProcess/test_configure_tmpdir.py`.

- **The report no longer says that a MAF of 1e-05 is a display floor.** The MAF note under the
  top-1000 tables and the legend's MAF entry said that a variant whose source allele frequency is
  0 is shown at a display floor of 1e-05 so that it renders on the log-scale plots, to be read as
  "present, frequency effectively 0"; the README's MAF column and METHODS sections 7 and 8 said
  so too. No release with these texts applied it: up to 2.3.3 `resultIntegrator.py` wrote 1e-05
  for an allele frequency of 0, and since 2.4.0, the first release with them, it writes `NA` for
  one that is 0 or missing, which the report leaves out of a row's MAF, the lowest frequency of
  the row's variants, blank when none of them has one; a frequency below 1e-05 is shown as it is,
  and the plots' log scale is their rank axis, a frequency only sizing a point. The note and the
  legend now say so, in that one sentence of each, in every report (the combined report of
  `assembly-search` too), as do the README, METHODS sections 7 and 8 and the website's example
  report (`docs/example-report.html`, a 2.8.0 report, which gets the sentence as 2.8.1 renders it
  there and nothing else). The rest of each report is unchanged, byte for byte, and no result
  changes. Covered by `PostProcess/test_generate_report.py`.
- **The report's MAF note no longer says that an indel's MAF is blank.** The note under the top-1000
  tables gave an indel-derived variant as one whose MAF is blank, the frequency registry being
  SNP-only. The registry is, but the indel post-analysis gives an indel the `AF` of the index's
  indel log, which CRISPRitz `add-variants` copied from the VCF's `INFO/AF`, in either mode and
  with or without the dict-less tiers, and the report shows it as it shows a SNP's frequency. In
  the self-test panel's `--per-sample` report with the tiers (the README's guide, `20bp-NGG-SpCas9`
  PAM, 4 mismatches, 1 DNA and 1 RNA bulge), the indel off-target at chr22:33985503 shows 3.00e-02
  and the indel-and-SNP one at chr22:21052642 1.80e-01, and the 3 of its 184 variant rows whose MAF
  is blank are indels whose `AF` reads 0 in the indel log (the VCF rounds it to two decimals); the
  183 indel rows of the website's example report all show one. The note now says that an indel has
  the `AF` of the index's indel log, and that its MAF is blank where the log gives that `AF` as 0,
  or not at all, in every report (the combined report of `assembly-search` too). The log gives
  none for a record that has `AF` at another place in `INFO` than the VCF's first record (the known
  limit on the place of `AF` under Changed), nor for any indel of the mega index, whose VCF has no
  `INFO/AF`. The README's MAF column, METHODS sections 3 and 7 and the example report
  (`docs/example-report.html`, which gets the note as 2.8.1 renders it there) say so too. Nothing
  else in a report changes, and no result changes. Pre-existing: the note has said so since 2.4.0,
  which wrote it. Covered by `PostProcess/test_generate_report.py`.
- **The report of a `--per-sample` search on an index without the dict-less tiers says what the
  search reported.** On such an index (one built without `--samplesID`, such as the one the
  benchmark `complete-test` builds), the report's Search mode row now says what
  `complete-search`'s banner says: the carrier samples are named, with the `AF` of the classic
  per-sample dictionaries (an indel's, of the index's indel log), which `add-variants` copied from
  the VCF, and no CONFIRMED/PUTATIVE cis phasing, SNP+SNP co-occurrence or exact joint allele
  frequency is reported, as they come from the tiers; a dense window can still be reported by its
  greedy worst-case alignment (the known limit on an index without the dict-less tiers, under
  Changed); and to get the tiers, rebuild the
  index with `build-index-only --samplesID` or download a prebuilt one. The texts that made the same
  claims elsewhere in the report now say what holds for those datasets: the Allele-frequency basis
  row and the MAF note under the top-1000 tables (a SNP's frequency is the `AF` of those
  dictionaries, not AC/AN, and one whose `AF` they give as 0, or not at all, has no frequency and is
  blank), the legend's MAF entry (no joint allele frequency is computed, and a variant whose `AF`
  the dictionaries or the indel log give as 0, or not at all, is left out of a row's minimum), the
  legend's Observed entry and the caption of the top-1000 tables (a dense window's greedy alignment
  lists the carriers of any of its variants, not its exact cis carriers). When no dataset of the
  search has the tiers, the Observed cell of such an alignment (a row of two or more SNPs flagged in
  `High_complexity_region`), in the report's tables and in every curated TSV that holds it, reads
  `N carriers (>=1 of K vars, poss. chimera)`, as in a population-level search, and links its site's
  other greedy alignments as "N other variant combinations". No column says which dataset a row
  comes from, and two datasets can name the same samples, so in a search that also has a dataset
  with the tiers, whose dense windows come as observed haplotypes flagged there too, every row keeps
  its cell, and the legend and the caption say why. The report read `.search_mode` alone, so every
  per-sample report said "Every per-sample haplotype is reported with CONFIRMED cis phasing, named
  carrier samples and exact joint allele frequency", AC/AN over the genotyped panel, an exact joint
  allele frequency and exact cis carriers. The self-test panel's search (the NGG index the README's
  self-test builds, `NGG_2_hg38_chr22+hg38_1000G`; the README's guide with the `20bp-NGG-SpCas9`
  PAM, 4 mismatches, 1 DNA and 1 RNA bulge) has no CONFIRMED line in any of its files; `Variant_MAF`
  is NA on 154 of its 189 alt rows (3 of 184 with the tiers), those whose variants' AF reads 0 in
  the VCF, which gives it to two decimals; four two-variant rows show the frequency of their
  commoner variant (0.09 to 0.86: chr22:16858220 shows 8.10e-01, from `Variant_MAF` "NA,0.81") for a
  haplotype that one to five individuals carry (9.81e-04 and 1.96e-04 with the tiers); and the 5
  greedy rows read "16, 40, 22, 18 and 18 cis carriers", three of them with a link to "N other
  observed haplotypes", though none of the 102 individuals they list (114 times) carries every
  variant of a row that lists it. The report cannot tell from the results, nor from the install,
  whose index can have the tiers by the time the report is built again, so `complete-search` now
  records the datasets its banner names in the result folder (`.untiered_datasets`), beside every
  dataset it searches (`.variant_datasets`), as it records `.search_mode`; a report built again from
  an earlier search's folder keeps the old texts. With the tiers this change leaves the report as it
  was, byte for byte (but for the order of tied rows, which the last entry below changes in every
  report), and no result changes. The `--per-sample` texts of `complete-search`'s help and of the
  README, the README's MAF column and METHODS sections 3 and 7 now say so too. The
  observed-haplotypes callout of Section 5 and the header of `observed_haplotypes.tsv`, as in every
  mode, and the download button of a greedy row's modal, which names the whole file, still call its
  rows observed haplotypes. The web search form offers per-sample only on an index with a genotype
  store, so a web search is never such a search. Pre-existing: the row has read so since 2.5.3,
  the allele-frequency texts since 2.4.0, the legend's entries since 2.6.0, the caption since
  2.7.10 and the Observed cell since 2.8.0. Covered by `PostProcess/test_generate_report.py` and
  `PostProcess/test_build_index_only_atomic.py`.
- **A population-level search on an index without the dict-less tiers says where its allele
  frequencies come from, and what `--per-sample` gives there.** On such an index both modes read
  the classic per-sample dictionaries, so a population-level search too reports each SNP with the
  `AF` that `add-variants` copied into them from the VCF: on the self-test panel above, 198 of its
  257 alt rows have a blank MAF. Its report said AC/AN over the panel all the same (the
  Allele-frequency basis row, the MAF note, the legend's MAF entry), as these texts did whatever the
  index. The Search mode row, as the banner, advised re-running with `--per-sample` for per-sample
  carriers, CONFIRMED cis phasing and exact joint allele frequency; the legend said that
  `--per-sample` computes the exact joint allele frequency and lists the exact cis carriers, and the
  caption of the top-1000 tables that it resolves the true cis haplotypes and carriers. On such a
  dataset `--per-sample` names the carriers only, with neither CONFIRMED cis phasing nor a joint
  allele frequency, and can still report a dense window by its greedy worst-case alignment.
  `complete-search` now records those datasets (`.untiered_datasets`, `.variant_datasets`) for any
  variant search, per-sample (the entry above) or population-level, and the web search job, which
  runs the search script itself, records them the same way. The population-level banner and the
  report's Search mode row say, after their `--per-sample` advice, what `--per-sample` gives on
  those datasets and how to get the tiers; the allele-frequency texts say what they say in a
  per-sample report; the legend's MAF and Observed entries and the caption say what `--per-sample`
  gives there. A sites-only dataset (one whose VCFs have no sample columns), to which a build
  without `--samplesID` writes no tiers either, is not one of those datasets, in the records, the
  banners or the report: with no carrier to name and no allele-frequency registry, its search
  reports none of its variants, in either mode. With the tiers this change leaves the banner and the
  report as they were, byte for byte (but for the order of tied rows, which the next entry changes
  in every report). Pre-existing: the allele-frequency texts have read so since 2.4.0, the banner
  and the Search mode row since 2.5.3, the legend's entries since 2.6.0 and the caption since 2.8.0.
  Covered by `PostProcess/test_generate_report.py`, `PostProcess/test_build_index_only_atomic.py`
  and `PostProcess/test_web_search_reuse.py`.
- **Rows of one score keep the table's order in the report's ranked tables and plots, and in the
  results page's plots.** The top-1000 by CRISPR-Bulge (`top1000_crispr_bulge.tsv` and its
  table), the top-1000 by CFD (`top1000.tsv` and its table) and the Section-2 plots' top-N were
  sorted with pandas' default sort, a quicksort, which is not stable: rows of exactly one score
  (common: identical aligned pairs at different sites score the same, and an out-of-domain
  CRISPR-Bulge score is -1) came out in an order that the whole table decided, so rows added
  elsewhere renumbered them; between two runs, 8 to 99 rows of `top1000_crispr_bulge.tsv` changed
  only their rank, each within an exact tie. These sorts are now stable, so tied rows keep the
  order of `integrated_results.tsv`. The self-test panel's `--per-sample` searches without and
  with the dict-less tiers (the searches above) differ by 5 rows (the known limit on an index
  without the dict-less tiers, under Changed):
  the 560 sites both reports hold came in another order at 282 places of
  `top1000_crispr_bulge.tsv` and 250 of `top1000.tsv`, and now come in the same order. On the
  search without the tiers, whose 565 rows hold 476 CRISPR-Bulge scores of -1 and 310 CFD scores
  of 0, the two tables put 477 and 483 rows at another rank than before, each inside its tie. Only
  the order of tied rows changes, and, where a tie spans the last row of a top-N, which of them it
  holds; no score changes, and no row leaves its tie. The recommended panel, sorted on seven keys
  at once (which pandas sorts stably), and the per-threshold tables, in the table's order, were
  already so and are unchanged. The results page's plots sorted the same way: the top-1000 plots of
  its Graphical Reports tab (by CFD, by CRISPR-Bulge, by fewest mismatches+bulges and by variant
  effect) and the Personal Risk Card's top-100 plots. They sort stably too now, so the top-1000
  plots put tied sites in the report's order: on the self-test panel's population-level search
  their rows of more than one mismatch+bulge, the rows the report keeps, came in another order than
  in the report's top-1000 tables at 547 (CFD) and 533 (CRISPR-Bulge) places, and now at none.
  Pre-existing since 2.4.0 (the report) and 2.2.0 (the results page's plots). Covered by
  `PostProcess/test_generate_report.py` and `PostProcess/test_top_plots.py`.

### Security
- **The web Settings page could be pointed outside `Settings/` and made to re-run a recorded
  command.** Retry, Show-full-log and the job-progress refresh built a path from the job id the
  browser sends back, so an id such as `../../x` reached outside `Settings/`, and Retry re-ran the
  `cmd.json` it found there through the shell. The page now accepts only ids of the exact form it
  issues and only directories inside `Settings/`, and Retry re-runs only a job its log shows as
  failed, never a maintainer-only publish while maintainer mode is off.
- **Genome archives are checked before they are extracted.** A genome `.tar.gz` uploaded in
  Settings, or fetched by `download --what genome --source url` (the Settings "Direct URL" source),
  was extracted with no check on its member names, so a crafted archive could write anywhere the app
  can (as root in the Docker image), including a job folder that Retry would then run. Every member
  is now checked first: an absolute or `..` path, a link, a device or a FIFO refuses the whole
  archive and nothing is written. Real genome archives (UCSC hg38 and mm39 chromFa) are unaffected.

- **The search form's genomic-sequence input passes nothing typed to a shell.** For a sequence
  given by BED coordinates, the sequence's `>` header and the sequence name on the coordinates'
  line went as typed into a shell command that ran `bedtools getfasta`: when Submit checked the
  form, as the base-editing window list was filled at each change, and when the search was
  launched. Characters a shell reads in either ran as commands, as the web app's user, which
  matters wherever the web app is served to other people. The header named the BED file written
  in the data folder, so a header holding a path wrote it elsewhere, and the index bedtools wrote
  beside the genome's FASTA was removed afterwards, one that was there before as well. bedtools
  now runs with no shell, on a BED file and a link to the FASTA in a private temporary folder,
  where it writes the index too: the header names no file, the sequence name must be that of a
  FASTA of the genome selected, and the genome an installed one. Pre-existing, in every earlier
  release. Covered by `PostProcess/test_web_shell_safety.py`.
- **Submit's launch of a search runs only what the form offers, and quotes each value it passes
  on.** A request can call the launch with any values, past Submit's check of the form: the
  genome, the PAM and the variant dataset name the folders and files the search reads, and the
  genome, the e-mail address and the number of mismatches went as sent into the shell command
  that starts the search, where characters a shell reads ran as commands. The launch now runs
  only a genome and a PAM that are installed, "Reference only" or a variant dataset the form
  offers for the genome, and whole numbers of mismatches and bulges, and it quotes each value of
  that command, which is unchanged for the values the form sends. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **A personal-assembly search takes only guides of the nucleotide codes a reference search keeps,
  and only a PAM and files the form offers.** Submit on the search form's Personal assembly tab
  wrote the guides of the PAM's length into the search's guide file as sent, upper-cased, and the
  search passes each guide of that file to a shell command: characters a shell reads in a guide
  ran as commands, as the web app's user, and a request can send that Submit any text. It now
  leaves out a guide holding a character other than those the Reference genome tab's Submit keeps
  of a guide (A, C, G, T and the IUPAC codes for two or three of them), as it leaves out one not
  of the PAM's length, rather than search a shorter one. A guide of those letters is searched as
  before; one holding an N is now left out too. That Submit searched, as sent, any PAM file and
  any folders and files for the individual's, outside the data folder too, and the search takes
  the PAM file's name, and that of a genome folder of `Genomes/`, into shell commands: it now
  searches only a PAM the form offers, with the files of the individual picked, as installed.
  Pre-existing, since 2.5.5. Covered by `PostProcess/test_web_shell_safety.py`.
- **The annotation a search applies is put together with no shell.** It is written as
  `Annotations/.active.<genome>.bed.gz`, and the name of a genome folder, which Settings took as
  typed and which can be put in `Genomes/` by hand, went into the commands that sort, compress,
  move and remove it, run through a shell, where characters a shell reads ran as commands. They
  now run with no shell. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **Settings refuses a name that CRISPRme would take into a shell command.** It took the name of
  a genome (typed, or that of the file uploaded), of a VCF dataset and of a nuclease with any
  character but a space, `+`, `/`, `\` and `..`, that of a file uploaded with any but `..`, and a
  genome download's URL as typed. CRISPRme takes these names into shell commands as they are (a
  search's post-analysis a genome's, a dataset's and a PAM file's, an index build a genome's and a
  dataset's, a genome download its own and the name its URL ends in), where characters a shell
  reads ran as commands, as the web app's user. Settings changes data in local mode only, which is
  how CRISPRme runs by default. It now refuses, and says why, a name that holds a blank, a
  character a sequence's name may not hold (`|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`, a
  backquote or `#`), `*`, `?` or `[`, or that starts with `-` or `.`: any name it takes, the name
  a URL ends in, and that of a file uploaded, but for a personal-assembly bundle, whose file's
  name only says how it is packed. Pre-existing. Covered by
  `PostProcess/test_settings_progress.py`.
- **Remove result, on a search's status page, removes only that search's folder.** It removed the
  folder the page's address named, through a shell: the job id there could run commands, or name
  a folder elsewhere, outside the data folder too, which was removed. It now removes, with no
  shell, only a folder of `Results/` whose id has the form the search form gives its jobs.
  Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **The results pages pass no value of a link or a table to a shell.** The guide, sample,
  chromosome and position a results link names, the guide and number of mismatches the results
  tables send back, a sample's name, as its VCF gives it, and the name of a search's results
  file, which holds the genome's, went into the shell commands that zip a table for download,
  draw the Graphical Reports and Personal Risk Card plots, and look up the targets at a
  position, where characters a shell reads ran as commands. Those commands now run with no
  shell, but for the look-up at a position, which now takes only a guide of the search, a
  position in digits and a chromosome name of letters, digits and `_.:,+=@%` (for another it
  says "The selected position does not exist"), and the results pages read only a folder of
  `Results/` whose id has the form the search form gives its jobs. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py` and `PostProcess/test_web_results_links.py`.
- **A job's status and results pages read, write and remove only that job's folders.** These
  pages take the job id from their address, which a client can set to any text, in local mode
  and on a website alike, and joined it as it was into the paths of the files they use: an id
  holding a path reached a folder elsewhere, outside the data folder too. Remove result of a
  personal-assembly search removed that folder, and the haplotype folders its `.Params.txt`
  names, joined the same way; Show full log, and the status of a personal-assembly search, read
  the log files there and showed them; the opening of a job's pages read its `.Params.txt`
  there; and the results pages wrote the table filter they keep (`.dropdown.json`) there, and
  their tables, sample and population lists and download links read the files there. Each now
  takes only a folder of `Results/` whose id has the form the search form gives its jobs, as
  Remove result of a search does (above), and otherwise changes nothing, shows no log, or says
  that there is no such job. History, which lists every search of `Results/`, shows one named
  otherwise (a command-line search's `--output` may be) with no link, and says why its pages do
  not open and what to do (`docs/crisprme_web_interface_user_guide.md`). Pre-existing. Covered
  by `PostProcess/test_web_shell_safety.py` and `PostProcess/test_history_page_assembly.py`.
- **Remove result removes only a search its status page offers it for: one that failed.** A
  search's status page offers Remove result once the search has failed, and a personal-assembly
  search's offered it once it had failed or was reconciled, but the button's request, which a
  client can send for any job, removed the search of the page's address whatever its state: a
  finished, running or queued search, and, through a personal-assembly search's button, a
  search of either kind. History lists every search in local mode, which is how CRISPRme runs by
  default, and on a website a search's link names it. Each now removes only a search in a state
  its page offers it for, and otherwise changes nothing, and a personal-assembly search's page
  offers it, as a search's does, once the search has failed: not once it is reconciled, which
  had anyone with a finished search's link remove it. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **Remove installed data, in Settings, removes only an item it lists.** It took the item as a
  request sent it, and a name of `..` named the data folder (or, through a folder of it that is a
  link to other storage, the folder there that holds the link's target), and `.` or no name the
  folder of every item of that kind: all of it was removed, every search's results with the data
  folder. It now measures and removes only an item the card lists, and says another is not
  found. Settings changes data in local mode only, which is how CRISPRme runs by default.
  Pre-existing, since 2.2.0. Covered by `PostProcess/test_settings_progress.py`.
- **A VCF upload in Settings goes only into a dataset uploaded or downloaded there, and takes a
  dataset's own files only.** An upload into a dataset registered as a server folder was written
  through the registration into that folder, under any file name, a new file or one replacing a
  file there: a request could write any file into any folder the web app's user can write. A
  registered folder is now read where it is and nothing is uploaded into it, and an upload takes
  `.vcf.gz` files and their `.tbi` or `.csi` index only (`docs/SETTINGS_DATA_MANAGER.md`).
  Pre-existing, since 2.2.0. Covered by `PostProcess/test_settings_progress.py`.
- **Settings publishes no index; an index is published from the command line.** With a Hugging
  Face token in the web app's environment (maintainer mode), Settings' Publish took any index
  name, and `publish-index` archives the folder that name leads to, beside it, and uploads the
  archive with that token to the shared repository, over a file of the same name there: a name
  holding a path had any folder of the server archived and uploaded. Taking only an index its
  list offers did not bound that, as a request has the list offer what it installs through the
  page: a VCF dataset of any server folder registered, or of a VCF uploaded, and its index
  built, whose archive holds the dataset's per-sample dictionaries and sample list. Settings now
  publishes nothing, and Retry runs no publish an earlier version started; in maintainer mode
  the page says to publish with `crisprme.py publish-index`, which is unchanged
  (`docs/SETTINGS_DATA_MANAGER.md`). Pre-existing, since 2.2.0. Covered by
  `PostProcess/test_settings_progress.py`.
- **Settings builds an index only of a genome, PAM and VCF dataset its card offers.** "Build an
  index locally" took them as a request sent them, and a path in one built an index of a folder
  outside the data folder. Pre-existing, since 2.2.0. Covered by
  `PostProcess/test_settings_progress.py`.
- **A personal-assembly bundle names its own files only.** The genome folder, chain and
  chromAlias file of each haplotype are paths its `metadata.json` gives, which an imported archive
  brings, and they were taken as they were: a bundle naming files and folders outside its own
  folder was imported, its individual offered as complete, and its search read them, its results
  holding what it found there. Such a bundle is now refused, and one imported before is listed
  as incomplete and not offered (`docs/crisprme_web_interface_user_guide.md`). Pre-existing,
  since 2.6.0. Covered by `PostProcess/test_personal_assembly.py` and
  `PostProcess/test_web_shell_safety.py`.
- **Remove of a personal assembly in Settings removes that bundle's folder, and of a pairing of
  the older layout its own markers only.** Settings listed a bundle under the individual its
  `metadata.json` names, which an archive imported there brings: one naming another individual
  took that one's place in the Individual picker and the Remove card, and one naming an
  individual with no folder was removed as a pairing of the older layout, by the files named
  after the assembly names the bundle gives, outside the data folder too. A bundle is now listed
  under the name of its folder, as `assembly-search --assembly-individual` takes it, and only a
  pairing of the older layout is removed as one (`docs/crisprme_web_interface_user_guide.md`).
  Pre-existing, since 2.6.0. Covered by `PostProcess/test_personal_assembly.py` and
  `PostProcess/test_settings_progress.py`.
- **The search form reads only what it offers.** The callbacks that fill in and check the form
  (the variant datasets offered for a genome, the analysis modes, the bulge options and the note
  under them, and Submit's check) took the genome, PAM and variant dataset a request sends into
  the paths of the folders and files they list and read: a path in them reached outside the data
  folder, and what was there changed what they answered, in local mode and on a website alike.
  They now answer for a value the form does not offer as for none, and Submit lists a variant
  dataset that is not installed for the genome, as it does a genome or PAM that is not installed.
  Submit's check of a dataset of one multi-chromosome VCF, new in 2.8.1, looks for its variant
  index with the PAM file of an installed PAM only, too. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **The indel steps of a search pass no chromosome name of a variant dataset to a shell.**
  A search with bulges takes each chromosome name of a variant dataset, read from the CHROM
  column of its VCFs (or from the per-chromosome indel logs named after it), into the shell
  commands that search and merge the indel off-targets and run the indel post-analysis, where
  characters a shell reads ran as commands, as the web app's user. A VCF dataset added in
  Settings carries its own CHROM column, which this did not check, and Settings changes data
  in local mode only, which is how CRISPRme runs by default. The indel search and its
  post-analysis now refuse a chromosome name that holds a blank or a character a sequence's
  name may not hold (the same characters Settings refuses in a name, above), and say so; a
  real chromosome name holds none, so a normal dataset is searched as before. Pre-existing.
  Covered by `PostProcess/test_indel_empty_chromosome.py`.
- **The post-analysis of a search passes no FASTA file's name to a shell, and a personal-assembly
  bundle brings only FASTA files named as a search can take them.** The post-analysis takes the
  name of each FASTA file of the genome folder searched, as its chromosome's, into shell
  commands as it is, where characters a shell reads ran as commands, as the web app's user. A
  genome added in Settings has its files named after their sequences, whose names are checked,
  but a personal-assembly bundle imported there brought its genome folders with the names its
  archive gives, which were not, and Settings changes data in local mode only, which is how
  CRISPRme runs by default. A bundle holding a FASTA file whose name holds a blank or a character
  a sequence's name may not hold (`|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`, a backquote or
  `#`) is now refused, and the import says why; the post-analysis refuses such a name too, and
  says so, and runs its step for each chromosome with no shell. A chromosome's off-targets are
  kept only from the file named after its sequence, whose name holds none of these, so a normal
  genome or bundle is searched as before. That step also looked a chromosome's off-targets up
  with its name as a shell word, where `*`, `?` and `[` (an HLA allele's name holds `*`) made a
  pattern of the files of the folder it runs in, outside the data folder, which were then read:
  it now looks the name up as it is. Pre-existing (bundles since 2.6.0). Covered by
  `PostProcess/test_personal_assembly.py`, `PostProcess/test_settings_progress.py` and
  `PostProcess/test_web_shell_safety.py`.
- **The results pages' tables bind the values of a link, a table and a form to their queries of
  a search's database.** The guide and sample of a results link, the guide and numbers of
  mismatches and bulges of a Show Targets link, the chromosome and positions typed in Query
  Genomic Regions, the guide of the summary table's row and a sample's name, as its VCF gives
  it, went as they were into the SQL that reads a search's results, and Custom Ranking wrote the
  columns it groups by, the order and the thresholds of its form into its SQL too: a quote in
  such a value ended the text it was written in, and the rest was read as SQL, which could widen
  the query, on the search's own database, or end it in an error. Each value is now bound to
  the query; Custom Ranking takes only the columns, orders and thresholds its form offers, and
  the positions of a region and the numbers of a link are whole numbers. A sample whose name
  holds a quote, whose table and Personal Risk Card failed, now has them. Pre-existing. Covered
  by `PostProcess/test_web_shell_safety.py`.
- **Two download links of the results pages look only for a file of their search's folder.**
  The Summary by Sample and Show Targets links offer their zip once it is there, under a name
  the page gives, which a client can set to any text, and they looked for it wherever that name
  led, outside the data folder too: whether a file was there showed in their answer, though
  nothing was served (the download route serves files of `Results/` only). They now look only in
  the search's folder, as that route would serve the file, and answer for any other name as for
  a zip not made yet. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **A saved e-mail password goes only to the mail server it was saved for.** Settings' Email
  notifications kept the saved password when a save named another server, port, SSL setting or
  sender address with the password field empty, as the page shows it, and the notification a
  search sends at its end logs in to the server saved, as the sender saved, with that password:
  a request, which can save any server in local mode (how CRISPRme runs by default), had the
  password sent to a server of its choosing. The password of `CRISPRME_SMTP_PSW` went the same
  way to a server the settings named with no password of their own. A save that changes any of
  the four now asks for the password again, and `CRISPRME_SMTP_PSW` is used with the default
  server (`smtp.gmail.com`) only (`docs/crisprme_web_interface_user_guide.md`). Pre-existing,
  since 2.2.0. Covered by `PostProcess/test_web_shell_safety.py`.
- **The archives the web app unpacks go in with safe modes, owned by the user unpacking
  them.** A genome `.tar.gz` uploaded in Settings or fetched from its Direct URL, a
  personal-assembly bundle, and a downloaded index, genotype store or build-inputs kit were
  unpacked with the mode bits and the owner each member of the archive gives: a setuid or setgid
  file, and, where the web app runs as root (the Docker image does), any owner the archive
  names, went into the data folder, which `docker run -v` mounts from the host, where a setuid
  file runs with its owner's rights. Each file now goes in with no setuid, setgid or sticky bit
  and no write for others than its owner, owned by the user unpacking it (tarfile's `data`
  filter, or the same done by hand on a Python before 3.11.4), and a member refused before (an
  absolute or `..` path, a link, a device, a FIFO) is still refused. A zip archive keeps no mode
  or owner when Python unpacks it, as before. Pre-existing. Covered by
  `PostProcess/test_settings_progress.py`, `PostProcess/test_genome_layout.py`,
  `PostProcess/test_personal_assembly.py` and `PostProcess/test_index_bundle.py`.
- **A search's e-mail notification goes to the one address typed, and links to the job by the
  server's own address.** Submit took any text holding an `@` as the address: a comma-separated
  list had one search notify every address on it, through the mail account saved in Settings,
  and a line break, which a request can send though the field holds one line, set the link of
  the e-mail sent. That link was made of the page's address as the request gave it, which a
  request sets to anything (and the one a browser gives lost its `//` there). Submit now takes
  one plain address, as `name@example.org`: no list, display name, comment, blank or line break,
  and no more than the 254 characters an address can have, counted before the address is read
  (the field is checked at each change of it, and reading a much longer value, which a request
  can send, took the server seconds). The field turns red for another one, and the search runs
  without a notification, as with no address. The link is the address the job's status page
  shows (in local mode `http://127.0.0.1:8080`) with the job's id. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **A job's link no longer gives the e-mail address its search was submitted with.** The
  `/Results/` route serves the files of a job's folder to anyone with the job's link, and in
  local mode, which is how CRISPRme runs by default, History lists every job. It served
  `email.txt`, which holds the address (and that of each search sent to the job as identical to
  it), and the search printed the address into its log (`MAIL: <address>`, then `EMAIL SENT to
  <address>`), which the status page's Show full log shows. The route now serves no `email.txt`,
  the search is not passed the address (its log says `MAIL: email.txt` for a notification,
  `MAIL: _` for none), and what the notification prints names no address, not even in a mail
  server's refusal. Notifications are sent as before; the log of a search run by an earlier
  version keeps its line. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **The Personal assembly tab's Submit takes whole numbers of mismatches and bulges, as its form
  offers them, before it writes anything.** It wrote the mismatches and bulges a request sends
  into the new search's `.Params.txt` as they were, a line break too, and only then converted
  them to numbers: a value that is not one left a job folder that never ran, listed as in queue,
  which no page removes, one for each such request; and a number the form does not offer was
  searched. It now takes 0 to 6 mismatches, the DNA and RNA bulges the dropdowns offer for the
  genome, PAM and variant dataset selected (0 to 2, more where an installed reference index is
  deeper, fewer where a variant dataset's index has fewer), and a Max edits of 1 to the
  slider's maximum where that governs (10, or the mismatches and bulges added where that is
  more), and otherwise says which is wrong, with nothing written. Pre-existing, since 2.5.5.
  Covered by `PostProcess/test_web_shell_safety.py`.
- **The Reference genome tab's Submit checks the Max edits value before it makes the job
  folder.** With Advanced options closed, where the Max edits slider governs, Submit converted
  the slider's value to a number only after it had made the job folder and written its first
  files: a value that is not a number, which a request can send, left a job folder that never
  ran, which History does not list, one for each such request. It now takes a whole number of
  the slider's range there, as the Personal assembly tab's Submit does, and runs nothing for
  another. Pre-existing, since 2.2.0. Covered by `PostProcess/test_web_shell_safety.py`.
- **What stays as it is.** A search's results TSVs, and their zips, hold the values they report
  as they are: an annotation label of a BED added in Settings, or a sample name or variant ID of
  a VCF dataset added there, that starts with `=`, `+`, `-` or `@` is a cell a spreadsheet may
  take for a formula. Escaping such cells would change the content of every TSV, so open them
  as text in a spreadsheet when that matters. Two features of Settings, which changes data in
  local mode only, are for whoever runs the server: Add a VCF dataset's "Register a server
  folder" registers any folder the web app can read as a dataset, and Add a reference genome's
  Direct URL has the server fetch any `http(s)` address. The online version of the web app
  disables both, with every other change Settings makes.

## [2.8.0] - 2026-10-24

### Added
- **Reference-sequence score columns + aligned spacer in the curated tables.** The top-1000 tables and
  every curated download TSV gain `CFD_ref` and `CRISPR-Bulge_ref` — the score of the REFERENCE sequence
  at each site, shown beside the ALT/worst-case score so a reviewer sees how much a variant moves each
  score — and an `Aligned_spacer+PAM` column placed immediately before `Aligned_protospacer+PAM`, so the
  guide and protospacer read as a pairwise alignment (bulge columns as `-`). `CRISPR-Bulge_ref` is gated
  like the CRISPR-Bulge column (shown only when CRISPR-Bulge was computed).
- **CRISPR-Bulge sensitive tier (benchmark 95%-recall operating point) in the Section-4 threshold
  table + a `crispr_bulge_gt_0.tsv` download.** The CRISPR-Bulge candidate-count table shows
  `≥ 0.5 / ≥ 0.2 / ≥ 1.9e-4 (95% recall)` — the orphan `≥ 0.1` tier is replaced by the **sensitive
  operating point** the recommended panel keys on: the threshold that recovers ~95% of validated
  off-targets on the 40-guide held-out benchmark (point 1.9e-4; guide-bootstrap 95% CI [5.8e-6, 3.1e-3];
  selectable via `CRISPRME_CB_RECALL_TARGET={90|95|99}`, single source of truth in `cb_operating_point`).
  The matching subset is bundled as a download, consistent with the risk-card sensitive row. The file
  key stays `crispr_bulge_gt_0.tsv` (stable download/web identifier).
- **CRISPR-Bulge scores are written at full precision (`{:.6g}`) instead of 3 decimals.** The sensitive
  default (1.9e-4) is below the old 3-decimal display floor (5e-4), which would have silently rounded
  ~7% of validated off-targets to `0.000` and dropped them from the sensitive tier (recall 0.93→1.00 on
  the benchmark). Full precision across the population/per-sample and indel scorers makes the operating
  point representable. CFD display is unchanged (3dp; its tiers are ≥0.2/≥0.5, far above the floor).

### Changed
- **CFD now penalizes bulges by position (Doench et al., 2016 gap matrix), instead of ignoring them.**
  CFD's score tables include position-weighted DNA-bulge (`r-:dX`) and RNA-bulge (`rX:d-`) factors, but
  the DNA-bulge factors were never applied: the scorer realigned a DNA-bulge target by trimming the 5'
  end, which collapsed the gap to a table position that has no entry, so every DNA bulge — and every
  bulge at the two undefined terminal columns — contributed a factor of 1.0 (no penalty). A bulge
  immediately adjacent to the PAM (maximally disruptive) therefore scored the same as no bulge
  (e.g. a chr19 site scored CFD 0.87 with a PAM-proximal DNA bulge, while CRISPR-Bulge scored it 0).
  CFD now scores the alignment as displayed — the extra genomic base of a DNA bulge is removed so the
  20 guide bases keep their true Cas9 positions (this also fixes a mispairing in the old 5'-trim), and
  the position-weighted gap factor is applied at the bulge's real position, with the two undefined
  terminal columns extrapolated (lenient at the PAM-distal 5' tip, abolished at the PAM-proximal
  junction). Mismatch-only sites and in-range RNA bulges are unchanged (bitwise); only previously
  unpenalized bulges move. The same position-aware scoring is applied across every CFD path (the
  population and per-sample analyzers, the indel analyzer, the two-pass worst-case screen, the
  standalone target scorer, and the factorized fast-mode scorer). A ≥2-bp bulge (opt-in `--bDNA`/
  `--bRNA ≥ 2`, outside the single-bulge calibration domain) applies the position factor per looped-out
  base; see METHODS §3. The deployment calibration (isotonic P(real), AUPR, precision/recall) was refit
  on the re-scored 2.94M-site benchmark and the report's PR curves regenerated, including a new
  bulge-stratified panel embedded in the risk-card modal. Measured effect: a per-site correctness fix,
  not a bulge-ranking gain — overall AUPR unchanged (0.168; the benchmark is mismatch-dominated) and
  CFD's bulge-only AP stays at the floor (≈0.01 vs CRISPR-Bulge ≈0.25), but removing the PAM-proximal
  over-call drops ~9,000 spurious high-CFD bulge calls at CFD≥0.2 (precision 0.006→0.007; at CFD≥0.5
  0.05→0.085). It reinforces the recommended-panel logic (CRISPR-Bulge ranks bulges; CFD rescues a
  single-bulge site only at high confidence) rather than changing it.
- **Context-aware default analysis mode: per-sample on a genotyped+phased index.** With no mode flag,
  `complete-search` now auto-selects **per-sample** (CONFIRMED cis phasing + named carrier samples +
  exact joint AF) when the variant index is genotyped **and phased** (e.g. the default HGDP+1000G
  panel), and **population-level** otherwise (sites-only / unphased / hybrid / reference-only). The
  index `data_type` is read from the Phase-2.1 registry manifest. New **`--population-level`** flag
  forces the fast worst-case screen (opt out of per-sample on a phased index — recommended on a very
  dense / high-edit / genome-wide search); `--per-sample` still forces resolution; the two are mutually
  exclusive. The web Analysis-mode selector mirrors this (defaults to per-sample on a phased panel). A
  one-line notice prints the auto-selected mode and how to override. Rationale: phased panels are built
  for exactly this resolution, so the richest, most honest result is the default where it is tractable.
- **The recommended validation panel (report Section 4) selects on a principled sensitive operating
  point (benchmark 95%-recall), not an ad-hoc `> 0`.** The panel hard-includes
  `(CRISPR-Bulge ≥ 1.9e-4) ∪ (mm+b ≤ 2) ∪ (CFD ≥ 0.2 at 0 bulges) ∪ (CFD ≥ 0.5 at 1 bulge) ∪
  (clean 0-mismatch bulge) ∪ must-includes` — making the shipped `panel_top100.tsv` *be* the recommended
  default panel rather than a CFD-centric subset. The CRISPR-Bulge sensitive threshold is the **95%-recall
  operating point** on the 40-guide held-out benchmark (CB cannot rank bulges: bulge=1 AP ≈ 0.02 vs ≈ 0.25,
  so it is the primary axis); the CFD rescue is sensitive off-bulge (≥ 0.2) but requires high confidence
  (≥ 0.5) on single-bulge sites where CFD is an extrapolation; ≥2-bulge CFD-only sites are kept-and-demoted
  (uncalibrated for both scores), never dropped. One consistent threshold message now runs across the risk
  card, the "how to read this" modal, the Section-4 methods note, the score glossary, and METHODS §3.5.
  Why a principled cutoff and not literal `> 0`: CRISPR-Bulge assigns a tiny non-zero score to ~every site
  (99.9% of benchmark negatives fall below 1.9e-4), so `> 0` would list essentially the whole genome.
  95% recall is the sensible trade — for the ACTGA/NRG worked example (189,957 candidate loci) it yields
  162 off-target validation loci (90% → 58; 99% → 1,158; the shipped crispr_bulge_gt_0.tsv tier); pushing from 95% to 99% recall is ~7× the sites because the
  last few percent of true off-targets sit in the score's noise floor. Framed as a risk/validation-size optimization, not a
  recall guarantee; recall above ~95% is reachable (selectable) but the list grows steeply.

### Fixed
- **Reference-only search on a read-only install no longer dies on the `vuoto.txt` placeholder.**
  The serial annotation-index pre-build (the flock-race fix) called `load_annotation_bed()` directly
  on the annotation file, bypassing the empty-annotation guard; for a reference-only search that file
  is the empty `vuoto.txt` placeholder, which `pysam.tabix_index()` bgzips in place and removes —
  consuming it, so the later gene-annotation step's `cp vuoto.txt` failed and aborted the search (and
  the write fails outright on a truly read-only install / Apptainer SIF). `load_annotation_bed` now
  returns `None` for `vuoto.txt` (never indexes the placeholder) and the pre-build skips it
  (`test_annotation.test_vuoto_placeholder_never_indexed_or_consumed`; release trap T-38).
- **The "Observed" cell no longer reports a misleading raw carrier count for population-level,
  multi-variant off-targets.** A dense-window site that stacks N variants into one worst-case
  protospacer previously showed e.g. `128 carriers`, which reads as "128 people carry this
  off-target" — but in the population-level (unphased) default that number is the **union** of
  samples carrying ≥1 of the N variants, and the stacked combination may exist on no single
  haplotype (a possible chimera; the row's own low joint AF gives it away). It now reads
  `128 carriers (>=1 of 12 vars, poss. chimera)`. A single-variant site still reads `N carriers`, and
  under `--per-sample` (phased, exact cis set) it reads `N cis carriers`. Applies to the in-report
  table and the download TSVs (one shared cell).
- **The per-site "Observed" modal/caption no longer calls population-level reconstructions
  "observed haplotypes."** In the population-level (unphased) default a site reached by several
  variant combinations now reads `N variant combinations (worst-case)` (modal), `N other variant
  combinations` (link), and the caption spells out that these are reconstructed from the observed
  variants and are unphased — "run `--per-sample` to resolve the true cis haplotypes and carriers."
  Under `--per-sample` the wording stays `observed haplotypes` (genuinely resolved/phased). Mirrors
  the carrier-count fix one level up.
- **`--per-sample` searches no longer fail at the annotation stage on Docker Desktop for Mac/Windows.**
  The post-search step annotates three target files in parallel, each building the same tabix index;
  `flock` is non-functional on VM-backed bind mounts (virtiofs/9p → ENOTTY), so the concurrent builds
  raced, corrupted the index read, and tripped the stderr-fatal gate — killing the whole search at the
  final stage with a misleading "file may be corrupt, re-download". Fix: pre-build the tabix index once,
  serially, before the parallel annotators (they then only read a complete index), and tolerate a
  `flock` that errors on such filesystems. Linux/cluster (working `flock`) was unaffected.
- **Transient truncated reads of the annotation no longer abort a search.** `load_annotation_bed` now
  retries the BGZF-EOF guard + index build + open with a short backoff (filesystem read-after-write lag
  on virtiofs/NFS), declaring the file corrupt only after the view fails to settle (+ 4 regression tests,
  now wired into CI).

### Documentation
- Troubleshooting entries for two CLI first-run messages a newcomer can hit: a guide needing `NNN`
  PAM-padding (23 nt), and the "output folder not empty" a validation-aborted run leaves behind.

### Removed
- Retired the standalone 1000G-2021 index (`NRG_3_hg38+hg38_1000G2021` + genotype companion) —
  superseded by the default HGDP1kGP (which includes the jointly re-phased 1000G-2021 samples);
  deleted from HuggingFace (26 GB reclaimed), rebuildable from the retained `hg38_1000G_2021`
  source. The `hg38_1000G_2021` source dataset remains (a mega + HGDP1kGP input).

## [2.7.10] - 2026-10-22

### Added
- **Top-of-report off-target risk summary (replaces the aggregated specificity score).** The report
  now opens with a per-scorer risk card: for CFD and for CRISPR-Bulge separately (never blended), the
  count of ≤1-gap off-target sites in each score band and a **calibrated expected number of true
  off-targets** with a coverage-valid 90% range, plus a qualitative note that a deployed guide's own
  off-target propensity (unknown, and not given a specific multiplier) can shift the count further.
  Estimates come from
  an isotonic score-to-probability calibration fit leave-one-guide-out on 359 experimentally confirmed
  off-targets across 46 guides (Tsai 2015, Refined TrueOT, Chen 2017, Listgarten 2018; Zenodo
  10.5281/zenodo.10902445); held-out AUPR 0.168 CFD / 0.613 CRISPR-Bulge; the range is a two-level
  guide-cluster bootstrap with a negative-binomial between-guide overdispersion term, validated
  per-held-out-guide to about 90% coverage. A CSS-only "How to read this" pop-up (no JavaScript, so the
  report stays self-contained) explains the two scorers and shows the precision-recall curves. New
  `PostProcess/offtarget_calibration.py`, `PostProcess/risk_card.py`, and `PostProcess/calibration/`
  (bundled isotonic maps + bootstrap ensemble + per-band overdispersion). The recommended validation
  panel reconciles with the card: its CFD tiers use the SAME cut points (≥ 0.5, ≥ 0.2), each showing
  the all-gaps assay count with the ≤1-gap calibrated subset (the card's domain) in parentheses, plus
  a bridge note; the wider CFD ≥ 0.05 net stays a bundled download referenced in the note. The MM/B
  matrix legend moved to a full-width band below the summary grid (flowing in columns) for a more
  compact page. Full method, provenance, and the statistician/clinician/genome-editor sign-off in
  `docs/offtarget_risk_summary_METHODS.md`.

- **Own implementation of the CFD score** (Doench et al., *Nat. Biotechnol.* 2016) in
  `PostProcess/cfd_score.py`, with the mismatch and PAM penalty tables as reviewable JSON
  (`cfd_mismatch_scores.json`, 392 keys; `cfd_pam_scores.json`, 16 keys) and a committed generator
  (`scripts/generate_cfd_tables.py`, with a `--check` table-validation mode). It is now **the**
  production CFD scorer: the SNP, INDEL, fast-mode and `targets-scoring` paths all score through
  `cfd_score.cfd` via the position-aware `bulge_cfd_penalty.cfd_score_bulge_aware`, and the previous
  per-module scorer (three `calc_cfd` copies + their helpers) and the binary `.pkl` penalty tables
  have been removed. Output is **bit-for-bit identical** to prior releases' CFD: a committed golden
  fixture (`cfd_golden.tsv.gz`, 6,132 cases spanning mismatches, DNA/RNA/terminal/PAM-proximal bulges,
  every gap key, IUPAC and non-canonical PAM) is reproduced exactly (raw IEEE-754) by
  `test_cfd_cleanroom_equiv`, which also checks clean-room internal consistency and the published
  Doench anchor values.

- **`crisprme.py targets-scoring` — score a supplied LIST of off-target sites** (no genome search).
  Takes Cas-OFFinder output (`--casoffinder`, schema auto-detected) or a generic TSV with an explicit
  column map (`--columns guide=..,dna=..,chrom=..,pos=..`) and annotates each site with the **same**
  CFD and CRISPR-Bulge calculations a normal CRISPRme+ report uses (the shared clean-room
  `cfd_score.cfd` via `cfd_score_bulge_aware`; CRISPR-Bulge via the `cbulge` env, `-1.0`
  when absent). With `--index-path` it adds a per-site **variant carrier overlay**: the variants
  overlapping each site's protospacer+PAM window and, per population dataset, their allele frequency +
  carrier count (read from the index's `Dictionaries/registry_<vcf>` tier). Enables scoring
  externally-nominated sites (e.g. non-canonical-PAM hits) without building a new index. New module
  `PostProcess/targets_scoring.py` + `PostProcess/test_targets_scoring.py` (wired into unit-tests CI).

- **Full-report opt-in surfaced on every entry point.** The slim report (above) can now be upgraded to
  the complete self-contained bundle from anywhere: `complete-search --report-full` makes the
  automatically-generated end-of-search report full (threaded via `CRISPRME_REPORT_FULL` → the shell's
  `generate_report.py --report-full`), the standalone `generate-report` CLI gained the `--report-full`
  flag (its `main()` now parses and threads it — previously only the Python `build_report(report_full=…)`
  API honored it), and the **web submission page** gained a *Report options → "Include the full report"*
  checkbox (present on every search tab, including Personal assembly) that sets the same env var. New
  `test_cli_main_report_full_flag_threads_to_bundle` covers the CLI→bundle path.

- **`observed_haplotypes.tsv` report companion — every in-budget haplotype per multi-haplotype site.**
  The per-cluster representative in `integrated_results` necessarily collapses a dense/multi-haplotype
  site to ONE headline row; this bundled companion lists **each** distinct observed variant haplotype
  at any site with ≥2 of them — its variant set, rsIDs, mm/bulges, **both** scores (CFD **and**
  CRISPR-Bulge), origin, AF, and named carriers — so a reviewer sees the full set a site contributes
  (serving the CFD-ranked **and** the CRISPR-Bulge-ranked views). Report-layer only (built in
  `generate_report` from the representatives + the alt-alignments file, merged per `(site, variant-set)`
  across axes); bundled in the SLIM report. New `PostProcess/observed_haplotypes_companion.py`.

### Changed
- **Performance: observed-haplotypes sidecar cuts the per-sample report render ~35%.** The report's
  clickable observed-haplotype modals were built by re-scanning the whole (30+ GB) alt-alignments
  table every render, even though only the alt-origin rows matter. The search now pre-reduces that
  table once (`make_oh_sidecar.py`, new STEP 8a) into a small `*_observed_haplotypes_raw.tsv` sidecar
  (alt-origin rows only — a strict superset of what the report keeps), and the report reads the sidecar
  instead. Measured on a 6,1,1 `--per-sample` run: 51.3M → 2.93M rows (31.5 GB → 8.9 GB), render
  ~25 → ~17 min, and `observed_haplotypes.tsv` is **byte-identical** (the report re-applies its exact
  per-row filter + clustering) with the prioritized top-1000 unchanged. A small zipped sidecar is kept
  after cleanup so standalone `generate-report` re-renders stay fast; older result dirs with no sidecar
  fall back to the full alt scan. New `test_make_oh_sidecar.py` + a byte-identity guard.

- **Recommended validation panel: domain-aware CFD hard-include with ≥2-bulge demotion.** The panel
  now hard-includes a site on CFD only within CFD's calibrated domain (CFD ≥ 0.5 **at ≤1 bulge**). A
  site that needs **≥2 bulges** and is high on CFD alone — where CFD is an uncalibrated over-call and
  CRISPR-Bulge is N/A — is **kept in the panel but *demoted*** (ranked after every reliable
  hard-include) rather than forced to the top on an untrustworthy score, and it is never dropped: a
  ≥2-bulge off-target can still be real (the benchmark contains no ≥2-bulge examples, so its risk is
  *unknown*, not zero) and only the user can decide what to assay. The count of such demoted sites is
  surfaced in the panel note even on a CFD-only run; the unconditional pure-bulge hard-include
  (0 mismatches, ≥1 bulge) is unchanged. The Section-7 CFD glossary tooltip and the in-report Methods
  note were updated to match, with regression tests pinning both the saturated-regime demote behaviour
  and the glossary wording.

- **`complete-search` now cleans up bulky intermediate files by default** (new `--keep-intermediates`
  to retain them). Once the finalized `*_integrated_results.tsv` and the report are produced, the
  pre-integration working files — the merged `bestMerge`/`altMerge` tables, the raw
  `all_results_with_alternative_alignments.tsv` (the `.zip` is kept), and the per-chromosome
  population summaries — are removed, reclaiming ~15 GB on a genome-wide run (measured on a 5+1+1
  `--per-sample` run). **Kept either way:** `integrated_results.tsv`/`.zip`, the `report.zip`, the
  alt-alignments `.zip`, the per-chromosome co-occurrence / high-variant-density sidecars (so a later
  `generate-report --report-full` can still rebuild those bundles), the per-guide/per-sample
  summaries, and the sqlite `.db` web results-page index (removing it would break browsing the job in
  the web interface). Gated by `CRISPRME_KEEP_INTERMEDIATES`.

- **`assembly-search` coordinate lifting now uses the `pyliftover` library** instead of an external
  liftOver binary. `pyliftover` replaces `ucsc-liftover` in `environment.yml`/`Dockerfile`, so the
  image no longer carries that binary and `assembly-search` works out of the box. Point-lifting
  behaviour is unchanged (per-endpoint target coordinate + strand; a point covered by more than one
  chain, or none, is left unmapped); chain files are still passed explicitly via `--chain-paternal` /
  `--chain-maternal` (nothing is fetched). Added a `THIRD_PARTY.md` attribution file for bundled
  components (CRISPR-Bulge, impg, pyliftover — all MIT).

- **The shareable report (`report.zip`) is now SLIM by default.** It bundles the human-facing layer
  — `report.html`, the curated top-1000 (both scorers), the validation panel, the per-tier CFD/
  CRISPR-Bulge/mmb subsets, the indel+SNP co-occurrence, and the run manifest/Params — which stays a
  few MB and easy to share. The bulky full-data files (the complete raw `integrated_results.tsv.gz`,
  the full `snp_snp_cooc`, `variant_created`, and `high_variant_density_regions`) are now bundled ONLY
  with the new **`--report-full`** flag (gzipped when bundled — e.g. a ~300 MB `snp_snp_cooc` → ~23 MB);
  they always remain in the result directory regardless, so nothing is lost. On a real genome-wide
  5+1+1 `--per-sample` run this took the default report from **274 MB → ~5 MB**; `--report-full`
  reproduces the complete bundle. The slim report degrades its HTML gracefully (no download links to
  the excluded files). `crisprme.py generate-report --report-full` opts into the full bundle.

- **The slim report now size-gates the observed-haplotypes companion.** The
  `observed_haplotypes.tsv` sidecar (one row per distinct observed variant haplotype at every
  multi-haplotype off-target site) is tiny on a normal run and keeps shipping in the slim bundle; on an
  extreme dense search it can grow to ~GB and would otherwise bloat the slim/shareable report. It is
  now bundled in the slim report only while it stays under a size threshold (default ~25 MB, overridable
  via `CRISPRME_OH_SLIM_MAX_MB`); over the threshold it moves to `--report-full`, exactly like the other
  bulky full-data files. The full file is always written to the result directory regardless, so nothing
  is lost. On a 6+1+1 `--per-sample` genome-wide run this takes the slim report from ~442 MB to ~20 MB.
  The clickable per-site top-1000 modals (both scorers, every haplotype for each shown site) are
  unchanged either way; only the top-level download link/callout/note and the modal's download button
  degrade (the per-site haplotype table and a short "available via `--report-full` / in the result
  directory" note remain) when the companion is not in the slim bundle.

- **Web Contacts page:** added Ann Cirincione under the Pinello Lab, and reworded "Please send any
  comments to Luca Pinello only" → "…to Luca Pinello" (dropped "only").

- **Web form hides the "Analysis mode" selector for a personal-assembly search.** Population-level vs
  per-sample resolves genotypes across a multi-sample variant *panel*; a personal-assembly search
  works on one individual's own two haplotypes (no panel to resolve), so the control is irrelevant
  there. The tab-toggle callback now hides the Analysis-mode block when the Personal-assembly tab is
  active and restores it for the reference-genome tab.

- **Web form shows an annotation note.** The per-search annotation selector was removed earlier, which
  obscured that searches are still annotated by default. The form now carries a short note that results
  are annotated with the functional regions enabled in Settings → Annotations (DHS, ENCODE SCREEN cCREs,
  GENCODE, COSMIC by default) plus the closest gene, and that tracks can be toggled — or a custom
  annotation added — there. (No new control; discoverability only.)

- **Install-walkthrough corrections (#92, Ann Cirincione).** Cluster SSH port-forwarding guidance
  (the server binds `0.0.0.0` with no auth; both tunnel forms given inline); per-OS Docker memory
  (macOS slider / Windows `.wslconfig` / Linux no cap); and a **backwards `--max-total-edits` rule** —
  the docs said it "must be ≥ mm+bDNA+bRNA," which is actually how you *disable* the cap (default is 6
  = 4mm+1DNA+1RNA; setting it below the sum prunes bulged sites, by design). Plus reconciled the
  reference-only disk figure (~25 GB), de-hardcoded a `:v2.7.9` pin, and refreshed stale
  `1000G2021` → `HGDP1kGP` references.

### Removed
- **Removed the stale "odd-length degenerate PAM + bulges can crash the search engine (issue #105)"
  runtime warning.** It fired on the *default* PAM (`NRG` is odd-length and carries the
  partially-degenerate `R`), so it cried wolf on essentially every default search — about a CRISPRitz
  C++ crash that was fixed upstream. CRISPRme+ hard-pins CRISPRitz 2.9.0 on every install surface
  (Docker `crispritz_ref=v2.9.0`, `install_from_source.sh`, Bioconda `crispritz >=2.9.0`), so no
  supported install can run a vulnerable engine. The `docs/` note is reworded to historical.

### Fixed
- **Observed-haplotypes companion no longer blows up memory on dense per-sample searches.**
  The build loaded the entire `all_results_with_alternative_alignments.tsv` (30+ GB on a
  dense run) and materialized every row, exploding to ~140 GB. It now streams the file in
  200k-row chunks and, for the slim report, keeps only rows near the top-1000 displayed
  sites (whose per-site modals are the only consumers of the haplotype map), with a
  vectorized pre-filter so the giant file is scanned, not materialized. Peak memory drops to
  tens of GB and the pass is ~2.5x faster; the full all-sites companion is built only under
  `--report-full`. Output is byte-identical to the unscoped build for the displayed sites.

- **A samplesID file with a trailing blank line (or any blank/comment/short row) no longer aborts the
  run.** (Reported independently by @ManuelTgn, PR #105, from a real 4091-sample HGDP+1kGP run.)
  The per-sample summary and specificity-radar steps (`process_summaries.py`,
  `radar_chart_dict_generator.py`) indexed samplesID columns without guarding against short rows, so an
  ordinary trailing newline in the `--samplesID` file crashed the job with `IndexError: list index out
  of range` at the "Creating images" stage — **after** the full multi-hour search had completed. The
  samplesID consumers now skip blank, `#`-comment, and short rows; `submit_job` strips blank lines when
  building the working `.sampleID.txt` (central fix for every consumer); and the pre-flight validator
  now checks the resolved samplesID association files up front and WARNs (never errors) if a *data* row
  has fewer than the expected 4 columns (`SAMPLE_ID`, `POPULATION_ID`, `SUPERPOPULATION_ID`, `SEX`),
  naming the file and an example line so a genuinely malformed samplesID is caught before the search
  rather than hours into it. +4 `test_validate_inputs` cases.

- **Indel off-targets now carry their cis-phased co-occurring SNP(s) in the `Variant` column.** When an
  indel off-target's true haplotype is an indel **plus** one or more co-occurring, cis-phased SNPs, the
  `Variant` column (and the `observed_haplotypes` companion, which reads the same field) previously
  showed the indel **only** — the co-occurring SNP(s) were recorded solely in the per-chromosome
  `*.indel_snp_cooc.tsv` sidecar, never merged into the off-target row. The INDEL post-analysis now
  appends the co-SNP forward-strand variant tokens (`chrom_pos_ref_alt`, matching the SNP-path format so
  `resultIntegrator` parses and strand-corrects them identically), with the `rsID` and `MAF` lists kept
  1:1 token-aligned, and narrows the `Samples` column to the **exact cis carriers** (individuals carrying
  the indel AND every required SNP on the same phased slot). **Gated to CONFIRMED genotype-phased
  co-occurrence only** — the gate is the co-occurrence *phase* (derived from the index's on-disk phased
  genotype stores), **not** the search mode, so the enrichment applies equally to a population-level and
  a `--per-sample` search on a genotyped-phased index, and cannot fire on a sites-only index (no
  genotype store → that row is byte-identical to before). A PUTATIVE (unphased / sites-only / mega) call
  is left sidecar-only and does not touch the row, preserving the anti-chimera contract (never assert a
  haplotype no individual is confirmed to carry). Verified end-to-end on a genome-wide 4+1+1 `--per-sample` run: 369 indel
  off-targets gained their co-SNP(s) with CFD byte-identical and the carrier roster narrowed to the
  exact cis set (e.g. chr10:105008040, carrier HG03711: `chr10_105008061_ATG_A` →
  `chr10_105008045_G_C,chr10_105008061_ATG_A`); SNP-only rows, PUTATIVE rows, and non-cis co-locations
  (indel and SNP carried by different individuals — e.g. the top-1000 rank-13 at chr18:12386282, whose
  worst-case alignment uses a SNP its indel carrier lacks) are left untouched. Because the corrected cis
  carrier count feeds the most-carriers representative tie-break, at a few dense multi-allele indel/STR
  loci a different — still real — observed haplotype becomes the cluster representative (no score
  change; this removes a pre-existing chimeric over-count where a worst-case indel alignment had used a
  cis SNP that most of its reported carriers do not carry).

- **`--per-sample` no longer over-states dense-window off-targets (the "chimera" over-call).** In a
  high-variant-density window (>`CRISPRME_IUPAC_CAP` co-located variants) the per-sample path emitted a
  synthetic **cross-sample** worst-case representative whose carrier roster was the UNION of different
  individuals' best-per-column alleles — a haplotype **no single individual carries** — which, being
  highest-CFD, became the per-cluster headline and over-stated both risk and carrier count (observed in
  the wild as a 15-SNP, "47-carrier", CFD-0.57 row at chr14:106770868 that nobody actually carries;
  the real observed haplotypes there are ≤2 variants, CFD ≈ 0). Per the per-sample contract — report an
  off-target only when ≥1 real individual carries it — that synthetic representative is no longer
  emitted; the dense region stays visible via the high-variant-density BED + the reference row. **Only**
  `--per-sample` is affected; the population-level default (which *should* report the conservative
  worst case) and sites-only panels are byte-identical. *Expected when diffing a `--per-sample` run
  against 2.7.9:* on dense-window variant sites the reported `highest_CFD` representative — and hence its
  `Start_coordinate` / `CFD_score` / `CFD_ref` and carrier roster — shifts to a **real** carrier
  haplotype (the chimera is gone), so a handful of such rows differ. This is the intended effect, not a
  scoring change: `CFD_alt` and all CFD arithmetic are byte-identical to 2.7.9 (validated on a
  genome-wide run — CFD identical on 189,119/189,132 shared alignments; the differing rows are exactly
  these double-bulge dense-window sites whose chimera representative was replaced by a real-carrier one).

- **Per-cluster representative tie-break → most carriers.** Among off-targets that tie EXACTLY on
  score and severity, the representative shown in `integrated_results` now prefers the haplotype carried
  by the most individuals (applies to the CFD, fewest-mm+b, and CRISPR-Bulge passes); a no-op for the
  unique-best ~99.7% of clusters.

- **CRISPR-Bulge aggregate specificity sanitized.** The CRISPR-Bulge `-1.0` out-of-domain sentinel
  (sites with ≥2 bulges, outside the model's domain) was summed into the aggregate `100/(100+Σ)`,
  driving the CRISPR-Bulge specificity negative; it is now excluded from the sum (no-op for CFD).

- **Indel-derived off-targets now get a `PAM_creation` value (#94, Ann Cirincione).** The field was
  only ever filled on the SNP path, so a PAM created by an insertion/deletion surfaced as `NA` even
  when it was the top-scoring site. The indel analyzer now mirrors the SNP predicate (if any reference
  base in the PAM window is incompatible with the motif, the variant created the PAM; store the ALT
  PAM). Writes the existing `PAM_gen` slot in place (no column shift); skipped when the window holds an
  IUPAC-undefined character (e.g. a bulge `-`). `calc_cfd` untouched.

- **`download --what index` no-ops on a *complete* install but still heals a partial one (#93, Ann
  Cirincione).** Re-running the documented download re-fetched the full tarball (~8.86 GB). It now
  skips the fetch when the install is complete, so a repeat costs no network call — but only when
  complete: a variant index needs its `_INDELS` companion, the shared reference index, and (unless
  `--no-genotypes`) its Tier-1 genotype store. Anything missing falls through to re-fetch just the
  absent piece, preserving the documented "no variant off-targets / Samples degraded → re-run
  download" recovery (and the `--no-genotypes`-then-add-genotypes workflow).

- **No more spurious `BiopythonWarning` at CLI startup.** `dash-bootstrap-components`' sdist drops a
  `pyproject.toml` into `site-packages/`, which Biopython misreads as "importing from a source tree"
  and warns about on every `crisprme.py` invocation. The Docker image already strips that stray file;
  conda/source installs did not, so `crisprme.py --version` printed the warning. Now silenced in two
  places: `crisprme.py` filters that one message before importing Bio (covers any install), and
  `install_from_source.sh` removes the stray file (matching the Dockerfile).

- **Web dropdown options no longer overlap for long labels.** The variant-index dropdown's long entries
  overflowed their fixed row height and collided with the next entry after repeated opens. The
  `variant-dataset` dropdown now reserves room for wrapped labels (`optionHeight`), and a CSS rule lets
  option text wrap and centers it within the row.

- **`targets-scoring` DNA-bulge CFD** now matches CRISPRme's reported CFD convention for in-place
  (Cas-OFFinder / report-style) alignments. The first cut front-trimmed the first `bulge_size`
  columns (a crispritz raw-format assumption), which misaligns an in-place gap that isn't at the
  5' end. It now DELETES the gap column(s) (the looped-out DNA base) from both strands and scores
  20-vs-20, so CFD's position weights stay guide-indexed. Verified against a real CRISPRme report:
  mismatch + RNA-bulge rows reproduce the reported CFD exactly; DNA-bulge CFD is the standard
  bulged-base-excluded value (CRISPRme's internal DNA-bulge gap placement isn't recoverable from an
  in-place alignment, so it stays a documented caveat). Updated `test_dna_bulge_delete_gap`.

- **Web "Load Example" no longer prunes bulged off-targets.** The Load Example button set the
  "Maximum edits" total-edit cap to 4 while loading 4 mm + 1 DNA + 1 RNA bulge, so `--max-total-edits=4`
  silently dropped every 4-mismatch-plus-bulge site (bulges had no headroom at the full mismatch
  budget). The cap is now 6 (= mm + bDNA + bRNA, the recommended default; matches the slider's own
  default and the CLI), and two inconsistent fallback defaults (5, 4) were aligned to 6. +1 regression
  test (`test/web/test_max_edits_default.py`, wired into web-e2e).

- **Web low-memory banner no longer nags a correctly-provisioned host.** A host set to the
  recommended ~64 GB reports a few GB less inside a container (GiB-vs-GB + VM reservation), so the
  banner fired even when Docker was set to 64 GB. Added a detection margin (`_DETECT_MARGIN_GB`), so
  the warning only fires for a clearly under-provisioned host; +1 unit test.

- **Removed the stale "download the offline version" note** from the web form footer. It was a
  hosted-web-service-era leftover — now that the service is retired, a user running the local web
  interface already *is* the full/offline version, so the note only caused confusion.


## [2.7.9] - 2026-10-11

### Added
- **Web low-memory pre-flight guard.** A new-user review found the clickable-Mac path dead-ends: the
  default first Submit is a genome-wide variant-aware search (peaks ~64 GB) that silently OOM-kills
  (exit 137) a typical 16-36 GB laptop. The web form now shows a **pre-Submit warning banner** when the
  host's detectable RAM is below the variant-aware whole-genome threshold (`memory_check.memory_status`,
  threshold `VARIANT_WGS_MIN_GB=64`, honoring `CRISPRME_MIN_GB`/`CRISPRME_SKIP_MEM_CHECK`), pointing to
  the escape hatches (reference-only index, single-chromosome/region, or raise Docker memory). Silent
  when RAM is undetectable; +7 unit tests. (A dedicated exit-137 message on the load page is a planned
  fast-follow; the search log already carries the `warn_low_memory` warning.)
- **`docs/RELEASE_CHECKLIST.md`** — a living release gate: the pre-tag consistency + functional +
  mechanics checks, plus a "regression traps" table that pairs every bug we've shipped (or nearly
  shipped) with the check that catches it. Appended to whenever a release pass finds a new problem.
- **`download --what assembly --sample <id>`** — fetch an HPRC individual (both haplotypes' genome +
  chain + chromAlias into `Assemblies/<id>/`) from the command line — the CLI equivalent of the web
  Data Manager's *Fetch from HPRC*. Closes the last CLI↔web gap; `docs/SETTINGS_DATA_MANAGER.md` now
  carries a full **CLI ↔ web parity** table (every web action maps to a `crisprme.py` command, so a
  CLI-built and a web-built install are interchangeable).
- **Build a variant-aware index from a single user-provided VCF.** `build-index-only --vcf` (and the
  web **Settings → Add a VCF dataset → Build index**) now **auto-splits** a single multi-chromosome
  VCF into the per-chromosome files the enricher (`crispritz add-variants`) requires — so a normal
  user VCF "just works". Previously a single multi-chromosome VCF failed with a cryptic
  empty-chromosome error (`Genomes/<g>/.fa`). (No-op for datasets already named per-chromosome, e.g.
  the shipped indexes.)
- **SNP-only custom VCFs now build a complete index.** A variant dataset with no indels previously
  aborted the whole build at the indel-indexing step (`Indels indexing failed`, empty fake-indel
  contigs); the indels index is now **skipped gracefully** and the SNP variant index is produced
  (SNP+indel co-occurrence simply does not apply to an indel-free dataset).
- **`docs/DATA_LAYOUT.md`** — a clear reference for *where your data lives*: the single user-owned
  data folder (CLI `--path`, Docker `-v …:/DATA`, Apptainer `--bind …:/DATA`), the **shared inputs**
  reused by every search (`Genomes/`, `genome_library/`, `VCFs/`, `Dictionaries/`, `samplesIDs/`,
  `Annotations/`, `PAMs/`, `Assemblies/`, `LiftoverFiles/`) vs. **where each run's outputs go**
  (`Results/<job>/`). Linked from the README quickstart, the Docker quickstart, and the Data-Manager
  guide.

### Changed
- **Final docs polish.** Reworded the web low-memory banner to actionable, web-achievable advice
  (switch to the reference-only index, or raise Docker memory — dropped the CLI-only "single
  chromosome" suggestion that has no web control). Restructured the densest METHODS sections into
  bullets/tables (§1 raw-vs-compressed encoding, §5 population-level comparison, §6 cancer
  annotations, §7 `Observed` column, §8 domain-of-validity + worst-case scoring). Reconciled the disk
  figures to one framing everywhere (~45 GB download → ~85 GB on disk → ~100 GB free; reference-only
  ~25 GB). Standardized the example guide on the full `CTAACAGTTGCTTTTATCAC` (what the web *Load
  Example* button loads). Added an SSH port-forward note to the website Singularity/Apptainer tab for
  reaching the web UI on a headless compute node.
- **CI + test hygiene.** Removed the redundant `docker-image.yml` workflow (single-arch, no-push,
  `checkout@v3` EOL — superseded by `docker-multiarch.yml`); wired the two Dash results-page unit
  tests (`test_results_page_assembly`, `test_results_summary_and_status`) into `web-e2e`; updated the
  `test_complete_test` mock index from `1000G2021` to the shipped default `HGDP1kGP`. The brute-force
  `validate-test` gate (exact set-equality vs an exhaustive aligner, fails the build on any mismatch)
  was audited and confirmed enforced on every PR + push + weekly cron.
- **complete-test docs (§2.2.2).** Download examples now use the default `NRG_3_hg38+hg38_HGDP1kGP`;
  added a note distinguishing the modern no-arg smoke (runs against whichever index you downloaded)
  from the legacy `--vcf_dataset` build-from-raw mode that feeds the `validate-test` gate (1000G only).
- **README/website readability pass (new-user experience).** Applied the docs-review nice-to-haves
  as surgical edits: Quickstart density trimmed (deduplicated the 3× RAM/default-index repeats, index
  pitch collapsed to the comparison table); Highlights five-index run-on → concise lead + table
  pointer; the mode/index decision bullet → a tight 2-column table; the report's `Observed`-column
  cascade and the §6 commercial-licence scenarios → tables; gnomAD-converter "pseudo-individuals"
  de-duplicated; `assembly-search` added to the function list; §1.1 reframed as the stable-2.1.x
  path with the CRISPRme+ update guidance first; Docker GPU wording qualified (CPU is the default,
  GPU optionally accelerates CRISPR-Bulge); `CRISPRITZ_REF` default (v2.9.0) + `http://localhost:8080`
  documented; a few internal-consistency nits (NGG→NRG example label, "three"→"four" commands,
  report list renumbered). No facts/commands/versions changed.

### Fixed
- **Final METHODS + linked-docs review pass (multi-agent, 2 blockers + consistency fixes):**
  - **`--max-total-edits -1`** — METHODS told users to pass `-1` to keep all off-targets, but the
    shipped validator rejects negatives (the search aborts). Dropped the `-1` option and fixed the
    stale in-code comment so doc + comment + validator agree (set it to `mm+bDNA+bRNA`).
  - **Mega-indel allele frequency** — METHODS said mega indels were "not yet frequency-annotated";
    they ARE (per-dataset `AF_<dataset>` + `AF_max` via the `indel_af` sidecar, read by the report).
    Scoped the "SNPs only" statement to the Tier-0 registry.
  - **`download --what vcf` deprecated in docs** — the raw per-dataset VCFs (`vcfs/1000G/`, `vcfs/HGDP/`)
    were never hosted on HF, so the command errored. Removed the misleading examples across README,
    DOCKER_QUICKSTART, the website, the web guide, DATA_LAYOUT and the legacy setup doc; clarified that
    the standard panels ship as prebuilt **indexes** (whose bundled genotype store covers sample-level
    analysis — raw VCFs are not needed) and that `--what vcf` is for custom/self-hosted repos. The CLI
    error now points at `--what index` / `build-inputs/` / `VCFs/<dataset>/`.
  - Consistency/accuracy: METHODS intro names the full default index + composition and the shipped
    version; stale present-tense "2.5.3 ships raw" / chrX "forthcoming" anchors reframed/resolved; a
    registry-cache speedup number reconciled to the source code (326→~89 s, 3.66×); the legacy-setup
    CRISPRitz-version self-contradiction collapsed to `v2.9.0`; three stale `INPUT_FORMATS` code line
    citations corrected; `PRECOMPUTED_INDEXES` `<vcf>` token defined once; `CRISPRME_IUPAC_CAP` default
    (10) stated. Also fixed a README download example that repeated `HGDP1kGP` under a "1000G-only" label.
- **Docs/website review pass (multi-agent section-by-section review, 3 blockers closed):**
  - README default-index per-source split said `HGDP (929) + 1000 Genomes (3162)` — contradicted the
    correct `925 HGDP + 3,166 1000 Genomes = 4,091` elsewhere in the same file; corrected.
  - README report-Summary example cited `3,477 individuals for 1000G+HGDP` (the retired-hybrid
    panel size) — corrected to `4,091` for `NRG_3_hg38+hg38_HGDP1kGP`.
  - Website Conda/Mamba path told pre-2.2.0 stock-Bioconda users to run CRISPRme+-only subcommands
    (`download`, `web-interface`) the stock package lacks — now routes them to the legacy `setup`
    and points the fast download/web path at the CRISPRme+ Docker/source install.
  - Also: §2.1 directory list was missing `genome_library` / `Dictionaries` / `Assemblies` and the
    one-working-folder framing (added, links DATA_LAYOUT.md); `sampleIDs` → `samplesIDs`; website
    hero now names the shipped default; the first-launch `~45 GB` figure is reconciled with the
    `~100 GB` disk budget; version-specific "unreleased 2.7.8 / in preparation" prose made
    version-agnostic; "dict-less" download jargon relabelled "variant-aware".
- **Web default variant panel now honours the re-phased default.** The first-paint variant-dataset
  seed (`_preferred_variant`) was called without the genome, so its genotyped-**phased** preference +
  sample-count tiebreak were skipped — with two or more phased panels installed it could seed a
  non-default panel. Both call sites now pass the genome, so the default page consistently proposes
  the SHAPEIT5 re-phased `NRG_3_hg38+hg38_HGDP1kGP` as the recommended variant default (still fully
  name-free: the preference is read from each index's `data_type` manifest).
- **No surface advertises a retired index.** Swept the repo for indexes removed from HuggingFace in
  v2.7.8 and repointed the few remaining references at the live catalog: the `seq_script/merge_panels/`
  README (retirement banner + repoint to `HGDP1kGP`), the `onekg2021_build.sh` header (reworded from
  "the default index" to an *alternative* single-source panel), a `pages_utils` docstring example, and
  the `test_variant_dropdown` docstring (clarifies the two-source `1000G_HGDP` fixture is a generic
  discovery-mechanism stand-in, not a shipped index). Archived 2.5.x docs already carry retirement
  banners and are unlinked from the live site.
- **Web Settings / Data-Manager (3 fresh-user-QA findings):**
  - Clicking **Download index** on an already-installed index re-downloaded all its GB and blocked
    the single-slot job queue; it now short-circuits with *"Index '…' is already installed — nothing
    to download."*
  - An in-progress download's `.stage_<genome>` staging dir was listed as an **installed genome** (and
    offered in the build / search / delete dropdowns); dot-dirs (staging + assembly markers) are now
    excluded from the installed-genome listing.
  - Launching a **Download / Build** left the inline feedback blank (job status appeared only in the
    Activity panel); every job launch now confirms *"Started — progress is shown in the Activity
    panel."*
- **Web-submitted `assembly-search` was broken** (`ModuleNotFoundError: No module named
  'validate_inputs'` at startup, so every personal-assembly search failed immediately). Two causes,
  both fixed: (1) crisprme.py's install-path correction hard-stripped a 3-char `"bin"` suffix, so when
  the web spawned its `opt/` copy (`<prefix>/opt/crisprme/crisprme.py`) the PostProcess path was
  garbled — it now detects whether PostProcess is a direct sibling (opt/ copy or source checkout) or
  reached via `opt/` (bin copy); (2) the job was launched with `--debug`, whose cwd-relative
  PostProcess path also missed (the job's cwd is the Results output dir) — `--debug` now only switches
  to that path when it exists, and the web assembly-search command no longer passes `--debug`. Verified
  live on a real HPRC individual (both haplotypes search + `--gene_annotation` closest-gene wired).
- **`assembly-search` (and `complete-search`) aborted on a read-only install** (Apptainer/Singularity
  SIF) with `truncate: … vuoto.txt: Read-only file system` (issue #143): the results-integration step
  stripped a trailing newline from the "no empirical data" placeholder, which is the shipped read-only
  `PostProcess/vuoto.txt`. It now only truncates a real, non-empty, writable file. (Docker's writable
  layer masked this; it affected HPC SIF runs.)
- **Web-submitted `assembly-search` crashed in the combine step** with `shutil.SameFileError`: the web
  writes the guide straight into `Results/<job>_combined/.guides.txt`, then `write_combined_guides_file`
  copied the guide onto that same path. It is now a no-op when source and destination are the same
  file. With this, a web assembly-search runs **end-to-end** — both haplotype searches → hg38
  reconciliation → the combined report with per-haplotype CFD + closest-gene columns (verified live on
  HPRC HG01255: 60,745 off-targets, `Annotation_closest_gene_*` present).
- **Docs/consistency precision pass (full release-gate QA, 161 PASS / 14 WARN / 2 FAIL).** Corrected
  the HGDP+1kGP cohort split to **925 HGDP + 3,166 1000 Genomes** (the shipped BCF header-derived
  split; still 4,091) across README / METHODS / DOCKER_QUICKSTART / PRECOMPUTED_INDEXES;
  `gnomAD-converter --help` now lists **v3.1, v4.0, and v4.1 (incl. joint VCFs)** (it already
  supported v4.1); documented the `scorer-env list` action; fixed the Dockerfile CRISPRitz comment
  (→ v2.9.0), the README `download --what` narrative (added `assembly`), the DESIGN_mega "two
  predefined indexes" framing (now the published set), the TensorFlow "2.13.1" → "2.13.x per-arch"
  wording (Apple-Silicon is 2.13.0), and drift-prone code-line citations in INPUT_FORMATS.
- **`complete-test` / `setup` help clarified.** `complete-test`'s "no prebuilt index" pointer now
  recommends the default `NRG_3_hg38+hg38_HGDP1kGP` (was `1000G2021`); the `--vcf_dataset`
  (1000G | HGDP | 1000G+HGDP) help now states it is the **legacy build-from-raw-VCF test panel**
  (1000 Genomes 2019 + HGDP 2019), explicitly **not** the shipped re-phased HGDP1kGP default and
  **not** the retired `NRG_3_hg38+hg38_1000G2021_HGDP` index. `setup --help` now points to the
  faster modern `download --what all` path (setup is the legacy from-scratch FTP download).

## [2.7.8] - 2026-10-11

### Added
- **Build a new-PAM index for a published dataset.** The PAM-independent enriched-genome "kits"
  for the re-phased `NRG_3_hg38+hg38_HGDP1kGP` and sites-only `NRG_3_hg38+hg38_mega` datasets are
  now published on HuggingFace (`build-inputs/hg38_HGDP1kGP.tar.gz`, `build-inputs/hg38_mega.tar.gz`,
  each with `SHA256SUMS`). A new `docs/PRECOMPUTED_INDEXES.md` section documents downloading a kit
  and indexing **any** nuclease/PAM (NGG, NNGRRT, TTTN, …) off it, plus the full source-to-kit
  reproduction (SHAPEIT5 re-phasing + data provenance for HGDP+1kGP; the five-source merge + MAF
  filters for mega). `METHODS.md` gains a re-phasing-provenance note and a "Re-indexing for a
  different PAM" subsection; `README.md` and the data-setup guide cross-link the new section.

### Fixed
- **Web `assembly-search` could not be submitted.** The hidden genome/chain/chromAlias holders
  auto-populated from the Individual picker were `dcc.Dropdown`s, which silently drop a value that
  is not one of their options; the picker stores cwd-relative *paths*, so the submit handler read
  `None` for all six and always reported "missing fields". They are now `dcc.Store` components
  (hold arbitrary strings), so a web-submitted assembly-search launches correctly.
- **Web `assembly-search` did not add closest-gene columns.** The submit handler discarded the
  gencode companion returned by `build_active_annotation`; it now passes `--gene_annotation`
  (the bundled `gencode.protein_coding.bed`) when the built-in annotation bundle is enabled, so a
  web-submitted job gets the `Annotation_closest_gene_*` columns (parity with complete-search).
- **`synth_sites_gt.py` wrote plain gzip instead of BGZF**, so `tabix` rejected the output
  ("not BGZF") and the documented sites-only (mega) fake-indel reproduction recipe failed. It now
  writes BGZF via pysam (with a BGZF-magic-byte test assertion).
- **`seq_script/merge_panels/mega_build_indels.sh` never placed the fake-indel genome** into the
  index directory — its `strip_and_place` glob (`*fakechr<N>*`) did not match the enricher's
  `fake_hg38_mega_<chr>` directory name. Search was unaffected (it uses the index bins), but the
  PAM-independent rebuild genome was missing; the glob is fixed.
- **Migration-guide license row was misleading.** The "What changed vs stock CRISPRme" table framed
  CRISPRme 2.1.14 as merely "AGPL-3.0 (copyleft)", implying AGPL alone permitted commercial use.
  Stock CRISPRme was dual-licensed — AGPL for academic/non-commercial use **plus** a separate
  commercial license for commercial use. The row now states that both CRISPRme and CRISPRme+ require
  a separate commercial license for commercial use; only the academic-use license changed
  (AGPL-3.0 → MGBOA 1.0).

### Changed
- **Documentation + website consistency sweep.** Version strings were unified to `2.7.8` across
  `README.md`, the website (`docs/index.html`), `docs/DOCKER_QUICKSTART.md`, the web-interface
  user guide, `docs/PRECOMPUTED_INDEXES.md`, and `install_from_source.sh` (several surfaces had
  drifted to `2.6.x` / `2.7.1` / `2.5.5`). The legacy data-setup guide now points to the modern
  `crisprme.py download` flow for first-time setup.
- **Fresh-user consistency pass.** A full install→run→build→interpret walkthrough surfaced and
  fixed: the **default index** is now consistently **HGDP + 1000 Genomes** across README, the
  website, the Docker quickstart, and the web-interface guide (several spots still named the old
  1000G-2021 default); the website RAM figure is corrected to **64 GB** for a genome-wide
  variant-aware search and now surfaces the **one-line installer**; a CLI-guide example that
  dropped the PAM `N`-padding (and would prune all off-targets) is fixed; the web guide gains an
  **"Understanding your results"** reference (CONFIRMED vs PUTATIVE, carriers, key columns); and
  a `docs/FUTURE_WORK.md` captures deferred enhancements.
- **One-click installer moved into this repo.** `install.sh` / `install.ps1` (previously served
  from the legacy `pinellolab/CRISPRme` repo) now live in `docs/` here and are served at
  `https://pinellolab.github.io/crisprme-plus/install.{sh,ps1}`, so the release owns its own
  installer. Their default first-run download is switched from the old `1000G-2021` index to the
  recommended **HGDP + 1000 Genomes** default, and the stated size corrected to **~45 GB** (the old
  "~75 GB" was wrong). The one-liner URLs in the README, website, and `install/README.md` are
  updated to match. (The old URL should be left as a redirect for existing bookmarks.)

### Removed
- **Retired ~71 GB of stale/superseded precomputed indexes from HuggingFace.** Deleted the
  superseded hybrid `NRG_3_hg38+hg38_1000G2021_HGDP` (+ its genotype companion), which is replaced
  by the fully re-phased `HGDP1kGP` default; the orphaned v2.3.0 dict-less lineage
  (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, `NGG_3_hg38-dictless+hg38_1000G_HGDP`,
  `genotypes_hg38_1000G_HGDP`); and the unused reference-only `NGG_3_hg38`. The hosted catalog is
  now the four production indexes (`HGDP1kGP` default, `1000G2021`, `HPRC`, `mega`) plus
  reference-only `NRG_3_hg38`; all user-facing docs were scrubbed to match (catalog entries,
  download commands, and every worked example repointed to a hosted index).

## [2.7.7] - 2026-10-10

### Added
- **Closest-gene annotation for `assembly-search`** (`--gene_annotation`, #66). A
  nearest-feature (BEDOPS `closest-features`) lookup adds the closest gene + signed
  distance to every hg38-mapped off-target of a personal diploid assembly — `Gene` +
  `Gene_distance_kb` in the curated report tables / web result page, and
  `Annotation_closest_gene_name` / `_ID` / `_distance_(kb)` in the full TSV
  (distance `0.0` = overlapping, negative = the other side; haplotype-non-mappable
  sites have no hg38 coordinate and so get no closest gene). Validated end-to-end on a
  real diploid assembly (HG01255 chr19 → 952/960 sites annotated) and in the web UI.
- **Variant-dataset retrieval recipes** under `paper/reproducibility/variant-datasets/`
  for 1000G 2019/2021, HGDP, gnomAD v4.1, TOPMed and All-of-Us (paper reproducibility,
  docs-only; #45).

### Changed
- Docs: README (`--gene_annotation` nearest-feature columns + assembly-search support +
  variant-datasets pointer) and the web-interface user guide (assembly result columns).

## [2.7.6] - 2026-10-09

Stability release candidate. One correctness hardening on top of v2.7.5; no change
to the search engine, default index, or Docker/Bioconda toolchain.

### Changed
- **`--samplesID` over-listing is now a hard ERROR (was a warning).** A `--samplesID`
  that lists samples genotyped in *none* of the provided VCFs silently inflates the
  panel AN (phantom samples are counted as hom-ref) and deflates every reported
  allele frequency (Samples / AF / MAF) — no crash, just wrong numbers (#46). Since
  there is no legitimate reason to list a sample absent from every VCF, pre-flight
  validation now refuses to run and points to the fix (a VCF-filtered samplesID, or
  the `registry_fix_an` / genotyped-panel transform). The shipped default indices
  (built with correct panels) and index-path searches are unaffected — this only
  guards fresh `--genome --vcf --samplesID` builds with a mis-listed panel.

### Fixed
- **Default 1000G samplesID over-listed its VCF panel.** The shipped
  `samplesIDs.1000G.txt` listed 3500 samples while the complete-test 1000G VCF
  (`ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019`) genotypes only 2548 —
  the 952 phantom samples were counted as hom-ref, deflating every reported AF
  ~1.37× (#46), and (with the hard ERROR above) blocked the default 1000G
  `complete-test`. The panel is now filtered to exactly the 2548 VCF-genotyped
  samples (HuggingFace data repo + in-repo fallback + `MD5SAMPLES` pin updated).
  HGDP was verified already correct (929 ↔ 929).

## [2.7.5] - 2026-10-09

Stability release: off-target scoring correctness (CFD audit), assembly-search
coordinate + report fixes, and test-suite hardening. No change to the search
engine, the default index, or the Docker/Bioconda toolchain.

### Fixed
- **CFD scorer — INDEL-path N-in-PAM divergence.** `analisi_indels_NNN.calc_cfd`
  had drifted from the canonical SNP-path `new_simple_analysis.calc_cfd`: an `N`
  in the scored 2 bp PAM kept a non-zero CFD on the INDEL path, whereas the SNP
  path zeroes a non-canonical PAM (`pam_scores.get(pam, 0.0)`, issue-#94 guard);
  it also carried a dead `if "N" == sl` branch. The INDEL twin is now byte-for-byte
  equivalent to the SNP path (an `N` genomic base is unknown → zeroing is
  consistent and conservative), pinned by a new regression test
  (`test_calc_cfd_twins_agree.py`). A full CFD correctness audit confirming the
  core is bit-identical to the Doench-2016/CRISPOR reference is documented in
  `docs/analyses/cfd_correctness_audit.md`. (#63)
- **assembly-search — hg38 off-target coordinates were off by one and lost
  strand.** Off-targets are now lifted via both end bases as stranded records,
  with orientation resolved and non-mappable sites refused rather than
  mis-placed; two edge-case crashes (contig-edge `impg`, empty lift) fixed. (#62)
- **assembly-search report — stale curated-output flags + inconsistent zip
  schema.** The module-level curation flags (`_DROP_MAF`, `_PRESENT_ANN_KINDS`,
  `_HAS_VARIANTS`, `_DROP_KINDS`) are now set on entry to *both* report builders,
  so an assembly report rendered after a complete-search can no longer inherit
  stale flags; and every table in an export zip now shares one schema (decided
  once from the pooled frame). Complete-search report output is unchanged
  (verified byte-identical). (#64)
- **Test suite — collection-order-dependent failures.** `test_bulge_cap.py` no
  longer leaks its import-time `app`/`dash`/library stubs into `sys.modules`, so a
  later-collected test's `from app import …` resolves the real module regardless
  of collection order (fixes two order-dependent `TestAnnotationLogNote` failures).

### Changed
- **assembly-search report columns.** The combined annotation column is split
  into per-annotation columns matching the complete-search report, and the new
  haplotype-context columns (`Haplotype_origin`, `hg38_orientation`,
  `hg38_lift_confidence`, `n_copies`/`copy_loci`, `hg38_end`) are surfaced;
  assembly-irrelevant (e.g. VCF-specific) columns are dropped. (#64)

## [2.7.4] - 2026-10-08

### Changed
- **Pinned CRISPRitz to v2.9.0** (`Dockerfile` + `install_from_source.sh`). CRISPRitz
  v2.9.0 is now a **pure enumeration engine** — all scoring, annotation and reporting
  were removed from CRISPRitz (they live in CRISPRme+), which drops its azimuth
  `scikit-learn==1.1.3 / numpy==1.24.4 / pandas==2.0.3 / scipy==1.10.1` pins. This has
  **no functional effect on CRISPRme+** (which only uses CRISPRitz's `index-genome` /
  `search` / `add-variants` engine — verified byte-identical, golden-search regression
  on x86_64 + arm64) and unblocks publishing CRISPRme+ on **Bioconda** (the `crisprme`
  package can now depend on `crispritz` without the old, unsolvable ML version pins).

## [2.7.3] - 2026-10-08

### Changed
- **Pinned CRISPRitz to v2.8.4** (`Dockerfile` + `install_from_source.sh`), which carries the
  CRISPRitz #44 enricher fix: on **multiallelic** records the `add-variants` enricher previously
  matched a sample to an ALT allele by substring-searching the raw genotype text and stopping at
  the first hit, so (1) a record with ≥10 ALTs could bind a sample to an allele it does not carry
  (`0|12` matched allele `1`) and (2) a sample carrying two different ALTs (`2|1`) had its second
  allele dropped. The genotype is now parsed into integer alleles and the sample is recorded under
  every ALT it carries. CRISPRme+'s shipped indices are `bcftools norm -m -any` pre-split to
  biallelic (unaffected); this corrects per-sample carrier assignment for **user-supplied
  multiallelic VCFs**. The biallelic path and the enriched FASTA are unchanged.
- **The CRISPR-Bulge scorer is now a guaranteed part of a source install, not best-effort.**
  `install_from_source.sh` builds the `cbulge` env **and compute-self-tests it** (`scorer-env
  create && scorer-env check`), and **aborts the install** if either fails — mirroring the Docker
  image (`build_scorer_envs=1`, failure fatal). Opt out deliberately with
  `CRISPRME_SKIP_SCORER_ENV=1`.

### Fixed
- **A missing/broken CRISPR-Bulge scorer env no longer silently degrades a search to CFD-only.**
  `complete-search` now **fails fast up front** (before the long search) with the one command that
  fixes it when the scorer env is not ready, instead of discovering it per-chromosome in
  post-analysis. A deliberate CFD-only run is still possible with `CRISPRME_SKIP_SCORER_ENV=1`,
  which proceeds with a single loud notice. The scorer-unavailable banner is written to **stdout,
  not stderr** (the post-analysis pipeline treats any stderr byte as fatal, so the old banner
  turned a deliberate CFD-only run into a failed run — supersedes #59). +1 regression test.

## [2.7.2] - 2026-10-08

### Fixed
- **#46 (co-occurrence path): Tier-1 genotype axis was not filtered to VCF-genotyped samples.**
  The phantom-sample filter (`genotyped_samples`) was wired on the Tier-0 registry
  (`compile_from_dict`) but not on the Tier-1 genotype store (`compile_genotypes_from_dict`), so an
  over-listing samplesID inflated the Tier-1 sample axis / AN and **deflated the per-sample
  CONFIRMED-cis joint-AF denominator**. The same `genotyped_samples` set is now threaded through to
  `tier0_compile.build_sample_meta` (already filter-capable); no-op when `None` (batteries install /
  source VCF absent). Build-path only — correctly-built shipped indices are unaffected (e.g.
  HGDP1kGP's header-derived 4,091-sample panel has no phantoms). Added
  `test_genotyped_samples_filter_drops_phantom_from_axis`.

## [2.7.1] - 2026-10-07

### Fixed
- **Report: `IntOGen_cancer_driver` column rendered `-` in every row.** `_curated_cell()` had no
  `intogen` branch, so the value fell through to the trailing `else: v = None` and displayed as `-`
  for every off-target. Added the missing `elif kind == "intogen"` branch (regenerating the curated
  top-1000 goes from 0/1000 to 45/1000 populated IntOGen cells; spot-checked on BRCA2 / EPHA3 /
  PTPRK coding-region off-targets). Added a guard that every `kind` declared in `CURATED_COLUMNS` is
  matched by a branch in `_curated_cell`. (#58, @anncir1)
- **CLI: `generate-report --no-maf` was silently ignored.** The documented top-level
  `crisprme.py generate-report` wrapper hand-parsed only its four value-flags, so the boolean
  `--no-maf` was dropped with no error and the MAF column was still emitted. It is now parsed
  (`"--no-maf" in args`) and threaded into `build_report(drop_maf=...)`, listed in the help, and
  covered by a test asserting the MAF column is actually removed.
- **Web: per-cluster view/download silently dropped all 24 CRISPR-Bulge columns.**
  `change_headers_bestMerge.py` reordered the 3-scorer `.bestMerge.txt` through a fixed `new_order`
  list that predates CRISPR-Bulge and contains no `CRISPR_BULGE_*` names, so the pandas column
  selection silently discarded the entire CRISPR-Bulge block (no crash — all listed names were
  present) from the genomic-region cluster download. The block is now preserved (no-op on a
  CFD-only input); new `test_change_headers_bestmerge.py` locks it in.
- **Web: the Reset button never cleared the threshold dropdown.** The reset callback's Output
  targeted the prop `"value "` (trailing space) instead of `"value"`, so Dash silently never wrote
  it. Fixed, plus `test_dash_output_props.py` guards every `pages/` Output/Input/State prop against
  stray-whitespace / empty names (the same silent-inert class as the IntOGen column).

### Added
- Coverage tests for the "declared-but-dead" bug class (an output field/flag that is declared but
  has no populating/wiring code, silently rendering a placeholder/no-op in every case): the
  CURATED_COLUMNS↔`_curated_cell` guard (#58), the `--no-maf` drop assertion, the CRISPR-Bulge
  per-cluster preservation test, and the Dash callback-prop whitespace guard.

## [2.7.0] - 2026-10-07

### Added
- **Scorer-env build is registry-driven + compute-verified (robust release).** The Docker image
  now provisions **and COMPUTE-self-tests every env in the `scorer_env.SCORER_ENVS` registry**
  (via `scorer_env.build_all`) instead of hardcoding `cbulge` — the modular scorer registry is the
  single source of truth, so adding a scorer env automatically builds it into the image. The build
  runs a real compute self-test (`scorer_env.selftest`: scores a perfect-match vs a mismatched
  target and asserts real, in-`[0,1]`, non-degenerate scores) rather than an import-only check, so
  a broken / miscomputing / weight-less scorer **fails the image build** instead of silently
  shipping (catches the tensorflow-metal GRU-miscompute class). `crisprme.py scorer-env
  check`/`doctor` now run the same compute self-test.
- **New recommended default variant index: `NRG_3_hg38+hg38_HGDP1kGP`** — the gnomAD HGDP + 1000
  Genomes callset (4,091 samples = 929 HGDP + 3,162 1000 Genomes), jointly **SHAPEIT5-re-phased** so
  **both** cohorts are phased: HGDP off-targets now get **CONFIRMED cis** + named per-sample carriers
  (not PUTATIVE-from-phasing), covering HGDP's seven genetic regions. Published to HuggingFace
  (`indexes/NRG_3_hg38+hg38_HGDP1kGP.tar.gz` + `genotypes_hg38_HGDP1kGP.tar.gz`); SNP+indel
  searchable; supersedes the hybrid `NRG_3_hg38+hg38_1000G2021_HGDP`.
- **PAR-aware chrX ploidy.** chrX allele frequencies now use the correct ploidy per region — diploid
  in the pseudoautosomal regions (PAR1/PAR2, AN=8182 for a 4,091-sample panel) and haploid-male in
  the non-PAR region (AN=5990) — in both the registry build and the search-time k≥2 co-occurrence
  denominators. Fixes a ~1.37× PAR allele-frequency inflation; autosomes/chrY unchanged.

### Changed
- **A missing/broken ML scorer env is now LOUD.** A search that degrades to CFD-only because the
  CRISPR-Bulge env is absent/broken prints an un-missable banner to **stdout + stderr** (was a
  one-line stderr note buried in `log_error`) — so a CFD-only run is never silent.
- **Default variant index flipped to `NRG_3_hg38+hg38_HGDP1kGP`** (was `NRG_3_hg38+hg38_1000G2021`).
  The web UI picks it automatically (richest genotyped-phased panel, by sample count — no hardcoded
  name); CLI help + docs recommend it first. `1000G2021` remains the 1000-Genomes-only alternative;
  `1000G2021_HGDP`, `HPRC`, and `mega` remain available.

## [2.6.2] - 2026-09-28

### Added
- **`assembly-search`: direct haplotype-vs-haplotype reconciliation (impg / minimap2).** The two
  personal haplotypes are now reconciled directly against each other (cached, threaded pairwise
  alignments), independent of `liftOver` — more accurate two-haplotype off-target reconciliation.
  `minimap2` / `impg` are checked at run start (`check_impg_available()`). (#53)
- **`assembly-search`: hg38 annotation of reconciled sites (`--annotation`).** Reconciled
  off-target sites can be annotated against the hg38 functional / gene tracks; a failed or missing
  annotation degrades gracefully and never discards an otherwise-finished two-haplotype search. (#55)

### Fixed
- **Non-human genome robustness.** Genome-layout detection and empty-result handling fixes
  surfaced while running CRISPRme on a non-human (pig) genome. (#56)
- **Restored the brute-force benchmark CI gate (`validate-benchmarks`).** The v2.6.0
  `complete-test` rewrite made `run_crisprme_test` always require a prebuilt index and run a
  single example-guide smoke, which silently dropped the CI's download + per-benchmark
  on-demand-build flow — so `validate-test` had no predictions to compare against the committed
  brute-force ground truth and the gate failed. `run_crisprme_test` now branches: the default
  (user) path keeps the no-auto-build prebuilt-index smoke, while the CI path (when
  `CRISPRME_ALLOW_ONDEMAND_BUILD` is set, as the workflow does) downloads the test data and runs
  every registered benchmark with on-demand index build, producing the per-benchmark output dirs
  `validate-test` needs. The brute-force comparison against CRISPRme's predicted off-targets is
  exercised again.
- **Guardrails so this class of regression can't recur silently.** Added `test_complete_test.py`
  — a fast, hermetic dispatch test (in the push-to-main unit suite) that pins the two-mode
  contract: CI mode (`CRISPRME_ALLOW_ONDEMAND_BUILD` set) downloads + builds each registered
  benchmark on demand with no `--index-path`, and user mode requires a prebuilt index and never
  downloads/builds. It fails against the broken rewrite. Also made `validate-benchmarks` (the
  brute-force gate) run on **direct pushes to `main`** touching the pipeline — not just PRs +
  the weekly cron — so a fast-forward release can no longer bypass it (the gap that let the
  breakage sit unnoticed since mid-September).

## [2.6.1] - 2026-09-28

### Added
- **The default image is now GPU-capable — one image, both modes.** `pinellolab/crisprme:latest`
  (and `:<tag>`) builds its **amd64** layer with the conda-forge CUDA TensorFlow
  (`tensorflow=2.13=cuda*`, validated correct + fast on an NVIDIA A100), so the CRISPR-Bulge scorer
  runs on the GPU under `docker run --gpus all …` and **falls back to CPU transparently** when no
  GPU is visible (the compute-backend device guard; byte-identical results) — no flag needed for
  CPU. The **arm64** layer stays CPU (CUDA has no ARM build), so Apple Silicon is unaffected. A
  single `docker pull pinellolab/crisprme:latest` therefore works everywhere and uses the GPU when
  present; the only cost is a larger amd64 image (bundles cudatoolkit/cudnn). Build the Dockerfile
  with `--build-arg scorer_backend=cpu` for a lean CPU-only image if size matters. Selection is via
  the Dockerfile `scorer_backend=auto` default keyed on `TARGETARCH`.
- **`Gene_region` column — exon/intron/UTR granularity for the gene annotation, on every
  report surface.** States in one word where an off-target sits relative to the nearest
  protein-coding gene: `CDS` (coding exon), `5'UTR`, `3'UTR`, `exon` (non-coding-transcript
  exon), `intron` (inside a gene body but not an exon), or `intergenic`. Computed from the
  existing `Annotation_GENCODE` feature-set (no new pipeline pass) and kept coherent with the
  `Gene` / `Gene_distance_kb` columns — `intergenic` is pinned to `Gene_distance_kb > 0`, so
  the region never contradicts the distance. Surfaced **consistently everywhere** from ONE
  shared helper (`utils.gene_region_class`): the raw `integrated_results.tsv`
  (`Annotation_gene_region`, written by the integrator), the in-report curated table, and
  every download TSV (all tiers, panel, top1000, variant_created); documented in the report's
  annotation legend. The raw `GENCODE` column still keeps the full overlapping feature-set.
  Report curation also derives it on the fly, so it renders on pre-existing
  `integrated_results.tsv` files too.

### Fixed
- **GPU image now actually builds its CUDA scorer env (and a broken scorer env can no longer
  ship silently).** The `cbulge` CUDA solve failed on GPU-less builders (e.g. GitHub CI) with
  `__cuda missing`, because the micromamba base activation exports `CONDA_OVERRIDE_CUDA=""` and
  `scorer_env` used `setdefault` (which will not replace an existing empty value). Now the
  override is force-set to `11.8` when empty/absent, so the CUDA TensorFlow build resolves on any
  builder (the image still falls back to CPU at runtime when no GPU is present). Additionally, the
  Dockerfile scorer-env build step is now **FATAL** — a failed or incomplete `cbulge` env fails
  the image build (verified by importing TensorFlow) instead of being swallowed by a warning and
  shipping an image with no scorer.

## [2.6.0] - 2026-09-28

### Changed
- **New default search: 4 mismatches + 1 DNA + 1 RNA bulge (6-edit budget).** The web form
  and CLI now default to `--mm 4 --bDNA 1 --bRNA 1` with `--max-total-edits 6` (was a 6/2/2
  ceiling with a 4-edit budget). This keeps the default in the CRISPR-Bulge scorer's validated
  **single-bulge** domain and, crucially, sizes the total-edit budget to the params so a bulged
  off-target at the full mismatch budget is no longer silently pruned. Simple mode (Advanced
  panel closed) fixes the per-type bulge cap at 1; the CLI likewise defaults omitted
  `--bDNA/--bRNA` to a single bulge of each type. Open the Advanced panel / pass explicit
  `--bDNA/--bRNA` (and raise `--max-total-edits` to their sum) for deeper searches.
- **Per-sample availability + label are index-aware on the web.** The "Per-sample" analysis
  option is auto-enabled only for genotyped panels (grayed out + reverted for sites-only like
  `mega`), and its label now states what it buys on the chosen index: **CONFIRMED cis + named
  carriers on a phased panel**, **carriers with PUTATIVE cis on an unphased panel**, or the
  hybrid mix.
- **CRISPR-Bulge is the machine-learning off-target scorer.**
  The legacy vendored RandomForest (scikit-learn) ML scorer is retired in favor of CRISPR-Bulge
  (Yaish & Orenstein, *NAR* 2024; TensorFlow GRU ensemble; MIT), which is more accurate
  — ~2.1× overall and ~5.4× on bulge-containing off-targets by AUPR on Refined-TrueOT —
  and faster per off-target on CPU. **CFD remains the primary score** (fast, in-process,
  pure-Python); CRISPR-Bulge is the secondary ML score. To keep the main environment free
  of a heavy TensorFlow pin, the scorer runs in its **own dedicated conda env** (`cbulge`),
  invoked over a persistent-worker batch protocol; if that env is absent the run completes
  on CFD with a single warning (the ML column is `-1`). Reports and the web UI now label
  the column/plots/legend "CRISPR-Bulge" and calibrate on its own [0,1] scale
  (0.5 / 0.2 / 0.1 bands); the physical `CRISPR_BULGE_score_*` column name and `.bestCRISPR_BULGE`
  file are retained as stable internal identifiers so older reports still render.
  The ensemble is a single-bulge model, so off-targets whose alignment needs ≥2 bulges (or a
  ≥2-bp bulge) are NOT ML-scored — they are nulled to `-1` and carried by CFD (the primary
  score) + mismatch/bulge counts. (An earlier build reduced such alignments to a single bulge;
  an adversarial review showed that reduction could corrupt the alignment and, on a non-monotonic
  model, was not provably conservative, so it was removed in favor of the honest null. Each
  off-target locus is still ML-scored on its best ≤1-bulge alignment when one exists.)
- **Modernized the main environment now that the legacy ML scorer's hard pins are gone.**
  `scikit-learn` is removed entirely (nothing in the main env imported it after the legacy
  RandomForest scorer was retired);
  `numpy` 1.24→≥1.26, `scipy` 1.10→≥1.11, `pandas` 2.0→≥2.1, and the `matplotlib-base` ceiling
  is dropped. `numpy` is capped `<2` and `pandas` `<3` as deliberate guards (the pysam /
  CRISPRitz C-ABI, and pandas 3.0 copy-on-write) — both are validated fast-follows. The
  full unit suite passes on the rebuilt, un-pinned environment (783 tests on the
  fresh-from-`environment.yml` source install).
- **The CRISPR-Bulge scorer env is architecture-aware, with identical scores on x86-64 and
  arm64.** conda-forge ships `tensorflow-cpu=2.13` only for x86-64 (aarch64 has just 2.18/2.19),
  so the `cbulge` env is now built per-arch behind the same name: conda `tensorflow-cpu=2.13` on
  x86-64, the conda CUDA build on GPU hosts, and — new — the **same TensorFlow 2.13.1** from the
  official PyPI aarch64 wheel on Linux ARM (base deps from conda, TF from pip, mirroring the
  Apple-Silicon Metal path). Same version + weights → validated **identical**: on a 254-pair panel
  (150 with bulges) x86-64 vs aarch64 scores agree to max |Δ| = 3.6e-7 with **zero** differences at
  display precision and **zero** threshold-band flips. The multi-arch Docker image therefore ships
  the full scorer on both arches (previously arm64 would have degraded to CFD-only). The scorer-env
  solve also honors the same `CONDA_CHANNEL_BASE` mirror as the main env (e.g. `https://prefix.dev`),
  so `scorer-env create` works on networks that block `conda.anaconda.org`.
- **Alternative-alignments output is now mode-driven (off by default in the
  population-level analysis).** The `..._all_results_with_alternative_alignments.tsv`
  dump — the *non-best* alignments per locus — grows combinatorially with the edit
  budget and its annotation/sort/integration was the single most expensive
  post-analysis step on a large search, yet the report, summary tables, mismatch/bulge
  matrix and top-100 panel are all built from the *best* alignment per locus and never
  read it. It is now emitted only under `--per-sample`; the default population-level
  run drops the non-best rows right after the merge, skipping the expensive step. New
  `--alt-alignments` / `--no-alt-alignments` flags override either way.
  `assembly-search` no longer forces it on either (its haplotype reconciliation reads
  only `*_integrated_results.tsv`), so its per-haplotype searches run lean too. (A worst-case
  6mm/2+2/max-8 search on the 5-source mega had ballooned this file past 160 GB and
  made the run take days; the default population path no longer pays that cost while
  the report is byte-identical.)
- **License: CRISPRme+ is now the MGB Open Access License 1.0 (MGBOA 1.0).** The
  software moves from AGPL-3.0 to MGBOA 1.0 (non-commercial, non-revenue-generating
  academic use; commercial use requires a license). `LICENSE` carries the full MGBOA
  text; the README, GitHub Pages footer, and the generated-report footer are updated
  to match (the report footer uses the MGB-approved "License and permitted use" +
  disclaimer language). The main `pinellolab/CRISPRme` repo is relicensed to match.
- **Report: the off-target Mismatch/Bulge matrix is now resolved by bulge TYPE
  (DNA vs RNA).** The single matrix splits each row into DNA-bulge and RNA-bulge
  counts (derived from the guide-aligned gap positions), so off-targets are countable
  by bulge type; the coarse total-bulge-only matrix is removed (the DNA/RNA table
  includes the 0-bulge row, so it is a strict superset and keeps the perfect-match
  highlight).
- **Report: variant identifiers are shown as `chrom;pos;ref;alt`** (with the rsID as a
  companion), instead of preferring the rsID — normalizing the previously mixed
  rsID / `pos:ref:alt` / `chr_pos_ref_alt` forms; multi-SNP haplotypes are joined.
- **Report: the high-variant-density column now reads "N in window"** (was
  "Yes (N var)"), and the legend spells out that this counts total variants in the
  surrounding window — distinct from the `Variant` column, which lists the variant(s)
  that create the individual off-target.
- **Personal-assembly storage is now one self-contained folder per individual with a
  single `metadata.json`.** `Assemblies/<individual>/{paternal,maternal}/genome/*.fa`
  plus each haplotype's `*_vs_GRCh38.chain.gz` and `*.chromAlias.txt`, described by a
  relocatable `metadata.json` (paths relative to the folder). This replaces the flat
  layout that scattered files across `Genomes/` + `LiftoverFiles/` with six hidden
  `.assembly_individual` marker files. Benefits: deleting an individual removes one
  folder and reclaims all its space (the old delete only removed markers), clear
  separation between individuals, and portable bring-your-own bundles. The old flat
  layout is still discovered for backward compatibility. New shared module
  `PostProcess/personal_assembly.py` (+ `test_personal_assembly.py`).

### Added
- **`assembly-search` annotates reconciled sites against hg38.** New optional
  `--annotation` flag (same file format, COSMIC licence gate and IntOGen handling as
  `complete-search`'s) adds an `Annotation` column to `combined_hg38.tsv` for every site
  with an hg38 coordinate (found on both haplotypes, or on one and mappable); sites with
  no hg38 equivalent get none. The web form applies the annotations enabled in
  Settings -> Annotations, like `complete-search`. In the results table the column is
  filterable and sortable like every other one. The annotation describes the *lifted*
  hg38 locus, not the haplotype's own sequence.
- **`scorer-env` command group + modular scorer environment.** `crisprme.py scorer-env
  {create,check,update,list,doctor}` manages the dedicated `cbulge` conda env (detects
  micromamba/mamba/conda, provisions the pinned CRISPR-Bulge source+weights, health-checks
  imports/weights, records state under `Annotations/.scorer_env.json`). The Docker image
  and `install_from_source.sh` create it at install time; set `build_scorer_envs=0` /
  `CRISPRME_SKIP_SCORER_ENV=1` to defer. `CRISPRME_SCORER_ENV` points CRISPRme at an
  existing/shared scorer env instead of creating its own.
- **`--compute-backend {cpu,gpu,auto,cuda,metal}` for the ML scorer** (default `cpu`;
  GPU is opt-in and always degrades gracefully to CPU). `cuda` uses a CUDA TensorFlow build
  on NVIDIA/Linux; `metal` uses tensorflow-metal on Apple Silicon **guarded by a load-time
  numerical self-test** — tensorflow-metal silently miscomputes this model's GRU kernel, so
  the scorer detects the miscompute and transparently falls back to a CPU device context
  (bit-identical scores) rather than emit wrong ones. A CPU thread cap
  (`CRISPRME_SCORER_THREADS`, default `min(16, cores)`) prevents TensorFlow oversubscription
  on many-core nodes. **CUDA is validated on an NVIDIA A100** (Linux): the GRU is numerically
  correct (~1e-5 of CPU) and scoring runs ~14.7k off-targets/s — ~1.6× the CPU rate, so a
  genome-wide worst-case (6.77M off-targets) scores in ~8 min. The Docker image accepts a
  `--build-arg scorer_backend=gpu` to bake in the CUDA scorer env (run with
  `docker run --gpus all … --compute-backend cuda`); CFD stays on the CPU alongside it.
- **New recommended default variant index: single-source, phased
  `NRG_3_hg38+hg38_1000G2021`** (published to HuggingFace). 1000 Genomes 2021,
  3,202 samples, fully phased (`data_type=genotyped-phased`) — so every reported
  co-occurrence is CONFIRMED cis with named per-sample carriers throughout (no
  PUTATIVE-from-phasing): a clean, simple default of real observed haplotypes.
  `NRG_3_hg38+hg38_1000G2021_HGDP` remains available as the broader-coverage option
  (adds HGDP's 929 individuals; hybrid confidence). The docs and CLI download hints
  now point at the single-source default first.
- **`--alt-alignments` / `--no-alt-alignments`** (`complete-search`) — explicitly
  emit or suppress the exhaustive alternative-alignments file, overriding the
  mode-driven default (off in population-level, on under `--per-sample`).
- **Index `data_type` / `phased` metadata in the registry manifest.** Each
  `reg_<chrom>.idx` now records `data_type` (`sites-only` / `genotyped-unphased` /
  `genotyped-phased` / `hybrid`) and a `phased` flag (per-database + overall),
  observed from the genotypes at build time — so tooling can identify an index's
  type without scanning the multi-GB per-sample dictionary. The variant-dataset
  selector in the web UI reads it to label each dataset (falling back to
  genotype-store presence for indexes whose manifest predates the flag). The three
  shipped indexes (`1000G2021_HGDP` = hybrid, `HPRC` = genotyped-phased, `mega` =
  sites-only) were backfilled in place on HuggingFace.
- **Report: an `Observed` column flagging off-targets supported by ≥1 real
  individual.** Values are `reference` (present in the reference genome, so carried
  by essentially every individual), `N carrier(s)` (a variant site carried by N
  named individuals in a genotyped panel), `observed` (a single AF>0 variant —
  carried by ≥1 individual by definition, even where carriers can't be named), or
  `putative` (a multi-variant combination whose cis co-occurrence cannot be
  confirmed). On a genotyped index the carrier roster is populated in both analysis
  modes; under `--per-sample` it is pruned to the EXACT cis carriers, so the column
  cleanly separates genuine haplotypes from worst-case reconstructions. The
  recommended validation panel now **prioritizes observed (real-carrier) sites**
  among equally-severe candidates — as the first tie-break *after* worst-case
  severity, so a strictly-worse reconstructed site is never displaced (worst-case
  coverage preserved), and the panel summary reports how many selected sites are
  observed. On a sites-only panel (e.g. `mega`) the `N carrier(s)` value and the
  tie-break are inactive (no per-sample roster) but the column is still shown; it is
  dropped only for a reference-only run.
- **`assembly-search --assembly-individual <name>`** — resolve both haplotypes'
  genome/chain/chromAlias from a registered `Assemblies/<name>/` bundle instead of
  passing the six paths explicitly (the six flags still work, incl. the legacy layout).
- **`assembly-search --max-total-edits <N>`** — override the total mismatch+bulge
  budget for the per-haplotype searches. The default is unchanged (the full
  `mm + bDNA + bRNA` budget — a personal assembly is a plain resolved genome with no
  IUPAC-variant lattice, so it does not need `complete-search`'s variant-density cap of
  4), but the override lets you tighten or widen it explicitly. Also fixes reuse
  detection so a re-run with the same parameters reuses the prior per-haplotype search.
- **`assembly-search` run history + self-describing jobs.** Each assembly-search job
  now records its parameters (individual, guide, PAM, edit budget, mode) alongside the
  results, and the web **History** page lists past assembly-search runs next to the
  standard searches so they can be reopened.
- **Bring-your-own personal assembly** — the "Add a personal assembly" Settings card has a
  file-picker to upload a `.zip`/`.tar.gz` bundle (one top-level `<individual>/` folder + a
  valid `metadata.json`); it is extracted, validated, and registered
  (`personal_assembly.import_archive`). Extraction is path-traversal-safe (rejects
  absolute / `..` members and tar symlink/device entries, and verifies every member stays
  within the staging directory).
- **Single canonical installer** — the clickable-app `install.sh`/`install.ps1` now live
  only in the served GitHub Pages repo (`pinellolab/CRISPRme` `docs/`); the drift-prone
  duplicate under `install/` was removed and replaced with a pointer README. Docs clarify
  that only the whole-genome variant-aware search needs ~64 GB RAM (many current laptops
  have it); reference-only / single-chromosome / region searches fit ~16 GB.
- **Direct haplotype-vs-haplotype reconciliation for `assembly-search`.** A
  subset of each haplotype's non-mappable (no hg38 equivalent) off-target
  predictions are now resolved by aligning the two haplotype genomes
  directly to each other (`minimap2`/`impg`), independent of `liftOver`.
  A prediction is only merged into the new `both_haplotype_private`
  category when both a single unambiguous alignment hit and mutual
  (bidirectional) reciprocal confirmation hold — real, measured directional
  asymmetry between haplotypes required this over a simpler one-directional
  check. `minimap2`/`impg` are new required dependencies for
  `assembly-search` (`check_impg_available()` mirrors the existing
  `check_liftover_available()` gate). 92/92 tests, real-data verified. The
  live results page and the static combined report both surface the new
  category: the "Both haplotypes" haplotype-coverage bar now has a real
  unmapped segment (previously hardcoded to 0, since the category didn't
  exist yet), and the live page's summary stats gained a "Found in both
  haplotypes, non-mappable to hg38" card. Reports from before this
  category existed render unaffected (the new count defaults to 0).

### Removed
- **The legacy vendored RandomForest ML scorer and azimuth are removed.** The 276 MB vendored
  RandomForest model (predictors archive, `dnaShape.pkl`, and its scorer module), the whole
  `PostProcess/azimuth/` tree (already inert / never invoked), and their now-unneeded
  scikit-learn 1.1.3 / numpy 1.24 pin are deleted. CRISPR-Bulge is now the ML scorer (see
  *Changed*); CFD's small lookup pickles (`mismatch_score.pkl`, `PAM_scores.pkl`) are
  version-agnostic and unaffected. `--scorer` is gone (the scorer is always CRISPR-Bulge); the
  internal `CRISPR_BULGE_*` column names are kept for back-compat only.
- **Legacy per-file personal-assembly upload in Settings.** The old "Add a personal
  assembly" card that took the six individual files (two genomes, two chains, two
  chromAlias) one at a time is gone; personal assemblies are now added as a single
  self-contained `.zip`/`.tar.gz` bundle (see *Bring-your-own personal assembly*
  above), matching the one-folder-per-individual storage layout. The CLI still accepts
  the six explicit paths for scripted use.

### Fixed
- **A genome delivered as one multi-sequence FASTA file no longer silently produces an
  EMPTY result.** CRISPRme identifies chromosomes by FASTA *file name*, one sequence per
  file (as the README requires), but CRISPRitz happily indexes and searches a single
  multi-sequence file -- so a UCSC assembly with no per-chromosome archive (e.g. the pig
  `susScr11`, `bigZips/susScr11.fa.gz`) ran to completion and then reported "No
  off-targets found" even with tens of thousands of raw hits, because the whole genome
  was treated as one "chromosome" named after the file. Now: `download` /
  the Data Manager split such a download into one `<sequence>.fa` per sequence; the
  input validator and the pipeline itself refuse a multi-sequence genome folder up front
  with the exact fix (`python PostProcess/genome_layout.py --split <genome_dir>`, which
  also repairs an existing folder -- rebuild any index made from the old layout); and the
  pipeline no longer reports success when the search found raw hits but none survived
  post-analysis (it fails with an explanatory error instead). New
  `PostProcess/genome_layout.py` (+ `test_genome_layout.py`); the layout check reads one
  header line per file, so it adds no measurable time to a search.
- **A search that genuinely finds nothing now says so on the job status page** (green
  "Finished: no off-targets were found ...", later steps marked "Not needed") instead of
  leaving "Merging Targets" and the rest as "To do" indefinitely.
- **The Result Summary matrix no longer shows bulge rows `0..5` for a reference-only
  search that only reaches 2 bulges.** The per-guide count file always carries a
  reference block and an all-zero variant block; the page only split the two blocks for
  searches with variants. A reference-only job now shows just its reference block.
- Removed two Dash dev-mode prop-type warnings on the results page (invalid extra keys in
  a `DataTable` `css` entry; a hidden radio's options passed as `[labels, values]`
  instead of one list of `{label, value}`). No visible change.
- **`complete-test` no longer relies on on-demand index building.** The smoke test ran a
  variant `complete-search` that expected an index to be built on demand — a path that was
  gated off (it produced incomplete dict-less tiers), so `complete-test` failed. CRISPRme does
  **not** build an index automatically: `complete-test` now checks for a prebuilt index and, if
  none is present, prints clear guidance to **download** one (`crisprme.py download --what index
  --index-name …`) and exits cleanly without downloading test data or building anything — the
  user decides which index(es) to fetch (as the web instructions already recommend). When an
  index is present it is used via `--index-path`.
- **Web interface now boots on the modernized environment.** The Dash cache was configured
  with the Flask-Caching **1.x** short-name `CACHE_TYPE: "filesystem"`; Flask-Caching **2.x**
  (pulled in by the un-pinned, modernized web stack) removed the short-name aliases and expects
  the backend class name, so `filesystem` raised `ImportStringError` at import time and the whole
  web app failed to start (`crisprme.py web-interface` and the served site). `CACHE_TYPE` is now
  `"FileSystemCache"`, which resolves on Flask-Caching 2.x (and remains valid on 1.10+).
- **Web results page no longer crashes (or silently drops a plot) when the CRISPR-Bulge
  filter is selected.** The filter dropdown value is the hyphen form `CRISPR-Bulge`, but the
  cached top-1000 / personal / private lolliplot files are written with the underscore token
  `CRISPR_BULGE`; the results page interpolated the hyphen value straight into the filename, so
  the personal-card page raised `ValueError('Personal and Private Lolliplots not found')` and
  the per-guide top-1000 plot vanished. The filename token is now normalized
  (`CRISPR-Bulge`→`CRISPR_BULGE`); `CFD`/`fewest` are unaffected.
- **The report now states, as a live count, how many off-targets are ≥2-bulge (N/A for
  CRISPR-Bulge).** The recommended-panel note reads "N of M off-targets are ≥2-bulge and
  therefore show `-` (N/A) in the CRISPR-Bulge column — ranked by CFD and edit distance
  instead, not dropped," so a mostly-empty ML column on a bulge-heavy run reads as expected
  behavior rather than missing data.
- **Variant search no longer aborts during report generation** (population-distribution
  step). `process_summaries.py` writes the per-criterion sidecars using the filter-criterion
  token (`..._CRISPR-Bulge.txt`, hyphen), but `submit_job` read the population-distribution
  input — and error-path cleanup globs — with the underscore identifier form
  (`..._CRISPR_BULGE.txt`); the mismatch raised `FileNotFoundError`, wrote to stderr, and the
  `[ -s $logerror ]` guard aborted the whole run with "population distribution plots creation
  failed" on **every** search. The consumers now use the hyphen form that `process_summaries`
  actually emits.
- **`build-index-only --vcf` now accepts the same config-file form as `complete-search`.**
  `complete-search --vcf` takes a file listing VCF folder names (under `VCFs/`), but
  `build-index-only --vcf` accepted only a bare directory path — so the same flag meant two
  different things. `build-index-only` now accepts **either** a VCF directory (unchanged) **or**
  a one-line config file naming a folder under `VCFs/`, with a clear error if the listing names
  more than one dataset (build each separately / merge first). Help text + README updated.
- **`complete-search`'s "build the index first" hint now shows the correct flags.** The error
  told users to run `build-index-only ... --bMax <N>`, but that command reads `--bDNA`/`--bRNA`
  (and ignores `--bMax`), so following the hint verbatim silently built nothing. The message now
  says `--bDNA <N> --bRNA <N>`.
- **Report ≥2-bulge N/A wording is precise.** The scores legend said such off-targets are
  "shown as -1"; the report tables actually show `–` (only the raw `integrated_results.tsv`
  records `-1`). The legend now states both.
- **`assembly-search` no longer crashes when a haplotype has zero off-targets.** A
  haplotype whose search finds no hits on the assembly is a legitimate result, not an
  error; reconciliation now treats it as an empty set — every locus on the other haplotype
  becomes heterozygous-equivalent (`<other>_only`), and the empty side contributes 0
  non-mappable — instead of aborting on the missing `*_integrated_results.tsv`. Unsupported
  multi-guide runs (more than one results file) still error.
- **`assembly-search` haplotype reconciliation reads only `*_integrated_results.tsv`.**
  It previously also outer-merged in `*_all_results_with_alternative_alignments.tsv`
  and re-clustered by CFD, on the theory that `integrated_results.tsv` alone
  was missing real distinct sites. Re-derived: that theory was wrong — the
  union's "duplicates" (up to 233 rows at one locus) were alt-file-only
  representations of an already-counted site, not additional real
  off-targets, and `integrated_results.tsv` already has exactly one row per
  real physical cluster by construction. The representative-selection
  criterion is also now consistently fewest-mismatches+bulges (matching
  `merge_contiguous_targets.py`'s own choice) rather than highest-CFD at
  every call site. Real genome-wide effect on a real dataset: net one fewer
  `both` row (a duplicate collapsing), zero effect elsewhere. 58/58 tests.

## [2.5.5] - 2026-09-13

### Added
- **HPRC pangenome production index** (`NRG_3_hg38+hg38_HPRC`) — a third prebuilt SpCas9
  NRG variant index, published on HuggingFace alongside the genotyped 1000G-2021 + HGDP and
  the sites-only mega. It enriches hg38 with the HPRC Minigraph-Cactus pangenome decomposed
  to a per-chromosome VCF (`vg deconstruct` → `vcfbub` → `vcfwave`; 232 assembly-derived
  genomes incl. CHM13). Genotypes are **phased**, so it supports CONFIRMED cis co-occurrence
  and named per-sample carriers, and it captures pangenome-specific variation absent from
  short-read panels. Reproducible recipe: `seq_script/merge_panels/hprc_build.sh` (RAW
  registry + SNP+indel, consistent with the other two indexes). METHODS §merge documents it.
- **`assembly-search` web interface** (PR #48) — the diploid personal-genome search
  (paternal + maternal assembly FASTAs, liftOver-reconciled to hg38; no VCF) is now available
  on the website: a "Personal assembly" tab under Select genome, a parallel status-polling
  `/load` page, and a dedicated results view (`result_page_assembly`) that mirrors
  complete-search's layout for two-haplotype data (haplotype-coverage stats, reconciled
  off-target table, per-haplotype non-mappable sites, Graphical Reports). A combined
  report.zip is produced at parity with complete-search (`build_combined_report`), and a
  Settings "Add a personal assembly" card fetches an HPRC Release-2 individual
  (`download_hprc_assembly.py`) or accepts uploaded assemblies.

### Changed
- Pinned **`dash-bootstrap-components<2`** in both `environment.yml` (PR #48) **and the
  `Dockerfile`** — the Docker build installs deps directly (not from `environment.yml`), so the
  image pin is what actually reaches users. dbc ≥ 2.0.0 requires React 18 / Dash ≥ 3; this app
  targets React 16 / Dash 2.x, and the assembly-search `dbc.Tabs` do not render under dbc 2.x.
  The published v2.5.5 image ships dbc 1.7.x so the web assembly-search feature renders (the CLI
  is unaffected).

### Fixed
- `assembly-search` now passes `--max-total-edits {mm+bDNA+bRNA}` into each haplotype's
  `complete-search` (PR #47) — previously it silently used the default cap of 4, pruning
  mismatch+bulge off-targets.
- `download_hprc_assembly.py`: the `--register` path used a bare `from pages_utils import …`
  (→ `ModuleNotFoundError`, the only mode the Settings card invokes) — fixed to
  `from pages.pages_utils import …`; and the multi-GB scratch `.fa.gz` is now removed after a
  successful split and on MD5-mismatch abort (was leaking ~2 GB per haplotype).

### Verified
- HPRC index validated end-to-end (`--per-sample`): 68,057 CONFIRMED phased-cis off-targets,
  855 SNP+indel and 140 SNP+SNP co-occurring rows, with named cis carriers ground-truth-verified
  against the source genotypes (e.g. HG03579 carries chr10:12896898 `TG>T` + chr10:12896900
  `T>A` in cis, joint AF 1/464). Assembly-search validated e2e (HG01255 chr19, merged code):
  reconciliation + combined report.zip. Full web click-through: tabs render (dbc<2), Load
  Example, Per-sample gray-out survives the assembly-tab round-trip, Settings cards, assembly
  results page. HPRC download card exercised against real HPRC S3 (92 contigs, register, leak-free).

## [2.5.4] - 2026-09-11

### Added
- **IntOGen (CC0) cancer-driver genes as the default cancer-gene annotation.** With COSMIC
  now licence-gated (below), CRISPRme+ ships **IntOGen** — a compendium of cancer *driver*
  genes whose current release is **CC0** (public domain; free for academic AND commercial
  use + redistribution) — **enabled by default** as the licence-free cancer-gene flag. New
  report column `Annotation_INTOGEN` (mirrors `_COSMIC`), a default-on `intogen_drivers.hg38.bed.gz`
  track (`seq_script/build_intogen_annotation.py` maps the CC0 driver symbols to GENCODE gene
  spans), report legend, Settings label. So a default search still flags cancer-driver
  off-targets — legally — and users with a COSMIC licence add the curated Census on top.
- **COSMIC cancer annotations are now licence-gated and EXCLUDED by default.** COSMIC
  (Genome Research Ltd / Wellcome Sanger) is free for academic / non-commercial research
  but **commercial use requires a licence** ([terms](https://www.cosmickb.org/terms/)).
  Because the built-in hg38 bundle bakes COSMIC in, CRISPRme+ now **strips COSMIC features
  from the annotation (and therefore the output) by default**, on both the web and the CLI.
  The user opts in **once** (the choice persists in `Annotations/.cosmic_license.json`):
  a **Settings → "COSMIC cancer annotations (licence)"** checkbox on the web (with a link to
  the terms), or `crisprme.py cosmic-license {enable|disable|status}` on the CLI (`enable`
  asks for confirmation; `--accept` for non-interactive/Docker use). New module
  `PostProcess/cosmic_license.py` (attestation flag + row-level `_COSMIC` filter);
  `build_active_annotation` and `complete-search` consult it. (`PostProcess/test_cosmic_license.py`.)

### Changed
- **Renamed the analysis-mode control: default is now "population-level"; `--full` → `--per-sample`.**
  The measurements showed the two modes differ in **genotype resolution**, not intrinsic speed
  (on a single guide the two take the same wall time; on a sites-only panel the opt-in mode does
  nothing). So "fast/full" was the wrong framing. The CLI now runs a **population-level analysis**
  by default (no flag) and exposes an opt-in **`--per-sample`** that resolves per-sample genotypes
  into observed haplotypes (CONFIRMED cis + named carriers + exact joint AF). The web form's mode
  control is relabeled **"Analysis mode: Population-level (default) / Per-sample"**, and the
  report's Search-mode note + the `.search_mode` marker use the new names. `--per-sample` requires
  a **genotyped** panel: on a sites-only index it is a **no-op with a warning**, and the web form
  **disables** the per-sample option when a sites-only index is selected. The internal
  `CRISPRME_FAST_MODE` env var (post-analysis plumbing) is unchanged and not user-facing.

### Removed
- **`--fast` / `--full` flags removed (clean break, no aliases).** They shipped only in 2.5.3 and
  are replaced by the population-level default + `--per-sample` above. Legacy `.search_mode`
  markers (`fast`/`full`) written by 2.5.3 are still recognized when rendering old result dirs.

### Verified
- **New-user clean-room on the released v2.5.3 (image + both published RAW indexes).**
  `apptainer pull docker://pinellolab/crisprme:v2.5.3` reports `v2.5.3`; a fresh
  `download` of **both** production indexes confirms the Tier-0 registry ships
  **uncompressed (`codec=0`, RAW)** as intended (`registry_hg38_1000G2021_HGDP` and
  `registry_hg38_mega`); a fast-default `complete-search` on **each** completes `rc=0`
  and produces a `report.zip` carrying the **Search-mode note + fast caveat**. Confirms
  the two-index design end-to-end: genotyped (CONFIRMED-capable) vs mega (sites-only,
  PUTATIVE + min-AF), and that the **mega now searches indels** (32-bin fake-indel
  genome; 387 SNP+indel co-occurrence + 393 `indel_af` rows with 5-source per-dataset AF).
- **Measured fast-vs-`--full` comparison on both indexes (adversarially verified).**
  Genotyped 1000G-2021+HGDP: detection is **lossless at the window level** (fast a superset;
  the 24 "full-only" loci are 1–10 bp anchor shifts of windows fast detected); **`fast_CFD ≥
  full_CFD` with 0 violations across 3,079 shared loci**; `--full` **tightens** carrier sets
  (fast over-lists a PUTATIVE union → 1,642 distinct samples; full prunes to observed cis →
  1,447) and turns 0 CONFIRMED co-occurrence rows into **106 CONFIRMED phased-cis** rows with
  exact joint AF + named carriers. Sites-only mega: `--full` yields **0 CONFIRMED / no
  carriers** in either mode and the `indel_snp_cooc` / `indel_af` companions are
  **byte-identical** — so `--full` is meaningful only on a genotyped panel. Single-guide
  runtimes are near-equal (fast ~13.8 min vs full ~13.4 min); fast's speedup is a
  dense/aggregate-panel effect, not a per-guide one.

### Fixed
- **New-user docs pinned stale artifacts.** `README.md` copy-paste commands referenced the
  previous image `v2.5.2` → bumped to `v2.5.3`; `docs/DOCKER_QUICKSTART.md` was a v2.4.0-era
  doc with a **dead index name** (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, whose download would
  fail) → fixed to `NRG_3_hg38+hg38_1000G2021_HGDP`, image tags bumped to `v2.5.3`, the mega
  index added, the "Alpha release" banner replaced, and a fast-vs-full explanation added.

### Documentation
- **METHODS §5/§8 + README** now document the *measured* fast-vs-`--full` behavior above
  (window-level losslessness, the exact worst-case CFD bound with 0 violations, carrier
  tightening, and the verified sites-only inertness of `--full`), so the fast/full contract
  is documented from data rather than design intent alone.

## [2.5.3] - 2026-09-11

### Performance
- **All indexes ship the Tier-0 registry uncompressed (raw) by default — "speed over space."**
  Registry lookups (millions per dense search, scattered across a chromosome) are now direct
  memory-mapped reads with **no per-lookup decompression** — measured **~2× faster** than the
  compressed path even with a warm cache. The genotyped build previously hard-coded zlib
  (`compress=True`); it now defaults to raw, matching the sites-only mega index, so both
  production indexes use **one consistent format**. Optional zlib block-compression (~3.5×
  smaller on disk) remains available via `CRISPRME_REGISTRY_COMPRESS=1` at build time (or
  `transcode_registry`); the reader reads either format byte-for-byte identically. The block
  cache + O(1) LRU below now serve that opt-in compressed path.
- **Tier-0 registry block-cache: ~3.7× faster dense variant post-analysis, genome-wide, byte-identical.**
  The genotyped v3 registry is zlib block-compressed (`~3.6×` smaller on disk); profiling a
  dense per-contig SNP post-analysis showed the cost had moved to lookup-time decompression —
  **`zlib.decompress` was 71% of the run** — because the decompressed-block LRU was only 8
  blocks while a dense IUPAC decomposition touches positions across the whole chromosome, so
  hot blocks were re-decompressed ~10⁶ times. The cache default is raised **8 → 4096 blocks**
  (tunable via `CRISPRME_REGISTRY_CACHE_BLOCKS`); the LRU only grows to the blocks actually
  touched, so a block being ~64 KB, this "caches the whole contig registry" for any single
  contig. **Measured** on a dense chr22 SNP post-analysis: **326 s → 89 s (3.66×), output
  byte-identical** (same md5). Sized from measured block counts — genotyped chr22 = 339
  blocks (~22 MB), chr1 (largest) = 1,868 blocks (~122 MB) — so **4096 holds every genotyped
  contig fully**, extending the speedup genome-wide (not just small chromosomes), and the cap
  bounds a pathological panel at ~256 MB/reader. The **mega sites-only index is uncompressed
  (`codec=RAW`) → unaffected** (harmless no-op). Also corrects the earlier "the ML scorer is the
  tail" assumption — the ML scorer is ~4% of the run; the registry decompression was the real tail.
- **O(1) registry block-cache LRU (`OrderedDict`).** Profiling the *indel* post-analysis with
  the enlarged cache exposed a second tail the small cache had hidden: the LRU was a plain
  list, so every cache hit did an O(cache-size) `list.remove()` — **43 s / 26 %** of a dense
  indel post-analysis at ~15 M hits. Replaced with an `OrderedDict` (`move_to_end` /
  `popitem`, all O(1)); byte-identical LRU semantics. **Indel post-analysis 167 s → 121 s**,
  and it also roughly halved the SNP path (cache=512: 89 s list-LRU → 43 s O(1)). A clean
  cache sweep shows the hot working set is tiny (plateau at ~128 blocks), so the combined
  cache + O(1)-LRU fix is **326 s → 43 s ≈ 7.6×** on a dense chr22 SNP post-analysis, and the
  4096 default is generous headroom that fills only to touched blocks (~8-40 MB in practice).

### Changed
- **Fast mode is now the DEFAULT search behavior; `--full` opts into exact enumeration.**
  The SNP variant post-analysis now reports one worst-possible representative per variant
  window by default (exact worst-case CFD; the ML scorer is a best-effort screen), instead of
  enumerating every observed haplotype — so a search stays tractable on dense / aggregate
  panels out of the box (the full enumeration was measured 49 h+ without completing on a
  4×-density panel). **Trade-off: per-sample carriers, CONFIRMED cis phasing and exact joint
  allele frequency are NOT computed in fast mode** — pass `--full` (CLI) or pick **Full** in
  the web search form for that per-sample resolution (recommended for genotyped panels /
  clinical validation). `--fast` is retained as a deprecated no-op alias. The web search form
  exposes a **Search mode** control with an inline explanation, and the shareable report now
  carries a **Search mode** row stating which mode was used and its implications. SNP+indel
  co-occurrence is unchanged by the search mode. A launch-time message and the report note
  make the carrier trade-off explicit (never silent). Mechanics unchanged: the resolved mode
  propagates through the post-analysis subprocess tree via `CRISPRME_FAST_MODE`
  (`0`=full, `1`=fast); a direct `new_simple_analysis.py` call still defaults to full, so unit
  tests and the byte-identical guarantee for the enumeration path are unaffected.

### Added
- **Parallel ML scoring (opt-in, default OFF)** — prototype
  that removes the per-contig ML-scoring tail (one large chromosome's worker running
  serially long after the others finish). The per-contig post-analysis pool is unchanged; this
  parallelizes the CPU-bound per-target **feature build** (`get_features`) — the genuinely
  serial part — at the finest seam across a small bounded **fork** pool.
  The model stays in the parent (children only build features). **Bit-identical feature
  matrix by construction** (a feature row is a pure per-target function; contiguous,
  order-preserving chunks) and default OFF ⇒ output + process count unchanged. Fork (not spawn)
  because the post-analysis entry points run their main at module top level. Measured
  feature-build speedup ≈ 3.8× at 8 workers. The model's
  `predict` is left untouched — the pickled predictors already carry `n_jobs=-1` (all cores), so
  it is already parallel (and, as a result, non-deterministic at the raw-float level: two serial
  runs differ ~5e-17, identical to the emitted 3 decimals — the ML scorer's output has always been
  reproducible only to the emitted precision, which this preserves).

### Fixed
- **HuggingFace download resilience against rate limits (HTTP 429).** A first-time,
  unauthenticated `crisprme download` of the many-small-file indices could hit HF's
  "Too Many Requests" and abort the whole download. The downloader now retries transient
  429 / 5xx responses with exponential backoff (honoring a sane `Retry-After` header) around
  both the snapshot fetch and the size-verification pass; `snapshot_download` resumes, so a
  retry only re-pulls what's missing. 4xx / auth errors still fail fast. Tunable via
  `CRISPRME_HF_MAX_RETRIES` (set to `1` to disable). Authenticating with an HF token still
  avoids the limit entirely.

## [2.5.2] - 2026-09-09

### Added
- **Two-pass fast mode (`complete-search --fast`, opt-in)** for dense / aggregate variant
  panels where per-haplotype enumeration is intractable (measured **49 h+ without
  completing** on a 4×-density 1000G+HGDP panel). Instead of enumerating the 2^k IUPAC
  haplotype lattice per window, it emits a small set of **worst-possible representatives**
  (reference + minimum-edit + maximum-CFD), propagated to the whole post-analysis via
  `CRISPRME_FAST_MODE`. Validated **lossless for locus detection** and **non-understating
  for the worst-case score** against the slow full-enumeration path on a real chr22
  1000G-2021+HGDP slice (0 CFD under-reports; surfaces *stronger* worst cases at 182 loci),
  and **confirmed genome-wide** (V2 non-`--fast` vs V3 `--fast` on the 2021 panel: V3 is a
  locus-level superset, **0 of V2's 1,458 CFD≥0.2 loci lost or demoted**, CFD exact-or-
  conservative, the ML scorer screen-grade only near the 0.2 line). The default (non-`--fast`) path
  is **byte-identical**. See `docs/DESIGN_2.5.1_two_pass_fast_mode.md` and METHODS §5/§8.
- **All-source "mega" sites-only index (5 datasets).** A new merged panel — 1000 Genomes
  2021, HGDP, gnomAD v4.1, TOPMed, All-of-Us — built directly from each source's aggregate
  `INFO/AF` (these carry no shared samples to reconstruct cross-source haplotypes), with
  per-dataset AF plus a per-site **`AF_max`** (the maximum across sources) as the reported
  global frequency. `merge_mega_sites.sh` + `build_mega_registry.py` compile the Tier-0
  registry directly from the frequencies (byte-identical binary; SNPs, with an indel-AF
  companion store), and the enricher tolerates sites-only VCFs. See METHODS §2/§3.
- **Report companions for aggregate panels.** The report now bundles the per-dataset
  **`indel_af.tsv`** companion (indel allele frequencies by source) alongside the SNP+indel
  co-occurrence TSV, and renders carrier / homozygote frequencies as **"NA"** for aggregate
  (sites-only) groups — where only allele frequencies exist and per-individual carrier counts
  are undefined — instead of fabricating numbers.
- **`download` always ships the shared reference index** alongside a variant index, so a
  reference-genome scan (or an on-demand variant-index rebuild) works immediately after a
  download with no separate build step.
- **SNP+SNP co-occurrence companion (`<output>.snp_snp_cooc.tsv`).** The SNP-side analogue
  of the SNP+indel co-occurrence companion: one joinable row per emitted variant off-target
  that USES ≥2 co-occurring SNP alt alleles, recorded from BOTH the observed-haplotype
  enumerator (genotyped panels → `CONFIRMED`/`PUTATIVE`, exact phase) and the registry-only /
  capped finalizer (sites-only panels → `PUTATIVE`). Reports the phase, a conservative
  `MinAF_bound` (the minimum participating marginal AF — a valid upper bound on the joint cis
  AF for both phases), and the observed carriers (`NA` on the sites-only path). Gated on a
  variant registry being present, so a legacy dict install is byte- and allocation-identical.
- **SNP+indel PUTATIVE co-occurrence on sites-only panels (mega).** The SNP+indel
  co-occurrence emission previously wrote a row only when cis carriers were countable, so on
  an aggregate sites-only panel (no per-sample genotypes) an indel + used SNP alt that
  co-occur at a locus produced no row. It now emits a `PUTATIVE` row with the conservative
  `min`-AF joint bound (over the indel `AF_max` + each used SNP's registry AF) and `NA`
  carrier fields. Gated on genotype ABSENCE, so a genotyped panel's real no-cis-carrier call
  is never turned into a phantom row.
- **Lossless dense-region worst-case haplotype (`CRISPRME_LOSSLESS_DENSE`, default ON for
  sites-only).** In a capped dense window the min-mismatch representative can be a strict subset
  of a genuine carried haplotype, so a sites-only off-target requiring ≥4 co-occurring variants
  on one haplotype could be dropped. The registry-only (sites-only) path therefore also emits
  the full co-located variant union — the maximal PUTATIVE haplotype — as an extra
  representative, bounded (not the 2^k lattice) and gated by the finalizer's own mm/PAM budget
  so nothing over-budget or PAM-invalid is emitted. The effect is **scoped to sites-only**
  (`registry_only_mode`), so it is **byte-identical for every genotyped / legacy / dict install**
  (those never enter the branch, and the genotyped path is already lossless via the observed
  enumerator). It is **on by default** so the sites-only panel fulfils "don't miss a region";
  set `CRISPRME_LOSSLESS_DENSE=0` to opt out (sites-only reverts to a greedy representative
  only). A per-sample union on the genotyped path is intentionally *not* done — it would risk
  phantom trans-as-cis haplotypes.

### Performance
- **ML scoring: load the model once + skip eager per-pentamer work.** The 276 MB
  RandomForest ensemble was re-`pickle.load()`ed on every scoring batch (both the SNP and
  INDEL post-analysis); it is now cached at module scope and loaded once per process.
  Separately, `get_features` no longer runs an eager per-pentamer `re.sub` + `random.choice`
  on N-free windows (guarded behind `if "N" in ...`), which also removes a latent
  nondeterminism from the hot path. Both are **byte-identical** at the reported score
  precision (validated against the shipped model); `get_features` micro-benchmarks ~1.6×
  faster. The dense-panel INDEL post-analysis remains dominated by single-threaded per-target
  ML-scorer compute (profiled scoring-bound, not enumeration-bound); parallelizing it is tracked
  separately.

### Fixed
- **Post-analysis no longer deadlocks at high thread counts.** A genome-wide variant
  post-analysis launched with a large `--thread` (e.g. 200 on a 200-core host) could hang
  indefinitely: a `multiprocessing.Pool` with hundreds of workers makes every worker inherit
  every other worker's internal sentinel pipe at fork, so the pool's join never finalizes
  (all per-chromosome tasks complete, the parent hangs in `pool.map` forever). Both the SNP
  and indel post-analysis pools now cap the worker count (default 32, override
  `CRISPRME_POSTPROC_MAX_WORKERS`) — no speed loss (post-analysis is I/O-bound, one
  subprocess per worker) and no deadlock.
- **Mega (all-source) index hardening.** Resolved 17 adversarial-review findings in the
  sites-only aggregate index path and de-duplicated colliding indel registry keys.
- **Clear pre-flight error when a variant search has no prebuilt index.** `complete-search
  --vcf` without a prebuilt index used to fall back to a legacy on-demand enrichment that
  builds a classic per-sample dict WITHOUT the dict-less tiers (registry/genotype) and
  could leave the indel coordinate map (`log_indels`) empty — crashing the variant/indel
  post-analysis ~30 min in with a cryptic `FileNotFoundError: log<chrom>.txt`. It now
  detects the missing tiers per dataset up front and stops immediately with an actionable
  message pointing to `download --what index` (recommended) or `build-index-only` first.
  Override for advanced/legacy on-demand builds with `CRISPRME_ALLOW_ONDEMAND_BUILD=1`.
- **Warn when SNP+indel co-occurrence is requested but the registry tier is absent.** The
  pre-flight guard accepts either the classic per-sample dict OR the dict-less registry, but
  co-occurrence needs the **registry** specifically. A classic-dict-only index (e.g. a local
  `build-index-only` run **without `--samplesID`**) previously ran cooc as a **silent no-op**;
  `complete-search` now prints a WARNING at launch naming the datasets missing the registry
  tier, so the missing co-occurrence output is no longer silent.
- **Fast-mode worst-case CFD is now exact on every path.** CFD is position-weighted, so the
  fewest-mismatch representative does not maximize CFD; fast mode now also emits the exact
  **maximum-CFD** representative (per-position argmax + bounded brute-force for the joint
  two-base PAM factor). On the **legacy dict / aggregate-panel (mega) path** this closes a
  **threshold-crossing** under-report (up to **0.14** CFD). Cross-checked against an
  independent factorized CFD oracle (4,000 random windows, joint-PAM case included) and a
  real CFD-scored multiallelic fixture.

### Notes
- In fast mode, **CFD is the exact worst case; the ML scorer is best-effort** (a non-factorizable
  RandomForest). On a chr22 1000G-2021+HGDP slice every ML-scorer ≥ 0.2 off-target was reported
  at full strength (under-reporting ≤ 0.04, sub-0.19 tail). **Genome-wide the ML-scorer tail is
  heavier:** ~5 % of ML-scorer ≥ 0.2 loci can drop below 0.2 under `--fast` (largest gap ~0.12),
  while **CFD had zero ≥ 0.2 losses** (V2-vs-V3 GW matrix). So `--fast` CFD is a safe
  actionable gate but **the ML scorer is a screen** — run **without** `--fast` for a guaranteed
  per-haplotype ML-scorer worst case / an ML-scorer action gate.
- **`--fast` scope + index/version compatibility.** `--fast` accelerates only the SNP
  post-analysis; the indel post-analysis is single-threaded and **unaffected** by it (and the
  `indel_snp_cooc.tsv` companion is byte-identical with/without `--fast`). The feature-on 2021
  index requires CRISPRme **≥ 2.5.0** — stock 2.4.0's indel post-analysis fails on it (2.4.0
  users use the 2019 / `-dictless` indices). The **mega** index currently searches **SNPs
  only** (its fake-indel genome is not built); searchable indels + PUTATIVE
  co-occurrence-without-genotypes (SNP+SNP and SNP+indel possible haplotypes from AF) are
  planned for **2.5.2**.
- **Multi-indel residual (documented):** the indel search materializes one indel per fake
  contig, so an off-target needing **≥ 2 co-occurring cis indels in one protospacer** is not
  a candidate (pre-existing single-indel-search property). Low-frequency: **~1–2%** of indel
  loci after repeat-masking/dedup; genuinely-missed off-targets **~0.1–0.2%**, all at the
  weakest (≈0-CFD) edit-budget ceiling.
- **`--fast` does NOT change SNP+indel co-occurrence.** `--fast` collapses only the *SNP*
  representative emission in `integrated_results.tsv`; the `indel_snp_cooc.tsv` companion is
  produced by the indel post-analysis' cis-phasing pass over the genotype tiers, which `--fast`
  does not touch. Measured on the completed genome-wide matrix (2021 panel, `--fast` vs
  non-`--fast`, guide TGCTTGGTCGGCACTGATAG): the two `indel_snp_cooc.tsv` files are
  **byte-identical** (same MD5, 2,729 rows, 843 CONFIRMED / 1,886 PUTATIVE, full per-sample
  `cis_samples` + joint-AF in both). Per-sample cis attribution is preserved in fast runs.
  (Corrects an earlier note that claimed `--fast` collapsed cooc to representatives only.)

## [2.5.0] - 2026-09-01

### Added
- **SNP+indel co-occurring off-targets — now enabled by default** (opt out with
  `CRISPRME_INDEL_SNP=0`; previously experimental/opt-in). Detects off-targets that
  require BOTH a nearby SNP **and** an indel on the same haplotype within a
  protospacer window — a class the classic two-pass search could not see. The build
  compiles a per-chromosome **phased indel genotype store** and overlays SNP IUPAC
  codes onto the fake-indel genome, so the indel search surfaces SNP+indel
  haplotypes; post-analysis reports **CONFIRMED-cis** (phased) and **PUTATIVE**
  (unphased) co-occurrences with per-sample carriers + joint allele frequency, and
  the shareable report bundles `indel_snp_cooc.tsv` and surfaces a
  "SNP + indel cis co-occurrences (N confirmed-cis)" section.
- **New 1000G-2021 + HGDP NRG variant index** published to HuggingFace
  (`NRG_3_hg38+hg38_1000G2021_HGDP`): 3,202 high-coverage phased 1000G samples +
  929 HGDP. `download --what index` yields a searchable, variant-ready dict-less
  install (registry + genotype + phased indel-genotype tiers bundled).
- Shareable report: publication-grade inline **vector SVG** figures.

### Changed
- **Dict-less observed-haplotype enumerator** replaces the legacy 2^k IUPAC lattice +
  greedy cap on the variant post-analysis path: every haplotype real individuals
  carry is emitted (sensitivity-first, **no truncation**), and dense/hypervariable
  windows are logged to `high_variant_density_regions.bed` for visibility rather
  than dropped — resolving the phased multi-SNP haplotype under-report on the
  default path. (Exhaustive worst-case coverage of *unobserved* haplotypes in dense
  windows remains a documented, flagged limitation; a lossless two-pass fix is planned.)
- `--max-total-edits` pruning now emits a WARNING and persists a `Pruning_note`
  into the run params so the "deep off-targets were pruned" signal reaches the
  report bundle.

### Fixed
- **Atomic per-chromosome build writes + integrity-checked resume** for the Tier-0
  registry, Tier-1 genotype, and phased indel-genotype stores: a SIGKILL/OOM
  mid-write no longer leaves a truncated file that a restart accepts as "done" (a
  silent-truncation class that had shipped a short chr1 indel-genotype store).
  Writers use tmp + `os.replace`; resume verifies integrity, not mere presence.
- Report: high-variant-density region count is **deduplicated** (was inflated by
  duplicate per-chromosome rows), and the co-occurrence TSV is deduplicated (the
  writer emits a row per alignment pass). Population-summary companion **degrades
  gracefully** (emits an NA-frequency row) instead of dropping multi-variant rows
  on a registry-only install.
- Pre-flight guard + build-time auto-fix so an over-listing `--samplesID` can no
  longer inflate the AN denominator (which deflated every reported allele frequency).

## [2.4.0] - 2026-08-20

### Added
- **Dict-less variant index is the default distribution.** The shipped SpCas9
  (NRG) variant index (`NRG_3_hg38-dictless+hg38_1000G_HGDP`, 1000G + HGDP)
  replaces the 152 GB per-sample SNP dictionaries with the compact Tier-0
  allele registry + Tier-1 genotype store — a ~16 GB download that searches out
  of the box with correct allele frequencies, rsIDs, per-individual Samples, and
  population summaries. See `METHODS.md`.
- **Shareable off-target report** (`generate-report`, auto-run at the end of
  `complete-search`): a single self-contained, branded `report.html` with an
  annotation / COSMIC legend plus a top-N table for validation. The web results
  page now shows a one-click **Download Report (.zip)** button in place of the
  three raw result-file links, and the CLI prints an "open the report" hint.
- **COSMIC cancer-gene annotation on by default.** The built-in functional-region
  bundle is now `dhs+encode_screenv4+gencode+cosmic.hg38.bed.gz` (ENCODE SCREEN v4
  cCREs + COSMIC Cancer Gene Census tiers); results carry an `Annotation_COSMIC`
  column. Legacy bundles remain supported via automatic fallback.
- **Observed-haplotype enumerator** (dict-less multi-variant path): multi-variant
  off-targets are enumerated only for haplotypes actually carried by sampled
  individuals (CONFIRMED = phased, PUTATIVE = unphased / mixed), removing phantom
  combinations. Gated on the Tier-1 genotype store.
- `download --what all` on hg38 also fetches the default reference index
  (`NRG_3_hg38`) so a reference-genome scan works with no on-demand index build.
- **Turn-key CLI search after a download.** `download --what index` now writes the
  `list_vcf.txt` / `list_samplesID.txt` list files (the ones the web form builds
  per-search) at the install root, so `complete-search --vcf list_vcf.txt
  --samplesID list_samplesID.txt` works out of the box — no hand-written config.
  A "Run a search from the command line" section is added to the Docker
  quickstart, README and the landing site, with the `generate-report` step.
- Manuscript-level `METHODS.md` (dict-less model, VCF homogenization / merging,
  AF estimation, phased / unphased / mixed haplotype scanning, variant-density
  control, annotations, report) and README "what's new" pointers.
- **Perfect-match (0-mismatch) multiplicity flagging.** A guide with more than one
  perfect genomic match (0 mismatch, 0 bulge) has no a-priori on-target -- each is
  an equally-efficient candidate cut site. The report **and** the web results page
  now force EVERY 0-mm site to the top of the validation panel (never truncated),
  add a `Perfect_match` column to the curated tables, and show a prominent red
  "no unambiguous on-target" banner listing the sites (an amber "presumed
  on-target" note when there is exactly one). Previously all 0-mm sites were
  bucketed as "the on-target" and excluded from the off-target set and the panel
  with no warning -- so a perfect-match off-target, the highest-risk class, was the
  one thing the report never surfaced.

### Fixed
- **Web job launch crashed on special characters in the job name.** A job name
  with shell / path metacharacters (e.g. `SBDS(T>C)`) was interpolated raw into a
  `mkdir` run through `/bin/sh`, aborting the launch (`Syntax error: "("
  unexpected`, HTTP 500 in `change_url`). The name is now sanitized to a
  filesystem/shell-safe charset, and the results directory + queue file are
  created without a shell (`os.makedirs` / `open`).
- **Combined allele-frequency denominator in the shipped dict-less index.** The
  registry's global panel size was built from an over-listed sample list
  (AN = 8858) rather than the VCF-resident panel (1000G 2548 + HGDP 929 = 3477 →
  AN = 6954), deflating combined MAF by ~27 %. The published registry now reports
  the correct combined AF (AC / 6954); allele counts were already correct. A
  singleton reads 1 / 6954 = 1.44e-4.
- **Annotation sorting is read-only-install + concurrent-job safe.**
  `_sort_annotation` / `_process_personal_annotation` build the sorted / merged
  annotation in a per-invocation temp dir instead of decompressing and
  overwriting files in the shared install `Annotations/` directory.
- **Dict-less unphased multi-variant off-targets: enumerate sub-combinations.**
  The observed-haplotype enumerator reported only the maximal-*union* variant set
  for an unphased individual, which hides a valid off-target when a variant
  *breaks* the target (e.g. one that disrupts the PAM) — the sub-combination that
  keeps the reference allele at that column (a real off-target under one possible
  phasing) was dropped. Unphased (and cross-phase-set) individuals now enumerate
  every non-empty subset of their own carried variants (never cross-individual
  chimeras); the mismatch/PAM gates prune out-of-budget subsets and the enumerator
  dedups shared subsets, bounded by `CRISPRME_IUPAC_CAP`. Recovers the missed
  off-targets with correct MAF + samples, with zero change to previously-reported
  rows or to the phased (confirmed-cis) path.
- **Registry-only (`--no-genotypes`) installs no longer drop variant off-target
  sites.** With a Tier-0 registry but no per-sample dict and no genotype tier, the
  variant finalizer's "≥1 carrier" guard silently discarded every variant-created
  off-target. Such sites are now emitted as degraded rows (Samples = `NA`) that
  still carry the creating-variant's registry rsID / AF, so the risk is surfaced
  without per-sample resolution. Gated on the registry-only mode, so legacy,
  registry+dict, and dict-less-with-genotypes installs are byte-identical.
- **`complete-test` sample-ID download.** The github fallback URL was built with a
  double slash (`…/data//samplesIDs/…`, a redirect the downloader could stumble
  on); it now normalizes to a single slash.
- **`download` prefers the corrected dict-less variant index.** Requesting a
  variant index by its bare (legacy dict-based) name now transparently resolves to
  the `-dictless` companion when one exists on the repo — both install to the same
  folder, so old links / muscle-memory get the smaller, combined-AF index.
- **Annotation sorting survives a small or full `/tmp`.** `_sort_annotation` /
  `_process_personal_annotation` now stage intermediates on a roomy writable disk
  (`$CRISPRME_TMPDIR`, else the working/output dir) instead of the default `/tmp`,
  which is tiny on many HPC/container setups.
- **Dict-less never misses a high-variant-density region.** In a capped dense
  window the observed path emitted only the real observed haplotypes (per-individual
  unions) and logged the region to `*.high_variant_density_regions.bed`; when those
  unions all exceeded the mismatch budget the region surfaced ONLY in the BED with no
  off-target ROW in the results (the dict path always reports a greedy representative
  there). The observed path now also emits the same greedy min-mismatch representative
  (argmin over the window's alleles) per capped region, so every dense region always
  surfaces ≥1 off-target row at parity with the dict path (verified genome-wide).
- **Dense-region transparency.** A new `High_variant_density_region` column in
  `integrated_results.tsv` flags each off-target that falls in a dense window,
  noting that a greedy worst-case alignment is reported and additional haplotype
  alignments may exist, and gives the FULL IUPAC protospacer (every variant column
  as its ref+alt ambiguity code, also added to the `.bed`) so a user can see all
  alleles in the window and dig into the other possible alignments.
- **Deterministic off-target representative.** When two alignments of a cluster
  tied exactly on CFD and mismatches+bulges, the reported "highest-CFD
  representative" fell to enumeration order, so the dict and dict-less paths (which
  enumerate a locus's alleles in different orders) could show a different
  representative *coordinate* for a small fraction of clusters (identical cluster,
  carriers and CFD). The greedy dense-region representative and the cluster
  representative selector now break exact ties by an intrinsic key (smaller alt
  allele; aligned-DNA + genomic position), so the representative is
  order-independent — a no-op on every already-unique cluster.
- **Dict-less correctly omits dict phantom off-targets (verified).** A genome-wide
  dict-vs-dict-less comparison leaves a handful of variant off-targets that the
  legacy dict path reports and the dict-less path does not. These are NOT missed
  regions: they are the dict's cross-variant *chimeras* — an off-target attributed
  to carriers who, in their actual (phased) haplotype, also carry an additional
  variant in the window that breaks the alignment, so no real individual carries the
  off-target. The dict-less observed enumerator (which reconstructs each individual's
  real haplotype and every sub-combination) correctly emits every off-target a real
  haplotype creates and drops these phantoms — it is *more* accurate than the dict
  here, not less (see §3/§4 of `METHODS.md`, and #118/#120).
- Report polish: combined-panel MAF is spelled out in the legend, the Section-3
  population bar no longer counts the on-target, duplicate reference rows are
  collapsed, and the section-numbering docstring matches the emitted sections.
- Doc/CLI accuracy: the Docker quick-start + web-interface guide use the v2.4.0
  image and the `-dictless` index name; `--max-total-edits` help states its real
  default (4); the web "Load Example" picks the installed default PAM; and a note
  documents pointing `HF_HOME` at a large disk on small-home systems.

## [2.3.3] - 2026-08-18

### Fixed
- **Bulge dropdowns could not be raised above 0 on a downloaded dict-less index.**
  The search form caps the selectable bulges at `min(index_max_bulges(reference),
  index_max_bulges(variant))`. A dict-less variant tarball ships the variant index
  but not the reference index, so the reference term was read as `0` and forced the
  cap — and therefore the bulge dropdowns — to `0`. The reference term is now
  **buildable-aware** (`reference_bulge_capacity`): the reference index is built
  on demand from the downloaded `Genomes/<ref>`, so its capacity is taken from the
  resolved index's `<PAM>_<N>` prefix instead of an absent folder. The default
  SpCas9 install now correctly offers up to 2 DNA / 2 RNA bulges.
- **Default search parameters restored to MM 6 / bDNA 2 / bRNA 2 with max-total-edits 4.**
  These are the intended out-of-the-box SpCas9 defaults (mismatches 6, both bulges 2,
  capped at 4 total edits); a prior change had drifted them.
- **A 0-bulge variant search failed on a dict-less / batteries install.** With
  `bDNA=0, bRNA=0` the indel branch ran a non-indexed CRISPRitz search of the
  **raw** fake-indel genome (`Genomes/<ref>+<vcf>_INDELS/`), which a dict-less
  install does not ship (build-time only) — producing missing `fakechr*.targets.txt`
  files and a `tail/head/rm/sed` cascade in `post_analisi_indel.sh`. The 0-bulge
  indel search now routes through the **indexed** indel genome
  (`genome_library/<idx>_INDELS/`, which the tarball does ship), matching the
  ≥1-bulge path. No behavior change for dict-based installs.
- **Population-distribution post-analysis crashed on out-of-range totals**
  (`max() arg is an empty sequence`). `populations_distribution.py` /
  `process_summaries.py` now guard the empty case (`max(..., default=0)`, an empty
  "no targets" figure) and the driver loop bound is computed from `mm + bDNA + bRNA`
  (the true edit budget) instead of the index `N`, so no row index falls outside the
  populated range.

## [2.3.2] - 2026-08-18

### Fixed
- **Self-built dict-less indexes are now self-complete for samplesIDs.**
  `build-index-only` (with `--samplesID`) writes the per-dataset and the combined
  `<vcf>.samplesID.txt` into the install (reusing the one union implementation), and
  `publish-index` bundles them into the index tarball — so a fresh `download
  --index-name <X>` + search works with no separate `--what samples`/`--what all`
  dependency (the download-side synthesis added in 2.3.1 becomes a pure fallback).
  No-op for single-dataset, reference, and dict-based indexes.

### Documentation
- **Rewrote the index-build docs for the dict-less flow.** `docs/PRECOMPUTED_INDEXES.md`
  now covers `build-index-only --vcf --samplesID` (**required** to emit the Tier-0
  registry + Tier-1 genotype tiers — without it you get a dicts-only index),
  `publish-index --dictless`, the separate `genotypes_<vcf>.tar.gz`, `download
  --no-genotypes`, the `-dictless` canonical install-name behavior, and the updated
  `manifest.json` fields (with a caveat that name-parsed fields are reference-index
  only). Added the variant-index build path to README §2.2.9, `crisprme_data_setup`
  §3.5, and `SETTINGS_DATA_MANAGER.md`; fixed an NGG→NRG example nit.

## [2.3.1] - 2026-08-18

### Fixed
- **Dict-less `download --what index` is now self-sufficient.** A standalone
  `download --what index` for a merged variant index (e.g.
  `NRG_3_hg38-dictless+hg38_1000G_HGDP`) now fetches the per-dataset `samplesID`
  lists it needs before synthesizing the combined `<vcf>.samplesID.txt`, so it no
  longer requires a prior `download --what all`/`--what samples`. No-op when those
  files are already present (the `--what all` → `--what index` path is unchanged).
- **Restored the `unit-tests` CI gate.** A colon in a workflow step name
  (`(#172: fakechrom match)`) made GitHub Actions reject `unit-tests.yml` at startup
  (a 0 s "workflow file issue"), so the unit-tests job had silently not run since the
  v2.3.0 merge (the tests themselves pass locally). Quoting the name fixes the parse;
  the full 22-step unit suite runs in CI again.

## [2.3.0] - 2026-08-18

The dict-less variant-analysis engine: replaces the ~152 GB per-sample SNP
dictionaries with compact, memory-mapped tiers, adds population-level off-target
summaries, and corrects the allele-frequency column. Backward-compatible — a
dict-based install still works unchanged.

### Added
- **Compact dict-less variant post-analysis.** A tiny, always-shipped **Tier-0
  registry** (per-(pos,alt) AC/AN + per-(db×subpop)/global allele counts + rsID,
  mmap+bisect) powers off-target detection and corrected allele frequencies out of
  the box; an optional **Tier-1 genotype store** reconstructs the exact per-sample
  Samples column. Together they replace the 152 GB per-sample dicts with ~7 GB
  (registry) + ~22 GB (genotype tier), random-access — SNP post-analysis no longer
  streams 152 GB of per-sample JSON. `build-index-only` emits the tiers; the search
  auto-detects and uses them (falling back to dicts when present).
- **Combination-aware population-summary output** — a `<output>.population_summary.tsv`
  companion per variant off-target with per-database, per-superpopulation and global
  allele/carrier frequencies, max-subpopulation AF (+label), homozygote counts and
  absolute Ns, with **dataset provenance preserved** (never conflates 1000G vs HGDP).
  Phased datasets use exact cis co-occurrence; unphased report an assume-cis upper
  bound + a labeled lower bound.
- **Dict-less HuggingFace distribution.** `publish-index --dictless` ships the
  registry in the main index tarball + the genotype tier as a separate optional
  `genotypes_<vcf>.tar.gz`; `download` fetches the genotype tier by default
  (`--no-genotypes` to skip). Published `NRG_3_hg38-dictless+hg38_1000G_HGDP` on
  `lucapinello/crisprme-data` (the dict-based `NRG_3`/`NGG_3` indexes remain the
  default). Fetch it with `download --what all` first, then `download --what index`,
  so the per-dataset `samplesID` lists are on disk when the combined
  `<vcf>.samplesID.txt` is auto-synthesized (a standalone `--what index` without a
  prior `samples` fetch leaves the combined file absent — hardening tracked as a
  follow-up).

### Fixed
- **Corrected allele frequencies.** The AF column — empty/mis-polarized for ~95 % of
  variants in the dict format (a documented 2.2.0 limitation) — is re-derived from
  AC/AN over the full panel.
- **Out-of-the-box dict-less install.** The skip-enrichment gates
  (`submit_job_automated_new_multiple_vcfs.sh`, `validate_inputs.py`) now accept a
  `registry_<vcf>/` tier in place of `dictionaries_<vcf>/`; `download` installs a
  `-dictless`-named index under its search-resolvable canonical name and synthesizes
  the combined `<vcf>.samplesID.txt` from the per-db files. Verified end-to-end: a
  fresh dict-less download finds the CPS1 off-target (chr2:210530658, CFD 0.947,
  rs114518452) with Samples reconstructed from the genotype tier and a populated
  population summary.

## [2.2.0] - 2026-08-18

First stable release of the CRISPRme+ 2.2.0 line (Python 3.11 / Dash 2.x). It
consolidates the `2.2.0-alpha.27`–`alpha.30` pre-releases (see those sections
below for the full history); the entries here are the changes since `alpha.30`.

### Fixed
- **Phased multi-SNP haplotype off-targets are no longer under-reported.** In the
  phased post-analysis path — the default for phased datasets such as 1000G/HGDP —
  `iupac_decomposition` failed to assemble the full-haplotype off-target when one
  sample carried **≥4 co-occurring variants in a single protospacer window**: the
  true worst target (fully substituted, lowest-mismatch/highest-CFD) was silently
  under-reported or dropped. The in-loop level-0 subtraction (a dedup device that
  starved the deeper combination layers) is replaced with a **deferred peel** that
  runs after the full lattice is built, so the maximal *cis* combination forms while
  per-haplotype attribution stays deduplicated. Validated byte-identical for every
  previously-correct case (unphased, phased ≤3-variant, single-SNP, >cap greedy)
  plus a real-data genome-wide no-regression diff. (#41; a long-standing defect
  inherited from classic CRISPRme, tracked there at pinellolab/CRISPRme#175.)

### Changed
- `environment.yml` now includes `ijson`, so from-source installs get the streaming
  low-RAM SNP-dict reader that the Docker image already had (without it, SNP
  post-analysis falls back to a whole-file `json.load`; results are identical, but
  RAM use and the OOM-guard estimate are higher).

## [2.2.0-alpha.30] - 2026-08-17

### Fixed
- **Indel off-targets were silently dropped since v2.1.9.** `post_analisi_indel.sh`
  subset the per-chromosome indel targets with `grep -F -w $chrom`, but the indel
  search runs on a per-chromosome *fake* genome (`pool_search_indels.py` searches and
  names its targets `fake<chrom>`), so every indel target row's Chromosome column is
  `fake<chrom>` (e.g. `fakechr22`) — which `grep -F -w chr22` can never match (`-w`'s
  left word boundary fails because `chr22` is preceded by the word character `e`). The
  per-chromosome subsets came back empty, the indel post-analysis processed nothing,
  and **all indel off-targets were dropped with a clean exit**. Both grep lines now
  match on `"$fakechrom"` (`-w` still prevents `fakechr2` matching `fakechr22`); the
  `NF >= 10` malformed-line guard that the same refactor had dropped is restored; and
  an empty subset now prints a WARNING to stdout (stderr is fatal in this pipeline) so
  the failure can never be silent again. A hermetic regression test
  (`test_indel_chrom_subset.py`) is added and wired into CI. SNP post-analysis is
  unaffected (SNP targets use real chromosome names). Thanks to @munchr-gene1 for the
  report and diagnosis (#172).

## [2.2.0-alpha.29] - 2026-08-16

### Fixed
- **Variant search without `--annotation` no longer crashes on a read-only install
  (apptainer SIF).** `post_process.sh` wrote the annotation intermediate BEDs
  (`<annot>.tmp.bed` / `.tmp.sorted.bed`) next to the annotation file; with no
  `--annotation`/`--gene_annotation` that file is the shipped empty placeholder
  `PostProcess/vuoto.txt` inside the (read-only) install dir, so the writes failed
  and the run died at "Integrating results". Temps now go to the run's writable
  output dir. As defense-in-depth, `resultIntegrator.py` guards the closest-gene
  distance parse so a missing/empty annotation field is treated as "no annotation"
  (columns stay `NA`) instead of raising `ValueError: could not convert string to
  float: ''`. The Docker/web path (real annotation, writable FS) is byte-identical.

### Added
- **Pre-flight check for unpadded guides.** A CLI `--guide` file whose guides lack
  the PAM's N-padding (e.g. the bare 20 nt `CTAACAGTTGCTTTTATCAC` instead of
  `CTAACAGTTGCTTTTATCACNNN`) previously ran the full multi-hour search then crashed
  deep in `radar_chart_dict_generator.py` (`IndexError`). The lightweight validator
  now catches this up front with an exact fix message (append/prepend N's), handling
  both 3' (SpCas9) and 5' (Cas12a) PAM orientations. The web is unaffected (it
  auto-pads guides). Covered by 10 new `check_guide_pam_consistency` tests; the whole
  `test_validate_inputs` suite is now run in the unit-tests CI (previously it wasn't).

## [2.2.0-alpha.28] - 2026-08-16

### Fixed
- **Batteries variant search no longer aborts at pre-flight validation.** The
  lightweight input validator (`validate_inputs.run_lightweight`, always-on before
  `complete-search`) hard-required `VCFs/<dataset>/` to exist and aborted otherwise —
  but a batteries-included install (`download --what index`) ships the precomputed
  variant index + per-sample dictionaries *without* the multi-GB source VCFs. The
  validator now recognizes that case: when the VCF dataset directory is absent but the
  dictionaries (`Dictionaries/dictionaries_<dataset>/` + `log_indels_<dataset>/`) are
  present, it passes the check (source VCFs genuinely not required), mirroring the
  skip-enrichment logic in `submit_job_automated_new_multiple_vcfs.sh`. A genuinely
  missing dataset (no dictionaries either) still fails loudly. Completes the alpha.27
  batteries fix: `download --index-name NRG_3_hg38+hg38_1000G_HGDP` + a variant search
  now runs end-to-end on a fresh install. Covered by two new `run_lightweight` tests.

## [2.2.0-alpha.27] - 2026-08-16

### Fixed
- **Batteries variant search is now self-sufficient (the index-only download can run
  a variant search).** The variant post-analysis needs per-sample SNP dictionaries
  (~152 GB for 1000G+HGDP) that weren't shipped, so a fresh `download --index-name
  <variant>` install failed at genome enrichment (no source VCFs). The dicts are now
  shipped **gzipped** (152 GB → ~28 GB; indel logs 18 GB → 3.3 GB) and read **on the
  fly**, mirroring the existing bgzipped-VCF pattern — no 150 GB on-disk decompress.
  - `new_simple_analysis.py` / `analisi_indels_NNN.py`: prefer `*.json.gz` / `*.txt.gz`
    (ijson streams through gzip; json.load + line reads via `gzip.open`), fall back to
    plain — backward-compatible. Bash wrappers unchanged (Python resolves the extension).
  - `pool_search_indels.py` / `pool_post_analisi_indel.py` + the `submit_job` batteries
    chrom glob: strip an optional `.gz` from `log*.txt(.gz)`.
  - `build-index-only` gzips the per-chromosome dicts it produces (pigz-preferred).
  - Verified end-to-end by `validate-benchmarks` (build gzips → search reads `.gz`) and
    a new `test_dict_gzip` unit test (publish bundles + download routes `.gz` dicts).

## [2.2.0-alpha.26] - 2026-08-16

### Changed
- **Default SpCas9 PAM is now NRG (NAG + NGG), not NGG.** The web form and the
  "Load Example" default, plus the CLI/download docs, now default to `NRG`
  (`20bp-NRG-SpCas9`) and the `NRG_3_*` indexes. NRG matches SpCas9's broad
  recognition, so **variant-created NAG off-targets are found out of the box** —
  e.g. the CPS1 off-target from the CRISPRme paper (rs114518452: chr2:210530658,
  2 mismatches + 1 RNA bulge, variant-created `TGG` PAM, CFD 0.947). Verified
  end-to-end on a from-scratch `NRG_3_hg38+hg38_1000G_HGDP` index. `_default_pam`
  prefers NRG, falls back to NGG, then the first available PAM. The `NGG_3_*`
  indexes remain available for canonical-NGG-only searches. New `NRG_3_*` indexes
  are published on `lucapinello/crisprme-data`.

### Fixed
- **Off-target search now uses all requested cores.** `submit_job` split
  `$ncpus` in half and ran reference + variant searches in parallel; the small
  reference genome finishes first, so the long variant-search tail ran on only
  half the cores. Now runs the two searches sequentially at full `$ncpus`
  (reference 44.6 → 25.6 min in testing; `searchTST` 24 → ~45 cores).
- **Freshly-built / downloaded indexes now search correctly.** `build-index-only`
  writes `.pam_build`/`.display_label` sidecars inside the index dir, but
  CRISPRitz's `-index` search requires the dir to contain only `.bin` files
  ("only .bin files" error). Added a symlink-only `.bin` view so the sidecars no
  longer break the search, regardless of how the index was built or downloaded.

## [2.2.0-alpha.25] - 2026-08-15

### Changed
- **Contacts page tidied.** Developers are now listed with the Pinello lab first
  (Luca Pinello), then the Verona/InfOmics group (Manuel Tognon et al.), then the
  Bauer lab (Linda Lin, Daniel Bauer). Comments are directed to Luca Pinello only
  (single email), and the page now says bugs should be filed as a GitHub issue.
  Also fixed a missing space after the comma in the per-lab author lists.

## [2.2.0-alpha.24] - 2026-08-15

### Changed
- **Submission page: consolidated to exactly three font sizes (Large / Medium /
  Small).** The form previously mixed ~7 inline sizes on a 10px root, which left
  helper notes rendering *larger* than the inputs they described. It now uses one
  clear hierarchy: Large `2.2rem` for the numbered section headers ("Select
  gRNA/genome/PAM/thresholds"), Medium `1.5rem` for every interactive control
  (spacer textarea, dropdowns, Load Example / Submit buttons, the "Maximum edits"
  label and form body), and Small `1.25rem` for helper text, notes and slider
  tick marks. Verified via full-page + high-DPI screenshots — no clipping or
  overflow. Cosmetic only; no change to search behavior.
- **"Max edits" slider maximum is now derived, not a magic number.** The slider
  topped out at a hardcoded 5; it now tops out at the app's mismatch ceiling
  (`MAX_MMS - 1` = 6). Rationale: mismatches are the non-index-limited edit
  dimension and the dominant search cost, so the total-edits knob should reach as
  far as the mismatch dropdown does. The installed index caps *bulges* only
  (`bMax = N - 1`, ≤2 in practice) — that limit is applied separately where
  simple-mode bulges are derived and must not shrink the mismatch-driven total.
  Default stays 4.

## [2.2.0-alpha.23] - 2026-08-15

### Changed
- **"Max edits" is now a self-sufficient single knob — bulges are derived from it.**
  Previously a variant off-target that needs a bulge (e.g. the CRISPRme paper's CPS1
  off-target rs114518452: 2 mismatches + 1 RNA bulge) was only found if bulges were
  explicitly enabled — the CLI defaulted to bulge-free, so `--max-total-edits` alone
  found no bulge targets. Now:
  - **CLI:** when `--bDNA/--bRNA` are omitted, they are derived from `--max-total-edits`,
    capped by the bulge depth the installed index supports (so `complete-search --vcf …
    --max-total-edits 4` finds up to 2mm+2bulge patterns with no bulge flags, reusing the
    installed index — never forcing a heavier build; explicit `--bDNA/--bRNA` still win).
  - **Web:** in simple mode the per-type bulge caps now track `min(max_edits, index_cap)`
    instead of a fixed 2/2, so the slider is genuinely the only control needed.
- **Consistent default cap of 4** across web slider, Load Example, CLI `--max-total-edits`,
  and the shell fallback (was 5). Four surfaces the CPS1 off-target (needs 3 total edits)
  with headroom and is the validated performance sweet spot (cap 6 was pathologically slow
  on genome-wide combined searches).
- Docs: `--bDNA/--bRNA` are documented as optional (derived from `--max-total-edits`), not
  "Required"; the variant-aware requirement for variant-created off-targets is noted.

## [2.2.0-alpha.22] - 2026-08-15

### Added
- **History page now reports the Max edits used.** The results-history table only showed
  Mismatches / DNA bulge / RNA bulge; it now has a "Max edits" column showing the total
  mismatches+bulges cap (the simple-mode slider value, or the per-type sum with an
  "(advanced)" note). Jobs predating this field show "-".

## [2.2.0-alpha.21] - 2026-08-15

### Fixed
- **Graphical Reports: every variant off-target now shows its blue "Alternative" point.**
  The top-1000 scatter drew a blue point only when a variant had a known allele frequency;
  a variant whose MAF was missing (e.g. the frequency is not annotated from the external
  AF source) was treated like a reference-only site and its blue point was suppressed. On
  the combined 1000G+HGDP panel that hid ~93% of the variant points, so the plot looked
  "all red." Now a blue point is drawn for any off-target carried by >= 1 sample, sized by
  allele frequency when known and a default size when not; only true reference-only sites
  get no blue point. This reproduces the CRISPRme Nature Genetics figure — e.g. the CPS1
  off-target rs114518452 (chr2:210530658, CFD 0.021 -> 0.947, a variant-created off-target)
  appears as the top blue point. Also fixed in the per-sample Personal Risk Card plots.

## [2.2.0-alpha.20] - 2026-08-15

### Added
- **Graphical Reports: a "top 1000 by variant effect" plot.** Alongside the existing
  score-ranked top-1000 scatter, the CFD and ML-scorer reports now also show a companion
  plot ranked by the size of the variant-induced score change (`|ALT - REF|`), so the
  variants that actually change the off-target score are foregrounded instead of being
  buried among the many that leave it unchanged. The score-ranked plot is unchanged; the
  fewest-mm+b report stays score-only. (`CRISPRme_plots.py`, shown in `results_page.py`.)

## [2.2.0-alpha.19] - 2026-08-15

### Fixed
- **Graphical Reports: the top-1000 scatter now shows the variant (ALT) points.** The blue
  ALT markers were drawn ~10x smaller than the red REF markers at the same x, so they were
  swallowed and the plot looked "all red." They now have a floored size + a black edge (and
  stay hidden for reference-only sites), so a variant that changes the CFD / ML-scorer score is
  visible at the tip of its arrow. (`CRISPRme_plots.py`, `CRISPRme_plots_personal.py`.)
- **Graphical Reports: the ML-scorer and fewest-mm+b top-1000 plots are generated again.** The
  plot generator sorted the ML-scorer plot on a non-existent column
  (`CFD_score_(highest_CRISPR_BULGE)`) and crashed with a silently-swallowed `KeyError` right after
  the CFD image, so only the CFD plot appeared. Fixed to `CRISPR_BULGE_score_(highest_CRISPR_BULGE)`.
- **Reference-only searches no longer break the Graphical Reports tab.** The population-barplot
  callback `Output` targeted a container the layout omitted for reference-only searches, which
  could fail the whole callback (radar chart included). The container is now always present
  (hidden for reference-only).
- **Personal Risk Card now shows the candidate sites.** The detail table was built only from
  *private* targets (unique to the sample), so it was empty for the common case of a sample
  with shared-but-not-private candidates. It now shows both a **Personal candidates** table
  (all off-targets the sample carries) and the **Private** table.
- `update_images_tabs` now guards an empty table selection (no crash on initial render) and
  applies the non-SpCas9 -> "fewest" criterion remap the rest of the page uses; fixed a radar
  chart `except` handler that assigned an unused variable.

## [2.2.0-alpha.18] - 2026-08-15

### Fixed
- **INDEL post-analysis: forward-ported the boundary/dense-variant crash guards.** The
  alpha.17 ML-scorer window guard was only in the SNP path (`new_simple_analysis.py`); the
  parallel INDEL path (`analisi_indels_NNN.py`) still assumed a full 29 nt window, so an
  out-of-range indel coordinate could crash post-analysis. The `!= 29` null-guard is now
  mirrored in both INDEL ML-scorer blocks.
- **INDEL CFD scoring no longer KeyErrors on an unexpected PAM.** `analisi_indels_NNN.py`
  used a raw `pam_scores[pam]` where the SNP path uses the guarded `pam_scores.get(pam,
  0.0)`; now matched (issue-#94 class).
- **CLI `--sorting-criteria` default now matches its own `--help` and the web.** It
  defaulted to `mm+bulges` but documented (and the web uses) `mm+bulges,mm`; the
  fewest/total tie-break now ranks tied targets identically on both paths.

### Changed
- **CLI and web now share the same default total-edits cap (5).** The web "Maximum edits"
  slider defaulted to 3 while the CLI `--max-total-edits` defaults to 5; the web slider
  (and Load Example) now default to 5 so a default web search and a default CLI search
  cover the same edit space.

## [2.2.0-alpha.17] - 2026-08-14

### Fixed
- **Post-analysis no longer crashes on variant-dense / chromosome-boundary targets.**
  On the combined 1000G+HGDP index the example search (and any dense-variant search)
  could abort in post-analysis with a `ZeroDivisionError` (ML scoring divided by the
  length of a 29 nt window that had been stripped empty of IUPAC codes) or an
  `IndexError` in `iupac_decomposition` (a reference window truncated at a chromosome
  boundary was shorter than the target). Both came from assuming a genome slice is always
  full length. Now an un-scoreable ML-scorer window is nulled (score −1, as for windows with
  `N`) and IUPAC positions past a truncated reference are skipped, with belt-and-suspenders
  divide-by-zero guards in the ML scorer module. Verified end-to-end on the batteries
  1000G+HGDP data: the example completes and its Personal Risk Card renders.

### Changed
- **Base editing is off by default after "Load Example."** The example now loads a
  standard (non-base-editing) search; the base-editor toggle stays on "No".
- **Larger helper text on the submission form.** The small gray guidance texts in the
  submission cards were ~10.5 px (the page root font is 10 px); enlarged to ~14 px so
  they are readable. The page's base font is unchanged.

## [2.2.0-alpha.16] - 2026-08-14

### Fixed
- **Result title is now mode-aware and no longer self-contradictory.** In simple mode
  the title listed the per-type mismatch/DNA/RNA caps *and* the total "Max edits" cap,
  which contradict when the total cap is tighter (e.g. `Mismatches 4 - DNA bulges 1 -
  RNA bulges 1 - Max edits 1`). Now a simple-mode search shows only `Max edits
  (mismatches + bulges) N` (the single control the user set), while an advanced search
  shows the explicit per-type caps. `change_url` records `Threshold_mode` (simple/
  advanced) in `.Params.txt`; older jobs are inferred (a total cap tighter than the
  per-type sum ⇒ simple) and pre-alpha.15 jobs fall back to the per-type title.
- **Dedup reuse wrote a malformed completion marker.** When a resubmission reused a
  prior clean result, the marker was written `Job\nDone` (newline) instead of
  `Job\tDone` (tab), so the dedup check and the status page's completion detection
  would not recognize the reused job as finished. Fixed to write `Job\tDone`.

## [2.2.0-alpha.15] - 2026-08-14

### Fixed
- **The result title now shows the "Max edits" cap that actually bound the search.**
  In simple mode the "Maximum edits (mismatches + bulges)" slider is the binding
  constraint, but the Result Summary title only listed the per-type mismatch/DNA/RNA
  caps — which can be looser than the total cap (and, after "Load Example", are values
  the user never explicitly chose), so the title overstated how broad the search was.
  The total-edits cap is now recorded in `.Params.txt` (`Max_total_edits`) and appended
  to the title as "Max edits (mismatches + bulges) N". Backward-compatible: results for
  jobs created before this field simply omit the clause.

### Changed
- **"Load Example" now sets the Max-edits slider (to 4) explicitly** instead of
  leaving it at the default, so the example is driven by the same total-edits cap the
  report displays and the per-type caps it fills in (4/1/1) are consistent with it.
- Docs: added mentions of the optional **email notification** (README, Docker
  quickstart, landing pages) and clarified that **Personal Risk Cards** appear only
  for variant-aware searches (landing pages, in-app help, quickstart, user guide).

## [2.2.0-alpha.14] - 2026-08-14

### Fixed
- **A resubmitted search no longer resurfaces a previously-FAILED job.** Results are
  deduplicated by search parameters; the dedup reused any matching prior job, so
  re-running a search whose earlier attempt had failed (e.g. left behind by a bug
  since fixed) returned that stale failure — surfaced in the UI as "The selected
  result encountered some errors, please remove it and try to submit again." Dedup
  now skips prior results that did not finish cleanly (non-empty `log_error.txt` or
  no `Job Done`) and runs the search fresh instead.

### Changed
- Docs: dropped the stale "peak memory is in post-processing and can exceed
  32 GB / ~64–100 GB" guidance (the alpha.12 streamed SNP-dict load removed that
  post-analysis spike); the 16 / 32 / 64 GB tiers stay.

## [2.2.0-alpha.13] - 2026-08-14

### Changed
- **Post-analysis performance** (follow-up to the streaming dict-load fix):
  - Added the `yajl` C library + `cffi` so `ijson` uses its fast `yajl2_cffi` C
    backend instead of the pure-python fallback (much faster streaming of the SNP
    dictionary).
  - The post-analysis worker cap no longer sizes per-worker RAM from the (large) SNP
    dictionary file: with streaming, per-worker peak is small, so it uses a modest
    fixed estimate and runs chromosomes in parallel again (the SNP path checks that
    `ijson` is importable; the INDEL path never loaded the big dict). Previously it
    over-estimated and serialized to a single worker.

## [2.2.0-alpha.12] - 2026-08-14

### Fixed
- **On-demand INDELS index built with a doubled `<PAM>_<N>` prefix** (regression
  from the alpha.10 indels-index-resolution change). crispritz `index-genome`
  already auto-prepends `<PAM>_<bMax>_` to the output dir, but the builder also
  added it manually, producing `NGG_2_NGG_2_..._INDELS`; detection/search look for
  the single-prefix name, so the indel search produced no output and the run failed
  with "off-targets search failed". Removed the manual prefix (let crispritz add
  it). This unblocks `validate-benchmarks` (complete-test + validate-test now match
  the brute-force ground truth: 2/2 light benchmarks). Batteries searches were
  unaffected (they use a correctly-named pre-built `_INDELS`).
- **Genome-wide variant post-analysis OOM ("Killed … EmptyDataError").** On a
  genome-wide 1000G+HGDP search, each per-chromosome worker `json.load()`ed the
  whole chromosome SNP dictionary — a chr2 dict is ~13 GB on disk and ~26 GB in
  RAM — so even one worker could exhaust a normal machine (and the concurrency
  guard defaulted to a fixed 64 GB budget / 4 GB-per-worker, over-subscribing on
  smaller hosts). `new_simple_analysis.py` now loads **only the dictionary entries
  the search's off-targets actually query** (streamed with `ijson`), dropping peak
  RAM from ~26 GB to <0.1 GB on the real chr2 dict (byte-identical results). The
  worker-cap guards now **detect the machine's real RAM** (host + cgroup limit),
  **estimate per-worker RAM from the actual dictionary sizes**, and print a clear
  "needs ~N GB / may be OOM-killed" warning instead of the cryptic cascade. Added
  `ijson` to the image.

### Changed
- Docs: corrected the in-app manual and web guide (removed the stale "up to 100
  spacers / first-100-sequences" copy; annotation is Settings-managed not a Step-3
  dropdown; email is configured in Settings, not "the server notifies").

## [2.2.0-alpha.11] - 2026-08-14

### Fixed
- **Graphical Reports population barplot never rendered.** The pipeline writes
  `populations_distribution_<guide>_<N>total_<score>_new.png`, but the results
  page read the name without the `_new` suffix, so every mismatch/bulge combo
  showed "No result found for this combination". Added the missing suffix.
- **Graphical Reports defaulted to an empty combo.** The "up to N edits" selector
  defaulted to `0` (on-target only), which is empty for typical sparse guides. It
  now defaults to the full budget (mm + bulges) — the cumulative barplot the
  pipeline only emits when it has data — so a plot shows immediately.
- **"Load Example" left the base-editing window blank.** The `be-window-start/stop`
  dropdowns start with no options, so the example's window values were dropped
  before `update_base_editing_dropdown` populated the options. Load Example now
  emits the window options atomically with the values (and passes ints, matching
  the option type), so it fills 4/8 correctly.

## [2.2.0-alpha.10] - 2026-08-14

### Added
- **Email notifications are now configurable in the app** (Settings → Email
  notifications, local mode): SMTP host/port/SSL, sender, and an app password are
  saved to a local `.email.json` (0600, password never echoed back to the
  browser). `send_mail.py` reads that config (falling back to the
  `CRISPRME_SMTP_*` env vars) and is now best-effort — it never writes to stderr
  or raises, so a mail-server problem can't fail an already-completed search.

### Changed
- **Search form polish**: the CRISPRme+ logo now sits above the first step; the
  *Load Example* button moved to the bottom of the *Select gRNA* box; larger,
  more legible fonts throughout the form; the spacer help text explains that
  multiple spacers can be entered (and points to the CLI for bulk runs); the
  *Notify me by email* box explains the one-time Settings setup.

### Fixed
- **Variant indel indexing/search now resolves the INDELS companion from the
  actually-resolved variant index** instead of reconstructing an exact
  `<true_pam>_<bMax>_..._INDELS` name. A batteries/precomputed index whose N or
  PAM differs from the request (e.g. the shipped `NGG_3` served against a
  `bMax=2` search, or a pamless `NNN` index) was not matched, so the pipeline
  attempted an impossible rebuild and crashed in `pool_index_indels.py`
  (`FileNotFoundError: .../Genomes/<ref>+<vcf>_INDELS/`). The INDELS index, its
  search path, and the on-demand builder now all key off the index's own
  `<PAM>_<N>` prefix (consistent with `build-index-only`/`publish`/`download`),
  so newly built-and-uploaded indexes work generally. When the companion is
  genuinely absent and cannot be rebuilt (no source genome), the run now fails
  with a clear "re-download the index" message instead of a raw traceback.

### Added
- **Annotation manager** in Settings → Data Manager: a per-genome
  *Manage annotations (enable / disable)* checklist backed by a persisted
  `.enabled.json` manifest, plus format validation on every uploaded annotation
  BED (≥4 tab columns, integer `start ≤ end`, whitespace-free label, size cap;
  malformed files are rejected with a clear message and nothing is written).
  Enabled tracks are merged on demand into a single cached active annotation and
  applied automatically to every search.

### Changed
- Annotations are now applied automatically from the enabled set rather than
  chosen per search. The annotation dropdown was removed from the search form
  (Step 3 is now just email + job name); the built-in
  ENCODE cCREs (SCREEN) + DHS + GENCODE bundle is enabled by default so
  annotations work out of the box.
- The **Maximum edits** slider is floored at 1 (the on-target / 0-edit hit is
  always reported at any setting), avoiding an empty-result edge case on large
  precomputed indexes.
- Search-submission form: consistent left indentation for the PAM and
  variant-dataset dropdowns and the submit/example buttons.

## [2.2.0] - 2026-08-06

### Added
- Graphical **Settings / Data Manager** page in the web interface: add reference
  genomes (UCSC by assembly name — e.g. the pig `susScr11` — HuggingFace, or a
  direct URL), precomputed indexes (download from HuggingFace or build locally
  from an installed genome + PAM), VCF datasets (HuggingFace or register an
  existing server folder), annotations (BED upload), and nucleases/PAMs (a small
  form). New data lands in the local data folder and is auto-discovered by the
  search form. Long operations run as detached jobs on a dedicated executor with
  live progress, so they never starve the search slots. Mutations are local-mode
  only; publishing an index to the shared HuggingFace repo is maintainer-only.
  `download --what genome` gained `--source {hf,ucsc,url}` (+ `--url`) so the CLI
  and web share one non-human-genome download path.
- Python 3.11 modernization: pipeline fixes for pandas 2.x / matplotlib 3.x and
  a Dash 1.x → 2.x web-app migration, plus a Python-3.11 Docker image built from
  source (CRISPRitz 2.8.1) (#131).
- `assembly-search` subcommand: off-target search on a personal diploid genome
  assembly (two haplotypes, no VCF), reconciled to hg38 via liftOver (#113).
- Reference-index UX: `build-index-only` pre-builds the reusable CRISPRitz
  reference index without running a search, and `complete-search --index-path`
  reuses a prebuilt/staged index library (a missing index is a hard error rather
  than a silent rebuild).
- HuggingFace data distribution: `download` fetches reference data (genome,
  annotations, PAMs, sample IDs, VCFs, precomputed indexes) from a HuggingFace
  dataset repository over its CDN, and `publish-index` uploads a locally built
  index for reuse. Default repo `lucapinello/crisprme-data`, overridable via
  `--hf-repo` / `CRISPRME_HF_REPO`. `setup`/`complete-test` also try HuggingFace
  first and fall back transparently to the original UCSC/EBI/Sanger sources
  (#140, #141).
- `complete-search --max-total-edits N`: cap the total edits (mismatches +
  bulges) per reported alignment; over-cap targets are dropped right after the
  search, shrinking intermediate files and post-analysis time (#107).

### CI
- New `unit tests` workflow: fast, hermetic byte-compile + network-free HF/index
  unit tests on every code PR.
- New `web e2e (playwright)` workflow: builds the py3.11 image, serves the web
  app, and drives Chromium to assert every Dash 2.x page renders (no blank pages
  / JS errors).
- `validate-benchmarks` gained a new-subcommand dispatch + unit-test smoke step.

### Changed
- Clearer failure reporting: when a search fails, CRISPRme now prints *which
  stage* failed (from the per-stage log) and the last lines of the error log,
  instead of only "run failed — see log_error.txt". Makes failures actionable
  for non-expert users.

### Fixed
- Web interface (Dash 2.x) hardening, from a full Playwright stress test of the
  running app:
  - The web server no longer crashes on a from-source install. Dash 2.x's
    `app.run()` lets the `HOST` environment variable override the host argument,
    and the from-source conda env sets `HOST` to a non-bindable compiler build
    triple (`x86_64-conda-linux-gnu`); the server now forces the intended host/
    port so it binds correctly.
  - The **Query Genomic Region** and **Personal Risk Cards** result tabs no
    longer return HTTP 500. Both callbacks type-checked their inputs before the
    "no click yet" guard, so Dash's initial (empty) render raised a `TypeError`;
    the guard now runs first, and Filter/Generate with nothing selected is a
    graceful no-op.
  - Removed the dead cross-origin "skeleton" stylesheet (blocked by browsers on
    every page); the layout already uses the Bootstrap grid.
  - The nuclease dropdown collapses case-variant duplicate PAM files so each
    nuclease is listed once.
- Zero-hit searches now complete cleanly with an empty result instead of
  aborting. A search that finds no off-targets (e.g. a very stringent
  guide/parameter combination) previously failed part-way through post-analysis
  ("off-targets post-analysis (reference) failed", then a cascade through the
  rsID / summary / integration steps, all of which assumed at least one target).
  Added a zero-target guard to the reference SNP post-analysis (mirroring the
  INDELs one), made `remove_n_and_dots.py` tolerate a header-only report, and
  added a high-level short-circuit that emits an empty-but-valid result and exits
  0 when no off-targets are found. Verified end-to-end on ml007 (empty result
  exits 0; a normal with-hits search is unaffected).

### Documentation
- New `docs/DOCKER_QUICKSTART.md` + a README quickstart callout: a few-command,
  no-conda / no-410 GB path to the web interface for non-experts — fast HF data
  download, a prebuilt index, then `web-interface` in the browser, with
  troubleshooting and "install more indexes as needed".
- Data-setup guide: documented what `setup` produces (including the combined
  1000G+HGDP config files), the HuggingFace fast-download path, and a new
  "Prebuild, reuse, and share the reference index" section.
- README: added a from-source install path and a reference-index /
  data-distribution commands section.
- Documented the variant PAM behaviour in `docs/INPUT_FORMATS.md`: PAM *creation*
  (variants that add an off-target) is reported; PAM *disruption* (variants that
  remove a reference PAM for carriers) is a known, intentional gap — with the
  rationale (disruption lowers rather than raises predicted risk and would need
  sample-specific, non-scorable semantics).

## [2.1.13] - 2026-08-04

### Added
- Robust `complete-test` downloads: HTTP retry with linear backoff plus
  checksum-verified resume, so a dropped connection or truncated transfer from a
  slow/flaky host retries (and skips already-verified files) instead of aborting
  a multi-hour run (#136).
- Memory-bounded post-analysis: per-chromosome SNP/indel workers are capped to a
  RAM budget (default 64 GB, overridable via `CRISPRME_MAX_MEM_GB` and
  `CRISPRME_POSTPROC_WORKER_GB`), preventing the peak-RAM spike observed on
  genome-wide 1000G runs (#136).

### Fixed
- Post-analysis worker-count diagnostic is now written to stdout (`log_verbose`)
  instead of stderr. The caller treats a non-empty stderr log as a fatal
  post-analysis failure, so the informational message previously aborted every
  run with a false "post-analysis (snps) failed" even when the analysis
  succeeded (#136).

## [2.1.12] - 2026-08-01

### Added
- Cas9 and Cas12a validation benchmark examples with precomputed brute-force
  ground-truth references and an extensible registry
  (`test/benchmark/benchmarks.json`); `validate-test` compares CRISPRme
  off-targets against them (#116).
- A fast, dependency-free Rust port of the brute-force ground-truth generator
  (`test/benchmark/rust/`) that produces output identical to the Python
  generator (#116).
- `validate-benchmarks` GitHub Actions CI that runs the full
  `complete-test` → `validate-test` round-trip on Linux for every registered
  benchmark (#116, #120).
- Native Apple Silicon (`linux/arm64`) Docker image via a multi-arch `buildx`
  workflow; the multi-arch image is published to Docker Hub on release tags
  (#121).
- Low-memory startup warning (`PostProcess/memory_check.py`) advising Docker
  Desktop users to raise the memory limit when total RAM is below ~32 GB (#121).
- Search-space budget estimator and guard (`PostProcess/search_budget.py`) that
  warns before resource-explosive runs, with a companion analysis in
  `docs/SCALABILITY_ANALYSIS.md` (#118).
- Release tooling: a `release-crisprme` Claude Code skill, `docs/RELEASING.md`,
  this `CHANGELOG.md`, and `scripts/prepare_release.py` (#117); plus docs for the
  fast Rust generator and using the release skill (#122).
- Input-format hardening: `--gene_annotation` is now auto-sorted like
  `--annotation`; clear errors for multi-line PAM files and malformed PAM
  filenames; and a non-fatal warning when a degenerate PAM motif is combined with
  bulges (the CRISPRitz #105 crash, fixed in CRISPRitz 2.7.1) (#125).
- Documentation: `docs/INPUT_FORMATS.md` (PAM/Cas12a/VCF/chromosome/annotation/
  bulge guidance, #124), `docs/ISSUE_AUDIT.md` (a review of all open + closed
  issues, #126), and `docs/ROADMAP.md` (the 2.1.12 → 2.1.13 plan, #128).

### Fixed
- Multiple bugs affecting new VCF dataset processing, concurrent runs, and
  gnomAD handling (#96).
- Guard the CFD PAM-score lookup against non-ATCG PAM bases present in real hg38,
  which previously crashed the run with a `KeyError` (#94, #125).
- Guard the post-analysis intermediate reads against non-UTF-8 bytes
  (`UnicodeDecodeError`) so a stray byte no longer aborts a run (#52, #125).
- `validate-test` now exits non-zero when a benchmark mismatches or the
  `complete-test` output is missing, instead of passing silently (#116).
- `complete-test` downloads 1000 Genomes VCFs over HTTPS instead of FTP,
  unblocking FTP-restricted and CI networks (#116).

### Changed
- Synced the `Dockerfile` `crisprme_version` pin to the released version, which
  had lagged behind `crisprme.py` and the Bioconda recipe (#117).

## [2.1.11] - 2026-07-06

### Fixed
- Consistent `.gz` suffix for compressed annotation files, removing a mismatch
  between file content and file name in the standard and personal-annotation
  pipelines (`crisprme.py`).
- `_check_personal_annotation` now strips a `.gz` suffix before processing so
  downstream steps see the expected filename.
- Removed double-compression of annotation files, preventing corrupted or
  unreadable annotation outputs.

### Changed
- Refactored per-chromosome test genome-directory handling: replaced
  `ensure_hg38_directory` with `_assign_genome_directory_name` (whole-genome
  builds in `<dest>/hg38`, single-chromosome builds in `<dest>/hg38_<chrom>`).
- `download_genome_data` now returns the resolved genome directory path directly.
- Renamed downloaded sample-ID files to include genome/dataset context
  (`hg38_{dataset}.samplesID.txt`) via a centralized mapping.
- Reduced default web-interface example-run parameters (max mismatches 6 → 4,
  DNA bulges 2 → 1, RNA bulges 2 → 1) for a faster example run.

### Added
- New `_write_sg1617_file` helper to generate an `sg1617` guide fixture during
  legacy database setup.
- New guides `docs/crisprme_data_setup_051826.md` (CLI setup/usage) and
  `docs/crisprme_web_interface_user_guide.md` (web interface walkthrough).
- "Setup Legacy Database" section in `README.md`.

### Removed
- Unused `_bgzip_ann_data` helper and stale typing imports/docstrings.

## [2.1.10] - 2026-05-29

### Added
- New `setup` command for automated download/installation of CRISPRme reference
  resources directly from the command line, removing the need for manual
  downloads from the CRISPRme website.

### Changed
- Updated documentation and container configuration to support the new
  installation workflow and improve reproducibility across deployments.

## [2.1.9] - 2026-01-16

### Added
- `validate-test` functionality: validates `complete-test` off-target
  predictions against brute-force ground-truth alignments derived from
  1000 Genomes variant data (PR #92).
- Chromosome-level validation support (single chromosome or genome-wide).

### Changed
- Strengthened the `complete-test` pipeline and refined benchmarking utilities
  for correctness and reproducibility with population-scale variant data.

## [2.1.8] - 2025-12-10

### Fixed
- Targeted fixes and stability improvements across annotation handling,
  off-target search execution, test coverage, and report generation.
- Fixed and updated the `complete-test` workflow to align with the latest
  internal logic and outputs (PR #86).

### Changed
- Updated Docker distribution and testing environment for consistent runtime
  behavior (PR #85).

## [2.1.7] - 2025-03-10

### Added
- Support for exome-based gnomAD VCF files.

### Fixed
- Corrected handling of indels in off-target reporting (previously missing or
  misreported).
- Fixed the visualization of the Personal Card tab in the GUI and website.
- Closed file streams after creating graphical reports to prevent memory issues.

### Changed
- Improved contiguous-target merge logic (always reports the leftmost off-target
  in a cluster).
- Improved error-tracing logic.

## [2.1.6] - 2024-11-27

### Added
- Support for gnomAD v4.1 VCFs, including joint variant files.
- New `complete-test` function for quick single-chromosome and full-genome
  off-target testing.

### Fixed
- Corrected indel handling in off-target reporting to eliminate false positives.
- Fixed the alternative-alignments report to list all possible alignments.
- Resolved Personal Card tab visualization issues in the GUI and website.
- Addressed a memory-overflow issue when merging contiguous targets.

### Changed
- Upgraded the DockerHub image with the latest fixes.

[Unreleased]: https://github.com/pinellolab/crisprme-plus/compare/v2.8.1...HEAD
[2.8.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.8.1
[2.8.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.8.0
[2.7.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.10
[2.7.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.9
[2.7.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.8
[2.7.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.7
[2.7.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.6
[2.7.5]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.5
[2.7.4]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.4
[2.7.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.3
[2.7.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.2
[2.7.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.1
[2.7.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.7.0
[2.6.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.2
[2.6.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.1
[2.6.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.6.0
[2.5.5]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.5
[2.5.4]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.4
[2.5.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.3
[2.5.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.2
[2.5.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.5.0
[2.4.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.4.0
[2.3.3]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.3
[2.3.2]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.2
[2.3.1]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.1
[2.3.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.0
[2.2.0]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0
[2.2.0-alpha.30]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.2.0-alpha.30
[2.1.13]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.13
[2.1.12]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.12
[2.1.11]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.11
[2.1.10]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.10
[2.1.9]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.9
[2.1.8]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.8
[2.1.7]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.7
[2.1.6]: https://github.com/pinellolab/crisprme-plus/releases/tag/v2.1.6

<!-- tierref track: [2.8.1] entries, merged by the docs track. -->

**For the docs track: three entries, all under `### Fixed`; one clause in the first holds only if
the wset track's Settings change ships.** The clause "and, in this release, every index built or
built again in the web Settings page" (in the sentence after the bold title) is true because the
wset track makes a variant build in the web Settings page pass the dataset's sample list as
`--samplesID` (the first entry under Fixed in `changelog/wset.md`). If that change does not ship
in 2.8.1, replace these two lines of the first entry

> `NRG_3_hg38+hg38_HGDP1kGP` among them), every index built with `build-index-only --samplesID`,
> and, in this release, every index built or built again in the web Settings page. A sites-only

with

> `NRG_3_hg38+hg38_HGDP1kGP` among them) and every index built with `build-index-only --samplesID`.
> A sites-only

With these entries in, the text that the preamble of `changelog/wset.md` says to add if this fix
does not land is not needed. The paragraph "Verified" of `docs/INDEX_PROVENANCE.md` says that a
search "by 2.8.1" matches the gate's report in every column but `Annotation_COSMIC` and gives its
md5: that search ran before these two fixes, so it should say a release candidate; a search by
2.8.1 as released differs from that report in what these entries name. So does this section's
entry on the rebuilt `NNN_3` archives of `CHANGELOG.md`, in its sentence "Searched genome-wide by
2.8.1 (NRG `--prune-pam`, 4 mm + 1 DNA + 1 RNA bulge, per-sample), the rebuilt pair matches the
gate's report in every column but `Annotation_COSMIC`": it should say a release candidate too. The
numbers in the first two entries are from the default index `NRG_3_hg38+hg38_HGDP1kGP` as
downloaded, NRG, 4 mismatches, 1 DNA and 1 RNA bulge: a per-sample search with guide
ACTGAAATCTGTAAGCAGGC, a population-level one with guide CTAACAGTTGCTTTTATCAC, and one of each with
guide CTCTCAGCTGGTACACGGCA, each run on aa9c54b and on these fixes alone: the fixes of the merge step
(`changelog/wset.md`) and of indel rows (`changelog/indel.md`) came after them, and the second gives
another count of indel rows for the same search. Their counts are not of the release as a whole.
The third entry's are with the merge step's fix, before that of indel rows.

### Fixed

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
  population-level ones since 2.5.2. On the default index, a per-sample search (guide
  ACTGAAATCTGTAAGCAGGC, 4 mismatches, 1 DNA and 1 RNA bulge) had 71,110 such rows of 189,949 for
  `highest_CFD` (68,348 for `fewest_mm+b`, 65,294 for `highest_CRISPR_BULGE`); for example
  chr7:66994206 (-) read `ACTGAAATCTGTAAGCAGGTGGG` with REF seed / non-seed 0 / 0 for its one
  mismatch. The tier path now writes these columns as the dictionary path does: each protospacer
  base that differs from the guide in lower case (`ACTGAAATCTGTAAGCAGGtGGG`, 1 / 0), except at a
  DNA-bulge column; the PAM stays in upper case.

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
  mismatches and bulges, so no row is written for a window on the minus strand now. Measured on the
  default index, guide CTCTCAGCTGGTACACGGCA, 4 mismatches, 1 DNA and 1 RNA bulge, genome-wide: per
  sample, 5 of 309,218 rows of alternative alignments go (chr1:147073620, chr16:12583323,
  chr16:85966126, chr6:5506134, chrX:150119503, none of them a site's best); `Other_motifs` goes
  down by one on the 5 rows of those sites in the table and the database, of 128,730, and on 76 rows
  of alternative alignments, one of which, at chr16:12583323, reads strand `+` for its CRISPR-Bulge
  block, where it paired the wrong row (`*_observed_haplotypes_raw` follows); nothing else changes
  in the table, in `observed_haplotypes.tsv` or in the tables of the report. At population level
  `Other_motifs` goes down by one on the same 5 rows, of 150,103, and one row's CRISPR-Bulge score
  moves in its last digit (chr16:1434667, `3.58341e-18` to `3.58338e-18`, as the scorer's batches
  shift); nothing else changes. With more mismatches and bulges such a row can be a site's best, and
  a site can have several. On chr19 and chr21 with 5 mismatches, 2 DNA and 2 RNA bulges, this guide
  and CTAACAGTTGCTTTTATCAC, per sample, 6 of 46,829 sites reported one as their best:
  chr19:28563886, with a CFD of 0.007, now reports a variant alignment with 0.006, and
  chr19:44589953, chr19:50205996, chr21:14062304, chr21:36964620 and chr21:45983672, with a CFD of
  0, now report a variant alignment and `Not_found_in_REF` `y`, as no reference alignment is left at
  the site (chr21:45983672 on strand `+`). `Other_motifs` goes down at 20 sites, by 1 to 10; 58 of
  316,189 rows of alternative alignments go, and at the sites that lost one the others pair
  differently across the three criteria; and the CRISPR-Bulge scores of 24 rows move in their last
  digit. At population level, with `hg38_mega` beside the default index, 5 sites change a best
  alignment so (4 for `highest_CFD`), and a sixth reads `y` for `Not_found_in_REF`. A search of the
  dictionaries, or of the reference alone, had no such row and is unchanged. Covered by
  `PostProcess/test_tier_ref_columns.py`.

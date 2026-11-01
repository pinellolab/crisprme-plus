<!-- tierref track: [2.8.1] entry, merged by the docs track. -->

**For the docs track: one entry, under `### Fixed`; one clause in it holds only if the wset
track's Settings change ships.** The clause "and, in this release, every index built or built
again in the web Settings page" (in the sentence after the bold title) is true because the wset
track makes a variant build in the web Settings page pass the dataset's sample list as
`--samplesID` (the first entry under Fixed in `changelog/wset.md`). If that change does not ship
in 2.8.1, replace these two lines of the entry

> `NRG_3_hg38+hg38_HGDP1kGP` among them), every index built with `build-index-only --samplesID`,
> and, in this release, every index built or built again in the web Settings page. A sites-only

with

> `NRG_3_hg38+hg38_HGDP1kGP` among them) and every index built with `build-index-only --samplesID`.
> A sites-only

With this entry in, the text that the preamble of `changelog/wset.md` says to add if this fix does
not land is not needed. The numbers in the entry are from the default index
`NRG_3_hg38+hg38_HGDP1kGP` as downloaded, NRG, 4 mismatches, 1 DNA and 1 RNA bulge: a per-sample
search with guide ACTGAAATCTGTAAGCAGGC and a population-level one with guide
CTAACAGTTGCTTTTATCAC, each run on aa9c54b and on this fix.

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
  `Non_seed_mismatches+bulges_REF` counted only the row's bulges and left out the reference's
  mismatches, in all three criteria (`highest_CFD`, `fewest_mm+b`, `highest_CRISPR_BULGE`). A
  search of an index's dictionaries marks them in lower case, as CRISPRme always has. The
  fewest mismatches and bulges plot (`imgs/CRISPRme_fewest_top_1000_log_for_main_text_<guide>`,
  `.png` and `.pdf`, and the same plot on a Personal Risk Card in the web) drew the reference
  point of such an off-target at its bulge count, 0 for one without a bulge, so a variant that
  removes mismatches looked like one that adds them. Per-sample searches have had this since
  2.4.0, population-level ones since 2.5.2. On the default index, a per-sample search
  (guide ACTGAAATCTGTAAGCAGGC, 4 mismatches, 1 DNA and 1 RNA bulge) had 71,110 such rows of
  189,949 for `highest_CFD` (68,348 for `fewest_mm+b`, 65,294 for `highest_CRISPR_BULGE`); for
  example chr7:66994206 (-) read `ACTGAAATCTGTAAGCAGGTGGG` with REF seed / non-seed 0 / 0 for its
  one mismatch. The tier path now writes these columns as the dictionary path does: each
  protospacer base that differs from the guide in lower case (`ACTGAAATCTGTAAGCAGGtGGG`, 1 / 0),
  except at a DNA-bulge column; the PAM stays in upper case. The pipeline breaks ties between
  rows by their whole line, which includes this column, so two more things can change. Rows with
  the same CFD score can come out in another order among themselves, in the results, in the
  report's tables and figures and in the perfect-match banner (on the per-sample search above,
  5,383 of its 189,949 rows are at another place in the results and otherwise the same; 4,598 of
  205,415 on a population-level search of the same index). And where two candidate alignments of
  a site tie exactly on start, score, mismatches and bulges, number of carriers and aligned
  sequence, the other one of the two can now be the one reported, with its own variant, rsID,
  allele frequency, carriers, PAM creation and REF and risk scores, the former one moving to the
  alternative alignments (which a per-sample search lists). That happened at 15 sites for
  `highest_CFD`, 17 for `fewest_mm+b` and 16 for `highest_CRISPR_BULGE` on the per-sample search
  above, and at 4, 5 and 7 on the population-level one; in every one of these cases the two were a
  row of SNPs and a row with an indel (alone or with a SNP) that give the same alignment, and the
  row with the indel is now the one reported. The counts made from the reported rows' carriers
  (`*.summary_by_samples.*`, the population distribution plots and the report's population panel)
  move by a few accordingly, and a PAM-creation count of `*.summary_by_guide.*` by one. The set of
  alignments reported (best and alternative), every score and every other column are otherwise
  unchanged, and the results of a search of an index without the tiers, or of a reference-only
  search, are the same byte for byte. To correct results made before, run the search again
  (`generate-report` alone reads the table as it is). Covered by
  `PostProcess/test_tier_ref_columns.py`.

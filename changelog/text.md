<!-- text track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the merge, not for the changelog.** The first Fixed entry completes the Changed entry
"`--per-sample` on an index without the dict-less tiers now says what it reports" (the banner,
new in 2.8.1), and can go right after it instead. That banner now also gives the known limit of
the Changed entry below; what the entry says of the banner still holds. Three references name an
entry by its place: the first Fixed entry's "the last entry below" and the second's "the next
entry" are the entry on tied rows, the last of this fragment, and the second's "the entry above"
is the first; keep them true where the entries land.

### Changed

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

### Fixed

- **The report of a `--per-sample` search on an index without the dict-less tiers says what the
  search reported.** On such an index (one built without `--samplesID`, such as the one the
  benchmark `complete-test` builds), the report's Search mode row now says what
  `complete-search`'s banner says: the carrier samples are named, with the allele frequencies the
  VCF gives, and no CONFIRMED/PUTATIVE cis phasing, SNP+SNP co-occurrence or exact joint allele
  frequency is reported, as they come from the tiers; a dense window can still be reported by its
  greedy worst-case alignment (the known limit under Changed); and to get the tiers, rebuild the
  index with `build-index-only --samplesID` or download a prebuilt one. The texts that made the
  same claims elsewhere in the report now say what holds for those datasets: the Allele-frequency
  basis row and the MAF note under the top-1000 tables (a variant's frequency is the `AF` its VCF
  gives, as rounded there, not AC/AN, and one whose `AF` reads 0 has no frequency and is blank),
  the legend's MAF entry (no joint allele frequency is computed, and a variant whose `AF` reads 0
  is left out of a row's minimum), the legend's Observed entry and the caption of the top-1000
  tables (a dense window's greedy alignment lists the carriers of any of its variants, not its
  exact cis carriers). When no dataset of the search has the tiers, the Observed cell of such an
  alignment (a row of two or more SNPs flagged in `High_complexity_region`), in the report's
  tables and in every curated TSV that holds it, reads `N carriers (>=1 of K vars, poss. chimera)`,
  as in a population-level search, and links its site's other greedy alignments as "N other
  variant combinations". No column says which dataset a row comes from, and two datasets can name
  the same samples, so in a search that also has a dataset with the tiers, whose dense windows come
  as observed haplotypes flagged there too, every row keeps its cell, and the legend and the
  caption say why. The report read `.search_mode` alone, so every per-sample report said "Every
  per-sample haplotype is reported with CONFIRMED cis phasing, named carrier samples and exact
  joint allele frequency", AC/AN over the genotyped panel, an exact joint allele frequency and
  exact cis carriers. The self-test panel's search (the NGG index the README's self-test builds,
  `NGG_2_hg38_chr22+hg38_1000G`; the README's guide with the `20bp-NGG-SpCas9` PAM, 4 mismatches,
  1 DNA and 1 RNA bulge) has no CONFIRMED line in any of its files; `Variant_MAF` is NA on 154 of
  its 189 alt rows (3 of 184 with the tiers), those whose variants' AF reads 0 in the VCF, which
  gives it to two decimals; four two-variant rows show the frequency of their commoner variant
  (0.09 to 0.86: chr22:16858220 shows 8.10e-01, from `Variant_MAF` "NA,0.81") for a haplotype
  that one to five individuals carry (9.81e-04 and 1.96e-04 with the tiers); and the 5 greedy
  rows read "16, 40, 22, 18 and 18 cis carriers", three of them with a link to "N other observed
  haplotypes", though none of the 102 individuals they list (114 times) carries every variant of
  a row that lists it. The report cannot tell from the results, nor from the install, whose index
  can have the tiers by the time the report is built again, so `complete-search` now records the
  datasets its banner names in the result folder (`.untiered_datasets`), beside every dataset it
  searches (`.variant_datasets`), as it records `.search_mode`; a report built again from an
  earlier search's folder keeps the old texts. With the tiers this change leaves the report as it
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
  the classic per-sample dictionaries, so a population-level search reports the `AF` the VCF gives
  too: on the self-test panel above, 198 of its 257 alt rows have a blank MAF. Its report said
  AC/AN over the panel all the same (the Allele-frequency basis row, the MAF note, the legend's MAF
  entry), as these texts did whatever the index. The Search mode row, as the banner, advised
  re-running with `--per-sample` for per-sample carriers, CONFIRMED cis phasing and exact joint
  allele frequency; the legend said that `--per-sample` computes the exact joint allele frequency
  and lists the exact cis carriers, and the caption of the top-1000 tables that it resolves the
  true cis haplotypes and carriers. On such a dataset `--per-sample` names the carriers only, with
  neither CONFIRMED cis phasing nor a joint allele frequency, and can still report a dense window
  by its greedy worst-case alignment. `complete-search` now records those datasets
  (`.untiered_datasets`, `.variant_datasets`) for any variant search, per-sample (the entry above)
  or population-level, and the web search job, which runs the search script itself, records them
  the same way. The population-level banner and the report's Search mode row say, after their
  `--per-sample` advice, what `--per-sample` gives on those datasets and how to get the tiers; the
  allele-frequency texts say what they say in a per-sample report; the legend's MAF and Observed
  entries and the caption say what `--per-sample` gives there. A sites-only dataset (one whose VCFs
  have no sample columns), to which a build without `--samplesID` writes no tiers either, is not
  one of those datasets, in the records, the banners or the report: with no carrier to name and no
  allele-frequency registry, its search reports none of its variants, in either mode. With the
  tiers this change leaves the banner and the report as they were, byte for byte (but for the order
  of tied rows, which the next entry changes in every report). Pre-existing: the allele-frequency
  texts have read so since 2.4.0, the banner and the Search mode row since 2.5.3, the legend's
  entries since 2.6.0 and the caption since 2.8.0. Covered by `PostProcess/test_generate_report.py`,
  `PostProcess/test_build_index_only_atomic.py` and `PostProcess/test_web_search_reuse.py`.
- **Rows of one score keep the table's order in the report's ranked tables and plots, and in the
  results page's plots.** The top-1000 by CRISPR-Bulge (`top1000_crispr_bulge.tsv` and its
  table), the top-1000 by CFD (`top1000.tsv` and its table) and the Section-2 plots' top-N were
  sorted with pandas' default sort, a quicksort, which is not stable: rows of exactly one score
  (common: identical aligned pairs at different sites score the same, and an out-of-domain
  CRISPR-Bulge score is -1) came out in an order that the whole table decided, so rows added
  elsewhere renumbered them; between two runs, 8 to 99 rows of `top1000_crispr_bulge.tsv` changed
  only their rank, each within an exact tie. These sorts are now stable, so tied rows keep the
  order of `integrated_results.tsv`. The self-test panel's `--per-sample` searches without and
  with the dict-less tiers (the searches above) differ by 5 rows (the known limit under Changed):
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

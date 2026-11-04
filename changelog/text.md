<!-- text track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track. -->

**For the merge, not for the changelog.** The first Fixed entry completes the Changed entry
"`--per-sample` on an index without the dict-less tiers now says what it reports" (the banner,
new in 2.8.1), and can go right after it instead. That banner now also gives the known limit of
the Changed entry below; what the entry says of the banner still holds.

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
  `High_variant_density_region`, with the carriers of any allele it takes. On the self-test's
  chr22 panel (the benchmark's phased 1000 Genomes VCF; the README's guide, 4 mismatches + 1
  DNA + 1 RNA bulge) that gives 5 rows (chr22:22410074, 22687072, 22687116, 22698190 and 22888934,
  all on the - strand, 8 to 13 variants each, 4 mismatches and 2 bulges, CFD at most 0.001), and
  each of the 114 carriers they list lacks one of its row's variants: 342 of the run's 56,730
  SNP-row carrier listings (a row lists its carriers once for each of its three alignments). After
  `build-index-only --samplesID` the same search returns exactly those 5 rows fewer, and every
  carrier it lists has all of its row's variants on one phased haplotype. The `complete-search`
  banner and the report of such a search (its Search mode row, the legend's Observed entry and the
  caption of the top-1000 tables), the `--per-sample` help, the README, METHODS sections 5 and 7
  and `docs/HGDP_1KGP_PHASED_INDEX_METHODS.md` now say so, and to rebuild the index with
  `build-index-only --samplesID` or use a prebuilt one. Results are unchanged. Pre-existing: 2.8.0
  behaves the same; it is fixed in the next release. Covered by
  `PostProcess/test_build_index_only_atomic.py` (the banner) and
  `PostProcess/test_generate_report.py` (the report).

### Fixed

- **The report of a `--per-sample` search on an index without the dict-less tiers says what the
  search reported, in every section.** On such an index (one built without `--samplesID`, such as
  the one the benchmark `complete-test` builds), the report's Search mode row now says what
  `complete-search`'s banner says: the carrier samples are named, with the allele frequencies the
  VCF gives, and no CONFIRMED/PUTATIVE cis phasing, SNP+SNP co-occurrence or exact joint allele
  frequency is reported, as they come from the tiers; a dense window can still be reported by its
  greedy worst-case alignment (the known limit under Changed); and to get the tiers, rebuild the
  index with `build-index-only --samplesID` or download a prebuilt one. The texts that made the
  same claims elsewhere in the report now say what holds for those datasets: the
  Allele-frequency basis row and the MAF note under the top-1000 tables (a variant's frequency is
  the `AF` its VCF gives, as rounded there, not AC/AN, and one whose `AF` reads 0 is blank), the
  legend's MAF entry (no joint allele frequency is computed, and a variant whose `AF` reads 0 is
  left out of a row's minimum), and the legend's Observed entry and the caption of the top-1000
  tables (a dense window's greedy alignment lists the carriers of any of its variants, not its
  exact cis carriers). The report read `.search_mode` alone, so every per-sample report said
  "Every per-sample haplotype is reported with CONFIRMED cis phasing, named carrier samples and
  exact joint allele frequency", AC/AN over the genotyped panel, an exact joint allele frequency
  and exact cis carriers. The self-test panel's search (the README's guide, 4 mismatches + 1 DNA +
  1 RNA bulge) has no CONFIRMED line in any of its files; `Variant_MAF` is NA on 154 of its 189
  alt rows (3 of 184 with the tiers), those whose variants' AF reads 0 in the VCF, which gives it
  to two decimals; four two-variant rows show the frequency of their commoner variant (0.09 to
  0.86: chr22:16858220 shows 8.10e-01, from `Variant_MAF` "NA,0.81") for a haplotype that one to
  five individuals carry (9.81e-04 and 1.96e-04 with the tiers); and the 5 greedy rows read "N cis
  carriers" though none of the 114 carriers they list carries every variant of its row. The
  report cannot tell from the results, nor from the install, whose index can have the tiers by
  the time the report is built again, so `complete-search` now records the datasets its banner
  names in the result folder (`.untiered_datasets`), as it records `.search_mode`; a report built
  again from an earlier search's folder keeps the old texts. With the tiers, and for a
  population-level search, the report is the same, byte for byte, and no result changes. The
  `--per-sample` texts of `complete-search`'s help and of the README, the README's MAF column and
  METHODS sections 3 and 7 now say so too. The web search form offers per-sample only on an index
  with a genotype store, so a web search is never such a search. Pre-existing: the row has read so
  since 2.5.3, the allele-frequency texts since 2.4.0, the legend's entries since 2.6.0 and the
  caption since 2.7.10. Covered by `PostProcess/test_generate_report.py` and
  `PostProcess/test_build_index_only_atomic.py`.
- **Rows of one score keep the table's order in the report's ranked tables and plots.** The
  top-1000 by CRISPR-Bulge (`top1000_crispr_bulge.tsv` and its table), the top-1000 by CFD
  (`top1000.tsv` and its table) and the Section-2 plots' top-N were sorted with pandas' default
  sort, a quicksort, which is not stable: rows of exactly one score (common: identical aligned
  pairs at different sites score the same, and an out-of-domain CRISPR-Bulge score is -1) came out
  in an order that the whole table decided, so rows added elsewhere renumbered them; between two
  runs, 8 to 99 rows of `top1000_crispr_bulge.tsv` changed only their rank, each within an exact
  tie. These sorts are now stable, so tied rows keep the order of `integrated_results.tsv`. The
  self-test panel's `--per-sample` searches without and with the dict-less tiers differ by 5 rows
  (the known limit under Changed): the 560 sites both reports hold came in another order at 282
  places of `top1000_crispr_bulge.tsv` and 250 of `top1000.tsv`, and now come in the same order.
  On the search without the tiers, whose 565 rows hold 476 CRISPR-Bulge scores of -1 and 310 CFD
  scores of 0, the two tables put 477 and 483 rows at another rank than before, each inside its
  tie. Only the order of tied rows changes, and, where a tie spans the last row of a top-N, which of
  them it holds; no score changes, and no row leaves its tie. The recommended panel, sorted on
  seven keys at once (which pandas sorts stably), and the per-threshold tables, in the table's
  order, were already so and are unchanged. Pre-existing since 2.4.0. Covered by
  `PostProcess/test_generate_report.py`.

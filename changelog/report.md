<!-- report track: entries for CHANGELOG.md [2.8.1], to be merged by the docs track. -->

### Fixed

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

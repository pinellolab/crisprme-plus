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

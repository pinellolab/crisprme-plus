<!-- web2 track: [2.8.1] entries for CHANGELOG.md, to be merged by the docs track (final review
WEB-1, WEB-2, WEB-5, WEB-10, REPORT-POP-DEFAULT-LABEL, REPORT-INDEL-ROSTER-LABEL). -->

**For the docs track: two amendments to entries already in `## [2.8.1]`, then four new entries.**
What the search form says of a dataset added as one multi-chromosome VCF file, or of a sites-only
one, is 2.8.1's own text: the bulge-note entry and the sample-list entry describe it, and with
these fixes they would contradict what 2.8.1 does. So those two are amendments: replace the quoted
text (line numbers are those of `CHANGELOG.md` at 9c6b404) with the text after it. Both entries
are under `### Fixed`.

1. **"The search form's bulge note and defaults did not match the installed indexes, and the docs
   said the web interface never builds an index"**, lines 1420-1423. Replace

   > the installed genome, once). For a panel with no variant index it says that only a 0-bulge search
   > is offered and that the first one prepares the variant genome and a 0-bulge index (hours on a
   > full genome), and to build the index in Settings ("Build an index locally") when the genome and
   > VCF files are installed, or to download one otherwise. The DNA and RNA bulges and the Max edits

   with

   > the installed genome, once). For a panel with no variant index it says that only a 0-bulge search
   > is offered and that the first one prepares the variant genome and a 0-bulge index (hours on a
   > full genome), and to build the index in Settings ("Build an index locally") when the genome and
   > VCF files are installed, or to download one otherwise. A search prepares the variant genome
   > from VCF files one per chromosome, and only a build of the index splits a single
   > multi-chromosome file, as the Settings card says: for a panel added as one such file (an upload
   > or a server folder) the note now says to build its index in Settings first, or to add the panel
   > again as one VCF file per chromosome, and Submit refuses its search with the same message. That
   > search was accepted and failed within seconds (`ERROR: cannot open .../Genomes/hg38/.fa`,
   > `Genome enrichment failed`), the status page saying only that the result encountered some
   > errors and to submit it again, which failed the same way (2.8.0 too). The DNA and RNA bulges and the Max edits

   and in lines 1427-1428 replace

   > and the data setup guide now say what a web search builds (a missing reference index; for a panel
   > without a variant index, a 0-bulge one), as does the message of the submit check. Covered by

   with

   > and the data setup guide now say what a web search builds (a missing reference index; for a panel
   > without a variant index, a 0-bulge one, from its VCF files one per chromosome), as does the
   > message of the submit check. Covered by

   and the entry's last line, 1429, `` `PostProcess/test_web_bulge_notes.py`. `` with
   `` `PostProcess/test_web_bulge_notes.py` and `test/web/test_settings_vcf_datasets.py`. ``

2. **"A VCF dataset added in the web Settings page can be searched from the browser, per sample
   too"**, lines 1488-1490. Replace

   > than letting it fail at its end (for a sites-only dataset, which has no samples to list, it says
   > so, and that a search reports such a VCF's variants only through a per-dataset allele-frequency
   > registry, as the shipped `hg38_mega` index has); and a genotyped dataset with no tiers yet reads

   with

   > than letting it fail at its end (for a sites-only dataset, which has no samples to list, it says
   > so, and that a search reports such a VCF's variants only through a per-dataset allele-frequency
   > registry, as the shipped `hg38_mega` index has; the note under the bulge options says the same of
   > such a dataset, an upload or the `hg38_mega` VCFs a Settings download installs, with or without
   > its index built in Settings, where 2.8.0's offered it a 0-bulge search and a build in Settings, or
   > the bulges of the index built; the web user guide says so, and that the published
   > `NRG_3_hg38+hg38_mega` index is the way to search the mega variants); and a genotyped dataset with no tiers yet reads

The new entries: the three `### Fixed` ones go with the web and report entries of `## [2.8.1]`
(after the bulge-note entry, line 1429, will do), the `### Changed` one at the end of that
section.

### Fixed

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

### Changed

- **In a population-level report an indel (or MNV) row's Observed cell reads `N carrier(s)`**, not
  `N carriers (>=1 of K vars, poss. chimera)` as a multi-variant SNP row does, now that such a row
  can name the SNP beside its indel: its carriers are resolved one by one, CONFIRMED where each is
  placed on its phased haplotype and PUTATIVE where the phase is unknown (an unphased call, an MNV).
  On an install without SNP genotypes (`download --what index --no-genotypes`) a SNP+indel row still
  lists every carrier of the indel, with or without the SNP (PUTATIVE in `indel_snp_cooc.tsv`), as
  in 2.8.0. The README and METHODS section 7 say so.

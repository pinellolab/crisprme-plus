# [2.8.1] entries from the wsearch track (web search form, job status and results pages)

All seven go under **### Fixed** of `## [2.8.1]` in CHANGELOG.md. None is new in 2.8.1: each
was found by the 2.8.1 web walkthrough or its review and is older.

### Fixed

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
  VCF files are installed, or to download one otherwise. The DNA and RNA bulges and the Max edits
  that the cap lowered come back as far as the cap allows them again, unless they were changed by
  hand in between: 1, 1 and 6 by default, on Reference only and on a panel with a 1-bulge index
  alike, also after several panels in a row. The web user guide, the Docker quickstart, the website
  and the data setup guide now say what a web search builds (a missing reference index; for a panel
  without a variant index, a 0-bulge one), as does the message of the submit check. Covered by
  `PostProcess/test_web_bulge_notes.py`.
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

<!-- web3 track: [2.8.1] entries for CHANGELOG.md, from the final validation of 38a7c4e (W1-W4). -->

### Fixed

- **A web search refused for an index a newer CRISPRitz built now says why on its status and
  results pages.** The search stops, as it should, with one ERROR (`<index> was built by a newer
  CRISPRitz (... is in index format 2, ...): this CRISPRme reads index format 1 only ...`, and what
  to do), but the job status page and its results link (`/load` and `/result`) knew only the other
  refusals' ERRORs, and said "The selected result encountered some errors, please remove it and try
  to submit again.": submitted again, the search was refused the same way, and the reason was shown
  only under Show full log. They now say "The search was not run: " and the ERROR, as they do for
  the other refusals. New in 2.8.1. Covered by `PostProcess/test_results_summary_and_status.py` and
  `PostProcess/test_web_results_links.py`.

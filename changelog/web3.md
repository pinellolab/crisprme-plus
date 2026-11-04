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
- **The search form says before Submit that a variant dataset has no sample list, and how to get
  one, also for a variant index built in Settings.** Removing a VCF dataset in Settings takes the
  sample list the page wrote for it (`samplesIDs/<dataset>.samplesID.txt`) and leaves its variant
  index, which the search form still offers, labelled as before. The note under the bulge options
  then offered its search, and Submit refused it, saying to download its index again, which
  installs the list, or to add the file yourself: an index built in Settings cannot be downloaded.
  The note offered the search of any dataset Submit refuses for want of its list (one whose VCFs
  name samples but which has no list yet, too). It now says what Submit says, as it does of a
  sites-only dataset or of one added as a single multi-chromosome VCF file, and for a dataset whose
  VCF folder is not there the message says to download its index again if it came from
  HuggingFace, to add the VCF dataset again in Settings under its name, which writes the list, or
  to add the file yourself. A merged panel, whose list Submit makes from those of its datasets, is
  offered as before. New in 2.8.1. Covered by `PostProcess/test_web_bulge_notes.py`.

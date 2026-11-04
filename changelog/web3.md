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
  one, whether its variant index was built in Settings or downloaded.** Removing a VCF dataset in
  Settings takes the sample list the page wrote for it (`samplesIDs/<dataset>.samplesID.txt`) and
  leaves its variant index, and the search form still offers the dataset if its name begins with
  its genome's and an underscore (as `hg38_<name>`; one named otherwise is not offered until it is
  added again). The note under the bulge options then offered its search, and
  Submit refused it, saying to download its index again, which installs the list, or to add the
  file yourself: an index built in Settings cannot be downloaded, and one downloaded is not
  downloaded again while it is installed (Settings answers "already installed — nothing to
  download", and `crisprme.py download --what index` skips a complete install whatever its list).
  The note offered the search of any dataset Submit refuses for want of its list (one whose VCFs
  name samples but which has no list yet, too). It now says what Submit says, as it does of a
  sites-only dataset or of one added as a single multi-chromosome VCF file, and for a dataset
  whose VCF folder is not there the message says, if its index came from HuggingFace, to remove
  the index in Settings and download it again there, which installs the list; to add the VCF
  dataset again in Settings under its name, which writes the list; or to add the file yourself. A
  merged panel, whose list Submit makes from those of its datasets, is offered as before. New in
  2.8.1. Covered by `PostProcess/test_web_bulge_notes.py`.
- **Submit with a genome or PAM selected that is not installed says what is missing and where to
  add it, instead of failing.** On a fresh data folder, which holds no genome or PAM until some
  are added in Settings / Data Manager, Load Example selects `hg38` and `20bp-NRG-SpCas9`, and
  Submit's input check opened `PAMs/20bp-NRG-SpCas9.txt` (the default one, too, with no PAM
  selected), which is not there: the check failed with a server error (HTTP 500), no "Missing or
  wrong input" dialog opened, and the page did nothing. It read the genome folder of a genomic
  sequence given by BED coordinates the same way, so a genome not installed failed too, and so did
  coordinates on a sequence whose FASTA the genome folder does not hold, or past the sequence's
  end. The dialog now lists the genome and the PAM that are not installed, with where to add one
  in Settings (Add a reference genome, Define a nuclease / PAM), and the BED coordinates it cannot
  read; the bulge check, which needs both, waits for them. Pre-existing: 2.8.0 did the same.
  Covered by `PostProcess/test_web_bulge_notes.py`.
- **The web interface no longer prints a `BiopythonWarning` as it starts on a Bioconda install.**
  A stray `site-packages/pyproject.toml`, which conda-forge's `dash-bootstrap-components` ships,
  makes Biopython warn at import that it is imported from inside its source tree. `crisprme.py`
  has filtered that warning since 2.7.10, but `web-interface` starts `index.py` in a new
  interpreter, which inherits no filter and whose pages import Biopython. `index.py` now sets the
  same filter before its first import. The Docker image and `install_from_source.sh` remove the
  stray file, so they never showed it. Pre-existing. Covered by `PostProcess/test_app_settings.py`:
  a check of `index.py`'s code, which runs everywhere, and its imports run as `web-interface` runs
  them, with Biopython beside a stray `pyproject.toml`, which runs where the web stack and
  Biopython are installed, as in the image.

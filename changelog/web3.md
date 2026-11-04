<!-- web3 track: [2.8.1] entries for CHANGELOG.md, from the final validation of 38a7c4e (W1-W4). -->

**For the merge, not for the changelog:** CHANGELOG.md's [2.8.1] entry "CI runs what it says it
runs" gives web-e2e's in-image step "440 tests, about 100 s". This track adds `test_app_settings`
(25 tests, the last Fixed entry below) and `test_web_shell_safety` (21 tests, the Security
block) to that step, and on this branch the step then ran 488 tests in about 100 s, with its
two allowed skips: correct the count there as these entries are merged.

### Fixed

- **A web search refused for an index a newer CRISPRitz built now says why on its status and
  results pages.** The search stops, as it should, with one ERROR (`<index> was built by a newer
  CRISPRitz (... is in index format 2, ...): this CRISPRme reads index format 1 only ...`, and what
  to do), but the job status page and its results link (`/load` and `/result`) knew only the other
  refusals' ERRORs, and said "The selected result encountered some errors, please remove it and try
  to submit again.": submitted again, the search was refused the same way, and the reason was shown
  only under Show full log. They now say "The search was not run: " and the ERROR, as they do for
  the other refusals. The search form gives no warning of such an index before Submit: which index
  a search reads is decided by the search script, and the form would need its own copy of that
  choice; after Submit, those pages say why. New in 2.8.1. Covered by
  `PostProcess/test_results_summary_and_status.py` and `PostProcess/test_web_results_links.py`.
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
  Biopython are installed: in the image, where `web-e2e.yml` runs it with the other tests that
  skip in the light job.

### Security

- **The search form's genomic-sequence input passes nothing typed to a shell.** For a sequence
  given by BED coordinates, the sequence's `>` header and the sequence name on the coordinates'
  line went as typed into a shell command that ran `bedtools getfasta`: when Submit checked the
  form, as the base-editing window list was filled at each change, and when the search was
  launched. Characters a shell reads in either ran as commands, as the web app's user, which
  matters wherever the web app is served to other people. The header named the BED file written
  in the data folder, so a header holding a path wrote it elsewhere, and the index bedtools wrote
  beside the genome's FASTA was removed afterwards, one that was there before as well. bedtools
  now runs with no shell, on a BED file and a link to the FASTA in a private temporary folder,
  where it writes the index too: the header names no file, the sequence name must be that of a
  FASTA of the genome selected, and the genome an installed one. Pre-existing, in every earlier
  release. Covered by `PostProcess/test_web_shell_safety.py`.
- **Submit's launch of a search runs only what the form offers, and quotes each value it passes
  on.** A request can call the launch with any values, past Submit's check of the form: the
  genome, the PAM and the variant dataset name the folders and files the search reads, and the
  genome, the e-mail address and the number of mismatches went as sent into the shell command
  that starts the search, where characters a shell reads ran as commands. The launch now runs
  only a genome and a PAM that are installed, "Reference only" or a variant dataset the form
  offers for the genome, and whole numbers of mismatches and bulges, and it quotes each value of
  that command, which is unchanged for the values the form sends. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **A personal-assembly search takes only guides of the nucleotide codes a reference search
  keeps.** Submit on the search form's Personal assembly tab wrote the guides of the PAM's length
  into the search's guide file as sent, upper-cased, and the search passes each guide of that
  file to a shell command: characters a shell reads in a guide ran as commands, as the web app's
  user, and a request can send that Submit any text. It now leaves out a guide holding a
  character other than those the Reference genome tab's Submit keeps of a guide (A, C, G, T and
  the IUPAC codes for two or three of them), as it leaves out one not of the PAM's length,
  rather than search a shorter one. A guide of those letters is searched as before; one holding
  an N is now left out too. Pre-existing, since 2.5.5. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **The annotation a search applies is put together with no shell.** It is written as
  `Annotations/.active.<genome>.bed.gz`, and the name of a genome folder, which Settings took as
  typed and which can be put in `Genomes/` by hand, went into the commands that sort, compress,
  move and remove it, run through a shell, where characters a shell reads ran as commands. They
  now run with no shell. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **Settings refuses a name that CRISPRme would take into a shell command.** It took the name of
  a genome (typed, or that of the file uploaded), of a VCF dataset and of a nuclease with any
  character but a space, `+`, `/`, `\` and `..`, that of a file uploaded with any but `..`, and a
  genome download's URL as typed. CRISPRme takes these names into shell commands as they are (a
  search's post-analysis a genome's, a dataset's and a PAM file's, an index build a genome's and a
  dataset's, a genome download its own and the name its URL ends in), where characters a shell
  reads ran as commands, as the web app's user. Settings changes data in local mode only, which is
  how CRISPRme runs by default. It now refuses, and says why, a name that holds a blank, a
  character a sequence's name may not hold (`|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`, a
  backquote or `#`), `*`, `?` or `[`, or that starts with `-` or `.`: any name it takes, the name
  a URL ends in, and that of a file uploaded, but for a personal-assembly bundle, whose file's
  name only says how it is packed. Pre-existing. Covered by
  `PostProcess/test_settings_progress.py`.
- **Remove result, on a search's status page, removes only that search's folder.** It removed the
  folder the page's address named, through a shell: the job id there could run commands, or name
  a folder elsewhere, outside the data folder too, which was removed. It now removes, with no
  shell, only a folder of `Results/` whose id has the form the search form gives its jobs.
  Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **The results pages pass no value of a link or a table to a shell.** The guide, sample,
  chromosome and position a results link names, the guide and number of mismatches the results
  tables send back, a sample's name, as its VCF gives it, and the name of a search's results
  file, which holds the genome's, went into the shell commands that zip a table for download,
  draw the Graphical Reports and Personal Risk Card plots, and look up the targets at a
  position, where characters a shell reads ran as commands. Those commands now run with no
  shell, but for the look-up at a position, which now takes only a guide of the search, a
  position in digits and a chromosome name of letters, digits and `_.:,+=@%` (for another it
  says "The selected position does not exist"), and the results pages read only a folder of
  `Results/` whose id has the form the search form gives its jobs. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py` and `PostProcess/test_web_results_links.py`.

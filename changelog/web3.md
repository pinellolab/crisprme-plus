<!-- web3 track: [2.8.1] entries for CHANGELOG.md, from the final validation of 38a7c4e (W1-W4). -->

**For the merge, not for the changelog:** CHANGELOG.md's [2.8.1] entry "CI runs what it says it
runs" gives web-e2e's in-image step "440 tests, about 100 s". This track adds `test_app_settings`
(25 tests, the last Fixed entry below) and `test_web_shell_safety` (49 tests, the Security
block) to that step, and on this branch the step then ran 516 tests in about 100 s, with its
two allowed skips: correct the count there as these entries are merged. The indel search's
chromosome-name check (its Security entry) adds three tests to `test_indel_empty_chromosome`,
and the Settings, personal-assembly and History entries add tests to `test_settings_progress`,
`test_personal_assembly` and `test_history_page_assembly`, in unit-tests' own job.

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
- **A personal-assembly search takes only guides of the nucleotide codes a reference search keeps,
  and only a PAM and files the form offers.** Submit on the search form's Personal assembly tab
  wrote the guides of the PAM's length into the search's guide file as sent, upper-cased, and the
  search passes each guide of that file to a shell command: characters a shell reads in a guide
  ran as commands, as the web app's user, and a request can send that Submit any text. It now
  leaves out a guide holding a character other than those the Reference genome tab's Submit keeps
  of a guide (A, C, G, T and the IUPAC codes for two or three of them), as it leaves out one not
  of the PAM's length, rather than search a shorter one. A guide of those letters is searched as
  before; one holding an N is now left out too. That Submit searched, as sent, any PAM file and
  any folders and files for the individual's, outside the data folder too, and the search takes
  the PAM file's name, and that of a genome folder of `Genomes/`, into shell commands: it now
  searches only a PAM the form offers, with the files of the individual picked, as installed.
  Pre-existing, since 2.5.5. Covered by `PostProcess/test_web_shell_safety.py`.
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
- **A job's status and results pages read, write and remove only that job's folders.** These
  pages take the job id from their address, which a client can set to any text, in local mode
  and on a website alike, and joined it as it was into the paths of the files they use: an id
  holding a path reached a folder elsewhere, outside the data folder too. Remove result of a
  personal-assembly search removed that folder, and the haplotype folders its `.Params.txt`
  names, joined the same way; Show full log, and the status of a personal-assembly search, read
  the log files there and showed them; the opening of a job's pages read its `.Params.txt`
  there; and the results pages wrote the table filter they keep (`.dropdown.json`) there, and
  their tables, sample and population lists and download links read the files there. Each now
  takes only a folder of `Results/` whose id has the form the search form gives its jobs, as
  Remove result of a search does (above), and otherwise changes nothing, shows no log, or says
  that there is no such job. History, which lists every search of `Results/`, shows one named
  otherwise (a command-line search's `--output` may be) with no link, and says why its pages do
  not open and what to do (`docs/crisprme_web_interface_user_guide.md`). Pre-existing. Covered
  by `PostProcess/test_web_shell_safety.py` and `PostProcess/test_history_page_assembly.py`.
- **Remove result removes only a search its status page offers it for.** A search's status page
  offers Remove result once the search has failed, and a personal-assembly search's once it has
  failed or is reconciled, but the button's request, which a client can send for any job,
  removed the search of the page's address whatever its state: a finished, running or queued
  search, and, through a personal-assembly search's button, a search of either kind. History
  lists every search in local mode, which is how CRISPRme runs by default, and on a website a
  search's link names it. Each now removes only a search in a state its page offers it for, and
  otherwise changes nothing. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **Remove installed data, in Settings, removes only an item it lists.** It took the item as a
  request sent it, and a name of `..` named the data folder (or, through a folder of it that is a
  link to other storage, the folder there that holds the link's target), and `.` or no name the
  folder of every item of that kind: all of it was removed, every search's results with the data
  folder. It now measures and removes only an item the card lists, and says another is not
  found. Settings changes data in local mode only, which is how CRISPRme runs by default.
  Pre-existing, since 2.2.0. Covered by `PostProcess/test_settings_progress.py`.
- **A VCF upload in Settings goes only into a dataset uploaded or downloaded there, and takes a
  dataset's own files only.** An upload into a dataset registered as a server folder was written
  through the registration into that folder, under any file name, a new file or one replacing a
  file there: a request could write any file into any folder the web app's user can write. A
  registered folder is now read where it is and nothing is uploaded into it, and an upload takes
  `.vcf.gz` files and their `.tbi` or `.csi` index only (`docs/SETTINGS_DATA_MANAGER.md`).
  Pre-existing, since 2.2.0. Covered by `PostProcess/test_settings_progress.py`.
- **Settings publishes no index; an index is published from the command line.** With a Hugging
  Face token in the web app's environment (maintainer mode), Settings' Publish took any index
  name, and `publish-index` archives the folder that name leads to, beside it, and uploads the
  archive with that token to the shared repository, over a file of the same name there: a name
  holding a path had any folder of the server archived and uploaded. Taking only an index its
  list offers did not bound that, as a request has the list offer what it installs through the
  page: a VCF dataset of any server folder registered, or of a VCF uploaded, and its index
  built, whose archive holds the dataset's per-sample dictionaries and sample list. Settings now
  publishes nothing, and Retry runs no publish an earlier version started; in maintainer mode
  the page says to publish with `crisprme.py publish-index`, which is unchanged
  (`docs/SETTINGS_DATA_MANAGER.md`). Pre-existing, since 2.2.0. Covered by
  `PostProcess/test_settings_progress.py`.
- **Settings builds an index only of a genome, PAM and VCF dataset its card offers.** "Build an
  index locally" took them as a request sent them, and a path in one built an index of a folder
  outside the data folder. Pre-existing, since 2.2.0. Covered by
  `PostProcess/test_settings_progress.py`.
- **A personal-assembly bundle names its own files only.** The genome folder, chain and
  chromAlias file of each haplotype are paths its `metadata.json` gives, which an imported archive
  brings, and they were taken as they were: a bundle naming files and folders outside its own
  folder was imported, its individual offered as complete, and its search read them, its results
  holding what it found there. Such a bundle is now refused, and one imported before is listed
  as incomplete and not offered (`docs/crisprme_web_interface_user_guide.md`). Pre-existing,
  since 2.6.0. Covered by `PostProcess/test_personal_assembly.py` and
  `PostProcess/test_web_shell_safety.py`.
- **Remove of a personal assembly in Settings removes that bundle's folder, and of a pairing of
  the older layout its own markers only.** Settings listed a bundle under the individual its
  `metadata.json` names, which an archive imported there brings: one naming another individual
  took that one's place in the Individual picker and the Remove card, and one naming an
  individual with no folder was removed as a pairing of the older layout, by the files named
  after the assembly names the bundle gives, outside the data folder too. A bundle is now listed
  under the name of its folder, as `assembly-search --assembly-individual` takes it, and only a
  pairing of the older layout is removed as one (`docs/crisprme_web_interface_user_guide.md`).
  Pre-existing, since 2.6.0. Covered by `PostProcess/test_personal_assembly.py` and
  `PostProcess/test_settings_progress.py`.
- **The search form reads only what it offers.** The callbacks that fill in and check the form
  (the variant datasets offered for a genome, the analysis modes, the bulge options and the note
  under them, and Submit's check) took the genome, PAM and variant dataset a request sends into
  the paths of the folders and files they list and read: a path in them reached outside the data
  folder, and what was there changed what they answered, in local mode and on a website alike.
  They now answer for a value the form does not offer as for none, and Submit lists a variant
  dataset that is not installed for the genome, as it does a genome or PAM that is not installed.
  Submit's check of a dataset of one multi-chromosome VCF, new in 2.8.1, looks for its variant
  index with the PAM file of an installed PAM only, too. Pre-existing. Covered by
  `PostProcess/test_web_shell_safety.py`.
- **The indel steps of a search pass no chromosome name of a variant dataset to a shell.**
  A search with bulges takes each chromosome name of a variant dataset, read from the CHROM
  column of its VCFs (or from the per-chromosome indel logs named after it), into the shell
  commands that search and merge the indel off-targets and run the indel post-analysis, where
  characters a shell reads ran as commands, as the web app's user. A VCF dataset added in
  Settings carries its own CHROM column, which this did not check, and Settings changes data
  in local mode only, which is how CRISPRme runs by default. The indel search and its
  post-analysis now refuse a chromosome name that holds a blank or a character a sequence's
  name may not hold (the same characters Settings refuses in a name, above), and say so; a
  real chromosome name holds none, so a normal dataset is searched as before. Pre-existing.
  Covered by `PostProcess/test_indel_empty_chromosome.py`.
- **The post-analysis of a search passes no FASTA file's name to a shell, and a personal-assembly
  bundle brings only FASTA files named as a search can take them.** The post-analysis takes the
  name of each FASTA file of the genome folder searched, as its chromosome's, into shell
  commands as it is, where characters a shell reads ran as commands, as the web app's user. A
  genome added in Settings has its files named after their sequences, whose names are checked,
  but a personal-assembly bundle imported there brought its genome folders with the names its
  archive gives, which were not, and Settings changes data in local mode only, which is how
  CRISPRme runs by default. A bundle holding a FASTA file whose name holds a blank or a character
  a sequence's name may not hold (`|`, `;`, `'`, `"`, `$`, `(`, `)`, `&`, `<`, `>`, a backquote or
  `#`) is now refused, and the import says why; the post-analysis refuses such a name too, and
  says so, and runs its step for each chromosome with no shell. A chromosome's off-targets are
  kept only from the file named after its sequence, whose name holds none of these, so a normal
  genome or bundle is searched as before. That step also looked a chromosome's off-targets up
  with its name as a shell word, where `*`, `?` and `[` (an HLA allele's name holds `*`) made a
  pattern of the files of the folder it runs in, outside the data folder, which were then read:
  it now looks the name up as it is. Pre-existing (bundles since 2.6.0). Covered by
  `PostProcess/test_personal_assembly.py`, `PostProcess/test_settings_progress.py` and
  `PostProcess/test_web_shell_safety.py`.
- **The results pages' tables bind the values of a link, a table and a form to their queries of
  a search's database.** The guide and sample of a results link, the guide and numbers of
  mismatches and bulges of a Show Targets link, the chromosome and positions typed in Query
  Genomic Regions, the guide of the summary table's row and a sample's name, as its VCF gives
  it, went as they were into the SQL that reads a search's results, and Custom Ranking wrote the
  columns it groups by, the order and the thresholds of its form into its SQL too: a quote in
  such a value ended the text it was written in, and the rest was read as SQL, which could widen
  the query, on the search's own database, or end it in an error. Each value is now bound to
  the query; Custom Ranking takes only the columns, orders and thresholds its form offers, and
  the positions of a region and the numbers of a link are whole numbers. A sample whose name
  holds a quote, whose table and Personal Risk Card failed, now has them. Pre-existing. Covered
  by `PostProcess/test_web_shell_safety.py`.
- **Two download links of the results pages look only for a file of their search's folder.**
  The Summary by Sample and Show Targets links offer their zip once it is there, under a name
  the page gives, which a client can set to any text, and they looked for it wherever that name
  led, outside the data folder too: whether a file was there showed in their answer, though
  nothing was served (the download route serves files of `Results/` only). They now look only in
  the search's folder, as that route would serve the file, and answer for any other name as for
  a zip not made yet. Pre-existing. Covered by `PostProcess/test_web_shell_safety.py`.
- **A saved e-mail password goes only to the mail server it was saved for.** Settings' Email
  notifications kept the saved password when a save named another server, port, SSL setting or
  sender address with the password field empty, as the page shows it, and the notification a
  search sends at its end logs in to the server saved, as the sender saved, with that password:
  a request, which can save any server in local mode (how CRISPRme runs by default), had the
  password sent to a server of its choosing. The password of `CRISPRME_SMTP_PSW` went the same
  way to a server the settings named with no password of their own. A save that changes any of
  the four now asks for the password again, and `CRISPRME_SMTP_PSW` is used with the default
  server (`smtp.gmail.com`) only (`docs/crisprme_web_interface_user_guide.md`). Pre-existing,
  since 2.2.0. Covered by `PostProcess/test_web_shell_safety.py`.

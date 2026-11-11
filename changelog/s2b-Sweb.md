Notes for the merge (not changelog text): the web app's links (S-R5-2 of the review of lane S,
round 5, which its refix documented, left to Luca: "ok to fix the symbolic link in the next
release"), on lane S's head c917347. The block goes at the end of Security in [Unreleased]. One
sentence of lane S's fragment (changelog/s2b-S.md, the "Three temporary names" entry), which
said the web app reads through such links, was corrected in place.

### Security
- **The web app follows no link another user put in its data folder.** A member of the data
  folder's group can write `Results/`, `Settings/` and the job folders in them under a umask of
  002, and the web app listens on every network interface. It followed a link put there: the
  Settings FAILED banner and Show full log, a search's log panel and the results downloads
  (`GET /Results/<path>`, a link in a job folder, or `Results/<name>` a link to a folder of the
  user's) showed and served any file of the user's to anyone who could reach the web app, and a
  Settings job wrote its log through one, over a file of the user's elsewhere (2.8.1 too). The
  inventory of what the web app reads, serves and writes found more of the kind, each fixed
  too: the status page's banner, the History table, the index labels and PAM files of the
  Settings page and the search form, the VCF headers, sample lists and markers, the results
  pages' tables, images, database and archives, and the email notification's files; the
  chunked upload (the browser's data went into a `.part` a link named), the annotation upload
  (into `<name>.uploading`), the email settings (the password into `.email.json.tmp`, or into
  that file when another user made it first), the annotation manifest and active annotation,
  the sample lists (through a temporary name the process id and thread made), the PAM files,
  the tables' filter (`.dropdown.json`), the mode file, the job's files a search writes at
  Submit, the logs of a search (the shell's `1>` and `2>`), the plots, zips and tables the
  results pages write (`zip`, `grep >`, `to_csv`), the COSMIC record and the
  concurrent-searches setting, and the download progress a Settings job writes; and the cache
  of the results tables, pickles Flask-Caching loads, kept in `Cache/` of the data folder,
  where another user's pickle ran as the web app's user when it was loaded.

  Each such file now goes through `PostProcess/safe_files.py`: every folder below the data
  folder is opened in turn with `O_NOFOLLOW` (`openat`), so no link is followed at any step;
  a file read must be a regular file with no other name (`fstat`: not a link, a FIFO, a device
  or, in a folder others can write, a hard link); a file written whole is made under a name of
  its own (`O_CREAT|O_EXCL|O_NOFOLLOW`) and renamed over its name, so a link there is replaced,
  not written through, and a log is opened with `O_NOFOLLOW` and truncated only once it is
  known to be a regular file; a folder is removed within the folder that holds it
  (`shutil.rmtree(dir_fd=)`), and an upload is moved only if it is a regular file. The
  download route refuses a link at any step of its path (403), and an empty, `.` or `..` step
  (404). A program the pages run on a job's files (the region page's `grep` and `awk`, the
  plots, the radar chart) reads copies of them or its standard input, with no shell, and
  writes in a folder of the web app's own whose files are then renamed into the job's folder;
  the archives are made in Python. A Settings job's folder is writable by the web app's user
  alone, and the tables' cache is in a folder of the server's own in the temporary folder
  (`tempfile.mkdtemp`, removed at exit; in memory when that folder is writable by others and
  not sticky). A link is followed only where no one else can have put it: a link of the
  user's (or of the data folder's owner, or root) in a folder only they can write, or a sticky
  one of theirs, such as a genome linked to a shared copy, a VCF dataset registered from a
  server folder, or `Results/` or `Settings/` linked onto a volume; in a folder the group can
  write the web app refuses those too, and below `Results/` and `Settings/` it follows none.

  Tests: `test_safe_files` (the unit-test job: each kind of file, link, FIFO and hard link, at
  each step), `test_web_no_links` (web-e2e's image step: the six cases of the live check, and
  each case the inventory added, through the real Dash server's Flask test client and the
  real functions that write, each refused and a legitimate file still read or served), and
  `test_settings_progress` (the Settings page's readers and writers, without Dash). Each fails
  on the code before. `docs/TROUBLESHOOTING.md` §7 says what this protects and what a shared
  server still needs: the search and the commands the web app starts still write through
  links in folders others can write, and the web app has no login.

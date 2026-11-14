"""Whether a variant index and the Dictionaries of its dataset, in a working directory, come from
one build of the dataset (CRISPRme 2.10.0).

A variant search reads the variant index it takes and, in the working directory, its dataset's
Dictionaries (``Dictionaries/``: the registry, the genotype stores, the indel logs): the alleles
come from the index, their carriers and frequencies from the Dictionaries, and both from one
enrichment of the dataset's VCFs. CRISPRme 2.10.0 installs a panel's Dictionaries once, as a
Dictionaries bundle, in place of those an archive of CRISPRme 2.8.1 installed, so a data folder
upgraded from 2.8.1 can hold an index of 2.8.1 (index format 1, ``NRG_3_hg38+hg38_HGDP1kGP``)
beside the Dictionaries of 2.10.0, which code every allele of a split multi-allelic site, where
2.8.1's enrichment coded the last record's alone (1.95 % of HGDP + 1000 Genomes' SNP records);
or, the other way, an index of 2.10.0 beside the Dictionaries an archive of 2.8.1 put back. The
download says so on stdout (crisprme_hf), and a search of such an index ran on: its SNP rows read
with the carriers and frequencies of the other build, until, hours later, its indel search
refused the INDELS index chromosome by chromosome (pool_search_indels' _window_mismatch).

The two builds are told apart by their fake-indel windows (indel_window): CRISPRitz 2.9.1 and
earlier, which built every index of format 1, wrote windows of 25/27, and 2.10.0 writes 30/30.
The indel logs record the windows of the Dictionaries' enrichment, chromosome by chromosome; an
index tells those of the enrichment it was built from (mismatch): each sub-index of the INDELS
index built with it records its own, or is of format 1, of 25/27 (indel_window.index_window).
A search of an index whose windows are not those of its dataset's indel logs is refused before
anything is written (refusal), by complete-search's pre-flight for each dataset it names, and by
the search script, which the web launches, before it searches the dataset: the check the indel
search makes of each chromosome (pool_search_indels' _window_mismatch), made first. Each
release's own pair runs, as does an index of 2.10.0 built from the enriched-genome kit of 2.8.1
beside that kit's Dictionaries (25/27 on both sides). Where a side holds no window (indel logs without a row, as a
dataset without indels has; an index without its INDELS index), nothing tells the builds apart,
and nothing is refused.

    python index_dictionaries.py <index folder> <working directory> <dataset>

prints the refusal, one line, and exits 3; or prints nothing and exits 0, an empty <index folder>
(no installed index) included. Stdlib, and the repository's indel_window (and crisprme_hf, for
the name of the index of 2.10.0 that replaces one of 2.8.1).
"""

import os
import shlex
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import indel_window  # noqa: E402

# The windows of the two enrichments a user meets, and whose data they are
_RELEASE = {indel_window.OLD_WINDOW: "CRISPRme 2.8.1 or earlier",
            indel_window.WINDOW_2100: "CRISPRme 2.10.0"}

# Where docs/TROUBLESHOOTING.md says what to do, with the commands
_SECTION = "docs/TROUBLESHOOTING.md section 16"


def sub_indexes(index_dir):
    """[(chromosome, folder)] of the sub-indexes of the INDELS index built with variant index
    ``index_dir`` (``<index>_INDELS``: ``<PAM>_<N>_fake<chromosome>`` a chromosome), sorted."""
    indels = os.path.normpath(index_dir) + "_INDELS"
    try:
        names = sorted(os.listdir(indels))
    except OSError:
        return []
    return [(name.rsplit("_fake", 1)[1], os.path.join(indels, name)) for name in names
            if "_fake" in name and os.path.isdir(os.path.join(indels, name))]


def mismatch(index_dir, workdir, dataset):
    """(the windows of the enrichment variant index ``index_dir`` was built from, those of the
    indel logs of ``dataset`` in ``workdir``), two sets, of the first chromosome where they
    differ; None where they agree, or nothing tells. Each sub-index of the INDELS index built
    with the variant index, from one enrichment, is read against its chromosome's log, as the
    indel search reads it (pool_search_indels' _window_mismatch): the window it records, or 25/27
    for one of index format 1, which only CRISPRitz 2.9.1 and earlier build
    (indel_window.index_window), against the log's. A variant index whose INDELS index is not
    there tells nothing here: the search script refuses it, or builds that INDELS index from the
    enrichment in place, whose fake-indel genome it reads against the logs first
    (indel_window.fake_genome_refusal)."""
    logs = os.path.join(workdir, "Dictionaries", "log_indels_%s" % dataset)
    for chrom, sub in sub_indexes(index_dir):
        window = indel_window.index_window(sub)
        if window is None:
            continue
        path = indel_window.log_path(logs, chrom)
        try:
            logged = indel_window.log_windows(path) if path else set()
        except (OSError, ValueError, EOFError):
            continue
        if logged and logged != {tuple(window)}:
            return {tuple(window)}, logged
    return None


def _built(windows):
    """What ``windows`` (a set of windows) say of the build: "25/27 (CRISPRme 2.8.1 or
    earlier)"."""
    parts = []
    for window in sorted(windows):
        who = _RELEASE.get(window)
        parts.append("%s (%s)" % (indel_window.describe(window), who) if who
                     else indel_window.describe(window))
    return " and ".join(parts)


def _successor(name):
    """The index of CRISPRme 2.10.0 that replaces the archive of 2.8.1 named ``name``, or None
    (crisprme_hf.successor_index; None too where it cannot be read)."""
    try:
        import crisprme_hf
        return crisprme_hf.successor_index(name)
    except Exception:  # noqa: BLE001  (no successor named: the advice says the other ways)
        return None


def refusal(index_dir, workdir, dataset):
    """Why a search must not read variant index ``index_dir`` with the Dictionaries of
    ``dataset`` in working directory ``workdir`` (one line, with what to do), or None: the
    windows of the enrichment the index was built from and those its dataset's indel logs
    record are both known and differ (mismatch)."""
    if not index_dir or not os.path.isdir(index_dir):
        return None
    found = mismatch(index_dir, workdir, dataset)
    if found is None:
        return None
    mine, theirs = found
    index_dir = os.path.normpath(index_dir)
    name = os.path.basename(index_dir)
    dicts = os.path.join(os.path.normpath(workdir), "Dictionaries")
    remove = "rm -rf %s %s" % (shlex.quote(index_dir), shlex.quote(index_dir + "_INDELS"))
    what = ("%s is an index of a build of %s whose enrichment wrote fake-indel windows of %s, and "
            "the Dictionaries of %s in %s are of another build, whose indel logs record windows "
            "of %s: a search of it would read the alleles of one build with the carriers and "
            "frequencies of the other." % (name, dataset, _built(mine), dataset, dicts,
                                            _built(theirs)))
    successor = _successor(name)
    if mine == {indel_window.OLD_WINDOW} and successor:
        todo = ("Search %s, the index of CRISPRme 2.10.0 that replaces it (crisprme.py download "
                "--what index --index-name %s), and remove this one (%s); or keep it in a data "
                "folder of its own, with its own Dictionaries" % (successor, successor, remove))
    else:
        todo = ("Put back the Dictionaries it was built with (of a downloaded index: download it "
                "again into a data folder of its own), or remove it (%s) and build it again "
                "from the enrichment those Dictionaries come from (build-index-only); keep the "
                "indexes of each build in a data folder of their own" % remove)
    return "%s %s (%s)." % (what, todo, _SECTION)


def main(argv):
    """``index_dictionaries.py <index folder> <working directory> <dataset>``: the refusal on
    stdout and exit 3, or nothing and exit 0 (the search script's check of a dataset)."""
    if len(argv) != 3:
        sys.stderr.write("usage: index_dictionaries.py <index folder> <working directory> "
                         "<dataset>\n")
        return 2
    why = refusal(*argv)
    if why:
        print(why)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

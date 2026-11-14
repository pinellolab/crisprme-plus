#!/usr/bin/env python3
"""A chromosome with no indel off-target possible must not abort a variant search.

The enricher writes one fake-indel contig per chromosome from the windows of the
indels some sample carries, so a chromosome with none gets a header-only FASTA and
a header-only indel log. buildTST cannot index that genome, which leaves the
chromosome's _INDELS sub-index without a .bin, and CRISPRitz 2.9.1
refuses to search such a directory -- so before this fix a single SNP-only
chromosome made the whole search fail, reported as an SNP post-analysis failure,
where 2.8.0 had returned the right answer (no indel off-target there). Reproduced
through the real CLI on chr21 (SNPs only) + chr22.

The same refusal hit a chromosome whose indels are there but whose windows hold no
site of the PAM (Cas12a's TTTV on a chromosome with one deletion, reproduced through
the real CLI: build-index-only --vcf exits 0, complete-search exits 1, where 2.8.0
returned the right report). buildTST writes no .bin there either, so
pool_index_indels.py records that build (.no_targets) and pool_search_indels.py
skips it; an index built before has no record, and the fake-indel genome it was
built from then tells, by buildTST's own rule (crispritz_bins).

A PAM-less index standing in for the search's PAM (NNN_3 serving NRG) holds the
targets of its own PAM: an index of it without the record is told by that PAM (the
PAM line its variant index's .pam_build records), and the search's own rebuild of a
missing INDELS index for it builds it with that PAM. Built with the search's PAM file,
that index went under NRG_3_..., where the search did not look, and a chromosome with
indels but no sub-index was searched in its fake-indel genome without bulges: exit 0,
with the bulged indel off-targets lost (measured through the real CLI). A search with
bulges now refuses such a chromosome.

complete-test --chrom runs its search from a per-run view of the working directory
(complete_test._chrom_view), gone once the run ends: an ERROR printed there gave the
view's paths, so that its command could not work, and the check of an index built
before the record read a fake-indel genome the view does not hold.

A contig with an underscore in its name (chrUn_KI270302v1) has never had a sub-index:
pool_index_indels.py takes what follows the last underscore for its chromosome. With an
indel on it, that refusal stopped every search with bulges, for good, with a rebuild
command that failed too: build-index-only took the sub-index it did not find for a full
disk, and exited 1, where 2.8.0 exited 0 (both reproduced through the real CLI). The
build now names such contigs and expects no sub-index of them, and the search goes on.

Pinned here: the empty chromosome is skipped on BOTH sides (not indexed, not
searched), the merge tolerates its missing part file and only its, an all-empty
dataset still yields a well-formed targets file, a build with no target is
recorded and skipped, an index without the record is skipped only where its
fake-indel genome shows it holds no target of the index's own PAM, any other index
that is empty although its log lists indels fails with an ERROR that says how to
rebuild it (no traceback, no merged file, and no command with a blank to fill in), a
bulge search of a chromosome with indels but no sub-index fails the same way while a
search without bulges still searches its fake-indel genome, a PAM-less stand-in's
INDELS index is built with its own PAM, a search from complete-test --chrom's view (the
real one) names and reads the working directory's paths, a refusal ends the run at its
own stage, a contig with an underscore in its name is named by the build, expected by no
check, and refused by no search, and the normal path is unchanged.

Every line the pool workers of either script print reaches the log. They printed
"Searching for INDELs in"/"Search ended for INDELs in" (and "Indexing ...") unflushed, and
the pool's exit terminates the workers it finds idle: 3 to 23 of those lines went missing
from a search's log_verbose.txt, a different number each run. A sitecustomize makes every
worker meet that end, so that no line gets through by chance.

Hermetic: a stub crispritz.py on PATH stands in for the engine and, like
CRISPRitz 2.9.1, refuses an index directory with no .bin; its index-genome writes
the layout the real one writes.
"""

import ast
import gzip
import json
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import complete_test as ct
import indel_window
import test_build_index_only_atomic as atomic
from test_index_resolution import _ShippedScript
from test_indel_window import paste_both, same_argv

HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = os.path.join(HERE, "pool_search_indels.py")
INDEX = os.path.join(HERE, "pool_index_indels.py")
POSTINDEL = os.path.join(HERE, "pool_post_analisi_indel.py")
HEADER = ("#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
          "Direction\tMismatches\tBulge Size\tTotal")
LOG_HEADER = "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq"

STUB = r'''#!/usr/bin/env python3
import os, shutil, sys
with open(os.environ["STUB_LOG"], "a") as fh:
    fh.write(" ".join(sys.argv[1:]) + "\n")
if sys.argv[1] == "search":
    index_dir, out_name = sys.argv[2], sys.argv[5]
    chrom = out_name.split("_pam.txt", 1)[0][len("fake"):]
    if "-index" in sys.argv:
        bins = [f for f in os.listdir(index_dir) if f.endswith(".bin")] if os.path.isdir(index_dir) else []
        if not bins:
            sys.stderr.write(f"ERROR: no .bin index files found in '{index_dir}'.\n")
            sys.exit(1)
        row = f"X\tACGT\tACGT\tfake{chrom}\t10\t10\t+\t1\t0\t1"
    else:
        # a search of the fake-indel genome itself, which takes no bulge (its row at 20)
        row = f"X\tACGT\tACGT\tfake{chrom}\t20\t20\t+\t1\t0\t1"
    with open(out_name + ".targets.txt", "w") as fh:
        # (STUB_FIRST_LINE: a first line another user put in the part file)
        fh.write(os.environ.get("STUB_FIRST_LINE", HEADER) + "\n")
        fh.write(row + "\n")
elif sys.argv[1] == "index-genome":
    # What crispritz.py index-genome writes: genome_library/<PAM>_<bMax>_<name>, made
    # anew, the PAM read from the PAM file, a .bin for each contig of a .fa there but
    # those $STUB_NO_BIN names (buildTST writes none for a contig without a target).
    # $STUB_RECORD_BLOCKED leaves a directory where the record of no target would go.
    # The log gets the PAM file's line too ("pam-line <line>"). A genome directory that
    # is not there writes nothing: the real one says so on stdout, and exits 0.
    name, gdir, pam_file = sys.argv[2], sys.argv[3], sys.argv[4]
    if not os.path.isdir(gdir):
        sys.exit(0)
    nb = sys.argv[sys.argv.index("-bMax") + 1]
    with open(pam_file) as fh:
        seq, size = fh.read().split()[:2]
    with open(os.environ["STUB_LOG"], "a") as fh:
        fh.write(f"pam-line {seq} {size}\n")
    pam = seq[:-int(size)] if int(size) < 0 else seq[-int(size):]
    out = os.path.join("genome_library", f"{pam}_{nb}_{name}")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    for fa in sorted(os.listdir(gdir)):
        if fa.split(".")[-1] != "fa":
            continue
        with open(os.path.join(gdir, fa)) as fh:
            contig = fh.readline().strip()[1:]
        if contig not in os.environ.get("STUB_NO_BIN", "").split(","):
            open(os.path.join(out, f"{pam}_{contig}_1.bin"), "w").close()
    if os.environ.get("STUB_RECORD_BLOCKED"):
        os.makedirs(os.path.join(out, ".no_targets"))
'''.replace("HEADER", repr(HEADER))
# The PAM files the tests use (NNN: the PAM line of a PAM-less index for NRG's guide),
# a fake-indel contig that holds a target of NRG, one that holds none (but targets of
# NNN), one too short for any, and the name of the record both scripts must agree on
NRG = "NNNNNNNNNNNNNNNNNNNNNRG 3"
NNN = "N" * 23 + " 3"
TTTV = "TTTVNNNNNNNNNNNNNNNNNNNNNNN -4"
WITH_NRG = "ACGTTGCAAG" * 6  # NAG at 27 (a window of 23 + 3 bases before it)
NO_NRG = "AT" * 30  # no G, no C: no NRG on either strand
SHORT = "ACGTACGTAC"  # shorter than a window (23 + 3 bases): no target of any PAM
NO_TARGETS = ".no_targets"
# what the stub's search writes for a chromosome through its sub-index, and through its
# fake-indel genome itself (a search without bulges)
ROW_A = "X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"
ROW = "X\tACGT\tACGT\tfake{}\t{}\t{}\t+\t1\t0\t1"
# the two indels the tests' logs list, and the chromosomes of their dataset's VCFs
LOGS_AB = {"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"], "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]}
ONLY_B = {"chrA": [], "chrB": LOGS_AB["chrB"]}  # (chrA has no indel: skipped, never searched)
AB = ("chrA", "chrB")
# contigs with an underscore in their name, which pool_index_indels.py does not index
UNPLACED, RANDOM = "chrUn_KI270302v1", "chr1_KI270706v1_random"
LEFT_OUT = ("No INDELS sub-index for {}: a contig with an underscore in its name gets none, so "
            "a search does not find the indel off-targets there that need a bulge")

# Imported first by each Python process of a run that has it on PYTHONPATH (sitecustomize):
# a pool worker that has taken the sentinel that ends it cleanly, flushing its stdout, waits
# instead for the SIGTERM with which the pool's exit (terminate()) ends a worker it finds
# idle. Which of the two a worker meets is a race; here every worker meets the second.
POOL_EXIT_SHIM = r'''
import multiprocessing.pool
import signal

_worker = multiprocessing.pool.worker


def worker(*args, **kwargs):
    _worker(*args, **kwargs)
    while True:
        signal.pause()


multiprocessing.pool.worker = worker
'''


def _gnu_sed():
    try:
        return subprocess.run(["sed", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


def window_rows(chrom, left, right):
    """Two indel log rows as CRISPRitz's add-variants writes them, of fake-indel windows of
    ``left`` reference bases before the record's REF and ``right`` after it (an insertion
    at 100 and a deletion at 200)."""
    out = []
    for k, (pos, ref, alt) in enumerate(((100, "A", "AT"), (200, "ACG", "A"))):
        out.append("%s_%d-%d_%d\tS1\trs%d\t0.1\t%s_%d_%s_%s\t%d,%d\tACGT" % (
            chrom, pos - 1 - left, pos - 1 + len(ref) + right, k + 1, k, chrom, pos, ref, alt,
            100 * k, 100 * k + left + len(alt) + right))
    return out


def window_fake(chrom, left, right):
    """The fake-indel FASTA whose windows lie where window_rows' rows place them (their
    FAKEPOS fields), an N after each and Ns between them, of A and T alone (no NRG)."""
    seq = ""
    for k, (pos, ref, alt) in enumerate(((100, "A", "AT"), (200, "ACG", "A"))):
        a, b = 100 * k, 100 * k + left + len(alt) + right
        seq += "N" * (a - len(seq)) + ("AT" * 50)[:b - a] + "N"
    return f">fake{chrom}\n{seq}\n"


def _write_log(path, records):
    with gzip.open(path, "wt") as fh:
        fh.write(LOG_HEADER + "\n")
        for r in records:
            fh.write(r + "\n")


class _Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = os.path.realpath(tempfile.mkdtemp(prefix="indel_empty_"))
        self.addCleanup(shutil.rmtree, self.tmp)
        self.bin = os.path.join(self.tmp, "bin")
        os.makedirs(self.bin)
        stub = os.path.join(self.bin, "crispritz.py")
        with open(stub, "w") as fh:
            fh.write(STUB)
        os.chmod(stub, 0o755)
        self.stub_log = os.path.join(self.tmp, "calls.log")
        open(self.stub_log, "w").close()
        self.env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                        STUB_LOG=self.stub_log)

    def calls(self):
        with open(self.stub_log) as fh:
            return [ln.split() for ln in fh if ln.strip()]

    def pool_exit_shim(self):
        """The PYTHONPATH that puts POOL_EXIT_SHIM in each Python process of a run."""
        shim = os.path.join(self.tmp, "shim")
        os.makedirs(shim, exist_ok=True)
        with open(os.path.join(shim, "sitecustomize.py"), "w") as fh:
            fh.write(POOL_EXIT_SHIM)
        return os.pathsep.join(filter(None, (shim, os.environ.get("PYTHONPATH"))))

    def worker_lines(self, out, *starts):
        """The lines of ``out`` that start with one of ``starts``, sorted (a pool's order)."""
        return sorted(ln for ln in out.splitlines() if ln.startswith(starts))


@unittest.skipUnless(_gnu_sed(), "pool_search_indels.py merges with GNU `sed -i`, as the pipeline image does")
class SearchSide(_Fixture):
    """pool_search_indels.py, chromosomes taken from the indel logs (no VCF dir)."""

    def run_search(self, logs, bins, records=(), fakes=None, true_pam="NRG", pam=NRG,
                   bulges=("1", "1"), vcfs=(), pam_files=None, pam_build=None, view=None,
                   listed="ds", no_trees=(), metas=None, enriched=False):
        """bins: chrom -> whether its sub-index holds a .bin (a chromosome left out has
        none); records: the chromosomes whose sub-index holds the record of a build with
        no target; fakes: chrom -> the sequence of its fake-indel contig, written where
        build-index-only leaves it (None: a header-only one); vcfs: the chromosomes of the
        dataset's VCF folder, whence the search then takes them (else from the logs);
        pam_files: name -> line of a PAM file in the working directory's PAMs/; pam_build:
        what the variant index's .pam_build says; view: a chromosome, to search as
        complete-test --chrom does, from complete_test's own view of the working directory;
        listed: how the search's --vcf list spells the dataset; no_trees: the chromosomes
        whose sub-index holds a directory x.bin and a file named .bin, neither of them a
        tree searchTST loads; metas: chrom -> the record of what its sub-index holds
        (.index_meta.json, index_fit.write_index_meta); enriched: the dataset's enriched
        genome in place (Genomes/ref+ds)."""
        cwd = self.work = os.path.join(self.tmp, "work")
        self.true_pam = true_pam
        logdir = os.path.join(cwd, "Dictionaries", "log_indels_ds")
        os.makedirs(logdir)
        for chrom, records_ in logs.items():
            _write_log(os.path.join(logdir, f"log{chrom}.txt.gz"), records_)
        idx = self.idx = os.path.join(cwd, "genome_library", f"{true_pam}_3_ref+ds_INDELS")
        for chrom, has_bin in bins.items():
            d = os.path.join(idx, f"{true_pam}_3_fake{chrom}")
            os.makedirs(d)
            if has_bin:
                # (the head of a tree of index format 1, as the releases these tests stand
                # for built: an empty sub-index is read by the rule of its siblings' format)
                with open(os.path.join(d, f"{true_pam}_fake{chrom}_0.bin"), "wb") as fh:
                    fh.write(struct.pack("=ii", 1, 20))
            if chrom in records:
                with open(os.path.join(d, NO_TARGETS), "w") as fh:
                    fh.write("built, 0 targets\n")
            if chrom in no_trees:
                os.makedirs(os.path.join(d, "x.bin"))
                open(os.path.join(d, ".bin"), "w").close()
            if chrom in (metas or {}):
                with open(os.path.join(d, ".index_meta.json"), "w") as fh:
                    json.dump(metas[chrom], fh)
        if pam_build is not None:
            os.makedirs(os.path.join(cwd, "genome_library", f"{true_pam}_3_ref+ds"))
            with open(os.path.join(cwd, "genome_library", f"{true_pam}_3_ref+ds", ".pam_build"), "w") as fh:
                fh.write(pam_build)
        if enriched:
            os.makedirs(os.path.join(cwd, "Genomes", "ref+ds"))
            with open(os.path.join(cwd, "Genomes", "ref+ds", "chrA.enriched.fa"), "w") as fh:
                fh.write(">chrA\nACGT\n")
        for chrom, seq in (fakes or {}).items():
            d = os.path.join(cwd, "Genomes", "ref+ds_INDELS", f"fake_ds_{chrom}")
            os.makedirs(d)
            with open(os.path.join(d, f"fake{chrom}.fa"), "w") as fh:
                fh.write(f">fake{chrom}\n{seq}N\n" if seq is not None else f">fake{chrom}\n")
        for name, line in (pam_files or {}).items():
            os.makedirs(os.path.join(cwd, "PAMs"), exist_ok=True)
            with open(os.path.join(cwd, "PAMs", name), "w") as fh:
                fh.write(line + "\n")
        self.vcf_dir = os.path.join(cwd, "VCFs", "ds")
        for chrom in vcfs:
            os.makedirs(self.vcf_dir, exist_ok=True)
            with gzip.open(os.path.join(self.vcf_dir, f"ds.{chrom}.vcf.gz"), "wt") as fh:
                fh.write(f"##fileformat=VCFv4.2\n#CHROM\tPOS\n{chrom}\t5\n")
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        open(os.path.join(self.tmp, "guide.txt"), "w").close()
        with open(os.path.join(self.tmp, "pam.txt"), "w") as fh:
            fh.write(pam + "\n")

        def search(root):  # from ``root``, as the search script calls it
            args = [os.path.join(root, "Genomes", "ref"), os.path.join(root, "VCFs", listed), "ds",
                    os.path.join(self.tmp, "guide.txt"), os.path.join(self.tmp, "pam.txt"),
                    "3", "4", *bulges, out, true_pam, root + "/", "2", "-1", "0"]
            return subprocess.run([sys.executable, SEARCH, *args], env=self.env,
                                  capture_output=True, text=True)

        if view:
            res = self.from_view(view, search)
        else:
            res = search(cwd)
        merged = os.path.join(out, f"indels_ref+ds_pam.txt_guide.txt_4_{bulges[0]}_{bulges[1]}.targets.txt")
        lines = None
        if os.path.isfile(merged):
            with open(merged) as fh:
                lines = fh.read().split("\n")
        return res, lines

    def from_view(self, chrom, search):
        """``search`` from complete_test's view of the working directory for ``chrom``,
        which is gone once it has returned, as it is when complete-test --chrom ends."""
        os.makedirs(os.path.join(self.work, "Genomes", "ref"), exist_ok=True)
        os.makedirs(self.idx[: -len("_INDELS")], exist_ok=True)
        vcf_list = os.path.join(self.tmp, "list_vcf.txt")
        with open(vcf_list, "w") as fh:
            fh.write("ds\n")
        here = os.getcwd()
        os.chdir(self.work)
        try:
            with ct._chrom_view(chrom, "ref", vcf_list) as view:
                return search(view)
        finally:
            os.chdir(here)

    def searched(self, index=True):
        """The chromosomes searched through their sub-index, or (index false) through
        their fake-indel genome itself."""
        return sorted(c[4].split("_pam.txt", 1)[0][len("fake"):] for c in self.calls()
                      if c[0] == "search" and ("-index" in c) == index)

    def fake(self, chrom):
        return f"{self.work}/Genomes/ref+ds_INDELS/fake_ds_{chrom}"

    def empty(self, why):
        """What the ERROR says of a sub-index with no .bin and no record."""
        return ("holds no .bin and no record that its build found no target there (an index "
                f"built before CRISPRme 2.8.1 has none), and {why}")

    def missing(self, chrom, shown):
        """What the ERROR says of a chromosome with indels and no sub-index."""
        return (f"is not there, though {chrom} has indels ({shown}): the INDELS index is "
                "incomplete, and a search with bulges cannot do without it")

    def chain(self, pam=None, enrich=False):
        """The rebuild command an ERROR gives where the dataset's VCFs are here, with the
        PAM file ``pam`` (the search's), and the enrichment and the download again when the
        fake-indel genome is not here."""
        gone = [self.idx] + ([f"{self.work}/Genomes/ref+ds"] if enrich else [])
        again = (f"; or, for an index you downloaded, download it again: rm -rf '{self.idx}' "
                 f"&& crisprme.py download --what index --index-name {self.true_pam}_3_ref+ds "
                 f"--path '{self.work}'" if enrich else "")
        return (": rm -rf " + " ".join("'%s'" % p for p in gone) + " && crisprme.py "
                f"build-index-only --genome '{self.work}/Genomes/ref' --pam "
                f"'{pam or self.tmp + '/pam.txt'}' --bDNA 2 --bRNA 2 --vcf '{self.vcf_dir}' "
                f"--path '{self.work}'{again}")

    def refused(self, res, what, advice, chrom="chrB"):
        """An ERROR for ``chrom`` that says ``what`` of its index, and how to rebuild it
        (``advice``, after "Rebuild it"), with nothing searched there and nothing merged."""
        self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
        self.assertNotIn("Traceback", res.stderr)
        self.assertEqual(res.stderr, (f"ERROR: the indel index of {chrom} ({self.idx}/"
                                      f"{self.true_pam}_3_fake{chrom}) {what}. Rebuild it{advice}\n"))
        # every command it gives runs as given: no blank to fill in, and the PAM file a
        # rebuild names is there
        self.assertNotIn("<", res.stderr)
        for named in re.findall(r"--pam '([^']*)'", res.stderr):
            self.assertTrue(os.path.isfile(named), named)
        self.assertNotIn(chrom, self.searched() + self.searched(index=False))
        # and each command it gives, rm -rf and the build or the download after it, pasted
        # into bash and into zsh, the default shell of macOS, runs the same programs with the
        # same arguments (stubs that record them: no build, no download; wave 1c-b of 2.10.0)
        for chain in re.findall(r"(rm -rf '[^']*'(?: '[^']*')* && crisprme\.py [^;\n]*)",
                                res.stderr):
            ran = same_argv(self, chain, self.tmp)
            if ran is not None:
                self.assertEqual([r[0] for r in ran], ["rm", "crisprme.py"], (chain, ran))
                self.assertEqual(ran[0][:2], ["rm", "-rf"], ran)

    def test_snp_only_chromosome_is_skipped_not_fatal(self):
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"], "chrB": []},
                                     {"chrA": True, "chrB": False})
        self.assertEqual(res.stderr, "")  # a stage's stderr is fatal in the pipeline
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn("No indels on chrB", res.stdout)
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([ln for ln in lines[1:] if ln], ["X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"])

    def test_a_first_line_a_part_file_holds_runs_nothing(self):
        # The merge takes the header from a chromosome's part file, a file of the output
        # folder, which others may write (its group, under a umask of 002), and put it into
        # its sed command as the shell's own text, where a "$(...)" or a backquote in it ran
        # (2.10.0, the review of lane S, round 4): it goes in as the line it is
        planted = "#Bulge type\t$(touch MARK_substitution)\t`touch MARK_backquotes`\tcrRNA"
        self.env["STUB_FIRST_LINE"] = planted
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": True})
        marks = [os.path.join(r, f) for r, _d, fs in os.walk(self.tmp) for f in fs
                 if f.startswith("MARK")]
        self.assertEqual(marks, [])
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(lines[0], planted)
        self.assertEqual(len([ln for ln in lines[1:] if ln]), 2)

    def test_no_chromosome_with_indels_still_writes_a_targets_file(self):
        res, lines = self.run_search({"chrA": [], "chrB": []}, {"chrA": False, "chrB": False})
        self.assertEqual(res.stderr, "")
        self.assertEqual(self.searched(), [])
        self.assertEqual(lines, [HEADER, ""])

    def test_an_empty_index_whose_log_lists_indels_still_fails_loudly(self):
        # a broken index must never pass as "no indels": with no record of a build that
        # found no target, and no fake-indel genome to check it against (as for a
        # downloaded one: the superseded NRR_3_hg38+hg38_mega archive holds no .bin), it
        # fails, and says how to rebuild it, with the enrichment that makes that genome,
        # or to download it again
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, vcfs=AB)
        self.refused(res, self.empty(f"the fake-indel genome it was built from "
                                     f"({self.fake('chrB')}) is not here to check it against"),
                     self.chain(enrich=True))
        self.assertIsNone(lines)

    def test_a_build_with_no_target_is_skipped(self):
        # the windows of chrB's indels hold no site of the PAM: its sub-index has no .bin,
        # and the record pool_index_indels.py leaves says why
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": False}, records={"chrB"})
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stderr, "")
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn("No NRG target in the indels of chrB: skipping its indel search "
                      "(its index was built with none)\n", res.stdout)
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([ln for ln in lines[1:] if ln], ["X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"])

    def test_every_chromosome_skipped_still_writes_a_targets_file(self):
        res, lines = self.run_search({"chrA": [], "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrB": False}, records={"chrB"})
        self.assertEqual(res.stderr, "")
        self.assertEqual(self.searched(), [])
        self.assertEqual(lines, [HEADER, ""])

    def test_an_index_built_before_the_record_is_told_by_its_fake_genome(self):
        # (one an earlier release built): no record, but the fake-indel genome it was
        # built from holds no target of the PAM, by buildTST's own rule
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": False},
                                     fakes={"chrA": WITH_NRG, "chrB": NO_NRG})
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stderr, "")
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn(f"No NRG target in the indels of chrB: skipping its indel search (its "
                      f"index holds no .bin, and {self.work}/Genomes/ref+ds_INDELS/fake_ds_chrB "
                      "no target)\n", res.stdout)
        self.assertEqual([ln for ln in lines[1:] if ln], ["X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"])

    def test_an_index_built_before_the_record_whose_genome_holds_a_target_fails(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, vcfs=AB,
                                     fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NRG: "
                                     "the index is incomplete"), self.chain())
        self.assertIsNone(lines)

    def test_an_index_built_before_the_record_whose_genome_cannot_be_read_fails(self):
        # no target found in a genome read otherwise than buildTST reads it says nothing of
        # the index: skipping that chromosome could lose its indel off-targets, so it fails
        for why, seq in (("CRLF line ends", f"{NO_NRG}\r\n{NO_NRG}"),
                         ("a byte no IUPAC base is", f"{NO_NRG}?{NO_NRG}")):
            with self.subTest(why=why):
                self.tearDown_dirs()
                res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, vcfs=AB,
                                             fakes={"chrA": WITH_NRG, "chrB": seq})
                self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa cannot be read "
                                             "as buildTST reads it"), self.chain())
                self.assertIsNone(lines)

    def test_a_name_that_ends_in_bin_is_no_tree(self):
        # a directory x.bin and a file named .bin: searchTST loads neither, so the sub-index
        # holds no .bin (with no target in its genome, the chromosome is skipped)
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, no_trees={"chrB"},
                                     fakes={"chrA": WITH_NRG, "chrB": NO_NRG})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn("No NRG target in the indels of chrB", res.stdout)
        self.assertEqual([ln for ln in lines[1:] if ln], [ROW_A])

    # A PAM-less index standing in for the search's PAM (NNN_3 serving NRG; NNNNNNNN_3
    # under --prune-pam): a build of it holds the targets of its own PAM, which the
    # search's does not tell, so its sub-index is checked with the PAM line it was built
    # with (from the variant index's .pam_build, else the search's guide around NNN).

    def test_a_pamless_index_standing_in_is_told_by_its_own_pam(self):
        # its windows on chrB are too short for a target of any PAM
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                     fakes={"chrA": WITH_NRG, "chrB": SHORT})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn(f"No NNN target in the indels of chrB: skipping its indel search (its "
                      f"index holds no .bin, and {self.fake('chrB')} no target)\n", res.stdout)
        self.assertEqual([ln for ln in lines[1:] if ln], [ROW_A])

    def test_a_pamless_index_standing_in_is_not_told_by_the_searchs_pam(self):
        # chrB's windows hold no NRG site but targets of NNN, so a build of the index
        # holds a .bin there: this one is incomplete, and the ERROR names the PAM file it
        # was built with, the one .pam_build names before any other of that line
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                     fakes={"chrA": WITH_NRG, "chrB": NO_NRG}, vcfs=AB,
                                     pam_files={"0-other.txt": NNN, "20bp-NNN-NO-PAM.txt": NNN},
                                     pam_build="20bp-NNN-NO-PAM 23 3\n")
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NNN: "
                                     "the index is incomplete"),
                     self.chain(pam=f"{self.work}/PAMs/20bp-NNN-NO-PAM.txt"))
        self.assertIsNone(lines)

    def test_a_pamless_index_standing_in_names_any_pam_file_of_its_line(self):
        # (no .pam_build: the line of the search's guide around NNN, in a file of any name)
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                 fakes={"chrA": WITH_NRG, "chrB": NO_NRG}, vcfs=AB,
                                 pam_files={"mine.txt": NNN, "x.txt": NRG})
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NNN: "
                                     "the index is incomplete"),
                     self.chain(pam=f"{self.work}/PAMs/mine.txt"))

    def test_a_pamless_index_standing_in_takes_the_geometry_of_its_pam_build(self):
        # built for a 21-nt guide: a window of 24 + 3 bases, where the search's 20-nt guide
        # takes 23 + 3 (so a run of 26 bases holds a target of the one, not the other)
        for seq, refused in ((SHORT + "A" * 17, True), (SHORT + "A" * 16, False)):
            with self.subTest(run=len(seq)):
                self.tearDown_dirs()
                res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                         fakes={"chrA": WITH_NRG, "chrB": seq}, vcfs=AB,
                                         pam_files={"20bp-NNN-NO-PAM.txt": NNN,
                                                    "21bp-NNN-NO-PAM.txt": "N" * 24 + " 3"},
                                         pam_build="21bp-NNN-NO-PAM 24 3\n")
                if refused:
                    self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a "
                                                 "target of NNN: the index is incomplete"),
                                 self.chain(pam=f"{self.work}/PAMs/21bp-NNN-NO-PAM.txt"))
                else:
                    self.assertEqual((res.returncode, res.stderr), (0, ""))
                    self.assertIn("No NNN target in the indels of chrB", res.stdout)

    def tearDown_dirs(self):
        for d in ("work", "out"):
            shutil.rmtree(os.path.join(self.tmp, d), ignore_errors=True)
        open(self.stub_log, "w").close()

    def test_a_pamless_index_standing_in_with_no_pam_file_of_its_line(self):
        # (none here: the one .pam_build names holds another line) with the fake-indel
        # genome here, a search rebuilds the INDELS index from it, once it is removed
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                 fakes={"chrA": WITH_NRG, "chrB": NO_NRG}, vcfs=AB,
                                 pam_files={"20bp-NNN-NO-PAM.txt": NRG},
                                 pam_build="20bp-NNN-NO-PAM 23 3\n")
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NNN: "
                                     "the index is incomplete"),
                     f": remove the INDELS index with rm -rf '{self.idx}' and run this search "
                     f"again, which builds it anew from the fake-indel genome in "
                     f"{self.work}/Genomes/ref+ds_INDELS")

    def test_a_pamless_index_standing_in_with_no_pam_file_and_no_fake_genome(self):
        # the steps in words, with the PAM line to write into a PAM file
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN", vcfs=AB,
                                 pam_build="20bp-NNN-NO-PAM 23 3\n")
        self.refused(res, self.empty(f"the fake-indel genome it was built from "
                                     f"({self.fake('chrB')}) is not here to check it against"),
                     " with build-index-only and a PAM file of the PAM it was built with "
                     f"(20bp-NNN-NO-PAM.txt, the line '{NNN}'), which is not in "
                     f"{self.work}/PAMs or {self.tmp}: remove '{self.idx}' and "
                     f"'{self.work}/Genomes/ref+ds', write that PAM file, and run crisprme.py "
                     "build-index-only with it as --pam, and --genome "
                     f"'{self.work}/Genomes/ref' --bDNA 2 --bRNA 2 --vcf '{self.vcf_dir}' "
                     f"--path '{self.work}'; or, for an index you downloaded, download it "
                     f"again: rm -rf '{self.idx}' && crisprme.py download --what index "
                     f"--index-name NNN_3_ref+ds --path '{self.work}'")
        self.assertNotIn(" && crisprme.py build-index-only", res.stderr)

    # A search with bulges whose chromosome has indels but no sub-index: the raw-genome
    # fallback takes no bulge, so it stops, where it lost those off-targets with no error.

    def test_a_bulge_search_refuses_a_chromosome_without_its_sub_index(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True}, vcfs=AB,
                                     fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, self.missing("chrB", "its indel log lists them"), self.chain())
        self.assertIsNone(lines)

    def test_a_bulge_search_refuses_one_without_its_sub_index_or_fake_genome(self):
        # (as from an index downloaded in part)
        res, lines = self.run_search(LOGS_AB, {"chrA": True}, vcfs=AB)
        self.refused(res, self.missing("chrB", "its indel log lists them"), self.chain(enrich=True))
        self.assertIsNone(lines)

    def test_a_bulge_search_refuses_a_stand_in_without_its_sub_index(self):
        # (a fake-indel genome that shows its indels, the log being gone)
        res, _ = self.run_search({"chrA": LOGS_AB["chrA"]}, {"chrA": True}, true_pam="NNN",
                                 vcfs=AB, fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, self.missing("chrB", f"{self.fake('chrB')} holds them"),
                     f": remove the INDELS index with rm -rf '{self.idx}' and run this search "
                     f"again, which builds it anew from the fake-indel genome in "
                     f"{self.work}/Genomes/ref+ds_INDELS")

    def test_a_search_without_bulges_still_searches_the_fake_genome(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True}, fakes={"chrA": WITH_NRG, "chrB": WITH_NRG},
                                     bulges=("0", "0"))
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual((self.searched(), self.searched(index=False)), (["chrA"], ["chrB"]))
        self.assertEqual([ln for ln in lines[1:] if ln], [ROW_A, ROW.format("chrB", 20, 20)])

    def test_a_missing_sub_index_is_not_refused_where_nothing_shows_indels(self):
        # chrB has no log: with a header-only fake-indel genome (pool_index_indels.py
        # indexes none), its search goes on as before; with none, the engine says what is
        # wrong, as before
        res, _ = self.run_search({"chrA": LOGS_AB["chrA"]}, {"chrA": True}, vcfs=("chrA", "chrB"),
                                 fakes={"chrA": WITH_NRG, "chrB": None})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual((self.searched(), self.searched(index=False)), (["chrA"], ["chrB"]))
        self.tearDown_dirs()
        res, _ = self.run_search({"chrA": LOGS_AB["chrA"]}, {"chrA": True}, vcfs=("chrA", "chrB"))
        self.assertNotIn("is not there", res.stderr)
        self.assertIn(f"ERROR: no .bin index files found in '{self.idx}/NRG_3_fakechrB/'.\n", res.stderr)

    # A contig with an underscore in its name: no build gives it a sub-index, so no rebuild
    # would help. Its search goes on, through its fake-indel genome and without bulges, as
    # on every release before the refusal, which took it for an index copied in part and
    # stopped each search with bulges for good.

    def contig_logs(self, contig):
        return dict(LOGS_AB, **{contig: [f"{contig}\t\trs3\t0.3\tA>AT\t5\tACGT"]})

    def test_a_bulge_search_goes_on_for_a_contig_no_build_gives_a_sub_index(self):
        for contig in (UNPLACED, RANDOM):
            for vcfs in ((), AB + (contig,)):  # (its chromosomes from the logs, or the VCFs)
                with self.subTest(contig=contig, vcfs=vcfs):
                    self.tearDown_dirs()
                    res, lines = self.run_search(
                        self.contig_logs(contig), {"chrA": True, "chrB": True}, vcfs=vcfs,
                        fakes={"chrA": WITH_NRG, "chrB": WITH_NRG, contig: WITH_NRG})
                    self.assertEqual((res.returncode, res.stderr), (0, ""))
                    self.assertEqual((self.searched(), self.searched(index=False)),
                                     (["chrA", "chrB"], [contig]))
                    self.assertEqual(sorted(ln for ln in lines[1:] if ln),
                                     sorted([ROW_A, ROW.format("chrB", 10, 10),
                                             ROW.format(contig, 20, 20)]))

    def test_a_bulge_search_leaves_such_a_contig_without_its_fake_genome_to_the_engine(self):
        # (as on every release before the refusal: there is nothing to search it in)
        res, _ = self.run_search(self.contig_logs(UNPLACED), {"chrA": True, "chrB": True})
        self.assertNotIn("is not there", res.stderr)
        self.assertNotIn("Rebuild it", res.stderr)
        self.assertEqual(self.searched(), ["chrA", "chrB", UNPLACED])
        self.assertIn(f"ERROR: no .bin index files found in '{self.idx}/NRG_3_fake{UNPLACED}/'.\n",
                      res.stderr)

    def test_a_chromosome_beside_such_a_contig_is_still_refused(self):
        # chrB is one a build gives a sub-index: this index was copied in part, and the ERROR
        # names it alone
        res, lines = self.run_search(self.contig_logs(UNPLACED), {"chrA": True},
                                     vcfs=AB + (UNPLACED,),
                                     fakes={"chrA": WITH_NRG, "chrB": WITH_NRG, UNPLACED: WITH_NRG})
        self.refused(res, self.missing("chrB", "its indel log lists them"), self.chain())
        self.assertIsNone(lines)

    # No VCFs of the dataset here (as for an index downloaded, or copied from where it was
    # built): build-index-only stops at its --vcf, after the rm -rf it was chained to.

    def words(self, search="this search again,"):
        """What an ERROR says to do with the fake-indel genome here."""
        return (f": remove the INDELS index with rm -rf '{self.idx}' and run {search} which "
                f"builds it anew from the fake-indel genome in {self.work}/Genomes/ref+ds_INDELS")

    def download(self):
        """What an ERROR says to do with neither those VCFs nor that genome here. It removes
        only the INDELS index, the one thing in the way of the download, so a download that
        fails (an index that is not on the hub) takes nothing more."""
        return (f": for an index you downloaded, download it again: rm -rf '{self.idx}' && "
                f"crisprme.py download --what index --index-name {self.true_pam}_3_ref+ds "
                f"--path '{self.work}'; for one you built, put its VCFs back in '{self.vcf_dir}' "
                "and run this search again, which then gives the command that rebuilds it")

    def test_without_the_vcfs_a_search_rebuilds_from_the_fake_genome(self):
        for bins, what in (({"chrA": True}, self.missing("chrB", "its indel log lists them")),
                           ({"chrA": True, "chrB": False},
                            self.empty("{}/fakechrB.fa holds a target of NRG: the index is incomplete"))):
            with self.subTest(what=what):
                self.tearDown_dirs()
                res, lines = self.run_search(LOGS_AB, bins, fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
                self.refused(res, what.format(self.fake("chrB")), self.words())
                self.assertNotIn("build-index-only", res.stderr)
                self.assertNotIn("Genomes/ref+ds'", res.stderr)
                self.assertIsNone(lines)

    def test_without_the_vcfs_or_the_fake_genome_there_is_the_download(self):
        for bins, what in (({"chrA": True}, self.missing("chrB", "its indel log lists them")),
                           ({"chrA": True, "chrB": False},
                            self.empty("the fake-indel genome it was built from ({}) is not here "
                                       "to check it against"))):
            with self.subTest(what=what):
                self.tearDown_dirs()
                res, lines = self.run_search(LOGS_AB, bins)
                self.refused(res, what.format(self.fake("chrB")), self.download())
                self.assertNotIn("build-index-only", res.stderr)
                self.assertNotIn("Genomes/ref+ds'", res.stderr)
                self.assertIsNone(lines)

    def test_a_list_that_spells_the_folder_is_read_as_the_dataset(self):
        # (--vcf list with VCFs/<dataset>, which the search's chromosomes survive: it takes
        # them from the logs) the VCFs are found, and the command names their folder ...
        for view in (None, "chrB"):
            with self.subTest(view=view):
                self.tearDown_dirs()
                res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, vcfs=AB,
                                         listed="VCFs/ds", view=view)
                self.refused(res, self.empty(f"the fake-indel genome it was built from "
                                             f"({self.fake('chrB')}) is not here to check it "
                                             "against"), self.chain(enrich=True))
                # ... and without them, it is that folder they go back in
                self.tearDown_dirs()
                res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, listed="VCFs/ds",
                                         view=view)
                self.refused(res, self.empty(f"the fake-indel genome it was built from "
                                             f"({self.fake('chrB')}) is not here to check it "
                                             "against"), self.download())

    def test_a_pamless_index_without_the_vcfs_needs_no_pam_file(self):
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                 fakes={"chrA": WITH_NRG, "chrB": NO_NRG})
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NNN: "
                                     "the index is incomplete"), self.words())
        self.tearDown_dirs()
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN")
        self.refused(res, self.empty(f"the fake-indel genome it was built from ({self.fake('chrB')}) "
                                     "is not here to check it against"), self.download())

    def test_a_vcf_folder_without_a_vcf_is_no_vcfs(self):
        # (a search refuses that folder at its input check, but for complete-test --chrom,
        # whose view has no VCFs/<dataset>) add-variants would enrich nothing: the command
        # then removes the enriched genome, exits 0, and leaves no INDELS index
        os.makedirs(os.path.join(self.tmp, "work", "VCFs", "ds"))
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, view="chrB")
        self.refused(res, self.empty(f"the fake-indel genome it was built from ({self.fake('chrB')}) "
                                     "is not here to check it against"), self.download())

    # complete-test --chrom searches from a view of the working directory, gone once the
    # run ends: an ERROR gives the working directory's paths, or its command cannot work,
    # and what tells of an index is read there, for the view leaves it out.

    def test_a_view_refuses_with_the_working_directorys_paths(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True}, vcfs=AB, view="chrB",
                                     fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, self.missing("chrB", "its indel log lists them"), self.chain())
        self.assertNotIn(ct.CHROMVIEWPREFIX, res.stderr)
        self.assertIsNone(lines)

    def test_a_view_refuses_with_the_enrichment_when_the_fake_genome_is_not_there(self):
        for bins, what in (({"chrA": True}, self.missing("chrB", "its indel log lists them")),
                           ({"chrA": True, "chrB": False},
                            "holds no .bin and no record that its build found no target there "
                            "(an index built before CRISPRme 2.8.1 has none), and the fake-indel "
                            "genome it was built from ({}) is not here to check it against")):
            with self.subTest(what=what):
                self.tearDown_dirs()
                res, _ = self.run_search(LOGS_AB, bins, vcfs=AB, view="chrB")
                self.refused(res, what.format(self.fake("chrB")), self.chain(enrich=True))
                self.assertNotIn(ct.CHROMVIEWPREFIX, res.stderr)

    def test_a_view_without_the_vcfs_sends_the_user_to_the_run_that_builds_the_index(self):
        # (its INDELS index, once removed, is not built by a run with --chrom, or goes
        # with the view)
        res, _ = self.run_search(LOGS_AB, {"chrA": True}, fakes={"chrA": WITH_NRG, "chrB": WITH_NRG},
                                 view="chrB")
        self.refused(res, self.missing("chrB", "its indel log lists them"),
                     self.words("complete-test without --chrom, or a complete-search of this index,"))
        self.assertNotIn(ct.CHROMVIEWPREFIX, res.stderr)

    def test_a_view_with_no_pam_file_of_a_stand_in_sends_the_user_there_too(self):
        res, _ = self.run_search(LOGS_AB, {"chrA": True}, true_pam="NNN", vcfs=AB, view="chrB",
                                 fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, self.missing("chrB", "its indel log lists them"),
                     self.words("complete-test without --chrom, or a complete-search of this index,"))

    def test_a_view_without_the_vcfs_or_the_fake_genome_gives_the_download(self):
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, view="chrB")
        self.refused(res, self.empty(f"the fake-indel genome it was built from ({self.fake('chrB')}) "
                                     "is not here to check it against"), self.download())
        self.assertNotIn(ct.CHROMVIEWPREFIX, res.stderr)

    def test_a_view_reads_the_fake_genome_in_the_working_directory(self):
        # an index built before the record whose genome holds no target: as outside a view
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, view="chrB",
                                     fakes={"chrA": WITH_NRG, "chrB": NO_NRG})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.searched(), [])
        self.assertIn(f"No NRG target in the indels of chrB: skipping its indel search (its "
                      f"index holds no .bin, and {self.fake('chrB')} no target)\n", res.stdout)
        self.assertNotIn(ct.CHROMVIEWPREFIX, res.stdout)
        self.assertEqual(lines, [HEADER, ""])

    def test_a_view_refuses_an_index_its_fake_genome_shows_incomplete(self):
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, vcfs=AB,
                                 fakes={"chrA": WITH_NRG, "chrB": WITH_NRG}, view="chrB")
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NRG: "
                                     "the index is incomplete"), self.chain())

    def test_a_view_reads_the_record_of_a_build_with_no_target(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, records={"chrB"}, view="chrB")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertIn("No NRG target in the indels of chrB: skipping its indel search (its index "
                      "was built with none)\n", res.stdout)
        self.assertEqual(lines, [HEADER, ""])

    def test_a_view_reads_the_pam_build_in_the_working_directory(self):
        # (the view leaves out that sidecar) the guide of a 21-nt build, and the PAM file it
        # names, in the working directory's PAMs/
        res, _ = self.run_search(LOGS_AB, {"chrA": True, "chrB": False}, true_pam="NNN",
                                 vcfs=AB, view="chrB",
                                 fakes={"chrA": WITH_NRG, "chrB": SHORT + "A" * 17},
                                 pam_files={"20bp-NNN-NO-PAM.txt": NNN,
                                            "21bp-NNN-NO-PAM.txt": "N" * 24 + " 3"},
                                 pam_build="21bp-NNN-NO-PAM 24 3\n")
        self.refused(res, self.empty(f"{self.fake('chrB')}/fakechrB.fa holds a target of NNN: "
                                     "the index is incomplete"),
                     self.chain(pam=f"{self.work}/PAMs/21bp-NNN-NO-PAM.txt"))

    def test_a_view_searches_a_whole_sub_index(self):
        res, lines = self.run_search(LOGS_AB, {"chrA": True, "chrB": True}, view="chrB")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.searched(), ["chrB"])
        self.assertEqual([ln for ln in lines[1:] if ln], [ROW.format("chrB", 10, 10)])

    def test_a_record_is_not_read_where_the_index_holds_a_bin(self):
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": True}, records={"chrB"})
        self.assertEqual(res.stderr, "")
        self.assertEqual(self.searched(), ["chrA", "chrB"])
        self.assertEqual([ln.split("\t")[3] for ln in lines[1:] if ln], ["fakechrA", "fakechrB"])

    def test_normal_path_is_unchanged(self):
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": True})
        self.assertEqual(res.stderr, "")
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([ln.split("\t")[3] for ln in lines[1:] if ln], ["fakechrA", "fakechrB"])
        self.assertNotIn("target in the indels", res.stdout)

    def test_a_sub_index_of_other_windows_than_its_log_is_refused(self):
        """A sub-index built from fake-indel windows of another width than its chromosome's
        indel log records (an index kept beside an enrichment made again with CRISPRitz
        2.10.0): the log would place every indel row of it wrong. One of index format 1,
        which only CRISPRitz 2.9.1 and earlier build (25/27), beside a log of 30/30, and
        one that records 30/30 beside a log of 25/27, are refused, saying both widths; one
        that records the log's own windows is searched."""
        for metas, rows, built, logged in (
                (None, window_rows("chrA", 30, 30), "25/27", "30/30"),
                ({"chrA": {"index_format": 2, "indel_window": [30, 30]}},
                 window_rows("chrA", 25, 27), "30/30", "25/27")):
            with self.subTest(built=built):
                shutil.rmtree(os.path.join(self.tmp, "work"), ignore_errors=True)
                shutil.rmtree(os.path.join(self.tmp, "out"), ignore_errors=True)
                open(self.stub_log, "w").close()
                res, lines = self.run_search({"chrA": rows}, {"chrA": True}, metas=metas)
                self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
                self.assertIsNone(lines)
                self.assertIn("ERROR: the indel index of chrA (%s) was built from fake-indel "
                              "windows of %s reference bases" % (
                                  os.path.join(self.idx, "NRG_3_fakechrA"), built), res.stderr)
                self.assertIn("holds windows of %s: the index and the log come from two "
                              "enrichments" % logged, res.stderr)
                self.assertEqual(self.searched(), [])
        shutil.rmtree(os.path.join(self.tmp, "work"))
        shutil.rmtree(os.path.join(self.tmp, "out"))
        res, lines = self.run_search({"chrA": window_rows("chrA", 30, 30)}, {"chrA": True},
                                     metas={"chrA": {"index_format": 2, "indel_window": [30, 30]}})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.searched(), ["chrA"])

    def test_the_rebuild_of_a_sub_index_of_other_windows_names_its_variant_index(self):
        """W1C-B2R2-2, the review of wave 1c-a: the refusal of a sub-index of other windows
        than its log named the INDELS index alone for its rebuild. A build makes the variant
        index of its name with it, from the same enrichment, and that one holds no window to
        tell: after the kit of 2.10.0 had replaced 2.8.1's beside the mega's SaCas9 indexes a
        search had built from the old one, the rebuild it gave ("Variant index already
        present") built the INDELS index alone, and the next search read the old kit's variant
        index with no word, 16 targets of chr22 lost. Where the variant index is there, each
        way the refusal gives removes it too: the build-index-only command (the VCFs and the
        PAM file here), this search again (the enrichment in place, with no VCFs), or, with no
        enriched genome either, the download, or the enriched genome put back first. Where it
        is not, the refusal reads as before."""
        rows, meta = window_rows("chrA", 30, 30), {"chrA": {"index_format": 2,
                                                          "indel_window": [25, 27]}}
        fake = {"chrA": window_fake("chrA", 30, 30)}
        variant = os.path.join(self.tmp, "work", "genome_library", "NRG_3_ref+ds")
        enriched = os.path.join(self.tmp, "work", "Genomes", "ref+ds")
        for case, kw, ends in (
                ("VCFs here", dict(vcfs=["chrA"]),
                 lambda: (": rm -rf '%s' '%s' && crisprme.py build-index-only --genome "
                          "'%s/Genomes/ref' --pam '%s/pam.txt' --bDNA 2 --bRNA 2 --vcf '%s' "
                          "--path '%s'\n" % (self.idx, variant, self.work, self.tmp,
                                              self.vcf_dir, self.work))),
                ("the enrichment in place", dict(enriched=True),
                 lambda: (": remove both with rm -rf '%s' '%s' and run this search again, "
                          "which builds them anew from the enrichment in place, %s and the "
                          "fake-indel genome in %s/Genomes/ref+ds_INDELS\n" % (
                              self.idx, variant, enriched, self.work))),
                ("neither", {},
                 lambda: (": for an index you downloaded, download it again: rm -rf '%s' && "
                          "crisprme.py download --what index --index-name NRG_3_ref+ds --path "
                          "'%s'; for one you built, put its enriched genome back in '%s' and run "
                          "this search again, which then says how to rebuild them\n" % (
                              self.idx, self.work, enriched)))):
            with self.subTest(case=case):
                shutil.rmtree(os.path.join(self.tmp, "work"), ignore_errors=True)
                shutil.rmtree(os.path.join(self.tmp, "out"), ignore_errors=True)
                open(self.stub_log, "w").close()
                res, lines = self.run_search({"chrA": rows}, {"chrA": True}, metas=meta,
                                             fakes=fake, pam_build="pam.txt 23 3\n", **kw)
                self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
                self.assertIsNone(lines)
                self.assertIn("The variant index of its name (%s) was built with it, from the "
                              "same enrichment, which it holds no window to tell: it goes too, "
                              "or it would be searched with the dataset's dictionaries, its "
                              "alleles another enrichment's. Rebuild them" % variant, res.stderr)
                self.assertTrue(res.stderr.endswith(ends()), res.stderr)
                self.assertEqual(self.searched(), [])
        # with no variant index of its name there, the INDELS index alone, as before
        shutil.rmtree(os.path.join(self.tmp, "work"))
        shutil.rmtree(os.path.join(self.tmp, "out"))
        res, lines = self.run_search({"chrA": rows}, {"chrA": True}, metas=meta, fakes=fake,
                                     vcfs=["chrA"])
        self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
        self.assertNotIn("variant index", res.stderr)
        self.assertTrue(res.stderr.endswith(
            "download. Rebuild it: rm -rf '%s' && crisprme.py build-index-only --genome "
            "'%s/Genomes/ref' --pam '%s/pam.txt' --bDNA 2 --bRNA 2 --vcf '%s' --path '%s'\n" % (
                self.idx, self.work, self.tmp, self.vcf_dir, self.work)), res.stderr)

    def test_a_search_whose_targets_need_more_than_the_windows_hold_says_so(self):
        """A 27-nt SaCas9 site (21-nt guide and NNGRRT) with a DNA bulge needs 28 bases,
        and a window of 25/27 holds up to 27 around an indel (26 at an MNV): the search
        says so on stdout (the run's log_verbose.txt). On windows of 30/30 (up to 31) it
        says nothing."""
        sa = "N" * 21 + "NNGRRT 6"
        res, _lines = self.run_search({"chrA": window_rows("chrA", 25, 27)}, {"chrA": True},
                                      true_pam="NNGRRT", pam=sa)
        self.assertEqual(res.stderr, "")
        self.assertIn("WARNING [indel search]: the fake-indel windows of ds hold 25 reference "
                      "bases before each indel and 27 after it, so a target that holds an "
                      "indel there holds up to 27 bases (26 at a record whose REF and ALT "
                      "differ at their first base, an MNV or complex record), guide, PAM and "
                      "DNA bulges counted; this search's targets need up to 28 (a 27-nt guide "
                      "and PAM with 1 DNA bulge)", res.stdout)
        self.assertIn("enrich ds again with it (build-index-only) to search them.", res.stdout)
        # a sites-only dataset (its VCFs name no sample, as the mega's) is not sent to
        # build-index-only, whose enrichment of them holds no window (W-R3-F1)
        shutil.rmtree(os.path.join(self.tmp, "work"))
        shutil.rmtree(os.path.join(self.tmp, "out"))
        res, _lines = self.run_search({"chrA": window_rows("chrA", 25, 27)}, {"chrA": True},
                                      true_pam="NNGRRT", pam=sa, vcfs=["chrA"])
        self.assertEqual(res.stderr, "")
        self.assertIn("To search them: ds is sites-only: its VCFs name no sample", res.stdout)
        self.assertIn("install its archive and enriched-genome kit of CRISPRme 2.10.0",
                      res.stdout)
        self.assertNotIn("enrich ds again", res.stdout)
        # and what goes first, as build-index-only's warning says: all that came from its kit
        # in place, each part by absolute path, its indel logs among it, and its indexes by a
        # find of their name patterns, which take those built after the note too, and the
        # three the recipe writes anew by absolute path too (the note named the three alone,
        # relative, and section 12's commands then left the kit in place: W1C-B1; then the
        # indexes in place alone, and its own flow left those its search built next:
        # W1C-B2R2-2; then the patterns given to rm -rf, which zsh, the default shell of
        # macOS, expanded itself, stopping the whole command where one matched nothing:
        # wave 1c-b of 2.10.0). Pasted into bash or into zsh, it does the same.
        self.assertIn("Either replaces its kit in place and all that was built from it",
                      res.stdout)
        gone = re.search(r"which go first \((.*?); and any copy", res.stdout)
        self.assertIsNotNone(gone, res.stdout)
        lib = os.path.dirname(self.idx)
        self.assertEqual(gone.group(1).split("; "), [
            "rm -rf '%s'" % os.path.join(self.work, "Dictionaries", "log_indels_ds")]
            + indel_window.index_removals(lib, "ref", "ds"))
        paste_both(self, gone.group(1), self.tmp)
        self.assertFalse(os.path.lexists(os.path.join(self.work, "Dictionaries", "log_indels_ds")))
        self.assertIn("(rm -rf '%s' '%s' '%s')" % (
            os.path.join(self.work, "VCFs", "ds"),
            os.path.join(self.work, "Dictionaries", "registry_ds"),
            os.path.join(self.work, "Dictionaries", "indel_af_ds")), res.stdout)
        # the same from complete-test --chrom's view of the working directory, which the run
        # removes when it ends: every path named is the working directory's, none the view's
        # (the dataset sites-only by its registry there, the view's VCFs/ being empty)
        shutil.rmtree(os.path.join(self.tmp, "work"))
        shutil.rmtree(os.path.join(self.tmp, "out"))
        registry = os.path.join(self.tmp, "work", "Dictionaries", "registry_ds")
        os.makedirs(registry)
        with open(os.path.join(registry, "reg_chrA.idx"), "w") as fh:
            fh.write('{"n_records": 1, "aggregation": "info_af", "databases": {}}')
        res, _lines = self.run_search({"chrA": window_rows("chrA", 25, 27)}, {"chrA": True},
                                      true_pam="NNGRRT", pam=sa, view="chrA")
        self.assertIn("To search them: ds is sites-only", res.stdout, res.stderr)
        gone = re.search(r"which go first \((.*?); and any copy", res.stdout)
        self.assertIsNotNone(gone, res.stdout)
        commands = gone.group(1).split("; ")
        self.assertIn(indel_window.index_removals(
            os.path.join(self.work, "genome_library"), "ref", "ds")[1], commands)
        # (each path the rm names, and the folder of each find)
        named = [w for c in commands if c.startswith("rm -rf ") for w in shlex.split(c)[2:]] + [
            shlex.split(c)[1] for c in commands if c.startswith("find ")]
        self.assertEqual(len(named), len(commands) - 1 + len(shlex.split(commands[0])[2:]))
        self.assertIn(os.path.join(self.work, "Dictionaries", "log_indels_ds"), named)
        self.assertIn(os.path.join(self.work, "genome_library"), named)
        self.assertNotIn(os.path.join(self.work, "Dictionaries", "registry_ds"), named)
        for path in named + re.findall(r"rm -rf '([^']*)'", res.stdout):
            self.assertTrue(path.startswith(self.work + os.sep), path)
            self.assertNotIn(ct.CHROMVIEWPREFIX, path)
        shutil.rmtree(os.path.join(self.tmp, "work"))
        shutil.rmtree(os.path.join(self.tmp, "out"))
        res, _lines = self.run_search({"chrA": window_rows("chrA", 30, 30)}, {"chrA": True},
                                      true_pam="NNGRRT", pam=sa,
                                      metas={"chrA": {"index_format": 2,
                                                      "indel_window": [30, 30]}})
        self.assertEqual(res.stderr, "")
        self.assertNotIn("fake-indel windows", res.stdout)

    def test_every_line_of_a_worker_reaches_the_log(self):
        chroms = [f"chr{c}" for c in "ABCDEF"]
        self.env["PYTHONPATH"] = self.pool_exit_shim()
        res, lines = self.run_search({c: [f"{c}\t\trs{i}\t0.1\tA>AT\t5\tACGT"]
                                      for i, c in enumerate(chroms)}, dict.fromkeys(chroms, True))
        self.assertEqual(res.stderr, "")
        self.assertEqual(sorted(ln.split("\t")[3] for ln in lines[1:] if ln),
                         [f"fake{c}" for c in chroms])
        self.assertEqual(self.worker_lines(res.stdout, "Search"), sorted(
            [f"Searching for INDELs in {c}" for c in chroms]
            + [f"Search ended for INDELs in {c}" for c in chroms]))


@unittest.skipUnless(_gnu_sed(), "pool_search_indels.py merges with GNU `sed -i`, as the pipeline image does")
class ItsCommandsWhereAPathHoldsAQuote(_Fixture):
    """IR2A-R2-4, the review of wave 2a of CRISPRme 2.10.0: the commands a refusal of
    pool_search_indels.py gives (_rebuild) quoted each path by hand, '<path>', and a path that
    holds a quote, a working directory "Luca's data", cut such a command there, pasted: the
    shell waited for a quote to close the rest, or ran it on other paths. Each path is quoted
    as a shell reads it back (indel_window.quote), and each command runs on those paths in
    bash and in zsh. (SearchSide's fixture, under a folder of that name; chrA has no indel, so
    that nothing is searched: the search's own commands take the paths unquoted.)"""

    run_search, from_view, searched = SearchSide.run_search, SearchSide.from_view, SearchSide.searched

    def setUp(self):
        super().setUp()
        self.top = self.tmp
        self.tmp = os.path.join(self.top, "Luca's data")
        os.makedirs(self.tmp)

    def advice(self, res):
        """What the refusal says to do: its ERROR, the one line on stderr, after "Rebuild it"
        (or "them")."""
        self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
        self.assertEqual(res.stderr.count("\n"), 1, res.stderr)
        self.assertEqual(self.searched() + self.searched(index=False), [])
        m = re.search(r"\. Rebuild (?:it|them)(.*)\n\Z", res.stderr, re.S)
        self.assertIsNotNone(m, res.stderr)
        return m.group(1)

    def download(self):
        return [["rm", "-rf", self.idx], ["crisprme.py", "download", "--what", "index",
                                          "--index-name", "NRG_3_ref+ds", "--path", self.work]]

    def test_the_rebuild_and_the_download_it_gives_run_on_those_paths(self):
        res, _ = self.run_search(ONLY_B, {"chrA": False, "chrB": False}, vcfs=AB)
        build, download = self.advice(res).split(
            "; or, for an index you downloaded, download it again: ")
        self.assertTrue(build.startswith(": rm -rf "), build)
        self.assertEqual(same_argv(self, build[2:], self.top), [
            ["rm", "-rf", self.idx, os.path.join(self.work, "Genomes", "ref+ds")],
            ["crisprme.py", "build-index-only", "--genome",
             os.path.join(self.work, "Genomes", "ref"), "--pam",
             os.path.join(self.tmp, "pam.txt"), "--bDNA", "2", "--bRNA", "2", "--vcf",
             self.vcf_dir, "--path", self.work]])
        self.assertEqual(same_argv(self, download, self.top), self.download())

    def test_the_rm_it_gives_removes_that_index_pasted_into_bash_or_zsh(self):
        # (the fake-indel genome here, no VCFs: the search builds it anew once it is gone)
        res, _ = self.run_search(ONLY_B, {"chrA": False}, fakes={"chrB": WITH_NRG})
        m = re.fullmatch(r": remove the INDELS index with (rm -rf .*) and run this search "
                         r"again, which builds it anew from the fake-indel genome in (.*)",
                         self.advice(res))
        self.assertIsNotNone(m, res.stderr)
        self.assertEqual(m.group(2), os.path.join(self.work, "Genomes", "ref+ds_INDELS"))
        paste_both(self, m.group(1), self.top)
        self.assertFalse(os.path.exists(self.idx))
        self.assertTrue(os.path.isdir(m.group(2)))

    def test_the_rm_of_both_removes_the_variant_index_too(self):
        # (an INDELS index of other windows than its log, the variant index of its name
        # there, and the enrichment in place)
        variant = os.path.join(self.tmp, "work", "genome_library", "NRG_3_ref+ds")
        res, _ = self.run_search({"chrA": window_rows("chrA", 30, 30)}, {"chrA": True},
                                 metas={"chrA": {"index_format": 2, "indel_window": [25, 27]}},
                                 fakes={"chrA": window_fake("chrA", 30, 30)},
                                 pam_build="pam.txt 23 3\n", enriched=True)
        m = re.match(r": remove both with (rm -rf .*) and run this search again, which builds "
                     r"them anew from the enrichment in place", self.advice(res))
        self.assertIsNotNone(m, res.stderr)
        paste_both(self, m.group(1), self.top)
        self.assertFalse(os.path.exists(self.idx) or os.path.exists(variant))
        self.assertTrue(os.path.isdir(os.path.join(self.work, "Genomes", "ref+ds")))

    def test_the_steps_in_words_name_those_paths_as_a_shell_reads_them(self):
        # (a PAM-less index standing in, with no PAM file of its PAM, and no fake-indel
        # genome): the arguments it gives for build-index-only, and the download
        res, _ = self.run_search(ONLY_B, {"chrA": False, "chrB": False}, true_pam="NNN",
                                 vcfs=AB, pam_build="20bp-NNN-NO-PAM 23 3\n")
        words, download = self.advice(res).split(
            "; or, for an index you downloaded, download it again: ")
        removed, args = re.fullmatch(r" with build-index-only and a PAM file of the PAM it was "
                                     r"built with \(.*\), which is not in .*: remove (.*), write "
                                     r"that PAM file, and run crisprme\.py build-index-only with "
                                     r"it as --pam, and (--genome .*)", words).groups()
        self.assertEqual(shlex.split(removed.replace(" and ", " ")),
                         [self.idx, os.path.join(self.work, "Genomes", "ref+ds")])
        self.assertEqual(shlex.split(args), [
            "--genome", os.path.join(self.work, "Genomes", "ref"), "--bDNA", "2", "--bRNA", "2",
            "--vcf", self.vcf_dir, "--path", self.work])
        self.assertEqual(same_argv(self, download, self.top),
                         [[w.replace("NRG_3", "NNN_3") for w in c] for c in self.download()])


class IndexSide(_Fixture):
    """pool_index_indels.py must not hand an empty fake genome to buildTST, and records
    a build that found no target."""

    def run_index(self, fastas, pam=NRG, true_pam=None, pam_build=None, dataset="ds", **env):
        """``true_pam``: the PAM of the variant index, when it is not that of ``pam``
        (a PAM-less index standing in for it, as a search passes it on), with what its
        .pam_build says."""
        self.dataset = dataset
        folder = os.path.join(self.tmp, "Genomes", f"ref+{dataset}_INDELS")
        for chrom, text in fastas.items():
            d = os.path.join(folder, f"fake_{dataset}_{chrom}")
            os.makedirs(d)
            if text is not None:
                with open(os.path.join(d, f"fake{chrom}.fa"),
                          "wb" if isinstance(text, bytes) else "w") as fh:
                    fh.write(text)
        with open(os.path.join(self.tmp, "pam.txt"), "w") as fh:
            fh.write(pam + "\n")
        true_pam = true_pam or (pam.split()[0][:4] if pam == TTTV else pam.split()[0][-3:])
        self.true_pam = true_pam
        if pam_build is not None:
            os.makedirs(os.path.join(self.tmp, "genome_library", f"{true_pam}_3_ref+ds"))
            with open(os.path.join(self.tmp, "genome_library", f"{true_pam}_3_ref+ds", ".pam_build"), "w") as fh:
                fh.write(pam_build)
        res = subprocess.run([sys.executable, INDEX, folder + "/", "pam.txt", true_pam, "ref", dataset, "3", "2"],
                             cwd=self.tmp, env=dict(self.env, **env), capture_output=True, text=True)
        indexed = sorted(c[2].rsplit("_", 1)[-1] for c in self.calls() if c[0] == "index-genome")
        return res, indexed

    def built_with(self):
        """The PAM file each index-genome call got, and the line it read there."""
        calls = self.calls()
        return ([c[3] for c in calls if c[0] == "index-genome"],
                sorted({" ".join(c[1:]) for c in calls if c[0] == "pam-line"}))

    def sub(self, chrom):
        return os.path.join(self.top(), f"{self.true_pam}_3_fake{chrom}")

    def top(self):
        return os.path.join(self.tmp, "genome_library", f"{self.true_pam}_3_ref+{self.dataset}_INDELS")

    def recorded(self):
        return sorted(d.rsplit("_fake", 1)[-1] for d in os.listdir(self.top())
                      if os.path.isfile(os.path.join(self.top(), d, NO_TARGETS)))

    def test_header_only_fake_genome_is_not_indexed(self):
        res, indexed = self.run_index({"chrA": ">fakechrA\nACGTNNACGT\n", "chrB": ">fakechrB\n"})
        self.assertEqual(res.stderr, "")
        self.assertEqual(indexed, ["chrA"])
        self.assertIn("No indels on chrB: nothing to index", res.stdout)

    def test_a_dir_without_any_fasta_is_not_mistaken_for_empty(self):
        _, indexed = self.run_index({"chrA": ">fakechrA\nACGT\n", "chrC": None})
        self.assertEqual(indexed, ["chrA", "chrC"])
        self.assertEqual(self.recorded(), [])  # (nor for a genome with no target)

    def test_a_build_with_no_target_is_recorded(self):
        # chrB's indels are there, but no window of one holds a site of the PAM: buildTST
        # writes no .bin, and the build says so in that sub-index
        for pam, no_target in ((NRG, NO_NRG), (TTTV, NO_NRG + "TTTA")):  # (one at the end)
            with self.subTest(pam=pam):
                self.tearDown_dirs()
                res, indexed = self.run_index(
                    {"chrA": f">fakechrA\n{WITH_NRG}TTTA{WITH_NRG}\n",
                     "chrB": f">fakechrB\n{no_target}\n"}, pam=pam, STUB_NO_BIN="fakechrB")
                self.assertEqual(res.returncode, 0, res.stdout + res.stderr)
                self.assertEqual(res.stderr, "")
                self.assertEqual(indexed, ["chrA", "chrB"])
                self.assertEqual(self.recorded(), ["chrB"])
                # (with the record of what the sub-index holds, its index format among it,
                # which a sub-index with no .bin records nowhere else: test_index_format2)
                self.assertEqual(sorted(os.listdir(self.sub("chrB"))),
                                 [".index_meta.json", NO_TARGETS])
                with open(os.path.join(self.sub("chrB"), NO_TARGETS)) as fh:
                    self.assertEqual(fh.read(), (
                        "built, 0 targets: no window of fake_ds_chrB holds a target of pam.txt "
                        "with -bMax 3, so buildTST wrote no .bin\n"))
                self.assertIn(
                    f"No {self.true_pam} target in the indels of chrB: its index holds no .bin "
                    f"(genome_library/{self.true_pam}_3_ref+ds_INDELS/{self.true_pam}_3_fakechrB/"
                    f"{NO_TARGETS} records it)\n", res.stdout)

    def tearDown_dirs(self):
        for d in ("Genomes", "genome_library"):
            shutil.rmtree(os.path.join(self.tmp, d), ignore_errors=True)
        open(self.stub_log, "w").close()

    def test_no_record_where_the_engine_left_out_the_bin_of_a_target(self):
        # (as buildTST does with the disk's inode table full): build-index-only fails on
        # that, and a search refuses the sub-index
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n",
                                 "chrB": f">fakechrB\n{WITH_NRG}\n"}, STUB_NO_BIN="fakechrB")
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(self.recorded(), [])
        self.assertEqual(os.listdir(self.sub("chrB")), [".index_meta.json"])  # (its format)
        self.assertNotIn("target in the indels", res.stdout)

    def test_no_record_for_a_genome_read_otherwise_than_buildtst_reads_it(self):
        res, _ = self.run_index({"chrB": f">fakechrB\r\n{NO_NRG}\r\n"}, STUB_NO_BIN="fakechrB")
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(self.recorded(), [])

    def test_a_record_that_cannot_be_written_fails_the_build(self):
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n",
                                 "chrB": f">fakechrB\n{NO_NRG}\n"},
                                STUB_NO_BIN="fakechrB", STUB_RECORD_BLOCKED="1")
        self.assertEqual(res.returncode, 1)
        self.assertEqual(res.stderr, "")  # (the exit status says it, as for a failed build)
        self.assertIn(f"Indexing INDELs in chrB: could not write {self.sub('chrB')[len(self.tmp) + 1:]}"
                      f"/{NO_TARGETS} (", res.stdout)
        self.assertIn("Indexing failed for INDELs in chrB\n", res.stdout)

    def test_each_sub_index_records_the_windows_its_log_records(self):
        """Each sub-index records the fake-indel windows its chromosome's indel log
        (Dictionaries/log_indels_<dataset> of the working directory that holds Genomes/)
        records, (left, right) reference bases around an indel, which a search reads beside
        the log (pool_search_indels' _window_mismatch); none where the log tells none."""
        logs = os.path.join(self.tmp, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        _write_log(os.path.join(logs, "logchrA.txt.gz"), window_rows("chrA", 30, 30))
        _write_log(os.path.join(logs, "logchrB.txt.gz"), window_rows("chrB", 25, 27))
        # (sub-indexes with no target, whose record this stub's empty .bin files would not
        # get: index_fit records an index of a format it reads)
        # (each fake-indel genome cut by the windows its log records, as add-variants cuts
        # it: one whose windows are another width is refused, _window_disagreement)
        res, indexed = self.run_index({"chrA": window_fake("chrA", 30, 30),
                                       "chrB": window_fake("chrB", 25, 27),
                                       "chrC": f">fakechrC\n{NO_NRG}\n"},
                                      STUB_NO_BIN="fakechrA,fakechrB,fakechrC")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(indexed, ["chrA", "chrB", "chrC"])
        for chrom, want in (("chrA", [30, 30]), ("chrB", [25, 27]), ("chrC", None)):
            with open(os.path.join(self.sub(chrom), ".index_meta.json")) as fh:
                self.assertEqual(json.load(fh).get("indel_window"), want, chrom)

    def test_a_fake_genome_of_other_windows_than_its_log_is_refused(self):
        """W-R2-F1's class: a fake-indel genome cut by windows of another width than its indel
        log records (a download installs an archive's logs beside an enrichment made here by
        another CRISPRitz) is not indexed: the sub-index would record the log's windows, and
        be read with that log, every indel row misplaced. The others are built."""
        logs = os.path.join(self.tmp, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        _write_log(os.path.join(logs, "logchrA.txt.gz"), window_rows("chrA", 30, 30))
        _write_log(os.path.join(logs, "logchrB.txt.gz"), window_rows("chrB", 30, 30))
        res, indexed = self.run_index({"chrA": window_fake("chrA", 30, 30),
                                       "chrB": window_fake("chrB", 25, 27)},
                                      STUB_NO_BIN="fakechrA,fakechrB")
        self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
        # (its first window's end falls among the Ns this fixture puts between two windows;
        # the next one's past the sequence's end)
        self.assertIn("ERROR: the fake-indel genome of chrB (%s) and its indel log, in %s, come "
                      "from two enrichments: the window of chrB_200_ACG_A, which the log places "
                      "at 100-161 of the fake sequence, ends elsewhere there (no separator after "
                      "it)." % (os.path.join(self.tmp, "Genomes", "ref+ds_INDELS", "fake_ds_chrB"),
                                logs), res.stderr)
        self.assertIn("Indexing failed for INDELs in chrB\n", res.stdout)
        self.assertEqual(indexed, ["chrA"])

    def test_a_normal_build_is_unchanged(self):
        res, indexed = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n",
                                       "chrB": f">fakechrB\n{WITH_NRG}\n"})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(indexed, ["chrA", "chrB"])
        self.assertEqual(self.recorded(), [])
        for chrom in ("chrA", "chrB"):
            self.assertEqual(os.listdir(self.sub(chrom)), [f"NRG_fake{chrom}_1.bin"])
        self.assertEqual(self.built_with(), (["pam.txt", "pam.txt"], [NRG]))
        self.assertNotIn("stands in", res.stdout)
        self.assertNotIn("No INDELS sub-index", res.stdout)

    def test_every_line_of_a_worker_reaches_the_log(self):
        chroms = [f"chr{c}" for c in "ABCDEF"]
        res, indexed = self.run_index({c: f">fake{c}\n{WITH_NRG}\n" for c in chroms},
                                      PYTHONPATH=self.pool_exit_shim())
        self.assertEqual((res.returncode, res.stderr, indexed), (0, "", chroms))
        self.assertEqual(self.worker_lines(res.stdout, "Indexing"), sorted(
            [f"Indexing INDELs in {c}" for c in chroms]
            + [f"Indexing ended for INDELs in {c}" for c in chroms]))

    # A contig with an underscore in its name gets no sub-index, on any release: the build
    # names them, in one line (a dataset's name may hold an underscore: it does not count).

    def left_out(self, res):
        return [ln for ln in res.stdout.splitlines() if ln.startswith("No INDELS sub-index")]

    def test_a_contig_with_an_underscore_gets_no_sub_index_and_the_build_names_it(self):
        for dataset in ("ds", "hg38_ds"):
            with self.subTest(dataset=dataset):
                self.tearDown_dirs()
                res, _ = self.run_index(
                    {"chrA": f">fakechrA\n{WITH_NRG}\n",
                     UNPLACED: f">fake{UNPLACED}\n{WITH_NRG}\n",
                     RANDOM: f">fake{RANDOM}\n{WITH_NRG}\n"}, dataset=dataset)
                self.assertEqual((res.returncode, res.stderr), (0, ""))
                self.assertEqual(self.left_out(res), [LEFT_OUT.format(f"{RANDOM}, {UNPLACED}")])
                self.assertEqual(sorted(os.listdir(self.top())), ["NRG_3_fakechrA"])

    def test_a_dataset_whose_only_indels_are_on_such_a_contig_gets_an_index_with_no_sub_index(self):
        for true_pam in (None, "NNN"):  # (of a PAM-less stand-in too, under its own PAM)
            with self.subTest(true_pam=true_pam):
                self.tearDown_dirs()
                res, _ = self.run_index({"chrA": ">fakechrA\n", UNPLACED: f">fake{UNPLACED}\n{WITH_NRG}\n"},
                                        true_pam=true_pam)
                self.assertEqual((res.returncode, res.stderr), (0, ""))
                self.assertEqual(self.left_out(res), [LEFT_OUT.format(UNPLACED)])
                self.assertEqual(os.listdir(self.top()), [])

    def test_a_contig_with_snps_only_is_left_out_of_nothing(self):
        # (its fake-indel genome is header-only: no indels, and so no sub-index to miss)
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n", UNPLACED: f">fake{UNPLACED}\n"})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.left_out(res), [])

    def test_a_dataset_whose_name_holds_an_underscore_leaves_out_no_chromosome(self):
        res, indexed = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n"}, dataset="hg38_ds")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.left_out(res), [])
        self.assertEqual(indexed, ["chrA"])
        self.assertEqual(os.listdir(self.sub("chrA")), ["NRG_fakechrA_1.bin"])

    def test_another_datasets_directory_is_named_by_no_line(self):
        # (it is no contig of this dataset, whatever its name holds past this one's prefix)
        d = os.path.join(self.tmp, "Genomes", "ref+ds_INDELS", "fake_other_set_chrB")
        os.makedirs(d)
        with open(os.path.join(d, "fakechrB.fa"), "w") as fh:
            fh.write(f">fakechrB\n{WITH_NRG}\n")
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n"})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.left_out(res), [])

    def test_a_contig_read_otherwise_is_named_and_does_not_stop_the_build(self):
        # (a genome that could not be read may hold indels, and naming it is all this does)
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n", UNPLACED: b"\xff\xfe\n"})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.left_out(res), [LEFT_OUT.format(UNPLACED)])

    # A search builds a missing INDELS index with pool_index_indels.py, handing it the
    # PAM file of the search and the PAM of the variant index it resolved, which may be a
    # PAM-less one standing in for the search's (NNN_3 serving NRG). Built with the
    # search's PAM file, the INDELS index went under NRG_3_..., where the search did not
    # look; built with the index's own PAM, it is the one build-index-only builds.

    def test_a_pamless_stand_in_is_indexed_with_its_own_pam(self):
        for pam_build, line in ((None, NNN), ("20bp-NNN-NO-PAM 23 3\n", NNN),
                                ("21bp-NNN-NO-PAM 24 3\n", "N" * 24 + " 3")):
            with self.subTest(pam_build=pam_build):
                self.tearDown_dirs()
                res, indexed = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n",
                                               "chrB": f">fakechrB\n{NO_NRG}\n"},
                                              true_pam="NNN", pam_build=pam_build)
                self.assertEqual((res.returncode, res.stderr), (0, ""))
                self.assertEqual(indexed, ["chrA", "chrB"])
                files, lines = self.built_with()
                self.assertEqual(lines, [line])
                # (a file of its own, gone with the run)
                self.assertEqual(len(set(files)), 1)
                self.assertFalse(os.path.exists(files[0]), files[0])
                self.assertEqual(sorted(os.listdir(os.path.join(self.tmp, "genome_library"))),
                                 (["NNN_3_ref+ds"] if pam_build else []) + ["NNN_3_ref+ds_INDELS"])
                for chrom in ("chrA", "chrB"):
                    self.assertEqual(os.listdir(self.sub(chrom)), [f"NNN_fake{chrom}_1.bin"])
                self.assertIn("The variant index NNN_3_ref+ds is PAM-less and stands in for NRG: "
                              f"indexing the INDELs with its own PAM ({line}"
                              + (f", as {pam_build.split()[0]})" if pam_build else ")"), res.stdout)

    def test_a_pamless_stand_in_with_no_target_is_recorded(self):
        res, _ = self.run_index({"chrA": f">fakechrA\n{WITH_NRG}\n", "chrB": f">fakechrB\n{SHORT}\n"},
                                true_pam="NNN", STUB_NO_BIN="fakechrB")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.recorded(), ["chrB"])
        with open(os.path.join(self.sub("chrB"), NO_TARGETS)) as fh:
            self.assertEqual(fh.read(), "built, 0 targets: no window of fake_ds_chrB holds a target "
                                        "of NNN.txt with -bMax 3, so buildTST wrote no .bin\n")


@unittest.skipUnless(_gnu_sed(), "pool_search_indels.py merges with GNU `sed -i`, as the pipeline image does")
class IndexThenSearch(_Fixture):
    """A build that found no target on chrB, then a search of it with the fake-indel
    genome gone (as from a copied or published index): the record alone lets it pass."""

    def test_the_search_skips_what_the_build_recorded(self):
        work = os.path.join(self.tmp, "work")
        indels = os.path.join(work, "Genomes", "ref+ds_INDELS")
        for chrom, seq in (("chrA", WITH_NRG), ("chrB", NO_NRG)):
            os.makedirs(os.path.join(indels, f"fake_ds_{chrom}"))
            with open(os.path.join(indels, f"fake_ds_{chrom}", f"fake{chrom}.fa"), "w") as fh:
                fh.write(f">fake{chrom}\n{seq}\n")
        logdir = os.path.join(work, "Dictionaries", "log_indels_ds")
        os.makedirs(logdir)
        _write_log(os.path.join(logdir, "logchrA.txt.gz"), ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"])
        _write_log(os.path.join(logdir, "logchrB.txt.gz"), ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"])
        pam = os.path.join(self.tmp, "pam.txt")
        with open(pam, "w") as fh:
            fh.write(NRG + "\n")
        res = subprocess.run([sys.executable, INDEX, indels + "/", pam, "NRG", "ref", "ds", "3", "2"],
                             cwd=work, env=dict(self.env, STUB_NO_BIN="fakechrB"),
                             capture_output=True, text=True)
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        shutil.rmtree(indels)
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        open(os.path.join(self.tmp, "guide.txt"), "w").close()
        res = subprocess.run(
            [sys.executable, SEARCH, os.path.join(work, "Genomes", "ref"),
             os.path.join(self.tmp, "no_vcfs"), "ds", os.path.join(self.tmp, "guide.txt"), pam,
             "3", "4", "1", "1", out, "NRG", work + "/", "2", "-1", "0"],
            env=self.env, capture_output=True, text=True)
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual([c[1].rstrip("/").rsplit("_fake", 1)[-1] for c in self.calls() if c[0] == "search"],
                         ["chrA"])
        with open(os.path.join(out, "indels_ref+ds_pam.txt_guide.txt_4_1_1.targets.txt")) as fh:
            self.assertEqual(fh.read().split("\n"),
                             [HEADER, "X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1", ""])

    def test_a_search_finds_the_indels_index_it_built_for_a_pamless_stand_in(self):
        # (the search's own rebuild of a missing INDELS index, for an NNN_3 index serving
        # NRG): built with the search's PAM file, it went under NRG_3_..., and the search
        # of NNN_3_... fell back to each fake-indel genome without bulges
        work = os.path.join(self.tmp, "work")
        indels = os.path.join(work, "Genomes", "ref+ds_INDELS")
        logdir = os.path.join(work, "Dictionaries", "log_indels_ds")
        os.makedirs(logdir)
        for chrom, seq in (("chrA", WITH_NRG), ("chrB", NO_NRG)):
            os.makedirs(os.path.join(indels, f"fake_ds_{chrom}"))
            with open(os.path.join(indels, f"fake_ds_{chrom}", f"fake{chrom}.fa"), "w") as fh:
                fh.write(f">fake{chrom}\n{seq}\n")
            _write_log(os.path.join(logdir, f"log{chrom}.txt.gz"), LOGS_AB[chrom])
        pam = os.path.join(self.tmp, "pam.txt")
        with open(pam, "w") as fh:
            fh.write(NRG + "\n")
        res = subprocess.run([sys.executable, INDEX, indels + "/", pam, "NNN", "ref", "ds", "3", "2"],
                             cwd=work, env=self.env, capture_output=True, text=True)
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        open(os.path.join(self.tmp, "guide.txt"), "w").close()
        res = subprocess.run(
            [sys.executable, SEARCH, os.path.join(work, "Genomes", "ref"),
             os.path.join(self.tmp, "no_vcfs"), "ds", os.path.join(self.tmp, "guide.txt"), pam,
             "3", "4", "1", "1", out, "NNN", work + "/", "2", "-1", "0"],
            env=self.env, capture_output=True, text=True)
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        idx = os.path.join(work, "genome_library", "NNN_3_ref+ds_INDELS")
        self.assertEqual(sorted(c[1] for c in self.calls() if c[0] == "search"),
                         [f"{idx}/NNN_3_fakechrA/", f"{idx}/NNN_3_fakechrB/"])
        self.assertTrue(all("-index" in c for c in self.calls() if c[0] == "search"))
        with open(os.path.join(out, "indels_ref+ds_pam.txt_guide.txt_4_1_1.targets.txt")) as fh:
            self.assertEqual(fh.read().split("\n"),
                             [HEADER, ROW_A, ROW.format("chrB", 10, 10), ""])


class TheSearchScriptAtARefusal(_ShippedScript):
    """The shipped script with a stub engine, up to the real pool_search_indels.py and past
    it. A refusal ends the run at its own stage, with its ERROR alone: it went on to read
    the merged file that was not written ("awk: cannot open"), to take the empty file awk's
    redirect left for the indel targets, and through the SNP post-analysis, the failure
    then read "Post-analysis SNPs"."""

    ENGINE = r"""#!/bin/bash
printf '%s\n' "$*" >> "$STUB_LOG"
if [ "$1" = search ]; then
	printf '%s\n' HEADER > "$5.targets.txt"
fi
exit 0
""".replace("HEADER", repr(HEADER))

    def setUp(self):
        super().setUp()
        (self.root / "bin" / "crispritz.py").write_text(self.ENGINE)
        (self.home / "Dictionaries" / f"log_indels_{self.VCF}" / "logchr1.txt").write_text(
            f"{LOG_HEADER}\nchr1\t\trs1\t0.1\tA>AT\t5\tACGT\n")
        self.install("NRG_3")
        self.sub = self.library / f"NRG_3_hg38+{self.VCF}_INDELS" / "NRG_3_fakechr1"

    def run_it(self):
        launched = self.launch()
        run = self.finish(launched)
        with open(launched.work / "out" / "log.txt") as fh:
            stages = [ln.split("\t")[:2] for ln in fh]
        out = launched.work / "out"
        targets = sorted(p.name for p in out.rglob("indels_*.targets.txt*"))
        return launched.proc.returncode, run, stages, targets

    def test_it_ends_the_run_at_its_own_stage(self):
        for why, state in (("is not there, though chr1 has indels", None),
                           ("holds no .bin and no record that its build found no target there",
                            self.sub)):
            with self.subTest(why=why):
                if state:
                    state.mkdir()
                code, run, stages, targets = self.run_it()
                self.assertEqual(code, 1)
                self.assertRegex(run.stderr, r"\AERROR: the indel index of chr1 \(\S+\) " + why)
                self.assertEqual(len(run.stderr.splitlines()), 1, run.stderr)
                self.assertIn(["Off-targets search", "Start"], stages)
                self.assertNotIn(["Off-targets search", "End"], stages)
                self.assertNotIn("Post-analysis SNPs", [s[0] for s in stages])
                self.assertEqual(targets, [])
                self.assertNotIn("Search INDELs completed", run.stdout)

    def test_a_search_of_the_index_goes_on_past_its_stage(self):
        self.sub.mkdir()
        (self.sub / "NRG_fakechr1_1.bin").write_bytes(b"")
        _, run, stages, targets = self.run_it()
        self.assertIn("Search INDELs completed", run.stdout)
        self.assertIn(["Off-targets search", "End"], stages)
        self.assertIn(["Post-analysis SNPs", "Start"], stages)
        self.assertEqual(len(targets), 1)
        self.assertNotIn("ERROR: the indel index", run.stderr)

    def test_a_contig_no_build_gives_a_sub_index_goes_on_past_its_stage(self):
        # (it stopped every search with bulges of the dataset, for good)
        self.sub.mkdir()
        (self.sub / "NRG_fakechr1_1.bin").write_bytes(b"")
        (self.home / "Dictionaries" / f"log_indels_{self.VCF}" / f"log{UNPLACED}.txt").write_text(
            f"{LOG_HEADER}\n{UNPLACED}\t\trs1\t0.1\tA>AT\t5\tACGT\n")
        _, run, stages, targets = self.run_it()
        self.assertNotIn("ERROR: the indel index", run.stderr)
        self.assertIn("Search INDELs completed", run.stdout)
        self.assertIn(["Off-targets search", "End"], stages)
        self.assertEqual(len(targets), 1)


class ABuildOfAContigWithAnUnderscore(atomic._Fixture):
    """The real build-index-only, with the stub engine of test_build_index_only_atomic, of a
    dataset with an indel on a contig with an underscore in its name. The check that each
    chromosome holding a target got its .bin looked for a sub-index pool_index_indels.py
    has never made, and exited 1 taking it for a full disk's, where 2.8.0 exited 0."""

    # as index-genome does (checkExistance): a genome directory that is not there is said on
    # stdout, and it exits 0
    ENGINE = r'''#!PYTHON
import os, sys
if sys.argv[1] == "index-genome" and not os.path.isdir(sys.argv[3]):
    print("ERROR! %s is not a directory." % sys.argv[3])
    sys.exit()
os.execv(__file__ + ".stub", [__file__ + ".stub"] + sys.argv[1:])
'''

    @classmethod
    def setUpClass(cls):
        atomic.setUpModule()
        cls.addClassCleanup(atomic.tearDownModule)

    def fresh(self):
        super().fresh()
        engine = os.path.join(self.bin, "crispritz.py")
        os.rename(engine, engine + ".stub")
        with open(engine, "w") as fh:
            fh.write(self.ENGINE.replace("#!PYTHON", "#!" + sys.executable, 1))
        os.chmod(engine, 0o755)

    def build(self, contigs=(), dataset="ds", **stub):
        """build-index-only --vcf of the inputs, their VCFs in VCFs/<dataset>, with each of
        contigs added, an insertion on it."""
        vcfs = self.path("VCFs", dataset)
        if dataset != "ds":
            os.rename(self.path("VCFs", "ds"), vcfs)
        bases = atomic.GENOME["chr1"]
        for contig in contigs:
            atomic._write(self.path("Genomes", "hg38", f"{contig}.fa"), f">{contig}\n{bases}\n")
            with gzip.open(os.path.join(vcfs, f"{dataset}.{contig}.vcf.gz"), "wt") as fh:
                fh.write("##fileformat=VCFv4.2\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                         f"{contig}\t26\trs_c_26\t{bases[25]}\t{bases[25]}T\t.\tPASS\tAF=0.25"
                         "\tGT\t0|0\t1|0\n")
        return self.run_bio(command=[
            "build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam", atomic.PAM,
            "--bDNA", "1", "--bRNA", "1", "--thread", "2", "--path", self.work,
            "--vcf", vcfs, "--samplesID", self.path("samples.txt")], **stub)

    def test_it_exits_0_and_names_the_contig(self):
        rc, out, err = self.build((UNPLACED, RANDOM))
        self.assertEqual((rc, err), (0, ""), out)
        self.assertEqual([ln for ln in out.splitlines() if ln.startswith("No INDELS sub-index")],
                         [LEFT_OUT.format(f"{RANDOM}, {UNPLACED}")])
        self.assertEqual(self.bins(atomic.INDELS), ["NRG_2_fakechr1/NRG_fakechr1_1.bin",
                                                    "NRG_2_fakechr2/NRG_fakechr2_1.bin"])
        self.assertIn(f"Variant index ready: {atomic.VAR} (+ _INDELS)", out)

    def test_it_exits_0_where_such_a_contig_alone_has_indels(self):
        # (no other chromosome makes the INDELS index's directory, and an index with no
        # directory is taken for a failed build, where 2.8.0 exited 0)
        for chrom in atomic.GENOME:
            os.remove(self.path("Genomes", "hg38", f"{chrom}.fa"))
            os.remove(self.path("VCFs", "ds", f"ds.{chrom}.vcf.gz"))
        rc, out, err = self.build((UNPLACED,))
        self.assertEqual((rc, err), (0, ""), out)
        self.assertIn(LEFT_OUT.format(UNPLACED) + "\n", out)
        self.assertEqual(os.listdir(self.lib(atomic.INDELS)), [])
        self.assertIn(f"Variant index ready: {atomic.VAR} (+ _INDELS)", out)

    def test_a_chromosome_without_its_bin_still_fails_it(self):
        # with no such contig, beside one, and in a dataset whose name holds an underscore
        for contigs, dataset in (((), "ds"), ((UNPLACED,), "ds"), ((UNPLACED,), "hg38_ds")):
            with self.subTest(contigs=contigs, dataset=dataset):
                self.fresh()
                rc, out, err = self.build(contigs, dataset, STUB_NO_BIN="NRG_2_fakechr2 fakechr2.fa")
                self.assertEqual(rc, 1, out + err)
                self.assertIn("wrote no .bin for fakechr2 (", err)
                self.assertFalse(os.path.exists(self.lib(f"NRG_2_hg38+{dataset}_INDELS")))


class TheRecordsName(unittest.TestCase):
    def test_both_scripts_name_the_record_alike(self):
        for script in (INDEX, SEARCH):
            with open(script) as fh:
                tree = ast.parse(fh.read())
            names = [n.value.value for n in tree.body if isinstance(n, ast.Assign)
                     and [t.id for t in n.targets if isinstance(t, ast.Name)] == ["NO_TARGETS"]]
            self.assertEqual(names, [NO_TARGETS], script)

    def test_the_search_knows_complete_tests_view_by_its_name(self):
        with open(SEARCH) as fh:
            tree = ast.parse(fh.read())
        names = [n.value.value for n in tree.body if isinstance(n, ast.Assign)
                 and [t.id for t in n.targets if isinstance(t, ast.Name)] == ["CHROMVIEWPREFIX"]]
        self.assertEqual(names, [ct.CHROMVIEWPREFIX])


class LogCriterion(unittest.TestCase):
    """_logged_no_indels, lifted out of the script (it reads sys.argv at import)."""

    @classmethod
    def setUpClass(cls):
        with open(SEARCH) as fh:
            tree = ast.parse(fh.read())
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_logged_no_indels")
        ns = {"os": os, "gzip": gzip}
        exec(compile(ast.Module(body=[fn], type_ignores=[]), SEARCH, "exec"), ns)
        cls.f = staticmethod(ns["_logged_no_indels"])

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="indel_log_")
        self.addCleanup(shutil.rmtree, self.d)

    def write(self, name, text):
        p = os.path.join(self.d, name)
        if name.endswith(".gz"):
            with gzip.open(p, "wt") as fh:
                fh.write(text)
        else:
            with open(p, "w") as fh:
                fh.write(text)

    def test_header_only_is_empty(self):
        self.write("logchr1.txt.gz", LOG_HEADER + "\n")
        self.assertTrue(self.f(self.d, "chr1"))

    def test_one_record_is_not_empty(self):
        self.write("logchr1.txt.gz", LOG_HEADER + "\nchr1\t\trs1\t0.1\tA>AT\t5\tACGT\n")
        self.assertFalse(self.f(self.d, "chr1"))

    def test_a_record_without_a_header_still_counts(self):
        self.write("logchr1.txt", "chr1\t\trs1\t0.1\tA>AT\t5\tACGT\n")
        self.assertFalse(self.f(self.d, "chr1"))

    def test_plain_text_log_is_read_too(self):
        self.write("logchr1.txt", LOG_HEADER + "\n\n")
        self.assertTrue(self.f(self.d, "chr1"))

    def test_missing_log_is_unknown_not_empty(self):
        self.assertFalse(self.f(self.d, "chr1"))


class TheChromColumn(unittest.TestCase):
    """_dataset_chroms (and its helpers), lifted out of both indel scripts, refuses a
    dataset whose CHROM holds a character a shell reads. The CHROM of a dataset VCF's
    first record is taken, unquoted, into the search, merge and post-analysis commands
    both scripts run through a shell (os.system / shell=True): a VCF uploaded with a
    crafted CHROM would run it. A real chromosome name holds none, so a normal dataset is
    read as before. The twin helper in each script is checked; both must agree."""

    @staticmethod
    def _dataset_chroms(script):
        with open(script) as fh:
            tree = ast.parse(fh.read())
        wanted = {"_chrom_from_vcf", "_normalize_chrom", "_checked_chrom", "_dataset_chroms"}
        body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
        ns = {"os": os, "gzip": gzip}
        exec(compile(ast.Module(body=body, type_ignores=[]), script, "exec"), ns)
        return ns["_dataset_chroms"]

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="indel_chrom_")
        self.addCleanup(shutil.rmtree, self.d)
        self.marker = os.path.join(self.d, "MARK_from_chrom")

    def _vcf_dir(self, chrom):
        vdir = os.path.join(self.d, "VCFs", "ds")
        os.makedirs(vdir, exist_ok=True)
        with gzip.open(os.path.join(vdir, "ds.vcf.gz"), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write(f"{chrom}\t1\t.\tA\tAT\t.\t.\t.\n")
        return vdir

    def _log_dir(self, chrom):
        ldir = os.path.join(self.d, "Dictionaries", "log_indels_ds")
        os.makedirs(ldir, exist_ok=True)
        with open(os.path.join(ldir, f"log{chrom}.txt"), "w") as fh:
            fh.write("CHR\n")
        return ldir

    # the CHROM a crafted VCF (or a bundled log named after it) would carry into a shell
    CRAFTED = ("chr1$(touch MARK_from_chrom)", "chr1;touch MARK;#", "chr1 extra", "chr1|x")

    def test_a_vcf_chrom_a_shell_reads_is_refused(self):
        for script in (SEARCH, POSTINDEL):
            dataset_chroms = self._dataset_chroms(script)
            for chrom in self.CRAFTED:
                with self.subTest(script=os.path.basename(script), chrom=chrom):
                    vdir = self._vcf_dir(chrom)
                    with self.assertRaises(ValueError):
                        dataset_chroms(vdir, "")
                    self.assertFalse(os.path.exists(self.marker))
                    shutil.rmtree(os.path.dirname(vdir))

    def test_a_bundled_log_chrom_a_shell_reads_is_refused(self):
        for script in (SEARCH, POSTINDEL):
            dataset_chroms = self._dataset_chroms(script)
            with self.subTest(script=os.path.basename(script)):
                ldir = self._log_dir("chr1;touch MARK;#")
                with self.assertRaises(ValueError):
                    dataset_chroms("", ldir)
                self.assertFalse(os.path.exists(self.marker))
                shutil.rmtree(os.path.dirname(ldir))

    def test_a_normal_dataset_is_read_as_before(self):
        for script in (SEARCH, POSTINDEL):
            dataset_chroms = self._dataset_chroms(script)
            with self.subTest(script=os.path.basename(script)):
                # a raw VCF whose CHROM has no 'chr' prefix is normalized, as before
                self.assertEqual(dataset_chroms(self._vcf_dir("1"), ""), ["chr1"])
                shutil.rmtree(os.path.join(self.d, "VCFs"))
                # the log fallback (no raw VCFs) reads chr1 and chrUn_KI270302v1 (an
                # underscore is allowed), in order
                ldir = os.path.join(self.d, "Dictionaries", "log_indels_ds")
                os.makedirs(ldir, exist_ok=True)
                for name in ("logchr1.txt", "logchrUn_KI270302v1.txt"):
                    with open(os.path.join(ldir, name), "w") as fh:
                        fh.write("CHR\n")
                self.assertEqual(dataset_chroms("", ldir), ["chr1", "chrUn_KI270302v1"])
                shutil.rmtree(os.path.dirname(ldir))


if __name__ == "__main__":
    unittest.main()

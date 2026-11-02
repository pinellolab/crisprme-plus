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

Hermetic: a stub crispritz.py on PATH stands in for the engine and, like
CRISPRitz 2.9.1, refuses an index directory with no .bin; its index-genome writes
the layout the real one writes.
"""

import ast
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import complete_test as ct
import test_build_index_only_atomic as atomic
from test_index_resolution import _ShippedScript

HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = os.path.join(HERE, "pool_search_indels.py")
INDEX = os.path.join(HERE, "pool_index_indels.py")
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
        fh.write(HEADER + "\n")
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
AB = ("chrA", "chrB")
# contigs with an underscore in their name, which pool_index_indels.py does not index
UNPLACED, RANDOM = "chrUn_KI270302v1", "chr1_KI270706v1_random"
LEFT_OUT = ("No INDELS sub-index for {}: a contig with an underscore in its name gets none, so "
            "a search does not find the indel off-targets there that need a bulge")


def _gnu_sed():
    try:
        return subprocess.run(["sed", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


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


@unittest.skipUnless(_gnu_sed(), "pool_search_indels.py merges with GNU `sed -i`, as the pipeline image does")
class SearchSide(_Fixture):
    """pool_search_indels.py, chromosomes taken from the indel logs (no VCF dir)."""

    def run_search(self, logs, bins, records=(), fakes=None, true_pam="NRG", pam=NRG,
                   bulges=("1", "1"), vcfs=(), pam_files=None, pam_build=None, view=None,
                   listed="ds"):
        """bins: chrom -> whether its sub-index holds a .bin (a chromosome left out has
        none); records: the chromosomes whose sub-index holds the record of a build with
        no target; fakes: chrom -> the sequence of its fake-indel contig, written where
        build-index-only leaves it (None: a header-only one); vcfs: the chromosomes of the
        dataset's VCF folder, whence the search then takes them (else from the logs);
        pam_files: name -> line of a PAM file in the working directory's PAMs/; pam_build:
        what the variant index's .pam_build says; view: a chromosome, to search as
        complete-test --chrom does, from complete_test's own view of the working directory;
        listed: how the search's --vcf list spells the dataset."""
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
                open(os.path.join(d, f"{true_pam}_fake{chrom}_0.bin"), "w").close()
            if chrom in records:
                with open(os.path.join(d, NO_TARGETS), "w") as fh:
                    fh.write("built, 0 targets\n")
        if pam_build is not None:
            os.makedirs(os.path.join(cwd, "genome_library", f"{true_pam}_3_ref+ds"))
            with open(os.path.join(cwd, "genome_library", f"{true_pam}_3_ref+ds", ".pam_build"), "w") as fh:
                fh.write(pam_build)
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

    def test_snp_only_chromosome_is_skipped_not_fatal(self):
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"], "chrB": []},
                                     {"chrA": True, "chrB": False})
        self.assertEqual(res.stderr, "")  # a stage's stderr is fatal in the pipeline
        self.assertEqual(self.searched(), ["chrA"])
        self.assertIn("No indels on chrB", res.stdout)
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([ln for ln in lines[1:] if ln], ["X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"])

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
                self.assertEqual(os.listdir(self.sub("chrB")), [NO_TARGETS])
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
        self.assertEqual(os.listdir(self.sub("chrB")), [])
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


if __name__ == "__main__":
    unittest.main()

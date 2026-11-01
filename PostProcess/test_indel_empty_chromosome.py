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

The search's own rebuild of a missing INDELS index for a PAM-less index standing in
for its PAM (NNN_3 serving NRG) builds it with the index's own PAM (the PAM line its
variant index's .pam_build records). Built with the search's PAM file, that index went
under NRG_3_..., where the search did not look, and each chromosome was searched in its
fake-indel genome without bulges: exit 0, with the bulged indel off-targets lost
(measured through the real CLI).

Pinned here: the empty chromosome is skipped on BOTH sides (not indexed, not
searched), the merge tolerates its missing part file and only its, an all-empty
dataset still yields a well-formed targets file, a build with no target is
recorded and skipped, an index without the record is skipped only where its
fake-indel genome shows it holds no target, any other index that is empty although
its log lists indels fails with an ERROR that names how to rebuild it (no
traceback, no merged file), a PAM-less stand-in's INDELS index is built with its
own PAM, and the normal path is unchanged.

Hermetic: a stub crispritz.py on PATH stands in for the engine and, like
CRISPRitz 2.9.1, refuses an index directory with no .bin; its index-genome writes
the layout the real one writes.
"""

import ast
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

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
    bins = [f for f in os.listdir(index_dir) if f.endswith(".bin")] if os.path.isdir(index_dir) else []
    if not bins:
        sys.stderr.write(f"ERROR: no .bin index files found in '{index_dir}'.\n")
        sys.exit(1)
    chrom = os.path.basename(os.path.normpath(index_dir)).split("_fake")[-1]
    with open(out_name + ".targets.txt", "w") as fh:
        fh.write(HEADER + "\n")
        fh.write(f"X\tACGT\tACGT\tfake{chrom}\t10\t10\t+\t1\t0\t1\n")
elif sys.argv[1] == "index-genome":
    # What crispritz.py index-genome writes: genome_library/<PAM>_<bMax>_<name>, made
    # anew, the PAM read from the PAM file, a .bin for each contig of a .fa there but
    # those $STUB_NO_BIN names (buildTST writes none for a contig without a target).
    # $STUB_RECORD_BLOCKED leaves a directory where the record of no target would go.
    # The log gets the PAM file's line too ("pam-line <line>").
    name, gdir, pam_file = sys.argv[2], sys.argv[3], sys.argv[4]
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
# what the stub's search writes for a chromosome through its sub-index
ROW_A = "X\tACGT\tACGT\tfakechrA\t10\t10\t+\t1\t0\t1"
ROW = "X\tACGT\tACGT\tfake{}\t{}\t{}\t+\t1\t0\t1"
# the two indels the tests' logs list
LOGS_AB = {"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"], "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]}


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
        self.tmp = tempfile.mkdtemp(prefix="indel_empty_")
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

    def run_search(self, logs, bins, records=(), fakes=None, true_pam="NRG", pam=NRG):
        """bins: chrom -> whether its sub-index holds a .bin; records: the chromosomes
        whose sub-index holds the record of a build with no target; fakes: chrom -> the
        sequence of its fake-indel contig, written where build-index-only leaves it."""
        cwd = self.work = os.path.join(self.tmp, "work")
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
        for chrom, seq in (fakes or {}).items():
            d = os.path.join(cwd, "Genomes", "ref+ds_INDELS", f"fake_ds_{chrom}")
            os.makedirs(d)
            with open(os.path.join(d, f"fake{chrom}.fa"), "w") as fh:
                fh.write(f">fake{chrom}\n{seq}N\n")
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        open(os.path.join(self.tmp, "guide.txt"), "w").close()
        with open(os.path.join(self.tmp, "pam.txt"), "w") as fh:
            fh.write(pam + "\n")
        args = [os.path.join(cwd, "Genomes", "ref"), os.path.join(self.tmp, "no_vcfs"), "ds",
                os.path.join(self.tmp, "guide.txt"), os.path.join(self.tmp, "pam.txt"),
                "3", "4", "1", "1", out, true_pam, cwd + "/", "2", "-1", "0"]
        res = subprocess.run([sys.executable, SEARCH, *args], env=self.env,
                             capture_output=True, text=True)
        merged = os.path.join(out, "indels_ref+ds_pam.txt_guide.txt_4_1_1.targets.txt")
        lines = None
        if os.path.isfile(merged):
            with open(merged) as fh:
                lines = fh.read().split("\n")
        return res, lines

    def searched(self):
        return sorted(c[1].rstrip("/").rsplit("_fake", 1)[-1] for c in self.calls() if c[0] == "search")

    def refused(self, res, why, removed=(), again=""):
        """An ERROR for chrB that says why and names the rebuild, with nothing merged."""
        self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
        self.assertNotIn("Traceback", res.stderr)
        gone = " ".join("'%s'" % p for p in (self.idx, *removed))
        rebuild = (f"rm -rf {gone} && crisprme.py "
                   f"build-index-only --genome '{self.work}/Genomes/ref' --pam "
                   f"'{self.tmp}/pam.txt' --bDNA 2 --bRNA 2 --vcf '{self.tmp}/no_vcfs' "
                   f"--path '{self.work}'")
        self.assertEqual(res.stderr, (
            f"ERROR: the indel index of chrB ({self.idx}/NRG_3_fakechrB) holds no .bin and no "
            "record that its build found no target there (an index built before CRISPRme "
            f"2.8.1 has none), and {why}. Rebuild it: {rebuild}{again}\n"))
        self.assertEqual(self.searched(), ["chrA"])  # (and nothing searched in chrB)

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
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": False})
        variant = self.idx[: -len("_INDELS")]
        self.refused(res, "the fake-indel genome it was built from "
                          f"({self.work}/Genomes/ref+ds_INDELS/fake_ds_chrB) is not here to check "
                          "it against", removed=[f"{self.work}/Genomes/ref+ds"],
                     again=(f"; or, for an index you downloaded, download it again: rm -rf "
                            f"'{variant}' '{self.idx}' && crisprme.py download --what index "
                            f"--index-name NRG_3_ref+ds --path '{self.work}'"))
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
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": False},
                                     fakes={"chrA": WITH_NRG, "chrB": WITH_NRG})
        self.refused(res, f"{self.work}/Genomes/ref+ds_INDELS/fake_ds_chrB/fakechrB.fa holds a "
                          "target of NRG: the index is incomplete")
        self.assertIsNone(lines)

    def test_a_pamless_index_searched_for_another_pam_is_not_checked(self):
        # a PAM-less index serving NRG (--prune-pam): the search's PAM does not tell
        # whether a build with the index's own holds a target
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": False}, true_pam="NNN",
                                     fakes={"chrA": WITH_NRG, "chrB": NO_NRG})
        self.assertEqual(res.returncode, 1, res.stderr)
        self.assertNotIn("Traceback", res.stderr)
        self.assertIn(f"ERROR: the indel index of chrB ({self.idx}/NNN_3_fakechrB) holds no .bin",
                      res.stderr)
        self.assertIn("and it is an index of NNN searched for NRG, which this cannot check. "
                      f"Rebuild it: rm -rf '{self.idx}' && crisprme.py build-index-only "
                      f"--genome '{self.work}/Genomes/ref' --pam <the PAM file of NNN> ",
                      res.stderr)
        self.assertIsNone(lines)

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

    def run_index(self, fastas, pam=NRG, true_pam=None, pam_build=None, **env):
        """``true_pam``: the PAM of the variant index, when it is not that of ``pam``
        (a PAM-less index standing in for it, as a search passes it on), with what its
        .pam_build says."""
        folder = os.path.join(self.tmp, "Genomes", "ref+ds_INDELS")
        for chrom, text in fastas.items():
            d = os.path.join(folder, f"fake_ds_{chrom}")
            os.makedirs(d)
            if text is not None:
                with open(os.path.join(d, f"fake{chrom}.fa"), "w") as fh:
                    fh.write(text)
        with open(os.path.join(self.tmp, "pam.txt"), "w") as fh:
            fh.write(pam + "\n")
        true_pam = true_pam or (pam.split()[0][:4] if pam == TTTV else pam.split()[0][-3:])
        self.true_pam = true_pam
        if pam_build is not None:
            os.makedirs(os.path.join(self.tmp, "genome_library", f"{true_pam}_3_ref+ds"))
            with open(os.path.join(self.tmp, "genome_library", f"{true_pam}_3_ref+ds", ".pam_build"), "w") as fh:
                fh.write(pam_build)
        res = subprocess.run([sys.executable, INDEX, folder + "/", "pam.txt", true_pam, "ref", "ds", "3", "2"],
                             cwd=self.tmp, env=dict(self.env, **env), capture_output=True, text=True)
        indexed = sorted(c[2].rsplit("_", 1)[-1] for c in self.calls() if c[0] == "index-genome")
        return res, indexed

    def built_with(self):
        """The PAM file each index-genome call got, and the line it read there."""
        calls = self.calls()
        return ([c[3] for c in calls if c[0] == "index-genome"],
                sorted({" ".join(c[1:]) for c in calls if c[0] == "pam-line"}))

    def sub(self, chrom):
        return os.path.join(self.tmp, "genome_library", f"{self.true_pam}_3_ref+ds_INDELS",
                            f"{self.true_pam}_3_fake{chrom}")

    def recorded(self):
        top = os.path.join(self.tmp, "genome_library", f"{self.true_pam}_3_ref+ds_INDELS")
        return sorted(d.rsplit("_fake", 1)[-1] for d in os.listdir(top)
                      if os.path.isfile(os.path.join(top, d, NO_TARGETS)))

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


class TheRecordsName(unittest.TestCase):
    def test_both_scripts_name_the_record_alike(self):
        for script in (INDEX, SEARCH):
            with open(script) as fh:
                tree = ast.parse(fh.read())
            names = [n.value.value for n in tree.body if isinstance(n, ast.Assign)
                     and [t.id for t in n.targets if isinstance(t, ast.Name)] == ["NO_TARGETS"]]
            self.assertEqual(names, [NO_TARGETS], script)


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

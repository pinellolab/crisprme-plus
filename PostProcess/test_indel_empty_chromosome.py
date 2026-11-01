#!/usr/bin/env python3
"""A chromosome with SNPs but no indels must not abort a variant search.

The enricher writes one fake-indel contig per chromosome from the windows of the
indels some sample carries, so a chromosome with none gets a header-only FASTA and
a header-only indel log. buildTST cannot index that genome, which leaves the
chromosome's _INDELS sub-index without a .bin, and CRISPRitz 2.9.1
refuses to search such a directory -- so before this fix a single SNP-only
chromosome made the whole search fail, reported as an SNP post-analysis failure,
where 2.8.0 had returned the right answer (no indel off-target there). Reproduced
through the real CLI on chr21 (SNPs only) + chr22.

Pinned here: the empty chromosome is skipped on BOTH sides (not indexed, not
searched), the merge tolerates its missing part file and only its, an all-empty
dataset still yields a well-formed targets file, an index that is empty although
its log lists indels still fails loudly, and the normal path is unchanged.

Hermetic: a stub crispritz.py on PATH stands in for the engine and, like
CRISPRitz 2.9.1, refuses an index directory with no .bin.
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
import os, sys
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
'''.replace("HEADER", repr(HEADER))


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

    def run_search(self, logs, bins):
        cwd = os.path.join(self.tmp, "work")
        logdir = os.path.join(cwd, "Dictionaries", "log_indels_ds")
        os.makedirs(logdir)
        for chrom, records in logs.items():
            _write_log(os.path.join(logdir, f"log{chrom}.txt.gz"), records)
        idx = os.path.join(cwd, "genome_library", "NRG_3_ref+ds_INDELS")
        for chrom, has_bin in bins.items():
            d = os.path.join(idx, f"NRG_3_fake{chrom}")
            os.makedirs(d)
            if has_bin:
                open(os.path.join(d, f"NRG_fake{chrom}_0.bin"), "w").close()
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        for name in ("guide.txt", "pam.txt"):
            open(os.path.join(self.tmp, name), "w").close()
        args = [os.path.join(cwd, "Genomes", "ref"), os.path.join(self.tmp, "no_vcfs"), "ds",
                os.path.join(self.tmp, "guide.txt"), os.path.join(self.tmp, "pam.txt"),
                "3", "4", "1", "1", out, "NRG", cwd + "/", "2", "-1", "0"]
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
        # a broken index must never pass as "no indels"
        res, _ = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                  "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                 {"chrA": True, "chrB": False})
        self.assertIn("no .bin index files found", res.stderr)
        self.assertEqual(self.searched(), ["chrA", "chrB"])

    def test_normal_path_is_unchanged(self):
        res, lines = self.run_search({"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"],
                                      "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]},
                                     {"chrA": True, "chrB": True})
        self.assertEqual(res.stderr, "")
        self.assertEqual(lines[0], HEADER)
        self.assertEqual([ln.split("\t")[3] for ln in lines[1:] if ln], ["fakechrA", "fakechrB"])


class IndexSide(_Fixture):
    """pool_index_indels.py must not hand an empty fake genome to buildTST."""

    def run_index(self, fastas):
        folder = os.path.join(self.tmp, "Genomes", "ref+ds_INDELS")
        for chrom, text in fastas.items():
            d = os.path.join(folder, f"fake_ds_{chrom}")
            os.makedirs(d)
            if text is not None:
                with open(os.path.join(d, f"fake{chrom}.fa"), "w") as fh:
                    fh.write(text)
        res = subprocess.run([sys.executable, INDEX, folder + "/", "pam.txt", "NRG", "ref", "ds", "3", "2"],
                             cwd=self.tmp, env=self.env, capture_output=True, text=True)
        indexed = sorted(c[2].rsplit("_", 1)[-1] for c in self.calls() if c[0] == "index-genome")
        return res, indexed

    def test_header_only_fake_genome_is_not_indexed(self):
        res, indexed = self.run_index({"chrA": ">fakechrA\nACGTNNACGT\n", "chrB": ">fakechrB\n"})
        self.assertEqual(res.stderr, "")
        self.assertEqual(indexed, ["chrA"])
        self.assertIn("No indels on chrB: nothing to index", res.stdout)

    def test_a_dir_without_any_fasta_is_not_mistaken_for_empty(self):
        _, indexed = self.run_index({"chrA": ">fakechrA\nACGT\n", "chrC": None})
        self.assertEqual(indexed, ["chrA", "chrC"])


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

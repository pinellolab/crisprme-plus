#!/usr/bin/env python3
"""validate-test where there is no complete-test benchmark output to check, and on a
chromosome that a benchmark's brute-force reference does not cover.

validate-test compares the off-targets that complete-test's BENCHMARK mode wrote
(CRISPRME_ALLOW_ONDEMAND_BUILD=1 ... --vcf_dataset 1000G, as validate-benchmarks.yml
runs it) with the committed brute-force ground truth, and starts by reading the
vcf.config.test.txt that mode writes. Run anywhere else -- a fresh folder, or after
the default prebuilt-index smoke, which writes no such file -- it ended in two Python
tracebacks: validate.py's FileNotFoundError, then crisprme.py's OSError. It now says
in one line what to run first, and exits 1.

That line names the pair for the chromosome asked for, but the Cas12a reference and
the single-"n edits" ones cover chr22 only ("scope" in benchmarks.json): compared on
chr21, or on a whole genome, each site CRISPRme found off chr22 counted as extra and
the benchmark FAILED. Such a benchmark is now skipped there.

Where a benchmark is compared, its brute-force sites are read by position, with no
pandas FutureWarning (and no KeyError on pandas 3).

validate.py is run as a subprocess, as crisprme.py runs it, and its
run_test_validation in-process with each comparison stubbed. crisprme.py's
validate_test is exec'd without the CLI dispatch (the truncate-and-exec loader of
test_bulge_cap), with its subprocess call stubbed.
"""
import contextlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATE = os.path.join(HERE, "validate.py")


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_bulge_cap._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    crisprme_py = os.path.join(os.path.dirname(HERE), "crisprme.py")
    sys.path.insert(0, HERE)
    with open(crisprme_py) as fh:
        src = fh.read()
    marker = "\nif len(sys.argv) < 2:"
    cut = src.index(marker)  # KeyError-loud if the dispatch shape ever changes
    mod = types.ModuleType("crisprme_under_test_validate")
    mod.__file__ = crisprme_py
    code = compile(src[:cut], crisprme_py, "exec")
    # Stubs exist ONLY so crisprme.py's top-level imports resolve while we exec its
    # definitions; they are removed again in the finally, so they cannot shadow the real
    # modules for a test collected later.
    _inserted = []
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        if "Bio" not in sys.modules:
            sys.modules["Bio"] = _bio
            _inserted.append("Bio")
        sys.modules["Bio.Seq"] = _bio_seq
        _inserted.append("Bio.Seq")
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISPR_BULGE_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = mock.MagicMock()
            _inserted.append(_m)
    try:
        exec(code, mod.__dict__)
    finally:
        for _m in _inserted:
            sys.modules.pop(_m, None)
    return mod


crisprme = _load_crisprme()
import validate  # noqa: E402 -- PostProcess is on sys.path now (_load_crisprme)


class ValidatePyWithoutBenchmarkOutput(unittest.TestCase):
    def test_one_line_on_stderr_and_exit_1(self):
        with tempfile.TemporaryDirectory() as d:  # a folder complete-test never ran in
            out = subprocess.run([sys.executable, VALIDATE, "chr22"], cwd=d,
                                 capture_output=True, text=True)
        self.assertEqual(out.returncode, 1, out.stderr)
        self.assertEqual(out.stdout, "")
        self.assertNotIn("Traceback", out.stderr)
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        # it says what is missing and what to run first: the documented benchmark mode
        for part in ("vcf.config.test.txt", "CRISPRME_ALLOW_ONDEMAND_BUILD=1",
                     "complete-test --chrom chr22 --vcf_dataset 1000G"):
            self.assertIn(part, lines[0])

    def hint(self, chrom):
        """validate.py's one line, run as validate-test runs it for `chrom`."""
        with tempfile.TemporaryDirectory() as d:
            out = subprocess.run([sys.executable, VALIDATE, chrom], cwd=d,
                                 capture_output=True, text=True)
        self.assertEqual(out.returncode, 1, out.stderr)
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        return lines[0]

    def test_the_commands_are_for_the_chromosome_asked_for(self):
        # validate-test --chrom chr21 reads hg38_chr21 reports, which a chr22 run never writes
        line = self.hint("chr21")
        self.assertIn("complete-test --chrom chr21 --vcf_dataset 1000G", line)
        self.assertIn("then crisprme.py validate-test --chrom chr21", line)
        self.assertNotIn("chr22", line)

    def test_without_chrom_they_are_the_documented_chr22_pair(self):
        # validate-test without --chrom ("all") checks a whole-genome run: rather than
        # that, the chr22 self-test, and the validate-test that then has something to check
        line = self.hint("all")
        self.assertIn("complete-test --chrom chr22 --vcf_dataset 1000G", line)
        self.assertIn("then crisprme.py validate-test --chrom chr22", line)

    def test_it_says_the_variable_goes_with_minus_e_in_docker(self):
        # set in front of `docker run`, it stays on the host: in the container
        # complete-test then runs the prebuilt-index smoke, which leaves nothing to check
        line = self.hint("chr22")
        self.assertIn("--vcf_dataset 1000G (Docker: pass it with -e), then crisprme.py "
                      "validate-test --chrom chr22", line)

    def test_it_runs_every_benchmark(self):
        # no CRISPRME_SKIP_HEAVY: the set fits 16 GiB, and a skipped benchmark checks nothing
        self.assertNotIn("CRISPRME_SKIP_HEAVY", self.hint("chr22"))


class ABenchmarkOffItsReferenceIsSkipped(unittest.TestCase):
    """After the `complete-test --chrom chr21 --vcf_dataset 1000G` that the line above
    names for chr21, `validate-test --chrom chr21` reported '[cas12a_hbg] FAILED:
    mismatch. CRISPRme sites: 223  Brute-force sites: 0': its reference covers chr22
    only, so every chr21 site counted as extra. A benchmark is now compared only on
    the chromosome its reference covers, and skipped elsewhere."""

    # what the benchmark mode runs with the CRISPRME_SKIP_HEAVY=1 opt-out
    LIGHT = ["cas9_sg1617", "cas12a_hbg"]

    def run_validation(self, chrom, ran, skip_heavy=False):
        """run_test_validation(chrom) on the shipped registry, with CRISPRME_SKIP_HEAVY
        unset as in the validate-test that line names (or set, the opt-out), where the
        benchmarks in `ran` left their output, each comparison stubbed to a match. Returns
        (exit code, stderr lines, the benchmarks compared)."""
        compared = []
        err = io.StringIO()
        env = dict(os.environ)
        env.pop("CRISPRME_SKIP_HEAVY", None)
        if skip_heavy:
            env["CRISPRME_SKIP_HEAVY"] = "1"

        def found(name, *_a):
            return (f"{name}.ref.txt", f"{name}.alt.txt") if name in ran else None

        def match(name, *_a):
            compared.append(name)
            return True

        with mock.patch.dict(os.environ, env, clear=True), \
             mock.patch.object(validate, "check_variant_dataset"), \
             mock.patch.object(validate, "find_crispritz_targets", side_effect=found), \
             mock.patch.object(validate, "download_brute_force_targets", return_value="bf.tsv"), \
             mock.patch.object(validate, "load_targets", return_value=(None, None)), \
             mock.patch.object(validate, "validate", side_effect=match), \
             contextlib.redirect_stderr(err):
            with self.assertRaises(SystemExit) as cm:
                validate.run_test_validation(chrom)
        return cm.exception.code, err.getvalue().splitlines(), compared

    def every(self):
        return [b["name"] for b in validate.load_benchmarks()["benchmarks"]]

    def test_chr21_compares_only_what_the_references_cover(self):
        every = self.every()
        code, lines, compared = self.run_validation("chr21", every)
        self.assertEqual(compared, ["cas9_sg1617"])
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not chr21.", lines)
        self.assertIn(f"Validation summary: 1 passed, 0 failed, {len(every) - 1} skipped.", lines)
        self.assertEqual(code, 0)

    def test_so_does_a_whole_genome_run(self):
        # validate-test without --chrom
        every = self.every()
        code, lines, compared = self.run_validation("all", every)
        self.assertEqual(compared, ["cas9_sg1617"])
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not the whole genome.", lines)
        self.assertIn(f"Validation summary: 1 passed, 0 failed, {len(every) - 1} skipped.", lines)
        self.assertEqual(code, 0)

    def test_chr22_compares_every_benchmark(self):
        # the documented self-test: every benchmark passed, none failed, none skipped
        every = self.every()
        code, lines, compared = self.run_validation("chr22", every)
        self.assertEqual(compared, every)
        self.assertIn(f"Validation summary: {len(every)} passed, 0 failed, 0 skipped.", lines)
        self.assertEqual(code, 0)

    def test_the_opt_out_leaves_out_the_bulge_2_benchmarks_saying_so(self):
        # CRISPRME_SKIP_HEAVY=1, for a very small machine, on complete-test and validate-test
        code, lines, compared = self.run_validation("chr22", self.LIGHT, skip_heavy=True)
        self.assertEqual(compared, self.LIGHT)
        heavy = [n for n in self.every() if n not in self.LIGHT]
        for name in heavy:
            self.assertIn(f"[{name}] SKIPPED: bulge-2 benchmark, left out by "
                          "CRISPRME_SKIP_HEAVY=1.", lines)
        self.assertIn(f"Validation summary: 2 passed, 0 failed, {len(heavy)} skipped.", lines)
        self.assertEqual(code, 0)

    def test_an_opt_out_on_complete_test_alone_is_said_too(self):
        # complete-test with it, validate-test without: what did not run is not a failure
        every = self.every()
        code, lines, compared = self.run_validation("chr22", self.LIGHT)
        self.assertEqual(compared, self.LIGHT)
        self.assertIn(f"Validation summary: 2 passed, 0 failed, {len(every) - 2} skipped.", lines)
        self.assertEqual(code, 0)

    def test_validate_py_no_longer_fails_chr21_sites_against_the_chr22_reference(self):
        # validate.py itself and the committed Cas12a reference, on chr21 output of
        # cas12a_hbg's: before, the real comparison counted each of these sites as extra
        registry = validate.load_benchmarks()
        bench = next(b for b in registry["benchmarks"] if b["name"] == "cas12a_hbg")
        th = registry["thresholds"]
        site = "X\tTTTVCCTTGTCAAGGCTATTGGTC\tTTTACCTTGTCAAGGCTATTGGTC\tchr21\t{}\t{}\t+\t0\t0\t0\n"
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "vcf.config.test.txt"), "w") as fout:
                fout.write("hg38_1000G\n")
            tdir = os.path.join(d, validate.targets_dir("cas12a_hbg"))
            os.makedirs(tdir)
            names = validate.target_report_names(bench["pam_name"], th["mm"], th["bDNA"],
                                                 th["bRNA"], "chr21")
            for fname, positions in zip(names, ((5010000, 5020000), (5030000,))):
                with open(os.path.join(tdir, fname), "w") as fout:
                    fout.writelines(site.format(p, p) for p in positions)
            out = subprocess.run([sys.executable, VALIDATE, "chr21"], cwd=d,
                                 capture_output=True, text=True)
        self.assertIn("[cas12a_hbg] SKIPPED: its brute-force reference covers chr22 only, "
                      "not chr21.", out.stderr.splitlines(), out.stderr)
        self.assertNotIn("FAILED", out.stderr)


REFERENCES = os.path.join(os.path.dirname(HERE), "test", "benchmark")


def registered(name):
    """A benchmark of the shipped registry, its thresholds resolved as complete-test does."""
    registry = validate.load_benchmarks()
    bench = next(b for b in registry["benchmarks"] if b["name"] == name)
    return bench, dict(registry["thresholds"], **bench.get("thresholds", {}))


class TheFiveEditBenchmarksCatchALostSite(unittest.TestCase):
    """cas12a_hbg_simple's reference is a header with no rows, and cas9_sg1617_simple's holds
    3 rows at 2 sites: a search that lost its sites still passed them. At 5 total edits the
    same guides, on the same bulge-2 indexes, have 381 and 92 sites. validate.py as
    crisprme.py runs it, on these committed references, with a complete-test output made
    from the reference itself, less what a test takes out."""

    SITES = {"cas9_sg1617_n5": 381, "cas12a_hbg_n5": 92}

    def rows(self, name):
        with open(os.path.join(REFERENCES, registered(name)[0]["reference"])) as fin:
            fin.readline()
            return [line.rstrip("\n").split("\t") for line in fin]

    def validate_py(self, lose=(), gain=()):
        """(exit status, stderr lines) of `validate.py chr22` where these two benchmarks ran,
        each site of `lose` missing from the search's output and each of `gain` added."""
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "vcf.config.test.txt"), "w") as fout:
                fout.write("hg38_1000G\n")
            for name in self.SITES:
                bench, th = registered(name)
                tdir = os.path.join(d, validate.targets_dir(name))
                os.makedirs(tdir)
                ref, alt = validate.target_report_names(bench["pam_name"], th["mm"], th["bDNA"],
                                                        th["bRNA"], "chr22")
                open(os.path.join(tdir, ref), "w").close()  # all in the variant report
                with open(os.path.join(tdir, alt), "w") as fout:
                    for chrom, rna, dna, strand, start, *_ in self.rows(name) + list(gain):
                        if (chrom, start, strand) not in lose:
                            fout.write(f"X\t{rna}\t{dna}\t{chrom}\t{start}\t{start}\t{strand}\t0\t0\t0\n")
            out = subprocess.run([sys.executable, VALIDATE, "chr22"], cwd=d, capture_output=True,
                                 text=True)
        return out.returncode, out.stderr.splitlines()

    def test_they_hold_what_they_are_there_for(self):
        for name, sites in self.SITES.items():
            with self.subTest(benchmark=name):
                rows = self.rows(name)
                self.assertEqual(len({(r[0], r[4], r[3]) for r in rows}), sites)
                self.assertEqual({r[0] for r in rows}, {"chr22"})
                bench, th = registered(name)
                self.assertEqual(th["max_total_edits"], 5)
                self.assertEqual(bench["scope"], "chr22")
                self.assertTrue(bench["heavy"])  # they search the bulge-2 index

    def test_every_site_found_passes(self):
        code, lines = self.validate_py()
        self.assertEqual(code, 0, lines)
        for name in self.SITES:
            self.assertIn(f"[{name}] PASSED: CRISPRme and brute-force off-target sets match.", lines)
        self.assertIn("Validation summary: 2 passed, 0 failed, 4 skipped.", lines)

    def test_one_lost_site_fails_it(self):
        for name in self.SITES:
            with self.subTest(benchmark=name):
                row = self.rows(name)[len(self.rows(name)) // 2]
                code, lines = self.validate_py(lose={(row[0], row[4], row[3])})
                self.assertEqual(code, 1, lines)
                self.assertIn(f"[{name}] FAILED: mismatch.", lines)
                self.assertIn("  Missing in CRISPRme: 1  Extra in CRISPRme: 0", lines)
                self.assertIn("Validation summary: 1 passed, 1 failed, 4 skipped.", lines)

    def test_one_extra_site_fails_it_too(self):
        extra = ("chr22", "CTAACAGTTGCTTTTATCACNGG", "CTAACAGTTGCTTTTATCACAGG", "+", "20000000")
        code, lines = self.validate_py(gain=[extra])
        self.assertEqual(code, 1, lines)
        self.assertIn("  Missing in CRISPRme: 0  Extra in CRISPRme: 1", lines)
        self.assertIn("Validation summary: 0 passed, 2 failed, 4 skipped.", lines)


class EveryInstallFindsTheFiveEditReferences(unittest.TestCase):
    """benchmarks.json's reference_base_url is the main branch, and validate.py goes there
    for a reference the install does not hold. Docker and source installs hold test/ whole,
    and CI copies the branch's over the image, so a reference committed with a change is found
    before that change is merged. The Bioconda recipe leaves test/benchmark/brute-force-* out
    of its package: a new reference there would be a 404 on a Bioconda build of a release
    until that release is on main, so these two live outside it."""

    # the recipe's own line (recipes/crisprme/build.sh in bioconda-recipes)
    RECIPE_LEAVES_OUT = "brute-force-*"

    def in_the_bioconda_layout(self, reference, md5, download):
        """download_brute_force_targets on a tree with what the recipe packages under
        test/benchmark, in a results folder of its own, `download` stubbed."""
        with tempfile.TemporaryDirectory() as d:
            shutil.copytree(REFERENCES, os.path.join(d, "test", "benchmark"),
                            ignore=shutil.ignore_patterns(self.RECIPE_LEAVES_OUT))
            os.makedirs(os.path.join(d, "PostProcess"))
            cwd = os.getcwd()
            os.chdir(d)  # BFDIR is relative to the working folder
            try:
                with mock.patch.object(validate, "__file__", os.path.join(d, "PostProcess", "validate.py")), \
                     mock.patch.object(validate, "download", side_effect=download) as fetched, \
                     contextlib.redirect_stderr(io.StringIO()):
                    path = validate.download_brute_force_targets("https://main.example/", reference, md5)
                return os.path.relpath(path, d) if path.startswith(d) else path, fetched.call_args
            finally:
                os.chdir(cwd)

    def test_they_need_no_network_there(self):
        def no_network(*_a, **_k):
            raise AssertionError("downloaded")
        for name in TheFiveEditBenchmarksCatchALostSite.SITES:
            with self.subTest(benchmark=name):
                bench, _ = registered(name)
                path, fetched = self.in_the_bioconda_layout(bench["reference"], bench["md5"], no_network)
                self.assertEqual(os.path.normpath(path),
                                 os.path.join("test", "benchmark", bench["reference"]))
                self.assertIsNone(fetched)

    def test_a_reference_in_a_brute_force_folder_is_downloaded_from_main_there(self):
        # what the recipe's line does to the others, so what a new one must stay out of
        bench, _ = registered("cas12a_hbg_simple")
        ref = os.path.join(REFERENCES, bench["reference"])

        def from_main(dest, http_url):
            copy = os.path.join(dest, os.path.basename(http_url))
            shutil.copyfile(ref, copy)
            return copy
        path, fetched = self.in_the_bioconda_layout(bench["reference"], bench["md5"], from_main)
        self.assertEqual(fetched.kwargs["http_url"], "https://main.example/" + bench["reference"])
        self.assertTrue(path.endswith(os.path.basename(bench["reference"])))

    def test_every_committed_reference_has_its_registered_md5(self):
        # a wrong one would send validate.py past the committed copy to the network
        for bench in validate.load_benchmarks()["benchmarks"]:
            with self.subTest(benchmark=bench["name"]):
                path = os.path.join(REFERENCES, bench["reference"])
                self.assertEqual(validate.compute_md5(path), bench["md5"])


class TheDriverRegeneratesThemWithTheirOwnBudgets(unittest.TestCase):
    """generate_references.py makes a benchmark's reference from its registry entry. It must
    ask the generator for what complete-test asks the search: a 5-edit cap binding over the
    wide-open per-type budgets. Run with a generator that records what it was asked."""

    STUB = ("import sys\n"
            "open(sys.argv[sys.argv.index('--output') + 1], 'w').write('CHR\\n')\n"
            "open(sys.argv[0] + '.args', 'a').write(' '.join(sys.argv[1:]) + '\\n')\n")

    def asked(self, name):
        with tempfile.TemporaryDirectory() as d:
            stub = os.path.join(d, "generator.py")
            with open(stub, "w") as fout:
                fout.write(self.STUB)
            subprocess.run([sys.executable, os.path.join(REFERENCES, "generate_references.py"),
                            "--enriched-fasta", "chr22.enriched.fa", "--chrom", "chr22",
                            "--only", name, "--generator", stub, "--out-dir", d],
                           capture_output=True, text=True, check=True)
            with open(stub + ".args") as fin:
                line = fin.read()
        return line

    def test_cas9(self):
        args = self.asked("cas9_sg1617_n5")
        for part in ("--rna CTAACAGTTGCTTTTATCACNGG", "--max-mismatches 5", "--max-dna-gaps 2",
                     "--max-rna-gaps 2", "--chrom chr22", "--pam GG", "--max-total-edits 5",
                     "cas9_sg1617_n5.chr22.tsv"):
            self.assertIn(part, args)
        self.assertNotIn("--pam-5prime", args)

    def test_cas12a(self):
        args = self.asked("cas12a_hbg_n5")
        for part in ("--rna TTTVCCTTGTCAAGGCTATTGGTC", "--max-mismatches 5", "--max-dna-gaps 2",
                     "--max-rna-gaps 2", "--chrom chr22", "--pam TTTV", "--max-total-edits 5",
                     "--pam-5prime", "cas12a_hbg_n5.chr22.tsv"):
            self.assertIn(part, args)


class TheBruteForceSitesAreReadByPosition(unittest.TestCase):
    """The brute-force reference is headered (CHR RNA DNA Strand Start ...), and its
    sites were read with x[0], x[4] and x[3] on rows labelled by column name: each
    benchmark validate-test compared printed a pandas FutureWarning ("Series.__getitem__
    treating keys as positions is deprecated"), and pandas 3 raises a KeyError there.
    The committed Cas12a reference, read as validate-test reads it."""

    def test_its_sites_are_chr_start_and_strand_with_no_warning(self):
        import warnings
        ref = os.path.join(os.path.dirname(HERE), "test", "benchmark",
                           "brute-force-cas12-1000G", "brute_force_cas12_1000G.tsv")
        targets = validate._load_targets(ref, "chr22")
        expected = [f"{c}_{p}_{s}" for c, p, s in
                    zip(targets["CHR"], targets["Start"], targets["Strand"])]
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            sites = validate.compute_sites(targets)
        self.assertEqual(sites, expected)
        self.assertGreater(len(set(sites)), 100)


class TheCliCallerAddsNoTraceback(unittest.TestCase):
    """validate.py has said why on stderr; crisprme.py only passes the failure on."""

    ARGV = ["crisprme.py", "validate-test", "--chrom", "chr22"]

    def test_a_failed_validation_exits_1(self):
        with mock.patch.object(crisprme, "input_args", self.ARGV), \
             mock.patch.object(crisprme, "subprocess") as msub:
            msub.call.return_value = 1
            with self.assertRaises(SystemExit) as cm:
                crisprme.validate_test()
        self.assertEqual(cm.exception.code, 1)
        (cmd,), _ = msub.call.call_args
        self.assertTrue(cmd.endswith("validate.py chr22"), cmd)

    def test_a_passed_validation_returns(self):
        with mock.patch.object(crisprme, "input_args", self.ARGV), \
             mock.patch.object(crisprme, "subprocess") as msub:
            msub.call.return_value = 0
            self.assertIsNone(crisprme.validate_test())


if __name__ == "__main__":
    unittest.main(verbosity=2)

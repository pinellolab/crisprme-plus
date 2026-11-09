"""Option P (2.10.0): a chromosome's post-analysis in pieces over one pool of workers,
merged into the files one process writes, byte for byte.

  * merge_pieces.score_cbulge scores the pieces' rows in exactly the batches one process
    forms (every cbulge_rows.BATCH_ROWS rows, cut after the input line that reaches it): the
    model's floats depend on which rows share a batch, so other batches would change the
    last printed digit of some CRISPR-Bulge scores (the cost judges measured it). Checked
    with a scorer whose output depends on the whole batch, against the loop of one process;
  * cbulge_rows.score_rows scores a batch whose rows it nearly all nulls at once (it read a
    list of them for each row);
  * the companions keep the first row of each identity, as one process does; the dense-window
    BED keeps one line a window, the one resultIntegrator reads first; the four-allele
    positions are the sorted union; the CFD graphs are summed;
  * postproc_pool cuts a chromosome into contiguous ranges, runs the stages of every
    chromosome on one pool, and stops at the first failure; a piece's output goes to its own
    file, and the log keeps one piece's start-up and what the others add;
  * through the real scripts: new_simple_analysis.py and analisi_indels_NNN.py run as one
    process and as pieces (CRISPRME_POSTPROC_RANGE, CRISPRME_POSTPROC_PIECE) then merged give
    the same files, byte for byte, with the pre-check on and off; the pieces write none of
    the notes a chromosome ends with, the merge writes one process's.

Run with:
    cd PostProcess && python3 -m unittest test_postproc_pieces -v
"""

import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cbulge_rows  # noqa: E402
import merge_pieces  # noqa: E402
import postproc_pool  # noqa: E402


def _batch_scorer(log):
    """A scorer whose scores depend on the whole batch (as the model's last bits do):
    each row gets the batch's number and size, and the batch is logged."""
    def score(rows):
        n = len(log)
        log.append([r[0] for r in rows])
        out = []
        for r in rows:
            t = list(r)
            t[-2] = "b%d/%d" % (n, len(rows))
            t.append("x")
            out.append(t)
        return out
    return score


def _one_process(lines, batch_rows, score):
    """The loop of new_simple_analysis.py / analisi_indels_NNN.py: the rows of each line
    gather, and after a line that brings them to batch_rows they are scored; the rest at
    the end. Returns the written rows."""
    out, buf = [], []
    for rows in lines:
        buf.extend(rows)
        if len(buf) >= batch_rows:
            out.extend(score(buf))
            buf = []
    if buf:
        out.extend(score(buf))
    return out


class ScoredInTheBatchesOfOneProcess(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="pieces_cbulge_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def write_pieces(self, lines, cuts):
        """The pieces' row files: line i's rows, [line, row] per row, cut at ``cuts``."""
        files = []
        bounds = [0] + cuts + [len(lines)]
        for k in range(len(bounds) - 1):
            path = os.path.join(self.tmp, "p%d.cbulge_rows.jsonl" % k)
            with open(path, "w") as fh:
                for i in range(bounds[k], bounds[k + 1]):
                    for row in lines[i]:
                        fh.write(json.dumps([i + 1, row]) + "\n")
            files.append(path)
        return files

    def test_random_lines_and_cuts(self):
        rng = random.Random(1)
        for trial in range(300):
            n_lines = rng.randint(0, 60)
            lines = []
            rid = 0
            for _ in range(n_lines):
                k = rng.choice([0, 0, 1, 1, 2, 3, 7])
                lines.append([["r%d" % (rid + j), "A", "C", 33, 0] for j in range(k)])
                rid += k
            batch = rng.randint(1, 12)
            cuts = sorted(rng.sample(range(1, n_lines), rng.randint(0, min(5, n_lines - 1))))\
                if n_lines > 1 else []
            want_log, got_log = [], []
            want = _one_process(lines, batch, _batch_scorer(want_log))
            files = self.write_pieces(lines, cuts)
            dst = os.path.join(self.tmp, "out%d.txt" % trial)
            n = merge_pieces.score_cbulge(dst, "#header\n", files, True, batch_rows=batch,
                                          score_fn=_batch_scorer(got_log))
            self.assertEqual(got_log, want_log, (batch, cuts, [len(x) for x in lines]))
            self.assertEqual(n, rid)
            with open(dst) as fh:
                written = fh.read().splitlines()
            self.assertEqual(written[0], "#header")
            expect = []
            for t in want:
                t = list(t)
                t.pop(-2)
                expect.append("\t".join(str(x) for x in t) + "\t0")
            self.assertEqual(written[1:], expect)

    def test_a_line_is_never_cut(self):
        """A line's rows past the batch size stay in its batch, as one process keeps them."""
        lines = [[["a%d" % j, "A", "C", 33, 0] for j in range(5)],
                 [["b%d" % j, "A", "C", 33, 0] for j in range(5)]]
        log = []
        merge_pieces.score_cbulge(os.path.join(self.tmp, "o.txt"), "#h\n",
                                  self.write_pieces(lines, [1]), True, batch_rows=3,
                                  score_fn=_batch_scorer(log))
        self.assertEqual([len(b) for b in log], [5, 5])

    def test_the_batch_is_the_scripts_one(self):
        self.assertEqual(cbulge_rows.BATCH_ROWS, 100000)
        for script in ("new_simple_analysis.py", "analisi_indels_NNN.py"):
            with open(os.path.join(HERE, script)) as fh:
                src = fh.read()
            self.assertIn("if len(cluster_to_save) >= _cbr.BATCH_ROWS:", src, script)
            self.assertIn("return _cbr.score_rows(cluster_targets, do_scores)", src, script)


class TheScorersBatchLoop(unittest.TestCase):
    """cbulge_rows.score_rows nulls the rows the model cannot score (two bulges, an N, ...):
    in a dense locus nearly every row of a batch, whose 100,000 rows it then read a list of
    for each row (about a minute a batch, which the merge of the pieces paid for every batch
    of a chromosome in one process). The same scores, in a fraction of a second."""

    def test_a_batch_of_nulled_rows_is_scored_at_once(self):
        import types
        old = sys.modules.get("scorer_runner")
        sys.modules["scorer_runner"] = types.SimpleNamespace(
            CRISPR_BULGE_predict_list=lambda sgs, offs: [0.25] * len(sgs))
        self.addCleanup(lambda: sys.modules.pop("scorer_runner", None) if old is None
                        else sys.modules.__setitem__("scorer_runner", old))
        rng = random.Random(5)
        rows, want = [], []
        for i in range(100000):
            if rng.random() < 0.9:  # a DNA and an RNA bulge: nulled
                rows.append(["X", "ACGTACGTAC-GTACGTACGTNNN", "ACGTACG-ACGTACGTACGTACGG",
                             "n", 33, 0])
                want.append("-1")
            else:
                rows.append(["X", "ACGTACGTACGTACGTACGTNNN", "ACGTACGTACGTACGTACGAAGG",
                             "n", 33, 0])
                want.append("0.25")
        t0 = time.time()
        got = cbulge_rows.score_rows(rows, True)
        took = time.time() - t0
        self.assertEqual([r[-1] for r in got], want)
        self.assertEqual([r[-3] for r in got], want)
        self.assertLess(took, 10.0)


class MergedFiles(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="pieces_merge_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def put(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_the_dense_window_bed_keeps_the_line_read_first(self):
        head = "#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\tiupac_protospacer\t" \
               "bits\tfour_allele\tgreedy\n"
        a = self.put("a.bed", head
                     + "chrT\t10\t33\tG\t3\tS1\tACRtYW\t3.00\t0\t1\n"
                     + "chrT\t50\t73\tG\t3\tS1\tAAA\t3.00\t0\t0\n")
        b = self.put("b.bed", head
                     + "chrT\t10\t33\tG\t3\tS1\tACRTYW\t3.00\t0\t1\n"
                     + "chrT\t9\t33\tG\t3\tS1\tACRTYW\t3.00\t0\t0\n")
        out = os.path.join(self.tmp, "m.bed")
        merge_pieces.merge_hvdr(out, [a, b])
        with open(out) as fh:
            got = fh.read()
        # the upper-case letters sort first; sorted by start, end, guide
        self.assertEqual(got, head
                         + "chrT\t9\t33\tG\t3\tS1\tACRTYW\t3.00\t0\t0\n"
                         + "chrT\t10\t33\tG\t3\tS1\tACRTYW\t3.00\t0\t1\n"
                         + "chrT\t50\t73\tG\t3\tS1\tAAA\t3.00\t0\t0\n")

    def test_companions_keep_the_first_row_of_an_identity(self):
        head = "# a note\n#Chromosome\tPosition\tDirection\tcrRNA\tDNA\tSNP\tPhase_Confirmation\n"
        a = self.put("a.tsv", head + "chrT\t1\t+\tG\tD\ts1\tCONFIRMED\n")
        b = self.put("b.tsv", head + "chrT\t1\t+\tG\tD\ts1\tPUTATIVE\n"
                     + "chrT\t2\t+\tG\tD\ts1\tPUTATIVE\n")
        out = os.path.join(self.tmp, "m.tsv")
        merge_pieces.merge_companion(out, [a, b], ("Chromosome", "Position", "Direction",
                                                    "crRNA", "DNA", "SNP"))
        with open(out) as fh:
            self.assertEqual(fh.read(), head + "chrT\t1\t+\tG\tD\ts1\tCONFIRMED\n"
                             + "chrT\t2\t+\tG\tD\ts1\tPUTATIVE\n")

    def test_four_allele_union_and_graph_sum(self):
        a = self.put("a.q", "#chrom\tstart\tend\nchrT\t9\t10\nchrT\t3\t4\n")
        b = self.put("b.q", "#chrom\tstart\tend\nchrT\t3\t4\nchrT\t100\t101\n")
        out = os.path.join(self.tmp, "m.q")
        merge_pieces.merge_four_allele(out, [a, b])
        with open(out) as fh:
            self.assertEqual(fh.read(), "#chrom\tstart\tend\nchrT\t3\t4\nchrT\t9\t10\n"
                             "chrT\t100\t101\n")
        g1 = self.put("g1", "ref\tvar\n0\t0\n1\t2\n")
        g2 = self.put("g2", "ref\tvar\n0\t5\n1\t0\n")
        merge_pieces.merge_cfd_graph(out, [g1, g2])
        with open(out) as fh:
            self.assertEqual(fh.read(), "ref\tvar\n0\t5\n2\t2\n")


class ThePool(unittest.TestCase):
    def test_pieces_and_ranges(self):
        self.assertEqual(postproc_pool.n_pieces(0, 100, 8), 1)
        self.assertEqual(postproc_pool.n_pieces(150, 100, 8), 1)
        self.assertEqual(postproc_pool.n_pieces(450, 100, 8), 4)
        self.assertEqual(postproc_pool.n_pieces(10 ** 6, 100, 8), 8)
        for lines in (1, 7, 100, 101, 1999):
            for n in (1, 2, 3, 7):
                if n > lines:
                    continue
                r = postproc_pool.ranges(lines, n)
                self.assertEqual(len(r), n)
                self.assertEqual(r[0][0], 1)
                self.assertEqual(r[-1][1], lines)
                for (a, b), (c, _d) in zip(r, r[1:]):
                    self.assertEqual(c, b + 1)
                sizes = [b - a + 1 for a, b in r]
                self.assertLessEqual(max(sizes) - min(sizes), 1)

    def test_the_switch(self):
        old = os.environ.get("CRISPRME_POSTPROC_PIECES")
        try:
            for v, on in (("0", False), ("no", False), ("1", True), ("", True)):
                os.environ["CRISPRME_POSTPROC_PIECES"] = v
                self.assertEqual(postproc_pool.pieces_on(), on, v)
        finally:
            if old is None:
                os.environ.pop("CRISPRME_POSTPROC_PIECES", None)
            else:
                os.environ["CRISPRME_POSTPROC_PIECES"] = old

    def test_tasks_queue_tasks_and_the_first_error_stops_the_pool(self):
        pool = postproc_pool.TaskPool(3)
        done = []
        lock = threading.Lock()

        def leaf(i):
            with lock:
                done.append(i)

        def parent(i):
            for j in range(4):
                pool.put(0, leaf, (i, j))

        for i in range(5):
            pool.put(1, parent, i)
        pool.run()
        self.assertEqual(sorted(done), sorted((i, j) for i in range(5) for j in range(4)))

        pool = postproc_pool.TaskPool(2)
        ran = []

        def boom():
            raise postproc_pool.StageError("piece failed")

        pool.put(0, boom)
        for i in range(10):
            pool.put(1, ran.append, i)
        with self.assertRaises(postproc_pool.StageError):
            pool.run()
        self.assertLess(len(ran), 10)

    def test_the_log_keeps_one_pieces_output_and_what_the_others_add(self):
        """Each piece starts up as one process does: the log keeps the first piece's output
        and the lines another piece adds (numbers aside), never a piece's own files or its
        run time."""
        tmp = tempfile.mkdtemp(prefix="pieces_log_")
        self.addCleanup(shutil.rmtree, tmp, True)
        pdir = os.path.join(tmp, "out.pieces")
        os.makedirs(pdir)
        logs = []
        for i, text in enumerate((
                "Loaded 12 SNP dictionary entr(y/ies) for chrT\n"
                "Read 3 indel(s) of the phased indel genotype store in 0.4s\n"
                "Wrote phase-confirmation companion (2 variant off-target row[s]) to "
                + os.path.join(pdir, "p0000.phase_confirmation.tsv") + "\n"
                "ANALYSIS COMPLETE IN 3.2\n",
                "Loaded 12 SNP dictionary entr(y/ies) for chrT\n"
                "Read 3 indel(s) of the phased indel genotype store in 0.7s\n"
                "WARNING [mega-indel-af]: emission error (further ones not shown)\n"
                "ANALYSIS COMPLETE IN 4.0\n",
                "Read 3 indel(s) of the phased indel genotype store in 0.5s\n"
                "WARNING [mega-indel-af]: emission error (further ones not shown)\n"
                "population-summary companion: 2 off-target row[s] on chrT could not be "
                "summarized (e.g. x)\n")):
            logs.append(os.path.join(pdir, "p%04d.log" % i))
            with open(logs[-1], "w") as fh:
                fh.write(text)
        self.assertEqual(postproc_pool.piece_notes(logs + [logs[0] + ".absent"], pdir), [
            "Loaded 12 SNP dictionary entr(y/ies) for chrT",
            "Read 3 indel(s) of the phased indel genotype store in 0.4s",
            "WARNING [mega-indel-af]: emission error (further ones not shown)",
            "population-summary companion: 2 off-target row[s] on chrT could not be "
            "summarized (e.g. x)"])

    def test_a_piece_writes_its_output_in_its_log_and_a_failure_shows_it(self):
        tmp = tempfile.mkdtemp(prefix="pieces_call_")
        self.addCleanup(shutil.rmtree, tmp, True)
        log = os.path.join(tmp, "p0000.log")
        ok = subprocess.run(
            [sys.executable, "-c", "import sys, postproc_pool; postproc_pool._call_logged("
             "[sys.executable, '-c', 'print(\"start-up\")'], sys.argv[1], what='piece 1')",
             log], cwd=HERE, capture_output=True, text=True)
        self.assertEqual((ok.returncode, ok.stdout, ok.stderr), (0, "", ""))
        with open(log) as fh:
            self.assertEqual(fh.read(), "start-up\n")
        with self.assertRaises(postproc_pool.StageError) as err:
            postproc_pool._call_logged(
                [sys.executable, "-c", "print('\\n'.join('line %d' % i for i in range(40)));"
                 " raise SystemExit(3)"], log, what="piece 2 of the post-analysis of chrT")
        msg = str(err.exception)
        self.assertIn("piece 2 of the post-analysis of chrT failed (exit 3)", msg)
        self.assertIn(log, msg)
        self.assertIn("line 39", msg)
        self.assertNotIn("line 19\n", msg)


def _gnu_wc():
    try:
        return subprocess.run(["wc", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


def _gnu_sed():
    """The scripts rewrite their best files' header with GNU ``sed -i`` (the pipeline
    image's): a BSD sed says so on stderr, and leaves the header as it was."""
    try:
        return subprocess.run(["sed", "--version"], capture_output=True).returncode == 0
    except OSError:
        return False


def _files(prefix, names):
    out = {}
    for name in names:
        path = prefix + name
        if os.path.exists(path):
            with open(path, "rb") as fh:
                out[name] = fh.read()
    return out


class PiecesThroughTheRealScripts(unittest.TestCase):
    """new_simple_analysis.py on a chromosome of dense and sparse windows, as one process
    and as pieces merged: the same files."""

    SNP_FILES = (".bestCFD.txt", ".bestmmblg.txt", ".bestCRISPR_BULGE.txt", ".CFDGraph.txt",
                 ".high_variant_density_regions.bed", ".four_allele_positions.bed",
                 ".population_summary.tsv", ".phase_confirmation.tsv", ".snp_snp_cooc.tsv")

    def setUp(self):
        from test_window_density import _PostAnalysis  # noqa: F401  (its constants below)
        import test_window_density as twd
        self.twd = twd
        self.tmp = tempfile.mkdtemp(prefix="pieces_snp_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        rng = random.Random(3)
        guide = twd.GUIDE
        sites, rows, entries = [], [], {}
        for k in range(14):
            site = list(guide + "TGG")
            # mismatches at three columns, and variant columns some samples carry
            for col in rng.sample(range(20), 3):
                site[col] = rng.choice([b for b in "ACGT" if b != guide[col]])
            site = "".join(site)
            sites.append(site)
            cols = rng.sample(range(23), rng.randint(2, 11))
            codes = {}
            for col in cols:
                ref = site[col]
                alt = guide[col] if col < 20 and site[col] != guide[col] else \
                    rng.choice([b for b in "ACGT" if b != ref])
                code = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
                        frozenset("AT"): "W", frozenset("GT"): "K",
                        frozenset("AC"): "M"}[frozenset((ref, alt))]
                codes[col] = code
                carriers = [("S%d" % s, rng.choice(["1|0", "0|1", "1|1"]))
                            for s in range(6) if rng.random() < 0.3]
                entries[(k, col)] = "%s;%s,%s;rs%d_%d;0.1" % (
                    ",".join("%s:%s" % c for c in carriers), ref, alt, k, col)
            window = twd._site_window(site, codes)
            rows.append((k, window, "+"))
            # a second alignment of the window, one letter in lower case
            low = list(window)
            low[0] = low[0].lower()
            rows.append((k, "".join(low), "+"))
        self.sites, self.rows, self.entries = sites, rows, entries

    def write(self, d):
        os.makedirs(d)
        pam = os.path.join(d, "pam.txt")
        with open(pam, "w") as fh:
            fh.write(self.twd.PAM_NRG)
        genome, starts = self.twd.FLANK, []
        for site in self.sites:
            starts.append(len(genome))
            genome += site + self.twd.FLANK
        fasta = os.path.join(d, "chrT.fa")
        with open(fasta, "w") as fh:
            fh.write(">chrT\n" + genome + "\n")
        dct = os.path.join(d, "my_dict_chrT.json")
        with open(dct, "w") as fh:
            json.dump({"chrT,%d" % (starts[i] + col + 1): v
                       for (i, col), v in self.entries.items()}, fh)
        targets = os.path.join(d, "chrT.targets.txt")
        with open(targets, "w") as fh:
            fh.write("#header\n")
            for i, window, strand in self.rows:
                fh.write("\t".join(
                    ["X", self.twd.CRRNA, window, "chrT", str(starts[i]), str(starts[i]),
                     strand, "0", "0", "0", "NA", "y", "NA", "NA", "NA", "NA", "NA",
                     "NA"]) + "\n")
        return fasta, dct, targets, pam

    def run_script(self, args, env):
        e = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
        # the dictionary path joins a row's carriers from a set, in the order of the string
        # hashes, which every process seeds anew (remove_contiguous_samples.py sorts them
        # downstream): one seed for every process here, so that one process and the pieces
        # can be compared byte for byte
        e.update(CRISPRME_SKIP_SCORER_ENV="1", PYTHONHASHSEED="0", **env)
        done = subprocess.run([sys.executable, os.path.join(HERE, "new_simple_analysis.py")]
                              + args, cwd=HERE, env=e, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout[-3000:] + done.stderr[-3000:])
        if _gnu_sed():
            self.assertEqual(done.stderr, "")
        return done.stdout

    def one_and_pieces(self, n_pieces, **env):
        d = os.path.join(self.tmp, "run_%d_%s" % (n_pieces, "_".join(sorted(env.values()))))
        fasta, dct, targets, pam = self.write(d)
        whole = os.path.join(d, "whole")
        whole_out = self.run_script([fasta, dct, targets, pam, whole, "4"], env)
        n_lines = len(self.rows)
        prefixes, piece_outs = [], []
        prep = os.path.join(d, "dict.pkl")
        self.run_script([fasta, dct, targets, pam, os.path.join(d, "prep"), "4"],
                        dict(env, CRISPRME_POSTPROC_PREPARE=prep))
        self.assertTrue(os.path.exists(prep))
        for i, (a, b) in enumerate(postproc_pool.ranges(n_lines, n_pieces)):
            p = os.path.join(d, "p%d" % i)
            prefixes.append(p)
            piece_outs.append(self.run_script(
                [fasta, dct, targets, pam, p, "4"],
                dict(env, CRISPRME_POSTPROC_PIECE="1", CRISPRME_POSTPROC_RANGE="%d:%d" % (a, b),
                     CRISPRME_POSTPROC_DICT=prep)))
        # the merged files, named as one process's (the notes name them)
        os.makedirs(os.path.join(d, "merged"))
        merged = os.path.join(d, "merged", "whole")
        merge_pieces.merge_snp(merged, prefixes)
        # what the runs said, for the notes
        self.said = {"whole": whole_out, "pieces": piece_outs, "merged": merged,
                     "prefixes": prefixes}
        return _files(whole, self.SNP_FILES), _files(merged, self.SNP_FILES)

    def test_per_sample_pieces_are_one_process(self):
        for n in (2, 5):
            for pre in ("1", "0"):
                with self.subTest(pieces=n, precheck=pre):
                    whole, merged = self.one_and_pieces(
                        n, CRISPRME_FAST_MODE="0", CRISPRME_HAPLOTYPE_PRECHECK=pre)
                    self.assertTrue(whole[".bestCFD.txt"].count(b"\n") > 2)
                    self.assertEqual(sorted(whole), sorted(merged))
                    for name in whole:
                        self.assertEqual(merged[name], whole[name], name)

    def test_population_level_pieces_are_one_process(self):
        whole, merged = self.one_and_pieces(3, CRISPRME_FAST_MODE="1")
        self.assertEqual(sorted(whole), sorted(merged))
        for name in whole:
            self.assertEqual(merged[name], whole[name], name)

    def test_precheck_on_and_off_give_the_same_files(self):
        on, _m = self.one_and_pieces(1, CRISPRME_FAST_MODE="0", CRISPRME_HAPLOTYPE_PRECHECK="1")
        off, _m = self.one_and_pieces(1, CRISPRME_FAST_MODE="0",
                                      CRISPRME_HAPLOTYPE_PRECHECK="0")
        self.assertEqual(sorted(on), sorted(off))
        for name in on:
            self.assertEqual(on[name], off[name], name)

    NOTES = ("Haplotype pre-check on ", "K-skip on ", "High-variant-density windows of ",
             "Four-allele positions ", "--max-degeneracy ", "snp-indel: ")

    def notes(self, text):
        return [line for line in text.splitlines() if line.startswith(self.NOTES)]

    def test_the_merge_writes_the_notes_of_one_process(self):
        """The notes a chromosome ends with (the pre-check or the K-skip, the dense windows,
        --max-degeneracy): no piece writes its own (each would count its lines only, in a
        log line that says it is the chromosome's), the merge writes one process's from the
        pieces' counts and the merged files."""
        for env in ({"CRISPRME_FAST_MODE": "0", "CRISPRME_IUPAC_CAP": "4"},
                    {"CRISPRME_FAST_MODE": "0", "CRISPRME_IUPAC_CAP": "4",
                     "CRISPRME_HAPLOTYPE_PRECHECK": "0"},
                    {"CRISPRME_FAST_MODE": "1", "CRISPRME_IUPAC_CAP": "4",
                     "CRISPRME_MAX_DEGENERACY": "1"}):
            with self.subTest(**env):
                self.one_and_pieces(3, **env)
                want = self.notes(self.said["whole"])
                self.assertTrue(
                    any(n.startswith("High-variant-density windows of chrT: ") for n in want),
                    want)
                self.assertTrue(any(n.startswith(("Haplotype pre-check on chrT: ",
                                                  "K-skip on chrT: ")) for n in want), want)
                for out in self.said["pieces"]:
                    self.assertEqual(self.notes(out), [])
                self.assertEqual(merge_pieces.snp_notes(self.said["merged"],
                                                        self.said["prefixes"]), want)


@unittest.skipUnless(_gnu_wc(), "analisi_indels_NNN.py reads GNU wc -l")
class IndelPiecesThroughTheRealScript(unittest.TestCase):
    """analisi_indels_NNN.py on test_indel_overlaid_rows' world, as one process and as
    pieces merged: the same files, with the carriers-first pre-check on and off."""

    FILES = (".bestCFD_INDEL.txt", ".bestmmblg_INDEL.txt", ".bestCRISPR_BULGE_INDEL.txt",
             ".indel_snp_cooc.tsv", ".indel_af.tsv")

    def test_pieces_are_one_process(self):
        import test_indel_overlaid_rows as tio
        for pre in ("1", "0"):
            with self.subTest(precheck=pre):
                tmp = tempfile.mkdtemp(prefix="pieces_indel_")
                self.addCleanup(shutil.rmtree, tmp, True)
                env = {"CRISPRME_HAPLOTYPE_PRECHECK": pre, "PYTHONHASHSEED": "0"}
                fx = tio._Fixture(tmp, tio.World(), env=env)
                self.assertEqual(fx.proc.returncode, 0, fx.proc.stderr[-2000:])
                whole = _files(fx.out, self.FILES)
                self.assertTrue(whole[".bestCFD_INDEL.txt"].count(b"\n") > 2)
                logd = os.path.join(tmp, "Dictionaries", "log_indels_" + tio.VCF)
                cluster = os.path.join(tmp, "targets.cluster.txt")
                n_lines = len(fx.rows_in)
                prefixes = []
                run_env = {k: v for k, v in os.environ.items()
                           if k not in ("MAMBA_EXE", "CONDA_EXE", "CRISPRME_CONDA_EXE")}
                run_env["PATH"] = os.path.dirname(sys.executable) + ":/usr/bin:/bin"
                run_env.update(env, PYTHONHASHSEED="0")
                for i, (a, b) in enumerate(postproc_pool.ranges(n_lines, 3)):
                    p = os.path.join(tmp, "out", "piece%d" % i)
                    prefixes.append(p)
                    argv = [sys.executable, tio.ANALYZER, "", cluster, p, logd,
                            os.path.join(tmp, "pam.txt"), "4",
                            os.path.join(tmp, tio.CHROM + ".fa"),
                            os.path.join(tmp, "guides.txt"), "1", "1"]
                    e = dict(run_env, CRISPRME_POSTPROC_PIECE="1",
                             CRISPRME_POSTPROC_RANGE="%d:%d" % (a, b))
                    done = subprocess.run(argv, cwd=HERE, env=e, capture_output=True,
                                          text=True, timeout=600)
                    self.assertEqual(done.returncode, 0, done.stderr[-2000:])
                    self.assertEqual(done.stderr, "")
                    # the notes the chromosome ends with are the merge's
                    self.assertEqual(self.notes(done.stdout), [])
                merged = os.path.join(tmp, "out", "merged")
                merge_pieces.merge_indel(merged, prefixes)
                got = _files(merged, self.FILES)
                self.assertEqual(sorted(whole), sorted(got))
                for name in whole:
                    self.assertEqual(got[name], whole[name], name)
                want = self.notes(fx.proc.stdout)
                self.assertTrue(any(" overlaid target(s) resolved per carrier haplotype "
                                    in n for n in want), want)
                self.assertEqual(merge_pieces.indel_notes(prefixes), want)

    @staticmethod
    def notes(text):
        return [line for line in text.splitlines()
                if line.startswith(("indel-snp: ", "WARNING [indel-snp]: ", "WARNING: "))]


if __name__ == "__main__":
    unittest.main()

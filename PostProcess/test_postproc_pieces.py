"""Option P (2.10.0): a chromosome's post-analysis in pieces over one pool of workers,
merged into the files one process writes, byte for byte.

  * merge_pieces.score_cbulge scores the pieces' rows in exactly the batches one process
    forms (every cbulge_rows.BATCH_ROWS rows, cut after the input line that reaches it): the
    model's floats depend on which rows share a batch, so other batches would change the
    last printed digit of some CRISPR-Bulge scores (the cost judges measured it). Checked
    with a scorer whose output depends on the whole batch, against the loop of one process;
  * cbulge_rows.score_rows scores a batch whose rows it nearly all nulls at once (it read a
    list of them for each row);
  * the companions keep the first row of each identity, as one process does, merge_snp too
    where two pieces hold one; the dense-window BED keeps one line a window, the one
    resultIntegrator reads first, whichever piece lists it first; the four-allele positions
    are the sorted union; the CFD graphs are summed;
  * postproc_pool cuts a chromosome into contiguous ranges, runs the stages of every
    chromosome on one pool, and stops at the first failure; a piece's output goes to its own
    file, and the log keeps one piece's start-up and what the others add;
  * the SNP pieces' copy of the dictionary keys sits in the run's output folder, which its
    group can write under a umask of 002: it is JSON (it was a pickle, and loading one runs
    what it names), made 0600 through a name nothing stood under, in a pieces' folder made
    0700; a pickle planted where a piece reads its copy runs nothing, and stops the piece;
  * through the real scripts: new_simple_analysis.py and analisi_indels_NNN.py run as one
    process and as pieces (CRISPRME_POSTPROC_RANGE, CRISPRME_POSTPROC_PIECE) then merged give
    the same files, byte for byte, with the pre-check on and off (the pieces skip and keep
    bulged raw targets too), and analisi_indels_NNN.py the same files with its carriers-first
    pre-check on as off, a target at exactly the budget among them; the pieces write none of
    the notes a chromosome ends with, the merge writes one process's, which say when the
    pre-check is off.

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
        # whichever piece lists it first (review r3 of lane S, its mutant D4: keeping a
        # window's last line passed when the smallest came last)
        for order in ((a, b), (b, a)):
            with self.subTest(order=[os.path.basename(p) for p in order]):
                merge_pieces.merge_hvdr(out, list(order))
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

    def test_merge_snp_keeps_the_first_row_of_an_identity_across_pieces(self):
        """merge_snp merges each SNP companion keeping the first row of each identity, as one
        process writes them: a window whose alignments a cut put in two pieces gives each a
        row of the same identity. (Review r3 of lane S, its mutant D3: merge_snp merging them
        with no identity kept both, and no test saw it, as the real-script fixtures write no
        companion: they need the dict-less tiers.)"""
        pieces = []
        for i in range(2):
            p = os.path.join(self.tmp, "p%d" % i)
            pieces.append(p)
            for name in (".bestCFD.txt", ".bestmmblg.txt"):
                self.put("p%d%s" % (i, name), "#Bulge_type\tcrRNA\tDNA\n")
            for name, identity in merge_pieces._SNP_COMPANIONS.items():
                rows = ["\t".join(["shared"] * len(identity) + ["piece %d" % i]),
                        "\t".join(["own %d" % i] * len(identity) + ["piece %d" % i])]
                self.put("p%d.%s" % (i, name), "# a note\n#" + "\t".join(identity + ("Note",))
                         + "\n" + "\n".join(rows) + "\n")
        merge_pieces.merge_snp(os.path.join(self.tmp, "m"), pieces)
        for name, identity in merge_pieces._SNP_COMPANIONS.items():
            with self.subTest(companion=name):
                with open(os.path.join(self.tmp, "m." + name)) as fh:
                    got = fh.read().splitlines()
                self.assertEqual(got, [
                    "# a note", "#" + "\t".join(identity + ("Note",)),
                    "\t".join(["shared"] * len(identity) + ["piece 0"]),
                    "\t".join(["own 0"] * len(identity) + ["piece 0"]),
                    "\t".join(["own 1"] * len(identity) + ["piece 1"])])

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


class _RunsCodeWhenLoaded(object):
    """What a pickle can hold: loading it calls the function it names, here os.mkdir(mark)
    (any function would do: os.system, exec)."""

    def __init__(self, mark):
        self.mark = mark

    def __reduce__(self):
        return (os.mkdir, (self.mark,))


class TheDictionaryCopy(unittest.TestCase):
    """The SNP pieces' copy of the chromosome's dictionary keys (postproc_pool.write_dict_copy,
    read_dict_copy): a file of the run's output folder, which other members of its group can
    write under a umask of 002, so it is JSON, which loads as data whoever wrote it (it was a
    pickle, which runs what it names: review r3 of lane S, S-R3-1), and it is made 0600
    through a name nothing stood under."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="pieces_copy_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.path = os.path.join(self.tmp, postproc_pool.DICT_COPY_NAME)

    def test_round_trip_json_0600(self):
        entries = {"chrT,5": "S1:1|0;A,G;rs1;0.1$S2:0|1;A,C;rs2;0.2", "chrT,9": "\u00e9;x"}
        old = os.umask(0o002)
        try:
            postproc_pool.write_dict_copy(self.path, 1, entries)
        finally:
            os.umask(old)
        self.assertEqual(os.stat(self.path).st_mode & 0o777, 0o600)
        with open(self.path) as fh:
            self.assertEqual(json.load(fh), {"haplo": True, "dict": entries})
        haplo, got = postproc_pool.read_dict_copy(self.path)
        self.assertIs(haplo, True)
        self.assertEqual(got, entries)
        self.assertEqual(list(got), list(entries))  # in the order of the dictionary
        self.assertEqual(os.listdir(self.tmp), [postproc_pool.DICT_COPY_NAME])

    def test_entries_json_cannot_hold_exactly_make_no_copy(self):
        import decimal
        for bad in ({"chrT,5": decimal.Decimal("0.1")}, {"chrT,5": float("nan")}):
            with self.subTest(bad=bad):
                with self.assertRaises((TypeError, ValueError)):
                    postproc_pool.write_dict_copy(self.path, False, bad)
                self.assertEqual(os.listdir(self.tmp), [])

    def test_never_through_a_name_already_there(self):
        """A link put where the copy is first written (another user's, to a file of the
        user's): the write fails, and the file it names is left as it was."""
        victim = os.path.join(self.tmp, "victim.txt")
        with open(victim, "w") as fh:
            fh.write("mine\n")
        os.symlink(victim, self.path + ".tmp")
        with self.assertRaises(FileExistsError):
            postproc_pool.write_dict_copy(self.path, True, {"chrT,5": "x"})
        with open(victim) as fh:
            self.assertEqual(fh.read(), "mine\n")
        self.assertFalse(os.path.exists(self.path))

    def test_a_file_that_is_no_copy_is_refused_and_runs_nothing(self):
        mark = os.path.join(self.tmp, "ran")
        with open(self.path, "wb") as fh:
            import pickle
            pickle.dump(_RunsCodeWhenLoaded(mark), fh, protocol=0)  # ASCII, so UTF-8 too
        for text in (None, "[1, 2]", '{"haplo": 1, "dict": {}}', '{"haplo": true}'):
            with self.subTest(text=text):
                if text is not None:
                    with open(self.path, "w") as fh:
                        fh.write(text)
                with self.assertRaises(ValueError) as err:
                    postproc_pool.read_dict_copy(self.path)
                self.assertIn(self.path, str(err.exception))
                self.assertFalse(os.path.exists(mark))


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
        """CRISPRME_POSTPROC_PIECES and CRISPRME_HAPLOTYPE_PRECHECK read alike: 0, false, no
        and off, in any case and with spaces around, turn each off; anything else leaves it
        on (haplotype_precheck.switch_on)."""
        import haplotype_precheck
        values = (("0", False), ("no", False), ("off", False), ("OFF", False),
                  (" Off ", False), ("False", False), ("FALSE", False), ("No", False),
                  ("1", True), ("", True), ("on", True), ("yes", True), ("2", True))
        old = os.environ.get("CRISPRME_POSTPROC_PIECES")
        try:
            for v, on in values:
                os.environ["CRISPRME_POSTPROC_PIECES"] = v
                self.assertEqual(postproc_pool.pieces_on(), on, v)
                self.assertEqual(haplotype_precheck.switch_on(
                    {"CRISPRME_HAPLOTYPE_PRECHECK": v}, "CRISPRME_HAPLOTYPE_PRECHECK"), on, v)
        finally:
            if old is None:
                os.environ.pop("CRISPRME_POSTPROC_PIECES", None)
            else:
                os.environ["CRISPRME_POSTPROC_PIECES"] = old
        self.assertTrue(haplotype_precheck.switch_on({}, "CRISPRME_HAPLOTYPE_PRECHECK"))

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
        run times (the indel script's "Done <chrom> <s>" too: the chromosome's is the
        pool's line)."""
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
                "Done chrT 3.1923\n"
                "ANALYSIS COMPLETE IN 3.2\n",
                "Loaded 12 SNP dictionary entr(y/ies) for chrT\n"
                "Read 3 indel(s) of the phased indel genotype store in 0.7s\n"
                "WARNING [mega-indel-af]: emission error (further ones not shown)\n"
                "Done chrT 3.9\n"
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


_FAKE_DRIVER = r"""#!/bin/bash
# a chromosome's shell driver for postproc_pool.run_stage: "cluster" writes a clustered file
# of $LINES lines and the analysis' arguments to $CRISPRME_POSTPROC_MANIFEST, as
# scriptAnalisiNNN_v3.sh does; "finish" leaves a mark
chrom=$1 stage=$2 d=$FAKE_DIR
echo "$stage $chrom ${CRISPRME_POSTPROC_MANIFEST:-}" >> "$d/driver.log"
if [ "$stage" = cluster ]; then
    { echo "#header"; seq 1 "$LINES"; } > "$d/$chrom.cluster.txt"
    printf '%s\n' genome dict "$d/$chrom.cluster.txt" pam "$d/out_$chrom" 4 \
        > "${CRISPRME_POSTPROC_MANIFEST:?}"
    exit 0
fi
touch "$d/finished_$chrom"
"""

_FAKE_SCRIPT = r"""#!/usr/bin/env python3
# the analysis script of a piece: its arguments and range in the run's log; on stdout what
# one process says, its start-up and its run times. The pool's start-up of a SNP stage
# (CRISPRME_POSTPROC_PREPARE) writes the dictionary copy as new_simple_analysis.py does; a
# piece given one writes in copies.log the modes of the copy and of its folder
import os, stat, sys
d = os.environ["FAKE_DIR"]
prepare = os.environ.get("CRISPRME_POSTPROC_PREPARE")
if prepare:
    sys.path.insert(0, os.environ["FAKE_HERE"])
    import postproc_pool
    postproc_pool.write_dict_copy(prepare, True, {"chrT,1": "S1:1|0;A,G;rs1;0.1"})
    raise SystemExit(0)
copy = os.environ.get("CRISPRME_POSTPROC_DICT")
if copy:
    with open(os.path.join(d, "copies.log"), "a") as fh:
        fh.write("%s %o %o\n" % (copy, stat.S_IMODE(os.stat(copy).st_mode),
                                 stat.S_IMODE(os.stat(os.path.dirname(copy)).st_mode)))
with open(os.path.join(d, "script.log"), "a") as fh:
    fh.write(" ".join(sys.argv[1:]) + " " + os.environ.get("CRISPRME_POSTPROC_RANGE", "-")
             + "\n")
print("Loaded 3 SNP dictionary entr(y/ies) for chrT")
print("Done chrT 1.25")
print("ANALYSIS COMPLETE IN 1.5")
"""


class TheNotes(unittest.TestCase):
    def test_the_switch_off_names_what_it_turns_off(self):
        """Per sample the switch turns off the pre-check; population-level with
        --max-degeneracy, the K-skip (there is no pre-check there); population-level
        without it, nothing (no note)."""
        import analysis_notes
        counts = analysis_notes.snp_counts()
        self.assertEqual(analysis_notes.snp_lines("chrT", counts, False, False, -1),
                         ["Haplotype pre-check on chrT: off (CRISPRME_HAPLOTYPE_PRECHECK=0)"])
        self.assertEqual(analysis_notes.snp_lines("chrT", counts, False, True, 4)[0],
                         "K-skip on chrT: off (CRISPRME_HAPLOTYPE_PRECHECK=0)")
        self.assertEqual(analysis_notes.snp_lines("chrT", counts, False, True, -1), [])

    def test_the_indel_side_says_its_pre_check_is_off(self):
        """The carriers-first pre-check of the indel side runs in either mode, so its notes say
        when the switch turned it off (after the count of the overlaid targets), and with it
        on how many targets it told; a piece's counts carry the switch to the merge."""
        import analysis_notes
        counts = analysis_notes.indel_counts()
        counts["overlaid"].update(targets=5, rows=2, dropped=3)
        off = "indel-snp: carriers-first pre-check off (CRISPRME_HAPLOTYPE_PRECHECK=0)"
        self.assertEqual(analysis_notes.indel_lines("chrT", counts, False)[-1], off)
        self.assertNotIn(off, analysis_notes.indel_lines("chrT", counts))
        counts["overlaid"]["prechecked"] = 3
        self.assertEqual(analysis_notes.indel_lines("chrT", counts)[-1],
                         "indel-snp: 3 of them told by the carriers-first pre-check, before the "
                         "second indels of the phased store were read")
        # no overlaid target, nothing to say
        self.assertEqual(analysis_notes.indel_lines("chrT", analysis_notes.indel_counts(),
                                                    False), [])
        tmp = tempfile.mkdtemp(prefix="indel_notes_")
        self.addCleanup(shutil.rmtree, tmp, True)
        counts["overlaid"]["prechecked"] = 0
        p = os.path.join(tmp, "p0")
        with open(p + ".piece.json", "w") as fh:
            json.dump({"do_scores": False, "lines": 5, "chrom": "chrT",
                       "notes": {"counts": counts, "precheck_on": False}}, fh)
        self.assertEqual(merge_pieces.indel_notes([p])[-1], off)


class TheStageThroughAFakeDriver(unittest.TestCase):
    """postproc_pool.run_stage over a fake driver and analysis script (the merge left out):
    the analysis' arguments come from a file only the run can write, and the log says the
    chromosome's run time once."""

    def setUp(self):
        import tempfile as _tf
        self.tmp = _tf.mkdtemp(prefix="pieces_stage_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.d = os.path.join(self.tmp, "run")
        self.tmpdir = os.path.join(self.tmp, "tmpdir")
        os.makedirs(self.d)
        os.makedirs(self.tmpdir)
        self.driver = os.path.join(self.tmp, "driver.sh")
        self.script = os.path.join(self.tmp, "script.py")
        for path, text in ((self.driver, _FAKE_DRIVER), (self.script, _FAKE_SCRIPT)):
            with open(path, "w") as fh:
                fh.write(text)
            os.chmod(path, 0o755)
        saved = {k: os.environ.get(k) for k in ("TMPDIR", "FAKE_DIR", "FAKE_HERE", "LINES",
                                                 "CRISPRME_POSTPROC_PIECES")}
        os.environ.update(TMPDIR=self.tmpdir, FAKE_DIR=self.d, FAKE_HERE=HERE, LINES="90",
                          CRISPRME_POSTPROC_PIECES="1")
        saved_tempdir = _tf.tempdir
        _tf.tempdir = None  # read TMPDIR again

        def restore():
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            _tf.tempdir = saved_tempdir
        self.addCleanup(restore)
        # the merge (merge_pieces.py, which needs a real piece's files) is left out
        real_call = postproc_pool._call
        self.merged = []

        def call(cmd, env=None, what=""):
            if any(str(c).endswith("merge_pieces.py") for c in cmd):
                self.merged.append(cmd)
                return
            real_call(cmd, env=env, what=what)
        postproc_pool._call = call
        self.addCleanup(setattr, postproc_pool, "_call", real_call)

    def run_stage(self, prepare_dict=False):
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            postproc_pool.run_stage("snp", ["chrT"], 3, lambda c, st: [self.driver, c, st],
                                    self.script, prepare_dict=prepare_dict, min_lines=30)
        return out.getvalue()

    def read(self, name):
        with open(os.path.join(self.d, name)) as fh:
            return fh.read().splitlines()

    def test_the_manifest_is_a_file_of_a_folder_only_the_run_writes(self):
        """The manifest was once $TMPDIR/crisprme_postproc_<kind>_<pid>_<chrom>.args:
        another user of a shared /tmp who made that file first, read-only to the run, had the
        driver's write fail and the run read that user's list. The name is made here as it
        was, with another list in it: the run reads its own."""
        planted = os.path.join(self.tmpdir, "crisprme_postproc_snp_%d_chrT.args" % os.getpid())
        with open(planted, "w") as fh:
            fh.write("genome\ndict\n%s\npam\n%s\n4\n" % (
                os.path.join(self.d, "chrT.cluster.txt"), os.path.join(self.d, "planted")))
        os.chmod(planted, 0o444)
        self.addCleanup(os.chmod, planted, 0o644)
        self.run_stage()
        ran = self.read("script.log")
        self.assertEqual(len(ran), 3, ran)
        for line in ran:
            self.assertIn(os.path.join(self.d, "out_chrT.pieces", "p000"), line)
            self.assertNotIn("planted", line)
        manifest = self.read("driver.log")[0].split()[2]
        folder = os.path.dirname(manifest)
        self.assertNotEqual(os.path.realpath(folder), os.path.realpath(self.tmpdir))
        self.assertEqual(os.path.realpath(os.path.dirname(folder)),
                         os.path.realpath(self.tmpdir))
        self.assertFalse(os.path.exists(folder))  # the run took its folder away
        self.assertEqual(sorted(os.listdir(self.tmpdir)), [os.path.basename(planted)])

    def test_the_folder_is_the_runs_alone(self):
        seen = []
        real = postproc_pool.subprocess.call

        def call(cmd, env=None, **kw):
            m = (env or {}).get("CRISPRME_POSTPROC_MANIFEST")
            if m:
                seen.append(os.stat(os.path.dirname(m)).st_mode & 0o777)
            return real(cmd, env=env, **kw)
        postproc_pool.subprocess.call = call
        try:
            self.run_stage()
        finally:
            postproc_pool.subprocess.call = real
        self.assertEqual(seen, [0o700])

    def test_the_pieces_folder_and_the_dictionary_copy_are_the_users_alone(self):
        """The pieces' folder sits in the run's output folder, which takes the umask (002 on
        shared nodes: its group can write it). The pool makes it 0700, and the copy of the
        dictionary its pieces read is made 0600 there, as JSON (review r3 of lane S, S-R3-1:
        group-writable, and a pickle, which ran what it named)."""
        with open(os.path.join(self.d, "dict"), "w") as fh:  # the chromosome's dictionary
            fh.write("{}")
        cwd = os.getcwd()
        os.chdir(self.d)  # the manifest names the dictionary "dict", as the driver writes it
        self.addCleanup(os.chdir, cwd)
        old = os.umask(0o002)
        try:
            self.run_stage(prepare_dict=True)
        finally:
            os.umask(old)
        copies = self.read("copies.log")
        self.assertEqual(len(copies), 3, copies)
        for line in copies:
            path, mode, folder_mode = line.split()
            self.assertEqual(os.path.basename(path), "dictionary.json")
            self.assertEqual(os.path.dirname(path), os.path.join(self.d, "out_chrT.pieces"))
            self.assertEqual((mode, folder_mode), ("600", "700"))
        self.assertFalse(os.path.exists(os.path.join(self.d, "out_chrT.pieces")))

    def test_the_chromosome_says_its_own_run_time(self):
        """The pieces' "Done chrT <s>" and "ANALYSIS COMPLETE IN <s>" are theirs: the log
        keeps neither, and the pool says the chromosome's, once its pieces are merged."""
        said = self.run_stage().splitlines()
        self.assertEqual(len(self.merged), 1)
        self.assertIn("Post-analysis SNP of chrT: 90 line(s) in 3 piece(s)", said)
        self.assertIn("Loaded 3 SNP dictionary entr(y/ies) for chrT", said)
        done = [line for line in said if line.startswith("Post-analysis SNP of chrT: done in ")]
        self.assertEqual(len(done), 1, said)
        self.assertTrue(done[0].endswith(" s (3 pieces)"), done)
        self.assertFalse([line for line in said
                          if line.startswith(("Done ", "ANALYSIS COMPLETE IN"))], said)
        self.assertTrue(os.path.exists(os.path.join(self.d, "finished_chrT")))


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
    and as pieces merged: the same files. Some windows have more mismatches than the budget,
    some of them an RNA bulge, so that the pre-check skips raw targets of the pieces (on the
    dictionary path and the observed one) and keeps others, which bulged rows reach only
    with the bulge taken off the mismatches."""

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
        # windows of 5 or 6 mismatches (a budget of 4 needs 1 or 2 fixed), half of them 5 and
        # an RNA bulge (a guide base the site lacks: the gap counts as a mismatch, which the
        # bulge count takes off, so 1 to fix), the ALTs that fix them carried by a sample or
        # two each
        for k in range(14, 30):
            bulge = k % 2 == 0
            gap = rng.randrange(3, 17) if bulge else None
            aligned = list(guide + "TGG")
            spacer = [c for c in range(20) if c != gap]
            mism = rng.sample(spacer, 5 if bulge else 5 + k % 4 // 2)
            for col in mism:
                aligned[col] = rng.choice([b for b in "ACGT" if b != guide[col]])
            if bulge:
                aligned[gap] = "-"
            site = "".join(ch for ch in aligned if ch != "-")
            sites.append(site)
            at = {c: c - (gap is not None and c > gap) for c in range(23) if c != gap}
            fixing = rng.sample(mism, rng.randint(1, 3))
            others = rng.sample([c for c in range(23) if c != gap and c not in mism],
                                rng.randint(1, 8))
            codes = {}
            for col in fixing + others:
                ref = aligned[col]
                alt = guide[col] if col in fixing else \
                    rng.choice([b for b in "ACGT" if b != ref])
                codes[col] = {frozenset("AG"): "R", frozenset("CT"): "Y",
                              frozenset("CG"): "S", frozenset("AT"): "W",
                              frozenset("GT"): "K", frozenset("AC"): "M"}[frozenset((ref, alt))]
                carriers = [("S%d" % s, rng.choice(["1|0", "0|1", "1|1"]))
                            for s in rng.sample(range(6), rng.randint(1, 2))]
                entries[(k, at[col])] = "%s;%s,%s;rs%d_%d;0.1" % (
                    ",".join("%s:%s" % c for c in carriers), ref, alt, k, col)
            dna = "".join(codes.get(c, ch) for c, ch in enumerate(aligned))
            rows.append((k, dna, "+", "RNA" if bulge else "X", 1 if bulge else 0))
        # windows of 6 mismatches and an RNA or a DNA bulge (the gap counts as a mismatch,
        # which the bulge count takes off: a budget of 4 needs 2 fixed), on the dictionary path
        # (4 two-allele columns) and the observed one (12): the two fixing ALTs carried by two
        # samples (skipped) or both by S2 (kept, S2's row); the other ALTs by S3 to S5
        self.bulged = {}
        for kind in ("RNA", "DNA"):
            for n_cols in (4, 12):
                for together in (False, True):
                    k = len(sites)
                    gap = 4 + len(self.bulged)
                    guide_al = list(guide)
                    t = list(guide)
                    if kind == "RNA":
                        t[gap] = "-"
                    else:
                        guide_al.insert(gap, "-")
                        t.insert(gap, "A")
                    n_spacer = len(t)
                    guide_al += ["N", "N", "N"]
                    t += ["T", "G", "G"]
                    mism = rng.sample([c for c in range(n_spacer) if c != gap], 6)
                    for col in mism:
                        t[col] = rng.choice([b for b in "ACGT" if b != guide_al[col]])
                    fixing = rng.sample(mism, 2)
                    others = rng.sample([c for c in range(len(t)) if c != gap and c not in mism],
                                        n_cols - 2)
                    site_at, n = {}, 0
                    for c, ch in enumerate(t):
                        if ch != "-":
                            site_at[c] = n
                            n += 1
                    codes = {}
                    for col in fixing + others:
                        ref = t[col]
                        alt = guide_al[col] if col in fixing else \
                            rng.choice([b for b in "ACGT" if b not in (ref, guide_al[col])])
                        codes[col] = {frozenset("AG"): "R", frozenset("CT"): "Y",
                                      frozenset("CG"): "S", frozenset("AT"): "W",
                                      frozenset("GT"): "K",
                                      frozenset("AC"): "M"}[frozenset((ref, alt))]
                        if col in fixing:
                            who = ["S2"] if together else ["S%d" % fixing.index(col)]
                        else:
                            who = rng.sample(["S3", "S4", "S5"], rng.randint(1, 2))
                        entries[(k, site_at[col])] = "%s;%s,%s;rs%d_%d;0.1" % (
                            ",".join("%s:1|0" % w for w in who), ref, alt, k, col)
                    sites.append("".join(ch for ch in t if ch != "-"))
                    dna = "".join(codes.get(c, ch) for c, ch in enumerate(t))
                    rows.append((k, dna, "+", kind, 1, "".join(guide_al)))
                    self.bulged[k] = (kind, n_cols, together)
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
            for row in self.rows:
                i, window, strand = row[:3]
                kind, bulges = row[3:5] if len(row) > 3 else ("X", 0)
                crrna = row[5] if len(row) > 5 else self.twd.CRRNA
                fh.write("\t".join(
                    [kind, crrna, window, "chrT", str(starts[i]), str(starts[i]),
                     strand, "0", str(bulges), str(bulges), "NA", "y", "NA", "NA", "NA", "NA",
                     "NA", "NA"]) + "\n")
        self.starts = starts
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
        prep = os.path.join(d, "dictionary.json")
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
                    if pre == "1":
                        # the pre-check skips raw targets of the pieces, and keeps others
                        note = [line for line in self.said["whole"].splitlines()
                                if line.startswith("Haplotype pre-check on chrT: ")]
                        skipped, seen = (int(x) for x in note[0].split(": ")[1].split()[:3:2])
                        self.assertGreater(skipped, 7, note)
                        self.assertGreater(seen - skipped, 7, note)
                        # bulged rows within the budget only with their bulge taken off
                        self.assertIn(b"\nRNA\t", whole[".bestCFD.txt"])
                        self.assertIn(b"\nDNA\t", whole[".bestCFD.txt"])
                    # a bulged window needing two fixes: none where two samples carry one
                    # each (skipped, with the pre-check on), S2's row where S2 carries both
                    rows = [line.split("\t") for line in
                            whole[".bestCFD.txt"].decode().splitlines()[1:]]
                    for k, (kind, n_cols, together) in self.bulged.items():
                        carriers = {s for r in rows if r[4] == str(self.starts[k])
                                    for s in r[12].split(",") if s not in ("n", "NA")}
                        self.assertEqual(carriers, {"S2"} if together else set(),
                                         (kind, n_cols, together))

    def test_the_dictionary_copy_is_json_and_the_users_alone(self):
        """The start-up of the pieces (CRISPRME_POSTPROC_PREPARE) writes the dictionary keys of
        the chromosome's lines as JSON, the dictionary's own format, 0600 whatever the umask:
        its entries those of the dictionary."""
        d = os.path.join(self.tmp, "copy")
        fasta, dct, targets, pam = self.write(d)
        prep = os.path.join(d, "dictionary.json")
        old = os.umask(0o002)
        try:
            self.run_script([fasta, dct, targets, pam, os.path.join(d, "prep"), "4"],
                            {"CRISPRME_FAST_MODE": "0", "CRISPRME_POSTPROC_PREPARE": prep})
        finally:
            os.umask(old)
        self.assertEqual(os.stat(prep).st_mode & 0o777, 0o600)
        with open(prep) as fh:
            copy = json.load(fh)
        with open(dct) as fh:
            dictionary = json.load(fh)
        self.assertEqual(sorted(copy), ["dict", "haplo"])
        self.assertIs(copy["haplo"], True)
        self.assertTrue(copy["dict"])
        self.assertEqual(copy["dict"], {k: dictionary[k] for k in copy["dict"]})

    def test_a_piece_runs_no_code_from_a_copy_others_can_write(self):
        """Review r3 of lane S (S-R3-1): the pieces read the copy from the pieces' folder of the
        run's output folder, and under a umask of 002 the folder's group could replace it
        between the start-up and a piece's load. A pickle was loaded there, and loading one
        runs what it names, as the user who runs the search. A pickle that makes a folder when
        loaded, group-writable as such a member would leave it, where a piece reads its copy:
        the piece runs none of it, and stops, naming the file."""
        d = os.path.join(self.tmp, "planted")
        fasta, dct, targets, pam = self.write(d)
        mark = os.path.join(d, "ran_code_from_the_copy")
        import pickle
        for name, protocol in (("dictionary.json", 0), ("dictionary.pkl", 4)):
            with self.subTest(name=name, protocol=protocol):
                copy = os.path.join(d, name)
                with open(copy, "wb") as fh:
                    pickle.dump(_RunsCodeWhenLoaded(mark), fh, protocol=protocol)
                os.chmod(copy, 0o664)
                e = {k: v for k, v in os.environ.items() if not k.startswith("CRISPRME_")}
                e.update(CRISPRME_SKIP_SCORER_ENV="1", CRISPRME_FAST_MODE="0",
                         CRISPRME_POSTPROC_PIECE="1",
                         CRISPRME_POSTPROC_RANGE="1:%d" % len(self.rows),
                         CRISPRME_POSTPROC_DICT=copy)
                done = subprocess.run(
                    [sys.executable, os.path.join(HERE, "new_simple_analysis.py"), fasta, dct,
                     targets, pam, os.path.join(d, "piece"), "4"],
                    cwd=HERE, env=e, capture_output=True, text=True)
                self.assertFalse(os.path.exists(mark), "a piece ran code its copy named")
                self.assertNotEqual(done.returncode, 0, done.stdout[-2000:])
                self.assertIn("%s is not a dictionary copy of the post-analysis" % copy,
                              done.stderr)

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


def _budget_world():
    """test_indel_overlaid_rows' world and two more sites, targets no real-script test held
    before (the review of lane S, round 2): BG, a deletion whose carrier's target has exactly
    the budget's 4 mismatches, beside the code of another sample's SNP it does not need; UN, a
    deletion whose target needs the guide's base at a code that does not hold the reference
    base (never expected: no genotype covers that column, so resolve_target takes the code's
    alleles there, and the pre-check must let every carrier fix it)."""
    import test_indel_overlaid_rows as tio
    w = tio.World()
    g = tio.GUIDE
    other = {"A": "C", "C": "A", "G": "T", "T": "G"}
    bg = list(g)
    for c in (1, 6, 13, 17):
        bg[c] = other[g[c]]
    w._put(5500, "".join(bg[:10]) + "G" + "".join(bg[10:]) + "TGG")
    w._indel("BG", 5500 + 10, bg[9] + "G", bg[9], {"S01": "1|0"})
    w._snp("BG", 5500 + 4 + 1, g[4], other[g[4]], {"S02": "1|0"})
    un = list(g)
    for c in (1, 7, 12, 18):
        un[c] = other[g[c]]
    assert g[3] == "G"
    un[3] = "T"  # the reference base; the code R there holds the guide's G, not T
    w._put(5700, "".join(un[:10]) + "C" + "".join(un[10:]) + "TGG")
    w._indel("UN", 5700 + 10, un[9] + "C", un[9], {"S03": "1|0"})
    w.extra_codes[5700 + 3 + 1] = "R"
    w.ref = "".join(w.seq)
    return w


@unittest.skipUnless(_gnu_wc(), "analisi_indels_NNN.py reads GNU wc -l")
class IndelPiecesThroughTheRealScript(unittest.TestCase):
    """analisi_indels_NNN.py on test_indel_overlaid_rows' world and _budget_world's two sites,
    as one process and as pieces merged: the same files, with the carriers-first pre-check on
    and off; and the same files with it on as off, the rows of a target at exactly the budget
    and of a column no genotype covers among them, the pre-check telling some targets (site
    C's chimera), and the run saying when it is off."""

    FILES = (".bestCFD_INDEL.txt", ".bestmmblg_INDEL.txt", ".bestCRISPR_BULGE_INDEL.txt",
             ".indel_snp_cooc.tsv", ".indel_af.tsv")
    OFF = "indel-snp: carriers-first pre-check off (CRISPRME_HAPLOTYPE_PRECHECK=0)"

    def one_and_pieces(self, pre):
        """One process and three pieces merged, with the pre-check on (pre "1") or off:
        (the one process's fixture, its files)."""
        import test_indel_overlaid_rows as tio
        tmp = tempfile.mkdtemp(prefix="pieces_indel_")
        self.addCleanup(shutil.rmtree, tmp, True)
        env = {"CRISPRME_HAPLOTYPE_PRECHECK": pre, "PYTHONHASHSEED": "0"}
        fx = tio._Fixture(tmp, _budget_world(), env=env)
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
        try:
            import intervaltree  # noqa: F401
        except ImportError:
            # CI's light image has no intervaltree: the pieces read the stub the fixture
            # wrote, as its one process did (test_indel_overlaid_rows._Fixture); without it
            # each piece stopped at its import, on lane S's own tree too
            run_env["PYTHONPATH"] = (os.path.join(tmp, "stub") + os.pathsep
                                     + run_env.get("PYTHONPATH", ""))
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
        return fx, whole

    def test_pieces_are_one_process(self):
        for pre in ("1", "0"):
            with self.subTest(precheck=pre):
                fx, _whole = self.one_and_pieces(pre)
                told = [n for n in self.notes(fx.proc.stdout)
                        if n.endswith(" of them told by the carriers-first pre-check, before "
                                      "the second indels of the phased store were read")]
                if pre == "1":
                    # site C's chimera at least: S05's deletion, the SNP only S06 carries
                    self.assertEqual(len(told), 1, self.notes(fx.proc.stdout))
                    self.assertGreater(int(told[0].split()[1]), 0, told)
                    self.assertNotIn(self.OFF, self.notes(fx.proc.stdout))
                else:
                    self.assertEqual(told, [])
                    self.assertIn(self.OFF, self.notes(fx.proc.stdout))

    def test_the_pre_check_changes_no_row(self):
        """With the pre-check on and off, one process each: the same files. Among the rows,
        BG's at exactly the budget (4 mismatches, S01's slot: one fewer allowed would skip it)
        and UN's, whose column no genotype covers (S03 takes the code's G there: read as known
        with nobody carrying an ALT, the target would be skipped)."""
        on_fx, on = self.one_and_pieces("1")
        _off_fx, off = self.one_and_pieces("0")
        self.assertEqual(sorted(on), sorted(off))
        for name in on:
            self.assertEqual(on[name], off[name], name)
        rows = on_fx.rows()
        bg = [r for r in rows if r["var"].split(",")[0].startswith("chrT_5510_")]
        self.assertEqual([(r["mm"], r["samples"]) for r in bg], [(4, "S01")], bg)
        un = [r for r in rows if r["var"].split(",")[0].startswith("chrT_5710_")]
        self.assertEqual([(r["mm"], r["samples"], r["dna"][3]) for r in un], [(4, "S03", "G")],
                         un)

    @staticmethod
    def notes(text):
        return [line for line in text.splitlines()
                if line.startswith(("indel-snp: ", "WARNING [indel-snp]: ", "WARNING: "))]


if __name__ == "__main__":
    unittest.main()

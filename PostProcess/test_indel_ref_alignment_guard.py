"""The reference-alignment CFD of an indel row whose reference alignment is shorter than
its alignment (analisi_indels_NNN.py).

alignRefFromVar builds a row's reference alignment by slicing the plain reference
chromosome from the row's start, as many bases as the target holds. Since CRISPRme
2.10.0 a row whose first base is one an insertion adds starts at the insertion's
anchor (indel_snp_cis.reference_start), so the slice holds that base and is never
empty; it is cut short only when an insertion lies within one target of the end of a
contig that ends in sequence, which hg38's primary chromosomes do not (each ends in
10 kb or more of N). Up to 2.8.1 such a row started at its place on the haplotype,
past the anchor, where the slice could be empty: Python returns '' for an
out-of-range slice instead of raising.

A short alignment breaks cfd_score_bulge_aware's precondition, two alignments of one
length. With a DNA bulge in the guide it indexes the protospacer by the crRNA's
positions and raises "IndexError: string index out of range": a complete-search
stops (the traceback lands in log_error.txt, which the search treats as a failure),
and only post_analisi_indel.sh run alone is silent, analisi_indels_NNN.sh not
checking the analyzer's exit status. Without a DNA bulge the scorer scores the bases
it has, a number that means nothing. The guard reports 0: there is no reference
target to score.

So that a loud failure does not become a silent 0, the rows the guard scores 0 are
counted (short_ref_rows), and the chromosome's notes say how many in one NOTE on
stdout, as they say how many rows of another contig were dropped (foreign_rows):
stderr would fail the search. A post-analysis cut into pieces leaves each piece's
count in its <piece>.piece.json, and merge_pieces.py writes the one NOTE from their
sum. TheShippedAnalyzerNearAContigEnd runs the real analyzer on rows 2.10.0 places
where the reference ends (pandas and numpy; intervaltree is stubbed when absent; GNU
`wc -l`, as the pipeline image has).

Following test_indel_pam_creation.py and test_indel_snp_overlay.py, the module's
top level reads sys.argv and opens the fake genome, so the functions under test
are extracted and exec'd against a synthetic namespace rather than imported.

Why this file exists, in the spirit of test_indel_pam_creation's note: the guard is
three lines and the interesting failure is not "does it return 0" but "does a row
whose reference alignment is well formed still get exactly the score it got
before". test_wellformed_* pin that, so a future change that widens the guard's
predicate and silently zeroes real reference CFDs is caught in milliseconds
instead of by a diff of a genome-wide run.
"""

import contextlib
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import analysis_notes
import merge_pieces
from bulge_cfd_penalty import cfd_score_bulge_aware
from cfd_score import cfd as cfd_clean, cfd_tables
from test_indel_overlaid_rows import _STUB_INTERVALTREE, _SKIP_WC, _gnu_wc

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "analisi_indels_NNN.py")
# the analyzer's count of the rows the guard scores 0, for the chromosome's notes
COUNTER = "short_ref_rows"

# A 20 nt spacer + NGG, with a DNA bulge ('-' in the crRNA) -- the geometry that
# turned a short reference alignment into an IndexError rather than a wrong score.
CRRNA = "GCTGAGCT-AGGTCAACTGGNGG"
ALT = "GCTGAGCTCAGGTCAACTGGCGG"
REF_OK = "GCTGAGCTCAGGTCAACTGGTGG"
# The same spacer without a bulge, where a short reference alignment raised nothing: the
# scorer scored the bases it had, reading the last two as the PAM.
CRRNA_NO_BULGE = "GCTGAGCTCAGGTCAACTGGNGG"

SENTINEL_REF = 33   # "fill this slot with the reference alignment's CFD"
SENTINEL_ALT = 55   # "fill this slot with the row's own CFD"


def _block():
    """Extract _cfd_from_alignment and preprocess_CFD_score."""
    with open(SRC) as fh:
        src = fh.read()
    a = src.index("def _cfd_from_alignment(")
    end = src.index("\ndef preprocess_CRISPR_BULGE_score(", a)
    return src[a:end]


def _ns(do_scores=True):
    mm_scores, pam_scores = cfd_tables()
    ns = {
        "_cfd_bulge_aware": cfd_score_bulge_aware,
        "_cfd_clean": cfd_clean,
        "mm_scores": mm_scores,
        "pam_scores": pam_scores,
        "do_scores": do_scores,
        COUNTER: 0,  # the module sets it before its loop scores a row
    }
    exec(compile(_block(), "<ref_alignment_guard>", "exec"), ns)
    return ns


def _row(ref_aligned, sentinel=SENTINEL_REF, bulge_type="RNA,DNA", crrna=CRRNA):
    """An indel cluster row as preprocess_CFD_score reads it.

    After the function appends the row's own CFD, the negative indices it uses are
    -1 the appended CFD, -3 the sentinel slot, -4 the reference alignment. So the
    row must end [..., ref_aligned, sentinel, filler].
    """
    return [bulge_type, crrna, ALT, "chr7", "159346029", "159346029", "+",
            "4", "1", "5", "n", "n", "SAMPLE1", "n", "GUIDE", "rs1", "0.01",
            ref_aligned, sentinel, "n"]


class RefAlignmentGuard(unittest.TestCase):
    def test_empty_reference_alignment_scores_zero(self):
        """An out-of-range slice: '' (a row that starts past the contig's end, which up
        to 2.8.1 a row inside a large insertion could). Used to raise IndexError."""
        ns = _ns()
        out = ns["preprocess_CFD_score"](_row(""))
        self.assertEqual(out[-3], "0.000")

    def test_dash_reference_alignment_scores_zero(self):
        """'-': an empty slice with a target's gap put back. Still short, still raised."""
        ns = _ns()
        out = ns["preprocess_CFD_score"](_row("-"))
        self.assertEqual(out[-3], "0.000")

    def test_short_reference_alignment_scores_zero(self):
        """Short but non-empty (2.10.0's case: the slice from an insertion's anchor
        reaches the end of a contig that ends in sequence), whatever the cause."""
        ns = _ns()
        out = ns["preprocess_CFD_score"](_row(REF_OK[:8]))
        self.assertEqual(out[-3], "0.000")

    def test_bulge_free_short_reference_alignment_scores_zero(self):
        """Without a DNA bulge a short alignment raised nothing: the scorer scored the
        bases it had, the last two read as the PAM. Here the reference ends 12 bases into
        a bulge-free target, and those 12 match the guide, so the score was 1.000, a
        perfect reference CFD for a target the reference does not hold. Now 0, counted,
        and the row keeps its own CFD."""
        ns = _ns()
        short = REF_OK[:12]
        self.assertNotEqual("{:.3f}".format(ns["_cfd_from_alignment"](CRRNA_NO_BULGE, short)),
                            "0.000", "the scorer must give the fixture a score, or this test "
                            "passes without the guard")
        out = ns["preprocess_CFD_score"](_row(short, bulge_type="X", crrna=CRRNA_NO_BULGE))
        self.assertEqual(out[-3], "0.000")
        self.assertEqual(out[-1], "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA_NO_BULGE,
                                                                           ALT)))
        self.assertEqual(ns[COUNTER], 1)

    def test_guard_does_not_disturb_the_rows_own_cfd(self):
        """Only the reference comparison is unavailable; the row keeps its ALT CFD."""
        ns = _ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, ALT))
        for ref in ("", "-", REF_OK[:8]):
            out = ns["preprocess_CFD_score"](_row(ref))
            self.assertEqual(out[-1], expected, msg=f"ref_aligned={ref!r}")

    def test_wellformed_reference_alignment_is_untouched(self):
        """THE non-regression pin: a valid row gets exactly the score it got before."""
        ns = _ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, REF_OK))
        out = ns["preprocess_CFD_score"](_row(REF_OK))
        self.assertEqual(out[-3], expected)
        self.assertNotEqual(out[-3], "0.000",
                            "fixture must have a non-zero reference CFD to be meaningful")

    def test_wellformed_equal_length_is_not_zeroed(self):
        """The predicate is len(ref) < len(frame); an equal-length ref must pass through."""
        ns = _ns()
        self.assertEqual(len(REF_OK), len(CRRNA))
        out = ns["preprocess_CFD_score"](_row(REF_OK))
        self.assertNotEqual(out[-3], "0.000")

    def test_alt_sentinel_path_unchanged(self):
        """Sentinel 55 fills the slot with the row's own CFD and never consults the ref."""
        ns = _ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, ALT))
        out = ns["preprocess_CFD_score"](_row("", sentinel=SENTINEL_ALT))
        self.assertEqual(out[-3], expected)

    def test_no_scores_path_unchanged(self):
        """With scoring off, every slot is -1 and the guard is never reached."""
        ns = _ns(do_scores=False)
        out = ns["preprocess_CFD_score"](_row(""))
        self.assertEqual(out[-1], "-1.000")
        self.assertEqual(out[-3], "-1.000")


class TheRowsTheGuardScoresZeroAreCounted(unittest.TestCase):
    def test_each_counts_once_and_no_other_row_counts(self):
        """Each row the guard scores 0 adds one to the count the chromosome's NOTE is written
        from; a well-formed row, a reference row (sentinel 55) and a row with scoring off
        add none."""
        ns = _ns()
        for ref in ("", "-", REF_OK[:8]):
            self.assertEqual(ns["preprocess_CFD_score"](_row(ref))[-3], "0.000")
        self.assertEqual(ns[COUNTER], 3)
        ns["preprocess_CFD_score"](_row(REF_OK))
        ns["preprocess_CFD_score"](_row("", sentinel=SENTINEL_ALT))
        self.assertEqual(ns[COUNTER], 3)
        off = _ns(do_scores=False)
        off["preprocess_CFD_score"](_row(""))
        self.assertEqual(off[COUNTER], 0)


class TheChromosomesNote(unittest.TestCase):
    """One NOTE per chromosome, on stdout, as the rows of another contig get their line:
    from one process's counts, and from the sum of a post-analysis's pieces' counts."""

    START = "NOTE: 3 indel target row(s) of chrE have a reference alignment shorter than "

    @staticmethod
    def notes(lines):
        return [line for line in lines
                if line.startswith("NOTE: ") and " have a reference alignment shorter " in line]

    def test_one_note_on_stdout(self):
        counts = analysis_notes.indel_counts()
        self.assertEqual(self.notes(analysis_notes.indel_lines("chrE", counts)), [])
        counts[COUNTER] = 3
        lines = analysis_notes.indel_lines("chrE", counts)
        notes = self.notes(lines)
        self.assertEqual(len(notes), 1, lines)
        self.assertTrue(notes[0].startswith(self.START), notes[0])
        self.assertTrue(notes[0].endswith("so their CFD_ref is 0."), notes[0])
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            analysis_notes.write(lines)
        self.assertIn(notes[0] + "\n", out.getvalue())
        self.assertEqual(err.getvalue(), "")

    def test_the_pieces_counts_make_one_note(self):
        """A piece writes no note: it leaves its count in its <piece>.piece.json, and the
        merge writes the chromosome's one NOTE from the pieces' sum."""
        tmp = tempfile.mkdtemp(prefix="ref_guard_notes_")
        self.addCleanup(shutil.rmtree, tmp, True)
        prefixes = []
        for i, n in enumerate((2, 0, 1)):
            counts = analysis_notes.indel_counts()
            counts[COUNTER] = n
            prefixes.append(os.path.join(tmp, "piece%d" % i))
            with open(prefixes[-1] + ".piece.json", "w") as fh:
                json.dump({"do_scores": True, "lines": 2, "chrom": "chrE",
                           "notes": {"counts": counts, "precheck_on": True}}, fh)
        notes = self.notes(merge_pieces.indel_notes(prefixes))
        self.assertEqual(len(notes), 1, notes)
        self.assertTrue(notes[0].startswith(self.START), notes[0])


# The world of TheShippedAnalyzerNearAContigEnd: the guide's perfect targets with a DNA
# bulge, none and an RNA bulge, (bulge type, crRNA, DNA, bulge size), planted inside two
# insertions of a contig that ends in sequence; the fake-indel windows cut as CRISPRitz
# 2.10.0 cuts them, FLANK reference bases each side of the REF, the right flank of the
# insertion near the end cut at the contig's end (its enricher's substr).
GUIDE = "CTAACAGTTGCTTTTATCAC"
E_CHROM = "chrE"
E_VCF = "synth"
FLANK = 30
PLANTED = (
    ("DNA", GUIDE[:10] + "-" + GUIDE[10:] + "NNN", GUIDE[:10] + "A" + GUIDE[10:] + "TGG", 1),
    ("X", GUIDE + "NNN", GUIDE + "TGG", 0),
    ("RNA", GUIDE + "NNN", GUIDE[:7] + "-" + GUIDE[8:] + "AGG", 1),
)
# a site of the guide with 2 mismatches, which the reference holds from the far insertion's
# anchor, so that the reference CFD of the rows starting there is one worth comparing
REF_SITE = GUIDE[:5] + "T" + GUIDE[6:15] + "G" + GUIDE[16:] + "AGG"


@unittest.skipUnless(_gnu_wc(), _SKIP_WC)
class TheShippedAnalyzerNearAContigEnd(unittest.TestCase):
    """The real analyzer on a 700-base contig that ends in sequence, with two insertions,
    anchored 400 and 3 bases before its end, each holding the three PLANTED targets. 2.10.0
    starts the rows of each at its anchor (indel_snp_cis.reference_start); from the anchor
    3 bases before the end the reference holds 4 bases, so those three rows' reference
    alignments are short: CFD_ref 0, and one NOTE on stdout counts the 3, also when the
    post-analysis runs in pieces. The other insertion's three rows keep their full reference
    alignment and the score of it. Before the guard the DNA-bulge row raised IndexError
    (the analyzer exited 1, its traceback on stderr) and the bulge-free row's CFD_ref was a
    score of the 4 bases."""

    N = 700
    ANCHORS_FROM_END = (400, 3)

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ref_guard_near_end_")
        rng = random.Random(113)
        ref = [rng.choice("ACGT") for _ in range(cls.N)]
        far = cls.N - cls.ANCHORS_FROM_END[0] - 1
        ref[far:far + len(REF_SITE)] = REF_SITE
        ref = "".join(ref)
        cls.ref = ref
        cls.inserts = []  # (pos1, anchor, near the end)
        rows, log, fake = [], [], ""
        for k, from_end in enumerate(cls.ANCHORS_FROM_END):
            pos1 = cls.N - from_end  # the insertion's VCF POS: its anchor, 1-based
            ins, offsets = "", []
            for _bt, _crrna, dna, _size in PLANTED:
                ins += "".join(rng.choice("ACGT") for _ in range(6))
                offsets.append(len(ins))
                ins += dna.replace("-", "")
            ins += "".join(rng.choice("ACGT") for _ in range(6))
            alt = ref[pos1 - 1] + ins
            w0 = pos1 - 1 - FLANK
            win = ref[w0:pos1 - 1] + alt + ref[pos1:pos1 + FLANK]
            if fake:
                fake += "N"
            f0 = len(fake)
            fake += win
            # the CHR field's end is not cut at the contig's end (indel_window.window_end)
            log.append("%s_%d-%d_%d\tS1,S2\trs%d\t0.05\t%s_%d_%s_%s\t%d,%d\t%s\n" % (
                E_CHROM, w0, pos1 + FLANK, k + 1, k + 1, E_CHROM, pos1, ref[pos1 - 1], alt,
                f0, f0 + len(win), ref[w0:w0 + len(win)]))
            cls.inserts.append((pos1, pos1 - 1, from_end < 23))
            for (bt, crrna, dna, size), off in zip(PLANTED, offsets):
                # the target's first base on the fake contig: the left flank, the anchor,
                # then its offset in the inserted bases; on the + strand with a 3' PAM
                # CRISPRitz's Cluster_Position is it less the RNA bulges, plus the DNA ones
                pos = f0 + FLANK + 1 + off
                cpos = pos - dna.count("-") + crrna.count("-")
                rows.append([bt, crrna, dna, "fake" + E_CHROM, str(pos), str(cpos), "+",
                             "0", str(size), str(size)])
        with open(os.path.join(cls.tmp, E_CHROM + ".fa"), "w") as fh:
            fh.write(">%s\n%s\n" % (E_CHROM, ref))
        cls.logd = os.path.join(cls.tmp, "Dictionaries", "log_indels_" + E_VCF)
        os.makedirs(cls.logd)
        with open(os.path.join(cls.logd, "log%s.txt" % E_CHROM), "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.writelines(log)
        with open(os.path.join(cls.tmp, "pam.txt"), "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRG 3\n")
        with open(os.path.join(cls.tmp, "guides.txt"), "w") as fh:
            fh.write(GUIDE + "NNN\n")
        cls.cluster = os.path.join(cls.tmp, "targets.cluster.txt")
        with open(cls.cluster, "w") as fh:
            fh.write("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
                     "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
                     "Annotation Type\tReal Guide\tID\tAF\tSNP\n")
            for r in rows:
                fh.write("\t".join(r + ["n", "n", "n", "n", GUIDE + "NNN", "n", "n", "n"])
                         + "\n")
        cls.n_rows = len(rows)
        cls.env = {k: v for k, v in os.environ.items()
                   if k not in ("MAMBA_EXE", "CONDA_EXE", "CRISPRME_CONDA_EXE")}
        cls.env["PATH"] = os.path.dirname(sys.executable) + ":/usr/bin:/bin"
        try:
            import intervaltree  # noqa: F401
        except ImportError:
            stub = os.path.join(cls.tmp, "stub")
            os.makedirs(stub)
            with open(os.path.join(stub, "intervaltree.py"), "w") as fh:
                fh.write(_STUB_INTERVALTREE)
            cls.env["PYTHONPATH"] = stub + os.pathsep + cls.env.get("PYTHONPATH", "")
        os.makedirs(os.path.join(cls.tmp, "out"))
        cls.out = os.path.join(cls.tmp, "out", "fake" + E_CHROM)
        cls.proc = cls.analyze(cls.out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def analyze(cls, out, env=None):
        here = os.path.dirname(os.path.abspath(__file__))
        argv = [sys.executable, os.path.join(here, "analisi_indels_NNN.py"), "", cls.cluster,
                out, cls.logd, os.path.join(cls.tmp, "pam.txt"), "4",
                os.path.join(cls.tmp, E_CHROM + ".fa"), os.path.join(cls.tmp, "guides.txt"),
                "1", "1"]
        return subprocess.run(argv, cwd=here, env=dict(cls.env, **(env or {})),
                              capture_output=True, text=True, timeout=600)

    def test_the_short_rows_score_zero_and_the_others_as_before(self):
        self.assertEqual(self.proc.returncode, 0,
                         self.proc.stdout[-2000:] + self.proc.stderr[-2000:])
        self.assertEqual(self.proc.stderr, "")
        mm_scores, pam_scores = cfd_tables()
        rows = {}
        with open(self.out + ".bestCFD_INDEL.txt") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                rows.setdefault(int(c[17].split("_")[1]), []).append(c)
        self.assertEqual(sorted(rows), sorted(p for p, _a, _near in self.inserts))
        for pos1, anchor, near in self.inserts:
            got = sorted(rows[pos1], key=lambda c: c[0])
            self.assertEqual([c[0] for c in got], sorted(p[0] for p in PLANTED))
            for c in got:
                where = "%s row of the insertion at %d" % (c[0], pos1)
                self.assertEqual(int(c[4]), anchor, where + ": starts at the anchor")
                if near:
                    # the reference from the anchor: its last 4 bases
                    self.assertEqual(c[18].replace("-", "").upper(), self.ref[anchor:], where)
                    self.assertLess(len(c[18]), len(c[1]), where)
                    self.assertEqual(c[19], "0.000", where + ": CFD_ref")
                else:
                    self.assertEqual(len(c[18]), len(c[1]), where)
                    want = cfd_score_bulge_aware(c[1], c[18], mm_scores, pam_scores)
                    self.assertEqual(c[19], "{:.3f}".format(want), where + ": CFD_ref")
            if not near:
                self.assertNotEqual([c[19] for c in got if c[0] == "X"], ["0.000"],
                                    "the far bulge-free row must have a reference CFD")

    def test_one_note_counts_them_on_stdout(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr[-2000:])
        notes = TheChromosomesNote.notes(self.proc.stdout.splitlines())
        self.assertEqual(len(notes), 1, self.proc.stdout[-2000:])
        self.assertTrue(notes[0].startswith(TheChromosomesNote.START), notes[0])
        self.assertEqual(self.proc.stderr, "")

    def test_in_pieces_the_merge_writes_the_one_note(self):
        """Two pieces (lines 1-4 and 5-6: the far insertion's rows and the DNA-bulge row
        near the end, then the other two near the end) count 1 and 2, write no note, and
        the merge writes the one NOTE of 3."""
        prefixes, counts = [], []
        for i, rng in enumerate(("1:4", "5:%d" % self.n_rows)):
            prefixes.append(os.path.join(self.tmp, "out", "piece%d" % i))
            done = self.analyze(prefixes[-1], {"CRISPRME_POSTPROC_PIECE": "1",
                                               "CRISPRME_POSTPROC_RANGE": rng})
            self.assertEqual(done.returncode, 0, done.stdout[-2000:] + done.stderr[-2000:])
            self.assertEqual(done.stderr, "")
            self.assertEqual(TheChromosomesNote.notes(done.stdout.splitlines()), [])
            with open(prefixes[-1] + ".piece.json") as fh:
                counts.append(json.load(fh)["notes"]["counts"][COUNTER])
        self.assertEqual(counts, [1, 2])
        notes = TheChromosomesNote.notes(merge_pieces.indel_notes(prefixes))
        self.assertEqual(len(notes), 1, notes)
        self.assertTrue(notes[0].startswith(TheChromosomesNote.START), notes[0])


if __name__ == "__main__":
    unittest.main()

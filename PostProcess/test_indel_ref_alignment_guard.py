"""The reference-alignment CFD of a row when there is no reference.

Both post-analysis paths score an alt row's REFERENCE alignment, and both can be
handed one shorter than the alignment frame.

alignRefFromVar builds a row's reference alignment by slicing the plain reference
chromosome at the row's coordinates. An indel row can sit past the end of its
chromosome -- the per-chromosome fake indel genome is longer than the real
chromosome, so a target found inside a large insertion has no reference position
at all -- and Python returns '' for an out-of-range slice instead of raising. The
empty alignment then reached cfd_score_bulge_aware, whose documented precondition
is that the two alignments are the same length; with a DNA bulge in the guide it
indexes the protospacer by the crRNA's positions and died with
"IndexError: string index out of range", taking the whole indel post-analysis --
and, because post_analisi_indel.sh exits 0 regardless, a whole chromosome's indel
targets -- with it.

The SNP path reaches the same call by a different route. refSeq is
genomeStr[pos : pos + len(replaceTarget)] and, near a chromosome boundary or from
an out-of-range position on the variant-enriched combined index, it comes back
truncated -- the hazard new_simple_analysis's own boundary guards document. Those
guards `continue` past the offending IUPAC POSITION so per-position indexing cannot
raise; they neither drop the row nor restore the missing bases, and _ref_with_bulges
carries the truncation into target[-4]. So the short alignment reaches the scorer
there too, and both paths are guarded and tested here together.

Following test_indel_pam_creation.py and test_indel_snp_overlay.py, each module's
top level reads sys.argv and opens its genome, so the functions under test are
extracted and exec'd against a synthetic namespace rather than imported.

Why this file exists, in the spirit of test_indel_pam_creation's note: the guard is
three lines and the interesting failure is not "does it return 0" but "does a row
whose reference alignment is well formed still get exactly the score it got
before". test_wellformed_* pin that, so a future change that widens the guard's
predicate and silently zeroes real reference CFDs is caught in milliseconds
instead of by a diff of a genome-wide run.
"""

import os
import unittest

from bulge_cfd_penalty import cfd_score_bulge_aware
from cfd_score import cfd as cfd_clean, cfd_tables

HERE = os.path.dirname(__file__)
INDEL_SRC = os.path.join(HERE, "analisi_indels_NNN.py")
SNP_SRC = os.path.join(HERE, "new_simple_analysis.py")

# A 20 nt spacer + NGG, with a DNA bulge ('-' in the crRNA) -- the geometry that
# turned a short reference alignment into an IndexError rather than a wrong score.
CRRNA = "GCTGAGCT-AGGTCAACTGGNGG"
ALT = "GCTGAGCTCAGGTCAACTGGCGG"
REF_OK = "GCTGAGCTCAGGTCAACTGGTGG"

SENTINEL_REF = 33   # "fill this slot with the reference alignment's CFD"
SENTINEL_ALT = 55   # "fill this slot with the row's own CFD"


def _block(src_path):
    """Extract _cfd_from_alignment and preprocess_CFD_score. Both modules define the
    two adjacently and follow them with preprocess_CRISPR_BULGE_score."""
    with open(src_path) as fh:
        src = fh.read()
    a = src.index("def _cfd_from_alignment(")
    end = src.index("\ndef preprocess_CRISPR_BULGE_score(", a)
    return src[a:end]


def _ns(src_path, do_scores=True):
    mm_scores, pam_scores = cfd_tables()
    ns = {
        "_cfd_bulge_aware": cfd_score_bulge_aware,
        "_cfd_clean": cfd_clean,
        "mm_scores": mm_scores,
        "pam_scores": pam_scores,
        "do_scores": do_scores,
    }
    exec(compile(_block(src_path), "<ref_alignment_guard>", "exec"), ns)
    return ns


def _row(ref_aligned, sentinel=SENTINEL_REF):
    """An indel cluster row as preprocess_CFD_score reads it.

    After the function appends the row's own CFD, the negative indices it uses are
    -1 the appended CFD, -3 the sentinel slot, -4 the reference alignment. So the
    row must end [..., ref_aligned, sentinel, filler].
    """
    return ["RNA,DNA", CRRNA, ALT, "chr7", "159346029", "159346029", "+",
            "4", "1", "5", "n", "n", "SAMPLE1", "n", "GUIDE", "rs1", "0.01",
            ref_aligned, sentinel, "n"]


class _RefAlignmentGuard:
    """Shared cases. SRC is set by the two concrete classes below."""

    SRC = None

    def _ns(self, do_scores=True):
        return _ns(self.SRC, do_scores=do_scores)

    def test_empty_reference_alignment_scores_zero(self):
        """The 2.8.0 out-of-range slice: '' . Used to raise IndexError."""
        ns = self._ns()
        out = ns["preprocess_CFD_score"](_row(""))
        self.assertEqual(out[-3], "0.000")

    def test_dash_reference_alignment_scores_zero(self):
        """2.8.1 rewrote alignRefFromVar to return '-' here. Still short, still raised."""
        ns = self._ns()
        out = ns["preprocess_CFD_score"](_row("-"))
        self.assertEqual(out[-3], "0.000")

    def test_short_reference_alignment_scores_zero(self):
        """Short but non-empty, whatever the cause."""
        ns = self._ns()
        out = ns["preprocess_CFD_score"](_row(REF_OK[:8]))
        self.assertEqual(out[-3], "0.000")

    def test_guard_does_not_disturb_the_rows_own_cfd(self):
        """Only the reference comparison is unavailable; the row keeps its ALT CFD."""
        ns = self._ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, ALT))
        for ref in ("", "-", REF_OK[:8]):
            out = ns["preprocess_CFD_score"](_row(ref))
            self.assertEqual(out[-1], expected, msg=f"ref_aligned={ref!r}")

    def test_wellformed_reference_alignment_is_untouched(self):
        """THE non-regression pin: a valid row gets exactly the score it got before."""
        ns = self._ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, REF_OK))
        out = ns["preprocess_CFD_score"](_row(REF_OK))
        self.assertEqual(out[-3], expected)
        self.assertNotEqual(out[-3], "0.000",
                            "fixture must have a non-zero reference CFD to be meaningful")

    def test_wellformed_equal_length_is_not_zeroed(self):
        """The predicate is len(ref) < len(frame); an equal-length ref must pass through."""
        ns = self._ns()
        self.assertEqual(len(REF_OK), len(CRRNA))
        out = ns["preprocess_CFD_score"](_row(REF_OK))
        self.assertNotEqual(out[-3], "0.000")

    def test_alt_sentinel_path_unchanged(self):
        """Sentinel 55 fills the slot with the row's own CFD and never consults the ref."""
        ns = self._ns()
        expected = "{:.3f}".format(ns["_cfd_from_alignment"](CRRNA, ALT))
        out = ns["preprocess_CFD_score"](_row("", sentinel=SENTINEL_ALT))
        self.assertEqual(out[-3], expected)

    def test_no_scores_path_unchanged(self):
        """With scoring off, every slot is -1 and the guard is never reached."""
        ns = self._ns(do_scores=False)
        out = ns["preprocess_CFD_score"](_row(""))
        self.assertEqual(out[-1], "-1.000")
        self.assertEqual(out[-3], "-1.000")


class IndelPathRefAlignmentGuard(_RefAlignmentGuard, unittest.TestCase):
    """analisi_indels_NNN.py -- where this was found, on 193 real chr7 rows."""

    SRC = INDEL_SRC


class SnpPathRefAlignmentGuard(_RefAlignmentGuard, unittest.TestCase):
    """new_simple_analysis.py -- the identical call, reachable via a truncated refSeq."""

    SRC = SNP_SRC


if __name__ == "__main__":
    unittest.main()

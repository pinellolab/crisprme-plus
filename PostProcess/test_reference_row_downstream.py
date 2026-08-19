"""DOWNSTREAM regression for the dict-less REFERENCE off-target row.

CONTEXT
-------
The dict-less observed path (``new_simple_analysis._iupac_decomposition_observed``)
emits, for every in-budget IUPAC candidate, a REFERENCE off-target row via
``_finalize_reference_entry`` so a candidate with NO productive observed haplotype
still keeps its locus (matches the separate reference-genome CRISPRitz search). A real
genome-wide e2e then CRASHED in ``resultIntegrator.py`` with::

    IndexError: list index out of range
    variantList = str(target[18]).split(",")  ->  elem.split("_")  ->  len(split[2])

ROOT CAUSE (reproduced here end-to-end)
---------------------------------------
The whole pipeline classifies a row origin=ref by the LITERAL marker ``"n"`` in the
Samples / rsID / AF / SNP columns:

  * ``remove_contiguous_samples.py`` / ``merge_contiguous_targets.py`` bin a row
    origin=ref ONLY when the SNP column (col18 post-adjust) == ``"n"``;
  * ``remove_n_and_dots.py`` then rewrites that literal ``"n"`` into the literal
    string ``"NA"`` (``chunk.replace({"n": "NA"})``), which ``resultIntegrator``
    reads as a REAL string -> ``target[18] != "NA"`` is FALSE (guard skips the
    variant parse) and ``target[13] == "NA"`` is TRUE (origin=ref).

The buggy ``_finalize_reference_entry`` wrote ``"NA"`` (not ``"n"``) into those
columns. That broke the invariant twice:

  1. SNP == "NA" (not "n") -> the merge scripts did NOT bin the row origin=ref; a
     ref-only candidate became ``var_only`` and was emitted as a spurious "best
     variant".
  2. The literal ``"NA"`` written to disk is read back by pandas
     (``remove_n_and_dots``) as NaN and re-serialized as an EMPTY string, so
     ``resultIntegrator`` saw SNP == "" -> ``"" != "NA"`` is TRUE -> it entered the
     variant parse -> ``"".split("_")`` == [""] -> ``len(split[2])`` -> IndexError.

WHAT THIS TEST DOES
-------------------
Drives the REAL ``_finalize_reference_entry`` + ``_finalize_observed_entry`` (via the
same AST harness ``test_phased_haplotype`` uses -- no genome/dict/arg files), serializes
the finalized rows into an on-disk POST-ADJUST ``bestCFD.txt`` (Reference at col3, SNP
at col18, exactly as ``adjust_cols.py`` reorders), then runs the REAL
``remove_contiguous_samples.py`` (pure stdlib) and -- when pandas is available --
``remove_n_and_dots.py``, before feeding each surviving row through the EXACT
``resultIntegrator`` col-18 variant-parse and the origin rule.

ASSERTS (fails on f48bb9e, passes after the "NA" -> "n" fix):
  * the reference locus survives the merge as origin=ref (SNP col18 == "n", NOT routed
    to ``.discarded_samples`` as a variant);
  * after ``remove_n_and_dots`` + the resultIntegrator parse: NO IndexError on any row;
  * the reference row classifies origin=ref (Samples reads literal "NA"), the variant
    row origin=alt;
  * negative control: an isolated SNP marker in {"", "n"} DOES raise IndexError in the
    resultIntegrator parse, proving the guard depends on the literal marker.

Run with:
    cd PostProcess && python3 -m unittest test_reference_row_downstream -v
"""

import io
import os
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN


# ----------------------------------------------------------------------------- #
# Column layouts (see the diagnosis).
#
# PRE-ADJUST finalized row (what _finalize_* build) = 18 split cols + a 3-elem tail
#   [ Reference, sentinel(33/55), tmp_pos ].
# The write loop pops tmp_pos (index -2) and appends the CFD; preprocess_CFD_score
# rewrote the sentinel (now at index 19) into CFD_ref. So the ON-DISK PRE-ADJUST
# bestCFD.txt row (22 cols) is:
#   0-17 split, 18 Reference, 19 CFD_ref, 20 CFD, 21 #Seq_in_cluster
PRE_ADJUST_HEADER = (
    "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\t"
    "Mismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\t"
    "Real_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster"
)
# POST-ADJUST (adjust_cols.py: cols[:3] + ["Reference"] + cols[3:] + ["CFD","CFD_ref"]).
# Reference moves to col3; SNP lands at col18; Samples at col13; CFD at col20. This is
# the layout every merge / integrator consumer sees.
SNP_COL = 18
SAMPLES_COL = 13
CFD_COL = 20


class _GtSentinel(object):
    """Truthy stand-in for a Tier-1 GenotypeReader (observed branch only checks
    ``mygt is not None``)."""

    def __len__(self):
        return 1


def _finalized_rows(gt_by_col, n_snps=None, split4="0", allowed_mms=6):
    """Drive the REAL observed path and return its finalized (pre-adjust) rows.

    ``gt_by_col[i]`` is the genotype token for ambiguity column i (e.g. "S:0|0" for
    ref/ref -> no productive variant, "S:1|0" -> S carries the alt)."""
    if n_snps is None:
        n_snps = len(gt_by_col)
    ref = ["A"] * TARGET_LEN
    for i in range(n_snps):
        ref[i] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"
    genome = "".join(ref)
    dna = list(genome)
    for i in range(n_snps):
        dna[i] = "W"
    dna = "".join(dna)
    split = ["X", GUIDE, dna, "chrT", split4, split4, "+",
             str(n_snps), "0", str(n_snps), "NGG", "y", "NA",
             "NA", "NA", "NA", "NA", "NA"]
    mydict = {}
    for i in range(n_snps):
        mydict["chrT," + str(i + 1)] = gt_by_col[i] + ";C,A;rs%d;0.10" % i
    overrides = dict(
        genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
        mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
        hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
        pam_begin=-3, pam_end=None, allowed_mms=allowed_mms,
    )
    ns = _load_pure_functions(overrides)
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    cluster = []
    ns["iupac_decomposition"](split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster


def _to_disk_pre_adjust(final_row, chrom, position, cfd="0.500"):
    """Apply the deterministic write-loop transform (pop tmp_pos at -2, append CFD;
    the sentinel at col19 becomes the CFD_ref score) to a finalized 21-col row,
    overriding Chromosome/Position/Cluster_Position so we can place rows into distinct
    clusters. Returns the 22-col on-disk PRE-ADJUST row (strings)."""
    row = [str(x) for x in final_row]
    # place the row at a controllable locus (Chromosome col3, Position col4,
    # Cluster_Position col5) so we can build separate clusters deterministically.
    row[3] = str(chrom)
    row[4] = str(position)
    row[5] = str(position)
    reference = row[18]           # "n" for a ref row, refSeq for a variant row
    sentinel = final_row[19]      # 55 (ref) or 33 (variant), int
    # tmp_pos is row[20]; the write loop pops index -2 AFTER the CFD append. Emulate:
    # pop tmp_pos, set CFD_ref (was the sentinel) to a numeric score, append CFD.
    cfd_ref = cfd if sentinel == 55 else "0.400"
    out = row[:18] + [reference, cfd_ref, cfd, "0"]
    return out


def _adjust_cols(row22):
    """Pure-python mirror of adjust_cols.py's reorder:
    cols[:3] + ["Reference"] + cols[3:] + ["CFD","CFD_ref"], dropping the original
    CFD/CFD_ref/Reference. Input is the 22-col pre-adjust row; output the 22-col
    post-adjust row (Reference at col3, SNP at col18, CFD at col20)."""
    # pre-adjust indices: 18 Reference, 19 CFD_ref, 20 CFD, 21 #Seq_in_cluster
    reference, cfd_ref, cfd = row22[18], row22[19], row22[20]
    cols = row22[:18] + [row22[21]]  # drop Reference/CFD_ref/CFD, keep #Seq_in_cluster
    return cols[:3] + [reference] + cols[3:] + [cfd, cfd_ref]


def _post_adjust_header():
    cols = PRE_ADJUST_HEADER.split("\t")
    reference, cfd_ref, cfd = cols[18], cols[19], cols[20]
    rest = cols[:18] + [cols[21]]
    return rest[:3] + [reference] + rest[3:] + [cfd, cfd_ref]


def _integrator_variant_parse(snp_cell):
    """The EXACT resultIntegrator col-18 variant-parse (lines ~409-414). Raises
    IndexError on the bad markers ("" / "n") and is a no-op on "NA" / real variants."""
    t18 = str(snp_cell)
    if t18 != "NA":
        for elem in t18.strip().split(","):
            split = str(elem).strip().split("_")
            _ = len(split[2])   # IndexError here for "" / "n"
            _ = len(split[3])


def _integrator_origin(samples_cell):
    """The EXACT resultIntegrator origin rule (lines 550-551): ref iff Samples == NA."""
    return "ref" if str(samples_cell) == "NA" else "alt"


try:
    import pandas as _pd  # noqa: F401
    _HAVE_PANDAS = True
except Exception:
    _HAVE_PANDAS = False


class ReferenceRowDownstream(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="refrow_down_")

    def _write_post_adjust(self, rows):
        path = os.path.join(self.d, "out.bestCFD.txt")
        with open(path, "w") as fh:
            fh.write("\t".join(_post_adjust_header()) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        return path

    def _build_mixed_and_ref_only(self):
        """Build TWO clusters (both real, from the patched module):
          cluster A (variant + its reference) at locus 100 -- a productive candidate;
          cluster B (reference only) at locus 100000 -- a no-productive-variant
        candidate (the exact e2e-crash case). Returns the post-adjust rows +
        bookkeeping."""
        # Cluster A: sample S carries the single alt -> one variant row + one ref row.
        rows_A = _finalized_rows(["S:1|0"], n_snps=1)
        # Cluster B: all ref/ref -> observed enumeration empty -> ONLY the ref row.
        rows_B = _finalized_rows(["S:0|0", "S:0|0", "S:0|0"], n_snps=3)
        self.assertTrue(rows_A, "cluster A must emit rows")
        self.assertEqual(len(rows_B), 1, "cluster B must be reference-only")

        post = []
        # place cluster A at position 100, cluster B (ref-only) at position 100000 so
        # the merge treats them as separate clusters (tau=3).
        for r in rows_A:
            post.append(_adjust_cols(_to_disk_pre_adjust(r, "chrT", 100)))
        for r in rows_B:
            post.append(_adjust_cols(_to_disk_pre_adjust(r, "chrT", 100000)))
        return post

    # --------------------------------------------------------------------- #
    # (1) MERGE leg -- pure stdlib, always runs (no pandas).
    # --------------------------------------------------------------------- #
    def test_merge_bins_reference_row_as_ref_not_variant(self):
        """The real remove_contiguous_samples.py must bin the ref-only candidate's
        reference row origin=ref (into the best file), NOT route it to
        .discarded_samples as a spurious variant. Pre-fix (SNP=="NA") it was wrongly
        binned var_only -> emitted as a best variant with SNP != "n"."""
        post = self._build_mixed_and_ref_only()
        infile = self._write_post_adjust(post)
        outfile = os.path.join(self.d, "merged.bestCFD.txt")
        # merge_close_targets_cfd.sh 1-based indices: chrom=5 pos=7 total=11
        # true_guide=16 snp_info=19 cfd=21; sort_order=score, criteria=mm+bulges.
        cmd = [sys.executable, os.path.join(PP, "remove_contiguous_samples.py"),
               infile, outfile, "3", "5", "7", "11", "16", "19", "21",
               "score", "mm+bulges"]
        r = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, "merge failed: %s" % r.stderr)

        def _read(path):
            with open(path) as fh:
                lines = [l.rstrip("\n") for l in fh]
            return [l.split("\t") for l in lines[1:] if l.strip()]

        best = _read(outfile)
        disc = _read(outfile + ".discarded_samples")
        # The ref-only candidate lives at chrT:100000. It must appear as a reference
        # row (SNP col18 == "n") in the BEST file, never as a variant.
        best_ref_only = [row for row in best if row[5] == "100000"]
        self.assertEqual(len(best_ref_only), 1,
                         "the ref-only candidate must survive the merge as ONE best "
                         "row; got %d (best=%r)"
                         % (len(best_ref_only), [r[5] for r in best]))
        self.assertEqual(best_ref_only[0][SNP_COL], "n",
                         "the ref-only candidate must be binned origin=ref (SNP col18 "
                         "== 'n'); pre-fix it was 'NA' -> misrouted as a variant")
        self.assertEqual(best_ref_only[0][SAMPLES_COL], "n",
                         "reference row Samples must be the legacy literal 'n'")
        # It must NOT have been routed to discarded_samples as a variant.
        disc_ref_only = [row for row in disc if row[5] == "100000"]
        self.assertFalse(disc_ref_only,
                         "the ref-only candidate must NOT be discarded as a variant")
        return outfile

    # --------------------------------------------------------------------- #
    # (2) CRASH GATE -- remove_n_and_dots (pandas) + resultIntegrator parse.
    # --------------------------------------------------------------------- #
    @unittest.skipUnless(_HAVE_PANDAS, "remove_n_and_dots.py requires pandas")
    def test_remove_n_then_integrator_parse_no_indexerror(self):
        """Full crash gate: merge -> remove_n_and_dots -> the EXACT resultIntegrator
        col-18 variant-parse + origin rule on every surviving row. Pre-fix the ref-only
        candidate's SNP became '' after the pandas NaN round-trip -> IndexError; now it
        is the literal 'NA' -> guarded, origin=ref."""
        post = self._build_mixed_and_ref_only()
        infile = self._write_post_adjust(post)
        outfile = os.path.join(self.d, "merged.bestCFD.txt")
        cmd = [sys.executable, os.path.join(PP, "remove_contiguous_samples.py"),
               infile, outfile, "3", "5", "7", "11", "16", "19", "21",
               "score", "mm+bulges"]
        r = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, "merge failed: %s" % r.stderr)

        # run the REAL remove_n_and_dots.py on BOTH the best and the discarded file.
        for f in (outfile, outfile + ".discarded_samples"):
            if os.stat(f).st_size <= 0:
                continue
            rc = subprocess.run(
                [sys.executable, os.path.join(PP, "remove_n_and_dots.py"), f],
                cwd=PP, capture_output=True, text=True)
            self.assertEqual(rc.returncode, 0, "remove_n_and_dots failed: %s" % rc.stderr)

        # feed every surviving data row through the EXACT resultIntegrator parse.
        def _rows(path):
            with open(path) as fh:
                lines = [l.rstrip("\n") for l in fh]
            return [l.split("\t") for l in lines[1:] if l.strip()]

        seen_ref = False
        seen_alt = False
        for path in (outfile, outfile + ".discarded_samples"):
            if not os.path.isfile(path) or os.stat(path).st_size <= 0:
                continue
            for row in _rows(path):
                # THIS is the line that raised IndexError before the fix.
                _integrator_variant_parse(row[SNP_COL])
                origin = _integrator_origin(row[SAMPLES_COL])
                if row[5] == "100000":
                    self.assertEqual(origin, "ref",
                                     "ref-only candidate must classify origin=ref "
                                     "(Samples must read literal 'NA'), got %r "
                                     "(Samples=%r SNP=%r)"
                                     % (origin, row[SAMPLES_COL], row[SNP_COL]))
                    seen_ref = True
                elif row[5] == "100" and origin == "alt":
                    seen_alt = True
        self.assertTrue(seen_ref, "the reference locus must be present as origin=ref")
        self.assertTrue(seen_alt, "the variant locus must be present as origin=alt")

    # --------------------------------------------------------------------- #
    # (3) NEGATIVE CONTROL -- the parse DOES crash on the bad markers.
    # --------------------------------------------------------------------- #
    def test_negative_control_bad_markers_raise(self):
        for bad in ("", "n"):
            with self.assertRaises(IndexError,
                                   msg="SNP %r must raise IndexError in the "
                                       "resultIntegrator parse" % bad):
                _integrator_variant_parse(bad)
        # the CORRECT markers do NOT raise.
        for ok in ("NA", "chrT_101_C_A"):
            _integrator_variant_parse(ok)


if __name__ == "__main__":
    unittest.main()

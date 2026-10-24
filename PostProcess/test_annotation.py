"""Regression test for CRISPRme's post-search annotation step
(``PostProcess/annotation.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_annotation.py' -v

``annotate_offtargets()`` used to split each off-target row with the
default, whitespace-collapsing ``str.split()`` instead of an explicit
``\\t``-split. Any row with a genuinely empty tab-separated column (e.g. a
blank AF value, reachable once CRISPRitz PR #36 lets AF-not-found return
``""`` instead of erroring) silently lost that column, shifting every later
column left by one and breaking the fixed-index readers downstream
(``add_risk_score.py``, ``resultIntegrator.py``) with an ``IndexError``.
"""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PostProcess import annotation as an  # noqa: E402


# 22-column header used by the primary bestCFD/bestmmblg/bestCRISPR_BULGE
# intermediates (submit_job_automated_new_multiple_vcfs.sh).
HEADER = (
    "#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\t"
    "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
    "Annotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref"
)


def _row(af=""):
    """A realistic 22-field off-target row; AF (index 17) blank by default,
    matching the "no AF field anywhere in the VCF" / AF-bounds-overflow case.
    """
    fields = [
        "DNA", "CTAAC-AGTTGCTTTTATCACNNN", "CTccCAAGTTGCTgTgATCACAGG",
        "CTccCAAGTTGCTgggATCACAGG", "chr22", "10730842", "10730843", "+",
        "4", "1", "5", "n", "y", "HG04185,HG04171", "n",
        "CTAACAGTTGCTTTTATCACNNN", ".", af, "chr22_10730857_G_T", "0",
        "0.0", "0.0",
    ]
    assert len(fields) == 22
    return "\t".join(fields)


class AnnotateOfftargetsTest(unittest.TestCase):
    def _run(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            infile = os.path.join(tmp, "in.txt")
            outfile = os.path.join(tmp, "out.txt")
            with open(infile, "w") as fh:
                fh.write(HEADER + "\n")
                for row in rows:
                    fh.write(row + "\n")
            an.annotate_offtargets(infile, None, outfile)  # annotation=None: vuoto.txt path
            with open(outfile) as fh:
                lines = [line.rstrip("\n") for line in fh]
        return lines

    def test_blank_af_field_preserved(self):
        """A blank AF value must survive as its own empty column, not vanish."""
        out = self._run([_row(af="")])
        header, data = out[0], out[1]
        self.assertEqual(len(header.split("\t")), 22)
        fields = data.split("\t")
        self.assertEqual(len(fields), 22, f"expected 22 fields, got {len(fields)}: {fields}")
        self.assertEqual(fields[17], "", "AF column should be blank, not dropped/shifted")
        # everything after AF must still be at its original index
        self.assertEqual(fields[18], "chr22_10730857_G_T")
        self.assertEqual(fields[20], "0.0")
        self.assertEqual(fields[21], "0.0")

    def test_populated_af_field_unaffected(self):
        """Normal (non-blank) rows must be unaffected by the parsing fix."""
        out = self._run([_row(af="0.02")])
        fields = out[1].split("\t")
        self.assertEqual(len(fields), 22)
        self.assertEqual(fields[17], "0.02")
        self.assertEqual(fields[21], "0.0")

    def test_annotation_type_column_still_gets_rewritten(self):
        """annotate_target() writes into fields[14] (Annotation_Type); confirm
        the explicit-split fix didn't disturb that in-place rewrite."""
        out = self._run([_row(af="0.1")])
        fields = out[1].split("\t")
        self.assertEqual(fields[14], "n")  # annotation=None -> always "n"

    def test_mixed_blank_and_populated_rows(self):
        """A blank-AF row and a normal row in the same file must both come out
        with 22 fields (guards against any per-file, not just per-row, state)."""
        out = self._run([_row(af=""), _row(af="0.5")])
        for data in out[1:]:
            self.assertEqual(len(data.split("\t")), 22)


class LoadAnnotationBedRetryTest(unittest.TestCase):
    """load_annotation_bed() must tolerate a *transient* truncated read of a valid
    annotation (a Docker Desktop virtiofs bind mount / NFS read-after-write lag can
    briefly present a complete file as short), retrying with backoff before declaring
    it corrupt -- while still raising an actionable error on a *persistent* failure.
    """

    def _gz(self, tmp):
        path = os.path.join(tmp, "annotation.sorted.bed.gz")
        open(path, "wb").close()  # real file so getmtime()/lock creation work
        return path

    def test_transient_read_then_success(self):
        """Two transient open failures then success -> returns the handle, having
        slept between attempts (no spurious 'corrupt' error)."""
        sentinel = object()
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            tbx = mock.Mock(side_effect=[OSError("Invalid BGZF header at offset 104282"),
                                         OSError("could not open index"), sentinel])
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex") as reindex, \
                 mock.patch.object(an.pysam, "TabixFile", tbx), \
                 mock.patch.object(an, "sleep") as slp:
                result = an.load_annotation_bed(gz)
            self.assertIs(result, sentinel)
            self.assertEqual(tbx.call_count, 3)       # failed twice, succeeded 3rd
            self.assertEqual(slp.call_count, 2)       # backed off before each retry
            self.assertTrue(reindex.called)

    def test_success_first_try_never_sleeps(self):
        """A healthy filesystem path must not pay any retry latency."""
        sentinel = object()
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex"), \
                 mock.patch.object(an.pysam, "TabixFile", return_value=sentinel), \
                 mock.patch.object(an, "sleep") as slp:
                result = an.load_annotation_bed(gz)
            self.assertIs(result, sentinel)
            slp.assert_not_called()

    def _valid_tbi(self, gz):
        """A non-empty .tbi, newer than the .gz, so the stale-check treats it as current
        (content is irrelevant here -- _reindex is mocked)."""
        tbi = gz + ".tbi"
        with open(tbi, "wb") as fh:
            fh.write(b"\x1f\x8b" + b"\0" * 64)
        os.utime(tbi, None)
        return tbi

    def test_corrupt_nonstale_index_is_rebuilt(self):
        """A present, non-empty, not-older .tbi that fails to OPEN is corrupt -> rebuilt on
        the next attempt. Regression guard: the stale-check alone would never rebuild it
        (this restores the base self-heal the retry loop dropped)."""
        sentinel = object()
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            self._valid_tbi(gz)
            tbx = mock.Mock(side_effect=[OSError("Invalid BGZF header"), sentinel])
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex") as reindex, \
                 mock.patch.object(an.pysam, "TabixFile", tbx), \
                 mock.patch.object(an, "sleep"):
                result = an.load_annotation_bed(gz)
            self.assertIs(result, sentinel)
            self.assertEqual(tbx.call_count, 2)   # corrupt open, rebuild, success
            self.assertTrue(reindex.called)       # the corrupt-not-stale index WAS rebuilt

    def test_valid_index_readonly_no_write(self):
        """A VALID index opens with ZERO writes (read-only-install safe): no reindex, and no
        .tbi.lock is created."""
        sentinel = object()
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            self._valid_tbi(gz)
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex") as reindex, \
                 mock.patch.object(an.pysam, "TabixFile", return_value=sentinel), \
                 mock.patch.object(an, "sleep"):
                result = an.load_annotation_bed(gz)
            self.assertIs(result, sentinel)
            reindex.assert_not_called()
            self.assertFalse(os.path.exists(gz + ".tbi.lock"))

    def test_readonly_dir_needing_rebuild_reports_permission_not_corruption(self):
        """Index missing/stale AND the directory not writable -> the error names the
        read-only directory, not a (false) 'the file may be corrupt. Re-download it'."""
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)  # no .tbi -> a rebuild is required
            real_open = open

            def fake_open(path, *a, **k):
                if str(path).endswith(".tbi.lock"):
                    raise PermissionError(13, "Permission denied")
                return real_open(path, *a, **k)

            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex"), \
                 mock.patch("builtins.open", side_effect=fake_open), \
                 mock.patch.object(an, "sleep"), \
                 self.assertRaises(an.SamtoolsError) as ctx:
                an.load_annotation_bed(gz)
            msg = str(ctx.exception).lower()
            self.assertIn("not writable", msg)
            self.assertNotIn("re-download", msg)

    def test_persistent_missing_eof_is_truncated_download(self):
        """A persistently absent BGZF EOF marker -> actionable 'truncated/incomplete
        download' message (not a generic corruption error), only after retries."""
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=False), \
                 mock.patch.object(an, "sleep") as slp, \
                 self.assertRaises(an.SamtoolsError) as ctx:
                an.load_annotation_bed(gz)
            self.assertIn("truncated or corrupt", str(ctx.exception))
            self.assertIn("incomplete download", str(ctx.exception))
            self.assertEqual(slp.call_count, an._ANNOT_LOAD_ATTEMPTS - 1)

    def test_persistent_load_failure_is_corrupt_after_attempts(self):
        """EOF present but every open fails -> 'corrupt even after N attempts'."""
        with tempfile.TemporaryDirectory() as tmp:
            gz = self._gz(tmp)
            with mock.patch.object(an, "_bgzf_eof_ok", return_value=True), \
                 mock.patch.object(an, "_reindex"), \
                 mock.patch.object(an.pysam, "TabixFile",
                                   side_effect=OSError("Invalid BGZF header")), \
                 mock.patch.object(an, "sleep") as slp, \
                 self.assertRaises(an.SamtoolsError) as ctx:
                an.load_annotation_bed(gz)
            self.assertIn("even after", str(ctx.exception))
            self.assertIn(str(an._ANNOT_LOAD_ATTEMPTS), str(ctx.exception))
            self.assertEqual(slp.call_count, an._ANNOT_LOAD_ATTEMPTS - 1)

    def test_vuoto_placeholder_never_indexed_or_consumed(self):
        """load_annotation_bed must NEVER index the empty-annotation placeholder vuoto.txt.
        pysam.tabix_index() bgzips in place and REMOVES the plain file, so indexing the
        placeholder consumes it -- which broke a later reference-only search's `cp vuoto.txt`
        (the serial pre-build in submit_job called this directly, bypassing annotation.py
        main()'s empty guard) and failed outright on a read-only install. It must return None
        (callers take the no-annotation path) and leave vuoto.txt untouched: no .gz / .tbi /
        .tbi.lock, and crucially pysam.tabix_index is never called."""
        with tempfile.TemporaryDirectory() as tmp:
            v = os.path.join(tmp, "vuoto.txt")
            open(v, "w").close()  # empty placeholder
            with mock.patch.object(an.pysam, "tabix_index",
                                   side_effect=AssertionError("must not index vuoto.txt")):
                self.assertIsNone(an.load_annotation_bed(v))
            self.assertTrue(os.path.isfile(v))                 # not consumed
            self.assertFalse(os.path.exists(v + ".gz"))        # not bgzipped
            self.assertFalse(os.path.exists(v + ".tbi"))
            self.assertFalse(os.path.exists(v + ".tbi.lock"))  # no lock written (read-only safe)


if __name__ == "__main__":
    unittest.main()

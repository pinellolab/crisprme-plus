"""Regression test: one indel target misattributed to EVERY chromosome.

The merged indel targets file is the only file in ``crispritz_targets/`` that
carries a header: ``pool_search_indels.py`` strips each per-contig header with
``tail -n +2`` while merging, then restores one with ``sed -i 1i``.
``post_analisi_indel.sh`` relies on that, copying line 1 of the merged file as
the header of each per-chromosome subset.

Two filters ran between those steps and both treated the header as a target:

  * ``pam_filter.py`` guarded headers with ``len(cols) < 3``, but the CRISPRitz
    header has all 10 tab-separated columns, so the guard never fired. The
    header's ``DNA`` column was read as a PAM region, ``_pam_can_match("DNA",
    "NRG")`` is False (``A`` vs ``G``), and the header was dropped as a
    "non-NRG target" -- the shipped log line even said ``dropped 1``.
  * the ``--max-total-edits`` cap in submit_job compared ``$10`` ("Total")
    against the numeric cap as strings, which is also false, dropping it too.

With the header gone, ``head -1`` returned the first DATA row, and
``sed -i 1i"$header"`` prepended that single row to ALL ~23 per-chromosome
subsets. Nothing downstream rejected it: every chromosome's fake genome starts
at coordinate 0 with near-uniform ~54 bp windows, so
``analisi_indels_NNN.py`` -- whose interval tree is keyed on the FAKE POSITION
ALONE -- resolved the foreign row against whichever chromosome it was running,
overwrote the Chromosome column with the ambient ``current_chr``, and stamped
that chromosome's real rsID, AF and carrier samples onto it.

One raw ``fakechr1`` search row therefore became 1 correct + N phantom report
rows on N different chromosomes, each looking like a genuine, well-annotated
variant off-target. Measured genome-wide on guide ACTGAAATCTGTAAGCAGGC at
NRG 4+1+1 over HGDP+1kGP: 8 phantoms on a dedicated NRG index (= published
2.8.0) and 11 on a K=8 pamless index, against an identical justified core of
13,253 indel rows. Which row gets injected depends on what lands at line 1, so
the count varied by index -- which is why the pamless/dedicated arms never
agreed.

These tests assert:
  * pam_filter.py keeps a ``#`` header while still dropping a genuinely
    non-matching target and keeping matching + IUPAC/PAM-creating ones;
  * the --max-total-edits cap awk keeps a ``#`` header while still capping;
  * prepending line 1 of a HEADERLESS file injects a foreign contig's row into
    another chromosome's subset (the bug), and the shipped guard refuses to;
  * the shipped scripts carry all four parts of the fix.

STDLIB only. Hermetic: no installed PAMs, index, or genome.

Run with:
    <env>/bin/python3 -m unittest discover -s PostProcess \
        -p 'test_indel_header_injection.py' -v
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
PAM_FILTER = os.path.join(PP, "pam_filter.py")
POST_INDEL = os.path.join(PP, "post_analisi_indel.sh")
ANALISI = os.path.join(PP, "analisi_indels_NNN.py")
SUBMIT = os.path.join(PP, "submit_job_automated_new_multiple_vcfs.sh")

# The real CRISPRitz header: 10 tab-separated columns, so a column-count guard
# does NOT catch it. Column 3 is "DNA", which is what got read as a PAM region.
HEADER = (
    "#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position\t"
    "Direction\tMismatches\tBulge Size\tTotal"
)
GUIDE = "ACTGAAATCTGTAAGCAGGCNNN"

# The row that was actually injected in the measured genome-wide run: a real
# fakechr1 hit whose IUPAC PAM region (GGG) satisfies NRG.
ROW_FAKECHR1 = (
    "RNA\t" + GUIDE + "\tgCrKAtATYTGT-RSYWGKkGGG\tfakechr1\t"
    "6100154\t6100153\t+\t4\t1\t5"
)
ROW_FAKECHR13 = (
    "X\t" + GUIDE + "\tGACTGACTGACTGACTGACTAGG\tfakechr13\t"
    "2000\t2000\t-\t1\t0\t1"
)
# PAM region TGC cannot satisfy NRG -> a genuine drop.
ROW_NON_NRG = (
    "X\t" + GUIDE + "\tGACTGACTGACTGACTGACTTGC\tfakechr13\t"
    "4000\t4000\t+\t0\t0\t0"
)
# PAM region NGG: variant-created / degenerate, must be KEPT (IUPAC-aware).
ROW_IUPAC_PAM = (
    "X\t" + GUIDE + "\tGACTGACTGACTGACTGACTNGG\tfakechr13\t"
    "5000\t5000\t+\t0\t0\t0"
)
# Total (col 10) = 9, above a cap of 6.
ROW_OVER_CAP = (
    "X\t" + GUIDE + "\tGACTGACTGACTGACTGACTCGG\tfakechr13\t"
    "6000\t6000\t+\t7\t2\t9"
)

NRG_PAM_FILE = "ACTGAAATCTGTAAGCAGGCNRG 3\n"  # 3' PAM, motif = last 3 = NRG


def _extract_shell_function(name):
    """Pull a tab-indented shell function verbatim out of submit_job.

    Extracted rather than copied, so the test fails if the shipped function and
    this file ever drift apart.
    """
    with open(SUBMIT, encoding="utf-8") as fh:
        text = fh.read()
    start = text.index(f"\t{name}() {{")
    end = text.index("\n\t}\n", start) + len("\n\t}\n")
    body = text[start:end]
    return "\n".join(
        line[1:] if line.startswith("\t") else line for line in body.splitlines()
    )


def _write(tmpdir, name, lines):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as fh:
        for ln in lines:
            fh.write(ln + "\n")
    return path


def _pam_file(tmpdir):
    path = os.path.join(tmpdir, "pam.txt")
    with open(path, "w") as fh:
        fh.write(NRG_PAM_FILE)
    return path


class PamFilterKeepsHeader(unittest.TestCase):
    def _run(self, lines):
        with tempfile.TemporaryDirectory() as td:
            targets = _write(td, "targets.txt", lines)
            r = subprocess.run(
                [sys.executable, "-I", PAM_FILTER, targets, _pam_file(td)],
                capture_output=True, text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(targets) as fh:
                return [ln.rstrip("\n") for ln in fh], r.stdout

    def test_header_survives_and_stays_first(self):
        """The 10-column header must pass through untouched, at line 1."""
        out, _ = self._run([HEADER, ROW_FAKECHR1, ROW_FAKECHR13])
        self.assertEqual(out[0], HEADER,
                         "pam_filter dropped the CRISPRitz header: downstream "
                         "post_analisi_indel.sh would then prepend a DATA row "
                         "to every chromosome's subset")
        self.assertEqual(len(out), 3, f"expected header + 2 rows, got {out!r}")

    def test_header_is_not_counted_as_a_kept_target(self):
        """The header must not inflate the kept/dropped tally."""
        _, log = self._run([HEADER, ROW_FAKECHR1])
        self.assertIn("kept 1", log)
        self.assertIn("dropped 0", log,
                      "a 'dropped 1' here is the header being discarded")

    def test_still_drops_a_genuinely_non_matching_pam(self):
        """The header fix must not turn pam_filter into a no-op."""
        out, log = self._run([HEADER, ROW_FAKECHR13, ROW_NON_NRG])
        self.assertEqual(out, [HEADER, ROW_FAKECHR13])
        self.assertIn("dropped 1", log)

    def test_keeps_iupac_and_variant_created_pam(self):
        """NGG/degenerate PAM regions stay (IUPAC-aware intersection)."""
        out, _ = self._run([HEADER, ROW_IUPAC_PAM, ROW_FAKECHR1])
        self.assertIn(ROW_IUPAC_PAM, out)
        self.assertIn(ROW_FAKECHR1, out)

    def test_headerless_input_is_unchanged(self):
        """Files that legitimately have no header still filter normally."""
        out, _ = self._run([ROW_FAKECHR13, ROW_NON_NRG])
        self.assertEqual(out, [ROW_FAKECHR13])


class EditsCapKeepsHeader(unittest.TestCase):
    """The --max-total-edits cap is the second path that dropped the header."""

    CAP_AWK = r"""awk -v m="$1" '/^#/ || NF < 10 || $10 <= m' "$2" """
    OLD_AWK = r"""awk -v m="$1" 'NF < 10 || $10 <= m' "$2" """

    def _run(self, snippet, cap, lines):
        with tempfile.TemporaryDirectory() as td:
            targets = _write(td, "targets.txt", lines)
            r = subprocess.run(["bash", "-c", snippet, "_", str(cap), targets],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            return r.stdout.splitlines()

    def test_fixed_awk_keeps_header_and_still_caps(self):
        out = self._run(self.CAP_AWK, 6, [HEADER, ROW_FAKECHR1, ROW_OVER_CAP])
        self.assertEqual(out[0], HEADER)
        self.assertIn(ROW_FAKECHR1, out)       # Total 5 <= 6
        self.assertNotIn(ROW_OVER_CAP, out)    # Total 9 > 6

    def test_old_awk_dropped_the_header(self):
        """Pins the defect: "Total" <= 6 is a false STRING comparison."""
        out = self._run(self.OLD_AWK, 6, [HEADER, ROW_FAKECHR1])
        self.assertNotIn(HEADER, out,
                         "if this now keeps the header, awk comparison "
                         "semantics changed and the fix needs revisiting")


class HeaderPrependIsGuarded(unittest.TestCase):
    """post_analisi_indel.sh must never prepend a DATA row to a subset."""

    # The shipped logic, equivalent in behaviour. The prepend itself is spelled
    # portably here ({ printf; cat; } instead of `sed -i 1i`) because BSD sed
    # reads `1i...` as -i's backup-extension argument, so the shipped form only
    # inserts under GNU sed and this suite also runs on macOS. What these tests
    # pin is the GUARD -- whether a data row is allowed through -- which is
    # independent of how the line is inserted. That the shipped script still
    # uses `sed -i 1i"$header"`, and only under a '#' case branch, is asserted
    # in ShippedScriptsCarryTheFix below.
    _PREPEND = r'{ printf "%s\n" "$header"; cat "$2"; } > "$2.n" && mv "$2.n" "$2"'
    GUARDED = r'''
header=$(head -1 "$1")
case "$header" in
    \#*) ''' + _PREPEND + r''' ;;
    *)   echo "WARNING: no header line; not prepending" ;;
esac
cat "$2"
'''
    UNGUARDED = r'''
header=$(head -1 "$1")
''' + _PREPEND + r'''
cat "$2"
'''

    def _run(self, snippet, merged_lines, subset_lines):
        with tempfile.TemporaryDirectory() as td:
            merged = _write(td, "merged.targets.txt", merged_lines)
            subset = _write(td, "merged.targets.txt.chr13", subset_lines)
            r = subprocess.run(["bash", "-c", snippet, "_", merged, subset],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            return r.stdout.splitlines()

    def test_unguarded_injects_a_foreign_contig_row(self):
        """The bug: a fakechr1 row lands in the fakechr13 subset."""
        out = self._run(self.UNGUARDED,
                        [ROW_FAKECHR1, ROW_FAKECHR13],   # headerless
                        [ROW_FAKECHR13])
        self.assertIn(ROW_FAKECHR1, out,
                      "expected the unguarded prepend to inject the foreign row")
        self.assertTrue(any("fakechr1\t" in ln for ln in out))

    def test_guarded_refuses_to_inject(self):
        out = self._run(self.GUARDED,
                        [ROW_FAKECHR1, ROW_FAKECHR13],   # headerless
                        [ROW_FAKECHR13])
        self.assertNotIn(ROW_FAKECHR1, out,
                         "a headerless merged file must not inject its first "
                         "data row into another chromosome's subset")
        self.assertFalse(any("fakechr1\t" in ln for ln in out))
        self.assertEqual([ln for ln in out if not ln.startswith("WARNING")],
                         [ROW_FAKECHR13])

    def test_guarded_still_prepends_a_real_header(self):
        out = self._run(self.GUARDED,
                        [HEADER, ROW_FAKECHR1, ROW_FAKECHR13],
                        [ROW_FAKECHR13])
        self.assertEqual(out[0], HEADER)
        self.assertIn(ROW_FAKECHR13, out)


class ShippedScriptsCarryTheFix(unittest.TestCase):
    def test_pam_filter_tests_for_a_comment_line(self):
        with open(PAM_FILTER) as fh:
            text = fh.read()
        self.assertRegex(
            text, r'row\.startswith\(\s*"#"\s*\)',
            "pam_filter.py must special-case '#' header lines; a column-count "
            "test alone does not catch the 10-column CRISPRitz header")

    def test_post_analisi_indel_guards_the_prepend(self):
        with open(POST_INDEL) as fh:
            text = fh.read()
        prepends = [ln for ln in text.splitlines()
                    if re.search(r'sed\s+-i\s+1i"\$header"', ln)]
        self.assertTrue(prepends, "expected the header prepend to still exist")
        self.assertRegex(
            text, r'case\s+"\$header"\s+in',
            "the prepend of $header must be guarded by a header test, or a "
            "dropped header injects a data row into every chromosome's subset")
        self.assertRegex(text, r'\\#\*\)',
                         "expected a '#'-prefix case branch guarding the prepend")

    def test_analisi_indels_rejects_foreign_contigs(self):
        with open(ANALISI) as fh:
            text = fh.read()
        self.assertRegex(
            text, r'if\s+line\[3\]\s*!=\s*"fake"\s*\+\s*current_chr',
            "analisi_indels_NNN.py must reject rows carrying another contig's "
            "label: its indel tree is keyed on the fake position alone, so a "
            "foreign row resolves to a REAL indel and is reported with that "
            "indel's rsID, AF and carrier samples")
        self.assertIn("foreign_rows", text,
                      "dropped foreign rows must be counted and reported")

    def test_edits_cap_awk_keeps_comment_lines(self):
        with open(SUBMIT) as fh:
            text = fh.read()
        cap = [ln for ln in text.splitlines()
               if "capped" in ln and "awk" in ln]
        self.assertTrue(cap, "expected the --max-total-edits cap awk")
        for ln in cap:
            self.assertIn("/^#/", ln,
                          f"cap awk must keep '#' header lines: {ln.strip()!r}")


class TargetRowCountIgnoresTheHeader(unittest.TestCase):
    """An empty result must not read as one hit.

    A CRISPRitz targets file carries a '#' header. While pam_filter was stripping
    it, `-s file` and `wc -l` were accidentally equivalent to counting data rows, so
    three call sites used them. With the header correctly preserved, a zero-hit
    search produced a header-only file, the zero-hit short-circuit saw a non-empty
    file, and a search that legitimately found nothing was reported as a failure --
    "the search found 1 raw off-target(s) but none survived post-analysis". Caught by
    the Cas12a brute-force benchmark, which has a genuinely empty chr22 result.
    """

    HEADER = ("#Bulge type\tcrRNA\tDNA\tChromosome\tPosition\tCluster Position"
              "\tDirection\tMismatches\tBulge Size\tTotal")
    ROW = "X\tACGT\tACGT\tchr22\t100\t100\t+\t0\t0\t0"

    def _rows(self, content):
        import subprocess, tempfile, os
        script = _extract_shell_function("_target_data_rows")
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write(content)
            path = fh.name
        self.addCleanup(os.unlink, path)
        out = subprocess.run(
            ["bash", "-c", script + '\n_target_data_rows "$1"\n', "bash", path],
            capture_output=True, text=True,
        )
        self.assertEqual(out.returncode, 0, out.stderr)
        return int(out.stdout.strip())

    def test_header_only_is_zero(self):
        self.assertEqual(self._rows(self.HEADER + "\n"), 0)

    def test_empty_file_is_zero(self):
        self.assertEqual(self._rows(""), 0)

    def test_header_plus_rows(self):
        self.assertEqual(self._rows(self.HEADER + "\n" + self.ROW + "\n"), 1)
        self.assertEqual(
            self._rows(self.HEADER + "\n" + self.ROW + "\n" + self.ROW + "\n"), 2
        )

    def test_headerless_file_still_counts_its_rows(self):
        """Robust either way: it must not blind-skip line 1."""
        self.assertEqual(self._rows(self.ROW + "\n"), 1)

    def test_missing_file_is_zero(self):
        import subprocess
        script = _extract_shell_function("_target_data_rows")
        out = subprocess.run(
            ["bash", "-c", script + '\n_target_data_rows "$1"\n', "bash",
             "/nonexistent/path/targets.txt"],
            capture_output=True, text=True,
        )
        self.assertEqual(out.stdout.strip(), "0")


class ShippedScriptCountsDataRows(unittest.TestCase):
    def test_no_site_sizes_a_targets_file_directly(self):
        with open(SUBMIT, encoding="utf-8") as fh:
            text = fh.read()
        self.assertNotIn('[ -s "$ref_targets_file" ]', text)
        self.assertNotIn('[ -s "$raw_ref_targets" ]', text)
        self.assertIn("_target_data_rows", text)

    def test_all_three_gates_use_the_helper(self):
        with open(SUBMIT, encoding="utf-8") as fh:
            text = fh.read()
        self.assertGreaterEqual(text.count("_target_data_rows "), 3)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""The search joins its three lists into one table line by line, with paste.

CONTEXT
-------
Each search writes three lists of its best alignments, one per criterion: CFD
(``.bestCFD.txt``), the fewest mismatches and bulges (``.bestmmblg.txt``) and CRISPR-Bulge
(``.bestCRISPR_BULGE.txt``), and three of the alternative alignments. Each list holds the
same clusters in the same order, one row each (``test_reference_rows_once``,
``test_representative_tiebreak``), and ``submit_job_automated_new_multiple_vcfs.sh`` joins
them into ``*.bestMerge.txt`` and ``*.altMerge.txt``: row N of the table is line N of the
three lists, joined by a tab, 3 x 24 columns.

It joined them with ``pr -m -t -J`` until 2.10.0. On a host whose ``pr`` is uutils'
(rust-coreutils 0.10.0, Ubuntu 26.04's default, which GitHub's ubuntu-latest runs since
November 2026), that ``pr`` fuses the lists with no tab between them and pads the table with
blank lines: 57 lines for three lists of 4 rows, each row with 2 columns fewer. The table
of a search run outside the release image on such a host (a Bioconda or source install) was
wrong. GNU ``pr`` changes rows too: it turns a run of two spaces or more into tabs where the
run crosses a tab stop, and takes a form feed for a new page, where it holds the rest of that
list back to the next page (up to 66 rows later), the table longer. ``paste`` copies each line
as it is, and GNU coreutils and uutils give the same bytes.

What such a table reached: the search stopped at "Integrating results", the stage that reads
it (docs/TROUBLESHOOTING.md section 17). The summaries and the plots of the stage before,
"Creating images", are made of the three lists, not of the table, and were right.

WHAT THIS TEST DOES
-------------------
Takes the two join lines from the REAL pipeline script and runs them with bash, on this
host's ``paste``: CI runs it on ubuntu-24.04 (GNU coreutils) and on ubuntu-26.04 (uutils).
``test_reference_rows_once`` and ``test_merge_sites_only_rows`` join their lists with these
lines too (``join_as_the_search_does``).

ASSERTS:
  * both joins are ``paste`` of the three lists, in their order, and no shell script of the
    search runs ``pr``;
  * row N is line N of each list joined by a tab, the line as it is: for rows of tab-
    separated columns, a header, a list shorter than the others or empty (its columns are
    empty), a line of 1 MiB, a last line with no newline, empty lines, and the bytes GNU
    ``pr`` changes (runs of spaces, a form feed, also in lists longer than its page) and a
    carriage return;
  * the summaries and the plots of "Creating images", which comes after the joins, are made of
    the three lists, and the stage names the table only to rename it, or to remove it on an
    error; the next stage, "Integrating results", reads it.
"""
import glob
import os
import re
import shutil
import subprocess
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
PIPELINE = os.path.join(PP, "submit_job_automated_new_multiple_vcfs.sh")
LISTS = ("bestCFD", "bestmmblg", "bestCRISPR_BULGE")


def _uncommented(path):
    with open(path) as fh:
        return [ln for ln in fh.read().splitlines() if not ln.lstrip().startswith("#")]


def pipeline_join_lines():
    """The search's two joins, of the best and of the alternative alignments, as the script
    has them."""
    picked = [ln for ln in _uncommented(PIPELINE)
              if re.search(r'>\s*"?\$final_res(_alt)?"?\s*$', ln) and ".bestCFD.txt" in ln]
    assert len(picked) == 2, "the search's joins of its three lists moved: %r" % picked
    return picked


def join_as_the_search_does(paths, alt=False):
    """The table the search makes of three files (the CFD, fewest mismatches and bulges and
    CRISPR-Bulge lists, in that order), with its own join line run by bash; its bytes."""
    var = "final_res_alt" if alt else "final_res"
    line = [ln for ln in pipeline_join_lines() if ("$final_res_alt" in ln) == alt][0]
    with tempfile.TemporaryDirectory(prefix="join_lists_") as d:
        final = os.path.join(d, "final")
        for name, path in zip(LISTS, paths):
            shutil.copyfile(path, "%s.%s.txt" % (final, name))
        subprocess.run(["bash", "-c", line], check=True, capture_output=True,
                       env=dict(os.environ, **{var: final}))
        with open(final, "rb") as fh:
            return fh.read()


def joined(*bodies):
    """Row N: line N of each body (bytes), joined by a tab; a body with fewer lines gives an
    empty field. A last line with no newline is a line."""
    lists = []
    for body in bodies:
        lines = body.split(b"\n")
        if lines[-1] == b"":
            lines.pop()
        lists.append(lines)
    rows = max(len(lines) for lines in lists)
    return b"".join(b"\t".join(lines[i] if i < len(lines) else b"" for lines in lists) + b"\n"
                    for i in range(rows))


def _row(i, tag):
    """A row of a list as the search writes it: 24 tab-separated columns."""
    return "\t".join(["X", "CTAACAGTTGCTTTTATCACNNN", "cTAACAGTTGCTTTTATCACAGG",
                      "CTAACAGTTGCTTTTATCACAGG", "chr22", str(10000 + i), str(10000 + i), "+",
                      "1", "0", "1", "NGG", "n", "HG00097,HG00099", "CTCF-only,CTCF-bound",
                      "CTAACAGTTGCTTTTATCACNNN", "rs123", "0.25", "chr22_10011_T_A", "2",
                      "0.4211", "0.1172", "0.3039", tag + str(i)]).encode()


def _rows(n, tag):
    return b"".join(_row(i, tag) + b"\n" for i in range(n))


class TheSearchJoinsItsListsWithPaste(unittest.TestCase):

    def test_both_joins_are_paste_of_the_three_lists_in_order(self):
        for line in pipeline_join_lines():
            with self.subTest(line=line):
                self.assertRegex(
                    line.strip(),
                    r'^paste "\$(final_res(?:_alt)?)\.bestCFD\.txt" "\$\1\.bestmmblg\.txt" '
                    r'"\$\1\.bestCRISPR_BULGE\.txt" >"\$\1"$')

    def test_no_shell_script_of_the_search_runs_pr(self):
        # pr -m -t -J joined the lists until 2.10.0: uutils' pr pads and fuses them
        command = re.compile(r"(?:^|[;&|(`{]|\$\(|\b(?:then|do|else|xargs|exec)\s)\s*pr\s")
        found = [(os.path.basename(path), ln.strip())
                 for path in sorted(glob.glob(os.path.join(PP, "*.sh")))
                 for ln in _uncommented(path) if command.search(ln)]
        self.assertEqual(found, [])


class PasteJoinsLineNOfEachList(unittest.TestCase):
    """Run on this host's paste: CI runs it with GNU coreutils and with uutils."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="join_cases_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def check(self, a, b, c, alt=False):
        paths = []
        for name, body in zip(LISTS, (a, b, c)):
            paths.append(os.path.join(self.d, name))
            with open(paths[-1], "wb") as fh:
                fh.write(body)
        table = join_as_the_search_does(paths, alt=alt)
        self.assertEqual(table, joined(a, b, c))
        return table

    def test_rows_of_tab_separated_columns(self):
        for alt in (False, True):
            with self.subTest(alt=alt):
                table = self.check(_rows(4, "a"), _rows(4, "b"), _rows(4, "c"), alt=alt)
                rows = table.decode().splitlines()
                self.assertEqual(len(rows), 4, "one row per line of the lists, no padding")
                self.assertEqual({len(r.split("\t")) for r in rows}, {72}, "3 x 24 columns")

    def test_a_header_and_more_rows_than_a_page_of_pr(self):
        header = ("\t".join("C%d" % i for i in range(24)) + "\n").encode()
        table = self.check(header + _rows(100, "a"), header + _rows(100, "b"),
                           header + _rows(100, "c"))
        self.assertEqual(len(table.splitlines()), 101)

    def test_a_list_shorter_than_the_others_or_empty(self):
        for lists in ((_rows(2, "a"), _rows(5, "b"), _rows(5, "c")),
                      (_rows(5, "a"), _rows(2, "b"), _rows(5, "c")),
                      (_rows(5, "a"), _rows(5, "b"), _rows(2, "c")),
                      (b"", _rows(3, "b"), _rows(3, "c")),
                      (_rows(3, "a"), b"", _rows(3, "c")),
                      (_rows(3, "a"), _rows(3, "b"), b""),
                      (b"", b"", b"")):
            with self.subTest(lines=[body.count(b"\n") for body in lists]):
                self.check(*lists)

    def test_a_long_line(self):
        long_line = b"\t".join(b"L" * 4096 for _ in range(256))   # 1 MiB
        self.check(long_line + b"\n" + _row(1, "a") + b"\n", _row(1, "b") + b"\n" + long_line + b"\n",
                   long_line + b"\n")

    def test_a_last_line_with_no_newline(self):
        self.check(_rows(3, "a").rstrip(b"\n"), _rows(3, "b"), _rows(3, "c"))
        self.check(_rows(3, "a"), _rows(3, "b"), _rows(3, "c").rstrip(b"\n"))
        self.check(_rows(1, "a").rstrip(b"\n"), _rows(3, "b"), _rows(2, "c").rstrip(b"\n"))

    def test_empty_lines(self):
        self.check(b"a1\n\na3\n", b"\nb2\n\n", b"c1\nc2\n\n")

    def test_each_line_as_it_is(self):
        # GNU pr turned a run of two spaces or more that crosses a tab stop into tabs: it made
        # "abcdef  gh" "abcdef<TAB>gh", a column more, and "a  b<TAB>        c"
        # "a  b<TAB><TAB>    c"
        self.check(b"abcdef  gh\n", b"b1\n", b"c1\n")
        self.check(b"a  b\t        c\n", b"two  spaces\n", b"eight        spaces\n")
        self.check(b"a\r\n\xc3\xa9\n", b"b\r\n\n", b"c\n\r\n")

    def test_a_form_feed(self):
        # GNU pr took a form feed for a new page and held the rest of that list back to the
        # next page, the table longer: here row 2 came out as row 4, of 5, and in the lists of
        # 100 rows, longer than its page of 66 lines, as row 67, of 165
        self.check(b"a1\n\fa2\na3\n", b"b1\nb2\nb3\n", b"c1\nc2\nc3\n")
        first = _rows(100, "a").split(b"\n")
        first[1] = b"\f" + first[1]
        table = self.check(b"\n".join(first), _rows(100, "b"), _rows(100, "c"))
        self.assertEqual(len(table.splitlines()), 100)


# the joined tables: $final_res and $final_res_alt with no list's suffix, and the names they take
TABLE = re.compile(r"\$\{?final_res(?:_alt)?\}?(?![\w.])|\.(?:best|alt)Merge\.txt\b")


class TheSummariesAndPlotsAreMadeOfTheListsNotOfTheTable(unittest.TestCase):
    """On a host whose pr was uutils', a search stopped at "Integrating results", where BEDOPS'
    sort-bed reads the table, and the summaries and the plots it had made were right: the stage
    that makes them, "Creating images", comes after the joins and reads the three lists
    (docs/TROUBLESHOOTING.md section 17 and the CHANGELOG say so)."""

    @classmethod
    def setUpClass(cls):
        with open(PIPELINE) as fh:
            cls.lines = fh.read().splitlines()

    def stage(self, name):
        """The numbers of the lines of the stage's Start and End marks."""
        marks = [[i for i, ln in enumerate(self.lines)
                  if ln.startswith("echo -e '%s\\t%s\\t'" % (name, edge))]
                 for edge in ("Start", "End")]
        self.assertEqual([len(found) for found in marks], [1, 1], name)
        return marks[0][0], marks[1][0]

    def commands(self, name):
        start, end = self.stage(name)
        return [ln.strip() for ln in self.lines[start:end] if not ln.lstrip().startswith("#")]

    def test_it_comes_after_the_joins_and_before_the_integration(self):
        joins = [i for i, ln in enumerate(self.lines) if ln in pipeline_join_lines()]
        self.assertEqual(len(joins), 2)
        images, integration = self.stage("Creating images"), self.stage("Integrating results")
        self.assertLess(max(joins), images[0])
        self.assertLess(images[1], integration[0])

    def test_it_names_the_table_only_to_rename_it_or_to_remove_it(self):
        named = [ln for ln in self.commands("Creating images") if TABLE.search(ln)]
        self.assertTrue(any(ln.startswith("mv $final_res ") for ln in named), named)
        self.assertEqual([ln for ln in named if not ln.startswith(("mv ", "rm "))], [])

    def test_its_summaries_and_plots_read_the_three_lists(self):
        runs = [ln for ln in self.commands("Creating images") if re.search(r"\.py\b", ln)]
        lists = {}
        for ln in runs:
            script = re.search(r"(\w+)\.py\b", ln).group(1)
            read = re.findall(r"\$final_res\.(best\w+)\.txt", ln)
            # populations_distribution.py reads what process_summaries.py wrote of the lists
            if not read:
                self.assertRegex(ln, r"\.PopulationDistribution_[\w-]+\.txt")
            lists.setdefault(script, set()).update(read)
        self.assertEqual(lists, {"process_summaries": set(LISTS),
                                 "populations_distribution": set(),
                                 "radar_chart_dict_generator": set(LISTS)})

    def test_the_integration_reads_the_table(self):
        self.assertTrue(any(
            "post_process.sh" in ln and "$(basename ${output_folder}).bestMerge.txt" in ln
            for ln in self.commands("Integrating results")))


if __name__ == "__main__":
    unittest.main()

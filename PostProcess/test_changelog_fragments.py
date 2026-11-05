#!/usr/bin/env python3
"""A changelog fragment has no heading but the sections it is merged into.

A release in the making keeps its CHANGELOG entries in changelog/<track>.md, which
merge_changelog.py merges into CHANGELOG.md's section for the release. It drops HTML comments,
leaves the text before the first heading out (it says that text is there, for whoever merges),
and appends each "### Added", "### Changed", "### Deprecated", "### Removed", "### Fixed" or
"### Security" block, line for line, to that subsection. Any other heading after the first is
text of the block above it: changelog/docs3.md ended with "### Notes for the merge", notes for
whoever merges, and the merge put them in CHANGELOG.md's [2.8.1], as a subsection between
Changed and Fixed, with no warning, while the amendment they asked for was never made. Notes for
the merge go before the first heading. A "### " heading of one word that is no section, or a
"## " heading, stops the merge.

The scan reads the fragments there are: none once a release has merged them. Fenced code and
HTML comments are not read for headings.
"""

import pathlib
import re
import unittest

REPO = pathlib.Path(__file__).resolve().parent.parent

SECTIONS = ("Added", "Changed", "Deprecated", "Removed", "Fixed", "Security")
HEADING = re.compile(r"#{1,6}(?:[ \t]|$)")
FENCE = re.compile(r"[ \t]*(`{3,}|~{3,})")


def misplaced_headings(text):
    """The numbers (from 1) and text of the headings of a fragment that the merge does not take
    for a section: after the first section, any heading but a section's; anywhere, a "## "
    heading or a one-word "### " heading that is no section, which stop it."""
    # blank the comments, keeping their lines, so that a number is the file's
    text = re.sub(r"<!--.*?-->", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S)
    found, fence, started = [], None, False
    for n, line in enumerate(text.split("\n"), 1):
        m = FENCE.match(line)
        if fence:
            if (m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence)
                    and not line.strip().strip(fence[0])):
                fence = None
            continue
        if m:
            fence = m.group(1)
            continue
        if not HEADING.match(line):
            continue
        word = re.fullmatch(r"### (\w+)\s*", line)
        if word and word.group(1) in SECTIONS:
            started = True
        elif word or line.startswith("## ") or started:
            found.append((n, line))
    return found


def fragments():
    return sorted((REPO / "changelog").glob("*.md"))


class TestAFragmentHasOnlyTheSectionsItIsMergedInto(unittest.TestCase):

    def test_it_finds_notes_after_the_sections(self):
        # changelog/docs3.md at 8cd4f4b, in short
        broken = ("<!-- docs3 track: ... The last\nsection holds notes for the merge -->\n\n"
                  "### Fixed\n\n- **An entry.**\n\n### Changed\n\n- **Another.**\n\n"
                  "### Notes for the merge\n\n- Amend the entry ...\n")
        self.assertEqual(misplaced_headings(broken), [(12, "### Notes for the merge")])
        self.assertEqual(misplaced_headings("### Fixed\n\n- a\n\n#### More\n\n# Notes\n"),
                         [(5, "#### More"), (7, "# Notes")])

    def test_but_not_the_same_notes_before_them(self):
        fixed = ("<!-- docs3 track: ... -->\n\n**For the merge, not for the changelog: two notes.**"
                 "\n\n1. Already corrected ...\n\n### Fixed\n\n- **An entry.**\n\n"
                 "### Changed\n\n- **Another.**\n")
        self.assertEqual(misplaced_headings(fixed), [])
        # the title some fragments open with
        self.assertEqual(misplaced_headings("# [2.8.1] entries from the asm track\n\n"
                                            "### Fixed\n\n- a\n"), [])

    def test_a_heading_that_stops_the_merge_is_found_anywhere(self):
        self.assertEqual(misplaced_headings("## [2.8.1]\n\n### Fixed\n- a\n"), [(1, "## [2.8.1]")])
        self.assertEqual(misplaced_headings("### Notes\n\n### Fixed\n- a\n"), [(1, "### Notes")])
        self.assertEqual(misplaced_headings("### Fixed\n- a\n### Fixes\n- b\n"), [(3, "### Fixes")])

    def test_comments_and_code_are_not_read_but_what_follows_them_is(self):
        self.assertEqual(misplaced_headings("### Fixed\n- a\n<!--\n### Notes for the merge\n-->\n"
                                            "- b\n"), [])
        self.assertEqual(misplaced_headings("### Fixed\n- a:\n\n```bash\n# 4b. the checks\n```\n"
                                            "#### After\n"), [(7, "#### After")])
        # a comment's lines are kept, so a number is the file's
        self.assertEqual(misplaced_headings("<!-- a\nb\nc -->\n### Fixed\n- a\n# Notes\n"),
                         [(6, "# Notes")])

    def test_every_fragment_has_only_its_sections(self):
        offenders = [f"{p.relative_to(REPO)}:{n}: {line}" for p in fragments()
                     for n, line in misplaced_headings(p.read_text(encoding="utf-8"))]
        self.assertEqual(offenders, [], "headings the merge does not take for a section of the "
                                        "release (notes for the merge go before the first "
                                        "heading): " + "; ".join(offenders))


if __name__ == "__main__":
    unittest.main()

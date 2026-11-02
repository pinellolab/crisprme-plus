#!/usr/bin/env python3
"""A wrapped line of the Markdown the release ships never opens a list.

CommonMark, and so GitHub, reads a line that starts with "+ " as a bullet, also in the
middle of a paragraph or of a list item. The --prune-pam entry of changelog/idx.md was
wrapped "... is reported on the" / "+ strand at another position; ...": the rest of the
entry rendered as a nested bullet, without the "+" that says which strand, and
merge_changelog.py copies a fragment into CHANGELOG.md line for line.
docs/PRECOMPUTED_INDEXES.md had such a line too, and seven older lines of docs/ did
("+ bulges)", "+ 1 RNA bulge)", ...). This repo never uses "+" as a bullet, so a line that
starts with one is a wrapped "plus".

The scan reads the README, METHODS, docs/ (the Hugging Face cards included) and the
changelog fragments. It skips fenced code blocks and HTML comments, which are not read as
Markdown text.
"""

import pathlib
import re
import unittest

REPO = pathlib.Path(__file__).resolve().parent.parent

# A "+" and then a space or a tab opens a list item at any depth of a list, so the
# indentation does not matter.
PLUS_BULLET = re.compile(r"[ \t]*\+[ \t]")
FENCE = re.compile(r"[ \t]*(`{3,}|~{3,})")


def plus_bullets(text):
    """The numbers (from 1) of the lines of ``text`` that start with a "+" bullet, outside
    fenced code blocks and HTML comments."""
    found, fence, comment = [], None, False
    for n, line in enumerate(text.split("\n"), 1):
        if comment:
            comment = "-->" not in line
            continue
        m = FENCE.match(line)
        if fence:
            if (m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence)
                    and not line.strip().strip(fence[0])):
                fence = None
            continue
        if m:
            fence = m.group(1)
        elif line.lstrip().startswith("<!--"):
            comment = "-->" not in line
        elif PLUS_BULLET.match(line):
            found.append(n)
    return found


def shipped_markdown():
    """The README, METHODS, docs/ and the changelog fragments."""
    paths = [REPO / "README.md", REPO / "METHODS.md"]
    paths += sorted((REPO / "docs").rglob("*.md"))
    paths += sorted((REPO / "changelog").glob("*.md"))
    return paths


class TestNoWrappedLineOpensAPlusBullet(unittest.TestCase):

    def test_the_scan_reads_the_docs_and_the_cards(self):
        """Guard the scan itself: if it reads nothing it proves nothing."""
        names = {str(p.relative_to(REPO)) for p in shipped_markdown()}
        for name in ("README.md", "METHODS.md", "docs/PRECOMPUTED_INDEXES.md",
                     "docs/hf/crisprme-data-general-index_README.md"):
            self.assertIn(name, names)

    def test_it_finds_the_wrapped_plus_strand(self):
        broken = ("- **Entry.** A - strand PAM K - L bases from the start is reported on the\n"
                  "  + strand at another position; K = L is unaffected.\n")
        self.assertEqual(plus_bullets(broken), [2])
        self.assertEqual(plus_bullets("A paragraph that ends with a\n+ bulges). Then\n"), [2])
        fixed = ("- **Entry.** A - strand PAM K - L bases from the start is reported on\n"
                 "  the + strand at another position; K = L is unaffected.\n")
        self.assertEqual(plus_bullets(fixed), [])
        self.assertEqual(plus_bullets("4 mismatches + 1 DNA +\n1 RNA bulge, or\n+1 RNA bulge\n"), [])

    def test_code_and_comments_are_not_read_but_what_follows_them_is(self):
        self.assertEqual(plus_bullets("```diff\n+ added\n```\n+ after\n"), [4])
        self.assertEqual(plus_bullets("  ```\n  + in a list\n  ```\n"), [])
        self.assertEqual(plus_bullets("````\n```\n+ still code\n````\n+ after\n"), [5])
        self.assertEqual(plus_bullets("<!-- a note\n+ more\n-->\n+ after\n"), [4])
        self.assertEqual(plus_bullets("<!-- one line -->\n+ after\n"), [2])

    def test_no_line_starts_with_a_plus_bullet(self):
        offenders = [f"{p.relative_to(REPO)}:{n}" for p in shipped_markdown()
                     for n in plus_bullets(p.read_text(encoding="utf-8"))]
        self.assertEqual(offenders, [], "lines that start with '+ ', which CommonMark reads "
                                        "as a bullet: " + ", ".join(offenders))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Tests for make_oh_sidecar.py -- the observed-haplotypes sidecar generator.

The sidecar keeps only alt-origin rows (highest-CFD OR highest-CRISPR_BULGE axis),
copied VERBATIM, as a strict superset of what the report's companion keeps.
"""
import os
import tempfile
import unittest

import make_oh_sidecar as m

_H = (
    ["Spacer+PAM", "Chromosome", "Start_coordinate_(highest_CFD)", "Strand_(highest_CFD)"]
    + [f"c{i}" for i in range(5, 17)]
    + ["REF/ALT_origin_(highest_CFD)"]
    + [f"d{i}" for i in range(18, 69)]
    + ["REF/ALT_origin_(highest_CRISPR_BULGE)"]
)


def _row(cfd_origin, cb_origin, tag="x"):
    cells = ["SP", "chr1", "100", "+"] + [tag] * 12 + [cfd_origin] + [tag] * 51 + [cb_origin]
    return "\t".join(cells)


class TestMakeOhSidecar(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def _run(self, lines):
        a = os.path.join(self.d, "alt.tsv")
        o = os.path.join(self.d, "side.tsv")
        with open(a, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        total, kept = m.write_sidecar(a, o)
        with open(o) as fh:
            out = [ln.rstrip("\n") for ln in fh]
        return total, kept, out, lines

    def test_keeps_alt_origin_either_axis_byte_identical(self):
        hdr = "\t".join(_H)
        src = [hdr,
               _row("alt", "ref", "a"),   # kept (CFD axis alt)
               _row("ref", "ref", "b"),   # dropped (neither alt)
               _row("ref", "alt", "c"),   # kept (CB axis alt)
               _row("ref", "ref", "d")]   # dropped
        total, kept, out, lines = self._run(src)
        self.assertEqual((total, kept), (4, 2))
        self.assertEqual(out[0], hdr)                 # header preserved verbatim
        self.assertEqual(out[1], lines[1])            # kept rows are BYTE-identical to source
        self.assertEqual(out[2], lines[3])
        self.assertEqual(len(out), 3)                 # header + 2 kept rows

    def test_case_and_space_insensitive_origin(self):
        hdr = "\t".join(_H)
        _, kept, _, _ = self._run([hdr, _row(" ALT ", "ref"), _row("Ref", " alt")])
        self.assertEqual(kept, 2)                     # both kept despite case/space

    def test_missing_origin_columns_copies_all(self):
        # a header WITHOUT the origin columns -> safe degrade: copy every row
        hdr = "\t".join(["Spacer+PAM", "Chromosome", "foo", "bar"])
        src = [hdr, "SP\tchr1\t1\t2", "SP\tchr2\t3\t4"]
        total, kept, out, _ = self._run(src)
        self.assertEqual((total, kept), (2, 2))
        self.assertEqual(out, src)                    # verbatim copy

    def test_empty_file(self):
        a = os.path.join(self.d, "empty.tsv")
        o = os.path.join(self.d, "o.tsv")
        open(a, "w").close()
        self.assertEqual(m.write_sidecar(a, o), (0, 0))


if __name__ == "__main__":
    unittest.main()

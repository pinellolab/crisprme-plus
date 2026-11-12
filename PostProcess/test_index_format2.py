#!/usr/bin/env python3
"""CRISPRme reads index format 2, the format CRISPRitz writes from 2.10.0.

An index of format 2 (CRISPRitz's docs/INDEX_FORMAT.md) opens with a header -- the magic
"CZI" 0xFE, the format version, the header's length, the flags, the guide length, the PAM
length, the max bulges and the PAM -- and ends with the magic again; each target carries a
flags byte where format 1 has a '0' or '_' marker; and it holds sites format 1 cannot: those
whose guide lies within a few bases of a contig end or an N run (the pad may run off the
contig or into an N), a PAM ending on a contig's last base, and, PAM-less, a flank partly off
the contig. CRISPRme 2.8.1 read none of it:

* crispritz_bins.bin_problem called every format-2 .bin "not a tree (numLeaves=-28747197)",
  and _leaves returned the magic, so the full-part check never ran;
* holds_target, count_targets and missing_bins applied the site rule of format 1 to any
  index, which on a format-2 index passes a lost .bin off as a contig with no target; the
  rule must be that of the format the index was built in, read from its own headers, not
  from the installed CRISPRitz; an index with no .bin records its format;
* index_fit read no guide length from a format-2 index, index_coverage checked nothing on
  one, and every search refused one ("this CRISPRme reads index format 1 only"); a format-2
  index must be refused only where the installed CRISPRitz cannot read it (2.9.1), naming
  both, and an index scan must prefer one it can read over one it would stop on;
* each index records what it holds (.index_meta.json: index format, build PAM, side, flank,
  guide length).

The trees here are made as buildTST 2.10.0 makes them (save_tst2, a port of saveTST), and the
site rule is checked against buildTST's own loops (buildtst2_targets, a port of
searchPAMonGenome and finishLeaf). test/index_format2/check_real_engine.py checks both
against the real engine, where one is installed (not in this job: run it locally).

Stdlib only.
"""

import json
import os
import random
import re
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import crispritz_bins  # noqa: E402
import index_coverage  # noqa: E402
import index_fit  # noqa: E402
import prune_flank  # noqa: E402
import test_build_index_only_atomic as atomic  # noqa: E402  (its stub engine and fixture)
import test_dictless_publish as dpub  # noqa: E402  (publish and download, offline)
import test_indel_empty_chromosome as empty  # noqa: E402  (its stub engine and harnesses)

MAGIC = b"CZI\xfe"
# the four bits a tree or a flank gives a base (writePair, saveTST); Q is 1111 in a flank
CODE = dict(zip("ACGTRYSWKMBDHVNQ", (1, 2, 4, 8, 5, 10, 6, 9, 12, 3, 14, 13, 11, 7, 0, 15)))


def save_tst2(targets, pam, guide_len=20, max_bulges=2, at_start=False, q=False):
    """A .bin of index format 2 as buildTST 2.10.0 writes one (insall, insert, saveTST and
    serialize in CRISPRitz's mainParallel.cpp, ported): ``targets`` are (window start,
    sequence, flank, leaf flags), sorted by sequence, all sequences of one length; the flags
    are PAM_LEFT 0x40 | FLANK_OFF 0x20 | pad_valid, HAS_NEXT being added here. ``q``: the
    file is flagged as holding the four-allele symbol Q."""
    nodes, nexts = [], [0] * len(targets)  # a node: [base, lokid, hikid, eqkid]

    def insert(seq, i):
        k, p = 0, 0
        while nodes:
            node = nodes[p]
            if seq[k] == node[0]:
                k += 1
                if k == len(seq):
                    nexts[i], node[3] = node[3], -(i + 1)
                    return
                if not node[3]:
                    node[3] = len(nodes)
                    break
                p = node[3]
            else:
                kid = 1 if seq[k] < node[0] else 2
                if not node[kid]:
                    node[kid] = len(nodes)
                    break
                p = node[kid]
        while True:
            nodes.append([seq[k], 0, 0, 0])
            k += 1
            if k == len(seq):
                nodes[-1][3] = -(i + 1)
                return
            nodes[-1][3] = len(nodes)

    def insall(lo, hi):
        if hi >= lo:
            mid = lo + (hi - lo) // 2
            insert(targets[mid][1], mid)
            insall(lo, mid - 1)
            insall(mid + 1, hi)

    insall(0, len(targets) - 1)
    flags = (1 if at_start else 0) | (2 if set(pam) == {"N"} else 0) | (4 if q else 0)
    out = bytearray(MAGIC + struct.pack("=6i", 2, 28 + len(pam), flags, guide_len, len(pam),
                                        max_bulges) + pam.encode())
    out += struct.pack("=ii", len(targets), guide_len)
    for i, (ws, _, flank, lflags) in enumerate(targets):
        out += struct.pack("=i", ws)
        nibble, k = 0, 0
        for j, base in enumerate(flank):
            nibble, k = nibble + CODE[base], k + 1
            if j == len(flank) - 1 or k == 2:
                out.append(nibble)
                nibble, k = 0, 0
            if k == 1:
                nibble <<= 4
        out.append(lflags | (0x80 if nexts[i] else 0))
        if nexts[i]:
            out += struct.pack("=i", nexts[i])
    out += struct.pack("=i", len(nodes))
    pair, first = ["0", "0"], [True]
    code = dict(CODE, **{"0": 0, "_": 15})

    def put(base):
        pair[0 if first[0] else 1] = base
        if base == "_" or not first[0]:
            out.append(0xF0 if pair[0] == "_" else code[pair[0]] << 4 | code[pair[1]])
            first[0] = True
        else:
            first[0] = False

    def serialize(p):
        base, lo, hi, eq = nodes[p]
        put(base)
        serialize(lo) if lo > 0 else put("0")
        serialize(hi) if hi > 0 else put("0")
        if eq > 0:
            serialize(eq)
        else:
            put("_")
            out.extend(struct.pack("=i", eq))

    serialize(0)
    return bytes(out + MAGIC)


def one_target_bin(path, pam, guide_len=20, max_bulges=2, at_start=False, fmt=2, leaves=1):
    """A whole .bin of one target (or of ``leaves`` alike), in index format ``fmt``."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if fmt == 2:
        data = save_tst2([(5, "ACGT", pam.replace("N", "A"), 0x02)] * leaves, pam, guide_len,
                         max_bulges, at_start)
    else:
        data = struct.pack("=ii", leaves, guide_len) + b"\0" * 8
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def revcomp_pam(pam):
    return pam.translate(crispritz_bins._REVERSE)[::-1]


def buildtst2_targets(seq, pam, pam_len, at_start):
    """The number of targets buildTST keeps in contig seq in index format 2:
    searchPAMonGenome (pamGenerator.cpp), then siteWindow and finishLeaf (mainParallel.cpp)
    of CRISPRitz 2.10.0, loop for loop, as slow as they are."""
    k, n = len(pam), len(seq)
    g = pam_len - k
    if g < 1 or k < 1 or n < g + 1:
        return 0
    right_ok = lambda nt: g <= nt <= n - 1  # noqa: E731 - as the engine's lambdas
    left_ok = lambda nt: nt + k - 1 >= 0 and nt + k + g <= n  # noqa: E731
    sites = []
    plus = [crispritz_bins._PAM_BITS[b] for b in pam]
    if all(b == 15 for b in plus):
        for nt in range(-(k - 1), n):
            for shape in (("R", "L") if not at_start else ("L", "R")):
                if (right_ok if shape == "R" else left_ok)(nt):
                    sites.append((shape, nt))
    else:
        if n < k:
            return 0
        genome = [crispritz_bins._BITS.get(b, 0) for b in seq]
        minus = [crispritz_bins._PAM_BITS[b] for b in revcomp_pam(pam)]
        for nt in range(0, n - k + 1):
            fp = all(genome[nt + j] & plus[j] for j in range(k))
            fn = all(genome[nt + j] & minus[j] for j in range(k))
            if not at_start:
                sites += [("R", nt)] * (fp and right_ok(nt)) + [("L", nt)] * (fn and left_ok(nt))
            else:
                sites += [("L", nt)] * (fp and left_ok(nt)) + [("R", nt)] * (fn and right_ok(nt))
    kept = 0
    for shape, nt in sites:
        start = nt - g if shape == "R" else nt + k
        guide = "".join(seq[p] if 0 <= p < n else "N" for p in range(start, start + g))
        kept += "N" not in guide
    return kept


def stub_engine(folder, version="2.10.0", index_format=2):
    """A crispritz.py on a PATH folder of its own that says ``version`` and, unless None,
    INDEX_FORMAT ``index_format``, as the real one does (CRISPRitz 2.9.1 states none)."""
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "crispritz.py")
    with open(path, "w") as fh:
        fh.write("#!/usr/bin/env python3\n# stub\n")
        fh.write(f'VERSION = "{version}"\n')
        if index_format is not None:
            fh.write(f"INDEX_FORMAT = {index_format}\n")
        fh.write("import sys\nsys.exit(0)\n")
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return folder


class _Engine(unittest.TestCase):
    """A temporary folder, and an installed CRISPRitz of a given index format on PATH."""

    def setUp(self):
        self.tmp = os.path.realpath(tempfile.mkdtemp(prefix="fmt2_"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.addCleanup(os.environ.__setitem__, "PATH", os.environ.get("PATH", ""))

    def engine(self, index_format=2, version=None):
        """Put on PATH a CRISPRitz whose crispritz.py states ``index_format`` (None: none,
        as 2.9.1's), or none at all with ``index_format`` False: then a folder of PATH that
        holds the installed one (the Docker image's /opt/conda/bin, a source or a Bioconda
        install's) gives way to one of links to all else it holds, so that bash, sed and
        python are found as before."""
        path = [p for p in os.environ.get("PATH", "").split(os.pathsep)
                if not p.startswith(os.path.join(self.tmp, "engine"))]
        if index_format is not False:
            folder = os.path.join(self.tmp, "engine%s" % len(os.listdir(self.tmp)))
            stub_engine(folder, version or ("2.9.1" if index_format is None else "2.10.0"),
                        index_format)
            path.insert(0, folder)
        else:
            path = [self.without_crispritz(p) if shutil.which("crispritz.py", path=p) else p
                    for p in path]
        os.environ["PATH"] = os.pathsep.join(path)
        return os.environ["PATH"]

    def without_crispritz(self, folder):
        """A folder of links to all ``folder`` holds but its crispritz.py."""
        links = tempfile.mkdtemp(prefix="bin_without_crispritz_", dir=self.tmp)
        for name in os.listdir(folder):
            if name != "crispritz.py":
                os.symlink(os.path.join(os.path.abspath(folder), name), os.path.join(links, name))
        return links


class TheFormat2TreeIsTold(unittest.TestCase):
    """bin_problem on trees as buildTST 2.10.0 saves them: of many targets, some of one
    sequence, with shared starts and IUPAC bases, for PAMs of odd and even length, dedicated
    and PAM-less, 3' and 5'. Each is whole, and each cut of it, at every length, is not."""

    def test_every_cut_of_a_whole_format_2_tree_is_told_apart(self):
        rng = random.Random(29)
        tmp = tempfile.mkdtemp(prefix="fmt2_bins_")
        self.addCleanup(shutil.rmtree, tmp, True)
        for pam, n, at_start in (("NRG", 40, False), ("TTTV", 25, True), ("NNNNNNNN", 3, False),
                                 ("NGG", 2, False), ("NNNN", 9, True), ("NGG", 1, False)):
            with self.subTest(pam=pam, targets=n):
                seqs = ["".join(rng.choice("ACGTACGTRY") for _ in range(6)) for _ in range(n)]
                seqs += seqs[: n // 4]
                pamless = set(pam) == {"N"}
                targets = sorted(((rng.randrange(-30, 1 << 20), s,
                                   "".join(rng.choice("ACGTN") for _ in pam),
                                   rng.choice((0, 0x40)) | rng.randrange(0, 3)
                                   | (0x20 if pamless and rng.random() < 0.3 else 0))
                                  for s in seqs), key=lambda t: t[1])
                data = save_tst2(targets, pam, 6, 2, at_start)
                path = os.path.join(tmp, f"{pam}_chr1_1.bin")
                with open(path, "wb") as fh:
                    fh.write(data)
                self.assertIsNone(crispritz_bins.bin_problem(path))
                self.assertEqual(crispritz_bins._leaves(path), len(targets))
                for size in range(len(data)):
                    with open(path, "wb") as fh:
                        fh.write(data[:size])
                    self.assertIsNotNone(crispritz_bins.bin_problem(path), size)

    def test_what_searchtst_refuses_or_buildtst_never_writes_is_told(self):
        tmp = tempfile.mkdtemp(prefix="fmt2_bad_")
        self.addCleanup(shutil.rmtree, tmp, True)
        good = [(3, "ACGTAC", "AGG", 0x02), (9, "ACGTTT", "TGG", 0x41)]
        path = os.path.join(tmp, "NGG_chr1_1.bin")

        def told(data, name=path):
            with open(name, "wb") as fh:
                fh.write(data)
            return crispritz_bins.bin_problem(name)

        self.assertIsNone(told(save_tst2(good, "NGG")))
        # the four-allele flag (bit 2) is one CRISPRitz 2.10.0 writes; bit 3 is none
        self.assertIsNone(told(save_tst2(good, "NGG", q=True)))
        flagged = bytearray(save_tst2(good, "NGG"))
        flagged[12] |= 8
        self.assertIn("malformed index header", told(bytes(flagged)))
        # more pad bases than the max bulges, and a flank off the contig of a dedicated PAM
        self.assertEqual(told(save_tst2([(3, "ACGTAC", "AGG", 0x03)] + good[1:], "NGG")),
                         "corrupt (target 0 of 2)")
        self.assertEqual(told(save_tst2([good[0], (9, "ACGTTT", "TGG", 0x21)], "NGG")),
                         "corrupt (target 1 of 2)")
        self.assertIsNone(told(save_tst2([(3, "ACGTAC", "NNN", 0x22)], "NNN"),
                               os.path.join(tmp, "NNN_chr1_1.bin")))
        # a header PAM other than the file name's
        self.assertIn("is named for 'NRG'", told(save_tst2(good, "NGG"),
                                                 os.path.join(tmp, "NRG_chr1_1.bin")))
        # a format newer than 2
        newer = bytearray(save_tst2(good, "NGG"))
        newer[4:8] = struct.pack("=i", 3)
        self.assertIn("index format 3", told(bytes(newer)))

    def test_the_header_reader_refuses_what_crispritz_2_10_refuses(self):
        """index_fit.read_index_header is CRISPRitz 2.10.0's own (crispritz.py index-info):
        the same refusal, in its words, of a file damaged in one way each, and of a good one
        the same build parameters, the four-allele flag among them
        (test/index_format2/check_real_engine.py compares the two on the real engine)."""
        tmp = tempfile.mkdtemp(prefix="fmt2_reader_")
        self.addCleanup(shutil.rmtree, tmp, True)
        data = save_tst2([(3, "ACGTAC", "AGG", 0x02), (9, "ACGTTT", "TGG", 0x41)], "NGG", 6, 2,
                         q=True)

        def read(blob, name="NGG_chr1_1.bin"):
            with open(os.path.join(tmp, name), "wb") as fh:
                fh.write(blob)
            return index_fit.read_index_header(os.path.join(tmp, name))

        def put(at, value):
            return data[:at] + struct.pack("=i", value) + data[at + 4:]

        self.assertEqual(read(data), {
            "format": 2, "pam_at_start": False, "pamless": False, "four_allele": True,
            "guide_len": 6, "pam_len": 3, "max_bulges": 2, "pam": "NGG", "complete": True,
            "num_leaves": 2})
        for blob, name, why in (
                (data[:20], "NGG_chr1_1.bin", "cut-short index header"),
                (put(4, 3), "NGG_chr1_1.bin", "index format 3, newer than this CRISPRitz reads"),
                (put(12, 12), "NGG_chr1_1.bin", "malformed index header (version 2, flags 12,"),
                (put(8, 32), "NGG_chr1_1.bin", "malformed index header (version 2, flags 4, "
                                               "PAM length 3, header length 32,"),
                (put(35, 7), "NGG_chr1_1.bin", "malformed index: its header gives guide length "
                                               "6, its body 7"),
                (data, "NRG_chr1_1.bin", "holds the PAM/flank region 'NGG' but is named for 'NRG'"),
                (put(31, 0), "NGG_chr1_1.bin", "not a valid CRISPRitz index (.bin) file "
                                               "(numLeaves=0)"),
                (put(31, 1 << 24), "NGG_chr1_1.bin", "damaged: its 16777216 leaves need at least "),
                (data[:-4], "NGG_chr1_1.bin", "cut short: it does not end in the end-of-index "
                                              "marker")):
            with self.subTest(why=why):
                self.assertTrue(read(blob, name).get("error", "").startswith(why),
                                read(blob, name))
        # format 1: a leaf count the file cannot hold is damage too
        self.assertTrue(read(struct.pack("=ii", 7, 20) + bytes(8)).get("error", "").startswith(
            "damaged: its 7 leaves need at least 61 bytes, the file holds 16"))

    def test_a_format_1_tree_is_read_as_before(self):
        """(the format-1 checks are TheIndexCheckTellsAWholeTreeFromOneCutShort's, in
        test_build_index_only_atomic: here, that _leaves still reads its first int)"""
        tmp = tempfile.mkdtemp(prefix="fmt1_bins_")
        self.addCleanup(shutil.rmtree, tmp, True)
        path = one_target_bin(os.path.join(tmp, "NGG_chr1_1.bin"), "NGG", fmt=1, leaves=7)
        self.assertEqual(crispritz_bins._leaves(path), 7)


class TheSiteRuleIsThatOfTheIndexsFormat(unittest.TestCase):
    """holds_target and count_targets by the rule of format 2 agree with buildTST 2.10.0's
    loops (buildtst2_targets) on random contigs of every kind (N runs, IUPAC codes, Q, sites
    at either end and beside N runs) and the shipped PAMs, 3' and 5', PAM-less among them; and
    a contig whose only site is one format 1 cannot hold is told apart."""

    PAMS = ("NNNNNNNNNNNNNNNNNNNNNRG 3", "NNNNNNNNNNNNNNNNNNNNNGG 3", "NNNNNNNNNNNNNNNNNNNNNNN 3",
            "TTTVNNNNNNNNNNNNNNNNNNNNNNN -4", "NNNNNNNNNNNNNNNNNNNNNNNGRRT 6",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNN -4", "NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8")

    def test_it_agrees_with_buildtst_2_10_on_random_contigs(self):
        r = random.Random(13)
        self.addCleanup(setattr, crispritz_bins, "_CHUNK", crispritz_bins._CHUNK)
        held = 0
        for case in range(3000):
            motif, limit = r.choice(self.PAMS).split()
            limit = int(limit)
            pam = motif[:-limit] if limit < 0 else motif[-limit:]
            bulges = r.randint(0, 3)
            n = r.choice((r.randint(0, 40), r.randint(20, 60), r.randint(40, 300)))
            mix = r.choice(("ACGT", "ACGT" * 6 + "N", "ACGTN", "ACGTNNNNRYSWKMBDHV",
                            "ACGTNRYSWKMBDHVQ", "ACGTQ"))
            seq = "".join(r.choice(mix) for _ in range(n))
            if n > 10 and r.random() < 0.3:
                a = r.randrange(n)
                seq = seq[:a] + "N" * min(n - a, r.randint(1, 30)) + seq[a + 30:]
            crispritz_bins._CHUNK = r.choice((1, 3, 1 << 20))
            want = buildtst2_targets(seq, pam, len(motif), limit < 0)
            held += want > 0
            self.assertEqual(crispritz_bins.count_targets(seq, pam, len(motif), limit < 0,
                                                          bulges, 2), want,
                             (seq, motif, limit))
            self.assertEqual(crispritz_bins.holds_target(seq, pam, len(motif), limit < 0, bulges,
                                                         2), want > 0, (seq, motif, limit))
        self.assertGreater(held, 600)
        self.assertLess(held, 2600)

    def test_a_site_only_format_2_holds(self):
        """INDEX_FORMAT section 8's two contigs one target long: a + strand NGG site whose
        PAM ends on the last base and whose pad runs off the start, and a - strand one whose
        pad runs off the end. Format 1 holds neither, format 2 both; the bulges change
        nothing in format 2."""
        guide = "CTAACAGTTGCTTTTATCAC"
        for seq in (guide + "TGG", "CCA" + guide):
            for bulges in (0, 1, 2):
                self.assertFalse(crispritz_bins.holds_target(seq, "NGG", 23, False, max(bulges, 1)))
                self.assertTrue(crispritz_bins.holds_target(seq, "NGG", 23, False, bulges, 2))
                self.assertEqual(crispritz_bins.count_targets(seq, "NGG", 23, False, bulges, 2), 1)
        # next to an N run, and a guide holding an N, which no format keeps
        self.assertTrue(crispritz_bins.holds_target("N" * 5 + guide + "TGG", "NGG", 23, False, 2, 2))
        self.assertFalse(crispritz_bins.holds_target("N" + guide[1:] + "TGG", "NGG", 23, False, 0, 2))
        # PAM-less: a window is kept when the flank base next to its guide is on the contig
        # (one more base, a window on each strand), the rest of the flank may be off it
        self.assertEqual(crispritz_bins.count_targets(guide, "NNN", 23, False, 2, 2), 0)
        self.assertEqual(crispritz_bins.count_targets(guide + "A", "NNN", 23, False, 2, 2), 2)
        self.assertEqual(crispritz_bins.count_targets("A" + guide + "A", "NNN", 23, False, 2, 2), 4)

    def test_q_is_any_base(self):
        guide = "CTAACAGTTGCTTTTATCAC"
        self.assertTrue(crispritz_bins.holds_target("Q" + guide[1:] + "TQG", "NGG", 23, False, 0, 2))
        self.assertEqual(crispritz_bins.count_targets("QQQ", "NGG", 3, False, 0, 2), None)


class MissingBinsReadsTheIndexsOwnFormat(_Engine):
    """missing_bins applies the site rule of the format the index was built in: from its own
    .bin files, whatever CRISPRitz is installed; for a folder with none, from its record,
    else from what its caller names (the format of the engine that has just built it, or of
    the index's other sub-indexes), else format 2's."""

    GUIDE = "CTAACAGTTGCTTTTATCAC"

    def genome(self, contigs):
        folder = os.path.join(self.tmp, "g")
        os.makedirs(folder, exist_ok=True)
        for name, seq in contigs.items():
            with open(os.path.join(folder, f"{name}.fa"), "w") as fh:
                fh.write(f">{name}\n{seq}\n")
        pam = os.path.join(self.tmp, "pam.txt")
        with open(pam, "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNGG 3\n")
        return folder, pam

    def test_a_lost_bin_of_a_site_only_format_2_holds_is_found(self):
        """chrE is one target long (a site at its end): a format-2 build holds it, a
        format-1 build has no site there. With its .bin lost from a format-2 index, the
        format-1 rule (2.8.1's) finds nothing missing; the index's own rule does, whichever
        CRISPRitz is installed."""
        genome, pam = self.genome({"chrA": "ACGT" * 10 + self.GUIDE + "AGG" + "ACGT" * 10,
                                   "chrE": self.GUIDE + "TGG"})
        index = os.path.join(self.tmp, "NGG_2_g")
        one_target_bin(os.path.join(index, "NGG_chrA_1.bin"), "NGG")
        for engine in (2, None, False):
            with self.subTest(engine=engine):
                self.engine(engine)
                self.assertEqual(crispritz_bins.missing_bins(index, genome, pam, 2),
                                 ([("chrE.fa", "chrE", 0, None)], []))
        # a format-1 index of the same genome, whole (chrE has no site in format 1)
        old = os.path.join(self.tmp, "NGG_2_old")
        one_target_bin(os.path.join(old, "NGG_chrA_1.bin"), "NGG", fmt=1)
        for engine in (2, None):
            with self.subTest(engine=engine, index="format 1"):
                self.engine(engine)
                self.assertEqual(crispritz_bins.missing_bins(old, genome, pam, 2), ([], []))

    def test_a_folder_with_no_bin_takes_its_record_then_its_callers_format(self):
        genome, pam = self.genome({"chrE": self.GUIDE + "TGG"})
        empty = os.path.join(self.tmp, "NGG_2_fakechrE")
        os.makedirs(empty)
        lost = ([("chrE.fa", "chrE", 0, None)], [])
        self.assertEqual(crispritz_bins.missing_bins(empty, genome, pam, 2), lost)  # format 2's
        self.assertEqual(crispritz_bins.missing_bins(empty, genome, pam, 2, 1), ([], []))
        self.assertEqual(crispritz_bins.missing_bins(empty, genome, pam, 2, 2), lost)
        for recorded, want in ((1, ([], [])), (2, lost)):
            with open(os.path.join(empty, index_fit.INDEX_META), "w") as fh:
                json.dump({"index_format": recorded}, fh)
            # (the record wins over the caller's format)
            self.assertEqual(crispritz_bins.missing_bins(empty, genome, pam, 2, 3 - recorded), want)

    def test_a_full_last_part_of_format_2_is_counted(self):
        """_leaves reads numLeaves after a format-2 header, so a contig whose last part is
        full is counted (here a tree of buildTST holds 1 target)."""
        genome, pam = self.genome({"chrA": self.GUIDE + "TGG" + "A" * 5 + "CCA" + self.GUIDE})
        index = os.path.join(self.tmp, "NGG_2_g")
        one_target_bin(os.path.join(index, "NGG_chrA_1.bin"), "NGG")
        self.addCleanup(setattr, crispritz_bins, "TARGETS_A_TREE", crispritz_bins.TARGETS_A_TREE)
        crispritz_bins.TARGETS_A_TREE = 1
        self.assertEqual(crispritz_bins.missing_bins(index, genome, pam, 2),
                         ([("chrA.fa", "chrA", 1, 2)], []))


def _functions(*names):
    """Those functions of the shipped search script, verbatim, as top-level functions."""
    with open(os.path.join(HERE, "submit_job_automated_new_multiple_vcfs.sh")) as fh:
        text = fh.read()
    out = []
    for name in names:
        start = text.index(f"\t{name}() {{")
        end = text.index("\n\t}\n", start) + len("\n\t}\n")
        out.append("\n".join(line[1:] if line.startswith("\t") else line
                             for line in text[start:end].splitlines()))
    return "\n".join(out)


def _shell():
    """The search script's functions these tests run."""
    return _functions("_engine_index_format", "_unreadable_bin", "_how_built_again",
                      "_check_index_format", "_index_guide_len", "_index_side", "_header_side",
                      "_index_fits_guide", "_resolve_index_dir", "_say_other_side",
                      "_index_name_parts", "_pam_holds", "_pam_narrower")


class _Index(_Engine):
    def index(self, name, fmt=2, guide=20, at_start=False, contigs=("chr1",), mb=None, lib=None):
        """An index folder ``name`` (<PAM>_<N>_<genome>) in library ``lib`` (tmp/lib), each
        contig a whole .bin of one target, in index format ``fmt``."""
        pam, n = name.split("_")[:2]
        folder = os.path.join(lib or os.path.join(self.tmp, "lib"), name)
        for c in contigs:
            one_target_bin(os.path.join(folder, f"{pam}_{c}_1.bin"), pam, guide,
                           int(n) if mb is None else mb, at_start, fmt)
        return folder

    def shell(self, call, folder="", pam_line="NNNNNNNNNNNNNNNNNNNNNGG 3", bulges=(1, 1),
              prune=0):
        """``call`` run after the search script's functions, with its variables for a PAM
        file of ``pam_line`` and these bulges, on the PATH the test set."""
        seq, pos = pam_line.split()
        script = (f"pos={pos}\nfullseqpam={seq}\nbDNA={bulges[0]}\nbRNA={bulges[1]}\n"
                  f"prune_pam={prune}\n{_shell()}\n{call}\n")
        return subprocess.run(["bash", "-c", script, "bash", folder], capture_output=True,
                              text=True, env=dict(os.environ))


class TheInstalledEngineSaysWhatItReads(_Engine):
    """index_fit reads the index format the installed CRISPRitz reads from its crispritz.py
    (INDEX_FORMAT, which CRISPRitz states from 2.10.0; 2.9.1 states none and reads format 1
    only), and the search script reads the same. Each case is told on this host's PATH, and
    on one whose only folder holds an installed CRISPRitz 2.9.1 beside the bash and sed the
    search script runs, as the Docker image's /opt/conda/bin holds its crispritz.py beside
    python. The case of none takes the installed one off PATH and leaves bash and sed: it
    was told on a PATH that still held it, and failed wherever CRISPRitz is installed (the
    release image, a source or a Bioconda install)."""

    def test_from_its_crispritz_py(self):
        installed = stub_engine(os.path.join(self.tmp, "installed"), "2.9.1", None)
        for tool in ("bash", "sed"):
            os.symlink(shutil.which(tool), os.path.join(installed, tool))
        for path in (os.environ.get("PATH", ""), installed):
            for stated, version, fmt, readable in ((2, "2.10.0", 2, (1, 2)),
                                                   (None, "2.9.1", 1, (1,)),
                                                   (False, None, None, (1, 2))):
                with self.subTest(path="installed" if path == installed else "host",
                                  stated=stated):
                    os.environ["PATH"] = path
                    self.engine(stated)
                    self.assertEqual(index_fit.engine_index_format(), fmt)
                    self.assertEqual(index_fit.engine_version(), version)
                    self.assertEqual(index_fit.readable_formats(), readable)
                    out = subprocess.run(["bash", "-c", f"{_shell()}\n_engine_index_format"],
                                         capture_output=True, text=True, env=dict(os.environ))
                    self.assertEqual((out.stdout, out.stderr), (f"{fmt or 2}\n", ""))


class AFormat2IndexIsReadNotRefused(_Index):
    """A format-2 index is read: its guide length and its side from its header. It is
    refused only where the installed CRISPRitz cannot read it, naming both; one of a format
    newer than this CRISPRme reads, always. The search script says the same."""

    def test_its_guide_length_and_side_come_from_its_header(self):
        self.engine(2)
        for guide, at_start in ((20, False), (21, False), (23, True)):
            with self.subTest(guide=guide, at_start=at_start):
                folder = self.index(f"{'TTTV' if at_start else 'NGG'}_2_g{guide}", guide=guide,
                                    at_start=at_start)
                self.assertEqual(index_fit.built_guide_len(folder), guide)
                self.assertEqual(index_fit.built_five_prime(folder), at_start)
                out = self.shell('_index_guide_len "$1"; echo; _index_side "$1"; echo "$?"', folder)
                self.assertEqual((out.stdout, out.stderr),
                                 (f"{guide}\n{5 if at_start else 3}\n0\n", ""))
        old = self.index("NGG_2_old", fmt=1, guide=21)
        self.assertEqual(index_fit.built_guide_len(old), 21)
        self.assertIsNone(index_fit.built_five_prime(old))  # (format 1 records no side)
        out = self.shell('_index_guide_len "$1"; echo; _index_side "$1"; echo "$?"', old)
        self.assertEqual(out.stdout, "21\n3\n1\n")

    def test_it_is_refused_only_where_the_installed_crispritz_cannot_read_it(self):
        folder = self.index("NGG_2_g")
        self.engine(2)
        self.assertIsNone(index_fit.unreadable_format_error(folder))
        out = self.shell('_check_index_format "$1"; echo went on', folder)
        self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "went on\n", ""))
        self.engine(None)  # CRISPRitz 2.9.1
        msg = index_fit.unreadable_format_error(folder)
        self.assertEqual(msg, (
            f"{folder} was built by a newer CRISPRitz (NGG_chr1_1.bin is in index format 2, "
            f"which CRISPRitz 2.10.0 and later write) than the one installed here, CRISPRitz "
            f"2.9.1, which reads index format 1 only. Install the CRISPRitz this CRISPRme "
            f"installs with (2.10.0 or later; the Docker image, install_from_source.sh and the "
            f"Bioconda package carry it), or build the index again with the CRISPRitz installed "
            f"here, with {folder} moved out of the way first"))
        out = self.shell('_check_index_format "$1"; echo went on', folder)
        self.assertEqual((out.returncode, out.stdout, out.stderr), (1, "", f"ERROR: {msg}\n"))
        # an INDELS index, one folder a chromosome
        indels = os.path.join(self.tmp, "lib", "NGG_2_g+ds_INDELS")
        one_target_bin(os.path.join(indels, "NGG_2_fakechr1", "NGG_fakechr1_1.bin"), "NGG", fmt=1)
        one_target_bin(os.path.join(indels, "NGG_2_fakechr2", "NGG_fakechr2_1.bin"), "NGG")
        self.assertIn("NGG_2_fakechr2/NGG_fakechr2_1.bin is in index format 2",
                      index_fit.unreadable_format_error(indels))
        out = self.shell('_check_index_format "$1"', indels)
        self.assertEqual(out.stderr.count("ERROR: "), 1)
        self.assertIn("NGG_2_fakechr2/NGG_fakechr2_1.bin is in index format 2", out.stderr)

    def test_a_format_newer_than_this_crisprme_reads_is_always_refused(self):
        folder = self.index("NGG_2_g")
        path = os.path.join(folder, "NGG_chr1_1.bin")
        with open(path, "r+b") as fh:
            fh.seek(4)
            fh.write(struct.pack("=i", 3))
        for engine in (2, False):
            with self.subTest(engine=engine):
                self.engine(engine)
                msg = index_fit.unreadable_format_error(folder)
                self.assertEqual(msg, (
                    f"{folder} was built by a newer CRISPRitz (NGG_chr1_1.bin is in index format "
                    f"3): this CRISPRme reads index formats 1 and 2 (which CRISPRitz writes up to "
                    f"2.9.1 and from 2.10.0), and needs a newer release to search it. Use a "
                    f"CRISPRme that reads index format 3, or build the index again with the "
                    f"CRISPRitz installed here, with {folder} moved out of the way first"))
                out = self.shell('_check_index_format "$1"', folder)
                self.assertEqual((out.returncode, out.stderr), (1, f"ERROR: {msg}\n"))
                self.assertIsNone(index_fit.built_guide_len(folder))

    def test_a_format_1_index_is_searched_with_a_note_once_the_engine_writes_format_2(self):
        # (of the working directory's genome_library; one of the --index-path library:
        # TheNoteMovesOnlyAnIndexThisInstallBuildsAgain)
        work = os.path.join(self.tmp, "work")
        old = self.index("NGG_2_g", fmt=1, lib=os.path.join(work, "genome_library"))
        self.engine(None)
        call = f'current_working_directory="{work}"\n_check_index_format "$1"; echo went on'
        out = self.shell(call, old)
        self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "went on\n", ""))
        self.engine(2)
        out = self.shell(call, old)
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        # (build-index-only keeps an index it finds in place, and download installs the
        # archives published until the rebuild, of format 1: AFormat1IndexIsBuiltAgain)
        self.assertEqual(out.stdout, (
            "NOTE: NGG_2_g is in index format 1 (built by CRISPRitz 2.9.1 or earlier): it lacks "
            "the sites whose guide lies within a few bases of a contig end or an N run, which an "
            "index built by CRISPRitz 2.10.0 or later holds (METHODS section 9). To search them, "
            f"move {old} out of the way and build it again with crisprme.py build-index-only "
            "(which keeps an index it finds in place).\n"
            "went on\n"))


def _replaced(folder, source=None, without=None):
    """The end of the note of index ``folder``, of format 1, of an --index-path library: replace
    it with one built again (from ``source``); with ``without``, or run the search without
    --index-path (under --prune-pam too, which needs no --index-path since 2.10.0)."""
    return ("A search builds no index in the --index-path library: to search them, replace "
            f"{folder} with one built again by crisprme.py build-index-only"
            + (f" from {source}" if source else "")
            + " (which builds into the genome_library folder of its --path, and keeps an index it "
            "finds there in place)"
            + (", or run this search without --index-path" if without else "") + ".")


class TheNoteMovesOnlyAnIndexThisInstallBuildsAgain(_Index):
    """The note of an index of format 1 said, of a variant and an INDELS index as of a
    reference one, to move it out of the way and build it again with build-index-only. On an
    install of downloaded indexes, which holds no VCFs and no enriched genome (the README's
    Quickstart), build-index-only refuses the dataset ("VCF folder ... not found") and the
    index stays moved out of the way; with the VCFs alone it enriches the whole dataset again
    first. The note now sends an index out of the way only where this install builds it again
    from what it holds -- the genome folder, or the enriched genomes in place and the dataset's
    VCFs -- and says elsewhere what building it again takes, and that the dataset's archives
    rebuilt in format 2 hold those sites.

    And only an index of the working directory's genome_library, where build-index-only builds
    it again. A search takes one of the --index-path library first -- its reference index from
    there alone, its INDELS index from beside its variant index -- and of an index there the
    note said the same: built again as it said, it was not searched, and the next search
    stopped ("no matching reference index under --index-path", "INDELS index ... is not beside
    its variant index"). Of an index there the note says to replace it with one built again,
    or, where the search then builds every index it needs, to run it without --index-path: a
    reference search, or a variant one with the dataset's VCFs, enriched genome and fake-indel
    genome in place (AFormat1IndexIsBuiltAgain follows each through the CLI). Under
    --prune-pam it said to run it without --prune-pam too, "which needs it": since 2.10.0
    --prune-pam needs no --index-path, and the search keeps it."""

    PREFIX = ("NOTE: {} is in index format 1 (built by CRISPRitz 2.9.1 or earlier): it lacks the "
              "sites whose guide lies within a few bases of a contig end or an N run, which an "
              "index built by CRISPRitz 2.10.0 or later holds (METHODS section 9). ")
    ARCHIVES = ("Where ds is published, its archives rebuilt in index format 2 hold them: "
                "building it again with crisprme.py build-index-only takes the dataset's ")

    def setUp(self):
        super().setUp()
        self.engine(2)
        self.work = os.path.join(self.tmp, "work")
        # the three in the working directory's genome_library, and in an --index-path library
        # outside it, of another name (tmp/lib)
        for lib, names in ((os.path.join(self.work, "genome_library"), ("ref", "var", "indels")),
                           (None, ("lref", "lvar", "lindels"))):
            ref = self.index("NGG_2_hg38", fmt=1, lib=lib)
            var = self.index("NGG_2_hg38+ds", fmt=1, lib=lib)
            indels = var + "_INDELS"
            one_target_bin(os.path.join(indels, "NGG_2_fakechr1", "NGG_fakechr1_1.bin"), "NGG",
                           fmt=1)
            for name, folder in zip(names, (ref, var, indels)):
                setattr(self, name, folder)

    def install(self, *parts):
        """The working directory of the search holds these of the dataset: "vcfs" (its VCFs),
        "enriched" (its enriched genome), "unfinished" (that genome, marked unfinished by a
        build-index-only run), "fake" (its fake-indel genome)."""
        for top in ("VCFs", "Genomes"):
            shutil.rmtree(os.path.join(self.work, top), ignore_errors=True)
        os.makedirs(os.path.join(self.work, "VCFs"))
        if "vcfs" in parts:
            os.makedirs(os.path.join(self.work, "VCFs", "ds"))
            open(os.path.join(self.work, "VCFs", "ds", "ds.chr1.vcf.gz"), "w").close()
        for part, folder in (("enriched", "hg38+ds"), ("unfinished", "hg38+ds/.building"),
                             ("fake", "hg38+ds_INDELS")):
            if part in parts:
                os.makedirs(os.path.join(self.work, "Genomes", folder), exist_ok=True)

    def note(self, folder, vcf="ds", prune=0):
        """The search script's note of index ``folder``, run as the search runs it, in the
        working directory install() made, a search of dataset ``vcf`` ("_": a reference search),
        under --prune-pam with ``prune``: on stdout, and nothing on stderr."""
        out = self.shell(f'ref_name=hg38\nvcf_name={vcf}\ncurrent_working_directory="{self.work}"\n'
                         f'vcf_folder="{self.work}/VCFs/{vcf}"\n_check_index_format "$1"', folder,
                         prune=prune)
        self.assertEqual((out.returncode, out.stderr), (0, ""))
        prefix = self.PREFIX.format(os.path.basename(folder))
        self.assertTrue(out.stdout.startswith(prefix) and out.stdout.endswith("\n"), out.stdout)
        return out.stdout[len(prefix):-1]

    def moved(self, folder, source=None):
        return (f"To search them, move {folder} out of the way and build it again with crisprme.py "
                f"build-index-only (which keeps an index it finds in place)"
                + (f", from {source}." if source else "."))

    def test_of_a_download_only_the_reference_index(self):
        # (the README's Quickstart: download --what index, VCFs/ empty, no enriched genome)
        self.install()
        self.assertEqual(self.note(self.ref), self.moved(self.ref))
        for folder in (self.var, self.indels, self.lvar, self.lindels):
            self.assertEqual(self.note(folder),
                             self.ARCHIVES + "VCFs (VCFs/ds), which this install does not hold.")
        # (a search without --index-path would find no variant index, and build none)
        self.assertEqual(self.note(self.lref), _replaced(self.lref))

    def test_with_the_vcfs_and_no_enriched_genome_the_whole_dataset_is_enriched_again(self):
        for parts in (("vcfs",), ("vcfs", "fake"), ("vcfs", "enriched", "unfinished", "fake")):
            with self.subTest(parts=parts):
                self.install(*parts)
                for folder in (self.var, self.indels, self.lvar, self.lindels):
                    self.assertEqual(self.note(folder), self.ARCHIVES + (
                        "enriched genome (Genomes/hg38+ds), which this install does not hold, or "
                        "else its VCFs (VCFs/ds), from which the whole dataset is enriched again "
                        "first."))
                self.assertEqual(self.note(self.ref), self.moved(self.ref))
                self.assertEqual(self.note(self.lref), _replaced(self.lref))

    def test_from_the_enriched_genomes_in_place_with_the_vcfs(self):
        self.install("vcfs", "enriched", "fake")
        self.assertEqual(self.note(self.var), self.moved(self.var, "Genomes/hg38+ds"))
        self.assertEqual(self.note(self.indels), self.moved(self.indels, "Genomes/hg38+ds_INDELS"))
        # of the --index-path library: replaced there, or the search run without it, which then
        # builds the three here
        self.assertEqual(self.note(self.lref), _replaced(self.lref, without=True))
        self.assertEqual(self.note(self.lvar),
                         _replaced(self.lvar, "Genomes/hg38+ds", without=True))
        self.assertEqual(self.note(self.lindels),
                         _replaced(self.lindels, "Genomes/hg38+ds_INDELS", without=True))
        # the INDELS index is built from the fake-indel genome, which build-index-only makes
        # only as it enriches the dataset: with the enriched genome in place, it does not (and a
        # search without --index-path builds no INDELS index then)
        self.install("vcfs", "enriched")
        self.assertEqual(self.note(self.var), self.moved(self.var, "Genomes/hg38+ds"))
        self.assertEqual(self.note(self.lvar), _replaced(self.lvar, "Genomes/hg38+ds"))
        self.assertEqual(self.note(self.lref), _replaced(self.lref))
        for folder in (self.indels, self.lindels):
            self.assertEqual(self.note(folder), self.ARCHIVES + (
                "fake-indel genome (Genomes/hg38+ds_INDELS), which this install does not hold."))
        # and build-index-only takes the dataset's VCFs, whatever is in place
        self.install("enriched", "fake")
        for folder in (self.var, self.indels, self.lvar, self.lindels):
            self.assertEqual(self.note(folder),
                             self.ARCHIVES + "VCFs (VCFs/ds), which this install does not hold.")

    def test_of_a_reference_search(self):
        # a reference search without --index-path builds its index in genome_library, under
        # --prune-pam too, which complete-search refused without --index-path before 2.10.0
        # (the note said to drop it as well)
        self.install()
        for prune in (0, 1):
            with self.subTest(prune=prune):
                self.assertEqual(self.note(self.ref, "_", prune), self.moved(self.ref))
                self.assertEqual(self.note(self.lref, "_", prune),
                                 _replaced(self.lref, without=True))


def _load_crisprme():
    """crisprme.py's definitions, as test_bulge_cap loads them."""
    import test_bulge_cap
    return test_bulge_cap.crisprme


class TheScanPrefersAnIndexTheEngineReads(_Index):
    """PKG1-scan-takes-format2: of a library holding a format-2 NGG_2_toy and a format-1
    NGG_3_toy, a 1/1 search took NGG_2_toy (the smaller N), which the installed CRISPRitz 2.9.1
    cannot read, and stopped on it, where NGG_3_toy served. The search script's scan and
    complete-search's (crisprme._index_the_search_takes) now take an index the installed
    CRISPRitz reads first, and the one they would stop on only where no other serves."""

    def test_the_search_script_and_complete_search_agree(self):
        new = self.index("NGG_2_toy")
        old = self.index("NGG_3_toy", fmt=1)
        lib = os.path.dirname(new)
        crisprme = _load_crisprme()
        for engine, want in ((None, old), (2, new)):
            with self.subTest(engine=engine):
                self.engine(engine)
                out = self.shell(f'_resolve_index_dir "$1" NGG toy 2', lib)
                self.assertEqual((out.stdout, out.stderr), (want + "\n", ""))
                self.assertEqual(crisprme._index_the_search_takes(
                    [lib], "NGG", "toy", 2, False, True, 20, 1, 1, False), want)
        shutil.rmtree(old)
        self.engine(None)  # (the only one: taken, and the search stops on it)
        self.assertEqual(self.shell(f'_resolve_index_dir "$1" NGG toy 2', lib).stdout, new + "\n")
        self.assertEqual(crisprme._index_the_search_takes(
            [lib], "NGG", "toy", 2, False, True, 20, 1, 1, False), new)


class EachIndexRecordsWhatItHolds(_Index):
    """An index records what it holds in .index_meta.json -- its index format, build PAM,
    side, flank and guide length, max bulges -- each read from its own headers where they
    hold it, else, for an index this run has just built, from its PAM file and the installed
    CRISPRitz; an index with no .bin records its format nowhere else."""

    def pam_file(self, line):
        path = os.path.join(self.tmp, "pam.txt")
        with open(path, "w") as fh:
            fh.write(line + "\n")
        return path

    def test_a_format_2_index_records_its_header(self):
        self.engine(2)
        folder = self.index("TTTV_3_g", guide=23, at_start=True, contigs=("chr1", "chr2"))
        # (a PAM file of another geometry does not overrule the header)
        self.assertIsNone(index_fit.write_index_meta(
            folder, self.pam_file("NNNNNNNNNNNNNNNNNNNNNNGG 3"), built=True))
        self.assertEqual(os.listdir(folder).count(".index_meta.json"), 1)
        self.assertEqual(index_fit.read_index_meta(folder), {
            "index_format": 2, "build_pam": "TTTV", "side": "5'", "flank_len": 4,
            "guide_len": 23, "max_bulges": 3, "pamless": False, "four_allele": False,
            "built_by": {"crispritz": "2.10.0", "crisprme": crisprme_version()}})

    def test_a_format_1_index_found_in_place_records_what_it_holds(self):
        self.engine(2)
        folder = self.index("NGG_2_g", fmt=1, guide=21)
        self.assertIsNone(index_fit.write_index_meta(folder, self.pam_file("N" * 20 + "NGG 3")))
        self.assertEqual(index_fit.read_index_meta(folder), {
            "index_format": 1, "build_pam": "NGG", "flank_len": 3, "guide_len": 21,
            "max_bulges": 2, "pamless": False, "four_allele": False})
        # its side only where its .pam_build records it (format 1 stores none)
        other = self.index("NNNN_2_h", fmt=1)
        with open(os.path.join(other, ".pam_build"), "w") as fh:
            fh.write("20bp-NNNN-x 24 -4\n")
        index_fit.write_index_meta(other)
        self.assertEqual(index_fit.read_index_meta(other)["side"], "5'")
        self.assertTrue(index_fit.read_index_meta(other)["pamless"])

    def test_one_with_no_bin_records_the_format_of_the_engine_that_built_it(self):
        empty = os.path.join(self.tmp, "lib", "NGG_2_fakechr9")
        os.makedirs(empty)
        pam = self.pam_file("N" * 20 + "NGG 3")
        self.assertIsNone(index_fit.write_index_meta(empty, pam))  # (found so: nothing to say)
        self.assertIsNone(index_fit.read_index_meta(empty))
        for engine, fmt in ((2, 2), (None, 1)):
            with self.subTest(engine=engine):
                self.engine(engine)
                self.assertIsNone(index_fit.write_index_meta(empty, pam, built=True))
                meta = index_fit.read_index_meta(empty)
                self.assertEqual((meta["index_format"], meta["build_pam"], meta["side"],
                                  meta["guide_len"], meta["max_bulges"]), (fmt, "NGG", "3'", 20, 2))
                self.assertEqual(index_fit.recorded_format(empty), fmt)
                self.assertEqual(crispritz_bins.index_format(empty), fmt)

    def test_a_record_is_kept_by_a_run_that_did_not_build_the_index(self):
        self.engine(2)
        folder = self.index("NGG_2_g")
        index_fit.write_index_meta(folder, built=True)
        with open(os.path.join(folder, index_fit.INDEX_META)) as fh:
            before = fh.read()
        index_fit.write_index_meta(folder)
        with open(os.path.join(folder, index_fit.INDEX_META)) as fh:
            self.assertEqual(fh.read(), before)
        self.assertIn("built_by", json.loads(before))

    def test_the_record_is_a_dot_file_and_the_search_script_writes_it(self):
        """crispritz.py search refuses an index folder that holds any other file but a .bin,
        unless its name starts with a dot (CRISPRitz 2.9.1 and 2.10.0: an ERROR on stdout,
        no targets file, exit 0), so the record is a dot file, as .pam_build is. The search
        script records an index it builds (_record_index_side) through index_fit's command."""
        self.assertTrue(index_fit.INDEX_META.startswith("."))
        self.engine(2)
        folder = self.index("NGG_2_g")
        out = subprocess.run([sys.executable, os.path.join(HERE, "index_fit.py"), "meta", folder,
                              self.pam_file("N" * 20 + "NGG 3")], capture_output=True, text=True)
        self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "", ""))
        self.assertEqual(index_fit.read_index_meta(folder)["index_format"], 2)
        self.assertIn('python "$starting_dir/index_fit.py" meta "$1" "$pam_file"',
                      _functions("_record_index_side"))


class TheSideOfAFormat2IndexIsItsHeaders(_Index):
    """An index of format 2 records the side of the guide it was built for in the header of
    its .bin files, and every reader takes it there before .pam_build. build-index-only read
    .pam_build alone: over a 5' PAM-less index with none (crispritz.py built it, or it was
    downloaded) it recorded the 3' PAM file it was given; the web read that record alone, a 3'
    geometry; and a 3' search refused the index as a 5' build "as its .pam_build records",
    which said 3'."""

    PAM_3P = os.path.join(os.path.dirname(HERE), "PAMs", "20bp-NNN-NO-PAM.txt")

    def five_prime(self, record=None, fmt=2):
        """A 5' PAM-less index for a 20-nt guide, in index format ``fmt``, with ``record`` for
        its .pam_build (none for None)."""
        folder = self.index(f"NNN_2_g{fmt}", fmt=fmt, at_start=True)
        if record:
            with open(os.path.join(folder, ".pam_build"), "w") as fh:
                fh.write(record + "\n")
        return folder

    def test_build_index_only_records_no_pam_file_of_the_other_side(self):
        crisprme = _load_crisprme()
        folder = self.five_prime()
        crisprme._write_pam_build(folder, self.PAM_3P)
        self.assertFalse(os.path.exists(os.path.join(folder, ".pam_build")))
        pam_5p = os.path.join(self.tmp, "23bp-NNN-5p.txt")
        with open(pam_5p, "w") as fh:
            fh.write("N" * 23 + " -3\n")
        crisprme._write_pam_build(folder, pam_5p)
        with open(os.path.join(folder, ".pam_build")) as fh:
            self.assertEqual(fh.read(), "23bp-NNN-5p 23 -3\n")

    def test_the_web_reads_the_geometry_in_the_header(self):
        import test_bulge_cap
        index_build_pam = test_bulge_cap.pages_utils.index_build_pam
        folder = self.five_prime()
        self.assertEqual(index_build_pam(folder), (None, 23, -3))
        # (a record of another geometry names a PAM file that does not describe the index)
        for record, name in (("20bp-NNN-NO-PAM 23 3", None), ("23bp-NNN-5p 23 -3", "23bp-NNN-5p")):
            with self.subTest(record=record):
                with open(os.path.join(folder, ".pam_build"), "w") as fh:
                    fh.write(record + "\n")
                self.assertEqual(index_build_pam(folder), (name, 23, -3))
        # an index of format 1 records no side in its headers: its .pam_build is read
        self.assertEqual(index_build_pam(self.five_prime("23bp-NNN-5p 23 -3", fmt=1)),
                         ("23bp-NNN-5p", 23, -3))

    def test_the_web_serves_a_search_of_the_side_in_the_header(self):
        """The web's rule of the side an index serves (pages_utils._index_side_serves, as
        crisprme.py's _side_serves and the search script's _index_side) read a name that says
        the side (2.10.0), then .pam_build, and not the header: a 5' PAM-less index of format 2
        with no .pam_build, of the published width 3, was offered to a 3' search, as a 3'
        build, and not to a 5' one. The header is read after the name, as the others read it."""
        import test_bulge_cap
        serves = test_bulge_cap.pages_utils._index_side_serves
        crisprme = _load_crisprme()
        folder = self.five_prime()
        for motif, three_prime in (("NGG", True), ("TTV", False)):
            with self.subTest(motif=motif):
                self.assertEqual(serves(folder, "NNN", motif, three_prime), not three_prime)
                self.assertEqual(crisprme._side_serves(folder, "NNN", motif, three_prime),
                                 not three_prime)

    def test_the_search_says_what_records_the_side(self):
        call = 'pam_name=20bp-NNN-NO-PAM.txt\ntrue_pam=NNN\n_say_other_side "$1"'
        for fmt, why in ((2, "the header of its .bin files"), (1, "its .pam_build")):
            with self.subTest(fmt=fmt):
                folder = self.five_prime("23bp-NNN-5p 23 -3" if fmt == 1 else
                                         "20bp-NNN-NO-PAM 23 3", fmt=fmt)
                out = self.shell(call, folder, pam_line="N" * 23 + " 3")
                self.assertEqual((out.stdout, out.stderr), ("", (
                    f"       {os.path.basename(folder)} is a PAM-less index for a PAM on the 5' "
                    f"side of the guide (as {why} records), and 20bp-NNN-NO-PAM.txt puts the PAM "
                    "on the 3' side, so a search of it would find nothing\n")))


def crisprme_version():
    with open(os.path.join(os.path.dirname(HERE), "crisprme.py")) as fh:
        for line in fh:
            if line.startswith("version = "):
                return line.split('"')[1]


class TheCoverageCheckReadsFormat2(_Index):
    """index_coverage, which the search script runs on an index before it searches it and
    build-index-only on one it finds in place, checked nothing on a format-2 index (it read
    no guide length from it): it applies the rule of the index's format now, reads its side
    from its header, and refuses a folder of two formats."""

    GUIDE = "CTAACAGTTGCTTTTATCAC"

    def genome(self):
        folder = os.path.join(self.tmp, "Genomes", "g")
        os.makedirs(folder)
        for name, seq in (("chrA", "ACGT" * 10 + self.GUIDE + "AGG" + "ACGT" * 10),
                          ("chrE", self.GUIDE + "TGG")):
            with open(os.path.join(folder, f"{name}.fa"), "w") as fh:
                fh.write(f">{name}\n{seq}\n")
        return folder

    def test_a_lost_tree_of_a_contig_only_format_2_holds_a_site_of(self):
        genome = self.genome()
        new = self.index("NGG_2_g", contigs=("chrA",))
        self.assertEqual(index_coverage.problems(new, genome), [("chrE.fa", "no tree")])
        old = self.index("NGG_2_h", fmt=1, contigs=("chrA",))  # (whole: chrE has no site)
        self.assertEqual(index_coverage.problems(old, genome), [])

    def test_its_side_is_its_headers(self):
        """A 5' index of TTTV that records no side in a .pam_build: chrF holds a site of it,
        read on the 5' side its header says, where a search's 3' PAM file would read none."""
        genome = self.genome()
        with open(os.path.join(genome, "chrF.fa"), "w") as fh:
            fh.write(f">chrF\nTTTA{self.GUIDE}ACGT\n")
        five = self.index("TTTV_2_five", contigs=("chrA",), at_start=True)
        self.assertEqual(index_coverage.problems(five, genome, five_prime_search=False),
                         [("chrF.fa", "no tree")])
        self.assertFalse(crispritz_bins.holds_target(f"TTTA{self.GUIDE}ACGT", "TTTV", 24, False,
                                                     2, 2))

    def test_a_folder_of_two_formats_is_refused(self):
        genome = self.genome()
        mixed = self.index("NGG_2_g", contigs=("chrA",))
        one_target_bin(os.path.join(mixed, "NGG_chrE_1.bin"), "NGG", fmt=1)
        found = index_coverage.problems(mixed, genome)
        self.assertEqual(found, [("", "holds trees of index formats 1 and 2, which no build "
                                      "writes together: it was put together from two builds")])
        self.assertEqual(index_coverage.message(mixed, found, "build it again."), (
            f"index 'NGG_2_g' ({mixed}) holds trees of index formats 1 and 2, which no build "
            f"writes together: it was put together from two builds: a tree of format 1 lacks "
            f"the sites at contig ends and N runs that one of format 2 holds, and either build "
            f"may lack trees of the other's FASTAs. Remove it (rm -rf '{mixed}') and build it "
            f"again."))


class AWideFlankOfFormat2NeedsNoWarning(_Index):
    """A PAM-less index of format 2 was built since CRISPRitz kept a target whose flank
    holds an N, and holds the PAMs within K - L nt of a contig end (B4): complete-search says
    neither of it under --prune-pam, as it does of one of format 1."""

    def test_no_note_and_no_warning(self):
        genome = os.path.join(self.tmp, "Genomes", "g")
        os.makedirs(genome)
        with open(os.path.join(genome, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 20 + "N" * 30 + "ACGT" * 20 + "\n")
        k8 = self.index("NNNNNNNN_3_g")
        self.assertIs(prune_flank.built_before_n_fix(k8, genome), False)
        indels = os.path.join(self.tmp, "lib", "NNNNNNNN_3_g+ds_INDELS")
        one_target_bin(os.path.join(indels, "NNNNNNNN_3_fakechr1", "NNNNNNNN_fakechr1_1.bin"),
                       "NNNNNNNN")
        self.assertIs(prune_flank.indels_built_before_n_fix(indels), False)
        crisprme = _load_crisprme()
        self.assertEqual(crisprme._wider_flank_notes([k8], "NRG", genome), [])
        old = self.index("NNNNNNNN_3_h", fmt=1)
        notes = crisprme._wider_flank_notes([k8, old], "NRG", genome)
        self.assertTrue(notes[0].startswith("NOTE [complete-search]: --prune-pam serves the 3-nt "
                                            "PAM NRG from NNNNNNNN_3_h"), notes)
        self.assertIn("nor has one built by CRISPRitz 2.10.0 or later (index format 2)", notes[0])


class TheNoteOfAWideFlankOfFormat1SaysWhatTheEngineMisreports(_Index):
    """A - strand PAM exactly K - L nt from a contig's start lies, in a PAM-less index of
    format 1 whose flank is K nt, in the window at index 0, which CRISPRitz 2.9.1's search read
    as a + strand one: it reported that PAM on the + strand at another position. 2.10.0's search
    reads index 0 of a 3' index right in either format (B1): of the published 8-nt-flank index
    it gives the dedicated index's row. complete-search's --prune-pam note of such an index said
    the misreport under 2.10.0 too, on stdout and in the run's Params.txt; it says it only under
    a CRISPRitz that reads format 1 alone (2.9.1), and the limits of the build under either."""

    def test_the_misreport_is_said_under_crispritz_2_9_1_alone(self):
        genome = os.path.join(self.tmp, "Genomes", "g")
        os.makedirs(genome)
        with open(os.path.join(genome, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 40 + "\n")
        k8 = self.index("NNNNNNNN_3_g", fmt=1)
        crisprme = _load_crisprme()
        limits = ("NOTE [complete-search]: --prune-pam serves the 3-nt PAM NGG from "
                  "NNNNNNNN_3_g, whose stored flank is 8 nt: a + strand PAM with 1 to 5 nt after "
                  "it to the end of a contig, and a - strand PAM in the first 5 nt of a contig, "
                  "are not indexed")
        misreport = ", and one 5 nt from the start is reported on the + strand at another position"
        rest = (" (METHODS section 9). An index whose flank is the PAM's width, or a dedicated "
                "one, has no such limit, nor has one built by CRISPRitz 2.10.0 or later (index "
                "format 2).")
        for stated, version, note in ((2, "2.10.0", limits + rest),
                                      (None, "2.9.1", limits + misreport + rest)):
            with self.subTest(engine=version):
                self.engine(stated)
                self.assertEqual(index_fit.engine_version(), version)
                self.assertEqual(crisprme._wider_flank_notes([k8], "NGG", genome), [note])


class AnEmptySubIndexRecordsItsFormat(empty._Fixture):
    """pool_index_indels.py records a sub-index buildTST built with no target, by buildTST's
    own rule: the rule of the format the installed CRISPRitz writes (the one that has just
    built it), and the sub-index records that format, the only record of it there is. With
    CRISPRitz 2.10.0 a fake-indel contig whose only site lies at its end holds a target: its
    .bin missing is a lost one, not a build with no target (2.8.1 recorded it as one, by the
    rule of format 1, and every search then skipped the chromosome's indels)."""

    run_index = empty.IndexSide.run_index
    sub = empty.IndexSide.sub
    top = empty.IndexSide.top
    recorded = empty.IndexSide.recorded
    tearDown_dirs = empty.IndexSide.tearDown_dirs
    # a fake-indel contig whose only NRG site's guide starts at its first base (its pad
    # off the contig): a site of format 2, none of format 1
    EDGE = "CTAACAGTTGCTTTTATCAC" + "AGG"

    def engine_states(self, index_format):
        """The stub engine's crispritz.py states INDEX_FORMAT ``index_format`` (None: none, as
        CRISPRitz 2.9.1's)."""
        shebang, _, rest = empty.STUB.partition("\n")
        with open(os.path.join(self.bin, "crispritz.py"), "w") as fh:
            fh.write(f"{shebang}\n" + (f"INDEX_FORMAT = {index_format}\n" if index_format else "")
                     + rest)

    def meta(self, chrom):
        with open(os.path.join(self.sub(chrom), index_fit.INDEX_META)) as fh:
            return json.load(fh)

    def test_a_build_with_no_target_records_the_format_of_its_engine(self):
        for fmt in (2, None):
            with self.subTest(engine=fmt):
                self.tearDown_dirs()
                self.engine_states(fmt)
                res, _ = self.run_index({"chrA": f">fakechrA\n{empty.WITH_NRG}\n",
                                         "chrB": f">fakechrB\n{empty.NO_NRG}\n"},
                                        STUB_NO_BIN="fakechrB")
                self.assertEqual((res.returncode, res.stderr), (0, ""))
                self.assertEqual(self.recorded(), ["chrB"])
                self.assertEqual(sorted(os.listdir(self.sub("chrB"))),
                                 [index_fit.INDEX_META, empty.NO_TARGETS])
                self.assertEqual(self.meta("chrB")["index_format"], fmt if fmt == 2 else 1)
                self.assertEqual((self.meta("chrB")["build_pam"], self.meta("chrB")["side"],
                                  self.meta("chrB")["guide_len"]), ("NRG", "3'", 20))

    def test_a_lost_bin_of_a_contig_only_format_2_holds_a_site_of_is_not_recorded(self):
        self.engine_states(2)
        res, _ = self.run_index({"chrA": f">fakechrA\n{empty.WITH_NRG}\n",
                                 "chrB": f">fakechrB\n{self.EDGE}\n"}, STUB_NO_BIN="fakechrB")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.recorded(), [])
        self.assertNotIn("target in the indels", res.stdout)
        # with CRISPRitz 2.9.1, which writes format 1, that contig holds no target
        self.tearDown_dirs()
        self.engine_states(None)
        res, _ = self.run_index({"chrA": f">fakechrA\n{empty.WITH_NRG}\n",
                                 "chrB": f">fakechrB\n{self.EDGE}\n"}, STUB_NO_BIN="fakechrB")
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertEqual(self.recorded(), ["chrB"])
        self.assertEqual(self.meta("chrB")["index_format"], 1)


@unittest.skipUnless(empty._gnu_sed(), "pool_search_indels.py merges with GNU `sed -i`, as the "
                                       "pipeline image does")
class ASearchReadsAnEmptySubIndexByItsFormat(empty._Fixture):
    """pool_search_indels.py, on a sub-index with no .bin and no record of a build with no
    target (one an earlier release built), reads its fake-indel genome by the rule of the
    format the sub-index records, else that of the INDELS index's other sub-indexes, else
    format 2's: a contig whose only site format 2 holds is a lost .bin in an index of
    format 2, and none in one of format 1."""

    run_search = empty.SearchSide.run_search
    from_view = empty.SearchSide.from_view
    searched = empty.SearchSide.searched
    fake = empty.SearchSide.fake
    LOGS = {"chrA": ["chrA\t\trs1\t0.1\tA>AT\t5\tACGT"], "chrB": ["chrB\t\trs2\t0.2\tC>CG\t5\tACGT"]}

    def sibling(self, fmt):
        """chrA's sub-index, holding a whole .bin of index format ``fmt``."""
        one_target_bin(os.path.join(self.tmp, "work", "genome_library", "NRG_3_ref+ds_INDELS",
                                    "NRG_3_fakechrA", "NRG_fakechrA_1.bin"), "NRG",
                       max_bulges=3, fmt=fmt)

    def test_by_the_format_of_the_other_sub_indexes(self):
        for fmt, refused in ((1, False), (2, True)):
            with self.subTest(sibling=fmt):
                shutil.rmtree(os.path.join(self.tmp, "work"), ignore_errors=True)
                shutil.rmtree(os.path.join(self.tmp, "out"), ignore_errors=True)
                self.sibling(fmt)
                res, _ = self.run_search(self.LOGS, {"chrB": False},
                                         fakes={"chrB": AnEmptySubIndexRecordsItsFormat.EDGE})
                if refused:
                    self.assertEqual(res.returncode, 1, res.stdout + res.stderr)
                    self.assertIn(f"{self.fake('chrB')}/fakechrB.fa holds a target of NRG: the "
                                  f"index is incomplete", res.stderr)
                else:
                    self.assertEqual((res.returncode, res.stderr), (0, ""))
                    self.assertIn("No NRG target in the indels of chrB: skipping its indel "
                                  "search", res.stdout)

    def test_by_the_other_sub_indexes_from_complete_tests_view(self):
        """complete-test --chrom searches from a view that links the chromosome's sub-index
        alone: its siblings are read where it really is."""
        for fmt, refused in ((1, False), (2, True)):
            with self.subTest(sibling=fmt):
                shutil.rmtree(os.path.join(self.tmp, "work"), ignore_errors=True)
                shutil.rmtree(os.path.join(self.tmp, "out"), ignore_errors=True)
                self.sibling(fmt)
                res, _ = self.run_search(self.LOGS, {"chrB": False}, view="chrB",
                                         fakes={"chrB": AnEmptySubIndexRecordsItsFormat.EDGE})
                self.assertEqual(res.returncode, 1 if refused else 0, res.stdout + res.stderr)
                self.assertEqual("holds a target of NRG: the index is incomplete" in res.stderr,
                                 refused, res.stderr)

    def test_by_its_own_record(self):
        self.sibling(2)
        sub = os.path.join(self.tmp, "work", "genome_library", "NRG_3_ref+ds_INDELS",
                           "NRG_3_fakechrB")
        os.makedirs(sub)
        with open(os.path.join(sub, index_fit.INDEX_META), "w") as fh:
            json.dump({"index_format": 1}, fh)
        res, _ = self.run_search(self.LOGS, {}, fakes={"chrB": AnEmptySubIndexRecordsItsFormat.EDGE})
        self.assertEqual((res.returncode, res.stderr), (0, ""))
        self.assertIn("No NRG target in the indels of chrB: skipping its indel search", res.stdout)


# The stub engine of test_build_index_only_atomic, as CRISPRitz 2.10.0: its crispritz.py
# states VERSION 2.10.0 and INDEX_FORMAT 2, and index-genome writes each .bin in index
# format 2 (its one target's flags byte, where format 1 has its '0', holds the whole pad).
STUB_2_10 = atomic.STUB_ENGINE.replace("\n", '\nVERSION = "2.10.0"\nINDEX_FORMAT = 2\n', 1).replace(
    '''def note(line):''', '''def tst2(pam, guide, seq, nb):
    """tst's .bin in index format 2: the header, then the body with the target's flags
    byte in place of its '0', then the magic again."""
    one, k = tst(pam, guide, seq), (len(pam) + 1) // 2
    head = b"CZI\\xfe" + struct.pack("=6i", 2, 28 + len(pam), 2 if set(pam) == {"N"} else 0,
                                    int(guide), len(pam), int(nb)) + pam.encode()
    return head + one[:12 + k] + bytes([min(int(nb), 31)]) + one[13 + k:] + b"CZI\\xfe"


def note(line):''', 1).replace(
    '''        data = tst(pam, len(motif) - abs(size), seq)
''', '''        data = tst2(pam, len(motif) - abs(size), seq, nb)
''', 1)


class BuildIndexOnlyTakesFormat2(atomic._Fixture):
    """build-index-only with a CRISPRitz that writes index format 2 and says so: 2.8.1 refused
    each index it built ("this CRISPRme reads index format 1 only"). The reference, variant
    and INDELS indexes now go into place, each recording what it holds, and a search takes
    them."""

    @classmethod
    def setUpClass(cls):
        atomic.setUpModule()  # (its stub Bio package, which crisprme.py imports)
        cls.addClassCleanup(atomic.tearDownModule)

    def fresh(self):
        super().fresh()
        path = os.path.join(self.bin, "crispritz.py")
        with open(path, "w") as fh:
            fh.write(STUB_2_10.replace("#!PYTHON", "#!" + sys.executable, 1))
        os.chmod(path, 0o755)

    def test_its_indexes_go_into_place_and_record_their_format(self):
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(err, "")
        for name in (atomic.REF, atomic.VAR, atomic.INDELS):
            bins = self.bins(name)
            self.assertTrue(bins, name)
            for b in bins:
                path = os.path.join(self.lib(name), b)
                self.assertEqual(index_fit.bin_format(path), 2, path)
                self.assertIsNone(crispritz_bins.bin_problem(path), path)
            folders = sorted({os.path.dirname(os.path.join(self.lib(name), b)) for b in bins})
            for folder in folders:
                meta = index_fit.read_index_meta(folder)
                self.assertEqual({k: meta[k] for k in ("index_format", "build_pam", "side",
                                                       "flank_len", "guide_len", "max_bulges")},
                                 {"index_format": 2, "build_pam": "NRG", "side": "3'",
                                  "flank_len": 3, "guide_len": 20, "max_bulges": 2}, folder)
                self.assertEqual(meta["built_by"],
                                 {"crispritz": "2.10.0", "crisprme": crisprme_version()})
        # a search takes them: it stops at the stub's search, past every check of an index
        rc, out, err = self.run_search()
        with open(self.path("Results", "search", "log_error.txt")) as fh:
            log = fh.read()
        self.assertNotIn("was built by a newer CRISPRitz", log)
        self.assertIn("stub crispritz.py: no search here", log)
        # and a run that finds them in place keeps them, and their records
        before = {n: index_fit.read_index_meta(self.lib(n)) for n in (atomic.REF, atomic.VAR)}
        rc, out, err = self.run_bio()
        self.assertEqual((rc, err), (0, ""))
        self.assertIn("Reference index already present", out)
        self.assertEqual({n: index_fit.read_index_meta(self.lib(n)) for n in before}, before)


class BuildIndexOnlyTakesTheSideInTheHeader(atomic._Fixture):
    """build-index-only read the side of the guide an index in place was built for from its
    .pam_build alone. Over a 5' PAM-less index of format 2 that had none (crispritz.py built
    it, or it was downloaded) it recorded the shipped 3' PAMs/20bp-NNN-NO-PAM.txt and said
    'already present', at exit 0; and it refused a 3' one of flank 4 with none, given a 3' PAM
    file of that flank, as 'a PAM-less index that records no side of the guide'. The header
    of its .bin files records the side, and is read first."""

    PAM_3P = os.path.join(atomic.REPO, "PAMs", "20bp-NNN-NO-PAM.txt")

    @classmethod
    def setUpClass(cls):
        atomic.setUpModule()  # (its stub Bio package, which crisprme.py imports)
        cls.addClassCleanup(atomic.tearDownModule)

    def fresh(self):
        super().fresh()
        path = os.path.join(self.bin, "crispritz.py")
        with open(path, "w") as fh:
            fh.write(STUB_2_10.replace("#!PYTHON", "#!" + sys.executable, 1))
        os.chmod(path, 0o755)

    def in_place(self, name, at_start):
        """A PAM-less index of format 2 for a 20-nt guide in genome_library, as crispritz.py
        builds it for a PAM on the 5' side of the guide (``at_start``) or the 3' side: a tree
        a chromosome, and no .pam_build."""
        pam = name.split("_")[0]
        for chrom in atomic.GENOME:
            one_target_bin(os.path.join(self.lib(name), f"{pam}_{chrom}_1.bin"), pam, 20, 2,
                           at_start)
        return self.lib(name)

    def build(self, pam):
        return self.finish(self.start(command=[
            "build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam", pam,
            "--bDNA", "1", "--bRNA", "1", "--thread", "2", "--path", self.work]))

    def test_a_pam_file_of_the_other_side_is_not_recorded_on_it(self):
        """Since the names of 2.10.0 (the defaults track's E2) a PAM file that says its guide
        length builds its own index beside it, under the name that says its side; one of its
        PAM alone, which says none and so takes the older name, is still refused, before the
        engine runs. Either says that the header records the other side, and neither records
        the PAM file on it."""
        index = self.in_place("NNN_2_hg38", at_start=True)
        alone = self.path("PAMs", "3bp-NNN-alone.txt")
        atomic._write(alone, "NNN 3\n")
        # (with no .pam_build, and with the 3' one a run of 2.8.1 recorded over it)
        for record in (None, "20bp-NNN-NO-PAM 23 3\n"):
            with self.subTest(record=record):
                if record:
                    atomic._write(os.path.join(index, ".pam_build"), record)
                shutil.rmtree(self.lib("20-NNN_2_hg38"), ignore_errors=True)
                rc, out, err = self.build(self.PAM_3P)
                self.assertEqual((rc, err), (0, ""), out)
                self.assertIn(
                    f"{index} is an index for a PAM on the 5' side of the guide (as the header of "
                    "its .bin files records): it is left as it is, and this index is built beside "
                    "it, as 20-NNN_2_hg38.\n", out)
                self.assertNotIn(f"already present: {index}", out)
                ran = len(self.calls())  # (the 3' build ran the engine, for its own index)
                rc, out, err = self.build(alone)
                self.assertEqual(rc, 1, out)
                self.assertEqual(self.calls()[ran:], [], "the refused build ran the engine")
                mv = f"mv '{index}' '{self.lib('NNN-20_2_hg38')}'"
                self.assertEqual(err, (
                    f"Error: {index} is already there, a PAM-less index for a PAM on the 5' side "
                    f"of the guide (as the header of its .bin files records), and {alone} puts "
                    "the PAM on the 3' side. 3bp-NNN-alone.txt says no guide length, so the index "
                    f"for it takes the same name: name this one for what it is ({mv}), which "
                    "records it, or build it in another working directory (--path).\n"))
                self.assertNotIn("already present", out)
                if record:
                    with open(os.path.join(index, ".pam_build")) as fh:
                        self.assertEqual(fh.read(), record)
                else:
                    self.assertFalse(os.path.exists(os.path.join(index, ".pam_build")))

    def test_one_of_the_side_its_header_records_is_kept_and_recorded(self):
        index = self.in_place("NNNN_2_hg38", at_start=False)
        pam = self.path("PAMs", "24bp-NNNN-3p.txt")
        atomic._write(pam, "N" * 24 + " 4\n")
        rc, out, err = self.build(pam)
        self.assertEqual((rc, err), (0, ""))
        self.assertIn(f"Reference index already present: {index}", out)
        with open(os.path.join(index, ".pam_build")) as fh:
            self.assertEqual(fh.read(), "24bp-NNNN-3p 24 4\n")
        self.assertEqual(index_fit.read_index_meta(index)["side"], "3'")
        self.assertEqual(self.calls(), [])


class AFormat1IndexIsBuiltAgain(atomic._Fixture):
    """A search by a CRISPRitz that writes index format 2 notes each index of format 1 it
    takes, as lacking the sites at contig ends and N runs: a variant search with bulges, its
    reference, variant and INDELS indexes. Each note said to build its index again with
    build-index-only, or to download it again, and neither gave format 2: build-index-only
    keeps an index it finds in place ('already present'), and download fetches the archives
    published until the rebuild, all of format 1. Each note now names its index, to be moved
    out of the way first, where this install builds it again; that done, build-index-only
    builds the three in format 2, and the next search notes none. An install of downloaded
    indexes holds no VCFs and no enriched genome, from which build-index-only builds the
    variant and INDELS indexes: their notes say so, and send neither out of the way. Of an
    index of an --index-path library other than the working directory's genome_library, which
    a search takes before one built again there, a note says to replace it with one built
    again, or to run the search without --index-path where the search then builds them."""

    PREFIX = TheNoteMovesOnlyAnIndexThisInstallBuildsAgain.PREFIX
    ARCHIVES = TheNoteMovesOnlyAnIndexThisInstallBuildsAgain.ARCHIVES

    @classmethod
    def setUpClass(cls):
        atomic.setUpModule()  # (its stub Bio package, which crisprme.py imports)
        cls.addClassCleanup(atomic.tearDownModule)

    def notes(self, *options, vcf=True):
        """The notes a search gives of an index of format 1, on its stdout (log_verbose.txt): a
        variant search with options, or without vcf a reference one."""
        if vcf:
            self.run_search(*options)
        else:
            shutil.rmtree(self.path("Results", "search"), ignore_errors=True)
            atomic._write(self.path("guides.txt"), atomic.GENOME["chr1"][:20] + "NNN\n")
            self.finish(self.start(command=[
                "complete-search", "--genome", self.path("Genomes", "hg38"), "--guide",
                self.path("guides.txt"), "--pam", atomic.PAM, "--mm", "4", "--bDNA", "1",
                "--bRNA", "1", "--output", "search", "--thread", "2", *options],
                CRISPRME_SKIP_SCORER_ENV="1"))
        return [line for line in self.search_log("log_verbose.txt").splitlines()
                if line.startswith("NOTE: ") and " is in index format 1 " in line]

    def search_log(self, name):
        """The last search's log ``name``, empty where complete-search launched none."""
        path = self.path("Results", "search", name)
        if not os.path.isfile(path):
            return ""
        with open(path) as fh:
            return fh.read()

    def stopped(self):
        """Where the last search stopped: '' at the stub's search, past every check of its
        indexes, else what it said (its log_error.txt, then complete-search's stderr)."""
        log = self.search_log("log_error.txt")
        return "" if "stub crispritz.py: no search here" in log else log + self.stderr()

    def formats(self, name, lib=None):
        top = os.path.join(lib, name) if lib else self.lib(name)
        return {index_fit.bin_format(os.path.join(d, f))
                for d, _, fs in os.walk(top) for f in fs if f.endswith(".bin")}

    def format1_install(self):
        """The three indexes as CRISPRme 2.8.1 built them, with CRISPRitz 2.9.1 (index format
        1), which CRISPRitz 2.10.0 then replaces."""
        rc, out, err = self.run_bio()  # by the stub of CRISPRitz 2.9.1: index format 1
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual([self.formats(n) for n in (atomic.REF, atomic.VAR, atomic.INDELS)],
                         [{1}] * 3)
        path = os.path.join(self.bin, "crispritz.py")
        with open(path, "w") as fh:
            fh.write(STUB_2_10.replace("#!PYTHON", "#!" + sys.executable, 1))
        os.chmod(path, 0o755)

    def move(self, *names):
        """Each index ``names`` names, moved out of the way, as its note says."""
        os.makedirs(os.path.join(self.root, "format1"), exist_ok=True)
        for name in names:
            os.rename(self.lib(name), os.path.join(self.root, "format1", name))

    def test_moved_out_of_the_way_as_its_note_says_it_is_built_in_format_2(self):
        names = (atomic.REF, atomic.VAR, atomic.INDELS)
        self.format1_install()
        notes = self.notes()
        self.assertEqual(len(notes), 3, notes)
        # build-index-only keeps each where it is, in format 1
        rc, out, err = self.run_bio()
        self.assertEqual((rc, err), (0, ""))
        for kept in ("Reference index", "Variant index", "Indels index"):
            self.assertIn(f"{kept} already present", out)
        self.assertEqual([self.formats(n) for n in names], [{1}] * 3)
        # so each note names its index, to be moved out of the way first, and what it is
        # built again from: the variant and INDELS indexes, the enriched genomes in place
        moved = [re.search(r"\. To search them, move (.+) out of the way and build it again with "
                           r"crisprme\.py build-index-only \(which keeps an index it finds in "
                           r"place\)(, from Genomes/hg38\+ds(_INDELS)?)?\.\Z", note)
                 for note in notes]
        self.assertTrue(all(moved), notes)
        self.assertEqual(sorted((os.path.realpath(m.group(1)), m.group(2) or "") for m in moved),
                         sorted((os.path.realpath(self.lib(n)), source) for n, source in (
                             (atomic.REF, ""), (atomic.VAR, ", from Genomes/hg38+ds"),
                             (atomic.INDELS, ", from Genomes/hg38+ds_INDELS"))))
        self.move(*(os.path.basename(m.group(1)) for m in moved))
        rc, out, err = self.run_bio()
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual([self.formats(n) for n in names], [{2}] * 3)
        self.assertEqual(self.notes(), [])

    def test_of_a_download_only_the_reference_index_is_moved_out_of_the_way(self):
        """The README's Quickstart installs the indexes and the dictionaries (download --what
        index), and no VCFs and no enriched genome: build-index-only then refuses the dataset,
        and a variant index the note had sent out of the way stayed there."""
        self.format1_install()
        vcfs = os.path.join(self.root, "VCFs.ds")
        shutil.move(self.path("VCFs", "ds"), vcfs)
        for genome in ("hg38+ds", "hg38+ds_INDELS"):
            shutil.rmtree(self.path("Genomes", genome))
        notes = self.notes()
        prefix = ("NOTE: {} is in index format 1 (built by CRISPRitz 2.9.1 or earlier): it lacks "
                  "the sites whose guide lies within a few bases of a contig end or an N run, "
                  "which an index built by CRISPRitz 2.10.0 or later holds (METHODS section 9). ")
        archives = ("Where ds is published, its archives rebuilt in index format 2 hold them: "
                    "building it again with crisprme.py build-index-only takes the dataset's ")
        self.assertEqual(notes, [
            prefix.format(atomic.REF) + f"To search them, move {self.lib(atomic.REF)} out of the "
            "way and build it again with crisprme.py build-index-only (which keeps an index it "
            "finds in place).",
            prefix.format(atomic.VAR) + archives + "VCFs (VCFs/ds), which this install does not "
            "hold.",
            prefix.format(atomic.INDELS) + archives + "VCFs (VCFs/ds), which this install does "
            "not hold."])
        # the reference index, moved out of the way, is built in format 2
        self.move(atomic.REF)
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual([self.formats(n) for n in (atomic.REF, atomic.VAR, atomic.INDELS)],
                         [{2}, {1}, {1}])
        self.assertEqual(self.notes(), notes[1:])
        # with the VCFs back and no enriched genome, building them again enriches the whole
        # dataset again: the notes say so, and send neither out of the way
        shutil.move(vcfs, self.path("VCFs", "ds"))
        self.assertEqual(self.notes(), [
            prefix.format(name) + archives + "enriched genome (Genomes/hg38+ds), which this "
            "install does not hold, or else its VCFs (VCFs/ds), from which the whole dataset is "
            "enriched again first." for name in (atomic.VAR, atomic.INDELS)])


    def library(self, *names):
        """An --index-path library outside the working directory, not named genome_library,
        holding the indexes ``names`` names, moved there."""
        lib = os.path.realpath(os.path.join(self.root, "staged"))
        os.makedirs(lib, exist_ok=True)
        for name in names:
            os.rename(self.lib(name), os.path.join(lib, name))
        return lib

    def test_of_an_index_path_library_each_is_replaced_there(self):
        """A search takes its reference index only from the --index-path library, and its INDELS
        index only from beside its variant index there: built again in the working directory, as
        the note said, neither was searched, and the next search stopped."""
        names = (atomic.REF, atomic.VAR, atomic.INDELS)
        self.format1_install()
        lib = self.library(*names)
        notes = self.notes("--index-path", lib)
        self.assertEqual(self.stopped(), "")
        self.assertEqual(notes, [self.PREFIX.format(name) + _replaced(
            os.path.join(lib, name), source, True) for name, source in (
                (atomic.REF, None), (atomic.VAR, "Genomes/hg38+ds"),
                (atomic.INDELS, "Genomes/hg38+ds_INDELS"))])
        # as the note said: each moved out of the way, and build-index-only run here, which
        # builds the three in format 2 into the working directory's genome_library
        away = os.path.join(self.root, "format1")
        os.makedirs(away)
        os.rename(os.path.join(lib, atomic.REF), os.path.join(away, atomic.REF))
        rc, out, err = self.run_bio()
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual([self.formats(n) for n in names], [{2}] * 3)
        self.notes("--index-path", lib)
        self.assertIn("ERROR: no matching reference index under --index-path", self.stopped())
        os.rename(os.path.join(away, atomic.REF), os.path.join(lib, atomic.REF))
        os.rename(os.path.join(lib, atomic.INDELS), os.path.join(away, atomic.INDELS))
        self.notes("--index-path", lib)
        self.assertIn(f"ERROR: INDELS index '{atomic.INDELS}' is not beside its variant index in "
                      f"'{lib}' (--index-path), and a search builds no index there.",
                      self.stopped())
        # as it says now: each replaced with the one built again
        os.rename(os.path.join(away, atomic.INDELS), os.path.join(lib, atomic.INDELS))
        for name in names:
            shutil.rmtree(os.path.join(lib, name))
            os.rename(self.lib(name), os.path.join(lib, name))
        self.assertEqual([self.formats(n, lib) for n in names], [{2}] * 3)
        self.assertEqual(self.notes("--index-path", lib), [])
        self.assertEqual(self.stopped(), "")

    def test_of_an_index_path_library_or_searched_without_it(self):
        names = (atomic.REF, atomic.VAR, atomic.INDELS)
        self.format1_install()
        lib = self.library(*names)
        notes = self.notes("--index-path", lib)
        self.assertEqual([note.endswith(", or run this search without --index-path.")
                          for note in notes], [True] * 3, notes)
        # without --index-path, the search builds the three in format 2 here, from the
        # genome folder, the enriched genome and the fake-indel genome
        self.assertEqual(self.notes(), [])
        self.assertEqual(self.stopped(), "")
        self.assertEqual([self.formats(n) for n in names], [{2}] * 3)

    def test_of_an_index_path_library_of_a_download_only_the_reference_index_is_replaced(self):
        """The README's staged library, of downloaded indexes: the working directory holds the
        dictionaries, and no VCFs and no enriched genome."""
        names = (atomic.REF, atomic.VAR, atomic.INDELS)
        self.format1_install()
        shutil.move(self.path("VCFs", "ds"), os.path.join(self.root, "VCFs.ds"))
        for genome in ("hg38+ds", "hg38+ds_INDELS"):
            shutil.rmtree(self.path("Genomes", genome))
        lib = self.library(*names)
        notes = self.notes("--index-path", lib)
        self.assertEqual(self.stopped(), "")
        self.assertEqual(notes, [
            self.PREFIX.format(atomic.REF) + _replaced(os.path.join(lib, atomic.REF)),
            self.PREFIX.format(atomic.VAR) + self.ARCHIVES + "VCFs (VCFs/ds), which this install "
            "does not hold.",
            self.PREFIX.format(atomic.INDELS) + self.ARCHIVES + "VCFs (VCFs/ds), which this "
            "install does not hold."])
        # (without --index-path the search has no variant index, and builds none)
        self.notes()
        self.assertIn("Error: no variant index here serves this search", self.stopped())
        # the reference index, replaced with one build-index-only builds again
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual((rc, err), (0, ""))
        shutil.rmtree(os.path.join(lib, atomic.REF))
        os.rename(self.lib(atomic.REF), os.path.join(lib, atomic.REF))
        self.assertEqual([self.formats(n, lib) for n in names], [{2}, {1}, {1}])
        self.assertEqual(self.notes("--index-path", lib), notes[1:])
        self.assertEqual(self.stopped(), "")

    def test_of_an_index_path_library_in_a_reference_search(self):
        self.format1_install()
        lib = self.library(atomic.REF)
        note = self.PREFIX.format(atomic.REF) + _replaced(os.path.join(lib, atomic.REF),
                                                          without=True)
        self.assertEqual(self.notes("--index-path", lib, vcf=False), [note])
        self.assertEqual(self.stopped(), "")
        # without --index-path, the search builds it in format 2 here
        self.assertEqual(self.notes(vcf=False), [])
        self.assertEqual(self.stopped(), "")
        self.assertEqual(self.formats(atomic.REF), {2})
        # (which the search with --index-path does not take)
        self.assertEqual(self.notes("--index-path", lib, vcf=False), [note])
        # or replaced with it
        shutil.rmtree(os.path.join(lib, atomic.REF))
        os.rename(self.lib(atomic.REF), os.path.join(lib, atomic.REF))
        self.assertEqual(self.notes("--index-path", lib, vcf=False), [])
        self.assertEqual(self.stopped(), "")

    def test_of_an_index_path_library_under_prune_pam_the_search_keeps_it(self):
        """Under --prune-pam the note said to run the search without --index-path "and
        --prune-pam (which needs it)": since 2.10.0 --prune-pam needs no --index-path, and the
        search without --index-path alone builds the index in format 2 here, and prunes."""
        self.format1_install()
        with open(os.path.join(self.bin, "crispritz.py"), "a") as fh:
            fh.write("# takes --prune-pam (the search refuses to prune with an engine whose source "
                     "never names it)\n")
        lib = self.library(atomic.REF)
        note = self.PREFIX.format(atomic.REF) + _replaced(os.path.join(lib, atomic.REF),
                                                          without=True)
        self.assertEqual(self.notes("--index-path", lib, "--prune-pam", vcf=False), [note])
        self.assertEqual(self.stopped(), "")
        # as it says: without --index-path, --prune-pam kept
        self.assertEqual(self.notes("--prune-pam", vcf=False), [])
        self.assertEqual(self.stopped(), "")
        self.assertEqual(self.formats(atomic.REF), {2})
        with open(self.log) as fh:  # (the engine is handed the index, or a view of its trees)
            searches = [line.split() for line in fh if line.startswith("search ")]
        self.assertEqual(os.path.basename(searches[-1][1].rstrip("/")), atomic.REF)
        self.assertIn("--prune-pam", searches[-1])


class AnArchiveSaysItsFormat(_Index):
    """publish-index puts the index format of an index's .bin files, and the record of what it
    holds, on the archive's manifest, and refuses an index of two formats; download says at
    once, on stdout, that an archive the installed CRISPRitz cannot read is installed but
    refused by every search."""

    NAME = "NGG_2_hg38+hg38_ds"

    def layout(self, fmt=2, mixed=False):
        root = os.path.join(self.tmp, "src")
        idx = os.path.join(root, "genome_library", self.NAME)
        one_target_bin(os.path.join(idx, "NGG_chr22_1.bin"), "NGG", fmt=fmt)
        if mixed:
            one_target_bin(os.path.join(idx, "NGG_chr21_1.bin"), "NGG", fmt=3 - fmt)
        one_target_bin(os.path.join(idx + "_INDELS", "NGG_2_fakechr22", "NGG_fakechr22_1.bin"),
                       "NGG", fmt=fmt)
        for name, text in (("log_indels_hg38_ds/logchr22.txt", "chr22\n"),
                           ("dictionaries_hg38_ds/my_dict_chr22.json", "{}")):
            os.makedirs(os.path.dirname(os.path.join(root, "Dictionaries", name)), exist_ok=True)
            with open(os.path.join(root, "Dictionaries", name), "w") as fh:
                fh.write(text)
        return idx

    def manifest(self, api):
        """The manifest of the index's own archive (an index of format 2 goes up after its
        panel's companions, its Dictionaries bundle first: crisprme_hf.publish_index)."""
        import tarfile
        archive = next(local for local, remote in api.uploads if remote.endswith(".tar.gz")
                       and not dpub.hf._companion_kind(os.path.basename(remote)[:-7]))
        with tarfile.open(archive) as tf:
            return json.load(tf.extractfile("manifest.json"))

    @staticmethod
    def publish(idx, sink):
        """test_dictless_publish's offline publish, its notes kept off the test's output."""
        import contextlib
        import io

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return dpub._publish(idx, sink)

    def test_publish_records_the_format_and_the_record(self):
        self.engine(2)
        idx = self.layout()
        index_fit.write_index_meta(idx, built=True)
        sink = os.path.join(self.tmp, "sink")
        os.makedirs(sink)
        manifest = self.manifest(self.publish(idx, sink))
        self.assertEqual(manifest["index_format"], 2)
        self.assertEqual(manifest["index_meta"]["index_format"], 2)
        self.assertEqual(manifest["index_meta"]["build_pam"], "NGG")
        with self.assertRaises(ValueError) as told:
            self.publish(self.layout_mixed(), sink)
        self.assertIn("holds .bin files of index formats 1 and 2: it was put together from two "
                      "builds", str(told.exception))

    def layout_mixed(self):
        shutil.rmtree(os.path.join(self.tmp, "src"))
        return self.layout(mixed=True)

    def test_publish_records_the_fake_indel_windows(self):
        """W-R1-N3.3: publish-index puts the fake-indel windows a variant index's indel logs
        record on the manifest (indel_window: [left, right] reference bases around each
        indel), the windows of the enrichment its INDELS index was built from; none where its
        logs hold no row, or windows of two widths (logs of two enrichments), and none for a
        reference index."""
        self.engine(2)
        idx = self.layout()
        logs = os.path.join(self.tmp, "src", "Dictionaries", "log_indels_hg38_ds")
        sink = os.path.join(self.tmp, "sink")
        os.makedirs(sink)

        def log(chrom, left, right):   # one insertion A>AT at POS 100, as add-variants logs it
            with open(os.path.join(logs, "log%s.txt" % chrom), "w") as fh:
                fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
                         "%s_%d-%d_1\tS1\trs1\t0.1\t%s_100_A_AT\t0,%d\tN\n"
                         % (chrom, 99 - left, 100 + right, chrom, left + 2 + right))

        self.assertNotIn("indel_window", self.manifest(self.publish(idx, sink)))  # no row
        log("chr22", 30, 30)
        self.assertEqual(self.manifest(self.publish(idx, sink))["indel_window"], [30, 30])
        log("chr22", 25, 27)
        self.assertEqual(self.manifest(self.publish(idx, sink))["indel_window"], [25, 27])
        log("chr21", 30, 30)
        self.assertNotIn("indel_window", self.manifest(self.publish(idx, sink)))
        ref = os.path.join(self.tmp, "src", "genome_library", "NGG_2_hg38")
        one_target_bin(os.path.join(ref, "NGG_chr22_1.bin"), "NGG", fmt=2)
        self.assertNotIn("indel_window", self.manifest(self.publish(ref, sink)))

    def test_download_says_at_once_an_archive_the_engine_cannot_read(self):
        import contextlib
        import io

        self.engine(2)
        idx = self.layout()
        sink = os.path.join(self.tmp, "sink")
        os.makedirs(sink)
        self.publish(idx, sink)
        orig = dpub.hf._hf_snapshot
        self.addCleanup(setattr, dpub.hf, "_hf_snapshot", orig)

        def fake_snapshot(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            os.makedirs(os.path.join(local_dir, "indexes"), exist_ok=True)
            for pat in allow_patterns:
                src = os.path.join(sink, os.path.basename(pat))
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(local_dir, "indexes", os.path.basename(pat)))
            return local_dir

        dpub.hf._hf_snapshot = fake_snapshot
        for engine, warned in ((None, True), (2, False)):
            with self.subTest(engine=engine):
                self.engine(engine)
                work = os.path.join(self.tmp, f"dst{engine}")
                os.makedirs(work)
                out, err = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                    dpub.hf.download_component("index", work, repo="local/test",
                                               index_name=self.NAME, token="x", genotypes=False)
                installed = os.path.join(work, "genome_library", self.NAME)
                self.assertTrue(os.path.isdir(installed))
                self.assertEqual(f"WARNING: {installed} was built by a newer CRISPRitz" in
                                 out.getvalue(), warned, out.getvalue())
                self.assertIn("index format 2", err.getvalue() if not warned else out.getvalue())


# The AppleDouble file ._<name> that macOS tar and cp write beside each file they copy onto a
# volume without extended attributes: its magic 0x00051607, its version 0x00020000, the filler
# "Mac OS X" and its entries, in 4096 bytes
APPLE_DOUBLE = bytes.fromhex("00051607000200004d6163204f5320582020202000020000000900000032000000"
                             "e00000000200000112000000000000000000").ljust(4096, b"\0")


def apple_double(path):
    """The ._<name> file macOS leaves beside the file ``path``."""
    with open(os.path.join(os.path.dirname(path), "._" + os.path.basename(path)), "wb") as fh:
        fh.write(APPLE_DOUBLE)


class AnAppleDoubleFileIsNoTree(_Index):
    """macOS tar and cp leave an AppleDouble file ._<name> beside each file they copy onto a
    volume without extended attributes. ._NGG_chr22_1.bin opens with no magic, and every reader
    of an index took it for one of its trees, of index format 1 and built for a 512-nt guide
    (its second int32): an index of format 2 beside one was refused as put together from two
    builds, with rm -rf as the remedy, by the search and by publish-index; under CRISPRitz 2.9.1
    its format went unseen; build-index-only kept an index built for too short a guide, and
    complete-search, the web and the search script took it as fitting any guide (CRISPRme 2.8.1
    did so too); its record was not written, or named the PAM "." of a 512-nt guide; and
    --prune-pam gave the note of a wide flank of format 1. The search script's probes read the
    first .bin find or sort lists, which the file may be, by the locale and the order of the
    file system (of the names' hashes on ext4 and APFS, else of their writing), so the 16 folders
    here hold trees of 16 names, each written before its ._ file in half of them and after it in
    the others. Every reader takes the index's trees alone now, the search hands the engine those
    alone (searchTST stops on a folder that holds such a file), and publish-index leaves the
    file out."""

    NAME = AnArchiveSaysItsFormat.NAME
    layout = AnArchiveSaysItsFormat.layout
    manifest = AnArchiveSaysItsFormat.manifest

    def trees(self, n=16, pam="TTTV", guide=21, at_start=True):
        """The trees of ``n`` indexes of one each, of index format 2, of a contig of a name of
        its own, the ._ file of each written before it in half of them and after it in the
        others."""
        found = []
        for i in range(n):
            tree = os.path.join(self.tmp, "lib", f"{pam}_2_g{i}", f"{pam}_chr{i + 1}_1.bin")
            os.makedirs(os.path.dirname(tree))
            if i % 2:
                apple_double(tree)
            one_target_bin(tree, pam, guide, 2, at_start)
            if not i % 2:
                apple_double(tree)
            found.append(tree)
        return found

    def test_index_fit_and_the_coverage_check_read_its_trees(self):
        genome = os.path.join(self.tmp, "Genomes", "g")
        os.makedirs(genome)
        with open(os.path.join(genome, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 10 + "\n")
        self.engine(2)
        five = self.index("TTTV_2_g", guide=21, at_start=True)
        old = self.index("NGG_2_g", fmt=1, guide=21)
        for folder, pam in ((five, "TTTV"), (old, "NGG")):
            apple_double(os.path.join(folder, f"{pam}_chr1_1.bin"))
        self.assertEqual((index_fit.bin_formats(five), index_fit.bin_formats(old)),
                         ({2: 1}, {1: 1}))
        self.assertEqual((index_fit.built_guide_len(five), index_fit.built_guide_len(old)),
                         (21, 21))
        self.assertIs(index_fit.built_five_prime(five), True)
        self.assertEqual(crispritz_bins.index_format(five), 2)
        self.assertEqual(index_coverage.problems(five, genome), [])
        self.assertEqual(crispritz_bins.incomplete_bins(five), [])
        self.assertIsNone(index_fit.write_index_meta(five))
        self.assertIsNone(index_fit.write_index_meta(old))
        self.assertEqual(index_fit.read_index_meta(five), {
            "index_format": 2, "build_pam": "TTTV", "side": "5'", "flank_len": 4,
            "guide_len": 21, "max_bulges": 2, "pamless": False, "four_allele": False})
        self.assertEqual(index_fit.read_index_meta(old), {
            "index_format": 1, "build_pam": "NGG", "flank_len": 3, "guide_len": 21,
            "max_bulges": 2, "pamless": False, "four_allele": False})
        self.engine(None)  # CRISPRitz 2.9.1, which reads format 1 alone
        self.assertEqual(index_fit.unreadable_bin(five), (os.path.join(five, "TTTV_chr1_1.bin"), 2))
        # a folder renamed: its trees are named for NGG alone
        renamed = os.path.join(self.tmp, "lib", "NRG_2_g")
        os.rename(old, renamed)
        self.assertEqual(index_coverage.problems(renamed, genome), [
            ("", "its trees are named for NGG, not NRG: the folder was renamed")])

    def test_the_search_script_reads_its_trees(self):
        # (a 23-nt guide with 1 + 1 bulges: an index built for 21 nt with N = 2 does not fit it)
        pam_line = "TTTV" + "N" * 23 + " -4"
        call = ('_index_guide_len "$1"; echo; _header_side "$1"; _index_fits_guide "$1" 2; '
                'echo "fits $?"; _check_index_format "$1"; echo went on')
        trees = self.trees()
        self.engine(2)
        for tree in trees:
            out = self.shell(call, os.path.dirname(tree), pam_line=pam_line)
            self.assertEqual((out.stdout, out.stderr), ("21\n5\nfits 1\nwent on\n", ""), tree)
        lib = os.path.dirname(os.path.dirname(trees[0]))
        for i in range(len(trees)):
            out = self.shell(f'_resolve_index_dir "$1" TTTV g{i} 2', lib, pam_line=pam_line)
            self.assertEqual((out.stdout, out.stderr), ("", ""), trees[i])
        self.engine(None)  # CRISPRitz 2.9.1, which reads format 1 alone
        for tree in trees:
            folder = os.path.dirname(tree)
            # (sorted as inside the container, where a ._ file sorts first)
            out = self.shell('export LC_ALL=C; _check_index_format "$1"', folder, pam_line=pam_line)
            self.assertEqual(out.returncode, 1, tree)
            self.assertIn(f"ERROR: {folder} was built by a newer CRISPRitz "
                          f"({os.path.basename(tree)} is in index format 2", out.stderr)

    def test_the_engine_is_handed_its_trees_alone(self):
        folder = self.index("NGG_2_g", contigs=("chr1", "chr2"))
        for c in ("chr1", "chr2"):
            apple_double(os.path.join(folder, f"NGG_{c}_1.bin"))
        out_dir = os.path.join(self.tmp, "out")
        script = f'output_folder="{out_dir}"\n{_functions("_bin_only_view")}\n_bin_only_view "$1"'
        out = subprocess.run(["bash", "-c", script, "bash", folder], capture_output=True, text=True)
        view = os.path.join(out_dir, ".binview", "NGG_2_g")
        self.assertEqual((out.stdout, out.stderr), (view, ""))
        self.assertEqual(sorted(os.listdir(view)), ["NGG_chr1_1.bin", "NGG_chr2_1.bin"])
        for c in ("chr1", "chr2"):
            os.remove(os.path.join(folder, f"._NGG_{c}_1.bin"))
        out = subprocess.run(["bash", "-c", script, "bash", folder], capture_output=True, text=True)
        self.assertEqual((out.stdout, out.stderr), (folder, ""))  # (its trees alone: as it is)

    def test_complete_search_and_the_web_read_its_trees(self):
        import test_bulge_cap

        crisprme = _load_crisprme()
        self.engine(2)
        short = self.index("NGG_2_toy")  # built for 20 nt: no 23-nt guide with 1 + 1 bulges
        apple_double(os.path.join(short, "NGG_chr1_1.bin"))
        self.assertIsNone(crisprme._index_the_search_takes(
            [os.path.dirname(short)], "NGG", "toy", 2, False, True, 23, 1, 1, False))
        five = self.index("TTTV_2_g", guide=21, at_start=True)  # (its geometry, in its header)
        apple_double(os.path.join(five, "TTTV_chr1_1.bin"))
        self.assertEqual(test_bulge_cap.pages_utils.index_build_pam(five), (None, 25, -4))
        genome = os.path.join(self.tmp, "Genomes", "g")
        os.makedirs(genome)
        with open(os.path.join(genome, "chr1.fa"), "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 20 + "N" * 30 + "ACGT" * 20 + "\n")
        k8 = self.index("NNNNNNNN_3_g")
        apple_double(os.path.join(k8, "NNNNNNNN_chr1_1.bin"))
        self.assertEqual(crisprme._wider_flank_notes([k8], "NRG", genome), [])

    def test_publish_index_takes_it_and_leaves_the_file_out(self):
        self.engine(2)
        idx = self.layout()
        apple_double(os.path.join(idx, "NGG_chr22_1.bin"))
        apple_double(os.path.join(idx + "_INDELS", "NGG_2_fakechr22", "NGG_fakechr22_1.bin"))
        sink = os.path.join(self.tmp, "sink")
        os.makedirs(sink)
        api = AnArchiveSaysItsFormat.publish(idx, sink)
        self.assertEqual(self.manifest(api)["index_format"], 2)
        # with tar and pigz, and with Python's tarfile alike (gzip stands in for pigz)
        import contextlib
        import io
        import tarfile
        from unittest import mock

        pigz = os.path.join(self.tmp, "pigz", "pigz")
        os.makedirs(os.path.dirname(pigz))
        with open(pigz, "w") as fh:
            fh.write('#!/bin/sh\nexec gzip "$@"\n')
        os.chmod(pigz, 0o755)
        which = shutil.which
        name = os.path.basename(idx)
        for found in (lambda program: pigz if program == "pigz" else which(program),
                      lambda program: None):
            tarball = os.path.join(self.tmp, f"{name}.tar.gz")
            with mock.patch.object(dpub.hf.shutil, "which", side_effect=found), \
                    contextlib.redirect_stderr(io.StringIO()):
                dpub.hf._make_index_tarball(tarball, idx, name, idx + "_INDELS", {"name": name})
            with tarfile.open(tarball) as tf:
                names = tf.getnames()
            self.assertIn(f"{name}/NGG_chr22_1.bin", names)
            self.assertIn(f"{name}_INDELS/NGG_2_fakechr22/NGG_fakechr22_1.bin", names)
            self.assertEqual([n for n in names if os.path.basename(n).startswith("._")], [])


class CheckBoundaryRowsChecksEachBlock(unittest.TestCase):
    """test/index_format2/check_boundary_rows.py checks that each row of a search's integrated
    results is its genome's bases at its coordinates, in each of the row's blocks: the
    alignments of the highest CFD score, of the fewest mismatches and bulges, and of the highest
    CRISPR-Bulge score, which a row reports apart. It named the third highest_CRISPR_Bulge,
    which the integrated results name highest_CRISPR_BULGE, and so never checked it. It checks
    each block resultIntegrator.py writes now, and a row that lacks one is wrong."""

    def test_a_wrong_alignment_in_any_block_is_told(self):
        import importlib.util

        path = os.path.join(os.path.dirname(HERE), "test", "index_format2",
                            "check_boundary_rows.py")
        spec = importlib.util.spec_from_file_location("check_boundary_rows", path)
        rows = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(rows)
        with open(os.path.join(HERE, "resultIntegrator.py")) as fh:
            blocks = list(dict.fromkeys(re.findall(r'"Start_coordinate_\(([^)]+)\)"', fh.read())))
        self.assertEqual(len(blocks), 3, blocks)
        seqs, recs = rows.genome()
        site = seqs["chrA"][2000:2023]  # the site planted inside chrA, on the + strand
        good = {"Chromosome": "chrA"}
        for block in blocks:
            good.update({f"Start_coordinate_({block})": "2000", f"Strand_({block})": "+",
                         f"Aligned_protospacer+PAM_REF_({block})": site,
                         f"Aligned_protospacer+PAM_ALT_({block})": "NA",
                         f"REF/ALT_origin_({block})": "ref",
                         f"Variant_info_genome_({block})": "NA"})
        self.assertEqual(rows.wrong_rows(seqs, recs, [good]), [])
        for block in blocks:
            with self.subTest(block=block):
                moved = dict(good, **{f"Start_coordinate_({block})": "2001"})  # a base off
                self.assertEqual([bad[0] for bad in rows.wrong_rows(seqs, recs, [moved])], [block])
                lacking = {k: v for k, v in good.items() if not k.endswith(f"_({block})")}
                self.assertEqual([bad[0] for bad in rows.wrong_rows(seqs, recs, [lacking])],
                                 [block])

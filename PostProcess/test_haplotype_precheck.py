"""The exact haplotype pre-check, its carriers-first twin on the indel side, and the K-skip
(2.10.0, haplotype_precheck).

A per-sample search skips a raw target no sample's genotypes can produce, and writes only its
reference row; a population-level search with --max-degeneracy K skips a raw target every row
of which rests on more than K ambiguous columns. Both must never skip a target a row comes
from, so each test here compares the module with a brute-force enumeration of every choice of
alleles, counted as the finalizer of new_simple_analysis.py and indel_snp_cis.resolve_target
count them:

  * sample_can_produce (over snp_columns) says True exactly when one sample, free to take at
    each variant column the reference base or any ALT it carries there, reaches the budget and
    the PAM: so it is True whenever a row of that sample's (a subset of those choices) can be;
  * degeneracy_floor is the fewest ambiguous columns (window_density.degeneracy) of any choice
    within the budget and the PAM, None when there is none;
  * indel_target_possible says True whenever one carrier has a choice resolve_target keeps;
  * and, through the real new_simple_analysis.py, the rows, the companions and both BEDs are
    the same with the pre-check on and off, a skipped target keeps its reference row, and the
    run says how many it skipped; a dense window two alignments share is one BED line.

Run with:
    cd PostProcess && python3 -m unittest test_haplotype_precheck -v
"""

import itertools
import os
import random
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import haplotype_precheck as hpc  # noqa: E402
import indel_snp_cis as isc  # noqa: E402
import window_density as wd  # noqa: E402

# the finalizer's PAM test (new_simple_analysis.iupac_code_set)
_SETS = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
_TAB = str.maketrans("ACGTacgt", "TGCAtgca")


def admits(letter, base):
    return base in _SETS[letter]


def revcomp(seq):
    return seq.translate(_TAB)[::-1]


def aligned(real_target, seq_prerevert, revert):
    """The finalizer's aligned target: reverse-complemented on the minus strand, the bulge
    gaps of ``real_target`` put back."""
    seq = revcomp(seq_prerevert) if revert else seq_prerevert
    t = list(seq)
    for pos, ch in enumerate(real_target):
        if ch == "-":
            t.insert(pos, "-")
    return t


def finalizer(real_target, seq_prerevert, revert, guide_no_pam, pam, sl):
    """(mismatch count, PAM ok) exactly as _finalize_observed_entry counts them."""
    pos_beg, pos_end, pam_begin, pam_end = sl
    t = aligned(real_target, seq_prerevert, revert)
    mm = sum(1 for i, ch in enumerate(t[pos_beg:pos_end]) if ch.upper() != guide_no_pam[i])
    ok = all(ch.upper() in _SETS[pam[i]] for i, ch in enumerate(t[pam_begin:pam_end]))
    return mm, ok, t


def random_case(rng, five_prime=False):
    """A raw target: its aligned guide and DNA (bulge gaps in either), the forward reference
    window, the strand, the PAM and slices, the variant columns and their ALTs, and samples
    carrying some of them."""
    if five_prime:
        pam, sl = "TTTV", (4, None, 0, 4)
    else:
        pam, sl = rng.choice(["NGG", "NRG", "NRR"]), (0, -3, -3, None)
    n_spacer = 20
    rna = rng.random() < 0.3  # a gap in the target (an RNA bulge)
    dna = rng.random() < 0.3  # a gap in the guide (a DNA bulge)
    spacer_guide = [rng.choice("ACGT") for _ in range(n_spacer)]
    spacer_target = ["x"] * n_spacer
    if rna:
        spacer_target[rng.randrange(2, n_spacer - 2)] = "-"
    if dna:
        while True:
            j = rng.randrange(2, n_spacer - 2)
            if spacer_target[j] != "-":
                spacer_guide.insert(j, "-")
                spacer_target.insert(j, "x")
                break
    if five_prime:
        guide = "TTTN" + "".join(spacer_guide)
        target_shape = ["x"] * 4 + spacer_target
    else:
        guide = "".join(spacer_guide) + "NNN"
        target_shape = spacer_target + ["x"] * 3
    n_bases = sum(1 for c in target_shape if c != "-")
    ref = "".join(rng.choice("ACGT") for _ in range(n_bases))
    revert = rng.random() < 0.5
    real_target = "".join("-" if c == "-" else "A" for c in target_shape)
    # the reference made close to the guide so that budgets are met often
    shown = aligned(real_target, ref, revert)
    guide_no_pam = guide[sl[0]:sl[1]]
    ps = list(range(len(shown)))[sl[0]:sl[1]]
    for t, j in enumerate(ps):
        if shown[j] != "-" and guide_no_pam[t] != "-" and rng.random() < 0.75:
            shown[j] = guide_no_pam[t]
    pm = list(range(len(shown)))[sl[2]:sl[3]]
    for t, j in enumerate(pm):
        if rng.random() < 0.8:
            shown[j] = rng.choice(_SETS[pam[t]])
    fwd = "".join(c for c in shown if c != "-")
    ref = revcomp(fwd) if revert else fwd
    cols = sorted(rng.sample(range(n_bases), rng.randint(1, min(7, n_bases))))
    variant_columns = []
    for c in cols:
        alts = [b for b in "ACGT" if b != ref[c]]
        rng.shuffle(alts)
        variant_columns.append((c, sorted(alts[:rng.randint(1, 3)])))
    samples = {}
    for s in range(rng.randint(0, 6)):
        carried = {}
        for c, alts in variant_columns:
            got = [a for a in alts if rng.random() < 0.35]
            if got:
                carried[c] = got
        samples["S%d" % s] = carried
    bulges = guide.count("-") + real_target.count("-")
    return dict(guide=guide, guide_no_pam=guide_no_pam, real_target=real_target, ref=ref,
                revert=revert, pam=pam, sl=sl, variant_columns=variant_columns,
                samples=samples, bulges=bulges, allowed=rng.randint(0, 4))


def choices_of(variant_columns, allowed_alts):
    """Every choice: at each column the reference (None) or one allowed ALT."""
    per = [[None] + allowed_alts.get(c, []) for c, _alts in variant_columns]
    for pick in itertools.product(*per):
        yield {c: a for (c, _alts), a in zip(variant_columns, pick) if a is not None}


def apply(ref, choice):
    seq = list(ref)
    for c, a in choice.items():
        seq[c] = a
    return "".join(seq)


class SampleBitsAndAtLeast(unittest.TestCase):
    def test_masks_round_trip(self):
        sb = hpc.SampleBits()
        m = sb.mask(["S3", "S1", "S2"])
        self.assertEqual(sorted(sb.samples(m)), ["S1", "S2", "S3"])
        self.assertEqual(sb.samples(sb.mask(["S2"])), ["S2"])
        self.assertEqual(sb.mask([]), 0)
        many = ["N%04d" % i for i in range(5000)]
        self.assertEqual(sorted(sb.samples(sb.mask(many))), many)

    def test_at_least_k_of_the_masks(self):
        rng = random.Random(7)
        for _ in range(500):
            masks = [rng.getrandbits(12) for _ in range(rng.randint(0, 9))]
            k = rng.randint(1, 6)
            want = 0
            for b in range(12):
                if sum(1 for m in masks if m >> b & 1) >= k:
                    want |= 1 << b
            self.assertEqual(hpc.at_least(masks, k), want)
        with self.assertRaises(ValueError):
            hpc.at_least([1], 0)


class SnpPrecheckAgainstEveryChoice(unittest.TestCase):
    """sample_can_produce over snp_columns, against every choice of every sample."""

    def check(self, case):
        sb = hpc.SampleBits()
        res = hpc.snp_columns(case["real_target"], case["ref"], case["revert"],
                              case["guide_no_pam"], case["pam"], *case["sl"],
                              case["variant_columns"], admits)
        self.assertIsNotNone(res)
        base_mm, fixes, pam_needs, blocked = res
        # the reference alignment's count is the finalizer's
        mm_ref, _ok, _t = finalizer(case["real_target"], case["ref"], case["revert"],
                                    case["guide_no_pam"], case["pam"], case["sl"])
        self.assertEqual(base_mm, mm_ref)
        carriers = {}
        for s, carried in case["samples"].items():
            for c, alts in carried.items():
                for a in alts:
                    carriers.setdefault((c, a), []).append(s)
        mask_of = lambda c, alts: sb.mask({s for a in alts for s in carriers.get((c, a), [])})  # noqa: E731
        window = 0
        for (c, a), ss in carriers.items():
            window |= sb.mask(ss)
        need = base_mm - case["bulges"] - case["allowed"]
        keep = not blocked and hpc.sample_can_produce(
            need, [mask_of(c, alts) for c, alts in fixes],
            [mask_of(c, alts) for c, alts in pam_needs], window)
        # brute force: a sample with a row (a non-empty choice of its own ALTs within the
        # budget and the PAM), and one whose union reaches them (the reference allowed)
        row, union = False, False
        for s, carried in case["samples"].items():
            if not carried:
                continue
            for choice in choices_of(case["variant_columns"], carried):
                mm, ok, _t = finalizer(case["real_target"], apply(case["ref"], choice),
                                       case["revert"], case["guide_no_pam"], case["pam"],
                                       case["sl"])
                if ok and mm - case["bulges"] <= case["allowed"]:
                    union = True
                    if choice:
                        row = True
        if row:
            self.assertTrue(keep, "the pre-check skipped a target a sample has: %r" % case)
        self.assertEqual(keep, union, case)
        return keep, row

    def test_random_windows_both_strands_3_prime(self):
        rng = random.Random(2100)
        kept = rows = 0
        for _ in range(3000):
            keep, row = self.check(random_case(rng))
            kept += keep
            rows += row
        # the cases exercise both outcomes
        self.assertGreater(rows, 100)
        self.assertGreater(3000 - kept, 100)

    def test_random_windows_5_prime_pam(self):
        rng = random.Random(12)
        for _ in range(1500):
            self.check(random_case(rng, five_prime=True))

    def test_a_pam_column_without_a_variant_blocks_every_row(self):
        # an NGG PAM whose last base is a T on the reference, no variant there
        guide = "ACGTACGTACGTACGTACGT" + "NNN"
        ref = "ACGTACGTACGTACGTACGT" + "AGT"
        res = hpc.snp_columns("A" * 23, ref, False, guide[:-3], "NGG", 0, -3, -3, None,
                              [(5, ["G"])], admits)
        self.assertTrue(res[3])
        # a variant at that base that the PAM admits makes it a PAM column to need
        res = hpc.snp_columns("A" * 23, ref, False, guide[:-3], "NGG", 0, -3, -3, None,
                              [(22, ["G", "C"])], admits)
        self.assertFalse(res[3])
        self.assertEqual(res[2], [(22, ["G"])])

    def test_a_cut_window_is_not_read(self):
        self.assertIsNone(hpc.snp_columns("A" * 23, "ACGT" * 5, False, "A" * 20, "NGG",
                                          0, -3, -3, None, [], admits))


class DegeneracyFloorAgainstEveryChoice(unittest.TestCase):
    def test_random_windows(self):
        rng = random.Random(99)
        seen_none = seen_value = 0
        for _ in range(2500):
            case = random_case(rng, five_prime=rng.random() < 0.3)
            res = hpc.snp_columns(case["real_target"], case["ref"], case["revert"],
                                  case["guide_no_pam"], case["pam"], *case["sl"],
                                  case["variant_columns"], admits)
            base_mm, fixes, pam_needs, blocked = res
            need = base_mm - case["bulges"] - case["allowed"]
            floor = hpc.degeneracy_floor(need, fixes, pam_needs, blocked)
            ref_al = "".join(aligned(case["real_target"], case["ref"], case["revert"]))
            best = None
            every = {c: alts for c, alts in case["variant_columns"]}
            for choice in choices_of(case["variant_columns"], every):
                mm, ok, t = finalizer(case["real_target"], apply(case["ref"], choice),
                                      case["revert"], case["guide_no_pam"], case["pam"],
                                      case["sl"])
                if not ok or mm - case["bulges"] > case["allowed"]:
                    continue
                d = wd.degeneracy(t, ref_al, case["guide_no_pam"], case["pam"], *case["sl"])
                best = d if best is None else min(best, d)
            self.assertEqual(floor, best, case)
            if best is None:
                seen_none += 1
            else:
                seen_value += 1
        self.assertGreater(seen_none, 50)
        self.assertGreater(seen_value, 50)


def random_indel_case(rng):
    """An overlaid indel target: aligned guide and DNA with code columns, the indel's
    carriers, and the SNP carriers of each code column's bases (or a column no genotype
    covers)."""
    sl = (0, -3)
    guide = "".join(rng.choice("ACGT") for _ in range(20)) + "NNN"
    target = list(guide[:20].lower() if rng.random() < 0.1 else guide[:20]) + list("AGG")
    for j in range(20):
        if rng.random() < 0.25:
            target[j] = rng.choice("ACGT")
    if rng.random() < 0.3:
        target[rng.randrange(1, 19)] = "-"
    for j in range(20, 23):
        if rng.random() < 0.2:
            target[j] = rng.choice("ACGT")
    codes = sorted(rng.sample([j for j in range(23) if target[j] != "-"], rng.randint(1, 5)))
    carriers = ["S%d" % i for i in range(rng.randint(1, 6))]
    others = ["T%d" % i for i in range(rng.randint(0, 3))]
    columns = []
    for j in codes:
        ref = target[j].upper()
        alts = [b for b in "ACGT" if b != ref]
        rng.shuffle(alts)
        alts = alts[:rng.randint(1, 3)]
        known = rng.random() < 0.85
        listed = {}
        for b in alts:
            listed[b] = [s for s in carriers + others if rng.random() < 0.4]
        columns.append(dict(col=j, ref=ref, known=known, code=frozenset([ref] + alts),
                            listed=listed))
        target[j] = "N"  # the code itself (any letter: the tests read the columns)
    return dict(guide=guide, target="".join(target), columns=columns, carriers=carriers,
                spacer=list(range(23))[sl[0]:sl[1]], motif={20: "N", 21: "G", 22: "G"},
                max_mm=rng.randint(0, 4), max_total=rng.choice([-1, 3, 4, 6]),
                bulge=rng.randint(0, 2))


class IndelPrecheckAgainstEveryChoice(unittest.TestCase):
    def test_random_targets(self):
        rng = random.Random(5)
        possible_seen = impossible_seen = 0
        for _ in range(3000):
            case = random_indel_case(rng)
            sb = hpc.SampleBits()
            everyone = sb.mask(case["carriers"])
            cols = [(c["col"], c["ref"], c["known"], c["code"],
                     {b: sb.mask(ss) for b, ss in c["listed"].items()})
                    for c in case["columns"]]
            got = hpc.indel_target_possible(
                case["guide"], case["target"], case["motif"], case["spacer"], case["max_mm"],
                case["max_total"], case["bulge"], cols, everyone, everyone,
                code_bases=lambda m: isc._CODE.get(m, ""))
            # brute force over every carrier's choices, counted as resolve_target counts
            snp_cols = []
            for c in case["columns"]:
                col = isc.SnpColumn(c["col"], 0, c["ref"], {}, "N", "+", known=c["known"])
                snp_cols.append(col)
            want = False
            for s in case["carriers"]:
                per = []
                for c in case["columns"]:
                    if not c["known"]:
                        per.append(sorted(c["code"] | {c["ref"]}))
                    else:
                        per.append(sorted({c["ref"]} | {b for b, ss in c["listed"].items()
                                                          if s in ss}))
                for bases in itertools.product(*per):
                    seq = isc.realize(case["target"], case["guide"], snp_cols, bases)
                    if not isc.pam_satisfied(seq, case["motif"]):
                        continue
                    mm = isc.count_mismatches(case["guide"], seq, case["spacer"])
                    if mm > case["max_mm"] or (case["max_total"] >= 0
                                               and mm + case["bulge"] > case["max_total"]):
                        continue
                    want = True
                    break
                if want:
                    break
            self.assertEqual(got, want, case)
            possible_seen += want
            impossible_seen += not want
        self.assertGreater(possible_seen, 200)
        self.assertGreater(impossible_seen, 200)

    def test_too_many_bulges_or_no_carrier(self):
        # the guide's T at column 3, a C there on the reference, the T an ALT sample 1 carries
        cols = [(3, "C", True, frozenset("CT"), {"T": 1})]
        guide, target = "ACGTACGTACGTACGTACGTNNN", "ACGYACGTACGTACGTACGTAGG"
        motif = {20: "N", 21: "G", 22: "G"}
        self.assertTrue(hpc.indel_target_possible(guide, target, motif, range(20), 0, -1, 0,
                                                  cols, 1, 1))
        self.assertFalse(hpc.indel_target_possible(guide, target, motif, range(20), 0, -1,
                                                   0, cols, 0, 0))
        self.assertFalse(hpc.indel_target_possible("ACGT-" + guide[4:], target, motif,
                                                   range(20), 4, -1, 1, cols, 1, 1,
                                                   max_dna_bulges=0))


from test_window_density import GUIDE as _G, _PostAnalysis, _site_window  # noqa: E402

_OUT_FILES = (".bestCFD.txt", ".bestmmblg.txt", ".bestCRISPR_BULGE.txt",
              ".high_variant_density_regions.bed", ".four_allele_positions.bed")


class ThroughThePostAnalysis(_PostAnalysis):
    """The real new_simple_analysis.py (classic dictionaries) on dense windows."""

    def files(self):
        out = {}
        for name in _OUT_FILES:
            path = os.path.join(self.tmp, "out" + name)
            if os.path.exists(path):
                with open(path) as fh:
                    out[name] = fh.read()
        return out

    def eleven_columns(self, carriers):
        """A window of eleven two-allele columns (2^11 sequences, over the default cap of
        10 bits) whose site is the guide + TGG but for a mismatch at column 0, 4 and 8
        (the ALTs there match the guide) and whose other columns' ALTs mismatch it.
        ``carriers``: {column: genotypes} of the ALTs."""
        site = list(_G + "TGG")
        alt = {}
        for col in (0, 4, 8):
            site[col] = {"A": "C", "C": "A", "G": "T", "T": "G"}[_G[col]]
            alt[col] = _G[col]
        for col in (2, 6, 10, 11, 12, 13, 14, 15):
            alt[col] = {"A": "G", "G": "A", "C": "T", "T": "C"}[site[col]]
        site = "".join(site)
        code = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
                frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M"}
        window = _site_window(site, {c: code[frozenset((site[c], a))] for c, a in alt.items()})
        entries = {(0, c): "%s;%s,%s;rs%d;0.1" % (carriers.get(c, ""), site[c], a, c)
                   for c, a in alt.items()}
        return site, window, entries

    def test_no_sample_reaches_the_budget_and_the_reference_row_stays(self):
        """Three mismatches, one ALT fixing each, every one carried by another sample; with
        a budget of 1 a row needs two of them on one sample: none has it, so the target is
        skipped. Within the budget of 3, its reference row is written, as without the
        pre-check."""
        site, window, entries = self.eleven_columns({0: "S1:1|0", 4: "S2:0|1", 8: "S3:1|1"})
        self.write_inputs([site], [(0, window, "+")], entries)
        _rows, said, _h, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="0")
        self.assertIn("Haplotype pre-check on chrT: 1 of 1 raw target(s) with a variant "
                      "letter skipped", said)
        on = self.files()
        _rows, said, _h, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="0",
                                           CRISPRME_HAPLOTYPE_PRECHECK="0")
        self.assertIn("Haplotype pre-check on chrT: off", said)
        self.assertEqual(self.files(), on)
        # within a budget of 3, the reference itself fits: a row of S1's or the reference
        rows, said, _h, _q = self.analyse(allowed_mms="3", CRISPRME_FAST_MODE="0")
        self.assertIn("0 of 1 raw target(s)", said)
        self.assertIn("n", [r[12] for r in rows])
        # nobody carries an ALT: skipped, and the reference row is still written
        site, window, entries = self.eleven_columns({})
        self.write_inputs([site], [(0, window, "+")], entries)
        rows, said, _h, _q = self.analyse(allowed_mms="3", CRISPRME_FAST_MODE="0")
        self.assertIn("1 of 1 raw target(s)", said)
        self.assertEqual([r[12] for r in rows], ["n"])
        on = self.files()
        self.analyse(allowed_mms="3", CRISPRME_FAST_MODE="0", CRISPRME_HAPLOTYPE_PRECHECK="0")
        self.assertEqual(self.files(), on)

    def test_one_sample_reaching_the_budget_keeps_the_target(self):
        site, window, entries = self.eleven_columns({0: "S1:1|0", 4: "S1:1|0", 8: "S3:1|1"})
        self.write_inputs([site], [(0, window, "+")], entries)
        rows, said, _h, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="0")
        self.assertIn("0 of 1 raw target(s)", said)
        self.assertIn("S1", [r[12] for r in rows])

    def test_a_window_two_alignments_share_is_one_bed_line(self):
        """Two raw targets of one window (the engine writes a mismatching base in lower
        case, so its alignments differ in case): one line, the letters resultIntegrator's
        lookup would read first (upper case sorts first); the samples carrying an ALT there."""
        site, window, entries = self.eleven_columns({0: "S1:1|0", 4: "S2:0|1"})
        low = window[:2] + window[2].lower() + window[3:]
        self.write_inputs([site], [(0, low, "+"), (0, window, "+")], entries)
        for fast in ("0", "1"):
            _rows, _said, hvdr, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE=fast)
            self.assertEqual(len(hvdr), 2, hvdr)
            self.assertEqual(hvdr[1][5:7], ["S1,S2", window])

    def test_the_k_skip_writes_what_the_cap_writes(self):
        """Population-level with --max-degeneracy 1: every row within a budget of 1 rests on
        two fixing columns or more, so the target is skipped before its genotypes are read;
        the files are those of the cap alone."""
        site, window, entries = self.eleven_columns({0: "S1:1|0", 4: "S2:0|1", 8: "S3:1|1"})
        self.write_inputs([site], [(0, window, "+")], entries)
        _rows, said, _h, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="1",
                                           CRISPRME_MAX_DEGENERACY="1")
        self.assertIn("K-skip on chrT: 1 of 1 raw target(s)", said)
        on = self.files()
        self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="1", CRISPRME_MAX_DEGENERACY="1",
                     CRISPRME_HAPLOTYPE_PRECHECK="0")
        self.assertEqual(self.files(), on)
        # at K = 2 the target is read: its worst case rests on two columns
        rows, said, _h, _q = self.analyse(allowed_mms="1", CRISPRME_FAST_MODE="1",
                                          CRISPRME_MAX_DEGENERACY="2")
        self.assertIn("K-skip on chrT: 0 of 1 raw target(s)", said)
        self.assertTrue([r for r in rows if r[12] not in ("n", "NA")])


class TheWideBudgetWarning(unittest.TestCase):
    """complete-search warns (stdout) above 1 DNA + 1 RNA bulge on a variant index that codes
    every allele (index format 2), or where none is installed and the installed CRISPRitz
    builds format 2; not at 1 + 1, nor on a format-1 index."""

    def setUp(self):
        import contextlib
        import io
        import shutil
        import tempfile
        from test_bulge_cap import _load_crisprme
        self.io, self.contextlib = io, contextlib
        self.cm = _load_crisprme()
        self.root = tempfile.mkdtemp(prefix="wide_budget_")
        self.addCleanup(shutil.rmtree, self.root, True)
        self.old_cwd = self.cm.current_working_directory
        self.cm.current_working_directory = self.root + "/"
        self.addCleanup(setattr, self.cm, "current_working_directory", self.old_cwd)
        self.vcf = os.path.join(self.root, "VCFs", "hg38_ds")
        os.makedirs(self.vcf)

    def index(self, fmt):
        d = os.path.join(self.root, "genome_library", "20-NRG_2_hg38+hg38_ds")
        os.makedirs(d, exist_ok=True)
        head = (b"CZI\xfe" + (2).to_bytes(4, sys.byteorder)) if fmt == 2 else \
            (5).to_bytes(4, sys.byteorder) + (23).to_bytes(4, sys.byteorder)
        with open(os.path.join(d, "NNNNNNNNNNNNNNNNNNNNNRG_chr22_1.bin"), "wb") as fh:
            fh.write(head + b"\0" * 64)

    def said(self, bdna, brna, engine=2):
        out = self.io.StringIO()
        old = self.cm.index_fit.engine_index_format
        self.cm.index_fit.engine_index_format = lambda: engine
        try:
            with self.contextlib.redirect_stdout(out):
                self.cm._warn_wide_bulge_budget([self.vcf], "hg38", "_", 4, bdna, brna)
        finally:
            self.cm.index_fit.engine_index_format = old
        return out.getvalue()

    def test_above_one_plus_one_on_a_format_2_index(self):
        self.index(2)
        for bdna, brna in ((2, 1), (1, 2), (2, 2), (2, 0)):
            self.assertIn("WARNING [complete-search]: 4 mismatch(es) + %d DNA + %d RNA "
                          "bulge(s) on a variant index that codes every allele" % (bdna, brna),
                          self.said(bdna, brna))
        for bdna, brna in ((1, 1), (0, 1), (1, 0), (0, 0)):
            self.assertEqual(self.said(bdna, brna), "")

    def test_not_on_a_format_1_index(self):
        self.index(1)
        self.assertEqual(self.said(2, 2), "")

    def test_no_index_yet_follows_the_installed_crispritz(self):
        self.assertIn("WARNING", self.said(2, 1, engine=2))
        self.assertEqual(self.said(2, 1, engine=1), "")


if __name__ == "__main__":
    unittest.main()

"""``scripts/haplotype_scores.py`` -- the per-haplotype scored companion.

The fixtures here are SYNTHETIC. The behaviour they pin was found on real delivered data,
but no real guide, sample identity or coordinate is committed: the repo rule is that those
never land in git, and a regression test does not need them to pin a rule.

What is pinned, and why each one matters:

  * ``row_reported_alignment_found`` fires when a site's own reported alignment is on NONE
    of the haplotypes written. On the delivering run this happened at exactly one site, of
    724 haplotypes -- the only signal distinguishing "all the haplotypes" from "the ones we
    could reconstruct". It must be a loud default, and a non-zero exit only under --strict.
  * A score is NEVER blank and the -1.0 sentinel is NEVER printed. A never-scored cell and
    an out-of-domain cell used to render alike, so a reader filtering on the score
    discarded both and called the site clean. That is the failure this file exists to stop.
  * Scores print ``%.6g``. At 3 decimal places the sensitive operating point (1.9e-4)
    rounds to 0.000 -- the threshold destroyed inside the file meant to carry it.
  * Counters are computed BEFORE any cap, so a capped file still reports true totals.
  * An unreadable budget, or a run with no haplotype source, REFUSES to write rather than
    emitting a complete-looking file that means "never enumerated".
"""

import io
import os
import shutil
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.normpath(os.path.join(_HERE, "..", "scripts"))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import haplotype_scores as hs  # noqa: E402


SPACER = "AAAACCCCGGGGTTTTAAAANNN"          # synthetic, not any real guide
REPORTED = "AAAACCCCGGGGTTTTAAAATGG"        # what the site's own row reports
HAP_ONE = "AAAACCCCGGGGTTTTAAACTGG"         # a haplotype that is NOT the reported one
HAP_TWO = "AAAACCCCGGGGTTTTAAAGTGG"

PARAMS = (
    "Genome_selected\thg38\n"
    "Pam\tNRR\n"
    "Max_bulges\t1\n"
    "Mismatches\t6\n"
    "DNA\t1\n"
    "RNA\t1\n"
    "Max_total_edits\t8\n"
    "Scorer\tcrispr-bulge\n"
)

_BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _header():
    """The REAL per-block schema, read off a real run rather than assumed.

    Measured on ml008, nrr281/N3 observed_haplotypes_raw.tsv:
        CFD_score_(highest_CFD)                  CFD_score_(fewest_mm+b)
        CRISPR_BULGE_score_(highest_CRISPR_BULGE)
    and NO CFD_score_(highest_CRISPR_BULGE). So NO observation ever carries both a CFD
    and a CRISPR-Bulge number: each block has exactly one score column.

    This fixture has now been wrong twice. First it wrote one start into all three blocks
    and invented a CRISPR-Bulge column on each, hiding five defects. Then it was
    "corrected" to put CFD on all three -- still a column no writer emits, which made the
    composite scored(cb-model-refused) status it was green on UNREACHABLE in production.
    A fixture that is not the real shape is not coverage; it is a second implementation
    nobody runs.
    """
    cols = ["Spacer+PAM", "Chromosome"]
    for b in _BLOCKS:
        cols += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                 "Aligned_spacer+PAM_(%s)" % b,
                 "Aligned_protospacer+PAM_REF_(%s)" % b,
                 "Aligned_protospacer+PAM_ALT_(%s)" % b,
                 "REF/ALT_origin_(%s)" % b]
        cols += (["CRISPR_BULGE_score_(%s)" % b] if b == "highest_CRISPR_BULGE"
                 else ["CFD_score_(%s)" % b])
        cols += ["Variant_info_genome_(%s)" % b, "Variant_MAF_(%s)" % b,
                 "Variant_samples_(%s)" % b]
    return cols


def _row(chrom, start, alt, cfd, cb, variants, samples="s1,s2", origin="alt",
         cb_start=None):
    """One results row. ``cb_start`` places the CRISPR-Bulge block at a DIFFERENT
    coordinate, which about 2% of real rows do (measured) and what the site-identity
    anchor has to survive."""
    out = [SPACER, chrom]
    for b in _BLOCKS:
        st = cb_start if (b == "highest_CRISPR_BULGE" and cb_start is not None) else start
        out += [str(st), "+", SPACER, REPORTED, alt, origin,
                str(cb) if b == "highest_CRISPR_BULGE" else str(cfd)]
        out += [variants, "0.01", samples]
    return out


class HaplotypeScoresMixin(object):
    """setUp and helpers ONLY -- never a TestCase.

    This used to be a TestCase carrying its own tests, and every class that reused the
    helpers inherited and re-ran them. Five subclasses turned ~35 real tests into 109
    collected, multiplying runtime and inflating the count.
    """

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="hapscores_")
        self.addCleanup(shutil.rmtree, self.dir, True)
        with io.open(os.path.join(self.dir, ".Params.txt"), "w", encoding="utf-8") as fh:
            fh.write(PARAMS)

    def _write_source(self, rows, name="x_observed_haplotypes_raw.tsv"):
        path = os.path.join(self.dir, name)
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        return path

    @staticmethod
    def _parse(text):
        """Body rows as dicts keyed by the header. NEVER index a column by position:
        adding rank_score_axis silently shifted every later index."""
        lines = [l for l in text.splitlines() if l and not l.startswith("#")]
        hdr = lines[0].split("\t")
        return [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]

    def _run(self, rows, **kw):
        self._write_source(rows)
        out = os.path.join(self.dir, "out.tsv")
        msgs = []
        res = hs.build(self.dir, out, log=msgs.append, **kw)
        with io.open(out, "r", encoding="utf-8") as fh:
            text = fh.read()
        return res, text, msgs

    def _targets(self, rows):
        """A curated table naming each site's OWN reported haplotype.

        Reconciliation keys on the VARIANT SET, not the alignment string: the curated
        alignment is the representative's own, while a per-haplotype row carries the
        haplotype's bases in the candidate's alignment shape. The curated notation is
        ``chr;pos;REF;ALT (rsID)`` against the per-row ``chr_pos_REF_ALT``, so both are
        reduced to a canonical key.
        """
        path = os.path.join(self.dir, "targets.tsv")
        guides = [r[4] if len(r) > 4 else None for r in rows]
        named = any(g is not None for g in guides)
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("Chromosome\tPosition\tVariant\tREF/ALT_origin\t"
                     "Aligned_protospacer+PAM\tCRISPR-Bulge"
                     + ("\tSpacer+PAM\n" if named else "\n"))
            for row, guide in zip(rows, guides):
                chrom, pos, variant, origin = row[:4]
                fh.write("%s\t%s\t%s\t%s\t%s\t0.01%s\n"
                         % (chrom, pos, variant, origin, REPORTED,
                            ("\t%s" % (guide or "")) if named else ""))
        return path


class HaplotypeScoresFixture(HaplotypeScoresMixin, unittest.TestCase):
    """The reconciliation signal and the file's invariants."""

    def test_a_site_whose_reported_alignment_is_on_no_haplotype_is_flagged(self):
        """The signal separating 'all the haplotypes' from 'the ones we could get'.

        Reproduces, synthetically, the one site in the delivered data where the published
        ALT protospacer appeared on none of the 724 haplotypes written.
        """
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.4, 0.003, "chr1_100_A_G")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;T (rs9)", "ALT")])  # on NEITHER
        (_n, _s, unreconciled), text, msgs = self._run(rows, targets=tgt)
        self.assertEqual(unreconciled, ["chr1:100:+:%s" % SPACER],
                         "the site must be reported as unreconciled")
        self.assertIn("row_reported_alignment_found", text)
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("row_reported_alignment_found")
        self.assertTrue(all(r[i] == "no" for r in data), "every row of the site says no")
        self.assertTrue(any("WARNING" in m for m in msgs), "it must warn loudly")

    def test_a_site_whose_reported_alignment_is_present_reconciles(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])  # IS one of them
        (_n, _s, unreconciled), text, _msgs = self._run(rows, targets=tgt)
        self.assertEqual(unreconciled, [], "it reconciles")
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("row_reported_alignment_found")
        self.assertTrue(all(r[i] == "yes" for r in data))

    def test_a_systematically_miskeyed_comparison_is_called_out_not_reported(self):
        """A signal that fires on nearly every site is a wrong key, not a finding.

        Measured on the delivering run: the curated Aligned_protospacer+PAM equals the
        per-haplotype ALT alignment on 0 of 244 rows, so a naive comparison flagged 106 of
        114 sites. Emitting that as 106 failures would be noise a reader cannot tell from
        signal, so it degrades to 'unknown' and says why.
        """
        rows, tgt_rows = [], []
        for i in range(12):
            pos = 1000 + i
            rows.append(_row("chr1", pos, HAP_ONE, 0.5, 0.004, "chr1_%d_A_C" % pos))
            tgt_rows.append(("chr1", str(pos), "chr1;%d;A;T (rs9)" % pos, "ALT"))
        tgt = self._targets(tgt_rows)
        (_n, _s, unreconciled), text, msgs = self._run(rows, targets=tgt)
        self.assertEqual(unreconciled, [], "a mis-keyed comparison reports no failures")
        self.assertTrue(any("LOOKS like a mis-keyed comparison" in m for m in msgs),
                        "it must say the key looks wrong")
        self.assertTrue(any("Variant column" in m for m in msgs),
                        "it must name the column actually compared, not the alignment")
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("row_reported_alignment_found")
        self.assertTrue(all(r[i] == "unknown(key-looks-wrong)" for r in data), 
                        set(r[i] for r in data))

    def test_one_genuine_failure_is_not_suppressed_by_that_guard(self):
        """The guard needs evidence: below its floor, a real failure still reports."""
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.4, 0.003, "chr1_100_A_G")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;T (rs9)", "ALT")])
        (_n, _s, unreconciled), _text, msgs = self._run(rows, targets=tgt)
        self.assertEqual(len(unreconciled), 1, "a lone genuine failure must survive")
        self.assertFalse(any("mis-keyed" in m for m in msgs))

    def test_without_targets_it_says_unknown_never_yes(self):
        """Silence must not be reported as success."""
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("row_reported_alignment_found")
        self.assertTrue(all(r[i] == "unknown" for r in data))

    def test_it_warns_loudly_and_still_exits_zero_by_default(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.4, 0.003, "chr1_100_A_G")]
        self._write_source(rows)
        tgt = self._targets([("chr1", "100", "chr1;100;A;T (rs9)", "ALT")])
        out = os.path.join(self.dir, "out.tsv")
        rc = hs.main(["--results", self.dir, "--out", out, "--targets", tgt])
        self.assertEqual(rc, 0, "default is warn, never fail")

    def test_strict_turns_the_same_condition_into_a_non_zero_exit(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.4, 0.003, "chr1_100_A_G")]
        self._write_source(rows)
        tgt = self._targets([("chr1", "100", "chr1;100;A;T (rs9)", "ALT")])
        out = os.path.join(self.dir, "out.tsv")
        rc = hs.main(["--results", self.dir, "--out", out, "--targets", tgt, "--strict"])
        self.assertEqual(rc, 1, "--strict opts in to failing")

    # -- never a blank, never -1.0 ------------------------------------------ #

    def test_the_null_sentinel_is_never_printed_as_a_number(self):
        """-1.0 means 'the model refused', not 'scored very low'."""
        self.assertEqual(hs.fmt_score(-1.0), "NA")
        self.assertEqual(hs.fmt_score(hs.CB_NULL), "NA")
        self.assertEqual(hs.fmt_score(None), "NA")

    def test_every_row_carries_a_score_status(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, -1.0, "chr1_100_A_C"),
                _row("chr2", 200, HAP_TWO, "NA", "NA", "chr2_200_A_G")]
        _res, text, _msgs = self._run(rows)
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("score_status")
        self.assertTrue(data, "fixture produced no rows")
        for row in data:
            self.assertTrue(row[i].strip(), "score_status must never be blank")

    def test_a_score_keeps_six_significant_digits(self):
        """At 3 dp the 1.9e-4 operating point rounds to 0.000."""
        self.assertEqual(hs.fmt_score(0.00019), "0.00019")
        self.assertNotEqual(hs.fmt_score(0.00019), "0.000")
        self.assertEqual(hs.fmt_score(0.001933604129590094), "0.0019336")

    # -- count before you cap ------------------------------------------------ #

    def test_a_cap_hides_a_sequence_but_never_a_number(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.9, 0.5, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.8, 0.4, "chr1_100_A_G"),
                _row("chr1", 100, "AAAACCCCGGGGTTTTAAATTGG", 0.7, 0.3, "chr1_100_A_T")]
        _res, text, _msgs = self._run(rows, cap=1)
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        n_i = head.index("n_haplotypes_at_site")
        totals = {r[n_i] for r in data}
        self.assertEqual(len(data), 1, "the cap limits rows written")
        self.assertEqual(totals, {"3"}, "counters are the TRUE totals, not the capped ones")

    # -- refuse rather than mislead ------------------------------------------ #

    def test_a_run_with_no_haplotype_source_refuses_to_write(self):
        """An empty file would read as 'nothing above threshold'. It must not."""
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.build(self.dir, os.path.join(self.dir, "out.tsv"))
        self.assertIn("population-level fast mode", str(cm.exception))

    def test_an_unreadable_budget_refuses_rather_than_guessing_one(self):
        os.remove(os.path.join(self.dir, ".Params.txt"))
        self._write_source([_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")])
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.build(self.dir, os.path.join(self.dir, "out.tsv"))
        self.assertIn("never passed as a flag", str(cm.exception))

    def test_the_budget_is_quoted_inside_the_status_not_invented(self):
        budget = hs.read_budget(self.dir)
        self.assertEqual(hs.budget_text(budget), "6mm/1DNA/1RNA/8total")

    # -- the header carries its own preconditions ---------------------------- #

    def test_the_header_states_the_threshold_its_ci_and_the_frame_caveat(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        for needle in ("score_source=", "budget=6mm/1DNA/1RNA/8total",
                       "[5.8e-06, 3.1e-03]", "candidate_frame", "COUNT BEFORE CAP",
                       "population-level fast mode", "targets-scoring"):
            self.assertIn(needle, text, "the header must carry %r" % needle)

    def test_carrier_identities_can_be_withheld(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C", samples="a,b,c")]
        _res, text, _msgs = self._run(rows, anonymise=True)
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        c_i, n_i = head.index("carriers"), head.index("n_carriers")
        self.assertEqual(data[0][c_i], "", "identities withheld under --anonymise")
        self.assertEqual(data[0][n_i], "3", "the count is still reported")




class CarrierReporting(unittest.TestCase):
    """Carriers are a reagent list: 1000G and HGDP lines can be bought, so a carrier id is
    something a reader acts on. Zygosity is REPORTED so they can plan the quantification --
    allele-specific quantification handles a heterozygote, so it is never a filter, a
    warning, or a reason to prefer one line over another."""

    def test_a_sample_on_both_slots_is_homozygous(self):
        n, hom, het = hs.carrier_breakdown("S1:0,S1:1,S2:0")[:3]
        self.assertEqual((n, hom, het), (2, "1", "1"))

    def test_zygosity_is_counted_not_ranked(self):
        """A haplotype carried only heterozygously is reported, never flagged."""
        n, hom, het = hs.carrier_breakdown("S1:0,S2:1,S3:0")[:3]
        self.assertEqual((n, hom, het), (3, "0", "3"))

    def test_a_source_without_slots_says_unknown_not_zero(self):
        """THE BUG THIS PINS. Only the SAMPLE:slot form records which chromosome copy
        carries the haplotype; the per-haplotype results columns give bare sample names.
        Printing 0 homozygous there reads as 'none are', when the truth is 'this source
        does not say' -- a not-measured cell wearing the clothes of a measured zero, which
        is the exact failure this file exists to prevent. Seen on real data: 690 rows all
        reporting 0 homozygous."""
        n, hom, het, panels = hs.carrier_breakdown("S1,S2,S3")[:4]
        self.assertEqual(n, 3, "the carriers are still counted")
        self.assertEqual(hom, "unknown")
        self.assertEqual(het, "unknown")
        self.assertNotEqual(hom, "0", "never a measured-looking zero")
        self.assertTrue(panels, "the panel breakdown still works without slots")

    def test_panels_are_named_so_nothing_implies_orderability(self):
        """Not every carrier id is a standard accession; say which collection it is from
        rather than implying every one can be bought."""
        panels = hs.carrier_breakdown(
            "NA19144:0,HG01052:1,HGDP00197:0,LP6005441-DNA_A11:1,SS6004474:0")[3]
        self.assertIn("1000G:2", panels)
        self.assertIn("HGDP:1", panels)
        self.assertIn("other:2", panels)

    def test_no_carriers_is_zero_not_a_crash(self):
        self.assertEqual(hs.carrier_breakdown("")[:4], (0, "0", "0", "-"))
        self.assertEqual(hs.carrier_breakdown(None)[:4], (0, "0", "0", "-"))

    # -- ancestry ----------------------------------------------------------- #

    def test_ancestry_is_summarised_per_carrier_set(self):
        """The run knows the ancestry of everyone it searched; no delivered file says so."""
        anc = {"S1": {"population": "Yoruba", "superpopulation": "AFR", "sex": "male"},
               "S2": {"population": "Tuscan", "superpopulation": "EUR", "sex": "female"},
               "S3": {"population": "Yoruba", "superpopulation": "AFR", "sex": "female"}}
        out = hs.carrier_breakdown("S1:0,S2:0,S3:1", ancestry=anc)
        self.assertIn("Yoruba:2", out[4])
        self.assertIn("Tuscan:1", out[4])
        self.assertIn("AFR:2", out[5])
        self.assertIn("EUR:1", out[5])

    def test_no_panel_file_says_unknown_not_an_undiverse_panel(self):
        """An absent ancestry file must not read as 'this panel has no diversity'."""
        out = hs.carrier_breakdown("S1:0,S2:1", ancestry={})
        self.assertEqual(out[4], "unknown")
        self.assertEqual(out[5], "unknown")

    def test_a_carrier_missing_from_the_panel_is_unlabelled_not_dropped(self):
        anc = {"S1": {"population": "Yoruba", "superpopulation": "AFR", "sex": ""}}
        out = hs.carrier_breakdown("S1:0,MYSTERY:1", ancestry=anc)
        self.assertIn("unlabelled:1", out[4])
        self.assertIn("Yoruba:1", out[4])


class TestSiteIdentityAcrossDivergentBlocks(HaplotypeScoresMixin, unittest.TestCase):
    """The three alignment blocks place the same physical off-target at DIFFERENT starts.

    Measured on real runs: 2.1% of 73,994 rows genome-wide per-sample, 2.6% of 31,502 on a
    three-chromosome run, almost all a 1-2 bp shift. Before the anchor, the site key was
    built from
    whichever block was being read, so one site split into up to three and every "per-site"
    roll-up was really per-FRAME. With --targets it was worse than cosmetic: the filter
    compared the per-block start against a curated Position that is always the highest-CFD
    start, so the CRISPR-Bulge observation was DROPPED -- and the only axis carrying a
    CRISPR-Bulge score disappeared, rendering a scored haplotype as unscored.

    The old fixture could not catch any of this: it wrote the same start into all three
    blocks, so divergence was never exercised.
    """

    def test_divergent_block_starts_are_one_site_not_three(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        (_n, n_sites, _u), text, _msgs = self._run(rows)
        self.assertEqual(n_sites, 1,
                         "one physical off-target must be one site, whatever each block "
                         "calls its start:\n" + text)
        self.assertEqual(len(set(r["site_key"] for r in self._parse(text))), 1)

    def test_each_row_still_reports_its_own_blocks_coordinate(self):
        """Anchoring site IDENTITY must not rewrite the coordinate a block actually used."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        _res, text, _msgs = self._run(rows)
        body = self._parse(text)
        starts = set(r["site_start"] for r in body)
        self.assertEqual(starts, {"100", "97"},
                         "each observation reports its own block's start: %s" % starts)

    def test_targets_does_not_drop_the_block_that_sits_elsewhere(self):
        """The false clean. --targets names the highest-CFD start; the CRISPR-Bulge block
        sits 3 bp away and is the ONLY axis carrying a CRISPR-Bulge score. If the filter
        keys on the per-block start, that score vanishes and a haplotype thousands of
        times over the operating point renders as never scored."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        _res, text, _msgs = self._run(rows, targets=tgt)
        body = self._parse(text)
        self.assertTrue(body, "the CRISPR-Bulge observation was dropped entirely")
        cb = set(r["crispr_bulge"] for r in body)
        self.assertIn("0.9", cb,
                      "the CRISPR-Bulge score was lost with the block it lives on: %s" % cb)
        self.assertEqual(set(r["site_max_crispr_bulge"] for r in body), {"0.9"},
                         "site_max_crispr_bulge must see the score too")

    def test_reconciliation_survives_the_divergence(self):
        """The reported haplotype is found even though one block sits at another start."""
        rows = [_row("chr1", 100, REPORTED, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        (_n, _s, unreconciled), text, _msgs = self._run(rows, targets=tgt)
        self.assertEqual(unreconciled, [], text)
        self.assertIn("\tyes", text)


class TestTargetsKeyIsPerGuide(HaplotypeScoresMixin, unittest.TestCase):
    """Two guides can report at the SAME locus. Keying --targets on (chrom, position)
    alone let the later curated row overwrite the earlier one, so one guide was quietly
    reconciled against the other guide's haplotype."""

    def test_two_guides_at_one_locus_are_kept_apart(self):
        rep = hs.read_reported(self._targets([
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT", "GUIDE_ONE"),
            ("chr1", "100", "chr1;100;A;T (rs9)", "ALT", "GUIDE_TWO"),
        ]))
        self.assertEqual(len(rep), 2, rep)
        self.assertNotEqual(rep[("chr1", "100", "GUIDE_ONE")].variants,
                            rep[("chr1", "100", "GUIDE_TWO")].variants)

    def test_lookup_prefers_the_matching_guide(self):
        rep = hs.read_reported(self._targets([
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT", "GUIDE_ONE"),
            ("chr1", "100", "chr1;100;A;T (rs9)", "ALT", "GUIDE_TWO"),
        ]))
        got = hs.lookup_reported(rep, "chr1", "100", "GUIDE_TWO")
        self.assertEqual(got.variants, ("chr1:100:A:T",))

    def test_a_table_with_no_guide_column_still_works(self):
        """Every curated table shipped so far has no guide column. It must keep working."""
        rep = hs.read_reported(self._targets(
            [("chr1", "100", "chr1;100;A;C (rs1)", "ALT")]))
        got = hs.lookup_reported(rep, "chr1", "100", SPACER)
        self.assertIsNotNone(got)
        self.assertEqual(got.variants, ("chr1:100:A:C",))

    def test_an_incompatible_duplicate_is_refused_not_guessed(self):
        """Same locus, same (or absent) guide, two DIFFERENT reported haplotypes. There is
        no right answer, so the answer is 'unknown', loudly -- never the last row read."""
        msgs = []
        rep = hs.read_reported(self._targets([
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
            ("chr1", "100", "chr1;100;A;T (rs9)", "ALT"),
        ]), say=msgs.append)
        self.assertIs(rep[("chr1", "100", "")], hs._AMBIGUOUS)
        self.assertTrue(any("named more than once" in m for m in msgs), msgs)

    def test_an_ambiguous_locus_never_claims_yes_or_no(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets([
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
            ("chr1", "100", "chr1;100;A;T (rs9)", "ALT"),
        ])
        (_n, _s, unreconciled), text, _msgs = self._run(rows, targets=tgt)
        self.assertIn("unknown(ambiguous-in-targets)", text)
        self.assertEqual(unreconciled, [],
                         "an ambiguous locus is not a reconciliation FAILURE")

    def test_an_identical_duplicate_row_is_not_ambiguous(self):
        """Curated tables repeat a site across guides/rows; identical repeats are fine."""
        msgs = []
        rep = hs.read_reported(self._targets([
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
            ("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
        ]), say=msgs.append)
        self.assertIsNot(rep[("chr1", "100", "")], hs._AMBIGUOUS)
        self.assertFalse([m for m in msgs if "named more than once" in m], msgs)


class TestTheFixtureMatchesTheRealSchema(unittest.TestCase):
    """The fixture is itself under test. It once invented a CRISPR-Bulge column on all
    three blocks -- a shape the pipeline never writes -- so the suite was green against a
    file that does not exist."""

    def test_crispr_bulge_is_only_on_the_crispr_bulge_block(self):
        cols = _header()
        cb = [c for c in cols if c.startswith("CRISPR_BULGE_score_")]
        self.assertEqual(cb, ["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"], cb)

    def test_cfd_is_not_on_the_crispr_bulge_block(self):
        """Measured on a real run: CFD is on highest_CFD and fewest_mm+b ONLY."""
        cols = _header()
        self.assertEqual(
            sorted(c for c in cols if c.startswith("CFD_score_")),
            sorted(["CFD_score_(highest_CFD)", "CFD_score_(fewest_mm+b)"]))

    def test_no_observation_can_carry_both_scores(self):
        """The consequence that matters: each block has exactly one score column, so a
        row with a CFD number AND a CRISPR-Bulge number is not a shape that exists."""
        cols = _header()
        for b in _BLOCKS:
            here = [c for c in cols if c.endswith("(%s)" % b) and "_score_" in c]
            self.assertEqual(len(here), 1, (b, here))


class TestAnEmptySourceNeverWins(HaplotypeScoresMixin, unittest.TestCase):
    """Found on real data, not by reading the code.

    The raw sidecar is gated on a Tier-1 dict-less genotype store. A run without that
    store writes the header and NO rows. probe_source preferred it unconditionally, so a
    run whose integrated table held 16,208 ALT rows with real carriers produced a 0-row
    companion and exited 0 -- a false clean, from the file written to stop false cleans.
    """

    def _write(self, name, rows):
        path = os.path.join(self.dir, name)
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        return path

    def test_a_header_only_raw_falls_through_to_the_alt_table(self):
        self._write("x_observed_haplotypes_raw.tsv", [])                      # header only
        self._write("x_all_results_with_alternative_alignments.tsv",
                    [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])
        path, kind = hs.probe_source(self.dir)
        self.assertEqual(kind, "all_results_with_alternative_alignments", path)

    def test_the_fallthrough_actually_produces_rows(self):
        self._write("x_observed_haplotypes_raw.tsv", [])
        self._write("x_all_results_with_alternative_alignments.tsv",
                    [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])
        out = os.path.join(self.dir, "out.tsv")
        (n_rows, n_sites, _u) = hs.build(self.dir, out, log=lambda m: None)
        self.assertTrue(n_rows > 0 and n_sites == 1, (n_rows, n_sites))

    def test_a_populated_raw_is_still_preferred(self):
        """The fall-through must not demote the richer source when it HAS rows."""
        self._write("x_observed_haplotypes_raw.tsv",
                    [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])
        self._write("x_all_results_with_alternative_alignments.tsv",
                    [_row("chr1", 100, HAP_TWO, 0.4, 0.9, "chr1_100_A_G")])
        _path, kind = hs.probe_source(self.dir)
        self.assertEqual(kind, "observed_haplotypes_raw")

    def test_all_sources_empty_refuses_rather_than_writing_a_clean_looking_file(self):
        self._write("x_observed_haplotypes_raw.tsv", [])
        self._write("x_all_results_with_alternative_alignments.tsv", [])
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.probe_source(self.dir)
        self.assertIn("header-only", str(cm.exception))

    def test_the_refusal_is_a_non_zero_exit(self):
        self._write("x_observed_haplotypes_raw.tsv", [])
        rc = hs.main(["--results", self.dir,
                      "--out", os.path.join(self.dir, "out.tsv")])
        self.assertNotEqual(rc, 0)
        self.assertFalse(os.path.exists(os.path.join(self.dir, "out.tsv")),
                         "a refusal must leave no file to be mistaken for an answer")


class TestTargetsGapsAreNamed(HaplotypeScoresMixin, unittest.TestCase):
    """A --targets site that produced nothing must be SAID, not silently absent.

    Measured on a real NRR run: of 200 curated top-CFD ALT sites, only 97 appear at that
    exact start in the haplotype source (14 more at another block's start, 47 within 5 bp,
    42 absent). Returning 97 rows with no further word lets a reader believe the other
    103 were examined and found clean.
    """

    def test_a_requested_site_that_produced_nothing_is_reported(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
                             ("chr9", "777", "chr9;777;G;T (rs7)", "ALT")])   # absent
        _res, _text, msgs = self._run(rows, targets=tgt)
        gap = [m for m in msgs if "produced" in m and "NO rows" in m]
        self.assertTrue(gap, msgs)
        self.assertIn("chr9:777", gap[0])
        self.assertIn("1 of the 2 site(s)", gap[0])

    def test_nothing_is_said_when_every_requested_site_produced_rows(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertFalse([m for m in msgs if "NO rows" in m], msgs)

    def test_a_divergent_block_is_not_miscounted_as_a_gap(self):
        """The gap is measured against the ANCHOR, not the block that created the site.
        Comparing against the creating block's own start would invent a gap on the ~2% of
        rows whose blocks diverge."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertFalse([m for m in msgs if "NO rows" in m], msgs)

    def test_the_gap_message_refuses_to_call_it_clean(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
                             ("chr9", "777", "chr9;777;G;T (rs7)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        gap = [m for m in msgs if "NO rows" in m][0]
        self.assertIn("not enumerated", gap)
        self.assertIn("NOT fuzzy-matched", gap)


class TestACapNeverDiscardsTheScoredAxis(HaplotypeScoresMixin, unittest.TestCase):
    """CFD and CRISPR-Bulge are different scales and were ordered in one sequence.

    CFD runs 0..1; the CRISPR-Bulge operating point is 1.9e-4. Only the CRISPR-Bulge
    block carries a CB column, so every CFD-only row sorted above every CB row -- and
    cap=1 on a site with CFD 0.6 and CB 0.004 (20x the operating point) wrote the CFD
    row and DISCARDED the only haplotype above the threshold, on a line that still said
    n_above_threshold=1.
    """

    def test_the_above_threshold_row_survives_a_cap_of_one(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.6, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows, cap=1)
        body = self._parse(text)
        self.assertEqual(len(body), 1)
        self.assertEqual(body[0]["crispr_bulge"], "0.004",
                         "the cap kept the CFD row and dropped the scored one")
        self.assertEqual(body[0]["above_threshold"], "yes")

    def test_the_file_never_claims_more_above_threshold_than_it_shows(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.6, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows, cap=1)
        body = self._parse(text)
        shown = len([r for r in body if r["above_threshold"] == "yes"])
        self.assertEqual(int(body[0]["n_above_threshold"]), shown)

    def test_rank_says_which_score_it_used(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.6, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        axes = set(r["rank_score_axis"] for r in self._parse(text))
        self.assertTrue(axes <= {"crispr_bulge", "cfd", "none"}, axes)
        self.assertIn("crispr_bulge", axes)

    def test_hidden_is_counted_in_haplotypes_not_rows(self):
        """One haplotype across three frames is ONE haplotype. The message said 2 were
        not written while the rows said n_haplotypes_at_site=1."""
        rows = [_row("chr1", 100, HAP_ONE, 0.6, 0.004, "chr1_100_A_C")]
        _res, text, msgs = self._run(rows, cap=1)
        note = [m for m in msgs if "not written" in m]
        self.assertTrue(note, msgs)
        self.assertIn("0 haplotype(s)", note[0])
        self.assertEqual(self._parse(text)[0]["n_haplotypes_at_site"], "1")


class TestARefusedCfdIsNotCalledUnalignable(HaplotypeScoresMixin, unittest.TestCase):
    """new_simple_analysis.py writes -1.000 into the CFD column when it computes no
    score. CFD is the one column on all three blocks, so on the two with no CB column a
    refused CFD left both None and typed the row 'needs-realignment' -- saying the frame
    could not hold it when the scorer reached it and declined."""

    def test_a_minus_one_cfd_types_as_refused_not_needs_realignment(self):
        rows = [_row("chr1", 100, HAP_ONE, -1.0, -1.0, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        statuses = set(r["score_status"] for r in self._parse(text))
        self.assertTrue(all("model-refused" in st for st in statuses), statuses)
        self.assertFalse([st for st in statuses if "needs-realignment" in st], statuses)

    def test_the_sentinel_is_still_never_printed(self):
        rows = [_row("chr1", 100, HAP_ONE, -1.0, -1.0, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        for r in self._parse(text):
            self.assertNotIn("-1.0", (r["cfd"], r["crispr_bulge"]))


class TestRefOnlySitesAreNotCalledAbsent(HaplotypeScoresMixin, unittest.TestCase):
    """A curated table is mostly REF-origin (115 of 187 sites on the delivering run).
    Those sites ARE in the source; they just carry no ALT haplotype. The gap warning
    skipped ref rows before noting the anchor, so it called them 'not enumerated here'."""

    def test_a_ref_only_site_is_reported_as_present_with_no_alt(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", origin="ref")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "REF")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        present = [m for m in msgs if "carry only reference rows" in m]
        self.assertTrue(present, msgs)
        self.assertFalse([m for m in msgs if "NO rows" in m],
                         "a ref-only site must not be called absent")

    def test_a_genuinely_absent_site_is_still_called_absent(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", origin="ref")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "REF"),
                             ("chr9", "777", "chr9;777;G;T (rs7)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        gap = [m for m in msgs if "NO rows" in m]
        self.assertTrue(gap, msgs)
        self.assertIn("chr9:777", gap[0])
        self.assertNotIn("chr1:100", gap[0])


class TestAGuideMismatchIsNotSilent(HaplotypeScoresMixin, unittest.TestCase):
    """A curated 'crRNA'/'Guide' column holds the 20-mer WITHOUT the PAM while the
    source's Spacer+PAM includes it. The locus passes the coarse filter, no guide-keyed
    entry matches, and reconciliation goes 'unknown' on every row with nothing said."""

    def _targets_named(self, guide):
        path = os.path.join(self.dir, "targets_g.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("Chromosome\tPosition\tVariant\tREF/ALT_origin\tcrRNA\n")
            fh.write("chr1\t100\tchr1;100;A;C (rs1)\tALT\t%s\n" % guide)
        return path

    def test_a_guide_that_does_not_match_the_source_is_warned_about(self):
        rows = [_row("chr1", 100, REPORTED, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets_named(SPACER[:20])          # the 20-mer, no PAM
        _res, text, msgs = self._run(rows, targets=tgt)
        warn = [m for m in msgs if "DIFFERENT guide" in m]
        self.assertTrue(warn, msgs)
        self.assertIn("WITHOUT the PAM", warn[0])
        self.assertEqual(set(r["row_reported_alignment_found"]
                             for r in self._parse(text)), {"unknown"})

    def test_a_matching_guide_says_nothing_and_reconciles(self):
        rows = [_row("chr1", 100, REPORTED, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets_named(SPACER)
        _res, text, msgs = self._run(rows, targets=tgt)
        self.assertFalse([m for m in msgs if "DIFFERENT guide" in m], msgs)
        self.assertEqual(set(r["row_reported_alignment_found"]
                             for r in self._parse(text)), {"yes"})


class TestAnUnreadableSourceIsNotCalledEmpty(HaplotypeScoresMixin, unittest.TestCase):

    def test_a_corrupt_archive_raises_instead_of_reading_as_header_only(self):
        path = os.path.join(self.dir, "x_observed_haplotypes_raw.zip")
        with io.open(path, "wb") as fh:
            fh.write(b"this is not a zip file")
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.probe_source(self.dir)
        msg = str(cm.exception)
        self.assertIn("cannot be read", msg)
        self.assertNotIn("header-only", msg)


class TestOneAxisRefusingIsNotTheRowRefusing(HaplotypeScoresMixin, unittest.TestCase):
    """Introduced by the fix for the CFD sentinel, caught on the next review pass.

    Marking the row out-of-domain whenever EITHER axis held the -1 sentinel made a
    haplotype with a real CRISPR-Bulge score of 0.004 -- twenty times the operating
    point -- count as never scored. The same row was above_threshold=yes AND part of
    n_not_scored, and a reader keeping only score_status=scored discarded the one row
    above the threshold.
    """

    def test_a_refused_cfd_and_a_real_cb_are_separate_rows_with_separate_answers(self):
        """In the REAL schema no observation carries both scores: each block has exactly
        one score column. So "CFD refused while CB scored" is two ROWS, and each must
        tell the truth about itself."""
        rows = [_row("chr1", 100, HAP_ONE, -1.0, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        body = self._parse(text)
        cb = [r for r in body if "CRISPR_BULGE" in r["haplotype_id"]]
        cfd = [r for r in body if "CRISPR_BULGE" not in r["haplotype_id"]]
        self.assertEqual(len(cb), 1)
        self.assertEqual(cb[0]["score_status"], "scored")
        self.assertEqual(cb[0]["above_threshold"], "yes")
        for r in cfd:
            self.assertIn("model-refused", r["score_status"],
                          "a -1 CFD must not read as never measured")

    def test_a_row_is_never_both_above_threshold_and_never_scored(self):
        rows = [_row("chr1", 100, HAP_ONE, -1.0, 0.004, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        for r in self._parse(text):
            if int(r["n_above_threshold"]) > 0:
                self.assertEqual(r["n_not_scored"], "0",
                                 "the one haplotype is counted as both")

    def test_nothing_surviving_is_still_out_of_domain(self):
        rows = [_row("chr1", 100, HAP_ONE, -1.0, -1.0, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        self.assertTrue(all("model-refused" in r["score_status"]
                            for r in self._parse(text)))


class TestAnUnreadableCandidateDoesNotKillTheFallback(HaplotypeScoresMixin,
                                                      unittest.TestCase):
    """Also introduced by the previous fix. Making has_data_row raise aborted
    probe_source before the alt table was looked at, so a corrupt raw sidecar beside a
    complete alt table refused a run that could be answered."""

    def _write(self, name, rows):
        path = os.path.join(self.dir, name)
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        return path

    def test_a_corrupt_raw_falls_through_to_a_good_alt_table(self):
        with io.open(os.path.join(self.dir, "x_observed_haplotypes_raw.zip"), "wb") as fh:
            fh.write(b"this is not a zip file")
        self._write("x_all_results_with_alternative_alignments.tsv",
                    [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])
        path, kind = hs.probe_source(self.dir)
        self.assertEqual(kind, "all_results_with_alternative_alignments", path)

    def test_it_still_refuses_when_nothing_is_usable(self):
        with io.open(os.path.join(self.dir, "x_observed_haplotypes_raw.zip"), "wb") as fh:
            fh.write(b"this is not a zip file")
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.probe_source(self.dir)
        self.assertIn("cannot be read", str(cm.exception))


class TestHiddenCountsTheSameHaplotypeTheCountersDo(HaplotypeScoresMixin,
                                                    unittest.TestCase):
    """``_hap`` tested ``variant_set`` for truthiness, but it is written as
    ``obs.variants or "-"`` and so is never empty -- the fallback to the aligned
    sequence was unreachable. Two variantless haplotypes with DIFFERENT sequences
    collapsed to one, so a cap hid a whole haplotype and reported 0 hidden."""

    def test_two_variantless_haplotypes_are_two(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, ""),
                _row("chr1", 100, HAP_TWO, 0.3, 0.8, "")]
        _res, text, msgs = self._run(rows, cap=1)
        body = self._parse(text)
        self.assertEqual(body[0]["n_haplotypes_at_site"], "2")
        note = [m for m in msgs if "not written" in m]
        self.assertTrue(note, msgs)
        self.assertNotIn("0 haplotype(s)", note[0],
                         "a whole haplotype was hidden and the message said 0")


class TestARefusalIsNeverErased(HaplotypeScoresMixin, unittest.TestCase):
    """The file's first invariant, re-broken twice while fixing it.

    Round 2 let a refusal on one axis mark the whole row unscored, so a haplotype 20x
    over the operating point counted as never measured. Round 3 let the surviving score
    erase the refusal, so a reached-and-declined axis rendered EXACTLY like one that has
    no column at all. The status now carries both.
    """

    def test_a_declined_cb_does_not_look_like_an_absent_column(self):
        """THE LIVE CASE: 139,656 rows of a real run. The CRISPR-Bulge axis was reached
        and declined; the other two blocks have no CRISPR-Bulge column at all. Both end
        with crispr_bulge=NA, so above_threshold is the only thing that can tell them
        apart -- and it read a bare 'unknown' on both, which the header defines as never
        measured."""
        declined = self._run([_row("chr1", 100, HAP_ONE, 0.6, -1.0, "chr1_100_A_C")])[1]
        body = self._parse(declined)
        cb_row = [r for r in body if "CRISPR_BULGE" in r["haplotype_id"]][0]
        other = [r for r in body if "CRISPR_BULGE" not in r["haplotype_id"]]
        self.assertIn("model-refused", cb_row["score_status"])
        self.assertEqual(cb_row["above_threshold"], "unknown(model-refused)")
        self.assertTrue(other)
        for r in other:
            self.assertEqual(r["above_threshold"], "unknown",
                             "never-measured must stay distinguishable from declined")

    def test_an_axis_with_no_column_still_says_plain_unknown(self):
        """The contrast that makes the value above mean something."""
        text = self._run([_row("chr1", 100, HAP_ONE, 0.6, 0.004, "chr1_100_A_C")])[1]
        plain = [r for r in self._parse(text)
                 if "CRISPR_BULGE" not in r["haplotype_id"]]
        self.assertTrue(plain)
        for r in plain:
            self.assertEqual(r["above_threshold"], "unknown")
            self.assertEqual(r["score_status"], "scored")

    def test_a_scored_with_refusal_row_counts_as_scored(self):
        text = self._run([_row("chr1", 100, HAP_ONE, -1.0, 0.004, "chr1_100_A_C")])[1]
        for r in self._parse(text):
            if r["crispr_bulge"] != "NA":
                self.assertEqual(r["n_not_scored"], "0")


class TestTheHaplotypeKeyIsCarriedNotRebuilt(HaplotypeScoresMixin, unittest.TestCase):
    """A source cell holding a literal "-" is truthy, so Site.add keyed on "-" while
    _hap mapped the written "-" back to empty and fell through to the sequence. A site
    then reported 1 haplotype, wrote 1, and said 1 was hidden."""

    def test_a_literal_dash_variant_cell_does_not_split_the_count(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "-"),
                _row("chr1", 100, HAP_TWO, 0.3, 0.8, "-")]
        _res, text, msgs = self._run(rows, cap=1)
        body = self._parse(text)
        total = int(body[0]["n_haplotypes_at_site"])
        note = [m for m in msgs if "not written" in m]
        hidden = int(note[0].split("covering ")[1].split(" haplotype")[0]) if note else 0
        self.assertEqual(total, len(body) + hidden,
                         "the site's own numbers do not add up: %d total, %d written, "
                         "%d hidden" % (total, len(body), hidden))

    def test_the_key_never_reaches_the_file(self):
        text = self._run([_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])[1]
        self.assertNotIn("_hap_key", text.splitlines()[
            next(i for i, l in enumerate(text.splitlines())
                 if l.startswith("site_key"))])


class TestADamagedPreferredSourceIsAnnounced(HaplotypeScoresMixin, unittest.TestCase):
    """An interrupted search leaves a truncated sidecar beside a complete alt table.
    Falling back silently leaves a damaged file sitting there unnoticed."""

    def test_the_fallback_says_the_preferred_source_is_damaged(self):
        with io.open(os.path.join(self.dir, "x_observed_haplotypes_raw.zip"), "wb") as fh:
            fh.write(b"this is not a zip file")
        path = os.path.join(self.dir, "x_all_results_with_alternative_alignments.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            fh.write("\t".join(_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")) + "\n")
        msgs = []
        _p, kind = hs.probe_source(self.dir, msgs.append)
        self.assertEqual(kind, "all_results_with_alternative_alignments")
        self.assertTrue([m for m in msgs if "cannot be read" in m], msgs)
        self.assertTrue([m for m in msgs if "should be" in m and "regenerated" in m], msgs)

    def test_a_clean_run_says_nothing(self):
        path = os.path.join(self.dir, "x_observed_haplotypes_raw.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            fh.write("\t".join(_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")) + "\n")
        msgs = []
        hs.probe_source(self.dir, msgs.append)
        self.assertEqual(msgs, [])


class TestTheThreeValuesAreEachExercised(HaplotypeScoresMixin, unittest.TestCase):
    """"THREE VALUES, NOT TWO" is the file's headline claim, and 'no' -- measured and
    BELOW -- had no positive assertion anywhere: flipping the branch to "yes" left the
    suite green. Same for needs-realignment and a non-zero n_not_scored, the module's
    stated contract for the Tier-2 case it deliberately did not build."""

    def test_measured_and_below_reads_no(self):
        thr = hs._threshold()[0]
        rows = [_row("chr1", 100, HAP_ONE, 0.5, thr / 190.0, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        cb = [r for r in self._parse(text) if "CRISPR_BULGE" in r["haplotype_id"]][0]
        self.assertEqual(cb["above_threshold"], "no")
        self.assertEqual(cb["score_status"], "scored")

    def test_at_the_threshold_exactly_reads_yes(self):
        thr = hs._threshold()[0]
        rows = [_row("chr1", 100, HAP_ONE, 0.5, thr, "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        cb = [r for r in self._parse(text) if "CRISPR_BULGE" in r["haplotype_id"]][0]
        self.assertEqual(cb["above_threshold"], "yes")

    def test_no_score_on_any_axis_is_needs_realignment_and_counted(self):
        rows = [_row("chr1", 100, HAP_ONE, "NA", "NA", "chr1_100_A_C")]
        _res, text, _msgs = self._run(rows)
        body = self._parse(text)
        self.assertTrue(body)
        for r in body:
            self.assertTrue(r["score_status"].startswith("needs-realignment:"),
                            r["score_status"])
            self.assertIn("chr1", r["score_status"], "the site key must be in it")
            self.assertEqual(r["n_not_scored"], "1")
            self.assertEqual(r["n_haplotypes_at_site"], "1")


class TestTheFourthReconciliationAbsence(HaplotypeScoresMixin, unittest.TestCase):
    """n/a asserts a MEASURED fact: "it reported a REFERENCE haplotype, so there is no
    ALT for it to appear among". Folding "the table has no origin column" and "the
    Variant cell did not parse" into it asserted a reason nobody read, left unreconciled
    empty, printed nothing, and made --strict exit 0 on a genuine failure."""

    def _targets_raw(self, header, row):
        path = os.path.join(self.dir, "t.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(header + "\n" + row + "\n")
        return path

    def test_an_absent_origin_column_is_not_n_a(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets_raw("Chromosome\tPosition\tVariant",
                                "chr1\t100\tchr1;100;A;T (rs9)")
        _res, text, msgs = self._run(rows, targets=tgt)
        vals = set(r["row_reported_alignment_found"] for r in self._parse(text))
        self.assertEqual(vals, {"unknown(origin-not-in-targets)"}, vals)
        self.assertTrue(any("no REF/ALT_origin column" in m for m in msgs), msgs)

    def test_an_unparsable_variant_on_an_alt_row_is_not_n_a(self):
        """generate_report._format_variant_id documents falling back to the rsID when no
        coordinate is available, so this is a shape the producer emits."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets_raw("Chromosome\tPosition\tVariant\tREF/ALT_origin",
                                "chr1\t100\trs10136\talt")
        _res, text, _msgs = self._run(rows, targets=tgt)
        vals = set(r["row_reported_alignment_found"] for r in self._parse(text))
        self.assertEqual(vals, {"unknown(unparsable-variant-in-targets)"}, vals)

    def test_a_real_ref_row_still_reads_n_a(self):
        """The contrast that makes n/a mean something."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "REF")])
        _res, text, _msgs = self._run(rows, targets=tgt)
        self.assertEqual(set(r["row_reported_alignment_found"]
                             for r in self._parse(text)), {"n/a"})


class TestAnEmptyTargetsTableIsNotNoTargets(HaplotypeScoresMixin, unittest.TestCase):
    """Branching on the truthiness of the parsed table made a --targets file that
    yielded zero usable rows indistinguishable from not passing --targets at all: the
    filter switched off and the WHOLE run was written, with none of the per-site
    accounting the flag exists for, and nothing said."""

    def test_a_targets_table_with_no_usable_rows_refuses(self):
        path = os.path.join(self.dir, "empty_targets.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("Chromosome\tPosition\tVariant\tREF/ALT_origin\n")
        self._write_source([_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C"),
                            _row("chr9", 900, HAP_TWO, 0.4, 0.9, "chr9_900_A_C")])
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.build(self.dir, os.path.join(self.dir, "o.tsv"),
                     log=lambda m: None, targets=path)
        self.assertIn("no usable rows", str(cm.exception))

    def test_it_does_not_quietly_write_the_whole_run(self):
        path = os.path.join(self.dir, "empty_targets.tsv")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("Chromosome\tPosition\tVariant\tREF/ALT_origin\n")
        out = os.path.join(self.dir, "o.tsv")
        self._write_source([_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C")])
        rc = hs.main(["--results", self.dir, "--out", out, "--targets", path])
        self.assertNotEqual(rc, 0)
        self.assertFalse(os.path.exists(out))


class TestAnAltSiteHoldingOnlyRefRowsIsAGap(HaplotypeScoresMixin, unittest.TestCase):
    """The curated row says ALT; this source holds only REF rows. The user's own table
    says a variant haplotype is there. Reporting "no variant haplotype here" is the worst
    answer available -- it is a clean, asserted against their own evidence."""

    def test_it_is_called_a_gap_not_a_clean(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", origin="ref")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertTrue([m for m in msgs if "ALT-origin in the curated table" in m], msgs)
        self.assertFalse([m for m in msgs if "no variant haplotype here" in m
                          and "NOT" not in m], msgs)

    def test_a_curated_ref_row_is_still_reported_as_present_with_no_alt(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", origin="ref")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "REF")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertTrue([m for m in msgs if "carry only reference rows" in m], msgs)
        self.assertFalse([m for m in msgs if "ALT-origin in the curated" in m], msgs)


class TestABudgetLabelIsNotAGate(HaplotypeScoresMixin, unittest.TestCase):
    """The budget is rendered into the header and quoted in score_status; nothing
    computes with it. Refusing the whole file over a missing label made the companion
    impossible for an assembly-search combined run, whose .Params.txt has no
    Max_total_edits by construction (pages/history_page.py carries a fallback for the
    same gap), and blamed search parameters that were all present."""

    def test_a_missing_max_total_edits_still_writes(self):
        with io.open(os.path.join(self.dir, ".Params.txt"), "w", encoding="utf-8") as fh:
            fh.write(PARAMS.replace("Max_total_edits\t8\n", ""))
        _res, text, _msgs = self._run([_row("chr1", 100, HAP_ONE, 0.4, 0.9,
                                            "chr1_100_A_C")])
        self.assertIn("?total", text, "the header must SHOW the gap, not hide it")
        self.assertTrue(self._parse(text))

    def test_a_missing_scorer_still_writes(self):
        with io.open(os.path.join(self.dir, ".Params.txt"), "w", encoding="utf-8") as fh:
            fh.write(PARAMS.replace("Scorer\tcrispr-bulge\n", ""))
        _res, text, _msgs = self._run([_row("chr1", 100, HAP_ONE, 0.4, 0.9,
                                            "chr1_100_A_C")])
        self.assertTrue(self._parse(text))

    def test_no_params_file_at_all_still_refuses(self):
        os.remove(os.path.join(self.dir, ".Params.txt"))
        with self.assertRaises(hs.HaplotypeScoreError):
            hs.build(self.dir, os.path.join(self.dir, "o.tsv"), log=lambda m: None)

    def test_the_budget_text_never_invents_a_number(self):
        self.assertEqual(hs.budget_text({"mismatches": 6, "dna": 1, "rna": 1}),
                         "6mm/1DNA/1RNA/?total")


class TestAContigNameWithAnUnderscoreKeepsItsAlleles(unittest.TestCase):
    """hg38 ships ~500 contigs whose NAME contains the field separator --
    chr1_KI270706v1_random, chrUn_KI270302v1. Splitting left-to-right took the first
    four pieces, so chr1_KI270706v1_random_12345_A_T reduced to
    chr1:KI270706v1:RANDOM:12345 with THE ALLELES DROPPED, and two different variants at
    that position compared EQUAL -- a site reconciling 'yes' against a haplotype it does
    not carry. A false yes is the one answer this file must never give."""

    def test_two_variants_at_one_position_on_an_alt_contig_do_not_collide(self):
        a = hs.canonical_variants("chr1_KI270706v1_random_12345_A_T")
        b = hs.canonical_variants("chr1_KI270706v1_random_12345_G_C")
        self.assertNotEqual(a, b, "different variants reduced to the same key")

    def test_the_alleles_survive(self):
        got = hs.canonical_variants("chr1_KI270706v1_random_12345_A_T")[0]
        self.assertTrue(got.endswith(":12345:A:T"), got)
        self.assertTrue(got.startswith("chr1_KI270706v1_random:"), got)

    def test_the_two_notations_still_agree_on_a_plain_contig(self):
        self.assertEqual(hs.canonical_variants("chr11_6482966_T_G"),
                         hs.canonical_variants("chr11;6482966;T;G (rs10136)"))

    def test_a_shape_it_cannot_read_is_kept_whole_not_invented(self):
        """Position must be numeric. If it is not, this is not the shape we think, and a
        made-up key could collide with a real one."""
        got = hs.canonical_variants("something_odd_here_now")
        self.assertEqual(got, ("SOMETHING_ODD_HERE_NOW",))


class TestTheTargetsAccountingKnowsTheGuide(HaplotypeScoresMixin, unittest.TestCase):
    """read_reported and site_key went per-guide when the multi-guide collision was
    fixed, but wanted/seen_anchors/ref_only_anchors stayed keyed on (chrom, start) -- so
    one guide's rows marked the coordinate 'seen' for every other guide curated at it,
    and a requested (locus, guide) that produced nothing vanished from a file whose
    header promises such sites are counted and named."""

    OTHER = "TTTTGGGGCCCCAAAATTTTNNN"

    def test_a_second_guide_at_the_same_locus_is_reported_as_a_gap(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT", SPACER),
                             ("chr1", "100", "chr1;100;A;T (rs9)", "ALT", self.OTHER)])
        _res, _text, msgs = self._run(rows, targets=tgt)
        gap = [m for m in msgs if "NO rows" in m]
        self.assertTrue(gap, "the second guide's site vanished silently: %s" % msgs)
        self.assertIn(self.OTHER, gap[0])

    def test_the_counts_in_the_two_log_lines_agree(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT", SPACER),
                             ("chr1", "100", "chr1;100;A;T (rs9)", "ALT", self.OTHER)])
        _res, _text, msgs = self._run(rows, targets=tgt)
        restricted = [m for m in msgs if "restricted to the" in m][0]
        self.assertIn("2 site(s)", restricted,
                      "two curated (locus,guide) pairs were counted as one")

    def test_a_guideless_table_is_unchanged(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertFalse([m for m in msgs if "NO rows" in m], msgs)


class TestTheAltVsRefVerdictDoesNotDependOnRowOrder(HaplotypeScoresMixin,
                                                    unittest.TestCase):
    """The ALT-in-targets / REF-only-in-source bucket used a lookup that fell back to
    'any curated row at this coordinate', so at a locus curated under two guides the
    verdict was whichever row the dict happened to yield: the SAME data gave 'GAP' or
    'no variant haplotype here' depending on the order of lines in the table."""

    OTHER = "TTTTGGGGCCCCAAAATTTTNNN"

    def _verdict(self, order):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C", origin="ref")]
        tgt = self._targets(order)
        _res, _text, msgs = self._run(rows, targets=tgt)
        return ("GAP" if [m for m in msgs if "ALT-origin in the curated table" in m]
                else ("CLEAN" if [m for m in msgs if "carry only reference rows" in m]
                      else "NONE"))

    def test_both_orders_give_the_same_answer(self):
        a = [("chr1", "100", "chr1;100;A;C (rs1)", "ALT", SPACER),
             ("chr1", "100", "chr1;100;A;T (rs9)", "REF", self.OTHER)]
        self.assertEqual(self._verdict(a), self._verdict(list(reversed(a))),
                         "the verdict flipped with the order of lines")

    def test_the_guide_that_was_curated_alt_is_the_one_reported(self):
        a = [("chr1", "100", "chr1;100;A;C (rs1)", "ALT", SPACER)]
        self.assertEqual(self._verdict(a), "GAP")

    def test_an_ambiguous_locus_is_in_neither_bucket(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C", origin="ref")]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT"),
                             ("chr1", "100", "chr1;100;A;T (rs9)", "REF")])
        _res, _text, msgs = self._run(rows, targets=tgt)
        self.assertTrue([m for m in msgs if "incompatibly" in m], msgs)
        self.assertFalse([m for m in msgs if "carry only reference rows" in m
                          and "incompatibly" not in m],
                         "an unknown was rendered as a measured clean")


class TestACutFileIsNotAnUnmeasuredRow(HaplotypeScoresMixin, unittest.TestCase):

    def test_a_row_cut_mid_line_refuses(self):
        path = os.path.join(self.dir, "x_observed_haplotypes_raw.tsv")
        full = "\t".join(_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"))
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(_header()) + "\n")
            fh.write(full + "\n")
            fh.write("\t".join(full.split("\t")[:9]) + "\n")       # cut mid-row
        with self.assertRaises(hs.HaplotypeScoreError) as cm:
            hs.build(self.dir, os.path.join(self.dir, "o.tsv"), log=lambda m: None)
        self.assertIn("truncated", str(cm.exception))


class TestANegativeCapMeansNoCap(HaplotypeScoresMixin, unittest.TestCase):
    """-1 is CRISPRme's own 'off' (--max-total-edits, --max-degeneracy), so a user asking
    for no cap got ranked[:-1]: every site truncated by one row, a single-row site
    emptied, silently, exit 0."""

    def test_minus_one_writes_every_row(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C")]
        _res, text, msgs = self._run(rows, cap=-1)
        self.assertTrue(self._parse(text), "a single-row site was emptied")
        self.assertTrue([m for m in msgs if "means NO cap" in m], msgs)

    def test_minus_one_matches_zero(self):
        rows = [_row("chr1", 100, HAP_ONE, 0.5, 0.004, "chr1_100_A_C"),
                _row("chr1", 100, HAP_TWO, 0.4, 0.003, "chr1_100_A_G")]
        a = len(self._parse(self._run(rows, cap=-1)[1]))
        self.setUp()
        b = len(self._parse(self._run(rows, cap=0)[1]))
        self.assertEqual(a, b)


class TestEveryCombinationUpholdsTheInvariants(HaplotypeScoresMixin, unittest.TestCase):
    """EXHAUSTIVE, not hand-picked.

    Five review rounds found the same class of defect over and over: two genuinely
    different situations rendering as the same cell, or a counter contradicting the rows
    beside it. Every one was found by a person noticing, which does not scale and did not
    converge -- rounds 4 and 5 each found defects inside the previous round's fixes.

    The input space is small and finite, so it is enumerated instead. Each case is built,
    run, and checked against the invariants the file states about itself. A defect of this
    class now fails here on the next run rather than waiting for a reviewer.
    """

    SCORES = ("0.5", "NA", "-1.0")          # measured / absent / model-refused sentinel
    ORIGINS = ("alt", "ref", "")
    CAPS = (0, 1, -1)

    def _one(self, cfd, cb, origin, cap, variants):
        self.setUp()
        rows = [_row("chr1", 100, HAP_ONE, cfd, cb, variants, origin=origin),
                _row("chr1", 100, HAP_TWO, cfd, cb, variants + ",chr1_101_G_T",
                     origin=origin)]
        try:
            _res, text, _msgs = self._run(rows, cap=cap)
        except hs.HaplotypeScoreError:
            return None                      # a refusal is always an allowed outcome
        return self._parse(text)

    def _cases(self):
        for cfd in self.SCORES:
            for cb in self.SCORES:
                for origin in self.ORIGINS:
                    for cap in self.CAPS:
                        for variants in ("chr1_100_A_C", "-", "NA"):
                            yield cfd, cb, origin, cap, variants

    def test_no_cell_is_ever_blank(self):
        n = 0
        for case in self._cases():
            body = self._one(*case)
            if not body:
                continue
            n += 1
            for r in body:
                blank = [k for k, v in r.items() if v == ""]
                self.assertFalse(blank, "blank %s on case %r" % (blank, (case,)))
        self.assertGreater(n, 50, "the enumeration produced almost no rows")

    def test_the_sentinel_is_never_printed(self):
        for case in self._cases():
            body = self._one(*case) or []
            for r in body:
                self.assertNotEqual(r["cfd"], "-1.0", case)
                self.assertNotEqual(r["crispr_bulge"], "-1.0", case)

    def test_every_written_value_is_declared_in_the_table(self):
        for case in self._cases():
            body = self._one(*case) or []
            for r in body:
                for col in ("score_status", "above_threshold",
                            "row_reported_alignment_found"):
                    hs._verdict(col, r[col])     # raises if undeclared

    def test_a_row_is_never_both_above_threshold_and_unscored(self):
        for case in self._cases():
            body = self._one(*case) or []
            for r in body:
                if r["above_threshold"] == "yes":
                    self.assertTrue(r["score_status"].startswith("scored"),
                                    "%s on case %r" % (r["score_status"], (case,)))

    def test_the_site_counters_add_up(self):
        """n_haplotypes_at_site is the TRUE total, so it must equal what was written
        plus what the cap hid -- the contradiction found twice by hand."""
        for case in self._cases():
            self.setUp()
            cfd, cb, origin, cap, variants = case
            rows = [_row("chr1", 100, HAP_ONE, cfd, cb, variants, origin=origin),
                    _row("chr1", 100, HAP_TWO, cfd, cb, variants + ",chr1_101_G_T",
                         origin=origin)]
            try:
                _res, text, msgs = self._run(rows, cap=cap)
            except hs.HaplotypeScoreError:
                continue
            body = self._parse(text)
            if not body:
                continue
            total = int(body[0]["n_haplotypes_at_site"])
            shown = len(set(r["variant_set"] for r in body))
            note = [m for m in msgs if "not written" in m]
            hidden = (int(note[0].split("covering ")[1].split(" haplotype")[0])
                      if note else 0)
            self.assertEqual(total, shown + hidden,
                             "case %r: total %d, shown %d, hidden %d"
                             % ((case,), total, shown, hidden))

    def test_above_threshold_agrees_with_the_number_beside_it(self):
        thr = hs._threshold()[0]
        for case in self._cases():
            body = self._one(*case) or []
            for r in body:
                if r["crispr_bulge"] == "NA":
                    self.assertTrue(r["above_threshold"].startswith("unknown"),
                                    "%s with no number, case %r"
                                    % (r["above_threshold"], (case,)))
                else:
                    want = "yes" if float(r["crispr_bulge"]) >= thr else "no"
                    self.assertEqual(r["above_threshold"], want, (case,))

    def test_the_composite_status_is_unreachable_today(self):
        """Mutation testing found no test could tell whether scored(<axis>-model-refused)
        existed. That is not a gap in the tests -- it is unreachable with the schema
        CRISPRme writes, because each block carries exactly one score column. Stated as a
        fact so it stays true: if a future schema puts both columns on one block, this
        fails and whoever does it has to decide deliberately."""
        cols = _header()
        for b in _BLOCKS:
            scores = [c for c in cols if "_score_" in c and c.endswith("(%s)" % b)]
            self.assertEqual(len(scores), 1,
                             "%s now has %s -- the composite status becomes REACHABLE "
                             "and needs a test of its own" % (b, scores))
        produced = set()
        for case in self._cases():
            for r in (self._one(*case) or []):
                produced.add(r["score_status"])
        self.assertFalse([v for v in produced if v.startswith("scored(")],
                         "the composite fired: %s" % produced)

    def test_declined_and_never_measured_never_render_alike(self):
        """The class behind four of the five rounds' findings.

        Asserted PER COLUMN, not on the row as a whole. A first version of this test
        compared the (score, above_threshold, status) tuple -- which still differed when
        above_threshold alone collapsed, so planting the exact round-3 regression left it
        green. A reader partitions on ONE column; the invariant has to hold there.
        """
        declined = self._one("0.5", "-1.0", "alt", 0, "chr1_100_A_C") or []
        absent = self._one("0.5", "NA", "alt", 0, "chr1_100_A_C") or []

        def cb_rows(body):
            return [r for r in body if "CRISPR_BULGE" in r["haplotype_id"]]

        d, a = cb_rows(declined), cb_rows(absent)
        self.assertTrue(d and a)
        # above_threshold is the column the header says tells these two apart.
        self.assertNotEqual(set(r["above_threshold"] for r in d),
                            set(r["above_threshold"] for r in a),
                            "above_threshold renders declined and never-measured alike")
        self.assertNotEqual(set(r["score_status"] for r in d),
                            set(r["score_status"] for r in a),
                            "score_status renders declined and never-measured alike")

    def test_each_verdict_column_distinguishes_what_its_table_says_it_does(self):
        """Generalises the above: for every pair of situations the VERDICTS table gives
        DIFFERENT meanings to, the column must actually differ between them."""
        cases = {
            "declined": self._one("0.5", "-1.0", "alt", 0, "chr1_100_A_C") or [],
            "never-measured": self._one("0.5", "NA", "alt", 0, "chr1_100_A_C") or [],
            "measured-below": self._one("0.5", "1e-09", "alt", 0, "chr1_100_A_C") or [],
            "measured-above": self._one("0.5", "0.5", "alt", 0, "chr1_100_A_C") or [],
        }
        seen = {}
        for name, body in cases.items():
            cb = [r for r in body if "CRISPR_BULGE" in r["haplotype_id"]]
            self.assertTrue(cb, name)
            val = cb[0]["above_threshold"]
            self.assertNotIn(val, seen,
                             "%s and %s both render above_threshold=%r"
                             % (name, seen.get(val), val))
            seen[val] = name


if __name__ == "__main__":
    unittest.main(verbosity=2)

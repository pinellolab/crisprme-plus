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
    """The REAL per-block schema.

    CFD_score_(...) exists on every block; CRISPR_BULGE_score_(...) exists ONLY on the
    CRISPR-Bulge block. The first version of this fixture invented a CRISPR-Bulge column
    on all three, so the suite exercised a shape the pipeline never writes -- and the
    rows that matter most, the ones carrying a score on one axis only, were never tested.
    """
    cols = ["Spacer+PAM", "Chromosome"]
    for b in _BLOCKS:
        cols += ["Start_coordinate_(%s)" % b, "Strand_(%s)" % b,
                 "Aligned_spacer+PAM_(%s)" % b,
                 "Aligned_protospacer+PAM_REF_(%s)" % b,
                 "Aligned_protospacer+PAM_ALT_(%s)" % b,
                 "REF/ALT_origin_(%s)" % b,
                 "CFD_score_(%s)" % b]
        if b == "highest_CRISPR_BULGE":
            cols += ["CRISPR_BULGE_score_(%s)" % b]
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
        out += [str(st), "+", SPACER, REPORTED, alt, origin, str(cfd)]
        if b == "highest_CRISPR_BULGE":
            out += [str(cb)]
        out += [variants, "0.01", samples]
    return out


class HaplotypeScoresFixture(unittest.TestCase):

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

    def _run(self, rows, **kw):
        self._write_source(rows)
        out = os.path.join(self.dir, "out.tsv")
        msgs = []
        res = hs.build(self.dir, out, log=msgs.append, **kw)
        with io.open(out, "r", encoding="utf-8") as fh:
            text = fh.read()
        return res, text, msgs

    # -- the reconciliation signal ------------------------------------------ #

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
        self.assertTrue(any("mis-keyed comparison" in m for m in msgs),
                        "it must say the key is wrong")
        body = [l for l in text.splitlines() if l and not l.startswith("#")]
        head, data = body[0].split("\t"), [l.split("\t") for l in body[1:]]
        i = head.index("row_reported_alignment_found")
        self.assertTrue(all(r[i] == "unknown" for r in data))

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


class TestSiteIdentityAcrossDivergentBlocks(HaplotypeScoresFixture):
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
        keys = set(l.split("\t")[0] for l in text.splitlines()
                   if l and not l.startswith("#") and not l.startswith("site_key"))
        self.assertEqual(len(keys), 1, keys)

    def test_each_row_still_reports_its_own_blocks_coordinate(self):
        """Anchoring site IDENTITY must not rewrite the coordinate a block actually used."""
        rows = [_row("chr1", 100, HAP_ONE, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        _res, text, _msgs = self._run(rows)
        body = [l.split("\t") for l in text.splitlines()
                if l and not l.startswith("#") and not l.startswith("site_key")]
        starts = set(r[2] for r in body)           # site_start
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
        body = [l.split("\t") for l in text.splitlines()
                if l and not l.startswith("#") and not l.startswith("site_key")]
        self.assertTrue(body, "the CRISPR-Bulge observation was dropped entirely")
        cb = set(r[11] for r in body)              # crispr_bulge
        self.assertIn("0.9", cb,
                      "the CRISPR-Bulge score was lost with the block it lives on: %s" % cb)
        self.assertEqual(set(r[28] for r in body), {"0.9"},
                         "site_max_crispr_bulge must see the score too")

    def test_reconciliation_survives_the_divergence(self):
        """The reported haplotype is found even though one block sits at another start."""
        rows = [_row("chr1", 100, REPORTED, 0.4, 0.9, "chr1_100_A_C", cb_start=97)]
        tgt = self._targets([("chr1", "100", "chr1;100;A;C (rs1)", "ALT")])
        (_n, _s, unreconciled), text, _msgs = self._run(rows, targets=tgt)
        self.assertEqual(unreconciled, [], text)
        self.assertIn("\tyes", text)


class TestTargetsKeyIsPerGuide(HaplotypeScoresFixture):
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

    def test_cfd_is_on_every_block(self):
        cols = _header()
        self.assertEqual(
            sorted(c for c in cols if c.startswith("CFD_score_")),
            sorted("CFD_score_(%s)" % b for b in _BLOCKS))


class TestAnEmptySourceNeverWins(HaplotypeScoresFixture):
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


class TestTargetsGapsAreNamed(HaplotypeScoresFixture):
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


if __name__ == "__main__":
    unittest.main(verbosity=2)

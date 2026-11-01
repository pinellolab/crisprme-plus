#!/usr/bin/env python3
"""The "Top 1000 candidate off-targets" plots (CRISPRme_plots.py, shown on the
results page's Graphical Reports tab) and the Personal Risk Card "Top 100" plots
(CRISPRme_plots_personal.py) draw EVERY site.

They dropped every site with mismatches+bulges <= 1 ("probably on-target introduced
by variants"), so all 1-mismatch-or-bulge off-targets and every perfect match, the
on-target among them, were missing, and nothing on the plot said so: a variant-created
perfect match, as strong a target as the on-target, never appeared. Now the 1-edit
sites are plotted like the others and the perfect matches are starred, with a legend
entry. The scripts still run from the command line as before, also for one site.

The shareable report (generate_report.py) draws its own top plots (Section 2) and
top-1000 tables (Section 6, top1000.tsv) and still leaves those sites out there; it
now says so under them, with how many and where they are (TestReportTopViews).

Run with:

    cd PostProcess && python -m unittest test_top_plots -v
"""

import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import CRISPRme_plots as plots  # noqa: E402
import CRISPRme_plots_personal as personal  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCACNNN"
CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _sites():
    """One table, as the pipeline writes it, of sites with 0 to 4 mismatches+bulges:
    the reference on-target, a variant-created perfect match (REF 1 mismatch), two
    1-edit sites (one variant), and three farther ones."""
    rows = [
        # mm+b, origin, CFD REF, CFD ALT, MAF, carriers, mismatches in the REF protospacer
        (0, "ref", 1.0, 1.0, None, None, "CTAACAGTTGCTTTTATCACTGG"),
        (0, "alt", 0.93, 1.0, 0.58, "S1,S2", "CTAACAGTTaCTTTTATCACTGG"),
        (1, "ref", 0.64, 0.64, None, None, "CTAACAGTTGCTTTTAaCACTGG"),
        (1, "alt", 0.01, 0.54, 0.25, "S4", "CTAACAGTTGCaTTTATtACTGG"),
        (2, "ref", 0.16, 0.16, None, None, "CTAACAGTTGCTTTTAatACTGG"),
        (3, "ref", 0.15, 0.15, None, None, "CTAACAGTTGCTTTaAatACTGG"),
        (4, "ref", 0.01, 0.01, None, None, "CTAACAGTTGCTTaaAatACTGG"),
    ]
    data = []
    for i, (mmb, origin, ref, alt, maf, samples, prot) in enumerate(rows):
        row = {"Spacer+PAM": GUIDE, "Chromosome": "chr1", "Start_coordinate": 1000 * (i + 1)}
        for crit in CRITERIA:
            row[f"Mismatches+bulges_({crit})"] = mmb
            row[f"Bulges_({crit})"] = 0
            row[f"REF/ALT_origin_({crit})"] = origin
            row[f"Variant_MAF_({crit})"] = maf
            row[f"Variant_samples_({crit})"] = samples
            row[f"Aligned_protospacer+PAM_REF_({crit})"] = prot
        row["CFD_score_(highest_CFD)"] = max(ref, alt)
        row["CFD_score_REF_(highest_CFD)"] = ref
        row["CFD_score_ALT_(highest_CFD)"] = alt
        row["CRISPR_BULGE_score_(highest_CRISPR_BULGE)"] = max(ref, alt)
        row["CRISPR_BULGE_score_REF_(highest_CRISPR_BULGE)"] = ref
        row["CRISPR_BULGE_score_ALT_(highest_CRISPR_BULGE)"] = alt
        data.append(row)
    return pd.DataFrame(data)


def _drawn(plot, df, *args, **kwargs):
    """Run a plot function and return, per scatter layer, the points it draws, plus
    the legend labels -- captured at savefig, before the figure is cleared."""
    seen = {}

    def _capture(*_a, **_k):
        if "layers" in seen:
            return
        ax = plt.gca()
        # REF points (zorder 1), ALT points (2), and the perfect-match stars (3)
        seen["layers"] = [(c.get_offsets().tolist(), c.get_zorder()) for c in ax.collections]
        seen["star_layers"] = [c.get_offsets().tolist() for c in ax.collections
                               if c.get_zorder() == 3]
        legend = ax.get_legend()
        seen["legend"] = [t.get_text() for t in legend.get_texts()] if legend else []

    with tempfile.TemporaryDirectory() as out, patch.object(plt, "savefig", _capture):
        plot(df, out + os.sep, *args, **kwargs)
    return seen


class TestTopPlots(unittest.TestCase):
    def test_filter_keeps_every_site(self):
        df = _sites()
        for term in ("cfd", "crispr_bulge", "mm+b"):
            with self.subTest(term=term):
                kept = plots.filter_table(df, term)
                self.assertEqual(len(kept), len(df))
                mmb = kept[{"cfd": "Mismatches+bulges_(highest_CFD)",
                            "crispr_bulge": "Mismatches+bulges_(highest_CRISPR_BULGE)",
                            "mm+b": "Mismatches+bulges_(fewest_mm+b)"}[term]]
                self.assertEqual(sorted(mmb.tolist()), [0, 0, 1, 1, 2, 3, 4])
        delta = plots.filter_table(df, "cfd", sort_by="delta")
        self.assertEqual(len(delta), len(df))
        # the variant-created perfect match and the variant 1-edit site lead the
        # variant-effect ranking (largest ALT-REF first)
        self.assertEqual(delta["CFD_score_ALT_(highest_CFD)"].tolist()[:2], [0.54, 1.0])

    def test_score_plots_draw_every_site_and_star_the_perfect_matches(self):
        df = _sites()
        for plot, term, score in (
            (plots.plot_with_CFD_score, "cfd", "CFD_score_(highest_CFD)"),
            (plots.plot_with_CRISPR_BULGE_score, "crispr_bulge",
             "CRISPR_BULGE_score_(highest_CRISPR_BULGE)"),
        ):
            with self.subTest(plot=plot.__name__):
                ranked = plots.filter_table(df, term)
                seen = _drawn(plot, ranked.copy(), GUIDE)
                red = seen["layers"][0][0]  # every site has a REF point
                self.assertEqual(len(red), len(df))
                # one star per perfect match, at its best score, at its rank
                stars = seen["star_layers"]
                self.assertEqual(len(stars), 1)
                ranks = [i + 1 for i, v in enumerate(ranked[
                    f"Mismatches+bulges_({'highest_CFD' if term == 'cfd' else 'highest_CRISPR_BULGE'})"
                ].tolist()) if v == 0]
                self.assertEqual(sorted(p[0] for p in stars[0]), ranks)
                self.assertEqual([p[1] for p in stars[0]], [1.0, 1.0])
                self.assertIn(plots.PERFECT_MATCH_LABEL, seen["legend"])

    def test_edit_count_plot_draws_every_site(self):
        ranked = plots.filter_table(_sites(), "mm+b")
        seen = _drawn(plots.plot_with_MMvBUL, ranked.copy(), GUIDE)
        self.assertEqual(len(seen["layers"][0][0]), 7)
        self.assertEqual(sorted(p[1] for p in seen["star_layers"][0]), [0.0, 0.0])
        self.assertIn(plots.PERFECT_MATCH_LABEL, seen["legend"])

    def test_no_star_and_no_legend_entry_without_a_perfect_match(self):
        df = _sites()
        df = df[df["Mismatches+bulges_(highest_CFD)"] > 0]
        seen = _drawn(plots.plot_with_CFD_score, plots.filter_table(df, "cfd").copy(), GUIDE)
        self.assertEqual(seen["star_layers"], [])
        self.assertNotIn(plots.PERFECT_MATCH_LABEL, seen["legend"])
        self.assertEqual(len(seen["layers"][0][0]), 5)  # the 1-edit sites stay

    def test_personal_card_plots_draw_every_site(self):
        df = _sites()
        for plot in (personal.plot_with_CFD_score, personal.plot_with_CRISPR_BULGE_score,
                     personal.plot_with_MMvBUL):
            with self.subTest(plot=plot.__name__):
                seen = _drawn(plot, df, GUIDE)
                self.assertEqual(len(seen["layers"][0][0]), 7)
                self.assertEqual(len(seen["star_layers"][0]), 2)
                self.assertIn(personal.PERFECT_MATCH_LABEL, seen["legend"])
        # the caller's table is left as it was (the three plots share it)
        self.assertEqual(df["Mismatches+bulges_(highest_CFD)"].tolist(), [0, 0, 1, 1, 2, 3, 4])
        self.assertNotIn("index", df.columns)

    def test_personal_card_plot_of_no_site(self):
        # a sample with no private target: an empty plot, not a failure
        empty = _sites().iloc[0:0]
        for plot in (personal.plot_with_CFD_score, personal.plot_with_CRISPR_BULGE_score,
                     personal.plot_with_MMvBUL):
            with self.subTest(plot=plot.__name__):
                seen = _drawn(plot, empty, GUIDE)
                self.assertEqual(seen["star_layers"], [])

    def test_a_single_site_is_plotted(self):
        # A search whose only site is the reference on-target (Reference only, 0
        # mismatches, 0 bulges), or one reference off-target: a site without a variant
        # has no ALT point (its size was NaN), and a plot of that one site alone failed
        # in savefig ("need at least one array to concatenate"). The script stopped
        # there without an image, and the shell dropped its error with warnings.txt.
        sites = _sites()
        for row in (0, 2, 4):  # 0, 1 and 2 mismatches+bulges, all reference
            one = sites.iloc[[row]]
            mmb = int(one["Mismatches+bulges_(highest_CFD)"].iloc[0])
            with self.subTest(mmb=mmb), tempfile.TemporaryDirectory() as tmp:
                tsv = os.path.join(tmp, "integrated.tsv")
                one.to_csv(tsv, sep="\t", index=False, na_rep="NA")
                out = os.path.join(tmp, "imgs") + os.sep
                os.makedirs(out)
                run = subprocess.run(
                    [sys.executable, os.path.join(HERE, "CRISPRme_plots.py"), tsv, out],
                    capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr[-600:])
                made = sorted(f for f in os.listdir(out) if f.endswith(".png"))
                self.assertEqual(len(made), 5, made)
                run = subprocess.run(
                    [sys.executable, os.path.join(HERE, "CRISPRme_plots_personal.py"), tsv,
                     out, f"{GUIDE}.S1.personal"], capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr[-600:])
                self.assertEqual(
                    len([f for f in os.listdir(out) if f.endswith(".personal.png")]), 3)
                # drawn: its REF point, and a star if it is a perfect match
                seen = _drawn(plots.plot_with_CFD_score, plots.filter_table(one, "cfd").copy(),
                              GUIDE)
                self.assertEqual(len(seen["layers"][0][0]), 1)
                self.assertEqual(len(seen["star_layers"]), 1 if mmb == 0 else 0)

    def test_scripts_still_run_from_the_command_line(self):
        df = _sites()
        with tempfile.TemporaryDirectory() as tmp:
            tsv = os.path.join(tmp, "integrated.tsv")
            df.to_csv(tsv, sep="\t", index=False, na_rep="NA")
            out = os.path.join(tmp, "imgs") + os.sep
            os.makedirs(out)
            subprocess.run([sys.executable, os.path.join(HERE, "CRISPRme_plots.py"), tsv, out],
                           check=True, capture_output=True)
            made = sorted(os.listdir(out))
            for kind in ("CFD", "CRISPR_BULGE", "fewest"):
                self.assertIn(f"CRISPRme_{kind}_top_1000_log_for_main_text_{GUIDE}.png", made)
            for kind in ("CFD", "CRISPR_BULGE"):
                self.assertIn(f"CRISPRme_{kind}_top_1000_by_variant_effect_{GUIDE}.png", made)
            subprocess.run([sys.executable, os.path.join(HERE, "CRISPRme_plots_personal.py"),
                            tsv, out, f"{GUIDE}.S1.personal"], check=True, capture_output=True)
            for kind in ("CFD", "CRISPR_BULGE", "fewest"):
                self.assertIn(f"CRISPRme_{kind}_top_1000_log_for_main_text_{GUIDE}.S1.personal.png",
                              os.listdir(out))


# A results table as the report reads it (columns by name): the on-target and a
# variant-created perfect match, a reference and a variant 1-edit off-target, and
# three farther sites.
_REPORT_COLUMNS = [
    "Spacer+PAM", "Chromosome", "Start_coordinate_(highest_CFD)", "Strand_(highest_CFD)",
    "Aligned_protospacer+PAM_REF_(highest_CFD)", "Aligned_protospacer+PAM_ALT_(highest_CFD)",
    "PAM_(highest_CFD)", "Mismatches_(highest_CFD)", "Bulges_(highest_CFD)",
    "Mismatches+bulges_(highest_CFD)", "REF/ALT_origin_(highest_CFD)",
    "PAM_creation_(highest_CFD)", "CFD_score_(highest_CFD)", "CFD_score_REF_(highest_CFD)",
    "CFD_score_ALT_(highest_CFD)", "Variant_info_genome_(highest_CFD)",
    "Variant_MAF_(highest_CFD)", "Variant_rsID_(highest_CFD)",
    "Variant_samples_(highest_CFD)", "Not_found_in_REF",
]
_REPORT_ROWS = [
    # chrom, pos, REF protospacer, ALT protospacer, mm, bulges, origin, CFD, CFD REF,
    # variant, MAF, rsID, carriers
    ("chr1", "5000", "CTAACAGTTGCTTTTATCACTGG", "NA", 0, 0, "ref", 1.0, 1.0,
     "NA", "NA", "NA", "NA"),
    ("chr3", "10000", "CTAACAGTTaCTTTTATCACTGG", "CTAACAGTTGCTTTTATCACTGG", 0, 0, "alt",
     1.0, 0.93, "chr3_10009_A_G", "0.25", "rs1", "S1,S2"),
    ("chr2", "30000", "CTAAtAGTTGCTTTTATCACAGG", "NA", 1, 0, "ref", 0.636, 0.636,
     "NA", "NA", "NA", "NA"),
    ("chr1", "40000", "CTAACAGTTGCTTTTATCtaTGG", "CTAACAGTTGCTTTTATCtCTGG", 1, 0, "alt",
     0.538, 0.012, "chr1_40019_A_C", "0.1", "rs2", "S4"),
    ("chr1", "30000", "CTAACAGTTGcTTTATCACTGG", "NA", 1, 1, "ref", 0.16, 0.16,
     "NA", "NA", "NA", "NA"),
    ("chr2", "5000", "tTAACAGgTGCTTcTATCAgTGG", "NA", 4, 0, "ref", 0.012, 0.012,
     "NA", "NA", "NA", "NA"),
    ("chr3", "30000", "tcACAGTTtCTTTTAaCACTAG", "tcACAGTTtCTTTTAgCACTAG", 4, 1, "alt",
     0.086, 0.01, "chr3_30016_A_G", "0.05", "rs3", "S2"),
]


@unittest.skipUnless(os.path.isfile(os.path.join(HERE, "generate_report.py")),
                     "generate_report.py not found")
class TestReportTopViews(unittest.TestCase):
    """The shareable report leaves the sites with 0 or 1 mismatches + bulges out of
    its Section-2 plots and Section-6 tables (and top1000.tsv), as before; it now
    says so under each, with how many and where they are: all in the validation
    panel (panel_top100.tsv), and the 1-edit off-targets also in mmb_le_1.tsv."""

    def setUp(self):
        import generate_report as gr
        self.gr = gr
        gr._DROP_MAF = False
        gr._PRESENT_ANN_KINDS = None
        self.tmp = tempfile.mkdtemp(prefix="top_views_")
        self.addCleanup(__import__("shutil").rmtree, self.tmp, ignore_errors=True)
        self.sid_dir = os.path.join(self.tmp, "samplesIDs")
        os.makedirs(self.sid_dir)
        with open(os.path.join(self.sid_dir, "hg38_walk.samplesID.txt"), "w") as h:
            h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
            for s in ("S1", "S2", "S4"):
                h.write(f"{s}\tPOPA\tSUPX\tmale\n")

    def build(self, rows):
        import zipfile
        tsv = os.path.join(self.tmp, f"{GUIDE}+NRG_hg38+hg38_walk_4+2_integrated_results.tsv")
        with open(tsv, "w") as h:
            h.write("\t".join(_REPORT_COLUMNS) + "\n")
            for (chrom, pos, ref, alt, mm, bul, origin, cfd, cfd_ref, var, maf, rsid,
                 samples) in rows:
                h.write("\t".join(str(v) for v in (
                    GUIDE, chrom, pos, "+", ref, alt, ref[-3:], mm, bul, mm + bul, origin,
                    "NA", cfd, cfd_ref, cfd, var, maf, rsid, samples,
                    "y" if origin == "alt" else "NA")) + "\n")
        out = self.gr.build_report(integrated_tsv=tsv, samplesid_dir=self.sid_dir,
                                   out_zip=os.path.join(self.tmp, "report.zip"))
        with zipfile.ZipFile(out) as z:
            files = {n: z.read(n).decode() for n in z.namelist()}
        return files

    @staticmethod
    def section(html, heading):
        """The text of a section from its heading to its first subsection."""
        return html.split(f"<h2>{heading}</h2>", 1)[1].split("<h3", 1)[0]

    @staticmethod
    def sites(tsv_text):
        rows = [r.split("\t") for r in tsv_text.strip().splitlines()]
        head = rows[0]
        return [(r[head.index("Chromosome")], r[head.index("Position")],
                 int(r[head.index("Mismatches+bulges")])) for r in rows[1:]]

    def test_the_plots_and_the_tables_say_what_they_leave_out(self):
        files = self.build(_REPORT_ROWS)
        html = files["report.html"]
        found = ("the sites with 0 or 1 mismatches + bulges, here 2 perfect matches "
                 "(0 mismatches and bulges) and 2 off-targets with 1 mismatch or bulge. "
                 "All 4 are in the recommended validation panel (Section 4, ")
        panel = '<a href="data/panel_top100.tsv" download><code>panel_top100.tsv</code></a>'
        one = '<a href="data/mmb_le_1.tsv" download><code>mmb_le_1.tsv</code></a>'
        for heading, view in (
            ("2. Key graphical report", "Not in these plots: "),
            ("6. Top 1000 putative off-targets",
             'Not in this table or its download (<a href="data/top1000.tsv" download>'
             "<code>top1000.tsv</code></a>): "),
        ):
            with self.subTest(section=heading):
                text = self.section(html, heading)
                self.assertIn(view + found + panel + "), the perfect matches first; the 2 "
                              "with 1 mismatch or bulge are also in " + one + ".", text)
        # the views themselves are as before: no site with 0 or 1 in top1000.tsv ...
        top = self.sites(files["data/top1000.tsv"])
        self.assertEqual(sorted(m for *_, m in top), [2, 4, 5])
        # ... and the caption is right: every one is in the panel, the 1-edit ones in
        # mmb_le_1.tsv
        left_out = {("chr1", "5000"), ("chr3", "10000"), ("chr2", "30000"), ("chr1", "40000")}
        in_panel = {(c, p) for c, p, _ in self.sites(files["data/panel_top100.tsv"])}
        self.assertLessEqual(left_out, in_panel)
        self.assertEqual({(c, p) for c, p, _ in self.sites(files["data/mmb_le_1.tsv"])},
                         {("chr2", "30000"), ("chr1", "40000")})

    def test_the_personal_assembly_report_says_it_too(self):
        # build_combined_report's Section-5 tables come from select_top too, and its
        # validation panel is Section 3 (checked end to end with assembly-search)
        import inspect
        src = inspect.getsource(self.gr.build_combined_report)
        self.assertIn("top_exclusion_note(", src)
        self.assertIn("panel_section=3", src)
        self.assertIn("{top_exclusion_html}", src)

    def test_nothing_left_out_nothing_said(self):
        files = self.build([r for r in _REPORT_ROWS if r[4] + r[5] > 1])
        self.assertNotIn("Not in these plots", files["report.html"])
        self.assertNotIn("Not in this table", files["report.html"])

    def test_note_wording(self):
        df = pd.DataFrame({"mmb": [0, 2, 3]})
        cols = {"mmb": "mmb"}
        note = self.gr.top_exclusion_note(df, cols, "these plots", "panel_top100.tsv")
        self.assertEqual(
            note, '<p class="caption">Not in these plots: the sites with 0 or 1 mismatches '
            "+ bulges, here 1 perfect match (0 mismatches and bulges). It is in the "
            "recommended validation panel (Section 4, <code>panel_top100.tsv</code>).</p>")
        df = pd.DataFrame({"mmb": [1, 2]})
        note = self.gr.top_exclusion_note(df, cols, "this table", None, 3, "mmb_le_1.tsv")
        self.assertIn("here 1 off-target with 1 mismatch or bulge. The one with 1 "
                      "mismatch or bulge is in <code>mmb_le_1.tsv</code>.", note)
        note = self.gr.top_exclusion_note(df, cols, "this table", "panel_top100.tsv", 3)
        self.assertIn("It is in the recommended validation panel (Section 3, ", note)
        self.assertEqual(self.gr.top_exclusion_note(df, {}, "this table"), "")
        df = pd.DataFrame({"mmb": [1, 0, 1]})
        note = self.gr.top_exclusion_note(df, cols, "these plots", "panel_top100.tsv")
        self.assertIn("here 1 perfect match (0 mismatches and bulges) and 2 off-targets "
                      "with 1 mismatch or bulge. All 3 are in the recommended validation "
                      "panel (Section 4, <code>panel_top100.tsv</code>), the perfect match "
                      "first.</p>", note)


if __name__ == "__main__":
    unittest.main()

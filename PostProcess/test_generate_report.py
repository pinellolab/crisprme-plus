#!/usr/bin/env python3
"""Unit test for the shareable-report generator (PostProcess/generate_report.py).

Builds a tiny integrated_results fixture (dict-based schema with CRISTA + PAM
creation columns, header names only for the columns the report reads) and
asserts that build_report produces a ZIP that flat-decompresses to exactly the
bundled files -- report.html + integrated_results.tsv.gz + top1000.tsv +
panel_top100.tsv + the per-tier "ready to go" TSVs (v2.2) -- and that the HTML
is a self-contained IND-briefing-book digest (report v2.2):

  * SECTION 1: the header summary card (guide, PAM, aggregated specificity score)
    AND the global Off-targets-by-MM-and-B matrix are present; the matrix
    REFERENCE + VARIANT totals reconcile to the canonical partition off-target
    total (REFERENCE + VARIANT + on-target == grand total),
  * the canonical partition invariant holds:
    variant-created + reference + on-target == total EXACTLY,
  * SECTION 2: exactly FOUR ref/alt scatter panels are embedded as inline base64
    PNGs when CRISTA is computed (by CFD, by CFD delta, by CRISTA, by CRISTA
    delta); 2 panels when CRISTA is absent,
  * SECTION 4: the suggested panel is the worst-case top-100 (a site worst by any
    single metric -- CFD, CRISTA, or mm+b -- is included; capped at 100),
  * SECTION 4/5 (v2.2): EACH validation tier (CFD>= {0.5,0.2,0.1,0.05},
    mm+b<= {1,2,3,4}, variant-created, + the worst-case top-100) is exported as
    its own ready-to-use TSV bundled in the zip, linked RELATIVELY in the
    section-4 Download column AND section-5; 0-row tiers (mm+b<=1 here) are
    skipped (no file, no link); each tier is column-identical to top1000.tsv,
  * SECTION 5/6: the top-N off-target table is rendered inline (rows sorted by
    CFD) with a PAM-creation column, and top1000.tsv + panel_top100.tsv are
    bundled in the zip; the concise MAF footnote (blank/em-dash explanation) is
    present under the summary and the top-1000 table,
  * SECTION 7: the fixed research-only disclaimer is in the footer,
  * there is no external dependency (<script>/<link>/http(s)), so it opens with
    file:// offline,
  * the download links are RELATIVE (resolve after unzip on any machine).

Dependency note: generate_report imports matplotlib + pandas at module top, and
the network-free CI env installs only requests + numpy. The whole test therefore
self-skips (never fails) when matplotlib/pandas are unavailable; it runs fully in
the real pipeline env (and locally), which is where the artifact is produced.
"""

import base64
import glob
import gzip
import os
import re
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:  # the report's real deps (present in the pipeline env, absent in light CI)
    import matplotlib  # noqa: F401
    import pandas  # noqa: F401

    import generate_report as gr

    _HAVE_DEPS = True
    _SKIP_REASON = ""
except Exception as exc:  # noqa: BLE001 - any import failure => skip, not fail
    _HAVE_DEPS = False
    _SKIP_REASON = f"matplotlib/pandas not available: {exc}"


# A minimal fixture that exercises the columns the report reads by NAME. Two
# reference off-targets and three variant-created ones (one multi-SNP haplotype
# with comma-joined rsID/MAF/samples), plus one on-target (mm+b == 0) that must
# be excluded from the top-N by the mm+b > 1 filter. The schema includes the
# highest_CFD projection (with CFD REF/ALT + PAM creation) and the highest_CRISTA
# projection (non-empty), so the report emits the CRISTA scatter + CRISTA column.
_HEADER = [
    "Spacer+PAM",
    "Chromosome",
    "Start_coordinate_(highest_CFD)",
    "Strand_(highest_CFD)",
    "Aligned_protospacer+PAM_REF_(highest_CFD)",
    "Aligned_protospacer+PAM_ALT_(highest_CFD)",
    "PAM_(highest_CFD)",
    "Mismatches_(highest_CFD)",
    "Bulges_(highest_CFD)",
    "Mismatches+bulges_(highest_CFD)",
    "REF/ALT_origin_(highest_CFD)",
    "PAM_creation_(highest_CFD)",
    "CFD_score_(highest_CFD)",
    "CFD_score_REF_(highest_CFD)",
    "CFD_score_ALT_(highest_CFD)",
    "Variant_info_genome_(highest_CFD)",
    "Variant_MAF_(highest_CFD)",
    "Variant_rsID_(highest_CFD)",
    "Variant_samples_(highest_CFD)",
    "CRISTA_score_(highest_CRISTA)",
    "CRISTA_score_REF_(highest_CRISTA)",
    "CRISTA_score_ALT_(highest_CRISTA)",
    "Variant_MAF_(highest_CRISTA)",
    "Variant_rsID_(highest_CRISTA)",
    "Variant_samples_(highest_CRISTA)",
    "Not_found_in_REF",
    "Annotation_closest_gene_name",
    "Annotation_closest_gene_distance_(kb)",
]

_GUIDE = "CTCTCAGCTGGTACACGGCANNN"

# columns: guide, chrom, pos, strand, aln_ref, aln_alt, pam, mm, bul, mmb,
#          origin, pam_creation, cfd, cfd_ref, cfd_alt, var_genome, maf, rsid,
#          samples, crista, crista_ref, crista_alt, crista_maf, crista_rsid,
#          crista_samples, not_in_ref, gene, dist
_ROWS = [
    # on-target (mm+b == 0) -> excluded from top-N; counted as 1 on-target
    [_GUIDE, "chr2", "100", "+", "CTCTCAGCTGGTACACGGCATGG", "NA", "TGG",
     "0", "0", "0", "ref", "NA", "1.0", "1.0", "1.0", "NA", "NA", "NA", "NA",
     "1.0", "1.0", "1.0", "NA", "NA", "NA", "NA", "GENE0", "0.0"],
    # reference off-targets
    [_GUIDE, "chr3", "200", "-", "cTCTCAGCTGGTACACGGCAAGG", "NA", "AGG",
     "2", "0", "2", "ref", "NA", "0.85", "0.85", "0.85", "NA", "NA", "NA", "NA",
     "0.80", "0.80", "0.80", "NA", "NA", "NA", "NA", "GENE1", "1.2"],
    [_GUIDE, "chr4", "300", "+", "ctCTCAGCTGGTACACGGCAcGG", "NA", "CGG",
     "3", "0", "3", "ref", "NA", "0.40", "0.40", "0.40", "NA", "NA", "NA", "NA",
     "0.35", "0.35", "0.35", "NA", "NA", "NA", "NA", "GENE2", "5.0"],
    # variant-created off-targets (1000G + HGDP carriers)
    [_GUIDE, "chr5", "400", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAAGG", "AGG",
     "2", "0", "2", "alt", "NA", "0.90", "0.20", "0.90", "chr5_400_T_A",
     "0.01", "rs111", "HG00096,HGDP00001",
     "0.88", "0.18", "0.88", "0.01", "rs111", "HG00096,HGDP00001",
     "y", "GENE3", "0.5"],
    [_GUIDE, "chr6", "500", "-", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAcGG", "CGG",
     "3", "1", "4", "alt", "pam_created", "0.55", "0.30", "0.55", "chr6_500_T_C",
     "0.02", "rs222", "NA18525",
     "0.50", "0.25", "0.50", "0.02", "rs222", "NA18525",
     "y", "GENE4", "2.0"],
    # multi-SNP haplotype: comma-joined rsID/MAF/samples (min-AF + first-rsID),
    # and a BLANK MAF -> em-dash + footnote path
    [_GUIDE, "chr7", "600", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAgGG", "GGG",
     "3", "0", "3", "alt", "NA", "0.30", "0.10", "0.30", "chr7_600_T_G,chr7_601_A_C",
     "NA,0.005,0.003", "NA,rs333,rs444", "HG00097,HGDP00003",
     "0.28", "0.08", "0.28", "NA,0.005,0.003", "NA,rs333,rs444",
     "HG00097,HGDP00003", "y", "GENE5", "10.0"],
]


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestGenerateReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gr_test_")
        self.addCleanup(lambda: __import__("shutil").rmtree(self.tmp, ignore_errors=True))
        # write the fixture under a canonical finalized filename so filename
        # fallback decodes guide/pam/genome/datasets/mm/bmax
        self.tsv = os.path.join(
            self.tmp,
            f"{_GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv",
        )
        with open(self.tsv, "w") as handle:
            handle.write("\t".join(_HEADER) + "\n")
            for row in _ROWS:
                handle.write("\t".join(row) + "\n")
        # a samplesID dir so the superpop plot resolves EUR/CSA
        self.sid_dir = os.path.join(self.tmp, "samplesIDs")
        os.makedirs(self.sid_dir)
        with open(os.path.join(self.sid_dir, "samplesIDs.1000G.txt"), "w") as h:
            h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
            h.write("HG00096\tGBR\tEUR\tmale\n")
            h.write("HG00097\tGBR\tEUR\tfemale\n")
            h.write("NA18525\tCHB\tEAS\tfemale\n")
        with open(os.path.join(self.sid_dir, "samplesIDs.HGDP.txt"), "w") as h:
            h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
            h.write("HGDP00001\tBrahui\tCSA\tmale\n")
            h.write("HGDP00003\tBrahui\tCSA\tmale\n")

    def _build_and_extract(self, **kwargs):
        out_zip = gr.build_report(
            integrated_tsv=self.tsv,
            out_zip=os.path.join(self.tmp, "fixture_report.zip"),
            samplesid_dir=self.sid_dir,
            **kwargs,
        )
        self.assertTrue(os.path.isfile(out_zip))
        extract = os.path.join(self.tmp, "extracted")
        with zipfile.ZipFile(out_zip) as zf:
            names = zf.namelist()
            zf.extractall(extract)
        return out_zip, names, extract

    @staticmethod
    def _read(path):
        with open(path, encoding="utf-8") as handle:
            return handle.read()

    def test_zip_has_exactly_the_flat_bundle(self):
        _, names, _ = self._build_and_extract()
        # report v2.2 bundles the worst-case panel (panel_top100.tsv) AND one
        # ready-to-go TSV per NON-EMPTY validation tier. The fixture's 5
        # off-targets yield: cfd_ge_{0.50,0.20,0.10,0.05}, mmb_le_{2,3,4} (NOT
        # mmb_le_1 -- 0 rows -> skipped), variant_created. All fixture tiers are
        # tiny -> plain .tsv (no .gz).
        self.assertEqual(
            sorted(names),
            [
                "cfd_ge_0.05.tsv",
                "cfd_ge_0.10.tsv",
                "cfd_ge_0.20.tsv",
                "cfd_ge_0.50.tsv",
                "integrated_results.tsv.gz",
                "mmb_le_2.tsv",
                "mmb_le_3.tsv",
                "mmb_le_4.tsv",
                "panel_top100.tsv",
                "report.html",
                "top1000.tsv",
                "variant_created.tsv",
            ],
        )
        # the 0-row tier (mm+b<=1 after the mm+b>1 base filter) is NOT bundled
        self.assertNotIn("mmb_le_1.tsv", names)
        # flat (no directory components)
        for name in names:
            self.assertNotIn("/", name)

    def test_bundled_tsv_gz_round_trips(self):
        _, _, extract = self._build_and_extract()
        gz = os.path.join(extract, "integrated_results.tsv.gz")
        self.assertTrue(os.path.isfile(gz))
        with gzip.open(gz, "rt") as handle:
            content = handle.read()
        self.assertIn("Spacer+PAM", content.splitlines()[0])
        self.assertEqual(len(content.strip().splitlines()), len(_ROWS) + 1)

    def test_top1000_tsv_bundled_with_table_rows(self):
        _, _, extract = self._build_and_extract()
        top = os.path.join(extract, "top1000.tsv")
        self.assertTrue(os.path.isfile(top))
        with open(top) as handle:
            lines = handle.read().strip().splitlines()
        # header + the 5 off-targets (on-target mm+b==0 filtered out)
        self.assertIn("Spacer+PAM", lines[0])
        self.assertEqual(len(lines), 5 + 1)

    def test_consistency_partition_sums_to_total(self):
        # variant-created + reference + on-target == total EXACTLY
        import pandas as pd

        df = pd.read_csv(self.tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))
        variant, reference, ontarget = gr.partition_masks(df, cols)
        n_v = int(variant.sum())
        n_r = int(reference.sum())
        n_o = int(ontarget.sum())
        self.assertEqual(n_v + n_r + n_o, len(df))
        # the fixture: 3 variant-created, 1 on-target, 2 reference
        self.assertEqual(n_v, 3)
        self.assertEqual(n_o, 1)
        self.assertEqual(n_r, 2)
        # masks are mutually exclusive
        self.assertEqual(int((variant & reference).sum()), 0)
        self.assertEqual(int((variant & ontarget).sum()), 0)
        self.assertEqual(int((reference & ontarget).sum()), 0)

    def test_section1_summary_card_and_global_matrix(self):
        _, _, extract = self._build_and_extract(
            params_override={
                "Nuclease": "SpCas9", "Pam": "NRG", "Genome_selected": "hg38",
                "Mismatches": "6", "DNA": "2", "RNA": "2", "Max_bulges": "2",
                "Max_total_edits": "4",
            }
        )
        html = self._read(os.path.join(extract, "report.html"))
        # header card
        self.assertIn("gRNA (spacer+PAM)", html)
        self.assertIn(_GUIDE, html)
        self.assertIn("Aggregated Specificity Score (0-100)", html)
        self.assertRegex(html, r"Nuclease</td><td>SpCas9<")
        # global MM/B matrix present, grouped REFERENCE vs VARIANT
        self.assertIn("Off-targets by Mismatch (MM) and Bulge (B)", html)
        self.assertIn('class="matrix"', html)
        self.assertIn(">REFERENCE", html)
        self.assertIn(">VARIANT", html)
        # per-mismatch columns up to the run's mm (6MM), on-target excluded
        self.assertIn(">6MM<", html)

    def test_section1_matrix_reconciles_to_grand_total(self):
        """Matrix REFERENCE+VARIANT totals == off-target partition total, and
        REFERENCE + VARIANT + on-target == grand total (every off-target lands
        in a cell). Bulge rows span 0..(bDNA+bRNA)=0..4, mm cols 0..6."""
        import pandas as pd

        df = pd.read_csv(self.tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))
        meta = gr.build_summary_meta(
            None, self.tsv, df, cols,
            params_override={
                "Nuclease": "SpCas9", "Pam": "NRG", "Genome_selected": "hg38",
                "Mismatches": "6", "DNA": "2", "RNA": "2", "Max_bulges": "2",
                "Max_total_edits": "4",
            },
        )
        matrix = gr.build_mmb_matrix(df, cols, meta)

        # bulge rows span 0..(bDNA+bRNA) == 0..4 (NOT Max_bulges=2)
        for label, mrows in matrix["groups"]:
            self.assertEqual([r[0] for r in mrows], [0, 1, 2, 3, 4])
        # mm columns span 0..6
        self.assertEqual(matrix["mm_cols"], [0, 1, 2, 3, 4, 5, 6])

        # matrix REFERENCE + VARIANT totals == canonical off-target totals
        variant, reference, ontarget = gr.partition_masks(df, cols)
        by_label = {lbl: sum(r[1] for r in rows) for lbl, rows in matrix["groups"]}
        self.assertEqual(by_label["REFERENCE"], int(reference.sum()))
        self.assertEqual(by_label["VARIANT"], int(variant.sum()))

        # every off-target lands in a cell: sum of ALL per-mm cells == off-target
        # count; and REFERENCE + VARIANT + on-target == grand total
        cell_sum = 0
        for _lbl, rows in matrix["groups"]:
            for _b, _tot, per_mm in rows:
                cell_sum += sum(per_mm)
        n_offtarget = int((~ontarget).sum())
        self.assertEqual(cell_sum, n_offtarget)
        self.assertEqual(
            by_label["REFERENCE"] + by_label["VARIANT"] + int(ontarget.sum()),
            len(df),
        )

    def test_section2_four_scatter_panels_when_crista_present(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        imgs = re.findall(r'src="data:image/png;base64,([A-Za-z0-9+/=]+)"', html)
        # 4 scatter panels + 1 population plot => 5 inline PNGs; each decodable
        self.assertEqual(len(imgs), 5)
        for encoded in imgs:
            raw = base64.b64decode(encoded)
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")
        # exactly the FOUR panel titles are present because CRISTA is computed
        self.assertIn("By CFD score", html)
        self.assertIn("By variant effect (CFD ALT - REF delta)", html)
        self.assertIn("By CRISTA score", html)
        self.assertIn("By variant effect (CRISTA)", html)

        # and the count is driven by crista_computed(): count scatter <h3> panels
        import pandas as pd

        df = pd.read_csv(self.tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))
        self.assertTrue(gr.crista_computed(df, cols))
        panels = gr.plot_scatter_panels(df, cols, n=1000, include_crista=True)
        self.assertEqual(len(panels), 4)

    def test_section2_two_scatter_panels_when_crista_absent(self):
        # blank out every CRISTA score -> crista_computed() False -> 2 panels
        import pandas as pd

        df = pd.read_csv(self.tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))
        df[cols["crista"]] = ""
        self.assertFalse(gr.crista_computed(df, cols))
        panels = gr.plot_scatter_panels(df, cols, n=1000, include_crista=False)
        self.assertEqual(len(panels), 2)
        titles = [p[0] for p in panels]
        self.assertEqual(
            titles,
            ["By CFD score", "By variant effect (CFD ALT - REF delta)"],
        )

    def test_section4_validation_panel_worstcase(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn("Recommended validation panel", html)
        # worst-case wording + cap (100) present
        self.assertIn("worst-case", html.lower())
        self.assertIn("worst-case top 100", html)
        # the FULL threshold table is kept (CFD>= {0.5,0.2,0.1,0.05},
        # mm+b<= {1,2,3,4})
        for t in ("0.5", "0.2", "0.1", "0.05"):
            self.assertIn(f"CFD &ge; {t}", html)
        for t in ("1", "2", "3", "4"):
            self.assertIn(f"mismatches + bulges &le; {t}", html)
        # of the selected panel, 3 are variant-created (the fixture has exactly 3)
        self.assertRegex(html, r"<strong>3</strong>[^<]*are\s+variant-created")

    def test_section4_worstcase_selection_and_cap(self):
        """A site worst by ANY single metric is included; cap honored."""
        import pandas as pd

        # constants are module-level
        self.assertEqual(gr.PANEL_WORSTCASE_CAP, 100)
        self.assertEqual(
            gr.PANEL_WORSTCASE_METRICS,
            (("cfd", "desc"), ("crista", "desc"), ("mmb", "asc")),
        )

        header = list(_HEADER)
        # Build a synthetic frame where each metric's single "worst" site is
        # NOT the leader on the other metrics -> only a worst-case (min-rank)
        # selection surfaces all three; a plain top-by-CFD would miss two.
        # Columns needed by select_worstcase_panel: mm, bulges, mmb, cfd, crista,
        # not_in_ref. We keep the on-target(mm+b==0) out of the panel.
        def _row(chrom, mm, b, cfd, crista, notref):
            r = ["G", chrom, "1", "+", "AAA", "AAA", "GGG", str(mm), str(b),
                 str(mm + b), notref and "alt" or "ref", "NA",
                 f"{cfd}", f"{cfd}", f"{cfd}", "NA", "NA", "NA", "NA",
                 f"{crista}", f"{crista}", f"{crista}", "NA", "NA", "NA",
                 "y" if notref else "NA", "GENE", "1.0"]
            return r

        rows = [
            _row("chrON", 0, 0, 1.0, 1.0, False),   # on-target -> excluded
            _row("chrCFD", 3, 2, 0.99, 0.10, True),  # worst by CFD only
            _row("chrCRI", 4, 3, 0.10, 0.99, True),  # worst by CRISTA only
            _row("chrMMB", 2, 0, 0.10, 0.10, False),  # worst (closest) by mm+b
            _row("chrLOW", 6, 4, 0.05, 0.05, False),  # low by every metric
        ]
        tsv = os.path.join(self.tmp, "wc.tsv")
        with open(tsv, "w") as h:
            h.write("\t".join(header) + "\n")
            for r in rows:
                h.write("\t".join(r) + "\n")
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))

        panel = gr.select_worstcase_panel(df, cols, cap=gr.PANEL_WORSTCASE_CAP)
        chroms = set(panel[cols["chrom"]])
        # on-target excluded
        self.assertNotIn("chrON", chroms)
        # each single-metric-worst site is surfaced by the worst-case rule
        self.assertIn("chrCFD", chroms)
        self.assertIn("chrCRI", chroms)
        self.assertIn("chrMMB", chroms)

        # top-3 by severity are exactly the three single-metric leaders (in some
        # order); the all-low site ranks last
        top3 = list(panel[cols["chrom"]].head(3))
        self.assertEqual({"chrCFD", "chrCRI", "chrMMB"}, set(top3))
        self.assertEqual(list(panel[cols["chrom"]])[-1], "chrLOW")

        # cap is honored: many low sites -> panel never exceeds the cap
        big = [_row("chrON", 0, 0, 1.0, 1.0, False)]
        for i in range(500):
            big.append(_row(f"c{i}", 5, 3, 0.20, 0.20, i % 2 == 0))
        tsv2 = os.path.join(self.tmp, "wc_big.tsv")
        with open(tsv2, "w") as h:
            h.write("\t".join(header) + "\n")
            for r in big:
                h.write("\t".join(r) + "\n")
        df2 = pd.read_csv(tsv2, sep="\t", dtype=str, na_filter=False)
        cols2 = gr._resolve(df2.columns, list(gr._COLS.keys()))
        panel2 = gr.select_worstcase_panel(df2, cols2, cap=gr.PANEL_WORSTCASE_CAP)
        self.assertEqual(len(panel2), gr.PANEL_WORSTCASE_CAP)
        # on-target never in the panel
        self.assertNotIn("chrON", set(panel2[cols2["chrom"]]))

    def test_section45_per_tier_files_bundled_and_linked_relative(self):
        """v2.2: each NON-EMPTY tier is bundled as its own ready-to-use TSV,
        column-identical to top1000.tsv, and linked RELATIVELY in the section-4
        Download column AND section-5; the 0-row tier is skipped (no file/link).
        """
        _, names, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))

        # the expected non-empty tier files for the fixture (mmb_le_1 skipped)
        expected = [
            "cfd_ge_0.50.tsv", "cfd_ge_0.20.tsv", "cfd_ge_0.10.tsv",
            "cfd_ge_0.05.tsv", "mmb_le_2.tsv", "mmb_le_3.tsv", "mmb_le_4.tsv",
            "variant_created.tsv",
        ]
        for name in expected:
            with self.subTest(tier=name):
                # bundled in the zip
                self.assertIn(name, names)
                self.assertTrue(os.path.isfile(os.path.join(extract, name)))
                # linked RELATIVELY (bare filename href, no scheme / no leading /)
                self.assertIn(f'href="{name}"', html)
                self.assertNotIn(f'href="/{name}"', html)
                self.assertNotIn(f"http://{name}", html)

        # the worst-case top-100 panel is ALSO surfaced in the tier table
        self.assertIn('href="panel_top100.tsv"', html)

        # 0-row tier (mm+b<=1) is skipped: no file, no link anywhere
        self.assertNotIn("mmb_le_1.tsv", names)
        self.assertNotIn("mmb_le_1", html)

        # a new Download column exists in the section-4 threshold tables
        self.assertIn("<th>Download</th>", html)

        # tier files are column-identical to top1000.tsv, sorted by CFD desc,
        # off-targets only (mm+b<=1 excluded)
        top_hdr = self._read(os.path.join(extract, "top1000.tsv")).splitlines()[0]
        for name in expected:
            with self.subTest(columns_of=name):
                lines = self._read(os.path.join(extract, name)).splitlines()
                self.assertEqual(lines[0], top_hdr)  # same column set/order
                # every data row has mm+b > 1 (on-/near-on-target excluded)
                import pandas as pd

                sub = pd.read_csv(
                    os.path.join(extract, name), sep="\t", dtype=str, na_filter=False
                )
                mmb = pd.to_numeric(
                    sub["Mismatches+bulges_(highest_CFD)"], errors="coerce"
                )
                self.assertTrue((mmb > 1).all())
                cfd = pd.to_numeric(
                    sub["CFD_score_(highest_CFD)"], errors="coerce"
                ).fillna(-1.0).tolist()
                self.assertEqual(cfd, sorted(cfd, reverse=True))  # CFD desc

        # exact per-tier row counts for the fixture (5 off-targets)
        counts = {
            "cfd_ge_0.50.tsv": 3, "cfd_ge_0.20.tsv": 5, "cfd_ge_0.10.tsv": 5,
            "cfd_ge_0.05.tsv": 5, "mmb_le_2.tsv": 2, "mmb_le_3.tsv": 4,
            "mmb_le_4.tsv": 5, "variant_created.tsv": 3, "panel_top100.tsv": 5,
        }
        for name, n in counts.items():
            with self.subTest(rowcount_of=name):
                lines = self._read(os.path.join(extract, name)).strip().splitlines()
                self.assertEqual(len(lines), n + 1)  # + header

    def test_tier_specs_are_module_level_and_editable(self):
        """The tier list is a module-level structure (easy to edit)."""
        self.assertTrue(hasattr(gr, "TIER_SPECS"))
        keys = [s["key"] for s in gr.TIER_SPECS]
        self.assertEqual(
            keys,
            [
                "cfd_ge_0.50", "cfd_ge_0.20", "cfd_ge_0.10", "cfd_ge_0.05",
                "mmb_le_1", "mmb_le_2", "mmb_le_3", "mmb_le_4",
                "variant_created",
            ],
        )
        # gzip threshold is a module-level constant (~2 MB)
        self.assertEqual(gr.TIER_GZIP_MAX_BYTES, 2 * 1024 * 1024)

    def test_large_tier_is_gzipped_and_labeled(self):
        """A tier whose plain TSV exceeds the ~2 MB threshold is gzipped, the
        bundled name gains .gz, and the link label reflects it (.tsv vs .tsv.gz).
        """
        import pandas as pd

        df = pd.read_csv(self.tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))
        frame = gr._offtarget_cfd_sorted(df, cols)

        staging = tempfile.mkdtemp(prefix="gr_tier_")
        self.addCleanup(
            lambda: __import__("shutil").rmtree(staging, ignore_errors=True)
        )
        # tiny threshold -> even the small fixture tier gzips
        name, path = gr.write_tier_tsv(
            frame, staging, "cfd_ge_0.05.tsv", gzip_max=10
        )
        self.assertEqual(name, "cfd_ge_0.05.tsv.gz")
        self.assertTrue(os.path.isfile(path))
        self.assertFalse(os.path.isfile(os.path.join(staging, "cfd_ge_0.05.tsv")))
        # round-trips + same columns
        with gzip.open(path, "rt") as h:
            content = h.read()
        self.assertIn("Spacer+PAM", content.splitlines()[0])

        # a large threshold keeps it plain
        name2, path2 = gr.write_tier_tsv(
            frame, staging, "cfd_ge_0.20.tsv", gzip_max=10 ** 9
        )
        self.assertEqual(name2, "cfd_ge_0.20.tsv")
        self.assertTrue(path2.endswith(".tsv"))

    def test_maf_footnote_present_under_summary_and_table(self):
        """v2.2: the concise Variant_MAF footnote (blank/em-dash explanation) is
        present, with the SNP-only-registry + AC/AN wording, and the existing
        em-dash rendering + rsID->genomic-key fallback are intact.
        """
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # the footnote text is a module-level constant, easy to audit/edit
        self.assertTrue(hasattr(gr, "MAF_FOOTNOTE"))
        for phrase in (
            "reference off-target (no variant)",
            "the allele-frequency registry is SNP-only",
            "a SNP variant not present in the frequency panel",
            "AC/AN over the genotyped panel",
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, gr.MAF_FOOTNOTE)
                self.assertIn(phrase, html)
        # rendered twice (under the summary AND under the top-1000 table)
        self.assertGreaterEqual(html.count("maf-footnote"), 2)
        # the em-dash blank-MAF rendering is still there (fixture has a blank MAF)
        self.assertIn("&mdash;", html)
        # rsID->genomic-key fallback intact: the multi-SNP row shows its rsID
        self.assertIn("rs333", html)

    def test_section6_table_has_pam_creation_and_crista(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn('class="ottable"', html)
        # PAM-creation column header + CRISTA column header (CRISTA computed)
        self.assertIn("<th>PAM creation</th>", html)
        self.assertIn("<th>CRISTA</th>", html)
        # rows sorted by CFD desc, on-target excluded => 5 rows
        tbody = html.split('class="ottable"')[-1].split("<tbody>")[-1].split("</tbody>")[0]
        n_rows = tbody.count("<tr>")
        self.assertEqual(n_rows, 5)
        cfds = re.findall(r"<td>(0\.\d{4})</td>", tbody)
        self.assertEqual(cfds, sorted(cfds, reverse=True))
        # the pam_created value from the fixture appears
        self.assertIn("pam_created", html)

    def test_section7_footer_disclaimer_and_provenance(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # exact disclaimer text (research-purposes-only wording)
        self.assertIn(
            "This report is provided for research purposes only.", html
        )
        self.assertIn(
            "are NOT a substitute for experimental validation", html
        )
        # provenance stamp: report generator version + source TSV basename
        self.assertIn("report generator v", html)
        self.assertIn(os.path.basename(self.tsv), html)

    def test_self_contained_offline_and_relative_links(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # NO external dependencies -> opens offline with file://
        self.assertNotIn("<script", html.lower())
        self.assertNotIn("<link", html.lower())
        self.assertNotIn("http://", html)
        self.assertNotIn("https://", html)
        # RELATIVE download links to both bundled siblings
        self.assertIn('href="integrated_results.tsv.gz"', html)
        self.assertIn('href="top1000.tsv"', html)

    def test_multi_snp_haplotype_min_maf_and_first_rsid(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # the multi-SNP row (rs333,rs444 / NA,0.005,0.003) must show rs333 and
        # the MIN maf 0.003 -> "3.00e-03"
        self.assertIn("rs333", html)
        self.assertIn("3.00e-03", html)


if __name__ == "__main__":
    unittest.main(verbosity=2)

#!/usr/bin/env python3
"""Unit test for the shareable-report generator (PostProcess/generate_report.py).

Builds a tiny integrated_results fixture (dict-based schema with CRISTA + PAM
creation columns, header names only for the columns the report reads) and
asserts that build_report produces a ZIP that flat-decompresses to exactly three
files -- report.html + integrated_results.tsv.gz + top1000.tsv -- and that the
HTML is a self-contained IND-briefing-book digest (report v2):

  * SECTION 1: the header summary card (guide, PAM, aggregated specificity score)
    AND the global Off-targets-by-MM-and-B matrix are present,
  * the canonical partition invariant holds:
    variant-created + reference + on-target == total EXACTLY,
  * SECTION 2: >= 2 CFD scatter panels are embedded as inline base64 PNGs
    (by CFD + by delta), plus a CRISTA panel because CRISTA is non-empty,
  * SECTION 5: the top-N off-target table is rendered inline (rows sorted by CFD)
    with a PAM-creation column, and top1000.tsv is bundled in the zip,
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

    def test_zip_has_exactly_three_flat_files(self):
        _, names, _ = self._build_and_extract()
        self.assertEqual(
            sorted(names),
            ["integrated_results.tsv.gz", "report.html", "top1000.tsv"],
        )
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

    def test_section2_scatter_panels_and_crista(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        imgs = re.findall(r'src="data:image/png;base64,([A-Za-z0-9+/=]+)"', html)
        # >= 2 CFD scatter panels + population plot + CRISTA panel => >= 3 images;
        # each must be a decodable PNG
        self.assertGreaterEqual(len(imgs), 2)
        for encoded in imgs:
            raw = base64.b64decode(encoded)
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")
        # both required scatter sort orders are present
        self.assertIn("By CFD score", html)
        self.assertIn("By variant effect (ALT - REF CFD)", html)
        # CRISTA panel included because CRISTA is non-empty in the fixture
        self.assertIn("By CRISTA score", html)

    def test_section4_validation_panel(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn("Recommended validation panel", html)
        self.assertIn("Suggested tiered panel", html)
        # threshold labels
        self.assertIn("CFD &ge; 0.5", html)
        self.assertIn("mismatches + bulges &le; 2", html)

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

#!/usr/bin/env python3
"""Unit test for the shareable-report generator (PostProcess/generate_report.py).

Builds a tiny integrated_results fixture (dict-less 85-col schema, header names
only for the columns the report reads) and asserts that build_report produces a
ZIP that flat-decompresses to exactly two files -- report.html +
integrated_results.tsv.gz -- and that the HTML is a self-contained digest:

  * the top-of-page summary table (guide, PAM, counts) is present,
  * all three plots are embedded as inline base64 PNG data URIs,
  * the top-N off-target table is rendered inline (rows sorted by CFD),
  * there is no external dependency (<script>/<link>/http(s)), so it opens with
    file:// offline,
  * the download link is RELATIVE (resolves after unzip on any machine).

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
# be excluded from the top-N by the mm+b > 1 filter. Only the header names the
# report resolves are meaningful; unused columns are omitted on purpose to prove
# the by-name resolution / graceful-degradation behaviour.
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
    "CFD_score_(highest_CFD)",
    "Variant_info_genome_(highest_CFD)",
    "Variant_MAF_(highest_CFD)",
    "Variant_rsID_(highest_CFD)",
    "Variant_samples_(highest_CFD)",
    "Not_found_in_REF",
    "Annotation_closest_gene_name",
    "Annotation_closest_gene_distance_(kb)",
]

_GUIDE = "CTCTCAGCTGGTACACGGCANNN"

# columns: guide, chrom, pos, strand, aln_ref, aln_alt, pam, mm, bul, mmb,
#          origin, cfd, var_genome, maf, rsid, samples, not_in_ref, gene, dist
_ROWS = [
    # on-target (mm+b == 0) -> excluded from top-N; counted as 1 on-target
    [_GUIDE, "chr2", "100", "+", "CTCTCAGCTGGTACACGGCATGG", "NA", "TGG",
     "0", "0", "0", "ref", "1.0", "NA", "NA", "NA", "NA", "NA", "GENE0", "0.0"],
    # reference off-targets
    [_GUIDE, "chr3", "200", "-", "cTCTCAGCTGGTACACGGCAAGG", "NA", "AGG",
     "2", "0", "2", "ref", "0.85", "NA", "NA", "NA", "NA", "NA", "GENE1", "1.2"],
    [_GUIDE, "chr4", "300", "+", "ctCTCAGCTGGTACACGGCAcGG", "NA", "CGG",
     "3", "0", "3", "ref", "0.40", "NA", "NA", "NA", "NA", "NA", "GENE2", "5.0"],
    # variant-created off-targets (1000G + HGDP carriers)
    [_GUIDE, "chr5", "400", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAAGG", "AGG",
     "2", "0", "2", "alt", "0.90", "chr5_400_T_A", "0.01", "rs111",
     "HG00096,HGDP00001", "y", "GENE3", "0.5"],
    [_GUIDE, "chr6", "500", "-", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAcGG", "CGG",
     "3", "1", "4", "alt", "0.55", "chr6_500_T_C", "0.02", "rs222",
     "NA18525", "y", "GENE4", "2.0"],
    # multi-SNP haplotype: comma-joined rsID/MAF/samples (min-AF + first-rsID)
    [_GUIDE, "chr7", "600", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAgGG", "GGG",
     "3", "0", "3", "alt", "0.30", "chr7_600_T_G,chr7_601_A_C",
     "NA,0.005,0.003", "NA,rs333,rs444", "HG00097,HGDP00003", "y",
     "GENE5", "10.0"],
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

    def test_zip_has_exactly_two_flat_files(self):
        _, names, _ = self._build_and_extract()
        self.assertEqual(
            sorted(names), ["integrated_results.tsv.gz", "report.html"]
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

    def test_html_is_self_contained_digest(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))

        # summary table present with the expected key fields
        self.assertIn('class="summary-table"', html)
        self.assertIn("Guide RNA (spacer+PAM)", html)
        self.assertIn(_GUIDE, html)
        self.assertIn("Total off-targets found", html)
        # 6 rows - 1 on-target = 5 off-targets
        self.assertRegex(html, r"Total off-targets found</td><td>5\b")
        self.assertRegex(html, r"On-target site\(s\)</td><td>1\b")
        # 3 variant-created rows (Not_found_in_REF == y)
        self.assertRegex(html, r"Variant-created off-targets</td><td>3\b")

        # at least one embedded base64 PNG (plots inline, not linked)
        imgs = re.findall(r'src="data:image/png;base64,([A-Za-z0-9+/=]+)"', html)
        self.assertGreaterEqual(len(imgs), 1)
        # and each is a decodable PNG
        for encoded in imgs:
            raw = base64.b64decode(encoded)
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")

        # top-N table present, rows sorted by CFD desc, on-target excluded
        self.assertIn('class="ottable"', html)
        tbody = html.split("<tbody>")[-1].split("</tbody>")[0]
        n_rows = tbody.count("<tr>")
        self.assertEqual(n_rows, 5)  # 5 off-targets, on-target filtered out
        # CFD column order: 0.90, 0.85, 0.55, 0.40, 0.30 (desc)
        cfds = re.findall(r"<td>(0\.\d{4})</td>", tbody)
        self.assertEqual(cfds, sorted(cfds, reverse=True))

        # NO external dependencies -> opens offline with file://
        self.assertNotIn("<script", html.lower())
        self.assertNotIn("<link", html.lower())
        self.assertNotIn("http://", html)
        self.assertNotIn("https://", html)

        # RELATIVE download link to the sibling gz
        self.assertIn('href="integrated_results.tsv.gz"', html)

    def test_params_override_populates_summary(self):
        _, _, extract = self._build_and_extract(
            params_override={
                "Nuclease": "SpCas9",
                "DNA": "2",
                "RNA": "2",
                "Max_total_edits": "4",
            }
        )
        html = self._read(os.path.join(extract, "report.html"))
        self.assertRegex(html, r"Nuclease</td><td>SpCas9<")
        self.assertRegex(html, r"Max total edits</td><td>4<")
        self.assertRegex(html, r"Bulges \(DNA / RNA\)</td><td>2 / 2")

    def test_multi_snp_haplotype_min_maf_and_first_rsid(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # the multi-SNP row (rs333,rs444 / NA,0.005,0.003) must show rs333 and
        # the MIN maf 0.003 -> "3.00e-03"
        self.assertIn("rs333", html)
        self.assertIn("3.00e-03", html)


if __name__ == "__main__":
    unittest.main(verbosity=2)

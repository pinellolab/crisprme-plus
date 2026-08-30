#!/usr/bin/env python3
"""Unit test for the shareable-report generator (PostProcess/generate_report.py).

Builds a tiny integrated_results fixture (dict-based schema with CRISTA + PAM
creation + annotation columns GENCODE/gene/distance/ENCODE/DHS, header names only
for the columns the report reads) and asserts that build_report produces a ZIP
that flat-decompresses to report.html + the RAW integrated_results.tsv.gz + the
CURATED top1000.tsv + panel_top100.tsv + the per-tier curated TSVs, and that the
HTML is a self-contained IND-briefing-book digest (report v2.4):

  * SECTION 1: the header summary card (guide, PAM, aggregated specificity score)
    AND the global Off-targets-by-MM-and-B matrix are present; the matrix
    REFERENCE + VARIANT totals reconcile to the canonical partition off-target
    total (REFERENCE + VARIANT + on-target == grand total),
  * the canonical partition invariant holds:
    variant-created + reference + on-target == total EXACTLY,
  * SECTION 2: exactly FOUR ref/alt scatter panels are embedded as inline base64
    SVGs when CRISTA is computed (by CFD, by CFD delta, by CRISTA, by CRISTA
    delta); 2 panels when CRISTA is absent,
  * CURATED COLUMNS (report v2.4): ONE curated column set is shared by BOTH the
    in-report top-1000 table AND every download file (top1000.tsv,
    panel_top100.tsv, per-tier tsvs). The table shows the annotation columns
    (Gene, Gene_distance_kb, GENCODE, ENCODE, DHS); the downloads use the SAME
    curated columns (header spot-checked on panel_top100.tsv),
  * SECTION 4: the HYBRID panel -- hard-include mm+b <= PANEL_FLOOR_MMB (2) OR
    CFD >= PANEL_FLOOR_CFD (0.5), fill to PANEL_CAP (100) by worst-case severity
    (no variant quota); the module-level floors are present, and the in-report
    methods note is present,
  * PER-TIER DOWNLOADS: the non-empty CFD/mm+b/variant-created subsets are bundled
    as curated TSVs and linked,
  * MAF FOOTNOTE: the v2.4 blank-MAF explanation is present,
  * SECTION 5/6: top1000.tsv + panel_top100.tsv are bundled in the zip,
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
    "Annotation_GENCODE",
    "Annotation_closest_gene_name",
    "Annotation_closest_gene_distance_(kb)",
    "Annotation_ENCODE",
    "Annotation_DHS",
    "Annotation_COSMIC",
]

_GUIDE = "CTCTCAGCTGGTACACGGCANNN"

# columns: guide, chrom, pos, strand, aln_ref, aln_alt, pam, mm, bul, mmb,
#          origin, pam_creation, cfd, cfd_ref, cfd_alt, var_genome, maf, rsid,
#          samples, crista, crista_ref, crista_alt, crista_maf, crista_rsid,
#          crista_samples, not_in_ref, GENCODE, gene, dist, ENCODE, DHS, COSMIC
_ROWS = [
    # on-target (mm+b == 0) -> excluded from top-N; counted as 1 on-target
    [_GUIDE, "chr2", "100", "+", "CTCTCAGCTGGTACACGGCATGG", "NA", "TGG",
     "0", "0", "0", "ref", "NA", "1.0", "1.0", "1.0", "NA", "NA", "NA", "NA",
     "1.0", "1.0", "1.0", "NA", "NA", "NA", "NA",
     "protein_coding", "GENE0", "0.0", "CTCF", "DHS_1", "-"],
    # reference off-targets
    [_GUIDE, "chr3", "200", "-", "cTCTCAGCTGGTACACGGCAAGG", "NA", "AGG",
     "2", "0", "2", "ref", "NA", "0.85", "0.85", "0.85", "NA", "NA", "NA", "NA",
     "0.80", "0.80", "0.80", "NA", "NA", "NA", "NA",
     "intron", "GENE1", "1.2", "enhancer", "DHS_2", "TSG"],
    [_GUIDE, "chr4", "300", "+", "ctCTCAGCTGGTACACGGCAcGG", "NA", "CGG",
     "3", "0", "3", "ref", "NA", "0.40", "0.40", "0.40", "NA", "NA", "NA", "NA",
     "0.35", "0.35", "0.35", "NA", "NA", "NA", "NA",
     "NA", "GENE2", "5.0", "NA", "NA", "-"],
    # variant-created off-targets (1000G + HGDP carriers)
    [_GUIDE, "chr5", "400", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAAGG", "AGG",
     "2", "0", "2", "alt", "NA", "0.90", "0.20", "0.90", "chr5_400_T_A",
     "0.01", "rs111", "HG00096,HGDP00001",
     "0.88", "0.18", "0.88", "0.01", "rs111", "HG00096,HGDP00001",
     "y", "exon", "GENE3", "0.5", "promoter", "DHS_3", "oncogene"],
    [_GUIDE, "chr6", "500", "-", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAcGG", "CGG",
     "3", "1", "4", "alt", "pam_created", "0.55", "0.30", "0.55", "chr6_500_T_C",
     "0.02", "rs222", "NA18525",
     "0.50", "0.25", "0.50", "0.02", "rs222", "NA18525",
     "y", "intergenic", "GENE4", "2.0", "CTCF", "DHS_4", "-"],
    # multi-SNP haplotype: comma-joined rsID/MAF/samples (min-AF + first-rsID),
    # and a BLANK MAF -> em-dash + footnote path
    [_GUIDE, "chr7", "600", "+", "CTCTCAGCTGGTACACGGCAtGG",
     "CTCTCAGCTGGTACACGGCAgGG", "GGG",
     "3", "0", "3", "alt", "NA", "0.30", "0.10", "0.30", "chr7_600_T_G,chr7_601_A_C",
     "NA,0.005,0.003", "NA,rs333,rs444", "HG00097,HGDP00003",
     "0.28", "0.08", "0.28", "NA,0.005,0.003", "NA,rs333,rs444",
     "HG00097,HGDP00003", "y", "lincRNA", "GENE5", "10.0", "enhancer", "DHS_5", "-"],
]


@unittest.skipUnless(_HAVE_DEPS, _SKIP_REASON)
class TestGenerateReport(unittest.TestCase):
    def setUp(self):
        # reset build_report's module-level run state so a prior test's run
        # (e.g. a drop_maf=True or partial-annotation build) cannot leak into a
        # standalone curated_headers()/_active_columns() call in this test
        gr._DROP_MAF = False
        gr._PRESENT_ANN_KINDS = None
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

    @staticmethod
    def _dpath(extract, name):
        """Path to a bundled file: report.html at the TOP level, everything else
        under the data/ subfolder (report v2.4 layout)."""
        if name == "report.html":
            return os.path.join(extract, name)
        return os.path.join(extract, "data", name)

    def test_zip_layout_report_at_root_rest_under_data(self):
        _, names, _ = self._build_and_extract()
        # ONLY report.html at the top level; every other bundled file under data/.
        # For this fixture (5 off-targets):
        #   CFD>= 0.5/0.2/0.05  -> all non-empty (0.1 tier removed)
        #   CRISTA>= 0.6/0.4/0.2 -> all non-empty (fixture has CRISTA)
        #   mm+b<= 1/2/3/4      -> 0/2/4/5 (mmb_le_1 is EMPTY -> not bundled)
        #   variant_created     -> 3
        self.assertEqual(
            sorted(names),
            sorted([
                "report.html",
                "data/integrated_results.tsv.gz",
                "data/run_manifest.json",  # IND traceability manifest
                "data/top1000.tsv",
                "data/top1000_crista.tsv",  # CRISTA-ranked companion
                "data/panel_top100.tsv",
                "data/cfd_ge_0.50.tsv",
                "data/cfd_ge_0.20.tsv",
                "data/cfd_ge_0.05.tsv",
                "data/crista_ge_0.60.tsv",  # CRISTA-appropriate tiers (not CFD's)
                "data/crista_ge_0.40.tsv",
                "data/crista_ge_0.20.tsv",
                "data/mmb_le_2.tsv",
                "data/mmb_le_3.tsv",
                "data/mmb_le_4.tsv",
                "data/variant_created.tsv",
            ]),
        )
        # report.html is the ONLY top-level entry (self-obvious what to open)
        top_level = [n for n in names if "/" not in n]
        self.assertEqual(top_level, ["report.html"])
        # the empty tier (mm+b <= 1) is NOT bundled
        self.assertNotIn("data/mmb_le_1.tsv", names)

    def test_pertier_downloads_are_curated_and_linked(self):
        _, names, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        expected_curated = gr.curated_headers(has_crista=True)
        # every bundled per-tier / panel / top1000 TSV uses the SAME curated header
        for fname in (
            "top1000.tsv", "panel_top100.tsv",
            "cfd_ge_0.50.tsv", "cfd_ge_0.20.tsv",
            "cfd_ge_0.05.tsv", "mmb_le_2.tsv", "mmb_le_3.tsv", "mmb_le_4.tsv",
            "variant_created.tsv",
        ):
            path = self._dpath(extract, fname)
            self.assertTrue(os.path.isfile(path), fname)
            with open(path) as handle:
                header = handle.readline().rstrip("\n").split("\t")
            self.assertEqual(header, expected_curated, fname)
            # each bundled tier file is linked from the HTML with a RELATIVE
            # href into the data/ subfolder
            self.assertIn(f'href="data/{fname}"', html)
        # curated header carries the annotation columns
        for ann in ("Gene", "Gene_distance_kb", "GENCODE", "ENCODE", "DHS",
                    "COSMIC_cancer_gene"):
            self.assertIn(ann, expected_curated)

    def test_curated_column_list_is_the_shared_set(self):
        # the ONE curated set, with CRISTA when computed
        self.assertEqual(
            gr.curated_headers(has_crista=True),
            [
                "rank", "Chromosome", "Position", "Strand",
                "Aligned_protospacer+PAM", "Mismatches", "Bulges",
                "Mismatches+bulges", "Perfect_match", "CFD", "CRISTA",
                "REF/ALT_origin",
                "PAM_creation", "Variant", "MAF", "Gene", "Gene_distance_kb",
                "GENCODE", "ENCODE", "DHS", "COSMIC_cancer_gene",
                "High_complexity_region",
            ],
        )
        # CRISTA is dropped when not computed (else identical)
        no_crista = gr.curated_headers(has_crista=False)
        self.assertNotIn("CRISTA", no_crista)
        self.assertEqual(
            [c for c in gr.curated_headers(has_crista=True) if c != "CRISTA"],
            no_crista,
        )

    def test_high_complexity_region_flag_projection(self):
        # the curated cell compacts the integrated_results note to "Yes (N var)"
        note = (
            "high_variant_density (12 variants): a greedy worst-case alignment is "
            "reported here; additional haplotype alignments may exist; full_IUPAC=ACGT"
        )
        cell = gr._curated_cell(
            "complex_region",
            {"High_variant_density_region": note},
            {"complex_region": "High_variant_density_region"},
        )
        self.assertEqual(cell, "Yes (12 var)")
        # "NA" -> the missing sentinel
        self.assertEqual(
            gr._curated_cell(
                "complex_region",
                {"High_variant_density_region": "NA"},
                {"complex_region": "High_variant_density_region"},
            ),
            gr.CURATED_MISSING,
        )

    def test_highly_complex_regions_bundled_and_surfaced(self):
        # a run with a high_variant_density_regions.bed in the result dir:
        # the report bundles the merged bed, links it, and shows the callout.
        rd = os.path.join(self.tmp, "hvdr_run")
        os.makedirs(rd)
        tsv = os.path.join(
            rd, f"{_GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv"
        )
        header = list(_HEADER) + ["High_variant_density_region"]
        note = (
            "high_variant_density (12 variants): a greedy worst-case alignment is "
            "reported here; additional haplotype alignments may exist; full_IUPAC=ACGTNNNN"
        )
        with open(tsv, "w") as h:
            h.write("\t".join(header) + "\n")
            for i, row in enumerate(_ROWS):
                h.write("\t".join(list(row) + [note if i == 0 else "NA"]) + "\n")
        with open(os.path.join(rd, "job.high_variant_density_regions.bed"), "w") as h:
            h.write(
                "#chrom\tstart\tend\tguide\tn_variants\tsamples_with_alt\tiupac_protospacer\n"
            )
            # overlaps the chr3:200 top-N off-target so it survives the top-N scoping.
            # Written 3x (per-chrom beds carry duplicate region rows) to exercise the
            # dedup in the merge -- the bundle must contain it exactly ONCE.
            _dupe = "chr3\t190\t213\tGUIDE\t12\tHG00096,HG00097\tACGTNNNNACGTACGTACGTNGG\n"
            h.write(_dupe)
            h.write(_dupe)
            h.write(_dupe)
        out_zip = gr.build_report(
            result_dir=rd,
            samplesid_dir=self.sid_dir,
            out_zip=os.path.join(self.tmp, "hvdr_report.zip"),
        )
        extract = os.path.join(self.tmp, "hvdr_extract")
        with zipfile.ZipFile(out_zip) as zf:
            names = zf.namelist()
            zf.extractall(extract)
        self.assertIn("data/high_variant_density_regions.bed", names)
        # the bundled bed keeps the IUPAC field
        bed = self._read(self._dpath(extract, "high_variant_density_regions.bed"))
        self.assertIn("iupac_protospacer", bed)
        self.assertIn("ACGTNNNNACGTACGTACGTNGG", bed)
        # dedup: the 3 identical per-chrom rows collapse to ONE in the bundle
        self.assertEqual(bed.count("chr3\t190\t213\tGUIDE"), 1)
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn("Highly complex (high-variant-density) regions", html)
        self.assertIn('href="data/high_variant_density_regions.bed"', html)
        # curated column is part of the shared set (header of top1000)
        top_header = self._read(
            self._dpath(extract, "top1000.tsv")
        ).splitlines()[0].split("\t")
        self.assertIn("High_complexity_region", top_header)

    def test_indel_snp_cooc_bundled_and_surfaced(self):
        # a run with per-chromosome *.indel_snp_cooc.tsv sidecars: the report
        # merges them (single header), bundles the merged TSV, links it, and shows
        # the confirmed-cis count (cis rows only, not trans/candidate rows).
        rd = os.path.join(self.tmp, "cooc_run")
        os.makedirs(rd)
        tsv = os.path.join(
            rd, f"{_GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv"
        )
        with open(tsv, "w") as h:
            h.write("\t".join(_HEADER) + "\n")
            for row in _ROWS:
                h.write("\t".join(row) + "\n")
        _cooc_header = (
            "chrom\tindel_pos\tindel_ref\tindel_alt\tofftarget_start\tstrand\t"
            "snp_dictpos\tsnp_rsid\tphase\tjoint_af\tn_cis\tcis_samples\n"
        )
        # phase vocabulary is CONFIRMED / PUTATIVE (indel_snp_cis.py) -- every row IS
        # a cis co-occurrence; CONFIRMED = proven phasing. chr3 sidecar: 1 CONFIRMED
        # + 1 PUTATIVE (PUTATIVE must NOT count toward confirmed-cis).
        with open(os.path.join(rd, "job_chr3.indel_snp_cooc.tsv"), "w") as h:
            h.write(_cooc_header)
            h.write("chr3\t190\tAT\tA\t200\t-\tchr3_195_C_T\trs900\tCONFIRMED\t0.0021\t2\tHG00096,HG00097\n")
            h.write("chr3\t250\tG\tGA\t260\t+\tchr3_255_A_G\trs901\tPUTATIVE\t0.0011\t1\tHG00097\n")
        # chr7 sidecar: 1 CONFIRMED row (has its OWN header -> dedup to one)
        with open(os.path.join(rd, "job_chr7.indel_snp_cooc.tsv"), "w") as h:
            h.write(_cooc_header)
            h.write("chr7\t600\tC\tCTT\t600\t+\tchr7_601_A_C\trs444\tCONFIRMED\t0.0033\t1\tHGDP00003\n")
        out_zip = gr.build_report(
            result_dir=rd,
            samplesid_dir=self.sid_dir,
            out_zip=os.path.join(self.tmp, "cooc_report.zip"),
        )
        extract = os.path.join(self.tmp, "cooc_extract")
        with zipfile.ZipFile(out_zip) as zf:
            names = zf.namelist()
            zf.extractall(extract)
        # (1) bundled into the ZIP under data/
        self.assertIn("data/indel_snp_cooc.tsv", names)
        merged = self._read(self._dpath(extract, "indel_snp_cooc.tsv"))
        lines = [ln for ln in merged.splitlines() if ln.strip()]
        # (2) exactly ONE header (dedup across the two per-chrom files) + both chroms
        self.assertEqual(lines[0].split("\t")[0], "chrom")
        self.assertEqual(sum(1 for ln in lines if ln.startswith("chrom\t")), 1)
        self.assertTrue(any(ln.startswith("chr3\t") for ln in lines))
        self.assertTrue(any(ln.startswith("chr7\t") for ln in lines))
        # (3) report.html links the file + names the CONFIRMED-CIS count (2 cis rows,
        #     NOT the 3 total candidate rows: the trans row is excluded)
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn('href="data/indel_snp_cooc.tsv"', html)
        self.assertIn("SNP + indel cis co-occurrences", html)
        self.assertIn("2 confirmed-cis", html)
        self.assertIn("of 3 candidate", html)

    def test_load_sample_dataset_native_provenance(self):
        """sample -> native per-db label, read from the files (no hardcoding)."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "hg38_1000G.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                h.write("HG00096\tGBR\tEUR\tmale\n")
                h.write("NA19017\tLWK\tAFR\tfemale\n")
            with open(os.path.join(d, "hg38_HGDP.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                h.write("HGDP00001\tBrahui\tCSA\tmale\n")
            # the combined superset must NOT override the per-db label
            with open(os.path.join(d, "hg38_1000G_HGDP.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                for s in ("HG00096", "NA19017", "HGDP00001"):
                    h.write(f"{s}\tX\tY\tZ\n")
            m = gr.load_sample_dataset(d)
            self.assertEqual(m["HG00096"], "1000G")
            self.assertEqual(m["NA19017"], "1000G")
            self.assertEqual(m["HGDP00001"], "HGDP")
        # a completely different database is labeled dynamically (never hardcoded)
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "hg38_gnomAD.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOP\tSUPERPOP\tSEX\n")
                h.write("SAMPLE_X\tp\ts\tm\n")
            self.assertEqual(gr.load_sample_dataset(d)["SAMPLE_X"], "gnomAD")
        self.assertEqual(gr.load_sample_dataset("/nonexistent"), {})

    def test_load_sample_superpop_reads_finalized_install_names(self):
        """Superpop must load from the REAL install layout hg38_<db>.samplesID.txt
        (not just the classic samplesIDs.<db>.txt), else the plot silently blanks."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "hg38_1000G.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                h.write("HG00096\tGBR\tEUR\tmale\n")
            with open(os.path.join(d, "hg38_HGDP.samplesID.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                h.write("HGDP00001\tBrahui\tCSA\tmale\n")
            # a config sidecar must be ignored, not parsed as samples
            with open(os.path.join(d, "samplesIDs.config.txt"), "w") as h:
                h.write("hg38\t1000G\tvcf\n")
            m = gr.load_sample_superpop(d)
            self.assertEqual(m.get("HG00096"), "EUR")
            self.assertEqual(m.get("HGDP00001"), "CSA")
            # classic download name also still works
            with open(os.path.join(d, "samplesIDs.EXTRA.txt"), "w") as h:
                h.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n")
                h.write("NA20000\tTSI\tEUR\tmale\n")
            self.assertEqual(gr.load_sample_superpop(d).get("NA20000"), "EUR")

    def test_genome_normalization_and_dataset_poison(self):
        """'_ref' suffix is stripped from the genome; a reference-only vcf token
        does not leak 'ref' as a dataset."""
        self.assertEqual(gr._normalize_genome("hg38_ref"), "hg38")
        self.assertEqual(gr._normalize_genome("mm10_reference"), "mm10")
        self.assertEqual(gr._normalize_genome("hg38"), "hg38")
        # reference-only web filename: <guide>+<pam>_hg38+hg38_ref_6+3
        fn = gr._parse_results_filename(
            "GUIDENNN+NGG_hg38+hg38_ref_6+3_integrated_results.tsv"
        )
        self.assertEqual(fn.get("genome"), "hg38")
        self.assertNotIn("datasets", fn)  # 'ref' must NOT become a dataset
        # a real variant run still decodes datasets
        fn2 = gr._parse_results_filename(
            "GUIDENNN+NGG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv"
        )
        self.assertEqual(fn2.get("datasets"), "1000G+HGDP")

    def test_bundled_tsv_gz_round_trips(self):
        _, _, extract = self._build_and_extract()
        gz = self._dpath(extract, "integrated_results.tsv.gz")
        self.assertTrue(os.path.isfile(gz))
        with gzip.open(gz, "rt") as handle:
            content = handle.read()
        self.assertIn("Spacer+PAM", content.splitlines()[0])
        self.assertEqual(len(content.strip().splitlines()), len(_ROWS) + 1)

    def test_top1000_tsv_bundled_with_table_rows(self):
        _, _, extract = self._build_and_extract()
        top = self._dpath(extract, "top1000.tsv")
        self.assertTrue(os.path.isfile(top))
        with open(top) as handle:
            lines = handle.read().strip().splitlines()
        # header is now the CURATED set (rank first), NOT the raw Spacer+PAM dump
        header = lines[0].split("\t")
        self.assertEqual(header, gr.curated_headers(has_crista=True))
        self.assertEqual(header[0], "rank")
        self.assertNotIn("Spacer+PAM", lines[0])
        # header + the 5 off-targets (on-target mm+b==0 filtered out)
        self.assertEqual(len(lines), 5 + 1)
        # rank column is 1..5 in order
        ranks = [ln.split("\t")[0] for ln in lines[1:]]
        self.assertEqual(ranks, ["1", "2", "3", "4", "5"])

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
        self.assertIn("Aggregated Specificity Score (0-100; higher = more specific)", html)
        self.assertRegex(html, r"Nuclease</td><td>SpCas9<")
        # global MM/B matrix present, grouped REFERENCE vs VARIANT
        self.assertIn("Off-targets by Mismatch (MM) and Bulge (B)", html)
        self.assertIn('class="matrix"', html)
        self.assertIn(">REFERENCE", html)
        self.assertIn(">VARIANT", html)
        # per-mismatch columns up to the run's mm (6MM), on-target excluded
        self.assertIn(">6MM<", html)

    def test_section1_matrix_reconciles_to_grand_total(self):
        """Matrix INCLUDES perfect matches (mm+b==0): REFERENCE + VARIANT totals
        == grand total (every site lands in a cell), the origin split covers all
        rows, and the perfect match(es) sit in the 0 MM / 0 B cell. Bulge rows
        span 0..(bDNA+bRNA)=0..4, mm cols 0..6."""
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

        # matrix now INCLUDES perfect matches (mm+b==0): origin split over ALL rows
        variant, reference, ontarget = gr.partition_masks(df, cols)
        if "not_in_ref" in cols:
            variant_all = df[cols["not_in_ref"]].astype(str).str.strip().str.lower().eq("y")
        else:
            variant_all = df[cols["origin"]].astype(str).str.strip().str.lower().eq("alt")
        reference_all = ~variant_all
        by_label = {lbl: sum(r[1] for r in rows) for lbl, rows in matrix["groups"]}
        self.assertEqual(by_label["REFERENCE"], int(reference_all.sum()))
        self.assertEqual(by_label["VARIANT"], int(variant_all.sum()))

        # every row (incl perfect matches) lands in a cell; REF+VAR == grand total
        cell_sum = sum(
            sum(per_mm) for _lbl, rows in matrix["groups"] for _b, _tot, per_mm in rows
        )
        self.assertEqual(cell_sum, len(df))
        self.assertEqual(by_label["REFERENCE"] + by_label["VARIANT"], len(df))

        # the perfect match(es) appear in the 0 MM / 0 B cell
        pm_cell = sum(
            per_mm[0]
            for _lbl, rows in matrix["groups"]
            for b, _tot, per_mm in rows
            if b == 0
        )
        self.assertEqual(pm_cell, int(ontarget.sum()))

    def test_section2_four_scatter_panels_when_crista_present(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # strip the branding logo <img> so we count only PLOT images
        plot_html = re.sub(r'<img class="logo"[^>]*>', "", html)
        imgs = re.findall(
            r'src="data:image/svg\+xml;base64,([A-Za-z0-9+/=]+)"', plot_html
        )
        # 4 scatter panels + 1 population plot => 5 inline SVGs; each decodable
        self.assertEqual(len(imgs), 5)
        for encoded in imgs:
            raw = base64.b64decode(encoded)
            self.assertTrue(raw[:5] == b"<?xml" or b"<svg" in raw[:256])
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

    def test_section4_hybrid_panel_floors_note_and_maf_footnote(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn("Recommended validation panel", html)
        # hybrid panel wording + cap (~100) present (title no longer says a flat
        # "top 100" since the panel can exceed the cap when many are hard-included)
        self.assertIn("hybrid", html.lower())
        self.assertIn("worst-case", html.lower())
        self.assertIn("~100 sites", html)
        # C) EXPLICIT IN-REPORT METHODS NOTE with the REAL constants:
        #    hard-include (mm+b <= 2 OR CFD >= 0.5), fill to 100 by worst-case
        self.assertIn("<strong>Methods.</strong>", html)
        self.assertIn("hard-include", html.lower())
        self.assertIn("mismatches+bulges &le; 2 OR CFD &ge; 0.5", html)
        self.assertIn("fill", html.lower())
        self.assertIn("worst by ANY single", html)
        self.assertIn("NO category quota", html)
        # the FULL threshold table is kept (CFD>= {0.5,0.2,0.05},
        # mm+b<= {1,2,3,4})
        for t in ("0.5", "0.2", "0.05"):
            self.assertIn(f"CFD &ge; {t}", html)
        self.assertNotIn("CFD &ge; 0.1", html)  # 0.1 tier removed
        for t in ("1", "2", "3", "4"):
            self.assertIn(f"mismatches + bulges &le; {t}", html)
        # of the selected panel, 3 are variant-created (the fixture has exactly 3)
        self.assertRegex(html, r"<strong>3</strong>[^<]*are\s+variant-created")
        # E) MAF FOOTNOTE (the v2.4 blank-MAF explanation)
        self.assertIn(
            "MAF blank (&mdash;) = reference off-target (no variant), "
            "an indel-derived variant", html,
        )
        self.assertIn("frequency registry is SNP-only", html)
        self.assertIn("AC/AN over the genotyped panel", html)

    def test_section4_hybrid_floors_are_module_level(self):
        # HYBRID panel floors + cap + metrics are module-level constants
        self.assertEqual(gr.PANEL_CAP, 100)
        self.assertEqual(gr.PANEL_FLOOR_MMB, 2)
        self.assertEqual(gr.PANEL_FLOOR_CFD, 0.5)
        self.assertEqual(
            gr.PANEL_WORSTCASE_METRICS,
            (("cfd", "desc"), ("crista", "desc"), ("mmb", "asc")),
        )

    def test_section4_hybrid_selection_and_cap(self):
        """Hard-include (mm+b<=2 OR CFD>=0.5), fill by worst-case severity, cap."""
        import pandas as pd

        header = list(_HEADER)

        def _row(chrom, mm, b, cfd, crista, notref):
            r = ["G", chrom, "1", "+", "AAA", "AAA", "GGG", str(mm), str(b),
                 str(mm + b), notref and "alt" or "ref", "NA",
                 f"{cfd}", f"{cfd}", f"{cfd}", "NA", "NA", "NA", "NA",
                 f"{crista}", f"{crista}", f"{crista}", "NA", "NA", "NA",
                 "y" if notref else "NA", "NA", "GENE", "1.0", "NA", "NA"]
            return r

        rows = [
            _row("chrON", 0, 0, 1.0, 1.0, False),    # perfect match -> forced to top
            _row("chrCFD", 3, 2, 0.99, 0.10, True),  # HARD (CFD>=0.5)
            _row("chrMMB", 2, 0, 0.10, 0.10, False), # HARD (mm+b<=2)
            _row("chrCRI", 4, 3, 0.10, 0.99, True),  # not hard; worst by CRISTA
            _row("chrLOW", 6, 4, 0.05, 0.05, False), # low by every metric
        ]
        tsv = os.path.join(self.tmp, "wc.tsv")
        with open(tsv, "w") as h:
            h.write("\t".join(header) + "\n")
            for r in rows:
                h.write("\t".join(r) + "\n")
        df = pd.read_csv(tsv, sep="\t", dtype=str, na_filter=False)
        cols = gr._resolve(df.columns, list(gr._COLS.keys()))

        panel = gr.select_worstcase_panel(df, cols, cap=gr.PANEL_CAP)
        chroms = list(panel[cols["chrom"]])
        # perfect match (mm+b==0) is FORCED to the top of the panel (candidate cut site)
        self.assertEqual(chroms[0], "chrON")
        # hard-includes always present; CRISTA-only site enters via the fill ranks
        self.assertIn("chrCFD", chroms)
        self.assertIn("chrMMB", chroms)
        self.assertIn("chrCRI", chroms)
        # after the perfect match, hard-includes come first (before the fill); all-low last
        self.assertEqual({"chrCFD", "chrMMB"}, set(chroms[1:3]))
        self.assertEqual(chroms[-1], "chrLOW")

        # HARD-INCLUDES EXCEED THE CAP -> they are ALL kept (panel > cap allowed).
        # Build many mm+b<=2 hard-include sites, cap at 3.
        hard = [_row("chrON", 0, 0, 1.0, 1.0, False)]
        for i in range(10):
            hard.append(_row(f"h{i}", 1, 1, 0.10, 0.10, False))  # mm+b=2 -> hard
        tsvh = os.path.join(self.tmp, "wc_hard.tsv")
        with open(tsvh, "w") as h:
            h.write("\t".join(header) + "\n")
            for r in hard:
                h.write("\t".join(r) + "\n")
        dfh = pd.read_csv(tsvh, sep="\t", dtype=str, na_filter=False)
        colsh = gr._resolve(dfh.columns, list(gr._COLS.keys()))
        panelh = gr.select_worstcase_panel(dfh, colsh, cap=3)
        # 10 hard-includes > cap(3) -> all 10 kept; perfect match prepended -> 11, at top
        self.assertEqual(len(panelh), 11)
        self.assertEqual(list(panelh[colsh["chrom"]])[0], "chrON")

        # cap is honored when hard-includes are few: many mid sites -> exactly cap
        big = [_row("chrON", 0, 0, 1.0, 1.0, False)]
        for i in range(500):
            # mm+b=4 (>floor) and CFD=0.20 (<floor) -> NOT hard-included
            big.append(_row(f"c{i}", 3, 1, 0.20, 0.20, i % 2 == 0))
        tsv2 = os.path.join(self.tmp, "wc_big.tsv")
        with open(tsv2, "w") as h:
            h.write("\t".join(header) + "\n")
            for r in big:
                h.write("\t".join(r) + "\n")
        df2 = pd.read_csv(tsv2, sep="\t", dtype=str, na_filter=False)
        cols2 = gr._resolve(df2.columns, list(gr._COLS.keys()))
        panel2 = gr.select_worstcase_panel(df2, cols2, cap=gr.PANEL_CAP)
        # PANEL_CAP off-targets + the prepended perfect match, which is at the top
        self.assertEqual(len(panel2), gr.PANEL_CAP + 1)
        self.assertEqual(list(panel2[cols2["chrom"]])[0], "chrON")

    def test_perfect_match_flag_and_banner(self):
        """0-mm/0-bulge sites are flagged and drive the warning banner."""
        mmb_col = {"mmb": "Mismatches+bulges_(highest_CFD)"}
        # curated cell: "Yes" for a perfect match, blank otherwise
        self.assertEqual(
            gr._curated_cell(
                "perfect_match", {"Mismatches+bulges_(highest_CFD)": "0"}, mmb_col
            ),
            "Yes",
        )
        self.assertEqual(
            gr._curated_cell(
                "perfect_match", {"Mismatches+bulges_(highest_CFD)": "2"}, mmb_col
            ),
            gr.CURATED_MISSING,
        )
        # >= 2 perfect matches -> red "no unambiguous on-target" banner, sites listed
        vp2 = {
            "n_perfect": 2,
            "perfect_sites": [
                {"chrom": "chr7", "pos": "66994206", "strand": "-"},
                {"chrom": "chr7", "pos": "72830803", "strand": "+"},
            ],
        }
        b2 = gr.render_perfect_match_banner(vp2)
        self.assertIn("Multiple perfect matches", b2)
        self.assertIn("chr7:66994206", b2)
        self.assertIn("chr7:72830803", b2)
        self.assertIn("#dc2626", b2)  # red border
        # exactly 1 -> amber presumed-on-target note, not the red warning
        b1 = gr.render_perfect_match_banner(
            {"n_perfect": 1, "perfect_sites": [{"chrom": "chr5", "pos": "1", "strand": "+"}]}
        )
        self.assertIn("presumed", b1.lower())
        self.assertNotIn("Multiple perfect matches", b1)
        # 0 -> no banner
        self.assertEqual(
            gr.render_perfect_match_banner({"n_perfect": 0, "perfect_sites": []}), ""
        )
        # MULTI-GUIDE: two guides, ONE perfect match each -> NOT ambiguous (amber),
        # even though n_perfect == 2. The red banner is per-guide (max 1 here).
        bmg = gr.render_perfect_match_banner({
            "n_perfect": 2, "max_perfect_per_guide": 1, "n_guides_with_perfect": 2,
            "perfect_sites": [
                {"guide": "GUIDE_A", "chrom": "chr5", "pos": "1", "strand": "+"},
                {"guide": "GUIDE_B", "chrom": "chr9", "pos": "2", "strand": "-"},
            ],
        })
        self.assertNotIn("Multiple perfect matches", bmg)
        self.assertIn("one per guide", bmg)
        self.assertIn("GUIDE_A", bmg)  # guide labelled when several contribute
        # ONE guide with 2 perfect matches -> red ambiguity even in a 2-guide vp
        bamb = gr.render_perfect_match_banner({
            "n_perfect": 2, "max_perfect_per_guide": 2, "n_guides_with_perfect": 1,
            "perfect_sites": [
                {"guide": "GUIDE_A", "chrom": "chr5", "pos": "1", "strand": "+"},
                {"guide": "GUIDE_A", "chrom": "chr9", "pos": "2", "strand": "-"},
            ],
        })
        self.assertIn("Multiple perfect matches", bamb)
        self.assertIn("#dc2626", bamb)

    def test_annotation_columns_and_legend_drop_when_absent(self):
        """A run WITHOUT a given annotation must NOT show its all-'-' column NOR a
        legend entry documenting a screen that never ran; a no-annotation run omits
        Section 7 entirely."""
        import tempfile
        # legend gates on the resolved column keys present in `cols`
        self.assertEqual(gr.build_annotation_legend_html({}), "")  # nothing present
        only_cosmic = gr.build_annotation_legend_html({"cosmic": "Annotation_COSMIC"})
        self.assertIn("COSMIC", only_cosmic)
        self.assertNotIn("DNase", only_cosmic)  # DHS entry absent
        # a run whose TSV carries NO annotation columns -> no COSMIC/ENCODE column,
        # and no "7. Annotation legend" section
        header = _HEADER[:26]  # drop GENCODE/gene/dist/ENCODE/DHS/COSMIC
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(
                d, f"{_GUIDE}+NRG_hg38+hg38_1000G_HGDP_6+3_integrated_results.tsv"
            )
            with open(tsv, "w") as h:
                h.write("\t".join(header) + "\n")
                for row in _ROWS:
                    h.write("\t".join(row[:26]) + "\n")
            out = gr.build_report(
                integrated_tsv=tsv, out_zip=os.path.join(d, "r.zip")
            )
            with zipfile.ZipFile(out) as zf:
                html = zf.read("report.html").decode()
                top = zf.read("data/top1000.tsv").decode().splitlines()[0].split("\t")
        self.assertNotIn("7. Annotation legend", html)
        self.assertNotIn("COSMIC_cancer_gene", top)
        self.assertNotIn("GENCODE", top)

    def test_section6_table_shows_curated_columns_incl_annotations(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertIn('class="ottable"', html)
        # the table now renders the CURATED headers (shared with downloads),
        # including PAM_creation + CRISTA (computed) AND the ANNOTATION columns
        for hdr in gr.curated_headers(has_crista=True):
            self.assertIn(f"<th>{hdr}</th>", html)
        for ann in ("Gene", "Gene_distance_kb", "GENCODE", "ENCODE", "DHS"):
            self.assertIn(f"<th>{ann}</th>", html)
        self.assertIn("<th>PAM_creation</th>", html)
        self.assertIn("<th>CRISTA</th>", html)
        # the annotation VALUES from the fixture render in the table body
        tbody = html.split('class="ottable"')[-1].split("<tbody>")[-1].split("</tbody>")[0]
        for val in ("promoter", "enhancer", "DHS_3", "DHS_5", "exon", "lincRNA"):
            self.assertIn(val, tbody)
        # rows sorted by CFD desc, on-target excluded => 5 rows
        n_rows = tbody.count("<tr>")
        self.assertEqual(n_rows, 5)
        cfds = re.findall(r"<td>(0\.\d{4})</td>", tbody)
        self.assertEqual(cfds, sorted(cfds, reverse=True))
        # the pam_created value from the fixture appears
        self.assertIn("pam_created", html)

    def test_table_and_downloads_share_curated_columns(self):
        """The in-report table AND every download expose the SAME curated set."""
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        curated = gr.curated_headers(has_crista=True)
        # table header order == curated order
        thead = html.split('class="ottable"')[-1].split("<thead>")[-1].split("</thead>")[0]
        table_headers = re.findall(r"<th>([^<]+)</th>", thead)
        self.assertEqual(table_headers, curated)
        # top1000.tsv + panel_top100.tsv headers == curated order (spot check)
        for fname in ("top1000.tsv", "panel_top100.tsv"):
            with open(self._dpath(extract, fname)) as handle:
                header = handle.readline().rstrip("\n").split("\t")
            self.assertEqual(header, curated, fname)

    def test_section7_footer_disclaimer_and_provenance(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # disclaimer: research-only, AS IS, outputs-may-change, no liability
        self.assertIn("research purposes only", html)
        self.assertIn("AS IS", html)
        self.assertIn("are NOT a substitute for experimental validation", html)
        self.assertIn("MAY CHANGE as CRISPRme", html)
        self.assertIn("accept no liability", html)
        # feedback / issue tracker at the very end
        self.assertIn("open an issue", html)
        self.assertIn(
            'href="https://github.com/pinellolab/crisprme-plus/issues"', html
        )
        # provenance stamp: report generator version + source TSV basename
        self.assertIn("report generator v", html)
        self.assertIn(os.path.basename(self.tsv), html)
        # the CRISPRme+ version must RESOLVE (not "n/a"): source/Docker installs
        # have no package metadata, so it is read from crisprme.py's canonical
        # `version`. Regression guard for the footer showing "version: n/a".
        self.assertNotIn("CRISPRme+ version: n/a", html)
        self.assertRegex(html, r"CRISPRme\+ version: \d+\.\d+")

    def test_package_version_resolves_without_metadata(self):
        """_package_version falls back to crisprme.py's canonical version when the
        package isn't pip/conda-installed (the source/Docker `COPY .` case), so the
        footer never reads 'n/a'."""
        v = gr._package_version()
        self.assertIsNotNone(v)
        self.assertRegex(v, r"^\d+\.\d+")

    def test_variants_included_panel_and_variant_count(self):
        """'Variants included': no AF-threshold claim, no impossible 1e-05 min-MAF
        line; panel size (registry sample_count preferred) + SNP variant count."""
        # cheap per-dataset individual counter (fixture: 1000G n=3, HGDP n=2)
        self.assertEqual(
            gr.dataset_individual_counts(self.sid_dir), {"1000G": 3, "HGDP": 2}
        )
        self.assertEqual(gr.panel_and_variants_note({}), "")
        self.assertEqual(gr.panel_and_variants_note(None), "")
        # registry sample_count OVERRIDES the samplesID-derived count (the latter
        # can over-list vs the AC/AN denominator); SNP count is surfaced
        # with an indel count (build manifest): both counts shown, "all searched"
        vc = {"n_records": 106_664_924, "n_indels": 14_255_298,
              "databases": {"1000G": 2548, "HGDP": 929}}
        note = gr.panel_and_variants_note({"1000G": 3, "HGDP": 2}, vc)
        self.assertIn("<strong>3,477</strong> individuals", note)   # 2548 + 929
        self.assertIn("1000G n=2,548", note)
        self.assertIn("106,664,924", note)
        self.assertIn("14,255,298", note)
        self.assertIn("indels, all searched", note)
        self.assertNotIn("excluded", note.lower())
        self.assertNotIn("indel pipeline", note.lower())
        # no indel count (legacy install / idx fallback): SNPs shown, indels stated
        # as also searched, still no "(indel pipeline)" nor "excluded"
        vc2 = {"n_records": 106_664_924, "n_indels": None,
               "databases": {"1000G": 2548, "HGDP": 929}}
        note_sc = gr.panel_and_variants_note({}, vc2)
        self.assertIn("106,664,924</strong> SNPs", note_sc)
        self.assertIn("also searched", note_sc)
        self.assertNotIn("indel pipeline", note_sc.lower())
        # no registry -> fall back to samplesID counts, no variant count clause
        note2 = gr.panel_and_variants_note({"1000G": 3, "HGDP": 2})
        self.assertIn("<strong>5</strong> individuals", note2)
        self.assertNotIn("SNPs", note2)

    def test_panel_note_is_database_agnostic(self):
        """The panel/variant-count note works for ANY database(s) + edge counts --
        nothing is tuned to 1000G/HGDP."""
        # arbitrary future databases, three of them, with an indel count
        vc = {"n_records": 500, "n_indels": 42,
              "databases": {"gnomAD": 76156, "TOPMed": 132345, "MyCohort": 88}}
        note = gr.panel_and_variants_note({}, vc)
        self.assertIn("<strong>208,589</strong> individuals", note)  # 76156+132345+88
        self.assertIn("MyCohort n=88", note)
        self.assertIn("TOPMed n=132,345", note)
        self.assertIn("500</strong> SNPs and <strong>42</strong> indels", note)
        # a SNP-only database: n_indels == 0 is KNOWN (not "unknown") -> "0 indels"
        note0 = gr.panel_and_variants_note({}, {"n_records": 9, "n_indels": 0,
                                               "databases": {"SnpOnlyDB": 10}})
        self.assertIn("9</strong> SNPs and <strong>0</strong> indels, all searched",
                      note0)
        # unknown indel count (n_indels None) -> "also searched", no invented number
        noteU = gr.panel_and_variants_note({}, {"n_records": 9, "n_indels": None,
                                               "databases": {"SomeDB": 10}})
        self.assertIn("also searched", noteU)
        self.assertNotIn("indels, all searched", noteU)
        # rendered report (fixture has no registry -> samplesID fallback)
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        self.assertNotIn("applies no allele-frequency threshold", html)
        self.assertNotIn("no variants are excluded by allele frequency", html)
        # inclusion statement explicitly covers indels (they ARE searched)
        self.assertIn("SNPs and insertions/deletions", html)
        self.assertIn("<strong>5</strong> individuals", html)  # 3 + 2 (samplesID)
        # the impossible 1e-05 "finest resolution" line is GONE
        self.assertNotIn("Finest allele-frequency resolution", html)
        # the MAF footnote now explains the 1e-05 DISPLAY floor
        self.assertIn("display floor of 1&times;10<sup>&minus;5</sup>", html)
        # captions no longer use the "paper-style" jargon; CRISTA panel references CFD
        self.assertNotIn("paper-style", html)
        self.assertIn("the same ref/alt scatter as for the cfd score", html.lower())

    def test_registry_variant_count_idx_and_sidecar(self):
        """Cheap variant count + genotyped sample counts from the Tier-0 registry
        (sums reg_*.idx n_records; sidecar takes precedence; never scans VCFs)."""
        import tempfile, json as _json
        with tempfile.TemporaryDirectory() as d:
            rdir = os.path.join(d, "Results", "job1"); os.makedirs(rdir)
            reg = os.path.join(d, "Dictionaries", "registry_hg38_1000G_HGDP")
            os.makedirs(reg)
            dbs = {"1000G": {"sample_count": 2548}, "HGDP": {"sample_count": 929}}
            for i, n in enumerate((100, 50)):
                with open(os.path.join(reg, f"reg_chr{i+1}.idx"), "w") as h:
                    _json.dump({"n_records": n, "databases": dbs}, h)
            meta = {"datasets": "1000G+HGDP"}
            vc = gr._registry_variant_count(rdir, meta)
            self.assertEqual(vc["n_records"], 150)  # summed across chroms
            self.assertEqual(vc["databases"], {"1000G": 2548, "HGDP": 929})
            self.assertIsNone(vc["n_indels"])  # idx headers carry no indel count
            # a build-time sidecar takes precedence + supplies the indel count
            with open(os.path.join(reg, "variant_count.json"), "w") as h:
                _json.dump(
                    {"n_records": 106664924, "n_indels": 14255298, "databases": dbs}, h
                )
            vc2 = gr._registry_variant_count(rdir, meta)
            self.assertEqual(vc2["n_records"], 106664924)
            self.assertEqual(vc2["n_indels"], 14255298)
            # no registry resolvable -> None (report simply omits the count)
            self.assertIsNone(gr._registry_variant_count(None, meta))
            self.assertIsNone(
                gr._registry_variant_count(os.path.join(d, "no", "such"), meta)
            )

    def test_review_fixes_bulge_bridge_emptytier_legend(self):
        """Careful-review fixes: bulge label shows total (not misleading '(max 2)');
        off-target line bridges to the matrix via the perfect-match subtraction;
        empty tiers read as '0 (none)'; Gene/Gene_distance_kb are in the legend."""
        # bulge label: per-type + total (direct render, needs int bdna/brna)
        meta = {"guide_display": "G", "nuclease": "SpCas9", "pam": "NRG",
                "genome": "hg38", "datasets": "1000G", "mm": "5",
                "bdna": 2, "brna": 2, "bmax": "2", "max_edits": "5"}
        hs = gr.render_summary_and_matrix(meta, "6.1", None)
        self.assertIn("up to 4 total", hs)
        self.assertNotIn("(max 2)", hs)
        # in the rendered report (fixture: 1 perfect match, empty mm+b<=1 tier)
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # off-target subtraction bridge (6 total - 1 perfect = 5)
        self.assertIn("total sites &minus; 1 perfect match", html)
        # empty edit-distance tier reads as a count, not a blank
        self.assertIn("0 <span class='caption'>(none)</span>", html)
        # Gene + Gene_distance_kb documented in the annotation legend
        self.assertIn("Signed distance in kilobases", html)
        self.assertIn("<div class=\"legend-term\">Gene</div>", html)
        # max-total-edits caveat present (cap bounds the search, not per-site edits)
        self.assertIn("individual reported alignments", html)

    def test_run_manifest_bundled_and_valid(self):
        """The ZIP carries a machine-readable run_manifest.json (IND traceability)."""
        import json as _json
        _, names, extract = self._build_and_extract()
        self.assertIn("data/run_manifest.json", names)
        m = _json.load(open(self._dpath(extract, "run_manifest.json")))
        self.assertEqual(m["crisprme_version"] or "2.4.0", m["crisprme_version"] or "2.4.0")
        for k in ("guides", "pam", "genome", "variant_datasets", "search", "counts",
                  "source_integrated_results", "crispritz_note"):
            self.assertIn(k, m)
        self.assertIn("max_total_edits", m["search"])
        self.assertIn("off_targets", m["counts"])
        self.assertIn("v2.8.2", m["crispritz_note"])  # engine version recorded

    def test_pipeline_no_readonly_or_bare_relative_writes(self):
        """Regression guard for the read-only-container bug: the search pipeline must
        not write a bare-relative dummy.txt nor append to $6 without a writability
        guard (both crashed a reference-only run under Apptainer)."""
        p = os.path.join(os.path.dirname(__file__),
                         "submit_job_automated_new_multiple_vcfs.sh")
        src = open(p).read()
        # the fixed bare-relative dummy.txt write must be gone
        self.assertNotIn("echo -e \"dummy_file\" >dummy.txt", src)
        self.assertNotIn(">dummy.txt", src.replace('"${output_folder}/dummy.txt"', ""))
        # appending to the samplefile $6 must be writability-guarded
        self.assertIn('[ -w "$6" ]', src)

    def test_persona_audit_report_additions(self):
        """Scores/columns legend (CFD/CRISTA + citations + alignment notation),
        'What to do next' box, per-guide max-edits caveat, and the perfect-match
        banner's REF/ALT origin tag."""
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # scores & columns legend with citations + lay gloss + notation
        self.assertIn("scores, columns", html.lower())          # section 7 retitle
        self.assertIn("Cutting Frequency Determination", html)   # CFD defined
        self.assertIn("Doench", html)                            # CFD citation
        self.assertIn("Abadi", html)                             # CRISTA citation
        self.assertIn("extrapolation beyond the model", html)    # bulge caveat
        self.assertIn("lowercase</b> = a mismatch", html)        # alignment notation
        # "What to do next" box
        self.assertIn("What to do next", html)
        self.assertIn("rhAMP-Seq", html)
        # max-total-edits caveat surfaces the observed max (fixture: cap 6, max mm+b 4
        # -> no caveat; force a case where obs>cap via render_summary_and_matrix)
        meta = {"guide_display": "G", "nuclease": "SpCas9", "pam": "NRG",
                "genome": "hg38", "datasets": "1000G", "mm": "5", "bdna": 2,
                "brna": 2, "bmax": "2", "max_edits": "5", "obs_max_mmb": 9}
        hs = gr.render_summary_and_matrix(meta, "6.1", None)
        self.assertIn("search cap", hs)
        self.assertIn("up to 9 total edits", hs)
        # perfect-match banner tags REF vs variant-created origin
        b = gr.render_perfect_match_banner({
            "n_perfect": 1, "max_perfect_per_guide": 1, "n_guides_with_perfect": 1,
            "perfect_sites": [{"chrom": "chr2", "pos": "100", "strand": "+",
                               "origin": "ref", "maf": None}],
        })
        self.assertIn("reference", b.lower())
        b2 = gr.render_perfect_match_banner({
            "n_perfect": 2, "max_perfect_per_guide": 2, "n_guides_with_perfect": 1,
            "perfect_sites": [
                {"chrom": "chr7", "pos": "1", "strand": "+", "origin": "ref", "maf": None},
                {"chrom": "chr7", "pos": "2", "strand": "-", "origin": "alt", "maf": 0.0035}],
        })
        self.assertIn("variant-created", b2.lower())
        self.assertIn("only in carriers", b2.lower())

    def test_inline_filename_mentions_are_clickable(self):
        """Any inline <code>FILE</code> reference to a bundled download (in prose,
        not just the download buttons) is itself a clickable link into data/."""
        h = "raw <code>integrated_results.tsv.gz</code>; also <code>other.tsv</code>."
        out = gr._linkify_bundled_filenames(h, {"integrated_results.tsv.gz"})
        self.assertIn(
            '<a href="data/integrated_results.tsv.gz" download>'
            "<code>integrated_results.tsv.gz</code></a>", out
        )
        self.assertIn("<code>other.tsv</code>", out)         # not bundled -> plain
        self.assertNotIn("<code>other.tsv</code></a>", out)  # not linked
        # in a real report, the prose mentions are links (single, not double-wrapped)
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        for name in ("integrated_results.tsv.gz", "panel_top100.tsv",
                     "variant_created.tsv"):
            self.assertIn(
                f'<a href="data/{name}" download><code>{name}</code></a>', html
            )
            self.assertNotIn(f"<code>{name}</code></a></a>", html)  # no double link

    def test_self_contained_offline_and_relative_links(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # Fully self-contained: renders offline with NO external RESOURCE loading
        # (no phone-home). Outbound <a href> links (the license contact mailto, the
        # issue tracker) are allowed -- they load nothing until the reader clicks.
        self.assertNotIn("<script", html.lower())
        self.assertNotIn("<link", html.lower())
        self.assertNotIn('src="http', html)  # no external images / iframes
        self.assertNotIn("url(http", html)    # no external CSS resources
        # every http(s) URL must sit inside an <a href="..."> (a plain link), never
        # a resource reference
        i = html.find("http")
        while i != -1:
            self.assertEqual(html[i - 6 : i], 'href="', f"non-link URL: {html[i-10:i+40]!r}")
            i = html.find("http", i + 1)
        # RELATIVE download links into the data/ subfolder (report.html is the
        # only top-level file, so it is self-evident what to open)
        self.assertIn('href="data/integrated_results.tsv.gz"', html)
        self.assertIn('href="data/top1000.tsv"', html)

    def test_multi_snp_haplotype_min_maf_and_first_rsid(self):
        _, _, extract = self._build_and_extract()
        html = self._read(os.path.join(extract, "report.html"))
        # the multi-SNP row (rs333,rs444 / NA,0.005,0.003) must show rs333 and
        # the MIN maf 0.003 -> "3.00e-03"
        self.assertIn("rs333", html)
        self.assertIn("3.00e-03", html)


if __name__ == "__main__":
    unittest.main(verbosity=2)

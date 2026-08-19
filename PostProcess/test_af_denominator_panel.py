"""Unit test for the AF panel denominator over VCF-FILTERED per-db samplesID lists
(#117/#121).

The registry / population-summary Panel builds AN (allele number) as sum-of-ploidy
over EXACTLY the samplesID rows handed to ``tier0_compile.build_sample_meta`` (via the
Tier-1 SampleAxis). If a per-db samplesID OVER-LISTS its VCF (1000G metadata lists
~3500 but the phased VCF has 2548), the panel counts phantom hom-ref individuals and
AN is inflated. The merge-script fix writes VCF-FILTERED per-db lists so AN reflects
the genotyped panel only.

This test drives that math directly: build a sample_meta from two synthetic per-db
samplesID files, build the SampleAxis + population_summary.Panel, and assert the GLOBAL
AN equals 2 * (genotyped-sample count) for autosomes -- and that over-listing would
inflate it. STDLIB ONLY (tier0_compile / tier1_genotypes / population_summary are all
stdlib-only via tier0_registry).

Run with:
    cd PostProcess && python3 -m unittest test_af_denominator_panel -v
"""

import os
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

import tier0_compile as t0c
import tier1_genotypes as t1g
import population_summary as ps


def _write_samplesid(path, sample_ids, subpop="EUR", sex="M"):
    with open(path, "w") as fh:
        fh.write("#SAMPLE_ID\tpopulation\tsuperpopulation\tsex\n")
        for sid in sample_ids:
            fh.write("%s\tPOP\t%s\t%s\n" % (sid, subpop, sex))


class AfPanelDenominator(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _panel_for(self, db_to_samplesid, chrom="chr1"):
        sample_meta, _overlaps = t0c.build_sample_meta(db_to_samplesid)
        axis = t1g.build_sample_axis(sample_meta)
        ploidy_of = t0c.ploidy_of_for_chrom(chrom)
        return ps.Panel(axis, ploidy_of), len(axis)

    def test_filtered_lists_give_genotyped_only_an(self):
        # Genotyped panel: 1000G has 3 genotyped, HGDP has 2 genotyped.
        f1 = os.path.join(self.tmp, "hg38_1000G.panel.samplesID.txt")
        f2 = os.path.join(self.tmp, "hg38_HGDP.panel.samplesID.txt")
        _write_samplesid(f1, ["A", "B", "C"], subpop="EUR")
        _write_samplesid(f2, ["D", "E"], subpop="AFR")
        panel, n = self._panel_for({"1000G": f1, "HGDP": f2})
        # 5 diploid autosomal samples -> GLOBAL AN = 2*5 = 10, n_called = 5.
        self.assertEqual(n, 5)
        self.assertEqual(panel.global_[0], 10, "AN must be 2 * genotyped samples")
        self.assertEqual(panel.global_[1], 5, "n_called must be genotyped samples")
        # per-db AN
        self.assertEqual(panel.db[ps.db_group_id("1000G")][0], 6)   # 2*3
        self.assertEqual(panel.db[ps.db_group_id("HGDP")][0], 4)    # 2*2

    def test_overlisting_inflates_an(self):
        # Demonstrate the bug the filter prevents: an over-listed samplesID (2 extra
        # phantom samples that are NOT in the VCF) inflates AN by 2 alleles each.
        f_over = os.path.join(self.tmp, "hg38_1000G.over.samplesID.txt")
        _write_samplesid(f_over, ["A", "B", "C", "PHANTOM1", "PHANTOM2"])
        panel_over, n_over = self._panel_for({"1000G": f_over})
        self.assertEqual(n_over, 5)
        self.assertEqual(panel_over.global_[0], 10)  # 2*5, inflated by the 2 phantoms

        f_filt = os.path.join(self.tmp, "hg38_1000G.filt.samplesID.txt")
        _write_samplesid(f_filt, ["A", "B", "C"])  # VCF-filtered
        panel_filt, n_filt = self._panel_for({"1000G": f_filt})
        self.assertEqual(n_filt, 3)
        self.assertEqual(panel_filt.global_[0], 6)   # 2*3, correct

        self.assertGreater(panel_over.global_[0], panel_filt.global_[0],
                           "over-listing must inflate AN vs the VCF-filtered panel")

    def test_combined_panel_an_is_sum_of_filtered_dbs(self):
        # The locked design's arithmetic shape: combined AN = 2 * (1000G_filtered +
        # HGDP_filtered). Use small stand-ins for 2548 + 929.
        f1 = os.path.join(self.tmp, "hg38_1000G.c.samplesID.txt")
        f2 = os.path.join(self.tmp, "hg38_HGDP.c.samplesID.txt")
        _write_samplesid(f1, ["S%d" % i for i in range(2548 % 50 + 4)], subpop="EUR")
        _write_samplesid(f2, ["H%d" % i for i in range(929 % 50 + 3)], subpop="AFR")
        n1 = 2548 % 50 + 4
        n2 = 929 % 50 + 3
        panel, n = self._panel_for({"1000G": f1, "HGDP": f2})
        self.assertEqual(n, n1 + n2)
        self.assertEqual(panel.global_[0], 2 * (n1 + n2))


if __name__ == "__main__":
    unittest.main()

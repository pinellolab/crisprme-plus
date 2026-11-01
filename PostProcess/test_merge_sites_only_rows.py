"""A sites-only dataset's rows stay at their sites in the merge step of a search.

CONTEXT
-------
A row of a sites-only dataset (``hg38_mega``, a registry with no genotypes) has no
carriers: its Samples is ``NA``, and ``adjust_cols.py`` reads that with pandas and writes
it back empty (its AF too, when it has none). The merge step's three sorts
(``submit_job_automated_new_multiple_vcfs.sh``, ``sort -k16,16 -k5,5 -k7,7n ...``) split
fields on blanks, so an empty column moves the later ones up: for such a row ``-k16,16``
read the rsID, not the guide, and all such rows sorted by rsID around the others.
``remove_contiguous_samples.py`` starts a new cluster only when the position rises by more
than the window, so rows in rsID order chained unrelated sites together: a site was listed
more than once, a variant row was marked ``Not_found_in_REF`` though its site has a
reference row, and a site could lose its row altogether. A row with only its AF empty (a
VCF with no ``AF``) stayed at its site but went by its REF score: the merge keeps the first
row of a variant, so a site could keep an alignment with a lower score.

With the fields tab-separated those rows reach their sites. There a row of the sites-only
dataset can meet the row another dataset gives for the same variant, and the merge joins
the two rows' lists: an empty column must add no empty item to them (the table would read
``,HG00096,HG00112``).

WHAT THIS TEST DOES
-------------------
Takes the three sort lines from the REAL pipeline script and runs them with bash; the rows
go through the REAL ``adjust_cols.py`` (so the empty columns come as they do in a search),
then the REAL ``remove_contiguous_samples.py`` runs with the indices
``merge_close_targets_cfd.sh`` passes.

ASSERTS (the ones the fix changes fail before it and pass after it):
  * the three sorts keep each row at its site, whatever its rsID;
  * the merge lists each site once, with the sites-only row among the site's rows and not
    marked as not found in the reference, and no site is lost;
  * the three lists of discarded rows, pasted by position as the pipeline does, give one
    site on each row;
  * a variant row two datasets give is one row, with the other dataset's carriers and no
    empty item in its lists, in either order;
  * of a variant's alignments with an empty AF, the merge keeps the one with the highest
    score, not the one with the highest REF score.
And these, on both sides of the fix: ``adjust_cols.py`` leaves a sites-only row's Samples
empty (the premise); a sites-only row with no other dataset's row keeps its empty lists;
rows with no empty column sort as they did, so a search with no sites-only dataset does
not change.

Run with:
    cd PostProcess && python3 -m unittest test_merge_sites_only_rows -v
"""

import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
PIPELINE = os.path.join(PP, "submit_job_automated_new_multiple_vcfs.sh")

GUIDE = "CTAACAGTTGCTTTTATCACNNN"
# as the post-analysis writes a row, before adjust_cols.py
PRE_ADJUST = ("#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\tDirection\t"
              "Mismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\tAnnotation_Type\t"
              "Real_Guide\trsID\tAF\tSNP\tReference\tCFD_ref\tCFD\t#Seq_in_cluster").split("\t")
LISTS = ("bestCFD", "bestmmblg", "bestCRISPR_BULGE")
TAB_OPTION = " -t $'\\t'"

# columns once adjust_cols.py has moved them (0-based)
POSITION, NOT_IN_REF, SAMPLES, RSID, AF, SNP, CFD = 6, 12, 13, 16, 17, 18, 20


def _row(position, dna="CTAACAGTTGCTTTTATCACTGG", samples="n", rsid="n", af="n", snp="n",
         cfd="0.5", total=2):
    values = {
        "#Bulge_type": "X", "crRNA": GUIDE, "DNA": dna, "Chromosome": "chr1",
        "Position": position, "Cluster_Position": position, "Direction": "+",
        "Mismatches": total, "Bulge_Size": 0, "Total": total, "PAM_gen": "n", "Var_uniq": "n",
        "Samples": samples, "Annotation_Type": "n", "Real_Guide": GUIDE, "rsID": rsid,
        "AF": af, "SNP": snp, "Reference": "n" if snp == "n" else "CTAACAGTTGCTTTTATCACTGG",
        "CFD_ref": cfd if snp == "n" else "0.01", "CFD": cfd, "#Seq_in_cluster": 0,
    }
    return "\t".join(str(values[c]) for c in PRE_ADJUST)


def _sites_only(position, rsid, snp, af="NA", cfd="0.9", total=1):
    """A variant row of a sites-only dataset as the post-analysis writes it."""
    return _row(position, dna="CTAACAGTTGCaTTTATCACTGG", samples="NA", rsid=rsid, af=af,
                snp=snp, cfd=cfd, total=total)


def _genotyped(position, samples, rsid, snp, af="0.25", cfd="0.9", total=1):
    return _row(position, dna="CTAACAGTTGCaTTTATCACTGG", samples=samples, rsid=rsid, af=af,
                snp=snp, cfd=cfd, total=total)


def _pipeline_sort_lines():
    """The merge step's own lines, as the script has them."""
    with open(PIPELINE) as fh:
        lines = fh.read().splitlines()
    picked = [ln for ln in lines
              if re.match(r"(head -1|tail -n \+2) \$final_res\.best\w+\.txt\b.*(\$final_res\.tmp)", ln)]
    assert len(picked) == 6, "the merge step's sort lines moved: %r" % picked
    assert all(k in "\n".join(picked) for k in LISTS)
    return picked


class SitesOnlyRowsInTheMergeStep(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="merge_sites_only_")
        self.addCleanup(shutil.rmtree, self.d, True)
        self.final = os.path.join(self.d, "final")

    def _adjusted(self, name, rows):
        """Write rows as the post-analysis does and run the REAL adjust_cols.py."""
        path = os.path.join(self.d, name + ".adj")
        with open(path, "w") as fh:
            fh.write("\t".join(PRE_ADJUST) + "\n")
            fh.write("".join(r + "\n" for r in rows))
        subprocess.run([sys.executable, "adjust_cols.py", path], cwd=PP, check=True,
                       capture_output=True, text=True)
        with open(path) as fh:
            return fh.read().splitlines()

    def _sorted_lists(self, rows):
        """The three lists, each of these rows, through adjust_cols.py and the pipeline's
        own sort lines. Returns {list name: path}."""
        for name in LISTS:
            with open(self.final + "." + name + ".txt", "w") as fh:
                fh.write("\n".join(self._adjusted(name, rows)) + "\n")
        subprocess.run(["bash", "-c", "\n".join(_pipeline_sort_lines())], check=True,
                       capture_output=True, text=True,
                       env=dict(os.environ, final_res=self.final, output_folder=self.d))
        return {name: self.final + "." + name + ".txt" for name in LISTS}

    @staticmethod
    def _read(path):
        with open(path) as fh:
            return [ln.rstrip("\n").split("\t") for ln in fh][1:]

    def _merge(self, path, pivot="score"):
        """The REAL merge with merge_close_targets_cfd.sh's indices; (best, discarded)."""
        out = path + ".trimmed"
        r = subprocess.run(
            [sys.executable, "remove_contiguous_samples.py", path, out, "3", "5", "7", "11",
             "16", "19", "21", pivot, "mm+bulges"],
            cwd=PP, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stderr, "", "a pipeline stage must write nothing on stderr")
        return self._read(out), self._read(out + ".discarded_samples")

    # a site with a reference row and a sites-only variant row, an rsID that sorts after
    # the guide, and another site: the review's three-row case
    ROWS = [_row(100), _sites_only(100, "rs999", "chr1_111_T_A"), _row(5000, cfd="0.2", total=3)]

    def test_adjust_cols_leaves_the_samples_of_a_sites_only_row_empty(self):
        """The premise: a real search's rows have an empty column here."""
        adjusted = self._adjusted("premise", [_sites_only(100, "rs999", "chr1_111_T_A")])
        self.assertEqual(adjusted[1].split("\t")[SAMPLES], "")
        self.assertEqual(adjusted[1].split("\t")[AF], "")

    def test_the_sorts_keep_each_row_at_its_site(self):
        paths = self._sorted_lists(self.ROWS)
        for name in LISTS:
            self.assertEqual([r[POSITION] for r in self._read(paths[name])],
                             ["100", "100", "5000"],
                             "%s: a row with an empty Samples went by its rsID" % name)

    def test_an_rsid_before_the_guide_keeps_its_row_at_its_site_too(self):
        """'.' and '22:...' sort before the guide, 'rs...' after it: both moved."""
        paths = self._sorted_lists([_row(100), _row(5000, cfd="0.2"),
                                    _sites_only(5000, ".", "chr1_5001_T_A")])
        for name in LISTS:
            self.assertEqual([r[POSITION] for r in self._read(paths[name])],
                             ["100", "5000", "5000"])

    def test_the_merge_lists_each_site_once_and_loses_none(self):
        paths = self._sorted_lists(self.ROWS)
        best, disc = self._merge(paths["bestCFD"])
        self.assertEqual([r[POSITION] for r in best], ["100", "5000"],
                         "each site has one best row, none is lost")
        self.assertEqual(best[0][SNP], "chr1_111_T_A")
        self.assertEqual(best[0][NOT_IN_REF], "n",
                         "the site has a reference row: the variant row is found in it")
        self.assertEqual([(r[POSITION], r[SNP]) for r in disc], [("100", "n")])

    def test_rows_in_rsid_order_do_not_chain_two_sites(self):
        """Two sites each with a reference row and a sites-only row whose rsIDs sort
        against their positions: the merge must give two clusters, one row each, and the
        lists of discarded rows, pasted as the pipeline does, one site on each row."""
        rows = [_row(100), _sites_only(100, "rs9", "chr1_111_T_A"),
                _row(5000, cfd="0.2", total=3), _sites_only(5000, "rs1", "chr1_5011_T_A", cfd="0.8"),
                _sites_only(5000, "rs2", "chr1_5012_T_A", cfd="0.7")]
        paths = self._sorted_lists(rows)
        disc_paths = []
        for name, pivot in zip(LISTS, ("score", "total", "score")):
            best, disc = self._merge(paths[name], pivot)
            self.assertEqual([r[POSITION] for r in best], ["100", "5000"], name)
            disc_paths.append(paths[name] + ".trimmed.discarded_samples")
        pasted = subprocess.run(["pr", "-m", "-t", "-J"] + disc_paths, check=True,
                                capture_output=True, text=True).stdout.splitlines()
        self.assertEqual(len(pasted), 4)  # the header, 1 at site 100, 2 at site 5000
        for line in pasted[1:]:
            cols = line.split("\t")
            self.assertEqual(len(cols), 66)
            self.assertEqual({cols[POSITION], cols[22 + POSITION], cols[44 + POSITION]},
                             {cols[POSITION]}, "a row's three blocks must be at one site")

    def _merged_variant_row(self, rows):
        best, disc = self._merge(self._sorted_lists(rows)["bestCFD"])
        self.assertEqual(disc, [], "one variant: one row")
        self.assertEqual(len(best), 1)
        return best[0]

    def test_a_variant_two_datasets_give_has_no_empty_item_in_its_lists(self):
        """Once the rows reach their site, the sites-only row meets the genotyped one."""
        genotyped = _genotyped(100, "S2,S1", "rs5", "chr1_111_T_A", af="0.25")
        sites_only = _sites_only(100, "rs5", "chr1_111_T_A")
        for rows in ([genotyped, sites_only], [sites_only, genotyped]):
            row = self._merged_variant_row(rows)
            self.assertEqual(row[SAMPLES], "S1,S2", "the carriers, and no empty item")
            self.assertEqual(row[AF], "0.25", "the allele frequency, and no empty item")
            self.assertEqual(row[RSID], "rs5")

    def test_a_sites_only_allele_frequency_is_still_listed(self):
        """Each dataset's allele frequency stays on the row, as two datasets' do."""
        row = self._merged_variant_row(
            [_genotyped(100, "S1", "rs5", "chr1_111_T_A", af="0.25"),
             _sites_only(100, ".", "chr1_111_T_A", af="0.1")])
        self.assertEqual(row[SAMPLES], "S1")
        self.assertEqual(row[AF], "0.1,0.25")
        self.assertEqual(row[RSID], ".,rs5")

    def test_a_sites_only_row_alone_keeps_its_empty_lists(self):
        row = self._merged_variant_row([_sites_only(100, ".", "chr1_111_T_A")])
        self.assertEqual((row[SAMPLES], row[AF]), ("", ""))
        self.assertEqual(row[NOT_IN_REF], "y", "no reference row at this site")

    def test_an_empty_af_does_not_make_the_merge_keep_a_lower_score(self):
        """A dataset whose VCF gives no AF: the rows stay at their site, but for a row
        with an empty AF the sort's -k21,21rg read the column after the score, the REF
        score. The merge keeps a variant's first row, so it kept the alignment with the
        higher REF score. It must keep the one with the higher score."""
        low = _genotyped(100, "S1", ".", "chr1_111_T_A", af="NA", cfd="0.2")
        high = _genotyped(100, "S1", ".", "chr1_111_T_A", af="NA", cfd="0.9").replace(
            "CTAACAGTTGCaTTTATCACTGG", "CTAACAGTTGCcTTTATCACTGG")
        # the REF scores run the other way
        low, high = low.replace("\t0.01\t", "\t0.7\t"), high.replace("\t0.01\t", "\t0.1\t")
        row = self._merged_variant_row([low, high])
        self.assertEqual(row[CFD], "0.9", "the variant's best alignment is lost")

    def _same_with_and_without_tab(self, sort_cmd, body):
        """``sort_cmd`` (the script's, whichever it holds) with and without the tab."""
        without = sort_cmd.replace(TAB_OPTION, "")
        self.assertIn("LC_ALL=C sort ", without)
        outputs = [
            subprocess.run(["bash", "-c", cmd + ' "$1"', "sort_line", body], check=True,
                           capture_output=True, text=True).stdout
            for cmd in (without, without.replace("sort ", "sort" + TAB_OPTION + " ", 1))
        ]
        self.assertTrue(outputs[0])
        self.assertEqual(outputs[0], outputs[1],
                         "with no empty column the tab must change nothing: %s" % sort_cmd)

    def test_rows_with_no_empty_column_sort_as_before(self):
        """A search with no sites-only dataset has no empty column, and the tab changes
        nothing for it: sort's blank-separated fields held the same columns, each with
        its one leading tab. Random rows, through the script's lines with and without it."""
        rng = random.Random(281)
        rows = []
        for _ in range(400):
            row = _row(rng.randrange(100, 140), dna=rng.choice("ACGT") * 23,
                       samples=rng.choice(["n", "S1", "S1,S2", "NA18597"]),
                       rsid=rng.choice(["n", ".", "rs5", ".,rs7"]),
                       af=rng.choice(["n", "0.1", "0.1,0.2"]),
                       snp=rng.choice(["n", "chr1_111_T_A"]),
                       cfd=rng.choice(["0.0", "0.25", "0.5", "1.0", "4.33906e-23", "1e-05"]),
                       total=rng.randrange(0, 7))
            # the guide and the chromosome vary too, in the columns the first keys read
            row = row.replace(GUIDE, rng.choice([GUIDE, "AAAACAGTTGCTTTTATCACNNN"]))
            rows.append(row.replace("chr1", rng.choice(["chr1", "chr2"]), 1))
        body = os.path.join(self.d, "body")
        with open(body, "w") as fh:
            fh.write("\n".join(self._adjusted("many", rows)[1:]) + "\n")
        sort_cmds = [m.group(0) for m in (re.search(r"LC_ALL=C sort.*?(?= -T )", ln)
                                          for ln in _pipeline_sort_lines()) if m]
        self.assertEqual(len(sort_cmds), 3)
        for sort_cmd in sort_cmds:
            self._same_with_and_without_tab(sort_cmd, body)


if __name__ == "__main__":
    unittest.main()

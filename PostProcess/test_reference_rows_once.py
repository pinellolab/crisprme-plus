"""Each reference alignment of a cluster is listed, and counted, once.

CONTEXT
-------
``remove_contiguous_samples.py`` (run by ``merge_close_targets_cfd.sh`` for the CFD,
fewest mismatches+bulges and CRISPR-Bulge files) keeps a cluster's best alignment and
writes the others to ``.discarded_samples`` (the alternative alignments), and it gives
every row the count of the others (``#Seq_in_cluster``, reported as ``Other_motifs``).
Two steps upstream hand it the same reference row more than once:

  * on the dict-less tiers, ``new_simple_analysis._finalize_reference_entry`` emits the
    reference alignment of every variant window (LOCUS-COVERAGE), which the reference
    search also reports, as a row identical to that one;
  * each variant dataset's post-analysis carries every row of the one reference search,
    so a search of two datasets hands the merge each reference row twice.

The merge already makes one row of a variant row that two datasets give, but it
appended every reference row as it came, so the reference alignments near a SNP were
listed twice on the tiers (``Other_motifs`` 444 where the dictionaries gave 222), and a
search of two datasets listed every reference alignment once per dataset.

WHAT THIS TEST DOES
-------------------
Builds the rows exactly as the pipeline does: the tier path's rows of one variant
window come from the REAL ``iupac_decomposition`` (AST harness of
``test_phased_haplotype``, with a genotype-store sentinel so the dict-less observed
branch runs), the reference search's row as the main loop of ``new_simple_analysis``
appends a row with no IUPAC code (``split + ["n", 55, 0]``), both with the columns
``scriptAnalisiNNN_v3.sh`` adds (``n n n n <guide> n n n``). They go through the
write-loop transform and ``adjust_cols.py``'s reorder (helpers of
``test_reference_row_downstream``), then the REAL ``remove_contiguous_samples.py``
runs on them with the indices ``merge_close_targets_cfd.sh`` passes.

ASSERTS (all but the first fail before the fix and pass after it; the first is their
premise and passes on both):
  * the tier path gives the reference search's row again, identical;
  * one variant window on the tiers: the reference alignment is listed once and the
    best row counts one other alignment, exactly as on the dictionaries, for the score
    and the fewest mismatches+bulges merges;
  * two variant datasets: each reference alignment is listed once, and the variant row
    both datasets give is still merged into one with the carriers of both;
  * two reference rows that differ (two alignments) are both kept and counted;
  * two copies of a reference alignment that differ only in the CRISPR-Bulge score (each
    dataset's post-analysis scores the reference rows itself) are one alignment: the
    three criteria's lists, which the pipeline pastes by position, stay the same length
    and name the same sites in the same order.

Run with:
    cd PostProcess && python3 -m unittest test_reference_rows_once -v
"""

import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PP)

from test_phased_haplotype import _load_pure_functions, GUIDE, GUIDE_NO_PAM, TARGET_LEN
from test_reference_row_downstream import (
    _GtSentinel,
    _adjust_cols,
    _post_adjust_header,
    _to_disk_pre_adjust,
)

# POST-ADJUST columns (see test_reference_row_downstream): DNA 2, Reference 3,
# Position 5, Cluster_Position 6, Total 10, Samples 13, SNP 18, #Seq_in_cluster 19,
# CFD 20.
DNA_COL = 2
SAMPLES_COL = 13
SNP_COL = 18
SEQ_IN_CLUSTER_COL = 19
CFD_COL = 20

# The 8 columns scriptAnalisiNNN_v3.sh appends to every search row: PAM_gen, Var_uniq,
# Samples, Annotation_Type, Real_Guide (the crRNA without gaps), rsID, AF, SNP.
ADDED_COLUMNS = ["n", "n", "n", "n", GUIDE.replace("-", ""), "n", "n", "n"]


def _window(n_snps):
    """The reference window of one candidate: 'C' (a mismatch against the all-'A'
    guide) at each SNP position, 'A' elsewhere, NGG PAM; its search rows on the
    variant genome (IUPAC 'W' at each SNP) and on the reference genome (CRISPRitz
    writes each mismatch in lower case)."""
    ref = ["A"] * TARGET_LEN
    for i in range(n_snps):
        ref[i] = "C"
    ref[20], ref[21], ref[22] = "A", "G", "G"
    genome = "".join(ref)
    dna_var = "W" * n_snps + genome[n_snps:]
    dna_ref = "c" * n_snps + genome[n_snps:]
    head = ["X", GUIDE, None, "chrT", "0", "0", "+", str(n_snps), "0", str(n_snps)]
    var_split = head[:2] + [dna_var] + head[3:] + ADDED_COLUMNS
    ref_split = head[:2] + [dna_ref] + head[3:] + ADDED_COLUMNS
    return genome, var_split, ref_split


def _tier_rows(gt_by_col, sample_prefix=""):
    """The rows the dict-less tiers give for one variant window: the REAL
    iupac_decomposition on the observed branch (a variant row per carried haplotype,
    plus the reference row of _finalize_reference_entry)."""
    n_snps = len(gt_by_col)
    genome, var_split, _ = _window(n_snps)
    mydict = {}
    for i in range(n_snps):
        mydict["chrT," + str(i + 1)] = (
            sample_prefix + gt_by_col[i] + ";C,A;rs%d;0.10" % i)
    ns = _load_pure_functions(dict(
        genomeStr=genome, current_chr="chrT", mydict=mydict, myreg=None,
        mygt=_GtSentinel(), haplotype_check=True, IUPAC_CAP=10,
        hvdr_bed=io.StringIO(), pam="NGG", pos_beg=0, pos_end=-3,
        pam_begin=-3, pam_end=None, allowed_mms=6,
    ))
    ns["_phase_confirmation_rows"].clear()
    ns["_phase_confirmation_keys"].clear()
    cluster = []
    ns["iupac_decomposition"](var_split, GUIDE.replace("-", ""), GUIDE_NO_PAM, cluster)
    return cluster


def _reference_search_row(n_snps):
    """The reference search's row for the same window, as new_simple_analysis's main
    loop appends a row with no IUPAC code."""
    _, _, ref_split = _window(n_snps)
    return ref_split + ["n", 55, 0]


def _on_disk(rows, position=100):
    """The write loop + adjust_cols.py: the post-adjust rows remove_contiguous_samples
    reads. A variant row scores above the reference row, as the alt allele matches."""
    out = []
    for r in rows:
        cfd = "0.500" if r[19] == 55 else "0.900"
        out.append(_adjust_cols(_to_disk_pre_adjust(r, "chrT", position, cfd=cfd)))
    return out


class ReferenceRowsListedOnce(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="refrows_once_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def _merge(self, rows, pivot="score"):
        """Run the REAL remove_contiguous_samples.py with merge_close_targets_cfd.sh's
        indices (chrom=5 pos=7 total=11 true_guide=16 snp_info=19 cfd=21) and return
        (best rows, discarded rows)."""
        infile = os.path.join(self.d, "in.%s.txt" % pivot)
        with open(infile, "w") as fh:
            fh.write("\t".join(_post_adjust_header()) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        outfile = os.path.join(self.d, "out.%s.txt" % pivot)
        cmd = [sys.executable, os.path.join(PP, "remove_contiguous_samples.py"),
               infile, outfile, "3", "5", "7", "11", "16", "19", "21", pivot,
               "mm+bulges"]
        r = subprocess.run(cmd, cwd=PP, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, "merge failed: %s" % r.stderr)
        self.assertEqual(r.stderr, "", "the merge must write nothing on stderr")

        def _read(path):
            with open(path) as fh:
                lines = [l.rstrip("\n") for l in fh]
            return [l.split("\t") for l in lines[1:] if l.strip()]

        return _read(outfile), _read(outfile + ".discarded_samples")

    def test_tier_path_gives_the_reference_search_row_again(self):
        """The premise: on the tiers, the reference row of a variant window is the
        reference search's row, column for column, once both are written."""
        tier = _tier_rows(["S:1|0"])
        tier_ref = [r for r in tier if r[19] == 55]
        self.assertEqual(len(tier_ref), 1, "one reference row per variant window")
        self.assertEqual(_on_disk(tier_ref), _on_disk([_reference_search_row(1)]))

    def test_one_window_on_the_tiers_lists_its_reference_alignment_once(self):
        """One dataset, one variant window: the tiers give the variant row and the
        window's reference row, and the reference search gives that reference row as
        well. The merge must list it once and count one other alignment, as it does on
        the dictionaries (where only the reference search gives it)."""
        tier = _on_disk(_tier_rows(["S:1|0"]) + [_reference_search_row(1)])
        variant = [r for r in tier if r[SNP_COL] != "n"]
        self.assertEqual(len(variant), 1)
        dicts = variant + _on_disk([_reference_search_row(1)])
        for pivot in ("score", "total"):
            best, disc = self._merge(tier, pivot)
            d_best, d_disc = self._merge(dicts, pivot)
            self.assertEqual(len(best), 1)
            self.assertEqual(best[0][SNP_COL], variant[0][SNP_COL],
                             "the carried variant stays the best alignment")
            refs = [r for r in best + disc if r[SNP_COL] == "n"]
            self.assertEqual(len(refs), 1,
                             "%s: the reference alignment is listed %d times"
                             % (pivot, len(refs)))
            self.assertEqual([r[SEQ_IN_CLUSTER_COL] for r in best + disc],
                             ["1", "0"],
                             "%s: the best row counts the one other alignment, the "
                             "alternative alignment none" % pivot)
            self.assertEqual((best, disc), (d_best, d_disc),
                             "%s: the tiers must give what the dictionaries give" % pivot)

    def test_two_datasets_list_each_reference_alignment_once(self):
        """Two variant datasets: each one's post-analysis gives the reference search's
        rows, so the merge gets every reference row twice. It must list it once; the
        variant row both datasets give is merged into one, with both carriers."""
        rows = _on_disk(
            _tier_rows(["S:1|0"]) + [_reference_search_row(1)]
            + _tier_rows(["S:1|0"], sample_prefix="B") + [_reference_search_row(1)])
        best, disc = self._merge(rows)
        self.assertEqual(len(best), 1)
        self.assertEqual(best[0][SAMPLES_COL], "BS,S",
                         "the variant row of both datasets is one row with both "
                         "carriers")
        self.assertEqual([r[SNP_COL] for r in disc], ["n"],
                         "the reference alignment is listed once, got %r"
                         % [r[DNA_COL] for r in disc])
        self.assertEqual([r[SEQ_IN_CLUSTER_COL] for r in best + disc], ["1", "0"])

    def test_two_reference_alignments_are_both_kept(self):
        """Only identical rows are one alignment: a reference row that differs in any
        column (here a 2-mismatch alignment beside the 1-mismatch one) is kept and
        counted."""
        other = _reference_search_row(2)
        rows = _on_disk([_reference_search_row(1), _reference_search_row(1), other])
        best, disc = self._merge(rows)
        self.assertEqual([r[DNA_COL] for r in best + disc],
                         ["c" + "A" * 20 + "GG", "cc" + "A" * 19 + "GG"])
        self.assertEqual([r[SEQ_IN_CLUSTER_COL] for r in best + disc], ["1", "0"])

    def _lists_of_two_datasets(self, position, cb_copies):
        """The rows of one site in the three criteria's lists, as the merge step's sort
        leaves them: the site's 1-mismatch alignment twice (each dataset's post-analysis
        gives it), then a 2-mismatch alignment twice and a 3-mismatch alignment twice.
        The CFD and the fewest lists' copies are equal. ``cb_copies`` gives the two
        CRISPR-Bulge values of the 1-mismatch alignment, highest first, as the score and
        the REF score of its two copies in the CRISPR-Bulge list."""
        def row(n_mismatches, score):
            r = _on_disk([_reference_search_row(n_mismatches)], position)[0]
            r[CFD_COL], r[CFD_COL + 1] = score, score
            return r

        same = [row(1, "0.500"), row(1, "0.500"), row(2, "0.300"), row(2, "0.300"),
                row(3, "0.100"), row(3, "0.100")]
        cb = [row(1, cb_copies[0]), row(1, cb_copies[1]), row(2, "3e-24"), row(2, "3e-24"),
              row(3, "1e-25"), row(3, "1e-25")]
        return same, same, cb

    def test_copies_scored_differently_are_one_alignment_in_the_three_lists(self):
        """Each variant dataset's post-analysis scores the reference search's rows
        itself, and CRISPR-Bulge can give one alignment two values (a real two-dataset
        search: 4.3390600000000004e-23 and 4.33904e-23). Keeping the copy in that list
        alone gives its discarded rows one row more, and as the pipeline pastes the
        three lists of discarded rows by position (pr -m), every row from there on
        pairs a CFD alignment with the next one's CRISPR-Bulge alignment, at the end of
        a cluster with another site's."""
        lists = [[], [], []]
        for position, cb_copies in ((100, ("4.3390600000000004e-23", "4.33904e-23")),
                                    (5000, ("2e-22", "2e-22"))):
            for rows, site_rows in zip(lists, self._lists_of_two_datasets(position, cb_copies)):
                rows.extend(site_rows)
        (cfd_best, cfd_disc), (few_best, few_disc), (cb_best, cb_disc) = (
            self._merge(lists[0], "score"),
            self._merge(lists[1], "total"),
            self._merge(lists[2], "score"),
        )
        self.assertEqual(len(cb_disc), len(cfd_disc),
                         "the CRISPR-Bulge list keeps a copy the others drop")
        self.assertEqual(len(few_disc), len(cfd_disc))
        self.assertEqual([r[DNA_COL] for r in cb_disc], [r[DNA_COL] for r in cfd_disc],
                         "row by row the lists must give the same alignment")

        # the pipeline's own paste: every row's blocks at one site
        paths = []
        for name, disc in (("cfd", cfd_disc), ("few", few_disc), ("cb", cb_disc)):
            paths.append(os.path.join(self.d, name + ".disc"))
            with open(paths[-1], "w") as fh:
                fh.write("".join("\t".join(r) + "\n" for r in disc))
        pasted = subprocess.run(["pr", "-m", "-t", "-J"] + paths, capture_output=True,
                                text=True, check=True).stdout.splitlines()
        self.assertEqual(len(pasted), 4)
        n = len(cfd_disc[0])
        for line in pasted:
            cols = line.split("\t")
            self.assertEqual(len(cols), 3 * n, "a row lacks a block: %r" % line)
            self.assertEqual(cols[5], cols[n + 5], "the fewest block is another site's")
            self.assertEqual(cols[5], cols[2 * n + 5], "the CRISPR-Bulge block is another site's")

        for best, disc in ((cfd_best, cfd_disc), (few_best, few_disc), (cb_best, cb_disc)):
            self.assertEqual([r[SEQ_IN_CLUSTER_COL] for r in best], ["2", "2"],
                             "each site's best row counts its two other alignments")
            self.assertEqual([r[SEQ_IN_CLUSTER_COL] for r in disc], ["1"] * 4)
        self.assertEqual(cb_best[0][CFD_COL], "4.3390600000000004e-23",
                         "the copy kept is the one with the higher score")

    def test_a_copy_of_another_alignment_is_still_another_alignment(self):
        """Two rows with one score are two alignments when any other column differs:
        leaving the scores out of the comparison merges nothing else."""
        one = _on_disk([_reference_search_row(1)])[0]
        other = list(one)
        other[DNA_COL] = "g" + other[DNA_COL][1:]  # another base at the mismatch
        best, disc = self._merge([one, other])
        self.assertEqual(sorted(r[DNA_COL] for r in best + disc),
                         sorted([one[DNA_COL], other[DNA_COL]]))


if __name__ == "__main__":
    unittest.main()

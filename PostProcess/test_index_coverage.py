#!/usr/bin/env python3
"""An index that lacks a chromosome's trees is not searched as a whole one.

A search took an index that holds one .bin for a whole one: a partial copy, a build an
earlier release left stopped between chromosomes, or an index of a genome folder that
gained a FASTA since, lost that chromosome's sites with exit 0 and no message (flows OF-6:
13 rows against 19 on the synthetic genome; 996 against 2731 on real chr21+chr22 without
chr22's trees). index_coverage.py says which FASTA lacks its trees, by buildTST's own rule
for a site, so that a healthy index -- whose FASTAs with no site of its PAM have no tree --
is never refused; the search script runs it before it opens an index, and build-index-only
on an index it finds in place.

Stdlib only. The trees are fabricated as buildTST names and heads them
(``<PAM>_<contig>_<part>.bin``: int numLeaves, int len_guide_used); the shipped search
script runs with a stub engine.
"""

import io
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import crispritz_bins  # noqa: E402
import index_coverage  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCAC"
SITE = GUIDE + "TGG"  # a site of NRG, and of NGG
FILL = "A" * 60  # (no site of NRG, NGG or TTTV on either strand: TheFixture)


def fasta(folder, name, seq, header=None):
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, f"{name}.fa"), "w") as fh:
        fh.write(f">{header or name}\n{seq}\n")


def tree(index_dir, pam, contig, part, leaves=1, guide=20):
    os.makedirs(index_dir, exist_ok=True)
    with open(os.path.join(index_dir, f"{pam}_{contig}_{part}.bin"), "wb") as fh:
        fh.write(struct.pack("=ii", leaves, guide) + b"\0" * 8)


class _Layout(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="coverage_")
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.genome = os.path.join(self.root, "Genomes", "g")
        self.lib = os.path.join(self.root, "genome_library")
        for c in ("chr1", "chr2", "chr3"):
            fasta(self.genome, c, FILL + SITE + FILL)

    def index(self, name, contigs=("chr1", "chr2", "chr3"), pam=None, record=None):
        folder = os.path.join(self.lib, name)
        for c in contigs:
            tree(folder, pam or name.split("_")[0], c, 1)
        if record:
            with open(os.path.join(folder, ".pam_build"), "w") as fh:
                fh.write(record)
        return folder


class TheFixture(_Layout):
    def test_the_filler_holds_no_site_and_the_site_one(self):
        for pam in ("NRG", "NGG"):
            self.assertFalse(crispritz_bins.holds_target(FILL * 3, pam, 23, False, 2))
            self.assertTrue(crispritz_bins.holds_target(FILL + SITE + FILL, pam, 23, False, 2))
        self.assertFalse(crispritz_bins.holds_target(FILL * 3, "TTTV", 24, True, 2))


class AReferenceIndex(_Layout):
    def test_a_whole_one_has_nothing_missing(self):
        self.assertEqual(index_coverage.problems(self.index("NRG_2_g"), self.genome), [])

    def test_a_fasta_with_sites_and_no_tree_is_missing(self):
        index = self.index("NRG_2_g", ("chr1", "chr3"))
        self.assertEqual(index_coverage.problems(index, self.genome), [("chr2.fa", "no tree")])

    def test_one_added_to_the_genome_folder_since_is_missing_too(self):
        index = self.index("NRG_2_g")
        fasta(self.genome, "chr4", FILL + SITE)
        self.assertEqual(index_coverage.problems(index, self.genome), [("chr4.fa", "no tree")])

    def test_a_fasta_with_no_site_of_its_pam_needs_none(self):
        # buildTST writes none for it: all N, shorter than one window, or no PAM
        index = self.index("NRG_2_g")
        fasta(self.genome, "chrN", "N" * 500)
        fasta(self.genome, "chrS", "ACGTACGTAC")
        fasta(self.genome, "chrA", "A" * 500)
        fasta(self.genome, "chrF", FILL * 4)
        self.assertEqual(index_coverage.problems(index, self.genome), [])

    def test_a_site_too_near_an_n_run_for_the_window_needs_none(self):
        index = self.index("NRG_2_g")
        fasta(self.genome, "chrE", "N" * 30 + SITE + "N" * 30)  # (window: guide + PAM + N)
        self.assertEqual(index_coverage.problems(index, self.genome), [])

    def test_a_gap_in_a_contigs_trees_is_missing(self):
        index = self.index("NRG_2_g", ("chr1", "chr3"))
        tree(index, "NRG", "chr2", 1)
        tree(index, "NRG", "chr2", 3)
        self.assertEqual(index_coverage.problems(index, self.genome),
                         [("chr2.fa", "trees 2 of 3")])

    def test_a_full_last_tree_is_missing_the_next_only_where_the_contig_holds_more(self):
        # (a tree holds TARGETS_A_TREE targets, made as many as chr2 holds here)
        full = crispritz_bins.count_targets(FILL + SITE + FILL, "NRG", 23, False, 2)
        self.assertGreater(full, 1)
        index = self.index("NRG_2_g", ("chr1", "chr3"))
        tree(index, "NRG", "chr2", 1, leaves=full)
        old = crispritz_bins.TARGETS_A_TREE
        self.addCleanup(setattr, crispritz_bins, "TARGETS_A_TREE", old)
        crispritz_bins.TARGETS_A_TREE = full
        self.assertEqual(index_coverage.problems(index, self.genome), [])
        fasta(self.genome, "chr2", FILL + SITE + FILL + SITE)
        got = index_coverage.problems(index, self.genome)
        self.assertEqual([f for f, _ in got], ["chr2.fa"])
        self.assertIn("its last tree is full", got[0][1])

    def test_a_contig_named_with_a_description_or_an_underscore(self):
        # buildTST names a tree after the whole header line
        index = self.index("NRG_2_g")
        fasta(self.genome, "chrUn_KI270302v1", FILL + SITE, header="chrUn_KI270302v1")
        fasta(self.genome, "chr5", FILL + SITE, header="chr5 Homo sapiens chromosome 5")
        tree(index, "NRG", "chrUn_KI270302v1", 1)
        tree(index, "NRG", "chr5 Homo sapiens chromosome 5", 1)
        self.assertEqual(index_coverage.problems(index, self.genome), [])

    def test_an_all_n_index_is_read_in_its_own_geometry(self):
        # every base is a site of NNN: a FASTA long enough for a window holds one
        index = self.index("NNN_3_g", ("chr1", "chr2"))
        self.assertEqual(index_coverage.problems(index, self.genome), [("chr3.fa", "no tree")])
        fasta(self.genome, "chrS", "ACGTACGTACGTACGTACGTACG")  # (shorter than 20 + 3 + 3)
        self.assertEqual([f for f, _ in index_coverage.problems(index, self.genome)], ["chr3.fa"])

    def test_a_5_prime_index_by_its_record_or_the_search(self):
        for c in ("chr1", "chr2", "chr3"):
            fasta(self.genome, c, FILL + "TTTA" + GUIDE + FILL)
        recorded = self.index("TTTV_2_g", ("chr1", "chr2"), record="23bp-TTTV-Cas12a 24 -4\n")
        self.assertEqual(index_coverage.problems(recorded, self.genome), [("chr3.fa", "no tree")])
        unrecorded = self.index("TTTV_3_g", ("chr1", "chr2"))
        self.assertEqual(index_coverage.problems(unrecorded, self.genome, True),
                         [("chr3.fa", "no tree")])
        fasta(self.genome, "chr3", FILL + "TTTA")  # (no guide after it: no 5' site)
        self.assertEqual(index_coverage.problems(recorded, self.genome), [])

    def test_a_folder_renamed_for_another_pam_is_said_so(self):
        index = self.index("NRG_2_g", pam="NGG")
        got = index_coverage.problems(index, self.genome)
        self.assertEqual(got, [("", "its trees are named for NGG, not NRG: the folder was "
                                    "renamed")])

    def test_an_index_whose_guide_length_cannot_be_read_is_left_alone(self):
        index = os.path.join(self.lib, "NRG_2_g")
        os.makedirs(index)
        with open(os.path.join(index, "NRG_chr1_1.bin"), "wb") as fh:
            fh.write(b"")
        self.assertEqual(index_coverage.problems(index, self.genome), [])


class AVariantIndex(_Layout):
    def test_it_answers_for_the_chromosomes_of_its_dataset_only(self):
        index = self.index("NRG_2_g+ds", ("chr1",))
        self.assertEqual(index_coverage.problems(index, self.genome, chroms=["chr1"]), [])
        self.assertEqual(index_coverage.problems(index, self.genome, chroms=["chr1", "chr2"]),
                         [("chr2.fa", "no tree")])

    def test_a_chromosome_the_genome_folder_has_no_fasta_of_is_left_alone(self):
        index = self.index("NRG_2_g+ds", ("chr1",))
        self.assertEqual(index_coverage.problems(index, self.genome, chroms=["chr1", "chrY"]), [])

    def test_its_chromosomes_are_taken_from_the_vcf_names_as_the_search_takes_them(self):
        vcfs = os.path.join(self.root, "VCFs", "ds")
        os.makedirs(vcfs)
        for f in ("ALL.chr22.phased.vcf.gz", "ds.chrX.vcf.gz", "ds.chrX.vcf.gz.tbi", "notes.txt"):
            open(os.path.join(vcfs, f), "w").close()
        self.assertEqual(index_coverage.vcf_chroms(vcfs), ["chr22", "chrX"])


class TheCommandLine(_Layout):
    def run_main(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            rc = index_coverage.main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_a_whole_index_says_nothing(self):
        self.assertEqual(self.run_main(self.index("NRG_2_g"), "3", self.genome), (0, "", ""))

    def test_one_that_is_not_says_which_fasta_and_what_to_do_on_stderr(self):
        index = self.index("NRG_2_g", ("chr1",))
        rc, out, err = self.run_main(index, "3", self.genome)
        self.assertEqual((rc, out), (1, ""))
        self.assertEqual(err, f"ERROR: index 'NRG_2_g' ({index}) lacks trees of FASTAs that hold "
                              "sites of its PAM: chr2.fa (no tree); chr3.fa (no tree): a search of "
                              "it would lose them with no error. It is a partial copy, or a build "
                              "stopped partway, or the genome folder has changed since it was "
                              f"built. Remove it (rm -rf '{index}') and run this search again, "
                              "which builds it anew where it can, or download it again or build "
                              "it with crisprme.py build-index-only\n")

    def test_one_it_cannot_read_is_a_note_on_stdout(self):
        index = self.index("NRG_2_g")
        rc, out, err = self.run_main(index, "3", os.path.join(self.root, "no_such_genome"))
        self.assertEqual((rc, err), (0, ""))
        self.assertIn("Note: could not check that", out)


@unittest.skipUnless(shutil.which("bash") and shutil.which("python"),
                     "the shipped script needs bash and python on PATH")
class TheSearchScript(unittest.TestCase):
    """The shipped search script, with a stub engine: it checks the index it is about to
    search with bulges, and goes on without them."""

    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="coverage_script_")
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.home = os.path.join(self.root, "home")
        self.genome = os.path.join(self.home, "Genomes", "g")
        for c in ("chr1", "chr2"):
            fasta(self.genome, c, FILL + SITE + FILL)
        self.lib = os.path.join(self.home, "genome_library")
        tree(os.path.join(self.lib, "NRG_2_g"), "NRG", "chr1", 1)
        stub = os.path.join(self.root, "bin")
        os.makedirs(stub)
        with open(os.path.join(stub, "crispritz.py"), "w") as fh:
            fh.write("#!/bin/sh\necho \"$*\" >>\"$STUB_LOG\"\n"
                     "printf '#Bulge type\\tcrRNA\\n' >\"$5.targets.txt\"\nexit 0\n")
        os.chmod(os.path.join(stub, "crispritz.py"), 0o755)
        self.env = dict(os.environ, PATH=stub + os.pathsep + os.environ.get("PATH", ""),
                        STUB_LOG=os.path.join(self.root, "calls.txt"), LC_ALL="C",
                        CRISPRME_NO_THREAD_GUARD="1")

    def run_script(self, bulges):
        out = os.path.join(self.root, f"out{bulges}")
        os.makedirs(out)
        pam = os.path.join(self.root, "20bp-NRG-SpCas9.txt")
        with open(pam, "w") as fh:
            fh.write("N" * 20 + "NRG 3\n")
        guides = os.path.join(self.root, "guides.txt")
        with open(guides, "w") as fh:
            fh.write(GUIDE + "NNN\n")
        empty = os.path.join(self.root, "empty.txt")
        open(empty, "w").close()
        args = [self.genome, "_", guides, pam, empty, "_", bulges + 1, 4, bulges, bulges, 3,
                out, HERE, 1, self.home, empty, "_", 1, 0, "_", "mm+bulges", "mm+bulges,mm",
                "False", "PASS,.", "_", 6, 0]
        with open(os.path.join(out, "log_error.txt"), "w") as err:
            done = subprocess.run(["bash", os.path.join(HERE,
                                                        "submit_job_automated_new_multiple_vcfs.sh"),
                                   *map(str, args)], cwd=self.root, env=self.env,
                                  stdout=subprocess.PIPE, stderr=err, text=True, timeout=300)
        with open(os.path.join(out, "log_error.txt")) as fh:
            return done.returncode, done.stdout, fh.read()

    def calls(self):
        try:
            with open(self.env["STUB_LOG"]) as fh:
                return fh.read().splitlines()
        except OSError:
            return []

    def test_a_search_with_bulges_of_an_index_lacking_a_chromosome_stops_before_it(self):
        rc, out, err = self.run_script(1)
        self.assertEqual(rc, 1, out)
        self.assertIn("ERROR: index 'NRG_2_g'", err)
        self.assertIn("chr2.fa (no tree)", err)
        self.assertFalse(any(c.startswith("search") for c in self.calls()), self.calls())

    def test_a_search_without_bulges_scans_the_genome_and_goes_on(self):
        rc, out, err = self.run_script(0)
        self.assertNotIn("ERROR: index", err)
        self.assertTrue(any(c.startswith("search") and " -r " in f" {c} " for c in self.calls()),
                        self.calls())


if __name__ == "__main__":
    unittest.main(verbosity=2)

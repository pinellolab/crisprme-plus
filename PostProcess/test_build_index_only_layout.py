#!/usr/bin/env python3
"""build-index-only refuses a genome folder it cannot index as a search reads it, and
an index it built that holds no .bin (Gate B walk P5).

index-genome indexes the ``*.fa`` files at the top of the genome folder only. A UCSC
chromFa archive unpacked with its chroms/ folder (as the web Settings genome upload
left one) gave an index with no .bin, which build-index-only put in place, and the web
showed as done, and every search then refused; a FASTA file of several sequences gave
an index of a genome every search refuses to use (genome_layout.py), and a FASTA file of
one sequence named otherwise one whose off-targets every search found and then dropped
(round 3 review, F2). Now the run ends with an ERROR that says what to do, before the
engine runs, and a reference or variant index the engine wrote no .bin into is never
put in place; nor is one with no .bin that an earlier run left taken for present
(final review round 1, F2).

These run the real crisprme.py build-index-only against the stub engine of
test_build_index_only_atomic (which writes a .bin for each contig that holds a site,
none for a contig of N alone), the way its own tests do.
"""

import os
import unittest

import test_build_index_only_atomic as bio

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule


class AGenomeFolderItCannotIndexIsRefused(bio._Fixture):
    def genome(self, *parts):
        return self.path("Genomes", "hg38", *parts)

    def assert_refused(self, rc, out, err, *says):
        self.assertEqual(rc, 1, out + err)
        for text in says:
            self.assertIn(text, err)
        self.assertTrue(err.startswith("Error: "), err)
        self.assertEqual(self.calls(), [])  # the engine never ran
        self.assertFalse(os.path.exists(self.path("genome_library")))
        self.assertEqual(self.leftovers(), [])

    def test_fasta_files_in_a_folder_of_the_genome(self):
        # a UCSC hg38.chromFa.tar.gz unpacked as it is: chroms/chr1.fa ...
        os.makedirs(self.genome("chroms"))
        for name in ("chr1.fa", "chr2.fa"):
            os.rename(self.genome(name), self.genome("chroms", name))
        rc, out, err = self.run_bio(vcf=False)
        self.assert_refused(rc, out, err, "holds no .fa file at its top",
                            "(it holds chroms/chr1.fa, chroms/chr2.fa)",
                            "move the .fa files up out of their folder")

    def test_fasta_files_index_genome_does_not_read(self):
        for label, rename in ((".fasta", lambda n: n[:-3] + ".fasta"),
                              ("gzipped", lambda n: n + ".gz")):
            with self.subTest(label):
                self.fresh()
                for name in ("chr1.fa", "chr2.fa"):
                    os.rename(self.genome(name), self.genome(rename(name)))
                rc, out, err = self.run_bio(vcf=False)
                self.assert_refused(rc, out, err, "holds no .fa file at its top",
                                    rename("chr1.fa"))

    def test_fasta_file_of_several_sequences(self):
        # a genome uploaded as one FASTA (the web Settings upload kept it so)
        with open(self.genome("hg38.fa"), "w") as out:
            for name in ("chr1.fa", "chr2.fa"):
                with open(self.genome(name)) as fh:
                    out.write(fh.read())
                os.remove(self.genome(name))
        rc, out, err = self.run_bio(vcf=False)
        # and to remove an index built from it before: it holds the sites of chr2 under
        # chr1's name, and a later run would take it as it is (final review round 1, F2)
        self.assert_refused(rc, out, err, "hg38.fa: holds 2 sequences (chr1, chr2)",
                            "python PostProcess/genome_layout.py --split",
                            "\nThen remove each index built from the folder while it held "
                            "that file (genome_library/*_hg38, genome_library/*_hg38+*, and "
                            "the variant genomes Genomes/hg38+*)",
                            "and run build-index-only again.",
                            "remove this one there, with the indexes built on it")

    def test_fasta_file_not_named_after_its_sequence(self):
        # one sequence in a file named otherwise (the web Settings upload of walkone.fa
        # holding >chr1 kept it so): the index was built, and every search of the folder
        # found its off-targets and then dropped them all, a search taking each file for
        # the sequence it is named after (round 3 review, F2)
        misnamed, named = self.genome("hg38_chr2.fa"), self.genome("chr2.fa")
        os.rename(named, misnamed)
        rc, out, err = self.run_bio(vcf=False)
        self.assert_refused(rc, out, err,
                            f"{misnamed} holds the sequence chr2, but a search takes each "
                            "FASTA file for the sequence it is named after (hg38_chr2), so "
                            "every search of this folder finds the off-targets of chr2 and "
                            "then drops them all.",
                            f"(mv '{misnamed}' '{named}', or python PostProcess/genome_layout.py "
                            f"--split {self.genome()}, which renames them all)",
                            "remove this one there, with the indexes built on it")
        # named after its sequence, the folder builds
        os.rename(misnamed, named)
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(sorted(f for f in os.listdir(self.lib(bio.REF)) if f.endswith(".bin")),
                         ["NRG_chr1_1.bin", "NRG_chr2_1.bin"])

    def test_an_index_with_no_bin_is_not_put_in_place(self):
        # no sequence holds a site (a contig of N alone): the engine writes no .bin
        # and exits 0, and every search refuses an index with none
        for name in ("chr1", "chr2"):
            with open(self.genome(name + ".fa"), "w") as fh:
                fh.write(f">{name}\n{'N' * 120}\n")
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 1, out + err)
        self.assertIn("Reference genome indexing failed: the engine exited 0 but wrote no "
                      ".bin", err)
        self.assertIn(f"nothing was put in place (expected {self.lib(bio.REF)})", err)
        self.assertFalse(os.path.exists(self.lib(bio.REF)))
        self.assertEqual(self.leftovers(), [])

    def test_a_reference_index_with_no_bin_is_not_taken_for_present(self):
        # one built from a folder with no .fa at its top before it was refused (an
        # upload of a chromFa archive kept chroms/): the run said 'Reference index
        # already present' and the web showed it done, though no search uses it
        os.makedirs(self.lib(bio.REF))
        with open(os.path.join(self.lib(bio.REF), ".pam_build"), "w") as fh:
            fh.write("20bp-NRG-SpCas9.txt\n")
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 1, out + err)
        self.assertEqual(err, f"Error: The reference index {self.lib(bio.REF)} holds no .bin "
                         "file, as one built from a genome folder with no .fa file at its "
                         "top does, and no search uses it: remove it (rm -rf "
                         f"'{self.lib(bio.REF)}', or remove the index in the web Settings "
                         "page) and run build-index-only again.\n")
        self.assertNotIn("already present", out)
        self.assertEqual(self.calls(), [])  # the engine never ran
        self.assertEqual(os.listdir(self.lib(bio.REF)), [".pam_build"])
        # removed, it is built
        os.remove(os.path.join(self.lib(bio.REF), ".pam_build"))
        os.rmdir(self.lib(bio.REF))
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, out + err)
        self.assertTrue(any(f.endswith(".bin") for f in os.listdir(self.lib(bio.REF))))

    def test_a_genome_of_one_fasta_a_sequence_builds_as_before(self):
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(err, "")
        self.assertEqual(sorted(f for f in os.listdir(self.lib(bio.REF)) if f.endswith(".bin")),
                         ["NRG_chr1_1.bin", "NRG_chr2_1.bin"])


if __name__ == "__main__":
    unittest.main()

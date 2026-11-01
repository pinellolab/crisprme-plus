#!/usr/bin/env python3
"""build-index-only refuses a genome folder it cannot index as a search reads it, and
an index it built that holds no .bin (Gate B walk P5).

index-genome indexes the ``*.fa`` files at the top of the genome folder only. A UCSC
chromFa archive unpacked with its chroms/ folder (as the web Settings genome upload
left one) gave an index with no .bin, which build-index-only put in place, and the web
showed as done, and every search then refused; a FASTA file of several sequences gave
an index of a genome every search refuses to use (genome_layout.py). Now the run ends
with an ERROR that says what to do, before the engine runs, and a reference or
variant index the engine wrote no .bin into is never put in place.

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
        self.assert_refused(rc, out, err, "hg38.fa: holds 2 sequences (chr1, chr2)",
                            "python PostProcess/genome_layout.py --split",
                            "Then run build-index-only again.")

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

    def test_a_genome_of_one_fasta_a_sequence_builds_as_before(self):
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, out + err)
        self.assertEqual(err, "")
        self.assertEqual(sorted(f for f in os.listdir(self.lib(bio.REF)) if f.endswith(".bin")),
                         ["NRG_chr1_1.bin", "NRG_chr2_1.bin"])


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""The search script's listing of the reference genome's chromosomes (real_chroms).

submit_job_automated_new_multiple_vcfs.sh lists the genome folder's FASTA files and takes
each file's name, less its extension, as a chromosome's. It took the name with
``$(basename $file_chr)``: a basename process per file, about 1 s for hg38's 455 contigs
in the image, and, the path unquoted, word splitting and pathname expansion on it, so
that a folder or a file whose name holds a space, or a file whose name is a pattern that
matches another file's, gave wrong names. PR 107 (Tim Barry) takes the name with
``${file_chr##*/}``: the same name for every other path, from the shell itself.

The loop is run here verbatim from the script, in bash, on a folder of such names, with a
``basename`` on PATH that fails, so a loop that runs it again lists wrong names.
"""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parent / "submit_job_automated_new_multiple_vcfs.sh"


def _listing_loop():
    """The script's loop that fills real_chroms, verbatim, from its `unset real_chroms`
    to its `done`."""
    text = _SCRIPT.read_text()
    start = text.index("\tunset real_chroms\n")
    end = text.index("\n\tdone\n", start) + len("\n\tdone\n")
    loop = text[start:end]
    assert 'for file_chr in "$ref_folder"/*.fa; do' in loop, loop
    return loop


# hg38's kinds of contig, and names with dots, dashes, spaces and pattern characters
NAMES = [
    "chr1.fa", "chr10.fa", "chrX.fa", "chrM.fa", "chrEBV.fa",
    "chr1_KI270706v1_random.fa", "chrUn_GL000195v1.fa", "chr19_KI270938v1_alt.fa",
    "KI270728.1.fa", "GCA_000001405.15_GRCh38.fa", "chr12.fa.fa",
    "chr-1.fa", "-chr2.fa", "--chr2b.fa", "HLA-A.fa",
    "chr 3.fa", " chr4.fa", "chr  6.fa",
    "HLA-A*01:01:01:01.fa", "HLA-A01:01:01:01.fa", "chr[8].fa", "chr8.fa", "chr?.fa",
]


class TheChromosomeListing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="fasta_listing_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        bin_dir = os.path.join(self.tmp, "bin")
        os.makedirs(bin_dir)
        fake = os.path.join(bin_dir, "basename")
        with open(fake, "w") as fh:
            fh.write("#!/bin/sh\necho basename-was-run\nexit 1\n")
        os.chmod(fake, 0o755)
        self.path = bin_dir + os.pathsep + os.environ.get("PATH", "")

    def listing(self, folder):
        """(real_chroms, what the loop printed) for ref_folder=folder."""
        out = os.path.join(self.tmp, "real_chroms")
        script = _listing_loop() + '\nprintf "%s\\0" "${real_chroms[@]}" >"$out"\n'
        proc = subprocess.run(
            ["bash", "-c", script, "bash"], capture_output=True, text=True,
            env=dict(os.environ, PATH=self.path, LC_ALL="C", ref_folder=folder, out=out))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            listed = fh.read().split("\0")[:-1]
        return listed, proc.stdout

    def check(self, folder):
        os.makedirs(folder)
        for name in NAMES:
            Path(folder, name).touch()
        Path(folder, "chr1.fa.fai").touch()   # not a FASTA: not listed
        listed, echoed = self.listing(folder)
        want = sorted(n[: -len(".fa")] for n in NAMES)
        self.assertEqual(sorted(listed), want)
        self.assertEqual(len(listed), len(NAMES))
        self.assertNotIn("basename-was-run", echoed)

    def test_names_with_dots_dashes_spaces_and_patterns(self):
        self.check(os.path.join(self.tmp, "hg38"))

    def test_in_a_folder_whose_path_holds_them_too(self):
        self.check(os.path.join(self.tmp, "my genome", "hg38.p14-alt [v2]"))

    def test_hg38_names_are_unchanged(self):
        # the names a genome folder of the pipeline holds (hg38 with alts: 455 contigs, each
        # [A-Za-z0-9_]): the chromosome is the file's name less its last extension
        folder = os.path.join(self.tmp, "hg38_alt")
        os.makedirs(folder)
        names = ["chr1", "chr22_KI270879v1_alt", "chrUn_KI270742v1", "chr17_GL000205v2_random"]
        for name in names:
            Path(folder, name + ".fa").touch()
        listed, echoed = self.listing(folder)
        self.assertEqual(listed, sorted(names))
        self.assertEqual(echoed.split("\n")[:-1], sorted(names))


if __name__ == "__main__":
    unittest.main()

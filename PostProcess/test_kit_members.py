#!/usr/bin/env python3
"""An enriched-genome kit's member list (scripts/kit_members.py; PP-4 of CRISPRme 2.10.0's release).

A kit, build-inputs/<dataset>.tar.gz, is packed by the rebuild with tar, and its tarball is not
byte-reproducible: tar records each file's time, owner and listing order, so the same files packed
twice give two tarballs, two SHA-256. The member list, the md5 of each file by its path in the kit,
is the same for both, read from the folders the kit packs or from either tarball, and checks a kit
once unpacked (md5sum -c); a changed file changes it.

Hermetic: a small kit packed here as the rebuild's row packs one (tar, then gzip for pigz).
"""

import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "kit_members.py")
DS = "hg38_ds"


def _write(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)


class AKitsMemberList(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="kitmembers_")
        self.addCleanup(shutil.rmtree, self.root, True)
        self.work = os.path.join(self.root, "work")
        enriched = os.path.join(self.work, "Genomes", f"hg38+{DS}")
        fake = os.path.join(self.work, "Genomes", f"hg38+{DS}_INDELS")
        _write(os.path.join(enriched, "chr1.enriched.fa"), b">chr1\nACGTRYACGT\n")
        _write(os.path.join(enriched, "chr2.enriched.fa"), b">chr2\nTTTTQAAAA\n")
        _write(os.path.join(fake, f"fake_{DS}_chr1", "fakechr1.fa"), b">fakechr1\nACGTNACGT\n")
        _write(os.path.join(fake, f"fake_{DS}_chr2", "fakechr2.fa"), b">fakechr2\n\n")

    def run_it(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)

    def pack(self, name, mtime):
        """The kit, as the rebuild's row packs it: tar of the two folders from the working
        directory, compressed (gzip standing for pigz); every file's time set to ``mtime``."""
        for folder, _, files in os.walk(os.path.join(self.work, "Genomes")):
            for f in files:
                os.utime(os.path.join(folder, f), (mtime, mtime))
        path = os.path.join(self.root, name)
        with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=mtime) as out:
            with tarfile.open(fileobj=out, mode="w|") as tf:
                for top in (f"Genomes/hg38+{DS}", f"Genomes/hg38+{DS}_INDELS"):
                    tf.add(os.path.join(self.work, top), arcname=top)
        return path

    def sha256(self, path):
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    def listed(self, path):
        with open(path) as fh:
            return fh.read()

    def test_one_list_for_two_packings_of_the_same_files(self):
        first, second = self.pack("a.tar.gz", 1_700_000_000), self.pack("b.tar.gz", 1_800_000_000)
        self.assertNotEqual(self.sha256(first), self.sha256(second))  # (not byte-reproducible)
        got = self.run_it(self.work, DS)
        self.assertEqual((got.returncode, got.stderr), (0, ""))
        mine = os.path.join(self.work, f"{DS}.members.md5")
        self.assertEqual(got.stdout, f"4 file(s) -> {mine}\n")
        lines = self.listed(mine).splitlines()
        self.assertEqual([ln.split("  ", 1)[1] for ln in lines], [
            f"Genomes/hg38+{DS}/chr1.enriched.fa", f"Genomes/hg38+{DS}/chr2.enriched.fa",
            f"Genomes/hg38+{DS}_INDELS/fake_{DS}_chr1/fakechr1.fa",
            f"Genomes/hg38+{DS}_INDELS/fake_{DS}_chr2/fakechr2.fa"])
        self.assertEqual(lines[0].split("  ")[0], hashlib.md5(b">chr1\nACGTRYACGT\n").hexdigest())
        for kit in (first, second):
            out = kit + ".members.md5"
            got = self.run_it("--kit", kit, "--out", out)
            self.assertEqual(got.returncode, 0, got.stderr)
            self.assertEqual(self.listed(out), self.listed(mine))

    def test_it_checks_a_kit_once_unpacked(self):
        kit = self.pack("a.tar.gz", 1_700_000_000)
        self.assertEqual(self.run_it(self.work, DS).returncode, 0)
        there = os.path.join(self.root, "unpacked")
        os.makedirs(there)
        with tarfile.open(kit) as tf:
            tf.extractall(there, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
        shutil.copy(os.path.join(self.work, f"{DS}.members.md5"), there)
        if shutil.which("md5sum"):
            got = subprocess.run(["md5sum", "-c", f"{DS}.members.md5"], cwd=there,
                                 capture_output=True, text=True)
            self.assertEqual(got.returncode, 0, got.stdout + got.stderr)
            self.assertEqual(got.stdout.count(": OK"), 4)
        # a file of it changed: its line no longer holds
        _write(os.path.join(there, "Genomes", f"hg38+{DS}", "chr2.enriched.fa"), b">chr2\nTTTTTAAAA\n")
        got = self.run_it(there, DS, "--out", os.path.join(self.root, "again.md5"))
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertNotEqual(self.listed(os.path.join(self.root, "again.md5")),
                            self.listed(os.path.join(self.work, f"{DS}.members.md5")))
        if shutil.which("md5sum"):
            got = subprocess.run(["md5sum", "-c", f"{DS}.members.md5"], cwd=there,
                                 capture_output=True, text=True)
            self.assertNotEqual(got.returncode, 0)
            self.assertIn(f"Genomes/hg38+{DS}/chr2.enriched.fa: FAILED", got.stdout)

    def test_a_member_that_is_no_regular_file_writes_no_list(self):
        os.symlink("chr1.enriched.fa",
                   os.path.join(self.work, "Genomes", f"hg38+{DS}", "chrZ.enriched.fa"))
        got = self.run_it(self.work, DS)
        self.assertEqual(got.returncode, 1)
        self.assertIn("1 member(s) that are no regular file (a kit holds none): "
                      f"Genomes/hg38+{DS}/chrZ.enriched.fa", got.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.work, f"{DS}.members.md5")))

    def test_a_folder_the_kit_packs_that_is_not_there(self):
        shutil.rmtree(os.path.join(self.work, "Genomes", f"hg38+{DS}_INDELS"))
        got = self.run_it(self.work, DS)
        self.assertEqual(got.returncode, 1)
        self.assertIn(f"hg38+{DS}_INDELS is not there: the kit packs it", got.stderr)

    def test_its_usage(self):
        for args, why in (((), "give a working directory and a dataset"),
                          ((self.work,), "give a working directory and a dataset"),
                          (("--kit", "x.tar.gz"), "--kit needs --out"),
                          ((self.work, DS, "--kit", "x"), "give a working directory and a dataset")):
            with self.subTest(args=args):
                got = self.run_it(*args)
                self.assertEqual(got.returncode, 2)
                self.assertIn("usage: kit_members.py", got.stderr)
                self.assertIn(why, got.stderr)


if __name__ == "__main__":
    unittest.main()

"""Unit tests for the personal-assembly folder+metadata layout
(``personal_assembly.py``): the one-folder-per-individual + ``metadata.json``
storage used by download_hprc_assembly.py, the assembly-search CLI, and the
web Data Manager.

Pure stdlib (no pandas / dash), so it runs anywhere:

    python -m unittest PostProcess.test_personal_assembly -v
"""
import os
import sys
import tarfile
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import personal_assembly as pa  # noqa: E402


def _make_haplotype(cwd, individual, hap, acc, n_contigs=3):
    """Create one haplotype's files + record it in metadata.json (mimics the
    downloader), returning the assembly_name."""
    hapdir = os.path.join(pa.individual_dir(cwd, individual), hap)
    gdir = os.path.join(hapdir, "genome")
    os.makedirs(gdir, exist_ok=True)
    for i in range(n_contigs):
        with open(os.path.join(gdir, f"chr{i + 1}.fa"), "w") as fh:
            fh.write(f">chr{i + 1}\nACGT\n")
    assembly_name = f"{individual}_{hap}_{acc}"
    chain = f"{assembly_name}_vs_GRCh38.chain.gz"
    ca = f"{assembly_name}.chromAlias.txt"
    with open(os.path.join(hapdir, chain), "w") as fh:
        fh.write("chain")
    with open(os.path.join(hapdir, ca), "w") as fh:
        fh.write("assembly\tucsc\tgenbank\n")
    pa.write_haplotype(
        cwd, individual, hap,
        {
            "assembly_name": assembly_name,
            "genome_dir": f"{hap}/genome",
            "chain": f"{hap}/{chain}",
            "chromalias": f"{hap}/{ca}",
            "n_contigs": n_contigs,
            "fasta_md5": "deadbeef",
        },
        source="HPRC release2",
    )
    return assembly_name


class TestWriteAndResolve(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()

    def test_write_then_resolve_complete(self):
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG01255", "maternal", "CM2")
        r = pa.resolve(self.cwd, "HG01255")
        self.assertIsNotNone(r)
        self.assertEqual(r["layout"], "bundle")
        self.assertTrue(r["complete"])
        self.assertEqual(r["source"], "HPRC release2")
        # cwd-relative resolved paths for the search engine
        self.assertEqual(
            r["paternal"]["genome_path"],
            os.path.join("Assemblies", "HG01255", "paternal", "genome"),
        )
        self.assertTrue(r["paternal"]["present"])
        self.assertTrue(r["maternal"]["present"])

    def test_single_haplotype_is_incomplete(self):
        _make_haplotype(self.cwd, "HG02", "paternal", "CM1")
        r = pa.resolve(self.cwd, "HG02")
        self.assertIsNotNone(r)
        self.assertFalse(r["complete"])
        self.assertIsNone(r["maternal"])

    def test_metadata_merges_haplotypes(self):
        _make_haplotype(self.cwd, "HG03", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG03", "maternal", "CM2")
        meta = pa.read_metadata(self.cwd, "HG03")
        self.assertEqual(set(meta["haplotypes"]), {"paternal", "maternal"})
        self.assertEqual(meta["format"], pa.BUNDLE_FORMAT)

    def test_missing_genome_files_marks_not_present(self):
        _make_haplotype(self.cwd, "HG04", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG04", "maternal", "CM2")
        # remove the paternal contigs -> not present -> incomplete
        gdir = os.path.join(pa.individual_dir(self.cwd, "HG04"), "paternal", "genome")
        for f in os.listdir(gdir):
            os.remove(os.path.join(gdir, f))
        r = pa.resolve(self.cwd, "HG04")
        self.assertFalse(r["paternal"]["present"])
        self.assertFalse(r["complete"])

    def test_resolve_missing_individual_is_none(self):
        self.assertIsNone(pa.resolve(self.cwd, "NOPE"))


class TestDiscoverAndDelete(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()

    def test_discover_lists_only_bundles(self):
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG01255", "maternal", "CM2")
        # a stray dir with no metadata.json must be ignored
        os.makedirs(os.path.join(self.cwd, pa.ASSEMBLIES_DIR, "junk"))
        found = pa.discover_bundles(self.cwd)
        self.assertEqual([f["individual"] for f in found], ["HG01255"])

    def test_delete_target_is_the_whole_folder(self):
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        self.assertEqual(
            pa.delete_target(self.cwd, "HG01255"),
            pa.individual_dir(self.cwd, "HG01255"),
        )

    def test_delete_target_none_for_non_bundle(self):
        self.assertIsNone(pa.delete_target(self.cwd, "NOPE"))


class TestBundleImport(unittest.TestCase):
    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG01255", "maternal", "CM2")

    def test_validate_good_bundle(self):
        self.assertIsNone(pa.validate_bundle_dir(pa.individual_dir(self.cwd, "HG01255")))

    def test_targz_roundtrip(self):
        arc = os.path.join(self.cwd, "b.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            t.add(pa.individual_dir(self.cwd, "HG01255"), arcname="HG01255")
        dest = tempfile.mkdtemp()
        name = pa.import_archive(arc, dest)
        self.assertEqual(name, "HG01255")
        self.assertTrue(pa.resolve(dest, "HG01255")["complete"])

    def test_zip_roundtrip(self):
        src = pa.individual_dir(self.cwd, "HG01255")
        arc = os.path.join(self.cwd, "b.zip")
        with zipfile.ZipFile(arc, "w") as z:
            for root, _dirs, files in os.walk(src):
                for f in files:
                    full = os.path.join(root, f)
                    z.write(full, os.path.join("HG01255", os.path.relpath(full, src)))
        dest = tempfile.mkdtemp()
        self.assertEqual(pa.import_archive(arc, dest), "HG01255")

    def test_reject_multiple_top_dirs(self):
        arc = os.path.join(self.cwd, "bad.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            for n in ("A/x", "B/y"):
                p = os.path.join(self.cwd, n)
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p, "w").write("z")
                t.add(p, arcname=n)
        with self.assertRaises(ValueError):
            pa.import_archive(arc, tempfile.mkdtemp())

    def test_reject_existing_destination(self):
        arc = os.path.join(self.cwd, "b.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            t.add(pa.individual_dir(self.cwd, "HG01255"), arcname="HG01255")
        # importing into the SAME cwd (already has HG01255) must refuse
        with self.assertRaises(ValueError):
            pa.import_archive(arc, self.cwd)


if __name__ == "__main__":
    unittest.main()

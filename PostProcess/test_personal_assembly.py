"""Unit tests for the personal-assembly folder+metadata layout
(``personal_assembly.py``): the one-folder-per-individual + ``metadata.json``
storage used by download_hprc_assembly.py, the assembly-search CLI, and the
web Data Manager.

Pure stdlib (no pandas / dash), so it runs anywhere:

    python -m unittest PostProcess.test_personal_assembly -v
"""
import io
import json
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import unittest
import warnings
import zipfile
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import personal_assembly as pa  # noqa: E402
import safe_tar  # noqa: E402


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


class TestHaplotypeLabel(unittest.TestCase):
    """The name a haplotype genome is shown under (results header, History): every
    bundle haplotype's folder is named "genome", so the folder name says nothing."""

    def setUp(self):
        self.cwd = tempfile.mkdtemp()

    def test_bundle_haplotype_is_its_assembly_name(self):
        name = _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        gdir = os.path.join(pa.individual_dir(self.cwd, "HG01255"), "paternal", "genome")
        self.assertEqual(pa.haplotype_label(gdir, "paternal"), name)
        self.assertEqual(pa.haplotype_label(gdir + os.sep), name)

    def test_bundle_without_assembly_name_is_individual_and_haplotype(self):
        _make_haplotype(self.cwd, "HG02", "maternal", "CM2")
        meta = pa.read_metadata(self.cwd, "HG02")
        del meta["haplotypes"]["maternal"]["assembly_name"]
        with open(pa.metadata_path(self.cwd, "HG02"), "w") as fh:
            json.dump(meta, fh)
        gdir = os.path.join(pa.individual_dir(self.cwd, "HG02"), "maternal", "genome")
        self.assertEqual(pa.haplotype_label(gdir, "maternal"), "HG02_maternal")

    def test_any_other_folder_is_its_name(self):
        legacy = os.path.join(self.cwd, "Genomes", "HG01255_paternal")
        os.makedirs(legacy)
        self.assertEqual(pa.haplotype_label(legacy, "paternal"), "HG01255_paternal")


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
                with open(p, "w") as fh:
                    fh.write("z")
                t.add(p, arcname=n)
        with self.assertRaises(ValueError):
            pa.import_archive(arc, tempfile.mkdtemp())

    def _tar_with_raw_names(self, path, names):
        """Build a .tar.gz whose members have EXACTLY the given (possibly
        malicious) stored names, bypassing tarfile.add's own normalization."""
        with tarfile.open(path, "w:gz") as t:
            for nm in names:
                data = b"x"
                ti = tarfile.TarInfo(name=nm)
                ti.size = len(data)
                t.addfile(ti, io.BytesIO(data))

    def test_reject_absolute_path_member(self):
        # A member with an absolute path (e.g. /etc/evil) must be refused BEFORE
        # extraction -- it must not survive an lstrip('/') and get written out.
        arc = os.path.join(self.cwd, "abs.tar.gz")
        self._tar_with_raw_names(arc, ["/malicious/metadata.json"])
        with self.assertRaises(ValueError):
            pa.import_archive(arc, tempfile.mkdtemp())

    def test_reject_dotdot_traversal_member(self):
        arc = os.path.join(self.cwd, "dd.tar.gz")
        self._tar_with_raw_names(arc, ["HG01255/../../evil.txt"])
        with self.assertRaises(ValueError):
            pa.import_archive(arc, tempfile.mkdtemp())

    def test_reject_existing_destination(self):
        arc = os.path.join(self.cwd, "b.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            t.add(pa.individual_dir(self.cwd, "HG01255"), arcname="HG01255")
        # importing into the SAME cwd (already has HG01255) must refuse
        with self.assertRaises(ValueError):
            pa.import_archive(arc, self.cwd)

    # the modes a bundle's archive gives its members (by path in the bundle; "" is its
    # folder): the folder setgid, sticky and writable by all, a FASTA setgid, a file setuid
    MODES = {"": 0o3777, "paternal/genome/chr1.fa": 0o2775, "tool": 0o4755}

    def test_files_go_in_with_safe_modes_and_none_of_the_archive_s_owners(self):
        # tarfile gave each file the mode its member has, setuid and setgid too, and, run
        # as root (as the Docker image runs the web app, whose data folder is mounted from
        # the host), the owner it names: a bundle put a setuid file of any owner there
        src = pa.individual_dir(self.cwd, "HG01255")
        with open(os.path.join(src, "tool"), "w") as fh:
            fh.write("#!/bin/sh\n")

        def owned(info):  # (each member owned by another user, 1234)
            info.mode = self.MODES.get(info.name.partition("/")[2], info.mode)
            info.uid = info.gid = 1234
            info.uname = info.gname = ""
            return info

        arc = os.path.join(self.cwd, "modes.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            t.add(src, arcname="HG01255", filter=owned)
        umask = os.umask(0)  # (the 'data' filter makes a folder with the mode umask gives)
        os.umask(umask)
        for fallback in (False, True):  # (True: Python before 3.11.4, with no 'data' filter)
            with self.subTest(fallback=fallback):
                dest = tempfile.mkdtemp()
                self.addCleanup(shutil.rmtree, dest, ignore_errors=True)
                owners = []
                record = lambda path, uid, gid, *a, **kw: owners.append((path, uid, gid))  # noqa: E731
                with mock.patch("os.geteuid", return_value=0), \
                        mock.patch("os.chown", record), mock.patch("os.lchown", record), \
                        mock.patch.object(safe_tar, "_DATA_FILTER",
                                          None if fallback else safe_tar._DATA_FILTER), \
                        warnings.catch_warnings():
                    warnings.simplefilter("ignore", DeprecationWarning)
                    self.assertEqual(pa.import_archive(arc, dest), "HG01255")
                self.assertTrue(pa.resolve(dest, "HG01255")["complete"])
                root = pa.individual_dir(dest, "HG01255")
                modes = {os.path.relpath(path, root): stat.S_IMODE(os.stat(path).st_mode)
                         for path in [root] + [os.path.join(top, name)
                                               for top, dirs, files in os.walk(root)
                                               for name in dirs + files]}
                self.assertEqual((modes.pop("."), modes.pop("paternal/genome/chr1.fa"),
                                  modes.pop("tool")),
                                 (0o755 if fallback else 0o777 & ~umask, 0o755, 0o755))
                self.assertEqual({path: oct(m) for path, m in modes.items() if m & 0o7022}, {})
                self.assertEqual([o for o in owners if o[1:] != (-1, -1)], [])

    def test_a_zip_s_files_go_in_with_no_mode_of_the_archive_s(self):
        # zipfile sets no mode bits (nor an owner, which a zip does not hold): as before
        src = pa.individual_dir(self.cwd, "HG01255")
        arc = os.path.join(self.cwd, "modes.zip")
        with zipfile.ZipFile(arc, "w") as z:
            for top, _dirs, files in os.walk(src):
                for name in files:
                    full = os.path.join(top, name)
                    info = zipfile.ZipInfo(os.path.join("HG01255", os.path.relpath(full, src)))
                    info.external_attr = (stat.S_IFREG | 0o4755) << 16
                    with open(full, "rb") as fh:
                        z.writestr(info, fh.read())
        dest = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, dest, ignore_errors=True)
        self.assertEqual(pa.import_archive(arc, dest), "HG01255")
        self.assertTrue(pa.resolve(dest, "HG01255")["complete"])
        root = pa.individual_dir(dest, "HG01255")
        self.assertEqual([os.path.join(top, name) for top, _dirs, files in os.walk(root)
                          for name in files
                          if os.stat(os.path.join(top, name)).st_mode & 0o7022], [])


class TestBundlePathsStayInsideIt(unittest.TestCase):
    """metadata.json is the bundle's own content, a client's for an archive imported in
    Settings, and its paths were joined to the bundle folder as they were: an absolute
    path or '../' named a genome folder, chain and chromAlias file anywhere on the server.
    Only their existence was checked, so the bundle was imported, the individual listed as
    complete, and its search read those folders and files (its results held their
    sequences). A bundle naming anything outside its folder is refused, and one already
    installed (imported before this check) is not resolved: it is incomplete, so neither
    the search form nor assembly-search --assembly-individual takes it."""

    def setUp(self):
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        self.cwd = os.path.join(self.top, "data")
        os.makedirs(self.cwd)
        # someone else's files beside the data folder, and another individual's bundle
        self.other = os.path.join(self.top, "other_user")
        for name in ("genome/chrS.fa", "genome2/chrS.fa", "private.chain",
                     "private.chromAlias.txt"):
            path = os.path.join(self.other, name)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write(">chrS\nACGT\n" if name.endswith(".fa") else "SECRET\n")
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG01255", "maternal", "CM2")

    def metadata(self, **paternal):
        """The metadata.json of bundle EVIL: its paternal entry with ``paternal``'s paths,
        the rest as a valid bundle gives them."""
        haps = {hap: {"genome_dir": f"{hap}/genome", "chain": f"{hap}/x.chain",
                      "chromalias": f"{hap}/x.chromAlias.txt"} for hap in pa.HAPLOTYPES}
        haps["paternal"].update(paternal)
        return {"format": pa.BUNDLE_FORMAT, "format_version": 1, "individual": "EVIL",
                "haplotypes": haps}

    def bundle(self, folder, meta):
        """``folder`` laid out as a valid bundle, with ``meta`` as its metadata.json."""
        for hap in pa.HAPLOTYPES:
            os.makedirs(os.path.join(folder, hap, "genome"), exist_ok=True)
            for name, text in (("genome/chr1.fa", ">chr1\nACGT\n"), ("x.chain", "chain\n"),
                               ("x.chromAlias.txt", "x\n")):
                with open(os.path.join(folder, hap, name), "w") as fh:
                    fh.write(text)
        with open(os.path.join(folder, pa.METADATA_NAME), "w") as fh:
            json.dump(meta, fh)

    def outside(self):
        """Paternal entries naming a file or folder outside the bundle, each one there
        once the bundle is Assemblies/EVIL/."""
        other = os.path.join("..", "..", "..", "other_user")  # from Assemblies/EVIL/
        return [
            {"genome_dir": os.path.join(self.other, "genome")},
            {"chain": os.path.join(self.other, "private.chain")},
            {"chromalias": os.path.join(self.other, "private.chromAlias.txt")},
            {"genome_dir": os.path.join(other, "genome2")},
            {"chain": os.path.join("paternal", "..", other, "private.chain")},
            {"genome_dir": os.path.join("..", "HG01255", "paternal", "genome")},
            {"chain": ["paternal/x.chain"]},  # not a path at all
        ]

    def test_an_archive_naming_paths_outside_its_folder_is_refused(self):
        src = os.path.join(self.top, "src")
        for paternal in self.outside():
            with self.subTest(**{k: str(v) for k, v in paternal.items()}):
                self.bundle(os.path.join(src, "EVIL"), self.metadata(**paternal))
                arc = os.path.join(self.top, "evil.zip")
                with zipfile.ZipFile(arc, "w") as z:
                    for root, _dirs, files in os.walk(src):
                        for f in files:
                            full = os.path.join(root, f)
                            z.write(full, os.path.relpath(full, src))
                with self.assertRaisesRegex(ValueError, "is not a path inside the bundle"):
                    pa.import_archive(arc, self.cwd)
                self.assertEqual(os.listdir(os.path.join(self.cwd, pa.ASSEMBLIES_DIR)),
                                 ["HG01255"])
                shutil.rmtree(src)
        # one naming its own files only is imported as before
        self.bundle(os.path.join(src, "EVIL"), self.metadata())
        arc = os.path.join(self.top, "good.tar.gz")
        with tarfile.open(arc, "w:gz") as t:
            t.add(os.path.join(src, "EVIL"), arcname="EVIL")
        self.assertEqual(pa.import_archive(arc, self.cwd), "EVIL")
        self.assertTrue(pa.resolve(self.cwd, "EVIL")["complete"])

    def test_an_installed_bundle_naming_paths_outside_its_folder_is_not_resolved(self):
        folder = pa.individual_dir(self.cwd, "EVIL")
        for paternal in self.outside():
            with self.subTest(**{k: str(v) for k, v in paternal.items()}):
                self.bundle(folder, self.metadata(**paternal))
                entry = pa.resolve(self.cwd, "EVIL")
                self.assertIsNone(entry["paternal"])
                self.assertIsNotNone(entry["maternal"])
                self.assertFalse(entry["complete"])
                self.assertIn("EVIL", [b["individual"] for b in pa.discover_bundles(self.cwd)])
                self.assertEqual(pa.delete_target(self.cwd, "EVIL"), folder)  # removable
        # the individual whose files one of them named is as it was
        self.assertTrue(pa.resolve(self.cwd, "HG01255")["complete"])
        self.bundle(folder, self.metadata())
        self.assertTrue(pa.resolve(self.cwd, "EVIL")["complete"])


class TestBundleFastaNamesAShellReads(unittest.TestCase):
    """A search takes the name of each FASTA file of a genome folder, less its extension,
    for a chromosome's, and its post-analysis took that name, unquoted, into shell
    commands. A bundle brings its genome folders with the names its archive gives, a
    client's for one imported in Settings, and only some .fa had to be there: a FASTA
    named with characters a shell reads, of a sequence named chr1 (the name the search's
    own check of the folder reads), was imported, the individual listed as complete, and
    its search ran them (pool_post_analisi_snp.py now refuses such a name too:
    test_web_shell_safety). A bundle holding one is refused, with nothing installed; one
    of the usual names goes in as before."""

    # a FASTA file's name, less .fa, holding what a shell reads
    CRAFTED = ("chr1;touch MARK;x", "chr1$(touch MARK)", "chr1`touch MARK`", "chr1|touch MARK",
               "chr1&&touch MARK", "chr1>MARK", "chr1<x", "chr1'x'", 'chr1"x"', "chr1#x",
               "chr1 x", "chr1\tx", "chr1\ntouch MARK")

    def setUp(self):
        self.cwd = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)

    def archive(self, names, kind):
        """Bundle IND1 as a ``kind`` (zip or tar.gz) archive, its paternal genome folder
        holding an entry of each name of ``names`` (a folder for one that ends in '/', else
        a FASTA of sequence chr1), its maternal one chr1.fa."""
        meta = {"format": pa.BUNDLE_FORMAT, "format_version": 1, "individual": "IND1",
                "haplotypes": {hap: {"genome_dir": f"{hap}/genome", "chain": f"{hap}/x.chain",
                                     "chromalias": f"{hap}/x.chromAlias.txt"}
                               for hap in pa.HAPLOTYPES}}
        entries = [("IND1/metadata.json", json.dumps(meta).encode())]
        for hap, fastas in (("paternal", names), ("maternal", ["chr1.fa"])):
            entries += [(f"IND1/{hap}/genome/{name}", b">chr1\nACGT\n") for name in fastas]
            entries += [(f"IND1/{hap}/x.chain", b"chain\n"), (f"IND1/{hap}/x.chromAlias.txt", b"x\n")]
        path = os.path.join(self.cwd, f"IND1.{kind}")
        if kind == "zip":
            with zipfile.ZipFile(path, "w") as z:
                for name, data in entries:
                    z.writestr(name, b"" if name.endswith("/") else data)
        else:
            with tarfile.open(path, "w:gz") as t:
                for name, data in entries:
                    info = tarfile.TarInfo(name.rstrip("/"))
                    if name.endswith("/"):
                        info.type, info.mode = tarfile.DIRTYPE, 0o755
                        t.addfile(info)
                    else:
                        info.size = len(data)
                        t.addfile(info, io.BytesIO(data))
        return path

    def test_a_bundle_holding_a_fasta_named_so_is_refused(self):
        for kind in ("zip", "tar.gz"):
            for chrom in self.CRAFTED:
                for name in (f"{chrom}.fa", f"{chrom}.fa/"):  # (the search lists a folder too)
                    with self.subTest(kind=kind, name=name):
                        shutil.rmtree(os.path.join(self.cwd, pa.ASSEMBLIES_DIR), ignore_errors=True)
                        arc = self.archive(["chr1.fa", name], kind)
                        with self.assertRaisesRegex(
                                ValueError, r"^invalid personal-assembly bundle: paternal genome "
                                r"dir 'paternal/genome' holds .*: a FASTA file's name must not "
                                r"hold a blank or any of "):
                            pa.import_archive(arc, self.cwd)
                        self.assertEqual(os.listdir(os.path.join(self.cwd, pa.ASSEMBLIES_DIR)),
                                         [])

    def test_a_bundle_of_the_usual_names_goes_in_as_before(self):
        # what a sequence's name may hold ('*' and ':' of an HLA allele's too), and the file
        # of attributes macOS tar writes beside a FASTA (._chr1.fa), which is not refused
        names = ["chr1.fa", "chrUn_KI270302v1.fa", "NC_000913.3.fa", "HLA-A*01:01.fa",
                 "._chr1.fa"]
        for kind in ("zip", "tar.gz"):
            with self.subTest(kind=kind):
                self.assertEqual(pa.import_archive(self.archive(names, kind), self.cwd), "IND1")
                self.assertTrue(pa.resolve(self.cwd, "IND1")["complete"])
                genome = os.path.join(pa.individual_dir(self.cwd, "IND1"), "paternal", "genome")
                self.assertEqual(sorted(os.listdir(genome)), sorted(names))
                shutil.rmtree(pa.individual_dir(self.cwd, "IND1"))


class TestABundleIsNamedAfterItsFolder(unittest.TestCase):
    """An individual is listed (discover_bundles, resolve) under the name of its folder of
    Assemblies/, as assembly-search --assembly-individual and Settings' Remove take it.
    metadata.json's "individual", the bundle's own content, a client's for an archive
    imported in Settings, was the name listed: a bundle naming another individual took that
    one's place in the Individual picker and the Remove card, and one naming an individual
    with no folder was taken by Remove for a pairing of the legacy layout, removing paths
    named after the bundle's assembly names, outside the data folder too
    (test_settings_progress)."""

    def setUp(self):
        self.top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.top, ignore_errors=True)
        self.cwd = os.path.join(self.top, "data")
        _make_haplotype(self.cwd, "HG01255", "paternal", "CM1")
        _make_haplotype(self.cwd, "HG01255", "maternal", "CM2")

    def test_an_imported_bundle_is_listed_under_its_folder_s_name(self):
        # archives of HG01255's bundle under other top folders (as the import's own error
        # says to rename one), their metadata naming HG01255, or GHOST, who has no folder
        src = pa.individual_dir(self.cwd, "HG01255")
        for folder, individual in (("EVIL", "HG01255"), ("EVIL2", "GHOST")):
            meta = pa.read_metadata(self.cwd, "HG01255")
            meta["individual"] = individual
            for hap in pa.HAPLOTYPES:
                meta["haplotypes"][hap]["assembly_name"] = "./../outside/victim"
            arc = os.path.join(self.top, f"{folder}.tar.gz")
            with tarfile.open(arc, "w:gz") as t:
                t.add(src, arcname=folder, filter=lambda info: None
                      if info.name.endswith(pa.METADATA_NAME) else info)
                data = json.dumps(meta).encode()
                info = tarfile.TarInfo(f"{folder}/{pa.METADATA_NAME}")
                info.size = len(data)
                t.addfile(info, io.BytesIO(data))
            self.assertEqual(pa.import_archive(arc, self.cwd), folder)
        found = {b["individual"]: b for b in pa.discover_bundles(self.cwd)}
        self.assertEqual(sorted(found), ["EVIL", "EVIL2", "HG01255"])
        for name in found:
            self.assertEqual(found[name]["root"], os.path.join(pa.ASSEMBLIES_DIR, name))
            self.assertTrue(found[name]["complete"])
            self.assertEqual(pa.resolve(self.cwd, name)["individual"], name)
            self.assertEqual(pa.delete_target(self.cwd, name), pa.individual_dir(self.cwd, name))
        # the downloader's bundle, its folder named after its individual: as before
        self.assertEqual(pa.read_metadata(self.cwd, "HG01255")["individual"], "HG01255")
        self.assertEqual(found["HG01255"]["paternal"]["assembly_name"], "HG01255_paternal_CM1")


if __name__ == "__main__":
    unittest.main()

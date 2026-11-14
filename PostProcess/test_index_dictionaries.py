#!/usr/bin/env python3
"""A variant index and the Dictionaries of its dataset come from one build, or no search reads
them together (PostProcess/index_dictionaries.py; PP-14 of CRISPRme 2.10.0's release).

A download of an index of 2.10.0 installs its panel's Dictionaries in place of those an archive
of 2.8.1 installed, and a data folder upgraded from 2.8.1 then holds the old index beside the new
Dictionaries. The download said so on stdout; complete-search, and the search the web launches,
ran it on, its SNP rows read with the other build's carriers, until its indel search refused the
INDELS index chromosome by chromosome, hours later. Each now refuses it before anything is
written, in one line that says what to do: the index of 2.8.1 (index format 1, its enrichment's
fake-indel windows 25/27) beside Dictionaries of 2.10.0 (30/30), and the other way. Each
release's own pair runs, and so does an index of 2.10.0 built from 2.8.1's enriched-genome kit
beside that kit's Dictionaries.

TheRefusal and TheCommandLine read fixtures of each case. ASearchRefusesItFirst runs the real
complete-search and the search script as the web launches it, with the stub engine of
test_build_index_only_atomic, on a dataset build-index-only built there.
"""

import contextlib
import io
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import index_dictionaries  # noqa: E402
import index_fit  # noqa: E402
import test_build_index_only_atomic as atomic  # noqa: E402  (its stub engine and fixture)
from test_indel_window import paste_both  # noqa: E402  (a printed command, in bash and zsh)

DS = "hg38_HGDP1kGP"
HEADER = "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq"


def _write(path, text, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode) as fh:
        fh.write(text)


def log_rows(chrom, left, right, n=3):
    """Rows of an indel log whose windows hold ``left`` bases before each indel and ``right``
    after it (one-base REF, an insertion), as add-variants writes them."""
    rows = []
    for i in range(n):
        pos = 1000 + 100 * i
        start, end = pos - 1 - left, pos + right
        rows.append(f"{chrom}_{start}-{end}_{i + 1}\tS1\trs{i}\t0.1\t{chrom}_{pos}_A_AT\t0,1\tACGT")
    return rows


def write_logs(workdir, window, chroms=("chr22",), dataset=DS):
    """Dictionaries/log_indels_<dataset>/log<chrom>.txt of rows of ``window`` (header-only for
    None)."""
    for chrom in chroms:
        rows = log_rows(chrom, *window) if window else []
        _write(os.path.join(workdir, "Dictionaries", f"log_indels_{dataset}", f"log{chrom}.txt"),
               "\n".join([HEADER] + rows) + "\n")


def _tree_bytes(fmt):
    """The first bytes of a .bin of index format ``fmt``: a format-1 file opens with its leaf
    count and the guide length; a format-2 one with the magic and the version."""
    if fmt == 1:
        return struct.pack("=ii", 5, 20) + b"\0" * 32
    return b"CZI\xfe" + struct.pack("=6i", 2, 31, 0, 20, 3, 3) + b"NRG" + b"\0" * 32


def make_index(library, name, fmt, sub_windows=None, sub_format=None, indels=True,
               chroms=("chr22",)):
    """A variant index ``name`` of index format ``fmt`` in ``library``, with its INDELS index:
    a sub-index a chromosome, of format ``sub_format`` (``fmt`` by default), recording
    ``sub_windows`` (its .index_meta.json's indel_window; none for None). Returns its path."""
    path = os.path.join(library, name)
    pam = index_fit.parse_index_name(name)[0]
    for chrom in chroms:
        _write(os.path.join(path, f"{pam}_{chrom}_1.bin"), _tree_bytes(fmt), "wb")
    if indels:
        n = index_fit.parse_index_name(name)[1]
        for chrom in chroms:
            sub = os.path.join(path + "_INDELS", f"{pam}_{n}_fake{chrom}")
            _write(os.path.join(sub, f"{pam}_fake{chrom}_1.bin"),
                   _tree_bytes(sub_format or fmt), "wb")
            if sub_windows:
                _write(os.path.join(sub, index_fit.INDEX_META),
                       json.dumps({"index_format": sub_format or fmt,
                                   "indel_window": list(sub_windows)}))
    return path


@contextlib.contextmanager
def engine(index_format):
    """A crispritz.py on PATH that states ``index_format`` (None: one that states none, as
    CRISPRitz 2.9.1's)."""
    with tempfile.TemporaryDirectory() as bindir:
        text = "#!/bin/sh\nVERSION = '%s'\n" % ("2.10.0" if index_format else "2.9.1")
        if index_format:
            text += f"INDEX_FORMAT = {index_format}\nINDEL_LEFT = 30\nINDEL_RIGHT = 30\n"
        _write(os.path.join(bindir, "crispritz.py"), text)
        os.chmod(os.path.join(bindir, "crispritz.py"), 0o755)
        with mock.patch.dict(os.environ, PATH=bindir + os.pathsep + os.environ.get("PATH", "")):
            yield


class _Folder(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="idxdicts_")
        self.addCleanup(shutil.rmtree, self.work, True)
        self.lib = os.path.join(self.work, "genome_library")
        os.makedirs(self.lib)


class TheRefusal(_Folder):
    """index_dictionaries.refusal: the windows of the enrichment an index was built from
    against those its dataset's indel logs record."""

    def test_an_index_of_2_8_1_beside_the_dictionaries_of_2_10_0_is_refused(self):
        # (a published archive of 2.8.1: index format 1, its INDELS index too, no record)
        idx = make_index(self.lib, f"NRG_3_hg38+{DS}", 1)
        write_logs(self.work, (30, 30))
        with engine(2):
            why = index_dictionaries.refusal(idx, self.work, DS)
        self.assertIsNotNone(why)
        self.assertEqual(len(why.splitlines()), 1)
        self.assertIn(f"NRG_3_hg38+{DS} is an index of a build of {DS} whose enrichment wrote "
                      "fake-indel windows of 25/27 (CRISPRme 2.8.1 or earlier)", why)
        self.assertIn(f"the Dictionaries of {DS} in {os.path.join(self.work, 'Dictionaries')} are "
                      "of another build, whose indel logs record windows of 30/30 (CRISPRme "
                      "2.10.0)", why)
        # what to do: its replacement of 2.10.0 (the index set's), and this one removed
        self.assertIn(f"Search 20-NRG_3_hg38+{DS}, the index of CRISPRme 2.10.0 that replaces it "
                      f"(crisprme.py download --what index --index-name 20-NRG_3_hg38+{DS})", why)
        self.assertIn(f"rm -rf {idx} {idx}_INDELS", why)
        self.assertIn("docs/TROUBLESHOOTING.md section 16", why)

    def test_the_other_way_an_index_of_2_10_0_beside_the_dictionaries_of_2_8_1(self):
        idx = make_index(self.lib, f"20-NRR_3_hg38+{DS}", 2, sub_windows=(30, 30))
        write_logs(self.work, (25, 27))
        with engine(2):
            why = index_dictionaries.refusal(idx, self.work, DS)
        self.assertIn("whose enrichment wrote fake-indel windows of 30/30 (CRISPRme 2.10.0)", why)
        self.assertIn("whose indel logs record windows of 25/27 (CRISPRme 2.8.1 or earlier)", why)
        self.assertIn("Put back the Dictionaries it was built with", why)
        self.assertEqual(len(why.splitlines()), 1)

    def test_each_release_with_its_own_dictionaries_runs(self):
        cases = (
            # the upgrade path of an index of 2.8.1: 2.10.0 reads format 1, with its own logs
            ("NRG_3_hg38+" + DS, 1, None, (25, 27)),
            # 2.10.0's index with its own
            ("20-NRR_3_hg38+" + DS, 2, (30, 30), (30, 30)),
            # an index of 2.10.0 built from 2.8.1's enriched-genome kit, beside that kit's
            ("20-NRR_3_hg38+" + DS, 2, (25, 27), (25, 27)),
        )
        for name, fmt, subs, logged in cases:
            with self.subTest(name=name, logged=logged):
                shutil.rmtree(self.lib)
                shutil.rmtree(os.path.join(self.work, "Dictionaries"), ignore_errors=True)
                idx = make_index(self.lib, name, fmt, sub_windows=subs)
                write_logs(self.work, logged)
                with engine(2):
                    self.assertIsNone(index_dictionaries.refusal(idx, self.work, DS))

    def test_nothing_tells_the_builds_apart_nothing_is_refused(self):
        idx = make_index(self.lib, "20-NRR_3_hg38+" + DS, 2, sub_windows=(30, 30))
        write_logs(self.work, None)  # a dataset without indels: logs of a header alone
        with engine(2):
            self.assertIsNone(index_dictionaries.refusal(idx, self.work, DS))
        # an index of format 2 without its INDELS index, which records no window
        shutil.rmtree(self.lib)
        idx = make_index(self.lib, "20-NRR_3_hg38+" + DS, 2, indels=False)
        write_logs(self.work, (25, 27))
        with engine(2):
            self.assertIsNone(index_dictionaries.refusal(idx, self.work, DS))
        self.assertIsNone(index_dictionaries.refusal("", self.work, DS))
        self.assertIsNone(index_dictionaries.refusal(os.path.join(self.lib, "none"), self.work, DS))

    def test_an_index_without_its_indels_index_tells_nothing_here(self):
        """Its windows are those of its INDELS index, built with it: without one, nothing here
        tells them (the search script refuses such an index, or builds its INDELS index from
        the enrichment in place, whose fake-indel genome it reads against the logs first)."""
        idx = make_index(self.lib, f"NRG_3_hg38+{DS}", 1, indels=False)
        write_logs(self.work, (30, 30))
        with engine(2):
            self.assertIsNone(index_dictionaries.refusal(idx, self.work, DS))

    def test_each_chromosome_is_read_against_its_own_log(self):
        """As the indel search reads them (pool_search_indels' _window_mismatch): a sub-index of
        chr1 built from 25/27 beside a log of 30/30 is refused, though chr2's agree."""
        idx = make_index(self.lib, "20-NRR_3_hg38+" + DS, 2, sub_windows=(30, 30),
                         chroms=("chr1", "chr2"))
        meta = os.path.join(idx + "_INDELS", "NRR_3_fakechr1", index_fit.INDEX_META)
        _write(meta, json.dumps({"index_format": 2, "indel_window": [25, 27]}))
        write_logs(self.work, (30, 30), chroms=("chr1", "chr2"))
        with engine(2):
            self.assertIsNotNone(index_dictionaries.refusal(idx, self.work, DS))

    def test_its_rm_command_removes_that_index_pasted_into_bash_or_zsh(self):
        # (a folder whose name holds a space and an apostrophe: the command is quoted)
        work = os.path.join(self.work, "Luca's data")
        lib = os.path.join(work, "genome_library")
        idx = make_index(lib, f"NRG_3_hg38+{DS}", 1)
        make_index(lib, f"20-NRG_3_hg38+{DS}", 2, sub_windows=(30, 30))
        write_logs(work, (30, 30))
        with engine(2):
            why = index_dictionaries.refusal(idx, work, DS)
        command = why.split("and remove this one (", 1)[1].split("); or keep it", 1)[0]
        self.assertTrue(command.startswith("rm -rf "), command)
        paste_both(self, command, self.work)
        self.assertEqual(sorted(os.listdir(lib)),
                         [f"20-NRG_3_hg38+{DS}", f"20-NRG_3_hg38+{DS}_INDELS"])


class TheCommandLine(_Folder):
    """``index_dictionaries.py <index> <workdir> <dataset>``, as the search script runs it: the
    refusal on stdout, exit 3; nothing, exit 0; an empty index (none installed), exit 0."""

    def run_it(self, *args):
        return subprocess.run([sys.executable, os.path.join(HERE, "index_dictionaries.py"),
                               *args], capture_output=True, text=True)

    def test_its_status_says_it(self):
        old = make_index(self.lib, f"NRG_3_hg38+{DS}", 1)
        new = make_index(self.lib, f"20-NRR_3_hg38+{DS}", 2, sub_windows=(30, 30))
        write_logs(self.work, (30, 30))
        with engine(2):
            refused = self.run_it(old, self.work, DS)
            fine = self.run_it(new, self.work, DS)
            none = self.run_it("", self.work, DS)
        self.assertEqual(refused.returncode, 3, refused.stderr)
        self.assertEqual(len(refused.stdout.splitlines()), 1)
        self.assertIn("are of another build", refused.stdout)
        self.assertEqual((refused.stderr, fine.returncode, fine.stdout, fine.stderr),
                         ("", 0, "", ""))
        self.assertEqual((none.returncode, none.stdout, none.stderr), (0, "", ""))
        self.assertEqual(self.run_it("a", "b").returncode, 2)


class CompleteSearchRefusesItFirst(_Folder):
    """crisprme.py's complete-search pre-flight (_refuse_mixed_builds) takes, for each dataset,
    the index the search script would take (_index_the_search_takes), and stops on one of
    another build than the dataset's Dictionaries, before the run writes anything."""

    @classmethod
    def setUpClass(cls):
        import test_bulge_cap
        cls.crisprme = test_bulge_cap._load_crisprme()

    def check(self, motif="NRG", vcf=DS):
        crisprme = self.crisprme
        err = io.StringIO()
        with mock.patch.object(crisprme, "current_working_directory", self.work + "/"):
            try:
                with contextlib.redirect_stderr(err):
                    crisprme._refuse_mixed_builds([os.path.join(self.work, "VCFs", vcf)], "hg38",
                                                  motif, 1, 1, False, True, 20, "_")
            except SystemExit as stop:
                self.assertEqual(stop.code, 1)
                return err.getvalue()
        self.assertEqual(err.getvalue(), "")
        return None

    def test_the_index_the_search_takes_is_the_one_read(self):
        # an upgraded data folder: 2.8.1's NRG index, 2.10.0's NRR index and its Dictionaries
        make_index(self.lib, f"NRG_3_hg38+{DS}", 1)
        make_index(self.lib, f"20-NRR_3_hg38+{DS}", 2, sub_windows=(30, 30))
        write_logs(self.work, (30, 30))
        with engine(2):
            # NRG: the search takes NRG_3 (its own PAM, ahead of NRR pruned): refused
            said = self.check("NRG")
            # NRR: the search takes 20-NRR_3, of 2.10.0 like the Dictionaries: it runs
            self.assertIsNone(self.check("NRR"))
        self.assertTrue(said.startswith(f"Error: NRG_3_hg38+{DS} is an index of a build of"), said)
        self.assertEqual(len(said.splitlines()), 1)

    def test_a_folder_of_2_8_1_runs_as_before(self):
        make_index(self.lib, f"NRG_3_hg38+{DS}", 1)
        write_logs(self.work, (25, 27))
        with engine(2):
            self.assertIsNone(self.check("NRG"))

    def test_it_comes_before_the_run_writes_anything(self):
        """complete-search calls it beside _refuse_no_variant_index, before it writes the
        output folder's first file (.command_line.txt)."""
        with open(os.path.join(REPO, "crisprme.py")) as fh:
            body = fh.read().split("\ndef complete_search() -> None:", 1)[1]
        body = body.split("\ndef ", 1)[0]
        self.assertIn("_refuse_mixed_builds(", body)
        self.assertLess(body.index("_refuse_mixed_builds("),
                        body.index('outputfolder + "/.command_line.txt"'))
        self.assertLess(body.index("_refuse_no_variant_index("),
                        body.index("_refuse_mixed_builds("))


def setUpModule():
    atomic.setUpModule()


def tearDownModule():
    atomic.tearDownModule()


@unittest.skipUnless(sys.platform.startswith("linux") and shutil.which("bash"),
                     "the stub engine's harness runs where the search script runs: Linux")
class ASearchRefusesItFirst(atomic._Fixture):
    """The real complete-search, and the search script as the web launches it, on a dataset
    build-index-only built here with the stub engine (whose fake-indel windows are 25/34), then
    on the same with its Dictionaries' indel logs another build's (30/30): the search goes past
    its checks of the index on the first, and is refused on the second, before anything is
    written or searched."""

    def other_builds_logs(self):
        """The dataset's indel logs replaced by another build's, of windows of 30/30."""
        logs = self.path("Dictionaries", "log_indels_ds")
        for name in os.listdir(logs):
            os.remove(os.path.join(logs, name))
        for chrom, seq in atomic.GENOME.items():
            _write(os.path.join(logs, f"log{chrom}.txt"),
                   f"{HEADER}\n{chrom}_9-70_1\tS2\trs_{chrom}_40\t0.25\t{chrom}_40_{seq[39]}_"
                   f"{seq[39]}T\t0,62\t{seq[9:70]}\n")

    def searches(self):
        return [c for c in self.calls() if c.split()[:1] == ["search"]]

    def test_its_own_dictionaries_are_searched(self):
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assert_searched()

    def test_another_builds_are_refused_before_anything_is_written(self):
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.other_builds_logs()
        os.remove(self.log)
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {atomic.VAR} is an index of a build of ds whose enrichment wrote "
                      f"fake-indel windows of 25/34, and the Dictionaries of ds in "
                      f"{self.path('Dictionaries')} are of another build, whose indel logs "
                      f"record windows of 30/30 (CRISPRme 2.10.0)", err)
        self.assertEqual(len([ln for ln in err.splitlines() if ln.startswith("Error: ")]), 1)
        self.assertNotIn("Launching job ", out)
        self.assertEqual(os.listdir(self.path("Results", "search")), [])
        self.assertEqual(self.searches(), [])

    def test_the_search_the_web_launches_is_refused_before_it_searches(self):
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.other_builds_logs()
        os.remove(self.log)
        rc, out, err = self.run_web_search()
        self.assertEqual(rc, 1, out + err)
        self.assertTrue(err.startswith(f"ERROR: {atomic.VAR} is an index of a build of ds whose "
                                       f"enrichment wrote fake-indel windows of 25/34"), err)
        self.assertEqual(len(err.splitlines()), 1, err)
        self.assertEqual(self.searches(), [])


if __name__ == "__main__":
    unittest.main()

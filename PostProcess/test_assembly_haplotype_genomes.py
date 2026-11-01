#!/usr/bin/env python3
"""Each assembly-search haplotype is searched on its own sequence.

complete-search names a genome -- and the index it builds or reuses,
genome_library/<PAM>_<N>_<name> -- by the genome folder's name. A personal-assembly
bundle keeps every haplotype in Assemblies/<individual>/<hap>/genome, so both
haplotypes of every individual were searched as "genome": with bulges, the paternal
search built NRG_2_genome, and the maternal search -- and every later individual's --
reused it, searched the paternal sequence, and the reconciliation then reported every
site as found on both haplotypes. A 0-bulge search scanned Genomes/genome/, which does
not exist, and came back empty.

These tests pin:
  * the name a haplotype genome is searched under (haplotype_search_name) and the
    folder it is searched with (haplotype_genome_view): unique per haplotype and per
    individual, for --assembly-individual and the explicit-folder form alike, stable
    across runs, and left out for a legacy Genomes/<name> folder; never the name of
    the legacy layout's Genomes/<individual>_<haplotype> folder, so neither one
    takes the other's index; made without a hitch by runs that start together;
  * assembly_search itself (crisprme.py run in-process; each haplotype's complete-search
    replaced by a stand-in that names, builds and reuses indexes by genome name as the
    pipeline does, and records which sequence it searched): each haplotype is searched
    on its own sequence, and a result written by the colliding search is not reused;
  * the record of which genome an index holds (.genome_source): an index of the same
    name built from another folder, or before the folder's FASTA files changed (names,
    sizes and modification times are compared, not the sequence), is removed and built
    again, one that records nothing (built by an earlier CRISPRme, or downloaded) is kept;
  * the pipeline keeps a linked genome's own name, picks only an index named for the
    genome itself (not one whose genome name merely ends in it), and its 0-bulge scan
    reads the genome folder it is given.
"""
import contextlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import types
import unittest
from unittest import mock

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CRISPRME_PY = os.path.join(os.path.dirname(HERE), "crisprme.py")
SHELL_SCRIPT = os.path.join(HERE, "submit_job_automated_new_multiple_vcfs.sh")
PAMS_DIR = os.path.join(os.path.dirname(HERE), "PAMs")
sys.path.insert(0, HERE)

import assembly_reconcile as ar  # noqa: E402
import personal_assembly as pa  # noqa: E402


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_annotation_tmpdirs._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    with open(CRISPRME_PY) as fh:
        src = fh.read()
    cut = src.index("\nif len(sys.argv) < 2:")  # KeyError-loud if the dispatch shape ever changes
    mod = types.ModuleType("crisprme_under_test_haplotype_genomes")
    mod.__file__ = CRISPRME_PY
    code = compile(src[:cut], CRISPRME_PY, "exec")
    # Stubs exist ONLY so crisprme.py's top-level imports resolve while we exec its
    # definitions; they are removed again in the finally, so they cannot shadow the real
    # modules for a test collected later.
    _inserted = []
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        if "Bio" not in sys.modules:
            sys.modules["Bio"] = _bio
            _inserted.append("Bio")
        sys.modules["Bio.Seq"] = _bio_seq
        _inserted.append("Bio.Seq")
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISPR_BULGE_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = mock.MagicMock()
            _inserted.append(_m)
    try:
        exec(code, mod.__dict__)
    finally:
        for _m in _inserted:
            sys.modules.pop(_m, None)
    return mod


crisprme = _load_crisprme()


# --------------------------------------------------------------------------- #
# Synthetic bundles
# --------------------------------------------------------------------------- #
def _write_genome(gdir, seq, mtime=None):
    os.makedirs(gdir, exist_ok=True)
    path = os.path.join(gdir, "chr1.fa")
    with open(path, "w") as fh:
        fh.write(f">chr1\n{seq}\n")
    if mtime is not None:
        os.utime(path, (mtime, mtime))


def _make_bundle(cwd, individual, seqs, assembly_names=True):
    """Assemblies/<individual>/ with a paternal and a maternal haplotype, each in a
    folder named "genome", as the HPRC fetch and an imported bundle store them."""
    root = os.path.join(cwd, pa.ASSEMBLIES_DIR, individual)
    meta = {"format": pa.BUNDLE_FORMAT, "format_version": 1, "individual": individual,
            "haplotypes": {}}
    for hap in pa.HAPLOTYPES:
        _write_genome(os.path.join(root, hap, "genome"), seqs[hap])
        for kind, fname in (("chain", f"{individual}_{hap}_vs_GRCh38.chain"),
                            ("chromalias", f"{individual}_{hap}.chromAlias.txt")):
            with open(os.path.join(root, hap, fname), "w") as fh:
                fh.write(f"{kind}\n")
        entry = {"genome_dir": f"{hap}/genome",
                 "chain": f"{hap}/{individual}_{hap}_vs_GRCh38.chain",
                 "chromalias": f"{hap}/{individual}_{hap}.chromAlias.txt"}
        if assembly_names:
            entry["assembly_name"] = f"{individual}_{hap}_CM00000{len(meta['haplotypes'])}.1"
        meta["haplotypes"][hap] = entry
    with open(os.path.join(root, pa.METADATA_NAME), "w") as fh:
        json.dump(meta, fh)
    return root


class _Workdir(unittest.TestCase):
    def setUp(self):
        self.cwd = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.cwd, ignore_errors=True)

    def genome(self, individual, hap):
        return os.path.join(self.cwd, pa.ASSEMBLIES_DIR, individual, hap, "genome")


# --------------------------------------------------------------------------- #
# The name and the folder each haplotype is searched with
# --------------------------------------------------------------------------- #
class EachHaplotypeHasItsOwnSearchName(_Workdir):
    def setUp(self):
        super().setUp()
        _make_bundle(self.cwd, "W1", {"paternal": "ACGT", "maternal": "TTTT"})
        _make_bundle(self.cwd, "W2", {"paternal": "GGGG", "maternal": "CCCC"})

    def test_bundle_haplotype_found_from_the_genome_folder_alone(self):
        hap = pa.bundle_haplotype(self.genome("W1", "maternal"), "maternal")
        self.assertEqual(hap["folder"], "W1")
        self.assertEqual(hap["haplotype"], "maternal")
        self.assertEqual(hap["assembly_name"], "W1_maternal_CM000001.1")
        # identity follows the folder, not the flag it was passed with
        self.assertEqual(pa.bundle_haplotype(self.genome("W1", "maternal"), "paternal")["haplotype"],
                         "maternal")
        self.assertIsNone(pa.bundle_haplotype(self.cwd))

    def test_names_differ_between_haplotypes_and_individuals(self):
        names = {(i, h): ar.haplotype_search_name(self.cwd, self.genome(i, h), h)
                 for i in ("W1", "W2") for h in pa.HAPLOTYPES}
        self.assertEqual(len(set(names.values())), 4)
        for (ind, hap), name in names.items():
            self.assertRegex(name, rf"^{ind}_{hap}_[0-9a-f]{{8}}$")

    def test_no_name_is_that_of_the_legacy_layouts_folder(self):
        # before 2.6.0 the same haplotype was kept as Genomes/W1_paternal, indexed as
        # NRG_2_W1_paternal: a bundle searched as W1_paternal took that index, and
        # handed its own to that folder's searches
        for hap in pa.HAPLOTYPES:
            name = ar.haplotype_search_name(self.cwd, self.genome("W1", hap), hap)
            self.assertNotEqual(name, f"W1_{hap}")

    def test_folders_whose_names_read_alike_get_different_names(self):
        # both are written HG_002 in a file name
        for ind in ("HG.002", "HG_002"):
            _make_bundle(self.cwd, ind, {"paternal": "ACGT", "maternal": "TTTT"})
        a = ar.haplotype_search_name(self.cwd, self.genome("HG.002", "paternal"), "paternal")
        b = ar.haplotype_search_name(self.cwd, self.genome("HG_002", "paternal"), "paternal")
        self.assertNotEqual(a, b)
        self.assertRegex(a, r"^HG_002_paternal_[0-9a-f]{8}$")
        self.assertRegex(b, r"^HG_002_paternal_[0-9a-f]{8}$")

    def test_other_folders_of_the_same_name_get_different_names(self):
        a = os.path.join(self.cwd, "elsewhere", "a", "genome")
        b = os.path.join(self.cwd, "elsewhere", "b", "genome")
        _write_genome(a, "ACGT")
        _write_genome(b, "ACGT")
        na = ar.haplotype_search_name(self.cwd, a, "paternal")
        nb = ar.haplotype_search_name(self.cwd, b, "paternal")
        self.assertNotEqual(na, nb)
        self.assertTrue(na.startswith("genome_") and nb.startswith("genome_"))
        self.assertRegex(na, r"^[A-Za-z0-9_-]+$")

    def test_name_keeps_only_characters_safe_in_file_names(self):
        _make_bundle(self.cwd, "odd name+v2", {"paternal": "A", "maternal": "C"})
        self.assertRegex(
            ar.haplotype_search_name(self.cwd, self.genome("odd name+v2", "paternal"), "paternal"),
            r"^odd_name_v2_paternal_[0-9a-f]{8}$")

    def test_view_is_named_after_the_haplotype_and_reads_its_genome(self):
        for ind in ("W1", "W2"):
            for hap in pa.HAPLOTYPES:
                view = ar.haplotype_genome_view(self.cwd, self.genome(ind, hap), hap)
                self.assertRegex(os.path.basename(view), rf"^{ind}_{hap}_[0-9a-f]{{8}}$")
                self.assertEqual(os.path.realpath(view), os.path.realpath(self.genome(ind, hap)))
                self.assertTrue(view.startswith(os.path.join(self.cwd, ar.HAPLOTYPE_GENOMES_DIRNAME)))

    def test_view_is_the_same_path_on_every_run(self):
        first = ar.haplotype_genome_view(self.cwd, self.genome("W1", "paternal"), "paternal")
        again = ar.haplotype_genome_view(self.cwd + "/", self.genome("W1", "paternal") + "/", "paternal")
        self.assertEqual(first, again)

    def test_view_link_is_relative_inside_the_working_directory(self):
        view = ar.haplotype_genome_view(self.cwd, self.genome("W1", "paternal"), "paternal")
        self.assertFalse(os.path.isabs(os.readlink(view)))
        # the working directory mounted elsewhere (Docker's /DATA) still resolves it
        moved = self.cwd + "_moved"
        os.rename(self.cwd, moved)
        try:
            self.assertTrue(os.path.isfile(os.path.join(
                moved, os.path.relpath(view, self.cwd), "chr1.fa")))
        finally:
            os.rename(moved, self.cwd)

    def test_a_stale_view_is_pointed_back_at_its_folder(self):
        view = ar.haplotype_genome_view(self.cwd, self.genome("W1", "paternal"), "paternal")
        os.remove(view)
        os.symlink(self.genome("W2", "paternal"), view)
        self.assertEqual(ar.haplotype_genome_view(self.cwd, self.genome("W1", "paternal"), "paternal"), view)
        self.assertEqual(os.path.realpath(view), os.path.realpath(self.genome("W1", "paternal")))

    def test_a_legacy_genomes_folder_is_searched_as_it_is(self):
        legacy = os.path.join(self.cwd, "Genomes", "HG01255_paternal")
        _write_genome(legacy, "ACGT")
        self.assertEqual(ar.haplotype_genome_view(self.cwd, legacy, "paternal"), legacy)


# --------------------------------------------------------------------------- #
# Which genome an index holds
# --------------------------------------------------------------------------- #
class IndexRecordsTheGenomeItHolds(_Workdir):
    def setUp(self):
        super().setUp()
        _make_bundle(self.cwd, "W1", {"paternal": "ACGT", "maternal": "TTTT"})
        self.view = ar.haplotype_genome_view(self.cwd, self.genome("W1", "paternal"), "paternal")
        self.name = os.path.basename(self.view)
        self.lib = os.path.join(self.cwd, "genome_library")

    def index(self, name, stamp=None):
        d = os.path.join(self.lib, name)
        os.makedirs(d)
        with open(os.path.join(d, "NRG_chr1_1.bin"), "w") as fh:
            fh.write("bin")
        if stamp is not None:
            with open(os.path.join(d, ar.INDEX_SOURCE_FILENAME), "w") as fh:
                json.dump(stamp, fh)
        return d

    def clear(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            present = ar.clear_stale_haplotype_indexes(self.cwd, self.view)
        return present, out.getvalue()

    def test_an_index_built_from_other_content_is_removed(self):
        d = self.index(f"NRG_2_{self.name}",
                       {"genome": "Assemblies/W1/paternal/genome", "fingerprint": "0" * 40})
        present, said = self.clear()
        self.assertFalse(os.path.exists(d))
        self.assertEqual(present, [])
        self.assertIn("has changed since it was built", said)

    def test_an_index_that_records_nothing_is_kept(self):
        # built by a search of this haplotype that has not ended yet (a run of another
        # guide, say), or was stopped: removing it could pull it from under that search
        d = self.index(f"NRG_2_{self.name}")
        present, said = self.clear()
        self.assertTrue(os.path.isdir(d))
        self.assertEqual(present, [d])
        self.assertEqual(said, "")

    def test_an_index_of_this_genome_is_kept(self):
        d = self.index(f"NRG_2_{self.name}",
                       {"genome": os.path.join("Assemblies", "W1", "paternal", "genome"),
                        "fingerprint": ar.genome_fingerprint(self.view)})
        self.assertEqual(self.clear(), ([d], ""))

    def test_an_index_built_from_another_folder_is_removed(self):
        # the same file names, sizes and times -- the maternal genome of the same
        # bundle, say -- but another folder
        d = self.index(f"NRG_2_{self.name}",
                       {"genome": os.path.join("Genomes", "W1_paternal"),
                        "fingerprint": ar.genome_fingerprint(self.view)})
        present, said = self.clear()
        self.assertFalse(os.path.exists(d))
        self.assertIn("it was built from Genomes/W1_paternal, not Assemblies/W1/paternal/genome", said)

    def test_an_index_still_being_built_is_left_to_its_builder(self):
        d = self.index(f"NRG_2_{self.name}", {"genome": "x", "fingerprint": "0" * 40})
        os.makedirs(d + ".lock")
        present, _ = self.clear()
        self.assertTrue(os.path.isdir(d))
        self.assertEqual(present, [])

    def test_other_genomes_indexes_are_never_touched(self):
        # NRG_2_W1_paternal: the index of Genomes/W1_paternal, the layout before 2.6.0
        others = [self.index(n, {"genome": "x", "fingerprint": "0" * 40})
                  for n in ("NRG_2_hg38", f"NRG_2_X_{self.name}", f"NRG_2_{self.name}+hg38_1000G",
                            "NRG_2_W1_paternal")]
        self.clear()
        for d in others:
            self.assertTrue(os.path.isdir(d), d)

    def test_only_indexes_built_by_the_search_are_stamped(self):
        kept = self.index(f"NRG_3_{self.name}")
        present, _ = self.clear()
        built = self.index(f"NRG_2_{self.name}")
        ar.stamp_built_haplotype_indexes(self.cwd, self.view, present)
        self.assertIsNone(ar._read_index_source(kept))
        record = ar._read_index_source(built)
        self.assertEqual(record["genome"], os.path.join("Assemblies", "W1", "paternal", "genome"))
        self.assertEqual(record["fingerprint"], ar.genome_fingerprint(self.view))
        self.assertEqual(sorted(os.listdir(built)), [ar.INDEX_SOURCE_FILENAME, "NRG_chr1_1.bin"])

    def test_fingerprint_follows_the_fasta_files(self):
        before = ar.genome_fingerprint(self.view)
        _write_genome(self.genome("W1", "paternal"), "ACGA", mtime=time.time() - 3600)
        self.assertNotEqual(ar.genome_fingerprint(self.view), before)


class RunsStartingTogetherDoNotTripEachOtherUp(_Workdir):
    """Searches that start together on one data folder (several guides launched at
    once, each in its own `docker run`) all make the haplotype's link, and may all
    record the index they found built. Their temporary files were named after the
    PID, which is 1 for crisprme.py in every container: one run removed or renamed
    the other's and stopped with FileNotFoundError. Here the other run starts at the
    worst moment, right after this one made its temporary file, with os.getpid() 1
    in both."""

    def setUp(self):
        super().setUp()
        _make_bundle(self.cwd, "W1", {"paternal": "ACGT", "maternal": "TTTT"})
        self.genome_dir = self.genome("W1", "paternal")
        patcher = mock.patch.object(ar.os, "getpid", return_value=1)
        patcher.start()
        self.addCleanup(patcher.stop)

    def view(self):
        return ar.haplotype_genome_view(self.cwd, self.genome_dir, "paternal")

    def test_both_make_the_haplotypes_link(self):
        real_symlink = os.symlink
        started, other = [], []

        def symlink(*args, **kwargs):
            real_symlink(*args, **kwargs)
            if not started:  # the other run starts now, and goes all the way
                started.append(True)
                other.append(self.view())

        with mock.patch.object(ar.os, "symlink", symlink):
            view = self.view()
        self.assertEqual(other, [view])
        self.assertEqual(os.path.realpath(view), os.path.realpath(self.genome_dir))
        self.assertEqual(os.listdir(os.path.dirname(view)), [os.path.basename(view)])

    def test_a_link_the_other_run_made_meanwhile_is_taken(self):
        # on a file system that will not rename over the other run's link
        view = self.view()
        target = os.readlink(view)
        os.remove(view)

        def replace(src, dst):
            os.symlink(target, dst)
            raise FileExistsError(dst)

        with mock.patch.object(ar.os, "replace", replace):
            self.assertEqual(self.view(), view)
        self.assertEqual(os.listdir(os.path.dirname(view)), [os.path.basename(view)])

    def test_any_other_failure_still_stops_the_run(self):
        view = self.view()
        os.remove(view)
        with mock.patch.object(ar.os, "replace", side_effect=PermissionError(view)):
            with self.assertRaises(PermissionError):
                self.view()
        self.assertEqual(os.listdir(os.path.dirname(view)), [])  # no temporary link left

    def test_both_record_the_index_they_found_built(self):
        view = self.view()
        index = os.path.join(self.cwd, "genome_library", f"NRG_2_{os.path.basename(view)}")
        os.makedirs(index)
        open(os.path.join(index, "NRG_chr1_1.bin"), "w").close()
        real_dump = json.dump
        started = []

        def dump(*args, **kwargs):
            real_dump(*args, **kwargs)
            if not started:  # the other run records it now
                started.append(True)
                ar.stamp_built_haplotype_indexes(self.cwd, view, [])

        with mock.patch.object(ar.json, "dump", dump), \
                contextlib.redirect_stdout(io.StringIO()) as out:
            ar.stamp_built_haplotype_indexes(self.cwd, view, [])
        self.assertEqual(out.getvalue(), "")  # no "Note: could not record ..."
        self.assertEqual(ar._read_index_source(index)["fingerprint"], ar.genome_fingerprint(view))
        self.assertEqual(sorted(os.listdir(index)), [ar.INDEX_SOURCE_FILENAME, "NRG_chr1_1.bin"])


# --------------------------------------------------------------------------- #
# assembly_search, run in-process
# --------------------------------------------------------------------------- #
class _FakeCompleteSearch:
    """Stands in for `subprocess.call` of a haplotype's complete-search, keeping the
    one thing that decides which sequence a search reads: the index is found by the
    genome folder's NAME (genome_library/<PAM>_<N>_<name>) and built from the folder
    only when none of that name exists. The output records the sequence searched."""

    def __init__(self, cwd):
        self.cwd = cwd
        self.calls = []

    def call(self, cmd, shell=True, cwd=None):
        argv = shlex.split(cmd)
        opt = {argv[i]: argv[i + 1] for i in range(len(argv) - 1) if argv[i].startswith("--")}
        name = os.path.basename(os.path.normpath(opt["--genome"]))
        lib = os.path.join(self.cwd, "genome_library")
        os.makedirs(lib, exist_ok=True)
        found = sorted(d for d in os.listdir(lib) if re.fullmatch(r"[^_]+_[0-9]+_" + re.escape(name), d))
        index = os.path.join(lib, found[0] if found else f"NRG_2_{name}")
        if not found:
            os.makedirs(index)
            with open(os.path.join(os.path.realpath(opt["--genome"]), "chr1.fa")) as fh:
                seq = fh.read().split("\n")[1]
            with open(os.path.join(index, "NRG_chr1_1.bin"), "w") as fh:
                fh.write(seq)
        with open(os.path.join(index, "NRG_chr1_1.bin")) as fh:
            searched = fh.read()
        out = os.path.join(self.cwd, "Results", opt["--output"])
        os.makedirs(out)
        with open(os.path.join(out, ".command_line.txt"), "w") as fh:
            fh.write("input_command\t" + " ".join(argv[1:]) + "\n")
        with open(os.path.join(out, ".Params.txt"), "w") as fh:
            fh.write("Pam\tNRG\n")
        with open(os.path.join(out, "g+NRG_genome+__4+2_integrated_results.tsv"), "w") as fh:
            fh.write(f"searched\n{searched}\n")
        open(os.path.join(out, ar.LOG_ERROR_NO_CHECK_FILENAME), "w").close()
        self.calls.append((opt["--output"], opt["--genome"], os.path.basename(index)))
        return 0


class AssemblySearchSearchesEachHaplotypeOnItsOwnGenome(_Workdir):
    SEQS = {"W1": {"paternal": "CTTACAGTTGCTTTTACCACTGG", "maternal": "CTTACTGTTGCCTTTACCACTGG"},
            "W2": {"paternal": "GGGGAGTTGCTTTTATCACTGG", "maternal": "AAAACAGTTGCTTTTATCACAGG"}}

    def setUp(self):
        super().setUp()
        for ind, seqs in self.SEQS.items():
            _make_bundle(self.cwd, ind, seqs)
        self.guide = os.path.join(self.cwd, "guide.txt")
        with open(self.guide, "w") as fh:
            fh.write("CTAACAGTTGCTTTTATCACNNN\n")
        self.pam = os.path.join(self.cwd, "20bp-NRG-SpCas9.txt")
        shutil.copy(os.path.join(PAMS_DIR, "20bp-NRG-SpCas9.txt"), self.pam)
        self.liftover = {}
        for hap in pa.HAPLOTYPES:
            for kind, ext in (("chain", "chain"), ("chrom-alias", "chromAlias.txt")):
                path = os.path.join(self.cwd, "LiftoverFiles", f"{hap}.{ext}")
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as fh:
                    fh.write("x\n")
                self.liftover[f"--{kind}-{hap}"] = path
        self.engine = _FakeCompleteSearch(self.cwd)
        for target, value in (
            ("current_working_directory", self.cwd + "/"),
            ("subprocess", types.SimpleNamespace(call=self.engine.call)),
            ("check_liftover_available", lambda: None),
            ("check_impg_available", lambda: None),
            ("reconcile_haplotypes", self._reconcile),
            ("build_combined_report", lambda *a, **k: "report.zip"),
        ):
            patcher = mock.patch.object(crisprme, target, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.reconciled = []

    def _reconcile(self, haplotypes, combined_output, **kwargs):
        self.reconciled.append({h: v["genome_dir"] for h, v in haplotypes.items()})
        return pd.DataFrame(), {}

    def run_search(self, output, individual=None, explicit=None, bulges=1):
        args = ["--guide", self.guide, "--pam", self.pam, "--mm", "4",
                "--bDNA", str(bulges), "--bRNA", str(bulges), "--output", output]
        if individual:
            args += ["--assembly-individual", individual]
        else:
            pat, mat = explicit
            args += ["--genome-paternal", pat, "--genome-maternal", mat]
            for flag, path in self.liftover.items():
                args += [flag, path]
        with mock.patch.object(crisprme, "input_args", ["crisprme.py", "assembly-search", *args]), \
                contextlib.redirect_stdout(io.StringIO()):
            crisprme.assembly_search()

    def searched(self, output):
        with open(os.path.join(self.cwd, "Results", output,
                               "g+NRG_genome+__4+2_integrated_results.tsv")) as fh:
            return fh.read().split("\n")[1]

    def assertOwnSequences(self, output, individual):
        for hap in pa.HAPLOTYPES:
            self.assertEqual(self.searched(f"{output}_{hap}"), self.SEQS[individual][hap],
                             f"{output}_{hap} searched another haplotype's sequence")

    def test_each_haplotype_of_a_registered_individual(self):
        self.run_search("job", individual="W1")
        self.assertOwnSequences("job", "W1")
        self.assertEqual(len(self.engine.calls), 2)
        for (_, _, index), hap in zip(self.engine.calls, pa.HAPLOTYPES):
            self.assertRegex(index, rf"^NRG_2_W1_{hap}_[0-9a-f]{{8}}$")
        # the reconciliation still reads the genome folders themselves
        self.assertEqual(self.reconciled[0], {h: self.genome("W1", h) for h in pa.HAPLOTYPES})

    def test_the_explicit_folder_form_as_the_web_submits_it(self):
        self.run_search("job", explicit=(self.genome("W1", "paternal"), self.genome("W1", "maternal")))
        self.assertOwnSequences("job", "W1")

    def test_a_second_individual_does_not_reuse_the_first_ones_index(self):
        self.run_search("one", individual="W1")
        self.run_search("two", individual="W2")
        self.assertOwnSequences("two", "W2")

    def test_a_result_written_by_the_colliding_search_is_not_reused(self):
        # what 2.8.0 and earlier left: a completed maternal search recorded with the
        # bare bundle folder, which found its index by the name "genome"
        self.run_search("job", individual="W1")
        out = os.path.join(self.cwd, "Results", "job_maternal")
        with open(os.path.join(out, ".command_line.txt")) as fh:
            line = fh.read()
        view = re.search(r"--genome (\S+)", line).group(1)
        with open(os.path.join(out, ".command_line.txt"), "w") as fh:
            fh.write(line.replace(view, self.genome("W1", "maternal")))
        os.remove(os.path.join(self.cwd, "Results", "job_combined", "job_combined_hg38.tsv"))
        self.engine.calls.clear()
        self.run_search("job", individual="W1")
        self.assertEqual([c[0] for c in self.engine.calls], ["job_maternal"])

    def test_a_finished_haplotype_search_is_still_reused(self):
        self.run_search("job", individual="W1")
        os.remove(os.path.join(self.cwd, "Results", "job_combined", "job_combined_hg38.tsv"))
        self.engine.calls.clear()
        self.run_search("job", individual="W1")
        self.assertEqual(self.engine.calls, [])

    def test_a_bundle_imported_again_with_other_content_is_indexed_again(self):
        self.run_search("one", individual="W1")
        _write_genome(self.genome("W1", "maternal"), "CTAACAGTTGCTTTTATCACTGG", mtime=time.time() - 3600)
        self.SEQS = {"W1": dict(self.SEQS["W1"], maternal="CTAACAGTTGCTTTTATCACTGG")}
        self.run_search("two", individual="W1")
        self.assertOwnSequences("two", "W1")

    def test_a_legacy_genomes_pair_is_searched_as_before(self):
        pat = os.path.join(self.cwd, "Genomes", "L_paternal")
        mat = os.path.join(self.cwd, "Genomes", "L_maternal")
        _write_genome(pat, "ACGTA")
        _write_genome(mat, "TTTTA")
        self.run_search("job", explicit=(pat, mat))
        self.assertEqual([c[1] for c in self.engine.calls], [pat, mat])
        self.assertEqual([c[2] for c in self.engine.calls], ["NRG_2_L_paternal", "NRG_2_L_maternal"])

    def test_the_legacy_layouts_index_of_the_individual_is_not_used(self):
        # what an earlier CRISPRme left for Genomes/W1_paternal and Genomes/W1_maternal
        # (the layout before 2.6.0, or a custom upload): indexes of other content under
        # those names, recording nothing
        legacy = {hap: os.path.join(self.cwd, "genome_library", f"NRG_2_W1_{hap}")
                  for hap in pa.HAPLOTYPES}
        for hap, index in legacy.items():
            os.makedirs(index)
            with open(os.path.join(index, "NRG_chr1_1.bin"), "w") as fh:
                fh.write(f"legacy {hap}")
        self.run_search("job", individual="W1")
        self.assertOwnSequences("job", "W1")
        for hap, index in legacy.items():  # left to those folders
            with open(os.path.join(index, "NRG_chr1_1.bin")) as fh:
                self.assertEqual(fh.read(), f"legacy {hap}")

    def test_the_legacy_layouts_folder_does_not_find_the_bundles_index(self):
        # complete-search of Genomes/W1_paternal (from the Reference genome tab, say)
        # picks an index by that name alone, reading no record
        self.run_search("job", individual="W1")
        lib = os.path.join(self.cwd, "genome_library")
        for hap in pa.HAPLOTYPES:
            own = ar.haplotype_search_name(self.cwd, self.genome("W1", hap), hap)
            self.assertEqual(_shipped_resolve(lib, own, 2), f"NRG_2_{own}")
            self.assertEqual(_shipped_resolve(lib, f"W1_{hap}", 2), "")


# --------------------------------------------------------------------------- #
# The pipeline: the name a linked genome keeps, and the 0-bulge scan
# --------------------------------------------------------------------------- #
class APipelineGenomeKeepsTheNameItIsGivenUnder(_Workdir):
    """The pipeline names the search's files and index after the genome folder it is
    given, and its post-analysis takes that name back from `$ref_folder`. Resolved
    in full, a link named after the haplotype became "genome" there: the
    post-analysis then looked for genome_<...>.targets.txt, which the search had
    written as <haplotype>_<...>.targets.txt. Runs the shipped assignment."""

    def ref_folder(self, given):
        with open(SHELL_SCRIPT) as fh:
            line = next(ln for ln in fh if ln.startswith("ref_folder="))
        return subprocess.run(["bash", "-c", line + 'printf "%s" "$ref_folder"', "bash", given],
                              capture_output=True, text=True, check=True).stdout

    def test_a_link_keeps_its_own_name(self):
        _make_bundle(self.cwd, "W1", {"paternal": "ACGT", "maternal": "TTTT"})
        view = os.path.join(self.cwd, "views", "W1_paternal")
        os.makedirs(os.path.dirname(view))
        os.symlink(self.genome("W1", "paternal"), view)
        got = self.ref_folder(view)
        self.assertEqual(os.path.basename(got), "W1_paternal")
        self.assertEqual(os.path.realpath(got), os.path.realpath(self.genome("W1", "paternal")))

    def test_a_folder_is_resolved_as_before(self):
        genome = os.path.join(self.cwd, "Genomes", "hg38")
        _write_genome(genome, "ACGT")
        for given in (genome, genome + "/", os.path.join(self.cwd, "Genomes", ".", "hg38")):
            self.assertEqual(self.ref_folder(given), os.path.realpath(genome))


def _shipped_resolver():
    """`_resolve_index_dir`, extracted from the shipped script (not copied), as
    test_index_resolution.py extracts it."""
    with open(SHELL_SCRIPT) as fh:
        text = fh.read()
    start = text.index("\t_resolve_index_dir() {")
    end = text.index("\n\t}\n", start) + len("\n\t}\n")
    return "\n".join(ln[1:] if ln.startswith("\t") else ln for ln in text[start:end].splitlines())


def _shipped_resolve(lib, suffix, bmax):
    """The index the shipped `_resolve_index_dir` picks in `lib` for the genome named
    `suffix` and the NRG PAM (its name, "" for none)."""
    script = f"prune_pam=0\npos=3\n{_shipped_resolver()}\n_resolve_index_dir \"$1\" NRG \"$2\" \"$3\"\n"
    out = subprocess.run(["bash", "-c", script, "bash", lib, suffix, str(bmax)],
                         capture_output=True, text=True, env=dict(os.environ, LC_ALL="C"))
    return os.path.basename(out.stdout.strip())


class TheIndexOfAnotherGenomeIsNeverPicked(_Workdir):
    """An index is picked by the genome's name, and the scan globbed
    <PAM>_<N>_<name>: that also matches <PAM>_<N>_<prefix>_<name>, the index of a
    genome whose name merely ends in _<name> -- individual "trio_P1"'s haplotype for
    "P1"'s, mini_hg38's for hg38 -- and a smaller budget won over the genome's own."""

    def resolve(self, suffix, bmax, *names):
        lib = os.path.join(self.cwd, "genome_library")
        for name in names:
            os.makedirs(os.path.join(lib, name), exist_ok=True)
            open(os.path.join(lib, name, "chr1.bin"), "w").close()
        return _shipped_resolve(lib, suffix, bmax)

    def test_the_genomes_own_index_wins(self):
        self.assertEqual(self.resolve("P1_paternal", 2, "NRG_2_trio_P1_paternal", "NRG_3_P1_paternal"),
                         "NRG_3_P1_paternal")

    def test_none_is_picked_when_the_genome_has_none(self):
        self.assertEqual(self.resolve("hg38", 2, "NRG_2_mini_hg38"), "")

    def test_a_variant_index_is_still_found(self):
        self.assertEqual(self.resolve("hg38+hg38_1000G", 2, "NRG_2_hg38+hg38_1000G", "NRG_2_hg38"),
                         "NRG_2_hg38+hg38_1000G")


_STUB_ENGINE = r"""#!/bin/bash
printf '%s\n' "$*" >> "$STUB_LOG"
# the reference search is the last step a reference-only run needs here
[ "$1" = search ] && exit 1
exit 0
"""


def _gnu_realpath():
    """The script runs `realpath` on the "_" placeholder of a reference-only search,
    which only GNU realpath resolves (BSD's writes to stderr, a failed stage)."""
    try:
        return subprocess.run(["realpath", "_no_such_path_"], capture_output=True).returncode == 0
    except OSError:
        return False


@unittest.skipUnless(shutil.which("python") and _gnu_realpath(),
                     "the shipped script needs python and GNU realpath on PATH")
class ZeroBulgeScanReadsTheGenomeItIsGiven(_Workdir):
    """The 0-bulge search scanned Genomes/<genome name>/, whatever folder it was
    given: for a haplotype genome, a folder that does not exist, so the search came
    back empty. Runs the shipped script the way complete-search launches it, with a
    stub engine that records its calls and stops the run at the search. The genome is
    handed over as haplotype_genome_view does: a link, named after the haplotype, to
    the bundle's genome folder."""

    def calls(self, bulges):
        _make_bundle(self.cwd, "W1", {"paternal": "ACGTACGT", "maternal": "TTTT"})
        view = os.path.join(self.cwd, "views", "W1_paternal")
        os.makedirs(os.path.dirname(view), exist_ok=True)
        if not os.path.islink(view):
            os.symlink(self.genome("W1", "paternal"), view)
        work = os.path.join(self.cwd, f"run{bulges}")
        out = os.path.join(work, "out")
        os.makedirs(out)
        bindir = os.path.join(work, "bin")
        os.makedirs(bindir)
        with open(os.path.join(bindir, "crispritz.py"), "w") as fh:
            fh.write(_STUB_ENGINE)
        os.chmod(os.path.join(bindir, "crispritz.py"), 0o755)
        pam = os.path.join(work, "20bp-NRG-SpCas9.txt")
        shutil.copy(os.path.join(PAMS_DIR, "20bp-NRG-SpCas9.txt"), pam)
        guides = os.path.join(work, "guides.txt")
        with open(guides, "w") as fh:
            fh.write("CTAACAGTTGCTTTTATCACNNN\n")
        empty = os.path.join(work, "empty.txt")
        open(empty, "w").close()
        log = os.path.join(work, "calls.txt")
        args = [view, "_", guides, pam, empty, "_", bulges + 1, 4, bulges, bulges, 3, out,
                HERE, 1, self.cwd, empty, "_", 1, 0, "_", "mm+bulges", "mm+bulges,mm",
                "False", "PASS,.", "_", 6, 0]
        env = dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}",
                   STUB_LOG=log, LC_ALL="C", LANG="C", CRISPRME_NO_THREAD_GUARD="1")
        with open(os.path.join(out, "log_error.txt"), "w") as err:
            subprocess.run(["bash", SHELL_SCRIPT, *map(str, args)], cwd=work, env=env,
                           stdout=subprocess.DEVNULL, stderr=err, timeout=300)
        with open(log) as fh:
            return [line.split() for line in fh]

    def test_without_bulges_the_given_folder_is_scanned(self):
        search = [c for c in self.calls(0) if c[0] == "search"]
        self.assertEqual(len(search), 1, search)
        self.assertIn("-r", search[0])
        self.assertEqual(os.path.realpath(search[0][1]), os.path.realpath(self.genome("W1", "paternal")))

    def test_with_bulges_the_index_is_the_haplotypes_own(self):
        calls = self.calls(1)
        self.assertEqual([c[:2] for c in calls if c[0] == "index-genome"], [["index-genome", "W1_paternal"]])
        search = [c for c in calls if c[0] == "search"]
        self.assertEqual(os.path.basename(search[0][1].rstrip("/")), "NRG_2_W1_paternal")

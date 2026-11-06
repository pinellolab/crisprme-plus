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
    on its own sequence, and a result written by the colliding search is not reused; a
    haplotype search that fails ends it in one line, exit status 1, not a traceback;
  * the record of which genome an index holds (.genome_source): an index of the same
    name built from another folder, or before the folder's FASTA files changed (names,
    sizes and modification times are compared, not the sequence), is removed and built
    again, one that records nothing (built by an earlier CRISPRme, or downloaded) is kept;
  * the pipeline keeps a linked genome's own name, picks only an index named for the
    genome itself (not one whose genome name merely ends in it), and its 0-bulge scan
    reads the genome folder it is given;
  * complete-search itself (run in-process as the CLI dispatch runs it, the pipeline
    replaced by a stand-in that takes or builds the index by the genome's name as the
    pipeline does): a search with bulges of a folder outside Genomes/ named like a
    genome in Genomes/, or like another folder whose index is here, is refused, the
    index it builds records its folder, and a genome in Genomes/ is searched as before;
    one haplotype of a bundle searched by its own folder is refused the index of its
    folder's name that records no folder (the one assembly-search shared before 2.8.1)
    and could serve the search, but not the link assembly-search searches it through;
    each refusal, followed as it says, runs the search with the bulges it was refused
    with (given none, a search takes the bulges the indexes installed allow), and
    under --index-path, which the pipeline builds no index in, without it (--prune-pam,
    which needed it before 2.9.0, stays);
  * the record names a folder of the machine that built the index: publish-index
    leaves it out of the archive, with and without pigz, and download drops it;
  * runs that start together: a stale index both find is taken away from its name in
    one step, and the data folders the pipeline makes are made without an error;
  * deleting an individual in Settings removes its haplotypes' indexes and links too.
"""
import ast
import contextlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import time
import types
import typing
import unittest
from unittest import mock

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CRISPRME_PY = os.path.join(os.path.dirname(HERE), "crisprme.py")
SHELL_SCRIPT = os.path.join(HERE, "submit_job_automated_new_multiple_vcfs.sh")
PAMS_DIR = os.path.join(os.path.dirname(HERE), "PAMs")
sys.path.insert(0, HERE)

import assembly_reconcile as ar  # noqa: E402
import crisprme_hf as hf  # noqa: E402
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
    in both. After the individual's FASTA files change, they all find its recorded
    index stale and remove it: the one that came second stopped with
    FileNotFoundError, as the index was gone."""

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

    def stale_index(self, view, bins=1):
        """An index of the haplotype recorded before its FASTA files changed (an
        individual imported again with other files): both runs find it stale."""
        index = os.path.join(self.cwd, "genome_library", f"NRG_2_{os.path.basename(view)}")
        os.makedirs(index)
        for i in range(bins):
            open(os.path.join(index, f"NRG_chr{i}_1.bin"), "w").close()
        with open(os.path.join(index, ar.INDEX_SOURCE_FILENAME), "w") as fh:
            json.dump({"genome": os.path.join("Assemblies", "W1", "paternal", "genome"),
                       "fingerprint": "0" * 40}, fh)
        return index

    def test_a_stale_index_the_other_run_removed_first_is_no_failure(self):
        # the second removal stopped the run with "An index of ... holds another
        # genome and could not be removed (FileNotFoundError ...)", the index gone
        view = self.view()
        self.stale_index(view)
        real_read, started = ar._read_index_source, []

        def read(index_dir):
            record = real_read(index_dir)
            if not started:  # the other run, which found it stale too, removes it now
                started.append(True)
                ar.clear_stale_haplotype_indexes(self.cwd, view)
            return record

        with mock.patch.object(ar, "_read_index_source", read), \
                contextlib.redirect_stdout(io.StringIO()):
            present = ar.clear_stale_haplotype_indexes(self.cwd, view)
        self.assertEqual(present, [])
        self.assertEqual(os.listdir(os.path.join(self.cwd, "genome_library")), [])

    def test_a_stale_index_leaves_its_name_before_it_is_emptied(self):
        # a search starting meanwhile finds no index of the haplotype, never part of one
        view = self.view()
        index = self.stale_index(view, bins=3)
        real_rmtree, removing = shutil.rmtree, []

        def rmtree(path, *args, **kwargs):
            removing.append((path, os.path.lexists(index)))
            real_rmtree(path, *args, **kwargs)

        with mock.patch.object(ar.shutil, "rmtree", rmtree), \
                contextlib.redirect_stdout(io.StringIO()) as out:
            ar.clear_stale_haplotype_indexes(self.cwd, view)
        self.assertEqual(len(removing), 1)
        self.assertNotEqual(removing[0][0], index)
        self.assertFalse(removing[0][1], "the index was emptied under its own name")
        self.assertEqual(os.listdir(os.path.join(self.cwd, "genome_library")), [])
        self.assertIn(f"Removing index {os.path.basename(index)}", out.getvalue())


class RunsStartingTogetherBothMakeTheDataFolders(_Workdir):
    """The pipeline makes Results/, Dictionaries/, Genomes/ and genome_library/ when
    they are missing: a test, then a plain mkdir. Two runs starting together on a data
    folder without genome_library/ (the first searches after an install, say) both
    found it missing, and the one that made it second wrote "mkdir: cannot create
    directory '/DATA/genome_library': File exists" to stderr, which fails the run
    ("ERROR: Genome enrichment failed on _"). Runs the shipped block, with the other
    run making each of those folders right between the test and the mkdir."""

    SHARED = ("Results", "Dictionaries", "Genomes", "genome_library")

    def block(self):
        with open(SHELL_SCRIPT) as fh:
            text = fh.read()
        start = text.index("\t# ensure that fundamental folders are available")
        return text[start:text.index("\t# START STEP 1 - Genome enrichment", start)]

    def test_a_folder_the_other_run_made_meanwhile_is_no_error(self):
        real_mkdir = shutil.which("mkdir")
        bindir = os.path.join(self.cwd, "bin")
        os.makedirs(bindir)
        with open(os.path.join(bindir, "mkdir"), "w") as fh:
            fh.write("#!/bin/bash\n"
                     'for a in "$@"; do\n'
                     f'\tcase "${{a##*/}}" in {"|".join(self.SHARED)}) "{real_mkdir}" -p "$a" ;; esac\n'
                     "done\n"
                     f'exec "{real_mkdir}" "$@"\n')
        os.chmod(os.path.join(bindir, "mkdir"), 0o755)
        data, out = os.path.join(self.cwd, "data"), os.path.join(self.cwd, "out")
        os.makedirs(data)
        os.makedirs(out)
        env = dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}",
                   current_working_directory=data, output_folder=out)
        run = subprocess.run(["bash", "-c", self.block()], env=env, capture_output=True, text=True)
        self.assertEqual(run.stderr, "")
        self.assertEqual(run.returncode, 0)
        self.assertEqual(sorted(os.listdir(data)), sorted(self.SHARED))
        self.assertEqual(os.listdir(out), ["crispritz_targets"])

    def test_crisprme_py_takes_a_folder_the_other_run_made_meanwhile(self):
        # crisprme.py makes Genomes/, Results/, ... first, the same way (FileExistsError)
        made, real_exists = [], os.path.exists

        def exists(path):
            if os.path.dirname(path) == self.cwd and os.path.basename(path) in crisprme.CRISPRMEDIRS:
                os.mkdir(path)  # made by the other run, once this one has looked
                made.append(path)
                return False
            return real_exists(path)

        with mock.patch.object(crisprme, "current_working_directory", self.cwd + os.sep), \
                mock.patch.object(crisprme.os.path, "exists", exists):
            crisprme.check_crisprme_dirtree()
        self.assertEqual(len(made), len(crisprme.CRISPRMEDIRS))
        self.assertEqual(sorted(os.listdir(self.cwd)), sorted(crisprme.CRISPRMEDIRS))


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

    def test_a_failed_haplotype_search_ends_in_one_line(self):
        # its complete-search has said why: assembly-search then ended in an OSError
        # traceback; the web's job reads this line as the failure, as it read that
        def call(cmd, shell=True, cwd=None):
            return 1 if "--output job_paternal " in cmd else self.engine.call(cmd, shell, cwd)

        err = io.StringIO()
        with mock.patch.object(crisprme, "subprocess", types.SimpleNamespace(call=call)), \
                contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as stop:
            self.run_search("job", individual="W1")
        self.assertEqual(stop.exception.code, 1)
        log = os.path.join(self.cwd, "Results", "job_paternal", "log_error.txt")
        self.assertEqual(err.getvalue(), f"Error: Haplotype search failed for --output "
                                         f"job_paternal! See {log} for details\n")
        self.assertEqual(self.reconciled, [])


# --------------------------------------------------------------------------- #
# complete-search, run in-process as the CLI dispatch runs it
# --------------------------------------------------------------------------- #
_DISPATCH_MARKER = "\nif len(sys.argv) < 2:"


def _dispatch(subcommand):
    """The statement crisprme.py's CLI dispatch runs for `subcommand`, verbatim (as
    test_annotation_tmpdirs takes it)."""
    with open(CRISPRME_PY) as fh:
        tail = fh.read().split(_DISPATCH_MARKER, 1)[1]
    m = re.search(rf'^elif sys\.argv\[1\] == "{re.escape(subcommand)}":[^\n]*\n'
                  r'((?:    [^\n]*\n)+)', tail, re.M)
    return textwrap.dedent(m.group(1))


class _FakePipeline:
    """Stands in for `subprocess.call` of submit_job_automated_new_multiple_vcfs.sh,
    keeping what decides which sequence a reference search reads: with bulges, the
    index of the genome's NAME (<PAM>_<N>_<name>, or 20-<PAM>_<N>_<name> as it builds one
    since 2.9.0, the guide length and the side of the PAM in its name; N at least the
    budget, the smallest first, then the newer name), in --index-path or genome_library/,
    built from the folder only when there is none, with the log line the pipeline prints
    then -- under --index-path, which it builds none in, it fails instead, with the
    pipeline's error; without bulges, the folder itself. The bulges each search ran with
    are kept in `bulges`."""

    def __init__(self, cwd):
        self.cwd = cwd
        self.calls, self.built, self.bulges = [], [], {}

    def call(self, cmd, shell=True, stderr=None, stdout=None, **kwargs):
        argv = cmd.split()
        genome, want, bdna, out, lib = argv[1], int(argv[7]), int(argv[9]), argv[12], argv[25]
        self.bulges[os.path.basename(out)] = (bdna, int(argv[10]))
        index_path = lib != "_"
        lib = lib if index_path else os.path.join(self.cwd, "genome_library")
        name = os.path.basename(genome)
        os.makedirs(lib, exist_ok=True)
        fits = []
        for d in os.listdir(lib):
            m = re.fullmatch(r"(20-)?NRG_([0-9]+)_" + re.escape(name), d)
            if m and int(m.group(2)) >= want:
                fits.append((int(m.group(2)), not m.group(1), d))
        found = [d for _n, _older, d in sorted(fits)]
        if found:
            index = os.path.join(lib, found[0])
            stdout.write(f"Reference index present ({found[0]})\n")
        elif index_path:
            stderr.write(f"ERROR: no matching reference index under --index-path '{lib}'\n")
            stderr.flush()
            return 1
        else:
            index = os.path.join(lib, f"20-NRG_{want}_{name}")
            stdout.write("Indexing reference genome\n")  # the pipeline's, checked in ZeroBulgeScan...
            os.makedirs(index)
            shutil.copy(os.path.join(genome, "chr1.fa"), os.path.join(index, "NRG_chr1_1.bin"))
            self.built.append(os.path.basename(index))
        stdout.flush()
        source = os.path.join(index, "NRG_chr1_1.bin") if bdna else os.path.join(genome, "chr1.fa")
        with open(source) as fh, open(os.path.join(out, "searched.txt"), "w") as dst:
            dst.write(fh.read().split("\n")[1])
        self.calls.append(os.path.basename(out))
        return 0


class ACompleteSearchTakesOnlyTheIndexOfItsOwnFolder(_Workdir):
    """complete-search names a genome's index (genome_library/<PAM>_<N>_<name>) after
    the genome folder's name, and takes any index of that name. Searched with
    bulges, projB/ref took the index projA/ref's search had built and reported
    projA/ref's sites; other/hg38 took Genomes/hg38's; one haplotype of a bundle,
    searched by its own folder, took the index that assembly-search had shared
    between all haplotypes before 2.8.1 (NRG_2_genome). Without bulges a search
    scans the folder it is given."""

    SEQ = {"paternal": "CTTACAGTTGCTTTTACCACTGG", "maternal": "CTTACTGTTGCCTTTACCACTGG",
           "reference": "CTAACAGTTGCTTTTATCACTGG"}

    def setUp(self):
        super().setUp()
        here = os.getcwd()
        os.chdir(self.cwd)
        self.addCleanup(os.chdir, here)
        for folder, seq in (("projA/ref", "paternal"), ("projB/ref", "maternal"),
                            ("Genomes/hg38", "reference"), ("other/hg38", "paternal")):
            _write_genome(os.path.join(self.cwd, folder), self.SEQ[seq])
        os.makedirs("PAMs")
        shutil.copy(os.path.join(PAMS_DIR, "20bp-NRG-SpCas9.txt"), "PAMs")
        with open("guide.txt", "w") as fh:
            fh.write("CTAACAGTTGCTTTTATCACNNN\n")
        self.lib = os.path.join(self.cwd, "genome_library")
        self.pipeline = _FakePipeline(self.cwd)
        env = dict(os.environ, CRISPRME_SKIP_SCORER_ENV="1", CRISPRME_SKIP_MEM_CHECK="1")
        for patcher in (mock.patch.dict(os.environ, env, clear=True),
                        mock.patch.object(crisprme, "current_working_directory", self.cwd + os.sep)):
            patcher.start()
            self.addCleanup(patcher.stop)

    def search(self, output, genome, bulges=1, *extra):
        """`crisprme.py complete-search ...` as the CLI dispatch runs it (with no
        --bDNA/--bRNA for `bulges` None); the sequence searched. Its stderr is kept
        as self.err."""
        flags = [] if bulges is None else ["--bDNA", str(bulges), "--bRNA", str(bulges)]
        argv = ["crisprme.py", "complete-search", "--genome", genome, "--guide", "guide.txt",
                "--pam", "PAMs/20bp-NRG-SpCas9.txt", "--mm", "4", *flags,
                "--output", output, "--thread", "1", *extra]
        err = io.StringIO()
        try:
            with mock.patch.object(crisprme, "input_args", argv), \
                    mock.patch.object(crisprme.subprocess, "call", side_effect=self.pipeline.call), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
                exec(_dispatch("complete-search"), crisprme.__dict__)
        finally:
            self.err = err.getvalue()
        with open(os.path.join(self.cwd, "Results", output, "searched.txt")) as fh:
            return fh.read()

    def record(self, index):
        return ar._read_index_source(os.path.join(self.lib, index))

    def follow(self, advice, output, genome, bulges, *extra):
        """Runs a refused search again as its message, `advice`, says to: without
        --index-path, when it says so, with the --bDNA/--bRNA it gives, and otherwise as
        it was run (`bulges`, `extra`), --prune-pam included. The sequence searched."""
        again = re.search(r"run this search again( without --index-path"
                          r"( and --prune-pam \(which needs it\))?)?(?: and)?"
                          r"(?: with --bDNA (\d+) --bRNA (\d+))?", advice)
        self.assertIsNotNone(again, f"no 'run this search again' in: {advice}")
        if again.group(1):
            at = extra.index("--index-path")
            extra = extra[:at] + extra[at + 2:]
        if again.group(2):
            extra = tuple(a for a in extra if a != "--prune-pam")
        if again.group(3) is not None:
            self.assertEqual(again.group(3), again.group(4))  # (search() gives both alike)
            bulges = int(again.group(3))
        return self.search(output, genome, bulges, *extra)

    def test_another_folder_of_the_same_name_is_refused(self):
        self.assertEqual(self.search("a", "projA/ref"), self.SEQ["paternal"])
        with self.assertRaises(SystemExit):
            self.search("b", "projB/ref")
        self.assertIn(f"{self.lib}/20-NRG_2_ref holds the genome in {self.cwd}/projA/ref, "
                      f"not {self.cwd}/projB/ref", self.err)
        self.assertEqual(self.pipeline.calls, ["a"])
        self.assertEqual(os.listdir(os.path.join(self.cwd, "Results", "b")), [])  # rerun as it is

    def test_the_folder_it_was_built_from_takes_it_again(self):
        self.search("a", "projA/ref")
        self.assertEqual(self.record("20-NRG_2_ref"),
                         {"genome": os.path.join("projA", "ref"),
                          "fingerprint": ar.genome_fingerprint(os.path.join(self.cwd, "projA", "ref"))})
        self.assertEqual(self.search("a2", "projA/ref"), self.SEQ["paternal"])
        self.assertEqual(self.pipeline.built, ["20-NRG_2_ref"])

    def test_a_folder_named_like_a_genome_in_genomes_is_refused(self):
        # whether that genome's index is built yet or not: else this one would build
        # it, for the searches of Genomes/hg38 (the web's too) to take
        for _ in range(2):
            with self.assertRaises(SystemExit):
                self.search("o", "other/hg38")
            self.assertIn(f"--genome {self.cwd}/other/hg38 has the name of another genome "
                          f"folder, {self.cwd}/Genomes/hg38", self.err)
            self.search(f"h{len(self.pipeline.calls)}", "Genomes/hg38")
        self.assertEqual(self.pipeline.calls, ["h0", "h1"])

    def test_without_bulges_such_a_folder_is_searched(self):
        # (it reads the folder itself, and takes nothing by its name)
        self.search("h", "Genomes/hg38")
        self.assertEqual(self.search("o0", "other/hg38", 0), self.SEQ["paternal"])

    def test_a_genome_in_genomes_is_searched_as_before(self):
        self.assertEqual(self.search("h", "Genomes/hg38"), self.SEQ["reference"])
        self.assertEqual(self.search("h2", "Genomes/hg38"), self.SEQ["reference"])
        self.assertEqual(self.pipeline.built, ["20-NRG_2_hg38"])
        self.assertEqual(os.listdir(os.path.join(self.lib, "20-NRG_2_hg38")), ["NRG_chr1_1.bin"])  # no record

    def test_an_index_that_records_nothing_is_taken_as_before(self):
        # (built by an earlier CRISPRme, or downloaded: which folder it holds is unknown)
        os.makedirs(os.path.join(self.lib, "NRG_2_ref"))
        shutil.copy(os.path.join(self.cwd, "projB", "ref", "chr1.fa"),
                    os.path.join(self.lib, "NRG_2_ref", "NRG_chr1_1.bin"))
        self.assertEqual(self.search("b", "projB/ref"), self.SEQ["maternal"])
        self.assertIsNone(self.record("NRG_2_ref"))  # this run did not build it

    def test_an_index_built_before_its_folder_changed_is_built_again(self):
        self.search("a", "projA/ref")
        _write_genome(os.path.join(self.cwd, "projA", "ref"), self.SEQ["maternal"],
                      mtime=time.time() - 3600)
        self.assertEqual(self.search("a2", "projA/ref"), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.built, ["20-NRG_2_ref", "20-NRG_2_ref"])

    def test_a_genome_in_genomes_is_refused_an_index_built_from_another_folder(self):
        # projA/ref was searched before a Genomes/ref was there
        self.search("a", "projA/ref")
        _write_genome(os.path.join(self.cwd, "Genomes", "ref"), self.SEQ["maternal"])
        with self.assertRaises(SystemExit):
            self.search("g", "Genomes/ref")
        self.assertIn(f"holds the genome in {self.cwd}/projA/ref, not {self.cwd}/Genomes/ref", self.err)

    def test_the_working_directorys_library_given_as_index_path_is_looked_at(self):
        self.search("a", "projA/ref")
        with self.assertRaises(SystemExit):
            self.search("b", "projB/ref", 1, "--index-path", "genome_library")
        self.assertIn("20-NRG_2_ref holds the genome in", self.err)
        # and, the pipeline building no index under --index-path, none is removed either
        _write_genome(os.path.join(self.cwd, "projA", "ref"), self.SEQ["maternal"],
                      mtime=time.time() - 3600)
        self.search("a2", "projA/ref", 1, "--index-path", "genome_library")
        self.assertEqual(self.pipeline.built, ["20-NRG_2_ref"])

    def test_another_index_path_is_not_looked_at(self):
        # no search builds there, and its records name folders of another working directory
        shared = os.path.join(self.cwd, "shared_library")
        os.makedirs(os.path.join(shared, "NRG_2_ref"))
        shutil.copy(os.path.join(self.cwd, "projB", "ref", "chr1.fa"),
                    os.path.join(shared, "NRG_2_ref", "NRG_chr1_1.bin"))
        with open(os.path.join(shared, "NRG_2_ref", ar.INDEX_SOURCE_FILENAME), "w") as fh:
            json.dump({"genome": "projB/ref", "fingerprint": "0" * 40}, fh)
        self.assertEqual(self.search("b", "projA/ref", 1, "--index-path", shared), self.SEQ["maternal"])

    # Each refusal says how to run the search again. Given no --bDNA/--bRNA, a search
    # takes the bulges the indexes of its folder's name allow: refused with the 1 + 1
    # that an index in the way gave it, the same command run again once that index was
    # removed, or the folder renamed, took none. Under --index-path the pipeline builds
    # no index, so the same command run again once the index was removed failed; under
    # --prune-pam, which complete-search refuses without --index-path, the command run
    # again without --index-path alone stopped at once.

    def test_followed_as_it_says_a_refused_search_keeps_its_bulges(self):
        self.search("a", "projA/ref")
        with self.assertRaises(SystemExit):
            self.search("b", "projB/ref", None)  # 1 + 1, from projA's 20-NRG_2_ref
        self.assertIn(f"If {self.cwd}/projA/ref is this same folder under another path, is no "
                      f"longer needed, or gets another name, remove {self.lib}/20-NRG_2_ref when no "
                      f"search is using it; otherwise give this folder another name. Then run "
                      f"this search again with --bDNA 1 --bRNA 1 (given no --bDNA/--bRNA, this "
                      f"search may then take no bulges): it builds the index of its name from "
                      f"this folder.", self.err)
        advice = self.err
        shutil.rmtree(os.path.join(self.lib, "20-NRG_2_ref"))
        self.assertEqual(self.follow(advice, "b2", "projB/ref", None), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.bulges["b2"], (1, 1))
        self.assertEqual(self.record("20-NRG_2_ref")["genome"], os.path.join("projB", "ref"))

    def test_under_index_path_it_says_to_run_again_without_it(self):
        self.search("a", "projA/ref")
        ip = ("--index-path", "genome_library")
        with self.assertRaises(SystemExit):
            self.search("b", "projB/ref", 1, *ip)
        self.assertIn("Then run this search again without --index-path and with --bDNA 1 "
                      "--bRNA 1 (given no --bDNA/--bRNA, this search may then take no bulges): "
                      "it builds the index of its name from this folder.", self.err)
        advice = self.err
        shutil.rmtree(os.path.join(self.lib, "20-NRG_2_ref"))
        with self.assertRaises(SystemExit) as stop:  # the same command: none there, none built
            self.search("b2", "projB/ref", 1, *ip)
        self.assertEqual(stop.exception.code, 1)
        with open(os.path.join(self.cwd, "Results", "b2", "log_error.txt")) as fh:
            self.assertIn("ERROR: no matching reference index under --index-path", fh.read())
        self.assertEqual(self.follow(advice, "b3", "projB/ref", 1, *ip), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.bulges["b3"], (1, 1))

    def test_under_prune_pam_it_says_to_run_again_without_index_path_alone(self):
        # without --index-path, the search builds the index of its name for the PAM
        # itself, which --prune-pam leaves as it is (it said to drop --prune-pam too,
        # which needed --index-path before 2.9.0)
        self.search("a", "projA/ref")
        ip = ("--index-path", "genome_library", "--prune-pam")
        with self.assertRaises(SystemExit):
            self.search("b", "projB/ref", 1, *ip)
        # (--prune-pam no longer needs --index-path, since 2.9.0, and stays)
        self.assertIn("Then run this search again without --index-path and with --bDNA 1 "
                      "--bRNA 1 (given no --bDNA/--bRNA, this search may then take no bulges): "
                      "it builds the index of its name from this folder.", self.err)
        advice = self.err
        shutil.rmtree(os.path.join(self.lib, "20-NRG_2_ref"))
        self.assertEqual(self.follow(advice, "b3", "projB/ref", 1, *ip), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.bulges["b3"], (1, 1))
        self.assertEqual(self.record("20-NRG_2_ref")["genome"], os.path.join("projB", "ref"))

    def test_followed_as_it_says_a_renamed_folder_under_prune_pam_runs(self):
        os.makedirs(self.lib)
        ip = ("--index-path", "genome_library", "--prune-pam")
        with self.assertRaises(SystemExit):
            self.search("i", "other/hg38", 1, *ip)
        self.assertIn(f"otherwise give this one a name no folder in {self.cwd}/Genomes has, and "
                      f"run this search again without --index-path and with --bDNA 1 --bRNA 1 "
                      f"(given no --bDNA/--bRNA, this search may then take no bulges).", self.err)
        advice = self.err
        os.rename(os.path.join(self.cwd, "other", "hg38"), os.path.join(self.cwd, "other", "hg38x"))
        self.assertEqual(self.follow(advice, "i2", "other/hg38x", 1, *ip), self.SEQ["paternal"])
        self.assertEqual(self.pipeline.bulges["i2"], (1, 1))

    def test_followed_as_it_says_a_renamed_folder_keeps_its_bulges(self):
        # under --index-path, it is to run without it: no index of the new name is there
        os.makedirs(self.lib)
        with self.assertRaises(SystemExit):
            self.search("i", "other/hg38", 1, "--index-path", "genome_library")
        self.assertIn(f"otherwise give this one a name no folder in {self.cwd}/Genomes has, and "
                      f"run this search again without --index-path and with --bDNA 1 --bRNA 1 "
                      f"(given no --bDNA/--bRNA, this search may then take no bulges).", self.err)
        with self.assertRaises(SystemExit):
            self.search("o", "other/hg38", None)  # 1 + 1, Genomes/hg38 being there
        self.assertIn(f"otherwise give this one a name no folder in {self.cwd}/Genomes has, and "
                      f"run this search again with --bDNA 1 --bRNA 1 (given no --bDNA/--bRNA, "
                      f"this search may then take no bulges).", self.err)
        advice = self.err
        os.rename(os.path.join(self.cwd, "other", "hg38"), os.path.join(self.cwd, "other", "hg38x"))
        self.assertEqual(self.follow(advice, "o2", "other/hg38x", None), self.SEQ["paternal"])
        self.assertEqual(self.pipeline.bulges["o2"], (1, 1))

    def test_an_index_it_cannot_remove_is_named_with_the_bulges_to_run_again_with(self):
        self.search("a", "projA/ref")
        with mock.patch.object(crisprme, "clear_stale_haplotype_indexes",
                               side_effect=OSError("Permission denied")), \
                self.assertRaises(SystemExit):
            self.search("a2", "projA/ref", None)
        self.assertIn("could not be removed (Permission denied). Remove it and run this search "
                      "again with --bDNA 1 --bRNA 1 (given no --bDNA/--bRNA, this search may then "
                      "take no bulges).", self.err)

    def test_a_search_without_bulges_is_given_none_to_run_again_with(self):
        # (a --vcf search without bulges is refused the same way, and takes none again)
        self.search("a", "projA/ref")
        err = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(err):
            crisprme._check_genome_name_is_its_own(os.path.join(self.cwd, "projB", "ref"), "_",
                                                   0, "NRG")
        self.assertIn("otherwise give this folder another name. Then run this search again: it "
                      "builds the index of its name from this folder.", err.getvalue())
        self.assertNotIn("--bDNA", err.getvalue())
        err = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(err):
            crisprme._check_genome_name_is_its_own(os.path.join(self.cwd, "projB", "ref"),
                                                   self.lib, 0, "NRG", True)
        self.assertIn("Then run this search again without --index-path: it builds the index "
                      "of its name from this folder.", err.getvalue())
        self.assertNotIn("--bDNA", err.getvalue())

    # One haplotype of a bundle searched by its own folder. Before 2.8.1 assembly-search
    # searched every haplotype under its folder's name, "genome", all on the index the
    # first one built, NRG_2_genome, which records no folder and is left in every data
    # folder that ran a bundle search with bulges.

    def bundle(self):
        _make_bundle(self.cwd, "W1", {hap: self.SEQ[hap] for hap in pa.HAPLOTYPES})
        return {hap: os.path.join("Assemblies", "W1", hap, "genome") for hap in pa.HAPLOTYPES}

    def old_index(self, name="NRG_2_genome"):
        """The index of the folder name as assembly-search left it before 2.8.1: built
        from the paternal haplotype, the first searched, and recording nothing."""
        index = os.path.join(self.lib, name)
        os.makedirs(index)
        shutil.copy(os.path.join(self.cwd, "Assemblies", "W1", "paternal", "genome", "chr1.fa"),
                    os.path.join(index, "NRG_chr1_1.bin"))
        return index

    def test_a_haplotype_is_refused_the_index_assembly_search_shared(self):
        haps = self.bundle()
        self.old_index()
        # without bulge flags the search takes 1 + 1 bulges, an index being there; under
        # --index-path the pipeline builds no index (--prune-pam needs it no more)
        for bulges, extra, again in ((1, (), "again with"), (None, (), "again with"),
                                     (1, ("--index-path", "genome_library"),
                                      "again without --index-path and with"),
                                     (1, ("--index-path", "genome_library", "--prune-pam"),
                                      "again without --index-path and with")):
            with self.assertRaises(SystemExit):
                self.search("m", haps["maternal"], bulges, *extra)
            self.assertIn(f"--genome {self.cwd}/{haps['maternal']} is the maternal haplotype of "
                          f"the personal assembly W1, and {self.lib}/NRG_2_genome records no "
                          f"genome folder", self.err)
            # given no --bDNA/--bRNA, assembly-search takes no bulges
            self.assertIn("Search the individual with assembly-search --bDNA 1 --bRNA 1, which "
                          "searches each haplotype on an index of its own (given no "
                          "--bDNA/--bRNA, assembly-search takes no bulges).", self.err)
            self.assertIn(f"remove it and run this search {again} --bDNA 1 --bRNA 1 (given no "
                          f"--bDNA/--bRNA, this search may then take no bulges): it then builds "
                          f"the index from {self.cwd}/{haps['maternal']}.", self.err)
        self.assertEqual(self.pipeline.calls, [])
        self.assertEqual(os.listdir(os.path.join(self.cwd, "Results", "m")), [])  # rerun as it is

    def test_followed_as_it_says_the_haplotypes_search_keeps_its_bulges(self):
        # given no --bDNA/--bRNA it takes 1 + 1 bulges, the index being there; the same
        # command run again once it is removed took none (and built NRG_1_genome)
        haps = self.bundle()
        index = self.old_index()
        with self.assertRaises(SystemExit):
            self.search("m", haps["maternal"], None)
        advice = self.err
        shutil.rmtree(index)
        self.assertEqual(self.follow(advice, "m2", haps["maternal"], None), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.bulges["m2"], (1, 1))
        self.assertEqual(self.pipeline.built, ["20-NRG_2_genome"])

    def test_followed_as_it_says_the_haplotypes_search_under_prune_pam_runs(self):
        haps = self.bundle()
        index = self.old_index()
        ip = ("--index-path", "genome_library", "--prune-pam")
        with self.assertRaises(SystemExit):
            self.search("m", haps["maternal"], 1, *ip)
        advice = self.err
        shutil.rmtree(index)
        self.assertEqual(self.follow(advice, "m2", haps["maternal"], 1, *ip), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.bulges["m2"], (1, 1))
        self.assertEqual(self.pipeline.built, ["20-NRG_2_genome"])

    def test_without_bulges_the_haplotype_is_searched(self):
        # (it scans its folder, and takes no index)
        haps = self.bundle()
        self.old_index()
        self.assertEqual(self.search("m0", haps["maternal"], 0), self.SEQ["maternal"])
        genome = os.path.join(self.cwd, haps["maternal"])
        # nor is a --vcf search without bulges, which reads the folder's variant genomes
        self.assertTrue(crisprme._check_genome_name_is_its_own(genome, "_", 0, "NRG"))

    def test_an_index_this_search_cannot_take_is_no_reason_to_refuse_it(self):
        # of no bulges (a search without bulges builds one, and never uses it), or of
        # another PAM: this search builds its own
        haps = self.bundle()
        self.old_index("NRG_1_genome")
        self.old_index("NGG_2_genome")
        self.assertEqual(self.search("m", haps["maternal"]), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.built, ["20-NRG_2_genome"])

    def test_the_link_assembly_search_searches_the_haplotype_through_is_not_refused(self):
        # its own name: an index of it records nothing only while a search of that
        # haplotype runs (another guide on the web, say), or after one was killed
        haps = self.bundle()
        self.old_index()
        view = ar.haplotype_genome_view(self.cwd, os.path.join(self.cwd, haps["maternal"]),
                                        "maternal")
        own = os.path.join(self.lib, f"NRG_2_{os.path.basename(view)}")
        os.makedirs(own)
        shutil.copy(os.path.join(self.cwd, haps["maternal"], "chr1.fa"),
                    os.path.join(own, "NRG_chr1_1.bin"))
        self.assertEqual(self.search("v", view), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.built, [])

    def test_with_no_such_index_the_haplotype_builds_and_takes_its_own(self):
        haps = self.bundle()
        self.assertEqual(self.search("m", haps["maternal"]), self.SEQ["maternal"])
        self.assertEqual(self.record("20-NRG_2_genome")["genome"], haps["maternal"])
        self.assertEqual(self.search("m2", haps["maternal"]), self.SEQ["maternal"])
        self.assertEqual(self.pipeline.built, ["20-NRG_2_genome"])

    def test_the_other_haplotype_is_then_pointed_to_assembly_search(self):
        # a bundle's folders keep their names: renaming one is no way out
        haps = self.bundle()
        self.search("m", haps["maternal"])
        with self.assertRaises(SystemExit):
            self.search("p", haps["paternal"])
        self.assertIn(f"is the paternal haplotype of the personal assembly W1, and "
                      f"{self.lib}/20-NRG_2_genome holds the genome in {self.cwd}/{haps['maternal']}",
                      self.err)
        self.assertIn("Search the individual with assembly-search", self.err)
        self.assertNotIn("another name", self.err)
        self.assertEqual(self.pipeline.calls, ["m"])

    def test_a_folder_named_genome_that_is_no_haplotype_takes_it_as_before(self):
        # which folder an index that records nothing holds is unknown
        _write_genome(os.path.join(self.cwd, "other", "genome"), self.SEQ["maternal"])
        self.bundle()
        self.old_index()
        self.search("o", "other/genome")
        self.assertEqual(self.pipeline.built, [])


# --------------------------------------------------------------------------- #
# A published index carries no record of the folder it was built from
# --------------------------------------------------------------------------- #
class APublishedIndexRecordsNoFolder(_Workdir):
    """An index records the genome folder it was built from, a folder of the machine
    that built it. publish-index put the record in the archive, and download installed
    it: a search of the downloader's Genomes/hg38 with bulges was then refused, and
    told to remove the index, naming a folder that is not there."""

    NAME = "NRG_2_hg38"

    def setUp(self):
        super().setUp()
        _write_genome(os.path.join(self.cwd, "ext", "hg38"), "ACGT")
        self.index = os.path.join(self.cwd, "genome_library", self.NAME)
        os.makedirs(self.index)
        with open(os.path.join(self.index, "NRG_chr1_1.bin"), "w") as fh:
            fh.write("bin")
        # as a search with bulges of ext/hg38 leaves it
        ar.stamp_index_source(self.cwd, os.path.join(self.cwd, "ext", "hg38"), self.index)
        self.assertEqual(sorted(os.listdir(self.index)), [ar.INDEX_SOURCE_FILENAME, "NRG_chr1_1.bin"])

    def archive(self, which):
        """The members of the archive publish-index makes of the index, with
        `shutil.which` answered by `which`, and what it wrote to stderr."""
        tarball = os.path.join(self.cwd, f"{self.NAME}.tar.gz")
        with mock.patch.object(hf.shutil, "which", side_effect=which), \
                contextlib.redirect_stderr(io.StringIO()) as err:
            hf._make_index_tarball(tarball, self.index, self.NAME, self.index + "_INDELS",
                                   {"name": self.NAME})
        with tarfile.open(tarball) as tf:
            return tf.getnames(), err.getvalue()

    def assertNoRecord(self, names):
        self.assertIn(f"{self.NAME}/NRG_chr1_1.bin", names)
        self.assertIn("manifest.json", names)
        self.assertEqual([n for n in names if os.path.basename(n) == ar.INDEX_SOURCE_FILENAME], [])

    def test_the_name_is_the_records(self):
        self.assertEqual(hf._INDEX_SOURCE_FILENAME, ar.INDEX_SOURCE_FILENAME)

    def test_the_archive_made_without_pigz_has_no_record(self):
        names, said = self.archive(lambda program: None)
        self.assertIn("pigz/tar not found", said)  # the tarfile path
        self.assertNoRecord(names)

    @unittest.skipUnless(shutil.which("tar") and shutil.which("gzip"), "needs tar and gzip")
    def test_the_archive_made_with_tar_and_pigz_has_no_record(self):
        # gzip stands in for pigz: tar runs either as the compressor the same way
        pigz = os.path.join(self.cwd, "bin", "pigz")
        os.makedirs(os.path.dirname(pigz))
        with open(pigz, "w") as fh:
            fh.write('#!/bin/sh\nexec gzip "$@"\n')
        os.chmod(pigz, 0o755)
        which = shutil.which
        names, said = self.archive(lambda program: pigz if program == "pigz" else which(program))
        self.assertEqual(said, "")  # the tar path, with no fallback
        self.assertNoRecord(names)

    def test_a_downloaded_index_records_no_folder(self):
        # an archive that carries the record, as publish-index made it before
        tarball = os.path.join(self.cwd, f"{self.NAME}.tar.gz")
        with tarfile.open(tarball, "w:gz") as tf:
            tf.add(self.index, arcname=self.NAME)
        workdir = os.path.join(self.cwd, "downloader")
        _write_genome(os.path.join(workdir, "Genomes", "hg38"), "ACGT")

        def fetch(repo, allow_patterns, local_dir, token=None, reporter=None, **_kw):
            os.makedirs(os.path.join(local_dir, "indexes"), exist_ok=True)
            shutil.copy(tarball, os.path.join(local_dir, "indexes", f"{self.NAME}.tar.gz"))
            return local_dir

        with mock.patch.object(hf, "_hf_snapshot", side_effect=fetch), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            installed = hf.download_component("index", workdir, repo="local/test",
                                              index_name=self.NAME, token="x")
        self.assertEqual(os.listdir(installed), ["NRG_chr1_1.bin"])
        # so a search of the downloader's own genome of that name takes it, as any download
        self.assertEqual(ar.foreign_index_sources(workdir, os.path.join(workdir, "Genomes", "hg38")), [])


# --------------------------------------------------------------------------- #
# Settings > Delete: an individual goes with what its searches left
# --------------------------------------------------------------------------- #
_SETTINGS_PAGE = os.path.join(os.path.dirname(HERE), "pages", "settings_page.py")


def _settings_functions(names, **namespace):
    """The named functions of pages/settings_page.py, extracted from its source and
    exec'd in `namespace` (it imports a live Dash app, as test_settings_progress
    explains)."""
    with open(_SETTINGS_PAGE) as fh:
        tree = ast.parse(fh.read())
    body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace.update(os=os, List=typing.List)
    exec(compile(ast.Module(body=body, type_ignores=[]), _SETTINGS_PAGE, "exec"), namespace)
    return namespace


class DeletingAnIndividualRemovesItsHaplotypeIndexes(_Workdir):
    """Each searched haplotype has its own index, the size of a reference index, and
    a link it is searched through. Settings > Delete > Personal assembly removed only
    Assemblies/<individual>/: both indexes stayed in Installed indexes, two more for
    every individual searched and deleted, and the links pointed nowhere."""

    def setUp(self):
        super().setUp()
        _make_bundle(self.cwd, "W1", {"paternal": "ACGT", "maternal": "TTTT"})
        _make_bundle(self.cwd, "W2", {"paternal": "GGGG", "maternal": "CCCC"})
        self.lib = os.path.join(self.cwd, "genome_library")
        self.built = {}
        for ind in ("W1", "W2"):
            for hap in pa.HAPLOTYPES:
                view = ar.haplotype_genome_view(self.cwd, self.genome(ind, hap), hap)
                for budget in ("NRG_2", "NGG_3"):
                    index = os.path.join(self.lib, f"{budget}_{os.path.basename(view)}")
                    os.makedirs(index)
                    with open(os.path.join(index, "NRG_chr1_1.bin"), "w") as fh:
                        fh.write("x" * 1000)
                    self.built.setdefault(ind, []).append(index)
                self.built[ind].append(os.path.dirname(view))
        self.others = [os.path.join(self.lib, n) for n in ("NRG_2_hg38", "NRG_2_genome", "NRG_2_W1_paternal")]
        for index in self.others:
            os.makedirs(index)

    def test_its_indexes_and_links_are_found(self):
        self.assertEqual(sorted(ar.bundle_search_leftovers(self.cwd, "W1")), sorted(self.built["W1"]))
        self.assertEqual(ar.bundle_search_leftovers(self.cwd, "not_there"), [])

    def settings(self):
        return _settings_functions(
            ("_delete_targets", "_delete_summary", "_dir_size", "_fmt_size"),
            current_working_directory=self.cwd + os.sep, personal_assembly=pa,
            assembly_reconcile=ar, installed_assemblies=lambda: [],
            GENOMES_DIR="Genomes", LIFTOVER_DIR="LiftoverFiles")

    def test_settings_deletes_them_with_the_individual(self):
        ns = self.settings()
        targets = ns["_delete_targets"]("assembly", "W1")
        self.assertEqual(targets[0], os.path.join(self.cwd, "Assemblies", "W1"))
        self.assertEqual(sorted(targets[1:]), sorted(self.built["W1"]))
        summary = ns["_delete_summary"]("assembly:W1")
        for index in self.built["W1"]:
            if os.path.dirname(index) == self.lib:
                self.assertIn(os.path.basename(index), summary)
        size = ns["_dir_size"](targets[0]) + 4 * 1000  # the four indexes counted
        self.assertIn(f"frees about {ns['_fmt_size'](size)}", summary)
        for path in targets:  # as do_delete removes them
            shutil.rmtree(path)
        self.assertEqual(sorted(os.listdir(self.lib)),
                         sorted(os.path.basename(i) for i in self.built["W2"] + self.others
                                if os.path.dirname(i) == self.lib))
        self.assertTrue(all(os.path.exists(p) for p in self.built["W2"] + self.others))

    def test_an_individual_never_searched_has_nothing_more(self):
        _make_bundle(self.cwd, "W3", {"paternal": "ACGT", "maternal": "TTTT"})
        ns = self.settings()
        self.assertEqual(ns["_delete_targets"]("assembly", "W3"),
                         [os.path.join(self.cwd, "Assemblies", "W3")])
        self.assertNotIn("indexes its searches built", ns["_delete_summary"]("assembly:W3"))


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
    """`_resolve_index_dir`, extracted from the shipped script (not copied) with the
    helpers it calls, which the script defines just before it, as test_index_resolution.py
    extracts it."""
    with open(SHELL_SCRIPT) as fh:
        text = fh.read()
    start = text.index("\t_index_side() {")
    end = text.index("\n\t}\n", text.index("\t_resolve_index_dir() {", start)) + len("\n\t}\n")
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
        with open(os.path.join(out, "log_error.txt"), "w") as err, \
                open(os.path.join(out, "log_verbose.txt"), "w") as verbose:
            subprocess.run(["bash", SHELL_SCRIPT, *map(str, args)], cwd=work, env=env,
                           stdout=verbose, stderr=err, timeout=300)
        with open(os.path.join(out, "log_verbose.txt")) as fh:
            self.said = fh.read().splitlines()
        with open(log) as fh:
            return [line.split() for line in fh]

    def test_without_bulges_the_given_folder_is_scanned(self):
        search = [c for c in self.calls(0) if c[0] == "search"]
        self.assertEqual(len(search), 1, search)
        self.assertIn("-r", search[0])
        self.assertEqual(os.path.realpath(search[0][1]), os.path.realpath(self.genome("W1", "paternal")))

    def test_with_bulges_the_index_is_the_haplotypes_own(self):
        calls = self.calls(1)
        # the line complete-search reads to know that this run built the index it records
        self.assertIn(crisprme._BUILT_REFERENCE_INDEX_LINE, self.said)
        self.assertEqual([c[:2] for c in calls if c[0] == "index-genome"], [["index-genome", "W1_paternal"]])
        search = [c for c in calls if c[0] == "search"]
        self.assertEqual(os.path.basename(search[0][1].rstrip("/")), "20-NRG_2_W1_paternal")

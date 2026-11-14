#!/usr/bin/env python
"""Tests for indel_window: the fake-indel window read from the data.

CRISPRitz 2.10.0 writes 30 reference bases on either side of an indel, 2.9.1 and earlier 25
before and 27 after; the indel log's CHR field (<chrom>_<start>-<end>_<id>) holds each
window's start and end, so a reader takes each window's flanks from its own row and never
assumes a width. A target of L bases (guide and PAM) with b DNA bulges holds a record's
change only where L + b <= LEFT + 2 and L + b <= RIGHT + 1 (LEFT + 1 at an MNV or complex
record).

Run: python -m unittest test_indel_window -v
"""

import fnmatch
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

import index_fit
import indel_window as iw


def row(chrom, pos, ref, alt, left, right, k=1):
    """An indel log row as CRISPRitz's add-variants writes it, of a window of left/right."""
    return "%s_%d-%d_%d\tS1\trs%d\t0.1\t%s_%d_%s_%s\t0,%d\tN" % (
        chrom, pos - 1 - left, pos - 1 + len(ref) + right, k, k, chrom, pos, ref, alt,
        left + len(alt) + right)


# ---- a command CRISPRme prints, or a block of its docs, pasted into bash and into zsh -------
# (wave 1c-b of 2.10.0: zsh, the default shell of macOS, expanded the index patterns of the
# advice's rm -rf itself, and stopped it whole where one matched nothing. Imported by
# test_build_index_only_atomic, test_mega_rebuild and test_indel_empty_chromosome.)

PASTE_SHELLS = ("bash", "zsh")
# what a shell says of a pasted line it cannot run as written
SHELL_ERRORS = ("no matches found", "command not found", "event not found", "parse error",
                "unbound variable", "bad pattern", "unmatched", "syntax error",
                "no such file or directory: #")


def pasted(test, shell, text, cwd, env=None, timeout=600):
    """``text`` pasted into an interactive ``shell`` ("bash" or "zsh") with its default
    options, no startup file of the user's (bash --norc --noprofile -i, zsh -f -i), its lines
    read on its standard input as a paste reaches it: zsh reads a '#' as a word, not a comment
    (INTERACTIVE_COMMENTS is off), and stops a command at a pattern that matches nothing
    (NOMATCH), where bash passes the pattern on as it is; either expands a '!' (history). zsh
    is not on every Linux host: the unit-tests job installs it, and without it the zsh half
    of a test is skipped ("zsh is not installed", test.skipTest: in a subTest, that one
    alone). Returns (the exit status of the text's last command, None where the shell never
    reached the line after it; the text's stdout; the shell's stderr, its prompts and what
    the text's commands wrote there)."""
    if shutil.which(shell) is None:
        test.skipTest("%s is not installed" % shell)
    argv = {"bash": [shell, "--norc", "--noprofile", "-i"], "zsh": [shell, "-f", "-i"]}[shell]
    e = {k: v for k, v in (os.environ if env is None else env).items()
         if k not in ("BASH_ENV", "ENV", "ZDOTDIR", "PROMPT_COMMAND", "TERM_PROGRAM",
                      "TERM_SESSION_ID")}
    # no history file written (an interactive bash writes ~/.bash_history), and macOS's
    # notes off (its bash's deprecation banner, its Terminal's session files)
    e.update(HISTFILE="", SHELL_SESSIONS_DISABLE="1", BASH_SILENCE_DEPRECATION_WARNING="1")
    mark = "__pasted_status__"
    done = subprocess.run(argv, input=text.rstrip("\n") + '\necho "%s=$?"\n' % mark, cwd=cwd,
                          env=e, capture_output=True, text=True, timeout=timeout)
    m = re.search(r"^%s=(\d+)$" % mark, done.stdout, re.M)
    return ((int(m.group(1)) if m else None), done.stdout[:m.start()] if m else done.stdout,
            done.stderr)


def shell_errors(err):
    """The lines of a pasted text's stderr in which the shell says it could not run a line as
    written (SHELL_ERRORS)."""
    return [ln for ln in err.splitlines() if any(s in ln.lower() for s in SHELL_ERRORS)]


def tree(root, *same):
    """{path relative to ``root``: ("dir",), ("file", its sha1), ("link", its target, ``root``
    and each of ``same`` in it read as <root>)} of everything under ``root``: what two runs
    left, compared."""
    out = {}
    for d, dirs, files in os.walk(root):
        for n in dirs + files:
            p = os.path.join(d, n)
            rel = os.path.relpath(p, root)
            if os.path.islink(p):
                target = os.readlink(p)
                for prefix in (root,) + same:
                    target = target.replace(prefix, "<root>")
                out[rel] = ("link", target)
            elif os.path.isdir(p):
                out[rel] = ("dir",)
            else:
                with open(p, "rb") as fh:
                    out[rel] = ("file", hashlib.sha1(fh.read()).hexdigest())
    return out


def paste_both(test, command, root, cwd=None, env=None):
    """``command``, as CRISPRme prints it (its paths absolute, under ``root``, or relative to
    ``cwd``, by default ``root``), pasted into zsh, on a copy of ``root`` (the command's paths
    there), then into bash, on ``root`` itself: each exits 0 and says no error of the shell's,
    and the two leave the same tree. Without zsh, bash alone (its subTest skipped). Returns
    bash's (status, stdout, stderr)."""
    cwd = root if cwd is None else cwd
    with tempfile.TemporaryDirectory(dir=os.path.dirname(root)) as tmp:
        copy = os.path.join(tmp, os.path.basename(root))
        shutil.copytree(root, copy, symlinks=True)
        theirs = None
        with test.subTest(shell="zsh", command=command):
            rc, _out, err = pasted(test, "zsh", command.replace(root, copy),
                                   cwd.replace(root, copy), env)
            test.assertEqual((rc, shell_errors(err)), (0, []), err)
            theirs = tree(copy, root)   # (a link copied as it was points into root)
        rc, out, err = pasted(test, "bash", command, cwd, env)
        test.assertEqual((rc, shell_errors(err)), (0, []), (command, err))
        if theirs is not None:
            test.assertEqual(tree(root), theirs, command)
    return rc, out, err


def same_argv(test, command, cwd, programs=("rm", "find", "mv", "crisprme.py")):
    """The programs ``command`` (as CRISPRme prints it) runs and their arguments, pasted into
    bash and into zsh with each of ``programs`` a stub on PATH that records them and does
    nothing else (exit 0): the same in either shell, and the shell says no error. For a
    command whose programs a test cannot run (a download, a build). Returns them, [[program,
    argument, ...], ...] in their order, bash's."""
    stubs = tempfile.mkdtemp(prefix="pasted_stubs_")
    test.addCleanup(shutil.rmtree, stubs, True)
    log = os.path.join(stubs, "argv.jsonl")
    for name in programs:
        path = os.path.join(stubs, name)
        with open(path, "w") as fh:
            fh.write("#!%s\nimport json, os, sys\nwith open(os.environ['STUB_ARGV'], 'a') as fh:\n"
                     "    fh.write(json.dumps([os.path.basename(sys.argv[0])] + sys.argv[1:]) "
                     "+ '\\n')\n" % sys.executable)
        os.chmod(path, 0o755)
    got = {}
    for shell in PASTE_SHELLS:
        with test.subTest(shell=shell, command=command):
            open(log, "w").close()
            env = dict(os.environ, PATH=stubs + os.pathsep + os.environ.get("PATH", ""),
                       STUB_ARGV=log)
            rc, _out, err = pasted(test, shell, command, cwd, env)
            test.assertEqual((rc, shell_errors(err)), (0, []), err)
            with open(log) as fh:
                got[shell] = [json.loads(ln) for ln in fh]
    if len(got) == len(PASTE_SHELLS):
        test.assertEqual(got["zsh"], got["bash"], command)
    return got.get("bash")


class TheRowsWindow(unittest.TestCase):
    def test_flanks_of_either_release(self):
        # chr22_10524814_TTAG_T, a row of CRISPRitz 2.9.1 (25/27) and one of 2.10.0 (30/30)
        self.assertEqual(iw.row_flanks("chr22_10524788-10524844_1", "chr22_10524814_TTAG_T"),
                         (25, 27))
        self.assertEqual(iw.row_flanks("chr22_10524783-10524847_1", "chr22_10524814_TTAG_T"),
                         (30, 30))

    def test_a_chromosome_name_holding_underscores_and_dashes(self):
        chr_field = "chrUn_KI2707-42v1_100-160_5"
        self.assertEqual((iw.window_start(chr_field), iw.window_end(chr_field)), (100, 160))
        self.assertEqual(iw.row_flanks(chr_field, "chrUn_KI2707-42v1_131_A_AT"), (30, 29))

    def test_a_start_past_the_pos_is_no_window(self):
        with self.assertRaises(ValueError):
            iw.left_flank(10, 20)


class ALogsWindows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def write(self, name, rows, gz=True):
        path = os.path.join(self.tmp, name)
        opener = gzip.open if gz else open
        with opener(path, "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
            fh.write("".join(r + "\n" for r in rows))
        return path

    def test_a_log_of_one_width(self):
        path = self.write("logchr1.txt.gz", [row("chr1", 100, "A", "AT", 30, 30, 1),
                                             row("chr1", 200, "ACG", "A", 30, 30, 2),
                                             row("chr1", 300, "AC", "TG", 30, 30, 3)])
        self.assertEqual(iw.log_windows(path), {(30, 30)})
        self.assertEqual(iw.log_path(self.tmp, "chr1"), path)

    def test_a_plain_log_a_header_alone_and_no_log(self):
        self.write("logchr2.txt", [row("chr2", 100, "A", "AT", 25, 27)], gz=False)
        self.write("logchr3.txt.gz", [])
        self.assertEqual(iw.dataset_windows(self.tmp),
                         {"chr2": {(25, 27)}, "chr3": set()})
        self.assertEqual(iw.dataset_windows(self.tmp, ["chr2", "chr9"]),
                         {"chr2": {(25, 27)}, "chr9": set()})
        self.assertEqual(iw.dataset_window(self.tmp), (25, 27))
        self.assertEqual(iw.log_chroms(self.tmp), ["chr2", "chr3"])

    def test_logs_of_two_enrichments(self):
        self.write("logchr1.txt.gz", [row("chr1", 100, "A", "AT", 30, 30)])
        self.write("logchr2.txt.gz", [row("chr2", 100, "A", "AT", 25, 27)])
        with self.assertRaises(ValueError):
            iw.dataset_window(self.tmp)

    def test_the_first_rows_say_the_width(self):
        rows = [row("chr1", 100 + 10 * k, "A", "AT", 30, 30, k) for k in range(5)]
        path = self.write("logchr1.txt.gz",
                          rows + ["chr1_x-y_9\tS1\t.\t0.1\tchr1_900_A_AT\t0,1\tN"])
        self.assertEqual(iw.log_windows(path, rows=5), {(30, 30)})
        with self.assertRaises(ValueError):
            iw.log_windows(path, rows=None)


class WhatAWindowHolds(unittest.TestCase):
    def test_reach(self):
        self.assertEqual((iw.reach((25, 27)), iw.reach((25, 27), mnv=True)), (27, 26))
        self.assertEqual((iw.reach((30, 30)), iw.reach((30, 30), mnv=True)), (31, 31))

    def test_a_search_that_fits_says_nothing(self):
        # every published index's motif at up to 2 DNA bulges: 29 + 2 on 30/30
        self.assertIsNone(iw.target_note((30, 30), 29, 2, "complete-search", "ds"))
        self.assertIsNone(iw.target_note((25, 27), 23, 3, "complete-search", "ds"))
        self.assertIsNone(iw.target_note(None, 40, 5, "complete-search", "ds"))

    def test_a_search_past_the_window_says_so(self):
        note = iw.target_note((30, 30), 29, 3, "complete-search", "ds")
        self.assertTrue(note.startswith("WARNING [complete-search]: the fake-indel windows of ds "
                                        "hold 30 reference bases before each indel and 30 "
                                        "after it"), note)
        self.assertIn("holds up to 31 bases, guide, PAM and DNA bulges counted; this search's "
                      "targets need up to 32 (a 29-nt guide and PAM with 3 DNA bulges)", note)
        self.assertNotIn("CRISPRitz 2.10.0 writes", note)
        # an archive of 2.9.1's windows: an MNV's limit too, and how to search them
        note = iw.target_note((25, 27), 27, 0, "indel search", "ds")
        self.assertIn("holds up to 27 bases (26 at a record whose REF and ALT differ at their "
                      "first base, an MNV or complex record)", note)
        self.assertIn("this search's targets need up to 27 (a 27-nt guide and PAM)", note)
        self.assertIn("enrich ds again with it (build-index-only)", note)


class _EngineCase(unittest.TestCase):
    """A temporary folder, and engine(): a crispritz.py of the lines given on PATH."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        index_fit._engine_cache.clear()
        self.addCleanup(index_fit._engine_cache.clear)

    def engine(self, lines):
        bindir = os.path.join(self.tmp, "bin")
        os.makedirs(bindir, exist_ok=True)
        exe = os.path.join(bindir, "crispritz.py")
        with open(exe, "w") as fh:
            fh.write("".join(x + "\n" for x in lines))
        os.chmod(exe, 0o755)
        index_fit._engine_cache.clear()
        return mock.patch.dict(os.environ, {"PATH": bindir})


class TheEngineAndTheIndex(_EngineCase):
    def test_the_engine_states_its_window(self):
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertEqual(iw.engine_window(), (30, 30))
        with self.engine(['VERSION = "2.9.1"']):
            self.assertEqual(iw.engine_window(), iw.OLD_WINDOW)
        with mock.patch.dict(os.environ, {"PATH": os.path.join(self.tmp, "none")}):
            index_fit._engine_cache.clear()
            self.assertEqual(iw.engine_window(), iw.WINDOW_2100)

    def test_an_indels_sub_index_says_its_window(self):
        sub = os.path.join(self.tmp, "NRG_2_fakechr1")
        os.makedirs(sub)
        self.assertIsNone(iw.index_window(sub))       # nothing to tell by
        with open(os.path.join(sub, "NRG_fakechr1_1.bin"), "wb") as fh:
            fh.write(b"\x01\x00\x00\x00\x14\x00\x00\x00")   # format 1: 2.9.1's, 25/27
        self.assertEqual(iw.index_window(sub), iw.OLD_WINDOW)
        with open(os.path.join(sub, index_fit.INDEX_META), "w") as fh:
            fh.write('{"index_format": 1, "indel_window": [30, 30]}')
        self.assertEqual((iw.recorded_window(sub), iw.index_window(sub)), ((30, 30), (30, 30)))


class ADatasetEnrichedAnew(_EngineCase):
    """W-R2-F1: what an enrichment left, which a new one of another CRISPRitz must not be
    mixed with: its genomes, dictionaries, indel logs, dict-less tiers and indexes, by
    absolute path (enrichment_stores); the windows its records tell (earlier_windows); and
    the refusal of build-index-only and of a search that would enrich the dataset anew beside
    the stores of an earlier enrichment of other windows (reenrichment_refusal)."""

    def put(self, *rel, text="x", data=None):
        path = os.path.join(self.tmp, "work", *rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb" if data is not None else "w") as fh:
            fh.write(data if data is not None else text)
        return path

    def dataset(self, window=(25, 27), data_type="genotyped-phased", samples=True):
        """A dataset ds of genome hg38 as build-index-only leaves it, its indel log of
        ``window``, its variant index of format 1 (2.9.1's), beside what is not its own; its
        registry's manifest of ``data_type`` (None: none, and the aggregation of a registry
        built from INFO AF fields, as every published mega's), its VCF naming a sample or,
        without ``samples``, none."""
        w = os.path.join(self.tmp, "work")
        shutil.rmtree(w, ignore_errors=True)
        self.put("Genomes", "hg38", "chr1.fa")
        self.put("Genomes", "hg38+ds", "chr1.enriched.fa")
        self.put("Genomes", "hg38+ds_INDELS", "fake_ds_chr1", "fakechr1.fa")
        self.put("Dictionaries", "dictionaries_ds", "my_dict_chr1.json.gz")
        with gzip.open(self.put("Dictionaries", "log_indels_ds", "logchr1.txt.gz", data=b""),
                       "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n"
                     + row("chr1", 100, "A", "AT", *window) + "\n")
        self.put("Dictionaries", "registry_ds", "reg_chr1.idx",
                 text='{"n_records": 2, "data_type": "%s"}' % data_type if data_type else
                 '{"n_records": 2, "aggregation": "info_af", "databases": {}}')
        with gzip.open(self.put("VCFs", "ds", "ds.chr1.vcf.gz", data=b""), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                     + ("\tFORMAT\tS1" if samples else "") + "\n")
        self.put("Dictionaries", "registry_ds", "reg_chr1.bin")
        self.put("Dictionaries", "genotypes_ds", "gt_chr1.idx")
        self.put("Dictionaries", "indel_genotypes_ds", "gt_indel_chr1.tsv.gz")
        fmt1 = b"\x01\x00\x00\x00\x14\x00\x00\x00"
        for name in ("NRG_3_hg38+ds", "20-NRG_2_hg38+ds", "20-NRG_2_hg38",
                     "20-NRG_2_x_hg38+ds", "20-NRG_2_hg38+ds2"):
            self.put("genome_library", name, "NRG_chr1_1.bin", data=fmt1)
        self.put("genome_library", "20-NRG_2_hg38+ds_INDELS", "NRG_2_fakechr1",
                 "NRG_fakechr1_1.bin", data=fmt1)
        os.makedirs(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds.lock"))
        os.makedirs(os.path.join(w, "genome_library", ".build-index-only-x", "NRG_3_hg38+ds"))
        return w

    def test_what_comes_from_an_enrichment(self):
        w = self.dataset()
        g, d, lib = (os.path.join(w, x) for x in ("Genomes", "Dictionaries", "genome_library"))
        indexes = [os.path.join(lib, n) for n in ("20-NRG_2_hg38+ds", "20-NRG_2_hg38+ds_INDELS",
                                                 "NRG_3_hg38+ds")]
        self.assertEqual(iw.enrichment_stores(w, "hg38", "ds"), [
            os.path.join(g, "hg38+ds"), os.path.join(g, "hg38+ds_INDELS"),
            os.path.join(d, "dictionaries_ds"), os.path.join(d, "log_indels_ds"),
            os.path.join(d, "registry_ds"), os.path.join(d, "genotypes_ds")] + indexes)
        # another library a search reads, each index once; the tiers there or not
        other = os.path.join(self.tmp, "staged")
        os.makedirs(os.path.join(other, "NRR_3_hg38+ds"))
        os.symlink(os.path.join(lib, "NRG_3_hg38+ds"), os.path.join(other, "NRG_3_hg38+ds"))
        shutil.rmtree(os.path.join(d, "genotypes_ds"))
        os.makedirs(os.path.join(d, "genotypes_ds", ".building"))
        self.assertEqual(iw.enrichment_stores(w, "hg38", "ds", [other])[4:],
                         [os.path.join(d, "registry_ds")] + indexes
                         + [os.path.join(other, "NRR_3_hg38+ds")])
        self.assertIn(os.path.join(d, "genotypes_ds"),
                      iw.enrichment_stores(w, "hg38", "ds", present=False))
        # a sites-only registry is its VCFs' (the mega's), not the enrichment's
        w = self.dataset(data_type="sites-only")
        self.assertNotIn(os.path.join(d, "registry_ds"), iw.enrichment_stores(w, "hg38", "ds"))
        # W-R3-F1: nor one whose manifest has no data_type, of INFO AF fields, as every
        # published mega's registry
        w = self.dataset(data_type=None)
        self.assertEqual(iw.registry_data_type(os.path.join(d, "registry_ds")), "sites-only")
        self.assertNotIn(os.path.join(d, "registry_ds"), iw.enrichment_stores(w, "hg38", "ds"))
        self.assertNotIn(os.path.join(d, "registry_ds"),
                         iw.enrichment_stores(w, "hg38", "ds", present=False))

    def test_the_windows_an_earlier_enrichment_tells(self):
        w = self.dataset((25, 27))
        lib = os.path.join(w, "genome_library")
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            # each told alone: its log, its INDELS index (format 1: 2.9.1's), its variant index
            os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            shutil.rmtree(os.path.join(lib, "20-NRG_2_hg38+ds_INDELS"))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), {(25, 27)})
            for n in ("NRG_3_hg38+ds", "20-NRG_2_hg38+ds"):
                shutil.rmtree(os.path.join(lib, n))
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), set())
        # a variant index of format 1 tells nothing where the installed CRISPRitz writes it
        w = self.dataset((30, 30))
        os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
        shutil.rmtree(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds_INDELS"))
        with self.engine(['VERSION = "2.9.1"']):
            self.assertEqual(iw.earlier_windows(w, "hg38", "ds"), set())

    def test_a_dataset_enriched_anew_beside_another_crispritzs_stores_is_refused(self):
        w = self.dataset((25, 27))
        g, d, lib = (os.path.join(w, x) for x in ("Genomes", "Dictionaries", "genome_library"))
        shutil.rmtree(os.path.join(g, "hg38+ds"))
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            why = iw.reenrichment_refusal(w, "hg38", "ds")
        built = [os.path.join(d, "registry_ds"), os.path.join(d, "genotypes_ds")] + [
            os.path.join(lib, n) for n in ("20-NRG_2_hg38+ds", "20-NRG_2_hg38+ds_INDELS",
                                           "NRG_3_hg38+ds")]
        self.assertTrue(why.startswith(
            "no finished enrichment of ds is in place (%s), so build-index-only would enrich "
            "it anew with the installed CRISPRitz, whose fake-indel windows are 30/30; but an "
            "earlier enrichment of it, whose records tell windows of 25/27 (CRISPRitz 2.9.1 "
            "and earlier wrote 25/27), left what was built from it: %s." % (
                os.path.join(g, "hg38+ds"), ", ".join(built))), why)
        self.assertIn("(rm -rf %s), then run build-index-only again" % " ".join(
            "'%s'" % p for p in [os.path.join(g, "hg38+ds_INDELS"),
                                 os.path.join(d, "dictionaries_ds"),
                                 os.path.join(d, "log_indels_ds")] + built), why)
        self.assertNotIn("indel_genotypes_ds", why)
        # of the installed CRISPRitz's windows: the same enrichment, nothing refused
        with self.engine(['VERSION = "2.9.1"']):
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))
        # nothing built from it: nothing to refuse
        for p in built:
            shutil.rmtree(p)
        with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                          "INDEL_RIGHT = 30"]):
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))

    def test_the_search_scripts_check(self):
        """indel_window.py reenrichment, as the search script runs it before it enriches a
        dataset itself (STEP 1): the refusal on stdout, exit 3; nothing, exit 0."""
        import subprocess
        import sys

        w = self.dataset((25, 27))
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "indel_window.py")
        for lines, rc in (([ 'VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                            "INDEL_RIGHT = 30"], 3), (['VERSION = "2.9.1"'], 0)):
            with self.engine(lines):
                done = subprocess.run([sys.executable, script, "reenrichment", w, "hg38", "ds",
                                       os.path.join(self.tmp, "staged")],
                                      capture_output=True, text=True, env=dict(os.environ))
            self.assertEqual((done.returncode, done.stderr), (rc, ""), done.stdout)
            self.assertEqual(done.stdout.startswith(
                "no finished enrichment of ds is in place (%s), so this search would enrich "
                "it anew" % os.path.join(w, "Genomes", "hg38+ds")), rc == 3, done.stdout)
        # --vcf: the dataset's VCF folder, which tells a sites-only one (W-R3-F1), beside the
        # windows of the installed CRISPRitz too
        w = self.dataset((30, 30), samples=False)   # (its registry a genotyped dataset's)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        for n in os.listdir(os.path.join(w, "genome_library")):   # (2.9.1's indexes)
            shutil.rmtree(os.path.join(w, "genome_library", n))
        for vcf, rc in ((["--vcf", os.path.join(w, "VCFs", "ds")], 3), ([], 0)):
            with self.engine(['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                              "INDEL_RIGHT = 30"]):
                done = subprocess.run([sys.executable, script, "reenrichment", *vcf, w, "hg38",
                                       "ds"], capture_output=True, text=True,
                                      env=dict(os.environ))
            self.assertEqual((done.returncode, done.stderr), (rc, ""), done.stdout)
            self.assertEqual("so this search would enrich it anew from its VCFs; but ds is "
                             "sites-only" in done.stdout, rc == 3, done.stdout)


    def test_a_registry_says_its_data_type(self):
        """registry_data_type reads a manifest as crisprme.py's _registry_data_type does: its
        data_type, else sites-only for an aggregation of INFO AF fields."""
        reg = os.path.join(self.tmp, "registry_x")
        os.makedirs(reg)
        self.assertIsNone(iw.registry_data_type(reg))
        self.assertIsNone(iw.registry_data_type(os.path.join(self.tmp, "none")))
        for text, kind in (('{"aggregation": "info_af"}', "sites-only"),
                           ('{"aggregation": "carriers"}', None),
                           ('{"aggregation": "info_af", "data_type": "sites-only"}', "sites-only"),
                           ('{"aggregation": "panel", "data_type": "genotyped-phased"}',
                            "genotyped-phased"),
                           ("not json", None), ("[1, 2]", None)):
            with open(os.path.join(reg, "reg_chr1.idx"), "w") as fh:
                fh.write(text)
            self.assertEqual(iw.registry_data_type(reg), kind, text)
        # the first manifest that tells
        with open(os.path.join(reg, "reg_chr2.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        self.assertEqual(iw.registry_data_type(reg), "sites-only")

    def test_a_sites_only_dataset(self):
        """sites_only: its VCFs name no sample, where there are any to read; else its
        registry says so. A VCF of the synthetic MEGA sample (synth_sites_gt.py) names one."""
        w = self.dataset(data_type=None, samples=False)
        vcfs = os.path.join(w, "VCFs", "ds")
        self.assertTrue(iw.sites_only(w, "ds", vcfs))
        self.assertTrue(iw.sites_only(w, "ds"))            # (by its registry)
        w = self.dataset(data_type=None, samples=True)
        self.assertFalse(iw.sites_only(w, "ds", vcfs))     # (the VCFs tell first)
        self.assertTrue(iw.sites_only(w, "ds", os.path.join(w, "VCFs", "none")))
        w = self.dataset(data_type="genotyped-phased", samples=False)
        self.assertTrue(iw.sites_only(w, "ds", vcfs))
        self.assertFalse(iw.sites_only(w, "ds"))

    def test_a_sites_only_dataset_is_not_enriched_anew_beside_its_kits_windows(self):
        """W-R3-F1: the mega's windows come from its kit (a synthetic sample carries its
        indels), and an enrichment of its VCFs, which name no sample, holds none. Enriched
        anew beside the windows of its logs or INDELS indexes, whatever their width, it
        would lose every indel off-target: refused, pointing to the kit, its registry (no
        data_type, as the published mega's) not named, nor build-index-only."""
        engine30 = ['VERSION = "2.10.0"', "INDEX_FORMAT = 2", "INDEL_LEFT = 30",
                    "INDEL_RIGHT = 30"]
        for window in ((25, 27), (30, 30)):
            w = self.dataset(window, data_type=None, samples=False)
            g, d = os.path.join(w, "Genomes"), os.path.join(w, "Dictionaries")
            vcfs = os.path.join(w, "VCFs", "ds")
            for p in (os.path.join(g, "hg38+ds"), os.path.join(g, "hg38+ds_INDELS"),
                      os.path.join(d, "genotypes_ds"), os.path.join(d, "dictionaries_ds")):
                shutil.rmtree(p)
            for n in os.listdir(os.path.join(w, "genome_library")):   # (MG1: logs alone)
                shutil.rmtree(os.path.join(w, "genome_library", n))
            with self.engine(engine30):
                why = iw.reenrichment_refusal(w, "hg38", "ds", vcf_dir=vcfs)
            self.assertTrue(why.startswith(
                "no finished enrichment of ds is in place (%s), so build-index-only would "
                "enrich it anew from its VCFs; but ds is sites-only: its VCFs name no sample, "
                "and CRISPRitz writes a fake-indel window only for an indel a sample carries, "
                "so an enrichment of them holds none, while its indel logs or its INDELS "
                "indexes record windows of %s" % (os.path.join(g, "hg38+ds"),
                                                  iw.describe(window))), why)
            self.assertIn("Put that kit back in place (%s and %s_INDELS: of a published "
                          "dataset, its enriched-genome kit" % ((os.path.join(g, "hg38+ds"),) * 2),
                          why)
            # (the removal of what came from the kit names no registry: it comes from the
            # VCFs alone; the recipe's own replacement of it, IR1-F2, is after)
            self.assertNotIn("registry_ds", why.split("or build its kit again that way")[0])
            self.assertNotIn("build-index-only again", why)
            self.assertEqual("rm -rf" in why, window != (30, 30), why)
            if window != (30, 30):
                self.assertIn("remove all that came from that kit (rm -rf '%s') and install "
                              "the dataset's archive and enriched-genome kit of CRISPRme 2.10.0 "
                              "in its place, or build its kit again that way with "
                              "seq_script/merge_panels/mega_rebuild.sh in a new working "
                              "directory of its own" % os.path.join(d, "log_indels_ds"), why)
                self.assertIn(iw.rebuild_kit_way("ds", w, "hg38"), why)
            else:
                self.assertNotIn("registry_ds", why)
            # by its registry alone, its VCFs not there
            with self.engine(engine30):
                self.assertTrue(iw.reenrichment_refusal(w, "hg38", "ds").startswith(
                    "no finished enrichment of ds is in place"))
        # its INDELS index alone tells it had windows
        w = self.dataset((30, 30), data_type=None, samples=False)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        os.remove(os.path.join(w, "Dictionaries", "log_indels_ds", "logchr1.txt.gz"))
        with self.engine(engine30):
            self.assertIn("is sites-only", iw.reenrichment_refusal(w, "hg38", "ds"))
            # nothing that holds a window, a variant index of format 1 alone (2.9.1's): no
            # window to lose, refused as any dataset beside another CRISPRitz's stores is,
            # its registry not named
            shutil.rmtree(os.path.join(w, "genome_library", "20-NRG_2_hg38+ds_INDELS"))
            why = iw.reenrichment_refusal(w, "hg38", "ds")
            self.assertTrue(why.startswith(
                "no finished enrichment of ds is in place (%s), so build-index-only would "
                "enrich it anew with the installed CRISPRitz" % os.path.join(
                    w, "Genomes", "hg38+ds")), why)
            self.assertNotIn("registry_ds", why)
            for n in os.listdir(os.path.join(w, "genome_library")):
                shutil.rmtree(os.path.join(w, "genome_library", n))
            self.assertIsNone(iw.reenrichment_refusal(w, "hg38", "ds"))
        # VCFs that name a sample (the synthetic MEGA sample) write windows: enriched anew as
        # any dataset beside the stores of the installed CRISPRitz's windows; by its registry
        # alone (its VCFs not given), refused
        w = self.dataset((30, 30), data_type=None, samples=True)
        shutil.rmtree(os.path.join(w, "Genomes", "hg38+ds"))
        for n in os.listdir(os.path.join(w, "genome_library")):
            shutil.rmtree(os.path.join(w, "genome_library", n))
        with self.engine(engine30):
            self.assertIsNone(iw.reenrichment_refusal(
                w, "hg38", "ds", vcf_dir=os.path.join(w, "VCFs", "ds")))
            self.assertIn("is sites-only", iw.reenrichment_refusal(w, "hg38", "ds"))

    def test_the_advice_on_a_sites_only_dataset(self):
        """W-R3-F1: a search's window note and the fake-genome refusal point a sites-only
        dataset to its archive and kit, or its kit built again, not to build-index-only."""
        note = iw.target_note((25, 27), 27, 1, "complete-search", "hg38_mega", sites=True)
        self.assertIn("CRISPRitz 2.10.0 writes windows of 30 bases on either side, which hold "
                      "up to 31. To search them: hg38_mega is sites-only: its VCFs name no "
                      "sample, and CRISPRitz writes a fake-indel window only for an indel a "
                      "sample carries, so an enrichment of them holds none (so build-index-only "
                      "would leave it with no indel off-target); its windows come from its "
                      "enriched-genome kit, whose indels a synthetic sample carries: install its "
                      "archive and enriched-genome kit of CRISPRme 2.10.0 "
                      "(docs/PRECOMPUTED_INDEXES.md), or build its kit again that way with "
                      "seq_script/merge_panels/mega_rebuild.sh in a new working directory of its "
                      "own, whose Genomes/<genome> holds the FASTAs of the chromosomes to build "
                      "and no other", note)
        self.assertTrue(note.endswith(iw.rebuild_kit_way("hg38_mega") + ". " + iw.kit_removal(
            "its kit in place", "hg38_mega")), note)
        self.assertNotIn("enrich hg38_mega again", note)
        self.assertIn("enrich ds again with it (build-index-only) to search them.",
                      iw.target_note((25, 27), 27, 1, "complete-search", "ds"))


class TheKitBuiltAgainGoesInANewWorkingDirectory(unittest.TestCase):
    """IR1-F2, the integration review r1 of 2.10.0: the advice on a sites-only dataset said to
    build its kit again with seq_script/merge_panels/mega_rebuild.sh, in place of what its rm
    named, which leaves the dataset's VCFs, registry and indel frequencies (they come from its
    VCFs alone); the recipe, run in that folder with the merged VCFs as its source, refused it
    ("... VCFs/hg38_mega already holds VCFs"). The advice gives the flow that works: the recipe
    in a new working directory of its own (it refuses one that holds a part it builds), then
    the dataset's VCFs, registry and indel frequencies removed too, by absolute path, and each
    part the recipe built moved in; and docs/TROUBLESHOOTING.md section 12, which it names, has
    those commands (test_mega_rebuild runs them)."""

    PARTS = ("Genomes/hg38+hg38_mega and Genomes/hg38+hg38_mega_INDELS",
             "the indexes of hg38_mega of its genome_library (not its reference index)",
             "Dictionaries/log_indels_hg38_mega, registry_hg38_mega and indel_af_hg38_mega, "
             "and VCFs/hg38_mega")

    def check(self, said, workdir):
        self.assertIn("build its kit again that way with seq_script/merge_panels/mega_rebuild.sh"
                      " in a new working directory of its own, whose Genomes/hg38 holds the "
                      "FASTAs of the chromosomes to build and no other, beside the PAM file (the "
                      "recipe refuses one that holds a part of hg38_mega it builds, its VCFs "
                      "among them, and makes nothing there), from the merged VCFs "
                      "mega.<chrom>.afmax.vcf.gz "
                      "(download --what vcf --dataset hg38_mega)", said)
        three = ("VCFs/hg38_mega", "Dictionaries/registry_hg38_mega",
                 "Dictionaries/indel_af_hg38_mega")
        self.assertIn("remove there the dataset's VCFs, registry and indel frequencies too, "
                      "which the recipe wrote anew (rm -rf %s)" % " ".join(
                          iw.quote(os.path.join(workdir, p)) if workdir else p for p in three),
                      said)
        for part in self.PARTS:
            self.assertIn(part, said)
        self.assertIn("(docs/TROUBLESHOOTING.md section 12 has the commands)", said)

    def test_each_message_gives_the_flow(self):
        self.check("or " + iw.sites_only_way("hg38_mega", "/data/w", "hg38"), "/data/w")
        note = iw.target_note((25, 27), 27, 1, "complete-search", "hg38_mega", sites=True)
        self.check(note.replace("<genome>", "hg38"), None)

    def test_a_searchs_note_names_the_removal_of_the_kit_in_place(self):
        """W1C-B1, the review of wave 1c-a: a search's note on a sites-only dataset named only
        the removal of its VCFs, registry and indel frequencies, which the recipe writes anew
        (rebuild_kit_way's, relative to the working directory), and none of the kit in place:
        followed with section 12, the kit stayed, every part the recipe built failed to move
        in, and build-index-only and every search then refused the dataset. As
        build-index-only's warning does, the note names that removal now (kit_removal): the
        kit's genomes and indel logs, each by absolute path, and every index of the dataset,
        in the working directory's genome_library and in a library the search reads besides,
        by a find of their name patterns (index_removals); not its registry, its indel
        frequencies or its VCFs (they come from its VCFs alone), nor the reference index.
        W1C-B2R2-2, the review of the fix: the note named the indexes in place when it was
        printed, and the search then built its own indexes of the dataset from the old kit (a
        SaCas9 search, beside 2.8.1's mega), which the note's own flow left in place, and a
        later search read: pasted into bash or into zsh, its removal takes those built after
        the note too, in either library, and not an index of another genome or dataset, the
        reference index or a hidden folder."""
        top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, top, True)
        w, lib = os.path.join(top, "work"), os.path.join(top, "lib")
        kit = [("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
               ("Dictionaries", "log_indels_hg38_mega")]
        indexes = [("genome_library", "NGG_2_hg38+hg38_mega"),
                   ("genome_library", "20-NRR_2_hg38+hg38_mega_INDELS")]
        kept = [("Dictionaries", "registry_hg38_mega"), ("VCFs", "hg38_mega"),
                ("Dictionaries", "indel_af_hg38_mega"), ("genome_library", "20-NRR_2_hg38")]
        for part in kit + indexes + kept:
            os.makedirs(os.path.join(w, *part))
        with open(os.path.join(w, "Dictionaries", "registry_hg38_mega", "reg_chr1.idx"), "w") as fh:
            fh.write('{"n_records": 1, "aggregation": "info_af", "databases": {}}')
        os.makedirs(os.path.join(lib, "20-NRR_2_hg38+hg38_mega"))
        # (built after the note, by its search: the indexes of its PAM, from the kit in place)
        later = [os.path.join(w, "genome_library", "21-NNGRRT_2_hg38+hg38_mega"),
                 os.path.join(w, "genome_library", "21-NNGRRT_2_hg38+hg38_mega_INDELS"),
                 os.path.join(lib, "21-NNGRRT_2_hg38+hg38_mega_INDELS")]
        others = [os.path.join(w, "genome_library", n) for n in (
            "21-NNGRRT_2_hg38", "20-NRR_2_x_hg38+hg38_mega", "20-NRR_2_hg38+hg38_HGDP1kGP",
            "20-NRR_2_hg38+hg38_mega.lock", ".20-NRR_2_hg38+hg38_mega")]
        for who in ("complete-search", "indel search"):
            note = iw.target_note((25, 27), 27, 1, who, "hg38_mega", sites=True, workdir=w,
                                  genome="hg38", libraries=[lib])
            self.check(note, w)
            m = re.search(r"\. Either replaces its kit in place and all that was built from it, "
                          r"as a search would read stores of the two together: its genomes, "
                          r"dictionaries, indel logs and dict-less tiers, and every index of the "
                          r"dataset, of every PAM, those a search builds from it after this note "
                          r"among them, which go first \((.*?); and any copy of those "
                          r"indexes in a library you search with --index-path\)\.$", note)
            self.assertIsNotNone(m, note)
            self.assertEqual(m.group(1).split("; "), [
                "rm -rf " + " ".join(iw.quote(os.path.join(w, *p)) for p in kit)]
                + iw.index_removals(os.path.join(w, "genome_library"), "hg38", "hg38_mega")
                + iw.index_removals(lib, "hg38", "hg38_mega"))
            self.assertIn("; find %s -mindepth 1 -maxdepth 1 -name '*_[0-9]_hg38+hg38_mega_INDELS' "
                          "-not -name '.*' -exec rm -rf -- {} +" % lib, m.group(1))
            for part in kept:
                self.assertNotIn(iw.quote(os.path.join(w, *part)), m.group(1))
            if who != "complete-search":
                continue
            # the removal pasted into zsh and into bash, after the search built its own indexes
            for d in later + others:
                os.makedirs(d)
            paste_both(self, m.group(1), top)
            for part in kit + indexes:
                self.assertFalse(os.path.lexists(os.path.join(w, *part)), part)
            for d in later + [os.path.join(lib, "20-NRR_2_hg38+hg38_mega")]:
                self.assertFalse(os.path.lexists(d), d)
            for d in others + [os.path.join(w, *p) for p in kept]:
                self.assertTrue(os.path.isdir(d), d)
            for part in kit + indexes:
                os.makedirs(os.path.join(w, *part))
            os.makedirs(os.path.join(lib, "20-NRR_2_hg38+hg38_mega"))

    def test_the_notes_removal_does_the_same_pasted_into_bash_or_zsh(self):
        """Wave 1c-b of 2.10.0: the note gave the dataset's indexes to rm -rf as shell
        patterns, which zsh, the default shell of macOS, expands itself: where one matched
        nothing (no index of the dataset of one kind, or in one library), zsh stopped the
        whole command (no matches found) and removed nothing, the kit in place among it,
        where bash passed the pattern on to rm as it was and removed the rest. A find of each
        pattern (index_removals) does the same in either shell: here with an index of each
        kind in both libraries (matched), with none (unmatched), and with a variant index in
        one and an INDELS index in the other (partly matched), the kit's genomes and indel
        logs and every index of the dataset go, in both shells alike, and nothing else."""
        states = {
            "matched": [("work", "genome_library", "NGG_2_hg38+hg38_mega"),
                        ("work", "genome_library", "20-NRR_2_hg38+hg38_mega_INDELS"),
                        ("lib", "20-NRR_2_hg38+hg38_mega"),
                        ("lib", "21-NNGRRT_2_hg38+hg38_mega_INDELS")],
            "unmatched": [],
            "partly matched": [("work", "genome_library", "NGG_2_hg38+hg38_mega"),
                               ("lib", "21-NNGRRT_2_hg38+hg38_mega_INDELS")]}
        kit = [("work", "Genomes", "hg38+hg38_mega"), ("work", "Genomes", "hg38+hg38_mega_INDELS"),
               ("work", "Dictionaries", "log_indels_hg38_mega")]
        kept = [("work", "Dictionaries", "registry_hg38_mega"), ("work", "VCFs", "hg38_mega"),
                ("work", "Dictionaries", "indel_af_hg38_mega"),
                ("work", "genome_library", "20-NRR_2_hg38"), ("lib", "20-NRR_2_hg38"),
                ("work", "genome_library", "20-NRR_2_x_hg38+hg38_mega"),
                ("work", "genome_library", "20-NRR_2_hg38+hg38_HGDP1kGP"),
                ("work", "genome_library", "20-NRR_2_hg38+hg38_mega.lock"),
                ("work", "genome_library", ".20-NRR_2_hg38+hg38_mega"),
                ("lib", ".21-NNGRRT_2_hg38+hg38_mega_INDELS")]
        for state, indexes in sorted(states.items()):
            with self.subTest(state=state):
                top = tempfile.mkdtemp()
                self.addCleanup(shutil.rmtree, top, True)
                w, lib = os.path.join(top, "work"), os.path.join(top, "lib")
                for part in kit + indexes + kept:
                    os.makedirs(os.path.join(top, *part))
                with open(os.path.join(w, "Dictionaries", "registry_hg38_mega", "reg_chr1.idx"),
                          "w") as fh:
                    fh.write('{"n_records": 1, "aggregation": "info_af", "databases": {}}')
                note = iw.target_note((25, 27), 27, 1, "complete-search", "hg38_mega",
                                      sites=True, workdir=w, genome="hg38", libraries=[lib])
                m = re.search(r"which go first \((.*?); and any copy of those indexes in a "
                              r"library you search with --index-path\)\.$", note)
                self.assertIsNotNone(m, note)
                paste_both(self, m.group(1), top)
                for part in kit + indexes:
                    self.assertFalse(os.path.lexists(os.path.join(top, *part)), part)
                for part in kept:
                    self.assertTrue(os.path.isdir(os.path.join(top, *part)), part)

    def test_section_12_has_its_commands(self):
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(os.path.dirname(here), "docs", "TROUBLESHOOTING.md")) as fh:
            text = fh.read()
        section = text[text.index("## 12. "):text.index("## 13. ")]
        start = section.index("**Building the mega's kit again.**")
        block = section[section.index("```bash\n", start):]
        block = block[:block.index("```", 8)].replace("\\\n", " ")
        lines = [ln.split("#")[0].strip() for ln in block.splitlines()]
        lines = [" ".join(ln.split()) for ln in lines if ln]
        # the recipe in a new working directory of the mega's chromosomes alone, with the PAM
        self.assertIn('mkdir -p "$NEW/Genomes/hg38" "$NEW/PAMs"', lines)
        self.assertIn('PAM=PAMs/20bp-NRR-SpCas9.txt bash '
                      '"$REPO/seq_script/merge_panels/mega_rebuild.sh" "$NEW" "$SRC" < /dev/null',
                      lines)
        # in the data folder, every part of the dataset the recipe builds, each it refuses in
        # its working directory, removed by name, whatever a message's own rm left: the kit in
        # place and all that was built from it, then the three the recipe writes anew. The
        # removal named the three alone, and a search's note named no other (W1C-B1)
        (rm,) = [ln + " " for ln in lines if ln.startswith("rm -rf ")]
        with open(os.path.join(os.path.dirname(here), "seq_script", "merge_panels",
                               "mega_rebuild.sh")) as fh:
            refused = set(re.findall(r'"\$D/([a-z_]+)_\$DATASET"', fh.read()))
        self.assertEqual(refused, {"dictionaries", "log_indels", "registry", "indel_af",
                                   "genotypes", "indel_genotypes"})
        for part in (['"$DATA/Genomes/hg38+hg38_mega" ', '"$DATA/Genomes/hg38+hg38_mega_INDELS" ',
                      '"$DATA/VCFs/hg38_mega" ']
                     + ['"$DATA/Dictionaries/%s_hg38_mega" ' % k for k in sorted(refused)]):
            self.assertIn(part, rm)
        # and the dataset's indexes of every PAM and either form of the name, by a find of the
        # name patterns the messages' own removal takes (indel_window.index_removals), not by
        # patterns given to rm -rf, which zsh, the default shell of macOS, expanded itself
        # and stopped the whole command where one matched nothing (wave 1c-b of 2.10.0)
        finds = [ln for ln in lines if ln.startswith("find ")]
        self.assertEqual(finds, ["find \"$DATA/genome_library\" -mindepth 1 -maxdepth 1 -name %s "
                                 "-not -name '.*' -exec rm -rf -- {} +" % shlex.quote(p)
                                 for p in iw.index_name_patterns("hg38", "hg38_mega")])
        self.assertFalse([ln for ln in lines if "*" in ln and ln.startswith("rm ")])
        # (the dataset's indexes of every PAM and either form of the name, not an index of
        # another genome whose name ends in _hg38, nor the reference index)
        for name, theirs in (("NGG_2_hg38+hg38_mega", True), ("20-NRR_2_hg38+hg38_mega", True),
                             ("20-NRG_2_x_hg38+hg38_mega", False), ("20-NRR_2_hg38", False)):
            self.assertEqual(fnmatch.fnmatchcase(name, "*_[0-9]_hg38+hg38_mega"), theirs, name)
        # each part the messages name, moved in (and the recipe's reference index not)
        moved = " ".join(ln for ln in lines if ln.startswith(("mv ", "[ -e")))
        for part in ('"$NEW/Genomes/hg38+hg38_mega"', '"$NEW/Genomes/hg38+hg38_mega_INDELS"',
                     '"$NEW"/genome_library/*_hg38+hg38_mega ',
                     '"$NEW"/genome_library/*_hg38+hg38_mega_INDELS',
                     '"$NEW/Dictionaries/log_indels_hg38_mega"',
                     '"$NEW/Dictionaries/registry_hg38_mega"',
                     '"$NEW/Dictionaries/indel_af_hg38_mega"', '"$NEW/VCFs/hg38_mega"',
                     '"$NEW/samplesIDs/hg38_mega.samplesID.txt"'):
            self.assertIn(part, moved)
        self.assertNotIn("genome_library/*_hg38 ", moved)
        # and the section says, as the recipe does, that it makes nothing in a folder it refuses
        self.assertIn("it refuses one that holds the dataset's VCFs, or a part it builds, and "
                      "makes nothing there", " ".join(section.split()))
        # W1C-B2R2-1, the review of wave 1c-a: the commands ran on whatever failed before them,
        # and took the PAM file from the data folder, which one that CRISPRme 2.8.1 set up
        # lacks: the recipe built nothing, and the removal took the dataset whole. They run in
        # bash, which stops at the first command that fails, an unset or empty variable among
        # them, whichever shell they are pasted into (wave 1c-b of 2.10.0: a subshell,
        # ( set -eu ... ), ran in zsh where a user pasted it into zsh), the PAM file the
        # repository's, and the data folder is changed only after the recipe, once every part
        # the moves take is there (test_mega_rebuild pastes them into bash and into zsh)
        self.assertEqual((lines[0], lines[1], lines[2], lines[-1]), (
            "```bash", 'bash -euo pipefail -s "$DATA" "$NEW" "$SRC" "$REPO" <<\'EOF\'',
            "DATA=${1:?DATA is not set} NEW=${2:?NEW is not set} SRC=${3:?SRC is not set} "
            "REPO=${4:?REPO is not set}", "EOF"))
        self.assertFalse([ln for ln in lines[1:] if "!" in ln])   # (no history expansion)
        self.assertIn('cp "$REPO/PAMs/20bp-NRR-SpCas9.txt" "$NEW/PAMs/"', lines)
        self.assertFalse([ln for ln in lines if "$DATA/PAMs" in ln])
        recipe = next(i for i, ln in enumerate(lines) if "mega_rebuild.sh" in ln)
        checked = next(i for i, ln in enumerate(lines) if ln.startswith('[ -d "$NEW/$p" ] ||'))
        changes = [i for i, ln in enumerate(lines) if "$DATA" in ln and
                   ln.startswith(("rm ", "mv ", "[ -e", "mkdir ", "find "))]
        self.assertTrue(changes)
        self.assertLess(recipe, checked)
        self.assertLess(checked, min(changes))
        self.assertIn("stops at the first command that fails", " ".join(section.split()))


class TheIndexesOfADataset(unittest.TestCase):
    """index_name_patterns and index_removals, the find of each name pattern an advice removes
    a dataset's indexes by, whenever they were built; and with_its_variant_index, an INDELS
    index of another enrichment than the logs with the variant index of its name, built with
    it (W1C-B2R2-2, the review of wave 1c-a)."""

    def test_the_finds_take_the_datasets_indexes_and_no_other(self):
        self.assertEqual(iw.index_name_patterns("hg38", "hg38_mega"),
                         ["*_[0-9]_hg38+hg38_mega", "*_[0-9]_hg38+hg38_mega_INDELS"])
        self.assertEqual(iw.index_removals("/data/w/genome_library", "hg38", "hg38_mega"), [
            "find /data/w/genome_library -mindepth 1 -maxdepth 1 -name '*_[0-9]_hg38+hg38_mega' "
            "-not -name '.*' -exec rm -rf -- {} +",
            "find /data/w/genome_library -mindepth 1 -maxdepth 1 -name "
            "'*_[0-9]_hg38+hg38_mega_INDELS' -not -name '.*' -exec rm -rf -- {} +"])
        # a folder or a name the shell reads is quoted, and a character find -name reads (a
        # '*' or a '[' of a dataset's name) escaped, so that it reads as itself
        self.assertEqual(iw.index_name_patterns("hg38", "a*b[1]"),
                         ["*_[0-9]_hg38+a\\*b\\[1]", "*_[0-9]_hg38+a\\*b\\[1]_INDELS"])
        self.assertEqual(iw.index_removals("/data/my lib", "hg38", "a b")[0],
                         "find '/data/my lib' -mindepth 1 -maxdepth 1 -name '*_[0-9]_hg38+a b' "
                         "-not -name '.*' -exec rm -rf -- {} +")
        # pasted into bash and into zsh, with an index of each kind (matched), with none
        # (unmatched), with one kind alone (partly matched): every index of the dataset, of
        # every PAM and either form of the name, and no other folder
        cases = {"hg38_mega": (["NGG_2_hg38+hg38_mega", "20-NRR_2_hg38+hg38_mega",
                                "21-NNGRRT_2_hg38+hg38_mega_INDELS", "TTTV-23_3_hg38+hg38_mega"],
                               ["20-NRR_2_hg38", "20-NRR_2_x_hg38+hg38_mega",
                                "20-NRR_2_hg38+hg38_HGDP1kGP", "20-NRR_2_hg38+hg38_mega.lock",
                                "20-NRR_2_hg38+x_hg38_mega", ".20-NRR_2_hg38+hg38_mega",
                                ".build-index-only-x_2_hg38+hg38_mega_INDELS"]),
                 "a b": (["NGG_2_hg38+a b", "NGG_2_hg38+a b_INDELS"], ["NGG_2_hg38+a bc"]),
                 "a*b": (["NGG_2_hg38+a*b", "NGG_2_hg38+a*b_INDELS"],
                         ["NGG_2_hg38+axyb", "NGG_2_hg38+a*b2"]),
                 "o'k": (["NGG_2_hg38+o'k_INDELS"], ["NGG_2_hg38+ok"])}
        for dataset, (ours, others) in sorted(cases.items()):
            for state, made in (("matched", ours), ("unmatched", []), ("partly matched", ours[:1])):
                with self.subTest(dataset=dataset, state=state):
                    top = tempfile.mkdtemp()
                    self.addCleanup(shutil.rmtree, top, True)
                    lib = os.path.join(top, "my lib")
                    for name in made + others:
                        os.makedirs(os.path.join(lib, name))
                    paste_both(self, "; ".join(iw.index_removals(lib, "hg38", dataset)), top)
                    self.assertEqual(sorted(os.listdir(lib)), sorted(others))

    def test_an_indels_index_goes_with_its_variant_index(self):
        top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, top, True)
        indels = os.path.join(top, "21-NNGRRT_2_hg38+hg38_mega_INDELS")
        os.makedirs(indels)
        self.assertEqual(iw.with_its_variant_index(indels), [indels])
        os.makedirs(indels[: -len("_INDELS")])
        self.assertEqual(iw.with_its_variant_index(indels), [indels, indels[: -len("_INDELS")]])


STUB_RECIPE = """#!/bin/bash
# a stub of seq_script/merge_panels/mega_rebuild.sh <working directory> <merged VCFs>: it reads
# its standard input to its end, as a program of the recipe that reads it would, then builds
# each part the recipe builds, a file of this build in each (or fails, with $STUB_FAIL)
cat > /dev/null
[ -z "${STUB_FAIL:-}" ] || { echo "the recipe failed" >&2; exit 1; }
for p in Genomes/hg38+hg38_mega Genomes/hg38+hg38_mega_INDELS Dictionaries/log_indels_hg38_mega \
    Dictionaries/registry_hg38_mega Dictionaries/indel_af_hg38_mega VCFs/hg38_mega \
    genome_library/20-NRR_2_hg38+hg38_mega genome_library/20-NRR_2_hg38+hg38_mega_INDELS \
    genome_library/20-NRR_2_hg38; do
  mkdir -p "$1/$p" && echo built > "$1/$p/new.txt"
done
mkdir -p "$1/samplesIDs" && echo MEGA > "$1/samplesIDs/hg38_mega.samplesID.txt"
echo "SUMMARY: ALL PASS"
"""


class TheSection12BlockPastedIntoBashOrZsh(unittest.TestCase):
    """docs/TROUBLESHOOTING.md section 12's commands, pasted into bash and into zsh with a stub
    of the recipe that builds each part the recipe builds (test_mega_rebuild pastes them with
    the recipe itself on a stub engine, on Linux: its fixture does not run on macOS, where
    this one runs in macOS's own bash and zsh). The stub reads its standard input to its end:
    bash reads the commands from its own (bash -s), and a program of the recipe that read it
    would take the commands after the recipe, so that none of them ran (the recipe's input is
    /dev/null). In a data folder holding the dataset's indexes of both kinds (matched), none
    (unmatched) or its variant index alone (partly matched), they put the recipe's parts in
    place of the dataset's, each part of the earlier kit gone and nothing else of the folder;
    where the recipe fails, the data folder is as it was; in either shell alike (wave 1c-b of
    2.10.0: they were a subshell, ( set -eu ... ), which zsh ran itself and stopped at its
    first comment, a word in an interactive zsh)."""

    PARTS = [("Genomes", "hg38+hg38_mega"), ("Genomes", "hg38+hg38_mega_INDELS"),
             ("Dictionaries", "dictionaries_hg38_mega"), ("Dictionaries", "log_indels_hg38_mega"),
             ("Dictionaries", "genotypes_hg38_mega"), ("Dictionaries", "indel_genotypes_hg38_mega"),
             ("Dictionaries", "registry_hg38_mega"), ("Dictionaries", "indel_af_hg38_mega"),
             ("VCFs", "hg38_mega")]
    KEPT = [("Genomes", "hg38"), ("genome_library", "20-NRR_2_hg38"),
            ("genome_library", "20-NRR_2_hg38+hg38_HGDP1kGP"),
            ("genome_library", ".20-NRR_2_hg38+hg38_mega"), ("samplesIDs",)]
    STATES = {"matched": ["NGG_2_hg38+hg38_mega", "NGG_2_hg38+hg38_mega_INDELS"],
              "unmatched": [], "partly matched": ["NGG_2_hg38+hg38_mega"]}

    @staticmethod
    def block():
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(os.path.dirname(here), "docs", "TROUBLESHOOTING.md")) as fh:
            text = fh.read()
        start = text.index("**Building the mega's kit again.**", text.index("## 12. "))
        first = text.index("```bash\n", start) + len("```bash\n")
        return text[first:text.index("```", first)]

    def folder(self, indexes):
        """A data folder, the earlier kit's parts and ``indexes`` of the dataset in it, beside
        what is no part of it, a new folder's name, the merged VCFs and a repository whose
        recipe is the stub. Returns (top, env)."""
        top = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, top, True)
        data, repo = os.path.join(top, "data"), os.path.join(top, "repo")
        for part in self.PARTS + [("genome_library", n) for n in indexes] + self.KEPT:
            os.makedirs(os.path.join(data, *part))
            if part[0] != "samplesIDs":
                with open(os.path.join(data, *part, "earlier.txt"), "w") as fh:
                    fh.write("the earlier kit's\n")
        with open(os.path.join(data, "samplesIDs", "hg38_mega.samplesID.txt"), "w") as fh:
            fh.write("its own\n")
        os.makedirs(os.path.join(repo, "PAMs"))
        os.makedirs(os.path.join(repo, "seq_script", "merge_panels"))
        os.makedirs(os.path.join(top, "merged"))
        with open(os.path.join(repo, "PAMs", "20bp-NRR-SpCas9.txt"), "w") as fh:
            fh.write("NNNNNNNNNNNNNNNNNNNNNRR 3\n")
        recipe = os.path.join(repo, "seq_script", "merge_panels", "mega_rebuild.sh")
        with open(recipe, "w") as fh:
            fh.write(STUB_RECIPE)
        os.chmod(recipe, 0o755)
        env = dict(os.environ, DATA=data, NEW=os.path.join(top, "new"),
                   SRC=os.path.join(top, "merged"), REPO=repo)
        return top, env

    def test_matched_unmatched_and_partly_matched_in_bash_and_zsh(self):
        block = self.block()
        for state, indexes in sorted(self.STATES.items()):
            left = {}
            for shell in PASTE_SHELLS:
                with self.subTest(state=state, shell=shell):
                    top, env = self.folder(indexes)
                    data = env["DATA"]
                    rc, out, err = pasted(self, shell, block, top, env)
                    self.assertEqual((rc, shell_errors(err)), (0, []), out + err)
                    self.assertIn("SUMMARY: ALL PASS", out)
                    got = tree(data)
                    # the recipe's parts in place, each part of the earlier kit gone
                    for part in [p for p in self.PARTS if p[1] not in (
                            "dictionaries_hg38_mega", "genotypes_hg38_mega",
                            "indel_genotypes_hg38_mega")] + [
                            ("genome_library", "20-NRR_2_hg38+hg38_mega"),
                            ("genome_library", "20-NRR_2_hg38+hg38_mega_INDELS")]:
                        self.assertEqual(got.get(os.path.join(*part, "new.txt")),
                                         ("file", hashlib.sha1(b"built\n").hexdigest()), part)
                    self.assertFalse([k for k in got if k.endswith("earlier.txt") and (
                        "hg38_mega" in k and not k.startswith("genome_library/.")
                        and "HGDP" not in k)], got)
                    for name in indexes:
                        self.assertNotIn(os.path.join("genome_library", name), got)
                    # what is no part of the dataset stays, the reference index and the
                    # samplesID file already there among it
                    for part in self.KEPT:
                        self.assertIn(os.path.join(*part), got)
                    self.assertEqual(got[os.path.join("genome_library", "20-NRR_2_hg38",
                                                      "earlier.txt")][0], "file")
                    with open(os.path.join(data, "samplesIDs", "hg38_mega.samplesID.txt")) as fh:
                        self.assertEqual(fh.read(), "its own\n")
                    left[shell] = got
            if len(left) == len(PASTE_SHELLS):
                self.assertEqual(left["zsh"], left["bash"], state)

    def test_a_recipe_that_fails_changes_nothing_in_bash_or_zsh(self):
        block = self.block()
        for shell in PASTE_SHELLS:
            with self.subTest(shell=shell):
                top, env = self.folder(self.STATES["matched"])
                before = tree(env["DATA"])
                rc, out, err = pasted(self, shell, block, top, dict(env, STUB_FAIL="1"))
                self.assertIsNotNone(rc, out + err)
                self.assertNotEqual(rc, 0, out + err)
                self.assertIn("the recipe failed", err)
                self.assertEqual(shell_errors(err), [], err)
                self.assertEqual(tree(env["DATA"]), before)
                # and with a variable not set, it stops before anything is made
                top, env = self.folder([])
                before = tree(env["DATA"])
                rc, out, err = pasted(self, shell, block, top, dict(env, SRC=""))
                self.assertNotEqual(rc, 0, out + err)
                self.assertIn("SRC is not set", err)
                self.assertEqual(tree(env["DATA"]), before)
                self.assertFalse(os.path.lexists(env["NEW"]))


class TheDocsRemovalOfHiddenBuildFolders(unittest.TestCase):
    """docs/TROUBLESHOOTING.md section 10 and docs/DOCKER_QUICKSTART.md give the removal of the
    hidden folders a build leaves after a container stop. Up to wave 1c-b of 2.10.0 it was one
    rm -rf of three patterns, which zsh, the default shell of macOS, expands itself: where one
    of them matched nothing (a build leaves its folder in one of the three folders alone, say)
    it stopped the whole command and removed none. A find of each pattern does the same pasted
    into bash or into zsh: with a hidden build folder of each kind (matched), with none
    (unmatched) and with one kind alone (partly matched), each goes, and nothing else (in
    Genomes and Dictionaries a build's own kind alone, as the patterns took them)."""

    KINDS = [("genome_library", ".build-index-only-x1"), ("Genomes", ".build-index-only-enrich-x2"),
             ("Dictionaries", ".build-index-only-dicts-x3")]
    KEPT = [("genome_library", "20-NRR_3_hg38"), ("genome_library", "20-NRR_3_hg38.lock"),
            ("Genomes", "hg38"), ("Genomes", ".build-index-only-x4"),
            ("Dictionaries", "dictionaries_ds"), ("Dictionaries", ".build-index-only-x5")]

    def blocks(self):
        """{doc: its block} of the two docs."""
        here = os.path.dirname(os.path.abspath(__file__))
        out = {}
        for doc, start in (("TROUBLESHOOTING.md", "## 10. "),
                           ("DOCKER_QUICKSTART.md", "Hidden `.build-index-only-*` folders")):
            with open(os.path.join(os.path.dirname(here), "docs", doc)) as fh:
                text = fh.read()
            first = text.index("```bash\n", text.index(start)) + len("```bash\n")
            out[doc] = textwrap.dedent(text[first:text.index("```", first)])
        return out

    def test_each_doc_removes_them_alike_in_bash_and_zsh(self):
        blocks = self.blocks()
        self.assertEqual(blocks["TROUBLESHOOTING.md"], blocks["DOCKER_QUICKSTART.md"])
        self.assertNotIn("rm -rf genome_library/", blocks["TROUBLESHOOTING.md"])
        for doc, block in sorted(blocks.items()):
            for state, made in (("matched", self.KINDS), ("unmatched", []),
                                ("partly matched", self.KINDS[1:2])):
                with self.subTest(doc=doc, state=state):
                    top = tempfile.mkdtemp()
                    self.addCleanup(shutil.rmtree, top, True)
                    for part in made + self.KEPT:
                        os.makedirs(os.path.join(top, *part))
                    paste_both(self, block, top)
                    for part in made:
                        self.assertFalse(os.path.lexists(os.path.join(top, *part)), part)
                    for part in self.KEPT:
                        self.assertTrue(os.path.isdir(os.path.join(top, *part)), part)


class AFakeGenomeAndItsLogs(unittest.TestCase):
    """The fake-indel genome of a chromosome and its indel log must come from one enrichment
    (fake_window_disagreement): each row's window at the offsets its FAKEPOS field gives,
    an N after it, its left flank the reference its refseq field records (a SNP code of the
    overlay standing for its base, Q for any). A download installs an archive's logs beside an
    enrichment made by another CRISPRitz: a 25/27 fake-indel genome read with 30/30 logs, or
    the other way round, disagrees at its first window."""

    REF = (lambda r: "".join(r.choice("ACGT") for _ in range(400)))(
        __import__("random").Random(7))

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="indel_window_fake_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def cut(self, left, right, name, overlay=False):
        """The fake-indel FASTA and indel log an enricher of windows left/right writes of an
        insertion at 100, a deletion at 200 and an MNV at 300 (overlay: every flank base
        coded, as the SNP overlay codes it, the first as Q)."""
        recs = [(100, self.REF[99], self.REF[99] + "TT"), (200, self.REF[199:202], self.REF[199]),
                (300, self.REF[299:301], "".join({"A": "C"}.get(b, "A") for b in self.REF[299:301]))]
        code = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
        flat, rows = "", []
        for k, (pos, ref, alt) in enumerate(recs):
            start = pos - 1 - left
            lf = self.REF[start:pos - 1]
            if overlay:
                lf = "Q" + "".join(code[b] for b in lf[1:])
            win = lf + alt + self.REF[pos - 1 + len(ref):pos - 1 + len(ref) + right]
            rows.append("chrT_%d-%d_%d\tS1\trs%d\t0.1\tchrT_%d_%s_%s\t%d,%d\t%s" % (
                start, pos - 1 + len(ref) + right, k + 1, k, pos, ref, alt, len(flat),
                len(flat) + len(win), self.REF[start:start + len(win)]))
            flat += win + "N"
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        fa = os.path.join(d, "fakechrT.fa")
        with open(fa, "w") as fh:
            fh.write(">fakechrT\n" + "".join(w + "\nN\n" for w in flat.split("N")[:-1]))
        log = os.path.join(d, "logchrT.txt.gz")
        with gzip.open(log, "wt") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n" + "\n".join(rows) + "\n")
        return fa, log

    def test_one_enrichment_agrees(self):
        for left, right in ((25, 27), (30, 30)):
            for overlay in (False, True):
                fa, log = self.cut(left, right, "w%d_%d_%s" % (left, right, overlay), overlay)
                self.assertIsNone(iw.fake_window_disagreement(fa, log), (left, overlay))
                self.assertIsNone(iw.fake_window_disagreement(fa, log, rows=None))

    def test_two_enrichments_disagree(self):
        fa27, log27 = self.cut(25, 27, "w25", overlay=True)
        fa30, log30 = self.cut(30, 30, "w30", overlay=True)
        how = iw.fake_window_disagreement(fa27, log30)
        self.assertEqual(how, "the window of chrT_100_%s_%sTT, which the log places at 0-63 of "
                              "the fake sequence, ends elsewhere there (no separator after it)"
                         % (self.REF[99], self.REF[99]))
        self.assertIn("ends elsewhere there", iw.fake_window_disagreement(fa30, log27))
        # the right length, but another flank: a window of 30 cut 3 bases later
        with open(fa30) as fh:
            text = fh.read()
        lines = text.split("\n")
        lines[1] = self.REF[99 - 27:99 - 27 + len(lines[1])]
        with open(fa30, "w") as fh:
            fh.write("\n".join(lines))
        self.assertIn("does not hold there the 30 reference bases before the record",
                      iw.fake_window_disagreement(fa30, log30))
        # nothing to compare: no log, a log of no row, no FASTA
        self.assertIsNone(iw.fake_window_disagreement(fa27, None))
        empty = os.path.join(self.tmp, "logempty.txt")
        with open(empty, "w") as fh:
            fh.write("CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq\n")
        self.assertIsNone(iw.fake_window_disagreement(fa27, empty))
        self.assertIsNone(iw.fake_window_disagreement(os.path.join(self.tmp, "none.fa"), log27))

    def test_a_dataset_of_two_enrichments_is_refused(self):
        """fake_genome_refusal, as build-index-only and the search script take it."""
        w = os.path.join(self.tmp, "work")
        fa, log = self.cut(25, 27, "w25", overlay=True)
        _, log30 = self.cut(30, 30, "w30")
        fake = os.path.join(w, "Genomes", "hg38+ds_INDELS", "fake_ds_chrT")
        os.makedirs(fake)
        shutil.copy(fa, fake)
        logs = os.path.join(w, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        shutil.copy(log, logs)
        self.assertIsNone(iw.fake_genome_refusal(w, "hg38", "ds"))
        shutil.copy(log30, logs)
        why = iw.fake_genome_refusal(w, "hg38", "ds", who="this search")
        self.assertTrue(why.startswith(
            "the fake-indel genome %s_INDELS and the indel logs %s of ds come from two "
            "enrichments: on chrT, the window of chrT_100_" % (
                os.path.join(w, "Genomes", "hg38+ds"), logs)), why)
        self.assertIn("so an index that this search built from the enrichment in place", why)
        self.assertIn("(rm -rf '%s' '%s_INDELS'; of a published index" % (
            os.path.join(w, "Genomes", "hg38+ds"), os.path.join(w, "Genomes", "hg38+ds")), why)
        self.assertIn("and run build-index-only again, which enriches it anew.", why)
        # W-R3-F1: a sites-only dataset (its registry of INFO AF fields, no data_type) is not
        # sent to build-index-only, which would enrich it with no window
        reg = os.path.join(w, "Dictionaries", "registry_ds")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chrT.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        why = iw.fake_genome_refusal(w, "hg38", "ds")
        self.assertNotIn("build-index-only again", why)
        # (the removal of what comes from an enrichment names no registry; the recipe's own
        # replacement of it, IR1-F2, is after)
        self.assertNotIn(reg, why.split("and build that kit again")[0])
        self.assertTrue(why.endswith(
            "or remove all of the dataset that comes from an enrichment (rm -rf '%s' '%s') and "
            "build that kit again, its indels carried by a synthetic sample, %s: ds is "
            "sites-only: its VCFs name no sample, and CRISPRitz writes a fake-indel window only "
            "for an indel a sample carries, so an enrichment of them holds none, so "
            "build-index-only would enrich it anew with no indel off-target." % (
                os.path.join(w, "Genomes", "hg38+ds_INDELS"), logs,
                iw.rebuild_kit_way("ds", w, "hg38"))), why)

    def test_the_recipe_the_advice_names_builds_the_whole_kit(self):
        """The advice on a sites-only dataset says to build its kit again with a synthetic
        sample, in place of all of it that comes from an enrichment (the enriched genome and
        its indexes among them). The recipe it names must build the whole kit of 2.10.0, the
        one the index set builds the mega with (seq_script/rebuild/index_set_2.10.0.tsv):
        lane M's mega_rebuild.sh, one build-index-only of the enriched genome, the SNP index,
        the fake-indel genome and its INDELS index. It named mega_build_indels.sh, up to
        2.8.1's recipe of the fake-indel half alone, after which the next search found no
        enriched genome and was refused, the dataset sites-only and its new logs holding
        windows (the integration of stage 2b)."""
        import csv
        import re
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        named = set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh",
                               iw.sites_only_way("hg38_mega")))
        # the fake-genome refusal of a sites-only dataset (test_a_dataset_of_two_enrichments_...)
        w = os.path.join(self.tmp, "work")
        fa, log = self.cut(25, 27, "w25", overlay=True)
        _, log30 = self.cut(30, 30, "w30")
        fake = os.path.join(w, "Genomes", "hg38+ds_INDELS", "fake_ds_chrT")
        os.makedirs(fake)
        shutil.copy(fa, fake)
        logs = os.path.join(w, "Dictionaries", "log_indels_ds")
        os.makedirs(logs)
        shutil.copy(log30, logs)
        reg = os.path.join(w, "Dictionaries", "registry_ds")
        os.makedirs(reg)
        with open(os.path.join(reg, "reg_chrT.idx"), "w") as fh:
            fh.write('{"aggregation": "info_af"}')
        why = iw.fake_genome_refusal(w, "hg38", "ds")
        self.assertIn("and build that kit again, its indels carried by a synthetic sample, "
                      "with seq_script/merge_panels/mega_rebuild.sh", why)
        named |= set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh", why))
        manifest = os.path.join(repo, "seq_script", "rebuild", "index_set_2.10.0.tsv")
        with open(manifest) as fh:
            lines = [x for x in fh if x.strip() and not x.startswith("#")]
        mega = [r for r in csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)
                if r["panel"] == "hg38_mega" and r["kind"] == "variant"]
        self.assertTrue(mega)
        builds = set(re.findall(r"seq_script/merge_panels/[A-Za-z0-9_]+\.sh",
                                " ".join(r["build"] for r in mega)))
        self.assertEqual(named, builds)
        for script in named:
            self.assertTrue(os.path.isfile(os.path.join(repo, script)), script)

    def test_the_overlay_refuses_a_genome_its_log_does_not_describe(self):
        """W-R3-F2: the SNP overlay places each code by the log's windows. Given the logs of
        another enrichment (30/30 beside a 25/27 fake-indel genome) it wrote 641,930 codes
        at the wrong bases of a real chr22 genome; it raises NotItsLog now, the FASTA left as
        it was, and the search's overlay step says so on stderr, exit 1. With the genome's
        own log it overlays as before."""
        import io
        from contextlib import redirect_stderr, redirect_stdout

        import overlay_indel_snps
        import search_indel_snp

        fa27, log27 = self.cut(25, 27, "w25")
        _, log30 = self.cut(30, 30, "w30")
        enriched = os.path.join(self.tmp, "chrT.enriched.fa")
        code = {"A": "R", "C": "Y", "G": "R", "T": "Y"}
        with open(enriched, "w") as fh:
            fh.write(">chrT\n" + "".join(code[b] for b in self.REF) + "\n")
        with open(fa27) as fh:
            before = fh.read()
        with self.assertRaises(overlay_indel_snps.NotItsLog) as said:
            overlay_indel_snps.overlay_fake_chromosome(fa27, enriched, log30)
        self.assertIn("come from two enrichments: the window of chrT_100_", str(said.exception))
        with open(fa27) as fh:
            self.assertEqual(fh.read(), before)
        # the search's step (search_indel_snp.py overlay <genome> <enriched> <logs>)
        top = os.path.join(self.tmp, "search")
        for d, f, src in (("genome/fake_ds_chrT", "fakechrT.fa", fa27),
                          ("enriched", "chrT.enriched.fa", enriched),
                          ("logs", "logchrT.txt.gz", log30)):
            os.makedirs(os.path.join(top, d))
            shutil.copy(src, os.path.join(top, d, f))
        err = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(err):
            rc = search_indel_snp.main(["search_indel_snp.py", "overlay",
                                        os.path.join(top, "genome"), os.path.join(top, "enriched"),
                                        os.path.join(top, "logs")])
        self.assertEqual(rc, 1)
        self.assertTrue(err.getvalue().startswith("ERROR: the fake-indel genome "), err.getvalue())
        self.assertIn("come from two enrichments", err.getvalue())
        with open(os.path.join(top, "genome", "fake_ds_chrT", "fakechrT.fa")) as fh:
            self.assertEqual(fh.read(), before)
        # its own log: overlaid
        self.assertGreater(overlay_indel_snps.overlay_fake_chromosome(fa27, enriched, log27), 0)
        self.assertIsNone(iw.fake_window_disagreement(fa27, log27))


if __name__ == "__main__":
    unittest.main(verbosity=2)

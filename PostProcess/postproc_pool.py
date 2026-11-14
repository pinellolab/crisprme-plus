#!/usr/bin/env python
"""Option P (2.10.0): a post-analysis stage run in pieces within each chromosome, over one
pool of workers for every chromosome of the stage.

Up to 2.8.1 the SNP and the indel post-analysis ran one process a chromosome
(pool_post_analisi_snp.py, pool_post_analisi_indel.py), so a dense locus (IGL on chr22,
IGH on chr14) ran on one core while the others waited: on chr22 at 4 mismatches + 1 DNA
+ 1 RNA bulge the SNP stage took 100 minutes, nearly all of it IGL. Here each chromosome
goes through three steps on the pool's workers:

1. its shell driver's "cluster" stage: the extraction and the clustering of its raw
   targets, as before; the driver writes the analysis script's arguments to a manifest, a
   file of a folder the pool makes in $TMPDIR (tempfile.mkdtemp: a new name, 0700), which no
   one else can replace where $TMPDIR is sticky, as /tmp is, or the user's own;
2. the analysis in pieces: the clustered file's lines cut into ``n`` contiguous ranges
   (``n_pieces``), the script run once a range (CRISPRME_POSTPROC_PIECE=1 and
   CRISPRME_POSTPROC_RANGE), each piece starting up exactly as one process does over the
   whole file, so each line is analysed as one process analyses it. The SNP side first
   reads the dictionary keys of the whole chromosome once (CRISPRME_POSTPROC_PREPARE), for
   the pieces to read them from a small copy (CRISPRME_POSTPROC_DICT) instead of streaming
   the whole dictionary each: a JSON file (write_dict_copy, read_dict_copy) in the pieces'
   folder, which loads as data whoever wrote it. Then merge_pieces.py joins the pieces' files
   into the chromosome's, byte for byte those of one process (the CRISPR-Bulge scores in the
   batches one process forms);
3. the driver's "finish" stage (adjust_cols.py and, for indels, remove_bad_indel_targets.py).

A chromosome with too few lines for two pieces runs its analysis as one process, as
before. CRISPRME_POSTPROC_PIECES=0 runs every chromosome as one process (the drivers' "all"
stage), as up to 2.8.1. CRISPRME_POSTPROC_PIECE_LINES sets the fewest lines a piece gets
(default: 16,000 for the SNP stage, 50,000 for the indel stage, whose lines cost less).

Stdout carries the notes (log_verbose.txt). A piece's own output goes to a file of the
pieces' folder (each piece starts up as one process does, so the pieces would repeat it):
once the pieces of a chromosome are done, the log takes the first piece's output and any
line another piece wrote that the first did not, but for the pieces' run times (piece_notes),
then merge_pieces.py writes the notes one process ends with, from the pieces' counts, and the
pool the chromosome's own run time ("Post-analysis SNP of chr22: done in 101.2 s (16
pieces)"). A failed step writes its error
on stderr, with the last lines of its output, and the stage stops (the search treats any
stderr write as a failure).
"""

import heapq
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import haplotype_precheck as _hpc


def pieces_on():
    """Option P is on unless CRISPRME_POSTPROC_PIECES says off: 0, false, no or off, in any
    case (haplotype_precheck.switch_on, the rule of CRISPRME_HAPLOTYPE_PRECHECK too)."""
    return _hpc.switch_on(os.environ, "CRISPRME_POSTPROC_PIECES")


def piece_lines(default):
    """The fewest lines a piece gets (CRISPRME_POSTPROC_PIECE_LINES, else ``default``)."""
    try:
        v = int(os.environ.get("CRISPRME_POSTPROC_PIECE_LINES", default))
        return v if v > 0 else default
    except ValueError:
        return default


def n_pieces(lines, min_lines, workers):
    """How many pieces a chromosome of ``lines`` data lines is cut into: as many as the
    pool has workers, each of at least ``min_lines`` lines; one for a small chromosome."""
    if lines <= 0:
        return 1
    return max(1, min(workers, lines // min_lines))


def ranges(lines, n):
    """``n`` contiguous ranges (first, last), 1-based and inclusive, of ``lines`` lines,
    as even as they can be."""
    out = []
    first = 1
    for i in range(n):
        size = lines // n + (1 if i < lines % n else 0)
        out.append((first, first + size - 1))
        first += size
    return out


def data_lines(path):
    """The lines of a clustered targets file after its header line."""
    n = 0
    with open(path, "rb") as fh:
        fh.readline()
        for _ in fh:
            n += 1
    return n


# The SNP side's dictionary copy: the dictionary keys a chromosome's lines need and its
# phasing, which the start-up (CRISPRME_POSTPROC_PREPARE) reads once for the chromosome's
# pieces (CRISPRME_POSTPROC_DICT). It is a file of the run's output folder, which other
# members of that folder's group can write where the umask lets them (002 on many shared
# nodes), so loading it must run no code whoever wrote it: JSON, the dictionary's own
# format, which loads as data only. (The first copy, a pickle, ran whatever its file named
# when a piece loaded it: 2.10.0 lane S review r3.) It is created 0600, never through a name
# that already stood (O_EXCL), in the pieces' folder (made 0700). That keeps others from
# reading or writing this copy, not from changing what a piece reads: in an output folder
# they can write, they can rename the pieces' folder and put one of their own in its place
# (a rename needs only write access to the folder that holds it), with a copy of their own
# and links the pieces then write through (docs/TROUBLESHOOTING.md section 7).
DICT_COPY_NAME = "dictionary.json"


def write_dict_copy(path, haplo, entries):
    """Write the dictionary copy ``path``: {"haplo": the phasing, "dict": the entries}, as
    JSON, created 0600 beside it and then renamed onto it. Entries JSON cannot hold exactly
    (a number other than an int or a float, a float that is not finite) raise TypeError or
    ValueError before any file is made; the pieces then read the dictionary themselves, as
    one process does."""
    text = json.dumps({"haplo": bool(haplo), "dict": entries}, allow_nan=False)
    tmp = path + ".tmp"
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def read_dict_copy(path):
    """The dictionary copy ``path`` (write_dict_copy) as (haplo, entries); a file that is not
    one raises ValueError, so a piece never runs without it."""
    try:
        with open(path, encoding="utf-8") as fh:
            copy = json.load(fh)
    except ValueError as err:  # not UTF-8 JSON (json.JSONDecodeError, UnicodeDecodeError)
        raise ValueError(f"{path} is not a dictionary copy of the post-analysis: {err}") \
            from None
    if (not isinstance(copy, dict) or not isinstance(copy.get("haplo"), bool)
            or not isinstance(copy.get("dict"), dict)):
        raise ValueError(f"{path} is not a dictionary copy of the post-analysis")
    return copy["haplo"], copy["dict"]


class TaskPool(object):
    """``workers`` threads running tasks from a priority queue (lower first); a task may
    queue more. ``run`` returns when no task is queued or running, or raises the first
    error (no new task starts once one failed)."""

    def __init__(self, workers):
        self.workers = max(1, workers)
        self._heap = []
        self._seq = 0
        self._running = 0
        self._error = None
        self._cv = threading.Condition()

    def put(self, priority, fn, *args):
        with self._cv:
            heapq.heappush(self._heap, (priority, self._seq, fn, args))
            self._seq += 1
            self._cv.notify()

    def _worker(self):
        while True:
            with self._cv:
                while not self._heap and self._running and self._error is None:
                    self._cv.wait()
                if self._error is not None or not self._heap:
                    self._cv.notify_all()
                    return
                _p, _s, fn, args = heapq.heappop(self._heap)
                self._running += 1
            try:
                fn(*args)
            except BaseException as err:  # noqa: BLE001 - reported by run()
                with self._cv:
                    if self._error is None:
                        self._error = err
            finally:
                with self._cv:
                    self._running -= 1
                    self._cv.notify_all()

    def run(self):
        threads = [threading.Thread(target=self._worker, daemon=True)
                   for _ in range(self.workers)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        if self._error is not None:
            raise self._error


class StageError(RuntimeError):
    pass


_OUT = threading.Lock()


def _say(lines):
    """Write ``lines`` on stdout at once (the pool's threads share it)."""
    if lines:
        with _OUT:
            sys.stdout.write("".join(line + "\n" for line in lines))
            sys.stdout.flush()


def _call(cmd, env=None, what=""):
    code = subprocess.call(cmd, env=env)
    if code != 0:
        raise StageError(f"{what} failed (exit {code}): {' '.join(cmd)}")


def _tail(path, n=20):
    try:
        with open(path, errors="replace") as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return ""


def _call_logged(cmd, log, env=None, what=""):
    """Run ``cmd`` with its stdout in the file ``log``; a failure raises StageError with
    the last lines of it."""
    with open(log, "w") as out:
        code = subprocess.call(cmd, env=env, stdout=out)
    if code != 0:
        raise StageError(f"{what} failed (exit {code}): {' '.join(cmd)}; the last lines of "
                         f"its output ({log}):\n{_tail(log)}")


_NUMBERS = re.compile(r"[0-9]+")
# A piece's run time, which its script prints last: "ANALYSIS COMPLETE IN <s>" (both
# scripts), after "Done <chrom> <s>" (analisi_indels_NNN.py). The chromosome's own is the
# pool's "... done in <s> s" line (run_stage).
_RUN_TIME = re.compile(r"ANALYSIS COMPLETE IN |Done \S+ [0-9.eE+-]+$")


def piece_notes(logs, pdir):
    """The lines of a chromosome's pieces' outputs (``logs``, in piece order) its log keeps:
    the first piece's, which say what one process says of its start-up (each piece starts up
    as one process does), and any line another piece wrote that the first did not, numbers
    aside (one of each); never a line that names the pieces' folder ``pdir`` (a piece's own
    files) or a piece's run time (_RUN_TIME: its last lines)."""
    out, seen = [], set()
    for i, log in enumerate(logs):
        try:
            with open(log, errors="replace") as fh:
                lines = fh.read().splitlines()
        except OSError:
            continue
        for line in lines:
            if pdir in line or _RUN_TIME.match(line):
                continue
            key = _NUMBERS.sub("#", line)
            if i and key in seen:
                continue
            seen.add(key)
            out.append(line)
    return out


def run_stage(kind, chroms, workers, driver_args, script, prepare_dict=False,
              min_lines=16000):
    """Run one post-analysis stage over ``chroms`` with ``workers`` workers.

    ``kind``: "snp" or "indel" (merge_pieces.py's). ``driver_args(chrom, stage)``: the
    command of the chromosome's shell driver for a stage ("all", "cluster", "finish").
    ``script``: the analysis script ("./new_simple_analysis.py", "./analisi_indels_NNN.py");
    its manifest arguments are its command line, the clustered file at index 1 (SNP: 2) and
    the output prefix at index 4 (SNP) or 2 (indel). ``prepare_dict``: the SNP side's
    dictionary copy, read once a chromosome."""
    pool = TaskPool(workers)
    lock = threading.Lock()
    if not pieces_on():
        for chrom in chroms:
            pool.put(1, lambda c=chrom: _call(driver_args(c, "all"),
                                              what=f"post-analysis of {c}"))
        pool.run()
        return
    cluster_at, out_at = (2, 4) if kind == "snp" else (1, 2)
    base_env = dict(os.environ)
    started = {}  # chromosome -> when its cluster stage started (its "done in" line)

    def finish(chrom, pdir, cluster, n=1):
        _call(driver_args(chrom, "finish"), what=f"post-analysis of {chrom} (finish)")
        if pdir is not None:
            shutil.rmtree(pdir, ignore_errors=True)
        if os.path.exists(cluster):
            os.remove(cluster)
        if n > 1:
            # the chromosome's run time, which no piece's own line tells (piece_notes)
            _say([f"Post-analysis {kind.upper()} of {chrom}: done in "
                  f"{time.time() - started[chrom]:.1f} s ({n} pieces)"])

    def merge(chrom, args, pdir, prefixes):
        _say(piece_notes([p + ".log" for p in prefixes], pdir))
        # it writes the notes one process ends with, from the pieces' counts
        _call([sys.executable, "merge_pieces.py", kind, args[out_at]] + prefixes,
              what=f"merge of the pieces of {chrom}")
        finish(chrom, pdir, args[cluster_at], len(prefixes))

    def piece(chrom, args, pdir, prefixes, i, rng, dict_copy, left):
        a = list(args)
        a[out_at] = prefixes[i]
        env = dict(base_env)
        env["CRISPRME_POSTPROC_PIECE"] = "1"
        env["CRISPRME_POSTPROC_RANGE"] = "%d:%d" % rng
        if dict_copy is not None and os.path.exists(dict_copy):
            env["CRISPRME_POSTPROC_DICT"] = dict_copy
        _call_logged([script] + a, prefixes[i] + ".log", env=env,
                     what=f"piece {i + 1} of the post-analysis of {chrom}")
        with lock:
            left[0] -= 1
            done = left[0] == 0
        if done:
            pool.put(0, merge, chrom, args, pdir, prefixes)

    def pieces(chrom, args, pdir, rngs, dict_copy):
        prefixes = [os.path.join(pdir, "p%04d" % i) for i in range(len(rngs))]
        left = [len(rngs)]
        for i, rng in enumerate(rngs):
            pool.put(1, piece, chrom, args, pdir, prefixes, i, rng, dict_copy, left)

    def prepare(chrom, args, pdir, rngs):
        a = list(args)
        a[out_at] = os.path.join(pdir, "prepare")
        env = dict(base_env)
        dict_copy = os.path.join(pdir, DICT_COPY_NAME)
        env["CRISPRME_POSTPROC_PREPARE"] = dict_copy
        # its output is the start-up's, which the first piece's repeats
        _call_logged([script] + a, os.path.join(pdir, "prepare.log"), env=env,
                     what=f"start-up of the post-analysis of {chrom}")
        pieces(chrom, args, pdir, rngs, dict_copy)

    def cluster(chrom, i):
        started[chrom] = time.time()
        env = dict(base_env)
        # the driver writes the analysis' arguments there, one a line: a file of the folder
        # the pool made (tempfile.mkdtemp: a new name, mode 0700), so no other user of a
        # shared $TMPDIR can name it before the driver does, nor, where $TMPDIR is sticky (as
        # /tmp) or the user's own, put another list there
        tmp = os.path.join(manifests, "%s_%d.args" % (kind, i))
        env["CRISPRME_POSTPROC_MANIFEST"] = tmp
        code = subprocess.call(driver_args(chrom, "cluster"), env=env)
        if code != 0:
            raise StageError(f"post-analysis of {chrom} (cluster) failed (exit {code})")
        if not os.path.exists(tmp):
            raise StageError(f"post-analysis of {chrom} (cluster) wrote no arguments for "
                             f"its analysis ({tmp})")
        with open(tmp) as fh:
            manifest = [ln.rstrip("\n") for ln in fh]
        os.remove(tmp)
        args = manifest
        cluster_file = args[cluster_at]
        lines = data_lines(cluster_file)
        n = n_pieces(lines, min_lines, workers)
        if n <= 1:
            _call([script] + args, what=f"post-analysis of {chrom}")
            pool.put(0, finish, chrom, None, cluster_file)
            return
        pdir = args[out_at] + ".pieces"
        shutil.rmtree(pdir, ignore_errors=True)
        # the pieces' files, the dictionary copy among them, in a folder made here, 0700:
        # others cannot read or write the files in it, but in an output folder they can
        # write (a umask of 002 lets the group) they can rename it and put their own in its
        # place, whose files the pieces then read and write through (DICT_COPY_NAME's note).
        # Anything left under this name that the line above could not remove (another
        # user's folder, a link) stops the stage (FileExistsError), and is never used.
        os.makedirs(pdir, mode=0o700)
        rngs = ranges(lines, n)
        _say([f"Post-analysis {kind.upper()} of {chrom}: {lines} line(s) in {n} piece(s)"])
        dict_path = args[1] if kind == "snp" else None
        if prepare_dict and dict_path and (os.path.exists(dict_path)
                                           or os.path.exists(dict_path + ".gz")):
            pool.put(1, prepare, chrom, args, pdir, rngs)
        else:
            pieces(chrom, args, pdir, rngs, None)

    manifests = tempfile.mkdtemp(prefix="crisprme_postproc_%s_" % kind)
    try:
        for i, chrom in enumerate(chroms):
            pool.put(2, cluster, chrom, i)
        pool.run()
    finally:
        shutil.rmtree(manifests, ignore_errors=True)

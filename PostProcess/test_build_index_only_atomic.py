#!/usr/bin/env python3
"""build-index-only never leaves a half-built index or enriched genome behind.

Every later run -- build-index-only and complete-search alike -- takes an index
directory, or an enriched genome, that is there for a finished one. build-index-only
used to run the engine straight into genome_library/<name>, and to move the enriched
genome into place before the rest of the enrichment, so a run that ended early
(Ctrl-C, an OOM kill, a container stop, an engine error) left a partial one there,
which every later search then used, with no error and fewer off-targets. It also
went ahead while a search held the build lock of what it was about to use or build.

These tests run the real crisprme.py build-index-only, and the real
pool_index_indels.py under it, with a stub engine on PATH. The stub crispritz.py
writes the layout the real one writes, each .bin a tree as buildTST saves one, and
can be told to fail, or to wait, before one chromosome, or to cut its .bin short;
the stub pigz compresses as pigz does and can be told to wait. A
wait ends on a kill, a signal or a release file, so every case is deterministic.
SameTreeAsBefore pins the successful builds to what the code before the fix
produced with the same stub, file for file.

Stdlib only: crisprme.py imports Bio.Seq at module scope, so a stub Bio package goes
on PYTHONPATH (build-index-only never uses it).
"""

import errno
import gzip
import hashlib
import json
import os
import random
import resource
import shutil
import signal
import stat
import struct
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CRISPRME = os.path.join(REPO, "crisprme.py")
PAM = os.path.join(REPO, "PAMs", "20bp-NRG-SpCas9.txt")
# --bDNA 1 --bRNA 1 -> indexes built with -bMax 2
REF = "NRG_2_hg38"
VAR = "NRG_2_hg38+ds"
INDELS = VAR + "_INDELS"
LOG_HEADER = "CHR\tSAMPLES\trsID\tAF\tindel\tFAKEPOS\trefseq"
GENOME = {"chr1": "ACGTTGCAAG" * 12, "chr2": "AGCTTCAGGT" * 12}

STUB_ENGINE = r'''#!PYTHON
"""Stub crispritz.py: the layout index-genome and add-variants write, and no more.

Each output file is a step "<command> <name> <input file>". A step whose key holds
$STUB_FAIL fails (exit 1) before it is written; one whose key holds $STUB_HANG
waits there, once $STUB_MARK is made, until $STUB_RELEASE exists. index-genome
writes each .bin as buildTST does, named by its contig's header, a tree of one
target (the whole sequence), and none for a contig of N alone; one whose step's key
holds $STUB_CUT it cuts to $STUB_CUT_TO bytes (half by default), and one whose key
holds $STUB_NO_BIN it does not write, and still exits 0, as buildTST does on a full
disk, or with the disk's inode table full. add-variants pads the fake
contig of chromosome $STUB_BIG_FAKE with a megabase of N (and this stub lifts any
file-size limit of its own), and with $STUB_GZ_BLOCKED leaves a directory where the
last dictionary's .gz would go.
"""
import json, os, resource, shutil, struct, sys, time

if os.environ.get("STUB_BIG_FAKE"):
    resource.setrlimit(resource.RLIMIT_FSIZE, (resource.getrlimit(resource.RLIMIT_FSIZE)[1],) * 2)

# the four bits a tree gives a base (writePair, in CRISPRitz's mainParallel.cpp), N's 0000
CODE = dict(zip("ACGTRYSWKMBDHV", (1, 2, 4, 8, 5, 10, 6, 9, 12, 3, 14, 13, 11, 7)))


def tst(pam, nb, seq):
    """A .bin as saveTST writes one, of a single target: the sequence seq, with pam
    for its PAM and nb in the guide-length field; its tree is a chain of one node a
    base, which ends with the target's index (crispritz_bins has the layout)."""
    nibbles = [n for base in seq for n in (CODE.get(base, 0), 0, 0)]  # a base, no lo/hi kid
    tree = bytes(a << 4 | b for a, b in zip(nibbles[::2], nibbles[1::2]))
    tree += bytes([nibbles[-1] << 4 | 0xF if len(nibbles) % 2 else 0xF0])  # "_": the end
    codes = [CODE.get(base, 0) for base in pam] + [0]
    target = struct.pack("=i", 0) + bytes(codes[i] << 4 | codes[i + 1]
                                          for i in range(0, len(pam), 2)) + b"0"
    return (struct.pack("=ii", 1, int(nb)) + target + struct.pack("=i", len(seq)) + tree
            + struct.pack("=i", -1))


def note(line):
    with open(os.environ["STUB_LOG"], "a") as fh:
        fh.write(line + "\n")


def step(key):
    if os.environ.get("STUB_FAIL") and os.environ["STUB_FAIL"] in key:
        sys.stderr.write("stub: %s failed\n" % key)
        sys.exit(1)
    if os.environ.get("STUB_HANG") and os.environ["STUB_HANG"] in key:
        open(os.environ["STUB_MARK"], "w").close()
        while not os.path.exists(os.environ["STUB_RELEASE"]):
            time.sleep(0.02)


def fastas(d):
    return sorted(f for f in os.listdir(d) if f.split(".")[-1] == "fa")


def sequence(path):
    with open(path) as fh:
        fh.readline()
        return "".join(line.strip() for line in fh)


def header(path):
    with open(path) as fh:
        return fh.readline().rstrip("\n")[1:]


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


note(" ".join(sys.argv[1:]))
if sys.argv[1] == "index-genome":
    name, gdir = sys.argv[2], os.path.realpath(sys.argv[3])
    nb = sys.argv[sys.argv.index("-bMax") + 1]
    with open(sys.argv[4]) as fh:
        motif, size = fh.read().split()[:2]
    size = int(size)
    pam = motif[:-size] if size < 0 else motif[-size:]
    out = os.path.join("genome_library", "%s_%s_%s" % (pam, nb, name))
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    for fa in fastas(gdir):
        key = "index-genome %s %s" % (name, fa)
        step(key)
        seq = sequence(os.path.join(gdir, fa))
        if not seq.strip("N"):  # (buildTST writes no tree where there is no target)
            continue
        if os.environ.get("STUB_NO_BIN") and os.environ["STUB_NO_BIN"] in key:
            continue
        data = tst(pam, nb, seq)
        if os.environ.get("STUB_CUT") and os.environ["STUB_CUT"] in key:
            data = data[: int(os.environ.get("STUB_CUT_TO", len(data) // 2))]
        with open(os.path.join(out, "%s_%s_1.bin" % (pam, header(os.path.join(gdir, fa)))),
                  "wb") as fh:
            fh.write(data)
elif sys.argv[1] == "add-variants":
    vdir, gdir = os.path.realpath(sys.argv[2]), os.path.realpath(sys.argv[3])
    vcf, genome = os.path.basename(vdir), os.path.basename(gdir)
    snps = os.path.join("variants_genome", "SNPs_genome")
    enriched = os.path.join(snps, genome + "_enriched")
    os.makedirs(enriched)
    os.makedirs(os.path.join("variants_genome", "INDELs_genome"))
    write(os.path.join("variants_genome", "INDELs_genome", "change_version.txt"), "obsolete\n")
    for fa in fastas(gdir):
        chrom, seq = fa[: -len(".fa")], sequence(os.path.join(gdir, fa))
        step("add-variants %s %s" % (vcf, fa))
        # a SNP at 11 (A>G), written as its IUPAC code, and an insertion at 26,
        # cut out with its flanks as a fake contig
        write(os.path.join(enriched, chrom + ".enriched.fa"),
              ">%s\n%sR%s\n" % (chrom, seq[:10], seq[11:]))
        write(os.path.join(snps, "my_dict_%s.json" % chrom),
              json.dumps({"%s,11" % chrom: "S1:0|1;A,G;rs_%s_11;0.25" % chrom}))
        fake = seq[:26] + "T" + seq[26:60]
        write(os.path.join(snps, "log%s.txt" % chrom),
              "%s\n%s_0-60_1\tS2\trs_%s_26\t0.25\t%s_26_%s_%sT\t0,%d\t%s\n"
              % (LOG_HEADER, chrom, chrom, chrom, seq[25], seq[25], len(fake), seq[:60]))
        if os.environ.get("STUB_BIG_FAKE") == chrom:
            fake += "N" * 1000000
        os.makedirs(os.path.join("variants_genome", "fake_%s_%s" % (vcf, chrom)))
        write(os.path.join("variants_genome", "fake_%s_%s" % (vcf, chrom), "fake%s.fa" % chrom),
              ">fake%s\n%s\n" % (chrom, fake))
    if os.environ.get("STUB_GZ_BLOCKED"):
        os.makedirs(os.path.join(snps, "my_dict_%s.json.gz" % chrom))
else:
    sys.exit("stub crispritz.py: no %s here" % sys.argv[1])
'''.replace("LOG_HEADER", repr(LOG_HEADER))

STUB_PIGZ = r'''#!PYTHON
"""Stub pigz: compresses each file to <file>.gz and removes it. With
$STUB_PIGZ_HANG set, it first waits, once $STUB_MARK is made, until $STUB_RELEASE
exists; with $STUB_PIGZ_FAIL, it fails after the first file, as on a write error."""
import gzip, os, sys, time

if os.environ.get("STUB_PIGZ_HANG"):
    open(os.environ["STUB_MARK"], "w").close()
    while not os.path.exists(os.environ["STUB_RELEASE"]):
        time.sleep(0.02)
args, files = iter(sys.argv[1:]), []
for a in args:
    if a == "-p":
        next(args)
    elif not a.startswith("-"):
        files.append(a)
for i, f in enumerate(files):
    if i and os.environ.get("STUB_PIGZ_FAIL"):
        sys.exit("pigz: abort: write error on %s.gz" % f)
    with open(f, "rb") as src, open(f + ".gz", "wb") as dst:
        dst.write(gzip.compress(src.read(), mtime=0))
    os.remove(f)
'''

# Imported first by each Python process of a run started with shim=True
# (sitecustomize), it stands in for what a test cannot set up for real.
# STUB_MOUNT=<dir>: <dir> is a mount point of its own, so a rename into or out of it
# fails with EXDEV, as between two bind mounts of one file system, whose st_dev is
# the same. STUB_SIGNAL_AT_CLEANUP=before: a SIGTERM comes just as the build's
# cleanup starts, before it holds signals back (as it first replaces the build's
# SIGTERM handler); =after: one comes just after that. STUB_OVERLAY_HANG=1: the SNP
# overlay, about to write a FASTA's temporary file, waits there once $STUB_MARK is
# made, until $STUB_RELEASE exists. STUB_SLOW_CLEANUP=<file>: the cleanup, about to
# remove what a build directory holds, makes <file> and waits until $STUB_RELEASE
# exists, as removing a large partial index takes a while. STUB_THREAD=1: a thread
# that does nothing, as numpy's OpenBLAS starts at import, for the kernel to give a
# signal to when the one that runs Python has it blocked. STUB_WRITE_HANG=<text>: a
# file opened for writing whose path holds <text> (a chromosome's dict-less tier or
# indel genotype store) waits there, once $STUB_MARK is made, until $STUB_RELEASE
# exists; STUB_WRITE_FAIL=<text>: that open fails, as on a full disk, or with the errno
# $STUB_WRITE_FAIL_ERRNO names (MemoryError: with that).
# STUB_SIGNAL_AT_RENAME=<path>: the first rename to <path>, or into it, is preceded by
# a signal to the process itself, SIGTERM or the one $STUB_SIGNAL_AT_RENAME_SIG names;
# STUB_RENAME_HANG=<path>: it waits instead, once $STUB_MARK is made, until
# $STUB_RELEASE exists. STUB_PROBE_FULL=<errno name>: the file build-index-only makes
# and removes in each directory of an index it built cannot be made, with that errno,
# as with the disk's inode table or a file-count quota full.
SHIM = r'''
import builtins, errno, os, signal, time

_mount = os.environ.get("STUB_MOUNT")
if _mount:
    _mount = os.path.realpath(_mount)

    def _inside(path):
        parent = os.path.realpath(os.path.dirname(os.path.abspath(path)))
        return parent == _mount or parent.startswith(_mount + os.sep)

    def _within_mounts(rename):
        def renamed(src, dst, *args, **kwargs):
            if _inside(src) != _inside(dst):
                raise OSError(errno.EXDEV, os.strerror(errno.EXDEV), src, None, dst)
            return rename(src, dst, *args, **kwargs)
        return renamed

    os.rename = _within_mounts(os.rename)
    os.replace = _within_mounts(os.replace)

if os.environ.get("STUB_SIGNAL_AT_CLEANUP"):
    _signal, _sent = signal.signal, []

    def _signal_at_cleanup(sig, handler):
        current = signal.getsignal(signal.SIGTERM)
        if not _sent and sig == signal.SIGTERM and callable(current):  # the build's own
            _sent.append(1)
            if os.environ["STUB_SIGNAL_AT_CLEANUP"] == "before":
                current(signal.SIGTERM, None)  # what a SIGTERM that comes now runs
            else:
                old = _signal(sig, handler)
                os.kill(os.getpid(), signal.SIGTERM)
                return old
        return _signal(sig, handler)

    signal.signal = _signal_at_cleanup

if os.environ.get("STUB_OVERLAY_HANG"):
    _open = builtins.open

    def open(file, mode="r", *args, **kwargs):
        if ("w" in mode and str(file).endswith(".fa.tmp")
                and os.path.basename(os.path.dirname(str(file))).startswith(".overlay-")):
            _open(os.environ["STUB_MARK"], "w").close()
            while not os.path.exists(os.environ["STUB_RELEASE"]):
                time.sleep(0.02)
        return _open(file, mode, *args, **kwargs)

    builtins.open = open

if os.environ.get("STUB_SLOW_CLEANUP"):
    import shutil
    _rmtree = shutil.rmtree

    def rmtree(path, *args, **kwargs):
        parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
        if parent.startswith(".build-index-only-"):
            builtins.open(os.environ["STUB_SLOW_CLEANUP"], "w").close()
            while not os.path.exists(os.environ["STUB_RELEASE"]):
                time.sleep(0.02)
        return _rmtree(path, *args, **kwargs)

    shutil.rmtree = rmtree

if os.environ.get("STUB_THREAD"):
    import threading
    threading.Thread(target=time.sleep, args=(3600,), daemon=True).start()

if os.environ.get("STUB_WRITE_HANG") or os.environ.get("STUB_WRITE_FAIL"):
    _open_for = builtins.open

    def open(file, mode="r", *args, **kwargs):
        if "w" in mode and os.environ.get("STUB_WRITE_HANG", "\0") in str(file):
            _open_for(os.environ["STUB_MARK"], "w").close()
            while not os.path.exists(os.environ["STUB_RELEASE"]):
                time.sleep(0.02)
        if "w" in mode and os.environ.get("STUB_WRITE_FAIL", "\0") in str(file):
            why = os.environ.get("STUB_WRITE_FAIL_ERRNO", "ENOSPC")
            if why == "MemoryError":
                raise MemoryError()
            raise OSError(getattr(errno, why), os.strerror(getattr(errno, why)), str(file))
        return _open_for(file, mode, *args, **kwargs)

    builtins.open = open

if os.environ.get("STUB_SIGNAL_AT_RENAME") or os.environ.get("STUB_RENAME_HANG"):
    _to = os.environ.get("STUB_SIGNAL_AT_RENAME") or os.environ["STUB_RENAME_HANG"]
    _stopped = []

    def _stop_before(rename):
        def renamed(src, dst, *args, **kwargs):
            if not _stopped and (str(dst) == _to or str(dst).startswith(_to + os.sep)):
                _stopped.append(1)
                if os.environ.get("STUB_SIGNAL_AT_RENAME"):
                    os.kill(os.getpid(), getattr(signal, "SIG" + os.environ.get(
                        "STUB_SIGNAL_AT_RENAME_SIG", "TERM")))
                else:
                    builtins.open(os.environ["STUB_MARK"], "w").close()
                    while not os.path.exists(os.environ["STUB_RELEASE"]):
                        time.sleep(0.02)
            return rename(src, dst, *args, **kwargs)
        return renamed

    os.rename = _stop_before(os.rename)

if os.environ.get("STUB_PROBE_FULL"):
    _os_open = os.open

    def _open_probe(path, flags, *args, **kwargs):
        if os.path.basename(str(path)) == ".build-index-only-probe":
            why = getattr(errno, os.environ["STUB_PROBE_FULL"])
            raise OSError(why, os.strerror(why), str(path))
        return _os_open(path, flags, *args, **kwargs)

    os.open = _open_probe
'''

SHARED = None  # the stub Bio package and a warmed matplotlib cache, for every run


def setUpModule():
    global SHARED
    SHARED = tempfile.mkdtemp(prefix="bio_atomic_shared_")
    os.makedirs(os.path.join(SHARED, "py", "Bio"))
    open(os.path.join(SHARED, "py", "Bio", "__init__.py"), "w").close()
    with open(os.path.join(SHARED, "py", "Bio", "Seq.py"), "w") as fh:
        fh.write("Seq = object\n")
    os.makedirs(os.path.join(SHARED, "shim"))
    with open(os.path.join(SHARED, "shim", "sitecustomize.py"), "w") as fh:
        fh.write(SHIM)
    # crisprme.py imports matplotlib, which builds its font cache (and says so on
    # stderr) the first time: do it here, once
    subprocess.run([sys.executable, "-c", "import matplotlib.pyplot"],
                   env=dict(os.environ, MPLCONFIGDIR=os.path.join(SHARED, "mpl")),
                   capture_output=True)


def tearDownModule():
    shutil.rmtree(SHARED, ignore_errors=True)


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def make_inputs(work):
    """A two-chromosome genome, a phased per-chromosome VCF dataset with a SNP and
    an insertion on each, and its samplesID listing."""
    for chrom, seq in GENOME.items():
        _write(os.path.join(work, "Genomes", "hg38", chrom + ".fa"), f">{chrom}\n{seq}\n")
        vcf = os.path.join(work, "VCFs", "ds", f"ds.{chrom}.vcf.gz")
        os.makedirs(os.path.dirname(vcf), exist_ok=True)
        with gzip.open(vcf, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                     f"{chrom}\t11\trs_{chrom}_11\tA\tG\t.\tPASS\tAF=0.25\tGT\t0|1\t0|0\n"
                     f"{chrom}\t26\trs_{chrom}_26\t{seq[25]}\t{seq[25]}T\t.\tPASS\tAF=0.25"
                     "\tGT\t0|0\t1|0\n")
    _write(os.path.join(work, "samplesIDs", "hg38_ds.samplesID.txt"),
           "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\n"
           "S1\tGBR\tEUR\tmale\nS2\tYRI\tAFR\tfemale\n")
    _write(os.path.join(work, "samples.txt"), "hg38_ds.samplesID.txt\n")


def make_stubs(bindir, pigz=True):
    for name, text in (("crispritz.py", STUB_ENGINE), ("pigz", STUB_PIGZ))[: 2 if pigz else 1]:
        path = os.path.join(bindir, name)
        _write(path, text.replace("#!PYTHON", "#!" + sys.executable, 1))
        os.chmod(path, 0o755)


def proc_start(pid):
    """Field 22 of /proc/<pid>/stat, the start time an owner record holds."""
    with open(f"/proc/{pid}/stat") as fh:
        return fh.read().rpartition(") ")[2].split()[19]


def proc_state(pid):
    """Field 3 of /proc/<pid>/stat: R, S, T (stopped), Z (a zombie)..."""
    with open(f"/proc/{pid}/stat") as fh:
        return fh.read().rpartition(") ")[2][:1]


def owner_record(owner, **fields):
    """An owner record of a build directory, as build-index-only writes it, for
    process owner on this host and PID namespace; fields replace (None drops) any."""
    record = {"host": os.uname().nodename, "pidns": os.readlink("/proc/self/ns/pid"),
              "pid": str(owner), "started": "2026-10-31T01:02:03Z"}
    try:
        record["pstart"] = proc_start(owner)
    except OSError:  # no such process
        record["pstart"] = "12345"
    record.update(fields)
    return "".join(f"{k}={v}\n" for k, v in record.items() if v is not None)


def tree(top, normalize=()):
    """{path: digest} of everything under top, a directory as "dir". A gzip file is
    taken by its content (the header holds a time), and every (path, label) in
    normalize is replaced by its label inside each file."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(top):
        dirnames.sort()
        for d in dirnames:
            out[os.path.relpath(os.path.join(dirpath, d), top) + "/"] = "dir"
        for f in sorted(filenames):
            p = os.path.join(dirpath, f)
            with open(p, "rb") as fh:
                data = fh.read()
            if f.endswith(".gz"):
                data = gzip.decompress(data)
            for path, label in normalize:
                data = data.replace(path.encode(), label.encode())
            out[os.path.relpath(p, top)] = hashlib.sha256(data).hexdigest()[:16]
    return out


class _Fixture(unittest.TestCase):
    def setUp(self):
        self.fresh()

    def fresh(self):
        """A new working directory with the inputs, and the stubs."""
        self.root = tempfile.mkdtemp(prefix="bio_atomic_")
        self.addCleanup(shutil.rmtree, self.root, True)
        self.work = os.path.join(self.root, "work")
        self.bin = os.path.join(self.root, "bin")
        os.makedirs(self.bin)
        make_stubs(self.bin)
        make_inputs(self.work)
        self.log = os.path.join(self.root, "calls.log")
        self.mark = os.path.join(self.root, "mark")
        self.release = os.path.join(self.root, "release")
        self.runs = 0

    def path(self, *parts):
        return os.path.join(self.work, *parts)

    def lib(self, name):
        return self.path("genome_library", name)

    def start(self, vcf=True, name=None, keep_hup=False, pigz=True, fsize=None, shim=False,
              command=None, **stub):
        """Launch build-index-only in a session of its own (so its whole process
        group can be killed, as an OOM kill or a container stop would). Without
        pigz, no pigz is on PATH; fsize limits the size of the files it writes, as
        a full disk would; shim loads SHIM in each Python process of the run;
        command runs that crisprme.py command line instead."""
        args = ["build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam", PAM,
                "--bDNA", "1", "--bRNA", "1", "--thread", "2", "--path", self.work]
        if vcf:
            args += ["--vcf", self.path("VCFs", "ds"), "--samplesID", self.path("samples.txt")]
        if name:
            args += ["--name", name]
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("STUB_", "CRISPRME_"))}
        path = os.environ.get("PATH", "").split(os.pathsep)
        if not pigz:
            bindir = os.path.join(self.root, "bin_without_pigz")
            if not os.path.isdir(bindir):
                make_stubs(bindir, pigz=False)
            path = [bindir] + [d for d in path if not shutil.which("pigz", path=d)]
        else:
            path = [self.bin] + path
        pythonpath = [os.path.join(SHARED, "py")] + [os.path.join(SHARED, "shim")] * shim
        env.update(PATH=os.pathsep.join(path), PYTHONPATH=os.pathsep.join(pythonpath),
                   MPLCONFIGDIR=os.path.join(SHARED, "mpl"),
                   STUB_LOG=self.log, STUB_MARK=self.mark, STUB_RELEASE=self.release,
                   **stub)
        for f in (self.mark, self.release):
            if os.path.exists(f):
                os.remove(f)
        self.runs += 1
        self.out = os.path.join(self.root, f"stdout{self.runs}.txt")
        self.err = os.path.join(self.root, f"stderr{self.runs}.txt")

        def signals():
            for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
                signal.signal(sig, signal.SIG_DFL)
            if keep_hup:  # as under nohup
                signal.signal(signal.SIGHUP, signal.SIG_IGN)
            if fsize:
                resource.setrlimit(resource.RLIMIT_FSIZE,
                                   (fsize, resource.getrlimit(resource.RLIMIT_FSIZE)[1]))

        with open(self.out, "w") as out, open(self.err, "w") as err:
            proc = subprocess.Popen([sys.executable, CRISPRME, *(command or args)],
                                    cwd=self.work, env=env, stdout=out, stderr=err,
                                    start_new_session=True, preexec_fn=signals)
        self.addCleanup(self._reap, proc)
        return proc

    @staticmethod
    def _reap(proc):
        try:
            os.killpg(proc.pid, signal.SIGKILL)  # whatever of the group is left
        except (ProcessLookupError, PermissionError):
            pass
        proc.wait()

    def finish(self, proc, timeout=120):
        proc.wait(timeout=timeout)
        return proc.returncode, self.stdout(), self.stderr()

    def run_bio(self, **kw):
        return self.finish(self.start(**kw))

    def run_search(self):
        """complete-search --vcf here, as a user runs it next (CFD only, as no
        CRISPR-Bulge env is here), into an output folder of its own."""
        shutil.rmtree(self.path("Results", "search"), ignore_errors=True)
        _write(self.path("guides.txt"), GENOME["chr1"][:20] + "NNN\n")
        _write(self.path("vcfs.txt"), "ds\n")
        command = ["complete-search", "--genome", self.path("Genomes", "hg38"),
                   "--vcf", self.path("vcfs.txt"), "--samplesID", self.path("samples.txt"),
                   "--guide", self.path("guides.txt"), "--pam", PAM, "--mm", "4",
                   "--bDNA", "1", "--bRNA", "1", "--output", "search", "--thread", "2"]
        return self.finish(self.start(command=command, CRISPRME_SKIP_SCORER_ENV="1"))

    def run_web_search(self):
        """The search the web UI launches (pages/main_page.py): the search script run
        straight, not complete-search, with the job directory, environment and
        arguments it gives it, into a job directory of its own. Returns its exit
        status, stdout and log_error.txt."""
        job = self.path("Results", "web")
        shutil.rmtree(job, ignore_errors=True)
        for name, text in ((".list_vcfs.txt", "ds\n"), (".samplesID.txt", "hg38_ds.samplesID.txt\n"),
                           (".guides.txt", GENOME["chr1"][:20] + "NNN\n"), ("vuoto.txt", "")):
            _write(os.path.join(job, name), text)
        empty = os.path.join(job, "vuoto.txt")
        args = [self.path("Genomes", "hg38"), os.path.join(job, ".list_vcfs.txt"),
                os.path.join(job, ".guides.txt"), PAM, empty, os.path.join(job, ".samplesID.txt"),
                2, 4, 1, 1, 3, job, HERE, 4, self.work + os.sep, empty, "_", 1, 0, "none",
                "mm+bulges", "mm+bulges,mm", "False", "PASS,.", "_", 6]
        env = {k: v for k, v in os.environ.items() if not k.startswith(("STUB_", "CRISPRME_"))}
        env.update(PATH=os.pathsep.join([self.bin, os.environ.get("PATH", "")]),
                   PYTHONPATH=os.path.join(SHARED, "py"), MPLCONFIGDIR=os.path.join(SHARED, "mpl"),
                   STUB_LOG=self.log, CRISPRME_FAST_MODE="1", CRISPRME_EMIT_ALT_ALIGNMENTS="0",
                   CRISPRME_REPORT_FULL="False")
        with open(os.path.join(job, "log_verbose.txt"), "w") as out, \
                open(os.path.join(job, "log_error.txt"), "w") as err:
            rc = subprocess.call(
                ["bash", os.path.join(HERE, "submit_job_automated_new_multiple_vcfs.sh"),
                 *map(str, args)], cwd=self.work, env=env, stdout=out, stderr=err, timeout=120)
        with open(os.path.join(job, "log_verbose.txt")) as out, \
                open(os.path.join(job, "log_error.txt")) as err:
            return rc, out.read(), err.read()

    def assert_searched(self):
        """complete-search goes past its checks of the variant index (the stub
        engine then fails it)."""
        rc, out, err = self.run_search()
        self.assertNotIn(" is unfinished: ", err)
        self.assertNotIn(" is being written", err)
        self.assertIn("Launching job ", out)

    def outputs(self):
        """What build-index-only builds, as tree() has it."""
        return {top: tree(self.path(top), ((self.work, "<W>"),))
                for top in ("Dictionaries", "Genomes", "genome_library")}

    def stdout(self):
        with open(self.out) as fh:
            return fh.read()

    def stderr(self):
        with open(self.err) as fh:
            return fh.read()

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return fh.read().splitlines()

    def wait_for_mark(self, proc, timeout=60):
        deadline = time.monotonic() + timeout
        while not os.path.exists(self.mark):
            if proc.poll() is not None:
                self.fail(f"build-index-only ended (rc {proc.returncode}) before the "
                          f"stub got there:\n{self.stdout()}\n{self.stderr()}")
            if time.monotonic() > deadline:
                self.fail("the stub never got there")
            time.sleep(0.02)

    def kill_group(self, proc):
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait(timeout=30)

    def leftovers(self):
        """Build directories left anywhere in the working directory."""
        found = []
        for dirpath, dirnames, _ in os.walk(self.work):
            found += [os.path.join(dirpath, d) for d in dirnames
                      if d.startswith(".build-index-only-")]
        return found

    def bins(self, name):
        top = self.lib(name)
        return sorted(os.path.relpath(os.path.join(d, f), top)
                      for d, _, fs in os.walk(top) for f in fs if f.endswith(".bin"))

    def lock(self, output, record=True):
        os.makedirs(output + ".lock")
        if record:
            _write(os.path.join(output + ".lock", "owner"),
                   "host=otherhost\npidns=pid:[4026531836]\npid=4242\n"
                   "started=2026-10-31T00:00:00Z\npstart=12345\n")

    def kit(self, pad=0):
        """An enrichment in place that came another way (a build-inputs kit, a run of
        an earlier release): its fake-indel genome as add-variants wrote it, without
        the SNP codes (chr2's padded with pad N), and no INDELS index yet."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        shutil.rmtree(self.lib(INDELS))
        for chrom, seq in GENOME.items():
            _write(self.path("Genomes", "hg38+ds_INDELS", f"fake_ds_{chrom}", f"fake{chrom}.fa"),
                   f">fake{chrom}\n{seq[:26]}T{seq[26:60]}{'N' * pad * (chrom == 'chr2')}\n")
        os.remove(self.log)


# What each index looks like once built, and the step of the stub's that builds its
# second chromosome (so stopping there leaves the first one built).
WHOLE = {
    REF: (False, "index-genome hg38 chr2.fa", ["NRG_chr1_1.bin", "NRG_chr2_1.bin"]),
    VAR: (True, "index-genome hg38+ds chr2.enriched.fa", ["NRG_chr1_1.bin", "NRG_chr2_1.bin"]),
    INDELS: (True, "NRG_2_fakechr2 fakechr2.fa",
             ["NRG_2_fakechr1/NRG_fakechr1_1.bin", "NRG_2_fakechr2/NRG_fakechr2_1.bin"]),
}


class AnIndexIsPutInPlaceWholeOrNotAtAll(_Fixture):
    """(a)-(c): a build killed partway, or one that fails, leaves no index."""

    def test_a_build_killed_partway_leaves_no_index_and_the_next_run_builds_it_whole(self):
        for name, (vcf, chr2, whole) in WHOLE.items():
            with self.subTest(index=name):
                self.fresh()
                proc = self.start(vcf=vcf, STUB_HANG=chr2)
                self.wait_for_mark(proc)
                self.kill_group(proc)
                self.assertFalse(os.path.lexists(self.lib(name)),
                                 f"a partial {name} was left in place")
                # a SIGKILL leaves the hidden build directory, which nothing reads, and
                # which the next build removes, as the run that left it is gone
                left = self.leftovers()
                self.assertEqual(len(left), 1)
                rc, out, err = self.run_bio(vcf=vcf)
                self.assertEqual(rc, 0, err)
                self.assertEqual(self.bins(name), whole)
                self.assertIn(f"Removing {left[0]}: the build-index-only run that made it "
                              f"(pid {proc.pid} on {os.uname().nodename}, started ", out)
                self.assertEqual(self.leftovers(), [])

    def test_an_index_is_built_in_a_hidden_directory_of_genome_library(self):
        """Where the index is renamed from: no rename then crosses a mount point, and
        the index is made with what genome_library gives a new directory."""
        if not os.path.exists("/proc/self/stat"):
            self.skipTest("no /proc here")
        proc = self.start(vcf=False, STUB_HANG=WHOLE[REF][1])
        self.wait_for_mark(proc)
        left = self.leftovers()
        self.assertEqual([os.path.dirname(d) for d in left], [self.path("genome_library")])
        # with the owner record that lets a later run tell whether it is still in use
        with open(os.path.join(left[0], "owner")) as fh:
            record = dict(line.split("=", 1) for line in fh.read().splitlines())
        self.assertEqual(sorted(record), ["host", "pid", "pidns", "pstart", "started"])
        self.assertEqual((record["host"], record["pid"], record["pstart"]),
                         (os.uname().nodename, str(proc.pid), proc_start(proc.pid)))
        self.assertEqual(record["pidns"], os.readlink("/proc/self/ns/pid"))
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.bins(REF), WHOLE[REF][2])
        self.assertEqual(self.leftovers(), [])

    def test_a_genome_library_mounted_on_its_own_gets_its_indexes(self):
        """genome_library as a bind mount of its own (docker -v, apptainer --bind) on
        the working directory's file system: a rename between the two fails with
        EXDEV though their st_dev is the same. SHIM makes os.rename and os.replace do
        that here, as a test cannot mount."""
        os.makedirs(self.path("genome_library"))
        rc, out, err = self.run_bio(shim=True, STUB_MOUNT=self.path("genome_library"))
        self.assertEqual(rc, 0, err)
        for name, (_, _, whole) in WHOLE.items():
            self.assertEqual(self.bins(name), whole)
        self.assertEqual(self.leftovers(), [])

    def test_the_indexes_take_the_setgid_group_of_genome_library(self):
        """A shared genome_library, setgid to its lab's group: each index and its
        files take that group, and the setgid bit, as they do built in place."""
        if not sys.platform.startswith("linux"):
            self.skipTest("setgid directories pass their group on as on Linux only")
        if os.geteuid() == 0:
            gid = 4242
        else:
            gid = next((g for g in os.getgroups() if g != os.getegid()), None)
            if gid is None:
                self.skipTest("no group to give genome_library but this user's own")
        lib = self.path("genome_library")
        os.makedirs(lib)
        os.chown(lib, -1, gid)
        os.chmod(lib, 0o2775)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        for name in WHOLE:
            for top, dirs, files in os.walk(self.lib(name)):
                self.assertTrue(os.stat(top).st_mode & stat.S_ISGID, top)
                for entry in [top] + [os.path.join(top, f) for f in files]:
                    self.assertEqual(os.stat(entry).st_gid, gid, entry)

    def test_a_failed_build_exits_non_zero_and_leaves_neither_index_nor_build_directory(self):
        for name, (vcf, chr2, _) in WHOLE.items():
            with self.subTest(index=name):
                self.fresh()
                rc, out, err = self.run_bio(vcf=vcf, STUB_FAIL=chr2)
                self.assertNotEqual(rc, 0)
                self.assertIn("failed (expected " + self.lib(name) + ")", err)
                if name == INDELS:  # pool_index_indels.py no longer passes it over
                    self.assertIn("Indexing failed for INDELs in chr2\n", out)
                self.assertFalse(os.path.lexists(self.lib(name)))
                self.assertEqual(self.leftovers(), [])

    def test_a_build_the_engine_left_cut_short_is_not_put_in_place(self):
        """buildTST does not check its writes: on a full disk it left a .bin empty or
        cut short and still exited 0 (the stub does so at $STUB_CUT). That index went
        into place ('Index built', 'Variant index ready'), every later run printed
        'already present' and kept it, and every search failed on it ('numLeaves=',
        'numNodes=0'), until it was removed by hand. Each .bin is now checked first."""
        for name, (vcf, chr2, whole) in WHOLE.items():
            for cut, size in (("empty", {"STUB_CUT_TO": "0"}), ("cut short", {})):
                with self.subTest(index=name, bin=cut):
                    self.fresh()
                    rc, out, err = self.run_bio(vcf=vcf, STUB_CUT=chr2, **size)
                    self.assertEqual(rc, 1, out + err)
                    # (for the INDELS index, the file of the chromosome it is)
                    self.assertIn(f" indexing failed: the engine exited 0, but {whole[1]} is "
                                  f"{cut}", err)
                    self.assertIn(", as when the disk fills up: nothing was put in place "
                                  f"(expected {self.lib(name)})\n", err)
                    self.assertFalse(os.path.lexists(self.lib(name)))
                    self.assertEqual(self.leftovers(), [])
                    # and once there is room, the next run builds it whole
                    rc, out, err = self.run_bio(vcf=vcf)
                    self.assertEqual(rc, 0, err)
                    self.assertEqual(self.bins(name), whole)

    def test_a_build_the_engine_left_without_a_bin_is_not_put_in_place(self):
        """buildTST does not check that a .bin opened either: with the disk's inode
        table or a file-count quota full, it wrote no .bin for the contigs left and
        exited 0, every .bin there whole (the stub does so at $STUB_NO_BIN). That index
        went into place, every later run kept it, and every search lacked the
        off-targets of those contigs, with no error. A contig that holds a target and
        has no .bin now fails the build. Where a contig has no .bin and cannot be read
        as buildTST reads it, so does a directory of the index that cannot take a file
        more as the engine ends. A contig with no target gets no .bin, and an index
        whole goes into place as before, though its build filled the inode table."""
        failure = {REF: "Reference genome", VAR: "Variant genome", INDELS: "Indels"}
        fasta = {REF: ("Genomes", "hg38", "chr2.fa"), VAR: ("Genomes", "hg38+ds", "chr2.enriched.fa"),
                 INDELS: ("Genomes", "hg38+ds_INDELS", "fake_ds_chr2", "fakechr2.fa")}
        contig = {REF: "chr2", VAR: "chr2", INDELS: "fakechr2"}
        for name, (vcf, chr2, whole) in WHOLE.items():
            for how, stub in (("no .bin", dict(STUB_NO_BIN=chr2)),
                              ("the inode table full as it ends", dict(shim=True, STUB_PROBE_FULL="ENOSPC"))):
                with self.subTest(index=name, how=how):
                    self.fresh()
                    if name != REF:  # (only that index to build)
                        rc, out, err = self.run_bio()
                        self.assertEqual(rc, 0, err)
                        shutil.rmtree(self.lib(name))
                    rc, out, err = self.run_bio(vcf=vcf, **stub)
                    if how != "no .bin":  # (every contig has its .bin)
                        self.assertEqual(rc, 0, err)
                        self.assertEqual(self.bins(name), whole)
                        continue
                    self.assertEqual(rc, 1, out + err)
                    self.assertEqual(err, f"Error: {failure[name]} indexing failed: the engine "
                                          f"exited 0, but wrote no .bin for {contig[name]} "
                                          f"({self.path(*fasta[name])}), which holds a target, as "
                                          "when the disk's inode table or a file-count quota is "
                                          f"full: nothing was put in place (expected {self.lib(name)})\n")
                    self.assertFalse(os.path.lexists(self.lib(name)))
                    self.assertEqual(self.leftovers(), [])
                    rc, out, err = self.run_bio(vcf=vcf)
                    self.assertEqual(rc, 0, err)
                    self.assertEqual(self.bins(name), whole)
        # a contig it cannot read as buildTST does (a byte no IUPAC base), with no .bin
        for full in (None, "ENOSPC", "EDQUOT"):
            with self.subTest(unread=True, full=full):
                self.fresh()
                _write(self.path("Genomes", "hg38", "chrUn.fa"), ">chrUn\n" + GENOME["chr1"] + "*\n")
                rc, out, err = self.run_bio(vcf=False, STUB_NO_BIN="chrUn.fa",
                                            **(dict(shim=True, STUB_PROBE_FULL=full) if full else {}))
                if not full:  # (as before: not told from one with no target)
                    self.assertEqual(rc, 0, err)
                    self.assertEqual(self.bins(REF), WHOLE[REF][2])
                    continue
                why = getattr(errno, full)
                self.assertEqual(rc, 1, out + err)
                self.assertTrue(err.startswith("Error: Reference genome indexing failed: the engine "
                                               "exited 0, but the disk or its inode table is full "
                                               f"([Errno {why}] {os.strerror(why)}: "), err)
                self.assertTrue(err.endswith(".build-index-only-probe'), so it may have left out "
                                             "the .bin of a contig: nothing was put in place "
                                             f"(expected {self.lib(REF)})\n"), err)
                self.assertFalse(os.path.lexists(self.lib(REF)))
                self.assertEqual(self.leftovers(), [])
        # a contig with no target (N alone) gets no .bin, as before
        self.fresh()
        _write(self.path("Genomes", "hg38", "chr3.fa"), ">chr3\n" + "N" * 120 + "\n")
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.bins(REF), WHOLE[REF][2])

    def test_a_genome_library_on_another_file_system_still_gets_its_index(self):
        """No rename crosses file systems: the build runs beside the indexes."""
        if not os.access("/dev/shm", os.W_OK):
            self.skipTest("no /dev/shm here")
        elsewhere = tempfile.mkdtemp(prefix="bio_atomic_lib_", dir="/dev/shm")
        self.addCleanup(shutil.rmtree, elsewhere, True)
        if os.stat(elsewhere).st_dev == os.stat(self.work).st_dev:
            self.skipTest("/dev/shm is on the same file system here")
        os.symlink(elsewhere, self.path("genome_library"))
        rc, out, err = self.run_bio(vcf=False)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.bins(REF), WHOLE[REF][2])
        self.assertEqual(os.listdir(elsewhere), [REF])
        self.assertEqual(self.leftovers(), [])


class AnEnrichmentIsDoneWholeOrRedone(_Fixture):
    """(d): nothing of an enrichment is in place before all of it is: the
    dictionaries are compressed in its build directory, and the enriched genome is
    moved into place last."""

    def enrichment(self):
        return (self.path("Genomes", "hg38+ds"), self.path("Genomes", "hg38+ds_INDELS"),
                self.path("Dictionaries", "dictionaries_ds"),
                self.path("Dictionaries", "log_indels_ds"))

    def test_an_enrichment_killed_before_its_end_leaves_nothing_of_it_and_is_redone(self):
        enriched = self.path("Genomes", "hg38+ds")
        proc = self.start(STUB_PIGZ_HANG="1")  # add-variants is done; compressing
        self.wait_for_mark(proc)
        self.kill_group(proc)
        for output in self.enrichment():
            self.assertFalse(os.path.lexists(output), f"{output} was left in place")
        # the dictionaries alone would pass complete-search's check for a prebuilt
        # variant index, and its own enrichment would then fail on the fake-indel
        # genomes in its way: it refuses up front instead
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn("Error: No prebuilt variant index found for: ds.\n", err)
        self.assertIn("crisprme.py build-index-only --genome ", err)
        left = self.leftovers()
        self.assertEqual([os.path.dirname(d) for d in left], [self.path("Genomes")])
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"Removing {left[0]}: the build-index-only run that made it "
                      f"(pid {proc.pid} on {os.uname().nodename}, started ", out)
        self.assertEqual(self.leftovers(), [])
        self.assertIn("Enriching hg38 with ds (add-variants)...", out)
        self.assertEqual(sum(c.startswith("add-variants ") for c in self.calls()), 2)
        self.assertEqual(sorted(os.listdir(enriched)), ["chr1.enriched.fa", "chr2.enriched.fa"])
        self.assertEqual(sorted(os.listdir(self.path("Dictionaries", "dictionaries_ds"))),
                         ["my_dict_chr1.json.gz", "my_dict_chr2.json.gz"])
        self.assertEqual(sorted(os.listdir(self.path("Dictionaries", "log_indels_ds"))),
                         ["logchr1.txt.gz", "logchr2.txt.gz"])
        self.assertEqual(sorted(os.listdir(self.path("Genomes", "hg38+ds_INDELS"))),
                         ["fake_ds_chr1", "fake_ds_chr2"])
        for name, (_, _, whole) in WHOLE.items():
            self.assertEqual(self.bins(name), whole)

    def test_a_failed_compression_fails_the_run_and_moves_nothing_into_place(self):
        """pigz failing partway, or the Python gzip that stands in for it without
        pigz failing to write, as on a full disk."""
        for compressor, kw in (("pigz", dict(STUB_PIGZ_FAIL="1")),
                               ("gzip", dict(pigz=False, STUB_GZ_BLOCKED="1"))):
            with self.subTest(compressor=compressor):
                self.fresh()
                rc, out, err = self.run_bio(**kw)
                self.assertEqual(rc, 1, out + err)
                self.assertIn("Error: Compressing the ds dictionaries failed (", err)
                self.assertIn("): nothing of the enrichment was put in place.\n", err)
                for output in self.enrichment():
                    self.assertFalse(os.path.lexists(output), f"{output} was put in place")
                self.assertEqual(self.leftovers(), [])
                self.assertEqual([c.split()[:2] for c in self.calls()],
                                 [["index-genome", "hg38"],
                                  ["add-variants", self.path("VCFs", "ds") + "/"]])

    def test_a_signal_as_it_goes_into_place_waits_for_its_last_rename(self):
        """Its moves take moments, but a Ctrl-C, SIGTERM or SIGHUP among them (SHIM
        sends one as the dictionaries go into place) left the dictionaries in place
        without the enriched genome: a search passed its check for them, enriched the
        genome itself, and failed on the fake-indel genomes in its way, leaving its own
        dictionaries, uncompressed, which every later run kept. The signal now waits
        for the last rename, then ends the run by it."""
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            with self.subTest(signal=sig.name):
                self.fresh()
                rc, out, err = self.run_bio(
                    shim=True, STUB_SIGNAL_AT_RENAME=self.path("Dictionaries", "dictionaries_ds"),
                    STUB_SIGNAL_AT_RENAME_SIG=sig.name[3:])
                self.assertEqual(rc, -sig, err)
                for output in self.enrichment():
                    self.assertTrue(os.path.isdir(output), f"{output} is not in place")
                self.assertEqual(sorted(os.listdir(self.path("Dictionaries", "dictionaries_ds"))),
                                 ["my_dict_chr1.json.gz", "my_dict_chr2.json.gz"])
                self.assertEqual(sorted(os.listdir(self.path("Dictionaries", "log_indels_ds"))),
                                 ["logchr1.txt.gz", "logchr2.txt.gz"])
                self.assertEqual(self.leftovers(), [])
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                self.assertIn(f"Enriched genome already present: {self.path('Genomes', 'hg38+ds')}\n",
                              out)

    def test_a_kill_between_its_renames_leaves_nothing_a_later_run_uses(self):
        """A SIGKILL can still land between two of them (SHIM stops the run before the
        enriched genome's, the dictionaries' or the indel logs' rename). What is in
        place then passes neither a search's check, which needs both sets of
        dictionaries, nor the next build-index-only's: the enriched genome goes into
        place marked unfinished until the rest is, so that one enriches again."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = self.outputs()
        for where in (("Genomes", "hg38+ds"), ("Dictionaries", "dictionaries_ds"),
                      ("Dictionaries", "log_indels_ds")):
            with self.subTest(stopped_before=where[1]):
                self.fresh()
                proc = self.start(shim=True, STUB_RENAME_HANG=self.path(*where))
                self.wait_for_mark(proc)
                self.kill_group(proc)
                self.assertFalse(os.path.lexists(self.path(*where)))
                # refused up front: by then the tiers the run was to write next are
                # marked unfinished too (AStoppedRunsTiersAreNotSearched)
                rc, out, err = self.run_search()
                self.assertEqual(rc, 1, out + err)
                self.assertIn(f"Error: {self.path('Dictionaries', 'registry_ds')} is unfinished: ",
                              err)
                self.assertNotIn("Launching job ", out)
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                self.assertIn("Enriching hg38 with ds (add-variants)...", out)
                self.assertEqual(self.leftovers(), [])
                self.assertEqual(self.outputs(), clean)

    def test_a_kill_while_its_dictionaries_are_copied_to_a_mount_of_their_own_leaves_none_of_it(self):
        """Dictionaries mounted on its own (docker -v, apptainer --bind; SHIM makes a
        rename into it fail with EXDEV): the dictionaries are copied there, which takes
        minutes at scale. A run killed in that copy left the fake-indel genomes and
        part of the dictionaries in place without the enriched genome, which let a
        search start and enrich the genome itself. Each set is now copied whole into a
        hidden directory of Dictionaries first, and renamed into place with the rest."""
        def mount():
            os.makedirs(self.path("Dictionaries"), exist_ok=True)
            return dict(shim=True, STUB_MOUNT=self.path("Dictionaries"))

        rc, out, err = self.run_bio(**mount())
        self.assertEqual(rc, 0, err)
        clean = self.outputs()
        self.fresh()
        proc = self.start(STUB_WRITE_HANG="dictionaries_ds/my_dict_chr2.json.gz", **mount())
        self.wait_for_mark(proc)
        self.kill_group(proc)
        for output in self.enrichment():
            self.assertFalse(os.path.lexists(output), f"{output} was left in place")
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn("Error: No prebuilt variant index found for: ds.\n", err)
        left = self.leftovers()
        self.assertEqual(sorted(os.path.dirname(d) for d in left),
                         [self.path("Dictionaries"), self.path("Genomes")])
        rc, out, err = self.run_bio(**mount())
        self.assertEqual(rc, 0, err)
        for d in left:
            self.assertIn(f"Removing {d}: the build-index-only run that made it (pid {proc.pid} "
                          f"on {os.uname().nodename}, started ", out)
        self.assertEqual(self.leftovers(), [])
        self.assertEqual(self.outputs(), clean)

    def test_an_enriched_genome_in_place_without_its_dictionaries_is_used_as_before(self):
        """A build-inputs kit ships the enriched genome and the fake-indel genome only,
        for build-index-only to index a new PAM from: it is taken for a finished
        enrichment, as before (only an enriched genome still marked unfinished is
        not)."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        for d in ("genome_library/" + VAR, "genome_library/" + INDELS, "Dictionaries"):
            shutil.rmtree(self.path(d))
        os.remove(self.log)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"Enriched genome already present: {self.path('Genomes', 'hg38+ds')}\n", out)
        self.assertEqual([c for c in self.calls() if c.startswith("add-variants ")], [])
        self.assertEqual(self.bins(VAR), WHOLE[VAR][2])
    """An index or enriched genome another run put in place while this one built it
    is kept, and this run's build discarded (a rename onto it would fail, and a move
    would nest this one inside it)."""

    def test_an_index_another_run_put_in_place_meanwhile_is_kept(self):
        proc = self.start(vcf=False, STUB_HANG=WHOLE[REF][1])
        self.wait_for_mark(proc)
        _write(os.path.join(self.lib(REF), "NRG_theirs_1.bin"), "theirs\n")
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 0, err)
        self.assertIn(f"{self.lib(REF)} appeared while this run was building it (another "
                      "run built it): keeping that one", out)
        self.assertEqual(self.bins(REF), ["NRG_theirs_1.bin"])
        self.assertEqual(self.leftovers(), [])

    def test_an_enriched_genome_another_run_put_in_place_meanwhile_is_kept(self):
        enriched = self.path("Genomes", "hg38+ds")
        proc = self.start(STUB_HANG="add-variants ds chr2.fa")
        self.wait_for_mark(proc)
        for chrom, seq in GENOME.items():
            _write(os.path.join(enriched, chrom + ".enriched.fa"), f">{chrom}\n{seq}\n")
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 0, err)
        self.assertIn(f"{enriched} appeared while this run was enriching it (another run "
                      "built it): keeping that one", out)
        self.assertEqual(sorted(os.listdir(enriched)), ["chr1.enriched.fa", "chr2.enriched.fa"])
        with open(os.path.join(enriched, "chr1.enriched.fa")) as fh:
            self.assertNotIn("R", fh.read().split("\n", 1)[1], "this run's enrichment went in")
        self.assertEqual([d for d in os.listdir(self.path("Genomes")) if d.startswith("run_")], [])


class SignalsRemoveTheBuildDirectory(_Fixture):
    """Ctrl-C, SIGTERM and SIGHUP end the run by that signal, as before, and leave
    neither the index nor its build directory; a hangup it ignores (nohup) stays
    ignored."""

    def test_ctrl_c_sigterm_and_a_hangup_leave_nothing(self):
        for sig, group in ((signal.SIGINT, True), (signal.SIGTERM, False),
                           (signal.SIGTERM, True), (signal.SIGHUP, False)):
            with self.subTest(signal=sig.name, to_the_group=group):
                self.fresh()
                proc = self.start(vcf=False, STUB_HANG=WHOLE[REF][1])
                self.wait_for_mark(proc)
                if group:
                    os.killpg(proc.pid, sig)
                else:
                    os.kill(proc.pid, sig)
                proc.wait(timeout=60)
                self.assertEqual(proc.returncode, -sig, self.stderr())
                self.assertFalse(os.path.lexists(self.lib(REF)))
                self.assertEqual(self.leftovers(), [])

    def test_a_sigterm_just_as_the_build_ends_still_ends_the_run_by_it(self):
        """One that comes once the index is in place, as the cleanup begins, before it
        holds signals back or just after (SHIM sends it there): the run removes its
        build directory and ends by that signal, with no traceback."""
        for when in ("before", "after"):
            with self.subTest(when=when):
                self.fresh()
                rc, out, err = self.run_bio(vcf=False, shim=True, STUB_SIGNAL_AT_CLEANUP=when)
                self.assertEqual(rc, -signal.SIGTERM, err)
                self.assertNotIn("Traceback", err)
                self.assertEqual(self.leftovers(), [])
                self.assertEqual(self.bins(REF), WHOLE[REF][2])

    def test_a_second_signal_while_it_cleans_up_waits_for_the_cleanup(self):
        """Ctrl-C pressed twice, or scancel and then Ctrl-C: the second signal comes
        while the cleanup removes the partial index, which takes seconds for a large
        one (SHIM holds the cleanup there), and the kernel gives it to another thread
        than the one cleaning up, as to OpenBLAS's in a run that imported numpy (SHIM
        starts one). It does not cut the cleanup short: the run removes its build
        directory, then ends by that signal."""
        cases = ((signal.SIGINT, signal.SIGINT), (signal.SIGTERM, signal.SIGINT),
                 (signal.SIGINT, signal.SIGTERM))
        for first, second in cases:
            with self.subTest(first=first.name, second=second.name):
                self.fresh()
                cleaning = os.path.join(self.root, "cleaning")
                proc = self.start(vcf=False, shim=True, STUB_HANG=WHOLE[REF][1],
                                  STUB_SLOW_CLEANUP=cleaning, STUB_THREAD="1")
                self.wait_for_mark(proc)
                os.killpg(proc.pid, first)
                deadline = time.monotonic() + 60
                while not os.path.exists(cleaning):
                    self.assertIsNone(proc.poll(), f"the run ended before its cleanup:\n"
                                                   f"{self.stderr()}")
                    self.assertLess(time.monotonic(), deadline, "the cleanup never began")
                    time.sleep(0.01)
                os.killpg(proc.pid, second)
                time.sleep(0.3)  # (one that cuts the cleanup short does so by then)
                open(self.release, "w").close()
                proc.wait(timeout=60)
                self.assertEqual(proc.returncode, -second, self.stderr())
                self.assertFalse(os.path.lexists(self.lib(REF)))
                self.assertEqual(self.leftovers(), [])

    def test_a_hangup_ignored_as_under_nohup_does_not_stop_the_build(self):
        proc = self.start(vcf=False, keep_hup=True, STUB_HANG=WHOLE[REF][1])
        self.wait_for_mark(proc)
        os.kill(proc.pid, signal.SIGHUP)
        time.sleep(0.3)
        self.assertIsNone(proc.poll(), "an ignored hangup stopped the build")
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.bins(REF), WHOLE[REF][2])


class ASearchsBuildLockStopsIt(_Fixture):
    """(e): while complete-search holds "<output>.lock", build-index-only neither
    takes that output as present nor builds it: it stops with an ERROR."""

    OUTPUTS = (("genome_library", REF), ("Genomes", "hg38+ds"),
               ("genome_library", VAR), ("genome_library", INDELS))

    def check_refused(self, output, before, owner=True, vcf=True):
        rc, out, err = self.run_bio(vcf=vcf)
        self.assertEqual(rc, 1, out + err)
        lock = output + ".lock"
        if owner:
            self.assertIn(f"Error: {output} is locked by the run that builds it: {lock} is "
                          "held by pid 4242 on otherhost (a run started "
                          "2026-10-31T00:00:00Z). Wait for that run to finish", err)
        else:
            self.assertIn(f"Error: {output} is locked by the run that builds it: {lock} "
                          "has no owner record", err)
        guarded = " '" + self.path("Genomes", "hg38+ds_INDELS") + "'" \
            if output.endswith("Genomes/hg38+ds") else ""
        # (not complete-search's takeover: an enrichment that makes has no SNP codes)
        self.assertIn(f"Wait for that run to finish, then run build-index-only again. If no "
                      f"run is building it, remove the lock and what it guards, then run "
                      f"build-index-only again: rm -rf '{lock}' '{output}'{guarded}\n", err)
        self.assertEqual(self.calls(), [], "it built something")
        self.assertEqual(tree(self.work), before, "it changed the working directory")

    def test_a_lock_beside_each_output_stops_the_run_before_it_builds_anything(self):
        for where, name in self.OUTPUTS:
            for present in (False, True):
                with self.subTest(output=name, present=present):
                    self.fresh()
                    if present:
                        rc, out, err = self.run_bio()
                        self.assertEqual(rc, 0, err)
                        os.remove(self.log)
                    output = self.path(where, name)
                    self.lock(output)
                    self.check_refused(output, tree(self.work))
        with self.subTest(output=REF, present=False, vcf=False):
            self.fresh()
            self.lock(self.lib(REF))
            self.check_refused(self.lib(REF), tree(self.work), vcf=False)

    def test_a_lock_without_an_owner_record_stops_it_too(self):
        output = self.lib(VAR)
        self.lock(output, record=False)
        self.check_refused(output, tree(self.work), owner=False)

    def test_a_lock_taken_while_an_index_is_built_discards_the_build(self):
        for name, (vcf, chr2, _) in WHOLE.items():
            with self.subTest(index=name):
                self.fresh()
                proc = self.start(vcf=vcf, STUB_HANG=chr2)
                self.wait_for_mark(proc)
                self.lock(self.lib(name))
                open(self.release, "w").close()
                rc, out, err = self.finish(proc)
                self.assertEqual(rc, 1, out + err)
                self.assertIn(f"Error: {self.lib(name)} is locked by the run that builds it", err)
                self.assertFalse(os.path.lexists(self.lib(name)))
                self.assertEqual(self.leftovers(), [])

    def test_a_lock_taken_after_the_run_began_stops_it_at_that_output(self):
        """Each step looks again: here the lock appears while the reference index is
        built, beside an output a later step would otherwise take as present."""
        for where, name in self.OUTPUTS[1:]:
            with self.subTest(output=name):
                self.fresh()
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                shutil.rmtree(self.lib(REF))
                os.remove(self.log)
                proc = self.start(STUB_HANG=WHOLE[REF][1])
                self.wait_for_mark(proc)
                output = self.path(where, name)
                self.lock(output)
                before = tree(output)
                open(self.release, "w").close()
                rc, out, err = self.finish(proc)
                self.assertEqual(rc, 1, out + err)
                self.assertIn(f"Error: {output} is locked by the run that builds it", err)
                self.assertNotIn("already present: " + output, out)
                self.assertEqual(tree(output), before)
                self.assertEqual([c.split()[1] for c in self.calls()], ["hg38"])

    def test_a_lock_taken_while_the_genome_is_enriched_moves_nothing_into_place(self):
        enriched = self.path("Genomes", "hg38+ds")
        proc = self.start(STUB_HANG="add-variants ds chr2.fa")
        self.wait_for_mark(proc)
        self.lock(enriched)
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {enriched} is locked by the run that builds it", err)
        self.assertEqual(sorted(os.listdir(self.path("Genomes"))), ["hg38", "hg38+ds.lock"])
        self.assertFalse(os.path.lexists(self.path("Dictionaries", "dictionaries_ds")))

    def test_a_lock_taken_while_the_dictionaries_are_compressed_moves_nothing_into_place(self):
        """Compressing takes minutes at scale: a search that starts meanwhile takes
        the free lock and enriches the genome itself."""
        enriched = self.path("Genomes", "hg38+ds")
        proc = self.start(STUB_PIGZ_HANG="1")
        self.wait_for_mark(proc)
        self.lock(enriched)
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {enriched} is locked by the run that builds it", err)
        self.assertEqual(sorted(os.listdir(self.path("Genomes"))), ["hg38", "hg38+ds.lock"])
        for d in ("dictionaries_ds", "log_indels_ds"):
            self.assertFalse(os.path.lexists(self.path("Dictionaries", d)), d)
        self.assertEqual([c.split()[1] for c in self.calls()[1:]],
                         [self.path("VCFs", "ds") + "/"])

    def test_a_lock_taken_while_the_fake_indel_genome_is_overlaid_moves_nothing_into_place(self):
        """The SNP overlay, after the compression, takes minutes at scale too (it
        reads each enriched chromosome whole)."""
        enriched = self.path("Genomes", "hg38+ds")
        proc = self.start(shim=True, STUB_OVERLAY_HANG="1")
        self.wait_for_mark(proc)
        self.lock(enriched)
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {enriched} is locked by the run that builds it", err)
        self.assertEqual(sorted(os.listdir(self.path("Genomes"))), ["hg38", "hg38+ds.lock"])
        for d in ("dictionaries_ds", "log_indels_ds"):
            self.assertFalse(os.path.lexists(self.path("Dictionaries", d)), d)
        self.assertEqual([c.split()[1] for c in self.calls()[1:]],
                         [self.path("VCFs", "ds") + "/"])


class TheBuildDirectoriesOfRunsThatAreGoneAreReclaimed(_Fixture):
    """A SIGKILL leaves a build directory, holding a partial index or enrichment as
    large as the finished one. Each build first removes those whose owner record
    names a run that is gone on this host and PID namespace, and only those."""

    def child(self, state=None):
        """A process of this test's: running, stopped ("T"), or a zombie ("Z")."""
        proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
        self.addCleanup(proc.wait)
        self.addCleanup(proc.kill)
        if state:
            os.kill(proc.pid, signal.SIGSTOP if state == "T" else signal.SIGKILL)
            deadline = time.monotonic() + 30
            while proc_state(proc.pid) != state:
                self.assertLess(time.monotonic(), deadline, f"pid {proc.pid} never got {state}")
                time.sleep(0.01)
        return proc.pid

    def test_only_those_of_runs_that_are_gone_are_removed(self):
        if not os.path.exists("/proc/self/stat"):
            self.skipTest("no /proc here")
        dead = subprocess.Popen([sys.executable, "-c", "pass"])
        dead.wait()
        live = self.child()
        lib, genomes = self.path("genome_library"), self.path("Genomes")
        dicts = self.path("Dictionaries")
        dirs = {  # path: (its owner's PID, other fields of its record or None for no record, removed?)
            os.path.join(lib, ".build-index-only-dead"): (dead.pid, {}, True),
            os.path.join(lib, ".build-index-only-zombie"): (self.child("Z"), {}, True),
            os.path.join(lib, ".build-index-only-reused"): (live, {"pstart": "1"}, True),
            os.path.join(genomes, ".build-index-only-enrich-dead"): (dead.pid, {}, True),
            os.path.join(dicts, ".build-index-only-dicts-dead"): (dead.pid, {}, True),
            os.path.join(lib, ".build-index-only-live"): (live, {}, False),
            os.path.join(lib, ".build-index-only-stopped"): (self.child("T"), {}, False),
            os.path.join(genomes, ".build-index-only-enrich-live"): (live, {}, False),
            os.path.join(dicts, ".build-index-only-dicts-live"): (live, {}, False),
            os.path.join(lib, ".build-index-only-elsewhere"): (dead.pid, {"host": "otherhost"}, False),
            os.path.join(lib, ".build-index-only-container"): (dead.pid, {"pidns": "pid:[1]"}, False),
            os.path.join(lib, ".build-index-only-no-pid"): (dead.pid, {"pid": "x"}, False),
            os.path.join(lib, ".build-index-only-no-record"): (dead.pid, None, False),
            # not of this kind: a search's scratch, and other hidden directories
            os.path.join(genomes, "run_a1b2c3"): (dead.pid, {}, False),
            os.path.join(genomes, ".stage_hg19"): (dead.pid, {}, False),
            os.path.join(lib, ".extract_NRG_2_hg19"): (dead.pid, {}, False),
        }
        for path, (pid, fields, _) in dirs.items():
            _write(os.path.join(path, "genome_library", "NRG_2_hg38", "NRG_chr1_1.bin"), "partial\n")
            if fields is not None:
                _write(os.path.join(path, "owner"), owner_record(pid, **fields))
        # one whose record is still being written (owner.<pid>, before its rename), and
        # a symlink by the name of a build directory, to a directory that must stay
        _write(os.path.join(lib, ".build-index-only-torn", f"owner.{dead.pid}"),
               owner_record(dead.pid))
        _write(self.path("elsewhere", "owner"), owner_record(dead.pid))
        os.symlink(self.path("elsewhere"), os.path.join(lib, ".build-index-only-link"))
        before = tree(self.work)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        gone = sorted(path for path, (_, _, removed) in dirs.items() if removed)
        self.assertEqual(
            sorted(line for line in out.splitlines() if line.startswith("Removing ")),
            [f"Removing {path}: the build-index-only run that made it (pid {dirs[path][0]} "
             f"on {os.uname().nodename}, started 2026-10-31T01:02:03Z) is gone." for path in gone])
        # where the owner cannot be checked, a line says how to remove it once it is gone
        self.assertEqual(
            sorted(line for line in out.splitlines() if line.startswith("Leaving ")),
            [f"Leaving {os.path.join(lib, name)}: the build-index-only run that made it (pid "
             f"{dead.pid} on {host}, started 2026-10-31T01:02:03Z) ran in another PID namespace "
             f"(another container, or this one before a restart) or on another host, so this "
             f"run cannot tell whether it still runs. If it does not, remove it: rm -rf "
             f"'{os.path.join(lib, name)}'"
             for name, host in ((".build-index-only-container", os.uname().nodename),
                                (".build-index-only-elsewhere", "otherhost"))])
        for path in gone:
            self.assertFalse(os.path.lexists(path), path)
        # and all else as it was
        after = tree(self.work)
        gone = tuple(os.path.relpath(path, self.work) + "/" for path in gone)
        self.assertEqual({k: v for k, v in before.items() if not (k + "/").startswith(gone)},
                         {k: v for k, v in after.items() if k in before})
        for name, (_, _, whole) in WHOLE.items():
            self.assertEqual(self.bins(name), whole)

    def test_one_whose_output_another_run_has_put_in_place_since_is_removed_too(self):
        """After a kill, the next run is often a complete-search, which builds the
        missing index or enriches the genome itself, in place (here the output is put
        back as it leaves it). The next build-index-only then builds nothing there
        that needs a build directory, and still removes the one the kill left."""
        kinds = (("index", ("genome_library", REF), dict(vcf=False, STUB_HANG=WHOLE[REF][1])),
                 ("enrichment", ("Genomes", "hg38+ds"), dict(STUB_PIGZ_HANG="1")))
        for kind, where, killed in kinds:
            with self.subTest(kind=kind):
                self.fresh()
                vcf = killed.get("vcf", True)
                rc, out, err = self.run_bio(vcf=vcf)
                self.assertEqual(rc, 0, err)
                output, aside = self.path(*where), os.path.join(self.root, "aside")
                os.rename(output, aside)
                proc = self.start(**killed)
                self.wait_for_mark(proc)
                self.kill_group(proc)
                left = self.leftovers()
                self.assertEqual(len(left), 1)
                os.rename(aside, output)
                rc, out, err = self.run_bio(vcf=vcf)
                self.assertEqual(rc, 0, err)
                self.assertIn(f" already present: {output}\n", out)
                self.assertIn(f"Removing {left[0]}: the build-index-only run that made it "
                              f"(pid {proc.pid} on {os.uname().nodename}, started ", out)
                self.assertEqual(self.leftovers(), [])


# Runs overlay_indel_snps.overlay_fake_chromosome on argv[2:5] (a fake-indel FASTA,
# the enriched chromosome, the indel log), with argv[1] the PostProcess directory.
# Given a mark file as argv[5], it stops at the first write of the overlaid FASTA,
# once the mark is made, for a kill -- whichever file object it writes through.
OVERLAY = r"""
import os, sys, time
sys.path.insert(0, sys.argv[1])
import overlay_indel_snps

if sys.argv[5:]:
    def stop_at_the_first_write(frame, event, arg):
        if (event == "c_call" and getattr(arg, "__name__", "") == "write"
                and getattr(arg, "__self__", None) not in (sys.stdout, sys.stderr)):
            open(sys.argv[5], "w").close()
            time.sleep(600)
    sys.setprofile(stop_at_the_first_write)
try:
    overlay_indel_snps.overlay_fake_chromosome(*sys.argv[2:5])
except OSError as exc:
    print("OSError", exc.errno)
    sys.exit(3)
"""


class TheSNPOverlayRewritesAFakeGenomeWholeOrNotAtAll(_Fixture):
    """Before indexing the fake-indel genome, build-index-only overlays SNP codes onto
    each of its FASTA files, rewriting it in place. A run stopped while it wrote one,
    or a write that failed, left that FASTA cut short, and the next run took it for
    a whole one: its INDELS index then had no target past the cut, with no error.
    The fake-indel genome of an enrichment the run makes is overlaid in that
    enrichment's build directory (TheFakeIndelGenomeGoesIntoPlaceWithItsSNPCodes); one
    in place, that came another way, is overlaid in place."""

    def overlays(self):
        """The overlay's hidden directories anywhere in the working directory."""
        return [d for d, _, _ in os.walk(self.work) if os.path.basename(d).startswith(".overlay-")]

    def fixture(self, pad=6000):
        """A fake-indel FASTA holding one indel's window (then N), the enriched
        chromosome with a SNP in that window's left flank, and the indel log; and
        that FASTA's text, before and after the overlay."""
        seq = GENOME["chr1"]
        fake = seq[:26] + "T" + seq[26:60]
        files = (os.path.join(self.root, "fake_ds_chr1", "fakechr1.fa"),
                 os.path.join(self.root, "chr1.enriched.fa"),
                 os.path.join(self.root, "logchr1.txt"))
        _write(files[0], f">fakechr1\n{fake}{'N' * pad}\n")
        _write(files[1], f">chr1\n{seq[:10]}R{seq[11:]}\n")
        _write(files[2], f"{LOG_HEADER}\nchr1_0-60_1\tS2\trs1\t0.25\tchr1_26_{seq[25]}_"
                         f"{seq[25]}T\t0,{len(fake)}\t{seq[:60]}\n")
        whole = f">fakechr1\n{fake}{'N' * pad}\n"
        return files, whole, f">fakechr1\n{fake[:10]}R{fake[11:]}{'N' * pad}\n"

    def overlay(self, files, *mark, fsize=None):
        def limit():
            if fsize:
                resource.setrlimit(resource.RLIMIT_FSIZE,
                                   (fsize, resource.getrlimit(resource.RLIMIT_FSIZE)[1]))

        proc = subprocess.Popen([sys.executable, "-c", OVERLAY, HERE, *files, *mark],
                                stdout=subprocess.PIPE, text=True, preexec_fn=limit,
                                env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        self.addCleanup(proc.wait)
        self.addCleanup(proc.kill)
        return proc

    def listed(self, fasta):
        """What a FASTA listing of its directory sees beside it: index-genome and a
        FASTA search take every file there."""
        d = os.path.dirname(fasta)
        return sorted(f for f in os.listdir(d) if not os.path.isdir(os.path.join(d, f)))

    def test_a_kill_while_it_writes_leaves_the_fasta_whole(self):
        files, whole, _ = self.fixture()
        mark = os.path.join(self.root, "mark")
        proc = self.overlay(files, mark)
        deadline = time.monotonic() + 60
        while not os.path.exists(mark):
            self.assertIsNone(proc.poll(), "the overlay ended before it wrote")
            self.assertLess(time.monotonic(), deadline, "the overlay never wrote")
            time.sleep(0.01)
        proc.kill()
        proc.communicate()
        with open(files[0]) as fh:
            self.assertEqual(fh.read(), whole)
        self.assertEqual(self.listed(files[0]), ["fakechr1.fa"])

    def test_a_write_that_fails_leaves_the_fasta_whole(self):
        files, whole, _ = self.fixture()
        proc = self.overlay(files, fsize=4096)  # less than the FASTA: a write fails partway
        out, _ = proc.communicate(timeout=60)
        self.assertEqual((proc.returncode, out), (3, f"OSError {errno.EFBIG}\n"))
        with open(files[0]) as fh:
            self.assertEqual(fh.read(), whole)
        self.assertEqual(os.listdir(os.path.dirname(files[0])), ["fakechr1.fa"])

    def test_the_fasta_it_writes_is_the_same_with_the_same_mode_and_through_a_symlink(self):
        files, _, overlaid = self.fixture()
        os.chmod(files[0], 0o640)
        # (as root, an owner that is not root, as a shared install can have)
        owner = (4242, 4242) if os.geteuid() == 0 else tuple(os.stat(files[0])[4:6])
        os.chown(files[0], *owner)
        proc = self.overlay(files)
        proc.communicate(timeout=60)
        self.assertEqual(proc.returncode, 0)
        with open(files[0]) as fh:
            self.assertEqual(fh.read(), overlaid)
        self.assertEqual(stat.S_IMODE(os.stat(files[0]).st_mode), 0o640)
        self.assertEqual(tuple(os.stat(files[0])[4:6]), owner)
        self.assertEqual(os.listdir(os.path.dirname(files[0])), ["fakechr1.fa"])
        # a FASTA that is a symlink: the file it names is overlaid, and it stays a link
        files, _, overlaid = self.fixture()
        store = os.path.join(self.root, "store", "fakechr1.fa")
        os.makedirs(os.path.dirname(store))
        os.rename(files[0], store)
        os.symlink(store, files[0])
        proc = self.overlay(files)
        proc.communicate(timeout=60)
        self.assertEqual(proc.returncode, 0)
        self.assertTrue(os.path.islink(files[0]))
        with open(store) as fh:
            self.assertEqual(fh.read(), overlaid)
        self.assertEqual(os.listdir(os.path.dirname(store)), ["fakechr1.fa"])

    def test_a_failed_overlay_fails_the_build_and_indexes_nothing(self):
        """A write inside the overlay that fails (here past a file-size limit, as on
        a full disk): the run fails rather than index what is left. Of the enrichment
        it makes, it puts nothing in place; of one in place, each fake FASTA is left
        whole -- chr1's overlaid, chr2's as add-variants wrote it."""
        with self.subTest(enrichment="this run's"):
            indels = self.path("Genomes", "hg38+ds_INDELS")
            rc, out, err = self.run_bio(fsize=1 << 19, STUB_BIG_FAKE="chr2")
            self.assertEqual(rc, 1, out + err)
            self.assertIn("Error: The SNP overlay of the ds fake-indel genome failed (", err)
            self.assertIn("): nothing of the enrichment was put in place.\n", err)
            for output in (self.path("Genomes", "hg38+ds"), indels,
                           self.path("Dictionaries", "dictionaries_ds"),
                           self.path("Dictionaries", "log_indels_ds")):
                self.assertFalse(os.path.lexists(output), f"{output} was put in place")
            self.assertEqual(self.leftovers(), [])
            # once the cause is gone, the next run does it all
            rc, out, err = self.run_bio()
            self.assertEqual(rc, 0, err)
            self.assertEqual(self.bins(INDELS), WHOLE[INDELS][2])
        with self.subTest(enrichment="in place"):
            self.fresh()
            indels = self.path("Genomes", "hg38+ds_INDELS")
            self.kit(pad=1000000)
            rc, out, err = self.run_bio(fsize=1 << 19)
            self.assertEqual(rc, 1, out + err)
            self.assertIn(f"Error: The SNP overlay of the fake-indel genome {indels} failed (", err)
            self.assertIn(f"), so {INDELS} is not built.", err)
            self.assertFalse(os.path.lexists(self.lib(INDELS)))
            self.assertEqual([c for c in self.calls() if "fakechr" in c], [])
            built = dict(reversed(line.split("  ", 1)) for line in BEFORE["variant"]["tree"])
            chr1 = "Genomes/hg38+ds_INDELS/fake_ds_chr1/fakechr1.fa"
            self.assertEqual(tree(self.work)[chr1], built[chr1])
            seq = GENOME["chr2"]
            with open(os.path.join(indels, "fake_ds_chr2", "fakechr2.fa")) as fh:
                self.assertEqual(fh.read(), f">fakechr2\n{seq[:26]}T{seq[26:60]}{'N' * 1000000}\n")
            for chrom in ("chr1", "chr2"):
                fasta = os.path.join(indels, f"fake_ds_{chrom}", f"fake{chrom}.fa")
                self.assertEqual(os.listdir(os.path.dirname(fasta)), [f"fake{chrom}.fa"])
            # once the cause is gone, the next run finishes it
            rc, out, err = self.run_bio()
            self.assertEqual(rc, 0, err)
            self.assertEqual(self.bins(INDELS), WHOLE[INDELS][2])

    def test_a_sigterm_while_it_overlays_ends_the_run_and_leaves_no_temporary_file(self):
        """scancel, docker stop: the run ends by the signal, as before, with no
        temporary file left: of the enrichment it makes, nothing is in place; of one
        in place, each fake FASTA is as add-variants wrote it."""
        for kit in (False, True):
            with self.subTest(enrichment="in place" if kit else "this run's"):
                self.fresh()
                if kit:
                    self.kit()
                proc = self.start(shim=True, STUB_OVERLAY_HANG="1")
                self.wait_for_mark(proc)
                os.kill(proc.pid, signal.SIGTERM)
                proc.wait(timeout=60)
                self.assertEqual(proc.returncode, -signal.SIGTERM, self.stderr())
                self.assertFalse(os.path.lexists(self.lib(INDELS)))
                self.assertEqual(self.leftovers(), [])
                self.assertEqual(self.overlays(), [])
                if not kit:
                    for output in (self.path("Genomes", "hg38+ds"),
                                   self.path("Genomes", "hg38+ds_INDELS")):
                        self.assertFalse(os.path.lexists(output), f"{output} was put in place")
                    continue
                for chrom, seq in GENOME.items():
                    fasta = self.path("Genomes", "hg38+ds_INDELS", f"fake_ds_{chrom}",
                                      f"fake{chrom}.fa")
                    self.assertEqual(os.listdir(os.path.dirname(fasta)), [f"fake{chrom}.fa"])
                    with open(fasta) as fh:
                        self.assertEqual(fh.read(), f">fake{chrom}\n{seq[:26]}T{seq[26:60]}\n")

    def test_a_kill_while_it_overlays_leaves_its_temporary_file_to_the_next_run(self):
        """An OOM kill or a container stop while the overlay writes a FASTA's
        temporary file leaves that file. Of the enrichment the run makes, it is in
        that enrichment's build directory, which the next run removes before it
        enriches again. Of one in place, it is in a hidden directory beside the FASTA,
        which holds the run's owner record, and the next run removes it, whether it
        overlays again or finds the INDELS index built since (by a complete-search,
        here put back as that leaves it)."""
        with self.subTest(enrichment="this run's"):
            fake = self.path("Genomes", "hg38+ds_INDELS", "fake_ds_chr1")
            proc = self.start(shim=True, STUB_OVERLAY_HANG="1")
            self.wait_for_mark(proc)
            self.kill_group(proc)
            self.assertFalse(os.path.lexists(fake), "the fake-indel genome was put in place")
            left = self.leftovers()
            self.assertEqual([os.path.dirname(d) for d in left], [self.path("Genomes")])
            # its temporary file is there with the rest of the enrichment, in no listing
            self.assertEqual([d.startswith(left[0] + os.sep) for d in self.overlays()], [True])
            rc, out, err = self.run_bio()
            self.assertEqual(rc, 0, err)
            self.assertIn(f"Removing {left[0]}: the build-index-only run that made it "
                          f"(pid {proc.pid} on {os.uname().nodename}, started ", out)
            self.assertIn("Enriching hg38 with ds (add-variants)...", out)
            self.assertEqual(self.leftovers(), [])
            self.assertEqual(self.overlays(), [])
            self.assertEqual(os.listdir(fake), ["fakechr1.fa"])
            self.assertEqual(self.bins(INDELS), WHOLE[INDELS][2])
        for since in (False, True):
            with self.subTest(enrichment="in place", indels_index_built_since=since):
                self.fresh()
                fake = self.path("Genomes", "hg38+ds_INDELS", "fake_ds_chr1")
                if since:
                    rc, out, err = self.run_bio()
                    self.assertEqual(rc, 0, err)
                    aside = os.path.join(self.root, "aside")
                    os.rename(self.lib(INDELS), aside)
                else:
                    self.kit()
                proc = self.start(shim=True, STUB_OVERLAY_HANG="1")
                self.wait_for_mark(proc)
                self.kill_group(proc)
                left = [os.path.join(fake, d) for d in os.listdir(fake) if d.startswith(".overlay-")]
                self.assertEqual(len(left), 1)
                with open(os.path.join(left[0], "owner")) as fh:
                    self.assertIn(f"\npid={proc.pid}\n", fh.read())
                if since:
                    os.rename(aside, self.lib(INDELS))
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                self.assertIn(f"Removing {left[0]}: the build-index-only run that made it "
                              f"(pid {proc.pid} on {os.uname().nodename}, started ", out)
                self.assertEqual(os.listdir(fake), ["fakechr1.fa"])
                self.assertEqual(self.bins(INDELS), WHOLE[INDELS][2])


class TheFakeIndelGenomeGoesIntoPlaceWithItsSNPCodes(_Fixture):
    """build-index-only added the SNP codes to the fake-indel genome only just before
    it indexed it, so a run stopped after the enrichment (in the tier emission, the
    indel genotype store, the variant index) left it in place without them. A
    complete-search then built the INDELS index from it, with pool_index_indels.py,
    which adds none, and every later build-index-only kept that one: the off-targets
    that need both an indel and a nearby SNP were missing or scored too low, with no
    error. The codes now go on in the enrichment's build directory, before any of it
    is in place."""

    def search_builds_the_indels_index(self):
        """What complete-search does when it finds the fake-indel genome but no
        INDELS index: pool_index_indels.py, in the working directory, with the
        variant index's PAM and N."""
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("STUB_", "CRISPRME_"))}
        env.update(PATH=os.pathsep.join([self.bin, os.environ.get("PATH", "")]),
                   STUB_LOG=self.log)
        subprocess.run([sys.executable, os.path.join(HERE, "pool_index_indels.py"),
                        self.path("Genomes", "hg38+ds_INDELS") + "/", PAM, "NRG", "hg38", "ds",
                        "2", "2"], cwd=self.work, env=env, check=True, capture_output=True)

    def test_an_indels_index_a_search_builds_from_it_is_the_one_build_index_only_builds(self):
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = tree(self.lib(INDELS)), tree(self.path("Genomes", "hg38+ds_INDELS"))
        # what the overlay prints comes where it did, as it overlaid just before indexing
        lines = out.splitlines()
        at = lines.index(f"Building indels index {INDELS}...")
        self.assertEqual(lines[at + 1:at + 5], [
            "overlay_indel_snps: chr1 -> 1 SNP overlays",
            "overlay_indel_snps: chr2 -> 1 SNP overlays",
            "overlay_indel_snps: 2 total SNP overlays",
            "[indel-snp] SNP-overlaid the fake-indel genome (2 flank bases IUPAC-coded).",
        ])
        cases = (("killed while it builds the variant index", WHOLE[VAR][1], True),
                 ("killed while it builds the INDELS index", WHOLE[INDELS][1], True),
                 ("a search builds the INDELS index meanwhile", WHOLE[VAR][1], False))
        for case, step, killed in cases:
            with self.subTest(case):
                self.fresh()
                proc = self.start(STUB_HANG=step)
                self.wait_for_mark(proc)
                if killed:
                    self.kill_group(proc)
                self.search_builds_the_indels_index()
                if killed:
                    rc, out, err = self.run_bio()
                else:
                    open(self.release, "w").close()
                    rc, out, err = self.finish(proc)
                self.assertEqual(rc, 0, err)
                self.assertIn(f"Indels index already present: {self.lib(INDELS)}\n", out)
                self.assertEqual(tree(self.lib(INDELS)), clean[0])
                self.assertEqual(tree(self.path("Genomes", "hg38+ds_INDELS")), clean[1])
                self.assertEqual(self.leftovers(), [])

    def test_one_in_place_without_them_gets_them_as_the_run_begins(self):
        """One that came another way (a kit, an earlier release, a search's own) got
        them at STEP 3 only, so a run stopped before then (the tiers, the indel
        genotype store, the variant index: hours at scale) left it without them, and
        a search built the INDELS index from it, SNP-blind. It now gets them right
        after the run finds it in place; STEP 3 prints what that printed, as before."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = tree(self.lib(INDELS)), tree(self.path("Genomes", "hg38+ds_INDELS"))
        self.fresh()
        self.kit()
        shutil.rmtree(self.lib(VAR))
        proc = self.start(STUB_HANG=WHOLE[VAR][1])
        self.wait_for_mark(proc)
        self.kill_group(proc)
        self.assertEqual(tree(self.path("Genomes", "hg38+ds_INDELS")), clean[1])
        self.search_builds_the_indels_index()
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"Indels index already present: {self.lib(INDELS)}\n", out)
        self.assertEqual(tree(self.lib(INDELS)), clean[0])
        # and a run to its end prints the overlay where it did
        self.fresh()
        self.kit()
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        lines = out.splitlines()
        at = lines.index(f"Building indels index {INDELS}...")
        self.assertEqual(lines[at + 1:at + 5], [
            "overlay_indel_snps: chr1 -> 1 SNP overlays",
            "overlay_indel_snps: chr2 -> 1 SNP overlays",
            "overlay_indel_snps: 2 total SNP overlays",
            "[indel-snp] SNP-overlaid the fake-indel genome (2 flank bases IUPAC-coded).",
        ])
        self.assertEqual(tree(self.lib(INDELS)), clean[0])

    def test_a_rerun_says_the_fake_indel_genome_already_carries_its_codes(self):
        """A re-run, after a kill while the variant index was built, overlays this
        release's own enrichment, coded before it went into place: no base changes,
        and the lines say why, rather than read as an empty or failed overlay."""
        proc = self.start(STUB_HANG=WHOLE[VAR][1])
        self.wait_for_mark(proc)
        self.kill_group(proc)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        lines = out.splitlines()
        at = lines.index(f"Building indels index {INDELS}...")
        self.assertEqual(lines[at + 1:at + 5], [
            "overlay_indel_snps: chr1 -> 0 SNP overlays (1 already coded)",
            "overlay_indel_snps: chr2 -> 0 SNP overlays (1 already coded)",
            "overlay_indel_snps: 0 total SNP overlays (2 already coded)",
            "[indel-snp] SNP-overlaid the fake-indel genome (0 flank bases IUPAC-coded; it "
            "already carries its SNP codes, 2 of them).",
        ])
        self.assertEqual(self.bins(INDELS), WHOLE[INDELS][2])

    def test_a_search_refuses_it_while_a_run_stopped_partway_has_left_it_half_coded(self):
        """A run stopped while it overlays one in place leaves some FASTAs overlaid
        and some not (here chr1's and chr2's): an INDELS index a search builds from
        them is blind to the SNPs on the others. The next run overlays the rest."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = tree(self.lib(INDELS))
        self.fresh()
        self.kit()
        indels = self.path("Genomes", "hg38+ds_INDELS")
        proc = self.start(shim=True, STUB_WRITE_HANG="fakechr2.fa.tmp")
        self.wait_for_mark(proc)
        self.kill_group(proc)
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {indels} is unfinished: the build-index-only run that wrote the "
                      f"SNP codes of the ds fake-indel genome there (pid {proc.pid} on "
                      f"{os.uname().nodename}, started ", err)
        self.assertIn(") stopped before it was done, and a search would use what is there as "
                      "if it were whole. Run build-index-only again with ds, which finishes "
                      "them, then run this search again.\n", err)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertEqual(sorted(os.listdir(indels)), ["fake_ds_chr1", "fake_ds_chr2"])
        self.assertEqual(tree(self.lib(INDELS)), clean)
        self.assert_searched()

    def test_an_indels_index_built_while_it_is_overlaid_in_place_is_built_again(self):
        """A search that got past its check of the mark before the run set it, or one
        that does not check it, builds the INDELS index from FASTAs not all overlaid
        yet, blind to the SNPs on the others, and every later run kept it ('Indels
        index already present'). The run removes it once its overlay ends, and builds
        it again; while a search still builds it, the run stops, its mark left for the
        next run to do so."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = tree(self.lib(INDELS))
        for case in ("built as it overlays", "still being built as it ends",
                     "built beside a stopped run's mark"):
            with self.subTest(case=case):
                self.fresh()
                self.kit()
                removing = (f"Removing {self.lib(INDELS)}: it was built from the fake-indel "
                            "genome while its SNP codes were going on, blind to those not on "
                            "yet. This run builds it again.\n")
                if case == "built beside a stopped run's mark":
                    proc = self.start(shim=True, STUB_WRITE_HANG="fakechr2.fa.tmp")
                    self.wait_for_mark(proc)
                    self.kill_group(proc)
                    self.search_builds_the_indels_index()  # chr1's FASTA coded, chr2's not
                    rc, out, err = self.run_bio()
                else:
                    proc = self.start(shim=True, STUB_OVERLAY_HANG="1")
                    self.wait_for_mark(proc)  # (before the first FASTA is overlaid)
                    if case == "built as it overlays":
                        self.search_builds_the_indels_index()
                    else:
                        self.lock(self.lib(INDELS))
                    open(self.release, "w").close()
                    rc, out, err = self.finish(proc)
                    if case == "still being built as it ends":
                        self.assertEqual(rc, 1, out + err)
                        self.assertIn(f"Error: {self.lib(INDELS)} is locked by the run that "
                                      "builds it", err)
                        self.assertTrue(os.path.isdir(self.path("Genomes", "hg38+ds_INDELS",
                                                                ".building")))
                        # that search ends, its INDELS index in place
                        shutil.rmtree(self.lib(INDELS) + ".lock")
                        self.search_builds_the_indels_index()
                        rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                self.assertIn(removing, out)
                self.assertIn(f"Building indels index {INDELS}...\n", out)
                self.assertEqual(tree(self.lib(INDELS)), clean)
                self.assertFalse(os.path.lexists(self.path("Genomes", "hg38+ds_INDELS",
                                                           ".building")))
                self.assertEqual(self.leftovers(), [])


class AStoppedRunsTiersAreNotSearched(_Fixture):
    """Once the enrichment is in place, build-index-only writes the dict-less tiers
    (registry_<vcf>/, genotypes_<vcf>/) and the phased indel genotype store
    (indel_genotypes_<vcf>/) in place, chromosome by chromosome, and a search takes
    them chromosome by chromosome. A run stopped there (an OOM kill: the tiers take
    ~150-180 GB a chromosome at scale) left the part written so far, and every
    complete-search used it as if it were whole, with no error: no SNP+indel
    co-occurrence on the chromosomes not written yet, other best alignments. Each is
    now marked unfinished while it is written, and complete-search refuses it until a
    build-index-only run has finished it."""

    # each step, the write of its second chromosome (held or failed there), and the
    # directory it writes (one worker, so the first chromosome is written by then)
    STEPS = (("tiers", "/registry_ds/reg_chr2.bin", "registry_ds",
              dict(CRISPRME_TIER_WORKERS="1")),
             ("indel genotype store", "/indel_genotypes_ds/gt_indel_chr2",
              "indel_genotypes_ds", dict(CRISPRME_INDELGT_WORKERS="1")))

    def test_a_run_stopped_while_it_writes_them_leaves_them_refused_until_one_finishes_them(self):
        self.run_bio()
        clean = self.outputs()
        for step, held, where, workers in self.STEPS:
            with self.subTest(step=step):
                self.fresh()
                proc = self.start(shim=True, STUB_WRITE_HANG=held, **workers)
                self.wait_for_mark(proc)
                self.kill_group(proc)
                directory = self.path("Dictionaries", where)
                rc, out, err = self.run_search()
                self.assertEqual(rc, 1, out + err)
                self.assertIn(f"Error: {directory} is unfinished: the build-index-only run that "
                              f"wrote the ", err)
                self.assertIn(f" there (pid {proc.pid} on {os.uname().nodename}, started ", err)
                self.assertIn("Run build-index-only again with ds", err)
                self.assertNotIn("Launching job ", out)
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                if step == "tiers":
                    self.assertIn("  1 chromosome(s) already emitted; resuming 1 remaining\n", out)
                self.assertEqual(self.outputs(), clean)
                self.assert_searched()

    def test_a_run_stopped_as_its_enrichment_goes_into_place_leaves_them_refused(self):
        """A SIGTERM as the enrichment goes into place (SHIM sends one there) waits for
        the last rename and then ends the run, before the tiers and the indel genotype
        store are begun: a whole enrichment without them, which a search took for one
        built without --samplesID. They are marked unfinished before the enrichment
        goes into place, until the steps that write them end."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = self.outputs()
        self.fresh()
        rc, out, err = self.run_bio(
            shim=True, STUB_SIGNAL_AT_RENAME=self.path("Dictionaries", "dictionaries_ds"))
        self.assertEqual(rc, -signal.SIGTERM, err)
        for where in ("registry_ds", "indel_genotypes_ds"):
            self.assertTrue(os.path.isdir(self.path("Dictionaries", where, ".building")), where)
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {self.path('Dictionaries', 'registry_ds')} is unfinished: the "
                      "build-index-only run that wrote the dict-less tiers of ds there", err)
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"Enriched genome already present: {self.path('Genomes', 'hg38+ds')}\n", out)
        self.assertEqual(self.outputs(), clean)
        self.assert_searched()

    def test_a_search_launched_from_the_web_refuses_what_a_stopped_run_left(self):
        """The web UI runs the search script straight (pages/main_page.py), past
        complete-search's checks, and that search used the part a stopped run left as
        if it were whole, with exit status 0. The script refuses it too: the job ends
        with the ERROR in its log_error.txt, which the results page shows."""
        cases = [(step, where, dict(shim=True, STUB_WRITE_HANG=held, **workers))
                 for step, held, where, workers in self.STEPS]
        cases.append(("SNP codes of the fake-indel genome", "../Genomes/hg38+ds_INDELS",
                      dict(shim=True, STUB_WRITE_HANG="fakechr2.fa.tmp")))
        for step, where, stub in cases:
            with self.subTest(step=step):
                self.fresh()
                if step.startswith("SNP codes"):
                    self.kit()
                proc = self.start(**stub)
                self.wait_for_mark(proc)
                self.kill_group(proc)
                calls = len(self.calls())
                rc, out, err = self.run_web_search()
                self.assertEqual(rc, 1, out + err)
                self.assertEqual(err, f"ERROR: {os.path.normpath(self.path('Dictionaries', where))} "
                                      "is unfinished: a build-index-only run is writing it, or "
                                      "stopped before it was done. Let that run finish, or run "
                                      "build-index-only again to its end, then run this search "
                                      "again\n")
                self.assertEqual(self.calls()[calls:], [])  # nothing searched or built
                rc, out, err = self.run_bio()
                self.assertEqual(rc, 0, err)
                rc, out, err = self.run_web_search()
                self.assertIn("Using precomputed variant index for hg38+ds; skipping genome "
                              "enrichment", out)
                self.assertNotIn(" is unfinished: ", err)

    def test_a_search_while_a_run_writes_them_is_told_to_wait_for_it(self):
        step, held, where, workers = self.STEPS[0]
        proc = self.start(shim=True, STUB_WRITE_HANG=held, **workers)
        self.wait_for_mark(proc)
        rc, out, err = self.run_search()
        self.assertEqual(rc, 1, out + err)
        self.assertIn(f"Error: {self.path('Dictionaries', where)} is being written: a "
                      f"build-index-only run (pid {proc.pid} on {os.uname().nodename}, started ",
                      err)
        self.assertIn(") is writing the dict-less tiers of ds there, and a search would use "
                      "the part written so far as if it were whole. Wait for that run to "
                      "finish, then run this search again.\n", err)
        open(self.release, "w").close()
        rc, out, err = self.finish(proc)
        self.assertEqual(rc, 0, err)
        self.assert_searched()

    def test_a_step_that_ends_with_a_guarded_failure_leaves_them_refused(self):
        """A chromosome whose tiers or indel genotypes failed to be written was skipped
        with a WARNING on stdout, and the step took its mark away as it ended: a search
        then used the part written as if it were whole, and after a full disk (the
        run then failed at the variant index) once space was freed. The mark now stays,
        and the line after the WARNING says so. A full disk or quota, an I/O error or no
        memory left, which every chromosome after would hit too, ends the run with an
        ERROR, the mark kept. Either way, the next run finishes them."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        clean = self.outputs()
        cases = [(why, True, why in ("EACCES", "ENOSPC"))  # (errno, one worker, then the rest)
                 for why in ("EACCES", "ENOSPC", "EDQUOT", "EIO", "MemoryError")]
        cases += [(why, False, False) for why in ("EACCES", "ENOSPC")]  # (a pool of two)
        for step, failed, where, workers in self.STEPS:
            for why, alone, then in cases:
                with self.subTest(step=step, why=why, alone=alone):
                    self.fresh()
                    directory = self.path("Dictionaries", where)
                    rc, out, err = self.run_bio(shim=True, STUB_WRITE_FAIL=failed,
                                                STUB_WRITE_FAIL_ERRNO=why, **(workers if alone else {}))
                    self.assertTrue(os.path.isdir(os.path.join(directory, ".building")), out + err)
                    if why == "EACCES":
                        self.assertEqual(rc, 0, err)
                        lines = out.splitlines()
                        warned = [i for i, line in enumerate(lines) if line.startswith(
                            ("WARNING [build_dictless_tiers]: could not emit dictless tiers for ",
                             "WARNING [indel-snp]: compile_indel_genotypes failed for "))]
                        self.assertEqual(len(warned), 1, out)
                        if step == "tiers":
                            kept = (f"WARNING [build-index-only]: the dict-less tiers of 1 of 2 "
                                    f"chromosome(s) of ds were not emitted (above), so {directory} "
                                    "stays marked unfinished and a search refuses it. Run "
                                    "build-index-only again with ds and its --samplesID, which "
                                    "emits them, or remove the tiers to search without them (rm -rf "
                                    f"'{directory}' '{self.path('Dictionaries', 'genotypes_ds')}').")
                        else:
                            kept = (f"WARNING [indel-snp]: the phased indel genotype store of 1 of 2 "
                                    f"chromosome(s) of ds was not built (above), so {directory} "
                                    "stays marked unfinished and a search refuses it. Run "
                                    "build-index-only again with ds, which builds it, or remove it to "
                                    f"search without it (rm -rf '{directory}').")
                        self.assertIn(kept, lines[warned[0] + 1:], out)
                    else:
                        self.assertEqual(rc, 1, out + err)
                        what = ("Emitting the dict-less tiers" if step == "tiers"
                                else "Building the phased indel genotype store")
                        cause = ("MemoryError" if why == "MemoryError" else
                                 f"[Errno {getattr(errno, why)}] {os.strerror(getattr(errno, why))}")
                        self.assertTrue(err.startswith(f"Error: {what} of ds failed ({cause}"), err)
                        self.assertTrue(err.endswith(
                            f"): {directory} stays marked unfinished, and a search refuses it. "
                            "Run build-index-only again once that is fixed, which "
                            f"{'emits' if step == 'tiers' else 'builds'} the rest.\n"), err)
                        self.assertEqual(err.count("\n"), 1, err)
                    if then:
                        rc, out, err = self.run_search()
                        self.assertEqual(rc, 1, out + err)
                        self.assertIn(f"Error: {directory} is unfinished: ", err)
                        rc, out, err = self.run_bio()
                        self.assertEqual(rc, 0, err)
                        self.assertEqual(self.outputs(), clean)
                        self.assert_searched()
        # and a chromosome's indel genotypes whole, but not the .done marker after them:
        # left empty by a full disk, it was taken for done, and no run wrote it again
        self.fresh()
        marker = self.path("Dictionaries", "indel_genotypes_ds", "gt_indel_chr2.tsv.gz.done")
        rc, out, err = self.run_bio(shim=True, STUB_WRITE_FAIL=os.path.basename(marker),
                                    CRISPRME_INDELGT_WORKERS="1")
        self.assertEqual(rc, 1, out + err)
        self.assertTrue(err.startswith("Error: Building the phased indel genotype store of ds "
                                       f"failed ([Errno 28] No space left on device: '{marker}')"),
                        err)
        self.assertFalse(os.path.lexists(marker[: -len(".done")]))
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.outputs(), clean)

    def test_a_mark_whose_run_cannot_be_checked_is_refused_too(self):
        """One a run left as it wrote its owner record, and one made in another PID
        namespace (another container, or this one before a restart)."""
        if not os.path.exists("/proc/self/stat"):
            self.skipTest("no /proc here")
        marker = self.path("Dictionaries", "registry_ds", ".building")
        for case, record, says in (
                ("no record", None, "is unfinished: the build-index-only run that wrote the "
                                    "dict-less tiers of ds there stopped before it was done"),
                ("another PID namespace", owner_record(os.getpid(), pidns="pid:[1]"),
                 "is being written, or was left unfinished: the build-index-only run that "
                 f"writes the dict-less tiers of ds there (pid {os.getpid()} on "
                 f"{os.uname().nodename}, started 2026-10-31T01:02:03Z) ran in another PID "
                 "namespace")):
            with self.subTest(case=case):
                os.makedirs(marker, exist_ok=True)
                if record:
                    _write(os.path.join(marker, "owner"), record)
                rc, out, err = self.run_search()
                self.assertEqual(rc, 1, out + err)
                self.assertIn(f"Error: {os.path.dirname(marker)} {says}", err)
                self.assertIn(" or remove them to search without them (rm -rf "
                              f"'{os.path.dirname(marker)}' "
                              f"'{self.path('Dictionaries', 'genotypes_ds')}'), then run this "
                              "search again.\n", err)

    def test_publish_index_refuses_them_too(self):
        """Published, they would be searched as if whole wherever they are downloaded."""
        rc, out, err = self.run_bio()
        self.assertEqual(rc, 0, err)
        for where in ("registry_ds", "indel_genotypes_ds"):
            with self.subTest(where=where):
                marker = self.path("Dictionaries", where, ".building")
                _write(os.path.join(marker, "owner"), "host=otherhost\npid=4242\n")
                rc, out, err = self.finish(self.start(
                    command=["publish-index", "--index", self.lib(VAR)]))
                self.assertEqual(rc, 1, out + err)
                self.assertEqual(err, f"Error: {os.path.dirname(marker)} is unfinished, or still "
                                      "being written: the build-index-only run that writes it "
                                      "has not finished it. Run build-index-only again to its "
                                      "end, then publish the index.\n")
                shutil.rmtree(marker)


def save_tst(targets, guide_len=23):
    """A .bin as buildTST writes one (insall, insert, saveTST and serialize in
    CRISPRitz's mainParallel.cpp, ported): targets are (guideIndex, sequence, PAM),
    sorted by sequence, all sequences of one length."""
    code = dict(zip("ACGTRYSWKMBDHVN0_", (1, 2, 4, 8, 5, 10, 6, 9, 12, 3, 14, 13, 11, 7, 0, 0, 15)))
    nodes, nexts = [], [0] * len(targets)  # a node: [base, lokid, hikid, eqkid]

    def insert(seq, i):
        k, p = 0, 0
        while nodes:
            node = nodes[p]
            if seq[k] == node[0]:
                k += 1
                if k == len(seq):  # another target with this sequence
                    nexts[i], node[3] = node[3], -(i + 1)
                    return
                if not node[3]:
                    node[3] = len(nodes)
                    break
                p = node[3]
            else:
                kid = 1 if seq[k] < node[0] else 2
                if not node[kid]:
                    node[kid] = len(nodes)
                    break
                p = node[kid]
        while True:
            nodes.append([seq[k], 0, 0, 0])
            k += 1
            if k == len(seq):
                nodes[-1][3] = -(i + 1)
                return
            nodes[-1][3] = len(nodes)

    def insall(lo, hi):
        if hi >= lo:
            mid = lo + (hi - lo) // 2
            insert(targets[mid][1], mid)
            insall(lo, mid - 1)
            insall(mid + 1, hi)

    insall(0, len(targets) - 1)
    out = bytearray(struct.pack("=ii", len(targets), guide_len))
    for i, (index, _, pam) in enumerate(targets):
        out += struct.pack("=i", index)
        nibble, k = 0, 0
        for j, base in enumerate(pam):
            nibble, k = nibble + code[base], k + 1
            if j == len(pam) - 1 or k == 2:
                out.append(nibble)
                nibble, k = 0, 0
            if k == 1:
                nibble <<= 4
        out += b"_" + struct.pack("=i", nexts[i]) if nexts[i] else b"0"
    out += struct.pack("=i", len(nodes))
    pair, first = ["0", "0"], [True]

    def put(base):
        pair[0 if first[0] else 1] = base
        if base == "_" or not first[0]:
            out.append(0xF0 if pair[0] == "_" else code[pair[0]] << 4 | code[pair[1]])
            first[0] = True
        else:
            first[0] = False

    def serialize(p):
        base, lo, hi, eq = nodes[p]
        put(base)
        serialize(lo) if lo > 0 else put("0")
        serialize(hi) if hi > 0 else put("0")
        if eq > 0:
            serialize(eq)
        else:
            put("_")
            out.extend(struct.pack("=i", eq))

    serialize(0)
    return bytes(out)


class TheIndexCheckTellsAWholeTreeFromOneCutShort(unittest.TestCase):
    """crispritz_bins, on trees as buildTST saves them: of many targets, some of one
    sequence, with shared starts and IUPAC bases, for PAMs of odd and even length.
    Each is whole, and each cut of it, at every length, is not."""

    def test_every_cut_of_a_whole_tree_is_told_apart(self):
        sys.path.insert(0, HERE)
        import crispritz_bins
        rng = random.Random(7)
        tmp = tempfile.mkdtemp(prefix="bio_atomic_bins_")
        self.addCleanup(shutil.rmtree, tmp, True)
        for pam, n in (("NRG", 40), ("TTTV", 25), ("NNNNNNNN", 1), ("NGG", 2)):
            with self.subTest(pam=pam, targets=n):
                seqs = ["".join(rng.choice("ACGTACGTRY") for _ in range(6)) for _ in range(n)]
                seqs += seqs[: n // 4]  # some sequences more than once
                targets = sorted((rng.randrange(1 << 20), s, "".join(rng.choice("ACGTN")
                                                                     for _ in pam))
                                 for s in seqs)
                targets = sorted(targets, key=lambda t: t[1])
                data = save_tst(targets)
                path = os.path.join(tmp, f"{pam}_chr1_1.bin")
                with open(path, "wb") as fh:
                    fh.write(data)
                self.assertIsNone(crispritz_bins.bin_problem(path))
                for size in range(len(data)):
                    with open(path, "wb") as fh:
                        fh.write(data[:size])
                    self.assertIsNotNone(crispritz_bins.bin_problem(path), size)


def buildtst_targets(seq, pam, pam_len, at_start, max_bulges):
    """The number of targets buildTST keeps in contig seq: searchPAMonGenome,
    guideRegionHasN and main in CRISPRitz 2.9.1 (pamGenerator.cpp, mainParallel.cpp),
    loop for loop, as slow as they are. It gives the real buildTST's count."""
    import crispritz_bins

    size, width, n = len(pam), pam_len + max_bulges, len(seq)
    if n < width:
        return 0
    genome = {b: crispritz_bins._BITS.get(b, 0) for b in "ACGTNRYSWKMBDHV"}
    plus = [crispritz_bins._PAM_BITS[b] for b in pam]
    minus = [crispritz_bins._PAM_BITS[b] for b in pam.translate(crispritz_bins._REVERSE)[::-1]]
    free = all(b == 15 for b in plus)
    indices = []
    for nt in range(n - (pam_len if at_start else size)):
        on_plus = free or all(genome[seq[nt + j]] & plus[j] for j in range(size))
        on_minus = free or all(genome[seq[nt + j]] & minus[j] for j in range(size))
        upstream, downstream = nt + size - width >= 0, nt <= n - width
        if not at_start:
            indices += [nt + size - width] * (on_plus and upstream) + [-nt] * (on_minus and downstream)
        else:
            indices += [-nt] * (on_plus and downstream) + [nt + size - width] * (on_minus and upstream)

    def has_n(target, pam_in_front):
        if pam_in_front:
            return target.find("N", size) >= 0
        return 0 <= target.find("N") < len(target) - size

    kept = 0
    for i in indices:
        if (i > 0) if not at_start else (i >= 0):
            kept += not has_n(seq[i:i + width], False)
        else:
            kept += not has_n(seq[-i:-i + width], True)
    return kept


class TheContigCheckFindsATargetWhereBuildTSTDoes(unittest.TestCase):
    """crispritz_bins.holds_target, which tells a contig a full inode table left
    without a .bin from one with no target, against buildTST's own loops
    (buildtst_targets): random contigs of every kind (N runs, IUPAC codes, lower case,
    targets at either end and windows at 0) and the shipped PAMs, 3' and 5', PAM-less
    among them, at 0 to 3 bulges. And missing_bins, on the .bin names buildTST
    writes: <PAM>_<the whole header>_<part>.bin."""

    PAMS = ("NNNNNNNNNNNNNNNNNNNNNRG 3", "NNNNNNNNNNNNNNNNNNNNNGG 3", "NNNNNNNNNNNNNNNNNNNNNNN 3",
            "TTTVNNNNNNNNNNNNNNNNNNNNNNN -4", "NNNNNNNNNNNNNNNNNNNNNNNGRRT 6",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNN -4", "NNNNNNNNNNNNNNNNNNNNNNNNNNNN 8")

    def test_it_agrees_with_buildtst_on_random_contigs(self):
        import crispritz_bins

        r = random.Random(7)
        held = 0
        for case in range(2500):
            motif, limit = r.choice(self.PAMS).split()
            limit = int(limit)
            pam = motif[:-limit] if limit < 0 else motif[-limit:]
            bulges = r.randint(0, 3)
            n = r.choice((r.randint(0, 40), r.randint(20, 60), r.randint(40, 300)))
            mix = r.choice(("ACGT", "ACGT" * 6 + "N", "ACGTN", "ACGTNNNNRYSWKMBDHV", "ACGTNRYSWKMBDHV"))
            seq = "".join(r.choice(mix) for _ in range(n))
            if n > 10 and r.random() < 0.3:
                a = r.randrange(n)
                seq = seq[:a] + "N" * min(n - a, r.randint(1, 30)) + seq[a + 30:]
            want = buildtst_targets(seq, pam, len(motif), limit < 0, bulges) > 0
            held += want
            self.assertEqual(crispritz_bins.holds_target(seq, pam, len(motif), limit < 0, bulges),
                             want, (seq, motif, limit, bulges))
        self.assertGreater(held, 500)
        self.assertLess(held, 2000)

    def test_missing_bins_names_only_a_contig_that_holds_a_target(self):
        import crispritz_bins

        top = tempfile.mkdtemp(prefix="bio_contigs_")
        self.addCleanup(shutil.rmtree, top, True)
        genome, index = os.path.join(top, "genome"), os.path.join(top, "index")
        os.makedirs(index)
        site = "ACGTTGCAAGACGTTGCAAGACGTAGG"  # a target for NRG at 0-2 bulges
        contigs = {
            "chr1.fa": ">chr1\n" + site + "\n",  # its .bin there: NRG_chr1_1.bin
            "chr1_random.fa": ">chr1_random\n" + site + "\n",  # no .bin of its own
            "chr2.fa": ">chr2 desc here\n" + site.lower() + "\n",  # no .bin
            "chr3.fa": ">chr3\n" + "N" * 60 + "\n",  # no target
            "chr4.fa": ">chr4\n" + site[:20] + "\n",  # shorter than a target
            "chr5.fa": ">chr5\n" + site,  # its only line has no end: buildTST drops it
            "chr6.fa": ">chr6\r\n" + site + "\r\n",  # read otherwise: unread, not missing
            "chr7.fa": ">chr7\n" + site + "U\n",  # as is a byte no IUPAC base
            "chr8.fa": ">chr8\n" + site + "\n",  # its .bin of part 2 alone: still there
            "notes.txt": ">chr9\n" + site + "\n",  # not indexed
        }
        for name, text in contigs.items():
            _write(os.path.join(genome, name), text)
        for name in ("NRG_chr1_1.bin", "NRG_chr8_2.bin", "NGG_chr2 desc here_1.bin"):
            open(os.path.join(index, name), "w").close()
        unread = ["chr6.fa", "chr7.fa"]
        self.assertEqual(crispritz_bins.missing_bins(index, genome, PAM, 2),
                         ([("chr1_random.fa", "chr1_random"), ("chr2.fa", "chr2 desc here")], unread))
        # with no index directory at all, every contig that holds a target
        self.assertEqual(crispritz_bins.missing_bins(os.path.join(top, "none"), genome, PAM, 2),
                         ([("chr1.fa", "chr1"), ("chr1_random.fa", "chr1_random"),
                           ("chr2.fa", "chr2 desc here"), ("chr8.fa", "chr8")], unread))
        # and with a PAM file it cannot read, none is told: each without a .bin is unread
        pam = os.path.join(top, "pam.txt")
        _write(pam, "NNNNNNNNNNNNNNNNNNNNNRG\n")
        self.assertEqual(crispritz_bins.missing_bins(index, genome, pam, 2)[0], [])


# What the code before this fix (379b5b0) did with the same stub and inputs, for
# SameTreeAsBefore: its exit status, the tree it left as `digest  path` (tree(), with
# the run's own paths written <W>, <R> and <REPO>), the lines it printed, sorted, and
# the engine calls it made. A change to the stub or the inputs changes these: record
# them again with that code, by running this file with --record there.
BEFORE = {
    'reference': {
        "rc": 0,
        "tree": [
            "dir  Genomes/",
            "dir  Genomes/hg38/",
            "9fb230dc8fb2ed4d  Genomes/hg38/chr1.fa",
            "c5336dd0168ef6b3  Genomes/hg38/chr2.fa",
            "dir  VCFs/",
            "dir  VCFs/ds/",
            "25cf25d958babc38  VCFs/ds/ds.chr1.vcf.gz",
            "26a5b9b8e3b7cfc5  VCFs/ds/ds.chr2.vcf.gz",
            "dir  genome_library/",
            "dir  genome_library/NRG_2_hg38/",
            "4983ff2823805684  genome_library/NRG_2_hg38/.display_label",
            "425fad713f395627  genome_library/NRG_2_hg38/.pam_build",
            "dee1e9b039f12fe0  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "338fbfb5efc7044c  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
            "c442c9cdd4d80f67  samples.txt",
            "dir  samplesIDs/",
            "036a130b3364d826  samplesIDs/hg38_ds.samplesID.txt",
        ],
        "stdout": [
            "Building reference index NRG_2_hg38 (bMax 2, 2 thread(s))...",
            "Index built: <W>/genome_library/NRG_2_hg38",
            "It will be reused automatically by complete-search runs launched from this working directory (or pointed here with --index-path) that use the same genome, PAM and bulge settings.",
        ],
        "stderr": [],
        "calls": [
            "index-genome hg38 <W>/Genomes/hg38/ <REPO>/PAMs/20bp-NRG-SpCas9.txt -bMax 2 -th 2",
        ],
    },
    'variant': {
        "rc": 0,
        "tree": [
            "dir  Dictionaries/",
            "dir  Dictionaries/dictionaries_ds/",
            "d1b69ccc5c7ccb40  Dictionaries/dictionaries_ds/my_dict_chr1.json.gz",
            "497425130c0deb8d  Dictionaries/dictionaries_ds/my_dict_chr2.json.gz",
            "dir  Dictionaries/genotypes_ds/",
            "c62654c4153c1576  Dictionaries/genotypes_ds/gt_chr1.bin",
            "acc7f66ace5d0084  Dictionaries/genotypes_ds/gt_chr1.idx",
            "c62654c4153c1576  Dictionaries/genotypes_ds/gt_chr2.bin",
            "acc7f66ace5d0084  Dictionaries/genotypes_ds/gt_chr2.idx",
            "dir  Dictionaries/indel_genotypes_ds/",
            "a3f726ac65dbc9e4  Dictionaries/indel_genotypes_ds/gt_indel_chr1.tsv.gz",
            "4355a46b19d348dc  Dictionaries/indel_genotypes_ds/gt_indel_chr1.tsv.gz.done",
            "8e3f9e931ce698b8  Dictionaries/indel_genotypes_ds/gt_indel_chr2.tsv.gz",
            "4355a46b19d348dc  Dictionaries/indel_genotypes_ds/gt_indel_chr2.tsv.gz.done",
            "dir  Dictionaries/log_indels_ds/",
            "5717b5bb7539d559  Dictionaries/log_indels_ds/logchr1.txt.gz",
            "0004c31cfbfba66d  Dictionaries/log_indels_ds/logchr2.txt.gz",
            "dir  Dictionaries/registry_ds/",
            "f57a3f78f5a31941  Dictionaries/registry_ds/reg_chr1.bin",
            "0923cbdf7545ac50  Dictionaries/registry_ds/reg_chr1.idx",
            "df1aff17472b46fc  Dictionaries/registry_ds/reg_chr2.bin",
            "0923cbdf7545ac50  Dictionaries/registry_ds/reg_chr2.idx",
            "3d344497202997f9  Dictionaries/registry_ds/variant_count.json",
            "dir  Genomes/",
            "dir  Genomes/hg38+ds/",
            "5228e2469ba9b172  Genomes/hg38+ds/chr1.enriched.fa",
            "a2bd555711c6e7ba  Genomes/hg38+ds/chr2.enriched.fa",
            "dir  Genomes/hg38+ds_INDELS/",
            "dir  Genomes/hg38+ds_INDELS/fake_ds_chr1/",
            "348a0d73b239b4a3  Genomes/hg38+ds_INDELS/fake_ds_chr1/fakechr1.fa",
            "dir  Genomes/hg38+ds_INDELS/fake_ds_chr2/",
            "0f767a34da5df5a5  Genomes/hg38+ds_INDELS/fake_ds_chr2/fakechr2.fa",
            "dir  Genomes/hg38/",
            "9fb230dc8fb2ed4d  Genomes/hg38/chr1.fa",
            "c5336dd0168ef6b3  Genomes/hg38/chr2.fa",
            "dir  VCFs/",
            "dir  VCFs/ds/",
            "25cf25d958babc38  VCFs/ds/ds.chr1.vcf.gz",
            "26a5b9b8e3b7cfc5  VCFs/ds/ds.chr2.vcf.gz",
            "dir  genome_library/",
            "dir  genome_library/NRG_2_hg38+ds/",
            "4983ff2823805684  genome_library/NRG_2_hg38+ds/.display_label",
            "425fad713f395627  genome_library/NRG_2_hg38+ds/.pam_build",
            "1dafb9f88a997704  genome_library/NRG_2_hg38+ds/NRG_chr1_1.bin",
            "fdb7a83262a278ee  genome_library/NRG_2_hg38+ds/NRG_chr2_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/",
            "270e5e2a33dbf6ed  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/NRG_fakechr1_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/",
            "5086fa72ea1b3534  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/NRG_fakechr2_1.bin",
            "dir  genome_library/NRG_2_hg38/",
            "425fad713f395627  genome_library/NRG_2_hg38/.pam_build",
            "dee1e9b039f12fe0  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "338fbfb5efc7044c  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
            "c442c9cdd4d80f67  samples.txt",
            "dir  samplesIDs/",
            "036a130b3364d826  samplesIDs/hg38_ds.samplesID.txt",
        ],
        "stdout": [
            "  emitting tiers with 2 parallel worker(s) (RAM-bounded; set CRISPRME_TIER_WORKERS to tune)",
            "Building indels index NRG_2_hg38+ds_INDELS...",
            "Building reference index NRG_2_hg38 (bMax 2, 2 thread(s))...",
            "Building variant index NRG_2_hg38+ds...",
            "Dictless tier emission complete: 2/2 chromosome(s) emitted (registry_<vcf>/ + genotypes_<vcf>/ siblings of the dicts).",
            "Emitting dictless tiers for 2 chromosome(s) from dictionaries_ds (databases: ds)...",
            "Enriching hg38 with ds (add-variants)...",
            "Index built: <W>/genome_library/NRG_2_hg38",
            "Variant index ready: NRG_2_hg38+ds (+ _INDELS). A later variant-aware search on hg38 with ds reuses it.",
            "Wrote variant-count manifest: 2 SNPs + 2 indels.",
            "[indel-snp] SNP-overlaid the fake-indel genome (2 flank bases IUPAC-coded).",
            "[indel-snp] building phased indel GT for 2 chromosome(s) with 2 parallel worker(s)",
            "[indel-snp] phased indel genotype store: 2 chromosome(s) -> indel_genotypes_ds",
            "overlay_indel_snps: 2 total SNP overlays",
            "overlay_indel_snps: chr1 -> 1 SNP overlays",
            "overlay_indel_snps: chr2 -> 1 SNP overlays",
        ],
        "stderr": [],
        "calls": [
            "add-variants <W>/VCFs/ds/ <W>/Genomes/hg38/ true",
            "index-genome hg38 <W>/Genomes/hg38/ <REPO>/PAMs/20bp-NRG-SpCas9.txt -bMax 2 -th 2",
            "index-genome hg38+ds <W>/Genomes/hg38+ds/ <REPO>/PAMs/20bp-NRG-SpCas9.txt -bMax 2 -th 2",
            "index-genome hg38+ds_INDELS/NRG_2_fakechr1 <W>/Genomes/hg38+ds_INDELS//fake_ds_chr1 <REPO>/PAMs/20bp-NRG-SpCas9.txt -bMax 2 -th 1",
            "index-genome hg38+ds_INDELS/NRG_2_fakechr2 <W>/Genomes/hg38+ds_INDELS//fake_ds_chr2 <REPO>/PAMs/20bp-NRG-SpCas9.txt -bMax 2 -th 1",
        ],
    },
    'rerun': {
        "rc": 0,
        "tree": [
            "dir  Dictionaries/",
            "dir  Dictionaries/dictionaries_ds/",
            "d1b69ccc5c7ccb40  Dictionaries/dictionaries_ds/my_dict_chr1.json.gz",
            "497425130c0deb8d  Dictionaries/dictionaries_ds/my_dict_chr2.json.gz",
            "dir  Dictionaries/genotypes_ds/",
            "c62654c4153c1576  Dictionaries/genotypes_ds/gt_chr1.bin",
            "acc7f66ace5d0084  Dictionaries/genotypes_ds/gt_chr1.idx",
            "c62654c4153c1576  Dictionaries/genotypes_ds/gt_chr2.bin",
            "acc7f66ace5d0084  Dictionaries/genotypes_ds/gt_chr2.idx",
            "dir  Dictionaries/indel_genotypes_ds/",
            "a3f726ac65dbc9e4  Dictionaries/indel_genotypes_ds/gt_indel_chr1.tsv.gz",
            "4355a46b19d348dc  Dictionaries/indel_genotypes_ds/gt_indel_chr1.tsv.gz.done",
            "8e3f9e931ce698b8  Dictionaries/indel_genotypes_ds/gt_indel_chr2.tsv.gz",
            "4355a46b19d348dc  Dictionaries/indel_genotypes_ds/gt_indel_chr2.tsv.gz.done",
            "dir  Dictionaries/log_indels_ds/",
            "5717b5bb7539d559  Dictionaries/log_indels_ds/logchr1.txt.gz",
            "0004c31cfbfba66d  Dictionaries/log_indels_ds/logchr2.txt.gz",
            "dir  Dictionaries/registry_ds/",
            "f57a3f78f5a31941  Dictionaries/registry_ds/reg_chr1.bin",
            "0923cbdf7545ac50  Dictionaries/registry_ds/reg_chr1.idx",
            "df1aff17472b46fc  Dictionaries/registry_ds/reg_chr2.bin",
            "0923cbdf7545ac50  Dictionaries/registry_ds/reg_chr2.idx",
            "3d344497202997f9  Dictionaries/registry_ds/variant_count.json",
            "dir  Genomes/",
            "dir  Genomes/hg38+ds/",
            "5228e2469ba9b172  Genomes/hg38+ds/chr1.enriched.fa",
            "a2bd555711c6e7ba  Genomes/hg38+ds/chr2.enriched.fa",
            "dir  Genomes/hg38+ds_INDELS/",
            "dir  Genomes/hg38+ds_INDELS/fake_ds_chr1/",
            "348a0d73b239b4a3  Genomes/hg38+ds_INDELS/fake_ds_chr1/fakechr1.fa",
            "dir  Genomes/hg38+ds_INDELS/fake_ds_chr2/",
            "0f767a34da5df5a5  Genomes/hg38+ds_INDELS/fake_ds_chr2/fakechr2.fa",
            "dir  Genomes/hg38/",
            "9fb230dc8fb2ed4d  Genomes/hg38/chr1.fa",
            "c5336dd0168ef6b3  Genomes/hg38/chr2.fa",
            "dir  VCFs/",
            "dir  VCFs/ds/",
            "25cf25d958babc38  VCFs/ds/ds.chr1.vcf.gz",
            "26a5b9b8e3b7cfc5  VCFs/ds/ds.chr2.vcf.gz",
            "dir  genome_library/",
            "dir  genome_library/NRG_2_hg38+ds/",
            "4983ff2823805684  genome_library/NRG_2_hg38+ds/.display_label",
            "425fad713f395627  genome_library/NRG_2_hg38+ds/.pam_build",
            "1dafb9f88a997704  genome_library/NRG_2_hg38+ds/NRG_chr1_1.bin",
            "fdb7a83262a278ee  genome_library/NRG_2_hg38+ds/NRG_chr2_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/",
            "270e5e2a33dbf6ed  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/NRG_fakechr1_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/",
            "5086fa72ea1b3534  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/NRG_fakechr2_1.bin",
            "dir  genome_library/NRG_2_hg38/",
            "425fad713f395627  genome_library/NRG_2_hg38/.pam_build",
            "dee1e9b039f12fe0  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "338fbfb5efc7044c  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
            "c442c9cdd4d80f67  samples.txt",
            "dir  samplesIDs/",
            "036a130b3364d826  samplesIDs/hg38_ds.samplesID.txt",
        ],
        "stdout": [
            "  2 chromosome(s) already emitted; resuming 0 remaining",
            "Dictless tier emission complete: 2/2 chromosome(s) emitted (registry_<vcf>/ + genotypes_<vcf>/ siblings of the dicts).",
            "Emitting dictless tiers for 2 chromosome(s) from dictionaries_ds (databases: ds)...",
            "Enriched genome already present: <W>/Genomes/hg38+ds",
            "Indels index already present: <W>/genome_library/NRG_2_hg38+ds_INDELS",
            "Reference index already present: <W>/genome_library/NRG_2_hg38",
            "Variant index already present: <W>/genome_library/NRG_2_hg38+ds",
            "Variant index ready: NRG_2_hg38+ds (+ _INDELS). A later variant-aware search on hg38 with ds reuses it.",
            "Wrote variant-count manifest: 2 SNPs + 2 indels.",
            "[indel-snp] phased indel genotype store: 2 chromosome(s) -> indel_genotypes_ds",
        ],
        "stderr": [],
        "calls": [],
    },
}


class SameTreeAsBefore(_Fixture):
    """(f): a successful build, and a re-run with everything in place, produce the
    same tree, file for file, and print the same lines, as before the fix."""

    maxDiff = None

    def describe(self, rc, out, err):
        normalize = ((self.work, "<W>"), (self.root, "<R>"), (REPO, "<REPO>"))

        def lines(text):
            for path, label in normalize:
                text = text.replace(path, label)
            # Not the INDELS pool's workers' lines: they print unflushed, and the pool's
            # exit terminates a worker that is still idle, so a line can go missing
            # (before this fix as after; its order is the pool's too, hence sorted).
            return sorted(ln for ln in text.splitlines()
                          if not ln.startswith(("Indexing INDELs in ",
                                                "Indexing ended for INDELs in ")))

        files = tree(self.work, normalize)
        return {"rc": rc, "tree": [f"{files[p]}  {p}" for p in sorted(files)],
                "stdout": lines(out), "stderr": lines(err), "calls": lines("\n".join(self.calls()))}

    def check(self, key, got):
        self.assertEqual(got, BEFORE[key])

    def test_a_reference_index(self):
        self.check("reference", self.describe(*self.run_bio(vcf=False, name="My index")))

    def test_a_variant_index_and_a_rerun_with_everything_in_place(self):
        self.check("variant", self.describe(*self.run_bio(name="My index")))
        os.remove(self.log)
        self.check("rerun", self.describe(*self.run_bio(name="My index")))


if __name__ == "__main__":
    if sys.argv[1:] == ["--record"]:  # print BEFORE for the code under test
        setUpModule()
        t = SameTreeAsBefore("test_a_reference_index")
        t.fresh()
        rec = {"reference": t.describe(*t.run_bio(vcf=False, name="My index"))}
        t.fresh()
        rec["variant"] = t.describe(*t.run_bio(name="My index"))
        os.remove(t.log)
        rec["rerun"] = t.describe(*t.run_bio(name="My index"))
        json.dump(rec, sys.stdout, indent=1)
        tearDownModule()
    else:
        unittest.main(verbosity=2)

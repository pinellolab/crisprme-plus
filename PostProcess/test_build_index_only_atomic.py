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
writes the layout the real one writes and can be told to fail, or to wait, before
one chromosome; the stub pigz compresses as pigz does and can be told to wait. A
wait ends on a kill, a signal or a release file, so every case is deterministic.
SameTreeAsBefore pins the successful builds to what the code before the fix
produced with the same stub, file for file.

Stdlib only: crisprme.py imports Bio.Seq at module scope, so a stub Bio package goes
on PYTHONPATH (build-index-only never uses it).
"""

import gzip
import hashlib
import json
import os
import shutil
import signal
import stat
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
waits there, once $STUB_MARK is made, until $STUB_RELEASE exists. add-variants,
with $STUB_GZ_BLOCKED, leaves a directory where the last dictionary's .gz would go.
"""
import json, os, shutil, sys, time


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
        step("index-genome %s %s" % (name, fa))
        write(os.path.join(out, "%s_%s_1.bin" % (pam, fa[: -len(".fa")])),
              "%s %s %s\n" % (pam, nb, sequence(os.path.join(gdir, fa))))
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
# the same. STUB_SIGNAL_AT_CLEANUP=1: a SIGTERM comes just as the build's cleanup
# starts, before its first pthread_sigmask holds signals back.
SHIM = r'''
import errno, os, signal

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
    _sigmask, _sent = signal.pthread_sigmask, []

    def pthread_sigmask(how, mask):
        handler = signal.getsignal(signal.SIGTERM)
        if not _sent and how == signal.SIG_BLOCK and callable(handler):
            _sent.append(1)
            handler(signal.SIGTERM, None)  # what a SIGTERM that comes now runs
        return _sigmask(how, mask)

    signal.pthread_sigmask = pthread_sigmask
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

    def start(self, vcf=True, name=None, keep_hup=False, pigz=True, shim=False, command=None,
              **stub):
        """Launch build-index-only in a session of its own (so its whole process
        group can be killed, as an OOM kill or a container stop would). Without
        pigz, no pigz is on PATH; shim loads SHIM in each Python process of the
        run; command runs that crisprme.py command line instead."""
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
        CRISPR-Bulge env is here)."""
        _write(self.path("guides.txt"), GENOME["chr1"][:20] + "NNN\n")
        _write(self.path("vcfs.txt"), "ds\n")
        command = ["complete-search", "--genome", self.path("Genomes", "hg38"),
                   "--vcf", self.path("vcfs.txt"), "--samplesID", self.path("samples.txt"),
                   "--guide", self.path("guides.txt"), "--pam", PAM, "--mm", "4",
                   "--bDNA", "1", "--bRNA", "1", "--output", "search", "--thread", "2"]
        return self.finish(self.start(command=command, CRISPRME_SKIP_SCORER_ENV="1"))

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


# What each index looks like once built, and the step of the stub's that builds its
# second chromosome (so stopping there leaves the first one built).
WHOLE = {
    REF: (False, "index-genome hg38 chr2.fa", ["NRG_chr1_1.bin", "NRG_chr2_1.bin"]),
    VAR: (True, "index-genome hg38+ds chr2.enriched.fa",
          ["NRG_chr1.enriched_1.bin", "NRG_chr2.enriched_1.bin"]),
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
                # a SIGKILL leaves the hidden build directory, which nothing reads
                left = self.leftovers()
                self.assertEqual(len(left), 1)
                rc, out, err = self.run_bio(vcf=vcf)
                self.assertEqual(rc, 0, err)
                self.assertEqual(self.bins(name), whole)
                self.assertEqual(self.leftovers(), left)

    def test_an_index_is_built_in_a_hidden_directory_of_genome_library(self):
        """Where the index is renamed from: no rename then crosses a mount point, and
        the index is made with what genome_library gives a new directory."""
        proc = self.start(vcf=False, STUB_HANG=WHOLE[REF][1])
        self.wait_for_mark(proc)
        left = self.leftovers()
        self.assertEqual([os.path.dirname(d) for d in left], [self.path("genome_library")])
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
        self.assertEqual(self.leftovers(), left)  # (the build directory a SIGKILL left)
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


class AnotherRunsBuildIsKept(_Fixture):
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
        """One that comes once the index is in place, as the cleanup begins and before
        it holds signals back (SHIM sends it there): the run removes its build
        directory and ends by that signal, with no traceback."""
        rc, out, err = self.run_bio(vcf=False, shim=True, STUB_SIGNAL_AT_CLEANUP="1")
        self.assertEqual(rc, -signal.SIGTERM, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(self.leftovers(), [])
        self.assertEqual(self.bins(REF), WHOLE[REF][2])

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
            self.assertIn("run complete-search, which takes a dead run's lock over", err)
        else:
            self.assertIn(f"Error: {output} is locked by the run that builds it: {lock} "
                          "has no owner record", err)
        guarded = " '" + self.path("Genomes", "hg38+ds_INDELS") + "'" \
            if output.endswith("Genomes/hg38+ds") else ""
        self.assertIn(f"rm -rf '{lock}' '{output}'{guarded}\n", err)
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
            "836398cfb63675bc  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "400e714b6d577f71  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
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
            "04fd219b01b57b8e  genome_library/NRG_2_hg38+ds/NRG_chr1.enriched_1.bin",
            "8e3da51d843e800e  genome_library/NRG_2_hg38+ds/NRG_chr2.enriched_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/",
            "e3c79f97ace30cbd  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/NRG_fakechr1_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/",
            "82ec6931a4d9ff3d  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/NRG_fakechr2_1.bin",
            "dir  genome_library/NRG_2_hg38/",
            "425fad713f395627  genome_library/NRG_2_hg38/.pam_build",
            "836398cfb63675bc  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "400e714b6d577f71  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
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
            "04fd219b01b57b8e  genome_library/NRG_2_hg38+ds/NRG_chr1.enriched_1.bin",
            "8e3da51d843e800e  genome_library/NRG_2_hg38+ds/NRG_chr2.enriched_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/",
            "e3c79f97ace30cbd  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr1/NRG_fakechr1_1.bin",
            "dir  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/",
            "82ec6931a4d9ff3d  genome_library/NRG_2_hg38+ds_INDELS/NRG_2_fakechr2/NRG_fakechr2_1.bin",
            "dir  genome_library/NRG_2_hg38/",
            "425fad713f395627  genome_library/NRG_2_hg38/.pam_build",
            "836398cfb63675bc  genome_library/NRG_2_hg38/NRG_chr1_1.bin",
            "400e714b6d577f71  genome_library/NRG_2_hg38/NRG_chr2_1.bin",
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

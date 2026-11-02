#!/usr/bin/env python3
"""A build-index-only worker that dies ends the run, with a message; it does not hang it.

The dict-less tiers and the phased indel genotype store are written by forked workers,
one chromosome each. They went through a multiprocessing.Pool, whose starmap does not
see a worker die: the OOM killer took one (the tiers take ~150-180 GB a chromosome at
panel scale), the pool started an idle worker in its place and waited for the result for
good. A web Settings variant build, which writes the tiers since 2.8.1, then read
"running" for ever, with every later job queued behind it. Now a dead worker sends the
chromosomes left to one worker at a time, and one that dies then too ends the run with
an ERROR that says memory is the likely cause.

These run the real crisprme.py build-index-only against the stub engine of
test_build_index_only_atomic, with a worker that kills itself (SIGKILL, as the OOM killer
does) as it opens its output, and a bound on how long the run may take.
"""

import os
import shutil
import signal
import time
import unittest

import test_build_index_only_atomic as bio

bio.SHIM += r'''
if os.environ.get("STUB_WRITE_KILL"):
    _open_killed = builtins.open

    def open(file, mode="r", *args, **kwargs):
        if "w" in mode and os.environ["STUB_WRITE_KILL"] in str(file):
            _once = os.environ.get("STUB_WRITE_KILL_ONCE")
            if not _once or not os.path.exists(_once):
                if _once:
                    _open_killed(_once, "w").close()
                os.kill(os.getpid(), signal.SIGKILL)
        return _open_killed(file, mode, *args, **kwargs)

    builtins.open = open
'''

setUpModule = bio.setUpModule
tearDownModule = bio.tearDownModule

# (a hang is what this catches: far more than the whole run takes)
BOUND = 90

# each step: the write of its second chromosome, its directory, the variable that sets its
# workers, and what its ERROR says
STEPS = (
    ("tiers", "/registry_ds/reg_chr2.bin", "registry_ds", "CRISPRME_TIER_WORKERS",
     "Emitting the dict-less tiers of ds failed: a worker process died"),
    ("indel genotype store", "/indel_genotypes_ds/gt_indel_chr2", "indel_genotypes_ds",
     "CRISPRME_INDELGT_WORKERS",
     "Building the phased indel genotype store of ds failed: a worker process died"),
)


def group_members(pgid):
    """The running processes of a process group (a worker that has ended stays a zombie
    until whoever took it on reaps it, which a container's first process may never do)."""
    found = []
    for pid in filter(str.isdigit, os.listdir("/proc")):
        try:
            with open(f"/proc/{pid}/stat") as fh:
                state, _ppid, group = fh.read().rsplit(")", 1)[1].split()[:3]
        except OSError:
            continue
        if int(group) == pgid and state != "Z":
            found.append(int(pid))
    return found


class AWorkerThatDies(bio._Fixture):
    def run_within_bound(self, **stub):
        return self.finish(self.start(shim=True, **stub), timeout=BOUND)

    def test_on_every_try_ends_the_run_with_a_message_and_the_mark_left(self):
        self.run_bio()
        clean = self.outputs()
        for step, write, where, knob, says in STEPS:
            for workers in ("2", "1"):
                with self.subTest(step=step, workers=workers):
                    self.fresh()
                    rc, out, err = self.run_within_bound(STUB_WRITE_KILL=write,
                                                         **{knob: workers})
                    self.assertEqual(rc, 1, out + err)
                    self.assertTrue(err.startswith(f"Error: {says}"), err)
                    self.assertIn("most likely killed for want of memory", err)
                    self.assertIn(knob, err)
                    directory = self.path("Dictionaries", where)
                    self.assertIn(f"{directory} stays marked unfinished", err)
                    self.assertIn("Run build-index-only again with more memory", err)
                    self.assertTrue(os.path.isdir(os.path.join(directory, ".building")))
                    self.assertFalse(os.path.exists(self.lib(bio.VAR)))  # it ended there
                    # (one try, or the one-at-a-time one after it; stdout says so)
                    self.assertEqual(f"{knob} below" in out, workers == "2", out)
                    # a search refuses what it left
                    rc, out, err = self.run_search()
                    self.assertEqual(rc, 1, out + err)
                    self.assertIn(f"Error: {directory} is unfinished: ", err)
                    # with memory again, a run ends what it began, as a run from nothing does
                    rc, out, err = self.run_bio()
                    self.assertEqual(rc, 0, out + err)
                    self.assertEqual(self.outputs(), clean)

    def test_the_way_the_message_gives_to_go_without_what_it_could_not_write(self):
        args = ["build-index-only", "--genome", self.path("Genomes", "hg38"), "--pam", bio.PAM,
                "--bDNA", "1", "--bRNA", "1", "--thread", "2", "--path", self.work,
                "--vcf", self.path("VCFs", "ds")]
        # the tiers: removed, and a run with no --samplesID: the search reads the dictionaries
        rc, out, err = self.run_within_bound(STUB_WRITE_KILL=STEPS[0][1])
        self.assertEqual(rc, 1, out + err)
        registry, genotypes = (self.path("Dictionaries", d) for d in ("registry_ds", "genotypes_ds"))
        self.assertIn(f"remove them (rm -rf '{registry}' '{genotypes}') and run build-index-only "
                      "from the command line with no --samplesID", err)
        shutil.rmtree(registry)
        shutil.rmtree(genotypes)
        rc, out, err = self.finish(self.start(command=args))
        self.assertEqual(rc, 0, out + err)
        self.assertFalse(os.path.exists(registry))
        self.assert_searched()
        # the indel genotype store: removed, and a run with the indel + SNP pass off
        self.fresh()
        rc, out, err = self.run_within_bound(STUB_WRITE_KILL=STEPS[1][1])
        self.assertEqual(rc, 1, out + err)
        store = self.path("Dictionaries", "indel_genotypes_ds")
        self.assertIn(f"remove it (rm -rf '{store}') and run build-index-only with "
                      "CRISPRME_INDEL_SNP=0", err)
        shutil.rmtree(store)
        rc, out, err = self.run_bio(CRISPRME_INDEL_SNP="0")
        self.assertEqual(rc, 0, out + err)
        self.assertFalse(os.path.exists(store))
        self.assert_searched()

    def test_once_leaves_the_chromosomes_left_to_one_worker_at_a_time(self):
        self.run_bio()
        clean = self.outputs()
        for step, write, where, knob, _says in STEPS:
            with self.subTest(step=step):
                self.fresh()
                rc, out, err = self.run_within_bound(
                    STUB_WRITE_KILL=write, STUB_WRITE_KILL_ONCE=self.path("killed_once"),
                    **{knob: "2"})
                self.assertEqual(rc, 0, out + err)
                self.assertEqual(err, "")
                self.assertIn("a worker process died, most likely killed for want of memory; "
                              f"set {knob} below 2 to start with fewer", out)
                self.assertTrue(os.path.exists(self.path("killed_once")))
                self.assertEqual(self.outputs(), clean)
                self.assertFalse(os.path.exists(self.path("Dictionaries", where, ".building")))
                self.assert_searched()

    def test_no_worker_outlives_a_run_killed_outright(self):
        # one worker writes chr1, ends it and waits for its next job; the other holds at
        # its write of chr2 (a run killed by SIGKILL, or by a SIGTERM with no handler)
        if not os.path.exists("/proc/self/stat"):
            self.skipTest("no /proc here")
        proc = self.start(shim=True, STUB_WRITE_HANG="/registry_ds/reg_chr2.bin",
                          CRISPRME_TIER_WORKERS="2")
        self.wait_for_mark(proc)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(timeout=30)
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline:
            if not group_members(proc.pid):
                return
            time.sleep(0.2)
        self.fail(f"workers of the run still there 30 s after it was killed: "
                  f"{group_members(proc.pid)}")


if __name__ == "__main__":
    unittest.main()

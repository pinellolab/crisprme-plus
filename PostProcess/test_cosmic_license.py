"""Tests for the COSMIC licence gate (PostProcess/cosmic_license.py)."""

import contextlib
import gzip
import json
import os
import stat
import tempfile
import unittest
from unittest import mock

import cosmic_license as cl


_BED_ROWS = [
    "chr1\t100\t200\tENCODE_cCRE_1\n",
    "chr1\t150\t250\tTP53_COSMIC\n",
    "chr1\t300\t400\tDHS_peak_3\n",
    "chr2\t500\t600\tBRCA1_COSMIC\n",
    "chr2\t700\t800\tgene_X_gencode\n",
]


@contextlib.contextmanager
def a_link_at_each_new_temporary_name(target, wanted):
    """Another user's step, which anyone who can write the folder can take: right after this
    process makes a file whose name ``wanted`` takes (os.open with O_CREAT|O_EXCL, as mkstemp
    and safe_files make a temporary file of their own), that file renamed away
    (<name>.moved) and a link to ``target`` put under its name, so every later step by that
    name meets the link. Taken where it would win the race, deterministically (os.open
    patched). Yields the names so taken."""
    real_open = os.open
    taken = []

    def open_(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        name = os.path.basename(os.fsdecode(path))
        if flags & os.O_CREAT and flags & os.O_EXCL and wanted(name):
            os.rename(path, os.fsdecode(path) + ".moved", src_dir_fd=dir_fd, dst_dir_fd=dir_fd)
            os.symlink(target, path, dir_fd=dir_fd)
            taken.append(name)
        return fd

    with mock.patch.object(os, "open", open_):
        yield taken


class TestCosmicFlag(unittest.TestCase):
    def test_default_disabled_and_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            # missing flag file -> disabled (licence-safe default)
            self.assertFalse(cl.cosmic_enabled(d))
            cl.set_cosmic_license(d, True, "attested via test")
            self.assertTrue(cl.cosmic_enabled(d))
            self.assertEqual(cl.get_cosmic_license(d).get("attestation"), "attested via test")
            # a regular file, readable by all as the record was before (mkstemp makes 0600)
            st = os.lstat(cl.license_path(d))
            self.assertTrue(stat.S_ISREG(st.st_mode))
            self.assertEqual(stat.S_IMODE(st.st_mode), 0o644)
            cl.set_cosmic_license(d, False)
            self.assertFalse(cl.cosmic_enabled(d))
            self.assertEqual(os.listdir(d), [cl.COSMIC_LICENSE_FILE])  # (no temporary left)

    def test_a_link_put_at_the_temporary_name_is_neither_changed_nor_followed(self):
        """The review of wave 1c-a, A1: Annotations/ may be writable by a group, and anyone
        who can write it can rename the record's temporary file away as soon as it is made and
        put a link in its place. The save set that file's mode by its name (os.chmod), which
        followed the link: a file of the user's it named, the email settings with their app
        password say (0600), became readable by all (0644). Its mode is set through its
        descriptor now: the file the link names keeps its mode and its content, and the record
        is in the file made for it (renamed away), readable by all."""
        with tempfile.TemporaryDirectory() as top:
            ann = os.path.join(top, "Annotations")
            os.makedirs(ann)
            os.chmod(ann, 0o775)  # (writable by the group, as under a umask of 002)
            secret = os.path.join(top, "email.json")
            with open(secret, "w") as fh:
                fh.write('{"password": "APP-PASSWORD"}\n')
            os.chmod(secret, 0o600)
            with a_link_at_each_new_temporary_name(
                    secret, lambda n: n.startswith(cl.COSMIC_LICENSE_FILE + ".")) as taken:
                cl.set_cosmic_license(ann, True, "attested via test")
            self.assertEqual(len(taken), 1, taken)  # (the other user's step was taken)
            self.assertEqual(stat.S_IMODE(os.stat(secret).st_mode), 0o600)
            with open(secret) as fh:
                self.assertEqual(fh.read(), '{"password": "APP-PASSWORD"}\n')
            moved = os.path.join(ann, taken[0] + ".moved")
            self.assertEqual(stat.S_IMODE(os.stat(moved).st_mode), 0o644)
            with open(moved) as fh:
                self.assertEqual(json.load(fh), {"version": 1, "enabled": True,
                                                 "attestation": "attested via test"})

    def test_the_record_text(self):
        # (the web Settings write this text through safe_files: pages_utils.set_cosmic_license)
        self.assertEqual(json.loads(cl.license_record(True, "x")),
                         {"version": 1, "enabled": True, "attestation": "x"})
        self.assertEqual(json.loads(cl.license_record(0)),
                         {"version": 1, "enabled": False, "attestation": ""})

    def test_corrupt_flag_file_is_disabled(self):
        with tempfile.TemporaryDirectory() as d:
            with open(cl.license_path(d), "w") as fh:
                fh.write("not json{")
            self.assertFalse(cl.cosmic_enabled(d))


class TestCosmicRowDetection(unittest.TestCase):
    def test_is_cosmic_row(self):
        self.assertTrue(cl.is_cosmic_row("chr1\t150\t250\tTP53_COSMIC\n"))
        self.assertTrue(cl.is_cosmic_row("chr2\t500\t600\tBRCA1_COSMIC"))
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tENCODE_cCRE_1\n"))
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tgene_gencode\n"))
        # the tag only counts in the name column, not coordinates/other cols
        self.assertFalse(cl.is_cosmic_row("chr1\t100\t200\tDHS\tnote_COSMIC_ish\n"))


class TestStripCosmic(unittest.TestCase):
    def _write(self, path, rows, gz=False):
        op = gzip.open if gz else open
        with op(path, "wt") as fh:
            fh.writelines(rows)

    def test_strip_plain_bed(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed")
            out = os.path.join(d, "ann.nocosmic.bed")
            self._write(src, _BED_ROWS)
            self.assertTrue(cl.bed_has_cosmic(src))
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual(dropped, 2)   # TP53_COSMIC + BRCA1_COSMIC
            self.assertEqual(kept, 3)
            with open(out) as fh:
                body = fh.read()
            self.assertNotIn("_COSMIC", body)
            self.assertIn("ENCODE_cCRE_1", body)
            self.assertIn("gene_X_gencode", body)
            self.assertFalse(cl.bed_has_cosmic(out))

    def test_strip_gzip_bed(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed.gz")
            out = os.path.join(d, "ann.nocosmic.bed.gz")
            self._write(src, _BED_ROWS, gz=True)
            self.assertTrue(cl.bed_has_cosmic(src))
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual((kept, dropped), (3, 2))
            with gzip.open(out, "rt") as fh:
                body = fh.read()
            self.assertNotIn("_COSMIC", body)

    def test_strip_keeps_comments(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed")
            out = os.path.join(d, "out.bed")
            self._write(src, ["#track cosmic annotations\n"] + _BED_ROWS)
            kept, dropped = cl.strip_cosmic(src, out)
            self.assertEqual(dropped, 2)
            with open(out) as fh:
                self.assertIn("#track cosmic annotations", fh.read())


def _concurrent_strip_worker(payload):
    import cosmic_license as _cl
    src, out, iters = payload
    for _ in range(iters):
        _cl.strip_cosmic(src, out)
    return True


class TestStripCosmicConcurrency(unittest.TestCase):
    def test_concurrent_strip_same_output_does_not_race(self):
        """Two+ complete-search jobs stripping the SAME annotation to the SAME out_bed
        must not race on the temp file. Pre-fix both wrote '<out>.tmp' and os.replace
        intermittently hit FileNotFoundError; a per-process temp name fixes it (the
        final os.replace to out_bed is atomic and the content identical, so
        last-writer-wins)."""
        import multiprocessing as mp
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed.gz")
            out = os.path.join(d, ".nocosmic.ann.bed.gz")
            with gzip.open(src, "wt") as fh:
                fh.writelines(_BED_ROWS * 300)  # enough size that writes overlap
            ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() \
                else mp.get_context()
            with ctx.Pool(4) as pool:
                res = pool.map(_concurrent_strip_worker, [(src, out, 20)] * 4)
            self.assertTrue(all(res))  # no FileNotFoundError in any worker
            with gzip.open(out, "rt") as fh:
                body = fh.read()
            self.assertNotIn("_COSMIC", body)      # output still valid + cosmic-free
            self.assertIn("ENCODE_cCRE_1", body)

    def test_never_through_a_name_already_there(self):
        """The strip writes beside its output, in Annotations/ of a data folder others may
        write (a lab's), under a name made of its process id. Another user could have put a
        link there first, to a file of the user's: the strip then overwrote that file with
        the annotation's lines (2.10.0). It makes its own file now (O_EXCL): a name that
        already stands fails it (the caller then strips into a folder of the run's own), and
        the file the link names is left as it was."""
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "ann.bed")
            with open(src, "w") as fh:
                fh.writelines(_BED_ROWS)
            victim = os.path.join(d, "victim.txt")
            with open(victim, "w") as fh:
                fh.write("mine\n")
            for out in (os.path.join(d, ".nocosmic.ann.bed"),
                        os.path.join(d, ".nocosmic.ann.bed.gz")):
                with self.subTest(out=os.path.basename(out)):
                    link = out + ".%d.tmp" % os.getpid()
                    os.symlink(victim, link)
                    with self.assertRaises(FileExistsError):
                        cl.strip_cosmic(src, out)
                    with open(victim) as fh:
                        self.assertEqual(fh.read(), "mine\n")
                    self.assertFalse(os.path.exists(out))
                    os.remove(link)
                    # and with no such name, the strip as before, its mode the umask's
                    old = os.umask(0o022)
                    try:
                        self.assertEqual(cl.strip_cosmic(src, out), (3, 2))
                    finally:
                        os.umask(old)
                    self.assertEqual(os.stat(out).st_mode & 0o777, 0o644)
                    opener = gzip.open if out.endswith(".gz") else open
                    with opener(out, "rt") as fh:
                        self.assertEqual(fh.read(), "".join(
                            r for r in _BED_ROWS if "_COSMIC" not in r))
                    self.assertEqual(sorted(os.listdir(d)), sorted(
                        ["ann.bed", "victim.txt"]
                        + [os.path.basename(p) for p in
                           (os.path.join(d, ".nocosmic.ann.bed"),
                            os.path.join(d, ".nocosmic.ann.bed.gz")) if os.path.exists(p)]))

    def test_tmp_name_is_not_the_racy_fixed_suffix(self):
        """Guard the fix: the temp must not be the fixed '<out_bed>.tmp' that races."""
        import inspect
        src = inspect.getsource(cl.strip_cosmic)
        self.assertNotIn('out_bed + ".tmp"', src)


if __name__ == "__main__":
    unittest.main()

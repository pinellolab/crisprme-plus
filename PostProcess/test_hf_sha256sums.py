#!/usr/bin/env python3
"""``scripts/hf_sha256sums.py`` must write a SHA256SUMS that ``sha256sum -c`` accepts.

The checksums come from the Hub's LFS metadata, so the failure modes worth pinning
are the ones that would publish a SHA256SUMS which verifies nothing or the wrong
thing: an empty listing read as a pass, a second page of the listing dropped, a
non-archive or a subfolder entry written as a line, or a staged replacement that
silently fails to take effect.

Hermetic: ``urllib.request.urlopen`` is replaced, no network.
"""

import hashlib
import importlib.util
import io
import json
import os
import subprocess
import tempfile
import unittest
from unittest import mock

_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "hf_sha256sums.py")
_spec = importlib.util.spec_from_file_location("hf_sha256sums", _SCRIPT)
sums_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sums_mod)

A = "a" * 64
B = "b" * 64
C = "c" * 64


class _Resp(io.BytesIO):
    def __init__(self, payload, link=None):
        super().__init__(json.dumps(payload).encode())
        self.headers = {"Link": link} if link else {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _file(path, oid):
    return {"type": "file", "path": path, "size": 1, "lfs": {"oid": oid, "size": 1}}


def _pages(*pages):
    """urlopen stand-in that serves ``pages`` in order, linking each to the next."""
    calls = []

    def fake(req, timeout=None):
        i = len(calls)
        calls.append(req.full_url)
        link = f'<https://huggingface.co/next/{i + 1}>; rel="next"' if i + 1 < len(pages) else None
        return _Resp(pages[i], link)

    return fake, calls


def _run(argv, pages):
    fake, calls = _pages(*pages)
    out, err = io.StringIO(), io.StringIO()
    with mock.patch.object(sums_mod.urllib.request, "urlopen", fake), \
            mock.patch("sys.stdout", out), mock.patch("sys.stderr", err):
        rc = sums_mod.main(argv)
    return rc, out.getvalue(), err.getvalue(), calls


class Listing(unittest.TestCase):
    def test_sha256sum_format_sorted_by_name(self):
        rc, out, _, _ = _run(["o/r"], [[_file("indexes/b.tar.gz", B), _file("indexes/a.tar.gz", A)]])
        self.assertEqual(rc, 0)
        self.assertEqual(out, f"{A}  a.tar.gz\n{B}  b.tar.gz\n")

    def test_every_page_is_read(self):
        _, out, _, calls = _run(["o/r"], [[_file("indexes/a.tar.gz", A)], [_file("indexes/b.tar.gz", B)]])
        self.assertEqual(len(calls), 2)
        self.assertIn("b.tar.gz", out)

    def test_only_archives_directly_under_the_folder(self):
        _, out, _, _ = _run(["o/r"], [[
            _file("indexes/a.tar.gz", A),
            _file("indexes/README.md", B),
            _file("indexes/SHA256SUMS", B),
            {"type": "directory", "path": "indexes/sub.tar.gz"},
        ]])
        self.assertEqual(out, f"{A}  a.tar.gz\n")

    def test_folder_is_part_of_the_request(self):
        _, _, _, calls = _run(["o/r", "--path", "build-inputs"], [[_file("build-inputs/k.tar.gz", A)]])
        self.assertIn("/tree/main/build-inputs?expand=true", calls[0])


class Refusals(unittest.TestCase):
    def test_empty_listing_is_an_error_not_an_empty_file(self):
        with self.assertRaises(SystemExit) as cm:
            _run(["o/r"], [[]])
        self.assertIn("empty SHA256SUMS", str(cm.exception))

    def test_archive_without_lfs_checksum_is_an_error(self):
        entry = {"type": "file", "path": "indexes/a.tar.gz", "size": 1}
        with self.assertRaises(SystemExit) as cm:
            _run(["o/r"], [[entry]])
        self.assertIn("a.tar.gz", str(cm.exception))

    def test_malformed_staged_line_is_an_error(self):
        with tempfile.NamedTemporaryFile("w", suffix=".sha256", delete=False) as fh:
            fh.write("deadbeef  a.tar.gz\n")
        try:
            with self.assertRaises(SystemExit):
                sums_mod.read_staged(fh.name)
        finally:
            os.unlink(fh.name)


class Staged(unittest.TestCase):
    def _staged(self, text):
        fh = tempfile.NamedTemporaryFile("w", suffix=".sha256", delete=False)
        fh.write(text)
        fh.close()
        self.addCleanup(os.unlink, fh.name)
        return fh.name

    def test_staged_entry_replaces_and_adds(self):
        staged = self._staged(f"{C}  a.tar.gz\n{C} *c.tar.gz\n")
        _, out, err, _ = _run(["o/r", "--staged", staged],
                              [[_file("indexes/a.tar.gz", A), _file("indexes/b.tar.gz", B)]])
        self.assertEqual(out, f"{C}  a.tar.gz\n{B}  b.tar.gz\n{C}  c.tar.gz\n")
        self.assertIn("replaces a.tar.gz", err)
        self.assertIn("adds c.tar.gz", err)

    def test_staged_path_prefix_is_reduced_to_the_name(self):
        staged = self._staged(f"{C}  /scratch/stage/indexes/a.tar.gz\n")
        _, out, _, _ = _run(["o/r", "--staged", staged], [[_file("indexes/a.tar.gz", A)]])
        self.assertEqual(out, f"{C}  a.tar.gz\n")


class Sha256sumAcceptsIt(unittest.TestCase):
    def test_round_trip_through_sha256sum_or_shasum(self):
        tool = None
        for cand in (["sha256sum"], ["shasum", "-a", "256"]):
            try:
                subprocess.run(cand + ["--version"], capture_output=True, check=True)
                tool = cand
                break
            except (OSError, subprocess.CalledProcessError):
                continue
        if tool is None:
            self.skipTest("neither sha256sum nor shasum available")
        with tempfile.TemporaryDirectory() as d:
            sums = {}
            for name, data in (("x.tar.gz", b"one"), ("y+z.tar.gz", b"two")):
                with open(os.path.join(d, name), "wb") as fh:
                    fh.write(data)
                sums[name] = hashlib.sha256(data).hexdigest()
            with open(os.path.join(d, "SHA256SUMS"), "w") as fh:
                fh.write(sums_mod.render(sums))
            res = subprocess.run(tool + ["-c", "SHA256SUMS"], cwd=d, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stdout + res.stderr)


if __name__ == "__main__":
    unittest.main()

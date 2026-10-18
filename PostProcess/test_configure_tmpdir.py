"""Unit test for crisprme.py::_configure_tmpdir (the $CRISPRME_TMPDIR knob).

crisprme.py runs its command dispatch at import time, so it cannot be imported in a
unit test; following test_indel_pam_creation, the function's source is extracted and
exec'd in a controlled namespace instead.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
CRISPRME = os.path.join(HERE, "..", "crisprme.py")


def _extract_func(name):
    src = open(CRISPRME).read()
    out, grab = [], False
    for ln in src.splitlines(keepends=True):
        if ln.startswith(f"def {name}("):
            grab = True
        elif grab and ln.startswith("def ") and not ln.startswith(f"def {name}("):
            break
        if grab:
            out.append(ln)
    if not out:
        raise AssertionError(f"could not extract {name} from crisprme.py")
    return "".join(out)


class TestConfigureTmpdir(unittest.TestCase):
    def setUp(self):
        ns = {"os": os, "sys": sys}
        exec(_extract_func("_configure_tmpdir"), ns)
        self.fn = ns["_configure_tmpdir"]
        self._saved = (os.environ.get("CRISPRME_TMPDIR"),
                       os.environ.get("TMPDIR"), tempfile.tempdir)

    def tearDown(self):
        c, t, td = self._saved
        for k, v in (("CRISPRME_TMPDIR", c), ("TMPDIR", t)):
            os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
        tempfile.tempdir = td

    def test_sets_tmpdir_and_creates_missing_dir(self):
        with tempfile.TemporaryDirectory() as d:
            sub = os.path.join(d, "crisprme_tmp")  # does not exist yet -> created
            os.environ["CRISPRME_TMPDIR"] = sub
            os.environ.pop("TMPDIR", None)
            self.fn()
            self.assertTrue(os.path.isdir(sub))
            self.assertEqual(os.environ.get("TMPDIR"), sub)
            self.assertEqual(tempfile.tempdir, sub)

    def test_noop_when_unset(self):
        os.environ.pop("CRISPRME_TMPDIR", None)
        os.environ.pop("TMPDIR", None)
        tempfile.tempdir = None
        self.fn()
        self.assertIsNone(os.environ.get("TMPDIR"))
        self.assertIsNone(tempfile.tempdir)

    def test_unwritable_dir_warns_and_falls_back(self):
        os.environ["CRISPRME_TMPDIR"] = "/proc/nonexistent/cannot_create"
        os.environ.pop("TMPDIR", None)
        tempfile.tempdir = None
        self.fn()  # must not raise; must not set TMPDIR to an unusable dir
        self.assertIsNone(os.environ.get("TMPDIR"))


if __name__ == "__main__":
    unittest.main()

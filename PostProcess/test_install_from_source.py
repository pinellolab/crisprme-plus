#!/usr/bin/env python3
"""install_from_source.sh refuses macOS before it does anything else.

The source build is Linux-only: environment.yml has packages with no macOS build, and
CRISPRitz's buildTST needs GCC's parallel STL, which no macOS toolchain provides. The
script used to start regardless and fail later, on a solver or compiler error. It now
stops at once on Darwin, with one line pointing to Docker.

Run with a stub `uname` on PATH and no conda env, so that without the guard the script
stops at its next check (the activated env), never at the build; a stub `git` would
fail loudly if it ever got that far.
"""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "install_from_source.sh"


class TheSourceBuildRefusesMacOS(unittest.TestCase):
    def run_script(self, system):
        with tempfile.TemporaryDirectory() as d:
            for name, body in (("uname", f"echo {system}"), ("git", "echo git >&2; exit 97")):
                stub = Path(d) / name
                stub.write_text(f"#!/bin/sh\n{body}\n")
                stub.chmod(0o755)
            env = {k: v for k, v in os.environ.items() if k != "CONDA_PREFIX"}
            env["PATH"] = f"{d}{os.pathsep}{env.get('PATH', '')}"
            return subprocess.run(["bash", str(SCRIPT)], env=env, capture_output=True,
                                  text=True)

    def test_darwin_stops_with_one_line_pointing_to_docker(self):
        out = self.run_script("Darwin")
        self.assertEqual(out.returncode, 1, out.stderr)
        self.assertEqual(out.stdout, "")
        lines = out.stderr.splitlines()
        self.assertEqual(len(lines), 1, out.stderr)
        self.assertIn("Docker", lines[0])
        self.assertIn("README 1.2", lines[0])

    def test_linux_goes_on_to_the_conda_env_check(self):
        out = self.run_script("Linux")
        self.assertNotEqual(out.returncode, 0)
        self.assertIn("CONDA_PREFIX", out.stderr)
        self.assertNotIn("Docker", out.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""Results-page links, found by the 2.8.1 web walkthrough:

* the job status page held the "view-results" link twice (one id, two elements).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_web_results_links.py' -v
"""

import os
import shutil
import sqlite3
import sys
import tempfile
import time
import unittest
import warnings
from unittest.mock import patch

import pandas as pd

try:  # page modules need the web stack (Dash); the light unit-test CI does not have it
    import dash  # noqa: F401
    from dash.development.base_component import Component
except ImportError:  # pragma: no cover
    raise unittest.SkipTest("dash not installed: web page tests skipped")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["crisprme.py", "--version"]  # importing the pages starts the app object

import pages.load_page as lp  # noqa: E402
import pages.results_page as rp  # noqa: E402

GUIDE = "CTAACAGTTGCTTTTATCACNNN"
CRITERIA = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")


def _walk(component):
    """Every component in a Dash layout."""
    yield component
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for child in children:
            yield from _walk(child)
    elif isinstance(children, Component):
        yield from _walk(children)


class TestStatusPageIds(unittest.TestCase):
    def test_no_id_twice(self):
        for name, layout in (("load_page", lp.load_page("x")),
                             ("load_page_assembly", lp.load_page_assembly("x"))):
            with self.subTest(page=name):
                ids = [c.id for c in _walk(layout) if getattr(c, "id", None)]
                self.assertEqual(sorted(i for i in set(ids) if ids.count(i) > 1), [])


if __name__ == "__main__":
    unittest.main()

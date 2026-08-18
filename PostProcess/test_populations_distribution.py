#!/usr/bin/env python3
"""Regression test for the DICT-LESS population-distribution crash (v2.3.3).

THE BUG (pre-existing; blocks the dict-less default + any 0-bulge search): a
dict-less complete-search died EXIT=1 at post-analysis in
``populations_distribution.py`` -- ``if max(lower_counts) > 0:`` raised
``ValueError: max() arg is an empty sequence`` (crisprme reports it as
"population distribution plots creation failed ... (mm+bulges: N)"). crisprme
treats the stderr traceback as fatal (``[ -s $logerror ]``), so no
integrated_results.tsv was produced.

Root cause was TWO-fold and both halves are fixed on this branch:

  (A) populations_distribution.py -- guard every max()-on-possibly-empty site so
      an out-of-range ``total`` can never crash (defense in depth): draw_barplot
      returns an empty "no targets" figure when there are no population rows, and
      the two ``max(lower_counts)`` / ``max(lower_barplot_values.values())`` calls
      use ``default=0``.

  (B) submit_job_automated_new_multiple_vcfs.sh -- the pop-dist totals loop is
      bounded by ``mm + bDNA + bRNA`` (the ACTUAL search bulge budget, matching
      process_summaries.py's ``bulge = bDNA + bRNA`` that sizes the file groups),
      NOT ``mm + bMax`` (the reused-index depth N). Verified here as a source check.

This test pins (A): running populations_distribution at a total BEYOND the file's
group count must NOT crash. The PopulationDistribution file for a 0-bulge search
(mm=6, bulge=0) has 7 groups (indices 0..6); asking for total=7,8,9 -- exactly
what the OLD ``mm + bMax`` loop did on an N=3 index -- previously threw. It now
produces an empty plot and returns cleanly.

Loading strategy: matplotlib is NOT installed in the network-free unit-test env
(unit-tests.yml installs only requests + numpy), and populations_distribution.py
imports matplotlib.pyplot / matplotlib.colors at module top and calls
matplotlib.use / plt.rcParams / plt.style at import time. We therefore install a
TARGETED, minimal matplotlib stub in sys.modules (never a catch-all -- a catch-all
would hand a MagicMock to optional-import-fallback targets and break them) and
exec the module source. numpy is present on CI and imported for real.
"""
import os
import sys
import tempfile
import types
import unittest
from unittest.mock import MagicMock


def _install_matplotlib_stub():
    """Register a minimal, real-enough matplotlib so populations_distribution.py
    imports and its plot calls no-op. If a REAL matplotlib is importable, use it and
    do NOT clobber. Otherwise install our package-shaped stub -- even if a prior test
    left a MagicMock under 'matplotlib' in sys.modules (test_bulge_cap /
    test_samplesid_emit stub scientific absents that way): a MagicMock is NOT a
    package, so `from matplotlib.legend_handler import HandlerTuple` would fail. We
    detect that (no real __path__ / __file__) and replace it with the proper stub.
    numpy is imported for real (present on CI)."""
    existing = sys.modules.get("matplotlib")
    if existing is not None:
        # a genuine matplotlib module has a real __path__ pointing at a package dir;
        # a MagicMock / bare ModuleType placeholder does not -> replace it.
        real_path = getattr(existing, "__path__", None)
        if isinstance(real_path, list) and real_path and getattr(existing, "__file__", None):
            return  # real matplotlib already present -> keep it
    else:
        try:
            import matplotlib  # noqa: F401  (present on some dev machines -> use real)
            return
        except Exception:
            pass
    # Drop any incomplete placeholder submodules so our fresh stub wins.
    for _stale in ("matplotlib", "matplotlib.pyplot", "matplotlib.colors",
                   "matplotlib.legend_handler"):
        sys.modules.pop(_stale, None)

    mpl = types.ModuleType("matplotlib")
    # A real `from matplotlib.pyplot import ...` / `from matplotlib.legend_handler
    # import ...` needs matplotlib to look like a PACKAGE (have __path__), else the
    # import machinery raises "'matplotlib' is not a package" before consulting the
    # submodules we register in sys.modules. Empty __path__ is enough.
    mpl.__path__ = []  # type: ignore[attr-defined]
    mpl.use = lambda *a, **k: None

    # matplotlib.pyplot -- every drawing call is a harmless no-op. rcParams must be a
    # real mutable mapping (the module does plt.rcParams["figure.dpi"] = 400 etc.).
    plt = types.ModuleType("matplotlib.pyplot")
    plt.rcParams = {}

    class _Style:
        # module does: `"seaborn-v0_8-poster" in plt.style.available` -> real list.
        available = ["seaborn-poster", "seaborn-v0_8-poster"]

        @staticmethod
        def use(*a, **k):
            return None

    plt.style = _Style()

    def _bar(*a, **k):
        return [MagicMock()]  # bars[i][0] indexing must work downstream

    plt.figure = lambda *a, **k: MagicMock()
    plt.bar = _bar
    plt.legend = lambda *a, **k: None
    plt.title = lambda *a, **k: None
    plt.annotate = lambda *a, **k: None
    plt.xticks = lambda *a, **k: None
    plt.yticks = lambda *a, **k: None
    plt.tight_layout = lambda *a, **k: None
    plt.savefig = lambda *a, **k: None

    mcolors = types.ModuleType("matplotlib.colors")
    mcolors.cnames = {}
    mcolors.to_rgb = lambda c: (0.0, 0.0, 0.0)
    mcolors.TABLEAU_COLORS = {f"tab:{i}": f"#00000{i}" for i in range(9)}

    legend_handler = types.ModuleType("matplotlib.legend_handler")

    class _HandlerTuple:
        def __init__(self, *a, **k):
            pass

    legend_handler.HandlerTuple = _HandlerTuple

    mpl.pyplot = plt
    mpl.colors = mcolors
    mpl.legend_handler = legend_handler
    sys.modules["matplotlib"] = mpl
    sys.modules["matplotlib.pyplot"] = plt
    sys.modules["matplotlib.colors"] = mcolors
    sys.modules["matplotlib.legend_handler"] = legend_handler


def _load_populations_distribution():
    """Import populations_distribution.py as a module object (its __main__ guard
    keeps import side-effect-free apart from the matplotlib setup the stub absorbs)."""
    _install_matplotlib_stub()
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import importlib
    return importlib.import_module("populations_distribution")


pd_mod = _load_populations_distribution()


# A minimal PopulationDistribution file for one guide, 0-bulge search (mm=6):
# each row has mm + 2*bulge + 1 = 7 groups (indices 0..6). Asking for total >= 7
# is out of range -- the crash the fix removes.
_GUIDE = "CCATCGGTGGCCGTTTGCCCNNN"
_POPDIST_7GROUPS = (
    f"-Summary_{_GUIDE}\n"
    "REF\t0,0\t1,0\t2,0\t0,0\t0,0\t0,0\t0,0\n"
    "EAS\t0,0\t0,0\t1,0\t0,0\t0,0\t0,0\t0,0\n"
    "EUR\t0,0\t0,0\t0,0\t0,0\t0,0\t0,0\t0,0\n"
    "AFR\t0,0\t1,0\t3,0\t0,0\t0,0\t0,0\t0,0\n"
)


class TestOutOfRangeTotalDoesNotCrash(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="popdist_")
        self.popfile = os.path.join(self.workdir, "PopulationDistribution_CFD.txt")
        with open(self.popfile, "w") as fh:
            fh.write(_POPDIST_7GROUPS)
        self._cwd = os.getcwd()
        os.chdir(self.workdir)  # draw_barplot writes the PNG into cwd

    def tearDown(self):
        os.chdir(self._cwd)
        import shutil
        shutil.rmtree(self.workdir, ignore_errors=True)

    def _run_total(self, total):
        """Reproduce the full pipeline of create_population_dist_plot for one total,
        exactly as the shell invokes it (argv-driven) but without spawning python."""
        barplot_values, max_value, number_bars = pd_mod.read_population_distribution(
            self.popfile, _GUIDE, total
        )
        lower_barplot_values = pd_mod.compute_lower_bar_values(
            barplot_values, total, number_bars
        )
        final_barplot_values = {
            pop: data[total] for pop, data in barplot_values.items()
        }
        populations = list(final_barplot_values.keys())
        # THE crash site: with an out-of-range total this used to raise
        # ValueError: max() arg is an empty sequence.
        pd_mod.draw_barplot(
            populations,
            final_barplot_values,
            lower_barplot_values,
            max_value,
            total,
            number_bars,
            _GUIDE,
            "CFD",
        )

    def test_out_of_range_totals_do_not_crash(self):
        """totals 7,8,9 exceed the 7-group file (the OLD mm+bMax=6+3 loop) -> must
        NOT raise (empty 'no targets' plot instead)."""
        for total in (7, 8, 9):
            with self.subTest(total=total):
                try:
                    self._run_total(total)
                except Exception as exc:  # noqa: BLE001 -- the regression is ANY crash
                    self.fail(
                        f"populations_distribution crashed at out-of-range total "
                        f"{total}: {type(exc).__name__}: {exc}"
                    )

    def test_in_range_totals_still_work(self):
        """In-range totals (0..6) must still run without error (no regression to the
        happy path)."""
        for total in range(0, 7):
            with self.subTest(total=total):
                try:
                    self._run_total(total)
                except Exception as exc:  # noqa: BLE001
                    self.fail(
                        f"populations_distribution regressed at in-range total "
                        f"{total}: {type(exc).__name__}: {exc}"
                    )

    def test_yrange_empty_lower_values_defaults_to_zero(self):
        """(A) unit: compute_plot_yrange never max()es an empty sequence -- an empty
        lower_barplot_values dict returns the default range, not a ValueError."""
        y_range, nodata = pd_mod.compute_plot_yrange({}, 0)
        self.assertTrue(nodata)
        self.assertIsNotNone(y_range)


class TestShellLoopBoundIsSearchBulgeBudget(unittest.TestCase):
    """(B) source check: the pop-dist totals loop is bounded by the SEARCH bulge
    budget mm + bDNA + bRNA, never mm + bMax (the reused-index depth N)."""

    def setUp(self):
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "submit_job_automated_new_multiple_vcfs.sh")) as fh:
            self.src = fh.read()

    def test_loop_uses_search_bulge_budget(self):
        self.assertIn("pop_dist_total_max=$(expr $mm + $bDNA + $bRNA)", self.src,
                      "pop-dist loop must bound totals by mm + bDNA + bRNA")
        self.assertIn("for total in $(seq 0 $pop_dist_total_max); do", self.src,
                      "pop-dist loop must iterate 0..pop_dist_total_max")

    def test_loop_no_longer_uses_bmax(self):
        self.assertNotIn("seq 0 $(expr $mm + $bMax)", self.src,
                         "pop-dist loop must not use the index-depth bMax bound")


if __name__ == "__main__":
    unittest.main()

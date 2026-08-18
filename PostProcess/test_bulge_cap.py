#!/usr/bin/env python3
"""Regression tests for the DICT-LESS bulge-cap fixes (v2.3.3).

A real user on the web could not select >0 DNA/RNA bulges with the dict-less index
(NRG_3_hg38+hg38_1000G_HGDP): the bulge dropdown was forced to 0, and the shipped
DNA=2/RNA=2 defaults were clobbered to 0. Root cause: the bulge cap is
``min(reference_term, variant_term)``, and a variant HF tarball ships NO reference-only
index (folder ``NRG_*_hg38`` with no ``+``). So on a genuinely-fresh install the
reference term (``index_max_bulges(..., None)``) returned 0, forcing the min -> 0 even
though the variant index supports bmax=2 AND the raw genome ``Genomes/hg38/`` IS shipped
(so a reference index is buildable on demand, exactly as the search shell already does).

These tests pin the three linked fixes, all pure-stdlib / CI-light:

  (1) BULGE-CAP -- the reference term is now BUILDABLE-aware:
      * ``crisprme.py._installed_index_bulge_cap`` (CLI mirror, no web deps): on a
        dict-less layout (variant index + raw genome, NO reference index) returns 2 for a
        variant search, not 0; and stays 0 when neither an index nor a raw genome exists.
        The VARIANT term stays STRICTLY installed-index-based (never buildable dict-less).
      * ``pages_utils.reference_bulge_capacity`` (web): the same asymmetric behavior.

  (3) 0-BULGE INDEL ROUTING -- ``pool_search_indels.py`` routes the 0-bulge case through
      the SHIPPED indexed indel genome (``-index``), not the raw ``Genomes/<ref>+<vcf>_INDELS/``
      (absent on batteries). Verified by running the script against a batteries layout with
      a FAKE ``crispritz.py`` that records its argv and writes the expected per-chrom
      ``fake<chrom>...targets.txt`` -- so post-processing never tail/head/rm/sed's a missing
      file. Also covers the raw-only fallback (non-batteries install).

Defaults (2) are static literals in pages/main_page.index_page (mms=6, dna=2, rna=2,
slider=4) whose SURVIVAL depends on (1); a targeted check that they are the shipped
literals lives in test_dictless_defaults below (source-level, no Dash import).

``crisprme.py`` is loaded via the truncate-and-exec pattern from test_samplesid_emit
(Bio.Seq + absent-scientific-module stubs). ``pages_utils`` is loaded with TARGETED stubs
for its web-only top-level imports (app, dash, pandas) -- never a catch-all.
"""
import os
import shutil
import sys
import tempfile
import types
import unittest


# --------------------------------------------------------------------------- #
# Load crisprme.py DEFINITIONS only (no CLI dispatch), Bio.Seq + scientific
# absents stubbed -- identical to test_samplesid_emit._load_crisprme.
# --------------------------------------------------------------------------- #
def _load_crisprme():
    here = os.path.dirname(os.path.abspath(__file__))
    crisprme_py = os.path.join(os.path.dirname(here), "crisprme.py")
    sys.path.insert(0, here)
    with open(crisprme_py) as fh:
        src = fh.read()
    marker = "\nif len(sys.argv) < 2:"
    cut = src.index(marker)  # KeyError-loud if the dispatch shape ever changes
    defs_only = src[:cut]
    mod = types.ModuleType("crisprme_under_test_bulge")
    mod.__file__ = crisprme_py
    code = compile(defs_only, crisprme_py, "exec")
    if "Bio.Seq" not in sys.modules:
        _bio = types.ModuleType("Bio")
        _bio_seq = types.ModuleType("Bio.Seq")
        _bio_seq.Seq = object
        _bio.Seq = _bio_seq
        sys.modules.setdefault("Bio", _bio)
        sys.modules["Bio.Seq"] = _bio_seq
    from unittest.mock import MagicMock
    for _m in ("pandas", "scipy", "sklearn", "matplotlib", "seaborn",
               "statsmodels", "intervaltree", "CRISTA_score"):
        if _m in sys.modules:
            continue
        try:
            __import__(_m)
        except Exception:
            sys.modules[_m] = MagicMock()
    exec(code, mod.__dict__)
    return mod


# --------------------------------------------------------------------------- #
# Load pages_utils with TARGETED stubs for its web-only top-level imports.
# --------------------------------------------------------------------------- #
def _load_pages_utils():
    from unittest.mock import MagicMock
    # `from app import operators, current_working_directory`
    if "app" not in sys.modules:
        _app = types.ModuleType("app")
        _app.operators = {}
        _app.current_working_directory = os.getcwd() + "/"
        sys.modules["app"] = _app
    # `from dash import html`
    if "dash" not in sys.modules:
        _dash = types.ModuleType("dash")
        _dash.html = MagicMock()
        sys.modules["dash"] = _dash
    # `import pandas as pd` -- present on CI (installed for numpy tests) but stub if absent
    if "pandas" not in sys.modules:
        try:
            __import__("pandas")
        except Exception:
            sys.modules["pandas"] = MagicMock()
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)
    if repo not in sys.path:
        sys.path.insert(0, repo)
    # import the module object so we can set its module-level current_working_directory
    import importlib
    pu = importlib.import_module("pages.pages_utils")
    return pu


crisprme = _load_crisprme()
pages_utils = _load_pages_utils()


# --------------------------------------------------------------------------- #
# Shared dict-less install fixtures.
# --------------------------------------------------------------------------- #
_PAM_VALUE = "20bp-NRG-SpCas9"
_PAM_LINE = "NNNNNNNNNNNNNNNNNNNNNRG 3\n"  # motif NRG, N=3 index depth
_MOTIF = "NRG"
_REF = "hg38"
_VCF = "hg38_1000G_HGDP"            # dropdown value 1000G_HGDP -> enriched hg38_1000G_HGDP
_DATASET = "1000G_HGDP"
_VARIANT_IDX = f"{_MOTIF}_3_{_REF}+{_VCF}"          # N=3 -> bmax 2
_VARIANT_INDELS = _VARIANT_IDX + "_INDELS"


def _mkfile(path, content="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(content)


def _make_dictless_install(root, with_raw_genome=True, with_ref_index=False):
    """The genuinely-fresh dict-less layout: variant index + its _INDELS companion, the
    raw genome Genomes/hg38 (per-chrom .fa), and NO reference-only index -- exactly what a
    variant HF tarball ships. Optionally add the raw genome / a reference index to probe
    the asymmetric cap."""
    os.makedirs(os.path.join(root, "PAMs"), exist_ok=True)
    _mkfile(os.path.join(root, "PAMs", _PAM_VALUE + ".txt"), _PAM_LINE)
    _mkfile(os.path.join(root, "genome_library", _VARIANT_IDX, "TSTgenome.NRG.bin"))
    _mkfile(os.path.join(root, "genome_library", _VARIANT_INDELS,
                         f"{_MOTIF}_3_fakechr22", "TSTgenome.NRG.bin"))
    if with_raw_genome:
        _mkfile(os.path.join(root, "Genomes", _REF, "chr22.fa"), ">chr22\nACGT\n")
    if with_ref_index:  # a reference index that a prior 1-bulge run auto-built (N=2 -> b1)
        _mkfile(os.path.join(root, "genome_library", f"{_MOTIF}_2_{_REF}",
                             "TSTgenome.NRG.bin"))


# --------------------------------------------------------------------------- #
# (1) BULGE-CAP -- CLI mirror _installed_index_bulge_cap
# --------------------------------------------------------------------------- #
class TestCliBulgeCap(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="bulgecap_cli_")
        self._orig_cwd = crisprme.current_working_directory
        crisprme.current_working_directory = self.root + "/"

    def tearDown(self):
        crisprme.current_working_directory = self._orig_cwd
        shutil.rmtree(self.root, ignore_errors=True)

    def test_dictless_variant_search_reaches_bmax2_not_zero(self):
        """THE BUG: dict-less variant search (variant index N=3 + raw genome, NO reference
        index) must reach 2 bulges, not 0. Old code returned min(ref=0, var=2)=0."""
        _make_dictless_install(self.root, with_raw_genome=True, with_ref_index=False)
        cap = crisprme._installed_index_bulge_cap(_MOTIF, _REF, is_variant=True)
        self.assertEqual(cap, 2, "dict-less variant cap must be 2 (min(ref_buildable=2, var=2))")

    def test_reference_only_search_uses_buildable_ceiling(self):
        """A reference-only search with only the raw genome shipped reaches the app
        ceiling (2) because the reference index is buildable on demand."""
        _make_dictless_install(self.root, with_raw_genome=True, with_ref_index=False)
        cap = crisprme._installed_index_bulge_cap(_MOTIF, _REF, is_variant=False)
        self.assertEqual(cap, 2, "reference cap must reflect the buildable-from-raw ceiling")

    def test_variant_term_is_strictly_installed_never_buildable(self):
        """(risk) the VARIANT term must NOT be relaxed to a buildable ceiling: a variant
        index absent (can't be built dict-less) forces the cap to 0 even with a raw
        genome present -- a 2/2 search must NOT try to build a variant index."""
        os.makedirs(os.path.join(self.root, "PAMs"), exist_ok=True)
        _mkfile(os.path.join(self.root, "PAMs", _PAM_VALUE + ".txt"), _PAM_LINE)
        _mkfile(os.path.join(self.root, "Genomes", _REF, "chr22.fa"), ">chr22\nA\n")
        # NO variant index installed
        cap = crisprme._installed_index_bulge_cap(_MOTIF, _REF, is_variant=True)
        self.assertEqual(cap, 0, "variant term must stay installed-only (no buildable relax)")

    def test_no_index_no_raw_genome_is_zero(self):
        """No reference index AND no raw genome -> 0 (only an index-free 0-bulge search)."""
        os.makedirs(os.path.join(self.root, "genome_library"), exist_ok=True)
        cap = crisprme._installed_index_bulge_cap(_MOTIF, _REF, is_variant=False)
        self.assertEqual(cap, 0)

    def test_installed_deeper_reference_index_wins_over_ceiling(self):
        """An already-built DEEPER reference index (N=4 -> b3) is honored above the app
        ceiling: max(installed, buildable)."""
        _mkfile(os.path.join(self.root, "PAMs", _PAM_VALUE + ".txt"), _PAM_LINE)
        _mkfile(os.path.join(self.root, "genome_library", f"{_MOTIF}_4_{_REF}",
                             "TSTgenome.NRG.bin"))
        _mkfile(os.path.join(self.root, "Genomes", _REF, "chr22.fa"), ">c\nA\n")
        cap = crisprme._installed_index_bulge_cap(_MOTIF, _REF, is_variant=False)
        self.assertEqual(cap, 3, "an installed deeper reference index must win")


# --------------------------------------------------------------------------- #
# (1) BULGE-CAP -- web pages_utils.reference_bulge_capacity
# --------------------------------------------------------------------------- #
class TestWebReferenceBulgeCapacity(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="bulgecap_web_")
        self._orig_cwd = pages_utils.current_working_directory
        pages_utils.current_working_directory = self.root

    def tearDown(self):
        pages_utils.current_working_directory = self._orig_cwd
        shutil.rmtree(self.root, ignore_errors=True)

    def test_raw_genome_present_no_ref_index_gives_buildable_ceiling(self):
        """THE BUG (web): raw genome shipped, no reference index -> capacity 2, not 0.
        (index_max_bulges(...,None)=0 previously forced the dropdown/defaults to 0.)"""
        _make_dictless_install(self.root, with_raw_genome=True, with_ref_index=False)
        self.assertEqual(pages_utils.index_max_bulges(_REF, _PAM_VALUE, None), 0,
                         "sanity: no reference index installed -> installed depth 0")
        self.assertEqual(pages_utils.index_max_bulges(_REF, _PAM_VALUE, _DATASET), 2,
                         "sanity: variant index N=3 -> 2 usable bulges")
        self.assertEqual(pages_utils.reference_bulge_capacity(_REF, _PAM_VALUE), 2,
                         "reference capacity must be the buildable ceiling (2)")

    def test_no_raw_genome_no_index_is_zero(self):
        """No raw genome and no reference index -> 0 (index-free only)."""
        _make_dictless_install(self.root, with_raw_genome=False, with_ref_index=False)
        self.assertEqual(pages_utils.reference_bulge_capacity(_REF, _PAM_VALUE), 0)

    def test_dictless_effective_cap_is_two(self):
        """End-to-end cap the dropdown uses: min(reference_bulge_capacity, variant term)
        = min(2, 2) = 2 on the dict-less install -> dropdown offers 0,1,2 and DNA/RNA=2
        survive (the whole point of the fix)."""
        _make_dictless_install(self.root, with_raw_genome=True, with_ref_index=False)
        ref = pages_utils.reference_bulge_capacity(_REF, _PAM_VALUE)
        var = pages_utils.index_max_bulges(_REF, _PAM_VALUE, _DATASET)
        self.assertEqual(min(ref, var), 2)


# --------------------------------------------------------------------------- #
# (2) DEFAULTS -- static literals in main_page.index_page survive (source check)
# --------------------------------------------------------------------------- #
class TestDictlessDefaults(unittest.TestCase):
    """The intended form defaults are mismatches=6, DNA bulge=2, RNA bulge=2,
    max-total-edits=4. These are static literals in pages/main_page.index_page; assert
    they are present as-shipped (no Dash import needed -- a targeted source scan)."""

    def setUp(self):
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(os.path.dirname(here), "pages", "main_page.py")) as fh:
            self.src = fh.read()

    def _dropdown_value(self, comp_id):
        """The `value=<int>` literal of the dcc.Dropdown whose id is `comp_id`."""
        import re
        # find `id="<comp_id>"` then walk back to the nearest `value=<int>` in the block
        m = re.search(rf'value=(\d+),[^)]*?id="{comp_id}"', self.src, re.S)
        self.assertIsNotNone(m, f"could not locate value=... id={comp_id} block")
        return int(m.group(1))

    def test_mismatch_default_is_6(self):
        self.assertEqual(self._dropdown_value("mms"), 6)

    def test_dna_bulge_default_is_2(self):
        self.assertEqual(self._dropdown_value("dna"), 2)

    def test_rna_bulge_default_is_2(self):
        self.assertEqual(self._dropdown_value("rna"), 2)

    def test_max_edits_slider_default_is_4(self):
        import re
        m = re.search(r'id="max-edits-slider".*?value=(\d+)', self.src, re.S)
        self.assertIsNotNone(m, "could not locate the max-edits-slider value")
        self.assertEqual(int(m.group(1)), 4)

    def test_cli_max_total_edits_default_is_4(self):
        """CLI parity: complete-search --max-total-edits default is 4."""
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(os.path.dirname(here), "crisprme.py")) as fh:
            csrc = fh.read()
        self.assertRegex(csrc, r"max_total_edits\s*=\s*4",
                         "crisprme.py complete-search default max_total_edits must be 4")


# --------------------------------------------------------------------------- #
# (3) 0-BULGE INDEL ROUTING -- exercise pool_search_indels.search_indels directly.
#
# The script parses sys.argv at module top, then defines search_indels, then runs a
# multiprocessing Pool. We exec ONLY the part up to the Pool (identical truncate-and-exec
# pattern as _load_crisprme) with a controlled argv and os.system monkeypatched to record
# the command -- so the routing decision is tested deterministically, with no Pool (which
# under macOS `spawn` would re-import this test module) and no real crispritz.
# --------------------------------------------------------------------------- #
def _load_pool_search_indels(argv, os_system_sink):
    """Exec pool_search_indels.py's DEFS + globals (argv-driven) but NOT the Pool block.
    Returns the module namespace with search_indels bound; os.system is patched to append
    each command to `os_system_sink` (so the routed crispritz command is captured)."""
    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, "pool_search_indels.py")
    with open(script) as fh:
        src = fh.read()
    marker = "\nchrs = _dataset_chroms("
    cut = src.index(marker)  # loud if the script's shape changes
    defs_only = src[:cut]
    ns = {"__name__": "pool_search_indels_under_test", "__file__": script}
    real_argv = sys.argv
    sys.argv = list(argv)
    try:
        exec(compile(defs_only, script, "exec"), ns)
    finally:
        sys.argv = real_argv

    def _fake_system(cmd):
        os_system_sink.append(cmd)
        return 0

    ns["os"].system = _fake_system
    return ns


class TestIndelZeroBulgeRouting(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="indel_route_")
        self.out = os.path.join(self.root, "out")
        os.makedirs(self.out, exist_ok=True)
        self.guide = os.path.join(self.root, "guides.txt")
        _mkfile(self.guide, "N" * 20 + "\n")
        self.pam = os.path.join(self.root, "pam.txt")
        _mkfile(self.pam, _PAM_LINE)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def _argv(self, bDNA, bRNA):
        #  1 ref_folder 2 vcf_dir 3 vcf_name 4 guide 5 pam 6 bMax 7 mm 8 bDNA 9 bRNA
        #  10 output_folder 11 true_pam 12 cwd 13 threads 14 max_edits
        return [
            "pool_search_indels.py",
            os.path.join(self.root, "Genomes", _REF),   # ref_folder -> ref_name hg38
            "_", _VCF, self.guide, self.pam,
            "3",                                          # bMax == idx_n
            "6", str(bDNA), str(bRNA),
            self.out,
            _MOTIF,                                       # true_pam == idx_pam
            self.root + "/",                              # current_working_directory
            "1", "4",
        ]

    def _search(self, bDNA, bRNA, chrom="chr22"):
        sink = []
        ns = _load_pool_search_indels(self._argv(bDNA, bRNA), sink)
        ns["search_indels"](chrom)
        self.assertEqual(len(sink), 1, f"expected exactly one crispritz call, got {sink}")
        return sink[0]

    def test_zero_bulge_uses_indexed_genome_when_only_indexed_shipped(self):
        """THE FIX: batteries/dict-less layout ships ONLY the indexed _INDELS genome. A
        0-bulge indel search must route through it with -index (not the absent raw
        Genomes/<ref>+<vcf>_INDELS/), so the per-chrom targets file is produced and the
        combined merge (tail/head/rm/sed) never hits a missing file."""
        _mkfile(os.path.join(self.root, "genome_library", _VARIANT_INDELS,
                             f"{_MOTIF}_3_fakechr22", "TSTgenome.NRG.bin"))
        cmd = self._search("0", "0")
        self.assertIn("genome_library", cmd,
                      f"0-bulge indel search must use the INDEXED genome: {cmd}")
        self.assertIn(_VARIANT_INDELS, cmd)
        self.assertIn(f"{_MOTIF}_3_fakechr22", cmd)
        self.assertIn("-index", cmd, "indexed search must pass -index for 0 bulges")
        self.assertNotIn("/Genomes/", cmd,
                         "0-bulge search must NOT touch the raw Genomes/<ref>+<vcf>_INDELS/")

    def test_bulge_search_still_uses_indexed_genome(self):
        """A >0-bulge indel search likewise uses the indexed genome with -index (unchanged
        behavior, now guarded by isdir checks)."""
        _mkfile(os.path.join(self.root, "genome_library", _VARIANT_INDELS,
                             f"{_MOTIF}_3_fakechr22", "TSTgenome.NRG.bin"))
        cmd = self._search("2", "2")
        self.assertIn("genome_library", cmd)
        self.assertIn(_VARIANT_INDELS, cmd)
        self.assertIn("-index", cmd)
        self.assertIn("-bDNA 2", cmd)
        self.assertIn("-bRNA 2", cmd)

    def test_raw_only_fallback_used_when_no_indexed_genome(self):
        """Fallback: a non-batteries install with ONLY the raw fake-indel genome (no
        indexed one) runs a 0-bulge search against the raw genome (brute-force, no
        -index)."""
        _mkfile(os.path.join(self.root, "Genomes", f"{_REF}+{_VCF}_INDELS",
                             f"fake_{_VCF}_chr22", "fakechr22.fa"), ">f\nA\n")
        cmd = self._search("0", "0")
        self.assertIn(f"{_REF}+{_VCF}_INDELS", cmd)
        self.assertIn("/Genomes/", cmd, "raw fallback must search under Genomes/")
        self.assertNotIn("-index", cmd, "raw fallback must NOT pass -index")

    def test_indexed_preferred_even_when_raw_also_present(self):
        """When BOTH the indexed and the raw indel genome exist, the indexed one is
        preferred (it always serves any bulge depth)."""
        _mkfile(os.path.join(self.root, "genome_library", _VARIANT_INDELS,
                             f"{_MOTIF}_3_fakechr22", "TSTgenome.NRG.bin"))
        _mkfile(os.path.join(self.root, "Genomes", f"{_REF}+{_VCF}_INDELS",
                             f"fake_{_VCF}_chr22", "fakechr22.fa"), ">f\nA\n")
        cmd = self._search("0", "0")
        self.assertIn("genome_library", cmd)
        self.assertIn("-index", cmd)

    def test_missing_both_genomes_still_targets_indexed_path(self):
        """(safety) if neither genome dir exists (should not happen -- the shell guarantees
        the indexed one or errors out first), the routing still targets the indexed -index
        path so the failure is a clear crispritz error, not a raw-genome brute-force on a
        path that also does not exist."""
        cmd = self._search("0", "0")
        self.assertIn("genome_library", cmd)
        self.assertIn("-index", cmd)


if __name__ == "__main__":
    unittest.main()

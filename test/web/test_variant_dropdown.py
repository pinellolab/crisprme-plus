#!/usr/bin/env python3
"""Regression test: the batteries-included variant INDEX must be selectable, and the
combined variant search must have the inputs it needs.

Guards the alpha bug where, after the batteries download (index only, no raw VCFs/),
the variant dropdown showed only "Reference only" -- get_variant_dataset_options read
only VCFs/ + a built-in marker that did not match the combined enriched index name, so
the shipped NGG_3_hg38+hg38_1000G_HGDP index was never offered. Also covers the VCF-folder
resolver and the combined-samplesID synthesis (the combined search dies at the sample-ID
step without the combined list, which the batteries download does not ship).

Standalone (no pytest); builds its own mock install and exits non-zero on failure, so it
can run directly in CI inside the full image env:

    docker run --rm -v "$PWD/test/web:/t" crisprme:web python /t/test_variant_dropdown.py
"""
import os
import sys
import tempfile
import shutil

sys.path.insert(0, "/opt/conda/opt/crisprme")

_fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL "), msg)
    if not cond:
        _fails.append(msg)


# Mock the exact failing layout: the combined index is installed, VCFs/ is EMPTY (the
# batteries-included download fetches the index, not the raw VCFs), and only the
# per-dataset sample lists are present.
d = tempfile.mkdtemp(prefix="crisprme_batt_")
for sub in (
    "VCFs",
    "samplesIDs",
    "Genomes/hg38",
    "genome_library/NGG_3_hg38",
    "genome_library/NGG_3_hg38+hg38_1000G_HGDP",
    "genome_library/NGG_3_hg38+hg38_1000G_HGDP_INDELS",
    "Annotations",
):
    os.makedirs(os.path.join(d, sub), exist_ok=True)
with open(os.path.join(d, "samplesIDs/hg38_1000G.samplesID.txt"), "w") as fh:
    fh.write("#SAMPLE\tPOP\tSUPERPOP\tSEX\nNA1\tGBR\tEUR\tM\n")
with open(os.path.join(d, "samplesIDs/hg38_HGDP.samplesID.txt"), "w") as fh:
    fh.write("#SAMPLE\tPOP\tSUPERPOP\tSEX\nHG1\tFrench\tEUR\tM\n")
# the batteries download ships the annotation bundle bgzipped (.bed.gz)
with open(os.path.join(d, "Annotations/dhs+encode+gencode.hg38.bed.gz"), "w") as fh:
    fh.write("")

# 1) discovery: the installed combined index must be offered as a variant option
from pages import pages_utils  # noqa: E402

pages_utils.current_working_directory = d
vals = [o["value"] for o in pages_utils.get_variant_dataset_options("hg38")]
check("1000G_HGDP" in vals,
      f"variant dropdown surfaces the installed combined index (got {vals})")

# 1b) annotation: the shipped bundle is bgzipped (.bed.gz) -- the option must still
# surface (checking only ".bed" left every default search with "No annotation"), and
# the callback defaults the value to it so annotations are ON by default.
avals = [o["value"] for o in pages_utils.get_annotation_options("hg38")]
check("EN" in avals,
      f"annotation option surfaces for the bgzipped bundle -> ON by default (got {avals})")

# 2) resolver + 3) combined-samplesID synthesis (main_page helpers)
from pages import main_page  # noqa: E402

main_page.current_working_directory = d
check(main_page._resolve_vcf_folder("hg38", "1000G_HGDP") == "hg38_1000G_HGDP",
      "resolver returns the enriched folder name for the combined dataset")
# genome dir name differing from the VCF folder prefix must still resolve
os.makedirs(os.path.join(d, "VCFs/hg38_1000G"), exist_ok=True)
check(main_page._resolve_vcf_folder("hg38_chr22", "1000G") == "hg38_1000G",
      "resolver is robust to a genome dir whose name differs from the VCF folder prefix")
main_page._ensure_samplesid("hg38", "hg38_1000G_HGDP")
combined = os.path.join(d, "samplesIDs/hg38_1000G_HGDP.samplesID.txt")
check(os.path.isfile(combined),
      "combined samplesID is synthesized from the per-dataset lists")
if os.path.isfile(combined):
    n = len([ln for ln in open(combined) if ln.strip() and not ln.startswith("#")])
    check(n == 2, f"combined samplesID unions both datasets (got {n} rows)")

shutil.rmtree(d, ignore_errors=True)
print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
sys.stdout.flush()
sys.stderr.flush()
# force-exit: importing the Dash app can leave non-daemon helper threads that would
# otherwise stall interpreter shutdown; this keeps the CI step's exit code clean.
# (os._exit skips buffer flushing, hence the explicit flush above.)
os._exit(1 if _fails else 0)

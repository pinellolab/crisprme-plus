#!/usr/bin/env python3
"""Annotation-manager foundation tests: BED validation, the per-genome enabled-set
manifest (incl. backward-compat default), and active-annotation assembly (bundle-only
fast path, empty case, multi-bed merge). Standalone; exits non-zero on failure. Runs
inside the image (the merge path needs sort-bed/bgzip)."""
import os
import sys
import tempfile
import shutil
import subprocess

sys.path.insert(0, "/opt/conda/opt/crisprme")
from pages import pages_utils as pu  # noqa: E402

_fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL "), msg)
    if not cond:
        _fails.append(msg)


d = tempfile.mkdtemp(prefix="annmgr_")
ann = os.path.join(d, "Annotations")
os.makedirs(ann)
pu.current_working_directory = d + "/"

BUILTIN = "dhs+encode+gencode.hg38.bed.gz"


def _write_gz(gz_path, text):
    plain = gz_path[:-3]
    with open(plain, "w") as fh:
        fh.write(text)
    subprocess.run(["bgzip", "-f", plain], check=True)  # -> gz_path


_write_gz(os.path.join(ann, BUILTIN), "chr1\t100\t200\tccre\nchr2\t300\t400\tdhs\n")
with open(os.path.join(ann, "gencode.protein_coding.bed"), "w") as fh:
    fh.write("chr1\t50\t60\tgene\n")
with open(os.path.join(ann, "mytrack.hg38.bed"), "w") as fh:
    fh.write("chr1\t10\t20\tfoo\nchr3\t5\t9\tbar\n")

# 1) validate_annotation_bed
check(pu.validate_annotation_bed(os.path.join(ann, "mytrack.hg38.bed")) is None,
      "valid .bed passes")
check(pu.validate_annotation_bed(os.path.join(ann, BUILTIN)) is None,
      "valid .bed.gz passes")
bad = os.path.join(ann, "bad.hg38.bed")
open(bad, "w").write("chr1\tNOTNUM\t20\tx\n")
check(pu.validate_annotation_bed(bad) is not None, "non-integer start rejected")
sp = os.path.join(ann, "space.hg38.bed")
open(sp, "w").write("chr1\t1\t2\tlabel with space\n")
check(pu.validate_annotation_bed(sp) is not None, "whitespace label rejected")
two = os.path.join(ann, "two.hg38.bed")
open(two, "w").write("chr1\t1\n")
check(pu.validate_annotation_bed(two) is not None, "<4 columns rejected")
rev = os.path.join(ann, "rev.hg38.bed")
open(rev, "w").write("chr1\t900\t100\tx\n")
check(pu.validate_annotation_bed(rev) is not None, "start>end rejected")
# clean the invalid scratch beds so they don't pollute the later listing check
for f in (bad, sp, two, rev):
    os.remove(f)

# 2) enabled manifest + backward-compat default
check(pu.read_enabled_annotations("hg38") == [BUILTIN],
      "no manifest -> built-in enabled for hg38 (annotations ON by default)")
check(pu.read_enabled_annotations("susScr11") == [],
      "no manifest -> none for a non-hg38 genome")
pu.write_enabled_annotations("hg38", [BUILTIN, "mytrack.hg38.bed"])
check(set(pu.read_enabled_annotations("hg38")) == {BUILTIN, "mytrack.hg38.bed"},
      "write/read round-trip")

# 3) build_active_annotation
pu.write_enabled_annotations("hg38", [BUILTIN])
a, g = pu.build_active_annotation("hg38")
check(a == BUILTIN and g == "gencode.protein_coding.bed",
      f"bundle-only -> fast path (got {a}, {g})")
pu.write_enabled_annotations("hg38", [])
check(pu.build_active_annotation("hg38") == ("vuoto.txt", "vuoto.txt"),
      "empty enabled -> vuoto.txt")
pu.write_enabled_annotations("hg38", [BUILTIN, "mytrack.hg38.bed"])
a, g = pu.build_active_annotation("hg38")
check(a.startswith(".active.hg38") and a.endswith(".bed.gz"),
      f"multi-enabled -> merged active annotation (got {a})")
check(g == "gencode.protein_coding.bed", "multi with built-in keeps gencode")
n = int(subprocess.run(["bash", "-lc", f"zcat {os.path.join(ann, a)} | wc -l"],
                       capture_output=True, text=True).stdout.strip())
check(n == 4, f"merged annotation has all intervals (2 built-in + 2 custom = 4, got {n})")
vals = [o["value"] for o in pu.get_custom_annotations()]
check(a not in vals and BUILTIN not in vals and "mytrack.hg38.bed" in vals,
      f"active/built-in excluded from the user list, custom listed (got {vals})")

shutil.rmtree(d, ignore_errors=True)
print(f"\n{'FAILED' if _fails else 'OK'}: {len(_fails)} failure(s)")
sys.stdout.flush()
os._exit(1 if _fails else 0)

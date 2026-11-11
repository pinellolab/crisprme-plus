#!/usr/bin/env python3
"""The panel recipes the rebuild of 2.10.0 calls (seq_script/merge_panels/hgdp1kgp_build.sh and
hprc_build.sh, from seq_script/rebuild/index_set_2.10.0.tsv) stop at a step that fails, so that
the rebuild's script, which stops at the first command that fails, stops there.

They printed BUILD_EXIT and went on: after a failed build-index-only the HGDP + 1000 Genomes
recipe exited 0 wherever an earlier index of the panel was in place (its last command lists
them), so the rebuild went on for hours, published the reference indexes and the panel's first
index, and stopped only at the missing index's publish-index (the stage-2b review of lane D,
D-R3-4). The steps before were silent too: a pipeline without pipefail took the status of its
last stage, so a VCF cut short by a stage that failed got its .tbi and was built on (and a run
again skipped it); a chromosome of the genome folder whose source was missing was skipped; and
a samples file of no sample was built on. Each now stops the recipe with a non-zero status. A
chromosome the genome folder does not hold (a test genome of chr22) is still not built, and a
build that succeeds still exits 0 (those two pass on the recipes before too).

Hermetic: the recipes run on stubs (apptainer, bcftools, tabix, zcat, crisprme.py) in a
temporary directory; nothing of the image runs.
"""

import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RECIPES = os.path.join(ROOT, "seq_script", "merge_panels")

# apptainer exec [options] <image> bash -c "<command>": the command run here, on the stub tools,
# with the image's /tmp a folder of the test's and its /opt/conda/bin the stubs' folder
_APPTAINER = r'''#!/usr/bin/env bash
cmd=
while [ $# -gt 0 ]; do
  if [ "$1" = bash ] && [ "${2:-}" = -c ]; then cmd=$3; break; fi
  shift
done
[ -n "$cmd" ] || { echo "stub apptainer: no bash -c" >&2; exit 64; }
cmd=${cmd//\/tmp\/sort_/$STUB_TMP/sort_}
cmd=${cmd//TMPDIR=\/tmp;/TMPDIR=$STUB_TMP;}
cmd=${cmd//\/opt\/conda\/bin/$STUB_BIN}
exec bash -c "$cmd"
'''
# a stage named in STUB_FAIL writes part of its output, then fails (as one killed between two
# records does: the stages after it read a stream that ends there)
_BCFTOOLS = r'''#!/usr/bin/env bash
sub=$1; shift
case "$sub" in
  query)
    [ -e "${2:-}" ] || { echo "stub bcftools: no input ${2:-}" >&2; exit 1; }
    cat "$2" ;;
  view)
    if [ "${STUB_FAIL:-}" = view ]; then head -c 3 "$1"; echo "stub bcftools: view failed" >&2; exit 1; fi
    [ -s "$1" ] || { echo "stub bcftools: no input $1" >&2; exit 1; }
    cat "$1" ;;
  norm)
    fa=; inp=
    while [ $# -gt 0 ]; do
      case "$1" in -f) fa=$2; shift 2 ;; -m) shift 2 ;; -O*) shift ;; *) inp=$1; shift ;; esac
    done
    [ -e "$fa" ] || { echo "stub bcftools: no reference $fa" >&2; exit 1; }
    if [ "${STUB_FAIL:-}" = norm ]; then head -c 3 "${inp:-/dev/stdin}"; echo "stub bcftools: norm failed" >&2; exit 1; fi
    if [ -n "$inp" ]; then cat "$inp"; else cat; fi ;;
  sort)
    out=
    while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift 2 ;; *) shift ;; esac; done
    cat > "$out"
    # (as bcftools does, nothing to read is an error; a stream cut short between two records
    # is not)
    [ -s "$out" ] || { rm -f "$out"; echo "stub bcftools: Failed to read from standard input" >&2; exit 1; } ;;
  *) cat ;;
esac
'''
_TABIX = r'''#!/usr/bin/env bash
for a in "$@"; do last=$a; done
touch "$last.tbi"
'''
_ZCAT = r'''#!/usr/bin/env bash
cat "$@"
'''
_CRISPRME = r'''#!/usr/bin/env bash
echo "$*" >> "$STUB_LOG"
path=
while [ $# -gt 0 ]; do [ "$1" = --path ] && path=$2; shift; done
if [ "${STUB_BUILD_RC:-0}" != 0 ]; then echo "stub: build-index-only failed" >&2; exit "$STUB_BUILD_RC"; fi
mkdir -p "$path/genome_library/$STUB_INDEX" "$path/genome_library/${STUB_INDEX}_INDELS"
'''


class _RecipeCases:
    """The cases, the same for both recipes (a mixin of each recipe's TestCase)."""

    RECIPE = None         # its file name in seq_script/merge_panels/
    DATASET = None        # the panel it builds
    INDEX = None          # the index its build writes, for the PAM of the test
    FAILING_STAGE = None  # a stage of its VCF pipeline that fails in one case

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.t = os.path.realpath(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self.bin = os.path.join(self.t, "bin")
        os.mkdir(self.bin)
        for name, text in (("apptainer", _APPTAINER), ("bcftools", _BCFTOOLS),
                           ("tabix", _TABIX), ("zcat", _ZCAT), ("crisprme.py", _CRISPRME)):
            path = os.path.join(self.bin, name)
            with open(path, "w") as fh:
                fh.write(text)
            os.chmod(path, 0o755)
        for d in ("B", "G", "SRC", "stubtmp"):
            os.mkdir(os.path.join(self.t, d))
        self.log = os.path.join(self.t, "builds.log")
        open(self.log, "w").close()
        self.genome("chr22")
        self.env = self.make_inputs()

    def genome(self, *chroms):
        for c in chroms:
            with open(os.path.join(self.t, "G", f"{c}.fa"), "w") as fh:
                fh.write(f">{c}\nACGT\n")

    def run_recipe(self, **env):
        base = {"PATH": os.pathsep.join([self.bin, "/usr/bin", "/bin"]), "HOME": self.t,
                "LC_ALL": "C", "STUB_BIN": self.bin, "STUB_TMP": os.path.join(self.t, "stubtmp"),
                "STUB_LOG": self.log, "STUB_INDEX": self.INDEX, "SIF": "image.sif",
                "B": os.path.join(self.t, "B"), "GENOME": os.path.join(self.t, "G"),
                "PAM": os.path.join(self.t, "pam.txt"), "THREADS": "1"}
        base.update(self.env)
        base.update(env)
        return subprocess.run([shutil.which("bash"), os.path.join(RECIPES, self.RECIPE)],
                              env=base, capture_output=True, text=True, timeout=120)

    def builds(self):
        with open(self.log) as fh:
            return [line for line in fh if line.startswith("build-index-only")]

    def tbi(self, chrom):
        return os.path.join(self.t, "B", "VCFs", self.DATASET, self.vcf_name(chrom) + ".tbi")

    def test_a_build_that_succeeds_exits_0(self):
        # (on a genome folder of chr22 alone: the recipe's other chromosomes are not built)
        got = self.run_recipe()
        self.assertEqual(got.returncode, 0, got.stdout + got.stderr)
        self.assertEqual(len(self.builds()), 1)
        self.assertTrue(os.path.exists(self.tbi("chr22")))
        self.assertFalse(os.path.exists(self.tbi("chr21")))
        self.assertIn(f"index: {self.t}/B/genome_library/{self.INDEX} (+ _INDELS)", got.stdout)

    def test_a_failed_build_stops_it(self):
        # the panel's earlier index in its working directory, as the rebuild leaves it
        for folder in (self.INDEX.replace("NRR", "NGG"),
                       self.INDEX.replace("NRR", "NGG") + "_INDELS"):
            os.makedirs(os.path.join(self.t, "B", "genome_library", folder))
        got = self.run_recipe(STUB_BUILD_RC="3")
        self.assertEqual(got.returncode, 3, got.stdout + got.stderr)
        self.assertIn("build-index-only failed (exit 3)", got.stderr)
        self.assertNotIn("index: ", got.stdout)

    def test_a_failed_stage_of_a_chromosome_s_vcf_stops_it_before_the_build(self):
        got = self.run_recipe(STUB_FAIL=self.FAILING_STAGE)
        self.assertNotEqual(got.returncode, 0, got.stdout)
        self.assertIn("[error] [chr22]", got.stderr)
        self.assertEqual(self.builds(), [], "it built on a VCF cut short")
        self.assertFalse(os.path.exists(self.tbi("chr22")),
                         "a .tbi says the VCF is whole, and a run again skips it")

    def test_a_chromosome_of_the_genome_without_its_source_stops_it(self):
        self.genome("chr21")
        got = self.run_recipe()
        self.assertNotEqual(got.returncode, 0, got.stdout)
        self.assertIn("[error] [chr21]", got.stderr)
        self.assertEqual(self.builds(), [])


class HgdpRecipe(_RecipeCases, unittest.TestCase):
    RECIPE = "hgdp1kgp_build.sh"
    DATASET = "hg38_HGDP1kGP"
    INDEX = "20-NRR_3_hg38+hg38_HGDP1kGP"
    FAILING_STAGE = "view"

    def vcf_name(self, chrom):
        return f"hgdp1kgp.{chrom}.norm.vcf.gz"

    def make_inputs(self):
        """chr22's re-phased BCF (to the stub, its samples, one a line) and gnomAD's metadata
        (the sample, then is_female in column 80, the population and region in 176 and
        177)."""
        self.bcf = os.path.join(self.t, "SRC", "hgdp1kg.chr22.rephased.bcf")
        with open(self.bcf, "w") as fh:
            fh.write("S1\nS2\n")
        meta = os.path.join(self.t, "meta.tsv")
        with open(meta, "w") as fh:
            for row in [["s"] + ["x"] * 176] + [[s] + ["x"] * 78 + ["false"] + ["x"] * 95
                                                 + ["POP", "SUP"] for s in ("S1", "S2")]:
                fh.write("\t".join(row) + "\n")
        return {"RP": os.path.join(self.t, "SRC"), "META": meta}

    def test_the_samples_of_the_bcf_unread_stop_it(self):
        os.remove(self.bcf)
        got = self.run_recipe()
        self.assertNotEqual(got.returncode, 0, got.stdout)
        self.assertIn("cannot read the samples", got.stderr)
        self.assertEqual(self.builds(), [])

    def test_a_bcf_of_no_sample_stops_it(self):
        open(self.bcf, "w").close()
        got = self.run_recipe()
        self.assertNotEqual(got.returncode, 0, got.stdout)
        self.assertIn("no sample in", got.stderr)
        self.assertEqual(self.builds(), [])


class HprcRecipe(_RecipeCases, unittest.TestCase):
    RECIPE = "hprc_build.sh"
    DATASET = "hg38_HPRC"
    INDEX = "20-NRR_3_hg38+hg38_HPRC"
    FAILING_STAGE = "norm"

    def vcf_name(self, chrom):
        return f"hprc.{chrom}.norm.vcf.gz"

    def make_inputs(self):
        """chr22's VCF (plain text: the stub zcat reads it) and the 1000 Genomes samples."""
        self.vcf = os.path.join(self.t, "SRC", "ALL.chr22.HPRC.vcfwave.vcf.gz")
        self.write_vcf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\n")
        p1kg = os.path.join(self.t, "p1kg.txt")
        with open(p1kg, "w") as fh:
            fh.write("#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\nHG1\tGBR\tEUR\tmale\n")
        return {"SRC": os.path.join(self.t, "SRC"), "P1KG": p1kg}

    def write_vcf(self, header):
        with open(self.vcf, "w") as fh:
            fh.write(f"##fileformat=VCFv4.2\n{header}chr22\t1\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|1\n")

    def test_the_samples_of_the_vcf_unread_stop_it(self):
        self.write_vcf("")  # no #CHROM line: no sample to read
        got = self.run_recipe()
        self.assertNotEqual(got.returncode, 0, got.stdout)
        self.assertIn("cannot read the samples", got.stderr)
        self.assertEqual(self.builds(), [])


if __name__ == "__main__":
    unittest.main()

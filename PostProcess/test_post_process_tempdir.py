"""Regression test: post_process.sh must NOT write annotation temp BEDs into the
install dir (PostProcess/). In a no-annotation run the annotation arg is the
empty placeholder PostProcess/vuoto.txt, which is read-only on an apptainer SIF;
the old script wrote vuoto.txt.tmp.bed/.tmp.sorted.bed next to it and failed.
This test makes PostProcess/ read-only, runs the annotation-prep prefix of
post_process.sh with $2=vuoto.txt, and asserts the temps land under the (writable)
output dir instead. Requires bedops (sort-bed / closest-features) on PATH; skipped
if absent.

Run with:
    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_post_process_tempdir.py' -v
"""
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))


class PostProcessTempDir(unittest.TestCase):
    @unittest.skipUnless(shutil.which("sort-bed") and shutil.which("closest-features"),
                         "bedops (sort-bed/closest-features) not on PATH")
    def test_no_annotation_temps_go_to_output_dir(self):
        vuoto = os.path.join(PP, "vuoto.txt")
        self.assertTrue(os.path.isfile(vuoto))
        out = tempfile.mkdtemp()
        try:
            # minimal bestMerge-like file: awk in post_process.sh reads cols $5,$7,$3
            best = os.path.join(out, os.path.basename(out) + ".bestMerge.txt")
            with open(best, "w") as fh:
                fh.write("#h\n")
                # cols 1..7 with col5=chr1 col7=1000 col3=GGGG
                fh.write("\t".join(["x","x","GGGG","x","chr1","x","1000"]) + "\n")
            # Run only the annotation-prep prefix of post_process.sh (through the
            # sort-bed of the annotation) with $2 = vuoto.txt, PostProcess/ read-only.
            mode = os.stat(PP).st_mode
            os.chmod(PP, mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH)
            try:
                script = (
                    'set -e; dir=$(dirname "$1"); annot_tmp="$dir/$(basename "$2")"; '
                    'if [ "$(basename "$2")" != "vuoto.txt" ]; then gunzip -k -c "$2" > "$annot_tmp.tmp.bed"; '
                    'else cp "$2" "$annot_tmp.tmp.bed"; fi; '
                    'sort-bed "$annot_tmp.tmp.bed" > "$annot_tmp.tmp.sorted.bed"; rm "$annot_tmp.tmp.bed"'
                )
                r = subprocess.run(["bash", "-c", script, "_", best, vuoto],
                                   capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)
            finally:
                os.chmod(PP, mode)
            # temp landed under the output dir, NOT in PostProcess/
            self.assertTrue(os.path.isfile(os.path.join(out, "vuoto.txt.tmp.sorted.bed")))
            self.assertFalse(os.path.exists(os.path.join(PP, "vuoto.txt.tmp.bed")))
            self.assertFalse(os.path.exists(os.path.join(PP, "vuoto.txt.tmp.sorted.bed")))
        finally:
            shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

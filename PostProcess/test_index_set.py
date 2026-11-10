#!/usr/bin/env python3
"""The index set of CRISPRme 2.10.0 (the rebuild manifest, seq_script/rebuild/index_set_2.10.0.tsv)
holds against the code that builds, names, publishes and downloads its archives
(seq_script/rebuild/check_index_set.py, which the D7 readiness audit runs before any build),
and the check finds what it is for: a name build-index-only would not give, an N that is not
the build's bulges + 1, a panel without its bundle, a genotyped panel without its store, an
archive of 2.8.1 nothing replaces (NNN_3's pointer to the general index), an index the web's
2 + 2 bulges do not fit. And the rebuild its commands make runs: the script --commands prints,
run dry (check_index_set.dry_run: every program stubbed, in temporary working directories), runs
to its end from the repository's root, and the check finds each way a command fails there.

Hermetic: reads the repository's files, and runs the rebuild's script dry in a temporary
directory (bash, and this Python for the stubs).
"""

import copy
import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_spec = importlib.util.spec_from_file_location(
    "check_index_set", os.path.join(ROOT, "seq_script", "rebuild", "check_index_set.py"))
cis = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cis)


def _row(rows, archive):
    return next(r for r in rows if r["archive"] == archive)


class TheSetHolds(unittest.TestCase):
    def setUp(self):
        self.rows = cis.read_manifest()

    def test_no_problem(self):
        self.assertEqual(cis.check(self.rows), [])

    def test_it_is_the_approved_set(self):
        indexes = sorted(r["archive"] for r in self.rows if r["kind"] in cis.INDEX_KINDS)
        hg = "+hg38_HGDP1kGP"
        self.assertEqual(indexes, sorted(
            [f"20-NRR_3_hg38{p}" for p in ("", hg, "+hg38_HPRC", "+hg38_mega")]
            + [f"{name}{p}" for p in ("", hg) for name in (
                "20-NGG_3_hg38", "20-NRG_3_hg38", "20-NNNNNNNN_3_hg38", "22-NNGRRN_3_hg38",
                "TTTN-20_3_hg38", "TTTN-23_3_hg38")]))
        self.assertFalse([a for a in indexes if a.startswith(("NNN_", "20-NNN_"))],
                         "NNN_3 is not built again for 2.10.0")
        self.assertEqual({r["bmax"] for r in self.rows if r["kind"] in cis.INDEX_KINDS}, {"3"})

    def test_each_variant_panel_ships_its_dictionaries_once(self):
        for panel, companions in (("hg38_HGDP1kGP", ("bundle", "genotypes", "indel-genotypes")),
                                  ("hg38_HPRC", ("bundle", "genotypes", "indel-genotypes")),
                                  ("hg38_mega", ("bundle",))):
            got = sorted(r["kind"] for r in self.rows if r["panel"] == panel
                         and r["kind"] not in cis.INDEX_KINDS + ("build-inputs",))
            self.assertEqual(got, sorted(companions), panel)

    def test_nnn_3_points_to_the_general_index(self):
        for rest in ("hg38", "hg38+hg38_HGDP1kGP"):
            self.assertIn(f"NNN_3_{rest}",
                          _row(self.rows, f"20-NNNNNNNN_3_{rest}")["replaces"].split(","))

    def test_the_commands_print(self):
        import contextlib
        import io
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(cis.main(["--commands"]), 0)
        text = out.getvalue()
        self.assertLess(text.index("# build"), text.index("# publish"))
        self.assertIn("$RUN crisprme.py publish-index --index "
                      "$HGDP/genome_library/20-NRR_3_hg38+hg38_HGDP1kGP", text)

    def test_what_commands_prints_is_a_script_bash_reads(self):
        """--commands prints every build command, then every publish command, for the rebuild
        to run. The build column of the kits and of mega's index ended in a note in brackets,
        and the companions' columns held prose, so bash refused the script (the stage-2b
        review of lane D, D-R1-8): a note now has its own column, and a companion, which its
        panel's first build writes and first publish uploads, has no command ("-"). (Its
        syntax alone: TheRebuildRuns runs it.)"""
        import contextlib
        import io
        import subprocess
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(cis.main(["--commands"]), 0)
        script = out.getvalue()
        got = subprocess.run(["bash", "-n"], input=script, text=True, capture_output=True)
        self.assertEqual((got.returncode, got.stderr), (0, ""))
        for r in self.rows:
            companion = r["kind"] in ("bundle", "genotypes", "indel-genotypes")
            for what in ("build", "publish"):
                self.assertEqual(r[what] == "-", companion, (r["archive"], what))
                if not companion:
                    one = subprocess.run(["bash", "-n"], input=r[what], text=True,
                                         capture_output=True)
                    self.assertEqual((one.returncode, one.stderr), (0, ""), r["archive"])
        self.assertIn("lane M's recipe", _row(self.rows, "20-NRR_3_hg38+hg38_mega")["notes"])
        self.assertIn("written by its first build",
                      _row(self.rows, "dictionaries-bundle_hg38_HGDP1kGP")["notes"])


class TheSearchPamsOfTheSetShip(unittest.TestCase):
    """An index of the set serves the PAMs its own holds, pruned: SaCas9's 22-nt NNGRRN index
    the CEP290 guides' NNGRRT, Cas12a's TTTN indexes TTTV. CRISPRme shipped the PAM files the
    set builds from, but none for a 22-nt NNGRRT or a 20-nt TTTV search, so the web, which
    lists the working directory's PAMs/, offered neither (the stage-2b review of lane D,
    D-R1-7). Each ships, of its index's guide length, a PAM that index serves."""

    SEARCHES = {"22-NNGRRN_3_hg38": ("22bp-NNGRRT-SaCas9.txt", "21bp-NNGRRT-SaCas9.txt"),
                "TTTN-20_3_hg38": ("20bp-TTTV-Cas12a.txt",),
                "TTTN-23_3_hg38": ("23bp-TTTV-Cas12a.txt",)}

    def test_each_ships_and_its_index_serves_it(self):
        rows = {r["archive"]: r for r in cis.read_manifest()}
        codes = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "V": "ACG", "N": "ACGT"}
        for archive, pam_files in self.SEARCHES.items():
            r = rows[archive]
            for name in pam_files:
                path = os.path.join(ROOT, "PAMs", name)
                self.assertTrue(os.path.isfile(path), name)
                pam, side, guide = cis.read_pam_file(path)
                self.assertEqual(side, r["side"], name)
                self.assertTrue(all(set(codes[c]) <= set(codes[i])
                                    for c, i in zip(pam, r["pam"])) and len(pam) == len(r["pam"]),
                                (name, pam, r["pam"]))
                # (a guide of the index's length, or one it serves within the web's bulges)
                self.assertTrue(cis.index_fit.fits(int(r["spacer"]), int(r["bmax"]), guide, 1, 1),
                                name)


class TheCheckFindsWhatItIsFor(unittest.TestCase):
    """Each mutation of the set is caught (and the set itself passes: TheSetHolds)."""

    def setUp(self):
        self.rows = cis.read_manifest()

    def mutated(self, which, **changes):
        rows = copy.deepcopy(self.rows)
        _row(rows, which).update(changes)
        return cis.check(rows)

    def test_a_name_build_index_only_would_not_give(self):
        problems = self.mutated("TTTN-23_3_hg38", archive="23-TTTN_3_hg38")
        self.assertTrue([p for p in problems if "build-index-only names this build TTTN-23_3_hg38"
                         in p], problems)

    def test_an_n_that_is_not_the_bulges_plus_one(self):
        problems = self.mutated("20-NGG_3_hg38", bdna="1", brna="1")
        self.assertTrue([p for p in problems if "is not max(bDNA 1, bRNA 1) + 1" in p], problems)

    def test_the_web_s_two_bulges_must_fit(self):
        # a 23-nt guide on an index of N = 2 is refused at 2 + 2 (and its name says N = 2)
        problems = self.mutated("TTTN-23_3_hg38", bmax="2", bdna="1", brna="1",
                                archive="TTTN-23_2_hg38")
        self.assertTrue([p for p in problems if "does not serve its own 23-nt guide" in p],
                        problems)

    def test_a_pam_file_that_says_otherwise(self):
        problems = self.mutated("22-NNGRRN_3_hg38", spacer="21")
        self.assertTrue([p for p in problems if "holds NNGRRN 3' 22 nt" in p], problems)

    def test_a_panel_without_its_bundle(self):
        rows = [r for r in self.rows if r["archive"] != "dictionaries-bundle_hg38_HPRC"]
        self.assertIn("panel hg38_HPRC: 0 bundle archive(s), not 1", cis.check(rows))

    def test_a_genotyped_panel_without_its_store(self):
        rows = [r for r in self.rows if r["archive"] != "genotypes_hg38_HGDP1kGP"]
        self.assertIn("panel hg38_HGDP1kGP: 0 genotypes archive(s), not 1", cis.check(rows))

    def test_an_archive_of_2_8_1_nothing_replaces(self):
        problems = self.mutated("20-NNNNNNNN_3_hg38+hg38_HGDP1kGP",
                                replaces="NNNNNNNN_3_hg38+hg38_HGDP1kGP")
        self.assertTrue([p for p in problems if p.startswith("NNN_3_hg38+hg38_HGDP1kGP (2.8.1): "
                                                             "replaced by 0")], problems)

    def test_a_repo_of_2_8_1(self):
        problems = self.mutated("20-NRR_3_hg38", repo="lucapinello/crisprme-data")
        self.assertTrue([p for p in problems if "repo lucapinello/crisprme-data" in p], problems)

    def test_a_variant_index_without_its_reference_index(self):
        rows = [r for r in self.rows if r["archive"] != "22-NNGRRN_3_hg38"]
        self.assertTrue([p for p in cis.check(rows)
                         if "its reference index 22-NNGRRN_3_hg38" in p])

    def test_a_note_in_a_command(self):
        kit = _row(self.rows, "build-inputs/hg38_HPRC")
        problems = self.mutated("build-inputs/hg38_HPRC",
                                build=kit["build"] + " (the enriched genome of the panel's first "
                                "build, 20-NRR_3_hg38+hg38_HPRC)")
        self.assertTrue([p for p in problems if p.startswith(
            "build-inputs/hg38_HPRC: its build command is not one bash can run")], problems)

    def test_a_command_for_a_companion(self):
        problems = self.mutated("genotypes_hg38_HPRC", publish="crisprme.py publish-index")
        self.assertTrue([p for p in problems if p.startswith(
            "genotypes_hg38_HPRC: a companion has no publish command of its own")], problems)



class TheRebuildRuns(unittest.TestCase):
    """The script --commands prints is the rebuild: bash runs it from the repository's root,
    every build command, then every publish command. bash -n took it, but it did not run: each
    variant build cd'd into its panel's working directory, then called its recipe by its
    repository path, which was no longer there; each reference build read its PAM file from a
    working directory nothing had put it in; a kit was packed after a cd that moved every
    command after it, and uploaded from the repository's root, where it is not; and mega's
    recipe read its source VCFs from a folder of the cluster it was written on (the stage-2b
    review of lane D, D-R2-1). The check now runs it dry (check_index_set.dry_run: every
    program stubbed, in temporary working directories), and reads each command for the ways it
    failed."""

    def setUp(self):
        self.rows = cis.read_manifest()

    def mutated(self, which, what, old, new):
        rows = copy.deepcopy(self.rows)
        row = _row(rows, which)
        self.assertIn(old, row[what])
        row[what] = row[what].replace(old, new)
        return rows

    def test_its_script_runs_dry_to_its_end(self):
        self.assertEqual(cis.dry_run(self.rows), [])

    def test_the_script_stops_before_its_first_command_where_a_variable_is_not_set(self):
        import shutil
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "HOME": t,
                   "RUN": "", "HGDP": t, "HPRC": t, "MEGA": t}
            got = subprocess.run([shutil.which("bash"), "-c", cis.commands_script(self.rows)],
                                 cwd=ROOT, env=env, capture_output=True, text=True)
            self.assertNotEqual(got.returncode, 0)
            self.assertIn("REF: the working directory of the reference indexes", got.stderr)
            self.assertEqual(os.listdir(t), [], "a command ran")

    def test_a_recipe_called_by_its_path_after_a_cd(self):
        # (the manifest's form before: cd $HGDP && ... bash seq_script/merge_panels/...)
        rows = self.mutated("20-NGG_3_hg38+hg38_HGDP1kGP", "build", "SIF=$SIF",
                            "cd $HGDP && SIF=$SIF")
        problems = cis.check(rows)
        self.assertIn("20-NGG_3_hg38+hg38_HGDP1kGP: its build command names the repository's "
                      "seq_script/merge_panels/hgdp1kgp_build.sh after a cd, where it is not",
                      problems)
        self.assertTrue([p for p in problems if p.startswith(
            "20-NGG_3_hg38+hg38_HGDP1kGP: its build command moves the directory of every "
            "command after it")], problems)
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the script seq_script/merge_panels/hgdp1kgp_build.sh "
                         "is not there" in p], dry)

    def test_a_build_that_reads_a_pam_file_it_did_not_put_there(self):
        rows = self.mutated("20-NRR_3_hg38", "build",
                            "mkdir -p $REF/PAMs && cp PAMs/20bp-NRR-SpCas9.txt $REF/PAMs/ && ", "")
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "20-NRR_3_hg38: its build command does not first copy its PAM file into its "
            "working directory")])
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the PAM file PAMs/20bp-NRR-SpCas9.txt is not there"
                         in p], dry)

    def test_a_kit_uploaded_from_the_repository_s_root(self):
        """(and the script stops there: set -e did not stop it at the first command of an &&
        list, which this upload is, so each command of the script ends in || exit)"""
        rows = self.mutated("build-inputs/hg38_HPRC", "publish", "$HPRC/hg38_HPRC.tar.gz",
                            "hg38_HPRC.tar.gz")
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the file to upload hg38_HPRC.tar.gz is not there"
                         in p], dry)

    def test_a_kit_packed_after_a_cd_that_moves_the_next_commands(self):
        rows = self.mutated("build-inputs/hg38_HGDP1kGP", "build", "(cd $HGDP && ", "cd $HGDP && ")
        _row(rows, "build-inputs/hg38_HGDP1kGP")["build"] = \
            _row(rows, "build-inputs/hg38_HGDP1kGP")["build"][:-1]
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "build-inputs/hg38_HGDP1kGP: its build command moves the directory of every command "
            "after it")])
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the script seq_script/merge_panels/hprc_build.sh is "
                         "not there" in p], dry)

    def test_an_input_a_recipe_reads_from_the_cluster_it_was_written_on(self):
        # (as mega's recipe read its source VCFs, SRC, from /srv/local/lp698/mega_gw)
        rows = self.mutated("20-NRR_3_hg38+hg38_HPRC", "build", "P1KG=$P1KG ", "")
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "20-NRR_3_hg38+hg38_HPRC: seq_script/merge_panels/hprc_build.sh reads P1KG from "
            "/srv/local/")])

    def test_a_program_of_the_release_image_run_outside_it(self):
        rows = self.mutated("20-NGG_3_hg38", "publish", "$RUN crisprme.py", "crisprme.py")
        self.assertIn("20-NGG_3_hg38: its publish command runs crisprme.py outside the release "
                      "image: $RUN crisprme.py", cis.check(rows))
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "runs outside the release image" in p], dry)

    def test_a_variable_the_rebuild_does_not_set(self):
        rows = self.mutated("20-NRR_3_hg38+hg38_mega", "build", "$MEGA_SRC", "$CK")
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "20-NRR_3_hg38+hg38_mega: its build command reads $CK, which the rebuild does not "
            "set")])

    def test_the_working_directory_of_another_panel(self):
        rows = self.mutated("20-NGG_3_hg38+hg38_HGDP1kGP", "publish", "$HGDP/", "$HPRC/")
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "20-NGG_3_hg38+hg38_HGDP1kGP: its publish command uses $HPRC, the working directory "
            "of another panel")])

    def test_a_script_the_repository_does_not_have(self):
        rows = self.mutated("TTTN-23_3_hg38+hg38_HGDP1kGP", "build", "hgdp1kgp_build.sh",
                            "hgdp1kgp.sh")
        self.assertIn("TTTN-23_3_hg38+hg38_HGDP1kGP: its build command names "
                      "seq_script/merge_panels/hgdp1kgp.sh, which the repository does not have",
                      cis.check(rows))
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the script seq_script/merge_panels/hgdp1kgp.sh is "
                         "not there" in p], dry)

    def test_a_pending_recipe_the_tree_holds(self):
        """mega's build calls lane M's mega_rebuild.sh, which the branch of lane D does not
        hold: the check says so, and the dry run takes its call as done, until the tree holds
        it; then its entry in PENDING_SCRIPTS has to go."""
        import tempfile
        self.assertIn("seq_script/merge_panels/mega_rebuild.sh",
                      _row(self.rows, "20-NRR_3_hg38+hg38_mega")["build"])
        with tempfile.TemporaryDirectory() as t:
            os.symlink(os.path.join(ROOT, "PAMs"), os.path.join(t, "PAMs"))
            panels = os.path.join(t, "seq_script", "merge_panels")
            os.makedirs(panels)
            for n in os.listdir(os.path.join(ROOT, "seq_script", "merge_panels")):
                os.symlink(os.path.join(ROOT, "seq_script", "merge_panels", n),
                           os.path.join(panels, n))
            for script in cis.PENDING_SCRIPTS:
                if not os.path.exists(os.path.join(t, script)):
                    open(os.path.join(t, script), "w").close()
            problems = cis.check(self.rows, root=t)
            self.assertEqual(cis.pending(root=t), [])
        for script in cis.PENDING_SCRIPTS:
            self.assertIn(f"{script}: this tree holds it, so it is no longer pending: take it out "
                          f"of check_index_set.PENDING_SCRIPTS", problems)

    def test_a_text_the_script_cannot_hold(self):
        import unittest.mock
        with unittest.mock.patch.dict(cis.VARIABLES, {"META": "gnomAD's metadata file"}):
            self.assertTrue([p for p in cis.check(self.rows) if p.startswith(
                "VARIABLES['META'] holds one of")])
            dry = cis.dry_run(self.rows)
        self.assertTrue([p for p in dry if "stops" in p], dry)

if __name__ == "__main__":
    unittest.main()

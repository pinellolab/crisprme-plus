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


# hg38's chromosomes chr1-22, X and Y, what the working directories of HPRC and mega hold
# (seq_script/rebuild/README.md, rule 1), and, with the mitochondrion and two of the 431 alt,
# random and chrUn contigs, all of hg38 as `download --what genome` installs it, what the
# reference's and HGDP + 1000 Genomes' hold: as 2.8.1's default index of each panel was built
PRIMARY = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
HG38 = PRIMARY + ("chrM", "chr10_GL383545v1_alt", "chrUn_KI270302v1")
GENOMES = {"REF": HG38, "HGDP": HG38, "HPRC": PRIMARY, "MEGA": PRIMARY}


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

    def test_an_archive_listed_twice(self):
        # (the check found it before; nothing made it happen: the stage-2b review of lane D,
        # D-R3-5)
        rows = self.rows + [copy.deepcopy(_row(self.rows, "TTTN-20_3_hg38"))]
        self.assertIn("TTTN-20_3_hg38: listed 2 times", cis.check(rows))

    # The ways the review of lane D (D-R3-5) found a manifest wrong that the check and its dry
    # run took: each is found now.

    def test_a_kit_uploaded_to_a_repo_of_2_8_1(self):
        # it would replace 2.8.1's kit in the repo 2.8.1 reads (decision C2)
        kit = _row(self.rows, "build-inputs/hg38_HGDP1kGP")
        problems = self.mutated("build-inputs/hg38_HGDP1kGP", publish=kit["publish"].replace(
            "upload lucapinello/crisprme-indexes-format2", "upload lucapinello/crisprme-data", 1))
        self.assertIn("build-inputs/hg38_HGDP1kGP: its publish command uploads to "
                      "lucapinello/crisprme-data, not lucapinello/crisprme-indexes-format2: the "
                      "set goes to the repo of index format 2 alone (decision C2), and an upload "
                      "elsewhere replaces what 2.8.1 reads", problems)

    def test_a_publish_that_names_another_repo(self):
        r = _row(self.rows, "20-NGG_3_hg38")
        problems = self.mutated("20-NGG_3_hg38",
                                publish=r["publish"] + " --hf-repo lucapinello/crisprme-data")
        self.assertIn("20-NGG_3_hg38: its publish command names the repo "
                      "lucapinello/crisprme-data, not lucapinello/crisprme-indexes-format2",
                      problems)

    def test_a_token_in_a_command(self):
        r = _row(self.rows, "20-NGG_3_hg38")
        problems = self.mutated("20-NGG_3_hg38", publish=r["publish"] + " --token hf_x")
        self.assertIn("20-NGG_3_hg38: its publish command passes --token: a token goes in "
                      "HF_TOKEN, never in the manifest", problems)

    def test_a_variant_index_published_with_its_per_sample_dictionaries(self):
        # its bundle would carry dictionaries_<panel> (152 GB of HGDP + 1000 Genomes), and
        # differ in fingerprint from the one the panel's other indexes name
        r = _row(self.rows, "20-NGG_3_hg38+hg38_HGDP1kGP")
        problems = self.mutated("20-NGG_3_hg38+hg38_HGDP1kGP",
                                publish=r["publish"].replace(" --dictless", ""))
        self.assertTrue([p for p in problems if p.startswith(
            "20-NGG_3_hg38+hg38_HGDP1kGP: its publish command does not pass --dictless")],
            problems)

    def test_a_recipe_input_left_to_a_default_under_the_working_directory(self):
        # hgdp1kgp_build.sh's RP=${RP:-$B}: the re-phased BCFs from the working directory
        r = _row(self.rows, "TTTN-23_3_hg38+hg38_HGDP1kGP")
        problems = self.mutated("TTTN-23_3_hg38+hg38_HGDP1kGP",
                                build=r["build"].replace("RP=$RP ", ""))
        self.assertIn("TTTN-23_3_hg38+hg38_HGDP1kGP: seq_script/merge_panels/hgdp1kgp_build.sh "
                      "reads RP from $B unless it is given it: its build command gives RP=",
                      problems)

    def test_an_input_no_command_reads(self):
        # mega's merged VCFs, dropped from its recipe's call: the recipe would stop on its
        # usage, after hours of the other panels' builds
        r = _row(self.rows, "20-NRR_3_hg38+hg38_mega")
        problems = self.mutated("20-NRR_3_hg38+hg38_mega",
                                build=r["build"].replace(" $MEGA_SRC", ""))
        self.assertIn("VARIABLES['MEGA_SRC']: no command reads $MEGA_SRC: a command lost an "
                      "input it needs, or the entry goes", problems)



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

    def test_megas_recipe_is_the_trees(self):
        """mega's build calls lane M's mega_rebuild.sh, which the branch of lane D did not
        hold, so the check named it pending and the dry run took its call as done. The
        integration of stage 2b brought it: no script is pending, and the dry run checks
        mega's call (the script there, its working directory and its source VCFs) as it
        checks every other recipe's."""
        build = _row(self.rows, "20-NRR_3_hg38+hg38_mega")["build"]
        self.assertIn("bash seq_script/merge_panels/mega_rebuild.sh $MEGA $MEGA_SRC", build)
        self.assertTrue(os.path.isfile(
            os.path.join(ROOT, "seq_script", "merge_panels", "mega_rebuild.sh")))
        self.assertEqual(cis.PENDING_SCRIPTS, {})
        self.assertEqual(cis.pending(), [])
        rows = self.mutated("20-NRR_3_hg38+hg38_mega", "build", "mega_rebuild.sh",
                            "mega_rebuilt.sh")
        dry = cis.dry_run(rows)
        self.assertTrue([p for p in dry if "the script seq_script/merge_panels/mega_rebuilt.sh "
                         "is not there" in p], dry)

    def test_a_pending_recipe_the_tree_holds(self):
        """A recipe a branch not merged yet brings (PENDING_SCRIPTS): the check names it, and
        the dry run takes its call as done, until the tree holds it; then its entry has to
        go."""
        import tempfile
        from unittest import mock
        script = "seq_script/merge_panels/mega_rebuild.sh"
        with mock.patch.dict(cis.PENDING_SCRIPTS, {script: "a recipe of a branch"}):
            with tempfile.TemporaryDirectory() as t:
                os.symlink(os.path.join(ROOT, "PAMs"), os.path.join(t, "PAMs"))
                panels = os.path.join(t, "seq_script", "merge_panels")
                os.makedirs(panels)
                for n in os.listdir(os.path.join(ROOT, "seq_script", "merge_panels")):
                    if n != "mega_rebuild.sh":
                        os.symlink(os.path.join(ROOT, "seq_script", "merge_panels", n),
                                   os.path.join(panels, n))
                self.assertEqual(cis.pending(root=t), [
                    f"{script}: a recipe of a branch, which this tree does not hold yet (the "
                    f"merge of its branch brings it): the dry run takes its call as done"])
                self.assertNotIn(f"{script}: this tree holds it, so it is no longer pending: "
                                 f"take it out of check_index_set.PENDING_SCRIPTS",
                                 cis.check(self.rows, root=t))
                open(os.path.join(t, script), "w").close()
                problems = cis.check(self.rows, root=t)
                self.assertEqual(cis.pending(root=t), [])
        self.assertIn(f"{script}: this tree holds it, so it is no longer pending: take it out "
                      f"of check_index_set.PENDING_SCRIPTS", problems)

    def run_on_genomes(self, genomes=None, mega_src=None):
        """The rebuild's script run from the repository's root on working directories whose
        Genomes/hg38 hold the FASTAs ``genomes`` names ({variable: [contig, ...]}; GENOMES's
        where it names none), mega's merged VCFs of ``mega_src`` (chr1-22, X and Y if None),
        every other variable set, and no program of the image on PATH: (exit status, stderr,
        whether a command ran: its first makes $REF/PAMs)."""
        import shutil
        import subprocess
        import tempfile
        genomes = genomes or {}
        with tempfile.TemporaryDirectory() as t:
            env = {"PATH": "/usr/bin:/bin", "HOME": t, "LC_ALL": "C", "RUN": ""}
            for v in ("REF", "HGDP", "HPRC", "MEGA"):
                env[v] = os.path.join(t, v)
                os.makedirs(os.path.join(env[v], "Genomes", "hg38"))
                for c in genomes.get(v, GENOMES[v]):
                    open(os.path.join(env[v], "Genomes", "hg38", c + ".fa"), "w").close()
            for v in ("SIF", "RP", "META", "HPRC_SRC", "P1KG", "MEGA_SRC"):
                env[v] = os.path.join(t, v)
                os.makedirs(env[v])
            for c in PRIMARY if mega_src is None else mega_src:
                open(os.path.join(env["MEGA_SRC"], f"mega.{c}.afmax.vcf.gz"), "w").close()
            got = subprocess.run([shutil.which("bash"), "-c", cis.commands_script(self.rows)],
                                 cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
            return got.returncode, got.stderr, os.path.exists(os.path.join(env["REF"], "PAMs"))

    def test_a_genome_folder_of_all_of_hg38_for_hprc_or_mega_stops_it_first(self):
        """The README said to put all of hg38 in every working directory, as download
        installs it; mega's recipe refuses any FASTA but the chromosomes of its merged VCFs,
        so the rebuild stopped at mega's build after hours of the other panels' (the stage-2b
        review of lane D, D-R3-1), and HPRC's indexes would have held 455 contigs where 2.8.1's
        hold 24. The script stops on it before its first command."""
        for v, panel in (("HPRC", "hg38_HPRC"), ("MEGA", "hg38_mega")):
            with self.subTest(panel=panel):
                rc, err, ran = self.run_on_genomes({v: HG38})
                self.assertEqual(rc, 2, err)
                self.assertFalse(ran, "a command ran")
                self.assertIn(f"/{v}/Genomes/hg38 (${v}) holds chr10_GL383545v1_alt chrM "
                              f"chrUn_KI270302v1: the working directory of {panel} holds the "
                              f"chromosomes chr1-22, X and Y of hg38 alone", err)

    def test_hgdp_s_genome_folder_of_chr1_22_x_and_y_alone_stops_it_first(self):
        # its indexes would hold 24 contigs, where 2.8.1's default index of the panel (NRR, and
        # NGG and NNN) holds all 455 of hg38 (2.8.1's NRG index of the panel held 24)
        rc, err, ran = self.run_on_genomes({"HGDP": PRIMARY})
        self.assertEqual(rc, 2, err)
        self.assertFalse(ran, "a command ran")
        self.assertIn("/HGDP/Genomes/hg38 ($HGDP) holds 24 FASTA(s), not the 27 of $REF's: the "
                      "working directory of hg38_HGDP1kGP holds all of hg38, as the reference's "
                      "does", err)

    def test_a_panel_without_a_chromosome_of_the_reference_s_stops_it_first(self):
        rc, err, ran = self.run_on_genomes({"HPRC": PRIMARY[:-1]})
        self.assertEqual(rc, 2, err)
        self.assertFalse(ran, "a command ran")
        self.assertIn("/HPRC/Genomes/hg38 ($HPRC) holds chr1 chr10", err)
        self.assertIn("the working directory of hg38_HPRC holds the chromosomes chr1-22, X and Y "
                      "of the reference's", err)

    def test_mega_s_merged_vcfs_of_other_chromosomes_stop_it_first(self):
        rc, err, ran = self.run_on_genomes(mega_src=PRIMARY[:-1])
        self.assertEqual(rc, 2, err)
        self.assertFalse(ran, "a command ran")
        self.assertIn("($MEGA_SRC) holds mega.<chr>.afmax.vcf.gz of chr1 chr10", err)
        self.assertIn("the recipe builds every chromosome of its sources, and refuses a genome "
                      "folder that holds the FASTA of another or lacks one of theirs", err)

    def test_a_test_genome_of_one_chromosome_runs(self):
        # (as the D7 readiness audit's chr22 dry run of every index type: its first command
        # runs, and stops here, where no crisprme.py is on PATH)
        rc, err, ran = self.run_on_genomes(
            {"REF": ("chr22", "chrM"), "HGDP": ("chr22", "chrM"), "HPRC": ("chr22",),
             "MEGA": ("chr22",)}, mega_src=("chr22",))
        self.assertTrue(ran, err)
        self.assertIn("rebuild: the build of 20-NRR_3_hg38 failed", err)

    def test_a_kit_packed_into_the_repository_s_root(self):
        # (the dry run found it before; nothing made it happen: D-R3-5)
        kit = _row(self.rows, "build-inputs/hg38_HGDP1kGP")["build"]
        packed = kit.replace("(cd $HGDP && tar -cf - Genomes/", "tar -cf - $HGDP/Genomes/") \
            .replace(" Genomes/hg38+hg38_HGDP1kGP_INDELS", " $HGDP/Genomes/hg38+hg38_HGDP1kGP_INDELS") \
            .rstrip(")")
        self.assertNotEqual(packed, kit)
        rows = copy.deepcopy(self.rows)
        _row(rows, "build-inputs/hg38_HGDP1kGP")["build"] = packed
        dry = cis.dry_run(rows)
        self.assertIn("the rebuild's script, run dry, writes into the repository's root: "
                      "hg38_HGDP1kGP.SHA256SUMS, hg38_HGDP1kGP.tar.gz", dry)

    def test_a_text_the_script_cannot_hold(self):
        import unittest.mock
        with unittest.mock.patch.dict(cis.VARIABLES, {"META": "gnomAD's metadata file"}):
            self.assertTrue([p for p in cis.check(self.rows) if p.startswith(
                "VARIABLES['META'] holds one of")])
            dry = cis.dry_run(self.rows)
        self.assertTrue([p for p in dry if "stops" in p], dry)

if __name__ == "__main__":
    unittest.main()

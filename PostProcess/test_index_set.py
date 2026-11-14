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
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_spec = importlib.util.spec_from_file_location(
    "check_index_set", os.path.join(ROOT, "seq_script", "rebuild", "check_index_set.py"))
cis = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cis)


def _row(rows, archive):
    return next(r for r in rows if r["archive"] == archive)


def printed(*args):
    """check_index_set's main run with ``args``: (exit status, stdout, stderr)."""
    import contextlib
    import io
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            rc = cis.main(list(args))
        except SystemExit as stop:
            rc = stop.code
    return rc, out.getvalue(), err.getvalue()


# hg38's chromosomes chr1-22, X and Y, what the working directories of HPRC and mega hold
# (seq_script/rebuild/README.md, rule 1), and, with the mitochondrion and two of the 431 alt,
# random and chrUn contigs, all of hg38 as `download --what genome` installs it, what the
# reference's and HGDP + 1000 Genomes' hold: as 2.8.1's default index of each panel was built
PRIMARY = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
HG38 = PRIMARY + ("chrM", "chr10_GL383545v1_alt", "chrUn_KI270302v1")
GENOMES = {"REF": HG38, "HGDP": HG38, "HPRC": PRIMARY, "MEGA": PRIMARY}
# What mega's recipe (seq_script/merge_panels/mega_rebuild.sh) leaves in its working directory
# once it has written every part of the dataset, before its last step, its checks
# (verify_mega_chrom.py): a build that failed them leaves all of them, its staging folder gone,
# beside the build VCF of each chromosome of its sources (VCFs/hg38_mega/mega.<chr>.build.vcf.gz)
MEGA = "20-NRR_3_hg38+hg38_mega"
MEGA_PARTS = (f"genome_library/{MEGA}", f"genome_library/{MEGA}_INDELS",
              "Genomes/hg38+hg38_mega", "Genomes/hg38+hg38_mega_INDELS",
              "Dictionaries/log_indels_hg38_mega", "Dictionaries/registry_hg38_mega",
              "Dictionaries/indel_af_hg38_mega", "samplesIDs/hg38_mega.samplesID.txt")


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
        out = {}
        for half in ("build", "publish"):
            out[half] = io.StringIO()
            with contextlib.redirect_stdout(out[half]):
                self.assertEqual(cis.main(["--commands", half]), 0)
        self.assertIn("# build\n", out["build"].getvalue())
        self.assertIn("$RUN crisprme.py publish-index --index "
                      "$HGDP/genome_library/20-NRR_3_hg38+hg38_HGDP1kGP",
                      out["publish"].getvalue())

    def test_what_commands_prints_is_a_script_bash_reads(self):
        """--commands prints the rebuild's build commands, or its publish commands, for the
        rebuild to run. The build column of the kits and of mega's index ended in a note in
        brackets, and the companions' columns held prose, so bash refused the script (the
        stage-2b review of lane D, D-R1-8): a note now has its own column, and a companion,
        which its panel's first build writes and first publish uploads, has no command ("-").
        (Its syntax alone: TheRebuildRuns runs it.)"""
        import contextlib
        import io
        import subprocess
        for half in ("build", "publish"):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(cis.main(["--commands", half]), 0)
            script = out.getvalue()
            got = subprocess.run(["bash", "-n"], input=script, text=True, capture_output=True)
            self.assertEqual((got.returncode, got.stderr), (0, ""), half)
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

    def run_on_genomes(self, genomes=None, mega_src=None, panels=None, unset=(), preamble=False,
                       setup=None):
        """The rebuild's script run from the repository's root on working directories whose
        Genomes/hg38 hold the FASTAs ``genomes`` names ({variable: [contig, ...]}; GENOMES's
        where it names none), mega's merged VCFs of ``mega_src`` (chr1-22, X and Y if None),
        every other variable set but those ``unset`` names, and no program of the image on PATH:
        the build half, or its part of the working directories ``panels`` names (--panels), or,
        with ``preamble``, what of it comes before its first row: (exit status, stderr, whether a
        command ran: each first makes the PAMs/ of its working directory). ``setup``, called with
        the environment once the folders are made, changes them first."""
        import shutil
        import subprocess
        import tempfile
        genomes = genomes or {}
        script = cis.commands_script(self.rows, panels=panels)
        if preamble:
            script = script[:script.index("\n# build\n") + len("\n# build\n")]
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
            folders = [env[v] for v in ("REF", "HGDP", "HPRC", "MEGA")]
            if setup:
                setup(env)
            for v in unset:
                del env[v]
            got = subprocess.run([shutil.which("bash"), "-c", script],
                                 cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
            return got.returncode, got.stderr, any(os.path.exists(os.path.join(w, "PAMs"))
                                                   for w in folders)

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

    def test_a_node_s_part_of_the_build_half_checks_its_genome_folders_as_the_whole_does(self):
        """A node runs its part of a half (--panels), but a part of the build half read only
        the variables of its rows, and so no $REF where it named none of the reference's: a part
        of HGDP + 1000 Genomes checked only that $HGDP's genome folder held a FASTA, and built on
        one of chr22 alone, which the whole half refuses; a part of HPRC or mega checked only
        that theirs held chr1-22, X or Y and no other; and a part of the reference's checked
        nothing (the PP review, round 2: A2-PP2-1, B2-1). A part checks them as the whole half
        does, against $REF's, which it asks for, and stops before its first command."""
        chr22 = ("chr22",)
        for panels, genomes, said in (
                (["HGDP"], {"HGDP": chr22},
                 "/HGDP/Genomes/hg38 ($HGDP) holds 1 FASTA(s), not the 27 of $REF's: the working "
                 "directory of hg38_HGDP1kGP holds all of hg38, as the reference's does"),
                (["HPRC"], {"HPRC": chr22},
                 "/HPRC/Genomes/hg38 ($HPRC) holds chr22, not chr1 chr10"),
                (["MEGA"], {"MEGA": chr22},
                 "/MEGA/Genomes/hg38 ($MEGA) holds chr22, not chr1 chr10"),
                (["HPRC", "MEGA"], {"HPRC": chr22, "MEGA": chr22},
                 "/HPRC/Genomes/hg38 ($HPRC) holds chr22, not chr1 chr10"),
                (["REF"], {"REF": ("chrM",)},
                 "/REF/Genomes/hg38 ($REF) holds no FASTA of chr1-22, X or Y")):
            with self.subTest(panels=panels):
                # (mega's merged VCFs those of its genome folder, which its own check reads)
                rc, err, ran = self.run_on_genomes(genomes, mega_src=genomes.get("MEGA"),
                                                   panels=panels)
                self.assertEqual(rc, 2, err)
                self.assertFalse(ran, "a command ran")
                self.assertIn(said, err)
        # $REF is asked for, and read alone: a part of HGDP + 1000 Genomes without it stops first
        rc, err, ran = self.run_on_genomes(panels=["HGDP"], unset=("REF",))
        self.assertNotEqual(rc, 0)
        self.assertFalse(ran, "a command ran")
        self.assertIn("REF: " + cis.REF_READ, err)
        self.assertIsNone(cis._VARIABLE_TEXT.search(cis.REF_READ))
        script = cis.commands_script(self.rows, panels=["HGDP"])
        self.assertNotIn("$REF", cis.half_lines(script))
        # on genome folders as rule 1 says, each node's part goes past its checks
        for panels in (["HGDP"], ["REF", "HPRC", "MEGA"], ["HPRC"], ["MEGA"], ["REF"]):
            with self.subTest(panels=panels, folders="as rule 1 says"):
                self.assertEqual(self.run_on_genomes(panels=panels, preamble=True), (0, "", False))

    def test_a_ref_of_links_into_a_folder_it_checks_stops_it_first(self):
        """The README said that on the node of HGDP + 1000 Genomes a $REF of links to $HGDP's
        own FASTAs would do, as the D7 plan's stub was (ln -sf "$HGDP/Genomes/hg38/"*.fa
        "$REF/Genomes/hg38/"): against it $HGDP's folder was checked against itself, a FASTA
        missing there missing from $REF's too (its link dangles, and a FASTA is listed only where
        it is there), and the node's part went past its checks on a genome folder of chr1-22, X
        and Y alone, as $REF=$HGDP did (the PP review, round 3: B3-1). $REF holds hg38 of its own:
        each half stops on a $REF that is a folder it checks, holds a link into one, or a link to
        nothing; one of links to another folder of hg38 goes past, as a copy does."""
        import shutil

        def into_hgdp(env, stub_while_whole=True):
            """$REF/Genomes/hg38 made the D7 plan's stub of links to $HGDP's FASTAs, while $HGDP
            held all of hg38 (or after it lost the others), then $HGDP left with chr1-22, X and
            Y."""
            ref, hgdp = (os.path.join(env[v], "Genomes", "hg38") for v in ("REF", "HGDP"))
            if not stub_while_whole:
                self.drop(hgdp)
            shutil.rmtree(ref)
            os.makedirs(ref)
            for name in os.listdir(hgdp):
                os.symlink(os.path.join(hgdp, name), os.path.join(ref, name))
            self.drop(hgdp)

        def into_hgdp_later(env):
            into_hgdp(env, stub_while_whole=False)

        def is_hgdp(env):
            self.drop(os.path.join(env["HGDP"], "Genomes", "hg38"))
            env["REF"] = env["HGDP"]

        def dangling(env):
            os.remove(os.path.join(env["REF"], "Genomes", "hg38", "chrM.fa"))
            os.symlink(os.path.join(env["REF"], "gone", "chrM.fa"),
                       os.path.join(env["REF"], "Genomes", "hg38", "chrM.fa"))

        def into_another(env):
            """$REF of links to a folder of all of hg38 of its own (on the node, the one $HGDP's
            were copied from), $HGDP whole."""
            other, ref = os.path.join(env["REF"], "..", "other"), os.path.join(
                env["REF"], "Genomes", "hg38")
            shutil.move(ref, other)
            os.makedirs(ref)
            for name in os.listdir(other):
                os.symlink(os.path.join(os.path.abspath(other), name), os.path.join(ref, name))

        said_into = ("/REF/Genomes/hg38 ($REF) is {}/Genomes/hg38 ($HGDP), or holds links to its "
                     "FASTAs: that folder would be checked against itself, a FASTA missing there "
                     "missing from $REF's too. $REF's Genomes/hg38 holds hg38 of its own")
        for panels in (["HGDP"], None):
            for how, setup in (("the D7 stub", into_hgdp), ("a stub made later", into_hgdp_later),
                               ("$HGDP itself", is_hgdp)):
                with self.subTest(panels=panels, ref=how):
                    rc, err, ran = self.run_on_genomes(panels=panels, setup=setup, preamble=True)
                    self.assertEqual(rc, 2, err)
                    self.assertFalse(ran, "a command ran")
                    self.assertIn("/Genomes/hg38 ($REF) is ", err)
                    self.assertIn(said_into.split("{}")[1], err)
            with self.subTest(panels=panels, ref="a link to nothing"):
                rc, err, ran = self.run_on_genomes(panels=panels, setup=dangling, preamble=True)
                self.assertEqual(rc, 2, err)
                self.assertIn("/REF/Genomes/hg38/chrM.fa ($REF) is a link to nothing", err)
            with self.subTest(panels=panels, ref="links to another folder of hg38"):
                self.assertEqual(self.run_on_genomes(panels=panels, setup=into_another,
                                                     preamble=True), (0, "", False))
        # (ml008's parts: $REF the reference's working directory, $HGDP's folder links to its
        # FASTAs, which is no $REF of links into a folder it checks)
        def hgdp_into_ref(env):
            ref, hgdp = (os.path.join(env[v], "Genomes", "hg38") for v in ("REF", "HGDP"))
            shutil.rmtree(hgdp)
            os.makedirs(hgdp)
            for name in os.listdir(ref):
                os.symlink(os.path.join(ref, name), os.path.join(hgdp, name))
        for panels in (["REF", "HPRC", "MEGA"], None):
            with self.subTest(panels=panels, hgdp="links to $REF's"):
                self.assertEqual(self.run_on_genomes(panels=panels, setup=hgdp_into_ref,
                                                     preamble=True), (0, "", False))

    @staticmethod
    def drop(folder):
        """Leaves genome folder ``folder`` with its FASTAs of chr1-22, X and Y alone, a folder
        the build half refuses for HGDP + 1000 Genomes (rule 1)."""
        for name in os.listdir(folder):
            if name[: -len(".fa")] not in PRIMARY:
                os.remove(os.path.join(folder, name))
        assert sorted(n[: -len(".fa")] for n in os.listdir(folder)) == sorted(PRIMARY)

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

class TheRebuildHasTwoHalves(unittest.TestCase):
    """The script --commands printed ran every build command, then every publish command, as
    the README ran it (bash rebuild.sh): about 47 hours of builds, then an upload of every
    archive before any gate or the PI's go, where a token or the cached login of the node held
    one; a re-run, the way back from a failed upload, stopped at mega's build, whose recipe
    refuses a working directory that holds any part of the dataset (the D7 readiness audit,
    R-PUB-2, CAP-2, CAP-5; PP-2 of the release). Its two halves are scripts of their own now:
    --commands build, which publishes nothing, and --commands publish, run after the gates and
    the go; each skips a row it ran to its end before (its mark), and the build half passes a kit
    whose SHA256SUMS holds. It passed mega's build too where every part its recipe makes was
    there, but the recipe writes them all before its checks, its last step: a build that failed
    them was taken for done, and a re-run exited 0 (the PP-2 review, A1-PP2-1). mega's build is
    passed on its mark alone, and its recipe refuses a working directory that holds a part of
    the dataset. The lines of the build half are the manifest's, as the rebuild ran them."""

    def setUp(self):
        self.rows = cis.read_manifest()

    def printed(self, *args):
        return printed(*args)

    def test_a_half_is_asked_for_by_name(self):
        for args in (["--commands"], ["--commands", "all"]):
            with self.subTest(args=args):
                rc, out, err = self.printed(*args)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn("--commands build (every build command; it publishes nothing) or "
                              "--commands publish", err)

    def test_the_build_half_publishes_nothing_and_the_publish_half_builds_nothing(self):
        rc, build, _ = self.printed("--commands", "build")
        self.assertEqual(rc, 0)
        rc, publish, _ = self.printed("--commands", "publish")
        self.assertEqual(rc, 0)
        for word in ("publish-index", "hf upload", "\n# publish\n"):
            self.assertNotIn(word, build)
        for word in ("build-index-only", "seq_script/merge_panels/", "tar -cf", "\n# build\n"):
            self.assertNotIn(word, publish)
        for r in self.rows:
            for half, text in (("build", build), ("publish", publish)):
                if r[half] != cis.NO_COMMAND:
                    self.assertIn(f"\n{r[half]} || {{ echo \"rebuild: the {half} of {r['archive']} "
                                  f"failed\"", text)

    def test_its_build_lines_are_those_the_rebuild_ran(self):
        """The lines of the build half, the if, else and mark around each left out, are the 19
        the D7 readiness audit checked and the rebuild of 2.10.0's indexes ran (md5
        7389f7e234d4e617860f8f9b34abdc49, as stage.sh cut them from the one script between its
        "# build" and "# publish" lines): the split changes no build. A change to a build
        command changes this md5, and is a change to the rebuild."""
        import hashlib
        lines = cis.half_lines(cis.commands_script(self.rows, half="build"))
        self.assertEqual(len(lines.splitlines()), 19)
        self.assertEqual(hashlib.md5(lines.encode()).hexdigest(),
                         "7389f7e234d4e617860f8f9b34abdc49")

    def test_a_build_command_that_publishes_is_found(self):
        rows = copy.deepcopy(self.rows)
        row = _row(rows, "20-NGG_3_hg38")
        row["build"] += " && $RUN crisprme.py publish-index --index $REF/genome_library/20-NGG_3_hg38"
        self.assertTrue([p for p in cis.check(rows) if p.startswith(
            "20-NGG_3_hg38: its build command publishes (crisprme.py publish-index): the build "
            "half of the rebuild publishes nothing")], cis.check(rows))
        self.assertTrue([p for p in cis.dry_run(rows) if p.startswith(
            "the rebuild's build half, run dry, publishes: publish-index")])

    def test_a_re_run_that_does_again_what_it_did_is_found(self):
        from unittest import mock
        wrapped = cis._row_lines

        def unwrapped(r, half):  # (each line alone, as the one script ran it)
            return wrapped(r, half)[1:2]
        with mock.patch.object(cis, "_row_lines", unwrapped):
            dry = cis.dry_run(self.rows)
        self.assertTrue([p for p in dry if p.startswith(
            "the rebuild's build half, run dry again, does again what it did")], dry)

    def test_a_row_is_done_by_its_mark_or_a_kit_by_its_list_and_by_nothing_else(self):
        """What a recipe leaves does not tell that it ran to its end: the build half's test of
        each row is its mark (with its index for an index), and for a kit its SHA256SUMS."""
        script = cis.commands_script(self.rows, half="build")
        tests = [line for line in script.splitlines() if line.startswith("if ")]
        self.assertEqual(len(tests), sum(r["build"] != cis.NO_COMMAND for r in self.rows))
        for r in self.rows:
            if r["build"] == cis.NO_COMMAND:
                continue
            w = "\"$" + cis.WORKDIR[r["panel"]] + "\""
            done = f"_done {w} build.{r['archive'].replace('/', '_')}"
            if r["kind"] in cis.INDEX_KINDS:
                done += f" {r['archive']}"
            else:
                done += f" || _kit_done {w} {r['archive'].split('/', 1)[1]}"
            self.assertIn(f"if {done}; then _skip \"build of {r['archive']}\"; else", tests)

    def run_build(self, done=(), mega_parts=None, mega_mark=False):
        """The build half run from the repository's root with no program of the image on PATH
        (RUN empty), every row's mark there but mega's build's (unless ``mega_mark``), and
        ``mega_parts`` in $MEGA: "all", what its recipe leaves once it has written every part of
        the dataset, as after a build that failed its checks, its last step (MEGA_PARTS and the
        build VCFs), or "some", the first two parts: (exit status, stdout, stderr)."""
        import shutil
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            env = {"PATH": "/usr/bin:/bin", "HOME": t, "LC_ALL": "C", "RUN": ""}
            for v in ("REF", "HGDP", "HPRC", "MEGA"):
                env[v] = os.path.join(t, v)
                os.makedirs(os.path.join(env[v], "Genomes", "hg38"))
                for c in GENOMES[v]:
                    open(os.path.join(env[v], "Genomes", "hg38", c + ".fa"), "w").close()
            for v in ("SIF", "RP", "META", "HPRC_SRC", "P1KG", "MEGA_SRC"):
                env[v] = os.path.join(t, v)
                os.makedirs(env[v])
            for c in PRIMARY:
                open(os.path.join(env["MEGA_SRC"], f"mega.{c}.afmax.vcf.gz"), "w").close()
            mega = MEGA
            for r in self.rows:
                if r["build"] == cis.NO_COMMAND or (r["archive"] == mega and not mega_mark):
                    continue
                w = env[cis.WORKDIR[r["panel"]]]
                os.makedirs(os.path.join(w, cis.MARKS), exist_ok=True)
                open(os.path.join(w, cis.MARKS, cis._mark_name("build", r["archive"])),
                     "w").close()
                if r["kind"] in cis.INDEX_KINDS:
                    os.makedirs(os.path.join(w, "genome_library", r["archive"]))
            parts = MEGA_PARTS
            for part in parts[:len(parts) if mega_parts == "all" else 2] if mega_parts else ():
                path = os.path.join(env["MEGA"], part)
                os.makedirs(os.path.dirname(path) if part.endswith(".txt") else path,
                            exist_ok=True)
                if part.endswith(".txt"):
                    open(path, "w").close()
            if mega_parts == "all":
                for c in PRIMARY:
                    path = os.path.join(env["MEGA"], "VCFs", "hg38_mega", f"mega.{c}.build.vcf.gz")
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    open(path, "w").close()
            got = subprocess.run([shutil.which("bash"), "-c",
                                  cis.commands_script(self.rows, half="build")],
                                 cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
            return got.returncode, got.stdout, got.stderr

    def test_megas_build_is_passed_on_its_mark(self):
        """A re-run once the build half has run mega's recipe to its end (its mark), every
        other row done before too: nothing runs, exit 0."""
        rc, out, err = self.run_build(mega_parts="all", mega_mark=True)
        self.assertEqual((rc, err), (0, ""))
        self.assertIn("rebuild: the build of 20-NRR_3_hg38+hg38_mega: done before, skipped", out)
        self.assertEqual(out.count(": done before, skipped"), 19, out)

    def test_a_mega_build_that_failed_its_checks_is_not_taken_for_done(self):
        """Every part the recipe makes there, its staging folder gone, the build VCFs there,
        and no mark: what a build that failed its checks (verify_mega_chrom.py, its last step)
        leaves, as a build of 2.10.0 run row by row leaves no mark either. The half took it for
        done and ended with exit 0; its line now runs the recipe, which refuses the working
        directory, and the half stops there, exit 1 (the PP-2 review, A1-PP2-1)."""
        rc, out, err = self.run_build(mega_parts="all")
        self.assertEqual(rc, 1, out + err)
        self.assertNotIn("the build of 20-NRR_3_hg38+hg38_mega: done before, skipped", out)
        self.assertIn("part(s) of hg38_mega, and nothing was made there", err)
        self.assertIn("rebuild: the build of 20-NRR_3_hg38+hg38_mega failed", err)

    def test_and_where_a_part_is_there_its_recipe_refuses_it(self):
        # (the recipe itself, which refuses the working directory before it makes anything)
        rc, out, err = self.run_build(mega_parts="some")
        self.assertEqual(rc, 1, out + err)
        self.assertIn("part(s) of hg38_mega, and nothing was made there", err)
        self.assertIn("rebuild: the build of 20-NRR_3_hg38+hg38_mega failed", err)


class ThePublishChecksItsTokenFirst(unittest.TestCase):
    """The kits' uploads, `$RUN hf upload ...`, failed as the README documented them, after the
    indexes' had gone up: the image sets HF_HUB_DISABLE_IMPLICIT_TOKEN=1, so its hf sends no
    token (only APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0 overrides it under apptainer, a host
    variable of that name does not), and from a login shell `$RUN hf` ran the host's
    ~/.local/bin/hf (the D7 readiness audit, R-PUB-1, G-P5; PP-1 of the release). The publish
    half checks first, before any upload: HF_TOKEN in the environment, exported for the image's
    hf with the override, not refused there, the hf it runs the one beside its crisprme.py, the
    token authenticating and the repo there. Here hf and crisprme.py are stubs on PATH, $RUN
    empty, so that this host stands for the image."""

    def setUp(self):
        self.rows = cis.read_manifest()

    def run_publish(self, token=True, hf_apart=False, whoami=0, info=0, panels=None,
                    trace=None, **env_extra):
        """The publish half run from the repository's root, of the rows of ``panels``'
        working directories alone where given (and only those set), under bash -x where
        ``trace`` is "-x", with SHELLOPTS=xtrace exported where it is "SHELLOPTS": (exit status,
        stderr, the stubs' log, one line a call: the program, then what the image's environment
        held, then its arguments)."""
        import shutil
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            image = os.path.join(t, "image", "bin")
            local = os.path.join(t, "home", ".local", "bin")
            os.makedirs(image)
            os.makedirs(local)
            log = os.path.join(t, "log")
            stub = ("#!/bin/sh\n"
                    "echo \"${0##*/} implicit=${APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN:-} "
                    "usersite=${PYTHONNOUSERSITE:-} $*\" >> %s\n"
                    "case \"$1 $2\" in \"auth whoami\") exit %d;; \"datasets info\") exit %d;; "
                    "esac\nexit 0\n" % (log, whoami, info))
            for folder, name in ((local if hf_apart else image, "hf"), (image, "crisprme.py")):
                with open(os.path.join(folder, name), "w") as fh:
                    fh.write(stub)
                os.chmod(os.path.join(folder, name), 0o755)
            env = {"PATH": os.pathsep.join([local, image, "/usr/bin", "/bin"]),
                   "HOME": os.path.join(t, "home"), "LC_ALL": "C", "RUN": ""}
            for v in panels or ("REF", "HGDP", "HPRC", "MEGA"):
                env[v] = os.path.join(t, v)
                os.makedirs(env[v])
            if token:
                env["HF_TOKEN"] = "test-not-a-token"
            if trace == "SHELLOPTS":
                env["SHELLOPTS"] = "xtrace"
            env.update(env_extra)
            script = cis.commands_script(self.rows, half="publish",
                                         **({"panels": panels} if panels else {}))
            got = subprocess.run([shutil.which("bash")] + (["-x"] if trace == "-x" else [])
                                 + ["-c", script],
                                 cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
            calls = []
            if os.path.exists(log):
                with open(log) as fh:
                    calls = fh.read().splitlines()
            return got.returncode, got.stderr, calls

    def uploads(self, calls):
        return [c for c in calls if " upload " in c or " publish-index " in c]

    def test_with_its_token_and_the_image_s_hf_it_checks_then_publishes(self):
        rc, err, calls = self.run_publish()
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual([c.split()[0] + " " + " ".join(c.split()[3:5]) for c in calls[:2]],
                         ["hf auth whoami", "hf datasets info"])
        self.assertIn("hf implicit=0 usersite=1 datasets info lucapinello/crisprme-indexes-format2",
                      calls[1])
        # every archive then, each command with the override exported for the image's hf
        self.assertEqual(len(self.uploads(calls)),
                         sum(r["kind"] in cis.INDEX_KINDS for r in self.rows) + 3 * 2)
        for c in self.uploads(calls):
            self.assertIn(" implicit=0 usersite=1 ", c)
        kits = [c for c in calls if c.startswith("hf ") and " upload " in c]
        self.assertTrue(kits)
        for c in kits:
            self.assertTrue(c.endswith(" --repo-type dataset"), c)

    def test_no_token_in_the_environment_publishes_nothing(self):
        """It stops with exit 2, the preflight's, as the README says: bash's own check of an
        unset variable (: "${HF_TOKEN:?...}") exited 1, a failed row's code (the PP-1 review,
        A1-PP1-1)."""
        for value in (None, ""):
            with self.subTest(HF_TOKEN=value):
                extra = {} if value is None else {"HF_TOKEN": value}
                rc, err, calls = self.run_publish(token=False, **extra)
                self.assertEqual(rc, 2, err)
                self.assertIn("publish: HF_TOKEN is not set: the token that publishes, a write "
                              "token of lucapinello/crisprme-indexes-format2, exported in the "
                              "environment: never on a command line", err)
                self.assertEqual(calls, [])

    def test_its_checks_never_write_the_token_into_a_trace(self):
        """Under bash -x, or with SHELLOPTS=xtrace exported, bash wrote '+ : <the token>', the
        argument of its check of the token, to stderr (the PP-1 review, A1-PP1-2): no line of
        the preflight expands the token, so a trace shows none (the hf of the image reads it
        from its environment)."""
        for trace in ("-x", "SHELLOPTS"):
            with self.subTest(trace=trace):
                rc, err, calls = self.run_publish(trace=trace)
                self.assertEqual(rc, 0, err)
                self.assertIn("+ export APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0", err)
                self.assertNotIn("test-not-a-token", err)
                self.assertTrue(self.uploads(calls))

    def test_each_node_publishes_the_rows_of_its_working_directories_alone(self):
        """The working directories are on two nodes, HGDP + 1000 Genomes' on one and the
        others on another, and each row reads its own: a half run once on either stopped at the
        first row of the other's, part of the set published (the PP review B-1). --panels
        prints the rows of the working directories named alone, which need only theirs set: run
        on each node, the two publish every archive of the set once."""
        published = []
        for panels in (["HGDP"], ["REF", "HPRC", "MEGA"]):
            with self.subTest(panels=panels):
                rc, err, calls = self.run_publish(panels=panels)
                self.assertEqual((rc, err), (0, ""))
                mine = [r for r in self.rows if r["publish"] != cis.NO_COMMAND
                        and cis.WORKDIR[r["panel"]] in panels]
                self.assertEqual(len(self.uploads(calls)),
                                 sum(1 if r["kind"] in cis.INDEX_KINDS else 2 for r in mine))
                for c in self.uploads(calls):
                    published.append(c)
                    self.assertTrue(any(f"/{p}/" in c for p in panels), c)
        indexes = [c.split(" --index ", 1)[1].split()[0] for c in published
                   if " publish-index " in c]
        self.assertEqual(sorted(i.rsplit("/", 1)[1] for i in indexes),
                         sorted(r["archive"] for r in self.rows if r["kind"] in cis.INDEX_KINDS))
        kits = [c.split()[-3] for c in published if " upload " in c]
        self.assertEqual(sorted(kits), sorted(
            f"{r['archive']}{ext}" for r in self.rows if r["kind"] == "build-inputs"
            for ext in (".SHA256SUMS", ".tar.gz")))

    def test_panels_are_named_by_their_working_directories(self):
        got, out, err = printed("--commands", "publish", "--panels", "HGDP")
        self.assertEqual(got, 0, err)
        self.assertIn("# The rows of $HGDP alone (--panels), for the node that holds them.", out)
        self.assertIn(": \"${HGDP:?", out)
        for other in ("REF", "HPRC", "MEGA"):
            self.assertNotIn(f"${other}", out)
        # (and $REF, read for its genome folder alone, which HPRC's is checked against)
        build = cis.commands_script(self.rows, half="build", panels=["HPRC"])
        self.assertEqual(
            [line for line in build.splitlines() if line.startswith(": \"${")],
            [line for line in build.splitlines() if line.startswith((": \"${HPRC:?",
                                                                     ": \"${RUN?", ": \"${SIF:?",
                                                                     ": \"${HPRC_SRC:?",
                                                                     ": \"${P1KG:?",
                                                                     f": \"${{REF:?{cis.REF_READ}"))])
        self.assertIn(f": \"${{REF:?{cis.REF_READ}}}\"", build)
        for args in (["--panels", "HGDP"], ["--commands", "publish", "--panels", "HGDP,FOO"],
                     ["--commands", "publish", "--panels", ","]):
            with self.subTest(args=args):
                got, out, err = printed(*args)
                self.assertEqual((got, out), (2, ""))
                self.assertIn("--panels", err)

    def test_an_image_that_would_send_no_token_publishes_nothing(self):
        for value in ("1", "true", "YES", "on"):
            with self.subTest(value=value):
                rc, err, calls = self.run_publish(HF_HUB_DISABLE_IMPLICIT_TOKEN=value)
                self.assertEqual(rc, 2, err)
                self.assertIn("publish: the release image would upload with no token", err)
                self.assertIn("docker run -e HF_TOKEN -e HF_HUB_DISABLE_IMPLICIT_TOKEN=0", err)
                self.assertEqual(self.uploads(calls), [])
        rc, err, calls = self.run_publish(HF_HUB_DISABLE_IMPLICIT_TOKEN="0")
        self.assertEqual((rc, err), (0, ""))

    def test_an_hf_that_is_not_the_image_s_own_publishes_nothing(self):
        rc, err, calls = self.run_publish(hf_apart=True)
        self.assertEqual(rc, 2, err)
        self.assertIn("/home/.local/bin/hf is the hf $RUN would run, not the release image's own, "
                      "installed beside its crisprme.py", err)
        self.assertIn("PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin bash", err)
        self.assertEqual(self.uploads(calls), [])

    def test_a_token_that_does_not_authenticate_or_no_repo_publishes_nothing(self):
        rc, err, calls = self.run_publish(whoami=1)
        self.assertEqual(rc, 2, err)
        self.assertIn("publish: HF_TOKEN does not authenticate on the Hub", err)
        self.assertEqual(self.uploads(calls), [])
        rc, err, calls = self.run_publish(info=1)
        self.assertEqual(rc, 2, err)
        self.assertIn("publish: lucapinello/crisprme-indexes-format2 is not on the Hub, or "
                      "HF_TOKEN cannot read it: its owner creates it", err)
        self.assertEqual(self.uploads(calls), [])

    def test_a_publish_half_that_does_not_check_first_is_found(self):
        from unittest import mock
        with mock.patch.object(cis, "_publish_checks", lambda: []):
            dry = cis.dry_run(self.rows)
        self.assertIn("the rebuild's publish half, run dry, publishes before it checks that the "
                      "token authenticates and the repo is there", dry)



def _read(*parts):
    with open(os.path.join(ROOT, *parts)) as fh:
        return fh.read()


def _flat(text):
    return " ".join(text.split())


class TheReleaseDocsOfTheSet(unittest.TestCase):
    """What the docs say of the publish of the set, and of the gates on it, holds (the PP review
    B, round 1): the publish runs on each node for the working directories it holds (B-1); the
    format-2 repo's card lists every file the publish uploads (B-2); the size and SHA-256 of an
    archive that publish-index packs at the publish are filled in after it, from the Hub (B-3);
    the timing gate names HPRC's runs and their baseline (B-4); the equivalence gate the release
    order names is a gate of the checklist (B-5); the archive gate reads each sub-index's window
    where it is recorded (B-6); and the rebuild's README points to sections the checklist has
    (B-7). And (round 2) the gates name what they compare with, and the release order fills in
    the figures known after the publish (B2-2 to B2-5)."""

    def setUp(self):
        self.rows = cis.read_manifest()
        self.readme = _read("seq_script", "rebuild", "README.md")
        self.checklist = _read("docs", "RELEASE_CHECKLIST.md")

    def gate(self, title):
        """The checklist's item whose bold title starts with ``title``, flattened."""
        text = _flat(self.checklist)
        start = text.index("**" + title)
        end = text.find("- [ ]", start)
        return text[start:end if end > 0 else None]

    def test_the_publish_half_runs_on_each_node_for_its_working_directories(self):
        """B-1: the working directories are on two nodes, and the docs ran the whole publish
        half once, which stopped at the first row of the other node's."""
        for name, text in (("README", self.readme), ("RELEASE_CHECKLIST", self.checklist)):
            calls = re.findall(r"check_index_set\.py --commands publish[^\n`]*", text)
            self.assertTrue(calls, name)
            for call in calls:
                self.assertIn(" --panels ", call, f"{name}: {call}")
        publish = self.readme[self.readme.index("**The publish.**"):]
        block = publish[publish.index("```bash\n"):]
        block = block[:block.index("\n```", 8)]
        named = [p for line in block.splitlines() if "--commands publish --panels " in line
                 for p in line.split("--panels ", 1)[1].split()[0].split(",")]
        self.assertEqual(sorted(named), sorted(set(cis.WORKDIR.values())))
        # the gates on ml008 read the indexes of HGDP + 1000 Genomes, built on ml007
        gate = self.gate("The real-data checker on the rebuilt indexes")
        self.assertIn("HGDP + 1000 Genomes' working directory, `$HGDP`, is on ml007", gate)
        self.assertIn("copy to ml008, into a folder of its own, what the install links of it, "
                      "from ml007 with `rsync -a`", gate)
        self.assertIn("The runs of HGDP + 1000 Genomes are made on ml008 on the copy of its "
                      "indexes the checker's gate makes", self.gate("The timing gate"))

    def test_the_card_lists_every_file_the_publish_uploads(self):
        """B-2: the publish uploads each kit's member list and HPRC's normalised VCFs (E, step
        4), which the card's layout did not list: a file on the Hub the card does not document
        fails the card's own gate."""
        card = _read("docs", "hf", "crisprme-indexes-format2_README.md")
        layout = card[card.index("## Layout"):]
        layout = layout[layout.index("```"):layout.index("```", layout.index("```") + 3)]
        for r in self.rows:
            if r["kind"] != "build-inputs":
                continue
            ds = r["archive"].split("/", 1)[1]
            line = next(ln for ln in layout.splitlines() if f"── {ds}.tar.gz" in ln)
            self.assertIn(".SHA256SUMS", line)
            self.assertTrue(".members.md5" in line or f"{ds}.members.md5" in layout, ds)
        # (what the release order uploads: E, step 4)
        self.assertIn("build-inputs/hg38_HGDP1kGP.members.md5 --repo-type dataset` (and HPRC's "
                      "and mega's alike)", _flat(self.checklist))
        hprc = layout[layout.index("hg38_HPRC/"):]
        self.assertIn("hprc.chr1.norm.vcf.gz (+ .tbi)", hprc)
        self.assertIn("SHA256SUMS", hprc)
        self.assertIn("build-inputs/hg38_HPRC --repo-type dataset", _flat(self.checklist))
        verifying = _flat(card[card.index("## Verifying"):])
        self.assertIn("`md5sum -c <dataset>.members.md5`", verifying)
        self.assertIn("`build-inputs/hg38_HPRC/`, `sha256sum -c --ignore-missing SHA256SUMS`",
                      verifying)

    def test_each_archive_s_size_and_sha256_are_filled_in_where_they_exist(self):
        """B-3: publish-index packs each index and companion archive at the publish, its
        manifest dated then, and deletes it once uploaded: their size and SHA-256 are the Hub's,
        known after the publish; a kit's, which its build row writes, after the build."""
        prov = _read("docs", "INDEX_PROVENANCE.md")
        section = prov[prov.index("### Each archive"):
                       prov.index("## The archives of CRISPRme 2.8.1")]
        header = next(ln for ln in section.splitlines() if ln.startswith("| Archive |"))
        cols = [c.strip() for c in header.strip().strip("|").split("|")]
        size, sha = cols.index("Size"), cols.index("SHA-256")
        seen = set()
        for line in section.splitlines():
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if not line.startswith("| `") or len(cells) != len(cols):
                continue  # (the kits' member lists have a table of their own)
            names = re.findall(r"`([^`]+)`", cells[0])
            want = "*after the build*" if names[0].startswith("build-inputs/") \
                else "*after the publish*"
            self.assertEqual((cells[size], cells[sha]), (want, want), line)
            seen.update(names)
        self.assertEqual(seen, {r["archive"] for r in self.rows})
        self.assertIn("after the publish, the size and SHA-256 of each of the 23 index and "
                      "companion archives, from the Hub",
                      self.gate("Published-index provenance recorded"))
        self.assertEqual(23, sum(r["kind"] != "build-inputs" for r in self.rows))

    def test_the_timing_gate_names_hprc_s_runs_and_their_baseline(self):
        """B-4: gate38a, the timing gate's baseline, has no run of HPRC."""
        gate = self.gate("The timing gate")
        self.assertIn("gate38a has no run of HPRC, so its configurations are timed now, both "
                      "sides", gate)
        self.assertIn("the README guide with NRR, per sample and population-level, on the "
                      "rebuilt `20-NRR_3_hg38+hg38_HPRC`", gate)
        self.assertIn("on 2.8.1's `NRR_3_hg38+hg38_HPRC` (format 1) with 2.8.1's image", gate)
        self.assertIn("NRR_3_hg38+hg38_HPRC", cis.ARCHIVES_281)
        self.assertIn("20-NRR_3_hg38+hg38_HPRC", {r["archive"] for r in self.rows})

    def test_the_equivalence_gate_the_release_order_names_is_one_of_the_checklist(self):
        """B-5: the release order's step 1 named "the G5 equivalence (D)", which no section of
        the checklist defined (it was the readiness audit's section D)."""
        text = _flat(self.checklist)
        order = text[text.index("### E. The order of 2.10.0's release"):]
        order = order[:order.index("2. **Luca's go #1**")]
        self.assertIn("the equivalence gate, G5 (D)", order)
        d = text[text.index("### D. Figures the release measures"):
                 text.index("### E. The order of 2.10.0's release")]
        self.assertIn("- [ ] **The equivalence gate (G5)**", d)
        g5 = d[d.index("**The equivalence gate (G5)**"):]
        g5 = g5[:g5.index("- [ ]")]
        for same in ("`20-NRG_3_hg38+hg38_HGDP1kGP` = `20-NNNNNNNN_3_hg38+hg38_HGDP1kGP` pruned "
                     "to NRG = `20-NRR_3_hg38+hg38_HGDP1kGP` pruned to NRG",
                     "`20-NRR_3_hg38+hg38_HGDP1kGP` = `20-NNNNNNNN_3_hg38+hg38_HGDP1kGP` pruned "
                     "to NRR",
                     "`20-NGG_3_hg38+hg38_HGDP1kGP` = `20-NRR_3_hg38+hg38_HGDP1kGP` pruned to NGG"):
            self.assertIn(same, g5)
        archives = {r["archive"] for r in self.rows}
        for name in re.findall(r"`(\d+-[A-Z]+_3_hg38\+hg38_HGDP1kGP)`", g5):
            self.assertIn(name, archives)

    def test_the_archive_gate_reads_each_sub_index_s_window_in_its_own_record(self):
        """B-6: the manifest gives one indel_window, the panel's logs', and no sub-index's:
        each sub-index records its own in its .index_meta.json, which the archive gate reads."""
        gate = self.gate("Every published variant archive can find indels")
        self.assertNotIn("`indel_window` [30, 30] for each sub-index", gate)
        self.assertIn("`indel_window` [30, 30] (one value, the window of the panel's indel logs",
                      gate)
        self.assertIn("in `index_meta` (the variant index's own record, which holds no window), "
                      "`built_by` {crispritz 2.10.0, crisprme 2.10.0}", gate)
        self.assertIn("each sub-index's own record in the stream, "
                      "`<index>_INDELS/<sub-index>/.index_meta.json`, says `indel_window` [30, 30]",
                      gate)
        hf = _read("PostProcess", "crisprme_hf.py")
        self.assertIn('manifest["indel_window"] = list(_w)', hf)
        self.assertIn("indel_window=_log_window(chrom)",
                      _read("PostProcess", "pool_index_indels.py"))

    def test_the_gates_name_what_they_compare_with_and_what_is_filled_in_after(self):
        """The PP review B, round 2: G5 diffed "each" of its runs against gate38a's runs of
        2.8.1, which hold none of NGG (B2-2); the 2.8.1 side of HPRC's timing run had no NRR PAM
        file, which neither 2.8.1's tree nor the data repo's pams/ holds before the publish
        (B2-3); no step of the release order filled in the provenance figures known only after
        the publish (B2-4); and docs/hf/README.md said that the repo holds the set's archives
        "and nothing else", not the files the publish puts beside them (B2-5)."""
        g5 = self.gate("The equivalence gate (G5)")
        self.assertNotIn("Then diff each against gate38a's runs", g5)
        self.assertIn("Then diff the NRG and NRR runs against gate38a's runs of 2.8.1 of the same "
                      "PAM", g5)
        self.assertIn("gate38a has no run of NGG", g5)
        self.assertIn("its search given 2.10.0's NRR PAM file, `PAMs/20bp-NRR-SpCas9.txt` "
                      "(`NNNNNNNNNNNNNNNNNNNNNRR 3`)", self.gate("The timing gate"))
        self.assertEqual(_read("PAMs", "20bp-NRR-SpCas9.txt").strip(),
                         "NNNNNNNNNNNNNNNNNNNNNRR 3")
        text = _flat(self.checklist)
        step5 = text[text.index("5. **`indexes/SHA256SUMS`**"):
                     text.index("6. **The archive gates, from the Hub**")]
        self.assertIn("Then fill in `docs/INDEX_PROVENANCE.md`'s size and SHA-256 of the 16 "
                      "indexes and 7 companions", step5)
        self.assertIn("before the tags (Gate A)", step5)
        self.assertIn("what the publish puts beside them (`indexes/SHA256SUMS`, each kit's "
                      "`.SHA256SUMS` and member list, `build-inputs/<dataset>.members.md5`, and "
                      "HPRC's normalised VCFs, `build-inputs/hg38_HPRC/`), and nothing else",
                      _flat(_read("docs", "hf", "README.md")))

    def test_the_readme_points_to_sections_the_checklist_has(self):
        """B-7: the README sent the uploads it does not make to "Gate G", which the checklist
        does not have (the readiness audit's Hugging Face section was)."""
        sections = re.findall(r"^#{2,3} ([A-Z0-9]+)\. ", self.checklist, re.M)
        refs = re.findall(r"`docs/RELEASE_CHECKLIST\.md`,\s+(?:Gate\s+)?([A-Z0-9]+)\b",
                          self.readme)
        self.assertTrue(refs)
        for ref in refs:
            self.assertIn(ref, sections)
        self.assertIn("(`docs/RELEASE_CHECKLIST.md`, E, step 4)", _flat(self.readme))


if __name__ == "__main__":
    unittest.main()

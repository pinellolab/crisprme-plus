#!/usr/bin/env python3
"""check_index_set.py - check the index set of CRISPRme 2.10.0 (the rebuild manifest) against the
code that builds, names, publishes and downloads its archives, and run its commands dry.

The manifest, ``seq_script/rebuild/index_set_2.10.0.tsv``, lists every archive the rebuild for
2.10.0 builds and publishes in ``lucapinello/crisprme-indexes-format2`` (decision C2), with its
build and publish commands. The rebuild follows it, and the D7 readiness audit runs this check
before any build starts: every index's name must be the one ``build-index-only`` gives its PAM
file (``index_fit.index_name``: guide length and side, E2), every archive must route to the
format-2 repo (``crisprme_hf.resolve_index_repo``), every index must serve its own guide at the
web's ceiling of 2 DNA + 2 RNA bulges (``index_fit.fits``), every variant panel must have its
Dictionaries bundle (and, genotyped, its genotype stores) once, and every index archive of 2.8.1
must be replaced by exactly one archive of the set (NNN_3, not built again, by the general
index), as a download of it says (``crisprme_hf.successor_index``).

Its commands are what the rebuild runs, and ``--commands build`` and ``--commands publish`` print
them as the two halves of the rebuild, a script each (HALVES): the build half publishes nothing,
and the publish half, which builds nothing, is run once the indexes have passed the release's
gates and the publish is approved, on each node for the working directories it holds
(``--panels``: the rows of those alone; a part of the build half checks their genome folders as
the whole does, against $REF's, REF_READ). Each command a bash command (``bash -n``; a note has
a column of its own) run from the repository's root on its own, so that none moves the next (a
``cd`` only inside ``( )``); a recipe called by its repository path, and given every input it
would otherwise read from a path of the cluster it was written on; each build first copying its
PAM file from the repository's ``PAMs/`` into its working directory's, and reading it there;
``crisprme.py``, ``hf`` and ``pigz`` run in the release image (``$RUN``); each command in its own
panel's working directory, with the variables of VARIABLES alone; a panel companion, which its
panel's first build writes and first publish uploads, with none ("-"). Each half skips a row it
ran to its end before (its mark in the working directory, MARKS), and the build half a kit whose
SHA256SUMS holds, written once its tarball is whole; nothing else, as what a recipe leaves does
not tell that it ran to its end (mega's writes every part of its dataset before its checks). The
check reads every command so, then runs the two halves dry, then each again (dry_run): every
program they call stubbed, in working directories that hold what the builds leave, from a copy
of the repository's root made of links. Each runs to its end only where every script it calls is
there, every PAM file a build reads is in place, and every file a publish uploads is where the
command says; run again, each does nothing again.

    python seq_script/rebuild/check_index_set.py                       # the check and the dry run; exit 1 on a problem
    python seq_script/rebuild/check_index_set.py --commands build      # the build half: a script bash runs from the repository's root
    python seq_script/rebuild/check_index_set.py --commands publish    # the publish half, after the gates and the go
    python seq_script/rebuild/check_index_set.py --commands publish --panels HGDP   # its rows of $HGDP alone (a node's)

Standard library and the repository's PostProcess modules only; no network.
"""

import argparse
import csv
import functools
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MANIFEST = os.path.join(HERE, "index_set_2.10.0.tsv")
sys.path.insert(0, os.path.join(ROOT, "PostProcess"))

import crisprme_hf  # noqa: E402
import index_fit  # noqa: E402

COLUMNS = ("archive", "kind", "panel", "pam_file", "pam", "side", "spacer", "bmax", "bdna",
           "brna", "repo", "replaces", "build", "publish", "notes")
INDEX_KINDS = ("reference", "variant")
# The panel companions, which the panel's first build writes (Dictionaries/) and its first
# publish-index uploads: no command of their own ("-" in both columns)
COMPANION_KINDS = ("bundle", "genotypes", "indel-genotypes")
NO_COMMAND = "-"
# The index archives of 2.8.1 (index format 1), as `scripts/hf_sha256sums.py` listed the repos
# of 2.8.1 on 2026-10-30: lucapinello/crisprme-data and lucapinello/crisprme-data-general-index.
ARCHIVES_281 = (
    "NGG_3_hg38", "NGG_3_hg38+hg38_HGDP1kGP",
    "NNN_3_hg38", "NNN_3_hg38+hg38_HGDP1kGP",
    "NRG_3_hg38", "NRG_3_hg38+hg38_HGDP1kGP", "NRG_3_hg38+hg38_HPRC", "NRG_3_hg38+hg38_mega",
    "NRR_3_hg38", "NRR_3_hg38+hg38_HGDP1kGP", "NRR_3_hg38+hg38_HPRC", "NRR_3_hg38+hg38_mega",
    "NNNNNNNN_3_hg38", "NNNNNNNN_3_hg38+hg38_HGDP1kGP",
)
# The panels whose stores hold genotypes (mega is sites-only), and the enriched-genome kits
# of 2.8.1 in lucapinello/crisprme-data
GENOTYPED = ("hg38_HGDP1kGP", "hg38_HPRC")
KITS_281 = ("build-inputs/hg38_HGDP1kGP", "build-inputs/hg38_mega")
# The ceiling of the web's Advanced panel: an index of the set must serve it for its own guide
WEB_BULGES = 2

# The variables the commands use, each what the rebuild sets it to (seq_script/rebuild/README.md);
# --commands' script stops before its first command where one is not set. No text here may hold
# a '}', a '$', a '`', a '"' or a "'": each goes into the script's ${NAME:?text} (_VARIABLE_TEXT).
VARIABLES = {
    "REF": "the working directory of the reference indexes, holding all of hg38 in Genomes/hg38, "
           "as crisprme.py download --what genome installs it",
    "HGDP": "the working directory of hg38_HGDP1kGP, holding all of hg38 in Genomes/hg38, as the "
            "reference working directory does",
    "HPRC": "the working directory of hg38_HPRC, holding the chromosomes chr1-22, X and Y of "
            "hg38 in Genomes/hg38, one FASTA each, and no other",
    "MEGA": "the working directory of hg38_mega, holding the chromosomes chr1-22, X and Y of "
            "hg38 in Genomes/hg38, one FASTA each, and no other",
    "RUN": "the prefix that runs a command in the release image, empty where the release is "
           "installed: apptainer exec --bind <the folders> <the image> "
           "/usr/local/bin/_entrypoint.sh",
    "SIF": "the release image, which the HGDP + 1000 Genomes and HPRC recipes run",
    "RP": "the folder of the re-phased HGDP + 1000 Genomes BCFs, hgdp1kg.<chr>.rephased.bcf",
    "META": "the metadata file of HGDP + 1000 Genomes by gnomAD, gnomad_meta_updated.tsv",
    "HPRC_SRC": "the folder of the HPRC VCFs, ALL.<chr>.HPRC.vcfwave.vcf.gz",
    "P1KG": "the samplesID file of 1000 Genomes 2021, which names the populations of the HPRC "
            "samples",
    "MEGA_SRC": "the folder of the merged VCFs of mega, mega.<chr>.afmax.vcf.gz",
}
# The ones that may be set empty (RUN: the release installed where the script runs)
MAY_BE_EMPTY = ("RUN",)
# What $REF is to a part of the build half (--panels) that names the rows of HGDP + 1000
# Genomes', HPRC's or mega's working directory and not the reference's: the genome folder theirs
# are checked against (README, rule 1), read and not written; on a node without the reference
# working directory, a folder that holds hg38 in Genomes/hg38 as it does, of its own: a copy, or
# links to the FASTAs of another folder of hg38, never to those of a folder it checks. A part
# that named HGDP + 1000 Genomes' alone checked only that its folder held a FASTA, and one of
# HPRC's or mega's only their names (the PP review, round 2: A2-PP2-1, B2-1); and the README then
# offered a $REF of links to $HGDP's own FASTAs, as the D7 plan's stub on ml007 is, against which
# $HGDP's folder was checked against itself: a FASTA missing there was missing from $REF's too,
# and the check passed (the PP review, round 3: B3-1). A text of VARIABLES' kind (_VARIABLE_TEXT).
REF_READ = ("a folder whose Genomes/hg38 holds all of hg38 of its own, as the working directory "
            "of the reference indexes does: that one, a copy, or links to the FASTAs of another "
            "such folder, never to those of a genome folder checked against it, which would then "
            "be checked against itself. The genome folders of the rows here are checked against "
            "it, and nothing is written there")
# What each working directory's genome folder, Genomes/hg38, holds (README, rule 1), as 2.8.1's
# default index of each panel was built (the contigs of the published archives): all of hg38
# (455 contigs, as `download --what genome` installs it) for the reference indexes and for HGDP
# + 1000 Genomes (2.8.1's NRR, NGG and NNN indexes of the panel; its NRG index held chr1-22, X
# and Y), its chromosomes chr1-22, X and Y alone (PRIMARY) for HPRC and mega (both of 2.8.1's
# indexes of each); mega's recipe builds every chromosome of its merged VCFs,
# mega.<chr>.afmax.vcf.gz (those 24 in vcfs/hg38_mega), and refuses a genome folder that holds
# any other FASTA. "REF": the reference working directory's FASTAs; "PRIMARY": those of them
# that are chr1-22, X and Y. The rebuild's script checks the genome folders of the working
# directories it uses before its first command (commands_script).
PRIMARY = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
_PRIMARY_RE = "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"
GENOME_FOLDERS = {"HGDP": "REF", "HPRC": "PRIMARY", "MEGA": "PRIMARY"}
# A panel's recipe that builds every chromosome of its sources: the variable of their folder,
# and how a file of one chromosome is named there (<prefix><chr><suffix>)
SOURCE_CHROMOSOMES = {"MEGA": ("MEGA_SRC", "mega.", ".afmax.vcf.gz")}
# What a text of VARIABLES may not hold
_VARIABLE_TEXT = re.compile(r"[}$`\"']")
# Each panel's working directory, the variable that names it
WORKDIR = {"-": "REF", "hg38_HGDP1kGP": "HGDP", "hg38_HPRC": "HPRC", "hg38_mega": "MEGA"}
# The programs of the release image, which a command runs through $RUN
IMAGE_PROGRAMS = ("crisprme.py", "hf", "huggingface-cli", "pigz", "python", "python3")
# The programs that upload to the Hub, whose repo a command names (hf upload <repo> ...)
UPLOADERS = ("hf", "huggingface-cli")
# The repository's folders a command names a path in, read from the repository's root (PAMs/
# names the working directory's after a cd into it)
REPOSITORY_FOLDERS = ("seq_script", "scripts", "PostProcess", "PAMs")
# A recipe the set calls that a branch not merged yet brings, which this tree does not hold
# yet: the merge brings it, and its entry here goes then (the check says so once the tree
# holds it). Each is said, not counted as a problem, while it is missing. None now: mega's
# recipe, seq_script/merge_panels/mega_rebuild.sh (lane M's, decision D2), came with the
# integration of stage 2b, and the dry run checks its call as it checks the others'.
PENDING_SCRIPTS = {}
# The two halves of the rebuild, each a script of its own (--commands build, --commands publish):
# the build half publishes nothing, and the publish half is run after the gates and the go
# (README, "Commands"; PP-2 of the release: the script of both halves published straight after
# its builds, before any gate)
HALVES = ("build", "publish")
# The folder of each working directory where a half of the rebuild marks a row it ran to its end
# (<half>.<archive>, a "/" as "_"), so that a re-run skips it (README, "A re-run"). A dot folder,
# which no archive packs: an index's is its genome_library folder, a kit's its Genomes/ folders.
MARKS = ".rebuild-2.10.0"


def read_manifest(path=MANIFEST):
    """The manifest's rows, each a dict of COLUMNS (comment lines, starting with '#', skipped)."""
    with open(path, newline="") as fh:
        lines = [line for line in fh if line.strip() and not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)
    if tuple(reader.fieldnames or ()) != COLUMNS:
        raise ValueError(f"{path}: the columns are {reader.fieldnames}, not {list(COLUMNS)}")
    return [dict(row) for row in reader]


def read_pam_file(path):
    """(PAM, side "3'" or "5'", guide length) of a PAM file, as build-index-only reads it."""
    with open(path) as fh:
        seq, pos = fh.readline().split()
    pos = int(pos)
    pam = seq[:-pos] if pos < 0 else seq[-pos:]
    return pam, ("5'" if pos < 0 else "3'"), len(seq) - abs(pos)


@functools.lru_cache(maxsize=None)
def bash_syntax_error(command):
    """What ``bash -n`` says of ``command`` (its syntax only: nothing runs), "" when bash takes
    it."""
    got = subprocess.run(["bash", "-n"], input=command, text=True, capture_output=True)
    return "" if got.returncode == 0 else (got.stderr.strip() or f"exit {got.returncode}")


def variables(command):
    """The names of the variables ``command`` reads ($NAME or ${NAME}), in order, once each."""
    seen = []
    for name in re.findall(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)", command):
        if name not in seen:
            seen.append(name)
    return seen


def simple_commands(command):
    """The simple commands of ``command``, in order: (assignments, words, depth, moved) each,
    its leading NAME=value words, the rest (a redirection's target left out), how many ( ) it is
    inside, and whether a cd before it, in its ( ) or one around it, moved its directory."""
    lexer = shlex.shlex(command, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    found, moved, assignments, words = [], [False], [], []
    redirect = False

    def end():
        if assignments or words:
            found.append((list(assignments), list(words), len(moved) - 1, moved[-1]))
            if words and words[0] == "cd":
                moved[-1] = True
        assignments.clear()
        words.clear()

    for token in lexer:
        if redirect:
            redirect = False
        elif token == "(":
            end()
            moved.append(moved[-1])
        elif token == ")":
            end()
            if len(moved) > 1:
                moved.pop()
        elif token in ("&&", "||", "|", ";", "&", ";;", "|&"):
            end()
        elif token[0] in "<>":
            redirect = True
        elif not words and re.match(r"[A-Za-z_][A-Za-z0-9_]*=", token):
            assignments.append(token)
        else:
            words.append(token)
    end()
    return found


@functools.lru_cache(maxsize=None)
def script_inputs(path):
    """The inputs recipe ``path`` reads from a path it was written with unless it is given
    them: (NAME, default) for each NAME=${NAME:-<default>} whose NAME it reads after, its
    default an absolute path (of the cluster it was written on), a command that names one, or
    a path under another of its variables ($B, the working directory: HGDP + 1000 Genomes'
    RP=${RP:-$B} reads the re-phased BCFs from there). A default that is a command of no
    absolute path (mega's PP, the repository's PostProcess/), a relative path or no path is not
    one."""
    with open(path) as fh:
        text = fh.read()
    found = []
    for m in re.finditer(r"^\s*([A-Za-z_][A-Za-z0-9_]*)=\$\{\1:-(/[^}]*|\$\([^)]*\s/[^)]*\)"
                         r"|\$\{?[A-Za-z_][A-Za-z0-9_]*\}?(?:/[^}]*)?)\}", text, re.M):
        name = m.group(1)
        if re.search(r"\$\{?%s\b" % name, text[m.end():]):
            found.append((name, m.group(2)))
    return tuple(found)


def check_commands(r, root=ROOT):
    """The problems of row ``r``'s build and publish commands, read as the rebuild runs them."""
    problems = []
    a, kind = r["archive"], r["kind"]
    own = WORKDIR.get(r["panel"])
    for what in ("build", "publish"):
        command = r[what]
        if kind in COMPANION_KINDS:
            if command != NO_COMMAND:
                problems.append(f"{a}: a companion has no {what} command of its own (its "
                                f"panel's first build writes it, its first publish uploads "
                                f"it): {NO_COMMAND!r}, not {command!r}")
            continue
        if command in ("", NO_COMMAND):
            problems.append(f"{a}: no {what} command")
            continue
        error = bash_syntax_error(command)
        if error:
            problems.append(f"{a}: its {what} command is not one bash can run ({error}); a "
                            f"note goes in the notes column")
            continue
        for name in variables(command):
            if name not in VARIABLES:
                problems.append(f"{a}: its {what} command reads ${name}, which the rebuild does "
                                f"not set (check_index_set.VARIABLES)")
        dirs = {name for name in variables(command) if name in WORKDIR.values()}
        if own is None or own not in dirs:
            problems.append(f"{a}: its {what} command does not name its working directory "
                            f"(${own})")
        for other in sorted(dirs - {own}):
            problems.append(f"{a}: its {what} command uses ${other}, the working directory of "
                            f"another panel: a panel's indexes share its working directory's "
                            f"Dictionaries alone")
        for assignments, words, depth, moved in simple_commands(command):
            program, in_image = (words[:1] or [""])[0], False
            if program == "$RUN":
                program, in_image = (words[1:2] or [""])[0], True
                if not program:
                    problems.append(f"{a}: its {what} command runs nothing through $RUN")
            if program in IMAGE_PROGRAMS and not in_image:
                problems.append(f"{a}: its {what} command runs {program} outside the release "
                                f"image: $RUN {program}")
            rest = words[2:] if in_image else words[1:]
            publishes = ((program in UPLOADERS and rest[:1] == ["upload"])
                         or (program == "crisprme.py" and rest[:1] == ["publish-index"]))
            builds = ((program == "crisprme.py" and rest[:1] == ["build-index-only"])
                      or program in ("bash", "sh"))
            if what == "build" and publishes:
                problems.append(f"{a}: its build command publishes ({program} {rest[0]}): the "
                                f"build half of the rebuild publishes nothing, its publish half "
                                f"runs once, after the gates and the go (--commands publish)")
            if what == "publish" and builds:
                problems.append(f"{a}: its publish command builds ({program}): the publish half "
                                f"of the rebuild runs once the build half has built the set")
            if program in UPLOADERS and rest[:1] == ["upload"]:
                repo = (rest[1:2] or [""])[0]
                if repo != crisprme_hf.FORMAT2_HF_REPO:
                    problems.append(f"{a}: its {what} command uploads to {repo or 'no repo'}, "
                                    f"not {crisprme_hf.FORMAT2_HF_REPO}: the set goes to the "
                                    f"repo of index format 2 alone (decision C2), and an "
                                    f"upload elsewhere replaces what 2.8.1 reads")
            if "--hf-repo" in rest:
                repo = (rest[rest.index("--hf-repo") + 1:][:1] or [""])[0]
                if repo != crisprme_hf.FORMAT2_HF_REPO:
                    problems.append(f"{a}: its {what} command names the repo {repo or '(none)'}"
                                    f", not {crisprme_hf.FORMAT2_HF_REPO}")
            if "--token" in rest:
                problems.append(f"{a}: its {what} command passes --token: a token goes in "
                                f"HF_TOKEN, never in the manifest")
            if program == "cd" and depth == 0:
                problems.append(f"{a}: its {what} command moves the directory of every command "
                                f"after it: its cd goes inside ( ), the commands running from "
                                f"the repository's root")
            for word in words:
                folder = word.split("/", 1)[0]
                if folder not in REPOSITORY_FOLDERS or "/" not in word:
                    continue
                if moved:
                    if folder != "PAMs":
                        problems.append(f"{a}: its {what} command names the repository's "
                                        f"{word} after a cd, where it is not")
                elif word in PENDING_SCRIPTS:
                    pass
                elif not os.path.exists(os.path.join(root, word)):
                    problems.append(f"{a}: its {what} command names {word}, which the "
                                    f"repository does not have")
            if program in ("bash", "sh"):
                script = next((w for w in words[1:] if not w.startswith("-")), "")
                if script.split("/", 1)[0] not in ("seq_script", "scripts"):
                    problems.append(f"{a}: its {what} command runs {script or 'no script'}, not "
                                    f"a recipe of the repository's seq_script/ or scripts/")
                elif os.path.isfile(os.path.join(root, script)) and not moved:
                    given = {x.split("=", 1)[0] for x in assignments}
                    for name, default in script_inputs(os.path.join(root, script)):
                        if name not in given:
                            problems.append(f"{a}: {script} reads {name} from {default} unless "
                                            f"it is given it: its {what} command gives "
                                            f"{name}=")
    if kind in INDEX_KINDS and r["build"] not in ("", NO_COMMAND):
        base = os.path.basename(r["pam_file"])
        copy = f"mkdir -p ${own}/PAMs && cp {r['pam_file']} ${own}/PAMs/ && "
        if not r["build"].startswith(copy):
            problems.append(f"{a}: its build command does not first copy its PAM file into its "
                            f"working directory, which the setup of a working directory does not "
                            f"hold: it starts {copy!r}")
        elif f"PAMs/{base}" not in r["build"][len(copy):]:
            problems.append(f"{a}: its build command does not read the PAM file it copies, "
                            f"PAMs/{base}")
    if kind in INDEX_KINDS and f"--index ${own}/genome_library/{a}" not in r["publish"]:
        problems.append(f"{a}: its publish command does not publish ${own}/genome_library/{a}")
    if kind == "variant" and "--dictless" not in r["publish"].split():
        problems.append(f"{a}: its publish command does not pass --dictless: the panel's "
                        f"Dictionaries bundle would carry its per-sample dictionaries "
                        f"(dictionaries_{r['panel']}), which the other indexes of the panel do "
                        f"not name, and differ in fingerprint from theirs")
    return problems


def check(rows, root=ROOT):
    """The problems of the index set ``rows`` (a list of strings; none: it holds)."""
    problems = []
    for r in rows:
        problems += check_commands(r, root)
    names = [r["archive"] for r in rows]
    for dup in sorted({n for n in names if names.count(n) > 1}):
        problems.append(f"{dup}: listed {names.count(dup)} times")
    texts = [(f"VARIABLES[{name!r}]", name, text) for name, text in VARIABLES.items()]
    for what, name, text in texts + [("REF_READ", "REF", REF_READ)]:
        if _VARIABLE_TEXT.search(text):
            problems.append(f"{what} holds one of }} $ ` \" ', which the rebuild's script cannot "
                            f"hold in ${{{name}:?...}}")
    read = {name for r in rows for what in ("build", "publish") if r[what] != NO_COMMAND
            for name in variables(r[what])}
    for name in VARIABLES:
        if name not in read:
            problems.append(f"VARIABLES[{name!r}]: no command reads ${name}: a command lost an "
                            f"input it needs, or the entry goes")
    for script in PENDING_SCRIPTS:
        if os.path.exists(os.path.join(root, script)):
            problems.append(f"{script}: this tree holds it, so it is no longer pending: take it "
                            f"out of check_index_set.PENDING_SCRIPTS")
    panels_with_index = set()
    replaced = {}
    for r in rows:
        a, kind = r["archive"], r["kind"]
        if r["repo"] != crisprme_hf.FORMAT2_HF_REPO:
            problems.append(f"{a}: repo {r['repo']}, not {crisprme_hf.FORMAT2_HF_REPO}")
        for old in [x for x in r["replaces"].split(",") if x and x != "-"]:
            replaced.setdefault(old, []).append(a)
        if kind in INDEX_KINDS:
            pam_file = os.path.join(root, r["pam_file"])
            if not os.path.isfile(pam_file):
                problems.append(f"{a}: no PAM file {r['pam_file']}")
                continue
            pam, side, spacer = read_pam_file(pam_file)
            if (pam, side, spacer) != (r["pam"], r["side"], int(r["spacer"])):
                problems.append(f"{a}: {r['pam_file']} holds {pam} {side} {spacer} nt, the row "
                                f"says {r['pam']} {r['side']} {r['spacer']} nt")
            bmax, bdna, brna = int(r["bmax"]), int(r["bdna"]), int(r["brna"])
            if bmax != max(bdna, brna) + 1:
                problems.append(f"{a}: bmax {bmax} is not max(bDNA {bdna}, bRNA {brna}) + 1")
            rest = "hg38" if kind == "reference" else f"hg38+{r['panel']}"
            want = index_fit.index_name(pam, max(bdna, brna) + 1, rest, spacer, side == "3'")
            if a != want:
                problems.append(f"{a}: build-index-only names this build {want}")
            if not crisprme_hf.is_named_for_guide(a):
                problems.append(f"{a}: not a name of the form of 2.10.0")
            if crisprme_hf.resolve_index_repo(a) != crisprme_hf.FORMAT2_HF_REPO:
                problems.append(f"{a}: a download routes it to "
                                f"{crisprme_hf.resolve_index_repo(a)}")
            if not index_fit.fits(spacer, bmax, spacer, WEB_BULGES, WEB_BULGES):
                problems.append(f"{a}: does not serve its own {spacer}-nt guide at "
                                f"{WEB_BULGES} DNA + {WEB_BULGES} RNA bulges")
            if kind == "variant":
                panels_with_index.add(r["panel"])
                if not a.endswith(f"_hg38+{r['panel']}"):
                    problems.append(f"{a}: not an index of panel {r['panel']}")
                base = a.partition("+")[0]
                if base not in names:
                    problems.append(f"{a}: its reference index {base}, which a download of it "
                                    f"installs, is not in the set")
            if kind == "reference" and f"--bDNA {bdna} --bRNA {brna}" not in r["build"]:
                problems.append(f"{a}: its build command does not give --bDNA {bdna} --bRNA {brna}")
        elif kind in ("bundle", "genotypes", "indel-genotypes"):
            prefix = {"bundle": crisprme_hf.BUNDLE_PREFIX, "genotypes": "genotypes_",
                      "indel-genotypes": "indel_genotypes_"}[kind]
            if a != prefix + r["panel"]:
                problems.append(f"{a}: the {kind} of {r['panel']} is {prefix}{r['panel']}")
            if crisprme_hf._companion_kind(a) != prefix:
                problems.append(f"{a}: the catalog would offer it as an index")
        elif kind == "build-inputs":
            if a != f"build-inputs/{r['panel']}":
                problems.append(f"{a}: the kit of {r['panel']} is build-inputs/{r['panel']}")
        else:
            problems.append(f"{a}: unknown kind {kind}")
    by_kind = {}
    for r in rows:
        by_kind.setdefault((r["kind"], r["panel"]), []).append(r["archive"])
    for panel in sorted(panels_with_index):
        want = ["bundle", "build-inputs"] + (["genotypes", "indel-genotypes"]
                                             if panel in GENOTYPED else [])
        for kind in want:
            got = by_kind.get((kind, panel), [])
            if len(got) != 1:
                problems.append(f"panel {panel}: {len(got)} {kind} archive(s), not 1")
        if panel not in GENOTYPED and by_kind.get(("genotypes", panel)):
            problems.append(f"panel {panel}: sites-only, with a genotype store")
    for old in ARCHIVES_281:
        got = replaced.get(old, [])
        if len(got) != 1:
            problems.append(f"{old} (2.8.1): replaced by {len(got)} archive(s) of the set "
                            f"({', '.join(got) or 'none'}), not 1")
    for kit in KITS_281:
        if kit not in replaced:
            problems.append(f"{kit} (2.8.1): replaced by no kit of the set")
    for old in sorted(replaced):
        if old in KITS_281:
            continue
        if old not in ARCHIVES_281:
            problems.append(f"{old}: replaced, but no archive of 2.8.1 has that name")
        elif crisprme_hf.successor_index(old) != replaced[old][0]:
            problems.append(f"{old}: a download of it names "
                            f"{crisprme_hf.successor_index(old)}, the set {replaced[old][0]}")
    return problems


def pending(root=ROOT):
    """What the check says of each script of PENDING_SCRIPTS this tree does not hold yet."""
    return [f"{script}: {why}, which this tree does not hold yet (the merge of its branch "
            f"brings it): the dry run takes its call as done"
            for script, why in PENDING_SCRIPTS.items()
            if not os.path.exists(os.path.join(root, script))]


def _genome_checks(used):
    """The lines of the rebuild's script that check, before its first command, the genome
    folders of the working directories it uses (README, rule 1; GENOME_FOLDERS): $REF's holds a
    FASTA of chr1-22, X or Y, and is hg38 of its own (REF_READ): no link to nothing, and neither
    the genome folder of a panel it checks nor links that lead into it at any hop, against which
    that folder would be checked against itself, a FASTA missing there missing from $REF's too
    (the PP review, round 3: B3-1; the PP fix review, round 1: B1-1, links into a folder of
    links); HGDP + 1000 Genomes' Genomes/hg38 holds the FASTAs of $REF's (all of hg38); HPRC's
    and mega's those of them that are chr1-22, X and Y, and no other (given no $REF, FASTAs of
    chr1-22, X and Y alone, the same in both: commands_script gives it $REF wherever it checks a
    panel's folder, a part of the build half too, REF_READ); and a panel whose recipe builds
    every chromosome of its sources (SOURCE_CHROMOSOMES: mega's), the chromosomes of those
    sources. Otherwise a panel's indexes would hold other contigs than 2.8.1's of the panel, and
    mega's recipe stops on a FASTA it does not build, after hours of the other panels' builds."""
    panels = [v for v in GENOME_FOLDERS if v in used]
    if not panels and "REF" not in used:
        return []
    lines = [
        "# the genome folders, before any command (seq_script/rebuild/README.md, rule 1)",
        "_chroms() { local f; for f in \"$1/${3:-}\"*\"$2\"; do [ -e \"$f\" ] || continue; "
        "f=${f##*/}; f=${f#\"${3:-}\"}; printf '%s\\n' \"${f%\"$2\"}\"; done | LC_ALL=C sort; }",
        "_stop() { echo \"rebuild: $*\" >&2; exit 2; }",
        f"_primary='{_PRIMARY_RE}'",
    ]
    if "REF" in used:
        lines += [
            "_all=$(_chroms \"$REF/Genomes/hg38\" .fa)",
            "_want=$(printf '%s\\n' $_all | grep -E \"$_primary\" || true)",
            "[ -n \"$_want\" ] || _stop \"$REF/Genomes/hg38 (\\$REF) holds no FASTA of chr1-22, "
            "X or Y: the working directory of the reference indexes holds all of hg38 there, as "
            "crisprme.py download --what genome installs it\"",
        ]
        if panels:
            # (0 where folder $1 is folder $2, or a FASTA of $1 is a link that leads into $2:
            # each hop of it read, a relative one from its own folder, and the folder it is in
            # taken as it is (cd -P). Its last target alone (readlink -f) passed a $REF of links
            # into a $2 whose own FASTAs are links to another folder, as ml008's $HGDP is: a
            # FASTA lost there was lost from $REF's too; the PP fix review, round 1: B1-1)
            lines.append(
                "_into() { local f t h d o n; o=$(cd -P \"$2\" 2>/dev/null && pwd) || return 1; "
                "[ \"$1\" -ef \"$2\" ] && return 0; for f in \"$1\"/*.fa; do t=$f; n=0; "
                "while [ -L \"$t\" ] && [ $n -lt 40 ]; do h=$(readlink -- \"$t\"); "
                "case $h in /*) ;; *) h=$(dirname -- \"$t\")/$h ;; esac; "
                "d=$(cd -P \"$(dirname -- \"$h\")\" 2>/dev/null && pwd); "
                "case \"$d/\" in \"$o\"/*) return 0 ;; esac; t=$h; n=$((n + 1)); done; done; "
                "return 1; }")
        for v in panels:
            lines.append(
                f"! _into \"$REF/Genomes/hg38\" \"${v}/Genomes/hg38\" || _stop \"$REF/Genomes/hg38 "
                f"(\\$REF) is ${v}/Genomes/hg38 (\\${v}), or holds links to its FASTAs: that "
                f"folder would be checked against itself, a FASTA missing there missing from "
                f"\\$REF's too. \\$REF's Genomes/hg38 holds hg38 of its own: the reference's "
                f"working directory, a copy of hg38 (crisprme.py download --what genome), or "
                f"links to the FASTAs of another such folder\"")
        lines.append(
            "for _f in \"$REF/Genomes/hg38\"/*.fa; do [ -e \"$_f\" ] || ! [ -L \"$_f\" ] || "
            "_stop \"$_f (\\$REF) is a link to nothing: \\$REF's Genomes/hg38 holds all of hg38, "
            "and the genome folders here are checked against it\"; done")
    else:
        lines.append("_all=; _want=")
    panel_of = {name: panel for panel, name in WORKDIR.items()}
    for v in panels:
        panel = panel_of[v]
        lines += [f"_got=$(_chroms \"${v}/Genomes/hg38\" .fa)",
                  f"[ -n \"$_got\" ] || _stop \"${v}/Genomes/hg38 (\\${v}) holds no FASTA\""]
        if GENOME_FOLDERS[v] == "REF":
            lines.append(
                f"[ -z \"$_all\" ] || [ \"$_got\" = \"$_all\" ] || _stop \"${v}/Genomes/hg38 "
                f"(\\${v}) holds $(set -- $_got; echo $#) FASTA(s), not the $(set -- $_all; echo $#) of "
                f"\\$REF's: the working directory of {panel} holds all of hg38, as the "
                f"reference's does, as 2.8.1's default index of the panel was built\"")
            continue
        refuses = "; its recipe refuses any other FASTA" if v in SOURCE_CHROMOSOMES else ""
        lines += [
            "_bad=$(printf '%s\\n' $_got | grep -vE \"$_primary\" || true)",
            f"[ -z \"$_bad\" ] || _stop \"${v}/Genomes/hg38 (\\${v}) holds $(echo $_bad): the "
            f"working directory of {panel} holds the chromosomes chr1-22, X and Y of hg38 alone, "
            f"one FASTA each, as 2.8.1's indexes of the panel were built{refuses}\"",
            "[ -n \"$_want\" ] || _want=$_got",
            f"[ \"$_got\" = \"$_want\" ] || _stop \"${v}/Genomes/hg38 (\\${v}) holds $(echo $_got),"
            f" not $(echo $_want): the working directory of {panel} holds the chromosomes chr1-22, "
            f"X and Y of the reference's\"",
        ]
    for v, (src, prefix, suffix) in SOURCE_CHROMOSOMES.items():
        if v in panels and src in used:
            lines += [
                f"_src=$(_chroms \"${src}\" {suffix} {prefix})",
                f"[ \"$_src\" = \"$(_chroms \"${v}/Genomes/hg38\" .fa)\" ] || _stop \"${src} "
                f"(\\${src}) holds {prefix}<chr>{suffix} of $(echo $_src), and ${v}/Genomes/hg38 "
                f"the FASTAs of $(echo $(_chroms \"${v}/Genomes/hg38\" .fa)): the recipe builds "
                f"every chromosome of its sources, and refuses a genome folder that holds the "
                f"FASTA of another or lacks one of theirs\"",
            ]
    return lines


def _mark_name(half, archive):
    """The name of the mark a half of the rebuild leaves for row ``archive`` once its command
    has run (in the row's working directory's MARKS folder): <half>.<archive>, a '/' as '_'."""
    return f"{half}.{archive.replace('/', '_')}"


def _done_checks(half):
    """The lines of a half of the rebuild that tell whether a row ran to its end before, so that
    a re-run skips it (README, "A re-run"): the mark the half leaves in the row's working
    directory once its command has run (MARKS: _mark), with, for an index, the index still in
    genome_library (_done); for the build of a kit, the kit's SHA256SUMS holding, written once
    its tarball is whole, the last thing its command does (_kit_done). Nothing else: what a
    recipe leaves does not tell that it ran to its end. mega's recipe writes every part of its
    dataset before its last step, its checks (verify_mega_chrom.py), so a build that failed them
    left every part there, and a re-run that took its parts for done passed it with exit 0; its
    line now runs the recipe, which refuses a working directory that holds any part of the
    dataset, and the half stops there (the PP-2 review, A1-PP2-1)."""
    lines = [
        f"# a re-run skips each row done before (seq_script/rebuild/README.md, \"A re-run\"): "
        f"its mark in its working directory ({MARKS}/), with its index",
        f"_mark() {{ mkdir -p \"$1/{MARKS}\" && : > \"$1/{MARKS}/$2\"; }}",
        f"_done() {{ [ -e \"$1/{MARKS}/$2\" ] && {{ [ -z \"${{3:-}}\" ] || "
        f"[ -d \"$1/genome_library/$3\" ]; }}; }}",
        "_skip() { echo \"rebuild: the $1: done before, skipped\"; }",
    ]
    if half == "build":
        lines += [
            "# a kit whose SHA256SUMS holds, which its command writes last, once the tarball is "
            "whole",
            "_kit_done() { [ -f \"$1/$2.SHA256SUMS\" ] && [ -f \"$1/$2.tar.gz\" ] "
            "&& (cd \"$1\" && sha256sum -c --status \"$2.SHA256SUMS\"); }",
        ]
    return lines


# What a clean environment's PATH is: the publish half runs the image's own hf only with it (a login
# shell's PATH puts the hf of ~/.local/bin ahead of the image's, which apptainer passes in)
CLEAN_PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"


def _publish_checks():
    """The publish half's preflight, before its first upload (README, "The publish"; PP-1 and
    PP-2 of the release). The image sets HF_HUB_DISABLE_IMPLICIT_TOKEN=1, so a bare `hf upload`
    in it sends no token, and from a login shell `$RUN hf` ran the host's ~/.local/bin/hf: the
    kits' uploads failed as documented, after the indexes' had gone up (the D7 readiness audit,
    R-PUB-1, G-P5). So: HF_TOKEN in the environment (never on a command line: no command of the
    set passes one), a stop with exit 2 where it is not set, as for each check here, and no line
    that expands it, so that a trace (bash -x, SHELLOPTS=xtrace) writes none (the PP-1 review,
    A1-PP1-1 and A1-PP1-2: bash's own check of an unset variable exited 1, and its argument was
    the token); its override for the image's hf, exported for apptainer
    (APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0; Docker's RUN gives -e HF_TOKEN -e
    HF_HUB_DISABLE_IMPLICIT_TOKEN=0), and PYTHONNOUSERSITE=1, so that the image's python reads
    no package of ~/.local; then, in the image, the token there and not refused by that setting,
    the hf it runs installed beside its crisprme.py (the image's own), the token authenticating
    (hf auth whoami), and the repo of index format 2 there (hf datasets info): its owner creates
    it, and `hf upload` would create one that is not there."""
    repo = crisprme_hf.FORMAT2_HF_REPO
    return [
        "# the preflight (seq_script/rebuild/README.md, \"The publish\"): the token in the "
        "environment, sent by the image's own hf, and the repo there",
        "_stop() { echo \"publish: $*\" >&2; exit 2; }",
        f'[ -n "${{HF_TOKEN:+x}}" ] || _stop "HF_TOKEN is not set: the token that publishes, a '
        f'write token of {repo}, exported in the environment: never on a command line, which ps '
        f'shows"',
        "export APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0 "
        "SINGULARITYENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0 PYTHONNOUSERSITE=1",
        "$RUN sh -c '[ -n \"${HF_TOKEN:+x}\" ] && case $(printf %s "
        "\"${HF_HUB_DISABLE_IMPLICIT_TOKEN:-}\" | tr a-z A-Z) in 1|ON|YES|TRUE) exit 1;; esac' "
        "|| _stop \"the release image would upload with no token: HF_TOKEN does not reach it, or "
        "its HF_HUB_DISABLE_IMPLICIT_TOKEN=1 holds there (apptainer takes this script's "
        "APPTAINERENV_HF_HUB_DISABLE_IMPLICIT_TOKEN=0 unless run with --cleanenv; a RUN of Docker "
        "gives docker run -e HF_TOKEN -e HF_HUB_DISABLE_IMPLICIT_TOKEN=0)\"",
        "_hf=$($RUN sh -c 'command -v hf') && _me=$($RUN sh -c 'command -v crisprme.py') "
        "|| _stop \"no hf, or no crisprme.py, where \\$RUN runs them\"",
        "[ \"${_hf%/*}\" = \"${_me%/*}\" ] || _stop \"$_hf is the hf \\$RUN would run, not the "
        "release image's own, installed beside its crisprme.py ($_me): run this script with a "
        f"clean PATH, PATH={CLEAN_PATH} bash <this script> (a login shell's PATH puts the hf of "
        "~/.local/bin first)\"",
        "$RUN hf auth whoami > /dev/null || _stop \"HF_TOKEN does not authenticate on the Hub "
        "(hf auth whoami)\"",
        f"$RUN hf datasets info {repo} > /dev/null || _stop \"{repo} is not on the Hub, or "
        f"HF_TOKEN cannot read it: its owner creates it, with its card, before the first publish "
        f"(docs/hf/README.md), and gives the token write access to it\"",
    ]


def _row_lines(r, half):
    """Row ``r``'s command of ``half`` as the rebuild runs it: the line itself, unchanged (its
    command, then || exit with what failed, then its archive and note as a comment), run where
    the row is not done before (_done_checks), and marked once it has run, the if, else and fi
    on lines of their own (half_lines leaves them out)."""
    note = "" if r["notes"] == NO_COMMAND else f": {r['notes']}"
    line = (f"{r[half]} || {{ echo \"rebuild: the {half} of {r['archive']} failed\" >&2; "
            f"exit 1; }}    # {r['archive']}{note}")
    w = "\"$" + WORKDIR[r["panel"]] + "\""
    mark = _mark_name(half, r["archive"])
    done = f"_done {w} {mark}" + (f" {r['archive']}" if r["kind"] in INDEX_KINDS else "")
    if half == "build" and r["kind"] == "build-inputs":
        done += f" || _kit_done {w} {r['archive'].split('/', 1)[1]}"
    return [f"if {done}; then _skip \"{half} of {r['archive']}\"; else", line,
            f"_mark {w} {mark}; fi"]


def panel_rows(rows, panels):
    """The rows of ``rows`` whose working directory is one of ``panels`` (the variables that name
    them, WORKDIR's: REF, HGDP, HPRC, MEGA), in the manifest's order; every row for None. A
    ValueError for a name that is no working directory's, or none."""
    if panels is None:
        return list(rows)
    known = list(dict.fromkeys(WORKDIR.values()))
    bad = [p for p in panels if p not in known]
    if bad or not panels:
        raise ValueError(f"--panels takes the working directories of the rows to print, by the "
                         f"variables that name them, comma-separated: {', '.join(known)} "
                         f"(not {', '.join(bad) or 'none'})")
    return [r for r in rows if WORKDIR[r["panel"]] in panels]


def commands_script(rows, manifest=MANIFEST, half="build", panels=None):
    """A half of the rebuild, ``half`` "build" or "publish", as one script bash runs from the
    repository's root. The build half runs every build command, and publishes nothing; the
    publish half, after the gates and the go, every publish command; each in the manifest's
    order, each line as the manifest gives it (_row_lines). With ``panels`` (the variables of
    working directories: panel_rows), the rows of those working directories alone, as a node
    that holds them runs its part of the half (README, "The publish": HGDP + 1000 Genomes' is on
    one node, the others on another, so that no one run of the whole half reads them all). It
    stops before its first command where a variable its commands read is not set, or, the build
    half, where a working directory's genome folder is not what the rebuild builds on
    (_genome_checks), a part of the build half as the whole does: against $REF's, which it asks
    for where it names no row of the reference's working directory (REF_READ); it skips a row
    done before (_done_checks); and it stops at the first command that fails: each ends in ||
    exit, as set -e does not stop at a command of an && list but its last."""
    if half not in HALVES:
        raise ValueError(f"no half {half!r} of the rebuild: {' or '.join(HALVES)}")
    rows = panel_rows(rows, panels)
    used = []
    for r in rows:
        if r[half] != NO_COMMAND:
            used += [v for v in variables(r[half]) if v not in used]
    ref_read = half == "build" and "REF" not in used and any(v in used for v in GENOME_FOLDERS)
    if ref_read:
        used.append("REF")
    lines = ["#!/usr/bin/env bash",
             f"# The rebuild of the index set of CRISPRme 2.10.0 ({os.path.basename(manifest)}), "
             f"its {half} half, as check_index_set.py --commands {half}"
             + (f" --panels {','.join(panels)}" if panels is not None else "")]
    if half == "build":
        lines += ["# prints it: every build command, in the manifest's order "
                  "(seq_script/rebuild/README.md). It publishes nothing: the publish half",
                  "# (--commands publish) is a script of its own, run once the indexes have "
                  "passed the release's gates and the publish is approved."]
    else:
        lines += ["# prints it: every publish command, in the manifest's order "
                  "(seq_script/rebuild/README.md), each into the repo of index format 2.",
                  "# Run it once, after the build half has built the set and the indexes have "
                  "passed the release's gates, with the publish approved."]
    if panels is not None:
        lines.append(f"# The rows of {', '.join('$' + p for p in panels)} alone (--panels), "
                     f"for the node that holds them.")
    if ref_read:
        lines.append("# $REF is read, not written: the genome folders of these rows' working "
                     "directories are checked against its Genomes/hg38, as the whole half checks "
                     "them (seq_script/rebuild/README.md, rule 1).")
    lines += ["# Run it from the repository's root, with the variables below set. A re-run skips "
              "each row done before.",
              "set -euo pipefail",
              "[ -f seq_script/rebuild/check_index_set.py ] || { echo \"run it from the "
              "repository's root\" >&2; exit 2; }"]
    for name in sorted(used, key=lambda v: list(VARIABLES).index(v) if v in VARIABLES else 99):
        text = REF_READ if ref_read and name == "REF" else \
            VARIABLES.get(name, "a variable of the manifest")
        lines.append(f': "${{{name}{"?" if name in MAY_BE_EMPTY else ":?"}{text}}}"')
    if half == "build":
        lines += _genome_checks(used)
    else:
        lines += _publish_checks()
    lines += _done_checks(half)
    lines.append(f"# {half}")
    for r in rows:
        if r[half] != NO_COMMAND:
            lines += _row_lines(r, half)
    return "\n".join(lines) + "\n"


def half_lines(script):
    """The rows' lines of a half of the rebuild (commands_script), as the manifest gives them:
    those after its "# build" or "# publish" line, the if, else and mark-and-fi lines around each
    left out. The build half's are what the D7 readiness audit checked (their md5)."""
    body = script.split("\n# build\n", 1)[-1].split("\n# publish\n", 1)[-1]
    return "".join(line + "\n" for line in body.splitlines()
                   if not (line.startswith("if ") and line.endswith("; else"))
                   and not (line.startswith("_mark ") and line.endswith("; fi")))


# The stub of every program the rebuild's commands run, for the dry run: one Python script, which
# does what its name (the link's) says it does, logs it, and stops the script, saying why, where a
# path it is handed is not there. bash runs a recipe as a stub (its own steps are its own tests');
# any other bash is the real one.
_STUB = r'''#!{python} -I
import os, sys
name = os.path.basename(sys.argv[0])
args = sys.argv[1:]
tmp = os.environ["CRISPRME_DRY_RUN_TMP"]


def log(*parts):
    with open(os.environ["CRISPRME_DRY_RUN_LOG"], "a") as fh:
        fh.write("\t".join([name, os.getcwd()] + [str(p) for p in parts]) + "\n")


def fail(why):
    log("FAIL", why)
    sys.stderr.write("dry run: %s: %s\n" % (name, why))
    sys.exit(1)


def need(path, what, test=os.path.exists):
    if not test(path):
        fail("%s %s is not there (the command runs in %s)" % (what, path, os.getcwd()))
    return os.path.abspath(path)


def opt(flag):
    if flag not in args[:-1]:
        fail("no %s" % flag)
    return args[args.index(flag) + 1]


if name in ("crisprme.py", "hf", "huggingface-cli", "pigz") \
        and os.environ.get("CRISPRME_REBUILD_IMAGE") != "1":
    fail("runs outside the release image: a command runs it through $RUN")
if name == "crisprme.py":
    if args[:1] == ["build-index-only"]:
        need(opt("--genome"), "the genome folder", os.path.isdir)
        log("build-index-only", need(opt("--pam"), "the PAM file", os.path.isfile),
            need(opt("--path"), "the working directory", os.path.isdir))
    elif args[:1] == ["publish-index"]:
        log("publish-index", need(opt("--index"), "the index", os.path.isdir))
    else:
        fail("crisprme.py %s is not a command of the rebuild" % " ".join(args[:1]))
elif name == "hf" and args[:2] == ["auth", "whoami"]:
    log("whoami")
elif name == "hf" and args[:2] == ["datasets", "info"] and len(args) > 2:
    log("repo", args[2])
elif name in ("hf", "huggingface-cli"):
    if args[:1] != ["upload"] or len(args) < 4:
        fail("not hf upload <repo> <file> <path in the repo>")
    log("upload", args[1], need(args[2], "the file to upload", os.path.isfile), args[3])
elif name == "pigz":
    sys.stdout.buffer.write(sys.stdin.buffer.read())
elif name == "tar":
    if args[:2] != ["-cf", "-"] or len(args) < 3:
        fail("not tar -cf - <folders>")
    log("pack", *[need(p, "the folder to pack", os.path.isdir) for p in args[2:]])
elif name == "sha256sum" and args[:1] == ["-c"]:
    # (a check of a list: whole where every file it names is there, as a dry run takes them)
    listed = [a for a in args[1:] if not a.startswith("-")]
    if len(listed) != 1 or not os.path.isfile(listed[0]):
        sys.exit(1)
    with open(listed[0]) as fh:
        names = [line.split(None, 1)[1].strip() for line in fh if line.strip()]
    log("check", os.path.abspath(listed[0]))
    sys.exit(0 if names and all(os.path.isfile(n) for n in names) else 1)
elif name == "sha256sum":
    for p in args:
        need(p, "the file to hash", os.path.isfile)
        print("0" * 64 + "  " + p)
elif name == "bash" and args and args[0].endswith(".sh") and not args[0].startswith("-"):
    script = args[0]
    if script in os.environ.get("CRISPRME_DRY_RUN_PENDING", "").split(":"):
        log("pending", script)
        sys.exit(0)
    need(script, "the script", os.path.isfile)
    repo = os.environ["CRISPRME_DRY_RUN_REPO"]
    if not os.path.realpath(script).startswith(repo + os.sep):
        fail("%s is not a script of the repository" % script)
    for key, value in sorted(os.environ.items()):
        if key != "PATH" and not key.startswith("CRISPRME_DRY_RUN_") \
                and value.startswith(tmp + os.sep):
            need(value, "%s=" % key)
    log("recipe", script, *[need(x, "the argument") if x.startswith(tmp + os.sep) else x
                            for x in args[1:]])
elif name == "bash":
    os.execv({bash!r}, [{bash!r}] + args)
else:
    fail("runs the image itself: a command runs the release image through $RUN")
'''
_STUBBED = ("crisprme.py", "hf", "huggingface-cli", "pigz", "tar", "sha256sum", "bash",
            "apptainer", "singularity")


def dry_run(rows, root=ROOT, manifest=MANIFEST):
    """The problems of running the rebuild dry, its build half, then its publish half
    (commands_script), then each again: every program their commands run stubbed (_STUB), from a
    copy of the repository's root made of links, with each variable set to a folder or file of a
    temporary directory: working directories that hold hg38's folder and what the builds leave
    (each index of the set, each panel's enriched genomes), and the recipes' inputs. Each half
    runs to its end only where every script it calls is there (a pending one, PENDING_SCRIPTS,
    taken as done), every PAM file a build reads is in place, every path a command hands a
    program or a recipe is there, crisprme.py, hf and pigz run in the release image, and every
    index of the set is published from its working directory, by the publish half alone (the
    build half publishes nothing, the publish half builds nothing); neither writes into the
    repository's root. Run again, each skips every row it ran (README, "A re-run"), and the build
    half skips a kit whose SHA256SUMS holds with no mark (as after a build run row by row)."""
    bash = shutil.which("bash")
    if not bash:
        return ["no bash to run the rebuild's script with"]
    problems = []
    with tempfile.TemporaryDirectory(prefix="check_index_set.") as t:
        t = os.path.realpath(t)
        top = os.path.join(t, "root")
        os.mkdir(top)
        links = sorted(n for n in os.listdir(root) if n != ".git")
        for n in links:
            os.symlink(os.path.join(root, n), os.path.join(top, n))
        stubs = os.path.join(t, "bin")
        os.mkdir(stubs)
        stub = os.path.join(t, "stub")
        with open(stub, "w") as fh:
            fh.write(_STUB.replace("{python}", sys.executable).replace("{bash!r}", repr(bash)))
        os.chmod(stub, 0o755)
        for n in _STUBBED:
            os.symlink(stub, os.path.join(stubs, n))
        env = {"PATH": os.pathsep.join([stubs, "/usr/bin", "/bin", "/usr/sbin", "/sbin"]),
               "HOME": os.path.join(t, "home"), "LC_ALL": "C",
               "RUN": "env CRISPRME_REBUILD_IMAGE=1", "HF_TOKEN": "dry-run-not-a-token",
               "CRISPRME_DRY_RUN_TMP": t, "CRISPRME_DRY_RUN_REPO": os.path.realpath(root),
               "CRISPRME_DRY_RUN_LOG": os.path.join(t, "log"),
               "CRISPRME_DRY_RUN_PENDING": ":".join(
                   s for s in PENDING_SCRIPTS if not os.path.exists(os.path.join(root, s)))}
        os.mkdir(env["HOME"])
        open(env["CRISPRME_DRY_RUN_LOG"], "w").close()
        work = {}
        for panel, name in WORKDIR.items():
            work[panel] = env[name] = os.path.join(t, "work", name)
            genome = os.path.join(work[panel], "Genomes", "hg38")
            os.makedirs(genome)
            # (README, rule 1: all of hg38, as download installs it, a mitochondrion and an alt
            # contig among the rest, for the reference and HGDP + 1000 Genomes; chr1-22, X and Y
            # for HPRC and mega)
            every = GENOME_FOLDERS.get(name, "REF") == "REF"
            for c in PRIMARY + (("chrM", "chr10_GL383545v1_alt") if every else ()):
                open(os.path.join(genome, f"{c}.fa"), "w").close()
        for r in rows:
            w = work.get(r["panel"])
            if w and r["kind"] in INDEX_KINDS:
                os.makedirs(os.path.join(w, "genome_library", r["archive"]), exist_ok=True)
            if w and r["panel"] != "-":
                for g in (f"hg38+{r['panel']}", f"hg38+{r['panel']}_INDELS"):
                    os.makedirs(os.path.join(w, "Genomes", g), exist_ok=True)
        inputs = os.path.join(t, "inputs")
        for name in ("RP", "HPRC_SRC", "MEGA_SRC"):
            env[name] = os.path.join(inputs, name)
            os.makedirs(env[name])
        for v, (src, prefix, suffix) in SOURCE_CHROMOSOMES.items():
            for c in PRIMARY:
                open(os.path.join(env[src], f"{prefix}{c}{suffix}"), "w").close()
        for name in ("SIF", "META", "P1KG"):
            env[name] = os.path.join(inputs, name)
            open(env[name], "w").close()
        # (a variable VARIABLES does not hold, which the check says, a folder: the run goes on)
        for r in rows:
            for what in HALVES:
                for name in variables(r[what]):
                    if name not in env:
                        env[name] = os.path.join(inputs, name)
                        os.makedirs(env[name])

        def run(half, label):
            """(exit status, the log lines it wrote) of a run of ``half``; a stop, a command
            that failed and went on, said as problems."""
            script = os.path.join(t, f"rebuild_{half}.sh")
            with open(script, "w") as fh:
                fh.write(commands_script(rows, manifest, half))
            with open(env["CRISPRME_DRY_RUN_LOG"]) as fh:
                before = len(fh.readlines())
            try:
                got = subprocess.run([bash, script], cwd=top, env=env, capture_output=True,
                                     text=True, timeout=600)
            except subprocess.TimeoutExpired:
                problems.append(f"the rebuild's {half} half, run dry{label}, did not end in 600 s")
                return None, []
            with open(env["CRISPRME_DRY_RUN_LOG"]) as fh:
                log = [line.rstrip("\n").split("\t") for line in fh][before:]
            if got.returncode != 0:
                why = [x for x in got.stderr.splitlines() if x.strip()][-3:]
                problems.append(f"the rebuild's {half} half, run dry{label} from the repository's "
                                f"root, stops (exit {got.returncode}): "
                                f"{' | '.join(why) or 'no message'}")
            failed = [x[3] for x in log if len(x) > 3 and x[2] == "FAIL"]
            if failed and got.returncode == 0:
                problems.append(f"the rebuild's {half} half, run dry{label}, goes on past a "
                                f"command that failed: {' | '.join(failed)}")
            return got.returncode, log

        def acts(log, *kinds):
            """The entries of ``log`` that build (a build, a recipe, a pack) or publish (an
            upload, a publish-index), as ``kinds`` names them."""
            what = {"build": ("build-index-only", "recipe", "pending", "pack"),
                    "publish": ("publish-index", "upload")}
            wanted = {k for kind in kinds for k in what[kind]}
            return [x[2] if x[2] in wanted else x[0] for x in log
                    if len(x) > 2 and (x[2] in wanted or (x[0] == "tar" and "build" in kinds))]

        rc, log = run("build", "")
        if rc == 0 and acts(log, "publish"):
            problems.append(f"the rebuild's build half, run dry, publishes: "
                            f"{', '.join(sorted(set(acts(log, 'publish'))))}")
        rc, log = run("publish", "")
        if rc == 0:
            first = next((i for i, x in enumerate(log) if len(x) > 2
                          and x[2] in ("publish-index", "upload")), None)
            checked = {x[2] for x in log[:first] if len(x) > 2}
            if first is not None and not {"whoami", "repo"} <= checked:
                problems.append("the rebuild's publish half, run dry, publishes before it checks "
                                "that the token authenticates and the repo is there")
            if acts(log, "build"):
                problems.append(f"the rebuild's publish half, run dry, builds: "
                                f"{', '.join(sorted(set(acts(log, 'build'))))}")
            published = {x[3] for x in log if x[0] == "crisprme.py" and x[2] == "publish-index"}
            for r in rows:
                if r["kind"] in INDEX_KINDS:
                    want = os.path.join(work[r["panel"]], "genome_library", r["archive"])
                    if want not in published:
                        problems.append(f"{r['archive']}: the rebuild's publish half, run dry, "
                                        f"does not publish {want}")
        if not problems:
            for half, kind in (("build", "build"), ("publish", "publish")):
                rc, log = run(half, " again")
                if rc == 0 and acts(log, kind):
                    problems.append(f"the rebuild's {half} half, run dry again, does again what "
                                    f"it did: {', '.join(sorted(set(acts(log, kind))))}")
        if not problems:
            # (a build run row by row leaves no mark: a kit whose SHA256SUMS holds is done)
            for r in rows:
                w = work.get(r["panel"])
                if w and r["kind"] == "build-inputs":
                    os.remove(os.path.join(w, MARKS, _mark_name("build", r["archive"])))
            rc, log = run("build", " again, without the marks of the kits")
            if rc == 0 and acts(log, "build"):
                problems.append(f"the rebuild's build half, run dry again without the marks of "
                                f"the kits, builds a kit whose SHA256SUMS holds: "
                                f"{', '.join(sorted(set(acts(log, 'build'))))}")
        extra = sorted(set(os.listdir(top)) - set(links))
        if extra:
            problems.append(f"the rebuild's script, run dry, writes into the repository's root: "
                            f"{', '.join(extra)}")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--commands", nargs="?", const="", metavar="build|publish",
                    help="print a half of the rebuild, a script bash runs from the repository's "
                         "root: build, every build command (it publishes nothing); publish, "
                         "every publish command, run once the gates have passed and the publish "
                         "is approved")
    ap.add_argument("--panels", metavar="REF,HGDP,HPRC,MEGA",
                    help="with --commands: the rows of these working directories alone (the "
                         "variables that name them, comma-separated), in the manifest's order, "
                         "for the node that holds them: --panels REF,HPRC,MEGA on one node and "
                         "--panels HGDP on the other (seq_script/rebuild/README.md, The publish). "
                         "A part of the build half asks for $REF too, whose genome folder theirs "
                         "are checked against (README, Commands)")
    args = ap.parse_args(argv)
    if args.commands is not None and args.commands not in HALVES:
        ap.error("--commands takes the half of the rebuild to print: --commands build (every "
                 "build command; it publishes nothing) or --commands publish (every publish "
                 "command, run once the indexes have passed the gates and the publish is "
                 "approved)")
    panels = None
    if args.panels is not None:
        if not args.commands:
            ap.error("--panels chooses the rows a half prints: give it with --commands build or "
                     "--commands publish")
        panels = list(dict.fromkeys(p.strip() for p in args.panels.split(",") if p.strip()))
        try:
            panel_rows([], panels)
        except ValueError as exc:
            ap.error(str(exc))
    rows = read_manifest(args.manifest)
    if args.commands:
        sys.stdout.write(commands_script(rows, args.manifest, args.commands, panels))
        return 0
    problems = check(rows)
    problems += dry_run(rows, manifest=args.manifest)
    for note in pending():
        print(f"NOTE: {note}")
    for p in problems:
        print(f"PROBLEM: {p}")
    n_idx = sum(r["kind"] in INDEX_KINDS for r in rows)
    print(f"index set {os.path.basename(args.manifest)}: {len(rows)} archives "
          f"({n_idx} indexes), {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())

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

Its commands are what the rebuild runs, and ``--commands`` prints them as one script: each a
bash command (``bash -n``; a note has a column of its own) run from the repository's root on its
own, so that none moves the next (a ``cd`` only inside ``( )``); a recipe called by its repository
path, and given every input it would otherwise read from a path of the cluster it was written on;
each build first copying its PAM file from the repository's ``PAMs/`` into its working
directory's, and reading it there; ``crisprme.py``, ``hf`` and ``pigz`` run in the release image
(``$RUN``); each command in its own panel's working directory, with the variables of VARIABLES
alone; a panel companion, which its panel's first build writes and first publish uploads, with
none ("-"). The check reads every command so, then runs the script dry (dry_run): every program
it calls stubbed, in working directories that hold what the builds leave, from a copy of the
repository's root made of links. The script runs to its end only where every script it calls is
there, every PAM file a build reads is in place, and every file a publish uploads is where the
command says.

    python seq_script/rebuild/check_index_set.py              # the check and the dry run; exit 1 on a problem
    python seq_script/rebuild/check_index_set.py --commands   # the rebuild: a script bash runs from the repository's root

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
# A recipe the set calls that another lane of stage 2b brings, which this tree does not hold
# yet: the integration of stage 2b brings it, and its entry here goes then (the check says so
# once the tree holds it). Each is said, not counted as a problem, while it is missing.
PENDING_SCRIPTS = {
    "seq_script/merge_panels/mega_rebuild.sh":
        "lane M's recipe of the mega rebuild (decision D2: PASS/'.' records, AF = AF_max), "
        "on branch s2b-M",
}


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
    for name, text in VARIABLES.items():
        if _VARIABLE_TEXT.search(text):
            problems.append(f"VARIABLES[{name!r}] holds one of }} $ ` \" ', which the rebuild's "
                            f"script cannot hold in ${{{name}:?...}}")
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
    return [f"{script}: {why}, which this tree does not hold yet (the integration of stage 2b "
            f"brings it): the dry run takes its call as done"
            for script, why in PENDING_SCRIPTS.items()
            if not os.path.exists(os.path.join(root, script))]


def _genome_checks(used):
    """The lines of the rebuild's script that check, before its first command, the genome
    folders of the working directories it uses (README, rule 1; GENOME_FOLDERS): HGDP + 1000
    Genomes' Genomes/hg38 holds the FASTAs of $REF's (all of hg38); HPRC's and mega's those of
    them that are chr1-22, X and Y, and no other; where the script does not use $REF, HPRC's
    and mega's are FASTAs of chr1-22, X and Y alone, the same in both (on a test genome of chr22,
    chr22); and a panel whose recipe builds every chromosome of its sources (SOURCE_CHROMOSOMES:
    mega's), the chromosomes of those sources. Otherwise a panel's indexes would hold other
    contigs than 2.8.1's of the panel, and mega's recipe stops on a FASTA it does not build,
    after hours of the other panels' builds."""
    panels = [v for v in GENOME_FOLDERS if v in used]
    if not panels:
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


def commands_script(rows, manifest=MANIFEST):
    """The rebuild as one script bash runs from the repository's root: it stops before its first
    command where a variable the commands read is not set, or a working directory's genome
    folder is not what the rebuild builds on (_genome_checks), then runs every build command,
    then every publish command, in the manifest's order, and stops at the first that fails:
    each ends in || exit, as set -e does not stop at a command of an && list but its last."""
    used = []
    for r in rows:
        for what in ("build", "publish"):
            if r[what] != NO_COMMAND:
                used += [v for v in variables(r[what]) if v not in used]
    lines = ["#!/usr/bin/env bash",
             f"# The rebuild of the index set of CRISPRme 2.10.0 ({os.path.basename(manifest)}), "
             f"as check_index_set.py --commands",
             "# prints it: every build command, then every publish command, in the manifest's "
             "order (seq_script/rebuild/README.md).",
             "# Run it from the repository's root, with the variables below set.",
             "set -euo pipefail",
             "[ -f seq_script/rebuild/check_index_set.py ] || { echo \"run it from the "
             "repository's root\" >&2; exit 2; }"]
    for name in sorted(used, key=lambda v: list(VARIABLES).index(v) if v in VARIABLES else 99):
        text = VARIABLES.get(name, "a variable of the manifest")
        lines.append(f': "${{{name}{"?" if name in MAY_BE_EMPTY else ":?"}{text}}}"')
    lines += _genome_checks(used)
    for what in ("build", "publish"):
        lines.append(f"# {what}")
        for r in rows:
            if r[what] != NO_COMMAND:
                note = "" if r["notes"] == NO_COMMAND else f": {r['notes']}"
                lines.append(f"{r[what]} || {{ echo \"rebuild: the {what} of {r['archive']} "
                             f"failed\" >&2; exit 1; }}    # {r['archive']}{note}")
    return "\n".join(lines) + "\n"


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
    """The problems of running the rebuild (commands_script) dry: every program its commands
    run stubbed (_STUB), from a copy of the repository's root made of links, with each variable
    set to a folder or file of a temporary directory: working directories that hold hg38's folder
    and what the builds leave (each index of the set, each panel's enriched genomes), and the
    recipes' inputs. It runs to its end only where every script it calls is there (a pending one,
    PENDING_SCRIPTS, taken as done), every PAM file a build reads is in place, every path a
    command hands a program or a recipe is there, crisprme.py, hf and pigz run in the release
    image, and every index of the set is published from its working directory; and it writes
    nothing into the repository's root."""
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
               "RUN": "env CRISPRME_REBUILD_IMAGE=1",
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
            for what in ("build", "publish"):
                for name in variables(r[what]):
                    if name not in env:
                        env[name] = os.path.join(inputs, name)
                        os.makedirs(env[name])
        script = os.path.join(t, "rebuild.sh")
        with open(script, "w") as fh:
            fh.write(commands_script(rows, manifest))
        try:
            got = subprocess.run([bash, script], cwd=top, env=env, capture_output=True,
                                 text=True, timeout=600)
        except subprocess.TimeoutExpired:
            return ["the rebuild's script, run dry, did not end in 600 s"]
        if got.returncode != 0:
            why = [x for x in got.stderr.splitlines() if x.strip()][-3:]
            problems.append(f"the rebuild's script, run dry from the repository's root, stops "
                            f"(exit {got.returncode}): {' | '.join(why) or 'no message'}")
        log = []
        if os.path.exists(env["CRISPRME_DRY_RUN_LOG"]):
            with open(env["CRISPRME_DRY_RUN_LOG"]) as fh:
                log = [line.rstrip("\n").split("\t") for line in fh]
        failed = [x[3] for x in log if len(x) > 3 and x[2] == "FAIL"]
        if failed and got.returncode == 0:
            problems.append(f"the rebuild's script, run dry, goes on past a command that "
                            f"failed: {' | '.join(failed)}")
        if got.returncode == 0:
            published = {x[3] for x in log if x[0] == "crisprme.py" and x[2] == "publish-index"}
            for r in rows:
                if r["kind"] in INDEX_KINDS:
                    want = os.path.join(work[r["panel"]], "genome_library", r["archive"])
                    if want not in published:
                        problems.append(f"{r['archive']}: the rebuild's script, run dry, does "
                                        f"not publish {want}")
        extra = sorted(set(os.listdir(top)) - set(links))
        if extra:
            problems.append(f"the rebuild's script, run dry, writes into the repository's root: "
                            f"{', '.join(extra)}")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--commands", action="store_true",
                    help="print the rebuild, a script bash runs from the repository's root: "
                         "every build command, then every publish command, in order")
    args = ap.parse_args(argv)
    rows = read_manifest(args.manifest)
    if args.commands:
        sys.stdout.write(commands_script(rows, args.manifest))
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

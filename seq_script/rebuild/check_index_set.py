#!/usr/bin/env python3
"""check_index_set.py - check the index set of CRISPRme 2.10.0 (the rebuild manifest) against the
code that builds, names, publishes and downloads its archives.

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

    python seq_script/rebuild/check_index_set.py              # the check; exit 1 on a problem
    python seq_script/rebuild/check_index_set.py --commands   # the build, then publish, commands

Standard library and the repository's PostProcess modules only; no network.
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MANIFEST = os.path.join(HERE, "index_set_2.10.0.tsv")
sys.path.insert(0, os.path.join(ROOT, "PostProcess"))

import crisprme_hf  # noqa: E402
import index_fit  # noqa: E402

COLUMNS = ("archive", "kind", "panel", "pam_file", "pam", "side", "spacer", "bmax", "bdna",
           "brna", "repo", "replaces", "build", "publish")
INDEX_KINDS = ("reference", "variant")
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


def check(rows, root=ROOT):
    """The problems of the index set ``rows`` (a list of strings; none: it holds)."""
    problems = []
    names = [r["archive"] for r in rows]
    for dup in sorted({n for n in names if names.count(n) > 1}):
        problems.append(f"{dup}: listed {names.count(dup)} times")
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
            if r["pam_file"] not in r["build"] and os.path.basename(r["pam_file"]) not in r["build"]:
                problems.append(f"{a}: its build command does not name {r['pam_file']}")
            if kind == "reference" and f"--bDNA {bdna} --bRNA {brna}" not in r["build"]:
                problems.append(f"{a}: its build command does not give --bDNA {bdna} --bRNA {brna}")
            if f"--index genome_library/{a}" not in r["publish"]:
                problems.append(f"{a}: its publish command does not publish genome_library/{a}")
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--commands", action="store_true",
                    help="print the build commands, then the publish commands, in order")
    args = ap.parse_args(argv)
    rows = read_manifest(args.manifest)
    if args.commands:
        for what in ("build", "publish"):
            print(f"# {what}")
            for r in rows:
                print(f"{r[what]}    # {r['archive']}")
        return 0
    problems = check(rows)
    for p in problems:
        print(f"PROBLEM: {p}")
    n_idx = sum(r["kind"] in INDEX_KINDS for r in rows)
    print(f"index set {os.path.basename(args.manifest)}: {len(rows)} archives "
          f"({n_idx} indexes), {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())

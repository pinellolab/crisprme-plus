#!/usr/bin/env python3
"""What a CRISPRme run records about its own inputs.

A results folder carries the sidecars the search wrote: ``.Params.txt`` (the genome, the
PAM, the budget it ran with), ``.command_line.txt`` (the command as it was given) and
``.sampleID.txt`` (the samples it read, gathered from the samplesID files of its
datasets). Reading them is how a tool handed a results folder alone can find the inputs
the run used, instead of asking the user to remember them -- or, worse, guessing one.

Shared so the two readers cannot drift: ``scripts/check_indel_rows.py``, the independent
check of CRISPRme's indel rows, which REFUSES a run whose records contradict what the
check was given, and ``PostProcess/extract_flanks.py`` (``crisprme.py extract-flanks``),
which resolves its inputs from the records where it is not given them. Standard library
only, as the checker is.
"""

import gzip
import os


def run_params(resdir):
    """The run's .Params.txt as {key: value} (empty when it has none)."""
    out = {}
    for name in (".Params.txt", "Params.txt"):
        p = os.path.join(resdir, name)
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) == 2:
                        out[f[0]] = f[1]
            break
    return out


def read_panel_meta(paths):
    """``{sample: {"population", "superpopulation", "sex"}}`` from samplesID files.

    ADDITIVE SIBLING of ``read_panel``, which returns a bare set and has callers that do
    set arithmetic on it -- changing that type to carry metadata would be a gratuitous
    risk, so this is a separate reader.

    WHY IT EXISTS. A samplesID file is a 4-column SAMPLE_ID / POPULATION_ID /
    SUPERPOPULATION_ID / SEX table written for exactly the panel a run searched, and every
    reader in the delivery path throws columns 2-4 away: the report keeps only the
    superpopulation, for one embedded chart, and ``read_panel`` keeps only column 1. So a
    run knows the ancestry of every individual it names and tells the reader none of it --
    which is the column an ancestry-diverse claim, or picking a representative carrier,
    actually rests on.

    A missing field is "" rather than a guess: the files are not guaranteed to carry all
    four columns, and an invented label is worse than an absent one.
    """
    out = {}
    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                sid = parts[0]
                if sid == "SAMPLE_ID":
                    continue
                out[sid] = {
                    "population": parts[1] if len(parts) > 1 else "",
                    "superpopulation": parts[2] if len(parts) > 2 else "",
                    "sex": parts[3] if len(parts) > 3 else "",
                }
    return out


def read_panel(paths):
    """The sample IDs of samplesID files (the first column of each line; '#' lines and a
    SAMPLE_ID header are skipped)."""
    out = set()
    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                sid = s.split()[0]
                if sid != "SAMPLE_ID":
                    out.add(sid)
    return out

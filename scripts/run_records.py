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

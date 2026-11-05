"""What is built on a genome's name, which a genome put in under that name must not replace.

A search, and a build of an index, take a genome's indexes (``genome_library/<PAM>_<N>_<name>``)
and its variant genomes (``Genomes/<name>+<dataset>``) by the genome folder's name alone, and
use each one that is there as it is. A genome put into ``Genomes/<name>`` over another one of
that name -- uploaded again with its sequences swapped, or downloaded from another URL -- was
then searched with the sequences and positions of the genome those were built from, at exit 0
and with no message. So the web Settings page refuses an upload, and ``crisprme.py download
--what genome --source url`` (the Settings page's Direct URL) a download, under the name of a
genome something installed is built on. ``--source ucsc`` and the HuggingFace download give
the same assembly by name, and are not refused. Both callers take the rule from here:
``pages/settings_page.py`` (``_built_on_genome``, ``_genome_upload_problem``) and
``crisprme.py`` ``download_data``.
"""

import os
from typing import List, Optional


def built_on(cwd: str, name: str) -> List[str]:
    """What is installed under ``cwd`` built on a genome named ``name``, as paths under it:
    its indexes (genome_library/<PAM>_<N>_<name>[+<dataset>], with their _INDELS companions)
    and its variant genomes (Genomes/<name>+<dataset>[_INDELS]). Hidden folders (a build in
    progress, an index set aside) are left out."""
    found = []
    for top in ("genome_library", "Genomes"):
        root = os.path.join(cwd, top)
        try:
            names = sorted(os.listdir(root))
        except OSError:
            continue
        for d in names:
            if top == "genome_library":  # (as the Settings page reads an index name)
                built = len(d.split("_", 2)) == 3 and d.split("_", 2)[2].split("+")[0] == name
            else:
                built = d.startswith(name + "+")
            if built and not d.startswith(".") and os.path.isdir(os.path.join(root, d)):
                found.append(f"{top}/{d}")
    return found


def replace_problem(cwd: str, name: str, how: str = "upload",
                    web: bool = True) -> Optional[str]:
    """Why a genome must not be put in as ``Genomes/<name>`` (``how``: "upload", or
    "download" for a Direct URL download; ``web``: said on the Settings page, else by
    crisprme.py download): what is listed by ``built_on``. None when nothing is built on it."""
    built = built_on(cwd, name)
    if not built:
        return None
    shown = ", ".join(built[:4]) + (", ..." if len(built) > 4 else "")
    if web:
        remove = (f"an index under 'Remove installed data'; a Genomes/{name}+<dataset> folder "
                  "from your data folder")
        again = f"{how} this genome under another name"
    else:
        remove = (f"each of them from {cwd.rstrip(os.sep) or os.sep}, or an index under 'Remove "
                  "installed data' in the web Settings page")
        again = "download it under another name (--ref)"
    return (
        f"Indexes or variant genomes built on a genome named {name!r} are installed "
        f"({shown}): searches of a genome {how}ed under that name would go on using "
        "them, with the sequences and positions of the genome they were built from. "
        f"Remove them first ({remove}), or {again}."
    )

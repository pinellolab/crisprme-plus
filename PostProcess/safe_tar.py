"""Unpacking a tar archive with safe modes: no setuid, setgid or sticky bit, and none of the
archive's owners.

A tar archive keeps each member's mode bits and owner, and tarfile's extractall, which trusts
an archive fully on Python 3.11, applies them as they are: run as root, as the Docker image runs
CRISPRme, it gives each file the owner the archive names, and keeps a setuid or setgid bit. The
web app unpacks into the data folder archives a request brings (a genome .tar.gz uploaded in
Settings or fetched from its Direct URL, a personal-assembly bundle) and archives it downloads
(an index, its genotype store, a build-inputs kit), and ``docker run -v`` mounts that folder
from the host: a setuid file put there by a request runs, from the host, with its owner's
rights, root's too. Each of them unpacks through ``extract``. A zip archive needs nothing of
the kind: zipfile sets no mode bits and no owner.
"""

import copy
import tarfile
from typing import List, Optional

# tarfile's 'data' filter (Python 3.11.4 and later)
_DATA_FILTER = getattr(tarfile, "data_filter", None)


def _safe_mode(member: tarfile.TarInfo) -> int:
    """The mode the 'data' filter gives ``member``, a regular file: no setuid, setgid or sticky
    bit, no write for the group or others, read and write for the owner, and execution for
    nobody unless for the owner; a folder, the same with all of the owner's."""
    mode = member.mode & 0o755
    if member.isdir():
        return mode | 0o700
    if not mode & 0o100:
        mode &= ~0o111
    return mode | 0o600


def extract(tar: tarfile.TarFile, dest: str,
            members: Optional[List[tarfile.TarInfo]] = None) -> None:
    """Unpacks ``members`` of ``tar`` (every member by default) into ``dest`` with safe modes,
    owned by the user unpacking it: through tarfile's 'data' filter where Python has it, which
    also refuses (tarfile.FilterError) a member with an absolute path or one leaving ``dest``, a
    link leading out of it, and a device or FIFO. Before 3.11.4 each member gets the mode that
    filter gives it, and no owner (an owner of -1 leaves the file's as it is), by hand: each
    caller unpacking a request's archive has refused any member but a regular file or a folder
    inside ``dest`` first."""
    if _DATA_FILTER is not None:
        tar.extractall(dest, members=members, filter="data")
        return
    safe = []
    for member in tar.getmembers() if members is None else members:
        member = copy.copy(member)
        member.mode = _safe_mode(member)
        member.uid = member.gid = -1
        member.uname = member.gname = ""
        safe.append(member)
    tar.extractall(dest, members=safe)

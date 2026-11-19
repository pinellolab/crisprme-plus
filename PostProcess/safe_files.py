"""Reading, writing and serving the data folder's files without following a link that someone
else put there.

The web app reads, serves and writes files in its data folder (``Results/``, ``Settings/`` and
the job folders in them, the email settings, the annotation manifest, ...) on paths a request
names (a job id, a file to download, a sample) in folders that every member of the folder's
group can write under a umask of 002. Such a member could put a link there, to a file of the
user's elsewhere, and the web app, which listens on every network interface, then read it out
to anyone who could reach it, serve it for download, or write through it (the review of lane S,
round 5, S-R5-2). Each such file goes through this module:

- every folder below the data folder is opened in turn, relative to the one before it
  (``openat``), with ``O_NOFOLLOW``: no link is followed at any step, and none can be put in
  the way between a check and the open;
- a file read is opened with ``O_NOFOLLOW`` and ``O_NONBLOCK`` and must be a regular file
  (``fstat``): not a link, a FIFO (which would hold the request) or a device;
- a file written whole is made anew under a name of this call's own
  (``O_CREAT|O_EXCL|O_NOFOLLOW``) in its folder and renamed over its name: a link there is
  replaced, never written through; a file that grows (a job's log) is opened with
  ``O_APPEND|O_NOFOLLOW``, and one written from its start with ``O_NOFOLLOW`` and truncated
  only once it is known to be a regular file;
- a file with another name (a hard link) is refused, read or written, in a folder others can
  write: a hard link to a file of the user's can be made where the system lets anyone link a
  file they cannot read (``fs.protected_hardlinks=0``);
- a folder the web app makes there for files of its own (a program's outputs, an upload
  unpacked) is used through its descriptor alone (``private_dir``), never by its name again:
  whoever can write the folder that holds it could rename it away and put a link in its
  place.

``links``, the rule for a link below the data folder's own folders: ``"none"`` refuses every
one; ``"own"`` follows one only where no one else can have put it: a link owned by the user
(the web app's user, the data folder's owner or root) in a folder that only such a user can
write, or in a sticky folder one of them owns. By default (``None``) the rule goes by the
folder: ``"none"`` in ``Results/`` and ``Settings/`` (``STRICT_FOLDERS``: the jobs' folders,
where CRISPRme makes and supports no link), ``"own"`` in the others, where links are a layout
it supports (a genome linked to a shared copy, a VCF dataset registered from a server folder,
an index's alias). A name in the data folder itself (``Results/`` or ``Settings/`` linked onto
a volume) always goes by ``"own"``. In a folder the group can write, a link is refused by
either rule.

Paths are given whole, with the data folder they must be under (``root``), which is trusted as
it is: a path that leaves it is refused. Every refusal is a ``RefusedPath``, an ``OSError``, so
a caller that takes an unreadable file in its stride takes a refused one the same way.
Standard library only: the light unit-test job imports it without Dash.
"""

import contextlib
import errno
import fcntl
import fnmatch
import os
import secrets
import shutil
import stat
import zipfile
from typing import IO, Iterator, List, NamedTuple, Optional

_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_NOFOLLOW = os.O_NOFOLLOW
_DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | _CLOEXEC
# the errno an open with O_NOFOLLOW gives on a link (FreeBSD: EMLINK), or with O_DIRECTORY on a
# link or a file
_LINK_ERRNOS = (errno.ELOOP, errno.EMLINK, errno.ENOTDIR)
LINKS_NONE = "none"
LINKS_OWN = "own"
# the data folder's folders where no link is a layout CRISPRme makes or supports: the jobs'
# folders (a search's, a Settings job's) and what is in them
STRICT_FOLDERS = ("Results", "Settings")


class RefusedPath(OSError):
    """A path refused: a link at some step below the data folder, a file that is not a regular
    one (or a folder that is not one), a file with another name, or a path outside the data
    folder."""

    def __init__(self, path: str, why: str):
        super().__init__(errno.ELOOP, why, path)
        self.why = why

    def __str__(self) -> str:
        return f"{self.filename}: {self.why}"


def job_root(job_folder: str) -> str:
    """The data folder of a job folder (``<data>/Results/<job>``, ``<data>/Settings/<job>``):
    the ``root`` to open its files under, for a function given the job folder alone."""
    return os.path.dirname(os.path.dirname(os.path.abspath(job_folder)))


def _split(path: str, root: str) -> List[str]:
    """The names of ``path``'s steps below ``root``: refused when it is ``root`` itself, leaves
    it, or has a ``..`` step anywhere. Taken as written, not normalised: a name a request
    sends, ``x/../../y``, never goes up through a folder ``x`` that is not there, as it would
    once ``..`` were taken out of it (empty and ``.`` steps are dropped)."""
    if not os.path.isabs(path):
        path = os.path.join(os.getcwd(), path)
    steps = [s for s in path.split(os.sep) if s not in ("", os.curdir)]
    top = [s for s in os.path.abspath(root).split(os.sep) if s]
    if os.pardir in steps or len(steps) <= len(top) or steps[:len(top)] != top:
        raise RefusedPath(path, "not a path in the data folder")
    return steps[len(top):]


def _trusted(root_fd: int) -> set:
    """The users whose links ``"own"`` follows: the web app's, the data folder's owner, root."""
    return {os.geteuid(), os.fstat(root_fd).st_uid, 0}


def _only_theirs(dir_fd: int, trusted: set) -> bool:
    """Whether no one but the users ``trusted`` can put or replace a name in the folder
    ``dir_fd``: one of them owns it, and the group and others cannot write it, or it is
    sticky (no one removes or renames a name another user owns)."""
    st = os.fstat(dir_fd)
    if st.st_uid not in trusted:
        return False
    return not st.st_mode & 0o022 or bool(st.st_mode & stat.S_ISVTX)


def _link_followable(dir_fd: int, name: str, links: str, trusted: set) -> bool:
    """Whether the link ``name`` in the folder ``dir_fd`` may be followed (``links`` "own")."""
    if links != LINKS_OWN:
        return False
    st = os.stat(name, dir_fd=dir_fd, follow_symlinks=False)
    return st.st_uid in trusted and _only_theirs(dir_fd, trusted)


def _is_link(dir_fd: int, name: str) -> bool:
    try:
        return stat.S_ISLNK(os.stat(name, dir_fd=dir_fd, follow_symlinks=False).st_mode)
    except OSError:
        return False


class _Walk:
    """The folder holding ``path``'s last step, opened step by step from ``root`` (``fd``), and
    that step's name (``name``); ``trusted``, the users of ``_trusted``; ``followed``, whether a
    link was followed on the way. A context manager that closes the folder."""

    def __init__(self, path: str, root: str, links: Optional[str] = None, *, whole: bool = False):
        parts = _split(path, root)
        if links is None:
            links = LINKS_NONE if parts[0] in STRICT_FOLDERS else LINKS_OWN
        if links not in (LINKS_NONE, LINKS_OWN):
            raise ValueError(f"links: {links!r}")
        self.path, self.links = path, links
        steps, self.name = (parts, None) if whole else (parts[:-1], parts[-1])
        fd = os.open(os.path.abspath(root), _DIR_FLAGS)
        try:
            self.trusted = _trusted(fd)
            self.followed = False
            self.depth = 0  # of the folder fd, below root
            for step in steps:
                if step in ("", os.curdir, os.pardir):
                    raise RefusedPath(path, "not a path in the data folder")
                nfd = self._step(fd, step)
                os.close(fd)
                fd = nfd
                self.depth += 1
        except BaseException:
            os.close(fd)
            raise
        self.fd = fd

    def rule(self) -> str:
        """The rule for a link in the folder ``fd``: "own" in the data folder itself."""
        return LINKS_OWN if self.depth == 0 else self.links

    def _step(self, fd: int, step: str) -> int:
        try:
            return os.open(step, _DIR_FLAGS | _NOFOLLOW, dir_fd=fd)
        except OSError as e:
            if e.errno not in _LINK_ERRNOS or not _is_link(fd, step):
                raise
            if not _link_followable(fd, step, self.rule(), self.trusted):
                raise RefusedPath(self.path, f"{step} is a link, which the web app does not "
                                             "follow there") from None
            self.followed = True
            return os.open(step, _DIR_FLAGS, dir_fd=fd)

    def close(self) -> None:
        if self.fd >= 0:
            os.close(self.fd)
            self.fd = -1

    def __enter__(self) -> "_Walk":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def _check_regular(fd: int, walk: _Walk, followed: bool) -> os.stat_result:
    """``fd``'s stat when it is a regular file with no other name where that matters, else
    RefusedPath (``fd`` is then closed)."""
    st = os.fstat(fd)
    why = None
    if not stat.S_ISREG(st.st_mode):
        why = "not a regular file"
    elif st.st_nlink > 1 and not followed and not _only_theirs(walk.fd, walk.trusted):
        why = "a file with another name (a hard link), in a folder others can write"
    if why:
        os.close(fd)
        raise RefusedPath(walk.path, why)
    return st


def _open_leaf(walk: _Walk, flags: int, perm: int = 0o666, follow: bool = True) -> int:
    """The last step opened with ``flags`` and O_NOFOLLOW, a regular file (see
    _check_regular); a link there is followed only when ``follow`` and the walk's rule let
    it."""
    followed = False
    try:
        fd = os.open(walk.name, flags | _NOFOLLOW | _CLOEXEC, perm, dir_fd=walk.fd)
    except OSError as e:
        if e.errno not in _LINK_ERRNOS or not _is_link(walk.fd, walk.name):
            raise
        if not (follow and _link_followable(walk.fd, walk.name, walk.rule(), walk.trusted)):
            raise RefusedPath(walk.path, f"{walk.name} is a link, which the web app does not "
                                         "follow there") from None
        fd = os.open(walk.name, flags | _CLOEXEC, perm, dir_fd=walk.fd)
        followed = True
    _check_regular(fd, walk, followed or walk.followed)
    if flags & os.O_NONBLOCK:
        fcntl.fcntl(fd, fcntl.F_SETFL, fcntl.fcntl(fd, fcntl.F_GETFL) & ~os.O_NONBLOCK)
    return fd


def _fdopen(fd: int, mode: str, encoding: Optional[str], errors: Optional[str],
            newline: Optional[str]) -> IO:
    try:
        if "b" in mode:
            return os.fdopen(fd, mode)
        return os.fdopen(fd, mode, encoding=encoding, errors=errors, newline=newline)
    except BaseException:
        os.close(fd)
        raise


# ---------------------------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------------------------
def open_fd(path: str, root: str, links: Optional[str] = None) -> int:
    """A read-only descriptor of the regular file ``path`` (the caller closes it)."""
    with _Walk(path, root, links) as walk:
        return _open_leaf(walk, os.O_RDONLY | os.O_NONBLOCK)


def open_read(path: str, root: str, mode: str = "r", *, links: Optional[str] = None,
              encoding: Optional[str] = None, errors: Optional[str] = None,
              newline: Optional[str] = None) -> IO:
    """The regular file ``path`` opened to read (``mode`` "r" or "rb")."""
    if mode not in ("r", "rb", "rt"):
        raise ValueError(f"open_read: mode {mode!r}")
    return _fdopen(open_fd(path, root, links), mode, encoding, errors, newline)


def read_text(path: str, root: str, *, links: Optional[str] = None, errors: Optional[str] = None,
              encoding: Optional[str] = None) -> str:
    with open_read(path, root, links=links, errors=errors, encoding=encoding) as fh:
        return fh.read()


def read_bytes(path: str, root: str, *, links: Optional[str] = None) -> bytes:
    with open_read(path, root, "rb", links=links) as fh:
        return fh.read()


def stat_file(path: str, root: str, links: Optional[str] = None) -> os.stat_result:
    """The stat of the regular file ``path`` (``os.stat``'s, with the rules above), without
    opening it."""
    with _Walk(path, root, links) as walk:
        st = os.stat(walk.name, dir_fd=walk.fd, follow_symlinks=False)
        if stat.S_ISLNK(st.st_mode):
            if not _link_followable(walk.fd, walk.name, walk.rule(), walk.trusted):
                raise RefusedPath(path, f"{walk.name} is a link, which the web app does not "
                                        "follow there")
            st = os.stat(walk.name, dir_fd=walk.fd)
        elif (stat.S_ISREG(st.st_mode) and st.st_nlink > 1 and not walk.followed
              and not _only_theirs(walk.fd, walk.trusted)):
            raise RefusedPath(path, "a file with another name (a hard link), in a folder "
                                    "others can write")
        if not stat.S_ISREG(st.st_mode):
            raise RefusedPath(path, "not a regular file")
        return st


def getsize(path: str, root: str, links: Optional[str] = None) -> int:
    return stat_file(path, root, links).st_size


def is_file(path: str, root: str, links: Optional[str] = None) -> bool:
    """``os.path.isfile`` with the rules above: False for a refused path, never raises."""
    try:
        stat_file(path, root, links)
    except (OSError, ValueError):
        return False
    return True


def is_dir(path: str, root: str, links: Optional[str] = None) -> bool:
    """``os.path.isdir`` with the rules above (every step of ``path`` a folder reached so)."""
    try:
        _Walk(path, root, links, whole=True).close()
    except (OSError, ValueError):
        return False
    return True


def check_dir(path: str, root: str, links: Optional[str] = None) -> None:
    """Raises (RefusedPath, or the OSError of a missing step) unless ``path`` is a folder
    reached with the rules above."""
    _Walk(path, root, links, whole=True).close()


def listdir(path: str, root: str, links: Optional[str] = None) -> List[str]:
    """The names in the folder ``path`` (``os.listdir``'s, in its order)."""
    with _Walk(path, root, links, whole=True) as walk:
        return os.listdir(walk.fd)


def glob(folder: str, pattern: str, root: str, links: Optional[str] = None) -> List[str]:
    """``glob.glob(os.path.join(folder, pattern))`` for a pattern of one name, in the folder
    ``folder`` reached with the rules above, sorted (glob gives the folder's order); a pattern
    that does not start with a dot leaves out the hidden names, as glob does. [] when the
    folder is not there or is refused."""
    try:
        names = listdir(folder, root, links)
    except OSError:
        return []
    hidden = pattern.startswith(".")
    return [os.path.join(folder, n) for n in sorted(names)
            if fnmatch.fnmatchcase(n, pattern) and (hidden or not n.startswith("."))]


# ---------------------------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------------------------
def _temp_name(name: str) -> str:
    return f".{name}.{secrets.token_hex(8)}.tmp"


@contextlib.contextmanager
def replace(path: str, root: str, mode: str = "w", *, perm: int = 0o666,
            encoding: Optional[str] = None, newline: Optional[str] = None) -> Iterator[IO]:
    """A new file to write ``path``'s content into, renamed over ``path`` when the block ends
    without an exception (removed when it raises): made under a name of its own in ``path``'s
    folder with O_CREAT|O_EXCL|O_NOFOLLOW and mode ``perm`` (less the umask), so a link at
    ``path`` is replaced, never written through, and the file is the user's whatever was
    there."""
    if mode not in ("w", "wb", "wt"):
        raise ValueError(f"replace: mode {mode!r}")
    with _Walk(path, root) as walk:
        tmp = _temp_name(walk.name)
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW | _CLOEXEC, perm,
                     dir_fd=walk.fd)
        try:
            with _fdopen(fd, mode, encoding, None, newline) as fh:
                yield fh
            os.replace(tmp, walk.name, src_dir_fd=walk.fd, dst_dir_fd=walk.fd)
        except BaseException:
            with contextlib.suppress(OSError):
                os.unlink(tmp, dir_fd=walk.fd)
            raise


def write_text(path: str, root: str, text: str, *, perm: int = 0o666,
               encoding: Optional[str] = None) -> None:
    with replace(path, root, perm=perm, encoding=encoding) as fh:
        fh.write(text)


def write_bytes(path: str, root: str, data: bytes, *, perm: int = 0o666) -> None:
    with replace(path, root, "wb", perm=perm) as fh:
        fh.write(data)


def open_append(path: str, root: str, mode: str = "a", *, perm: int = 0o666,
                encoding: Optional[str] = None, errors: Optional[str] = None) -> IO:
    """``path`` opened to add to its end (made when it is not there): O_APPEND|O_NOFOLLOW, a
    regular file with no other name (a link at ``path`` is refused)."""
    if mode not in ("a", "ab", "at"):
        raise ValueError(f"open_append: mode {mode!r}")
    with _Walk(path, root) as walk:
        fd = _open_leaf(walk, os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_NONBLOCK, perm,
                        follow=False)
    return _fdopen(fd, mode, encoding, errors, None)


def open_write(path: str, root: str, mode: str = "w", *, perm: int = 0o666,
               encoding: Optional[str] = None, errors: Optional[str] = None) -> IO:
    """``path`` opened to write from its start (made when it is not there), for a file
    another process writes as it goes (a job's log): O_NOFOLLOW without O_TRUNC, then
    truncated once it is known to be a regular file with no other name (a link at ``path``
    is refused, and the file a hard link names is never truncated)."""
    if mode not in ("w", "wb", "wt"):
        raise ValueError(f"open_write: mode {mode!r}")
    with _Walk(path, root) as walk:
        fd = _open_leaf(walk, os.O_WRONLY | os.O_CREAT | os.O_NONBLOCK, perm, follow=False)
    try:
        os.ftruncate(fd, 0)
    except BaseException:
        os.close(fd)
        raise
    return _fdopen(fd, mode, encoding, errors, None)


def create(path: str, root: str, mode: str = "w", *, perm: int = 0o666,
           encoding: Optional[str] = None) -> IO:
    """``path`` made anew (O_CREAT|O_EXCL|O_NOFOLLOW): FileExistsError when anything, a link
    too, is there."""
    if mode not in ("w", "wb", "wt"):
        raise ValueError(f"create: mode {mode!r}")
    with _Walk(path, root) as walk:
        fd = os.open(walk.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW | _CLOEXEC,
                     perm, dir_fd=walk.fd)
    return _fdopen(fd, mode, encoding, None, None)


def touch(path: str, root: str) -> None:
    """Makes ``path`` when it is not there and sets its times to now (a link refused)."""
    with open_append(path, root, "ab") as fh:
        os.utime(fh.fileno())


def mkdir(path: str, root: str, perm: int = 0o777, *, exist_ok: bool = False) -> None:
    """``os.mkdir`` in a folder reached with the rules above; with ``exist_ok``, one there
    already must be a folder reached so too (a link there only by the rule for the data
    folder's own folders)."""
    with _Walk(path, root) as walk:
        try:
            os.mkdir(walk.name, perm, dir_fd=walk.fd)
            return
        except FileExistsError:
            if not exist_ok:
                raise
    check_dir(path, root)


def makedirs(path: str, root: str, perm: int = 0o777) -> None:
    """``os.makedirs(path, exist_ok=True)`` below ``root``, each step with the rules above."""
    parts = _split(path, root)
    for i in range(1, len(parts) + 1):
        mkdir(os.path.join(os.path.abspath(root), *parts[:i]), root, perm, exist_ok=True)


def remove(path: str, root: str, links: Optional[str] = None) -> None:
    """Removes the name ``path`` (a link itself, never what it names) in a folder reached with
    the rules above."""
    with _Walk(path, root, links) as walk:
        os.unlink(walk.name, dir_fd=walk.fd)


def rmtree(path: str, root: str, links: Optional[str] = None, *,
           ignore_errors: bool = False) -> None:
    """``shutil.rmtree`` of the folder ``path`` within the folder that holds it, reached with
    the rules above (a link at ``path`` itself is refused, whatever ``links`` says): no
    folder another user renames into a step of ``path`` meanwhile is removed in its
    place."""
    try:
        with _Walk(path, root, links) as walk:
            st = os.stat(walk.name, dir_fd=walk.fd, follow_symlinks=False)
            if not stat.S_ISDIR(st.st_mode):
                raise RefusedPath(path, "not a folder")
            shutil.rmtree(walk.name, ignore_errors=ignore_errors, dir_fd=walk.fd)
    except OSError:
        if not ignore_errors:
            raise


def rename(src: str, dst: str, root: str, *, dst_links: Optional[str] = None) -> None:
    """Renames the regular file ``src`` over ``dst``, each in a folder reached with the rules
    above (``dst_links``: ``dst``'s rule), as ``os.replace`` does: a link at ``dst`` is
    replaced, and ``src`` must be a regular file with no other name (a link there is
    refused). One put in place of ``src`` while it is moved is taken back out of ``dst``."""
    with _Walk(src, root) as swalk, _Walk(dst, root, dst_links) as dwalk:
        _move_regular(swalk, src, dwalk.fd, dwalk.name)


def _move_regular(swalk: _Walk, src: str, dst_fd: int, dst_name: str) -> None:
    """``rename``'s move of the file ``swalk`` reached to ``dst_name`` in the folder
    ``dst_fd``."""
    fd = _open_leaf(swalk, os.O_RDONLY | os.O_NONBLOCK, follow=False)
    try:
        st = os.fstat(fd)
        os.replace(swalk.name, dst_name, src_dir_fd=swalk.fd, dst_dir_fd=dst_fd)
        now = os.stat(dst_name, dir_fd=dst_fd, follow_symlinks=False)
        if (now.st_dev, now.st_ino) != (st.st_dev, st.st_ino):
            with contextlib.suppress(OSError):
                os.unlink(dst_name, dir_fd=dst_fd)
            raise RefusedPath(src, "replaced by another file while it was moved")
    finally:
        os.close(fd)


def move_into(src: str, root: str, dst: str) -> None:
    """Renames the regular file ``src`` (as ``rename`` takes it) to ``dst``, a path in a
    folder of this process's own (``private_dir``'s ``path``, or a folder made below it),
    which is taken as it is: no one else can put a name there (and ``path`` is itself a
    link, to the folder's descriptor)."""
    folder = os.open(os.path.dirname(dst), _DIR_FLAGS)
    try:
        with _Walk(src, root) as swalk:
            _move_regular(swalk, src, folder, os.path.basename(dst))
    finally:
        os.close(folder)


def replace_from(src: str, dst: str, root: str) -> None:
    """Renames the file ``src``, a path of a folder the caller made and only it writes (a
    private temporary folder), over ``dst``, whose folder is reached with the rules above:
    a link at ``dst`` is replaced, never written through. Copied there (``replace``) when
    the two are on different file systems."""
    with _Walk(dst, root) as walk:
        try:
            os.replace(src, walk.name, dst_dir_fd=walk.fd)
            return
        except OSError as e:
            if e.errno != errno.EXDEV:
                raise
    with open(src, "rb") as fin, replace(dst, root, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.unlink(src)


# ---------------------------------------------------------------------------------------------
# what other code opens by name
# ---------------------------------------------------------------------------------------------
def fd_path(fd: int) -> str:
    """A name that opens the very file ``fd`` is open on, whatever is renamed meanwhile (Linux's
    /proc/self/fd, else /dev/fd), for a library or a program that takes a file name; a child
    process must inherit ``fd`` (subprocess's ``pass_fds``)."""
    proc = f"/proc/self/fd/{fd}"
    return proc if os.path.isdir("/proc/self/fd") else f"/dev/fd/{fd}"


@contextlib.contextmanager
def opened_path(path: str, root: str, links: Optional[str] = None) -> Iterator[str]:
    """The regular file ``path``, opened with the rules above, as a name (``fd_path``) for a
    library that takes a file name: it opens that very file, whatever is put at ``path``
    meanwhile. The descriptor is closed when the block ends."""
    fd = open_fd(path, root, links)
    try:
        yield fd_path(fd)
    finally:
        os.close(fd)


class PrivateDir(NamedTuple):
    """A folder of ``private_dir``'s: ``path``, the name of its descriptor
    (``/proc/self/fd/<fd>``), and ``fd``, which a program given a path below it must
    inherit (subprocess's ``pass_fds=(d.fd,)``)."""
    path: str
    fd: int


def _empty_dir(fd: int) -> None:
    """Removes everything in the folder ``fd`` (each name in it, never what a link there
    names), through its descriptor."""
    try:
        names = os.listdir(fd)
    except OSError:
        return
    for name in names:
        with contextlib.suppress(OSError):
            if stat.S_ISDIR(os.stat(name, dir_fd=fd, follow_symlinks=False).st_mode):
                shutil.rmtree(name, ignore_errors=True, dir_fd=fd)
            else:
                os.unlink(name, dir_fd=fd)


@contextlib.contextmanager
def private_dir(folder: str, root: str, prefix: str = ".crisprme-", *,
                links: Optional[str] = None) -> Iterator[PrivateDir]:
    """A folder of this process's own, made in the folder ``folder`` (reached with the rules
    above) under a name no one could make first, mode 0700, and used through its descriptor
    alone: every path below ``path`` (``/proc/self/fd/<fd>``) reaches the very folder made
    here, whatever is renamed or put in ``folder`` meanwhile. Removed, with what is in it,
    when the block ends.

    ``folder`` may be one others can write (a job's folder, ``Annotations/``, ``Genomes/``,
    ``Assemblies/``): a folder made there by ``tempfile.mkdtemp`` and then used by its name
    could be renamed away within ``folder`` by anyone who can write it (a rename there needs
    no more), and a link to a folder of the user's put under its name, which the web app
    then wrote into, unpacked an upload into, and moved files out of by that name. A
    program given a path below it inherits ``fd`` (``pass_fds``). Linux only (there is no
    other way to reach a folder below a descriptor by name): elsewhere OSError."""
    proc = "/proc/self/fd"
    if not os.path.isdir(proc):
        raise OSError(errno.ENOTSUP, "no /proc/self/fd here, to reach a folder by its "
                                     "descriptor", folder)
    with _Walk(folder, root, links, whole=True) as walk:
        for _ in range(100):
            name = f"{prefix}{secrets.token_hex(8)}"
            try:
                os.mkdir(name, 0o700, dir_fd=walk.fd)
                break
            except FileExistsError:
                continue
        else:  # pragma: no cover
            raise FileExistsError(errno.EEXIST, "no free name for a private folder", folder)
        try:
            fd = os.open(name, _DIR_FLAGS | _NOFOLLOW, dir_fd=walk.fd)
        except OSError as e:
            if e.errno in _LINK_ERRNOS:
                raise RefusedPath(os.path.join(folder, name), "replaced by a link as soon as "
                                                              "it was made") from None
            raise
        try:
            # the folder made here, not one another user put under its name meanwhile (a
            # folder of the user's that the group can write can be moved there, and one of
            # theirs, which root can write into): only an empty, private one, owned by
            # whoever a file this process makes in it is owned by, passes (another such
            # folder of this process's, then). Checked before anything in it is used, or
            # removed.
            st = os.fstat(fd)
            if not stat.S_ISDIR(st.st_mode) or st.st_mode & 0o077 or os.listdir(fd):
                raise RefusedPath(os.path.join(folder, name), "replaced by another folder as "
                                                              "soon as it was made")
            try:
                probe = os.open(".owner", os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW
                                | _CLOEXEC, 0o600, dir_fd=fd)
            except OSError:
                raise RefusedPath(os.path.join(folder, name), "replaced by another folder as "
                                                              "soon as it was made") from None
            try:
                mine = os.fstat(probe).st_uid
            finally:
                os.close(probe)
                os.unlink(".owner", dir_fd=fd)
            if st.st_uid != mine:
                raise RefusedPath(os.path.join(folder, name), "replaced by another folder as "
                                                              "soon as it was made")
            parent = os.dup(walk.fd)
        except BaseException:
            os.close(fd)
            raise
    try:
        yield PrivateDir(f"{proc}/{fd}", fd)
    finally:
        _empty_dir(fd)
        with contextlib.suppress(OSError):  # its name, when that is still this folder
            now, mine = os.stat(name, dir_fd=parent, follow_symlinks=False), os.fstat(fd)
            if (now.st_dev, now.st_ino) == (mine.st_dev, mine.st_ino):
                os.rmdir(name, dir_fd=parent)
        os.close(fd)
        os.close(parent)


def sqlite_connect(path: str, root: str):
    """A read-only connection (``mode=ro&immutable=1``) to the SQLite database ``path``, a
    regular file reached with the rules above: SQLite opens the database by name, so it is
    given the name of the descriptor opened here (``fd_path``)."""
    import sqlite3

    fd = open_fd(path, root)
    try:
        conn = sqlite3.connect(f"file:{fd_path(fd)}?mode=ro&immutable=1", uri=True)
        try:
            conn.execute("PRAGMA schema_version").fetchone()  # SQLite has opened it now
        except BaseException:
            conn.close()
            raise
    finally:
        os.close(fd)
    return conn


_COMPRESSION = {".gz": "gzip", ".bz2": "bz2", ".zip": "zip", ".xz": "xz", ".zst": "zstd"}


def read_csv(path: str, root: str, *, links: Optional[str] = None, **kwargs):
    """``pandas.read_csv(path, **kwargs)`` of the regular file ``path``, opened with the rules
    above (its compression, unless given, by its extension, as read_csv infers it from a
    name)."""
    import pandas as pd

    kwargs.setdefault("compression", _COMPRESSION.get(os.path.splitext(path)[1].lower()))
    with open_read(path, root, "rb", links=links) as fh:
        return pd.read_csv(fh, **kwargs)


def to_csv(frame, path: str, root: str, **kwargs) -> None:
    """``frame.to_csv(path, **kwargs)`` written through ``replace``."""
    with replace(path, root, newline="") as fh:
        frame.to_csv(fh, **kwargs)


def zip_file(src: str, dst: str, root: str, arcname: Optional[str] = None) -> None:
    """The zip archive ``dst`` holding the regular file ``src`` under ``arcname`` (its name by
    default), as ``zip -j dst src`` makes it, from a file opened with the rules above and
    written through ``replace`` (``zip`` follows a link at ``src``, and reads as an archive to
    add to whatever a link at ``dst`` names)."""
    import time

    with open_read(src, root, "rb") as fin:
        st = os.fstat(fin.fileno())
        stamp = time.localtime(max(st.st_mtime, 315532800))[:6]  # (zip dates start in 1980)
        info = zipfile.ZipInfo(arcname or os.path.basename(src), date_time=stamp)
        info.external_attr = (st.st_mode & 0xFFFF) << 16
        info.compress_type = zipfile.ZIP_DEFLATED
        info.file_size = st.st_size
        with replace(dst, root, "wb") as fout, zipfile.ZipFile(fout, "w") as zf, \
                zf.open(info, "w", force_zip64=st.st_size >= zipfile.ZIP64_LIMIT) as member:
            shutil.copyfileobj(fin, member, 1 << 20)


def send_file(path: str, root: str, *, as_attachment: bool = True,
              download_name: Optional[str] = None):
    """A Flask response serving the regular file ``path``, opened with the rules above
    (``flask.send_from_directory`` opens the name it joins, links and all)."""
    import flask

    fh = open_read(path, root, "rb")
    try:
        st = os.fstat(fh.fileno())
        response = flask.send_file(fh, as_attachment=as_attachment,
                                   download_name=download_name or os.path.basename(path),
                                   last_modified=st.st_mtime)
    except BaseException:
        fh.close()
        raise
    if response.status_code == 200:
        response.content_length = st.st_size
    return response

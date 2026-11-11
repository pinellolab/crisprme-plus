"""Small persisted web-app settings (currently: max concurrent search jobs).

The number of search jobs that run at once used to be hard-coded to 2 (an old
CRISPRme web-server constraint). It is now configurable so an operator can raise
it on a bigger machine — but the DEFAULT stays a conservative **2** so the public
server / a laptop is never accidentally oversubscribed.

Resolution order (first that applies wins):
  1. the ``CRISPRME_MAX_CONCURRENT_JOBS`` environment variable (operator pin at launch)
  2. the persisted value set from the web Settings page (``<workdir>/.crisprme_app_config.json``)
  3. the built-in default (2)

Values are clamped to [1, MAX_JOBS_CAP]. The executor is created once at startup,
so a change made in Settings takes effect on the next app restart.
"""

import json
import os
import tempfile

CONFIG_FILE = ".crisprme_app_config.json"
DEFAULT_MAX_JOBS = 2
MAX_JOBS_CAP = 16  # sane upper bound; a search job is memory-heavy (~tens of GB each)


def _config_path(workdir=None):
    base = workdir or os.environ.get("CRISPRME_CONFIG_DIR") or os.getcwd()
    return os.path.join(base, CONFIG_FILE)


def _clamp(n):
    return max(1, min(int(n), MAX_JOBS_CAP))


def env_override():
    """The max concurrent search jobs CRISPRME_MAX_CONCURRENT_JOBS pins (clamped),
    overriding the persisted value; None when it is unset or not a whole number in
    ASCII digits, which is ignored. Never raises."""
    env = os.environ.get("CRISPRME_MAX_CONCURRENT_JOBS", "").strip()
    # ASCII only: int() rejects some other str.isdigit() characters, such as
    # superscript or circled digits.
    if env.isascii() and env.isdigit():
        try:
            return _clamp(env)
        except ValueError:  # more digits than int() converts (4300 by default)
            pass
    return None


def get_max_concurrent_jobs(workdir=None, default=DEFAULT_MAX_JOBS):
    """Resolve the configured max concurrent search jobs (env > file > default)."""
    pinned = env_override()
    if pinned is not None:
        return pinned
    try:
        with open(_config_path(workdir)) as fh:
            return _clamp(json.load(fh).get("max_concurrent_jobs", default))
    except (OSError, ValueError, TypeError):
        return _clamp(default)


def set_max_concurrent_jobs(n, workdir=None):
    """Persist the max concurrent search jobs (clamped). Returns the stored value.

    Takes effect at the next app restart (the process pool is built at startup).
    """
    n = _clamp(n)
    path = _config_path(workdir)
    cfg = {}
    try:
        with open(path) as fh:
            cfg = json.load(fh)
    except (OSError, ValueError):
        cfg = {}
    if not isinstance(cfg, dict):
        cfg = {}
    cfg["max_concurrent_jobs"] = n
    # a temporary name of its own (mkstemp: O_EXCL, O_NOFOLLOW) renamed over the file: open()
    # wrote through a link another user put at its name (the web Settings page saves it in
    # the data folder, which a group may write)
    fd, tmp = tempfile.mkstemp(prefix=CONFIG_FILE + ".", suffix=".tmp",
                               dir=os.path.dirname(path) or ".")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(cfg, fh, indent=2)
        os.chmod(tmp, 0o644)  # (readable, as open() made it: mkstemp makes it 0600)
        os.replace(tmp, path)
    except BaseException:
        if os.path.lexists(tmp):
            os.remove(tmp)
        raise
    return n

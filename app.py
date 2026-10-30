"""Main module for the CRISPRme web application.

This module sets up the Flask and Dash applications, configures the server, and initializes the web application. It also defines stylesheets, caching, and various application settings.

Attributes:
    WEBADDRESS (str): The web address for accessing the application.
    IPADDRESS (str): The IP address and port for local access to the web application.
    URL (str): The server URL.
    external_stylesheets (list): List of external stylesheets for the web application.
    server (Flask): The Flask server instance.
    app (Dash): The Dash application instance.
    app_directory (str): The directory of the current application.
    current_working_directory (str): The current working directory.
    operators (list): List of filtering operators used for querying tables.
    ONLINE (bool): Flag indicating if the application is online or offline.
    DISPLAY_OFFLINE (str): CSS display property for offline mode.
    DISPLAY_ONLINE (str): CSS display property for online mode.
    pool_executor (ProcessPoolExecutor): Executor for running multiple jobs 
        concurrently.
    CACHE_CONFIG (dict): Configuration settings for caching.
    cache (Cache): Cache instance for the application.
"""

import concurrent.futures
import sys
import os

try:
    from flask_caching import Cache

    import dash_bootstrap_components as dbc

    import dash
    import flask
except ImportError as e:
    sys.stderr.write(
        f"Error: Web interface dependencies not installed: {e}\n"
        "Install them with: pip install crisprme[web]\n"
        "or: pip install dash flask flask-caching dash-bootstrap-components\n"
    )
    sys.exit(1)


# define website address -> DEPRECATED <-
WEBADDRESS = "http://crisprme.di.univr.it"

# define server IP address (local deployment)
IPADDRESS = "127.0.0.1:8080"

# define server URL
URL = ""

# define current working directory
current_working_directory = f"{os.getcwd()}/"


def start_message() -> None:
    """Prints a startup message to the standard error stream.

    This function outputs a message indicating that the server has started and
    provides the URL to access the web application.

    Args:
        None

    Returns:
        None
    """

    sys.stderr.write("SERVER STARTED\n")
    sys.stderr.write(f"GO TO http://{IPADDRESS} TO USE THE WEB APP\n\n")


# --> entry point <-- #
# The old chriddyp/bWLwgP "skeleton" stylesheet was dropped: modern browsers
# block it cross-origin (ERR_BLOCKED_BY_RESPONSE.NotSameOrigin), so it failed to
# load on every page and only produced a console error. The layout uses the
# dash-bootstrap-components grid (dbc.Row/dbc.Col), so Bootstrap is all we need.
external_stylesheets = [dbc.themes.BOOTSTRAP]  # CSS used to style the website
server = flask.Flask(__name__)  # define flask app.server
app = dash.Dash(
    __name__,
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True,
    server=server, # type: ignore
)  # initialize server app
app_directory = os.path.dirname(os.path.realpath(__file__))  # current location
start_message()  # print server start message
app.title = "CRISPRme"  # assign flask app name
# necessary if update element in a callback generated in another callback
# app.config['suppress_callback_exceptions'] = True
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True
# define filtering operators used when querying tables
operators = [
    ["ge ", ">="],
    ["le ", "<="],
    ["lt ", "<"],
    ["gt ", ">"],
    ["ne ", "!="],
    ["eq ", "="],
    ["contains "],
]
ONLINE = False  # NOTE change to True for online version, False for offline
DISPLAY_OFFLINE = ""
DISPLAY_ONLINE = ""
if ONLINE:
    DISPLAY_OFFLINE = "none"
else:
    DISPLAY_ONLINE = "none"
# Max concurrent search jobs. Configurable (env CRISPRME_MAX_CONCURRENT_JOBS or the
# web Settings page), but the DEFAULT stays a conservative 2 so a laptop / the public
# server is never oversubscribed; a search job is memory-heavy (~tens of GB each). The
# pool is built once at startup, so a Settings change applies on the next restart.
# app_settings lives in PostProcess/, which is NOT on sys.path yet: `crisprme.py
# web-interface` starts index.py as a fresh subprocess (sys.path[0] is just this
# directory), and index.py imports this module before pages/pages_utils.py inserts
# PostProcess/. Insert it here, from this file's own location like pages_utils
# does, so the setting is found whatever the cwd.
_POSTPROCESS_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "PostProcess"
)
if _POSTPROCESS_PATH not in sys.path:
    sys.path.insert(0, _POSTPROCESS_PATH)
try:
    import app_settings
    _MAX_CONCURRENT_JOBS = app_settings.get_max_concurrent_jobs()
except Exception as e:
    _MAX_CONCURRENT_JOBS = 2
    sys.stdout.write(
        "WARNING: cannot read the concurrent search jobs setting (app_settings: "
        f"{type(e).__name__}: {e}); using the default of {_MAX_CONCURRENT_JOBS} "
        "concurrent search jobs\n"
    )
    sys.stdout.flush()
else:
    # app_settings ignores a CRISPRME_MAX_CONCURRENT_JOBS that is not a whole number
    # (e.g. "5" with the quotes `docker run --env-file` keeps) and clamps one into
    # [1, MAX_JOBS_CAP]; say so rather than silently run another pool size.
    _env_jobs = os.environ.get("CRISPRME_MAX_CONCURRENT_JOBS", "").strip()
    if _env_jobs and not (_env_jobs.isdigit() and int(_env_jobs) == _MAX_CONCURRENT_JOBS):
        sys.stdout.write(
            f"WARNING: CRISPRME_MAX_CONCURRENT_JOBS={_env_jobs!r} is not a whole number "
            f"from 1 to {app_settings.MAX_JOBS_CAP}; using {_MAX_CONCURRENT_JOBS} "
            "concurrent search job(s)\n"
        )
        sys.stdout.flush()
pool_executor = concurrent.futures.ProcessPoolExecutor(max_workers=_MAX_CONCURRENT_JOBS)
# A separate single-slot executor for Settings / Data-Manager operations
# (genome/index/VCF downloads and local index builds). Keeping these off the
# search `pool_executor` means a long download or index build never starves the
# search slots (_MAX_CONCURRENT_JOBS of them: 1-16, 2 by default).
settings_executor = concurrent.futures.ProcessPoolExecutor(max_workers=1)
# Maintainer-only actions (e.g. publishing an index to the shared HuggingFace
# repo) are gated behind this flag: enabled only when a HuggingFace write token
# is configured or CRISPRME_MAINTAINER is set. Never exposed to normal users.
MAINTAINER_MODE = bool(
    os.environ.get("CRISPRME_MAINTAINER")
    or os.environ.get("HF_TOKEN")
    or os.environ.get("HUGGING_FACE_HUB_TOKEN")
)
# Guardrail for the small-file (annotation) browser uploads; large genomes/VCFs
# never go through an in-memory upload (server-side fetch / chunked upload only).
app.server.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # 200 MB
# configure caching
CACHE_CONFIG = {
    # "FileSystemCache" (a Flask-Caching backend class name) if you don't want to
    # set up redis. NB: use the class name, not the legacy lowercase "filesystem"
    # short-name -- Flask-Caching >=2.0 dropped the old short-name aliases, so
    # "filesystem" raises ImportStringError on the modernized web stack.
    "CACHE_TYPE": "FileSystemCache",
    "CACHE_DIR": ("Cache"),  # os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()  # initialize cache
cache.init_app(app.server, config=CACHE_CONFIG)  # start web-app

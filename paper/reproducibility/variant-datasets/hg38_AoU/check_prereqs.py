"""The Python packages retrieve_AoU.sh's steps import (requests, pandas, pysam), imported as
those steps import them: run as a script, so Python's path is that of the steps (this
folder, PYTHONPATH, the user's site-packages) and not the working directory, which
"python -c" would put first (a module there would run on its first import)."""

import pandas  # noqa: F401
import pysam  # noqa: F401
import requests  # noqa: F401

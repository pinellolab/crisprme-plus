"""Provides functions for displaying and managing CRISPRme job history.

This module implements the history page of the CRISPRme web application,
allowing users to view and interact with past job results. It includes
functions for retrieving job parameters, displaying a history table, and
handling user interactions such as row selection and filtering.
"""

from .pages_utils import RESULTS_DIR, PARAMS_FILE, LOG_FILE, GUIDES_FILE

from typing import Dict, List, Optional, Tuple
from dash import Input, Output, State
from dash.exceptions import PreventUpdate
from app import app, current_working_directory, URL

from dash import html
import pandas as pd
import numpy as np

from datetime import datetime

import math
import os
import pathlib

SUMMARYTABCOLS = [
    "Job", "Genome", "Variants", "Mismatches", "DNA bulge", "RNA bulge", "Max edits", "PAM", "Number of Guides", "Start"
]


def support_filter_history(
    result_df: pd.DataFrame, genome_filter: str, pam_filter: str
) -> Tuple[pd.DataFrame, int]:
    """Filter the results history dataframe based on genome and PAM criteria.

    This function filters the input dataframe to include only jobs matching the
    specified genome build and PAM sequence(s). It also calculates the maximum
    number of pages for displaying the filtered results.

    Args:
        result_df: The dataframe containing the job history.
        genome_filter: The genome build to filter by.
        pam_filter: The PAM sequence(s) to filter by.

    Returns:
        A tuple containing the filtered dataframe and the maximum number of pages.
    """
    if not isinstance(result_df, pd.DataFrame):
        raise TypeError(
            f"Expected {type(pd.DataFrame).__name__}, got {type(result_df).__name__}"
        )
    if not isinstance(genome_filter, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome_filter).__name__}")
    if not isinstance(pam_filter, str):
        raise TypeError(f"Exepcted {str.__name__}, got {type(pam_filter).__name__}")
    if genome_filter is not None:
        # keep only rows related to the requested genome type
        drop_rows = result_df[result_df["Genome"] != genome_filter].index
        result_df.drop(drop_rows, axis=0, inplace=True)
    if pam_filter is not None:
        drop_rows = result_df[~result_df["PAM"].isin(pam_filter)].index
        result_df.drop(drop_rows, axis=0, inplace=True)
    # allow also empty results history?
    max_page = result_df.shape[0]
    max_page = math.floor(max_page / 1000000) + 1  # max size
    return result_df, max_page


# trigger row highlighting
@app.callback(
    Output("results-table", "style_data_conditional"),
    [Input("results-table", "selected_cells")],
    [State("results-table", "data")],
)
def highlight_row(sel_cel: List, all_guides: List) -> List[Dict[str, str]]:
    """Highlight the selected row in the results table.

    This callback function dynamically styles the results table to highlight the
    row corresponding to the selected cell.

    Args:
        sel_cel: A list representing the selected cell(s).
        all_guides: The data of the results table.

    Returns:
        A list of dictionaries specifying the styling for the highlighted row.

    Raises:
        PreventUpdate: If no cell is selected or the table data is empty.
    """
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    job_name = all_guides[int(sel_cel[0]["row"])]["Job"]
    assert isinstance(job_name, str)
    return [
        {
            "if": {"filter_query": "".join(["{Job} eq ", f'"{job_name}"'])},
            "background-color": "rgba(0, 0, 255,0.15)",  # highlighting color
        }
    ]


def retrieve_resultsdirs(resultsdir: str) -> List[str]:
    """Retrieve a list of valid result directories.

    This function scans the specified results directory and returns a list of
    subdirectories that represent valid CRISPRme job results. A directory is
    considered valid if it contains a PARAMS_FILE.

    Args:
        resultsdir: The path to the results directory.

    Returns:
        A list of strings, where each string is the path to a valid result
        directory.
    """
    return [
        d
        for d in os.listdir(resultsdir)
        if os.path.isdir(os.path.join(resultsdir, d)) and os.path.isfile(os.path.join(resultsdir, d, PARAMS_FILE))
    ]

def read_params(paramsfile: str, jobid: str) -> Dict[str, str]:
    """Read job parameters from a file.

    This function reads the parameters used for a specific CRISPRme job from
    the given parameters file.

    Args:
        paramsfile: The path to the parameters file.
        jobid: The ID of the job.

    Returns:
        A dictionary containing the job parameters, where keys are parameter
        names and values are parameter values.

    Raises:
        IOError: If an error occurs while reading the parameters file.
    """
    try:
        with open(paramsfile, mode="r") as infile:
            params = {
                fields[0]: fields[1] 
                for line in infile 
                for fields in [line.strip().split()]
            }  # read search parameters
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e
    return params

def read_job_info(logfile: str, jobid: str) -> str:
    """Read job start time from the log file.

    This function extracts the job start time from the specified log file.

    Args:
        logfile: The path to the log file.
        jobid: The ID of the job.

    Returns:
        A string representing the job start time.

    Raises:
        IOError: If an error occurs while reading the log file.
    """
    try:
        with open(logfile, mode="r") as infile:
            log = infile.read()  # read job log data
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e
    return (next(s for s in log.split("\n") if "Job\tStart" in s)).split("\t")[-1]

def count_guides(guidesfile: str, jobid: str) -> int:
    """Count the number of guides in a guides file.

    This function reads the specified guides file and returns the number of
    guides found within it.

    Args:
        guidesfile: The path to the guides file.
        jobid: The ID of the job.

    Returns:
        The number of guides found in the file.

    Raises:
        IOError: If an error occurs while reading the guides file.
    """
    try:
        with open(guidesfile, mode="r") as infile:
            return len(infile.read().strip().split("\n"))
    except OSError as e:
        raise IOError(f"An error occurred while collecting results in {jobid}") from e

def process_genome(genome_selected: str, genome_idx: str) -> Tuple[str, str]:
    """Process genome information from job parameters.

    This function extracts the genome build and variant information from the
    given genome selection and index parameters.

    Args:
        genome_selected: The selected genome string.
        genome_idx: The genome index string.

    Returns:
        A tuple containing the genome build and a comma-separated string of
        variants.
    """
    genome = genome_selected.split("+")[0] if "+" in genome_selected else genome_selected
    variants = "Reference"
    if "+" in genome_idx:
        variants = ",".join([e.split("+")[-1] for e in genome_idx.split(",")])
    return genome, variants


# =============================================================================
# Assembly-search jobs in the History table.
#
# An assembly-search combined job's own .Params.txt is a different SHAPE from
# complete-search's (3-column index/key/value vs. 2-column key/value, and a
# different key set entirely -- Genome_paternal/Genome_maternal instead of
# Genome_selected/Genome_idx, no Genome_ref/Genome_idx at all), and its
# log.txt is a raw merged stdout/stderr dump, not complete-search's
# structured per-stage "stage\tStart"/"stage\tEnd" log read_job_info()
# expects. Feeding an assembly job through read_params()/read_job_info()/
# construct_history_summary()'s complete-search-only column logic mis-parses
# the file (read_params()'s whitespace .split() maps the leading index
# number as the key) and then KeyErrors on params["Genome_selected"] --
# confirmed against every real assembly-search job already in Results/.
#
# A small independent set of helpers here, rather than threading
# genome_type-awareness through the shared complete-search parsing, mirrors
# the same call already made for /load's own status polling (see
# refresh_assembly_search in load_page.py: "a fully independent parallel...
# NOT an extension of the existing callback").
# =============================================================================
def _is_assembly_params(paramsfile: str) -> bool:
    """Same content-substring check index.py's _is_assembly_job() uses,
    applied to a .Params.txt path retrieve_resultsdirs() already confirmed
    exists."""
    try:
        with open(paramsfile) as f:
            return "Genome_type\tassembly" in f.read()
    except OSError:
        return False


def _read_assembly_params_file(paramsfile: str) -> Dict[str, str]:
    """Parses an assembly-search combined job's .Params.txt: 3-column
    (index, key, value) -- unlike read_params()'s 2-column (key, value)
    complete-search format."""
    params: Dict[str, str] = {}
    with open(paramsfile) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                params[fields[1]] = fields[2]
    return params


def _assembly_genome_display(paternal_label: str, maternal_label: str) -> str:
    """A single Genome-column label for a two-haplotype job. Both known
    Genome_paternal/Genome_maternal label shapes (main_page.py's
    _assembly_genome_label()) end in a haplotype suffix -- strip it and
    collapse to "<individual> (assembly)" when both sides name the same
    individual; otherwise fall back to showing both labels plainly rather
    than guessing."""
    def _individual(label: str) -> Optional[str]:
        for suffix in ("_paternal", "_maternal", " (paternal)", " (maternal)"):
            if label.endswith(suffix):
                return label[: -len(suffix)]
        return None

    pat, mat = _individual(paternal_label), _individual(maternal_label)
    if pat and pat == mat:
        return f"{pat} (assembly)"
    return f"{paternal_label} / {maternal_label}"


def _assembly_max_total_edits(params: Dict[str, str], results_directory: str) -> str:
    """Assembly's own combined .Params.txt has no Max_total_edits field --
    mirrors results_page.py's own fallback for the same gap: read it back
    from either haplotype's own (complete-search-written) .Params.txt, both
    run with the same mm/bDNA/bRNA/max-total-edits by construction."""
    for dirname_key in ("Paternal_dir", "Maternal_dir"):
        hap_dir = params.get(dirname_key)
        if not hap_dir:
            continue
        hap_params_path = os.path.join(results_directory, hap_dir, PARAMS_FILE)
        if not os.path.isfile(hap_params_path):
            continue
        with open(hap_params_path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0] == "Max_total_edits":
                    return fields[1]
    return "-"


def _int_or_dash(value: Optional[str]):
    try:
        return int(value)
    except (TypeError, ValueError):
        return "-"


def _assembly_history_row(jobid: str, params: Dict[str, str], results_directory: str) -> Dict[str, object]:
    """Builds one History-table row for an assembly-search combined job:
    complete-search's columns where they apply, sane assembly-specific
    values where they don't (no single reference genome, no VCF variants,
    no structured log.txt timestamp)."""
    genome = _assembly_genome_display(
        params.get("Genome_paternal", "?"), params.get("Genome_maternal", "?")
    )
    guides_file = os.path.join(results_directory, jobid, GUIDES_FILE)
    try:
        guidesnum = count_guides(guides_file, jobid)
    except IOError:
        guidesnum = 0
    start = params.get("Job_start")
    if not start:
        # Pre-fix jobs (created before assembly_search() wrote its own
        # .Params.txt, or web-submitted before Job_start was added) have no
        # Job_start field -- the file's own mtime is a reasonable stand-in
        # so these still sort/display sensibly rather than erroring.
        params_path = os.path.join(results_directory, jobid, PARAMS_FILE)
        try:
            start = datetime.fromtimestamp(
                os.path.getmtime(params_path)
            ).strftime("%Y-%m-%d %H:%M:%S")
        except OSError:
            start = "-"
    return {
        "Job": jobid,
        "Genome": genome,
        "Variants": "-",
        "Mismatches": _int_or_dash(params.get("Mismatches")),
        "DNA bulge": _int_or_dash(params.get("DNA")),
        "RNA bulge": _int_or_dash(params.get("RNA")),
        "Max edits": _assembly_max_total_edits(params, results_directory),
        "PAM": params.get("Pam", "?"),
        "Number of Guides": guidesnum,
        "Start": start,
    }


def construct_history_summary(results: List[str]) -> pd.DataFrame:
    """Construct a summary dataframe of CRISPRme job history.

    This function reads job information from parameter and log files for each
    completed job ID and compiles a summary dataframe. The dataframe includes
    job parameters, start time, and number of guides.

    Args:
        results: A list of job IDs.

    Returns:
        A pandas DataFrame summarizing the job history.
    """
    summary = {c: [] for c in SUMMARYTABCOLS}  # initialize table
    results_directory = os.path.join(current_working_directory, RESULTS_DIR)
    for jobid in results:
        paramsfile = os.path.join(results_directory, jobid, PARAMS_FILE)
        if _is_assembly_params(paramsfile):
            row = _assembly_history_row(
                jobid, _read_assembly_params_file(paramsfile), results_directory
            )
            for col in SUMMARYTABCOLS:
                summary[col].append(row[col])
            continue
        params = read_params(os.path.join(current_working_directory, RESULTS_DIR, jobid, PARAMS_FILE), jobid)
        jobinfo = read_job_info(os.path.join(current_working_directory, RESULTS_DIR, jobid, LOG_FILE), jobid)
        guidesnum = count_guides(os.path.join(current_working_directory, RESULTS_DIR, jobid, GUIDES_FILE), jobid)
        genome, variants = process_genome(params["Genome_selected"], params["Genome_idx"])
        summary[SUMMARYTABCOLS[0]].append(jobid)  # job id
        summary[SUMMARYTABCOLS[1]].append(genome)  # genome 
        summary[SUMMARYTABCOLS[2]].append(variants)  # variants
        summary[SUMMARYTABCOLS[3]].append(int(params["Mismatches"]))  # mms
        summary[SUMMARYTABCOLS[4]].append(int(params["DNA"]))  # dna bulges
        summary[SUMMARYTABCOLS[5]].append(int(params["RNA"]))  # rna bulges
        # Max edits (total mismatches + bulges cap). Simple-mode searches are governed
        # by this slider value; advanced-mode searches use the per-type caps and record
        # their sum here. Older jobs (before this was recorded) show "-".
        _mte = params.get("Max_total_edits")
        _mode = params.get("Threshold_mode")
        if _mte in (None, "", "None"):
            _maxedits = "-"
        elif _mode == "advanced":
            _maxedits = f"{_mte} (advanced)"
        else:
            _maxedits = str(_mte)
        summary[SUMMARYTABCOLS[6]].append(_maxedits)  # max total edits
        summary[SUMMARYTABCOLS[7]].append(params["Pam"])  # pam
        summary[SUMMARYTABCOLS[8]].append(guidesnum)  # number of guides
        summary[SUMMARYTABCOLS[9]].append(jobinfo)  # job start time
        print(jobinfo)
    summary = pd.DataFrame(summary)
    # format="mixed": complete-search's own Start values are a ctime-style
    # string ("Fri 28 Aug 2026 04:09:43 PM UTC"); assembly-search rows above
    # use "%Y-%m-%d %H:%M:%S" (see _assembly_history_row). Once both job
    # types coexist in the same Results/ dir, pandas' default format
    # inference locks onto whichever shape the first row happens to be and
    # raises on every later row of the other shape -- format="mixed" parses
    # each value independently instead (pandas' own suggested fix for
    # exactly this error). utc=True: complete-search's ctime string carries
    # an explicit "UTC" (tz-aware once parsed); assembly's naive string has
    # none -- sort_values() can't compare tz-aware and tz-naive datetimes in
    # the same column, and both times genuinely are UTC (Job_start is
    # datetime.now() on this same server; the .Params.txt mtime fallback is
    # too), so coercing both to UTC-aware is correct, not just convenient.
    summary[SUMMARYTABCOLS[9]] = pd.to_datetime(
        summary[SUMMARYTABCOLS[9]], format="mixed", utc=True
    )
    summary = summary.sort_values([SUMMARYTABCOLS[9]], ascending=False)
    return summary


def get_available_results() -> pd.DataFrame:
    """Retrieve available CRISPRme job results.

    This function scans the results directory for valid job result directories
    and constructs a summary dataframe of the available results.

    Returns:
        A pandas DataFrame summarizing the available job results.
    """
    # retrieve results stored in Results folder
    results_directory = os.path.join(current_working_directory, RESULTS_DIR)
    results_dirs = retrieve_resultsdirs(results_directory) 
    assert len(results_dirs) <= len(os.listdir(results_directory))  # empty results are skipped
    return construct_history_summary(results_dirs)


def table_header_(results: pd.DataFrame) -> html.Thead:
    return html.Thead(
        html.Tr(
            [
                html.Th(c, style={"vertical-align": "middle", "text-align": "center"})
                if c not in ["Load", "Delete"]
                else html.Th("", style={"vertical-align": "middle", "text-align": "center"})
                for c in results.columns.tolist()
            ]
        )
    )

def history_table_body_(results: pd.DataFrame, page: int, max_rows: int, remaining_rows: int) -> List[html.Tr]:
    """Generate the body of the history table.

    This function creates the rows for the history table, displaying job
    information for each result. It handles pagination by displaying only a
    subset of rows based on the current page and maximum rows per page.

    Args:
        results: The DataFrame containing job history data.
        page: The current page number.
        max_rows: The maximum number of rows to display per page.
        remaining_rows: The number of rows remaining to be displayed.

    Returns:
        A list of html.Tr elements representing the table rows.
    """
    return [
        html.Tr([
            html.Td(
                html.A(
                    str(results.at[i + (page - 1) * max_rows, "Job"]),
                    target="_blank",
                    href=os.path.join(URL, f"load?job={results.at[i + (page - 1) * max_rows, 'Job']}")
                ),
                style={"vertical-align": "middle", "text-align": "center"},
            ) if col == "Job" else 
            html.Td(
                results.at[i + (page - 1) * max_rows, col],
                style={"vertical-align": "middle", "text-align": "center"},
            )
            for col in results.columns
        ])
        for i in range(min(remaining_rows, max_rows))
    ]


def display_history_table(
    results: pd.DataFrame, page: int, max_rows: Optional[int] = 1000000
) -> html.Div:
    """Display the history table with pagination.

    This function creates the HTML representation of the history table, including
    pagination controls. It generates the table header and body, and handles
    displaying only a subset of rows based on the current page and maximum rows
    per page.

    Args:
        results: The DataFrame containing job history data.
        page: The current page number.
        max_rows: The maximum number of rows to display per page.

    Returns:
        An html.Div element containing the history table.
    """
    remaining_rows = results.shape[0] - (page - 1) * max_rows
    hist_tab_body = history_table_body_(results, page, max_rows, remaining_rows)
    return [
        html.Table(
            [table_header_(results), html.Tbody(hist_tab_body)],
            style={"display": "inline-block"},
        ),
    ] + [
        html.Button(
            f"{i}",
            id=f"button-delete-history-{i}",
            **{"data-jobid": "None"},
            style={"display": "none"},
        )
        for i in range(min(remaining_rows, max_rows) - 1, 10)
    ]


def retrieve_mode() -> str:
    """Retrieve the CRISPRme running mode from the mode file.

    Returns:
        "server" or "local". Defaults to "local" if the file is absent or unreadable.
    """
    mode_file = pathlib.Path(current_working_directory) / ".mode_type.txt"
    try:
        content = mode_file.read_text(encoding="utf-8").strip()
        return content if content in ("server", "local") else "local"
    except (FileNotFoundError, OSError):
        return "local"  # safe default

def history_header() -> html.Div:
    """Create the header for the history page.

    This function generates the HTML header section for the history page,
    including the title and a brief description.

    Returns:
        An html.Div element containing the header content.
    """
    return html.Div(
        [
            html.H3("Results History"),
            html.P(
                "List of available results. Click on a link to open the " 
                "corresponding results page in a new tab."
            )
        ]
    )

def history_table(results: pd.DataFrame, mode: str) -> html.Div:
    """Create the history table or display a message if not available.

    This function generates the HTML for the history table if the CRISPRme
    running mode allows it (local mode). If running in server mode, it displays
    a message indicating that history is not available.

    Args:
        results: The DataFrame containing job history data.
        mode: The CRISPRme running mode ("server" or "local").

    Returns:
        An html.Div element containing the history table or a message.
    """
    if mode != "server":
        return html.Div(
            display_history_table(results, 1), id="div-history-table", style={"text-align": "center"}
        )
    return html.Div("History is not available while using website mode")


def history_page() -> html.Div:
    """Create the layout for the history page.

    This function generates the HTML structure for the history page, displaying
    a table of previous CRISPRme job results if available. If running in server
    mode, the history table is not displayed.

    Returns:
        An html.Div element containing the history page layout.
    """
    results = get_available_results()  # retrieve available results 
    if retrieve_mode() == "server":  # running online -> do not disply history
        results = pd.DataFrame()
    html_divs = [
        history_header(), 
        html.Div("None,None", id="div-history-filter-query", style={"display": "none"}),
        history_table(results, retrieve_mode()),
        html.Div(id="div-remove-jobid", style={"display": "none"}),
        html.Div(html.Br(), style={"text-align": "center"})
    ]
    max_page = results.shape[0]
    max_page = np.floor(max_page / 1000000) + 1
    return html.Div(html_divs, style={"margin": "1%"})

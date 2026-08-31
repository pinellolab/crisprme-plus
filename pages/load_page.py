"""Define the layout of the webpage displayed while CRISPRme is running the
analysis.

The webpage shows the status of each step of CRISPRme analysis, e.g. "done", 
"queued", etc. Moreover, the page provide the user the opportunity to check 
dinamically the state of the submitted job. 

The analysis results are kept in storage for 3 days. After 3 days the results
are automatically deleted and could not be accessed anymore.
"""

from app import app, current_working_directory, URL
from .pages_utils import RESULTS_DIR, GUIDES_FILE, LOG_FILE, PARAMS_FILE, QUEUE_FILE

from dash import Input, Output, State
from dash.exceptions import PreventUpdate
from typing import List, Optional, Tuple

from dash import dcc
from dash import html
import dash_bootstrap_components as dbc

import subprocess
import shutil
import os


# Check job completion
@app.callback(
    [
        Output("view-results", "style"),
        Output("index-status", "children"),
        Output("search-status", "children"),
        Output("post-process-status", "children"),
        Output("merge-status", "children"),
        Output("images-status", "children"),
        # NB: order matches the return tuples below (integrate_status then
        # database_status) AND the step rows in the layout. Keep these two in sync
        # -- they used to be swapped here and swapped again in the layout, which
        # only happened to cancel out.
        Output("integrate-status", "children"),
        Output("database-status", "children"),
        Output("view-results", "href"),
        Output("no-directory-error", "children"),
        Output("button-remove-result", "hidden"),
    ],
    [Input("load-page-check", "n_intervals")],
    [State("url", "search")],
)
def refresh_search(n: int, dir_name: str) -> Tuple:
    """Check the job status and refresh the webpage if the a job is nearly
    finished. The function is called every 3 seconds.

    Once completed the search, display the link to switch to results page.
    If the selected job does not exist, it is displayed a warning message.

    ...

    Parameters
    ----------
    n : int
    dir_name : str
        Job directory name

    Returns
    -------
    Tuple
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(dir_name, str):
        raise TypeError(f"Expected {str.__name__}, gopt {type(dir_name).__name__}")
    if n is None:
        raise PreventUpdate
    # recover job directories
    job_data = [
        d
        for d in os.listdir(os.path.join(current_working_directory, RESULTS_DIR))
        if os.path.isdir(os.path.join(current_working_directory, RESULTS_DIR, d))
    ]
    current_job_directory = os.path.join(
        current_working_directory, RESULTS_DIR, dir_name.split("=")[-1]
    )
    if dir_name.split("=")[-1] in job_data:
        job_data = [
            f
            for f in os.listdir(current_job_directory)
            if (
                not f.startswith(".")
                and os.path.isfile(os.path.join(current_job_directory, f))
            )
        ]
        if os.path.exists(os.path.join(current_job_directory, GUIDES_FILE[1:])):
            try:
                with open(
                    os.path.join(current_job_directory, GUIDES_FILE[1:])
                ) as handle_guides:
                    n_guides = len(handle_guides.read().strip().split("\n"))
            except OSError as e:
                raise e
        else:
            n_guides = -1
        if LOG_FILE in job_data:
            try:
                with open(os.path.join(current_job_directory, LOG_FILE)) as handle_log:
                    done = 0
                    index_status = html.P("To do", style={"color": "red"})
                    search_status = html.P("To do", style={"color": "red"})
                    post_process_status = html.P("To do", style={"color": "red"})
                    merge_status = html.P("To do", style={"color": "red"})
                    images_status = html.P("To do", style={"color": "red"})
                    database_status = html.P("To do", style={"color": "red"})
                    integrate_status = html.P("To do", style={"color": "red"})
                    current_log = handle_log.read()
                    # check if variants are required
                    variants = False
                    try:
                        with open(
                            os.path.join(current_job_directory, PARAMS_FILE)
                        ) as handle_params:
                            if "Ref_comp\tTrue" in handle_params.read():
                                variants = True
                    except OSError as e:
                        raise e
                    if variants:
                        if "Index-genome Variant\tEnd" in current_log:
                            index_status = html.P("Done", style={"color": "green"})
                            done += 1
                        elif "Index-genome Variant\tStart" in current_log:
                            index_status = html.P(
                                "Preparing enriched-genome index... [4/4]",
                                style={"color": "orange"},
                            )
                        elif "Index-genome Reference\tStart" in current_log:
                            index_status = html.P(
                                "Preparing reference-genome index... [3/4]",
                                style={"color": "orange"},
                            )
                        elif "Indexing Indels\tStart" in current_log:
                            index_status = html.P(
                                "Preparing indel-genome index... [2/4]",
                                style={"color": "orange"},
                            )
                        elif "Add-variants\tStart" in current_log:
                            index_status = html.P(
                                "Enriching genome with variants... [1/4]",
                                style={"color": "orange"},
                            )
                        elif "Search Reference\tStart" in current_log:
                            html.P("Done", style={"color": "green"})
                            done += 1
                    else:
                        if "Index-genome Reference\tEnd" in current_log:
                            index_status = html.P("Done", style={"color": "green"})
                            done += 1
                        elif "Index-genome Reference\tStart" in current_log:
                            index_status = html.P(
                                "Preparing reference-genome index... [1/1]",
                                style={"color": "orange"},
                            )
                        elif "Search Reference\tStart" in current_log:
                            index_status = html.P("Done", style={"color": "green"})
                            done += 1
                    # The pipeline logs the search stage as "Off-targets search"
                    # (see submit_job_automated_new_multiple_vcfs.sh), not "Search
                    # Reference/Variant/INDELs" — the old strings never matched, so
                    # this step never flipped to Done and the VIEW RESULTS button
                    # (gated on done == 7) never appeared.
                    if "Off-targets search\tEnd" in current_log:
                        search_status = html.P("Done", style={"color": "green"})
                        done += 1
                    elif "Off-targets search\tStart" in current_log:
                        search_status = html.P(
                            "Searching...", style={"color": "orange"}
                        )
                    if variants:
                        if (
                            "Post-analysis SNPs\tEnd" in current_log
                            and "Post-analysis INDELs\tEnd" in current_log
                        ):
                            post_process_status = html.P(
                                "Done", style={"color": "green"}
                            )
                            done += 1
                        elif "Post-analysis SNPs\tEnd" in current_log:
                            post_process_status = html.P(
                                "Post-analysis on INDELs... Step [2/2]",
                                style={"color": "orange"},
                            )
                        elif "Post-analysis SNPs\tStart" in current_log:
                            post_process_status = html.P(
                                "Post-analysis on SNPs... Step [1/2]",
                                style={"color": "orange"},
                            )
                    else:
                        if "Post-analysis\tEnd" in current_log:
                            post_process_status = html.P(
                                "Done", style={"color": "green"}
                            )
                            done += 1
                        elif "Post-analysis\tStart" in current_log:
                            post_process_status = html.P(
                                "Post-analysis... Step [1/1]", style={"color": "orange"}
                            )
                    if "Merging Targets\tEnd" in current_log:
                        merge_status = html.P("Done", style={"color": "green"})
                        done += 1
                    elif "Merging Targets\tStart" in current_log:
                        merge_status = html.P(
                            "Processing... Step [1/1]", style={"color": "orange"}
                        )
                    if "Annotating results\tStart" in current_log:
                        images_status = html.P(
                            "Annotating... Step[1/2]", style={"color": "orange"}
                        )
                    if "Creating images\tEnd" in current_log:
                        images_status = html.P("Done", style={"color": "green"})
                        done += 1
                    elif "Creating images\tStart" in current_log:
                        images_status = html.P(
                            "Generating images... Step [2/2]", style={"color": "orange"}
                        )
                    if "Integrating results\tEnd" in current_log:
                        integrate_status = html.P("Done", style={"color": "green"})
                        done += 1
                    elif "Integrating results\tStart" in current_log:
                        integrate_status = html.P(
                            "Processing... Step [1/1]", style={"color": "orange"}
                        )
                    if "Creating database\tEnd" in current_log:
                        database_status = html.P("Done", style={"color": "green"})
                        done += 1
                    elif "Creating database\tStart" in current_log:
                        database_status = html.P(
                            "Inserting data... Step [1/1]", style={"color": "orange"}
                        )
                    if (
                        os.path.isfile(
                            os.path.join(current_job_directory, "log_error.txt")
                        )
                        and os.path.getsize(
                            os.path.join(current_job_directory, "log_error.txt")
                        )
                        > 0
                    ):
                        return (
                            {"visibility": "hidden"},
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            html.P("Not available", style={"color": "red"}),
                            "",
                            dbc.Alert(
                                str(
                                    "The selected result encountered some errors, "
                                    "please remove it and try to submit again."
                                ),
                                color="danger",
                            ),
                            False,
                        )
                    # Show VIEW RESULTS as soon as the pipeline's last stage is
                    # done. Keying only on done==7 / "Job\tDone" was fragile: the
                    # search step didn't count (marker mismatch, fixed above) and
                    # the final marker is written as "Job\nDone" (newline), not a
                    # tab — so the button could stay hidden on a finished job.
                    if (
                        done == 7
                        or "Creating database\tEnd" in current_log
                        or "Job\tDone" in current_log
                        or "Job\nDone" in current_log
                    ):
                        return (
                            {"visibility": "visible"},
                            index_status,
                            search_status,
                            post_process_status,
                            merge_status,
                            images_status,
                            integrate_status,
                            database_status,
                            os.path.join(URL, f"result?job={dir_name.split('=')[-1]}"),
                            "",
                            True,
                        )
                    else:
                        return (
                            {"visibility": "hidden"},
                            index_status,
                            search_status,
                            post_process_status,
                            merge_status,
                            images_status,
                            integrate_status,
                            database_status,
                            "",  # no link to results (unfinished job)
                            "",
                            True,
                        )
            except OSError as e:
                raise e
        # job has been queued
        elif "queue.txt" in job_data:
            return (
                {"visibility": "hidden"},
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                html.P("Queued", style={"color": "red"}),
                "",
                dbc.Alert("Job submitted. Current status: IN QUEUE", color="info"),
                True,
            )
    # job data not found
    return (
        {"visibility": "hidden"},
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        html.P("Not available", style={"color": "red"}),
        "",
        dbc.Alert("The selected result does not exist", color="danger"),
        True,
    )


# =============================================================================
# Assembly-search (personal diploid genome) status polling
#
# A fully independent parallel of refresh_search/remove_result/show_full_log
# above, with its own component ids and its own dcc.Interval -- NOT an
# extension of those existing callbacks. refresh_search's parser is hardwired
# to complete-search's 7 pipeline stage strings (Index/Search/Post-process/
# Merge/Images/Integrate/Database) read via a "Ref_comp\tTrue" check in
# .Params.txt; an assembly-search job's log.txt never contains any of those
# strings (assembly_search() prints its own, completely different stage
# transitions -- see assembly_search_web_plan.md component B), so it would
# just show "To do" forever, never reaching the done gate. Building a
# separate, small parser for the 4 real stage-transition prints
# assembly_search() already emits is far lower-risk than threading
# genome_type-awareness through the existing, working complete-search
# parser.
#
# remove_result also can't be reused as-is: it deletes exactly one directory
# (the one named in ?job=), but an assembly-search job's canonical URL points
# at only <job_id>_combined -- reusing it would silently orphan the
# <job_id>_paternal/_maternal directories on "Remove result".
# =============================================================================
def _assembly_status_cell(started: bool, done: bool, progress_label: str) -> html.P:
    if done:
        return html.P("Done", style={"color": "green"})
    if started:
        return html.P(progress_label, style={"color": "orange"})
    return html.P("To do", style={"color": "red"})


@app.callback(
    [
        Output("assembly-view-results", "style"),
        Output("assembly-paternal-status", "children"),
        Output("assembly-maternal-status", "children"),
        Output("assembly-reconcile-status", "children"),
        Output("assembly-view-results", "href"),
        Output("assembly-no-directory-error", "children"),
        Output("assembly-button-remove-result", "hidden"),
    ],
    [Input("assembly-load-page-check", "n_intervals")],
    [State("url", "search")],
)
def refresh_assembly_search(n: int, dir_name: str) -> Tuple:
    """Polls an assembly-search job's Results/<job_id>_combined/log.txt,
    matching assembly_search()'s own real stage-transition prints (paternal
    search start/reuse, maternal search start/reuse, reconciliation
    start, reconciliation complete) -- no new instrumentation needed on the
    crisprme.py side. Same visual conventions (red/orange/green
    To-do/in-progress/Done) as refresh_search, for a consistent look."""
    if n is None:
        raise PreventUpdate
    not_avail = html.P("Not available", style={"color": "red"})
    job_id = dir_name.split("=")[-1]
    job_dir = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    if not os.path.isdir(job_dir):
        return (
            {"visibility": "hidden"}, not_avail, not_avail, not_avail, "",
            dbc.Alert("The selected result does not exist", color="danger"),
            True,
        )
    log_path = os.path.join(job_dir, LOG_FILE)
    if not os.path.isfile(log_path):
        if os.path.isfile(os.path.join(job_dir, QUEUE_FILE)):
            queued = html.P("Queued", style={"color": "red"})
            return (
                {"visibility": "hidden"}, queued, queued, queued, "",
                dbc.Alert("Job submitted. Current status: IN QUEUE", color="info"),
                True,
            )
        return ({"visibility": "hidden"}, not_avail, not_avail, not_avail, "", "", True)

    with open(log_path, errors="replace") as handle_log:
        current_log = handle_log.read()

    paternal_started = (
        "Running paternal haplotype search" in current_log
        or "Found existing completed paternal search results" in current_log
    )
    maternal_started = (
        "Running maternal haplotype search" in current_log
        or "Found existing completed maternal search results" in current_log
    )
    reconcile_started = "Reconciling paternal and maternal predictions" in current_log
    # Only printed AFTER combined.to_csv() succeeds (assembly_search(),
    # crisprme.py:2653-2654) -- seeing this string is itself proof the
    # combined TSV was written, no separate file-existence check needed.
    reconcile_done = "Reconciliation complete. Wrote" in current_log

    paternal_status = _assembly_status_cell(
        paternal_started, maternal_started or reconcile_started, "Searching paternal haplotype..."
    )
    maternal_status = _assembly_status_cell(
        maternal_started, reconcile_started, "Searching maternal haplotype..."
    )
    reconcile_status = _assembly_status_cell(
        reconcile_started, reconcile_done, "Reconciling against hg38..."
    )

    # Mirrors refresh_search's log_error.txt-non-empty check, adapted for a
    # single merged log.txt (stdout+stderr combined -- see
    # _run_assembly_search_job): a clean crisprme.py error() call always
    # writes "Error: ..." to stderr before exiting; an unhandled exception
    # writes a Python traceback; the web launcher's own defensive except
    # clause prefixes with "[web] launcher error:". Not exhaustive -- a raw
    # shell-level failure (e.g. a command not found) won't match any of
    # these and would just appear to hang at "To do" -- known gap, not
    # solved here.
    failed = (
        "\nError: " in current_log
        or current_log.startswith("Error: ")
        or "Traceback (most recent call last)" in current_log
        or "[web] launcher error:" in current_log
    )
    if failed and not reconcile_done:
        return (
            {"visibility": "hidden"}, not_avail, not_avail, not_avail, "",
            dbc.Alert(
                "The selected result encountered some errors, please remove "
                "it and try to submit again.",
                color="danger",
            ),
            False,
        )
    if reconcile_done:
        return (
            {"visibility": "visible"},
            paternal_status,
            maternal_status,
            reconcile_status,
            os.path.join(URL, f"result?job={job_id}"),
            "",
            False,
        )
    return (
        {"visibility": "hidden"},
        paternal_status,
        maternal_status,
        reconcile_status,
        "",
        "",
        True,
    )


@app.callback(
    Output("assembly-result-deleted", "children"),
    [Input("assembly-button-remove-result", "n_clicks")],
    [State("url", "search")],
)
def remove_assembly_result(n: int, dir_name: str) -> Optional[html.P]:
    """Deletes all three sibling directories a completed assembly-search job
    left behind (<id>_paternal, _maternal, _combined) -- not just the single
    <id>_combined directory named in the URL. Reads Paternal_dir/
    Maternal_dir from .Params.txt (written explicitly at submit time by
    submit_assembly_search_job, see assembly_search_web_plan.md) rather than
    re-deriving the sibling names by string convention here too."""
    if not n:
        raise PreventUpdate
    combined_dir_name = dir_name.split("=")[-1]
    combined_dir = os.path.join(current_working_directory, RESULTS_DIR, combined_dir_name)
    dirs_to_remove = [combined_dir]
    params_path = os.path.join(combined_dir, PARAMS_FILE)
    if os.path.isfile(params_path):
        params = {}
        with open(params_path) as pf:
            for line in pf:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    params[fields[1]] = fields[2]
        for key in ("Paternal_dir", "Maternal_dir"):
            name = params.get(key)
            if name:
                sibling_dir = os.path.join(current_working_directory, RESULTS_DIR, name)
                if os.path.isdir(sibling_dir):
                    dirs_to_remove.append(sibling_dir)
    for d in dirs_to_remove:
        shutil.rmtree(d, ignore_errors=True)
    return html.P("Results deleted")


@app.callback(
    [Output("assembly-full-log-container", "hidden"), Output("assembly-full-log-text", "value")],
    [Input("assembly-button-show-log", "n_clicks")],
    [State("url", "search"), State("assembly-full-log-container", "hidden")],
)
def show_assembly_full_log(n: int, dir_name: str, currently_hidden: bool) -> Tuple[bool, str]:
    """Same collapsible copy-pasteable log pattern as show_full_log, reading
    log.txt (assembly-search's single merged stdout+stderr log) instead of
    the separate log_error.txt/log_verbose.txt complete-search writes."""
    if not n:
        raise PreventUpdate
    if not currently_hidden:
        return True, ""
    job_dir = os.path.join(current_working_directory, RESULTS_DIR, dir_name.split("=")[-1])
    path = os.path.join(job_dir, LOG_FILE)
    if not os.path.isfile(path):
        return False, "No log file found for this job yet."
    try:
        with open(path, errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        return False, f"(could not read {LOG_FILE}: {exc})"
    tail = 400
    if len(lines) > tail:
        lines = [f"... (showing last {tail} lines) ...\n"] + lines[-tail:]
    text = "".join(lines).rstrip() or "(empty)"
    return False, text


def load_page_assembly(job_link: str = "link") -> html.Div:
    """Assembly-search's own /load page layout -- a parallel of load_page()
    below (not a branch inside it), since the pipeline stages, "remove
    result" semantics (3 directories, not 1), and log file shape (one
    merged log.txt, not log_error.txt/log_verbose.txt) all differ enough
    that sharing the layout/callbacks would mean threading conditional
    logic through code Luca already built and maintains for
    complete-search."""
    final_list = []
    final_list.append(
        html.Div(
            html.Div(
                html.Div(
                    [
                        html.P(
                            "Job submitted. Copy this link to view the status "
                            "and the result page:"
                        ),
                        html.Div(
                            html.P(
                                job_link,
                                style={"margin-top": "0.75rem", "font-size": "large"},
                            ),
                            style={
                                "border-radius": "5px",
                                "border": "2px solid",
                                "border-color": "blue",
                                "width": "100%",
                                "display": "inline-block",
                                "margin": "5px",
                            },
                        ),
                        html.P("Results will be kept available for 3 days"),
                    ],
                    style={"display": "inline-block"},
                ),
                style={
                    "display": "inline-block",
                    "background-color": "rgba(154, 208, 150, 0.39)",
                    "border-radius": "10px",
                    "border": "1px solid black",
                    "width": "70%",
                },
            ),
            style={"text-align": "center"},
        )
    )
    view_results = dcc.Link(
        html.Button(
            "View Results",
            style={
                "font-size": "large",
                "width": "700 px",
                "margin-top": "0.75rem",
                "border-radius": "5px",
                "border": "2px solid",
            },
        ),
        style={"visibility": "hidden"},
        id="assembly-view-results",
        href=URL,
    )
    final_list.append(
        html.Div(
            [
                html.H4("Status report"),
                html.Div(
                    [
                        cell
                        for label, status_id in [
                            ("Paternal haplotype search", "assembly-paternal-status"),
                            ("Maternal haplotype search", "assembly-maternal-status"),
                            ("Reconciling against hg38", "assembly-reconcile-status"),
                        ]
                        for cell in (
                            html.Div(label, className="status-step-label"),
                            html.Div(
                                "To do",
                                id=status_id,
                                className="status-step-value",
                                style={"color": "red"},
                            ),
                        )
                    ],
                    className="status-grid",
                ),
                html.Div(
                    [
                        html.Div([view_results]),
                        html.Div(id="assembly-no-directory-error"),
                        html.Div(
                            [
                                html.Button(
                                    "Remove result",
                                    id="assembly-button-remove-result",
                                    n_clicks=0,
                                    hidden=True,
                                ),
                                html.Button(
                                    "Show full log",
                                    id="assembly-button-show-log",
                                    n_clicks=0,
                                    style={"margin-left": "8px"},
                                ),
                            ]
                        ),
                        html.Div(id="assembly-result-deleted"),
                        html.Div(
                            dcc.Textarea(
                                id="assembly-full-log-text",
                                value="",
                                readOnly=True,
                                style={
                                    "width": "100%",
                                    "height": "320px",
                                    "fontFamily": "monospace",
                                    "fontSize": "12px",
                                    "whiteSpace": "pre",
                                    "overflow": "auto",
                                },
                            ),
                            id="assembly-full-log-container",
                            hidden=True,
                            style={"marginTop": "8px"},
                        ),
                    ]
                ),
            ],
            id="div-assembly-status-report",
        )
    )
    final_list.append(html.Div([view_results], style={"text-align": "center"}))
    final_list.append(html.P("", id="assembly-done"))
    final_list.append(dcc.Interval(id="assembly-load-page-check", interval=(3 * 1000)))
    return html.Div(final_list, style={"margin": "1%"})


# remove job results
@app.callback(
    Output("result-deleted", "children"),
    [Input("button-remove-result", "n_clicks")],
    [State("url", "search")],
)
def remove_result(n: int, dir_name: str) -> html.P:
    """Delete results obtained running the current CRISPRme job.

    ...

    Parameters
    ----------
    n : int
    dir_name : str
        Job directory name

    Returns
    -------
    html.P
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(dir_name, str):
        raise TypeError(f"Expected {str.__name__}, got {type(dir_name).__name__}")
    if n == 0 or n is None:
        raise PreventUpdate  # do not do anything
    elif n == 1:
        current_job_directory = os.path.join(
            current_working_directory, RESULTS_DIR, dir_name.split("=")[-1]
        )
        # remove job data
        cmd = f"rm -rf {current_job_directory}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurrend while running {cmd}")
        return html.P("Results deleted")
    return None


# Show / hide the full job log (errors + verbose tail) so users can copy-paste
# the complete error text into a bug report instead of only seeing the generic
# "encountered some errors" banner.
@app.callback(
    [Output("full-log-container", "hidden"), Output("full-log-text", "value")],
    [Input("button-show-log", "n_clicks")],
    [State("url", "search"), State("full-log-container", "hidden")],
)
def show_full_log(n: int, dir_name: str, currently_hidden: bool) -> Tuple[bool, str]:
    """Toggle the log panel; when opening, load log_error.txt + a tail of
    log_verbose.txt from the job directory into a read-only, copyable text box."""
    if not n:
        raise PreventUpdate
    if not currently_hidden:  # it's open -> collapse it
        return True, ""
    job_dir = os.path.join(
        current_working_directory, RESULTS_DIR, dir_name.split("=")[-1]
    )
    sections: List[str] = []
    # (display label, filename, tail_lines or None for full)
    for label, fname, tail in (
        ("log_error.txt", "log_error.txt", None),
        ("log_verbose.txt (tail)", "log_verbose.txt", 400),
    ):
        path = os.path.join(job_dir, fname)
        if not os.path.isfile(path):
            continue
        try:
            with open(path, "r", errors="replace") as fh:
                lines = fh.readlines()
        except OSError as exc:
            sections.append(f"===== {label} =====\n(could not read: {exc})")
            continue
        if tail and len(lines) > tail:
            lines = [f"... (showing last {tail} lines) ...\n"] + lines[-tail:]
        body = "".join(lines).rstrip() or "(empty)"
        sections.append(f"===== {label} =====\n{body}")
    text = (
        "\n\n".join(sections)
        if sections
        else "No log files found for this job yet."
    )
    return False, text


# Load Page
def load_page_no_job() -> List:
    """Empty-state shown when /load is opened without a valid job id.

    A bare visit to ``/load`` (e.g. a stale bookmark, or clicking the nav link
    directly) previously rendered a misleading "Job submitted" status box for a
    job that does not exist. Instead, tell the user plainly and point them back
    to the search page.

    Returns
    -------
    List
        Empty-state load page layout
    """

    return [
        html.Div(
            html.Div(
                [
                    html.H3("No job to display"),
                    html.P(
                        str(
                            "This page shows the status and results of a specific "
                            "CRISPRme job. Open it from the link you were given when "
                            "you submitted a search, or start a new search."
                        )
                    ),
                    dcc.Link(
                        html.Button(
                            "Start a new search",
                            style={
                                "font-size": "large",
                                "margin-top": "0.75rem",
                                "border-radius": "5px",
                                "border": "2px solid",
                            },
                        ),
                        href="/",
                    ),
                ],
                style={"display": "inline-block", "width": "70%"},
            ),
            style={"text-align": "center", "margin-top": "3rem"},
        )
    ]


def load_page(job_link: str = "link") -> List:
    """Construct the layout of the results load page. The page is displayed
    while CRISPRme analysis is running, and show the user the status of each
    analysis step.

    ...

    Parameters
    ----------
    job_link : str
        URL the user can copy to check the job status/results. Rendered inline
        so the load page no longer needs its own ``job-link`` component (which
        would duplicate the persistent ``job-link`` placeholder in the base
        layout and break the page in the browser).

    Returns
    -------
    List
        Results load page
    """

    # begin construction of results loading page
    final_list = []
    # construct box with the link to results
    final_list.append(
        html.Div(
            html.Div(
                html.Div(
                    [
                        html.P(
                            str(
                                "Job submitted. Copy this link to view the status "
                                "and the result page:"
                            )
                        ),
                        html.Div(
                            html.P(
                                job_link,
                                style={"margin-top": "0.75rem", "font-size": "large"},
                            ),
                            style={
                                "border-radius": "5px",
                                "border": "2px solid",
                                "border-color": "blue",
                                "width": "100%",
                                "display": "inline-block",
                                "margin": "5px",
                            },
                        ),
                        html.P("Results will be kept available for 3 days"),
                    ],
                    style={"display": "inline-block"},
                ),
                style={
                    "display": "inline-block",
                    "background-color": "rgba(154, 208, 150, 0.39)",
                    "border-radius": "10px",
                    "border": "1px solid black",
                    "width": "70%",
                },
            ),
            style={"text-align": "center"},
        )
    )
    # button to view results
    view_results = dcc.Link(
        html.Button(
            "View Results",
            style={
                "font-size": "large",
                "width": "700 px",
                "margin-top": "0.75rem",
                "border-radius": "5px",
                "border": "2px solid",
            },
        ),
        style={"visibility": "hidden"},
        id="view-results",
        href=URL,
    )
    # report status
    final_list.append(
        html.Div(
            [
                html.H4("Status report"),
                # One row per step, laid out as a 2-column grid (label | live
                # status) so each status always lines up with its step and the
                # statuses form a clean column. The previous layout put labels and
                # statuses in two separate <ul>s side by side, so a long label that
                # wrapped -- or the <p> margins the callback injects into the status
                # cells -- pushed the two columns out of vertical sync.
                html.Div(
                    [
                        cell
                        for label, status_id in [
                            (
                                "Preparing genome index (instant when precomputed)",
                                "index-status",
                            ),
                            ("Searching off-targets", "search-status"),
                            ("Post processing", "post-process-status"),
                            ("Merging targets", "merge-status"),
                            ("Annotating and generating images", "images-status"),
                            ("Integrating results", "integrate-status"),
                            ("Populating database", "database-status"),
                        ]
                        for cell in (
                            html.Div(label, className="status-step-label"),
                            html.Div(
                                "To do",
                                id=status_id,
                                className="status-step-value",
                                style={"color": "red"},
                            ),
                        )
                    ],
                    className="status-grid",
                ),
                html.Div(
                    [
                        html.Div([view_results]),  # hidden till analysis is completed
                        html.Div(id="no-directory-error"),
                        html.Div(
                            [
                                html.Button(
                                    "Remove result",
                                    id="button-remove-result",
                                    n_clicks=0,
                                    hidden=True,
                                ),
                                html.Button(
                                    "Show full log",
                                    id="button-show-log",
                                    n_clicks=0,
                                    style={"margin-left": "8px"},
                                ),
                            ]
                        ),
                        html.Div(id="result-deleted"),
                        # collapsible, copy-pasteable job log (errors + verbose
                        # tail) so users can grab the full error text for a report
                        html.Div(
                            dcc.Textarea(
                                id="full-log-text",
                                value="",
                                readOnly=True,
                                style={
                                    "width": "100%",
                                    "height": "320px",
                                    "fontFamily": "monospace",
                                    "fontSize": "12px",
                                    "whiteSpace": "pre",
                                    "overflow": "auto",
                                },
                            ),
                            id="full-log-container",
                            hidden=True,
                            style={"marginTop": "8px"},
                        ),
                    ]
                ),
            ],
            id="div-status-report",
        )
    )
    # view results button
    final_list.append(
        html.Div([view_results], style={"text-align": "center"}),
    )
    final_list.append(html.P("", id="done"))
    final_list.append(dcc.Interval(id="load-page-check", interval=(3 * 1000)))
    load_page = html.Div(final_list, style={"margin": "1%"})
    return load_page

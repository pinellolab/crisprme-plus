"""Script to manage and print the results page. The results page allows the user 
to navigate through the CRISPRme analysis results, providing several different
filtering options to visualize the results.

The result page consists of 2 (?) layers the ... on the top of the page and 
tab selection layer to select the way to visualize the results. 
Options (tabs) to visualize the results:
    - Custom Ranking
        Display the best targets for each input guide, according to the scoring
        criterion selected by the user.
    - Summary by Mismatch/Bulges
    - Summary by Sample
        Display the best targets for each individual sample considered during 
        CRISPRme analysis.
    - Query Genomic Region
        Display the best targets for each input guide in a specific genomic
        region. The results are sorted according to the scoring criterion 
        selected by the user.
    - Graphical Reports
        Plot information regarding the target guides (plots computed at execution
        time -> could require few secs to complete).
    - Personal Risk Card (displayed only if used individual data)

The results could be sorted and filtered according to 3 criteria:
    - CFD score
    - CRISTA score
    - Number of mismacthes and bulges

TODO: complete doc string with missing info --> read paper carefully
"""

from .pages_utils import (
    GUIDES_FILE,
    PAGE_SIZE,
    BARPLOT_LEN,
    COL_BOTH,
    COL_BOTH_TYPE,
    COL_BOTH_RENAME,
    GUIDE_COLUMN,
    CHR_COLUMN,
    VARIANTS_CRISTA,
    VARIANTS_CFD,
    VARIANTS_FEWEST,
    RESULTS_DIR,
    DATA_DIR,
    IMGS_DIR,
    FILTERING_CRITERIA,
    PARAMS_FILE,
    SAMPLES_ID_FILE,
    CAS9,
    PANDAS_OPERATORS,
    drop_columns,
    write_json,
    read_json,
    get_query_column,
    split_filter_part,
    generate_table,
    generate_table_samples,
)
from app import (
    app,
    cache,
    app_directory,
    current_working_directory,
    URL,
)
from PostProcess.supportFunctions.loadSample import associateSample
from PostProcess import CFDGraph, query_manager
from PostProcess.assembly_reconcile import (
    PRED_COLS,
    find_results_prefix,
    load_crisprme_predictions,
    load_unlifted_ids,
)
from PostProcess.generate_report import read_specificity_score

from dash.exceptions import PreventUpdate
from dash import Input, Output, State
from typing import Dict, List, Optional, Tuple
from glob import glob

from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import pandas as pd

import subprocess
import math
import base64  # for decoding upload content
import plotly.graph_objects as go
from dash import dash_table
import sqlite3
import flask
import errno
import re
import os


# -------------------------------------------------------------------------------
# Result page layout
#


def _perfect_match_sites(integrated_tsv):
    """Distinct perfect-match (0 mismatch + 0 bulge) sites in an integrated TSV.

    A guide with >1 perfect genomic match has no a-priori on-target -- each is an
    equally-efficient candidate cut site (a perfect-match off-target is the
    highest-risk class). Reads ONLY the locus + mm+b columns so it is fast even on
    a genome-wide result. Fully guarded: returns [] on any error (never blocks the
    results page). Each entry is a (chrom, pos, strand) tuple.
    """
    try:
        import pandas as pd

        header = pd.read_csv(integrated_tsv, sep="\t", nrows=0).columns.tolist()

        def _col(*cands):
            for c in cands:
                if c in header:
                    return c
            return None

        mmb = _col("Mismatches+bulges_(highest_CFD)", "Mismatches+bulges")
        chrom = _col("Chromosome")
        pos = _col("Start_coordinate_(highest_CFD)", "Start_coordinate")
        strand = _col("Strand_(highest_CFD)", "Strand")
        if not mmb:
            return []
        usecols = [c for c in (mmb, chrom, pos, strand) if c]
        df = pd.read_csv(integrated_tsv, sep="\t", usecols=usecols, dtype=str)
        mask = pd.to_numeric(df[mmb], errors="coerce") == 0
        seen, out = set(), []
        for _, r in df[mask].iterrows():
            key = (
                str(r.get(chrom, "?")),
                str(r.get(pos, "?")),
                str(r.get(strand, "")),
            )
            if key not in seen:
                seen.add(key)
                out.append(key)
        return out
    except Exception:
        return []


def result_page(job_id: str) -> html.Div:
    """Print the results page layout (guides table + images).
    The guides table contains the research profile found during
    target search. Creates 10 buttons (mismatch number + 2), the
    remaining ones are set to style = {"display":None}, in order
    to have the right number of buttons, based on mismatches required
    in input during the target search. This choice solves some
    callback issues that have in input elements not created. In this
    case, all the possible buttons are created, but are shown only
    those correct based on the selected number of mismatches.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier

    Returns
    -------
    html.Div
        Result page layout
    """

    # check input function arguments
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    # start result page creation code
    value = job_id
    job_directory = os.path.join(current_working_directory, "Results", f"{job_id}")
    # check existance and zip integrated file
    integrated_file_name = glob(
        os.path.join(current_working_directory, "Results", f"{job_id}", "*integrated*")
    )[
        0
    ]  # take the first list element
    assert isinstance(integrated_file_name, str)
    integrated_file_name_zip = integrated_file_name.replace("tsv", "zip")
    # check existence and zip alt_merge file
    alt_merge_file_name = glob(
        os.path.join(
            current_working_directory,
            "Results",
            f"{job_id}",
            "*all_results_with_alternative_alignments*",
        )
    )[
        0
    ]  # take the first list element
    assert isinstance(alt_merge_file_name, str)
    alt_merge_file_name_zip = alt_merge_file_name.replace("tsv", "zip")
    # check job directory existence to avoid crush
    if not os.path.isdir(job_directory):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    count_guides = 0
    guides_file = os.path.join(current_working_directory, "Results", f"{value}", GUIDES_FILE)
    assert os.path.isfile(guides_file)
    try:
        with open(guides_file) as handle:
            for line in handle:
                count_guides += 1
    except:
        raise IOError(f"Unable to read {guides_file}.")
    finally:
        handle.close()
    # Load mismatches
    try:
        with open(os.path.join(current_working_directory, RESULTS_DIR, value, PARAMS_FILE)) as p:
            all_params = p.read()
            real_genome_name = (
                next(s for s in all_params.split("\n") if "Genome_idx" in s)
            ).split("\t")[-1]
            mms = (next(s for s in all_params.split("\n") if "Mismatches" in s)).split(
                "\t"
            )[-1]
            bulge_dna = (next(s for s in all_params.split("\n") if "DNA" in s)).split(
                "\t"
            )[-1]
            bulge_rna = (next(s for s in all_params.split("\n") if "RNA" in s)).split(
                "\t"
            )[-1]
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
            max_bulges = (
                next(s for s in all_params.split("\n") if "Max_bulges" in s)
            ).split("\t")[-1]
            # DNA and RNA bulges are INDEPENDENT per-type caps and a single alignment
            # can carry BOTH, so the total bulge span (and the mm+bulge column extent)
            # is their SUM -- not Max_bulges (the older single/max value). Fall back to
            # Max_bulges for jobs created before the DNA/RNA fields existed.
            def _param_val(key):
                return next(
                    (s.split("\t")[-1] for s in all_params.split("\n")
                     if s.startswith(key + "\t")),
                    None,
                )
            try:
                total_bulges = int(_param_val("DNA")) + int(_param_val("RNA"))
            except (TypeError, ValueError):
                total_bulges = int(max_bulges)
            pam_name = (next(s for s in all_params.split("\n") if "Pam" in s)).split(
                "\t"
            )[-1]
            # Optional: the total mismatches+bulges cap that actually bound the search
            # (the "Max edits" the user set). next(..., None) keeps result pages for jobs
            # created before this field existed loadable.
            max_total_edits = next(
                (
                    s.split("\t")[-1]
                    for s in all_params.split("\n")
                    if s.startswith("Max_total_edits")
                ),
                None,
            )
            # Which threshold control governed the search: 'simple' (the single
            # "Max edits" slider) or 'advanced' (explicit per-type mm/DNA/RNA caps).
            # Absent for jobs created before this field -> inferred below.
            threshold_mode = next(
                (
                    s.split("\t")[-1]
                    for s in all_params.split("\n")
                    if s.startswith("Threshold_mode")
                ),
                None,
            )
    except OSError as e:
        raise e
    finally:
        p.close()
    # recover genome name
    genome_name = genome_type_f
    if "+" in real_genome_name:
        genome_name = [genome_name] + [
            name.split("+")[1] for name in real_genome_name.strip().split(",")
        ]
        genome_name = "+".join(genome_name)
    if "True" in ref_comp:
        genome_type = "both"
    else:
        genome_type = "ref"
    mms = int(mms[0])
    # load acfd for each guide
    acfd_file = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        "".join([".", job_id, ".acfd_CFD.txt"]),
    )
    if not os.path.isfile(acfd_file):
        # something went wrong
        raise FileNotFoundError(f"Unable to locate {acfd_file}")
    try:
        with open(acfd_file) as handle:
            all_scores = handle.read().strip().split("\n")
    except OSError as e:
        raise e
    finally:
        handle.close()
    guides_error_file = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, "guides_error.txt"
    )
    list_error_guides = []
    if os.path.exists(guides_error_file):
        try:
            with open(guides_error_file) as handle_error_g:
                for e_g in handle_error_g:
                    list_error_guides.append(e_g.strip())
        except OSError as e:
            raise e
        finally:
            handle_error_g.close()
    col_targetfor = "("
    for i in range(1, (mms + total_bulges)):
        col_targetfor = "".join([col_targetfor, str(i), "-"])
    col_targetfor = "".join([col_targetfor, str(mms + total_bulges)])
    col_targetfor = " ".join([col_targetfor, "Mismatches + Bulges)"])
    # Column of headers. Remove the entries accordingly when checking genome type
    columns_profile_table = [
        {"name": ["", "gRNA (spacer+PAM)"], "id": "Guide", "type": "text"},
        {"name": ["", "Nuclease", ""], "id": "Nuclease", "type": "text"},
        {
            "name": ["", "Aggregated Specificity Score (0-100)"],
            "id": "CFD",
            "type": "text",
        },
        {
            "name": ["Off-targets for Mismatch (MM) and Bulge (B) Value", "Total"],
            "id": "Total",
            "type": "text",
        },
    ]
    columns_profile_table.append(
        {
            "name": ["Off-targets for Mismatch (MM) and Bulge (B) Value", "# Bulges"],
            "id": "# Bulges",
            "type": "text",
        }
    )
    for i in range(mms + 1):
        columns_profile_table.append(
            {
                "name": [
                    "Off-targets for Mismatch (MM) and Bulge (B) Value",
                    "".join([str(i), "MM"]),
                ],
                "id": "".join([str(i), "MM"]),
                "type": "text",
            }
        )
    remove_indices = set()
    if "NO SCORES" in all_scores:
        # remove CFD and Doench header from table
        remove_indices.add("CFD", "Doench 2016", "Reference", "Enriched")
    if genome_type == "ref":
        # remove reference header
        remove_indices.update(["Reference", "Enriched"])
    else:
        # remove reference and reference target headers
        remove_indices.update(
            [
                "Reference",
                "Enriched",
                "On-Targets Reference",
                "Samples in Class 0 - 0+ - 1 - 1+",
            ]
        )
    # Remove headers not used in the selected search results
    columns_profile_table = [
        i
        for j, i in enumerate(columns_profile_table)
        if columns_profile_table[j]["id"] not in remove_indices
    ]
    # build final result page with corresponding fields
    final_list = []
    if list_error_guides:
        final_list.append(
            dbc.Alert(
                [
                    "Warning: Some guides have too many targets! ",
                    html.A(
                        "Click here",
                        href=os.path.join(URL, DATA_DIR, job_id, "guides_error.txt"),
                        className="alert-link",
                    ),
                    " to view them",
                ],
                color="warning",
            )
        )
    # Multiple perfect matches (0 mm + 0 bulge) => no a-priori on-target. A
    # perfect-match OFF-target cuts as efficiently as the intended site, so warn
    # loudly and point to the report (where all are forced into the panel + flagged).
    _perfect_sites = _perfect_match_sites(integrated_file_name)
    if len(_perfect_sites) >= 2:
        _site_txt = ", ".join(f"{c}:{p}" for c, p, _s in _perfect_sites[:10])
        if len(_perfect_sites) > 10:
            _site_txt += ", …"
        final_list.append(
            dbc.Alert(
                [
                    html.Strong(
                        f"⚠ Multiple perfect matches — no unambiguous "
                        f"on-target. "
                    ),
                    f"This guide matches {len(_perfect_sites)} genomic sites with "
                    f"0 mismatches and 0 bulges ({_site_txt}). Each is a candidate "
                    f"cut site — a perfect-match off-target cuts as efficiently "
                    f"as the intended one. All are placed at the top of the "
                    f"validation panel and flagged ",
                    html.Code("Perfect_match = Yes"),
                    " in the downloadable report.",
                ],
                color="danger",
            )
        )
    # Present the title according to which threshold control governed the search.
    # SIMPLE mode: show ONLY the single "Max edits" (mismatches + bulges) cap the user
    # set. The per-type mm/DNA/RNA caps are wide internal defaults (or the ones Load
    # Example fills in), and showing them next to a tighter total cap is contradictory
    # (e.g. "Mismatches 4 - DNA bulges 1 - RNA bulges 1 - Max edits 1"). ADVANCED mode:
    # show the explicit per-type caps the user chose (the total cap equals their sum).
    try:
        _mte_int = int(max_total_edits) if max_total_edits is not None else None
    except (TypeError, ValueError):
        _mte_int = None
    try:
        _per_type_sum = int(mms) + int(bulge_dna) + int(bulge_rna)
    except (TypeError, ValueError):
        _per_type_sum = None
    if threshold_mode == "simple":
        _is_simple = True
    elif threshold_mode == "advanced":
        _is_simple = False
    else:
        # Pre-existing job without the field: a total cap tighter than the per-type
        # sum means the "Max edits" slider governed (simple mode). If the cap wasn't
        # recorded at all (pre-alpha.15 job), fall back to the per-type title.
        _is_simple = (
            _mte_int is not None and _per_type_sum is not None and _mte_int < _per_type_sum
        )
    _summary_parts = ["Result Summary", "-", genome_name, "-", pam_name]
    if _is_simple and _mte_int is not None:
        _summary_parts += ["-", "Max edits (mismatches + bulges)", str(_mte_int)]
    else:
        _summary_parts += [
            "-", "Mismatches", str(mms),
            "-", "DNA bulges", bulge_dna,
            "-", "RNA bulges", bulge_rna,
        ]
    final_list.append(html.H3(" ".join(_summary_parts)))
    # The detailed explanation of this matrix is rendered as a caption BELOW the table
    # (see `matrix_explanation` after the DataTable), so it reads right next to the cells
    # it describes and stays consistent, word-for-word in message, with the downloadable
    # report's matrix caption.
    # define upper page box
    final_list.append(
        # Single PROMINENT, CENTERED "Download Report" action (replaces the former
        # three raw-file links). The self-contained report ZIP bundles the integrated
        # results, all tiers, the complete raw table, and the high-variant-density
        # regions BED, so it supersedes them. A polling interval reveals the button
        # once the report has been generated at the end of the job. Boxed + centered
        # so a user cannot miss it.
        html.Div(
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            html.H4(
                                "📄 Off-target report",
                                style={
                                    "margin": "0 0 8px 0",
                                    "fontWeight": "700",
                                    "color": "#1a365d",
                                },
                            ),
                            html.P(
                                "Preparing the off-target report for download…",
                                id="download-link-report",
                                style={
                                    "margin": "0",
                                    "fontSize": "1.05rem",
                                    "textAlign": "center",
                                },
                            ),
                            dcc.Interval(interval=2 * 1000, id="interval-report"),
                            html.Div(
                                os.path.join(
                                    current_working_directory,
                                    RESULTS_DIR,
                                    job_id,
                                    job_id + "_report.zip",
                                ),
                                style={"display": "none"},
                                id="div-info-report",
                            ),
                        ],
                        style={"textAlign": "center"},
                    ),
                    width={"size": 10},
                ),
                justify="center",
            ),
            style={
                "textAlign": "center",
                "margin": "16px auto 24px auto",
                "padding": "22px",
                "background": "#eef6ff",
                "border": "2px solid #2b6cb0",
                "borderRadius": "12px",
                "boxShadow": "0 2px 10px rgba(43,108,176,0.18)",
            },
        )
    )
    # results table (middle of page layout)
    final_list.append(
        html.Div(
            html.Div(
                dash_table.DataTable(
                    id="general-profile-table",
                    # page_size=PAGE_SIZE,
                    columns=columns_profile_table,
                    merge_duplicate_headers=True,
                    # fixed_rows={ 'headers': True, 'data': 0 },
                    # data = profile.to_dict('records'),
                    selected_cells=[{"row": 0, "column": 0}],
                    # layout CSS style
                    css=[
                        {
                            "selector": ".row",
                            "rule": "margin: 0",
                            "selector": "td.cell--selected, td.focused",
                            "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                        },
                        {
                            "selector": "td.cell--selected *, td.focused *",
                            "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                        },
                    ],
                    page_current=0,
                    page_size=10,
                    # The Result Summary has one row per guide (usually one), each a tall
                    # multi-line REF/VAR breakdown. Show every row at once -- no pager,
                    # no vertical scroll -- so the whole summary fits on the page instead
                    # of a weird pagination bar + scrollbar on a one-row table. Only wide
                    # tables scroll horizontally (overflowX auto).
                    page_action="none",
                    # virtualization = True,
                    filter_action="custom",
                    filter_query="",
                    sort_action="custom",
                    sort_mode="multi",
                    sort_by=[],
                    style_table={
                        # 'margin-left': "10%",
                        "overflowX": "auto",
                    },
                    style_data={
                        "whiteSpace": "pre",
                        "height": "auto",
                        "font-size": "1.30rem",
                    },
                    # style_cell={
                    #    'width':f'{1/len(columns_profile_table)*100}%'
                    # },
                    style_data_conditional=[
                        {
                            "if": {"column_id": "Genome"},
                            "font-weight": "bold",
                            "textAlign": "center",
                        },
                        # {'if': {'column_id': 'Guide'},
                        #                    'width': '10%',
                        #                    }
                    ],
                    style_cell_conditional=[
                        {
                            "if": {"column_id": "Guide"},
                            "width": "20%",
                        },
                        {
                            "if": {"column_id": "Total"},
                            "width": "15%",
                        },
                        {
                            "if": {"column_id": "Doench 2016"},
                            "width": "5%",
                        },
                        {
                            "if": {"column_id": "# Bulges"},
                            "width": "5%",
                        },
                        {
                            "if": {"column_id": "Nuclease"},
                            "width": "5%",
                        },
                    ],
                    #                        {'if': {'column_id': 'Reference'},
                    #                        'width': '10%',
                    #                        }],
                ),
                id="div-general-profile-table",
                style={"margin-left": "5%", "margin-right": "5%"},
            )
        )
    )
    # Explanation BELOW the table, consistent in message with the downloadable report's
    # matrix caption: what a cell is, the origin split, the perfect-match cell, and that
    # every site is counted once (one row carrying both alleles as columns).
    _cap_style = {
        "margin": "0 5% 0 5%",
        "fontSize": "1.02rem",
        "color": "#334155",
        "lineHeight": "1.5",
    }
    _budget_note = []
    if _mte_int is not None:
        _budget_note = [
            " Some entries combine mismatches and bulges beyond the ",
            html.Strong(f"Max edits total ({_mte_int})"),
            ": the search caps mismatches and bulges independently, so an alignment's "
            "total can exceed that budget while the same site stays within scope on its "
            "other (reference or variant) allele — the same distinct sites, not extra "
            "off-target risk.",
        ]
    final_list.append(
        html.Div(
            [
                html.P(
                    [
                        "Each entry counts ",
                        html.Strong("distinct putative off-target sites"),
                        ", grouped by mismatch count and bulge size and split by origin: ",
                        html.Strong("REFERENCE"),
                        " — the target exists in the reference genome (even if a variant "
                        "also alters it) — vs ",
                        html.Strong("VARIANT"),
                        " — the target exists only because a variant creates it (",
                        html.Code("Not_found_in_REF"),
                        ").",
                    ],
                    style=_cap_style,
                ),
                html.P(
                    [
                        "The ",
                        html.Strong("0 MM / 0 B"),
                        " entry holds the guide's ",
                        html.Strong("perfect genomic match(es)"),
                        " — the candidate on-target(s). The intended on-target cannot be "
                        "told from a perfect-match off-target by sequence alone, so every "
                        "perfect match is reported and forced to the top of the validation "
                        "panel in the downloadable report.",
                    ],
                    style={**_cap_style, "marginTop": "0.5em"},
                ),
                html.P(
                    [
                        "Every site is counted ",
                        html.Strong("once"),
                        ": a site is a single row that carries both its reference-genome "
                        "alignment and its variant-carrier alignment as side-by-side "
                        "columns, so it is never double-counted across cells.",
                        *_budget_note,
                    ],
                    style={**_cap_style, "marginTop": "0.5em"},
                ),
            ],
            style={"marginTop": "10px", "marginBottom": "6px"},
        )
    )
    final_list.append(html.Br())  # add space between HTML lines
    # drop-down bar (filetring criterion selection)
    final_list.append(
        html.Div(
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            html.H4("Select filter criteria for targets"),
                            dcc.Dropdown(
                                options=[
                                    {"label": "CFD score", "value": "CFD"},
                                    {"label": "CRISTA Score", "value": "CRISTA"},
                                    {
                                        "label": "Fewest Mismatches and Bulges",
                                        "value": "fewest",
                                    },
                                ],
                                value="CFD",
                                id="target_filter_dropdown",
                            ),
                            dcc.Store(id="store"),
                        ]
                    )
                )
            ),
        )
    )
    final_list.append(html.Br())  # add space between HTML lines
    if genome_type == "ref":
        final_list.append(
            dcc.Tabs(
                id="tabs-reports",
                value="tab-query-table",
                children=[
                    dcc.Tab(label="Custom Ranking", value="tab-query-table"),
                    dcc.Tab(
                        label="Summary by Mismatches/Bulges",
                        value="tab-summary-by-guide",
                    ),
                    dcc.Tab(
                        label="Query Genomic Region", value="tab-summary-by-position"
                    ),
                    dcc.Tab(label="Graphical Reports", value="tab-summary-graphical"),
                ],
            )
        )
    else:
        # Barplot for population distributions
        final_list.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Button(
                                    "Show/Hide Target Distribution in SuperPopulations",
                                    id="btn-collapse-populations",
                                )
                            ),
                        ]
                    ),
                    dbc.Collapse(
                        dbc.Card(
                            dbc.CardBody(html.Div(id="content-collapse-population"))
                        ),
                        id="collapse-populations",
                    ),
                ],
                hidden=True,
            )
        )
        final_list.append(html.Br())  # add space between HTML lines
        # define results page tabs
        final_list.append(
            dcc.Tabs(
                id="tabs-reports",
                value="tab-query-table",
                children=[
                    dcc.Tab(label="Custom Ranking", value="tab-query-table"),
                    dcc.Tab(
                        label="Summary by Mismatches/Bulges",
                        value="tab-summary-by-guide",
                    ),
                    dcc.Tab(label="Summary by Sample", value="tab-summary-by-sample"),
                    dcc.Tab(
                        label="Query Genomic Region", value="tab-summary-by-position"
                    ),
                    dcc.Tab(label="Graphical Reports", value="tab-summary-graphical"),
                    dcc.Tab(
                        label="Personal Risk Cards", value="tab-graphical-sample-card"
                    ),
                ],
            )
        )
    final_list.append(html.Div(id="div-tab-content"))

    final_list.append(
        html.Div(genome_type, style={"display": "none"}, id="div-genome-type")
    )
    result_page = html.Div(final_list, style={"margin": "1%"})
    return result_page


def _encode_png(path: str) -> Optional[str]:
    """Base64 data-URI for a PNG file, or None if it can't be read -- same
    inline-embedding pattern already used throughout this file for
    complete-search's own Graphical Reports plots (e.g. the radar-chart and
    top-N images), reused as-is here rather than reinvented."""
    try:
        with open(path, "rb") as fh:
            return "data:image/png;base64," + base64.b64encode(fh.read()).decode()
    except OSError:
        return None


def result_page_assembly(job_id: str) -> html.Div:
    """Results page for an assembly-search (personal diploid genome) job --
    fully independent from result_page() above, not a branch inside it.
    result_page() crashes on this job type: it does unguarded .Params.txt
    lookups (Genome_idx, Ref_comp, Genome_selected -- no fallback, so
    StopIteration) and assumes complete-search's exact output filenames
    (glob for *integrated*, a per-job SQLite .db, .acfd_CFD.txt) from its
    very first lines, well before its own genome_type tab-set branch --
    none of that exists for an assembly-search job. See
    assembly_search_web_plan.md component D for the full investigation
    (confirmed by directly reading result_page() and generate_sample_card,
    not assumed).

    First-pass scope, deliberately: the reconciled off-target table (the 3
    real `origin` categories reconcile_haplotypes() can actually produce --
    paternal_only, maternal_only, both; a 4th, "both_haplotype_private", is
    defined in assembly_reconcile.py but not called anywhere in the current
    pipeline -- confirmed by grepping for its call sites -- so it can't
    appear in real data yet) plus a haplotype-coverage summary, including
    the two non-mappable ("haplotype-private") site COUNTS.
    reconcile_haplotypes() doesn't persist per-site detail for non-mappable
    predictions anywhere, only a count, so a detailed haplotype-private
    table isn't buildable from current pipeline output -- not a UI choice,
    a real upstream data gap. A per-haplotype (maternal vs. paternal)
    comparison view, in the spirit of complete-search's "Personal Risk
    Cards" (a UX-shape precedent only -- none of its actual code, built on
    a per-job SQLite database and VCF-sample files, applies here), is
    deferred pending design discussion.
    """
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    if not os.path.isdir(job_directory):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))

    params = {}
    params_path = os.path.join(job_directory, PARAMS_FILE)
    if os.path.isfile(params_path):
        with open(params_path) as pf:
            for line in pf:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    params[fields[1]] = fields[2]

    combined_tsv_matches = glob(os.path.join(job_directory, "*_combined_hg38.tsv"))
    if not combined_tsv_matches:
        return html.Div(
            html.Div(
                [
                    html.H3(f"Personal Assembly Search Results — {job_id}"),
                    dbc.Alert(
                        "This job's reconciled results file was not found -- "
                        "the run may not have finished yet.",
                        color="warning",
                    ),
                ],
                style={"margin": "1%"},
            )
        )
    df = pd.read_csv(combined_tsv_matches[0], sep="\t")

    # Display-only column trimming -- doesn't touch combined_hg38.tsv on
    # disk or the assembly_reconcile.py PRED_COLS that produced it, just
    # what this table renders:
    #  - Spacer+PAM_paternal/_maternal are always identical (checked
    #    against real diverse test data: 0/542 "both" rows differ -- it's
    #    the same guide regardless of haplotype), so merge into one
    #    "Spacer+PAM" column, falling back to whichever haplotype actually
    #    has a value on paternal_only/maternal_only rows.
    #  - The *_ALT_(fewest_mm+b)_* columns are always "NA": there's no VCF
    #    in an assembly-search haplotype search, so there's no alternate
    #    allele for them to hold (REF/ALT_origin_(fewest_mm+b) is "ref" for
    #    every row, checked directly).
    #  - off_target_id_paternal/_maternal are internal join/dedup keys, not
    #    something a user needs to read.
    if "Spacer+PAM_paternal" in df.columns and "Spacer+PAM_maternal" in df.columns:
        df["Spacer+PAM"] = df["Spacer+PAM_paternal"].combine_first(
            df["Spacer+PAM_maternal"]
        )
    _HIDDEN_DISPLAY_COLS = {
        "Spacer+PAM_paternal",
        "Spacer+PAM_maternal",
        "off_target_id_paternal",
        "off_target_id_maternal",
        "Aligned_protospacer+PAM_ALT_(fewest_mm+b)_paternal",
        "Aligned_protospacer+PAM_ALT_(fewest_mm+b)_maternal",
    }
    display_cols = ["Spacer+PAM"] if "Spacer+PAM" in df.columns else []
    display_cols += [
        c for c in df.columns if c not in _HIDDEN_DISPLAY_COLS and c != "Spacer+PAM"
    ]

    # Haplotype coverage summary, parsed from reconcile_haplotypes()'s own
    # structured log (log_verbose.txt, NOT log.txt -- the latter only has
    # assembly_search()'s coarser stage-transition prints, see component B)
    # -- the only place the non-mappable counts are available at all; there
    # is no per-site detail for them anywhere in current pipeline output.
    summary_counts = {}
    log_verbose_path = os.path.join(job_directory, "log_verbose.txt")
    if os.path.isfile(log_verbose_path):
        with open(log_verbose_path, errors="replace") as lf:
            lines = lf.read().splitlines()
        if "Reconciliation complete:" in lines:
            start = lines.index("Reconciliation complete:") + 1
            for line in lines[start:]:
                line = line.strip()
                if not line or ":" not in line:
                    break
                key, _, val = line.partition(":")
                try:
                    summary_counts[key.strip()] = int(val.strip())
                except ValueError:
                    break

    origin_counts = df["origin"].value_counts().to_dict() if "origin" in df.columns else {}

    # Per-haplotype non-mappable ("private") sites: dropped from `df`'s rows
    # entirely (reconcile_haplotypes() only carries successfully-lifted
    # predictions into the combined table), but their per-site detail is
    # NOT actually lost -- it's sitting untouched in each haplotype's own
    # complete-search output, and reconcile_haplotypes() already writes the
    # exact non-mappable off_target_id list to
    # {paternal,maternal}_offtargets_not_lifted.bed. Cross-referencing the
    # two recovers real per-site detail (coordinates/sequence/CFD in the
    # haplotype's OWN coordinate system, not hg38 -- these sites have none)
    # with no core/CLI change (corrected an earlier wrong finding here --
    # see assembly_search_web_plan.md component D). Only splits into
    # "private to paternal" / "private to maternal" -- NOT further resolved
    # against whether the same site is also present (under a different,
    # also-unmappable representation) in the other haplotype; that
    # resolution needs impg and is deferred (see design doc).
    # Loaded once per haplotype and reused for two different purposes below:
    # (1) private_tables (non-mappable per-site detail, as before) and (2)
    # the new per-haplotype guide summary table + CFD distribution plot.
    # Previously this only loaded when a haplotype actually had non-mappable
    # sites (an early `continue`) -- widened to load unconditionally since
    # the summary table needs each haplotype's full prediction set
    # regardless of non-mappable status.
    private_tables: Dict[str, pd.DataFrame] = {}
    hap_predictions_all: Dict[str, pd.DataFrame] = {}
    hap_results_dirs: Dict[str, str] = {}
    for hap, dirname_key in (("paternal", "Paternal_dir"), ("maternal", "Maternal_dir")):
        hap_dir_name = params.get(dirname_key)
        if not hap_dir_name:
            continue
        hap_results_dir = os.path.join(current_working_directory, RESULTS_DIR, hap_dir_name)
        hap_results_dirs[hap] = hap_results_dir
        try:
            prefix = find_results_prefix(hap_results_dir)
            # merge_bp hardcoded to 3: matches the --merge value
            # submit_assembly_search_job's cmd always passes today (not yet
            # user-configurable). Must match what the job actually ran with
            # -- cluster_collapse's row order (and therefore the
            # off_target_id -> row mapping) depends on it, same as
            # reconcile_haplotypes()'s own call.
            # Widened past the default PRED_COLS (web-side only -- doesn't
            # touch assembly_reconcile.py or what the CLI itself writes) to
            # also pull Bulge_type_(fewest_mm+b), needed for the mismatch +
            # bulge breakdown table below to show the same "Bulge Type"
            # column complete-search's own summary table has.
            hap_predictions = load_crisprme_predictions(
                hap_results_dir, prefix, merge_bp=3,
                cols=PRED_COLS + ["Bulge_type_(fewest_mm+b)"],
            )
        except (FileNotFoundError, OSError):
            continue
        hap_predictions["off_target_id"] = hap_predictions["off_target_id"].astype(str)
        hap_predictions_all[hap] = hap_predictions
        not_lifted_path = os.path.join(job_directory, f"{hap}_offtargets_not_lifted.bed")
        if not os.path.isfile(not_lifted_path):
            continue
        unlifted_ids = load_unlifted_ids(not_lifted_path)
        if not unlifted_ids:
            continue
        private = hap_predictions[hap_predictions["off_target_id"].isin(unlifted_ids)]
        if not private.empty:
            private_tables[hap] = private.drop(columns=["off_target_id"])

    # Per-haplotype guide summary, in the spirit of complete-search's own
    # Result Summary table (one row per guide: Guide/Nuclease/aggregate
    # specificity score/Total) -- adapted to show paternal vs. maternal
    # side by side instead of a single reference-genome value.
    # read_specificity_score() is complete-search's own aggregate-score
    # reader (PostProcess/generate_report.py, unmodified, read-only reuse)
    # -- same "100/(100+sum_cfds)" definition already shown on every
    # complete-search results page, not a new formula invented here.
    guide_nuclease = "?"
    for hap in ("paternal", "maternal"):
        hap_params_path = os.path.join(hap_results_dirs.get(hap, ""), PARAMS_FILE)
        if os.path.isfile(hap_params_path):
            with open(hap_params_path) as hp:
                for line in hp:
                    # each haplotype's .Params.txt is complete-search's own
                    # native format: "key\tvalue" (2 columns, no leading
                    # index) -- NOT the combined job's own .Params.txt
                    # format parsed into `params` above ("index\tkey\tvalue",
                    # 3 columns). Different files, different shapes.
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 2 and fields[0] == "Nuclease":
                        guide_nuclease = fields[1]
                        break
            break
    guides_seen: List[str] = []
    for hap in ("paternal", "maternal"):
        for g in hap_predictions_all.get(hap, pd.DataFrame(columns=["Spacer+PAM"]))["Spacer+PAM"]:
            if g not in guides_seen:
                guides_seen.append(g)
    guide_summary_rows = []
    for guide in guides_seen:
        row = {"Guide": guide, "Nuclease": guide_nuclease}
        for hap, label in (("paternal", "Paternal"), ("maternal", "Maternal")):
            hap_df = hap_predictions_all.get(hap)
            hap_dir_name = params.get(f"{label}_dir")
            if hap_df is not None:
                row[f"Total off-targets ({label})"] = int(
                    (hap_df["Spacer+PAM"] == guide).sum()
                )
            else:
                row[f"Total off-targets ({label})"] = "?"
            if hap_dir_name and hap in hap_results_dirs:
                row[f"Specificity score ({label})"] = read_specificity_score(
                    hap_results_dirs[hap], hap_dir_name, [guide]
                )
            else:
                row[f"Specificity score ({label})"] = "CFD score not available"
        guide_summary_rows.append(row)

    # Mismatch+bulge breakdown, MAPPED sites only, same column shape as
    # complete-search's own "Summary by Mismatches/Bulges" table (Bulge
    # Type / Mismatches / Bulge Size / Reference / Variant / Combined) --
    # here Paternal/Maternal stand in for Reference/Variant, and Combined is
    # their sum, same as complete-search's own definition (checked directly
    # against a real summary_by_guide.*.txt file: Combined = Reference +
    # Variant, not a distinct-site count). One real difference worth being
    # upfront about: complete-search's Reference/Variant are mutually
    # exclusive categories of the SAME distinct-site list, so summing them
    # can't double-count a site. Paternal/Maternal here are two independent
    # per-genome searches, so a site reconciled as `origin: both` DOES
    # contribute to both the Paternal and the Maternal tally (and so to
    # Combined) -- Combined is "total (site, haplotype) observations at
    # this combo", not "distinct reconciled sites".
    #
    # "Mapped only" = restricted to each haplotype's own off_target_ids
    # that appear in the reconciled table `df` (df only ever contains
    # successfully-lifted-to-hg38 predictions) -- the haplotype-private
    # (non-mappable) sites already have their own tables above and are
    # deliberately excluded here, per request.
    mapped_ids: Dict[str, set] = {}
    for hap in ("paternal", "maternal"):
        id_col = f"off_target_id_{hap}"
        if id_col in df.columns:
            mapped_ids[hap] = set(
                df[id_col].dropna().astype(float).astype(int).astype(str)
            )
        else:
            mapped_ids[hap] = set()

    def _mmb_counts(hap: str) -> "pd.Series":
        hap_df = hap_predictions_all.get(hap)
        btype_col = "Bulge_type_(fewest_mm+b)"
        if hap_df is None or btype_col not in hap_df.columns:
            return pd.Series(dtype=int)
        mapped = hap_df[hap_df["off_target_id"].isin(mapped_ids.get(hap, set()))]
        if mapped.empty:
            return pd.Series(dtype=int)
        mm = pd.to_numeric(mapped["Mismatches_(fewest_mm+b)"], errors="coerce")
        b = pd.to_numeric(mapped["Bulges_(fewest_mm+b)"], errors="coerce")
        bt = mapped[btype_col]
        valid = mm.notna() & b.notna()
        if not valid.any():
            return pd.Series(dtype=int)
        return pd.Series(
            list(zip(bt[valid], mm[valid].astype(int), b[valid].astype(int)))
        ).value_counts()

    _pat_mmb, _mat_mmb = _mmb_counts("paternal"), _mmb_counts("maternal")
    mmb_rows = [
        {
            "Bulge Type": key[0],
            "Mismatches": key[1],
            "Bulge Size": key[2],
            "Paternal": int(_pat_mmb.get(key, 0)),
            "Maternal": int(_mat_mmb.get(key, 0)),
            "Combined": int(_pat_mmb.get(key, 0)) + int(_mat_mmb.get(key, 0)),
        }
        for key in sorted(set(_pat_mmb.index) | set(_mat_mmb.index), key=lambda k: (k[1], k[2], k[0]))
    ]
    # Same html.Table + rowSpan/colSpan grouped-header structure as
    # complete-search's own Summary by Mismatches/Bulges table
    # (pages_utils.py generate_table(): Bulge type/Mismatches/Bulge
    # Size/PAM Creation each row-span both header rows, "Targets found in
    # Genome" col-spans 3 sub-columns Reference/Variant/Combined) -- not a
    # dash_table.DataTable at all, which is the real structural difference
    # behind the "different visual layout" (generate_table() itself sets no
    # header color -- grepped, no style_header/CSS rule anywhere styles it,
    # so the layout, not a specific color, is the verifiable thing to
    # mirror here).
    _cell_style = {"vertical-align": "middle", "text-align": "center", "padding": "6px 10px"}
    mmb_summary_block = []
    if mmb_rows:
        mmb_header = [
            html.Tr(
                [
                    html.Th("Bulge Type", rowSpan="2", style=_cell_style),
                    html.Th("Mismatches", rowSpan="2", style=_cell_style),
                    html.Th("Bulge Size", rowSpan="2", style=_cell_style),
                    html.Th("Off-targets found (mapped)", colSpan="3", style=_cell_style),
                ]
            ),
            html.Tr(
                [html.Th(x, style=_cell_style) for x in ["Paternal", "Maternal", "Combined"]]
            ),
        ]
        mmb_body = [
            html.Tr(
                [
                    html.Td(r["Bulge Type"], style=_cell_style),
                    html.Td(r["Mismatches"], style=_cell_style),
                    html.Td(r["Bulge Size"], style=_cell_style),
                    html.Td(r["Paternal"], style=_cell_style),
                    html.Td(r["Maternal"], style=_cell_style),
                    html.Td(r["Combined"], style=_cell_style),
                ]
            )
            for r in mmb_rows
        ]
        mmb_summary_block = [
            html.H4("Mismatch + bulge breakdown (mapped sites only)"),
            html.P(
                "Same shape as complete-search's own Summary by "
                "Mismatches/Bulges table, with Paternal/Maternal standing "
                "in for Reference/Variant. Excludes the haplotype-private "
                "(non-mappable) sites shown separately above.",
                style={"font-size": "0.95rem", "color": "#777"},
            ),
            html.Table(
                mmb_header + mmb_body,
                id="assembly-mmb-summary-table",
                style={"display": "inline-block", "borderCollapse": "collapse"},
            ),
            html.Br(),
        ]

    def _stat(label: str, value) -> html.Div:
        return html.Div(
            [
                html.Div(str(value), style={"font-size": "1.8rem", "font-weight": "700"}),
                html.Div(label, style={"font-size": "1.0rem", "color": "#555"}),
            ],
            style={"display": "inline-block", "margin-right": "36px", "textAlign": "center"},
        )

    # Report downloads: each haplotype already has its own self-contained
    # report.zip (built automatically by that haplotype's complete-search
    # run, same generate_report.py used everywhere else) -- link both. A
    # combined (both-haplotype) report needs its own report-generation
    # logic, since the reconciled table's paternal/maternal-paired column
    # schema isn't something build_report() already understands -- deferred
    # pending a separate scoping pass, not implemented here. Positioned near
    # the top of the page, boxed like complete-search's own "Off-target
    # report" button, so it occupies the same prominent slot even though --
    # unlike that single combined button -- it's two links for now.
    zip_links = []
    for hap, label, dirname_key in (
        ("paternal", "Paternal", "Paternal_dir"),
        ("maternal", "Maternal", "Maternal_dir"),
    ):
        hap_dir_name = params.get(dirname_key)
        if not hap_dir_name:
            continue
        zip_name = f"{hap_dir_name}_report.zip"
        zip_path = os.path.join(current_working_directory, RESULTS_DIR, hap_dir_name, zip_name)
        if os.path.isfile(zip_path):
            zip_links.append(
                html.A(
                    f"⬇ {label} report (.zip)",
                    href=os.path.join(URL, RESULTS_DIR, hap_dir_name, zip_name),
                    target="_blank",
                    style={"marginRight": "24px", "fontWeight": "600", "color": "#fff"},
                )
            )
    report_box = []
    if zip_links:
        report_box = [
            html.Div(
                html.Div(
                    [
                        html.H4(
                            "📄 Off-target reports (per haplotype)",
                            style={"margin": "0 0 8px 0", "fontWeight": "700", "color": "#1a365d"},
                        ),
                        html.P(
                            "No single combined report yet -- each haplotype's own "
                            "full complete-search report, before reconciliation "
                            "against the other haplotype.",
                            style={"margin": "0 0 10px 0", "fontSize": "0.95rem", "color": "#334155"},
                        ),
                        html.Div(
                            zip_links,
                            style={"background": "#2b6cb0", "borderRadius": "8px", "padding": "10px"},
                        ),
                    ],
                    style={"textAlign": "center"},
                ),
                style={
                    "textAlign": "center",
                    "margin": "16px auto 24px auto",
                    "padding": "22px",
                    "background": "#eef6ff",
                    "border": "2px solid #2b6cb0",
                    "borderRadius": "12px",
                    "boxShadow": "0 2px 10px rgba(43,108,176,0.18)",
                },
            )
        ]

    # Per-haplotype guide summary table, complete-search's own Result
    # Summary table (Guide/Nuclease/aggregate specificity score/Total)
    # adapted to show paternal vs. maternal side by side -- same slot in
    # the layout (right after the report box, before the main table).
    # Nested/grouped header, mirroring complete-search's own Result Summary
    # table exactly (general-profile-table: merge_duplicate_headers=True,
    # columns given as [top-level, sub-level] name pairs, e.g.
    # ["Off-targets for Mismatch (MM) and Bulge (B) Value", "Total"]) --
    # here "Specificity score"/"Total off-targets" are the top-level groups,
    # Paternal/Maternal the sub-columns underneath, instead of their single
    # reference-genome value. (Their table's blue header tint isn't backed
    # by any style_header/CSS rule anywhere in this codebase -- grepped for
    # both and found neither, so the exact color can't be confirmed from
    # source; applying this page's own already-established blue instead of
    # guessing theirs.)
    guide_summary_block = []
    if guide_summary_rows:
        _gcol_defs = (
            [("", "Guide"), ("", "Nuclease")]
            + [("Specificity score", l) for l in ("Paternal", "Maternal")]
            + [("Total off-targets", l) for l in ("Paternal", "Maternal")]
        )
        guide_summary_block = [
            html.H4("Guide summary"),
            dash_table.DataTable(
                id="assembly-guide-summary-table",
                columns=[
                    {"name": [top, sub], "id": f"{top} {sub}".strip()}
                    for top, sub in _gcol_defs
                ],
                data=[
                    {f"{top} {sub}".strip(): row[f"{sub}" if not top else f"{top} ({sub})"]
                     for top, sub in _gcol_defs}
                    for row in guide_summary_rows
                ],
                merge_duplicate_headers=True,
                page_action="none",
                style_table={"overflowX": "auto"},
                style_header={
                    "backgroundColor": "#2b6cb0",
                    "color": "#fff",
                    "fontWeight": "700",
                    "textAlign": "center",
                },
                style_data={"whiteSpace": "normal", "height": "auto", "font-size": "1.15rem"},
                style_cell={"textAlign": "center", "padding": "6px"},
            ),
        ]
        # Plain-language burden callout -- "which haplotype carries more
        # off-target risk", computed directly from the row above (no new
        # data), one sentence per guide.
        _callouts = []
        for row in guide_summary_rows:
            pat_n = row.get("Total off-targets (Paternal)")
            mat_n = row.get("Total off-targets (Maternal)")
            if not isinstance(pat_n, int) or not isinstance(mat_n, int) or pat_n == mat_n:
                continue
            more, fewer, diff = (
                ("Maternal", "Paternal", mat_n - pat_n)
                if mat_n > pat_n
                else ("Paternal", "Maternal", pat_n - mat_n)
            )
            _callouts.append(
                html.Li(
                    f"{row['Guide']}: {more} haplotype has {diff} more off-target(s) "
                    f"than {fewer} ({pat_n} vs {mat_n})."
                )
            )
        if _callouts:
            guide_summary_block.append(
                html.Ul(_callouts, style={"font-size": "0.95rem", "color": "#334155"})
            )
        guide_summary_block.append(html.Br())

    # Shared layout applied to every plotly figure on this page -- "plotly_white"
    # instead of plotly's default (grey plot background, purple/blue default
    # trace colors) for a cleaner look closer to the rest of the page, plus a
    # font size matching the fontsize=17 convention complete-search's own
    # matplotlib images use (generate_img_radar_chart.py) so the interactive
    # and static plots on this page don't look like they're from two
    # different tools.
    _PLOTLY_LAYOUT = dict(template="plotly_white", font=dict(size=14), margin=dict(t=50))

    # Origin-split visual: a proportional horizontal bar rather than a
    # geometrically-precise Venn diagram -- a true Venn's circle-overlap
    # areas would need to be drawn to scale to not visually mislead about
    # the real proportions, which is real design/implementation work of its
    # own; a proportional bar shows the same three real counts accurately
    # with much less risk of that, so it's the safer honest default. Not
    # mapped-only -- deliberately covers the origin split which BY
    # DEFINITION only exists among mapped sites (non-mappable sites have no
    # origin category), so no filtering question applies here.
    origin_chart_block = []
    _origin_total = sum(origin_counts.get(k, 0) for k in ("both", "paternal_only", "maternal_only"))
    if _origin_total:
        _origin_fig = go.Figure()
        for key, label, color in (
            ("both", "Both haplotypes", "#2b6cb0"),
            ("paternal_only", "Paternal-only", "#63b3ed"),
            ("maternal_only", "Maternal-only", "#f6ad55"),
        ):
            _origin_fig.add_trace(
                go.Bar(
                    y=["Origin"],
                    x=[origin_counts.get(key, 0)],
                    name=f"{label} ({origin_counts.get(key, 0)})",
                    orientation="h",
                    marker_color=color,
                )
            )
        _origin_layout = {**_PLOTLY_LAYOUT, "margin": dict(t=70, b=40)}
        _origin_fig.update_layout(
            **_origin_layout,
            barmode="stack",
            height=220,
            xaxis_title="Reconciled off-target sites",
            yaxis=dict(visible=False),
            # legend above the plot, not below -- below collided with the
            # x-axis title in the same cramped margin (the previous cause of
            # the overlap).
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        )
        origin_chart_block = [dcc.Graph(figure=_origin_fig, id="assembly-origin-split-graph")]

    final_list = [
        html.H3(f"Personal Assembly Search Results — {job_id}"),
        html.P(
            f"Paternal genome: {params.get('Genome_paternal', '?')}  |  "
            f"Maternal genome: {params.get('Genome_maternal', '?')}  |  "
            f"PAM: {params.get('Pam', '?')}  |  "
            f"Mismatches: {params.get('Mismatches', '?')}  |  "
            f"DNA bulges: {params.get('DNA', '?')}  |  "
            f"RNA bulges: {params.get('RNA', '?')}",
            style={"color": "#555"},
        ),
        html.Hr(),
        *report_box,
        *guide_summary_block,
        html.H4("Haplotype coverage"),
        html.Div(
            [
                _stat("Found in both haplotypes", origin_counts.get("both", 0)),
                _stat("Paternal-only", origin_counts.get("paternal_only", 0)),
                _stat("Maternal-only", origin_counts.get("maternal_only", 0)),
                _stat(
                    "Paternal non-mappable to hg38",
                    summary_counts.get("paternal_non_mappable", "?"),
                ),
                _stat(
                    "Maternal non-mappable to hg38",
                    summary_counts.get("maternal_non_mappable", "?"),
                ),
            ]
        ),
        *origin_chart_block,
        html.P(
            "Non-mappable sites have no hg38 equivalent -- invisible to any "
            "reference-based search. Listed below (coordinates are in each "
            "haplotype's own assembly, not hg38, since these sites have no "
            "hg38 equivalent).",
            style={
                "font-size": "1.0rem",
                "color": "#777",
                "font-style": "italic",
                "marginTop": "6px",
            },
        ),
    ]
    for hap, label in (("paternal", "Paternal"), ("maternal", "Maternal")):
        private = private_tables.get(hap)
        if private is None or private.empty:
            continue
        final_list.append(
            html.Details(
                [
                    html.Summary(
                        f"{label}-private off-targets, no hg38 equivalent "
                        f"({len(private)} sites)",
                        style={"cursor": "pointer", "fontWeight": "600", "marginTop": "10px"},
                    ),
                    dash_table.DataTable(
                        id=f"assembly-private-table-{hap}",
                        columns=[{"name": c, "id": c, "hideable": True} for c in private.columns],
                        data=private.to_dict("records"),
                        page_size=25,
                        sort_action="native",
                        filter_action="native",
                        style_table={"overflowX": "auto"},
                        style_data={
                            "whiteSpace": "normal",
                            "height": "auto",
                            "font-size": "1.05rem",
                        },
                        style_cell={"textAlign": "left", "padding": "4px"},
                    ),
                ]
            )
        )
    # ---- "Custom Ranking" tab content: the main reconciled table ----
    custom_ranking_tab = [
        html.Br(),
        html.P(
            "hg38_chr/hg38_start (and strand) are the coordinates the two "
            "haplotypes were matched on, so they're shared between them by "
            "construction. hg38_end is lifted independently per haplotype "
            "and can differ by a few bp when an indel private to one "
            "haplotype shifts where its alignment ends in hg38 space -- "
            "hg38_end_paternal/hg38_end_maternal are both kept for that "
            "reason, rather than merged into one column.",
            style={
                "font-size": "1.0rem",
                "color": "#777",
                "font-style": "italic",
                "marginTop": "6px",
            },
        ),
        dash_table.DataTable(
            id="assembly-results-table",
            columns=[{"name": c, "id": c, "hideable": True} for c in display_cols],
            data=df.to_dict("records"),
            page_size=25,
            sort_action="native",
            filter_action="native",
            style_table={"overflowX": "auto"},
            style_data={"whiteSpace": "normal", "height": "auto", "font-size": "1.05rem"},
            style_cell={"textAlign": "left", "padding": "4px"},
        ),
    ]

    # ---- "Graphical Reports" tab content ----
    # Top-1000 image: reuse each haplotype's own pre-rendered CFD plot
    # (already generated by that haplotype's underlying complete-search
    # run, nothing new computed here). Only CFD -- not the CRISTA or fewest
    # mm+b variants (per request), and not the "_by_variant_effect" sibling,
    # which splits/colors by REF-vs-ALT variant origin -- always degenerate
    # here (no VCF, so always "ref"), same reasoning as dropping the ALT
    # table columns above.
    graphical_blocks = []
    for hap, label, dirname_key in (
        ("paternal", "Paternal", "Paternal_dir"),
        ("maternal", "Maternal", "Maternal_dir"),
    ):
        hap_dir_name = params.get(dirname_key)
        if not hap_dir_name:
            continue
        imgs_dir = os.path.join(current_working_directory, RESULTS_DIR, hap_dir_name, "imgs")
        img_paths = sorted(glob(os.path.join(imgs_dir, "CRISPRme_CFD_top_1000_log_for_main_text_*.png")))
        encoded = [p for p in ((path, _encode_png(path)) for path in img_paths) if p[1]]
        if not encoded:
            continue
        graphical_blocks.append(
            html.Div(
                [html.H6(label)]
                + [
                    html.Img(src=src, style={"maxWidth": "420px", "margin": "6px"})
                    for _, src in encoded
                ],
                style={"display": "inline-block", "verticalAlign": "top", "marginRight": "24px"},
            )
        )

    # Paternal-vs-maternal CFD-score distribution: same area-plot shape as
    # complete-search's own CFDGraph (PostProcess/CFDGraph.py, REF vs VAR
    # series), but that module's legend labels are hardcoded strings inside
    # createGraph() (not parameters), so it can't be reused as-is for a
    # paternal/maternal split without editing that file -- which per-branch
    # convention this work doesn't touch. New, independent figure instead:
    # same plotly area-plot approach, correct labels, built from data this
    # page already loaded above (hap_predictions_all), not a new file format.
    # CFD is natively 0-1 (same scale shown everywhere else on this page,
    # e.g. the CFD_score_(fewest_mm+b) column) -- bucketed at 0.01 resolution
    # (101 buckets) for a reasonably smooth curve, but the axis stays on the
    # native 0-1 scale rather than rescaled to 0-100.
    def _cfd_bucket_counts(hap: str) -> List[int]:
        counts = [0] * 101
        hap_df = hap_predictions_all.get(hap)
        if hap_df is None or "CFD_score_(fewest_mm+b)" not in hap_df.columns:
            return counts
        for v in pd.to_numeric(hap_df["CFD_score_(fewest_mm+b)"], errors="coerce").dropna():
            counts[min(100, max(0, int(round(v * 100))))] += 1
        return counts

    if hap_predictions_all:
        cfd_fig = go.Figure()
        for hap, label in (("paternal", "Paternal"), ("maternal", "Maternal")):
            cfd_fig.add_trace(
                go.Scatter(
                    x=[i / 100 for i in range(101)],
                    y=_cfd_bucket_counts(hap),
                    fill="tozeroy",
                    name=label,
                )
            )
        cfd_fig.update_layout(
            **_PLOTLY_LAYOUT,
            xaxis_title="CFD score (0-1)",
            # dtick=1 on a log axis means "one tick per decade" (1, 10, 100,
            # ...) -- without it plotly's auto-ticking on this data (many
            # zero/near-zero buckets next to a handful of large ones) mixes
            # tick spacings inconsistently, which is what looked "weird".
            yaxis=dict(title="Number of off-targets (log scale)", type="log", dtick=1),
            hovermode="x",
        )
        graphical_blocks.append(
            html.Div(
                [
                    html.H6("CFD-score distribution, paternal vs. maternal"),
                    dcc.Graph(figure=cfd_fig, id="assembly-cfd-distribution-graph"),
                ]
            )
        )

    # Per-position mismatch/bulge distribution -- what complete-search's own
    # Graphical Reports tab actually shows for "what basepair it is
    # mismatched/bulged to around the guide": the bottom bar-chart panel of
    # generate_img_radar_chart.py's combined ENCODE/GENCODE+motif figure
    # (PostProcess/generate_img_radar_chart.py:258-283, driven by a
    # motifDict built in PostProcess/radar_chart_dict_generator.py from the
    # raw per-alignment aligned strings: lowercase letter = mismatch-to-that-
    # base, "-" = bulge). Not reused directly -- that combined figure bundles
    # an ENCODE/GENCODE annotation radar chart in the SAME image, which
    # would render degenerate here (this run's annotation is "vuoto.txt" /
    # empty, so every Annotation_* column is "NA" -- see the annotation
    # column decision elsewhere on this page), and its motifDict is built
    # from raw CRISPRitz target-file columns this page doesn't load. New,
    # independent tally instead, against the SAME aligned-sequence columns
    # this page already has (Aligned_protospacer+PAM_REF_(fewest_mm+b) for
    # mismatches/RNA-bulges, Aligned_spacer+PAM_(fewest_mm+b) for DNA-bulges
    # -- confirmed directly against real rows: lowercase = mismatch base,
    # "-" on the REF side = RNA bulge, "-" on the spacer side = DNA bulge,
    # same convention). Position is the index within the aligned strings
    # (which can be one longer than the guide itself when there's a DNA
    # bulge), not re-derived against the original script's guide-relative
    # edge-case handling for N-padded PAMs -- labeled "alignment position",
    # not literal guide bases, to not overclaim exactness. Built from each
    # haplotype's full prediction set (not mapped-only): this is a per-
    # haplotype sequence-composition statistic, not tied to hg38 mappability.
    def _position_tallies(hap: str) -> Dict[str, List[int]]:
        hap_df = hap_predictions_all.get(hap)
        ref_col = "Aligned_protospacer+PAM_REF_(fewest_mm+b)"
        spacer_col = "Aligned_spacer+PAM_(fewest_mm+b)"
        if hap_df is None or ref_col not in hap_df.columns or spacer_col not in hap_df.columns:
            return {}
        ref_seqs = hap_df[ref_col].dropna().astype(str)
        spacer_seqs = hap_df[spacer_col].dropna().astype(str)
        width = max([len(s) for s in ref_seqs] + [len(s) for s in spacer_seqs] + [0])
        if width == 0:
            return {}
        tallies = {k: [0] * width for k in ("A", "C", "G", "T", "RNA bulge", "DNA bulge")}
        for seq in ref_seqs:
            for i, ch in enumerate(seq):
                if ch.islower() and ch.upper() in "ACGT":
                    tallies[ch.upper()][i] += 1
                elif ch == "-":
                    tallies["RNA bulge"][i] += 1
        for seq in spacer_seqs:
            for i, ch in enumerate(seq):
                if ch == "-":
                    tallies["DNA bulge"][i] += 1
        return tallies

    # x-axis labels: bare guide letters, matching the original exactly
    # (plt.xticks(ticks=ind, labels=list(guide))) -- an earlier pass here
    # added a "position:" prefix to disambiguate the guide's repeated
    # letters, which is NOT what the original does; reverted to match.
    # y-axis: the original normalizes each position's stacked values by
    # that position's own total-across-all-series divided into the GLOBAL
    # max such total (generate_img_radar_chart.py:243-252: `maxmax =
    # max(totalMotif)`, `motifDict[nuc][count] /= maxmax`) -- i.e. every
    # bar's height is relative to the single most-covered position, 0-1,
    # not a raw count. Re-read that code and initially didn't carry this
    # through; fixed here to match.
    _guide_seq = guides_seen[0] if guides_seen else ""
    position_figs = []
    for hap, label in (("paternal", "Paternal"), ("maternal", "Maternal")):
        tallies = _position_tallies(hap)
        if not tallies:
            continue
        width = len(next(iter(tallies.values())))
        x = [_guide_seq[i] if i < len(_guide_seq) else "" for i in range(width)]
        totals_per_position = [sum(tallies[s][i] for s in tallies) for i in range(width)]
        maxmax = max(totals_per_position) if totals_per_position else 0
        normalized = (
            {s: [v / maxmax for v in tallies[s]] for s in tallies}
            if maxmax
            else tallies
        )
        fig = go.Figure()
        for series in ("A", "C", "G", "T", "RNA bulge", "DNA bulge"):
            fig.add_trace(go.Bar(x=x, y=normalized[series], name=series))
        fig.update_layout(
            **_PLOTLY_LAYOUT,
            barmode="stack",
            title=label,
            xaxis_title="Guide position",
            yaxis=dict(title="Fraction of most-covered position's total", range=[0, 1]),
        )
        # stacked in one column (not side by side) so each is bigger/more
        # readable, per request.
        position_figs.append(
            html.Div(
                dcc.Graph(figure=fig, id=f"assembly-position-mmb-barplot-{hap}"),
                style={"width": "100%", "maxWidth": "900px"},
            )
        )
    if position_figs:
        graphical_blocks.append(
            html.Div(
                [
                    html.H6("Mismatch/bulge distribution by alignment position"),
                    html.P(
                        "Per-position tally of which base an off-target "
                        "mismatches to, and where RNA/DNA bulges occur, "
                        "across all off-targets found in that haplotype.",
                        style={"font-size": "0.95rem", "color": "#777"},
                    ),
                ]
                + position_figs
            )
        )

    graphical_reports_tab = (
        [html.Br(), html.Div(graphical_blocks)]
        if graphical_blocks
        else [html.Br(), html.P("No graphical reports available for this job.")]
    )

    # Own tab, matching complete-search's actual "Summary by
    # Mismatches/Bulges" tab (moved out of the always-visible top section --
    # an earlier pass put it there since only 2 tabs were planned at the
    # time; now a 3rd tab matching complete-search's real tab set instead).
    mmb_tab = (
        [html.Br()] + mmb_summary_block
        if mmb_summary_block
        else [html.Br(), html.P("No mismatch+bulge summary available for this job.")]
    )

    final_list.append(html.Hr())
    final_list.append(
        dcc.Tabs(
            [
                dcc.Tab(label="Custom Ranking", children=custom_ranking_tab),
                dcc.Tab(label="Summary by Mismatches/Bulges", children=mmb_tab),
                dcc.Tab(label="Graphical Reports", children=graphical_reports_tab),
            ]
        )
    )

    return html.Div(final_list, style={"margin": "1%"})


# store drop-down value in auxiliary file
@app.callback(
    Output("store", "data"),
    [Input("target_filter_dropdown", "value")],
    [State("url", "search")],
)
def sendto_write_json(filter_criterion: str, search: str) -> None:
    """Write auxiliary file to store the table filtering criterion
    (received from the drop-down) and filter the tables displayed in
    Summary by Mismatches/Bulges accordingly.

    The function is triggered by the user, when choosing the filtering
    criterion from the drop-down bar.

    ...

    Parameters
    ----------
    filter_criterion : str
        Table filtering criterion
    search : str
        Target search name

    Returns
    -------
    None
    """
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not filter_criterion in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    job_id = search.split("=")[-1]
    write_json(filter_criterion, job_id)


# -------------------------------------------------------------------------------
# Download links generation and actions definition
#


# Generate download link summary_by_sample
@app.callback(
    [
        Output("download-link-summary_by_sample", "children"),
        Output("interval-summary_by_sample", "disabled"),
    ],
    [Input("interval-summary_by_sample", "n_intervals")],
    [State("div-info-summary_by_sample", "children"), State("url", "search")],
)
def download_link_sample(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load =
    """Create the link to download CRISPRme result files.

    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate  # nothing to do
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "txt"])
    file_to_load = file_to_load.strip().split("/")[-1]
    # print(file_to_load)
    if os.path.exists(os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)):
        return (
            html.A(
                "Download file",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    # The interval fires once per second; if the file has not appeared after
    # ~30s it is not going to, so stop polling instead of showing the
    # "Generating download link..." banner forever.
    if n >= 30:
        return (
            "Download link unavailable (the file could not be generated).",
            True,
        )
    return "Generating download link, Please wait...", False


# download the self-contained off-target report (replaces the former three
# raw-file links: general table, integrated results, alt-alignments). The report
# ZIP already bundles the integrated results, per-tier subsets, and the complete
# raw table, so a single button supersedes all three. The report is generated at
# the end of the job; the polling interval reveals the button once it exists.
@app.callback(
    [
        Output("download-link-report", "children"),
        Output("interval-report", "disabled"),
    ],
    [Input("interval-report", "n_intervals")],
    [State("div-info-report", "children"), State("url", "search")],
)
def download_report(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:
    """Reveal the 'Download Report' button once <job>_report.zip exists."""
    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = file_to_load.split("/")[-1]
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)
    ):
        return (
            html.A(
                "⬇  Download Report (.zip)",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
                style={
                    "display": "inline-block",
                    "background": "#2f855a",
                    "color": "#fff",
                    "padding": "16px 44px",
                    "borderRadius": "8px",
                    "textDecoration": "none",
                    "fontWeight": "700",
                    "fontSize": "1.4rem",
                    "boxShadow": "0 3px 10px rgba(0,0,0,0.25)",
                    "marginTop": "4px",
                },
            ),
            True,
        )
    return "Preparing the off-target report for download…", False


# Generate download link sumbysample
@app.callback(
    [
        Output("download-link-sumbysample", "children"),
        Output("interval-sumbysample", "disabled"),
    ],
    [Input("interval-sumbysample", "n_intervals")],
    [State("div-info-sumbysample-targets", "children"), State("url", "search")],
)
def download_link_sample(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load = job_id.HG001.guide
    """Create the link to download CRISPRme results by sample table.

    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "zip"])
    if os.path.exists(os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)):
        return (
            html.A(
                "Download zip",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# Generate download link sumbyguide
@app.callback(
    [
        Output("download-link-sumbyguide", "children"),
        Output("interval-sumbyguide", "disabled"),
    ],
    [Input("interval-sumbyguide", "n_intervals")],
    [State("div-info-sumbyguide-targets", "children"), State("url", "search")],
)
def downloadLinkGuide(
    n: int, file_to_load: str, search: str
) -> Tuple[str, bool]:  # file to load = job_id.RNA.1.0.guide
    """Create the link to download CRISPRme results by sample table.

    ...

    Parameters
    ----------
    n : int
    file_to_load : str
        File to download
    search : str
        Target search name

    Returns
    -------
    str
    bool
    """

    if not isinstance(file_to_load, str):
        raise TypeError(f"Expected {str.__name__}, got {type(file_to_load).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if n is None:
        raise PreventUpdate
    job_id = search.split("=")[-1]
    file_to_load = ".".join([file_to_load, "zip"])
    if os.path.exists(os.path.join(current_working_directory, RESULTS_DIR, job_id, file_to_load)):
        return (
            html.A(
                "Download zip",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
        )
    return "Generating download link, Please wait...", False


# trigger file download
@app.server.route("/Results/<path:path>")
def download_file(path: str) -> flask.Response:
    """Download the chosen file.

    ...

    Parameters
    ----------
    path : str
        Path to file location

    Returns
    -------
    flask.Response
    """

    if not isinstance(path, str):
        raise TypeError(f"Expected {str.__name__}, got {type(path).__name__}")
    # print(current_working_directory)
    # print('test', path)
    return flask.send_from_directory(
        os.path.join(current_working_directory, "Results/"), path, as_attachment=True
    )


# Filter/sort IUPAC decomposition table for cluster page
@app.callback(
    Output("table-scomposition-cluster", "data"),
    [
        Input("table-scomposition-cluster", "page_current"),
        Input("table-scomposition-cluster", "page_size"),
        Input("table-scomposition-cluster", "sort_by"),
        Input("table-scomposition-cluster", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_iupac_decomposition_table_cluster(
    page_current: int,
    page_size: int,
    filter_criterion: str,
    search: str,
    hash_term: str,
) -> Dict[str, str]:
    """

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    filter_criterion : str
        Data table filter
    search : str
        Unique search ID
    hash_term : str
        Hashing

    Returns
    -------
    Dict[str, str]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    hash_term = hash_term.split("#")[1]
    guide = hash_term[: hash_term.find("-Pos-")]
    chr_pos = hash_term[(hash_term.find("-Pos-") + 5) :]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    try:
        with open(os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    if genome_type == "ref":
        raise PreventUpdate
    filtering_expressions = filter_criterion.split(" && ")
    decomp_fname = (
        job_id + "." + chromosome + "_" + position + "." + guide + ".scomposition.txt"
    )
    # load data and cache the data table (in pd.DataFrame)
    df_cached = global_store_general(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, decomp_fname)
    )
    if df_cached is None:  # nothing to display and do not update the page
        raise PreventUpdate
    df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    # filter data table
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in PANDAS_OPERATORS:
            # these operators match pandas series operator method names
            df_cached = df_cached.loc[
                getattr(df_cached[col_name], operator)(filter_value)
            ]
        elif operator == "contains":
            df_cached = df_cached.loc[df_cached[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            df_cached = df_cached.loc[df_cached[col_name].str.startswith(filter_value)]
    # Calculate sample count
    data_to_send = df_cached.iloc[
        page_current * page_size : (page_current + 1) * page_size
    ].to_dict("records")
    return data_to_send


@app.callback(
    Output("table-position-target", "data"),
    [
        Input("table-position-target", "page_current"),
        Input("table-position-target", "page_size"),
        Input("table-position-target", "sort_by"),
        Input("table-position-target", "filter_query"),
        Input("hide-reference-targets", "value"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table_cluster(
    page_current: int,
    page_size: int,
    sort_by: List[str],
    filter_criterion: str,
    hide_reference: str,
    search: str,
    hash_term: str,
) -> Dict[str, str]:
    """

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    sort_by : List[str]
        Columns used while sorting the data table
    filter_criterion : str
        Data table filter
    hide_reference : str
        Hide reference data
    search : str
        Unique search ID
    has_term : str
        Hashing

    Returns
    -------
    Dict[str, str]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Exepcted {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(sort_by, list):
        raise TypeError(f"Expected {list.__name__}, got {type(sort_by).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Exepcted {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not isinstance(hide_reference, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hide_reference).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    hash_term = hash_term.split("#")[1]
    guide = hash_term[: hash_term.find("-Pos-")]
    chr_pos = hash_term[hash_term.find("-Pos-") + 5 :]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    try:
        with open(os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    filtering_expressions = filter_criterion.split(" && ")
    guide_fname = job_id + "." + chromosome + "_" + position + "." + guide + ".txt"
    # cache guide data table
    df_cached = global_store_general(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, guide_fname)
    )
    if df_cached is None:  # empty file -> nothing cached and nothing to do
        raise PreventUpdate
    if genome_type == "ref":
        df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    else:
        df_cached.rename(columns=COL_BOTH_RENAME, inplace=True)
    # drop unused columns
    if "hide-ref" in hide_reference or genome_type == "var":
        df_cached.drop(df_cached[(df_cached["Samples"] == "n")].index, inplace=True)
    # hide reference data
    if "hide-cluster" in hide_reference:
        df_cached = df_cached.head(1)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in PANDAS_OPERATORS:
            # these operators match pandas series operator method names
            df_cached = df_cached.loc[
                getattr(df_cached[col_name], operator)(filter_value)
            ]
        elif operator == "contains":
            df_cached = df_cached.loc[df_cached[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            df_cached = df_cached.loc[df_cached[col_name].str.startswith(filter_value)]
    # sort data table by the defined columns
    if bool(sort_by):
        df_cached = df_cached.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )
    # Calculate sample count
    data_to_send = df_cached.iloc[
        (page_current * page_size) : ((page_current + 1) * page_size)
    ].to_dict("records")
    if genome_type != "ref":
        (
            dict_sample_to_pop,
            dict_pop_to_superpop,
        ) = associateSample.loadSampleAssociation(
            os.path.join(job_directory, SAMPLES_ID_FILE)
        )[
            :2
        ]
        for row in data_to_send:
            summarized_sample_cell = {}
            for s in row["Samples"].split(","):
                if s == "n":
                    break
                try:
                    summarized_sample_cell[
                        dict_pop_to_superpop[dict_sample_to_pop[s]]
                    ] += 1
                except:
                    summarized_sample_cell[
                        dict_pop_to_superpop[dict_sample_to_pop[s]]
                    ] = 1
            if summarized_sample_cell:
                row["Samples Summary"] = ", ".join(
                    [
                        str(summarized_sample_cell[sp]) + " " + sp
                        for sp in summarized_sample_cell
                    ]
                )
            else:
                row["Samples Summary"] = "n"
    return data_to_send


def cluster_page(job_id: str, hash_term: str) -> html.Div:
    """Recover CRISPR targets for the selected cluster.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    hash_term : str
        Hashing

    Returns
    -------
    html.Div
        Sample page layout
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    guide = hash_term[: hash_term.find("-Pos-")]
    chr_pos = hash_term[(hash_term.find("-Pos-") + 5) :]
    chromosome = chr_pos.split("-")[0]
    position = chr_pos.split("-")[1]
    if not os.path.isdir(os.path.join(current_working_directory, RESULTS_DIR, job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (next(s for s in params.split("\n") if "Ref_comp" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    style_hide_reference = {"display": "none"}  # display reference data
    value_hide_reference = []
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
        style_hide_reference = {}
        value_hide_reference = ["hide-ref", "hide-cluster"]  # hide reference data
    # begin page body construction
    final_list = []  # HTML page handler
    assert isinstance(chromosome, str)
    assert isinstance(position, str)
    final_list.append(html.H3(f"Selected Position: {chromosome} - {position}"))
    if genome_type == "ref":
        cols = [
            {"name": i, "id": i, "type": t, "hideable": True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)
        ]
        file_to_grep = ".bestMerge.txt"
    else:
        cols = [
            {"name": i, "id": i, "type": t, "hideable": True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)
        ]
        file_to_grep = ".bestMerge.txt"
    cluster_grep_result = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        ".".join([job_id, f"{chromosome}_{position}", guide, "txt"]),
    )
    put_header_cmd = " ".join(
        [
            "head -1",
            os.path.join(
                current_working_directory,
                RESULTS_DIR,
                job_id,
                f".{job_id}{file_to_grep}",
            ),
            f"> {cluster_grep_result} ; ",
        ]
    )
    # Example    job_id.chr3_100.guide.txt
    if not os.path.exists(cluster_grep_result):
        # os.system(f'touch {cluster_grep_result}')
        # Grep annotation for ref
        cmd = f"head -1 {file_to_grep} > {cluster_grep_result}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
        if genome_type == "ref":  # NOTE HEADER NOT SAVED
            cmd = " ".join(
                [
                    "grep -F",
                    guide,
                    os.path.join(
                        current_working_directory,
                        RESULTS_DIR,
                        job_id,
                        f"{job_id}.Annotation.targets.txt",
                    ),
                    "|",
                    f"awk '$6=={position} && $4==\"{chromosome}\"'",
                ]
            )
            get_annotation = subprocess.Popen(
                [cmd],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            out, err = get_annotation.communicate()
            annotation_type = out.decode("UTF-8").strip().split("\t")[-1]
            os.popen(
                put_header_cmd
                + " grep -F "
                + guide
                + " "
                + current_working_directory
                + "Results/"
                + job_id
                + "/"
                + job_id
                + file_to_grep
                + " | awk '$6=="
                + position
                + ' && $4=="'
                + chromosome
                + '" {###print $0"\\t'
                + annotation_type
                + "\"}' >> "
                + cluster_grep_result
            ).read()
        else:  # NOTE HEADER NOT SAVED
            os.popen(
                " ".join(
                    [
                        put_header_cmd,
                        "grep -F",
                        guide,
                        os.path.join(
                            current_working_directory,
                            RESULTS_DIR,
                            job_id,
                            f"{job_id}{file_to_grep}",
                        ),
                        "|",
                        f"awk '$6=={position} && $4==\"{chromosome}\"'",
                        ">>",
                        cluster_grep_result,
                    ]
                )
            ).read()
            # NOTE top1 will have sample and annotation, other targets will
            # have '.'-> 18/03 all samples and annotation are already writter
            # for all targets

        # TODO: review this part
        os.system(
            f"python {app_directory}/PostProcess/change_headers_bestMerge.py {cluster_grep_result} {cluster_grep_result}.tmp"
        )
        os.system(
            f"mv -f {cluster_grep_result}.tmp {cluster_grep_result} > /dev/null 2>&1"
        )
        # zip cluster results
        cmd = f"zip -j {cluster_grep_result.replace('txt', 'zip')} {cluster_grep_result} &"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
    final_list.append(
        html.Div(
            f"{job_id}.{chromosome}_{position}.{guide}",
            style={"display": "none"},
            id="div-info-sumbyposition-targets",
        )
    )
    decomp_fname = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.{chromosome}_{position}.{guide}.scomposition.txt",
    )
    iupac_decomp_visibility = {"display": "none"}
    if genome_type != "ref":
        iupac_decomp_visibility = {}
        # Example    job_id.chr_pos.guide.scomposition.txt
        # if not os.path.exists(scomposition_file):
        # os.system(f'touch {scomposition_file}')
        cmd = " ".join(
            [
                "grep -F",
                guide,
                os.path.join(
                    current_working_directory,
                    RESULTS_DIR,
                    job_id,
                    f".{job_id}{file_to_grep}",
                ),
                "|",
                f'awk \'$6=={position} && $4=="{chromosome}" && $13!="n"\'',
                ">",
                decomp_fname,
            ]
        )
        os.popen(cmd).read()
    final_list.append(
        html.P(
            [
                html.P(
                    "List of all the configurations for the target in the selected position.",
                    style=iupac_decomp_visibility,
                ),
                dcc.Checklist(
                    options=[
                        {"label": "Hide Reference Targets", "value": "hide-ref"},
                        {"label": "Show only TOP1 Target", "value": "hide-cluster"},
                    ],
                    id="hide-reference-targets",
                    value=value_hide_reference,
                    style=style_hide_reference,
                ),
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbyposition",
                        ),
                        dcc.Interval(interval=5 * 1000, id="interval-sumbyposition"),
                    ]
                ),
            ]
        )
    )
    cols_for_decomp = cols.copy()
    cols_for_decomp.append(
        {"name": "Samples", "id": "Samples", "type": "text", "hideable": True}
    )
    final_list.append(
        html.Div(
            dash_table.DataTable(
                # Table storing IUPAC decomposition of the selected target
                # rows are recovered from top1.samples.txt
                id="table-scomposition-cluster",
                columns=cols_for_decomp,
                virtualization=True,
                fixed_rows={"headers": True, "data": 0},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                sort_action="custom",
                sort_mode="multi",
                sort_by=[],
                filter_action="custom",
                filter_query="",
                style_table={"max-height": "600px"},
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            style=iupac_decomp_visibility,
        )
    )
    final_list.append(html.Hr())
    # Build cluster Table
    final_list.append(
        # if rows are highlighted in red, the target was found only in
        # non-reference genome (enriched with variants)
        str(
            "List of Targets found for the selected position. Other possible "
            "configurations of the target are listed in the table above, along "
            "with the corresponding samples list."
        ),
    )
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-position-target",
                columns=cols,
                virtualization=True,
                fixed_rows={"headers": True, "data": 0},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                sort_action="custom",
                sort_mode="multi",
                sort_by=[],
                filter_action="custom",
                filter_query="",
                style_table={
                    "max-height": "600px",
                    "overflowY": "scroll",
                },
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            id="div-result-table",
        )
    )
    return html.Div(final_list, style={"margin": "1%"})


# -------------------------------------------------------------------------------
# Summary by Sample tab
#


def global_get_sample_targets(
    job_id: str, sample: str, guide: str, page: int
) -> pd.DataFrame:
    """Recover CRISPRme analysis report regarding the selected sample.
    The sample related report can be filtered using the criteria available
    in the drop-down bar, above the report tabs:
    - CFD score
    - CRISTA score
    - Fewest Mismatches and Bulges

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    sample : str
        Sample identifier
    guide : str
        CRISPR guide
    page : int
        Current page

    Returns
    -------
    pd.DataFrame
        Data table reporting CRISPRme analysis results related to the
        selected sample
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(sample, str):
        raise TypeError(f"Expected {str.__name__}, got {type(sample).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if job_id is None:
        return ""
    db_path = glob(os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db"))[0]
    assert isinstance(db_path, str)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    filter_criterion = read_json(job_id)  # recover filter criterion selected
    query_cols = get_query_column(filter_criterion)
    # query the db
    result = pd.read_sql_query(
        "SELECT * FROM final_table WHERE \"{}\"='{}' AND \"{}\" LIKE '%{}%' LIMIT {} OFFSET {}".format(
            GUIDE_COLUMN,
            guide,
            query_cols["samples"],
            sample,
            PAGE_SIZE,
            page * PAGE_SIZE,
        ),
        conn,
    )
    return result


# callback to update the samples table
@app.callback(
    [Output("table-sample-target", "data"), Output("table-sample-target", "columns")],
    [
        Input("table-sample-target", "page_current"),
        Input("table-sample-target", "page_size"),
        Input("table-sample-target", "sort_by"),
        Input("table-sample-target", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table_sample(
    page_current: int,
    page_size: int,
    sort_by: str,
    filter_criterion: str,
    search: str,
    hash_term: str,
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Update the sample table accordingly to the filtering criterion
    selected in the drop-down bar.

    ...

    Parameters
    ----------
    page_current : int
        Current webpage
    page_size : int
        Webpage size
    sort_by : str
        Data table sorting criterion
    filter_criterion : str
        Data table filtering criterion
    search : str
        Search identifier
    hash_term : str
        Hashing term

    Returns
    -------
    Tuple[Dict[str, str], pd.DataFrame]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    job_id = search.split("=")[-1]
    filter_criterion = read_json(job_id)  # recover filter criterion
    assert isinstance(filter_criterion, str)
    assert filter_criterion in FILTERING_CRITERIA
    hash_term = hash_term.split("#")[1]
    guide = hash_term[: hash_term.find("-Sample-")]
    sample = str(hash_term[hash_term.rfind("-") + 1 :])
    try:
        with open(os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)) as handle:
            all_params = handle.read()
            genome_type_f = (
                next(s for s in all_params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (
                next(s for s in all_params.split("\n") if "Ref_comp" in s)
            ).split("\t")[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    # populate the sample table
    sample_df = global_get_sample_targets(job_id, sample, guide, page_current)
    # filter the sample table
    drop_cols = drop_columns(sample_df, filter_criterion)
    sample_df.drop(drop_cols, inplace=True, axis=1)
    # personal targets report filename
    integrated_sample_personal_fname = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        ".".join([job_id, sample, guide, "personal_targets.tsv"]),
    )
    # store sample table to personal targets file
    sample_df.to_csv(
        integrated_sample_personal_fname, sep="\t", na_rep="NA", index=False
    )
    # personal targets report ZIP
    integrated_sample_personal_zip_fname = integrated_sample_personal_fname.replace(
        "tsv", "zip"
    )
    # zip operation, non blocking
    cmd = f"zip -j {integrated_sample_personal_zip_fname} {integrated_sample_personal_fname} &"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f'An error occurred while running "{cmd}"')
    columns_df = [
        {"name": i, "id": i, "hideable": True}
        for col, i in enumerate(sample_df.columns)
    ]
    return sample_df.to_dict("records"), columns_df


# Return the targets found for the selected sample


def sample_page(job_id: str, hash_term: str) -> html.Div:
    """Build the sample webpage.
    The sample page contains the CRISPR targets found for the selected
    sample.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    hash_term : str
        Hashing

    Returns
    -------
    html.Div
        Sample webpage
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(hash_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_term).__name__}")
    guide = hash_term[: hash_term.find("-Sample-")]
    sample = str(hash_term[(hash_term.rfind("-") + 1) :])
    if not os.path.isdir(os.path.join(current_working_directory, RESULTS_DIR, job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (next(s for s in params.split("\n") if "Ref_comp" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    # begin sample page construction
    final_list = []  # HTML page handler
    final_list.append(html.H3(f"Selected Sample: {sample}"))  # page header
    final_list.append(
        html.P(
            [
                # if rows are highlghted in red, the CRISPR target was found
                # only in non reference genome (enriched with variants)
                "List of Targets found for the selected sample.",
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbysample",
                        ),
                        dcc.Interval(interval=(5 * 1000), id="interval-sumbysample"),
                    ]
                ),
            ]
        )
    )
    # header file
    header = os.path.join(current_working_directory, RESULTS_DIR, job_id, "header.txt")
    # file_to_grep = current_working_directory + 'Results/' + \
    #     job_id + '/.' + job_id + '.bestMerge.txt'
    integrated_fname = glob(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, "*integrated*")
    )[
        0
    ]  # take the first element
    assert isinstance(integrated_fname, str)
    file_to_grep = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.bestMerge.txt.integrated_results.tsv",
    )
    sample_grep_result = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, f"{job_id}.{sample}.{guide}.txt"
    )
    integrated_sample_personal = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.{sample}.{guide}.personal_targets.tsv",
    )
    integrated_sample_personal_zip = integrated_sample_personal.replace("tsv", "zip")
    final_list.append(
        html.Div(
            f"{job_id}.{sample}.{guide}.personal_targets",
            style={"display": "none"},
            id="div-info-sumbysample-targets",
        )
    )
    # define path to db
    db_path = glob(os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db"))[0]
    assert isinstance(db_path, str)
    # initialize db for queries
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    total_private_sample = f"SELECT * FROM final_table LIMIT 1"
    rows = c.execute(total_private_sample)
    header = [description[0] for description in rows.description]
    conn.commit()
    conn.close()  # close db connection
    # define columns
    cols = [{"name": i, "id": i, "hideable": True} for i in header]
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-sample-target",
                columns=cols,
                style_cell={"textAlign": "left"},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                style_table={
                    "overflowX": "scroll",
                    "overflowY": "scroll",
                    "max-height": "300px",
                },
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
            ),
            id="div-result-table",
        )
    )
    return html.Div(final_list, style={"margin": "1%"})


# TODO: move auxiliary functions close to each other in this file


@cache.memoize()
def global_store_general(path_file_to_load: str) -> pd.DataFrame:
    """Cache target files to improve results visualization and get better
    performances.

    ...

    Parameters
    ----------
    path_file_to_load : str
        Path to file to cache

    Returns
    -------
    pandas.DataFrame
        Results table
    """

    if not isinstance(path_file_to_load, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(path_file_to_load).__name__}"
        )
    if path_file_to_load is not None and not os.path.isfile(path_file_to_load):
        raise FileNotFoundError(f"Unable to locate {path_file_to_load}")
    if path_file_to_load is None:
        return ""  # do not cache anything
    if "scomposition" in path_file_to_load:
        rows_to_skip = 1
    else:
        rows_to_skip = 1  # Skip header
    # make sure file to cache is not empty
    if os.path.getsize(path_file_to_load) > 0:
        # TSV format -> sep="\t"
        df = pd.read_csv(path_file_to_load, sep="\t", index_col=False, na_filter=False)
    else:
        df = None  # empty file, no need for caching
    return df


# -------------------------------------------------------------------------------
# Summary by Mismatches/Bulges tab
#


# Update primary table of 'Show targets' of Summary by Mismatches/Bulges
@app.callback(
    [
        Output("table-subset-target", "data"),
        Output("table-subset-target", "columns"),
    ],
    [
        Input("table-subset-target", "page_current"),
        Input("table-subset-target", "page_size"),
        Input("table-subset-target", "sort_by"),
        Input("table-subset-target", "filter_query"),
        Input("hide-reference-targets", "value"),
    ],
    [
        State("url", "search"),
        State("url", "hash"),
    ],
)
def update_table_subset(
    page_current: int,
    page_size: int,
    sort_by: str,
    filter_term: str,
    hide_reference: str,
    search: str,
    hash_guide: str,
) -> List:
    """The function splits the results according to user-defined filtering
    or sorting criteria.

    The function also updates the visualized results when the user clicks
    the button next/prev page.

    The function loads the CRISPR targets/scores files if available and use
    them to create a pandas DataFrame. The DataFrame column names are changed
    accordingly to those used as IDs of webpage datatable columns.

    If no target is available, the function returns an error message.

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Current page size
    sort_by : str
        Sorting criterion
    filter_term : str
        Filtering criterion
    hide_reference : bool
        Displays only non reference data
    search : str
        Search
    hash_guide : str
        Guide hashing

    Returns
    -------
    List
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(hide_reference, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(hide_reference).__name__}"
        )
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_guide).__name__}")
    # recover job identifier
    job_id = search.split("=")[-1]
    # recover the filtering criterion from drop-down bar
    filter_criterion = read_json(job_id)
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (next(s for s in params.split("\n") if "Ref_comp" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    # value = job_id
    if search is None:
        raise PreventUpdate  # do not do anything
    if filter_term is not None:
        filtering_expressions = filter_term.split(" && ")
    # filtering_expressions.append(['{crRNA} = ' + guide])
    # recover guide, mismatches and bulges
    guide = hash_guide[1 : hash_guide.find("new")]
    mms = hash_guide[-1:]
    bulge_s = hash_guide[-2:-1]
    if "DNA" in hash_guide:
        bulge_t = "DNA"
    elif "RNA" in hash_guide:
        bulge_t = "RNA"
    else:
        bulge_t = "X"
    # choose if hide reference data
    if "hide-ref" in hide_reference or genome_type == "var":
        result = global_store_subset_no_ref(
            job_id, bulge_t, bulge_s, mms, guide, page_current
        )
    else:
        result = global_store_subset(job_id, bulge_t, bulge_s, mms, guide, page_current)
    drop_cols = drop_columns(result, filter_criterion)
    result.drop(drop_cols, axis=1, inplace=True)
    # name of target file filtered with bul-type, mm and bul
    targets_with_mm_bul = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.{bulge_t}.{mms}.{bulge_s}.{guide}.targets.tsv",
    )
    # save df to tsv with filtered data
    result.to_csv(targets_with_mm_bul, sep="\t", na_rep="NA", index=False)
    # change name to zip file
    targets_with_mm_bul_zip = targets_with_mm_bul.replace("tsv", "zip")
    # zip operation, non blocking
    cmd = f"zip -j {targets_with_mm_bul_zip} {targets_with_mm_bul} &"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f"An error occurred while running {cmd}")
    columns_result = [
        {"name": i, "id": i, "hideable": True}
        for col, i in enumerate(result.columns.tolist())
    ]
    data_to_send = result.to_dict("records")
    return [data_to_send, columns_result]


def guidePagev3(job_id, hash):
    guide = hash[: hash.find("new")]
    mms = hash[-1:]
    bulge_s = hash[-2:-1]
    if "DNA" in hash:
        bulge_t = "DNA"
    elif "RNA" in hash:
        bulge_t = "RNA"
    else:
        bulge_t = "X"
    add_header = " - Mismatches " + str(mms)
    if bulge_t != "X":
        add_header += " - " + str(bulge_t) + " " + str(bulge_s)
    value = job_id
    print(job_id)
    print("------")
    if not os.path.isdir(os.path.join(current_working_directory, "Results/", job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    with open(os.path.join(current_working_directory, "Results/", job_id, ".Params.txt")) as p:
        all_params = p.read()
        genome_type_f = (
            next(s for s in all_params.split("\n") if "Genome_selected" in s)
        ).split("\t")[-1]
        ref_comp = (next(s for s in all_params.split("\n") if "Ref_comp" in s)).split(
            "\t"
        )[-1]
        pam = (next(s for s in all_params.split("\n") if "Pam" in s)).split("\t")[-1]

    job_directory = os.path.join(current_working_directory, "Results", job_id)
    genome_type = "ref"
    style_hide_reference = {"display": "none"}
    value_hide_reference = []
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
        style_hide_reference = {}
        value_hide_reference = ["hide-ref"]

    pam_at_start = False
    if str(guide)[0] == "N":
        pam_at_start = True

    final_list = []
    if pam_at_start:
        final_list.append(
            html.H3(
                "Selected Guide: " + str(pam) + str(guide).replace("N", "") + add_header
            )
        )
    else:
        final_list.append(
            html.H3(
                "Selected Guide: " + str(guide).replace("N", "") + str(pam) + add_header
            )
        )
    final_list.append(
        html.P(
            [
                # 'Select a row to view the target IUPAC character scomposition. The rows highlighted in red indicates that the target was found only in the genome with variants.',
                "List of Targets found for the selected guide.",
                dcc.Checklist(
                    options=[{"label": "Hide Reference Targets", "value": "hide-ref"}],
                    id="hide-reference-targets",
                    value=value_hide_reference,
                    style=style_hide_reference,
                ),
                html.Div(
                    [
                        html.P(
                            "Generating download link, Please wait...",
                            id="download-link-sumbyguide",
                        ),
                        dcc.Interval(interval=5 * 1000, id="interval-sumbyguide"),
                    ]
                ),
            ]
        )
    )
    integrated_file_name = glob(
        current_working_directory + "/Results/" + job_id + "/" + "*integrated*"
    )[0]
    integrated_file_name = str(integrated_file_name)
    file_to_grep = os.path.join(job_directory, job_id, ".bestMerge.txt.integrated_results.tsv")
    # file_to_grep_alt = job_directory + job_id + '.altMerge.txt'

    guide_grep_result = (
        job_directory
        + job_id
        + "."
        + bulge_t
        + "."
        + bulge_s
        + "."
        + mms
        + "."
        + guide
        + ".txt"
    )
    # put_header = 'head -1 ' + job_directory + job_id + file_to_grep + ' > ' + guide_grep_result + ' ; '

    final_list.append(
        html.Div(
            f"{job_id}/"
            + "."
            + str(bulge_t)
            + "."
            + str(mms)
            + "."
            + str(bulge_s)
            + "."
            + guide
            + ".targets",
            style={"display": "none"},
            id="div-info-sumbyguide-targets",
        )
    )

    path_db = glob(current_working_directory + "/Results/" + job_id + "/.*.db")[0]
    path_db = str(path_db)
    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    total_private_sample = f"SELECT * FROM final_table LIMIT 1"
    rows = c.execute(total_private_sample)
    header = [description[0] for description in rows.description]
    conn.commit()
    conn.close()

    cols = [{"name": i, "id": i, "hideable": True} for i in header]
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id="table-subset-target",
                # columns=cols,
                style_cell={"textAlign": "left"},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action="custom",
                style_table={"max-height": "600px", "overflowX": "scroll"},
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                css=[
                    {"selector": ".row", "rule": "margin: 0"},
                    {
                        "selector": "td.cell--selected, td.focused",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                    {
                        "selector": "td.cell--selected *, td.focused *",
                        "rule": "background-color: rgba(0, 0, 255,0.15) !important;",
                    },
                ],
                # style_data_conditional=[ ]
            ),
            id="div-result-table",
        )
    )
    final_list.append(html.Br())

    return html.Div(final_list, style={"margin": "1%"})


# @cache.memoize()
def global_store_subset_no_ref(
    job_id: str, bulge_t: str, bulge_s: str, mms: str, guide: str, page: int
) -> pd.DataFrame:
    """Cache targets files to improve visualization performance.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    bulge_t : str
        Bulge type
    bulge_s : str
    mms : str
        Mismatches
    guide : str
        Guide
    page : int
        Current page

    Returns
    -------
    pd.DataFrame
        Results table
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(bulge_t, str):
        raise TypeError(f"Expected {str.__name__}, got {type(bulge_t).__name__}")
    if not isinstance(bulge_s, str):
        raise TypeError(f"Expected {str.__name__}, got {type(bulge_s).__name__}")
    if not isinstance(mms, str):
        raise TypeError(f"Expected {str.__name__}, got {type(mms).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if job_id is None:
        return ""  # do not do anything
    # recover path to db file
    db_path = glob(os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db"))[
        0
    ]  # take the first element
    assert isinstance(db_path, str)
    # initialize db
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # recover the filtering criterion from drop-down bar
    filter_criterion = read_json(job_id)
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    query_cols = get_query_column(filter_criterion)
    # perform query on db
    result = pd.read_sql_query(
        'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}"=\'{}\' AND "{}"={} AND "{}"={} AND "{}"<>\'NA\' LIMIT {} OFFSET {}'.format(
            GUIDE_COLUMN,
            guide,
            query_cols["bul_type"],
            bulge_t,
            query_cols["bul"],
            bulge_s,
            query_cols["mm"],
            mms,
            query_cols["samples"],
            PAGE_SIZE,
            page * PAGE_SIZE,
        ),
        conn,
    )
    return result


# @cache.memoize()
def global_store_subset(
    job_id: str, bulge_t: str, bulge_s: str, mms: str, guide: str, page: int
) -> pd.DataFrame:
    """Cache targets files to improve visualization performance.

    ...

    Parameters
    ----------
    job_id : str
        Unique job identifier
    bulge_t : str
        Bulge type
    bulge_s : str
    mms : str
        Mismatches
    guide : str
        Guide
    page : int
        Current page

    Returns
    -------
    pd.DataFrame
        Res
    """

    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(bulge_t, str):
        raise TypeError(f"Expected {str.__name__}, got {type(bulge_t).__name__}")
    if not isinstance(bulge_s, str):
        raise TypeError(f"Expected {str.__name__}, got {type(bulge_s).__name__}")
    if not isinstance(mms, str):
        raise TypeError(f"Expected {str.__name__}, got {type(mms).__name__}")
    if not isinstance(guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide).__name__}")
    if not isinstance(page, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page).__name__}")
    if job_id is None:
        return ""  # do not do anything
    # recover path to db
    db_path = glob(os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db"))[0]
    assert isinstance(db_path, str)
    # initialize db
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # recover filtering criterion from drop-down bar
    filter_criterion = read_json(job_id)
    if not filter_criterion in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    query_cols = get_query_column(filter_criterion)
    # perform query on db
    result = pd.read_sql_query(
        'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}"=\'{}\' AND "{}"={} AND "{}"={} LIMIT {} OFFSET {}'.format(
            GUIDE_COLUMN,
            guide,
            query_cols["bul_type"],
            bulge_t,
            query_cols["bul"],
            bulge_s,
            query_cols["mm"],
            mms,
            PAGE_SIZE,
            page * PAGE_SIZE,
        ),
        conn,
    )
    # add mismatches and bulges
    targets_with_mm_bul = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        f"{job_id}.{bulge_t}.{mms}.{bulge_s}.{guide}.targets.tsv",
    )
    # store query results in TSV file
    result.to_csv(targets_with_mm_bul, sep="\t", na_rep="NA", index=False)
    return result


# Load barplot of population distribution for selected guide


@app.callback(
    Output("content-collapse-population", "children"),
    [Input("general-profile-table", "selected_cells")],
    [State("general-profile-table", "data"), State("url", "search")],
)
def load_distribution_populations(
    sel_cel: List, all_guides: List[str], job_id: str
) -> List[html.Div]:
    """Load targets distribution by superpopulation and display
    them in the corresponding webpage.

    ...

    Parameters
    ----------
    sel_cel : List
    all_guides : List[str]
        CRISPR guides
    job_id : str
        Unique job identifier

    Returns
    -------
    List[html.Div]
        Webpage with target distribution plots
    """

    if not isinstance(sel_cel, list):
        raise TypeError(f"Expected {list.__name__}, got {type(sel_cel).__name__}")
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate  # do not do anything
    # get the guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    job_id = job_id.split("=")[-1]  # job identifier
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            all_params = handle_params.read()
            mms = (next(s for s in all_params.split("\n") if "Mismatches" in s)).split(
                "\t"
            )[-1]
            mms = int(mms)
            max_bulges = (
                next(s for s in all_params.split("\n") if "Max_bulges" in s)
            ).split("\t")[-1]
            max_bulges = int(max_bulges)
            # total bulge span = DNA + RNA (independent per-type caps, both may
            # co-occur in one alignment); fall back to Max_bulges for old jobs.
            try:
                total_bulges = int(
                    next(s.split("\t")[-1] for s in all_params.split("\n") if s.startswith("DNA\t"))
                ) + int(
                    next(s.split("\t")[-1] for s in all_params.split("\n") if s.startswith("RNA\t"))
                )
            except (StopIteration, ValueError):
                total_bulges = max_bulges
    except OSError as e:
        raise e
    # begin page construction
    distributions = [
        dbc.Row(
            html.P(
                str(
                    "On- and Off-Targets distributions in the Reference and "
                    "Variant Genome. For the Variant Genome, the targets are "
                    "divided into SuperPopulations."
                ),
                style={"margin-left": "0.75rem"},
            )
        )
    ]
    # compute plots
    for i in range(math.ceil((mms + total_bulges + 1) / BARPLOT_LEN)):
        all_images = []
        for mm in range(i * BARPLOT_LEN, (i + 1) * BARPLOT_LEN):
            if mm < (mms + total_bulges + 1):
                try:
                    all_images.append(
                        dbc.Col(
                            [
                                html.A(
                                    html.Img(
                                        src="data:image/png;base64,{}".format(
                                            base64.b64encode(
                                                open(
                                                    os.path.join(
                                                        current_working_directory,
                                                        RESULTS_DIR,
                                                        job_id,
                                                        "imgs",
                                                        "_".join(
                                                            [
                                                                "populations",
                                                                "distribution",
                                                                guide,
                                                                f"{mm}total.png",
                                                            ]
                                                        ),
                                                    ),
                                                    mode="rb",
                                                ).read(),
                                            ).decode()
                                        ),
                                        id=f"distribution-population{mm}",
                                        width="100%",
                                        height="auto",
                                    ),
                                    target="_blank",
                                    href=os.path.join(
                                        RESULTS_DIR,
                                        job_id,
                                        "imgs",
                                        "_".join(
                                            [
                                                "populations",
                                                "distribution",
                                                guide,
                                                f"{mm}total.png",
                                            ]
                                        ),
                                    ),
                                ),
                                html.Div(
                                    html.P(
                                        f"Distribution {mm} Mismatches + Bulges ",
                                        style={"display": "inline-block"},
                                    ),
                                    style={"text-align": "center"},
                                ),
                            ]
                        )
                    )
                except:
                    all_images.append(
                        dbc.Col(
                            [
                                html.Div(
                                    html.P(
                                        f"No Targets found with {mm} Mismatches + Bulges",
                                        style={"display": "inline-block"},
                                    ),
                                    style={"text-align": "center"},
                                ),
                            ],
                            align="center",
                        )
                    )
            else:
                all_images.append(dbc.Col(html.P("")))
        distributions.append(html.Div([dbc.Row(all_images)]))
    return distributions


# Open/close barplot for population distribution
@app.callback(
    Output("collapse-populations", "is_open"),
    [Input("btn-collapse-populations", "n_clicks")],
    [State("collapse-populations", "is_open")],
)
def toggle_collapse_distribution_populations(n, is_open):
    if n:
        return not is_open
    return is_open


# -------------------------------------------------------------------------------
# Custom Ranking tab
#


# trigger guides table construction
@app.callback(
    [
        Output("general-profile-table", "data"),
        Output("general-profile-table", "selected_cells"),
    ],
    [
        Input("general-profile-table", "page_current"),
        Input("general-profile-table", "page_size"),
        Input("general-profile-table", "sort_by"),
        Input("general-profile-table", "filter_query"),
        Input("target_filter_dropdown", "value"),
    ],
    [State("url", "search")],
)
def update_table_general_profile(
    page_current: int,
    page_size: int,
    sort_by: List[str],
    filter_term: str,
    filter_criterion: str,
    search: str,
) -> Tuple[Dict, List]:
    """Construct the custom ranking tab page.
    The tab displays a table summarizing the CRISPRme analysis results
    for each input guide.

    The displayed table columns are filtered according to the filter
    criterion selected by the user through the drop-down bar.

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    sort_by : List[str]
        Sorting criterion
    filter_term : str
        Filter
    filter_criterion : str
        Filter criterion
    search : str
        Search

    Returns
    -------
    Tuple[Dict, List]
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(sort_by, list):
        raise TypeError(f"Expected {list.__name__}, got {type(sort_by).__name__}")
    if not isinstance(filter_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_term).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filter criterion ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    # recover job identifier
    job_id = search.split("=")[-1]
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (next(s for s in params.split("\n") if "Ref_comp" in s)).split(
                "\t"
            )[-1]
            mms = (next(s for s in params.split("\n") if "Mismatches" in s)).split(
                "\t"
            )[-1]
            mms = int(mms)
            max_bulges = (
                next(s for s in params.split("\n") if "Max_bulges" in s)
            ).split("\t")[-1]
            max_bulges = int(max_bulges)
            nuclease = (next(s for s in params.split("\n") if "Nuclease" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    filtering_expressions = filter_term.split(" && ")
    # Get error guides
    error_guides = []
    if os.path.exists(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, "guides_error.txt")
    ):
        try:
            with open(
                os.path.join(current_working_directory, RESULTS_DIR, job_id, "guides_error.txt")
            ) as handle_guide_error:
                for e_g in handle_guide_error:
                    error_guides.append(e_g.strip())
        except OSError as e:
            raise e
    # Get guide from .guide.txt
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, GUIDES_FILE)
        ) as handle_guides:
            guides = handle_guides.read().strip().split("\n")
            guides.sort()
    except OSError as e:
        raise e
    acfd_file = os.path.join(
        current_working_directory,
        RESULTS_DIR,
        job_id,
        "".join([".", job_id, ".acfd_", filter_criterion, ".txt"]),
    )
    if not os.path.isfile(acfd_file):
        raise FileNotFoundError(f"Unable to locate {acfd_file}")
    # load acfd for each guide
    try:
        with open(acfd_file) as handle_acfd:
            all_scores = handle_acfd.read().strip().split("\n")
    except OSError as e:
        raise e
    # Load scores
    if "NO SCORES" not in all_scores:
        all_scores.sort()
        acfd = [
            float(a.split("\t")[1])
            for a in all_scores
            if a.split("\t")[0] not in error_guides
        ]
        doench = [
            a.split("\t")[2] for a in all_scores if a.split("\t")[0] not in error_guides
        ]
        if genome_type == "both":
            doench_enr = [
                a.split("\t")[3]
                for a in all_scores
                if a.split("\t")[0] not in error_guides
            ]
        # acfd = [int(round((100/(100 + x))*100)) for x in acfd]
        acfd = [
            (
                float("{:.3f}".format(x * 100))
                if x < 1 and x >= 0
                else "CFD score not available"
            )
            for x in acfd
        ]
    df = []
    table_to_file = []
    for i, g in enumerate(guides):
        table_to_file.append(g)  # append guide to table
        # append nuclease to table
        table_to_file.append(f"Nuclease: {nuclease}")
        data_general_count = pd.read_csv(
            os.path.join(
                current_working_directory,
                RESULTS_DIR,
                job_id,
                f".{job_id}.general_target_count.{g}_{filter_criterion}.txt",
            ),
            sep="\t",
            na_filter=False,
        )
        data_guides = {}
        data_guides["Guide"] = g
        data_guides["Nuclease"] = nuclease
        data_general_count_copy = data_general_count.copy()

        n_rows = len(data_general_count_copy)
        half = n_rows // 2
        origin_concat = ["REF"] * half + ["VAR"] * (n_rows - half)
        count_bulges_concat = list(range(half)) + list(range(n_rows - half))

        data_general_count_copy.insert(0, "Genome", origin_concat, True)
        data_general_count_copy.insert(1, "Bulges", count_bulges_concat, True)
        if "NO SCORES" not in all_scores:
            data_guides["CFD"] = acfd[i]
            table_to_file.append(f"Score: {acfd[i]}")  # append CFD to table
            # append CFD to table
            table_to_file.append(f"Filter_criterion: {filter_criterion}")
            table_to_file.append("\t\t\t\tMismatches")
            table_to_file.append(data_general_count_copy.to_string(index=False))
            if genome_type == "both":
                data_guides["Doench 2016"] = doench[i]
            else:
                data_guides["Doench 2016"] = doench[i]
        # Bulge ROWS span 0..(bDNA+bRNA); the count file (built by process_summaries
        # with bDNA+bRNA+1 rows per origin) is the source of truth, so derive the row
        # labels/placement from len(data_general_count), NOT Max_bulges (= max(bDNA,bRNA)).
        # The old code hardcoded max_bulges in {0,1,2} and silently dropped the 3-/4-bulge
        # rows of e.g. a 2/2 search from this summary table.
        nb = len(data_general_count)
        per_origin = (nb // 2) if genome_type == "both" else nb
        if genome_type == "both":
            labels = [str(j) for j in range(per_origin)] * 2
            labels.insert(len(labels) // 2, "")
            data_guides["# Bulges"] = "\n".join(labels)
        else:
            data_guides["# Bulges"] = "\n".join(str(j) for j in range(per_origin))
        data_guides["Total"] = []
        # Per-row Totals with REFERENCE/VARIANT labels centered in each origin block,
        # aligned to the # Bulges / MM cells (which carry a blank separator between the
        # two halves). ref_mid/var_mid reproduce the old hardcoded placement exactly for
        # per_origin in {1,2,3} (max_bulges 0/1/2) and generalize to any bulge count.
        totals = [str(sum(data_general_count.iloc[j, :])) for j in range(nb)]
        if genome_type == "both":
            ref_mid = per_origin // 2
            var_mid = per_origin + (per_origin // 2)
            for j, t in enumerate(totals):
                if j == ref_mid:
                    data_guides["Total"].append("REFERENCE\t" + t)
                elif j == var_mid:
                    data_guides["Total"].append("VARIANT\t\t" + t)
                else:
                    data_guides["Total"].append("\t" + t)
            data_guides["Total"].insert(len(data_guides["Total"]) // 2, "")
        else:
            mid = nb // 2
            for j, t in enumerate(totals):
                data_guides["Total"].append(
                    ("REFERENCE\t" + t) if j == mid else ("\t" + t)
                )
        for mm_i in range(mms + 1):
            col = list(data_general_count.iloc[:, mm_i].values.astype(str))
            if genome_type == "both":
                col.insert(len(col) // 2, "")  # blank separator between REF/VAR halves
            # NO truncation: show every bulge row the data contains (was tmp[:max_bulges+1],
            # which dropped rows with total bulges > max(bDNA,bRNA)).
            data_guides[str(mm_i) + "MM"] = "\n".join(col)
        data_guides["Total"] = "\n".join(data_guides["Total"])
        df.append(data_guides)
    dff = pd.DataFrame(df)  # create data table
    table_to_file_save_dest = os.path.join(
        current_working_directory, RESULTS_DIR, job_id, f"{job_id}.general_table.txt"
    )
    try:
        outfile = open(table_to_file_save_dest, "w")
        for elem in table_to_file:
            outfile.write(elem + "\n")
    except OSError as e:
        raise e
    finally:
        outfile.close()
    # zip integrated results
    integrated_fname = glob(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, "*integrated*")
    )[0]
    assert isinstance(integrated_fname, str)
    # integrated_file = integrated_fname
    # zip integrated file
    integrated_to_zip = integrated_fname.replace("tsv", "zip")
    if not os.path.exists(integrated_to_zip):
        cmd = f"zip -j {integrated_to_zip} {integrated_fname} &"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
    # zip alt_merge results
    alt_merge_fname = glob(
        os.path.join(
            current_working_directory, RESULTS_DIR, job_id, "*all_results_with_alternative_alignments*"
        )
    )[0]
    assert isinstance(alt_merge_fname, str)
    # integrated_file = alt_merge_fname
    # zip integrated file
    alt_merge_to_zip = alt_merge_fname.replace("tsv", "zip")
    if not os.path.exists(alt_merge_to_zip):
        cmd = f"zip -j {alt_merge_to_zip} {alt_merge_fname} &"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
    # score checking
    if "NO SCORES" not in all_scores:
        try:
            dff = dff.sort_values(["CFD", "Doench 2016"], ascending=[False, False])
        except:  # for BOTH
            dff = dff.sort_values(["CFD", "Enriched"], ascending=[False, False])
    else:
        try:
            dff = dff.sort_values("On-Targets Reference", ascending=True)
        except:
            dff = dff.sort_values("On-Targets Enriched", ascending=True)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if operator in PANDAS_OPERATORS:
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == "contains":
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]
    if bool(sort_by):
        dff = dff.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )
    # Return every guide row (page_action="none" on this table): the Result Summary
    # must fit on the page in full, without a pager or vertical scroll. Custom
    # filter/sort above still apply; there just is no page slice.
    data_to_send = dff.to_dict("records")
    return data_to_send, [{"row": 0, "column": 0}]


# Update color on selected row
@app.callback(
    Output("general-profile-table", "style_data_conditional"),
    [Input("general-profile-table", "selected_cells")],
    [State("general-profile-table", "data")],
)
def color_selected_row(sel_cel: List, all_guides: List) -> List:
    """Color the selected row of the data table.

    ...

    Parameters
    ----------
    sel_cel : List
        Selected row
    all_guides : List
        Guides list

    Returns
    -------
    List
    """

    # check if the table has to be updated or not
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate  # do not do anything
    # recover the guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # color the selected row
    return [
        {
            "if": {
                "filter_query": '{Guide} eq "' + guide + '"',
            },
            "background-color": "rgba(0, 0, 255,0.15)",  # rgb(255, 102, 102)
        },
        {"if": {"column_id": "Genome"}, "font-weight": "bold", "textAlign": "center"},
    ]


# ------------------------------------------------------------------------------
# Query genomic region tab
#


# trigger filtering table by genomic coordinates
@app.callback(
    [
        Output("div-table-position", "children"),
        Output("div-current-page-table-position", "children"),
    ],
    [Input("div-position-filter-query", "children")],
    [
        State("button-filter-position", "n_clicks_timestamp"),
        State("target_filter_dropdown", "value"),
        State("url", "search"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("div-current-page-table-position", "children"),
    ],
)
def filter_position_table(
    filter_q: List[str],
    n: int,
    filter_criterion: str,
    search: str,
    sel_cel: List[int],
    all_guides: List[int],
    current_page: str,
) -> Tuple[List[html.P], str]:
    """Filter result table by genomic region. The table is filtered in order to
    display only those targets falling within the genomic interval defined
    by the user.

    The results can be furtherly filtered by scoring criterion. The available
    criteria are CFD score, CRISTA score and the number of mismatches and bulges.

    ...

    Parameters
    ----------
    filter_q : List[str]
        Filtering query (coordinates, filtering criterion)
    n : int
        Input click listened
    filter_criterion : str
        result table filtering criterion
    search : str
        Search ID
    sel_cel : List[int]
    all_guides : List[int]
        List of the guides
    current_page : str
        Current table page number

    Returns
    -------
    List[html.P]
        HTML result page
    str
        Page numeration

    """

    if n is not None and not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if not filter_criterion in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if sel_cel is None:
        raise PreventUpdate
    if n is None:
        raise PreventUpdate
    # recover filter query fields;
    # if there are NULL fields -> prevent table update
    # query structure: chrom,start,stop
    if isinstance(filter_q, str):  # simple regular query
        filter_q = filter_q.split(",")
        assert isinstance(filter_q, list)  # it should be list of fields
    elif isinstance(filter_q, list):  # updated by callback
        assert len(filter_q) == 2  # we should have just two elements
        filter_criterion = filter_q[1]  # recover table filtering criterion
        filter_q = filter_q[0].split(",")  # query genomic coordinates
    assert filter_criterion in FILTERING_CRITERIA
    chrom = filter_q[0]
    if chrom == "None":
        raise PreventUpdate  # invalid chromosome
    start = filter_q[1]
    if start == "None":
        raise PreventUpdate  # invalid start
    end = filter_q[2]
    if end == "None":
        raise PreventUpdate  # invalid stop
    current_page = int(current_page.split("/")[0])
    job_id = search.split("=")[-1]
    # recover the guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # recover db file
    db_path = glob(os.path.join(current_working_directory, RESULTS_DIR, job_id, ".*.db"))[0]
    assert isinstance(db_path, str)
    assert os.path.isfile(db_path)
    # connect db with sqlite3
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = 'SELECT * FROM final_table WHERE "{}"=\'{}\' AND "{}">={} AND "{}"<={} AND "{}"=\'{}\''
    # recover filtering criterion
    filter_criterion = read_json(job_id)
    query_cols = get_query_column(filter_criterion)  # recover query columns
    # query the db
    result = pd.read_sql_query(
        query.format(
            GUIDE_COLUMN,
            guide,
            query_cols["start"],
            start,
            query_cols["start"],
            end,
            CHR_COLUMN,
            chrom,
        ),
        conn,
    )
    conn.commit()
    conn.close()  # close db connection
    assert isinstance(result, pd.DataFrame)
    if result.empty:  # no guides found ?
        df_check = False
    else:  # check table fit to page
        df_check = True
    # filter diplayed column using filtering criterion
    drop_cols = drop_columns(result, filter_criterion)
    result.drop(drop_cols, inplace=True, axis=1)  # remove columns from table
    # check table characteristics to fit it into html page
    if df_check:
        out_1 = [
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="table-position",
                export_format="xlsx",
                columns=[
                    {"name": i, "id": i, "hideable": True}
                    for count, i in enumerate(result.columns)
                ],
                data=result.to_dict("records"),
                style_cell={"textAlign": "left"},
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                style_table={
                    "overflowX": "scroll",
                },
                page_size=PAGE_SIZE,
            )
        ]
    else:
        out_1 = [html.P("No results found with this genomic coordinates")]
    return out_1, "/".join([str(1), str(1)])


# update filters to filter results by position
@app.callback(
    Output("div-position-filter-query", "children"),
    [Input("button-filter-position", "n_clicks")],
    [
        State("target_filter_dropdown", "value"),
        State("dropdown-chr-table-position", "value"),
        State("input-position-start", "value"),
        State("input-position-end", "value"),
    ],
)
def update_position_filter(
    n: int, filter_criterion: str, chrom: str, pos_start: str, pos_end: str
) -> Tuple[str, str, int]:
    """Callback to update the result table filtered by genomic location.

    ...

    Parameters
    ----------
    n : int
        Number of clicks listened
    filter_criterion : str
        Filtering citerion to apply to the table
    chrom : str
        Chromosome
    pos_start : str
        Start position
    pos_end : str
        Stop position

    Returns
    -------
    Tuple[str, str, int]
        New genomic locations and potential filtering criterion
    """

    # Dash fires this callback once on tab render with n=None (no click) and the
    # State values still None; guard *before* the type checks so an initial
    # render is a no-op instead of a 500 (TypeError: Expected str, got NoneType).
    if n is None:  # no click -> no page update
        raise PreventUpdate
    if not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    # The chromosome/position State values are None until the user picks them.
    # Normalize empties to the sentinel "None" so a Filter click with nothing
    # selected is a graceful no-op downstream (chrom == "None" -> PreventUpdate)
    # rather than a crash.
    chrom = "None" if chrom in (None, "") else chrom
    pos_start = "None" if pos_start in (None, "") else pos_start
    pos_end = "None" if pos_end in (None, "") else pos_end
    coords = ",".join([chrom, pos_start, pos_end])
    return coords, filter_criterion


# ------------------------------------------------------------------------------
# Summary by Sample tab
#


# View next/prev page on sample table
@app.callback(
    [
        Output("div-table-samples", "children"),
        Output("div-current-page-table-samples", "children"),
    ],
    [
        Input("prev-page-sample", "n_clicks_timestamp"),
        Input("next-page-sample", "n_clicks_timestamp"),
        Input("div-sample-filter-query", "children"),
        Input("target_filter_dropdown", "value"),
    ],
    [
        State("button-filter-population-sample", "n_clicks_timestamp"),
        State("url", "search"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("div-current-page-table-samples", "children"),
    ],
)
def filter_sample_table(
    n_prev: int,
    n_next: int,
    filter_q: str,
    filter_criterion: str,
    n: int,
    search: str,
    sel_cel: List,
    all_guides: List,
    current_page: str,
) -> Tuple[html.Table, str]:
    """Filter summary by sample table according to the filtering crietrion
    selected by the user.

    The adopted filtering criterion is selected by the user through the
    drop-down bar available for all webpage result tabs.

    ...

    Parameters
    ----------
    n_prev : int
        Previous pages number
    n_next :
        Next pages number
    filter_q : str
        Filter query
    filter_criterion : str
        Filtering criterion
    n : int
        Clicks
    search : str
        Search
    sel_cel : List
        Selected table rows
    all_guides : List
        Guides list
    current_page : str
        Current webpage

    Returns
    -------
    Tuple[html.Table, str]
        Updated samples table
    """

    if n_prev is not None and not isinstance(n_prev, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n_prev).__name__}")
    if n_next is not None and not isinstance(n_next, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n_next).__name__}")
    if n is not None and not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(current_page, str):
        raise TypeError(f"Expected {str.__name__}, got {type(current_page).__name__}")
    if sel_cel is None:
        raise PreventUpdate  # do not do anything
    if n_prev is None and n_next is None and n is None:
        raise PreventUpdate  # do not do anything
    if n_prev is None:
        n_prev = 0
    if n_next is None:
        n_next = 0
    if n is None:
        n = 0
    # get superpopulation name
    sup_pop = filter_q.split(",")[0]
    # get population name
    pop = filter_q.split(",")[1]
    # get sample
    sample = str(filter_q.split(",")[2])
    if sup_pop == "None":
        sup_pop = None
    if pop == "None":
        pop = None
    if sample == "None" or sample == "NONE":
        sample = None
    current_page = int(current_page.split("/")[0])
    btn_sample_section = [n, n_prev, n_next]
    # get job identifier
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    population_1000gp = associateSample.loadSampleAssociation(
        os.path.join(job_directory, SAMPLES_ID_FILE)
    )[2]
    # read CRISPRme run parameters
    try:
        with open(
            os.path.join(current_working_directory, RESULTS_DIR, job_id, PARAMS_FILE)
        ) as handle_params:
            params = handle_params.read()
            genome_type_f = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            ref_comp = (next(s for s in params.split("\n") if "Ref_comp" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    genome_type = "ref"
    if "+" in genome_type_f:
        genome_type = "var"
    if "True" in ref_comp:
        genome_type = "both"
    # recover the guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    if genome_type == "both":
        col_names_sample = [
            "Sample",
            "Sex",
            "Population",
            "Super Population",  # 'Targets in Reference',
            "Targets in Sample",
            "Targets in Population",
            "Targets in Super Population",
            "PAM Creation",
        ]  # , 'Class']
    else:
        col_names_sample = [
            "Sample",
            "Sex",
            "Population",
            "Super Population",  # 'Targets in Reference',
            "Targets in Sample",
            "Targets in Population",
            "Targets in Super Population",
            "PAM Creation",
        ]  # , 'Class']
    # Last button pressed is filtering, return the first page of the
    # filtered table
    if max(btn_sample_section) == n:
        df = pd.read_csv(
            os.path.join(
                job_directory,
                f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt",
            ),
            sep="\t",
            names=col_names_sample,
            skiprows=2,
            na_filter=False,
        )
        df = df.sort_values("Targets in Sample", ascending=False)
        more_info_col = ["Show Targets" for _ in range(df.shape[0])]
        df[""] = more_info_col
        if (
            (sup_pop is not None and sup_pop != "")
            or (pop is not None and pop != "")
            or (sample is not None and sample != "")
        ):
            if sample is None or sample == "":
                if pop is None or pop == "":
                    df.drop(
                        df[
                            (~(df["Population"].isin(population_1000gp[sup_pop])))
                        ].index,
                        inplace=True,
                    )
                else:
                    df.drop(df[(df["Sample"] != sample)].index, inplace=True)
            else:
                df.drop(df[(df["Sample"] != sample)].index, inplace=True)
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        return (
            generate_table_samples(df, "table-samples", 1, guide, job_id),
            f"1/{max_page}",
        )
    else:
        if max(btn_sample_section) == n_next:  # go to next page
            current_page += 1
            df = pd.read_csv(
                os.path.join(
                    job_directory,
                    f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt",
                ),
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            if genome_type == "both":
                df = df.sort_values("Targets in Sample", ascending=False)
            else:
                df = df.sort_values("Targets in Reference", ascending=False)
            more_info_col = ["Show Targets" for _ in range(df.shape[0])]
            df[""] = more_info_col
            # Active filter
            if pop or sup_pop or sample:
                if sample is None or sample == "":
                    if pop is None or pop == "":
                        df.drop(
                            df[
                                (~(df["Population"].isin(population_1000gp[sup_pop])))
                            ].index,
                            inplace=True,
                        )
                    else:
                        df.drop(df[(df["Population"] != pop)].index, inplace=True)
                else:
                    df.drop(df[(df["Sample"] != sample)].index, inplace=True)
            if ((current_page - 1) * 10) > len(df):
                current_page = current_page - 1
                if current_page < 1:
                    current_page = 1
        else:  # go to previous page
            current_page -= 1
            if current_page < 1:
                current_page = 1
            df = pd.read_csv(
                os.path.join(
                    job_directory,
                    f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt",
                ),
                sep="\t",
                names=col_names_sample,
                skiprows=2,
                na_filter=False,
            )
            if genome_type == "both":
                df = df.sort_values("Targets in Sample", ascending=False)
            else:
                df = df.sort_values("Targets in Sample", ascending=False)
            more_info_col = ["Show Targets" for _ in range(df.shape[0])]
            df[""] = more_info_col
            if pop or sup_pop or sample:
                if sample is None or sample == "":
                    if pop is None or pop == "":
                        df.drop(
                            df[
                                (~(df["Population"].isin(population_1000gp[sup_pop])))
                            ].index,
                            inplace=True,
                        )
                    else:
                        df.drop(df[(df["Population"] != pop)].index, inplace=True)
                else:
                    df.drop(df[(df["Sample"] != sample)].index, inplace=True)
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        return (
            generate_table_samples(df, "table-samples", current_page, guide, job_id),
            f"{current_page}/{max_page}",
        )


# Callback to update the hidden div filter
@app.callback(
    Output("div-sample-filter-query", "children"),
    [Input("button-filter-population-sample", "n_clicks")],
    [
        State("dropdown-superpopulation-sample", "value"),
        State("dropdown-population-sample", "value"),
        State("input-sample", "value"),
    ],
)
def update_sample_filter(
    n: int, superpopulation: str, population: str, sample: str
) -> str:
    """Update the filter for the samples table.

    ...

    Parameters
    ----------
    n : int
        Clicks
    superpopulation : str
        Superpopulation code
    population : str
        Population code
    sample : str
        Sample identifier

    Returns
    -------
    str
        New filter query
    """

    if n is not None and not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if superpopulation is not None and not isinstance(superpopulation, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(superpopulation).__name__}"
        )
    if population is not None and not isinstance(population, str):
        raise TypeError(f"Expected {str.__name__}, got {type(population).__name__}")
    if sample is not None and not isinstance(sample, str):
        raise TypeError(f"Expected {str.__name__}, got {type(sample).__name__}")
    if n is None:
        raise PreventUpdate
    # prevent page updates when at least one filter element is none
    if any([field is None for field in [superpopulation, population, sample]]):
        raise PreventUpdate
    filter_new = ",".join(
        [superpopulation, population, sample.replace(" ", "").upper()]
    )
    return filter_new


# Callback to update the sample based on population selected
@app.callback(
    [Output("dropdown-sample", "options"), Output("dropdown-sample", "value")],
    [Input("dropdown-population-sample", "value")],
    [State("url", "search")],
)
def update_sample_drop(pop: str, search: str) -> Tuple[List, None]:
    """Update Summary by Sample data table accordingly to the selected
    population.

    ...

    Parameters
    ----------
    pop : str
        Population
    search: str
        Search

    Returns
    -------
    Tuple[List, None]
    """

    if pop is not None and not isinstance(pop, str):
        raise TypeError(f"Expected {str.__name__}, got {type(pop).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if pop is None or pop == "":
        return [], None  # no update required
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    pop_dict = associateSample.loadSampleAssociation(
        os.path.join(job_directory, SAMPLES_ID_FILE)
    )[3]
    return [{"label": sample, "value": sample} for sample in pop_dict[pop]], None


# Callback to update the population tab based on superpopulation selected
@app.callback(
    [
        Output("dropdown-population-sample", "options"),
        Output("dropdown-population-sample", "value"),
    ],
    [Input("dropdown-superpopulation-sample", "value")],
    [State("url", "search")],
)
def update_population_drop(superpop: str, search: str) -> Tuple[Dict, None]:
    """Update Summary by Sample data table accordingly to the selected
    superpopulation.

    ...

    Parameters
    ----------
    pop : str
        Population
    search: str
        Search

    Returns
    -------
    Tuple[List, None]
    """

    if superpop is not None and not isinstance(superpop, str):
        raise TypeError(f"Expected {str.__name__}, got {type(superpop).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if superpop is None or superpop == "":
        raise PreventUpdate  # no update required
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    population_1000gp = associateSample.loadSampleAssociation(
        os.path.join(job_directory, SAMPLES_ID_FILE)
    )[2]
    return [{"label": i, "value": i} for i in population_1000gp[superpop]], None


def check_existance_sample(job_directory: str, job_id: str, sample: str) -> bool:
    """Check if the selected sample exists in the dataset.

    ...

    Parameters
    ----------
    job_directory : str
        Path to job results
    job_id : str
        Unique job identifier
    sample : str
        Sample

    Returns
    -------
    bool
    """

    if not isinstance(job_directory, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_directory).__name__}")
    if not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if not isinstance(sample, str):
        raise TypeError(f"Expected {str.__name__}, got {type(sample).__name__}")
    dataset = pd.read_csv(
        os.path.join(job_directory, job_id, SAMPLES_ID_FILE), sep="\t", na_filter=False
    )
    samples = dataset.iloc[:, 0].tolist()
    if sample in samples:
        return True
    return False


# -------------------------------------------------------------------------------
# Graphical Reports tab
#


# Select figures on mms value, sample value
@app.callback(
    [
        Output("div-radar-chart-encode_gencode", "children"),
        Output("div-population-barplot", "children"),
    ],
    [
        Input("mm-dropdown", "value"),
        Input("general-profile-table", "selected_cells"),
        Input("target_filter_dropdown", "value"),
    ],
    [State("url", "search"), State("general-profile-table", "data")],
)
def update_images_tabs(
    mm: int, sel_cel: List, filter_criterion: str, search: str, all_guides: List
) -> Tuple[List, List]:
    """Compute the plots displayed when watching at the Graphical Reports
    tab in the main CRISPRme results webpage.

    The plots are computed at execution time.

    ...

    Parameters
    ----------
    mm : str
        Mismatches
    sel_cel : List
        Selected table cells
    filter_criterion : str
        Filter criterion selected by the user via the global drop-down bar
    search : str
        Search
    all_guides : List
        All CRISPR guides

    Returns
    -------
    Tuple[List, List]
        The ENCODE+GENCODE radar chart and the population barplot containers
    """

    if not isinstance(mm, str):
        raise TypeError(f"Expected {str.__name__}, got {type(mm).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    # Dash fires this on tab render before a row is selected; guard before indexing
    # sel_cel (mirrors update_content_tab / generate_sample_card) so the initial render
    # is a no-op instead of a crash on sel_cel[0] / all_guides[...].
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    bulge = 0
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # For a non-SpCas9 nuclease, CFD/CRISTA are not computed, so the barplot + radar
    # images are only produced for the "fewest mm+bulges" criterion. Remap here (as
    # update_content_tab does) so this tab reads the files that actually exist instead
    # of rendering "No result found".
    try:
        with open(os.path.join(job_directory, PARAMS_FILE)) as _pf:
            _params = _pf.read()
        _nuclease = (
            next(s for s in _params.split("\n") if "Nuclease" in s)
        ).split("\t")[-1]
        if _nuclease != CAS9:
            filter_criterion = FILTERING_CRITERIA[0]  # fewest mm + bulges
    except (OSError, StopIteration):
        pass
    # define plot containers
    # radar_chart_images = list()
    radar_chart_encode_gencode = []
    # radar_chart_gencode = list()
    population_barplots = []
    # begin graphical reports page construction
    try:
        # population barplot
        population_barplots.extend(
            [
                html.A(
                    html.Img(
                        src="data:image/png;base64,{}".format(
                            base64.b64encode(
                                open(
                                    os.path.join(
                                        current_working_directory,
                                        RESULTS_DIR,
                                        job_id,
                                        IMGS_DIR,
                                        str(
                                            f"populations_distribution_{guide}"
                                            f"_{int(mm) + int(bulge)}total_{filter_criterion}_new.png"
                                        ),
                                    ),
                                    mode="rb",
                                ).read()
                            ).decode()
                        ),
                        id=f"distribution-population{int(mm) + int(bulge)}",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                    href=os.path.join(
                        f"/{RESULTS_DIR}",
                        job_id,
                        IMGS_DIR,
                        str(
                            f"populations_distribution_{guide}_"
                            f"{int(mm) + int(bulge)}total_{filter_criterion}_new.png"
                        ),
                    ),
                ),
            ]
        )
    except:
        population_barplots = [
            html.Div(
                html.H2("No result found for this combination of mismatches and bulges")
            )
        ]
    # radar chart
    radar_img_encode_gencode = os.path.join(
        f"{IMGS_DIR}",
        str(
            f"summary_single_guide_{guide}_{mm}."
            f"{bulge}_TOTAL_{filter_criterion}.ENCODE+GENCODE.png"
        ),
    )
    # TODO: do not call python script, rather define functions
    cmd = f"python {app_directory}/PostProcess/generate_img_radar_chart.py {guide} {job_directory}/.guide_dict_{guide}_{filter_criterion}.json {job_directory}/.motif_dict_{guide}_{filter_criterion}.json {mm} {bulge} TOTAL_{filter_criterion} {job_directory}/imgs/"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f'An error occurred while running "{cmd}"')
    img_found = False  # look for radar chart image
    try:
        radar_src_encode_gencode = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    os.path.join(
                        current_working_directory,
                        RESULTS_DIR,
                        job_id,
                        radar_img_encode_gencode,
                    ),
                    mode="rb",
                ).read()
            ).decode()
        )
        img_found = True
    except:
        pass  # ignore
    try:
        radar_href_encode_gencode = (
            "/Results/" + job_id + "/" + radar_img_encode_gencode
        )
    except:
        # assign the variable actually used below (href=radar_href_encode_gencode);
        # the old code set an unused 'radar_href' -> NameError if this ever failed
        radar_href_encode_gencode = ""
    if img_found:
        radar_chart_encode_gencode.append(
            html.A(
                html.Img(
                    src=radar_src_encode_gencode,
                    id="radar-img-guide",
                    width="100%",
                    height="auto",
                ),
                target="_blank",
                href=radar_href_encode_gencode,
            )
        )
    if len(radar_chart_encode_gencode) == 0:  # no radar chart
        radar_chart_encode_gencode.append(
            html.H2("No result found for this combination of mismatches and bulges")
        )
    return (
        radar_chart_encode_gencode,
        population_barplots,
    )


@app.callback(
    [
        Output("download-link-personal-card", "children"),
        Output("download-link-personal-card", "hidden"),
        Output("div-personal-plot", "children"),
        Output("div-private-plot", "children"),
        Output("div-table-sample-card", "children"),
        Output("div-top-target-sample-card", "children"),
    ],
    [Input("button-sample-card", "n_clicks")],
    [
        State("target_filter_dropdown", "value"),
        State("dropdown-sample-card", "value"),
        State("general-profile-table", "selected_cells"),
        State("general-profile-table", "data"),
        State("url", "search"),
    ],
)
def generate_sample_card(
    n: int,
    filter_criterion: str,
    sample: str,
    sel_cel: List,
    all_guides: List,
    search: str,
) -> List:
    """Generate the sample risk card for each CRISPR guide analyzed.

    The webpage plots the top 1000 personal and private targets, showing their
    allelic frequency. The results can be filtered by the user selecting one
    criterion from the general drop-down bar.

    ...

    Parameters
    ----------
    n : int
        Clicks
    filter_criterion : str
        Filtering criterion
    sample : str
        Sample ID
    sel_cel : List
        Selected cells
    all_guides : List
        All CRISPR guides
    search : str
        Search

    Returns
    -------
    List
        Sample card webpage
    """

    # Dash fires this callback once on tab render with n=None (no click) and a
    # None sample; guard *before* the type checks so the initial render is a
    # no-op instead of a 500 (TypeError: Expected str, got NoneType).
    if n is None:
        raise PreventUpdate  # do not do anything
    if not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion ({filter_criterion})")
    # no sample chosen yet -> don't try to build a card (graceful no-op instead
    # of a crash when the user clicks Generate without selecting a sample)
    if sample is None or sample == "":
        raise PreventUpdate
    if not isinstance(sample, str):
        raise TypeError(f"Expected {str.__name__}, got {type(sample).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    # recover guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # recover job id
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    # directory holding the personal/private target plots for this job
    imgsdir = os.path.join(job_directory, "imgs")
    # read summary by sample data
    samples_summary = pd.read_csv(
        os.path.join(
            job_directory, f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt"
        ),
        sep="\t",
        skiprows=2,  # skip first two rows (guide + header)
        index_col=0,  # use sample ids as index
        header=None,  # header has been skipped
        na_filter=False,  # keep nan values
    )
    # select the number of personal targets and pam creation events for sample
    # these data are in the 4th and 7th column of the summary by sample files,
    # respectively
    targets_sample_num = samples_summary.loc[sample, 4]
    pam_creation_sample_num = samples_summary.loc[sample, 7]
    # output filenames to store personal and private targets for  sample
    targets_personal_fname = os.path.join(
        job_directory, f"{job_id}.{sample}.{guide}.personal_targets.tsv"
    )
    targets_private_fname = os.path.join(
        job_directory, f"{job_id}.{sample}.{guide}.private_targets.tsv"
    )
    # template for the cached plot filenames -> .format(criterion, guide, sample, kind)
    plotfname = os.path.join(
        imgsdir,
        "CRISPRme_{0}_top_1000_log_for_main_text_{1}.{2}.{3}.png",
    )
    if os.path.isfile(targets_private_fname) and (
        all(
            os.path.isfile(plotfname.format(c, guide, sample, "personal"))
            for c in FILTERING_CRITERIA
        )
        and all(
            os.path.isfile(plotfname.format(c, guide, sample, "private"))
            for c in FILTERING_CRITERIA
        )
    ):  # cached private targets for sample
        targets_private = pd.read_csv(targets_private_fname, sep="\t")
        # sort private targets according the filtering criterion
        order = filter_criterion == FILTERING_CRITERIA[0]
        criterion_cname = get_query_column(filter_criterion)["sort"]
        targets_private = targets_private.sort_values(
            [criterion_cname], ascending=order
        )
        # also load the full personal-candidate set (written alongside the private file
        # on first query) so the card can show ALL candidate sites, not just the private
        # subset; fall back to the private set if an older cache lacks the personal file
        if os.path.isfile(targets_personal_fname):
            targets_personal = pd.read_csv(targets_personal_fname, sep="\t").sort_values(
                [criterion_cname], ascending=order
            )
        else:
            targets_personal = targets_private
    else:  # first query for this sample
        # sql database to be used for queries
        db_path = os.path.join(job_directory, f".{job_id}.db")
        if not os.path.isfile(db_path):
            raise FileNotFoundError(f"Cannot locate database file {db_path}")
        try:
            conn = sqlite3.connect(db_path)  # connect to database
            c = conn.cursor()  # sql database cursor to execute queries
            query_cols = get_query_column(filter_criterion)  # get query columns
            # recover colname for sample column and filtering/sorting criterion
            samples_cname, criterion_cname = query_cols["samples"], query_cols["sort"]
            # define and perform sql database to retrieve sample targets
            sqlquery = (
                f"SELECT * FROM final_table WHERE \"{GUIDE_COLUMN}\"='{guide}' "
                f"AND \"{samples_cname}\" LIKE '%{sample}%'"
            )
            targets_personal = pd.read_sql_query(sqlquery, conn)  # perform query
        except sqlite3.Error as e:
            raise sqlite3.Error(f"Database error ({db_path})") from e
        except KeyError as e:
            raise KeyError(
                "Key error - possibly missing column in query columns"
            ) from e
        except Exception as e:
            # sourcery skip: raise-specific-error
            raise Exception(
                f"Unexpected error while querying {db_path} for sample {sample}"
            ) from e
        finally:
            if "conn" in locals():  # close connection to database
                conn.commit()
                conn.close()
        if targets_personal.shape[0] != targets_sample_num:
            raise ValueError(
                f"Mismatching sample targets number (expected: {targets_sample_num}, got: {targets_personal.shape[0]})"
            )
        # sort personal targets - sort in ascending order if fewest is the criterion
        # top targets have less mm+bulges
        order = filter_criterion == FILTERING_CRITERIA[0]
        targets_personal = targets_personal.sort_values(
            [criterion_cname], ascending=order
        )
        # retrieve sample private targets from personal targets
        targets_private = targets_personal[targets_personal[samples_cname] == sample]
        # plot images to be displayed in the personal risk card tab
        try:
            targets_personal.to_csv(targets_personal_fname, sep="\t", index=False)
            targets_private.to_csv(targets_private_fname, sep="\t", index=False)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"File path not found {e.filename}") from e
        except PermissionError as e:
            raise PermissionError(f"Permission denied: {e.filename}") from e
        except OSError as e:
            if e.errno == errno.ENOSPC:
                raise OSError("No space left on device") from e
            else:
                raise OSError(f"Cannot write {e.filename}") from e
        except Exception as e:
            # sourcery skip: raise-specific-error
            raise Exception(
                "An unexpected error occurred while saving personal and integrated targets"
            ) from e
        # compute sample's personal and private targets with lolliplots
        crisprme_plots_personal = (
            f"python {app_directory}/PostProcess/CRISPRme_plots_personal.py"
        )
        crisprme_plots_personal_cmd = (
            f"{crisprme_plots_personal} {targets_personal_fname} {imgsdir}/"
        )
        code = subprocess.call(
            f"{crisprme_plots_personal_cmd} {guide}.{sample}.personal", shell=True
        )
        if code != 0:
            raise subprocess.SubprocessError(
                f"Failed personal targets plot generation on sample {sample}"
            )
        code = subprocess.call(
            f"{crisprme_plots_personal_cmd} {guide}.{sample}.private", shell=True
        )
        if code != 0:
            raise subprocess.SubprocessError(
                f"Failed private targets plot generation on sample {sample}"
            )
    # compress private targets file for download
    targets_private_zip = os.path.join(
        job_directory,
        os.path.basename(f"{os.path.splitext(targets_private_fname)[0]}.zip"),
    )
    if not os.path.isfile(targets_private_zip):
        code = subprocess.call(
            f"zip -j {targets_private_zip} {targets_private_fname}", shell=True
        )
        if code != 0:
            raise subprocess.SubprocessError(
                f"Failed to compress {targets_private_fname}"
            )
    # compute targets stats for sample -> displayed on top of the page
    sample_stats = pd.DataFrame(
        {
            "Personal": [targets_sample_num],
            "PAM creation": [pam_creation_sample_num],
            "Private": [targets_private.shape[0]],
        }
    ).astype(str)
    ans = targets_private  # private targets = the subset unique to this sample
    # Build the candidate-sites tables. Show BOTH the full personal-candidate set (all
    # off-targets the sample carries) and the private subset (unique to the sample).
    # Previously only the private table was shown, so a sample with personal-but-not-
    # private candidates (the common case) saw an empty "candidates" table.
    _risk_sample_cols = [
        "Variant_samples_(highest_CFD)",
        "Variant_samples_(fewest_mm+b)",
        "Variant_samples_(highest_CRISTA)",
    ]

    def _risk_table(df, tid):
        return dash_table.DataTable(
            css=[{"selector": ".row", "rule": "margin: 0"}],
            id=tid,
            columns=[{"name": i, "id": i, "hideable": True} for i in df.columns],
            data=df.to_dict("records"),
            style_cell_conditional=[
                {
                    "if": {"column_id": c},
                    "textAlign": "left",
                    "minWidth": "180px",
                    "width": "180px",
                    "maxWidth": "180px",
                    "overflow": "hidden",
                }
                for c in _risk_sample_cols
            ],
            style_table={
                "overflowX": "scroll",
                "overflowY": "scroll",
                "max-height": "300px",
            },
        )

    candidate_tables = [
        html.H5(f"Personal candidate off-targets ({targets_personal.shape[0]})"),
        _risk_table(targets_personal, "results-table-risk-personal"),
        html.Br(),
        html.H5(f"Private off-targets — unique to this sample ({ans.shape[0]})"),
        _risk_table(ans, "results-table-risk"),
    ]
    # put images for personal and private targets in HTML
    try:
        image_personal_top = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    os.path.join(
                        imgsdir,
                        f"CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.{sample}.personal.png",
                    ),
                    mode="rb",
                ).read()
            ).decode()
        )
        image_private_top = "data:image/png;base64,{}".format(
            base64.b64encode(
                open(
                    os.path.join(
                        imgsdir,
                        f"CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.{sample}.private.png",
                    ),
                    mode="rb",
                ).read()
            ).decode()
        )
    except:
        raise ValueError("Personal and Private Lolliplots not found")
    # recover filtering criterion selected via drop-down bar
    filter_criterion = read_json(job_id)
    assert filter_criterion in FILTERING_CRITERIA
    # create personal risk card page
    try:
        file_to_load = f"{job_id}.{sample}.{guide}.private_targets.zip"
        out_1 = [
            html.A(
                "Download private targets",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            False,
            [
                html.P(f"Top 100 Personal Targets ordered by {filter_criterion}"),
                html.A(
                    html.Img(
                        src=image_personal_top,
                        id="sample-personal-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            [
                html.P(f"Top 100 Private Targets ordered by {filter_criterion}"),
                html.A(
                    html.Img(
                        src=image_private_top,
                        id="sample-private-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="results-table",
                columns=[{"name": i, "id": i} for i in sample_stats.columns],
                data=sample_stats.to_dict("records"),
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
            ),
            candidate_tables,
        ]
    except:
        out_1 = [
            html.A(
                "Download private targets",
                href=os.path.join(URL, RESULTS_DIR, job_id, file_to_load),
                target="_blank",
            ),
            True,
            [
                html.P(f"Top 100 Personal Targets ordered by {filter_criterion}"),
                html.A(
                    html.Img(
                        src=image_personal_top,
                        id="sample-personal-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            [
                html.P(f"Top 100 Private Targets ordered by {filter_criterion}"),
                html.A(
                    html.Img(
                        src=image_private_top,
                        id="sample-private-top",
                        width="100%",
                        height="auto",
                    ),
                    target="_blank",
                ),
            ],
            dash_table.DataTable(
                css=[{"selector": ".row", "rule": "margin: 0"}],
                id="results-table",
                style_cell_conditional=[
                    {
                        "if": {"column_id": "Variant_samples_(highest_CFD)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(fewest_mm+b)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                    {
                        "if": {"column_id": "Variant_samples_(highest_CRISTA)"},
                        "textAlign": "left",
                        "minWidth": "180px",
                        "width": "180px",
                        "maxWidth": "180px",
                        "overflow": "hidden",
                    },
                ],
                columns=[{"name": i, "id": i} for i in sample_stats.columns],
                data=sample_stats.to_dict("records"),
            ),
            [],
        ]
    return out_1


# ------------------------------------------------------------------------------
# main page layout


# update the main content table
@app.callback(
    Output("div-tab-content", "children"),
    [
        Input("tabs-reports", "value"),
        Input("general-profile-table", "selected_cells"),
        Input("target_filter_dropdown", "value"),
    ],
    [
        State("general-profile-table", "data"),
        State("url", "search"),
        State("div-genome-type", "children"),
    ],
)
def update_content_tab(
    value: str,
    sel_cel: List,
    filter_criterion: str,
    all_guides: List,
    search: str,
    genome_type: str,
) -> List:
    """Build and update the layout of the results page.

    ...

    Parameters
    ----------
    value : str
        HTML identifier
    sel_cel : List
        Selected cells
    filetr_criterion : str
        Filtering criterion
    all_guides : List
        All CRISPR guides
    search : str
        Search
    genome_type : str
        Selected genome type

    Returns
    -------
    List
        Results page layout
    """

    if value is not None:
        if not isinstance(value, str):
            raise TypeError(f"Expected {str.__name__}, got {type(value).__name__}")
    if not isinstance(filter_criterion, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_criterion).__name__}"
        )
    if filter_criterion not in FILTERING_CRITERIA:
        raise ValueError(f"Forbidden filtering criterion selected ({filter_criterion})")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(genome_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome_type).__name__}")
    if value is None or sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate  # do not do anything
    # recover current guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # recover job ID
    job_id = search.split("=")[-1]
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    # read parameters file
    try:
        with open(os.path.join(job_directory, PARAMS_FILE)) as handle_params:
            params = handle_params.read()
            mms = (next(s for s in params.split("\n") if "Mismatches" in s)).split(
                "\t"
            )[-1]
            genome_selected = (
                next(s for s in params.split("\n") if "Genome_selected" in s)
            ).split("\t")[-1]
            max_bulges = (
                next(s for s in params.split("\n") if "Max_bulges" in s)
            ).split("\t")[-1]
            pam = (next(s for s in params.split("\n") if "Pam" in s)).split("\t")[-1]
            nuclease = (next(s for s in params.split("\n") if "Nuclease" in s)).split(
                "\t"
            )[-1]
    except OSError as e:
        raise e
    # initialize page layout list
    fl = []
    fl.append(html.Br())
    # check the selected nuclease
    if nuclease != "SpCas9":
        CFD_notification = html.Div(
            "CFD score is not calculated if the used nuclease is not SpCas9"
        )
        # if nuclease is not SpCas9 filter only by fewest mm + bulges
        filter_criterion = FILTERING_CRITERIA[0]  # fewest mm + b
    else:
        CFD_notification = html.Div("", hidden=True)
    if nuclease != CAS9 and filter_criterion != FILTERING_CRITERIA[0]:
        raise ValueError(f"Wrong filtering criterion selected for nuclease {nuclease}")
    # PAM(s)
    pam_at_start = False
    assert isinstance(guide, str)
    if guide[0] == "N":
        pam_at_start = True
    if pam_at_start:
        fl.append(html.H5(f"Focus on: {pam}{guide.replace('N', '')}"))
    else:
        fl.append(html.H5(f"Focus on: {guide.replace('N', '')}{pam}"))
    # BUG changing the selected guide two times, drop mms to 0
    # TODO: add hidden div ?
    if value == "tab-summary-by-guide":
        # Show Summary by Mismatches/Bulges
        # NOTE use old id (summary-by-guide -> summary by mm+b)
        fl.append(
            html.P(
                [
                    str(
                        "Summary table counting the number of targets found in "
                        "the Reference and Variant Genome for each combination "
                        'of Bulge Type, Bulge Size and Mismatch. Select "Show '
                        'Targets" to view the corresponding list of targets.'
                    ),
                ]
            )
        )
        fl.append(html.Br())
        # load summary by guide table
        guides_summary = pd.read_csv(
            os.path.join(
                job_directory,
                f"{job_id}.summary_by_guide.{guide}_{filter_criterion}.txt",
            ),
            sep="\t",
            na_filter=False,
        )
        more_info_col = []
        total_col = []
        for _ in range(guides_summary.shape[0]):
            more_info_col.append("Show Targets")
            total_col.append(guides_summary["Bulge Size"])
        guides_summary[""] = more_info_col
        fl.append(
            html.Div(
                generate_table(
                    guides_summary,
                    "table-summary-by-guide",
                    genome_type,
                    guide,
                    job_id,
                ),
                style={"text-align": "center"},
            )
        )
        return fl
    elif value == "tab-summary-by-sample":
        # Show Summary by Sample table
        fl.append(
            html.P(
                str(
                    "Summary table counting the number of targets found in the "
                    "Variant Genome for each sample. Filter the table by "
                    "selecting the Population or Superpopulation desired from "
                    "the dropdowns."
                )
            )
        )
        if genome_type == "both":
            col_names_sample = [
                "Sample",
                "Sex",
                "Population",
                "Super Population",  # 'Targets in Reference',
                "Targets in Sample",
                "Targets in Population",
                "Targets in Super Population",
                "PAM Creation",
            ]
        else:
            col_names_sample = [
                "Sample",
                "Sex",
                "Population",
                "Super Population",  # 'Targets in Reference',
                "Targets in Sample",
                "Targets in Population",
                "Targets in Super Population",
                "PAM Creation",
            ]
        # load summary by samples table
        samples_summary = pd.read_csv(
            os.path.join(
                job_directory,
                f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt",
            ),
            sep="\t",
            names=col_names_sample,
            skiprows=2,
            na_filter=False,
        )
        samples_summary = samples_summary.sort_values(
            "Targets in Sample", ascending=False
        )
        more_info_col = ["Show Targets" for _ in range(samples_summary.shape[0])]
        samples_summary[""] = more_info_col

        population_1000gp = associateSample.loadSampleAssociation(
            os.path.join(job_directory, SAMPLES_ID_FILE)
        )[2]
        super_populations = [{"label": i, "value": i} for i in population_1000gp.keys()]
        populations = []
        for pop in population_1000gp.keys():
            for i in population_1000gp[pop]:
                populations.append({"label": i, "value": i})
        fl.append(
            html.Div(
                [
                    html.Div(
                        os.path.join(
                            job_directory, f"{job_id}.summary_by_samples.{guide}"
                        ),
                        style={"display": "none"},
                        id="div-info-summary_by_sample",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=super_populations,
                                        id="dropdown-superpopulation-sample",
                                        placeholder="Select a Super Population",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=populations,
                                        id="dropdown-population-sample",
                                        placeholder="Select a Population",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        id="input-sample",
                                        placeholder="Select a Sample",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button(
                                        "Filter",
                                        id="button-filter-population-sample",
                                    )
                                )
                            ),
                        ]
                    ),
                    dbc.Row(
                        dbc.Col(
                            html.Div(
                                [
                                    html.P(
                                        "Generating download link, Please wait...",
                                        id="download-link-summary_by_sample",
                                    ),
                                    dcc.Interval(
                                        interval=(1 * 1000),
                                        id="interval-summary_by_sample",
                                    ),
                                ]
                            )
                        )
                    ),
                ],
                style={"width": "50%"},
            )
        )
        fl.append(
            html.Div(
                "None,None,None",
                id="div-sample-filter-query",
                style={"display": "none"},
            )
        )  # keep current filter:  Superpop + Pop
        fl.append(
            html.Div(
                generate_table_samples(
                    samples_summary, "table-samples", 1, guide, job_id
                ),
                style={"text-align": "center"},
                id="div-table-samples",
            )
        )
        fl.append(
            html.Div(
                [
                    html.Button("Prev", id="prev-page-sample"),
                    html.Button("Next", id="next-page-sample"),
                ],
                style={"text-align": "center"},
            )
        )
        max_page = samples_summary.shape[0]
        max_page = math.floor(max_page / 10) + 1
        fl.append(html.Div(f"1/{max_page}", id="div-current-page-table-samples"))
        return fl
    elif value == "tab-summary-by-position":
        # Show Summary by position table (Query Genomic regions tab)
        fl.append(
            html.P(
                str(
                    "Summary table containing all the targets found in a "
                    "specific range of positions (chr, start, end) of the "
                    "genome."
                )
            )
        )
        fl.append(
            html.P(
                str(
                    "Filter the table by selecting the chromosome of interest, "
                    "and writing the start and end position of the region to "
                    "view."
                )
            )
        )
        # Dropdown chromosomes
        try:
            # read chromosomes from FASTA files
            onlyfile = [
                f
                for f in os.listdir(
                    os.path.join(current_working_directory, "Genomes", genome_selected)
                )
                if (
                    os.path.isfile(
                        os.path.join(current_working_directory, "Genomes", genome_selected, f)
                    )
                    and (f.endswith(".fa") or f.endswith(".fasta"))
                )
            ]
        except:
            # guess chromosomes
            onlyfile = [f"chr{i}.fa" for i in range(1, 23)]
            # NOTE in case no chr in "Genomes", put 22 chr + X Y M
            onlyfile += ["chrX.fa", "chrY.fa", "chrM"]
        # remove file extension (.fa)
        onlyfile = [x[: x.rfind(".")] for x in onlyfile]
        chr_file = []
        chr_file_unset = []
        for chr_name in onlyfile:
            chr_name = chr_name.replace(".enriched", "")
            if "_" in chr_name:
                chr_file_unset.append(chr_name)
            else:
                chr_file.append(chr_name)
        chr_file.sort(
            key=lambda s: [
                int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)
            ]
        )
        chr_file_unset.sort(
            key=lambda s: [
                int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)
            ]
        )
        chr_file += chr_file_unset
        chr_file = [{"label": chr_name, "value": chr_name} for chr_name in chr_file]
        # TODO: insert failsafe if no chromosome is found
        fl.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        options=chr_file,
                                        id="dropdown-chr-table-position",
                                        placeholder="Select a chromosome",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        placeholder="Start Position",
                                        id="input-position-start",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    dcc.Input(
                                        placeholder="End Position",
                                        id="input-position-end",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button("Filter", id="button-filter-position")
                                )
                            ),
                            html.Br(),
                        ]
                    ),
                ],
                style={"width": "50%"},
            )
        )
        # keep current filter:  chr,pos_start,pos_end
        fl.append(
            html.Div(
                "None,None,None",
                id="div-position-filter-query",
                style={"display": "none"},
            )
        )
        fl.append(html.Br())
        fl.append(html.Div(style={"text-align": "center"}, id="div-table-position"))
        max_page = 1  # maximum one single page
        fl.append(html.Div(f"1/{max_page}", id="div-current-page-table-position"))
        fl.append(
            html.Div(
                f"{mms}-{max_bulges}",
                id="div-mms-bulges-position",
                style={"display": "none"},
            )
        )
        return fl
    elif value == "tab-graphical-sample-card":
        # Show Personal Risk Cards table
        samples_summary = pd.read_csv(
            os.path.join(
                job_directory,
                f"{job_id}.summary_by_samples.{guide}_{filter_criterion}.txt",
            ),
            skiprows=2,
            sep="\t",
            header=None,
            na_filter=False,
        )
        samples = samples_summary.iloc[:, 0]  # samples on 1st column
        fl.append(
            html.P(
                str(
                    "Summary page containing the single Personal Risk card to "
                    "be inspected and downloaded"
                )
            )
        )
        fl.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    dcc.Dropdown(
                                        id="dropdown-sample-card",
                                        options=[
                                            {"label": sam, "value": sam}
                                            for sam in samples
                                        ],
                                        placeholder="Select a Sample",
                                    )
                                )
                            ),
                            dbc.Col(
                                html.Div(
                                    html.Button("Generate", id="button-sample-card")
                                )
                            ),
                            dbc.Col(
                                html.Div(id="download-link-personal-card", hidden=True)
                            ),
                        ]
                    ),
                ],
                style={"width": "50%"},
            )
        )
        fl.append(
            html.Div(
                [
                    html.Br(),
                    dbc.Row(
                        [
                            dbc.Col(html.Div("", id="div-personal-plot")),
                            dbc.Col(html.Div("", id="div-private-plot")),
                        ]
                    ),
                ]
            )
        )
        fl.append(
            html.Div(
                "",
                id="div-table-sample-card",
                style={
                    "text-align": "center",
                    "margin-left": "1%",
                    "margin-right": "1%",
                },
            )
        )
        fl.append(
            html.Div(
                "",
                id="div-top-target-sample-card",
                style={
                    "text-align": "center",
                    "margin-left": "1%",
                    "margin-right": "1%",
                },
            )
        )
        fl.append(html.Div("", id="div-sample-card"))
        return fl
    elif value == "tab-query-table":
        # Show Custom Ranking table
        fl.append(
            html.P(
                str(
                    "Summary page to query the final result file selecting "
                    "one/two column to group by the table and extract "
                    "requested targets"
                )
            )
        )
        all_value = {
            "Target1 :with highest CFD": [
                "Mismatches",
                "Bulges",
                "Mismatches+bulges",
                "CFD_score",
                "CFD_risk_Score",
            ],  # , 'Highest_CFD_Absolute_Risk_Score'
            "Target2 :with lowest Mismatches + Bulge Count": [
                "Mismatches",
                "Bulges",
                "Mismatches+bulges",
                "CFD_score",
                "CFD_risk_Score",
            ],
        }  # , 'CFD_Absolute_Risk_Score'
        all_options = {
            "Target1 :with highest CFD": [
                " Mismatches",
                " Bulges",
                " Mismatch+Bulges",
                " CFD",
                " Risk Score",
            ],  # , ' Absolute Risk Score'
            "Target2 :with lowest Mismatches + Bulges Count": [
                " Mismatches",
                " Bulges",
                " Mismatch+Bulges",
                " CFD",
                " Risk Score",
            ],
        }  # , ' Absolute Risk Score'
        label = [{"label": lab} for lab in all_options.keys()]
        value = [{"value": val} for val in all_value.keys()]
        target_opt = [label, value]
        query_tab_content = html.Div(
            [
                # row with the first and second group by and thresholds
                dbc.Row(
                    [
                        dbc.Col(  # col0 phantom target select
                            [
                                html.Div(
                                    [
                                        html.H4("Order by"),
                                        dcc.RadioItems(
                                            id="target",
                                            options=target_opt,
                                            value="Target1 :with highest CFD",
                                        ),
                                    ]
                                )
                            ],
                            style={"display": "none"},
                        ),
                        dbc.Col(  # col1 main group by
                            html.Div(
                                [
                                    html.H4("Group by"),
                                    dcc.RadioItems(id="order", value="CFD_score"),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(  # col2 second group by
                            html.Div(
                                [
                                    html.H4("And group by"),
                                    html.P(
                                        "First select the left group by value",
                                        id="secondtext",
                                    ),
                                    dcc.RadioItems(id="multiorder"),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(  # select threshold
                            html.Div(
                                [
                                    html.H4("Select thresholds"),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6("Min"),
                                                            dcc.Dropdown(
                                                                id="thresh_drop"
                                                            ),
                                                        ]
                                                    ),
                                                ]
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6("Max"),
                                                            dcc.Dropdown(id="maxdrop"),
                                                        ]
                                                    )
                                                ]
                                            ),
                                        ]
                                    ),
                                ]
                            ),
                            width=3,
                        ),
                        dbc.Col(
                            html.Div(
                                [
                                    html.H4("Select ordering"),
                                    dcc.RadioItems(
                                        id="Radio-asc-1",
                                        options=[
                                            {"label": " Ascending", "value": "ASC"},
                                            {"label": " Descending", "value": "DESC"},
                                        ],
                                        value="DESC",
                                        labelStyle={
                                            "display": "inline-block",
                                            "margin": "10px",
                                        },
                                    ),
                                ]
                            ),
                            width=3,
                        ),
                    ],
                    justify="center",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button(
                                        "Submit",
                                        id="submit-val",
                                        n_clicks=0,
                                    ),
                                )
                            ],
                            width={"size": 1},
                        ),
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button(
                                        "Reset",
                                        id="reset-val",
                                        n_clicks=0,
                                    )
                                )
                            ],
                            width={"size": 1, "offset": 1},
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Br(),
                                        html.Hr(),
                                    ]
                                )
                            ]
                        ),
                    ],
                    justify="center",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                CFD_notification,
                                html.P(
                                    str(
                                        "Export will download 1000 lines "
                                        "contained in the current view of the "
                                        "table"
                                    )
                                ),
                                html.Div(
                                    dash_table.DataTable(
                                        css=[
                                            {
                                                "word-break": "break-all",
                                                "line-break": "anywhere",
                                                "overflow-wrap": "break-word",
                                                "selector": ".row",
                                                "rule": "margin: 0; overflow: inherit; word-break: break-all; overflow-wrap: break-word; line-break: anywhere;",
                                            }
                                        ],
                                        style_cell={
                                            "height": "auto",
                                            "textAlign": "left",
                                        },
                                        export_format="xlsx",
                                        id="live_table",
                                        style_cell_conditional=[
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(highest_CFD)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(fewest_mm+b)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                            {
                                                "if": {
                                                    "column_id": "Variant_samples_(highest_CRISTA)"
                                                },
                                                "textAlign": "left",
                                                "minWidth": "180px",
                                                "width": "180px",
                                                "maxWidth": "180px",
                                                "overflow": "hidden",
                                            },
                                        ],
                                        style_table={
                                            "overflowX": "scroll",
                                            "overflowY": "scroll",
                                            "max-height": "300px",
                                        },
                                        page_current=0,
                                        page_size=1000,
                                        page_action="custom",
                                        tooltip_delay=0,
                                        tooltip_duration=None,
                                    ),
                                    id="div-query-table",
                                ),
                            ],
                        ),
                    ],
                ),
                html.Div(
                    [
                        dbc.Row(
                            dbc.Col(
                                [
                                    dbc.Alert(
                                        "Select a main order before submitting the query",
                                        id="message-alert",
                                        color="danger",
                                        dismissable=True,
                                        fade=True,
                                        is_open=False,
                                        duration=4000,
                                    ),
                                ]
                            ),
                        )
                    ],
                    style={"display": "inline-block"},
                ),
            ]
        )
        fl.append(query_tab_content)
        return fl
    else:
        # Show Graphical Report images
        samp_style = {}
        if genome_type == "ref":
            samp_style = {"display": "none"}
        fl.append(
            html.P(
                str(
                    "Summary Graphical report collecting all the plots and "
                    "images produced during the search"
                )
            )
        )
        total = int(mms) + int(max_bulges)
        opt_mm = [{"label": str(i), "value": str(i)} for i in range(total + 1)]
        opt_blg = [
            {"label": str(i), "value": str(i)} for i in range(int(max_bulges) + 1)
        ]
        if genome_type != "ref":
            population_1000gp = associateSample.loadSampleAssociation(
                os.path.join(job_directory, SAMPLES_ID_FILE)
            )[2]
            super_populations = [
                {"label": i, "value": i} for i in population_1000gp.keys()
            ]
            populations = []
            for k in population_1000gp.keys():
                for i in population_1000gp[k]:
                    populations.append({"label": i, "value": i})
        else:
            super_populations = []
            populations = []
        # embed a top-1000 plot PNG if it exists (returns None if absent, e.g. the
        # variant-effect plot is only produced for the CFD/CRISTA score criteria)
        def _embed_top1000(fname, img_id):
            try:
                with open(
                    os.path.join(
                        current_working_directory, RESULTS_DIR, job_id, IMGS_DIR, fname
                    ),
                    mode="rb",
                ) as _fh:
                    _src = "data:image/png;base64,{}".format(
                        base64.b64encode(_fh.read()).decode()
                    )
            except Exception:
                return None
            return html.A(
                html.Img(src=_src, id=img_id, width="80%", height="auto"),
                target="_blank",
            )

        _score_img = _embed_top1000(
            f"CRISPRme_{filter_criterion}_top_1000_log_for_main_text_{guide}.png",
            "top-1000-score",
        )
        _delta_img = _embed_top1000(
            f"CRISPRme_{filter_criterion}_top_1000_by_variant_effect_{guide}.png",
            "top-1000-delta",
        )
        _t1000_parts = []
        if _score_img is not None:
            _t1000_parts += [
                html.H5("Top 1000 candidate off-targets (ranked by score)"),
                _score_img,
            ]
        if _delta_img is not None:
            _t1000_parts += [
                html.Br(),
                html.H5("Top 1000 by variant effect (largest |ALT-REF| score change)"),
                html.P(
                    "Ranked by how much each variant changes the off-target score "
                    "(largest first), so the impactful variants are not buried among "
                    "the many that leave the score unchanged.",
                    style={"font-size": "1.05rem", "color": "#555"},
                ),
                _delta_img,
            ]
        top1000_image = html.Div(_t1000_parts) if _t1000_parts else html.Div("")
        total_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P(
                            str(
                                "Select total number of mismatches and/or "
                                "bulges to consider, up to"
                            )
                        ),
                        dcc.Dropdown(
                            id="mm-dropdown",
                            # default to the full budget (mm + bulges): the population
                            # barplot is cumulative ("up to N"), and the pipeline only
                            # emits a PNG for a total that has data, so the max total is
                            # the all-inclusive view and is non-empty whenever the search
                            # found any off-target (vs the old "0" = on-target-only, which
                            # showed "No result" for typical sparse guides).
                            options=opt_mm,
                            value=str(total),
                            clearable=False,
                        ),
                    ]
                )
            )
        ]
        sample_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Superpopulation", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                options=super_populations,
                                id="dropdown-superpopulation-sample",
                                placeholder="SuperPopulation",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Population", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                options=populations,
                                id="dropdown-population-sample",
                                placeholder="Population",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Sample", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                id="dropdown-sample",
                                placeholder="Sample",
                                style=samp_style,
                            ),
                        ),
                    ]
                ),
                md=4,
            ),
        ]
        fl.append(
            html.Div(
                [
                    CFD_notification,
                    dbc.Row(dbc.Col(top1000_image, width={"size": 10, "offset": 2})),
                    dbc.Row(total_buttons, justify="center"),
                    html.Br(),
                ]
            )
        )
        radar_chart_encode_gencode = dbc.Col(
            html.Div(id="div-radar-chart-encode_gencode")
        )
        # Always include the population-barplot container so the callback Output
        # ("div-population-barplot") always resolves. For reference-only searches
        # (genome_type == "ref") there is no per-population variant data, so hide it
        # instead of dropping it from the layout -- dropping it made the whole
        # update_images_tabs callback fail (the radar chart died too) on ref searches.
        populations_barplots = dbc.Col(
            html.Div(id="div-population-barplot"),
            style=({"display": "none"} if genome_type == "ref" else {}),
        )
        graph_summary_both = [populations_barplots, radar_chart_encode_gencode]
        fl.append(html.Div([dbc.Row(graph_summary_both)]))
        fl.append(
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            "The GENCODE and ENCODE annotations are defined in detail ",
                            html.A(
                                "here",
                                target="_blank",
                                href="https://www.gencodegenes.org/human/",
                            ),
                            " and ",
                            html.A(
                                "here",
                                target="_blank",
                                href="https://screen.encodeproject.org/",
                            ),
                        ]
                    ),
                    width={"size": 6, "offset": 6},
                )
            )
        )
        cfd_path = os.path.join(job_directory, f"{job_id}.CFDGraph.txt")
        if not os.path.isfile(cfd_path):  # No file found to display CFD graph
            return fl
        fl.extend(CFDGraph.CFDGraph(cfd_path))
        return fl
    raise PreventUpdate


# TODO: move auxiliary functions close to each other in this file
# Perform expensive loading of a dataframe and save result into 'global store'
# Cache are in the Cache directory


@cache.memoize()
def global_store(job_id: str) -> pd.DataFrame:
    """Perform once dataframe loading and cache data in Cache directory.

    ...

    Parameters
    ----------
    value : str
        Job ID

    Returns
    -------
    pd.DataFrame
    """

    if job_id is not None and not isinstance(job_id, str):
        raise TypeError(f"Expected {str.__name__}, got {type(job_id).__name__}")
    if job_id is None:
        return ""  # nothing to return
    target = [
        f
        for f in os.listdir(os.path.join(current_working_directory, RESULTS_DIR, job_id))
        if os.path.isfile(os.path.join(current_working_directory, RESULTS_DIR, job_id, f))
        and f.endswith("scores.txt")
    ]
    # use targets file
    if not target:
        target = [
            f
            for f in os.listdir(os.path.join(current_working_directory, RESULTS_DIR, job_id))
            if os.path.isfile(os.path.join(current_working_directory, RESULTS_DIR, job_id, f))
            and f.endswith("targets.txt")
        ]
    targets_summary = pd.read_csv(
        os.path.join(current_working_directory, RESULTS_DIR, job_id, target[0]),
        sep="\t",
        usecols=range(38),
        na_filter=False,
    )
    targets_summary.rename(
        columns={
            "#Bulge type": "BulgeType",
            "#Bulge_type": "BulgeType",
            "Bulge Size": "BulgeSize",
            "Bulge_Size": "BulgeSize",
            "Doench 2016": "Doench2016",
            "Doench_2016": "Doench2016",
        },
        inplace=True,
    )
    return targets_summary


# trigger tables update
@app.callback(
    Output("result-table", "data"),
    [
        Input("result-table", "page_current"),
        Input("result-table", "page_size"),
        Input("result-table", "sort_by"),
        Input("result-table", "filter_query"),
    ],
    [State("url", "search"), State("url", "hash")],
)
def update_table(
    page_current: int,
    page_size: int,
    sort_by: str,
    filter_term: str,
    search: str,
    hash_guide: str,
) -> Dict:
    """Split the results according to a filtering or sorting criterion selected
    by the user.

    Update the shown results once the user clicks on the "next page"/"prev page"
    buttons.

    Load the targets or scores (if available) files, and store it in a pandas
    DataFrame object. The column names are changed in order to match those
    of the table displayed within the webpage.

    If no targets are found a warning message is returned.

    ...

    Parameters
    ----------
    page_current : int
        Current page
    page_size : int
        Page size
    sort_by : List
        Sorting criterion
    filter_term : str
        Filtering
    search : str
        Search
    hash_guide : str
        Guide hashing

    Returns
    -------
    Dict
    """

    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(sort_by, list):
        raise TypeError(f"Expected {list.__name__}, got {type(sort_by).__name__}")
    if not isinstance(filter_term, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_term).__name__}")
    if search is not None:
        if not isinstance(search, str):
            raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_guide).__name__}")
    if search is None:
        raise PreventUpdate  # do not do anything
    # recover job ID
    job_id = search.split("=")[-1]
    # recover job directory
    job_directory = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    # recover guide
    guide = hash_guide.split("#")[1]
    filtering_expressions = filter_term.split(" && ")
    # keep data for the current guide
    df = global_store(job_id)
    df_filtered = df[df["crRNA"] == guide]
    # insert sorting criteria
    sort_by.insert(0, {"column_id": "Mismatches", "direction": "asc"})
    sort_by.insert(1, {"column_id": "BulgeSize", "direction": "asc"})
    # get filters
    for f in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(f)
        if operator in PANDAS_OPERATORS:
            # these operators match pandas series operator method names
            df_filtered = df_filtered.loc[
                getattr(df_filtered[col_name], operator)(filter_value)
            ].sort_values(
                [col["column_id"] for col in sort_by],
                ascending=[col["direction"] == "asc" for col in sort_by],
                inplace=False,
            )
        elif operator == "contains":
            df_filtered = df_filtered.loc[
                df_filtered[col_name].str.contains(filter_value)
            ]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            df_filtered = df_filtered.loc[
                df_filtered[col_name].str.startswith(filter_value)
            ]
    if sort_by:
        df_filtered = df_filtered.sort_values(
            [
                "Samples" if col["column_id"] == "Samples Summary" else col["column_id"]
                for col in sort_by
            ],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )
    # Check if we have some results
    warning_no_res = ""
    try:
        with open(
            os.path.join(job_directory, f"{job_id}.targets.txt")
        ) as handle_targets:
            no_result = False
            handle_targets.readline()  # consume buffer
            line = handle_targets.readline().strip()  # last line
            if not line:
                no_result = True
    except OSError as e:
        raise e
    if no_result:  # display warning message
        warning_no_res = dbc.Alert(
            "No results were found with the given parameters", color="warning"
        )
    return df_filtered.iloc[
        page_current * page_size : (page_current + 1) * page_size
    ].to_dict("records")


# ------------------------------------------------------------------------------
# Callbacks for querying part


# Return the table with the query's result
@app.callback(
    # [Output('live_table', 'data'),
    [
        Output("live_table", "columns"),
        Output("live_table", "data"),
        Output("live_table", "tooltip_data"),
        Output("message-alert", "is_open"),
    ],
    [
        Input("submit-val", "n_clicks"),
        Input("live_table", "page_current"),
        Input("target_filter_dropdown", "value"),
    ],  # take this value (as state)
    [
        State("live_table", "page_size"),
        State("general-profile-table", "selected_cells"),
        State("target", "value"),
        State("order", "value"),
        State("general-profile-table", "data"),
        State("multiorder", "value"),
        State("thresh_drop", "value"),
        State("Radio-asc-1", "value"),
        State("maxdrop", "value"),
        State("url", "search"),
        State("message-alert", "is_open"),
    ],
)
def update_output(
    n_clicks: int,
    page_current: int,
    filter_target_value: str,
    page_size: int,
    sel_cel: List,
    target: str,
    radio_order: str,
    all_guides: List,
    order_drop: str,
    thresh_drop: str,
    asc1: str,
    maxdrop: int,
    url: str,
    alert: bool,
) -> Tuple:
    """Update the dispalyed table according to the query performed by the user.

    ...

    Paramters
    ---------
    n_clicks : int
        Clicks
    page_current : int
        Current page
    filter_target_values : str
        Targets filter
    page_size : int
        Page size
    sel_cel : List
        Selected cells
    target : str
        Current target
    radio_order : str
        First group by criterion
    all_guides : List
        All CRISPR guides
    order_drop : str
        Second group by criterion
    thresh_drop : str
        Threshold value
    asc1 : str
        Sorting in ascending or descending order
    maxdrop : int
        Maximum number of dropped rows
    alert : bool
        Alert

    Returns
    -------
    Tuple
    """

    if not isinstance(n_clicks, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n_clicks).__name__}")
    if not isinstance(page_current, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_current).__name__}")
    if not isinstance(filter_target_value, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(filter_target_value).__name__}"
        )
    if not isinstance(page_size, int):
        raise TypeError(f"Expected {int.__name__}, got {type(page_size).__name__}")
    if not isinstance(target, str):
        raise TypeError(f"Expected {str.__name__}, got {type(target).__name__}")
    # prevent update on None inputs
    if radio_order is None or (
        order_drop is None and thresh_drop is None and asc1 is None
    ):
        raise PreventUpdate  # do not do anything
    # recover guide
    guide = all_guides[int(sel_cel[0]["row"])]["Guide"]
    # target is the filter value to query on the db
    target = filter_target_value
    if n_clicks > 0:
        # no input by user
        if radio_order == None:
            data = []
            tooltip_data = []
            return (data, tooltip_data, not alert)
        else:
            # perform queries on data
            if thresh_drop != None:
                alert = False
                data = query_manager.shold(
                    target,
                    n_clicks,
                    page_current,
                    page_size,
                    radio_order,
                    order_drop,
                    thresh_drop,
                    maxdrop,
                    asc1,
                    url,
                    guide,
                    current_working_directory,
                )
            else:
                data = query_manager.noshold(
                    target,
                    n_clicks,
                    page_current,
                    page_size,
                    radio_order,
                    order_drop,
                    asc1,
                    url,
                    guide,
                    current_working_directory,
                )
            # find columns to drop (use user's filter)
            drop_cols = drop_columns(data, filter_target_value)
            # drop column from datatable to show
            data.drop(drop_cols, inplace=True, axis=1)
            # extract cols for datatable
            columns = [
                {"name": i, "id": i, "hideable": True}
                for count, i in enumerate(data.columns)
            ]
            # select SNPs columns to filter
            if filter_target_value == FILTERING_CRITERIA[0]:  # fewest
                snps = pd.DataFrame(data[VARIANTS_FEWEST]).to_dict("records")
            if filter_target_value == FILTERING_CRITERIA[1]:
                snps = pd.DataFrame(data[VARIANTS_CFD]).to_dict("records")
            if filter_target_value == FILTERING_CRITERIA[2]:
                snps = pd.DataFrame(data[VARIANTS_CRISTA]).to_dict("records")
            # extract data and list datas
            data = data.to_dict("records")
            tooltip_data = [
                {
                    column: {"value": str(value), "type": "markdown"}
                    for column, value in row.items()
                }
                for row in snps
            ]
    else:
        raise PreventUpdate  # do not do anything
    return (columns, data, tooltip_data, alert)


# trigger page number reset
@app.callback(Output("live_table", "page_current"), [Input("submit-val", "n_clicks")])
def reset_pagenumber(n: int) -> int:
    """Reset page number.

    ...

    Paramters
    ---------
    n : int
        Current page number

    Returns
    -------
    Reset page number
    """

    if not isinstance(n, int):
        raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if n > 0:
        number_reset = 0
        return number_reset
    else:
        raise PreventUpdate  # page number already reset


# trigger columns options selection
@app.callback(Output("order", "options"), [Input("target", "value")])
def set_columns_options(selected_target: str) -> List[Dict]:
    """Set options to be selected by the user.

    ...

    Parameters
    ----------
    selected_target : str
        Selected columns

    Returns
    -------
    List[Dict]
    """

    if not isinstance(selected_target, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(selected_target).__name__}"
        )
    # all possible column values
    all_value = {
        "Target1 :with highest CFD": [
            "Mismatches",
            "Bulges",
            "Mismatches+bulges",
            "CFD_score",
            "CFD_risk_score",
        ],  # , 'Highest_CFD_Absolute_Risk_Score'
        "Target2 :with lowest Mismatches + Bulge Count": [
            "Mismatches",
            "Bulges",
            "Mismatches+bulges",
            "CFD_score",
            "CFD_risk_score",
        ],
    }  # , 'CFD_Absolute_Risk_Score'
    # all possible columns options
    all_options = {
        "Target1 :with highest CFD": [
            " Mismatches",
            " Bulges",
            " Mismatch+Bulges",
            " Score",
            " Risk Score",
        ],  # , ' Absolute Risk Score'
        "Target2 :with lowest Mismatches + Bulges Count": [
            " Mismatches",
            " Bulges",
            " Mismatch+Bulges",
            " CFD",
            " Risk Score",
        ],
    }  # , ' Absolute Risk Score'
    gi = [
        {
            "label": all_options[selected_target][count],
            "value": all_value[selected_target][count],
        }
        for count in range(len(all_value[selected_target]))
    ]
    return gi


# callback to return the parameters in the various cases
@app.callback(
    [
        Output("multiorder", "options"),
        Output("thresh_drop", "options"),
        Output(component_id="secondtext", component_property="style"),
    ],
    [Input("order", "value")],
)
def set_display_children(selected_order: str) -> Tuple:
    """Display table options.

    ...

    Parameters
    ----------
    selected_order : str
        Selected ordering

    Returns
    -------
    Tuple
    """

    if selected_order is not None:
        if not isinstance(selected_order, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(selected_order).__name__}"
            )
    target_value = {
        "Mismatches": ["Bulges", "Mismatches+bulges", "CFD"],
        "Bulges": ["Mismatches", "Mismatches+bulges", "CFD_score"],
        "Mismatches+bulges": ["Mismatches", "Bulges", "CFD_score"],
        "CFD_score": ["Mismatches", "Bulges", "Mismatches+bulges"],
        "CFD_risk_score": [],
    }
    target_label = {
        "Mismatches": [" Bulges", " Mismatch+Bulges", " Score"],
        "Bulges": [" Mismatches", " Mismatch+Bulges", " Score"],
        "Mismatches+bulges": [" Mismatches", " Bulges", " Score"],
        "CFD_score": [" Mismatches", " Bulges", " Mismatch+Bulges"],
        "Highest CFD Risk Score": [],
        "Highest CFD Absolute Risk Score": [],
        "CFD_risk_score": [],
        "CFD Absolute Risk Score": [],
    }
    if selected_order is None:
        gi = []
        data = []
    else:  # selected order is not None
        gi = [
            {
                "label": target_label[selected_order][count],
                "value": target_value[selected_order][count],
            }
            for count in range(len(target_value[selected_order]))
        ]
        if selected_order == "Mismatches":
            data = [
                {"label": "0", "value": "0"},
                {"label": "1", "value": "1"},
                {"label": "2", "value": "2"},
                {"label": "3", "value": "3"},
                {"label": "4", "value": "4"},
                {"label": "5", "value": "5"},
                {"label": "6", "value": "6"},
            ]
        elif selected_order == "CFD_score":
            data = [
                {"label": "0.01", "value": "0.01"},
                {"label": "0.1", "value": "0.1"},
                {"label": "0.2", "value": "0.2"},
                {"label": "0.3", "value": "0.3"},
                {"label": "0.4", "value": "0.4"},
                {"label": "0.5", "value": "0.5"},
                {"label": "0.6", "value": "0.6"},
                {"label": "0.7", "value": "0.7"},
                {"label": "0.8", "value": "0.8"},
                {"label": "0.9", "value": "0.9"},
            ]
        elif selected_order == "Mismatches+bulges":
            data = [
                {"label": "0", "value": "0"},
                {"label": "1", "value": "1"},
                {"label": "2", "value": "2"},
                {"label": "3", "value": "3"},
                {"label": "4", "value": "4"},
                {"label": "5", "value": "5"},
                {"label": "6", "value": "6"},
                {"label": "7", "value": "7"},
                {"label": "8", "value": "8"},
            ]
        elif selected_order == "Bulges":
            data = [
                {"label": "0", "value": "0"},
                {"label": "1", "value": "1"},
                {"label": "2", "value": "2"},
            ]
        else:
            gi = []
            data = []
    return gi, data, {"display": "none"}


# drop columns according to threshold
@app.callback(
    Output("maxdrop", "options"),
    [Input("thresh_drop", "value"), Input("order", "value")],
)
def maxdrop(thresh_drop: str, order: str) -> List:
    """Filter the targets table, using the selected threshold value on the
    scores.

    ...

    Parameters
    ----------
    thresh_drop : str
        Threshold value
    order : str
        Ordering criterion

    Returns
    -------
    List
        Filtered data
    """

    if thresh_drop is not None:
        if not isinstance(thresh_drop, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(thresh_drop).__name__}"
            )
    if order is not None:
        if not isinstance(order, str):
            raise TypeError(f"Expected {str.__name__}, got {type(order).__name__}")
    if order == "Mismatches":
        if thresh_drop:
            start_value = int(thresh_drop)
            data = [{"label": str(i), "value": str(i)} for i in range(start_value, 7)]
        else:
            data = []
    elif order == "CFD_score":
        if thresh_drop:
            start_value = int(float(thresh_drop) * 10)
            small_value = False
            if start_value < 1:
                small_value = True
                start_value = 1
            if start_value < 10:
                data = [
                    {"label": f"0.{i}", "value": f"0.{i}"}
                    for i in range(start_value, 10)
                ]
                data.append({"label": "1", "value": "1"})
                if small_value:
                    data.insert(0, {"label": "0.01", "value": "0.01"})
            else:
                data = []
        else:
            data = []
    elif order == "Bulges":
        if thresh_drop:
            start_value = int(thresh_drop)
            data = [{"label": str(i), "value": str(i)} for i in range(start_value, 3)]
        else:
            data = []
    elif order == "Mismatches+bulges":
        if thresh_drop:
            start_value = int(thresh_drop)
            data = [{"label": str(i), "value": str(i)} for i in range(start_value, 9)]
        else:
            data = []
    else:
        data = []
    return data


# reset buttons
@app.callback(
    [
        Output("order", "value"),
        Output("multiorder", "value"),
        Output("maxdrop", "value"),
        Output("thresh_drop", "value "),
        Output("Radio-asc-1", "value"),
    ],
    [Input("reset-val", "n_clicks")],
)
def resetbutton(n_clicks: int) -> Tuple:
    if n_clicks > 0:
        return None, None, None, None, None
    else:
        return None, None, None, None, None
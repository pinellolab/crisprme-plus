"""
This file contains two main functions to deal with querying of the database
in all the specific cases we want.
The two functions are called by the main app file and take in input the parameters 
that are chosen by the user. 
If the parameters are linked to a sort query without shold only the parameters needed are in input and passed to noshold().
Otherwise, we have all the parameters including sholds and the other stuff we need for querying. shold()
The parameters in input are modified with a prefix and a suffix to be coherent with the database structure (column names).
So, if the name of columns in the input file are modified, also these parameters inside here need to be modified.
In output the two functions return a dict that is then included in the table in the main app file.
"""

import dash
import plotly
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash import dash_table
import pandas as pd
import sqlite3

import safe_files  # (PostProcess/: the job's database opened read only, with no link followed)

guide_column = "Spacer+PAM"


def _column(name):
    """``name`` back, to be written between double quotes into a query as a column's:
    one holding a double quote would end the name there and have the rest read as SQL.
    A column's name cannot be bound to a query as a value is."""
    if '"' in str(name):
        raise ValueError(f"Not a column name: {name!r}")
    return name


def _order(asc1):
    """``asc1`` back, if it is the order ASC or DESC, written into a query as it is."""
    if asc1 not in ("ASC", "DESC"):
        raise ValueError(f"Not an order: {asc1!r}")
    return asc1


def shold(
    target,
    n_clicks,
    page_current,
    page_size,
    radio_order,
    orderdrop,
    sholddrop,
    maxdrop,
    asc1,
    url_job,
    guide,
    current_working_directory,
):
    # query with threshold
    # print('entro shold')
    if asc1 == None:
        asc1 = "DESC"

    url = url_job[5:]

    path = current_working_directory + "/Results/" + url + "/." + url + ".db"

    # read only, with no link followed (safe_files): url is the job id of the address, and
    # a link another user put at the job's .db (or a step to it) was opened
    conn = safe_files.sqlite_connect(path, current_working_directory)
    c = conn.cursor()
    param = [
        guide,
        page_size,
        page_current * page_size,
    ]

    if target == "CFD":
        radio_order = str(radio_order) + "_(highest_CFD)"
        orderdrop = str(orderdrop) + "_(highest_CFD)"
    elif target == "fewest":
        radio_order = str(radio_order) + "_(fewest_mm+b)"
        orderdrop = str(orderdrop) + "_(fewest_mm+b)"
    elif target == "CRISPR-Bulge":
        if "CFD" in radio_order:
            radio_order = str(radio_order).replace("CFD", "CRISPR_BULGE")
        try:
            if "CFD" in orderdrop:
                orderdrop = str(orderdrop).replace("CFD", "CRISPR_BULGE")
        except:
            pass
        radio_order = str(radio_order) + "_(highest_CRISPR_BULGE)"
        orderdrop = str(orderdrop) + "_(highest_CRISPR_BULGE)"
    # the columns and the order are written into the query, the values bound to it
    radio_order, orderdrop, asc1 = _column(radio_order), _column(orderdrop), _order(asc1)

    # print(radio_order, orderdrop)

    if orderdrop != "None_(highest_CFD)":  # ordinamento doppio
        if maxdrop == None:
            maxdrop = 1000
        # if maxdrop == None:  # min
        #     df = pd.read_sql_query("SELECT * FROM final_table WHERE \"{}\"=? AND {}>={} ORDER BY \"{}\" {},\"{}\" {} LIMIT ? OFFSET ?".format(
        #         guide_column, radio_order, sholddrop, radio_order, asc1, orderdrop, asc1), conn, params=param)
        # else:  # min e max
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? AND "{}" BETWEEN ? AND ? ORDER BY "{}" {},"{}" {}  LIMIT ? OFFSET ?'.format(
                guide_column,
                radio_order,
                radio_order,
                asc1,
                orderdrop,
                asc1,
            ),
            conn,
            params=[guide, sholddrop, maxdrop] + param[1:],
        )
    else:  # ordinamento singolo
        if maxdrop == None:  # min
            maxdrop = 1000
            # df = pd.read_sql_query("SELECT * FROM final_table WHERE \"{}\"=? AND \"{}\">={} ORDER BY \"{}\" {} LIMIT ? OFFSET ?".format(
            #     guide_column, radio_order, sholddrop, radio_order, asc1), conn, params=param)  # ok
        # else:  # min e max
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? AND "{}" BETWEEN ? AND ? ORDER BY "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, radio_order, asc1
            ),
            conn,
            params=[guide, sholddrop, maxdrop] + param[1:],
        )  # ok

    conn.commit()
    conn.close()
    data = df
    return data


def noshold(
    target,
    n_clicks,
    page_current,
    page_size,
    radio_order,
    orderdrop,
    asc1,
    url_job,
    guide,
    current_working_directory,
):
    # print(asc1)
    # print('entro noshold')
    if asc1 == None:
        asc1 = "DESC"
    url = url_job[5:]
    path = current_working_directory + "/Results/" + url + "/." + url + ".db"

    # read only, with no link followed (safe_files): url is the job id of the address, and
    # a link another user put at the job's .db (or a step to it) was opened
    conn = safe_files.sqlite_connect(path, current_working_directory)
    c = conn.cursor()
    param = [
        guide,
        page_size,
        page_current * page_size,
    ]
    # print("noshold")

    if target == "CFD":
        radio_order = str(radio_order) + "_(highest_CFD)"
        orderdrop = str(orderdrop) + "_(highest_CFD)"
    elif target == "fewest":
        radio_order = str(radio_order) + "_(fewest_mm+b)"
        orderdrop = str(orderdrop) + "_(fewest_mm+b)"
    elif target == "CRISPR-Bulge":
        if "CFD" in radio_order:
            radio_order = str(radio_order).replace("CFD", "CRISPR_BULGE")
        try:
            if "CFD" in orderdrop:
                orderdrop = str(orderdrop).replace("CFD", "CRISPR_BULGE")
        except:
            pass
        radio_order = str(radio_order) + "_(highest_CRISPR_BULGE)"
        orderdrop = str(orderdrop) + "_(highest_CRISPR_BULGE)"
    # the columns and the order are written into the query, the values bound to it
    radio_order, orderdrop, asc1 = _column(radio_order), _column(orderdrop), _order(asc1)

    # print(radio_order, orderdrop)

    if orderdrop != "None_(highest_CFD)":
        # query con multiordinamento
        # print('double ordering')
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? ORDER BY "{}" {}, "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, asc1, orderdrop, asc1
            ),
            conn,
            params=param,
        )
    else:
        # query con ordinamento singolo
        # print('single ordering')
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? ORDER BY "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, asc1
            ),
            conn,
            params=param,
        )

    # query senza soglia
    conn.commit()
    conn.close()
    data = df
    return data

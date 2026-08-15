"""Build the layout of the main webpage of CRISPRme.
The main webpage reads the input data and manages the analysis.
"""

from seq_script import extract_seq, convert_pam
from .pages_utils import (
    ANNOTATIONS_DIR,
    DNA_ALPHABET,
    EMAIL_FILE,
    GENOMES_DIR,
    GITHUB_LINK,
    GUIDES_FILE,
    JOBID_ITERATIONS_MAX,
    JOBID_MAXLEN,
    LOG_FILE,
    PAMS_DIR,
    PAMS_FILE,
    PARAMS_FILE,
    POSTPROCESS_DIR,
    PAPER_LINK,
    QUEUE_FILE,
    RESULTS_DIR,
    SAMPLES_FILE_LIST,
    VALID_CHARS,
    VARIANTS_DATA,
    select_same_len_guides,
    get_available_PAM,
    get_available_CAS,
    index_max_bulges,
    variant_dataset_present,
    has_variant_index,
    get_variant_dataset_options,
    build_active_annotation,
    get_pam_options,
    get_custom_VCF,
    get_available_genomes,
    get_custom_annotations,
    sort_annotation,
    compress_file,
)
from app import (
    URL,
    app,
    operators,
    current_working_directory,
    app_directory,
    DISPLAY_OFFLINE,
    ONLINE,
    pool_executor,
)

from dash.exceptions import PreventUpdate
from dash import Input, Output, State, no_update
from typing import Dict, List, Optional, Tuple, Union
from datetime import datetime

import dash_bootstrap_components as dbc
from dash import html
from dash import dcc

import collections
import subprocess
import filecmp
import random
import string
import os


def _resolve_vcf_folder(genome_selected: str, tok: str) -> str:
    """Resolve the on-disk ``VCFs/`` folder (and enriched-genome suffix) for a
    (genome, dataset-token) pair.

    By convention the folder is ``<genome>_<tok>`` (minus a trailing ``_ref``), and
    that name is returned whenever it exists on disk -- so behaviour is UNCHANGED when
    the genome and its VCF folders are named consistently (the shipped hg38 case, and
    any custom genome named consistently). The fallback makes the web robust to a
    genome directory whose name differs from the VCF folder's genome prefix -- e.g. a
    per-chromosome ``hg38_chr22`` whose variants live in ``hg38_1000G`` -- by picking
    the installed ``VCFs/*_<tok>`` folder whose genome prefix is the longest *token*
    prefix of the selected genome. Returns the convention name if nothing matches, so
    the downstream .list_vcfs / .samplesID / genome_idx stay self-consistent."""
    gpref = genome_selected[:-4] if genome_selected.endswith("_ref") else genome_selected
    conv = f"{gpref}_{tok}"
    vcf_dir = os.path.join(current_working_directory, "VCFs")
    if os.path.isdir(os.path.join(vcf_dir, conv)):
        return conv
    best = None  # (folder_name, prefix)
    if os.path.isdir(vcf_dir):
        suffix = f"_{tok}"
        for d in sorted(os.listdir(vcf_dir)):
            if not os.path.isdir(os.path.join(vcf_dir, d)):
                continue
            if d == tok:
                pfx = ""
            elif d.endswith(suffix):
                pfx = d[: -len(suffix)]
            else:
                continue
            # accept only a whole-token genome prefix of the selected genome
            if pfx == "" or gpref == pfx or gpref.startswith(pfx + "_"):
                if best is None or len(pfx) > len(best[1]):
                    best = (d, pfx)
    return best[0] if best else conv


def _ensure_samplesid(genome_selected: str, vcf_folder: str) -> None:
    """Make sure ``samplesIDs/<vcf_folder>.samplesID.txt`` exists.

    A combined/merged variant dataset (e.g. the batteries-included 1000G+HGDP index,
    folder ``hg38_1000G_HGDP``) needs a single combined sample-ID list, but the data
    repo ships only the per-dataset lists (``hg38_1000G``, ``hg38_HGDP``). Without the
    combined file the search dies at the sample-ID step (submit_job greps a missing
    file -> non-empty log_error -> exit 1). If it's missing, synthesize it by unioning
    the component per-dataset lists (dataset suffix split on ``_``). Best effort: if the
    components aren't all present, leave things unchanged."""
    sdir = os.path.join(current_working_directory, "samplesIDs")
    target = os.path.join(sdir, f"{vcf_folder}.samplesID.txt")
    if os.path.isfile(target):
        return
    gpref = genome_selected[:-4] if genome_selected.endswith("_ref") else genome_selected
    dataset = vcf_folder[len(gpref) + 1:] if vcf_folder.startswith(gpref + "_") else vcf_folder
    if "_" not in dataset:
        return  # a single (non-merged) dataset; nothing to synthesize
    comps = [os.path.join(sdir, f"{gpref}_{c}.samplesID.txt") for c in dataset.split("_")]
    if not comps or not all(os.path.isfile(c) for c in comps):
        return  # cannot synthesize from components
    seen, rows = set(), []
    for c in comps:
        with open(c) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                key = line.split("\t", 1)[0]
                if key not in seen:
                    seen.add(key)
                    rows.append(line.rstrip("\n"))
    try:
        with open(target, "w") as out:
            out.write("\n".join(rows) + "\n")
    except OSError:
        pass


MAX_BULGES = 3  # max allowed bulges
MAX_MMS = 7  # max allowed mismatches
# mismatches, bulges and guides values
AV_MISMATCHES = [{"label": i, "value": i} for i in range(MAX_MMS)]
AV_BULGES = [{"label": i, "value": i} for i in range(MAX_BULGES)]
AV_GUIDE_SEQUENCE = [{"label": i, "value": i} for i in range(15, 26)]
# base editing options
BE_NTS = [{"label": nt, "value": nt} for nt in DNA_ALPHABET]


def split_filter_part(filter_part: str) -> Tuple:
    """Recover filtering operator.

    ...

    Paramters
    --------
    filter_part : str
        Filter

    Returns
    -------
    Tuple
    """

    if not isinstance(filter_part, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_part).__name__}")
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[(name_part.find("{") + 1) : name_part.rfind("}")]
                value_part = value_part.strip()
                v0 = value_part[0]
                if v0 == value_part[-1] and v0 in ("'", '"', "`"):
                    value = value_part[1:-1].replace(f"\\{v0}", v0)
                else:
                    try:
                        value = float(value_part)
                    except:
                        value = value_part
                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value
    return [None] * 3


# load example data
@app.callback(
    [
        Output("text-guides", "value"),
        Output("available-pam", "value"),
        Output("available-genome", "value"),
        # variant-dataset.value is also driven by change_variant_dataset_options
        # (genome-change); load-example is the secondary writer.
        Output("variant-dataset", "value", allow_duplicate=True),
        Output("mms", "value"),
        Output("dna", "value"),
        Output("rna", "value"),
        Output("be-window-start", "value"),
        Output("be-window-stop", "value"),
        Output("be-nts", "value"),
        Output("radio-base_editor", "value"),
        # also emit the window OPTIONS here so the values below actually stick: the
        # be-window dropdowns start with no options, and update_base_editing_dropdown
        # (Input=text-guides) only fills them on the next callback round -- by which
        # point Dash has already dropped a value that matched no option. Setting both
        # atomically fixes the race (duplicate of update_base_editing_dropdown's output).
        Output("be-window-start", "options", allow_duplicate=True),
        Output("be-window-stop", "options", allow_duplicate=True),
        # Set the "Max edits" slider explicitly so the example is driven by the same
        # total-edits cap the report displays -- rather than leaving the slider at its
        # default (which could be tighter than the per-type caps below and make the
        # result title look inconsistent with what the user "chose".)
        Output("max-edits-slider", "value", allow_duplicate=True),
    ],
    [Input("load-example-button", "n_clicks")],
    prevent_initial_call=True,
)
def load_example_data(load_button_click: int) -> List[Union[str, List[str]]]:
    """Load data for CRISPRme example run.

    ...

    Parameters
    ----------
    load_button_click : int
        Click on "Load Example" button

    Returns
    -------
    List
        Example parameters
    """

    # Pick an installed variant panel dynamically (the richest one) so the example never
    # selects a dataset that isn't installed -- e.g. when only the combined 1000G+HGDP
    # panel is present, not a standalone 1000G. Falls back to reference-only.
    example_variant = _preferred_variant(
        [o["value"] for o in get_variant_dataset_options("hg38")]
    )
    example_guide = "CTAACAGTTGCTTTTATCAC"  # 20 nt
    # window options span 1..guide_len (ints), matching update_base_editing_dropdown
    be_window_options = [
        {"label": i, "value": i} for i in range(1, len(example_guide) + 1)
    ]
    return [
        example_guide,  # guide to use
        "20bp-NGG-SpCas9",  # PAM/enzyme to use
        "hg38",  # ref genome to use
        example_variant,  # variant dataset (installed panel, chosen dynamically)
        4,  # MM (int, to match the dropdown option values)
        1,  # DNA bulges
        1,  # RNA bulges
        4,  # start window in base editor (int: the be-window options are ints)
        8,  # stop window in base editor
        "A",  # nt to check in base editor (pre-filled, but base editing stays OFF)
        "N",  # base editing OFF by default -- the example runs a standard search
        be_window_options,  # be-window-start options (set atomically with the value)
        be_window_options,  # be-window-stop options
        5,  # Max edits (mismatches + bulges): the binding total-edits cap (== the default)
    ]


# Job submission and results URL definition
@app.callback(
    [Output("url", "pathname"), Output("url", "search")],
    [Input("submit-job", "n_clicks")],
    [
        State("url", "href"),
        State("available-genome", "value"),
        State("variant-dataset", "value"),
        State("available-pam", "value"),
        State("radio-guide", "value"),
        State("text-guides", "value"),
        State("mms", "value"),
        State("dna", "value"),
        State("rna", "value"),
        State("radio-base_editor", "value"),
        State("be-window-start", "value"),
        State("be-window-stop", "value"),
        State("be-nts", "value"),
        State("checklist-mail", "value"),
        State("example-email", "value"),
        State("job-name", "value"),
        State("max-edits-slider", "value"),
        State("advanced-thresholds-collapse", "is_open"),
    ],
)
def change_url(
    n: int,
    href: str,
    genome_selected: str,
    variant_choice: str,
    pam: str,
    guide_type: str,
    text_guides: List[str],
    mms: int,
    dna: int,
    rna: int,
    radio_be_value: str,
    be_start: int,
    be_stop: int,
    be_nt: str,
    adv_opts: List,
    dest_email: str,
    job_name: str,
    max_edits_val: int,
    advanced_open: bool,
) -> Tuple[str, str]:
    """Launch the targets search and generates the input files for
    post-processing operations, and results visualization.

    It manages the input data given by the user in the main webpage of CRISPRme
    and run the search, notify the user by sending an email when the job is
    completed, and produce the link to the webpage used to visualize the results.

    ** Further details **

    Perform checks on input parameters' consistency.

    To each received job is assigned a different identifier (or job name). This
    allows to easily recognize different job submissions. The job IDs consist
    in alpha-numeric strings of 10 characters (A-z 0-9). The IDs are randomly
    generated. If the generated ID is already assigned to some other job,
    compute another ID. Every 7 iterations, add +1 to ID length (up to 20 chars
    as max length). Once generated the ID, create the job directory. Within the
    job directory, create the `queue.txt` file (for job queueing).

    If the input parameters match those of an already processed search, the
    current job ID is modified to match that of the available analysis (even if
    completed/currently submitted/in queue). Update the email address associated
    to the job and reset the 3 days availability of the results.

    The current policy of CRISPRme allows up to 2 jobs to run concurrently. The
    others are put in queue.

    ...

    Parameters
    ----------
    n : int
        Clicks
    href : str
        URL
    genome_selected : str
        Selected genome
    variant_choice : str
        Selected variant dataset ("ref" / "1000G" / "1000G+HGDP")
    pam : str
        Selected PAM
    guide_type : str
        RNA guide type
    text_guides : str
        Input guides
    mms : int
        Number of mismatches
    dna : int
        Number of DNA bulges
    rna : int
        Number of RNA bulges
    adv_opts : List
        Selected advanced options
    dest_email : str
        User mail address
    job_name : str
        Submitted job ID

    Returns
    -------
    Tuple[str, str]
        URL to retrieve CRISPRme analysis results
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(href, str):
        raise TypeError(f"Expected {str.__name__}, got {type(href).__name__}")
    if genome_selected is not None:
        if not isinstance(genome_selected, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(genome_selected).__name__}"
            )
    # The variant selector is now a single genome-driven dropdown whose value is a
    # scalar token: "ref" (reference only) or a dataset / "+"-joined combo (e.g.
    # "1000G", "1000G+HGDP"). Normalize it to the dataset list the search wiring
    # expects. Custom-VCF ("PV") selection was removed from the form.
    if variant_choice in (None, "", "ref"):
        ref_var = []
    else:
        ref_var = str(variant_choice).split("+")
    vcf_input = None  # personal/custom VCF path removed from the search form
    if pam is not None:
        if not isinstance(pam, str):
            raise TypeError(f"Exepcted {str.__name__}, got {type(pam).__name__}")
    if text_guides is not None:
        if not isinstance(text_guides, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(text_guides).__name__}"
            )
    # if rna is not None:
    #     if not isinstance(rna, int):
    #         raise TypeError(f"Expected {str.__name__}, got {type(rna).__name__}")
    # if dna is not None:
    #     if not isinstance(dna, int):
    #         raise TypeError(f"Expected {str.__name__}, got {type(dna).__name__}")
    if adv_opts is not None:
        if not isinstance(adv_opts, list):
            raise TypeError(f"Expected {list.__name__}, got {type(adv_opts).__name__}")
    if dest_email is not None:
        if not isinstance(dest_email, str):
            raise TypeError(f"Expected {str.__name__}, got {type(dest_email).__name__}")
    if job_name is not None:
        if not isinstance(job_name, str):
            raise TypeError(f"Expected {str.__name__}, got {type(job_name).__name__}")
    if n is None:
        raise PreventUpdate  # do not update the page
    # job start
    print("Launching JOB")
    # ---- Check input. If fails, give simple input
    if (genome_selected is None) or (not genome_selected):
        genome_selected = "hg38_ref"  # use hg38 by default
    if (pam is None) or (not pam):
        pam = "20bp-NGG-SpCas9"  # use Cas9 PAM
        guide_seqlen = 20  # set guide length to 20
    else:
        for c in pam.split("-"):
            if "bp" in c:  # use length specified in PAM
                guide_seqlen = int(c.replace("bp", ""))
    if (text_guides is None) or (not text_guides):
        text_guides = "A" * guide_seqlen
    elif guide_type != "GS":
        text_guides = text_guides.strip()
        if not all(
            [
                len(guide) == len(text_guides.split("\n")[0])
                for guide in text_guides.split("\n")
            ]
        ):
            text_guides = select_same_len_guides(text_guides)
    # remove Ns from guides
    guides_tmp = "\n".join(
        [guide.replace("N", "") for guide in text_guides.split("\n")]
    )
    text_guides = guides_tmp.strip()
    # ---- Generate random job ids
    id_len = 10
    for i in range(JOBID_ITERATIONS_MAX):
        # get already assigned job ids
        assigned_ids = [
            d
            for d in os.listdir(os.path.join(current_working_directory, RESULTS_DIR))
            if (
                os.path.isdir(os.path.join(current_working_directory, RESULTS_DIR, d))
                and not d.startswith(".")  # avoid hidden files/directories
            )
        ]
        job_id = "".join(
            random.choices(string.ascii_uppercase + string.digits, k=id_len)
        )
        if job_id not in assigned_ids:  # suitable job id
            break
        if i > 7:
            i = 0  # restart
            id_len += 1  # increase ID length
            if id_len > JOBID_MAXLEN:  # reached maximum length
                break
    if job_name and job_name != "None":
        assert isinstance(job_name, str)
        job_id = f"{job_name}_{job_id}"
    result_dir = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    # create results directory
    cmd = f"mkdir {result_dir}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f"An error occurred while running {cmd}")
    # NOTE test command for queue
    cmd = f"touch {os.path.join(current_working_directory, RESULTS_DIR, job_id, QUEUE_FILE)}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f"An error occurred while running {cmd}")
    # ---- Set search parameters
    # ANNOTATION: there is no per-search annotation selector anymore. The search
    # applies the annotations ENABLED for this genome (managed in Settings ->
    # Annotations); build_active_annotation assembles them into a single annotation
    # bed (the built-in ENCODE cCREs + DHS + GENCODE bundle is enabled by default,
    # multiple enabled beds are merged, none -> "vuoto.txt"). gencode is the built-in
    # bundle's gene-annotation companion and rides with its enabled state.
    annotation_name, gencode_name = build_active_annotation(genome_selected)
    # GENOME TYPE CHECK
    ref_comparison = False
    genome_type = "ref"  # search is 'ref' or 'both'
    if len(ref_var) > 0:
        ref_comparison = True
        genome_type = "both"
    search_index = True
    genome_selected = genome_selected.replace(" ", "_")
    genome_ref = genome_selected
    # NOTE indexed genomes names format:
    # PAM + _ + bMax + _ + genome_selected
    # VCF CHECK
    # TODO: check here
    if genome_type == "ref":
        sample_list = None
    sample_list = []
    try:
        with open(os.path.join(result_dir, ".list_vcfs.txt"), mode="w") as handle_vcf:
            if not ref_var:
                vcf_folder = "_"
                handle_vcf.write(f"{vcf_folder}\n")
            # One dropdown value == one dataset == one enriched-genome/VCF folder, whose
            # name is derived DYNAMICALLY as "<genome>_<dataset>" (no hardcoded 1000G/
            # HGDP or hg38 literals). This works for any current or future genome/dataset
            # (e.g. a pig susScr11 + a custom VCF) and for a merged panel whose dropdown
            # token is itself "1000G_HGDP" -> folder "hg38_1000G_HGDP". The legacy
            # free-text custom-VCF box (VARIANTS_DATA[2]) still maps to the typed name.
            for _tok in ref_var:
                if not _tok or _tok == "ref":
                    continue
                if _tok == VARIANTS_DATA[2] and vcf_input:  # free-text custom VCF
                    vcf_folder = vcf_input
                else:
                    # resolve the actual on-disk VCFs/ folder (robust to a genome dir
                    # whose name differs from the VCF folder's genome prefix)
                    vcf_folder = _resolve_vcf_folder(genome_selected, _tok)
                # a batteries-included combined index ships without a combined sample-ID
                # list; synthesize it from the per-dataset lists so the search doesn't die
                _ensure_samplesid(genome_selected, vcf_folder)
                sample_list.append(f"{vcf_folder}.samplesID.txt")
                handle_vcf.write(f"{vcf_folder}\n")
    except OSError as e:
        raise e
    try:
        with open(
            os.path.join(result_dir, ".samplesID.txt"), mode="w"
        ) as handle_samples:
            for e in sample_list:
                handle_samples.write(f"{e}\n")
    except OSError as e:
        raise e
    # manage email sending
    send_email = False
    if adv_opts is None:
        adv_opts = []
    if "email" in adv_opts and check_mail_address(dest_email):
        send_email = True
        try:
            with open(os.path.join(result_dir, EMAIL_FILE), mode="w") as handle_mail:
                handle_mail.write(f"{dest_email}\n")
                handle_mail.write(
                    f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                )
                handle_mail.write(
                    f"{datetime.utcnow().strftime('%m/%d/%Y, %H:%M:%S')}\n"
                )
        except OSError as e:
            raise e
    else:
        dest_email = "_"  # null value
    # manage PAM
    pam_len = 0
    pam_begin = False
    try:
        with open(
            os.path.join(current_working_directory, PAMS_DIR, f"{pam}.txt")
        ) as handle_pam:
            pam_char = handle_pam.readline()
            index_pam_value = pam_char.split()[-1]
            if int(index_pam_value) < 0:
                end_idx = -int(index_pam_value)
                pam_char = pam_char.split()[0][:end_idx]
                pam_begin = True
            else:
                end_idx = int(index_pam_value)
                pam_char = pam_char.split()[0][-end_idx:]
            pam_len = end_idx
    except OSError as e:
        raise e
    # manage guide type
    if guide_type == "GS":
        # text_sequence = text_guides
        # Extract sequence and create the guides
        guides = list()
        for seqname_and_seq in text_guides.split(">"):
            if not seqname_and_seq:
                continue
            seqname = seqname_and_seq[: seqname_and_seq.find("\n")]
            seq = seqname_and_seq[seqname_and_seq.find("\n") :]
            seq = seq.strip()  # remove endline
            if "chr" in seq:
                for line in seq.split("\n"):
                    if not line:
                        continue
                    line_split = line.strip().split()
                    # line_split = re.split(r";|,|.|:|-| ", line.strip())
                    # print(line_split)
                    seq_read = f"{line_split[0]}:{line_split[1]}-{line_split[2]}"
                    assert bool(seqname)
                    assert bool(seq_read)
                    assert bool(genome_ref)
                    seq_read = extract_seq.extractSequence(
                        seqname, seq_read, genome_ref.replace(" ", "_")
                    )
            else:
                seq_read = "".join(seq.split()).strip()
            guides.extend(
                convert_pam.getGuides(seq_read, pam_char, guide_seqlen, pam_begin)
            )
        guides = list(set(guides))  # remove potential duplicate guides
        # create new guides dataset
        if not guides:
            guides = "A" * guide_seqlen
        text_guides = "\n".join(guides).strip()
        assert bool(guides)
    # force guides to be upper case
    text_guides = text_guides.upper()
    text_guides_tmp = [
        guide for guide in text_guides.split("\n") if len(guide) == guide_seqlen
    ]
    if not text_guides_tmp:  # no suitable guide found
        text_guides_tmp.append("A" * guide_seqlen)
    text_guides = "\n".join(text_guides_tmp)
    for guide in text_guides.split("\n"):
        for nt in guide:
            if nt not in VALID_CHARS:
                # remove forbidden characters from guide
                text_guides = text_guides.replace(nt, "")
    # set limit to 100 guides per run in the website
    if len(text_guides.split("\n")) > 100:
        text_guides = "\n".join(text_guides.split("\n")[:100]).strip()
    # Adjust guides by adding Ns (compatible with Crispritz)
    if pam_begin:
        pam_to_file = pam_char + ("N" * guide_seqlen) + " " + index_pam_value
        pam_to_indexing = pam_char + ("N" * 25) + " " + index_pam_value
    else:
        pam_to_file = ("N" * guide_seqlen) + pam_char + " " + index_pam_value
        pam_to_indexing = ("N" * 25) + pam_char + " " + index_pam_value
    # store PAMs file
    try:
        with open(os.path.join(result_dir, PAMS_FILE), mode="w") as handle_pam:
            handle_pam.write(pam_to_file)
    except OSError as e:
        raise e
    pams_file = os.path.join(result_dir, PAMS_FILE)
    guides_file = os.path.join(result_dir, GUIDES_FILE)
    if text_guides:
        try:
            with open(guides_file, mode="w") as handle_guides:
                if pam_begin:
                    text_guides = "N" * pam_len + text_guides.replace(
                        "\n", "\n" + "N" * pam_len
                    )
                else:
                    text_guides = (
                        text_guides.replace("\n", "N" * pam_len + "\n") + "N" * pam_len
                    )
                handle_guides.write(text_guides)
        except OSError as e:
            raise e
    # bulges
    dna = int(dna)
    rna = int(rna)
    # Index name budget: mirror the validated CLI (crisprme.py complete-search) EXACTLY.
    # There, bMax = max(bDNA, bRNA) and the TST index is named/built as
    # "<pam>_<bMax+1>_<genome>" (the +1 is for alignments starting with a gap). So the
    # folder budget N = max(dna,rna)+1 and the index supports up to N-1 = max(dna,rna)
    # bulges of each type. Using dna+rna here would mislabel the index and diverge from
    # the CLI/shell, which resolve by max(dna,rna)+1.
    max_bulges = max(dna, rna) + 1
    assert isinstance(dna, int)
    assert isinstance(rna, int)
    # base editing
    if be_start is None or not bool(be_start) or radio_be_value == "N":
        be_start = 1
    else:
        be_start = int(be_start)
    if be_stop is None or not bool(be_stop) or radio_be_value == "N":
        be_stop = 0
    else:
        be_stop = int(be_stop)
    if be_nt is None or not bool(be_nt) or radio_be_value == "N":
        be_nt = "none"
    else:
        assert be_nt in DNA_ALPHABET
    assert isinstance(be_start, int)
    assert isinstance(be_stop, int)
    assert isinstance(be_nt, str)
    if search_index:
        search = False
    # Check if index exists, otherwise set generate_index to true
    genome_idx_list = []
    if genome_type == "ref":
        genome_idx_list.append(f"{pam_char}_{max_bulges}_{genome_selected}")
    else:
        # Variant index name == "<pam>_<budget>_<genome>+<enriched>", where the enriched
        # genome name is derived DYNAMICALLY as "<genome>_<dataset>" (no hardcoded 1000G/
        # HGDP or hg38 literals) -- mirroring the VCF-folder mapping above so it matches
        # the built index (e.g. "NNN_3_hg38+hg38_1000G_HGDP") for any genome/dataset.
        # Display-only: the shell resolves the actual precomputed index by scanning for a
        # sufficient bulge budget, so this string need not match the on-disk budget.
        for _tok in ref_var:
            if not _tok or _tok == "ref":
                continue
            if _tok == VARIANTS_DATA[2] and vcf_input:  # free-text custom VCF
                _enriched = vcf_input
            else:
                # same resolver as the VCF-folder block above, so the index name
                # matches the enriched genome/index actually on disk
                _enriched = _resolve_vcf_folder(genome_selected, _tok)
            genome_idx_list.append(f"{pam_char}_{max_bulges}_{genome_selected}+{_enriched}")
    genome_idx = ",".join(genome_idx_list)
    # Total-edits cap (submit_job arg 26 + recorded in .Params.txt for the report).
    # This is the SINGLE source of truth: the same value is written to Params below and
    # passed to the search below. Simple mode: the "Max edits" slider governs (it is the
    # constraint the user actually set). Advanced ("old") mode: the per-type mm/bulge caps
    # govern, so the total cap is their sum (== unbounded, 2.1.x behavior). Floored at 1
    # (a cap of 0 hits a crispritz --max-edits 0 empty-set bug; 1 still keeps the 0-edit
    # on-target).
    if advanced_open:
        max_total_edits = int(mms) + int(dna) + int(rna)
    else:
        max_total_edits = int(max_edits_val) if max_edits_val is not None else 5
    max_total_edits = max(1, max_total_edits)
    # Create .Params.txt file
    try:
        with open(os.path.join(result_dir, PARAMS_FILE), mode="w") as handle_params:
            handle_params.write(f"Genome_selected\t{genome_selected}\n")
            handle_params.write(f"Genome_ref\t{genome_ref}\n")
            if search_index:
                handle_params.write(f"Genome_idx\t{genome_idx}\n")
            else:
                handle_params.write(f"Genome_idx\tNone\n")
            handle_params.write(f"Pam\t{pam_char}\n")
            handle_params.write(f"Max_bulges\t{max_bulges}\n")
            handle_params.write(f"Mismatches\t{mms}\n")
            handle_params.write(f"DNA\t{dna}\n")
            handle_params.write(f"RNA\t{rna}\n")
            # The actual binding constraint for the search: total mismatches + bulges.
            # In simple mode this is the "Max edits" slider value the user set (which can
            # be tighter than the per-type caps above); recording it lets the report show
            # the real cap rather than only the looser per-type numbers.
            handle_params.write(f"Max_total_edits\t{max_total_edits}\n")
            # Which threshold control the user used: 'advanced' (explicit per-type
            # mm/DNA/RNA caps) or 'simple' (the single "Max edits" slider). The result
            # title uses this to avoid showing per-type caps that contradict a tighter
            # total cap in simple mode.
            handle_params.write(
                f"Threshold_mode\t{'advanced' if advanced_open else 'simple'}\n"
            )
            handle_params.write(f"Annotation\t{annotation_name}\n")
            # nuclease is derived from the PAM token (<len>bp-<motif>-<enzyme>),
            # since the separate Cas-protein selector was removed
            _pam_parts = str(pam).split("-")
            nuclease = "-".join(_pam_parts[2:]) if len(_pam_parts) >= 3 else str(pam)
            handle_params.write(f"Nuclease\t{nuclease}\n")
            handle_params.write(f"Ref_comp\t{ref_comparison}\n")
            handle_params.write(f"BE_nucleotide\t{be_nt}\n")
            handle_params.write(f"BE_start\t{be_start}\n")
            handle_params.write(f"BE_stop\t{be_stop}\n")
    except OSError as e:
        raise e
    # ---- Check if input parameters (mms, bulges, pam, guides, genome) match
    # those of previous searches
    computed_results_dirs = [
        d
        for d in os.listdir(os.path.join(current_working_directory, RESULTS_DIR))
        if (
            os.path.isdir(os.path.join(current_working_directory, RESULTS_DIR, d))
            and not d.startswith(".")  # ignore hidden directories
        )
    ]
    computed_results_dirs.remove(job_id)  # remove current job results
    for res_dir in computed_results_dirs:
        if os.path.exists(
            os.path.join(current_working_directory, RESULTS_DIR, res_dir, PARAMS_FILE)
        ):
            if filecmp.cmp(
                os.path.join(
                    current_working_directory, RESULTS_DIR, res_dir, PARAMS_FILE
                ),
                os.path.join(result_dir, PARAMS_FILE),
            ):
                try:
                    # old job guides
                    guides_old = (
                        open(
                            os.path.join(
                                current_working_directory,
                                RESULTS_DIR,
                                res_dir,
                                GUIDES_FILE,
                            )
                        )
                        .read()
                        .split("\n")
                    )
                    # current job guides
                    guides_current = (
                        open(
                            os.path.join(
                                current_working_directory,
                                RESULTS_DIR,
                                job_id,
                                GUIDES_FILE,
                            )
                        )
                        .read()
                        .split("\n")
                    )
                except OSError as e:
                    raise e
                if collections.Counter(guides_old) == collections.Counter(
                    guides_current
                ):
                    if os.path.exists(
                        os.path.join(
                            current_working_directory, RESULTS_DIR, res_dir, LOG_FILE
                        )
                    ):  # log file found
                        adj_date = False
                        try:
                            with open(
                                os.path.join(
                                    current_working_directory,
                                    RESULTS_DIR,
                                    res_dir,
                                    LOG_FILE,
                                )
                            ) as handle_log:
                                log_data = handle_log.read().strip()
                                if "Job\tDone" in log_data:
                                    adj_date = True
                                    log_data = log_data.split("\n")
                                    date_new = subprocess.Popen(
                                        ["echo $(date)"],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        shell=True,
                                    )
                                    out, err = date_new.communicate()
                                    log_to_write = "\n".join(log_data[:-1])
                                    date_write = str(
                                        f"{log_to_write}\nJob\tDone\t"
                                        f"{out.decode('UTF-8').strip()}"
                                    )
                        except OSError as e:
                            raise e
                        # Only dedup onto a prior run that FINISHED CLEANLY. A
                        # failed/incomplete prior result (no "Job\tDone", or a
                        # non-empty log_error.txt -- e.g. left behind by a bug that
                        # has since been fixed) must NOT be resurfaced: reusing it
                        # shows "The selected result encountered some errors, please
                        # remove it and try to submit again". Skip it so THIS
                        # submission runs fresh instead of inheriting the old error.
                        _prev_err = os.path.join(
                            current_working_directory, RESULTS_DIR, res_dir, "log_error.txt"
                        )
                        _prev_failed = (
                            os.path.exists(_prev_err) and os.path.getsize(_prev_err) > 0
                        )
                        if (not adj_date) or _prev_failed:
                            continue  # unhealthy prior result -> do not dedup onto it
                        if adj_date:
                            try:
                                with open(
                                    os.path.join(
                                        current_working_directory,
                                        RESULTS_DIR,
                                        res_dir,
                                        LOG_FILE,
                                    ),
                                    mode="w+",
                                ) as handle_log:
                                    assert date_write
                                    handle_log.write(date_write)
                            except OSError as e:
                                raise e
                            if send_email:
                                # Send mail with file in job_id dir with link to
                                # job already done, note that job_id directory
                                # will be deleted
                                try:
                                    with open(
                                        os.path.join(
                                            current_working_directory,
                                            RESULTS_DIR,
                                            res_dir,
                                            EMAIL_FILE,
                                        ),
                                        mode="w+",
                                    ) as handle_email:
                                        handle_email.write(f"{dest_email}\n")
                                        handle_email.write(
                                            f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                                        )
                                        handle_email.write(
                                            "".join(
                                                [
                                                    datetime.utcnow().strftime(
                                                        "%m/%d/%Y, %H:%M:%S"
                                                    ),
                                                    "\n",
                                                ]
                                            )
                                        )
                                except OSError as e:
                                    raise e
                        elif send_email:
                            # Job is not finished yet. Add current user's email
                            # to email.txt
                            if os.path.exists(
                                os.path.join(
                                    current_working_directory,
                                    RESULTS_DIR,
                                    res_dir,
                                    EMAIL_FILE,
                                )
                            ):
                                try:
                                    with open(
                                        os.path.join(
                                            current_working_directory,
                                            RESULTS_DIR,
                                            res_dir,
                                            EMAIL_FILE,
                                        ),
                                        mode="a+",
                                    ) as handle_email:
                                        handle_email.write("--OTHEREMAIL--")
                                        handle_email.write(f"{dest_email}\n")
                                        handle_email.write(
                                            f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                                        )
                                        handle_email.write(
                                            "".join(
                                                [
                                                    datetime.utcnow().strftime(
                                                        "%m/%d/%Y, %H:%M:%S"
                                                    ),
                                                    "\n",
                                                ]
                                            )
                                        )
                                except OSError as e:
                                    raise e
                            else:
                                try:
                                    with open(
                                        os.path.join(
                                            current_working_directory,
                                            RESULTS_DIR,
                                            res_dir,
                                            EMAIL_FILE,
                                        ),
                                        mode="w+",
                                    ) as handle_email:
                                        handle_email.write(f"{dest_email}\n")
                                        handle_email.write(
                                            f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                                        )
                                        handle_email.write(
                                            "".join(
                                                [
                                                    datetime.utcnow().strftime(
                                                        "%m/%d/%Y, %H:%M:%S"
                                                    ),
                                                    "\n",
                                                ]
                                            )
                                        )
                                except OSError as e:
                                    raise e
                        current_job_dir = os.path.join(
                            current_working_directory, RESULTS_DIR, job_id
                        )
                        cmd = f"rm -r {current_job_dir}"
                        code = subprocess.call(cmd, shell=True)
                        if code != 0:
                            raise ValueError(f"An error occurred while running {cmd}")
                        return "/load", f"?job={res_dir}"
                    else:
                        # log file not found
                        # we may have entered a job directory that was in queue
                        if os.path.exists(
                            os.path.join(
                                current_working_directory,
                                RESULTS_DIR,
                                res_dir,
                                QUEUE_FILE,
                            )
                        ):
                            if send_email:
                                if os.path.exists(
                                    os.path.join(
                                        current_working_directory,
                                        RESULTS_DIR,
                                        res_dir,
                                        EMAIL_FILE,
                                    )
                                ):
                                    try:
                                        with open(
                                            os.path.join(
                                                current_working_directory,
                                                RESULTS_DIR,
                                                res_dir,
                                                EMAIL_FILE,
                                            ),
                                            mode="a+",
                                        ) as handle_email:
                                            handle_email.write("--OTHEREMAIL--")
                                            handle_email.write(f"{dest_email}\n")
                                            handle_email.write(
                                                f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                                            )
                                            handle_email.write(
                                                "".join(
                                                    [
                                                        datetime.utcnow().strftime(
                                                            "%m/%d/%Y, %H:%M:%S"
                                                        ),
                                                        "\n",
                                                    ]
                                                )
                                            )
                                    except OSError as e:
                                        raise e
                                else:
                                    try:
                                        with open(
                                            os.path.join(
                                                current_working_directory,
                                                RESULTS_DIR,
                                                res_dir,
                                                EMAIL_FILE,
                                            ),
                                            mode="w+",
                                        ) as handle_email:
                                            handle_email.write(f"{dest_email}\n")
                                            handle_email.write(
                                                f"{''.join(href.split('/')[:-1])}/load?job={job_id}\n"
                                            )
                                            handle_email.write(
                                                "".join(
                                                    [
                                                        datetime.utcnow().strftime(
                                                            "%m/%d/%Y, %H:%M:%S"
                                                        ),
                                                        "\n",
                                                    ]
                                                )
                                            )
                                    except OSError as e:
                                        raise e
                            return ("/load", f"?job={res_dir}")
    # merge default is 3 nt wide
    merge_default = 3
    print(
        str(
            f"Submitted JOB {job_id}. STDOUT > log_verbose.txt; STDERR > "
            "log_error.txt"
        )
    )
    # set sorting criteria for score and fewest
    sorting_criteria_scoring = "mm+bulges"
    sorting_criteria = "mm+bulges,mm"
    # TODO: use functions rather than calling scripts
    run_job_sh = os.path.join(
        app_directory, POSTPROCESS_DIR, "submit_job_automated_new_multiple_vcfs.sh"
    )
    genome = os.path.join(current_working_directory, GENOMES_DIR, genome_ref)
    vcfs = os.path.join(result_dir, ".list_vcfs.txt")
    annotation = os.path.join(
        current_working_directory, ANNOTATIONS_DIR, annotation_name
    )
    pam_file = os.path.join(current_working_directory, PAMS_DIR, f"{pam}.txt")
    samples_ids = os.path.join(result_dir, SAMPLES_FILE_LIST)
    postprocess = os.path.join(app_directory, POSTPROCESS_DIR)
    gencode = os.path.join(current_working_directory, ANNOTATIONS_DIR, gencode_name)
    log_verbose = os.path.join(result_dir, "log_verbose.txt")
    log_error = os.path.join(result_dir, "log_error.txt")
    assert isinstance(dna, int)
    assert isinstance(rna, int)

    # if annotation requested, compress and index bed 
    if annotation_name != "vuoto.txt":
        annotation = sort_annotation(annotation)
    else:
        if not os.path.isfile(annotation):
            code = subprocess.call(f"touch {annotation}", shell=True)

    if gencode_name != "vuoto.txt":
        gencode = compress_file(gencode)
    else:
        if not os.path.isfile(gencode):
            code = subprocess.call(f"touch {gencode}", shell=True)

    # max_total_edits (submit_job arg 26) was computed above as the single source of
    # truth and written to .Params.txt; reuse it here so Params and the actual search
    # can never disagree.
    # args 23-25 keep submit_job's defaults (cicd_test, vcf-filter-pass-values,
    # index_path) so that arg 26 (max_total_edits) lands in the right position.
    cmd = f"{run_job_sh} {genome} {vcfs} {guides_file} {pam_file} {annotation} {samples_ids} {max_bulges} {mms} {dna} {rna} {merge_default} {result_dir} {postprocess} {4} {current_working_directory} {gencode} {dest_email} {be_start} {be_stop} {be_nt} {sorting_criteria_scoring} {sorting_criteria} False PASS,. _ {max_total_edits} 1> {log_verbose} 2>{log_error}"
    # run job
    pool_executor.submit(subprocess.run, cmd, shell=True)
    return ("/load", f"?job={job_id}")


# Check input presence
@app.callback(
    [
        Output("submit-job", "n_clicks"),
        Output("modal", "is_open"),
        Output("available-genome", "className"),
        Output("available-pam", "className"),
        Output("text-guides", "style"),
        Output("mms", "className"),
        Output("dna", "className"),
        Output("rna", "className"),
        Output("warning-list", "children"),
    ],
    [Input("check-job", "n_clicks"), Input("close", "n_clicks")],
    [
        State("available-genome", "value"),
        State("available-pam", "value"),
        State("radio-guide", "value"),
        State("text-guides", "value"),
        State("mms", "value"),
        State("dna", "value"),
        State("rna", "value"),
        State("variant-dataset", "value"),
        State("modal", "is_open"),
    ],
)
# len_guide_seq, active_tab ,
def check_input(
    n: int,
    n_close: int,
    genome_selected: str,
    pam: str,
    guide_type: str,
    text_guides: List[str],
    mms: int,
    dna: int,
    rna: int,
    variant_choice: str,
    is_open: bool,
) -> Tuple:
    """Check the correctness of input data and fields. If the input data are
    missing or wrong the borders of the corresponding box are colored in red.
    If input are missing, a Modal element is displayed listing the missing
    elements. The callback is triggered when clicking on the "Submit" button or
    when the Modal object is closed ("Close" button or clicking on-screen when
    the Modal object is open).

    ...

    Parameters
    ----------
    n : int
        Clicks
    n_close : int
        Clicks
    genome_selected : str
        Selected genome
    pam : str
        PAM
    guide_type : str
        Guide type
    text_guides : List[str]
        List of selected guides
    mms : str
        Number of mismatches
    dna : str
        Number of DNA bulges
    rna : str
        Number of RNA bulges
    is_open : bool
        True if Modal object is open

    Returns
    -------
    Tuple
        Input data used during CRISPRme analysis
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if is_open is not None:
        if not isinstance(is_open, bool):
            raise TypeError(f"Expected {bool.__name__}, got {type(is_open).__name__}")
    print("Check input for JOB")
    if n is None:
        raise PreventUpdate  # do not check data --> no trigger
    if is_open is None:
        is_open = False
    classname_red = "missing-input"
    genome_update = None
    pam_update = None
    text_update = {"width": "300px", "height": "30px"}
    mms_update = None
    dna_update = None
    rna_update = None
    be_start_update = None
    be_stop_update = None
    update_style = False
    miss_input_list = []  # recover missing inputs
    # display missing genome
    if genome_selected is None or not bool(genome_selected):
        genome_update = classname_red
        update_style = True
        miss_input_list.append("Genome")
    if genome_selected is None or not bool(genome_selected):
        genome_selected = "hg38_ref"
    genome_ref = genome_selected
    if pam is None or not bool(pam):
        pam_update = classname_red
        update_style = True
        miss_input_list.append("PAM")
    if mms is None:
        mms_update = classname_red
        update_style = True
        miss_input_list.append("Allowed Mismatches")
    if dna is None:
        dna_update = classname_red
        update_style = True
        miss_input_list.append("Bulge DNA size")
    if rna is None:
        rna_update = classname_red
        update_style = True
        miss_input_list.append("Bulge RNA size")
    if pam is None or not bool(pam):
        pam = "20bp-NGG-SpCas9"
        len_guide_sequence = 20
    else:
        for e in pam.split("-"):
            if "bp" in e:
                len_guide_sequence = int(e.replace("bp", ""))
    no_guides = False
    if text_guides is None or not bool(text_guides):
        text_guides = "A" * len_guide_sequence
        no_guides = True
    elif guide_type != "GS":
        text_guides = text_guides.strip()
        if not all(
            [len(g) == len(text_guides.split("\n")[0]) for g in text_guides.split("\n")]
        ):
            text_guides = select_same_len_guides(text_guides)
    # check PAM
    try:
        with open(
            os.path.join(current_working_directory, PAMS_DIR, f"{pam}.txt")
        ) as handle_pam:
            pam_char = handle_pam.readline()
            index_pam_value = int(pam_char.split()[-1])
            if index_pam_value < 0:
                end_idx = index_pam_value * (-1)
                pam_char = pam_char.split()[0][:end_idx]
                pam_begin = True
            else:
                end_idx = index_pam_value
                pam_char = pam_char.split()[0][(end_idx * (-1)) :]
                pam_begin = False
    except OSError as e:
        raise e
    if guide_type == "GS":
        # Extract sequence and create the guides
        guides = []
        for seqname_and_seq in text_guides.split(">"):
            if not seqname_and_seq:
                continue
            seqname = seqname_and_seq[: seqname_and_seq.find("\n")]
            seq = seqname_and_seq[seqname_and_seq.find("\n") :]
            seq = seq.strip()
            if "chr" in seq:
                for line in seq.split("\n"):
                    if not line.strip():
                        continue
                    line_split = line.strip().split()
                    # check suitable BED-like input
                    if len(line_split) < 3:  # chr start stop (minimal)
                        miss_input_list.append(
                            str(
                                "Wrong guides BED coordinates read. Please input "
                                "genomic coordinates as 'chr    start   stop'"
                            )
                        )
                        guides = []  # reset guides
                        break
                    if not line_split[1].isdigit():
                        miss_input_list.append(
                            str("The start coordinate must contain only digits")
                        )
                        guides = []  # reset guides
                        break
                    if not line_split[2].isdigit():
                        miss_input_list.append(
                            str("The stop coordinate must contain only digits")
                        )
                        guides = []  # reset guides
                        break
                    if int(line_split[1]) > int(line_split[2]):
                        miss_input_list.append(
                            str(
                                "Wrong genomic coordinates. The stop coordinate "
                                "seems larger than the start coordinate."
                            )
                        )
                        guides = []  # reset guides
                        break
                    # line_split = re.split(r";|,|.|:|-| ", line.strip())
                    # print(line_split)
                    seq_read = f"{line_split[0]}:{line_split[1]}-{line_split[2]}"
                    assert bool(seqname)
                    assert bool(seq_read)
                    assert bool(genome_ref)
                    seq_read = extract_seq.extractSequence(
                        seqname, seq_read, genome_ref.replace(" ", "_")
                    )
                    guides.extend(
                        convert_pam.getGuides(
                            seq_read, pam_char, len_guide_sequence, pam_begin
                        )
                    )
            else:
                seq_read = "".join(seq.split()).strip()
                guides.extend(
                    convert_pam.getGuides(
                        seq_read, pam_char, len_guide_sequence, pam_begin
                    )
                )
        guides = list(set(guides))  # remove potential duplicates
        if not guides:
            guides = "A" * len_guide_sequence
            no_guides = True
        text_guides = "\n".join(guides).strip()
    text_guides = text_guides.upper()
    text_guides_tmp = [
        guide.replace("N", "")
        for guide in text_guides.split("\n")
        if len(guide.replace("N", "")) == len_guide_sequence
    ]
    if not text_guides_tmp:  # no guide found
        text_guides_tmp.append("A" * len_guide_sequence)
        no_guides = True
    text_guides = "\n".join(text_guides_tmp)
    # remove forbidden characters from guides
    for guide in text_guides.split("\n"):
        for nt in guide:
            if nt not in VALID_CHARS:
                text_guides = text_guides.replace(nt, "")
    # set limit to 1000000000 guides per run
    if len(text_guides.split("\n")) > 1000000000:
        text_guides = "\n".join(text_guides.split("\n")[:1000000000]).strip()
    if no_guides:
        text_update = {"width": "300px", "height": "30px", "border": "1px solid red"}
        update_style = True
        miss_input_list.append(
            str(
                "Input at least one correct guide, correct guides must have the "
                "length requested for the selected PAM sequence (e.g., 20bp, "
                "21bp, etc)"
            )
        )
    # WEB-ONLY guard: never build an index on the fly. Block the web submit when no
    # installed index covers the requested bulge depth. This only withholds the browser
    # launch (update_style=True -> submit-job returns None below, so change_url never
    # fires); the CLI complete-search build path in submit_job is deliberately intact.
    # The bulge depth an index must support is max(dna, rna) -- exactly matching the CLI
    # (bMax = max(bDNA, bRNA); index folder N = bMax+1) and the shell scan (N >= bMax+1).
    # index_max_bulges returns N-1, so "index_max_bulges >= max(dna, rna)" is precisely
    # the shell's index-availability test (incl. NNN-pamless + combined-dataset matching).
    if genome_selected and pam and dna is not None and rna is not None:
        need = max(int(dna), int(rna))
        _ok = index_max_bulges(genome_selected, pam, None) >= need
        _sel = (
            []
            if variant_choice in (None, "", "ref")
            else [v for v in str(variant_choice).split("+") if v]
        )
        for _d in _sel:
            _ok = _ok and (index_max_bulges(genome_selected, pam, _d) >= need)
        if not _ok:
            update_style = True
            miss_input_list.append(
                "No installed index supports %d bulge(s) (DNA %s / RNA %s) for this "
                "genome / PAM / variant selection. The web app never builds an index on "
                "the fly — reduce the DNA/RNA bulges, or download/build a matching index "
                "first (Settings → Data manager, or 'crisprme.py download --what index' "
                "/ 'crisprme.py build-index-only')."
                % (need, dna, rna)
            )
    miss_input = html.Div(
        [
            html.P("The following inputs are wrong or missing:"),
            html.Ul([html.Li(x) for x in miss_input_list]),
            html.P("Please fill in the values before submitting the job"),
        ]
    )
    if not update_style:
        print("All input read correctly")
        return (
            1,
            False,
            genome_update,
            pam_update,
            text_update,
            mms_update,
            dna_update,
            rna_update,
            miss_input,
        )
    return (
        None,
        (not is_open),
        genome_update,
        pam_update,
        text_update,
        mms_update,
        dna_update,
        rna_update,
        miss_input,
    )


@app.callback(
    Output("fade-len-guide", "is_in"),
    [Input("tabs", "active_tab")],
    [State("fade-len-guide", "is_in")],
)
def reset_tab(current_tab: str, is_in: bool) -> bool:
    """Manages the fading of the dropdown bar for the guide length, when the tab
    'Sequence' is active.

    ...

    Parameters
    ----------
    current_tab : str
        Current active tab
    is_in : bool
        True if dropdown's guide length is displayed

    Returns
    -------
    bool
    """

    if current_tab is not None:
        if not isinstance(current_tab, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(current_tab).__name__}"
            )
    if current_tab is None:
        raise PreventUpdate  # do not do anything
    if current_tab == "guide-tab":
        return False
    return True


# Check if email address is valid
@app.callback(Output("example-email", "style"), [Input("example-email", "value")])
def is_email_valid(email: str) -> Dict[str, str]:
    """Check if the provided mail address is valid.
    Change the mail box borders to green or red, accordingly.

    ...

    Parameters
    ----------
    email: str
        Email address

    Returns
    -------
    Dict[str, str]
        Email box borders color
    """
    if email is not None:
        if not isinstance(email, str):
            raise TypeError(f"Expected {str.__name__}, got {type(email).__name__}")
    if email is None:
        raise PreventUpdate  # do not do anything
    if ("@" in email) and (len(email.split("@")) == 2):
        # mail address should be valid
        return {"border": "1px solid #94f033", "outline": "0"}
    return {"border": "1px solid red"}


# Fade in/out email
@app.callback(Output("example-email", "disabled"), [Input("checklist-mail", "value")])
def disabled_mail(checklist_value: List) -> bool:
    """Disable email if not in the checklist.

    ...

    Parameters
    ----------
    checklist_value : List

    Returns
    -------
    bool
    """

    if not isinstance(checklist_value, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(checklist_value).__name__}"
        )
    if "email" not in checklist_value:
        return True
    return False


# disable job ID
@app.callback(Output("job-name", "disabled"), [Input("checklist-job-name", "value")])
def disable_job_name(checklist_value: List) -> bool:
    """Disable job name if not in the checklist.

    ...

    Parameters
    ----------
    checklist_value : List

    Returns
    -------
    bool
    """

    if not isinstance(checklist_value, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(checklist_value).__name__}"
        )
    if "job_name" not in checklist_value:
        return True
    return False


# (change_disabled_vcf_dropdown removed with the personal-variant VCF picker: the
# variant selector is now a single genome-driven dataset dropdown.)


# (The annotation selector was removed from the search form: searches now apply the
# annotations enabled in Settings -> Annotations, assembled by build_active_annotation.)


# select Cas protein from dropdown
# (select_cas_pam_dropdown removed: the Cas-protein selector was dropped, so the PAM
# dropdown lists all available PAMs directly with an enzyme-aware label.)


# add place holder to guide box
@app.callback([Output("text-guides", "placeholder")], [Input("radio-guide", "value")])
def change_placeholder_guide_textbox(guide_type: str) -> List:
    """Add place holders to guides text box.

    ...

    Parameters
    ----------
    guide_type : str
        Guide type

    Returns
    -------
    List
    """

    if not isinstance(guide_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide_type).__name__}")
    place_holder_text = ""
    if guide_type == "IP":  # individual spacers
        place_holder_text = str("GAGTCCGAGCAGAAGAAGAA\n" "CCATCGGTGGCCGTTTGCCC")
    elif guide_type == "GS":  # genomic sequences
        place_holder_text = str(
            ">sequence1\n"
            "AAGTCCCAGGACTTCAGAAGagctgtgagaccttggc\n"
            ">sequence_bed\n"
            "chr1 11130540 11130751\n"
            "chr1 1023000 1024000"
        )
    else:
        raise ValueError(f"Forbidden guide type ({guide_type})")
    assert bool(place_holder_text)
    return [place_holder_text]


def _preferred_variant(option_values: List[str]) -> str:
    """Default variant selection, chosen DYNAMICALLY from the installed datasets —
    NO hardcoded 1000G/HGDP or genome names, so it works for any future genome or
    dataset (e.g. a pig susScr11 + a custom VCF). Prefer a variant-aware search over
    reference-only, and among installed panels prefer the 'richest' one: the panel
    whose dataset tokens are a superset of the most other panels (e.g. a combined
    1000G+HGDP panel over either alone). Falls back to the first listed variant,
    then 'ref' when none are installed."""
    variants = [v for v in option_values if v and v != "ref"]
    if not variants:
        return "ref"

    def _covers(v: str) -> int:
        toks = set(v.split("_"))
        return sum(1 for o in variants if o != v and set(o.split("_")) <= toks)

    # most-covering (combined) first; ties broken by original listed order
    return max(variants, key=lambda v: (_covers(v), -variants.index(v)))


# change variants options
@app.callback(
    [Output("variant-dataset", "options"), Output("variant-dataset", "value")],
    [Input("available-genome", "value")],
)
def change_variant_dataset_options(genome_value: str) -> List:
    """Repopulate the variant-dataset dropdown for the selected genome.

    Genome-driven: only datasets actually installed for this genome are offered
    (built-in 1000G/HGDP for hg38, a combined entry when both are present), always
    with "Reference only" first. A genome with no variant data shows only
    "Reference only". The value is reset to a still-valid option so a stale
    selection from a previous genome cannot leak into the search.
    """

    if genome_value is not None and not isinstance(genome_value, str):
        raise TypeError(f"Expected {str.__name__}, got {type(genome_value).__name__}")
    options = get_variant_dataset_options(genome_value)
    # default to a variant-aware search, preferring the richest installed panel —
    # chosen dynamically (no hardcoded dataset/genome names).
    value = _preferred_variant([o["value"] for o in options])
    return [options, value]


# Limit the DNA/RNA bulge options to what a built index supports (hard cap).
# Bulge searches need a per-PAM TST index; 0-bulge searches run index-free, so 0
# is always available. A variant search also needs the variant index for each
# selected dataset, so the cap is the min across the reference and variant
# indexes. Mismatches are never index-limited, so they are left untouched.
@app.callback(
    [
        Output("dna", "options"),
        Output("rna", "options"),
        Output("dna", "value", allow_duplicate=True),
        Output("rna", "value", allow_duplicate=True),
        Output("bulge-guard-note", "children"),
    ],
    [
        Input("available-genome", "value"),
        Input("available-pam", "value"),
        Input("variant-dataset", "value"),
    ],
    [State("dna", "value"), State("rna", "value")],
    prevent_initial_call=True,
)
def limit_bulges_to_index(genome, pam, variant_choice, cur_dna, cur_rna):
    if not genome or not pam:
        return AV_BULGES, AV_BULGES, no_update, no_update, ""
    maxb = index_max_bulges(genome, pam, None)  # reference index (always needed)
    # scalar dropdown value -> list of datasets ("ref" -> none)
    # ANY selected dataset (built-in OR a custom VCF registered via Settings) — not a
    # hardcoded 1000G/HGDP whitelist, which previously skipped custom datasets and left
    # their variant-index bulge ceiling unchecked (overstating the bulge dropdown).
    selected = (
        []
        if variant_choice in (None, "", "ref")
        else [v for v in str(variant_choice).split("+") if v]
    )
    for v in selected:  # a variant bulge search also needs the variant index
        maxb = min(maxb, index_max_bulges(genome, pam, v))
    opts = [{"label": i, "value": i} for i in range(0, maxb + 1)]
    dna_v = cur_dna if (isinstance(cur_dna, int) and cur_dna <= maxb) else 0
    rna_v = cur_rna if (isinstance(cur_rna, int) and cur_rna <= maxb) else 0
    if maxb == 0:
        note = (
            "No bulge index for this genome/PAM"
            + (" + selected variant set" if selected else "")
            + " yet — only a fast 0-bulge search is available. Build an index in "
            "Settings to enable bulges."
        )
    else:
        note = f"Up to {maxb} DNA/RNA bulge(s) available (limited by the built index)."
    return opts, opts, dna_v, rna_v, note


@app.callback(
    [
        Output("available-pam", "options"),
        Output("available-pam", "value", allow_duplicate=True),
    ],
    [Input("available-genome", "value"), Input("variant-dataset", "value")],
    [State("available-pam", "value")],
    prevent_initial_call=True,
)
def update_pam_options(genome, variant_choice, current_pam):
    """Restrict the PAM list to what is searchable for the current genome + variant
    selection: all PAMs for a reference-only search, but only PAMs with a variant
    index (pamless NNN counts for all) when a variant dataset is included. Keeps the
    current PAM if it is still valid, else falls back to a sensible default."""
    options = get_pam_options(genome, variant_choice)
    valid = {o["value"] for o in options}
    if current_pam in valid:
        value = current_pam
    else:
        value = _default_pam(_default_cas()) if any(
            o["value"] == _default_pam(_default_cas()) for o in options
        ) else (options[0]["value"] if options else None)
    return options, value


def _default_genome() -> Optional[str]:
    """Sensible default genome: hg38 if installed, else the first available."""
    gs = [g["value"] for g in get_available_genomes()]
    return "hg38" if "hg38" in gs else (gs[0] if gs else None)


def _default_cas() -> Optional[str]:
    """Default nuclease: SpCas9 if installed, else the first available."""
    cs = [c["value"] for c in get_available_CAS()]
    return "SpCas9" if "SpCas9" in cs else (cs[0] if cs else None)


def _default_pam(cas: Optional[str]) -> Optional[str]:
    """Default PAM for a nuclease: prefer an NGG PAM, else the first for that Cas."""
    if not cas:
        return None
    pams = [
        p["value"]
        for p in get_available_PAM()
        if "-".join(p["value"].split(".")[0].split("-")[2:]) == cas
    ]
    if not pams:
        return None
    for p in pams:
        if "-NGG-" in p:
            return p
    return pams[0]


def index_page() -> html.Div:
    """Construct the layout of CRISPRme main page.
    When a new genome is added to /Genomes directory, reload genomes and PAMs
    dropdowns (via page reloading).

    ...

    Parameters
    ----------
    None

    Returns
    -------
    html.Div
    """

    # begin main page construction
    final_list = []
    # smart defaults, based on what is actually installed: hg38 + SpCas9/NGG +
    # 1000G variants + the standard 4/1/1 thresholds, so a non-expert can submit
    # a sensible search without configuring everything from scratch.
    _def_genome = _default_genome()
    _def_cas = _default_cas()
    _def_pam = _default_pam(_def_cas)
    # Default variant selection, chosen dynamically (mirrors change_variant_dataset_options
    # so the seed and the genome-change callback agree; no hardcoded dataset names).
    _def_variants = _preferred_variant(
        [o["value"] for o in get_variant_dataset_options(_def_genome)]
    )
    # seed the PAM dropdown options for the default nuclease so the default PAM
    # value is valid on first render (an empty options list makes Dash drop the
    # preset value before the cas->pam callback can populate it)
    _def_pam_options = [
        p
        for p in get_available_PAM()
        if _def_cas and "-".join(p["value"].split(".")[0].split("-")[2:]) == _def_cas
    ]
    # page intro — the CRISPRme+ wordmark logo above the first step, then the form
    # starts directly at "Select gRNA".
    introduction_content = html.Div(
        html.Img(
            src="assets/crisprme-logo.svg",
            style={"height": "80px"},
            alt="CRISPRme+",
        ),
        style={"textAlign": "center", "margin": "4px 0 20px"},
    )
    # warnings
    modal = html.Div(
        [
            dbc.Modal(
                [
                    dbc.ModalHeader("WARNING! Missing or wrong input"),
                    dbc.ModalBody(
                        str(
                            "The following inputs are missing, please select "
                            "values before submitting the job"
                        ),
                        id="warning-list",
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="close", className="modal-button")
                    ),
                ],
                id="modal",
                centered=True,
            ),
        ]
    )
    # guides table
    tab_guides_content = html.Div(
        [
            html.H4("Select gRNA"),
            dcc.RadioItems(
                id="radio-guide",
                options=[
                    {"label": " Input individual spacer(s)", "value": "IP"},
                    {"label": " Input genomic sequence(s)", "value": "GS"},
                ],
                value="IP",
            ),
            dcc.Textarea(
                id="text-guides",
                placeholder=str("GAGTCCGAGCAGAAGAAGAA\n" "CCATCGGTGGCCGTTTGCCC"),
                style={"width": "300px", "height": "60px", "fontSize": "1rem"},
            ),
            dbc.FormText(
                str(
                    "Spacer must be provided as a DNA sequence without a PAM. "
                    "Multiple spacers can be entered, one per line. To process "
                    "many spacers at once we recommend the command-line version "
                    "of CRISPRme+ (see the offline instructions below)."
                ),
                color="secondary",
                style={"fontSize": "1.4rem"},
            ),
            # Load-example lives at the bottom of this box so a new user can populate
            # the whole form with one click right where they start.
            html.Div(
                html.Button(
                    "Load Example",
                    id="load-example-button",
                    style={"background-color": "#E6E6E6", "width": "300px"},
                ),
                style={"textAlign": "left", "marginTop": "12px"},
            ),
        ],
        style={"width": "300px"},  # NOTE same as text-area
    )
    # cas protein dropdown
    # PAM dropdown. The Cas-protein selector was removed as redundant: the PAM value
    # already encodes the enzyme (e.g. 20bp-NGG-SpCas9) and the label now shows it
    # (e.g. "SpCas9 · NGG"), so a single self-describing PAM dropdown suffices.
    pam_content = html.Div(
        [
            html.H4("Select PAM"),
            html.Div(
                dcc.Dropdown(
                    options=get_available_PAM(),
                    value=_def_pam,
                    clearable=False,
                    id="available-pam",
                    style={"width": "300px"},
                )
            ),
        ],
    )
    personal_data_management_content = html.Div(
        [
            html.Br(),
            html.A(
                html.Button(
                    "Settings / Data Manager",
                    id="add-genome",
                    style={"display": DISPLAY_OFFLINE},
                ),
                href=os.path.join(URL, "settings"),
                target="",
                style={"text-decoration": "none", "color": "#555"},
            ),
        ]
    )
    # genome dropdown
    genome_content = html.Div(
        [
            html.H4("Select genome"),
            html.Div(
                dcc.Dropdown(
                    options=get_available_genomes(),
                    value=_def_genome,
                    clearable=False,
                    id="available-genome",
                ),
                style={"width": "300px"},
            ),
            html.P("Variants", style={"margin": "8px 0 2px"}),
            html.Div(
                dcc.Dropdown(
                    options=get_variant_dataset_options(_def_genome),
                    value=_def_variants,
                    clearable=False,
                    id="variant-dataset",
                    style={"width": "300px"},
                ),
            ),
        ]
    )
    # thresholds boxes
    thresholds_content = html.Div(
        [
            html.H4("Select thresholds"),
            # PRIMARY control: a single "maximum total edits" slider (mismatches +
            # bulges). This is the simple, non-expert knob; the per-type mismatch /
            # bulge limits live under "Advanced" below and stay wide open by default
            # so the slider is the governing constraint (CRISPRme issue #107).
            html.Div(
                [
                    html.P(
                        "Maximum edits (mismatches + bulges)",
                        style={"margin-bottom": "2px", "font-weight": "600"},
                    ),
                    dcc.Slider(
                        id="max-edits-slider",
                        min=1,
                        max=5,
                        step=1,
                        value=5,
                        marks={i: str(i) for i in range(1, 6)},
                        tooltip={"placement": "bottom", "always_visible": False},
                    ),
                    html.P(
                        "Total number of differences (mismatches + DNA/RNA bulges) "
                        "allowed between a guide and an off-target. The default is 5 (matching the "
                        "command line); lower it for a faster, narrower search. "
                        "The on-target (0 edits) is always reported at any setting.",
                        style={"font-size": "1.45rem", "color": "#555"},
                    ),
                    html.P(
                        "Note: this limit is applied during the search against the "
                        "variant-enriched genome. A variant off-target can therefore "
                        "appear with a slightly higher mismatch+bulge count in the "
                        "results (its count is reported against the reference); "
                        "reference off-targets always stay within the limit.",
                        style={"font-size": "1.4rem", "color": "#777", "font-style": "italic"},
                    ),
                ],
                style={"max-width": "420px", "margin-bottom": "12px"},
            ),
            # ADVANCED: per-type mismatch / bulge caps, hidden by default.
            dbc.Button(
                "Advanced options ▾",
                id="advanced-thresholds-toggle",
                color="link",
                n_clicks=0,
                style={"padding": "0", "font-size": "1rem"},
            ),
            dbc.Collapse(
                html.Div(
                    [
                        html.P(
                            "Per-type caps. Left at their maxima the slider above "
                            "governs; lower them to further restrict a single type.",
                            style={"font-size": "1.4rem", "color": "#555"},
                        ),
                        html.Div(  # mismatches box
                            [
                                html.P("Mismatches"),
                                dcc.Dropdown(
                                    options=AV_MISMATCHES,
                                    value=6,
                                    clearable=False,
                                    id="mms",
                                    style={"width": "60px"},
                                ),
                            ],
                            style={"display": "inline-block", "margin-right": "20px"},
                        ),
                        html.Div(  # DNA bulges box
                            [
                                html.P(["DNA", html.Br(), "Bulges"]),
                                dcc.Dropdown(
                                    options=AV_BULGES,
                                    value=2,
                                    clearable=False,
                                    id="dna",
                                    style={"width": "60px"},
                                ),
                            ],
                            style={"display": "inline-block", "margin-right": "20px"},
                        ),
                        html.Div(  # RNA bulges box
                            [
                                html.P(["RNA", html.Br(), "Bulges"]),
                                dcc.Dropdown(
                                    options=AV_BULGES,
                                    value=2,
                                    clearable=False,
                                    id="rna",
                                    style={"width": "60px"},
                                ),
                            ],
                            style={"display": "inline-block"},
                        ),
                        html.Div(
                            id="bulge-guard-note",
                            style={
                                "font-size": "1.4rem",
                                "color": "#555",
                                "margin-top": "6px",
                            },
                        ),
                    ],
                    style={"margin-top": "8px"},
                ),
                id="advanced-thresholds-collapse",
                is_open=False,
            ),
        ],
    )
    # base editing boxes
    base_editing_content = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        html.H4("Base editing?"),
                        style={"display": "inline-block", "margin-right": "20px"},
                    ),
                    html.Div(
                        dcc.RadioItems(
                            id="radio-base_editor",
                            options=[
                                {"label": "Yes", "value": "Y"},
                                {"label": "No", "value": "N"},
                            ],
                            value="N",
                            labelStyle={
                                "margin-right": "5px",
                                "display": "inline-block",
                            },
                        ),
                        style={"display": "inline-block"},
                    ),
                ]
            ),
            html.Div(
                [
                    html.Div(  # BE window start dropdown
                        [
                            html.P("Window start"),
                            dcc.Dropdown(
                                clearable=False,
                                id="be-window-start",
                                style={"width": "60px"},
                            ),
                        ],
                        style={"display": "inline-block", "margin-right": "20px"},
                    ),
                    html.Div(  # BE window stop dropdown
                        [
                            html.P("Window stop"),
                            dcc.Dropdown(
                                clearable=False,
                                id="be-window-stop",
                                style={"width": "60px"},
                            ),
                        ],
                        style={"display": "inline-block", "margin-right": "20px"},
                    ),
                    html.Div(  # BE nucleotides dropdown
                        [
                            html.P(["Nucleotide"]),
                            dcc.Dropdown(
                                options=BE_NTS,
                                clearable=False,
                                id="be-nts",
                                style={"width": "60px"},
                            ),
                        ],
                        style={"display": "inline-block", "margin-right": "20px"},
                    ),
                ],
                id="div-base-editor-dropdowns",
                style={"display": "none"},
            ),
        ],
        style={"margin-top": "16px", "border-top": "1px solid #eef2f4", "padding-top": "12px"},
    )
    # annotations dropdown
    # (annotation selector removed: searches apply the annotations enabled in
    # Settings -> Annotations; nothing to pick on the form.)
    # mail box
    mail_content = html.Div(
        [
            dcc.Checklist(
                options=[
                    {
                        "label": " Notify me by email",
                        "value": "email",
                        "disabled": False,
                    }
                ],
                id="checklist-mail",
                value=[],
            ),
            html.Div(
                dbc.Input(
                    type="email",
                    id="example-email",
                    placeholder="name@mail.com",
                    className="exampleEmail",
                    disabled=True,
                    style={"width": "300px"},
                )
            ),
            dbc.FormText(
                [
                    "Tick the box and enter your address to get an email when the "
                    "job finishes. Sending requires a one-time mail-server setup in ",
                    html.A(
                        "Settings → Email notifications",
                        href=os.path.join(URL, "settings"),
                    ),
                    " (SMTP host, sender address and an app password). Until that is "
                    "configured the job still runs — it just won't email you.",
                ],
                color="secondary",
                style={"fontSize": "1.4rem"},
            ),
        ]
    )
    # job name box
    job_name_content = html.Div(
        [
            dcc.Checklist(
                options=[
                    {"label": " Job name", "value": "job_name", "disabled": False}
                ],
                id="checklist-job-name",
                value=[],
            ),
            html.Div(
                dbc.Input(
                    type="text",
                    id="job-name",
                    placeholder="my_job",
                    className="jobName",
                    disabled=True,
                    style={"width": "300px"},
                )
            ),
        ]
    )
    # submit button
    submit_content = html.Div(
        [
            html.Button(
                "Submit",
                id="check-job",
                style={"background-color": "#E6E6E6", "width": "300px"},
            ),
            html.Button("", id="submit-job", style={"display": "none"}),
        ],
        style={"textAlign": "left"},  # left-align to match the inputs/dropdowns
    )
    # (Load Example button moved to the bottom of the "Select gRNA" box above.)
    # terms and conditions link
    terms_and_conditions_content = html.Div(
        [
            html.Div("By clicking submit you are agreeing to the"),
            html.Div(
                html.A(
                    "Terms and Conditions.",
                    target="_blank",
                    href=f"{GITHUB_LINK}/blob/main/LICENSE",
                )
            ),
        ]
    )
    # insert introduction in the page layout
    final_list.append(introduction_content)
    # add other content
    final_list.append(
        html.Div(
            [modal]
            + [
                # one numbered step-card per stage: a single centered column that
                # reads top-to-bottom (Guide -> Genome+variants -> PAM -> Thresholds
                # -> Annotation -> Run). Content blocks are reused verbatim, so all
                # component ids / callbacks are unchanged.
                html.Div(
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Div(
                                    str(_n),
                                    style={
                                        "width": "34px",
                                        "height": "34px",
                                        "borderRadius": "50%",
                                        "backgroundColor": "#388396",
                                        "color": "white",
                                        "fontWeight": "700",
                                        "fontSize": "1.1rem",
                                        "display": "flex",
                                        "alignItems": "center",
                                        "justifyContent": "center",
                                    },
                                ),
                                width="auto",
                            ),
                            dbc.Col(_content),
                        ],
                        align="start",
                    ),
                    style={
                        "backgroundColor": "white",
                        "border": "1px solid #e2e8ec",
                        "borderRadius": "10px",
                        "padding": "16px 22px",
                        "marginBottom": "14px",
                        "boxShadow": "0 1px 4px rgba(0,0,0,0.05)",
                    },
                )
                for _n, _content in enumerate(
                    [
                        tab_guides_content,
                        genome_content,
                        pam_content,
                        html.Div([thresholds_content, base_editing_content]),
                        html.Div(
                            [
                                mail_content,
                                job_name_content,
                                html.Br(),
                                submit_content,
                                terms_and_conditions_content,
                            ]
                        ),
                    ],
                    start=1,
                )
            ],
            style={"width": "100%"},
            id="steps-background",
        )
    )
    final_list.append(html.Br())
    final_list.append(
        html.P(
            str(
                "*The offline version of CRISPRme can be downloaded from GitHub "
                "and offers additional functionalities, including the option to "
                "input personal data (such as genetic variants, annotations, "
                "and/or empirical off-target results) as well as custom PAMs and "
                "genomes. There is no limit on the length or number of spacers, "
                "mismatches, and/or bulges used in the offline search."
            )
        )
    )
    # constrain the whole page to the same centered width as the step-cards so the
    # intro text and footer line up with the submitting form
    index_page = html.Div(
        final_list,
        style={
            "maxWidth": "680px",
            "margin": "0 auto",
            "padding": "24px 12px",
            # larger, more legible base font for the whole submission form (headings,
            # labels, help text and radio/checkbox labels all inherit from here)
            "fontSize": "1.05rem",
        },
    )
    return index_page


@app.callback(
    Output("div-base-editor-dropdowns", "style"), [Input("radio-base_editor", "value")]
)
def update_visibility_base_editor_dropdowns(radio_value: str) -> Dict:
    """Update visibilyt of base editing dropdowns.
    default is display none.
    ...

    Parameters
    ----------
    radio_value : str

    Returns
    -------
    Dict
    """

    if radio_value == "Y":
        return {"display": ""}
    else:
        return {"display": "none"}


@app.callback(
    [
        Output("advanced-thresholds-collapse", "is_open"),
        Output("max-edits-slider", "disabled"),
        Output("advanced-thresholds-toggle", "children"),
    ],
    [Input("advanced-thresholds-toggle", "n_clicks")],
    [State("advanced-thresholds-collapse", "is_open")],
    prevent_initial_call=True,
)
def toggle_advanced_thresholds(n_clicks: int, is_open: bool) -> Tuple:
    """Toggle the Advanced (per-type mismatch/bulge) panel. Opening it switches to
    the explicit per-type mode and GRAYS OUT the max-edits slider to make clear the
    total-edits cap no longer governs; closing it restores the simple slider mode."""
    new_open = not is_open
    label = "Advanced options ▴" if new_open else "Advanced options ▾"
    # slider disabled (grayed out) exactly when the advanced panel is open
    return new_open, new_open, label


@app.callback(
    [Output("be-window-start", "options"), Output("be-window-stop", "options")],
    [Input("text-guides", "value")],
    [State("radio-guide", "value"), State("available-genome", "value")],
)
def update_base_editing_dropdown(
    text_guides: str, guide_type: str, genome: str
) -> Tuple:
    """Update base editing dropdown dinamically. The start and stop values for
    base editing are changed accordingly to the guides provided in input by
    the user.

    ...

    Parameters
    ----------
    text_guides : str
        Guides
    guide_type : str
        Guide type
    genome : str
        Reference genome

    Returns
    -------
    Tuple
    """

    if text_guides is not None:
        if not isinstance(text_guides, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(text_guides).__name__}"
            )
    if not isinstance(guide_type, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guide_type).__name__}")
    dropdown_options = [{"label": "", "value": ""}]
    if text_guides is None:
        return dropdown_options, dropdown_options
    if guide_type == "IP":  # individual spacers
        guides = text_guides.strip()
    elif guide_type == "GS":  # genomic sequences
        guides = list()
        for seqname_and_seq in text_guides.split(">"):
            if not seqname_and_seq:
                continue
            seqname = seqname_and_seq[: seqname_and_seq.find("\n")]
            seq = seqname_and_seq[seqname_and_seq.find("\n") :].strip()
            if "chr" in seq:  # BED regions
                for line in seq.split("\n"):
                    if not line:
                        continue
                    line_split = line.strip().split()
                    # line_split = re.split(r";|,|.|:|-| ", line.strip())
                    # print(line_split)
                    seq_read = f"{line_split[0]}:{line_split[1]}-{line_split[2]}"
                    seq_read = extract_seq.extractSequence(
                        seqname, seq_read, genome.replace(" ", "_")
                    )
            else:
                seq_read = "".join(seq.split()).strip()
            guides.append(seq_read)
        guides = "\n".join(list(set(guides)))
    if not all(
        [len(guide) == len(guides.split("\n")[0]) for guide in guides.split("\n")]
    ):
        guides = select_same_len_guides(guides)
    guides = guides.split("\n")
    dropdown_options = [{"label": i, "value": i} for i in range(1, len(guides[0]) + 1)]
    return dropdown_options, dropdown_options


def check_mail_address(mail_address: str) -> bool:
    """Check mail address consistency.

    ...

    Parameters
    ----------
    mail_address : str
        Mail address

    Returns
    -------
    bool
    """

    if not mail_address:  # check wether is None or empty
        return False
    assert mail_address is not None
    if not isinstance(mail_address, str):
        raise TypeError(f"Expected {str.__name__}, got {type(mail_address).__name__}")
    mail_address_fields = mail_address.split("@")
    if len(mail_address_fields) > 1 and bool(mail_address_fields[-1]):
        return True
    return False
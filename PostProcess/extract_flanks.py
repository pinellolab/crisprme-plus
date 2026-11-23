#!/usr/bin/env python3
"""``crisprme.py extract-flanks``: the genomic sequence around an off-target site, per
observed haplotype, to order validation oligos from.

A shortlist of off-targets (the report's top-N download, or an integrated results table)
says WHERE each site is. To validate one in the lab you need WHAT IS THERE: the sequence
around it, to place amplicon primers in and to check your oligos against. For a site in a
variant region that is not one sequence. The reference has one; each distinct haplotype
observed in the search's panel has another; and the haplotype that carries the variant
which MAKES the off-target is usually not the reference one. Ordering primers off the
reference for a site that exists only on an ALT haplotype is how a validation assay
quietly measures the wrong locus.

So this writes, per input row:

  * the REFERENCE flank, and
  * one flank per DISTINCT observed haplotype sequence, deduplicated by sequence, each
    with the samples (and, where the VCF is phased, the haplotype side) that carry it.

A site with no variant in its window yields one entry. A site whose haplotypes all read
alike yields one entry. The target's own coordinates are marked in every entry -- the
offset and length of the protospacer+PAM inside that flank -- so nobody has to count
bases, and so an indel that shifts the site is visible as a shifted offset.

HOW A HAPLOTYPE IS WRITTEN. By ``scripts/variant_haplotypes.py``: the reference with
every ALT allele a sample's genotype puts on that slot, each in place of its REF bases, as
VCF 4.x defines it. That is the same code ``scripts/check_indel_rows.py`` is built on --
the independent check of CRISPRme's indel rows, hardened against real data. Two readers
of one VCF that disagreed would be two different answers to "what does this individual
hold here?", and one of them would send a user's oligos to the wrong sequence.

WHAT IT NEVER DOES. It never invents a haplotype. A sample whose call is unphased, has a
missing allele, or whose alleles overlap in the window (the VCF gives no sequence there)
is REPORTED, under ``unresolved``, with the reason and the records at fault -- not guessed
and not silently dropped. A target near a contig end is truncated and said so. An input it
cannot resolve is named exactly, as the checker names it, never guessed.

THE WINDOW is the target plus a padding at each of its two ends, named for what each
sits beside: ``--upstream`` beyond the PROTOSPACER's outer end and ``--downstream``
beyond the PAM's outer end, both read in the GUIDE's orientation. Which GENOMIC side
each lands on follows from the row's strand and from which end of the target the PAM is
-- read from the run's own records, never assumed, because on a 5' PAM (Cas12a) the two
ends swap. ``--window W`` is the symmetric form, and is exactly ``--upstream W
--downstream W``. A ``--preset`` names a configuration, so another assay is a line in
``flank_presets.conf`` rather than a change to this file.

OUTPUT, three ways of reading one computation.
  * ``<out>.fa`` -- one record per emitted sequence, no sample identities in any name.
  * ``<out>.tsv`` -- one line per emitted sequence, joining back to the input row by its
    line number and, where the input has one, its rank. Deterministic: rows in input
    order, the reference entry first, then the haplotypes that HOLD the off-target, then
    by carrier count and by the sequence itself; carriers and variants sorted.
  * ``--annotate-table`` -- the INPUT table written back out with the sequences added as
    columns (the reference sequence, the primary ALT, how many distinct ALTs the row
    has), so a results table becomes an oligo order without a join. A row with several
    distinct ALTs keeps them in the two files above, reached by ``flank_row_id``: the
    ``--max-haplotypes`` most useful (DEFAULT 10, and the cap is on those two files,
    not only on the cell), every one of them at ``--max-haplotypes 0``, with
    ``n_haplotypes`` / ``flank_alt_count`` always the true total and the notes counting
    what was left out. None is crammed into a cell unless ``--alt-list`` asks for it.

``--filter`` chooses which rows are done at all, as an expression over the table's OWN
columns, with the CRISPR-Bulge sensitive operating point available by name
(``CB_SENSITIVE``) rather than as a number to copy.

COORDINATES. CRISPRme's ``Start_coordinate_*`` / curated ``Position`` is the 0-based
reference coordinate of the target's leftmost base, on either strand. Flank coordinates
in the TSV are printed 1-based inclusive (as a genome browser reads them) and are always
GENOMIC; the target offset inside a flank is 0-based, as a Python slice reads it.
ORIENTATION: by default a flank is written as the GUIDE reads it, so the protospacer in
the file runs 5'->3' and an oligo can be ordered off it directly; ``--orientation
genome`` writes the reference's plus strand instead, whatever the row's strand. Which
one was used is in every output line (``orientation``, ``emitted_strand``), and the
TSV's ``protospacer_on_strand`` always gives the target read on the row's own strand.
"""

import ast
import contextlib
import csv
import glob
import gzip
import os
import re
import sys
import zipfile
from collections import defaultdict

# the VCF reader, the reference FASTA reader and the haplotype writer are shared with
# scripts/check_indel_rows.py (which they were built for and hardened in): a user's
# validation sequence and the release gate's notion of what a haplotype holds must not be
# two different pieces of code. scripts/ ships beside PostProcess/ in every install.
_SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "scripts"))
if _SCRIPTS not in sys.path:
    sys.path.append(_SCRIPTS)
try:
    from run_records import read_panel, run_params  # noqa: E402  (its folder, above)
    from variant_haplotypes import (  # noqa: E402  (its folder is settled just above)
        PASS_FILTERS, Allele, Reference, Vcf, make_hap, on_slot, parse_calls, rc,
        reads_at, span_of, touches,
    )
except ImportError as _exc:  # pragma: no cover - a broken install, named not guessed
    raise ImportError(
        "extract-flanks needs variant_haplotypes.py and run_records.py (the shared VCF, "
        "reference and run-sidecar readers it writes haplotypes with), looked for in %s. "
        "Every CRISPRme install ships scripts/ beside PostProcess/." % _SCRIPTS) from _exc


class FlankError(Exception):
    """An input the command cannot resolve, or a run it will not model. Its text is the
    whole message: it names what is missing and what to pass, never a traceback."""


# THE WINDOW. Two forms of one thing. A window is the target plus a padding at each of
# its two ends, and the ends are named for what they sit beside:
#
#   --upstream N    N bases beyond the PROTOSPACER's outer end (its PAM-distal side)
#   --downstream M  M bases beyond the PAM's outer end
#   --window W      the symmetric form: the same as --upstream W --downstream W
#
# Both are read in the GUIDE's orientation, so they do not change meaning with the row's
# strand or with where the nuclease's PAM sits; which genomic side each lands on is
# worked out per row (Geometry.edges). For a 3' PAM (SpCas9: protospacer then PAM)
# --upstream is literally "before the protospacer's first base" and --downstream
# "after the PAM's last base". For a 5' PAM (Cas12a: PAM then protospacer) the two ends
# of the target swap, so the SAME two flags land on the opposite genomic sides -- which
# is the point of naming them for the protospacer and the PAM rather than for left and
# right.
#
# THE DEFAULT IS THE ANCHORED FORM, 3 and 6. Three reasons. (1) An anchored window is
# the general one: its meaning does not depend on how long the target is, while a
# symmetric window is only the special case where the two paddings are equal, and is
# still reachable in one flag (--window W). (2) The assay these numbers come from -- a
# short oligo read straight off the site, no amplicon -- is the case where getting the
# padding wrong is unrecoverable, because there is no primer slack to absorb it. (3) A
# small default keeps a 1,000-row shortlist a file you can open; a user who wants the
# amplicon-scale context asks for it, and gets it from one flag or the `amplicon` preset.
# The length is never left to be inferred: every output says the window it used and the
# bases it actually emitted (window_spec, pad_upstream, pad_downstream, flank_len).
DEFAULT_UPSTREAM = 3
DEFAULT_DOWNSTREAM = 6
# the symmetric alternative. An amplicon for validation sequencing (rhAmpSeq, CRISPResso2
# and the like) is typically 200-300 bp with its primers 50-150 bp away from the cut, and
# a Sanger/TIDE amplicon 400-700 bp: 200 + target + 200 covers the first comfortably and
# the short end of the second. It is the `amplicon` preset and the value --window takes
# when it is given no number of its own.
DEFAULT_WINDOW = 200
# how many haplotypes of a site are written beside the reference. How many a site HAS
# depends strongly on the window: at --window 200 a panel of thousands gives dozens (on
# the 1000-site shortlist of a genome-wide HGDP+1000G run: a median of 23, and one site
# with 4264), nearly all differing from each other far from the target, while a 32-mer
# spans too little to hold many. Writing them all turns a shortlist into megabytes, so
# the most useful come first -- the sequences that HOLD the off-target, then by how many
# haplotypes carry them -- and the rest are counted in the notes rather than written.
# --max-haplotypes 0 writes every one. n_haplotypes always gives the true total.
DEFAULT_MAX_HAPLOTYPES = 10
# the alignment blocks an integrated results table carries, and the one taken by default
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
DEFAULT_BLOCK = "highest_CFD"
_NO_VALUE = ("", "NA", "n", "-", ".")

# ---- where the nuclease's PAM sits, relative to the protospacer, in the GUIDE's
# orientation. Read from the run's own records, never assumed: a 5' PAM puts the two
# paddings on the opposite genomic sides from a 3' one, so guessing SpCas9 on a Cas12a
# run would write every window back to front.
PAM_3PRIME = "3prime"     # SpCas9 and the rest of the NGG/NRG family: spacer then PAM
PAM_5PRIME = "5prime"     # Cas12a (TTTV), and anything else that reads PAM then spacer
PAM_SIDES = (PAM_3PRIME, PAM_5PRIME)

# ---- what to do when the target's genomic footprint is not the nominal spacer+PAM
# length, which is what a bulge in the alignment means: a DNA bulge is an extra genomic
# base (the footprint is longer), an RNA bulge is a base of the guide with no genomic
# partner (shorter). See Geometry.edges.
POLICY_ANCHORED = "anchored"        # keep the paddings, let the total vary (default)
POLICY_FIXED_TOTAL = "fixed-total"  # keep the total, move a padding
LENGTH_POLICIES = (POLICY_ANCHORED, POLICY_FIXED_TOTAL)

# ---- which orientation the emitted sequence is written in
ORIENT_GUIDE = "guide"    # as the guide reads it: the protospacer 5'->3' (default)
ORIENT_GENOME = "genome"  # the reference's plus strand, whatever the row's strand
ORIENTATIONS = (ORIENT_GUIDE, ORIENT_GENOME)
DEFAULT_ORIENTATION = ORIENT_GUIDE

# ---- named window configurations. A preset is a NAME for a set of the flags below, so
# that an assay with its own geometry is a line of configuration rather than a change to
# this file: --preset-file <path> and $CRISPRME_FLANK_PRESETS add to (and override) these
# by name, and --help lists whatever is loaded. Nothing here is special-cased anywhere
# else in the module: a preset only supplies defaults for flags the user did not give.
PRESET_KEYS = ("upstream", "downstream", "window", "length_policy", "orientation",
               "max_haplotypes", "filter", "description")
BUILTIN_PRESETS = {
    "oligo": {
        "upstream": "3", "downstream": "6",
        "description": "a short oligo read straight off the site: 3 nt beyond the "
                       "protospacer, 6 nt beyond the PAM (a 32-mer for a 20-nt spacer "
                       "with a 3-nt PAM). The default geometry.",
    },
    "amplicon": {
        "window": "200",
        "description": "amplicon-scale context: 200 bases beyond each end of the "
                       "target (a 423 bp interval for a 23-nt target), room to place "
                       "primers 50-150 bp from the cut.",
    },
}
# a preset file shipped beside this module, so the built-in set itself is editable
# without touching code; then $CRISPRME_FLANK_PRESETS, then --preset-file (last wins).
PRESET_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "flank_presets.conf")
PRESET_ENV = "CRISPRME_FLANK_PRESETS"

# What the REFERENCE entry's note may say, and the three different things it means. The
# reference sequence is the same in all three; what differs is how much was CHECKED, and
# only the first is a statement about the panel's samples. Writing it on a row nobody was
# read for would be a positive claim no VCF was opened to make.
REF_NOTE_CHECKED = "reference-of-every-sample-without-a-variant-here"
REF_NOTE_NO_CONTIG = "no-vcf-covers-this-contig"   # VCFs read, none holds this contig
REF_NOTE_NO_VCF = "reference-only"                 # no VCF was read at all
# a VCF holds the contig but has no sample columns: there were no genotypes to check
REF_NOTE_NO_GENOTYPES = "no-genotyped-vcf-for-this-contig"

TSV_COLUMNS = (
    "row_id",            # the input row's 1-based data line number: the join key
    "rank",              # the input's own rank column, where it has one
    "guide",             # the input's Spacer+PAM, where it has one
    "chrom",
    "target_start",      # 1-based inclusive, on the reference
    "target_end",        # 1-based inclusive, on the reference
    "target_strand",
    "haplotype",         # ref | hap1, hap2, ... | unresolved
    "fasta_id",
    "flank_start",       # 1-based inclusive, on the reference (ALWAYS genomic)
    "flank_end",         # 1-based inclusive, on the reference (ALWAYS genomic)
    "flank_len",         # the bases of THIS haplotype between them (an indel moves it)
    "target_offset",     # 0-based, into the flank sequence of this entry AS EMITTED
    "target_len",        # bases of this entry's flank the target covers
    "target_found",      # yes: this haplotype reads the row's protospacer+PAM here
    "target_sequence",   # those bases, plus strand
    "protospacer_on_strand",   # the same read on the row's strand
    "n_haplotypes",      # distinct sequences this row gave (the reference counted),
                         # whether or not --max-haplotypes wrote them all
    "n_samples",         # panel samples with a record in this flank (every other sample
                         # reads as the reference here)
    "n_carriers",        # reconstructed haplotype slots with this sequence
    "carriers",          # sample:slot, sorted; blank under --anonymise. On an
                         # unresolved line: sample:slot where only that copy is
                         # unknown (overlapping records), the bare sample where
                         # neither copy is (an unphased or missing genotype)
    "variants",          # the records written into this flank
    # ---- the geometry this line was cut with, so nobody has to infer it from the
    # command line that produced the file (or from a length that a bulge or an indel
    # moved). Every one of these is a property of the ROW, not of the run.
    "orientation",       # guide | genome: which way the sequence above is written
    "emitted_strand",    # + | -: the genomic strand the emitted sequence reads
    "window_spec",       # e.g. up3+down6, or window200 for the symmetric form
    "pad_upstream",      # bases emitted beyond the protospacer's outer end, AS APPLIED
    "pad_downstream",    # bases emitted beyond the PAM's outer end, AS APPLIED
    "pam_side",          # 3prime | 5prime | unknown: which end of the target is the PAM
    "length_policy",     # anchored | fixed-total
    "nominal_target_len",   # spacer+PAM length, gaps taken out: 23 for a 20-nt SpCas9
    "target_span_ref",   # reference coordinates the target covers (a bulge moves it)
    "bulge_type",        # none | DNA:n | RNA:n | DNA:n,RNA:n | unknown
    "notes",
)

# a row whose genomic footprint is not the nominal spacer+PAM length: the class a fixed
# total-length assumption gets wrong. Greppable, and the three columns above say by how
# much and why.
NOTE_FOOTPRINT = "target-span-not-nominal"
NOTE_LENGTH_ADJUSTED = "padding-moved-to-hold-the-total"
# A haplotype that DELETES the whole window has no sequence at all -- a real case at a
# short window (one entry of 4,017 on a 1,000-row shortlist at 3/6, none of 10,354 of the
# same shortlist at --window 200, because a 32 bp window fits inside a deletion that a
# 423 bp one does not). It is a fact about that genome, so the line stays and says this;
# but a FASTA record with no bases, or an empty element in a '|'-delimited cell, is how a
# reader ends up pairing every later sequence with the wrong haplotype, so the record is
# left out and the cell says SEQ_DELETED where the sequence would be.
NOTE_WINDOW_DELETED = "window-deleted-on-this-haplotype"
SEQ_DELETED = "(window-deleted)"


# --------------------------------------------------------------------------- #
# the input table
# --------------------------------------------------------------------------- #
@contextlib.contextmanager
def open_table(path):
    """A results table as a stream of lines: a .tsv, a .tsv.gz or the one .tsv inside a
    .zip (the shapes CRISPRme's downloads come in). Closed on the way out, whatever
    happens: a top-1000 run opens one file, a scripted one opens thousands."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist()
                     if n.endswith(".tsv") or n.endswith(".txt")]
            if len(names) != 1:
                raise FlankError(
                    "%s holds %d .tsv/.txt members, not one, so which is the table is "
                    "not ours to choose: unzip it and pass one of %s" % (
                        path, len(names), ", ".join(sorted(names)[:8])
                        + (", ..." if len(names) > 8 else "")))
            with zf.open(names[0]) as fh:
                yield (line.decode("utf-8", "replace").rstrip("\r\n") for line in fh)
        return
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        # "\r\n" and not "\n": a table written on Windows would otherwise leave a
        # carriage return on the end of its last column, in the notes of every line
        # and in the middle of a line the column mode writes back out
        yield (line.rstrip("\r\n") for line in fh)


def table_reader(lines):
    """The one way this module parses a results table: tab-separated fields, ONE ROW PER
    PHYSICAL LINE.

    QUOTE_NONE, and that is the point. A TSV is not a CSV: a results table has no quoting
    convention, a double quote in a cell is a double quote, and csv's default dialect
    would instead swallow it and -- where a cell opens a quote that closes on a later line
    (Excel's tab-delimited export writes exactly that for a cell holding a newline) --
    JOIN several lines into one row. The row number this returns is the key every output
    joins on (``flank_row_id``, the FASTA name) and the key ``write_annotated_table``
    reads off the raw text of the same file, so a reader that counted lines differently
    from the writer would put one site's sequence on another site's row, silently, with
    the right-looking number of columns. write_annotated_table refuses if the two ever
    disagree; this is why they cannot.
    """
    return csv.reader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)


# Where the guide itself is named, in the order preferred. ``Aligned_spacer+PAM`` is the
# report's curated download -- the input of both README examples -- which has no plain
# ``Spacer+PAM`` column at all, so without it the guide column of the primary documented
# input is empty on every line. It is the guide written against the row's alignment, so
# on a DNA-bulge row it carries the alignment's gap (81 of the 1,000 rows of a real
# top-N): the gap is a property of that row, not of the guide the file is about, and
# read_rows takes it out, so several guides' files merge into one order on this column.
_GUIDE_COLUMNS = ("Spacer+PAM", "crRNA", "Aligned_spacer+PAM")


def _column(head, names):
    for n in names:
        if n in head:
            return head[n]
    return None


def table_columns(header, block):
    """(column indexes, the shape's name) of a results table, or raise. Two shapes are
    read: the integrated results table (per-block columns, one block chosen) and the
    curated download the report writes (one projection; ``rank``, ``Position``,
    ``Aligned_protospacer+PAM``). Both are keyed by name, never by position."""
    head = {n: i for i, n in enumerate(header)}
    chrom = _column(head, ("Chromosome", "chrom"))
    if chrom is None:
        raise FlankError("the table has no Chromosome column: is it a CRISPRme results "
                         "table or the report's download?")
    integrated = "Start_coordinate_(%s)" % block
    if integrated in head:
        return {
            "chrom": chrom,
            "start": head[integrated],
            "strand": head["Strand_(%s)" % block],
            "aligned": head["Aligned_protospacer+PAM_ALT_(%s)" % block],
            "aligned_ref": _column(head, ("Aligned_protospacer+PAM_REF_corrected_(%s)" % block,
                                          "Aligned_protospacer+PAM_REF_(%s)" % block)),
            "rank": _column(head, ("rank",)),
            "guide": _column(head, _GUIDE_COLUMNS),
            # the guide side of the alignment AS WRITTEN, gaps and all: the PAM's side
            # reads off its run of N's and a DNA bulge off its gaps, and both are lost
            # once the gaps are taken out for the ``guide`` column
            "spacer": _column(head, ("Aligned_spacer+PAM_(%s)" % block,
                                     "Aligned_spacer+PAM")),
        }, "integrated (%s)" % block
    start = _column(head, ("Position", "Start_coordinate"))
    aligned = _column(head, ("Aligned_protospacer+PAM",))
    strand = _column(head, ("Strand",))
    if start is None or aligned is None or strand is None:
        missing = [n for n, v in (("Position/Start_coordinate", start),
                                  ("Aligned_protospacer+PAM", aligned),
                                  ("Strand", strand)) if v is None]
        raise FlankError(
            "the table is neither a CRISPRme integrated results table (no "
            "'Start_coordinate_(%s)') nor the report's curated download (no %s)"
            % (block, ", ".join(missing)))
    return {"chrom": chrom, "start": start, "strand": strand, "aligned": aligned,
            "aligned_ref": None, "rank": _column(head, ("rank",)),
            "guide": _column(head, _GUIDE_COLUMNS),
            "spacer": _column(head, ("Aligned_spacer+PAM",))}, "curated"


def read_rows(path, block=DEFAULT_BLOCK, max_rows=0, row_filter=None):
    """([(row_id, dict)], the shape's name, the header, {filter counts}) of a results
    table: its sites, in file order. ``row_id`` is the 1-based data line number, the key
    the output joins back on. A line whose start or alignment is unreadable is kept with
    a note rather than dropped: a shortlist the user handed us must come back the same
    length.

    ``row_filter`` (a RowFilter) selects which rows are reconstructed. A row it does not
    keep is dropped HERE and counted, so nothing downstream has to know about it; the
    annotate-table writer reads the input again and marks those lines, so the user's own
    table still comes back whole.
    """
    csv.field_size_limit(1 << 30)
    out = []
    counts = {"read": 0, "kept": 0, "filtered_out": 0, "filter_missing": 0,
              "truncated": False, "last_line": 0}
    with open_table(path) as lines:
        reader = table_reader(lines)
        try:
            header = next(reader)
        except StopIteration:
            raise FlankError("%s is empty" % path)
        cols, shape = table_columns(header, block)
        for line_no, row in enumerate(reader, start=1):
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if max_rows and len(out) >= max_rows:
                # stop READING, do not just skip: a genome-wide integrated table has
                # millions of lines and --max-rows is how a user says "a shortlist"
                counts["truncated"] = True
                break
            counts["read"] += 1
            counts["last_line"] = line_no

            def cell(key, _row=row):
                i = cols.get(key)
                return _row[i].strip() if i is not None and i < len(_row) else ""

            if row_filter is not None:
                keep, missing = row_filter.keep(row)
                if missing:
                    counts["filter_missing"] += 1
                if not keep:
                    counts["filtered_out"] += 1
                    continue

            rec = {"chrom": cell("chrom"), "strand": cell("strand"),
                   "rank": cell("rank"),
                   # the guide, without the gap an aligned column writes into it
                   # (see _GUIDE_COLUMNS): the same text on every row of one guide
                   "guide": cell("guide").replace("-", ""), "note": ""}
            aligned = cell("aligned")
            if aligned in _NO_VALUE:
                aligned = cell("aligned_ref")
            rec["aligned"] = aligned.upper()
            # the guide side of the alignment as written. Its gaps are the row's DNA
            # bulges and its run of N's names the PAM's side; both are read per row, so
            # a table mixing alignments (or guides) is still read right.
            spacer = cell("spacer").upper()
            rec["spacer"] = "" if spacer in ("NA", ".") else spacer
            try:
                rec["start"] = int(cell("start"))
            except ValueError:
                rec["start"] = None
                rec["note"] = "unreadable-start"
            if rec["strand"] not in ("+", "-"):
                rec["note"] = ((rec["note"] + ";" if rec["note"] else "")
                               + "unreadable-strand")
            if not rec["aligned"] or set(rec["aligned"].replace("-", "")) - set("ACGTN"):
                rec["note"] = ((rec["note"] + ";" if rec["note"] else "")
                               + "unreadable-alignment")
            out.append((line_no, rec))
    counts["kept"] = len(out)
    return out, shape, header, counts


# --------------------------------------------------------------------------- #
# named window configurations
# --------------------------------------------------------------------------- #
def load_presets(paths=()):
    """{name: {key: value}}: BUILTIN_PRESETS, then the file shipped beside this module,
    then $CRISPRME_FLANK_PRESETS, then each --preset-file, every one overriding the ones
    before it BY NAME.

    An unknown key is refused rather than ignored: a preset that silently does nothing
    is a window the user believes they asked for and did not get, and this command's
    output goes to an oligo order.
    """
    import configparser
    out = {name: dict(body) for name, body in BUILTIN_PRESETS.items()}
    files = [PRESET_FILE] if os.path.isfile(PRESET_FILE) else []
    env = os.environ.get(PRESET_ENV)
    if env:
        for p in env.split(os.pathsep):
            if p.strip():
                files.append(p.strip())
    files.extend(paths)
    for path in files:
        if not os.path.isfile(path):
            raise FlankError("the preset file does not exist: %s" % path)
        parser = configparser.ConfigParser()
        parser.optionxform = str   # keys as written
        try:
            with open(path, errors="replace") as fh:
                parser.read_file(fh, source=path)
        except configparser.Error as exc:
            raise FlankError("%s is not a readable preset file (%s). A preset is an INI "
                             "section: [name] then one key = value per line, from: %s"
                             % (path, exc, ", ".join(PRESET_KEYS)))
        for name in parser.sections():
            body = {}
            for key, value in parser.items(name):
                if key not in PRESET_KEYS:
                    raise FlankError(
                        "preset [%s] of %s sets %r, which is not a preset key. The keys "
                        "are: %s" % (name, path, key, ", ".join(PRESET_KEYS)))
                body[key] = value.strip()
            out[name] = body
    return out


def preset_help(presets):
    """The presets, as --help lists them."""
    lines = []
    for name in sorted(presets):
        body = presets[name]
        how = []
        if body.get("window"):
            how.append("--window %s" % body["window"])
        if body.get("upstream"):
            how.append("--upstream %s" % body["upstream"])
        if body.get("downstream"):
            how.append("--downstream %s" % body["downstream"])
        for key in ("length_policy", "orientation", "max_haplotypes", "filter"):
            if body.get(key):
                how.append("--%s %s" % (key.replace("_", "-"), body[key]))
        lines.append("      %-14s %s" % (name, " ".join(how) or "(nothing)"))
        if body.get("description"):
            import textwrap
            lines.extend(textwrap.wrap(body["description"], 70,
                                       initial_indent=" " * 21,
                                       subsequent_indent=" " * 21))
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# the nuclease: which end of the target is the PAM
# --------------------------------------------------------------------------- #
def read_pam_file(path):
    """(motif, side, length) of a CRISPRme PAM file: one line, ``<sequence> <n>``, where
    the sequence is the guide written as N's with the PAM in place and ``n`` is the PAM's
    length, NEGATIVE when the PAM leads (5', Cas12a) and positive when it trails (3').
    The same reading crisprme.py complete-search does."""
    with open(path, errors="replace") as fh:
        first = fh.readline().split()
    if len(first) < 2:
        raise FlankError("%s is not a CRISPRme PAM file (its first line is not "
                         "'<sequence> <signed PAM length>')" % path)
    try:
        signed = int(first[-1])
    except ValueError:
        raise FlankError("%s is not a CRISPRme PAM file (its first line does not end in "
                         "a signed PAM length)" % path)
    seq = first[0].upper()
    if signed < 0:
        return seq[:-signed], PAM_5PRIME, -signed
    if signed == 0:
        raise FlankError("%s gives a PAM length of 0" % path)
    return seq[-signed:], PAM_3PRIME, signed


def pam_side_from_spacer(spacer, motif_len=None):
    """The PAM's side read off the guide side of a row's alignment, or None.

    The guide is written with an N at every PAM position, so a trailing run of N's is a
    3' PAM and a leading one a 5' PAM -- the same reading ``scripts/check_alt_pam.py``
    does. Gaps are taken out first: a DNA bulge in the PAM region writes one into the
    guide. Returns (side, length), or None where the two ends cannot be told apart (a
    guide with no N's, or one with N's at both ends -- a degenerate spacer).
    """
    s = (spacer or "").replace("-", "").upper()
    if not s:
        return None
    lead = len(s) - len(s.lstrip("N"))
    trail = len(s) - len(s.rstrip("N"))
    if lead >= len(s):       # all N: a PAM-less guide tells us nothing about a side
        return None
    want = motif_len
    if want:
        if trail >= want and lead < want:
            return PAM_3PRIME, want
        if lead >= want and trail < want:
            return PAM_5PRIME, want
        return None
    if trail and not lead:
        return PAM_3PRIME, trail
    if lead and not trail:
        return PAM_5PRIME, lead
    return None


def resolve_pam(resdir, cmd, workdir, pam_file, pam_motif, pam_side, spacers):
    """(side, motif, length, how it was found). Never assumes SpCas9.

    In order: what the user gave; the PAM FILE the run's own command line names (which
    carries the side in the sign of its length, so it is the one unambiguous source);
    the ``Pam`` line of the run's .Params.txt for the motif, with the side read off the
    alignments in the table; the alignments alone. Where none of them answers, the side
    is None and the caller decides whether that matters -- it does not for a symmetric
    window, and run() refuses an asymmetric one rather than guess.
    """
    if pam_side and pam_side not in PAM_SIDES:
        raise FlankError("--pam-side must be one of: %s" % ", ".join(PAM_SIDES))
    if pam_file:
        motif, side, length = read_pam_file(pam_file)
        return side, motif, length, "PAM: %s %s (--pam-file %s)" % (
            motif, side, pam_file)
    if pam_motif and pam_side:
        return pam_side, pam_motif.upper(), len(pam_motif), (
            "PAM: %s %s (--pam-motif/--pam-side)" % (pam_motif.upper(), pam_side))
    motif = pam_motif.upper() if pam_motif else None
    if resdir and not motif:
        named = _flag_value(cmd, "--pam")
        if named:
            path = named if os.path.isabs(named) or workdir is None else os.path.join(
                workdir, named)
            if os.path.isfile(path):
                motif, side, length = read_pam_file(path)
                # an explicit --pam-side still wins: the user may be reading a table
                # the run's own PAM does not describe
                side = pam_side or side
                return side, motif, length, (
                    "PAM: %s %s (%s)" % (motif, side, "--pam-side, with the motif of "
                                         "the run's --pam %s" % named if pam_side else
                                         "the run's --pam %s" % named))
        params = run_params(resdir) if resdir else {}
        if params.get("Pam"):
            motif = params["Pam"].strip().upper()
    if pam_side:
        return pam_side, motif, len(motif) if motif else None, (
            "PAM: %s %s (--pam-side)" % (motif or "motif not recorded", pam_side))
    # the motif's length is known but not its side (or neither): ask the table's own
    # alignments, which carry the guide written with an N at every PAM position
    seen = set()
    for spacer in spacers:
        answer = pam_side_from_spacer(spacer, len(motif) if motif else None)
        if answer:
            seen.add(answer)
    if len(seen) == 1:
        side, length = seen.pop()
        return side, motif or ("N" * length), length, (
            "PAM: %s %s (read off the table's Aligned_spacer+PAM%s)"
            % (motif or "N" * length, side,
               ", its length from the run's Params.txt" if motif else ""))
    if motif:
        return None, motif, len(motif), (
            "PAM: %s, side unknown (the run records the motif but not which end it is "
            "on, and the table's alignments do not say)" % motif)
    return None, None, None, "PAM: not recorded by the run and not readable from the table"


# --------------------------------------------------------------------------- #
# the window
# --------------------------------------------------------------------------- #
def bulge_counts(spacer, protospacer):
    """(DNA bulges, RNA bulges) of a row's alignment, from its two aligned strings.

    CRISPRme writes a bulge as a gap in the side that has no base there: a gap in the
    GUIDE is a DNA bulge (a genomic base with no partner -- the footprint is one base
    LONGER per bulge) and a gap in the PROTOSPACER an RNA bulge (a guide base with no
    genomic partner -- one base SHORTER). Checked against a real run's
    ``Bulge_type`` column, which agrees on all four cases (X, RNA, DNA, RNA,DNA).
    """
    return (spacer or "").count("-"), (protospacer or "").count("-")


def bulge_label(dna, rna, known):
    if not known:
        return "unknown"
    parts = [("DNA:%d" % dna) if dna else "", ("RNA:%d" % rna) if rna else ""]
    return ",".join(p for p in parts if p) or "none"


class Geometry(object):
    """How far beyond each end of the target to reach, and in which direction.

    ``upstream`` is the padding beyond the PROTOSPACER's outer end and ``downstream`` the
    padding beyond the PAM's outer end, both in the guide's orientation. For a 3' PAM
    those are literally "before the protospacer's first base" and "after the PAM's last
    base"; for a 5' PAM the target's two ends swap, so the same two numbers land on the
    opposite genomic sides. ``Geometry.edges`` does that mapping per row.
    """

    __slots__ = ("upstream", "downstream", "pam_side", "policy", "symmetric")

    def __init__(self, upstream, downstream, pam_side=None, policy=POLICY_ANCHORED,
                 symmetric=False):
        self.upstream = upstream
        self.downstream = downstream
        self.pam_side = pam_side
        self.policy = policy
        self.symmetric = symmetric

    @property
    def spec(self):
        """What the output calls this window."""
        if self.symmetric:
            return "window%d" % self.upstream
        return "up%d+down%d" % (self.upstream, self.downstream)

    def edges(self, strand, span, nominal):
        """(left, right, upstream applied, downstream applied, adjustment).

        ``span`` is the reference coordinates the target covers and ``nominal`` the
        spacer+PAM length with its gaps taken out (None where the table does not say).
        They differ exactly where the alignment carries a bulge, or where an indel of
        the haplotype the target was located on sits inside it.

        Under the default ANCHORED policy the paddings are kept and the emitted length
        varies: the bases beyond the target are the bases the user asked for, and the
        sequence is what that piece of genome really is. Under FIXED-TOTAL the reference
        interval is held at upstream + nominal + downstream, and the difference comes
        out of the padding -- beyond the PAM first, then beyond the protospacer, never
        below zero.
        """
        up, down, adjust = self.upstream, self.downstream, 0
        if self.policy == POLICY_FIXED_TOTAL and nominal is not None:
            adjust = nominal - span
            down += adjust
            if down < 0:
                up, down = max(0, up + down), 0
            up = max(0, up)
        five, three = ((down, up) if self.pam_side == PAM_5PRIME else (up, down))
        # the guide's 5' end is at the genomic left on a + row and at the right on a -
        left, right = (five, three) if strand != "-" else (three, five)
        return left, right, up, down, adjust


# --------------------------------------------------------------------------- #
# selecting rows: a general expression over the table's own columns
# --------------------------------------------------------------------------- #
class _Missing(object):
    """A cell with no value where the expression wanted one: not zero, not a string to
    compare against. Any comparison it takes part in is False, and the row is counted."""

    __slots__ = ()

    def __repr__(self):
        return "<missing>"


MISSING = _Missing()
# What a CRISPRme table writes where a number was not measured. '-' and '.' are on this
# list and are ALSO ordinary values: '-' is the minus strand, an intergenic Gene_region,
# a PAM_creation that did not happen; '.' is a VCF's own empty field. Nothing in the file
# says which a given cell is, so _Cell lets the EXPRESSION say: read as a number, one of
# these is missing; compared with a string, a cell is its own text.
_MISSING_CELLS = frozenset(("", "NA", "N/A", "-", ".", "NaN", "nan", "None", "null"))


class _Cell(object):
    """One cell as read, before the expression says how to read it.

    ``CFD > 0.5`` and ``Strand == "-"`` ask different questions of the same characters.
    A score column's '-' means "not measured" and a row chosen on it would be a site put
    into an assay for no reason; a Strand column's '-' is half the table. So the cell
    carries both readings and the comparison picks one: against a string (or a list of
    them, for ``in``) the cell is its TEXT, and anywhere a number is wanted -- a numeric
    comparison, arithmetic, a function -- it is its number, or MISSING where it holds one
    of the markers above. Only that second case counts as a missing value, so
    ``Strand == "-"`` selects the minus-strand rows and reports nothing missing, while
    ``CRISPR_Bulge >= CB_SENSITIVE`` still skips the rows nobody measured and says how
    many there were.
    """

    __slots__ = ("text", "number")

    def __init__(self, text):
        self.text = text
        try:
            self.number = float(text)
        except ValueError:
            self.number = None

    def __repr__(self):                      # pragma: no cover - debugging only
        return "<cell %r>" % self.text

    def as_text(self):
        return self.text

    def as_value(self, state):
        """The cell where a value is wanted and no string says otherwise: its number,
        else MISSING where it holds a marker, else the word itself (which ``_compare``
        then refuses against a number, naming both)."""
        if self.number is not None:
            return self.number
        if self.text in _MISSING_CELLS:
            state["missing"] = True
            return MISSING
        return self.text


def _is_text(value):
    """Does this side of a comparison make it a comparison of TEXT? A string, or a
    list/tuple of them -- the right-hand side of ``Strand in ("+", "-")``."""
    if isinstance(value, str):
        return True
    if isinstance(value, (list, tuple)):
        return bool(value) and all(isinstance(v, str) for v in value)
    return False


def _read_pair(left, right, state):
    """The two sides of one comparison, each read in the light of the other."""
    if isinstance(left, _Cell) and not isinstance(right, _Cell) and _is_text(right):
        return left.as_text(), right
    if isinstance(right, _Cell) and not isinstance(left, _Cell) and _is_text(left):
        return left, right.as_text()
    return _read_value(left, state), _read_value(right, state)


def _read_value(value, state):
    return value.as_value(state) if isinstance(value, _Cell) else value


def column_alias(name):
    """The name a filter expression can write a column by: every run of characters that
    cannot be in a Python name becomes one underscore (``CRISPR-Bulge`` ->
    ``CRISPR_Bulge``, ``CFD_score_(highest_CFD)`` -> ``CFD_score_highest_CFD``). The
    exact name is always reachable as ``col("...")``."""
    return re.sub(r"\W+", "_", name).strip("_")


def cb_sensitive_point():
    """(threshold, recall target) of the CRISPR-Bulge sensitive operating point, from
    PostProcess/cb_operating_point.py -- the single source of truth the report, the risk
    card and the panel selection all key on. Taken BY NAME rather than written as a
    literal, so a filter written today still means "the sensitive threshold" after the
    operating point is re-fitted or CRISPRME_CB_RECALL_TARGET is set."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.append(here)
    try:
        import cb_operating_point as cbop
    except ImportError as exc:   # pragma: no cover - a broken install, named not guessed
        raise FlankError(
            "--filter uses CB_SENSITIVE, the CRISPR-Bulge sensitive operating point, "
            "which lives in cb_operating_point.py beside this command in %s, and it is "
            "not there: write the number itself, or repair the install (%s)"
            % (here, exc))
    return cbop.sensitive_threshold(), cbop.recall_target()


class RowFilter(object):
    """``--filter`` : which rows of the table get flanks, as an expression over the
    table's OWN columns.

    The expression is Python's grammar cut down to what a row test needs -- comparisons,
    and/or/not, arithmetic, ``in`` -- evaluated over a whitelist of names: one per column
    (``column_alias``), ``col("exact name")`` for any column, and the named constants
    below. No attribute access, no indexing, no call but the few named here, so a table
    or an expression cannot run anything.

    Names are checked against the header ONCE, before a single row is read, so a typo is
    an error that names the columns rather than a run that quietly selects nothing.
    """

    FUNCTIONS = {"abs": abs, "min": min, "max": max, "len": len, "int": int,
                 "float": float}
    _NODES = (ast.Expression, ast.BoolOp, ast.UnaryOp, ast.BinOp, ast.Compare,
              ast.Constant, ast.Name, ast.Call, ast.Tuple, ast.List, ast.Load,
              ast.And, ast.Or, ast.Not, ast.USub, ast.UAdd, ast.Add, ast.Sub,
              ast.Mult, ast.Div, ast.Mod, ast.Eq, ast.NotEq, ast.Lt, ast.LtE,
              ast.Gt, ast.GtE, ast.In, ast.NotIn)

    def __init__(self, expression, header):
        self.expression = expression
        self.header = list(header)
        threshold, recall = cb_sensitive_point()
        self.constants = {"CB_SENSITIVE": threshold, "CB_RECALL": recall,
                          "MISSING": MISSING}
        self.index = {}
        clash = set()
        for i, name in enumerate(header):
            alias = column_alias(name)
            if not alias:
                continue
            if alias in self.index and self.index[alias] != i:
                clash.add(alias)
            self.index[alias] = i
        for alias in clash:
            del self.index[alias]
        self.exact = {name: i for i, name in enumerate(header)}
        self.clash = clash
        try:
            self.tree = ast.parse(expression.strip(), mode="eval")
        except SyntaxError as exc:
            raise FlankError("--filter %r is not an expression: %s" % (expression, exc))
        self._check(self.tree)

    # ---- validation, once
    def _check(self, tree):
        # the Name node of ``max(a, b)`` is the function, not a column: told apart by
        # identity, so a column really named max is still read as a column
        called = {id(n.func) for n in ast.walk(tree)
                  if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
        for node in ast.walk(tree):
            if not isinstance(node, self._NODES):
                raise FlankError(
                    "--filter may not use %s. It is a row test over the table's "
                    "columns: comparisons, and/or/not, arithmetic, 'in', and the "
                    "functions %s." % (type(node).__name__,
                                       ", ".join(sorted(self.FUNCTIONS))))
            if isinstance(node, ast.Call):
                if not isinstance(node.func, ast.Name) or node.func.id not in (
                        set(self.FUNCTIONS) | {"col", "has"}):
                    raise FlankError(
                        "--filter may only call col(), has() and %s"
                        % ", ".join(sorted(self.FUNCTIONS)))
                if node.func.id in ("col", "has"):
                    if (len(node.args) != 1
                            or not isinstance(node.args[0], ast.Constant)
                            or not isinstance(node.args[0].value, str)):
                        raise FlankError("--filter: %s() takes one column name in "
                                         "quotes" % node.func.id)
                    self._known_column(node.args[0].value)
            elif isinstance(node, ast.Name):
                if node.id in self.constants or node.id in self.index:
                    continue
                if id(node) in called:
                    continue   # the name of a function being called, checked above
                if node.id in self.FUNCTIONS or node.id in ("col", "has"):
                    raise FlankError("--filter: %s is a function, not a column of this "
                                     "table: write %s(...)" % (node.id, node.id))
                if node.id in self.clash:
                    raise FlankError(
                        "--filter: %s is the name of more than one column of this "
                        "table; write the one you mean as col(\"...\")" % node.id)
                if node.id not in self.index:
                    raise FlankError(
                        "--filter: this table has no column %s. Its columns are: %s. "
                        "A name with a character Python cannot hold is written with "
                        "underscores (CRISPR-Bulge is CRISPR_Bulge) or as "
                        "col(\"CRISPR-Bulge\"). The named constants are: %s."
                        % (node.id, ", ".join(self.header),
                           ", ".join(sorted(self.constants))))

    def _known_column(self, name):
        if name not in self.exact:
            raise FlankError('--filter: this table has no column "%s". Its columns are: '
                             "%s" % (name, ", ".join(self.header)))

    # ---- one row
    def keep(self, row):
        """(keep it?, was a value missing?)

        ONE rule, so there is nothing subtle to get wrong: a row is kept when the
        expression is true AND every cell the expression read AS A NUMBER had one. A
        cell the table does not fill is not zero and not an empty string, and a row
        selected on a value nobody measured is a site put into an assay for no reason;
        a negated test would otherwise let exactly those rows through. A cell compared
        with a STRING is a different question -- ``Strand == "-"`` asks for the text,
        and the minus-strand rows are not rows with no strand -- so there the cell is
        its own characters and nothing is missing (see _Cell). Evaluation
        short-circuits, so ``has("CRISPR-Bulge")`` is the way to ask about presence
        without reading the value: ``not has("X") or X > 1`` keeps a row with no X.
        """
        state = {"missing": False}
        value = _read_value(self._eval(self.tree.body, row, state), state)
        if state["missing"]:
            return False, True
        return bool(value) and value is not MISSING, False

    def _cell(self, i, row, state):
        return _Cell(row[i].strip() if i < len(row) else "")

    def _eval(self, node, row, state):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id in self.constants:
                return self.constants[node.id]
            if node.id not in self.index:
                # a function's name written where a value belongs ("abs > 1"): said as
                # a message, never as a KeyError out of the middle of a run
                raise FlankError("--filter: %s is a function, not a column of this "
                                 "table: write %s(...)" % (node.id, node.id))
            return self._cell(self.index[node.id], row, state)
        if isinstance(node, ast.Call):
            name = node.func.id
            if name == "col":
                return self._cell(self.exact[node.args[0].value], row, state)
            if name == "has":
                # "does this cell hold a measured value": the markers are absence here
                # whatever the column, because that is the question has() is for. To
                # ask about the characters instead, compare with them: Strand == "-"
                i = self.exact[node.args[0].value]
                text = row[i].strip() if i < len(row) else ""
                return text not in _MISSING_CELLS
            args = [_read_value(self._eval(a, row, state), state) for a in node.args]
            if any(a is MISSING for a in args):
                return MISSING
            return self.FUNCTIONS[name](*args)
        if isinstance(node, (ast.Tuple, ast.List)):
            return [_read_value(self._eval(e, row, state), state) for e in node.elts]
        if isinstance(node, ast.BoolOp):
            if isinstance(node.op, ast.And):
                for sub in node.values:
                    v = _read_value(self._eval(sub, row, state), state)
                    if v is MISSING or not v:
                        return False
                return True
            for sub in node.values:
                v = _read_value(self._eval(sub, row, state), state)
                if v is not MISSING and v:
                    return True
            return False
        if isinstance(node, ast.UnaryOp):
            v = _read_value(self._eval(node.operand, row, state), state)
            if isinstance(node.op, ast.Not):
                return v is MISSING or not v
            if v is MISSING:
                return MISSING
            return -v if isinstance(node.op, ast.USub) else +v
        if isinstance(node, ast.BinOp):
            a, b = self._eval(node.left, row, state), self._eval(node.right, row, state)
            a, b = _read_value(a, state), _read_value(b, state)
            if a is MISSING or b is MISSING:
                return MISSING
            return self._arith(node.op, a, b)
        if isinstance(node, ast.Compare):
            left = self._eval(node.left, row, state)
            for op, right_node in zip(node.ops, node.comparators):
                right = self._eval(right_node, row, state)
                # each side read in the light of the other (_read_pair): against a
                # string a cell is its own text, so Strand == "-" works; where a number
                # is wanted a cell the table does not fill is not zero and not a string,
                # so the row does not pass and is counted, and the user is told how many
                # there were rather than left to wonder
                a, b = _read_pair(left, right, state)
                if a is MISSING or b is MISSING:
                    return False
                if not self._compare(op, a, b):
                    return False
                left = right
            return True
        raise FlankError("--filter may not use %s" % type(node).__name__)

    def _arith(self, op, a, b):
        try:
            if isinstance(op, ast.Add):
                return a + b
            if isinstance(op, ast.Sub):
                return a - b
            if isinstance(op, ast.Mult):
                return a * b
            if isinstance(op, ast.Div):
                return a / b
            return a % b
        except (TypeError, ZeroDivisionError) as exc:
            raise FlankError("--filter %r: %s" % (self.expression, exc))

    def _compare(self, op, a, b):
        try:
            if isinstance(op, ast.Eq):
                return a == b
            if isinstance(op, ast.NotEq):
                return a != b
            if isinstance(op, ast.In):
                return a in b
            if isinstance(op, ast.NotIn):
                return a not in b
            if isinstance(op, ast.Lt):
                return a < b
            if isinstance(op, ast.LtE):
                return a <= b
            if isinstance(op, ast.Gt):
                return a > b
            return a >= b
        except TypeError:
            # comparing a word with a number: the expression asks something of a column
            # that does not hold numbers. Said once, with both values, rather than
            # quietly dropping every row of the table
            raise FlankError(
                "--filter %r compares %r with %r, which have no order. Check the column "
                "names: a text column (an origin, a gene) cannot be compared with a "
                "number." % (self.expression, a, b))


# --------------------------------------------------------------------------- #
# the run's own records: where its reference, VCFs and samples list are
# --------------------------------------------------------------------------- #
def run_command_line(resdir):
    """The search's command line (its .command_line.txt), or ''."""
    p = os.path.join(resdir, ".command_line.txt")
    if not os.path.isfile(p):
        return ""
    with open(p, errors="replace") as fh:
        return fh.read()


def _flag_value(command_line, flag):
    m = re.search(r"%s\s+(\S+)" % re.escape(flag), command_line)
    return m.group(1) if m else None


def run_workdir(resdir):
    """The working directory the search ran in: a results folder is <workdir>/Results/<job>,
    and the search's command line names its inputs relative to it."""
    parent = os.path.dirname(os.path.normpath(os.path.abspath(resdir)))
    if os.path.basename(parent) == "Results":
        return os.path.dirname(parent)
    return None


def result_table_of(resdir):
    """The run's integrated results table, or raise naming what was looked for."""
    for pattern in ("*integrated_results.tsv", "*integrated_results.tsv.gz",
                    "*integrated_results.zip"):
        found = sorted(glob.glob(os.path.join(resdir, pattern)))
        if found:
            return found[0]
    raise FlankError(
        "no integrated results table under %s (looked for *integrated_results.tsv, "
        ".tsv.gz and .zip): pass the table with --targets" % resdir)


def _config_entries(path):
    """The non-blank, non-comment lines of a CRISPRme .config.txt list."""
    with open(path, errors="replace") as fh:
        return [line.strip() for line in fh
                if line.strip() and not line.strip().startswith("#")]


def resolve_inputs(resdir, reference, vcfs, vcf_dirs, samples, reference_only=False):
    """(reference path, [VCF paths], [samplesID paths], [how each was found]) for a run.

    What the user gives wins. What the user does not give is read from the run's OWN
    records -- its .command_line.txt (--genome, --vcf, --samplesID, as the search was
    given them, resolved against the working directory a results folder sits under) and
    its .sampleID.txt -- and never guessed: an input the records name but that is not
    there is raised by name, with what to pass instead, as check_indel_rows does.
    """
    how = []
    cmd = run_command_line(resdir) if resdir else ""
    workdir = run_workdir(resdir) if resdir else None

    def under(value):
        """A path the search's command line gives, where it gave it."""
        if os.path.isabs(value) or workdir is None:
            return value
        return os.path.join(workdir, value)

    # ---- the reference genome
    if not reference:
        if not resdir:
            raise FlankError("--reference is needed (no --results folder to read the "
                             "run's own record of its genome from)")
        named = _flag_value(cmd, "--genome")
        if named is None:
            params = run_params(resdir)
            raise FlankError(
                "the run does not record its reference genome (%s has no .command_line.txt "
                "with --genome%s): pass --reference"
                % (resdir, "; its .Params.txt says Genome_selected %s" % params["Genome_selected"]
                   if params.get("Genome_selected") else ""))
        reference = under(named)
        if not os.path.exists(reference):
            raise FlankError(
                "the run searched the genome %r, which resolves to %s and is not there: "
                "pass --reference with the genome the run used"
                % (named, os.path.abspath(reference)))
        how.append("reference: %s (the run's --genome)" % reference)
    else:
        how.append("reference: %s (--reference)" % reference)
    if not os.path.exists(reference):
        raise FlankError("the reference given does not exist: %s" % reference)

    # ---- the VCFs
    paths = list(vcfs)
    if reference_only:
        how.append("VCFs: none read (--reference-only): the reference flank alone")
        vcf_dirs = []
        paths = []
    for d in vcf_dirs:
        if not os.path.isdir(d):
            raise FlankError("--vcf-dir is not a folder: %s" % d)
        found = sorted(glob.glob(os.path.join(d, "*.vcf.gz")))
        if not found:
            raise FlankError("no *.vcf.gz under --vcf-dir %s" % d)
        paths.extend(found)
    if reference_only:
        pass
    elif paths:
        how.append("VCFs: %d given (--vcf / --vcf-dir)" % len(paths))
    elif resdir:
        named = _flag_value(cmd, "--vcf")
        if named is None:
            how.append("VCFs: none -- the run's command line names no --vcf (a "
                       "reference-only search): reference flanks only")
        else:
            listed = under(named)
            if not os.path.isfile(listed):
                raise FlankError(
                    "the run searched the VCF list %r, which resolves to %s and is not "
                    "there: pass --vcf-dir with the dataset folder(s) the run used"
                    % (named, os.path.abspath(listed)))
            missing = []
            for dataset in _config_entries(listed):
                folder = dataset if os.path.isabs(dataset) else (
                    os.path.join(workdir, "VCFs", dataset) if workdir else
                    os.path.join("VCFs", dataset))
                found = sorted(glob.glob(os.path.join(folder, "*.vcf.gz")))
                if not found:
                    missing.append((dataset, os.path.abspath(folder)))
                paths.extend(found)
            if missing:
                raise FlankError(
                    "the run searched the variant dataset(s) %s, whose VCFs are not where "
                    "the run read them (%s: no *.vcf.gz). A run made with --index-path "
                    "does not need its VCFs at search time, so they are often gone. Pass "
                    "--vcf-dir with the dataset folder(s) the run was built from, or "
                    "--reference-only for the reference flank alone."
                    % (", ".join(d for d, _ in missing),
                       "; ".join(p for _, p in missing)))
            how.append("VCFs: %d from the run's --vcf %s" % (len(paths), named))
    else:
        how.append("VCFs: none given: reference flanks only")
    seen = set()
    paths = [p for p in paths
             if not (os.path.realpath(p) in seen or seen.add(os.path.realpath(p)))]
    for p in paths:
        if not os.path.isfile(p):
            raise FlankError("the VCF given does not exist: %s" % p)

    # ---- the samples list
    sample_paths = list(samples)
    if sample_paths:
        how.append("samples: %d file(s) given (--samples)" % len(sample_paths))
    elif resdir:
        own = os.path.join(resdir, ".sampleID.txt")
        if os.path.isfile(own) and read_panel([own]):
            sample_paths = [own]
            how.append("samples: the run's own .sampleID.txt")
        else:
            named = _flag_value(cmd, "--samplesID")
            if named is None:
                how.append("samples: the run records none: every sample of the VCF(s)")
            else:
                listed = under(named)
                if not os.path.isfile(listed):
                    raise FlankError(
                        "the run searched the samples list %r, which resolves to %s and is "
                        "not there: pass --samples" % (named, os.path.abspath(listed)))
                for entry in _config_entries(listed):
                    f = entry if os.path.isabs(entry) else (
                        os.path.join(workdir, "samplesIDs", entry) if workdir else
                        os.path.join("samplesIDs", entry))
                    if not os.path.isfile(f):
                        raise FlankError(
                            "the run searched the samples file %r, which resolves to %s "
                            "and is not there: pass --samples"
                            % (entry, os.path.abspath(f)))
                    sample_paths.append(f)
                how.append("samples: the run's --samplesID %s" % named)
    else:
        how.append("samples: none given: every sample of the VCF(s)")
    for p in sample_paths:
        if not os.path.isfile(p):
            raise FlankError("the samples file given does not exist: %s" % p)
    return reference, paths, sample_paths, how


# --------------------------------------------------------------------------- #
# the VCFs, by chromosome
# --------------------------------------------------------------------------- #
class Variants(object):
    """The VCF(s) a run read, opened once and asked by region.

    Two genotyped VCFs that cover one chromosome and name a sample in common are
    REFUSED: one sample would then have two sets of calls there and no one haplotype, and
    guessing which to write is exactly the thing this command must not do.
    """

    def __init__(self, paths, panel=None):
        self.vcfs = [Vcf(p) for p in paths]
        self.panel = panel
        self._for = {}

    def samples(self):
        out = set()
        for v in self.vcfs:
            out.update(v.samples)
        return out

    def covering(self, chrom):
        """The VCFs that hold ``chrom`` (checked once per chromosome)."""
        if chrom in self._for:
            return self._for[chrom]
        use = [v for v in self.vcfs if v.holds(chrom)]
        seen = {}
        for v in use:
            for s in v.samples:
                if self.panel is not None and s not in self.panel:
                    continue
                if s in seen and seen[s] != v.path:
                    raise FlankError(
                        "two VCFs of %s name the sample in common (%s and %s): one "
                        "sample cannot have two sets of calls at one place, so its "
                        "haplotype is not written. Pass the one dataset the site should "
                        "be read from." % (chrom, os.path.basename(seen[s]),
                                           os.path.basename(v.path)))
                seen[s] = v.path
        self._for[chrom] = use
        return use

    def alleles(self, chrom, beg, end, include_filtered=False):
        """([Allele], {line id: calls}, [non-PASS records in the window]) over
        [beg, end), 0-based.

        A record outside PASS and '.' is left out by default -- that is what CRISPRitz's
        enricher reads, so it is what the search modelled, and the flank then matches the
        row's own alignment -- but it is RETURNED in the third list, so the caller can say
        'there are also these calls in your amplicon window' rather than hide them.

        The third list is the records outside PASS WHETHER OR NOT they were written:
        with ``include_filtered`` their bases are in the sequence the user will order
        oligos against, which is the more important of the two cases to say out loud.
        The written ones also carry their FILTER into the variants column (token below).
        """
        out, calls, filtered = [], {}, []
        for v in self.covering(chrom):
            names = [s if self.panel is None or s in self.panel else None
                     for s in v.samples]
            for r in v.records(chrom, beg, end):
                passing = r["filter"] in PASS_FILTERS
                if not passing:
                    for alt in r["alts"]:
                        if alt and alt not in (".", "*") and not set(alt) - set("ACGTN"):
                            filtered.append("%s_%d_%s_%s|%s" % (chrom, r["pos1"], r["ref"],
                                                                alt, r["filter"]))
                    if not include_filtered:
                        continue
                lid = len(calls)
                here = []
                for k, alt in enumerate(r["alts"], 1):
                    if (not alt or alt in (".", "*") or alt.startswith("<")
                            or set(alt) - set("ACGTN")):
                        continue
                    a = Allele(r["pos1"], r["ref"], alt, lid, k)
                    a.filtered = None if passing else r["filter"]
                    here.append(a)
                if not here:
                    continue
                calls[lid] = parse_calls(r["format"], r["tail"], names)
                out.extend(here)
        out.sort(key=lambda a: (a.lo, a.hi, a.key))
        return out, calls, sorted(set(filtered))


# --------------------------------------------------------------------------- #
# the haplotypes of one site
# --------------------------------------------------------------------------- #
def token_of(allele, chrom):
    """The record as the TSV names it: ``chrom_POS_REF_ALT``, and its FILTER behind a
    '|' for one outside PASS that ``--include-filtered`` wrote into the sequence.

    Without the FILTER the variants column cannot tell a low-quality call from a PASS
    one, and the bases of that call are in the sequence a user orders oligos against.
    The same text the notes' ``filtered-records-*`` list uses, so the two agree.
    """
    return allele.token(chrom) + ("|%s" % allele.filtered if allele.filtered else "")


def slot_of(call, number, h, ploidy):
    """Whether slot ``h`` of a call of ``ploidy`` copies holds allele ``number``:
    'yes' / 'no' / '?'. This is variant_haplotypes.on_slot asked about a slot that IS
    known -- a phased reference call of that ploidy -- rather than about another record's
    slot, so the three answers are its own: a phased call says which slot, a homozygous
    call holds the allele on every slot however it is written, and an unphased
    heterozygous call, a missing allele or a call of another ploidy says '?'.
    """
    return on_slot(((-1,) * max(ploidy, 1), True), h, call, number)


class Entry(object):
    """One emitted sequence of one input row: a flank, who carries it, and where the
    target sits inside it."""

    __slots__ = ("label", "seq", "carriers", "variants", "notes", "offset", "length",
                 "found")

    def __init__(self, label, seq):
        self.label = label
        self.seq = seq
        self.carriers = []      # (sample, slot); empty for a reference nobody differs from
        self.variants = set()
        self.notes = set()
        self.offset = 0
        self.length = 0
        self.found = False


def locate_target(hap, t_lo, t_hi, probe):
    """(offset, length, found) of the target inside a flank haplotype.

    Where the haplotype READS the row's protospacer+PAM at the site (reads_at: the bases
    start at a base of reference coordinate ``t_lo``), that is the answer and ``found`` is
    True -- on a haplotype carrying an insertion inside the target, its exact bases.
    Otherwise the target's reference span [t_lo, t_hi) is mapped onto the haplotype
    through its places: the stretch a reader of this haplotype must look at to see the
    same piece of genome. A coordinate this haplotype deletes maps to the junction, so
    the stretch can be empty -- said by a length of 0, never by a crash.
    """
    if probe:
        here = reads_at(hap, t_lo, probe)
        if here:
            return here[0], len(probe), True
    i0 = hap.span_at(t_lo)[0]
    j1 = hap.span_at(t_hi)[0]
    return i0, max(0, j1 - i0), False


def slot_sets(alleles, calls, chrom):
    """({slot key: [the slot's alleles]}, {slot key: [(sample, slot)]}, [(who, why)]).

    One sample's slot is the alleles its genotype puts on that copy. Samples that come to
    the SAME set of alleles share a key, so the window is written once for all of them:
    in a panel of thousands most slots are the reference, or one common variant.

    A sample whose call in this window is unphased and heterozygous, has a missing
    allele, or has fewer copies than another of its calls has no known slot -- which of
    its alleles share a haplotype is not in the VCF. It is REPORTED with the records at
    fault, and neither of its haplotypes is written: a guess here is a wrong oligo.
    """
    by_line = {}
    for a in alleles:
        by_line.setdefault(a.line, []).append(a)
    lines_of = {}
    ploidy = {}
    for lid, per_sample in calls.items():
        for s, (alle, _phased) in per_sample.items():
            lines_of.setdefault(s, []).append(lid)
            ploidy[s] = max(ploidy.get(s, 0), len(alle))
    sets, members, unresolved = {}, {}, []
    for sample in sorted(lines_of):
        lids = sorted(lines_of[sample])
        unknown = set()
        slots = []
        for h in range(ploidy[sample]):
            mine = []
            for lid in lids:
                call = calls[lid].get(sample)
                if call is None:
                    continue
                for a in by_line.get(lid, ()):
                    answer = slot_of(call, a.number, h, ploidy[sample])
                    if answer == "yes":
                        mine.append(a)
                    elif answer == "?":
                        unknown.add(token_of(a, chrom))
            slots.append(mine)
        if unknown:
            unresolved.append((sample, "unphased-or-missing-genotype:%s"
                               % ",".join(sorted(unknown))))
            continue
        for h, mine in enumerate(slots):
            key = tuple(sorted((a.lo, a.hi, a.key) for a in mine))
            sets.setdefault(key, mine)
            members.setdefault(key, []).append((sample, h))
    return sets, members, unresolved


def target_span(refseq, chrom, row, variants, probe, include_filtered):
    """([t_lo, t_hi), a note or None): the REFERENCE coordinates the row's target covers.

    The row's start is its leftmost base; how many reference coordinates the target
    covers is not the protospacer's length where an insertion inside it adds bases. So
    the span is read off the first sequence that READS the protospacer at the row's
    start -- the reference where the reference reads it (every REF-origin row), else the
    most-carried haplotype of the few records around the target that does. Where nothing
    reads it (a table from another reference, or a REF-fallback alignment of a site that
    exists only on an ALT haplotype) the protospacer's own length stands in, and the note
    says so: every entry of the row then also says protospacer-not-on-this-haplotype.

    A slot whose records overlap has no sequence where they do (make_hap leaves the
    reference standing in a ``clash`` and says so). Only a clash the READING ITSELF
    crosses can make that reading a guess: this span is the row's own target, a couple of
    dozen bases, while the window read here is wider, and a clash elsewhere in it changes
    neither which bases the target covers nor where they are. Dropping the slot for a
    clash anywhere cost the row its protospacer and shifted the target -- and with it
    every flank coordinate the row prints -- by the length of an indel it does not hold.
    """
    start, n = row["start"], len(probe)
    fallback = (start, min(start + max(n, 1), len(refseq)))
    if not n:
        return fallback, "unreadable-alignment"
    # only the records AROUND the target can change which bases it covers, so this reads
    # a window of a few dozen bases rather than the whole flank
    lo = max(0, start - n - 2)
    hi = min(len(refseq), start + 2 * n + 2)
    alleles, calls, _filtered = variants.alleles(chrom, lo, hi, include_filtered)
    sets, members, _unresolved = slot_sets(alleles, calls, chrom)
    order = [[]] + [sets[k] for k in sorted(sets, key=lambda k: (-len(members[k]), k))]
    for mine in order:
        hap = make_hap(refseq, mine, lo, hi)
        for i0 in reads_at(hap, start, probe):
            c0, c1 = span_of(hap.place[i0:i0 + n])
            if hap.clashes(c0, c1):
                continue   # read off bases this slot may not hold: not an answer
            return (c0, min(max(c1 + 1, c0 + 1), len(refseq))), None
    return fallback, "protospacer-not-located"


def site_entries(refseq, chrom, row, variants, geometry, include_filtered,
                 max_haplotypes=0, ref_note=REF_NOTE_CHECKED):
    """([Entry], [(who, reason)], {facts}) for one input row.

    The flank is cut on the REFERENCE -- the same genomic interval for every haplotype of
    a site -- so the entries of a row line up and the TSV can print one interval. An
    indel inside the window then makes one haplotype's flank longer or shorter than
    another's, which is what that individual's genome really holds; ``flank_len`` says so.

    Only samples with a record in the flank are reconstructed: every other sample of the
    panel reads as the reference here, and listing thousands of them would say nothing.

    ``ref_note`` is what the reference entry may claim (REF_NOTE_*): the caller knows
    whether a VCF was read for this contig at all, and the default claim is about every
    sample of the panel.

    ``geometry`` says how far beyond each end of the target to reach. Which genomic side
    each padding lands on is worked out HERE, from the row's own strand and the PAM's
    side, because both can differ from row to row.
    """
    notes = set(row["note"].split(";")) if row["note"] else set()
    probe_bases = row["aligned"].replace("-", "")
    probe = probe_bases if row["strand"] != "-" else rc(probe_bases)

    (t_lo, t_hi), note = target_span(refseq, chrom, row, variants, probe,
                                     include_filtered)
    if note:
        notes.add(note)
    # the spacer+PAM length a window of a fixed total assumes: the guide side of the
    # alignment with its gaps taken out. A bulge makes the genomic footprint differ
    # from it -- longer for a DNA bulge, shorter for an RNA one -- and that is the
    # class a 32-mer assumption gets wrong, so it is said per row either way.
    dna_b, rna_b = bulge_counts(row.get("spacer"), row["aligned"])
    nominal = len(row["spacer"].replace("-", "")) if row.get("spacer") else None
    span = t_hi - t_lo
    pad_left, pad_right, up_applied, down_applied, adjust = geometry.edges(
        row["strand"], span, nominal)
    if nominal is not None and span != nominal:
        notes.add("%s:%d/%d" % (NOTE_FOOTPRINT, span, nominal))
    elif nominal is None and geometry.policy == POLICY_FIXED_TOTAL:
        # fixed-total is measured against the spacer+PAM length, and this table does
        # not carry it: say so rather than quietly behave as if anchored
        notes.add("no-nominal-length(the table has no Aligned_spacer+PAM: "
                  "--length-policy fixed-total has nothing to hold the total to)")
    if adjust:
        notes.add("%s:%+d" % (NOTE_LENGTH_ADJUSTED, adjust))
    flo, fhi = t_lo - pad_left, t_hi + pad_right
    if flo < 0:
        notes.add("truncated-at-contig-start")
        flo = 0
    if fhi > len(refseq):
        notes.add("truncated-at-contig-end")
        fhi = len(refseq)

    alleles, calls, filtered = variants.alleles(chrom, flo, fhi, include_filtered)
    if filtered:
        # WRITTEN: their bases are in these sequences (and their tokens carry the FILTER,
        # token_of). IN-WINDOW: they are not written, as the search did not model them,
        # but a low-quality call in the user's amplicon is still worth a word.
        notes.add(("filtered-records-written:%d" if include_filtered
                   else "filtered-records-in-window:%d") % len(filtered))
    sets, members, unresolved = slot_sets(alleles, calls, chrom)

    def entry_of(mine):
        hap = make_hap(refseq, mine, flo, fhi)
        e = Entry(None, hap.seq)
        e.offset, e.length, e.found = locate_target(hap, t_lo, t_hi, probe)
        return e, hap

    ref_entry, _ref_hap = entry_of([])
    ref_entry.label = "ref"
    entries = {ref_entry.seq: ref_entry}
    for key in sorted(sets, key=lambda k: (-len(members[k]), k)):
        made, hap = entry_of(sets[key])
        here = [z for z in hap.clash if z[1] > flo and z[0] < fhi]
        if here:
            # the VCF describes those bases twice, so THIS slot has no sequence. The
            # sample's other slot may still be known, and is written: a clash takes off
            # the copy it is on, not the individual.
            for sample, h in members[key]:
                unresolved.append(("%s:%d" % (sample, h), "overlapping-records:%s"
                                   % ",".join("%d-%d" % (z[0] + 1, z[1]) for z in here)))
            continue
        e = entries.get(made.seq)
        if e is None:
            entries[made.seq] = e = made
        # an allele is listed for a flank only where it CHANGES a base of it (touches:
        # one whose REF reaches into the window while its edit lies outside writes the
        # reference here, and naming it would read as a difference that is not there)
        here_vars = {token_of(a, chrom) for a in sets[key] if touches(a, flo, fhi - 1)}
        if e.carriers and e.variants != here_vars:
            # two different sets of records that spell one sequence (an indel written at
            # two anchors of a repeat, say). Deduplication is BY SEQUENCE, which is what
            # an oligo order cares about, so they are one entry, and the entry says that
            # its variant list is a union rather than one haplotype's record set.
            e.notes.add("one-sequence-from-more-than-one-variant-set")
        e.carriers.extend(members[key])
        e.variants.update(here_vars)

    alts = [e for e in entries.values() if e is not ref_entry]
    # the sequences that actually HOLD the row's off-target come first (that is what a
    # validation assay is for), each group then by how many haplotypes carry it, and the
    # sequence itself breaks any tie: a deterministic order, most useful first
    alts.sort(key=lambda e: (not e.found, -len(e.carriers), e.seq))
    for i, e in enumerate(alts, 1):
        e.label = "hap%d" % i
    if max_haplotypes and len(alts) > max_haplotypes:
        left = alts[max_haplotypes:]
        notes.add("haplotypes-not-written:%d(%d carrier haplotypes, "
                  "--max-haplotypes 0 writes them all)"
                  % (len(left), sum(len(e.carriers) for e in left)))
        alts = alts[:max_haplotypes]
    facts = {"t_lo": t_lo, "t_hi": t_hi, "flo": flo, "fhi": fhi, "probe": probe,
             "notes": notes, "filtered": filtered, "ref_note": ref_note,
             # the panel samples with a record here: the UNION of the two lists, because
             # a sample can be in both -- one copy written, the other unresolvable (an
             # overlap is reported per slot, while the sample's key stays in members).
             # Adding the lengths counted such a sample twice, and the column then ran
             # past the size of the panel (two rows of a 4,091-sample run printed 4,104
             # and 4,124).
             "n_samples": len({s for v in members.values() for s, _h in v}
                              | {w.split(":")[0] for w, _r in unresolved}),
             "n_haplotypes": len(entries),
             # the geometry of THIS row, so no reader has to infer it from a length a
             # bulge or an indel moved
             "window_spec": geometry.spec, "pad_upstream": up_applied,
             "pad_downstream": down_applied, "pam_side": geometry.pam_side or "unknown",
             "length_policy": geometry.policy, "nominal_target_len": nominal,
             "target_span_ref": span,
             "bulge_type": bulge_label(dna_b, rna_b, bool(row.get("spacer")))}
    return [ref_entry] + alts, sorted(set(unresolved)), facts


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #
def fasta_id(row_id, label):
    """The FASTA name of one entry: the input row's line number and the entry's label.
    Unambiguous (one per emitted sequence) and free of sample identities by construction
    -- carrier names live only in the TSV, which --anonymise blanks."""
    return "crisprme_row%d_%s" % (row_id, label)


def wrap(seq, width=60):
    return [seq[i:i + width] for i in range(0, len(seq), width)] or [""]


def has_sequence(line):
    """Does this output line have bases for a FASTA record? An unresolved sample has no
    sequence by definition, and nor has a haplotype whose window is wholly deleted."""
    return line["haplotype"] != "unresolved" and bool(line.get("_seq"))


def write_outputs(results, out_prefix, anonymise):
    """Write <out_prefix>.fa and <out_prefix>.tsv. Returns (sequences, rows).

    EVERY FASTA RECORD HOLDS BASES. A line with no sequence -- an unresolved sample, or
    a haplotype that deletes the whole window (NOTE_WINDOW_DELETED) -- keeps its TSV
    line, which is where that fact belongs, and writes no record: a header followed by a
    blank line is not a sequence, and the only thing it can do is shift every record
    after it by one for a reader pairing names with sequences. Such a line's ``fasta_id``
    is blank, so nothing in the TSV points at a record that is not there.
    """
    fa_path, tsv_path = out_prefix + ".fa", out_prefix + ".tsv"
    n_seq = 0
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("\t".join(TSV_COLUMNS) + "\n")
        for line in results:
            write_it = has_sequence(line)
            seq = line.pop("_seq", "")
            if write_it:
                n_seq += 1
                fa.write(">%s %s\n" % (line["fasta_id"], " ".join(
                    "%s=%s" % (k, line[k]) for k in
                    ("chrom", "flank_start", "flank_end", "target_offset", "target_len",
                     "target_strand", "target_found", "n_carriers", "row_id",
                     # which way round these bases are written, and in what window:
                     # the FASTA is used on its own, so it must not need the TSV to be
                     # read right
                     "orientation", "emitted_strand", "window_spec"))))
                for piece in wrap(seq):
                    fa.write(piece + "\n")
            if anonymise:
                line["carriers"] = ""
            tsv.write("\t".join(str(line.get(c, "")) for c in TSV_COLUMNS) + "\n")
    return n_seq, len(results)


# --------------------------------------------------------------------------- #
# the column mode: the input table, with the sequences beside the row
# --------------------------------------------------------------------------- #
# The REF and the primary ALT sequence ON THE ROW, which is what turns a results table
# into an oligo order without a join. A row can have SEVERAL distinct ALT sequences, and
# a column is one cell wide: so the count is a column of its own and the sequences stay
# in the companion .fa/.tsv, reachable by flank_row_id -- the --max-haplotypes most useful
# (default 10) of them, every one at --max-haplotypes 0, with flank_alt_count the true
# total either way. Nothing is dropped silently, and nothing is crammed into a cell unless
# --alt-list asks for it.
ANNOTATE_COLUMNS = (
    "flank_row_id",          # the join key into <out>.tsv and <out>.fa
    "flank_window",          # up3+down6 / window200: the geometry these cells were cut with
    "flank_orientation",     # guide | genome
    "flank_strand",          # + | -: the genomic strand the two sequences below read
    "flank_ref",             # the REFERENCE sequence over the window
    "flank_ref_len",
    "flank_ref_offset",      # 0-based offset of the target inside flank_ref
    "flank_ref_target_len",  # bases of flank_ref the target covers
    "flank_alt",             # the PRIMARY ALT sequence (see primary_alt), or SEQ_DELETED
    "flank_alt_len",
    "flank_alt_offset",
    "flank_alt_target_len",
    "flank_alt_count",       # how many DISTINCT ALT sequences this row has, all of them
    "flank_alt_carriers",    # carrier haplotypes of the primary ALT
    "flank_notes",
)
ANNOTATE_ALT_LIST = "flank_alt_all"   # --alt-list: every ALT sequence, delimited
ALT_LIST_SEP = "|"


def annotate_columns(alt_list=False):
    return list(ANNOTATE_COLUMNS) + ([ANNOTATE_ALT_LIST] if alt_list else [])


def annotation_rows(out_lines):
    """{row_id: {column: value}} for the column mode, read off the same lines the TSV
    and the FASTA are written from -- one computation, three outputs, so the three can
    never disagree.

    The PRIMARY ALT is ``hap1``: the first of the row's ALT sequences in the order the
    rest of the command already uses -- the sequences that HOLD the off-target first,
    then by how many haplotypes carry them, then by the sequence itself. That is the one
    a validation assay wants if it only has room for one.

    ``flank_alt_all`` is hap1, hap2, ... in that order, one element per ALT entry
    written: element k IS hapK of ``<prefix>.tsv``, and the first element IS
    ``flank_alt``. A haplotype that deletes the whole window has no sequence, and its
    element is SEQ_DELETED rather than nothing, so that correspondence survives reading
    the cell -- an empty element would be dropped by the obvious split and every later
    sequence would then be read as the wrong haplotype. For the same reason an empty
    ``flank_alt`` means "this row has no ALT entry" and nothing else.
    """
    by_row = defaultdict(list)
    for line in out_lines:
        by_row[line["row_id"]].append(line)
    out = {}
    for row_id, lines in by_row.items():
        ref = next((d for d in lines if d["haplotype"] == "ref"), None)
        alts = [d for d in lines if d["haplotype"].startswith("hap")]
        alts.sort(key=lambda d: int(d["haplotype"][3:]))
        first = alts[0] if alts else None
        notes = ref["notes"] if ref is not None else (
            lines[0]["notes"] if lines else "")
        # n_haplotypes counts the reference, and it is the TRUE total whether or not
        # --max-haplotypes wrote them all: the count must not shrink because the file
        # was kept readable
        total = 0
        if ref is not None and str(ref.get("n_haplotypes", "")).isdigit():
            total = max(0, int(ref["n_haplotypes"]) - 1)
        # a haplotype whose window is wholly deleted: named here too, because the cells
        # below are the only part of the output some readers ever see
        gone = [d["haplotype"] for d in alts if not d.get("_seq")]
        if gone:
            notes = ";".join(filter(None, [notes, "%s:%s" % (NOTE_WINDOW_DELETED,
                                                             ",".join(gone))]))
        cells = {
            "flank_row_id": row_id,
            "flank_window": (ref or first or {}).get("window_spec", ""),
            "flank_orientation": (ref or first or {}).get("orientation", ""),
            "flank_strand": (ref or first or {}).get("emitted_strand", ""),
            "flank_ref": (ref or {}).get("_seq", ""),
            "flank_ref_len": (ref or {}).get("flank_len", ""),
            "flank_ref_offset": (ref or {}).get("target_offset", ""),
            "flank_ref_target_len": (ref or {}).get("target_len", ""),
            "flank_alt": ((first or {}).get("_seq", "") or
                          (SEQ_DELETED if first is not None else "")),
            "flank_alt_len": (first or {}).get("flank_len", ""),
            "flank_alt_offset": (first or {}).get("target_offset", ""),
            "flank_alt_target_len": (first or {}).get("target_len", ""),
            "flank_alt_count": total,
            "flank_alt_carriers": (first or {}).get("n_carriers", ""),
            "flank_notes": notes,
            ANNOTATE_ALT_LIST: ALT_LIST_SEP.join(d["_seq"] or SEQ_DELETED
                                                 for d in alts),
        }
        out[row_id] = cells
    return out


def write_annotated_table(table, path, annot, header, alt_list=False,
                          filter_drop=False, filtered=False, counts=None):
    """Write the input table back out with the sequence columns added, line for line.

    Every original column is copied VERBATIM, in its own order: the user's table comes
    back as their table, with columns on the end. A row the command did not reconstruct
    -- one the filter did not keep, one past --max-rows, one whose coordinates could not
    be read -- keeps its line with the added cells empty and ``flank_notes`` saying
    which, so the table never comes back shorter than it went in (``--filter-drop``
    asks for the shortened form).

    ``counts`` is read_rows' own tally. The join key is the data line's number, counted
    here off the RAW TEXT and there off ``table_reader``, so the two must see the same
    number of data lines; if they ever do not, no table is left behind and the run
    stops, because the failure it would otherwise commit -- one site's sequence on
    another site's row, with every column present and every number plausible -- is not
    one a user can see.
    """
    columns = annotate_columns(alt_list)
    clash = [c for c in columns if c in header]
    if clash:
        raise FlankError(
            "the table already has the column(s) %s, and the column mode would have to "
            "overwrite them: rename them, or write the companion .tsv instead"
            % ", ".join(clash))
    op = gzip.open if path.endswith(".gz") else open
    counts = counts or {}
    last_read = counts.get("last_line", 0)
    written = kept = 0
    with open_table(table) as lines, op(path, "wt") as fh:
        # the ORIGINAL text, not a csv round trip: a cell holding a quote character
        # would come back without it, and the user's table must come back as theirs.
        # The blank-line rule is read_rows', so the line numbers are the same key.
        stream = iter(lines)
        try:
            fh.write(next(stream) + "\t" + "\t".join(columns) + "\n")
        except StopIteration:   # pragma: no cover - run() refuses an empty table
            raise FlankError("%s is empty" % table)
        for line_no, raw in enumerate(stream, start=1):
            row = raw.split("\t")
            if len(row) == 1 and not row[0].strip():
                continue
            written += 1
            cells = annot.get(line_no)
            if cells is None:
                past = last_read and line_no > last_read
                if filter_drop and filtered and not past:
                    continue
                cells = {c: "" for c in columns}
                cells["flank_row_id"] = line_no
                cells["flank_notes"] = ("not-read(--max-rows)" if past
                                        else ("filtered-out" if filtered
                                              else "no-flank-written"))
            kept += 1
            # the line as it was read, then the added columns on the end
            fh.write(raw + "\t" + "\t".join(str(cells.get(c, "")) for c in columns)
                     + "\n")
    expected = counts.get("read")
    if expected is not None and not counts.get("truncated") and written != expected:
        # the join key has gone out of step with the file: every cell added below would
        # be on the wrong line, and nothing in the output would look wrong
        try:
            os.remove(path)
        except OSError:                      # pragma: no cover - whatever left it there
            pass
        raise FlankError(
            "%s has %d data line(s) but was read as %d row(s), so the line number the "
            "sequences are joined on no longer means the same line in both: the "
            "annotated table is NOT written, because every sequence in it could be on "
            "another site's row. Check the table for a stray quote character or an "
            "embedded newline (export it as plain tab-separated text), and report this."
            % (table, written, expected))
    return {"path": path, "rows_in": written, "rows_out": kept, "columns": columns}


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def resolve_window(upstream, downstream, window, preset, presets):
    """(upstream, downstream, symmetric?) from what the user gave and the preset chosen.

    Precedence, highest first: the flag itself, the preset, the built-in default -- the
    rule the --help, the README and flank_presets.conf all state, and it holds between
    the two FORMS of a window as well as within one. ``--window W`` IS ``--upstream W
    --downstream W``, so a preset's ``window = 200`` is 200 at both ends and
    ``--preset amplicon --upstream 5`` is 5 and 200: the flag replaces the end it names
    and nothing else. ``--preset <anything> --window W`` is W at both ends.

    The two forms are refused only where they come from the SAME source and so say two
    different things at once: both on the command line, or both in one preset. Refusing
    the combination of a preset with a flag would make the one preset a user is told to
    write (upstream/downstream, as in flank_presets.conf) unusable with the one flag the
    documentation says always wins.
    """
    body = {}
    if preset:
        if preset not in presets:
            raise FlankError("there is no preset %r. The presets are: %s. Add one with a "
                             "--preset-file (or $%s): a section [name] with any of %s."
                             % (preset, ", ".join(sorted(presets)), PRESET_ENV,
                                ", ".join(PRESET_KEYS)))
        body = presets[preset]

    def of(key):
        if body.get(key):
            try:
                return int(body[key])
            except ValueError:
                raise FlankError("preset [%s] gives %s = %r, which is not a whole "
                                 "number of bases" % (preset, key, body[key]))
        return None

    if window is not None and (upstream is not None or downstream is not None):
        raise FlankError(
            "--window %d and --upstream/--downstream say two different things about one "
            "window: --window W IS --upstream W --downstream W, so give one form or the "
            "other." % window)
    p_up, p_down, p_win = of("upstream"), of("downstream"), of("window")
    if p_win is not None and (p_up is not None or p_down is not None):
        raise FlankError(
            "preset [%s] gives window = %d AND upstream/downstream: they are two ways of "
            "writing one window (window W IS upstream W and downstream W), so the preset "
            "must give one form or the other." % (preset, p_win))
    if p_win is not None:
        # the symmetric form of a preset, read as the two ends it is: a flag can then
        # replace one of them without the preset having to be rewritten
        p_up = p_down = p_win
    if window is not None:
        if window < 0:
            raise FlankError("--window must be 0 or more")
        return window, window, True
    # the symmetric LABEL belongs to a window written in the symmetric form, not to one
    # that happens to have equal ends: a flag naming an end makes it the anchored form
    symmetric = p_win is not None and upstream is None and downstream is None
    up = upstream if upstream is not None else p_up
    down = downstream if downstream is not None else p_down
    up = DEFAULT_UPSTREAM if up is None else up
    down = DEFAULT_DOWNSTREAM if down is None else down
    if up < 0 or down < 0:
        raise FlankError("--upstream and --downstream must be 0 or more")
    return up, down, symmetric


def run(targets=None, results=None, reference=None, vcfs=(), vcf_dirs=(), samples=(),
        out=None, window=None, upstream=None, downstream=None, preset=None,
        preset_files=(), length_policy=None, orientation=None, pam_file=None,
        pam_motif=None, pam_side=None, row_filter=None, annotate_table=None,
        alt_list=False, filter_drop=False, block=DEFAULT_BLOCK, anonymise=False,
        include_filtered=False, max_rows=0, max_haplotypes=None,
        reference_only=False, log=None):
    """Write the flanking sequences of a results table. Returns a summary dict.

    Raises FlankError with the whole message for anything it cannot resolve.
    """
    say = log if log is not None else (lambda _m: None)
    if not out:
        raise FlankError("--out is needed (the prefix of the .fa and .tsv written)")
    if max_haplotypes is not None and max_haplotypes < 0:
        raise FlankError("--max-haplotypes must be 0 (every one) or more")
    if block not in BLOCKS:
        raise FlankError("--block must be one of: %s" % ", ".join(BLOCKS))
    presets = load_presets(preset_files)
    body = presets.get(preset, {}) if preset else {}
    up, down, symmetric = resolve_window(upstream, downstream, window, preset, presets)
    policy = length_policy or body.get("length_policy") or POLICY_ANCHORED
    if policy not in LENGTH_POLICIES:
        raise FlankError("--length-policy must be one of: %s"
                         % ", ".join(LENGTH_POLICIES))
    orientation = orientation or body.get("orientation") or DEFAULT_ORIENTATION
    if orientation not in ORIENTATIONS:
        raise FlankError("--orientation must be one of: %s" % ", ".join(ORIENTATIONS))
    if row_filter is None and body.get("filter"):
        row_filter = body["filter"]
    if max_haplotypes is None:
        # not given: the preset's, else the default. (None rather than the default as
        # the sentinel, so "--max-haplotypes 10" still beats a preset that says 0.)
        if body.get("max_haplotypes"):
            try:
                max_haplotypes = int(body["max_haplotypes"])
            except ValueError:
                raise FlankError("preset [%s] gives max_haplotypes = %r, which is not "
                                 "a whole number" % (preset, body["max_haplotypes"]))
            if max_haplotypes < 0:
                raise FlankError("preset [%s] gives a negative max_haplotypes" % preset)
        else:
            max_haplotypes = DEFAULT_MAX_HAPLOTYPES
    resdir = os.path.abspath(results) if results else None
    if resdir and not os.path.isdir(resdir):
        raise FlankError("the folder given for --results does not exist: %s" % resdir)
    if targets:
        table = os.path.abspath(targets)
        if not os.path.isfile(table):
            raise FlankError("the file given for --targets does not exist: %s" % table)
    elif resdir:
        table = result_table_of(resdir)
    else:
        raise FlankError("give --targets <results table> or --results <results folder>")
    if annotate_table and os.path.abspath(annotate_table) == table:
        raise FlankError(
            "--annotate-table would write over the input table (%s): give it a path of "
            "its own. The input is never modified in place." % table)

    ref_path, vcf_paths, sample_paths, how = resolve_inputs(
        resdir, reference, list(vcfs), list(vcf_dirs), list(samples), reference_only)
    for line in how:
        say("extract-flanks: " + line)

    selector = None
    if row_filter:
        # the header first, so a name the table does not have is an error before a
        # single site is reconstructed rather than an empty output an hour later
        with open_table(table) as lines:
            try:
                header = next(table_reader(lines))
            except StopIteration:
                raise FlankError("%s is empty" % table)
        selector = RowFilter(row_filter, header)
        say("extract-flanks: --filter %s (CB_SENSITIVE = %g, the %d%%-recall "
            "CRISPR-Bulge operating point)"
            % (row_filter, selector.constants["CB_SENSITIVE"],
               selector.constants["CB_RECALL"]))

    rows, shape, header, counts = read_rows(table, block=block, max_rows=max_rows,
                                            row_filter=selector)
    say("extract-flanks: %d row(s) from %s (%s table)" % (len(rows), table, shape))
    if selector is not None:
        say("extract-flanks: the filter kept %d of the %d row(s) read%s"
            % (counts["kept"], counts["read"],
               "; %d had no value for a column it asks about, and a row whose value is "
               "missing is not kept" % counts["filter_missing"]
               if counts["filter_missing"] else ""))
    if not rows:
        if selector is not None and counts["read"]:
            raise FlankError(
                "the filter %r kept none of the %d rows of %s. Check the column names "
                "and the direction of the comparison; the columns are: %s"
                % (row_filter, counts["read"], table, ", ".join(header)))
        raise FlankError("%s has a header but no data rows" % table)

    # which end of the target is the PAM. Read from the run's own records, never
    # assumed: a 5' PAM (Cas12a) puts the two paddings on the opposite genomic sides
    # from a 3' one, so a guess here writes every asymmetric window back to front.
    spacers = []
    for _i, r in rows:                 # the first 200 that have one, not the first 200
        if r.get("spacer"):
            spacers.append(r["spacer"])
            if len(spacers) >= 200:
                break
    side, motif, pam_len, pam_how = resolve_pam(
        resdir, run_command_line(resdir) if resdir else "",
        run_workdir(resdir) if resdir else None, pam_file, pam_motif, pam_side,
        spacers)
    say("extract-flanks: " + pam_how)
    if side is None and up != down:
        raise FlankError(
            "--upstream %d and --downstream %d put different numbers of bases beyond "
            "the two ends of the target, so which end is the PAM decides which genomic "
            "side each lands on -- and %s. Pass --pam-file <the PAM file the search "
            "used> (its signed length carries the side), or --pam-motif <MOTIF> "
            "--pam-side %s, or use the symmetric --window, for which the side does not "
            "matter." % (up, down, pam_how[5:] if pam_how.startswith("PAM: ") else
                         pam_how, "|".join(PAM_SIDES)))
    geometry = Geometry(up, down, side, policy, symmetric)
    say("extract-flanks: window %s, %s length, written in the %s orientation%s"
        % (geometry.spec, policy, orientation,
           " (%s %s PAM)" % (motif, side) if side else ""))

    panel = read_panel(sample_paths) if sample_paths else None
    variants = Variants(vcf_paths, panel)
    if panel is not None and vcf_paths:
        unseen = panel - variants.samples()
        if len(unseen) == len(panel):
            raise FlankError(
                "none of the %d samples of the samples list is in the VCF(s) given: the "
                "list and the VCFs are of different datasets" % len(panel))
        if unseen:
            say("extract-flanks: %d of %d listed samples are in no VCF given"
                % (len(unseen), len(panel)))

    sites_only = [os.path.basename(v.path) for v in variants.vcfs if not v.samples]
    if sites_only:
        # such a dataset has no genotypes, so it names no carrier and no haplotype can be
        # written from it: said once, not left for the user to infer from empty output
        say("extract-flanks: %d VCF(s) have no sample columns (%s): they list no carrier, "
            "so no haplotype is written from them" % (len(sites_only),
                                                      ", ".join(sorted(sites_only)[:3])))
    reference_reader = Reference(ref_path)
    reference_reader.keep_all = False
    out_lines = []
    multi = 0
    # how many DISTINCT ALT sequences each row gave, keyed by row: the number an assay
    # sizes its extra targets from, and the column mode's flank_alt_count
    alt_counts = {}
    by_chrom = defaultdict(list)
    for row_id, row in rows:
        if not row["chrom"]:
            # a line with no chromosome at all is not a site, and it must not be read as
            # a contig the reference is then blamed for lacking: it comes back as a line
            # saying so, like any other row this cannot place (a table exported with a
            # newline inside a cell leaves exactly such a line behind)
            row["note"] = ((row["note"] + ";" if row["note"] else "") + "no-chromosome")
            out_lines.append(_unreadable_line(row_id, row))
            continue
        by_chrom[row["chrom"]].append((row_id, row))
    # Which contigs a VCF was actually read for. A contig no VCF holds -- chrY, an _alt,
    # a _random, a chrUn, 46 of them under the 23-file HGDP+1000G set on one real
    # shortlist -- is NOT a contig where every sample reads as the reference: nothing was
    # read there at all, and the reference entry must say which of the two it is.
    uncovered = []
    for chrom in sorted(by_chrom):
        if not reference_reader.has(chrom):
            raise FlankError(
                "the reference %s has no contig %s, which row(s) %s name: pass the "
                "reference the run searched"
                % (ref_path, chrom, ",".join(str(i) for i, _ in by_chrom[chrom][:5])))
        covering = variants.covering(chrom) if vcf_paths else []
        if not vcf_paths:
            ref_note = REF_NOTE_NO_VCF
        elif not covering:
            ref_note = REF_NOTE_NO_CONTIG
            uncovered.append(chrom)
        elif not any(v.samples for v in covering):
            ref_note = REF_NOTE_NO_GENOTYPES   # a sites-only dataset: no genotype to read
        else:
            ref_note = REF_NOTE_CHECKED
        refseq = reference_reader.seq(chrom)
        for row_id, row in sorted(by_chrom[chrom], key=lambda x: x[0]):
            if row["start"] is None or not row["aligned"]:
                out_lines.append(_unreadable_line(row_id, row))
                continue
            if not 0 <= row["start"] < len(refseq):
                # a coordinate off the end of the contig: another reference, or another
                # assembly's coordinates. Said as a line of its own, never read as a site
                row["note"] = "start-past-the-contig-end(%s is %d bases)" % (
                    chrom, len(refseq))
                out_lines.append(_unreadable_line(row_id, row))
                continue
            entries, unresolved, facts = site_entries(
                refseq, chrom, row, variants, geometry, include_filtered,
                max_haplotypes, ref_note=ref_note)
            if facts["n_haplotypes"] > 1:
                multi += 1
            # n_haplotypes counts the reference entry, so the distinct ALT sequences of
            # the row are one fewer: that is "how many ALT 32-mers does this site have",
            # the number a user adds assay targets for
            alt_counts[row_id] = facts["n_haplotypes"] - 1
            out_lines.extend(_lines_for(row_id, row, chrom, entries, unresolved, facts,
                                        orientation))
    if uncovered:
        say("extract-flanks: %d of the %d contig(s) the table names are in none of the "
            "%d VCF(s) read (%s): no sample was checked on them, so their rows are the "
            "reference sequence alone, noted %s"
            % (len(uncovered), len(by_chrom), len(vcf_paths), ", ".join(uncovered),
               REF_NOTE_NO_CONTIG))
    out_lines.sort(key=lambda d: (d["row_id"], d["_order"]))
    for d in out_lines:
        d.pop("_order", None)
    # read off the lines BEFORE they are written: write_outputs takes the sequences out
    # of them as it writes the FASTA, and the column mode needs the same bytes
    annot = annotation_rows(out_lines) if annotate_table else None
    # entries with a line but no bases: a deletion covers their whole window. Said here
    # because it is why the FASTA holds fewer records than the TSV has entry lines
    deleted = sum(1 for d in out_lines
                  if d["haplotype"] != "unresolved" and not d.get("_seq"))
    n_seq, n_lines = write_outputs(out_lines, out, anonymise)
    multi_alt = sum(1 for n in alt_counts.values() if n > 1)
    say("extract-flanks: %d sequence(s) for %d row(s); %d row(s) have more than one "
        "distinct haplotype (the reference counted), %d have more than one distinct "
        "ALT sequence, AT THE %s WINDOW"
        % (n_seq, len(rows), multi, multi_alt, geometry.spec))
    if deleted:
        say("extract-flanks: %d haplotype(s) have a deletion over the WHOLE window and "
            "so no sequence: their line says %s, and no FASTA record is written for "
            "them" % (deleted, NOTE_WINDOW_DELETED))
    say("extract-flanks: wrote %s.fa and %s.tsv" % (out, out))
    annotated = None
    if annotate_table:
        annotated = write_annotated_table(
            table, annotate_table, annot, header, alt_list=alt_list,
            filter_drop=filter_drop, filtered=(selector is not None), counts=counts)
        say("extract-flanks: wrote %s: the %d line(s) of the input table with %d column(s) "
            "added (%s)" % (annotate_table, annotated["rows_out"],
                            len(annotated["columns"]),
                            ", ".join(annotated["columns"])))
    return {"rows": len(rows), "sequences": n_seq, "lines": n_lines,
            "window_deleted": deleted,
            "multi_haplotype_rows": multi, "multi_alt_rows": multi_alt,
            "alt_counts": alt_counts, "window": geometry.spec,
            "upstream": up, "downstream": down, "length_policy": policy,
            "orientation": orientation, "pam_side": side, "pam_motif": motif,
            "pam_len": pam_len, "filter": row_filter, "filter_counts": counts,
            "annotated_table": annotated,
            "fasta": out + ".fa", "tsv": out + ".tsv",
            "table": table, "shape": shape, "reference": ref_path,
            "vcfs": vcf_paths, "samples": sample_paths}


def _unreadable_line(row_id, row):
    d = {c: "" for c in TSV_COLUMNS}
    d.update({"row_id": row_id, "rank": row["rank"], "guide": row["guide"],
              "chrom": row["chrom"], "target_strand": row["strand"],
              "haplotype": "unresolved", "notes": row["note"] or "unreadable-row",
              "_order": 0, "_seq": ""})
    return d


def _geometry_cells(row, facts, orientation):
    """The columns that say how this row's window was cut, and which way round it is
    written. The same on every entry of a row (a flank's geometry is the row's), and
    written on the unresolved lines too, so a reader never has to go back to the command
    line that made the file."""
    flip = orientation == ORIENT_GUIDE and row["strand"] == "-"
    return {
        "orientation": orientation,
        "emitted_strand": "-" if flip else "+",
        "window_spec": facts["window_spec"],
        "pad_upstream": facts["pad_upstream"],
        "pad_downstream": facts["pad_downstream"],
        "pam_side": facts["pam_side"],
        "length_policy": facts["length_policy"],
        "nominal_target_len": ("" if facts["nominal_target_len"] is None
                               else facts["nominal_target_len"]),
        "target_span_ref": facts["target_span_ref"],
        "bulge_type": facts["bulge_type"],
    }


def _lines_for(row_id, row, chrom, entries, unresolved, facts,
               orientation=DEFAULT_ORIENTATION):
    """The TSV lines (and the sequences the FASTA gets) of one input row.

    ORIENTATION is applied here and nowhere else: everything upstream reconstructs on
    the reference's plus strand, which is where the VCF puts its records, and a
    minus-strand row is turned round once, at the end, with its offset moved to match.
    The coordinates stay genomic whatever the orientation -- a coordinate has only one
    meaning -- and ``emitted_strand`` says which way the bases of that interval were
    written.
    """
    out = []
    n_haps = facts["n_haplotypes"]
    notes_site = sorted(facts["notes"])
    geom = _geometry_cells(row, facts, orientation)
    flip = geom["emitted_strand"] == "-"
    for k, e in enumerate(entries):
        seq, offset = e.seq, e.offset
        if flip:
            seq = rc(seq)
            offset = len(seq) - e.offset - e.length
        target = seq[offset:offset + e.length] if e.length else ""
        notes = sorted(set(notes_site) | e.notes)
        if not e.found:
            notes.append("protospacer-not-on-this-haplotype")
        if not seq:
            # this genome holds NOTHING between the flank's two coordinates: a deletion
            # covers the whole window. flank_len 0 is the only other sign of it, and a
            # length is read as a length, not as a warning
            notes.append(NOTE_WINDOW_DELETED)
        if e.label == "ref":
            # what the reference entry may claim: that every sample reads as the
            # reference here is a statement about the panel, and only a contig a VCF
            # covers was read to make it (REF_NOTE_*)
            notes.append(facts.get("ref_note", REF_NOTE_CHECKED))
        d = {
            "row_id": row_id, "rank": row["rank"], "guide": row["guide"], "chrom": chrom,
            "target_start": facts["t_lo"] + 1, "target_end": facts["t_hi"],
            "target_strand": row["strand"], "haplotype": e.label,
            # blank where no record is written for it (see write_outputs): a name that
            # points at nothing is worse than no name
            "fasta_id": fasta_id(row_id, e.label) if seq else "",
            "flank_start": facts["flo"] + 1, "flank_end": facts["fhi"],
            "flank_len": len(seq), "target_offset": offset, "target_len": e.length,
            "target_found": "yes" if e.found else "no", "target_sequence": target,
            # the target read on the row's own strand, to check against the row's
            # alignment. Under the guide orientation the emitted sequence already is
            # that, so the two columns agree; under the genome orientation they differ
            # on a minus-strand row, which is the point of keeping both.
            "protospacer_on_strand": (target if flip else
                                      (rc(target) if row["strand"] == "-" else target)),
            "n_haplotypes": n_haps, "n_samples": facts["n_samples"],
            "n_carriers": len(e.carriers),
            "carriers": ",".join("%s:%d" % (s, h) for s, h in sorted(e.carriers)),
            "variants": ",".join(sorted(e.variants)),
            "notes": ";".join(notes),
            "_order": k, "_seq": seq,
        }
        d.update(geom)
        out.append(d)
    by_reason = defaultdict(list)
    for sample, reason in unresolved:
        by_reason[reason].append(sample)
    for j, reason in enumerate(sorted(by_reason), start=len(entries)):
        d = {c: "" for c in TSV_COLUMNS}
        d.update({
            "row_id": row_id, "rank": row["rank"], "guide": row["guide"], "chrom": chrom,
            "target_start": facts["t_lo"] + 1, "target_end": facts["t_hi"],
            "target_strand": row["strand"], "haplotype": "unresolved",
            "flank_start": facts["flo"] + 1, "flank_end": facts["fhi"],
            "n_haplotypes": n_haps, "n_samples": facts["n_samples"],
            "n_carriers": len(by_reason[reason]),
            "carriers": ",".join(sorted(by_reason[reason])), "notes": reason,
            "_order": j, "_seq": "",
        })
        d.update(geom)
        out.append(d)
    return out


HELP_TEMPLATE = """crisprme.py extract-flanks --targets <table> | --results <folder>
    --out <prefix>
    [--upstream 3] [--downstream 6] [--window W] [--preset <name>]
    [--length-policy anchored|fixed-total] [--orientation guide|genome]
    [--annotate-table <file.tsv>] [--alt-list] [--filter <expression>] [--filter-drop]
    [--reference <genome dir or FASTA>] [--vcf <bgzipped VCF>] [--vcf-dir <folder>]
    [--samples <samplesID file>] [--pam-file <PAM file>] [--pam-motif M --pam-side S]
    [--block highest_CFD] [--anonymise] [--include-filtered] [--max-rows N]
    [--max-haplotypes 10] [--reference-only] [--preset-file <file>]

Writes the genomic sequence around each off-target site -- the REFERENCE flank and one
flank per distinct observed haplotype, with the samples that carry each -- so a site on a
shortlist can be turned into an oligo order. By default the ten most useful haplotypes of
a row are written and the rest are counted in the notes (--max-haplotypes, default 10;
--max-haplotypes 0 writes every one, and n_haplotypes is always the true total). The
protospacer+PAM is marked in every flank by its offset and length, so an indel that
shifts it is visible rather than miscounted.

THE WINDOW is the target plus a padding at each of its two ends, named for what each
sits beside: --upstream beyond the PROTOSPACER's outer end, --downstream beyond the
PAM's outer end, both read in the GUIDE's orientation. THE DEFAULT IS 3 AND 6 -- for a
20-nt spacer with a 3-nt PAM, a 32-mer -- because an anchored window is the general
form (a symmetric one is the case where the two are equal) and because a short oligo
read straight off the site has no primer slack to absorb a padding error. For
amplicon-scale context use --window 200 (which IS --upstream 200 --downstream 200), or
--preset amplicon. Every output says the window it used and the bases it emitted.

    --targets          a results table: the report's top-N download (curated columns) or
                       an integrated_results.tsv (.tsv, .tsv.gz or .zip)
    --results          a CRISPRme results folder: its own records give the reference,
                       the VCFs, the samples list and the PAM, and -- unless --targets
                       names another table -- its integrated table is the input. Give
                       both to read a shortlist (the report's top-N) with the run's
                       own inputs
    --out              the prefix of the two files written: <prefix>.fa and <prefix>.tsv
    --upstream         bases beyond the protospacer's outer end (default 3)
    --downstream       bases beyond the PAM's outer end (default 6)
    --window           the symmetric form: the same as --upstream W --downstream W.
                       Giving both forms in one place -- both on the command line, or
                       both in one preset -- is refused, not resolved
    --preset           a named window configuration (below). A flag given on the command
                       line always wins over the preset, including over its other form
                       of the window: a preset's window W IS upstream W + downstream W,
                       so --preset <name> --window 200 is 200 at both ends and
                       --preset <name> --upstream 5 replaces that end alone
    --length-policy    what to do where the target's genomic footprint is not the
                       nominal spacer+PAM length -- which is what a bulge means (a DNA
                       bulge is one genomic base MORE, an RNA bulge one LESS):
                         anchored     (default) keep the paddings, let the total vary.
                                      The bases beyond the target are the bases asked
                                      for and the sequence is genomically true
                         fixed-total  hold the reference interval at
                                      upstream+nominal+downstream and take the
                                      difference out of the padding (beyond the PAM
                                      first), for an assay that needs one length
                       Either way every row says its own: flank_len, target_span_ref,
                       nominal_target_len, bulge_type, and the note
                       target-span-not-nominal on the rows affected
    --orientation      guide (default: the protospacer reads 5'->3' as the guide does,
                       which is what an oligo order wants) or genome (always the
                       reference's plus strand). Coordinates are genomic either way;
                       the orientation and emitted_strand columns say which was used
    --annotate-table   also write the INPUT TABLE back out with the sequences added as
                       columns: the reference sequence, the primary ALT sequence, how
                       many distinct ALT sequences the row has, and where the target
                       sits in each. A row with several ALTs keeps them all in the
                       companion <prefix>.tsv/.fa, joined by flank_row_id
    --alt-list         add one more column holding the ALT sequences of the row,
                       '|'-delimited, for a user who wants a single file. It holds the
                       ones WRITTEN: flank_alt_count is always the true total and
                       flank_notes says how many --max-haplotypes left out, so pass
                       --max-haplotypes 0 to have every one in the cell. Element k IS
                       hapK of <prefix>.tsv; a haplotype that deletes the whole window
                       has no sequence and its element is %s, never empty
    --filter           an expression over the TABLE'S OWN columns; only the rows it
                       keeps get flanks. Column names with a character a Python name
                       cannot hold are written with underscores (CRISPR-Bulge is
                       CRISPR_Bulge) or as col("CRISPR-Bulge"). CB_SENSITIVE is the
                       CRISPR-Bulge sensitive operating point, taken from
                       cb_operating_point rather than written as a number, so a filter
                       written today still means the same thing after a re-fit. E.g.
                       both the ALT and the REF alignment above it:
                         --filter 'CRISPR_Bulge >= CB_SENSITIVE and
                                   CRISPR_Bulge_ref >= CB_SENSITIVE'
                       A row whose value for a column the expression read AS A NUMBER
                       is missing (blank, NA, '-', '.') is NOT kept, whichever way
                       round the test is written, and how many there were is said.
                       Compared with a STRING a cell is its own text, because the same
                       characters are data in another column: Strand == "-" selects the
                       minus-strand rows, Gene_region == "-" the intergenic ones, and
                       neither is a missing value. Evaluation short-circuits, so
                       has("col") asks about presence without reading the value:
                       not has("X") or X > 1 keeps a row with no X. has() reads the
                       markers as absence in every column, which is the question it is
                       for; to ask about the characters, compare with them
    --filter-drop      leave the rows the filter did not keep OUT of --annotate-table
                       (by default they keep their line, with the added cells empty and
                       flank_notes filtered-out, so your table comes back whole)
    --reference        the reference genome: a folder of per-chromosome FASTAs, or one
                       FASTA (default: the run's --genome, from its .command_line.txt)
    --vcf              a bgzipped VCF the run searched (repeatable). INDEX THEM FIRST
                       (tabix -p vcf): without a .tbi each file is read end to end once
                       to find its contigs -- minutes per file on a panel-sized VCF,
                       and nothing is printed while it happens
    --vcf-dir          a folder of them (repeatable). Default: the dataset folder(s) the
                       run's --vcf list names, where they are still there
    --samples          a samplesID file (repeatable). Default: the run's .sampleID.txt
    --pam-file         the CRISPRme PAM file the search used: its signed length says
                       which end of the target the PAM is on. Default: the run's own
                       --pam, else the Pam of its Params.txt with the side read off the
                       table's Aligned_spacer+PAM. SpCas9 is never assumed -- an
                       asymmetric window on a run whose PAM side cannot be read is
                       refused rather than guessed
    --pam-motif        the PAM motif (NGG, NRG, TTTV, ...), with --pam-side
    --pam-side         3prime (spacer then PAM: SpCas9) or 5prime (PAM then spacer:
                       Cas12a). On a 5' PAM the two ends of the target swap, so
                       --upstream and --downstream land on the opposite genomic sides
    --block            which alignment of an integrated row to take the site from:
                       highest_CFD (default), fewest_mm+b or highest_CRISPR_BULGE
    --anonymise        write no sample identity anywhere: the TSV's carriers column is
                       blank and only the carrier counts are kept (FASTA names never
                       carry a sample identity)
    --include-filtered also write the records whose FILTER is outside PASS and '.'
                       (by default they are left out, as the search left them out, and
                       counted in the notes: filtered-records-in-window:N). Written,
                       they are marked: the note says filtered-records-written:N and
                       each such record carries its FILTER in the variants column
    --max-rows         stop after N rows of the table (a genome-wide integrated table has
                       millions; this command is meant for a shortlist)
    --max-haplotypes   how many haplotypes to write per row beside the reference
                       (default 10, the most useful first: the ones that hold the
                       off-target, then by carrier count; 0 writes every one). The
                       n_haplotypes column always gives the true total, and the notes
                       say how many were left out
    --reference-only   do not read any VCF: the reference flank alone
    --preset-file      an INI file of presets (repeatable), also read from
                       $%s. A section [name] with any of
                       %s

PRESETS (--preset <name>):
%s
""" % (SEQ_DELETED, PRESET_ENV, "\n                       ".join(
    __import__("textwrap").wrap(", ".join(PRESET_KEYS), 55)), "%s")


def help_text():
    """The help, with whatever presets this install actually has."""
    try:
        presets = load_presets()
    except FlankError as exc:
        return HELP_TEMPLATE % ("      (the preset file could not be read: %s)" % exc)
    return HELP_TEMPLATE % preset_help(presets)


# the module-level name the CLI and crisprme.py have always written to stderr
HELP = help_text()


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--targets")
    p.add_argument("--results")
    p.add_argument("--out")
    p.add_argument("--reference")
    p.add_argument("--vcf", action="append", default=[])
    p.add_argument("--vcf-dir", dest="vcf_dir", action="append", default=[])
    p.add_argument("--samples", action="append", default=[])
    # None, not a default: "not given" has to be told from "given the default", so a
    # preset can supply the value and the flag can still win over it
    p.add_argument("--window", type=int, default=None)
    p.add_argument("--upstream", type=int, default=None)
    p.add_argument("--downstream", type=int, default=None)
    p.add_argument("--preset")
    p.add_argument("--preset-file", dest="preset_file", action="append", default=[])
    p.add_argument("--length-policy", dest="length_policy")
    p.add_argument("--orientation")
    p.add_argument("--annotate-table", dest="annotate_table")
    p.add_argument("--alt-list", dest="alt_list", action="store_true")
    p.add_argument("--filter", dest="row_filter")
    p.add_argument("--filter-drop", dest="filter_drop", action="store_true")
    p.add_argument("--pam-file", dest="pam_file")
    p.add_argument("--pam-motif", dest="pam_motif")
    p.add_argument("--pam-side", dest="pam_side")
    p.add_argument("--block", default=DEFAULT_BLOCK)
    p.add_argument("--anonymise", action="store_true")
    p.add_argument("--anonymize", dest="anonymise", action="store_true")
    p.add_argument("--include-filtered", dest="include_filtered", action="store_true")
    p.add_argument("--max-rows", dest="max_rows", type=int, default=0)
    p.add_argument("--max-haplotypes", dest="max_haplotypes", type=int, default=None)
    p.add_argument("--reference-only", dest="reference_only", action="store_true")
    p.add_argument("--help", "-h", action="store_true")
    args = p.parse_args(argv)
    if args.help or (not args.targets and not args.results):
        sys.stderr.write(help_text())
        return 1
    try:
        run(targets=args.targets, results=args.results, reference=args.reference,
            vcfs=args.vcf, vcf_dirs=args.vcf_dir, samples=args.samples, out=args.out,
            window=args.window, upstream=args.upstream, downstream=args.downstream,
            preset=args.preset, preset_files=args.preset_file,
            length_policy=args.length_policy, orientation=args.orientation,
            annotate_table=args.annotate_table, alt_list=args.alt_list,
            row_filter=args.row_filter, filter_drop=args.filter_drop,
            pam_file=args.pam_file, pam_motif=args.pam_motif, pam_side=args.pam_side,
            block=args.block, anonymise=args.anonymise,
            include_filtered=args.include_filtered, max_rows=args.max_rows,
            max_haplotypes=args.max_haplotypes,
            reference_only=args.reference_only, log=lambda m: print(m))
    except FlankError as e:
        sys.stderr.write("ERROR: %s\n" % e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

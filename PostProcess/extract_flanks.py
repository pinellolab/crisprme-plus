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

OUTPUT. ``<out>.fa`` (one record per emitted sequence, no sample identities in any name)
and ``<out>.tsv`` (one line per emitted sequence, joining back to the input row by its
line number and, where the input has one, its rank). Deterministic: rows in input order,
the reference entry first, then the haplotypes that HOLD the off-target, then by carrier
count and by the sequence itself; carriers and variants sorted.

COORDINATES. CRISPRme's ``Start_coordinate_*`` / curated ``Position`` is the 0-based
reference coordinate of the target's leftmost base, on either strand. Flank coordinates in
the TSV are printed 1-based inclusive (as a genome browser reads them); the target offset
inside a flank is 0-based, as a Python slice reads it. Every flank is written on the PLUS
strand of the reference, whatever the row's strand: one orientation, stated once, is what
an oligo order needs. The TSV's ``protospacer_on_strand`` gives the target read on the
row's own strand, for checking against the row's alignment.
"""

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


# 200 bases either side of the target. An amplicon for validation sequencing (rhAmpSeq,
# CRISPResso2 and the like) is typically 200-300 bp with its primers 50-150 bp away from
# the cut, and a Sanger/TIDE amplicon 400-700 bp: 200 + target + 200 covers the first
# comfortably and the short end of the second, while staying small enough that a top-1000
# list is a file you can open. Change it with --window.
DEFAULT_WINDOW = 200
# how many haplotypes of a site are written beside the reference. A 400 bp window over a
# panel of thousands holds dozens of distinct haplotypes (on the 1000-site shortlist of a
# genome-wide HGDP+1000G run: a median of 23, and one site with 4264), nearly all differing from
# each other far from the target. Writing them all turns a shortlist into megabytes, so
# the most useful come first -- the sequences that HOLD the off-target, then by how many
# haplotypes carry them -- and the rest are counted in the notes rather than written.
# --max-haplotypes 0 writes every one. n_haplotypes always gives the true total.
DEFAULT_MAX_HAPLOTYPES = 10
# the alignment blocks an integrated results table carries, and the one taken by default
BLOCKS = ("highest_CFD", "fewest_mm+b", "highest_CRISPR_BULGE")
DEFAULT_BLOCK = "highest_CFD"
_NO_VALUE = ("", "NA", "n", "-", ".")

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
    "flank_start",       # 1-based inclusive, on the reference
    "flank_end",         # 1-based inclusive, on the reference
    "flank_len",         # the bases of THIS haplotype between them (an indel moves it)
    "target_offset",     # 0-based, into the flank sequence of this entry
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
    "notes",
)


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
                yield (line.decode("utf-8", "replace").rstrip("\n") for line in fh)
        return
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        yield (line.rstrip("\n") for line in fh)


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
            "guide": _column(head, ("Spacer+PAM", "crRNA")),
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
            "guide": _column(head, ("Spacer+PAM", "crRNA"))}, "curated"


def read_rows(path, block=DEFAULT_BLOCK, max_rows=0):
    """[(row_id, dict)] of a results table: its sites, in file order. ``row_id`` is the
    1-based data line number, the key the output joins back on. A line whose start or
    alignment is unreadable is kept with a note rather than dropped: a shortlist the user
    handed us must come back the same length."""
    csv.field_size_limit(1 << 30)
    out = []
    with open_table(path) as lines:
        reader = csv.reader(lines, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise FlankError("%s is empty" % path)
        cols, shape = table_columns(header, block)
        for line_no, row in enumerate(reader, start=1):
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if max_rows and len(out) >= max_rows:
                break

            def cell(key, _row=row):
                i = cols.get(key)
                return _row[i].strip() if i is not None and i < len(_row) else ""

            rec = {"chrom": cell("chrom"), "strand": cell("strand"),
                   "rank": cell("rank"), "guide": cell("guide"), "note": ""}
            aligned = cell("aligned")
            if aligned in _NO_VALUE:
                aligned = cell("aligned_ref")
            rec["aligned"] = aligned.upper()
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
    return out, shape


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
        """([Allele], {line id: calls}, [filtered-out records]) over [beg, end), 0-based.

        A record outside PASS and '.' is left out by default -- that is what CRISPRitz's
        enricher reads, so it is what the search modelled, and the flank then matches the
        row's own alignment -- but it is RETURNED in the third list, so the caller can say
        'there are also these calls in your amplicon window' rather than hide them.
        """
        out, calls, filtered = [], {}, []
        for v in self.covering(chrom):
            names = [s if self.panel is None or s in self.panel else None
                     for s in v.samples]
            for r in v.records(chrom, beg, end):
                passing = r["filter"] in PASS_FILTERS
                if not passing and not include_filtered:
                    for alt in r["alts"]:
                        if alt and alt not in (".", "*") and not set(alt) - set("ACGTN"):
                            filtered.append("%s_%d_%s_%s|%s" % (chrom, r["pos1"], r["ref"],
                                                                alt, r["filter"]))
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
                        unknown.add(a.token(chrom))
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
        if hap.clash:
            continue
        here = reads_at(hap, start, probe)
        if here:
            c0, c1 = span_of(hap.place[here[0]:here[0] + n])
            return (c0, min(max(c1 + 1, c0 + 1), len(refseq))), None
    return fallback, "protospacer-not-located"


def site_entries(refseq, chrom, row, variants, window, include_filtered,
                 max_haplotypes=0):
    """([Entry], [(who, reason)], {facts}) for one input row.

    The flank is cut on the REFERENCE -- the same genomic interval for every haplotype of
    a site -- so the entries of a row line up and the TSV can print one interval. An
    indel inside the window then makes one haplotype's flank longer or shorter than
    another's, which is what that individual's genome really holds; ``flank_len`` says so.

    Only samples with a record in the flank are reconstructed: every other sample of the
    panel reads as the reference here, and listing thousands of them would say nothing.
    """
    notes = set(row["note"].split(";")) if row["note"] else set()
    probe_bases = row["aligned"].replace("-", "")
    probe = probe_bases if row["strand"] != "-" else rc(probe_bases)

    (t_lo, t_hi), note = target_span(refseq, chrom, row, variants, probe,
                                     include_filtered)
    if note:
        notes.add(note)
    flo, fhi = t_lo - window, t_hi + window
    if flo < 0:
        notes.add("truncated-at-contig-start")
        flo = 0
    if fhi > len(refseq):
        notes.add("truncated-at-contig-end")
        fhi = len(refseq)

    alleles, calls, filtered = variants.alleles(chrom, flo, fhi, include_filtered)
    if filtered:
        notes.add("filtered-records-in-window:%d" % len(filtered))
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
        here_vars = {a.token(chrom) for a in sets[key] if touches(a, flo, fhi - 1)}
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
             "notes": notes, "filtered": filtered,
             "n_samples": len({s for v in members.values() for s, _h in v})
                          + len({w.split(":")[0] for w, _r in unresolved}),
             "n_haplotypes": len(entries)}
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


def write_outputs(results, out_prefix, anonymise):
    """Write <out_prefix>.fa and <out_prefix>.tsv. Returns (sequences, rows)."""
    fa_path, tsv_path = out_prefix + ".fa", out_prefix + ".tsv"
    n_seq = 0
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("\t".join(TSV_COLUMNS) + "\n")
        for line in results:
            if line["haplotype"] != "unresolved":
                n_seq += 1
                fa.write(">%s %s\n" % (line["fasta_id"], " ".join(
                    "%s=%s" % (k, line[k]) for k in
                    ("chrom", "flank_start", "flank_end", "target_offset", "target_len",
                     "target_strand", "target_found", "n_carriers", "row_id"))))
                for piece in wrap(line.pop("_seq")):
                    fa.write(piece + "\n")
            else:
                line.pop("_seq", None)
            if anonymise:
                line["carriers"] = ""
            tsv.write("\t".join(str(line.get(c, "")) for c in TSV_COLUMNS) + "\n")
    return n_seq, len(results)


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #
def run(targets=None, results=None, reference=None, vcfs=(), vcf_dirs=(), samples=(),
        out=None, window=DEFAULT_WINDOW, block=DEFAULT_BLOCK, anonymise=False,
        include_filtered=False, max_rows=0, max_haplotypes=DEFAULT_MAX_HAPLOTYPES,
        reference_only=False, log=None):
    """Write the flanking sequences of a results table. Returns a summary dict.

    Raises FlankError with the whole message for anything it cannot resolve.
    """
    say = log if log is not None else (lambda _m: None)
    if not out:
        raise FlankError("--out is needed (the prefix of the .fa and .tsv written)")
    if window < 0:
        raise FlankError("--window must be 0 or more")
    if max_haplotypes < 0:
        raise FlankError("--max-haplotypes must be 0 (every one) or more")
    if block not in BLOCKS:
        raise FlankError("--block must be one of: %s" % ", ".join(BLOCKS))
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

    ref_path, vcf_paths, sample_paths, how = resolve_inputs(
        resdir, reference, list(vcfs), list(vcf_dirs), list(samples), reference_only)
    for line in how:
        say("extract-flanks: " + line)

    rows, shape = read_rows(table, block=block, max_rows=max_rows)
    say("extract-flanks: %d row(s) from %s (%s table)" % (len(rows), table, shape))
    if not rows:
        raise FlankError("%s has a header but no data rows" % table)

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
    by_chrom = defaultdict(list)
    for row_id, row in rows:
        by_chrom[row["chrom"]].append((row_id, row))
    for chrom in sorted(by_chrom):
        if not reference_reader.has(chrom):
            raise FlankError(
                "the reference %s has no contig %s, which row(s) %s name: pass the "
                "reference the run searched"
                % (ref_path, chrom, ",".join(str(i) for i, _ in by_chrom[chrom][:5])))
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
                refseq, chrom, row, variants, window, include_filtered, max_haplotypes)
            if facts["n_haplotypes"] > 1:
                multi += 1
            out_lines.extend(_lines_for(row_id, row, chrom, entries, unresolved, facts))
    out_lines.sort(key=lambda d: (d["row_id"], d["_order"]))
    for d in out_lines:
        d.pop("_order", None)
    n_seq, n_lines = write_outputs(out_lines, out, anonymise)
    say("extract-flanks: %d sequence(s) for %d row(s); %d row(s) have more than one "
        "distinct haplotype" % (n_seq, len(rows), multi))
    say("extract-flanks: wrote %s.fa and %s.tsv" % (out, out))
    return {"rows": len(rows), "sequences": n_seq, "lines": n_lines,
            "multi_haplotype_rows": multi, "fasta": out + ".fa", "tsv": out + ".tsv",
            "table": table, "shape": shape, "reference": ref_path,
            "vcfs": vcf_paths, "samples": sample_paths}


def _unreadable_line(row_id, row):
    d = {c: "" for c in TSV_COLUMNS}
    d.update({"row_id": row_id, "rank": row["rank"], "guide": row["guide"],
              "chrom": row["chrom"], "target_strand": row["strand"],
              "haplotype": "unresolved", "notes": row["note"] or "unreadable-row",
              "_order": 0, "_seq": ""})
    return d


def _lines_for(row_id, row, chrom, entries, unresolved, facts):
    out = []
    n_haps = facts["n_haplotypes"]
    notes_site = sorted(facts["notes"])
    for k, e in enumerate(entries):
        seq = e.seq
        target = seq[e.offset:e.offset + e.length] if e.length else ""
        notes = sorted(set(notes_site) | e.notes)
        if not e.found:
            notes.append("protospacer-not-on-this-haplotype")
        if e.label == "ref":
            notes.append("reference-of-every-sample-without-a-variant-here")
        d = {
            "row_id": row_id, "rank": row["rank"], "guide": row["guide"], "chrom": chrom,
            "target_start": facts["t_lo"] + 1, "target_end": facts["t_hi"],
            "target_strand": row["strand"], "haplotype": e.label,
            "fasta_id": fasta_id(row_id, e.label),
            "flank_start": facts["flo"] + 1, "flank_end": facts["fhi"],
            "flank_len": len(seq), "target_offset": e.offset, "target_len": e.length,
            "target_found": "yes" if e.found else "no", "target_sequence": target,
            "protospacer_on_strand": (rc(target) if row["strand"] == "-" else target),
            "n_haplotypes": n_haps, "n_samples": facts["n_samples"],
            "n_carriers": len(e.carriers),
            "carriers": ",".join("%s:%d" % (s, h) for s, h in sorted(e.carriers)),
            "variants": ",".join(sorted(e.variants)),
            "notes": ";".join(notes),
            "_order": k, "_seq": seq,
        }
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
        out.append(d)
    return out


HELP = """crisprme.py extract-flanks --targets <table> | --results <folder> --out <prefix>
    [--reference <genome dir or FASTA>] [--vcf <bgzipped VCF>] [--vcf-dir <folder>]
    [--samples <samplesID file>] [--window 200] [--block highest_CFD]
    [--anonymise] [--include-filtered] [--max-rows N] [--max-haplotypes 10]
    [--reference-only]

Writes the genomic sequence around each off-target site -- the REFERENCE flank and one
flank per distinct observed haplotype, with the samples that carry each -- so a site on a
shortlist can be turned into an oligo order. The protospacer+PAM is marked in every flank
by its offset and length, so an indel that shifts it is visible rather than miscounted.

    --targets          a results table: the report's top-N download (curated columns) or
                       an integrated_results.tsv (.tsv, .tsv.gz or .zip)
    --results          a CRISPRme results folder: its own records give the reference,
                       the VCFs and the samples list, and -- unless --targets names
                       another table -- its integrated table is the input. Give both to
                       read a shortlist (the report's top-N) with the run's own inputs
    --out              the prefix of the two files written: <prefix>.fa and <prefix>.tsv
    --reference        the reference genome: a folder of per-chromosome FASTAs, or one
                       FASTA (default: the run's --genome, from its .command_line.txt)
    --vcf              a bgzipped VCF the run searched (repeatable)
    --vcf-dir          a folder of them (repeatable). Default: the dataset folder(s) the
                       run's --vcf list names, where they are still there
    --samples          a samplesID file (repeatable). Default: the run's .sampleID.txt
    --window           bases of flank on either side of the target (default 200: an
                       amplicon of ~250-300 bp with room to place primers 50-150 bp out)
    --block            which alignment of an integrated row to take the site from:
                       highest_CFD (default), fewest_mm+b or highest_CRISPR_BULGE
    --anonymise        write no sample identity anywhere: the TSV's carriers column is
                       blank and only the carrier counts are kept (FASTA names never
                       carry a sample identity)
    --include-filtered also write the records whose FILTER is outside PASS and '.'
                       (by default they are left out, as the search left them out, and
                       counted in the notes column)
    --max-rows         stop after N rows of the table (a genome-wide integrated table has
                       millions; this command is meant for a shortlist)
    --max-haplotypes   how many haplotypes to write per row beside the reference
                       (default 10, the most useful first: the ones that hold the
                       off-target, then by carrier count; 0 writes every one). The
                       n_haplotypes column always gives the true total, and the notes
                       say how many were left out
    --reference-only   do not read any VCF: the reference flank alone
"""


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
    p.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    p.add_argument("--block", default=DEFAULT_BLOCK)
    p.add_argument("--anonymise", action="store_true")
    p.add_argument("--anonymize", dest="anonymise", action="store_true")
    p.add_argument("--include-filtered", dest="include_filtered", action="store_true")
    p.add_argument("--max-rows", dest="max_rows", type=int, default=0)
    p.add_argument("--max-haplotypes", dest="max_haplotypes", type=int,
                   default=DEFAULT_MAX_HAPLOTYPES)
    p.add_argument("--reference-only", dest="reference_only", action="store_true")
    p.add_argument("--help", "-h", action="store_true")
    args = p.parse_args(argv)
    if args.help or (not args.targets and not args.results):
        sys.stderr.write(HELP)
        return 1
    try:
        run(targets=args.targets, results=args.results, reference=args.reference,
            vcfs=args.vcf, vcf_dirs=args.vcf_dir, samples=args.samples, out=args.out,
            window=args.window, block=args.block, anonymise=args.anonymise,
            include_filtered=args.include_filtered, max_rows=args.max_rows,
            max_haplotypes=args.max_haplotypes,
            reference_only=args.reference_only, log=lambda m: print(m))
    except FlankError as e:
        sys.stderr.write("ERROR: %s\n" % e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

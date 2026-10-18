#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""targets-scoring: score a user-supplied LIST of candidate off-target sites.

Takes an externally-produced list of off-target sites -- e.g. CasOFFinder output,
or any TSV carrying the aligned guide/off-target columns -- and annotates each with
the SAME calculations a normal CRISPRme+ report uses, with the genome search
replaced by the supplied list:

  * CFD score            -- the canonical CRISPRme CFD (``calc_cfd`` below is kept
                            byte-identical to new_simple_analysis.calc_cfd /
                            analisi_indels_NNN.calc_cfd; pinned by
                            test_calc_cfd_twins_agree.py).
  * CRISPR-Bulge score   -- the CRISPRme+ ML scorer (dedicated ``cbulge`` env);
                            -1.0 when the env is absent (never fails the run).
  * carrier overlay      -- OPTIONAL, when an installed variant index is given
                            (--index-path): for each site, the variants overlapping
                            its protospacer+PAM window and, per population dataset,
                            their allele frequency + carrier count ("in how many
                            people the site is polymorphic").

No index is needed for pure scoring; the carrier overlay reads the index's
``Dictionaries/registry_<vcf>/`` tier (the same registry the report uses).

Scope note (v1): mismatch-only sites (the common case, and CasOFFinder's default
NNN scan) are scored exactly. Bulge rows are supported when the aligned
guide/off-target follow CRISPRme's alignment convention (bulge columns as '-',
DNA-bulge offset = bulge size). CRISPR-Bulge natively scores <=1-bulge alignments;
rows it cannot represent (>1 gap column, N in the off-target, length mismatch,
length > 24) are nulled to -1.0 and carried by CFD + edit distance -- identical to
the main pipeline's behavior.
"""

import argparse
import glob
import os
import pickle
import sys
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# CFD scorer -- VERBATIM copy of the canonical calc_cfd/revcom/get_mm_pam_scores
# (new_simple_analysis.py / analisi_indels_NNN.py). Those two modules parse
# sys.argv at import time, so they cannot be imported here; the codebase already
# keeps byte-identical copies guarded by test_calc_cfd_twins_agree.py. That test
# is extended to pin THIS copy to the canonical one as well. Do not edit without
# updating the twin (the test will fail otherwise).
# ---------------------------------------------------------------------------


def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(
            open(
                os.path.dirname(os.path.realpath(__file__)) + "/mismatch_score.pkl",
                "rb",
            )
        )
        pam_scores = pickle.load(
            open(os.path.dirname(os.path.realpath(__file__)) + "/PAM_scores.pkl", "rb")
        )
        return (mm_scores, pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")


def revcom(s):
    basecomp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}
    letters = list(s[::-1])
    try:
        letters = [basecomp[base] for base in letters]
    except:
        return None  # If some IUPAC were not translated
    return "".join(letters)


def calc_cfd(guide_seq, sg, pam, mm_scores, pam_scores, do_scores):
    score = 1
    sg = sg.replace("T", "U")
    guide_seq = guide_seq.replace("T", "U")
    s_list = list(sg)
    guide_seq_list = list(guide_seq)

    for i, sl in enumerate(s_list):
        if guide_seq_list[i] == sl:
            score *= 1
        else:
            try:  # Catch exception if IUPAC character
                key = "r" + guide_seq_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            except Exception as e:
                score = 0
                break
            try:
                score *= mm_scores[key]
            except (
                Exception
            ) as e:  # If '-' is in first position, i do not have the score for that position
                pass

    # Guard against non-ATCG PAM bases (issue #94): a non-canonical PAM contributes
    # 0 to the CFD product instead of raising KeyError.
    score *= pam_scores.get(pam, 0.0)
    return score


# ---------------------------------------------------------------------------
# Site model + scoring
# ---------------------------------------------------------------------------


class Site(object):
    """One candidate off-target site in CRISPRme's aligned representation.

    crrna / dna are the aligned guide / off-target strings INCLUDING the PAM, same
    length, with '-' for bulge columns (RNA bulge -> '-' in dna, DNA bulge -> '-'
    in crrna). ``pam_len`` is the PAM length in bp (3 for Sp/NRG).
    """

    __slots__ = ("crrna", "dna", "chrom", "pos", "strand", "bulge_type",
                 "bulge_size", "mismatches", "pam_len", "extra")

    def __init__(self, crrna, dna, chrom=None, pos=None, strand=None,
                 bulge_type="X0", bulge_size=0, mismatches=None, pam_len=3,
                 extra=None):
        self.crrna = (crrna or "").upper()
        self.dna = (dna or "").upper()
        self.chrom = chrom
        self.pos = pos               # 1-based genomic start of the aligned window
        self.strand = strand
        self.bulge_type = bulge_type  # "X0" | "DNA" | "RNA"
        self.bulge_size = int(bulge_size or 0)
        self.mismatches = mismatches
        self.pam_len = int(pam_len)
        self.extra = extra or {}

    def genomic_span(self):
        """Number of reference bases the aligned off-target consumes (gaps excluded)."""
        return len(self.dna) - self.dna.count("-")


def cfd_for_site(site, mm_scores, pam_scores):
    """CFD for one site from its (in-place) aligned crRNA/DNA.

    Inputs follow the standard display / Cas-OFFinder convention: gaps are placed
    IN-PLACE at the bulge position (DNA bulge -> '-' in crRNA, RNA bulge -> '-' in
    DNA), crRNA and DNA are the same length.

    * mismatch-only / RNA bulge: calc_cfd(crRNA, DNA[:-pam], DNA[-2:]). The RNA-bulge
      gap sits in DNA, so the 20 guide bases keep CFD positions 1..20 and the bulged
      guide base (vs DNA '-') contributes factor 1 -- verified to reproduce CRISPRme's
      reported CFD exactly on mismatch + RNA-bulge rows.
    * DNA bulge: the extra DNA base has no guide partner, so DELETE the gap column(s)
      (where crRNA has '-') from BOTH strands, then score the 20-vs-20. CFD is
      position-WEIGHTED, so keeping the in-place gap would shift positions 1..20 by
      the bulge width and mis-weight every downstream mismatch; dropping the looped-out
      base restores the correct positions. (This is the textbook mismatch-only CFD with
      the bulged base excluded; it does not try to reproduce CRISPRme's internal
      DNA-bulge gap placement, which is not recoverable from an in-place alignment and
      can differ -- a DNA-bulge caveat the caller should surface.)
    """
    pam = site.pam_len
    dna = site.dna
    guide = site.crrna
    if site.bulge_type == "DNA" and "-" in guide:
        keep = [i for i, ch in enumerate(guide) if ch != "-"]
        guide = "".join(guide[i] for i in keep)
        dna = "".join(dna[i] for i in keep)
    return calc_cfd(guide, dna[:-pam], dna[-2:], mm_scores, pam_scores, True)


# CRISPR-Bulge null rules -- identical to preprocess_CRISPR_BULGE_score._pair in
# analisi_indels_NNN.py: the seq_len=24 model scores a protospacer+PAM alignment
# with AT MOST ONE 1-bp bulge. Score only what the model represents; NULL (-1.0)
# everything else (>1 gap column, N in the off-target, length mismatch, >24).
_CB_MAXLEN = 24
_CB_SG_DUMMY = "A" * 20 + "NGG"
_CB_OFF_DUMMY = "A" * 20 + "AGG"


def _cb_pair(site, index, null_index):
    # mask the guide PAM to NGG exactly as the pipeline does (model scores spacer;
    # PAM normalized), keep the off-target's real PAM + any bulge gaps.
    sg = site.crrna[: len(site.crrna) - site.pam_len] + "NGG"
    off = site.dna
    total_gaps = sg.count("-") + off.count("-")
    if (("N" in off) or (len(off) != len(sg)) or total_gaps > 1
            or len(sg) > _CB_MAXLEN):
        null_index.append(index)
        return _CB_SG_DUMMY, _CB_OFF_DUMMY
    return sg, off


def crispr_bulge_for_sites(sites, device=None):
    """Batch CRISPR-Bulge scores for the sites. Returns list[float] (len==sites),
    -1.0 for rows the 1-bulge model cannot represent or when the cbulge env is
    unavailable (graceful; never raises)."""
    try:
        import scorer_runner  # lazy: only when scoring is requested
    except Exception as e:  # pragma: no cover - env issue, degrade gracefully
        sys.stderr.write(
            "WARNING: CRISPR-Bulge scorer unavailable (%s); emitting -1.0.\n" % e
        )
        return [-1.0] * len(sites)
    null_index = []
    pairs = [_cb_pair(s, i, null_index) for i, s in enumerate(sites)]
    scores = scorer_runner.CRISPR_BULGE_predict_list(
        [p[0] for p in pairs], [p[1] for p in pairs], device=device
    )
    null = set(null_index)
    return [(-1.0 if i in null else scores[i]) for i in range(len(sites))]


# ---------------------------------------------------------------------------
# Carrier overlay -- "in how many people is the site polymorphic"
# ---------------------------------------------------------------------------


def _find_registry_dir(index_path):
    """Locate the registry_<vcf> dir under an index install's Dictionaries/."""
    cands = []
    for base in (index_path, os.path.join(index_path, "Dictionaries")):
        cands += glob.glob(os.path.join(base, "registry_*"))
    cands = [c for c in cands if os.path.isdir(c)]
    if not cands:
        raise FileNotFoundError(
            "no Dictionaries/registry_<vcf> found under %r -- pass an installed "
            "variant index via --index-path (the carrier overlay needs the "
            "registry tier)." % index_path
        )
    return sorted(cands)[0]


class _RegistryOverlay(object):
    """Opens per-chromosome RegistryReaders lazily and answers window queries."""

    def __init__(self, registry_dir):
        import tier0_registry  # lazy: only when overlay is requested
        self._t0 = tier0_registry
        self._dir = registry_dir
        self._readers = {}  # chrom -> RegistryReader | None

    def _reader(self, chrom):
        if chrom not in self._readers:
            binp = os.path.join(self._dir, "reg_%s.bin" % chrom)
            idxp = os.path.join(self._dir, "reg_%s.idx" % chrom)
            if os.path.exists(binp) and os.path.exists(idxp):
                self._readers[chrom] = self._t0.RegistryReader(binp, idxp)
            else:
                self._readers[chrom] = None
        return self._readers[chrom]

    def overlay(self, site):
        """Return a dict summarising variants overlapping the site's window."""
        out = {
            "n_overlapping_variants": 0,
            "max_global_AF": "",
            "max_carrier_count": "",
            "variants": "",
            "per_dataset_AF": "",
        }
        if site.chrom is None or site.pos is None:
            return out
        reader = self._reader(site.chrom)
        if reader is None:
            return out
        start = int(site.pos)
        end = start + site.genomic_span() - 1
        variants = []
        max_af, max_carriers = 0.0, 0
        per_db = {}
        for pos in range(start, end + 1):
            for alt in reader.alts_at(pos):
                groups = reader.lookup(pos, alt)
                g = groups.get(self._t0.GLOBAL_GROUP_ID)
                af = g.allele_freq() if g else 0.0
                carriers = g.n_carrier_indiv if g else 0
                ref = reader.ref(pos, alt) or ""
                rsid = reader.rsid(pos, alt) or "."
                variants.append(
                    "%s_%d_%s_%s:%s:%.3g:%d" % (site.chrom, pos, ref, alt, rsid, af, carriers)
                )
                if af > max_af:
                    max_af = af
                if carriers > max_carriers:
                    max_carriers = carriers
                for gid, cnt in groups.items():
                    if gid == self._t0.GLOBAL_GROUP_ID or self._t0.SEP in gid:
                        continue  # keep per-dataset (not per-subpop) groups only
                    per_db[gid] = max(per_db.get(gid, 0.0), cnt.allele_freq())
        out["n_overlapping_variants"] = len(variants)
        if variants:
            out["max_global_AF"] = "%.6g" % max_af
            out["max_carrier_count"] = max_carriers
            out["variants"] = ";".join(variants)
            out["per_dataset_AF"] = ";".join(
                "%s=%.6g" % (k, v) for k, v in sorted(per_db.items())
            )
        return out

    def close(self):
        for r in self._readers.values():
            if r is not None:
                try:
                    r.close()
                except Exception:
                    pass


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------


def _norm_chrom(c):
    c = str(c).strip()
    return c if c.startswith("chr") else ("chr" + c if c else c)


def parse_casoffinder(path, pam_len=3, position_base=0):
    """Parse a Cas-OFFinder output file (tab-separated).

    Supports both the modern header schema (``#Id/Bulge Type/crRNA/DNA/Chromosome/
    Location/Direction/Mismatches/Bulge Size``) and the legacy 6-column schema
    (``crRNA/Chromosome/Location/DNA/Direction/Mismatches``), auto-detected from the
    header / column count. Cas-OFFinder encodes bulges as '-' in crRNA (DNA bulge)
    or DNA (RNA bulge); mismatch-only rows have no gaps. Positions are 0-based by
    default (Cas-OFFinder convention); set position_base=1 for 1-based input.
    """
    sites = []
    with open(path) as fh:
        header = None
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if header is None and (raw.startswith("#") or _looks_like_header(fields)):
                header = [h.strip().lstrip("#").strip().lower() for h in fields]
                continue
            if header is not None:
                rec = dict(zip(header, fields))
                crrna = rec.get("crrna") or rec.get("crrna (with bulges)")
                dna = rec.get("dna") or rec.get("dna (with bulges)")
                chrom = rec.get("chromosome") or rec.get("chr")
                loc = rec.get("location") or rec.get("position")
                strand = rec.get("direction") or rec.get("strand")
                mm = rec.get("mismatches")
                btype = (rec.get("bulge type") or rec.get("bulge_type") or "").strip()
                bsize = rec.get("bulge size") or rec.get("bulge_size") or 0
            else:
                # legacy 6-col, no header
                if len(fields) < 6:
                    continue
                crrna, chrom, loc, dna, strand, mm = fields[:6]
                btype, bsize = "", 0
            sites.append(
                _make_site(crrna, dna, chrom, loc, strand, mm, btype, bsize,
                           pam_len, position_base)
            )
    return [s for s in sites if s is not None]


def _looks_like_header(fields):
    joined = "\t".join(fields).lower()
    return ("crrna" in joined) and ("dna" in joined)


def _make_site(crrna, dna, chrom, loc, strand, mm, btype, bsize, pam_len, position_base):
    if not crrna or not dna:
        return None
    crrna = crrna.strip().upper()
    dna = dna.strip().upper()
    try:
        pos1 = int(str(loc).strip()) + (1 if position_base == 0 else 0)
    except (TypeError, ValueError):
        pos1 = None
    # derive bulge type/size from the alignment when not given explicitly
    bt = (btype or "").strip().upper()
    if bt in ("", "X", "X0", "0"):
        if crrna.count("-") and not dna.count("-"):
            bt, bsize = "DNA", crrna.count("-")
        elif dna.count("-") and not crrna.count("-"):
            bt, bsize = "RNA", dna.count("-")
        else:
            bt, bsize = "X0", 0
    try:
        bsize = int(str(bsize).strip() or 0)
    except ValueError:
        bsize = crrna.count("-") + dna.count("-")
    try:
        mm = int(str(mm).strip()) if mm not in (None, "") else None
    except ValueError:
        mm = None
    return Site(
        crrna=crrna, dna=dna, chrom=_norm_chrom(chrom), pos=pos1,
        strand=(strand.strip() if strand else None),
        bulge_type=bt, bulge_size=bsize, mismatches=mm, pam_len=pam_len,
        extra={},
    )


def parse_tsv(path, guide_col, dna_col, chrom_col=None, pos_col=None,
              strand_col=None, mm_col=None, btype_col=None, bsize_col=None,
              pam_len=3, position_base=1):
    """Parse a generic named-column TSV into Sites (explicit column mapping)."""
    sites = []
    with open(path) as fh:
        header = None
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if header is None:
                header = [h.strip().lstrip("#").strip() for h in fields]
                continue
            rec = dict(zip(header, fields))
            sites.append(
                _make_site(
                    rec.get(guide_col), rec.get(dna_col), rec.get(chrom_col),
                    rec.get(pos_col), rec.get(strand_col), rec.get(mm_col),
                    rec.get(btype_col) if btype_col else "",
                    rec.get(bsize_col) if bsize_col else 0,
                    pam_len, position_base,
                )
            )
    return [s for s in sites if s is not None]


# ---------------------------------------------------------------------------
# Orchestration + output
# ---------------------------------------------------------------------------

_BASE_COLUMNS = [
    "crRNA_aligned", "DNA_aligned", "Chromosome", "Position", "Strand",
    "Bulge_type", "Bulge_size", "Mismatches", "PAM", "CFD", "CRISPR_BULGE",
]
_OVERLAY_COLUMNS = [
    "n_overlapping_variants", "max_global_AF", "max_carrier_count",
    "variants", "per_dataset_AF",
]


def score_sites(sites, index_path=None, device=None):
    """Score every site with CFD + CRISPR-Bulge and (if index_path) the carrier
    overlay. Returns (header, rows) where rows are lists of strings sorted by CFD
    descending."""
    mm_scores, pam_scores = get_mm_pam_scores()
    cfds = [cfd_for_site(s, mm_scores, pam_scores) for s in sites]
    cbs = crispr_bulge_for_sites(sites, device=device)

    overlay = None
    header = list(_BASE_COLUMNS)
    if index_path:
        overlay = _RegistryOverlay(_find_registry_dir(index_path))
        header += _OVERLAY_COLUMNS

    rows = []
    try:
        for i, s in enumerate(sites):
            pam_seq = s.dna[-s.pam_len:] if len(s.dna) >= s.pam_len else s.dna
            row = [
                s.crrna, s.dna, s.chrom or "", "" if s.pos is None else str(s.pos),
                s.strand or "", s.bulge_type, str(s.bulge_size),
                "" if s.mismatches is None else str(s.mismatches),
                pam_seq, "%.3f" % cfds[i], "%.3f" % cbs[i],
            ]
            if overlay is not None:
                ov = overlay.overlay(s)
                row += [
                    str(ov["n_overlapping_variants"]), str(ov["max_global_AF"]),
                    str(ov["max_carrier_count"]), ov["variants"], ov["per_dataset_AF"],
                ]
            rows.append(row)
    finally:
        if overlay is not None:
            overlay.close()

    cfd_idx = header.index("CFD")
    rows.sort(key=lambda r: float(r[cfd_idx]), reverse=True)
    return header, rows


def write_tsv(header, rows, out_path):
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def write_xlsx(header, rows, out_path):
    try:
        import openpyxl
    except Exception as e:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "openpyxl is required for --xlsx output (%s); use the default TSV." % e
        )
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "targets_scoring"
    ws.append(header)
    for r in rows:
        ws.append(r)
    wb.save(out_path)


def run(targets, output, casoffinder=False, index_path=None, device=None,
        xlsx=False, pam_len=3, position_base=None, columns=None):
    """Top-level entry used by the CLI. ``columns`` is a dict for generic TSV
    parsing (keys guide/dna/chrom/pos/strand/mm/btype/bsize)."""
    if casoffinder:
        pb = 0 if position_base is None else position_base
        sites = parse_casoffinder(targets, pam_len=pam_len, position_base=pb)
    elif columns:
        pb = 1 if position_base is None else position_base
        sites = parse_tsv(
            targets, columns["guide"], columns["dna"], columns.get("chrom"),
            columns.get("pos"), columns.get("strand"), columns.get("mm"),
            columns.get("btype"), columns.get("bsize"), pam_len=pam_len,
            position_base=pb,
        )
    else:
        raise ValueError("specify --casoffinder or a generic --columns mapping")

    if not sites:
        raise ValueError("no sites parsed from %r" % targets)

    header, rows = score_sites(sites, index_path=index_path, device=device)
    if xlsx:
        write_xlsx(header, rows, output)
    else:
        write_tsv(header, rows, output)
    return len(rows)


def _build_argparser():
    p = argparse.ArgumentParser(
        prog="targets-scoring",
        description="Score a supplied list of off-target sites (CFD + CRISPR-Bulge) "
                    "with an optional variant carrier overlay.",
    )
    p.add_argument("--targets", required=True, help="input site list (CasOFFinder or TSV)")
    p.add_argument("--output", required=True, help="output path (.tsv, or .xlsx with --xlsx)")
    p.add_argument("--casoffinder", action="store_true",
                   help="parse the input as Cas-OFFinder output (auto-detect schema)")
    p.add_argument("--index-path", default=None,
                   help="installed variant index (adds the carrier overlay)")
    p.add_argument("--device", default=None,
                   help="CRISPR-Bulge scorer device: cpu|gpu|auto|cuda|metal (default cpu)")
    p.add_argument("--xlsx", action="store_true", help="write .xlsx instead of TSV")
    p.add_argument("--pam-length", type=int, default=3, help="PAM length in bp (default 3)")
    p.add_argument("--position-base", type=int, choices=(0, 1), default=None,
                   help="0- or 1-based positions (default: 0 for CasOFFinder, 1 for TSV)")
    # generic TSV column mapping: the key=val,... form (same syntax as
    # `crisprme.py targets-scoring --columns`) OR the individual --*-column flags.
    p.add_argument("--columns", default=None,
                   help="generic-TSV column map, e.g. "
                        "guide=crRNA,dna=DNA,chrom=Chr,pos=Pos,strand=Strand,mm=Mismatches")
    p.add_argument("--guide-column")
    p.add_argument("--dna-column")
    p.add_argument("--chrom-column")
    p.add_argument("--pos-column")
    p.add_argument("--strand-column")
    p.add_argument("--mismatches-column")
    p.add_argument("--bulge-type-column")
    p.add_argument("--bulge-size-column")
    return p


def main(argv=None):
    parser = _build_argparser()
    args = parser.parse_args(argv)
    if args.device is not None and args.device.lower() not in (
        "cpu", "gpu", "auto", "cuda", "metal"
    ):
        parser.error("--device must be one of: cpu | gpu | auto | cuda | metal")
    columns = None
    if not args.casoffinder:
        # Prefer the --columns key=val,... form (same syntax the crisprme.py wrapper
        # documents); fall back to the individual --*-column flags.
        if args.columns:
            columns = {}
            for kv in args.columns.split(","):
                if "=" not in kv:
                    parser.error("bad --columns entry %r; use key=value pairs" % kv)
                k, v = kv.split("=", 1)
                columns[k.strip()] = v.strip()
            if not (columns.get("guide") and columns.get("dna")):
                parser.error("--columns must map at least guide=<col> and dna=<col>")
        elif args.guide_column and args.dna_column:
            columns = {
                "guide": args.guide_column, "dna": args.dna_column,
                "chrom": args.chrom_column, "pos": args.pos_column,
                "strand": args.strand_column, "mm": args.mismatches_column,
                "btype": args.bulge_type_column, "bsize": args.bulge_size_column,
            }
        else:
            parser.error(
                "for a generic TSV, pass --columns guide=..,dna=.. (or the "
                "--guide-column/--dna-column flags), or --casoffinder for "
                "Cas-OFFinder output."
            )
    n = run(
        args.targets, args.output, casoffinder=args.casoffinder,
        index_path=args.index_path, device=args.device, xlsx=args.xlsx,
        pam_len=args.pam_length, position_base=args.position_base, columns=columns,
    )
    sys.stderr.write("targets-scoring: wrote %d scored sites to %s\n" % (n, args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())

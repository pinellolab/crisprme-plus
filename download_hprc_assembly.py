#!/usr/bin/env python3
"""Download an HPRC release-2 individual's assembly-search input files
directly from the public human-pangenomics S3 bucket, and (optionally)
register them with the Settings page's Data Manager -- the CLI-side
mirror of the existing "download a VCF dataset" precedent
(settings_page.py's add_vcf), for the meeting's requested "HPRC download
script/example."

WHERE THE URLS COME FROM (verified directly against the real S3 bucket
and the real index CSV, 2026-09-02 -- not guessed):
  - The authoritative sample -> file mapping is
    https://github.com/human-pangenomics/hprc_intermediate_assembly's
    data_tables/assemblies_release2_v1.0.index.csv (columns include
    sample_id, haplotype [1=paternal, 2=maternal], assembly_name, source,
    and the FASTA's own s3:// URI). This script always reads the FASTA
    URL from that CSV directly -- it does NOT reconstruct it from a
    naming template, because not every row follows one: e.g. HG002's
    rows have source="extramural" and live under a completely different
    prefix (HPRC_PLUS/.../Q100/pansn/) with different assembly-tool
    provenance (verkko, not hifiasm). Only source == "hprc" rows (432 of
    466 in the current table) are supported by this script; others are
    reported, not guessed at.
  - For a supported ("hprc"-source) row, the chain (vs GRCh38) and
    chromAlias files live as siblings of the FASTA's own
    ".../assemblies/release2/" directory, under
    "annotation/chains/minigraph-cactus/<assembly_name>_vs_GRCh38.chain.gz"
    and "annotation/chrom_assignment/<assembly_name>.chromAlias.txt" --
    confirmed present (real HTTP 200) for two independently-checked
    samples (HG01255, HG00408) before trusting this pattern. Every URL is
    HEAD-checked before use regardless; a missing file is reported
    clearly, never silently skipped.

FASTA SPLITTING: HPRC's release FASTA is one gzipped multi-record file
with PanSN-style headers (e.g. ">HG01255#1#CM086702.1"), but CRISPRme's
search engine expects one plain .fa file per contig, named by its
UCSC-style chromosome name. The downloaded chromAlias.txt's own "assembly"
(PanSN header) -> "ucsc" columns are the exact, already-correct mapping for
this -- confirmed directly against the real file, including compound names
like "chr14_JBHDTB010000001.1_random" for non-canonical contigs. No
heuristic/guessing involved: this script looks up each record's real PanSN
header in that mapping and uses the ucsc name verbatim as the output name.

ON-DISK LAYOUT: everything for one individual lands in a single
self-contained folder with a metadata.json describing every file (see
personal_assembly.py for the full rationale)::

    <data-dir>/Assemblies/HG01255/
      metadata.json
      paternal/ genome/ chr1.fa … chrN.fa
                HG01255_paternal_<acc>_vs_GRCh38.chain.gz
                HG01255_paternal_<acc>.chromAlias.txt
      maternal/ …

so the individual can be removed with a single `rm -rf`, and shared/imported
as a .zip/.tar.gz of the same structure. Writing metadata.json IS the
registration -- the Settings Data Manager and the search forms discover the
individual from it, no separate marker files.

Usage:
    python download_hprc_assembly.py HG01255
    python download_hprc_assembly.py HG01255 --haplotypes paternal
    python download_hprc_assembly.py HG01255 --data-dir /path/to/crisprme-data
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import sys
import urllib.request
from typing import Dict, List, Optional, Tuple

# personal_assembly lives in PostProcess/ (shared with crisprme.py + the web app)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "PostProcess"))
import personal_assembly  # noqa: E402

INDEX_CSV_URL = (
    "https://raw.githubusercontent.com/human-pangenomics/"
    "hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv"
)
S3_HTTPS_BASE = "https://s3-us-west-2.amazonaws.com/human-pangenomics/"
HAPLOTYPE_CODE = {"1": "paternal", "2": "maternal"}


def _s3_to_https(uri: str) -> str:
    """s3://human-pangenomics/<key> -> https://s3-us-west-2.amazonaws.com/human-pangenomics/<key>"""
    prefix = "s3://human-pangenomics/"
    if not uri.startswith(prefix):
        raise ValueError(f"Unexpected S3 URI shape: {uri!r}")
    return S3_HTTPS_BASE + uri[len(prefix):]


_INDEX_CSV_CACHE: Optional[List[Dict[str, str]]] = None


def _fetch_index_csv_rows() -> List[Dict[str, str]]:
    """All rows of the release-2 index CSV, fetched once and cached
    in-process (module-level -- lives for the running process's lifetime,
    e.g. a dev-server worker; dies on restart, which is fine here since
    this is a small convenience cache, not a correctness-critical one).
    """
    global _INDEX_CSV_CACHE
    if _INDEX_CSV_CACHE is None:
        with urllib.request.urlopen(INDEX_CSV_URL, timeout=30) as resp:
            text = resp.read().decode("utf-8")
        _INDEX_CSV_CACHE = list(csv.DictReader(io.StringIO(text)))
    return _INDEX_CSV_CACHE


def fetch_index_rows(sample_id: str) -> List[Dict[str, str]]:
    """Real per-haplotype rows for one sample from the release-2 index CSV."""
    return [row for row in _fetch_index_csv_rows() if row["sample_id"] == sample_id]


def fetch_all_sample_ids() -> List[str]:
    """Sorted, distinct sample_ids from the release-2 index, restricted to
    source == "hprc" rows -- the only rows resolve_haplotype_urls() will
    actually accept, so a picker built from this list never offers a
    choice that's guaranteed to fail at submit time (benchmark/extramural
    entries like HG002 are excluded here, not merely rejected later).
    """
    ids = {row["sample_id"] for row in _fetch_index_csv_rows() if row["source"] == "hprc"}
    return sorted(ids)


def _url_exists(url: str) -> bool:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return resp.status == 200
    except Exception:
        return False


def resolve_haplotype_urls(row: Dict[str, str]) -> Dict[str, str]:
    """{'fasta','fasta_md5','chain','chromalias'} https URLs for one haplotype
    row, built from its own FASTA S3 URI (never guessed from sample_id
    alone -- see module docstring). Raises ValueError with a clear reason
    if this row isn't a supported "hprc"-source, standard-layout row.
    """
    if row["source"] != "hprc":
        raise ValueError(
            f"{row['sample_id']} haplotype {row['haplotype']} has source="
            f"{row['source']!r}, not 'hprc' -- this script only supports the "
            "standard HPRC release-2 layout (verified for 'hprc'-source rows "
            "only). Check its files manually on the HPRC data portal."
        )
    assembly_name = row["assembly_name"]
    fasta_url = _s3_to_https(row["assembly"])
    # .../<sample>/assemblies/release2/<assembly_name>.fa.gz
    #   -> annotation dir is the "release2" dir's own "annotation/" sibling
    release_dir = fasta_url.rsplit("/", 1)[0]
    annotation_dir = release_dir + "/annotation"
    urls = {
        "fasta": fasta_url,
        "fasta_md5": _s3_to_https(row["assembly_md5"]),
        "chain": f"{annotation_dir}/chains/minigraph-cactus/{assembly_name}_vs_GRCh38.chain.gz",
        "chromalias": f"{annotation_dir}/chrom_assignment/{assembly_name}.chromAlias.txt",
    }
    missing = [k for k, u in urls.items() if not _url_exists(u)]
    if missing:
        raise ValueError(
            f"{assembly_name}: expected files not found at the predicted URL(s) "
            f"for {', '.join(missing)} -- the standard layout may not hold for "
            "this sample. Check the HPRC data portal manually rather than "
            "trust a guessed path."
        )
    return urls


def _download(url: str, dest: str, note: str = "") -> None:
    print(f"  downloading {note or url} -> {dest}")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=60) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    os.replace(tmp, dest)


def _verify_md5(path: str, md5_url: str) -> Optional[str]:
    """Returns the expected MD5 hex string if the file matches it, else None."""
    with urllib.request.urlopen(md5_url, timeout=15) as resp:
        expected = resp.read().decode("utf-8").split()[0].strip()
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return expected if h.hexdigest() == expected else None


def split_fasta_by_chromalias(fasta_gz_path: str, chromalias_path: str, dest_dir: str) -> List[str]:
    """Splits a gzipped multi-record HPRC FASTA into one plain .fa file per
    record under dest_dir, named by chromAlias.txt's own "ucsc" column for
    that record's exact PanSN header -- the same mapping (including
    compound "chrN_<accession>_random" names for non-canonical contigs)
    already used by the existing Genomes/<name>/ folders. A header with no
    chromAlias entry is skipped with a warning, never silently renamed by
    guesswork.

    Returns the list of written filenames.
    """
    header_to_ucsc: Dict[str, str] = {}
    with open(chromalias_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                header_to_ucsc[parts[0]] = parts[1]

    os.makedirs(dest_dir, exist_ok=True)
    written = []
    out_fh = None
    skipped = []
    with gzip.open(fasta_gz_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if out_fh:
                    out_fh.close()
                header = line[1:].split()[0].strip()
                ucsc_name = header_to_ucsc.get(header)
                if ucsc_name is None:
                    skipped.append(header)
                    out_fh = None
                    continue
                out_path = os.path.join(dest_dir, f"{ucsc_name}.fa")
                out_fh = open(out_path, "w")
                out_fh.write(f">{ucsc_name}\n")
                written.append(out_path)
            elif out_fh:
                out_fh.write(line)
        if out_fh:
            out_fh.close()
    if skipped:
        print(
            f"  WARNING: {len(skipped)} record(s) had no chromAlias entry, "
            f"skipped (not written): {', '.join(skipped[:5])}"
            + ("…" if len(skipped) > 5 else "")
        )
    return written


def process_individual(
    sample_id: str,
    haplotypes: List[str],
    data_dir: str,
    work_dir: str,
) -> None:
    rows = fetch_index_rows(sample_id)
    if not rows:
        print(f"No rows found for {sample_id!r} in the release-2 index. Typo, or not in this freeze?")
        sys.exit(1)

    by_hap = {}
    for row in rows:
        hap = HAPLOTYPE_CODE.get(row["haplotype"])
        if hap in haplotypes:
            by_hap[hap] = row
    missing_hap = [h for h in haplotypes if h not in by_hap]
    if missing_hap:
        print(f"No {', '.join(missing_hap)} row found for {sample_id} in the index.")
        sys.exit(1)

    # One self-contained folder per individual: <data-dir>/Assemblies/<sample>/
    ind_dir = personal_assembly.individual_dir(data_dir, sample_id)
    source = "HPRC release2 (human-pangenomics)"

    registered = {}
    for hap, row in by_hap.items():
        assembly_name = row["assembly_name"]
        print(f"\n{sample_id} {hap} ({assembly_name}):")
        urls = resolve_haplotype_urls(row)  # raises clearly if this row can't be trusted

        hap_dir = os.path.join(ind_dir, hap)
        chromalias_name = f"{assembly_name}.chromAlias.txt"
        chromalias_dest = os.path.join(hap_dir, chromalias_name)
        _download(urls["chromalias"], chromalias_dest, "chromAlias")

        chain_name = f"{assembly_name}_vs_GRCh38.chain.gz"
        chain_dest = os.path.join(hap_dir, chain_name)
        _download(urls["chain"], chain_dest, "liftOver chain")

        fasta_tmp = os.path.join(work_dir, f"{assembly_name}.fa.gz")
        _download(urls["fasta"], fasta_tmp, "genome FASTA (large)")
        print("  verifying FASTA checksum...")
        fasta_md5 = _verify_md5(fasta_tmp, urls["fasta_md5"])
        if fasta_md5 is None:
            print(f"  MD5 MISMATCH for {fasta_tmp} -- aborting, not registering a possibly-corrupt download.")
            try:
                os.remove(fasta_tmp)  # don't leave a multi-GB corrupt scratch file behind
            except OSError:
                pass
            sys.exit(1)

        genome_dest_dir = os.path.join(hap_dir, "genome")
        print(f"  splitting FASTA into {genome_dest_dir}/ (one file per contig)...")
        written = split_fasta_by_chromalias(fasta_tmp, chromalias_dest, genome_dest_dir)
        print(f"  wrote {len(written)} contig files")
        # Remove the multi-GB scratch FASTA now that per-contig files are written
        # (otherwise ~2GB/haplotype accumulates under the scratch dir).
        try:
            os.remove(fasta_tmp)
        except OSError:
            pass

        # Record this haplotype in the individual's metadata.json (paths relative
        # to the individual folder). Writing metadata IS the registration -- the
        # Data Manager and search forms discover the individual from it.
        personal_assembly.write_haplotype(
            data_dir,
            sample_id,
            hap,
            {
                "assembly_name": assembly_name,
                "genome_dir": f"{hap}/genome",
                "chain": f"{hap}/{chain_name}",
                "chromalias": f"{hap}/{chromalias_name}",
                "n_contigs": len(written),
                "fasta_md5": fasta_md5,
            },
            source=source,
        )
        registered[hap] = assembly_name

    print(
        f"\nRegistered {sample_id} ({', '.join(registered)}) -> "
        f"{os.path.relpath(ind_dir, data_dir)}/ (metadata.json written). "
        "It will appear in the Data Manager and the assembly-search form."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_id", help="HPRC sample id, e.g. HG01255")
    ap.add_argument(
        "--haplotypes", nargs="+", choices=["paternal", "maternal"],
        default=["paternal", "maternal"], help="Which haplotype(s) to fetch (default: both)",
    )
    ap.add_argument(
        "--data-dir", default=".",
        help="CRISPRme data folder; the assembly is written to <data-dir>/Assemblies/<sample>/ (default: .)",
    )
    ap.add_argument(
        "--register", action="store_true",
        help="Accepted for backward compatibility; registration is now automatic (metadata.json is written)",
    )
    ap.add_argument(
        "--work-dir", default=None,
        help="Scratch dir for the raw .fa.gz download (default: <data-dir>/Assemblies/.tmp)",
    )
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    work_dir = args.work_dir or os.path.join(data_dir, personal_assembly.ASSEMBLIES_DIR, ".tmp")
    process_individual(args.sample_id, args.haplotypes, data_dir, work_dir)


if __name__ == "__main__":
    main()

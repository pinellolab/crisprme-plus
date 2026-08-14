#!/usr/bin/env python

from multiprocessing import Pool

import gzip
import subprocess
import sys
import os


def _chrom_from_vcf(vcf_path: str) -> str:
    """Return the chromosome name from the first data line of a VCF (plain or gzipped)."""
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t")[0]
    raise ValueError(f"No data lines found in VCF: {vcf_path}")


def _normalize_chrom(chrom: str) -> str:
    """Ensure chromosome name has chr prefix (e.g. '22' -> 'chr22').
    Some VCF datasets (e.g. 1000G GRCh38) store chromosomes without the
    prefix in the CHROM field, while the genome indices use 'chr'-prefixed names.
    """
    return chrom if chrom.startswith("chr") else "chr" + chrom


def _dataset_chroms(vcf_dir: str, log_indels_dir: str):
    """Chromosome list for a variant dataset. Prefer the raw VCFs; fall back to the
    per-chromosome indel logs bundled with a precomputed index
    (``Dictionaries/log_indels_<vcf>/log<chrom>.txt``) when the multi-GB source VCFs
    were not downloaded. See the twin helper in ``pool_search_indels.py``."""
    if vcf_dir and os.path.isdir(vcf_dir):
        vcfs = sorted(f for f in os.listdir(vcf_dir) if f.endswith("vcf.gz"))
        if vcfs:
            return [_normalize_chrom(_chrom_from_vcf(os.path.join(vcf_dir, f))) for f in vcfs]
    if log_indels_dir and os.path.isdir(log_indels_dir):
        chroms = [
            _normalize_chrom(f[len("log"):-len(".txt")])
            for f in sorted(os.listdir(log_indels_dir))
            if f.startswith("log") and f.endswith(".txt")
        ]
        if chroms:
            return chroms
    raise ValueError(
        f"No VCFs found in '{vcf_dir}' and no indel logs under '{log_indels_dir}'; "
        "cannot determine chromosomes for the indel post-analysis."
    )


# post-analysis script name
POSTANALYSIS = "./post_analisi_indel.sh"

# read input arguments
output_folder = sys.argv[1]
ref_folder = sys.argv[2]
vcf_folder = sys.argv[3]
guide_file = sys.argv[4]
mm = sys.argv[5]
bDNA = sys.argv[6]
bRNA = sys.argv[7]
annotation_file = sys.argv[8]
pam_file = sys.argv[9]
dict_folder = sys.argv[10]
final_res = sys.argv[11]
final_res_alt = sys.argv[12]
ncpus = int(sys.argv[13])


def start_analysis(chrom: str) -> None:
    code = subprocess.call(
        f"{POSTANALYSIS} {output_folder} {ref_folder} {vcf_folder} {guide_file} "
        f"{mm} {bDNA} {bRNA} {annotation_file} {pam_file} {dict_folder} "
        f"{final_res} {final_res_alt} {chrom}",
        shell=True
    )
    if code != 0:
        raise subprocess.SubprocessError(
            f"Post-analysis on indels failed on chromosome {chrom}"
        )


def _detect_budget_gb():
    """Auto-detect the memory budget (GB): CRISPRME_MAX_MEM_GB override, else the
    SMALLEST of host available/total RAM (/proc/meminfo) and any cgroup limit (so a
    memory-limited container is respected); 64 GB fallback. See the SNP twin
    `pool_post_analisi_snp.py` for the rationale (fixes OOM on < 64 GB machines)."""
    env = os.environ.get("CRISPRME_MAX_MEM_GB")
    if env:
        try:
            v = float(env)
            if v > 0:
                return v
        except ValueError:
            pass
    cands = []
    try:
        info = {}
        with open("/proc/meminfo") as fh:
            for line in fh:
                k, _, rest = line.partition(":")
                info[k.strip()] = rest.strip()
        for key in ("MemAvailable", "MemTotal"):
            if key in info:
                cands.append(int(info[key].split()[0]) / (1024.0 ** 2))
                break
    except Exception:
        pass
    for p in ("/sys/fs/cgroup/memory.max",
              "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(p) as fh:
                raw = fh.read().strip()
            if raw and raw != "max":
                b = int(raw)
                if 0 < b < (1 << 62):
                    cands.append(b / (1024.0 ** 3))
        except Exception:
            pass
    return min(cands) if cands else 64.0


def _estimate_worker_gb(dict_folder):
    """Estimate per-worker peak RAM (GB) from the largest .json in the dictionary
    folder (json.load ~2.2x on-disk size + ~1 GB working set); see the SNP twin.
    Overridable via CRISPRME_POSTPROC_WORKER_GB; floors at 4 GB, 6 GB fallback."""
    env = os.environ.get("CRISPRME_POSTPROC_WORKER_GB")
    if env:
        try:
            v = float(env)
            if v > 0:
                return v
        except ValueError:
            pass
    biggest_bytes = 0
    try:
        for f in os.listdir(dict_folder):
            if f.endswith(".json"):
                b = os.path.getsize(os.path.join(dict_folder, f))
                if b > biggest_bytes:
                    biggest_bytes = b
    except Exception:
        pass
    if biggest_bytes <= 0:
        return 6.0
    return max(4.0, (biggest_bytes / (1024.0 ** 3)) * 2.2 + 1.0)


def memory_capped_workers(requested, n_tasks, dict_folder):
    """Bound concurrent post-analysis workers to the machine's ACTUAL memory (see
    the SNP twin `pool_post_analisi_snp.py`). Budget auto-detected minus headroom;
    per-worker estimated from dictionary sizes. Returns
    (workers, budget_gb, per_worker_gb)."""
    budget_gb = _detect_budget_gb()
    usable_gb = max(1.0, budget_gb - max(2.0, budget_gb * 0.15))
    per_worker_gb = _estimate_worker_gb(dict_folder)
    cap = max(1, int(usable_gb // per_worker_gb))
    workers = max(1, min(requested, cap, n_tasks))
    return workers, budget_gb, per_worker_gb


# chromosome list (raw VCFs when present; else the bundled indel logs)
vcf_name = os.path.basename(vcf_folder.rstrip("/"))
chrs = _dataset_chroms(vcf_folder, os.path.join(dict_folder, "log_indels_" + vcf_name))
workers, budget_gb, per_worker_gb = memory_capped_workers(ncpus, len(chrs), dict_folder)
# NOTE: write these diagnostics to STDOUT (log_verbose.txt), never STDERR.
# The caller treats a non-empty stderr log (`[ -s $logerror ]`) as a fatal
# post-analysis failure, so informational text on stderr aborts the run.
sys.stdout.write(
    f"Post-analysis INDELs: {workers} concurrent worker(s) "
    f"(cores={ncpus}, detected memory budget {budget_gb:.1f} GB, "
    f"est. per-worker {per_worker_gb:.1f} GB)\n"
)
usable_gb = max(1.0, budget_gb - max(2.0, budget_gb * 0.15))
if per_worker_gb > usable_gb:
    sys.stdout.write(
        f"WARNING: the largest per-chromosome dictionary needs ~{per_worker_gb:.0f} GB "
        f"of RAM to load, but only ~{usable_gb:.0f} GB is available; this step may be "
        f"OOM-killed. Use a machine with more RAM (>= ~{per_worker_gb + 4:.0f} GB for this "
        f"dataset), give the container more memory, or set CRISPRME_MAX_MEM_GB.\n"
    )
sys.stdout.flush()
with Pool(processes=workers) as pool:  # run chrom-wise post-analysis in parallel
    pool.map(start_analysis, chrs)

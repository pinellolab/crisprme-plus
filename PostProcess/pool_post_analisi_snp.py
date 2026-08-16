#!/usr/bin/env python

from multiprocessing import Pool

import subprocess
import os
import sys

output_folder = sys.argv[1]
ref_folder = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
mm = sys.argv[5]
bDNA = sys.argv[6]
bRNA = sys.argv[7]
annotation_file = sys.argv[8]
pam_file = sys.argv[9]
# sampleID=sys.argv[10]
dict_folder = sys.argv[10]
final_res = sys.argv[11]
final_res_alt = sys.argv[12]
ncpus = int(sys.argv[13])


def start_analysis(chrom):
    cmd = f'./post_analisi_snp.sh "{output_folder}" "{ref_folder}" "{vcf_name}" "{guide_file}" "{mm}" "{bDNA}" "{bRNA}" {annotation_file} {pam_file} {dict_folder} {final_res} {final_res_alt} {chrom}'
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise OSError(f"Post-analysis SNP failed on chromsomes {chrom}")


def _detect_budget_gb():
    """Best-effort detection of the memory budget (GB) for post-analysis.

    Order: explicit CRISPRME_MAX_MEM_GB override, else the SMALLEST of the signals
    we can read -- the host available/total RAM (/proc/meminfo) and any cgroup
    memory limit (so a `docker run --memory=8g` container is respected even though
    /proc/meminfo reports the host's RAM). Falls back to 64 GB only if nothing is
    readable. This is the key fix for OOM kills on machines with < 64 GB: the old
    code assumed a fixed 64 GB budget and over-subscribed workers on small hosts.
    """
    env = os.environ.get("CRISPRME_MAX_MEM_GB")
    if env:
        try:
            v = float(env)
            if v > 0:
                return v
        except ValueError:
            pass
    cands = []
    # host memory (prefer available, fall back to total)
    try:
        info = {}
        with open("/proc/meminfo") as fh:
            for line in fh:
                k, _, rest = line.partition(":")
                info[k.strip()] = rest.strip()
        for key in ("MemAvailable", "MemTotal"):
            if key in info:
                cands.append(int(info[key].split()[0]) / (1024.0 ** 2))  # kB -> GB
                break
    except Exception:
        pass
    # cgroup memory limit (container). v2 then v1; "max"/unlimited sentinels ignored.
    for p in ("/sys/fs/cgroup/memory.max",
              "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(p) as fh:
                raw = fh.read().strip()
            if raw and raw != "max":
                b = int(raw)
                if 0 < b < (1 << 62):  # skip the "unlimited" sentinel
                    cands.append(b / (1024.0 ** 3))
        except Exception:
            pass
    return min(cands) if cands else 64.0


def _estimate_worker_gb(dict_folder):
    """Estimate per-worker peak RAM (GB). Each per-chromosome worker json.load()s
    that chromosome's SNP dictionary, which balloons to ~2.2x its on-disk size in
    RAM (measured: a 13 GB 1000G+HGDP chr2 dict -> ~26 GB), plus ~1 GB for the
    chromosome genome string + working set. Estimated from the LARGEST .json in the
    dictionary folder so the cap reflects the worst-case chromosome. Overridable via
    CRISPRME_POSTPROC_WORKER_GB. Floors at 4 GB; 6 GB fallback if nothing measurable.
    """
    env = os.environ.get("CRISPRME_POSTPROC_WORKER_GB")
    if env:
        try:
            v = float(env)
            if v > 0:
                return v
        except ValueError:
            pass
    # When ijson is available, new_simple_analysis STREAMS only the queried dict
    # entries, so per-worker peak RAM is ~the chromosome genome string + a small
    # subset + scoring working set (<~3 GB), independent of the (up to ~26 GB)
    # dictionary. A modest fixed estimate then lets workers run in parallel again.
    try:
        import ijson  # noqa: F401
        return 3.0
    except ImportError:
        pass
    # No ijson -> new_simple_analysis json.load()s the whole per-chromosome dict
    # (~2.2x its on-disk size), so size the estimate from the largest dictionary.
    biggest_bytes = 0
    biggest_is_gz = False
    try:
        for f in os.listdir(dict_folder):
            if f.endswith(".json") or f.endswith(".json.gz"):
                b = os.path.getsize(os.path.join(dict_folder, f))
                if b > biggest_bytes:
                    biggest_bytes = b
                    biggest_is_gz = f.endswith(".gz")
    except Exception:
        pass
    if biggest_bytes <= 0:
        return 6.0
    # A gzipped dict is ~3.5x smaller on disk but json.load expands it to the full
    # uncompressed size in RAM, so scale the on-disk size up before the 2.2x factor.
    gz_factor = 3.5 if biggest_is_gz else 1.0
    return max(4.0, (biggest_bytes / (1024.0 ** 3)) * gz_factor * 2.2 + 1.0)


def memory_capped_workers(requested, n_tasks, dict_folder):
    """Bound concurrent post-analysis workers to the machine's ACTUAL memory.

    Each per-chromosome worker loads that chromosome's genome string + SNP
    dictionary (json.load), so peak RAM scales with the worker count; on a
    genome-wide 1000G+HGDP run this can spike to ~100 GB at one-worker-per-core.
    Budget is auto-detected (host RAM / cgroup limit; see _detect_budget_gb) minus
    headroom; per-worker RAM is estimated from the dictionary sizes (see
    _estimate_worker_gb). Overridable via CRISPRME_MAX_MEM_GB /
    CRISPRME_POSTPROC_WORKER_GB. Returns (workers, budget_gb, per_worker_gb) so the
    caller can warn if even a single worker will not fit.
    """
    budget_gb = _detect_budget_gb()
    usable_gb = max(1.0, budget_gb - max(2.0, budget_gb * 0.15))
    per_worker_gb = _estimate_worker_gb(dict_folder)
    cap = max(1, int(usable_gb // per_worker_gb))
    workers = max(1, min(requested, cap, n_tasks))
    return workers, budget_gb, per_worker_gb


chroms = [
    os.path.splitext(os.path.basename(f))[0]
    for f in os.listdir(ref_folder)
    if f.endswith(".fa") and not f.endswith(".fai")
]

workers, budget_gb, per_worker_gb = memory_capped_workers(ncpus, len(chroms), dict_folder)
# NOTE: write these diagnostics to STDOUT (log_verbose.txt), never STDERR.
# The caller (submit_job_automated_new_multiple_vcfs.sh) treats a non-empty
# stderr log (`[ -s $logerror ]`) as a fatal post-analysis failure, so any
# informational text on stderr here would abort the run with a false error.
sys.stdout.write(
    f"Post-analysis SNPs: {workers} concurrent worker(s) "
    f"(cores={ncpus}, detected memory budget {budget_gb:.1f} GB, "
    f"est. per-worker {per_worker_gb:.1f} GB)\n"
)
# Warn (loudly, but to STDOUT) when a single worker's estimated peak exceeds the
# budget: the per-chromosome SNP dictionary is simply too large for this machine's
# RAM, so the worker will likely be OOM-killed no matter the concurrency. This
# turns the otherwise-cryptic "Killed ... EmptyDataError: No columns to parse"
# cascade into an actionable message.
usable_gb = max(1.0, budget_gb - max(2.0, budget_gb * 0.15))
if per_worker_gb > usable_gb:
    sys.stdout.write(
        f"WARNING: the largest per-chromosome SNP dictionary needs ~{per_worker_gb:.0f} GB "
        f"of RAM to load, but only ~{usable_gb:.0f} GB is available. This genome-wide "
        f"variant post-analysis may be killed (out of memory). Run on a machine with more "
        f"RAM (>= ~{per_worker_gb + 4:.0f} GB for this dataset), give the container more "
        f"memory, or set CRISPRME_MAX_MEM_GB to your real limit.\n"
    )
sys.stdout.flush()
with Pool(processes=workers) as pool:
    pool.map(start_analysis, chroms)

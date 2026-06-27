#!/usr/bin/env bash
# Submit the two Quantum ESPRESSO mold jobs prepared under qe/glcn and qe/glcnac.
# Run this on the HPC system from the STMFit repository root.

set -euo pipefail

DIRS=("qe/glcn" "qe/glcnac")
MAX_TOTAL_TASKS=8
REPORT="hpc/qe_molds/qe_input_preflight.tsv"
JOBS_TSV="hpc/qe_molds/qe_jobs.tsv"
DRY_RUN=0
WATCH=0
WATCH_POLL=60
JULIA_MODULE_VERSION=${JULIA_MODULE_VERSION:-}
QE_COMPILER_MODULE=${QE_COMPILER_MODULE:-intel/2024.0}
QE_MPI_MODULE=${QE_MPI_MODULE:-impi/2021.11}
QE_MODULE=${QE_MODULE:-qe/7.4.1}
QE_QOS=${QE_QOS:-}
QE_MIN_MEM_MB=${QE_MIN_MEM_MB:-0}
QE_DEPENDENCY=${QE_DEPENDENCY:-}
SEQUENTIAL=0

usage() {
    cat <<'EOF'
Usage: hpc/submit_qe_molds.sh [options]

Options:
  --dir PATH              QE run directory. Repeatable. Default: qe/glcn qe/glcnac
  --max-total-tasks INT   Max simultaneous Slurm ntasks (sum if parallel, max-per-job if --sequential) [8]
  --min-mem-mb INT        Optional minimum #SBATCH --mem per job [QE_MIN_MEM_MB or 0]
  --report PATH           Preflight report TSV [hpc/qe_molds/qe_input_preflight.tsv]
  --jobs-tsv PATH         Submitted job-id report [hpc/qe_molds/qe_jobs.tsv]
  --watch                 Poll squeue until submitted jobs finish, then show sacct
  --watch-poll SEC        Poll interval for --watch [60]
  --dependency SPEC       Slurm dependency for queued jobs, e.g. afterany:12345
  --sequential            Submit each dir with afterok dependency on the previous job
  --dry-run               Validate and print sbatch commands without submitting
  -h, --help              Show this help

Run from the STMFit repository root on the HPC system. The script does not sync
files; it assumes QE run directories are already present on the cluster. Full
runs submit run_qe_mold.sbatch; SCF+PP-only preliminary dirs may instead provide
run_scf_pp.sbatch.
EOF
}

run_script_for_dir() {
    local dir="$1"
    if [[ -f "$dir/run_qe_mold.sbatch" ]]; then
        printf '%s\n' "run_qe_mold.sbatch"
    elif [[ -f "$dir/run_scf_pp.sbatch" ]]; then
        printf '%s\n' "run_scf_pp.sbatch"
    else
        echo "ERROR: missing run_qe_mold.sbatch or run_scf_pp.sbatch in $dir" >&2
        exit 1
    fi
}

custom_dirs=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)
            if (( custom_dirs == 0 )); then DIRS=(); custom_dirs=1; fi
            DIRS+=("$2"); shift 2 ;;
        --dir=*)
            if (( custom_dirs == 0 )); then DIRS=(); custom_dirs=1; fi
            DIRS+=("${1#*=}"); shift ;;
        --max-total-tasks) MAX_TOTAL_TASKS="$2"; shift 2 ;;
        --max-total-tasks=*) MAX_TOTAL_TASKS="${1#*=}"; shift ;;
        --min-mem-mb) QE_MIN_MEM_MB="$2"; shift 2 ;;
        --min-mem-mb=*) QE_MIN_MEM_MB="${1#*=}"; shift ;;
        --report) REPORT="$2"; shift 2 ;;
        --report=*) REPORT="${1#*=}"; shift ;;
        --jobs-tsv) JOBS_TSV="$2"; shift 2 ;;
        --jobs-tsv=*) JOBS_TSV="${1#*=}"; shift ;;
        --watch) WATCH=1; shift ;;
        --watch-poll) WATCH_POLL="$2"; shift 2 ;;
        --watch-poll=*) WATCH_POLL="${1#*=}"; shift ;;
        --dependency) QE_DEPENDENCY="$2"; shift 2 ;;
        --dependency=*) QE_DEPENDENCY="${1#*=}"; shift ;;
        --sequential) SEQUENTIAL=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ ! -f Project.toml || ! -d test || ! -d hpc ]]; then
    echo "ERROR: run from the STMFit repository root." >&2
    exit 1
fi

preflight_args=()
for dir in "${DIRS[@]}"; do
    preflight_args+=(--dir "$dir")
done
preflight_args+=(--out "$REPORT" --max-total-tasks "$MAX_TOTAL_TASKS")
[[ "$QE_MIN_MEM_MB" != "0" ]] && preflight_args+=(--min-mem-mb "$QE_MIN_MEM_MB")
(( SEQUENTIAL )) && preflight_args+=(--sequential)

echo "Preflight: ${DIRS[*]}"
if ! command -v julia >/dev/null 2>&1; then
    if type module >/dev/null 2>&1; then
        module purge
        module load "julia${JULIA_MODULE_VERSION:+/$JULIA_MODULE_VERSION}"
    fi
fi
command -v julia >/dev/null 2>&1 || { echo "ERROR: julia not found. Load a Julia module or set JULIA_MODULE_VERSION." >&2; exit 1; }
julia --project=. test/preflight_qe_mold_inputs.jl "${preflight_args[@]}"

if (( DRY_RUN )); then
    echo "Dry run: no jobs submitted. Commands would be:"
    prev_label=""
    for dir in "${DIRS[@]}"; do
        script=$(run_script_for_dir "$dir")
        dep=""
        if (( SEQUENTIAL )) && [[ -n "$prev_label" ]]; then
            dep=" --dependency=afterok:$prev_label"
        elif [[ -n "$QE_DEPENDENCY" ]]; then
            dep=" --dependency=$QE_DEPENDENCY"
        fi
        echo "  (cd '$dir' && sbatch$dep $script)"
        prev_label="<previous-job-id>"
    done
    exit 0
fi

command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found. Run on the HPC login node." >&2; exit 1; }

mkdir -p "$(dirname "$JOBS_TSV")"
printf "run_dir\tjob_id\tsubmit_output\n" > "$JOBS_TSV"
job_ids=()
prev_job_id=""
for dir in "${DIRS[@]}"; do
    echo "Submitting $dir"
    script=$(run_script_for_dir "$dir")
    sbatch_args=(--export=ALL,JULIA_MODULE_VERSION="$JULIA_MODULE_VERSION",QE_COMPILER_MODULE="$QE_COMPILER_MODULE",QE_MPI_MODULE="$QE_MPI_MODULE",QE_MODULE="$QE_MODULE")
    [[ -n "$QE_QOS" ]] && sbatch_args=(--qos="$QE_QOS" "${sbatch_args[@]}")
    if (( SEQUENTIAL )) && [[ -n "$prev_job_id" ]]; then
        sbatch_args=(--dependency="afterok:$prev_job_id" "${sbatch_args[@]}")
    elif [[ -n "$QE_DEPENDENCY" ]]; then
        sbatch_args=(--dependency="$QE_DEPENDENCY" "${sbatch_args[@]}")
    fi
    submit_out=$(cd "$dir" && sbatch "${sbatch_args[@]}" "$script")
    job_id=$(printf "%s\n" "$submit_out" | grep -oE '[0-9]+' | head -1)
    [[ -n "$job_id" ]] || { echo "ERROR: could not parse job id from: $submit_out" >&2; exit 1; }
    job_ids+=("$job_id")
    (( SEQUENTIAL )) && prev_job_id="$job_id"
    printf "%s\t%s\t%s\n" "$dir" "$job_id" "$submit_out" >> "$JOBS_TSV"
    echo "  $dir -> $job_id"
done

echo "Submitted jobs: ${job_ids[*]}"
echo "Job report: $JOBS_TSV"

if (( WATCH )); then
    command -v squeue >/dev/null 2>&1 || { echo "ERROR: squeue not found; cannot watch." >&2; exit 1; }
    ids_csv=$(IFS=,; echo "${job_ids[*]}")
    echo "Watching jobs $ids_csv"
    while :; do
        n=$(squeue -h -j "$ids_csv" 2>/dev/null | wc -l)
        [[ "$n" =~ ^[0-9]+$ ]] || n=0
        (( n == 0 )) && break
        printf "  [%s] %s jobs/tasks still queued/running...\r" "$(date +%H:%M:%S)" "$n"
        sleep "$WATCH_POLL"
    done
    echo
    echo "Jobs finished. Accounting summary:"
    command -v sacct >/dev/null 2>&1 && sacct -j "$ids_csv" --format=JobID,JobName,State,Elapsed,ExitCode,MaxRSS -n || true
fi

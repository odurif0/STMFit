#!/usr/bin/env bash
#SBATCH --job-name=stmfit-local-sigma
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
# Pass all four variables explicitly in sbatch --export; no login-node fits.
[[ $# -eq 0 || ( $# -eq 1 && "$1" = --dry-run ) ]] || exit 2
: "${STMFIT_PROJECT_DIR:?Set source checkout}"
: "${STMFIT_RAW_DIR:?Set complete raw cohort}"
: "${STMFIT_OUTDIR:?Set new output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_RAW_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
args=(--data-dir "$STMFIT_RAW_DIR" --count-config config/chitosan.toml --settings config/local_sigma_counting.toml)
if [[ $# -eq 1 ]]; then
    for chunk in 1 2 3 4; do
        "$JULIA_BIN" --startup-file=no --threads=1 --project=. test/diagnose_local_sigma_counting.jl \
            "${args[@]}" --outdir "$STMFIT_OUTDIR/chunk$chunk" --chunk "$chunk/4" --dry-run
    done
    exit 0
fi
: "${SLURM_JOB_ID:?Run scientific comparison in Slurm}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || exit 2
mkdir "$STMFIT_OUTDIR"
pids=()
for chunk in 1 2 3 4; do
    "$JULIA_BIN" --startup-file=no --threads=1 --project=. test/diagnose_local_sigma_counting.jl \
        "${args[@]}" --outdir "$STMFIT_OUTDIR/chunk$chunk" --chunk "$chunk/4" > "$STMFIT_OUTDIR/chunk$chunk.log" 2>&1 &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do wait "$pid" || status=1; done
[[ "$status" -eq 0 ]] || { echo "Incomplete fits; keep all outputs, no retry or grade" >&2; exit 1; }
exec "$JULIA_BIN" --startup-file=no --threads=1 --project=. test/verify_local_sigma_counting.jl "$STMFIT_OUTDIR" "$STMFIT_RAW_DIR"

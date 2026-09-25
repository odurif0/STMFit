#!/usr/bin/env bash
#SBATCH --job-name=stmfit-background
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
[[ $# -eq 0 || ( $# -eq 1 && "$1" = --dry-run ) ]] || exit 2
: "${STMFIT_PROJECT_DIR:?Set source checkout}"
: "${STMFIT_INPUT_DIR:?Set saved counting run}"
: "${STMFIT_OUTDIR:?Set fresh output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for input_path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$input_path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
args=(--input-root "$STMFIT_INPUT_DIR" --settings config/background_conditioning.toml)
if [[ $# -eq 1 ]]; then
    exec "$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_background_conditioning.jl \
        "${args[@]}" --outdir "$STMFIT_OUTDIR" --dry-run
fi
: "${SLURM_JOB_ID:?Run scientific fits on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || exit 2
mkdir "$STMFIT_OUTDIR"
pids=()
for chunk in 1 2 3 4; do
    "$JULIA_BIN" --startup-file=no -t1 --project=. test/diagnose_background_conditioning.jl \
        "${args[@]}" --outdir "$STMFIT_OUTDIR/chunk$chunk" --chunk "$chunk/4" > "$STMFIT_OUTDIR/chunk$chunk.log" 2>&1 &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do wait "$pid" || status=1; done
[[ "$status" -eq 0 ]] || { echo "Incomplete comparison; keep all outputs, no retry or grade" >&2; exit 1; }
exec "$JULIA_BIN" --startup-file=no -t1 --project=. test/verify_background_conditioning.jl "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR"

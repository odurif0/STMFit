#!/usr/bin/env bash
#SBATCH --job-name=stmfit-complete-observation
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --no-requeue
set -euo pipefail
# Four one-thread workers; each verifies every saved bank before completion.
# Viper requires explicit command-line sbatch --export for these variables.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set committed source checkout}"
: "${STMFIT_INPUT_DIR:?Set first complete observed-fit output directory}"
: "${STMFIT_OUTDIR:?Set new diagnostic output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 GKSwstype=100
julia=("$JULIA_BIN" --startup-file=no -t 1 --project=.)
command=("${julia[@]}" test/diagnose_scan_exclusion.jl --input "$STMFIT_INPUT_DIR"
    --settings config/scan_exclusion_complete_observed.toml)
if [[ $# -eq 1 ]]; then
    for arm in control observed; do
        "${command[@]}" --outdir "$STMFIT_OUTDIR/repeat1/$arm" --arm "$arm" --dry-run
    done
    exit
fi
: "${SLURM_JOB_ID:?Run scientific diagnostic on a compute node}"
[[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
mkdir "$STMFIT_OUTDIR"
pids=()
for repetition in 1 2; do
    mkdir "$STMFIT_OUTDIR/repeat$repetition"
    for arm in control observed; do
        (
            output="$STMFIT_OUTDIR/repeat$repetition/$arm"
            "${command[@]}" --outdir "$output" --arm "$arm"
            "${julia[@]}" test/verify_scan_exclusion.jl "$STMFIT_INPUT_DIR" "$output"
        ) > "$STMFIT_OUTDIR/repeat${repetition}_${arm}.log" 2>&1 &
        pids+=("$!")
    done
done
status=0
for pid in "${pids[@]}"; do
    if wait "$pid"; then :; else status=1; fi
done
[[ "$status" -eq 0 ]] || exit "$status"
for repetition in 1 2; do
    "${julia[@]}" test/verify_scan_exclusion.jl --complete-pair "$STMFIT_INPUT_DIR" \
        "$STMFIT_OUTDIR/repeat$repetition/control" "$STMFIT_OUTDIR/repeat$repetition/observed"
done

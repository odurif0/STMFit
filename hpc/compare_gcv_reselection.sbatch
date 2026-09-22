#!/usr/bin/env bash
#SBATCH --job-name=stmfit-gcv-reselect
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail

# One job; four one-thread fit shards, then whole-cohort assignment.
# On Viper supply all five variables below in sbatch's explicit --export list.
# No grading inputs, geometry caches, or batch policy substitution.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set the source checkout}"
: "${STMFIT_INPUT_DIR:?Set raw/template input directory}"
: "${STMFIT_SELECTED_SUMMARY:?Set original cached-N summary}"
: "${STMFIT_OUTDIR:?Set a new output directory}"
: "${JULIA_BIN:?Set the Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_SELECTED_SUMMARY" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required: $path" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || {
    echo "Output already exists: $STMFIT_OUTDIR" >&2; exit 2;
}
[[ -f "$STMFIT_SELECTED_SUMMARY" && -f "$STMFIT_INPUT_DIR/templates_cc.tsv" \
   && -d "$STMFIT_INPUT_DIR/full146_raw" ]] || { echo "Missing common inputs" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
if [[ $# -eq 0 ]]; then
    : "${SLURM_JOB_ID:?Scientific comparison must run in a Slurm job}"
    [[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
fi
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=4
exec "$JULIA_BIN" --startup-file=no -t 4 --project=. test/run_gcv_reselection.jl \
    --data-dir "$STMFIT_INPUT_DIR/full146_raw" --count-config config/chitosan.toml \
    --config config/unit_assignment_patch_support.toml \
    --selected-summary "$STMFIT_SELECTED_SUMMARY" --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" \
    --outdir "$STMFIT_OUTDIR" "$@"

#!/bin/bash
#SBATCH --job-name=stmfit-labelfree-diagnostic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=02:00:00

# Four independent one-thread cases; no fitting on login nodes.
# Required paths are explicit and no completed run may be overwritten.
# Viper exports SBATCH_EXPORT=NONE in the submitting environment. The caller
# must override it on the sbatch command line (a #SBATCH directive loses to it):
# --export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN
# A shell --dry-run or sbatch --test-only does not test compute-node export.
set -euo pipefail
: "${STMFIT_PROJECT_DIR:?required}"
: "${STMFIT_INPUT_DIR:?required}"
: "${STMFIT_OUTDIR:?required}"
: "${JULIA_BIN:?Julia 1.13 executable required}"
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--dry-run" ) ]]; then
    echo "Only optional --dry-run is accepted" >&2
    exit 2
fi
export GKSwstype=100
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$STMFIT_PROJECT_DIR"
exec "$JULIA_BIN" --startup-file=no --threads=1 --project="$STMFIT_PROJECT_DIR"     test/run_label_free_exploration.jl     --input-dir "$STMFIT_INPUT_DIR"     --config config/chitosan_10_20mer_adaptive_support_rescue.toml     --settings config/label_free_exploration.toml     --outdir "$STMFIT_OUTDIR" "$@"

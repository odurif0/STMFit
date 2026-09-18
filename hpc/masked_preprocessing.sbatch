#!/bin/bash
#SBATCH --job-name=stmfit-masked-background
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=00:30:00

# One bounded four-file diagnostic, no count/classifier fits or grades.
# Viper sets SBATCH_EXPORT=NONE. The caller MUST override it on the command line:
# --export=STMFIT_PROJECT_DIR,STMFIT_INPUT_DIR,STMFIT_OUTDIR,JULIA_BIN
# A #SBATCH directive does not override that submitting environment variable.
# Shell --dry-run and sbatch --test-only do NOT verify compute-node export.
# This script never submits, retries, requeues or overwrites an existing run.
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
exec "$JULIA_BIN" --startup-file=no --threads=1 --project="$STMFIT_PROJECT_DIR" \
    test/run_masked_preprocessing.jl \
    --input-dir "$STMFIT_INPUT_DIR" \
    --config config/chitosan_10_20mer_adaptive_support_rescue.toml \
    --acquisition-settings config/label_free_exploration.toml \
    --settings config/masked_robust_preprocessing.toml \
    --outdir "$STMFIT_OUTDIR" "$@"

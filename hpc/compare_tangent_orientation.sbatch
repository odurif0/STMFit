#!/usr/bin/env bash
#SBATCH --job-name=stmfit-tangent-angle
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
# Explicitly export all five paths on Viper. No retries or login-node science.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set source checkout}"
: "${STMFIT_CACHE_DIR:?Set original symmetric geometry directory}"
: "${STMFIT_INPUT_DIR:?Set full146_raw and templates_cc.tsv directory}"
: "${STMFIT_OUTDIR:?Set new output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_CACHE_DIR" "$STMFIT_INPUT_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
if [[ $# -eq 0 ]]; then
    : "${SLURM_JOB_ID:?Run real comparison on a compute node}"
    [[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
fi
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=4
exec "$JULIA_BIN" --startup-file=no -t 4 --project=. test/run_tangent_mold_comparison.jl \
    --data-dir "$STMFIT_INPUT_DIR/full146_raw" --count-config config/chitosan.toml \
    --config config/unit_assignment_patch_support.toml \
    --reference-settings config/tangent_mold_projection.toml --settings config/tangent_mold_orientation.toml \
    --features "$STMFIT_CACHE_DIR/features.tsv" --split-features "$STMFIT_CACHE_DIR/features_split.tsv" \
    --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" --outdir "$STMFIT_OUTDIR" "$@"

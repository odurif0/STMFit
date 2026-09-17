#!/usr/bin/env bash
#SBATCH --job-name=stmfit-full146
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --output=results/reconstructed_cc_soft_v1/full146-%j.log
set -euo pipefail

# Fixed-cache full146 reconstruction only; no counting/base refit. The native
# runner derives N from base_geometry_full146.tsv and computes genuine split
# fits on all supplied files with at most four single-thread children, then fits
# classifiers on the whole cohort. No split cache or external labels enter here.
# Required: STMFIT_PROJECT_DIR, JULIA_BIN (1.13, checked by the native runner),
# STMFIT_INPUT_DIR (absolute directory containing base_geometry_full146.tsv,
# full146_raw/ and templates_cc.tsv), STMFIT_OUTDIR.
# Optional: STMFIT_SELECTED_SUMMARY explicitly names a compatible original
# summary. Empty/unset means features-only; no summary is auto-discovered.
# STMFIT_OUTDIR is the NEW direct full146 output directory, not a parent root;
# relative output paths resolve from STMFIT_PROJECT_DIR. No /full146 is appended.
# Metadata-only check: bash hpc/reconstructed_full146.sbatch --dry-run
# Submit from the synced project root; the Slurm log parent must already exist.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 2
fi
cd "${STMFIT_PROJECT_DIR:?Set the synced project directory}"
: "${JULIA_BIN:?Set the Julia 1.13 executable}"
: "${STMFIT_INPUT_DIR:?Set the prepared full146 input directory}"
: "${STMFIT_OUTDIR:?Set a new direct full146 output directory}"
export GKSwstype=100
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=4

summary_args=()
if [[ -n "${STMFIT_SELECTED_SUMMARY:-}" ]]; then
    summary_args=(--selected-summary "$STMFIT_SELECTED_SUMMARY")
fi
exec "$JULIA_BIN" -t "$JULIA_NUM_THREADS" --project=. test/run_reconstructed_chitosan.jl \
  --data-dir "$STMFIT_INPUT_DIR/full146_raw" \
  --count-config config/chitosan.toml \
  --config config/unit_assignment_reconstructed.toml \
  --features "$STMFIT_INPUT_DIR/base_geometry_full146.tsv" \
  --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" \
  --outdir "$STMFIT_OUTDIR" "${summary_args[@]}" "$@"

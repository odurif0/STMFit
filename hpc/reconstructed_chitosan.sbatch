#!/usr/bin/env bash
#SBATCH --job-name=stmfit-native
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
#SBATCH --output=results/reconstructed_cc_soft_v1/native-%j.log
set -euo pipefail

cd "${STMFIT_PROJECT_DIR:?Set the synced project directory}"
export GKSwstype=100
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
"${JULIA_BIN:?Set the Julia 1.13 executable}" --version

# Application first. No saved auto/manual count is selected or mixed in.
"$JULIA_BIN" -t "$JULIA_NUM_THREADS" --project=. test/run_reconstructed_chitosan.jl \
  --data-dir "${STMFIT_INPUT_DIR:?Set the prepared raw-input directory}/unknown_raw" \
  --count-config config/chitosan_10_20mer_adaptive_support_rescue.toml \
  --config config/unit_assignment_reconstructed.toml \
  --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" \
  --outdir "${STMFIT_OUTDIR:?Set a new output directory}/unknown25"

# Fixed saved selected-N geometry, all 146 files; labels are not inputs.
# Regenerate genuine split fits and patches. Grading/comparison happens locally.
"$JULIA_BIN" -t "$JULIA_NUM_THREADS" --project=. test/run_reconstructed_chitosan.jl \
  --data-dir "$STMFIT_INPUT_DIR/full146_raw" \
  --count-config config/chitosan.toml \
  --config config/unit_assignment_reconstructed.toml \
  --features "$STMFIT_INPUT_DIR/base_geometry_full146.tsv" \
  --templates "$STMFIT_INPUT_DIR/templates_cc.tsv" \
  --outdir "$STMFIT_OUTDIR/full146"

#!/usr/bin/env bash
#SBATCH --job-name=stmfit-mold-loo
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000MB
set -euo pipefail
# Pass all four paths explicitly with sbatch --export on Viper.
if [[ $# -gt 1 ]] || [[ $# -eq 1 && "$1" != "--dry-run" ]]; then
    echo "Usage: $0 [--dry-run]" >&2; exit 2
fi
: "${STMFIT_PROJECT_DIR:?Set committed source directory}"
: "${STMFIT_REFERENCE_DIR:?Set saved finite tangent reference directory}"
: "${STMFIT_OUTDIR:?Set new output directory}"
: "${JULIA_BIN:?Set Julia 1.13 executable}"
for path in "$STMFIT_PROJECT_DIR" "$STMFIT_REFERENCE_DIR" "$STMFIT_OUTDIR" "$JULIA_BIN"; do
    [[ "$path" = /* ]] || { echo "Absolute paths required" >&2; exit 2; }
done
[[ ! -e "$STMFIT_OUTDIR" && ! -L "$STMFIT_OUTDIR" ]] || { echo "Output exists" >&2; exit 2; }
cd "$STMFIT_PROJECT_DIR"
if [[ $# -eq 0 ]]; then
    : "${SLURM_JOB_ID:?Run real inference on a compute node}"
    [[ "${SLURM_CPUS_PER_TASK:-0}" -ge 4 ]] || { echo "Four CPUs required" >&2; exit 2; }
fi
export GKSwstype=100 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=4
exec "$JULIA_BIN" --startup-file=no -t 4 --project=. test/run_mold_loo_comparison.jl \
    --reference-dir "$STMFIT_REFERENCE_DIR" --config config/unit_assignment_patch_support.toml \
    --settings config/mold_leave_one_out.toml --outdir "$STMFIT_OUTDIR" "$@"


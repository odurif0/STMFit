#!/usr/bin/env bash
# Read-only Fourier/profile analysis, not a QE calculation.
#SBATCH --job-name=qe-vacuum-crossings
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_crossings_%j.out
#SBATCH --error=qe_crossings_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
STMFIT_CROSS_ROOT=$(cd ../.. && pwd)
STMFIT_GLCN_WFC=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save/wfc1.dat
STMFIT_GLCNAC_WFC=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save/wfc1.dat
test ! -e analysis
sha256sum "$STMFIT_GLCN_WFC" "$STMFIT_GLCNAC_WFC" > wavefunctions_before.sha256
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_CROSS_ROOT" \
    "$STMFIT_CROSS_ROOT/test/qe_vacuum_crossings.jl" analyze . "$STMFIT_GLCN_WFC" "$STMFIT_GLCNAC_WFC"
sha256sum -c --quiet wavefunctions_before.sha256

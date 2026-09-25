#!/usr/bin/env bash
# Interval-enclosed Fourier analysis; no QE calculation or new checkpoint.
#SBATCH --job-name=qe-vacuum-continuous
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_continuous_%j.out
#SBATCH --error=qe_continuous_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
STMFIT_CONT_ROOT=$(cd ../.. && pwd)
STMFIT_GLCN_WFC=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save/wfc1.dat
STMFIT_GLCNAC_WFC=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save/wfc1.dat
test ! -e glcn/continuous
test ! -e glcnac/continuous
sha256sum "$STMFIT_GLCN_WFC" "$STMFIT_GLCNAC_WFC" > wavefunctions_before.sha256
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_CONT_ROOT" \
    "$STMFIT_CONT_ROOT/test/qe_vacuum_continuous.jl" analyze . glcn "$STMFIT_GLCN_WFC"
"$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_CONT_ROOT" \
    "$STMFIT_CONT_ROOT/test/qe_vacuum_continuous.jl" analyze . glcnac "$STMFIT_GLCNAC_WFC"
sha256sum -c --quiet wavefunctions_before.sha256
echo 'VACUUM_CONTINUOUS_COMPLETE; all outcomes retained; no grade or production change'

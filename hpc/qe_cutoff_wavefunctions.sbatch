#!/usr/bin/env bash
# One read-only saved-orbital analysis job. No electronic calculation or retry.
#SBATCH --job-name=qe-cutoff-wfc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_cutoff_wfc_%j.out
#SBATCH --error=qe_cutoff_wfc_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_CPUS_PER_TASK:?} == 4 ]] || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
STMFIT_CW_ROOT=$(cd ../.. && pwd)
STMFIT_CW_RUN=$PWD
STMFIT_CW_GLCN=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save/wfc1.dat
STMFIT_CW_GLCNAC=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save/wfc1.dat
sha256sum -c --quiet input_hashes.sha256
test ! -e comparison
sha256sum "$STMFIT_CW_GLCN" "$STMFIT_CW_GLCNAC" \
    "$STMFIT_CW_ROOT/qe/density_cutoff_20260926/glcn/work/glcn_central.save/wfc1.dat" \
    "$STMFIT_CW_ROOT/qe/density_cutoff_20260926/glcnac/work/glcnac_central.save/wfc1.dat" > wavefunctions_before.sha256
trap 'sha256sum -c --quiet wavefunctions_before.sha256' EXIT
for molecule in glcn glcnac; do
    if [[ $molecule == glcn ]]; then STMFIT_CW_OLD=$STMFIT_CW_GLCN; else STMFIT_CW_OLD=$STMFIT_CW_GLCNAC; fi
    for cutoff in 360 720; do
        STMFIT_CW_FILE=$STMFIT_CW_OLD
        if [[ $cutoff == 720 ]]; then STMFIT_CW_FILE=$STMFIT_CW_ROOT/qe/density_cutoff_20260926/$molecule/work/${molecule}_central.save/wfc1.dat; fi
        "$JULIA_BIN" --startup-file=no -t 4 --project="$STMFIT_CW_ROOT" \
            "$STMFIT_CW_ROOT/test/qe_cutoff_wavefunctions.jl" analyze "$STMFIT_CW_RUN" "$molecule" "$cutoff" "$STMFIT_CW_FILE"
    done
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_CW_ROOT" \
    "$STMFIT_CW_ROOT/test/qe_cutoff_wavefunctions.jl" summarize "$STMFIT_CW_RUN"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_CW_ROOT" \
    "$STMFIT_CW_ROOT/test/test_qe_cutoff_wavefunctions.jl" "$STMFIT_CW_RUN" > independent_verification.out 2>&1
sha256sum -c --quiet input_hashes.sha256 wavefunctions_before.sha256
echo 'CUTOFF_WAVEFUNCTION_COMPARISON_COMPLETE; no mold, calibration, fit or benchmark adopted'

#!/usr/bin/env bash
# One state per job, at most two jobs sequentially. No retry or requeue.
#SBATCH --job-name=qe-rho720
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_cutoff_%j.out
#SBATCH --error=qe_cutoff_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_NTASKS:?} == 8 ]] || exit 1
STMFIT_RHO_STATE=$(basename "$PWD")
case "$STMFIT_RHO_STATE" in
    glcn) STMFIT_RHO_SOURCE=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save ;;
    glcnac) STMFIT_RHO_SOURCE=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save ;;
    *) exit 1 ;;
esac
STMFIT_RHO_ROOT=$(cd ../../.. && pwd)
STMFIT_RHO_RUN=$(cd .. && pwd)
STMFIT_RHO_PP=$STMFIT_RHO_ROOT/qe/xc_build_control_20260926/binaries/pp_stock.x
sha256sum -c --quiet input_hashes.sha256
(cd "$STMFIT_RHO_ROOT" && sha256sum -c --quiet "$STMFIT_RHO_RUN/stock_pp.sha256")
test ! -e work
cmp accepted.xml "$STMFIT_RHO_SOURCE/data-file-schema.xml"
sha256sum "$STMFIT_RHO_SOURCE/data-file-schema.xml" "$STMFIT_RHO_SOURCE/wfc1.dat" \
    "$STMFIT_RHO_SOURCE/charge-density.dat" > original_before.sha256
trap 'sha256sum -c --quiet original_before.sha256' EXIT
mkdir -p "work/${STMFIT_RHO_STATE}_central.save"
# Only wavefunctions plus their metadata. No old density or mixing history.
cp -p "$STMFIT_RHO_SOURCE/data-file-schema.xml" "$STMFIT_RHO_SOURCE/wfc1.dat" "work/${STMFIT_RHO_STATE}_central.save/"
sha256sum -c --quiet warmstart.sha256
module purge
module load intel/2024.0 impi/2021.11 qe/7.4.1 mkl/2024.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
sha256sum "$(command -v pw.x)" "$STMFIT_RHO_PP" > executables.sha256
srun -n 8 --cpu-bind=cores pw.x -in pw_scf.in > pw_scf.out
"$JULIA_BIN" --startup-file=no --project="$STMFIT_RHO_ROOT" \
    "$STMFIT_RHO_ROOT/test/qe_density_cutoff.jl" check-scf "$STMFIT_RHO_RUN" "$STMFIT_RHO_STATE"
find work -type f -print0 | sort -z | xargs -0 sha256sum > converged_before_pp.sha256
for stem in density total electrostatic repeat_total; do
    srun -n 8 --cpu-bind=cores "$STMFIT_RHO_PP" -in "$stem.in" > "$stem.out"
    grep -q 'JOB DONE.' "$stem.out"
done
sha256sum -c --quiet converged_before_pp.sha256
"$JULIA_BIN" --startup-file=no --project="$STMFIT_RHO_ROOT" \
    "$STMFIT_RHO_ROOT/test/qe_density_cutoff.jl" analyze "$STMFIT_RHO_ROOT" "$STMFIT_RHO_RUN" "$STMFIT_RHO_STATE"
"$JULIA_BIN" --startup-file=no --project="$STMFIT_RHO_ROOT" \
    "$STMFIT_RHO_ROOT/test/test_qe_density_cutoff.jl" "$STMFIT_RHO_ROOT" "$STMFIT_RHO_RUN" "$STMFIT_RHO_STATE" > independent_verification.out 2>&1
sha256sum -c --quiet input_hashes.sha256
echo 'DENSITY_CUTOFF_SENSITIVITY_COMPLETE; no mold or potential adopted'

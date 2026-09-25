#!/usr/bin/env bash
# One isolated build-and-PP experiment; no SCF or geometry optimization.
#SBATCH --job-name=qe-gamma-check
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_gamma_%j.out
#SBATCH --error=qe_gamma_%j.err
set -euo pipefail
module purge
module load "${QE_COMPILER_MODULE:-intel/2024.0}" "${QE_MPI_MODULE:-impi/2021.11}" "${QE_MODULE:-qe/7.4.1}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
QE_NTASKS=${SLURM_NTASKS:?Slurm allocation required}
[[ "$QE_NTASKS" == 8 ]] || exit 1
: "${JULIA_BIN:?Explicit Julia 1.13 binary required}"
STMFIT_GAMMA_ROOT=$(cd ../.. && pwd)
STMFIT_GAMMA_RUN=$PWD
STMFIT_QE_BUILD=$STMFIT_GAMMA_ROOT/qe/gamma_source_20260925
STMFIT_GLCN_SOURCE=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save
STMFIT_GLCNAC_SOURCE=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save
test ! -e work
test ! -e binaries
test ! -e "$STMFIT_QE_BUILD/make.inc"
test "$(git -C "$STMFIT_QE_BUILD" rev-parse HEAD)" = 500de340b820e1cb8c05f2d8bb8fced102f377c1
git -C "$STMFIT_QE_BUILD" diff --exit-code
git -C "$STMFIT_QE_BUILD" diff --cached --exit-code
test "$(git -C "$STMFIT_QE_BUILD/external/fox" rev-parse HEAD)" = 3453648e6837658b747b895bb7bef4b1ed2eac40
test "$(git -C "$STMFIT_QE_BUILD/external/mbd" rev-parse HEAD)" = 89a3cc199c0a200c9f0f688c3229ef6b9a8d63bd
test "$(git -C "$STMFIT_QE_BUILD/external/devxlib" rev-parse HEAD)" = a6b89ef77b1ceda48e967921f1f5488d2df9226d
sha256sum "$STMFIT_GLCN_SOURCE/wfc1.dat" "$STMFIT_GLCNAC_SOURCE/wfc1.dat" > original_wfc_before.sha256
mkdir -p work/glcn work/glcnac binaries build_logs
cp -a "$STMFIT_GLCN_SOURCE" work/glcn/
cp -a "$STMFIT_GLCNAC_SOURCE" work/glcnac/
sha256sum -c --quiet collected_wfc.sha256
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_before.sha256
for molecule in glcn glcnac; do
    cmp "work/$molecule/${molecule}_central.save/data-file-schema.xml" "$molecule/data-file-schema.xml"
    cd "$molecule"
    srun -n "$QE_NTASKS" --cpu-bind=cores pp.x -in pp_ldos.in > site_control.out
    grep -q 'JOB DONE.' site_control.out
    cmp site_control.cube reference.cube
    cd "$STMFIT_GAMMA_RUN"
done
cd "$STMFIT_QE_BUILD"
# Dependencies were fetched at their pinned commits before the job. Forbid
# transports here; a missing dependency must fail rather than fetch on compute.
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=protocol.allow GIT_CONFIG_VALUE_0=never
F90=ifort MPIF90=mpiifort CC=icx BLAS_LIBS='-mkl=sequential' \
    ./configure --enable-parallel --disable-openmp --with-scalapack=no --with-fox=yes \
    > "$STMFIT_GAMMA_RUN/build_logs/configure.out" 2>&1
make -j 8 pwlibs > "$STMFIT_GAMMA_RUN/build_logs/stock_libraries.out" 2>&1
make -j 8 -C PP/src pp.x > "$STMFIT_GAMMA_RUN/build_logs/stock_pp.out" 2>&1
cp -p PP/src/pp.x "$STMFIT_GAMMA_RUN/binaries/pp_stock.x"
cp make.inc "$STMFIT_GAMMA_RUN/build_logs/make.inc"
cp PP/src/local_dos.f90 "$STMFIT_GAMMA_RUN/build_logs/local_dos.stock.f90"
git apply --check "$STMFIT_GAMMA_ROOT/hpc/patches/qe-7.4.1-gamma-density.patch"
git apply "$STMFIT_GAMMA_ROOT/hpc/patches/qe-7.4.1-gamma-density.patch"
git diff -- PP/src/local_dos.f90 > "$STMFIT_GAMMA_RUN/build_logs/applied.patch"
make -j 8 -C PP/src pp.x > "$STMFIT_GAMMA_RUN/build_logs/patched_pp.out" 2>&1
cp -p PP/src/pp.x "$STMFIT_GAMMA_RUN/binaries/pp_gamma.x"
cp PP/src/local_dos.f90 "$STMFIT_GAMMA_RUN/build_logs/local_dos.gamma.f90"
cd "$STMFIT_GAMMA_RUN"
sha256sum binaries/pp_stock.x binaries/pp_gamma.x build_logs/local_dos.stock.f90 build_logs/local_dos.gamma.f90 > executable_sources.sha256
for molecule in glcn glcnac; do
    cd "$molecule"
    for observable in stock_control stock_ildos; do
        srun -n "$QE_NTASKS" --cpu-bind=cores ../binaries/pp_stock.x -in "$observable.in" > "$observable.out"
        grep -q 'JOB DONE.' "$observable.out"
    done
    for observable in patched_control patched_ildos repeat_ildos; do
        srun -n "$QE_NTASKS" --cpu-bind=cores ../binaries/pp_gamma.x -in "$observable.in" > "$observable.out"
        grep -q 'JOB DONE.' "$observable.out"
    done
    cd "$STMFIT_GAMMA_RUN"
done
sha256sum -c --quiet checkpoint_before.sha256
sha256sum -c --quiet original_wfc_before.sha256
for molecule in glcn glcnac; do
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_GAMMA_ROOT" \
        "$STMFIT_GAMMA_ROOT/test/qe_gamma_reconstruction.jl" compare "$molecule" settings.toml
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_GAMMA_ROOT" \
    "$STMFIT_GAMMA_ROOT/test/verify_qe_gamma_reconstruction.jl" "$STMFIT_GAMMA_RUN" \
    > independent_verification.out 2>&1
echo 'GAMMA_DIAGNOSTIC_COMPLETE; no fit, template, calibration or benchmark produced'

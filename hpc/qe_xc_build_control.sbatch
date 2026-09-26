#!/usr/bin/env bash
# One PP-only paired build. No SCF, threshold change, retry or requeue.
# First measure both control cases; components are analyzed only if both pass.
#SBATCH --job-name=qe-xc-build-pair
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=96000MB
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=qe_xc_pair_%j.out
#SBATCH --error=qe_xc_pair_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Compute allocation required}"
: "${JULIA_BIN:?Explicit Julia 1.13 required}"
[[ ${SLURM_NTASKS:?} == 8 ]] || exit 1
module purge
module load intel/2024.0 impi/2021.11 qe/7.4.1 mkl/2024.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
STMFIT_XCB_ROOT=$(cd ../.. && pwd)
STMFIT_XCB_RUN=$PWD
STMFIT_XCB_BUILD=$STMFIT_XCB_ROOT/qe/xc_paired_source_20260926
STMFIT_MKL_LIBDIR=${MKLROOT:?}/lib/intel64
STMFIT_MKL_LINK_FLAGS=("-L$STMFIT_MKL_LIBDIR" "-Wl,-rpath,$STMFIT_MKL_LIBDIR" -mkl=sequential)
STMFIT_GLCN_SOURCE=/u/oldu/code/STMFit/qe/glcn_restart5/qe_tmp/glcn_central.save
STMFIT_GLCNAC_SOURCE=/u/oldu/code/STMFit_spectral_20260925/STMFit/qe/spectral_window_20260925/work/glcnac/glcnac_central.save
test ! -e work
test ! -e binaries
test ! -e "$STMFIT_XCB_BUILD/make.inc"
test "$(git -C "$STMFIT_XCB_BUILD" rev-parse HEAD)" = 500de340b820e1cb8c05f2d8bb8fced102f377c1
git -C "$STMFIT_XCB_BUILD" diff --exit-code
git -C "$STMFIT_XCB_BUILD" diff --cached --exit-code
test "$(git -C "$STMFIT_XCB_BUILD/external/fox" rev-parse HEAD)" = 3453648e6837658b747b895bb7bef4b1ed2eac40
test "$(git -C "$STMFIT_XCB_BUILD/external/mbd" rev-parse HEAD)" = 89a3cc199c0a200c9f0f688c3229ef6b9a8d63bd
test "$(git -C "$STMFIT_XCB_BUILD/external/devxlib" rev-parse HEAD)" = a6b89ef77b1ceda48e967921f1f5488d2df9226d
mkdir build_logs binaries
mpiifort "$STMFIT_XCB_ROOT/test/fixtures/qe_mkl_link_probe.f90" \
    "${STMFIT_MKL_LINK_FLAGS[@]}" -o build_logs/mkl_link_probe > build_logs/mkl_link.out 2>&1
build_logs/mkl_link_probe > build_logs/mkl_smoke.out 2>&1
sha256sum "$STMFIT_GLCN_SOURCE/wfc1.dat" "$STMFIT_GLCNAC_SOURCE/wfc1.dat" \
    "$STMFIT_GLCN_SOURCE/charge-density.dat" "$STMFIT_GLCNAC_SOURCE/charge-density.dat" > original_state_before.sha256
mkdir -p work/glcn work/glcnac
cp -a "$STMFIT_GLCN_SOURCE" work/glcn/
cp -a "$STMFIT_GLCNAC_SOURCE" work/glcnac/
sha256sum -c --quiet accepted_states.sha256
find work -type f -print0 | sort -z | xargs -0 sha256sum > checkpoint_before.sha256
sha256sum "$(command -v pp.x)" > site_pp_binary.sha256
cd "$STMFIT_XCB_BUILD"
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=protocol.allow GIT_CONFIG_VALUE_0=never
F90=ifort MPIF90=mpiifort CC=icx BLAS_LIBS="${STMFIT_MKL_LINK_FLAGS[*]}" \
    ./configure --enable-parallel --disable-openmp --with-scalapack=no --with-fox=yes \
    > "$STMFIT_XCB_RUN/build_logs/configure.out" 2>&1
cp install/config.log "$STMFIT_XCB_RUN/build_logs/config.log"
grep -q '^configure: exit 0$' install/config.log
cp make.inc "$STMFIT_XCB_RUN/build_logs/make.inc"
grep -Eq '^LAPACK[[:space:]]*=[[:space:]]*$' make.inc
make -j 8 pwlibs > "$STMFIT_XCB_RUN/build_logs/libraries.out" 2>&1
make -j 8 -C PP/src pp.x > "$STMFIT_XCB_RUN/build_logs/stock_pp.out" 2>&1
cp -p PP/src/pp.x "$STMFIT_XCB_RUN/binaries/pp_stock.x"
# Only punch_plot.o and its containing libpp.a may change; the helper is new.
find . -type f \( -name '*.o' -o -name '*.a' -o -name '*.mod' \) \
    ! -path ./PP/src/punch_plot.o ! -path ./PP/src/libpp.a -print0 \
    | sort -z | xargs -0 sha256sum > "$STMFIT_XCB_RUN/build_logs/unchanged_libraries.sha256"
git apply --check "$STMFIT_XCB_ROOT/hpc/patches/qe-7.4.1-xc-components.patch"
git apply "$STMFIT_XCB_ROOT/hpc/patches/qe-7.4.1-xc-components.patch"
cp "$STMFIT_XCB_ROOT/hpc/qe_xc_components.f90" PP/src/stmfit_xc_components.f90
make -j 8 -C PP/src pp.x > "$STMFIT_XCB_RUN/build_logs/extraction_pp.out" 2>&1
sha256sum -c --quiet "$STMFIT_XCB_RUN/build_logs/unchanged_libraries.sha256"
cmp make.inc "$STMFIT_XCB_RUN/build_logs/make.inc"
cp -p PP/src/pp.x "$STMFIT_XCB_RUN/binaries/pp_xc.x"
cp PP/src/stmfit_xc_components.f90 "$STMFIT_XCB_RUN/build_logs/"
git diff -- PP/src/Makefile PP/src/punch_plot.f90 > "$STMFIT_XCB_RUN/build_logs/extraction.patch"
cd "$STMFIT_XCB_RUN"
sha256sum binaries/pp_stock.x binaries/pp_xc.x build_logs/stmfit_xc_components.f90 \
    build_logs/extraction.patch > binary_sources.sha256
for molecule in glcn glcnac; do
    cmp "work/$molecule/${molecule}_central.save/data-file-schema.xml" "$molecule/data-file-schema.xml"
    cd "$molecule/stock"
    for stem in density total electrostatic; do
        srun -n 8 --cpu-bind=cores ../../binaries/pp_stock.x -in "$stem.in" > "$stem.out"
        grep -q 'JOB DONE.' "$stem.out"
    done
    cd ../site
    srun -n 8 --cpu-bind=cores pp.x -in density.in > density.out
    grep -q 'JOB DONE.' density.out
    cd ..
    srun -n 8 --cpu-bind=cores ../binaries/pp_xc.x -in extract.in > extract.out
    grep -q 'XC_COMPONENTS_EXPORTED' extract.out
    grep -q 'JOB DONE.' extract.out
    cd repeat
    srun -n 8 --cpu-bind=cores ../../binaries/pp_xc.x -in extract.in > extract.out
    grep -q 'XC_COMPONENTS_EXPORTED' extract.out
    grep -q 'JOB DONE.' extract.out
    cd "$STMFIT_XCB_RUN"
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_XCB_ROOT" \
        "$STMFIT_XCB_ROOT/test/qe_xc_build_control.jl" report "$STMFIT_XCB_ROOT" . "$molecule"
done
sha256sum -c --quiet checkpoint_before.sha256
sha256sum -c --quiet original_state_before.sha256
"$JULIA_BIN" --startup-file=no --project="$STMFIT_XCB_ROOT" \
    "$STMFIT_XCB_ROOT/test/qe_xc_build_control.jl" check .
for molecule in glcn glcnac; do
    cd "$molecule"
    for field in rho_valence rho_xc_input grad2_xc_input gga_active lda gga_local gga_divergence; do
        srun -n 8 --cpu-bind=cores ../binaries/pp_stock.x -in "cube_$field.in" > "cube_$field.out"
        grep -q 'JOB DONE.' "cube_$field.out"
    done
    cd "$STMFIT_XCB_RUN"
    "$JULIA_BIN" --startup-file=no --project="$STMFIT_XCB_ROOT" \
        "$STMFIT_XCB_ROOT/test/qe_xc_build_control.jl" analyze "$STMFIT_XCB_ROOT" . "$molecule"
done
"$JULIA_BIN" --startup-file=no --project="$STMFIT_XCB_ROOT" \
    "$STMFIT_XCB_ROOT/test/test_qe_xc_build_control.jl" "$STMFIT_XCB_ROOT" "$STMFIT_XCB_RUN" > independent_verification.out 2>&1
sha256sum -c --quiet checkpoint_before.sha256
sha256sum -c --quiet original_state_before.sha256
echo 'XC_PAIRED_BUILD_COMPLETE; cross-build differences measured, no new physical potential adopted'

#!/bin/bash -l
#SBATCH --job-name=stmfit-jp-cal
#SBATCH --output=results/hpc_logs/joint_proxy_calibration_full_%j.out
#SBATCH --error=results/hpc_logs/joint_proxy_calibration_full_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=12000MB
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="${STMFIT_PROJECT_DIR:-/u/${USER}/code/STMFit}"
CONFIG="${STMFIT_JOINT_PROXY_CONFIG:-config/joint_proxy_molds.toml}"
OUT="${STMFIT_CALIBRATION_OUT:-results/joint_proxy_calibration_dft_m030_h050_v1_full_candidate.toml}"
SEED="${STMFIT_CALIBRATION_SEED:-20260721}"
CASES="${STMFIT_CALIBRATION_CASES:-100}"

cd "$PROJECT_DIR"
module purge
module load julia/${JULIA_MODULE_VERSION:-1.12}
# The root environment must be instantiated on the login node before submission;
# compute nodes have no internet access.
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

julia --project=. test/calibrate_joint_proxy_molds.jl \
    --config "$CONFIG" \
    --seed "$SEED" \
    --cases "$CASES" \
    --n-min 2 \
    --n-max 12 \
    --out "$OUT"

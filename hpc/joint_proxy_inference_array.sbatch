#!/bin/bash -l
#SBATCH --job-name=stmfit-jp
#SBATCH --output=results/hpc_logs/joint_proxy_%A_%a.out
#SBATCH --error=results/hpc_logs/joint_proxy_%A_%a.err
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=12000MB
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="${STMFIT_PROJECT_DIR:-/u/${USER}/code/STMFit}"
DATA_DIR="${STMFIT_DATA_DIR:-/ptmp/${USER}/stmfit/data}"
FILES_FROM="${STMFIT_FILES_FROM:-hpc/joint_proxy_files_all.txt}"
CONFIG="${STMFIT_JOINT_PROXY_CONFIG:-config/joint_proxy_molds.toml}"
CALIBRATION="${STMFIT_JOINT_PROXY_CALIBRATION:-config/joint_proxy_calibration_dft_m030_h050_v1.toml}"
OUTBASE="${STMFIT_JOINT_PROXY_OUTBASE:-results/joint_proxy_dft_m030_h050_v1_all}"
N_CHUNKS="${N_CHUNKS:-${SLURM_ARRAY_TASK_MAX:-8}}"

cd "$PROJECT_DIR"
module purge
module load julia/${JULIA_MODULE_VERSION:-1.12}
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

SHARD=$(printf '%s/shard_%02d' "$OUTBASE" "$SLURM_ARRAY_TASK_ID")
mkdir -p "$SHARD" results/hpc_logs
SHARD_CALIBRATION="$SHARD/joint_proxy_calibration.toml"
cp "$CALIBRATION" "$SHARD_CALIBRATION"

julia --project=. test/infer_joint_proxy_molds.jl \
    --config "$CONFIG" \
    --data-dir "$DATA_DIR" \
    --files-from "$FILES_FROM" \
    --calibration "$SHARD_CALIBRATION" \
    --outdir "$SHARD" \
    --chunk "${SLURM_ARRAY_TASK_ID}/${N_CHUNKS}"

julia --project=. test/validate_joint_proxy_predictions.jl \
    --artifacts "$SHARD" \
    --calibration "$SHARD_CALIBRATION"

#!/bin/bash -l
#SBATCH --job-name=stmfit-hua
#SBATCH --output=results/hpc_logs/hierarchical_unit_%A_%a.out
#SBATCH --error=results/hpc_logs/hierarchical_unit_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=8000MB
#SBATCH --mail-type=END,FAIL

set -euo pipefail

: "${STMFIT_PROJECT_DIR:?STMFIT_PROJECT_DIR is required}"
: "${STMFIT_HIERARCHICAL_FEATURES:?STMFIT_HIERARCHICAL_FEATURES is required}"
: "${STMFIT_HIERARCHICAL_CONFIG:?STMFIT_HIERARCHICAL_CONFIG is required}"
: "${STMFIT_HIERARCHICAL_OUTBASE:?STMFIT_HIERARCHICAL_OUTBASE is required}"
: "${STMFIT_HIERARCHICAL_FOLDS:?STMFIT_HIERARCHICAL_FOLDS is required}"
: "${STMFIT_HIERARCHICAL_CONFIG_SHA256:?STMFIT_HIERARCHICAL_CONFIG_SHA256 is required}"
: "${STMFIT_HIERARCHICAL_INPUT_SHA256:?STMFIT_HIERARCHICAL_INPUT_SHA256 is required}"
: "${STMFIT_HIERARCHICAL_SOURCE_SHA256:?STMFIT_HIERARCHICAL_SOURCE_SHA256 is required}"

cd "$STMFIT_PROJECT_DIR"
mkdir -p "$STMFIT_HIERARCHICAL_OUTBASE/shards" results/hpc_logs

module purge
module load julia/${JULIA_MODULE_VERSION:-1.12}
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

IFS=':' read -r -a DATE_FOLDS <<< "$STMFIT_HIERARCHICAL_FOLDS"
BLOCKS_PER_FOLD=11  # one unbootstrapped shard + ten fixed 50-seed blocks
TASK_ZERO=$((SLURM_ARRAY_TASK_ID - 1))
DATE_FOLD_INDEX=$((TASK_ZERO / BLOCKS_PER_FOLD))
BOOTSTRAP_BLOCK_INDEX=$((TASK_ZERO % BLOCKS_PER_FOLD))
(( DATE_FOLD_INDEX < ${#DATE_FOLDS[@]} )) || {
    echo "ERROR: array task does not map to a declared date_fold" >&2
    exit 1
}
DATE_FOLD="${DATE_FOLDS[$DATE_FOLD_INDEX]}"

if (( BOOTSTRAP_BLOCK_INDEX == 0 )); then
    BOOTSTRAP_BLOCK="unbootstrapped"
    BLOCK_TAG="unbootstrapped"
else
    BLOCK_START=$((50 * (BOOTSTRAP_BLOCK_INDEX - 1)))
    BLOCK_STOP=$((BLOCK_START + 49))
    BOOTSTRAP_BLOCK="${BLOCK_START}:${BLOCK_STOP}"
    BLOCK_TAG=$(printf '%03d-%03d' "$BLOCK_START" "$BLOCK_STOP")
fi

OUT=$(printf '%s/shards/fold_%s_block_%s.tsv' \
    "$STMFIT_HIERARCHICAL_OUTBASE" "$DATE_FOLD" "$BLOCK_TAG")

echo "date_fold=$DATE_FOLD bootstrap_block=$BOOTSTRAP_BLOCK task=$SLURM_ARRAY_TASK_ID"
echo "threads=$JULIA_NUM_THREADS output=$OUT"

julia --project=. test/evaluate_hierarchical_unit_assignment.jl \
    --features "$STMFIT_HIERARCHICAL_FEATURES" \
    --config "$STMFIT_HIERARCHICAL_CONFIG" \
    --fold "$DATE_FOLD" \
    --bootstrap-block "$BOOTSTRAP_BLOCK" \
    --config-sha256 "$STMFIT_HIERARCHICAL_CONFIG_SHA256" \
    --input-sha256 "$STMFIT_HIERARCHICAL_INPUT_SHA256" \
    --source-sha256 "$STMFIT_HIERARCHICAL_SOURCE_SHA256" \
    --out "$OUT"

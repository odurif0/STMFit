#!/bin/bash -l
#SBATCH --job-name=stmfit-type-merge
#SBATCH --output=results/hpc_logs/type_ablation_merge_%j.out
#SBATCH --error=results/hpc_logs/type_ablation_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4000MB
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="${STMFIT_PROJECT_DIR:-/u/${USER}/code/STMFit}"
OUTBASE="${STMFIT_TYPE_ABLATION_OUTBASE:-results/joint_proxy_type_ablation_v1}"
N_CHUNKS="${N_CHUNKS:-8}"

cd "$PROJECT_DIR"
module purge
module load julia/${JULIA_MODULE_VERSION:-1.12}
export JULIA_NUM_THREADS=1

ARGS=()
for i in $(seq 1 "$N_CHUNKS"); do
    ARGS+=(--shard "$(printf '%s/shard_%02d' "$OUTBASE" "$i")")
done

julia --project=. test/merge_joint_proxy_type_ablation.jl \
    "${ARGS[@]}" --outdir "$OUTBASE/merged"

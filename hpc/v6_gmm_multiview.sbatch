#!/usr/bin/env bash
# Slurm job: GMM multi-view + mold margin (v6). Last untested combination.
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V5, STMFIT_OUTPUT_DIR,
# JULIA_MODULE_VERSION.

#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module purge
module load "julia/${JULIA_MODULE_VERSION}"
cd "${STMFIT_PROJECT_DIR}"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1

V5="${STMFIT_FEATURES_V5}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"

# reference: 1-view + mold (77.3/35)
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v6_ref.tsv" \
    --view "v_avg=${BASE4},patch_u_asym,mold_margin_avg" --seeds 10 --interactions
# 2-view: add backward-moment view
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v6_2v.tsv" \
    --view "v_avg=${BASE4},patch_u_asym,mold_margin_avg" \
    --view "v_comt=${BASE4},bwd_neg_com_t" --seeds 10 --interactions
# 3-view: historical 3-view + mold
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v6_3v.tsv" \
    --view "v_avg=${BASE4},patch_u_asym,mold_margin_avg" \
    --view "v_comt=${BASE4},bwd_neg_com_t" \
    --view "v_diag45=${BASE4},bwd_neg_diag45" --seeds 10 --interactions
# 2-view with split instead of com_t
predict build_labelfree_gmm_predictions.jl "${V5}" "${OUT}/pred_gmm_v6_2v_split.tsv" \
    --view "v_avg=${BASE4},patch_u_asym,mold_margin_avg" \
    --view "v_split=${BASE4},split_log_skew" --seeds 10 --interactions

echo "V6 SWEEP COMPLETE"
ls -la "${OUT}"

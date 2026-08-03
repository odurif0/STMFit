#!/usr/bin/env bash
# Slurm job: mold-cost feature sweep (v3).
# Tests the DFT mold cost_margin/cost_GlcN/cost_GlcNAc (17x17 contrast
# scoring, independent physical signal) as features for the k-means 4-view
# and GMM 1-view baselines.
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V3, STMFIT_OUTPUT_DIR,
# JULIA_MODULE_VERSION.

#SBATCH --time=03:00:00
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

V3="${STMFIT_FEATURES_V3}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
MARGIN="mold_cost_margin"
MOLD3="mold_cost_margin,mold_cost_glcn,mold_cost_glcnac"
KM_VIEWS=(--view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" \
          --view "v_diag135=${BASE4},bwd_neg_diag135" \
          --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view + mold costs ---
predict build_labelfree_unit_predictions.jl "${V3}" "${OUT}/pred_km_v3_ref.tsv" "${KM_VIEWS[@]}" --interactions
predict build_labelfree_unit_predictions.jl "${V3}" "${OUT}/pred_km_v3_margin.tsv" "${KM_VIEWS[@]}" \
    --view "v_margin=${BASE4},${MARGIN}" --interactions
predict build_labelfree_unit_predictions.jl "${V3}" "${OUT}/pred_km_v3_mold3.tsv" "${KM_VIEWS[@]}" \
    --view "v_mold3=${BASE4},${MOLD3}" --interactions

# --- GMM 1-view + mold costs ---
predict build_labelfree_gmm_predictions.jl "${V3}" "${OUT}/pred_gmm_v3_ref.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V3}" "${OUT}/pred_gmm_v3_margin.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym,${MARGIN}" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V3}" "${OUT}/pred_gmm_v3_margin_ni.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym,${MARGIN}" --seeds 10 --no-interactions

echo "V3 SWEEP COMPLETE"
ls -la "${OUT}"

#!/usr/bin/env bash
# Slurm job: adaptive-contour (constant-current) coupling sweep (v8).
# - k-means 4-view + mold_cc_avg view (the k-means x adaptive-contour coupling)
# - GMM 1-view + mold_cc_avg refinements (selftrain, seeds, no-uasym)
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V7, STMFIT_OUTPUT_DIR,
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

V7="${STMFIT_FEATURES_V7}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
CC="mold_cc_avg"
CH="mold_margin_avg"
KM_VIEWS=(--view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" \
          --view "v_diag135=${BASE4},bwd_neg_diag135" \
          --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view x adaptive contour ---
predict build_labelfree_unit_predictions.jl "${V7}" "${OUT}/pred_km_v8_ref.tsv" "${KM_VIEWS[@]}" --interactions
predict build_labelfree_unit_predictions.jl "${V7}" "${OUT}/pred_km_v8_cc.tsv" "${KM_VIEWS[@]}" \
    --view "v_cc=${BASE4},${CC}" --interactions
predict build_labelfree_unit_predictions.jl "${V7}" "${OUT}/pred_km_v8_ccch.tsv" "${KM_VIEWS[@]}" \
    --view "v_cc=${BASE4},${CC}" --view "v_ch=${BASE4},${CH}" --interactions

# --- GMM 1-view x adaptive contour, refinements ---
predict build_labelfree_gmm_predictions.jl "${V7}" "${OUT}/pred_gmm_v8_cc.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,${CC}" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V7}" "${OUT}/pred_gmm_v8_cc_st2.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,${CC}" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V7}" "${OUT}/pred_gmm_v8_cc_s20.tsv" \
    --view "v_cc=${BASE4},patch_u_asym,${CC}" --seeds 20 --interactions
predict build_labelfree_gmm_predictions.jl "${V7}" "${OUT}/pred_gmm_v8_cc_nouasym.tsv" \
    --view "v_cc=${BASE4},${CC}" --seeds 10 --interactions

echo "V8 SWEEP COMPLETE"
ls -la "${OUT}"

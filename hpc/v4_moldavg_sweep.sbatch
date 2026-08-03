#!/usr/bin/env bash
# Slurm job: mold-averaged-cost feature sweep (v4).
# Tests mold_margin_avg (forward+backward mean of the 17x17 contrast
# cost_margin) as a feature for k-means 4-view and GMM 1-view.
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V4, STMFIT_OUTPUT_DIR,
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

V4="${STMFIT_FEATURES_V4}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
KM_VIEWS=(--view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" \
          --view "v_diag135=${BASE4},bwd_neg_diag135" \
          --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view + mold margins ---
predict build_labelfree_unit_predictions.jl "${V4}" "${OUT}/pred_km_v4_ref.tsv" "${KM_VIEWS[@]}" --interactions
predict build_labelfree_unit_predictions.jl "${V4}" "${OUT}/pred_km_v4_avg.tsv" "${KM_VIEWS[@]}" \
    --view "v_avg=${BASE4},mold_margin_avg" --interactions
predict build_labelfree_unit_predictions.jl "${V4}" "${OUT}/pred_km_v4_avg3.tsv" "${KM_VIEWS[@]}" \
    --view "v_avg=${BASE4},mold_margin_avg,mold_margin_fwd,mold_margin_bwd" --interactions

# --- GMM 1-view + mold margins ---
predict build_labelfree_gmm_predictions.jl "${V4}" "${OUT}/pred_gmm_v4_ref.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V4}" "${OUT}/pred_gmm_v4_avg.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym,mold_margin_avg" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V4}" "${OUT}/pred_gmm_v4_avg_ni.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym,mold_margin_avg" --seeds 10 --no-interactions

echo "V4 SWEEP COMPLETE"
ls -la "${OUT}"

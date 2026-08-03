#!/usr/bin/env bash
# Slurm job: physical-feature sweep v2b.
# - k-means: 4-view baseline WITH --interactions (fixes the v2 omission) plus
#   one reduced family view per config (3-5 features to avoid k-means dilution).
# - GMM 1-view: baseline with/without interactions; family views run
#   --no-interactions to avoid full-covariance over-parameterization.
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V2, STMFIT_OUTPUT_DIR,
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

V2="${STMFIT_FEATURES_V2}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"

predict() {
    local script="$1" feats="$2" out="$3"
    shift 3
    echo "=== $(basename "$out") :: $script :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}

BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
# reduced physical families (kept small: k-means views stay close to 15 dims)
BWD3="bwd_res_en,bwd_res_com_u,bwd_res_asym_u"
DIFF5="diff_res_en,diff_res_com_u,diff_res_asym_u,diff_res_sku,diff_res_com_t"
GEOM4="spacing_next_nm,spacing_asym,chain_curv_deg,elongation"
KM_VIEWS=(--view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" \
          --view "v_diag135=${BASE4},bwd_neg_diag135" \
          --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view (default views = historical reference) ---
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2b_ref.tsv" "${KM_VIEWS[@]}" --interactions
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2b_bwdres.tsv" "${KM_VIEWS[@]}" \
    --view "v_bwd3=${BASE4},${BWD3}" --interactions
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2b_diffres.tsv" "${KM_VIEWS[@]}" \
    --view "v_diff5=${BASE4},${DIFF5}" --interactions
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2b_geom.tsv" "${KM_VIEWS[@]}" \
    --view "v_geom4=${BASE4},${GEOM4}" --interactions
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2b_all.tsv" "${KM_VIEWS[@]}" \
    --view "v_bwd3=${BASE4},${BWD3}" --view "v_diff5=${BASE4},${DIFF5}" \
    --view "v_geom4=${BASE4},${GEOM4}" --interactions

# --- GMM 1-view ---
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2b_ref.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2b_ref_ni.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --seeds 10 --no-interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2b_diffres_ni.tsv" \
    --view "v_diff5=${BASE4},patch_u_asym,${DIFF5}" --seeds 10 --no-interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2b_bwdres_ni.tsv" \
    --view "v_bwd3=${BASE4},patch_u_asym,${BWD3}" --seeds 10 --no-interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2b_geom_ni.tsv" \
    --view "v_geom4=${BASE4},patch_u_asym,${GEOM4}" --seeds 10 --no-interactions

echo "V2B SWEEP COMPLETE"
ls -la "${OUT}"

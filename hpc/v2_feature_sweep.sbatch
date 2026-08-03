#!/usr/bin/env bash
# Slurm job: physical-feature sweep on the enriched table v2.
# Tests each new feature family added to the k-means 4-view and GMM 1-view
# baselines. Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V2,
# STMFIT_OUTPUT_DIR, JULIA_MODULE_VERSION.

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
BWD_RES="bwd_res_en,bwd_res_com_u,bwd_res_com_t,bwd_res_sku,bwd_res_skt,bwd_res_asym_u,bwd_res_asym_t"
DIFF_RES="diff_res_en,diff_res_com_u,diff_res_com_t,diff_res_sku,diff_res_skt,diff_res_asym_u,diff_res_asym_t,diff_res_kurt_u,diff_res_ring1,diff_res_ring2,diff_res_ring3"
DIFF_RAW="diff_raw_en,diff_raw_com_u,diff_raw_com_t,diff_raw_sku,diff_raw_skt,diff_raw_asym_u,diff_raw_asym_t"
GEOM="spacing_next_nm,spacing_asym,chain_curv_deg,elongation"
KM_VIEWS=(--view "v_base=${BASE4}" --view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view reference on v2 + each family ---
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_ref.tsv" "${KM_VIEWS[@]}"
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_bwdres.tsv" "${KM_VIEWS[@]}" \
    --view "v_bwdres=${BASE4},${BWD_RES}"
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_diffres.tsv" "${KM_VIEWS[@]}" \
    --view "v_diffres=${BASE4},${DIFF_RES}"
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_diffraw.tsv" "${KM_VIEWS[@]}" \
    --view "v_diffraw=${BASE4},${DIFF_RAW}"
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_geom.tsv" "${KM_VIEWS[@]}" \
    --view "v_geom=${BASE4},${GEOM}"
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2_all.tsv" "${KM_VIEWS[@]}" \
    --view "v_bwdres=${BASE4},${BWD_RES}" --view "v_diffres=${BASE4},${DIFF_RES}" \
    --view "v_geom=${BASE4},${GEOM}"

# --- GMM 1-view reference + families ---
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2_ref.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2_diffres.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --view "v_diffres=${BASE4},${DIFF_RES}" \
    --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2_bwdres.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --view "v_bwdres=${BASE4},${BWD_RES}" \
    --seeds 10 --interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2_geom.tsv" \
    --view "v_uasym=${BASE4},patch_u_asym" --view "v_geom=${BASE4},${GEOM}" \
    --seeds 10 --interactions

echo "V2 SWEEP COMPLETE"
ls -la "${OUT}"

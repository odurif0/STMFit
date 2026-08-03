#!/usr/bin/env bash
# Slurm job: physical-feature sweep v2c — envelope-residual family
# (amp_env_resid_linear, amp_env_resid_quadratic, amp_second_diff), the last
# untested physical features in the table.
# Expected env: STMFIT_PROJECT_DIR, STMFIT_FEATURES_V2, STMFIT_OUTPUT_DIR,
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
ENV3="amp_env_resid_linear,amp_env_resid_quadratic,amp_second_diff"
BWD3="bwd_res_en,bwd_res_com_u,bwd_res_asym_u"
KM_VIEWS=(--view "v_comt=${BASE4},bwd_neg_com_t" \
          --view "v_diag45=${BASE4},bwd_neg_diag45" \
          --view "v_diag135=${BASE4},bwd_neg_diag135" \
          --view "v_split=${BASE4},split_log_skew")

# --- k-means 4-view + envelope-residual family ---
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2c_env.tsv" "${KM_VIEWS[@]}" \
    --view "v_env3=${BASE4},${ENV3}" --interactions
predict build_labelfree_unit_predictions.jl "${V2}" "${OUT}/pred_km_v2c_envbwd.tsv" "${KM_VIEWS[@]}" \
    --view "v_env3=${BASE4},${ENV3}" --view "v_bwd3=${BASE4},${BWD3}" --interactions

# --- GMM 1-view + envelope-residual family (low-dim, no interactions) ---
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2c_env_ni.tsv" \
    --view "v_env3=${BASE4},patch_u_asym,${ENV3}" --seeds 10 --no-interactions
predict build_labelfree_gmm_predictions.jl "${V2}" "${OUT}/pred_gmm_v2c_envbwd_ni.tsv" \
    --view "v_envbwd=${BASE4},patch_u_asym,${ENV3},${BWD3}" --seeds 10 --no-interactions

echo "V2C SWEEP COMPLETE"
ls -la "${OUT}"

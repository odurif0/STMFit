#!/usr/bin/env bash
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
V="${STMFIT_FEATURES_V13}"
OUT="${STMFIT_OUTPUT_DIR}"
mkdir -p "${OUT}"
predict() {
    local script="$1" feats="$2" out="$3"; shift 3
    echo "=== $(basename "$out") :: $* ==="
    julia --project=. "test/${script}" --features "${feats}" --out "${out}" "$@"
}
BASE4="amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v13_margin.tsv" \
    --view "v_m=${BASE4},patch_u_asym,def_margin_fwd,def_margin_bwd" --seeds 10 --interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v13_shape_ni.tsv" \
    --view "v_s=${BASE4},patch_u_asym,def_st0_f,def_su0_f,def_st1_f,def_su1_f" --seeds 10 --no-interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v13_all_ni.tsv" \
    --view "v_a=${BASE4},patch_u_asym,def_margin_fwd,def_margin_bwd,def_st0_f,def_su0_f,def_st1_f,def_su1_f" \
    --seeds 10 --no-interactions --selftrain 2
predict build_labelfree_gmm_predictions.jl "${V}" "${OUT}/pred_gmm_v13_margin_b.tsv" \
    --view "v_m=${BASE4},patch_u_asym,def_margin_fwd,def_margin_bwd,def_st0_b,def_su0_b" \
    --seeds 10 --interactions --selftrain 2
echo "V13 SWEEP COMPLETE"
ls -la "${OUT}"

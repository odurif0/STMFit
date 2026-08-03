#!/usr/bin/env bash
#SBATCH --job-name=stmfit-frozen-contrast
#SBATCH --output=/ptmp/%u/stmfit/frozen_contrast/logs/%A_%a.out
#SBATCH --error=/ptmp/%u/stmfit/frozen_contrast/logs/%A_%a.err

set -euo pipefail

: "${STMFIT_REMOTE_PROJECT:?missing STMFIT_REMOTE_PROJECT}"
: "${STMFIT_DATA_DIR:?missing STMFIT_DATA_DIR}"
: "${STMFIT_WHOLE_ROI_MOLDS:?missing STMFIT_WHOLE_ROI_MOLDS}"
: "${STMFIT_WHOLE_ROI_SUMMARY:?missing STMFIT_WHOLE_ROI_SUMMARY}"
: "${STMFIT_WHOLE_ROI_CONFIG:?missing STMFIT_WHOLE_ROI_CONFIG}"
: "${STMFIT_WHOLE_ROI_FILE_LIST:?missing STMFIT_WHOLE_ROI_FILE_LIST}"
: "${STMFIT_WHOLE_ROI_OUTDIR:?missing STMFIT_WHOLE_ROI_OUTDIR}"

for var in \
    STMFIT_WHOLE_ROI_TRUTH STMFIT_WHOLE_ROI_GRADE \
    STMFIT_WHOLE_ROI_EXPECTED_N STMFIT_WHOLE_ROI_CONTROL_SEQUENCE \
    STMFIT_WHOLE_ROI_MANIFEST STMFIT_WHOLE_ROI_BENCHMARK_MANIFEST \
    STMFIT_WHOLE_ROI_FULL145 STMFIT_TRUTH STMFIT_GRADE \
    STMFIT_EXPECTED_N STMFIT_CONTROL_SEQUENCE STMFIT_MANIFEST
do
    if [[ -n "${!var:-}" ]]; then
        echo "forbidden benchmark/truth input at array boundary: $var" >&2
        exit 2
    fi
done

# SXM file-list entries are transported and used as local basenames. Keep this
# contract deliberately narrow and identical to the launcher: no whitespace,
# quotes, shell metacharacters, leading dashes, command substitutions, or path
# separators.
safe_sxm_basename_re='^[A-Za-z0-9][A-Za-z0-9._-]*[.][sS][xX][mM]$'
unsafe_sxm_basename() {
    echo "unsafe SXM basename contains characters outside the launcher transport contract: $1" >&2
    exit 2
}

files=()
while IFS= read -r raw || [[ -n "$raw" ]]; do
    line="${raw%%#*}"
    trimmed="$(printf '%s' "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    [[ -z "$trimmed" ]] && continue
    [[ "$line" == "$trimmed" ]] || unsafe_sxm_basename "$raw"
    [[ "$trimmed" =~ $safe_sxm_basename_re ]] || unsafe_sxm_basename "$trimmed"
    files+=("$trimmed")
done < "$STMFIT_WHOLE_ROI_FILE_LIST"
(( ${#files[@]} > 0 )) || {
    echo "empty STMFIT_WHOLE_ROI_FILE_LIST: $STMFIT_WHOLE_ROI_FILE_LIST" >&2
    exit 2
}

index=$((SLURM_ARRAY_TASK_ID - 1))
(( index >= 0 && index < ${#files[@]} )) || {
    echo "invalid array index: $SLURM_ARRAY_TASK_ID for ${#files[@]} files" >&2
    exit 2
}
file="${files[$index]}"
stem="${file%.sxm}"

julia_threads="${STMFIT_JULIA_THREADS_PER_TASK:-${SLURM_CPUS_PER_TASK:-4}}"
[[ "$julia_threads" =~ ^[0-9]+$ ]] || {
    echo "STMFIT_JULIA_THREADS_PER_TASK must be an integer" >&2
    exit 2
}
(( julia_threads >= 1 && julia_threads <= 4 )) || {
    echo "useful Julia threads must stay in 1..4/task, got $julia_threads" >&2
    exit 2
}
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    (( julia_threads <= SLURM_CPUS_PER_TASK )) || {
        echo "Julia threads ($julia_threads) exceed SLURM_CPUS_PER_TASK ($SLURM_CPUS_PER_TASK)" >&2
        exit 2
    }
fi

if [[ "$(basename "$STMFIT_WHOLE_ROI_CONFIG")" == "joint_proxy_whole_roi_constant_current.toml" ]]; then
    : "${STMFIT_CONSTANT_CURRENT_PROVENANCE:?missing STMFIT_CONSTANT_CURRENT_PROVENANCE}"
    : "${STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256:?missing STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256}"
    : "${STMFIT_CONSTANT_CURRENT_MAPS:?missing STMFIT_CONSTANT_CURRENT_MAPS}"
    : "${STMFIT_CONSTANT_CURRENT_INVALID_MASK:?missing STMFIT_CONSTANT_CURRENT_INVALID_MASK}"
    : "${STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE:?missing STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE}"
fi

module purge
if [[ -n "${JULIA_MODULE_VERSION:-}" ]]; then
    module load "julia/${JULIA_MODULE_VERSION}"
else
    module load julia
fi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GKSwstype=100

mkdir -p "$STMFIT_WHOLE_ROI_OUTDIR/$stem"
cd "$STMFIT_REMOTE_PROJECT"
julia --threads="$julia_threads" --project=. \
    test/plot_frozen_contrast_diagnostic.jl \
    --out-dir "$STMFIT_WHOLE_ROI_OUTDIR/$stem" \
    --selection-summary "$STMFIT_WHOLE_ROI_SUMMARY" \
    --molds "$STMFIT_WHOLE_ROI_MOLDS" \
    --config "$STMFIT_WHOLE_ROI_CONFIG" \
    --observable-provenance "${STMFIT_CONSTANT_CURRENT_PROVENANCE:-}" \
    --observable-provenance-sha256 "${STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256:-}" \
    --observable-map "${STMFIT_CONSTANT_CURRENT_MAPS:-}" \
    --validity-mask "${STMFIT_CONSTANT_CURRENT_INVALID_MASK:-}" \
    --observable-mold-provenance "${STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE:-}" \
    --half-nm 0.80 --step-nm 0.08 "$file"

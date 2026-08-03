#!/usr/bin/env bash
# Dedicated launcher for the heavyweight frozen-contrast diagnostic.
# Always run --dry-run before --watch.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_FILE="$SCRIPT_DIR/remote.env"

DRY_RUN=0
WATCH=0
case "${1:-}" in
    --dry-run) DRY_RUN=1 ;;
    --submit) ;;
    --watch) WATCH=1 ;;
    -h|--help)
        echo "Usage: $0 --dry-run|--submit|--watch"
        exit 0 ;;
    "") echo "Usage: $0 --dry-run|--submit|--watch" >&2; exit 2 ;;
    --truth|--truth=*|--grade|--grade=*|--expected-N|--expected-N=*|\
    --expected-n|--expected-n=*|--control-sequence|--control-sequence=*|\
    --manifest|--manifest=*|--benchmark-manifest|--benchmark-manifest=*|\
    --full145|--full145=*|--control|--control=*)
        echo "forbidden benchmark/truth input at launcher boundary: ${1:-}" >&2
        exit 2 ;;
    *) echo "Usage: $0 --dry-run|--submit|--watch" >&2; exit 2 ;;
esac

if [[ -f "$ENV_FILE" ]]; then
    # shellcheck source=/dev/null
    source "$ENV_FILE"
elif (( ! DRY_RUN )); then
    echo "missing $ENV_FILE" >&2
    exit 1
fi

for var in \
    STMFIT_WHOLE_ROI_TRUTH STMFIT_WHOLE_ROI_GRADE \
    STMFIT_WHOLE_ROI_EXPECTED_N STMFIT_WHOLE_ROI_CONTROL_SEQUENCE \
    STMFIT_WHOLE_ROI_MANIFEST STMFIT_WHOLE_ROI_BENCHMARK_MANIFEST \
    STMFIT_WHOLE_ROI_FULL145 STMFIT_TRUTH STMFIT_GRADE \
    STMFIT_EXPECTED_N STMFIT_CONTROL_SEQUENCE STMFIT_MANIFEST
do
    if [[ -n "${!var:-}" ]]; then
        echo "forbidden benchmark/truth input at launcher boundary: $var" >&2
        exit 2
    fi
done

if [[ -z "${STMFIT_SSH_HOST:-}" ]]; then
    if (( DRY_RUN )); then
        STMFIT_SSH_HOST="dry-run-no-remote"
    else
        echo "missing STMFIT_SSH_HOST" >&2
        exit 1
    fi
fi
if [[ -z "${STMFIT_REMOTE_USER:-}" ]]; then
    if (( DRY_RUN )); then
        STMFIT_REMOTE_USER="${USER:-$(id -un)}"
    else
        echo "missing STMFIT_REMOTE_USER" >&2
        exit 1
    fi
fi

: "${STMFIT_LOCAL_DATA:=${STMFIT_DATA_DIR:-}}"
: "${REMOTE_PROJECT_BASE:=/u/$STMFIT_REMOTE_USER/code}"
: "${REMOTE_DATA_BASE:=/ptmp/$STMFIT_REMOTE_USER/stmfit}"
: "${CPUS_PER_TASK:=4}"
: "${MEM_PER_CPU:=4000}"
: "${WALLTIME:=08:00:00}"
: "${JULIA_MODULE_VERSION:=}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"

[[ "$CPUS_PER_TASK" =~ ^[0-9]+$ ]] || { echo "CPUS_PER_TASK must be an integer" >&2; exit 2; }
(( CPUS_PER_TASK >= 1 && CPUS_PER_TASK <= 4 )) || {
    echo "CPUS_PER_TASK must be in 1..4 so useful Julia threads stay <=4/task" >&2
    exit 2
}
[[ "$MEM_PER_CPU" =~ ^[0-9]+$ ]] || { echo "MEM_PER_CPU must be an integer" >&2; exit 2; }
(( MEM_PER_CPU > 0 )) || { echo "MEM_PER_CPU must be positive" >&2; exit 2; }

safe_path_re='^[-_./+A-Za-z0-9=@]+$'
safe_host_re='^[-_.:@A-Za-z0-9]+$'
safe_token_re='^[-_./+A-Za-z0-9=@:]+$'
# SXM file-list entries are transported as remote basenames. Keep this
# deliberately narrower than POSIX filenames: no whitespace, quotes, shell
# metacharacters, leading dashes, command substitutions, or path separators.
safe_sxm_basename_re='^[A-Za-z0-9][A-Za-z0-9._-]*[.][sS][xX][mM]$'

unsafe_value() {
    local label=$1
    local value=$2
    echo "unsafe $label contains characters outside the launcher transport contract: $value" >&2
    exit 2
}

require_safe_path() {
    local label=$1
    local value=$2
    [[ -z "$value" ]] && return
    [[ "$value" =~ $safe_path_re ]] || unsafe_value "$label" "$value"
}

require_safe_token() {
    local label=$1
    local value=$2
    [[ -z "$value" ]] && return
    [[ "$value" =~ $safe_token_re ]] || unsafe_value "$label" "$value"
}

require_safe_host() {
    local label=$1
    local value=$2
    [[ "$value" =~ $safe_host_re ]] || unsafe_value "$label" "$value"
}

require_safe_sxm_basename() {
    local value=$1
    [[ "$value" =~ $safe_sxm_basename_re ]] || unsafe_value "SXM basename" "$value"
}

shell_quote() {
    local value=$1
    printf "'%s'" "${value//\'/\'\\\'\'}"
}

REMOTE_PROJECT="$REMOTE_PROJECT_BASE/STMFit"
REMOTE_DATA="$REMOTE_DATA_BASE/data"
REMOTE_INPUT="$REMOTE_DATA_BASE/frozen_contrast_input"
REMOTE_OUT="${STMFIT_WHOLE_ROI_REMOTE_OUT:-${STMFIT_WHOLE_ROI_OUTDIR:-$REMOTE_DATA_BASE/frozen_contrast}}"
LOCAL_OUT="${STMFIT_WHOLE_ROI_LOCAL_OUT:-/tmp/opencode/frozen_contrast_hpc}"
LOCAL_SUMMARY="${STMFIT_WHOLE_ROI_SUMMARY:-$REPO_ROOT/results/best_plots_240817_primary_rerun/summary_overlap060_hard.tsv}"
LOCAL_MOLDS="${STMFIT_WHOLE_ROI_MOLDS:-/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv}"
LOCAL_CONFIG="${STMFIT_WHOLE_ROI_CONFIG:-config/joint_proxy_whole_roi_diagnostic.toml}"
LOCAL_FILE_LIST="${STMFIT_WHOLE_ROI_FILE_LIST:-}"

if [[ -z "$LOCAL_FILE_LIST" ]]; then
    LOCAL_FILE_LIST="$LOCAL_OUT/predeclared_two_scan_file_list.txt"
    mkdir -p "$(dirname "$LOCAL_FILE_LIST")"
    printf '%s\n' "240817_007.sxm" "240817_050.sxm" > "$LOCAL_FILE_LIST"
fi

read_file_list() {
    local file_list=$1
    [[ -f "$file_list" ]] || { echo "file list not found: $file_list" >&2; exit 1; }
    FILES=()
    local raw line trimmed
    while IFS= read -r raw || [[ -n "$raw" ]]; do
        line="${raw%%#*}"
        trimmed="$(printf '%s' "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
        [[ -z "$trimmed" ]] && continue
        [[ "$line" == "$trimmed" ]] || unsafe_value "SXM basename" "$raw"
        require_safe_sxm_basename "$trimmed"
        FILES+=("$trimmed")
    done < "$file_list"
    (( ${#FILES[@]} > 0 )) || { echo "file list is empty: $file_list" >&2; exit 2; }
}
read_file_list "$LOCAL_FILE_LIST"
array_size=${#FILES[@]}
max_by_cpu=$((8 / CPUS_PER_TASK))
(( max_by_cpu >= 1 )) || { echo "CPU budget would exceed 8 total CPUs" >&2; exit 2; }
requested_concurrency="${STMFIT_WHOLE_ROI_MAX_CONCURRENT_TASKS:-$max_by_cpu}"
[[ "$requested_concurrency" =~ ^[0-9]+$ ]] || {
    echo "STMFIT_WHOLE_ROI_MAX_CONCURRENT_TASKS must be an integer" >&2
    exit 2
}
(( requested_concurrency >= 1 )) || {
    echo "STMFIT_WHOLE_ROI_MAX_CONCURRENT_TASKS must be positive" >&2
    exit 2
}
(( requested_concurrency <= max_by_cpu )) || {
    echo "requested concurrency $requested_concurrency exceeds 8-CPU budget with CPUS_PER_TASK=$CPUS_PER_TASK" >&2
    exit 2
}
max_concurrency=$requested_concurrency
(( max_concurrency <= array_size )) || max_concurrency=$array_size
JULIA_THREADS_PER_TASK="${STMFIT_JULIA_THREADS_PER_TASK:-$CPUS_PER_TASK}"
[[ "$JULIA_THREADS_PER_TASK" =~ ^[0-9]+$ ]] || {
    echo "STMFIT_JULIA_THREADS_PER_TASK must be an integer" >&2
    exit 2
}
(( JULIA_THREADS_PER_TASK >= 1 && JULIA_THREADS_PER_TASK <= 4 )) || {
    echo "STMFIT_JULIA_THREADS_PER_TASK must be in 1..4" >&2
    exit 2
}
(( JULIA_THREADS_PER_TASK <= CPUS_PER_TASK )) || {
    echo "STMFIT_JULIA_THREADS_PER_TASK cannot exceed CPUS_PER_TASK" >&2
    exit 2
}

remote_stage_path() {
    local path=$1
    [[ -z "$path" ]] && return 0
    printf '%s/%s' "$REMOTE_INPUT" "$(basename "$path")"
}

REMOTE_MOLDS="${STMFIT_REMOTE_WHOLE_ROI_MOLDS:-$(remote_stage_path "$LOCAL_MOLDS")}"
REMOTE_SUMMARY="$(remote_stage_path "$LOCAL_SUMMARY")"
REMOTE_FILE_LIST="$(remote_stage_path "$LOCAL_FILE_LIST")"
if [[ "$LOCAL_CONFIG" = /* ]]; then
    REMOTE_CONFIG="$(remote_stage_path "$LOCAL_CONFIG")"
else
    REMOTE_CONFIG="$LOCAL_CONFIG"
fi

LOCAL_CC_PROVENANCE="${STMFIT_CONSTANT_CURRENT_PROVENANCE:-}"
LOCAL_CC_MAPS="${STMFIT_CONSTANT_CURRENT_MAPS:-}"
LOCAL_CC_MASK="${STMFIT_CONSTANT_CURRENT_INVALID_MASK:-}"
LOCAL_CC_MOLD_PROVENANCE="${STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE:-}"
CC_PROVENANCE_SHA="${STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256:-}"
if [[ "$(basename "$LOCAL_CONFIG")" == "joint_proxy_whole_roi_constant_current.toml" && -z "$LOCAL_CC_MOLD_PROVENANCE" && -n "$LOCAL_MOLDS" ]]; then
    LOCAL_CC_MOLD_PROVENANCE="${LOCAL_MOLDS}.provenance.toml"
fi
REMOTE_CC_PROVENANCE="$(remote_stage_path "$LOCAL_CC_PROVENANCE")"
REMOTE_CC_MAPS="$(remote_stage_path "$LOCAL_CC_MAPS")"
REMOTE_CC_MASK="$(remote_stage_path "$LOCAL_CC_MASK")"
REMOTE_CC_MOLD_PROVENANCE="$(remote_stage_path "$LOCAL_CC_MOLD_PROVENANCE")"

LOCAL_CC_EXTRA_ARTIFACTS=()
if [[ "$(basename "$LOCAL_CONFIG")" == "joint_proxy_whole_roi_constant_current.toml" && -n "$LOCAL_CC_PROVENANCE" ]]; then
    cc_prefix="${LOCAL_CC_PROVENANCE%.provenance.toml}"
    shopt -s nullglob
    for artifact in "${cc_prefix}"_h*.tsv "${cc_prefix}"_h*.mask.tsv; do
        [[ -f "$artifact" ]] && LOCAL_CC_EXTRA_ARTIFACTS+=("$artifact")
    done
    shopt -u nullglob
fi

require_safe_host "STMFIT_SSH_HOST" "$STMFIT_SSH_HOST"
require_safe_token "STMFIT_REMOTE_USER" "$STMFIT_REMOTE_USER"
require_safe_path "REMOTE_PROJECT_BASE" "$REMOTE_PROJECT_BASE"
require_safe_path "REMOTE_DATA_BASE" "$REMOTE_DATA_BASE"
require_safe_path "REMOTE_PROJECT" "$REMOTE_PROJECT"
require_safe_path "REMOTE_DATA" "$REMOTE_DATA"
require_safe_path "REMOTE_INPUT" "$REMOTE_INPUT"
require_safe_path "REMOTE_OUT" "$REMOTE_OUT"
require_safe_path "LOCAL_OUT" "$LOCAL_OUT"
require_safe_path "LOCAL_SUMMARY" "$LOCAL_SUMMARY"
require_safe_path "LOCAL_MOLDS" "$LOCAL_MOLDS"
require_safe_path "LOCAL_CONFIG" "$LOCAL_CONFIG"
require_safe_path "LOCAL_FILE_LIST" "$LOCAL_FILE_LIST"
require_safe_path "STMFIT_LOCAL_DATA" "${STMFIT_LOCAL_DATA:-}"
require_safe_path "LOCAL_CC_PROVENANCE" "$LOCAL_CC_PROVENANCE"
require_safe_path "LOCAL_CC_MAPS" "$LOCAL_CC_MAPS"
require_safe_path "LOCAL_CC_MASK" "$LOCAL_CC_MASK"
require_safe_path "LOCAL_CC_MOLD_PROVENANCE" "$LOCAL_CC_MOLD_PROVENANCE"
require_safe_path "REMOTE_MOLDS" "$REMOTE_MOLDS"
require_safe_path "REMOTE_SUMMARY" "$REMOTE_SUMMARY"
require_safe_path "REMOTE_FILE_LIST" "$REMOTE_FILE_LIST"
require_safe_path "REMOTE_CONFIG" "$REMOTE_CONFIG"
require_safe_path "REMOTE_CC_PROVENANCE" "$REMOTE_CC_PROVENANCE"
require_safe_path "REMOTE_CC_MAPS" "$REMOTE_CC_MAPS"
require_safe_path "REMOTE_CC_MASK" "$REMOTE_CC_MASK"
require_safe_path "REMOTE_CC_MOLD_PROVENANCE" "$REMOTE_CC_MOLD_PROVENANCE"
require_safe_token "JULIA_MODULE_VERSION" "$JULIA_MODULE_VERSION"
require_safe_token "WALLTIME" "$WALLTIME"
require_safe_token "CC_PROVENANCE_SHA" "$CC_PROVENANCE_SHA"
for artifact in "${LOCAL_CC_EXTRA_ARTIFACTS[@]}"; do
    require_safe_path "constant-current bracket artifact" "$artifact"
done

require_file_for_submit() {
    local path=$1
    local label=$2
    if [[ ! -f "$path" ]]; then
        if (( DRY_RUN )); then
            echo "warning: dry-run did not find $label: $path" >&2
        else
            echo "missing $label: $path" >&2
            exit 1
        fi
    fi
}

require_data_for_submit() {
    if [[ -z "$STMFIT_LOCAL_DATA" || ! -d "$STMFIT_LOCAL_DATA" ]]; then
        if (( DRY_RUN )); then
            echo "warning: dry-run did not find local SXM directory: ${STMFIT_LOCAL_DATA:-<unset>}" >&2
            return
        fi
        echo "missing local SXM directory: ${STMFIT_LOCAL_DATA:-<unset>}" >&2
        exit 1
    fi
    local file
    for file in "${FILES[@]}"; do
        if [[ ! -f "$STMFIT_LOCAL_DATA/$file" ]]; then
            if (( DRY_RUN )); then
                echo "warning: dry-run did not find $STMFIT_LOCAL_DATA/$file" >&2
            else
                echo "missing $STMFIT_LOCAL_DATA/$file" >&2
                exit 1
            fi
        fi
    done
}

require_file_for_submit "$LOCAL_SUMMARY" "selection summary"
require_file_for_submit "$LOCAL_MOLDS" "converged mold TSV"
if [[ "$LOCAL_CONFIG" = /* ]]; then
    require_file_for_submit "$LOCAL_CONFIG" "diagnostic config"
else
    require_file_for_submit "$REPO_ROOT/$LOCAL_CONFIG" "diagnostic config"
fi
require_file_for_submit "$LOCAL_FILE_LIST" "SXM file list"
if [[ "$(basename "$LOCAL_CONFIG")" == "joint_proxy_whole_roi_constant_current.toml" ]]; then
    require_file_for_submit "$LOCAL_CC_PROVENANCE" "constant-current provenance"
    require_file_for_submit "$LOCAL_CC_MAPS" "constant-current map"
    require_file_for_submit "$LOCAL_CC_MASK" "constant-current validity mask"
    require_file_for_submit "$LOCAL_CC_MOLD_PROVENANCE" "constant-current mold binding provenance"
    if [[ -z "$CC_PROVENANCE_SHA" ]]; then
        if (( DRY_RUN )); then
            echo "warning: dry-run did not receive STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256" >&2
        else
            echo "missing STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256" >&2
            exit 1
        fi
    fi
fi
require_data_for_submit

ssh_cmd=(ssh -o BatchMode=no -o ConnectTimeout="$SSH_CONNECT_TIMEOUT"
    -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" "$STMFIT_SSH_HOST")
run() {
    if (( DRY_RUN )); then printf '  [dry-run]'; printf ' %q' "$@"; printf '\n'
    else "$@"; fi
}
remote() {
    if (( DRY_RUN )); then echo "  [dry-run] ssh $STMFIT_SSH_HOST $1"
    else "${ssh_cmd[@]}" "$1"; fi
}

echo "Frozen-contrast HPC launch"
echo "  target:  $STMFIT_SSH_HOST"
echo "  config:  $LOCAL_CONFIG"
echo "  molds:   $LOCAL_MOLDS"
echo "  summary: $LOCAL_SUMMARY"
echo "  file-list: $LOCAL_FILE_LIST"
echo "  files:   ${FILES[*]}"
echo "  output:  $REMOTE_OUT -> $LOCAL_OUT"
echo "  array size: $array_size"
echo "  max concurrency: $max_concurrency"
echo "  cpus per task: $CPUS_PER_TASK"
echo "  julia threads/task: $JULIA_THREADS_PER_TASK"
echo "  memory/task: $((CPUS_PER_TASK * MEM_PER_CPU)) MB"
echo "  walltime: $WALLTIME"
[[ -n "$LOCAL_CC_PROVENANCE" ]] && echo "  observable provenance: $LOCAL_CC_PROVENANCE"
[[ -n "$LOCAL_CC_MAPS" ]] && echo "  observable map: $LOCAL_CC_MAPS"
[[ -n "$LOCAL_CC_MASK" ]] && echo "  validity mask: $LOCAL_CC_MASK"
[[ -n "$LOCAL_CC_MOLD_PROVENANCE" ]] && echo "  mold provenance: $LOCAL_CC_MOLD_PROVENANCE"
[[ -n "$CC_PROVENANCE_SHA" ]] && echo "  observable provenance sha256: $CC_PROVENANCE_SHA"

remote "mkdir -p $(shell_quote "$REMOTE_PROJECT") $(shell_quote "$REMOTE_DATA") $(shell_quote "$REMOTE_INPUT") $(shell_quote "$REMOTE_OUT/logs")"
run rsync -avz --exclude=results/ --exclude=.git/ --exclude=Manifest.toml \
    --exclude=hpc/remote.env --exclude=.omo/ --exclude=.slim/ \
    --exclude=docs/build/ --exclude=qe/ "$REPO_ROOT/" "$STMFIT_SSH_HOST:$REMOTE_PROJECT/"
for file in "${FILES[@]}"; do
    run rsync -avz --update "${STMFIT_LOCAL_DATA:-/missing-local-sxm}/$file" \
        "$STMFIT_SSH_HOST:$REMOTE_DATA/$file"
done
run rsync -avz "$LOCAL_SUMMARY" "$STMFIT_SSH_HOST:$REMOTE_SUMMARY"
run rsync -avz "$LOCAL_FILE_LIST" "$STMFIT_SSH_HOST:$REMOTE_FILE_LIST"
if [[ -f "$LOCAL_MOLDS" || ! -e "$LOCAL_MOLDS" ]]; then
    run rsync -avz "$LOCAL_MOLDS" "$STMFIT_SSH_HOST:$REMOTE_MOLDS"
fi
if [[ "$LOCAL_CONFIG" = /* ]]; then
    run rsync -avz "$LOCAL_CONFIG" "$STMFIT_SSH_HOST:$REMOTE_CONFIG"
fi
[[ -n "$LOCAL_CC_PROVENANCE" ]] && run rsync -avz "$LOCAL_CC_PROVENANCE" "$STMFIT_SSH_HOST:$REMOTE_CC_PROVENANCE"
[[ -n "$LOCAL_CC_MAPS" ]] && run rsync -avz "$LOCAL_CC_MAPS" "$STMFIT_SSH_HOST:$REMOTE_CC_MAPS"
[[ -n "$LOCAL_CC_MASK" ]] && run rsync -avz "$LOCAL_CC_MASK" "$STMFIT_SSH_HOST:$REMOTE_CC_MASK"
[[ -n "$LOCAL_CC_MOLD_PROVENANCE" ]] && run rsync -avz "$LOCAL_CC_MOLD_PROVENANCE" "$STMFIT_SSH_HOST:$REMOTE_CC_MOLD_PROVENANCE"
for artifact in "${LOCAL_CC_EXTRA_ARTIFACTS[@]}"; do
    run rsync -avz "$artifact" "$STMFIT_SSH_HOST:$(remote_stage_path "$artifact")"
done

remote "test -s $(shell_quote "$REMOTE_MOLDS") || { echo $(shell_quote "missing converged mold: $REMOTE_MOLDS") >&2; exit 3; }"
if [[ "$(basename "$LOCAL_CONFIG")" == "joint_proxy_whole_roi_constant_current.toml" ]]; then
    if [[ -n "$REMOTE_CC_PROVENANCE" && -n "$REMOTE_CC_MAPS" && -n "$REMOTE_CC_MASK" && -n "$REMOTE_CC_MOLD_PROVENANCE" ]]; then
        remote "test -s $(shell_quote "$REMOTE_CC_PROVENANCE") && test -s $(shell_quote "$REMOTE_CC_MAPS") && test -s $(shell_quote "$REMOTE_CC_MASK") && test -s $(shell_quote "$REMOTE_CC_MOLD_PROVENANCE") || { echo 'missing constant-current observable artifacts' >&2; exit 3; }"
    elif (( DRY_RUN )); then
        echo "  [dry-run] constant-current observable artifact check deferred: provenance/map/mask/mold-provenance paths are unset"
    else
        echo "missing constant-current provenance/map/mask/mold-provenance inputs" >&2
        exit 1
    fi
fi
remote "cd $(shell_quote "$REMOTE_PROJECT") && module purge && module load $(shell_quote "julia${JULIA_MODULE_VERSION:+/$JULIA_MODULE_VERSION}") && julia --project=. -e 'using Pkg; Pkg.instantiate()'"

mem_mb=$((CPUS_PER_TASK * MEM_PER_CPU))
exports="ALL,STMFIT_REMOTE_PROJECT=$REMOTE_PROJECT,STMFIT_DATA_DIR=$REMOTE_DATA,STMFIT_WHOLE_ROI_MOLDS=$REMOTE_MOLDS,STMFIT_WHOLE_ROI_SUMMARY=$REMOTE_SUMMARY,STMFIT_WHOLE_ROI_CONFIG=$REMOTE_CONFIG,STMFIT_WHOLE_ROI_FILE_LIST=$REMOTE_FILE_LIST,STMFIT_WHOLE_ROI_OUTDIR=$REMOTE_OUT,STMFIT_JULIA_THREADS_PER_TASK=$JULIA_THREADS_PER_TASK,JULIA_MODULE_VERSION=$JULIA_MODULE_VERSION,STMFIT_CONSTANT_CURRENT_PROVENANCE=$REMOTE_CC_PROVENANCE,STMFIT_CONSTANT_CURRENT_PROVENANCE_SHA256=$CC_PROVENANCE_SHA,STMFIT_CONSTANT_CURRENT_MAPS=$REMOTE_CC_MAPS,STMFIT_CONSTANT_CURRENT_INVALID_MASK=$REMOTE_CC_MASK,STMFIT_CONSTANT_CURRENT_MOLD_PROVENANCE=$REMOTE_CC_MOLD_PROVENANCE"
submit="cd $(shell_quote "$REMOTE_PROJECT") && sbatch --parsable --array=1-$array_size%$max_concurrency --cpus-per-task=$CPUS_PER_TASK --mem=${mem_mb}MB --time=$WALLTIME --export=$(shell_quote "$exports") hpc/whole_roi_frozen_contrast_array.sbatch"
if (( DRY_RUN )); then
    remote "$submit"
    echo "Dry-run complete; no SSH connection or Slurm submission was performed. This records launch plumbing only."
    exit 0
fi

jobid=$("${ssh_cmd[@]}" "$submit")
jobid="${jobid%%;*}"
echo "submitted job $jobid"
if (( WATCH )); then
    while "${ssh_cmd[@]}" "squeue -h -j $(shell_quote "$jobid")" | grep -q .; do sleep 30; done
    mkdir -p "$LOCAL_OUT"
    rsync -avz "$STMFIT_SSH_HOST:$REMOTE_OUT/" "$LOCAL_OUT/"
    echo "fetched results to $LOCAL_OUT"
fi

#!/usr/bin/env bash
# Sync and submit the Quantum ESPRESSO mold jobs on MPCDF.
#
# This is separate from launch_remote.sh, which submits the STMFit image batch.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_FILE="$SCRIPT_DIR/remote.env"

if [[ ! -f "$ENV_FILE" ]]; then
    echo "ERROR: $ENV_FILE not found." >&2
    echo "       cp hpc/remote.env.example hpc/remote.env   then edit it." >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$ENV_FILE"

: "${STMFIT_SSH_HOST:?STMFIT_SSH_HOST not set in $ENV_FILE}"
: "${JULIA_MODULE_VERSION:=}"
: "${WATCH_POLL:=60}"
: "${RSYNC_EXTRA:=}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"
: "${QE_COMPILER_MODULE:=intel/2024.0}"
: "${QE_MPI_MODULE:=impi/2021.11}"
: "${QE_MODULE:=qe/7.4.1}"
: "${QE_QOS:=}"
: "${QE_MIN_MEM_MB:=0}"
: "${QE_DEPENDENCY:=}"

if [[ -z "${STMFIT_REMOTE_USER:-}" ]]; then
    echo "ERROR: STMFIT_REMOTE_USER is not set in $ENV_FILE." >&2
    exit 1
fi

: "${REMOTE_PROJECT_BASE:=/u/$STMFIT_REMOTE_USER/code}"
STMFIT_REMOTE_PROJECT="$REMOTE_PROJECT_BASE/STMFit"

DRY_RUN=0
WATCH=0
SYNC_CODE=1
SYNC_QE=1
MAX_TOTAL_TASKS=8
SEQUENTIAL=${QE_SEQUENTIAL:-1}
DIRS=("qe/glcn" "qe/glcnac")

usage() {
    cat <<'EOF'
Usage: hpc/launch_qe_molds_remote.sh [options]

Options:
  --dry-run             Show sync/preflight/submit commands without changing remote state
  --watch               Pass --watch to hpc/submit_qe_molds.sh on the cluster
  --no-sync-code        Skip repository code sync
  --no-sync-qe          Skip qe/ input-directory sync
  --max-total-tasks N   Max total Slurm tasks across QE jobs [8]
  --min-mem-mb N        Optional minimum #SBATCH --mem per job [$QE_MIN_MEM_MB]
  --dependency SPEC     Slurm dependency for submitted jobs, e.g. afterany:12345
  --sequential          Submit QE run dirs as afterok dependency chain [default]
  --parallel            Submit QE run dirs independently
  --dir PATH            QE run directory relative to repo. Repeatable. Default qe/glcn qe/glcnac
  -h, --help            Show this help

Requires hpc/remote.env. Run locally from the STMFit checkout after preparing
QE input directories. Full runs provide run_qe_mold.sbatch; preliminary SCF+PP
runs may provide run_scf_pp.sbatch.
EOF
}

custom_dirs=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --watch) WATCH=1; shift ;;
        --no-sync-code) SYNC_CODE=0; shift ;;
        --no-sync-qe) SYNC_QE=0; shift ;;
        --sequential) SEQUENTIAL=1; shift ;;
        --parallel) SEQUENTIAL=0; shift ;;
        --max-total-tasks) MAX_TOTAL_TASKS="$2"; shift 2 ;;
        --max-total-tasks=*) MAX_TOTAL_TASKS="${1#*=}"; shift ;;
        --min-mem-mb) QE_MIN_MEM_MB="$2"; shift 2 ;;
        --min-mem-mb=*) QE_MIN_MEM_MB="${1#*=}"; shift ;;
        --dependency) QE_DEPENDENCY="$2"; shift 2 ;;
        --dependency=*) QE_DEPENDENCY="${1#*=}"; shift ;;
        --dir)
            if (( custom_dirs == 0 )); then DIRS=(); custom_dirs=1; fi
            DIRS+=("$2"); shift 2 ;;
        --dir=*)
            if (( custom_dirs == 0 )); then DIRS=(); custom_dirs=1; fi
            DIRS+=("${1#*=}"); shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

module_load_cmd="module purge && module load julia${JULIA_MODULE_VERSION:+/$JULIA_MODULE_VERSION}"

run_local() {
    if (( DRY_RUN )); then
        echo "  [dry-run] $*"
    else
        "$@"
    fi
}

run_remote() {
    local cmd="$1"
    if (( DRY_RUN )); then
        echo "  [dry-run] ssh $STMFIT_SSH_HOST \"$cmd\""
    else
        ssh -o BatchMode=no \
            -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
            -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
            "$STMFIT_SSH_HOST" "$cmd"
    fi
}

echo "QE mold remote launch"
echo "  host:    $STMFIT_SSH_HOST"
echo "  project: $STMFIT_REMOTE_PROJECT"
echo "  dirs:    ${DIRS[*]}"
echo "  tasks:   max $MAX_TOTAL_TASKS"
echo "  mem:     min ${QE_MIN_MEM_MB}MB"
echo "  submit:  $([[ $SEQUENTIAL == 1 ]] && echo sequential || echo parallel)"
echo "  qe:      $QE_COMPILER_MODULE $QE_MPI_MODULE $QE_MODULE"
echo "  qos:     ${QE_QOS:-default}"
echo "  dep:     ${QE_DEPENDENCY:-none}"
echo "  ssh:     ConnectTimeout=$SSH_CONNECT_TIMEOUT ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL"
(( DRY_RUN )) && echo "  mode:    dry-run"

if (( ! DRY_RUN )); then
    ssh -o BatchMode=no \
        -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
        -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
        "$STMFIT_SSH_HOST" 'echo ok' >/dev/null \
        || { echo "ERROR: SSH to $STMFIT_SSH_HOST failed." >&2; exit 1; }
fi

if (( SYNC_CODE )); then
    echo "Syncing code to $STMFIT_REMOTE_PROJECT"
    run_remote "mkdir -p '$STMFIT_REMOTE_PROJECT'"
    run_local rsync -avz --delete \
        --exclude='.git/' \
        --exclude='results/' \
        --exclude='*.log' \
        --exclude='Manifest.toml' \
        --exclude='hpc/remote.env' \
        --exclude='docs/build/' \
        --exclude='qe/' \
        ${RSYNC_EXTRA} \
        -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
        "$REPO_ROOT/" "$STMFIT_SSH_HOST:$STMFIT_REMOTE_PROJECT/"
fi

if (( SYNC_QE )); then
    echo "Syncing prepared QE input directories"
    for dir in "${DIRS[@]}"; do
        [[ -d "$REPO_ROOT/$dir" ]] || { echo "ERROR: local QE dir not found: $dir" >&2; exit 1; }
        run_remote "mkdir -p '$STMFIT_REMOTE_PROJECT/$dir'"
        run_local rsync -avz --delete \
            --include='pw_relax.in' \
            --include='pw_scf.in' \
            --include='pp_ldos.in' \
            --include='run_qe_mold.sbatch' \
            --include='run_scf_pp.sbatch' \
            --include='pseudo/' \
            --include='pseudo/*.UPF' \
            --exclude='*' \
            ${RSYNC_EXTRA} \
            -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
            "$REPO_ROOT/$dir/" "$STMFIT_SSH_HOST:$STMFIT_REMOTE_PROJECT/$dir/"
    done
fi

preflight_args=()
submit_args=(--max-total-tasks "$MAX_TOTAL_TASKS")
[[ "$QE_MIN_MEM_MB" != "0" ]] && submit_args+=(--min-mem-mb "$QE_MIN_MEM_MB")
[[ -n "$QE_DEPENDENCY" ]] && submit_args+=(--dependency "$QE_DEPENDENCY")
for dir in "${DIRS[@]}"; do
    preflight_args+=(--dir "$dir")
    submit_args+=(--dir "$dir")
done
preflight_args+=(--out hpc/qe_molds/qe_input_preflight.tsv --max-total-tasks "$MAX_TOTAL_TASKS")
[[ "$QE_MIN_MEM_MB" != "0" ]] && preflight_args+=(--min-mem-mb "$QE_MIN_MEM_MB")
(( SEQUENTIAL )) && preflight_args+=(--sequential)
(( WATCH )) && submit_args+=(--watch --watch-poll "$WATCH_POLL")
(( SEQUENTIAL )) && submit_args+=(--sequential)

remote_preflight="cd '$STMFIT_REMOTE_PROJECT' && $module_load_cmd && julia --project=. test/preflight_qe_mold_inputs.jl"
for arg in "${preflight_args[@]}"; do
    remote_preflight+=" '$arg'"
done

echo "Running remote preflight"
run_remote "$remote_preflight"

remote_submit="cd '$STMFIT_REMOTE_PROJECT' && JULIA_MODULE_VERSION='$JULIA_MODULE_VERSION' QE_COMPILER_MODULE='$QE_COMPILER_MODULE' QE_MPI_MODULE='$QE_MPI_MODULE' QE_MODULE='$QE_MODULE' QE_QOS='$QE_QOS' QE_MIN_MEM_MB='$QE_MIN_MEM_MB' QE_DEPENDENCY='$QE_DEPENDENCY' bash hpc/submit_qe_molds.sh"
for arg in "${submit_args[@]}"; do
    remote_submit+=" '$arg'"
done

echo "Submitting remote QE jobs"
run_remote "$remote_submit"

echo "Done. Remote reports:"
echo "  $STMFIT_REMOTE_PROJECT/hpc/qe_molds/qe_input_preflight.tsv"
echo "  $STMFIT_REMOTE_PROJECT/hpc/qe_molds/qe_jobs.tsv"

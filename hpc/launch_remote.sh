#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# launch_remote.sh — push-button STMFit run on the MPCDF cluster (Raven/Viper).
#
# Workflow: sync code + data → instantiate Julia → submit Slurm array →
#           (optionally) watch → merge chunk summaries → fetch results back.
#
# Configure first:   cp hpc/remote.env.example hpc/remote.env   (then edit)
# Dry run:           ./hpc/launch_remote.sh --dry-run
# Full run:          ./hpc/launch_remote.sh
# Watch + fetch:     ./hpc/launch_remote.sh --watch
# Just fetch an already-finished job:  ./hpc/launch_remote.sh --fetch-only
#
# All remote commands go through your ~/.ssh/config host alias (ProxyJump +
# ControlMaster recommended so you only authenticate once). See hpc/README.md.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ── Locate ourselves & the repo root ─────────────────────────────────────────
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

# ── Defaults & arg parsing ──────────────────────────────────────────────────
: "${STMFIT_SSH_HOST:?STMFIT_SSH_HOST not set in $ENV_FILE}"
: "${N_CHUNKS:=4}"
: "${CPUS_PER_TASK:=4}"
: "${WALLTIME:=04:00:00}"
: "${MEM_PER_CPU:=4000}"
: "${JULIA_MODULE_VERSION:=}"
: "${STMFIT_CONFIG:=config/chitosan.toml}"
: "${STMFIT_OUTDIR:=results/best_plots}"
: "${STMFIT_BATCH_ARGS:=}"
: "${N_FILES:=}"
: "${STMFIT_MAIL_USER:=}"
: "${WATCH_POLL:=30}"
: "${RSYNC_EXTRA:=}"
: "${SSH_CONNECT_TIMEOUT:=180}"
: "${SSH_SERVER_ALIVE_INTERVAL:=60}"

# ── Remote username + paths (MPCDF convention) ───────────────────────────────
# STMFIT_REMOTE_USER is the short MPCDF account name used in filesystem paths
# (/u/<user>/..., /ptmp/<user>/...). NOTE: this is NOT necessarily the same as
# the SSH `User` in ~/.ssh/config (MPCDF often accepts your full email for SSH
# 2FA, but filesystem paths always use the short account name). Set it in
# hpc/remote.env — it's the only personal value the launcher needs.
if [[ -z "${STMFIT_REMOTE_USER:-}" ]]; then
    printf "\033[31m✗ STMFIT_REMOTE_USER is not set.\033[0m\n" >&2
    echo "  Set it in $ENV_FILE to your short MPCDF account name." >&2
    echo "  This is the name used in /u/<user>/... and /ptmp/<user>/... paths." >&2
    exit 1
fi

# Remote paths follow the MPCDF convention. Override via REMOTE_PROJECT_BASE /
# REMOTE_DATA_BASE / STMFIT_LOCAL_DATA only if your layout differs (see env file).
: "${REMOTE_PROJECT_BASE:=/u/$STMFIT_REMOTE_USER/code}"
: "${REMOTE_DATA_BASE:=/ptmp/$STMFIT_REMOTE_USER/stmfit}"
STMFIT_REMOTE_PROJECT="$REMOTE_PROJECT_BASE/STMFit"
STMFIT_REMOTE_DATA="$REMOTE_DATA_BASE/data"
# Local .sxm source: prefer $STMFIT_DATA_DIR (project convention), else error
# out — never hardcode a personal path. Set STMFIT_LOCAL_DATA in remote.env.
: "${STMFIT_LOCAL_DATA:=${STMFIT_DATA_DIR:-}}"
export STMFIT_REMOTE_PROJECT STMFIT_REMOTE_DATA STMFIT_LOCAL_DATA

DRY_RUN=0
SYNC_CODE=1
SYNC_DATA=1
WATCH=0          # default: submit and exit (don't block). Use --watch to wait.
INstantiate=1
MERGE_REMOTE=1
FETCH_ONLY=0

usage() {
    cat <<'EOF'
Usage: launch_remote.sh [options]
  --dry-run         Show what would happen; don't sync/submit/fetch.
  --no-sync-code    Skip rsync of the repo to the cluster.
  --no-sync-data    Skip rsync of .sxm data.
  --no-instantiate  Skip remote `Pkg.instantiate()` (assume already done).
  --no-merge        Skip remote chunk-summary merge after the array.
  --watch           Block until the array finishes, then merge + fetch.
  --fetch-only      Only fetch results from a job already submitted/finished.
  -h, --help        Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)        DRY_RUN=1; shift ;;
        --no-sync-code)   SYNC_CODE=0; shift ;;
        --no-sync-data)   SYNC_DATA=0; shift ;;
        --no-instantiate) INstantiate=0; shift ;;
        --no-merge)       MERGE_REMOTE=0; shift ;;
        --watch)          WATCH=1; shift ;;
        --fetch-only)     FETCH_ONLY=1; SYNC_CODE=0; SYNC_DATA=0; INstantiate=0; WATCH=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# ── Helpers ─────────────────────────────────────────────────────────────────
C_BOLD=$'\033[1m'; C_CYAN=$'\033[36m'; C_GREEN=$'\033[32m'; C_YELLOW=$'\033[33m'; C_RED=$'\033[31m'; C_RESET=$'\033[0m'
say()  { printf "%s▸ %s%s%s\n" "$C_CYAN" "$C_BOLD" "$*" "$C_RESET"; }
ok()   { printf "%s✓ %s%s\n" "$C_GREEN" "$*" "$C_RESET"; }
warn() { printf "%s⚠ %s%s\n" "$C_YELLOW" "$*" "$C_RESET" >&2; }
die()  { printf "%s✗ %s%s\n" "$C_RED" "$*" "$C_RESET" >&2; exit 1; }

runc() { # runc <label> <cmd...>  — runs locally, echoes under --dry-run
    if (( DRY_RUN )); then echo "  [dry-run] $*"; return; fi
    "$@"
}
runssh() { # runssh <cmd>  — runs on the cluster; echoes under --dry-run
    local cmd="$1"
    if (( DRY_RUN )); then echo "  [dry-run] ssh $STMFIT_SSH_HOST \"$cmd\""; return; fi
    # shellcheck disable=SC2086
    ssh -o BatchMode=no \
        -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
        -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
        "$STMFIT_SSH_HOST" "$cmd"
}

# ── 0. Pre-flight ───────────────────────────────────────────────────────────
say "STMFit remote launch — target: $STMFIT_SSH_HOST"
[[ -n "$N_FILES" ]] && echo "    files: $N_FILES   chunks: $N_CHUNKS   cpus/task: $CPUS_PER_TASK   time: $WALLTIME" \
                    || echo "    chunks: $N_CHUNKS   cpus/task: $CPUS_PER_TASK   time: $WALLTIME"
echo "    project: $STMFIT_REMOTE_PROJECT"
echo "    data:    $STMFIT_LOCAL_DATA  →  $STMFIT_REMOTE_DATA"
echo "    config:  $STMFIT_CONFIG   outdir: $STMFIT_OUTDIR"
echo "    ssh:     ConnectTimeout=$SSH_CONNECT_TIMEOUT ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL"
(( DRY_RUN )) && warn "DRY RUN — no changes will be made."

# A quick connectivity probe (skipped in dry-run to avoid prompting).
if (( ! DRY_RUN )); then
    ssh -o BatchMode=no \
        -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
        -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
        "$STMFIT_SSH_HOST" 'echo ok' >/dev/null 2>&1 \
        || die "SSH to '$STMFIT_SSH_HOST' failed. Check ~/.ssh/config (ProxyJump + 2FA)."
    ok "SSH reachable."
fi

# ── 1. Sync code ────────────────────────────────────────────────────────────
if (( SYNC_CODE )); then
    say "Syncing code → $STMFIT_REMOTE_PROJECT"
    runc rsync -avz --delete \
        --exclude='results/' --exclude='*.log' --exclude='.git/' \
        --exclude='hpc/remote.env' --exclude='Manifest.toml' \
        ${RSYNC_EXTRA} \
        -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
        "$REPO_ROOT/" "$STMFIT_SSH_HOST:$STMFIT_REMOTE_PROJECT/"
    ok "Code synced."
fi

# ── 2. Sync data ────────────────────────────────────────────────────────────
if (( SYNC_DATA )); then
    if [[ -z "$STMFIT_LOCAL_DATA" ]]; then
        warn "STMFIT_LOCAL_DATA not set — skipping data sync."
        warn "  Set it in $ENV_FILE (or export STMFIT_DATA_DIR) to push .sxm files."
    elif [[ ! -d "$STMFIT_LOCAL_DATA" ]]; then
        warn "Local data dir not found: $STMFIT_LOCAL_DATA — skipping data sync."
    else
        say "Syncing data → $STMFIT_REMOTE_DATA"
        runssh "mkdir -p '$STMFIT_REMOTE_DATA'"
        runc rsync -avz --update --include='*/' --include='*.sxm' --exclude='*' \
            ${RSYNC_EXTRA} \
            -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
            "$STMFIT_LOCAL_DATA/" "$STMFIT_SSH_HOST:$STMFIT_REMOTE_DATA/"
        ok "Data synced."
    fi
fi

# ── 3. Julia instantiate (on the login node — has internet) ─────────────────
if (( INstantiate )); then
    say "Instantiating Julia project on login node (downloads packages if needed)"
    local_julia_ver="${JULIA_MODULE_VERSION:-default}"
    runssh "cd '$STMFIT_REMOTE_PROJECT' && \
        module purge && module load julia/${JULIA_MODULE_VERSION} && \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); @info \"instantiate ok\"'"
    ok "Julia project ready."
fi

# ── fetch-only short-circuits here ──────────────────────────────────────────
if (( FETCH_ONLY )); then
    say "Fetch-only: skipping submission."
    WATCH=1
    JOBID=""
else
    # ── 4. Submit Slurm array ───────────────────────────────────────────────
    say "Submitting Slurm array ($N_CHUNKS chunks)..."
    # Build the export list. N_CHUNKS is passed so --chunk i/N_CHUNKS stays
    # consistent across partial re-submissions. STMFIT_* mirror the sbatch env.
    export_vars="ALL,N_CHUNKS=$N_CHUNKS"
    [[ -n "$STMFIT_CONFIG" ]]        && export_vars+=",STMFIT_CONFIG=$STMFIT_CONFIG"
    [[ -n "$STMFIT_OUTDIR" ]]        && export_vars+=",STMFIT_OUTDIR=$STMFIT_OUTDIR"
    [[ -n "$STMFIT_BATCH_ARGS" ]]    && export_vars+=",STMFIT_BATCH_ARGS=$STMFIT_BATCH_ARGS"
    [[ -n "$N_FILES" ]]              && export_vars+=",N_FILES=$N_FILES"
    [[ -n "$JULIA_MODULE_VERSION" ]] && export_vars+=",JULIA_MODULE_VERSION=$JULIA_MODULE_VERSION"
    # NOTE: STMFIT_MAIL_USER is passed via --mail-user below (the only place
    # Slurm reads it), so it's not added to --export.

    # Total job memory = MEM_PER_CPU × CPUS_PER_TASK, expressed in MB. Raven's
    # shared-job submit filter REQUIRES a total `--mem` (with explicit MB unit);
    # `--mem-per-cpu` is rejected ("memory limit must be provided for shared
    # jobs"). Viper accepts either, so --mem is the portable choice.
    MEM_TOTAL_MB=$(( MEM_PER_CPU * CPUS_PER_TASK ))

    sbatch_cmd=(
        sbatch
        --export="$export_vars"
        --array="1-$N_CHUNKS"
        --cpus-per-task="$CPUS_PER_TASK"
        --time="$WALLTIME"
        --mem="${MEM_TOTAL_MB}MB"
    )
    [[ -n "$STMFIT_MAIL_USER" ]] && sbatch_cmd+=(--mail-user="$STMFIT_MAIL_USER")
    sbatch_cmd+=("$STMFIT_REMOTE_PROJECT/hpc/batch_array.sbatch")

    if (( DRY_RUN )); then
        echo "  [dry-run] ssh $STMFIT_SSH_HOST \"cd $STMFIT_REMOTE_PROJECT && ${sbatch_cmd[*]}\""
        JOBID="<dry-run>"
    else
        # shellcheck disable=SC2086
        submit_out=$(ssh -o BatchMode=no \
            -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
            -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
            "$STMFIT_SSH_HOST" \
            "cd '$STMFIT_REMOTE_PROJECT' && ${sbatch_cmd[*]}") \
            || die "sbatch submission failed: $submit_out"
        JOBID=$(echo "$submit_out" | grep -oE '[0-9]+' | head -1)
        [[ -n "$JOBID" ]] || die "Could not parse job id from: $submit_out"
        ok "Submitted job array $JOBID."
        echo "    Monitor:  ssh $STMFIT_SSH_HOST 'squeue -j $JOBID'"
        echo "    Cancel:   ssh $STMFIT_SSH_HOST 'scancel $JOBID'"
    fi
fi

# ── 5. Watch (optional) ────────────────────────────────────────────────────
if (( WATCH )); then
    if [[ -z "${JOBID:-}" ]]; then
        warn "No job id to watch (--fetch-only or dry-run). Enter the array job id to watch,"
        warn "or blank to skip: " ; read -r JOBID
        [[ -z "$JOBID" ]] && { warn "No job id; finishing."; exit 0; }
    fi
    if (( DRY_RUN )); then
        echo "  [dry-run] would poll squeue -j $JOBID every ${WATCH_POLL}s until the array ends"
        exit 0
    fi
    say "Watching array $JOBID (poll every ${WATCH_POLL}s)..."
    # Poll until no task of this job is queued/running.
    while :; do
        n=$(ssh -o BatchMode=no \
            -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
            -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
            "$STMFIT_SSH_HOST" \
            "squeue -h -j '$JOBID' -t R,PD,CF,CA,S,ST 2>/dev/null | wc -l" || echo 0)
        if [[ "$n" =~ ^[0-9]+$ ]] && (( n == 0 )); then break; fi
        printf "  [%s] %s tasks still queued/running...\r" "$(date +%H:%M:%S)" "$n"
        sleep "$WATCH_POLL"
    done
    echo
    ok "Array $JOBID finished."
    ssh -o BatchMode=no \
        -o ConnectTimeout="$SSH_CONNECT_TIMEOUT" \
        -o ServerAliveInterval="$SSH_SERVER_ALIVE_INTERVAL" \
        "$STMFIT_SSH_HOST" "sacct -j '$JOBID' --format=JobID,ArrayTaskID,State,Elapsed,ExitCode,MaxRSS -n" || true

    # ── 6. Merge chunk summaries on the login node ─────────────────────────
    if (( MERGE_REMOTE )); then
        say "Merging chunk summaries on the cluster..."
        runssh "cd '$STMFIT_REMOTE_PROJECT' && \
            module purge && module load julia/${JULIA_MODULE_VERSION} && \
            julia --project=. hpc/merge_chunks.jl '$STMFIT_OUTDIR' --total '$N_CHUNKS'" \
            || warn "Remote merge reported issues (missing chunks?). See output above."
        ok "Merge step done."
    fi
fi

# ── 7. Fetch results ────────────────────────────────────────────────────────
say "Fetching results ← $STMFIT_REMOTE_PROJECT/$STMFIT_OUTDIR"
LOCAL_OUT="$REPO_ROOT/$STMFIT_OUTDIR"
runc rsync -avz --update \
    --include='*/' --include='*.png' --include='*.tsv' --include='*.txt' --exclude='*' \
    ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_SSH_HOST:$STMFIT_REMOTE_PROJECT/$STMFIT_OUTDIR/" "$LOCAL_OUT/"
ok "Results fetched to $LOCAL_OUT"

# Also pull the per-chunk Slurm logs for inspection.
runc rsync -avz --update \
    ${RSYNC_EXTRA} \
    -e "ssh -o ConnectTimeout=$SSH_CONNECT_TIMEOUT -o ServerAliveInterval=$SSH_SERVER_ALIVE_INTERVAL" \
    "$STMFIT_SSH_HOST:$STMFIT_REMOTE_PROJECT/results/hpc_logs/" "$REPO_ROOT/results/hpc_logs/" \
    || warn "No hpc_logs to fetch (ok if logs dir was empty)."

(( DRY_RUN )) && warn "DRY RUN complete — nothing was actually run." \
              || ok "Done. Summary: $LOCAL_OUT/summary_overlap060_hard.tsv"

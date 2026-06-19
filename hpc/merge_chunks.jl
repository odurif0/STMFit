#!/usr/bin/env julia
# Merge the per-chunk summary TSVs produced by `test/batch_full.jl --chunk i/n`
# back into a single summary file, mirroring the non-chunked output name.
#
# batch_full.jl writes one of:
#   summary_overlap060_hard.tsv                 (when CHUNK_TOTAL == 1)
#   summary_overlap060_hard_chunkNNofMM.tsv     (one per chunk, when sharding)
# This script concatenates the chunkNNofMM shards into a single TSV with the
# header written exactly once and data rows sorted by filepath. It reuses the
# read/concat pattern from test/summarize.jl (lines 42-51).
#
# Usage (run from the repo root, on a login node, after the array finishes):
#   julia --project=. hpc/merge_chunks.jl <outdir> [options]
#
# Examples:
#   julia --project=. hpc/merge_chunks.jl results/best_plots
#   julia --project=. hpc/merge_chunks.jl results/best_plots --total 4
#   julia --project=. hpc/merge_chunks.jl results/best_plots \
#       --summary-name summary_overlap060_hard.tsv --chunk-prefix summary_overlap060_hard_chunk
#
# Missing chunks are reported but do not abort the merge (so you can salvage a
# partial run). Check the printed warnings and re-submit only the missing array
# indices, e.g.  sbatch --array=2,3 ...   then re-run this script.

using DelimitedFiles
using Printf

# ── CLI parsing ─────────────────────────────────────────────────────────────
const DEFAULT_SUMMARY_NAME = "summary_overlap060_hard.tsv"
const DEFAULT_CHUNK_PREFIX = "summary_overlap060_hard_chunk"

function _parse_cli(args)
    outdir = ""
    total = 0
    summary_name = DEFAULT_SUMMARY_NAME
    chunk_prefix = DEFAULT_CHUNK_PREFIX
    dry_run = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--total"
            i < length(args) || error("--total requires an integer")
            total = parse(Int, args[i + 1]); i += 2; continue
        elseif startswith(arg, "--total=")
            total = parse(Int, split(arg, "=", limit=2)[2]); i += 1; continue
        elseif arg == "--summary-name"
            i < length(args) || error("--summary-name requires a value")
            summary_name = args[i + 1]; i += 2; continue
        elseif startswith(arg, "--summary-name=")
            summary_name = split(arg, "=", limit=2)[2]; i += 1; continue
        elseif arg == "--chunk-prefix"
            i < length(args) || error("--chunk-prefix requires a value")
            chunk_prefix = args[i + 1]; i += 2; continue
        elseif startswith(arg, "--chunk-prefix=")
            chunk_prefix = split(arg, "=", limit=2)[2]; i += 1; continue
        elseif arg == "--dry-run"
            dry_run = true; i += 1; continue
        elseif startswith(arg, "-")
            error("Unknown option: $arg")
        else
            outdir = arg; i += 1; continue
        end
    end
    isempty(outdir) && error("Usage: julia hpc/merge_chunks.jl <outdir> [--total N] [--summary-name NAME] [--chunk-prefix PREFIX] [--dry-run]")
    total >= 0 || error("--total must be >= 0")
    isdir(outdir) || error("outdir does not exist: $outdir")
    return outdir, total, summary_name, chunk_prefix, dry_run
end

# Match summary_overlap060_hard_chunkNNofMM.tsv and return (idx, total).
function _parse_chunk_filename(fn::AbstractString, prefix::AbstractString)
    b = basename(fn)
    startswith(b, prefix) || return nothing
    rest = b[nextind(b, sizeof(prefix)):end]   # "NNofMM.tsv"
    endswith(rest, ".tsv") || return nothing
    core = rest[1:prevind(rest, sizeof(rest) - 3)]  # strip ".tsv"
    m = match(r"^(\d+)of(\d+)$", core)
    m === nothing && return nothing
    return (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
end

# ── Main ────────────────────────────────────────────────────────────────────
outdir, total_hint, summary_name, chunk_prefix, dry_run = _parse_cli(ARGS)

const SUMMARY_PATH = joinpath(outdir, summary_name)
const SUMMARY_EXISTS = isfile(SUMMARY_PATH)

# Discover chunk shards.
shards = sort([f for f in readdir(outdir; join=true)
               if _parse_chunk_filename(f, chunk_prefix) !== nothing])

if isempty(shards)
    @warn "No chunk shards matching '$chunk_prefix*of*.tsv' found in $outdir"
    if SUMMARY_EXISTS
        println("A single (non-chunked) summary already exists: $SUMMARY_PATH — nothing to do.")
        exit(0)
    end
    exit(1)
end

idx_total = [_parse_chunk_filename(f, chunk_prefix) for f in shards]
inferred_total = maximum(t for (_, t) in idx_total)
total = total_hint > 0 ? total_hint : inferred_total
present_idx = sort(unique(i for (i, _) in idx_total))
expected_idx = 1:total

missing_idx = sort(setdiff(collect(expected_idx), present_idx))
if !isempty(missing_idx)
    @warn "Missing $(length(missing_idx)) of $total chunk(s): $(join(missing_idx, ","))"
    println("  → re-submit only those, e.g.:  sbatch --array=$(join(missing_idx, ",")) hpc/batch_array.sbatch")
    println("  → then re-run this merge. Continuing with $(length(present_idx)) present chunk(s).")
end

# Sanity: all shards agree on the declared total.
declared_totals = unique(t for (_, t) in idx_total)
length(declared_totals) > 1 &&
    @warn "Shards declare inconsistent totals: $(join(declared_totals, ",")); using --total/ inferred = $total"

@printf("Merging %d chunk shard(s) (of %d) → %s%s\n",
        length(shards), total, SUMMARY_PATH, dry_run ? "  [DRY-RUN]" : "")

# Read & concatenate inside a function so the locals (header, rows) have clean
# lexical scope (avoids the global soft-scope ambiguity). First shard
# contributes the header; the rest contribute data rows only (reuses the vcat
# pattern from test/summarize.jl:42-51).
function merge_shards(shards, outpath, n_missing::Int)
    header = String[]
    all_rows = Vector{Vector{String}}()
    for f in shards
        isfile(f) || (@warn "$f vanished mid-merge, skipping"; continue)
        raw = readdlm(f, '\t', String)
        isempty(header) && (header = vec(String.(raw[1, :])))
        for i in 2:size(raw, 1)
            push!(all_rows, String.(collect(raw[i, :])))
        end
    end
    # Sort by filepath (column 1) for a stable, readable output.
    sort!(all_rows; by=r -> isempty(r) ? "" : String(r[1]))

    n_data = length(all_rows)
    n_ok = count(r -> length(r) >= 2 && String(r[2]) == "ok", all_rows)
    @printf("  %d data rows (%d ok) across %d shard(s)\n", n_data, n_ok, length(shards))

    open(outpath, "w") do io
        isempty(header) || println(io, join(header, '\t'))
        for r in all_rows
            println(io, join(r, '\t'))
        end
    end
    println("Wrote $outpath")
    if n_missing > 0
        println("NOTE: $n_missing chunk(s) were missing — the merged summary is partial.")
        println("      Re-run after re-submitting the missing array indices to complete it.")
    end
    return n_data
end

if dry_run
    # Still count what we'd merge, but write nothing.
    local n = 0
    for f in shards
        isfile(f) || continue
        raw = readdlm(f, '\t', String)
        n += max(0, size(raw, 1) - 1)
    end
    @printf("  [DRY-RUN] would write %d data rows to %s\n", n, SUMMARY_PATH)
    exit(0)
end

n_missing = length(missing_idx)
merge_shards(shards, SUMMARY_PATH, n_missing)

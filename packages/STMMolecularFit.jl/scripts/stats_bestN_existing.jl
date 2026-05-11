#!/usr/bin/env julia
# Backward-compatible wrapper: summarize the newest enriched best-plot summary.
# Prefer using summarize_enriched.jl directly.

using Printf

function newest_summary(dir::String="results/best_plots")
    isdir(dir) || return nothing
    files = [joinpath(dir, f) for f in readdir(dir) if startswith(f, "summary_overlap060_hard") && endswith(f, ".tsv")]
    isempty(files) && (files = [joinpath(dir, f) for f in readdir(dir) if startswith(f, "summary_overlap060") && endswith(f, ".tsv")])
    isempty(files) && (files = [joinpath(dir, f) for f in readdir(dir) if startswith(f, "summary_") && endswith(f, ".tsv")])
    isempty(files) && return nothing
    return sort(files; by=mtime, rev=true)[1]
end

summary = length(ARGS) >= 1 ? ARGS[1] : newest_summary()
summary === nothing && error("No summary TSV found. Run centralize_best_plots.jl first, or pass a summary path.")

println("stats_bestN_existing.jl is deprecated; forwarding to summarize_enriched.jl")
println("summary: $summary")
const SUMMARY_OVERRIDE = summary
include(joinpath(@__DIR__, "summarize_enriched.jl"))

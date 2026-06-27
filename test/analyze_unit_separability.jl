#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# analyze_unit_separability.jl — Phase 1: test whether GlcNAc/GlcN are
# separable from Gaussian fit features (amplitude, sigma, integrated intensity).
#
# This is a DIAGNOSTIC script (post-hoc, label-free). It reads a per-lobe
# features TSV (produced by extract_lobe_features.jl) and tests:
#   1. Unimodal vs bimodal distribution of per-lobe features (kmeans k=1 vs
#      k=2, BIC comparison).
#   2. With --with-truth: cross with ground-truth sequences to measure the
#      actual separability (AUC per feature, clustering accuracy). This is
#      EXTERNAL EVALUATION ONLY — the truth is never used by the fitter.
#
# Features extracted per lobe (from the TSV):
#   amplitude             fitted Gaussian peak height
#   sigma_parallel_nm     fitted width along chain axis
#   sigma_perp_nm         fitted width perpendicular to axis
#   integrated            amplitude * sigma_parallel * sigma_perp (volume proxy)
#   spacing_prev_nm       distance to previous lobe
#   t_nm                  position along chain axis
#
# All features are per-file z-scored before pooling to remove file-to-file
# STM contrast variability.
#
# Examples:
#   # Label-free separability test
#   julia --project=. test/analyze_unit_separability.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --out results/unit_separability
#
#   # With truth cross-evaluation (external, post-hoc)
#   julia --project=. test/analyze_unit_separability.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --truth benchmarks/chitosan_240817_unit_sequences.tsv \
#       --out results/unit_separability
#
#   # Include residual features from Phase 1b
#   julia --project=. test/analyze_unit_separability.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --extra results/unit_separability/residual_features.tsv \
#       --truth benchmarks/chitosan_240817_unit_sequences.tsv \
#       --out results/unit_separability
# ──────────────────────────────────────────────────────────────────────────────

using Printf
using Dates
using Statistics
using StatsBase
using Clustering
using Random
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _read_tsv

const DEFAULT_FEATURES = "results/unit_separability/lobe_features.tsv"
const DEFAULT_OUT = "results/unit_separability"
const DEFAULT_TRUTH = "benchmarks/chitosan_240817_unit_sequences.tsv"

struct SepOptions
    features::String
    extra::Union{Nothing,String}
    out_dir::String
    truth::Union{Nothing,String}
    feat_names::Vector{String}
    primary_only::Bool
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    extra::Union{Nothing,String} = nothing
    out_dir = DEFAULT_OUT
    truth::Union{Nothing,String} = nothing
    feat_names = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]
    primary_only = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"
            features = args[i+1]; i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--extra"
            extra = args[i+1]; i += 2
        elseif startswith(arg, "--extra=")
            extra = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out_dir = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--truth"
            truth = args[i+1]; i += 2
        elseif startswith(arg, "--truth=")
            truth = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--features-list"
            feat_names = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--features-list=")
            feat_names = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--with-truth"
            truth = DEFAULT_TRUTH; i += 1
        elseif arg == "--primary-only"
            primary_only = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/analyze_unit_separability.jl [options]

            Options:
              --features PATH       Per-lobe features TSV from extract_lobe_features.jl [$(DEFAULT_FEATURES)]
              --extra PATH          Additional residual features TSV from extract_blob_residual_features.jl
              --out PATH            Output directory [$(DEFAULT_OUT)]
              --truth PATH          Ground-truth TSV for cross-evaluation
              --with-truth          Use default truth path: $(DEFAULT_TRUTH)
              --features-list LIST  Comma-separated feature names to use
                                     [amplitude,sigma_parallel_nm,sigma_perp_nm,integrated]
              --primary-only        Only use files with quality=clean or clean_target

            Reads:  <features>  (one row per lobe, from extract_lobe_features.jl)
                    <extra>     (optional, one row per lobe, from extract_blob_residual_features.jl)
            Writes: <out>/gaussian_features.tsv   (per-lobe features + z-scored + truth)
                    <out>/separability_report.txt  (analysis summary)
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    return SepOptions(features, extra, out_dir, truth, String.(strip.(feat_names)), primary_only)
end

_basename_file(s::AbstractString) = basename(strip(s))

function _parse_seq(s::AbstractString)
    t = strip(s)
    isempty(t) && return Int[]
    return [parse(Int, c) for c in t if c in "01"]
end

function _parse_float(s::AbstractString)
    t = strip(s)
    (isempty(t) || t in ("NA", "NaN", "nan")) && return NaN
    return parse(Float64, t)
end

# Compute AUC for a single feature vs binary labels (Mann-Whitney U).
function _roc_auc(scores::Vector{Float64}, labels::Vector{Int})
    n_pos = count(==(1), labels)
    n_neg = count(==(0), labels)
    (n_pos == 0 || n_neg == 0) && return NaN
    perm = sortperm(scores)
    ranks = zeros(Float64, length(scores))
    i = 1
    while i <= length(scores)
        j = i
        while j < length(scores) && scores[perm[j+1]] == scores[perm[i]]
            j += 1
        end
        avg_rank = (i + j) / 2.0
        for k in i:j
            ranks[perm[k]] = avg_rank
        end
        i = j + 1
    end
    sum_ranks_pos = sum(ranks[labels .== 1])
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return auc
end

# BIC for k-means. Clustering.jl expects observations in columns, so data is
# n_features × n_samples.
function _kmeans_bic(data::Matrix{Float64}, k::Int)
    nf, n = size(data)
    if k == 1
        center = mean(data, dims=2)
        sse = sum((data .- center).^2)
        n_params = nf
    else
        result = kmeans(data, k; maxiter=200, rng=MersenneTwister(42), display=:none)
        sse = result.totalcost
        n_params = k * nf + (k - 1)
    end
    bic = n * log(sse / n + 1e-30) + n_params * log(n)
    return bic, sse
end

function main()
    opt = _parse_cli(ARGS)
    mkpath(opt.out_dir)

    # Read truth if provided
    truth_by_file = Dict{String,Vector{Int}}()
    truth_quality = Dict{String,String}()
    if opt.truth !== nothing
        _, truth_rows = _read_tsv(opt.truth)
        for row in truth_rows
            file = _basename_file(row["file"])
            truth_by_file[file] = _parse_seq(row["sequence"])
            truth_quality[file] = get(row, "quality", "clean")
        end
    end

    # Read features
    _, feat_rows = _read_tsv(opt.features)
    isempty(feat_rows) && error("No feature rows in $(opt.features)")

    # Read extra features (residual) if provided
    extra_by_file_lobe = Dict{String,Dict{Int,Dict{String,Float64}}}()
    extra_cols = String[]
    if opt.extra !== nothing
        extra_header, extra_rows = _read_tsv(opt.extra)
        extra_cols = [c for c in extra_header if !(c in ["file", "lobe", "amplitude", "t_nm", "u_nm", "sigma_parallel_nm", "sigma_perp_nm"])]
        for row in extra_rows
            file = _basename_file(row["file"])
            lobe = parse(Int, row["lobe"])
            d = Dict{String,Float64}()
            for c in extra_cols
                d[c] = _parse_float(get(row, c, "NA"))
            end
            get!(extra_by_file_lobe, file, Dict{Int,Dict{String,Float64}}())[lobe] = d
        end
    end

    # Build feature table
    records = []
    for row in feat_rows
        file = _basename_file(row["file"])
        lobe = parse(Int, row["lobe"])
        if opt.primary_only && haskey(truth_quality, file) && truth_quality[file] == "poor_quality"
            continue
        end
        amp = _parse_float(row["amplitude"])
        spar = _parse_float(row["sigma_parallel_nm"])
        sperp = _parse_float(row["sigma_perp_nm"])
        integrated = amp * spar * sperp
        t_nm = _parse_float(row["t_nm"])
        spacing_prev = _parse_float(row["spacing_prev_nm"])

        feat = Dict{String,Float64}(
            "amplitude" => amp,
            "sigma_parallel_nm" => spar,
            "sigma_perp_nm" => sperp,
            "integrated" => integrated,
            "spacing_prev_nm" => spacing_prev,
            "t_nm" => t_nm,
        )
        for (k, v) in row
            k in ("file", "N", "lobe", "source") && continue
            haskey(feat, k) && continue
            parsed = tryparse(Float64, strip(v))
            parsed === nothing || (feat[k] = parsed)
        end
        # Merge extra features
        if haskey(extra_by_file_lobe, file) && haskey(extra_by_file_lobe[file], lobe)
            merge!(feat, extra_by_file_lobe[file][lobe])
        end
        push!(records, (file, lobe, feat))
    end

    isempty(records) && error("No lobe records after filtering")
    n_total = length(records)
    n_files = length(unique(r[1] for r in records))
    println("Collected ", n_total, " lobes from ", n_files, " files")

    # Select features (including extra if present)
    feat_names = copy(opt.feat_names)
    if !isempty(extra_cols)
        for c in extra_cols
            c in feat_names || push!(feat_names, c)
        end
    end
    nf = length(feat_names)

    # Per-file z-scoring
    files_seen = unique(r[1] for r in records)
    file_stats = Dict{String,Dict{String,Tuple{Float64,Float64}}}()
    for file in files_seen
        file_records = [r for r in records if r[1] == file]
        stats = Dict{String,Tuple{Float64,Float64}}()
        for fn in feat_names
            vals = [get(r[3], fn, NaN) for r in file_records]
            vals = filter(isfinite, vals)
            isempty(vals) && (stats[fn] = (0.0, 1.0); continue)
            mu = mean(vals)
            sigma = std(vals) > 0 ? std(vals) : 1.0
            stats[fn] = (mu, sigma)
        end
        file_stats[file] = stats
    end

    # Build z-scored feature matrix (skip rows with NaN in any selected feature)
    valid_mask = trues(n_total)
    X = zeros(n_total, nf)
    for (i, (file, lobe, feat)) in enumerate(records)
        for (j, fn) in enumerate(feat_names)
            v = get(feat, fn, NaN)
            if !isfinite(v)
                valid_mask[i] = false
                X[i, j] = 0.0
            else
                mu, sigma = file_stats[file][fn]
                X[i, j] = sigma > 0 ? (v - mu) / sigma : 0.0
            end
        end
    end

    n_valid = sum(valid_mask)
    X_valid = X[valid_mask, :]
    records_valid = records[valid_mask]
    println("Valid rows (no NaN in selected features): ", n_valid, "/", n_total)

    # Write per-lobe features TSV
    features_path = joinpath(opt.out_dir, "gaussian_features.tsv")
    open(features_path, "w") do io
        cols = ["file", "lobe", [fn for fn in feat_names]..., ["z_$(fn)" for fn in feat_names]..., "truth_label"]
        println(io, join(cols, '\t'))
        for (i, (file, lobe, feat)) in enumerate(records)
            truth_seq = get(truth_by_file, file, Int[])
            truth_label = lobe <= length(truth_seq) ? string(truth_seq[lobe]) : ""
            raw_vals = join([isnan(get(feat, fn, NaN)) ? "NA" : @sprintf("%.6f", get(feat, fn, NaN)) for fn in feat_names], '\t')
            zvals = valid_mask[i] ? join([@sprintf("%.6f", X[i, j]) for j in 1:nf], '\t') :
                                     join(["NA" for _ in 1:nf], '\t')
            println(io, join([file, string(lobe), raw_vals, zvals, truth_label], '\t'))
        end
    end
    println("Wrote features: ", features_path)

    # ── Label-free: test k=1 vs k=2 ──
    X_t = permutedims(X_valid)
    bic1, sse1 = _kmeans_bic(X_t, 1)
    bic2, sse2 = _kmeans_bic(X_t, 2)
    delta_bic = bic1 - bic2

    result_k2 = kmeans(X_t, 2; maxiter=200, rng=MersenneTwister(42), display=:none)
    assignments = result_k2.assignments
    c1_size = count(==(1), assignments)
    c2_size = count(==(2), assignments)
    centroid_dist = norm(result_k2.centers[:, 1] .- result_k2.centers[:, 2])

    # ── With truth: cross-evaluation ──
    truth_labels_all = Int[]
    truth_mask = falses(n_total)
    for (i, (file, lobe, feat)) in enumerate(records)
        seq = get(truth_by_file, file, Int[])
        if lobe <= length(seq)
            truth_labels_all = push!(truth_labels_all, seq[lobe])
            truth_mask[i] = true
        end
    end

    # Build report
    report_path = joinpath(opt.out_dir, "separability_report.txt")
    open(report_path, "w") do io
        println(io, "═")
        println(io, "Unit Separability Report")
        println(io, "═")
        println(io, "Date: ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
        println(io, "Features TSV: ", opt.features)
        opt.extra !== nothing && println(io, "Extra TSV:    ", opt.extra)
        println(io, "Features: ", join(feat_names, ", "))
        println(io, "Lobes: ", n_total, " from ", n_files, " files (", n_valid, " valid)")
        println(io, "Per-file z-scoring: applied")
        println(io)

        println(io, "── Label-free: unimodal vs bimodal test ──")
        println(io, "k=1 BIC: ", @sprintf("%.1f", bic1), "  (SSE=", @sprintf("%.4f", sse1), ")")
        println(io, "k=2 BIC: ", @sprintf("%.1f", bic2), "  (SSE=", @sprintf("%.4f", sse2), ")")
        println(io, "ΔBIC (k=1 - k=2): ", @sprintf("%.1f", delta_bic),
                  delta_bic > 10 ? "  → BIMODAL (k=2 strongly preferred)" :
                  delta_bic > 0 ? "  → weakly bimodal" :
                  "  → UNIMODAL (k=1 preferred)")
        println(io, "Cluster sizes: [", c1_size, ", ", c2_size, "]")
        println(io, "Centroid distance: ", @sprintf("%.3f", centroid_dist))
        println(io)

        if !isempty(truth_labels_all) && opt.truth !== nothing
            # Rebuild truth labels for valid rows only
            truth_valid = Int[]
            for (i, (file, lobe, feat)) in enumerate(records)
                if valid_mask[i]
                    seq = get(truth_by_file, file, Int[])
                    if lobe <= length(seq)
                        push!(truth_valid, seq[lobe])
                    end
                end
            end

            println(io, "── With truth: cross-evaluation (external, post-hoc) ──")
            println(io, "Labeled valid lobes: ", length(truth_valid), " / ", n_valid)
            n_pos = count(==(1), truth_valid)
            n_neg = count(==(0), truth_valid)
            println(io, "Truth: GlcNAc(1)=", n_pos, "  GlcN(0)=", n_neg)
            println(io)

            # AUC per feature
            println(io, "AUC per feature (pooled, z-scored):")
            X_labeled = X_valid[1:length(truth_valid), :]  # need to align properly
            # Rebuild X_labeled properly
            X_labeled = zeros(length(truth_valid), nf)
            idx = 0
            for (i, (file, lobe, feat)) in enumerate(records)
                if valid_mask[i]
                    seq = get(truth_by_file, file, Int[])
                    if lobe <= length(seq)
                        idx += 1
                        X_labeled[idx, :] = X[i, :]
                    end
                end
            end

            for (j, fn) in enumerate(feat_names)
                auc = _roc_auc(X_labeled[:, j], truth_valid)
                auc_sym = isnan(auc) ? "NA" : (auc < 0.5 ? @sprintf("%.3f (inv=%.3f)", auc, 1-auc) : @sprintf("%.3f", auc))
                println(io, "  ", fn, ": ", auc_sym)
            end
            println(io)

            # Clustering accuracy
            X_lab_t = permutedims(X_labeled)
            result_lab = kmeans(X_lab_t, 2; maxiter=200, rng=MersenneTwister(42), display=:none)
            clusters = result_lab.assignments

            # Physical convention: GlcNAc = higher amplitude
            amp_idx = findfirst(==("amplitude"), feat_names)
            glcnac_cluster = 2
            if amp_idx !== nothing
                amp_labeled = X_labeled[:, amp_idx]
                mean_amp_c1 = mean(amp_labeled[clusters .== 1])
                mean_amp_c2 = mean(amp_labeled[clusters .== 2])
                glcnac_cluster = mean_amp_c1 >= mean_amp_c2 ? 1 : 2
            end
            pred = [c == glcnac_cluster ? 1 : 0 for c in clusters]

            phys_correct = count(pred[i] == truth_valid[i] for i in 1:length(truth_valid))
            phys_acc = phys_correct / length(truth_valid)
            acc_flip = count((1 .- pred)[i] == truth_valid[i] for i in 1:length(truth_valid)) / length(truth_valid)
            oracle_acc = max(phys_acc, acc_flip)

            println(io, "Clustering accuracy (pooled 2-means):")
            println(io, "  Physical (GlcNAc=amp max):  ", @sprintf("%.1f%%", 100*phys_acc),
                      "  (", phys_correct, "/", length(truth_valid), ")")
            println(io, "  Oracle (best flip):         ", @sprintf("%.1f%%", 100*oracle_acc))
            println(io, "  Gap: ", @sprintf("%.1f%%", 100*(oracle_acc - phys_acc)))
            println(io)
            println(io, "Interpretation:")
            if phys_acc > 0.75
                println(io, "  → Good separability: features discriminate GlcNAc/GlcN.")
            elseif phys_acc > 0.60
                println(io, "  → Weak separability: signal present but noisy.")
            else
                println(io, "  → Poor separability: GlcNAc/GlcN not distinguishable from these features.")
                println(io, "    Consider Phase 1b (residual non-Gaussian features) or different STM conditions.")
            end
        else
            println(io, "── No truth provided: skipping cross-evaluation ──")
            println(io, "Run with --truth to measure AUC and clustering accuracy.")
        end
        println(io)
        println(io, "── Output files ──")
        println(io, "  ", features_path, " — per-lobe features (raw + z-scored + truth label)")
        println(io, "  ", report_path, " — this report")
    end

    println("Wrote report: ", report_path)

    # Summary
    println()
    println("Separability summary:")
    println("  ΔBIC (k=1 - k=2): ", @sprintf("%.1f", delta_bic),
              delta_bic > 10 ? "  → BIMODAL" :
              delta_bic > 0 ? "  → weakly bimodal" :
              "  → UNIMODAL")
    if opt.truth !== nothing && @isdefined(phys_acc)
        println("  Physical accuracy: ", @sprintf("%.1f%%", 100*phys_acc))
        println("  Oracle accuracy:   ", @sprintf("%.1f%%", 100*oracle_acc))
    end
end

main()

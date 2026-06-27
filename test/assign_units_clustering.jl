#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# assign_units_clustering.jl — Phase 3: assign GlcNAc/GlcN (0/1) to each lobe
# via unsupervised clustering of per-lobe features.
#
# This is the main unit-assignment script. It reads the per-lobe features TSV
# (from extract_lobe_features.jl) and optionally the residual features TSV
# (from extract_blob_residual_features.jl), then:
#   1. Pools all lobes, per-file z-scores the features.
#   2. Runs k-means k=1 vs k=2 (BIC) to test if two populations exist.
#   3. If k=2 is preferred (or forced), assigns each lobe to a cluster.
#   4. Maps clusters to {GlcN=0, GlcNAc=1} using the physical convention:
#      GlcNAc = higher mean amplitude (acetyl group is larger → brighter).
#   5. Writes the predicted sequence per file.
#
# The assignment is LABEL-FREE: it uses only the fitted features, never the
# ground-truth sequence. The physical mapping (amplitude) is a prior, not a
# label.
#
# Examples:
#   # Basic: Gaussian features only
#   julia --project=. test/assign_units_clustering.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --out results/unit_assignment/assigned_sequences.tsv
#
#   # With residual features (Phase 1b)
#   julia --project=. test/assign_units_clustering.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --extra results/unit_separability/residual_features.tsv \
#       --out results/unit_assignment/assigned_sequences.tsv
#
#   # Force k=2 even if BIC prefers k=1
#   julia --project=. test/assign_units_clustering.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --force-k2 --out results/unit_assignment/assigned_sequences.tsv
# ──────────────────────────────────────────────────────────────────────────────

using Printf
using Statistics
using StatsBase
using Clustering
using Random
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_FEATURES = "results/unit_separability/lobe_features.tsv"
const DEFAULT_OUT = "results/unit_assignment/assigned_sequences.tsv"

struct AssignOptions
    features::String
    extra::Union{Nothing,String}
    out_tsv::String
    feat_names::Vector{String}
    force_k2::Bool
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    extra::Union{Nothing,String} = nothing
    out_tsv = DEFAULT_OUT
    feat_names = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]
    force_k2 = false

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
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--features-list"
            feat_names = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--features-list=")
            feat_names = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--force-k2"
            force_k2 = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/assign_units_clustering.jl [options]

            Options:
              --features PATH       Per-lobe features TSV [$(DEFAULT_FEATURES)]
              --extra PATH          Additional residual features TSV
              --out PATH            Output TSV [$(DEFAULT_OUT)]
              --features-list LIST  Comma-separated feature names
              --force-k2            Force k=2 even if BIC prefers k=1
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    return AssignOptions(features, extra, out_tsv, String.(strip.(feat_names)), force_k2)
end

_basename_file(s::AbstractString) = basename(strip(s))

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
    mkpath(dirname(opt.out_tsv))

    # Read features
    _, feat_rows = _read_tsv(opt.features)
    isempty(feat_rows) && error("No feature rows in $(opt.features)")

    # Read extra features if provided
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
                d[c] = _parse_f(get(row, c, "NA"))
            end
            get!(extra_by_file_lobe, file, Dict{Int,Dict{String,Float64}}())[lobe] = d
        end
    end

    # Build records
    records = []
    for row in feat_rows
        file = _basename_file(row["file"])
        lobe = parse(Int, row["lobe"])
        amp = _parse_f(row["amplitude"])
        spar = _parse_f(row["sigma_parallel_nm"])
        sperp = _parse_f(row["sigma_perp_nm"])
        integrated = amp * spar * sperp

        feat = Dict{String,Float64}(
            "amplitude" => amp,
            "sigma_parallel_nm" => spar,
            "sigma_perp_nm" => sperp,
            "integrated" => integrated,
        )
        for (k, v) in row
            k in ("file", "N", "lobe", "source") && continue
            haskey(feat, k) && continue
            parsed = tryparse(Float64, strip(v))
            parsed === nothing || (feat[k] = parsed)
        end
        if haskey(extra_by_file_lobe, file) && haskey(extra_by_file_lobe[file], lobe)
            merge!(feat, extra_by_file_lobe[file][lobe])
        end
        push!(records, (file, lobe, feat))
    end

    # Feature names
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
            vals = filter(isfinite, [get(r[3], fn, NaN) for r in file_records])
            isempty(vals) && (stats[fn] = (0.0, 1.0); continue)
            mu = mean(vals)
            sigma = std(vals) > 0 ? std(vals) : 1.0
            stats[fn] = (mu, sigma)
        end
        file_stats[file] = stats
    end

    # Build z-scored matrix, filtering NaN rows
    valid_mask = trues(length(records))
    X = zeros(length(records), nf)
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
    println("Valid lobes: ", n_valid, "/", length(records))

    # k=1 vs k=2 BIC
    X_t = permutedims(X_valid)
    bic1, _ = _kmeans_bic(X_t, 1)
    bic2, _ = _kmeans_bic(X_t, 2)
    delta_bic = bic1 - bic2
    println("ΔBIC (k=1 - k=2): ", @sprintf("%.1f", delta_bic),
              delta_bic > 10 ? "  → bimodal" :
              delta_bic > 0 ? "  → weakly bimodal" :
              "  → unimodal")

    use_k2 = opt.force_k2 || delta_bic > 0

    if !use_k2
        println("k=1 preferred by BIC. Using --force-k2=false → all lobes assigned to same type.")
        # All same type — assign all to 1 (GlcNAc) by convention
        assignments = ones(Int, n_valid)
        result_k2 = nothing
    else
        result_k2 = kmeans(X_t, 2; maxiter=200, rng=MersenneTwister(42), display=:none)
        assignments = result_k2.assignments
    end

    # Physical mapping: GlcNAc (1) = higher mean raw amplitude. Use the raw
    # amplitude even when amplitude is not part of the clustering feature list.
    glcnac_cluster = 1
    if result_k2 !== nothing
        amp_valid = [get(r[3], "amplitude", NaN) for r in records_valid]
        mean_amp_c1 = mean(filter(isfinite, amp_valid[assignments .== 1]))
        mean_amp_c2 = mean(filter(isfinite, amp_valid[assignments .== 2]))
        glcnac_cluster = mean_amp_c1 >= mean_amp_c2 ? 1 : 2
        labels = [a == glcnac_cluster ? 1 : 0 for a in assignments]
    else
        labels = assignments  # all 1s if k=1
    end

    # Write per-lobe predictions
    open(opt.out_tsv, "w") do io
        println(io, join([
            "file", "lobe", "amplitude", "predicted", "physical_label",
            "cluster", "P_GlcNAc"
        ], '\t'))

        for (i, (file, lobe, feat)) in enumerate(records_valid)
            amp = get(feat, "amplitude", NaN)
            label = labels[i]
            cluster = assignments[i]
            # Rough posterior: distance to each centroid (softmax)
            p_glcnac = NaN
            if result_k2 !== nothing
                c1 = result_k2.centers[:, glcnac_cluster]
                c2 = result_k2.centers[:, 3 - glcnac_cluster]
                d1 = norm(X_valid[i, :] .- c1)
                d2 = norm(X_valid[i, :] .- c2)
                p_glcnac = 1.0 / (1.0 + exp(d1 - d2))
            end
            println(io, join([
                file, string(lobe),
                isnan(amp) ? "NA" : @sprintf("%.8e", amp),
                string(label),
                label == 1 ? "GlcNAc" : "GlcN",
                string(cluster),
                isnan(p_glcnac) ? "NA" : @sprintf("%.4f", p_glcnac),
            ], '\t'))
        end
    end

    println("Wrote: ", opt.out_tsv)

    # Print per-file sequences
    println("\n── Predicted sequences ──")
    by_file = Dict{String,Vector{Tuple{Int,Int,Float64}}}()
    for (i, (file, lobe, feat)) in enumerate(records_valid)
        push!(get!(by_file, file, Tuple{Int,Int,Float64}[]),
              (lobe, labels[i], get(feat, "amplitude", NaN)))
    end

    for file in sort(collect(keys(by_file)))
        entries = sort(by_file[file], by=x->x[1])
        seq = join([string(e[2]) for e in entries])
        n = length(entries)
        n_glcnac = count(e[2] == 1 for e in entries)
        n_glcn = n - n_glcnac
        println("  ", file, "  N=", n, "  seq=", seq,
                "  (GlcNAc=", n_glcnac, ", GlcN=", n_glcn, ")")
    end
end

main()

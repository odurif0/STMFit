#!/usr/bin/env julia

# Build label-free 0/1/? unit-assignment predictions via GMM (2-component full-covariance EM).
#
# This script is deliberately prediction-only: it never reads benchmark truth,
# expected sequences, or composition priors. Cluster labels are mapped by the
# physical convention already used by the grader: the higher-amplitude cluster is
# GlcNAc (1). Use grade_unit_assignment.jl or report_unit_assignment_benchmark.jl
# only after this TSV has been written.

using Clustering
using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv

const DEFAULT_FEATURES = "results/unit_separability/lobe_features_selectedN_primary_local.tsv"
const DEFAULT_SPLIT = ""
const DEFAULT_PATCHES = ""
const DEFAULT_OUT = "results/unit_assignment/labelfree_gmm_predictions.tsv"

struct Options
    features::String
    split_features::String
    patches::String
    out_tsv::String
    view_specs::Vector{Pair{String,Vector{String}}}
    first_seed::Int
    n_seeds::Int
    interactions::Bool
    selftrain::Int
end

mutable struct LobeRecord
    file::String
    lobe::Int
    amplitude::Float64
    features::Dict{String,Float64}
end

function _parse_view(s::AbstractString)
    parts = split(String(s), '='; limit=2)
    length(parts) == 2 || error("--view must be NAME=feature1,feature2,...")
    name = strip(parts[1])
    isempty(name) && error("--view name is empty")
    features = [strip(f) for f in split(parts[2], ',') if !isempty(strip(f))]
    isempty(features) && error("--view $name has no features")
    return name => features
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    split_features = DEFAULT_SPLIT
    patches = DEFAULT_PATCHES
    out_tsv = DEFAULT_OUT
    view_specs = Pair{String,Vector{String}}[]
    first_seed = 0
    n_seeds = 20
    interactions = false
    selftrain = 0

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"
            features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--split-features"
            split_features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--split-features=")
            split_features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--patches"
            patches = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--patches=")
            patches = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--view"
            push!(view_specs, _parse_view(_arg_value(args, i, arg))); i += 2
        elseif startswith(arg, "--view=")
            push!(view_specs, _parse_view(split(arg, "="; limit=2)[2])); i += 1
        elseif arg == "--first-seed"
            first_seed = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--first-seed=")
            first_seed = parse(Int, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--seeds"
            n_seeds = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--seeds=")
            n_seeds = parse(Int, split(arg, "="; limit=2)[2]); i += 1
        elseif arg == "--interactions"
            interactions = true; i += 1
        elseif arg == "--no-interactions"
            interactions = false; i += 1
        elseif arg == "--selftrain"
            selftrain = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--selftrain=")
            selftrain = parse(Int, split(arg, "="; limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_labelfree_gmm_predictions.jl [options]

            Options:
              --features PATH        Main per-lobe feature TSV [$(DEFAULT_FEATURES)]
              --split-features PATH  Optional split-width feature TSV; adds split_log_skew
              --patches PATH         Optional backward patch TSV; adds bwd_neg_* descriptors
              --out PATH             Output prediction TSV [$(DEFAULT_OUT)]
              --view NAME=LIST       Feature view to cluster. Repeatable. If omitted,
                                     sensible label-free defaults are chosen from
                                     available columns.
              --first-seed INT       First k-means seed [0]
              --seeds INT            Number of k-means seeds [20]
              --interactions         Append pairwise products within each view after
                                     per-file z-scoring.
              --no-interactions      Disable the pairwise-product expansion (overrides
                                     an earlier --interactions; keeps views low-dim).
              --selftrain INT        After each seed's EM fit, run INT hard
                                     reassignment iterations using Mahalanobis
                                     distance to the cluster means (GMM seeds ->
                                     Mahalanobis self-training). [0]

            Output columns: file, lobe, predicted, confidence, amplitude,
            probability_1, views_used, invalid_reason.

            Label-free constraint: this script reads no truth sequence, expected N,
            control motif, or composition count. Benchmark grading is a separate
            post-hoc step.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    n_seeds > 0 || error("--seeds must be positive")
    selftrain >= 0 || error("--selftrain must be >= 0")
    isfile(features) || error("Feature TSV not found: $features")
    !isempty(split_features) && !isfile(split_features) && error("Split feature TSV not found: $split_features")
    !isempty(patches) && !isfile(patches) && error("Patch TSV not found: $patches")
    return Options(features, split_features, patches, out_tsv, view_specs,
                   first_seed, n_seeds, interactions, selftrain)
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

_key(file::AbstractString, lobe::Integer) = (basename(strip(file)), Int(lobe))

function _row_file(row::Dict{String,String}, source::AbstractString)
    haskey(row, "file") || error("$source missing file column")
    file = basename(strip(row["file"]))
    isempty(file) && error("$source has empty file value")
    return file
end

function _row_lobe(row::Dict{String,String}, source::AbstractString)
    haskey(row, "lobe") || error("$source missing lobe column")
    lobe = tryparse(Int, strip(row["lobe"]))
    lobe === nothing && error("$source has invalid lobe index: $(row["lobe"])")
    lobe >= 1 || error("$source lobe index must be >= 1, got $lobe")
    return lobe
end

function _parse_f(s)
    t = strip(String(s))
    (isempty(t) || t in ("NA", "NaN", "nan", "?")) && return NaN
    return parse(Float64, t)
end

function _record_features(row::Dict{String,String})
    feat = Dict{String,Float64}()
    for (k, v) in row
        k in ("file", "N", "lobe", "source") && continue
        parsed = tryparse(Float64, strip(v))
        parsed === nothing && continue
        feat[k] = parsed
    end
    if !haskey(feat, "integrated") && all(haskey(feat, k) for k in ("amplitude", "sigma_parallel_nm", "sigma_perp_nm"))
        feat["integrated"] = feat["amplitude"] * feat["sigma_parallel_nm"] * feat["sigma_perp_nm"]
    end
    return feat
end

function _load_records(path::String)
    _, rows = _read_tsv(path)
    isempty(rows) && error("No rows in $path")
    records = LobeRecord[]
    seen = Set{Tuple{String,Int}}()
    for row in rows
        for col in ("file", "lobe", "amplitude")
            haskey(row, col) || error("Feature TSV missing column: $col")
        end
        file = _row_file(row, "Feature TSV $path")
        lobe = _row_lobe(row, "Feature TSV $path")
        key = _key(file, lobe)
        key in seen && error("Feature TSV has duplicate row for $(key[1]) lobe $(key[2]): $path")
        push!(seen, key)
        amp = _parse_f(row["amplitude"])
        isfinite(amp) || error("Feature TSV has non-finite amplitude for $file lobe $lobe: $(row["amplitude"])")
        push!(records, LobeRecord(file, lobe, amp, _record_features(row)))
    end
    sort!(records, by=r -> (r.file, r.lobe))
    return records
end

function _merge_split!(records::Vector{LobeRecord}, path::String)
    isempty(path) && return nothing
    _, rows = _read_tsv(path)
    by_key = Dict{Tuple{String,Int},Float64}()
    for row in rows
        haskey(row, "skew_ratio") || error("Split TSV missing skew_ratio column")
        file = _row_file(row, "Split TSV $path")
        lobe = _row_lobe(row, "Split TSV $path")
        key = _key(file, lobe)
        haskey(by_key, key) && error("Split TSV has duplicate row for $(key[1]) lobe $(key[2]): $path")
        skew = _parse_f(row["skew_ratio"])
        by_key[key] = isfinite(skew) && skew > 0 ? log(skew) : NaN
    end
    for rec in records
        rec.features["split_log_skew"] = get(by_key, _key(rec.file, rec.lobe), NaN)
    end
    return nothing
end

function _patch_columns(header::Vector{String}, prefix::String)
    cols = [c for c in header if startswith(c, prefix)]
    isempty(cols) && return String[]
    side = round(Int, sqrt(length(cols)))
    side * side == length(cols) || error("Patch columns for prefix $prefix are not a square grid: $(length(cols))")
    return cols
end

function _negative_moment(vals::Vector{Float64}, coords::Vector{Float64})
    side = length(coords)
    weights = max.(-vals, 0.0)
    total = sum(weights)
    total <= eps(Float64) && return (com_t=NaN, diag45=NaN, diag135=NaN)
    com_t = 0.0
    diag45 = 0.0
    diag135 = 0.0
    idx = 1
    for u in coords, t in coords
        w = weights[idx]
        com_t += w * t
        diag45 += w * ((t + u) / sqrt(2.0))
        diag135 += w * ((t - u) / sqrt(2.0))
        idx += 1
    end
    return (com_t=com_t / total, diag45=diag45 / total, diag135=diag135 / total)
end

function _merge_patches!(records::Vector{LobeRecord}, path::String)
    isempty(path) && return nothing
    header, rows = _read_tsv(path)
    cols = _patch_columns(String.(header), "bwd_res_p")
    isempty(cols) && error("Patch TSV has no bwd_res_pNNN columns: $path")
    side = round(Int, sqrt(length(cols)))
    coords = side == 1 ? [0.0] : collect(range(-1.0, 1.0; length=side))
    by_key = Dict{Tuple{String,Int},NamedTuple{(:com_t,:diag45,:diag135),Tuple{Float64,Float64,Float64}}}()
    for row in rows
        file = _row_file(row, "Patch TSV $path")
        lobe = _row_lobe(row, "Patch TSV $path")
        key = _key(file, lobe)
        haskey(by_key, key) && error("Patch TSV has duplicate row for $(key[1]) lobe $(key[2]): $path")
        vals = [_parse_f(row[c]) for c in cols]
        by_key[key] = _negative_moment(vals, coords)
    end
    for rec in records
        d = get(by_key, _key(rec.file, rec.lobe), (com_t=NaN, diag45=NaN, diag135=NaN))
        rec.features["bwd_neg_com_t"] = d.com_t
        rec.features["bwd_neg_diag45"] = d.diag45
        rec.features["bwd_neg_diag135"] = d.diag135
    end
    return nothing
end

function _available_features(records::Vector{LobeRecord})
    names = Set{String}()
    for rec in records
        union!(names, keys(rec.features))
    end
    return names
end

function _default_views(records::Vector{LobeRecord})
    available = _available_features(records)
    local_base = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"]
    gaussian_base = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]
    base = all(in(available), local_base) ? local_base : gaussian_base
    all(in(available), base) || error("No default feature base is available; pass --view explicitly")

    views = Pair{String,Vector{String}}[]
    pushed = false
    for (name, extra) in (("base_bwd_neg_com_t", "bwd_neg_com_t"),
                          ("base_bwd_neg_diag45", "bwd_neg_diag45"),
                          ("base_split_log_skew", "split_log_skew"))
        if extra in available
            push!(views, name => vcat(base, [extra]))
            pushed = true
        end
    end
    pushed || push!(views, "base" => base)
    return views
end

function _standardized_matrix(records::Vector{LobeRecord}, features::Vector{String}; interactions::Bool=false)
    n = length(records)
    p = length(features)
    raw = fill(NaN, n, p)
    valid = trues(n)
    for (i, rec) in enumerate(records)
        for (j, fn) in enumerate(features)
            v = get(rec.features, fn, NaN)
            raw[i, j] = v
            isfinite(v) || (valid[i] = false)
        end
    end

    z = fill(NaN, n, p)
    files = sort(unique(rec.file for rec in records))
    for file in files
        idxs = findall(i -> records[i].file == file, 1:n)
        for j in 1:p
            vals = raw[idxs, j]
            good = filter(isfinite, vals)
            if isempty(good)
                z[idxs, j] .= NaN
            else
                μ = mean(good)
                σ = std(good)
                σ = σ > 0 ? σ : 1.0
                for i in idxs
                    z[i, j] = isfinite(raw[i, j]) ? (raw[i, j] - μ) / σ : NaN
                end
            end
        end
    end

    if interactions && p > 1
        extra = Matrix{Float64}(undef, n, p * (p - 1) ÷ 2)
        col = 1
        for a in 1:(p-1), b in (a+1):p
            extra[:, col] = z[:, a] .* z[:, b]
            col += 1
        end
        z = hcat(z, extra)
    end

    for i in 1:n
        all(isfinite, z[i, :]) || (valid[i] = false)
    end
    return z, valid
end

function _gmm_log_density(x::AbstractVector, mu::AbstractVector, Sigma::AbstractMatrix)
    p = length(x)
    d = x - mu
    L = cholesky(Symmetric(Sigma) + 1e-8 * I).L
    z = L \ d
    return -0.5 * (p * log(2π) + 2 * sum(log.(diag(L))) + dot(z, z))
end

# EM fit of a 2-component full-covariance GMM, k-means initialized.
# Returns (means, covariances, weights).
function _gmm_fit(X::Matrix{Float64}, seed::Int; max_iter::Int=200, tol::Float64=1e-6)
    p, n = size(X)
    rng = MersenneTwister(seed)
    k = 2
    km = kmeans(X, k; maxiter=100, rng=rng, display=:none)
    means = zeros(p, k)
    covs = [Matrix{Float64}(I, p, p) for _ in 1:k]
    weights = fill(1.0 / k, k)
    for c in 1:k
        members = findall(km.assignments .== c)
        isempty(members) && continue
        means[:, c] = sum(X[:, members]; dims=2) / length(members)
        centered = X[:, members] .- means[:, c]
        covs[c] = centered * centered' / length(members)
        covs[c] += 1e-6 * I
    end
    prev_ll = -Inf
    for iter in 1:max_iter
        # E-step in log space
        log_resp = zeros(n, k)
        for j in 1:k
            for i in 1:n
                log_resp[i, j] = log(weights[j]) + _gmm_log_density(X[:, i], means[:, j], covs[j])
            end
        end
        ll = 0.0
        for i in 1:n
            m = maximum(log_resp[i, :])
            ll += m + log(sum(exp.(log_resp[i, :] .- m)))
        end
        abs(ll - prev_ll) < tol * (1 + abs(prev_ll)) && break
        prev_ll = ll
        for i in 1:n
            m = maximum(log_resp[i, :])
            log_resp[i, :] .-= m
            log_resp[i, :] .-= log(sum(exp.(log_resp[i, :])))
        end
        resp = exp.(log_resp)
        # M-step
        for j in 1:k
            nk = sum(resp[:, j])
            nk > 1e-12 || continue
            weights[j] = nk / n
            means[:, j] = X * resp[:, j] / nk
            centered = X .- means[:, j]
            covs[j] = (centered * Diagonal(resp[:, j]) * centered') / nk
            covs[j] += 1e-6 * I
        end
    end
    return means, covs, weights
end

# Hard self-training: reassign every point to the closest cluster by Mahalanobis
# distance, then re-estimate means/covariances from the hard memberships.
# Repeats `iters` times. Returns updated (means, covariances, weights).
function _mahalanobis_self_train(X::Matrix{Float64}, means::Matrix{Float64},
                                 covs::Vector{Matrix{Float64}}, weights::Vector{Float64};
                                 iters::Int=5)
    p, n = size(X)
    k = size(means, 2)
    for _ in 1:iters
        d2 = zeros(n, k)
        for j in 1:k
            L = cholesky(Symmetric(covs[j]) + 1e-8 * I).L
            for i in 1:n
                z = L \ (X[:, i] .- means[:, j])
                d2[i, j] = dot(z, z)
            end
        end
        assign = [argmin(d2[i, :]) for i in 1:n]
        for j in 1:k
            members = findall(==(j), assign)
            isempty(members) && continue
            nk = length(members)
            weights[j] = nk / n
            means[:, j] = sum(X[:, members]; dims=2) / nk
            centered = X[:, members] .- means[:, j]
            covs[j] = centered * centered' / nk + 1e-6 * I
        end
    end
    return means, covs, weights
end

function _view_probability(records::Vector{LobeRecord}, features::Vector{String}, opt::Options)
    X, valid = _standardized_matrix(records, features; interactions=opt.interactions)
    idxs = findall(valid)
    length(idxs) >= 2 || return fill(NaN, length(records))

    probs = fill(NaN, length(records))
    votes = zeros(Float64, length(records))
    counts = zeros(Int, length(records))
    data = permutedims(X[idxs, :])

    for seed in opt.first_seed:(opt.first_seed + opt.n_seeds - 1)
        means, covs, weights = _gmm_fit(data, seed)
        if opt.selftrain > 0
            means, covs, weights = _mahalanobis_self_train(data, means, covs, weights;
                                                           iters=opt.selftrain)
        end
        # Physical mapping: GlcNAc = higher-amplitude cluster.
        cluster_amp = Dict{Int,Vector{Float64}}(1 => Float64[], 2 => Float64[])
        log_resp = zeros(size(data, 2), 2)
        for (j, i) in enumerate(idxs)
            for c in 1:2
                if opt.selftrain > 0
                    # After self-training, score by Mahalanobis distance instead
                    # of the GMM responsibility.
                    L = cholesky(Symmetric(covs[c]) + 1e-8 * I).L
                    z = L \ (data[:, j] .- means[:, c])
                    log_resp[j, c] = log(weights[c]) - 0.5 * dot(z, z)
                else
                    log_resp[j, c] = log(weights[c]) + _gmm_log_density(data[:, j], means[:, c], covs[c])
                end
            end
            m = maximum(log_resp[j, :])
            resp = exp.(log_resp[j, :] .- m)
            resp ./= sum(resp)
            assigned = argmax(resp)
            push!(cluster_amp[assigned], records[i].amplitude)
        end
        mean_amp = Dict(c => mean(vals) for (c, vals) in cluster_amp if !isempty(vals))
        length(mean_amp) == 2 || continue
        high_cluster = first(sort(collect(keys(mean_amp)); by=c -> mean_amp[c], rev=true))
        for (j, i) in enumerate(idxs)
            m = maximum(log_resp[j, :])
            resp = exp.(log_resp[j, :] .- m)
            resp ./= sum(resp)
            label = argmax(resp) == high_cluster ? 1.0 : 0.0
            votes[i] += label
            counts[i] += 1
        end
    end

    for i in eachindex(records)
        counts[i] > 0 && (probs[i] = votes[i] / counts[i])
    end
    return probs
end

function _write_predictions(path::String, records::Vector{LobeRecord}, probs::Vector{Float64}, views_used::Vector{Int})
    _ensure_parent(path)
    islink(path) && error("Refusing to overwrite symlink: $path")
    n_forced = n_uncertain = 0
    open(path, "w") do io
        println(io, join(["file", "lobe", "predicted", "confidence", "amplitude",
                          "probability_1", "views_used", "invalid_reason"], '\t'))
        for (i, rec) in enumerate(records)
            p = probs[i]
            if isfinite(p) && views_used[i] > 0
                pred = p >= 0.5 ? 1 : 0
                conf = max(p, 1 - p)
                n_forced += 1
                println(io, join([rec.file, rec.lobe, pred, @sprintf("%.8f", conf),
                                  @sprintf("%.10g", rec.amplitude), @sprintf("%.8f", p),
                                  views_used[i], "ok"], '\t'))
            else
                n_uncertain += 1
                println(io, join([rec.file, rec.lobe, "?", "0.00000000",
                                  @sprintf("%.10g", rec.amplitude), "NA", views_used[i],
                                  "no_valid_view"], '\t'))
            end
        end
    end
    return n_forced, n_uncertain
end

function main(args=ARGS)
    opt = _parse_cli(args)
    records = _load_records(opt.features)
    _merge_split!(records, opt.split_features)
    _merge_patches!(records, opt.patches)
    views = isempty(opt.view_specs) ? _default_views(records) : opt.view_specs

    p_sum = zeros(Float64, length(records))
    p_count = zeros(Int, length(records))
    for (name, features) in views
        missing = [fn for fn in features if !(fn in _available_features(records))]
        isempty(missing) || error("View $name uses unavailable features: $(join(missing, ", "))")
        probs = _view_probability(records, features, opt)
        usable = count(isfinite, probs)
        @printf("view %-24s rows=%d/%d features=%s\n", name, usable, length(records), join(features, ","))
        for i in eachindex(records)
            if isfinite(probs[i])
                p_sum[i] += probs[i]
                p_count[i] += 1
            end
        end
    end

    probs = [p_count[i] > 0 ? p_sum[i] / p_count[i] : NaN for i in eachindex(records)]
    n_forced, n_uncertain = _write_predictions(opt.out_tsv, records, probs, p_count)
    files = length(unique(rec.file for rec in records))
    println("\nLabel-free unit predictions")
    println("  features:   ", opt.features)
    println("  split:      ", isempty(opt.split_features) ? "none" : opt.split_features)
    println("  patches:    ", isempty(opt.patches) ? "none" : opt.patches)
    println("  out:        ", opt.out_tsv)
    println("  files:      ", files)
    println("  rows:       ", length(records))
    println("  predicted:  ", n_forced)
    println("  uncertain:  ", n_uncertain)
    println("  views:      ", join(first.(views), ", "))
end

main()

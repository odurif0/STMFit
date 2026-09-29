#!/usr/bin/env julia

# Label-free two-component Gaussian mixture per feature view -> 0/1/? per lobe.
#
# Prediction-only: no benchmark truth, expected sequence or composition prior.
# Per view: per-scan standardization (+ pairwise interactions), EM with a full
# covariance (k-means start, ridge-regularized), hard Mahalanobis self-training,
# Mahalanobis scoring, and physical naming (the component with the higher mean
# raw amplitude is class 1, GlcNAc). Seeds are combined by a hard vote. With
# assignment_training_scans = "corroborated_counts", the mixture is learned only
# on scans whose count agrees with the repeated-scan molecule consensus; every
# scan is still standardized and assigned.

using Clustering
using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: load_config, load_training_scans, read_corroborated_scans

const DEFAULT_OUT = "results/unit_assignment/labelfree_gmm_predictions.tsv"

struct Options
    features::String
    out_tsv::String
    view_specs::Vector{Pair{String,Vector{String}}}
    first_seed::Int
    n_seeds::Int
    interactions::Bool
    selftrain::Int
    covariance_ridge::Float64
    scale_fallback::Float64
    maxiter::Int
    tolerance::Float64
    min_mass::Float64
    training_scans::Union{Nothing,Set{String}}
end

mutable struct LobeRecord
    file::String
    lobe::Int
    amplitude::Float64
    features::Dict{String,Float64}
end

_arg_value(args, i::Int, flag::String) = (i < length(args) || error("$flag requires a value"); args[i+1])

function _parse_view(s::AbstractString)
    parts = split(String(s), '='; limit=2)
    length(parts) == 2 || error("--view must be NAME=feature1,feature2,...")
    name = strip(parts[1]); isempty(name) && error("--view name is empty")
    features = [strip(f) for f in split(parts[2], ',') if !isempty(strip(f))]
    isempty(features) && error("--view $name has no features")
    return String(name) => String.(features)
end

"Reject any configuration outside the single supported learning path."
function _require_supported(cfg)
    m, s, p = cfg["model"], cfg["selection"], cfg["preprocessing"]
    expected = ((m, "gmm_learning_family", "gaussian"), (m, "gmm_covariance_scope", "final_only"),
                (m, "gmm_hard_assignment", "mahalanobis"), (m, "gmm_covariance_structure", "full"),
                (m, "gmm_final_score", "mahalanobis"), (m, "gmm_final_covariance", "ridge"),
                (s, "gmm_cluster_naming", "raw_amplitude"), (s, "gmm_resampling", "none"),
                (s, "gmm_seed_aggregation", "hard_vote"), (s, "gmm_training_weighting", "equal_lobes"),
                (s, "assignment_training_support", "all_admissible"),
                (p, "gmm_feature_normalization", "mean_sample_std"))
    for (section, key, value) in expected
        get(section, key, nothing) == value || error("Unsupported GMM setting $key=$(repr(get(section, key, nothing))); only $value is defined")
    end
end

function _parse_cli(args)
    out_tsv = DEFAULT_OUT; features = ""; config = ""; training_scans = ""
    view_specs = Pair{String,Vector{String}}[]
    first_seed = 0; n_seeds = 10; interactions = false; selftrain = 0
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"; features = _arg_value(args, i, arg); i += 2
        elseif arg == "--out"; out_tsv = _arg_value(args, i, arg); i += 2
        elseif arg == "--config"; config = _arg_value(args, i, arg); i += 2
        elseif arg == "--training-scans"; training_scans = _arg_value(args, i, arg); i += 2
        elseif arg == "--view"; push!(view_specs, _parse_view(_arg_value(args, i, arg))); i += 2
        elseif arg == "--first-seed"; first_seed = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif arg == "--seeds"; n_seeds = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif arg == "--selftrain"; selftrain = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif arg == "--interactions"; interactions = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/build_labelfree_gmm_predictions.jl --features TSV --config TOML
                   --view NAME=f1,f2,... [--view ...] [--seeds N] [--first-seed S] [--selftrain K]
                   [--interactions] [--training-scans consensus_summary.tsv] [--out TSV]

            Output columns: file, lobe, predicted, confidence, amplitude, probability_1,
            views_used, invalid_reason. Label-free: no truth, sequence or composition input.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(features) || error("Feature TSV not found: $features")
    isfile(config) || error("Config not found: $config")
    isempty(view_specs) && error("At least one --view is required")
    n_seeds > 0 && selftrain >= 0 || error("--seeds must be positive and --selftrain non-negative")
    cfg = load_config(config)
    _require_supported(cfg)
    m = cfg["model"]; pre = cfg["preprocessing"]
    corroborated = load_training_scans(cfg) == "corroborated_counts"
    corroborated == !isempty(training_scans) || error("--training-scans is required only with corroborated_counts")
    scans = isempty(training_scans) ? nothing : read_corroborated_scans(training_scans)
    return Options(features, out_tsv, view_specs, first_seed, n_seeds, interactions, selftrain,
        Float64(m["gmm_covariance_ridge"]), Float64(pre["gmm_scale_fallback"]),
        Int(m["gmm_learning_maxiter"]), Float64(m["gmm_learning_tolerance"]), Float64(m["gmm_learning_min_mass"]), scans)
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
        file = basename(strip(row["file"])); isempty(file) && error("Empty file value in $path")
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe !== nothing && lobe >= 1 || error("Invalid lobe index in $path: $(row["lobe"])")
        (file, lobe) in seen && error("Duplicate row for $file lobe $lobe: $path")
        push!(seen, (file, lobe))
        amp = _parse_f(row["amplitude"])
        isfinite(amp) || error("Non-finite amplitude for $file lobe $lobe")
        push!(records, LobeRecord(file, lobe, amp, _record_features(row)))
    end
    sort!(records, by=r -> (r.file, r.lobe))
    return records
end

"Per-scan mean/sample-std standardization, then pairwise interactions."
function _standardized_matrix(records::Vector{LobeRecord}, features::Vector{String}; interactions::Bool, scale_fallback::Float64)
    n = length(records); p = length(features)
    raw = fill(NaN, n, p); valid = trues(n)
    for (i, rec) in enumerate(records), (j, fn) in enumerate(features)
        v = get(rec.features, fn, NaN)
        raw[i, j] = v
        isfinite(v) || (valid[i] = false)
    end
    z = fill(NaN, n, p)
    for file in sort(unique(rec.file for rec in records))
        idxs = findall(i -> records[i].file == file, 1:n)
        for j in 1:p
            good = filter(isfinite, raw[idxs, j])
            isempty(good) && continue
            μ = mean(good); σ = std(good); σ = σ > 0 ? σ : scale_fallback
            for i in idxs
                z[i, j] = isfinite(raw[i, j]) ? (raw[i, j] - μ) / σ : NaN
            end
        end
    end
    if interactions && p > 1
        extra = Matrix{Float64}(undef, n, p * (p - 1) ÷ 2)
        col = 1
        for a in 1:(p-1), b in (a+1):p
            extra[:, col] = z[:, a] .* z[:, b]; col += 1
        end
        z = hcat(z, extra)
    end
    for i in 1:n
        all(isfinite, z[i, :]) || (valid[i] = false)
    end
    return z, valid
end

function _gmm_log_density(x::AbstractVector, mu::AbstractVector, Sigma::AbstractMatrix)
    d = x - mu
    L = cholesky(Symmetric(Sigma) + 1e-8 * I).L
    z = L \ d
    return -0.5 * (length(x) * log(2π) + 2 * sum(log.(diag(L))) + dot(z, z))
end

"Component score used for the final hard assignment."
function _final_component_score(x, mu, covariance, weight, selftrain::Int)
    if selftrain > 0
        L = cholesky(Symmetric(covariance) + 1e-8 * I).L
        z = L \ (x .- mu)
        return log(weight) - 0.5 * dot(z, z)
    end
    return log(weight) + _gmm_log_density(x, mu, covariance)
end

"Two-component EM with full, ridge-regularized covariances from a k-means start."
function _gmm_fit(X::Matrix{Float64}, seed::Int; ridge::Float64, max_iter::Int, tol::Float64, min_mass::Float64)
    p, n = size(X)
    k = 2
    km = kmeans(X, k; maxiter=100, rng=MersenneTwister(seed), display=:none)
    means = zeros(p, k)
    covs = [Matrix{Float64}(I, p, p) for _ in 1:k]
    weights = fill(1.0 / k, k)
    for c in 1:k
        members = findall(km.assignments .== c)
        isempty(members) && continue
        means[:, c] = sum(X[:, members]; dims=2) / length(members)
        centered = X[:, members] .- means[:, c]
        covs[c] = centered * centered' / length(members)
        covs[c] += ridge * I
    end
    prev_ll = -Inf
    for iter in 1:max_iter
        log_resp = zeros(n, k)
        for j in 1:k, i in 1:n
            log_resp[i, j] = log(weights[j]) + _gmm_log_density(X[:, i], means[:, j], covs[j])
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
        for j in 1:k
            nk = sum(resp[:, j])
            nk > min_mass || continue
            weights[j] = nk / n
            means[:, j] = X * resp[:, j] / nk
            centered = X .- means[:, j]
            covs[j] = (centered * Diagonal(resp[:, j]) * centered') / nk
            covs[j] += ridge * I
        end
    end
    return means, covs, weights
end

"Hard self-training: reassign by Mahalanobis distance, re-estimate `iters` times."
function _mahalanobis_self_train(X::Matrix{Float64}, means::Matrix{Float64}, covs::Vector{Matrix{Float64}},
                                 weights::Vector{Float64}; iters::Int, ridge::Float64)
    p, n = size(X); k = size(means, 2)
    for iteration in 1:iters
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
            covs[j] = iteration == iters ? centered * transpose(centered) / nk + ridge * I :
                centered * centered' / nk + ridge * I
        end
    end
    return means, covs, weights
end

function _view_probability(records::Vector{LobeRecord}, features::Vector{String}, opt::Options)
    X, valid = _standardized_matrix(records, features; interactions=opt.interactions, scale_fallback=opt.scale_fallback)
    eligible = opt.training_scans === nothing ? trues(length(records)) : [r.file in opt.training_scans for r in records]
    idxs = findall(valid .& eligible)
    score_idxs = findall(valid)
    length(idxs) >= 2 || return fill(NaN, length(records))
    votes = zeros(Float64, length(records)); counts = zeros(Int, length(records))
    data = permutedims(X[idxs, :]); scoring_data = permutedims(X[score_idxs, :])
    for seed in opt.first_seed:(opt.first_seed + opt.n_seeds - 1)
        means, covs, weights = _gmm_fit(data, seed; ridge=opt.covariance_ridge, max_iter=opt.maxiter,
                                        tol=opt.tolerance, min_mass=opt.min_mass)
        if opt.selftrain > 0
            means, covs, weights = _mahalanobis_self_train(data, means, covs, weights;
                iters=opt.selftrain, ridge=opt.covariance_ridge)
        end
        assignments = Int[]
        for j in eachindex(score_idxs)
            scores = [_final_component_score(scoring_data[:, j], means[:, c], covs[c], weights[c], opt.selftrain) for c in 1:2]
            push!(assignments, argmax(exp.(scores .- maximum(scores))))
        end
        # Physical naming on training rows only: higher mean raw amplitude is class 1.
        values = Dict(1 => Float64[], 2 => Float64[])
        for (j, i) in enumerate(score_idxs)
            eligible[i] && push!(values[assignments[j]], records[i].amplitude)
        end
        mean_amp = Dict(c => mean(v) for (c, v) in values if !isempty(v))
        length(mean_amp) == 2 || continue
        high = first(sort(collect(keys(mean_amp)); by=c -> mean_amp[c], rev=true))
        println("GMM seed ", seed, " component_weights=", join(weights, ','), " high_cluster=", high)
        for (j, i) in enumerate(score_idxs)
            votes[i] += assignments[j] == high ? 1.0 : 0.0
            counts[i] += 1
        end
    end
    return [counts[i] > 0 ? votes[i] / counts[i] : NaN for i in eachindex(records)]
end

function _write_predictions(path::String, records::Vector{LobeRecord}, probs::Vector{Float64}, views_used::Vector{Int})
    _ensure_parent(path)
    islink(path) && error("Refusing to overwrite symlink: $path")
    open(path, "w") do io
        println(io, join(["file", "lobe", "predicted", "confidence", "amplitude",
                          "probability_1", "views_used", "invalid_reason"], '\t'))
        for (i, rec) in enumerate(records)
            p = probs[i]
            if isfinite(p) && views_used[i] > 0
                println(io, join([rec.file, rec.lobe, p >= 0.5 ? 1 : 0, @sprintf("%.8f", max(p, 1 - p)),
                                  @sprintf("%.10g", rec.amplitude), @sprintf("%.8f", p), views_used[i], "ok"], '\t'))
            else
                println(io, join([rec.file, rec.lobe, "?", "0.00000000", @sprintf("%.10g", rec.amplitude),
                                  "NA", views_used[i], "no_valid_view"], '\t'))
            end
        end
    end
end

function main(args=ARGS)
    opt = _parse_cli(args)
    records = _load_records(opt.features)
    available = reduce(union, (keys(r.features) for r in records); init=Set{String}())
    p_sum = zeros(Float64, length(records)); p_count = zeros(Int, length(records))
    for (name, features) in opt.view_specs
        missing = [fn for fn in features if !(fn in available)]
        isempty(missing) || error("View $name uses unavailable features: $(join(missing, ", "))")
        probs = _view_probability(records, features, opt)
        @printf("view %-24s rows=%d/%d features=%s\n", name, count(isfinite, probs), length(records), join(features, ","))
        for i in eachindex(records)
            isfinite(probs[i]) && (p_sum[i] += probs[i]; p_count[i] += 1)
        end
    end
    probs = [p_count[i] > 0 ? p_sum[i] / p_count[i] : NaN for i in eachindex(records)]
    _write_predictions(opt.out_tsv, records, probs, p_count)
    println("Label-free GMM predictions: ", length(records), " rows, ",
            length(unique(r.file for r in records)), " scans -> ", opt.out_tsv)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

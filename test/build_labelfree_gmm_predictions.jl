#!/usr/bin/env julia

# Build label-free 0/1/? predictions via explicitly configured two-component learning.
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
using StatsBase: Weights, sample
using TOML
using SHA

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _ensure_parent, _read_tsv
include(joinpath(@__DIR__, "lib", "assignment_covariance.jl"))
using .AssignmentCovariance
include(joinpath(@__DIR__, "lib", "assignment_mixtures.jl"))
using .AssignmentMixtures
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: load_training_policy, load_training_mask, validate_training_mask,
    load_gmm_normalization, load_gmm_weighting, load_gmm_seed_aggregation, load_gmm_resampling,
    load_gmm_covariance_structure, load_gmm_cluster_naming, load_gmm_learning, load_gmm_covariance_scope, write_table

const DEFAULT_FEATURES = "results/unit_separability/lobe_features_selectedN_primary_local.tsv"
const DEFAULT_SPLIT = ""
const DEFAULT_PATCHES = ""
const DEFAULT_OUT = "results/unit_assignment/labelfree_gmm_predictions.tsv"
const DEFAULT_CONFIG = joinpath(dirname(@__DIR__), "config", "unit_assignment_reconstructed.toml")

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
    covariance_mode::String
    covariance_scope::String
    covariance_ridge::Float64
    final_score::String
    training_policy::String
    training_support::String
    normalization::String
    scale_fallback::Float64
    training_weighting::String
    seed_aggregation::String
    resampling::String
    bootstrap_replicates::Int
    bootstrap_seed::Int
    bootstrap_audit::String
    covariance_structure::String
    cluster_naming::String
    learning::NamedTuple
end

function _load_final_score(config::AbstractDict)
    mode = get(get(config, "model", Dict()), "gmm_final_score", nothing)
    mode in ("mahalanobis", "gaussian_density", "student_density") ||
        throw(ArgumentError("explicit gmm_final_score must be mahalanobis, gaussian_density or student_density"))
    return String(mode)
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
    config = DEFAULT_CONFIG
    training_support = ""
    bootstrap_audit = ""

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--bootstrap-audit"
            bootstrap_audit = _arg_value(args, i, arg); i += 2
        elseif arg == "--training-support"
            training_support = _arg_value(args, i, arg); i += 2
        elseif arg == "--config"
            config = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--config=")
            config = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--features"
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
              --config PATH          Explicit covariance structure, ridge and score
                                     [config/unit_assignment_reconstructed.toml].
                                     gmm_learning_family=gaussian preserves legacy
                                     learning. factor_analyzer and student_t change
                                     EM and the hard parameter updates. Free masses
                                     remain unconstrained. gmm_hard_assignment and
                                     gmm_final_score both set to student_density
                                     retain Student density through hard assignment
                                     and final scoring; otherwise the learning-only
                                     alternatives keep Mahalanobis decisions.
                                     Their rank/df and iteration settings are fixed
                                     in config; no automatic family selection.
                                     They require selftrain >= 1 and the uncombined
                                     equal-lobe, raw-naming, hard-vote policies.
                                     With gmm_covariance_scope=final_only,
                                     Ledoit-Wolf acts only on the last hard
                                     self-training covariance. all_updates also
                                     regularizes initialization, every Gaussian EM
                                     M-step and earlier hard updates. It requires
                                     --selftrain >= 1. gaussian_density adds
                                     the covariance-volume term only to final
                                     scoring after self-training (>= 1).
                                     Neither final-score option changes vote weights.
                                     gmm_covariance_structure=full preserves separate
                                     component covariances; tied pools within-group
                                     scatter at initialization, every EM M-step and
                                     every hard update. Means and masses stay free.
                                     Tied mode requires the unresampled ridge,
                                     hard-vote, equal-lobe support-control policies.
                                     gmm_cluster_naming=raw_amplitude preserves the
                                     global raw-amplitude anchor. within_scan_z uses
                                     training-only within-scan amplitude mean/sample
                                     std (same explicit scale fallback) to name the
                                     unchanged groups; no composition is imposed.
                                     Also declares per-file GMM feature scaling:
                                     mean_sample_std or median_iqr (Type-7 Q75-Q25),
                                     and its degenerate-scale fallback. Scaling
                                     uses finite observations before interactions;
                                     it neither clips nor imputes observations.
                                     gmm_training_weighting is equal_lobes or
                                     equal_scans: the latter weights each usable
                                     training row inversely by its scan's usable
                                     row count, in initialization, EM, hard updates
                                     and amplitude naming. Mixture weights stay free.
                                     gmm_seed_aggregation is hard_vote (legacy) or
                                     mean_membership: average the normalized scores
                                     of the physically named high-amplitude group
                                     across seeds, without hardening them to 0/1.
                                     These are not calibrated chemical probabilities;
                                     fitting, scores and group naming stay unchanged.
                                     gmm_resampling=whole_scans draws whole usable
                                     scans with replacement, once per configured
                                     bootstrap replicate. Each uses the same seed
                                     range and hard votes. Per-scan scaling is fixed;
                                     learning and amplitude naming use sampled rows
                                     with their multiplicities. All valid rows are
                                     predicted, including out-of-bag scans. This is
                                     bagging, not held-out validation or calibration.
              --bootstrap-audit PATH Required new TSV for whole_scans; records all
                                     scan multiplicities and valid seed counts.
              --split-features PATH  Optional split-width feature TSV; adds split_log_skew
              --training-support PATH Observed pixel-count TSV, required only for
                                     complete_patches training. Normalization,
                                     fits and amplitude naming use complete rows;
                                     admissible partial rows are prediction-only.
              --patches PATH         Optional backward patch TSV; adds bwd_neg_* descriptors
              --out PATH             Output prediction TSV [$(DEFAULT_OUT)]
              --view NAME=LIST       Feature view to cluster. Repeatable. If omitted,
                                     sensible label-free defaults are chosen from
                                     available columns.
              --first-seed INT       First k-means seed [0]
              --seeds INT            Number of k-means seeds [20]
              --interactions         Append pairwise products within each view after
                                     configured per-file scaling.
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
    covariance = load_covariance_config(config)
    covariance.mode == "ledoit_wolf" && selftrain < 1 &&
        error("Final Ledoit-Wolf covariance requires --selftrain >= 1")
    cfg = TOML.parsefile(config)
    covariance_structure = load_gmm_covariance_structure(cfg)
    cluster_naming = load_gmm_cluster_naming(cfg)
    learning = load_gmm_learning(cfg)
    covariance_scope = load_gmm_covariance_scope(cfg)
    learning.family != "gaussian" && selftrain < 1 &&
        error("Alternative learning requires --selftrain >= 1 for the configured hard updates")
    normalization = load_gmm_normalization(cfg)
    training_weighting = load_gmm_weighting(cfg)
    seed_aggregation = load_gmm_seed_aggregation(cfg)
    bootstrap = load_gmm_resampling(cfg)
    (bootstrap.mode == "whole_scans") == !isempty(bootstrap_audit) ||
        error("--bootstrap-audit is required only with whole_scans")
    if !isempty(bootstrap_audit)
        (ispath(bootstrap_audit) || islink(bootstrap_audit)) && error("Bootstrap audit already exists")
        abspath(bootstrap_audit) != abspath(out_tsv) || error("Bootstrap audit must differ from prediction output")
    end
    final_score = _load_final_score(cfg)
    training_policy = load_training_policy(cfg)
    (training_policy == "complete_patches") == !isempty(training_support) ||
        error("--training-support is required only with complete_patches")
    final_score == "gaussian_density" && selftrain < 1 &&
        error("Explicit final Gaussian score requires --selftrain >= 1")
    return Options(features, split_features, patches, out_tsv, view_specs,
                   first_seed, n_seeds, interactions, selftrain, covariance.mode, covariance_scope, covariance.ridge, final_score,
                   training_policy, training_support, normalization.mode, normalization.scale_fallback,
                   training_weighting, seed_aggregation, bootstrap.mode, bootstrap.replicates, bootstrap.seed,
                   bootstrap_audit, covariance_structure, cluster_naming, learning)
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

function _standardized_matrix(records::Vector{LobeRecord}, features::Vector{String}; interactions::Bool=false,
                              training_mask=nothing, normalization::String, scale_fallback::Real)
    normalization in ("mean_sample_std", "median_iqr") || throw(ArgumentError("unknown GMM normalization"))
    !(scale_fallback isa Bool) && isfinite(scale_fallback) && scale_fallback > 0 ||
        throw(ArgumentError("scale fallback must be positive and finite"))
    n = length(records)
    training_mask === nothing || (training_mask isa AbstractVector{Bool} && length(training_mask) == n) ||
        throw(ArgumentError("normalization training mask must match rows"))
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
        fit_idxs = training_mask === nothing ? idxs : filter(i -> training_mask[i], idxs)
        for j in 1:p
            vals = raw[fit_idxs, j]
            good = filter(isfinite, vals)
            if isempty(good)
                z[idxs, j] .= NaN
            else
                if normalization == "mean_sample_std"
                    # Preserve the original operations exactly for the control.
                    μ = mean(good)
                    σ = std(good)
                    σ = σ > 0 ? σ : scale_fallback
                else
                    μ = median(good)
                    q = quantile(good, [0.25, 0.75]; alpha=1, beta=1)
                    σ = q[2] - q[1]
                    σ = isfinite(σ) && σ > 0 ? σ : scale_fallback
                end
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

"Final hard-cluster scoring only: same fitted means, covariances and weights."
function _final_component_score(x, mu, covariance, weight, opt::Options)
    if opt.final_score == "student_density"
        return student_component_score(x, mu, covariance, weight;
            df=opt.learning.student_df, guard=opt.learning.cholesky_guard)
    end
    if opt.selftrain > 0 && opt.final_score == "mahalanobis"
        L = cholesky(Symmetric(covariance) + 1e-8 * I).L
        z = L \ (x .- mu)
        return log(weight) - 0.5 * dot(z, z)
    end
    return log(weight) + _gmm_log_density(x, mu, covariance)
end

"Mean-one weights; each represented scan has total n / number_of_scans."
function _observation_weights(records, idxs, mode::String)
    mode == "equal_lobes" && return nothing # Preserve legacy floating-point operations.
    mode == "equal_scans" || throw(ArgumentError("unknown GMM training weighting"))
    isempty(idxs) && return Float64[]
    counts = Dict{String,Int}()
    for i in idxs
        counts[records[i].file] = get(counts, records[i].file, 0) + 1
    end
    mass = length(idxs) / length(counts)
    return [mass / counts[records[i].file] for i in idxs]
end

function _check_observation_weights(w, n)
    w === nothing && return nothing
    w isa AbstractVector{<:Real} && length(w) == n &&
        all(x -> !(x isa Bool) && isfinite(x) && x > 0, w) && isfinite(sum(w)) ||
        throw(ArgumentError("observation weights must be positive finite values matching training rows"))
    return Float64.(w)
end

"Weighted k-means++ for the two GMM components, without a composition constraint."
function _weighted_seeds(X, w, rng)
    n = size(X, 2)
    first_seed = sample(rng, 1:n, Weights(w))
    costs = [w[i] * sum(abs2, X[:, i] .- X[:, first_seed]) for i in 1:n]
    if sum(costs) == 0
        costs = copy(w)
        costs[first_seed] = 0.0
    end
    return [first_seed, sample(rng, 1:n, Weights(costs))]
end

function _cluster_amplitude_means(records, score_idxs, assignments, eligible, observation_weights;
                                  mode::String="raw_amplitude", scale_fallback=nothing)
    mode in ("raw_amplitude", "within_scan_z") || throw(ArgumentError("unknown GMM cluster naming"))
    normalized = Dict{Int,Float64}()
    if mode == "within_scan_z"
        observation_weights === nothing || throw(ArgumentError("relative naming requires equal-lobe weights"))
        scale_fallback isa Real && !(scale_fallback isa Bool) && isfinite(scale_fallback) && scale_fallback > 0 ||
            throw(ArgumentError("relative naming requires an explicit positive finite scale fallback"))
        training = [i for i in score_idxs if eligible[i]]
        length(unique(training)) == length(training) || throw(ArgumentError("relative naming requires unresampled rows"))
        for file in sort(unique(records[i].file for i in training))
            members = [i for i in training if records[i].file == file]
            amplitudes = [records[i].amplitude for i in members]
            all(isfinite, amplitudes) || throw(ArgumentError("nonfinite training naming amplitude"))
            center, scale = mean(amplitudes), std(amplitudes)
            scale = isfinite(scale) && scale > 0 ? scale : scale_fallback
            for i in members
                normalized[i] = (records[i].amplitude - center) / scale
            end
        end
    end
    values = Dict(1 => Float64[], 2 => Float64[])
    masses = Dict(1 => Float64[], 2 => Float64[])
    for (j, i) in enumerate(score_idxs)
        eligible[i] || continue
        c = assignments[j]
        push!(values[c], mode == "raw_amplitude" ? records[i].amplitude : normalized[i])
        observation_weights === nothing || push!(masses[c], observation_weights[i])
    end
    return Dict(c => (observation_weights === nothing ? mean(vals) : dot(vals, masses[c]) / sum(masses[c]))
                for (c, vals) in values if !isempty(vals))
end

"Shared within-component ML covariance; responsibilities have one unit of mass per row."
function _tied_covariance(X::AbstractMatrix, means::AbstractMatrix, resp::AbstractMatrix; ridge::Real)
    p, n = size(X)
    n > 0 && p > 0 && size(means, 1) == p && size(resp) == (n, size(means, 2)) ||
        throw(ArgumentError("invalid tied covariance dimensions"))
    all(isfinite, X) && all(isfinite, means) && all(isfinite, resp) && all(>=(0), resp) ||
        throw(ArgumentError("invalid tied covariance observations or responsibilities"))
    all(s -> isapprox(s, 1; atol=1e-12, rtol=1e-12), vec(sum(resp; dims=2))) ||
        throw(ArgumentError("tied responsibilities must sum to one per row"))
    !(ridge isa Bool) && isfinite(ridge) && ridge > 0 || throw(ArgumentError("positive finite ridge required"))
    scatter = zeros(p, p)
    for c in axes(means, 2)
        centered = X .- means[:, c]
        scatter .+= centered * Diagonal(resp[:, c]) * centered'
    end
    sample = scatter / n
    all(isfinite, sample) || throw(ArgumentError("tied covariance overflow"))
    return (sample=sample, covariance=sample + ridge * I)
end

# EM fit of a 2-component GMM, k-means initialized. The legacy full route
# retains its exact operations; production always passes the explicit config mode.
# Returns (means, covariances, weights).
function _gmm_fit(X::Matrix{Float64}, seed::Int; ridge::Float64, max_iter::Int=200, tol::Float64=1e-6,
                  observation_weights=nothing, covariance_structure::String="full",
                  covariance_scope::String="final_only", min_mass::Float64=1e-12, trace=nothing)
    covariance_scope in ("final_only", "all_updates") || throw(ArgumentError("unknown covariance scope"))
    regularized = covariance_scope == "all_updates"
    isfinite(min_mass) && min_mass > 0 || throw(ArgumentError("positive finite component mass floor required"))
    regularized && (covariance_structure != "full" || observation_weights !== nothing) &&
        throw(ArgumentError("all-updates shrinkage requires full covariance and equal-lobe training"))
    covariance_structure in ("full", "tied") || throw(ArgumentError("unknown covariance structure"))
    covariance_structure == "tied" && observation_weights !== nothing &&
        throw(ArgumentError("tied covariance requires equal-lobe training"))
    p, n = size(X)
    w = _check_observation_weights(observation_weights, n)
    total = w === nothing ? n : sum(w)
    rng = MersenneTwister(seed)
    k = 2
    km = w === nothing ? kmeans(X, k; maxiter=100, rng=rng, display=:none) :
        kmeans(X, k; weights=w, init=_weighted_seeds(X, w, rng), maxiter=100, rng=rng, display=:none)
    means = zeros(p, k)
    covs = [Matrix{Float64}(I, p, p) for _ in 1:k]
    weights = fill(1.0 / k, k)
    shrinkages = zeros(k)
    all_shrinkages = Float64[]
    for c in 1:k
        members = findall(km.assignments .== c)
        isempty(members) && (regularized ? throw(MixtureFitError("empty_initial_cluster")) : continue)
        mass = w === nothing ? length(members) : sum(w[members])
        means[:, c] = w === nothing ? sum(X[:, members]; dims=2) / length(members) : X[:, members] * w[members] / mass
        centered = X[:, members] .- means[:, c]
        covs[c] = w === nothing ? centered * centered' / length(members) :
            centered * Diagonal(w[members]) * centered' / mass
        covs[c] += ridge * I
        if regularized
            estimate = final_covariance(centered; mode="ledoit_wolf", ridge)
            covs[c], shrinkages[c] = estimate.covariance, estimate.shrinkage
        end
    end
    if covariance_structure == "tied"
        resp = Float64.([a == c for a in km.assignments, c in 1:k])
        shared = _tied_covariance(X, means, resp; ridge).covariance
        for c in 1:k
            covs[c] = copy(shared)
        end
    end
    if regularized
        append!(all_shrinkages, shrinkages)
        trace === nothing || push!(trace, (stage="initial", iteration=0, means=copy(means), covs=deepcopy(covs),
            weights=copy(weights), responsibilities=Float64.([a == c for a in km.assignments, c in 1:k]),
            shrinkages=copy(shrinkages)))
    end
    prev_ll = -Inf
    converged, updates = false, 0
    ll_decreases, max_ll_drop = 0, 0.0
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
            term = m + log(sum(exp.(log_resp[i, :] .- m)))
            ll += w === nothing ? term : w[i] * term
        end
        if regularized && isfinite(prev_ll) && ll < prev_ll
            ll_decreases += 1
            max_ll_drop = max(max_ll_drop, prev_ll - ll)
        end
        if abs(ll - prev_ll) < tol * (1 + abs(prev_ll))
            converged = true
            break
        end
        prev_ll = ll
        for i in 1:n
            m = maximum(log_resp[i, :])
            log_resp[i, :] .-= m
            log_resp[i, :] .-= log(sum(exp.(log_resp[i, :])))
        end
        resp = exp.(log_resp)
        w === nothing || (resp .*= w)
        # M-step
        for j in 1:k
            nk = sum(resp[:, j])
            nk > min_mass || (regularized ? throw(MixtureFitError("collapsed_component")) : continue)
            weights[j] = nk / total
            means[:, j] = X * resp[:, j] / nk
            centered = X .- means[:, j]
            covs[j] = (centered * Diagonal(resp[:, j]) * centered') / nk
            covs[j] += ridge * I
            if regularized
                estimate = responsibility_covariance(centered, resp[:, j]; ridge)
                covs[j], shrinkages[j] = estimate.covariance, estimate.shrinkage
            end
        end
        if covariance_structure == "tied"
            shared = _tied_covariance(X, means, resp; ridge).covariance
            for j in 1:k
                covs[j] = copy(shared)
            end
        end
        updates += 1
        if regularized
            append!(all_shrinkages, shrinkages)
            trace === nothing || push!(trace, (stage="em", iteration=iter, means=copy(means), covs=deepcopy(covs),
                weights=copy(weights), responsibilities=copy(resp), shrinkages=copy(shrinkages)))
        end
    end
    regularized && println("GMM shrinkage seed ", seed, " scope=all_updates em_updates=", updates,
        " converged=", converged, " shrinkage_min=", minimum(all_shrinkages),
        " shrinkage_max=", maximum(all_shrinkages), " likelihood_decreases=", ll_decreases,
        " max_likelihood_drop=", max_ll_drop)
    return means, covs, weights
end

# Hard self-training: reassign every point to the closest cluster by Mahalanobis
# distance, then re-estimate means/covariances from the hard memberships.
# Repeats `iters` times. Returns updated (means, covariances, weights).
function _mahalanobis_self_train(X::Matrix{Float64}, means::Matrix{Float64},
                                 covs::Vector{Matrix{Float64}}, weights::Vector{Float64};
                                 iters::Int=5, covariance_mode::String, ridge::Float64,
                                 diagnostics=nothing, observation_weights=nothing,
                                 covariance_structure::String="full", covariance_scope::String="final_only",
                                 min_mass::Float64=1e-12, trace=nothing)
    covariance_scope in ("final_only", "all_updates") || throw(ArgumentError("unknown covariance scope"))
    regularized = covariance_scope == "all_updates"
    isfinite(min_mass) && min_mass > 0 || throw(ArgumentError("positive finite component mass floor required"))
    regularized && (covariance_mode != "ledoit_wolf" || covariance_structure != "full" || observation_weights !== nothing) &&
        throw(ArgumentError("all-updates hard shrinkage requires full Ledoit-Wolf and equal-lobe training"))
    covariance_structure in ("full", "tied") || throw(ArgumentError("unknown covariance structure"))
    covariance_structure == "tied" && (covariance_mode != "ridge" || observation_weights !== nothing) &&
        throw(ArgumentError("tied hard updates require ridge and equal-lobe training"))
    covariance_mode in ("ridge", "ledoit_wolf") || error("Unknown final covariance mode")
    covariance_mode == "ledoit_wolf" && iters < 1 && error("Final covariance requires hard memberships")
    p, n = size(X)
    w = _check_observation_weights(observation_weights, n)
    w !== nothing && covariance_mode != "ridge" && error("Weighted hard updates require ridge covariance")
    total = w === nothing ? n : sum(w)
    k = size(means, 2)
    for iteration in 1:iters
        shrinkages = zeros(k)
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
            isempty(members) && (regularized ? throw(MixtureFitError("empty_hard_cluster")) : continue)
            nk = w === nothing ? length(members) : sum(w[members])
            regularized && nk <= min_mass && throw(MixtureFitError("collapsed_hard_component"))
            weights[j] = nk / total
            means[:, j] = w === nothing ? sum(X[:, members]; dims=2) / nk : X[:, members] * w[members] / nk
            centered = X[:, members] .- means[:, j]
            if iteration == iters || regularized
                estimate = if w === nothing
                    final_covariance(centered; mode=covariance_mode, ridge=ridge)
                else
                    covariance = centered * Diagonal(w[members]) * centered' / nk
                    (sample=covariance, covariance=covariance + ridge * I, shrinkage=0.0)
                end
                covs[j] = estimate.covariance
                shrinkages[j] = estimate.shrinkage
                if diagnostics !== nothing && covariance_structure == "full" && iteration == iters
                    push!(diagnostics, (component=j, members=copy(members), weight=weights[j],
                        mean=copy(means[:, j]), sample=estimate.sample, covariance=copy(covs[j]),
                        shrinkage=estimate.shrinkage))
                end
            else
                covs[j] = (w === nothing ? centered * centered' / nk :
                    centered * Diagonal(w[members]) * centered' / nk) + ridge * I
            end
        end
        if regularized
            println("GMM hard shrinkage iteration=", iteration, " shrinkages=", join(shrinkages, ','))
            trace === nothing || push!(trace, (stage="hard", iteration, means=copy(means), covs=deepcopy(covs),
                weights=copy(weights), responsibilities=Float64.([a == c for a in assign, c in 1:k]),
                shrinkages=copy(shrinkages)))
        end
        if covariance_structure == "tied"
            resp = Float64.([a == c for a in assign, c in 1:k])
            shared = _tied_covariance(X, means, resp; ridge)
            for j in 1:k
                covs[j] = copy(shared.covariance)
                members = findall(==(j), assign)
                if diagnostics !== nothing && iteration == iters && !isempty(members)
                    push!(diagnostics, (component=j, members=copy(members), weight=weights[j],
                        mean=copy(means[:, j]), sample=copy(shared.sample), covariance=copy(covs[j]),
                        shrinkage=0.0))
                end
            end
        end
    end
    return means, covs, weights
end

"Only the contribution to the seed average changes; resp is already normalized."
function _seed_contribution(resp, high_cluster, mode::String)
    if mode == "hard_vote"
        return argmax(resp) == high_cluster ? 1.0 : 0.0
    end
    mode == "mean_membership" || throw(ArgumentError("unknown GMM seed aggregation"))
    return resp[high_cluster]
end

"Draw S whole usable scans with replacement, preserving every usable row per draw."
function _scan_bootstrap_sample(records, idxs, seed)
    groups = Vector{Int}[]
    locations = Dict{String,Int}()
    for i in idxs
        f = records[i].file
        if !haskey(locations, f)
            locations[f] = length(groups) + 1
            push!(groups, Int[])
        end
        push!(groups[locations[f]], i)
    end
    isempty(groups) && return Int[]
    rng = MersenneTwister(seed)
    draws = rand(rng, 1:length(groups), length(groups))
    return reduce(vcat, (groups[j] for j in draws))
end

function _bootstrap_view_probability(records, X, idxs, score_idxs, opt;
                                     diagnostics=nothing, audit=nothing, view_name="")
    probabilities = fill(NaN, length(records))
    score_data = permutedims(X[score_idxs, :])
    position = zeros(Int, length(records))
    position[score_idxs] = eachindex(score_idxs)
    total = zeros(length(score_idxs))
    usable_bags = 0
    files = unique(r.file for r in records)
    usable = Dict(f => count(i -> records[i].file == f, idxs) for f in files)
    eligible = trues(length(records)) # train_idxs below already contains only eligible valid rows.
    for bag in 1:opt.bootstrap_replicates
        draw_seed = opt.bootstrap_seed + bag - 1
        train_idxs = _scan_bootstrap_sample(records, idxs, draw_seed)
        data = permutedims(X[train_idxs, :])
        votes = zeros(length(score_idxs))
        accepted = 0
        for seed in opt.first_seed:(opt.first_seed + opt.n_seeds - 1)
            means, covs, weights = _gmm_fit(data, seed; ridge=opt.covariance_ridge,
                covariance_structure=opt.covariance_structure)
            clusters = diagnostics === nothing ? nothing : []
            if opt.selftrain > 0
                means, covs, weights = _mahalanobis_self_train(data, means, covs, weights;
                    iters=opt.selftrain, covariance_mode=opt.covariance_mode, ridge=opt.covariance_ridge,
                    diagnostics=clusters, covariance_structure=opt.covariance_structure)
            end
            assignments = Int[]
            for j in axes(score_data, 2)
                scores = [_final_component_score(score_data[:,j], means[:,c], covs[c], weights[c], opt) for c in 1:2]
                resp = exp.(scores .- maximum(scores)); resp ./= sum(resp)
                push!(assignments, argmax(resp))
            end
            # Duplicate training indices count each scan draw in physical naming,
            # while amplitudes from unsampled scans cannot rename either group.
            mean_amp = _cluster_amplitude_means(records, train_idxs,
                assignments[position[train_idxs]], eligible, nothing;
                mode=opt.cluster_naming, scale_fallback=opt.scale_fallback)
            high = length(mean_amp) == 2 ? first(sort(collect(keys(mean_amp)); by=c -> mean_amp[c], rev=true)) : 0
            diagnostics === nothing || push!(diagnostics, (replicate=bag, bootstrap_seed=draw_seed,
                seed=seed, indices=copy(train_idxs), clusters=clusters, amplitude_means=copy(mean_amp),
                high_cluster=high, means=copy(means), covariances=deepcopy(covs), weights=copy(weights)))
            high == 0 && continue # Same missing-group rule as the unresampled head; no retry.
            votes .+= assignments .== high
            accepted += 1
        end
        if accepted > 0
            total .+= votes ./ accepted
            usable_bags += 1
        end
        println("bootstrap replicate ", bag, "/", opt.bootstrap_replicates,
            " draw_seed=", draw_seed, " training_rows=", length(train_idxs),
            " unique_scans=", length(unique(records[i].file for i in train_idxs)),
            " valid_seeds=", accepted, "/", opt.n_seeds)
        if audit !== nothing
            for f in files
                rows = count(i -> records[i].file == f, train_idxs)
                push!(audit, Dict("view"=>view_name, "replicate"=>string(bag),
                    "bootstrap_seed"=>string(draw_seed), "file"=>f, "usable_rows"=>string(usable[f]),
                    "multiplicity"=>string(usable[f] == 0 ? 0 : rows ÷ usable[f]),
                    "training_rows"=>string(rows), "valid_seeds"=>string(accepted), "total_seeds"=>string(opt.n_seeds)))
            end
        end
    end
    usable_bags > 0 && (probabilities[score_idxs] = total / usable_bags)
    return probabilities
end

function _view_probability(records::Vector{LobeRecord}, features::Vector{String}, opt::Options;
                           diagnostics=nothing, training_mask=nothing, audit=nothing, view_name="")
    eligible = validate_training_mask(opt.training_policy, training_mask, length(records))
    X, valid = _standardized_matrix(records, features; interactions=opt.interactions, training_mask,
        normalization=opt.normalization, scale_fallback=opt.scale_fallback)
    idxs = findall(valid .& eligible)
    score_idxs = findall(valid)
    length(idxs) >= 2 || return fill(NaN, length(records))
    if opt.resampling == "whole_scans"
        return _bootstrap_view_probability(records, X, idxs, score_idxs, opt; diagnostics, audit, view_name)
    end

    probs = fill(NaN, length(records))
    votes = zeros(Float64, length(records))
    counts = zeros(Int, length(records))
    data = permutedims(X[idxs, :])
    scoring_data = permutedims(X[score_idxs, :])
    observation_weights = _observation_weights(records, idxs, opt.training_weighting)
    naming_weights = if observation_weights === nothing
        nothing
    else
        w = zeros(length(records))
        w[idxs] = observation_weights
        w
    end

    for seed in opt.first_seed:(opt.first_seed + opt.n_seeds - 1)
        seed_diagnostics = diagnostics === nothing ? nothing : []
        update_trace = diagnostics === nothing ? nothing : []
        if opt.learning.family == "gaussian"
            try
                means, covs, weights = _gmm_fit(data, seed; ridge=opt.covariance_ridge, observation_weights,
                    covariance_structure=opt.covariance_structure, covariance_scope=opt.covariance_scope,
                    max_iter=opt.learning.maxiter, tol=opt.learning.tolerance,
                    min_mass=opt.learning.min_mass, trace=update_trace)
                if opt.selftrain > 0
                    means, covs, weights = _mahalanobis_self_train(data, means, covs, weights;
                        iters=opt.selftrain, covariance_mode=opt.covariance_mode, ridge=opt.covariance_ridge,
                        diagnostics=seed_diagnostics, observation_weights, covariance_structure=opt.covariance_structure,
                        covariance_scope=opt.covariance_scope, min_mass=opt.learning.min_mass, trace=update_trace)
                end
            catch err
                opt.covariance_scope == "all_updates" && err isa MixtureFitError || rethrow()
                println("Mixture fit seed ", seed, " family=gaussian scope=all_updates status=unavailable reason=", err.reason)
                continue
            end
        else
            initial = _gmm_fit(data, seed; ridge=opt.covariance_ridge, max_iter=0)
            result = try
                fit_mixture(data, initial, opt.learning; ridge=opt.covariance_ridge,
                    hard_iterations=opt.selftrain, trace=seed_diagnostics)
            catch err
                err isa MixtureFitError || rethrow()
                println("Mixture fit seed ", seed, " family=", opt.learning.family, " status=unavailable reason=", err.reason)
                continue
            end
            means, covs, weights = result.state.means, result.state.covs, result.state.weights
            noise_min = isempty(result.state.noise) ? NaN : minimum(minimum, result.state.noise)
            factor_error = isempty(result.state.noise) ? NaN : maximum(maximum(abs.(covs[c] -
                result.state.loadings[c] * result.state.loadings[c]' - Diagonal(result.state.noise[c]))) for c in 1:2)
            println("Mixture fit seed ", seed, " family=", opt.learning.family,
                " status=ok em_updates=", result.updates, " converged=", result.converged,
                " hard_updates=", opt.selftrain, " rank=", opt.learning.factor_rank, " df=", opt.learning.student_df,
                " hard_assignment=", opt.learning.hard_assignment,
                " noise_min=", noise_min, " factor_error=", factor_error,
                " precision_min=", result.precision_range[1], " precision_max=", result.precision_range[2],
                " log_kernel_sum=", result.log_kernel_sum)
        end
        diagnostics === nothing || push!(diagnostics,
            (seed=seed, indices=copy(idxs), clusters=seed_diagnostics,
             observation_weights=observation_weights === nothing ? nothing : copy(observation_weights), updates=update_trace))
        # Physical mapping is fitted on training members only, then frozen.
        assignments = Int[]
        log_resp = zeros(size(scoring_data, 2), 2)
        for (j, i) in enumerate(score_idxs)
            for c in 1:2
                log_resp[j, c] = _final_component_score(scoring_data[:, j], means[:, c], covs[c], weights[c], opt)
            end
            m = maximum(log_resp[j, :])
            resp = exp.(log_resp[j, :] .- m)
            resp ./= sum(resp)
            assigned = argmax(resp)
            push!(assignments, assigned)
        end
        mean_amp = _cluster_amplitude_means(records, score_idxs, assignments, eligible, naming_weights;
            mode=opt.cluster_naming, scale_fallback=opt.scale_fallback)
        println("GMM seed ", seed, " covariance_structure=", opt.covariance_structure,
            " covariance_maxdiff=", maximum(abs.(covs[1] - covs[2])),
            " component_weights=", join(weights, ','), " named=", length(mean_amp) == 2)
        high_cluster = length(mean_amp) == 2 ? first(sort(collect(keys(mean_amp)); by=c -> mean_amp[c], rev=true)) : 0
        println("GMM naming seed ", seed, " mode=", opt.cluster_naming,
            " naming_means=", join((get(mean_amp,c,NaN) for c in 1:2), ','), " high_cluster=", high_cluster,
            " fit_sha256=", bytes2hex(sha256(reinterpret(UInt8,
                vcat(vec(means), reduce(vcat, vec.(covs)), weights)))))
        high_cluster == 0 && continue
        for (j, i) in enumerate(score_idxs)
            m = maximum(log_resp[j, :])
            resp = exp.(log_resp[j, :] .- m)
            resp ./= sum(resp)
            votes[i] += _seed_contribution(resp, high_cluster, opt.seed_aggregation)
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
    training_mask = load_training_mask(opt.training_policy, opt.training_support,
                                      [(r.file, r.lobe) for r in records], "gmm")
    training_mask === nothing || println("GMM complete-patch eligibility: $(count(training_mask))/$(length(records)) rows before feature validity")
    views = isempty(opt.view_specs) ? _default_views(records) : opt.view_specs

    p_sum = zeros(Float64, length(records))
    p_count = zeros(Int, length(records))
    bootstrap_audit = Dict{String,String}[]
    for (name, features) in views
        missing = [fn for fn in features if !(fn in _available_features(records))]
        isempty(missing) || error("View $name uses unavailable features: $(join(missing, ", "))")
        probs = _view_probability(records, features, opt; training_mask,
            audit=isempty(opt.bootstrap_audit) ? nothing : bootstrap_audit, view_name=name)
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
    if !isempty(opt.bootstrap_audit)
        write_table(opt.bootstrap_audit, ["view", "replicate", "bootstrap_seed", "file", "usable_rows",
            "multiplicity", "training_rows", "valid_seeds", "total_seeds"], bootstrap_audit)
    end
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
    println("  weighting:  ", opt.training_weighting)
    println("  covariance structure: ", opt.covariance_structure)
    println("  covariance scope: ", opt.covariance_scope)
    println("  hard assignment: ", opt.learning.hard_assignment)
    println("  learning family: ", opt.learning.family)
    println("  cluster naming: ", opt.cluster_naming)
    println("  seed aggregation: ", opt.seed_aggregation)
    println("  resampling: ", opt.resampling, " replicates=", opt.bootstrap_replicates,
        " bootstrap_seed=", opt.bootstrap_seed)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

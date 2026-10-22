# Frozen unit assignment: learn the cohort-level parts of the promoted unit
# assignment (k-means centroids of 4 views, GMM of the v_cc view, physical
# naming, fusion emission rates) on a reference set of scans, then apply them
# unchanged to any scan. The production driver (run_reconstructed_chitosan.jl)
# learns the same parts transductively on its own run cohort.
#
# Features stay as in production: standardized within each scan (mean, sample
# std, fallback 1) with pairwise interactions, so a lobe's feature vector does
# not depend on which other scans are present. The learning and scoring code
# repeats build_labelfree_unit_predictions.jl (trained path) and
# build_labelfree_gmm_predictions.jl step for step; training on a run's own
# corroborated scans reproduces that run's pred_kmeans.tsv, pred_gmm.tsv and
# predictions.tsv (test/test_frozen_assignment.jl).
#
# Label-free: no truth, expected N, sequence or composition is read. Fusion of
# a frozen model uses the reference emission rates and a fixed prior pi
# (config: fusion_prior_pi), so the reference cohort's composition is not
# carried to new molecules.
#
# A fully frozen model (train_reference_model, apply_reference_model) also
# freezes the empirical Fisher feature, the only other feature learned across
# scans: one Fisher model per lobe parity is fitted on the reference lobes and
# scores the lobes of the other parity, as the production lobe-parity CV does
# on a run. With the whole run as reference it reproduces the run's
# fisher_cv.tsv.
module FrozenAssignment

using Clustering, LinearAlgebra, Printf, Random, Statistics
using ..ReconstructedUnitAssignment: read_table, lobe_table
using ..MoleculeFusion: FusionSettings, fit_latent_class, load_fusion_settings
using ..EmpiricalFisherNative: load_fisher_config, load_patches, fit_fisher, FisherFitError, flip_u_disk,
                               PatchTable, FisherOptions

export AssignmentData, FrozenModel, FrozenFisher, FisherData, load_assignment_data, train_frozen, predict_frozen,
       fuse_calls, fit_emission_rates, write_frozen_model, read_frozen_model, KMEANS_VIEWS, GMM_VIEW,
       load_fisher_data, train_fisher, fisher_feature, with_emp_fisher, reference_files, load_frozen_settings,
       train_reference_model, apply_reference_model, load_run_assignment

const BASE4 = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"]
# Views of the production driver (scripts/pipeline/run_reconstructed_chitosan.jl).
const KMEANS_VIEWS = ["v_base" => BASE4, "v_split" => vcat(BASE4, "split_log_skew"),
                      "v_comt" => vcat(BASE4, "bwd_neg_com_t"), "v_diag45" => vcat(BASE4, "bwd_neg_diag45")]
const GMM_VIEW = "v_cc" => vcat(BASE4, ["patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher"])

"Standardized view: X is p x n (one column per lobe); valid marks complete rows."
struct ViewMatrix
    X::Matrix{Float64}
    valid::BitVector
end

"Lobes of a run, sorted by (file, lobe), with every view precomputed."
struct AssignmentData
    file::Vector{String}
    lobe::Vector{Int}
    amplitude::Vector{Float64}
    kmeans::Vector{ViewMatrix}
    gmm::ViewMatrix
    gmm_raw::Matrix{Float64}       # unstandardized GMM view columns (lobe x feature)
    scale_fallback::Float64
end

"Fisher direction and origin of one lobe parity (score = dot(x - mid, w_p))."
struct FrozenFisher
    w_p::Vector{Float64}
    mid::Vector{Float64}
end

struct GMMSeed
    means::Matrix{Float64}
    covs::Vector{Matrix{Float64}}
    weights::Vector{Float64}
    high::Int          # 0: physical naming failed, the seed does not vote
end

"Cohort-level parameters of the assignment; centers[v][s] is p x 2, high the class-1 (brighter) cluster."
struct FrozenModel
    kmeans_centers::Vector{Vector{Matrix{Float64}}}
    kmeans_high::Vector{Vector{Int}}
    gmm::Vector{GMMSeed}
    selftrain::Int
    theta0::Float64
    theta1::Float64
    pi::Float64
    # fisher[1]: fitted on even lobes, scores odd lobes; fisher[2]: the reverse.
    # Empty when the Fisher feature is taken from the run (not frozen).
    fisher::Vector{Union{Nothing,FrozenFisher}}
end

_value(s) = something(tryparse(Float64, strip(String(s))), NaN)

"Negative-lobe moments of a backward residual patch (as build_labelfree_unit_predictions.jl)."
function _negative_moment(vals::Vector{Float64}, coords::Vector{Float64})
    weights = max.(-vals, 0.0)
    total = sum(weights)
    total <= eps(Float64) && return (NaN, NaN)
    com_t = 0.0; diag45 = 0.0; idx = 1
    for u in coords, t in coords
        w = weights[idx]
        com_t += w * t
        diag45 += w * ((t + u) / sqrt(2.0))
        idx += 1
    end
    return (com_t / total, diag45 / total)
end

function _view(raw::Matrix{Float64}, file::Vector{String}, scale_fallback::Float64)
    n, p = size(raw)
    valid = trues(n)
    for i in 1:n, j in 1:p
        isfinite(raw[i, j]) || (valid[i] = false)
    end
    z = fill(NaN, n, p)
    start = 1
    while start <= n
        stop = start
        while stop < n && file[stop+1] == file[start]; stop += 1; end
        idxs = start:stop
        for j in 1:p
            good = filter(isfinite, raw[idxs, j])
            isempty(good) && continue
            μ = mean(good); σ = std(good); σ = σ > 0 ? σ : scale_fallback
            for i in idxs
                z[i, j] = isfinite(raw[i, j]) ? (raw[i, j] - μ) / σ : NaN
            end
        end
        start = stop + 1
    end
    if p > 1
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
    return ViewMatrix(permutedims(z), valid)
end

"""
    load_assignment_data(features_predictor, patches_bwd17; scale_fallback=1.0)

Reads a run's assignment/features_predictor.tsv and assignment/patches_bwd17.tsv.
"""
function load_assignment_data(features::AbstractString, patches::AbstractString; scale_fallback::Real=1.0)
    _, rows = lobe_table(features; required=vcat(["amplitude"], unique(vcat(BASE4, last(GMM_VIEW), "split_log_skew"))))
    pixels = [@sprintf("bwd_res_p%03d", i) for i in 1:289]
    _, prows = lobe_table(patches; required=pixels)
    keys_sorted = sort(collect(keys(rows)))
    Set(keys_sorted) == Set(keys(prows)) || error("Features and backward patches cover different lobes")
    coords = collect(range(-1.0, 1.0; length=17))
    moments = Dict(k => _negative_moment([_value(prows[k][c]) for c in pixels], coords) for k in keys_sorted)
    file = first.(keys_sorted); lobe = last.(keys_sorted)
    amplitude = [_value(rows[k]["amplitude"]) for k in keys_sorted]
    all(isfinite, amplitude) || error("Non-finite amplitude")
    column(name) = name == "bwd_neg_com_t" ? [moments[k][1] for k in keys_sorted] :
                   name == "bwd_neg_diag45" ? [moments[k][2] for k in keys_sorted] :
                   [_value(rows[k][name]) for k in keys_sorted]
    rawview(names) = reduce(hcat, [column(n) for n in names])
    kviews = [_view(rawview(names), file, 1.0) for (_, names) in KMEANS_VIEWS]
    graw = rawview(last(GMM_VIEW))
    gview = _view(graw, file, Float64(scale_fallback))
    return AssignmentData(file, lobe, amplitude, kviews, gview, graw, Float64(scale_fallback))
end

# --- GMM, identical arithmetic to build_labelfree_gmm_predictions.jl; the
# Cholesky factor of each component is computed once per iteration instead of
# once per lobe, and lobe vectors use buffers (same factor, same operations,
# same result; checked against the production script in the tests).

_chol(S) = cholesky(Symmetric(S) + 1e-8 * I).L

"d = X[:, i] - M[:, j] without allocation."
function _centered!(d, X, i, M, j)
    @inbounds for r in eachindex(d)
        d[r] = X[r, i] - M[r, j]
    end
    return d
end

function _gmm_fit(X::Matrix{Float64}, seed::Int; ridge, max_iter, tol, min_mass)
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
    d = Vector{Float64}(undef, p)
    centered = similar(X); scaled = similar(X); product = Matrix{Float64}(undef, p, p)
    for iter in 1:max_iter
        log_resp = zeros(n, k)
        for j in 1:k
            L = _chol(covs[j]); logdet2 = 2 * sum(log.(diag(L))); lw = log(weights[j])
            for i in 1:n
                _centered!(d, X, i, means, j)
                ldiv!(L, d)
                log_resp[i, j] = lw + -0.5 * (p * log(2π) + logdet2 + dot(d, d))
            end
        end
        # k = 2: row maxima and log-sum-exp written out (same operations as the
        # production script's row slices).
        ll = 0.0
        for i in 1:n
            a, b = log_resp[i, 1], log_resp[i, 2]
            m = max(a, b)
            ll += m + log(exp(a - m) + exp(b - m))
        end
        abs(ll - prev_ll) < tol * (1 + abs(prev_ll)) && break
        prev_ll = ll
        for i in 1:n
            a, b = log_resp[i, 1], log_resp[i, 2]
            m = max(a, b)
            a -= m; b -= m
            s = log(exp(a) + exp(b))
            log_resp[i, 1] = a - s; log_resp[i, 2] = b - s
        end
        resp = exp.(log_resp)
        for j in 1:k
            nk = sum(resp[:, j])
            nk > min_mass || continue
            weights[j] = nk / n
            means[:, j] = X * resp[:, j] / nk
            # (X .- mean) * Diagonal(resp) * (X .- mean)' in preallocated buffers
            centered .= X .- means[:, j]
            scaled .= centered .* transpose(resp[:, j])
            covs[j] = mul!(product, scaled, centered') / nk
            covs[j] += ridge * I
        end
    end
    return means, covs, weights
end

function _self_train(X::Matrix{Float64}, means, covs, weights; iters::Int, ridge)
    p, n = size(X); k = size(means, 2)
    for iteration in 1:iters
        d2 = zeros(n, k); d = Vector{Float64}(undef, p)
        for j in 1:k
            L = _chol(covs[j])
            for i in 1:n
                _centered!(d, X, i, means, j)
                ldiv!(L, d)
                d2[i, j] = dot(d, d)
            end
        end
        assign = [argmin(view(d2, i, :)) for i in 1:n]
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

"Final hard assignment of the columns of X (1 or 2)."
function _gmm_assign(X::AbstractMatrix{Float64}, g::GMMSeed, selftrain::Int)
    k = length(g.weights); p = size(X, 1)
    Ls = [_chol(g.covs[c]) for c in 1:k]
    logdet2 = [2 * sum(log.(diag(L))) for L in Ls]
    out = Vector{Int}(undef, size(X, 2))
    scores = zeros(k); e = zeros(k); d = Vector{Float64}(undef, p)
    for j in axes(X, 2)
        for c in 1:k
            _centered!(d, X, j, g.means, c)
            ldiv!(Ls[c], d)
            scores[c] = selftrain > 0 ? log(g.weights[c]) - 0.5 * dot(d, d) :
                log(g.weights[c]) + -0.5 * (p * log(2π) + logdet2[c] + dot(d, d))
        end
        e .= exp.(scores .- maximum(scores))
        out[j] = argmax(e)
    end
    return out
end

"""
    train_frozen(data, train_rows, cfg; fusion=nothing)

Learns k-means centroids, GMM components and their physical naming (higher
mean raw amplitude = class 1) on `train_rows` (indices into `data`). Settings
come from the assignment config `cfg` (seeds, selftrain, GMM numerics).
Returns a FrozenModel without emission rates (theta0 = theta1 = NaN). With
equal_weights, the GMM mixing weights are set to 1/2 after learning, for
naming and scoring (no class fraction carried from the training scans).
"""
function train_frozen(data::AssignmentData, train_rows::AbstractVector{<:Integer}, cfg::AbstractDict; equal_weights::Bool=false)
    sel = cfg["selection"]; m = cfg["model"]
    first_seed = Int(sel["first_seed"])
    kcenters = Vector{Vector{Matrix{Float64}}}(); khigh = Vector{Vector{Int}}()
    for V in data.kmeans
        train = [i for i in train_rows if V.valid[i]]
        cs = Matrix{Float64}[]; hs = Int[]
        if length(train) >= 2
            X = V.X[:, train]
            for seed in first_seed:(first_seed + Int(sel["kmeans_seeds"]) - 1)
                km = kmeans(X, 2; maxiter=200, rng=MersenneTwister(seed), display=:none)
                amp = [mean(data.amplitude[train[j]] for j in eachindex(train) if km.assignments[j] == c) for c in 1:2]
                push!(cs, copy(km.centers)); push!(hs, amp[1] >= amp[2] ? 1 : 2)
            end
        end
        push!(kcenters, cs); push!(khigh, hs)
    end
    V = data.gmm
    train = [i for i in train_rows if V.valid[i]]
    seeds = GMMSeed[]
    selftrain = Int(sel["gmm_selftrain"])
    if length(train) >= 2
        X = V.X[:, train]
        ridge = Float64(m["gmm_covariance_ridge"])
        for seed in first_seed:(first_seed + Int(sel["gmm_seeds"]) - 1)
            means, covs, weights = _gmm_fit(X, seed; ridge, max_iter=Int(m["gmm_learning_maxiter"]),
                                            tol=Float64(m["gmm_learning_tolerance"]), min_mass=Float64(m["gmm_learning_min_mass"]))
            selftrain > 0 && ((means, covs, weights) = _self_train(X, means, covs, weights; iters=selftrain, ridge))
            equal_weights && (weights = fill(1.0 / length(weights), length(weights)))
            g = GMMSeed(means, covs, weights, 0)
            assign = _gmm_assign(X, g, selftrain)
            values = Dict(1 => Float64[], 2 => Float64[])
            for (j, i) in enumerate(train)
                push!(values[assign[j]], data.amplitude[i])
            end
            mean_amp = Dict(c => mean(v) for (c, v) in values if !isempty(v))
            high = length(mean_amp) == 2 ? first(sort(collect(keys(mean_amp)); by=c -> mean_amp[c], rev=true)) : 0
            push!(seeds, GMMSeed(means, covs, weights, high))
        end
    end
    return FrozenModel(kcenters, khigh, seeds, selftrain, NaN, NaN, NaN, Union{Nothing,FrozenFisher}[])
end

_round8(p) = isfinite(p) ? parse(Float64, @sprintf("%.8f", p)) : NaN

"""
    predict_frozen(model, data, rows; vote_threshold=0.5)

Per-scan calls for `rows`: (call, p_kmeans, p_gmm, p_vote) with call in
Int8 0/1, or -1 for "?" (a component without a valid view). Component
probabilities are rounded to 8 decimals as the production TSVs are.
"""
function predict_frozen(model::FrozenModel, data::AssignmentData, rows::AbstractVector{<:Integer}; vote_threshold::Real=0.5)
    n = length(rows)
    psum = zeros(n); pcount = zeros(Int, n)
    for (v, V) in enumerate(data.kmeans)
        cs = model.kmeans_centers[v]; hs = model.kmeans_high[v]
        isempty(cs) && continue
        buf = zeros(size(V.X, 1)); d = zeros(2)
        for (r, i) in enumerate(rows)
            V.valid[i] || continue
            votes = 0.0
            for (C, h) in zip(cs, hs)
                for c in 1:2
                    _centered!(buf, V.X, i, C, c)
                    d[c] = sum(abs2, buf)
                end
                votes += argmin(d) == h ? 1.0 : 0.0
            end
            p = votes / length(cs)
            isfinite(p) && (psum[r] += p; pcount[r] += 1)
        end
    end
    pk = [pcount[r] > 0 ? _round8(psum[r] / pcount[r]) : NaN for r in 1:n]
    pg = fill(NaN, n)
    V = data.gmm
    scored = [r for r in 1:n if V.valid[rows[r]]]
    if !isempty(scored) && !isempty(model.gmm)
        X = V.X[:, rows[scored]]
        votes = zeros(length(scored)); counts = zeros(Int, length(scored))
        for g in model.gmm
            g.high == 0 && continue
            assign = _gmm_assign(X, g, model.selftrain)
            votes .+= (assign .== g.high); counts .+= 1
        end
        for (j, r) in enumerate(scored)
            counts[j] > 0 && (pg[r] = _round8(votes[j] / counts[j]))
        end
    end
    call = fill(Int8(-1), n); pv = fill(NaN, n)
    for r in 1:n
        isfinite(pk[r]) && isfinite(pg[r]) || continue
        p = (pk[r] + pg[r]) / 2
        pv[r] = p
        call[r] = p >= vote_threshold ? Int8(1) : Int8(0)
    end
    return (call=call, p_kmeans=pk, p_gmm=pg, p_vote=pv)
end

"""
    fit_emission_rates(calls, physical, settings) -> fit_latent_class result

`calls`: per-scan calls (Int8, -1 = "?"); `physical`: physical-lobe key per
call ("" when the lobe is not on a fused track). Uses lobes with at least
settings.min_scans classified calls, as build_molecule_fusion.jl does.
"""
function fit_emission_rates(calls::AbstractVector{Int8}, physical::AbstractVector{<:AbstractString}, st::FusionSettings)
    k, n = _physical_counts(calls, physical, st.min_scans)
    isempty(k) && error("No physical lobe observed in at least $(st.min_scans) scans")
    fused = sort(collect(keys(k)))
    return fit_latent_class([k[key] for key in fused], [n[key] for key in fused], st)
end

function _physical_counts(calls, physical, min_scans)
    k = Dict{String,Int}(); n = Dict{String,Int}()
    for (c, key) in zip(calls, physical)
        isempty(key) && continue
        c >= 0 || continue
        k[key] = get(k, key, 0) + c; n[key] = get(n, key, 0) + 1
    end
    for key in collect(keys(n))
        n[key] >= min_scans || (delete!(n, key); delete!(k, key))
    end
    return k, n
end

"""
    fuse_calls(calls, physical, theta0, theta1, pi, min_scans) -> Vector{Int8}

Posterior of the latent-class model with fixed parameters; a physical lobe
with at least min_scans classified calls gets the fused call (posterior >= 0.5),
every other row keeps its per-scan call.
"""
function fuse_calls(calls::AbstractVector{Int8}, physical::AbstractVector{<:AbstractString},
                    theta0::Real, theta1::Real, pi::Real, min_scans::Integer)
    k, n = _physical_counts(calls, physical, min_scans)
    post = Dict{String,Float64}()
    for (key, kk) in k
        nn = n[key]
        l1 = log(pi) + kk * log(theta1) + (nn - kk) * log(1 - theta1)
        l0 = log(1 - pi) + kk * log(theta0) + (nn - kk) * log(1 - theta0)
        mx = max(l0, l1)
        post[key] = exp(l1 - mx) / (exp(l1 - mx) + exp(l0 - mx))
    end
    return [haskey(post, key) ? (post[key] >= 0.5 ? Int8(1) : Int8(0)) : c for (c, key) in zip(calls, physical)]
end

with_fusion(m::FrozenModel, theta0, theta1, pi) =
    FrozenModel(m.kmeans_centers, m.kmeans_high, m.gmm, m.selftrain, Float64(theta0), Float64(theta1), Float64(pi), m.fisher)
with_fisher(m::FrozenModel, fisher) =
    FrozenModel(m.kmeans_centers, m.kmeans_high, m.gmm, m.selftrain, m.theta0, m.theta1, m.pi, fisher)

"""
    write_frozen_model(path, model)

Long TSV: part, view, seed, index, i, j, value (%.17g). part is kmeans_center,
kmeans_high, gmm_mean, gmm_cov, gmm_weight, gmm_high, fisher_present,
fisher_wp, fisher_mid (view = parity model) or scalar.
"""
function write_frozen_model(path::AbstractString, m::FrozenModel)
    (ispath(path) || islink(path)) && error("Output already exists: $path")
    open(path, "w") do io
        println(io, join(["part", "view", "seed", "index", "i", "j", "value"], '\t'))
        line(a...) = println(io, join(string.(a), '\t'))
        f(x) = @sprintf("%.17g", x)
        for (v, cs) in enumerate(m.kmeans_centers), (s, C) in enumerate(cs)
            for c in axes(C, 2), i in axes(C, 1); line("kmeans_center", v, s, c, i, 0, f(C[i, c])); end
            line("kmeans_high", v, s, 0, 0, 0, m.kmeans_high[v][s])
        end
        for (s, g) in enumerate(m.gmm)
            for c in eachindex(g.weights)
                for i in axes(g.means, 1); line("gmm_mean", 0, s, c, i, 0, f(g.means[i, c])); end
                for i in axes(g.covs[c], 1), j in axes(g.covs[c], 2); line("gmm_cov", 0, s, c, i, j, f(g.covs[c][i, j])); end
                line("gmm_weight", 0, s, c, 0, 0, f(g.weights[c]))
            end
            line("gmm_high", 0, s, 0, 0, 0, g.high)
        end
        for (k, ff) in enumerate(m.fisher)
            line("fisher_present", k, 0, 0, 0, 0, ff === nothing ? 0 : 1)
            ff === nothing && continue
            for i in eachindex(ff.w_p); line("fisher_wp", k, 0, 0, i, 0, f(ff.w_p[i])); end
            for i in eachindex(ff.mid); line("fisher_mid", k, 0, 0, i, 0, f(ff.mid[i])); end
        end
        for (name, x) in (("selftrain", m.selftrain), ("theta0", m.theta0), ("theta1", m.theta1), ("pi", m.pi))
            line("scalar_" * name, 0, 0, 0, 0, 0, x isa Integer ? x : f(x))
        end
    end
    return path
end

function read_frozen_model(path::AbstractString)
    lines = readlines(path)
    split(lines[1], '\t') == ["part", "view", "seed", "index", "i", "j", "value"] || error("Not a frozen model: $path")
    rows = [split(l, '\t') for l in lines[2:end] if !isempty(l)]
    ints(r) = parse.(Int, r[2:6])
    nview = maximum((ints(r)[1] for r in rows if startswith(r[1], "kmeans")); init=0)
    kc = [Matrix{Float64}[] for _ in 1:nview]; kh = [Int[] for _ in 1:nview]
    gm = Dict{Int,Dict{Symbol,Any}}(); scal = Dict{String,String}()
    fpresent = Dict{Int,Bool}(); fwp = Dict{Int,Dict{Int,Float64}}(); fmid = Dict{Int,Dict{Int,Float64}}()
    for r in rows
        v, s, c, i, j = ints(r)
        if r[1] == "fisher_present"
            fpresent[v] = r[7] == "1"
        elseif r[1] == "fisher_wp"
            get!(fwp, v, Dict{Int,Float64}())[i] = parse(Float64, r[7])
        elseif r[1] == "fisher_mid"
            get!(fmid, v, Dict{Int,Float64}())[i] = parse(Float64, r[7])
        elseif r[1] == "kmeans_center"
            while length(kc[v]) < s; push!(kc[v], zeros(0, 2)); end
            C = kc[v][s]
            size(C, 1) < i && (C = vcat(C, zeros(i - size(C, 1), 2)); kc[v][s] = C)
            C[i, c] = parse(Float64, r[7])
        elseif r[1] == "kmeans_high"
            while length(kh[v]) < s; push!(kh[v], 0); end
            kh[v][s] = parse(Int, r[7])
        elseif startswith(r[1], "gmm_")
            g = get!(gm, s, Dict{Symbol,Any}(:mean => Dict{Tuple{Int,Int},Float64}(), :cov => Dict{Tuple{Int,Int,Int},Float64}(),
                                             :weight => Dict{Int,Float64}(), :high => 0))
            r[1] == "gmm_mean" && (g[:mean][(i, c)] = parse(Float64, r[7]))
            r[1] == "gmm_cov" && (g[:cov][(c, i, j)] = parse(Float64, r[7]))
            r[1] == "gmm_weight" && (g[:weight][c] = parse(Float64, r[7]))
            r[1] == "gmm_high" && (g[:high] = parse(Int, r[7]))
        elseif startswith(r[1], "scalar_")
            scal[r[1][8:end]] = r[7]
        end
    end
    seeds = GMMSeed[]
    for s in sort(collect(keys(gm)))
        g = gm[s]; k = length(g[:weight]); p = maximum(first.(collect(keys(g[:mean]))))
        means = [g[:mean][(i, c)] for i in 1:p, c in 1:k]
        covs = [[g[:cov][(c, i, j)] for i in 1:p, j in 1:p] for c in 1:k]
        push!(seeds, GMMSeed(means, covs, [g[:weight][c] for c in 1:k], g[:high]))
    end
    vec_of(d) = [d[i] for i in 1:length(d)]
    fisher = Union{Nothing,FrozenFisher}[fpresent[k] ? FrozenFisher(vec_of(fwp[k]), vec_of(fmid[k])) : nothing
                                         for k in sort(collect(keys(fpresent)))]
    return FrozenModel(kc, kh, seeds, parse(Int, scal["selftrain"]), parse(Float64, scal["theta0"]),
                       parse(Float64, scal["theta1"]), parse(Float64, scal["pi"]), fisher)
end

# --- fully frozen model: Fisher feature, reference rule, training and application

"Forward 17x17 residual patches in the Fisher layout, with the production validity rule."
struct FisherData
    options::FisherOptions
    table::PatchTable
    valid::BitVector
end

"Frozen-model settings of a frozen_reference config, checked."
function load_frozen_settings(cfg::AbstractDict)
    sel = cfg["selection"]
    get(sel, "assignment_model", "") == "frozen_reference" || error("Config is not a frozen_reference assignment")
    get(sel, "fusion_emission_source", "") == "reference_sessions_all_scans" || error("Unsupported fusion_emission_source")
    get(sel, "reference_scans", "") == "corroborated_counts" || error("reference_scans must be corroborated_counts")
    cap = get(sel, "reference_max_scans_per_track", 0)
    cap isa Integer && cap > 0 || error("reference_max_scans_per_track must be a positive integer")
    get(sel, "reference_scan_order", "") == "file_name" || error("reference_scan_order must be file_name")
    get(sel, "frozen_fisher", "") == "reference_lobe_parity" || error("frozen_fisher must be reference_lobe_parity")
    weights = get(sel, "frozen_gmm_weights", "")
    weights in ("equal", "learned") || error("frozen_gmm_weights must be equal or learned")
    pi = Float64(sel["fusion_prior_pi"]); 0 < pi < 1 || error("fusion_prior_pi must be in (0, 1)")
    return (max_per_track=Int(cap), equal_weights=weights == "equal", pi=pi, fusion=load_fusion_settings(cfg))
end

"""
    load_fisher_data(patches_fwd17, cfg, data)

Reads assignment/patches_fwd17.tsv (prefix res) with the config's Fisher
settings; its lobes must be those of `data`, in the same order.
"""
function load_fisher_data(path::AbstractString, cfg::AbstractDict, data::AssignmentData)
    options = load_fisher_config(cfg)
    options.patch_projection == "none" || error("The frozen Fisher supports fisher_patch_projection = none only")
    options.cv_scheme == "lobe_parity" || error("The frozen Fisher follows fisher_cv_scheme = lobe_parity")
    options.training_support == "all_admissible" || error("The frozen Fisher needs assignment_training_support = all_admissible")
    table = load_patches(path, "res", options)
    table.keys == collect(zip(data.file, data.lobe)) || error("Forward patches and features cover different lobes")
    valid = BitVector([isempty(table.invalid_reasons[i]) && all(isfinite, table.X[i, :]) && isfinite(table.amplitudes[i])
                       for i in eachindex(table.keys)])
    return FisherData(options, table, valid)
end

"Two parity models fitted on the valid reference rows (nothing when a fit fails)."
function train_fisher(fd::FisherData, rows::AbstractVector{<:Integer})
    out = Union{Nothing,FrozenFisher}[]
    for parity in (0, 1)
        train = [i for i in sort(rows) if fd.valid[i] && mod(fd.table.keys[i][2], 2) == parity]
        model = try
            fit_fisher(fd.table.X[train, :], fd.table.amplitudes[train], fd.options)
        catch err
            err isa FisherFitError || rethrow()
            nothing
        end
        push!(out, model === nothing ? nothing : FrozenFisher(model.w_p, model.mid))
    end
    return out
end

"""
    fisher_feature(fd, fisher, i)

emp_fisher of lobe i: the mirrored Fisher margin of the opposite-parity
model, rounded and negated as the production files carry it (fisher_cv.tsv
%.6f, then features_predictor.tsv %.8g of the negated score). NaN if the
patch is invalid or the model is missing.
"""
function fisher_feature(fd::FisherData, fisher::AbstractVector, i::Integer)
    fd.valid[i] || return NaN
    m = fisher[mod(fd.table.keys[i][2], 2) == 1 ? 1 : 2]
    m === nothing && return NaN
    x = fd.table.X[i, :]
    v = max(dot(x - m.mid, m.w_p), dot(flip_u_disk(x, fd.table.grid) - m.mid, m.w_p))
    isfinite(v) || return NaN
    return parse(Float64, @sprintf("%.8g", -parse(Float64, @sprintf("%.6f", v))))
end

"`data` with the emp_fisher column replaced (GMM view rebuilt per scan)."
function with_emp_fisher(data::AssignmentData, values::AbstractVector{<:Real})
    length(values) == length(data.file) || error("One emp_fisher value per lobe is required")
    raw = copy(data.gmm_raw)
    raw[:, findfirst(==("emp_fisher"), last(GMM_VIEW))] .= values
    return AssignmentData(data.file, data.lobe, data.amplitude, data.kmeans, _view(raw, data.file, data.scale_fallback),
                          raw, data.scale_fallback)
end

"""
    reference_files(rows, max_per_track)

Reference scans from molecule-consensus summary rows (filepath, track,
count_rule): the corroborated scans (count_rule = agrees), at most
max_per_track per track, first in file-name (acquisition-number) order.
"""
function reference_files(rows, max_per_track::Integer)
    by_track = Dict{String,Vector{String}}()
    for r in rows
        r["count_rule"] == "agrees" || continue
        push!(get!(by_track, r["track"], String[]), basename(r["filepath"]))
    end
    return Set(f for v in values(by_track) for f in sort(v)[1:min(max_per_track, length(v))])
end

"""
    train_reference_model(data, fd, ref_rows, theta_rows, physical, cfg)

Fully frozen model: Fisher per parity on `ref_rows`, emp_fisher recomputed
for every lobe, k-means/GMM/naming on `ref_rows`, emission rates fitted on
the per-scan calls of `theta_rows`, prior pi from the config. Throws if the
emission rates cannot be fitted. freeze_fisher=false keeps the run's own
(transductive) emp_fisher; it exists for diagnostics, not production.
"""
function train_reference_model(data::AssignmentData, fd::FisherData, ref_rows, theta_rows, physical, cfg::AbstractDict;
                               freeze_fisher::Bool=true)
    st = load_frozen_settings(cfg)
    # freeze_fisher=false (diagnostics only): keep the run's transductive emp_fisher
    fisher = freeze_fisher ? train_fisher(fd, ref_rows) : Union{Nothing,FrozenFisher}[]
    d = freeze_fisher ? with_emp_fisher(data, [fisher_feature(fd, fisher, i) for i in eachindex(data.file)]) : data
    model = with_fisher(train_frozen(d, ref_rows, cfg; equal_weights=st.equal_weights), fisher)
    p = predict_frozen(model, d, theta_rows; vote_threshold=cfg["selection"]["vote_threshold"])
    fit = fit_emission_rates(p.call, physical[theta_rows], st.fusion)
    return with_fusion(model, fit.theta0, fit.theta1, st.pi)
end

"""
    apply_reference_model(model, data, fd, rows, physical, cfg) -> (fused, per_scan)

Per-scan calls of `rows` and their fusion within tracks (physical-lobe keys),
with the model's Fisher, clusters, emission rates and prior.
"""
function apply_reference_model(model::FrozenModel, data::AssignmentData, fd::FisherData, rows, physical, cfg::AbstractDict;
                               freeze_fisher::Bool=true)
    length(model.fisher) == (freeze_fisher ? 2 : 0) || error("Fisher parity models do not match freeze_fisher=$freeze_fisher")
    st = load_frozen_settings(cfg)
    d = freeze_fisher ? with_emp_fisher(data, [fisher_feature(fd, model.fisher, i) for i in eachindex(data.file)]) : data
    p = predict_frozen(model, d, rows; vote_threshold=cfg["selection"]["vote_threshold"])
    return fuse_calls(p.call, physical[rows], model.theta0, model.theta1, model.pi, st.fusion.min_scans), p
end

"""
    load_run_assignment(run, cfg) -> (data, fisher, physical, consensus)

Assignment inputs of a run of run_molecule_consensus_chitosan.jl: features
and patches (assignment/), physical-lobe keys (fusion/predictions_fused.tsv)
and the final consensus summary rows. Label-free.
"""
function load_run_assignment(run::AbstractString, cfg::AbstractDict)
    a = joinpath(run, "assignment")
    data = load_assignment_data(joinpath(a, "features_predictor.tsv"), joinpath(a, "patches_bwd17.tsv");
                                scale_fallback=cfg["preprocessing"]["gmm_scale_fallback"])
    fd = load_fisher_data(joinpath(a, "patches_fwd17.tsv"), cfg, data)
    _, fused = lobe_table(joinpath(run, "fusion", "predictions_fused.tsv"); required=["physical_lobe"])
    physical = [fused[(data.file[i], data.lobe[i])]["physical_lobe"] for i in eachindex(data.file)]
    _, cons = read_table(joinpath(run, "consensus", "consensus_summary_final.tsv"))
    Set(basename(r["filepath"]) for r in cons) == Set(data.file) || error("Consensus and assignment cover different scans")
    return (data=data, fisher=fd, physical=physical, consensus=cons)
end

end

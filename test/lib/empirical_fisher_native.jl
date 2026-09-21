"""
Native reconstruction of `empirical_fisher_mold.py` (Julia 1.13).

Only patch pixels and lobe parity enter the calculation. PCA dimension, both
regularizers, EM limit/tolerance, seed, and grid widths come from the config.
The two full-covariance mixture weights are estimated, never fixed.
Optional affine-disk projection removes the fixed [1,t,u] subspace before PCA
and before each held-out score. Original center amplitudes still order clusters;
the projection removes affine molecular signal as well as possible nuisance.

This is NOT a byte-identical sklearn implementation. Initialization uses Julia's
MersenneTwister, one two-center k-means++ draw, and Lloyd iterations until the
partition is unchanged (bounded by the configured GMM iteration limit). sklearn
uses NumPy's random stream, greedy k-means++ candidates, and its own k-means
stopping rule. EM uses exact responsibility sums instead of sklearn's machine-
epsilon mass addition. SVD/BLAS rounding can differ as well. No seed search,
model selection, score calibration, or reference predictions are used.

Historical conventions deliberately retained:
* `legacy_row_major` reads numbered pixels as Python reshape(side, side).
  The current extractor writes u-outer/t-inner, so the legacy [:, ::-1]
  mirror reverses physical t, despite its historical name `flip_u_disk`.
  The opt-in `physical_u_outer_t_inner` layout instead reverses physical u.
  Both layouts keep exactly the same disk pixels, PCA and fit conventions;
  only the reflected held-out patch changes. The legacy config is unchanged.
* PCA fits centered training patches. The noise covariance is the sample
  covariance of ALL latent training rows, not pooled within-cluster covariance.
* The two hard-cluster patch means are ordered by center-pixel amplitude.
  The weight points from low amplitude to high amplitude.
* `mid` is mean(X_centered), approximately zero, NOT the original training
  mean or the midpoint of the cluster means. Scores use raw held-out patches.
  The explicit `training_mean` alternative instead subtracts the training-patch
  mean at scoring, with unchanged PCA, mixture, weights and amplitude mapping.
* Scores are max(original, mirrored), with unchanged sign, rounded to six
  decimals only when written. The old negative-score => 1 rule is not applied
  here: this file produces a margin, not a chemical prediction.
"""
module EmpiricalFisherNative

using LinearAlgebra
using Statistics
using Random
using Printf
using TOML
include(joinpath(@__DIR__, "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment: load_training_policy, validate_training_mask

export FisherOptions, FisherGrid, PatchTable, FisherModel, FisherScore,
       load_fisher_config, fisher_grid, load_patches, flip_u_disk,
       fit_fisher, score, maxmirror_score, cv_scores, write_scores

struct FisherOptions
    pca_components::Int
    noise_regularization::Float64
    gmm_regularization::Float64
    gmm_maxiter::Int
    gmm_tolerance::Float64
    seed::Int
    half_nm::Float64
    step_nm::Float64
    layout::String
    patch_projection::String
    projection_zero_l1::Float64
    score_center::String
    patch_support::String
    training_support::String

    function FisherOptions(pca_components, noise_regularization,
                           gmm_regularization, gmm_maxiter, gmm_tolerance,
                           seed, half_nm, step_nm, layout, patch_projection, projection_zero_l1, score_center,
                           patch_support, training_support)
        for (name, value) in (("fisher_pca_components", pca_components),
                              ("fisher_gmm_maxiter", gmm_maxiter))
            value isa Integer && !(value isa Bool) && value > 0 ||
                throw(ArgumentError("$name must be a positive integer"))
        end
        seed isa Integer && !(seed isa Bool) && seed >= 0 ||
            throw(ArgumentError("fisher_seed must be a nonnegative integer"))
        for (name, value) in (("fisher_noise_regularization", noise_regularization),
                              ("fisher_gmm_regularization", gmm_regularization),
                              ("fisher_gmm_tolerance", gmm_tolerance),
                              ("mold_half_nm", half_nm), ("mold_step_nm", step_nm))
            value isa Real && !(value isa Bool) && isfinite(value) && value > 0 ||
                throw(ArgumentError("$name must be positive and finite"))
        end
        layout in ("legacy_row_major", "physical_u_outer_t_inner") ||
            throw(ArgumentError("unsupported fisher_layout: $layout"))
        patch_projection in ("none", "affine_disk") ||
            throw(ArgumentError("unsupported fisher_patch_projection: $patch_projection"))
        score_center in ("legacy_centered_mean", "training_mean") ||
            throw(ArgumentError("unsupported fisher_score_center: $score_center"))
        patch_support in ("full_square", "complete_disk_symmetric") ||
            throw(ArgumentError("unsupported assignment_patch_support: $patch_support"))
        training_support in ("all_admissible", "complete_patches") ||
            throw(ArgumentError("unsupported assignment_training_support: $training_support"))
        projection_zero_l1 isa Real && !(projection_zero_l1 isa Bool) &&
            isfinite(projection_zero_l1) && projection_zero_l1 >= 0 ||
            throw(ArgumentError("fisher_projection_zero_l1 must be finite and nonnegative"))
        ratio = half_nm / step_nm
        isfinite(ratio) && 1 <= ratio < (typemax(Int) - 1) / 2 ||
            throw(ArgumentError("mold_half_nm / mold_step_nm must give a finite grid"))
        nhalf = round(Int, ratio)
        # This tolerance only absorbs binary representation of an integer grid.
        isapprox(ratio, nhalf; rtol=8eps(Float64), atol=0) ||
            throw(ArgumentError("mold_half_nm must be an integer multiple of mold_step_nm"))
        isfinite(half_nm^2) && half_nm^2 > 0 ||
            throw(ArgumentError("mold_half_nm squared must be positive and finite"))
        new(Int(pca_components), Float64(noise_regularization),
            Float64(gmm_regularization), Int(gmm_maxiter), Float64(gmm_tolerance),
            Int(seed), Float64(half_nm), Float64(step_nm), String(layout), String(patch_projection),
            Float64(projection_zero_l1), String(score_center), String(patch_support), String(training_support))
    end
end

"Load all scientific settings explicitly; there are no fallback model values."
function load_fisher_config(config::AbstractDict)
    model = get(config, "model", Dict())
    pre = get(config, "preprocessing", Dict())
    required = ("fisher_pca_components", "fisher_noise_regularization",
                "fisher_gmm_regularization", "fisher_gmm_maxiter",
                "fisher_gmm_tolerance", "fisher_seed", "mold_half_nm", "mold_step_nm")
    for name in required
        haskey(model, name) || throw(ArgumentError("missing [model] $name"))
    end
    for name in ("fisher_layout", "fisher_patch_projection", "assignment_patch_support")
        haskey(pre, name) || throw(ArgumentError("missing [preprocessing] $name"))
    end
    haskey(model, "fisher_projection_zero_l1") ||
        throw(ArgumentError("missing [model] fisher_projection_zero_l1"))
    haskey(model, "fisher_score_center") ||
        throw(ArgumentError("missing [model] fisher_score_center"))
    options = FisherOptions((model[name] for name in required)...,
                            pre["fisher_layout"], pre["fisher_patch_projection"], model["fisher_projection_zero_l1"],
                            model["fisher_score_center"], pre["assignment_patch_support"], load_training_policy(config))
    grid = fisher_grid(options)
    options.pca_components <= length(grid.disk_indices) ||
        throw(ArgumentError("fisher_pca_components exceeds the number of disk pixels"))
    return options
end
load_fisher_config(path::AbstractString) = load_fisher_config(TOML.parsefile(path))

struct FisherGrid
    side::Int
    coords::Vector{Float64}
    disk_indices::Vector{Int}       # one-based row-major full-patch indices
    mirror_indices::Vector{Int}     # disk index, or 0 for zero padding
    center_index::Int               # full-patch index
end

"Keep the NumPy disk enumeration; select the explicitly configured mirror axis."
function fisher_grid(options::FisherOptions)
    nhalf = round(Int, options.half_nm / options.step_nm)
    side = 2nhalf + 1
    start = -options.half_nm
    # NumPy arange uses the representable (start + step) - start increment.
    increment = (start + options.step_nm) - start
    coords = [start + i * increment for i in 0:(side - 1)]
    disk = Int[]
    for row in 1:side, col in 1:side
        coords[row]^2 + coords[col]^2 <= options.half_nm^2 &&
            push!(disk, (row - 1) * side + col)
    end
    positions = Dict(index => i for (i, index) in enumerate(disk))
    mirror = Int[]
    for index in disk
        row, col = divrem(index - 1, side)
        reflected = if options.layout == "legacy_row_major"
            row * side + (side - col)              # physical t -> -t
        else
            (side - 1 - row) * side + col + 1       # physical u -> -u
        end
        push!(mirror, get(positions, reflected, 0))
    end
    if options.layout == "physical_u_outer_t_inner"
        # A maximum over original/reflected patches is reflection-invariant
        # only when the disk is closed under an involutive reflection. Do not
        # silently zero-pad a different grid in this scientific variant.
        all(>(0), mirror) && mirror[mirror] == collect(eachindex(disk)) ||
            throw(ArgumentError("Fisher disk is not closed under transverse reflection"))
    end
    return FisherGrid(side, coords, disk, mirror, (side^2 + 1) ÷ 2)
end

struct PatchTable
    keys::Vector{Tuple{String,Int}}
    X::Matrix{Float64}              # rows follow sorted keys, columns disk order
    amplitudes::Vector{Float64}     # center pixels of the original patches
    invalid_reasons::Vector{String}
    grid::FisherGrid
end

"""
    load_patches(path, prefix, options) -> PatchTable

Keep every unique (basename(file), positive lobe) key, including incomplete and
nonfinite patches. Bad keys and duplicate keys are fatal schema errors. Missing
pixel columns and bad pixel values invalidate rows, not their keys. All full-
patch pixels must be finite in full_square mode. The complete_disk_symmetric
mode requires only the actual scoring disk to be finite (including its center).
Nonfinite/nonnumeric values outside it are ignored, never imputed. Missing schema
columns remain invalid in both modes.
"""
function load_patches(path::AbstractString, prefix::AbstractString, options::FisherOptions)
    isempty(prefix) && throw(ArgumentError("patch prefix must not be empty"))
    grid = fisher_grid(options)
    required_pixels = options.patch_support == "full_square" ? Set(1:grid.side^2) : Set(grid.disk_indices)
    pixels = [@sprintf("%s_p%03d", prefix, i) for i in 1:grid.side^2]
    keys = Tuple{String,Int}[]
    data = Vector{Float64}[]
    amplitudes = Float64[]
    reasons = String[]
    seen = Set{Tuple{String,Int}}()
    open(path, "r") do io
        eof(io) && throw(ArgumentError("empty patch TSV: missing file/lobe header"))
        header = split(chomp(readline(io)), '\t'; keepempty=true)
        length(unique(header)) == length(header) ||
            throw(ArgumentError("duplicate column names in patch TSV"))
        index = Dict(name => i for (i, name) in enumerate(header))
        all(haskey(index, name) for name in ("file", "lobe")) ||
            throw(ArgumentError("patch TSV requires file and lobe columns"))
        for (line_index, line) in enumerate(eachline(io))
            isempty(strip(line)) && continue
            values = split(chomp(line), '\t'; keepempty=true)
            length(values) <= length(header) ||
                throw(ArgumentError("too many fields on patch TSV line $(line_index + 1)"))
            value(name) = get(index, name, 0) in 1:length(values) ? values[index[name]] : ""
            file = basename(strip(value("file")))
            lobe = tryparse(Int, strip(value("lobe")))
            !isempty(file) && !(file in (".", "..")) ||
                throw(ArgumentError("missing file on patch TSV line $(line_index + 1)"))
            lobe !== nothing && lobe > 0 ||
                throw(ArgumentError("lobe must be a positive integer on line $(line_index + 1)"))
            key = (String(file), lobe)
            key in seen && throw(ArgumentError("duplicate patch key: $key"))
            push!(seen, key)
            patch = fill(NaN, grid.side^2)
            reason = ""
            for (i, name) in enumerate(pixels)
                if !haskey(index, name)
                    reason = "missing_patch_column:$name"
                    break
                end
                text = strip(value(name))
                number = tryparse(Float64, text)
                if !(i in required_pixels)
                    continue
                end
                if number === nothing
                    reason = "invalid_patch_value:$name"
                    break
                elseif !isfinite(number)
                    reason = "nonfinite_patch_value:$name"
                    break
                end
                patch[i] = number
            end
            push!(keys, key)
            push!(data, isempty(reason) ? patch[grid.disk_indices] : fill(NaN, length(grid.disk_indices)))
            push!(amplitudes, isempty(reason) ? patch[grid.center_index] : NaN)
            push!(reasons, reason)
        end
    end
    order = sortperm(keys)
    X = Matrix{Float64}(undef, length(keys), length(grid.disk_indices))
    for (i, j) in enumerate(order)
        X[i, :] = data[j]
    end
    return PatchTable(keys[order], X, amplitudes[order], reasons[order], grid)
end

"Configured disk mirror (historical API name); legacy t or opt-in physical u."
function flip_u_disk(x::AbstractVector, grid::FisherGrid)
    length(x) == length(grid.disk_indices) || throw(DimensionMismatch("wrong disk vector length"))
    return [index == 0 ? 0.0 : Float64(x[index]) for index in grid.mirror_indices]
end

struct FisherFitError <: Exception
    reason::String
end
Base.showerror(io::IO, error::FisherFitError) = print(io, "invalid Fisher fit: ", error.reason)

struct GMMFit
    weights::Vector{Float64}
    means::Matrix{Float64}          # component x latent dimension
    covariances::Vector{Matrix{Float64}}
    assignments::Vector{Int}
    converged::Bool
    iterations::Int
    lower_bound::Float64
end

# k=2 is structural. No composition prior or amplitude enters initialization.
function _initial_partition(Z::Matrix{Float64}, options::FisherOptions)
    n, d = size(Z)
    rng = MersenneTwister(options.seed)
    first = rand(rng, 1:n)
    distances = [sum(abs2, view(Z, i, :) .- view(Z, first, :)) for i in 1:n]
    total = sum(distances)
    isfinite(total) && total > 0 || throw(FisherFitError("no_patch_variation"))
    draw = rand(rng) * total
    second = findfirst(>(draw), cumsum(distances))
    second === nothing && (second = findlast(>(0), distances))
    centers = Z[[first, second], :]
    previous = zeros(Int, n)
    for _ in 1:options.gmm_maxiter
        assignment = [sum(abs2, view(Z, i, :) .- view(centers, 1, :)) <=
                      sum(abs2, view(Z, i, :) .- view(centers, 2, :)) ? 1 : 2 for i in 1:n]
        all(any(==(c), assignment) for c in 1:2) ||
            throw(FisherFitError("empty_initial_component"))
        assignment == previous && return assignment
        for c in 1:2
            centers[c, :] = vec(mean(Z[assignment .== c, :]; dims=1))
        end
        previous = assignment
    end
    return previous
end

"Weighted maximum-likelihood covariance (denominator sum responsibilities)."
function _m_step(Z::Matrix{Float64}, responsibilities::Matrix{Float64}, regularization::Float64)
    n, d = size(Z)
    size(responsibilities) == (n, 2) || throw(DimensionMismatch("expected two responsibility columns"))
    mass = vec(sum(responsibilities; dims=1))
    all(x -> isfinite(x) && x > 0, mass) || throw(FisherFitError("empty_soft_component"))
    weights = mass ./ sum(mass)
    means = (transpose(responsibilities) * Z) ./ mass
    covariances = Matrix{Float64}[]
    for c in 1:2
        residual = Z .- transpose(means[c, :])
        covariance = transpose(residual) * (residual .* responsibilities[:, c]) / mass[c]
        covariance[diagind(covariance)] .+= regularization
        all(isfinite, covariance) || throw(FisherFitError("nonfinite_gmm_covariance"))
        push!(covariances, covariance)
    end
    return weights, means, covariances
end

function _e_step(Z::Matrix{Float64}, weights, means, covariances)
    n, d = size(Z)
    log_probability = Matrix{Float64}(undef, n, 2)
    for c in 1:2
        factor = try
            cholesky(Symmetric(covariances[c]))
        catch error
            error isa PosDefException || rethrow()
            throw(FisherFitError("non_positive_gmm_covariance"))
        end
        whitened = factor.L \ transpose(Z .- transpose(means[c, :]))
        logdet = 2sum(log, diag(factor.L))
        log_probability[:, c] = log(weights[c]) .-
            (d * log(2pi) + logdet .+ vec(sum(abs2, whitened; dims=1))) ./ 2
    end
    largest = maximum(log_probability; dims=2)
    lognorm = largest .+ log.(sum(exp.(log_probability .- largest); dims=2))
    all(isfinite, lognorm) || throw(FisherFitError("nonfinite_gmm_likelihood"))
    responsibilities = exp.(log_probability .- lognorm)
    return responsibilities, mean(lognorm)
end

function _fit_gmm(Z::Matrix{Float64}, options::FisherOptions)
    initial = _initial_partition(Z, options)
    responsibilities = Float64.([a == c for a in initial, c in 1:2])
    weights, means, covariances = _m_step(Z, responsibilities, options.gmm_regularization)
    lower_bound = -Inf
    converged = false
    iterations = 0
    for iteration in 1:options.gmm_maxiter
        responsibilities, current = _e_step(Z, weights, means, covariances)
        weights, means, covariances = _m_step(Z, responsibilities, options.gmm_regularization)
        iterations = iteration
        converged = abs(current - lower_bound) < options.gmm_tolerance
        lower_bound = current
        converged && break
    end
    responsibilities, _ = _e_step(Z, weights, means, covariances)
    assignments = [responsibilities[i, 1] >= responsibilities[i, 2] ? 1 : 2 for i in axes(Z, 1)]
    all(any(==(c), assignments) for c in 1:2) ||
        throw(FisherFitError("empty_hard_component"))
    return GMMFit(weights, means, covariances, assignments, converged, iterations, lower_bound)
end

struct FisherModel
    w_p::Vector{Float64}
    mid::Vector{Float64}
    g0::Vector{Float64}             # centered low-amplitude patch mean
    g1::Vector{Float64}             # centered high-amplitude patch mean
    amplitude_means::Tuple{Float64,Float64}
    gmm::GMMFit
    projection_basis::Union{Nothing,Matrix{Float64}} # fixed disk basis, not learned from rows
end

"Orthonormal affine basis on the actual disk; no molecular pixels or labels enter."
function _affine_disk_basis(grid::FisherGrid)
    u = [grid.coords[div(index - 1, grid.side) + 1] for index in grid.disk_indices]
    t = [grid.coords[mod(index - 1, grid.side) + 1] for index in grid.disk_indices]
    return Matrix(qr(hcat(ones(length(t)), t, u)).Q)[:, 1:3]
end

_project_disk(x::AbstractVector, basis::AbstractMatrix) = x - basis * (transpose(basis) * x)

# Separate the linear algebra from GMM for fixed-partition numerical tests.
function _fisher_linear(Xc, Z, V, amplitudes, assignments, regularization)
    all(any(==(c), assignments) for c in 1:2) ||
        throw(FisherFitError("empty_hard_component"))
    amplitude_means = [mean(amplitudes[assignments .== c]) for c in 1:2]
    all(isfinite, amplitude_means) || throw(FisherFitError("nonfinite_component_amplitude"))
    amplitude_means[1] != amplitude_means[2] || throw(FisherFitError("amplitude_tie"))
    order = sortperm(amplitude_means)
    g0 = vec(mean(Xc[assignments .== order[1], :]; dims=1))
    g1 = vec(mean(Xc[assignments .== order[2], :]; dims=1))
    mid = vec(mean(Xc; dims=1))
    Zc = Z .- mean(Z; dims=1)
    covariance = transpose(Zc) * Zc / (size(Z, 1) - 1)
    covariance[diagind(covariance)] .+= regularization
    z0 = transpose(V) * (g0 - mid)
    z1 = transpose(V) * (g1 - mid)
    w_p = V * (Symmetric(covariance) \ (z1 - z0))
    all(isfinite, w_p) || throw(FisherFitError("nonfinite_fisher_vector"))
    any(!iszero, w_p) || throw(FisherFitError("zero_fisher_contrast"))
    return (w_p=w_p, mid=mid, g0=g0, g1=g1,
            amplitude_means=(amplitude_means[order[1]], amplitude_means[order[2]]))
end

"PCA -> native GMM2 -> amplitude-ordered hard means -> regularized Fisher solve."
function fit_fisher(X::AbstractMatrix{<:Real}, amplitudes::AbstractVector{<:Real}, options::FisherOptions)
    n, d = size(X)
    length(amplitudes) == n || throw(DimensionMismatch("one center amplitude is required per patch"))
    n >= max(2, options.pca_components) || throw(FisherFitError("too_few_training_rows"))
    d >= options.pca_components || throw(FisherFitError("too_few_patch_pixels"))
    all(isfinite, X) && all(isfinite, amplitudes) || throw(FisherFitError("nonfinite_training_input"))
    basis = options.patch_projection == "affine_disk" ? _affine_disk_basis(fisher_grid(options)) : nothing
    if basis !== nothing
        size(basis, 1) == d || throw(DimensionMismatch("affine projection requires the configured Fisher disk"))
        X = Matrix{Float64}(X) - (X * basis) * transpose(basis)
        any(row -> sum(abs, row) > options.projection_zero_l1, eachrow(X)) ||
            throw(FisherFitError("zero_affine_patch_mass"))
    end
    Xc = Matrix{Float64}(X) .- mean(X; dims=1)
    all(isfinite, Xc) || throw(FisherFitError("nonfinite_centered_input"))
    decomposition = try
        svd(Xc; full=false)
    catch error
        error isa LAPACKException || rethrow()
        throw(FisherFitError("pca_failure"))
    end
    !isempty(decomposition.S) && decomposition.S[1] > 0 || throw(FisherFitError("no_patch_variation"))
    V = Matrix(decomposition.V[:, 1:options.pca_components])
    Z = Xc * V
    gmm = _fit_gmm(Z, options)
    linear = _fisher_linear(Xc, Z, V, Float64.(amplitudes), gmm.assignments,
                            options.noise_regularization)
    mid = options.score_center == "training_mean" ? vec(mean(X; dims=1)) : linear.mid
    return FisherModel(linear.w_p, mid, linear.g0, linear.g1, linear.amplitude_means, gmm, basis)
end

"Same projection at training and scoring, with the explicitly selected score origin."
function score(x::AbstractVector, model::FisherModel)
    prepared = model.projection_basis === nothing ? x : _project_disk(x, model.projection_basis)
    return dot(prepared - model.mid, model.w_p)
end
maxmirror_score(x::AbstractVector, model::FisherModel, grid::FisherGrid) =
    max(score(x, model), score(flip_u_disk(x, grid), model))

struct FisherScore
    file::String
    lobe::Int
    score::Float64                 # NaN represents invalid, never serialized as a number
    invalid_reason::String
end

"""
Fit even and odd lobe folds independently. Each row uses only the OTHER fold's
PCA, GMM, amplitude mapping, covariance, and Fisher vector. No full-data fit is
built. Missing/degenerate training folds invalidate their held-out rows without
removing keys or inventing a cluster. Finite nonconverged EM fits are retained,
like sklearn's max_iter behavior, with a warning rather than silent fallback.
"""
function cv_scores(patches::PatchTable, options::FisherOptions; training_mask=nothing)
    n = length(patches.keys)
    eligible = validate_training_mask(options.training_support, training_mask, n)
    size(patches.X, 1) == n && length(patches.amplitudes) == n && length(patches.invalid_reasons) == n ||
        throw(DimensionMismatch("patch table rows must match keys"))
    grid = fisher_grid(options)
    grid.disk_indices == patches.grid.disk_indices && grid.coords == patches.grid.coords ||
        throw(ArgumentError("patch table grid does not match Fisher config"))
    reasons = copy(patches.invalid_reasons)
    basis = options.patch_projection == "affine_disk" ? _affine_disk_basis(grid) : nothing
    for i in 1:n
        if isempty(reasons[i]) && !(all(isfinite, patches.X[i, :]) && isfinite(patches.amplitudes[i]))
            reasons[i] = "nonfinite_patch_input"
        elseif isempty(reasons[i]) && basis !== nothing &&
               sum(abs, _project_disk(patches.X[i, :], basis)) <= options.projection_zero_l1
            reasons[i] = "zero_affine_patch_mass"
        end
    end
    margins = fill(NaN, n)
    valid = isempty.(reasons)
    for parity in (0, 1)
        train = [i for i in 1:n if valid[i] && eligible[i] && mod(patches.keys[i][2], 2) == parity]
        held = [i for i in 1:n if valid[i] && mod(patches.keys[i][2], 2) != parity]
        isempty(held) && continue
        fold = parity == 0 ? "even" : "odd"
        model = try
            fit_fisher(patches.X[train, :], patches.amplitudes[train], options)
        catch error
            error isa FisherFitError || rethrow()
            for i in held
                reasons[i] = "training_$fold:" * error.reason
            end
            continue
        end
        model.gmm.converged || @warn "Fisher GMM reached configured iteration limit" fold iterations=model.gmm.iterations
        for i in held
            value = maxmirror_score(patches.X[i, :], model, grid)
            if isfinite(value)
                margins[i] = value
            else
                reasons[i] = "nonfinite_fisher_score"
            end
        end
    end
    return [FisherScore(key[1], key[2], margins[i], reasons[i]) for (i, key) in enumerate(patches.keys)]
end

"Write one CV margin table. No chemical prediction and no full-data artifact."
function write_scores(path::AbstractString, rows::AbstractVector{FisherScore})
    ispath(path) && throw(ArgumentError("output already exists: $path"))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "file\tlobe\tscore\tinvalid_reason")
        for row in rows
            if isempty(row.invalid_reason) && isfinite(row.score)
                @printf(io, "%s\t%d\t%.6f\t\n", row.file, row.lobe, row.score)
            else
                reason = isempty(row.invalid_reason) ? "nonfinite_fisher_score" : row.invalid_reason
                println(io, row.file, '\t', row.lobe, "\tNA\t", reason)
            end
        end
    end
    return path
end

end # module

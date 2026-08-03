#!/usr/bin/env julia

using Printf
using Random
using Statistics
using LinearAlgebra
using STMSXMIO
using GaussianFit2D

include(joinpath(@__DIR__, "lib", "joint_proxy", "inference_cli.jl"))

using Main.JointProxyRegistry: ProxyEnsemble, load_registry
using Main.JointProxyTypePosterior: TypePosteriorLobeEvidence, infer_type_posterior

const Inference = Main.JointProxyInferenceCLI
const REAL_LOBE_FIELDS = [:file, :n, :lobe, :patch_kind, :ablation, :p0, :p1,
    :margin_p1, :argmax_type, :confidence, :residual_corr, :geometric_mass, :stm_dft_v1_mass]
const REAL_SUMMARY_FIELDS = [:file, :status, :n, :patch_kind, :ablation, :lobe_count,
    :argmax0, :argmax1, :ties, :mean_p1, :mean_confidence, :log_evidence,
    :geometric_mass, :stm_dft_v1_mass]
const VIEW_TRANSFER_FIELDS = [:file, :status, :n, :lobe, :n_pairs, :corr, :slope,
    :intercept, :affine_nrmse, :best_dt_steps, :best_du_steps, :best_dt_nm,
    :best_du_nm, :best_corr, :corr_gain]

Base.@kwdef struct DiagnosticOptions
    config::String
    synthetic_out::Union{Nothing,String} = nothing
    data_dir::Union{Nothing,String} = nothing
    files_from::Union{Nothing,String} = nothing
    real_out::Union{Nothing,String} = nothing
    real_summary::Union{Nothing,String} = nothing
    view_out::Union{Nothing,String} = nothing
    view_max_shift_steps::Union{Nothing,Int} = nothing
    chunk::Union{Nothing,Tuple{Int,Int}} = nothing
    seed::Int = 20260721
    replicates::Int = 5
    noise_sigmas::Vector{Float64} = [0.0, 0.10]
end

function parse_options(args)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--config", "--synthetic-out", "--data-dir", "--files-from",
            "--real-out", "--real-summary", "--chunk", "--seed", "--replicates",
            "--noise-sigmas", "--view-out", "--view-max-shift-steps") ||
            error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        values[flag] = String(args[i + 1])
        i += 2
    end
    haskey(values, "--config") || error("--config is required")
    synthetic_out = get(values, "--synthetic-out", nothing)
    real_out = get(values, "--real-out", nothing)
    synthetic_out === nothing && real_out === nothing && error("request --synthetic-out and/or --real-out")
    real_requested = real_out !== nothing
    for flag in ("--data-dir", "--files-from", "--real-summary")
        real_requested && !haskey(values, flag) && error("$flag is required with --real-out")
    end
    view_out = get(values, "--view-out", nothing)
    view_shift = haskey(values, "--view-max-shift-steps") ?
        parse(Int, values["--view-max-shift-steps"]) : nothing
    view_out !== nothing && !real_requested && error("--view-out requires --real-out")
    view_out !== nothing && view_shift === nothing &&
        error("--view-max-shift-steps is required with --view-out")
    view_out === nothing && view_shift !== nothing &&
        error("--view-max-shift-steps requires --view-out")
    view_shift !== nothing && view_shift < 0 && error("--view-max-shift-steps must be non-negative")
    chunk = if haskey(values, "--chunk")
        parts = split(values["--chunk"], '/')
        length(parts) == 2 || error("--chunk must be i/n")
        parsed = (parse(Int, parts[1]), parse(Int, parts[2]))
        1 <= parsed[1] <= parsed[2] || error("invalid --chunk")
        parsed
    else
        nothing
    end
    noise = parse.(Float64, split(get(values, "--noise-sigmas", "0,0.10"), ','))
    all(x -> isfinite(x) && x >= 0, noise) || error("noise sigmas must be finite and non-negative")
    return DiagnosticOptions(
        config=values["--config"], synthetic_out=synthetic_out,
        data_dir=get(values, "--data-dir", nothing), files_from=get(values, "--files-from", nothing),
        real_out=real_out, real_summary=get(values, "--real-summary", nothing), chunk=chunk,
        view_out=view_out, view_max_shift_steps=view_shift,
        seed=parse(Int, get(values, "--seed", "20260721")),
        replicates=parse(Int, get(values, "--replicates", "5")), noise_sigmas=noise,
    )
end

function family_ensembles(registry::ProxyEnsemble)
    families = sort(unique(entry.source.family for entry in registry.entries))
    out = Pair{String,ProxyEnsemble}[]
    for family in families
        entries = filter(entry -> entry.source.family == family, registry.entries)
        push!(out, family => ProxyEnsemble(entries, registry.grid_half_nm, registry.grid_step_nm,
            registry.grid_n, registry.npix, registry.payload_sha256, ["diagnostic_family=$family"]))
    end
    push!(out, "combined" => registry)
    return out
end

function family_mass(result, family::String)
    pair = findfirst(item -> first(item) == family, result.proxy_family_sensitivity)
    return pair === nothing ? 0.0 : last(result.proxy_family_sensitivity[pair])
end

function backbone_residual(patch::AbstractVector{<:Real}, points,
    sigma_parallel_nm::Real, sigma_perp_nm::Real)
    length(patch) == length(points) || error("patch/nuisance-grid length mismatch")
    sigma_parallel_nm > 0 && sigma_perp_nm > 0 || error("nuisance sigmas must be positive")
    basis = [exp(-0.5 * ((point.t / sigma_parallel_nm)^2 +
        (point.u / sigma_perp_nm)^2)) for point in points]
    finite = [i for i in eachindex(patch) if isfinite(patch[i]) && isfinite(basis[i])]
    if length(finite) < 2
        values = [Float64(value) for value in patch if isfinite(value)]
        center = isempty(values) ? 0.0 : mean(values)
        return [isfinite(value) ? Float64(value) - center : 0.0 for value in patch]
    end
    design = hcat(ones(length(finite)), basis[finite])
    coefficients = design \ Float64.(patch[finite])
    return [isfinite(patch[i]) ? Float64(patch[i]) -
        (coefficients[1] + coefficients[2] * basis[i]) : 0.0 for i in eachindex(patch)]
end

function matched_type_posterior(candidate, patch_views, ensemble; rho=NaN, effective_factor=nothing)
    posterior = Main.JointProxyTypePosterior
    scale = posterior.effective_view_scale(rho, effective_factor)
    n = length(candidate.lobes)
    n == length(patch_views) || error("matched posterior lobe count mismatch")
    state_logweights = Float64[]
    state_marginals = Matrix{Float64}[]
    state_families = String[]
    for direction in 0:1, phase in 0:1, mirror in 0:1, entry in ensemble.entries
        unary = zeros(Float64, n, 2)
        for (i, (lobe, views)) in enumerate(zip(candidate.lobes, patch_views))
            parity = posterior._parity_for_lobe(i, n, direction, phase)
            observed = Dict(view => backbone_residual(patch, candidate.patch_tu,
                lobe.sigma_parallel_nm, lobe.sigma_perp_nm) for (view, patch) in views)
            for typ in 0:1
                template = posterior._template_for(entry, typ, parity, mirror)
                matched_template = backbone_residual(template, candidate.patch_tu,
                    lobe.sigma_parallel_nm, lobe.sigma_perp_nm)
                unary[i, typ + 1] = scale * sum(
                    posterior.log_ncc_score(patch, matched_template) for patch in values(observed)) +
                    log(0.5)
            end
        end
        all(isfinite, unary) || return nothing
        chain = posterior.binary_chain_forward_backward(unary)
        log_prior = -3 * log(2.0) + log(entry.source.weight)
        push!(state_logweights, log_prior + chain.log_evidence)
        push!(state_marginals, chain.state_marginals)
        push!(state_families, entry.source.family)
    end
    weights, log_evidence = posterior._normalize_logweights(state_logweights)
    marginals = zeros(Float64, n, 2)
    masses = Dict{String,Float64}()
    for (weight, lobe_marginals, family) in zip(weights, state_marginals, state_families)
        marginals .+= weight .* lobe_marginals
        masses[family] = get(masses, family, 0.0) + weight
    end
    return (; lobe_marginals=marginals, log_evidence,
        proxy_family_sensitivity=sort(collect(masses); by=first))
end

function _nonnegative_backbone_mold_fit(y::Vector{Float64}, gaussian::Vector{Float64},
    mold::Vector{Float64})
    yc = y .- mean(y)
    gc = gaussian .- mean(gaussian)
    mc = mold .- mean(mold)
    candidates = Tuple{Float64,Float64}[(0.0, 0.0)]
    gg, mm = dot(gc, gc), dot(mc, mc)
    gg > 0 && push!(candidates, (max(0.0, dot(gc, yc) / gg), 0.0))
    mm > 0 && push!(candidates, (0.0, max(0.0, dot(mc, yc) / mm)))
    gm = dot(gc, mc)
    determinant = gg * mm - gm^2
    if determinant > eps(Float64) * max(gg * mm, floatmin(Float64))
        gy, my = dot(gc, yc), dot(mc, yc)
        amplitude = (mm * gy - gm * my) / determinant
        coefficient = (gg * my - gm * gy) / determinant
        amplitude >= 0 && coefficient >= 0 && push!(candidates, (amplitude, coefficient))
    end
    best = first(candidates)
    best_rss = Inf
    for (amplitude, coefficient) in candidates
        residual = yc .- amplitude .* gc .- coefficient .* mc
        rss = dot(residual, residual)
        if rss < best_rss
            best = (amplitude, coefficient)
            best_rss = rss
        end
    end
    amplitude, coefficient = best
    background = mean(y .- amplitude .* gaussian .- coefficient .* mold)
    return (; background, amplitude, coefficient, rss=best_rss)
end

function joint_generative_score(patch::AbstractVector{<:Real}, template::AbstractVector{<:Real},
    points, sigma_parallel_nm::Real, sigma_perp_nm::Real)
    length(patch) == length(template) == length(points) ||
        error("patch/template/nuisance-grid length mismatch")
    sigma_parallel_nm > 0 && sigma_perp_nm > 0 || error("nuisance sigmas must be positive")
    gaussian = [exp(-0.5 * ((point.t / sigma_parallel_nm)^2 +
        (point.u / sigma_perp_nm)^2)) for point in points]
    finite = [i for i in eachindex(patch) if
        isfinite(patch[i]) && isfinite(template[i]) && isfinite(gaussian[i])]
    length(finite) >= 3 || return NaN
    y = Float64.(patch[finite])
    g = Float64.(gaussian[finite])
    m = Float64.(template[finite])
    null_fit = _nonnegative_backbone_mold_fit(y, g, zeros(length(y)))
    model_fit = _nonnegative_backbone_mold_fit(y, g, m)
    rss_model = min(model_fit.rss, null_fit.rss)
    centered_energy = sum(abs2, y .- mean(y))
    floor_rss = eps(Float64) * max(centered_energy, floatmin(Float64))
    null_fit.rss > floor_rss || return 0.0
    return clamp((null_fit.rss - rss_model) / null_fit.rss, 0.0, 1.0)
end

function generative_type_posterior(candidate, patch_views, ensemble; rho=NaN, effective_factor=nothing)
    posterior = Main.JointProxyTypePosterior
    scale = posterior.effective_view_scale(rho, effective_factor)
    n = length(candidate.lobes)
    n > 0 || error("generative posterior candidate has no lobes")
    n == length(patch_views) || error("generative posterior lobe count mismatch")
    labels = sort(collect(keys(first(patch_views))))
    isempty(labels) && error("generative posterior lobe has no view patches")
    all(views -> sort(collect(keys(views))) == labels, patch_views) ||
        error("generative posterior view-label mismatch")
    state_logweights = Float64[]
    state_marginals = Matrix{Float64}[]
    state_families = String[]
    for direction in 0:1, phase in 0:1, mirror in 0:1, entry in ensemble.entries
        unary = zeros(Float64, n, 2)
        for (i, (lobe, views)) in enumerate(zip(candidate.lobes, patch_views))
            parity = posterior._parity_for_lobe(i, n, direction, phase)
            for typ in 0:1
                template = posterior._template_for(entry, typ, parity, mirror)
                unary[i, typ + 1] = scale * sum(joint_generative_score(
                    views[label], template, candidate.patch_tu,
                    lobe.sigma_parallel_nm, lobe.sigma_perp_nm,
                ) for label in labels) + log(0.5)
            end
        end
        all(isfinite, unary) || return nothing
        chain = posterior.binary_chain_forward_backward(unary)
        log_prior = -3 * log(2.0) + log(entry.source.weight)
        push!(state_logweights, log_prior + chain.log_evidence)
        push!(state_marginals, chain.state_marginals)
        push!(state_families, entry.source.family)
    end
    weights, log_evidence = posterior._normalize_logweights(state_logweights)
    marginals = zeros(Float64, n, 2)
    masses = Dict{String,Float64}()
    for (weight, lobe_marginals, family) in zip(weights, state_marginals, state_families)
        marginals .+= weight .* lobe_marginals
        masses[family] = get(masses, family, 0.0) + weight
    end
    return (; lobe_marginals=marginals, log_evidence,
        proxy_family_sensitivity=sort(collect(masses); by=first))
end

function diagnostic_posterior(lobes, ensemble; rho=NaN, effective_factor=nothing)
    try
        return infer_type_posterior(lobes, ensemble; rho, effective_factor)
    catch error
        message = sprint(showerror, error)
        error isa ArgumentError && occursin("non-finite unary evidence", message) && return nothing
        rethrow()
    end
end

function synthetic_rows(registry::ProxyEnsemble, options::DiagnosticOptions)
    rng = MersenneTwister(options.seed)
    ablations = family_ensembles(registry)
    coords = collect(-registry.grid_half_nm:registry.grid_step_nm:registry.grid_half_nm)
    patch_tu = [(; t, u) for u in coords for t in coords]
    synthetic_candidate = (; lobes=[(; sigma_parallel_nm=0.20, sigma_perp_nm=0.14)], patch_tu)
    gaussian = [exp(-0.5 * ((point.t / 0.20)^2 + (point.u / 0.14)^2)) for point in patch_tu]
    rows = NamedTuple[]
    for entry in sort(registry.entries; by=e -> (e.source.family, e.source.name))
        for template in sort(entry.templates; by=t -> (t.type, t.parity, t.mirror))
            for sigma in options.noise_sigmas, replicate in 1:options.replicates
                noise = sigma .* randn(rng, length(template.pixels))
                patch = template.pixels .+ noise
                generative_patch = 0.3 .+ 2.0 .* gaussian .+ template.pixels .+ noise
                lobes = [TypePosteriorLobeEvidence(Dict("fwd" => patch))]
                for (ablation, ensemble) in ablations
                    result = diagnostic_posterior(lobes, ensemble; effective_factor=1.0)
                    result === nothing && error("template self-control produced non-finite evidence")
                    p0, p1 = result.lobe_marginals[1, 1], result.lobe_marginals[1, 2]
                    predicted = p1 > p0 ? 1 : p0 > p1 ? 0 : -1
                    ptrue = template.type == 0 ? p0 : p1
                    pother = template.type == 0 ? p1 : p0
                    push!(rows, (; injected_family=entry.source.family, injected_source=entry.source.name,
                        injected_type=template.type, parity=template.parity, mirror=template.mirror,
                        noise_sigma=sigma, replicate, ablation, p0, p1, predicted,
                        correct=predicted == template.type, true_margin=ptrue - pother,
                        log_evidence=result.log_evidence,
                        geometric_mass=family_mass(result, "geometric"),
                        stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
                    matched = matched_type_posterior(synthetic_candidate,
                        [Dict("fwd" => patch)], ensemble; effective_factor=1.0)
                    matched === nothing && error("matched template self-control produced non-finite evidence")
                    mp0, mp1 = matched.lobe_marginals[1, 1], matched.lobe_marginals[1, 2]
                    mpredicted = mp1 > mp0 ? 1 : mp0 > mp1 ? 0 : -1
                    mptrue = template.type == 0 ? mp0 : mp1
                    mpother = template.type == 0 ? mp1 : mp0
                    push!(rows, (; injected_family=entry.source.family, injected_source=entry.source.name,
                        injected_type=template.type, parity=template.parity, mirror=template.mirror,
                        noise_sigma=sigma, replicate, ablation="$(ablation)_matched", p0=mp0, p1=mp1,
                        predicted=mpredicted, correct=mpredicted == template.type,
                        true_margin=mptrue - mpother, log_evidence=matched.log_evidence,
                        geometric_mass=family_mass(matched, "geometric"),
                        stm_dft_v1_mass=family_mass(matched, "stm_dft_v1")))
                    generative = generative_type_posterior(synthetic_candidate,
                        [Dict("fwd" => generative_patch)], ensemble; effective_factor=1.0)
                    generative === nothing && error("generative template self-control produced non-finite evidence")
                    gp0, gp1 = generative.lobe_marginals[1, 1], generative.lobe_marginals[1, 2]
                    gpredicted = gp1 > gp0 ? 1 : gp0 > gp1 ? 0 : -1
                    gptrue = template.type == 0 ? gp0 : gp1
                    gpother = template.type == 0 ? gp1 : gp0
                    push!(rows, (; injected_family=entry.source.family, injected_source=entry.source.name,
                        injected_type=template.type, parity=template.parity, mirror=template.mirror,
                        noise_sigma=sigma, replicate, ablation="$(ablation)_generative", p0=gp0, p1=gp1,
                        predicted=gpredicted, correct=gpredicted == template.type,
                        true_margin=gptrue - gpother, log_evidence=generative.log_evidence,
                        geometric_mass=family_mass(generative, "geometric"),
                        stm_dft_v1_mass=family_mass(generative, "stm_dft_v1")))
                end
            end
        end
    end
    return rows
end

function read_files(options::DiagnosticOptions)
    files = [strip(line) for line in readlines(something(options.files_from)) if !isempty(strip(line))]
    all(file -> basename(file) == file && endswith(lowercase(file), ".sxm"), files) ||
        error("--files-from must contain bare .sxm basenames")
    length(unique(files)) == length(files) || error("--files-from contains duplicate basenames")
    if options.chunk !== nothing
        index, total = options.chunk
        files = files[index:total:end]
    end
    return files
end

function diagnostic_file_seed(base_seed::Integer, file::AbstractString)
    state = UInt64(base_seed)
    for byte in codeunits(file)
        state = state * UInt64(0x00000100000001b3) ⊻ UInt64(byte)
    end
    return Int(mod(state, UInt64(typemax(Int))))
end

function diagnostic_chain_config(seed::Integer)
    base = Inference._default_chain_config()
    names = fieldnames(typeof(base))
    values = NamedTuple{names}(ntuple(i -> getfield(base, names[i]), length(names)))
    return GaussianFit2D.ChainSweepConfig(;
        merge(values, (; skip_global=true, rng_seed=Int(seed)))...)
end

function image_lobe_evidence(candidate, views)
    first_lobe, last_lobe = first(candidate.lobes), last(candidate.lobes)
    dx = last_lobe.x_nm - first_lobe.x_nm
    dy = last_lobe.y_nm - first_lobe.y_nm
    norm_axis = hypot(dx, dy)
    ax, ay = norm_axis > eps(Float64) ? (dx / norm_axis, dy / norm_axis) : (1.0, 0.0)
    px, py = -ay, ax
    return [TypePosteriorLobeEvidence(Dict(
        view.label => [Main.JointProxyCandidateViews._interp(
            view.xs, view.ys, view.z,
            lobe.x_nm + offset.t * ax + offset.u * px,
            lobe.y_nm + offset.t * ay + offset.u * py,
        ) for offset in candidate.patch_tu]
        for view in views if view.present,
    )) for lobe in candidate.lobes]
end

function generative_patch_sets(lobes::Vector{TypePosteriorLobeEvidence})
    isempty(lobes) && return Pair{String,Vector{Dict{String,Vector{Float64}}}}[]
    patches = [lobe.view_patches for lobe in lobes]
    labels = sort(collect(keys(first(patches))))
    all(patch -> sort(collect(keys(patch))) == labels, patches) ||
        error("generative image view-label mismatch")
    sets = Pair{String,Vector{Dict{String,Vector{Float64}}}}["image" => patches]
    if length(labels) > 1
        for label in labels
            push!(sets, "image_$(label)" => [Dict(label => patch[label]) for patch in patches])
        end
    end
    return sets
end

generative_view_factor(patches, combined_factor::Real) =
    all(length(patch) == 1 for patch in patches) ? 1.0 : Float64(combined_factor)

function _affine_patch_metrics(fwd::AbstractVector{<:Real}, bwd::AbstractVector{<:Real})
    length(fwd) == length(bwd) || error("fwd/bwd patch length mismatch")
    finite = [i for i in eachindex(fwd) if isfinite(fwd[i]) && isfinite(bwd[i])]
    n_pairs = length(finite)
    n_pairs >= 3 || return (; n_pairs, corr=NaN, slope=NaN, intercept=NaN,
        affine_nrmse=NaN)
    x = Float64.(fwd[finite])
    y = Float64.(bwd[finite])
    xc = x .- mean(x)
    yc = y .- mean(y)
    xx = dot(xc, xc)
    yy = dot(yc, yc)
    denom = sqrt(xx * yy)
    corr = denom > eps(Float64) ? dot(xc, yc) / denom : NaN
    slope = xx > eps(Float64) ? dot(xc, yc) / xx : NaN
    intercept = isfinite(slope) ? mean(y) - slope * mean(x) : NaN
    rss = isfinite(slope) ? sum(abs2, y .- (slope .* x .+ intercept)) : NaN
    affine_nrmse = yy > eps(Float64) && isfinite(rss) ? sqrt(rss / yy) : NaN
    return (; n_pairs, corr, slope, intercept, affine_nrmse)
end

function _shifted_patch_pairs(fwd::AbstractVector{<:Real}, bwd::AbstractVector{<:Real},
    grid_n::Int, dt::Int, du::Int)
    length(fwd) == length(bwd) == grid_n^2 || error("patch/grid size mismatch")
    shifted_fwd = Float64[]
    shifted_bwd = Float64[]
    for iu in 1:grid_n, it in 1:grid_n
        jt, ju = it + dt, iu + du
        1 <= jt <= grid_n && 1 <= ju <= grid_n || continue
        push!(shifted_fwd, Float64(fwd[(iu - 1) * grid_n + it]))
        push!(shifted_bwd, Float64(bwd[(ju - 1) * grid_n + jt]))
    end
    return shifted_fwd, shifted_bwd
end

function view_transfer_metrics(fwd::AbstractVector{<:Real}, bwd::AbstractVector{<:Real},
    grid_n::Int, max_shift_steps::Int)
    max_shift_steps >= 0 || error("max shift must be non-negative")
    max_shift_steps < grid_n || error("max shift must be smaller than the patch grid")
    base = _affine_patch_metrics(fwd, bwd)
    best_dt, best_du = 0, 0
    best_corr = base.corr
    best_abs = isfinite(best_corr) ? abs(best_corr) : -Inf
    best_distance = 0
    for dt in -max_shift_steps:max_shift_steps, du in -max_shift_steps:max_shift_steps
        dt == 0 && du == 0 && continue
        shifted_fwd, shifted_bwd = _shifted_patch_pairs(fwd, bwd, grid_n, dt, du)
        shifted = _affine_patch_metrics(shifted_fwd, shifted_bwd)
        isfinite(shifted.corr) || continue
        candidate_abs = abs(shifted.corr)
        candidate_distance = abs(dt) + abs(du)
        if candidate_abs > best_abs + 1e-12 ||
            (abs(candidate_abs - best_abs) <= 1e-12 && candidate_distance < best_distance)
            best_dt, best_du = dt, du
            best_corr = shifted.corr
            best_abs = candidate_abs
            best_distance = candidate_distance
        end
    end
    corr_gain = isfinite(base.corr) && isfinite(best_corr) ?
        abs(best_corr) - abs(base.corr) : NaN
    return merge(base, (; best_dt_steps=best_dt, best_du_steps=best_du,
        best_corr, corr_gain))
end

function real_rows(registry::ProxyEnsemble, options::DiagnosticOptions)
    ablations = family_ensembles(registry)
    lobe_rows = NamedTuple[]
    summaries = NamedTuple[]
    view_rows = NamedTuple[]
    for file in read_files(options)
        path = joinpath(something(options.data_dir), file)
        image = STMSXMIO.read_sxm(path)
        pattern = Inference._default_pattern_config(path, dirname(something(options.real_out)))
        chain = diagnostic_chain_config(diagnostic_file_seed(options.seed, file))
        results, _best, context = redirect_stdout(devnull) do
            GaussianFit2D.chain_gaussian_sweep(image, pattern, chain)
        end
        views = Inference._build_views(image, pattern, Main.JointProxyCandidateViews.build_view_data)
        report = Main.JointProxyCandidateViews.extract_candidate_views(results, context, chain, views;
            patch_half_nm=registry.grid_half_nm, patch_step_nm=registry.grid_step_nm)
        if isempty(report.candidates)
            push!(summaries, (; file, status="no_candidates", n=0, patch_kind="none",
                ablation="none", lobe_count=0,
                argmax0=0, argmax1=0, ties=0, mean_p1=NaN, mean_confidence=NaN,
                log_evidence=NaN, geometric_mass=NaN, stm_dft_v1_mass=NaN))
            continue
        end
        candidate = argmin(item -> item.joint_gcv, report.candidates)
        matched_ablations = family_ensembles(registry)
        patch_sets = Pair{String,Vector{TypePosteriorLobeEvidence}}[
            "image" => image_lobe_evidence(candidate, views),
            "residual" => [TypePosteriorLobeEvidence(lobe.residual_patches) for lobe in candidate.lobes],
            "raw" => [TypePosteriorLobeEvidence(lobe.raw_patches) for lobe in candidate.lobes],
            "neg_residual" => [TypePosteriorLobeEvidence(Dict(
                view => -patch for (view, patch) in lobe.residual_patches)) for lobe in candidate.lobes],
        ]
        if options.view_out !== nothing
            image_lobes = first(patch_sets).second
            max_shift = something(options.view_max_shift_steps)
            for (lobe, evidence) in zip(candidate.lobes, image_lobes)
                if haskey(evidence.view_patches, "fwd") && haskey(evidence.view_patches, "bwd")
                    metrics = view_transfer_metrics(evidence.view_patches["fwd"],
                        evidence.view_patches["bwd"], registry.grid_n, max_shift)
                    push!(view_rows, (; file, status="ok", n=candidate.n, lobe=lobe.index,
                        metrics.n_pairs, metrics.corr, metrics.slope, metrics.intercept,
                        metrics.affine_nrmse, metrics.best_dt_steps, metrics.best_du_steps,
                        best_dt_nm=metrics.best_dt_steps * registry.grid_step_nm,
                        best_du_nm=metrics.best_du_steps * registry.grid_step_nm,
                        metrics.best_corr, metrics.corr_gain))
                else
                    push!(view_rows, (; file, status="bwd_missing", n=candidate.n,
                        lobe=lobe.index, n_pairs=0, corr=NaN, slope=NaN, intercept=NaN,
                        affine_nrmse=NaN, best_dt_steps=0, best_du_steps=0,
                        best_dt_nm=0.0, best_du_nm=0.0, best_corr=NaN, corr_gain=NaN))
                end
            end
        end
        for (patch_kind, lobes) in patch_sets, (ablation, ensemble) in ablations
                result = diagnostic_posterior(lobes, ensemble; rho=candidate.residual_corr,
                    effective_factor=candidate.effective_view_factor)
                if result === nothing
                    push!(summaries, (; file, status="nonfinite_evidence", n=candidate.n,
                        patch_kind, ablation, lobe_count=length(lobes), argmax0=0, argmax1=0,
                        ties=0, mean_p1=NaN, mean_confidence=NaN, log_evidence=NaN,
                        geometric_mass=NaN, stm_dft_v1_mass=NaN))
                    continue
                end
                predicted = Int[]
                for (lobe, probabilities) in zip(candidate.lobes, eachrow(result.lobe_marginals))
                    p0, p1 = probabilities
                    argmax_type = p1 > p0 ? 1 : p0 > p1 ? 0 : -1
                    push!(predicted, argmax_type)
                    push!(lobe_rows, (; file, n=candidate.n, lobe=lobe.index, patch_kind,
                        ablation, p0, p1, margin_p1=p1 - p0, argmax_type,
                        confidence=max(p0, p1), residual_corr=candidate.residual_corr,
                        geometric_mass=family_mass(result, "geometric"),
                        stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
                end
                push!(summaries, (; file, status="ok", n=candidate.n, patch_kind, ablation,
                    lobe_count=length(predicted), argmax0=count(==(0), predicted),
                    argmax1=count(==(1), predicted), ties=count(==(-1), predicted),
                    mean_p1=mean(result.lobe_marginals[:, 2]),
                    mean_confidence=mean(maximum(row) for row in eachrow(result.lobe_marginals)),
                    log_evidence=result.log_evidence,
                    geometric_mass=family_mass(result, "geometric"),
                    stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
            end
        for (patch_kind, lobes) in patch_sets, (ablation, ensemble) in matched_ablations
            patches = [lobe.view_patches for lobe in lobes]
            result = matched_type_posterior(candidate, patches, ensemble;
                rho=candidate.residual_corr, effective_factor=candidate.effective_view_factor)
            matched_name = "$(ablation)_matched"
            if result === nothing
                push!(summaries, (; file, status="nonfinite_evidence", n=candidate.n,
                    patch_kind, ablation=matched_name, lobe_count=length(lobes), argmax0=0,
                    argmax1=0, ties=0, mean_p1=NaN, mean_confidence=NaN,
                    log_evidence=NaN, geometric_mass=NaN, stm_dft_v1_mass=NaN))
                continue
            end
            predicted = Int[]
            for (lobe, probabilities) in zip(candidate.lobes, eachrow(result.lobe_marginals))
                p0, p1 = probabilities
                argmax_type = p1 > p0 ? 1 : p0 > p1 ? 0 : -1
                push!(predicted, argmax_type)
                push!(lobe_rows, (; file, n=candidate.n, lobe=lobe.index, patch_kind,
                    ablation=matched_name, p0, p1, margin_p1=p1 - p0, argmax_type,
                    confidence=max(p0, p1), residual_corr=candidate.residual_corr,
                    geometric_mass=family_mass(result, "geometric"),
                    stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
            end
            push!(summaries, (; file, status="ok", n=candidate.n, patch_kind,
                ablation=matched_name, lobe_count=length(predicted),
                argmax0=count(==(0), predicted), argmax1=count(==(1), predicted),
                ties=count(==(-1), predicted), mean_p1=mean(result.lobe_marginals[:, 2]),
                mean_confidence=mean(maximum(row) for row in eachrow(result.lobe_marginals)),
                log_evidence=result.log_evidence,
                geometric_mass=family_mass(result, "geometric"),
                stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
        end
        for (patch_kind, image_patches) in generative_patch_sets(first(patch_sets).second),
            (ablation, ensemble) in ablations
            view_factor = generative_view_factor(image_patches, candidate.effective_view_factor)
            result = generative_type_posterior(candidate, image_patches, ensemble;
                rho=candidate.residual_corr, effective_factor=view_factor)
            generative_name = "$(ablation)_generative"
            if result === nothing
                push!(summaries, (; file, status="nonfinite_evidence", n=candidate.n,
                    patch_kind, ablation=generative_name, lobe_count=length(image_patches),
                    argmax0=0, argmax1=0, ties=0, mean_p1=NaN, mean_confidence=NaN,
                    log_evidence=NaN, geometric_mass=NaN, stm_dft_v1_mass=NaN))
                continue
            end
            predicted = Int[]
            for (lobe, probabilities) in zip(candidate.lobes, eachrow(result.lobe_marginals))
                p0, p1 = probabilities
                argmax_type = p1 > p0 ? 1 : p0 > p1 ? 0 : -1
                push!(predicted, argmax_type)
                push!(lobe_rows, (; file, n=candidate.n, lobe=lobe.index, patch_kind,
                    ablation=generative_name, p0, p1, margin_p1=p1 - p0, argmax_type,
                    confidence=max(p0, p1), residual_corr=candidate.residual_corr,
                    geometric_mass=family_mass(result, "geometric"),
                    stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
            end
            push!(summaries, (; file, status="ok", n=candidate.n, patch_kind,
                ablation=generative_name, lobe_count=length(predicted),
                argmax0=count(==(0), predicted), argmax1=count(==(1), predicted),
                ties=count(==(-1), predicted), mean_p1=mean(result.lobe_marginals[:, 2]),
                mean_confidence=mean(maximum(row) for row in eachrow(result.lobe_marginals)),
                log_evidence=result.log_evidence,
                geometric_mass=family_mass(result, "geometric"),
                stm_dft_v1_mass=family_mass(result, "stm_dft_v1")))
        end
        end
    return lobe_rows, summaries, view_rows
end

format_tsv(value::Bool) = value ? "true" : "false"
format_tsv(value::Real) = isfinite(Float64(value)) ? @sprintf("%.15g", Float64(value)) : "NaN"
format_tsv(value) = string(value)

function write_tsv(path::AbstractString, rows; fields=nothing)
    fields === nothing && isempty(rows) && error("fields are required for an empty diagnostic table: $path")
    mkpath(dirname(path))
    selected_fields = fields === nothing ? collect(keys(first(rows))) : collect(fields)
    open(path, "w") do io
        println(io, join(selected_fields, '\t'))
        for row in rows
            println(io, join((format_tsv(getproperty(row, field)) for field in selected_fields), '\t'))
        end
    end
end

function run_diagnostic(options::DiagnosticOptions)
    registry = load_registry(options.config)
    options.replicates > 0 || error("--replicates must be positive")
    options.synthetic_out === nothing || write_tsv(options.synthetic_out, synthetic_rows(registry, options))
    if options.real_out !== nothing
        lobes, summaries, views = real_rows(registry, options)
        write_tsv(options.real_out, lobes; fields=REAL_LOBE_FIELDS)
        write_tsv(something(options.real_summary), summaries; fields=REAL_SUMMARY_FIELDS)
        options.view_out === nothing || write_tsv(options.view_out, views; fields=VIEW_TRANSFER_FIELDS)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_diagnostic(parse_options(ARGS))
end

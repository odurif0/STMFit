"""
Core computation module for multi-Gaussian fitting.

Pure math — no plotting, no CLI. Dependencies: Optimization.jl + NLopt, LsqFit, Statistics, Printf.

Uses typed structs (FitConfig, FitResult, FitMetrics) from types.jl.
"""

# ===========================================================================
# 1. Data loading
# ===========================================================================

function load_data(filepath::String)
    """Load x and intensity from a whitespace-delimited file.

    Handles 2 or 3 columns. Deduplicates, filters NaN/Inf, and sorts by x.
    """
    data = readdlm(filepath)
    ncols = size(data, 2)
    if ncols == 3
        x = vec(data[:, 2])
        y = vec(data[:, 3])
    elseif ncols == 2
        x = vec(data[:, 1])
        y = vec(data[:, 2])
    else
        error("Expected 2 or 3 columns, got $ncols")
    end

    # Filter NaN/Inf rows (corrupted data)
    finite_mask = isfinite.(x) .& isfinite.(y)
    if any(.!finite_mask)
        n_bad = count(.!finite_mask)
        x = x[finite_mask]
        y = y[finite_mask]
        @warn "Dropped $n_bad rows with NaN/Inf values."
    end

    # Deduplicate x, keeping first-occurrence y values
    ux = unique(x)
    idx = Int[findfirst(==(xi), x) for xi in ux]
    x = ux
    y = y[idx]

    # Sort
    p = sortperm(x)
    x = x[p]
    y = y[p]

    return x, y
end

# ===========================================================================
# 2. Gaussian model (delta-reparameterized)
# ===========================================================================
#
# params = [y0, A_0, mu_0, sigma_0, A_1, delta_1, sigma_1, ...]
# n_params = 1 + 3 * n_peaks

function _effective_min_spacing(cfg::FitConfig)
    """Return the hard adjacent-center spacing lower bound.

    In addition to the explicit physical/chemical `min_spacing`, enforce the
    minimum spacing needed so two widest allowed peaks overlap by at most
    `max_overlap` at each other's center.
    """
    return effective_spacing_min(cfg.min_spacing, cfg.max_spacing,
                                 sigma_from_fwhm(cfg.fwhm_max), cfg.max_overlap)
end

function _get_sigma(params::AbstractVector{<:Real}, i::Int)
    return params[4 + 3 * i]
end

function _get_amplitude(params::AbstractVector{<:Real}, i::Int; use_log::Bool=false)
    v = params[2 + 3 * i]
    return use_log ? exp(v) : v
end

function _params_to_centers(params::AbstractVector{<:Real}, n_peaks::Int)
    centers = [params[3]]  # mu_0
    for i in 1:(n_peaks - 1)
        push!(centers, centers[end] + params[3 + 3 * i])
    end
    return centers
end

function _center_errors(perr::AbstractVector{<:Real}, n_peaks::Int)
    """Propagate parameter errors to absolute center errors.

    center_i = mu_0 + sum(delta_1..delta_i)
    σ²(center_i) = σ²(mu_0) + Σ_j σ²(delta_j)
    perr is the full error vector (with y0=0 at index 1).
    """
    isempty(perr) && return Float64[]
    # Build cumulative variance: start with mu_0, add deltas in order
    cvar = perr[3]^2  # mu_0 variance (index 3 in full params)
    cerr = Float64[sqrt(cvar)]
    for i in 1:(n_peaks - 1)
        # delta_i is at index 3 + 3*i in full params
        di = 3 + 3 * i
        if length(perr) >= di
            cvar += perr[di]^2
        end
        push!(cerr, sqrt(cvar))
    end
    return cerr
end

function multi_gaussian(x::AbstractVector{<:Real}, params::AbstractVector{<:Real};
                        n_peaks::Int=0, asymmetric_edges::Bool=false,
                        use_log_amplitude::Bool=false,
                        y_buf::Union{Vector{Float64},Nothing}=nothing)
    """y0 + sum_i A_i * exp(-0.5 * ((x - mu_i) / sigma_i)^2)

    When asymmetric_edges=true, the first and last peaks use split Gaussians
    with different sigma on each side of the center.

    If y_buf is provided (same length as x), it is reused in-place instead
    of allocating a new output vector on each call.  The caller must not
    rely on y_buf's contents across calls.
    """
    y0 = params[1]
    nx = length(x)

    # Build centers — small allocation (n_peaks <= ~30), amortized cost negligible
    centers = _params_to_centers(params, n_peaks)

    # Pre-allocated output buffer or fresh allocation
    y_model = y_buf !== nothing ? y_buf : similar(x)
    fill!(y_model, y0)

    if asymmetric_edges && n_peaks >= 2
        sigma_outer_left  = params[end - 1]
        sigma_outer_right = params[end]

        if n_peaks > 2
            for i in 1:(n_peaks - 2)
                @inbounds mu    = centers[i+1]
                @inbounds sigma = params[4 + 3 * i]
                @inbounds amp_raw = params[2 + 3 * i]
                amp = use_log_amplitude ? exp(amp_raw) : amp_raw
                @. y_model += amp * exp(-0.5 * ((x - mu) / sigma)^2)
            end
        end

        # Edge peak 0: split
        @inbounds amp0_raw = params[2]
        amp0    = use_log_amplitude ? exp(amp0_raw) : amp0_raw
        c0     = centers[1]
        @inbounds sig0_in = params[4]
        @. y_model += amp0 * exp(-0.5 * ((x - c0) / ifelse(x < c0, sigma_outer_left, sig0_in))^2)

        # Edge peak N-1: split
        @inbounds ampN_raw = params[2 + 3 * (n_peaks - 1)]
        ampN    = use_log_amplitude ? exp(ampN_raw) : ampN_raw
        cN     = centers[end]
        @inbounds sigN_in = params[4 + 3 * (n_peaks - 1)]
        @. y_model += ampN * exp(-0.5 * ((x - cN) / ifelse(x < cN, sigN_in, sigma_outer_right))^2)

        return y_model
    else
        # Symmetric: broadcast-fused
        for i in 0:(n_peaks - 1)
            @inbounds mu    = centers[i+1]
            @inbounds sigma = params[4 + 3 * i]
            @inbounds amp_raw = params[2 + 3 * i]
            amp = use_log_amplitude ? exp(amp_raw) : amp_raw
            @. y_model += amp * exp(-0.5 * ((x - mu) / sigma)^2)
        end
        return y_model
    end
end

function predict_fit(x::AbstractVector{<:Real}, result::FitResult, cfg::FitConfig)
    """Evaluate a fitted model on an arbitrary x grid."""
    asymmetric = cfg.asymmetric_edges && result.n_peaks >= 2
    return multi_gaussian(x, result.popt; n_peaks=result.n_peaks,
                          asymmetric_edges=asymmetric,
                          use_log_amplitude=cfg.use_log_amplitude)
end

# ===========================================================================
# 3. Bounds construction
# ===========================================================================

function build_bounds(n_peaks::Int, x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig)
    """Build parameter bounds for the delta-parameterization.

    Baseline is fixed at 0 (data is offset so y.min() = 0 before fitting).

    Returns (lower, upper).
    """
    x_min, x_max = extrema(x)
    fwhm_min = cfg.fwhm_min
    fwhm_max = cfg.fwhm_max
    max_sp   = cfg.max_spacing
    sigma_min = fwhm_min / FWHM_TO_SIGMA
    sigma_max = fwhm_max / FWHM_TO_SIGMA
    margin    = cfg.edge_sigma_min * sigma_max
    amp_max   = cfg.amplitude_max !== nothing ? cfg.amplitude_max : maximum(y)
    amp_min   = cfg.amplitude_min !== nothing ? cfg.amplitude_min :
                cfg.amplitude_min_fraction * amp_max
    amp_low   = cfg.use_log_amplitude ? log(max(amp_min, 1e-10)) : amp_min
    amp_high  = cfg.use_log_amplitude ? log(max(amp_max, 1e-10)) : amp_max
    min_sp    = _effective_min_spacing(cfg)

    lower = Float64[]
    upper = Float64[]

    for i in 0:(n_peaks - 1)
        # Amplitude
        push!(lower, amp_low)
        push!(upper, amp_high)

        if i == 0
            # mu_0: first peak center
            push!(lower, x_min + margin)
            mu0_upper = x_max - margin - (n_peaks - 1) * min_sp
            if mu0_upper < lower[end]
                error("""
                    Impossible bounds for n_peaks=$n_peaks: mu_0 upper ($(@sprintf("%.4f", mu0_upper)))
                    is below lower ($(@sprintf("%.4f", lower[end])), x_min+margin).
                    Reduce effective min_spacing ($(@sprintf("%.4f", min_sp))) or n_peaks to fit the x range
                    [$(@sprintf("%.4f", x_min)), $(@sprintf("%.4f", x_max))].
                    Required range: $(@sprintf("%.4f", 2*margin + (n_peaks-1)*min_sp)),
                    available: $(@sprintf("%.4f", x_max - x_min)).
                    """)
            end
            push!(upper, mu0_upper)
        else
            # delta_i: spacing from peak i-1 to peak i, bounded by [min_sp, max_sp]
            push!(lower, min_sp)
            push!(upper, max_sp)
        end

        # Per-peak sigma
        push!(lower, sigma_min)
        push!(upper, sigma_max)
    end

    return lower, upper
end

function _build_edge_sigma_bounds(cfg::FitConfig)
    fwhm_min = cfg.fwhm_min
    fwhm_max = cfg.fwhm_max
    sigma_min = fwhm_min / FWHM_TO_SIGMA
    sigma_max = fwhm_max / FWHM_TO_SIGMA
    sigma_outer_max = 1.5 * sigma_max
    return [sigma_min, sigma_min], [sigma_outer_max, sigma_outer_max]
end

function _build_initial_guess(n_peaks::Int, x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig)
    lower, upper = build_bounds(n_peaks, x, y, cfg)
    return [(lo + hi) / 2.0 for (lo, hi) in zip(lower, upper)]
end

# ===========================================================================
# 4. Warm-start helpers
# ===========================================================================
# ===========================================================================
# 5. Fitting: NLopt global + Levenberg-Marquardt local
# ===========================================================================

function _make_objective_function(x, y, n_peaks, asymmetric_edges, use_log_amplitude=false;
                                   kappa_max=25.0, kappa_weight=1.0)
    """Create the RSS objective (1-arg, compatible with NLopt via wrapper).

    Baseline is always fixed at 0 (data is pre-offset so y.min() = 0).
    Uses pre-allocated buffers to avoid allocations on every call.
    When kappa_max > 0 and n_peaks > 1, adds a progressive κ penalty.
    """
    n_params_inner = 3 * n_peaks + (asymmetric_edges && n_peaks >= 2 ? 2 : 0)
    full_buf = zeros(n_params_inner + 1)  # +1 for y0=0
    y_buf = similar(x)  # pre-alloc output buffer for multi_gaussian
    # y0 at index 1 is always 0
    function objective(params_vec::AbstractVector{<:Real})
        @. full_buf[2:end] = params_vec
        residuals = y - multi_gaussian(x, full_buf; n_peaks=n_peaks,
                                        asymmetric_edges=asymmetric_edges,
                                        use_log_amplitude=use_log_amplitude,
                                        y_buf=y_buf)
        rss = sum(abs2, residuals)
        # Condition-number penalty for adjacent overlap
        if kappa_max > 0 && n_peaks > 1
            deltas = [params_vec[3 + 3*j] for j in 1:(n_peaks-1)]
            sigmas = [params_vec[4 + 3*k] for k in 0:(n_peaks-1)]
            κ = STMFitCore.adjacent_kappa_max(deltas, sigmas)
            rss *= (1.0 + STMFitCore.kappa_penalty(κ; kappa_max, weight=kappa_weight))
        end
        return rss
    end
    return objective
end

function fit_model(x::Vector{Float64}, y::Vector{Float64}, n_peaks::Int, cfg::FitConfig;
                     p0::Union{Vector{Float64},Nothing}=nothing,
                     maxtime_override::Union{Float64,Nothing}=nothing)
    """Two-stage fit: NLopt global + Levenberg-Marquardt (LsqFit).

    Baseline is fixed at 0 (data is pre-offset so y.min() = 0).
    If p0 is provided, it is used as the initial guess instead of the midpoint.
    maxtime_override (seconds) overrides cfg.global_maxtime for this call only.

    Returns a FitResult.
    """
    asymmetric = cfg.asymmetric_edges && n_peaks >= 2
    lo_edge = Float64[]
    hi_edge = Float64[]
    if asymmetric
        lo_edge, hi_edge = _build_edge_sigma_bounds(cfg)
    end

    lower, upper = build_bounds(n_peaks, x, y, cfg)

    # Append edge sigma bounds if asymmetric
    if asymmetric
        append!(lower, lo_edge)
        append!(upper, hi_edge)
    end

    if p0 === nothing
        p0 = _build_initial_guess(n_peaks, x, y, cfg)
        if asymmetric
            for (lo, hi) in zip(lo_edge, hi_edge)
                push!(p0, (lo + hi) / 2.0)
            end
        end
    else
        # Clip initial guess p0 to bounds
        p0 = copy(p0)
        for idx in eachindex(p0)
            if p0[idx] < lower[idx]
                p0[idx] = lower[idx] + 1e-6
            elseif p0[idx] > upper[idx]
                p0[idx] = upper[idx] - 1e-6
            end
        end
        if asymmetric && length(p0) < length(lower)
            for (lo, hi) in zip(lo_edge, hi_edge)
                push!(p0, (lo + hi) / 2.0)
            end
        end
    end

    objective_1arg = _make_objective_function(x, y, n_peaks, asymmetric, cfg.use_log_amplitude;
                                               kappa_max=cfg.kappa_max, kappa_weight=cfg.kappa_weight)
    # Optimization.jl expects f(u, p); the second argument p is unused
    objective_opt(u, _) = objective_1arg(u)

    # ---- Stage 1: NLopt global optimization ----
    algo_sym = cfg.nlopt_algorithm
    nlop = OptimizationNLopt.NLopt.Opt(algo_sym, length(lower))
    nlop.xtol_rel = cfg.global_tol
    nlop.ftol_rel = cfg.global_tol
    nlop.maxtime  = something(maxtime_override, cfg.global_maxtime)

    prob = OptimizationProblem(objective_opt, p0; lb=lower, ub=upper)

    sol = try
        solve(prob, nlop; maxiters=cfg.global_maxiter)
    catch e
        @warn "NLopt global optimization failed: $e"
        return FitResult(
            popt=Float64[], pcov=nothing, perr=Float64[],
            y_fit=Float64[], success=false,
            warnings=["Global optimization failed: $e"],
        )
    end

    global_params = sol.u

    # ---- Stage 2: Levenberg-Marquardt refinement (LsqFit) ----
    warnings_list = String[]

    # Build model function: always prepend y0 = 0 (reuse pre-allocated buffers)
    n_inner = length(lower)
    lm_buf = zeros(n_inner + 1)  # +1 for y0=0
    lm_y_buf = similar(x)        # output buffer for multi_gaussian
    model_func_inner = (xdata, p) -> begin
        @. lm_buf[2:end] = p
        return multi_gaussian(xdata, lm_buf; n_peaks=n_peaks,
                              asymmetric_edges=asymmetric,
                              use_log_amplitude=cfg.use_log_amplitude,
                              y_buf=lm_y_buf)
    end

    popt_inner = Float64[]
    perr_inner = Float64[]
    pcov = nothing
    success = true

    try
        fit = curve_fit(model_func_inner, x, y, global_params;
                        lower=lower, upper=upper,
                        maxIter=cfg.curve_fit_maxfev,
                        show_trace=false,
                        autodiff=:finite)
        popt_inner = fit.param
        # Estimate covariance; may fail for ill-conditioned problems
        try
            pcov = estimate_covar(fit)
            if any(isinf.(pcov)) || any(isnan.(pcov))
                push!(warnings_list,
                      "Covariance matrix ill-conditioned. Uncertainties unreliable.")
                perr_inner = fill(NaN, length(popt_inner))
            else
                perr_inner = sqrt.(diag(pcov))
            end
        catch cov_e
            push!(warnings_list,
                  "Covariance estimation failed ($(typeof(cov_e))). Uncertainties unreliable.")
            perr_inner = fill(NaN, length(popt_inner))
        end
    catch e
        push!(warnings_list, "curve_fit failed: $e.")
        popt_inner = global_params
        perr_inner = fill(NaN, length(popt_inner))
        success = false
    end

    # Reconstruct full params vector (y0 = 0 prepended)
    popt = vcat([0.0], popt_inner)
    perr = vcat([0.0], perr_inner)

    y_fit = multi_gaussian(x, popt; n_peaks=n_peaks, asymmetric_edges=asymmetric,
                            use_log_amplitude=cfg.use_log_amplitude)

    # Bound-at-limit warnings (full vector: y0=0 at index 1)
    x_unit = cfg.x_unit
    full_lower = vcat([0.0], lower)
    full_upper = vcat([0.0], upper)
    for (idx, (val, lo, hi)) in enumerate(zip(popt, full_lower, full_upper))
        # y0 is fixed at 0 by construction — never a real bound issue
        idx == 1 && continue
        if abs(val - lo) < 1e-6 || abs(val - hi) < 1e-6
            label = _param_label(idx - 1, n_peaks=n_peaks, asymmetric_edges=asymmetric)
            if occursin("sigma", lowercase(label))
                fwhm = FWHM_TO_SIGMA * val
                push!(warnings_list,
                      "Parameter $label is at its bound " *
                      "($(@sprintf("%.4f", val)) $x_unit, FWHM=$(@sprintf("%.4f", fwhm)) $x_unit).")
            else
                push!(warnings_list,
                      "Parameter $label is at its bound ($(@sprintf("%.6f", val))).")
            end
        end
    end

    return FitResult(
        popt=popt, pcov=pcov, perr=perr,
        y_fit=y_fit, success=success,
        warnings=warnings_list,
    )
end

function _param_label(idx::Int;
                      n_peaks::Int=0,
                      asymmetric_edges::Bool=false)
    """Human-readable parameter name from full-vector index.

    The idx always refers to the *full* parameter vector (with y0=0 prepended).
    When asymmetric_edges, the last 2 params are edge outer sigmas.
    """
    if idx == 0
        return "y0 (baseline)"
    end

    if asymmetric_edges && n_peaks >= 2
        n_main = 1 + 3 * n_peaks
        if idx == n_main
            return "sigma_outer_left (peak 1)"
        elseif idx == n_main + 1
            return "sigma_outer_right (peak N)"
        end
    end

    peak = (idx - 1) ÷ 3
    slot = (idx - 1) % 3
    if slot == 0
        return "A_$(peak + 1)"
    elseif slot == 1
        return peak == 0 ? "mu_0" : "delta_$(peak + 1)"
    else
        return "sigma_$(peak + 1)"
    end
end

# ===========================================================================
# 6. Model comparison & metrics
# ===========================================================================

function compute_metrics(y::Vector{Float64}, y_fit::Vector{Float64}, n_params::Int; student_nu::Float64=4.0, noise_estimate::Float64=NaN)
    n = length(y)
    rss = sum((y - y_fit).^2)
    tss = sum((y .- mean(y)).^2)
    dof = max(1, n - n_params)
    # Gaussian BIC (legacy)
    bic_g = n * log(rss / n) + n_params * log(n)
    # Student-t BIC (consistent with 2D)
    resid = y .- y_fit
    if isfinite(noise_estimate)
        noise_est = noise_estimate
    else
        noise_est = max(std(resid) * 0.1, median(abs.(resid)) * 1.4826, 1e-12)
    end
    nu = student_nu
    student_nll = sum(0.5 * (nu + 1) .* log1p.((resid ./ noise_est).^2 ./ nu))
    bic_s = 2 * student_nll + n_params * log(n)
    # AICc using Student-t NLL (consistent with 2D)
    aicc_s = 2 * student_nll + 2 * n_params + (2 * n_params * (n_params + 1)) / max(n - n_params - 1, 1)
    return FitMetrics(
        bic_g,                                      # bic (Gaussian)
        bic_s,                                      # student_bic
        n * log(rss / n) + 2 * n_params,             # aic
        aicc_s,                                      # aicc (Student-t NLL)
        1.0 - rss / tss,                             # r_squared
        rss / dof / max(noise_est^2, 1e-12),            # chi2_red (normalized by noise)
        dof,
        rss,
        n_params,
    )
end

function _determine_n_peaks_range(x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig;
                                  verbose::Bool=true)
    """Determine min/max n_peaks and center_n for the sweep."""
    x_range = maximum(x) - minimum(x)
    fwhm_max = cfg.fwhm_max
    sigma_max = fwhm_max / FWHM_TO_SIGMA
    edge_min = cfg.edge_sigma_min * sigma_max
    edge_max = cfg.edge_sigma_max * sigma_max
    max_sp   = cfg.max_spacing
    x_unit   = cfg.x_unit

    usable_conservative = x_range - 2 * edge_max
    usable_optimistic   = x_range - 2 * edge_min

    min_sp = _effective_min_spacing(cfg)
    min_n = max(2, Int(floor(usable_conservative / max_sp)))
    max_n = max(min_n, Int(floor(usable_optimistic / min_sp)))

    # Center estimate: use same margin as build_bounds (edge_sigma_min),
    # not the conservative edge_sigma_max used for n_peaks range above.
    mean_sp = (min_sp + max_sp) / 2.0
    edge_margin = cfg.edge_sigma_min * sigma_max
    usable = x_range - 2 * edge_margin
    center_n = Int(round(usable / mean_sp))
    center_n = max(min_n, min(max_n, center_n))

    if verbose
        @printf("  x range = %.2f %s\n", x_range, x_unit)
        @printf("  edge_margin: min=%.3f (bounds), max=%.3f (n_peaks estimate)\n", edge_min, edge_max)
        @printf("  spacing in [%.2f, %.2f] %s (configured min=%.2f, max_overlap=%.2f)\n",
                min_sp, max_sp, x_unit, cfg.min_spacing, cfg.max_overlap)
        @printf("  usable range: %.2f (conservative) to %.2f (optimistic) %s\n", usable_conservative, usable_optimistic, x_unit)
        @printf("  n_peaks range: %d to %d\n", min_n, max_n)
    end
    return min_n, max_n, center_n
end

function _fit_one(n_peaks::Int, x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig;
                 p0::Union{Vector{Float64},Nothing}=nothing,
                 maxtime_override::Union{Float64,Nothing}=nothing)
    n_extra = (cfg.asymmetric_edges && n_peaks >= 2) ? 2 : 0
    n_params = 1 + 3 * n_peaks + n_extra
    fit = fit_model(x, y, n_peaks, cfg; p0=p0, maxtime_override=maxtime_override)
    if !fit.success || isempty(fit.popt)
        return nothing
    end
    m = compute_metrics(y, fit.y_fit, n_params; student_nu=cfg.student_nu, noise_estimate=cfg.noise_estimate)
    # Post-fit: max adjacent condition number (from full params with y0 prepended)
    kappa_val = 1.0
    if n_peaks > 1
        deltas = [fit.popt[3 + 3*i] for i in 1:(n_peaks-1)]   # delta_i at full index 3+3i
        sigmas = [fit.popt[4 + 3*k] for k in 0:(n_peaks-1)]   # sigma_k at full index 4+3k
        kappa_val = STMFitCore.adjacent_kappa_max(deltas, sigmas)
    end
    result = FitResult(
        n_peaks=n_peaks,
        popt=fit.popt, pcov=fit.pcov, perr=fit.perr,
        y_fit=fit.y_fit, success=fit.success, warnings=fit.warnings,
        bic=m.bic, student_bic=m.student_bic, aic=m.aic, aicc=m.aicc, r_squared=m.r_squared,
        chi2_red=m.chi2_red, dof=m.dof, rss=m.rss, n_params=m.n_params,
        kappa_max_adj=kappa_val,
    )
    return result
end

function _sweep_direction(x, y, cfg::FitConfig, seed_result, start_n, extreme_n, step, direction,
                           expand_warm::Bool, early_stop_patience, early_stop_dbic;
                           progress_callback=nothing, verbose::Bool=true)
    """Sweep one direction (right=expand, left=shrink). Returns Dict{Int,FitResult}."""
    _bic(r) = cfg.use_student_bic ? r.student_bic : r.bic
    results = Dict{Int,FitResult}()
    best_bic = seed_result !== nothing ? _bic(seed_result) : Inf
    streak = 0
    n = start_n
    last_result = seed_result
    done = (step > 0 ? n > extreme_n : n < extreme_n)

    while !done
        if streak >= early_stop_patience
            if verbose
                @printf("    -> %s skip n=%d (BIC already increased %dx)\n",
                        direction, n, early_stop_patience)
            end
            done = true
        else
            if verbose
                @printf("  Fitting n_peaks = %d (%d params) [%s]...\n", n, 1 + 3 * n, direction)
            end
            r = try
                _fit_one(n, x, y, cfg; p0=nothing, maxtime_override=cfg.global_maxtime / 2)
            catch e
                if verbose
                    println("    FAILED n=$n: $e")
                end
                nothing
            end
            if r !== nothing
                if verbose
                    tag = ""
                    if !isempty(r.warnings)
                        tag = "  [$(length(r.warnings)) warning(s)]"
                    end
                    bic_label = cfg.use_student_bic ? "sBIC" : "BIC"
                    @printf("    DONE n=%d: %s = %.1f  |  R^2 = %.4f%s\n",
                            n, bic_label, _bic(r), r.r_squared, tag)
                end
                r.popt_inner = r.popt[2:end]
                results[n] = r
                last_result = r
                # Progress callback
                if progress_callback !== nothing
                    progress_callback(direction, n, length(results), abs(extreme_n - start_n) + 2, r)
                end
                # Check early stop
                if _bic(r) < best_bic
                    best_bic = _bic(r)
                    streak = 0
                elseif (_bic(r) - best_bic) > early_stop_dbic
                    if verbose
                        @printf("    -> %s early stop at n=%d (d%s = %.0f > %.0f)\n",
                                direction, n, bic_label, _bic(r) - best_bic, early_stop_dbic)
                    end
                    done = true
                else
                    streak += 1
                    if streak >= early_stop_patience
                        if verbose
                            @printf("    -> %s early stop at n=%d (%s increased %dx)\n",
                                    direction, n, bic_label, early_stop_patience)
                        end
                        done = true
                    end
                end
            end
        end
        n += step
        if step > 0 ? n > extreme_n : n < extreme_n
            done = true
        end
    end
    return results
end

function run_model_comparison(x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig;
                              progress_callback=nothing)
    """Bidirectional sweep from center n_peaks with early stop.
    Left and right directions run in parallel via Threads.@spawn."""
    early_stop_patience = cfg.early_stop_patience
    early_stop_dbic     = cfg.early_stop_dbic
    bic_threshold       = cfg.bic_competition_threshold

    min_n, max_n, center_n = _determine_n_peaks_range(x, y, cfg)

    println()
    @printf("--- Bidirectional sweep: n_peaks from %d to %d (start=%d, patience=%d, dBIC_stop=%.0f) ---\n\n",
            min_n, max_n, center_n, early_stop_patience, early_stop_dbic)

    results = Dict{Int,FitResult}()

    function _try_fit(n::Int; p0::Union{Vector{Float64},Nothing}=nothing)
        try
            r = _fit_one(n, x, y, cfg; p0=p0)
            if r !== nothing
                tag = ""
                if !isempty(r.warnings)
                    tag = "  [$(length(r.warnings)) warning(s)]"
                end
                bic_val = cfg.use_student_bic ? r.student_bic : r.bic
                bic_label = cfg.use_student_bic ? "sBIC" : "BIC"
                @printf("    DONE n=%d: %s = %.1f  |  R^2 = %.4f%s\n",
                        n, bic_label, bic_val, r.r_squared, tag)
            else
                println("    FAILED n=$n: fit returned no valid result")
            end
            return r
        catch e
            println("    FAILED n=$n: $e")
            return nothing
        end
    end

    # ---- 1. Fit center (midpoint initialization) ----
    @printf("  Fitting n_peaks = %d (%d params) [center]...\n", center_n, 1 + 3 * center_n)
    r = _try_fit(center_n)
    if r !== nothing
        results[center_n] = r
        if progress_callback !== nothing
            progress_callback("center", center_n, 1, max_n - min_n + 1, r)
        end
    end

    # ---- 2. Expand bidirectionally in parallel ----
    if r !== nothing
        t_right = Threads.@spawn _sweep_direction(
            x, y, cfg, r, center_n + 1, max_n, +1, "right", true,
            early_stop_patience, early_stop_dbic;
            progress_callback=progress_callback)
        t_left  = Threads.@spawn _sweep_direction(
            x, y, cfg, r, center_n - 1, min_n, -1, "left", false,
            early_stop_patience, early_stop_dbic;
            progress_callback=progress_callback)
        for (k, v) in fetch(t_right); results[k] = v; end
        for (k, v) in fetch(t_left);  results[k] = v; end
    end

    # Sort
    result_keys = sort(collect(keys(results)))
    results_list = [results[k] for k in result_keys]

    if isempty(results_list)
        return results_list
    end

    # Tag competitive
    use_sbic = cfg.use_student_bic
    best_val = use_sbic ? minimum(r.student_bic for r in results_list) : minimum(r.bic for r in results_list)
    for r in results_list
        delta = (use_sbic ? r.student_bic : r.bic) - best_val
        r.competitive = delta <= bic_threshold
        r.delta_bic = delta
    end

    n_total = length(results_list)
    n_comp = count(r -> r.competitive, results_list)
    @printf("\n  %d models fitted, %d competitive (delta_BIC <= %.0f)\n",
            n_total, n_comp, bic_threshold)

    return results_list
end

# ===========================================================================
# 7. Application API
# ===========================================================================

function _print_input_summary(x, y, cfg::FitConfig)
    x_unit = cfg.x_unit
    max_sp = cfg.max_spacing
    fwhm_min = cfg.fwhm_min
    fwhm_max = cfg.fwhm_max
    sigma_max = fwhm_max / FWHM_TO_SIGMA
    edge_min = cfg.edge_sigma_min * sigma_max
    edge_max = cfg.edge_sigma_max * sigma_max

    @printf("  %d data points, x in [%.3f, %.3f] %s\n", length(x), minimum(x), maximum(x), x_unit)
    @printf("  Intensity range: [%.3f, %.3f]\n", minimum(y), maximum(y))
    min_sp_eff = _effective_min_spacing(cfg)
    @printf("  spacing in [%.2f, %.2f] %s (configured min=%.2f, max_overlap=%.2f), FWHM in [%.2f, %.2f] %s\n",
            min_sp_eff, max_sp, x_unit, cfg.min_spacing, cfg.max_overlap, fwhm_min, fwhm_max, x_unit)
    @printf("  edge_margin: %.3f (%.1f sigma) to %.3f (%.1f sigma) %s\n",
            edge_min, cfg.edge_sigma_min, edge_max, cfg.edge_sigma_max, x_unit)
end

function run_fit(cfg::FitConfig=DEFAULT_CONFIG; save_cache::Bool=true, verbose::Bool=true,
                 progress_callback=nothing)
    """Load data, run model comparison, return FitRunResult."""
    if verbose
        println("Loading data from: $(cfg.filepath)")
    end
    x, y = load_data(cfg.filepath)
    return run_fit(x, y, cfg; save_cache=save_cache, verbose=verbose,
                   progress_callback=progress_callback)
end

function run_fit(x::Vector{Float64}, y::Vector{Float64}, cfg::FitConfig;
                 save_cache::Bool=true, verbose::Bool=true,
                 progress_callback=nothing)
    """Run model comparison on in-memory x/y data (no file I/O)."""
    if cfg.offset_to_zero
        y_offset = minimum(y)
        if y_offset != 0.0
            y = y .- y_offset
            if verbose
                @printf("  Data offset by %.4f (min set to 0)\n", -y_offset)
            end
        end
    end
    if verbose
        _print_input_summary(x, y, cfg)
    end

    all_results = run_model_comparison(x, y, cfg; progress_callback=progress_callback)
    if isempty(all_results)
        if verbose
            @warn "No models could be fitted. Check your data and parameter settings."
        end
        return FitRunResult(x, y, all_results, cfg)
    end

    fr = FitRunResult(x, y, all_results, cfg)
    if save_cache
        fr.cache_file = save_results(cfg.filepath, x, y, all_results, cfg)
    end
    return fr
end

# ===========================================================================
# 8. Serialization (JLD2)
# ===========================================================================

function save_results(filepath::String, x, y, all_results, cfg::FitConfig)
    out_dir = cfg.output_dir
    cache_file = _cache_path(filepath, out_dir)
    dir = dirname(cache_file)
    isempty(dir) || mkpath(dir)
    chash = _config_hash(cfg)
    JLD2.jldsave(cache_file; x=x, y=y, all_results=all_results, cfg=cfg,
                 mgf_version=MGF_VERSION, config_hash=chash)
    println("  Results cached to: $cache_file")
    return cache_file
end

function _cache_path(filepath::String, output_dir)
    base = splitext(basename(filepath))[1]
    filename = base * "_cache.jld2"
    if output_dir !== nothing && !isempty(output_dir)
        mkpath(output_dir)
        return joinpath(output_dir, filename)
    end
    cd = _cache_dir()
    mkpath(cd)
    return joinpath(cd, filename)
end

function load_results(cache_file::String)
    data = JLD2.load(cache_file)
    println("  Loaded cached results from: $cache_file")

    # Version check
    cached_ver = get(data, "mgf_version", nothing)
    if cached_ver !== nothing && cached_ver != MGF_VERSION
        @warn "Cache version ($cached_ver) differs from current version ($MGF_VERSION). Results may be stale."
    end

    # Config hash check
    cached_hash = get(data, "config_hash", nothing)
    cfg = data["cfg"]
    if cached_hash !== nothing && cfg isa FitConfig
        current_hash = _config_hash(cfg)
        if cached_hash != current_hash
            @warn "Cache config hash mismatch. The cached config differs from the current config. Results may be misleading."
        end
    end

    n_models = length(data["all_results"])
    n_list = [r.n_peaks for r in data["all_results"]]
    println("  $n_models models: n_peaks = $n_list")
    return data["x"], data["y"], data["all_results"], cfg
end

# ===========================================================================
# 9. Export results
# ===========================================================================

function export_results(x, y, all_results, cfg::FitConfig)
    """Export all models (parameters + residuals) to a text file."""
    x_unit = cfg.x_unit
    best = cfg.use_student_bic ? argmin(r -> r.student_bic, all_results) : argmin(r -> r.bic, all_results)
    out_file = output_path(cfg, "_results.txt")
    bic_label = cfg.use_student_bic ? "sBIC" : "BIC"

    open(out_file, "w") do f
        write(f, "# Multi-Gaussian fit results\n")
        @printf(f, "# Data: %s  (%d points, x in [%.3f, %.3f] %s)\n",
                cfg.filepath, length(x), minimum(x), maximum(x), x_unit)
        @printf(f, "# spacing_min_configured = %.6f %s\n", cfg.min_spacing, x_unit)
        @printf(f, "# spacing_min_effective = %.6f %s\n", _effective_min_spacing(cfg), x_unit)
        @printf(f, "# spacing_min_effective_source = max(min_spacing, sqrt(-2log(max_overlap))*sigma_max)\n")
        @printf(f, "# max_overlap = %.6f\n", cfg.max_overlap)
        @printf(f, "# kappa_max = %.1f\n", cfg.kappa_max)
        write(f, "#\n")
        @printf(f, "# %8s  %10s  %6s  %10s  %10s  %8s  %10s  %10s  %6s  %4s  %4s\n",
                "n_peaks", bic_label, "d"*bic_label, "AICc", "AIC", "R2", "chi2_red", "RSS", "kappa", "BEST", "COMP")
        for r in all_results
            bval = cfg.use_student_bic ? r.student_bic : r.bic
            tag = r === best ? "   *" : ""
            comp = r.competitive && r !== best ? "  +" : ""
            @printf(f, "# %8d  %10.1f  %6.1f  %10.1f  %10.1f  %8.4f  %10.4f  %10.4f  %6.1f%s%s\n",
                    r.n_peaks, bval, r.delta_bic, r.aicc, r.aic,
                    r.r_squared, r.chi2_red, r.rss, r.kappa_max_adj, tag, comp)
        end
        bic_thresh = cfg.bic_competition_threshold
        @printf(f, "# * = BEST, + = competitive (dBIC <= %.0f)\n#\n\n", bic_thresh)

        for r in all_results
            n_peaks = r.n_peaks
            popt = r.popt
            perr = r.perr
            y_fit = predict_fit(x, r, cfg)
            is_best = r === best
            tag = is_best ? " (BEST)" : ""

            write(f, repeat("=", 72), "\n")
            bval = cfg.use_student_bic ? r.student_bic : r.bic
            @printf(f, "# n_peaks = %d%s  |  %s = %.1f  |  R^2 = %.4f\n",
                    n_peaks, tag, bic_label, bval, r.r_squared)
            write(f, repeat("=", 72), "\n")

            centers = _params_to_centers(popt, n_peaks)
            cerr = _center_errors(perr, n_peaks)

            write(f, "peak,center($x_unit),center_err,amplitude,amplitude_err," *
                   "sigma($x_unit),sigma_err,fwhm($x_unit),fwhm_err\n")
            @printf(f, "baseline,%.6f,%.6f,,,,,,\n", popt[1], perr[1])

            for i in 0:(n_peaks - 1)
                A = _get_amplitude(popt, i; use_log=cfg.use_log_amplitude)
                sigma = _get_sigma(popt, i)
                A_e = _get_amplitude(perr, i; use_log=cfg.use_log_amplitude)
                s_e = _get_sigma(perr, i)
                fwhm = FWHM_TO_SIGMA * sigma
                fwhm_e = FWHM_TO_SIGMA * s_e
                ce = get(cerr, i + 1, NaN)
                @printf(f, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                        i + 1, centers[i+1], ce, A, A_e, sigma, s_e, fwhm, fwhm_e)
            end

            if n_peaks > 1
                spacings = [popt[3 + 3 * i] for i in 1:(n_peaks - 1)]
                parts = join([@sprintf("%.6f", s) for s in spacings], ",")
                write(f, "spacing($x_unit),$parts\n")
                @printf(f, "kappa_max_adj,%.4f\n", r.kappa_max_adj)
            end

            write(f, "\n# Residuals: x($x_unit), y_data, y_fit, residual\n")
            residuals = y - y_fit
            for (xi, yi, fi, ri) in zip(x, y, y_fit, residuals)
                @printf(f, "%.6f,%.6f,%.6f,%.6f\n", xi, yi, fi, ri)
            end
            write(f, "\n\n")
        end
    end

    println("  All results saved to: $out_file")
    return out_file
end

# ===========================================================================
# 10. Console summary
# ===========================================================================

function print_summary(all_results, best_res, cfg::FitConfig)
    x_unit = cfg.x_unit
    n_peaks = best_res.n_peaks
    popt = best_res.popt
    perr = best_res.perr
    bic_label = cfg.use_student_bic ? "sBIC" : "BIC"

    println("\n" * repeat("=", 90))
    println("MODEL COMPARISON")
    println(repeat("=", 90))
    @printf("  %8s  %10s  %6s  %10s  %10s  %8s  %10s  %10s  %9s\n",
            "n_peaks", bic_label, "d"*bic_label, "AICc", "AIC", "R^2", "chi2_red", "RSS", "n_params")
    @printf("  %8s  %10s  %6s  %10s  %10s  %8s  %10s  %10s  %9s\n",
            repeat("-", 6), repeat("-", 10), repeat("-", 6), repeat("-", 10),
            repeat("-", 10), repeat("-", 8), repeat("-", 10), repeat("-", 10), repeat("-", 9))

    for r in all_results
        bval = cfg.use_student_bic ? r.student_bic : r.bic
        tag = r === best_res ? " <-- BEST" : ""
        comp = r.competitive && r !== best_res ? " *" : ""
        @printf("  %8d  %10.1f  %6.1f  %10.1f  %10.1f  %8.4f  %10.4f  %10.4f  %9d%s%s\n",
                r.n_peaks, bval, r.delta_bic, r.aicc, r.aic,
                r.r_squared, r.chi2_red, r.rss, r.n_params, tag, comp)
    end

    if any(r -> r.competitive && r !== best_res, all_results)
        @printf("  * = competitive model (d%s <= %.0f)\n",
                bic_label, cfg.bic_competition_threshold)
    end

    println()
    best_bval = cfg.use_student_bic ? best_res.student_bic : best_res.bic
    println(repeat("=", 80))
    @printf("BEST MODEL: %d Gaussians  |  %s = %.1f  |  R^2 = %.4f  |  chi2_red = %.4f\n",
            n_peaks, bic_label, best_bval, best_res.r_squared, best_res.chi2_red)
    println(repeat("=", 80))

    # Baseline
    @printf("\n  Baseline: y0 = %.4f +/- %.4f\n", popt[1], perr[1])

    # Parameters header
    println()
    @printf("  %5s  %14s  %14s  %14s  %14s\n",
            "Peak", "Center", "Amplitude", "FWHM", "sigma")
    @printf("  %5s  %14s  %14s  %14s  %14s\n",
            "", "($x_unit)", "", "($x_unit)", "($x_unit)")
    @printf("  %5s  %14s  %14s  %14s  %14s\n",
            repeat("-", 5), repeat("-", 10), repeat("-", 10), repeat("-", 10), repeat("-", 10))

    asymmetric = cfg.asymmetric_edges && n_peaks >= 2
    centers = _params_to_centers(popt, n_peaks)
    cerr = _center_errors(perr, n_peaks)
    for i in 0:(n_peaks - 1)
        A = _get_amplitude(popt, i; use_log=cfg.use_log_amplitude)
        sigma = _get_sigma(popt, i)
        A_e = _get_amplitude(perr, i; use_log=cfg.use_log_amplitude)
        s_e = _get_sigma(perr, i)
        fwhm = FWHM_TO_SIGMA * sigma
        fwhm_e = FWHM_TO_SIGMA * s_e
        ce = get(cerr, i + 1, NaN)
        edge_tag = ""
        if asymmetric
            if i == 0
                so = popt[end-1]; so_e = perr[end-1]
                edge_tag = @sprintf("  [outer sigma=%.4f +/-%.4f]", so, so_e)
            elseif i == n_peaks - 1
                so = popt[end]; so_e = perr[end]
                edge_tag = @sprintf("  [outer sigma=%.4f +/-%.4f]", so, so_e)
            end
        end
        @printf("  %5d  %8.3f +/-%6.3f  %8.3f +/-%.3f  %8.3f +/-%.3f  %8.4f +/-%.4f%s\n",
                i + 1, centers[i+1], ce, A, A_e, fwhm, fwhm_e, sigma, s_e, edge_tag)
    end

    if n_peaks > 1
        spacings = [popt[3 + 3 * i] for i in 1:(n_peaks - 1)]
        parts = [@sprintf("%d->%d: %.3f %s", i, i + 1, s, x_unit)
                 for (i, s) in enumerate(spacings)]
        println("\n  SPACING BETWEEN ADJACENT PEAKS:")
        println("  " * join(parts, "  |  "))
        @printf("\n  Max adjacent κ = %.1f  (threshold = %.1f)%s\n",
                best_res.kappa_max_adj, cfg.kappa_max,
                best_res.kappa_max_adj > cfg.kappa_max ? "  ← EXCEEDS THRESHOLD" : "")
    end

    if !isempty(best_res.warnings)
        println("\n  WARNINGS:")
        for w in best_res.warnings
            println("    - $w")
        end
    end

    println(repeat("=", 80))
end

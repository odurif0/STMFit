"""
Opt-in counting diagnostic. Nothing in this module is loaded by production.

The linear subproblem uses the *finite mapped* amplitude/background boxes of
`GaussianFit2D._fit_chain_n`, not an ordinary nonnegative least-squares model.
The outer search refines raw RSS from a native result. It is not an equal-budget
replacement of the native global (geometry-penalized) + local optimization.
"""
module CountingVariableProjection

using LinearAlgebra
using GaussianFit2D
const G = GaussianFit2D
const NL = G.OptimizationNLopt.NLopt

export DiagnosticOptions, read_options, box_lsq, native_raw_bounds,
       linear_layout, design_matrix, linear_coefficients, encode_linear,
       mapping_check, fixed_geometry_profile, refine_profile, circular_to_elliptical,
       finalize_native!, full_gcv

struct DiagnosticOptions
    outer_maxeval::Int
    outer_maxtime_s::Float64
    outer_xtol_rel::Float64
    outer_ftol_rel::Float64
    linear_maxiter::Int
    linear_kkt_atol::Float64
    linear_kkt_rtol::Float64
    linear_bound_atol::Float64
    linear_svd_rtol::Float64
    mapping_atol::Float64
    mapping_rtol::Float64
    native_elliptical_maxiter::Int
end

"Read every numerical budget/tolerance explicitly; never supply defaults."
function read_options(table::AbstractDict)
    names = string.(fieldnames(DiagnosticOptions))
    Set(keys(table)) == Set(names) || error("Diagnostic keys must be exactly: " * join(names, ", "))
    vals = Any[]
    for (key, typ) in zip(names, fieldtypes(DiagnosticOptions))
        v = table[key]
        v isa Real && !(v isa Bool) && isfinite(v) && v > 0 ||
            error("counting_variable_projection.$key must be finite and positive")
        typ == Int && !(v isa Integer) && error("$key must be an integer")
        push!(vals, typ(v))
    end
    opts = DiagnosticOptions(vals...)
    opts.linear_svd_rtol < 1 || error("linear_svd_rtol must be < 1")
    return opts
end

function _svd_step(A, r, rtol)
    isempty(A) && return zeros(size(A, 2)), 0
    f = svd(A; full=false)
    isempty(f.S) && return zeros(size(A, 2)), 0
    keep = findall(>(rtol * maximum(f.S)), f.S)
    isempty(keep) && return zeros(size(A, 2)), 0
    return f.V[:, keep] * ((f.U[:, keep]' * r) ./ f.S[keep]), length(keep)
end

function _box_report(A, b, x, lo, hi, state, iterations, status, opts, rank_free)
    residual = A*x - b
    gradient = A' * residual
    violation = 0.0
    for i in eachindex(x)
        v = state[i] == 2 ? 0.0 : state[i] == -1 ? max(-gradient[i], 0.0) :
            state[i] == 1 ? max(gradient[i], 0.0) : abs(gradient[i])
        violation = max(violation, v)
    end
    scale = max(norm(A' * b, Inf), norm(A' * (A*x), Inf))
    ktol = opts.linear_kkt_atol + opts.linear_kkt_rtol * scale
    feasible = all(lo .<= x .<= hi)
    converged = feasible && isfinite(violation) && violation <= ktol
    return (x=copy(x), rss=sum(abs2, residual), residual_norm=norm(residual),
            gradient=gradient, kkt_violation=violation, kkt_tolerance=ktol,
            converged=converged, status=converged ? "kkt_satisfied" : status,
            iterations=iterations, active_lower=findall(abs.(x-lo) .<= opts.linear_bound_atol),
            active_upper=findall(abs.(hi-x) .<= opts.linear_bound_atol),
            rank_free=rank_free, feasible=feasible)
end

"""
Primal active-set box least squares, minimizing ||A*x-b||².

Free variables are solved by an SVD correction from the current iterate, so a
rank-deficient nullspace does not force an arbitrary reset to zero. A feasible
line step adds the first blocking bound. The largest wrong-sign active gradient
is released only after the free least-squares step. KKT uses A'*(A*x-b) (half
RSS gradient). Fixed variables impose no gradient sign condition. Nonconvergence
is returned, never silently reclassified as a solution.
"""
function box_lsq(A0::AbstractMatrix, b0::AbstractVector, lo0::AbstractVector,
                 hi0::AbstractVector; x0::AbstractVector, options::DiagnosticOptions)
    A, b, lo, hi = Matrix{Float64}(A0), Float64.(b0), Float64.(lo0), Float64.(hi0)
    m, n = size(A)
    m > 0 && n > 0 || throw(ArgumentError("empty least-squares problem"))
    length(b) == m && length(lo) == length(hi) == length(x0) == n || throw(DimensionMismatch())
    all(isfinite, A) && all(isfinite, b) && all(isfinite, lo) && all(isfinite, hi) &&
        all(isfinite, x0) || throw(ArgumentError("nonfinite least-squares input"))
    all(lo .<= hi) || throw(ArgumentError("inverted box"))
    x = clamp.(Float64.(x0), lo, hi)
    state = zeros(Int, n) # -1 lower, 0 free, +1 upper, 2 fixed
    for i in 1:n
        if lo[i] == hi[i]
            state[i] = 2; x[i] = lo[i]
        elseif x[i] - lo[i] <= options.linear_bound_atol
            state[i] = -1; x[i] = lo[i]
        elseif hi[i] - x[i] <= options.linear_bound_atol
            state[i] = 1; x[i] = hi[i]
        end
    end
    rank_free = 0
    for iteration in 1:options.linear_maxiter
        free = findall(==(0), state)
        if !isempty(free)
            d, rank_free = _svd_step(A[:, free], b - A*x, options.linear_svd_rtol)
            alpha, blocker, bound = 1.0, 0, 0
            for (k, i) in enumerate(free)
                if d[k] < 0 && x[i] + d[k] < lo[i]
                    a = (lo[i] - x[i]) / d[k]
                    if a <= alpha; alpha, blocker, bound = max(0.0, a), i, -1; end
                elseif d[k] > 0 && x[i] + d[k] > hi[i]
                    a = (hi[i] - x[i]) / d[k]
                    if a <= alpha; alpha, blocker, bound = max(0.0, a), i, 1; end
                end
            end
            x[free] .+= alpha .* d
            x .= clamp.(x, lo, hi) # only floating-point step roundoff
            if blocker != 0
                state[blocker] = bound
                x[blocker] = bound == -1 ? lo[blocker] : hi[blocker]
                continue
            end
        else
            rank_free = 0
        end
        report = _box_report(A, b, x, lo, hi, state, iteration, "stationarity_failure", options, rank_free)
        report.converged && return report
        worst, release = report.kkt_tolerance, 0
        for i in 1:n
            v = state[i] == -1 ? -report.gradient[i] : state[i] == 1 ? report.gradient[i] : 0.0
            if v > worst; worst, release = v, i; end
        end
        release == 0 && return report # e.g. excessively truncated free SVD
        state[release] = 0
    end
    return _box_report(A, b, x, lo, hi, state, options.linear_maxiter,
                       "iteration_limit", options, rank_free)
end

"""
Mirror the current raw optimizer box inside `_fit_chain_n` (not exported by
core). This is deliberate version coupling, not a new physical calibration.
Actual forward/decode equivalence and native raw endpoints are tested in the
companion synthetic suite, including compiled `_fit_chain_n` zero-iteration probes.
"""
function native_raw_bounds(n::Int, cfg)
    n > 0 || throw(ArgumentError("positive N required"))
    np = G._chain_nparams(n, cfg)
    lo, hi = fill(-10.0, np), fill(10.0, np)
    lo[1], hi[1] = -5.0, 5.0
    j = 2
    if cfg.chain_tilted_baseline
        lo[j:j+1] .= -1.0; hi[j:j+1] .= 1.0; j += 2
    end
    lo[j:j+n-1] .= -5.0; hi[j:j+n-1] .= 5.0; j += n
    ns = G._chain_spacing_param_count(n, cfg)
    lo[j], hi[j] = -4.0, 4.0; j += 1
    lo[j:j+ns-2] .= -5.0; hi[j:j+ns-2] .= 5.0; j += ns-1
    lo[j:j+n-1] .= -3.0; hi[j:j+n-1] .= 3.0; j += n
    nw = (cfg.chain_circular_sigmas ? 1 : 2) * G._chain_sigma_param_count(n, cfg)
    lo[j:j+nw-1] .= -5.0; hi[j:j+nw-1] .= 5.0; j += nw
    nk = G._chain_skew_param_count(n, cfg)
    lo[j:j+nk-1] .= -5.0; hi[j:j+nk-1] .= 5.0; j += nk
    j == np+1 || error("Native parameter layout changed")
    return lo, hi
end

function linear_layout(n, cfg)
    prefix = cfg.chain_tilted_baseline ? 3 : 1
    return (baseline=1:prefix, amplitudes=prefix+1:prefix+n,
            geometry=prefix+n+1:G._chain_nparams(n, cfg),
            names=vcat(cfg.chain_tilted_baseline ? ["baseline", "tilt_x", "tilt_y"] : ["baseline"],
                       ["amplitude_$i" for i in 1:n]))
end

function linear_coefficients(p, n, axisctx, cfg; amp_min, amp_range)
    b0, feats = G._decode_chain(p, n, axisctx, cfg; amp_min=amp_min, amp_range=amp_range)
    return vcat(cfg.chain_tilted_baseline ? [b0, p[2], p[3]] : [b0],
                [f.amplitude for f in feats])
end

"Actual decoded geometry, with exactly the Gaussian forward basis (no surrogate)."
function design_matrix(x, y, p, n, axisctx, cfg; amp_min, amp_range)
    G._chain_peak_profile(cfg) == :gaussian || error("First-pass counting diagnostic requires peak_profile=gaussian")
    _, feats, ts, us = G._decode_chain(p, n, axisctx, cfg; amp_min=amp_min, amp_range=amp_range)
    layout = linear_layout(n, cfg)
    A = Matrix{Float64}(undef, length(x), length(layout.names))
    A[:, 1] .= 1.0
    if cfg.chain_tilted_baseline; A[:, 2] .= x; A[:, 3] .= y; end
    peak_axes = G._chain_peak_axes(ts, us, axisctx, cfg)
    for (j, f) in enumerate(feats)
        ax, ay = peak_axes[j]
        A[:, length(layout.baseline)+j] .= @. exp(-0.5 *
            ((((x-f.x_nm)*ax + (y-f.y_nm)*ay)/f.sigma_x_nm)^2 +
             (((x-f.x_nm)*(-ay) + (y-f.y_nm)*ax)/f.sigma_y_nm)^2))
    end
    return A
end

function _linear_bounds(p, n, axisctx, cfg; amp_min, amp_range)
    rawlo, rawhi = native_raw_bounds(n, cfg)
    layout = linear_layout(n, cfg)
    pl, ph = copy(p), copy(p)
    linidx = vcat(collect(layout.baseline), collect(layout.amplitudes))
    pl[linidx] .= rawlo[linidx]; ph[linidx] .= rawhi[linidx]
    lo = linear_coefficients(pl, n, axisctx, cfg; amp_min, amp_range)
    hi = linear_coefficients(ph, n, axisctx, cfg; amp_min, amp_range)
    return lo, hi
end

function encode_linear(p, coeff, n, axisctx, cfg; amp_min, amp_range)
    layout = linear_layout(n, cfg)
    out = copy(p)
    out[layout.baseline] .= coeff[layout.baseline]
    mapped = !isnan(amp_min) && !isnan(amp_range) && amp_range > G.EPS
    for i in layout.amplitudes
        frac = mapped ? (coeff[i] - amp_min) / amp_range : NaN
        out[i] = mapped ? log(frac / (1-frac)) : log(coeff[i])
    end
    lo, hi = native_raw_bounds(n, cfg)
    out[layout.amplitudes] .= clamp.(out[layout.amplitudes], lo[layout.amplitudes], hi[layout.amplitudes])
    return out
end

function mapping_check(A, coeff, x, y, p, n, axisctx, cfg; amp_min, amp_range, options)
    expected = A * coeff
    actual = G._chain_model_values(x, y, p, n, axisctx, cfg; amp_min, amp_range)
    err = norm(actual - expected, Inf)
    tol = options.mapping_atol + options.mapping_rtol * max(norm(actual, Inf), norm(expected, Inf))
    return (passed=isfinite(err) && err <= tol, max_abs_error=err, tolerance=tol)
end

function fixed_geometry_profile(p0, n, x, y, z, axisctx, cfg; amp_min, amp_range,
                                options::DiagnosticOptions)
    lo, hi = native_raw_bounds(n, cfg)
    length(p0) == length(lo) || throw(DimensionMismatch("native parameter count"))
    all(isfinite, p0) && all(lo .<= p0 .<= hi) || error("Initial vector outside native finite box")
    G._chain_can_fit_support(n, axisctx, cfg) || error("Infeasible native support")
    A = design_matrix(x, y, p0, n, axisctx, cfg; amp_min, amp_range)
    coeff0 = linear_coefficients(p0, n, axisctx, cfg; amp_min, amp_range)
    check0 = mapping_check(A, coeff0, x, y, p0, n, axisctx, cfg; amp_min, amp_range, options)
    check0.passed || error("Initial native forward mapping failed: $check0")
    lower, upper = _linear_bounds(p0, n, axisctx, cfg; amp_min, amp_range)
    linear = box_lsq(A, z, lower, upper; x0=coeff0, options)
    p = encode_linear(p0, linear.x, n, axisctx, cfg; amp_min, amp_range)
    check = mapping_check(A, linear.x, x, y, p, n, axisctx, cfg; amp_min, amp_range, options)
    pred = G._chain_model_values(x, y, p, n, axisctx, cfg; amp_min, amp_range)
    rss = sum(abs2, pred-z)
    return (params=p, linear=linear, mapping=check, initial_mapping=check0,
            rss=rss, lower=lower, upper=upper,
            success=linear.converged && check.passed && isfinite(rss))
end

"Bounded local raw-RSS profiling from a native result; no selector is called."
function refine_profile(p0, n, x, y, z, axisctx, cfg; amp_min, amp_range,
                        options::DiagnosticOptions, initial=nothing)
    started = time_ns()
    first = initial === nothing ? fixed_geometry_profile(p0, n, x, y, z, axisctx, cfg;
                                                        amp_min, amp_range, options) : initial
    first.success || error("Initial profile failed: $(first.linear.status), mapping=$(first.mapping.passed)")
    geom = linear_layout(n, cfg).geometry
    lower, upper = native_raw_bounds(n, cfg)
    best, evaluations, failed = first, 0, 0
    last_error = ""
    objective = function (q, grad)
        isempty(grad) || error("Derivative-free outer algorithm required")
        evaluations += 1
        p = copy(p0); p[geom] .= q
        candidate = fixed_geometry_profile(p, n, x, y, z, axisctx, cfg; amp_min, amp_range, options)
        if !candidate.success
            failed += 1
            last_error = "linear=$(candidate.linear.status), mapping=$(candidate.mapping.passed)"
            return Inf
        end
        if candidate.rss < best.rss; best = candidate; end
        return candidate.rss
    end
    opt = NL.Opt(:LN_BOBYQA, length(geom))
    opt.lower_bounds = lower[geom]; opt.upper_bounds = upper[geom]
    opt.maxeval = options.outer_maxeval
    opt.maxtime = options.outer_maxtime_s
    opt.xtol_rel = options.outer_xtol_rel
    opt.ftol_rel = options.outer_ftol_rel
    opt.min_objective = objective
    status = "not_started"
    try
        _, _, ret = NL.optimize(opt, p0[geom])
        status = string(ret)
    catch err
        status = "exception"
        last_error = sprint(showerror, err)
    end
    # Evaluation/time limits yield usable iterates, not a convergence claim.
    converged = status in ("SUCCESS", "FTOL_REACHED", "XTOL_REACHED", "STOPVAL_REACHED") && failed == 0
    return (best=best, initial=first, status=status, converged=converged,
            evaluations=evaluations, failed_evaluations=failed, error=last_error,
            elapsed_s=(time_ns()-started)/1e9)
end

"Same circular→elliptical raw-width duplication as extract_lobe_features.jl."
function circular_to_elliptical(p, n, cfg_ell)
    cfg_ell.chain_circular_sigmas && error("Elliptical config required")
    prefix = 1 + (cfg_ell.chain_tilted_baseline ? 2 : 0)
    split = prefix+n+G._chain_spacing_param_count(n, cfg_ell)+n
    nw = G._chain_sigma_param_count(n, cfg_ell)
    sigma = p[split+1:split+nw]
    tail = p[split+nw+1:end]
    return vcat(p[1:split], sigma, sigma, tail)
end

full_gcv(rss, ndata, n, cfg) = ndata > G._chain_nparams(n, cfg) ?
    ndata / (ndata-G._chain_nparams(n, cfg))^2 * rss : Inf

function finalize_native!(r, data, cfg)
    cfg.cv_method == "gcv" || error("Diagnostic finalization requires cv_method=gcv (no hidden refits)")
    r.success || return r
    pred = G._chain_model_values(data.x, data.y, r.params, r.n, data.axisctx, cfg;
                                 amp_min=r.amp_min, amp_range=r.amp_range)
    G._finalize_chain_result!(r, data.z, pred, data.noise, r.n, max(10, length(data.z) ÷ 9),
        data.zfull, data.xs, data.ys, data.zimg, data.x, data.y, data.axisctx, cfg)
    return r
end

end # module

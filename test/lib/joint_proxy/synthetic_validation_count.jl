function _profile(view::ViewData, ax)
    ts = Float64[]; ws = Float64[]
    offset = median(view.z)
    for (iy, y) in enumerate(view.ys), (ix, x) in enumerate(view.xs)
        v = view.z[iy, ix]; isfinite(v) || continue
        w = max(v - offset, 0.0); w == 0 && continue
        push!(ts, (x - ax.cx) * ax.ax + (y - ax.cy) * ax.ay)
        push!(ws, w)
    end
    isempty(ws) && return Float64[], Float64[], 0.0
    order = sortperm(ts); ts, ws = ts[order], ws[order]
    bins = collect(range(first(ts), last(ts); length=128))
    prof = zeros(Float64, length(bins))
    for (t, w) in zip(ts, ws)
        i = clamp(searchsortedlast(bins, t), 1, length(bins))
        prof[i] += w
    end
    prof = [mean(prof[max(1, i - 1):min(length(prof), i + 1)]) for i in eachindex(prof)]
    signal = maximum(prof) - median(prof)
    return bins, prof, signal
end

function _support_count_estimate(ts::AbstractVector{<:Real}, prof::AbstractVector{<:Real}, n_min::Int, n_max::Int;
                                 spacing_nm::Float64 = Main.JointProxySimulator.SimulatorConfig().spacing_nm)
    isempty(ts) && return n_min, 0.0, false
    isempty(prof) && return n_min, 0.0, false
    base = median(Float64.(prof))
    signal = maximum(Float64.(prof)) - base
    isfinite(signal) && signal > 0 || return n_min, 0.0, false
    thresh = base + 0.25 * signal
    idx = findall(p -> isfinite(p) && p >= thresh, prof)
    support_nm = isempty(idx) ? 0.0 : max(0.0, Float64(ts[last(idx)]) - Float64(ts[first(idx)]))
    estimate = clamp(round(Int, support_nm / max(spacing_nm, eps(Float64))) + 1, n_min, n_max)
    return estimate, support_nm, true
end

function _gaussian_columns(ts::Vector{Float64}, centers::Vector{Float64}, sigma_nm::Float64)
    cols = Matrix{Float64}(undef, length(ts), length(centers))
    sigma_nm > 0 && isfinite(sigma_nm) || return fill(0.0, length(ts), length(centers))
    inv2s2 = 0.5 / (sigma_nm * sigma_nm)
    for (j, c) in enumerate(centers)
        @inbounds for i in eachindex(ts)
            d = ts[i] - c
            cols[i, j] = exp(-d * d * inv2s2)
        end
    end
    return cols
end

function _ols_gcv_score(ts::Vector{Float64}, prof::Vector{Float64}, centers::Vector{Float64}, sigma_nm::Float64)
    nobs = length(ts)
    nobs == length(prof) || return Inf, Float64[], 0.0, true
    nobs < 4 && return 1.0, centers, 0.0, true
    cols = _gaussian_columns(ts, centers, sigma_nm)
    X = hcat(ones(nobs), cols)
    β = X \ prof
    amps = Float64.(β[2:end])
    max_amp = maximum(abs.(amps); init=0.0)
    # Reject collapsed lobes while allowing profile-bin and center-grid mismatch
    # to attenuate an edge amplitude in otherwise physical long chains.
    if isempty(amps) || max_amp <= 0.0 || any(a -> !isfinite(a) || a <= 0.0 || a < 0.15 * max_amp, amps)
        return Inf, centers, 0.0, true
    end
    resid = prof .- X * β
    rss = sum(abs2, resid)
    # This profile fit estimates one amplitude per lobe plus a baseline, while
    # spacing, width, and offset are selected by the surrounding grid search.
    p = length(centers) + 4
    denom = (nobs - p)
    denom > 0 || return Inf, centers, 0.0, true
    gcv = (rss / nobs) / (max(1e-12, (denom / nobs)^2))
    return gcv, centers, sum(amps), false
end

function _candidate_geometry_scan(ts::Vector{Float64}, prof::Vector{Float64}, n::Int;
                                  spacing_nm::Float64=Main.JointProxySimulator.SimulatorConfig().spacing_nm,
                                  sigma_par_nm::Float64=Main.JointProxySimulator.SimulatorConfig().sigma_par_nm)
    isempty(ts) && return 1.0, Float64[], true
    isempty(prof) && return 1.0, Float64[], true
    base_spacing = spacing_nm
    base_sigma = sigma_par_nm
    spacing_grid = collect(range(0.925 * base_spacing, 1.075 * base_spacing; length=5))
    sigma_grid = collect(range(0.90 * base_sigma, 1.10 * base_sigma; length=5))
    offsets = collect(ts)
    best_gcv = Inf
    best_centers = Float64[]
    for spacing in spacing_grid, sigma in sigma_grid
        for offset in offsets
            centers = [offset + (k - 1 - (n - 1) / 2) * spacing for k in 1:n]
            centers[1] < first(ts) && continue
            centers[end] > last(ts) && continue
            gcv, good_centers, _, rejected = _ols_gcv_score(ts, prof, centers, sigma)
            rejected && continue
            gcv < best_gcv || continue
            best_gcv = gcv
            best_centers = good_centers
        end
    end
    if !isfinite(best_gcv)
        return 1.0, [mean(ts) + (k - 1 - (n - 1) / 2) * base_spacing for k in 1:n], true
    end
    return best_gcv, best_centers, false
end

function _gaussian_kernel1d(sigma_bins::Float64)
    sigma_bins > 0 && isfinite(sigma_bins) || return [1.0]
    radius = max(1, ceil(Int, 3 * sigma_bins))
    xs = collect(-radius:radius)
    weights = exp.(-0.5 .* (Float64.(xs) ./ sigma_bins).^2)
    s = sum(weights)
    s > 0 || return [1.0]
    return weights ./ s
end

function _smooth_profile(prof::AbstractVector{<:Real}, sigma_bins::Float64)
    vals = Float64.(prof)
    isempty(vals) && return Float64[]
    ker = _gaussian_kernel1d(sigma_bins)
    radius = (length(ker) - 1) ÷ 2
    out = similar(vals)
    for i in eachindex(vals)
        acc = 0.0
        wsum = 0.0
        for (k, w) in enumerate(ker)
            j = i + (k - 1) - radius
            1 <= j <= length(vals) || continue
            v = vals[j]
            isfinite(v) || continue
            acc += w * v
            wsum += w
        end
        out[i] = wsum > 0 ? acc / wsum : vals[i]
    end
    return out
end

function _local_maxima(ts::AbstractVector{<:Real}, prof::AbstractVector{<:Real}; min_prominence::Float64)
    isempty(ts) && return Float64[], Float64[]
    isempty(prof) && return Float64[], Float64[]
    xs = Float64.(ts)
    ys = Float64.(prof)
    base = median(ys)
    peaks_t = Float64[]
    peaks_v = Float64[]
    for i in 2:length(ys)-1
        y = ys[i]
        isfinite(y) || continue
        y >= ys[i - 1] || continue
        y > ys[i + 1] || continue
        y - base >= min_prominence || continue
        push!(peaks_t, xs[i])
        push!(peaks_v, y)
    end
    return peaks_t, peaks_v
end

function _nms_centers(peaks_t::Vector{Float64}, peaks_v::Vector{Float64}, min_sep_nm::Float64)
    isempty(peaks_t) && return Float64[]
    order = sortperm(eachindex(peaks_t); by=i -> peaks_v[i], rev=true)
    chosen = Float64[]
    for i in order
        t = peaks_t[i]
        all(abs(t - c) >= min_sep_nm for c in chosen) || continue
        push!(chosen, t)
    end
    sort!(chosen)
    return chosen
end

function _detect_lobe_centers(ts::AbstractVector{<:Real}, prof::AbstractVector{<:Real};
                              spacing_nm::Float64 = Main.JointProxySimulator.SimulatorConfig().spacing_nm,
                              sigma_par_nm::Float64 = Main.JointProxySimulator.SimulatorConfig().sigma_par_nm)
    isempty(ts) && return Float64[], 0, true
    isempty(prof) && return Float64[], 0, true
    dt = length(ts) > 1 ? abs(Float64(ts[2]) - Float64(ts[1])) : spacing_nm
    dt > 0 && isfinite(dt) || return Float64[], 0, true
    sigma_nm = max(0.5 * spacing_nm, 1.5 * sigma_par_nm)
    smoothed = _smooth_profile(prof, sigma_nm / dt)
    signal = maximum(smoothed) - median(smoothed)
    signal > 0 && isfinite(signal) || return Float64[], 0, true
    peaks_t, peaks_v = _local_maxima(ts, smoothed; min_prominence=0.12 * signal)
    centers = _nms_centers(peaks_t, peaks_v, max(0.65 * spacing_nm, 0.9 * sigma_par_nm))
    return centers, length(centers), isempty(centers)
end


module STMMolecularFit

using DelimitedFiles
using LinearAlgebra
using GaussianFit1D
using GaussianFit2D
using STMFitCore: FWHM_TO_SIGMA, sigma_from_fwhm
using Plots
using Printf
using Statistics

export SXMImage, SXMChannel, PreprocessConfig, SlideConfig, SlideResult, FitSlideConfig, SlideFitResult
export Fit2DResult
export ChainSweepConfig, ChainModelResult
export read_sxm, channel_names, get_channel, preprocess_channel, molecule_roi_mask
export extract_slide, write_slide_outputs, fit_slide, write_fit_outputs, extract_and_fit_slide
export fit_2d_chain, chain_gaussian_sweep
export fit_chain_1d_bootstrapped, compare_1d_2d, compare_2d_1d_by_N

const EPS = 1e-12

struct SXMChannel
    name::String
    unit::String
    direction::String
    data::Matrix{Float64}
end

struct SXMImage
    filepath::String
    header::Dict{String,String}
    width::Int
    height::Int
    range_nm::Tuple{Float64,Float64}
    offset_nm::Tuple{Float64,Float64}
    channels::Vector{SXMChannel}
end

Base.@kwdef mutable struct PreprocessConfig
    channel::String = "Z"
    direction::String = "fwd"
    stride::Int = 1
    flatten::String = "plane+rows"
    smooth_radius_px::Int = 1
    roi_channel::String = "Z"
    roi_threshold_fraction::Float64 = 0.35
    roi_noise_k::Float64 = 3.0  # noise-based threshold multiplier
    roi_dilate_px::Int = 10
end

Base.@kwdef mutable struct SlideConfig
    channel::String = "Z"
    direction::String = "fwd"
    stride::Int = 1
    flatten::String = "plane+rows"
    smooth_radius_px::Int = 1
    contrast::String = "bright" # bright | dark | auto
    axis::String = "auto"       # auto | manual
    axis_angle_deg::Float64 = 0.0
    origin_x_nm::Float64 = NaN
    origin_y_nm::Float64 = NaN
    width_nm::Float64 = 0.30
    slide_mode::Symbol = :ridge_mean   # :ridge_mean (crest+disc avg, recommended) | :mean | :ridge
    n_samples::Int = 900
    baseline_quantile::Float64 = 0.10
    support_noise_k::Float64 = 2.5
    support_smooth_radius::Int = 5
    support_padding_nm::Float64 = 0.20
    min_support_nm::Float64 = 1.0
    output_dir::String = "results/slide"
    no_plot::Bool = false
end

Base.@kwdef struct SlideResult
    x::Vector{Float64}
    y::Vector{Float64}
    full_x::Vector{Float64}
    full_y::Vector{Float64}
    counts::Vector{Int}
    origin::Tuple{Float64,Float64}
    axis::Vector{Float64}
    perp::Vector{Float64}
    baseline::Float64
    threshold::Float64
    support_start_nm::Float64
    support_end_nm::Float64
    support_length_nm::Float64
    unit::String
    noise::Float64
    noise_1d::Float64 = 0.0
    arc_ratio::Float64 = 1.0  # arc_length / axial_length (>1 for curved ridges)
end

Base.@kwdef mutable struct FitSlideConfig
    min_spacing::Float64 = 0.4
    max_spacing::Float64 = 0.675
    fwhm_min::Float64 = 0.45
    fwhm_max::Float64 = 1.2
    max_overlap::Float64 = 0.6
    kappa_max::Float64 = 10.0
    kappa_weight::Float64 = 1.0
    amplitude_min_fraction::Float64 = 0.3
    global_maxtime::Float64 = 8.0
    global_maxiter::Int = 5000
    asymmetric_edges::Bool = false
    peak_profile::Symbol = :gaussian
    output_dir::String = "results/slide_fit"
    no_plot::Bool = false
end

Base.@kwdef struct SlideFitResult
    fit_run
    best_model
    results_file::String
    plot_files::Vector{String}
end

function _parse_header(header_text::String)
    header = Dict{String,String}()
    current = nothing
    buf = IOBuffer()
    for line in split(header_text, '\n')
        if startswith(line, ":") && endswith(strip(line), ":")
            current !== nothing && (header[current] = strip(String(take!(buf))))
            current = strip(line, [':', ' ', '\r', '\n', '\t'])
        elseif current !== nothing
            println(buf, line)
        end
    end
    current !== nothing && (header[current] = strip(String(take!(buf))))
    return header
end

function _parse_pair(value::String, T=Float64)
    vals = split(strip(value))
    length(vals) >= 2 || error("Expected two values, got: $value")
    return parse(T, vals[1]), parse(T, vals[2])
end

function _parse_data_info(value::String)
    infos = NamedTuple[]
    for row in split(value, '\n')
        s = strip(row)
        isempty(s) && continue
        startswith(s, "Channel") && continue
        parts = split(s)
        length(parts) >= 6 || continue
        push!(infos, (name=parts[2], unit=parts[3], direction=parts[4]))
    end
    return infos
end

function _read_be_float32(bytes::Vector{UInt8}, offset::Int, n::Int)
    io = IOBuffer(bytes[offset:(offset + 4n - 1)])
    vals = Vector{Float64}(undef, n)
    for i in 1:n
        u = read(io, UInt32)
        ENDIAN_BOM == 0x04030201 && (u = bswap(u))
        vals[i] = Float64(reinterpret(Float32, u))
    end
    return vals
end

function read_sxm(filepath::String)
    bytes = read(filepath)
    marker = b":SCANIT_END:"
    pos = findfirst(marker, bytes)
    pos === nothing && error("Could not find :SCANIT_END: in $filepath")
    header = _parse_header(String(bytes[1:(first(pos)-1)]))
    nx, ny = _parse_pair(header["SCAN_PIXELS"], Int)
    rx_m, ry_m = _parse_pair(header["SCAN_RANGE"], Float64)
    ox_m, oy_m = haskey(header, "SCAN_OFFSET") ? _parse_pair(header["SCAN_OFFSET"], Float64) : (0.0, 0.0)

    expanded = Tuple{String,String,String}[]
    for info in _parse_data_info(header["DATA_INFO"])
        dirs = lowercase(info.direction) == "both" ? ["fwd", "bwd"] : [lowercase(info.direction)]
        for dir in dirs
            push!(expanded, (info.name, info.unit, dir))
        end
    end

    nvals_per_image = nx * ny
    nvals = nvals_per_image * length(expanded)
    data_offset = length(bytes) - 4nvals + 1
    data_offset > 0 || error("Invalid SXM data size")
    vals = _read_be_float32(bytes, data_offset, nvals)

    channels = SXMChannel[]
    k = 1
    for (name, unit, dir) in expanded
        raw = vals[k:(k+nvals_per_image-1)]
        k += nvals_per_image
        data = permutedims(reshape(raw, nx, ny))
        # Nanonis backward scan is stored in acquisition order; flip x so fwd/bwd share coordinates.
        lowercase(dir) == "bwd" && (data = reverse(data; dims=2))
        push!(channels, SXMChannel(name, unit, dir, data))
    end
    return SXMImage(filepath, header, nx, ny, (rx_m*1e9, ry_m*1e9), (ox_m*1e9, oy_m*1e9), channels)
end

channel_names(img::SXMImage) = unique(c.name for c in img.channels)

function get_channel(img::SXMImage, name::String; direction::String="fwd")
    lname, ldir = lowercase(name), lowercase(direction)
    for c in img.channels
        lowercase(c.name) == lname && lowercase(c.direction) == ldir && return c
    end
    for c in img.channels
        lowercase(c.name) == lname && return c
    end
    error("Channel '$name' direction '$direction' not found. Available: " * string([(c.name,c.direction) for c in img.channels]))
end

function _coordinate_vectors(img::SXMImage; stride::Int=1)
    xs = collect(range(0, img.range_nm[1]; length=img.width))[1:stride:end]
    ys = collect(range(0, img.range_nm[2]; length=img.height))[1:stride:end]
    return xs, ys
end

function _value_scale(unit::String)
    u = lowercase(strip(unit))
    u == "m" && return 1e9, "nm"
    u == "a" && return 1e12, "pA"
    return 1.0, unit
end

function _plane_fit(xs, ys, z)
    xflat = repeat(xs', length(ys), 1)[:]
    yflat = repeat(ys, 1, length(xs))[:]
    zflat = vec(z)
    coeff = hcat(ones(length(xflat)), xflat, yflat) \ zflat
    plane = similar(z)
    for iy in eachindex(ys), ix in eachindex(xs)
        plane[iy, ix] = coeff[1] + coeff[2]*xs[ix] + coeff[3]*ys[iy]
    end
    return plane
end

function _row_median_flatten(z)
    out = copy(z)
    for iy in axes(out, 1)
        out[iy, :] .-= median(out[iy, :])
    end
    return out
end

function _box_smooth(z::Matrix{Float64}, radius::Int)
    radius <= 0 && return copy(z)
    out = similar(z)
    ny, nx = size(z)
    for iy in 1:ny, ix in 1:nx
        out[iy, ix] = mean(@view z[max(1,iy-radius):min(ny,iy+radius), max(1,ix-radius):min(nx,ix+radius)])
    end
    return out
end

function preprocess_channel(img::SXMImage, ch::SXMChannel, cfg::PreprocessConfig=PreprocessConfig())
    scale, scaled_unit = _value_scale(ch.unit)
    stride = max(1, cfg.stride)
    xs, ys = _coordinate_vectors(img; stride=stride)
    z = ch.data[1:stride:end, 1:stride:end] .* scale
    finite_vals = z[isfinite.(z)]
    fill_value = isempty(finite_vals) ? 0.0 : median(finite_vals)
    z[.!isfinite.(z)] .= fill_value
    raw = copy(z)
    occursin("plane", lowercase(cfg.flatten)) && (z .-= _plane_fit(xs, ys, z))
    occursin("rows", lowercase(cfg.flatten)) && (z = _row_median_flatten(z))
    z_smooth = _box_smooth(z, cfg.smooth_radius_px)
    noise = 1.4826 * median(abs.(vec(z_smooth) .- median(vec(z_smooth))))
    noise = max(noise, std(vec(z_smooth))*0.1, EPS)
    return xs, ys, raw, z, z_smooth, scaled_unit, noise
end

function _largest_component(mask::BitMatrix)
    ny, nx = size(mask)
    seen = falses(ny, nx)
    best = Tuple{Int,Int}[]
    for iy in 1:ny, ix in 1:nx
        (!mask[iy, ix] || seen[iy, ix]) && continue
        comp = Tuple{Int,Int}[]
        stack = [(iy, ix)]
        seen[iy, ix] = true
        while !isempty(stack)
            y, x = pop!(stack)
            push!(comp, (y, x))
            for yy in max(1,y-1):min(ny,y+1), xx in max(1,x-1):min(nx,x+1)
                if mask[yy, xx] && !seen[yy, xx]
                    seen[yy, xx] = true
                    push!(stack, (yy, xx))
                end
            end
        end
        length(comp) > length(best) && (best = comp)
    end
    out = falses(ny, nx)
    for idx in best
        out[idx...] = true
    end
    return out
end

function _dilate_mask(mask::BitMatrix, radius::Int)
    radius <= 0 && return copy(mask)
    ny, nx = size(mask)
    out = falses(ny, nx)
    for iy in 1:ny, ix in 1:nx
        if mask[iy, ix]
            out[max(1,iy-radius):min(ny,iy+radius), max(1,ix-radius):min(nx,ix+radius)] .= true
        end
    end
    return out
end

function _otsu_threshold(signal::Matrix{Float64})
    """Otsu's automatic threshold: maximizes inter-class variance."""
    v = vec(signal)
    v = v[isfinite.(v) .& (v .> 0)]
    isempty(v) && return 0.0
    nbins = 128
    lo, hi = minimum(v), maximum(v)
    hi <= lo && return 0.0
    bin_edges = range(lo, hi, length=nbins+1)
    counts = zeros(Int, nbins)
    for x in v
        idx = min(nbins, max(1, Int(floor((x-lo)/(hi-lo)*nbins))+1))
        counts[idx] += 1
    end
    total = sum(counts)
    total == 0 && return 0.0
    sumB, wB, max_var, best = 0.0, 0.0, 0.0, 0
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2
    sumT = sum(bin_centers .* counts)
    for i in 1:nbins
        wB += counts[i]
        wB == 0 && continue
        wF = total - wB
        wF <= 0 && break
        sumB += bin_centers[i] * counts[i]
        mB = sumB / wB
        mF = (sumT - sumB) / wF
        var_between = wB * wF * (mB - mF)^2 / (total^2)
        if var_between > max_var
            max_var = var_between
            best = i
        end
    end
    return best > 0 ? bin_centers[best] : 0.0
end

function molecule_roi_mask(img::SXMImage, cfg::PreprocessConfig=PreprocessConfig())
    ch = get_channel(img, cfg.roi_channel; direction="fwd")
    xs, ys, _raw, _z, z_smooth, _unit, _noise = preprocess_channel(img, ch, cfg)
    finite_vals = z_smooth[isfinite.(z_smooth)]
    isempty(finite_vals) && return xs, ys, trues(size(z_smooth))
    signal = z_smooth .- minimum(finite_vals)
    signal[.!isfinite.(signal)] .= 0.0
    maxsig = maximum(signal)
    maxsig <= EPS && return xs, ys, trues(size(z_smooth))
    noise = 1.4826 * median(abs.(vec(z_smooth) .- median(vec(z_smooth))))
    noise = max(noise, std(vec(z_smooth))*0.1, EPS)
    threshold_otsu = _otsu_threshold(signal)
    threshold_frac = cfg.roi_threshold_fraction * maxsig
    threshold_noise = cfg.roi_noise_k * noise
    threshold = max(threshold_otsu, threshold_frac, threshold_noise)
    mask = signal .>= threshold
    mask = _largest_component(mask)
    mask = _dilate_mask(mask, max(0, cfg.roi_dilate_px ÷ max(1, cfg.stride)))
    return xs, ys, mask
end

function _weighted_axis(x, y, z)
    keep = isfinite.(x) .& isfinite.(y) .& isfinite.(z)
    x, y, z = x[keep], y[keep], z[keep]
    isempty(z) && error("Cannot estimate weighted axis: no finite points")
    w = z .- minimum(z) .+ EPS  # all pixels contribute, not just the brightest
    sw = sum(w)
    ox = sum(x .* w) / sw
    oy = sum(y .* w) / sw
    X = hcat(x .- ox, y .- oy)
    W = Diagonal(w ./ maximum(w))
    _, _, V = svd(W * X; full=false)
    axis = collect(V[:, 1])
    axis ./= max(norm(axis), EPS)
    (axis[2] < 0 || (abs(axis[2]) < EPS && axis[1] < 0)) && (axis .*= -1)
    return (origin=(ox, oy), axis=axis)
end

function _bilinear(xs, ys, z, x0, y0)
    ix = clamp(searchsortedlast(xs, x0), 1, length(xs)-1)
    iy = clamp(searchsortedlast(ys, y0), 1, length(ys)-1)
    x1, x2 = xs[ix], xs[ix+1]
    y1, y2 = ys[iy], ys[iy+1]
    tx = (x0-x1)/max(x2-x1, EPS)
    ty = (y0-y1)/max(y2-y1, EPS)
    return (1-tx)*(1-ty)*z[iy,ix] + tx*(1-ty)*z[iy,ix+1] + (1-tx)*ty*z[iy+1,ix] + tx*ty*z[iy+1,ix+1]
end

function _line_bounds(xs, ys, origin, axis)
    ox, oy = origin
    ax, ay = axis
    cand = Float64[]
    for xedge in (minimum(xs), maximum(xs)); abs(ax) > EPS && push!(cand, (xedge-ox)/ax); end
    for yedge in (minimum(ys), maximum(ys)); abs(ay) > EPS && push!(cand, (yedge-oy)/ay); end
    valid = [t for t in cand if minimum(xs)-1e-9 <= ox+t*ax <= maximum(xs)+1e-9 && minimum(ys)-1e-9 <= oy+t*ay <= maximum(ys)+1e-9]
    length(valid) >= 2 || error("Could not intersect slide axis with image bounds")
    return minimum(valid), maximum(valid)
end

function _trace_crest(xs, ys, z, ts, origin, axis, perp, halfw::Float64; step_u::Float64=0.02)
    """Robust crest tracing with continuity constraint.
    Returns (rx, ry, arc, vals_max) where rx,ry are crest positions, arc is arc length,
    vals_max are the maximum pixel values at each crest position."""
    n = length(ts)
    rx = Float64[]; ry = Float64[]
    vals_max = zeros(n)
    # Pass 1: find max perpendicular to axis at each t
    for i in 1:n
        t = ts[i]
        best_val = -Inf; bx = origin[1] + t*axis[1]; by = origin[2] + t*axis[2]
        for u in range(-halfw, halfw, step=step_u)
            x0 = origin[1] + t*axis[1] + u*perp[1]
            y0 = origin[2] + t*axis[2] + u*perp[2]
            if minimum(xs) <= x0 <= maximum(xs) && minimum(ys) <= y0 <= maximum(ys)
                v = _bilinear(xs, ys, z, x0, y0)
                if v > best_val; best_val = v; bx = x0; by = y0; end
            end
        end
        push!(rx, bx); push!(ry, by)
        vals_max[i] = best_val > -Inf ? best_val : _bilinear(xs, ys, z, origin[1]+t*axis[1], origin[2]+t*axis[2])
    end
    # Smooth and compute arc from crest positions
    rad = max(1, n ÷ 30)
    sx = [mean(rx[max(1,i-rad):min(n,i+rad)]) for i in 1:n]
    sy = [mean(ry[max(1,i-rad):min(n,i+rad)]) for i in 1:n]
    # Arc from smoothed crest
    arc = [0.0]
    for i in 2:n
        d = sqrt((sx[i]-sx[i-1])^2 + (sy[i]-sy[i-1])^2)
        push!(arc, arc[end] + d)
    end
    return (rx=rx, ry=ry, arc=arc, vals_max=vals_max, sx=sx, sy=sy)
end

function _smooth1d(v, radius::Int)
    radius <= 0 && return copy(v)
    out = similar(v)
    for i in eachindex(v)
        out[i] = mean(@view v[max(firstindex(v), i-radius):min(lastindex(v), i+radius)])
    end
    return out
end

_mad_std(v) = 1.4826 * median(abs.(v .- median(v)))

function _otsu_threshold_1d(signal::AbstractVector{Float64})
    v = signal[isfinite.(signal) .& (signal .> 0)]
    isempty(v) && return 0.0
    nbins = 128
    lo, hi = minimum(v), maximum(v)
    hi <= lo && return 0.0
    bin_edges = range(lo, hi, length=nbins+1)
    counts = zeros(Int, nbins)
    for x in v
        idx = min(nbins, max(1, Int(floor((x-lo)/(hi-lo)*nbins))+1))
        counts[idx] += 1
    end
    total = sum(counts)
    total == 0 && return 0.0
    sumB, wB, max_var, best = 0.0, 0.0, 0.0, 0
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2
    sumT = sum(bin_centers .* counts)
    for i in 1:nbins
        wB += counts[i]
        wB == 0 && continue
        wF = total - wB
        wF <= 0 && break
        sumB += bin_centers[i] * counts[i]
        mB = sumB / wB
        mF = (sumT - sumB) / wF
        var_between = wB * wF * (mB - mF)^2 / (total^2)
        if var_between > max_var
            max_var = var_between
            best = i
        end
    end
    return best > 0 ? bin_centers[best] : 0.0
end

function _support_from_baseline(dist, y, cfg::SlideConfig)
    ys = _smooth1d(y, cfg.support_smooth_radius)
    baseline = quantile(ys, cfg.baseline_quantile)
    peak = maximum(ys)
    low = ys[ys .<= baseline]
    noise = isempty(low) ? _mad_std(ys) : _mad_std(low)
    signal_above_baseline = ys .- baseline
    threshold_otsu = baseline + _otsu_threshold_1d(signal_above_baseline[signal_above_baseline .> 0])
    threshold_noise = baseline + cfg.support_noise_k * max(noise, EPS)
    # Use the LESS aggressive threshold (min) to avoid over-truncating support.
    # Otsu can be too aggressive on 1D profiles; noise alone can be too permissive.
    threshold = min(threshold_otsu, threshold_noise)
    active = findall(>(threshold), ys)
    isempty(active) && error("No active support found")
    runs = UnitRange{Int}[]
    s = active[1]; prev = active[1]
    for idx in active[2:end]
        if idx == prev + 1
            prev = idx
        else
            push!(runs, s:prev); s = prev = idx
        end
    end
    push!(runs, s:prev)
    # Filter by min length
    long_runs = filter(r -> dist[last(r)] - dist[first(r)] >= cfg.min_support_nm, runs)
    isempty(long_runs) && error("No active component longer than min_support_nm")
    # Select component containing global peak, else longest
    peak_idx = argmax(ys)
    containing = filter(r -> first(r) <= peak_idx <= last(r), long_runs)
    best = isempty(containing) ? argmax(r -> dist[last(r)] - dist[first(r)], long_runs) :
           argmax(r -> dist[last(r)] - dist[first(r)], containing)
    lo = max(first(dist), dist[first(best)] - cfg.support_padding_nm)
    hi = min(last(dist), dist[last(best)] + cfg.support_padding_nm)
    keep = findall(i -> lo <= dist[i] <= hi, eachindex(dist))
    return keep, baseline, threshold, lo, hi
end

function extract_slide(img::SXMImage, cfg::SlideConfig=SlideConfig())
    pcfg = PreprocessConfig(channel=cfg.channel, direction=cfg.direction, stride=cfg.stride, flatten=cfg.flatten,
                            smooth_radius_px=cfg.smooth_radius_px, roi_channel=cfg.channel)
    ch = get_channel(img, cfg.channel; direction=cfg.direction)
    xs, ys, _raw, z, _zs, unit, noise = preprocess_channel(img, ch, pcfg)
    contrast = lowercase(cfg.contrast)
    contrast == "dark" && (z = .-z)
    if contrast == "auto"
        hi, lo = quantile(vec(z), 0.995), quantile(vec(z), 0.005)
        abs(lo) > abs(hi) && (z = .-z)
    end
    # Fuse Z bwd if available (matches 2D model, improves axis centering)
    has_bwd = any(c -> lowercase(c.name) == lowercase(cfg.channel) && lowercase(c.direction) == "bwd", img.channels)
    if has_bwd
        ch_bwd = get_channel(img, cfg.channel; direction="bwd")
        _, _, _, z_bwd, _, _, _ = preprocess_channel(img, ch_bwd, pcfg)
        z = (z .+ z_bwd) ./ 2.0
    end
    if lowercase(cfg.axis) == "manual"
        if isnan(cfg.origin_x_nm) || isnan(cfg.origin_y_nm)
            _rx, _ry, mask = molecule_roi_mask(img, pcfg)
            mx = Float64[]; my = Float64[]
            for iy in eachindex(ys), ix in eachindex(xs)
                mask[iy, ix] || continue
                push!(mx, xs[ix]); push!(my, ys[iy])
            end
            origin = (sum(mx)/length(mx), sum(my)/length(my))
        else
            origin = (cfg.origin_x_nm, cfg.origin_y_nm)
        end
        θ = deg2rad(cfg.axis_angle_deg)
        axis = [cos(θ), sin(θ)]
    else
        # Use 2D model's ROI logic (fused Z, same mask as chain sweep) for axis
        img2d = GaussianFit2D.read_sxm(img.filepath)
        _, _, _, _, xr2, yr2, zr2, _ = GaussianFit2D._fused_roi_data(img2d,
            GaussianFit2D.PatternConfig(filepath=img.filepath, channel=cfg.channel, direction="fwd",
                stride=cfg.stride, flatten=cfg.flatten, smooth_radius_px=cfg.smooth_radius_px))
        axis2d = GaussianFit2D._weighted_roi_axis(xr2, yr2, zr2)
        origin = axis2d.origin
        axis_val = axis2d.axis
        # Still collect mask pixels for slide profile extraction
        _rx, _ry, mask = molecule_roi_mask(img, pcfg)
        xv = Float64[]; yv = Float64[]; zv = Float64[]
        for iy in eachindex(ys), ix in eachindex(xs)
            mask[iy, ix] || continue
            push!(xv, xs[ix]); push!(yv, ys[iy]); push!(zv, z[iy, ix])
        end
        if isempty(xv)
            @warn "Automatic ROI was empty; falling back to all image pixels for slide axis estimation"
            for iy in eachindex(ys), ix in eachindex(xs)
                push!(xv, xs[ix]); push!(yv, ys[iy]); push!(zv, z[iy, ix])
            end
        end
        axis = axis_val
    end
    axis ./= max(norm(axis), EPS)
    (axis[2] < 0 || (abs(axis[2]) < EPS && axis[1] < 0)) && (axis .*= -1)
    perp = [-axis[2], axis[1]]
    tmin, tmax = _line_bounds(xs, ys, origin, axis)
    ts = collect(range(tmin, tmax; length=cfg.n_samples))
    dist = ts .- first(ts)
    vals = zeros(length(ts)); counts = zeros(Int, length(ts))
    if cfg.width_nm <= 0
        for i in eachindex(ts)
            vals[i] = _bilinear(xs, ys, z, origin[1]+ts[i]*axis[1], origin[2]+ts[i]*axis[2])
            counts[i] = 1
    end
    else
        halfw = cfg.width_nm/2
        if cfg.slide_mode == :ridge
            halfw = cfg.width_nm / 2
            crest = _trace_crest(xs, ys, z, ts, origin, axis, perp, halfw; step_u=0.02)
            for i in eachindex(ts)
                vals[i] = crest.vals_max[i]
                counts[i] = 1
            end
            dist = crest.arc
        elseif cfg.slide_mode == :ridge_mean
            # Trace crest (shared with :ridge), then average in disc around each point
            halfw = cfg.width_nm / 2
            crest = _trace_crest(xs, ys, z, ts, origin, axis, perp, halfw; step_u=0.02)
            r_avg = 0.15
            for i in eachindex(ts)
                bx, by = crest.rx[i], crest.ry[i]
                acc, n = 0.0, 0
                for dx in range(-r_avg, r_avg, step=0.02), dy in range(-r_avg, r_avg, step=0.02)
                    dx*dx + dy*dy <= r_avg*r_avg || continue
                    x0 = bx + dx; y0 = by + dy
                    if minimum(xs) <= x0 <= maximum(xs) && minimum(ys) <= y0 <= maximum(ys)
                        acc += _bilinear(xs, ys, z, x0, y0); n += 1
                    end
                end
                vals[i] = n > 0 ? acc / n : _bilinear(xs, ys, z, origin[1]+ts[i]*axis[1], origin[2]+ts[i]*axis[2])
                counts[i] = 1
            end
            dist = crest.arc
        else  # :mean (default)
            for iy in eachindex(ys), ix in eachindex(xs)
                dx = xs[ix]-origin[1]; dy = ys[iy]-origin[2]
                t = dx*axis[1] + dy*axis[2]
                u = dx*perp[1] + dy*perp[2]
                abs(u) <= halfw || continue
                tmin <= t <= tmax || continue
                b = clamp(Int(floor((t-tmin)/max(tmax-tmin, EPS)*(length(ts)-1)))+1, 1, length(ts))
                vals[b] += z[iy, ix]; counts[b] += 1
            end
            for i in eachindex(vals)
                vals[i] = counts[i] > 0 ? vals[i]/counts[i] : _bilinear(xs, ys, z, origin[1]+ts[i]*axis[1], origin[2]+ts[i]*axis[2])
            end
        end
    end
    vals0 = vals .- quantile(vals, cfg.baseline_quantile)
    vals0 .-= minimum(vals0)
    keep, baseline, threshold, lo, hi = _support_from_baseline(dist, vals0, cfg)
    xcrop = dist[keep] .- first(dist[keep])
    ycrop = vals0[keep] .- minimum(vals0[keep])
    # 1D-specific noise estimate from adjacent differences (independent of fit)
    dy = diff(ycrop)
    noise_1d = isempty(dy) ? noise : 1.4826 * median(abs.(dy .- median(dy))) / sqrt(2)
    return SlideResult(x=xcrop, y=ycrop, full_x=dist, full_y=vals0, counts=counts, origin=(origin[1],origin[2]),
                       axis=axis, perp=perp, baseline=baseline, threshold=threshold, support_start_nm=lo,
                       support_end_nm=hi, support_length_nm=hi-lo, unit=unit, noise=noise, noise_1d=noise_1d)
end

function write_slide_outputs(slide::SlideResult, cfg::SlideConfig=SlideConfig())
    mkpath(cfg.output_dir)
    profile = joinpath(cfg.output_dir, "slide_profile.txt")
    full = joinpath(cfg.output_dir, "slide_full_profile.txt")
    meta = joinpath(cfg.output_dir, "slide_metadata.tsv")
    writedlm(profile, hcat(slide.x, slide.y))
    writedlm(full, hcat(slide.full_x, slide.full_y, slide.counts))
    open(meta, "w") do io
        println(io, "key\tvalue")
        for (k,v) in [("channel",cfg.channel),("direction",cfg.direction),("axis_mode",cfg.axis),("width_nm",cfg.width_nm),
                      ("origin_x_nm",slide.origin[1]),("origin_y_nm",slide.origin[2]),("axis_x",slide.axis[1]),("axis_y",slide.axis[2]),
                      ("baseline",slide.baseline),("threshold",slide.threshold),("support_start_nm",slide.support_start_nm),
                      ("support_end_nm",slide.support_end_nm),("support_length_nm",slide.support_length_nm),("profile",profile),("full_profile",full)]
            println(io, "$(k)\t$(v)")
        end
    end
    plotfile = nothing
    if !cfg.no_plot
        p = plot(slide.full_x, slide.full_y; label="full slide", xlabel="distance (nm)", ylabel="intensity", title="Extracted STM slide")
        hline!(p, [slide.threshold]; label="support threshold", linestyle=:dash)
        vline!(p, [slide.support_start_nm, slide.support_end_nm]; color=:red, label="support")
        plot!(p, slide.x .+ slide.support_start_nm, slide.y; label="exported profile", linewidth=3)
        plotfile = joinpath(cfg.output_dir, "slide_profile.png")
        savefig(p, plotfile)
    end
    return (profile=profile, full=full, metadata=meta, plot=plotfile)
end

function fit_slide(slide::SlideResult, cfg::FitSlideConfig=FitSlideConfig())
    mkpath(cfg.output_dir)
    profile_file = joinpath(cfg.output_dir, "slide_profile_for_fit.txt")
    writedlm(profile_file, hcat(slide.x, slide.y))
    return fit_slide(profile_file, cfg; noise_estimate=slide.noise_1d, arc_ratio=slide.arc_ratio)
end

function fit_slide(profile_file::String, cfg::FitSlideConfig=FitSlideConfig(); noise_estimate::Float64=NaN, arc_ratio::Float64=1.0)
    fit_dir = joinpath(cfg.output_dir, "fit_1d")
    mkpath(fit_dir)
    # Scale max_spacing when using arc length (avoids overfitting on stretched x-axis)
    max_sp = cfg.max_spacing
    mgf_cfg = GaussianFit1D.build_config(Dict{String,Any}(
        "filepath" => profile_file,
        "output_dir" => fit_dir,
        "min_spacing" => cfg.min_spacing,
        "max_spacing" => max_sp,
        "fwhm_min" => cfg.fwhm_min,
        "fwhm_max" => cfg.fwhm_max,
        "max_overlap" => cfg.max_overlap,
        "kappa_max" => cfg.kappa_max,
        "kappa_weight" => cfg.kappa_weight,
        "amplitude_min_fraction" => cfg.amplitude_min_fraction,
        "global_maxtime" => cfg.global_maxtime,
        "global_maxiter" => cfg.global_maxiter,
        "asymmetric_edges" => cfg.asymmetric_edges,
        "peak_profile" => cfg.peak_profile,
        "no_show" => true,
    ))
    if isfinite(noise_estimate)
        mgf_cfg.noise_estimate = noise_estimate
    end
    fr = GaussianFit1D.run_fit(mgf_cfg; save_cache=true, verbose=true)
    isempty(fr.all_results) && error("No multi-Gaussian models were fitted")
    GaussianFit1D.update_model_rankings(fr.all_results, mgf_cfg)
    best = GaussianFit1D.best_result(fr)
    results_file = GaussianFit1D.export_results(fr.x, fr.y, fr.all_results, mgf_cfg)
    plot_files = cfg.no_plot ? String[] : GaussianFit1D.plot_results(fr.x, fr.y, best, fr.all_results, mgf_cfg)
    return SlideFitResult(fit_run=fr, best_model=best, results_file=results_file, plot_files=plot_files)
end

function write_fit_outputs(fit::SlideFitResult, cfg::FitSlideConfig=FitSlideConfig())
    mkpath(cfg.output_dir)
    model_file = joinpath(cfg.output_dir, "model_selection.tsv")
    open(model_file, "w") do io
        println(io, "n_peaks\tn_params\tbic\tdelta_bic\taicc\tr_squared\tchi2_red\trss\tcompetitive\tselected")
        for r in sort(fit.fit_run.all_results; by=r -> r.n_peaks)
            println(io, "$(r.n_peaks)\t$(r.n_params)\t$(r.bic)\t$(r.delta_bic)\t$(r.aicc)\t$(r.r_squared)\t$(r.chi2_red)\t$(r.rss)\t$(r.competitive)\t$(r === fit.best_model)")
        end
    end
    best_file = joinpath(cfg.output_dir, "best_model.tsv")
    best_plot = ""
    idx = findfirst(p -> occursin("best", lowercase(p)), fit.plot_files)
    idx !== nothing && (best_plot = fit.plot_files[idx])
    open(best_file, "w") do io
        b = fit.best_model
        println(io, "key\tvalue")
        println(io, "best_n_peaks\t$(b.n_peaks)")
        println(io, "bic\t$(b.bic)")
        println(io, "aicc\t$(b.aicc)")
        println(io, "r_squared\t$(b.r_squared)")
        println(io, "chi2_red\t$(b.chi2_red)")
        println(io, "rss\t$(b.rss)")
        println(io, "n_params\t$(b.n_params)")
        println(io, "spacing_min_configured\t$(fit.fit_run.cfg.min_spacing)")
        println(io, "spacing_min_effective\t$(GaussianFit1D._effective_min_spacing(fit.fit_run.cfg))")
        println(io, "max_overlap\t$(fit.fit_run.cfg.max_overlap)")
        println(io, "results_file\t$(fit.results_file)")
        println(io, "best_plot\t$(best_plot)")
    end
    peaks_file = joinpath(cfg.output_dir, "best_peaks.tsv")
    centers = GaussianFit1D._params_to_centers(fit.best_model.popt, fit.best_model.n_peaks)
    open(peaks_file, "w") do io
        println(io, "peak\tcenter_nm\tamplitude\tsigma_nm\tfwhm_nm")
        for i in 0:(fit.best_model.n_peaks-1)
            A = GaussianFit1D._get_amplitude(fit.best_model.popt, i)
            σ = GaussianFit1D._get_sigma(fit.best_model.popt, i)
            println(io, "$(i+1)\t$(centers[i+1])\t$A\t$σ\t$(GaussianFit1D.FWHM_TO_SIGMA * σ)")
        end
    end
    return (model_selection=model_file, best_model=best_file, best_peaks=peaks_file)
end

include("Fit2D.jl")
using .Fit2D

include("selectors.jl")

function extract_and_fit_slide(filepath::String; slide_config::SlideConfig=SlideConfig(), fit_config::FitSlideConfig=FitSlideConfig(), write_outputs::Bool=true)
    img = read_sxm(filepath)
    slide = extract_slide(img, slide_config)
    slide_files = write_outputs ? write_slide_outputs(slide, slide_config) : nothing
    profile_for_fit = slide_files === nothing ? slide : slide_files.profile
    fit = fit_slide(profile_for_fit, fit_config)
    fit_files = write_outputs ? write_fit_outputs(fit, fit_config) : nothing
    return (slide=slide, slide_files=slide_files, fit=fit, fit_files=fit_files)
end

# ===========================================================================
# 1D-Bootstrapped 2D Chain Fit
# ===========================================================================

function fit_chain_1d_bootstrapped(filepath::String;
        chain_config::GaussianFit2D.ChainSweepConfig=GaussianFit2D.ChainSweepConfig(),
        slide_config::SlideConfig=SlideConfig(),
        fit_config::FitSlideConfig=FitSlideConfig(),
        boot_halfwidth::Int=2,
        no_plot::Bool=false)
    """Run 1D slide extraction + fit, then a narrow 2D chain sweep around the 1D BIC-best N.
    
    Compares N values N₁ - halfwidth … N₁ + halfwidth with 1D bootstrap init.
    This is comparable to a full sweep but ~3-5× faster."""
    # Step 1: Extract and fit 1D slide (using STMMolecularFit's own pipeline)
    img_own = read_sxm(filepath)
    slide = extract_slide(img_own, slide_config)
    slide_files = write_slide_outputs(slide, slide_config)
    fit = fit_slide(slide_files.profile, fit_config)
    fit_files = write_fit_outputs(fit, fit_config)

    # Step 2: Extract 1D bootstrap parameters
    best = fit.best_model
    N_1d = best.n_peaks
    centers_1d = GaussianFit1D._params_to_centers(best.popt, N_1d)
    amps_1d = [GaussianFit1D._get_amplitude(best.popt, i) for i in 0:(N_1d-1)]
    sigmas_1d = [GaussianFit1D._get_sigma(best.popt, i) for i in 0:(N_1d-1)]
    sigma_parallel_1d = mean(sigmas_1d)

    # Step 3: Build ChainSweepConfig with narrow sweep around N_1d, 1D bootstrap init
    ccfg = deepcopy(chain_config)
    ccfg.n_min = max(2, N_1d - boot_halfwidth)
    ccfg.n_max = N_1d + boot_halfwidth
    ccfg.init_centers_t = centers_1d
    ccfg.init_amplitudes = amps_1d
    ccfg.init_sigma_parallel = isfinite(ccfg.init_sigma_parallel) ? ccfg.init_sigma_parallel : sigma_parallel_1d
    ccfg.boot_sweep_halfwidth = boot_halfwidth

    # Step 4: Run narrow 2D chain sweep
    img_2d = GaussianFit2D.read_sxm(filepath)
    pcfg = GaussianFit2D.PatternConfig(
        filepath=filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1,
        output_dir=slide_config.output_dir)
    results, best_chain, ctx = GaussianFit2D.chain_gaussian_sweep(img_2d, pcfg, ccfg)

    # Write outputs
    mkpath(slide_config.output_dir)
    ccfg_write = deepcopy(ccfg)
    GaussianFit2D._write_chain_sweep(results, best_chain, ctx, pcfg, ccfg_write)
    no_plot || GaussianFit2D._plot_chain_sweep(results, best_chain, ctx, pcfg, ccfg_write)

    println("1D_best_N: $N_1d (BIC=$(best.bic))")
    println("2D_boot_sweep: N=$(ccfg.n_min)..$(ccfg.n_max)")
    println("2D_boot_best_N: $(best_chain.n) valid=$(best_chain.valid) BIC=$(best_chain.bic)")
    println("2D_boot_mean_spacing_nm: $(best_chain.mean_spacing_nm)")
    println("2D_boot_residual_peak_snr: $(best_chain.residual_peak_snr)")

    return (slide=slide, slide_files=slide_files, fit=fit, fit_files=fit_files,
            chain_results=results, chain_best=best_chain, chain_config=ccfg)
end

function compare_1d_2d(filepath::String;
        chain_config::GaussianFit2D.ChainSweepConfig=GaussianFit2D.ChainSweepConfig(),
        slide_config::SlideConfig=SlideConfig(),
        fit_config::FitSlideConfig=FitSlideConfig(),
        boot_halfwidth::Int=2,
        no_plot::Bool=false)
    """Compare 1D slice fitting with 2D chain fitting (both sweep and narrow 1D-bootstrapped sweep)."""
    t_boot = @elapsed bootstrap = fit_chain_1d_bootstrapped(filepath;
        chain_config=chain_config, slide_config=slide_config, fit_config=fit_config,
        boot_halfwidth=boot_halfwidth, no_plot=no_plot)

    ccfg_sweep = deepcopy(chain_config)
    ccfg_sweep.n_min = 2
    ccfg_sweep.n_max = 14
    pcfg = GaussianFit2D.PatternConfig(
        filepath=filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1,
        output_dir=slide_config.output_dir)
    img_2d = GaussianFit2D.read_sxm(filepath)
    t_sweep = @elapsed results_sweep, best_sweep, ctx_sweep = GaussianFit2D.chain_gaussian_sweep(img_2d, pcfg, ccfg_sweep)

    println("\n=== COMPARAISON 1D vs 2D ===")
    println("1D N=$(bootstrap.fit.best_model.n_peaks) BIC=$(bootstrap.fit.best_model.bic)")
    println("2D bootstrapped sweep N=$(bootstrap.chain_best.n) valid=$(bootstrap.chain_best.valid) BIC=$(bootstrap.chain_best.bic)")
    println("2D full sweep N=$(best_sweep.n) valid=$(best_sweep.valid) BIC=$(best_sweep.bic)")
    @printf("\nTiming:\n  bootstrapped narrow sweep: %.1f s (N=%d..%d)\n  full sweep:               %.1f s (N=2..14)\n",
        t_boot, bootstrap.chain_config.n_min, bootstrap.chain_config.n_max, t_sweep)
    @printf("  Speedup: %.1f×\n", t_sweep / t_boot)

    return (bootstrapped=bootstrap, sweep=(results=results_sweep, best=best_sweep, ctx=ctx_sweep),
            timing=(t_bootstrapped=t_boot, t_sweep=t_sweep))
end

function compare_2d_1d_by_N(filepath::String;
        output_dir::String = "results/comparison_2d_1d_by_N",
        slide_width_nm::Float64 = 0.30,
        support_noise_k::Float64 = 2.5,
        support_padding_nm::Float64 = 0.20,
        n_min::Int = 2, n_max::Int = 14,
        spacing_min_nm::Float64 = 0.35, spacing_max_nm::Float64 = 0.75,
        fit_width_nm::Float64 = 0.15,
        global_maxtime::Float64 = 20.0,
        skip_plots::Bool = false,
        min_spacing_1d::Float64 = 0.35, max_spacing_1d::Float64 = 0.75,
        max_overlap::Float64 = 0.6)
    """
    compare_2d_1d_by_N(filepath; ...)

    For each N in the 2D chain sweep, build a 2×2 comparison plot:
    - Top-left: 2D heatmap + FWHM ellipses + axis
    - Top-right: 1D slide profile + multi-Gaussian fit components
    - Bottom-left: 2D residuals / noise
    - Bottom-right: 1D residuals

    Prints a summary table of ΔsBIC for both methods.

    Returns `(plot_paths=..., best_1d=..., best_2d=..., results_2d=..., all_results_1d=...)`.
    """
    FWHM_SIGMA = 2.355
    function _ellipse!(p, x0, y0, a, b, angle; color=:cyan, alpha=0.3, label="")
        θ = range(0, 2π, length=72)
        cosθ, sinθ = cos.(θ), sin.(θ)
        ca, sa = cos(angle), sin(angle)
        xe = x0 .+ a .* cosθ .* ca .- b .* sinθ .* sa
        ye = y0 .+ a .* cosθ .* sa .+ b .* sinθ .* ca
        plot!(p, xe, ye; color=color, alpha=alpha, label=label, linewidth=1.5)
    end

    COLORMAP_RESID = cgrad([:blue, :lightgray, :red])
    mkpath(output_dir)

    # ━━━ 1D extraction + fit ━━━
    println("Extracting 1D slide...")
    slide_cfg = SlideConfig(
        width_nm=slide_width_nm,
        support_noise_k=support_noise_k, support_padding_nm=support_padding_nm,
        output_dir=output_dir, no_plot=true)
    img_own = read_sxm(filepath)
    slide = extract_slide(img_own, slide_cfg)
    write_slide_outputs(slide, slide_cfg)

    println("Fitting 1D profile...")
    fit_cfg = FitSlideConfig(
        min_spacing=min_spacing_1d, max_spacing=max_spacing_1d,
        max_overlap=max_overlap,
        output_dir=output_dir)
    fit_1d = fit_slide(slide, fit_cfg)
    write_fit_outputs(fit_1d, fit_cfg)

    # ━━━ 2D chain sweep ━━━
    println("Running 2D chain sweep...")
    img = GaussianFit2D.read_sxm(filepath)
    pcfg = GaussianFit2D.PatternConfig(
        filepath=filepath, channel="Z", direction="fwd",
        stride=1, flatten="plane+rows", smooth_radius_px=1,
        output_dir=output_dir, no_plot=false)
    ccfg = GaussianFit2D.ChainSweepConfig(
        n_min=n_min, n_max=n_max,
        spacing_min_nm=spacing_min_nm, spacing_max_nm=spacing_max_nm,
        fit_width_nm=fit_width_nm,
        support_noise_k=support_noise_k,
        support_padding_nm=support_padding_nm,
        max_overlap=max_overlap,
        sigma_parallel_min_nm=fit_cfg.fwhm_min / GaussianFit1D.FWHM_TO_SIGMA,
        sigma_parallel_max_nm=fit_cfg.fwhm_max / GaussianFit1D.FWHM_TO_SIGMA,
        sigma_perp_min_nm=fit_cfg.fwhm_min / GaussianFit1D.FWHM_TO_SIGMA,
        sigma_perp_max_nm=fit_cfg.fwhm_max / GaussianFit1D.FWHM_TO_SIGMA,
        global_maxtime=global_maxtime, global_maxiter=10000, cv_folds=3,
        intelligent_sweep=true, fuse_z_bwd=true)
    results_2d, best_2d, ctx = GaussianFit2D.chain_gaussian_sweep(img, pcfg, ccfg)

    # ━━━ Build comparison plot per N ━━━
    xs_2d = ctx.xs; ys_2d = ctx.ys
    zimg = ctx.zimg
    axctx = ctx.axisctx
    t_all = (ctx.x .- axctx.origin[1]) .* axctx.axis[1] .+ (ctx.y .- axctx.origin[2]) .* axctx.axis[2]
    tbins = range(minimum(t_all), maximum(t_all), length=200)
    prof_2d_data = zeros(length(tbins)-1)
    prof_counts = zeros(Int, length(tbins)-1)
    for i in eachindex(t_all)
        b = clamp(Int(floor((t_all[i] - minimum(t_all)) / max(maximum(t_all) - minimum(t_all), 1e-12) * (length(tbins)-1))) + 1, 1, length(tbins)-1)
        prof_2d_data[b] += ctx.z[i]
        prof_counts[b] += 1
    end
    for b in eachindex(prof_2d_data)
        prof_counts[b] > 0 && (prof_2d_data[b] /= prof_counts[b])
    end
    t_centers = (tbins[1:end-1] .+ tbins[2:end]) ./ 2

    # Get 1D results
    fr_1d = fit_1d.fit_run
    all_results_1d = fr_1d.all_results
    x_1d, y_1d = fr_1d.x, fr_1d.y
    cfg_1d = fr_1d.cfg
    best1d = GaussianFit1D.best_result(GaussianFit1D.FitRunResult(x_1d, y_1d, all_results_1d, cfg_1d, "", nothing, nothing, nothing))

    plot_paths = String[]

    if !skip_plots
        for r in results_2d
            r.success || continue
            n = r.n
            N_label = @sprintf("N=%02d", n)
            println("Generating plot for N=$n...")

            # ── 2D model profile along axis ──
            pred_2d = GaussianFit2D._chain_model_values(ctx.x, ctx.y, r.params, n, axctx, ccfg;
                amp_min=r.amp_min, amp_range=r.amp_range)
            prof_2d_model = zeros(length(tbins)-1)
            for i in eachindex(t_all)
                b = clamp(Int(floor((t_all[i] - minimum(t_all)) / max(maximum(t_all) - minimum(t_all), 1e-12) * (length(tbins)-1))) + 1, 1, length(tbins)-1)
                prof_2d_model[b] += pred_2d[i]
            end
            for b in eachindex(prof_2d_model)
                prof_counts[b] > 0 && (prof_2d_model[b] /= prof_counts[b])
            end

            # ── 2D heatmap ──
            pred_img = zeros(size(zimg))
            for iy in eachindex(ys_2d), ix in eachindex(xs_2d)
                pred_img[iy, ix] = GaussianFit2D._chain_model_values([xs_2d[ix]], [ys_2d[iy]], r.params, n, axctx, ccfg;
                    amp_min=r.amp_min, amp_range=r.amp_range)[1]
            end

            # ── Compute ROI bounds for 2D zoom ──
            roi_rows = [iy for iy in eachindex(ys_2d) if any(ctx.mask[iy, :])]
            roi_cols = [ix for ix in eachindex(xs_2d) if any(ctx.mask[:, ix])]
            roi_xmin = isempty(roi_cols) ? minimum(xs_2d) : xs_2d[minimum(roi_cols)] - 0.5
            roi_xmax = isempty(roi_cols) ? maximum(xs_2d) : xs_2d[maximum(roi_cols)] + 0.5
            roi_ymin = isempty(roi_rows) ? minimum(ys_2d) : ys_2d[minimum(roi_rows)] - 0.5
            roi_ymax = isempty(roi_rows) ? maximum(ys_2d) : ys_2d[maximum(roi_rows)] + 0.5

            # ── 2D data + FWHM ellipses (merged) ──
            z_clims = (quantile(vec(zimg), 0.10), quantile(vec(zimg), 0.995))
            p_top1 = heatmap(xs_2d, ys_2d, zimg; aspect_ratio=:equal,
                             title="2D data + model", xlabel="x (nm)", ylabel="y (nm)", colorbar=false,
                             colormap=:thermal, clims=z_clims)
            xlims!(p_top1, roi_xmin, roi_xmax); ylims!(p_top1, roi_ymin, roi_ymax)
            contour!(p_top1, xs_2d, ys_2d, Float64.(ctx.mask); levels=[0.5], color=:white, linewidth=2)
            ox, oy = axctx.origin; ax, ay = axctx.axis
            plot!(p_top1, [ox + minimum(t_all)*ax, ox + maximum(t_all)*ax], [oy + minimum(t_all)*ay, oy + maximum(t_all)*ay];
                  color=:yellow, linewidth=2, label="axis")
            # Ridge path: for each axial position t, take max pixel value in perpendicular direction
            perp_r = [-ay, ax]; t_ridge = range(minimum(t_all), maximum(t_all), length=120)
            rx2, ry2 = Float64[], Float64[]
            for tv in t_ridge
                bx, by, bv = 0.0, 0.0, -Inf
                for uv in range(-0.35, 0.35, step=0.02)
                    x0, y0 = ox + tv*ax + uv*perp_r[1], oy + tv*ay + uv*perp_r[2]
                    if minimum(xs_2d) <= x0 <= maximum(xs_2d) && minimum(ys_2d) <= y0 <= maximum(ys_2d)
                        v = _bilinear(xs_2d, ys_2d, zimg, x0, y0)
                        if v > bv; bv = v; bx = x0; by = y0; end
                    end
                end
                push!(rx2, bx); push!(ry2, by)
            end
            plot!(p_top1, rx2, ry2; color=:lime, linewidth=2, alpha=0.8, label="ridge")
            if n > 0 && r.success
                b, feats, ts, us, spars, sperps = GaussianFit2D._decode_chain(r.params, n, axctx, ccfg;
                    amp_min=r.amp_min, amp_range=r.amp_range)
                axis_angle = atan(ax, ay)
                for f in feats
                    a_ellipse = f.sigma_x_nm * FWHM_SIGMA / 2
                    b_ellipse = f.sigma_y_nm * FWHM_SIGMA / 2
                    _ellipse!(p_top1, f.x_nm, f.y_nm, a_ellipse, b_ellipse, axis_angle; color=:cyan, alpha=0.5, label="")
                end
            end

            # ── 2D residuals ──
            resid_img = (zimg .- pred_img) .* Float64.(ctx.mask) ./ max(ctx.noise, 1e-12)
            p_res2d = heatmap(xs_2d, ys_2d, resid_img; aspect_ratio=:equal,
                              title="2D residuals / noise", xlabel="x (nm)", ylabel="y (nm)", colorbar=false,
                              colormap=COLORMAP_RESID, clims=(-3, 3))
            xlims!(p_res2d, roi_xmin, roi_xmax); ylims!(p_res2d, roi_ymin, roi_ymax)

            # ── 1D slide + fit (top-right) ──
            t_shift_1d = 0.0
            if n > 0 && r.success
                b, feats, ts, us, spars, sperps = GaussianFit2D._decode_chain(r.params, n, axctx, ccfg;
                    amp_min=r.amp_min, amp_range=r.amp_range)
                t_shift_1d = ts[1]  # first 2D lobe position as reference
            end
            x_1d_t = x_1d .+ t_shift_1d

            # ── 1D slide + fit ──
            # Find the 1D result for this N, or use the closest
            r1d = nothing
            for r1 in all_results_1d
                if r1.n_peaks == n
                    r1d = r1
                    break
                end
            end

            if r1d !== nothing
                y1d_pred = GaussianFit1D.predict_fit(x_1d, r1d, cfg_1d)
                y1d_resid = y_1d .- y1d_pred
                p_1d = plot(x_1d_t, y_1d; color=:gray, alpha=0.7, label="1D data", linewidth=1)
                plot!(p_1d, x_1d_t, y1d_pred; color=:red, label="1D fit $N_label", linewidth=2)
                centers = GaussianFit1D._params_to_centers(r1d.popt, n)
                comp_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :lime, :teal, :gold]
                asymmetric = cfg_1d.asymmetric_edges && n >= 2
                y0 = r1d.popt[1]
                for (i, c) in enumerate(centers)
                    idx = i - 1
                    A = GaussianFit1D._get_amplitude(r1d.popt, idx)
                    σ_in = GaussianFit1D._get_sigma(r1d.popt, idx)
                    if asymmetric && (idx == 0 || idx == n - 1)
                        σ_out = idx == 0 ? r1d.popt[end-1] : r1d.popt[end]
                        z = x_1d .- c
                        s = idx == 0 ? (z .< 0) .* σ_out .+ (z .>= 0) .* σ_in :
                                       (z .< 0) .* σ_in .+ (z .>= 0) .* σ_out
                        y_comp = y0 .+ A .* exp.(-0.5 .* (z ./ s).^2)
                    else
                        y_comp = y0 .+ A .* exp.(-0.5 .* ((x_1d .- c) ./ max(σ_in, 1e-9)).^2)
                    end
                    col = comp_colors[mod1(i, length(comp_colors))]
                    plot!(p_1d, x_1d_t, y_comp; color=col, alpha=0.35, linestyle=:dash, linewidth=1,
                          label=(i==1 ? "components" : ""))
                end
                xlabel!("position (nm)"); ylabel!(p_1d, "intensity")
                dbic1d = round(r1d.student_bic - best1d.student_bic, digits=0)
                dbic1d_str = dbic1d == 0 ? "ΔsBIC=0 (best)" : @sprintf("ΔsBIC=%+.0f", dbic1d)
                title!(p_1d, "1D fit $N_label  $dbic1d_str  sBIC=$(round(r1d.student_bic, digits=0))  resid_σ=$(round(std(y1d_resid), digits=5))")
                # 1D residuals (bottom-right)
                p_res1d = plot(x_1d_t, y1d_resid; color=:red, label="1D resid", linewidth=1)
                hline!(p_res1d, [0]; color=:gray, linestyle=:dash, label="")
                xlabel!("position (nm)"); ylabel!(p_res1d, "residual")
                title!(p_res1d, "1D residuals  σ=$(round(std(y1d_resid), digits=5))")
            else
                best1d_local = GaussianFit1D.best_result(GaussianFit1D.FitRunResult(x_1d, y_1d, all_results_1d, cfg_1d, "", nothing, nothing, nothing))
                y1d_best = GaussianFit1D.predict_fit(x_1d, best1d_local, cfg_1d)
                y1d_resid = y_1d .- y1d_best
                p_1d = plot(x_1d_t, y_1d; color=:gray, alpha=0.7, label="1D data", linewidth=1)
                plot!(p_1d, x_1d_t, y1d_best; color=:red, alpha=0.5, label="1D best (N=$(best1d_local.n_peaks))", linewidth=2)
                xlabel!("position (nm)"); ylabel!(p_1d, "intensity")
                title!(p_1d, "1D best=$(best1d_local.n_peaks)  (no 1D fit for N=$n)")
                p_res1d = plot(x_1d_t, y1d_resid; color=:red, alpha=0.5, label="1D resid best", linewidth=1)
                hline!(p_res1d, [0]; color=:gray, linestyle=:dash, label="")
                xlabel!("position (nm)"); ylabel!(p_res1d, "residual")
                title!(p_res1d, "1D residuals  σ=$(round(std(y1d_resid), digits=5))")
            end

            # ── Metrics annotation (ΔBIC, best = 0) ──
            dbic_2d = r.bic - best_2d.bic
            dbic_label = isapprox(dbic_2d, 0, atol=0.1) ? "Δ=0" : @sprintf("Δ=%.0f", dbic_2d)
            title_str = @sprintf("2D chain vs 1D slide — N=%d\nsBIC=%d  %s  chi2=%.2f  spacing=%.3f  spar=%.3f  sperp=%.3f  %s",
                               n, round(Int, r.bic), dbic_label, r.chi2_reduced,
                               r.mean_spacing_nm, r.sigma_parallel_nm, r.sigma_perp_nm,
                               r.valid ? "" : "[INVALID]")

            # ── Combine: models top, residuals bottom ──
            fig = plot(p_top1, p_1d, p_res2d, p_res1d;
                       layout=(2, 2), size=(1800, 1200),
                       plot_title=title_str,
                       plot_titlefontsize=10)

            # Save
            out = joinpath(output_dir, @sprintf("comparison_N%02d.png", n))
            savefig(fig, out)
            push!(plot_paths, out)
            println("  -> $(out)")
        end
    end

    # ── Summary table ──
    # best1d already computed above

    # Build a lookup for 1D results
    d1d = Dict{Int, GaussianFit1D.FitResult}()
    for r1 in all_results_1d
        d1d[r1.n_peaks] = r1
    end

    # Collect all N tested by either method
    all_N = sort(unique(vcat([r.n for r in results_2d if r.success], collect(keys(d1d)))))

    println("\n" * "="^90)
    println("COMPARAISON 1D vs 2D — support $(round(slide.support_length_nm, digits=2)) nm, spacing [0.35, 0.75] nm")
    println("="^90)
    println(rpad("", 90))
    println(rpad("N", 4) * rpad("ΔsBIC_1D", 12) * rpad("χ²_1D", 10) *
            rpad("ΔsBIC_2D", 12) * rpad("χ²_2D", 10) * rpad("spacing_2D", 12) *
            rpad("valid_2D", 10) * rpad("σ∥_2D", 10) * rpad("σ⊥_2D", 10))
    println(repeat("-", 88))

    for n in all_N
        r1 = get(d1d, n, nothing)
        r2 = findfirst(r -> r.n == n && r.success, results_2d)
        r2 = r2 === nothing ? nothing : results_2d[r2]

        db1 = r1 !== nothing ? @sprintf("%.0f", r1.student_bic - best1d.student_bic) : "--"
        ch1 = r1 !== nothing ? @sprintf("%.2f", r1.chi2_red) : "--"
        db2 = r2 !== nothing ? @sprintf("%.0f", r2.bic - best_2d.bic) : "--"
        ch2 = r2 !== nothing ? @sprintf("%.2f", r2.chi2_reduced) : "--"
        sp2 = r2 !== nothing && r2.success ? @sprintf("%.3f", r2.mean_spacing_nm) : "--"
        vl2 = r2 !== nothing ? (r2.valid ? "✓" : "✗") : "--"
        sP2 = r2 !== nothing && r2.success ? @sprintf("%.3f", r2.sigma_parallel_nm) : "--"
        sT2 = r2 !== nothing && r2.success ? @sprintf("%.3f", r2.sigma_perp_nm) : "--"

        prefix = ((r1 !== nothing && r1 === best1d) || (r2 !== nothing && r2 === best_2d)) ? "→" : " "
        println(rpad(prefix * string(n), 4) * rpad(db1, 12) * rpad(ch1, 10) *
                rpad(db2, 12) * rpad(ch2, 10) * rpad(sp2, 12) * rpad(vl2, 8) *
                rpad(sP2, 10) * rpad(sT2, 10))
    end

    println()
    println("→ Best 1D: N=$(best1d.n_peaks)  (ΔsBIC = 0)")
    println("→ Best 2D: N=$(best_2d.n)  (ΔsBIC = 0)  valid=$(best_2d.valid)")
    println("  Consensus: N=$(best_2d.n)")
    println("\n1D: $(length(x_1d)) points, noise_1d=$(round(slide.noise_1d, digits=6))")
    println("2D: $(length(ctx.z)) pixels in fit mask, noise=$(round(ctx.noise, digits=4)), Z fwd+bwd fusion")
    println("\nPlots saved in: $output_dir")

    return (plot_paths=plot_paths, best_1d=best1d, best_2d=best_2d, results_2d=results_2d, all_results_1d=all_results_1d)
end

end

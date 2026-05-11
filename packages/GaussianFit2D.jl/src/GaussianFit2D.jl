module GaussianFit2D

using ArgParse
using LinearAlgebra
using LsqFit
using Optimization
using OptimizationNLopt
using Plots
using Printf
using Statistics
using STMFitCore: effective_spacing_min, endpoint_overrun

export SXMImage, SXMChannel, PatternConfig, PatternFitResult, MolecularFeature,
       MolecularChain
export read_sxm, channel_names, get_channel, preprocess_channel, detect_blobs,
       fit_molecular_pattern, detect_molecular_chains, main_cli
export chain_gaussian_sweep, chain_direct_fit, fit_chain_consensus, fit_chain_batch

const EPS = 1e-12

struct SXMChannel
    name::String
    unit::String
    direction::String
    data::Matrix{Float64} # [y, x], in physical unit stored by Nanonis
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

struct MolecularFeature
    amplitude::Float64
    x_nm::Float64
    y_nm::Float64
    sigma_x_nm::Float64
    sigma_y_nm::Float64
    score::Float64
end

struct MolecularChain
    id::Int
    features::Vector{MolecularFeature}
    score::Float64
    mean_spacing_nm::Float64
    spacing_cv::Float64
    max_turn_angle_deg::Float64
end

Base.@kwdef mutable struct ChainSweepConfig
    n_min::Int = 2
    n_max::Int = 14
    multistart::Int = 20
    cv_folds::Int = 5
    rng_seed::Int = 4321
    student_nu::Float64 = 4.0
    spacing_min_nm::Float64 = 0.35
    spacing_max_nm::Float64 = 0.75
    lateral_max_nm::Float64 = 0.35
    fit_width_nm::Float64 = 0.45
    support_baseline_quantile::Float64 = 0.10
    support_threshold_fraction::Float64 = 0.25
    support_noise_k::Float64 = 2.5
    support_padding_nm::Float64 = 0.20
    support_min_length_nm::Float64 = 1.0
    t_min_nm::Float64 = NaN
    t_max_nm::Float64 = NaN
    sigma_parallel_min_nm::Float64 = 0.12
    sigma_parallel_max_nm::Float64 = 0.45
    sigma_perp_min_nm::Float64 = 0.10
    sigma_perp_max_nm::Float64 = 0.55
    # 1D bootstrap / external initialization
    init_centers_t::Vector{Float64} = Float64[]
    init_amplitudes::Vector{Float64} = Float64[]
    init_laterals::Vector{Float64} = Float64[]
    init_sigma_parallel::Float64 = NaN
    init_sigma_perp::Float64 = NaN
    skip_global::Bool = false
    boot_sweep_halfwidth::Int = 2
    intelligent_sweep::Bool = true
    early_stop_patience::Int = 3
    early_stop_dbic::Float64 = 100.0
    fuse_z_bwd::Bool = true
    residual_peak_snr_threshold::Float64 = 3.5
    max_overlap::Float64 = 0.60
    min_amplitude_fraction::Float64 = 0.3  # reject models with any peak amplitude < 30% of max (matches 1D)
    max_iter::Int = 300
    global_maxtime::Float64 = 15.0
    global_maxiter::Int = 5000
    global_tol::Float64 = 1e-5
    chain_circular_sigmas::Bool = false  # circular gaussians (spar=sperp per peak, fewer params)
    selection_criterion::String = "bic"  # model selection: "bic" | "aicc" | "cv"
end

Base.@kwdef mutable struct ChainModelResult
    n::Int = 0
    params::Vector{Float64} = Float64[]
    success::Bool = false
    train_nll::Float64 = Inf
    cv_nll_mean::Float64 = Inf
    cv_nll_std::Float64 = Inf
    bic::Float64 = Inf
    aicc::Float64 = Inf
    mdl::Float64 = Inf
    residual_peak_snr::Float64 = Inf
    mean_spacing_nm::Float64 = Inf
    spacing_cv::Float64 = Inf
    max_lateral_nm::Float64 = Inf
    sigma_parallel_nm::Float64 = Inf
    sigma_perp_nm::Float64 = Inf
    overlap::Float64 = Inf
    endpoint_overrun_nm::Float64 = Inf
    bound_like::Int = 0
    valid::Bool = false
    reason::String = ""
    rss::Float64 = Inf
    chi2_reduced::Float64 = Inf
    mad::Float64 = Inf
    amp_min::Float64 = NaN
    amp_range::Float64 = NaN
end

Base.@kwdef mutable struct PatternConfig
    filepath::String = ""
    channel::String = "Z"
    direction::String = "fwd"
    output_dir::String = "results"

    # Preprocessing
    stride::Int = 2
    flatten::String = "plane+rows" # none | plane | rows | plane+rows
    smooth_radius_px::Int = 1
    ignore_border_px::Int = 8

    # Detection
    contrast::String = "auto" # auto | bright | dark
    max_features::Int = 50
    max_fit_features::Int = 8
    min_features::Int = 1
    threshold_sigma::Float64 = 2.5
    min_distance_px::Int = 10

    # Chain validation. All image channels/directions are views of one measure;
    # fusion + chain validation is the main recognition mode.
    fusion::Bool = true
    fusion_channels::String = "Z" # comma-separated channel names, or "all"
    chain_min_length::Int = 3
    chain_min_spacing_nm::Float64 = 0.35
    chain_max_spacing_nm::Float64 = 1.6
    chain_max_angle_deg::Float64 = 45.0
    chain_max_spacing_cv::Float64 = 0.45
    min_chain_score::Float64 = 5.0
    max_chains::Int = 8
    max_path_branches::Int = 6

    # Molecule ROI / principal-axis mode. The views are measurements of the same
    # molecule, so first isolate the molecular object in topography, then search
    # features inside that object instead of chasing background chains.
    roi::Bool = true
    roi_channel::String = "Z"
    roi_threshold_fraction::Float64 = 0.35
    roi_dilate_px::Int = 10
    axis_profile::Bool = true
    axis_peak_sigma::Float64 = 1.5
    axis_min_peak_distance_nm::Float64 = 0.45
    repeat_spacing_nm::Float64 = 0.52
    end_extension_nm::Float64 = 0.35

    # Fit constraints
    min_sigma_nm::Float64 = 0.05
    max_sigma_nm::Float64 = 1.5
    initial_sigma_nm::Union{Nothing,Float64} = nothing
    max_iter::Int = 300
    no_fit::Bool = false
    no_plot::Bool = false
    all_images::Bool = false
    inspect_images::Bool = false
end

Base.@kwdef mutable struct PatternFitResult
    features::Vector{MolecularFeature} = MolecularFeature[]
    raw_features::Vector{MolecularFeature} = MolecularFeature[]
    chains::Vector{MolecularChain} = MolecularChain[]
    evidence_map::Matrix{Float64} = zeros(0, 0)
    roi_mask::BitMatrix = falses(0, 0)
    axis_peaks::Vector{MolecularFeature} = MolecularFeature[]
    evidence_unit::String = "SNR"
    roi_length_nm::Float64 = NaN
    estimated_repeats::Int = 0
    estimated_repeat_range::Tuple{Int,Int} = (0, 0)
    params_unconstrained::Vector{Float64} = Float64[]
    success::Bool = false
    rss::Float64 = Inf
    r_squared::Float64 = -Inf
    bic::Float64 = Inf
    warnings::Vector{String} = String[]
    parameter_file::Union{String,Nothing} = nothing
    plot_file::Union{String,Nothing} = nothing
end

# -----------------------------------------------------------------------------
# SXM reader
# -----------------------------------------------------------------------------

function _parse_header(header_text::String)
    header = Dict{String,String}()
    current = nothing
    buf = IOBuffer()
    for line in split(header_text, '\n')
        if startswith(line, ":") && endswith(strip(line), ":")
            if current !== nothing
                header[current] = strip(String(take!(buf)))
            end
            current = strip(line, [':', ' ', '\r', '\n', '\t'])
        elseif current !== nothing
            println(buf, line)
        end
    end
    if current !== nothing
        header[current] = strip(String(take!(buf)))
    end
    return header
end

function _parse_pair(value::String, T=Float64)
    vals = split(strip(value))
    length(vals) >= 2 || error("Expected two values, got: $value")
    return parse(T, vals[1]), parse(T, vals[2])
end

function _parse_data_info(value::String)
    rows = split(value, '\n')
    infos = NamedTuple[]
    for row in rows
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
    marker = Vector{UInt8}(":SCANIT_END:")
    idx = findfirst(marker, bytes)
    idx === nothing && error("SXM marker :SCANIT_END: not found")
    header = _parse_header(String(bytes[1:(first(idx)-1)]))

    width, height = _parse_pair(header["SCAN_PIXELS"], Int)
    rx, ry = _parse_pair(header["SCAN_RANGE"], Float64)
    ox, oy = _parse_pair(get(header, "SCAN_OFFSET", "0 0"), Float64)
    range_nm = (rx * 1e9, ry * 1e9)
    offset_nm = (ox * 1e9, oy * 1e9)

    expanded = Tuple{String,String,String}[]
    for info in _parse_data_info(header["DATA_INFO"])
        dirs = lowercase(info.direction) == "both" ? ["fwd", "bwd"] : [lowercase(info.direction)]
        for dir in dirs
            push!(expanded, (info.name, info.unit, dir))
        end
    end

    nvals = width * height * length(expanded)
    data_offset = length(bytes) - 4nvals + 1
    data_offset > 0 || error("Invalid SXM data size")
    vals = _read_be_float32(bytes, data_offset, nvals)

    channels = SXMChannel[]
    k = 1
    for (name, unit, dir) in expanded
        raw = vals[k:(k + width * height - 1)]
        k += width * height
        mat = permutedims(reshape(raw, width, height))
        # Nanonis stores backward scans in acquisition order. Flip x so fwd/bwd
        # share the same spatial coordinates. For the example file, Z fwd vs
        # bwd correlation changes from ≈ -0.35 to ≈ 0.96 after this flip.
        lowercase(dir) == "bwd" && (mat = reverse(mat; dims=2))
        push!(channels, SXMChannel(name, unit, dir, mat))
    end
    return SXMImage(filepath, header, width, height, range_nm, offset_nm, channels)
end

channel_names(img::SXMImage) = ["$(c.name) ($(c.direction), $(c.unit))" for c in img.channels]

function get_channel(img::SXMImage, name::String; direction::String="fwd")
    lname, ldir = lowercase(name), lowercase(direction)
    for c in img.channels
        lowercase(c.name) == lname && lowercase(c.direction) == ldir && return c
    end
    error("Channel '$name' direction '$direction' not found. Available: $(join(channel_names(img), ", "))")
end

function _value_scale(unit::String)
    u = lowercase(strip(unit))
    u == "m" && return 1e9, "nm"
    u == "a" && return 1e12, "pA"
    return 1.0, isempty(unit) ? "a.u." : unit
end

_coordinate_vectors(img::SXMImage; stride::Int=1) = (
    collect(range(0, img.range_nm[1], length=img.width)[1:stride:end]),
    collect(range(0, img.range_nm[2], length=img.height)[1:stride:end]),
)

# -----------------------------------------------------------------------------
# Preprocessing and blob detection
# -----------------------------------------------------------------------------

function _plane_fit(xs, ys, z::Matrix{Float64})
    xflat = repeat(xs, inner=length(ys))
    yflat = repeat(ys, outer=length(xs))
    zflat = vec(z) # z[y,x] aligned with xflat/yflat: y varies fastest
    coeff = hcat(ones(length(xflat)), xflat, yflat) \ zflat
    return [coeff[1] + coeff[2] * x + coeff[3] * y for y in ys, x in xs]
end

function _row_median_flatten(z::Matrix{Float64})
    out = copy(z)
    global_med = median(vec(out))
    for iy in axes(out, 1)
        out[iy, :] .-= median(view(out, iy, :)) - global_med
    end
    return out
end

function _box_smooth(z::Matrix{Float64}, radius::Int)
    radius <= 0 && return copy(z)
    out = similar(z)
    ny, nx = size(z)
    for iy in 1:ny, ix in 1:nx
        ylo, yhi = max(1, iy-radius), min(ny, iy+radius)
        xlo, xhi = max(1, ix-radius), min(nx, ix+radius)
        out[iy, ix] = mean(@view z[ylo:yhi, xlo:xhi])
    end
    return out
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
            for yy in max(1, y-1):min(ny, y+1), xx in max(1, x-1):min(nx, x+1)
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
            out[max(1, iy-radius):min(ny, iy+radius), max(1, ix-radius):min(nx, ix+radius)] .= true
        end
    end
    return out
end

function preprocess_channel(img::SXMImage, ch::SXMChannel, cfg::PatternConfig)
    scale, scaled_unit = _value_scale(ch.unit)
    stride = max(1, cfg.stride)
    xs, ys = _coordinate_vectors(img; stride=stride)
    z = ch.data[1:stride:end, 1:stride:end] .* scale
    raw = copy(z)

    if occursin("plane", lowercase(cfg.flatten))
        z .-= _plane_fit(xs, ys, z)
    end
    if occursin("rows", lowercase(cfg.flatten))
        z = _row_median_flatten(z)
    end
    z_smooth = _box_smooth(z, cfg.smooth_radius_px)
    noise = 1.4826 * median(abs.(vec(z_smooth) .- median(vec(z_smooth))))
    noise = max(noise, std(vec(z_smooth)) * 0.1, EPS)
    return xs, ys, raw, z, z_smooth, scaled_unit, noise
end

function molecule_roi_mask(img::SXMImage, cfg::PatternConfig)
    ch = get_channel(img, cfg.roi_channel; direction="fwd")
    xs, ys, _raw, z, z_smooth, _unit, _noise = preprocess_channel(img, ch, cfg)
    signal = z_smooth .- minimum(z_smooth)
    maxsig = maximum(signal)
    if maxsig <= EPS
        return xs, ys, trues(size(z_smooth))
    end
    mask = signal .>= cfg.roi_threshold_fraction * maxsig
    mask = _largest_component(mask)
    mask = _dilate_mask(mask, max(0, cfg.roi_dilate_px ÷ max(1, cfg.stride)))
    return xs, ys, mask
end

function _inside_mask(f::MolecularFeature, xs, ys, mask::BitMatrix)
    ix = clamp(searchsortedfirst(xs, f.x_nm), 1, length(xs))
    iy = clamp(searchsortedfirst(ys, f.y_nm), 1, length(ys))
    return mask[iy, ix]
end

function _candidate_score(v::Float64, contrast::String)
    contrast == "bright" && return v
    contrast == "dark" && return -v
    return abs(v)
end

function detect_blobs(z::Matrix{Float64}, xs, ys, cfg::PatternConfig, noise::Float64)
    contrast = lowercase(cfg.contrast)
    if contrast == "auto"
        hi = quantile(vec(z), 0.995)
        lo = quantile(vec(z), 0.005)
        contrast = abs(hi) >= abs(lo) ? "bright" : "dark"
    end

    ny, nx = size(z)
    border = max(1, cfg.ignore_border_px ÷ max(1, cfg.stride))
    min_dist = max(1, cfg.min_distance_px ÷ max(1, cfg.stride))
    taken = falses(ny, nx)
    candidates = MolecularFeature[]
    threshold = cfg.threshold_sigma * noise

    for iy in (1+border):(ny-border), ix in (1+border):(nx-border)
        taken[iy, ix] && continue
        val = z[iy, ix]
        score = _candidate_score(val, contrast)
        score < threshold && continue
        ylo, yhi = max(1, iy-min_dist), min(ny, iy+min_dist)
        xlo, xhi = max(1, ix-min_dist), min(nx, ix+min_dist)
        local_scores = [_candidate_score(z[yy, xx], contrast) for yy in ylo:yhi, xx in xlo:xhi]
        score < maximum(local_scores) && continue
        amp = contrast == "dark" ? -abs(val) : abs(val)
        push!(candidates, MolecularFeature(amp, xs[ix], ys[iy], NaN, NaN, score / noise))
        taken[ylo:yhi, xlo:xhi] .= true
    end

    sort!(candidates; by=f -> f.score, rev=true)
    return candidates[1:min(length(candidates), cfg.max_features)]
end

function _channel_snr_map(img::SXMImage, ch::SXMChannel, cfg::PatternConfig)
    xs, ys, _raw, _z, z_smooth, _unit, noise = preprocess_channel(img, ch, cfg)
    if lowercase(cfg.contrast) == "bright"
        snr = z_smooth ./ noise
    elseif lowercase(cfg.contrast) == "dark"
        snr = .-z_smooth ./ noise
    else
        # Same object seen through different channels may invert contrast.
        # Use absolute robust SNR for fusion; chain geometry will reject isolated noise.
        snr = abs.(z_smooth ./ noise)
    end
    return xs, ys, snr
end

function fused_evidence_map(img::SXMImage, cfg::PatternConfig)
    evidence = nothing
    xs_final = Float64[]
    ys_final = Float64[]
    used = 0
    requested = Set(strip.(lowercase.(split(cfg.fusion_channels, ','))))
    for ch in img.channels
        if !("all" in requested) && !(lowercase(ch.name) in requested)
            continue
        end
        xs, ys, snr = _channel_snr_map(img, ch, cfg)
        if evidence === nothing
            evidence = zeros(size(snr))
            xs_final, ys_final = xs, ys
        end
        size(snr) == size(evidence) || continue
        evidence .= max.(evidence, snr)
        used += 1
    end
    used > 0 || error("No channels matching fusion_channels='$(cfg.fusion_channels)' available")
    return xs_final, ys_final, evidence::Matrix{Float64}
end

_dist(a::MolecularFeature, b::MolecularFeature) = hypot(a.x_nm - b.x_nm, a.y_nm - b.y_nm)

function _turn_angle_deg(a::MolecularFeature, b::MolecularFeature, c::MolecularFeature)
    v1 = [b.x_nm - a.x_nm, b.y_nm - a.y_nm]
    v2 = [c.x_nm - b.x_nm, c.y_nm - b.y_nm]
    n1, n2 = norm(v1), norm(v2)
    (n1 < EPS || n2 < EPS) && return 180.0
    cosang = clamp(dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return acosd(cosang)
end

function _component_indices(adj::Vector{Vector{Int}})
    n = length(adj)
    seen = falses(n)
    comps = Vector{Int}[]
    for i in 1:n
        seen[i] && continue
        stack = [i]
        seen[i] = true
        comp = Int[]
        while !isempty(stack)
            v = pop!(stack)
            push!(comp, v)
            for u in adj[v]
                if !seen[u]
                    seen[u] = true
                    push!(stack, u)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

function _order_component_by_pca(features::Vector{MolecularFeature})
    length(features) <= 2 && return features
    xs = [f.x_nm for f in features]
    ys = [f.y_nm for f in features]
    X = hcat(xs .- mean(xs), ys .- mean(ys))
    _, _, V = svd(X; full=false)
    axis = V[:, 1]
    proj = X * axis
    return features[sortperm(proj)]
end

function _chain_stats(features::Vector{MolecularFeature})
    d = [_dist(features[i], features[i+1]) for i in 1:(length(features)-1)]
    mean_d = mean(d)
    cv = length(d) <= 1 ? 0.0 : std(d) / max(mean_d, EPS)
    angles = length(features) <= 2 ? [0.0] : [_turn_angle_deg(features[i], features[i+1], features[i+2]) for i in 1:(length(features)-2)]
    return mean_d, cv, maximum(angles)
end

function _make_chain(id::Int, features::Vector{MolecularFeature}, cfg::PatternConfig)
    mean_d, cv, max_angle = _chain_stats(features)
    snr_score = sum(f.score for f in features)
    spacing_penalty = 10cv
    angle_penalty = max_angle / max(cfg.chain_max_angle_deg, EPS)
    # Favour coherent long chains: true weak chain members should not be dropped
    # just because a shorter high-SNR subchain has a slightly better average SNR.
    score = snr_score + 2.0 * length(features) - spacing_penalty - angle_penalty
    return MolecularChain(id, features, score, mean_d, cv, max_angle)
end

function extract_chains(candidates::Vector{MolecularFeature}, cfg::PatternConfig)
    n = length(candidates)
    n == 0 && return MolecularChain[]
    adj = [Int[] for _ in 1:n]
    for i in 1:n-1, j in i+1:n
        d = _dist(candidates[i], candidates[j])
        if cfg.chain_min_spacing_nm <= d <= cfg.chain_max_spacing_nm
            push!(adj[i], j)
            push!(adj[j], i)
        end
    end

    candidate_paths = Vector{Int}[]
    max_branches = max(1, cfg.max_path_branches)

    function dfs!(path::Vector{Int}, used::Set{Int})
        last = path[end]
        neigh = filter(u -> !(u in used), adj[last])
        if length(path) >= 2
            prev = path[end-1]
            neigh = filter(u -> _turn_angle_deg(candidates[prev], candidates[last], candidates[u]) <= cfg.chain_max_angle_deg, neigh)
        end
        sort!(neigh; by=u -> candidates[u].score, rev=true)
        length(neigh) > max_branches && (neigh = neigh[1:max_branches])

        if isempty(neigh)
            length(path) >= cfg.chain_min_length && push!(candidate_paths, copy(path))
            return
        end
        for u in neigh
            push!(path, u)
            push!(used, u)
            dfs!(path, used)
            delete!(used, u)
            pop!(path)
        end
    end

    for start in 1:n
        dfs!(Int[start], Set([start]))
    end

    seen = Set{String}()
    chains = MolecularChain[]
    for path in candidate_paths
        # Canonicalize reversed duplicates but preserve geometric order for stats.
        key = join(min(path, reverse(path)), ",")
        key in seen && continue
        push!(seen, key)
        ch = _make_chain(length(chains) + 1, candidates[path], cfg)
        if ch.spacing_cv <= cfg.chain_max_spacing_cv && ch.score >= cfg.min_chain_score
            push!(chains, ch)
        end
    end

    sort!(chains; by=c -> c.score, rev=true)

    # Greedy non-redundant selection: keep the best paths without reusing nodes.
    selected = MolecularChain[]
    occupied_features = MolecularFeature[]
    for ch in chains
        overlaps = any(_dist(f, g) < cfg.axis_min_peak_distance_nm for f in ch.features for g in occupied_features)
        overlaps && continue
        push!(selected, MolecularChain(length(selected) + 1, ch.features, ch.score,
                                       ch.mean_spacing_nm, ch.spacing_cv,
                                       ch.max_turn_angle_deg))
        append!(occupied_features, ch.features)
        length(selected) >= cfg.max_chains && break
    end
    return selected
end

function detect_molecular_chains(img::SXMImage, cfg::PatternConfig)
    xs, ys, evidence = fused_evidence_map(img, cfg)
    roi_mask = trues(size(evidence))
    if cfg.roi
        _rxs, _rys, roi_mask = molecule_roi_mask(img, cfg)
        size(roi_mask) == size(evidence) || (roi_mask = trues(size(evidence)))
        evidence = evidence .* roi_mask
    end
    candidates = detect_blobs(evidence, xs, ys, cfg, 1.0)
    cfg.roi && (candidates = [f for f in candidates if _inside_mask(f, xs, ys, roi_mask)])
    axis_peaks = cfg.axis_profile ? axis_profile_peaks(xs, ys, evidence, roi_mask, cfg) : MolecularFeature[]
    roi_length, estimated_n, estimated_range = estimate_repeat_count(xs, ys, evidence, roi_mask, cfg)
    # Axis peaks are not extra truth; they seed weak/merged lobes on the same molecule.
    if !isempty(axis_peaks)
        append!(candidates, axis_peaks)
        sort!(candidates; by=f -> f.score, rev=true)
        candidates = candidates[1:min(length(candidates), cfg.max_features)]
    end
    chains = extract_chains(candidates, cfg)
    accepted = MolecularFeature[]
    for chain in chains
        append!(accepted, chain.features)
    end
    return xs, ys, evidence, candidates, chains, accepted, roi_mask, axis_peaks, roi_length, estimated_n, estimated_range
end

function axis_profile_peaks(xs, ys, evidence::Matrix{Float64}, roi_mask::BitMatrix, cfg::PatternConfig)
    pts = Tuple{Float64,Float64,Float64}[]
    for iy in eachindex(ys), ix in eachindex(xs)
        roi_mask[iy, ix] || continue
        w = max(evidence[iy, ix], 0.0)
        w > 0 || continue
        push!(pts, (xs[ix], ys[iy], w))
    end
    length(pts) < 5 && return MolecularFeature[]
    ws = [p[3] for p in pts]
    xmean = sum(p[1] * p[3] for p in pts) / sum(ws)
    ymean = sum(p[2] * p[3] for p in pts) / sum(ws)
    X = hcat([p[1] - xmean for p in pts], [p[2] - ymean for p in pts])
    W = Diagonal(ws ./ maximum(ws))
    _, _, V = svd(W * X; full=false)
    axis = V[:, 1]
    ts = X * axis
    tmin, tmax = extrema(ts)
    dt = median(diff(xs))
    nb = max(8, Int(ceil((tmax - tmin) / max(dt, EPS))))
    prof = zeros(nb)
    counts = zeros(nb)
    xacc = zeros(nb)
    yacc = zeros(nb)
    for (idx, p) in enumerate(pts)
        b = clamp(Int(floor((ts[idx] - tmin) / max(tmax - tmin, EPS) * (nb - 1))) + 1, 1, nb)
        prof[b] += p[3]
        counts[b] += 1
        xacc[b] += p[1] * p[3]
        yacc[b] += p[2] * p[3]
    end
    prof ./= max.(counts, 1.0)
    prof_s = _box_smooth(reshape(prof, 1, :), max(1, cfg.smooth_radius_px))[1, :]
    min_sep_bins = max(1, Int(round(cfg.axis_min_peak_distance_nm / max(dt, EPS))))
    peaks = MolecularFeature[]
    thr = cfg.axis_peak_sigma
    for b in 2:(nb-1)
        prof_s[b] >= prof_s[b-1] && prof_s[b] >= prof_s[b+1] || continue
        prof_s[b] >= thr || continue
        if counts[b] > 0 && prof[b] > 0
            x0 = xacc[b] / (prof[b] * counts[b])
            y0 = yacc[b] / (prof[b] * counts[b])
            push!(peaks, MolecularFeature(prof_s[b], x0, y0, NaN, NaN, prof_s[b]))
        end
    end
    sort!(peaks; by=f -> f.score, rev=true)
    # Merge near-duplicates along the axis/spatially.
    kept = MolecularFeature[]
    for p in peaks
        any(_dist(p, q) < cfg.axis_min_peak_distance_nm for q in kept) && continue
        push!(kept, p)
    end
    return kept
end

function estimate_repeat_count(xs, ys, evidence::Matrix{Float64}, roi_mask::BitMatrix, cfg::PatternConfig)
    pts = Tuple{Float64,Float64,Float64}[]
    for iy in eachindex(ys), ix in eachindex(xs)
        roi_mask[iy, ix] || continue
        w = max(evidence[iy, ix], 0.0)
        w > 0 || continue
        push!(pts, (xs[ix], ys[iy], w))
    end
    length(pts) < 5 && return NaN, 0, (0, 0)
    ws = [p[3] for p in pts]
    cenx = sum(p[1] * p[3] for p in pts) / sum(ws)
    ceny = sum(p[2] * p[3] for p in pts) / sum(ws)
    X = hcat([p[1] - cenx for p in pts], [p[2] - ceny for p in pts])
    W = Diagonal(ws ./ maximum(ws))
    _, _, V = svd(W * X; full=false)
    axis = V[:, 1]
    t = X * axis
    # Use robust support instead of visible local maxima. The STM envelope blurs
    # several repeats into broad lobes, so maxima under-count the chain.
    qlo, qhi = quantile(t, 0.02), quantile(t, 0.98)
    apparent_length = qhi - qlo
    center_span = max(apparent_length - 2cfg.end_extension_nm, cfg.repeat_spacing_nm)
    n = max(1, Int(round(center_span / cfg.repeat_spacing_nm)) + 1)
    nlo = max(1, Int(floor((center_span - cfg.repeat_spacing_nm) / cfg.repeat_spacing_nm)) + 1)
    nhi = max(nlo, Int(ceil((center_span + cfg.repeat_spacing_nm) / cfg.repeat_spacing_nm)) + 1)
    return apparent_length, n, (nlo, nhi)
end

function _robust_roi_data(img::SXMImage, cfg::PatternConfig)
    ch = get_channel(img, cfg.roi_channel; direction="fwd")
    xs, ys, _raw, z, z_smooth, _unit, noise = preprocess_channel(img, ch, cfg)
    _rxs, _rys, mask = molecule_roi_mask(img, cfg)
    vals = z[mask]
    z0 = z .- quantile(vals, 0.05)
    xflat = Float64[]; yflat = Float64[]; zflat = Float64[]
    for iy in eachindex(ys), ix in eachindex(xs)
        mask[iy, ix] || continue
        push!(xflat, xs[ix]); push!(yflat, ys[iy]); push!(zflat, z0[iy, ix])
    end
    return xs, ys, z0, mask, xflat, yflat, zflat, max(noise, EPS)
end

function _flatten_roi(z0::AbstractMatrix{Float64}, mask::BitMatrix, xs::AbstractVector{Float64}, ys::AbstractVector{Float64})
    xflat = Float64[]; yflat = Float64[]; zflat = Float64[]
    for iy in eachindex(ys), ix in eachindex(xs)
        mask[iy, ix] || continue
        push!(xflat, xs[ix]); push!(yflat, ys[iy]); push!(zflat, z0[iy, ix])
    end
    return xflat, yflat, zflat
end

function _finalize_chain_result!(r::ChainModelResult, zfit::AbstractVector{Float64}, pred::AbstractVector{Float64},
                                  noise::Float64, n::Int, n_eff::Real, z::AbstractVector{Float64},
                                  xs::AbstractVector{Float64}, ys::AbstractVector{Float64},
                                  zimg::AbstractMatrix{Float64}, xfit::AbstractVector{Float64},
                                  yfit::AbstractVector{Float64}, axisctx, ccfg::ChainSweepConfig)
    r.train_nll = _student_nll(zfit .- pred, noise, ccfg.student_nu) / length(zfit)
    r.cv_nll_mean, r.cv_nll_std = _chain_cv_score(xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg)
    r.residual_peak_snr = _residual_peak_snr(xfit, yfit, zfit, pred, noise)
    pcount = _chain_nparams(n, ccfg)
    full_nll = r.train_nll * length(z)
    r.bic = 2full_nll + pcount * log(n_eff)
    r.aicc = 2full_nll + 2pcount + (2pcount*(pcount+1)) / max(n_eff-pcount-1, 1)
    resid = zfit .- pred
    r.rss = sum(abs2, resid)
    r.chi2_reduced = r.rss / max(1, length(zfit) - pcount) / max(noise^2, EPS)
    r.mad = median(abs.(resid))
    _chain_metrics!(r, axisctx, ccfg)
    r.mdl = 0.0  # backward compat
    checks = String[]
    isfinite(r.cv_nll_mean) || push!(checks, "CV failed")
    r.overlap <= ccfg.max_overlap || push!(checks, "overlap")
    r.endpoint_overrun_nm <= 1e-6 || push!(checks, "outside support")
    r.residual_peak_snr <= ccfg.residual_peak_snr_threshold || push!(checks, "residual high")
    r.valid = isempty(checks)
    r.reason = r.valid ? "ok" : join(checks, ",")
end

function _fused_roi_data(img::SXMImage, cfg::PatternConfig)
    """Fuse Z fwd + Z bwd for SNR boost. The bwd scan is already flipped in read_sxm."""
    # Preprocessing: share coordinate grid, plane fit, row flatten from fwd
    ch_fwd = get_channel(img, cfg.roi_channel; direction="fwd")
    xs, ys, raw_fwd, z_fwd, zs_fwd, _unit, noise_fwd = preprocess_channel(img, ch_fwd, cfg)
    # bwd: only extract and preprocess data, reuse same xs/ys grid
    ch_bwd = get_channel(img, cfg.roi_channel; direction="bwd")
    # Minimal preprocessing on bwd: stride, scale, flatten, smooth
    scale, _ = _value_scale(ch_bwd.unit)
    stride = max(1, cfg.stride)
    z_bwd = ch_bwd.data[1:stride:end, 1:stride:end] .* scale
    if occursin("plane", lowercase(cfg.flatten))
        z_bwd .-= _plane_fit(xs, ys, z_bwd)
    end
    if occursin("rows", lowercase(cfg.flatten))
        z_bwd = _row_median_flatten(z_bwd)
    end
    z_bwd = _box_smooth(z_bwd, cfg.smooth_radius_px)
    noise_bwd = 1.4826 * median(abs.(vec(z_bwd) .- median(vec(z_bwd))))
    noise_bwd = max(noise_bwd, std(vec(z_bwd)) * 0.1, EPS)
    # Fuse: average preprocessed Z fwd and Z bwd
    z_fused = (z_fwd .+ z_bwd) ./ 2.0
    zs_fused = (zs_fwd .+ z_bwd) ./ 2.0
    noise = max(noise_fwd, noise_bwd)
    # ROI mask from fused data
    _rxs, _rys, mask = molecule_roi_mask_fused(img, cfg, zs_fused)
    vals = z_fused[mask]
    z0 = z_fused .- quantile(vals, 0.05)
    xflat, yflat, zflat = _flatten_roi(z0, mask, xs, ys)
    return xs, ys, z0, mask, xflat, yflat, zflat, max(noise, EPS)
end

function _channel_roi_data(img::SXMImage, cfg::PatternConfig, z_mask::BitMatrix, channel_name::String)
    """Load fused fwd+bwd data from an alternate channel, using the same ROI mask as Z."""
    ch_fwd = get_channel(img, channel_name; direction="fwd")
    xs, ys, _, z_fwd, _, _, noise_fwd = preprocess_channel(img, ch_fwd, cfg)
    has_bwd = any(c -> lowercase(c.name) == lowercase(channel_name) && lowercase(c.direction) == "bwd", img.channels)
    if has_bwd
        ch_bwd = get_channel(img, channel_name; direction="bwd")
        _, _, _, z_bwd, _, _, noise_bwd = preprocess_channel(img, ch_bwd, cfg)
        z_fused = (z_fwd .+ z_bwd) ./ 2.0
        noise = max(noise_fwd, noise_bwd)
    else
        z_fused = z_fwd
        noise = noise_fwd
    end
    vals = z_fused[z_mask]
    z0 = z_fused .- quantile(vals, 0.05)
    xflat, yflat, zflat = _flatten_roi(z0, z_mask, xs, ys)
    return xs, ys, z0, z_mask, xflat, yflat, zflat, max(noise, EPS)
end

function _current_evidence_weights(img::SXMImage, cfg::PatternConfig, z_mask::BitMatrix)
    """Compute soft weights [0.3, 1.0] from Current channel evidence, aligned with Z mask.
    Where Current shows strong signal → weight ≈ 1.0. Where Current is flat/baseline → weight ≈ 0.3."""
    has_current = any(c -> lowercase(c.name) == "current", img.channels)
    !has_current && return nothing
    # Load Current fwd+bwd, same preprocessing as Z
    chc_fwd = get_channel(img, "Current"; direction="fwd")
    xs, ys, _, c_fwd, _, _, _ = preprocess_channel(img, chc_fwd, cfg)
    has_cbwd = any(c -> lowercase(c.name) == "current" && lowercase(c.direction) == "bwd", img.channels)
    if has_cbwd
        chc_bwd = get_channel(img, "Current"; direction="bwd")
        _, _, _, c_bwd, _, _, _ = preprocess_channel(img, chc_bwd, cfg)
        c_avg = (c_fwd .+ c_bwd) ./ 2.0
    else
        c_avg = c_fwd
    end
    # Soft sigmoid: normalize Current to [0,1] per-pixel
    cm = c_avg .- quantile(c_avg[z_mask], 0.05)
    c_range = maximum(cm[z_mask]) - minimum(cm[z_mask])
    if c_range <= EPS
        return ones(length(xs) * length(ys))
    end
    cn = clamp.((cm .- minimum(cm[z_mask])) ./ max(c_range, EPS), 0.0, 1.0)
    # Map to weight range [0.3, 1.0]: low current → downweight, high current → full weight
    w_map = 0.3 .+ 0.7 .* cn
    # Flatten for mask
    w = Float64[]
    for iy in eachindex(ys), ix in eachindex(xs)
        z_mask[iy, ix] || continue
        push!(w, w_map[iy, ix])
    end
    return w
end

function molecule_roi_mask_fused(img, cfg::PatternConfig, z_smooth)
    """ROI mask from fused (or single-view) preprocessed data."""
    signal = z_smooth .- minimum(z_smooth)
    maxsig = maximum(signal)
    if maxsig <= EPS
        return nothing, nothing, trues(size(z_smooth))
    end
    mask = signal .>= cfg.roi_threshold_fraction * maxsig
    mask = _largest_component(mask)
    mask = _dilate_mask(mask, max(0, cfg.roi_dilate_px ÷ max(1, cfg.stride)))
    size_z = size(z_smooth)
    # Compute xs, ys from img
    stride = max(1, cfg.stride)
    xv, yv = _coordinate_vectors(img; stride=stride)
    return xv, yv, mask
end

_rsigmoid(t) = 1 / (1 + exp(-clamp(t, -60.0, 60.0)))
_rlogit(u) = log(clamp(u, 1e-6, 1-1e-6) / (1 - clamp(u, 1e-6, 1-1e-6)))
_mad_std(v) = 1.4826 * median(abs.(v .- median(v)))

function _effective_spacing_min_nm(ccfg::ChainSweepConfig)
    """Hard physical center-spacing constraint derived from max_overlap and sigma bounds.

    Adjacent axial spacings are parameterized with this lower bound. Since lateral
    offsets only increase Euclidean peak distance, this guarantees that even at
    the largest allowed sigma the pairwise overlap cannot exceed `max_overlap`.
    """
    sigma_max = ccfg.chain_circular_sigmas ? ccfg.sigma_parallel_max_nm :
        max(ccfg.sigma_parallel_max_nm, ccfg.sigma_perp_max_nm)
    return effective_spacing_min(ccfg.spacing_min_nm, ccfg.spacing_max_nm, sigma_max, ccfg.max_overlap)
end

function _chain_support_length(axisctx)
    return max(axisctx.tmax - axisctx.tmin, 0.0)
end

function _chain_min_span(n::Int, ccfg::ChainSweepConfig)
    return max(n - 1, 0) * _effective_spacing_min_nm(ccfg)
end

function _chain_can_fit_support(n::Int, axisctx, ccfg::ChainSweepConfig)
    n <= 1 && return true
    return _chain_min_span(n, ccfg) <= _chain_support_length(axisctx) + 1e-9
end

function _student_nll(resid, noise::Float64, nu::Float64; weights=nothing)
    r = resid ./ max(noise, EPS)
    vals = 0.5 * (nu + 1) .* log1p.((r.^2) ./ nu)
    if weights !== nothing
        vals .*= weights
    end
    return sum(vals)
end

function _chain_nparams(n::Int, ccfg::Union{ChainSweepConfig,Nothing}=nothing)
    n == 0 && return 1
    n_sigmas = (ccfg !== nothing && ccfg.chain_circular_sigmas) ? n : 2n
    return 1 + n + 1 + (n-1) + n + n_sigmas
end

function _residual_peak_snr(x, y, z, pred, noise)
    return maximum(abs.(z .- pred)) / max(noise, EPS)
end

# -----------------------------------------------------------------------------
# Ordered 2D chain model selection: Gaussian lobes constrained along one axis
# -----------------------------------------------------------------------------

function _weighted_roi_axis(x, y, z)
    w = max.(z .- quantile(z, 0.10), 0.0) .+ EPS
    sw = sum(w)
    ox = sum(x .* w) / sw
    oy = sum(y .* w) / sw
    X = hcat(x .- ox, y .- oy)
    W = Diagonal(w ./ maximum(w))
    _, _, V = svd(W * X; full=false)
    axis = Vector{Float64}(V[:, 1])
    axis ./= max(norm(axis), EPS)
    # Deterministic orientation: increasing y, then x.
    (axis[2] < 0 || (abs(axis[2]) < EPS && axis[1] < 0)) && (axis .*= -1)
    perp = [-axis[2], axis[1]]
    t = (x .- ox) .* axis[1] .+ (y .- oy) .* axis[2]
    qlo, qhi = quantile(t, 0.01), quantile(t, 0.99)
    return (origin=(ox, oy), axis=axis, perp=perp, tmin=qlo, tmax=qhi)
end

function _chain_coordinates(x, y, axisctx)
    t = (x .- axisctx.origin[1]) .* axisctx.axis[1] .+ (y .- axisctx.origin[2]) .* axisctx.axis[2]
    u = (x .- axisctx.origin[1]) .* axisctx.perp[1] .+ (y .- axisctx.origin[2]) .* axisctx.perp[2]
    return t, u
end

function _active_t_support(t, z, ccfg::ChainSweepConfig)
    if length(t) < 5
        return minimum(t), maximum(t), (support_method="full_t_range_too_few_points",)
    end
    nb = max(20, min(200, Int(ceil((maximum(t) - minimum(t)) / max(0.02, ccfg.spacing_min_nm / 8)))))
    prof = zeros(nb); counts = zeros(Int, nb)
    tlo, thi = minimum(t), maximum(t)
    for i in eachindex(t)
        b = clamp(Int(floor((t[i] - tlo) / max(thi - tlo, EPS) * (nb - 1))) + 1, 1, nb)
        prof[b] += z[i]
        counts[b] += 1
    end
    prof ./= max.(counts, 1)
    baseline = quantile(prof, ccfg.support_baseline_quantile)
    peak = maximum(prof)
    low = prof[prof .<= baseline]
    noise = isempty(low) ? _mad_std(prof) : _mad_std(low)
    threshold_contrast = baseline + ccfg.support_threshold_fraction * max(peak - baseline, EPS)
    threshold_noise = baseline + ccfg.support_noise_k * max(noise, EPS)
    thr = max(threshold_contrast, threshold_noise)
    active = findall(i -> counts[i] > 0 && prof[i] >= thr, eachindex(prof))
    base_meta = (support_method="auto_axis_profile_support",
                 profile_bins=nb, baseline=baseline, peak=peak,
                 noise_sigma_profile=noise, threshold=thr,
                 threshold_contrast=threshold_contrast, threshold_noise=threshold_noise,
                 baseline_quantile=ccfg.support_baseline_quantile,
                 threshold_fraction=ccfg.support_threshold_fraction,
                 noise_k=ccfg.support_noise_k, padding_nm=ccfg.support_padding_nm,
                 min_support_nm=ccfg.support_min_length_nm, active_bins=length(active))
    isempty(active) && return tlo, thi, merge(base_meta, (fallback="no_active_bins", active_components_raw=0, active_components_after_min_length=0))
    # Use the active support containing the global peak when available; otherwise largest component.
    runs = UnitRange{Int}[]
    start = active[1]; prev = active[1]
    for idx in active[2:end]
        if idx == prev + 1
            prev = idx
        else
            push!(runs, start:prev)
            start = prev = idx
        end
    end
    push!(runs, start:prev)
    dt = (thi - tlo) / nb
    runs_all = copy(runs)
    runs = filter(r -> length(r) * dt >= ccfg.support_min_length_nm, runs)
    isempty(runs) && return tlo, thi, merge(base_meta, (fallback="no_component_long_enough", active_components_raw=length(runs_all), active_components_after_min_length=0))
    peak_bin = argmax(prof)
    containing = filter(r -> first(r) <= peak_bin <= last(r), runs)
    best = isempty(containing) ? argmax(r -> length(r), runs) : argmax(r -> length(r), containing)
    newlo = tlo + (first(best) - 1) * dt - ccfg.support_padding_nm
    newhi = tlo + last(best) * dt + ccfg.support_padding_nm
    final_lo, final_hi = max(tlo, newlo), min(thi, newhi)
    return final_lo, final_hi, merge(base_meta, (
        fallback="none", active_components_raw=length(runs_all),
        active_components_after_min_length=length(runs), selected_component_start_bin=first(best),
        selected_component_end_bin=last(best), peak_bin=peak_bin,
        peak_t_nm=tlo + (peak_bin - 0.5) * dt,
        support_start_t_nm=final_lo, support_end_t_nm=final_hi,
        support_length_nm=final_hi-final_lo))
end

function _chain_fit_data(x, y, z, axisctx, ccfg::ChainSweepConfig)
    t, u = _chain_coordinates(x, y, axisctx)
    tube = abs.(u) .<= ccfg.fit_width_nm
    sum(tube) >= 20 || return x, y, z, axisctx, trues(length(x)), (support_method="full_roi_tube_too_small",)
    tlo, thi, support_meta = _active_t_support(t[tube], z[tube], ccfg)
    isfinite(ccfg.t_min_nm) && (tlo = max(tlo, ccfg.t_min_nm))
    isfinite(ccfg.t_max_nm) && (thi = min(thi, ccfg.t_max_nm))
    tlo < thi || error("Invalid chain t window: t_min_nm=$(ccfg.t_min_nm), t_max_nm=$(ccfg.t_max_nm), active support=$(tlo)-$(thi)")
    keep = tube .& (t .>= tlo) .& (t .<= thi)
    sum(keep) >= 20 || (keep = tube)
    xf, yf, zf = x[keep], y[keep], z[keep]
    tfit, _ufit = _chain_coordinates(xf, yf, axisctx)
    axisctx_fit = (origin=axisctx.origin, axis=axisctx.axis, perp=axisctx.perp,
                   tmin=tlo, tmax=thi)
    support_meta = merge(support_meta, (user_t_min_nm=ccfg.t_min_nm, user_t_max_nm=ccfg.t_max_nm,
                                        final_t_min_nm=tlo, final_t_max_nm=thi,
                                        final_support_length_nm=thi-tlo,
                                        fit_width_nm=ccfg.fit_width_nm,
                                        tube_pixels=sum(tube), fit_mask_pixels=sum(keep)))
    return xf, yf, zf, axisctx_fit, keep, support_meta
end

function _decode_chain(p::AbstractVector, n::Int, axisctx, ccfg::ChainSweepConfig;
                       amp_min::Float64=NaN, amp_range::Float64=NaN)
    n == 0 && return p[1], MolecularFeature[], Float64[], Float64[], Float64[], Float64[]
    b0 = p[1]
    j = 2
    if !isnan(amp_min) && !isnan(amp_range) && amp_range > EPS
        amps = [amp_min + amp_range * clamp(_rsigmoid(p[j+k-1]), 1e-9, 1.0-1e-9) for k in 1:n]
    else
        amps = [exp(clamp(p[j+k-1], -30, 30)) for k in 1:n]
    end
    j += n
    spacing_min_eff = _effective_spacing_min_nm(ccfg)
    support_len = _chain_support_length(axisctx)
    t0_raw = p[j]; j += 1
    deltas = Float64[]
    used_span = 0.0
    for i in 1:(n-1)
        remaining_after = (n - 1 - i) * spacing_min_eff
        upper_i = min(ccfg.spacing_max_nm, support_len - used_span - remaining_after)
        upper_i = max(upper_i, spacing_min_eff)
        d = spacing_min_eff + (upper_i - spacing_min_eff) * _rsigmoid(p[j])
        push!(deltas, d)
        used_span += d
        j += 1
    end
    t0_hi = axisctx.tmax - used_span
    t0_span = max(t0_hi - axisctx.tmin, EPS)
    t0 = axisctx.tmin + t0_span * _rsigmoid(t0_raw)
    us = Float64[]
    for _ in 1:n
        push!(us, ccfg.lateral_max_nm * tanh(p[j]))
        j += 1
    end
    if ccfg.chain_circular_sigmas
        sigmas = [ccfg.sigma_parallel_min_nm + (ccfg.sigma_parallel_max_nm - ccfg.sigma_parallel_min_nm) * _rsigmoid(p[j+k-1]) for k in 1:n]
        spars = sigmas
        sperps = sigmas
    else
        spars = [ccfg.sigma_parallel_min_nm + (ccfg.sigma_parallel_max_nm - ccfg.sigma_parallel_min_nm) * _rsigmoid(p[j+k-1]) for k in 1:n]
        j += n
        sperps = [ccfg.sigma_perp_min_nm + (ccfg.sigma_perp_max_nm - ccfg.sigma_perp_min_nm) * _rsigmoid(p[j+k-1]) for k in 1:n]
    end
    ts = Float64[t0]
    for d in deltas
        push!(ts, ts[end] + d)
    end
    ox, oy = axisctx.origin
    ax, ay = axisctx.axis
    px, py = axisctx.perp
    feats = MolecularFeature[]
    for k in 1:n
        xk = ox + ts[k] * ax + us[k] * px
        yk = oy + ts[k] * ay + us[k] * py
        push!(feats, MolecularFeature(amps[k], xk, yk, spars[k], sperps[k], amps[k]))
    end
    return b0, feats, ts, us, spars, sperps
end

function _chain_model_values(x, y, p, n::Int, axisctx, ccfg::ChainSweepConfig;
                             amp_min::Float64=NaN, amp_range::Float64=NaN)
    b0, feats, _ts, _us, _spars, _sperps = _decode_chain(p, n, axisctx, ccfg;
                                                          amp_min=amp_min, amp_range=amp_range)
    pred = fill(b0, length(x))
    for f in feats
        @. pred += f.amplitude * exp(-0.5 * ((((x - f.x_nm)*axisctx.axis[1] + (y - f.y_nm)*axisctx.axis[2])/f.sigma_x_nm)^2 + (((x - f.x_nm)*axisctx.perp[1] + (y - f.y_nm)*axisctx.perp[2])/f.sigma_y_nm)^2))
    end
    return pred
end

function _nearest_values_on_grid(xs, ys, zimg, feats::Vector{MolecularFeature})
    vals = Float64[]
    for f in feats
        ix = clamp(searchsortedfirst(xs, f.x_nm), 1, length(xs))
        iy = clamp(searchsortedfirst(ys, f.y_nm), 1, length(ys))
        push!(vals, max(zimg[iy, ix], EPS))
    end
    return vals
end

function _pack_chain_initial(xs, ys, zimg, n::Int, axisctx, ccfg::ChainSweepConfig)
    have_1d_init = length(ccfg.init_centers_t) >= n && length(ccfg.init_amplitudes) >= n
    spacing_min_eff = _effective_spacing_min_nm(ccfg)
    support_len = _chain_support_length(axisctx)
    if !_chain_can_fit_support(n, axisctx, ccfg)
        error(@sprintf("N=%d cannot fit support %.4f nm with effective min spacing %.4f nm", n, support_len, spacing_min_eff))
    end

    p = Float64[quantile(vec(zimg), 0.05)]

    if have_1d_init
        # Use 1D bootstrap centers and amplitudes
        centers_t = ccfg.init_centers_t[1:n]
        amps = ccfg.init_amplitudes[1:n]
        medamp = max(median(amps), EPS)
        raw_deltas0 = n > 1 ? diff(centers_t) : Float64[]
    else
        spacing0 = n > 1 ? clamp(support_len / max(n - 1, 1), spacing_min_eff, ccfg.spacing_max_nm) : spacing_min_eff
        total = spacing0 * max(n - 1, 0)
        t0 = axisctx.tmin + 0.5 * max(support_len - total, 0.0)
        ts0 = [t0 + (k-1) * spacing0 for k in 1:n]
        ox, oy = axisctx.origin; ax, ay = axisctx.axis
        feats0 = [MolecularFeature(1.0, ox + t*ax, oy + t*ay, 0.2, 0.2, 1.0) for t in ts0]
        amps = _nearest_values_on_grid(xs, ys, zimg, feats0)
        medamp = max(median(amps), EPS)
        raw_deltas0 = fill(spacing0, max(n - 1, 0))
    end

    deltas0 = Float64[]
    used_span = 0.0
    for i in 1:(n-1)
        remaining_after = (n - 1 - i) * spacing_min_eff
        upper_i = min(ccfg.spacing_max_nm, support_len - used_span - remaining_after)
        upper_i = max(upper_i, spacing_min_eff)
        raw = i <= length(raw_deltas0) ? raw_deltas0[i] : spacing_min_eff
        d = clamp(raw, spacing_min_eff, upper_i)
        push!(deltas0, d)
        used_span += d
    end
    if have_1d_init
        t0 = centers_t[1]
    end
    total = sum(deltas0)
    t0_hi = axisctx.tmax - total
    t0 = clamp(t0, axisctx.tmin, t0_hi)

    amp_max_data = max(maximum(zimg), EPS)
    amp_min_val = ccfg.min_amplitude_fraction * amp_max_data
    amp_range_val = max(amp_max_data - amp_min_val, EPS)
    for a in amps
        # floor at 2% above amp_min to avoid extreme sigmoid gradients at initialization
        raw_frac = (a - amp_min_val) / amp_range_val
        frac = clamp(raw_frac, 0.02, 0.98)
        push!(p, _rlogit(frac))
    end
    t0u = clamp((t0 - axisctx.tmin) / max(t0_hi - axisctx.tmin, EPS), 1e-5, 1-1e-5)
    push!(p, _rlogit(t0u))
    used_span = 0.0
    for i in 1:(n-1)
        remaining_after = (n - 1 - i) * spacing_min_eff
        upper_i = min(ccfg.spacing_max_nm, support_len - used_span - remaining_after)
        upper_i = max(upper_i, spacing_min_eff)
        push!(p, _rlogit((deltas0[i] - spacing_min_eff) / max(upper_i - spacing_min_eff, EPS)))
        used_span += deltas0[i]
    end
    # lateral offsets
    if have_1d_init && length(ccfg.init_laterals) >= n
        for k in 1:n
            u0 = clamp(ccfg.init_laterals[k] / max(ccfg.lateral_max_nm, EPS), -1+1e-6, 1-1e-6)
            push!(p, atanh(u0))
        end
    else
        for _ in 1:n
            push!(p, 0.0)
        end
    end
    spacing0 = isempty(deltas0) ? spacing_min_eff : mean(deltas0)

    # sigma parallel
    spar0 = if isfinite(ccfg.init_sigma_parallel)
        clamp(ccfg.init_sigma_parallel, ccfg.sigma_parallel_min_nm, ccfg.sigma_parallel_max_nm)
    else
        clamp(0.5 * spacing0, ccfg.sigma_parallel_min_nm, ccfg.sigma_parallel_max_nm)
    end
    # sigma perp
    sperp0 = if isfinite(ccfg.init_sigma_perp)
        clamp(ccfg.init_sigma_perp, ccfg.sigma_perp_min_nm, ccfg.sigma_perp_max_nm)
    else
        clamp(0.35 * spacing0, ccfg.sigma_perp_min_nm, ccfg.sigma_perp_max_nm)
    end
    if ccfg.chain_circular_sigmas
        sigma_trans = _rlogit((spar0 - ccfg.sigma_parallel_min_nm) / max(ccfg.sigma_parallel_max_nm - ccfg.sigma_parallel_min_nm, EPS))
        for _ in 1:n
            push!(p, sigma_trans)
        end
    else
        spar_trans = _rlogit((spar0 - ccfg.sigma_parallel_min_nm) / max(ccfg.sigma_parallel_max_nm - ccfg.sigma_parallel_min_nm, EPS))
        sperp_trans = _rlogit((sperp0 - ccfg.sigma_perp_min_nm) / max(ccfg.sigma_perp_max_nm - ccfg.sigma_perp_min_nm, EPS))
        for _ in 1:n
            push!(p, spar_trans)
        end
        for _ in 1:n
            push!(p, sperp_trans)
        end
    end
    return p, amp_min_val, amp_range_val
end

function _fit_chain_n(xs, ys, zimg, x, y, z, noise, n::Int, axisctx, ccfg::ChainSweepConfig; starts::Int=ccfg.multistart)
    n == 0 && return ChainModelResult(n=0, params=[median(z)], success=true)
    if !_chain_can_fit_support(n, axisctx, ccfg)
        return ChainModelResult(n=n, success=false, valid=false,
            reason=@sprintf("infeasible support: min span %.4f > support %.4f", _chain_min_span(n, ccfg), _chain_support_length(axisctx)))
    end
    p0, amp_min, amp_range = _pack_chain_initial(xs, ys, zimg, n, axisctx, ccfg)
    xy = vcat(reshape(x, 1, :), reshape(y, 1, :))
    model = (xydata, p) -> _chain_model_values(view(xydata,1,:), view(xydata,2,:), p, n, axisctx, ccfg;
                                                amp_min=amp_min, amp_range=amp_range)

    np = _chain_nparams(n, ccfg)
    lower = fill(-10.0, np); upper = fill(10.0, np)
    lower[1] = -5.0; upper[1] = 5.0
    j = 2
    for _ in 1:n; lower[j] = -5.0; upper[j] = 5.0; j += 1 end  # sigmoid-transformed amplitude
    lower[j] = -4.0; upper[j] = 4.0; j += 1  # t0
    for _ in 1:n-1; lower[j] = -5.0; upper[j] = 5.0; j += 1 end  # deltas
    for _ in 1:n; lower[j] = -3.0; upper[j] = 3.0; j += 1 end  # u
    # spars (n) and sperps (n)
    if ccfg.chain_circular_sigmas
        for _ in 1:n
            lower[j] = -5.0; upper[j] = 5.0
            j += 1
        end
    else
        for _ in 1:(2n)
            lower[j] = -5.0; upper[j] = 5.0
            j += 1
        end
    end
    p0 = clamp.(p0, lower .+ 1e-9, upper .- 1e-9)

    p_global = p0
    global_success = true  # default true if we skip NLopt
    if !ccfg.skip_global
        objective = (u, _) -> sum(abs2, z .- model(xy, u))
        nlop = OptimizationNLopt.NLopt.Opt(:GN_DIRECT_L, np)
        nlop.xtol_rel = ccfg.global_tol
        nlop.ftol_rel = ccfg.global_tol
        nlop.maxtime  = ccfg.global_maxtime
        prob = OptimizationProblem(objective, p0; lb=lower, ub=upper)
        global_success = false
        try
            sol = solve(prob, nlop; maxiters=ccfg.global_maxiter)
            p_global = sol.u
            global_success = true
        catch e
            # If NLopt fails but we have good 1D init, fall through to LM
            length(ccfg.init_centers_t) >= n || return ChainModelResult(n=n, success=false, reason="NLopt failed: $e")
        end
    end

    p_final = p_global
    try
        fit = curve_fit(model, xy, z, p_global; lower=lower, upper=upper, maxIter=ccfg.max_iter, autodiff=:finite)
        p_final = fit.param
    catch
    end
    pred = model(xy, p_final)
    nll = _student_nll(z .- pred, noise, ccfg.student_nu)
    return ChainModelResult(n=n, params=p_final, success=true, train_nll=nll,
                            amp_min=amp_min, amp_range=amp_range)
end

function _chain_overlap(feats::Vector{MolecularFeature}, spar::Float64, sperp::Float64)
    length(feats) <= 1 && return 0.0
    s = max(spar, sperp, EPS)
    ov = 0.0
    for i in 1:length(feats)-1, j in i+1:length(feats)
        d = _dist(feats[i], feats[j])
        ov = max(ov, exp(-0.5 * (d/s)^2))
    end
    return ov
end

function _chain_metrics!(r::ChainModelResult, axisctx, ccfg::ChainSweepConfig)
    (!r.success || isempty(r.params)) && return
    if r.n == 0
        r.mean_spacing_nm = Inf; r.spacing_cv = 0.0; r.max_lateral_nm = 0.0
        r.sigma_parallel_nm = NaN; r.sigma_perp_nm = NaN; r.overlap = 0.0; r.endpoint_overrun_nm = 0.0
        return
    end
    _b, feats, ts, us, spars, sperps = _decode_chain(r.params, r.n, axisctx, ccfg;
                                                      amp_min=r.amp_min, amp_range=r.amp_range)
    ds = diff(ts)
    r.mean_spacing_nm = isempty(ds) ? Inf : mean(ds)
    r.spacing_cv = length(ds) <= 1 ? 0.0 : std(ds) / max(mean(ds), EPS)
    r.max_lateral_nm = maximum(abs.(us))
    r.sigma_parallel_nm = mean(spars)
    r.sigma_perp_nm = mean(sperps)
    r.overlap = _chain_overlap(feats, mean(spars), mean(sperps))
    r.endpoint_overrun_nm = endpoint_overrun(ts, axisctx.tmin, axisctx.tmax)
    near(v, lo, hi) = (v - lo) / max(hi - lo, EPS) < 0.03 || (hi - v) / max(hi - lo, EPS) < 0.03
    spacing_min_eff = _effective_spacing_min_nm(ccfg)
    r.bound_like = count(d -> near(d, spacing_min_eff, ccfg.spacing_max_nm), ds)
    for sp in spars
        r.bound_like += near(sp, ccfg.sigma_parallel_min_nm, ccfg.sigma_parallel_max_nm) ? 1 : 0
    end
    for sp in sperps
        r.bound_like += near(sp, ccfg.sigma_perp_min_nm, ccfg.sigma_perp_max_nm) ? 1 : 0
    end
end

function _chain_cv_score(xs, ys, zimg, x, y, z, noise, n, axisctx, ccfg::ChainSweepConfig)
    n == 0 && return _student_nll(z .- median(z), noise, ccfg.student_nu), 0.0
    folds = max(2, ccfg.cv_folds)
    t = (x .- axisctx.origin[1]) .* axisctx.axis[1] .+ (y .- axisctx.origin[2]) .* axisctx.axis[2]
    tmin, tmax = extrema(t)
    scores = Float64[]
    for fold in 1:folds
        val = [mod(floor(Int, (t[i]-tmin)/(tmax-tmin+EPS)*folds), folds) + 1 == fold for i in eachindex(t)]
        train = .!val
        sum(train) > 10 && sum(val) > 5 || continue
        r = _fit_chain_n(xs, ys, zimg, x[train], y[train], z[train], noise, n, axisctx, ccfg; starts=max(3, ccfg.multistart ÷ 4))
        r.success || continue
        pred = _chain_model_values(x[val], y[val], r.params, n, axisctx, ccfg;
                                   amp_min=r.amp_min, amp_range=r.amp_range)
        push!(scores, _student_nll(z[val] .- pred, noise, ccfg.student_nu) / sum(val))
    end
    isempty(scores) && return Inf, Inf
    return mean(scores), length(scores) > 1 ? std(scores)/sqrt(length(scores)) : 0.0
end

function _select_chain_model(results::Vector{ChainModelResult}, ccfg::ChainSweepConfig)
    # Select by the configured criterion among all successful fits
    succ = filter(r -> r.success, results)
    isempty(succ) && error("No chain model fit succeeded")
    criterion = lowercase(ccfg.selection_criterion)
    if criterion == "aicc"
        best = argmin(r -> r.aicc, succ)
    elseif criterion == "cv"
        best = argmin(r -> r.cv_nll_mean, succ)
    else
        best = argmin(r -> r.bic, succ)
    end
    # Warn if the selected model has quality flags
    if !best.valid
        @warn "Selected model N=$(best.n) has quality issues: $(best.reason)"
    end
    return best
end

function chain_gaussian_sweep(img::SXMImage, cfg::PatternConfig, ccfg::ChainSweepConfig;
                               override_data=nothing, override_axisctx=nothing)
    # ── Data loading ──
    if override_data !== nothing
        xs, ys, zimg, mask, x, y, z, noise = override_data
    else
        has_bwd = any(c -> lowercase(c.name) == lowercase(cfg.roi_channel) && lowercase(c.direction) == "bwd", img.channels)
        use_fusion = ccfg.fuse_z_bwd && has_bwd
        if use_fusion
            xs, ys, zimg, mask, x, y, z, noise = _fused_roi_data(img, cfg)
        else
            xs, ys, zimg, mask, x, y, z, noise = _robust_roi_data(img, cfg)
        end
    end
    # ── Axis / support ──
    axisctx_full = override_axisctx !== nothing ? override_axisctx : _weighted_roi_axis(x, y, z)
    xfit, yfit, zfit, axisctx, fit_keep, support_meta = _chain_fit_data(x, y, z, axisctx_full, ccfg)
    # Adaptive range over the hard axial support used by the parameterization.
    axis_length = _chain_support_length(axisctx)
    spacing_min_eff = _effective_spacing_min_nm(ccfg)
    n_max_data = max(1, Int(floor(axis_length / spacing_min_eff)) + 1)
    n_max_eff = min(ccfg.n_max, n_max_data)

    # effective sample size (pixels in fit mask ÷ typical spatial correlation factor)
    n_eff = max(10, length(zfit) ÷ 9)
    results = ChainModelResult[]

    if ccfg.intelligent_sweep
        # ── Intelligent centre-out sweep with early stopping (same logic as 1D) ──
        n_min_data = max(2, Int(floor(axis_length / ccfg.spacing_max_nm)))
        n_min_eff = max(ccfg.n_min, n_min_data)
        if n_min_eff > n_max_eff
            @printf("    -> requested/coverage N_min=%d exceeds feasible N_max=%d for %.2f nm support; clamping to feasible range\n",
                    n_min_eff, n_max_eff, axis_length)
            n_min_eff = n_max_eff
        end
        mean_sp = (spacing_min_eff + ccfg.spacing_max_nm) / 2.0
        center_n = clamp(Int(round(axis_length / mean_sp)), n_min_eff, n_max_eff)
        # User can anchor the center via n_min config (e.g. from 1D)
        if ccfg.n_min > 2 && ccfg.n_min == ccfg.n_max
            center_n = ccfg.n_min  # fixed N, no sweep
            n_min_eff = ccfg.n_min
            n_max_eff = ccfg.n_max
        elseif ccfg.n_min > n_min_eff + 2
            # User suggested a different floor, use it as hint
            center_n = max(center_n, ccfg.n_min)
        end
        if n_min_eff == n_max_eff
            center_n = n_min_eff
        end

        @printf("--- Sweep: N=%d..%d (adaptive from %.1f nm support, spacing %.2f-%.2f; effective min %.2f) ---\n",
                n_min_eff, n_max_eff, axis_length, ccfg.spacing_min_nm, ccfg.spacing_max_nm, spacing_min_eff)
        @printf("    center N=%d\n", center_n)

        # Fit center first
        @printf("  Fitting N=%d...\n", center_n)
        r_center = _fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, center_n, axisctx, ccfg)
        best_bic = Inf
        if r_center.success
            pred = _chain_model_values(xfit, yfit, r_center.params, center_n, axisctx, ccfg;
                                       amp_min=r_center.amp_min, amp_range=r_center.amp_range)
            _finalize_chain_result!(r_center, zfit, pred, noise, center_n, n_eff, z, xs, ys, zimg, xfit, yfit, axisctx, ccfg)
            best_bic = r_center.bic
        end
        push!(results, r_center)

        # Sweep right (N increasing)
        streak = 0
        for n in (center_n+1):n_max_eff
            if streak >= ccfg.early_stop_patience
                @printf("    -> right early stop at N=%d (BIC incr %dx)\n", n, ccfg.early_stop_patience)
                break
            end
            @printf("  Fitting N=%d...\n", n)
            r = _fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg)
            if r.success
                pred = _chain_model_values(xfit, yfit, r.params, n, axisctx, ccfg;
                                           amp_min=r.amp_min, amp_range=r.amp_range)
                _finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, axisctx, ccfg)
                # Immediate stop: model overflows the support (lobes outside visible range)
                if r.endpoint_overrun_nm > ccfg.spacing_max_nm
                    @printf("    -> right stop N=%d (overflows support by %.2f nm)\n", n, r.endpoint_overrun_nm)
                    push!(results, r)
                    break
                end
                # Immediate stop: BIC worsened and we're more than 1 step past centre — excessive lobes
                if r.bic > best_bic && n > center_n + 1
                    @printf("    -> right stop N=%d (BIC degraded past centre+1)\n", n)
                    push!(results, r)
                    break
                end
                if r.bic > best_bic + ccfg.early_stop_dbic
                    # Not competitive: increase streak
                    streak += 1
                else
                    streak = 0
                    best_bic = min(best_bic, r.bic)
                end
            end
            push!(results, r)
        end

        # Sweep left (N decreasing)
        streak = 0
        for n in (center_n-1):-1:n_min_eff
            if streak >= ccfg.early_stop_patience
                @printf("    -> left early stop at N=%d (BIC incr %dx)\n", n, ccfg.early_stop_patience)
                break
            end
            @printf("  Fitting N=%d...\n", n)
            r = _fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg)
            if r.success
                pred = _chain_model_values(xfit, yfit, r.params, n, axisctx, ccfg;
                                           amp_min=r.amp_min, amp_range=r.amp_range)
                _finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, axisctx, ccfg)
                # Left sweep (removing lobes): no overflow risk
                if r.bic > best_bic + ccfg.early_stop_dbic
                    streak += 1
                else
                    streak = 0
                    best_bic = min(best_bic, r.bic)
                end
                if r.bic > best_bic + ccfg.early_stop_dbic
                    streak += 1
                else
                    streak = 0
                    best_bic = min(best_bic, r.bic)
                end
            end
            push!(results, r)
        end
    else
        # ── Legacy linear sweep ──
        for n in ccfg.n_min:n_max_eff
            @printf("  Fitting N=%d...\n", n)
            r = _fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg)
            if r.success
                pred = n == 0 ? fill(r.params[1], length(zfit)) :
                    _chain_model_values(xfit, yfit, r.params, n, axisctx, ccfg;
                                        amp_min=r.amp_min, amp_range=r.amp_range)
                _finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, axisctx, ccfg)
            end
            push!(results, r)
        end
    end

    best = _select_chain_model(results, ccfg)
    return results, best, (xs=xs, ys=ys, zimg=zimg, mask=mask, x=xfit, y=yfit, z=zfit, noise=noise,
                            axisctx=axisctx, axisctx_full=axisctx_full, fit_keep=fit_keep,
                            fit_width_nm=ccfg.fit_width_nm, support_meta=support_meta)
end

function fit_chain_consensus(img::SXMImage, cfg::PatternConfig, ccfg::ChainSweepConfig)
    """Fit chain on Z and Current independently; compare for consensus."""
    # ── Z channel ──
    println("━━━ Consensus: fitting Z channel ━━━")
    results_z, best_z, ctx_z = chain_gaussian_sweep(img, cfg, ccfg)

    # ── Current channel ──
    has_current = any(c -> lowercase(c.name) == "current", img.channels)
    if !has_current
        @warn "No Current channel found in image; consensus skipped"
        return (z=best_z, current=nothing, consensus=false, agreement="no current channel")
    end

    ccfg_c = deepcopy(ccfg)
    ccfg_c.fuse_z_bwd = true  # fuse Current fwd+bwd
    # Narrow sweep: only test near Z's best N
    ccfg_c.n_min = max(2, best_z.n - 1)
    ccfg_c.n_max = best_z.n + 1

    # Use Z's ROI mask and axis for Current fit
    cfg_z = deepcopy(cfg)
    cfg_z.roi_channel = "Z"
    # Get Z mask and data for override
    _, _, _, z_mask, _, _, _, _ = _fused_roi_data(img, cfg_z)
    # Load Current data through Z's mask
    cfg_current = deepcopy(cfg)
    cfg_current.roi_channel = "Current"
    curr_data = _channel_roi_data(img, cfg_current, z_mask, "Current")
    ccfg_c.n_min = max(2, best_z.n - 1)
    ccfg_c.n_max = best_z.n + 1

    println("━━━ Consensus: fitting Current channel (Z mask, Z axis) ━━━")
    results_c, best_c, ctx_c = chain_gaussian_sweep(img, cfg_current, ccfg_c;
        override_data=curr_data, override_axisctx=ctx_z.axisctx_full)

    # ── Compare ──
    n_z = best_z.n; n_c = best_c.n
    s_z = best_z.mean_spacing_nm; s_c = best_c.mean_spacing_nm
    n_agree = n_z == n_c
    s_agree = isapprox(s_z, s_c, rtol=0.20)
    consensus = n_agree  # Primary: N agreement across channels
    if n_agree && s_agree
        agreement = "full — N and spacing match"
    elseif n_agree
        agreement = "N matches (validated) — spacing differs (expected: Z=$(@sprintf("%.3f", s_z)) vs C=$(@sprintf("%.3f", s_c)) nm)"
    else
        agreement = "N mismatch: Z=$n_z vs C=$n_c"
    end

    println("─── Consensus ───")
    println("Z:      N=$n_z  spacing=$(round(s_z, digits=4)) nm  BIC=$(round(best_z.bic, digits=1))")
    if best_c !== nothing
        println("Current:N=$n_c  spacing=$(round(s_c, digits=4)) nm  BIC=$(round(best_c.bic, digits=1))")
    end
    println("Agreement: $agreement")
    println("Consensus: $consensus")

    return (z=(results=results_z, best=best_z, ctx=ctx_z),
            current=(results=results_c, best=best_c, ctx=ctx_c),
            consensus=consensus, agreement=agreement,
            consensus_n=consensus ? n_z : nothing)
end

function fit_chain_batch(filepaths::Vector{String}, cfg::PatternConfig, ccfg::ChainSweepConfig;
                          consensus::Bool=true)
    """Process multiple SXM images; return summary table and per-image results."""
    results = []
    for fp in filepaths
        println("\n" * "="^50)
        println("Processing: $fp")
        println("="^50)
        img = read_sxm(fp)
        local_cfg = deepcopy(cfg)
        local_cfg.filepath = fp
        local_cfg.output_dir = joinpath(cfg.output_dir, replace(basename(fp), r"\.sxm$"i => ""))
        mkpath(local_cfg.output_dir)
        if consensus
            cres = fit_chain_consensus(img, local_cfg, ccfg)
            push!(results, (filepath=fp, z_n=cres.z.best.n, z_bic=cres.z.best.bic,
                           current_n=cres.current.best !== nothing ? cres.current.best.n : nothing,
                           consensus=cres.consensus, agreement=cres.agreement))
        else
            res, best, ctx = chain_gaussian_sweep(img, local_cfg, ccfg)
            _write_chain_sweep(res, best, ctx, local_cfg, ccfg)
            push!(results, (filepath=fp, n=best.n, bic=best.bic, valid=best.valid))
        end
    end
    # Summary table
    println("\n" * "="^70)
    println("BATCH SUMMARY")
    println("="^70)
    if consensus
        for r in results
            println("$(r.filepath): Z_N=$(r.z_n) C_N=$(r.current_n) consensus=$(r.consensus)")
        end
    else
        for r in results
            println("$(r.filepath): N=$(r.n) BIC=$(round(r.bic,digits=1)) valid=$(r.valid)")
        end
    end
    return results
end

function chain_direct_fit(img::SXMImage, cfg::PatternConfig, ccfg::ChainSweepConfig)
    """Fit a single chain model at N = ccfg.n_min (= ccfg.n_max) with optional 1D bootstrap init.
    No sweep.  Returns (result, context).  When skip_global=true, skips NLopt global."""
    ccfg.n_max == ccfg.n_min || error("chain_direct_fit requires n_min == n_max")
    has_bwd = any(c -> lowercase(c.name) == lowercase(cfg.roi_channel) && lowercase(c.direction) == "bwd", img.channels)
    use_fusion = ccfg.fuse_z_bwd && has_bwd
    if use_fusion
        xs, ys, zimg, mask, x, y, z, noise = _fused_roi_data(img, cfg)
    else
        xs, ys, zimg, mask, x, y, z, noise = _robust_roi_data(img, cfg)
    end
    axisctx_full = _weighted_roi_axis(x, y, z)
    xfit, yfit, zfit, axisctx, fit_keep, support_meta = _chain_fit_data(x, y, z, axisctx_full, ccfg)
    # effective sample size (pixels in fit mask ÷ typical spatial correlation factor)
    n_eff = max(10, length(zfit) ÷ 9)
    n = ccfg.n_min
    r = _fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, axisctx, ccfg)
    if r.success
        pred = n == 0 ? fill(r.params[1], length(zfit)) :
            _chain_model_values(xfit, yfit, r.params, n, axisctx, ccfg;
                                amp_min=r.amp_min, amp_range=r.amp_range)
        _finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, axisctx, ccfg)
    end
    return r, (xs=xs, ys=ys, zimg=zimg, mask=mask, x=xfit, y=yfit, z=zfit, noise=noise,
               axisctx=axisctx, axisctx_full=axisctx_full, fit_keep=fit_keep,
               fit_width_nm=ccfg.fit_width_nm, support_meta=support_meta)
end

# -----------------------------------------------------------------------------
# Constrained Gaussian refinement
# -----------------------------------------------------------------------------

_sigmoid(t) = 1 / (1 + exp(-clamp(t, -60.0, 60.0)))
_logit(u) = log(clamp(u, 1e-6, 1 - 1e-6) / (1 - clamp(u, 1e-6, 1 - 1e-6)))

function _pack_initial(features::Vector{MolecularFeature}, img::SXMImage, cfg::PatternConfig)
    p = Float64[0.0]
    sx0 = something(cfg.initial_sigma_nm, clamp(min(img.range_nm...) / 25, cfg.min_sigma_nm, cfg.max_sigma_nm))
    sy0 = sx0
    sigma_span = max(cfg.max_sigma_nm - cfg.min_sigma_nm, EPS)
    for f in features
        push!(p, f.amplitude)
        push!(p, _logit(f.x_nm / img.range_nm[1]))
        push!(p, _logit(f.y_nm / img.range_nm[2]))
        push!(p, _logit((sx0 - cfg.min_sigma_nm) / sigma_span))
        push!(p, _logit((sy0 - cfg.min_sigma_nm) / sigma_span))
    end
    return p
end

function _unpack_features(p::AbstractVector, img::SXMImage, cfg::PatternConfig)
    n = (length(p) - 1) ÷ 5
    features = MolecularFeature[]
    for i in 0:(n-1)
        j = 2 + 5i
        A = p[j]
        x0 = img.range_nm[1] * _sigmoid(p[j+1])
        y0 = img.range_nm[2] * _sigmoid(p[j+2])
        sigma_span = max(cfg.max_sigma_nm - cfg.min_sigma_nm, EPS)
        sx = cfg.min_sigma_nm + sigma_span * _sigmoid(p[j+3])
        sy = cfg.min_sigma_nm + sigma_span * _sigmoid(p[j+4])
        push!(features, MolecularFeature(A, x0, y0, sx, sy, abs(A)))
    end
    return features
end

function _model_values(x::AbstractVector, y::AbstractVector, p::AbstractVector, img::SXMImage, cfg::PatternConfig)
    isempty(p) && return zeros(length(x))
    z = fill(p[1], length(x))
    for f in _unpack_features(p, img, cfg)
        @. z += f.amplitude * exp(-0.5 * (((x - f.x_nm) / f.sigma_x_nm)^2 + ((y - f.y_nm) / f.sigma_y_nm)^2))
    end
    return z
end

function fit_molecular_pattern(img::SXMImage, cfg::PatternConfig)
    ch = get_channel(img, cfg.channel; direction=cfg.direction)
    xs, ys, raw, z, z_smooth, unit, noise = preprocess_channel(img, ch, cfg)
    roi_mask = falses(0, 0)
    axis_peaks = MolecularFeature[]
    roi_length = NaN
    estimated_n = 0
    estimated_range = (0, 0)
    if cfg.fusion
        _fxs, _fys, evidence, candidates, chains, accepted, roi_mask, axis_peaks, roi_length, estimated_n, estimated_range = detect_molecular_chains(img, cfg)
    else
        candidates = detect_blobs(z_smooth, xs, ys, cfg, noise)
        chains = extract_chains(candidates, cfg)
        accepted = MolecularFeature[]
        for chain in chains
            append!(accepted, chain.features)
        end
        evidence = zeros(0, 0)
    end
    if length(candidates) < cfg.min_features
        return PatternFitResult(raw_features=candidates, chains=chains, evidence_map=evidence,
                                roi_mask=roi_mask, axis_peaks=axis_peaks,
                                roi_length_nm=roi_length, estimated_repeats=estimated_n,
                                estimated_repeat_range=estimated_range,
                                warnings=["Not enough candidate features detected."])
    end
    if isempty(accepted)
        return PatternFitResult(raw_features=candidates, chains=chains, evidence_map=evidence,
                                roi_mask=roi_mask, axis_peaks=axis_peaks,
                                roi_length_nm=roi_length, estimated_repeats=estimated_n,
                                estimated_repeat_range=estimated_range,
                                warnings=["No valid molecular chain found; candidates are exported for threshold tuning."])
    end

    nfit = min(length(accepted), cfg.max_fit_features)
    initial = accepted[1:nfit]
    xflat = repeat(xs, inner=length(ys))
    yflat = repeat(ys, outer=length(xs))
    zflat = vec(z) # correct alignment: z[y,x]

    result = PatternFitResult(raw_features=candidates, chains=chains, evidence_map=evidence,
                              roi_mask=roi_mask, axis_peaks=axis_peaks,
                              roi_length_nm=roi_length, estimated_repeats=estimated_n,
                              estimated_repeat_range=estimated_range)
    p0 = _pack_initial(initial, img, cfg)
    if cfg.no_fit
        result.params_unconstrained = p0
        result.features = _unpack_features(p0, img, cfg)
        result.success = true
    else
        xydata = vcat(reshape(xflat, 1, :), reshape(yflat, 1, :))
        model = (xy, p) -> _model_values(view(xy, 1, :), view(xy, 2, :), p, img, cfg)
        try
            fit = curve_fit(model, xydata, zflat, p0; maxIter=cfg.max_iter, autodiff=:finite)
            result.params_unconstrained = fit.param
            result.features = _unpack_features(fit.param, img, cfg)
            result.success = true
        catch e
            result.params_unconstrained = p0
            result.features = _unpack_features(p0, img, cfg)
            result.success = false
            push!(result.warnings, "Constrained refinement failed; using detected blobs. Error: $e")
        end
    end

    pred = _model_values(xflat, yflat, result.params_unconstrained, img, cfg)
    result.rss = sum(abs2, zflat .- pred)
    tss = sum(abs2, zflat .- mean(zflat))
    result.r_squared = 1 - result.rss / max(tss, EPS)
    k = 1 + 5length(result.features)
    n_eff = max(1, length(zflat) ÷ max(1, cfg.smooth_radius_px + cfg.min_distance_px ÷ max(1, cfg.stride)))
    result.bic = n_eff * log(max(result.rss, EPS) / length(zflat)) + k * log(n_eff)
    return result
end

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

function _write_parameters(result::PatternFitResult, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "molecular_features.tsv")
    open(out, "w") do io
        println(io, "type\tfeature\tamplitude\tx_nm\ty_nm\tsigma_x_nm\tsigma_y_nm\tscore")
        for (i, f) in enumerate(result.raw_features)
            println(io, "detected\t$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(f.sigma_x_nm)\t$(f.sigma_y_nm)\t$(f.score)")
        end
        for (i, f) in enumerate(result.features)
            println(io, "fitted\t$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(f.sigma_x_nm)\t$(f.sigma_y_nm)\t$(f.score)")
        end
    end
    result.parameter_file = out
    estimate_out = joinpath(cfg.output_dir, "chain_estimate.tsv")
    open(estimate_out, "w") do io
        println(io, "roi_length_nm\trepeat_spacing_nm\tend_extension_nm\testimated_repeats\testimated_range_min\testimated_range_max\tvisible_chain_features")
        println(io, "$(result.roi_length_nm)\t$(cfg.repeat_spacing_nm)\t$(cfg.end_extension_nm)\t$(result.estimated_repeats)\t$(result.estimated_repeat_range[1])\t$(result.estimated_repeat_range[2])\t$(isempty(result.chains) ? 0 : length(result.chains[1].features))")
    end
    chain_out = joinpath(cfg.output_dir, "molecular_chains.tsv")
    open(chain_out, "w") do io
        println(io, "chain_id\tindex_in_chain\tx_nm\ty_nm\tscore\tmean_spacing_nm\tspacing_cv\tmax_turn_angle_deg\tchain_score")
        for chain in result.chains
            for (i, f) in enumerate(chain.features)
                println(io, "$(chain.id)\t$(i)\t$(f.x_nm)\t$(f.y_nm)\t$(f.score)\t$(chain.mean_spacing_nm)\t$(chain.spacing_cv)\t$(chain.max_turn_angle_deg)\t$(chain.score)")
            end
        end
    end
    return out
end

_safe_name(s::String) = replace(lowercase(strip(s)), r"[^a-z0-9]+" => "_")

function _write_all_images_summary(rows, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "all_images_summary.tsv")
    open(out, "w") do io
        println(io, "channel\tdirection\tsuccess\tn_detected\tn_fitted\tr_squared\trss\toutput_dir")
        for row in rows
            println(io, "$(row.channel)\t$(row.direction)\t$(row.success)\t$(row.n_detected)\t$(row.n_fitted)\t$(row.r_squared)\t$(row.rss)\t$(row.output_dir)")
        end
    end
    return out
end

function _plot_result(img::SXMImage, cfg::PatternConfig, result::PatternFitResult)
    ch = get_channel(img, cfg.channel; direction=cfg.direction)
    xs, ys, raw, z, z_smooth, unit, _noise = preprocess_channel(img, ch, cfg)
    xflat = repeat(xs, inner=length(ys))
    yflat = repeat(ys, outer=length(xs))
    pred = reshape(_model_values(xflat, yflat, result.params_unconstrained, img, cfg), length(ys), length(xs))

    p1 = heatmap(xs, ys, raw; aspect_ratio=:equal, title="Raw $(ch.name) $(ch.direction) [$unit]", xlabel="x (nm)", ylabel="y (nm)")
    evidence = result.evidence_map
    if cfg.fusion && !isempty(evidence)
        p2 = heatmap(xs, ys, evidence; aspect_ratio=:equal, title="Fused multi-view evidence + chains", xlabel="x (nm)", ylabel="y (nm)")
    else
        p2 = heatmap(xs, ys, z; aspect_ratio=:equal, title="Preprocessed + chains", xlabel="x (nm)", ylabel="y (nm)")
    end
    for f in result.raw_features
        scatter!(p2, [f.x_nm], [f.y_nm]; marker=:circle, markersize=3, label=false, color=:gray)
    end
    chain_colors = [:red, :cyan, :yellow, :lime, :magenta, :orange, :white, :blue]
    for (ci, chain) in enumerate(result.chains)
        color = chain_colors[mod1(ci, length(chain_colors))]
        xs_chain = [f.x_nm for f in chain.features]
        ys_chain = [f.y_nm for f in chain.features]
        plot!(p2, xs_chain, ys_chain; linewidth=2.5, color=color, label=false)
        scatter!(p2, xs_chain, ys_chain; marker=:circle, markersize=5, color=color, label=false)
    end
    p3 = heatmap(xs, ys, pred; aspect_ratio=:equal, title="Chain-constrained Gaussian refinement ($(length(result.features)) fitted)", xlabel="x (nm)", ylabel="y (nm)")
    for f in result.features
        scatter!(p3, [f.x_nm], [f.y_nm]; marker=:cross, markersize=5, label=false, color=:red)
    end
    p4 = heatmap(xs, ys, z .- pred; aspect_ratio=:equal, title="Residuals", xlabel="x (nm)", ylabel="y (nm)")
    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(1200, 950))
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "molecular_pattern_overview.png")
    savefig(fig, out)
    result.plot_file = out
    return out
end

function _plot_all_images(img::SXMImage, cfg::PatternConfig)
    mkpath(cfg.output_dir)
    panels = Any[]
    for ch in img.channels
        local_cfg = deepcopy(cfg)
        local_cfg.channel = ch.name
        local_cfg.direction = ch.direction
        xs, ys, raw, z, z_smooth, unit, noise = preprocess_channel(img, ch, local_cfg)
        candidates = detect_blobs(z_smooth, xs, ys, local_cfg, noise)
        p_raw = heatmap(xs, ys, raw; aspect_ratio=:equal,
                        title="Raw $(ch.name) $(ch.direction) [$unit]",
                        xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
        p_pre = heatmap(xs, ys, z; aspect_ratio=:equal,
                        title="Preprocessed + candidates ($(length(candidates)))",
                        xlabel="x (nm)", ylabel="y (nm)", colorbar=false)
        for f in candidates
            scatter!(p_pre, [f.x_nm], [f.y_nm]; marker=:circle, markersize=3,
                     label=false, color=:white)
        end
        push!(panels, p_raw, p_pre)
    end
    fig = plot(panels...; layout=(length(img.channels), 2), size=(1100, 280length(img.channels)))
    out = joinpath(cfg.output_dir, "all_images_overview.png")
    savefig(fig, out)
    return out
end

function _run_all_images(img::SXMImage, cfg::PatternConfig)
    rows = NamedTuple[]
    for ch in img.channels
        local_cfg = deepcopy(cfg)
        local_cfg.channel = ch.name
        local_cfg.direction = ch.direction
        local_cfg.output_dir = joinpath(cfg.output_dir, "$(_safe_name(ch.name))_$(ch.direction)")
        result = fit_molecular_pattern(img, local_cfg)
        _write_parameters(result, local_cfg)
        !cfg.no_plot && _plot_result(img, local_cfg, result)
        push!(rows, (channel=ch.name, direction=ch.direction, success=result.success,
                     n_detected=length(result.raw_features), n_fitted=length(result.features), r_squared=result.r_squared,
                     rss=result.rss, output_dir=local_cfg.output_dir))
        println("$(ch.name) $(ch.direction): detected=$(length(result.raw_features)) fitted=$(length(result.features)) R²=$(@sprintf("%.5f", result.r_squared)) dir=$(local_cfg.output_dir)")
    end
    summary = _write_all_images_summary(rows, cfg)
    overview = cfg.no_plot ? nothing : _plot_all_images(img, cfg)
    println("all_images_summary: $summary")
    overview !== nothing && println("all_images_overview: $overview")
end

function _parse_sweep(spec::String)
    parts = split(spec, ':')
    length(parts) == 3 || error("threshold sweep must be start:step:stop, e.g. 2.0:0.25:4.0")
    a, step, b = parse.(Float64, parts)
    step > 0 || error("threshold sweep step must be > 0")
    return collect(a:step:b)
end

function _run_threshold_sweep(img::SXMImage, cfg::PatternConfig, spec::String)
    mkpath(cfg.output_dir)
    out = joinpath(cfg.output_dir, "threshold_sweep.tsv")
    open(out, "w") do io
        println(io, "threshold_sigma\tn_candidates\tn_chains\tn_chain_features\tbest_chain_length\tbest_chain_score\troi_length_nm\testimated_repeats\testimated_range")
        for thr in _parse_sweep(spec)
            local_cfg = deepcopy(cfg)
            local_cfg.threshold_sigma = thr
            _xs, _ys, _evidence, candidates, chains, accepted, _roi, _axis, roi_len, est_n, est_range = detect_molecular_chains(img, local_cfg)
            best_len = isempty(chains) ? 0 : maximum(length(c.features) for c in chains)
            best_score = isempty(chains) ? -Inf : maximum(c.score for c in chains)
            println(io, "$(thr)\t$(length(candidates))\t$(length(chains))\t$(length(accepted))\t$(best_len)\t$(best_score)\t$(roi_len)\t$(est_n)\t$(est_range[1])-$(est_range[2])")
            println("threshold=$(thr): candidates=$(length(candidates)) chains=$(length(chains)) best_length=$(best_len) estimated_repeats=$(est_n)")
        end
    end
    println("threshold_sweep: $out")
    return out
end

function _write_chain_sweep(results::Vector{ChainModelResult}, best::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    mkpath(cfg.output_dir)
    summary = joinpath(cfg.output_dir, "chain_model_selection.tsv")
    open(summary, "w") do io
        println(io, "N\tn_params\tsuccess\tvalid\ttrain_nll\tcv_nll_mean\tcv_nll_std\tbic\taicc\tbic_per_pixel\taicc_per_pixel\trss\tchi2_reduced\tmad\tresidual_peak_snr\tmean_spacing_nm\tspacing_cv\tmax_lateral_nm\tsigma_parallel_nm\tsigma_perp_nm\toverlap\tendpoint_overrun_nm\tbound_like\treason\twarnings\tselected")
        for r in results
            npx = max(1, length(ctx.z))
            println(io, "$(r.n)\t$(_chain_nparams(r.n, ccfg))\t$(r.success)\t$(r.valid)\t$(r.train_nll)\t$(r.cv_nll_mean)\t$(r.cv_nll_std)\t$(r.bic)\t$(r.aicc)\t$(r.bic/npx)\t$(r.aicc/npx)\t$(r.rss)\t$(r.chi2_reduced)\t$(r.mad)\t$(r.residual_peak_snr)\t$(r.mean_spacing_nm)\t$(r.spacing_cv)\t$(r.max_lateral_nm)\t$(r.sigma_parallel_nm)\t$(r.sigma_perp_nm)\t$(r.overlap)\t$(r.endpoint_overrun_nm)\t$(r.bound_like)\t$(r.reason)\t$(r.valid ? "" : r.reason)\t$(r === best)")
        end
    end
    params = joinpath(cfg.output_dir, "chain_selected_lobes.tsv")
    open(params, "w") do io
        println(io, "lobe\tamplitude\tx_nm\ty_nm\tt_nm\tu_nm\tsigma_parallel_nm\tsigma_perp_nm\tspacing_prev_nm")
        if best.n > 0 && best.success
            _b, feats, ts, us, spars, sperps = _decode_chain(best.params, best.n, ctx.axisctx, ccfg;
                                                              amp_min=best.amp_min, amp_range=best.amp_range)
            for (i, f) in enumerate(feats)
                spacing_prev = i == 1 ? NaN : ts[i] - ts[i-1]
                println(io, "$(i)\t$(f.amplitude)\t$(f.x_nm)\t$(f.y_nm)\t$(ts[i])\t$(us[i])\t$(spars[i])\t$(sperps[i])\t$(spacing_prev)")
            end
        end
    end
    axisfile = joinpath(cfg.output_dir, "chain_axis.tsv")
    open(axisfile, "w") do io
        println(io, "origin_x_nm\torigin_y_nm\taxis_x\taxis_y\tperp_x\tperp_y\ttmin_nm\ttmax_nm")
        println(io, "$(ctx.axisctx.origin[1])\t$(ctx.axisctx.origin[2])\t$(ctx.axisctx.axis[1])\t$(ctx.axisctx.axis[2])\t$(ctx.axisctx.perp[1])\t$(ctx.axisctx.perp[2])\t$(ctx.axisctx.tmin)\t$(ctx.axisctx.tmax)")
    end
    supportfile = joinpath(cfg.output_dir, "chain_support_metadata.tsv")
    open(supportfile, "w") do io
        println(io, "key\tvalue")
        println(io, "pipeline\tauto ROI -> weighted PCA axis -> tube axial profile -> robust support -> adaptive n_max -> deterministic global+local fit -> BIC selection")
        println(io, "n_min_config\t$(ccfg.n_min)")
        println(io, "n_max_config\t$(ccfg.n_max)")
        println(io, "spacing_min_nm\t$(ccfg.spacing_min_nm)")
        println(io, "spacing_min_effective_nm\t$(_effective_spacing_min_nm(ccfg))")
        println(io, "spacing_min_effective_source\tmax(spacing_min_nm, sqrt(-2log(max_overlap))*sigma_max)")
        println(io, "spacing_max_nm\t$(ccfg.spacing_max_nm)")
        println(io, "max_overlap\t$(ccfg.max_overlap)")
        println(io, "sigma_parallel_min_nm\t$(ccfg.sigma_parallel_min_nm)")
        println(io, "sigma_parallel_max_nm\t$(ccfg.sigma_parallel_max_nm)")
        println(io, "sigma_perp_min_nm\t$(ccfg.sigma_perp_min_nm)")
        println(io, "sigma_perp_max_nm\t$(ccfg.sigma_perp_max_nm)")
        for k in keys(ctx.support_meta)
            println(io, "$(k)\t$(getfield(ctx.support_meta, k))")
        end
    end
    profile = joinpath(cfg.output_dir, "chain_axis_profile.tsv")
    pred = best.success ? (best.n == 0 ? fill(best.params[1], length(ctx.z)) :
        _chain_model_values(ctx.x, ctx.y, best.params, best.n, ctx.axisctx, ccfg;
                            amp_min=best.amp_min, amp_range=best.amp_range)) : fill(NaN, length(ctx.z))
    t = (ctx.x .- ctx.axisctx.origin[1]) .* ctx.axisctx.axis[1] .+ (ctx.y .- ctx.axisctx.origin[2]) .* ctx.axisctx.axis[2]
    nb = max(20, min(200, Int(ceil((maximum(t) - minimum(t)) / max(median(diff(ctx.xs)), EPS)))))
    sdata = zeros(nb); smodel = zeros(nb); sresid = zeros(nb); counts = zeros(Int, nb)
    for i in eachindex(t)
        b = clamp(Int(floor((t[i] - minimum(t)) / max(maximum(t) - minimum(t), EPS) * (nb - 1))) + 1, 1, nb)
        sdata[b] += ctx.z[i]; smodel[b] += pred[i]; sresid[b] += ctx.z[i] - pred[i]; counts[b] += 1
    end
    open(profile, "w") do io
        println(io, "t_nm\tdata_mean\tmodel_mean\tresidual_mean\tcount")
        for b in 1:nb
            tb = minimum(t) + (b - 0.5) / nb * (maximum(t) - minimum(t))
            c = max(counts[b], 1)
            println(io, "$(tb)\t$(sdata[b]/c)\t$(smodel[b]/c)\t$(sresid[b]/c)\t$(counts[b])")
        end
    end
    return summary, params, axisfile, profile, supportfile
end

function _plot_chain_sweep(results::Vector{ChainModelResult}, best::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    xs, ys, zimg, mask = ctx.xs, ctx.ys, ctx.zimg, ctx.mask
    pred_img = zeros(size(zimg))
    if best.success
        for iy in eachindex(ys), ix in eachindex(xs)
            pred_img[iy, ix] = best.n == 0 ? best.params[1] :
                _chain_model_values([xs[ix]], [ys[iy]], best.params, best.n, ctx.axisctx, ccfg;
                                    amp_min=best.amp_min, amp_range=best.amp_range)[1]
        end
    end
    p1 = heatmap(xs, ys, zimg; aspect_ratio=:equal, title="Chain ROI data", xlabel="x (nm)", ylabel="y (nm)")
    contour!(p1, xs, ys, Float64.(mask); levels=[0.5], color=:white, linewidth=2)
    ox, oy = ctx.axisctx.origin; ax, ay = ctx.axisctx.axis
    plot!(p1, [ox + ctx.axisctx.tmin*ax, ox + ctx.axisctx.tmax*ax], [oy + ctx.axisctx.tmin*ay, oy + ctx.axisctx.tmax*ay]; color=:yellow, linewidth=2, label="axis")
    if best.n > 0 && best.success
        _b, feats, _ts, _us, _spars, _sperps = _decode_chain(best.params, best.n, ctx.axisctx, ccfg;
                                                              amp_min=best.amp_min, amp_range=best.amp_range)
        scatter!(p1, [f.x_nm for f in feats], [f.y_nm for f in feats]; color=:cyan, marker=:cross, label="lobes")
    end
    p2 = heatmap(xs, ys, pred_img; aspect_ratio=:equal, title="Ordered chain model N=$(best.n)", xlabel="x (nm)", ylabel="y (nm)")
    p3 = heatmap(xs, ys, (zimg .- pred_img) .* mask ./ max(ctx.noise, EPS); aspect_ratio=:equal, title="Residuals / noise", xlabel="x (nm)", ylabel="y (nm)")
    Ns = [r.n for r in results]
    p4 = plot(Ns, [r.cv_nll_mean for r in results]; marker=:circle, label="CV NLL", xlabel="N", ylabel="score", title="Chain model selection")
    plot!(p4, Ns, [r.bic / max(1, length(ctx.z)) for r in results]; marker=:diamond, label="BIC / n")
    plot!(p4, Ns, [r.residual_peak_snr for r in results]; marker=:utriangle, label="resid SNR")
    vline!(p4, [best.n]; label="selected", color=:red)
    fig = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 950))
    out = joinpath(cfg.output_dir, "chain_model_selection.png")
    savefig(fig, out)
    return out
end

function _plot_chain_model_grid(r::ChainModelResult, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig, outdir::String)
    r.success || return nothing
    xs, ys, zimg, mask = ctx.xs, ctx.ys, ctx.zimg, ctx.mask
    pred_img = zeros(size(zimg))
    for iy in eachindex(ys), ix in eachindex(xs)
        pred_img[iy, ix] = r.n == 0 ? r.params[1] :
            _chain_model_values([xs[ix]], [ys[iy]], r.params, r.n, ctx.axisctx, ccfg;
                                amp_min=r.amp_min, amp_range=r.amp_range)[1]
    end
    resid = (zimg .- pred_img) .* mask ./ max(ctx.noise, EPS)
    p1 = heatmap(xs, ys, zimg; aspect_ratio=:equal, title="Data ROI", xlabel="x (nm)", ylabel="y (nm)")
    contour!(p1, xs, ys, Float64.(mask); levels=[0.5], color=:white, linewidth=2)
    ox, oy = ctx.axisctx.origin; ax, ay = ctx.axisctx.axis
    plot!(p1, [ox + ctx.axisctx.tmin*ax, ox + ctx.axisctx.tmax*ax], [oy + ctx.axisctx.tmin*ay, oy + ctx.axisctx.tmax*ay]; color=:yellow, linewidth=2, label=false)
    if r.n > 0
        _b, feats, _ts, _us, _spars, _sperps = _decode_chain(r.params, r.n, ctx.axisctx, ccfg;
                                                              amp_min=r.amp_min, amp_range=r.amp_range)
        scatter!(p1, [f.x_nm for f in feats], [f.y_nm for f in feats]; color=:cyan, marker=:cross, label=false)
    end
    p2 = heatmap(xs, ys, pred_img; aspect_ratio=:equal, title="Chain model N=$(r.n)", xlabel="x (nm)", ylabel="y (nm)")
    p3 = heatmap(xs, ys, resid; aspect_ratio=:equal, title="Residual / noise", xlabel="x (nm)", ylabel="y (nm)")
    p4 = plot(title="scores", legend=false, framestyle=:box)
    annotate!(p4, 0.05, 0.82, text("valid=$(r.valid)", :left, 10))
    annotate!(p4, 0.05, 0.68, text(@sprintf("CV %.4g ± %.3g", r.cv_nll_mean, r.cv_nll_std), :left, 10))
    annotate!(p4, 0.05, 0.54, text(@sprintf("BIC %.4g", r.bic), :left, 10))
    annotate!(p4, 0.05, 0.40, text(@sprintf("resid SNR %.3g", r.residual_peak_snr), :left, 10))
    annotate!(p4, 0.05, 0.26, text(@sprintf("spacing %.3g CV %.3g", r.mean_spacing_nm, r.spacing_cv), :left, 10))
    annotate!(p4, 0.05, 0.12, text(@sprintf("σ∥ %.3g σ⊥ %.3g", r.sigma_parallel_nm, r.sigma_perp_nm), :left, 10))
    xlims!(p4, 0, 1); ylims!(p4, 0, 1)
    fig = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 950))
    out = joinpath(outdir, @sprintf("chain_model_N%02d.png", r.n))
    savefig(fig, out)
    return out
end

function _plot_chain_models_by_n(results::Vector{ChainModelResult}, ctx, cfg::PatternConfig, ccfg::ChainSweepConfig)
    outdir = joinpath(cfg.output_dir, "models_by_N")
    mkpath(outdir)
    paths = String[]
    for r in results
        path = _plot_chain_model_grid(r, ctx, cfg, ccfg, outdir)
        path !== nothing && push!(paths, path)
    end
    return paths
end

function _run_chain_sweep(img::SXMImage, cfg::PatternConfig, ccfg::ChainSweepConfig)
    results, best, ctx = chain_gaussian_sweep(img, cfg, ccfg)
    summary, params, axisfile, profile, supportfile = _write_chain_sweep(results, best, ctx, cfg, ccfg)
    plot_path = cfg.no_plot ? nothing : _plot_chain_sweep(results, best, ctx, cfg, ccfg)
    model_plots = cfg.no_plot ? String[] : _plot_chain_models_by_n(results, ctx, cfg, ccfg)
    println("chain_selected_N: $(best.n)")
    println("chain_selected_valid: $(best.valid)")
    println("chain_selected_cv_nll: $(best.cv_nll_mean) ± $(best.cv_nll_std)")
    println("chain_selected_bic: $(best.bic)")
    println("chain_residual_peak_snr: $(best.residual_peak_snr)")
    println("chain_mean_spacing_nm: $(best.mean_spacing_nm)")
    println("chain_summary: $summary")
    println("chain_lobes: $params")
    println("chain_axis: $axisfile")
    println("chain_axis_profile: $profile")
    println("chain_support_metadata: $supportfile")
    plot_path !== nothing && println("chain_plot: $plot_path")
    !isempty(model_plots) && println("chain_models_by_N_dir: $(joinpath(cfg.output_dir, "models_by_N"))")
end

function _print_info(img::SXMImage)
    println("File: $(img.filepath)")
    println("Pixels: $(img.width) × $(img.height)")
    println("Range: $(@sprintf("%.3f", img.range_nm[1])) nm × $(@sprintf("%.3f", img.range_nm[2])) nm")
    println("Channels:")
    for c in img.channels
        println("  - $(c.name) direction=$(c.direction) unit=$(c.unit) min=$(minimum(c.data)) max=$(maximum(c.data))")
    end
end

function _arg_settings()
    s = ArgParseSettings(description="STM molecular pattern recognition by preprocessing, blob detection and constrained Gaussian refinement.")
    @add_arg_table! s begin
        "filepath"
            help = "Fichier Nanonis .sxm"
            required = true
        "--info"
            action = :store_true
        "--channel"
            default = "Z"
        "--direction"
            default = "fwd"
        "--output-dir", "-o"
            default = "results"
        "--stride"
            arg_type = Int
            default = 2
        "--flatten"
            help = "none | plane | rows | plane+rows"
            default = "plane+rows"
        "--smooth-radius-px"
            arg_type = Int
            default = 1
        "--contrast"
            help = "auto | bright | dark"
            default = "auto"
        "--max-features"
            help = "Nombre maximal de candidats détectés/exportés"
            arg_type = Int
            default = 50
        "--max-fit-features"
            help = "Nombre maximal de candidats raffinés par fit non linéaire"
            arg_type = Int
            default = 8
        "--threshold-sigma"
            arg_type = Float64
            default = 2.5
        "--min-distance-px"
            arg_type = Int
            default = 10
        "--ignore-border-px"
            arg_type = Int
            default = 8
        "--chain-min-length"
            arg_type = Int
            default = 3
        "--chain-min-spacing-nm"
            arg_type = Float64
            default = 0.35
        "--chain-max-spacing-nm"
            arg_type = Float64
            default = 1.6
        "--chain-max-angle-deg"
            arg_type = Float64
            default = 45.0
        "--chain-max-spacing-cv"
            arg_type = Float64
            default = 0.45
        "--min-chain-score"
            arg_type = Float64
            default = 5.0
        "--max-chains"
            arg_type = Int
            default = 8
        "--max-path-branches"
            arg_type = Int
            default = 6
        "--repeat-spacing-nm"
            help = "Pas moléculaire attendu pour estimer le nombre d'unités (chitosane ≈ 0.52 nm)"
            arg_type = Float64
            default = 0.52
        "--end-extension-nm"
            help = "Correction d'enveloppe STM aux deux extrémités avant longueur/pas"
            arg_type = Float64
            default = 0.35
        "--no-fusion"
            help = "Ne pas fusionner les vues; détecter seulement le canal/direction choisi"
            action = :store_true
        "--fusion-channels"
            help = "Canaux utilisés pour la carte d'évidence fusionnée: Z, Current, Z,Current, ou all"
            default = "Z"
        "--min-sigma-nm"
            arg_type = Float64
            default = 0.05
        "--max-sigma-nm"
            arg_type = Float64
            default = 1.5
        "--initial-sigma-nm"
            arg_type = Float64
            default = nothing
        "--no-fit"
            action = :store_true
        "--no-plot"
            action = :store_true
        "--all-images"
            help = "Traiter tous les couples canal/direction du fichier SXM"
            action = :store_true
        "--inspect-images"
            help = "Créer seulement une planche raw/preprocessed pour toutes les images"
            action = :store_true
        "--threshold-sweep"
            help = "Balayage du seuil sans fit, format start:step:stop, ex. 2.0:0.25:4.0"
            default = ""
        "--selection-criterion"
            help = "Model selection criterion: bic | mdl | cv (default: bic)"
            default = "bic"
        "--chain-sweep"
            help = "Sélection N par chaîne 2D ordonnée"
            action = :store_true
        "--consensus"
            help = "Multi-channel consensus: fit on Z then on Current, compare N/spacing"
            action = :store_true
        "--batch"
            help = "Consensus batch processing on all .sxm files in a directory"
            arg_type = String
            default = ""
        "--chain-mode"
            help = "Chain fitting mode: sweep (default) | 1d-bootstrapped | direct"
            arg_type = String
            default = "sweep"
        "--chain-skip-global"
            help = "Skip NLopt global optimization when 1D bootstrap init is available (faster)"
            action = :store_true
        "--chain-init-centers-t"
            help = "Comma-separated list of t-axis centers (nm) from 1D fit, for 1d-bootstrapped mode"
            arg_type = String
            default = ""
        "--chain-init-amplitudes"
            help = "Comma-separated list of amplitudes from 1D fit, for 1d-bootstrapped mode"
            arg_type = String
            default = ""
        "--chain-init-sigma-parallel"
            help = "Mean sigma from 1D fit, for sigma_parallel initialization (nm)"
            arg_type = Float64
            default = NaN
        "--chain-init-sigma-perp"
            help = "Mean sigma from 1D fit, for sigma_perp initialization (nm)"
            arg_type = Float64
            default = NaN
        "--multistart"
            arg_type = Int
            default = 20
        "--cv-folds"
            arg_type = Int
            default = 5
        "--chain-n-min"
            arg_type = Int
            default = 2
        "--chain-n-max"
            arg_type = Int
            default = 14
        "--chain-spacing-min-nm"
            arg_type = Float64
            default = 0.35
        "--chain-spacing-max-nm"
            arg_type = Float64
            default = 0.75
        "--chain-lateral-max-nm"
            arg_type = Float64
            default = 0.35
        "--chain-fit-width-nm"
            help = "Demi-largeur du tube fitté autour de l'axe moléculaire; rend le fit comparable à une slice 1D"
            arg_type = Float64
            default = 0.45
        "--chain-support-baseline-quantile"
            help = "Quantile robuste de baseline du profil axial pour support auto"
            arg_type = Float64
            default = 0.10
        "--chain-support-threshold-fraction"
            help = "Seuil relatif du profil axial pour garder seulement le support longitudinal actif"
            arg_type = Float64
            default = 0.25
        "--chain-support-noise-k"
            help = "Seuil bruit du support auto: baseline + k*MAD; combiné au seuil contraste par max()"
            arg_type = Float64
            default = 2.5
        "--chain-support-padding-nm"
            arg_type = Float64
            default = 0.20
        "--chain-support-min-length-nm"
            help = "Longueur minimale d'un composant actif longitudinal"
            arg_type = Float64
            default = 1.0
        "--chain-t-min-nm"
            help = "Absolute minimum t coordinate along chain axis; use slide_metadata support_start_t_nm to match a 1D slice"
            arg_type = Float64
            default = NaN
        "--chain-t-max-nm"
            help = "Absolute maximum t coordinate along chain axis; use slide_metadata support_end_t_nm to match a 1D slice"
            arg_type = Float64
            default = NaN
        "--chain-sigma-parallel-min-nm"
            arg_type = Float64
            default = 0.12
        "--chain-sigma-parallel-max-nm"
            arg_type = Float64
            default = 0.45
        "--chain-sigma-perp-min-nm"
            arg_type = Float64
            default = 0.10
        "--chain-sigma-perp-max-nm"
            arg_type = Float64
            default = 0.55
        "--fuse-z-bwd"
            help = "Fuse Z fwd + Z bwd for SNR boost (default: true)"
            action = :store_true
            default = true
        "--no-fuse-z-bwd"
            help = "Disable Z fwd+bwd fusion, use Z fwd only"
            action = :store_true
            default = false
        "--global-maxtime"
            help = "NLopt global optimization max time per model (seconds)"
            arg_type = Float64
            default = 15.0
        "--global-maxiter"
            help = "NLopt global optimization max evaluations per model"
            arg_type = Int
            default = 5000
        "--global-tol"
            help = "NLopt relative tolerance"
            arg_type = Float64
            default = 1e-5
    end
    return s
end

function main_cli()
    args = parse_args(_arg_settings())
    img = read_sxm(args["filepath"])
    if args["info"]
        _print_info(img)
        return
    end
    cfg = PatternConfig(
        filepath=args["filepath"], channel=args["channel"], direction=args["direction"],
        output_dir=args["output-dir"], stride=args["stride"], flatten=args["flatten"],
        smooth_radius_px=args["smooth-radius-px"], contrast=args["contrast"],
        max_features=args["max-features"], max_fit_features=args["max-fit-features"],
        threshold_sigma=args["threshold-sigma"],
        min_distance_px=args["min-distance-px"], ignore_border_px=args["ignore-border-px"],
        fusion=!args["no-fusion"], chain_min_length=args["chain-min-length"],
        fusion_channels=args["fusion-channels"],
        chain_min_spacing_nm=args["chain-min-spacing-nm"],
        chain_max_spacing_nm=args["chain-max-spacing-nm"],
        chain_max_angle_deg=args["chain-max-angle-deg"],
        chain_max_spacing_cv=args["chain-max-spacing-cv"], min_chain_score=args["min-chain-score"],
        max_chains=args["max-chains"], max_path_branches=args["max-path-branches"],
        repeat_spacing_nm=args["repeat-spacing-nm"], end_extension_nm=args["end-extension-nm"],
        min_sigma_nm=args["min-sigma-nm"], max_sigma_nm=args["max-sigma-nm"],
        initial_sigma_nm=args["initial-sigma-nm"], no_fit=args["no-fit"], no_plot=args["no-plot"])
    cfg.all_images = args["all-images"]
    cfg.inspect_images = args["inspect-images"]

    if !isempty(args["threshold-sweep"])
        _run_threshold_sweep(img, cfg, args["threshold-sweep"])
        return
    end

    if args["chain-sweep"]
        chain_mode = lowercase(args["chain-mode"])
        # Parse optional 1D bootstrap parameters
        parse_csv(s) = isempty(s) ? Float64[] : parse.(Float64, split(s, ','))
        init_centers_t = parse_csv(args["chain-init-centers-t"])
        init_amps      = parse_csv(args["chain-init-amplitudes"])
        init_sigpar    = args["chain-init-sigma-parallel"]
        init_sigperp   = args["chain-init-sigma-perp"]

        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_baseline_quantile=args["chain-support-baseline-quantile"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                t_min_nm=args["chain-t-min-nm"],
                                t_max_nm=args["chain-t-max-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"],
                                init_centers_t=init_centers_t,
                                init_amplitudes=init_amps,
                                init_sigma_parallel=init_sigpar,
                                init_sigma_perp=init_sigperp,
                                skip_global=args["chain-skip-global"],
                                selection_criterion=args["selection-criterion"])
        if chain_mode in ("1d-bootstrapped", "direct")
            # Force n_min == n_max if bootstrapped with init data
            if !isempty(init_centers_t)
                n_boot = length(init_centers_t)
                ccfg.n_min = n_boot
                ccfg.n_max = n_boot
            end
            r, ctx = chain_direct_fit(img, cfg, ccfg)
            _write_chain_sweep([r], r, ctx, cfg, ccfg)
            cfg.no_plot || _plot_chain_sweep([r], r, ctx, cfg, ccfg)
            println("chain_direct_N: $(r.n)")
            println("chain_direct_valid: $(r.valid)")
            println("chain_direct_bic: $(r.bic)")
            println("chain_direct_bic: $(r.bic)")
            println("chain_direct_residual_peak_snr: $(r.residual_peak_snr)")
            println("chain_direct_train_nll: $(r.train_nll)")
            r.success && println("chain_direct_mean_spacing_nm: $(r.mean_spacing_nm)")
        else
            _run_chain_sweep(img, cfg, ccfg)
        end
        return
    end

    if args["consensus"]
        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_baseline_quantile=args["chain-support-baseline-quantile"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                t_min_nm=args["chain-t-min-nm"],
                                t_max_nm=args["chain-t-max-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"])
        consensus = fit_chain_consensus(img, cfg, ccfg)
        # Write Z results
        _write_chain_sweep(consensus.z.results, consensus.z.best, consensus.z.ctx, cfg, ccfg)
        cfg.no_plot || _plot_chain_sweep(consensus.z.results, consensus.z.best, consensus.z.ctx, cfg, ccfg)
        if consensus.current.best !== nothing
            println("consensus_z_N: $(consensus.z.best.n)")
            println("consensus_current_N: $(consensus.current.best.n)")
            println("consensus_agreement: $(consensus.agreement)")
            println("consensus_match: $(consensus.consensus)")
        end
        return
    end

    if !isempty(args["batch"])
        dir = args["batch"]
        files = filter(f -> endswith(lowercase(f), ".sxm"), readdir(dir; join=true))
        isempty(files) && error("No .sxm files found in $dir")
        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"])
        fit_chain_batch(files, cfg, ccfg; consensus=true)
        return
    end

    if cfg.inspect_images
        path = _plot_all_images(img, cfg)
        println("all_images_overview: $path")
        return
    end

    if cfg.all_images
        _run_all_images(img, cfg)
        return
    end

    result = fit_molecular_pattern(img, cfg)
    param_file = _write_parameters(result, cfg)
    plot_file = cfg.no_plot ? nothing : _plot_result(img, cfg, result)

    println("success: $(result.success)")
    println("estimated_repeats: $(result.estimated_repeats) (range $(result.estimated_repeat_range[1])-$(result.estimated_repeat_range[2]), ROI length $(@sprintf("%.3f", result.roi_length_nm)) nm)")
    println("chains: $(length(result.chains))")
    println("detected_features: $(length(result.raw_features))")
    println("fitted_features: $(length(result.features))")
    println("RSS: $(@sprintf("%.6g", result.rss))")
    println("R²: $(@sprintf("%.6f", result.r_squared))")
    println("BIC-like score: $(@sprintf("%.6g", result.bic))")
    println("parameters: $param_file")
    plot_file !== nothing && println("plot: $plot_file")
    for w in result.warnings
        println("warning: $w")
    end
end

end # module

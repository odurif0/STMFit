# Label-free repeated-scan grouping and molecule-level count consensus.
#
# Scans of one session are ordered by acquisition time. Two consecutive scans
# are the same, unchanged molecule when their flattened images, resampled in
# the absolute (piezo) frame, correlate above a fixed threshold after a pure
# translation (thermal drift). A molecule observed at least `min_track_scans`
# times with a strict-majority selected count gets that count in every scan
# whose frame fully contains the registered molecule footprint.
#
# Inputs are raw SXM headers/pixels and label-free per-scan counts/geometry.
# No expected count, sequence, class proportion or benchmark membership exists
# in this module.
module MoleculeConsensus

using STMSXMIO, Statistics, Dates, Printf

export ScanFrame, scan_frame, absolute_xy, image_xy, ConsensusSettings,
       load_consensus_settings, registration_image, resample_absolute, highpass,
       register_grids, link_scans, consensus_counts

struct ScanFrame
    file::String
    acquired::DateTime
    session::String
    offset::NTuple{2,Float64}
    extent::NTuple{2,Float64}
    angle_deg::Float64
    width::Int
    height::Int
end

function scan_frame(img::SXMImage)
    h = img.header
    for key in ("REC_DATE", "REC_TIME")
        haskey(h, key) || error("SXM header lacks $key: $(img.filepath)")
    end
    date = strip(h["REC_DATE"]); time = strip(h["REC_TIME"])
    acquired = DateTime(date * " " * time, dateformat"dd.mm.yyyy HH:MM:SS")
    angle = parse(Float64, strip(get(h, "SCAN_ANGLE", "0")))
    return ScanFrame(basename(img.filepath), acquired, date, img.offset_nm, img.range_nm,
                     angle, img.width, img.height)
end

# Image coordinates (x right, y down, origin at the first pixel, nm) to the
# absolute piezo frame. SCAN_OFFSET is the frame centre; SCAN_ANGLE rotates the
# frame clockwise. The sign convention is the one under which chain axes agree
# across the rotated repeat scans (a label-free geometric check); its global
# mirror image is equivalent for registration.
function absolute_xy(f::ScanFrame, x, y)
    lx = x - f.extent[1] / 2; ly = -(y - f.extent[2] / 2)
    th = -deg2rad(f.angle_deg); c, s = cos(th), sin(th)
    return f.offset[1] + c * lx - s * ly, f.offset[2] + s * lx + c * ly
end

function image_xy(f::ScanFrame, X, Y)
    th = -deg2rad(f.angle_deg); c, s = cos(th), sin(th)
    dx = X - f.offset[1]; dy = Y - f.offset[2]
    lx = c * dx + s * dy; ly = -s * dx + c * dy
    return lx + f.extent[1] / 2, -ly + f.extent[2] / 2
end

struct ConsensusSettings
    window_half_nm::Float64
    grid_step_nm::Float64
    highpass_sigma_nm::Float64
    max_shift_nm::Float64
    coarse_step_px::Int
    min_overlap_px::Int
    ncc_min::Float64
    max_center_distance_nm::Float64
    min_track_scans::Int
    frame_margin_nm::Float64
    flatten::String
    channel::String
    roi_margin_nm::Float64
end

function load_consensus_settings(cfg::AbstractDict)
    m = cfg["model"]; s = cfg["selection"]; p = cfg["preprocessing"]
    get(s, "majority", "") == "strict" || error("Only the strict-majority consensus is defined")
    out = ConsensusSettings(Float64(m["window_half_nm"]), Float64(m["grid_step_nm"]),
        Float64(m["highpass_sigma_nm"]), Float64(m["max_shift_nm"]), Int(m["coarse_step_px"]),
        Int(s["min_overlap_px"]), Float64(s["ncc_min"]), Float64(s["max_center_distance_nm"]),
        Int(s["min_track_scans"]), Float64(s["frame_margin_nm"]), String(p["flatten"]), String(p["channel"]),
        Float64(m["roi_margin_nm"]))
    out.window_half_nm > 0 && out.grid_step_nm > 0 && out.highpass_sigma_nm > 0 &&
        out.max_shift_nm > 0 && out.coarse_step_px >= 1 && out.min_overlap_px > 0 &&
        -1 < out.ncc_min < 1 && out.max_center_distance_nm > 0 && out.min_track_scans >= 2 &&
        out.frame_margin_nm >= 0 && out.roi_margin_nm > 0 || error("Invalid molecule-consensus settings")
    return out
end

# Mean of the flattened forward and backward topography (nm). Missing pixels
# stay missing; no imputation.
function registration_image(img::SXMImage, st::ConsensusSettings)
    views = Matrix{Float64}[]
    for dir in ("fwd", "bwd")
        ch = nothing
        for c in img.channels
            lowercase(c.name) == lowercase(st.channel) && lowercase(c.direction) == dir && (ch = c)
        end
        ch === nothing && continue
        pre = STMSXMIO.preprocess_observed_channel(img, ch; stride=1, flatten=st.flatten,
                                          smooth_radius_px=0, plane_rank_rtol=1e-12)
        push!(views, pre.z)
    end
    isempty(views) && error("No $(st.channel) channel in $(img.filepath)")
    out = fill(NaN, size(views[1]))
    for i in eachindex(out)
        vals = [v[i] for v in views if isfinite(v[i])]
        isempty(vals) || (out[i] = mean(vals))
    end
    return out
end

function _bilinear(z::AbstractMatrix{Float64}, row::Float64, col::Float64)
    h, w = size(z)
    (1 <= row <= h && 1 <= col <= w) || return NaN
    r0 = min(floor(Int, row), h - 1); c0 = min(floor(Int, col), w - 1)
    r0 = max(r0, 1); c0 = max(c0, 1)
    fr = row - r0; fc = col - c0
    a = z[r0, c0]; b = z[r0, c0+1]; c = z[r0+1, c0]; d = z[r0+1, c0+1]
    (isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d)) || return NaN
    return (1 - fr) * ((1 - fc) * a + fc * b) + fr * ((1 - fc) * c + fc * d)
end

# Grid rows follow absolute Y, columns absolute X, both centred on `center`.
function resample_absolute(z::AbstractMatrix{Float64}, f::ScanFrame, center, st::ConsensusSettings)
    g = collect(-st.window_half_nm:st.grid_step_nm:st.window_half_nm + st.grid_step_nm / 2)
    out = fill(NaN, length(g), length(g))
    sx = (f.width - 1) / f.extent[1]; sy = (f.height - 1) / f.extent[2]
    for (i, gy) in enumerate(g), (j, gx) in enumerate(g)
        x, y = image_xy(f, center[1] + gx, center[2] + gy)
        out[i, j] = _bilinear(z, y * sy + 1, x * sx + 1)
    end
    return out
end

function _blur1(a::Matrix{Float64}, k::Vector{Float64}, dim::Int)
    r = (length(k) - 1) ÷ 2; out = zeros(size(a)); n = size(a, dim)
    for idx in CartesianIndices(a)
        i = idx[dim]; acc = 0.0
        for t in -r:r
            j = i + t
            1 <= j <= n || continue
            src = dim == 1 ? a[j, idx[2]] : a[idx[1], j]
            acc += k[t+r+1] * src
        end
        out[idx] = acc
    end
    return out
end

# Remove structure broader than sigma with a missing-aware Gaussian estimate.
function highpass(a::Matrix{Float64}, sigma_px::Float64)
    r = ceil(Int, 3sigma_px); k = [exp(-t^2 / (2sigma_px^2)) for t in -r:r]; k ./= sum(k)
    m = Float64.(isfinite.(a)); b = ifelse.(isfinite.(a), a, 0.0)
    num = _blur1(_blur1(b, k, 1), k, 2); den = _blur1(_blur1(m, k, 1), k, 2)
    out = fill(NaN, size(a))
    for i in eachindex(a)
        isfinite(a[i]) && den[i] > 1e-6 && (out[i] = a[i] - num[i] / den[i])
    end
    return out
end

function _ncc(A::Matrix{Float64}, B::Matrix{Float64}, dy::Int, dx::Int, min_overlap::Int)
    n1, n2 = size(A); sa = sb = saa = sbb = sab = 0.0; n = 0
    for i in max(1, 1 - dy):min(n1, n1 - dy), j in max(1, 1 - dx):min(n2, n2 - dx)
        a = A[i, j]; b = B[i+dy, j+dx]
        (isfinite(a) && isfinite(b)) || continue
        n += 1; sa += a; sb += b; saa += a * a; sbb += b * b; sab += a * b
    end
    n >= min_overlap || return -Inf
    ca = saa - sa^2 / n; cb = sbb - sb^2 / n
    (ca > 0 && cb > 0) || return -Inf
    return (sab - sa * sb / n) / sqrt(ca * cb)
end

# Translation t (nm, absolute frame) such that content at X in A appears at
# X + t in B, with its normalized cross-correlation.
function register_grids(A::Matrix{Float64}, B::Matrix{Float64}, st::ConsensusSettings)
    S = round(Int, st.max_shift_nm / st.grid_step_nm); c = st.coarse_step_px
    best = (-Inf, 0, 0)
    for dy in -S:c:S, dx in -S:c:S
        v = _ncc(A, B, dy, dx, st.min_overlap_px)
        v > best[1] && (best = (v, dy, dx))
    end
    _, y0, x0 = best
    for dy in max(-S, y0 - c):min(S, y0 + c), dx in max(-S, x0 - c):min(S, x0 + c)
        v = _ncc(A, B, dy, dx, st.min_overlap_px)
        v > best[1] && (best = (v, dy, dx))
    end
    v, dy, dx = best
    return (ncc=v, tx=dx * st.grid_step_nm, ty=dy * st.grid_step_nm)
end

"""
    link_scans(frames, centroids, images, st)

`frames`: Dict file => ScanFrame; `centroids`: Dict file => absolute chain
centroid; `images`: function file => registration image. Returns the
acquisition order and one row per consecutive pair.
"""
function link_scans(frames, centroids, images, st::ConsensusSettings)
    order = sort(collect(keys(frames)); by=f -> (frames[f].acquired, f))
    pairs = NamedTuple[]
    sigma_px = st.highpass_sigma_nm / st.grid_step_nm
    for k in 1:length(order)-1
        a, b = order[k], order[k+1]; fa, fb = frames[a], frames[b]
        dist = hypot((centroids[a] .- centroids[b])...)
        reason = fa.session != fb.session ? "different_session" :
                 dist > st.max_center_distance_nm ? "centroids_apart" : ""
        reg = (ncc=NaN, tx=NaN, ty=NaN)
        if isempty(reason)
            A = highpass(resample_absolute(images(a), fa, centroids[a], st), sigma_px)
            B = highpass(resample_absolute(images(b), fb, centroids[a], st), sigma_px)
            reg = register_grids(A, B, st)
            reason = isfinite(reg.ncc) && reg.ncc >= st.ncc_min ? "linked" : "low_correlation"
        end
        push!(pairs, (a=a, b=b, center_distance_nm=dist, ncc=reg.ncc, tx_nm=reg.tx, ty_nm=reg.ty,
                      linked=reason == "linked", reason=reason))
    end
    return order, pairs
end

"""
    consensus_counts(order, pairs, counts, lobes_abs, frames, st)

`counts`: Dict file => per-scan selected N; `lobes_abs`: Dict file => Vector of
absolute lobe positions (label-free per-scan geometry). Returns one row per scan.
"""
function consensus_counts(order, pairs, counts, lobes_abs, frames, st::ConsensusSettings)
    track = Dict{String,Int}(); tid = 0
    link = Dict((p.a, p.b) => p for p in pairs if p.linked)
    for (k, f) in enumerate(order)
        if k > 1 && haskey(link, (order[k-1], f))
            track[f] = track[order[k-1]]
        else
            tid += 1; track[f] = tid
        end
    end
    members = Dict{Int,Vector{String}}()
    for f in order; push!(get!(members, track[f], String[]), f); end
    pos = Dict(f => k for (k, f) in enumerate(order))
    rows = Dict{String,Any}[]
    for f in order
        mem = members[track[f]]; n = length(mem); ns = [counts[m] for m in mem]
        tally = Dict{Int,Int}(); for v in ns; tally[v] = get(tally, v, 0) + 1; end
        top = maximum(values(tally)); winners = [k for (k, v) in tally if v == top]
        nstar = n >= st.min_track_scans && length(winners) == 1 && top > n / 2 ? only(winners) : 0
        row = Dict{String,Any}("file" => f, "track" => track[f], "track_size" => n,
            "N_scan" => counts[f], "N_consensus" => nstar, "agreement" => top / n,
            "N_final" => counts[f], "reference_scan" => "", "rule" => "",
            "roi_x0_nm" => NaN, "roi_y0_nm" => NaN, "roi_x1_nm" => NaN, "roi_y1_nm" => NaN)
        if n < st.min_track_scans
            row["rule"] = "track_too_short"
        elseif nstar == 0
            row["rule"] = "no_strict_majority"
        elseif counts[f] == nstar
            row["rule"] = "agrees"
        else
            refs = [m for m in mem if counts[m] == nstar]
            ref = refs[argmin([(abs(pos[m] - pos[f]), pos[m]) for m in refs])]
            t = (0.0, 0.0); lo, hi = minmax(pos[ref], pos[f]); sgn = pos[ref] < pos[f] ? 1.0 : -1.0
            for k in lo:hi-1
                p = link[(order[k], order[k+1])]; t = (t[1] + sgn * p.tx_nm, t[2] + sgn * p.ty_nm)
            end
            fr = frames[f]; m = st.frame_margin_nm
            pts = [image_xy(fr, P[1] + t[1], P[2] + t[2]) for P in lobes_abs[ref]]
            inside = all(p -> m <= p[1] <= fr.extent[1] - m && m <= p[2] <= fr.extent[2] - m, pts)
            r = st.roi_margin_nm
            row["roi_x0_nm"] = max(0.0, minimum(first.(pts)) - r); row["roi_x1_nm"] = min(fr.extent[1], maximum(first.(pts)) + r)
            row["roi_y0_nm"] = max(0.0, minimum(last.(pts)) - r); row["roi_y1_nm"] = min(fr.extent[2], maximum(last.(pts)) + r)
            row["reference_scan"] = ref
            if inside
                row["N_final"] = nstar; row["rule"] = "consensus_applied"
            else
                row["rule"] = "footprint_outside_frame"
            end
        end
        push!(rows, row)
    end
    return rows
end

end

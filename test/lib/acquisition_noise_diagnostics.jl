module AcquisitionNoiseDiagnostics

# Opt-in, label-free acquisition evidence. No optimizer, selector, class score,
# effective sample size, or likelihood is defined here. July 22/26 diagnostics
# already established native x reversal, residual hysteresis and boundary failures.
using Statistics
using LinearAlgebra
using TOML
using SHA
using STMSXMIO
using GaussianFit2D

# Reuse the production artifact reader's adaptive-support decision, not N alone.
# extract_lobe_features.jl guards main; none of its fitting functions is called.
module SavedSupport
include(joinpath(@__DIR__, "..", "extract_lobe_features.jl"))
end

export settings_from_dict, load_settings, pair_metrics, common_support,
       registration_scan, lag_acf, background_diagnostic, block_feasibility,
       read_geometry, load_case, diagnose, write_diagnostic

const REQUIRED_SETTINGS = (
    "channel", "auxiliary_channel", "lag_x_min_px", "lag_x_max_px", "lag_y_min_px", "lag_y_max_px",
    "lag_bound_basis", "min_registration_pixels",
    "min_positive_correlation", "ambiguity_correlation_gap", "min_row_pixels",
    "footprint_sigma", "background_guard_nm", "min_background_pixels",
    "min_background_fraction", "acf_max_lag_px", "min_acf_pairs", "acf_threshold",
    "block_length_multiplier", "min_usable_blocks", "block_min_occupancy",
    "local_patch_radius_nm", "min_patch_pixels")

function settings_from_dict(d::AbstractDict)
    for k in REQUIRED_SETTINGS
        haskey(d, k) || error("Missing acquisition_noise.$k; diagnostics have no physical defaults")
    end
    ints = ("lag_x_min_px", "lag_x_max_px", "lag_y_min_px", "lag_y_max_px",
            "min_registration_pixels", "min_row_pixels", "min_background_pixels",
            "acf_max_lag_px", "min_acf_pairs", "min_usable_blocks", "min_patch_pixels")
    for k in ints
        d[k] isa Integer && !(d[k] isa Bool) || error("$k must be an integer")
    end
    for k in ("min_registration_pixels", "min_row_pixels", "min_background_pixels",
              "acf_max_lag_px", "min_acf_pairs", "min_usable_blocks", "min_patch_pixels")
        d[k] >= (k in ("acf_max_lag_px", "min_usable_blocks") ? 1 : 3) || error("$k is too small")
    end
    d["lag_x_min_px"] <= 0 <= d["lag_x_max_px"] || error("x lag grid must include zero")
    d["lag_y_min_px"] <= 0 <= d["lag_y_max_px"] || error("y lag grid must include zero")
    d["lag_bound_basis"] isa AbstractString && !isempty(strip(d["lag_bound_basis"])) ||
        error("lag_bound_basis must explain the fixed diagnostic window")
    for k in ("min_positive_correlation", "ambiguity_correlation_gap", "footprint_sigma",
              "background_guard_nm", "min_background_fraction", "acf_threshold",
              "block_length_multiplier", "block_min_occupancy", "local_patch_radius_nm")
        d[k] isa Real && !(d[k] isa Bool) && isfinite(d[k]) || error("$k must be finite")
    end
    0 < d["min_positive_correlation"] <= 1 || error("min_positive_correlation must be in (0,1]")
    0 <= d["ambiguity_correlation_gap"] < 1 || error("Invalid ambiguity_correlation_gap")
    0 < d["min_background_fraction"] <= 1 || error("Invalid min_background_fraction")
    0 < d["acf_threshold"] < 1 || error("Invalid acf_threshold")
    0 < d["block_min_occupancy"] <= 1 || error("Invalid block_min_occupancy")
    d["footprint_sigma"] > 0 && d["local_patch_radius_nm"] > 0 || error("Footprint/patch size must be positive")
    d["background_guard_nm"] >= 0 || error("background_guard_nm must be nonnegative")
    d["block_length_multiplier"] >= 1 || error("block_length_multiplier must be at least one")
    for k in ("channel", "auxiliary_channel")
        d[k] isa AbstractString && !isempty(strip(d[k])) || error("$k must be named")
    end
    return Dict{String,Any}(k => d[k] for k in REQUIRED_SETTINGS)
end

function load_settings(path::AbstractString)
    d = TOML.parsefile(path)
    haskey(d, "acquisition_noise") || error("Missing [acquisition_noise] settings")
    return settings_from_dict(d["acquisition_noise"])
end

"""Signed correlation and physical positive-affine discrepancy, b ≈ gain*a+offset.
NRMSE denominators are population centered RMS, not a fitted noise estimate.
A negative-correlation image is never made to agree by flipping its sign.
"""
function pair_metrics(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || throw(DimensionMismatch("Paired vectors differ"))
    n = length(a)
    empty = (n=n, status="insufficient_or_constant", correlation=NaN, gain=NaN,
             offset=NaN, direct_nrmse=NaN, positive_affine_nrmse=NaN,
             difference_rms=NaN, difference_median=NaN)
    n >= 3 || return empty
    all(isfinite, a) && all(isfinite, b) || error("pair_metrics requires finite pairs")
    ma, mb = mean(a), mean(b)
    ac, bc = a .- ma, b .- mb
    va, vb = sum(abs2, ac), sum(abs2, bc)
    va > 0 && vb > 0 || return empty
    c = clamp(dot(ac, bc) / sqrt(va * vb), -1.0, 1.0)
    gain = max(0.0, dot(ac, bc) / va)
    offset = mb - gain * ma
    return (n=n, status="ok", correlation=c, gain=gain, offset=offset,
            direct_nrmse=sqrt(sum(abs2, a .- b) / va),
            positive_affine_nrmse=sqrt(sum(abs2, b .- (gain .* a .+ offset)) / vb),
            difference_rms=sqrt(mean(abs2, a .- b)), difference_median=median(b .- a))
end

"""Forward locations valid at EVERY candidate backward sample (y+dy,x+dx).
Support is fixed before scoring. Never compare changing overlaps or fill edges.
"""
function common_support(a::AbstractMatrix, b::AbstractMatrix, mask::AbstractMatrix{Bool},
                        dxs, dys)
    size(a) == size(b) == size(mask) || throw(DimensionMismatch("Image/support sizes differ"))
    (isempty(dxs) || isempty(dys)) && error("Empty lag grid")
    out = falses(size(mask))
    ny, nx = size(mask)
    for I in findall(mask)
        y, x = Tuple(I)
        isfinite(a[I]) || continue
        valid = true
        for dy in dys, dx in dxs
            yy, xx = y + dy, x + dx
            if !(1 <= yy <= ny && 1 <= xx <= nx && isfinite(b[yy, xx]))
                valid = false
                break
            end
        end
        out[I] = valid
    end
    return out
end

function sample_pairs(a, b, support, dx::Int, dy::Int)
    ii = findall(support)
    return a[ii], [b[I[1]+dy, I[2]+dx] for I in ii]
end

function scan_fixed(a, b, support, dxs, dys, s; min_pixels::Int)
    candidates = NamedTuple[]
    n = count(support)
    if n < min_pixels
        return (status="insufficient_support", usable=false, n=n, dx=0, dy=0,
                best_correlation=NaN, correlation_gap=NaN, near_maxima=0,
                boundary=false, ambiguous=false, candidates=candidates)
    end
    indices = findall(support)
    forward_values = a[indices]
    for dy in dys, dx in dxs
        backward_values = [b[I[1]+dy, I[2]+dx] for I in indices]
        m = pair_metrics(forward_values, backward_values)
        push!(candidates, merge((dx_px=dx, dy_px=dy), m))
    end
    valid = filter(r -> isfinite(r.correlation), candidates)
    if isempty(valid)
        return (status="constant_signal", usable=false, n=n, dx=0, dy=0,
                best_correlation=NaN, correlation_gap=NaN, near_maxima=0,
                boundary=false, ambiguous=false, candidates=candidates)
    end
    # Deterministic tie order only; ambiguity is always reported, never hidden.
    sort!(valid; by=r -> (-r.correlation, abs(r.dx_px)+abs(r.dy_px), r.dy_px, r.dx_px))
    best = first(valid)
    gap = length(valid) > 1 ? best.correlation - valid[2].correlation : Inf
    near = count(r -> best.correlation-r.correlation <= s["ambiguity_correlation_gap"], valid)
    boundary = (length(dxs)>1 && best.dx_px in (first(dxs), last(dxs))) ||
               (length(dys)>1 && best.dy_px in (first(dys), last(dys)))
    ambiguous = near > 1
    reasons = String[]
    best.correlation <= 0 && push!(reasons, "nonpositive_correlation")
    0 < best.correlation < s["min_positive_correlation"] && push!(reasons, "weak_positive_correlation")
    ambiguous && push!(reasons, "ambiguous_maximum")
    boundary && push!(reasons, "boundary_maximum")
    return (status=isempty(reasons) ? "interior_positive_maximum" : join(reasons, ";"),
            usable=isempty(reasons), n=n, dx=best.dx_px, dy=best.dy_px,
            best_correlation=best.correlation, correlation_gap=gap, near_maxima=near,
            boundary=boundary, ambiguous=ambiguous, candidates=candidates)
end

function registration_scan(a, b, mask, s)
    dxs = s["lag_x_min_px"]:s["lag_x_max_px"]
    dys = s["lag_y_min_px"]:s["lag_y_max_px"]
    support = common_support(a, b, mask, dxs, dys)
    result = scan_fixed(a, b, support, dxs, dys, s; min_pixels=s["min_registration_pixels"])
    return merge(result, (support=support,))
end

function lag_acf(z::AbstractMatrix, mask::AbstractMatrix{Bool}, s)
    size(z) == size(mask) || throw(DimensionMismatch("Noise/support sizes differ"))
    rows = NamedTuple[]
    ny, nx = size(z)
    for axis in ("x", "y"), lag in 0:s["acf_max_lag_px"]
        a, b = Float64[], Float64[]
        dy, dx = axis == "x" ? (0, lag) : (lag, 0)
        for y in 1:max(0, ny-dy), x in 1:max(0, nx-dx)
            mask[y,x] && mask[y+dy,x+dx] || continue
            isfinite(z[y,x]) && isfinite(z[y+dy,x+dx]) || continue
            push!(a, z[y,x]); push!(b, z[y+dy,x+dx])
        end
        ok = length(a) >= s["min_acf_pairs"]
        m = ok ? pair_metrics(a, b) : nothing
        push!(rows, (axis=axis, lag_px=lag, pairs=length(a),
                     correlation=ok ? m.correlation : NaN,
                     status=ok ? m.status : "insufficient_pairs"))
    end
    return rows
end

function acf_length(rows, axis::String, threshold::Real)
    r = filter(r -> r.axis == axis && r.lag_px > 0, rows)
    all(row -> isfinite(row.correlation), r) ||
        return (pixels=NaN, status="insufficient_pairs_or_constant")
    # Negative correlation is dependence too. An early zero crossing in an
    # oscillatory scan artefact does not justify small validation blocks.
    first_crossing = findfirst(row -> abs(row.correlation) <= threshold, r)
    first_crossing === nothing && return (pixels=NaN, status="right_censored_at_max_lag")
    any(row -> abs(row.correlation) > threshold, r[first_crossing:end]) &&
        return (pixels=NaN, status="resurgent_or_oscillatory_acf_tail")
    return (pixels=Float64(r[first_crossing].lag_px), status="first_magnitude_threshold_crossing_no_observed_resurgence")
end

"""Off-footprint image diagnostics only. A fitted residual is not background noise.
MAD is descriptive: substrate structure/line artefacts need not be stochastic noise.
"""
function background_diagnostic(z, mask, s)
    valid = mask .& isfinite.(z)
    n = count(valid)
    fraction = n / length(z)
    adequate = n >= s["min_background_pixels"] && fraction >= s["min_background_fraction"]
    if !adequate
        return (status="insufficient_background", n=n, fraction=fraction,
                mad_sigma=NaN, standard_deviation=NaN, mean=NaN,
                length_x=(pixels=NaN, status="insufficient_background"),
                length_y=(pixels=NaN, status="insufficient_background"), acf=NamedTuple[])
    end
    v = z[valid]
    rows = lag_acf(z, valid, s)
    return (status="off_footprint_estimate_not_verified_pure_noise", n=n, fraction=fraction,
            mad_sigma=1.4826*median(abs.(v .- median(v))), standard_deviation=std(v), mean=mean(v),
            length_x=acf_length(rows, "x", s["acf_threshold"]),
            length_y=acf_length(rows, "y", s["acf_threshold"]), acf=rows)
end

"""Count disjoint occupied rectangles; this is NOT an effective sample size.
First |ACF| crossings and finite block counts do not establish independence.
"""
function block_feasibility(mask, length_x::Real, length_y::Real, s)
    if !(isfinite(length_x) && isfinite(length_y) && length_x > 0 && length_y > 0)
        return (status="blocked_unresolved_correlation_length", block_x_px=0,
                block_y_px=0, usable_blocks=0, feasible=false)
    end
    bx = ceil(Int, s["block_length_multiplier"]*length_x)
    by = ceil(Int, s["block_length_multiplier"]*length_y)
    ny, nx = size(mask)
    usable = 0
    for y in 1:by:(ny-by+1), x in 1:bx:(nx-bx+1)
        count(@view(mask[y:y+by-1, x:x+bx-1])) / (bx*by) >= s["block_min_occupancy"] && (usable += 1)
    end
    feasible = usable >= s["min_usable_blocks"]
    return (status=feasible ? "descriptive_blocks_available_not_independence" : "blocked_too_few_occupied_blocks",
            block_x_px=bx, block_y_px=by, usable_blocks=usable, feasible=feasible)
end

function read_geometry(path::AbstractString, file::AbstractString, n::Int)
    lines = readlines(path)
    isempty(lines) && error("Empty geometry TSV")
    header = String.(strip.(split(first(lines), '\t'; keepempty=true)))
    length(unique(header)) == length(header) || error("Duplicate geometry headers")
    required = ("file", "N", "lobe", "amplitude", "x_nm", "y_nm", "sigma_parallel_nm",
                "sigma_perp_nm", "axis_x", "axis_y", "origin_x_nm", "origin_y_nm",
                "baseline", "tilt_x", "tilt_y")
    all(k -> k in header, required) || error("Geometry requires actual saved lobe centers, widths, axis, baseline and tilt")
    forbidden = ("truth", "sequence", "expected_N", "target_N", "control_sequence", "label", "predicted")
    any(in(header), forbidden) && error("Geometry must not contain external labels or predictions")
    rows = Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = String.(strip.(split(line, '\t'; keepempty=true)))
        length(vals) == length(header) || error("Malformed geometry row")
        row = Dict(zip(header, vals))
        basename(row["file"]) == basename(file) && push!(rows, row)
    end
    length(rows) == n || error("Need exactly $n saved geometry rows for file; ambiguous/missing candidate geometry")
    all(r -> tryparse(Int, r["N"]) == n, rows) || error("Saved geometry N differs from selected summary")
    ids = [tryparse(Int, r["lobe"]) for r in rows]
    Set(ids) == Set(1:n) || error("Duplicate/missing/invalid lobe indices")
    sort!(rows; by=r -> parse(Int, r["lobe"]))
    numeric = setdiff(collect(required), ["file", "N", "lobe"])
    lobes = NamedTuple[]
    for r in rows
        values = Dict{String,Float64}()
        for k in numeric
            v = tryparse(Float64, r[k])
            v !== nothing && isfinite(v) || error("Non-finite geometry field $k")
            values[k] = v
        end
        values["sigma_parallel_nm"] > 0 && values["sigma_perp_nm"] > 0 || error("Nonpositive geometry width")
        values["amplitude"] >= 0 || error("Negative geometry amplitude")
        skew = parse(Float64, get(r, "skew_ratio", "1.0"))
        isfinite(skew) && skew > 0 || error("Invalid geometry skew")
        push!(lobes, (; (Symbol(k) => values[k] for k in numeric)..., lobe=parse(Int,r["lobe"]), skew_ratio=skew))
    end
    common = (:axis_x, :axis_y, :origin_x_nm, :origin_y_nm, :baseline, :tilt_x, :tilt_y)
    all(l -> all(k -> getproperty(l,k) == getproperty(first(lobes),k), common), lobes) ||
        error("Saved geometry has inconsistent axis/origin/baseline")
    l = first(lobes)
    abs(hypot(l.axis_x, l.axis_y)-1) <= 1e-6 || error("Saved geometry axis is not normalized")
    return lobes
end

function strict_channel(img, name, dir)
    matches = filter(c -> lowercase(c.name)==lowercase(name) && lowercase(c.direction)==dir, img.channels)
    length(matches) == 1 || error("Exactly one $name $dir channel is required; never reuse fallback direction")
    return only(matches)
end

# Auxiliary channel is descriptive only: no optimizer, registration search,
# nuisance projection, fit weights, noise model or independence assumption.
function load_auxiliary(img, pcfg, name, direction)
    matches = filter(c -> lowercase(c.name)==lowercase(name) && lowercase(c.direction)==direction, img.channels)
    empty = (available=length(matches)==1, matching_channels=length(matches), unit="",
             analysis_pixels=0, finite_pixels=0, raw_min=NaN, raw_max=NaN,
             raw_std=NaN, image=nothing)
    length(matches)==1 || return merge(empty,(status=isempty(matches) ? "missing_channel" : "duplicate_direction_channels",))
    ch=only(matches)
    scale,unit=STMSXMIO._value_scale(ch.unit)
    raw=ch.data[1:pcfg.stride:end,1:pcfg.stride:end].*scale
    finite=isfinite.(raw)
    values=raw[finite]
    base=merge(empty,(unit=unit,analysis_pixels=length(raw),finite_pixels=length(values)))
    isempty(values) && return merge(base,(status="all_nonfinite_channel",))
    _,_,_,z,_,_,_=GaussianFit2D.preprocess_channel(img,ch,pcfg)
    z[.!finite].=NaN # Native median imputation is not an observed Current sample.
    return merge(base,(status=all(finite) ? "available" : "available_with_nonfinite_samples",
        raw_min=minimum(values),raw_max=maximum(values),raw_std=length(values)>1 ? std(values) : NaN,image=z))
end

function auxiliary_report(case,reg,s)
    af,ab=case.auxiliary.forward,case.auxiliary.backward
    rows=NamedTuple[]
    n=reg.n
    blank=(channel=s["auxiliary_channel"],comparison="",condition="",source_unit="",target_unit="",
           dx_px=0,dy_px=0,z_support_pixels=n,observed_pairs=0,status="",
           correlation=NaN,linear_r2=NaN,positive_affine_nrmse=NaN,
           gain=NaN,offset=NaN,analysis_pixels=0,finite_pixels=0,
           raw_min=NaN,raw_max=NaN,raw_std=NaN)
    for (direction,a) in (("fwd",af),("bwd",ab))
        push!(rows,merge(blank,(comparison="direction_scale",condition=direction,status=a.status,
            target_unit=a.unit,analysis_pixels=a.analysis_pixels,finite_pixels=a.finite_pixels,
            raw_min=a.raw_min,raw_max=a.raw_max,raw_std=a.raw_std)))
    end
    comparisons=(("auxiliary_fwd_bwd","native",af.image,ab.image,0,0,af.unit,ab.unit),
                 ("auxiliary_fwd_bwd","at_fixed_Z_lag",af.image,ab.image,reg.dx,reg.dy,af.unit,ab.unit),
                 ("Z_auxiliary_fwd","native",case.forward,af.image,0,0,case.unit,af.unit),
                 ("Z_auxiliary_bwd","native",case.backward,ab.image,0,0,case.unit,ab.unit),
                 ("Z_auxiliary_bwd","at_fixed_Z_lag",case.backward,ab.image,reg.dx,reg.dy,case.unit,ab.unit))
    ii=findall(reg.support)
    for (comparison,condition,a,b,dx,dy,u1,u2) in comparisons
        row=merge(blank,(comparison=comparison,condition=condition,dx_px=dx,dy_px=dy,source_unit=u1,target_unit=u2))
        if a===nothing || b===nothing
            push!(rows,merge(row,(status="unavailable_direction",)))
            continue
        elseif comparison=="auxiliary_fwd_bwd" && u1!=u2
            push!(rows,merge(row,(status="incompatible_direction_units",)))
            continue
        end
        # For same-direction backward coupling, sample both channels at the same
        # backward location. For acquisition concordance, keep forward fixed.
        shift_a=comparison=="Z_auxiliary_bwd"
        va=[a[I[1]+(shift_a ? dy : 0),I[2]+(shift_a ? dx : 0)] for I in ii]
        vb=[b[I[1]+dy,I[2]+dx] for I in ii]
        observed=count(isfinite.(va).&isfinite.(vb))
        if observed!=n
            push!(rows,merge(row,(status="nonfinite_on_fixed_Z_support",observed_pairs=observed)))
        elseif n<s["min_registration_pixels"]
            push!(rows,merge(row,(status="insufficient_fixed_Z_support",observed_pairs=observed)))
        else
            m=pair_metrics(va,vb)
            coupled=comparison!="auxiliary_fwd_bwd"
            # Z and Current have different units: no direct difference/NRMSE is
            # meaningful between them. r² is descriptive signed-Pearson squared,
            # equivalent to one unconstrained linear predictor plus an intercept.
            push!(rows,merge(row,(status=m.status,observed_pairs=n,correlation=m.correlation,
                linear_r2=isfinite(m.correlation) ? m.correlation^2 : NaN,
                positive_affine_nrmse=coupled ? NaN : m.positive_affine_nrmse,
                gain=coupled ? NaN : m.gain,offset=coupled ? NaN : m.offset)))
        end
    end
    return rows
end

function load_case(raw::AbstractString, geometry::AbstractString, physical::AbstractString,
                   summary::AbstractString, s)
    cfg = TOML.parsefile(physical)
    model, pre = cfg["model"], cfg["preprocessing"]
    # Require actual physical/preprocessing values used by support reconstruction.
    for k in ("stride", "flatten", "smooth_radius_px")
        haskey(pre, k) || error("Frozen physical config missing preprocessing.$k")
    end
    pre["stride"] isa Integer && !(pre["stride"] isa Bool) && pre["stride"] > 0 || error("Invalid frozen stride")
    pre["smooth_radius_px"] isa Integer && !(pre["smooth_radius_px"] isa Bool) && pre["smooth_radius_px"] >= 0 || error("Invalid frozen smooth radius")
    lowercase(pre["flatten"]) in ("none", "plane", "rows", "plane+rows") || error("Unsupported frozen flatten mode")
    for k in ("fit_width_nm", "support_noise_k", "support_padding_nm", "support_min_length_nm", "support_baseline_quantile")
        haskey(model, k) || error("Frozen physical config missing model.$k")
    end
    header = String.(strip.(split(first(readlines(summary)), '\t'; keepempty=true)))
    any(k -> k in header, ("truth", "expected_N", "target_N", "sequence", "control_sequence")) &&
        error("Selected summary must be label-free")
    contexts = SavedSupport._read_selected_context(String(summary);
        selection_policy=String(get(model, "selection_policy", "")))
    file = basename(raw)
    haskey(contexts, file) || error("No selected summary for $file")
    selected = contexts[file]
    lobes = read_geometry(geometry, file, selected.n)
    img = STMSXMIO.read_sxm(String(raw)) # Native backward x flip happens exactly once.
    pcfg, ccfg, _ = SavedSupport._configs(model, pre, ""; selected_context=selected)
    # Selected chain geometry and support were built from Z; other channel
    # likelihoods need a separate, explicitly justified geometry model.
    lowercase(s["channel"]) == lowercase(pcfg.roi_channel) ||
        error("This saved geometry supports only its native ROI channel $(pcfg.roi_channel)")
    fch = strict_channel(img, s["channel"], "fwd")
    bch = strict_channel(img, s["channel"], "bwd")
    fch.unit == bch.unit || error("Forward/backward units differ")
    xs, ys, _, f, _, unit, _ = GaussianFit2D.preprocess_channel(img, fch, pcfg)
    _, _, _, b, _, _, _ = GaussianFit2D.preprocess_channel(img, bch, pcfg)
    length(xs) >= 2 && length(ys) >= 2 || error("Need at least two grid locations on each image axis")
    all(l -> first(xs) <= l.x_nm <= last(xs) && first(ys) <= l.y_nm <= last(ys), lobes) ||
        error("Saved lobe centers are outside the raw image coordinate frame")
    profile = String(get(model, "peak_profile", "gaussian"))
    profile == "gaussian" && any(l -> l.skew_ratio != 1.0, lobes) &&
        error("Asymmetric saved geometry does not match the frozen Gaussian profile")
    stride = pcfg.stride
    stride > 0 || error("Invalid frozen stride")
    valid_f = isfinite.(fch.data[1:stride:end,1:stride:end])
    valid_b = isfinite.(bch.data[1:stride:end,1:stride:end])
    # The fused fit now uses both unsmoothed views; neither direction spreads
    # an imputed sample through a smoothing kernel in the fit map.
    valid_fused = valid_f .& valid_b
    all(isfinite, f) && all(isfinite, b) || error("No finite preprocessed views")
    xg, yg, fused, roi, x, y, z, _ = GaussianFit2D._fused_roi_data(img, pcfg)
    xs == xg && ys == yg || error("Native preprocessing grids differ")
    l = first(lobes)
    axis = (l.axis_x, l.axis_y)
    t = (x .- l.origin_x_nm).*axis[1] .+ (y .- l.origin_y_nm).*axis[2]
    ac = (origin=(l.origin_x_nm,l.origin_y_nm), axis=axis, perp=(-axis[2],axis[1]),
          tmin=minimum(t), tmax=maximum(t))
    _, _, _, _, keep, support_meta = GaussianFit2D._chain_fit_data(x, y, z, ac, ccfg)
    # _flatten_roi is row-major, whereas Julia findall(mask) is column-major.
    support = falses(size(roi))
    k = 0
    for iy in eachindex(ys), ix in eachindex(xs)
        roi[iy,ix] || continue
        k += 1
        support[iy,ix] = keep[k]
    end
    k == length(keep) || error("Native fit-support order mismatch")
    f[.!valid_f] .= NaN; b[.!valid_b] .= NaN
    auxiliary=(forward=load_auxiliary(img,pcfg,s["auxiliary_channel"],"fwd"),
               backward=load_auxiliary(img,pcfg,s["auxiliary_channel"],"bwd"))
    return (file=file, xs=xs, ys=ys, forward=f, backward=b, fused=fused, auxiliary=auxiliary,
            fit_support=support, fused_valid=valid_fused, roi=roi, lobes=lobes, selected=selected,
            support_meta=support_meta, unit=unit, physical_config=cfg,
            effective_support_noise_k=ccfg.support_noise_k,
            effective_support_padding_nm=ccfg.support_padding_nm,
            imputed_fwd=count(.!valid_f), imputed_bwd=count(.!valid_b),
            profile=profile)
end

function geometry_footprint(case, s)
    mask = copy(case.roi) # Also exclude native detected image support, not only fit tube.
    ax, ay = first(case.lobes).axis_x, first(case.lobes).axis_y
    for l in case.lobes, iy in eachindex(case.ys), ix in eachindex(case.xs)
        dx, dy = case.xs[ix]-l.x_nm, case.ys[iy]-l.y_nm
        # Conservatively include both halves of split-width geometry, if present.
        rx = s["footprint_sigma"]*l.sigma_parallel_nm*max(sqrt(l.skew_ratio),1/sqrt(l.skew_ratio)) + s["background_guard_nm"]
        ry = s["footprint_sigma"]*l.sigma_perp_nm + s["background_guard_nm"]
        ((dx*ax+dy*ay)/rx)^2 + ((-dx*ay+dy*ax)/ry)^2 <= 1 && (mask[iy,ix]=true)
    end
    # Off-footprint in BOTH acquisition coordinates at every declared shift.
    # Conservative square dilation also protects against lag search uncertainty.
    r = maximum(abs(s[k]) for k in ("lag_x_min_px","lag_x_max_px","lag_y_min_px","lag_y_max_px"))
    return STMSXMIO._dilate_mask(mask, r)
end

function saved_prediction(case)
    l0 = first(case.lobes)
    pred = [l0.baseline+l0.tilt_x*x+l0.tilt_y*y for y in case.ys, x in case.xs]
    case.profile in ("gaussian", "split") || error("Unsupported saved peak profile")
    for l in case.lobes, iy in eachindex(case.ys), ix in eachindex(case.xs)
        dx, dy = case.xs[ix]-l.x_nm, case.ys[iy]-l.y_nm
        t, u = dx*l.axis_x+dy*l.axis_y, -dx*l.axis_y+dy*l.axis_x
        spar = case.profile == "split" ? l.sigma_parallel_nm*(t<0 ? 1/sqrt(l.skew_ratio) : sqrt(l.skew_ratio)) : l.sigma_parallel_nm
        pred[iy,ix] += l.amplitude*exp(-0.5*((t/spar)^2+(u/l.sigma_perp_nm)^2))
    end
    return pred
end

function diagnose(case, s)
    f, b = case.forward, case.backward
    reg = registration_scan(f, b, case.fit_support, s)
    before = pair_metrics(sample_pairs(f,b,reg.support,0,0)...)
    after = pair_metrics(sample_pairs(f,b,reg.support,reg.dx,reg.dy)...)
    dx_nm, dy_nm = case.xs[2]-case.xs[1], case.ys[2]-case.ys[1]
    row_lags = NamedTuple[]
    for y in axes(f,1)
        row_mask = falses(size(reg.support))
        row_mask[y,:] .= reg.support[y,:]
        count(row_mask) == 0 && continue
        rr = scan_fixed(f,b,row_mask,s["lag_x_min_px"]:s["lag_x_max_px"],reg.dy:reg.dy,s;
                        min_pixels=s["min_row_pixels"])
        push!(row_lags,(row=y,y_nm=case.ys[y],pairs=rr.n,dx_px=rr.dx,dy_px=rr.dy,
                       correlation=rr.best_correlation,gap=rr.correlation_gap,
                       status=rr.status,usable=rr.usable))
    end
    qualified = Float64[r.dx_px for r in row_lags if r.usable]
    row_median = isempty(qualified) ? NaN : median(qualified)
    row_mad = isempty(qualified) ? NaN : 1.4826*median(abs.(qualified.-row_median))
    local_rows = NamedTuple[]
    for l in case.lobes
        local_mask = copy(reg.support)
        for I in findall(local_mask)
            (case.xs[I[2]]-l.x_nm)^2+(case.ys[I[1]]-l.y_nm)^2 <= s["local_patch_radius_nm"]^2 || (local_mask[I]=false)
        end
        for (condition,dx,dy) in (("native",0,0),("best_grid_shift_diagnostic_only",reg.dx,reg.dy))
            m = pair_metrics(sample_pairs(f,b,local_mask,dx,dy)...)
            status = count(local_mask) >= s["min_patch_pixels"] ? m.status : "insufficient_patch_support"
            push!(local_rows,merge((lobe=l.lobe,condition=condition,dx_px=dx,dy_px=dy,
                                  registration_usable=reg.usable),m,(status=status,)))
        end
    end
    footprint = geometry_footprint(case,s)
    background_mask = .!footprint .& isfinite.(f) .& isfinite.(b)
    nf = background_diagnostic(f,background_mask,s)
    nb = background_diagnostic(b,background_mask,s)
    noise_rows = NamedTuple[]
    for (view,noise) in (("fwd",nf),("bwd",nb)), row in noise.acf
        push!(noise_rows,merge((source="off_footprint_image",view=view),row))
    end
    # Genuine fitted residual diagnostics are kept separate; never used to supply
    # a background scale or an ACF-based block length when background is absent.
    pred = saved_prediction(case)
    residual = case.fused .- pred
    rm = case.fit_support .& case.fused_valid .& isfinite.(residual)
    residual_values = residual[rm]
    resid_mad = isempty(residual_values) ? NaN : 1.4826*median(abs.(residual_values.-median(residual_values)))
    for row in lag_acf(residual,rm,s)
        push!(noise_rows,merge((source="fitted_residual_not_noise",view="native_fused"),row))
    end
    lx = max(nf.length_x.pixels,nb.length_x.pixels)
    ly = max(nf.length_y.pixels,nb.length_y.pixels)
    blocks = block_feasibility(reg.support,lx,ly,s)
    bg_pair_support = common_support(f,b,background_mask,[0,reg.dx],[0,reg.dy])
    bg_pair = pair_metrics(sample_pairs(f,b,bg_pair_support,reg.dx,reg.dy)...)
    if count(bg_pair_support) < s["min_background_pixels"] ||
       count(bg_pair_support)/length(f) < s["min_background_fraction"]
        bg_pair = merge(bg_pair, (status="insufficient_background_crossview_support",))
    end
    # Finite sensitivity suggestions only. A boundary/ambiguous maximum is not a
    # measured correction and yields no proposed correction grid.
    perturb_x = Float64[]; perturb_y = Float64[]
    if reg.usable
        delta = isfinite(row_mad) ? max(1,ceil(Int,row_mad)) : 1
        perturb_x = unique(Float64[clamp(reg.dx+d,s["lag_x_min_px"],s["lag_x_max_px"])*dx_nm for d in (-delta,0,delta)])
        perturb_y = unique(Float64[clamp(reg.dy+d,s["lag_y_min_px"],s["lag_y_max_px"])*dy_nm for d in (-1,0,1)])
    end
    reasons = String["No stationary joint noise/covariance model validated; no calibrated likelihood/probabilities"]
    !reg.usable && push!(reasons,"Registration unresolved: "*reg.status)
    !blocks.feasible && push!(reasons,"Spatial validation "*blocks.status)
    auxiliary_rows=auxiliary_report(case,reg,s)
    return (registration=reg,before=before,after=after,row_lags=row_lags,local_rows=local_rows,
            auxiliary_rows=auxiliary_rows,
            noise_rows=noise_rows,noise_fwd=nf,noise_bwd=nb,background_crossview=bg_pair,
            background_pixels=count(background_mask),footprint_pixels=count(footprint),
            residual_mad=resid_mad,blocks=blocks,pixel_x_nm=dx_nm,pixel_y_nm=dy_nm,
            qualified_rows=length(qualified),row_median_dx_px=row_median,row_mad_dx_px=row_mad,
            perturb_dx_nm=perturb_x,perturb_dy_nm=perturb_y,blocked_reasons=reasons)
end

function write_tsv(path, rows, fields)
    open(path,"w") do io
        println(io,join(string.(fields),'\t'))
        for row in rows
            vals = [getproperty(row,k) for k in fields]
            println(io,join([v isa AbstractFloat && !isfinite(v) ? "NA" : replace(string(v),'\t'=>' ','\n'=>' ') for v in vals],'\t'))
        end
    end
end

plain(x::NamedTuple) = Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractDict) = Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractVector) = plain.(x)
plain(x) = x

function write_diagnostic(outdir, case, d, s; inputs::AbstractDict)
    (ispath(outdir) || islink(outdir)) && error("Exclusive diagnostic output directory already exists: $outdir")
    mkpath(dirname(abspath(outdir)))
    mkdir(outdir) # exclusive even if another invocation raced us
    reg=d.registration
    write_tsv(joinpath(outdir,"registration_candidates.tsv"),reg.candidates,
        (:dx_px,:dy_px,:n,:status,:correlation,:gain,:offset,:direct_nrmse,:positive_affine_nrmse,:difference_rms,:difference_median))
    write_tsv(joinpath(outdir,"row_lags.tsv"),d.row_lags,
        (:row,:y_nm,:pairs,:dx_px,:dy_px,:correlation,:gap,:status,:usable))
    write_tsv(joinpath(outdir,"local_views.tsv"),d.local_rows,
        (:lobe,:condition,:dx_px,:dy_px,:registration_usable,:n,:status,:correlation,:gain,:offset,:direct_nrmse,:positive_affine_nrmse,:difference_rms,:difference_median))
    write_tsv(joinpath(outdir,"noise_acf.tsv"),d.noise_rows,(:source,:view,:axis,:lag_px,:pairs,:correlation,:status))
    write_tsv(joinpath(outdir,"auxiliary_views.tsv"),d.auxiliary_rows,
        (:channel,:comparison,:condition,:source_unit,:target_unit,:dx_px,:dy_px,:z_support_pixels,
         :observed_pairs,:status,:correlation,:linear_r2,:positive_affine_nrmse,:gain,:offset,
         :analysis_pixels,:finite_pixels,:raw_min,:raw_max,:raw_std))
    summary=(file=case.file,status="diagnostic_only_joint_likelihood_blocked",selected_n=case.selected.n,
        registration_status=reg.status,registration_usable=reg.usable,pixels=reg.n,
        dx_px=reg.dx,dy_px=reg.dy,dx_nm=reg.dx*d.pixel_x_nm,dy_nm=reg.dy*d.pixel_y_nm,
        native_correlation=d.before.correlation,shifted_correlation=d.after.correlation,
        native_nrmse=d.before.positive_affine_nrmse,shifted_nrmse=d.after.positive_affine_nrmse,
        background_pixels=d.background_pixels,background_fwd_mad=d.noise_fwd.mad_sigma,
        background_bwd_mad=d.noise_bwd.mad_sigma,fitted_residual_mad_not_noise=d.residual_mad,
        block_status=d.blocks.status,usable_blocks=d.blocks.usable_blocks)
    write_tsv(joinpath(outdir,"summary.tsv"),[summary],keys(summary))
    meta=Dict{String,Any}(
        "schema"=>"acquisition_noise_diagnostic_v1", "julia_version"=>string(VERSION),
        "file"=>case.file, "unit"=>case.unit, "selected_context"=>plain(case.selected),
        "settings"=>plain(s),"support"=>plain(case.support_meta),
        "effective_support_noise_k"=>case.effective_support_noise_k,
        "effective_support_padding_nm"=>case.effective_support_padding_nm,
        "image_rows"=>size(case.fit_support,1),"image_columns"=>size(case.fit_support,2),
        "fit_support_pixels"=>count(case.fit_support),"constant_registration_pixels"=>reg.n,
        "fit_support_sha256_column_major_bytes"=>bytes2hex(sha256(UInt8.(vec(case.fit_support)))),
        "constant_registration_support_sha256_column_major_bytes"=>bytes2hex(sha256(UInt8.(vec(reg.support)))),
        "footprint_exclusion_pixels"=>d.footprint_pixels,
        "imputed_fwd_excluded"=>case.imputed_fwd,"imputed_bwd_excluded"=>case.imputed_bwd,
        "registration"=>plain(Base.structdiff(reg,(candidates=reg.candidates,support=reg.support))),
        "pixel_x_nm"=>d.pixel_x_nm,"pixel_y_nm"=>d.pixel_y_nm,
        "before"=>plain(d.before),"after"=>plain(d.after),
        "noise_fwd"=>plain(Base.structdiff(d.noise_fwd,(acf=d.noise_fwd.acf,))),
        "noise_bwd"=>plain(Base.structdiff(d.noise_bwd,(acf=d.noise_bwd.acf,))),
        "background_crossview"=>plain(d.background_crossview), "blocks"=>plain(d.blocks),
        "qualified_row_count"=>d.qualified_rows,"row_median_dx_px"=>d.row_median_dx_px,
        "row_mad_dx_px"=>d.row_mad_dx_px,"finite_registration_perturbations_dx_nm"=>d.perturb_dx_nm,
        "finite_registration_perturbations_dy_nm"=>d.perturb_dy_nm,
        "blocked_reasons"=>d.blocked_reasons,
        "auxiliary_diagnostic"=>Dict(
            "channel"=>s["auxiliary_channel"],"level"=>"direct_fixed_Z_support_and_lag_only",
            "nuisance_projection"=>"not_performed",
            "forward"=>plain(Base.structdiff(case.auxiliary.forward,(image=case.auxiliary.forward.image,))),
            "backward"=>plain(Base.structdiff(case.auxiliary.backward,(image=case.auxiliary.backward.image,))),
            "limit"=>"Coupled acquisition structure, not independent chemical evidence; no fit weights or likelihood changed"),
        "conventions"=>[
            "Native STMSXMIO.read_sxm backward x reversal already applied; never repeat it",
            "Shift sign: compare forward[y,x] with backward[y+dy,x+dx]; integer analysis-grid pixels after stride",
            "Native preprocess_channel unsmoothed flattened views used symmetrically for acquisition/noise",
            "Native fused target retained only for selected support replay and saved fitted residual; fwd unsmoothed, bwd smoothed",
            "Saved axis/origin and actual geometry fixed; adaptive support chosen from original refined_policy, not N",
            "Positive signed-correlation objective on constant finite support, no label/class agreement",
            "Off-footprint estimates can include substrate/scan structure; no assumption that views are independent",
            "First |ACF| threshold crossing is descriptive; censoring, missing pairs or observed resurgence block block-size suggestions",
            "Disjoint occupied rectangles are not independent observations or effective sample size",
            "Row MAD and grid perturbations are sensitivity ranges, not confidence intervals",
            "Auxiliary channel uses exactly the fixed Z support/lag; missing or nonfinite pairs are flagged, not filled or silently dropped",
            "Auxiliary raw scale is full-image finite analysis-grid dynamic range/std before flattening, not a noise estimate",
            "Z-auxiliary r squared is descriptive one-predictor linear variance explained, not calibrated evidence or nuisance-removed contrast",
            "No N selection, class refitting, production modification, likelihood or probability calibration"])
    meta["inputs"]=Dict(string(k)=>Dict("path"=>abspath(v),"sha256"=>bytes2hex(sha256(read(v)))) for (k,v) in inputs)
    open(joinpath(outdir,"metadata.toml"),"w") do io
        TOML.print(io,meta;sorted=true)
    end
    return outdir
end

end # module

module AcquisitionRegistration

using TOML, Statistics
export load_settings, scan, common_support

const FIELDS = Dict(
    "model" => ["method","reference_frame"],
    "selection" => ["bands","min_band_pixels","min_band_rows","min_row_pixels",
        "min_correlation","min_distant_peak_gap","peak_neighborhood_nm",
        "max_band_disagreement_nm","unresolved_policy"],
    "preprocessing" => ["max_lag_nm","max_lag_width_fraction","missing_policy"])

function load_settings(path)
    cfg = TOML.parsefile(path)
    Set(keys(cfg)) == Set(keys(FIELDS)) || error("Unexpected registration sections")
    for (section,fields) in FIELDS
        Set(keys(cfg[section])) == Set(fields) || error("Incomplete/unknown registration settings")
    end
    s = merge(values(cfg)...)
    s["method"] == "signed_row_centered_integer_x_translation" || error("Unknown registration method")
    s["reference_frame"] == "forward" || error("Only forward reference frame is supported")
    s["unresolved_policy"] == "zero_shift_keep_all_files" || error("Unknown unresolved policy")
    s["missing_policy"] == "restore_raw_mask_before_smoothing" || error("Unknown missing policy")
    for k in ("bands","min_band_pixels","min_band_rows","min_row_pixels")
        s[k] isa Integer && !(s[k] isa Bool) && s[k] >= 2 || error("Invalid integer $k")
    end
    for k in ("min_correlation","min_distant_peak_gap","peak_neighborhood_nm",
              "max_band_disagreement_nm","max_lag_nm","max_lag_width_fraction")
        s[k] isa Real && !(s[k] isa Bool) && isfinite(s[k]) && s[k] > 0 || error("Invalid $k")
    end
    s["min_correlation"] <= 1 && s["min_distant_peak_gap"] < 1 || error("Invalid correlation gates")
    s["max_lag_width_fraction"] < 0.5 || error("Lag window leaves no fixed interior")
    s["max_lag_nm"] > s["peak_neighborhood_nm"] || error("Lag window smaller than peak neighborhood")
    return s
end

"Forward pixels observed at every candidate backward x; never variable overlaps."
function common_support(a,b,limit,min_row_pixels)
    size(a) == size(b) || throw(DimensionMismatch("View shapes differ"))
    limit >= 0 || error("Negative search limit")
    ny,nx = size(a); mask = falses(ny,nx)
    for y in 1:ny
        badprefix = vcat(0,cumsum(.!isfinite.(@view b[y,:])))
        for x in (1+limit):(nx-limit)
            mask[y,x] = isfinite(a[y,x]) && badprefix[x+limit+1] == badprefix[x-limit]
        end
        count(@view mask[y,:]) >= min_row_pixels || (mask[y,:] .= false)
    end
    return mask
end

function peak_result(lags, correlations, npixels, nrows, s, tolerance)
    enough = npixels >= s["min_band_pixels"] && nrows >= s["min_band_rows"]
    valid = findall(isfinite,correlations)
    bestidx = isempty(valid) ? nothing : first(sort(valid;by=i->(-correlations[i],abs(lags[i]),lags[i])))
    best = bestidx === nothing ? 0 : lags[bestidx]
    corr = bestidx === nothing ? NaN : correlations[bestidx]
    distant = [correlations[i] for i in valid if abs(lags[i]-best) > tolerance]
    gap = isempty(distant) ? NaN : corr-maximum(distant)
    boundary = bestidx !== nothing && best in (first(lags),last(lags))
    reasons = String[]
    enough || push!(reasons,"insufficient_support")
    isfinite(corr) || push!(reasons,"constant_signal")
    isfinite(corr) && corr < s["min_correlation"] && push!(reasons,"weak_or_negative_correlation")
    (!isfinite(gap) || gap < s["min_distant_peak_gap"]) && push!(reasons,"ambiguous_peak")
    boundary && push!(reasons,"boundary_peak")
    return (best_dx_px=best,correlation=corr,distant_peak_gap=gap,boundary=boundary,
        pixels=npixels,rows=nrows,usable=isempty(reasons),status=isempty(reasons) ? "identified" : join(reasons,";"))
end

"""One fixed signed translation, accepted only if all contiguous bands identify it.
Row-centering removes row offsets from the comparison, without fitting gain or
changing the production flattening. Bands are stability checks, NOT independent
replicates or a calibrated confidence interval. No molecular/template/label input.
"""
function scan(a::Matrix{Float64},b::Matrix{Float64},pixel_nm::Real,s)
    size(a) == size(b) || throw(DimensionMismatch("View shapes differ"))
    isfinite(pixel_nm) && pixel_nm > 0 || error("Invalid pixel scale")
    ny,nx = size(a); nb = s["bands"]
    limit = min(floor(Int,s["max_lag_nm"]/pixel_nm),floor(Int,s["max_lag_width_fraction"]*nx))
    tolerance = max(1,ceil(Int,s["peak_neighborhood_nm"]/pixel_nm))
    agreement = max(1,ceil(Int,s["max_band_disagreement_nm"]/pixel_nm))
    lags = collect(-limit:limit)
    support = common_support(a,b,limit,s["min_row_pixels"])
    # Sufficient statistics pool within-row centered contrast, not row baselines.
    aa = zeros(nb); npix=zeros(Int,nb); nrows=zeros(Int,nb)
    ab = zeros(nb,length(lags)); bb = zeros(nb,length(lags))
    for y in 1:ny
        xs = findall(@view support[y,:]); isempty(xs) && continue
        band = min(nb,1+((y-1)*nb)÷ny)
        av = a[y,xs]; av .-= mean(av)
        aa[band] += sum(abs2,av); npix[band] += length(xs); nrows[band] += 1
        for (i,lag) in enumerate(lags)
            bv = b[y,xs .+ lag]; bv .-= mean(bv)
            bb[band,i] += sum(abs2,bv)
            ab[band,i] += sum(av .* bv)
        end
    end
    correlations = fill(NaN,nb+1,length(lags))
    for i in eachindex(lags), band in 0:nb
        va = band == 0 ? sum(aa) : aa[band]
        vb = band == 0 ? sum(@view bb[:,i]) : bb[band,i]
        cab = band == 0 ? sum(@view ab[:,i]) : ab[band,i]
        va > 0 && vb > 0 && (correlations[band+1,i] = clamp(cab/sqrt(va*vb),-1.,1.))
    end
    peaks = [peak_result(lags,vec(correlations[band+1,:]),band==0 ? sum(npix) : npix[band],
        band==0 ? sum(nrows) : nrows[band],s,tolerance) for band in 0:nb]
    global_peak = first(peaks)
    reasons = String[]
    limit > tolerance+1 || push!(reasons,"insufficient_lag_span")
    global_peak.usable || push!(reasons,"global:"*global_peak.status)
    for band in 1:nb
        p = peaks[band+1]
        p.usable || push!(reasons,"band$band:"*p.status)
        abs(p.best_dx_px-global_peak.best_dx_px) <= agreement || push!(reasons,"band$(band)_disagreement")
    end
    accepted = isempty(reasons)
    return (accepted=accepted,applied_dx_px=accepted ? global_peak.best_dx_px : 0,
        status=accepted ? "identified" : join(reasons,";"),limit_px=limit,tolerance_px=tolerance,
        agreement_px=agreement,lags=lags,correlations=correlations,peaks=peaks,support=support)
end

end

"""
Opt-in four-view preprocessing comparison. No count fit, prediction, residual
model, Current stage, calibrated noise, or accepted registration is produced.
Frozen native geometry/support remain affected by the original preprocessing.
"""
module MaskedPreprocessingComparison

using Statistics, LinearAlgebra, TOML, SHA, Serialization
using STMSXMIO, GaussianFit2D
include(joinpath(@__DIR__, "masked_robust_preprocessing.jl"))
include(joinpath(@__DIR__, "acquisition_noise_diagnostics.jl"))
const MRP = MaskedRobustPreprocessing
const AD = AcquisitionNoiseDiagnostics
const METHODS = ("native_reference", "finite_only_ols", "guarded_ols", "guarded_huber")
const DIRECTIONS = ("fwd", "bwd")

export METHODS, metadata_context, load_real_context, build_views, compare_views,
       write_comparison, require_new_output, plain, require_julia

require_julia() = v"1.13.0" <= VERSION < v"1.14.0" || error("This comparison requires Julia 1.13")

function require_new_output(path)
    isempty(strip(path)) && error("Output path is empty")
    out = abspath(path)
    (ispath(out) || islink(out)) && error("Output already exists: $out")
    parent = dirname(out)
    while true
        islink(parent) && error("Output ancestor is a symlink: $parent")
        ispath(parent) && !isdir(parent) && error("Output ancestor is not a directory: $parent")
        parent == dirname(parent) && break
        parent = dirname(parent)
    end
    return out
end

function comparison_settings(s)
    AD.settings_from_dict(s)
    Tuple(s[k] for k in ("lag_x_min_px", "lag_x_max_px", "lag_y_min_px", "lag_y_max_px")) ==
        (-8, 8, -2, 2) || error("Comparison retains the original ±8 x / ±2 y window (85 lags)")
    lowercase(s["channel"]) == "z" || error("Comparison requires the native Z channel")
    return s
end

"""Metadata only: never reads raw file bytes or creates output, even for fake SXM."""
function metadata_context(raw, geometry, physical, summary, acquisition, settings)
    for path in (raw, geometry, physical, summary, acquisition, settings)
        isfile(path) || error("Input is not a file: $path")
    end
    endswith(lowercase(raw), ".sxm") || error("Raw input must name an SXM file")
    cfg = TOML.parsefile(physical)
    all(k -> haskey(cfg, k), ("model", "preprocessing")) || error("Incomplete physical config")
    model, pre = cfg["model"], cfg["preprocessing"]
    for k in ("stride", "flatten", "smooth_radius_px")
        haskey(pre, k) || error("Frozen preprocessing.$k is required")
    end
    lowercase(pre["flatten"]) == "plane+rows" || error("Native reference must use original plane+rows")
    pre["stride"] isa Integer && !(pre["stride"] isa Bool) && pre["stride"] > 0 || error("Invalid stride")
    pre["smooth_radius_px"] isa Integer && !(pre["smooth_radius_px"] isa Bool) && pre["smooth_radius_px"] >= 0 || error("Invalid native support smoothing")
    for k in ("fit_width_nm", "support_noise_k", "support_padding_nm", "support_min_length_nm", "support_baseline_quantile")
        haskey(model, k) || error("Frozen model.$k is required")
        model[k] isa Real && !(model[k] isa Bool) && isfinite(model[k]) || error("Invalid model.$k")
    end
    lines = readlines(summary)
    isempty(lines) && error("Empty selected summary")
    header = split(first(lines), '\t') .|> strip
    any(k -> k in header, ("truth", "expected_N", "target_N", "sequence", "control_sequence", "label", "predicted")) &&
        error("Selected summary must be label-free")
    contexts = AD.SavedSupport._read_selected_context(String(summary);
        selection_policy=String(get(model, "selection_policy", "")))
    file = basename(raw)
    haskey(contexts, file) || error("No selected summary for $file")
    selected = contexts[file]
    lobes = AD.read_geometry(geometry, file, selected.n)
    profile = String(get(model, "peak_profile", "gaussian"))
    profile in ("gaussian", "split") || error("Unsupported saved peak profile")
    profile == "gaussian" && any(l -> l.skew_ratio != 1.0, lobes) && error("Saved split widths contradict Gaussian config")
    # Construction only: no optimizer, ROI, preprocessing or raw reader is called.
    _, ccfg, _ = AD.SavedSupport._configs(model, pre, ""; selected_context=selected)
    s = comparison_settings(AD.load_settings(acquisition))
    mask_config = TOML.parsefile(settings)
    options = MRP.read_options(mask_config)
    return (; file, selected, lobes, physical_config=cfg, s, options, mask_config,
        effective_support_noise_k=ccfg.support_noise_k,
        effective_support_padding_nm=ccfg.support_padding_nm)
end

"""Real execution only. Second native read restores scaled raw Z, never a second flip."""
function load_real_context(raw, geometry, physical, summary, context)
    # Reuse the exact frozen support replay with UNCHANGED settings. load_case
    # incidentally preprocesses auxiliary Current arrays, then discards them.
    # No Current evidence, whole acquisition diagnose or residual stage is used.
    loaded = AD.load_case(raw, geometry, physical, summary, context.s)
    case = Base.structdiff(loaded, (auxiliary=loaded.auxiliary,))
    case.unit == "nm" || error("Native Z height must be in nm")
    case.selected == context.selected || error("Selected context changed during load")
    case.lobes == context.lobes || error("Saved geometry changed during load")
    case.physical_config == context.physical_config || error("Physical config changed during load")
    img = STMSXMIO.read_sxm(String(raw))
    stride = context.physical_config["preprocessing"]["stride"]
    xs, ys = STMSXMIO._coordinate_vectors(img; stride=stride)
    xs == case.xs && ys == case.ys || error("Second native read grid mismatch")
    fch, bch = (AD.strict_channel(img, context.s["channel"], d) for d in DIRECTIONS)
    raw_views = map((fch, bch), (case.forward, case.backward)) do ch, native
        scale, unit = STMSXMIO._value_scale(ch.unit)
        unit == "nm" || error("Raw Z scale must be nm")
        full = ch.data .* scale
        z = full[1:stride:end, 1:stride:end]
        size(z) == size(native) || error("Raw/native size mismatch")
        isfinite.(z) == isfinite.(native) || error("Raw/native observed-mask mismatch")
        return (z=z, observed=BitMatrix(isfinite.(z)), fullimage_pixels=length(full),
            fullimage_observed=count(isfinite, full), fullimage_rows=size(full,1), fullimage_columns=size(full,2))
    end
    rawdata = (forward=raw_views[1], backward=raw_views[2])
    footprint = AD.geometry_footprint(case, context.s)
    return (; case, rawdata, footprint)
end

function native_view(z, raw)
    size(z) == size(raw.z) || throw(DimensionMismatch("Native/raw image sizes differ"))
    isfinite.(z) == raw.observed == isfinite.(raw.z) || error("Native observed support differs from raw")
    ny = size(z, 1)
    return (corrected=copy(z), valid=BitMatrix(isfinite.(z)), observed=copy(raw.observed),
        background_mask="not_exported_native_uses_imputed_full_grid", plane=fill(NaN,size(z)),
        background=fill(NaN,size(z)), row_offsets=fill(NaN,ny),
        row_status=fill("native_row_diagnostics_not_exported",ny), row_background_counts=fill(NaN,ny),
        status="REFERENCE", reason="native_diagnostics_not_exported", converged="unknown", iterations=NaN,
        coefficients=(intercept_nm=NaN, slope_x=NaN, slope_y=NaN), robust_scale=NaN,
        coefficient_change=NaN, weight_change=NaN, scale_floor_active="unknown", design_rank=NaN,
        loss="native", initial_objective_nm2=NaN, objective_nm2=NaN, stationarity_inf_nm=NaN)
end

"""Exactly eight method/direction records. Numerical failures are not removed."""
function build_views(case, rawdata, footprint, options)
    size(footprint) == size(case.fit_support) == size(case.forward) == size(case.backward) ||
        throw(DimensionMismatch("Frozen image/mask sizes differ"))
    views = Dict{String,Any}()
    for (direction, raw, native) in zip(DIRECTIONS, (rawdata.forward, rawdata.backward), (case.forward, case.backward))
        views["native_reference/" * direction] = native_view(native, raw)
        for method in METHODS[2:end]
            bgmask = method == "finite_only_ols" ? trues(size(footprint)) : .!footprint
            loss = method == "guarded_huber" ? :huber : :ols
            views[method * "/" * direction] = MRP.masked_background(case.xs, case.ys, raw.z, bgmask, options; loss=loss)
        end
    end
    return views
end

viewkey(method, direction) = method * "/" * direction
image(views, method, direction) = views[viewkey(method, direction)].corrected
maskhash(mask) = bytes2hex(sha256(UInt8.(vec(mask))))
safe_fraction(n, d) = d > 0 ? n / d : NaN

function comparison_scopes(case, views, footprint, s)
    Set(keys(views)) == Set(viewkey(m,d) for m in METHODS for d in DIRECTIONS) || error("Need all eight method/view records")
    dxs, dys = s["lag_x_min_px"]:s["lag_x_max_px"], s["lag_y_min_px"]:s["lag_y_max_px"]
    own = Dict(m => AD.common_support(image(views,m,"fwd"), image(views,m,"bwd"),
        case.fit_support, dxs, dys) for m in METHODS)
    all_common = reduce((a,b) -> a .& b, [own[m] for m in METHODS])
    definitions = [(id="own/"*m, kind="own", methods=[m], support=own[m]) for m in METHODS]
    push!(definitions, (id="all_common", kind="all_common", methods=collect(METHODS), support=all_common))
    for m in METHODS[2:end]
        push!(definitions, (id="native_pair/"*m, kind="native_pair", methods=[METHODS[1],m],
            support=own[METHODS[1]] .& own[m]))
    end
    scopes = NamedTuple[]
    for scope in definitions
        # One common zero-lag background anchor for EVERY view in this scope.
        # It is independent of fit support (which is inside the footprint),
        # never adapted per lag and never chosen using any truth/signal array.
        anchor = .!footprint
        for m in scope.methods, direction in DIRECTIONS
            v = views[viewkey(m,direction)]
            size(v.corrected) == size(case.fit_support) || error("View grid mismatch")
            v.valid == isfinite.(v.corrected) || error("Corrected availability mismatch")
            any(v.valid .& .!v.observed) && error("Imputed pixels cannot become observations")
            anchor .&= v.valid .& v.observed
        end
        push!(scopes, merge(scope, (background_support=anchor, support_sha256=maskhash(scope.support),
                                   background_sha256=maskhash(anchor))))
    end
    return (; own, all_common, scopes)
end

function unavailable_metrics(n, status)
    return (n=n, status=status, correlation=NaN, gain=NaN, offset=NaN,
        direct_nrmse=NaN, positive_affine_nrmse=NaN, difference_rms=NaN, difference_median=NaN)
end

function metrics_on(a, b, mask, dx, dy, minimum_pixels)
    n = count(mask)
    n >= minimum_pixels || return unavailable_metrics(n, "insufficient_support")
    return AD.pair_metrics(AD.sample_pairs(a,b,mask,dx,dy)...)
end

function fixed_gauge(a, b, scope, s)
    n = count(scope.background_support)
    adequate = n >= s["min_background_pixels"] && n/length(a) >= s["min_background_fraction"]
    if !adequate || !any(scope.support)
        return (status=!adequate ? "insufficient_shared_observed_background" : "no_comparison_support",
                anchor_pixels=n, forward_level_nm=NaN, backward_level_nm=NaN)
    end
    return (status="observable_shared_background_median", anchor_pixels=n,
        forward_level_nm=median(a[scope.background_support]), backward_level_nm=median(b[scope.background_support]))
end

function level_rms(a, b, mask, dx, dy, gauge, minimum_pixels)
    count(mask) >= minimum_pixels || return (level_status="insufficient_support", level_aligned_rms_nm=NaN)
    gauge.status == "observable_shared_background_median" || return (level_status=gauge.status, level_aligned_rms_nm=NaN)
    va,vb = AD.sample_pairs(a,b,mask,dx,dy)
    return (level_status="observable_constant_only", level_aligned_rms_nm=
        sqrt(mean(abs2, (va .- gauge.forward_level_nm) .- (vb .- gauge.backward_level_nm))))
end

function background_rows(views, masks, s)
    rows = NamedTuple[]
    for scope in masks.scopes, method in scope.methods, direction in DIRECTIONS
        mask = scope.background_support
        n = count(mask)
        fraction = n/length(mask)
        adequate = n >= s["min_background_pixels"] && fraction >= s["min_background_fraction"]
        z, native = image(views,method,direction), image(views,METHODS[1],direction)
        all(isfinite,native[mask]) || error("Native lacks observations on fixed background support")
        v, reference = z[mask], native[mask]
        level = adequate ? median(v) : NaN
        native_level = adequate ? median(reference) : NaN
        push!(rows, (scope=scope.id, scope_kind=scope.kind, method=method, direction=direction,
            status=adequate ? "descriptive_postfit_background_not_noise_calibration" : "insufficient_shared_observed_background",
            background_sha256=scope.background_sha256, n=n, fraction=fraction,
            comparison_support_pixels=count(scope.support), mean_nm=adequate ? mean(v) : NaN,
            mad_sigma_nm=adequate ? 1.4826*median(abs.(v .- median(v))) : NaN,
            standard_deviation_nm=adequate ? std(v) : NaN,
            level_nm=level, native_level_nm=native_level,
            change_level_aligned_rms_vs_native_nm=adequate ?
                sqrt(mean(abs2,(v .- level) .- (reference .- native_level))) : NaN))
    end
    return rows
end

function local_support(case, support, lobe, s)
    mask = copy(support)
    for I in findall(mask)
        (case.xs[I[2]]-lobe.x_nm)^2+(case.ys[I[1]]-lobe.y_nm)^2 <= s["local_patch_radius_nm"]^2 || (mask[I]=false)
    end
    return mask
end

function view_tables(case, rawdata, views, footprint)
    view_rows, row_rows = NamedTuple[], NamedTuple[]
    for m in METHODS, (direction,raw) in zip(DIRECTIONS,(rawdata.forward,rawdata.backward))
        v = views[viewkey(m,direction)]
        n = length(v.corrected)
        push!(view_rows, (method=m, direction=direction, status=string(v.status), reason=string(v.reason),
            converged=v.converged, iterations=v.iterations, observed_pixels=count(v.observed),
            valid_pixels=count(v.valid), analysis_pixels=n, raw_missing_pixels=count(.!v.observed),
            corrected_unavailable_observed=count(v.observed .& .!v.valid),
            observed_fraction=count(v.observed)/n, valid_fraction=count(v.valid)/n,
            fullimage_pixels=raw.fullimage_pixels, fullimage_observed=raw.fullimage_observed,
            fullimage_missing_pixels=raw.fullimage_pixels-raw.fullimage_observed,
            fullimage_observed_fraction=raw.fullimage_observed/raw.fullimage_pixels,
            supplied_background_pixels=m == "native_reference" ? NaN : count(v.background_mask),
            plane_intercept_nm=v.coefficients.intercept_nm, plane_slope_x=v.coefficients.slope_x,
            plane_slope_y=v.coefficients.slope_y, robust_scale_nm=v.robust_scale,
            coefficient_change=v.coefficient_change, weight_change=v.weight_change,
            scale_floor_active=v.scale_floor_active, design_rank=v.design_rank,
            initial_objective_nm2=v.initial_objective_nm2, objective_nm2=v.objective_nm2,
            stationarity_inf_nm=v.stationarity_inf_nm))
        for y in eachindex(case.ys)
            push!(row_rows, (method=m, direction=direction, row=y, y_nm=case.ys[y],
                status=string(v.row_status[y]), observed_pixels=count(v.observed[y,:]),
                valid_pixels=count(v.valid[y,:]), background_pixels=v.row_background_counts[y],
                guarded_observed_pixels=count(v.observed[y,:] .& .!footprint[y,:]), row_offset_nm=v.row_offsets[y]))
        end
    end
    return (; view_rows, row_rows)
end

"""
Each method has whole-lag-valid own support, then all-four and native-pair
intersections. Empty/insufficient scopes still produce all 85 lag rows.
Grid-best shifts are diagnostics only, also when the native helper marks usable.
"""
function compare_views(case, rawdata, views, footprint, s)
    comparison_settings(s)
    masks = comparison_scopes(case, views, footprint, s)
    tables = view_tables(case, rawdata, views, footprint)
    summaries, candidates, local_rows = NamedTuple[], NamedTuple[], NamedTuple[]
    dxs, dys = s["lag_x_min_px"]:s["lag_x_max_px"], s["lag_y_min_px"]:s["lag_y_max_px"]
    for scope in masks.scopes, method in scope.methods
        f,b = image(views,method,"fwd"), image(views,method,"bwd")
        reg = AD.registration_scan(f,b,scope.support,s)
        reg.support == scope.support || error("Already-common registration support was not idempotent")
        has_best = isfinite(reg.best_correlation)
        gauge = fixed_gauge(f,b,scope,s)
        base = (scope=scope.id, scope_kind=scope.kind, method=method,
            forward_status=string(views[viewkey(method,"fwd")].status),
            backward_status=string(views[viewkey(method,"bwd")].status),
            fit_support_pixels=count(case.fit_support), own_support_pixels=count(masks.own[method]),
            all_common_pixels=count(masks.all_common), scope_pixels=count(scope.support),
            scope_fraction_fit=safe_fraction(count(scope.support),count(case.fit_support)),
            scope_fraction_own=safe_fraction(count(scope.support),count(masks.own[method])),
            support_sha256=scope.support_sha256, background_sha256=scope.background_sha256)
        flags = (registration_status=reg.status, registration_usable_diagnostic=reg.usable,
            registration_accepted=false, boundary=reg.boundary, ambiguous=reg.ambiguous,
            insufficient_support=reg.n<s["min_registration_pixels"],
            correlation_gap=reg.correlation_gap, near_maxima=reg.near_maxima)
        levels = (anchor_status=gauge.status, anchor_pixels=gauge.anchor_pixels,
            forward_level_nm=gauge.forward_level_nm, backward_level_nm=gauge.backward_level_nm)
        # The native helper omits insufficient candidates. Keep the DECLARED
        # grid here without inventing scores, correlations, or a best shift.
        lagrows = isempty(reg.candidates) ?
            [merge((dx_px=dx,dy_px=dy), unavailable_metrics(reg.n,"insufficient_support")) for dy in dys for dx in dxs] : reg.candidates
        length(lagrows) == length(dxs)*length(dys) || error("Incomplete lag grid")
        for r in lagrows
            lr = level_rms(f,b,scope.support,r.dx_px,r.dy_px,gauge,s["min_registration_pixels"])
            push!(candidates, merge(base, flags, levels, r, lr))
        end
        for condition in ("zero_lag", "grid_best_diagnostic_only")
            known = condition == "zero_lag" || has_best
            dx,dy = condition == "zero_lag" ? (0,0) : (reg.dx,reg.dy)
            m = known ? metrics_on(f,b,scope.support,dx,dy,s["min_registration_pixels"]) :
                unavailable_metrics(reg.n,"no_grid_maximum")
            lr = known ? level_rms(f,b,scope.support,dx,dy,gauge,s["min_registration_pixels"]) :
                (level_status="no_grid_maximum",level_aligned_rms_nm=NaN)
            shift = (condition=condition, dx_px=known ? dx : NaN, dy_px=known ? dy : NaN)
            push!(summaries, merge(base, flags, levels, shift, m, lr))
            for l in case.lobes
                lm = local_support(case,scope.support,l,s)
                metric = known ? metrics_on(f,b,lm,dx,dy,s["min_patch_pixels"]) :
                    unavailable_metrics(count(lm),"no_grid_maximum")
                push!(local_rows, merge(base, flags, shift,
                    (lobe=l.lobe, x_nm=l.x_nm, y_nm=l.y_nm, local_support_sha256=maskhash(lm)), metric))
            end
        end
    end
    return merge(masks,tables,(; summaries,candidates,local_rows, background_rows=background_rows(views,masks,s)))
end

# File format contains only built-in scalar/string/array/dictionary values.
# NaN is preserved; unavailable symbols/nothing/missing become explicit strings.
plain(x::NamedTuple) = Dict{String,Any}(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractDict) = Dict{String,Any}(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractArray) = map(plain, x)
plain(x::Symbol) = string(x)
plain(::Nothing) = "unavailable"
plain(::Missing) = "NA"
plain(x::Union{Bool,Integer,AbstractFloat,AbstractString}) = x
plain(x) = error("Unsupported serialized type $(typeof(x)); save explicit plain fields")

function snapshot(case, rawdata, views, footprint, d, context)
    frozen = (file=case.file, xs=collect(case.xs), ys=collect(case.ys), unit=case.unit,
        fit_support=case.fit_support, roi=case.roi, footprint=footprint,
        lobes=case.lobes, selected=case.selected, support_meta=case.support_meta,
        effective_support_noise_k=case.effective_support_noise_k,
        effective_support_padding_nm=case.effective_support_padding_nm)
    options = Dict{String,Any}(string(k)=>getfield(context.options,k) for k in fieldnames(typeof(context.options)))
    return plain(Dict("schema"=>"masked_preprocessing_comparison_arrays_v1",
        "julia_version"=>string(VERSION), "frozen"=>frozen, "raw_scaled_analysis_nm"=>rawdata,
        "views"=>views, "own_supports"=>d.own, "all_common_support"=>d.all_common,
        "scopes"=>d.scopes, "acquisition_settings"=>context.s,
        "masked_options"=>options, "physical_config"=>context.physical_config))
end

function source_paths()
    root = normpath(joinpath(@__DIR__, "..", ".."))
    paths = [joinpath(root,p) for p in ("Project.toml", "Manifest.toml",
        "test/lib/masked_preprocessing_comparison.jl", "test/diagnose_masked_preprocessing_real.jl",
        "test/lib/masked_robust_preprocessing.jl", "test/lib/acquisition_noise_diagnostics.jl",
        "test/extract_lobe_features.jl")]
    for package in ("STMFitCore", "STMSXMIO", "GaussianFit1D", "GaussianFit2D", "STMMolecularFit")
        for (dir, _, files) in walkdir(joinpath(root,"packages",package*".jl","src")), file in files
            endswith(file,".jl") && push!(paths,joinpath(dir,file))
        end
    end
    return Dict(relpath(p,root)=>p for p in paths)
end

filehash(path) = open(io -> bytes2hex(sha256(io)), path, "r")
hash_records(paths) = Dict(string(k)=>Dict("path"=>abspath(p),"sha256"=>filehash(p)) for (k,p) in paths)

function write_comparison(outdir, case, rawdata, views, footprint, d, context; inputs)
    require_julia()
    out = require_new_output(outdir)
    # Convert before creating output; never serialize SXMImage, options or models.
    arrays = snapshot(case,rawdata,views,footprint,d,context)
    input_hashes, source_hashes = hash_records(inputs), hash_records(source_paths())
    mkpath(dirname(out)); mkdir(out)
    outputs = Dict("view_status.tsv"=>d.view_rows, "row_status.tsv"=>d.row_rows,
        "comparison_summary.tsv"=>d.summaries, "lag_candidates.tsv"=>d.candidates,
        "local_comparisons.tsv"=>d.local_rows, "background_summary.tsv"=>d.background_rows)
    for (name, rows) in outputs
        isempty(rows) && error("Expected retained rows for $name")
        AD.write_tsv(joinpath(out,name), rows, keys(first(rows)))
    end
    serialize(joinpath(out,"arrays.jls"),arrays)
    artifacts = Dict(name=>filehash(joinpath(out,name)) for name in [collect(keys(outputs)); "arrays.jls"])
    meta = Dict{String,Any}(
        "schema"=>"masked_preprocessing_comparison_v1", "julia_version"=>string(VERSION),
        "file"=>case.file, "height_unit"=>"nm", "methods"=>collect(METHODS),
        "selected_context"=>plain(case.selected), "support"=>plain(case.support_meta),
        "physical_config"=>plain(context.physical_config), "acquisition_settings"=>plain(context.s),
        "masked_options"=>arrays["masked_options"], "inputs"=>input_hashes, "sources"=>source_hashes,
        "artifacts_sha256"=>artifacts, "fit_support_pixels"=>count(case.fit_support),
        "all_common_pixels"=>count(d.all_common), "footprint_pixels"=>count(footprint),
        "fit_support_sha256"=>maskhash(case.fit_support), "roi_sha256"=>maskhash(case.roi),
        "footprint_sha256"=>maskhash(footprint), "all_common_sha256"=>maskhash(d.all_common),
        "effective_support_noise_k"=>case.effective_support_noise_k,
        "effective_support_padding_nm"=>case.effective_support_padding_nm,
        "lag_candidates_per_scope_method"=>85, "view_records"=>length(d.view_rows),
        "row_records"=>length(d.row_rows), "summary_records"=>length(d.summaries),
        "lag_records"=>length(d.candidates), "local_records"=>length(d.local_rows),
        "background_records"=>length(d.background_rows),
        "array_format"=>"Julia 1.13 Serialization; plain built-in scalars, arrays and string-key dictionaries only; z[y,x]; no custom objects. Load trusted files only.",
        "raw_array_scope"=>"Raw scaled nm arrays after original stride; full-image counts before stride also saved. No filling, flattening or second backward flip in these arrays.",
        "conventions"=>[
            "Frozen native ROI, fit support and saved geometry inherited from old preprocessing, including imputation effects; not an independently validated mask.",
            "The unchanged geometry_footprint excludes native ROI plus guarded saved ellipses, then conservatively square-dilates by max lag radius 8, including y beyond ±2.",
            "Native reference is actual preprocess_channel unsmoothed plane+rows output; its raw missing mask is restored. Unknown native estimator diagnostics are NA/unknown, never zero.",
            "Own support requires forward observations and backward observations at every one of 85 lags. All-common and native-pair supports intersect own masks before any score.",
            "A failed method empties all-common; blocked rows and unaffected pairwise evidence remain separately identified. Coverage is never silently restored.",
            "Forward[y,x] versus backward[y+dy,x+dx], after original stride; no additional smoothing, interpolation, local registration or sign flip.",
            "All grid-best shifts are diagnostic-only and never accepted. Signed correlations, ambiguity, boundary and insufficiency remain explicit.",
            "For level-aligned physical RMS, one median per view uses exactly the shared finite observed guarded background mask of that scope at zero lag. Constants remain fixed for every lag.",
            "Correlation is offset-invariant. Positive-affine metrics are separately named diagnostics, not a correction or gain warp.",
            "Background mean/MAD/std and level-aligned change versus native are descriptive post-fit structure, not calibrated noise; shared support and original minimum background counts/fraction apply.",
            "Frozen load_case incidentally preprocesses auxiliary Current arrays, then discards them. No Current evidence is computed or exported by this comparison. No whole acquisition diagnose, saved-fit residual, ACF/block sizing, n_eff, calibrated noise or joint likelihood; no winner selection, count or chemical fit."])
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,meta;sorted=true)
    end
    return out
end

end # module

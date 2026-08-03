#!/usr/bin/env julia

include(joinpath(@__DIR__, "diagnose_joint_proxy_view_preprocessing.jl"))

const ROW_DRIFT_FIELDS = [:file, :flatten_mode, :status, :n, :row_index, :y_nm,
    :x_start_index, :x_end_index, :base_corr, :best_lag_px, :best_lag_nm,
    :best_corr, :corr_gain, :n_pairs, :slope, :intercept, :affine_nrmse]

Base.@kwdef struct RowDriftOptions
    config::String
    data_dir::String
    file::String
    out::String
    flatten_mode::String
    max_lag_px::Int
    min_pairs::Int
    seed::Int = 20260721
end

function parse_row_drift_options(args)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        flag = String(args[i])
        flag in ("--config", "--data-dir", "--file", "--out", "--flatten-mode",
            "--max-lag-px", "--min-pairs", "--seed") || error("unknown argument: $flag")
        i < length(args) || error("$flag requires a value")
        values[flag] = String(args[i + 1])
        i += 2
    end
    for flag in ("--config", "--data-dir", "--file", "--out", "--flatten-mode",
        "--max-lag-px", "--min-pairs")
        haskey(values, flag) || error("$flag is required")
    end
    file = values["--file"]
    basename(file) == file && endswith(lowercase(file), ".sxm") ||
        error("--file must be a bare .sxm basename")
    flatten_mode = lowercase(strip(values["--flatten-mode"]))
    flatten_mode in VALID_FLATTEN_MODES ||
        error("--flatten-mode must use: $(join(VALID_FLATTEN_MODES, ','))")
    max_lag_px = parse(Int, values["--max-lag-px"])
    min_pairs = parse(Int, values["--min-pairs"])
    max_lag_px >= 0 || error("--max-lag-px must be non-negative")
    min_pairs >= 3 || error("--min-pairs must be at least 3")
    return RowDriftOptions(
        config=values["--config"], data_dir=values["--data-dir"], file=file,
        out=values["--out"], flatten_mode=flatten_mode, max_lag_px=max_lag_px,
        min_pairs=min_pairs, seed=parse(Int, get(values, "--seed", "20260721")),
    )
end

function _lagged_row_pairs(fwd::AbstractVector{<:Real}, bwd::AbstractVector{<:Real},
    lag::Int, max_lag_px::Int)
    length(fwd) == length(bwd) || error("fwd/bwd row length mismatch")
    max_lag_px < length(fwd) ÷ 2 || error("max lag leaves no common row support")
    abs(lag) <= max_lag_px || error("lag exceeds common support bound")
    paired_fwd = Float64[]
    paired_bwd = Float64[]
    for i in (firstindex(fwd) + max_lag_px):(lastindex(fwd) - max_lag_px)
        j = i + lag
        push!(paired_fwd, Float64(fwd[i]))
        push!(paired_bwd, Float64(bwd[j]))
    end
    return paired_fwd, paired_bwd
end

function row_lag_metrics(fwd::AbstractVector{<:Real}, bwd::AbstractVector{<:Real},
    max_lag_px::Int, min_pairs::Int)
    max_lag_px >= 0 || error("max lag must be non-negative")
    min_pairs >= 3 || error("minimum pairs must be at least 3")
    length(fwd) == length(bwd) || error("fwd/bwd row length mismatch")
    if 2max_lag_px + min_pairs > length(fwd)
        base = _affine_patch_metrics(fwd, bwd)
        return (; status="insufficient_pairs", base_corr=base.corr,
            best_lag_px=0, best_corr=NaN, corr_gain=NaN, n_pairs=base.n_pairs,
            slope=NaN, intercept=NaN, affine_nrmse=NaN)
    end
    base_fwd, base_bwd = _lagged_row_pairs(fwd, bwd, 0, max_lag_px)
    base = _affine_patch_metrics(base_fwd, base_bwd)
    best_lag = 0
    best = nothing
    best_abs = -Inf
    for lag in -max_lag_px:max_lag_px
        paired_fwd, paired_bwd = _lagged_row_pairs(fwd, bwd, lag, max_lag_px)
        candidate = _affine_patch_metrics(paired_fwd, paired_bwd)
        candidate.n_pairs >= min_pairs && isfinite(candidate.corr) || continue
        candidate_abs = abs(candidate.corr)
        if candidate_abs > best_abs + 1e-12 ||
            (abs(candidate_abs - best_abs) <= 1e-12 && abs(lag) < abs(best_lag))
            best_lag = lag
            best = candidate
            best_abs = candidate_abs
        end
    end
    best === nothing && return (; status="insufficient_pairs", base_corr=base.corr,
        best_lag_px=0, best_corr=NaN, corr_gain=NaN, n_pairs=base.n_pairs,
        slope=NaN, intercept=NaN, affine_nrmse=NaN)
    corr_gain = isfinite(base.corr) ? abs(best.corr) - abs(base.corr) : NaN
    return (; status="ok", base_corr=base.corr, best_lag_px=best_lag,
        best_corr=best.corr, corr_gain, n_pairs=best.n_pairs, slope=best.slope,
        intercept=best.intercept, affine_nrmse=best.affine_nrmse)
end

function roi_bounds(mask::AbstractMatrix{Bool})
    indices = findall(mask)
    isempty(indices) && error("fixed ROI mask is empty")
    rows = first.(Tuple.(indices))
    cols = last.(Tuple.(indices))
    return minimum(rows):maximum(rows), minimum(cols):maximum(cols)
end

function compute_row_drift_rows(case, file::String, flatten_mode::String,
    max_lag_px::Int, min_pairs::Int)
    mask = case.context.mask
    mask === nothing && error("fixed-geometry fit produced no ROI mask")
    size(mask) == (length(case.context.ys), length(case.context.xs)) ||
        error("fixed ROI mask/grid size mismatch")
    row_range, col_range = roi_bounds(mask)
    length(col_range) >= 2max_lag_px + min_pairs ||
        error("ROI is too narrow for requested lag and minimum pairs")

    pattern = pattern_with_flatten(case.base_pattern, flatten_mode)
    views = Inference._build_views(case.image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    by_label = Dict(view.label => view for view in views if view.present)
    haskey(by_label, "fwd") && haskey(by_label, "bwd") || error("both fwd and bwd views are required")
    fwd, bwd = by_label["fwd"], by_label["bwd"]
    fwd.xs == bwd.xs && fwd.ys == bwd.ys || error("fwd/bwd preprocessing grids differ")
    x_step_nm = length(fwd.xs) > 1 ? median(diff(fwd.xs)) : NaN

    rows = NamedTuple[]
    for iy in row_range
        metrics = row_lag_metrics(view(fwd.z, iy, col_range), view(bwd.z, iy, col_range),
            max_lag_px, min_pairs)
        push!(rows, (; file, flatten_mode,
            metrics.status, n=case.candidate.n, row_index=iy, y_nm=fwd.ys[iy],
            x_start_index=first(col_range), x_end_index=last(col_range),
            metrics.base_corr, metrics.best_lag_px,
            best_lag_nm=metrics.best_lag_px * x_step_nm, metrics.best_corr,
            metrics.corr_gain, metrics.n_pairs, metrics.slope, metrics.intercept,
            metrics.affine_nrmse))
    end
    return rows
end

function run_row_drift_diagnostic(options::RowDriftOptions)
    case = fixed_geometry_case(options.config, options.data_dir, options.file,
        dirname(options.out), options.seed)
    rows = compute_row_drift_rows(case, options.file, options.flatten_mode,
        options.max_lag_px, options.min_pairs)
    write_tsv(options.out, rows; fields=ROW_DRIFT_FIELDS)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_row_drift_diagnostic(parse_row_drift_options(ARGS))
end

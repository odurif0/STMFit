#!/usr/bin/env julia

# Extract aligned per-lobe STM patches from an existing lobe-features TSV.
# The lobe order/labels are not used. Output is a wide TSV with normalized raw
# and residual patch pixels for downstream PCA/template diagnostics.

using Printf
using Statistics
using GaussianFit2D: read_sxm, preprocess_channel, PatternConfig, get_channel

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

include(joinpath(@__DIR__, "lib", "patch_preprocessing.jl"))
using .PatchPreprocessing: PreprocessingSettings, load_patch_preprocessing,
    load_patch_residual_filter, patch_residual
include(joinpath(@__DIR__, "lib", "patch_acquisition.jl"))
using .PatchAcquisition

const DEFAULT_FEATURES = "results/unit_separability/lobe_features_selectedN_primary.tsv"
const DEFAULT_OUT = "results/unit_separability/lobe_patches_selectedN_primary.tsv"

struct Options
    features::String
    out_tsv::String
    data_dir::String
    half_nm::Float64
    half_u_nm::Float64
    step_nm::Float64
    preprocessing::PreprocessingSettings
    residual_filter::String
    acquisition_shifts::Union{Nothing,String}
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    out_tsv = DEFAULT_OUT
    data_dir = get(ENV, "STMFIT_DATA_DIR", "")
    config_path::Union{Nothing,String} = nothing
    assignment_path::Union{Nothing,String} = nothing
    acquisition_path::Union{Nothing,String} = nothing
    half_nm = 0.32
    half_u_nm = 0.32
    step_nm = 0.08
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"; features = args[i+1]; i += 2
        elseif startswith(arg, "--features="); features = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--data-dir"; data_dir = args[i+1]; i += 2
        elseif startswith(arg, "--data-dir="); data_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--config"
            config_path === nothing || error("Duplicate --config option")
            i < length(args) && !startswith(args[i+1], "--") || error("Missing value for --config")
            config_path = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config_path === nothing || error("Duplicate --config option")
            config_path = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--acquisition-shifts"
            acquisition_path === nothing || error("Duplicate --acquisition-shifts")
            i < length(args) && !startswith(args[i+1], "--") || error("Missing acquisition table")
            acquisition_path = args[i+1]; i += 2
        elseif arg == "--assignment-config"
            assignment_path === nothing || error("Duplicate --assignment-config option")
            i < length(args) && !startswith(args[i+1], "--") || error("Missing value for --assignment-config")
            assignment_path = args[i+1]; i += 2
        elseif startswith(arg, "--assignment-config=")
            assignment_path === nothing || error("Duplicate --assignment-config option")
            assignment_path = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--half-nm"; half_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-nm="); half_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--half-u-nm"; half_u_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--half-u-nm="); half_u_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: STMFIT_DATA_DIR=/data julia --project=. test/extract_lobe_patches.jl [options]

            Options:
              --features PATH   Lobe feature TSV [$(DEFAULT_FEATURES)]
              --out PATH        Output patch TSV [$(DEFAULT_OUT)]
              --data-dir PATH   SXM data directory [\$STMFIT_DATA_DIR]
              --config PATH     Optional count TOML: all three [preprocessing] fields
                                (omitted: stride=1, flatten=plane+rows, smooth_radius_px=1)
              --acquisition-shifts PATH  Opt-in backward shift table; forward only restores raw masks
              --assignment-config PATH  TOML with [preprocessing] patch_residual_filter
                                        (omitted: legacy smooth_data_only)
              --half-nm FLOAT   Patch half-size along the chain t [0.32]
              --half-u-nm FLOAT Patch half-size transverse u [0.32]
              --step-nm FLOAT   Patch grid spacing [0.08]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isempty(data_dir) && error("No data directory: set STMFIT_DATA_DIR or pass --data-dir")
    isdir(data_dir) || error("Data directory not found: $data_dir")
    isfile(features) || error("Features TSV not found: $features")
    preprocessing = load_patch_preprocessing(config_path)
    residual_filter = load_patch_residual_filter(assignment_path)
    return Options(features, out_tsv, data_dir, half_nm, half_u_nm, step_nm, preprocessing, residual_filter, acquisition_path)
end

function _eval_peak(x, y, cx, cy, ax, ay, A, spar, sperp, skew_ratio)
    dt = (x - cx) * ax + (y - cy) * ay
    du = (x - cx) * (-ay) + (y - cy) * ax
    r = max(skew_ratio, eps(Float64))
    sigma_t = dt < 0 ? spar / sqrt(r) : spar * sqrt(r)
    return A * exp(-0.5 * (dt/sigma_t)^2 - 0.5 * (du/sperp)^2)
end

function _interp(xs, ys, z, x, y)
    ix = searchsortedlast(xs, x)
    iy = searchsortedlast(ys, y)
    if ix < 1 || iy < 1 || ix >= length(xs) || iy >= length(ys)
        return NaN
    end
    x1, x2 = xs[ix], xs[ix+1]
    y1, y2 = ys[iy], ys[iy+1]
    tx = (x - x1) / (x2 - x1)
    ty = (y - y1) / (y2 - y1)
    z11 = z[iy, ix]; z21 = z[iy, ix+1]
    z12 = z[iy+1, ix]; z22 = z[iy+1, ix+1]
    return (1-tx)*(1-ty)*z11 + tx*(1-ty)*z21 + (1-tx)*ty*z12 + tx*ty*z22
end

function _normalize_patch(vals)
    good = filter(isfinite, vals)
    isempty(good) && return fill(NaN, length(vals))
    μ = median(good)
    σ = std(good)
    σ = σ > 0 ? σ : 1.0
    return [isfinite(v) ? (v - μ) / σ : NaN for v in vals]
end

function main(args=ARGS)
    opt = _parse_cli(args)
    _, rows = _read_tsv(opt.features)
    isempty(rows) && error("No feature rows in $(opt.features)")
    by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in rows
        push!(get!(by_file, basename(row["file"]), Dict{String,String}[]), row)
    end
    shifts = read_shifts(opt.acquisition_shifts,keys(by_file))
    coords = collect(-opt.half_nm:opt.step_nm:opt.half_nm)
    coords_u = collect(-opt.half_u_nm:opt.step_nm:opt.half_u_nm)
    pix_names = [@sprintf("%03d", i) for i in 1:(length(coords) * length(coords_u))]

    mkpath(dirname(opt.out_tsv))
    open(opt.out_tsv, "w") do io
        println(io, join(vcat(["file", "lobe", "t_nm", "u_nm", "amplitude"],
                              ["raw_p$(p)" for p in pix_names],
                              ["res_p$(p)" for p in pix_names]), '\t'))
        for (idx, file) in enumerate(sort(collect(keys(by_file))))
            sxm_path = joinpath(opt.data_dir, file)
            isfile(sxm_path) || (@warn "SXM not found, skipping" file; continue)
            rs = sort(by_file[file], by=r -> parse(Int, r["lobe"]))
            try
                img = read_sxm(sxm_path)
                ch = get_channel(img, "Z"; direction="fwd")
                pcfg = PatternConfig(filepath=sxm_path, channel="Z", direction="fwd",
                    stride=opt.preprocessing.stride, flatten=opt.preprocessing.flatten,
                    smooth_radius_px=opt.preprocessing.smooth_radius_px,
                    output_dir=dirname(opt.out_tsv), no_plot=true)
                xs, ys, raw, z, z_smooth, scaled_unit, noise = preprocess_channel(img, ch, pcfg)
                if shifts !== nothing
                    require_direction(ch,"fwd")
                    z,z_smooth = observed_shift(z,ch,opt.preprocessing,0)
                end
                nx, ny = length(xs), length(ys)
                ax = _parse_f(rs[1]["axis_x"]); ay = _parse_f(rs[1]["axis_y"])
                b0 = haskey(rs[1], "baseline") ? _parse_f(rs[1]["baseline"]) : 0.0
                bx = haskey(rs[1], "tilt_x") ? _parse_f(rs[1]["tilt_x"]) : 0.0
                by = haskey(rs[1], "tilt_y") ? _parse_f(rs[1]["tilt_y"]) : 0.0
                model = zeros(ny, nx)
                for iy in 1:ny, ix in 1:nx
                    model[iy, ix] = b0 + bx * xs[ix] + by * ys[iy]
                end
                for row in rs
                    A = _parse_f(row["amplitude"]); cx = _parse_f(row["x_nm"]); cy = _parse_f(row["y_nm"])
                    spar = _parse_f(row["sigma_parallel_nm"]); sperp = _parse_f(row["sigma_perp_nm"])
                    skew = haskey(row, "skew_ratio") ? _parse_f(row["skew_ratio"]) : 1.0
                    isfinite(skew) || (skew = 1.0)
                    for iy in 1:ny, ix in 1:nx
                        model[iy, ix] += _eval_peak(xs[ix], ys[iy], cx, cy, ax, ay, A, spar, sperp, skew)
                    end
                end
                residual = patch_residual(z, z_smooth, model, opt.preprocessing, opt.residual_filter)

                for row in rs
                    cx = _parse_f(row["x_nm"]); cy = _parse_f(row["y_nm"])
                    raw_vals = Float64[]; res_vals = Float64[]
                    for u in coords_u, t in coords
                        x = cx + t * ax + u * (-ay)
                        y = cy + t * ay + u * ax
                        push!(raw_vals, _interp(xs, ys, z_smooth, x, y))
                        push!(res_vals, _interp(xs, ys, residual, x, y))
                    end
                    raw_norm = _normalize_patch(raw_vals)
                    res_norm = _normalize_patch(res_vals)
                    vals = vcat([file, row["lobe"], row["t_nm"], row["u_nm"], row["amplitude"]],
                                [isfinite(v) ? @sprintf("%.7g", v) : "NA" for v in raw_norm],
                                [isfinite(v) ? @sprintf("%.7g", v) : "NA" for v in res_norm])
                    println(io, join(vals, '\t'))
                end
                @printf("[%d/%d] %-24s patches=%d grid=%dx%d\n", idx, length(by_file), file, length(rs), length(coords), length(coords_u))
            catch e
                @warn "Failed" file reason=sprint(showerror, e)
            end
        end
    end
    println("Wrote: ", opt.out_tsv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

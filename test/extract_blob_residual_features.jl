#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# extract_blob_residual_features.jl — Phase 1b: extract non-Gaussian features
# from the fit residual around each lobe center.
#
# The Gaussian model is symmetric per lobe. If GlcNAc and GlcN differ in
# intra-blob asymmetry (e.g. the C2 acetyl group creates a "shoulder" offset
# along the chain axis), this signal is lost in the symmetric Gaussian but
# visible in the residual (data - model).
#
# This script reads the per-lobe features TSV from extract_lobe_features.jl
# (which contains the fitted Gaussian parameters + chain axis direction) and
# the original SXM files, then:
#   1. Preprocesses each SXM (same config as the batch).
#   2. Reconstructs the Gaussian model (sum of 2D Gaussians, no baseline).
#   3. Computes the residual = data_smoothed - gaussian_sum.
#   4. For each lobe center, extracts:
#      - skewness_axial:     third standardized moment of the residual profile
#                            along the chain axis (±2σ∥ window)
#      - kurtosis_axial:     fourth standardized moment (excess kurtosis)
#      - shoulder_left:      mean residual at t = t_k - δ
#      - shoulder_right:     mean residual at t = t_k + δ
#      - shoulder_asymmetry: (shoulder_right - shoulder_left) / (|R| + |L|)
#      - lr_asymmetry:       (|R| - |L|) / (|R| + |L|) for |residual| sums
#      - residual_peak_snr:  max |residual| / noise in the local window
#
# The shoulder offset δ (default 0.15 nm) is guided by the pyranose ring
# geometry: the C2 substituent is ~1/3 of the monomer spacing from the ring
# center. This is a PHYSICAL PRIOR (not a label).
#
# Examples:
#   STMFIT_DATA_DIR=/path/to/data julia --project=. \
#       test/extract_blob_residual_features.jl \
#       --features results/unit_separability/lobe_features.tsv \
#       --out results/unit_separability/residual_features.tsv
# ──────────────────────────────────────────────────────────────────────────────

using Printf
using Statistics
using StatsBase
using DelimitedFiles
using GaussianFit2D: read_sxm, preprocess_channel, PatternConfig, get_channel
using GaussianFit2D: ChainSweepConfig

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_FEATURES = "results/unit_separability/lobe_features.tsv"
const DEFAULT_OUT = "results/unit_separability/residual_features.tsv"
const DEFAULT_DELTA_NM = 0.15
const DEFAULT_WINDOW_SIGMA = 2.0

struct ResOptions
    features::String
    out_tsv::String
    data_dir::String
    delta_nm::Float64
    window_sigma::Float64
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    out_tsv = DEFAULT_OUT
    data_dir = get(ENV, "STMFIT_DATA_DIR", "")
    delta_nm = DEFAULT_DELTA_NM
    window_sigma = DEFAULT_WINDOW_SIGMA

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"
            features = args[i+1]; i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--data-dir"
            data_dir = args[i+1]; i += 2
        elseif startswith(arg, "--data-dir=")
            data_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--delta-nm"
            delta_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--delta-nm=")
            delta_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--window-sigma"
            window_sigma = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--window-sigma=")
            window_sigma = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: STMFIT_DATA_DIR=/data julia --project=. test/extract_blob_residual_features.jl [options]

            Options:
              --features PATH       Lobe features TSV from extract_lobe_features.jl [$(DEFAULT_FEATURES)]
              --out PATH            Output TSV [$(DEFAULT_OUT)]
              --data-dir PATH       SXM data directory [\$STMFIT_DATA_DIR]
              --delta-nm FLOAT      Shoulder offset along axis (nm) [$(DEFAULT_DELTA_NM)]
              --window-sigma FLOAT  Window half-width in σ∥ units [$(DEFAULT_WINDOW_SIGMA)]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isempty(data_dir) && error("No data directory: set STMFIT_DATA_DIR or pass --data-dir")
    isdir(data_dir) || error("Data directory not found: $data_dir")
    return ResOptions(features, out_tsv, data_dir, delta_nm, window_sigma)
end

_basename_file(s::AbstractString) = basename(strip(s))

# Evaluate a single 2D Gaussian at (x, y) given center, axis, widths.
function _eval_gaussian(x, y, cx, cy, ax, ay, A, spar, sperp)
    dt = (x - cx) * ax + (y - cy) * ay
    du = (x - cx) * (-ay) + (y - cy) * ax
    return A * exp(-0.5 * (dt/spar)^2 - 0.5 * (du/sperp)^2)
end

function _weighted_moments(values::Vector{Float64}, weights::Vector{Float64})
    sw = sum(weights)
    sw <= 0 && return NaN, NaN
    mu = sum(values .* weights) / sw
    d = values .- mu
    m2 = sum(weights .* d.^2) / sw
    m3 = sum(weights .* d.^3) / sw
    m4 = sum(weights .* d.^4) / sw
    m2 <= 0 && return NaN, NaN
    return m3 / m2^1.5, m4 / m2^2 - 3.0
end

function main()
    opt = _parse_cli(ARGS)
    mkpath(dirname(opt.out_tsv))

    # Read features TSV (group by file)
    _, feat_rows = _read_tsv(opt.features)
    isempty(feat_rows) && error("No feature rows in $(opt.features)")

    by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in feat_rows
        file = _basename_file(row["file"])
        push!(get!(by_file, file, Dict{String,String}[]), row)
    end

    files = sort(collect(keys(by_file)))
    println("Processing ", length(files), " files")

    open(opt.out_tsv, "w") do io
        println(io, join([
            "file", "lobe", "amplitude", "t_nm", "u_nm", "sigma_parallel_nm", "sigma_perp_nm",
            "skewness_axial", "kurtosis_axial",
            "shoulder_left", "shoulder_right", "shoulder_asymmetry",
            "lr_asymmetry", "residual_peak_snr"
        ], '\t'))

        for (i_file, file) in enumerate(files)
            sxm_path = joinpath(opt.data_dir, file)
            if !isfile(sxm_path)
                @warn "SXM not found, skipping" file
                continue
            end

            rows = by_file[file]
            # Axis from first row (same for all lobes of a file)
            ax = _parse_f(rows[1]["axis_x"])
            ay = _parse_f(rows[1]["axis_y"])
            ox = _parse_f(rows[1]["origin_x_nm"])
            oy = _parse_f(rows[1]["origin_y_nm"])

            try
                img = read_sxm(sxm_path)
                ch = get_channel(img, "Z"; direction="fwd")
                pcfg = PatternConfig(filepath=sxm_path, channel="Z", direction="fwd",
                    stride=1, flatten="plane+rows", smooth_radius_px=1,
                    output_dir=dirname(opt.out_tsv), no_plot=true)
                xs, ys, raw, z, z_smooth, scaled_unit, noise = preprocess_channel(img, ch, pcfg)

                nx, ny = length(xs), length(ys)

                # Reconstruct full model if baseline/tilt are available from
                # extract_lobe_features.jl. Older TSVs do not have these
                # columns; in that case fall back to Gaussian-only residual.
                b0 = haskey(rows[1], "baseline") ? _parse_f(rows[1]["baseline"]) : 0.0
                bx = haskey(rows[1], "tilt_x") ? _parse_f(rows[1]["tilt_x"]) : 0.0
                by = haskey(rows[1], "tilt_y") ? _parse_f(rows[1]["tilt_y"]) : 0.0
                model = zeros(ny, nx)
                for iy in 1:ny, ix in 1:nx
                    model[iy, ix] = b0 + bx * xs[ix] + by * ys[iy]
                end
                for row in rows
                    A = _parse_f(row["amplitude"])
                    cx = _parse_f(row["x_nm"])
                    cy = _parse_f(row["y_nm"])
                    spar = _parse_f(row["sigma_parallel_nm"])
                    sperp = _parse_f(row["sigma_perp_nm"])
                    for iy in 1:ny, ix in 1:nx
                        model[iy, ix] += _eval_gaussian(xs[ix], ys[iy], cx, cy, ax, ay, A, spar, sperp)
                    end
                end

                residual = z_smooth .- model

                for (i, row) in enumerate(rows)
                    A = _parse_f(row["amplitude"])
                    spar = _parse_f(row["sigma_parallel_nm"])
                    sperp = _parse_f(row["sigma_perp_nm"])
                    cx = _parse_f(row["x_nm"])
                    cy = _parse_f(row["y_nm"])
                    t_k = _parse_f(row["t_nm"])
                    u_k = _parse_f(row["u_nm"])

                    wt = opt.window_sigma * spar
                    wu = 2.0 * sperp

                    # Extract local residual profile along axis
                    t_vals = Float64[]
                    r_vals = Float64[]
                    for iy in 1:ny, ix in 1:nx
                        x_px, y_px = xs[ix], ys[iy]
                        dt = (x_px - cx) * ax + (y_px - cy) * ay
                        du = (x_px - cx) * (-ay) + (y_px - cy) * ax
                        if abs(dt) <= wt && abs(du) <= wu
                            push!(t_vals, dt)
                            push!(r_vals, residual[iy, ix])
                        end
                    end

                    if length(t_vals) < 5
                        println(io, join([file, string(i), @sprintf("%.6f", A),
                                          @sprintf("%.6f", t_k), @sprintf("%.6f", u_k),
                                          @sprintf("%.6f", spar), @sprintf("%.6f", sperp),
                                          "NA", "NA", "NA", "NA", "NA", "NA", "NA"], '\t'))
                        continue
                    end

                    # Weighted skewness and kurtosis
                    weights = abs.(r_vals) .+ 1e-30
                    skew, kurt = _weighted_moments(t_vals, weights)

                    # Shoulder at ±δ
                    delta = opt.delta_nm
                    sw = spar * 0.5
                    sr = [r_vals[j] for j in 1:length(t_vals) if abs(t_vals[j] - delta) <= sw]
                    sl = [r_vals[j] for j in 1:length(t_vals) if abs(t_vals[j] + delta) <= sw]
                    shoulder_r = isempty(sr) ? NaN : mean(sr)
                    shoulder_l = isempty(sl) ? NaN : mean(sl)
                    denom = abs(shoulder_r) + abs(shoulder_l)
                    shoulder_asym = denom > 1e-30 ? (shoulder_r - shoulder_l) / denom : NaN

                    # L/R asymmetry
                    L = sum(abs.(r_vals[t_vals .< 0]))
                    R = sum(abs.(r_vals[t_vals .> 0]))
                    lr_asym = (L + R) > 1e-30 ? (R - L) / (R + L) : NaN

                    # Residual peak SNR
                    local_peak = isempty(r_vals) ? 0.0 : maximum(abs.(r_vals))
                    res_snr = noise > 1e-30 ? local_peak / noise : NaN

                    println(io, join([file, string(i), @sprintf("%.6f", A),
                                      @sprintf("%.6f", t_k), @sprintf("%.6f", u_k),
                                      @sprintf("%.6f", spar), @sprintf("%.6f", sperp),
                                      isnan(skew) ? "NA" : @sprintf("%.6f", skew),
                                      isnan(kurt) ? "NA" : @sprintf("%.6f", kurt),
                                      isnan(shoulder_l) ? "NA" : @sprintf("%.6f", shoulder_l),
                                      isnan(shoulder_r) ? "NA" : @sprintf("%.6f", shoulder_r),
                                      isnan(shoulder_asym) ? "NA" : @sprintf("%.6f", shoulder_asym),
                                      isnan(lr_asym) ? "NA" : @sprintf("%.6f", lr_asym),
                                      isnan(res_snr) ? "NA" : @sprintf("%.6f", res_snr)], '\t'))
                end

                @printf("[%d/%d] %-24s  lobes=%d  noise=%.4e\n", i_file, length(files), file, length(rows), noise)
            catch e
                @warn "Failed" file reason=sprint(showerror, e)
            end
        end
    end

    println("\nWrote: ", opt.out_tsv)
end

main()

#!/usr/bin/env julia

# Diagnostic-only: whole-ROI type comparison on a real 240817 scan using
# converged ±0.80 nm periodic DFT molds.
#
# This script fits a 2D chain on the scan, decodes the lobe geometry, crops
# the molecular ROI from the preprocessed forward image, assembles the
# converged DFT molds at decoded positions, and runs a whole-ROI nuisance
# fit + exhaustive sequence search. It uses NO benchmark labels, expected N,
# control sequence, or composition prior.
#
# Prerequisites:
#   /tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv
#   STMFIT_DATA_DIR pointing to the 240817 data
#
# Usage:
#   STMFIT_DATA_DIR=/path/to/data julia --project=. \
#       test/diagnose_joint_proxy_whole_roi_real.jl 240817_001.sxm

using Printf
using LinearAlgebra
using STMSXMIO
using GaussianFit2D

include(joinpath(@__DIR__, "lib", "joint_proxy", "proxy_registry.jl"))
include(joinpath(@__DIR__, "lib", "joint_proxy", "inference_cli.jl"))
include(joinpath(@__DIR__, "diagnose_joint_proxy_whole_roi.jl"))

const WR_REG = Main.JointProxyRegistry
const Inference = Main.JointProxyInferenceCLI

function load_converged_entry(path::String)
    isfile(path) || error("converged mold TSV not found: $path")
    templates = WR_REG.load_stm_template_tsv(path; npix=441)
    source = WR_REG.ProxySource("stm_dft_periodic080", "converged_diagnostic",
        path, WR_REG.sha256_file(path), 0.50, 1.0, true)
    return WR_REG.ProxyEntry(source, templates)
end

function _parse_args(args)
    isempty(args) || return String(args[1])
    error("usage: julia --project=. test/diagnose_joint_proxy_whole_roi_real.jl <file.sxm>")
end

function _chain_theta(axisctx)
    ax, ay = axisctx.axis
    return atan(ay, ax)
end

function _crop_roi(xs, ys, z, centers, margin_nm)
    tmin = minimum(c[1] for c in centers) - margin_nm
    tmax = maximum(c[1] for c in centers) + margin_nm
    umin = minimum(c[2] for c in centers) - margin_nm
    umax = maximum(c[2] for c in centers) + margin_nm
    ix = findall(t -> tmin <= t <= tmax, xs)
    iy = findall(u -> umin <= u <= umax, ys)
    isempty(ix) && error("ROI crop produced empty x range")
    isempty(iy) && error("ROI crop produced empty y range")
    return xs[ix], ys[iy], z[iy, ix], ix, iy
end

function main()
    filename = _parse_args(ARGS)
    data_dir = something(get(ENV, "STMFIT_DATA_DIR", nothing), ".")
    path = joinpath(data_dir, filename)
    isfile(path) || error("SXM file not found: $path")

    mold_path = "/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv"
    entry = load_converged_entry(mold_path)
    half_nm = 0.80
    step_nm = 0.08

    config_path = "config/joint_proxy_molds.toml"
    registry = WR_REG.load_registry(config_path)

    println("Real-ROI whole-ROI diagnostic")
    println("  file:  ", filename)
    println("  molds: ", mold_path)
    println()

    image = STMSXMIO.read_sxm(path)
    pattern = Inference._default_pattern_config(path, "/tmp/opencode/real_roi_diagnostic")
    chain = Inference._default_chain_config()

    results, best, context = redirect_stdout(devnull) do
        GaussianFit2D.chain_gaussian_sweep(image, pattern, chain)
    end

    views = Inference._build_views(image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    report = Main.JointProxyCandidateViews.extract_candidate_views(
        results, context, chain, views;
        patch_half_nm=registry.grid_half_nm, patch_step_nm=registry.grid_step_nm)

    isempty(report.candidates) && error("no valid candidates from fit")
    candidate = argmin(item -> item.joint_gcv, report.candidates)
    n = candidate.n
    println("  fit:   N_selected=", n, "  joint_gcv=", @sprintf("%.6f", candidate.joint_gcv))
    println("  lobes: ", length(candidate.lobes))

    fwd = only(v for v in views if v.label == "fwd" && v.present)
    theta = _chain_theta(context.axisctx)
    centers = [(lobe.x_nm, lobe.y_nm) for lobe in candidate.lobes]

    roi_margin = half_nm + 0.1
    rxs, rys, rz, ixx, iyy = _crop_roi(fwd.xs, fwd.ys, fwd.z, centers, roi_margin)
    println("  ROI:   ", length(rxs), "x", length(rys), " pixels")
    println("  theta: ", @sprintf("%.4f rad", theta))

    # Build the Gaussian backbone from the fitted model on the ROI grid
    model_full = reshape(GaussianFit2D._chain_model_values(
        vec(repeat(reshape(fwd.xs, 1, :), length(fwd.ys), 1)),
        vec(repeat(reshape(fwd.ys, :, 1), 1, length(fwd.xs))),
        candidate.source_gcv > 0 ? results[min(end, n)].params : results[1].params,
        n, context.axisctx, chain),
        length(fwd.ys), length(fwd.xs))
    backbone_roi = model_full[iyy, ixx]

    # Build backbone from decoded parameters directly (robust to which ChainModelResult)
    best_result = results[findfirst(r -> r.n == n && r.success && r.valid, results)]
    model_vec = GaussianFit2D._chain_model_values(
        vec(repeat(reshape(rxs, 1, :), length(rys), 1)),
        vec(repeat(reshape(rys, :, 1), 1, length(rxs))),
        best_result.params, n, context.axisctx, chain;
        amp_min=best_result.amp_min, amp_range=best_result.amp_range)
    backbone_roi = reshape(model_vec, length(rys), length(rxs))

    observed = rz

    # Assemble molds for each candidate sequence and find best
    best_seq = nothing
    best_fit = nothing
    println()
    println("Searching ", 2^n, " sequences (N=", n, ")...")
    for bits in 0:(2^n - 1)
        sequence = [Int((bits >> (lobe - 1)) & 1) for lobe in 1:n]
        mold = assemble_whole_roi(entry, sequence, centers, rxs, rys;
            half_nm=half_nm, step_nm=step_nm, theta=theta,
            direction=0, phase=0, mirror=0)
        fit = fit_whole_roi_nuisance(observed, backbone_roi, mold)
        if best_fit === nothing || fit.sse < best_fit.sse
            best_seq = sequence
            best_fit = fit
        end
    end

    # Also search over direction/phase/mirror states for the best sequence
    best_state = (0, 0, 0)
    for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
        for bits in 0:(2^n - 1)
            sequence = [Int((bits >> (lobe - 1)) & 1) for lobe in 1:n]
            mold = assemble_whole_roi(entry, sequence, centers, rxs, rys;
                half_nm=half_nm, step_nm=step_nm, theta=theta,
                direction=direction, phase=phase, mirror=mirror)
            fit = fit_whole_roi_nuisance(observed, backbone_roi, mold)
            if best_fit === nothing || fit.sse < best_fit.sse
                best_seq = sequence
                best_fit = fit
                best_state = (direction, phase, mirror)
            end
        end
    end

    # Compute null (backbone-only) for comparison
    good = isfinite.(observed) .& isfinite.(backbone_roi)
    null_design = ones(Float64, count(good), 2)
    null_design[:, 2] .= backbone_roi[good]
    null_coeff = null_design \ observed[good]
    null_coeff[2] = max(null_coeff[2], 0.0)
    null_residual = observed[good] .- null_design * null_coeff
    null_sse = sum(abs2, null_residual)

    # Compute relative SSE reduction
    sse_reduction = 1.0 - best_fit.sse / null_sse

    println()
    println("Results:")
    println("  best sequence:      ", best_seq)
    println("  direction/phase/mirror: ", best_state)
    println("  null SSE (backbone only):  ", @sprintf("%.6e", null_sse))
    println("  best SSE (backbone+mold):  ", @sprintf("%.6e", best_fit.sse))
    println("  SSE reduction:      ", @sprintf("%.4f%%", 100 * sse_reduction))
    println("  background:         ", @sprintf("%.6f", best_fit.background))
    println("  backbone coeff:     ", @sprintf("%.6f", best_fit.backbone_coefficient))
    println("  mold coeff:         ", @sprintf("%.6f", best_fit.mold_coefficient))
    println("  n_pixels:           ", best_fit.n_pixels)

    # Honest assessment
    println()
    if best_fit.mold_coefficient < 1e-8
        println("Decision: mold contribution is negligible; honest abstention.")
    elseif sse_reduction < 0.001
        println("Decision: SSE reduction < 0.1%; insufficient evidence for type assignment.")
    else
        println("Decision: SSE reduction = ", @sprintf("%.2f%%", 100 * sse_reduction),
                "; mold coefficient = ", @sprintf("%.4f", best_fit.mold_coefficient))
        println("  NOTE: this is a real-ROI comparison, not a validated chemical",
                " identity. The scan-domain gap (constant-height LDOS vs",
                " constant-current topography) remains unmodeled.")
    end

    # Per-sequence SSE summary for the top candidates
    println()
    println("Per-sequence SSE (direction=0, phase=0, mirror=0):")
    seq_sses = NamedTuple[]
    for bits in 0:(2^n - 1)
        sequence = [Int((bits >> (lobe - 1)) & 1) for lobe in 1:n]
        mold = assemble_whole_roi(entry, sequence, centers, rxs, rys;
            half_nm=half_nm, step_nm=step_nm, theta=theta,
            direction=0, phase=0, mirror=0)
        fit = fit_whole_roi_nuisance(observed, backbone_roi, mold)
        push!(seq_sses, (sequence=sequence, sse=fit.sse,
            mold_coeff=fit.mold_coefficient))
    end
    sort!(seq_sses, by=s -> s.sse)
    for (rank, item) in enumerate(seq_sses[1:min(end, 8)])
        println("  #", rank, ": seq=", item.sequence,
            "  SSE=", @sprintf("%.6e", item.sse),
            "  mold_coeff=", @sprintf("%.4f", item.mold_coeff))
    end
end

main()

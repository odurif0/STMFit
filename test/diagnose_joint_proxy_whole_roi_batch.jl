#!/usr/bin/env julia

# Optimized batch real-ROI whole-ROI diagnostic.
# Precomputes per-lobe per-type mold images, then combines for each sequence.
#
# Usage:
#   STMFIT_DATA_DIR=/path/to/data julia --project=. \
#       test/diagnose_joint_proxy_whole_roi_batch.jl --out /tmp/opencode/batch.tsv \
#       240817_064.sxm ...

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
    source = WR_REG.ProxySource("stm_dft_periodic080", "diagnostic",
        path, WR_REG.sha256_file(path), 0.50, 1.0, true)
    return WR_REG.ProxyEntry(source, templates)
end

_template_pixels(entry, typ, parity, mirror) =
    only(t for t in entry.templates if t.type == typ && t.parity == parity && t.mirror == mirror).pixels

function _template_matrix(pixels, grid_n)
    length(pixels) == grid_n^2 || throw(ArgumentError("template/grid mismatch"))
    return permutedims(reshape(Float64.(pixels), grid_n, grid_n))
end

function _sample_template(matrix, t, u, half_nm, step_nm)
    tol = 32eps(Float64) * max(1.0, half_nm)
    (t < -half_nm - tol || t > half_nm + tol ||
     u < -half_nm - tol || u > half_nm + tol) && return 0.0
    n = size(matrix, 1)
    qt = clamp((t + half_nm) / step_nm + 1.0, 1.0, Float64(n))
    qu = clamp((u + half_nm) / step_nm + 1.0, 1.0, Float64(n))
    it = min(floor(Int, qt), n - 1)
    iu = min(floor(Int, qu), n - 1)
    ft, fu = qt - it, qu - iu
    return (1-ft)*(1-fu)*matrix[iu,it] + ft*(1-fu)*matrix[iu,it+1] +
           (1-ft)*fu*matrix[iu+1,it] + ft*fu*matrix[iu+1,it+1]
end

function precompute_lobe_images(entry, centers, n, xs, ys, theta,
        half_nm, step_nm, direction, phase, mirror)
    grid_n = length(collect(-half_nm:step_nm:half_nm))
    cos_th, sin_th = cos(theta), sin(theta)
    images = Vector{NTuple{2,Matrix{Float64}}}()
    for lobe in 1:n
        parity = mod((direction == 0 ? lobe : n - lobe + 1) - 1 + phase, 2)
        mat0 = _template_matrix(_template_pixels(entry, 0, parity, mirror), grid_n)
        mat1 = _template_matrix(_template_pixels(entry, 1, parity, mirror), grid_n)
        img0 = zeros(Float64, length(ys), length(xs))
        img1 = zeros(Float64, length(ys), length(xs))
        cx, cy = centers[lobe]
        @inbounds for iy in eachindex(ys), ix in eachindex(xs)
            dx, dy = xs[ix] - cx, ys[iy] - cy
            t = dx * cos_th + dy * sin_th
            u = -dx * sin_th + dy * cos_th
            img0[iy, ix] = _sample_template(mat0, t, u, half_nm, step_nm)
            img1[iy, ix] = _sample_template(mat1, t, u, half_nm, step_nm)
        end
        push!(images, (img0, img1))
    end
    return images
end

function _active_fit(y, predictors, active)
    design = ones(Float64, length(y), length(active) + 1)
    for (column, predictor) in enumerate(active)
        design[:, column + 1] .= predictors[predictor]
    end
    coefficients = design \ y
    any(coefficients[2:end] .< -sqrt(eps(Float64))) && return nothing
    coefficients[2:end] .= max.(coefficients[2:end], 0.0)
    residual = y .- design * coefficients
    full = zeros(Float64, length(predictors))
    for (column, predictor) in enumerate(active)
        full[predictor] = coefficients[column + 1]
    end
    return (; sse=sum(abs2, residual), coefficients=full)
end

function _nnls_fit(y, backbone_v, mold_v)
    predictors = [backbone_v, mold_v]
    candidates = filter(!isnothing, [_active_fit(y, predictors, active)
        for active in (Int[], [1], [2], [1, 2])])
    best = candidates[argmin(getproperty.(candidates, :sse))]
    return best.sse, best.coefficients[1], best.coefficients[2]
end

function _null_sse(y, backbone_v)
    predictors = [backbone_v]
    candidates = filter(!isnothing, [_active_fit(y, predictors, active)
        for active in (Int[], [1])])
    return minimum(getproperty.(candidates, :sse))
end

function process_file(filename, selected_n, data_dir, entry, half_nm, step_nm)
    path = joinpath(data_dir, filename)
    isfile(path) || return (; file=filename, status="not_found", n=0, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")

    image = STMSXMIO.read_sxm(path)
    pattern = Inference._default_pattern_config(path, "/tmp/opencode/batch")
    chain = Inference._default_chain_config()
    results, _, context = redirect_stdout(devnull) do
        GaussianFit2D.chain_gaussian_sweep(image, pattern, chain)
    end
    views = Inference._build_views(image, pattern,
        Main.JointProxyCandidateViews.build_view_data)
    report = Main.JointProxyCandidateViews.extract_candidate_views(
        results, context, chain, views; patch_half_nm=0.32, patch_step_nm=0.08)
    isempty(report.candidates) && return (; file=filename, status="no_cand", n=0, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")
    candidate_index = findfirst(item -> item.n == selected_n, report.candidates)
    candidate_index === nothing && return (; file=filename, status="selected_n_missing", n=selected_n, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")
    candidate = report.candidates[candidate_index]
    n = candidate.n

    fwd = nothing
    for v in views; v.label == "fwd" && v.present && (fwd = v; break) end
    fwd === nothing && return (; file=filename, status="no_fwd", n=n, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")

    ax, ay = context.axisctx.axis
    theta = atan(ay, ax)
    centers = [(l.x_nm, l.y_nm) for l in candidate.lobes]
    margin = half_nm + 0.1
    ix = findall(t -> minimum(c[1] for c in centers)-margin <= t <= maximum(c[1] for c in centers)+margin, fwd.xs)
    iy = findall(u -> minimum(c[2] for c in centers)-margin <= u <= maximum(c[2] for c in centers)+margin, fwd.ys)
    isempty(ix) && return (; file=filename, status="empty_roi", n=n, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")
    rxs, rys, rz = fwd.xs[ix], fwd.ys[iy], fwd.z[iy, ix]

    best_result = results[findfirst(r -> r.n == n && r.success && r.valid, results)]
    best_result === nothing && return (; file=filename, status="no_fit", n=n, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")
    mv = GaussianFit2D._chain_model_values(
        vec(repeat(reshape(rxs, 1, :), length(rys), 1)),
        vec(repeat(reshape(rys, :, 1), 1, length(rxs))),
        best_result.params, n, context.axisctx, chain;
        amp_min=best_result.amp_min, amp_range=best_result.amp_range)
    backbone = reshape(mv, length(rys), length(rxs))
    observed = rz
    good = isfinite.(observed) .& isfinite.(backbone)
    count(good) < 20 && return (; file=filename, status="few_px", n=n, sse_red=NaN, gap_pct=NaN, mc=NaN, seq="")

    y = observed[good]; bv = backbone[good]
    null_sse = _null_sse(y, bv)

    sequence_sse = Dict{String,Float64}()
    sequence_mc = Dict{String,Float64}()
    for direction in (0,1), phase in (0,1), mirror in (0,1)
        images = precompute_lobe_images(entry, centers, n, rxs, rys, theta,
            half_nm, step_nm, direction, phase, mirror)
        for bits in 0:(2^n-1)
            mold = zeros(Float64, length(rys), length(rxs))
            for lobe in 1:n
                mold .+= (bits >> (lobe-1)) & 1 == 0 ? images[lobe][1] : images[lobe][2]
            end
            mv_good = mold[good]
            sse, bc, mc = _nnls_fit(y, bv, mv_good)
            sequence = join([Int((bits >> (l-1)) & 1) for l in 1:n], "")
            if sse < get(sequence_sse, sequence, Inf)
                sequence_sse[sequence] = sse
                sequence_mc[sequence] = mc
            end
        end
    end

    ranked = sort(collect(sequence_sse); by=last)
    best_seq, best_sse = first(ranked)
    second_sse = length(ranked) > 1 ? last(ranked[2]) : Inf
    best_mc = sequence_mc[best_seq]

    return (; file=filename, status="ok", n=n,
        sse_red=null_sse > 0 ? 100.0*(1.0 - best_sse/null_sse) : NaN,
        gap_pct=null_sse > 0 ? 100.0*(second_sse-best_sse)/null_sse : NaN,
        mc=best_mc, seq=best_seq)
end

function read_selected_counts(path)
    lines = readlines(path)
    isempty(lines) && error("empty selection summary: $path")
    header = split(first(lines), '\t')
    file_col = findfirst(==("filepath"), header)
    selected_col = findfirst(==("N_selected"), header)
    file_col === nothing && error("selection summary missing filepath")
    selected_col === nothing && error("selection summary missing N_selected")
    counts = Dict{String,Int}()
    for line in lines[2:end]
        fields = split(line, '\t')
        length(fields) >= max(file_col, selected_col) || continue
        parsed = tryparse(Int, fields[selected_col])
        parsed === nothing || (counts[basename(fields[file_col])] = parsed)
    end
    return counts
end

function main()
    out = "/tmp/opencode/real_roi_batch.tsv"; selection_summary = ""; files = String[]
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--out"
            out = ARGS[i+1]; i += 2
        elseif ARGS[i] == "--selection-summary"
            selection_summary = ARGS[i+1]; i += 2
        else
            push!(files, String(ARGS[i])); i += 1
        end
    end
    isempty(files) && error("at least one file required")
    isempty(selection_summary) && error("--selection-summary is required (label-free production output)")
    selected_counts = read_selected_counts(selection_summary)
    data_dir = something(get(ENV, "STMFIT_DATA_DIR", nothing), ".")
    mold_path = "/tmp/opencode/chitosan_connected_molds_dft_m030_h050_periodic080_half080.tsv"
    entry = load_converged_entry(mold_path)

    println("Batch real-ROI whole-ROI diagnostic (optimized)")
    println("  files: ", length(files))
    println("  selection summary: ", selection_summary)
    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "file\tstatus\tn\tsse_reduction_pct\tsequence_gap_pct\tmold_coeff\tbest_sequence")
        flush(io)
        for file in files
            haskey(selected_counts, file) || error("file missing from selection summary: $file")
            r = process_file(file, selected_counts[file], data_dir, entry, 0.80, 0.08)
            println(io, join([r.file, r.status, r.n,
                @sprintf("%.4f", r.sse_red), @sprintf("%.6f", r.gap_pct),
                @sprintf("%.6f", r.mc), r.seq], '\t'))
            flush(io)
            status_str = r.status == "ok" ?
                " SSE↓=$(@sprintf("%.2f%%", r.sse_red)) gap=$(@sprintf("%.4f%%", r.gap_pct)) mold=$(@sprintf("%.4f", r.mc)) seq=$(r.seq)" : ""
            println("  ", r.file, ": N=", r.n, " status=", r.status, status_str)
        end
    end
    println("  output: ", out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

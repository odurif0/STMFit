#!/usr/bin/env julia
# Joint shape fit of every molecule track of a run of
# run_molecule_consensus_chitosan.jl on its registered repeat scans (campaign J,
# stage 2; docs/src/journal.md, 2026-10-26, for the design, fixed before this
# script was run). Label-free: scan geometry, physical lobes, tracks and raw
# topography; no benchmark label is read.
#
# Usage:
#   julia -t 8 --project=. scripts/pipeline/fit_joint_shape.jl --run RUN_DIR --data-dir RAW_DIR \
#       --outdir NEW_DIR [--cap 8] [--step 0.02] [--margin 0.8] [--tracks T1,T2,...] [--save-stacks]
#
# For each track of at least 3 scans (at most --cap, evenly spaced in file
# order): rigid registration of the scans onto a reference scan by their shared
# physical lobes, stack on a --step grid of the molecule frame (--margin nm
# across and beyond the end lobes), and one elliptical Gaussian per unit; the
# same on the odd and even scans when each half has at least 2 scans.
# Outputs: units.tsv (one row per unit and half), tracks.tsv (fit status), and
# with --save-stacks the stacks (stacks/<track>_<half>.tsv: first row t, first
# column u).
using STMFit.JointShape, STMFit.MoleculeConsensus, STMFit.ReconstructedUnitAssignment
using STMFit.UnitReliability: capped_files
import STMFit.MoleculeConsensus: _bilinear
using STMSXMIO
using Printf, Statistics

function shape_fit_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--cap" => "8", "--step" => "0.02", "--margin" => "0.8"); i = 1
    while i <= length(args)
        k = args[i]
        if k == "--save-stacks"
            o[k] = "yes"; i += 1; continue
        end
        k in ("--run", "--data-dir", "--outdir", "--cap", "--step", "--margin", "--tracks") && i < length(args) ||
            error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    for k in ("--run", "--data-dir", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

function raw_index(dir)
    index = Dict{String,String}()
    for (d, _, names) in walkdir(dir; follow_symlinks=true), n in names
        endswith(lowercase(n), ".sxm") && (index[n] = joinpath(d, n))
    end
    return index
end

physical_index(p) = parse(Int, split(p, '_')[end])

"Tracks (file order) and, per scan, its lobes as (k, x_nm, y_nm, physical lobe)."
function load_tracks(run)
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    _, geo = read_table(joinpath(run, "scan_geometry.tsv"))
    _, pred = read_table(joinpath(run, "predictions.tsv"))
    physical = Dict((r["file"], parse(Int, r["lobe"])) => r["physical_lobe"] for r in pred)
    lobes = Dict{String,Vector{Tuple{Int,Float64,Float64,String}}}()
    for r in geo
        f = r["file"]; k = parse(Int, r["lobe"])
        push!(get!(lobes, f, Tuple{Int,Float64,Float64,String}[]),
              (k, parse(Float64, r["x_nm"]), parse(Float64, r["y_nm"]), get(physical, (f, k), "")))
    end
    foreach(v -> sort!(v; by=first), values(lobes))
    tracks = Dict{String,Vector{String}}()
    for r in cons
        haskey(lobes, r["file"]) && push!(get!(tracks, r["track"], String[]), r["file"])
    end
    foreach(sort!, values(tracks))
    return tracks, lobes
end

modal(v) = first(sort(unique(v); by=x -> (-count(==(x), v), x)))

"Scan data for the stack: frame, flattened forward topography, absolute positions of its physical lobes."
function scan_data(path, lobes)
    img = read_sxm(path); fr = scan_frame(img)
    ch = nothing
    for c in img.channels
        lowercase(c.name) == "z" && lowercase(c.direction) == "fwd" && (ch = c)
    end
    ch === nothing && error("No forward Z channel: $path")
    z = STMSXMIO.preprocess_observed_channel(img, ch; stride=1, flatten="plane+rows", smooth_radius_px=0, plane_rank_rtol=1e-12).z
    pos = Dict(l[4] => absolute_xy(fr, l[2], l[3]) for l in lobes if !isempty(l[4]))
    return fr, z, pos
end

function sampler(F, T, fr, z)
    inv = invert_rigid(T)
    sx = (fr.width - 1) / fr.extent[1]; sy = (fr.height - 1) / fr.extent[2]
    return (t, u) -> begin
        q = apply_rigid(inv, from_frame(F, t, u))
        x, y = image_xy(fr, q...)
        _bilinear(z, y * sy + 1, x * sx + 1)
    end
end

function fit_track(track, files, lobes, raw, cap, step, margin)
    rows = Dict{String,Any}[]; status = Dict{String,Any}[]
    used = capped_files(files, cap)
    nmode = modal([length(lobes[f]) for f in used])
    iref = findfirst(f -> length(lobes[f]) == nmode && all(!isempty(l[4]) for l in lobes[f]), used)
    if iref === nothing
        push!(status, Dict("track" => track, "half" => "all", "scans" => 0, "units" => nmode, "converged" => false,
                           "iterations" => 0, "rms_rel" => NaN, "reason" => "no reference scan"))
        return rows, status, nothing
    end
    data = Dict(f => scan_data(raw[f], lobes[f]) for f in used)
    ref = used[iref]; refpos = data[ref][3]
    names = sort(collect(keys(refpos)); by=physical_index)
    F = molecule_frame([refpos[p] for p in names])
    init = [to_frame(F, refpos[p]) for p in names]
    samplers = Dict{String,Any}()
    for f in used
        shared = [p for p in names if haskey(data[f][3], p)]
        length(shared) >= 3 || continue
        T = rigid_fit([data[f][3][p] for p in shared], [refpos[p] for p in shared])
        samplers[f] = sampler(F, T, data[f][1], data[f][2])
    end
    reg = [f for f in used if haskey(samplers, f)]
    ts = collect(minimum(first.(init)) - margin:step:maximum(first.(init)) + margin + step / 2)
    us = collect(-margin:step:margin + step / 2)
    chain_t = (minimum(first.(init)), maximum(first.(init)))
    halves = ["all" => reg]
    length(reg) >= 4 && append!(halves, ["odd" => reg[1:2:end], "even" => reg[2:2:end]])
    stacks = Dict{String,Any}()
    for (half, fs) in halves
        length(fs) >= (half == "all" ? 3 : 2) || continue
        img, _, nused = stack_scans([samplers[f] for f in fs], ts, us; chain_t)
        st = Dict{String,Any}("track" => track, "half" => half, "scans" => nused, "units" => length(init))
        try
            fit, bg, info = fit_lobes(ts, us, img, init)
            amed = median(l.amplitude for l in fit); n = length(fit)
            for (r, (p, l)) in enumerate(zip(names, fit))
                push!(rows, Dict{String,Any}("track" => track, "half" => half, "scans" => nused, "unit" => p,
                    "rank" => r, "n_units" => n, "d" => min(r, n + 1 - r), "t_nm" => l.t, "u_nm" => l.u,
                    "amplitude" => l.amplitude, "amplitude_rel" => l.amplitude / amed,
                    "sigma_par_nm" => l.sigma_par, "sigma_perp_nm" => l.sigma_perp, "theta_deg" => rad2deg(l.theta)))
            end
            merge!(st, Dict("converged" => info.converged, "iterations" => info.iterations, "rms_rel" => info.rms_rel, "reason" => "ok"))
        catch e
            merge!(st, Dict("converged" => false, "iterations" => 0, "rms_rel" => NaN, "reason" => sprint(showerror, e)))
        end
        push!(status, st)
        stacks[half] = img
    end
    return rows, status, (ts, us, stacks)
end

fmtv(x) = x isa AbstractFloat ? (isfinite(x) ? @sprintf("%.6f", x) : "NaN") : string(x)

function shape_fit_main(args=ARGS)
    o = shape_fit_options(args)
    o === nothing && return println("fit_joint_shape.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR [--cap 8] [--step 0.02] [--margin 0.8] [--tracks T1,T2] [--save-stacks]")
    cap, step, margin = parse(Int, o["--cap"]), parse(Float64, o["--step"]), parse(Float64, o["--margin"])
    tracks, lobes = load_tracks(o["--run"])
    raw = raw_index(o["--data-dir"])
    names = sort([t for (t, fs) in tracks if length(fs) >= 3]; by=t -> (length(t), t))
    haskey(o, "--tracks") && (keep = Set(split(o["--tracks"], ',')); names = [t for t in names if t in keep])
    results = Vector{Any}(undef, length(names))
    Threads.@threads for i in eachindex(names)
        t = names[i]
        results[i] = try
            fit_track(t, tracks[t], lobes, raw, cap, step, margin)
        catch e
            (Dict{String,Any}[], [Dict{String,Any}("track" => t, "half" => "all", "scans" => 0, "units" => 0, "converged" => false,
                                                   "iterations" => 0, "rms_rel" => NaN, "reason" => sprint(showerror, e))], nothing)
        end
    end
    out = o["--outdir"]; mkpath(out)
    ucols = ["track", "half", "scans", "unit", "rank", "n_units", "d", "t_nm", "u_nm", "amplitude", "amplitude_rel",
             "sigma_par_nm", "sigma_perp_nm", "theta_deg"]
    scols = ["track", "half", "scans", "units", "converged", "iterations", "rms_rel", "reason"]
    open(joinpath(out, "units.tsv"), "w") do io
        println(io, join(ucols, '\t'))
        for r in results, row in r[1]
            println(io, join((fmtv(row[c]) for c in ucols), '\t'))
        end
    end
    open(joinpath(out, "tracks.tsv"), "w") do io
        println(io, join(scols, '\t'))
        for r in results, row in r[2]
            println(io, join((replace(fmtv(row[c]), '\t' => ' ', '\n' => ' ') for c in scols), '\t'))
        end
    end
    if haskey(o, "--save-stacks")
        mkpath(joinpath(out, "stacks"))
        for (t, r) in zip(names, results)
            r[3] === nothing && continue
            ts, us, stacks = r[3]
            for (half, img) in stacks
                open(joinpath(out, "stacks", "$(t)_$(half).tsv"), "w") do io
                    println(io, join(vcat(["u\\t"], [@sprintf("%.3f", x) for x in ts]), '\t'))
                    for (i, u) in enumerate(us)
                        println(io, join(vcat([@sprintf("%.3f", u)], [fmtv(v) for v in img[i, :]]), '\t'))
                    end
                end
            end
        end
    end
    ok = count(r -> any(s -> s["half"] == "all" && s["reason"] == "ok", r[2]), results)
    println("Joint shape fit: $(length(names)) tracks, $ok fitted (all scans); outputs in $out")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && shape_fit_main()

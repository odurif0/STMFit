#!/usr/bin/env julia
# Metadata and tip-position integrity check of the molecule tracks of a run of
# run_molecule_consensus_chitosan.jl. Label-free: SXM headers and fitted lobe
# positions only; no benchmark label, no grade; nothing in the run is changed.
#
# Usage:
#   julia --project=. scripts/data/check_track_integrity.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR
#
# For every track (>= 2 scans) of RUN_DIR/consensus/consensus.tsv it recomputes,
# in the absolute piezo frame given by the headers (offset, range, angle):
#   sessions          distinct acquisition dates in the track (should be 1);
#   span_min, max_gap_min   acquisition time span and largest gap between consecutive members;
#   centroid_spread_nm      largest distance of a scan's chain centroid from the track median;
#   axis_spread_deg   largest difference between the absolute chain-axis directions
#                     (a molecule does not turn between scans), modulo 180 degrees;
#   n_range           max - min of the per-scan counts;
#   interleaved       scans of other tracks acquired between the first and last member.
# scans.tsv lists every scan (track, time, absolute centroid, absolute axis direction, count).
# A track is flagged when sessions > 1, max_gap_min > 30, centroid_spread_nm > 6,
# axis_spread_deg > 30 or n_range > 2. These thresholds are diagnostic only (they
# follow the collection limits of config/molecule_consensus.toml) and select nothing.
using STMFit.MoleculeConsensus, STMFit.ReconstructedUnitAssignment
using STMSXMIO, Dates, Printf, Statistics

function integrity_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--data-dir", "--outdir") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--data-dir", "--outdir")) || error("Required: --run --data-dir --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

circ_diff_deg(a, b) = (d = abs(mod(a - b + 90, 180) - 90); d)

function raw_files(dir)
    index = Dict{String,String}()
    for (d, _, names) in walkdir(dir; follow_symlinks=true), name in names
        endswith(lowercase(name), ".sxm") && (index[name] = joinpath(d, name))
    end
    return index
end

function integrity_main(args=ARGS)
    o = integrity_options(args)
    o === nothing && return println("check_track_integrity.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR")
    run = o["--run"]
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    _, lobes = lobe_table(joinpath(run, "scan_geometry.tsv"); required=["x_nm", "y_nm", "axis_x", "axis_y"])
    raw = raw_files(o["--data-dir"])
    track_of = Dict(r["file"] => parse(Int, r["track"]) for r in cons)
    n_of = Dict(r["file"] => parse(Int, r["N_scan"]) for r in cons)
    frames = Dict{String,ScanFrame}()
    centroid = Dict{String,NTuple{2,Float64}}(); axis_deg = Dict{String,Float64}()
    for f in sort(collect(keys(track_of)))
        frames[f] = scan_frame(read_sxm(raw[f]))
        keys_f = sort([k for k in keys(lobes) if first(k) == f]; by=last)
        x = [parse(Float64, lobes[k]["x_nm"]) for k in keys_f]; y = [parse(Float64, lobes[k]["y_nm"]) for k in keys_f]
        cx, cy = mean(x), mean(y)
        ax, ay = parse(Float64, lobes[first(keys_f)]["axis_x"]), parse(Float64, lobes[first(keys_f)]["axis_y"])
        X0, Y0 = absolute_xy(frames[f], cx, cy); X1, Y1 = absolute_xy(frames[f], cx + ax, cy + ay)
        centroid[f] = (X0, Y0); axis_deg[f] = mod(rad2deg(atan(Y1 - Y0, X1 - X0)), 180)
    end
    order = sort(collect(keys(frames)); by=f -> (frames[f].acquired, f))
    rows = Dict{String,String}[]
    for t in sort(unique(values(track_of)))
        members = sort([f for f in keys(track_of) if track_of[f] == t]; by=f -> (frames[f].acquired, f))
        length(members) >= 2 || continue
        times = [frames[f].acquired for f in members]
        gaps = [Dates.value(times[k+1] - times[k]) / 60000 for k in 1:length(members)-1]
        cx = median(first.(getindex.(Ref(centroid), members))); cy = median(last.(getindex.(Ref(centroid), members)))
        spread = maximum(hypot(centroid[f][1] - cx, centroid[f][2] - cy) for f in members)
        axspread = maximum(circ_diff_deg(axis_deg[a], axis_deg[b]) for a in members, b in members)
        ns = [n_of[f] for f in members]
        i1 = findfirst(==(first(members)), order); i2 = findfirst(==(last(members)), order)
        inter = count(f -> track_of[f] != t, order[i1:i2])
        sessions = length(unique(frames[f].session for f in members))
        flags = String[]
        sessions > 1 && push!(flags, "sessions")
        maximum(gaps) > 30 && push!(flags, "time_gap")
        spread > 6 && push!(flags, "centroid_spread")
        axspread > 30 && push!(flags, "axis_spread")
        maximum(ns) - minimum(ns) > 2 && push!(flags, "n_range")
        push!(rows, Dict("track" => string(t), "scans" => string(length(members)), "sessions" => string(sessions),
            "span_min" => @sprintf("%.1f", Dates.value(times[end] - times[1]) / 60000), "max_gap_min" => @sprintf("%.1f", maximum(gaps)),
            "centroid_spread_nm" => @sprintf("%.2f", spread), "axis_spread_deg" => @sprintf("%.1f", axspread),
            "n_range" => string(maximum(ns) - minimum(ns)), "interleaved" => string(inter),
            "flags" => join(flags, ","), "first_scan" => first(members)))
    end
    mkpath(o["--outdir"])
    per_scan = [Dict("track" => string(track_of[f]), "file" => f, "acquired" => string(frames[f].acquired),
                     "X_nm" => @sprintf("%.3f", centroid[f][1]), "Y_nm" => @sprintf("%.3f", centroid[f][2]),
                     "axis_deg" => @sprintf("%.1f", axis_deg[f]), "N_scan" => string(n_of[f])) for f in order]
    write_table(joinpath(o["--outdir"], "scans.tsv"), ["track", "file", "acquired", "X_nm", "Y_nm", "axis_deg", "N_scan"], per_scan)
    cols = ["track", "scans", "sessions", "span_min", "max_gap_min", "centroid_spread_nm", "axis_spread_deg", "n_range", "interleaved", "flags", "first_scan"]
    write_table(joinpath(o["--outdir"], "tracks.tsv"), cols, rows)
    flagged = [r for r in rows if !isempty(r["flags"])]
    q(k) = sort(parse.(Float64, [r[k] for r in rows]))
    open(joinpath(o["--outdir"], "report.md"), "w") do io
        println(io, "# Track integrity: $run\n")
        println(io, "$(length(rows)) tracks of at least 2 scans ($(sum(parse(Int, r["scans"]) for r in rows)) scans); flagged: $(length(flagged)).\n")
        println(io, "| check | median | q95 | max | tracks flagged |\n|---|---:|---:|---:|---:|")
        for (k, fl) in (("span_min", ""), ("max_gap_min", "time_gap"), ("centroid_spread_nm", "centroid_spread"), ("axis_spread_deg", "axis_spread"), ("n_range", "n_range"), ("interleaved", ""))
            v = q(k); println(io, "| $k | $(@sprintf("%.1f", median(v))) | $(@sprintf("%.1f", quantile(v, 0.95))) | $(@sprintf("%.1f", maximum(v))) | ",
                              isempty(fl) ? "-" : string(count(r -> occursin(fl, r["flags"]), rows)), " |")
        end
        println(io, "\nTracks with more than one session: $(count(r -> r["sessions"] != "1", rows)).")
        println(io, "\n## Flagged tracks\n\n| track | scans | flags | gap min | centroid nm | axis deg | N range |\n|---|---:|---|---:|---:|---:|---:|")
        for r in sort(flagged; by=r -> -parse(Int, r["scans"]))
            println(io, "| $(r["track"]) | $(r["scans"]) | $(r["flags"]) | $(r["max_gap_min"]) | $(r["centroid_spread_nm"]) | $(r["axis_spread_deg"]) | $(r["n_range"]) |")
        end
    end
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && integrity_main()

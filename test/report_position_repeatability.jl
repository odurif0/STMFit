#!/usr/bin/env julia
# Label-free lobe-position repeatability across repeated scans of one molecule.
#
# Usage:
#   julia --project=. test/report_position_repeatability.jl --run RUN_DIR \
#       --data-dir RAW_DIR --outdir NEW_DIR
#
# RUN_DIR is a run_molecule_consensus_chitosan.jl output. For every molecule
# track, the scans with the track's modal final count are matched lobe by lobe
# (rank along the chain axis, as in the fusion stage). Lobe positions are
# mapped to the absolute piezo frame and aligned on the track's first scan:
#   image  : drift from the consecutive-scan image registration only
#            (independent of the lobe fits);
#   rigid  : additionally, the per-scan translation that best superposes the
#            lobes (least squares) is removed;
#   procrustes : per-scan translation and rotation (no scaling) onto the mean
#            shape, iterated; separates rotation between scans from noise.
# The scatter of each physical lobe around its mean is split along the chain
# axis and across it. No label, expected count or sequence is read.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
include(joinpath(@__DIR__, "lib", "molecule_consensus.jl"))
include(joinpath(@__DIR__, "lib", "molecule_fusion.jl"))
using .ReconstructedUnitAssignment, .MoleculeConsensus, .MoleculeFusion
using STMSXMIO, Statistics, Printf

function repeat_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--run", "--data-dir", "--outdir") && i < length(args) || error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--data-dir", "--outdir")) || error("Required: --run --data-dir --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

"Pooled per-component SD over groups: sqrt(sum of squared residuals / sum of (n - 1))."
pooled_sd(ss, dof) = dof > 0 ? sqrt(ss / dof) : NaN

const Shape = Vector{NTuple{2,Float64}}

lobe_means(Q, scans, n) = [(mean(Q[f][r][1] for f in scans), mean(Q[f][r][2] for f in scans)) for r in 1:n]

"Remove each scan's mean offset from the across-scan lobe means."
function align_translation(P, scans, n)
    m = lobe_means(P, scans, n)
    Dict(f => (o = (mean(P[f][r][1] - m[r][1] for r in 1:n), mean(P[f][r][2] - m[r][2] for r in 1:n));
               [(p[1] - o[1], p[2] - o[2]) for p in P[f]]) for f in scans)
end

"Generalized Procrustes without scaling: per-scan rotation and translation onto the mean shape."
function align_procrustes(P, scans, n; iterations=5)
    Q = Dict(f => copy(P[f]) for f in scans); rot = Dict(f => 0.0 for f in scans)
    for _ in 1:iterations
        m = lobe_means(Q, scans, n); mc = (mean(p[1] for p in m), mean(p[2] for p in m))
        for f in scans
            q = P[f]; qc = (mean(p[1] for p in q), mean(p[2] for p in q))
            sxy = sum((q[r][1]-qc[1])*(m[r][2]-mc[2]) - (q[r][2]-qc[2])*(m[r][1]-mc[1]) for r in 1:n)
            sxx = sum((q[r][1]-qc[1])*(m[r][1]-mc[1]) + (q[r][2]-qc[2])*(m[r][2]-mc[2]) for r in 1:n)
            th = atan(sxy, sxx); c, sn = cos(th), sin(th); rot[f] = th
            Q[f] = [(mc[1] + c*(p[1]-qc[1]) - sn*(p[2]-qc[2]), mc[2] + sn*(p[1]-qc[1]) + c*(p[2]-qc[2])) for p in q]
        end
    end
    return Q, rot
end

function repeatability(run, rawdir)
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    _, final = read_table(joinpath(run, "consensus", "consensus_summary_final.tsv"))
    _, pairs = read_table(joinpath(run, "consensus", "pairs.tsv"))
    tracks = Dict(r["file"] => parse(Int, r["track"]) for r in cons)
    counts = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in final)
    _, geo = lobe_table(joinpath(run, "assignment", "features.tsv"); required=["x_nm", "y_nm"])
    raw = Dict{String,String}()
    for (d, _, names) in walkdir(rawdir), n in names
        endswith(lowercase(n), ".sxm") && (raw[n] = joinpath(d, n))
    end
    frames = Dict(f => scan_frame(read_sxm(raw[f])) for f in keys(counts))
    lobes_abs = Dict{String,Vector{NTuple{2,Float64}}}()
    for key in sort(collect(keys(geo)))
        r = geo[key]
        push!(get!(lobes_abs, key[1], NTuple{2,Float64}[]),
              absolute_xy(frames[key[1]], parse(Float64, r["x_nm"]), parse(Float64, r["y_nm"])))
    end
    order = sort(collect(keys(counts)); by=f -> (frames[f].acquired, f))
    phys = physical_lobes(order, tracks, counts, lobes_abs)
    link = Dict((p["a"], p["b"]) => (parse(Float64, p["tx_nm"]), parse(Float64, p["ty_nm"])) for p in pairs if p["linked"] == "true")
    pos = Dict(f => k for (k, f) in enumerate(order))
    byphys = Dict{String,Vector{Tuple{String,Int}}}()
    for (key, pk) in phys; push!(get!(byphys, pk, Tuple{String,Int}[]), key); end
    trackscans = Dict{Int,Vector{String}}()
    for key in keys(phys)
        f = key[1]; t = tracks[f]
        v = get!(trackscans, t, String[]); f in v || push!(v, f)
    end
    rows = Dict{String,Any}[]; spacing_rows = Dict{String,Any}[]; rotation_rows = Dict{String,Any}[]
    for (t, fs) in sort(collect(trackscans))
        length(fs) >= 2 || continue
        sort!(fs; by=f -> pos[f]); ref = fs[1]; n = counts[ref]
        # cumulative image-registration drift ref -> f (consecutive links of the track)
        drift = Dict{String,NTuple{2,Float64}}(ref => (0.0, 0.0))
        for f in fs[2:end]
            d = (0.0, 0.0); ok = true
            for k in pos[ref]:pos[f]-1
                l = get(link, (order[k], order[k+1]), nothing)
                l === nothing && (ok = false; break)
                d = (d[1] + l[1], d[2] + l[2])
            end
            ok && (drift[f] = d)
        end
        used = [f for f in fs if haskey(drift, f)]
        length(used) >= 2 || continue
        # positions per scan, indexed by chain rank
        P = Dict{String,Vector{NTuple{2,Float64}}}()
        for f in used
            v = Vector{NTuple{2,Float64}}(undef, n)
            for lobe in 1:n
                r = parse(Int, split(phys[(f, lobe)], "_")[end])
                a = lobes_abs[f][lobe]
                v[r] = (a[1] - drift[f][1], a[2] - drift[f][2])
            end
            P[f] = v
        end
        ax = (P[ref][end][1] - P[ref][1][1], P[ref][end][2] - P[ref][1][2]); L = hypot(ax...)
        ax = (ax[1] / L, ax[2] / L); pe = (-ax[2], ax[1])
        for mode in ("image", "rigid", "procrustes")
            if mode == "image"
                Q = P
            elseif mode == "procrustes"
                Q, rot = align_procrustes(P, used, n)
                for f in used
                    push!(rotation_rows, Dict("track" => string(t), "file" => f,
                        "scan_angle_deg" => @sprintf("%.2f", frames[f].angle_deg),
                        "rotation_deg" => @sprintf("%.3f", rad2deg(rot[f] - rot[ref]))))
                end
            else
                Q = align_translation(P, used, n)
            end
            for r in 1:n
                m = (mean(Q[f][r][1] for f in used), mean(Q[f][r][2] for f in used))
                res = [(Q[f][r][1] - m[1], Q[f][r][2] - m[2]) for f in used]
                along = [d[1] * ax[1] + d[2] * ax[2] for d in res]
                across = [d[1] * pe[1] + d[2] * pe[2] for d in res]
                dof = length(used) - 1 - (mode == "rigid" ? (length(used) - 1) / n : mode == "procrustes" ? 1.5 * (length(used) - 1) / n : 0)
                push!(rows, Dict("track" => t, "rank" => r, "N" => n, "scans" => length(used), "mode" => mode,
                    "position" => r in (1, n) ? "end" : "interior",
                    "ss_along" => sum(abs2, along), "ss_across" => sum(abs2, across), "dof" => dof))
            end
        end
        for f in used, r in 1:n-1
            push!(spacing_rows, Dict("track" => t, "file" => f, "gap" => r,
                "spacing_nm" => hypot(P[f][r+1][1] - P[f][r][1], P[f][r+1][2] - P[f][r][2])))
        end
    end
    return rows, spacing_rows, rotation_rows
end

function main(args=ARGS)
    o = repeat_options(args)
    o === nothing && return println("report_position_repeatability.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR")
    rows, spacing, rotations = repeatability(abspath(o["--run"]), o["--data-dir"])
    out = o["--outdir"]; mkpath(out)
    write_table(joinpath(out, "scan_rotation.tsv"), ["track", "file", "scan_angle_deg", "rotation_deg"], rotations)
    write_table(joinpath(out, "lobe_scatter.tsv"), ["track", "rank", "N", "scans", "mode", "position", "sd_along_nm", "sd_across_nm"],
        [Dict("track" => string(r["track"]), "rank" => string(r["rank"]), "N" => string(r["N"]), "scans" => string(r["scans"]),
              "mode" => r["mode"], "position" => r["position"],
              "sd_along_nm" => @sprintf("%.4f", pooled_sd(r["ss_along"], r["dof"])),
              "sd_across_nm" => @sprintf("%.4f", pooled_sd(r["ss_across"], r["dof"]))) for r in rows])
    summary = Dict{String,String}[]
    for mode in ("image", "rigid", "procrustes"), pos in ("all", "end", "interior")
        sel = [r for r in rows if r["mode"] == mode && (pos == "all" || r["position"] == pos)]
        isempty(sel) && continue
        dof = sum(r["dof"] for r in sel)
        push!(summary, Dict("mode" => mode, "position" => pos, "lobes" => string(length(sel)),
            "tracks" => string(length(unique(r["track"] for r in sel))),
            "sd_along_nm" => @sprintf("%.4f", pooled_sd(sum(r["ss_along"] for r in sel), dof)),
            "sd_across_nm" => @sprintf("%.4f", pooled_sd(sum(r["ss_across"] for r in sel), dof))))
    end
    s = [r["spacing_nm"] for r in spacing]
    gaps = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in spacing; push!(get!(gaps, (r["track"], r["gap"]), Float64[]), r["spacing_nm"]); end
    within = [std(v) for v in values(gaps) if length(v) >= 2]
    push!(summary, Dict("mode" => "spacing", "position" => "all", "lobes" => string(length(s)), "tracks" => string(length(unique(r["track"] for r in spacing))),
        "sd_along_nm" => @sprintf("mean %.4f sd %.4f", mean(s), std(s)), "sd_across_nm" => @sprintf("within-gap sd %.4f", sqrt(mean(abs2, within)))))
    write_table(joinpath(out, "summary.tsv"), ["mode", "position", "lobes", "tracks", "sd_along_nm", "sd_across_nm"], summary)
    for r in summary
        @printf("%-8s %-9s lobes %-5s tracks %-4s along %-24s across %s\n", r["mode"], r["position"], r["lobes"], r["tracks"], r["sd_along_nm"], r["sd_across_nm"])
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

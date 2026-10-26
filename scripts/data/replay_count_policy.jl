#!/usr/bin/env julia
# Offline replay of the per-scan count layers of batch_full.jl on a finished run,
# with each layer on or off, followed by the molecule consensus (campaign K;
# docs/src/journal.md, 2026-10-26, for the variants and the decision rule, fixed
# before this script was run). Label-free: counting summary, config bounds,
# tracks and, for the split halves, SXM acquisition times.
#
# Usage:
#   julia --project=. scripts/data/replay_count_policy.jl --run RUN_DIR --config COUNT_TOML \
#       --outdir NEW_DIR [--data-dir RAW_DIR]
#
# Variants of the per-scan count, from the counting_summary.tsv columns N_eff,
# robust_aicc_N, ambiguous_eff, runnerup_N_eff, delta_GCV_rel_eff and
# support_2D_ell_nm, with the rules of batch_full.jl (support_midpoint_hybrid):
#   V0  robust-AICc guard, then support midpoint (the current 6-mer policy)
#   V1  support midpoint only        V2  robust-AICc guard only
#   V3  neither (N_eff)
# The replay of V0 is compared with N_selected (it must match on a run of that
# policy). Each variant is followed by the consensus: strict majority over the
# scans of a track of at least 3 scans, else the per-scan count.
# Outputs: variant_counts.tsv, split_half.tsv, report.md.
using STMFit.ReconstructedUnitAssignment, STMFit.MoleculeConsensus, STMFitCore
using STMSXMIO: read_sxm
using Printf, Statistics, TOML

function replay_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--config", "--outdir", "--data-dir") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--run", "--config", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Bounds and thresholds the layers use, read from a count config (engine defaults where the config is silent, as batch_full.jl)."
function policy_settings(path)
    c = TOML.parsefile(path); m = c["model"]; s = get(c, "selection", Dict{String,Any}())
    return (spacing_min=Float64(get(m, "spacing_min_nm", 0.35)), spacing_max=Float64(get(m, "spacing_max_nm", 0.75)),
            sigma_max=Float64(get(m, "sigma_parallel_max_nm", 0.509)), max_overlap=Float64(get(m, "max_overlap", 0.60)),
            n_min=Int(get(m, "n_min", 2)), n_max=Int(get(m, "n_max", 14)),
            amb=Float64(get(s, "gcv_ambiguity_rel_threshold", 0.05)), up=Float64(get(s, "support_midpoint_up_gcv_rel_threshold", 0.30)))
end

"Robust-AICc guard of batch_full.jl (_refined_selection): down when robust AICc wants fewer lobes, up by one when GCV is ambiguous."
function guard(n_eff, robust, amb, runner, dgcv, P)
    robust === nothing && return n_eff
    robust < n_eff && return robust
    robust == n_eff + 1 && amb && isfinite(dgcv) && dgcv <= P.amb && runner == n_eff + 1 && return robust
    return n_eff
end

"Support midpoint of batch_full.jl (_support_midpoint_hybrid_selection)."
function midpoint(current, n_eff, support, dgcv, P)
    (isfinite(support) && support > 0) || return current
    smin = effective_spacing_min(P.spacing_min, P.spacing_max, P.sigma_max, P.max_overlap)
    lo, hi = support_n_bounds(support, smin, P.spacing_max; n_min_config=P.n_min, n_max_config=P.n_max)
    mid = round(Int, (lo + hi) / 2)
    current > mid + 1 && return mid
    current > mid && return current - 1
    current < mid && n_eff > current && isfinite(dgcv) && dgcv <= P.up && return current + 1
    return current
end

intor(s) = (v = tryparse(Int, s); v)
floator(s) = (v = tryparse(Float64, s); v === nothing ? NaN : v)

"Per-scan counts of the four variants."
function variants(r, P)
    n = parse(Int, r["N_eff"]); robust = intor(get(r, "robust_aicc_N", "NA"))
    amb = lowercase(get(r, "ambiguous_eff", "false")) == "true"; runner = intor(get(r, "runnerup_N_eff", "NA"))
    dgcv = floator(get(r, "delta_GCV_rel_eff", "NA")); support = floator(get(r, "support_2D_ell_nm", "NA"))
    g = guard(n, robust, amb, runner, dgcv, P)
    return (V0=midpoint(g, n, support, dgcv, P), V1=midpoint(n, n, support, dgcv, P), V2=g, V3=n)
end

"Strict-majority count of a list, or nothing."
function strict_majority(ns)
    t = Dict{Int,Int}()
    for x in ns
        t[x] = get(t, x, 0) + 1
    end
    n, k = first(sort(collect(t); by=p -> (-last(p), first(p))))
    return 2k > length(ns) ? n : nothing
end
vote(ns) = (m = strict_majority(ns); m === nothing ? sort(ns)[(length(ns) + 1) ÷ 2] : m)

function raw_index(dir)
    index = Dict{String,String}()
    for (d, _, names) in walkdir(dir; follow_symlinks=true), n in names
        endswith(lowercase(n), ".sxm") && (index[n] = joinpath(d, n))
    end
    return index
end

const VARIANTS = (:V0, :V1, :V2, :V3)

function replay_main(args=ARGS)
    o = replay_options(args)
    o === nothing && return println("replay_count_policy.jl --run RUN_DIR --config COUNT_TOML --outdir NEW_DIR [--data-dir RAW_DIR]")
    P = policy_settings(o["--config"])
    _, summary = read_table(joinpath(o["--run"], "counting_summary.tsv"))
    rows = Dict(basename(r["filepath"]) => r for r in summary if get(r, "status", "ok") == "ok")
    per = Dict(f => variants(r, P) for (f, r) in rows)
    _, cons = read_table(joinpath(o["--run"], "consensus", "consensus.tsv"))
    track = Dict(r["file"] => r["track"] for r in cons if haskey(per, r["file"]))
    nfinal = Dict(r["file"] => r["N_final"] for r in cons)
    members = Dict{String,Vector{String}}()
    for (f, t) in track
        push!(get!(members, t, String[]), f)
    end
    if haskey(o, "--data-dir")
        raw = raw_index(o["--data-dir"])
        for fs in values(members)
            length(fs) >= 4 || continue
            times = Dict(f => scan_frame(read_sxm(raw[f])).acquired for f in fs)
            sort!(fs; by=f -> (times[f], f))
        end
    else
        foreach(sort!, values(members))
    end
    consensus = Dict{Symbol,Dict{String,Int}}(v => Dict{String,Int}() for v in VARIANTS)
    for (t, fs) in members, v in VARIANTS
        m = length(fs) >= 3 ? strict_majority([getproperty(per[f], v) for f in fs]) : nothing
        for f in fs
            consensus[v][f] = m === nothing ? getproperty(per[f], v) : m
        end
    end
    out = o["--outdir"]; mkpath(out)
    files = sort(collect(keys(track)))
    open(joinpath(out, "variant_counts.tsv"), "w") do io
        println(io, join(vcat(["file", "track", "track_size", "N_selected", "N_final"], [string(v) for v in VARIANTS],
                              [string(v) * "_consensus" for v in VARIANTS]), '\t'))
        for f in files
            println(io, join(vcat([f, track[f], length(members[track[f]]), rows[f]["N_selected"], get(nfinal, f, "NA")],
                                  [getproperty(per[f], v) for v in VARIANTS], [consensus[v][f] for v in VARIANTS]), '\t'))
        end
    end
    halves = Tuple{String,String,Int,Int}[]
    for (t, fs) in sort(collect(members); by=first)
        length(fs) >= 4 || continue
        a, b = fs[1:2:end], fs[2:2:end]
        for v in VARIANTS
            push!(halves, (t, string(v), vote([getproperty(per[f], v) for f in a]), vote([getproperty(per[f], v) for f in b])))
        end
        push!(halves, (t, "own", vote([parse(Int, rows[f]["N_selected"]) for f in a]), vote([parse(Int, rows[f]["N_selected"]) for f in b])))
    end
    open(joinpath(out, "split_half.tsv"), "w") do io
        println(io, "track\tvariant\tN_odd\tN_even")
        foreach(h -> println(io, join(h, '\t')), halves)
    end
    md = IOBuffer()
    match0 = count(f -> per[f].V0 == parse(Int, rows[f]["N_selected"]), files)
    println(md, "# Offline replay of the per-scan count layers (campaign K)\n")
    println(md, "Run `$(o["--run"])`, config `$(o["--config"])` (n_max $(P.n_max)). $(length(files)) scans, $(count(fs -> length(fs) >= 3, values(members))) tracks of at least 3 scans.\n")
    println(md, "Replay check: V0 equals N_selected for $match0 of $(length(files)) scans (it must be all of them on a run of the support_midpoint_hybrid policy).\n")
    println(md, "| variant | per-scan counts changed against V0 | consensus counts changed against V0 |\n|---|---:|---:|")
    for v in VARIANTS
        println(md, "| $v | $(count(f -> getproperty(per[f], v) != per[f].V0, files)) | $(count(f -> consensus[v][f] != consensus[:V0][f], files)) |")
    end
    println(md, "\n## Split-half agreement (tracks of at least 4 scans; votes: strict majority, else lower median)\n")
    println(md, "| variant | tracks | halves agree | mean abs difference | difference > 2 |\n|---|---:|---:|---:|---:|")
    for v in vcat([string(v) for v in VARIANTS], ["own"])
        d = [abs(h[3] - h[4]) for h in halves if h[2] == v]
        isempty(d) && continue
        println(md, @sprintf("| %s | %d | %d (%.0f %%) | %.2f | %d |", v, length(d), count(==(0), d), 100 * count(==(0), d) / length(d), mean(d), count(>(2), d)))
    end
    println(md, "\n## Count range across the scans of a track (tracks of at least 3 scans)\n")
    println(md, "| variant | range 0 | range ≤ 1 | range > 2 |\n|---|---:|---:|---:|")
    for v in VARIANTS
        r = [maximum(getproperty(per[f], v) for f in fs) - minimum(getproperty(per[f], v) for f in fs) for fs in values(members) if length(fs) >= 3]
        println(md, "| $v | $(count(==(0), r)) | $(count(<=(1), r)) | $(count(>(2), r)) |")
    end
    write(joinpath(out, "report.md"), take!(md))
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && replay_main()

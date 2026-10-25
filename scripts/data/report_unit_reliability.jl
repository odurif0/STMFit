#!/usr/bin/env julia
# Benchmark-free reliability of per-unit signals over repeated scans of one
# molecule, for one or more runs of run_molecule_consensus_chitosan.jl.
# Label-free: scan_geometry.tsv, predictions.tsv (physical lobes), consensus.tsv
# and, for acquisition times, the SXM headers. No benchmark label is read.
#
# Usage:
#   julia --project=. scripts/data/report_unit_reliability.jl \
#       --run NAME=RUN_DIR [--run NAME=RUN_DIR ...] [--raw NAME=RAW_DIR ...] \
#       [--control NAME] --outdir NEW_DIR [--cap 8] [--dmin 2] [--draws 1000] [--seed 20261025]
#
# --raw gives the folder of the raw scans of a run (acquisition times, used for
# the order of the scans and the time-lag table). --control names the run that
# serves as negative control for sequence (default: the first run): a run whose
# molecules all have the same sequence, so that unit-level variation which
# repeats across scans cannot be sequence once the distance-to-end profile is
# removed. Every other run is compared with it by the ratio of between-lobe
# variances. See docs/src/journal.md (campaign R, 2026-10-25) for the
# definitions and the decision rules, fixed before this script was run.
using STMFit.UnitReliability, STMFit.MoleculeConsensus
using STMSXMIO: read_sxm
using Printf, Random, Statistics, Dates

function reliability_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,Any}("--run" => Pair{String,String}[], "--raw" => Pair{String,String}[],
                         "--cap" => "8", "--dmin" => "2", "--draws" => "1000", "--seed" => "20261025"); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--run", "--raw", "--control", "--outdir", "--cap", "--dmin", "--draws", "--seed") && i < length(args) ||
            error("Unknown option or missing value: $k")
        if k in ("--run", "--raw")
            name, path = split(args[i+1], "="; limit=2)
            push!(o[k], String(name) => String(path))
        else
            o[k] = args[i+1]
        end
        i += 2
    end
    isempty(o["--run"]) && error("Required: --run NAME=RUN_DIR")
    haskey(o, "--outdir") || error("Required: --outdir")
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

"Acquisition time of every scan of the tracks of at least 3 scans."
function acquisition_times(ru, raw_dir)
    raw = raw_index(raw_dir); times = Dict{String,DateTime}()
    for fs in values(ru.tracks), f in fs
        length(fs) >= 3 && haskey(raw, f) && (times[f] = scan_frame(read_sxm(raw[f])).acquired)
    end
    return times
end

modal(v) = (c = Dict{Int,Int}(); for x in v; c[x] = get(c, x, 0) + 1; end; first(sort(collect(c); by=p -> (-last(p), first(p)))).first)

const LENGTH_GROUPS = (("all", 0, 99), ("N<=6", 0, 6), ("N=7-9", 7, 9), ("N>=10", 10, 99))

"Per-track moments of the groups, for the tracks whose modal final count lies in [lo, hi]."
function group_moments(ru, kind; cap, dmin, lo, hi, profile)
    tg = track_groups(ru, kind; cap, dmin, profile)
    out = Moments[]
    for (t, g) in tg
        n = modal([ru.n_final[f] for f in ru.tracks[t]])
        lo <= n <= hi && push!(out, moments(g))
    end
    return out
end

fmt(x; d=3) = isnan(x) ? "NA" : @sprintf("%.*f", d, x)
sd(v) = isnan(v) ? NaN : sqrt(max(v, 0.0))

function report_main(args=ARGS)
    o = reliability_options(args)
    o === nothing && return println("report_unit_reliability.jl --run NAME=RUN_DIR [--run ...] [--raw NAME=RAW_DIR ...] [--control NAME] --outdir NEW_DIR")
    cap, dmin, draws, seed = parse(Int, o["--cap"]), parse(Int, o["--dmin"]), parse(Int, o["--draws"]), parse(Int, o["--seed"])
    raw = Dict(o["--raw"])
    runs = Dict{String,RunUnits}(); times = Dict{String,Dict{String,DateTime}}()
    for (name, path) in o["--run"]
        ru = load_run_units(path; name)
        if haskey(raw, name)
            times[name] = acquisition_times(ru, raw[name])
            ru = load_run_units(path; name, times=times[name])
        end
        runs[name] = ru
    end
    names = [n for (n, _) in o["--run"]]
    control = get(o, "--control", first(names))
    haskey(runs, control) || error("--control $control is not one of the runs")
    rng = MersenneTwister(seed)
    mkpath(o["--outdir"])
    tsv = open(joinpath(o["--outdir"], "variance_components.tsv"), "w")
    println(tsv, join(["run", "signal", "dmin", "group", "tracks", "lobes", "observations", "sd_between", "sd_within", "icc", "icc_lo", "icc_hi"], '\t'))
    md = IOBuffer()
    println(md, "# Reliability of per-unit signals\n")
    println(md, "Residual = log signal minus scan median minus the mean at the distance class to the nearest end. Physical lobes seen in at least 3 scans of a track (at most $cap scans per track), d ≥ $dmin unless stated. Intervals: bootstrap over tracks ($draws draws, seed $seed). Negative control for sequence: `$control`.\n")
    println(md, "| run | tracks with ≥ 3 scans | scans in them |\n|---|---:|---:|")
    for n in names
        ru = runs[n]; big = [t for (t, fs) in ru.tracks if length(fs) >= 3]
        println(md, "| $n | $(length(big)) | $(sum(length(ru.tracks[t]) for t in big)) |")
    end
    ratios = String[]
    for kind in SIGNALS
        println(md, "\n## Signal: $kind\n")
        println(md, "| run | group | tracks | lobes | σ between | σ within | ICC | 95% interval |\n|---|---|---:|---:|---:|---:|---:|---|")
        tm = Dict{Tuple{String,Int},Vector{Moments}}()
        for n in names
            ru = runs[n]; profile = distance_profile(ru, kind)
            variants = vcat([(dmin, g) for g in LENGTH_GROUPS], [(dm, LENGTH_GROUPS[1]) for dm in (1, 2, 3) if dm != dmin])
            for (dm, (gname, lo, hi)) in variants
                v = group_moments(ru, kind; cap, dmin=dm, lo, hi, profile)
                length(v) >= 3 || continue
                m = sum(v; init=Moments()); b, w, i = variance_components(m)
                _, _, ib = bootstrap_components(v; draws, rng)
                lo_i, hi_i = isempty(ib) ? (NaN, NaN) : (quantile(ib, 0.025), quantile(ib, 0.975))
                label = dm == dmin ? gname : "$gname (d ≥ $dm)"
                println(md, "| $n | $label | $(length(v)) | $(m.groups) | $(fmt(sd(b))) | $(fmt(sqrt(w))) | $(fmt(i, d=2)) | $(fmt(lo_i, d=2)) to $(fmt(hi_i, d=2)) |")
                println(tsv, join([n, kind, dm, gname, length(v), m.groups, m.n, sd(b), sqrt(w), i, lo_i, hi_i], '\t'))
                dm == dmin && gname == "all" && (tm[(n, dm)] = v)
            end
        end
        for n in names
            n == control && continue
            (haskey(tm, (n, dmin)) && haskey(tm, (control, dmin))) || continue
            r, lo_r, hi_r = ratio_interval(tm[(n, dmin)], tm[(control, dmin)]; draws, rng)
            push!(ratios, "| $n | $kind | $(fmt(r, d=2)) | $(fmt(lo_r, d=2)) to $(fmt(hi_r, d=2)) | $(lo_r > 1 ? "above 1: extra structure" : "no evidence of extra structure") |")
        end
    end
    println(md, "\n## Between-lobe variance against the negative control\n")
    println(md, "Ratio σ²_between(run) / σ²_between($control), all tracks, d ≥ $dmin. The rule fixed in advance asks for an interval above 1.\n")
    println(md, "| run | signal | ratio | 95% interval | reading |\n|---|---|---:|---|---|")
    foreach(r -> println(md, r), ratios)
    # time lag
    lagrows = String[]
    for n in names
        haskey(times, n) || continue
        ru = runs[n]; profile = distance_profile(ru, :amplitude)
        ps = pair_correlations(ru, :amplitude, times[n]; dmin, profile)
        for (lo, hi, r, l, h, nt) in lag_summary(ps, [(0.0, 3.0), (3.0, 10.0), (10.0, 30.0), (30.0, Inf)]; draws, rng)
            push!(lagrows, "| $n | $(Int(lo))–$(isinf(hi) ? "∞" : string(Int(hi))) | $(fmt(r, d=2)) | $(fmt(l, d=2)) to $(fmt(h, d=2)) | $nt |")
        end
    end
    if !isempty(lagrows)
        println(md, "\n## Correlation of the residual profile between two scans of a molecule, by time lag (amplitude)\n")
        println(md, "| run | lag (min) | mean r | 95% interval | tracks |\n|---|---|---:|---|---:|")
        foreach(r -> println(md, r), lagrows)
    end
    # order along the chain
    println(md, "\n## Order along the chain of the persistent profile (amplitude)\n")
    println(md, "Mean autocorrelation of each molecule's mean residual profile, against the same profiles shuffled along the chain.\n")
    println(md, "| run | lag | molecules | observed | shuffled mean | shuffled 95% range |\n|---|---:|---:|---:|---:|---|")
    for n in names
        ru = runs[n]; P = collect(values(persistent_profiles(ru, :amplitude; dmin, profile=distance_profile(ru, :amplitude))))
        length(P) >= 3 || continue
        for lag in (1, 2)
            obs, null = autocorrelation_null(P, lag; draws=2000, rng)
            println(md, "| $n | $lag | $(length(P)) | $(fmt(obs)) | $(fmt(mean(null))) | $(fmt(quantile(null, 0.025))) to $(fmt(quantile(null, 0.975))) |")
        end
    end
    # counts
    println(md, "\n## Reproducibility of the per-scan counts\n")
    println(md, "Range (max − min) of the counts of the scans of one molecule, tracks of at least 2 scans.\n")
    println(md, "| run | tracks | all scans agree | range ≤ 1 | range ≤ 2 | range > 2 | largest range |\n|---|---:|---:|---:|---:|---:|---:|")
    for n in names
        cd = count_dispersion(runs[n]); r = [x[3] for x in cd]
        isempty(r) && continue
        println(md, "| $n | $(length(r)) | $(count(==(0), r)) | $(count(<=(1), r)) | $(count(<=(2), r)) | $(count(>(2), r)) | $(maximum(r)) |")
    end
    close(tsv)
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && report_main()

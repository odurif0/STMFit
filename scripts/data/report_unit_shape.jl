#!/usr/bin/env julia
# Shape of the units from the lobe patches of runs of
# run_molecule_consensus_chitosan.jl, over repeated scans of one molecule
# (campaign J, stage 2a; docs/src/journal.md, 2026-10-26, for the descriptors and
# the decision rule, fixed before this script was run). Label-free: patches,
# physical lobes and tracks only; no benchmark label is read.
#
# Usage:
#   julia --project=. scripts/data/report_unit_shape.jl --run NAME=RUN_DIR --run NAME=RUN_DIR \
#       --control NAME --outdir NEW_DIR [--cap 8] [--dmin 2] [--draws 1000] [--seed 20261026]
#
# Patches: RUN_DIR/assignment/patches_fwd17.tsv (raw_p columns, 17 x 17, u-major,
# t along the chain, normalised per patch). Within a track, the patches of scans
# whose physical lobes run against the lobe index are turned by 180°, so a track
# has one orientation. For the end analysis and the averaged patches, t points
# away from the nearest chain end (the middle lobe of an odd chain keeps its
# orientation). --control names the negative control for sequence (a run whose
# molecules share one sequence); every other run is compared with it.
# Outputs: report.md, variance_components.tsv, mean_patches.tsv, ends.tsv.
using STMFit.UnitReliability
using Printf, Random, Statistics

function shape_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,Any}("--run" => Pair{String,String}[], "--cap" => "8", "--dmin" => "2",
                         "--draws" => "1000", "--seed" => "20261026"); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--run", "--control", "--outdir", "--cap", "--dmin", "--draws", "--seed") && i < length(args) ||
            error("Unknown option or missing value: $k")
        if k == "--run"
            name, path = split(args[i+1], "="; limit=2)
            push!(o[k], String(name) => String(path))
        else
            o[k] = args[i+1]
        end
        i += 2
    end
    length(o["--run"]) >= 2 || error("Required: at least two --run NAME=RUN_DIR")
    haskey(o, "--outdir") || error("Required: --outdir")
    haskey(o, "--control") || error("Required: --control NAME")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Raw patches of a run: (file, lobe) => 289 values (NaN for NA)."
function read_patches(run_dir)
    lines = readlines(joinpath(run_dir, "assignment", "patches_fwd17.tsv"))
    h = split(lines[1], '\t')
    cols = findall(c -> startswith(c, "raw_p"), h)
    length(cols) == 289 || error("Expected 289 raw_p columns in $run_dir")
    fi, li = findfirst(==("file"), h), findfirst(==("lobe"), h)
    out = Dict{Tuple{String,Int},Vector{Float64}}()
    for l in lines[2:end]
        isempty(strip(l)) && continue
        v = split(l, '\t')
        out[(String(v[fi]), parse(Int, v[li]))] = [x == "NA" ? NaN : parse(Float64, x) for x in v[cols]]
    end
    return out
end

"Patch of lobe `l` of a scan of `n` lobes, turned so that t points away from the nearest end."
inward(P, l, n) = 2l.k > n + 1 ? rotate180(P) : P

const DCLASS = 4

struct RunShape
    name::String
    ru::RunUnits
    patches::Dict{Tuple{String,Int},Vector{Float64}}
    used::Dict{String,Vector{String}}           # track => used scans (at least 3 scans, capped)
end

function run_shape(name, path, cap)
    ru = load_run_units(path; name)
    patches = read_patches(path)
    used = Dict(t => capped_files(fs, cap) for (t, fs) in ru.tracks if length(fs) >= 3)
    return RunShape(name, ru, patches, used)
end

"Descriptor values of every lobe of a scan in the track orientation (NaN where unavailable)."
function scan_values(rs::RunShape, f, kind)
    lobes = rs.ru.lobes[f]; o = scan_orientation(lobes)
    return [begin
                P = get(rs.patches, (f, l.k), nothing)
                (P === nothing || o == 0) ? NaN : getproperty(patch_descriptors(o < 0 ? rotate180(P) : P), kind)
            end for l in lobes]
end

"Residual groups by physical lobe per track, as campaign R: minus scan median, minus the distance-class mean."
function shape_groups(rs::RunShape, kind; dmin)
    vals = Dict(f => scan_values(rs, f, kind) for fs in values(rs.used) for f in fs)
    centred = Dict(f => (ok = filter(!isnan, v); isempty(ok) ? v : v .- median(ok)) for (f, v) in vals)
    sums = zeros(DCLASS); cnt = zeros(Int, DCLASS)                   # every track counts once
    for fs in values(rs.used)
        s = zeros(DCLASS); c = zeros(Int, DCLASS)
        for f in fs, (l, v) in zip(rs.ru.lobes[f], centred[f])
            isnan(v) && continue
            j = min(l.d, DCLASS); s[j] += v; c[j] += 1
        end
        for j in 1:DCLASS
            c[j] > 0 && (sums[j] += s[j] / c[j]; cnt[j] += 1)
        end
    end
    profile = [cnt[j] > 0 ? sums[j] / cnt[j] : NaN for j in 1:DCLASS]
    out = Moments[]
    for fs in values(rs.used)
        groups = Dict{String,Vector{Float64}}()
        for f in fs, (l, v) in zip(rs.ru.lobes[f], centred[f])
            (isnan(v) || l.d < dmin || isempty(l.physical)) && continue
            push!(get!(groups, l.physical, Float64[]), v - profile[min(l.d, DCLASS)])
        end
        g = [x for x in values(groups) if length(x) >= 3]
        isempty(g) || push!(out, moments(g))
    end
    return out, profile
end

"Per track: mean end asymmetry of its end lobes and summed inward patches by distance class."
function end_and_means(rs::RunShape)
    E = Float64[]
    sums = Dict(j => zeros(289) for j in 1:DCLASS); cnt = Dict(j => 0 for j in 1:DCLASS)
    for fs in values(rs.used)
        e = Float64[]
        for f in fs
            lobes = rs.ru.lobes[f]; n = length(lobes)
            for l in lobes
                P = get(rs.patches, (f, l.k), nothing)
                (P === nothing || !all(isfinite, P)) && continue
                Q = inward(P, l, n)
                j = min(l.d, DCLASS); sums[j] .+= Q; cnt[j] += 1
                l.d == 1 && push!(e, end_asymmetry(Q))
            end
        end
        isempty(e) || push!(E, mean(e))
    end
    return E, Dict(j => cnt[j] > 0 ? sums[j] ./ cnt[j] : fill(NaN, 289) for j in 1:DCLASS), cnt
end

fmt(x; d=3) = isnan(x) ? "NA" : @sprintf("%.*f", d, x)
sdv(b) = isnan(b) ? NaN : sqrt(max(b, 0.0))

function shape_main(args=ARGS)
    o = shape_options(args)
    o === nothing && return println("report_unit_shape.jl --run NAME=RUN_DIR --run NAME=RUN_DIR --control NAME --outdir NEW_DIR")
    cap, dmin, draws = parse(Int, o["--cap"]), parse(Int, o["--dmin"]), parse(Int, o["--draws"])
    rng = MersenneTwister(parse(Int, o["--seed"]))
    runs = [run_shape(n, p, cap) for (n, p) in o["--run"]]
    ctrl = findfirst(r -> r.name == o["--control"], runs)
    ctrl === nothing && error("--control $(o["--control"]) is not one of the runs")
    mkpath(o["--outdir"])
    md = IOBuffer()
    println(md, "# Shape of the units from averaged lobe patches (campaign J, stage 2a)\n")
    println(md, "Descriptors of the normalised forward patches (17 × 17, ±0.32 nm), residuals minus scan median and distance-class mean, physical lobes with d ≥ $dmin seen in at least 3 of at most $cap scans per track. Negative control for sequence: `$(runs[ctrl].name)`. Intervals: bootstrap over tracks ($draws draws).\n")
    println(md, "| run | tracks used | scans used |\n|---|---:|---:|")
    for r in runs
        println(md, "| $(r.name) | $(length(r.used)) | $(sum(length, values(r.used))) |")
    end
    tsv = open(joinpath(o["--outdir"], "variance_components.tsv"), "w")
    println(tsv, join(["run", "descriptor", "tracks", "lobes", "sd_between", "sd_within", "icc", "profile_d1", "profile_d2", "profile_d3", "profile_d4"], '\t'))
    println(md, "\n## Between-lobe variance of each descriptor\n")
    println(md, "| run | descriptor | tracks | lobes | σ between | σ within | ICC | d = 2 minus d = 3 | its d′ |\n|---|---|---:|---:|---:|---:|---:|---:|---:|")
    tm = Dict{Tuple{String,Symbol},Vector{Moments}}()
    for kind in SHAPE_DESCRIPTORS, r in runs
        m, prof = shape_groups(r, kind; dmin)
        tm[(r.name, kind)] = m
        b, w, icc = variance_components(sum(m; init=Moments()))
        Δ = prof[2] - prof[3]
        println(md, "| $(r.name) | $kind | $(length(m)) | $(sum(x -> x.groups, m; init=0)) | $(fmt(sdv(b), d=4)) | $(fmt(sqrt(w), d=4)) | $(fmt(icc, d=2)) | $(fmt(Δ, d=4)) | $(fmt(abs(Δ) / sdv(b), d=2)) |")
        println(tsv, join([r.name, kind, length(m), sum(x -> x.groups, m; init=0), sdv(b), sqrt(w), icc, prof...], '\t'))
    end
    close(tsv)
    println(md, "\n## Rule (a): between-lobe variance against the negative control (99% intervals)\n")
    println(md, "| run | descriptor | ratio | 99% interval | reading |\n|---|---|---:|---|---|")
    hits_a = String[]
    for r in runs, kind in SHAPE_DESCRIPTORS
        r.name == runs[ctrl].name && continue
        q, lo, hi = ratio_interval(tm[(r.name, kind)], tm[(runs[ctrl].name, kind)]; draws, rng, level=0.99)
        lo > 1 && push!(hits_a, "$(r.name) $kind")
        println(md, "| $(r.name) | $kind | $(fmt(q, d=2)) | $(fmt(lo, d=2)) to $(fmt(hi, d=2)) | $(lo > 1 ? "above 1" : "not above 1") |")
    end
    println(md, "\n## Rule (b): chain ends (end lobes oriented inward)\n")
    println(md, "E = mean of the patch over the outer part (t ≤ −0.12 nm, |u| ≤ 0.16 nm) minus the inner part (t ≥ 0.12 nm); one value per track (mean over its end lobes), then the mean over tracks.\n")
    ends = Dict(r.name => end_and_means(r) for r in runs)
    E0 = ends[runs[ctrl].name][1]
    println(md, "| run | tracks | mean E | difference with the control | 99% interval | half of the control's abs E | reading |\n|---|---:|---:|---:|---|---:|---|")
    hits_b = String[]
    open(joinpath(o["--outdir"], "ends.tsv"), "w") do io
        println(io, "run\ttrack_mean_E")
        for r in runs, e in ends[r.name][1]
            println(io, "$(r.name)\t$e")
        end
    end
    for r in runs
        E = ends[r.name][1]
        if r.name == runs[ctrl].name
            println(md, "| $(r.name) | $(length(E)) | $(fmt(mean(E))) | | | | control |")
            continue
        end
        Δ = mean(E) - mean(E0)
        boot = [mean(E[rand(rng, 1:length(E), length(E))]) - mean(E0[rand(rng, 1:length(E0), length(E0))]) for _ in 1:draws]
        lo, hi = quantile(boot, 0.005), quantile(boot, 0.995)
        clear = (lo > 0 || hi < 0) && abs(Δ) >= abs(mean(E0)) / 2
        clear && push!(hits_b, r.name)
        println(md, "| $(r.name) | $(length(E)) | $(fmt(mean(E))) | $(fmt(Δ)) | $(fmt(lo)) to $(fmt(hi)) | $(fmt(abs(mean(E0)) / 2)) | $(clear ? "clearly different" : "not clearly different") |")
    end
    open(joinpath(o["--outdir"], "mean_patches.tsv"), "w") do io
        println(io, join(vcat(["run", "class", "patches"], [@sprintf("p%03d", i) for i in 1:289]), '\t'))
        for r in runs, j in 1:DCLASS
            _, means, cnt = ends[r.name]
            cnt[j] > 0 && println(io, join(vcat([r.name, j == DCLASS ? "d>=$DCLASS" : "d=$j", cnt[j]], [@sprintf("%.5f", v) for v in means[j]]), '\t'))
        end
    end
    build = !isempty(hits_a) || !isempty(hits_b)
    println(md, "\n## Decision (rule fixed on 2026-10-26)\n")
    println(md, "(a) $(isempty(hits_a) ? "no descriptor above the control" : "above the control: " * join(hits_a, ", ")); (b) $(isempty(hits_b) ? "ends not clearly different" : "ends clearly different: " * join(hits_b, ", ")). ",
            build ? "**The joint fit is built.**" : "**The joint fit is not built; the campaign stops.**")
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && shape_main()

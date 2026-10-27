#!/usr/bin/env julia
# Two unit types on an application run, read without position (campaign T;
# docs/src/journal.md 2026-10-27, for the features, the model and the criteria,
# fixed before this script was run). Label-free: assignment features, lobe
# patches, physical lobes and tracks of runs of run_molecule_consensus_chitosan.jl.
#
# Usage:
#   julia --project=. scripts/data/type_units.jl --run RUN_DIR --control CONTROL_RUN_DIR \
#       --outdir NEW_DIR [--cap 16] [--seed 20261027]
#
# The same pipeline runs on the application run and on the fixed-sequence
# control. Outputs: units.tsv (both runs), sequences.tsv (application run,
# molecules whose per-scan counts vary by at most one unit), report.md with the
# criteria C1–C4 and, when they all hold, the sample-level predictions.
using STMFit.UnitTyping, STMFit.UnitReliability, STMFit.ReconstructedUnitAssignment
using Printf, Statistics, SHA

const FEATURES = ["log_amplitude", "amp_rel", "amp_prominence", "integrated_prominence", "amp_neighbor_ratio",
                  "split_log_skew", "patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd",
                  "width_t", "width_u", "central_fraction"]

function typing_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--cap" => "16", "--seed" => "20261027"); i = 1
    while i <= length(args)
        args[i] in ("--run", "--control", "--outdir", "--cap", "--seed") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--run", "--control", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

num(s) = (v = tryparse(Float64, s); v === nothing ? NaN : v)

"Raw patches of a run: (file, lobe) => values."
function read_patches(run)
    lines = readlines(joinpath(run, "assignment", "patches_fwd17.tsv"))
    h = split(lines[1], '\t'); cols = findall(c -> startswith(c, "raw_p"), h)
    out = Dict{Tuple{String,Int},Vector{Float64}}()
    for l in lines[2:end]
        v = split(l, '\t')
        out[(String(v[1]), parse(Int, v[2]))] = [x == "NA" ? NaN : parse(Float64, x) for x in v[cols]]
    end
    return out
end

"Units of a run: per physical unit, its per-scan residual feature vectors, halves and distance class."
function run_units(run; cap)
    _, feat = read_table(joinpath(run, "assignment", "features_predictor.tsv"))
    _, pred = read_table(joinpath(run, "predictions.tsv"))
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    physical = Dict((r["file"], parse(Int, r["lobe"])) => r["physical_lobe"] for r in pred)
    patches = read_patches(run)
    byfile = Dict{String,Vector{Dict{String,String}}}()
    for r in feat
        push!(get!(byfile, r["file"], Dict{String,String}[]), r)
    end
    track = Dict(r["file"] => r["track"] for r in cons); nscan = Dict(r["file"] => parse(Int, r["N_scan"]) for r in cons)
    members = Dict{String,Vector{String}}()
    for f in keys(byfile)
        push!(get!(members, track[f], String[]), f)
    end
    # per-scan robust z of every feature, in track orientation
    obs = Dict{String,Vector{Tuple{String,Int,Int,Vector{Float64}}}}()      # file => (physical, k, d, z)
    for (f, rows) in byfile
        sort!(rows; by=r -> parse(Int, r["lobe"])); n = length(rows)
        lobes = [Lobe(parse(Int, r["lobe"]), min(parse(Int, r["lobe"]), n + 1 - parse(Int, r["lobe"])),
                      get(physical, (f, parse(Int, r["lobe"])), ""), num(r["amplitude"]), 1.0, 1.0) for r in rows]
        o = scan_orientation(lobes)
        M = fill(NaN, n, length(FEATURES))
        for (i, r) in enumerate(rows)
            a = num(r["amplitude"]); M[i, 1] = a > 0 ? log(a) : NaN
            for (j, c) in enumerate(FEATURES[2:9])
                M[i, j+1] = num(get(r, c, "NA"))
            end
            P = get(patches, (f, lobes[i].k), nothing)
            if P !== nothing && o != 0
                d = patch_descriptors(o < 0 ? rotate180(P) : P)
                M[i, 10], M[i, 11], M[i, 12] = d.width_t, d.width_u, d.central_fraction
            end
        end
        for j in axes(M, 2)
            M[:, j] = robust_z(M[:, j])
        end
        obs[f] = [(lobes[i].physical, lobes[i].k, lobes[i].d, M[i, :]) for i in 1:n]
    end
    # distance-class means, every track once
    p = length(FEATURES); prof = zeros(4, p)
    for j in 1:4, c in 1:p
        tm = Float64[]
        for fs in values(members)
            v = [x[4][c] for f in fs for x in obs[f] if min(x[3], 4) == j && isfinite(x[4][c])]
            isempty(v) || push!(tm, mean(v))
        end
        prof[j, c] = isempty(tm) ? 0.0 : mean(tm)
    end
    units = Dict{Tuple{String,String},Any}()
    for (t, fs) in members
        length(fs) >= 4 || continue
        used = capped_files(sort(fs), cap)
        for (h, f) in enumerate(used), (ph, k, d, z) in obs[f]
            isempty(ph) && continue
            u = get!(units, (t, ph), (odd=Vector{Float64}[], even=Vector{Float64}[], d=Int[]))
            push!(isodd(h) ? u.odd : u.even, z .- prof[min(d, 4), :]); push!(u.d, d)
        end
    end
    return units, members, nscan
end

meanfinite(V) = [(x = filter(isfinite, [v[c] for v in V]); isempty(x) ? NaN : mean(x)) for c in eachindex(V[1])]
modal(v) = first(sort(unique(v); by=x -> (-count(==(x), v), x)))
physical_index(p) = parse(Int, split(p, '_')[end])

"Typing of one run: unit table, mixtures and criteria values."
function type_run(units; seed)
    keys_ = sort([k for (k, u) in units if length(u.odd) >= 2 && length(u.even) >= 2])
    A = [meanfinite(vcat(units[k].odd, units[k].even)) for k in keys_]
    Od = [meanfinite(units[k].odd) for k in keys_]; Ev = [meanfinite(units[k].even) for k in keys_]
    ok = [all(isfinite, A[i]) && all(isfinite, Od[i]) && all(isfinite, Ev[i]) for i in eachindex(keys_)]
    keys_, A, Od, Ev = keys_[ok], reduce(hcat, A[ok])', reduce(hcat, Od[ok])', reduce(hcat, Ev[ok])'
    μ = mean(A; dims=1); σ = std(A; dims=1); σ[σ .== 0] .= 1
    Z, Zo, Ze = (A .- μ) ./ σ, (Od .- μ) ./ σ, (Ev .- μ) ./ σ
    V, vals = reproducible_subspace(Zo, Ze)
    S, So, Se = Z * V, Zo * V, Ze * V
    M1 = fit_mixture(S, 1; seed); M2 = fit_mixture(S, 2; seed)
    dbic = mixture_bic(M1, size(S, 1)) - mixture_bic(M2, size(S, 1))
    P = posterior(M2, S)
    bright = argmax([sum(P[:, j] .* Z[:, 1]) / sum(P[:, j]) for j in 1:2])        # higher log-amplitude residual
    type1 = P[:, bright] .> 0.5
    to = posterior(M2, So)[:, bright] .> 0.5; te = posterior(M2, Se)[:, bright] .> 0.5
    d = [min(modal(units[k].d), 4) for k in keys_]
    shares = [count(i -> d[i] == j && type1[i], eachindex(d)) / max(1, count(==(j), d)) for j in 1:4]
    return (keys=keys_, Z=Z, V=V, vals=vals, M2=M2, dbic=dbic, post=P[:, bright], type1=type1, type_odd=to, type_even=te,
            kappa=cohen_kappa(to, te), d=d, shares=shares, nscans=[length(units[k].odd) + length(units[k].even) for k in keys_])
end

fmt(x; d=2) = isnan(x) ? "NA" : @sprintf("%.*f", d, x)

function typing_main(args=ARGS)
    o = typing_options(args)
    o === nothing && return println("type_units.jl --run RUN_DIR --control CONTROL_RUN_DIR --outdir NEW_DIR [--cap 16] [--seed 20261027]")
    cap, seed = parse(Int, o["--cap"]), parse(Int, o["--seed"])
    ua, membersa, nscana = run_units(o["--run"]; cap)
    uc, _, _ = run_units(o["--control"]; cap)
    ra = type_run(ua; seed); rc = type_run(uc; seed)
    out = o["--outdir"]; mkpath(out)
    open(joinpath(out, "units.tsv"), "w") do io
        println(io, "run\ttrack\tphysical_lobe\td_class\tscans\ttype\tposterior_type1\ttype_odd\ttype_even")
        for (name, r) in (("application", ra), ("control", rc)), i in eachindex(r.keys)
            println(io, join([name, r.keys[i][1], r.keys[i][2], r.d[i], r.nscans[i], r.type1[i] ? 1 : 0, fmt(r.post[i]; d=4),
                              r.type_odd[i] ? 1 : 0, r.type_even[i] ? 1 : 0], '\t'))
        end
    end
    c1 = ra.dbic >= 6; c2 = ra.kappa >= 0.5; c3 = rc.dbic < 6; c4 = maximum(ra.shares) - minimum(ra.shares) < 0.15
    md = IOBuffer()
    println(md, "# Two unit types without position (campaign T)\n")
    println(md, "Application `$(o["--run"])`, control `$(o["--control"])`. Features: $(join(FEATURES, ", ")); robust z within scan, distance-class means removed, unit means over odd, even and all scans (at most $cap scans per molecule).\n")
    println(md, "| | application | control |\n|---|---:|---:|")
    println(md, "| units (≥ 2 scans in each half) | $(length(ra.keys)) | $(length(rc.keys)) |")
    println(md, "| molecules | $(length(unique(first.(ra.keys)))) | $(length(unique(first.(rc.keys)))) |")
    println(md, "| reproducible dimensions kept | $(size(ra.V, 2)) | $(size(rc.V, 2)) |")
    println(md, "| ΔBIC (1 minus 2 components) | $(fmt(ra.dbic; d=1)) | $(fmt(rc.dbic; d=1)) |")
    println(md, "| type-1 share (brighter family) | $(fmt(100mean(ra.type1); d=0)) % | $(fmt(100mean(rc.type1); d=0)) % |")
    println(md, "| κ, odd against even scans | $(fmt(ra.kappa)) | $(fmt(rc.kappa)) |")
    println(md, "| type-1 share at d = 1 / 2 / 3 / ≥ 4 | $(join([fmt(100s; d=0) for s in ra.shares], " / ")) % | $(join([fmt(100s; d=0) for s in rc.shares], " / ")) % |\n")
    println(md, "| criterion | value | holds |\n|---|---|---|")
    println(md, "| C1 two families in the application (ΔBIC ≥ 6) | $(fmt(ra.dbic; d=1)) | $(c1 ? "yes" : "no") |")
    println(md, "| C2 types repeat (κ ≥ 0.5) | $(fmt(ra.kappa)) | $(c2 ? "yes" : "no") |")
    println(md, "| C3 no two families in the control (ΔBIC < 6) | $(fmt(rc.dbic; d=1)) | $(c3 ? "yes" : "no") |")
    println(md, "| C4 position does not decide (spread of type-1 shares < 15 points) | $(fmt(100(maximum(ra.shares) - minimum(ra.shares)); d=0)) points | $(c4 ? "yes" : "no") |\n")
    if c1 && c2 && c3 && c4
        # sequences of the molecules with stable counts
        post = Dict(ra.keys[i] => ra.post[i] for i in eachindex(ra.keys))
        sc = Dict(k => length(u.odd) + length(u.even) for (k, u) in ua)
        rows = String[]; reads = Vector{Char}[]
        for (t, fs) in sort(collect(membersa); by=first)
            length(fs) >= 4 || continue
            counts = [nscana[f] for f in fs]; maximum(counts) - minimum(counts) <= 1 || continue
            phys = sort(unique(k[2] for k in keys(ua) if k[1] == t); by=physical_index)
            s = [haskey(post, (t, p)) && sc[(t, p)] >= 3 ? (post[(t, p)] >= 0.8 ? '1' : post[(t, p)] <= 0.2 ? '0' : '?') : '?' for p in phys]
            push!(reads, s); push!(rows, "$t\t$(length(fs))\t$(join(s))")
        end
        open(joinpath(out, "sequences.tsv"), "w") do io
            println(io, "track\tscans\ttypes (1 = brighter family, ? = not read)")
            foreach(r -> println(io, r), rows)
        end
        typed = [c for s in reads for c in s if c != '?']
        runs = Int[]; unlike = 0; pairs = 0
        for s in reads
            cur = 0; prev = '?'
            for c in s
                if c == '?'
                    cur > 0 && push!(runs, cur); cur = 0; prev = '?'; continue
                end
                if prev != '?'
                    pairs += 1; unlike += c != prev
                    c == prev ? (cur += 1) : (push!(runs, cur); cur = 1)
                else
                    cur = 1
                end
                prev = c
            end
            cur > 0 && push!(runs, cur)
        end
        h = bytes2hex(sha256(read(joinpath(out, "sequences.tsv"))))
        println(md, "## Sample-level predictions (registered before any comparison)\n")
        println(md, "$(length(reads)) molecules with stable counts; $(length(typed)) units read of $(sum(length, reads)). Type-1 share $(fmt(100count(==('1'), typed) / max(1, length(typed)); d=0)) %; mean run length of one type $(fmt(isempty(runs) ? NaN : mean(runs))); share of unlike neighbours $(fmt(100unlike / max(1, pairs); d=0)) %. `sequences.tsv` sha256 $h.")
    else
        println(md, "A criterion fails: no sequence is written, as declared.")
    end
    write(joinpath(out, "report.md"), take!(md))
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && typing_main()

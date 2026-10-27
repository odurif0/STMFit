#!/usr/bin/env julia
# Two unit types along the direction where the application run varies more,
# reproducibly, than a fixed-sequence control (campaign T2; docs/src/journal.md
# 2026-10-27, for the measurements, the direction, the folds and the criteria,
# fixed before this script was run). Label-free.
#
# Usage:
#   julia --project=. scripts/data/type_units_contrastive.jl --run RUN_DIR --control CONTROL_RUN_DIR \
#       --outdir NEW_DIR [--cap 16] [--seed 20261027]
#
# Units, features and residuals as scripts/data/type_units.jl, with the patch
# descriptors averaged over the forward and backward patches. Molecules of each
# run are split in two folds; the direction learned on one fold of both runs
# scores the other fold, so every score is out of fold. Outputs: units.tsv,
# report.md and, when C1–C4 hold, sequences.tsv with sample-level predictions.
include(joinpath(@__DIR__, "type_units.jl"))
using LinearAlgebra, Random

function contrastive_options(args)
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

"Unit means of a run (all, odd and even scans), distance class and scans per unit, for units with ≥ 2 scans in each half."
function unit_matrices(units)
    keys_ = sort([k for (k, u) in units if length(u.odd) >= 2 && length(u.even) >= 2])
    A = [meanfinite(vcat(units[k].odd, units[k].even)) for k in keys_]
    Od = [meanfinite(units[k].odd) for k in keys_]; Ev = [meanfinite(units[k].even) for k in keys_]
    ok = [all(isfinite, A[i]) && all(isfinite, Od[i]) && all(isfinite, Ev[i]) for i in eachindex(keys_)]
    k = keys_[ok]
    return (keys=k, A=Matrix(reduce(hcat, A[ok])'), Od=Matrix(reduce(hcat, Od[ok])'), Ev=Matrix(reduce(hcat, Ev[ok])'),
            d=[min(modal(units[x].d), 4) for x in k], scans=[length(units[x].odd) + length(units[x].even) for x in k])
end

"Fold (1 or 2) of every track: tracks sorted, shuffled with `seed`, alternated."
function track_folds(tracks, seed)
    t = shuffle(MersenneTwister(seed), sort(unique(tracks)))
    return Dict(x => isodd(i) ? 1 : 2 for (i, x) in enumerate(t))
end

repcov(Zo, Ze) = (Zo' * Ze + Ze' * Zo) ./ (2size(Zo, 1))

function mixtures(s; seed)
    X = reshape(s, :, 1)
    M1 = fit_mixture(X, 1; seed); M2 = fit_mixture(X, 2; seed)
    return M2, mixture_bic(M1, length(s)) - mixture_bic(M2, length(s))
end

function contrastive_main(args=ARGS)
    o = contrastive_options(args)
    o === nothing && return println("type_units_contrastive.jl --run RUN_DIR --control CONTROL_RUN_DIR --outdir NEW_DIR [--cap 16] [--seed 20261027]")
    cap, seed = parse(Int, o["--cap"]), parse(Int, o["--seed"])
    ua, membersa, nscana = run_units(o["--run"]; cap, bwd=true)
    uc, _, _ = run_units(o["--control"]; cap, bwd=true)
    a = unit_matrices(ua); c = unit_matrices(uc)
    σ = std(vcat(a.A .- mean(a.A; dims=1), c.A .- mean(c.A; dims=1)); dims=1); σ[σ .== 0] .= 1
    scale(R) = (μ = mean(R.A; dims=1); ((R.A .- μ) ./ σ, (R.Od .- μ) ./ σ, (R.Ev .- μ) ./ σ))
    Za, Zoa, Zea = scale(a); Zc, Zoc, Zec = scale(c)
    fa = track_folds(first.(a.keys), seed); fc = track_folds(first.(c.keys), seed + 1)
    ia = [fa[k[1]] for k in a.keys]; ic = [fc[k[1]] for k in c.keys]
    sa = zeros(length(a.keys)); soa = similar(sa); sea = similar(sa)
    sc = zeros(length(c.keys)); soc = similar(sc); sec = similar(sc)
    dirs = Vector{Vector{Float64}}(); lams = Float64[]
    for k in 1:2
        tra = ia .== k; trc = ic .== k
        v, λ = contrastive_direction(repcov(Zoa[tra, :], Zea[tra, :]), repcov(Zoc[trc, :], Zec[trc, :]))
        !isempty(dirs) && dot(v, dirs[1]) < 0 && (v = -v)
        push!(dirs, v); push!(lams, λ)
        tea = .!tra; tec = .!trc
        sa[tea] = Za[tea, :] * v; soa[tea] = Zoa[tea, :] * v; sea[tea] = Zea[tea, :] * v
        sc[tec] = Zc[tec, :] * v; soc[tec] = Zoc[tec, :] * v; sec[tec] = Zec[tec, :] * v
    end
    lo, hi = quantile(sa, 0.025), quantile(sa, 0.975)
    wa = clamp.(sa, lo, hi)
    M2, dbic = mixtures(wa; seed)
    P = posterior(M2, reshape(wa, :, 1))
    bright = argmax([sum(P[:, j] .* Za[:, 1]) / sum(P[:, j]) for j in 1:2])
    post = P[:, bright]; type1 = post .> 0.5
    to = posterior(M2, reshape(clamp.(soa, lo, hi), :, 1))[:, bright] .> 0.5
    te = posterior(M2, reshape(clamp.(sea, lo, hi), :, 1))[:, bright] .> 0.5
    κ = cohen_kappa(to, te)
    shares = [count(i -> a.d[i] == j && type1[i], eachindex(a.d)) / max(1, count(==(j), a.d)) for j in 1:4]
    Mc, dbicc = mixtures(winsorize(sc); seed)
    c1 = dbic >= 6 && minimum(M2.weights) >= 0.15
    c2 = κ >= 0.5
    c3 = !(dbicc >= 6 && minimum(Mc.weights) >= 0.15)
    c4 = maximum(shares) - minimum(shares) < 0.15
    out = o["--outdir"]; mkpath(out)
    open(joinpath(out, "units.tsv"), "w") do io
        println(io, "run\ttrack\tphysical_lobe\tfold\td_class\tscans\tscore\ttype\tposterior_type1\ttype_odd\ttype_even")
        for i in eachindex(a.keys)
            println(io, join(["application", a.keys[i]..., ia[i], a.d[i], a.scans[i], fmt(sa[i]; d=4), type1[i] ? 1 : 0, fmt(post[i]; d=4), to[i] ? 1 : 0, te[i] ? 1 : 0], '\t'))
        end
        for i in eachindex(c.keys)
            println(io, join(["control", c.keys[i]..., ic[i], c.d[i], c.scans[i], fmt(sc[i]; d=4), "", "", "", ""], '\t'))
        end
    end
    md = IOBuffer()
    println(md, "# Contrastive unit types (campaign T2)\n")
    println(md, "Application `$(o["--run"])`, control `$(o["--control"])`; $(length(a.keys)) and $(length(c.keys)) units (≥ 2 scans in each half). Leading λ by fold: $(join(fmt.(lams), ", ")). Direction loadings (fold 1): " *
                join(["$(FEATURES[j]) $(fmt(dirs[1][j]))" for j in sortperm(abs.(dirs[1]); rev=true)[1:5]], ", ") * ".\n")
    println(md, "| criterion | value | holds |\n|---|---|---|")
    println(md, "| C1 application: ΔBIC ≥ 6 and weights ≥ 0.15 | ΔBIC $(fmt(dbic; d=1)), weights $(join(fmt.(sort(M2.weights)), " / ")) | $(c1 ? "yes" : "no") |")
    println(md, "| C2 types repeat (κ ≥ 0.5) | $(fmt(κ)) | $(c2 ? "yes" : "no") |")
    println(md, "| C3 control without two families | ΔBIC $(fmt(dbicc; d=1)), weights $(join(fmt.(sort(Mc.weights)), " / ")) | $(c3 ? "yes" : "no") |")
    println(md, "| C4 type-1 shares by distance class within 15 points | $(join([fmt(100s; d=0) for s in shares], " / ")) % | $(c4 ? "yes" : "no") |\n")
    if c1 && c2 && c3 && c4
        postd = Dict(a.keys[i] => post[i] for i in eachindex(a.keys))
        sc_ = Dict(k => length(u.odd) + length(u.even) for (k, u) in ua)
        rows = String[]; reads = Vector{Char}[]
        for (t, fs) in sort(collect(membersa); by=first)
            length(fs) >= 4 || continue
            counts = [nscana[f] for f in fs]; maximum(counts) - minimum(counts) <= 1 || continue
            phys = sort(unique(k[2] for k in keys(ua) if k[1] == t); by=physical_index)
            s = [haskey(postd, (t, p)) && sc_[(t, p)] >= 3 ? (postd[(t, p)] >= 0.8 ? '1' : postd[(t, p)] <= 0.2 ? '0' : '?') : '?' for p in phys]
            push!(reads, s); push!(rows, "$t\t$(length(fs))\t$(join(s))")
        end
        open(joinpath(out, "sequences.tsv"), "w") do io
            println(io, "track\tscans\ttypes (1 = brighter family, ? = not read)")
            foreach(r -> println(io, r), rows)
        end
        typed = [ch for s in reads for ch in s if ch != '?']
        pairs = [(s[i], s[i+1]) for s in reads for i in 1:length(s)-1 if s[i] != '?' && s[i+1] != '?']
        h = bytes2hex(sha256(read(joinpath(out, "sequences.tsv"))))
        println(md, "## Sample-level predictions (registered before any comparison)\n")
        println(md, "$(length(reads)) molecules with stable counts; $(length(typed)) of $(sum(length, reads; init=0)) units read; type-1 share $(fmt(100count(==('1'), typed) / max(1, length(typed)); d=0)) %; unlike neighbours $(fmt(100count(p -> p[1] != p[2], pairs) / max(1, length(pairs)); d=0)) % of $(length(pairs)) read pairs. `sequences.tsv` sha256 $h.")
    else
        println(md, "A criterion fails: no sequence is written, as declared.")
    end
    write(joinpath(out, "report.md"), take!(md))
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && contrastive_main()

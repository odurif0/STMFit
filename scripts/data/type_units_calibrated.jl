#!/usr/bin/env julia
# Two unit types or a continuum along the 10–20mer-specific direction, with the
# type width calibrated on the fixed-sequence control, measurement noise and
# neighbour contamination removed (campaign T3; docs/src/journal.md 2026-10-27,
# for every choice, fixed before this script was run). Label-free.
#
# Usage:
#   julia --project=. scripts/data/type_units_calibrated.jl --run RUN_DIR --control CONTROL_RUN_DIR \
#       --outdir NEW_DIR [--cap 16] [--seed 20261027]
#
# Direction, folds and measurements as scripts/data/type_units_contrastive.jl.
# Outputs: units.tsv, report.md and, when C1–C4 hold, sequences.tsv with
# sample-level predictions.
include(joinpath(@__DIR__, "type_units_contrastive.jl"))

"Per-scan scores of the units of one run along the fold directions, corrected for neighbour contamination; returns per-unit score lists, β and σ_w²."
function corrected_scores(units, R, μ, σ, folds, dirs)
    raw = Dict{Tuple{String,Int},Float64}(); owner = Dict{Tuple{String,Int},Tuple{Tuple{String,String},Bool}}()
    for k in R.keys
        u = units[k]; v = dirs[3 - folds[k[1]]]                       # direction learned on the other fold
        for (obs, w, isodd_) in ((u.odd, u.where_odd, true), (u.even, u.where_even, false)), (x, fk) in zip(obs, w)
            z = (x .- vec(μ)) ./ vec(σ)
            all(isfinite, z) || continue
            raw[fk] = dot(z, v); owner[fk] = (k, isodd_)
        end
    end
    nb(fk) = (vals = [raw[(fk[1], fk[2] + d)] for d in (-1, 1) if haskey(raw, (fk[1], fk[2] + d))]; isempty(vals) ? NaN : mean(vals))
    xs = Float64[]; ys = Float64[]
    for (fk, y) in raw
        x = nb(fk); isfinite(x) && (push!(xs, x); push!(ys, y))
    end
    β = sum((xs .- mean(xs)) .* (ys .- mean(ys))) / sum(abs2, xs .- mean(xs))
    scores = Dict(k => (odd=Float64[], even=Float64[]) for k in R.keys)
    for (fk, y) in raw
        x = nb(fk); k, o = owner[fk]
        push!(o ? scores[k].odd : scores[k].even, isfinite(x) ? y - β * x : y)
    end
    num = 0.0; den = 0
    for sc in values(scores)
        a = vcat(sc.odd, sc.even); length(a) >= 2 || continue
        num += sum(abs2, a .- mean(a)); den += length(a) - 1
    end
    return scores, β, num / den
end

function calibrated_main(args=ARGS)
    o = contrastive_options(args)
    o === nothing && return println("type_units_calibrated.jl --run RUN_DIR --control CONTROL_RUN_DIR --outdir NEW_DIR [--cap 16] [--seed 20261027]")
    cap, seed = parse(Int, o["--cap"]), parse(Int, o["--seed"])
    ua, membersa, nscana = run_units(o["--run"]; cap, bwd=true)
    uc, _, _ = run_units(o["--control"]; cap, bwd=true)
    a = unit_matrices(ua); c = unit_matrices(uc)
    σ = std(vcat(a.A .- mean(a.A; dims=1), c.A .- mean(c.A; dims=1)); dims=1); σ[σ .== 0] .= 1
    μa = mean(a.A; dims=1); μc = mean(c.A; dims=1)
    Zoa, Zea = (a.Od .- μa) ./ σ, (a.Ev .- μa) ./ σ; Zoc, Zec = (c.Od .- μc) ./ σ, (c.Ev .- μc) ./ σ
    fa = track_folds(first.(a.keys), seed); fc = track_folds(first.(c.keys), seed + 1)
    ia = [fa[k[1]] for k in a.keys]; ic = [fc[k[1]] for k in c.keys]
    dirs = Vector{Vector{Float64}}()
    for k in 1:2
        v, _ = contrastive_direction(repcov(Zoa[ia .== k, :], Zea[ia .== k, :]), repcov(Zoc[ic .== k, :], Zec[ic .== k, :]))
        !isempty(dirs) && dot(v, dirs[1]) < 0 && (v = -v)
        push!(dirs, v)
    end
    sca, βa, w2a = corrected_scores(ua, a, μa, σ, fa, dirs)
    scc, βc, w2c = corrected_scores(uc, c, μc, σ, fc, dirs)
    stats(sc, keys_) = (m = [mean(vcat(sc[k].odd, sc[k].even)) for k in keys_], n = [length(sc[k].odd) + length(sc[k].even) for k in keys_],
                        mo = [mean(sc[k].odd) for k in keys_], no = [length(sc[k].odd) for k in keys_],
                        me = [mean(sc[k].even) for k in keys_], ne = [length(sc[k].even) for k in keys_])
    A = stats(sca, a.keys); C = stats(scc, c.keys)
    τ2 = max(var(C.m) - mean(w2c ./ C.n), 1e-6)
    function test(S, w2, τ2)
        T = fit_two_types(S.m, τ2 .+ w2 ./ S.n; seed); Cn = fit_continuum(S.m, w2 ./ S.n); n = length(S.m)
        return T, Cn, (-2Cn.loglik + 2log(n)) - (-2T.loglik + 3log(n))
    end
    Ta, Ca, dbica = test(A, w2a, τ2)
    Tc, Cc, dbicc = test(C, w2c, τ2)
    pa = type_posterior(A.m, τ2 .+ w2a ./ A.n, Ta.π, Ta.μ1, Ta.μ2)
    za1 = (a.A[:, 1] .- μa[1]) ./ σ[1]
    one_is_bright = sum(pa .* za1) / sum(pa) >= sum((1 .- pa) .* za1) / sum(1 .- pa)
    post = one_is_bright ? pa : 1 .- pa; type1 = post .> 0.5
    half(m, n) = (p = type_posterior(m, τ2 .+ w2a ./ n, Ta.π, Ta.μ1, Ta.μ2); (one_is_bright ? p : 1 .- p) .> 0.5)
    to = half(A.mo, A.no); te = half(A.me, A.ne); κ = cohen_kappa(to, te)
    shares = [count(i -> a.d[i] == j && type1[i], eachindex(a.d)) / max(1, count(==(j), a.d)) for j in 1:4]
    small(T) = min(T.π, 1 - T.π)
    c1 = dbica >= 6 && small(Ta) >= 0.10
    c2 = κ >= 0.5
    c3 = !(dbicc >= 6 && small(Tc) >= 0.10)
    c4 = maximum(shares) - minimum(shares) < 0.15
    out = o["--outdir"]; mkpath(out)
    open(joinpath(out, "units.tsv"), "w") do io
        println(io, "track\tphysical_lobe\td_class\tscans\tscore\tposterior_type1\ttype\ttype_odd\ttype_even")
        for i in eachindex(a.keys)
            println(io, join([a.keys[i]..., a.d[i], A.n[i], fmt(A.m[i]; d=4), fmt(post[i]; d=4), type1[i] ? 1 : 0, to[i] ? 1 : 0, te[i] ? 1 : 0], '\t'))
        end
    end
    md = IOBuffer()
    println(md, "# Calibrated types or continuum (campaign T3)\n")
    println(md, "Application $(length(a.keys)) units, control $(length(c.keys)) units. Neighbour slope β: application $(fmt(βa; d=3)), control $(fmt(βc; d=3)). Per-scan noise σ_w²: $(fmt(w2a; d=3)) and $(fmt(w2c; d=3)). Type width τ² from the control: $(fmt(τ2; d=3)). Application continuum width τa²: $(fmt(Ca.τ2; d=3)); control: $(fmt(Cc.τ2; d=3)).\n")
    println(md, "| criterion | value | holds |\n|---|---|---|")
    println(md, "| C1 application: ΔBIC ≥ 6, smaller type ≥ 10% | ΔBIC $(fmt(dbica; d=1)), types $(fmt(100small(Ta); d=0)) / $(fmt(100(1 - small(Ta)); d=0)) %, means $(fmt(Ta.μ1)) and $(fmt(Ta.μ2)) | $(c1 ? "yes" : "no") |")
    println(md, "| C2 types repeat (κ ≥ 0.5) | $(fmt(κ)) | $(c2 ? "yes" : "no") |")
    println(md, "| C3 control does not meet C1 | ΔBIC $(fmt(dbicc; d=1)), types $(fmt(100small(Tc); d=0)) % | $(c3 ? "yes" : "no") |")
    println(md, "| C4 type-1 shares by distance class within 15 points | $(join([fmt(100s; d=0) for s in shares], " / ")) % | $(c4 ? "yes" : "no") |\n")
    if c1 && c2 && c3 && c4
        postd = Dict(a.keys[i] => post[i] for i in eachindex(a.keys)); nd = Dict(a.keys[i] => A.n[i] for i in eachindex(a.keys))
        rows = String[]; reads = Vector{Char}[]
        for (t, fs) in sort(collect(membersa); by=first)
            length(fs) >= 4 || continue
            counts = [nscana[f] for f in fs]; maximum(counts) - minimum(counts) <= 1 || continue
            phys = sort(unique(k[2] for k in keys(ua) if k[1] == t); by=physical_index)
            s = [haskey(postd, (t, p)) && nd[(t, p)] >= 3 ? (postd[(t, p)] >= 0.8 ? '1' : postd[(t, p)] <= 0.2 ? '0' : '?') : '?' for p in phys]
            push!(reads, s); push!(rows, "$t\t$(length(fs))\t$(join(s))")
        end
        open(joinpath(out, "sequences.tsv"), "w") do io
            println(io, "track\tscans\ttypes (1 = brighter type, ? = not read)")
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

abspath(PROGRAM_FILE) == abspath(@__FILE__) && calibrated_main()

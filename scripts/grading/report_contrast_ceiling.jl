#!/usr/bin/env julia
# Reading ceiling of a per-unit K against N call from a lobe signal (grading
# side: it reads the control sequence; the method never does).
#
# Usage:
#   julia --project=. scripts/grading/report_contrast_ceiling.jl --run NAME=RUN_DIR \
#       [--run NAME=RUN_DIR ...] --outdir NEW_DIR [--control 010010] [--draws 1000] [--seed 20261025]
#
# On a run of a fixed control molecule (default NKNNKN written 010010, 1 = K,
# palindromic so the orientation of a scan does not matter) it measures, per
# signal, the contrast between the K lobes and the core N lobes (positions 3 and
# 4), averaged over scans of each molecule and then over molecules. In the
# control the sequence is a function of the position, so this contrast contains
# the position effect (K lobes are the second from either end): it is an upper
# bound for the chemical contrast. It is set against the unit-level variation
# that repeats across scans of one molecule (between-lobe sd, after removing the
# distance-to-end profile, from STMFit.UnitReliability) and against the
# within-lobe sd. d' = |contrast| / sd and the accuracy of the best single
# threshold between two equally likely classes is Φ(d'/2). See docs/src/journal.md
# (campaign R, 2026-10-25) for the definitions and the decision rule.
using STMFit.UnitReliability
using Printf, Random, Statistics

function ceiling_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,Any}("--run" => Pair{String,String}[], "--control" => "010010", "--draws" => "1000", "--seed" => "20261025"); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--run", "--outdir", "--control", "--draws", "--seed") && i < length(args) || error("Unknown option or missing value: $k")
        if k == "--run"
            name, path = split(args[i+1], "="; limit=2)
            push!(o[k], String(name) => String(path))
        else
            o[k] = args[i+1]
        end
        i += 2
    end
    isempty(o["--run"]) && error("Required: --run NAME=RUN_DIR")
    haskey(o, "--outdir") || error("Required: --outdir")
    o["--control"] == reverse(o["--control"]) || error("The control must read the same from both ends (orientation is unknown)")
    all(c -> c in "01", o["--control"]) || error("--control is a 0/1 string, 1 = K")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Mean scan-centred signal of the K, core and terminal lobes of every molecule (scans with the control's length)."
function class_means(ru, kind, control)
    n = length(control)
    cls(k) = control[k] == '1' ? :K : (k == 1 || k == n ? :terminal : :core)
    out = Dict{String,NamedTuple{(:K, :core, :terminal),NTuple{3,Float64}}}()
    for (t, fs) in ru.tracks
        acc = Dict(:K => Float64[], :core => Float64[], :terminal => Float64[])
        for f in fs
            length(ru.lobes[f]) == n || continue
            v = centred(ru.lobes[f], kind)
            any(isnan, v) && continue
            for k in 1:n
                push!(acc[cls(k)], v[k])
            end
        end
        all(!isempty(acc[c]) for c in (:K, :core, :terminal)) &&
            (out[t] = (K=mean(acc[:K]), core=mean(acc[:core]), terminal=mean(acc[:terminal])))
    end
    return out
end

fmt(x; d=3) = isnan(x) ? "NA" : @sprintf("%.*f", d, x)
ci(x) = isempty(x) ? (NaN, NaN) : (quantile(x, 0.025), quantile(x, 0.975))

function ceiling_main(args=ARGS)
    o = ceiling_options(args)
    o === nothing && return println("report_contrast_ceiling.jl --run NAME=RUN_DIR [--run ...] --outdir NEW_DIR [--control 010010]")
    control = o["--control"]; draws = parse(Int, o["--draws"]); rng = MersenneTwister(parse(Int, o["--seed"]))
    md = IOBuffer()
    println(md, "# Reading ceiling of K against N from one lobe signal\n")
    println(md, "Control $control (1 = K). Contrast = K lobes minus core N lobes, on the scan-centred log signal, per molecule then over molecules; it contains the position effect, so it is an upper bound for chemistry. Between-lobe sd = unit-level variation that repeats across scans of a molecule (distance-to-end profile removed, d ≥ 2); within-lobe sd = what does not repeat. Intervals: bootstrap over molecules ($draws draws).\n")
    println(md, "| run | signal | molecules | K − core N | 95% interval | sd between | sd within | d′ (mean of scans) | ceiling | d′ (one scan) | ceiling | contrast needed for 80% |\n|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for (name, path) in o["--run"]
        ru = load_run_units(path; name)
        for kind in SIGNALS
            cm = class_means(ru, kind, control)
            isempty(cm) && continue
            c = [v.K - v.core for v in values(cm)]
            cb = [mean(c[rand(rng, 1:length(c), length(c))]) for _ in 1:draws]
            tm = [moments(g) for g in values(track_groups(ru, kind; profile=distance_profile(ru, kind)))]
            m = sum(tm; init=Moments()); b, w, _ = variance_components(m)
            σb = sqrt(max(b, 0.0)); σw = sqrt(w)
            Bb, _, _ = bootstrap_components(tm; draws, rng)
            dp = abs(mean(c)) / σb; ds = abs(mean(c)) / sqrt(σb^2 + σw^2)
            dpb = [abs(cb[i]) / sqrt(max(Bb[i], 1e-12)) for i in 1:min(length(cb), length(Bb))]
            dsb = [abs(cb[i]) / sqrt(max(Bb[i], 1e-12) + w) for i in 1:min(length(cb), length(Bb))]
            lo, hi = ci(cb)
            needed = 1.6821 * σb                       # Φ(d′/2) = 0.80  ->  d′ = 1.6821
            println(md, "| $name | $kind | $(length(c)) | $(fmt(mean(c))) | $(fmt(lo)) to $(fmt(hi)) | $(fmt(σb)) | $(fmt(σw)) | ",
                    "$(fmt(dp, d=2)) | $(fmt(100 * normal_cdf(dp / 2), d=0)) % ($(fmt(100 * normal_cdf(ci(dpb)[1] / 2), d=0)) to $(fmt(100 * normal_cdf(ci(dpb)[2] / 2), d=0))) | ",
                    "$(fmt(ds, d=2)) | $(fmt(100 * normal_cdf(ds / 2), d=0)) % | $(fmt(needed)) |")
        end
    end
    println(md, "\nThe ceiling is the accuracy of the best threshold on that signal alone, for two equally likely classes: it cannot be beaten by a better classifier of the same signal, only by a signal with a larger contrast or a smaller between-lobe sd.")
    mkpath(o["--outdir"])
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && ceiling_main()

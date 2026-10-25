#!/usr/bin/env julia
# Dependence of the frozen, label-free unit assignment of an application run
# (10-20mers) on the subset of the 6-mer scans it is learned on. Grading side:
# it reads the 6-mer control sequence and the tables that define the subsets;
# the models themselves read no label.
#
# Usage:
#   julia -t 4 --project=. scripts/grading/report_training_stability.jl \
#       --run6 RUN_6MER --run-app RUN_APPLICATION --config config/unit_assignment_frozen_labelfree.toml \
#       --manifest benchmarks/chitosan_6mer_benchmark_v2.toml --quality benchmarks/chitosan_6mer_cu100_quality.tsv \
#       --outdir NEW_DIR [--random-size 167] [--random-seeds 5]
#
# For every training subset (see `subsets`) a fully frozen model is learned on
# the scans of RUN_6MER in the subset (reference rule of the config) and applied
# unchanged to RUN_APPLICATION. Reported: the call of every application lobe under every model
# (application_calls.tsv, with the share of models calling it bright), agreement of the application calls
# between models, class-1 fraction by distance to the nearest end, and the
# accuracy by position class on the 6-mer scans outside the subset. Declared in
# docs/src/journal.md (campaign D, 2026-10-25); nothing is selected on it.
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment
using LinearAlgebra, Printf, Random, Statistics, TOML

const CONTROL = "010010"       # NKNNKN, 1 = K

function stability_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--random-size" => "167", "--random-seeds" => "5"); i = 1
    while i <= length(args)
        args[i] in ("--run6", "--run-app", "--config", "--manifest", "--quality", "--outdir", "--random-size", "--random-seeds") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run6", "--run-app", "--config", "--manifest", "--quality", "--outdir")) ||
        error("Required: --run6 --run-app --config --manifest --quality --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

function tsv_rows(path)
    lines = [l for l in readlines(path) if !isempty(strip(l)) && !startswith(l, "#")]
    header = String.(split(lines[1], '\t'; keepempty=true))
    return [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

"Scan facts of the run (manifest key = run file name): session, cohort, quality-table flags."
function scan_facts(manifest_path, quality_path, files)
    m = TOML.parsefile(manifest_path)["files"]
    q = Dict(r["raw"] => r for r in tsv_rows(quality_path))
    facts = Dict{String,NamedTuple}()
    for f in files
        haskey(m, f) || continue
        r = get(q, m[f]["raw_path"], nothing)
        r === nothing && continue
        facts[f] = (session=m[f]["session"], cohort=m[f]["cohort"], in_full=r["in_full"] == "yes",
                    rule=r["rule_accepted"] == "yes", elite=r["elite"] == "yes")
    end
    return facts
end

"Calendar year of a session name (20240307_LHe_Cu100 or Cu100LHe01052025)."
function session_year(s)
    m = match(r"^(20\d\d)\d{4}_", s); m !== nothing && return parse(Int, m[1])
    m = match(r"(20\d\d)$", s); m !== nothing && return parse(Int, m[1])
    return 0
end

function subsets(facts, consensus, size, seeds)
    full = sort([f for (f, x) in facts if x.in_full])
    track = Dict(basename(r["filepath"]) => r["track"] for r in consensus)
    out = Pair{String,Tuple{Vector{String},Int}}[]       # name => (files, cap per track)
    push!(out, "full" => (full, 5))
    for cap in (1, 3, 10, 100000)
        push!(out, (cap == 100000 ? "full_nocap" : "full_cap$cap") => (full, cap))
    end
    push!(out, "elite" => (sort([f for (f, x) in facts if x.elite]), 5))
    push!(out, "rule419" => (sort([f for (f, x) in facts if x.rule]), 5))
    push!(out, "original145" => (sort([f for (f, x) in facts if x.cohort == "full145"]), 5))
    push!(out, "sessions2024" => (sort([f for f in full if session_year(facts[f].session) == 2024]), 5))
    push!(out, "sessions2025" => (sort([f for f in full if session_year(facts[f].session) == 2025]), 5))
    for s in 1:seeds
        rng = MersenneTwister(1000 + s)
        push!(out, "rand_scans_s$s" => (sort(shuffle(rng, full)[1:min(size, length(full))]), 5))
    end
    tracks = sort(unique(track[f] for f in full))
    for s in 1:seeds
        rng = MersenneTwister(2000 + s); picked = String[]
        for t in shuffle(rng, tracks)
            append!(picked, [f for f in full if track[f] == t])
            length(picked) >= size && break
        end
        push!(out, "rand_molecules_s$s" => (sort(picked), 5))
    end
    return out
end

function train_on(R, cfg, files::Set{String}, cap)
    cfg2 = deepcopy(cfg); cfg2["selection"]["reference_max_scans_per_track"] = cap
    st = load_frozen_settings(cfg2)
    ref = reference_files([r for r in R.consensus if basename(r["filepath"]) in files], st.max_per_track)
    ref_rows = findall(in(ref), R.data.file); theta_rows = findall(in(files), R.data.file)
    model = train_reference_model(R.data, R.fisher, ref_rows, theta_rows, R.physical, cfg2)
    return model, length(ref)
end

"Distance to the nearest chain end of every lobe of a run table."
function distance_to_end(data)
    n = Dict{String,Int}(); for f in data.file; n[f] = get(n, f, 0) + 1; end
    return [min(data.lobe[i], n[data.file[i]] + 1 - data.lobe[i]) for i in eachindex(data.file)]
end

"Accuracy by position class on the six-lobe scans of `rows` (calls aligned with `rows`)."
function position_accuracy(data, rows, calls)
    byfile = Dict{String,Vector{Tuple{Int,Int8}}}()
    for (r, c) in zip(rows, calls)
        push!(get!(byfile, data.file[r], Tuple{Int,Int8}[]), (data.lobe[r], c))
    end
    ok = Dict(:terminal => [0, 0], :K => [0, 0], :core => [0, 0])
    for (f, v) in byfile
        length(v) == 6 || continue
        for (k, c) in v
            cls = CONTROL[k] == '1' ? :K : (k == 1 || k == 6 ? :terminal : :core)
            ok[cls][2] += 1; ok[cls][1] += (c == parse(Int, CONTROL[k]))
        end
    end
    tot = [sum(ok[c][1] for c in keys(ok)), sum(ok[c][2] for c in keys(ok))]
    return (terminal=ok[:terminal], K=ok[:K], core=ok[:core], all=tot)
end

pct(a) = a[2] == 0 ? NaN : 100 * a[1] / a[2]
f1(x) = isnan(x) ? "NA" : @sprintf("%.1f", x)

function stability_main(args=ARGS)
    o = stability_options(args)
    o === nothing && return println("report_training_stability.jl --run6 RUN --run-app RUN --config TOML --manifest TOML --quality TSV --outdir NEW_DIR")
    BLAS.set_num_threads(1)
    cfg = load_config(o["--config"])
    R6 = load_run_assignment(o["--run6"], cfg); RA = load_run_assignment(o["--run-app"], cfg)
    facts = scan_facts(o["--manifest"], o["--quality"], unique(R6.data.file))
    subs = subsets(facts, R6.consensus, parse(Int, o["--random-size"]), parse(Int, o["--random-seeds"]))
    d = distance_to_end(RA.data); dclass = min.(d, 4)
    rows_app = collect(eachindex(RA.data.file))
    calls = Dict{String,Vector{Int8}}(); info = Dict{String,Dict{String,String}}()
    for (name, (files, cap)) in subs
        fset = Set(files)
        model, nref = train_on(R6, cfg, fset, cap)
        fused, _ = apply_reference_model(model, RA.data, RA.fisher, rows_app, RA.physical, cfg)
        calls[name] = Int8.(fused)
        held = findall(f -> !(f in fset), R6.data.file)
        fh, _ = apply_reference_model(model, R6.data, R6.fisher, held, R6.physical, cfg)
        acc = position_accuracy(R6.data, held, Int8.(fh))
        bright = [begin sel = [i for i in rows_app if dclass[i] == c && fused[i] >= 0]; isempty(sel) ? NaN : 100 * count(i -> fused[i] == 1, sel) / length(sel) end for c in 1:4]
        info[name] = Dict("model" => name, "scans" => string(length(files)), "reference_scans" => string(nref), "theta0" => @sprintf("%.3f", model.theta0),
                          "theta1" => @sprintf("%.3f", model.theta1), "pi" => @sprintf("%.3f", model.pi),
                          "heldout_scans" => string(length(unique(R6.data.file[held]))), "acc_all" => f1(pct(acc.all)), "acc_terminal" => f1(pct(acc.terminal)),
                          "acc_K" => f1(pct(acc.K)), "acc_core" => f1(pct(acc.core)),
                          "app_bright_d1" => f1(bright[1]), "app_bright_d2" => f1(bright[2]), "app_bright_d3" => f1(bright[3]), "app_bright_d4plus" => f1(bright[4]),
                          "app_bright_all" => f1(100 * count(==(1), calls[name]) / max(1, count(>=(0), calls[name]))))
        @printf("%-22s scans %3d reference %3d  theta0 %.3f  held-out acc %s  app bright %s%%\n", name, length(files), nref, model.theta0,
                info[name]["acc_all"], info[name]["app_bright_all"]); flush(stdout)
    end
    names = [first(p) for p in subs]
    agree(a, b) = (m = (calls[a] .>= 0) .& (calls[b] .>= 0); count(calls[a][m] .== calls[b][m]) / max(1, count(m)))
    A = [agree(a, b) for a in names, b in names]
    mkpath(o["--outdir"])
    cols = ["model", "scans", "reference_scans", "theta0", "theta1", "pi", "heldout_scans", "acc_all", "acc_terminal", "acc_K", "acc_core",
            "app_bright_d1", "app_bright_d2", "app_bright_d3", "app_bright_d4plus", "app_bright_all"]
    write_table(joinpath(o["--outdir"], "models.tsv"), cols, [info[n] for n in names])
    write_table(joinpath(o["--outdir"], "agreement.tsv"), vcat(["model"], names),
                [merge(Dict("model" => a), Dict(b => @sprintf("%.4f", A[i, j]) for (j, b) in enumerate(names))) for (i, a) in enumerate(names)])
    # per-lobe calls of every model, and how many models call the lobe bright
    out = Dict{String,String}[]
    for (j, i) in enumerate(rows_app)
        votes = [calls[n][j] for n in names if calls[n][j] >= 0]
        row = Dict("file" => RA.data.file[i], "lobe" => string(RA.data.lobe[i]), "physical_lobe" => RA.physical[i], "d" => string(d[i]),
                   "models_classified" => string(length(votes)),
                   "bright_share" => isempty(votes) ? "NA" : @sprintf("%.3f", count(==(1), votes) / length(votes)))
        for n in names
            row[n] = calls[n][j] < 0 ? "?" : string(calls[n][j])
        end
        push!(out, row)
    end
    write_table(joinpath(o["--outdir"], "application_calls.tsv"),
                vcat(["file", "lobe", "physical_lobe", "d", "models_classified", "bright_share"], names), out)
    rs = [n for n in names if startswith(n, "rand_scans")]; rm = [n for n in names if startswith(n, "rand_molecules")]
    pairs(group) = [agree(a, b) for (i, a) in enumerate(group) for b in group[i+1:end]]
    vs_full(group) = [agree("full", n) for n in group]
    floor_ = vs_full(vcat(rs, rm)); q05 = quantile(floor_, 0.05)
    open(joinpath(o["--outdir"], "report.md"), "w") do io
        println(io, "# Training subset and the frozen assignment of the application run\n")
        println(io, "6-mer run `$(o["--run6"])`, application run `$(o["--run-app"])`, config `$(basename(o["--config"]))`. Agreement = share of application lobes with the same fused call (lobes classified by both models).\n")
        println(io, "| model | scans | reference | θ0 | θ1 | π | 6-mer held-out accuracy % (all / terminal / K / core N) | application: % bright at d = 1 / 2 / 3 / ≥4 | all |\n|---|---:|---:|---:|---:|---:|---|---|---:|")
        for n in names
            x = info[n]
            println(io, "| $n | $(x["scans"]) | $(x["reference_scans"]) | $(x["theta0"]) | $(x["theta1"]) | $(x["pi"]) | $(x["acc_all"]) / $(x["acc_terminal"]) / $(x["acc_K"]) / $(x["acc_core"]) | $(x["app_bright_d1"]) / $(x["app_bright_d2"]) / $(x["app_bright_d3"]) / $(x["app_bright_d4plus"]) | $(x["app_bright_all"]) |")
        end
        println(io, "\n## Agreement of the application calls\n")
        println(io, "- Noise floor: agreement with the full-set model of the $(length(floor_)) random subsets: mean $(f1(100 * mean(floor_))) %, 5th percentile $(f1(100 * q05)) %, range $(f1(100 * minimum(floor_))) to $(f1(100 * maximum(floor_))) %.")
        println(io, "- Among random scan subsets: mean pairwise $(f1(100 * mean(pairs(rs)))) % (range $(f1(100 * minimum(pairs(rs)))) to $(f1(100 * maximum(pairs(rs)))) %); among random molecule subsets: $(f1(100 * mean(pairs(rm)))) % (range $(f1(100 * minimum(pairs(rm)))) to $(f1(100 * maximum(pairs(rm)))) %).\n")
        println(io, "| structured subset | agreement with the full-set model % | below the noise floor 5th percentile |\n|---|---:|---|")
        for n in names
            (startswith(n, "rand_") || n == "full") && continue
            a = agree("full", n)
            println(io, "| $n | $(f1(100 * a)) | $(a < q05 ? "yes" : "no") |")
        end
        println(io, "\nFull matrix: `agreement.tsv`.")
    end
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && stability_main()

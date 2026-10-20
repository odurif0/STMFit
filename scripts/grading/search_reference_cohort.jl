#!/usr/bin/env julia
# Campaign B: search of the reference cohort of a frozen unit assignment by its
# benchmark score. Grading side: reads the benchmark labels. Protocol and
# parameters: benchmarks/reference_cohort_search_v1.toml.
#
# Usage:
#   julia -t N --project=. scripts/grading/search_reference_cohort.jl \
#       --protocol benchmarks/reference_cohort_search_v1.toml --mode MODE --outdir NEW_DIR [options]
# Modes (one search per process, so a cluster job can run them side by side):
#   nested --outer SESSION  search on the other development sessions, then grade
#                           the selected cohort and the full pool on SESSION
#                           (trace.tsv, outer.tsv, outer_files.tsv)
#   final                   search on all development sessions (trace.tsv,
#                           cohort.tsv, frozen models model_selected.tsv and
#                           model_pool.tsv with their .units lists)
#   null --seed K           the final search with per-track random truths
#                           (protocol [null]; trace.tsv)
#   report --campaign DIR   merge the nested_*, final and null_* outputs of DIR
#                           (outer.tsv, outer_files.tsv, null_summary.tsv, report.md)
#   sealed --models DIR     apply the frozen models of a final output to the
#                           sealed sessions and grade them once (sealed.tsv,
#                           sealed_files.tsv, predictions_<model>.tsv)
# Multistart protocols ([search] strategy = "multistart_flip_swap"): nested,
# final and null take --start K (0 = the pool, K >= 1 a seeded random subset)
# and write only trace.tsv, start.units and end.units; report picks the best
# start of every target, grades the held-out sessions and writes the final
# models (report.md, outer.tsv, outer_files.tsv, null_summary.tsv, cohort.tsv).
#
# The model, its training and its application are label-free
# (src/frozen_assignment.jl); labels enter only the objective here. Trainings
# are cached by (held-out session, training sessions, units). Results do not
# depend on the thread count.
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment, STMFit.MoleculeFusion
using LinearAlgebra, Printf, Random, TOML

const ROOT = dirname(dirname(@__DIR__))

function search_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--protocol", "--mode", "--outdir", "--outer", "--seed", "--start", "--campaign", "--models") && i < length(args) ||
            error("Unknown option or missing value: $k")
        haskey(o, k) && error("Repeated option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--protocol", "--mode", "--outdir")) || error("Required: --protocol --mode --outdir")
    need = Dict("nested" => "--outer", "final" => "", "null" => "--seed", "report" => "--campaign", "sealed" => "--models")
    haskey(need, o["--mode"]) || error("--mode must be nested, final, null, report or sealed")
    isempty(need[o["--mode"]]) || haskey(o, need[o["--mode"]]) || error("--mode $(o["--mode"]) needs $(need[o["--mode"]])")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Run tables and benchmark split, indexed by the lobes of the assignment data."
struct Setup
    data::AssignmentData
    cfg::Dict{String,Any}
    fusion::FusionSettings
    pi::Float64
    vote::Float64
    session::Vector{String}
    physical::Vector{String}
    file_rows::Dict{String,Vector{Int}}
    file_track::Dict{String,Int}
    session_rows::Dict{String,Vector{Int}}
    session_files::Dict{String,Vector{String}}
    transductive::Vector{Int8}           # final calls of the run itself (promoted method, run cohort)
    units::Dict{Int,Vector{Int}}          # track => rows of its corroborated scans
    unit_session::Dict{Int,String}
    development::Vector{String}
    sealed::Vector{String}
end

function load_setup(protocol)
    d = protocol["data"]; path(p) = isabspath(p) ? p : joinpath(ROOT, p)
    run = path(d["run"])
    cfg = load_config(path(d["config"]))
    get(cfg["selection"], "assignment_model", "") == "frozen_reference" || error("Config is not a frozen_reference assignment")
    get(cfg["selection"], "fusion_emission_source", "") == "reference_sessions_all_scans" || error("Unsupported fusion_emission_source")
    data = load_assignment_data(joinpath(run, "assignment", "features_predictor.tsv"),
                                joinpath(run, "assignment", "patches_bwd17.tsv");
                                scale_fallback=cfg["preprocessing"]["gmm_scale_fallback"])
    manifest = TOML.parsefile(path(d["manifest"]))["files"]
    _, final = read_table(joinpath(run, "consensus", "consensus_summary_final.tsv"))
    rule = Dict(basename(r["filepath"]) => r["count_rule"] for r in final)
    track = Dict(basename(r["filepath"]) => parse(Int, r["track"]) for r in final)
    _, fused = lobe_table(joinpath(run, "fusion", "predictions_fused.tsv"); required=["physical_lobe"])
    files = unique(data.file)
    Set(files) == Set(keys(manifest)) == Set(keys(rule)) || error("Run cohort and manifest differ")
    session = [manifest[f]["session"] for f in data.file]
    physical = [fused[(data.file[i], data.lobe[i])]["physical_lobe"] for i in eachindex(data.file)]
    call(x) = x == "?" ? Int8(-1) : parse(Int8, x)
    transductive = [call(fused[(data.file[i], data.lobe[i])]["predicted"]) for i in eachindex(data.file)]
    file_rows = Dict(f => Int[] for f in files); session_rows = Dict{String,Vector{Int}}()
    for i in eachindex(data.file)
        push!(file_rows[data.file[i]], i); push!(get!(session_rows, session[i], Int[]), i)
    end
    session_files = Dict(s => sort(unique(data.file[r])) for (s, r) in session_rows)
    units = Dict{Int,Vector{Int}}(); unit_session = Dict{Int,String}()
    for f in files
        s = manifest[f]["session"]; t = track[f]
        get!(unit_session, t, s) == s || error("Track $t crosses sessions")
        rule[f] == "agrees" && append!(get!(units, t, Int[]), file_rows[f])
    end
    foreach(sort!, values(units))
    cohort(names) = sort(unique(manifest[f]["session"] for f in files if manifest[f]["cohort"] in names))
    dev = cohort(protocol["split"]["development_cohorts"]); sealed = cohort(protocol["split"]["sealed_cohorts"])
    isempty(intersect(dev, sealed)) || error("Development and sealed sessions overlap")
    fusion = load_fusion_settings(cfg)
    return Setup(data, cfg, fusion, Float64(cfg["selection"]["fusion_prior_pi"]), Float64(cfg["selection"]["vote_threshold"]),
                 session, physical, file_rows, track, session_rows, session_files, transductive, units, unit_session, dev, sealed)
end

"Units (tracks with corroborated scans) of the given sessions."
pool_units(S::Setup, sessions) = sort([u for (u, s) in S.unit_session if s in sessions && haskey(S.units, u)])

"""Train on `units`, fit emission rates on every scan of `train_sessions`.
Returns (model with theta and pi, emission fit ok)."""
function frozen_model(S::Setup, units, train_sessions)
    all(u -> S.unit_session[u] in train_sessions, units) || error("Training unit outside the training sessions")
    rows = isempty(units) ? Int[] : sort(reduce(vcat, [S.units[u] for u in units]))
    model = train_frozen(S.data, rows, S.cfg)
    ref = sort(reduce(vcat, [S.session_rows[s] for s in train_sessions]))
    pr = predict_frozen(model, S.data, ref; vote_threshold=S.vote)
    fit = try
        fit_emission_rates(pr.call, S.physical[ref], S.fusion)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    theta0, theta1 = fit === nothing ? (NaN, NaN) : (fit.theta0, fit.theta1)
    return FrozenAssignment.with_fusion(model, theta0, theta1, S.pi), fit !== nothing
end

"Final calls of `rows` (fused within tracks when emission rates exist)."
function final_calls(S::Setup, model::FrozenModel, rows)
    p = predict_frozen(model, S.data, rows; vote_threshold=S.vote)
    isfinite(model.theta0) || return p.call, p
    return fuse_calls(p.call, S.physical[rows], model.theta0, model.theta1, model.pi, S.fusion.min_scans), p
end

# --- cache of held-out calls: key (held-out session, training sessions, units)
const CACHE = Dict{Tuple{String,Vector{String},Vector{Int}},Vector{Int8}}()
const CACHE_LOCK = ReentrantLock()
const TRAININGS = Threads.Atomic{Int}(0)

function heldout_calls(S::Setup, units::Vector{Int}, train_sessions::Vector{String}, heldout::String)
    key = (heldout, train_sessions, units)
    hit = lock(() -> get(CACHE, key, nothing), CACHE_LOCK)
    hit === nothing || return hit
    model, _ = frozen_model(S, units, train_sessions)
    calls, _ = final_calls(S, model, S.session_rows[heldout])
    Threads.atomic_add!(TRAININGS, 1)
    lock(() -> (CACHE[key] = calls), CACHE_LOCK)
    return calls
end

# --- grading (rule of grade_unit_assignment.jl: own N, first min(N, length)
# positions, identity or reverse alignment chosen by accuracy)

function grade_file(calls::AbstractVector{Int8}, truth::Vector{Int8})
    best = (acc=-1.0, correct=0, classified=0, exact=false)
    for a in (calls, reverse(calls))
        n = min(length(a), length(truth)); correct = classified = 0
        for i in 1:n
            a[i] < 0 && continue
            classified += 1; correct += a[i] == truth[i]
        end
        acc = classified == 0 ? 0.0 : correct / classified
        acc > best.acc && (best = (acc=acc, correct=correct, classified=classified,
            exact=length(a) == length(truth) && classified == length(truth) && correct == length(truth)))
    end
    return (n=length(calls), correct=best.correct, classified=best.classified, exact=best.exact)
end

"Per-file grades of one held-out session's calls."
function grade_session(S::Setup, session::String, calls::Vector{Int8}, truth)
    rows = S.session_rows[session]; pos = Dict(r => k for (k, r) in enumerate(rows))
    return [(file=f, track=S.file_track[f], session=session,
             grade_file(calls[[pos[r] for r in S.file_rows[f]]], truth[f])...) for f in S.session_files[session]]
end

"Totals and the per-molecule criterion over per-file grades."
function summarize(grades, truth)
    files = length(grades); correct = sum(g.correct for g in grades; init=0)
    classified = sum(g.classified for g in grades; init=0)
    possible = sum(length(truth[g.file]) for g in grades; init=0)
    tracks = Dict{Int,Vector{Any}}()
    for g in grades; push!(get!(tracks, g.track, Any[]), g); end
    acc = [sum(g.correct for g in v) / sum(length(truth[g.file]) for g in v) for v in values(tracks)]
    maj = count(v -> 2count(g -> g.exact, v) > length(v), values(tracks))
    return (files=files, possible=possible, correct=correct, classified=classified, errors=classified - correct,
            unclassified=sum(min(g.n, length(truth[g.file])) for g in grades; init=0) - classified,
            exact=count(g -> g.exact, grades), tracks=length(tracks),
            molecule_accuracy=isempty(acc) ? NaN : sum(acc) / length(acc), tracks_majority_exact=maj)
end

"Leave-one-session-out grades of every unit set, computed in parallel."
function loso_grades(S::Setup, sets::Vector{Vector{Int}}, sessions::Vector{String}, truth)
    tasks = [(k, s) for k in eachindex(sets), s in sessions]
    calls = Vector{Vector{Int8}}(undef, length(tasks))
    Threads.@threads :dynamic for t in eachindex(tasks)
        k, s = tasks[t]
        train = filter(!=(s), sessions)
        calls[t] = heldout_calls(S, filter(u -> S.unit_session[u] in train, sets[k]), train, s)
    end
    out = [Any[] for _ in sets]
    for (t, (k, s)) in enumerate(tasks)
        append!(out[k], grade_session(S, s, calls[t], truth))
    end
    return out
end

"""
    floating_search(pool, objective; min_gain=1, floating=true)

`objective(sets)` returns one integer per unit set. Backward elimination from
`pool` with optional re-additions; every accepted step gains >= min_gain, so
the search ends. Ties go to the smallest unit id. Returns (set, value, trace).
"""
function floating_search(pool::Vector{Int}, objective; min_gain::Integer=1, floating::Bool=true, log=nothing)
    C = sort(pool); value = only(objective([C]))
    trace = [(step=0, action="start", unit=0, value=value, gain=0, size=length(C))]
    note(t) = (push!(trace, t); log === nothing || log(t))
    log === nothing || log(trace[1])
    while !isempty(C)
        values = objective([filter(!=(u), C) for u in C])
        gain, k = findmax(values .- value)
        gain >= min_gain || break
        u = C[k]; C = filter(!=(u), C); value = values[k]
        note((step=length(trace), action="remove", unit=u, value=value, gain=gain, size=length(C)))
        while floating
            removed = sort(setdiff(pool, C))
            isempty(removed) && break
            values = objective([sort(vcat(C, r)) for r in removed])
            gain, k = findmax(values .- value)
            gain >= min_gain || break
            C = sort(vcat(C, removed[k])); value = values[k]
            note((step=length(trace), action="add", unit=removed[k], value=value, gain=gain, size=length(C)))
        end
    end
    return C, value, trace
end

"""
    local_search(pool, start, objective; min_gain=1, swap_top=5)

Steepest ascent from `start`: every step evaluates flipping each unit of
`pool` (remove it if in the set, add it otherwise) and takes the best flip.
When no flip gains >= min_gain, it evaluates swaps (remove one unit, add
one) among the swap_top best removals and the swap_top best additions of
that flip round. Stops when neither gains. Ties go to the smallest ids.
Returns (set, value, trace); swaps record the added unit in `added`.
"""
function local_search(pool::Vector{Int}, start::Vector{Int}, objective; min_gain::Integer=1, swap_top::Integer=5, log=nothing)
    pool = sort(pool); C = sort(start); value = only(objective([C]))
    trace = [(step=0, action="start", unit=0, added=0, value=value, gain=0, size=length(C))]
    note(t) = (push!(trace, t); log === nothing || log(t))
    log === nothing || log(trace[1])
    while true
        flips = [u in C ? filter(!=(u), C) : sort(vcat(C, u)) for u in pool]
        values = objective(flips)
        gains = values .- value
        gain, k = findmax(gains)
        if gain >= min_gain
            action = pool[k] in C ? "remove" : "add"
            C = flips[k]; value = values[k]
            note((step=length(trace), action=action, unit=pool[k], added=0, value=value, gain=gain, size=length(C)))
            continue
        end
        inside = [i for i in eachindex(pool) if pool[i] in C]; outside = [i for i in eachindex(pool) if !(pool[i] in C)]
        (swap_top > 0 && !isempty(inside) && !isempty(outside)) || break
        best(ix) = sort(sort(ix; by=i -> (-gains[i], pool[i]))[1:min(swap_top, length(ix))])
        pairs = [(pool[i], pool[j]) for i in best(inside) for j in best(outside)]
        sets = [sort(vcat(filter(!=(u), C), r)) for (u, r) in pairs]
        values = objective(sets)
        gain, k = findmax(values .- value)
        gain >= min_gain || break
        C = sets[k]; value = values[k]
        note((step=length(trace), action="swap", unit=pairs[k][1], added=pairs[k][2], value=value, gain=gain, size=length(C)))
    end
    return C, value, trace
end

"Start set k of a multistart search: 0 is the pool, k >= 1 a seeded random subset."
function start_set(pool::Vector{Int}, sp, k::Integer)
    k == 0 && return sort(pool)
    1 <= k <= length(sp["start_seeds"]) || error("Start $k is not declared in the protocol")
    rng = MersenneTwister(sp["start_seeds"][k])
    C = [u for u in sort(pool) if rand(rng) < sp["start_keep_probability"]]
    isempty(C) && error("Empty start set $k")
    return C
end

multistart(protocol) = get(protocol["search"], "strategy", "floating_backward") == "multistart_flip_swap"

function benchmark_truth(S::Setup, sequence::AbstractString)
    t = Int8[c == '1' ? 1 : 0 for c in sequence]
    return Dict(f => t for f in keys(S.file_rows))
end

"Per-track random palindromic truths of the null (seeded, sorted track order)."
function null_truth(S::Setup, patterns, seed)
    rng = MersenneTwister(seed)
    pat = Dict(t => patterns[rand(rng, 1:length(patterns))] for t in sort(unique(values(S.file_track))))
    return Dict(f => Int8[c == '1' ? 1 : 0 for c in pat[S.file_track[f]]] for f in keys(S.file_rows))
end

"One search on `sessions`; writes the trace with the totals of every step."
function run_search(S::Setup, protocol, sessions, truth; label="", start=nothing)
    sp = protocol["search"]
    pool = pool_units(S, sessions)
    objective(sets) = [summarize(g, truth).correct for g in loso_grades(S, sets, sessions, truth)]
    t0 = time()
    logger(t) = (@printf("%s step %d %s unit=%d%s correct=%d gain=%d units=%d trainings=%d %.0f s\n", label, t.step,
                         t.action, t.unit, get(t, :added, 0) == 0 ? "" : "+$(t.added)", t.value, t.gain, t.size,
                         TRAININGS[], time() - t0); flush(stdout))
    if multistart(protocol)
        start === nothing && error("A multistart protocol needs --start")
        first_set = start_set(pool, sp, start)
        C, value, trace = local_search(pool, first_set, objective; min_gain=sp["min_gain"], swap_top=sp["swap_top"], log=logger)
    else
        start === nothing || error("--start is only defined for a multistart protocol")
        first_set = pool
        C, value, trace = floating_search(pool, objective; min_gain=sp["min_gain"], floating=sp["floating"], log=logger)
    end
    rows = Dict{String,String}[]
    current = copy(first_set)
    for t in trace
        t.action == "remove" && filter!(!=(t.unit), current)
        t.action == "add" && (current = sort(vcat(current, t.unit)))
        t.action == "swap" && (current = sort(vcat(filter(!=(t.unit), current), t.added)))
        s = summarize(only(loso_grades(S, [current], sessions, truth)), truth)
        push!(rows, Dict("step" => string(t.step), "action" => t.action, "unit" => string(t.unit),
            "unit_added" => get(t, :added, 0) == 0 ? "" : string(t.added),
            "unit_session" => t.unit == 0 ? "" : S.unit_session[t.unit],
            "unit_scans" => t.unit == 0 ? "" : string(length(unique(S.data.file[S.units[t.unit]]))),
            "gain" => string(t.gain), "units" => string(t.size), stat_strings(s)...))
    end
    current == C || error("Trace replay does not reach the search result")
    return C, pool, rows, first_set
end

stat_strings(s) = ["files" => string(s.files), "possible" => string(s.possible), "correct" => string(s.correct),
    "classified" => string(s.classified), "errors" => string(s.errors), "unclassified" => string(s.unclassified),
    "exact" => string(s.exact), "tracks" => string(s.tracks),
    "molecule_accuracy" => @sprintf("%.4f", s.molecule_accuracy), "tracks_majority_exact" => string(s.tracks_majority_exact)]
const STAT_COLUMNS = ["files", "possible", "correct", "classified", "errors", "unclassified", "exact", "tracks",
                      "molecule_accuracy", "tracks_majority_exact"]
const TRACE_COLUMNS = vcat(["step", "action", "unit", "unit_added", "unit_session", "unit_scans", "gain", "units"], STAT_COLUMNS)

file_rows_out(grades, extra) = [merge(Dict("file" => g.file, "session" => g.session, "track" => string(g.track),
    "N" => string(g.n), "correct" => string(g.correct), "classified" => string(g.classified),
    "exact" => string(g.exact)), extra) for g in grades]
const FILE_COLUMNS = ["file", "session", "track", "N", "correct", "classified", "exact"]

function run_nested(S::Setup, protocol, out, o)
    truth = benchmark_truth(S, protocol["data"]["control_sequence"])
    o in S.development || error("--outer must be a development session")
    inner = filter(!=(o), S.development)
    C, pool, trace, _ = run_search(S, protocol, inner, truth; label=o)
    write_table(joinpath(out, "trace.tsv"), TRACE_COLUMNS, trace)
    summary = Dict{String,String}[]; files = Dict{String,String}[]
    for (name, units) in (("pool", pool), ("selected", C))
        model, ok = frozen_model(S, units, inner)
        calls, _ = final_calls(S, model, S.session_rows[o])
        g = grade_session(S, o, calls, truth)
        push!(summary, Dict("outer" => o, "cohort" => name, "units" => string(length(units)),
            "removed" => join(setdiff(pool, units), ","), "emission_fit" => string(ok),
            "theta0" => @sprintf("%.6f", model.theta0), "theta1" => @sprintf("%.6f", model.theta1),
            stat_strings(summarize(g, truth))...))
        append!(files, file_rows_out(g, Dict("outer" => o, "cohort" => name)))
    end
    write_table(joinpath(out, "outer.tsv"), OUTER_COLUMNS, summary)
    write_table(joinpath(out, "outer_files.tsv"), vcat(["outer", "cohort"], FILE_COLUMNS), files)
    @printf("outer %s: pool %s / selected %s correct units\n", o, summary[1]["correct"], summary[2]["correct"])
end
const OUTER_COLUMNS = vcat(["outer", "cohort", "units", "removed", "emission_fit", "theta0", "theta1"], STAT_COLUMNS)

function run_final(S::Setup, protocol, out)
    truth = benchmark_truth(S, protocol["data"]["control_sequence"])
    C, pool, trace, _ = run_search(S, protocol, S.development, truth; label="final")
    write_table(joinpath(out, "trace.tsv"), TRACE_COLUMNS, trace)
    removed_at = Dict{Int,String}()
    for r in trace
        r["action"] == "remove" && (removed_at[parse(Int, r["unit"])] = r["step"])
        r["action"] == "add" && delete!(removed_at, parse(Int, r["unit"]))
    end
    write_table(joinpath(out, "cohort.tsv"), ["track", "session", "status", "removed_at_step", "corroborated_scans", "files"],
        [Dict("track" => string(u), "session" => S.unit_session[u], "status" => u in C ? "selected" : "removed",
              "removed_at_step" => get(removed_at, u, ""), "corroborated_scans" => string(length(unique(S.data.file[S.units[u]]))),
              "files" => join(unique(S.data.file[S.units[u]]), ",")) for u in pool])
    for (name, units) in (("pool", pool), ("selected", C))
        model, ok = frozen_model(S, units, S.development)
        ok || error("Emission rates could not be fitted for the $name model")
        write_frozen_model(joinpath(out, "model_$(name).tsv"), model)
        open(joinpath(out, "model_$(name).units"), "w") do io; foreach(u -> println(io, u), units); end
    end
end

function run_null(S::Setup, protocol, out, seed)
    np = protocol["null"]
    seed in np["seeds"] || error("Seed $seed is not declared in the protocol")
    C, _, trace, _ = run_search(S, protocol, S.development, null_truth(S, np["patterns"], seed); label="null$seed")
    write_table(joinpath(out, "trace.tsv"), TRACE_COLUMNS, trace)
end

"Merge one campaign directory: nested_<session>/, final/, null_<seed>/."
function run_report(S::Setup, protocol, out, dir)
    truth = benchmark_truth(S, protocol["data"]["control_sequence"])
    outer = Dict{String,String}[]; ofiles = Dict{String,String}[]
    for o in S.development
        d = joinpath(dir, "nested_$o")
        isfile(joinpath(d, "outer.tsv")) || error("Missing nested output for $o")
        append!(outer, last(read_grading(joinpath(d, "outer.tsv")))); append!(ofiles, last(read_grading(joinpath(d, "outer_files.tsv"))))
        # reference: the promoted transductive method of the run (learned on every scan of the run, label-free)
        g = grade_session(S, o, S.transductive[S.session_rows[o]], truth)
        push!(outer, Dict("outer" => o, "cohort" => "transductive_run", "units" => "", "removed" => "", "emission_fit" => "",
                          "theta0" => "", "theta1" => "", stat_strings(summarize(g, truth))...))
        append!(ofiles, file_rows_out(g, Dict("outer" => o, "cohort" => "transductive_run")))
    end
    write_table(joinpath(out, "outer.tsv"), OUTER_COLUMNS, outer)
    write_table(joinpath(out, "outer_files.tsv"), vcat(["outer", "cohort"], FILE_COLUMNS), ofiles)
    final = last(read_grading(joinpath(dir, "final", "trace.tsv")))
    nulls = Dict{String,String}[]
    for seed in protocol["null"]["seeds"]
        p = joinpath(dir, "null_$seed", "trace.tsv")
        isfile(p) || continue
        tr = last(read_grading(p))
        push!(nulls, Dict("seed" => string(seed), "units" => tr[end]["units"], "steps" => string(length(tr) - 1),
            "start_correct" => tr[1]["correct"], "end_correct" => tr[end]["correct"],
            "gain" => string(parse(Int, tr[end]["correct"]) - parse(Int, tr[1]["correct"])),
            "start_molecule_accuracy" => tr[1]["molecule_accuracy"], "end_molecule_accuracy" => tr[end]["molecule_accuracy"]))
    end
    write_table(joinpath(out, "null_summary.tsv"), ["seed", "units", "steps", "start_correct", "end_correct", "gain",
        "start_molecule_accuracy", "end_molecule_accuracy"], nulls)
    cohorts = ("pool", "selected", "transductive_run")
    tot(c, k) = sum(parse(Int, r[k]) for r in outer if r["cohort"] == c)
    gains = [parse(Int, r["gain"]) for r in nulls]
    fgain = parse(Int, final[end]["correct"]) - parse(Int, final[1]["correct"])
    open(joinpath(out, "report.md"), "w") do io
        println(io, "# Campaign B: reference-cohort search\n")
        println(io, "Development sessions: $(length(S.development)); pool: $(length(pool_units(S, S.development))) units.\n")
        println(io, "## Nested leave-one-session-out (honest estimate)\n")
        println(io, "Held-out sessions only. `transductive_run` is the promoted method of the run itself (learned without labels " *
                    "on all 493 scans, held-out session included), for reference.\n")
        println(io, "| cohort | correct | classified | errors | exact chains | molecules | molecule accuracy | molecules majority exact |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for c in cohorts
            m = molecule_stats([r for r in ofiles if r["cohort"] == c], truth)
            @printf(io, "| %s | %d/%d | %d | %d | %d | %d | %.3f | %d |\n", c, tot(c, "correct"), tot(c, "possible"),
                    tot(c, "classified"), tot(c, "errors"), tot(c, "exact"), m.tracks, m.accuracy, m.majority_exact)
        end
        println(io, "\nGain of the search on held-out sessions: $(tot("selected", "correct") - tot("pool", "correct")) correct units.\n")
        println(io, "| held-out session | files | pool | selected | transductive_run | units kept / pool |")
        println(io, "|---|---|---|---|---|---|")
        for o in S.development
            r = Dict(x["cohort"] => x for x in outer if x["outer"] == o)
            pool_n = parse(Int, r["pool"]["units"]); kept = parse(Int, r["selected"]["units"])
            println(io, "| $o | $(r["pool"]["files"]) | $(r["pool"]["correct"]) | $(r["selected"]["correct"]) | " *
                        "$(r["transductive_run"]["correct"]) | $kept / $pool_n |")
        end
        println(io, "\n## Final search ($(length(S.development)) sessions, in-sample objective)\n")
        println(io, "Start $(final[1]["correct"]) -> end $(final[end]["correct"]) correct units (gain $fgain); " *
                    "$(final[end]["units"]) of $(final[1]["units"]) units kept; per-molecule accuracy " *
                    "$(final[1]["molecule_accuracy"]) -> $(final[end]["molecule_accuracy"]).\n")
        println(io, "## Null searches (random palindromic truths)\n")
        isempty(gains) ? println(io, "No null output.") :
            println(io, "$(length(gains)) seeds; in-sample gain: min $(minimum(gains)), median $(sort(gains)[cld(length(gains), 2)]), " *
                        "max $(maximum(gains)); real gain $fgain ranks above $(count(<(fgain), gains)) of $(length(gains)) null gains.")
    end
    print(read(joinpath(out, "report.md"), String))
end

"One vote per physical molecule over file-level rows (tracks pooled across sessions)."
function molecule_stats(rows, truth)
    tracks = Dict{String,Vector{Dict{String,String}}}()
    for r in rows; push!(get!(tracks, r["track"], Dict{String,String}[]), r); end
    acc = [sum(parse(Int, r["correct"]) for r in v) / sum(length(truth[r["file"]]) for r in v) for v in values(tracks)]
    maj = count(v -> 2count(r -> r["exact"] == "true", v) > length(v), values(tracks))
    return (tracks=length(tracks), accuracy=isempty(acc) ? NaN : sum(acc) / length(acc), majority_exact=maj)
end

"One start of a multistart search on a nested, final or null target: trace and cohorts only."
function run_start(S::Setup, protocol, out, mode, o, k)
    if mode == "nested"
        o in S.development || error("--outer must be a development session")
        sessions = filter(!=(o), S.development); truth = benchmark_truth(S, protocol["data"]["control_sequence"]); label = "$(o)_$k"
    elseif mode == "final"
        sessions = S.development; truth = benchmark_truth(S, protocol["data"]["control_sequence"]); label = "final_$k"
    else
        seed = parse(Int, o)
        seed in protocol["null"]["seeds"] || error("Seed $seed is not declared in the protocol")
        sessions = S.development; truth = null_truth(S, protocol["null"]["patterns"], seed); label = "null$(seed)_$k"
    end
    C, _, trace, first_set = run_search(S, protocol, sessions, truth; label, start=k)
    write_table(joinpath(out, "trace.tsv"), TRACE_COLUMNS, trace)
    for (name, units) in (("start", first_set), ("end", C))
        open(joinpath(out, "$name.units"), "w") do io; foreach(u -> println(io, u), units); end
    end
end

"Best end of the starts <prefix>_0 ... <prefix>_n (ties: lowest start; start 0 is the pool)."
function best_start(dir, prefix, n)
    ends = Int[]; units = Vector{Int}[]; pool_value = 0
    for k in 0:n
        d = joinpath(dir, "$(prefix)_$k")
        isfile(joinpath(d, "end.units")) || error("Missing start output: $d")
        tr = last(read_grading(joinpath(d, "trace.tsv")))
        push!(ends, parse(Int, tr[end]["correct"]))
        push!(units, parse.(Int, filter(!isempty, readlines(joinpath(d, "end.units")))))
        k == 0 && (pool_value = parse(Int, tr[1]["correct"]))
    end
    v, i = findmax(ends)
    return (start=i - 1, units=units[i], value=v, ends=ends, pool_value=pool_value)
end

"Per-molecule paired difference b - a of correct units, with a bootstrap over molecules."
function paired_molecules(ofiles, a, b; draws=20000)
    sums = Dict{Tuple{String,String},Int}()
    for r in ofiles
        r["cohort"] in (a, b) || continue
        k = (r["track"], r["cohort"]); sums[k] = get(sums, k, 0) + parse(Int, r["correct"])
    end
    tracks = sort(unique(first.(collect(keys(sums)))))
    d = [sums[(t, b)] - sums[(t, a)] for t in tracks]
    rng = MersenneTwister(1)
    boot = sort([sum(d[rand(rng, 1:length(d), length(d))]) for _ in 1:draws])
    return (molecules=length(d), better=count(>(0), d), worse=count(<(0), d), same=count(==(0), d), total=sum(d),
            lo=boot[round(Int, 0.025draws)], hi=boot[round(Int, 0.975draws)],
            largest=sort(collect(zip(tracks, d)); by=x -> -abs(x[2]))[1:min(5, length(d))])
end

const OUTER2_COLUMNS = vcat(["outer", "cohort", "start", "in_sample", "in_sample_pool", "starts_end", "units", "removed",
                             "emission_fit", "theta0", "theta1"], STAT_COLUMNS)

"Merge a multistart campaign: nested_<session>_<k>/, final_<k>/, null_<seed>_<k>/."
function run_report_multistart(S::Setup, protocol, out, dir)
    truth = benchmark_truth(S, protocol["data"]["control_sequence"])
    n = length(protocol["search"]["start_seeds"])
    outer = Dict{String,String}[]; ofiles = Dict{String,String}[]
    for o in S.development
        b = best_start(dir, "nested_$o", n)
        inner = filter(!=(o), S.development); pool = pool_units(S, inner)
        for (name, units) in (("pool", pool), ("selected", b.units))
            model, ok = frozen_model(S, units, inner)
            calls, _ = final_calls(S, model, S.session_rows[o])
            g = grade_session(S, o, calls, truth)
            push!(outer, Dict("outer" => o, "cohort" => name, "start" => name == "pool" ? "" : string(b.start),
                "in_sample" => string(name == "pool" ? b.pool_value : b.value), "in_sample_pool" => string(b.pool_value),
                "starts_end" => join(b.ends, ","), "units" => string(length(units)), "removed" => join(setdiff(pool, units), ","),
                "emission_fit" => string(ok), "theta0" => @sprintf("%.6f", model.theta0), "theta1" => @sprintf("%.6f", model.theta1),
                stat_strings(summarize(g, truth))...))
            append!(ofiles, file_rows_out(g, Dict("outer" => o, "cohort" => name)))
        end
        g = grade_session(S, o, S.transductive[S.session_rows[o]], truth)
        push!(outer, Dict("outer" => o, "cohort" => "transductive_run", stat_strings(summarize(g, truth))...))
        append!(ofiles, file_rows_out(g, Dict("outer" => o, "cohort" => "transductive_run")))
    end
    write_table(joinpath(out, "outer.tsv"), OUTER2_COLUMNS, outer)
    write_table(joinpath(out, "outer_files.tsv"), vcat(["outer", "cohort"], FILE_COLUMNS), ofiles)
    f = best_start(dir, "final", n); pool = pool_units(S, S.development)
    write_table(joinpath(out, "cohort.tsv"), ["track", "session", "status", "corroborated_scans"],
        [Dict("track" => string(u), "session" => S.unit_session[u], "status" => u in f.units ? "selected" : "removed",
              "corroborated_scans" => string(length(unique(S.data.file[S.units[u]])))) for u in pool])
    for (name, units) in (("pool", pool), ("selected", f.units))
        model, ok = frozen_model(S, units, S.development)
        ok || error("Emission rates could not be fitted for the $name model")
        write_frozen_model(joinpath(out, "model_$(name).tsv"), model)
        open(joinpath(out, "model_$(name).units"), "w") do io; foreach(u -> println(io, u), units); end
    end
    nulls = Dict{String,String}[]
    for seed in protocol["null"]["seeds"]
        b = best_start(dir, "null_$seed", n)
        push!(nulls, Dict("seed" => string(seed), "start" => string(b.start), "pool_correct" => string(b.pool_value),
                          "best_correct" => string(b.value), "gain" => string(b.value - b.pool_value), "starts_end" => join(b.ends, ",")))
    end
    write_table(joinpath(out, "null_summary.tsv"), ["seed", "start", "pool_correct", "best_correct", "gain", "starts_end"], nulls)
    tot(c, k) = sum(parse(Int, r[k]) for r in outer if r["cohort"] == c)
    pm = paired_molecules(ofiles, "pool", "selected")
    gains = [parse(Int, r["gain"]) for r in nulls]; fgain = f.value - f.pool_value
    dec = get(protocol, "decision", Dict())
    open(joinpath(out, "report.md"), "w") do io
        println(io, "# Campaign B, option 1: multistart flip/swap search\n")
        println(io, "Development sessions: $(length(S.development)); pool: $(length(pool)) units; starts: pool + $n random.\n")
        println(io, "## Nested leave-one-session-out (honest estimate)\n")
        println(io, "| cohort | correct | classified | errors | exact chains | molecules | molecule accuracy | molecules majority exact |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for c in ("pool", "selected", "transductive_run")
            m = molecule_stats([r for r in ofiles if r["cohort"] == c], truth)
            @printf(io, "| %s | %d/%d | %d | %d | %d | %d | %.3f | %d |\n", c, tot(c, "correct"), tot(c, "possible"),
                    tot(c, "classified"), tot(c, "errors"), tot(c, "exact"), m.tracks, m.accuracy, m.majority_exact)
        end
        gain = tot("selected", "correct") - tot("pool", "correct")
        println(io, "\nHeld-out gain over the pool: $gain correct units. Per molecule (selected - pool): $(pm.better) better, " *
                    "$(pm.worse) worse, $(pm.same) unchanged of $(pm.molecules); bootstrap 95% interval over molecules " *
                    "[$(pm.lo), $(pm.hi)]; largest changes $(pm.largest).\n")
        println(io, "| held-out session | files | pool | selected | transductive_run | chosen start | in-sample pool -> best | ends per start |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for o in S.development
            r = Dict(x["cohort"] => x for x in outer if x["outer"] == o)
            println(io, "| $o | $(r["pool"]["files"]) | $(r["pool"]["correct"]) | $(r["selected"]["correct"]) | " *
                        "$(r["transductive_run"]["correct"]) | $(r["selected"]["start"]) | $(r["pool"]["in_sample"]) -> " *
                        "$(r["selected"]["in_sample"]) | $(r["selected"]["starts_end"]) |")
        end
        println(io, "\n## Final search ($(length(S.development)) sessions, in-sample objective)\n")
        println(io, "Pool $(f.pool_value) -> best $(f.value) correct units (gain $fgain, start $(f.start)); ends per start " *
                    "$(join(f.ends, ", ")); $(length(f.units)) of $(length(pool)) units kept.\n")
        println(io, "## Null searches (random palindromic truths, same multistart search)\n")
        println(io, isempty(gains) ? "No null output." :
            "$(length(gains)) seeds; in-sample gains $(join(gains, ", ")); the real gain $fgain ranks above " *
            "$(count(<(fgain), gains)) of $(length(gains)).")
        if !isempty(dec)
            println(io, "\n## Pre-registered decision (protocol [decision])\n")
            c1 = gain > dec["v1_nested_gain"]; c2 = pm.lo > 0; c3 = pm.better > pm.worse
            println(io, "1. Held-out gain $gain > v1 gain $(dec["v1_nested_gain"]): $(c1 ? "yes" : "no")")
            println(io, "2. Bootstrap lower bound $(pm.lo) > 0: $(c2 ? "yes" : "no")")
            println(io, "3. More molecules better ($(pm.better)) than worse ($(pm.worse)): $(c3 ? "yes" : "no")")
            println(io, "\nThe stronger search supports the hypothesis that a better optimum generalises: $(c1 && c2 && c3 ? "yes" : "no").")
        end
    end
    print(read(joinpath(out, "report.md"), String))
end

"Grading-side TSV reader (the method reader refuses label-like columns)."
function read_grading(path)
    lines = filter(!isempty, readlines(path))
    header = String.(split(lines[1], '\t'))
    return header, [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

function run_sealed(S::Setup, protocol, out, models)
    truth = benchmark_truth(S, protocol["data"]["control_sequence"])
    summary = Dict{String,String}[]; files = Dict{String,String}[]
    for name in ("pool", "selected")
        model = read_frozen_model(joinpath(models, "model_$(name).tsv"))
        grades = Any[]; preds = Dict{String,String}[]
        for s in S.sealed
            rows = S.session_rows[s]
            calls, p = final_calls(S, model, rows)
            append!(grades, grade_session(S, s, calls, truth))
            for (k, r) in enumerate(rows)
                c(x) = x < 0 ? "?" : string(x)
                push!(preds, Dict("file" => S.data.file[r], "lobe" => string(S.data.lobe[r]), "predicted" => c(calls[k]),
                    "scan_predicted" => c(p.call[k]), "probability_1" => isfinite(p.p_vote[k]) ? @sprintf("%.8f", p.p_vote[k]) : "NA",
                    "physical_lobe" => S.physical[r], "model" => S.cfg["model"]["name"] * "+" * name))
            end
        end
        sort!(preds; by=r -> (r["file"], parse(Int, r["lobe"])))
        write_table(joinpath(out, "predictions_$(name).tsv"),
                    ["file", "lobe", "predicted", "scan_predicted", "probability_1", "physical_lobe", "model"], preds)
        push!(summary, Dict("cohort" => name, "theta0" => @sprintf("%.6f", model.theta0), "theta1" => @sprintf("%.6f", model.theta1),
                            "pi" => @sprintf("%.3f", model.pi), stat_strings(summarize(grades, truth))...))
        append!(files, file_rows_out(grades, Dict("cohort" => name)))
    end
    grades = reduce(vcat, [grade_session(S, s, S.transductive[S.session_rows[s]], truth) for s in S.sealed])
    push!(summary, Dict("cohort" => "transductive_run", "theta0" => "", "theta1" => "", "pi" => "",
                        stat_strings(summarize(grades, truth))...))
    append!(files, file_rows_out(grades, Dict("cohort" => "transductive_run")))
    write_table(joinpath(out, "sealed.tsv"), vcat(["cohort", "theta0", "theta1", "pi"], STAT_COLUMNS), summary)
    write_table(joinpath(out, "sealed_files.tsv"), vcat(["cohort"], FILE_COLUMNS), files)
end

function search_main(args=ARGS)
    o = search_options(args)
    o === nothing && return println("search_reference_cohort.jl --protocol TOML --mode nested|final|null|report|sealed --outdir NEW_DIR " *
                                    "[--outer SESSION] [--seed K] [--start K] [--campaign DIR] [--models DIR]")
    BLAS.set_num_threads(1)
    protocol = TOML.parsefile(o["--protocol"])
    S = load_setup(protocol)
    out = o["--outdir"]; mkpath(out)
    cp(o["--protocol"], joinpath(out, "protocol.toml"))
    @printf("Julia %s, %d threads; development %d sessions, sealed %d; pool %d units\n", VERSION, Threads.nthreads(),
            length(S.development), length(S.sealed), length(pool_units(S, S.development)))
    mode = o["--mode"]
    if multistart(protocol) && mode in ("nested", "final", "null")
        haskey(o, "--start") || error("A multistart protocol needs --start K")
        run_start(S, protocol, out, mode, mode == "nested" ? o["--outer"] : get(o, "--seed", ""), parse(Int, o["--start"]))
    else
        haskey(o, "--start") && error("--start is only defined for a multistart protocol")
        mode == "nested" && run_nested(S, protocol, out, o["--outer"])
        mode == "final" && run_final(S, protocol, out)
        mode == "null" && run_null(S, protocol, out, parse(Int, o["--seed"]))
        mode == "report" && (multistart(protocol) ? run_report_multistart : run_report)(S, protocol, out, o["--campaign"])
    end
    mode == "sealed" && run_sealed(S, protocol, out, o["--models"])
    println("Done: $(TRAININGS[]) trainings -> $out")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && search_main()

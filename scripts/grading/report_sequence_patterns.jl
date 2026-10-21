#!/usr/bin/env julia
# Sequence patterns of an application run (10-20mers): one sequence per
# molecule, structure statistics against a composition-preserving null, and a
# comparison with manual readings. Grading side: reads the manual readings.
#
# Usage:
#   julia --project=. scripts/grading/report_sequence_patterns.jl --run RUN_DIR \
#       --subsets TSV --manual benchmarks/chitosan_10_20mer_manual_sequences.tsv \
#       --config config/unit_assignment_frozen_labelfree.toml --outdir NEW_DIR \
#       [--campaign-gap-days 7] [--permutations 10000] [--model NAME=MODEL.tsv ...]
#
# Assignments compared:
#   transductive  the run's own final (fused) calls (promoted method);
#   frozen        the fully frozen label-free model built per measurement
#                 campaign (src/frozen_assignment.jl, config reference rule);
#                 a campaign is a run of sessions whose acquisition dates are
#                 at most --campaign-gap-days apart;
#   NAME          (--model, repeatable) a frozen model built elsewhere, for
#                 example on the 6-mer benchmark, applied unchanged.
# A molecule is a consensus track; its sequence takes, for every physical lobe
# (rank along the chain axis), the majority of the final calls of its scans
# ("?" if none or tied). Orientation is unknown, so every statistic is either
# orientation-free or takes the better orientation, and every null applies the
# same procedure to sequences shuffled within each chain (composition kept).
# Statistics: GlcNAc (1) fraction; blockiness P = f(01+10) / (2 f0 f1)
# (1 Bernoulli, < 1 blocks, > 1 alternation); mean run length; ends vs
# interior; longest alternating end segment; agreement with the manual
# template "1010" + "1"s; profile along the chain oriented by its longer
# alternating end. --subsets (file, subset, session, date dd.mm.yyyy)
# describes the run's scans (subset, for example selection / short_chain /
# target; a track takes its scans' majority subset). Manual comparison: N and units (better of the two
# orientations, unreadable manual units excluded) for the scans read by hand.
using STMFit.FrozenAssignment, STMFit.ReconstructedUnitAssignment
using Dates, LinearAlgebra, Printf, Random, Statistics

const TEMPLATE_HEAD = Int8[1, 0, 1, 0]

function pattern_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--campaign-gap-days" => "7", "--permutations" => "10000", "--model" => ""); i = 1
    while i <= length(args)
        args[i] in ("--run", "--subsets", "--manual", "--config", "--outdir", "--campaign-gap-days", "--permutations", "--model") &&
            i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i] == "--model" && !isempty(o["--model"]) ? o["--model"] * ";" * args[i+1] : args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--subsets", "--manual", "--config", "--outdir")) ||
        error("Required: --run --subsets --manual --config --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

"Grading-side TSV reader (comment lines skipped)."
function tsv(path)
    lines = filter(l -> !isempty(l) && !startswith(l, "#"), readlines(path))
    h = String.(split(lines[1], '\t'))
    return [Dict(zip(h, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

call(x) = x == "?" ? Int8(-1) : parse(Int8, x)
seqstr(s) = join(c < 0 ? "?" : string(c) for c in s)

"Campaign index per session: sessions sorted by date, a gap > gap_days starts a new campaign."
function campaigns(session_date, gap_days)
    sessions = sort(collect(keys(session_date)); by=s -> (session_date[s], s))
    out = Dict{String,Int}(); c = 0; prev = nothing
    for s in sessions
        (prev === nothing || Dates.value(session_date[s] - prev) > gap_days) && (c += 1)
        out[s] = c; prev = session_date[s]
    end
    return out
end

"Track sequences from final calls per row: majority per physical lobe (rank order)."
function track_sequences(data, physical, calls)
    votes = Dict{String,Dict{Int,Vector{Int8}}}()
    for i in eachindex(data.file)
        k = physical[i]; isempty(k) && continue
        t, r = rsplit(k, "_"; limit=2)
        push!(get!(get!(votes, t, Dict{Int,Vector{Int8}}()), parse(Int, r), Int8[]), calls[i])
    end
    out = Dict{String,Vector{Int8}}()
    for (t, byrank) in votes
        n = maximum(keys(byrank))
        out[t] = [begin
                      v = filter(>=(0), get(byrank, r, Int8[]))
                      ones = count(==(1), v); zeros = length(v) - ones
                      isempty(v) || ones == zeros ? Int8(-1) : Int8(ones > zeros ? 1 : 0)
                  end for r in 1:n]
    end
    return out
end

alt_prefix(s) = (k = 1; while k < length(s) && s[k] != s[k+1]; k += 1; end; k)
alt_end(s) = max(alt_prefix(s), alt_prefix(reverse(s)))
template(n) = vcat(TEMPLATE_HEAD[1:min(4, n)], fill(Int8(1), max(0, n - 4)))
agree(s, t) = count(i -> s[i] == t[i], eachindex(s)) / length(s)
template_agreement(s) = max(agree(s, template(length(s))), agree(reverse(s), template(length(s))))
oriented(s) = alt_prefix(reverse(s)) > alt_prefix(s) ? reverse(s) : s

"Statistics over complete sequences (no unknown unit), each returning one number."
function seq_stats(seqs)
    isempty(seqs) && return Dict{String,Float64}()
    all1 = reduce(vcat, seqs)
    f1 = mean(all1 .== 1)
    pairs = [(s[i], s[i+1]) for s in seqs for i in 1:length(s)-1]
    mixed = count(p -> p[1] != p[2], pairs) / max(length(pairs), 1)
    runs = [count(i -> i == 1 || s[i] != s[i-1], eachindex(s)) for s in seqs]
    ends = [x for s in seqs for x in (s[1], s[end])]
    interior = reduce(vcat, [s[2:end-1] for s in seqs if length(s) > 2]; init=Int8[])
    return Dict("glcnac_fraction" => f1,
                "blockiness_P" => f1 in (0.0, 1.0) ? NaN : mixed / (2 * f1 * (1 - f1)),
                "mean_run_length" => mean(length.(seqs) ./ runs),
                "ends_glcnac" => mean(ends .== 1), "interior_glcnac" => isempty(interior) ? NaN : mean(interior .== 1),
                "alternating_end" => mean(alt_end.(seqs)),
                "template_agreement" => mean(template_agreement.(seqs)))
end

"Two-sided permutation p-values: shuffle within each chain, recompute every statistic."
function permutation_test(seqs, nperm, rng)
    obs = seq_stats(seqs)
    names = sort(collect(keys(obs)))
    null = Dict(k => Float64[] for k in names)
    for _ in 1:nperm
        st = seq_stats([shuffle(rng, s) for s in seqs])
        for k in names; push!(null[k], st[k]); end
    end
    return [(stat=k, observed=obs[k], null_mean=mean(filter(isfinite, null[k])),
             p_low=mean(filter(isfinite, null[k]) .<= obs[k]), p_high=mean(filter(isfinite, null[k]) .>= obs[k])) for k in names]
end

"Fraction of GlcNAc by position (from the oriented end), positions 1..maxpos."
function profile(seqs, maxpos)
    os = oriented.(seqs)
    [(pos=p, n=count(s -> length(s) >= p, os),
      glcnac=(v = [s[p] for s in os if length(s) >= p]; isempty(v) ? NaN : mean(v))) for p in 1:maxpos]
end

"Units of a manual reading against calls: better orientation, unreadable manual units excluded."
function manual_agreement(calls, manual)
    m = [c == '?' ? Int8(-1) : parse(Int8, c) for c in manual]
    best = (agree=0, compared=0)
    for c in (calls, reverse(calls))
        n = min(length(c), length(m))
        idx = [i for i in 1:n if m[i] >= 0 && c[i] >= 0]
        a = count(i -> c[i] == m[i], idx)
        (isempty(idx) ? 0.0 : a / length(idx)) > (best.compared == 0 ? -1.0 : best.agree / best.compared) &&
            (best = (agree=a, compared=length(idx)))
    end
    return best
end

function pattern_main(args=ARGS)
    o = pattern_options(args)
    o === nothing && return println("report_sequence_patterns.jl --run DIR --subsets TSV --manual TSV --config TOML --outdir NEW_DIR [--campaign-gap-days 7] [--permutations 10000]")
    BLAS.set_num_threads(1)
    run = o["--run"]; out = o["--outdir"]; mkpath(out)
    cfg = load_config(o["--config"]); st = load_frozen_settings(cfg)
    R = load_run_assignment(run, cfg)
    meta = Dict(r["file"] => r for r in tsv(o["--subsets"]))
    subset = Dict(f => r["subset"] for (f, r) in meta)
    cons = Dict(r["file"] => r for r in tsv(joinpath(run, "consensus", "consensus.tsv")))
    final_n = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in R.consensus)
    _, fused = lobe_table(joinpath(run, "predictions.tsv"); required=["predicted"])
    trans = [call(fused[(R.data.file[i], R.data.lobe[i])]["predicted"]) for i in eachindex(R.data.file)]
    trans_scan = [call(fused[(R.data.file[i], R.data.lobe[i])]["scan_predicted"]) for i in eachindex(R.data.file)]
    # sessions -> campaigns (acquisition dates from --subsets)
    files = sort(unique(R.data.file))
    all(f -> haskey(meta, f), files) || error("--subsets misses scans of the run")
    session = Dict(f => meta[f]["session"] for f in files)
    date = Dict(f => Date(meta[f]["date"], dateformat"dd.mm.yyyy") for f in files)
    sdate = Dict{String,Date}()
    for f in files; sdate[session[f]] = min(get(sdate, session[f], date[f]), date[f]); end
    camp = campaigns(sdate, parse(Int, o["--campaign-gap-days"]))
    # frozen model per campaign
    frozen = fill(Int8(-1), length(R.data.file)); cinfo = Dict{String,String}[]
    for c in sort(unique(values(camp)))
        cf = Set(f for f in files if camp[session[f]] == c)
        rows = findall(in(cf), R.data.file)
        ref = reference_files([r for r in R.consensus if basename(r["filepath"]) in cf], st.max_per_track)
        status = "ok"; th = (NaN, NaN)
        try
            model = train_reference_model(R.data, R.fisher, findall(in(ref), R.data.file), rows, R.physical, cfg)
            f, _ = apply_reference_model(model, R.data, R.fisher, rows, R.physical, cfg)
            frozen[rows] .= f; th = (model.theta0, model.theta1)
        catch err
            err isa InterruptException && rethrow()
            status = "failed: " * first(replace(sprint(showerror, err), '\n' => ' '), 120)
        end
        ss = sort(unique(session[f] for f in cf))
        push!(cinfo, Dict("campaign" => string(c), "sessions" => join(ss, ","), "first_date" => string(minimum(sdate[s] for s in ss)),
                          "scans" => string(length(cf)), "reference_scans" => string(length(ref)), "status" => status,
                          "theta0" => @sprintf("%.4f", th[1]), "theta1" => @sprintf("%.4f", th[2])))
    end
    write_table(joinpath(out, "campaigns.tsv"), ["campaign", "sessions", "first_date", "scans", "reference_scans", "status", "theta0", "theta1"], cinfo)
    # molecule sequences
    seqs = Dict("transductive" => track_sequences(R.data, R.physical, trans), "frozen" => track_sequences(R.data, R.physical, frozen))
    extra = Dict{String,Vector{Int8}}()
    for spec in split(o["--model"], ";"; keepempty=false)
        name, path = split(spec, "="; limit=2)
        f, _ = apply_reference_model(read_frozen_model(path), R.data, R.fisher, collect(eachindex(R.data.file)), R.physical, cfg)
        extra[name] = f; seqs[name] = track_sequences(R.data, R.physical, f)
    end
    assignments = vcat(["transductive", "frozen"], sort(collect(keys(extra))))
    trackfiles = Dict{String,Vector{String}}()
    for f in files; push!(get!(trackfiles, "track" * cons[f]["track"], String[]), f); end
    majority(v) = (c = Dict{String,Int}(); for x in v; c[x] = get(c, x, 0) + 1; end; first(sort(collect(c); by=p -> (-p[2], p[1])))[1])
    tracks = sort(collect(keys(seqs["transductive"])); by=t -> parse(Int, t[6:end]))
    trows = [Dict("track" => t, "scans" => string(length(trackfiles[t])), "subset" => majority([get(subset, f, "?") for f in trackfiles[t]]),
                  "campaign" => string(camp[session[first(trackfiles[t])]]), "N" => string(length(seqs["transductive"][t])),
                  "files" => join(trackfiles[t], ","), "seq_transductive" => seqstr(seqs["transductive"][t]),
                  "seq_frozen" => seqstr(get(seqs["frozen"], t, Int8[])),
                  ["seq_" * a => seqstr(get(seqs[a], t, Int8[])) for a in keys(extra)]...) for t in tracks]
    write_table(joinpath(out, "tracks.tsv"), vcat(["track", "scans", "subset", "campaign", "N"], "seq_" .* assignments, ["files"]), trows)
    # statistics with permutation nulls, per assignment and molecule group
    rng = MersenneTwister(1); nperm = parse(Int, o["--permutations"])
    srows = Dict{String,String}[]
    groups = [("all", t -> true), ("N>=7", t -> parse(Int, t["N"]) >= 7), ("selection+target", t -> t["subset"] in ("selection", "target")),
              ("short_chain", t -> t["subset"] == "short_chain")]
    for a in assignments, (gname, keep) in groups
        S = [seqs[a][t["track"]] for t in trows if keep(t) && haskey(seqs[a], t["track"]) && length(seqs[a][t["track"]]) >= 4 &&
             all(>=(0), seqs[a][t["track"]])]
        length(S) >= 3 || continue
        for r in permutation_test(S, nperm, rng)
            push!(srows, Dict("assignment" => a, "group" => gname, "molecules" => string(length(S)), "stat" => r.stat,
                              "observed" => @sprintf("%.4f", r.observed), "null_mean" => @sprintf("%.4f", r.null_mean),
                              "p_low" => @sprintf("%.4f", r.p_low), "p_high" => @sprintf("%.4f", r.p_high)))
        end
    end
    write_table(joinpath(out, "stats.tsv"), ["assignment", "group", "molecules", "stat", "observed", "null_mean", "p_low", "p_high"], srows)
    prows = Dict{String,String}[]
    for a in assignments
        S = [seqs[a][t["track"]] for t in trows if haskey(seqs[a], t["track"]) && length(seqs[a][t["track"]]) >= 4 && all(>=(0), seqs[a][t["track"]])]
        for p in profile(S, 20)
            push!(prows, Dict("assignment" => a, "position" => string(p.pos), "molecules" => string(p.n), "glcnac" => @sprintf("%.3f", p.glcnac)))
        end
    end
    write_table(joinpath(out, "profile.tsv"), ["assignment", "position", "molecules", "glcnac"], prows)
    # manual readings
    mrows = Dict{String,String}[]
    for m in tsv(o["--manual"])
        f = m["file"]
        row = Dict("file" => f, "status" => m["status"], "manual" => m["sequence"], "manual_N" => string(length(m["sequence"])),
                   "in_run" => string(f in files))
        if f in files
            idx = findall(==(f), R.data.file)
            t = "track" * cons[f]["track"]
            row["N_scan"] = cons[f]["N_scan"]; row["N_final"] = string(final_n[f]); row["track"] = t
            row["subset"] = get(subset, f, "?")
            for (name, c) in vcat([("scan", trans_scan[idx]), ("transductive", trans[idx]), ("frozen", frozen[idx])],
                                  [(a, extra[a][idx]) for a in sort(collect(keys(extra)))])
                row["calls_" * name] = seqstr(c)
                if m["status"] == "read"
                    a = manual_agreement(c, m["sequence"]); row["agree_" * name] = "$(a.agree)/$(a.compared)"
                end
            end
            row["track_seq_transductive"] = seqstr(seqs["transductive"][t])
        end
        push!(mrows, row)
    end
    mcols = vcat(["file", "status", "manual", "manual_N", "in_run", "subset", "N_scan", "N_final", "track", "calls_scan", "agree_scan",
                  "calls_transductive", "agree_transductive", "calls_frozen", "agree_frozen"],
                 vcat([["calls_" * a, "agree_" * a] for a in sort(collect(keys(extra)))]...), ["track_seq_transductive"])
    write_table(joinpath(out, "manual_comparison.tsv"), mcols, mrows)
    # report
    open(joinpath(out, "report.md"), "w") do io
        println(io, "# Sequence patterns: $(basename(abspath(run)))\n")
        println(io, "$(length(files)) scans, $(length(tracks)) molecules (consensus tracks). Campaigns (gap > $(o["--campaign-gap-days"]) days):\n")
        println(io, "| campaign | first date | scans | reference scans | frozen model |"); println(io, "|---|---|---|---|---|")
        for c in cinfo; println(io, "| $(c["campaign"]) | $(c["first_date"]) | $(c["scans"]) | $(c["reference_scans"]) | $(c["status"]) |"); end
        Ns = [parse(Int, t["N"]) for t in trows]
        println(io, "\nMolecule lengths N: median $(median(Ns)), range $(minimum(Ns))–$(maximum(Ns)); N ≥ 7: $(count(>=(7), Ns)) molecules.\n")
        println(io, "## Structure against the composition-preserving null ($nperm shuffles)\n")
        println(io, "| assignment | group | molecules | statistic | observed | null mean | p (obs ≤ null) | p (obs ≥ null) |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in srows
            println(io, "| $(r["assignment"]) | $(r["group"]) | $(r["molecules"]) | $(r["stat"]) | $(r["observed"]) | $(r["null_mean"]) | $(r["p_low"]) | $(r["p_high"]) |")
        end
        println(io, "\n## Manual readings\n")
        println(io, "| file | manual | N final | per-scan calls | agree | final (transductive) | agree | frozen | agree |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        for r in mrows
            r["status"] == "read" && r["in_run"] == "true" || continue
            println(io, "| $(r["file"]) | $(r["manual"]) | $(r["N_final"]) | $(r["calls_scan"]) | $(r["agree_scan"]) | $(r["calls_transductive"]) | " *
                        "$(r["agree_transductive"]) | $(r["calls_frozen"]) | $(r["agree_frozen"]) |")
        end
    end
    print(read(joinpath(out, "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && pattern_main()

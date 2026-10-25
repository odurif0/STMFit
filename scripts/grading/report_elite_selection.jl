#!/usr/bin/env julia
# Elite subset from the frozen image-quality score (campaign E). Grading side: it
# reads the earlier human triage as weak labels; the score itself reads none.
#
# Usage:
#   julia --project=. scripts/grading/report_elite_selection.jl --score SCORE_TSV \
#       --triage benchmarks/chitosan_6mer_cu100_triage.tsv \
#       --quality benchmarks/chitosan_6mer_cu100_quality.tsv --outdir NEW_DIR [--precision 0.98]
#
# Scans are matched by sha256. Human kept = non-empty triage set; human rejected =
# empty set (overview rejections have no score and drop out). Reports:
#   * the AUC of the score against the triage;
#   * the leave-one-session-out precision: for each session, the threshold is the
#     smallest score above which at least --precision of the training scans are
#     human-kept, then applied to the held-out session;
#   * the final threshold on all sessions, and the elite subset = scans of the
#     full set (quality table, in_full = yes) with score >= threshold.
# Writes elite.tsv (file, session, sha256, score, elite), by_session.tsv, report.md.
using STMFit.ReconstructedUnitAssignment: write_table
using Printf, Statistics

function elite_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--precision" => "0.98"); i = 1
    while i <= length(args)
        args[i] in ("--score", "--triage", "--quality", "--outdir", "--precision") && i < length(args) || error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--score", "--triage", "--quality", "--outdir")) || error("Required: --score --triage --quality --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function tsv(path)
    lines = [l for l in readlines(path) if !isempty(strip(l)) && !startswith(l, "#")]
    header = String.(split(lines[1], '\t'; keepempty=true))
    return [Dict(zip(header, String.(split(l, '\t'; keepempty=true)))) for l in lines[2:end]]
end

"Smallest score t such that at least `prec` of the scans scoring >= t are human-kept (NaN if none)."
function threshold(scores, kept, prec)
    idx = sortperm(scores; rev=true)
    best = NaN; ck = 0
    for (rank, i) in enumerate(idx)
        ck += kept[i]
        ck / rank >= prec && (best = scores[i])
    end
    return best
end

"Mann-Whitney AUC: probability that a kept scan outscores a rejected one (ties count half)."
function auc(scores, kept)
    pos = scores[kept .== 1]; neg = scores[kept .== 0]
    (isempty(pos) || isempty(neg)) && return NaN
    return mean((p > n) + 0.5 * (p == n) for p in pos, n in neg)
end

function elite_main(args=ARGS)
    o = elite_options(args)
    o === nothing && return println("report_elite_selection.jl --score TSV --triage TSV --quality TSV --outdir NEW_DIR [--precision 0.98]")
    prec = parse(Float64, o["--precision"])
    score = Dict(r["sha256"] => r for r in tsv(o["--score"]) if isfinite(something(tryparse(Float64, r["score"]), NaN)))
    rows = NamedTuple[]
    for t in tsv(o["--triage"])
        haskey(score, t["sha256"]) || continue
        push!(rows, (file=t["file"], session=t["session"], sha=t["sha256"], kept=isempty(t["set"]) ? 0 : 1,
                     score=parse(Float64, score[t["sha256"]]["score"])))
    end
    s = [r.score for r in rows]; k = [r.kept for r in rows]; sess = [r.session for r in rows]
    a = auc(s, k)
    t_all = threshold(s, k, prec)
    # leave-one-session-out
    pooled_n = 0; pooled_kept = 0; by_session = Dict{String,String}[]
    for session in sort(unique(sess))
        train = sess .!= session; test = .!train
        t = threshold(s[train], k[train], prec)
        sel = isnan(t) ? falses(sum(test)) : s[test] .>= t
        n = count(sel); kk = sum(k[test][sel]); pooled_n += n; pooled_kept += kk
        push!(by_session, Dict("session" => session, "scored" => string(sum(test)), "human_kept" => string(sum(k[test])),
            "threshold" => isnan(t) ? "NaN" : @sprintf("%.3f", t), "predicted_elite" => string(n), "of_which_human_kept" => string(kk)))
    end
    quality = Dict(r["sha256"] => r for r in tsv(o["--quality"]))
    out = Dict{String,String}[]
    for r in rows
        q = get(quality, r.sha, nothing)
        q !== nothing && q["in_full"] == "yes" || continue
        push!(out, Dict("file" => r.file, "session" => r.session, "sha256" => r.sha, "score" => @sprintf("%.4f", r.score),
                        "elite" => !isnan(t_all) && r.score >= t_all ? "yes" : "no"))
    end
    mkpath(o["--outdir"])
    write_table(joinpath(o["--outdir"], "elite.tsv"), ["file", "session", "sha256", "score", "elite"], out)
    write_table(joinpath(o["--outdir"], "by_session.tsv"), ["session", "scored", "human_kept", "threshold", "predicted_elite", "of_which_human_kept"], by_session)
    n_elite = count(r -> r["elite"] == "yes", out)
    open(joinpath(o["--outdir"], "report.md"), "w") do io
        println(io, "# Elite subset from the frozen quality score\n")
        println(io, "$(length(rows)) scored scans with a triage decision ($(sum(k)) human-kept, $(length(k) - sum(k)) rejected).\n")
        println(io, "- AUC of the score against the triage: **$(@sprintf("%.3f", a))**.")
        println(io, "- Target precision: $prec (share of human-kept among scans above the threshold).")
        println(io, "- Final threshold (all sessions): **$(@sprintf("%.3f", t_all))**.")
        println(io, "- Leave-one-session-out: $pooled_n scans predicted elite, $pooled_kept human-kept (precision ",
                @sprintf("%.3f", pooled_n == 0 ? NaN : pooled_kept / pooled_n), ").")
        println(io, "- **Elite subset: $n_elite of the $(length(out)) full-set scans.**\n")
        println(io, "| session | scored | human kept | threshold | predicted elite | of which kept |\n|---|---:|---:|---:|---:|---:|")
        for b in by_session
            println(io, "| $(b["session"]) | $(b["scored"]) | $(b["human_kept"]) | $(b["threshold"]) | $(b["predicted_elite"]) | $(b["of_which_human_kept"]) |")
        end
    end
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && elite_main()

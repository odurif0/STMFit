#!/usr/bin/env julia
# Joint count of the repeated scans of one molecule from the per-scan GCV curves
# of a run of run_molecule_consensus_chitosan.jl (campaign J, stage 1; see
# docs/src/journal.md, 2026-10-26, for the rule and the comparisons, fixed
# before this script was run). Label-free: score tables, counting summary,
# molecule tracks and SXM acquisition times; no benchmark label is read.
#
# Usage:
#   julia --project=. scripts/pipeline/select_joint_count.jl --run RUN_DIR \
#       --data-dir RAW_DIR --outdir NEW_DIR [--scores SCORES_DIR]
#
# --scores   folder holding <scan>/ell_scores.tsv and circ_scores.tsv at any depth
#            (the counting_chunk*/ folders of the run; default RUN_DIR).
# --data-dir raw scans of the run, for the acquisition order of the split halves.
#
# Outputs:
#   joint_counts.tsv  per scan: track, N_selected (per scan), N_eff (Level 1),
#                     N_final (current consensus), N_vote_eff (vote of N_eff),
#                     N_joint, J_margin (J of the runner-up minus J of N_joint)
#   split_half.tsv    per track of at least 4 scans and per rule: the count of
#                     the odd and of the even scans in time order
#   report.md         counts changed and split-half agreement by rule
using STMFit.JointCount, STMFit.ReconstructedUnitAssignment, STMFit.MoleculeConsensus
using STMSXMIO: read_sxm
using Printf, Statistics

function joint_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--run", "--data-dir", "--outdir", "--scores") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--run", "--data-dir", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
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

const RULES = ("selected_vote", "eff_vote", "joint")

"Count of a group of scans by one rule."
function group_count(rule, fs, n_selected, n_eff, curves)
    rule == "selected_vote" && return vote_count([n_selected[f] for f in fs])
    rule == "eff_vote" && return vote_count([n_eff[f] for f in fs])
    n, _ = joint_count([get(curves, f, Dict{Int,Float64}()) for f in fs])
    return n == 0 ? vote_count([n_eff[f] for f in fs]) : n
end

function joint_main(args=ARGS)
    o = joint_options(args)
    o === nothing && return println("select_joint_count.jl --run RUN_DIR --data-dir RAW_DIR --outdir NEW_DIR [--scores SCORES_DIR]")
    run = o["--run"]
    curves = read_score_curves(get(o, "--scores", run))
    _, summary = read_table(joinpath(run, "counting_summary.tsv"))
    n_selected = Dict(basename(r["filepath"]) => parse(Int, r["N_selected"]) for r in summary)
    n_eff = Dict(basename(r["filepath"]) => parse(Int, r["N_eff"]) for r in summary)
    _, cons = read_table(joinpath(run, "consensus", "consensus.tsv"))
    track = Dict(r["file"] => r["track"] for r in cons)
    n_final = Dict(r["file"] => parse(Int, r["N_final"]) for r in cons)
    members = Dict{String,Vector{String}}()
    for r in cons
        push!(get!(members, r["track"], String[]), r["file"])
    end
    raw = raw_index(o["--data-dir"])
    acquired(f) = haskey(raw, f) ? scan_frame(read_sxm(raw[f])).acquired : error("Raw scan missing: $f")
    for (t, fs) in members
        length(fs) >= 4 || continue
        times = Dict(f => acquired(f) for f in fs)
        sort!(fs; by=f -> (times[f], f))
    end
    rows = Dict{String,Any}[]
    for (t, fs) in sort(collect(members); by=first)
        big = length(fs) >= 3
        nj, J = big ? joint_count([get(curves, f, Dict{Int,Float64}()) for f in fs]) : (0, Dict{Int,Float64}())
        margin = length(J) > 1 ? (s = sort(collect(values(J))); s[2] - s[1]) : NaN
        nv = big ? vote_count([n_eff[f] for f in fs]) : 0
        for f in fs
            push!(rows, Dict("file" => f, "track" => t, "track_size" => length(fs), "N_selected" => n_selected[f],
                             "N_eff" => n_eff[f], "N_final" => n_final[f],
                             "N_vote_eff" => big ? nv : n_selected[f],
                             "N_joint" => big && nj > 0 ? nj : (big ? nv : n_selected[f]),
                             "J_margin" => isnan(margin) ? "NA" : @sprintf("%.4f", margin)))
        end
    end
    mkpath(o["--outdir"])
    cols = ["file", "track", "track_size", "N_selected", "N_eff", "N_final", "N_vote_eff", "N_joint", "J_margin"]
    open(joinpath(o["--outdir"], "joint_counts.tsv"), "w") do io
        println(io, join(cols, '\t'))
        foreach(r -> println(io, join((r[c] for c in cols), '\t')), sort(rows; by=r -> r["file"]))
    end
    halves = Tuple{String,Int,String,Int,Int}[]
    for (t, fs) in sort(collect(members); by=first)
        length(fs) >= 4 || continue
        a, b = split_halves(fs)
        for rule in RULES
            push!(halves, (t, length(fs), rule, group_count(rule, a, n_selected, n_eff, curves),
                           group_count(rule, b, n_selected, n_eff, curves)))
        end
    end
    open(joinpath(o["--outdir"], "split_half.tsv"), "w") do io
        println(io, join(["track", "scans", "rule", "N_odd", "N_even"], '\t'))
        foreach(h -> println(io, join(h, '\t')), halves)
    end
    md = IOBuffer()
    println(md, "# Joint count of repeated scans (campaign J, stage 1)\n")
    println(md, "Run `$(run)`: $(length(rows)) scans, $(count(fs -> length(fs) >= 3, values(members))) tracks of at least 3 scans; $(length(curves)) score curves, $(count(isempty, values(curves))) without a valid candidate.\n")
    println(md, "| count | scans changed against N_final |\n|---|---:|")
    for c in ("N_selected", "N_vote_eff", "N_joint")
        println(md, "| $c | $(count(r -> r[c] != r["N_final"], rows)) |")
    end
    println(md, "\n## Split-half agreement (tracks of at least 4 scans, odd against even scans in time order)\n")
    println(md, "| rule | tracks | halves agree | mean abs difference | difference > 2 |\n|---|---:|---:|---:|---:|")
    for rule in RULES
        d = [abs(h[4] - h[5]) for h in halves if h[3] == rule]
        isempty(d) && continue
        println(md, @sprintf("| %s | %d | %d (%.0f %%) | %.2f | %d |", rule, length(d), count(==(0), d), 100 * count(==(0), d) / length(d), mean(d), count(>(2), d)))
    end
    println(md, "\nRules: `selected_vote` = strict majority of the per-scan N_selected, else the lower median (the current consensus rule on a half); `eff_vote` = the same on N_eff; `joint` = argmin of J.")
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && joint_main()

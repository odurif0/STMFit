#!/usr/bin/env julia
# Sequences of an application run (10–20mer) with the robustness of every unit.
# Label-free: reads the per-lobe calls table written by
# scripts/grading/report_training_stability.jl (application_calls.tsv) and, if
# given, a table of scans. No benchmark label is read.
#
# Usage:
#   julia --project=. scripts/data/list_application_sequences.jl --calls application_calls.tsv \
#       --outdir NEW_DIR [--model full] [--scans TABLE.tsv] [--threshold 0.8]
#
# --model       column of the calls table to show (default `full`: the frozen
#               assignment learned on the full 6-mer set).
# --scans       tab-separated table whose first column is the scan file name (a
#               quality table works); only those scans are listed, in table order.
# --threshold   a unit is robust when at least this share of the models calls it
#               the same way (default 0.8, descriptive, not tuned).
# Output: sequences.tsv and report.md. Calls: 0 = class 0 (dimmer), 1 = class 1
# (brighter), ? = unclassified. Robustness: . dim, # bright, ~ unstable, ? unknown.
# A robust call is a call the method repeats across training subsets; it is not
# a verified K or N (docs/src/journal.md, campaigns R and D, 2026-10-25).
using STMFit.UnitReliability, STMFit.ReconstructedUnitAssignment
using Printf

function sequences_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}("--model" => "full", "--threshold" => "0.8"); i = 1
    while i <= length(args)
        args[i] in ("--calls", "--outdir", "--model", "--scans", "--threshold") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    for k in ("--calls", "--outdir")
        haskey(o, k) || error("Required: $k")
    end
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

const ROBUST_MARK = Dict(:dim => '.', :bright => '#', :unstable => '~', :unknown => '?')

"Calls and robustness strings of one scan, units in chain order."
function scan_sequence(rows, model, threshold)
    rows = sort(rows; by=r -> parse(Int, r["lobe"]))
    calls = join((get(r, model, "") in ("0", "1") ? r[model] : "?" for r in rows))
    marks = join((ROBUST_MARK[call_robustness(parse(Float64, replace(r["bright_share"], "NA" => "NaN")); threshold)] for r in rows))
    return calls, marks
end

function sequences_main(args=ARGS)
    o = sequences_options(args)
    o === nothing && return println("list_application_sequences.jl --calls application_calls.tsv --outdir NEW_DIR [--model full] [--scans TABLE.tsv] [--threshold 0.8]")
    threshold, model = parse(Float64, o["--threshold"]), o["--model"]
    header, rows = read_table(o["--calls"])
    model in header || error("No column $model in $(o["--calls"])")
    by_scan = Dict{String,Vector{Dict{String,String}}}()
    for r in rows
        push!(get!(by_scan, r["file"], Dict{String,String}[]), r)
    end
    scans = if haskey(o, "--scans")
        h, t = read_table(o["--scans"])
        [r[first(h)] for r in t]
    else
        sort(collect(keys(by_scan)))
    end
    mkpath(o["--outdir"])
    tsv = open(joinpath(o["--outdir"], "sequences.tsv"), "w")
    println(tsv, join(["file", "units", "calls", "robustness", "bright", "robust_bright", "unstable"], '\t'))
    md = IOBuffer()
    println(md, "# Sequences of the application run, with the robustness of every unit\n")
    println(md, "Calls of the model `$model` (0 = dimmer class, 1 = brighter class, ? = unclassified); robustness over the models of `$(basename(o["--calls"]))`: `.` dim, `#` bright, `~` unstable, `?` unknown, at a threshold of $(threshold) of the models. A robust call is a call the method repeats across training subsets, not a verified K or N.\n")
    println(md, "| scan | units | calls | robustness |\n|---|---:|---|---|")
    tally = Dict(:dim => 0, :bright => 0, :unstable => 0, :unknown => 0)
    for f in scans
        if !haskey(by_scan, f)
            println(md, "| $f | | absent from the calls table | |")
            continue
        end
        calls, marks = scan_sequence(by_scan[f], model, threshold)
        for c in marks
            tally[first(k for (k, v) in ROBUST_MARK if v == c)] += 1
        end
        println(md, "| $f | $(length(calls)) | `$calls` | `$marks` |")
        println(tsv, join([f, length(calls), calls, marks, count(==('1'), calls), count(==('#'), marks), count(==('~'), marks)], '\t'))
    end
    close(tsv)
    n = sum(values(tally))
    n > 0 && println(md, "\nUnits listed: $n; dim in at least $(round(Int, 100threshold))% of the models: $(tally[:dim]), bright: $(tally[:bright]), unstable: $(tally[:unstable]), unknown: $(tally[:unknown]).")
    write(joinpath(o["--outdir"], "report.md"), take!(md))
    print(read(joinpath(o["--outdir"], "report.md"), String))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && sequences_main()

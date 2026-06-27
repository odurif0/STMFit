#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# grade_unit_assignment.jl — Grade GlcNAc/GlcN (0/1) unit assignment per lobe
#
# This script evaluates predicted unit assignments against the ground-truth
# sequence. It is INTENTIONALLY OUTSIDE the fitting/model-selection path:
# it may use the true sequence to grade candidate methods, but those labels
# must not be used by the fitter itself (same label-free rule as
# grade_chitosan_benchmark.jl).
#
# Two conventions are reported:
#   1. Physical (label-free): GlcNAc = highest-amplitude cluster. The 0↔1
#      flip is fixed by physics (amplitude), not by the truth. Only the
#      sequence orientation (identity vs reverse) is tested (2 alignments).
#   2. Oracle (supervised upper bound): best of 4 alignments (identity,
#      reverse, flip, reverse+flip). Uses the truth to choose the best
#      alignment — clearly labeled as supervised.
#
# The 4 alignments handle the orientation ambiguity of the PCA axis:
#   identity:     pred[i]          vs truth[i]
#   reverse:      pred[N-i+1]      vs truth[i]    (chain read backwards)
#   flip:         1 - pred[i]      vs truth[i]    (0↔1 swapped)
#   reverse+flip: 1 - pred[N-i+1]  vs truth[i]    (both)
#
# Examples:
#   julia --project=. test/grade_unit_assignment.jl \
#       --predictions results/unit_assignment/assigned_sequences.tsv \
#       --truth benchmarks/chitosan_240817_unit_sequences.tsv \
#       --out results/benchmark_grades/unit_assignment.tsv
# ──────────────────────────────────────────────────────────────────────────────

using Printf
using DelimitedFiles

const DEFAULT_TRUTH = "benchmarks/chitosan_240817_unit_sequences.tsv"
const DEFAULT_PRED = "results/unit_assignment/assigned_sequences.tsv"
const DEFAULT_OUT = "results/benchmark_grades/unit_assignment.tsv"

struct GradeOptions
    predictions::String
    truth::String
    out_tsv::String
    file_column::String
    lobe_column::String
    predicted_column::String
    physical_column::String
    amplitude_column::String
    primary_only::Bool
end

function _parse_cli(args)
    predictions = DEFAULT_PRED
    truth = DEFAULT_TRUTH
    out_tsv = DEFAULT_OUT
    file_column = "file"
    lobe_column = "lobe"
    predicted_column = "predicted"
    physical_column = "physical_label"
    amplitude_column = "amplitude"
    primary_only = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--predictions"
            i < length(args) || error("--predictions requires a path")
            predictions = args[i+1]; i += 2
        elseif startswith(arg, "--predictions=")
            predictions = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--truth"
            i < length(args) || error("--truth requires a path")
            truth = args[i+1]; i += 2
        elseif startswith(arg, "--truth=")
            truth = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a path")
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--file-column"
            i < length(args) || error("--file-column requires a name")
            file_column = args[i+1]; i += 2
        elseif startswith(arg, "--file-column=")
            file_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--lobe-column"
            i < length(args) || error("--lobe-column requires a name")
            lobe_column = args[i+1]; i += 2
        elseif startswith(arg, "--lobe-column=")
            lobe_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--predicted-column"
            i < length(args) || error("--predicted-column requires a name")
            predicted_column = args[i+1]; i += 2
        elseif startswith(arg, "--predicted-column=")
            predicted_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--physical-column"
            i < length(args) || error("--physical-column requires a name")
            physical_column = args[i+1]; i += 2
        elseif startswith(arg, "--physical-column=")
            physical_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--amplitude-column"
            i < length(args) || error("--amplitude-column requires a name")
            amplitude_column = args[i+1]; i += 2
        elseif startswith(arg, "--amplitude-column=")
            amplitude_column = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--primary-only"
            primary_only = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/grade_unit_assignment.jl [options]

            Options:
              --predictions PATH       Predicted assignment TSV [$(DEFAULT_PRED)]
              --truth PATH             Ground-truth sequence TSV [$(DEFAULT_TRUTH)]
              --out PATH               Per-file grade TSV [$(DEFAULT_OUT)]
              --file-column NAME       File column in predictions [file]
              --lobe-column NAME       Lobe index column [lobe]
              --predicted-column NAME  Raw predicted label column [predicted]
              --physical-column NAME   Physical-mapped label column [physical_label]
              --amplitude-column NAME  Amplitude column for physical convention [amplitude]
              --primary-only           Only grade files with quality=clean or clean_target

            Predictions TSV: one row per lobe, with columns:
              $(file_column), $(lobe_column), $(predicted_column), [$(physical_column)], [$(amplitude_column)]

            Truth TSV: one row per file, with columns:
              file, sequence, [quality], [target_N], [notes]
              sequence = ordered string of 0/1 (t_nm increasing, PCA convention)
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    return GradeOptions(predictions, truth, out_tsv, file_column, lobe_column,
                        predicted_column, physical_column, amplitude_column, primary_only)
end

# Read TSV, skipping comment lines (starting with #) and empty lines.
function _read_tsv_skip_comments(path::AbstractString)
    lines = readlines(path)
    data_lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '#'), lines)
    isempty(data_lines) && error("No data lines in TSV: $path")
    header = split(data_lines[1], '\t'; keepempty=true)
    rows = Vector{Dict{String,String}}()
    for line in data_lines[2:end]
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (j, name) in enumerate(header)
            row[name] = j <= length(vals) ? vals[j] : ""
        end
        push!(rows, row)
    end
    return header, rows
end

_basename_file(s::AbstractString) = basename(strip(s))

# Parse a sequence string "010110" into Vector{Int} [0,1,0,1,1,0].
function _parse_seq(s::AbstractString)
    t = strip(s)
    isempty(t) && return Int[]
    return [parse(Int, c) for c in t if c in "01"]
end

# Compute per-position accuracy between two equal-length Int vectors.
function _seq_accuracy(pred::Vector{Int}, truth::Vector{Int})
    n = min(length(pred), length(truth))
    n == 0 && return 0.0, 0
    correct = count(pred[i] == truth[i] for i in 1:n)
    return correct / n, correct
end

# Levenshtein edit distance between two Int vectors.
function _edit_distance(a::Vector{Int}, b::Vector{Int})
    la, lb = length(a), length(b)
    la == 0 && return lb
    lb == 0 && return la
    prev = collect(0:lb)
    curr = zeros(Int, lb + 1)
    for i in 1:la
        curr[1] = i
        for j in 1:lb
            cost = a[i] == b[j] ? 0 : 1
            curr[j+1] = min(prev[j+1] + 1, curr[j] + 1, prev[j] + cost)
        end
        prev, curr = curr, prev
    end
    return prev[lb+1]
end

# Apply alignment transform to a predicted sequence.
# alignment ∈ {"identity", "reverse", "flip", "reverse+flip"}
function _align(pred::Vector{Int}, alignment::String)
    if alignment == "identity"
        return copy(pred)
    elseif alignment == "reverse"
        return reverse(pred)
    elseif alignment == "flip"
        return [1 - x for x in pred]
    elseif alignment == "reverse+flip"
        return [1 - x for x in reverse(pred)]
    else
        error("Unknown alignment: $alignment")
    end
end

const ALL_ALIGNMENTS = ["identity", "reverse", "flip", "reverse+flip"]
const PHYSICAL_ALIGNMENTS = ["identity", "reverse"]

# Compute grade for one file under a given alignment.
function _grade_one(pred_seq::Vector{Int}, truth_seq::Vector{Int}, alignment::String)
    aligned = _align(pred_seq, alignment)
    acc, correct = _seq_accuracy(aligned, truth_seq)
    ed = _edit_distance(aligned, truth_seq)
    n = min(length(aligned), length(truth_seq))
    # Confusion matrix: [truth\pred] = [[TP(1,1), FN(1,0)], [FP(0,1), TN(0,0)]]
    tp = fn = fp = tn = 0
    for i in 1:n
        t, p = truth_seq[i], aligned[i]
        if t == 1 && p == 1; tp += 1
        elseif t == 1 && p == 0; fn += 1
        elseif t == 0 && p == 1; fp += 1
        else tn += 1
        end
    end
    seq_match = (length(aligned) == length(truth_seq)) && all(aligned .== truth_seq)
    return (accuracy=acc, correct=correct, n=n, edit_dist=ed,
            tp=tp, fn=fn, fp=fp, tn=tn, seq_match=seq_match)
end

# Determine physical mapping: which label (0 or 1) corresponds to the
# higher-amplitude cluster (GlcNAc). Returns a remapping function.
# This is label-free: it uses only the amplitude, not the truth.
function _physical_remap(pred_rows::Vector{Dict{String,String}}, opt::GradeOptions)
    isempty(pred_rows) && return identity
    has_amplitude = haskey(pred_rows[1], opt.amplitude_column)
    if !has_amplitude
        return identity  # no amplitude → no physical remap possible
    end
    # Mean amplitude per predicted label
    amps_by_label = Dict{Int,Vector{Float64}}()
    for row in pred_rows
        label = parse(Int, row[opt.predicted_column])
        amp = tryparse(Float64, strip(row[opt.amplitude_column]))
        amp === nothing && continue
        isfinite(amp) || continue
        push!(get!(amps_by_label, label, Float64[]), amp)
    end
    mean_amp = Dict(k => sum(v) / length(v) for (k, v) in amps_by_label if !isempty(v))
    # GlcNAc (label 1 in truth) = higher amplitude
    labels = collect(keys(mean_amp))
    length(labels) < 2 && return identity
    high_label = first(sort(labels; by=k -> mean_amp[k], rev=true))
    # If high_label == 1, no remap needed. If high_label == 0, flip.
    if high_label == 1
        return identity
    else
        return x -> 1 - x
    end
end

function main(args=ARGS)
    opt = _parse_cli(args)

    # Read truth
    _, truth_rows = _read_tsv_skip_comments(opt.truth)
    truth_by_file = Dict{String,Vector{Int}}()
    truth_quality = Dict{String,String}()
    for row in truth_rows
        file = _basename_file(row["file"])
        seq = _parse_seq(row["sequence"])
        truth_by_file[file] = seq
        truth_quality[file] = get(row, "quality", "clean")
    end

    # Read predictions (one row per lobe)
    _, pred_rows = _read_tsv_skip_comments(opt.predictions)
    pred_by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in pred_rows
        haskey(row, opt.file_column) || error("Missing file column: $(opt.file_column)")
        haskey(row, opt.predicted_column) || error("Missing predicted column: $(opt.predicted_column)")
        file = _basename_file(row[opt.file_column])
        push!(get!(pred_by_file, file, Dict{String,String}[]), row)
    end

    # Sort lobes by lobe index within each file
    for (file, rows) in pred_by_file
        haskey(rows[1], opt.lobe_column) && sort!(rows, by=r -> parse(Int, r[opt.lobe_column]))
    end
    global_remap = _physical_remap(pred_rows, opt)

    mkpath(dirname(opt.out_tsv))

    # Aggregate stats
    files = sort(collect(keys(pred_by_file)))
    physical_total = physical_correct = 0
    oracle_total = oracle_correct = 0
    physical_seq_match = 0
    oracle_seq_match = 0
    n_files_graded = 0

    open(opt.out_tsv, "w") do io
        println(io, join([
            "file", "N_pred", "N_truth", "quality",
            "phys_alignment", "phys_accuracy", "phys_correct", "phys_edit_dist", "phys_seq_match",
            "phys_tp", "phys_fn", "phys_fp", "phys_tn",
            "oracle_alignment", "oracle_accuracy", "oracle_correct", "oracle_edit_dist", "oracle_seq_match",
            "oracle_tp", "oracle_fn", "oracle_fp", "oracle_tn",
        ], '\t'))

        for file in files
            truth_seq = get(truth_by_file, file, Int[])
            isempty(truth_seq) && continue  # no ground truth for this file

            quality = get(truth_quality, file, "clean")
            if opt.primary_only && quality == "poor_quality"
                continue
            end

            rows = pred_by_file[file]
            pred_raw = [parse(Int, r[opt.predicted_column]) for r in rows]

            # Physical convention: one global amplitude-based remap, then test identity/reverse.
            pred_phys = [global_remap(p) for p in pred_raw]

            best_phys = nothing
            for al in PHYSICAL_ALIGNMENTS
                g = _grade_one(pred_phys, truth_seq, al)
                if best_phys === nothing || g.accuracy > best_phys.accuracy
                    best_phys = (alignment=al, g...)
                end
            end

            # Oracle convention: best of all 4 alignments on raw predictions
            best_oracle = nothing
            for al in ALL_ALIGNMENTS
                g = _grade_one(pred_raw, truth_seq, al)
                if best_oracle === nothing || g.accuracy > best_oracle.accuracy
                    best_oracle = (alignment=al, g...)
                end
            end

            n_phys = best_phys.n
            n_ora = best_oracle.n
            physical_total += n_phys
            physical_correct += best_phys.correct
            oracle_total += n_ora
            oracle_correct += best_oracle.correct
            physical_seq_match += best_phys.seq_match ? 1 : 0
            oracle_seq_match += best_oracle.seq_match ? 1 : 0
            n_files_graded += 1

            println(io, join([
                file,
                string(length(pred_raw)),
                string(length(truth_seq)),
                quality,
                best_phys.alignment,
                @sprintf("%.4f", best_phys.accuracy),
                string(best_phys.correct),
                string(best_phys.edit_dist),
                string(best_phys.seq_match),
                string(best_phys.tp), string(best_phys.fn),
                string(best_phys.fp), string(best_phys.tn),
                best_oracle.alignment,
                @sprintf("%.4f", best_oracle.accuracy),
                string(best_oracle.correct),
                string(best_oracle.edit_dist),
                string(best_oracle.seq_match),
                string(best_oracle.tp), string(best_oracle.fn),
                string(best_oracle.fp), string(best_oracle.tn),
            ], '\t'))
        end
    end

    # Summary
    println("\nUnit assignment grade")
    println("  predictions:    ", opt.predictions)
    println("  truth:          ", opt.truth)
    println("  output:         ", opt.out_tsv)
    println("  files graded:   ", n_files_graded)
    println()
    println("  Physical convention (label-free, GlcNAc=amp max):")
    println("    per-blob accuracy: ", physical_correct, "/", physical_total,
            " (", physical_total > 0 ? @sprintf("%.1f%%", 100*physical_correct/physical_total) : "NA", ")")
    println("    sequence match:    ", physical_seq_match, "/", n_files_graded)
    println()
    println("  Oracle convention (supervised upper bound, best of 4 alignments):")
    println("    per-blob accuracy: ", oracle_correct, "/", oracle_total,
            " (", oracle_total > 0 ? @sprintf("%.1f%%", 100*oracle_correct/oracle_total) : "NA", ")")
    println("    sequence match:    ", oracle_seq_match, "/", n_files_graded)
    println()
    phys_acc = physical_total > 0 ? 100*physical_correct/physical_total : NaN
    ora_acc = oracle_total > 0 ? 100*oracle_correct/oracle_total : NaN
    gap = ora_acc - phys_acc
    println("  Gap (oracle - physical): ", @sprintf("%.1f%%", gap),
            gap < 5.0 ? "  [physical convention validated]" :
            gap < 15.0 ? "  [moderate gap — some supervision gain available]" :
            "  [large gap — physical mapping imperfect]")
end

main()

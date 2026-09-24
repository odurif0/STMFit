#!/usr/bin/env julia
# Two complete, independent executions of the frozen counting + assignment path.
# Historical counts/predictions and external grading are deliberately not inputs.
module HybridReproduction
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
using SHA

const REQUIRED = Set(["--data-dir", "--count-config", "--config", "--templates", "--outdir"])

function options(args)
    args == ["--help"] && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        haskey(o, k) && error("Repeated option $k")
        if k == "--dry-run"; o[k] = "true"; i += 1; continue; end
        k in REQUIRED && i < length(args) && !startswith(args[i+1], "--") || error("Missing/forbidden option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in REQUIRED) || error("All reproduction inputs required")
    isdir(o["--data-dir"]) || error("Missing raw directory")
    all(isfile(o[k]) for k in ("--count-config", "--config", "--templates")) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function merge_counts(paths, files, output)
    header = String[]; rows = Dict{String,String}[]; seen = Set{String}()
    for path in paths
        h, rs = read_table(path)
        isempty(header) && (header = h)
        h == header || error("Counting shard schemas differ")
        counts = selected_counts(path)
        isempty(intersect(seen, keys(counts))) || error("Duplicate scan across counting shards")
        union!(seen, keys(counts))
        all(get(r, "selection_policy", "") == "support_midpoint_hybrid" for r in rs) || error("Wrong counting policy")
        any(startswith(get(r, k, ""), "adaptive_support") for r in rs for k in
            ("refined_policy", "refined_source", "selection_source")) && error("Adaptive support is outside this reproduction")
        append!(rows, rs)
    end
    seen == Set(files) || error("Counting omitted or added scans")
    sort!(rows; by=r -> basename(r["filepath"]))
    write_table(output, header, rows)
    return selected_counts(output)
end

function execute(o; runner=run_stage, thread_budget=min(Threads.nthreads(), 4))
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
    thread_budget > 0 || error("Positive thread budget required")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    cfg = TOML.parsefile(o["--count-config"])
    cfg["model"]["selection_criterion"] == "gcv" || error("GCV required")
    cfg["model"]["selection_policy"] == "support_midpoint_hybrid" || error("Frozen hybrid policy required")
    assignment = load_config(o["--config"])
    assignment["preprocessing"]["patch_residual_filter"] == "smooth_residual" || error("Symmetric residuals required")
    rawfiles = raw_index(o["--data-dir"]); files = sort(collect(keys(rawfiles)))
    nchunks = min(length(files), thread_budget, 4)
    tag = @sprintf("%03d", round(Int, 100cfg["model"]["max_overlap"]))
    println("Two fresh repetitions; $(length(files)) raw scans; $nchunks counting shards each")
    haskey(o, "--dry-run") && return println("No pixels, output, count/geometry caches, labels or grading are used by this dry-run")
    out = abspath(o["--outdir"]); mkpath(joinpath(out, "logs")); stage = "inputs"
    try
        raw = prepare_raw(o["--data-dir"], out, Set(files))
        write_table(joinpath(out, "input_hashes.tsv"), ["input", "sha256"],
            [Dict("input"=>k, "sha256"=>bytes2hex(sha256(read(o[k])))) for k in ("--count-config", "--config", "--templates")])
        write_table(joinpath(out, "raw_hashes.tsv"), ["file", "sha256"],
            [Dict("file"=>f, "sha256"=>bytes2hex(sha256(read(rawfiles[f])))) for f in files])
        for repeat in 1:2
            rep = joinpath(out, "repeat$repeat"); mkpath(joinpath(rep, "logs"))
            stage = "repeat$(repeat)_counting"
            paths = String[]
            for i in 1:nchunks
                name = nchunks == 1 ? "summary_overlap$(tag)_hard.tsv" :
                    @sprintf("summary_overlap%s_hard_chunk%02dof%02d.tsv", tag, i, nchunks)
                push!(paths, joinpath(rep, "counting_chunk$i", name))
            end
            @sync for i in 1:nchunks
                Threads.@spawn runner(rep, "counting_chunk$i", "batch_full.jl",
                    [string(length(files)), "--config", abspath(o["--count-config"]), "--data-dir", raw,
                     "--outdir", dirname(paths[i]), "--tsv", joinpath(rep, "no_triage_input.tsv"),
                     "--skip-1d", "--chunk", "$i/$nchunks"]; threads=1)
            end
            selected = joinpath(rep, "counting_summary.tsv")
            counts = merge_counts(paths, files, selected)
            stage = "repeat$(repeat)_assignment"
            dest = joinpath(rep, "assignment")
            runner(rep, "assignment", "run_reconstructed_chitosan.jl",
                ["--data-dir", raw, "--count-config", abspath(o["--count-config"]),
                 "--config", abspath(o["--config"]), "--templates", abspath(o["--templates"]),
                 "--selected-summary", selected, "--outdir", dest]; threads=min(thread_budget, 4))
            for name in ("features.tsv", "features_split.tsv", "predictions.tsv")
                check_counts(joinpath(dest, name), counts)
            end
        end
    catch e
        write_table(joinpath(out, "failures.tsv"), ["stage", "reason"],
            [Dict("stage"=>stage, "reason"=>replace(sprint(showerror, e), '\n'=>' ', '\t'=>' '))])
        rethrow()
    end
    println("Both raw-to-prediction repetitions complete; no best-repeat selection. External grading remains separate.")
end

function main(args=ARGS)
    o = options(args)
    o === nothing ? println("run_hybrid_reproduction.jl --data-dir DIR --count-config TOML --config TOML --templates TSV --outdir NEW_DIR [--dry-run]") : execute(o)
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && HybridReproduction.main()

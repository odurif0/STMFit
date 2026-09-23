#!/usr/bin/env julia
# Two complete, label-free arms; only the per-lobe patch sampling frame changes.
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
using SHA
const ORIENTATION_OPTIONS = Set(["--data-dir", "--count-config", "--config", "--settings",
    "--features", "--split-features", "--templates", "--outdir"])

function orientation_comparison(args=ARGS; runner=run_stage)
    if "--help" in args || "-h" in args
        println("run_local_orientation_comparison.jl: ", join(sort(collect(ORIENTATION_OPTIONS)), " VALUE "), " VALUE [--dry-run]")
        println("Global versus local-tangent patches, unchanged saved geometry/N/templates/classifier settings. No labels.")
        return
    end
    opts = Dict{String,String}(); filtered = String[]; i = 1
    while i <= length(args)
        key = args[i]
        haskey(opts, key) && error("Repeated comparison option: $key")
        if key == "--dry-run"
            opts[key] = "true"; push!(filtered, key); i += 1; continue
        end
        key in ORIENTATION_OPTIONS && i < length(args) && !startswith(args[i+1], "--") ||
            error("Missing/forbidden comparison option: $key")
        opts[key] = args[i+1]
        key == "--settings" || append!(filtered, args[i:i+1])
        i += 2
    end
    all(haskey(opts, k) for k in ORIENTATION_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered)
    settings = PatchFrames.load_frame_settings(opts["--settings"])
    load_config(opts["--config"])
    _, base = lobe_table(opts["--features"]); _, splitrows = lobe_table(opts["--split-features"])
    require_same_keys(base, splitrows, "split geometry")
    files = Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"]))) == files || error("Raw/cache cohorts differ")
    counts = Dict(f => count(k -> first(k) == f, keys(base)) for f in files)
    for table in (base, splitrows), (key, row) in table
        parse(Int, row["N"]) == counts[first(key)] || error("Saved N differs from lobe count")
    end
    # This small geometry-only calculation also validates the full cohort in dry-run.
    frames = PatchFrames.local_frame_rows(values(base), settings)
    haskey(opts, "--dry-run") && return println("Dry run: $(length(files)) files, $(length(base)) keys; two frozen-geometry arms; no SXM reading, fitting or output.")
    out = abspath(opts["--outdir"]); mkpath(joinpath(out, "logs")); stage = "inputs"
    try
        inputs = [Dict("input"=>k, "sha256"=>bytes2hex(sha256(read(opts[k]))))
            for k in sort(collect(ORIENTATION_OPTIONS)) if isfile(opts[k])]
        write_table(joinpath(out, "input_hashes.tsv"), ["input", "sha256"], inputs)
        write_table(joinpath(out, "raw_hashes.tsv"), ["file", "sha256"],
            [Dict("file"=>f, "sha256"=>bytes2hex(sha256(read(p)))) for (f,p) in sort(collect(raw_index(opts["--data-dir"])))])
        basepath = joinpath(out, "cached_features.tsv"); splitpath = joinpath(out, "cached_split.tsv")
        cp(opts["--features"], basepath); cp(opts["--split-features"], splitpath)
        framepath = PatchFrames.write_frames(joinpath(out, "patch_frames.tsv"), frames)
        PatchFrames.read_frames(framepath, values(base))
        common = ["--data-dir", abspath(opts["--data-dir"]), "--count-config", abspath(opts["--count-config"]),
            "--config", abspath(opts["--config"]), "--templates", abspath(opts["--templates"]),
            "--features", basepath, "--split-features", splitpath]
        for mode in ("reference", "local_frame")
            stage = mode
            extra = mode == "reference" ? String[] : ["--patch-frames", framepath]
            runner(out, mode, "run_reconstructed_chitosan.jl", vcat(common, extra, ["--outdir", joinpath(out, mode)]))
            cp(basepath, joinpath(out, mode, "features.tsv")); cp(splitpath, joinpath(out, mode, "features_split.tsv"))
            check_counts(joinpath(out, mode, "predictions.tsv"), counts)
        end
    catch err
        write_table(joinpath(out, "failures.tsv"), ["stage", "reason"],
            [Dict("stage"=>stage, "reason"=>replace(sprint(showerror, err), '\n'=>' ', '\t'=>' '))])
        rethrow()
    end
    println("Completed local patch-frame comparison. External benchmark grading remains separate.")
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && orientation_comparison()

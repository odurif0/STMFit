#!/usr/bin/env julia
# Bounded experiment: cached own-N GCV counts versus a fresh native GCV sweep.
# Both arms subsequently refit at fixed N; sweep geometry is NEVER reused for
# assignment. No production selector, inference setting, or benchmark is changed.
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))

const RESELECTION_OPTIONS = Set(["--data-dir", "--count-config", "--config",
    "--selected-summary", "--templates", "--outdir"])

function parse_reselection(args)
    if "--help" in args || "-h" in args
        println("""
        Usage: julia -t 4 --project=. test/run_gcv_reselection.jl [options]
        Required: --data-dir DIR --count-config TOML --config TOML
          --selected-summary TSV --templates TSV --outdir NEW_DIR
        --dry-run checks metadata only, without reading SXM pixels or writing.
        The selected summary supplies the control's label-free counts and cohort.
        Treatment repeats extract_lobe_features.jl's raw GCV selection, including
        its unchanged intelligent sweep. It does NOT use the batch hybrid guard.
        Both assignments regenerate base/split fits, patches and learned heads.
        Julia 1.13 required; real multi-file execution belongs on Viper.
        """)
        return nothing
    end
    i = 1
    while i <= length(args)
        if args[i] == "--dry-run"
            i += 1
        else
            args[i] in RESELECTION_OPTIONS || error("Unknown or forbidden option: $(args[i])")
            i += 2
        end
    end
    opts = parse_options(args)
    all(k -> haskey(opts, k), RESELECTION_OPTIONS) || error("Missing experiment input")
    return opts
end

function reselection_summary(features, raw, files, output)
    _, rows = lobe_table(features; required=["N", "gcv", "source"])
    Set(first.(collect(keys(rows)))) == Set(files) || error("GCV sweep omitted or added files")
    summary = Dict{String,String}[]
    for file in sort(collect(files))
        chain = [rows[k] for k in sort(collect(keys(rows))) if first(k) == file]
        n = length(chain)
        all(r -> tryparse(Int, r["N"]) == n, chain) || error("Inconsistent selected N: $file")
        gcv = unique(r["gcv"] for r in chain)
        source = unique(r["source"] for r in chain)
        length(gcv) == length(source) == 1 || error("Inconsistent selected fit: $file")
        score = something(tryparse(Float64, only(gcv)), NaN)
        isfinite(score) && score >= 0 || error("Invalid selected GCV: $file")
        only(source) in ("ell", "circ") || error("Invalid selected geometry: $file")
        push!(summary, Dict("filepath"=>joinpath(raw, file), "N_selected"=>string(n),
            "gcv"=>only(gcv), "source"=>only(source), "selection_source"=>"extractor_raw_gcv"))
    end
    write_table(output, ["filepath", "N_selected", "gcv", "source", "selection_source"], summary)
    return selected_counts(output)
end

function execute_reselection(opts; exporter=export_features, runner=run_stage)
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
    cfg = TOML.parsefile(opts["--count-config"])
    get(cfg["model"], "selection_criterion", "gcv") == "gcv" || error("GCV required")
    startswith(get(cfg["model"], "selection_policy", ""), "adaptive_support") &&
        error("This experiment requires nonadaptive base support")
    assignment = load_config(opts["--config"])
    counts = selected_counts(opts["--selected-summary"])
    files = Set(keys(counts))
    Set(keys(raw_index(opts["--data-dir"]))) == files || error("Raw/control cohorts differ")
    outdir = abspath(opts["--outdir"])
    (ispath(outdir) || islink(outdir)) && error("Output already exists")
    if haskey(opts, "--dry-run")
        println("Julia $VERSION; files=$(length(files)); method=$(assignment["model"]["name"])")
        println("Control: cached N -> fresh fixed-N geometry and assignment")
        println("Treatment: unchanged extractor raw-GCV sweep -> fresh fixed-N geometry and assignment")
        println("No geometry cache, benchmark, batch hybrid policy or exhaustive-sweep override. Output: $outdir")
        return
    end
    mkpath(joinpath(outdir, "logs"))
    stage = "inputs"
    try
        raw = prepare_raw(opts["--data-dir"], outdir, files)
        control = joinpath(outdir, "control_counts.tsv")
        write_table(control, ["filepath", "N_selected"],
            [Dict("filepath"=>joinpath(raw, f), "N_selected"=>string(counts[f])) for f in sort(collect(files))])
        common = ["--data-dir", raw, "--count-config", abspath(opts["--count-config"]),
            "--config", abspath(opts["--config"]), "--templates", abspath(opts["--templates"])]
        stage = "control"
        runner(outdir, stage, "run_reconstructed_chitosan.jl",
            vcat(common, ["--selected-summary", control, "--outdir", joinpath(outdir, stage)]))
        stage = "gcv_sweep"
        sweep = joinpath(outdir, "sweep_features.tsv")
        exporter(outdir, stage, "extract_lobe_features.jl",
            ["--data-dir", raw, "--config", abspath(opts["--count-config"])], sweep;
            nfiles=length(files))
        selected = joinpath(outdir, "reselected_counts.tsv")
        newcounts = reselection_summary(sweep, raw, files, selected)
        write_table(joinpath(outdir, "count_changes.tsv"), ["file", "control_N", "reselected_N", "delta"],
            [Dict("file"=>f, "control_N"=>string(counts[f]), "reselected_N"=>string(newcounts[f]),
                "delta"=>string(newcounts[f]-counts[f])) for f in sort(collect(files))])
        stage = "reselected"
        runner(outdir, stage, "run_reconstructed_chitosan.jl",
            vcat(common, ["--selected-summary", selected, "--outdir", joinpath(outdir, stage)]))
        for (arm, expected) in (("control", counts), ("reselected", newcounts))
            for name in ("features.tsv", "features_split.tsv", "predictions.tsv")
                check_counts(joinpath(outdir, arm, name), expected)
            end
        end
    catch err
        write_table(joinpath(outdir, "failures.tsv"), ["stage", "status", "reason"],
            [Dict("stage"=>stage, "status"=>"incomplete", "reason"=>replace(sprint(showerror, err), '\n'=>' ', '\t'=>' '))])
        rethrow()
    end
    println("Completed both fixed-method arms: $outdir. External grading is separate.")
end

function reselection_main(args=ARGS)
    opts = parse_reselection(args)
    opts === nothing || execute_reselection(opts)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && reselection_main()

#!/usr/bin/env julia
# One native Julia production command: raw SXM/selected-N artifacts -> 0/1/? + QC + maps.
# This is an explicit reconstructed method, not a new name for the frozen champion.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using TOML, Printf

const ROOT = dirname(@__DIR__)
const BASE4 = "amp_prominence,amp_neighbor_ratio,integrated_prominence,amp_rel"
const VALUE_OPTIONS = Set(["--data-dir", "--count-config", "--config", "--outdir",
    "--selected-summary", "--features", "--split-features", "--patches-fwd",
    "--patches-bwd", "--descriptor-patches", "--templates",
    "--cube0", "--cube1", "--frame0", "--frame1"])

function parse_options(args)
    if "--help" in args || "-h" in args
        println("""
        Usage: julia -t 4 --project=. test/run_reconstructed_chitosan.jl [options]
        Required: --data-dir DIR --count-config TOML --config TOML --outdir NEW_DIR
        Mold inputs: --cube0 PATH --cube1 PATH --frame0 PATH --frame1 PATH
          OR --templates PATH to reuse native template output.
        Optional cached label-free inputs:
          --selected-summary PATH (filepath/N_selected; no external labels)
          --features PATH (selected-N Gaussian geometry; skips counting/base refit)
          --split-features PATH --patches-fwd PATH --patches-bwd PATH
          --descriptor-patches PATH (backward residual 9x9; mold patches are 17x17)
          Matched-residual mode regenerates all patches; patch caches are rejected.
        --dry-run: check supplied input paths and print the stages without computing.

        Production only: no benchmark labels, expected count, control sequence,
        comparison file, grading, composition prior, or profile selection.
        Outputs: predictions.tsv, summary.tsv, review_queue.tsv, plots/, stage logs.
        Missing or failed files are reported in failures.tsv; coverage gaps are errors.
        Julia 1.13 is required. No Python subprocess is used.
        """)
        return nothing
    end
    opts = Dict{String,String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--dry-run"
            opts[arg] = "true"; i += 1; continue
        end
        arg in VALUE_OPTIONS || error("Unknown or forbidden option: $arg")
        i < length(args) || error("$arg needs a value")
        haskey(opts, arg) && error("Repeated option: $arg")
        opts[arg] = args[i+1]; i += 2
    end
    for key in ("--data-dir", "--count-config", "--config", "--outdir")
        haskey(opts, key) || error("Missing required option: $key")
    end
    isdir(opts["--data-dir"]) || error("Data directory not found")
    for (key, value) in opts
        key in ("--data-dir", "--outdir", "--dry-run") && continue
        isfile(value) || error("$key input does not exist: $value")
    end
    if !haskey(opts, "--templates")
        all(k -> haskey(opts, k), ("--cube0", "--cube1", "--frame0", "--frame1")) ||
            error("Provide --templates or all four cube/frame inputs")
    end
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Output directory already exists")
    return opts
end

function run_stage(outdir, name, script, args; threads=min(Threads.nthreads(), 4))
    cmd = `$(Base.julia_cmd()) --threads=$threads --project=$ROOT $(joinpath(@__DIR__, script)) $args`
    println("[$name] ", cmd); flush(stdout)
    open(joinpath(outdir, "logs", name * ".log"), "w") do io
        run(pipeline(addenv(cmd, "GKSwstype" => "100", "OPENBLAS_NUM_THREADS" => "1"); stdout=io, stderr=io))
    end
end

function selected_counts(path)
    header, rows = read_table(path)
    all(c -> c in header, ("filepath", "N_selected")) || error("Summary requires filepath and N_selected")
    counts = Dict{String,Int}()
    for row in rows
        file = basename(row["filepath"])
        isempty(file) && error("Empty selected-summary file")
        haskey(counts, file) && error("Duplicate selected-summary file: $file")
        get(row, "status", "ok") == "ok" || error("Counting did not succeed: $file")
        n = tryparse(Int, row["N_selected"])
        n !== nothing && n > 0 || error("Invalid N_selected for $file")
        counts[file] = n
    end
    return counts
end

function check_counts(features, counts)
    _, rows = lobe_table(features)
    actual = Dict{String,Int}()
    for (file, _) in keys(rows)
        actual[file] = get(actual, file, 0) + 1
    end
    if actual != counts
        missing = sort(collect(setdiff(keys(counts), keys(actual))))
        extra = sort(collect(setdiff(keys(actual), keys(counts))))
        wrong = sort(["$f: $(actual[f])/$(counts[f])" for f in intersect(keys(actual), keys(counts)) if actual[f] != counts[f]])
        error("Feature/count coverage mismatch; missing=$(join(missing, ',')); extra=$(join(extra, ',')); rows/N_selected=$(join(wrong, ','))")
    end
    return rows
end

function raw_index(data_dir)
    index = Dict{String,String}()
    for (dir, _, files) in walkdir(data_dir), name in files
        endswith(lowercase(name), ".sxm") || continue
        haskey(index, name) && error("Duplicate raw SXM basename: $name")
        index[name] = abspath(joinpath(dir, name))
    end
    isempty(index) && error("No SXM files under $data_dir")
    return index
end

function prepare_raw(data_dir, outdir, files)
    index = raw_index(data_dir)
    all(f -> haskey(index, f), files) || error("Raw SXMs missing for supplied selected-N keys")
    path = joinpath(outdir, "raw_inputs")
    mkdir(path)
    for file in sort(collect(files))
        symlink(index[file], joinpath(path, file))
    end
    return path
end

function merge_feature_chunks(paths, output)
    header = String[]
    rows = Dict{String,String}[]
    seen = Set{Tuple{String,Int}}()
    for path in paths
        chunk_header, chunk = lobe_table(path)
        isempty(header) && (header = chunk_header)
        header == chunk_header || error("Feature chunk headers differ")
        for key in sort(collect(keys(chunk)))
            key in seen && error("Duplicate lobe across feature chunks: $key")
            push!(seen, key)
            push!(rows, chunk[key])
        end
    end
    isempty(rows) && error("No exported feature rows")
    sort!(rows; by=row -> (row["file"], parse(Int, row["lobe"])))
    return write_table(output, header, rows)
end

function export_features(outdir, name, script, args, output; nfiles,
                         thread_budget=min(Threads.nthreads(), 4), runner=run_stage)
    nfiles > 0 && thread_budget > 0 || error("Invalid feature-export workload")
    nchunks = min(nfiles, thread_budget, 4)
    if nchunks == 1
        runner(outdir, name, script, vcat(args, ["--out", output]); threads=thread_budget)
        return output
    end
    paths = [joinpath(outdir, "$(name)_chunk$(i).tsv") for i in 1:nchunks]
    # Each child gets one Julia/BLAS thread; total CPU use stays within the
    # parent's budget. Counting and downstream cohort-wide fitting stay intact.
    @sync for (i, path) in enumerate(paths)
        Threads.@spawn runner(outdir, "$(name)_chunk$(i)", script,
            vcat(args, ["--chunk", "$i/$nchunks", "--out", path]); threads=1)
    end
    return merge_feature_chunks(paths, output)
end

function cached_or_run(opts, key, output, name, script, args, outdir; fit_files=0)
    if haskey(opts, key)
        println("[$name] reuse ", opts[key])
        return abspath(opts[key])
    end
    if fit_files > 0
        export_features(outdir, name, script, args, output; nfiles=fit_files)
    else
        run_stage(outdir, name, script, vcat(args, ["--out", output]))
    end
    return output
end

function join_predictor_features(base_path, split_path, fwd_score, bwd_score, fisher_path, out)
    header, base = lobe_table(base_path)
    _, splitrows = lobe_table(split_path; required=["skew_ratio"])
    _, fwd = lobe_table(fwd_score; required=["cost_margin"])
    _, bwd = lobe_table(bwd_score; required=["cost_margin"])
    _, fisher = lobe_table(fisher_path; required=["score"])
    for (name, table) in (("split", splitrows), ("forward score", fwd), ("backward score", bwd), ("Fisher", fisher))
        require_same_keys(base, table, name)
    end
    extras = ["split_log_skew", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher"]
    any(c -> c in header, extras) && error("Recomputed predictor columns already present")
    rows = Dict{String,String}[]
    for key in sort(collect(keys(base)))
        row = copy(base[key])
        skew = something(tryparse(Float64, splitrows[key]["skew_ratio"]), NaN)
        row["split_log_skew"] = isfinite(skew) && skew > 0 ? @sprintf("%.8g", log(skew)) : "NaN"
        for (col, value) in (("mold_cc_fwd", fwd[key]["cost_margin"]), ("mold_cc_bwd", bwd[key]["cost_margin"]))
            p = something(tryparse(Float64, value), NaN)
            row[col] = isfinite(p) ? @sprintf("%.8g", p) : "NaN"
        end
        f = something(tryparse(Float64, fisher[key]["score"]), NaN)
        row["emp_fisher"] = isfinite(f) ? @sprintf("%.8g", -f) : "NaN"
        # Native predictors must see unavailable columns as nonfinite, even if
        # every lobe is unavailable (not as an absent feature definition).
        row["patch_u_asym_reconstructed"] == "NA" && (row["patch_u_asym_reconstructed"] = "NaN")
        push!(rows, row)
    end
    write_table(out, vcat(header, extras), rows)
end

function write_summary(predictions, out)
    _, rows = lobe_table(predictions)
    summary = Dict{String,String}[]
    for file in sort(unique(first.(collect(keys(rows)))))
        chain = [rows[k] for k in sort(collect(keys(rows))) if first(k) == file]
        labels = [r["predicted"] for r in chain]
        push!(summary, Dict("file" => file, "N_selected" => string(length(chain)),
            "assignment" => join(labels), "predicted_0" => string(count(==("0"), labels)),
            "predicted_1" => string(count(==("1"), labels)), "uncertain" => string(count(==("?"), labels))))
    end
    write_table(out, ["file", "N_selected", "assignment", "predicted_0", "predicted_1", "uncertain"], summary)
end

function execute_pipeline(opts)
    VERSION.major == 1 && VERSION.minor == 13 || error("This reconstruction requires Julia 1.13")
    cfg = load_config(opts["--config"])
    if cfg["preprocessing"]["patch_residual_filter"] == "smooth_residual"
        any(haskey(opts, k) for k in ("--patches-fwd", "--patches-bwd", "--descriptor-patches")) &&
            error("Matched-residual mode requires fresh patches; remove cached patch inputs")
    end
    outdir = abspath(opts["--outdir"])
    selected = get(opts, "--selected-summary", "")
    counts = isempty(selected) ? Dict{String,Int}() : selected_counts(selected)
    if haskey(opts, "--features")
        _, supplied = lobe_table(opts["--features"])
        if isempty(counts)
            for (file, _) in keys(supplied)
                counts[file] = get(counts, file, 0) + 1
            end
        end
        check_counts(opts["--features"], counts)
    end
    files = isempty(counts) ? Set(keys(raw_index(opts["--data-dir"]))) : Set(keys(counts))
    if haskey(opts, "--dry-run")
        println("Julia ", VERSION, "; method=", cfg["model"]["name"], "; files=", length(files))
        println("Stages: selected-N counting/reuse -> base/split geometry -> 9x9 descriptor + 17x17 patches -> native CC/Fisher -> native predictors -> vote -> validation/QC/maps")
        println("No computation or output writing in dry-run. Output: ", outdir)
        return
    end
    mkpath(joinpath(outdir, "logs"))
    stage = "input"
    try
        raw = prepare_raw(opts["--data-dir"], outdir, files)
        if isempty(counts)
            stage = "counting"
            countdir = joinpath(outdir, "counting")
            # --tsv is a triage INPUT. A guaranteed nonexistent path disables
            # the driver's legacy benchmark-shaped default input.
            unused_triage = joinpath(outdir, "no_triage_input.tsv")
            run_stage(outdir, stage, "batch_full.jl", [string(length(files)), "--config", abspath(opts["--count-config"]),
                "--data-dir", raw, "--outdir", countdir, "--tsv", unused_triage, "--skip-1d"])
            selected = joinpath(countdir, "summary_overlap060_hard.tsv")
            counts = selected_counts(selected)
            Set(keys(counts)) == files || error("Counting omitted input files")
        elseif isempty(selected)
            selected = joinpath(outdir, "selected_from_features.tsv")
            write_table(selected, ["filepath", "N_selected"],
                [Dict("filepath" => joinpath(raw, f), "N_selected" => string(counts[f])) for f in sort(collect(files))])
        end
        stage = "base_features"
        geometry = cached_or_run(opts, "--features", joinpath(outdir, "features.tsv"), stage,
            "extract_lobe_features.jl", ["--config", abspath(opts["--count-config"]), "--data-dir", raw,
            "--selected-summary", abspath(selected)], outdir; fit_files=length(files))
        basekeys = check_counts(geometry, counts)
        stage = "split_features"
        split_cfg = TOML.parsefile(opts["--count-config"])
        split_cfg["model"]["peak_profile"] = "split"
        split_cfg["model"]["skew_ratio_max"] = cfg["model"]["split_skew_ratio_max"]
        split_config_path = joinpath(outdir, "split_config.toml")
        open(io -> TOML.print(io, split_cfg), split_config_path, "w")
        split_features = cached_or_run(opts, "--split-features", joinpath(outdir, "features_split.tsv"), stage,
            "extract_lobe_features.jl", ["--config", split_config_path, "--data-dir", raw,
            "--selected-summary", abspath(selected)], outdir; fit_files=length(files))
        check_counts(split_features, counts)
        stage = "local_features"
        local_features = joinpath(outdir, "features_local.tsv")
        run_stage(outdir, stage, "augment_lobe_local_features.jl", ["--features", geometry, "--out", local_features])
        check_counts(local_features, counts)
        stage = "patches"
        paths = Dict{String,String}()
        for (key, name, script, half, step, prefix, side) in (
            ("--patches-fwd", "patches_fwd17", "extract_lobe_patches.jl", cfg["model"]["mold_half_nm"], cfg["model"]["mold_step_nm"], "res_p", 17),
            ("--patches-bwd", "patches_bwd17", "extract_lobe_patches_bwd.jl", cfg["model"]["mold_half_nm"], cfg["model"]["mold_step_nm"], "bwd_res_p", 17),
            ("--descriptor-patches", "patches_bwd9", "extract_lobe_patches_bwd.jl", cfg["model"]["descriptor_half_nm"], cfg["model"]["descriptor_step_nm"], "bwd_res_p", 9))
            path = cached_or_run(opts, key, joinpath(outdir, name * ".tsv"), name, script,
                ["--features", geometry, "--data-dir", raw, "--config", abspath(opts["--count-config"]),
                 "--assignment-config", abspath(opts["--config"]),
                 "--half-nm", string(half), "--step-nm", string(step)], outdir)
            header, patchkeys = lobe_table(path; required=[prefix * lpad(string(i), 3, '0') for i in 1:side^2])
            length(filter(c -> startswith(c, prefix), header)) == side^2 || error("Unexpected patch dimensions: $path")
            require_same_keys(basekeys, patchkeys, name)
            paths[key] = path
        end
        stage = "descriptor"
        descriptor_features = joinpath(outdir, "features_descriptor.tsv")
        augment_descriptor(local_features, paths["--descriptor-patches"], descriptor_features, opts["--config"])
        stage = "molds"
        templates = get(opts, "--templates", joinpath(outdir, "templates_cc.tsv"))
        if !haskey(opts, "--templates")
            run_stage(outdir, stage, "build_cc_molds_native.jl", ["--cube0", abspath(opts["--cube0"]), "--cube1", abspath(opts["--cube1"]),
                "--frame0", abspath(opts["--frame0"]), "--frame1", abspath(opts["--frame1"]), "--config", abspath(opts["--config"]), "--out", templates])
        end
        score_paths = String[]
        stage = "mold_scores"
        for (key, name, prefix) in (("--patches-fwd", "score_fwd", "res"), ("--patches-bwd", "score_bwd", "bwd_res"))
            path = joinpath(outdir, name * ".tsv")
            run_stage(outdir, name, "score_connected_mold_templates.jl", ["--patches", paths[key], "--templates", templates,
                "--template-mode", "contrast", "--prefix", prefix, "--out", path])
            push!(score_paths, path)
        end
        stage = "fisher"
        fisher = joinpath(outdir, "fisher_cv.tsv")
        run_stage(outdir, stage, "build_empirical_fisher_native.jl", ["--patches", paths["--patches-fwd"], "--prefix", "res",
            "--config", abspath(opts["--config"]), "--out", fisher])
        stage = "predictor_features"
        table = joinpath(outdir, "features_predictor.tsv")
        join_predictor_features(descriptor_features, split_features, score_paths[1], score_paths[2], fisher, table)
        stage = "gmm"
        gmm = joinpath(outdir, "pred_gmm.tsv")
        sel = cfg["selection"]
        common = ["--features", table, "--first-seed", string(sel["first_seed"])]
        sel["interactions"] && push!(common, "--interactions")
        run_stage(outdir, stage, "build_labelfree_gmm_predictions.jl", vcat(common,
            ["--out", gmm, "--view", "v_cc=$BASE4,patch_u_asym_reconstructed,mold_cc_fwd,mold_cc_bwd,emp_fisher",
             "--seeds", string(sel["gmm_seeds"]), "--selftrain", string(sel["gmm_selftrain"])]))
        stage = "kmeans"
        km = joinpath(outdir, "pred_kmeans.tsv")
        run_stage(outdir, stage, "build_labelfree_unit_predictions.jl", vcat(common,
            ["--out", km, "--patches", paths["--patches-bwd"], "--view", "v_base=$BASE4",
             "--view", "v_split=$BASE4,split_log_skew", "--view", "v_comt=$BASE4,bwd_neg_com_t",
             "--view", "v_diag45=$BASE4,bwd_neg_diag45", "--seeds", string(sel["kmeans_seeds"])]))
        stage = "vote"
        predictions = joinpath(outdir, "predictions.tsv")
        write_soft_vote(table, km, gmm, predictions, opts["--config"])
        run_stage(outdir, "validation", "validate_unit_predictions.jl", ["--predictions", predictions, "--features", geometry])
        stage = "maps"
        plots_dir = joinpath(outdir, "plots")
        run_stage(outdir, stage, "plot_reconstructed_unit_assignment.jl", ["--features", geometry, "--predictions", predictions, "--outdir", plots_dir])
        stage = "qc"
        run_stage(outdir, stage, "summarize_unknown_unit_qc.jl", ["--predictions", predictions,
            "--plots-dir", plots_dir, "--out", joinpath(outdir, "review_queue.tsv")])
        write_summary(predictions, joinpath(outdir, "summary.tsv"))
        println("Completed ", cfg["model"]["name"], ": ", length(counts), " chains; ", outdir)
    catch err
        reason = replace(sprint(showerror, err), '\n' => ' ', '\t' => ' ', '\r' => ' ')
        write_table(joinpath(outdir, "failures.tsv"), ["file", "N_selected", "stage", "status", "reason"],
            [Dict("file" => f, "N_selected" => (haskey(counts, f) ? string(counts[f]) : "NA"), "stage" => stage,
                  "status" => "incomplete", "reason" => reason) for f in sort(collect(files))])
        rethrow()
    end
end

function main(args=ARGS)
    opts = parse_options(args)
    opts === nothing || execute_pipeline(opts)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

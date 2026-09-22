#!/usr/bin/env julia
# One shared initialization, matched joint/profiled optimization, unchanged heads.
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
include(joinpath(@__DIR__, "refine_fixed_n_geometry.jl"))
const GEOMETRY_OPTIONS = Set(["--data-dir", "--count-config", "--config", "--settings",
    "--features", "--split-features", "--templates", "--outdir"])

function geometry_comparison(args=ARGS; runner=run_stage, exporter=export_features)
    if "--help" in args || "-h" in args
        println("run_geometry_profile_comparison.jl: ", join(sort(collect(GEOMETRY_OPTIONS)), " VALUE "), " VALUE [--dry-run]")
        println("Reference replays saved support. Joint and profiled arms share fresh native initialization at saved N.")
        return
    end
    filtered = String[]; settings_path = ""; i = 1
    while i <= length(args)
        key = args[i]
        key == "--dry-run" && (push!(filtered,key); i += 1; continue)
        key in GEOMETRY_OPTIONS && i < length(args) || error("Missing or forbidden comparison option")
        if key == "--settings"
            isempty(settings_path) || error("Repeated settings"); settings_path = args[i+1]
        else
            append!(filtered,args[i:i+1])
        end
        i += 2
    end
    opts = parse_options(filtered)
    isempty(settings_path) && error("Explicit comparison settings required")
    all(haskey(opts,k) for k in setdiff(GEOMETRY_OPTIONS,Set(["--settings"]))) || error("Missing input")
    _, base = lobe_table(opts["--features"])
    _, splitrows = lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files = Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"]))) == files || error("Raw/cache cohorts differ")
    counts = Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for (key,r) in splitrows
        parse(Int,r["N"]) == counts[first(key)] || error("Split N differs")
    end
    out = abspath(opts["--outdir"])
    refine_args = ["--features",abspath(opts["--features"]),"--data-dir",abspath(opts["--data-dir"]),
        "--config",abspath(opts["--count-config"]),"--settings",abspath(settings_path)]
    FixedNGeometryRefinement.main(vcat(refine_args,["--out",joinpath(out,"geometry.tsv"),"--dry-run"]))
    load_config(opts["--config"])
    haskey(opts,"--dry-run") && return println("Saved reference plus matched joint/profiled arms. Same N, split cache and assignment settings; no computation.")
    mkpath(joinpath(out,"logs")); stage = "inputs"
    try
        basepath = joinpath(out,"cached_features.tsv"); splitpath = joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        common = ["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),"--split-features",splitpath]
        stage = "reference"
        runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",basepath,"--outdir",joinpath(out,stage)]))
        cp(basepath,joinpath(out,stage,"features.tsv")); cp(splitpath,joinpath(out,stage,"features_split.tsv"))
        stage = "geometry"
        mainpath = joinpath(out,"geometry.tsv")
        exporter(out,stage,"refine_fixed_n_geometry.jl",refine_args,mainpath;nfiles=length(files))
        nchunks = min(Threads.nthreads(),4,length(files))
        shards = nchunks == 1 ? [mainpath] : [joinpath(out,"geometry_chunk$i.tsv") for i in 1:nchunks]
        paths = Dict{String,String}()
        for mode in ("native","joint","profiled")
            path = joinpath(out,"$(mode)_features.tsv")
            merge_feature_chunks([p * ".$mode.tsv" for p in shards],path)
            check_counts(path,counts)
            paths[mode] = path
        end
        read(paths["profiled"]) == read(mainpath) || error("Profile shard product differs from companion")
        for mode in ("joint","profiled")
            stage = mode
            runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",paths[mode],"--outdir",joinpath(out,mode)]))
            cp(paths[mode],joinpath(out,mode,"features.tsv")); cp(splitpath,joinpath(out,mode,"features_split.tsv"))
        end
        for mode in ("reference","joint","profiled")
            check_counts(joinpath(out,mode,"predictions.tsv"),counts)
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed fixed-N geometry comparison. External grading remains separate.")
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && geometry_comparison()

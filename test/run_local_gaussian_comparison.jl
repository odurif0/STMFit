#!/usr/bin/env julia
# Saved reference, matched global refit, and valid-GCV expanded orientation pool.
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
include(joinpath(@__DIR__,"refit_local_gaussian_orientation.jl"))
using SHA
const LOCAL_GAUSSIAN_OPTIONS = Set(["--data-dir","--count-config","--config","--settings",
    "--features","--split-features","--templates","--outdir"])

function local_gaussian_comparison(args=ARGS; runner=run_stage, exporter=export_features)
    if "--help" in args || "-h" in args
        println("run_local_gaussian_comparison.jl: ",join(sort(collect(LOCAL_GAUSSIAN_OPTIONS))," VALUE ")," VALUE [--dry-run]")
        println("Saved reference plus global/GCV-local Gaussian refits at fixed N. Split cache, sampling frame and classifier unchanged.")
        return
    end
    opts=Dict{String,String}(); filtered=String[]; i=1
    while i<=length(args)
        key=args[i]; haskey(opts,key) && error("Repeated comparison option")
        if key=="--dry-run"
            opts[key]="true"; push!(filtered,key); i+=1; continue
        end
        key in LOCAL_GAUSSIAN_OPTIONS && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden comparison option")
        opts[key]=args[i+1]; key=="--settings" || append!(filtered,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in LOCAL_GAUSSIAN_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered); load_config(opts["--config"])
    _,base=lobe_table(opts["--features"]); _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files=Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"])))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for (k,row) in splitrows; parse(Int,row["N"])==counts[first(k)] || error("Split N mismatch"); end
    out=abspath(opts["--outdir"])
    refit_args=["--features",abspath(opts["--features"]),"--data-dir",abspath(opts["--data-dir"]),
        "--config",abspath(opts["--count-config"]),"--settings",abspath(opts["--settings"])]
    LocalGaussianOrientation.main(vcat(refit_args,["--out",joinpath(out,"refit.tsv"),"--dry-run"]))
    haskey(opts,"--dry-run") && return println("Three complete arms, $(length(files)) files / $(length(base)) keys; no image read, fit or output.")
    mkpath(joinpath(out,"logs")); stage="inputs"
    try
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(LOCAL_GAUSSIAN_OPTIONS)) if isfile(opts[k])])
        write_table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],
            [Dict("file"=>f,"sha256"=>bytes2hex(sha256(read(p)))) for (f,p) in sort(collect(raw_index(opts["--data-dir"])))])
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),"--split-features",splitpath]
        stage="reference"
        runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",basepath,"--outdir",joinpath(out,stage)]))
        cp(basepath,joinpath(out,stage,"features.tsv")); cp(splitpath,joinpath(out,stage,"features_split.tsv"))
        stage="refit"; mainpath=joinpath(out,"refit.tsv")
        exporter(out,stage,"refit_local_gaussian_orientation.jl",refit_args,mainpath;nfiles=length(files))
        nchunks=min(Threads.nthreads(),4,length(files))
        shards=nchunks==1 ? [mainpath] : [joinpath(out,"refit_chunk$i.tsv") for i in 1:nchunks]
        for (arm,mode) in (("global_refit","global"),("gcv_orientation","gcv"))
            path=joinpath(out,arm*"_features.tsv")
            merge_feature_chunks([p*".$mode.tsv" for p in shards],path)
            selected=check_counts(path,counts)
            PatchFrames.read_model_axes(values(selected)) === nothing && error("Fitted model axes missing")
            for (key,row) in selected
                parse(Int,row["N"])==counts[first(key)] || error("Refit changed N")
            end
            mode=="gcv" && read(path)!=read(mainpath) && error("Main shard product differs from companion")
            stage=arm
            runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",path,"--outdir",joinpath(out,arm)]))
            cp(path,joinpath(out,arm,"features.tsv")); cp(splitpath,joinpath(out,arm,"features_split.tsv"))
        end
        for arm in ("reference","global_refit","gcv_orientation")
            check_counts(joinpath(out,arm,"predictions.tsv"),counts)
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed Gaussian orientation comparison. External grading remains separate.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && local_gaussian_comparison()

#!/usr/bin/env julia
# Three complete inference arms; external grading is deliberately separate.
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
include(joinpath(@__DIR__,"profile_cross_view_residuals.jl"))
using SHA
const CROSS_VIEW_OPTIONS=Set(["--data-dir","--count-config","--config","--settings",
    "--features","--split-features","--templates","--outdir"])

function cross_view_comparison(args=ARGS;runner=run_stage,exporter=export_features)
    if "--help" in args || "-h" in args
        println("run_cross_view_residual_comparison.jl: ",join(sort(collect(CROSS_VIEW_OPTIONS))," VALUE ")," VALUE [--dry-run]")
        println("Saved reference, same-view, opposite-view subtraction. N, geometry, main/split features and classifier fixed.")
        return
    end
    opts=Dict{String,String}(); filtered=String[]; i=1
    while i<=length(args)
        key=args[i]; haskey(opts,key) && error("Repeated comparison option")
        if key=="--dry-run"
            opts[key]="true"; push!(filtered,key); i+=1; continue
        end
        key in CROSS_VIEW_OPTIONS && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden comparison option")
        opts[key]=args[i+1]; key=="--settings" || append!(filtered,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in CROSS_VIEW_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered); load_config(opts["--config"])
    _,base=lobe_table(opts["--features"]); _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files=Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"])))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for (k,row) in splitrows; parse(Int,row["N"])==counts[first(k)] || error("Split N mismatch"); end
    out=abspath(opts["--outdir"])
    profile_args=["--features",abspath(opts["--features"]),"--data-dir",abspath(opts["--data-dir"]),
        "--config",abspath(opts["--count-config"]),"--settings",abspath(opts["--settings"])]
    CrossViewResidual.main(vcat(profile_args,["--out",joinpath(out,"profile_fwd.tsv"),"--dry-run"]))
    haskey(opts,"--dry-run") && return println("Three complete arms, $(length(files)) files / $(length(base)) keys; metadata only, no output.")
    mkpath(joinpath(out,"logs")); stage="inputs"
    try
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(CROSS_VIEW_OPTIONS)) if isfile(opts[k])])
        write_table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],
            [Dict("file"=>f,"sha256"=>bytes2hex(sha256(read(p)))) for (f,p) in sort(collect(raw_index(opts["--data-dir"])))])
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        stage="profile"; fwd=joinpath(out,"profile_fwd.tsv"); bwd=joinpath(out,"profile_bwd.tsv")
        exporter(out,stage,"profile_cross_view_residuals.jl",profile_args,fwd;nfiles=length(files))
        nchunks=min(Threads.nthreads(),4,length(files))
        shards=nchunks==1 ? [fwd] : [joinpath(out,"profile_chunk$i.tsv") for i in 1:nchunks]
        merge_feature_chunks([p*".bwd.tsv" for p in shards],bwd)
        check_counts(fwd,counts); check_counts(bwd,counts)
        ResidualFeatures.read_residual_models(fwd,bwd,values(base))
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),
            "--features",basepath,"--split-features",splitpath]
        for (arm,extras) in (("reference",String[]),
            ("same_view",["--residual-features-fwd",fwd,"--residual-features-bwd",bwd]),
            ("cross_view",["--residual-features-fwd",bwd,"--residual-features-bwd",fwd]))
            stage=arm
            runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,extras,["--outdir",joinpath(out,arm)]))
            cp(basepath,joinpath(out,arm,"features.tsv")); cp(splitpath,joinpath(out,arm,"features_split.tsv"))
            check_counts(joinpath(out,arm,"predictions.tsv"),counts)
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed cross-view residual comparison. External grading remains separate.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && cross_view_comparison()

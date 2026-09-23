#!/usr/bin/env julia
# One saved-support replay and one tangent-CC variant. No external labels.
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
using SHA
const TANGENT_OPTIONS=Set(["--data-dir","--count-config","--config","--settings",
    "--features","--split-features","--templates","--outdir"])

function tangent_comparison(args=ARGS;runner=run_stage)
    if "--help" in args
        println("run_tangent_mold_comparison.jl: ",join(sort(collect(TANGENT_OPTIONS))," VALUE ")," VALUE [--dry-run]")
        return
    end
    opts=Dict{String,String}(); filtered=String[]; i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated comparison option")
        if k=="--dry-run"; opts[k]="true"; push!(filtered,k); i+=1; continue; end
        k in TANGENT_OPTIONS && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden comparison option")
        opts[k]=args[i+1]; k=="--settings" || append!(filtered,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in TANGENT_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered); cfg=load_config(opts["--config"])
    TangentMoldProjection.settings(opts["--settings"])
    cfg["preprocessing"]["patch_residual_filter"]=="smooth_residual" || error("Matched residual required")
    _,base=lobe_table(opts["--features"]); _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    for r in values(base); TangentMoldProjection.validate_geometry(r); end
    files=Set(first.(collect(keys(base)))); raw=raw_index(opts["--data-dir"])
    Set(keys(raw))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for (k,r) in splitrows; parse(Int,r["N"])==counts[first(k)] || error("Split N mismatch"); end
    check_counts(opts["--features"],counts)
    haskey(opts,"--dry-run") && return println("Two complete arms; $(length(files)) files / $(length(base)) keys; metadata only, no output")
    out=abspath(opts["--outdir"]); mkpath(joinpath(out,"logs")); stage="inputs"
    try
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(TANGENT_OPTIONS)) if isfile(opts[k])])
        write_table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],
            [Dict("file"=>f,"sha256"=>bytes2hex(sha256(read(p)))) for (f,p) in sort(collect(raw))])
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),
            "--features",basepath,"--split-features",splitpath]
        for (arm,extras) in (("reference",String[]),("tangent",["--mold-tangent-settings",abspath(opts["--settings"])]))
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
    println("Complete tangent-CC comparison; external grading remains separate")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && tangent_comparison()

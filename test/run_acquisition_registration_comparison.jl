#!/usr/bin/env julia
# Same frozen geometry and classifier in all three arms. No count selection.
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
include(joinpath(@__DIR__,"estimate_acquisition_shifts.jl"))
using SHA
const REGISTRATION_OPTIONS = Set(["--data-dir","--count-config","--config","--settings",
    "--features","--split-features","--templates","--outdir"])

function estimate_shards(out,opts,base;runner=run_stage)
    files=sort(unique(first.(collect(keys(base)))))
    nchunks=min(Threads.nthreads(),4,length(files))
    header,_=lobe_table(opts["--features"])
    paths=String[]
    @sync for i in 1:nchunks
        subset=Set(files[i:nchunks:end]); part=joinpath(out,"registration_features_chunk$i.tsv")
        write_table(part,header,[base[k] for k in sort(collect(keys(base))) if first(k) in subset])
        path=joinpath(out,"shifts_chunk$i.tsv"); push!(paths,path)
        args=["--features",part,"--data-dir",abspath(opts["--data-dir"]),
            "--config",abspath(opts["--count-config"]),"--settings",abspath(opts["--settings"]),"--out",path]
        @async runner(out,"registration_chunk$i","estimate_acquisition_shifts.jl",args;threads=1)
    end
    for suffix in ("",".summary.tsv",".peaks.tsv",".scores.tsv")
        merged=Dict{String,String}[]; expected=nothing
        for p in paths
            h,rs=read_table(p*suffix)
            expected===nothing ? (expected=h) : (expected==h || error("Registration shard headers differ"))
            append!(merged,rs)
        end
        sort!(merged;by=r->(r["file"],parse(Int,get(r,"band","0")),parse(Int,get(r,"dx_px","0"))))
        write_table(joinpath(out,"shifts.tsv")*suffix,expected,merged)
    end
    EstimateAcquisitionShifts.PatchAcquisition.read_shifts(joinpath(out,"shifts.tsv"),files)
end

function registration_comparison(args=ARGS;runner=run_stage,estimator=estimate_shards)
    if "--help" in args || "-h" in args
        println("run_acquisition_registration_comparison.jl: ",join(sort(collect(REGISTRATION_OPTIONS))," VALUE ")," VALUE [--dry-run]")
        println("Saved reference, raw-observation-mask control, and gated integer-x registration; fixed N/geometry/classifier.")
        return
    end
    opts=Dict{String,String}(); filtered=String[]; i=1
    while i<=length(args)
        k=args[i]
        if k=="--dry-run"
            haskey(opts,k) && error("Repeated dry run"); opts[k]="true"; push!(filtered,k); i+=1; continue
        end
        k in REGISTRATION_OPTIONS && i<length(args) || error("Missing/forbidden registration option")
        haskey(opts,k) && error("Repeated registration option")
        opts[k]=args[i+1]; k=="--settings" || append!(filtered,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in REGISTRATION_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered)
    EstimateAcquisitionShifts.AcquisitionRegistration.load_settings(opts["--settings"])
    load_config(opts["--config"])
    _,base=lobe_table(opts["--features"]); _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files=Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"])))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for table in (base,splitrows), (k,row) in table
        parse(Int,row["N"])==counts[first(k)] || error("Saved N differs from lobe count")
    end
    haskey(opts,"--dry-run") && return println("Dry run: $(length(files)) files, $(length(base)) keys; no compute/output. Three frozen-geometry arms, at most four registration workers.")
    out=abspath(opts["--outdir"]); mkpath(joinpath(out,"logs")); stage="inputs"
    try
        inputs=[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(REGISTRATION_OPTIONS)) if isfile(opts[k])]
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],inputs)
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),
            "--features",basepath,"--split-features",splitpath]
        stage="registration"; estimator(out,opts,base;runner)
        zero=joinpath(out,"zero_shifts.tsv")
        write_table(zero,["file","bwd_sample_dx_px"],[Dict("file"=>f,"bwd_sample_dx_px"=>"0") for f in sort(collect(files))])
        for mode in ("reference","observed_control","registered")
            stage=mode
            extra=mode=="reference" ? String[] : ["--acquisition-shifts",mode=="observed_control" ? zero : joinpath(out,"shifts.tsv")]
            runner(out,mode,"run_reconstructed_chitosan.jl",vcat(common,extra,["--outdir",joinpath(out,mode)]))
            cp(basepath,joinpath(out,mode,"features.tsv")); cp(splitpath,joinpath(out,mode,"features_split.tsv"))
            check_counts(joinpath(out,mode,"predictions.tsv"),counts)
        end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed acquisition comparison. External benchmark grading remains separate.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && registration_comparison()

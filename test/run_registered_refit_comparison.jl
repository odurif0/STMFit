#!/usr/bin/env julia
# One registered-refit hypothesis with matched zero-shift and saved controls.
include(joinpath(@__DIR__,"run_acquisition_registration_comparison.jl"))
include(joinpath(@__DIR__,"refit_registered_geometry.jl"))
const REFIT_OPTIONS=union(REGISTRATION_OPTIONS,Set(["--refit-settings"]))

function registered_refit_comparison(args=ARGS;runner=run_stage,estimator=estimate_shards,exporter=export_features)
    if "--help" in args || "-h" in args
        return println("run_registered_refit_comparison.jl: ",join(sort(collect(REFIT_OPTIONS))," VALUE ")," VALUE [--dry-run]")
    end
    opts=Dict{String,String}(); filtered=String[]; i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option")
        if k=="--dry-run"; opts[k]="true"; push!(filtered,k); i+=1; continue; end
        k in REFIT_OPTIONS && i<length(args) || error("Missing/forbidden refit comparison option")
        opts[k]=args[i+1]; k in ("--settings","--refit-settings") || append!(filtered,args[i:i+1]); i+=2
    end
    all(haskey(opts,k) for k in REFIT_OPTIONS) || error("All comparison inputs required")
    parse_options(filtered)
    EstimateAcquisitionShifts.AcquisitionRegistration.load_settings(opts["--settings"])
    refit_options=RegisteredRefit.settings(opts["--refit-settings"])
    load_config(opts["--config"])
    _,base=lobe_table(opts["--features"];required=RegisteredRefit.F.REQUIRED)
    _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files=Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"])))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->first(k)==f,keys(base)) for f in files)
    for f in files
        RegisteredRefit.F.validate_chain([base[k] for k in sort(collect(keys(base))) if first(k)==f])
    end
    all(parse(Int,r["N"])==counts[first(k)] for (k,r) in splitrows) || error("Split count differs")
    haskey(opts,"--dry-run") && return println("Dry run: $(length(files)) files / $(length(base)) keys; native Gaussian/split refits, saved N, original_support=$(refit_options.original_support), unchanged assignment; no pixels/output.")
    out=abspath(opts["--outdir"]); mkpath(joinpath(out,"logs")); stage="inputs"
    try
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(REFIT_OPTIONS)) if isfile(opts[k])])
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        cp(opts["--refit-settings"],joinpath(out,"refit_settings.toml"))
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"])]
        stage="registration"; estimator(out,opts,base;runner)
        zero=joinpath(out,"zero_shifts.tsv")
        write_table(zero,["file","bwd_sample_dx_px"],[Dict("file"=>f,"bwd_sample_dx_px"=>"0") for f in sort(collect(files))])
        stage="reference"
        runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",basepath,"--split-features",splitpath,"--outdir",joinpath(out,stage)]))
        cp(basepath,joinpath(out,stage,"features.tsv")); cp(splitpath,joinpath(out,stage,"features_split.tsv"))
        stage="refit"; mainpath=joinpath(out,"refit.tsv")
        refit_args=["--data-dir",abspath(opts["--data-dir"]),"--features",basepath,
            "--config",abspath(opts["--count-config"]),"--assignment-config",abspath(opts["--config"]),
            "--settings",abspath(opts["--refit-settings"]),"--shifts",joinpath(out,"shifts.tsv")]
        exporter(out,stage,"refit_registered_geometry.jl",refit_args,mainpath;nfiles=length(files))
        nchunks=min(Threads.nthreads(),4,length(files))
        shards=nchunks==1 ? [mainpath] : [joinpath(out,"refit_chunk$i.tsv") for i in 1:nchunks]
        for arm in RegisteredRefit.ARMS
            paths=Dict{String,String}()
            for profile in RegisteredRefit.PROFILES
                path=joinpath(out,"$arm.$profile.tsv")
                merge_feature_chunks([p*".$arm.$profile.tsv" for p in shards],path)
                check_counts(path,counts); paths[profile]=path
            end
            stage=arm
            runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",paths["gaussian"],
                "--split-features",paths["split"],"--acquisition-shifts",arm=="control" ? zero : joinpath(out,"shifts.tsv"),
                "--outdir",joinpath(out,arm)]))
            cp(paths["gaussian"],joinpath(out,arm,"features.tsv")); cp(paths["split"],joinpath(out,arm,"features_split.tsv"))
        end
        read(joinpath(out,"registered.gaussian.tsv"))==read(mainpath) || error("Refit shard output differs")
        for arm in ("reference","control","registered"); check_counts(joinpath(out,arm,"predictions.tsv"),counts); end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed registered native refit comparison. External grading remains separate.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && registered_refit_comparison()

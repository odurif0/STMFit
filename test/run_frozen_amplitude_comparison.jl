#!/usr/bin/env julia
# Two assignments on the identical saved N/geometry; one bounded linear change.
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
include(joinpath(@__DIR__,"profile_frozen_amplitudes.jl"))
const PROFILE_OPTIONS=Set(["--data-dir","--count-config","--config","--settings",
    "--features","--split-features","--templates","--outdir"])

function profile_comparison(args=ARGS; runner=run_stage, exporter=export_features)
    if "--help" in args
        println("run_frozen_amplitude_comparison.jl: ",join(sort(collect(PROFILE_OPTIONS))," VALUE ")," VALUE [--dry-run]")
        return
    end
    filtered=String[]; settings=""; i=1
    while i<=length(args)
        key=args[i]
        key=="--dry-run" && (push!(filtered,key); i+=1; continue)
        key in PROFILE_OPTIONS && i<length(args) || error("Missing/forbidden comparison option")
        if key=="--settings"
            isempty(settings) || error("Repeated --settings"); settings=args[i+1]
        else
            append!(filtered,args[i:i+1])
        end
        i+=2
    end
    opts=parse_options(filtered)
    isempty(settings) && error("Explicit linear settings required")
    all(haskey(opts,k) for k in setdiff(PROFILE_OPTIONS,Set(["--settings"]))) || error("Missing experiment input")
    _,base=lobe_table(opts["--features"]); _,splitrows=lobe_table(opts["--split-features"])
    require_same_keys(base,splitrows,"split geometry")
    files=Set(first.(collect(keys(base))))
    Set(keys(raw_index(opts["--data-dir"])))==files || error("Raw/cache cohorts differ")
    counts=Dict(f=>count(k->k[1]==f,keys(base)) for f in files)
    for (key,r) in splitrows; parse(Int,r["N"])==counts[key[1]] || error("Split N mismatch"); end
    out=abspath(opts["--outdir"])
    profile_args=["--features",abspath(opts["--features"]),"--data-dir",abspath(opts["--data-dir"]),
        "--config",abspath(opts["--count-config"]),"--settings",abspath(settings)]
    FrozenAmplitudeProfile.main(vcat(profile_args,["--out",joinpath(out,"profiled_features.tsv"),"--dry-run"]))
    load_config(opts["--config"])
    haskey(opts,"--dry-run") && return println("Two arms, frozen N/geometry, no nonlinear fit; output $out")
    mkpath(joinpath(out,"logs")); stage="inputs"
    try
        basepath=joinpath(out,"cached_features.tsv"); splitpath=joinpath(out,"cached_split.tsv")
        cp(opts["--features"],basepath); cp(opts["--split-features"],splitpath)
        common=["--data-dir",abspath(opts["--data-dir"]),"--count-config",abspath(opts["--count-config"]),
            "--config",abspath(opts["--config"]),"--templates",abspath(opts["--templates"]),"--split-features",splitpath]
        stage="control"
        runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",basepath,"--outdir",joinpath(out,stage)]))
        cp(basepath,joinpath(out,stage,"features.tsv")); cp(splitpath,joinpath(out,stage,"features_split.tsv"))
        stage="profile"
        profiled=joinpath(out,"profiled_features.tsv")
        exporter(out,stage,"profile_frozen_amplitudes.jl",profile_args,profiled;nfiles=length(files))
        updated=check_counts(profiled,counts)
        for k in keys(base),c in setdiff(keys(base[k]),FrozenAmplitudeProfile.MUTABLE_COLUMNS)
            base[k][c]==updated[k][c] || error("Frozen geometry changed: $k $c")
        end
        stage="profiled"
        runner(out,stage,"run_reconstructed_chitosan.jl",vcat(common,["--features",profiled,"--outdir",joinpath(out,stage)]))
        cp(profiled,joinpath(out,stage,"features.tsv")); cp(splitpath,joinpath(out,stage,"features_split.tsv"))
        for arm in ("control","profiled"); check_counts(joinpath(out,arm,"predictions.tsv"),counts); end
    catch err
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Completed frozen-amplitude comparison. External grading is separate.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && profile_comparison()

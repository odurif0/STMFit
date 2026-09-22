#!/usr/bin/env julia
# Merge the eight numerical cases, never a benchmark or classifier.
include(joinpath(@__DIR__,"diagnose_paired_convergence.jl"))
const D=PairedConvergenceDiagnostic
function main(args=ARGS)
    length(args)==1 || error("Usage: merge_paired_convergence.jl NEW_RUN_DIR")
    out=args[1]; dirs=[joinpath(out,"chunk$i") for i in 1:4]
    all(isdir,dirs) || error("Four shards required")
    fits=Dict{String,String}[]
    for dir in dirs
        isfile(joinpath(dir,"failures.tsv")) && error("A diagnostic shard failed")
        append!(fits,D.table(joinpath(dir,"fits.tsv")))
    end
    expected=Set((c.profile,c.family,c.mode,s) for c in D.CASES for s in ("control","extended"))
    keys=[(r["profile"],r["family"],r["mode"],r["stage"]) for r in fits]
    length(keys)==length(Set(keys))==16 && Set(keys)==expected || error("Incomplete/duplicate cases")
    for name in ("input_hashes.tsv","settings.toml","model_settings.toml")
        all(read(joinpath(dir,name))==read(joinpath(first(dirs),name)) for dir in dirs) || error("Shard inputs differ")
        cp(joinpath(first(dirs),name),joinpath(out,name))
    end
    sort!(fits;by=r->(r["profile"],r["family"],r["mode"],r["stage"]))
    D.write_records(joinpath(out,"fits.tsv"),fits)
    eligibility=Dict{String,String}[]
    for stage in ("control","extended"),mode in ("fused","paired"),profile in ("gaussian","split")
        valid=sort(filter(r->r["stage"]==stage && r["mode"]==mode && r["profile"]==profile && r["valid"]=="true",fits);by=r->D.value(r,"gcv"))
        push!(eligibility,Dict("stage"=>stage,"mode"=>mode,"profile"=>profile,"valid_families"=>string(length(valid)),
            "selected_family"=>isempty(valid) ? "none" : first(valid)["family"],
            "selected_lm_converged"=>isempty(valid) ? "NA" : first(valid)["converged"],
            "selected_stationary"=>isempty(valid) ? "NA" : first(valid)["stationarity_passed"]))
    end
    D.write_records(joinpath(out,"eligibility.tsv"),eligibility)
    println("All eight cases / sixteen fits accounted for. Convergence and validity remain separate; no grade.")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

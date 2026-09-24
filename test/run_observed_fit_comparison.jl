#!/usr/bin/env julia
# One comparison; all files must fit before either cohort-wide assignment runs.
module ObservedFitComparison
include(joinpath(@__DIR__,"refit_observed_geometry.jl"))
using SHA, TOML
const P=ObservedGeometry.Pipeline
const REQUIRED=Set(["--data-dir","--count-config","--config","--templates","--selected-summary","--settings","--outdir"])

function options(args)
    args==["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in REQUIRED && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option $k")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in REQUIRED) || error("All comparison inputs required")
    isdir(o["--data-dir"]) || error("Raw directory missing")
    all(isfile(o[k]) for k in setdiff(REQUIRED,Set(["--data-dir","--outdir"]))) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    o
end

function candidate_config(raw, settings)
    cfg=deepcopy(raw)
    cfg["preprocessing"]["missing_pixel_policy"]="observed_only"
    cfg["preprocessing"]["plane_rank_rtol"]=settings.plane_rank_rtol
    cfg
end

function execute(o; runner=P.run_stage, exporter=P.export_features)
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    input=ObservedGeometry.inputs(Dict("--data-dir"=>o["--data-dir"],"--selected-summary"=>o["--selected-summary"],
        "--config"=>o["--count-config"],"--assignment-config"=>o["--config"]))
    settings=ObservedGeometry.settings(o["--settings"])
    counts=input.counts; files=sort(collect(keys(counts)))
    haskey(o,"--dry-run") && return println("Dry run: $(length(files)) files / $(sum(values(counts))) lobes; fixed saved N; native versus observed-only fit and patches; no pixels or output")
    out=abspath(o["--outdir"]); mkpath(joinpath(out,"logs")); stage="inputs"
    try
        P.write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(o[k])))) for k in sort(collect(REQUIRED)) if isfile(o[k])])
        P.write_table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],
            [Dict("file"=>f,"sha256"=>bytes2hex(sha256(read(input.paths[f])))) for f in files])
        for (key,name) in (("--selected-summary","selected_summary.tsv"),("--count-config","control_count.toml"),
                ("--config","assignment.toml"),("--settings","settings.toml"),("--templates","templates.tsv"))
            cp(o[key],joinpath(out,name))
        end
        observed_cfg=joinpath(out,"observed_count.toml")
        open(io->TOML.print(io,candidate_config(input.raw,settings)),observed_cfg,"w")
        mainpath=joinpath(out,"refit.tsv"); stage="refit"
        args=["--data-dir",abspath(o["--data-dir"]),"--selected-summary",abspath(o["--selected-summary"]),
            "--config",abspath(o["--count-config"]),"--assignment-config",abspath(o["--config"]),"--settings",abspath(o["--settings"])]
        exporter(out,stage,"refit_observed_geometry.jl",args,mainpath;nfiles=length(files))
        chunks=min(Threads.nthreads(),4,length(files))
        shards=chunks==1 ? [mainpath] : [joinpath(out,"refit_chunk$i.tsv") for i in 1:chunks]
        for arm in ObservedGeometry.ARMS
            paths=Dict{String,String}()
            for profile in ObservedGeometry.PROFILES
                path=joinpath(out,"$arm.$profile.tsv")
                P.merge_feature_chunks([s*".$arm.$profile.tsv" for s in shards],path)
                P.check_counts(path,counts); paths[profile]=path
            end
            stage=arm
            runner(out,stage,"run_reconstructed_chitosan.jl",[
                "--data-dir",abspath(o["--data-dir"]),"--count-config",arm=="control" ? abspath(o["--count-config"]) : observed_cfg,
                "--config",abspath(o["--config"]),"--templates",abspath(o["--templates"]),
                "--selected-summary",abspath(o["--selected-summary"]),"--features",paths["gaussian"],
                "--split-features",paths["split"],"--outdir",joinpath(out,arm)])
            for (profile,name) in (("gaussian","features.tsv"),("split","features_split.tsv"))
                cp(paths[profile],joinpath(out,arm,name))
            end
            P.check_counts(joinpath(out,arm,"predictions.tsv"),counts)
        end
        read(mainpath)==read(joinpath(out,"observed.gaussian.tsv")) || error("Merged refit differs")
    catch err
        P.write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Complete observed-fit comparison; grading remains external. No calibration or champion change.")
end
function main(args=ARGS)
    o=options(args)
    o===nothing ? println("run_observed_fit_comparison.jl ",join(sort(collect(REQUIRED))," VALUE ")," VALUE [--dry-run]") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ObservedFitComparison.main()

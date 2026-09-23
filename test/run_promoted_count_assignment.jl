#!/usr/bin/env julia
# Compare two saved label-free count policies, rebuilding both geometries.
# This does not revalidate/recompute the historical counting-policy benchmark.
module PromotedCountAssignment
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
using SHA
const REQUIRED = Set(["--data-dir","--count-config","--config","--templates",
    "--control-summary","--promoted-summary","--outdir"])

function options(args)
    args == ["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in REQUIRED && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option $k")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in REQUIRED) || error("All comparison inputs required")
    isdir(o["--data-dir"]) || error("Missing raw directory")
    all(isfile(o[k]) for k in setdiff(REQUIRED,Set(["--data-dir","--outdir"]))) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    o
end

function inputs(o)
    cfg=TOML.parsefile(o["--count-config"])
    cfg["model"]["selection_criterion"]=="gcv" || error("GCV required")
    cfg["model"]["selection_policy"]=="support_midpoint_hybrid" || error("Unchanged promoted count config required")
    assignment=load_config(o["--config"])
    assignment["preprocessing"]["patch_residual_filter"]=="smooth_residual" || error("Symmetric residual processing required")
    control=selected_counts(o["--control-summary"])
    promoted=selected_counts(o["--promoted-summary"])
    raw=raw_index(o["--data-dir"])
    Set(keys(control))==Set(keys(promoted))==Set(keys(raw)) || error("Count/raw cohorts differ")
    header,rows=read_table(o["--promoted-summary"])
    "selection_policy" in header || error("Promoted summary must name its policy")
    all(r["selection_policy"]=="support_midpoint_hybrid" for r in rows) || error("Mixed/wrong promoted policy")
    for key in ("--control-summary","--promoted-summary")
        _,rs=read_table(o[key])
        any(startswith(get(r,k,""),"adaptive_support") for r in rs for k in
            ("refined_policy","refined_source","selection_policy","selection_source")) &&
            error("Adaptive-support counts are outside this fixed-support comparison")
    end
    (;control,promoted,raw)
end

function execute(o;runner=run_stage)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    r=inputs(o); files=sort(collect(keys(r.control)))
    changed=count(f->r.control[f]!=r.promoted[f],files)
    println("Saved counts only: $(length(files)) files, $changed changed; $(sum(values(r.control))) -> $(sum(values(r.promoted))) lobes")
    haskey(o,"--dry-run") && return println("Two fresh fixed-N arms; no counting sweep, cache trimming, pixels or output")
    out=abspath(o["--outdir"]); mkpath(joinpath(out,"logs")); stage="inputs"
    try
        raw=prepare_raw(o["--data-dir"],out,Set(files))
        write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
            [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(o[k])))) for k in sort(collect(REQUIRED)) if isfile(o[k])])
        write_table(joinpath(out,"raw_hashes.tsv"),["file","sha256"],
            [Dict("file"=>f,"sha256"=>bytes2hex(sha256(read(r.raw[f])))) for f in files])
        cp(o["--control-summary"],joinpath(out,"original_control_summary.tsv"))
        cp(o["--promoted-summary"],joinpath(out,"original_promoted_summary.tsv"))
        write_table(joinpath(out,"count_changes.tsv"),["file","control_N","promoted_N","delta"],
            [Dict("file"=>f,"control_N"=>string(r.control[f]),"promoted_N"=>string(r.promoted[f]),
                "delta"=>string(r.promoted[f]-r.control[f])) for f in files])
        for (arm,counts) in (("control",r.control),("hybrid",r.promoted))
            selected=joinpath(out,arm*"_counts.tsv")
            write_table(selected,["filepath","N_selected"],
                [Dict("filepath"=>joinpath(raw,f),"N_selected"=>string(counts[f])) for f in files])
            stage=arm
            runner(out,stage,"run_reconstructed_chitosan.jl",[
                "--data-dir",raw,"--count-config",abspath(o["--count-config"]),
                "--config",abspath(o["--config"]),"--templates",abspath(o["--templates"]),
                "--selected-summary",selected,"--outdir",joinpath(out,arm)])
            for file in ("features.tsv","features_split.tsv","predictions.tsv")
                check_counts(joinpath(out,arm,file),counts)
            end
        end
    catch e
        write_table(joinpath(out,"failures.tsv"),["stage","reason"],
            [Dict("stage"=>stage,"reason"=>replace(sprint(showerror,e),'\n'=>' ','\t'=>' '))])
        rethrow()
    end
    println("Complete comparison; grading is external. Timed optimizer variability must be reported.")
end
function main(args=ARGS)
    o=options(args)
    o===nothing ? println("run_promoted_count_assignment.jl --data-dir DIR --count-config TOML --config TOML --templates TSV --control-summary TSV --promoted-summary TSV --outdir NEW_DIR [--dry-run]") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && PromotedCountAssignment.main()

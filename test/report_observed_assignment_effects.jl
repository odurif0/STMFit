#!/usr/bin/env julia
# Saved-output attribution only: no labels, fitting or model selection.
module ObservedAssignmentEffects
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
const TABLES=("features","features_split","features_local","patches_fwd17","patches_bwd17","patches_bwd9",
    "features_descriptor","features_predictor","training_support","score_fwd","score_bwd","fisher_cv","pred_gmm","pred_kmeans","predictions")
const LOCAL_INPUTS=Set(["features","features_split","features_local","patches_fwd17","patches_bwd17","patches_bwd9"])

function verify_vote(g,k,p,threshold,model)
    Set(keys(g))==Set(keys(k))==Set(keys(p)) || error("Vote component keys differ")
    for (key,r) in p
        r["model"]==model || error("Unexpected assignment model")
        if g[key]["predicted"]=="?" || k[key]["predicted"]=="?"
            r["predicted"]=="?" && r["probability_1"]=="NA" || error("Missing component did not abstain")
        else
            a=parse(Float64,g[key]["probability_1"]); b=parse(Float64,k[key]["probability_1"])
            0<=a<=1 && 0<=b<=1 || error("Invalid component probability")
            v=(a+b)/2
            isapprox(parse(Float64,r["probability_1"]),v;atol=5.1e-9,rtol=0) || error("Incorrect mean probability")
            isapprox(parse(Float64,r["confidence"]),2abs(v-.5);atol=5.1e-9,rtol=0) || error("Incorrect vote confidence")
            r["predicted"]==(v>=threshold ? "1" : "0") || error("Incorrect threshold decision")
        end
    end
    length(p)
end

function differences(a,b,files; prediction=false)
    Set(keys(a))==Set(keys(b)) || error("Stage cohorts differ")
    keys0=sort([key for key in keys(a) if first(key) in files])
    changed=[key for key in keys0 if prediction ? a[key]["predicted"]!=b[key]["predicted"] : a[key]!=b[key]]
    (;lobes=length(keys0),changed_lobes=length(changed),changed_scans=length(unique(first.(changed))))
end

function analyze(root)
    isfile(joinpath(root,"failures.tsv")) && error("Comparison failed; no partial report")
    counts=selected_counts(joinpath(root,"selected_summary.tsv")); files=Set(keys(counts))
    audits=vcat([last(read_table(p)) for p in readdir(root;join=true) if startswith(basename(p),"refit") && endswith(p,".fits.tsv")]...)
    slots=Set((r["file"],r["arm"],r["profile"]) for r in audits)
    slots==Set((f,a,p) for f in files for a in ("control","observed") for p in ("gaussian","split")) || error("Incomplete fit slots")
    length(slots)==length(audits) && all(r["valid"]=="true" for r in audits) || error("Duplicate or invalid fits")
    reused=Set(r["file"] for r in audits if r["reused_control"]=="true")
    for f in files
        candidate=[r for r in audits if r["file"]==f && r["arm"]=="observed"]
        all((r["reused_control"]=="true")== (f in reused) for r in candidate) || error("Mixed reuse within file")
        all((parse(Int,r["observed_pixels"])==parse(Int,r["total_pixels"]))==(f in reused) for r in candidate) || error("Reuse/mask mismatch")
    end
    groups=("fully_observed"=>reused,"partial"=>setdiff(files,reused))
    rows=Dict{String,String}[]
    cfg=load_config(joinpath(root,"assignment.toml"))
    for arm in ("control","observed")
        dir=joinpath(root,arm)
        isfile(joinpath(dir,"failures.tsv")) && error("Assignment failed")
        tables=Dict(t=>check_counts(joinpath(dir,t*".tsv"),counts) for t in ("pred_gmm","pred_kmeans","predictions"))
        verify_vote(tables["pred_gmm"],tables["pred_kmeans"],tables["predictions"],cfg["selection"]["vote_threshold"],cfg["model"]["name"])
        _,summary=read_table(joinpath(dir,"summary.tsv"))
        length(summary)==length(counts) && Dict(r["file"]=>parse(Int,r["N_selected"]) for r in summary)==counts || error("Final N changed")
        length(readdir(joinpath(dir,"plots","standalone")))==length(counts) || error("Missing maps")
    end
    for stage in TABLES
        a=check_counts(joinpath(root,"control",stage*".tsv"),counts)
        b=check_counts(joinpath(root,"observed",stage*".tsv"),counts)
        for (group,subset) in groups
            d=differences(a,b,subset)
            stage in LOCAL_INPUTS && group=="fully_observed" && d.changed_lobes!=0 && error("Unchanged local inputs differ: $stage")
            labels=stage in ("pred_gmm","pred_kmeans","predictions") ? differences(a,b,subset;prediction=true) : nothing
            push!(rows,Dict("stage"=>stage,"group"=>group,"scans"=>string(length(subset)),"lobes"=>string(d.lobes),
                "changed_rows"=>string(d.changed_lobes),"changed_scans"=>string(d.changed_scans),
                "changed_assignments"=>labels===nothing ? "NA" : string(labels.changed_lobes),
                "changed_assignment_scans"=>labels===nothing ? "NA" : string(labels.changed_scans)))
        end
    end
    rows
end

function main(args=ARGS)
    args==["--help"] && return println("report_observed_assignment_effects.jl RUN_DIR NEW_OUTDIR")
    length(args)==2 || error("Run directory and new report directory required")
    (ispath(args[2]) || islink(args[2])) && error("Report exists")
    rows=analyze(args[1]); mkpath(args[2])
    write_table(joinpath(args[2],"stage_changes.tsv"),["stage","group","scans","lobes","changed_rows","changed_scans","changed_assignments","changed_assignment_scans"],rows)
    open(joinpath(args[2],"report.md"),"w") do io
        println(io,"# Assignment sensitivity at exactly reused local inputs\n\nBoth complete cohorts and independent vote arithmetic pass. No benchmark labels were read.\n")
        println(io,"| Stage | Raw group | Changed rows / lobes | Changed scans | Changed assignments |\n|---|---|---:|---:|---:|")
        for r in rows
            println(io,"| $(r["stage"]) | $(r["group"]) | $(r["changed_rows"])/$(r["lobes"]) | $(r["changed_scans"])/$(r["scans"]) | $(r["changed_assignments"]) |")
        end
        println(io,"\nThe fully observed group's Gaussian/split fits, local features and three patch tables are exactly unchanged. Changes later in this group therefore arise downstream of those local inputs when the rest of the cohort changes. This identifies cohort dependence, not one uniquely responsible learning stage, chemical accuracy or an optimum. Independent learner repeats would be needed to separate cohort sensitivity from any learner nondeterminism.")
    end
    for r in rows
        r["stage"]=="predictions" && println(r)
    end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ObservedAssignmentEffects.main()

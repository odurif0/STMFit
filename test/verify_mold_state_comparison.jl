#!/usr/bin/env julia
# Independent saved-output arithmetic, before and without external grading.
using Test, SHA
include(joinpath(@__DIR__,"lib","reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
length(ARGS)==2 || error("SAVED_TANGENT_DIR COMPARISON_DIR required")
ref,out=ARGS
_,base=lobe_table(joinpath(ref,"features.tsv"))
files=sort(unique(first.(keys(base))))
costs=Dict()
for view in ("fwd","bwd")
    for path in sort(readdir(ref;join=true))
        occursin(Regex("^score_"*view*"(?:_chunk[1-9][0-9]*)?\\.tsv\\.audit\\.tsv\$"),basename(path)) || continue
        _,rows=read_table(path)
        for r in rows
            k=(view,r["file"],parse(Int,r["lobe"]),parse(Int,r["parity"]),parse(Int,r["mirror"]))
            haskey(costs,k) && error("Duplicate input state")
            costs[k]=parse.(Float64,[r["cost0"],r["cost1"]])
        end
    end
end
totals=Dict()
for view in ("fwd","bwd"),f in files,p in 0:1,m in 0:1
    keys_f=sort([k for k in keys(base) if k[1]==f])
    total=0.0
    for (_,i) in keys_f
        c=costs[(view,f,i,mod(i-1+p,2),m)]
        if all(isfinite,c); total+=min(c...); else; @assert c==[Inf,Inf]; end
    end
    totals[(view,f,p,m)]=total
end
old_tables=Dict(name=>lobe_table(joinpath(ref,name))[2] for name in
    ("score_fwd.tsv","score_bwd.tsv","features_predictor.tsv","predictions.tsv"))
@testset "Whole-cohort decoder comparison, no benchmark data" begin
    @test length(costs)==8length(base)
    _,hashes=read_table(joinpath(out,"input_hashes.tsv"))
    # Input paths refer to the compute host; compare cached files by basename.
    for r in hashes
        p=joinpath(ref,basename(r["input"]))
        isfile(p) || continue
        @test bytes2hex(sha256(read(p)))==r["sha256"]
    end
    for mode in ("reference","finite","shared")
        dir=joinpath(out,mode)
        _,pred=lobe_table(joinpath(dir,"predictions.tsv"))
        @test Set(keys(pred))==Set(keys(base))
        old=old_tables["predictions.tsv"]
        @test [k for k in sort(collect(keys(pred))) if pred[k]["predicted"]=="?"]==
              [k for k in sort(collect(keys(old))) if old[k]["predicted"]=="?"]
        for name in ("features.tsv","features_split.tsv","features_local.tsv","features_descriptor.tsv",
            "patches_fwd17.tsv","patches_bwd17.tsv","patches_bwd9.tsv","fisher_cv.tsv","training_support.tsv","pred_kmeans.tsv")
            @test read(joinpath(dir,name))==read(joinpath(ref,name))
        end
        if mode=="reference"
            for name in ("score_fwd.tsv","score_bwd.tsv","features_predictor.tsv","pred_gmm.tsv","predictions.tsv","summary.tsv")
                @test read(joinpath(dir,name))==read(joinpath(ref,name))
            end
        end
        h,features=lobe_table(joinpath(dir,"features_predictor.tsv"))
        @test Set(keys(features))==Set(keys(base))
        for k in keys(base),c in setdiff(h,["mold_cc_fwd","mold_cc_bwd"])
            @test features[k][c]==old_tables["features_predictor.tsv"][k][c]
        end
        scores=Dict(view=>lobe_table(joinpath(dir,"score_"*view*".tsv"))[2] for view in ("fwd","bwd"))
        for f in files
            keys_f=sort([k for k in keys(base) if k[1]==f]); n=length(keys_f)
            selected=Dict()
            for view in ("fwd","bwd")
                s=scores[view][(f,1)]
                p=mod(parse(Int,s["global_phase"])+parse(Int,s["global_direction"])*(n-1),2)
                m=parse(Int,s["global_mirror"]); selected[view]=(p,m)
                for k in keys_f
                    r=scores[view][k]
                    @test r["global_phase"]==s["global_phase"] && r["global_mirror"]==s["global_mirror"]
                    c=costs[(view,f,k[2],mod(k[2]-1+p,2),m)]
                    # Saved costs use eight significant digits; compare within serialization error.
                    for (col,v) in zip(("cost_GlcN","cost_GlcNAc"),c)
                        x=parse(Float64,r[col]); @test isequal(x,v) || isapprox(x,v;rtol=1e-7,atol=0)
                    end
                    margin=abs(c[1]-c[2]); got=parse(Float64,r["cost_margin"])
                    @test isequal(got,margin) || isapprox(got,margin;rtol=1e-7,atol=0)
                    if mode!="reference" && !all(isfinite,c)
                        @test r["predicted"]=="-1" && r["physical_label"]=="?"
                    end
                end
                if mode!="reference"
                    @test isapprox(parse(Float64,s["file_cost"]),totals[(view,f,p,m)];rtol=1e-7,atol=1e-15)
                    if mode=="finite"
                        @test totals[(view,f,p,m)]≈minimum(totals[(view,f,q,t)] for q in 0:1 for t in 0:1) atol=1e-12
                        if all(all(isfinite,costs[(view,f,k[2],0,0)]) for k in keys_f)
                            @test all(scores[view][k]==old_tables["score_"*view*".tsv"][k] for k in keys_f)
                        end
                    end
                end
            end
            if mode=="shared"
                @test selected["fwd"]==selected["bwd"]
                p,m=selected["fwd"]
                @test sum(totals[(v,f,p,m)] for v in ("fwd","bwd"))≈
                    minimum(sum(totals[(v,f,q,t)] for v in ("fwd","bwd")) for q in 0:1 for t in 0:1) atol=1e-12
            end
        end
        changed=[k for k in keys(base) if pred[k]["predicted"]!=old[k]["predicted"]]
        println(mode,": ",length(files)," scans, ",length(base)," keys, ",length(changed)," final changes on ",length(unique(first.(changed)))," scans")
    end
end

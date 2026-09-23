#!/usr/bin/env julia
# Independent saved-cost arithmetic; no production decoder or benchmark labels.
using Test, SHA, Printf
include(joinpath(@__DIR__,"lib/reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
length(ARGS)==2 || error("FINITE_TANGENT_DIR COMPARISON_DIR required")
ref,out=ARGS
_,base=lobe_table(joinpath(ref,"features.tsv")); keys_all=sort(collect(keys(base)))
files=sort(unique(first.(keys_all)))
num(r,c)=parse(Float64,r[c]); integer(r,c)=parse(Int,r[c])
rounded(x)=parse(Float64,@sprintf("%.8g",x))
costs=Dict()
for view in ("fwd","bwd"), path in sort(readdir(ref;join=true))
    occursin(Regex("^score_"*view*"(?:_chunk[1-9][0-9]*)?\\.tsv\\.audit\\.tsv\$"),basename(path)) || continue
    for r in last(read_table(path))
        k=(view,r["file"],integer(r,"lobe"),integer(r,"parity"),integer(r,"mirror"))
        haskey(costs,k) && error("Duplicate input state")
        costs[k]=[num(r,"cost0"),num(r,"cost1")]
    end
end
unchanged=("features.tsv","features_split.tsv","features_local.tsv","features_descriptor.tsv",
    "patches_fwd17.tsv","patches_bwd17.tsv","patches_bwd9.tsv","fisher_cv.tsv","training_support.tsv","pred_kmeans.tsv")
@testset "Complete fixed cohort and exact finite control" begin
    @test length(costs)==8length(base)
    @test !isfile(joinpath(out,"failures.tsv"))
    for r in last(read_table(joinpath(out,"input_hashes.tsv")))
        name=basename(r["input"])
        p=name=="mold_leave_one_out.toml" ? joinpath(out,"decoder_settings.toml") :
          name=="unit_assignment_patch_support.toml" ? joinpath(out,"assignment_config.toml") : joinpath(ref,name)
        @test isfile(p)
        @test bytes2hex(sha256(read(p)))==r["sha256"]
    end
    h,oldfeatures=lobe_table(joinpath(ref,"features_predictor.tsv"))
    _,old=lobe_table(joinpath(ref,"predictions.tsv"))
    for arm in ("reference","loo")
        dir=joinpath(out,arm)
        for name in unchanged
            @test read(joinpath(dir,name))==read(joinpath(ref,name))
        end
        for name in (unchanged...,"score_fwd.tsv","score_bwd.tsv","features_predictor.tsv","pred_gmm.tsv","predictions.tsv")
            @test Set(keys(last(lobe_table(joinpath(dir,name)))))==Set(keys_all)
        end
        _,features=lobe_table(joinpath(dir,"features_predictor.tsv"))
        for k in keys_all,c in setdiff(h,["mold_cc_fwd","mold_cc_bwd"])
            @test features[k][c]==oldfeatures[k][c]
        end
        _,summary=read_table(joinpath(dir,"summary.tsv"))
        @test length(summary)==length(files)
        @test Dict(r["file"]=>integer(r,"N_selected") for r in summary)==
            Dict(f=>count(k->k[1]==f,keys_all) for f in files)
        if arm=="reference"
            for name in ("score_fwd.tsv","score_bwd.tsv","features_predictor.tsv","pred_gmm.tsv","predictions.tsv","summary.tsv")
                @test read(joinpath(dir,name))==read(joinpath(ref,name))
            end
        end
        _,pred=lobe_table(joinpath(dir,"predictions.tsv"))
        println(arm," keys=",length(pred)," scans=",length(files)," unavailable=",count(r->r["predicted"]=="?",values(pred)),
            " decision changes=",count(k->pred[k]["predicted"]!=old[k]["predicted"],keys_all))
    end
end

@testset "Independent target-excluded objectives and serialized margins" begin
    _,ar=read_table(joinpath(out,"loo/state_selection.tsv"))
    audits=Dict((r["view"],r["file"],integer(r,"lobe"),integer(r,"phase"),integer(r,"mirror"))=>r for r in ar)
    @test length(audits)==length(ar)==8length(base)
    for view in ("fwd","bwd")
        h,rows=lobe_table(joinpath(out,"loo","score_"*view*".tsv"))
        @test !any(startswith(c,"global_") || c=="file_cost" for c in h)
        _,old=lobe_table(joinpath(ref,"score_"*view*".tsv"))
        changed_states=0; changed_margins=0; empty_training=0
        for f in files
            ks=filter(k->k[1]==f,keys_all)
            for target in ks
                objectives=Float64[]; counts=Int[]
                for p in 0:1,m in 0:1
                    total=0.0; ntrain=0
                    for other in ks
                        other==target && continue
                        c=costs[(view,f,other[2],mod(other[2]-1+p,2),m)]
                        @test all(isfinite,c) || c==[Inf,Inf]
                        if all(isfinite,c); total+=min(c...); ntrain+=1; end
                    end
                    push!(objectives,total); push!(counts,ntrain)
                    a=audits[(view,target...,p,m)]
                    @test num(a,"training_cost")≈total atol=2e-14
                    @test integer(a,"training_lobes")==ntrain
                end
                @test all(==(first(counts)),counts)
                ntrain=first(counts); r=rows[target]; selected=argmin(objectives)
                p=ntrain==0 ? -1 : div(selected-1,2); m=ntrain==0 ? -1 : mod(selected-1,2)
                @test (integer(r,"state_phase"),integer(r,"state_mirror"))==(p,m)
                @test integer(r,"training_lobes")==ntrain
                @test num(r,"training_cost")≈(ntrain==0 ? 0.0 : objectives[selected]) atol=2e-14
                c=ntrain==0 ? [Inf,Inf] : costs[(view,target...,mod(target[2]-1+p,2),m)]
                reason=ntrain==0 ? "no_remaining_observations" : all(isfinite,c) ? "ok" : "insufficient_native_support"
                @test r["state_reason"]==reason
                for q in 0:1,t in 0:1
                    a=audits[(view,target...,q,t)]
                    @test a["selected"]==string((q,t)==(p,m)) && a["reason"]==reason
                end
                for (column,x) in zip(("cost_GlcN","cost_GlcNAc","cost_margin"),(c...,abs(c[1]-c[2])))
                    @test isequal(num(r,column),rounded(x))
                end
                label=all(isfinite,c) ? Int(c[2]<c[1]) : -1
                @test integer(r,"predicted")==label
                @test r["physical_label"]==(label<0 ? "?" : label==1 ? "GlcNAc" : "GlcN")
                prior=old[target]; phase=mod(integer(prior,"global_phase")+integer(prior,"global_direction")*(length(ks)-1),2)
                changed_states+=((p,m)!=(phase,integer(prior,"global_mirror")))
                changed_margins+=(r["cost_margin"]!=prior["cost_margin"])
                empty_training+=(ntrain==0)
            end
        end
        println(view," changed per-lobe states=",changed_states," changed margins=",changed_margins," no other evidence=",empty_training)
    end
end

#!/usr/bin/env julia
# Saved-output integrity, independent vote arithmetic and refit variability.
# No labels, benchmark membership, fitting or selected-N replacement.
using Test, SHA
include(joinpath(@__DIR__,"run_reconstructed_chitosan.jl"))
length(ARGS)==3 || error("RUN SAVED_SUPPORT SAVED_GEOMETRY required")
run,saved,geometry=ARGS
control=selected_counts(joinpath(run,"original_control_summary.tsv"))
hybrid=selected_counts(joinpath(run,"original_promoted_summary.tsv"))
files=Set(keys(control)); sameN=Set(f for f in files if control[f]==hybrid[f])
tables=("features","features_split","features_local","features_descriptor","features_predictor",
    "patches_fwd17","patches_bwd17","patches_bwd9","training_support","score_fwd","score_bwd",
    "fisher_cv","pred_gmm","pred_kmeans","predictions")
num(r,k)=parse(Float64,r[k])
@testset "Complete selected-count cohorts and invariant settings" begin
    @test !isfile(joinpath(run,"failures.tsv"))
    @test Set(keys(hybrid))==files
    @test control==selected_counts(joinpath(run,"control_counts.tsv"))
    @test hybrid==selected_counts(joinpath(run,"hybrid_counts.tsv"))
    _,changes=read_table(joinpath(run,"count_changes.tsv"))
    @test length(changes)==length(files)
    for r in changes
        @test parse(Int,r["control_N"])==control[r["file"]]
        @test parse(Int,r["promoted_N"])==hybrid[r["file"]]
        @test parse(Int,r["delta"])==hybrid[r["file"]]-control[r["file"]]
    end
    for (arm,counts) in (("control",control),("hybrid",hybrid))
        dir=joinpath(run,arm)
        @test !isfile(joinpath(dir,"failures.tsv"))
        for t in tables
            rows=check_counts(joinpath(dir,t*".tsv"),counts)
            if t in ("features","features_split")
                for (k,r) in rows
                    @test parse(Int,r["N"])==counts[k[1]]
                    @test isfinite(num(r,"gcv")) && num(r,"gcv")>=0
                    @test r["source"] in ("circ","ell")
                    @test all(isfinite(num(r,c)) for c in ("x_nm","y_nm","sigma_parallel_nm","sigma_perp_nm"))
                end
            end
        end
        _,summary=read_table(joinpath(dir,"summary.tsv"))
        @test Dict(r["file"]=>parse(Int,r["N_selected"]) for r in summary)==counts
        @test length(summary)==length(files)
        @test length(readdir(joinpath(dir,"plots","standalone")))==length(files)
        log=read(joinpath(dir,"logs/gmm.log"),String)
        for setting in ("learning family: gaussian","covariance structure: full","covariance scope: final_only",
            "hard assignment: mahalanobis","resampling: none replicates=1 bootstrap_seed=0")
            @test occursin(setting,log)
        end
        names=filter(l->startswith(l,"GMM naming seed "),split(log,'\n'))
        @test length(names)==10
        @test [parse(Int,split(l)[4]) for l in names]==collect(0:9)
        for l in names
            m=match(r"mode=(\w+) naming_means=(\S+) high_cluster=(\d+)",l)
            @test m!==nothing && m[1]=="raw_amplitude"
            means=parse.(Float64,split(m[2],',')); high=parse(Int,m[3])
            @test high in (0,1,2)
            @test high==0 ? !all(isfinite,means) : all(isfinite,means) && means[high]==maximum(means)
        end
    end
end
@testset "Independent vote and missing-component arithmetic" begin
    for arm in ("control","hybrid")
        dir=joinpath(run,arm)
        _,g=lobe_table(joinpath(dir,"pred_gmm.tsv")); _,k=lobe_table(joinpath(dir,"pred_kmeans.tsv"))
        _,p=lobe_table(joinpath(dir,"predictions.tsv"))
        for (key,r) in p
            @test r["model"]=="cc_soft_patch_support_v1"
            if g[key]["predicted"]=="?" || k[key]["predicted"]=="?"
                @test r["predicted"]=="?" && r["probability_1"]=="NA"
            else
                a,b=num(g[key],"probability_1"),num(k[key],"probability_1")
                @test 0<=a<=1 && 0<=b<=1
                v=(a+b)/2
                @test isapprox(num(r,"probability_1"),v;atol=5.1e-9,rtol=0)
                @test isapprox(num(r,"confidence"),2abs(v-.5);atol=5.1e-9,rtol=0)
                @test r["predicted"]==(v>=.5 ? "1" : "0")
            end
        end
        println(arm," unavailable=",sort([key for (key,r) in p if r["predicted"]=="?"]))
        println(arm," prediction_sha256=",bytes2hex(sha256(read(joinpath(dir,"predictions.tsv")))))
    end
end
# Exact replay is NOT presumed for fresh timed fits. Emit measurable differences.
for t in tables
    ref=t in ("features","features_split") ? geometry : saved
    println("control vs saved byte identity ",t,": ",read(joinpath(run,"control",t*".tsv"))==read(joinpath(ref,t*".tsv")))
end
for (a,b,label,subset) in ((geometry,joinpath(run,"control"),"saved -> control",files),
    (joinpath(run,"control"),joinpath(run,"hybrid"),"unchanged N: control -> hybrid",sameN))
    for t in ("features","features_split")
        _,ra=lobe_table(joinpath(a,t*".tsv")); _,rb=lobe_table(joinpath(b,t*".tsv"))
        ks=sort([k for k in intersect(keys(ra),keys(rb)) if k[1] in subset])
        changed=[k for k in ks if ra[k]!=rb[k]]
        shift(k)=hypot(num(ra[k],"x_nm")-num(rb[k],"x_nm"),num(ra[k],"y_nm")-num(rb[k],"y_nm"))
        relative(k)=abs(num(rb[k],"gcv")/num(ra[k],"gcv")-1)
        println(label," ",t,": rows=",length(ks)," changed_rows=",length(changed),
            " changed_scans=",length(unique(first.(changed)))," max_position_nm=",maximum(shift,ks;init=0.),
            " max_relative_GCV=",maximum(relative,ks;init=0.))
    end
end
println("Counts control/hybrid: ",sum(values(control))," / ",sum(values(hybrid)),
    "; changed scans ",length(files)-length(sameN))

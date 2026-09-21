using Test, Random, Statistics, LinearAlgebra, Printf, TOML
include(joinpath(@__DIR__, "lib", "empirical_fisher_native.jl"))
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .EmpiricalFisherNative, .ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_complete_training.toml")
const BASE = load_fisher_config(CONTROL)
const COMPLETE = load_fisher_config(CANDIDATE)
const GRID = fisher_grid(BASE)

function support_fixture(path, rowkeys, mask; family="fwd17")
    rows = [Dict("file"=>key[1], "lobe"=>string(key[2]), "fwd17_observed"=>"289",
                 "bwd17_observed"=>"289", "bwd9_observed"=>"81") for key in rowkeys]
    for i in findall(.!mask)
        rows[i][family*"_observed"] = family == "bwd9" ? "80" : "288"
    end
    write_table(path, ["file","lobe","fwd17_observed","bwd17_observed","bwd9_observed"], rows)
    return rows
end

@testset "Explicit policy changes no other parameter" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    @test BASE.training_support == "all_admissible"
    @test COMPLETE.training_support == "complete_patches"
    a,b = load_config(CONTROL),load_config(CANDIDATE)
    b["selection"]["assignment_training_support"] = a["selection"]["assignment_training_support"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    for value in (nothing,"guess",true,1)
        c = deepcopy(a)
        value === nothing ? delete!(c["selection"],"assignment_training_support") :
            (c["selection"]["assignment_training_support"] = value)
        @test_throws ArgumentError load_training_policy(c)
        @test_throws ArgumentError load_fisher_config(c)
    end
    @test validate_training_mask("all_admissible",nothing,3) == trues(3)
    @test_throws ArgumentError validate_training_mask("all_admissible",trues(3),3)
    for mask in (nothing,[1,0,1],trues(2))
        @test_throws ArgumentError validate_training_mask("complete_patches",mask,3)
    end
end

@testset "Observed pixels determine eligibility with exact, order-independent keys" begin
    mktempdir() do dir
        paths = String[]
        rowkeys = [("synthetic.sxm",i) for i in 1:4]
        for (family,prefix,side) in ReconstructedUnitAssignment.TRAINING_PATCH_FAMILIES
            rows = Dict{String,String}[]
            columns = [prefix*lpad(string(i),3,'0') for i in 1:side^2]
            for i in 1:4
                row = Dict("file"=>rowkeys[i][1], "lobe"=>string(i))
                for col in columns; row[col]="0.0"; end
                i == 2 && family == "fwd17" && (row[first(columns)]="NaN")
                i == 3 && family == "bwd17" && (row[last(columns)]="NA")
                i == 4 && family == "bwd9" && (row[first(columns)]="missing")
                push!(rows,row)
            end
            path=joinpath(dir,family*".tsv")
            write_table(path,vcat(["file","lobe"],columns),reverse(rows)); push!(paths,path)
        end
        output=joinpath(dir,"support.tsv")
        write_training_support(paths...,output)
        @test load_training_mask("complete_patches",output,rowkeys,"fisher") == [true,false,true,true]
        @test load_training_mask("complete_patches",output,rowkeys,"gmm") == [true,false,false,false]
        @test load_training_mask("complete_patches",output,reverse(rowkeys),"gmm") == [false,false,false,true]
        @test load_training_mask("all_admissible","",rowkeys,"gmm") === nothing
        @test_throws ArgumentError load_training_mask("all_admissible",output,rowkeys,"gmm")
        @test_throws ArgumentError load_training_mask("complete_patches","",rowkeys,"gmm")
        @test_throws ErrorException load_training_mask("complete_patches",output,rowkeys[1:3],"gmm")
        @test_throws ErrorException load_training_mask("complete_patches",output,vcat(rowkeys,[rowkeys[1]]),"gmm")
        @test_throws ErrorException write_training_support(paths...,output)
        header,rows=read_table(output)
        for (name,value) in (("fwd17_observed","290"),("bwd9_observed","-1"),("bwd17_observed","NaN"))
            bad=deepcopy(rows); bad[1][name]=value
            path=joinpath(dir,name*"_bad.tsv"); write_table(path,header,bad)
            @test_throws ErrorException load_training_mask("complete_patches",path,rowkeys,"gmm")
        end
        path=joinpath(dir,"forbidden.tsv")
        write_table(path,vcat(header,["expected_N"]),rows)
        @test_throws ErrorException load_training_mask("complete_patches",path,rowkeys,"gmm")
    end
end

@testset "Fisher uses complete opposite folds only, still scores admissible partials" begin
    rng=MersenneTwister(27); n=64
    full=0.08randn(rng,n,GRID.side^2)
    pattern=[exp(-((r-9)^2+(c-9)^2)/18)+0.03(c-9) for r in 1:GRID.side for c in 1:GRID.side]
    for i in 1:n; full[i,:] .+= (i <= 44 ? -3.0 : 3.0).*pattern; end
    keys=[("synthetic.sxm",i) for i in 1:n]
    table=PatchTable(keys,full[:,GRID.disk_indices],full[:,GRID.center_index],fill("",n),GRID)
    a=cv_scores(table,BASE)
    b=cv_scores(table,COMPLETE;training_mask=trues(n))
    @test [s.score for s in a] == [s.score for s in b]
    @test_throws ArgumentError cv_scores(table,COMPLETE)
    mask=trues(n); mask[[10,21,30,47]].=false
    scores=cv_scores(table,COMPLETE;training_mask=mask)
    @test all(s->isfinite(s.score),scores) && length(scores)==n
    for parity in (0,1)
        train=findall(i->mask[i] && mod(i,2)==parity,1:n)
        model=fit_fisher(table.X[train,:],table.amplitudes[train],COMPLETE)
        for i in findall(i->mod(i,2)!=parity,1:n)
            @test scores[i].score == maxmirror_score(table.X[i,:],model,GRID)
        end
    end
    changedX,changedA=copy(table.X),copy(table.amplitudes)
    changedX[.!mask,:] .*= -1000; changedA[.!mask] .= 1e9
    changed=PatchTable(keys,changedX,changedA,table.invalid_reasons,GRID)
    other=cv_scores(changed,COMPLETE;training_mask=mask)
    @test [s.score for s in scores[mask]] == [s.score for s in other[mask]]
    order=randperm(rng,n)
    permuted=PatchTable(keys[order],table.X[order,:],table.amplitudes[order],table.invalid_reasons[order],GRID)
    # Row permutations can change floating-point SVD signs/rounding; compare tightly.
    @test [s.score for s in cv_scores(permuted,COMPLETE;training_mask=mask[order])] ≈
          [s.score for s in scores[order]] rtol=1e-9 atol=1e-9
    noeven=mask .& isodd.(1:n)
    missing=cv_scores(table,COMPLETE;training_mask=noeven)
    @test all(s->isnan(s.score) && startswith(s.invalid_reason,"training_even:"),missing[1:2:n])
    @test all(s->isfinite(s.score),missing[2:2:n])
    badX=copy(table.X); badX[10,1]=NaN
    bad=cv_scores(PatchTable(keys,badX,table.amplitudes,table.invalid_reasons,GRID),COMPLETE;training_mask=mask)
    @test isnan(bad[10].score) && bad[10].invalid_reason=="nonfinite_patch_input"
    mktempdir() do dir
        path=joinpath(dir,"patch.tsv")
        open(path,"w") do io
            println(io,join(vcat(["file","lobe"],[@sprintf("res_p%03d",i) for i in 1:289]),'\t'))
            for i in 1:n
                row=string.(full[i,:]); !mask[i] && (row[1]="NaN")
                println(io,join(vcat([keys[i][1],string(i)],row),'\t'))
            end
        end
        support=joinpath(dir,"support.tsv"); support_fixture(support,keys,mask)
        out=joinpath(dir,"scores.tsv"); expected=joinpath(dir,"expected.tsv")
        loaded=load_patches(path,"res",COMPLETE)
        write_scores(expected,cv_scores(loaded,COMPLETE;training_mask=mask))
        args=["--patches",path,"--prefix","res","--config",CANDIDATE,"--training-support",support,"--out",out]
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_empirical_fisher_native.jl")) $args`)
        @test read(out)==read(expected)
    end
end

@testset "GMM freezes normalization, fitting, free weights and physical naming before prediction" begin
    mktempdir() do dir
        rng=MersenneTwister(971); names=["f$i" for i in 1:4]
        records=GMM.LobeRecord[]
        for i in 1:80
            high=mod1(i,10)>7
            amp=(high ? 0.3 : 0.1)+0.005randn(rng)
            vals=0.2randn(rng,4).+(high ? 2.0 : -2.0); vals[1]=amp
            push!(records,GMM.LobeRecord("synthetic_$(cld(i,10)).sxm",mod1(i,10),amp,Dict(zip(names,vals))))
        end
        keys=[(r.file,r.lobe) for r in records]
        mask=trues(80); mask[10:10:80].=false
        path=joinpath(dir,"features.tsv")
        rows=[Dict("file"=>r.file,"lobe"=>string(r.lobe),"amplitude"=>string(r.amplitude),
                   (f=>string(r.features[f]) for f in names)...) for r in records]
        write_table(path,vcat(["file","lobe","amplitude"],names),rows)
        support=joinpath(dir,"support.tsv"); support_fixture(support,keys,mask;family="bwd9")
        args=["--features",path,"--seeds","2","--selftrain","2","--interactions","--view","v="*join(names,',')]
        old=GMM._parse_cli(vcat(args,["--config",CONTROL]))
        opt=GMM._parse_cli(vcat(args,["--config",CANDIDATE,"--training-support",support]))
        @test_throws ErrorException GMM._parse_cli(vcat(args,["--config",CANDIDATE]))
        @test_throws ErrorException GMM._parse_cli(vcat(args,["--config",CONTROL,"--training-support",support]))
        @test_throws ArgumentError GMM._view_probability(records,names,opt)
        @test isequal(GMM._view_probability(records,names,old),
                      GMM._view_probability(records,names,opt;training_mask=trues(80)))
        Z,valid=GMM._standardized_matrix(records,names;training_mask=mask)
        @test all(valid)
        for file in unique(r.file for r in records), j in 1:4
            idxs=findall(r->r.file==file,records); train=filter(i->mask[i],idxs)
            raw=[records[i].features[names[j]] for i in train]
            @test Z[idxs,j] == [(records[i].features[names[j]]-mean(raw))/std(raw) for i in idxs]
        end
        da,db=[],[]
        a=GMM._view_probability(records,names,opt;training_mask=mask,diagnostics=da)
        changed=deepcopy(records)
        for i in findall(.!mask)
            changed[i].amplitude=-1e12
            for f in names; changed[i].features[f] *= -100; end
        end
        b=GMM._view_probability(changed,names,opt;training_mask=mask,diagnostics=db)
        @test length(a)==80 && all(isfinite,a)
        @test a[mask]==b[mask] # Includes frozen high-amplitude group naming.
        @test da==db && length(da)==2
        @test all(d->d.indices==findall(mask),da)
        @test all(d->length(d.clusters)==2 && d.clusters[1].weight != d.clusters[2].weight,da)
        z,_=GMM._standardized_matrix(changed,names;training_mask=mask)
        @test z[mask,:]==Z[mask,:]
        extra=GMM.LobeRecord("new_partial_file.sxm",1,1e99,Dict(f=>1e10 for f in names))
        extended=vcat(records,[extra]); dc=[]
        c=GMM._view_probability(extended,names,opt;training_mask=vcat(mask,false),diagnostics=dc)
        @test isequal(c[1:80],a) && isnan(c[81]) && dc==da
        nofile=copy(mask); nofile[1:10].=false
        _,finite=GMM._standardized_matrix(records,names;training_mask=nofile)
        @test !any(finite[1:10]) && all(finite[11:80])
        @test all(isnan,GMM._view_probability(records,names,opt;training_mask=falses(80)))
        invalid=deepcopy(records); invalid[10].features["f2"]=NaN
        d=GMM._view_probability(invalid,names,opt;training_mask=mask)
        @test isnan(d[10]) && isequal(d[mask],a[mask])
        expected=joinpath(dir,"expected.tsv"); out=joinpath(dir,"out.tsv")
        GMM._write_predictions(expected,records,a,Int.(isfinite.(a)))
        cli=vcat(args,["--config",CANDIDATE,"--training-support",support,"--out",out])
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $cli`)
        @test read(out)==read(expected)
        renamed=deepcopy(records)
        for r in renamed; r.file="renamed_"*r.file; end
        @test isequal(a,GMM._view_probability(renamed,names,opt;training_mask=mask))
    end
end

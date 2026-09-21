using Test, Statistics, Random, LinearAlgebra, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_scan_bagging.toml")

@testset "One explicit whole-scan bootstrap; no combined variant" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_gmm_resampling(a) == (mode="none", replicates=1, seed=0)
    @test load_gmm_resampling(b) == (mode="whole_scans", replicates=20, seed=0)
    b["model"]["name"] = a["model"]["name"]
    for field in ("gmm_resampling", "gmm_bootstrap_replicates")
        b["selection"][field] = a["selection"][field]
    end
    @test a == b
    @test_throws ErrorException ReconstructedRepresentationDiagnostics.load_inputs(Dict(), CANDIDATE)
    for name in readdir(joinpath(ROOT,"config"))
        startswith(name,"unit_assignment_") && endswith(name,".toml") || continue
        cfg = TOML.parsefile(joinpath(ROOT,"config",name))
        haskey(get(cfg,"selection",Dict()), "gmm_seed_aggregation") || continue
        @test load_gmm_resampling(cfg).mode == (name == "unit_assignment_scan_bagging.toml" ? "whole_scans" : "none")
    end
    for (field, values) in (("gmm_resampling", (nothing,"lobe","class_stratified",true,1)),
                          ("gmm_bootstrap_replicates", (nothing,0,-1,true,2.5)),
                          ("gmm_bootstrap_seed", (nothing,-1,true,0.5,typemax(Int))))
        for value in values
            cfg = load_config(CANDIDATE)
            value === nothing ? delete!(cfg["selection"],field) : (cfg["selection"][field]=value)
            @test_throws ArgumentError load_gmm_resampling(cfg)
        end
    end
    bad = deepcopy(a); bad["selection"]["gmm_bootstrap_replicates"] = 2
    @test_throws ArgumentError load_gmm_resampling(bad)
    for (section,key,value) in (("selection","gmm_training_weighting","equal_scans"),
        ("selection","gmm_seed_aggregation","mean_membership"),
        ("selection","assignment_training_support","complete_patches"),
        ("preprocessing","gmm_feature_normalization","median_iqr"),
        ("model","gmm_final_covariance","ledoit_wolf"))
        cfg=load_config(CANDIDATE); cfg[section][key]=value
        @test_throws ArgumentError load_gmm_resampling(cfg)
    end
end

function fixture_records()
    rng = MersenneTwister(922)
    records = GMM.LobeRecord[]
    for file in 1:6, lobe in 1:(7+file)
        high = lobe > .65(7+file)
        push!(records, GMM.LobeRecord("synthetic_$file.sxm", lobe,
            1 + 2high + .1randn(rng),
            Dict("f1"=>(high ? 1. : -.5) + randn(rng) + file,
                 "f2"=>(high ? .4 : -.2) + randn(rng) - file)))
    end
    records[2].features["f2"] = NaN
    push!(records,GMM.LobeRecord("unavailable.sxm",1,1.,Dict("f1"=>NaN,"f2"=>NaN)))
    records
end

@testset "Literal whole-scan draws, multiplicities and identity-only grouping" begin
    records=fixture_records()
    ids=setdiff(1:length(records),[2,length(records)])
    files=unique(records[i].file for i in ids)
    groups=[[i for i in ids if records[i].file==f] for f in files]
    for seed in 0:19
        expected=reduce(vcat,groups[rand(MersenneTwister(seed),1:length(groups),length(groups))])
        drawn=GMM._scan_bootstrap_sample(records,ids,seed)
        @test drawn==expected
        @test all(i->i in ids,drawn)
        @test sum(count(==(first(g)),drawn) for g in groups)==length(groups)
        @test all(length(unique(count(==(i),drawn) for i in g))==1 for g in groups)
    end
    renamed=deepcopy(records)
    for (i,r) in enumerate(renamed)
        r.file="opaque_"*r.file; r.lobe+=100; r.features["N"]=-999
    end
    @test GMM._scan_bootstrap_sample(records,ids,17)==GMM._scan_bootstrap_sample(renamed,ids,17)
    @test isempty(GMM._scan_bootstrap_sample(records,Int[],0))
    @test GMM._scan_bootstrap_sample(records,[1,3,4],0)==[1,3,4]
end

@testset "Independent bag arithmetic, training-only naming and unchanged output rules" begin
    records=fixture_records()
    mktempdir() do dir
        features=joinpath(dir,"features.tsv")
        write_table(features,["file","lobe","amplitude","f1","f2"],
            [Dict("file"=>r.file,"lobe"=>r.lobe,"amplitude"=>r.amplitude,
                  "f1"=>r.features["f1"],"f2"=>r.features["f2"]) for r in records])
        cfg=load_config(CANDIDATE); cfg["selection"]["gmm_bootstrap_replicates"]=4
        config=joinpath(dir,"synthetic.toml")
        open(io->TOML.print(io,cfg),config,"w")
        auditpath=joinpath(dir,"audit.tsv")
        args=["--features",features,"--config",config,"--view","toy=f1,f2",
              "--seeds","3","--selftrain","2","--interactions","--bootstrap-audit",auditpath]
        opt=GMM._parse_cli(args)
        @test opt.resampling=="whole_scans" && opt.bootstrap_replicates==4
        diagnostics=[]; audit=Dict{String,String}[]
        actual=GMM._view_probability(records,["f1","f2"],opt;diagnostics,audit,view_name="toy")
        Z,valid=GMM._standardized_matrix(records,["f1","f2"];interactions=true,
            normalization="mean_sample_std",scale_fallback=1.)
        ids=findall(valid); expected=zeros(length(ids)); accepted_bags=0
        @test length(diagnostics)==12
        @test length(audit)==4*7
        for bag in 1:4
            sums=zeros(length(ids)); accepted=0
            ds=filter(d->d.replicate==bag,diagnostics)
            idx=GMM._scan_bootstrap_sample(records,ids,bag-1)
            for d in ds
                @test d.indices==idx && d.bootstrap_seed==bag-1
                fit=GMM._gmm_fit(permutedims(Z[idx,:]),d.seed;ridge=opt.covariance_ridge)
                fit=GMM._mahalanobis_self_train(permutedims(Z[idx,:]),fit...;
                    iters=2,covariance_mode="ridge",ridge=opt.covariance_ridge)
                @test fit==(d.means,d.covariances,d.weights)
                @test sum(d.weights)≈1
                labels=Int[]
                for i in ids
                    scores=[log(d.weights[c])-sum(abs2,
                        cholesky(Symmetric(d.covariances[c])+1e-8I).L \ (Z[i,:]-d.means[:,c]))/2 for c in 1:2]
                    push!(labels,argmax(scores))
                end
                byrow=Dict(zip(ids,labels))
                amps=Dict(c=>[records[i].amplitude for i in idx if byrow[i]==c] for c in 1:2)
                means=Dict(c=>mean(v) for (c,v) in amps if !isempty(v))
                @test means==d.amplitude_means
                length(means)==2 || continue
                high=argmax([means[1],means[2]])
                @test high==d.high_cluster
                sums .+= labels.==high; accepted+=1
            end
            rows=filter(r->r["replicate"]==string(bag),audit)
            @test sum(parse(Int,r["multiplicity"]) for r in rows)==6
            @test sum(parse(Int,r["training_rows"]) for r in rows)==length(idx)
            @test all(parse(Int,r["training_rows"])==parse(Int,r["multiplicity"])*parse(Int,r["usable_rows"]) for r in rows)
            @test all(parse(Int,r["valid_seeds"])==accepted for r in rows)
            accepted>0 && (expected .+= sums/accepted; accepted_bags+=1)
        end
        @test actual[ids]==expected/accepted_bags
        @test all(isnan,actual[.!valid]) && count(.!valid)==2
        @test all(p->0<=p<=1,actual[ids])
        @test isequal(actual,GMM._view_probability(records,["f1","f2"],opt))
        renamed=deepcopy(records)
        for r in renamed
            r.file="opaque_"*r.file; r.lobe+=100; r.features["N"]=-999
        end
        @test isequal(actual,GMM._view_probability(renamed,["f1","f2"],opt))
        # An out-of-bag scan remains prediction-eligible but cannot rename groups.
        one=GMM.Options((f==:bootstrap_replicates ? 1 : getfield(opt,f) for f in fieldnames(GMM.Options))...)
        sampled=Set(records[i].file for i in first(diagnostics).indices)
        @test length(sampled)<6
        perturbed=deepcopy(records)
        for r in perturbed
            r.file in sampled || (r.amplitude=1e12)
        end
        @test isequal(GMM._view_probability(records,["f1","f2"],one),
                      GMM._view_probability(perturbed,["f1","f2"],one))
        @test any(i->!(records[i].file in sampled) && isfinite(actual[i]),ids)
        output,reference=joinpath.(dir,("out.tsv","reference.tsv"))
        GMM._write_predictions(reference,records,actual,Int.(valid))
        cli=vcat(args,["--out",output])
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $cli`)
        @test read(output)==read(reference)
        _,saved_audit=read_table(auditpath)
        @test saved_audit==audit
        @test_throws ErrorException GMM._parse_cli(cli) # Never overwrite the audit.
        for flag in ("--truth","--sequence","--expected-N","--benchmark-manifest")
            @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
        end
        @test_throws ErrorException GMM._parse_cli(args[1:end-2])
        @test_throws ErrorException GMM._parse_cli(["--features",features,"--config",CONTROL,"--bootstrap-audit",joinpath(dir,"bad.tsv")])
    end
end

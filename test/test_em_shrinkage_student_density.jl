using Test, Random, LinearAlgebra, Statistics, TOML
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const RU = GMM.ReconstructedUnitAssignment
const AC = GMM.AssignmentCovariance
const AM = GMM.AssignmentMixtures
const ROOT = dirname(@__DIR__)
config(name) = joinpath(ROOT, "config", "unit_assignment_" * name * ".toml")

@testset "Explicit, isolated covariance scope and Student decisions" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    control = RU.load_config(config("patch_support"))
    shrink = RU.load_config(config("em_shrinkage"))
    student = RU.load_config(config("student_density"))
    for (cfg, changes) in ((shrink, ["gmm_covariance_scope", "gmm_final_covariance"]),
                          (student, ["gmm_learning_family", "gmm_hard_assignment", "gmm_final_score"]))
        matched = deepcopy(cfg)
        for key in vcat(changes, ["name"])
            matched["model"][key] = control["model"][key]
        end
        @test matched == control
    end
    @test RU.load_gmm_covariance_scope(shrink) == "all_updates"
    @test RU.load_gmm_learning(student).hard_assignment == "student_density"
    for key in ("gmm_covariance_scope", "gmm_hard_assignment"), bad in (nothing, "auto", true, 2)
        cfg = deepcopy(control)
        bad === nothing ? delete!(cfg["model"], key) : (cfg["model"][key] = bad)
        @test_throws ArgumentError (key == "gmm_covariance_scope" ?
            RU.load_gmm_covariance_scope(cfg) : RU.load_gmm_learning(cfg))
    end
    for (section, key, bad) in (("model", "gmm_learning_family", "student_t"),
            ("model", "gmm_final_covariance", "ridge"), ("model", "gmm_covariance_structure", "tied"),
            ("model", "gmm_final_score", "gaussian_density"), ("model", "gmm_hard_assignment", "student_density"),
            ("selection", "fisher_cv_scheme", "scan_hash_twofold"), ("selection", "gmm_cluster_naming", "within_scan_z"),
            ("selection", "gmm_training_weighting", "equal_scans"), ("selection", "gmm_seed_aggregation", "mean_membership"),
            ("selection", "gmm_resampling", "whole_scans"), ("selection", "assignment_training_support", "complete_patches"),
            ("preprocessing", "gmm_feature_normalization", "median_iqr"))
        cfg = deepcopy(shrink); cfg[section][key] = bad
        @test_throws ArgumentError RU.load_gmm_covariance_scope(cfg)
    end
    for (key, bad) in (("gmm_learning_family", "gaussian"), ("gmm_learning_family", "factor_analyzer"),
            ("gmm_final_score", "mahalanobis"), ("gmm_hard_assignment", "mahalanobis"),
            ("gmm_covariance_scope", "all_updates"))
        cfg = deepcopy(student); cfg["model"][key] = bad
        @test_throws ArgumentError RU.load_gmm_learning(cfg)
    end
    for name in ("reconstructed", "patch_support", "student_t", "factor_analyzer", "shrunk_gmm")
        cfg = RU.load_config(config(name))
        @test RU.load_gmm_covariance_scope(cfg) == "final_only"
        @test RU.load_gmm_learning(cfg).hard_assignment == "mahalanobis"
    end
end

# Independent dense outer-product calculation, deliberately not the fast formula.
function reference_covariance(Y, r; ridge=1e-6)
    p, n = size(Y); a = r / sum(r)
    outer = [Y[:, i] * Y[:, i]' for i in 1:n]
    S = sum(a[i] * outer[i] for i in 1:n)
    T = tr(S) / p * Matrix{Float64}(I, p, p)
    denominator = sum(abs2, S-T)
    numerator = sum(a[i]^2 * sum(abs2, outer[i]-S) for i in 1:n)
    alpha = denominator > 0 ? clamp(numerator/denominator, 0., 1.) : 0.
    return (covariance=(1-alpha)*S+alpha*T+ridge*I, sample=S, shrinkage=alpha)
end

@testset "Responsibility shrinkage algebra and hard-membership limit" begin
    rng = MersenneTwister(2026092201)
    for (p,n) in ((1,12), (4,30), (15,7)), seed in 1:3
        X = randn(rng,p,n); r = rand(rng,n); r[1]=0
        Y = X .- X*r/sum(r)
        got = AC.responsibility_covariance(Y,r;ridge=1e-6)
        expected = reference_covariance(Y,r)
        @test got.sample ≈ expected.sample
        @test got.covariance ≈ expected.covariance
        @test got.shrinkage ≈ expected.shrinkage atol=1e-12
        @test 0 <= got.shrinkage <= 1 && isposdef(Symmetric(got.covariance))
        @test AC.responsibility_covariance(Y,7r;ridge=1e-6).covariance ≈ got.covariance
        order = randperm(rng,n)
        @test AC.responsibility_covariance(Y[:,order],r[order];ridge=1e-6).covariance ≈ got.covariance
        centered = X .- mean(X;dims=2)
        hard = AC.final_covariance(centered;mode="ledoit_wolf",ridge=1e-6)
        weighted = AC.responsibility_covariance(centered,ones(n);ridge=1e-6)
        @test weighted.covariance ≈ hard.covariance
        @test weighted.shrinkage ≈ hard.shrinkage atol=1e-12
        binary = Float64.(isodd.(1:n)); ids = findall(==(1),binary)
        centered = X .- mean(X[:,ids];dims=2)
        @test AC.responsibility_covariance(centered,binary;ridge=1e-6).covariance ≈
            AC.final_covariance(centered[:,ids];mode="ledoit_wolf",ridge=1e-6).covariance
        Q = Matrix(qr(randn(rng,p,p)).Q)
        @test AC.responsibility_covariance(Q*Y,r;ridge=1e-6).covariance ≈ Q*got.covariance*Q'
    end
    for X in (zeros(3,1), zeros(3,5), [1. -1. 0. 0.;0. 0. 1. -1.])
        got=AC.responsibility_covariance(X,ones(size(X,2));ridge=1e-6)
        @test got.shrinkage == 0
        @test isposdef(Symmetric(got.covariance))
    end
    for (X,r) in ((zeros(2,0),Float64[]), (zeros(0,2),ones(2)),
            (zeros(2,3),ones(2)), (zeros(2,3),zeros(3)), (zeros(2,3),[1.,-1.,1.]),
            (zeros(2,3),[NaN,1.,1.]), (fill(Inf,2,3),ones(3)))
        @test_throws ArgumentError AC.responsibility_covariance(X,r;ridge=1e-6)
    end
    for ridge in (0.,-1.,Inf,NaN,true)
        @test_throws ArgumentError AC.responsibility_covariance(zeros(2,3),ones(3);ridge)
    end
end

@testset "Shrinkage at initialization, every EM M-step and both hard updates" begin
    rng=MersenneTwister(92201)
    X=hcat(randn(rng,5,63).-1,randn(rng,5,27).+2)
    X[1,:] .+= .8X[2,:]
    for seed in 0:2
        trace=[]
        mu,C,w=GMM._gmm_fit(X,seed;ridge=1e-6,max_iter=4,tol=1e-30,covariance_scope="all_updates",trace)
        GMM._mahalanobis_self_train(X,mu,C,w;iters=2,covariance_mode="ledoit_wolf",ridge=1e-6,
            covariance_scope="all_updates",trace)
        @test [s.stage for s in trace] == ["initial", "em", "em", "em", "em", "hard", "hard"]
        previous=nothing
        for step in trace
            if previous !== nothing
                d=[dot(X[:,i]-previous.means[:,c],inv(previous.covs[c]+1e-8I)*(X[:,i]-previous.means[:,c]))
                    for i in axes(X,2),c in 1:2]
                if step.stage == "em"
                    score=[log(previous.weights[c])-log(det(previous.covs[c]+1e-8I))/2-d[i,c]/2
                        for i in axes(X,2),c in 1:2]
                    prob=exp.(score.-maximum(score;dims=2)); prob ./= sum(prob;dims=2)
                    @test step.responsibilities ≈ prob atol=1e-11
                else
                    a=[argmin(d[i,:]) for i in axes(X,2)]
                    @test step.responsibilities == Float64.([a[i]==c for i in axes(X,2),c in 1:2])
                end
            end
            for c in 1:2
                r=step.responsibilities[:,c]; mean=X*r/sum(r)
                expected=reference_covariance(X.-mean,r)
                @test step.means[:,c] ≈ mean
                @test step.covs[c] ≈ expected.covariance
                @test step.shrinkages[c] ≈ expected.shrinkage atol=1e-11
                @test isposdef(Symmetric(step.covs[c]))
            end
            @test sum(step.weights) ≈ 1
            @test step.covs[1] != step.covs[2]
            if step.stage != "initial"
                @test step.weights ≈ vec(sum(step.responsibilities;dims=1))/size(X,2)
            end
            previous=step
        end
        @test trace[1].weights == [.5,.5] # Unchanged initializer, not a final composition quota.
    end
    @test_throws ArgumentError GMM._gmm_fit(X,0;ridge=1e-6,covariance_scope="all_updates",covariance_structure="tied")
    @test_throws ArgumentError GMM._gmm_fit(X,0;ridge=1e-6,covariance_scope="all_updates",observation_weights=ones(90))
end

@testset "Student density agrees between E-step, hard assignment and final scoring" begin
    rng=MersenneTwister(92202)
    X=hcat(randn(rng,5,63).-1,randn(rng,5,27).+2)
    s=merge(RU.load_gmm_learning(RU.load_config(config("student_density"))),(maxiter=4,tolerance=1e-30))
    for seed in 0:2
        initial=GMM._gmm_fit(X,seed;ridge=1e-6,max_iter=0)
        state=AM.initialize_state(initial...,s;ridge=1e-6)
        trace=[]; fit=AM.fit_mixture(X,initial,s;ridge=1e-6,hard_iterations=2,trace)
        @test length(trace)==6 && !fit.converged && fit.updates==4
        for step in trace
            e=AM.expectation(X,state,s)
            scores=[AM.student_component_score(X[:,i],state.means[:,c],state.covs[c],state.weights[c];
                df=s.student_df,guard=s.cholesky_guard) for i in axes(X,2),c in 1:2]
            densities=[state.weights[c]/sqrt(det(state.covs[c]+s.cholesky_guard*I)) *
                (1+dot(X[:,i]-state.means[:,c],inv(state.covs[c]+s.cholesky_guard*I)*(X[:,i]-state.means[:,c]))/s.student_df)^(-(s.student_df+5)/2)
                for i in axes(X,2),c in 1:2]
            @test exp.(scores) ≈ densities
            @test e.responsibilities ≈ densities ./ sum(densities;dims=2)
            r=e.responsibilities
            if step.stage == "hard"
                a=[argmax(densities[i,:]) for i in axes(X,2)]
                r=Float64.([a[i]==c for i in axes(X,2),c in 1:2])
            end
            for c in 1:2
                u=(s.student_df+5)./(s.student_df.+e.distances[:,c]); ru=r[:,c].*u
                mean=sum(ru[i]*X[:,i] for i in axes(X,2))/sum(ru)
                C=sum(ru[i]*(X[:,i]-mean)*(X[:,i]-mean)' for i in axes(X,2))/sum(r[:,c])+1e-6I
                @test step.state.means[:,c] ≈ mean
                @test step.state.covs[c] ≈ C
            end
            @test step.state.weights ≈ vec(sum(r;dims=1))/size(X,2)
            state=step.state
        end
    end
    # Deliberately separate distance-only assignment from a mass/volume-aware score.
    scores=[AM.student_component_score([0.,0.],m,C,w;df=5.,guard=1e-8)
        for (m,C,w) in (([0.,0.],100Matrix{Float64}(I,2,2),.1),([1.,0.],Matrix{Float64}(I,2,2),.9))]
    @test argmax(scores)==2 # Pure distance would choose the first mean.
end

@testset "CLI determinism, final vote reconstruction and label-free row retention" begin
    rng=MersenneTwister(92203); records=GMM.LobeRecord[]; features=["f$i" for i in 1:5]
    for f in 1:6, lobe in 1:(8+f)
        high=lobe>.7(8+f)
        push!(records,GMM.LobeRecord("synthetic_$f.sxm",lobe,1+2high+.1randn(rng),
            Dict(name=>3high+.5randn(rng)+f for name in features)))
    end
    records[2].features["f2"]=NaN
    push!(records,GMM.LobeRecord("unavailable.sxm",1,1.,Dict(name=>NaN for name in features)))
    mktempdir() do dir
        file=joinpath(dir,"features.tsv")
        RU.write_table(file,vcat(["file","lobe","amplitude"],features),
            [merge(Dict{String,Any}(r.features),Dict("file"=>r.file,"lobe"=>r.lobe,"amplitude"=>r.amplitude)) for r in records])
        for variant in ("em_shrinkage","student_density","patch_support")
            args=["--features",file,"--config",config(variant),"--view","test="*join(features,','),
                "--seeds","3","--selftrain","2","--interactions"]
            opt=GMM._parse_cli(args); trace=[]
            probs=GMM._view_probability(records,features,opt;diagnostics=trace)
            @test isequal(probs,GMM._view_probability(records,features,opt))
            Z,valid=GMM._standardized_matrix(records,features;interactions=true,
                normalization="mean_sample_std",scale_fallback=1.)
            ids=findall(valid); votes=zeros(length(ids))
            @test length(trace)==3
            for d in trace
                if variant == "student_density"
                    state=last(d.clusters).state; mu,C,w=state.means,state.covs,state.weights
                else
                    mu=hcat([c.mean for c in d.clusters]...); C=[c.covariance for c in d.clusters]; w=[c.weight for c in d.clusters]
                end
                scores=[GMM._final_component_score(Z[i,:],mu[:,c],C[c],w[c],opt) for i in ids,c in 1:2]
                a=[argmax(scores[j,:]) for j in eachindex(ids)]
                high=argmax([mean(records[ids[j]].amplitude for j in eachindex(ids) if a[j]==c) for c in 1:2])
                votes .+= a.==high
                if variant == "student_density"
                    state=(means=mu,covs=C,weights=w)
                    @test a == [argmax(r) for r in eachrow(AM.expectation(permutedims(Z[ids,:]),state,opt.learning).responsibilities)]
                elseif variant == "em_shrinkage"
                    @test count(t->t.stage=="hard",d.updates)==2
                    @test first(d.updates).stage=="initial"
                end
            end
            @test probs[ids]==votes/3
            @test count(.!valid)==2 && all(isnan,probs[.!valid])
            renamed=deepcopy(records)
            for r in renamed
                r.file="renamed_"*r.file; r.lobe+=100
                r.features["N"]=-999; r.features["expected_N"]=1e20
            end
            @test isequal(probs,GMM._view_probability(renamed,features,opt))
            expected,out=joinpath.(dir,(variant*"_expected.tsv",variant*"_cli.tsv"))
            GMM._write_predictions(expected,records,probs,Int.(valid))
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args --out $out`)
            @test read(expected)==read(out)
            @test length(last(RU.lobe_table(out)))==length(records)
            for flag in ("--expected-N","--truth","--sequence","--benchmark-manifest")
                @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
            end
            variant=="patch_support" || @test_throws ErrorException GMM._parse_cli(vcat(args,["--selftrain","0"]))
            if variant != "patch_support"
                failed=GMM.Options((f==:learning ? merge(opt.learning,(min_mass=1e20,)) : getfield(opt,f)
                    for f in fieldnames(GMM.Options))...)
                missing=GMM._view_probability(records,features,failed)
                @test all(isnan,missing)
                path=joinpath(dir,variant*"_failed.tsv")
                GMM._write_predictions(path,records,missing,zeros(Int,length(records)))
                _,rows=RU.lobe_table(path)
                @test length(rows)==length(records)
                @test all(r["predicted"]=="?" && r["invalid_reason"]=="no_valid_view" for r in values(rows))
            end
        end
    end
end

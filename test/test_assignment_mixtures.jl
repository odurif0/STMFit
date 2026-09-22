using Test, Statistics, Random, LinearAlgebra, TOML
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
using .GMM.AssignmentMixtures
const RU = GMM.ReconstructedUnitAssignment
const ROOT = dirname(@__DIR__)
config(family) = joinpath(ROOT, "config", "unit_assignment_" * family * ".toml")
settings(family) = RU.load_gmm_learning(RU.load_config(config(family)))

@testset "Explicit isolated learning families, unchanged support policies" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    control = RU.load_config(config("patch_support"))
    for family in ("factor_analyzer", "student_t")
        cfg = RU.load_config(config(family))
        @test settings(family).family == family
        @test settings(family).factor_rank == 4 && settings(family).student_df == 5.
        cfg["model"]["name"] = control["model"]["name"]
        cfg["model"]["gmm_learning_family"] = "gaussian"
        @test cfg == control
        for (section,key,bad) in (("model","gmm_covariance_structure","tied"),
                ("model","gmm_final_covariance","ledoit_wolf"),("model","gmm_final_score","gaussian_density"),
                ("selection","fisher_cv_scheme","scan_hash_twofold"),("selection","gmm_cluster_naming","within_scan_z"),
                ("selection","gmm_training_weighting","equal_scans"),("selection","gmm_seed_aggregation","mean_membership"),
                ("selection","gmm_resampling","whole_scans"),("selection","assignment_training_support","complete_patches"),
                ("preprocessing","gmm_feature_normalization","median_iqr"),("model","gmm_learning_cholesky_guard",1e-7))
            cfg = RU.load_config(config(family)); cfg[section][key] = bad
            @test_throws ArgumentError RU.load_gmm_learning(cfg)
        end
    end
    for file in readdir(joinpath(ROOT,"config");join=true)
        startswith(basename(file),"unit_assignment_") && endswith(file,".toml") || continue
        cfg = TOML.parsefile(file)
        haskey(get(cfg,"model",Dict()),"gmm_final_covariance") || continue
        @test RU.load_gmm_learning(cfg).family == (file == config("factor_analyzer") ? "factor_analyzer" :
            file in (config("student_t"), config("student_density")) ? "student_t" : "gaussian")
    end
    for key in ("gmm_factor_rank", "gmm_learning_maxiter", "gmm_student_df", "gmm_learning_tolerance",
                "gmm_learning_min_mass", "gmm_learning_cholesky_guard"), bad in (nothing, 0, -1, true, Inf, NaN, "5")
        cfg = deepcopy(control)
        bad === nothing ? delete!(cfg["model"],key) : (cfg["model"][key] = bad)
        @test_throws ArgumentError RU.load_gmm_learning(cfg)
    end
    for bad in (nothing,"auto","gaussian_factor_student",true,2)
        cfg = deepcopy(control)
        bad === nothing ? delete!(cfg["model"],"gmm_learning_family") : (cfg["model"]["gmm_learning_family"] = bad)
        @test_throws ArgumentError RU.load_gmm_learning(cfg)
    end
end

# Independent observation-wise equations use full covariance inverse and outer products.
function independent_expectation(X, state, s)
    p,n = size(X); delta=zeros(n,2); kernels=zeros(n,2)
    for c in 1:2
        C = state.covs[c] + s.cholesky_guard * I
        for i in 1:n
            d = X[:,i] - state.means[:,c]
            delta[i,c] = dot(d, inv(C)*d)
            kernels[i,c] = state.weights[c]/sqrt(det(C)) * (s.family == "student_t" ?
                (1 + delta[i,c]/s.student_df)^(-(s.student_df+p)/2) : exp(-delta[i,c]/2))
        end
    end
    return kernels ./ sum(kernels;dims=2), delta
end

function independent_update(X, state, r, delta, s; ridge)
    p,n = size(X); mu=similar(state.means); covs=Matrix{Float64}[]
    for c in 1:2
        mass=sum(r[:,c])
        if s.family == "student_t"
            u = [(s.student_df+p)/(s.student_df+delta[i,c]) for i in 1:n]
            mu[:,c] = sum(r[i,c]*u[i]*X[:,i] for i in 1:n)/sum(r[:,c].*u)
            push!(covs,sum(r[i,c]*u[i]*(X[:,i]-mu[:,c])*(X[:,i]-mu[:,c])' for i in 1:n)/mass+ridge*I)
        else
            L=state.loadings[c]; C=state.covs[c]+s.cholesky_guard*I; q=size(L,2)
            beta=L'*inv(C); V=Matrix{Float64}(I,q,q)-beta*L
            z=[beta*(X[:,i]-state.means[:,c]) for i in 1:n]
            zz=[V+z[i]*z[i]' for i in 1:n]
            T=sum(r[i,c]*[zz[i] z[i]; z[i]' 1.] for i in 1:n)
            B=sum(r[i,c]*X[:,i]*vcat(z[i],1.)' for i in 1:n)
            W=B*inv(T); nextL=W[:,1:q]; mu[:,c]=W[:,end]
            # Full residual expression independently includes both cross terms.
            residual=sum(r[i,c]*(X[:,i]*X[:,i]' - X[:,i]*vcat(z[i],1.)'*W' -
                W*vcat(z[i],1.)*X[:,i]' + W*[zz[i] z[i];z[i]' 1.]*W') for i in 1:n)/mass
            psi=max.(diag(residual).-s.cholesky_guard,ridge)
            push!(covs,nextL*nextL'+Diagonal(psi))
        end
    end
    return mu,covs,vec(sum(r;dims=1))/n
end

@testset "Latent and Student equations, free masses, EM and hard updates" begin
    rng=MersenneTwister(20260922)
    X=hcat(randn(rng,6,55).-1,randn(rng,6,25).+2)
    X[1,:] .+= .4X[2,:]; ridge=1e-6
    for family in ("factor_analyzer","student_t"), seed in 0:2
        s=merge(settings(family),(maxiter=4,tolerance=1e-30))
        initial=GMM._gmm_fit(X,seed;ridge,max_iter=0)
        state=initialize_state(initial...,s;ridge)
        trace=[]
        fit=fit_mixture(X,initial,s;ridge,hard_iterations=2,trace)
        @test length(trace)==6 && !fit.converged && fit.updates==4
        @test isequal(fit,fit_mixture(X,initial,s;ridge,hard_iterations=2))
        for step in trace
            e=expectation(X,state,s); r,d=independent_expectation(X,state,s)
            @test e.responsibilities ≈ r atol=1e-10
            @test e.distances ≈ d atol=1e-10
            if step.stage == "hard"
                a=[argmin(d[i,:]) for i in axes(X,2)]
                r=Float64.([a[i]==c for i in axes(X,2),c in 1:2])
            end
            mu,C,w=independent_update(X,state,r,d,s;ridge)
            @test step.state.means ≈ mu atol=1e-9
            @test all(step.state.covs[c] ≈ C[c] for c in 1:2)
            @test step.state.weights ≈ w atol=1e-10
            @test sum(w) ≈ 1 && w[1] != w[2]
            @test all(C -> isapprox(C,C';atol=1e-12) && isposdef(Symmetric(C)),step.state.covs)
            if family == "factor_analyzer"
                for c in 1:2
                    L,psi=step.state.loadings[c],step.state.noise[c]
                    @test size(L)==(6,4) && minimum(psi)>=ridge
                    @test step.state.covs[c] ≈ L*L'+Diagonal(psi)
                end
            end
            state=step.state
        end
        @test_throws MixtureFitError update_state(X,state,hcat(ones(80),zeros(80)),zeros(80,2),s;ridge)
        @test_throws ArgumentError update_state(X,state,fill(NaN,80,2),zeros(80,2),s;ridge)
        @test_throws ArgumentError fit_mixture(fill(NaN,6,80),initial,s;ridge,hard_iterations=2)
        @test_throws MixtureFitError expectation(X,merge(state,(covs=[-Matrix{Float64}(I,6,6),state.covs[2]],)),s)
    end
    @test_throws ArgumentError initialize_state(zeros(4,2),[Matrix{Float64}(I,4,4) for _ in 1:2],
        [.5,.5],settings("factor_analyzer");ridge)
end

@testset "Student outlier influence, scale denominator and Gaussian limit" begin
    X=[0. 1. -1. 50.; 0. .5 -.5 0.]; p,n=size(X)
    s=settings("student_t"); ridge=1e-6
    state=initialize_state(zeros(2,2),[Matrix{Float64}(I,2,2) for _ in 1:2],[.5,.5],s;ridge)
    e=expectation(X,state,s); next=update_state(X,state,e.responsibilities,e.distances,s;ridge)
    u=(s.student_df+p)./(s.student_df.+e.distances[:,1])
    @test u[end]<u[1]/100
    @test abs(next.means[1,1])<abs(mean(X[1,:]))/50
    c=X.-next.means[:,1]
    @test next.covs[1] ≈ c*Diagonal(u)*c'/n+ridge*I
    @test !isapprox(next.covs[1],c*Diagonal(u)*c'/sum(u)+ridge*I)
    large=merge(s,(student_df=1e12,))
    e=expectation(X,state,large); gaussian=update_state(X,state,e.responsibilities,e.distances,large;ridge)
    @test gaussian.means[:,1] ≈ vec(mean(X;dims=2)) rtol=1e-8
    @test gaussian.covs[1] ≈ cov(X;dims=2,corrected=false)+ridge*I rtol=1e-8
end

@testset "Label-free CLI, final vote reconstruction, unavailable rows and failures" begin
    rng=MersenneTwister(922)
    records=GMM.LobeRecord[]; features=["f$i" for i in 1:5]
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
        for family in ("factor_analyzer","student_t","patch_support")
            args=["--features",file,"--config",config(family),"--view","test="*join(features,','),
                "--seeds","3","--selftrain","2","--interactions"]
            opt=GMM._parse_cli(args); trace=[]
            probs=GMM._view_probability(records,features,opt;diagnostics=trace)
            @test isequal(probs,GMM._view_probability(records,features,opt))
            Z,valid=GMM._standardized_matrix(records,features;interactions=true,
                normalization="mean_sample_std",scale_fallback=1.)
            ids=findall(valid); votes=zeros(length(ids))
            @test length(trace)==3
            for d in trace
                if family == "patch_support"
                    mu=hcat([c.mean for c in d.clusters]...); C=[c.covariance for c in d.clusters]; w=[c.weight for c in d.clusters]
                else
                    state=last(d.clusters).state; mu,C,w=state.means,state.covs,state.weights
                    @test count(t->t.stage=="hard",d.clusters)==2
                end
                scores=[log(w[c])-dot(Z[i,:]-mu[:,c],inv(C[c]+1e-8I)*(Z[i,:]-mu[:,c]))/2 for i in ids,c in 1:2]
                a=[argmax(scores[j,:]) for j in eachindex(ids)]
                high=argmax([mean(records[ids[j]].amplitude for j in eachindex(ids) if a[j]==c) for c in 1:2])
                votes .+= a.==high
            end
            @test probs[ids]==votes/3
            @test count(.!valid)==2 && all(isnan,probs[.!valid])
            renamed=deepcopy(records)
            for r in renamed
                r.file="renamed_"*r.file; r.lobe+=100; r.features["N"]=-999
            end
            @test isequal(probs,GMM._view_probability(renamed,features,opt))
            expected,out=joinpath.(dir,(family*"_expected.tsv",family*"_cli.tsv"))
            GMM._write_predictions(expected,records,probs,Int.(valid))
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args --out $out`)
            @test read(expected)==read(out)
            @test length(last(RU.lobe_table(out)))==length(records)
            for flag in ("--expected-N","--truth","--sequence","--benchmark-manifest")
                @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
            end
            if family != "patch_support"
                @test_throws ErrorException GMM._parse_cli(vcat(args,["--selftrain","0"]))
                failed=GMM.Options((f==:learning ? merge(opt.learning,(min_mass=1e20,)) : getfield(opt,f) for f in fieldnames(GMM.Options))...)
                missing=GMM._view_probability(records,features,failed)
                @test all(isnan,missing)
                path=joinpath(dir,family*"_failed.tsv")
                GMM._write_predictions(path,records,missing,zeros(Int,length(records)))
                _,rows=RU.lobe_table(path)
                @test length(rows)==length(records)
                @test all(r["predicted"]=="?" && r["invalid_reason"]=="no_valid_view" for r in values(rows))
            end
        end
    end
end

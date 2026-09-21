using Test, Statistics, Random, LinearAlgebra, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_tied_covariance.toml")

@testset "Explicit covariance sharing only, with no combined experiment" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_gmm_covariance_structure(a) == "full"
    @test load_gmm_covariance_structure(b) == "tied"
    b["model"]["gmm_covariance_structure"] = a["model"]["gmm_covariance_structure"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    @test_throws ErrorException ReconstructedRepresentationDiagnostics.load_inputs(Dict(), CANDIDATE)
    for name in readdir(joinpath(ROOT, "config"))
        startswith(name, "unit_assignment_") && endswith(name, ".toml") || continue
        cfg = TOML.parsefile(joinpath(ROOT, "config", name))
        haskey(get(cfg, "model", Dict()), "gmm_final_covariance") || continue
        @test load_gmm_covariance_structure(cfg) == (name == basename(CANDIDATE) ? "tied" : "full")
    end
    mktempdir() do dir
        for value in (nothing, "auto", "diag", true, 1, 0.5, NaN)
            cfg = deepcopy(a)
            value === nothing ? delete!(cfg["model"], "gmm_covariance_structure") :
                (cfg["model"]["gmm_covariance_structure"] = value)
            @test_throws ArgumentError load_gmm_covariance_structure(cfg)
            path = joinpath(dir, "bad.toml")
            open(io -> TOML.print(io,cfg), path, "w")
            @test_throws ArgumentError load_config(path)
        end
    end
    for (section,key,value) in (("model","gmm_final_covariance","ledoit_wolf"),
            ("model","gmm_final_score","gaussian_density"),
            ("selection","gmm_training_weighting","equal_scans"),
            ("selection","gmm_seed_aggregation","mean_membership"),
            ("selection","gmm_resampling","whole_scans"),
            ("selection","assignment_training_support","complete_patches"),
            ("preprocessing","gmm_feature_normalization","median_iqr"))
        cfg = load_config(CANDIDATE); cfg[section][key] = value
        @test_throws ArgumentError load_gmm_covariance_structure(cfg)
    end
end

# Independent outer-product arithmetic, not the production covariance helper.
function pooled(X, means, R)
    sum(R[i,c] * ((X[:,i]-means[:,c])*(X[:,i]-means[:,c])')
        for i in axes(X,2), c in axes(means,2)) / size(X,2)
end
function independent_em(X, means, covariance, weights; ridge)
    n, k = size(X,2), length(weights)
    L = cholesky(Symmetric(covariance) + 1e-8I).L
    R = zeros(n,k)
    for i in 1:n
        logs = [log(weights[c]) - sum(abs2,L\(X[:,i]-means[:,c]))/2 for c in 1:k]
        R[i,:] = exp.(logs .- maximum(logs)); R[i,:] ./= sum(R[i,:])
    end
    masses = vec(sum(R;dims=1))
    mu = hcat([sum(R[i,c]*X[:,i] for i in 1:n)/masses[c] for c in 1:k]...)
    mu, pooled(X,mu,R) + ridge*I, masses/n
end

@testset "Pooled within-component scatter, not between-group or equal-group variance" begin
    X = [0. 1. 2. 10. 14.; 0. 2. 1. 12. 10.]
    R = [1. 0.; 1. 0.; 1. 0.; 0. 1.; 0. 1.]
    mu = hcat(vec(mean(X[:,1:3];dims=2)),vec(mean(X[:,4:5];dims=2)))
    S = pooled(X,mu,R); result = GMM._tied_covariance(X,mu,R;ridge=1e-6)
    @test result.sample ≈ S
    @test result.covariance ≈ S + 1e-6I
    a,b = X[:,1:3].-mu[:,1], X[:,4:5].-mu[:,2]
    @test S ≈ (a*a' + b*b')/5
    @test !isapprox(S,(a*a'/3 + b*b'/2)/2)
    @test !isapprox(S,cov(X;dims=2,corrected=false))
    shift = [100.,-23.]
    @test GMM._tied_covariance(X.+shift,mu.+shift,R;ridge=1e-6).sample ≈ S
    @test GMM._tied_covariance(zeros(2,5),zeros(2,2),R;ridge=1e-6).covariance == 1e-6Matrix(I,2,2)
    for bad in (R.*2, fill(NaN,5,2), -R, zeros(4,2))
        @test_throws ArgumentError GMM._tied_covariance(X,mu,bad;ridge=1e-6)
    end
    for bad in (0.,-1.,Inf,NaN,true)
        @test_throws ArgumentError GMM._tied_covariance(X,mu,R;ridge=bad)
    end
end

@testset "Shared initialization and every EM/hard update match independent equations" begin
    rng = MersenneTwister(922)
    X = hcat(randn(rng,3,25).-1.5, randn(rng,3,11).+1.5)
    ridge = 1e-6
    for seed in 0:2
        km = GMM.kmeans(X,2;maxiter=100,rng=MersenneTwister(seed),display=:none)
        R = Float64.([a==c for a in km.assignments,c in 1:2])
        mu = hcat([vec(mean(X[:,findall(==(c),km.assignments)];dims=2)) for c in 1:2]...)
        S = pooled(X,mu,R)+ridge*I
        fit = GMM._gmm_fit(X,seed;ridge,max_iter=0,covariance_structure="tied")
        @test fit[1] ≈ mu
        @test fit[2][1] == fit[2][2]
        @test fit[2][1] !== fit[2][2]
        @test fit[2][1] ≈ S
        @test fit[3] == [.5,.5] # Initialization only, not a composition constraint.
        weights = copy(fit[3])
        for step in 1:3
            mu,S,weights = independent_em(X,mu,S,weights;ridge)
            fit = GMM._gmm_fit(X,seed;ridge,max_iter=step,tol=0.,covariance_structure="tied")
            @test fit[1] ≈ mu atol=1e-12
            @test fit[2][1] == fit[2][2]
            @test fit[2][1] ≈ S atol=1e-12
            @test fit[3] ≈ weights atol=1e-12
            @test sum(fit[3]) ≈ 1
            @test fit[3][1] != fit[3][2]
        end
        start = deepcopy(fit)
        for step in 1:2
            L = cholesky(Symmetric(S)+1e-8I).L
            assignments = [argmin([sum(abs2,L\(X[:,i]-mu[:,c])) for c in 1:2]) for i in axes(X,2)]
            R = Float64.([a==c for a in assignments,c in 1:2])
            mu = hcat([vec(mean(X[:,findall(==(c),assignments)];dims=2)) for c in 1:2]...)
            weights = vec(sum(R;dims=1))/size(X,2)
            S = pooled(X,mu,R)+ridge*I
            diagnostic = []
            result = GMM._mahalanobis_self_train(X,deepcopy(start)...;iters=step,
                covariance_mode="ridge",ridge,covariance_structure="tied",diagnostics=diagnostic)
            @test result[1] ≈ mu atol=1e-12
            @test result[2][1] == result[2][2]
            @test result[2][1] ≈ S atol=1e-12
            @test result[3] == weights
            @test length(diagnostic) == 2
            @test diagnostic[1].covariance == diagnostic[2].covariance
            @test all(d.sample ≈ S-ridge*I for d in diagnostic)
        end
        legacy = GMM._gmm_fit(X,seed;ridge)
        @test isequal(legacy,GMM._gmm_fit(X,seed;ridge,covariance_structure="full"))
        @test isequal(GMM._mahalanobis_self_train(X,deepcopy(legacy)...;iters=2,covariance_mode="ridge",ridge),
            GMM._mahalanobis_self_train(X,deepcopy(legacy)...;iters=2,covariance_mode="ridge",ridge,covariance_structure="full"))
        @test_throws ArgumentError GMM._gmm_fit(X,seed;ridge,covariance_structure="bad")
        @test_throws ArgumentError GMM._gmm_fit(X,seed;ridge,covariance_structure="tied",observation_weights=ones(size(X,2)))
        @test_throws ArgumentError GMM._mahalanobis_self_train(X,deepcopy(start)...;iters=2,
            covariance_mode="ledoit_wolf",ridge,covariance_structure="tied")
    end
end

@testset "Label-free naming, unavailable rows, deterministic vote and CLI replay" begin
    rng = MersenneTwister(913)
    records = GMM.LobeRecord[]
    for file in 1:6, lobe in 1:(7+file)
        high = lobe > .7(7+file)
        push!(records,GMM.LobeRecord("synthetic_$file.sxm",lobe,1+2high+.1randn(rng),
            Dict("f1"=>(high ? 2. : -1.)+.5randn(rng)+file,
                 "f2"=>(high ? 1. : -.5)+.5randn(rng)-file)))
    end
    records[2].features["f2"] = NaN
    push!(records,GMM.LobeRecord("unavailable.sxm",1,1.,Dict("f1"=>NaN,"f2"=>NaN)))
    mktempdir() do dir
        features = joinpath(dir,"features.tsv")
        write_table(features,["file","lobe","amplitude","f1","f2"],
            [Dict("file"=>r.file,"lobe"=>r.lobe,"amplitude"=>r.amplitude,
                  "f1"=>r.features["f1"],"f2"=>r.features["f2"]) for r in records])
        args = ["--features",features,"--config",CANDIDATE,"--view","test=f1,f2",
            "--seeds","3","--selftrain","2","--interactions"]
        opt = GMM._parse_cli(args)
        @test opt.covariance_structure == "tied"
        @test opt.seed_aggregation == "hard_vote" && opt.resampling == "none"
        diagnostic = []; probs = GMM._view_probability(records,["f1","f2"],opt;diagnostics=diagnostic)
        @test isequal(probs,GMM._view_probability(records,["f1","f2"],opt))
        Z,valid = GMM._standardized_matrix(records,["f1","f2"];interactions=true,
            normalization="mean_sample_std",scale_fallback=1.)
        ids = findall(valid); votes = zeros(length(ids))
        @test length(diagnostic)==3
        for d in diagnostic
            @test d.indices == ids && d.observation_weights === nothing
            @test length(d.clusters) == 2
            @test d.clusters[1].covariance == d.clusters[2].covariance
            @test sum(c.weight for c in d.clusters) ≈ 1
            scores = zeros(length(ids),2)
            for (j,i) in enumerate(ids), c in d.clusters
                L = cholesky(Symmetric(c.covariance)+1e-8I).L
                scores[j,c.component] = log(c.weight)-sum(abs2,L\(Z[i,:]-c.mean))/2
            end
            assignment = [argmax(scores[j,:]) for j in eachindex(ids)]
            means = [mean(records[ids[j]].amplitude for j in eachindex(ids) if assignment[j]==c) for c in 1:2]
            votes .+= assignment .== argmax(means)
        end
        @test probs[ids] == votes/3
        @test count(.!valid)==2 && all(isnan,probs[.!valid])
        @test all(p->isfinite(p)&&0<=p<=1,probs[ids])
        renamed = deepcopy(records)
        for r in renamed
            r.file="renamed_"*r.file; r.lobe+=100; r.features["N"]=-999
        end
        @test isequal(probs,GMM._view_probability(renamed,["f1","f2"],opt))
        expected,out = joinpath.(dir,("expected.tsv","cli.tsv"))
        GMM._write_predictions(expected,records,probs,Int.(valid))
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args --out $out`)
        @test read(expected)==read(out)
        @test length(last(lobe_table(out)))==length(records)
        for flag in ("--expected-N","--truth","--sequence","--benchmark-manifest")
            @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
        end
    end
end

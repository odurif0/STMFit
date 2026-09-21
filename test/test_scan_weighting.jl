using Test, Statistics, Random, LinearAlgebra, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_scan_weighting.toml")

@testset "Only the explicit observation-weighting policy changes" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_gmm_weighting(a) == "equal_lobes"
    @test load_gmm_weighting(b) == "equal_scans"
    b["selection"]["gmm_training_weighting"] = a["selection"]["gmm_training_weighting"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    @test_throws ErrorException ReconstructedRepresentationDiagnostics.load_inputs(Dict(), CANDIDATE)
    for value in (nothing, "guess", true, 1, 0.5, NaN)
        cfg = deepcopy(a)
        value === nothing ? delete!(cfg["selection"], "gmm_training_weighting") :
            (cfg["selection"]["gmm_training_weighting"] = value)
        @test_throws ArgumentError load_gmm_weighting(cfg)
    end
    cfg = load_config(CANDIDATE)
    cfg["model"]["gmm_final_covariance"] = "ledoit_wolf"
    @test_throws ArgumentError load_gmm_weighting(cfg)
end

@testset "Equal scan totals use usable rows, never N, labels or ordinal position" begin
    records = [GMM.LobeRecord(f, i, 1., Dict("N" => 999.))
        for (f, n) in (("a.sxm", 4), ("b.sxm", 2), ("unavailable.sxm", 1)) for i in 1:n]
    idxs = [1, 2, 3, 5, 6] # One inadmissible row and an entirely unavailable scan.
    @test GMM._observation_weights(records, idxs, "equal_lobes") === nothing
    w = GMM._observation_weights(records, idxs, "equal_scans")
    @test w ≈ [5/6, 5/6, 5/6, 5/4, 5/4]
    @test sum(w) ≈ length(idxs)
    @test sum(w[1:3]) ≈ sum(w[4:5])
    @test isempty(GMM._observation_weights(records, Int[], "equal_scans"))
    @test_throws ArgumentError GMM._observation_weights(records, idxs, "unknown")
    renamed = deepcopy(records)
    for r in renamed
        r.file = "renamed_" * r.file; r.lobe += 100; r.features["N"] = -1
    end
    @test GMM._observation_weights(renamed, idxs, "equal_scans") == w
    perm = [3, 1, 5, 2, 4]
    @test GMM._observation_weights(records, idxs[perm], "equal_scans") == w[perm]
    @test GMM._observation_weights(records, [1, 2, 5, 6], "equal_scans") == ones(4)
    @test GMM._observation_weights(records, [5], "equal_scans") == [1.]
    for bad in ([1., 2.], [1., 0., 1.], [1., -1., 1.], [1., NaN, 1.], [1., Inf, 1.], [true, true, true])
        @test_throws ArgumentError GMM._check_observation_weights(bad, 3)
    end
end

@testset "Weighted initialization and one EM update match explicit arithmetic" begin
    rng = MersenneTwister(410)
    X = hcat(randn(rng, 2, 17) .- 2, randn(rng, 2, 9) .+ 3)
    w = vcat(fill(0.5, 17), fill(2., 9)); w .*= length(w)/sum(w)
    ridge = 1e-6
    for seed in 0:2
        srng = MersenneTwister(seed)
        seeds = GMM._weighted_seeds(X, w, srng)
        @test length(unique(seeds)) == 2
        @test seeds == GMM._weighted_seeds(X, w, MersenneTwister(seed))
        km = GMM.kmeans(X, 2; weights=w, init=seeds, maxiter=100, rng=srng, display=:none)
        means, covs, mixing = GMM._gmm_fit(X, seed; ridge, max_iter=0, observation_weights=w)
        @test mixing == [0.5, 0.5] # Starting values only; never a constraint.
        for c in 1:2
            ids = findall(==(c), km.assignments)
            mass = sum(w[ids])
            mu = sum(w[i] .* X[:,i] for i in ids) / mass
            covariance = sum(w[i] .* ((X[:,i]-mu)*(X[:,i]-mu)') for i in ids) / mass
            @test means[:,c] ≈ mu
            @test covs[c] ≈ covariance + ridge*I
        end
        resp = zeros(size(X,2), 2)
        for i in axes(X,2)
            logs = [log(mixing[c]) + GMM._gmm_log_density(X[:,i],means[:,c],covs[c]) for c in 1:2]
            resp[i,:] = exp.(logs .- maximum(logs)); resp[i,:] ./= sum(resp[i,:])
        end
        after = GMM._gmm_fit(X, seed; ridge, max_iter=1, observation_weights=w)
        for c in 1:2
            wr = w .* resp[:,c]; mass = sum(wr)
            mu = sum(wr[i] .* X[:,i] for i in axes(X,2)) / mass
            covariance = sum(wr[i] .* ((X[:,i]-mu)*(X[:,i]-mu)') for i in axes(X,2)) / mass
            @test after[1][:,c] ≈ mu
            @test after[2][c] ≈ covariance + ridge*I
            @test after[3][c] ≈ mass/sum(w)
        end
        @test sum(after[3]) ≈ 1
        @test after[3][1] != after[3][2]
        diag = []
        fit = GMM._mahalanobis_self_train(X, deepcopy(after)...; iters=2,
            covariance_mode="ridge", ridge, observation_weights=w, diagnostics=diag)
        @test length(diag) == 2
        for c in diag
            ids = c.members; mass = sum(w[ids])
            mu = sum(w[i] .* X[:,i] for i in ids) / mass
            covariance = sum(w[i] .* ((X[:,i]-mu)*(X[:,i]-mu)') for i in ids) / mass
            @test c.mean ≈ mu
            @test c.sample ≈ covariance
            @test c.covariance ≈ covariance + ridge*I
            @test c.weight ≈ mass/sum(w)
            @test isposdef(Symmetric(c.covariance))
        end
        @test sum(fit[3]) ≈ 1
        @test_throws ErrorException GMM._mahalanobis_self_train(X, deepcopy(after)...;
            iters=2, covariance_mode="ledoit_wolf", ridge, observation_weights=w)
        # Explicit absence of weights is exactly the legacy route.
        @test isequal(GMM._gmm_fit(X, seed; ridge), GMM._gmm_fit(X, seed; ridge, observation_weights=nothing))
    end
    @test length(unique(GMM._weighted_seeds(zeros(2,4), ones(4), MersenneTwister(1)))) == 2
end

@testset "Physical naming is weighted; prediction-only amplitudes cannot rename groups" begin
    records = [GMM.LobeRecord(i <= 4 ? "a.sxm" : "b.sxm", i, a, Dict{String,Float64}())
        for (i,a) in enumerate([9.,9.,9.,6.,0.,6.,1e12])]
    assignments = [1,1,1,2,1,2,1]; eligible = [trues(6);false]
    w = [GMM._observation_weights(records, 1:6, "equal_scans");0.]
    unweighted = GMM._cluster_amplitude_means(records, 1:7, assignments, eligible, nothing)
    weighted = GMM._cluster_amplitude_means(records, 1:7, assignments, eligible, w)
    @test unweighted[1] == 6.75 && unweighted[2] == 6.
    @test weighted[1] ≈ 5.4
    @test weighted[2] == 6.
    records[7].amplitude = -1e12
    @test GMM._cluster_amplitude_means(records, 1:7, assignments, eligible, w) == weighted
end

@testset "CLI/API deterministic identity, per-view eligibility and unchanged scaling" begin
    rng = MersenneTwister(840)
    records = GMM.LobeRecord[]
    for file in 1:6, lobe in 1:(6+file)
        high = lobe > 0.7(6+file)
        push!(records, GMM.LobeRecord("synthetic_$file.sxm", lobe, high ? 4. : 1.,
            Dict("f1" => (high ? 4. : -2.) + .3randn(rng) + file,
                 "f2" => (high ? 3. : -1.) + .4randn(rng) - file)))
    end
    records[2].features["f2"] = NaN
    push!(records, GMM.LobeRecord("unavailable.sxm",1,1.,Dict("f1"=>NaN,"f2"=>NaN)))
    mktempdir() do dir
        path = joinpath(dir,"features.tsv")
        write_table(path,["file","lobe","amplitude","f1","f2"],
            [Dict("file"=>r.file,"lobe"=>r.lobe,"amplitude"=>r.amplitude,
                  "f1"=>r.features["f1"],"f2"=>r.features["f2"]) for r in records])
        args = ["--features",path,"--config",CANDIDATE,"--view","test=f1,f2","--seeds","2","--selftrain","2","--interactions"]
        opt = GMM._parse_cli(args)
        @test opt.training_weighting == "equal_scans"
        @test opt.normalization == "mean_sample_std" && opt.scale_fallback == 1.
        diag = []; p = GMM._view_probability(records,["f1","f2"],opt;diagnostics=diag)
        @test isequal(p,GMM._view_probability(records,["f1","f2"],opt))
        @test length(p) == length(records) && isnan(p[2]) && isnan(p[end])
        @test all(isfinite,p[setdiff(1:length(p),[2,length(p)])])
        for d in diag
            @test !(2 in d.indices) && !(length(records) in d.indices)
            @test d.observation_weights == GMM._observation_weights(records,d.indices,"equal_scans")
            totals = [sum(d.observation_weights[j] for (j,i) in enumerate(d.indices) if records[i].file==f)
                for f in unique(records[i].file for i in d.indices)]
            @test maximum(totals) ≈ minimum(totals)
            @test sum(c.weight for c in d.clusters) ≈ 1
        end
        expected, out = joinpath.(dir,("expected.tsv","cli.tsv"))
        GMM._write_predictions(expected,records,p,Int.(isfinite.(p)))
        @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args --out $out`)
        @test read(expected) == read(out)
        for flag in ("--expected-N","--truth","--sequence")
            @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
        end
    end
end

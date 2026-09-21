using Test, Statistics, Random, LinearAlgebra, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
include(joinpath(@__DIR__, "lib", "reconstructed_representation_diagnostics.jl"))
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_continuous_vote.toml")

@testset "One explicit aggregation change, no training or threshold change" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_gmm_seed_aggregation(a) == "hard_vote"
    @test load_gmm_seed_aggregation(b) == "mean_membership"
    b["selection"]["gmm_seed_aggregation"] = a["selection"]["gmm_seed_aggregation"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    @test_throws ErrorException ReconstructedRepresentationDiagnostics.load_inputs(Dict(), CANDIDATE)
    for name in readdir(joinpath(ROOT, "config"))
        startswith(name, "unit_assignment_") && endswith(name, ".toml") || continue
        cfg = TOML.parsefile(joinpath(ROOT, "config", name))
        haskey(get(cfg, "selection", Dict()), "gmm_training_weighting") || continue
        @test load_gmm_seed_aggregation(cfg) ==
            (name == "unit_assignment_continuous_vote.toml" ? "mean_membership" : "hard_vote")
    end
    mktempdir() do dir
        for value in (nothing, "guess", "posterior", true, 1, 0.5, NaN)
            cfg = deepcopy(a)
            value === nothing ? delete!(cfg["selection"], "gmm_seed_aggregation") :
                (cfg["selection"]["gmm_seed_aggregation"] = value)
            @test_throws ArgumentError load_gmm_seed_aggregation(cfg)
            path = joinpath(dir, "bad.toml")
            open(io -> TOML.print(io, cfg), path, "w")
            @test_throws ArgumentError load_config(path)
        end
    end
end

@testset "Continuous contribution, legacy ties and no confidence recalibration" begin
    for resp in ([0.1,0.9], [0.51,0.49], [0.5,0.5], [0.,1.], [1.,0.]), high in 1:2
        @test GMM._seed_contribution(resp,high,"hard_vote") == (argmax(resp)==high ? 1. : 0.)
        @test GMM._seed_contribution(resp,high,"mean_membership") == resp[high]
    end
    scores = [[0.51,0.49], [0.51,0.49], [0.01,0.99]]
    @test mean(GMM._seed_contribution(r,1,"hard_vote") for r in scores) == 2/3
    @test mean(GMM._seed_contribution(r,1,"mean_membership") for r in scores) ≈ 1.03/3
    @test_throws ArgumentError GMM._seed_contribution([0.1,0.9],1,"unknown")
end

function fixture_records()
    rng = MersenneTwister(921)
    records = GMM.LobeRecord[]
    for file in 1:6, lobe in 1:(7+file)
        high = lobe > 0.65(7+file)
        push!(records, GMM.LobeRecord("synthetic_$file.sxm", lobe,
            1 + 2high + .1randn(rng),
            Dict("f1" => (high ? 1. : -.5) + randn(rng) + file,
                 "f2" => (high ? .4 : -.2) + randn(rng) - file)))
    end
    records[2].features["f2"] = NaN
    push!(records,GMM.LobeRecord("unavailable.sxm",1,1.,Dict("f1"=>NaN,"f2"=>NaN)))
    records
end

@testset "Identical learned parameters, independent aggregation and label-free boundaries" begin
    records = fixture_records()
    mktempdir() do dir
        features = joinpath(dir, "features.tsv")
        write_table(features,["file","lobe","amplitude","f1","f2"],
            [Dict("file"=>r.file,"lobe"=>r.lobe,"amplitude"=>r.amplitude,
                  "f1"=>r.features["f1"],"f2"=>r.features["f2"]) for r in records])
        common = ["--features",features,"--view","test=f1,f2","--seeds","3","--selftrain","2","--interactions"]
        options = [GMM._parse_cli(vcat(common,["--config",cfg])) for cfg in (CONTROL,CANDIDATE)]
        @test options[1].seed_aggregation == "hard_vote"
        @test options[2].seed_aggregation == "mean_membership"
        diagnostics = [[],[]]
        probs = [GMM._view_probability(records,["f1","f2"],options[j];diagnostics=diagnostics[j]) for j in 1:2]
        @test isequal(diagnostics[1],diagnostics[2]) # Every seed: same members, weights, means and covariances.
        @test length(diagnostics[1]) == 3
        Z, valid = GMM._standardized_matrix(records,["f1","f2"];interactions=true,
            normalization="mean_sample_std",scale_fallback=1.)
        ids = findall(valid)
        sums = zeros(length(ids),2)
        for d in diagnostics[1]
            @test d.indices == ids && d.observation_weights === nothing
            @test length(d.clusters) == 2
            @test sum(c.weight for c in d.clusters) ≈ 1
            @test all(c.weight == length(c.members)/length(ids) for c in d.clusters)
            # Reconstruct the unchanged Mahalanobis scores from fitted moments,
            # not via the production score/contribution helpers.
            scores = zeros(length(ids),2)
            for (j,i) in enumerate(ids), c in d.clusters
                z = cholesky(Symmetric(c.covariance) + 1e-8I).L \ (Z[i,:] - c.mean)
                scores[j,c.component] = log(c.weight) - sum(abs2,z)/2
            end
            memberships = exp.(scores .- maximum(scores;dims=2))
            memberships ./= sum(memberships;dims=2)
            assignment = [argmax(scores[j,:]) for j in axes(scores,1)]
            amplitude = [mean(records[ids[j]].amplitude for j in eachindex(ids) if assignment[j]==c) for c in 1:2]
            high = argmax(amplitude)
            sums[:,1] .+= assignment .== high
            sums[:,2] .+= memberships[:,high]
        end
        @test any(d.clusters[1].weight != d.clusters[2].weight for d in diagnostics[1])
        @test probs[1][ids] == sums[:,1]/3
        @test probs[2][ids] ≈ sums[:,2]/3 atol=2e-15 rtol=2e-15
        @test any(p -> abs(3p-round(3p)) > 1e-6, probs[2][ids]) # Genuinely continuous on this fixture.
        @test probs[1][ids] != probs[2][ids]
        for (j,opt) in enumerate(options)
            @test length(probs[j]) == length(records)
            @test all(isnan,probs[j][.!valid]) && count(.!valid)==2
            @test all(p -> isfinite(p) && 0<=p<=1,probs[j][valid])
            @test isequal(probs[j],GMM._view_probability(records,["f1","f2"],opt))
            renamed = deepcopy(records)
            for r in renamed
                r.file="renamed_"*r.file; r.lobe+=100; r.features["N"]=-999
            end
            @test isequal(probs[j],GMM._view_probability(renamed,["f1","f2"],opt))
            expected, output = joinpath.(dir,("expected_$j.tsv","cli_$j.tsv"))
            GMM._write_predictions(expected,records,probs[j],Int.(valid))
            args=vcat(common,["--config",j==1 ? CONTROL : CANDIDATE,"--out",output])
            @test success(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"build_labelfree_gmm_predictions.jl")) $args`)
            @test read(expected)==read(output)
            @test length(last(lobe_table(output)))==length(records)
            for flag in ("--expected-N","--truth","--sequence","--benchmark-manifest")
                @test_throws ErrorException GMM._parse_cli(vcat(args,[flag,"forbidden"]))
            end
        end
    end
end

@testset "Final vote threshold and abstention are unchanged" begin
    mktempdir() do dir
        features, km, gmm, out = joinpath.(dir,("features.tsv","km.tsv","gmm.tsv","out.tsv"))
        rows = [Dict("file"=>"toy.sxm","lobe"=>i) for i in 1:4]
        write_table(features,["file","lobe"],rows)
        header=["file","lobe","predicted","probability_1"]
        write_table(km,header,[merge(r,Dict("predicted"=>"0","probability_1"=>0.)) for r in rows])
        write_table(gmm,header,[merge(r,Dict("predicted"=>(p===nothing ? "?" : "1"),
            "probability_1"=>(p===nothing ? "NA" : p))) for (r,p) in zip(rows,(1.,0.999999,0.6,nothing))])
        write_soft_vote(features,km,gmm,out,CANDIDATE)
        _,table=lobe_table(out)
        @test [table[("toy.sxm",i)]["predicted"] for i in 1:4]==["1","0","0","?"]
        @test table[("toy.sxm",1)]["confidence"]=="0.00000000"
        @test table[("toy.sxm",4)]["invalid_reason"]=="unavailable_gmm"
        rounded=joinpath(dir,"rounded.tsv")
        records=[GMM.LobeRecord("toy.sxm",i,1.,Dict{String,Float64}()) for i in 1:4]
        GMM._write_predictions(rounded,records,[.499999999,.500000001,.5,NaN],[1,1,1,0])
        _,saved=lobe_table(rounded)
        @test [saved[("toy.sxm",i)]["predicted"] for i in 1:4]==["0","1","1","?"]
        @test all(saved[("toy.sxm",i)]["probability_1"]=="0.50000000" for i in 1:3)
    end
end

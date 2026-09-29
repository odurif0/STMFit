# Molecule-balanced cohort learning: weights from label-free repeated-scan tracks.
#   julia --project=. test/test_molecule_weighting.jl
using Test, TOML
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
module GMM
include(joinpath(@__DIR__, "build_labelfree_gmm_predictions.jl"))
end
module KM
include(joinpath(@__DIR__, "build_labelfree_unit_predictions.jl"))
end
const ROOT = dirname(@__DIR__)
const CONTROL = joinpath(ROOT, "config", "unit_assignment_patch_support.toml")
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_molecule_balanced.toml")

@testset "only the declared training weighting differs from the control" begin
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_gmm_weighting(b) == "equal_molecules"
    @test b["selection"]["kmeans_training_weighting"] == "equal_molecules"
    delete!(b["selection"], "kmeans_training_weighting")
    b["selection"]["gmm_training_weighting"] = a["selection"]["gmm_training_weighting"]
    b["model"]["name"] = a["model"]["name"]
    @test a == b
    cfg = load_config(CANDIDATE); cfg["model"]["gmm_final_covariance"] = "ledoit_wolf"
    @test_throws ArgumentError load_gmm_weighting(cfg)
end

@testset "each molecule track gets equal total mass" begin
    mktempdir() do dir
        groups = joinpath(dir, "groups.tsv")
        write(groups, "file\tgroup\na.sxm\ttrack1\nb.sxm\ttrack1\nc.sxm\ttrack2\n")
        g = read_training_groups(groups)
        @test g == Dict("a.sxm" => "track1", "b.sxm" => "track1", "c.sxm" => "track2")
        recs = [GMM.LobeRecord(f, i, 1.0, Dict{String,Float64}()) for (f, n) in (("a.sxm", 6), ("b.sxm", 6), ("c.sxm", 6)) for i in 1:n]
        idxs = collect(eachindex(recs))
        w = GMM._observation_weights(recs, idxs, "equal_molecules", g)
        @test sum(w[1:12]) ≈ sum(w[13:18]) && sum(w) ≈ length(idxs)
        @test GMM._observation_weights(recs, idxs, "equal_lobes", g) === nothing
        @test_throws ArgumentError GMM._observation_weights(recs, idxs, "equal_molecules", Dict("a.sxm" => "t"))
        krecs = [KM.LobeRecord(r.file, r.lobe, r.amplitude, r.features) for r in recs]
        kw = KM.group_weights(krecs, idxs, groups)
        @test kw ≈ w
        @test KM.group_weights(krecs, idxs, "") === nothing
        write(groups, "file\tgroup\na.sxm\ttrack1\na.sxm\ttrack2\n")
        @test_throws ErrorException read_training_groups(groups)
    end
end

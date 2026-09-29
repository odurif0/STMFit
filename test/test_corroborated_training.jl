# Cohort learning restricted to scans whose count is corroborated by repeated scans.
#   julia --project=. test/test_corroborated_training.jl
using Test, TOML, Random
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
const CANDIDATE = joinpath(ROOT, "config", "unit_assignment_corroborated_training.toml")

function consensus_summary(dir, rules)
    path = joinpath(dir, "consensus_summary_final.tsv")
    open(path, "w") do io
        println(io, "filepath\tN_selected\tcount_rule")
        for (f, r) in rules; println(io, "$f\t6\t$r"); end
    end
    return path
end

@testset "only the declared training-scan policy differs from the control" begin
    a, b = load_config(CONTROL), load_config(CANDIDATE)
    @test load_training_scans(a) == "all" && load_training_scans(b) == "corroborated_counts"
    delete!(b["selection"], "assignment_training_scans"); b["model"]["name"] = a["model"]["name"]
    @test a == b
    bad = load_config(CONTROL); bad["selection"]["assignment_training_scans"] = "labels"
    @test_throws ArgumentError load_training_scans(bad)
end

@testset "corroborated scans and CLI pairing" begin
    mktempdir() do dir
        s = consensus_summary(dir, ["a.sxm" => "agrees", "b.sxm" => "consensus_applied", "c.sxm" => "track_too_short"])
        @test read_corroborated_scans(s) == Set(["a.sxm"])
        @test KM.corroborated_scans(s) == Set(["a.sxm"])
        feats = joinpath(dir, "f.tsv"); write(feats, "file\tlobe\tamplitude\tx\na.sxm\t1\t1.0\t0.1\n")
        base = ["--features", feats, "--out", joinpath(dir, "o.tsv"), "--view", "v=x"]
        @test GMM._parse_cli(vcat(base, ["--config", CANDIDATE, "--training-scans", s])).training_scans == Set(["a.sxm"])
        @test_throws ErrorException GMM._parse_cli(vcat(base, ["--config", CANDIDATE]))
        @test_throws ErrorException GMM._parse_cli(vcat(base, ["--config", CONTROL, "--training-scans", s]))
        @test GMM._parse_cli(vcat(base, ["--config", CONTROL])).training_scans === nothing
    end
end

@testset "k-means learns on corroborated scans and still assigns every row" begin
    mktempdir() do dir
        rng = MersenneTwister(3)
        recs = KM.LobeRecord[]
        for s in 1:20, l in 1:6   # corroborated: two well-separated types per chain
            hi = l in (2, 5)
            push!(recs, KM.LobeRecord("s$s.sxm", l, hi ? 2.0 : 1.0, Dict("x" => (hi ? 3.0 : 0.0) + 0.1randn(rng))))
        end
        for l in 1:6              # uncorroborated scan with an extreme artifact
            push!(recs, KM.LobeRecord("bad.sxm", l, 1.0, Dict("x" => 50.0 + l)))
        end
        rules = vcat(["s$s.sxm" => "agrees" for s in 1:20], ["bad.sxm" => "consensus_fit_failed"])
        s = consensus_summary(dir, rules)
        X = reshape([r.features["x"] for r in recs], :, 1)
        opt = KM.Options("", "", "", "", ["v" => ["x"]], 0, 5, false, s)
        p = KM._view_probability_trained(recs, X, collect(eachindex(recs)), ["x"], opt)
        @test all(isfinite, p)
        truth = [r.lobe in (2, 5) && r.file != "bad.sxm" for r in recs]
        @test all(p[i] == (truth[i] ? 1.0 : 0.0) for i in eachindex(recs) if recs[i].file != "bad.sxm")
        @test all(p[i] == 1.0 for i in eachindex(recs) if recs[i].file == "bad.sxm")  # assigned, not learned from
    end
end

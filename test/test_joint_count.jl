# Joint count of repeated scans from per-scan score curves (campaign J, stage 1).
#   julia --project=. test/test_joint_count.jl
using Test
using STMFit.JointCount

@testset "score tables and effective curves" begin
    mktempdir() do tmp
        d = mkpath(joinpath(tmp, "chunk1", "scanA"))
        head = "N\tGCV\tdelta_raw\tdelta_valid\tBIC\tGCV\tAICc\tvalid\traw_best\tselected\treason\tchi2\tcv_nll_mean\tcv_nll_std"
        row(n, s, v) = join([n, s, 0, 0, 1, s, 1, v, false, false, "ok", 0.1, "NaN", "NaN"], '\t')
        write(joinpath(d, "ell_scores.tsv"), join([head, row(4, 3.0, true), row(5, 1.0, true), row(6, 2.0, true), row(7, 0.5, false)], "\n") * "\n")
        write(joinpath(d, "circ_scores.tsv"), join([head, row(5, 1.5, true), row(6, 1.2, true)], "\n") * "\n")
        @test read_score_table(joinpath(d, "ell_scores.tsv")) == Dict(4 => 3.0, 5 => 1.0, 6 => 2.0)   # invalid N=7 dropped
        curves = read_score_curves(tmp)
        @test collect(keys(curves)) == ["scanA.sxm"]
        @test curves["scanA.sxm"] == Dict(4 => 3.0, 5 => 1.0, 6 => 1.2)                           # min(ell, circ) per N
        write(joinpath(d, "bad.tsv"), "N\tGCV\n1\t2\n")
        @test_throws ErrorException read_score_table(joinpath(d, "bad.tsv"))
    end
end

@testset "joint count" begin
    # two scans prefer 6 sharply, one prefers 7 weakly: the joint count is 6
    a = Dict(5 => 4.0, 6 => 1.0, 7 => 3.0)
    b = Dict(5 => 5.0, 6 => 1.0, 7 => 2.0)
    c = Dict(6 => 1.1, 7 => 1.0, 8 => 2.0)
    n, J = joint_count([a, b, c])
    @test n == 6
    @test isapprox(J[6], log(1.1))
    @test isapprox(J[8], log(4.0) + log(5.0) + log(2.0))     # a and b lack N=8: each adds its largest log ratio
    @test isapprox(J[5], log(4.0) + log(5.0) + log(2.0))     # c lacks N=5
    # the scale of a curve does not matter, only its shape
    @test joint_count([Dict(k => 1e-6v for (k, v) in a), b, c])[1] == 6
    # ties go to the smaller N; empty curves are ignored; no curve gives 0
    @test joint_count([Dict(5 => 1.0, 6 => 1.0)])[1] == 5
    @test joint_count([Dict{Int,Float64}(), a])[1] == 6
    @test joint_count(Dict{Int,Float64}[])[1] == 0
end

@testset "votes and halves" begin
    @test vote_count([6, 6, 7]) == 6
    @test vote_count([6, 7]) == 6                 # no strict majority: lower median
    @test vote_count([5, 6, 7, 8]) == 6
    @test vote_count([7, 7, 6, 6, 9]) == 7        # no strict majority: median of 6 6 7 7 9
    @test_throws ArgumentError vote_count(Int[])
    @test split_halves(["a", "b", "c", "d", "e"]) == (["a", "c", "e"], ["b", "d"])
end

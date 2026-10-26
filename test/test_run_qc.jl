# Helpers of the label-free run quality page.
#   julia --project=. test/test_run_qc.jl
using Test
include(joinpath(@__DIR__, "..", "scripts", "pipeline", "report_run_qc.jl"))

@testset "votes, split halves and count ranges" begin
    @test vote([6, 6, 7]) == 6 && vote([6, 7]) == 6 && vote([5, 7, 9, 9]) == 7
    n = Dict("a1" => 6, "a2" => 6, "a3" => 6, "a4" => 6, "b1" => 8, "b2" => 10, "b3" => 8, "b4" => 10, "c1" => 5, "c2" => 6)
    members = Dict("A" => ["a1", "a2", "a3", "a4"], "B" => ["b1", "b2", "b3", "b4"], "C" => ["c1", "c2"])
    @test sort(split_half(members, n)) == [0, 2]          # B: odd scans count 8, even scans 10; C has fewer than 4 scans
    @test sort(count_range(members, n)) == [0, 2]         # molecules of at least 3 scans
    @test sort(count_range(members, n; min_scans=2)) == [0, 1, 2]
end

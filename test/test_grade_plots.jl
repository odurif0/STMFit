# Benchmark categories used to sort molecule plots (external grading only).
#   julia --project=. test/test_grade_plots.jl
using Test
include(joinpath(@__DIR__, "grade_consensus_run.jl"))

@testset "plot categories" begin
    grades = [Dict("file" => "a.sxm", "phys_alignment" => "identity"),
              Dict("file" => "dir/b.sxm", "phys_alignment" => "reverse"),
              Dict("file" => "c.sxm", "phys_alignment" => "identity"),
              Dict("file" => "d.sxm", "phys_alignment" => "identity")]
    truth = Dict("a.sxm" => "010010", "b.sxm" => "011000", "c.sxm" => "010010", "d.sxm" => "010010")
    pred = Dict("a.sxm" => "010010", "b.sxm" => "000110", "c.sxm" => "01?000", "d.sxm" => "0100100", "e.sxm" => "0101")
    cats = plot_categories(grades, truth, pred, sort(collect(keys(pred))))
    @test cats["a.sxm"][1] == "1_all_exact" && isempty(cats["a.sxm"][3])
    @test cats["b.sxm"][1] == "1_all_exact"                     # reverse orientation honoured
    @test cats["c.sxm"][1] == "2_N_exact_calls_wrong" && cats["c.sxm"][3] == [3, 5]   # '?' counts as not correct
    @test cats["d.sxm"][1] == "3_N_wrong" && occursin("found 7", cats["d.sxm"][2])
    @test cats["e.sxm"][1] == "4_not_in_benchmark"
end

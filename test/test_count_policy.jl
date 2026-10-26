# Offline replay of the per-scan count layers (campaign K).
#   julia --project=. test/test_count_policy.jl
using Test
include(joinpath(@__DIR__, "..", "scripts", "data", "replay_count_policy.jl"))

P = (spacing_min=0.35, spacing_max=0.75, sigma_max=0.509, max_overlap=0.60, n_min=2, n_max=14, amb=0.05, up=0.30)

@testset "robust-AICc guard" begin
    @test guard(7, 6, false, nothing, NaN, P) == 6                 # down whenever robust AICc wants fewer
    @test guard(6, 7, true, 7, 0.03, P) == 7                       # up by one when GCV is ambiguous
    @test guard(6, 7, true, 7, 0.08, P) == 6                       # not ambiguous enough
    @test guard(6, 8, true, 8, 0.01, P) == 6                       # never up by two
    @test guard(6, nothing, true, 7, 0.01, P) == 6                 # guard failed: N_eff
end

@testset "support midpoint" begin
    # spacing_min_eff = sqrt(-2 log 0.6) * 0.509 = 0.5145; support 3.0 nm: lo = max(2, floor(3/0.75)) = 4,
    # hi = floor(3/0.5145) + 1 = 6, midpoint 5
    @test midpoint(8, 8, 3.0, 0.5, P) == 5                          # two or more above: to the midpoint
    @test midpoint(6, 6, 3.0, 0.5, P) == 5                          # one above: down by one
    @test midpoint(4, 6, 3.0, 0.2, P) == 5                          # below, N_eff higher and close GCV: up by one
    @test midpoint(4, 4, 3.0, 0.2, P) == 4                          # below but N_eff not higher: unchanged
    @test midpoint(5, 5, NaN, 0.2, P) == 5                          # no support: unchanged
end

@testset "variants and consensus" begin
    r = Dict("N_eff" => "7", "robust_aicc_N" => "6", "ambiguous_eff" => "false", "runnerup_N_eff" => "NA",
             "delta_GCV_rel_eff" => "0.4", "support_2D_ell_nm" => "3.0")
    v = variants(r, P)
    @test v == (V0=5, V1=5, V2=6, V3=7)
    @test strict_majority([6, 6, 7]) == 6 && strict_majority([6, 7]) === nothing
    @test vote([6, 7]) == 6 && vote([7, 7, 6]) == 7
end

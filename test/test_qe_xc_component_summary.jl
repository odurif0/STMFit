#!/usr/bin/env julia
using Test, Statistics
include(joinpath(@__DIR__, "summarize_qe_xc_components.jl"))
const M = SummarizeQEXCComponents

@testset "XC small-table summaries and centered triangle bound" begin
    @test M.residual_bound(2., .5, .25, .125) == .4375
    @test isnan(M.residual_bound(0., 0., 0., 0.))
    @test_throws ErrorException M.residual_bound(-1., 0., 0., 0.)
    @test_throws ErrorException M.residual_bound(1., NaN, 0., 0.)
    @test all(isnan, M.extent(Float64[]))
    sd(x) = sqrt(mean(abs2, x .- mean(x)))
    lda = [.1, .2, .4, .3]; local_gga = [.2, -.1, .3, .4]
    divergence = [-1., 2., -3., 4.]; errors = [1e-7, -1e-7, 0., 1e-8]
    xc = lda .+ local_gga .+ divergence .+ errors
    @test sd(xc .- divergence)/sd(xc) <= M.residual_bound(sd(xc), sd(lda), sd(local_gga), maximum(abs, errors))
    counts = [Dict("k"=>string(k), "z_nm"=>string(k/10), "paw_gap"=>string(k>0),
        "lower_half"=>string(k==1), "samples"=>"4", "gga_active"=>string(k==2 ? 4 : 2),
        "negative_valence"=>"1", "negative_xc_input"=>"0", "below_lda"=>"0", "below_gga_density"=>"1") for k in 0:2]
    rows = Dict{String,String}[]
    for r in counts, subset in M.SUBSETS, component in M.COMPONENTS
        empty = r["k"]=="2" && subset=="gga_inactive"
        value = empty ? "NaN" : "0.5"
        push!(rows, merge(r, Dict("subset"=>subset, "component"=>component,
            "samples"=>empty ? "0" : subset=="all" ? "4" : r["k"]=="2" ? "4" : "2",
            "mean"=>value, "minimum"=>value, "maximum"=>value, "spatial_sd"=>value)))
    end
    result = M.summarize_tables("fixture", rows, counts, 0.)
    @test length(result.components)==36 && length(result.branches)==3 && length(result.bounds)==9
    @test [r.planes for r in result.branches]==[3,2,1]
    @test all(r -> r.negative_valence_fraction_min==r.negative_valence_fraction_max==.25, result.branches)
    @test all(r -> r.mean_min_ev==r.mean_max_ev==.5M.B.S.RY_EV, result.components)
    @test all(r -> r.centered_residual_relative_l2_upper_bound_min==r.centered_residual_relative_l2_upper_bound_max==2, result.bounds)
    inactive = only(filter(r -> r.domain=="paw_gap" && r.subset=="gga_inactive", result.bounds))
    @test inactive.defined_planes==1 && inactive.empty_planes==1 && inactive.zero_xc_sd_planes==0
    zerorows = deepcopy(rows)
    for r in zerorows
        r["component"]=="xc_total" && r["samples"]!="0" && (r["spatial_sd"]="0")
    end
    zeroresult = M.summarize_tables("fixture", zerorows, counts, 0.)
    @test all(r -> r.defined_planes==0 && r.zero_xc_sd_planes+r.empty_planes==r.planes, zeroresult.bounds)
    @test all(r -> isnan(r.centered_residual_relative_l2_upper_bound_min) && isnan(r.centered_residual_relative_l2_upper_bound_max), zeroresult.bounds)
    @test_throws ErrorException M.summarize_tables("fixture", rows[2:end], counts, 0.)
    @test_throws ErrorException M.summarize_tables("fixture", rows, [counts; counts[1]], 0.)
end

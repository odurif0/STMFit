#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "audit_joint_proxy_mold_support.jl"))

@testset "joint-proxy mold support audit" begin
    coords = [(t, u) for u in (-1.0, 0.0, 1.0) for t in (-1.0, 0.0, 1.0)]
    centered = Dict(coord => (coord == (0.0, 0.0) ? 2.0 : 0.0) for coord in coords)
    centered_metrics = support_metrics("centered", centered)
    @test centered_metrics.edge_l2_fraction == 0.0
    @test !centered_metrics.peak_on_edge
    @test (centered_metrics.peak_t_nm, centered_metrics.peak_u_nm) == (0.0, 0.0)

    edged = Dict(coord => (coord == (1.0, 1.0) ? -3.0 : 0.0) for coord in coords)
    edge_metrics = support_metrics("edge", edged)
    @test edge_metrics.edge_l2_fraction == 1.0
    @test edge_metrics.edge_absmax_fraction == 1.0
    @test edge_metrics.peak_on_edge

    maps = Dict(0 => centered, 1 => edged)
    rows = audit_map_support(maps)
    @test Set(row.name for row in rows) == Set(["type_0", "type_1", "type_1_minus_type_0"])
    diff = only(row for row in rows if row.name == "type_1_minus_type_0")
    @test diff.peak_on_edge

    mktempdir() do directory
        maps_path = joinpath(directory, "maps.tsv")
        out_path = joinpath(directory, "audit.tsv")
        open(maps_path, "w") do io
            println(io, "type\tt_nm\tu_nm\tvalue")
            for typ in (0, 1), coordinate in coords
                println(io, join((typ, coordinate..., maps[typ][coordinate]), '\t'))
            end
        end
        parsed = read_typed_maps(maps_path)
        @test parsed == maps
        main(["--maps", maps_path, "--out", out_path])
        output = readlines(out_path)
        @test length(output) == 4
        @test startswith(output[1], "name\tn_points")
    end
end

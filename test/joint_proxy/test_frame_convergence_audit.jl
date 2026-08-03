#!/usr/bin/env julia

using Test

include(joinpath(dirname(@__DIR__), "audit_joint_proxy_frame_convergence.jl"))

function _write_map(path, half; decay=0.5)
    coords = collect(-half:1.0:half)
    open(path, "w") do io
        println(io, "type\tt_nm\tu_nm\tvalue")
        for typ in (0, 1), u in coords, t in coords
            r = sqrt(t^2 + u^2)
            value = (typ == 0 ? 1.0 : 2.0) * exp(-r / decay)
            println(io, join((typ, t, u, value), '\t'))
        end
    end
end

@testset "joint-proxy frame convergence audit" begin
    mktempdir() do dir
        ref = joinpath(dir, "ref.tsv")
        wide = joinpath(dir, "wide.tsv")
        _write_map(ref, 1.0)
        _write_map(wide, 2.0)

        rows = audit_frame_convergence(ref, [wide])
        @test length(rows) == 1
        @test rows[1].grid_n == 5
        @test rows[1].central_max_abs_difference == 0.0
        @test rows[1].finite_fraction == 1.0
        @test rows[1].edge_l2_fraction_type0 < 1.0

        out = joinpath(dir, "convergence.tsv")
        write_convergence_audit(out, rows)
        lines = readlines(out)
        @test length(lines) == 2
        @test startswith(lines[1], "grid_n")
    end
end

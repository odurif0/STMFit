#!/usr/bin/env julia

using Test
using SHA

include(joinpath(dirname(@__DIR__), "audit_joint_proxy_frame_compatibility.jl"))

function write_map(path, half; central=nothing, nonfinite_edge=false)
    coords = collect(-half:1.0:half)
    open(path, "w") do io
        println(io, "type\tt_nm\tu_nm\tvalue")
        for typ in (0, 1), u in coords, t in coords
            value = central !== nothing && abs(t) <= 1 && abs(u) <= 1 ?
                central[(typ, t, u)] : typ + 3t - 2u
            nonfinite_edge && typ == 1 && t == half && u == half && (value = "NA")
            println(io, join((typ, t, u, value), '\t'))
        end
    end
end

function write_provenance(path, maps_path; half, glcn_hash="a"^64)
    maps_hash = bytes2hex(open(SHA.sha256, maps_path))
    open(path, "w") do io
        println(io, "schema = \"stmfit-qe-mold-provenance-v1\"")
        println(io, "provider = \"test\"")
        println(io, "maps_path = \"$maps_path\"")
        println(io, "maps_sha256 = \"$maps_hash\"")
        println(io, "glcn_cube_path = \"missing-glcn.cube\"")
        println(io, "glcn_cube_sha256 = \"$glcn_hash\"")
        println(io, "glcnac_cube_path = \"missing-glcnac.cube\"")
        println(io, "glcnac_cube_sha256 = \"$("b"^64)\"")
        println(io, "half_nm = $half")
        println(io, "step_nm = 1.0")
        println(io, "height_nm = 0.5")
        println(io, "sample_bias_ev = -0.3")
        println(io, "cube_units = \"bohr\"")
    end
end

@testset "joint-proxy frame compatibility audit" begin
    mktempdir() do directory
        reference_map = joinpath(directory, "reference.tsv")
        wider_map = joinpath(directory, "wider.tsv")
        reference_provenance = joinpath(directory, "reference.toml")
        wider_provenance = joinpath(directory, "wider.toml")
        write_map(reference_map, 1.0)
        central = read_typed_map(reference_map).values
        write_map(wider_map, 2.0; central=central)
        write_provenance(reference_provenance, reference_map; half=1.0)
        write_provenance(wider_provenance, wider_map; half=2.0)

        compatible = audit_frame_compatibility(reference_map, reference_provenance,
            wider_map, wider_provenance)
        @test compatible.compatible
        @test compatible.reason == "compatible_nested_frame"
        @test compatible.reference.grid_n == 3
        @test compatible.candidate.grid_n == 5
        @test compatible.central_max_abs_difference == 0.0
        @test !compatible.reference.cubes_present

        missing = audit_frame_compatibility(reference_map, reference_provenance,
            wider_map, nothing)
        @test !missing.compatible
        @test missing.reason == "candidate_provenance_missing"

        mismatched_provenance = joinpath(directory, "mismatched.toml")
        write_provenance(mismatched_provenance, wider_map; half=2.0, glcn_hash="c"^64)
        mismatched = audit_frame_compatibility(reference_map, reference_provenance,
            wider_map, mismatched_provenance)
        @test !mismatched.compatible
        @test mismatched.reason == "cube_hash_mismatch"

        nonfinite_map = joinpath(directory, "nonfinite.tsv")
        nonfinite_provenance = joinpath(directory, "nonfinite.toml")
        write_map(nonfinite_map, 2.0; central=central, nonfinite_edge=true)
        write_provenance(nonfinite_provenance, nonfinite_map; half=2.0)
        nonfinite = audit_frame_compatibility(reference_map, reference_provenance,
            nonfinite_map, nonfinite_provenance)
        @test !nonfinite.compatible
        @test nonfinite.reason == "nonfinite_map_values"

        out = joinpath(directory, "audit.tsv")
        write_frame_compatibility(out, missing)
        @test length(readlines(out)) == 2
    end
end

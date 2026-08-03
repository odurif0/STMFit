using Test
using TOML
using SHA

include(joinpath(@__DIR__, "lib", "qe_mold_provenance.jl"))
using .QEMoldProvenance: write_qe_mold_provenance

_sha256_bytes(path::AbstractString) = open(path, "r") do io
    SHA.sha256(io)
end
_sha256_hex(path::AbstractString) = bytes2hex(_sha256_bytes(path))

# Exact v1 (constant-height) key set. The writer must remain byte-identical to
# the frozen stm_dft_v1 sidecar when no observable is requested.
const V1_KEYS = Set([
    "schema", "provider", "sample_bias_ev", "height_nm", "half_nm", "step_nm",
    "cube_units", "glcn_cube_path", "glcnac_cube_path", "maps_path",
    "templates_path", "glcn_cube_sha256", "glcnac_cube_sha256",
    "maps_sha256", "templates_sha256",
])

@testset "QE mold provenance binds physical inputs and generated outputs" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "maps.tsv")
        templates = joinpath(tmp, "templates.tsv")
        out = joinpath(tmp, "provenance.toml")
        write(glcn, "glcn\n")
        write(glcnac, "glcnac\n")
        write(maps, "maps\n")
        write(templates, "templates\n")

        write_qe_mold_provenance(out; provider="stm_dft_v1", glcn_cube=glcn,
            glcnac_cube=glcnac, maps=maps, templates=templates,
            sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr")

        manifest = TOML.parsefile(out)
        @test manifest["schema"] == "stmfit-qe-mold-provenance-v1"
        @test manifest["provider"] == "stm_dft_v1"
        @test manifest["sample_bias_ev"] == -0.3
        @test manifest["height_nm"] == 0.5
        @test all(occursin(r"^[0-9a-f]{64}$", manifest[key]) for key in
                  ("glcn_cube_sha256", "glcnac_cube_sha256", "maps_sha256", "templates_sha256"))
    end
end

@testset "constant-height provenance is byte-identical to v1 (no observable leak)" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "maps.tsv")
        templates = joinpath(tmp, "templates.tsv")
        write(glcn, "glcn\n"); write(glcnac, "glcnac\n")
        write(maps, "maps\n"); write(templates, "templates\n")
        out = joinpath(tmp, "prov_ch.toml")
        write_qe_mold_provenance(out; provider="stm_dft_v1", glcn_cube=glcn,
            glcnac_cube=glcnac, maps=maps, templates=templates,
            sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr")
        manifest = TOML.parsefile(out)
        # No observable/cc fields leak into the constant-height schema.
        @test Set(keys(manifest)) == V1_KEYS
        @test !haskey(manifest, "observable")
        @test !haskey(manifest, "isovalue")
        @test !haskey(manifest, "type_isovalues")
        @test !haskey(manifest, "invalid_mask_sha256")
        # Byte-stable across two writes (deterministic sorted TOML).
        out2 = joinpath(tmp, "prov_ch2.toml")
        write_qe_mold_provenance(out2; provider="stm_dft_v1", glcn_cube=glcn,
            glcnac_cube=glcnac, maps=maps, templates=templates,
            sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr")
        @test read(out, String) == read(out2, String)
    end
end

@testset "constant-current provenance binds the full observable chain" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        write(glcn, "glcn-bytes\n"); write(glcnac, "glcnac-bytes\n")
        write(maps, "cc-map\n"); write(mask, "cc-mask\n")
        bracket_artifacts = Tuple{Float64,String,String}[]
        for height in (0.40, 0.45, 0.55, 0.60)
            bracket_map = joinpath(tmp, "cc_map_$(height).tsv")
            bracket_mask = joinpath(tmp, "cc_map_$(height).mask.tsv")
            write(bracket_map, "map $(height)\n")
            write(bracket_mask, "mask $(height)\n")
            push!(bracket_artifacts, (height, bracket_map, bracket_mask))
        end
        out = joinpath(tmp, "cc_prov.toml")
        write_qe_mold_provenance(out;
            provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr",
            observable="constant-current",
            nominal_height_nm=0.50,
            bracket_heights_nm=[0.40, 0.45, 0.50, 0.55, 0.60],
            isovalue=1.234e-4,
            z_spacing_nm=0.0212,
            crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=mask,
            frame_origin_nm=[0.0, 0.0, 0.0],
            frame_t_axis=[1.0, 0.0, 0.0],
            frame_u_axis=[0.0, 1.0, 0.0],
            bracket_artifacts=bracket_artifacts)
        m = TOML.parsefile(out)
        @test m["schema"] == "stmfit-qe-mold-provenance-v1"
        @test m["observable"] == "constant-current"
        @test m["provider"] == "stm_dft_cc_diag"
        @test m["nominal_height_nm"] == 0.50
        @test m["bracket_heights_nm"] == [0.40, 0.45, 0.50, 0.55, 0.60]
        @test !haskey(m, "isovalue")
        isos = Dict(item["type"] => item for item in m["type_isovalues"])
        @test Set(keys(isos)) == Set([0, 1])
        @test isos[0]["isovalue"] == 1.234e-4
        @test isos[1]["isovalue"] == 1.234e-4
        @test m["z_spacing_nm"] == 0.0212
        @test m["crossing_policy"] == "vacuum_first_bracketed_linear"
        @test !haskey(m, "reference_plane")
        @test length(m["type_frames"]) == 2
        frames = Dict(frame["type"] => frame for frame in m["type_frames"])
        @test Set(keys(frames)) == Set([0, 1])
        @test frames[0]["origin_nm"] == [0.0, 0.0, 0.0]
        @test frames[1]["origin_nm"] == [0.0, 0.0, 0.0]
        @test frames[0]["t_axis"] == [1.0, 0.0, 0.0]
        @test frames[1]["u_axis"] == [0.0, 1.0, 0.0]
        @test m["invalid_mask_sha256"] == _sha256_hex(mask)
        @test m["maps_sha256"] == _sha256_hex(maps)
        @test m["glcn_cube_sha256"] == _sha256_hex(glcn)
        @test m["glcnac_cube_sha256"] == _sha256_hex(glcnac)
        @test occursin(r"^[0-9a-f]{64}$", m["invalid_mask_sha256"])
        # No templates block: the constant-current build emits maps+mask only.
        @test !haskey(m, "templates_path")
        @test !haskey(m, "templates_sha256")
    end
end

@testset "constant-current provenance binds distinct typed nominal isovalues" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        for path in (glcn, glcnac, maps, mask)
            write(path, "x\n")
        end
        out = joinpath(tmp, "cc_typed_iso.toml")
        write_qe_mold_provenance(out;
            provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50],
            type_isovalues=[
                (type=0, isovalue=1.0e-4, source="calibrated"),
                (type=1, isovalue=2.0e-4, source="calibrated"),
            ],
            z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear", invalid_mask=mask,
            type_frames=[
                (type=0, origin_nm=[0.0, 0.0, 0.0],
                 t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0]),
                (type=1, origin_nm=[0.1, 0.2, 0.3],
                 t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0]),
            ],
            bracket_artifacts=Tuple{Float64,String,String}[])

        m = TOML.parsefile(out)
        @test !haskey(m, "isovalue")
        isos = Dict(item["type"] => item for item in m["type_isovalues"])
        @test isos[0]["isovalue"] == 1.0e-4
        @test isos[1]["isovalue"] == 2.0e-4
        @test isos[0]["source"] == "calibrated"
        @test isos[1]["source"] == "calibrated"
    end
end

@testset "constant-current provenance binds distinct typed frames" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        low_map = joinpath(tmp, "low.tsv")
        low_mask = joinpath(tmp, "low.mask.tsv")
        high_map = joinpath(tmp, "high.tsv")
        high_mask = joinpath(tmp, "high.mask.tsv")
        for path in (glcn, glcnac, maps, mask, low_map, low_mask, high_map, high_mask)
            write(path, "x\n")
        end
        out = joinpath(tmp, "cc_typed_prov.toml")
        write_qe_mold_provenance(out;
            provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr",
            observable="constant-current",
            nominal_height_nm=0.50,
            bracket_heights_nm=[0.40, 0.50, 0.60],
            isovalue=1.0e-4,
            z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=mask,
            type_frames=[
                (type=0, origin_nm=[1.0, 2.0, 3.0],
                 t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0]),
                (type=1, origin_nm=[1.1, 2.2, 3.3],
                 t_axis=[0.0, 1.0, 0.0], u_axis=[0.0, 0.0, 1.0]),
            ],
            bracket_artifacts=[(0.40, low_map, low_mask), (0.60, high_map, high_mask)])

        m = TOML.parsefile(out)
        @test !haskey(m, "reference_plane")
        @test length(m["type_frames"]) == 2
        frames = Dict(frame["type"] => frame for frame in m["type_frames"])
        @test Set(keys(frames)) == Set([0, 1])
        @test frames[0]["origin_nm"] == [1.0, 2.0, 3.0]
        @test frames[1]["origin_nm"] == [1.1, 2.2, 3.3]
        @test frames[0]["t_axis"] == [1.0, 0.0, 0.0]
        @test frames[1]["u_axis"] == [0.0, 0.0, 1.0]
    end
end

@testset "constant-current provenance rejects incomplete typed frames" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        for path in (glcn, glcnac, maps, mask)
            write(path, "x\n")
        end
        base_kw = (provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50],
            isovalue=1.0e-4, z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear", invalid_mask=mask,
            bracket_artifacts=Tuple{Float64,String,String}[])
        base_no_iso = NamedTuple{filter(!=(:isovalue), keys(base_kw))}(base_kw)
        good0 = (type=0, origin_nm=[0.0, 0.0, 0.0],
            t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0])
        good1 = (type=1, origin_nm=[0.1, 0.2, 0.3],
            t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0])
        malformed = (type=1, origin_nm=[0.1, 0.2],
            t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0])

        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "one.toml");
            base_kw..., type_frames=[good0])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "dup.toml");
            base_kw..., type_frames=[good0, good0])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "badtype.toml");
            base_kw..., type_frames=[good0, merge(good1, (type=2,))])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "malformed.toml");
            base_kw..., type_frames=[good0, malformed])
    end
end

@testset "constant-current provenance rejects Bool numeric fields" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        for path in (glcn, glcnac, maps, mask)
            write(path, "x\n")
        end
        base_kw = (provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50],
            isovalue=1.0e-4, z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear", invalid_mask=mask,
            bracket_artifacts=Tuple{Float64,String,String}[])
        base_no_iso = NamedTuple{filter(!=(:isovalue), keys(base_kw))}(base_kw)
        good0 = (type=0, origin_nm=[0.0, 0.0, 0.0],
            t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0])
        good1 = (type=1, origin_nm=[0.1, 0.2, 0.3],
            t_axis=[1.0, 0.0, 0.0], u_axis=[0.0, 1.0, 0.0])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "bool_type.toml");
            base_kw..., type_frames=[merge(good0, (type=false,)), good1])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "bool_frame.toml");
            base_kw..., type_frames=[merge(good0, (origin_nm=Any[false, 0.0, 0.0],)), good1])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "bool_iso.toml");
            base_no_iso..., isovalue=true, type_frames=[good0, good1])
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "bool_typed_iso.toml");
            base_no_iso...,
            type_isovalues=[(type=0, isovalue=true), (type=1, isovalue=1.0e-4)],
            type_frames=[good0, good1])
    end
end

@testset "constant-current provenance rejects unbound bracket heights" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube"); glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv"); mask = joinpath(tmp, "cc.mask.tsv")
        write(glcn, "x\n"); write(glcnac, "x\n"); write(maps, "x\n"); write(mask, "x\n")
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "x.toml");
            provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.40, 0.50, 0.60],
            isovalue=1.0e-4, z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear", invalid_mask=mask,
            frame_origin_nm=[0.0,0.0,0.0], frame_t_axis=[1.0,0.0,0.0],
            frame_u_axis=[0.0,1.0,0.0])
    end
end

@testset "provenance atomically replaces an output symlink without touching its target" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube"); glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "maps.tsv"); templates = joinpath(tmp, "templates.tsv")
        victim = joinpath(tmp, "victim.toml"); out = joinpath(tmp, "provenance.toml")
        write(glcn, "x\n"); write(glcnac, "x\n"); write(maps, "x\n")
        write(templates, "x\n"); write(victim, "sentinel\n"); symlink(victim, out)

        write_qe_mold_provenance(out; provider="stm_dft_v1", glcn_cube=glcn,
            glcnac_cube=glcnac, maps=maps, templates=templates,
            sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr")

        @test read(victim, String) == "sentinel\n"
        @test !islink(out)
        @test TOML.parsefile(out)["provider"] == "stm_dft_v1"
    end
end

@testset "constant-current provenance requires every observable field" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube")
        glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        mask = joinpath(tmp, "cc_maps.mask.tsv")
        write(glcn, "x\n"); write(glcnac, "x\n")
        write(maps, "x\n"); write(mask, "x\n")
        low_map = joinpath(tmp, "low.tsv"); low_mask = joinpath(tmp, "low.mask.tsv")
        high_map = joinpath(tmp, "high.tsv"); high_mask = joinpath(tmp, "high.mask.tsv")
        for path in (low_map, low_mask, high_map, high_mask); write(path, "x\n"); end
        base_kw = (provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.40, 0.50, 0.60],
            isovalue=1.0e-4, z_spacing_nm=0.02,
            crossing_policy="vacuum_first_bracketed_linear", invalid_mask=mask,
            frame_origin_nm=[0.0,0.0,0.0], frame_t_axis=[1.0,0.0,0.0],
            frame_u_axis=[0.0,1.0,0.0],
            bracket_artifacts=[(0.40, low_map, low_mask), (0.60, high_map, high_mask)])
        # Complete call succeeds.
        out = joinpath(tmp, "ok.toml")
        write_qe_mold_provenance(out; base_kw...)
        @test TOML.parsefile(out)["observable"] == "constant-current"
        # Dropping any required cc field must error with a named diagnostic.
        for drop in (:nominal_height_nm, :bracket_heights_nm, :isovalue,
                     :z_spacing_nm, :crossing_policy, :invalid_mask,
                     :frame_origin_nm, :frame_t_axis, :frame_u_axis)
            keep = filter(!=(drop), keys(base_kw))
            reduced = NamedTuple{keep}(base_kw)
            @test_throws ErrorException write_qe_mold_provenance(
                joinpath(tmp, "bad.toml"); reduced...)
        end
    end
end

@testset "constant-current provenance: tampered mask changes recorded hash" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube"); glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv"); mask = joinpath(tmp, "cc.mask.tsv")
        write(glcn, "x\n"); write(glcnac, "x\n"); write(maps, "x\n"); write(mask, "v1\n")
        cc_kw = (provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50], isovalue=1.0e-4,
            z_spacing_nm=0.02, crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=mask, frame_origin_nm=[0.0,0.0,0.0],
            frame_t_axis=[1.0,0.0,0.0], frame_u_axis=[0.0,1.0,0.0],
            bracket_artifacts=Tuple{Float64,String,String}[])
        out = joinpath(tmp, "p.toml")
        write_qe_mold_provenance(out; cc_kw...)
        h1 = TOML.parsefile(out)["invalid_mask_sha256"]
        @test h1 == _sha256_hex(mask)
        # Tamper the mask; the recomputed provenance must record a different hash.
        write(mask, "v2-tampered\n")
        write_qe_mold_provenance(out; cc_kw...)
        h2 = TOML.parsefile(out)["invalid_mask_sha256"]
        @test h2 == _sha256_hex(mask)
        @test h1 != h2
    end
end

@testset "constant-current provenance: missing mask file errors cleanly" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube"); glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv")
        write(glcn, "x\n"); write(glcnac, "x\n"); write(maps, "x\n")
        @test_throws ErrorException write_qe_mold_provenance(joinpath(tmp, "x.toml");
            provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50], isovalue=1.0e-4,
            z_spacing_nm=0.02, crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=joinpath(tmp, "missing.mask.tsv"),
            frame_origin_nm=[0.0,0.0,0.0], frame_t_axis=[1.0,0.0,0.0],
            frame_u_axis=[0.0,1.0,0.0],
            bracket_artifacts=Tuple{Float64,String,String}[])
    end
end

@testset "constant-current provenance binds isovalue scan policy and hash" begin
    mktempdir() do tmp
        glcn = joinpath(tmp, "glcn.cube"); glcnac = joinpath(tmp, "glcnac.cube")
        maps = joinpath(tmp, "cc_maps.tsv"); mask = joinpath(tmp, "cc.mask.tsv")
        for path in (glcn, glcnac, maps, mask)
            write(path, "x\n")
        end
        base_kw = (provider="stm_dft_cc_diag", glcn_cube=glcn, glcnac_cube=glcnac,
            maps=maps, sample_bias_ev=-0.3, height_nm=0.5, half_nm=0.32,
            step_nm=0.08, cube_units="bohr", observable="constant-current",
            nominal_height_nm=0.50, bracket_heights_nm=[0.50], isovalue=1.0e-4,
            z_spacing_nm=0.02, crossing_policy="vacuum_first_bracketed_linear",
            invalid_mask=mask, frame_origin_nm=[0.0,0.0,0.0],
            frame_t_axis=[1.0,0.0,0.0], frame_u_axis=[0.0,1.0,0.0],
            bracket_artifacts=Tuple{Float64,String,String}[])
        default_out = joinpath(tmp, "default.toml")
        explicit_out = joinpath(tmp, "explicit.toml")
        changed_out = joinpath(tmp, "changed.toml")
        write_qe_mold_provenance(default_out; base_kw...)
        write_qe_mold_provenance(explicit_out; base_kw..., isovalue_scan_intervals=1024)
        write_qe_mold_provenance(changed_out; base_kw..., isovalue_scan_intervals=2048)
        @test TOML.parsefile(default_out)["isovalue_scan_intervals"] == 1024
        @test read(default_out) == read(explicit_out)
        @test TOML.parsefile(changed_out)["isovalue_scan_intervals"] == 2048
        @test _sha256_hex(default_out) != _sha256_hex(changed_out)
        @test_throws ErrorException write_qe_mold_provenance(
            joinpath(tmp, "invalid.toml"); base_kw..., isovalue_scan_intervals=0)
    end
end

@testset "constant-current diagnostic config preserves frozen registration" begin
    diagnostic_path = joinpath(@__DIR__, "..", "config",
        "joint_proxy_whole_roi_diagnostic.toml")
    constant_current_path = joinpath(@__DIR__, "..", "config",
        "joint_proxy_whole_roi_constant_current.toml")
    diagnostic = TOML.parsefile(diagnostic_path)
    constant_current = TOML.parsefile(constant_current_path)
    @test constant_current["registration"] == diagnostic["registration"]
    @test constant_current["numerical_tie"] == diagnostic["numerical_tie"]
    @test Set(keys(constant_current)) == Set(["registration", "numerical_tie", "observable"])

    observable = constant_current["observable"]
    @test observable["mode"] == "constant-current"
    @test observable["provider"] == "stm_dft_cc_diag"
    @test observable["sample_bias_ev"] == -0.300
    @test observable["nominal_height_nm"] == 0.50
    @test observable["bracket_heights_nm"] == [0.40, 0.45, 0.50, 0.55, 0.60]
    @test observable["bracket_policy"] == "diagnostics_only"
    @test observable["crossing_policy"] == "vacuum_first_bracketed_linear"
    @test observable["provenance_schema"] == "stmfit-qe-mold-provenance-v1"
    @test observable["maps_sha256_key"] == "maps_sha256"
    @test observable["invalid_mask_sha256_key"] == "invalid_mask_sha256"
    @test observable["bracket_artifacts_key"] == "bracket_artifacts"
    @test all(occursin(r"^[0-9a-f]{64}$", observable[key]) for key in
        ("glcn_cube_sha256", "glcnac_cube_sha256"))
end

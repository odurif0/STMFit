#!/usr/bin/env julia

# End-to-end smoke test for the synthetic Quantum ESPRESSO -> STM mold helper
# workflow. This does not run QE; it validates the local file handoffs and formats.

using Printf
using Test

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_vec3, _read_key_tsv

const ROOT = normpath(joinpath(@__DIR__, ".."))

function _run_script(script::String, args::Vector{String})
    project = "--project=$(ROOT)"
    run(`$(Base.julia_cmd()) $project $(joinpath(ROOT, script)) $args`)
end

function _write_synthetic_trimer(path::String)
    open(path, "w") do io
        println(io, "4")
        println(io, "synthetic oriented trimer proxy for smoke testing")
        println(io, "C  0.0  0.0  0.0")
        println(io, "C  1.0  0.0  0.0")
        println(io, "O  1.0  1.0  0.0")
        println(io, "H  0.0  0.0  1.0")
    end
end

function _parse_vec3_checked(s::String)
    vals = _parse_vec3(s)
    @test length(vals) == 3
    return vals
end

function _read_qe_positions(path::String)
    lines = readlines(path)
    idx = findfirst(line -> startswith(strip(line), "ATOMIC_POSITIONS"), lines)
    idx === nothing && error("No ATOMIC_POSITIONS in $path")
    atoms = Tuple{String,Float64,Float64,Float64}[]
    for line in lines[(idx+1):end]
        s = strip(line)
        isempty(s) && break
        startswith(s, "K_POINTS") && break
        parts = split(s)
        length(parts) >= 4 || break
        push!(atoms, (parts[1], parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])))
    end
    isempty(atoms) && error("No atom rows in $path")
    return atoms
end

function _write_fake_qe_relax(path::String, metadata::String, atoms)
    md = _read_key_tsv(metadata)
    cell = [_parse_vec3_checked(md[k]) for k in ("cell_a", "cell_b", "cell_c")]
    open(path, "w") do io
        println(io, "synthetic QE output")
        println(io, "CELL_PARAMETERS (angstrom)")
        for v in cell
            println(io, @sprintf("  %.10f %.10f %.10f", v[1], v[2], v[3]))
        end
        println(io, "ATOMIC_POSITIONS (angstrom)")
        for (elem, x, y, z) in atoms
            println(io, @sprintf("%s %.10f %.10f %.10f", elem, x, y, z + 0.01))
        end
    end
end

function _write_cube(path::String; offset::Float64=0.0)
    n = 12
    step = 0.1
    vals = Float64[]
    for ix in 0:(n-1), iy in 0:(n-1), iz in 0:(n-1)
        x = ix * step
        y = iy * step
        z = iz * step
        push!(vals, offset + x + 2y + 3z)
    end
    open(path, "w") do io
        println(io, "synthetic cube")
        println(io, "value=x+2y+3z")
        println(io, "0 0.0 0.0 0.0")
        println(io, @sprintf("%d %.10f 0.0 0.0", n, step))
        println(io, @sprintf("%d 0.0 %.10f 0.0", n, step))
        println(io, @sprintf("%d 0.0 0.0 %.10f", n, step))
        for (i, v) in enumerate(vals)
            print(io, @sprintf(" %.8e", v))
            i % 6 == 0 && println(io)
        end
        length(vals) % 6 == 0 || println(io)
    end
end

function _n_data_rows(path::String)
    return count(line -> !isempty(strip(line)), readlines(path)) - 1
end

@testset "QE STM mold helper smoke test" begin
    mktempdir() do dir
        molecule = joinpath(dir, "trimer.xyz")
        slab = joinpath(dir, "trimer_slab.xyz")
        slab_meta = joinpath(dir, "trimer_slab_meta.tsv")
        qe_dir = joinpath(dir, "qe")
        fake_relax = joinpath(qe_dir, "smoke_relax.out")
        relaxed_xyz = joinpath(qe_dir, "smoke_relaxed.xyz")
        relaxed_meta = joinpath(qe_dir, "smoke_relaxed_meta.tsv")
        scf_relaxed = joinpath(qe_dir, "pw_scf_relaxed.in")
        frame = joinpath(qe_dir, "frame.tsv")
        cube0 = joinpath(dir, "type0.cube")
        cube1 = joinpath(dir, "type1.cube")
        maps = joinpath(dir, "stm_maps.tsv")
        molds = joinpath(dir, "molds.tsv")
        bonds = joinpath(dir, "bonds.tsv")

        _write_synthetic_trimer(molecule)
        _run_script("test/build_qe_slab_trimer_xyz.jl", [
            "--molecule", molecule,
            "--out", slab,
            "--metadata", slab_meta,
            "--nx", "3",
            "--ny", "3",
            "--layers", "2",
            "--height-above-top", "2.6",
            "--vacuum", "8.0",
        ])
        @test isfile(slab)
        @test isfile(slab_meta)

        _run_script("test/prepare_qe_mold_inputs.jl", [
            "--xyz", slab,
            "--cell-metadata", slab_meta,
            "--out-dir", qe_dir,
            "--prefix", "smoke",
            "--fix-below-z", "1.0",
            "--emin-ev", "-0.5",
            "--emax-ev", "0.1",
            "--ntasks", "2",
            "--walltime", "00:10:00",
        ])
        @test all(isfile, joinpath.(Ref(qe_dir), ["pw_relax.in", "pw_scf.in", "pp_ldos.in", "run_qe_mold.sbatch"]))
        @test occursin("CELL_PARAMETERS angstrom", read(joinpath(qe_dir, "pw_relax.in"), String))

        atoms = _read_qe_positions(joinpath(qe_dir, "pw_scf.in"))
        _write_fake_qe_relax(fake_relax, slab_meta, atoms)
        _run_script("test/extract_qe_relaxed_xyz.jl", [
            "--qe-out", fake_relax,
            "--out", relaxed_xyz,
            "--metadata", relaxed_meta,
        ])
        @test isfile(relaxed_xyz)
        @test isfile(relaxed_meta)
        @test occursin("cell_found\ttrue", read(relaxed_meta, String))

        _run_script("test/update_qe_positions_from_xyz.jl", [
            "--input", joinpath(qe_dir, "pw_scf.in"),
            "--xyz", relaxed_xyz,
            "--out", scf_relaxed,
        ])
        @test occursin("ATOMIC_POSITIONS angstrom", read(scf_relaxed, String))

        n_slab = 3 * 3 * 2
        _run_script("test/extract_qe_mold_frame.jl", [
            "--xyz", relaxed_xyz,
            "--origin-indices", "$(n_slab + 1),$(n_slab + 2)",
            "--axis-from", string(n_slab + 1),
            "--axis-to", string(n_slab + 2),
            "--plane-index", string(n_slab + 3),
            "--height-nm", "0.05",
            "--out", frame,
        ])
        @test occursin("height_nm\t0.05", read(frame, String))

        _write_cube(cube0; offset=0.0)
        _write_cube(cube1; offset=1.0)
        _run_script("test/cube_to_stm_maps.jl", [
            "--cube", "0:$(cube0)",
            "--frame", "0:$(frame)",
            "--cube", "1:$(cube1)",
            "--frame", "1:$(frame)",
            "--cube-units", "nm",
            "--half-nm", "0.16",
            "--step-nm", "0.08",
            "--out", maps,
        ])
        @test _n_data_rows(maps) == 2 * 5 * 5

        _run_script("test/import_stm_mold_maps.jl", [
            "--maps", maps,
            "--out", molds,
            "--bond-out", bonds,
            "--half-nm", "0.16",
            "--step-nm", "0.08",
            "--normalize", "none",
        ])
        @test _n_data_rows(molds) == 8
        @test _n_data_rows(bonds) == 16
    end
end

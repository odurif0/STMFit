#!/usr/bin/env julia
# Julia 1.13 focused tests. Python is used ONLY here as a numerical reference.
# Synthetic: julia --project=. test/test_cc_mold_native.jl
# One real pair (no STM fit or labels):
#   julia --project=. test/test_cc_mold_native.jl --real-cubes \
#     CUBE0 CUBE1 FRAME0 FRAME1 NEW_OUT_DIR
using Test
using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics
using TOML

include(joinpath(@__DIR__, "lib", "cc_mold_native.jl"))
using .CCMoldNative

const ROOT = dirname(@__DIR__)
const CONFIG = joinpath(ROOT, "config", "unit_assignment_reconstructed.toml")
const SETTINGS = load_settings(CONFIG)

function changed_settings(; kwargs...)
    names = fieldnames(MoldSettings)
    values = NamedTuple{names}(Tuple(getfield(SETTINGS, n) for n in names))
    return MoldSettings(merge(values, (; kwargs...))...)
end

function write_cube(path, dims, origin, axes, density; signed_axes=false)
    open(path, "w") do io
        println(io, "synthetic scalar LDOS\nlegacy token-order fixture")
        @printf(io, "0 %.17g %.17g %.17g\n", (origin ./ CCMoldNative.BOHR_NM)...)
        for a in 1:3
            @printf(io, "%d %.17g %.17g %.17g\n", signed_axes ? -dims[a] : dims[a],
                    (axes[:, a] ./ CCMoldNative.BOHR_NM)...)
        end
        # First axis fast, deliberately testing the reference convention.
        for k in 0:(dims[3]-1), j in 0:(dims[2]-1), i in 0:(dims[1]-1)
            text = @sprintf("%.17E", density(i, j, k))
            println(io, replace(text, 'E' => 'D'))
        end
    end
    return path
end

function write_frame(path, origin, t, u)
    open(path, "w") do io
        println(io, "field\tvalue")
        println(io, "origin_nm\t", join(origin, ','))
        println(io, "t_axis\t", join(t, ','))
        println(io, "u_axis\t", join(u, ','))
        println(io, "height_nm\t0.55") # not a replacement for configured target
        println(io, "normal_axis\tignored reference metadata")
    end
    return path
end

function synthetic_pair(dir)
    dims = (21, 23, 71)
    origin = [-0.5, -0.55, -0.6]
    axes = Matrix(Diagonal([0.05, 0.05, 0.05]))
    cubes = String[]
    for typ in 0:1
        density = function(i, j, k)
            x, y, z = origin + axes * [i, j, k]
            x > 0.15 && y > 0.15 && return 0.0
            peak = typ == 0 ? 0.25 + 0.17x - 0.12y : 0.21 - 0.19x + 0.24y
            amplitude = 0.025 * exp(-((x + 0.04)^2 / 0.12 + (y - 0.09)^2 / 0.2))
            amplitude *= typ == 0 ? (1 + 0.28x + 0.1y) : (1 + 0.3sin(4x + 2y))
            return amplitude * exp(-((z - peak) / 0.27)^2)
        end
        push!(cubes, write_cube(joinpath(dir, "type$typ.cube"), dims, origin, axes, density))
    end
    f0 = write_frame(joinpath(dir, "frame0.tsv"), [0.0, 0.0, 0.0], [0.97, 0.24, 0.02], [-0.2, 0.94, 0.01])
    f1 = write_frame(joinpath(dir, "frame1.tsv"), [0.015, -0.01, 0.0], [0.96, -0.20, 0.01], [0.25, 0.93, -0.02])
    return cubes[1], cubes[2], f0, f1
end

# Calls the original Python implementation, not a rewritten reference. All
# reference defaults used by surface_grid/legacy calibration are checked against
# the frozen explicit config before it runs. This script never reads predictions.
const PYTHON_REFERENCE = raw"""
import importlib.util, os, sys, tomllib
import numpy as np
root, c0, c1, f0, f1, config, outdir = sys.argv[1:]
spec = importlib.util.spec_from_file_location("reference_cc", os.path.join(root, "test/lib/cc_mold_builder.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.os = os  # reference imports os only from its __main__ block
with open(config, "rb") as f:
    cfg = tomllib.load(f)["model"]
assert cfg["mold_z_min_nm"] == -0.5 and cfg["mold_z_max_nm"] == 2.6
assert cfg["mold_isovalue_min_log10"] == -5.5
assert cfg["mold_isovalue_max_fraction"] == 0.8 and cfg["mold_isovalue_count"] == 80
assert cfg["mold_calibration"] == "first_below_target"
assert cfg["mold_z_step_nm"] == 0.005  # build_templates has this reference default
half, step, target = (cfg[k] for k in ("mold_half_nm", "mold_step_nm", "mold_target_height_nm"))
xy = np.arange(-half, half + 0.001, step)
for typ, cube, framepath, label in ((0,c0,f0,"GlcN"), (1,c1,f1,"GlcNAc")):
    d = {}
    for line in open(framepath):
        p = line.strip().split("\t")
        if len(p) == 2 and p[0] in ("origin_nm", "t_axis", "u_axis"):
            d[p[0]] = np.array([float(x) for x in p[1].split(",")])
    frame = tuple(d[k] for k in ("origin_nm", "t_axis", "u_axis"))
    zs, vc = m.surface_grid(cube, frame, xy, xy, cfg["mold_z_step_nm"])
    iso, mh, nv, heights = m.iso_for_mean_height_legacy(vc, zs, target)
    zf = np.where(np.isfinite(heights), heights, np.nan)
    normalized = np.where(np.isfinite(zf), (zf - np.nanmean(zf)) / np.nanstd(zf), 0.0)
    prefix = os.path.join(outdir, "python_type" + str(typ))
    for name, data in (("zs",zs), ("vc",vc), ("selection",[iso,mh,nv]), ("heights",heights), ("normalized",normalized)):
        np.savetxt(prefix + "_" + name + ".tsv", data, delimiter="\t", fmt="%.17g")
    m.build_templates(cube, frame, label, os.path.join(outdir,"python_templates.tsv"), target, half, half, step, legacy=True)
"""

function same_numerics(actual, expected; rtol=2e-12, atol=2e-14)
    @test size(actual) == size(expected)
    @test isnan.(actual) == isnan.(expected)
    good = isfinite.(expected)
    @test all(isfinite, actual[good])
    @test all(isapprox.(actual[good], expected[good]; rtol=rtol, atol=atol))
end

function compare_python(cube0, cube1, frame0, frame1, outdir)
    native = build_molds(cube0=cube0, cube1=cube1, frame0=frame0, frame1=frame1,
                         config=CONFIG, out=joinpath(outdir, "native_templates.tsv"))
    # The command name is intentional: no REPL interpreter or Python production path.
    run(`python3 -c $PYTHON_REFERENCE $ROOT $cube0 $cube1 $frame0 $frame1 $CONFIG $outdir`)
    for (typ, cube, frame, result) in ((0,cube0,frame0,native.type0), (1,cube1,frame1,native.type1))
        zs, vc = surface_grid(read_cube(cube), read_frame(frame), SETTINGS)
        prefix = joinpath(outdir, "python_type$(typ)_")
        ref_zs = vec(readdlm(prefix * "zs.tsv", '\t', Float64))
        ref_vc = readdlm(prefix * "vc.tsv", '\t', Float64)
        ref_selection = vec(readdlm(prefix * "selection.tsv", '\t', Float64))
        ref_heights = vec(readdlm(prefix * "heights.tsv", '\t', Float64))
        ref_normalized = vec(readdlm(prefix * "normalized.tsv", '\t', Float64))
        same_numerics(zs, ref_zs)
        same_numerics(vc, ref_vc)
        same_numerics([result.iso, result.mean_height, result.nvalid], ref_selection)
        @test isequal(result.heights, ref_heights) # exact grid-sample decisions
        same_numerics(result.normalized, ref_normalized)
        finite = isfinite.(ref_vc)
        println("type$typ parity: max |LDOS Δ|=", maximum(abs.(vc[finite] - ref_vc[finite])),
                "; max |zscore Δ|=", maximum(abs.(result.normalized - ref_normalized)),
                "; iso=", result.iso, "; nvalid=", result.nvalid)
    end
    a, b = readlines(joinpath(outdir,"native_templates.tsv")), readlines(joinpath(outdir,"python_templates.tsv"))
    @test length(a) == length(b) == 9
    @test a[1] == b[1]
    for i in 2:9
        ac, bc = split(a[i], '\t'), split(b[i], '\t')
        @test ac[1:4] == bc[1:4]
        # Both writers use seven significant digits. Permit one final digit of
        # rounding at an exact decimal tie, not a height-map/indexing change.
        same_numerics(parse.(Float64, ac[5:end]), parse.(Float64, bc[5:end]); rtol=2e-6, atol=1e-7)
    end
    println("template TSV byte equality: ", a == b)
    return native
end

function synthetic_tests()
    @testset "explicit mold config" begin
        @test VERSION >= v"1.13.0" && VERSION < v"1.14.0"
        @test SETTINGS.target_height_nm == 0.50
        @test SETTINGS.half_nm == 0.32
        @test SETTINGS.isovalue_count == 80
        @test_throws ArgumentError CCMoldNative._validate(changed_settings(step_nm=0.03))
        @test_throws ArgumentError CCMoldNative._validate(changed_settings(z_step_nm=0.0))
        @test_throws ArgumentError CCMoldNative._validate(changed_settings(target_height_nm=NaN))
        @test_throws ArgumentError CCMoldNative._validate(changed_settings(isovalue_count=1))
        @test_throws ArgumentError CCMoldNative._validate(changed_settings(calibration="nearest"))
        mktempdir() do dir
            cfg = TOML.parsefile(CONFIG)
            delete!(cfg["model"], "mold_z_step_nm")
            path = joinpath(dir, "missing.toml")
            open(io -> TOML.print(io, cfg), path, "w")
            @test_throws ArgumentError load_settings(path)
        end
    end
    mktempdir() do dir
        @testset "cube indexing, skew axes, boundary weights and invalid input" begin
            dims = (4, 3, 5)
            origin = [-0.1, 0.2, -0.3]
            axes = [0.13 0.01 0.02; 0.02 0.17 0.01; 0.0 0.03 0.19]
            path = write_cube(joinpath(dir,"index.cube"), dims, origin, axes, (i,j,k) -> 1+i+10j+100k; signed_axes=true)
            cube = read_cube(path)
            @test cube.dims == dims
            @test cube.origin_nm ≈ origin
            @test cube.axes_nm ≈ axes
            @test cube.values[1:5] == [1,2,3,4,11]
            c = [1.25 0.5 2.75; 0.0 0.0 0.0; -0.5 0.0 0.0; -1.1 0.0 0.0; 3.5 2.0 4.0]
            points = c * transpose(cube.axes_nm) .+ transpose(cube.origin_nm)
            sampled = sample_volume(cube, points)
            @test sampled[1] ≈ 1 + 1.25 + 10*0.5 + 100*2.75
            @test sampled[2] ≈ 1.0
            @test sampled[3] ≈ 0.5 # not renormalized to 1
            @test isnan(sampled[4])
            @test sampled[5] ≈ 0.5 * (1+3+20+400)
            @test_throws ArgumentError sample_volume(cube, fill(NaN, 1, 3))
            lines = readlines(path)
            for (name, bad) in (("short", lines[1:end-1]), ("extra", vcat(lines,["1"])),
                                ("nonnumeric", vcat(lines[1:end-1],["broken"])),
                                ("nonfinite", vcat(lines[1:end-1],["Inf"])))
                p = joinpath(dir,"$name.cube"); write(p, join(bad,'\n') * "\n")
                @test_throws ArgumentError read_cube(p)
            end
            singular = copy(lines); singular[5] = singular[4]
            p = joinpath(dir,"singular.cube"); write(p, join(singular,'\n') * "\n")
            @test_throws ArgumentError read_cube(p)
            p = joinpath(dir,"empty.cube"); write(p, "")
            @test_throws ArgumentError read_cube(p)
        end
        @testset "frame validation and non-normalized t/u" begin
            f = write_frame(joinpath(dir,"frame.tsv"), [0,0,0], [2,0,0], [0,3,0])
            frame = read_frame(f)
            @test frame.t_axis == [2,0,0]
            @test frame.u_axis == [0,3,0]
            @test CCMoldNative._normal(frame) == [0,0,1]
            @test_throws ArgumentError CCMoldNative._normal(Frame([0.,0,0], [1.,0,0], [2.,0,0]))
            write(f, read(f,String) * "t_axis\t1,0,0\n")
            @test_throws ArgumentError read_frame(f)
            write(f, "origin_nm\t0,0,0\nu_axis\t0,1,0\n")
            @test_throws ArgumentError read_frame(f)
        end
        @testset "highest occupied height, first below target, no support minimum" begin
            s = changed_settings(isovalue_min_log10=0.0, isovalue_max_fraction=1.0, isovalue_count=3)
            vc = [2.0 0.0 2.0 0.0; zeros(4,4)]
            zs = [0.0,0.2,0.4,0.6]
            chosen = iso_for_mean_height_legacy(vc, zs, s)
            @test chosen.iso == 1.0
            @test chosen.nvalid == 1 # 20% support is allowed by the legacy path
            @test chosen.heights[1] == 0.4 # choose the upper island, not the first crossing
            @test all(isnan, chosen.heights[2:end])
            @test_throws ArgumentError iso_for_mean_height_legacy(vc,zs,changed_settings(
                target_height_nm=0.4,isovalue_min_log10=0.0,isovalue_max_fraction=1.0,isovalue_count=3))
            @test_throws ArgumentError iso_for_mean_height_legacy(zeros(2,4),zs,s)
            @test_throws ArgumentError iso_for_mean_height_legacy(fill(NaN,2,4),zs,s)
            @test_throws ArgumentError iso_for_mean_height_legacy(fill(Inf,2,4),zs,s)
            @test_throws ArgumentError iso_for_mean_height_legacy(vc,reverse(zs),s)
            normalized = CCMoldNative._zscore([1.,2.,3.,NaN])
            @test normalized ≈ [-sqrt(1.5),0,sqrt(1.5),0]
            @test_throws ArgumentError CCMoldNative._zscore([0.1,0.1,NaN])
            @test_throws ArgumentError CCMoldNative._zscore([NaN,NaN])
        end
        cube0,cube1,frame0,frame1 = synthetic_pair(dir)
        @testset "reference Python numerical intermediates and TSV" begin
            result = compare_python(cube0,cube1,frame0,frame1,dir)
            for r in (result.type0,result.type1)
                @test size(r.variants) == (4,289)
                @test r.variants[1,:] == r.normalized
                @test r.variants[2,1:17] == reverse(r.normalized[1:17])
                @test r.variants[3,1:17] == r.normalized[(end-16):end]
                @test r.variants[4,:] == reverse(r.normalized)
                @test r.nvalid < 289 # synthetic absent columns
            end
        end
        @testset "output safety and standalone CLI" begin
            existing = joinpath(dir,"existing.tsv"); write(existing,"keep exactly\n")
            @test_throws ArgumentError build_molds(cube0=cube0,cube1=cube1,frame0=frame0,frame1=frame1,config=CONFIG,out=existing)
            @test read(existing,String) == "keep exactly\n"
            link = joinpath(dir,"link.tsv"); symlink(existing,link)
            @test_throws ArgumentError build_molds(cube0=cube0,cube1=cube1,frame0=frame0,frame1=frame1,config=CONFIG,out=link)
            dangling = joinpath(dir,"dangling.tsv"); symlink(joinpath(dir,"absent"),dangling)
            @test_throws ArgumentError build_molds(cube0=cube0,cube1=cube1,frame0=frame0,frame1=frame1,config=CONFIG,out=dangling)
            failedout = joinpath(dir,"failed.tsv")
            badframe = joinpath(dir,"badframe.tsv"); write(badframe,"")
            @test_throws ArgumentError build_molds(cube0=cube0,cube1=cube1,frame0=frame0,frame1=badframe,config=CONFIG,out=failedout)
            @test !ispath(failedout)
            cli = joinpath(@__DIR__,"build_cc_molds_native.jl")
            @test occursin("All six options", read(`$(Base.julia_cmd()) --project=$ROOT $cli --help`,String))
            out = joinpath(dir,"cli.tsv")
            args = ["--cube0",cube0,"--cube1",cube1,"--frame0",frame0,"--frame1",frame1,"--config",CONFIG,"--out",out]
            run(pipeline(`$(Base.julia_cmd()) --project=$ROOT $cli $args`; stdout=devnull))
            @test read(out) == read(joinpath(dir,"native_templates.tsv"))
            @test !success(pipeline(`$(Base.julia_cmd()) --project=$ROOT $cli $args`; stdout=devnull,stderr=devnull))
            @test !success(pipeline(`$(Base.julia_cmd()) --project=$ROOT $cli --cube0`; stdout=devnull,stderr=devnull))
            @test !success(pipeline(`$(Base.julia_cmd()) --project=$ROOT $cli --unknown=x`; stdout=devnull,stderr=devnull))
        end
    end
end

if isempty(ARGS)
    synthetic_tests()
elseif length(ARGS) == 6 && ARGS[1] == "--real-cubes"
    @test VERSION >= v"1.13.0" && VERSION < v"1.14.0"
    outdir = ARGS[6]
    (ispath(outdir) || islink(outdir)) && error("refusing existing test output directory: $outdir")
    mkpath(outdir)
    @testset "one DFT pair: native/reference numerical parity (no fits)" begin
        compare_python(ARGS[2],ARGS[3],ARGS[4],ARGS[5],outdir)
    end
else
    error("usage: test/test_cc_mold_native.jl [--real-cubes CUBE0 CUBE1 FRAME0 FRAME1 NEW_OUT_DIR]")
end

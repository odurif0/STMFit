#!/usr/bin/env julia
# Focused native checks: julia --project=. test/test_patch_preprocessing.jl
# Synthetic SXM only. No fitting, benchmark inputs, or real scans.
using Test
using TOML

module ForwardPatches
include(joinpath(@__DIR__, "extract_lobe_patches.jl"))
end
module BackwardPatches
include(joinpath(@__DIR__, "extract_lobe_patches_bwd.jl"))
end

include(joinpath(@__DIR__, "lib", "patch_preprocessing.jl"))
using .PatchPreprocessing

const ROOT = dirname(@__DIR__)
const COUNT_CONFIG = joinpath(ROOT, "config", "chitosan.toml")
const LEGACY = (stride=1, flatten="plane+rows", smooth_radius_px=1)
settings_tuple(p) = (stride=p.stride, flatten=p.flatten, smooth_radius_px=p.smooth_radius_px)

function write_config(path, fields)
    open(path, "w") do io
        TOML.print(io, Dict("preprocessing" => fields))
    end
    return path
end

function synthetic_sxm(path)
    n = 33
    xs = range(0, 2.4; length=n)
    fwd = [0.4 + 0.17x + 0.08y + 2exp(-((x-1.2)^2/0.07 + (y-1.2)^2/0.04)) +
           0.25cos(8x+3y) + 0.14sin(9y) for y in xs, x in xs]
    bwd = [0.4 + 0.17x + 0.08y + 1.8exp(-((x-1.2)^2/0.07 + (y-1.2)^2/0.04)) +
           0.16cos(6x+2y) - 0.11sin(7y) for y in xs, x in xs]
    open(path, "w") do io
        println(io, ":SCAN_PIXELS:\n$n $n\n:SCAN_RANGE:\n2.4e-9 2.4e-9")
        println(io, ":DATA_INFO:\nChannel Name Unit Direction Calibration Offset\n0 Z nm both 1 0")
        println(io, ":SCANIT_END:")
        for image in (fwd, reverse(bwd; dims=2)), row in axes(image,1), col in axes(image,2)
            bits = reinterpret(UInt32, Float32(image[row,col]))
            write(io, hton(bits))
        end
    end
    return path
end

function extract_quietly(mod, args)
    redirect_stdout(devnull) do
        mod.main(args)
    end
end

@testset "patch preprocessing validation" begin
    @test VERSION >= v"1.13.0" && VERSION < v"1.14.0"
    @test settings_tuple(load_patch_preprocessing()) == LEGACY
    @test settings_tuple(load_patch_preprocessing(COUNT_CONFIG)) == LEGACY
    @test_throws ArgumentError load_patch_preprocessing("")
    mktempdir() do dir
        path = joinpath(dir,"config.toml")
        @test_throws ArgumentError load_patch_preprocessing(path)
        for flatten in ("none", "plane", "rows", "plane+rows", "PLANE+ROWS")
            fields = Dict("stride"=>2, "flatten"=>flatten, "smooth_radius_px"=>0)
            write_config(path,fields)
            @test settings_tuple(load_patch_preprocessing(path)) == (stride=2,flatten=flatten,smooth_radius_px=0)
        end
        fields = Dict{String,Any}("stride"=>1,"flatten"=>"plane+rows","smooth_radius_px"=>1)
        for key in keys(fields)
            incomplete = copy(fields); delete!(incomplete,key)
            write_config(path,incomplete)
            @test_throws ArgumentError load_patch_preprocessing(path)
        end
        for (key, invalid) in (("stride",0), ("stride",-1), ("stride",1.0), ("stride",true),
                               ("smooth_radius_px",-1), ("smooth_radius_px",0.5), ("smooth_radius_px",false),
                               ("flatten","plnae"), ("flatten",""), ("flatten",1))
            bad = copy(fields); bad[key] = invalid
            write_config(path,bad)
            @test_throws ArgumentError load_patch_preprocessing(path)
        end
        write(path,"[model]\nname = \"missing preprocessing\"\n")
        @test_throws ArgumentError load_patch_preprocessing(path)
        write(path,"preprocessing = 3\n")
        @test_throws ArgumentError load_patch_preprocessing(path)
    end
end

mktempdir() do dir
    feature_path = joinpath(dir,"features.tsv")
    write(feature_path, "file\tlobe\tt_nm\tu_nm\tamplitude\taxis_x\taxis_y\tx_nm\ty_nm\tsigma_parallel_nm\tsigma_perp_nm\tbaseline\ttilt_x\ttilt_y\tskew_ratio\n" *
                        "synthetic.sxm\t1\t0\t0\t2\t1\t0\t1.2\t1.2\t0.2\t0.15\t0.4\t0.17\t0.08\t1\n")
    synthetic_sxm(joinpath(dir,"synthetic.sxm"))
    config = joinpath(dir,"custom.toml")
    defaults = Dict{String,Any}("stride"=>1,"flatten"=>"plane+rows","smooth_radius_px"=>1)
    write_config(config,Dict("stride"=>2,"flatten"=>"none","smooth_radius_px"=>0))
    for (name, mod, cli) in (("forward",ForwardPatches,"extract_lobe_patches.jl"),
                             ("backward",BackwardPatches,"extract_lobe_patches_bwd.jl"))
        @testset "$name parsing and synthetic preprocessing wiring" begin
            output = joinpath(dir,"$name.tsv")
            args = ["--features",feature_path,"--data-dir",dir,"--out",output]
            @test settings_tuple(mod._parse_cli(args).preprocessing) == LEGACY
            for suffix in (["--config",config],["--config=$config"])
                @test settings_tuple(mod._parse_cli(vcat(args,suffix)).preprocessing) ==
                      (stride=2,flatten="none",smooth_radius_px=0)
            end
            @test_throws ErrorException mod._parse_cli(vcat(args,["--config"]))
            @test_throws ErrorException mod._parse_cli(vcat(args,["--config","--out"]))
            @test_throws ErrorException mod._parse_cli(vcat(args,["--config",config,"--config",config]))
            @test_throws ArgumentError mod._parse_cli(vcat(args,["--config="]))
            @test_throws ArgumentError mod._parse_cli(vcat(args,["--config",joinpath(dir,"absent.toml")]))
            help = read(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,cli)) --help`,String)
            @test occursin("--config PATH",help)
            @test occursin("smooth_radius_px=1",help)

            extract_quietly(mod,args)
            legacy_bytes = read(output)
            @test length(readlines(output)) == 2
            extract_quietly(mod,vcat(args,["--config",COUNT_CONFIG]))
            @test read(output) == legacy_bytes # identical behavior on explicit old settings
            for (key,value) in (("stride",2),("flatten","none"),("smooth_radius_px",0))
                config_fields = copy(defaults); config_fields[key] = value
                cfg = write_config(joinpath(dir,"$(name)_$(key).toml"),config_fields)
                extract_quietly(mod,vcat(args,["--config",cfg]))
                @test length(readlines(output)) == 2
                @test read(output) != legacy_bytes # all three fields reach actual preprocessing
            end

            missing = joinpath(dir,"missing_features.tsv")
            write(missing,"file\tlobe\nnot_here.sxm\t1\n")
            missing_args = ["--features",missing,"--data-dir",dir,"--out",output,"--config",COUNT_CONFIG]
            @test_logs (:warn,r"SXM not found, skipping") extract_quietly(mod,missing_args)
            @test length(readlines(output)) == 1 # missing-file behavior remains a skipped row

            write_config(config,Dict("stride"=>0,"flatten"=>"none","smooth_radius_px"=>0))
            sentinel = "do not touch\n"; write(output,sentinel)
            @test_throws ArgumentError extract_quietly(mod,vcat(args,["--config",config]))
            @test read(output,String) == sentinel # reject invalid config before output creation
            write_config(config,Dict("stride"=>2,"flatten"=>"none","smooth_radius_px"=>0))
        end
    end
end

using Test, LinearAlgebra, Statistics, TOML, SHA
include(joinpath(@__DIR__, "lib", "patch_frames.jl"))
using .PatchFrames
const FRAME_CONFIG = joinpath(dirname(@__DIR__), "config", "local_patch_orientation.toml")

function geometry_rows(ts; theta=0.0, norm=1.0, offset=(3.0, -2.0), curve=(0.15, 0.2, -0.07))
    ax, ay = norm*cos(theta), norm*sin(theta)
    [Dict("file"=>"synthetic.sxm", "lobe"=>string(i), "N"=>string(length(ts)),
        "axis_x"=>string(ax), "axis_y"=>string(ay),
        "x_nm"=>string(offset[1] + t*cos(theta) - (curve[1]+curve[2]*t+curve[3]*t^2)*sin(theta)),
        "y_nm"=>string(offset[2] + t*sin(theta) + (curve[1]+curve[2]*t+curve[3]*t^2)*cos(theta))) for (i,t) in enumerate(ts)]
end

@testset "Local tangent from frozen geometry, no label input" begin
    settings = load_frame_settings(FRAME_CONFIG)
    @test settings == (degree=2,)
    for theta in (0.0, 0.3, 1.5, 3.1), n in (1,2,3,6,12), norm in (1.0, 0.9999996)
        ts = collect(range(-1, 1; length=max(n,2)))[1:n]
        rows = geometry_rows(ts; theta, norm)
        before = deepcopy(rows)
        frames = local_frame_rows(rows, settings)
        @test rows == before
        @test length(frames) == n
        for (i,r) in enumerate(frames)
            slope = n == 1 ? 0.0 : n == 2 ? .2-.07sum(ts) : .2-.14ts[i]
            expected = norm .* [cos(theta+atan(slope)), sin(theta+atan(slope))]
            axis = parse.(Float64, [r["frame_axis_x"],r["frame_axis_y"]])
            @test axis ≈ expected atol=3e-14
            @test hypot(axis...) ≈ norm atol=2e-15
            @test parse(Float64,r["tangent_slope"]) ≈ slope atol=3e-14
        end
        mktempdir() do dir
            path = write_frames(joinpath(dir,"frames.tsv"),frames)
            restored = read_frames(path, rows)
            @test length(restored) == n
            for r in rows
                @test patch_axis(restored,r,0,0) == Tuple(parse.(Float64,
                    [frames[parse(Int,r["lobe"])][k] for k in ("frame_axis_x","frame_axis_y")]))
                @test patch_axis(nothing,r,0.7,-0.3) == (0.7,-0.3)
            end
        end
    end
    # An exactly straight global control is not implicitly normalized.
    rows = geometry_rows(-2:2; norm=.9999996, curve=(0.0,0.0,0.0))
    frames = local_frame_rows(rows,settings)
    @test all(parse(Float64,r["frame_axis_x"]) == .9999996 && parse(Float64,r["frame_axis_y"]) == 0 for r in frames)
    # Chain reversal rotates both axes by pi, not a one-axis reflection.
    rows = geometry_rows(collect(-2:2))
    original = local_frame_rows(rows,settings)
    reversed = deepcopy(reverse(rows))
    for (i,r) in enumerate(reversed)
        r["lobe"] = string(i)
        for k in ("axis_x","axis_y"); r[k] = string(-parse(Float64,r[k])); end
    end
    reverse_frames = local_frame_rows(reversed,settings)
    for (a,b) in zip(original,reverse(reverse_frames)), k in ("frame_axis_x","frame_axis_y")
        @test parse(Float64,a[k]) ≈ -parse(Float64,b[k]) atol=2e-15
    end
    # Label-shaped additions cannot affect the centers-only calculation.
    extra = deepcopy(rows)
    for r in extra
        r["expected_N"] = "999"; r["class"] = "ignored"; r["t_nm"] = "NaN"
    end
    @test local_frame_rows(extra,settings) == original
    for mutate! in (r->push!(r,deepcopy(r[1])), r->(r[2]["lobe"]="99"),
        r->(r[2]["x_nm"]=r[1]["x_nm"]), r->(r[2]["x_nm"]="NaN"),
        r->(r[2]["axis_x"]="0.9"), r->(r[1]["axis_x"]="0"))
        bad = deepcopy(rows); mutate!(bad)
        @test_throws ErrorException local_frame_rows(bad,settings)
    end
    mktempdir() do dir
        path = joinpath(dir,"frames.tsv")
        for mutate! in (r->pop!(r), r->push!(r,deepcopy(r[1])),
            r->(r[1]["x_nm"]="987"), r->(r[1]["frame_axis_y"]="NaN"),
            r->(r[1]["frame_axis_x"]="-1"), r->(r[1]["centerline_rms_nm"]="-1"))
            bad = deepcopy(original); mutate!(bad); write_frames(path,bad)
            @test_throws ErrorException read_frames(path,rows)
        end
        write_frames(path,original)
        write(path, replace(read(path,String), "tangent_slope"=>"label"))
        @test_throws ErrorException read_frames(path,rows)
        cfg = TOML.parsefile(FRAME_CONFIG)
        for (section,key,value) in (("model","centerline_degree",true), ("model","centerline_degree",3),
            ("selection","centerline_fit","label_search"), ("preprocessing","patch_orientation","global"),
            ("model","expected_N",6))
            bad = deepcopy(cfg); bad[section][key] = value
            p = joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),p,"w")
            @test_throws ErrorException load_frame_settings(p)
        end
    end
end

module Comparison
include(joinpath(@__DIR__, "run_local_orientation_comparison.jl"))
end

@testset "Two-arm comparison wiring and failures" begin
    @test_throws ErrorException Comparison.orientation_comparison(["--truth","labels.tsv"])
    mktempdir() do dir
        raw = joinpath(dir,"raw"); mkdir(raw); write(joinpath(raw,"synthetic.sxm"),"fixture, no SXM compute")
        rows = geometry_rows(-2:2)
        features = joinpath(dir,"features.tsv")
        Comparison.write_table(features,sort(collect(keys(rows[1]))),rows)
        templates = joinpath(dir,"templates.tsv"); write(templates,"fixture")
        out = joinpath(dir,"out")
        args = ["--data-dir",raw,"--count-config",joinpath(Comparison.ROOT,"config","chitosan.toml"),
            "--config",joinpath(Comparison.ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",FRAME_CONFIG,"--features",features,"--split-features",features,
            "--templates",templates,"--outdir",out]
        Comparison.orientation_comparison(vcat(args,["--dry-run"]))
        @test !ispath(out)
        calls = String[]
        function runner(output,name,script,cli)
            @test script == "run_reconstructed_chitosan.jl"
            opts = Comparison.parse_options(cli)
            push!(calls,name)
            @test haskey(opts,"--patch-frames") == (name == "local_frame")
            if name == "local_frame"
                @test length(Comparison.PatchFrames.read_frames(opts["--patch-frames"],rows)) == 5
                # Production rejects patch caches and unbound geometry before writing.
                for k in ("--patches-fwd","--patches-bwd","--descriptor-patches")
                    @test_throws ErrorException Comparison.execute_pipeline(merge(opts,Dict(k=>features)))
                end
                missing_base = copy(opts); delete!(missing_base,"--features")
                @test_throws ErrorException Comparison.execute_pipeline(missing_base)
            end
            mkdir(opts["--outdir"])
            cp(features,joinpath(opts["--outdir"],"predictions.tsv"))
        end
        Comparison.orientation_comparison(args;runner)
        @test calls == ["reference","local_frame"]
        for mode in calls
            @test read(joinpath(out,mode,"features.tsv")) == read(features)
            @test read(joinpath(out,mode,"features_split.tsv")) == read(features)
        end
        @test !isfile(joinpath(out,"failures.tsv"))
        @test_throws ErrorException Comparison.orientation_comparison(args;runner)
        failed_args = copy(args); failed_args[end] = joinpath(dir,"failed")
        @test_throws ErrorException Comparison.orientation_comparison(failed_args;runner=(a...)->error("fixture failure"))
        _,failed = Comparison.read_table(joinpath(dir,"failed","failures.tsv"))
        @test only(failed)["stage"] == "reference"
    end
end

#!/usr/bin/env julia
# Focused native checks: julia --project=. test/test_patch_preprocessing.jl
# Synthetic SXM only. No fitting, benchmark inputs, or real scans.
using Test
using TOML
using Statistics
using STMSXMIO: _box_smooth

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
const ASSIGNMENT_CONFIG = joinpath(ROOT, "config", "unit_assignment_reconstructed.toml")
const MATCHED_CONFIG = joinpath(ROOT, "config", "unit_assignment_matched_residual.toml")
const LEGACY = (stride=1, flatten="plane+rows", smooth_radius_px=1)
settings_tuple(p) = (stride=p.stride, flatten=p.flatten, smooth_radius_px=p.smooth_radius_px)

function write_config(path, fields)
    open(path, "w") do io
        TOML.print(io, Dict("preprocessing" => fields))
    end
    return path
end

@testset "Explicit residual policy and matched filtering" begin
    @test load_patch_residual_filter() == "smooth_data_only"
    @test load_patch_residual_filter(ASSIGNMENT_CONFIG) == "smooth_data_only"
    @test load_patch_residual_filter(MATCHED_CONFIG) == "smooth_residual"
    candidate, reference = TOML.parsefile(MATCHED_CONFIG), TOML.parsefile(ASSIGNMENT_CONFIG)
    @test candidate["model"]["name"] == "cc_soft_matched_residual_v1"
    candidate["model"]["name"] = reference["model"]["name"]
    candidate["preprocessing"]["patch_residual_filter"] = reference["preprocessing"]["patch_residual_filter"]
    @test candidate == reference # Exactly one scientific setting differs.
    @test_throws ArgumentError load_patch_residual_filter("")
    mktempdir() do dir
        cfg = joinpath(dir,"assignment.toml")
        for value in ("unknown", "", 1, true)
            write_config(cfg, Dict("patch_residual_filter"=>value))
            @test_throws ArgumentError load_patch_residual_filter(cfg)
        end
        write_config(cfg, Dict{String,String}())
        @test_throws ArgumentError load_patch_residual_filter(cfg)
        write(cfg, "[model]\nname = \"synthetic\"\n")
        @test_throws ArgumentError load_patch_residual_filter(cfg)
    end
    xs, ys = collect(range(-1.2, 1.2; length=31)), collect(range(-0.9, 0.9; length=25))
    model = [0.4 + 0.07x - 0.03y + 2exp(-0.5((x/.20)^2+(y/.14)^2)) for y in ys, x in xs]
    shoulder = [.08exp(-0.5(((x-.09)/.11)^2+((y-.22)/.09)^2)) for y in ys, x in xs]
    for radius in (0,1,2)
        settings = PreprocessingSettings(1,"none",radius)
        smooth = _box_smooth(model,radius)
        null = patch_residual(model,smooth,model,settings,"smooth_residual")
        legacy = patch_residual(model,smooth,model,settings,"smooth_data_only")
        @test all(iszero,null)
        @test legacy == smooth-model
        if radius > 0
            @test maximum(abs,legacy) > 0.01
        else
            @test all(iszero,legacy)
        end
        z = model+shoulder
        snapshot = (copy(z),copy(model),copy(smooth))
        matched = patch_residual(z,_box_smooth(z,radius),model,settings,"smooth_residual")
        @test matched ≈ _box_smooth(shoulder,radius) atol=5e-16 rtol=5e-14
        @test matched ≈ _box_smooth(z,radius)-_box_smooth(model,radius) atol=1e-15 rtol=1e-13
        # Independent finite-window calculation, including edges/corners.
        expected = [mean((z-model)[max(1,y-radius):min(end,y+radius),max(1,x-radius):min(end,x+radius)])
                    for y in axes(z,1), x in axes(z,2)]
        # A copied window and a strided view can reduce in different orders.
        @test maximum(abs,matched-expected) <= 8eps(maximum(abs,expected))
        @test (z,model,smooth) == snapshot
        forward, backward = model+shoulder, model-0.4shoulder
        rf = patch_residual(forward,_box_smooth(forward,radius),model,settings,"smooth_residual")
        rb = patch_residual(backward,_box_smooth(backward,radius),model,settings,"smooth_residual")
        @test rf-rb ≈ _box_smooth(forward-backward,radius) atol=5e-16
        if radius == 0
            @test matched == patch_residual(z,z,model,settings,"smooth_data_only")
        end
        bad = copy(z); bad[12,16] = NaN
        rb = patch_residual(bad,_box_smooth(bad,radius),model,settings,"smooth_residual")
        old = patch_residual(bad,_box_smooth(bad,radius),model,settings,"smooth_data_only")
        @test isfinite.(rb) == isfinite.(old) # No missing-sample repair or support relaxation.
        @test count(!isfinite,rb) == (2radius+1)^2
    end
    @test_throws DimensionMismatch patch_residual(model,model,zeros(2,2),PreprocessingSettings(1,"none",1),"smooth_residual")
    @test_throws ArgumentError patch_residual(model,model,model,PreprocessingSettings(1,"none",-1),"smooth_residual")
    @test_throws ArgumentError patch_residual(model,model,model,PreprocessingSettings(1,"none",1),"unknown")
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
            @test mod._parse_cli(args).residual_filter == "smooth_data_only"
            for suffix in (["--assignment-config",MATCHED_CONFIG],["--assignment-config=$MATCHED_CONFIG"])
                @test mod._parse_cli(vcat(args,suffix)).residual_filter == "smooth_residual"
            end
            @test_throws ErrorException mod._parse_cli(vcat(args,["--assignment-config"]))
            @test_throws ErrorException mod._parse_cli(vcat(args,["--assignment-config","--out"]))
            @test_throws ErrorException mod._parse_cli(vcat(args,["--assignment-config",MATCHED_CONFIG,"--assignment-config",MATCHED_CONFIG]))
            @test_throws ArgumentError mod._parse_cli(vcat(args,["--assignment-config="]))
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
            extract_quietly(mod,vcat(args,["--config",COUNT_CONFIG,"--assignment-config",ASSIGNMENT_CONFIG]))
            @test read(output) == legacy_bytes
            extract_quietly(mod,vcat(args,["--config",COUNT_CONFIG,"--assignment-config",MATCHED_CONFIG]))
            @test read(output) != legacy_bytes
            header, actual = mod.ScriptUtils._read_tsv(output)
            legacy_lines = split(chomp(String(copy(legacy_bytes))), '\n')
            old = Dict(zip(split(legacy_lines[1],'\t'), split(legacy_lines[2],'\t')))
            @test length(actual) == 1
            for column in header
                if startswith(column,"raw_p") || startswith(column,"bwd_raw_p") || startswith(column,"diff_raw_p")
                    @test actual[1][column] == old[column]
                end
            end
            # End-to-end extraction agrees with an independent model and S(z-M),
            # before the existing interpolation, normalization and TSV rounding.
            sxm = joinpath(dir,"synthetic.sxm")
            img = mod.read_sxm(sxm)
            direction = name == "forward" ? "fwd" : "bwd"
            cfg = mod.PatternConfig(filepath=sxm,channel="Z",direction=direction,
                stride=1,flatten="plane+rows",smooth_radius_px=1,output_dir=dir,no_plot=true)
            xs, ys, _, z, _, _, _ = mod.preprocess_channel(img,mod.get_channel(img,"Z"; direction),cfg)
            model = [.4+.17x+.08y+2exp(-.5(((x-1.2)/.2)^2+((y-1.2)/.15)^2)) for y in ys, x in xs]
            residual = _box_smooth(z-model,1)
            coords = collect(-.32:.08:.32)
            values = [mod._interp(xs,ys,residual,1.2+t,1.2+u) for u in coords for t in coords]
            expected = (values .- median(values)) ./ std(values)
            prefix = name == "forward" ? "res_p" : "bwd_res_p"
            measured = [parse(Float64,actual[1][prefix*lpad(string(i),3,'0')]) for i in eachindex(expected)]
            @test maximum(abs,measured-expected) < 5e-6 # seven significant TSV digits
            # Opt-in observation-mask zero transform is byte-identical when raw
            # pixels are finite. Only backward moves, before model subtraction.
            matched_bytes=read(output)
            shifts=joinpath(dir,"$(name)_shifts.tsv")
            write(shifts,"file\tbwd_sample_dx_px\nsynthetic.sxm\t0\n")
            registered_args=vcat(args,["--config",COUNT_CONFIG,"--assignment-config",MATCHED_CONFIG,
                "--acquisition-shifts",shifts])
            extract_quietly(mod,registered_args)
            @test read(output)==matched_bytes
            write(shifts,"file\tbwd_sample_dx_px\nsynthetic.sxm\t2\n")
            extract_quietly(mod,registered_args)
            if name=="forward"
                @test read(output)==matched_bytes
            else
                @test read(output)!=matched_bytes
                moved=fill(NaN,size(z)); moved[:,1:end-2] .= z[:,3:end]
                independent=_box_smooth(moved-model,1)
                vals=[mod._interp(xs,ys,independent,1.2+t,1.2+u) for u in coords for t in coords]
                want=mod._normalize_patch(vals)
                _,registered=mod.ScriptUtils._read_tsv(output)
                got=[parse(Float64,registered[1][prefix*lpad(string(i),3,'0')]) for i in eachindex(want)]
                @test maximum(abs,got-want)<5e-6
                # Translating the already-subtracted residual is NOT this model.
                wrong=fill(NaN,size(z)); wrong[:,1:end-2] .= residual[:,3:end]
                wrongvals=[mod._interp(xs,ys,wrong,1.2+t,1.2+u) for u in coords for t in coords]
                @test maximum(abs,want-mod._normalize_patch(wrongvals))>1e-3
            end
            @test_throws ErrorException mod._parse_cli(vcat(args,["--acquisition-shifts"]))
            @test_throws ErrorException mod._parse_cli(vcat(registered_args,["--acquisition-shifts",shifts]))

            # Local axes affect sampling only, including the difference views.
            _, base_rows = mod.ScriptUtils._read_tsv(feature_path)
            frame_rows = mod.PatchFrames.local_frame_rows(base_rows,(degree=2,))
            frames = mod.PatchFrames.write_frames(joinpath(dir,"$(name)_frames.tsv"),frame_rows)
            frame_args = vcat(args,["--config",COUNT_CONFIG,"--assignment-config",MATCHED_CONFIG,"--patch-frames",frames])
            extract_quietly(mod,frame_args)
            @test read(output) == matched_bytes # identity is byte-identical
            @test_throws ErrorException mod._parse_cli(vcat(args,["--patch-frames"]))
            @test_throws ErrorException mod._parse_cli(vcat(frame_args,["--patch-frames",frames]))
            angle = pi/6
            frame_rows[1]["frame_axis_x"] = string(cos(angle))
            frame_rows[1]["frame_axis_y"] = string(sin(angle))
            frame_rows[1]["tangent_slope"] = string(tan(angle))
            mod.PatchFrames.write_frames(frames,frame_rows)
            for policy in (ASSIGNMENT_CONFIG,MATCHED_CONFIG)
                extract_quietly(mod,vcat(args,["--config",COUNT_CONFIG,"--assignment-config",policy,"--patch-frames",frames]))
                _, measured_rows = mod.ScriptUtils._read_tsv(output)
                _, _, _, zf, sf, _, _ = mod.preprocess_channel(img,mod.get_channel(img,"Z";direction="fwd"),cfg)
                _, _, _, zb, sb, _, _ = mod.preprocess_channel(img,mod.get_channel(img,"Z";direction="bwd"),cfg)
                matched = policy == MATCHED_CONFIG
                rf = matched ? _box_smooth(zf-model,1) : sf-model
                rb = matched ? _box_smooth(zb-model,1) : sb-model
                sources = name == "forward" ? [("raw_p",sf),("res_p",rf)] :
                    [("bwd_raw_p",sb),("bwd_res_p",rb),("diff_raw_p",sf-sb),("diff_res_p",rf-rb)]
                for (pixel_prefix,source) in sources
                    samples = [mod._interp(xs,ys,source,1.2+t*cos(angle)-u*sin(angle),
                        1.2+t*sin(angle)+u*cos(angle)) for u in coords for t in coords]
                    want = mod._normalize_patch(samples)
                    got = [parse(Float64,measured_rows[1][pixel_prefix*lpad(string(i),3,'0')]) for i in eachindex(want)]
                    @test maximum(abs,got-want) < 5e-6
                end
                rotated_model = [begin
                    dt = (x-1.2)*cos(angle)+(y-1.2)*sin(angle)
                    du = -(x-1.2)*sin(angle)+(y-1.2)*cos(angle)
                    .4+.17x+.08y+2exp(-.5*(dt/.2)^2-.5*(du/.15)^2)
                end for y in ys,x in xs]
                wrong = matched ? _box_smooth(z-rotated_model,1) : _box_smooth(z,1)-rotated_model
                wrong_samples = [mod._interp(xs,ys,wrong,1.2+t*cos(angle)-u*sin(angle),
                    1.2+t*sin(angle)+u*cos(angle)) for u in coords for t in coords]
                actual_res = [parse(Float64,measured_rows[1][prefix*lpad(string(i),3,'0')]) for i in eachindex(wrong_samples)]
                @test maximum(abs,actual_res-mod._normalize_patch(wrong_samples)) > 1e-3
            end
            # Bad frame coverage fails before an existing output can be clobbered.
            write(output,"frame sentinel")
            mod.PatchFrames.write_frames(frames,Dict{String,String}[])
            @test_throws ErrorException extract_quietly(mod,frame_args)
            @test read(output,String) == "frame sentinel"
            # A fitted MODEL rotation changes subtraction, NOT patch sampling.
            model_feature_path = joinpath(dir,"$(name)_model_features.tsv")
            model_rows = deepcopy(base_rows)
            original_header = first(mod.ScriptUtils._read_tsv(feature_path))
            model_header = vcat(original_header,["model_orientation","model_axis_x","model_axis_y"])
            function save_model_rows()
                open(model_feature_path,"w") do io
                    println(io,join(model_header,'\t'))
                    for r in model_rows; println(io,join([r[k] for k in model_header],'\t')); end
                end
            end
            merge!(model_rows[1],Dict("model_orientation"=>"global","model_axis_x"=>"1","model_axis_y"=>"0"))
            save_model_rows()
            model_args = ["--features",model_feature_path,"--data-dir",dir,"--out",output,
                "--config",COUNT_CONFIG,"--assignment-config",MATCHED_CONFIG]
            extract_quietly(mod,model_args)
            @test read(output)==matched_bytes
            merge!(model_rows[1],Dict("model_orientation"=>"local_tangent",
                "model_axis_x"=>string(cos(angle)),"model_axis_y"=>string(sin(angle))))
            save_model_rows(); extract_quietly(mod,model_args)
            _, actual_model = mod.ScriptUtils._read_tsv(output)
            rotated_model = [begin
                dt = (x-1.2)*cos(angle)+(y-1.2)*sin(angle)
                du = -(x-1.2)*sin(angle)+(y-1.2)*cos(angle)
                .4+.17x+.08y+2exp(-.5*(dt/.2)^2-.5*(du/.15)^2)
            end for y in ys,x in xs]
            rotated_residual = _box_smooth(z-rotated_model,1)
            want = mod._normalize_patch([mod._interp(xs,ys,rotated_residual,1.2+t,1.2+u) for u in coords for t in coords])
            got = [parse(Float64,actual_model[1][prefix*lpad(string(i),3,'0')]) for i in eachindex(want)]
            @test maximum(abs,got-want)<5e-6
            @test maximum(abs,got-expected)>1e-3 # not the original global subtraction
            for column in original_header
                @test model_rows[1][column] == base_rows[1][column]
            end
            for column in keys(actual_model[1])
                startswith(column,name=="forward" ? "raw_p" : "bwd_raw_p") || continue
                @test actual_model[1][column] == old[column]
            end
            for (field,value) in (("model_axis_x","NaN"),("model_axis_x","2"),("model_orientation","unknown"))
                broken=deepcopy(model_rows); broken[1][field]=value
                @test_throws ErrorException mod.PatchFrames.read_model_axes(broken)
            end
            broken=deepcopy(model_rows); delete!(broken[1],"model_axis_y")
            @test_throws ErrorException mod.PatchFrames.read_model_axes(broken)
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

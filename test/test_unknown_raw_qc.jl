#!/usr/bin/env julia
# Synthetic, saved-geometry-only regression. No repository raw data is needed.
using Test, STMSXMIO
ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__,"plot_unknown_raw_qc.jl"))
const Q = UnknownRawQC
using Plots

function fixture_sxm(path; both=true, nonfinite=false)
    nx,ny = 5,4
    fwd = [Float64(3iy+ix)*1e-9 for iy in 1:ny, ix in 1:nx]
    backward = [Float64(7iy+2ix)*1e-9 for iy in 1:ny, ix in 1:nx]
    nonfinite && (fwd[2,3]=NaN; fwd[3,4]=Inf)
    directions = both ? "both" : "fwd"
    header = """
    :NANONIS_VERSION:
    2
    :SCAN_PIXELS:
    $nx $ny
    :SCAN_RANGE:
    4e-9 6e-9
    :SCAN_OFFSET:
    31e-9 -44e-9
    :SCAN_ANGLE:
    90
    :SCAN_DIR:
    down
    :SCAN_TIME:
    0.1 0.2
    :ACQ_TIME:
    3.0
    :REC_DATE:
    18.09.2026
    :REC_TIME:
    12:30:00
    :BIAS:
    -0.3
    :Lock-in>Lock-in status:
    OFF
    :DATA_INFO:
    Channel Name Unit Direction Calibration Offset
    14 Z m $directions 1 0
    0 Current A $directions 1 0
    :SCANIT_END:
    """
    # Nanonis x-fast data stream. Backward x acquisition order is reversed.
    matrices = both ? [fwd,reverse(backward;dims=2),fill(1e-11,ny,nx),fill(2e-11,ny,nx)] :
                      [fwd,fill(1e-11,ny,nx)]
    open(path,"w") do io
        write(io,header); write(io,UInt8[0x1a,0x04])
        for m in matrices, v in vec(permutedims(m))
            write(io,hton(reinterpret(UInt32,Float32(v))))
        end
    end
    return fwd,backward
end

function fixture_inputs(dir; raw=true, both=true)
    data = joinpath(dir,"raw"); mkpath(data)
    raw && fixture_sxm(joinpath(data,"synthetic.sxm");both)
    config = joinpath(dir,"config.toml")
    write(config,"[preprocessing]\nstride=1\nflatten=\"plane+rows\"\nsmooth_radius_px=1\n")
    geom = [Dict("file"=>"synthetic.sxm","N"=>"2","lobe"=>string(i),
        "x_nm"=>string(i),"y_nm"=>string(2i),"t_nm"=>string(i),"u_nm"=>string(2i),
        "axis_x"=>"1","axis_y"=>"0","origin_x_nm"=>"0","origin_y_nm"=>"0","source"=>"ell") for i in 1:2]
    counts = [Dict("filepath"=>"synthetic.sxm","N_selected"=>"2","refined_policy"=>"adaptive_support_rescue_keep",
        "support_2D_ell_nm"=>"2.5","support_2D_circ_nm"=>"2.5","artifact_fwd_bwd_corr"=>"0.4",
        "artifact_fwd_bwd_nrmse"=>"1.5")]
    # Deliberately include meaningless chemical labels: renderer must not read/use them.
    summaries = [Dict("file"=>"synthetic.sxm","N_selected"=>"2","assignment"=>"SECRET_LABELS")]
    reviews = [Dict("file"=>"synthetic.sxm","N_predicted_lobes"=>"2","review_status"=>"review",
        "review_reasons"=>"low_mean_confidence","mean_confidence"=>"0.2","uncertain_fraction"=>"0")]
    paths=String[]
    for (name,rows) in (("features",geom),("selected",counts),("summary",summaries),("review",reviews))
        p = joinpath(dir,name*".tsv"); push!(paths,p)
        Q.write_tsv(p,sort(collect(keys(first(rows)))),rows)
    end
    args = ["--features",paths[1],"--selected-summary",paths[2],"--summary",paths[3],
            "--review-queue",paths[4],"--data-dir",data,"--config",config,"--out-dir",joinpath(dir,"output")]
    return (; paths,config,data,args,geom)
end

@testset "unknown raw QC: coordinates, native IO, isolation" begin
    @test VERSION.major == 1 && VERSION.minor == 13
    mktempdir() do dir
        path=joinpath(dir,"x.sxm")
        f,b=fixture_sxm(path)
        img=read_sxm(path)
        @test size(img.channels[1].data)==(4,5)
        @test img.channels[1].data ≈ f rtol=1e-6
        @test img.channels[2].data ≈ b rtol=1e-6
        @test all(isapprox.(img.range_nm, (4.0,6.0)))
        @test all(isapprox.(img.offset_nm, (31.0,-44.0)))
        @test img.header["SCAN_DIR"] == "down"
        @test img.header["SCAN_ANGLE"] == "90"
        cfg=Q.PatchPreprocessing.PreprocessingSettings(2,"none",0)
        v=Q.prepare_view(img,"fwd",cfg)
        @test v.raw_x == collect(0.0:1.0:4.0)
        @test v.raw_y == collect(0.0:2.0:6.0)
        @test v.xs == [0.0,2.0,4.0]
        @test v.ys == [0.0,4.0]
        @test size(v.raw)==(4,5) # raw remains full resolution
        @test size(v.smooth)==(2,3)
        @test v.raw ≈ f*1e9 rtol=1e-6
        @test v.smooth ≈ f[1:2:end,1:2:end]*1e9 rtol=1e-6
        @test v.unit == "nm"
        @test v.raw_nonfinite == 0
        vb=Q.prepare_view(img,"bwd",cfg)
        @test vb.raw ≈ b*1e9 rtol=1e-6
        @test vb.raw != reverse(b;dims=2)*1e9 # no accidental second flip
        inv=Q.header_inventory(img)
        @test inv.scan["bias_V"] == -0.3
        @test inv.scan["native_dx_nm"] ≈ 1.0
        @test inv.scan["native_dy_nm"] ≈ 2.0
        @test inv.scan["scan_angle_header"] == "90"
        @test inv.scan["lockin_status"] == "OFF"
        @test length(inv.channels)==4
        @test [r["recorded_direction_expanded"] for r in inv.channels] == ["fwd","bwd","fwd","bwd"]
        @test Set(r["recorded_unit"] for r in inv.channels)==Set(["m","A"])
        @test inv.channels[2]["reader_transform"]=="x-reversed from acquisition bytes"
        @test inv.channels[1]["finite_pixels"] == 20
        @test inv.channels[1]["median_recorded_unit"] == Q.median(vec(img.channels[1].data))
        @test inv.channels[1]["sample_sd_recorded_unit"] == Q.std(vec(img.channels[1].data);corrected=true)
        @test inv.channels[3]["min_recorded_unit"] == inv.channels[3]["median_recorded_unit"] == inv.channels[3]["max_recorded_unit"]
        @test inv.channels[3]["sample_sd_recorded_unit"] == 0
        @test any(r->r["key"]=="DATA_INFO" && occursin("Current",r["value"]),inv.headers)
        @test !any(r->occursin("Lock",r["channel_name"]),inv.channels)
        @test Q.exact_channel(img,"Z","anything")===nothing

        # Preprocessing must be the native function (shared helpers), not a clone.
        cfg2=Q.PatchPreprocessing.PreprocessingSettings(1,"plane+rows",1)
        v2=Q.prepare_view(img,"fwd",cfg2)
        native_cfg=Q.PatternConfig(filepath=path,channel="Z",direction="fwd",stride=1,
            flatten="plane+rows",smooth_radius_px=1,no_plot=true)
        nx,ny,nraw,nflat,nsmooth,nu,nn=Q.preprocess_channel(img,img.channels[1],native_cfg)
        @test v2.xs==nx && v2.ys==ny
        @test v2.flat==nflat && v2.smooth==nsmooth
        @test v2.raw==nraw
        @test v2.noise==nn

        # Unequal x/y scales, rotated axes, nonzero origin and signed u.
        # Plot saved xy directly; t/u inverse is only an audit.
        rows=[Dict("lobe"=>"1","x_nm"=>"1.4","y_nm"=>"3.3","t_nm"=>"1.2","u_nm"=>"-0.8",
            "axis_x"=>"0.6","axis_y"=>"0.8","origin_x_nm"=>"0.04","origin_y_nm"=>"2.82")]
        @test Q.geometry_check(rows,img).max_projection_error_nm < 1e-14
        @test Q.geometry_check(rows,img).centers_outside==0
        rows[1]["x_nm"]="4.5"
        @test Q.geometry_check(rows,img).centers_outside==1
        rows[1]["x_nm"]="1.4"
        centers=Q.saved_centers(rows)
        p=Q.image_panel(v2,:raw,img.range_nm;centers)
        @test p.series_list[1][:x] == v2.raw_x
        @test p.series_list[1][:y] == v2.raw_y
        @test p.series_list[1][:z].surf == v2.raw
        @test p.series_list[2][:x] == [1.4]
        @test p.series_list[2][:y] == [3.3]
        @test p.subplots[1][:aspect_ratio] == :equal
        @test p.subplots[1][:colorbar] == :right
        @test p.subplots[1][:left_margin] == 12Plots.mm
        @test p.subplots[1][:right_margin] == 16Plots.mm
        @test p.subplots[1][:yaxis][:flip] == false
        blind=Q.image_panel(v2,:raw,img.range_nm)
        @test length(blind.series_list)==1
        @test isempty(blind.subplots[1][:annotations])
        smooth=Q.image_panel(v2,:smooth,img.range_nm)
        @test smooth.series_list[1][:z].surf==v2.smooth
        @test Q.display_limits(fill(1.0,2,2))[1]<1<Q.display_limits(fill(1.0,2,2))[2]
        @test_throws ErrorException Q.display_limits(fill(NaN,2,2))
        closeall()

        fixture_sxm(path;both=false)
        img=read_sxm(path)
        @test get_channel(img,"Z";direction="bwd").direction=="fwd" # native fallback exists
        @test Q.prepare_view(img,"bwd",cfg)===nothing # QC explicitly refuses fallback
        @test length(Q.header_inventory(img).channels)==2

        fixture_sxm(path;nonfinite=true)
        img=read_sxm(path)
        v=Q.prepare_view(img,"fwd",cfg2)
        @test isnan(v.raw[2,3]) && isinf(v.raw[3,4])
        @test all(isfinite,v.smooth)
        @test v.raw_nonfinite==2
        inv_nan=Q.header_inventory(img)
        @test inv_nan.channels[1]["finite_pixels"]==18
        @test isfinite(inv_nan.channels[1]["sample_sd_recorded_unit"])
        @test Q.centers_without_finite_raw_pixel(v,[(x=2.0,y=2.0),(x=1.0,y=0.0),(x=-1.0,y=0.0)])==2
        @test Q.centers_without_finite_raw_pixel(nothing,[(x=2.0,y=2.0)])==1
        nanraw=Q.image_panel(v,:raw,img.range_nm)
        nanprep=Q.image_panel(v,:smooth,img.range_nm)
        @test occursin("Missing raw pixels: 2",nanraw.subplots[1][:title])
        @test occursin("native imputation of 10.00%",nanprep.subplots[1][:title])
        closeall()
    end
end

@testset "unknown raw QC: keyed metadata and all-file failure retention" begin
    mktempdir() do dir
        f=fixture_inputs(dir)
        data=Q.load_inputs(f.paths...)
        @test data.files==["synthetic.sxm"]
        @test length(data.byfile["synthetic.sxm"])==2
        @test !haskey(data.counts["synthetic.sxm"],"assignment")
        @test_throws ErrorException Q.parse_cli(["--bogus","x"])
        @test_throws ErrorException Q.parse_cli(["--features"])
        @test_throws ErrorException Q.parse_cli(vcat(f.args,["--config",f.config]))
        @test Q.parse_cli(["--help"])===nothing
        @test_throws ErrorException Q.main(vcat(f.args,["--file","not_selected.sxm"]))
        @test !isdir(joinpath(dir,"output"))
        source=read(f.paths[1],String)
        bad=deepcopy(f.geom); bad[2]["lobe"]="1"
        Q.write_tsv(f.paths[1],sort(collect(keys(first(bad)))),bad)
        @test_throws ErrorException Q.load_inputs(f.paths...)
        write(f.paths[1],source)
        before=Dict(p=>Q.filehash(p) for p in vcat(f.paths,[f.config,joinpath(f.data,"synthetic.sxm")]))
        result=Q.main(f.args)
        out=joinpath(dir,"output")
        @test length(result.qc_rows)==1
        @test result.qc_rows[1]["status"]=="ok"
        @test result.qc_rows[1]["saved_N_algorithmic"]=="2"
        @test result.qc_rows[1]["support_policy_saved"]=="adaptive_support_rescue_keep"
        @test result.qc_rows[1]["max_projection_error_nm"]==0
        @test all(Q.filehash(p)==h for (p,h) in before)
        for name in ("synthetic_blind.png","synthetic_overlay.png")
            bytes=read(joinpath(out,name))
            @test bytes[1:8]==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a]
            width=foldl((x,y)->256x+Int(y),bytes[17:20];init=0)
            height=foldl((x,y)->256x+Int(y),bytes[21:24];init=0)
            @test 1798<=width<=1800 && 1598<=height<=1600
            @test length(bytes)>10_000
        end
        @test all(!occursin("SECRET_LABELS",read(joinpath(out,p),String)) for p in readdir(out) if !endswith(p,".png"))
        @test occursin("algorithmic",read(joinpath(out,"README.md"),String))
        @test occursin("not the exact",read(joinpath(out,"README.md"),String))
        @test occursin("no extra registration",read(joinpath(out,"qc_index.tsv"),String))
        @test length(readlines(joinpath(out,"display_limits.tsv")))==5
        @test_throws ErrorException Q.main(f.args) # no output overwrite
        @test Q.escape_tsv("a\nb\tc\\d")=="a\\nb\\tc\\\\d"
    end
    mktempdir() do dir
        f=fixture_inputs(dir;raw=false)
        @test_throws ErrorException Q.main(f.args)
        out=joinpath(dir,"output")
        qc=Q.read_columns(joinpath(out,"qc_index.tsv"),["file","status","error"])
        @test only(qc)["file"]=="synthetic.sxm"
        @test only(qc)["status"]=="error"
        @test !isempty(only(qc)["error"])
        @test isfile(joinpath(out,"synthetic_blind.png"))
        @test isfile(joinpath(out,"synthetic_overlay.png"))
        @test occursin("MISSING",read(joinpath(out,"source_hashes.tsv"),String))
    end
    mktempdir() do dir
        f=fixture_inputs(dir;both=false)
        @test_throws ErrorException Q.main(f.args)
        out=joinpath(dir,"output")
        qc=Q.read_columns(joinpath(out,"qc_index.tsv"),["file","status","error"])
        @test only(qc)["status"]=="missing_view"
        @test occursin("bwd not recorded",only(qc)["error"])
        @test length(readlines(joinpath(out,"channels.tsv")))==3
        @test length(readlines(joinpath(out,"display_limits.tsv")))==3
    end
end

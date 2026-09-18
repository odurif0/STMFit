#!/usr/bin/env julia
# Synthetic arrays/artifacts only. No optimizer, labels, or real image files.
# julia --startup-file=no --threads=1 --project=. test/test_acquisition_noise_diagnostics.jl
using Test, Random, Statistics, TOML, STMSXMIO
module AcquisitionRunner
include(joinpath(@__DIR__,"diagnose_acquisition_noise.jl"))
end
const AD=AcquisitionRunner.AcquisitionNoiseDiagnostics

function synthetic_settings(; kwargs...)
    d=Dict{String,Any}(
        "channel"=>"Z","auxiliary_channel"=>"Current","lag_x_min_px"=>-4,"lag_x_max_px"=>4,
        "lag_y_min_px"=>-2,"lag_y_max_px"=>2,
        "lag_bound_basis"=>"Fixed small synthetic search; no automatic expansion",
        "min_registration_pixels"=>100,"min_positive_correlation"=>0.2,
        "ambiguity_correlation_gap"=>0.001,"min_row_pixels"=>20,
        "footprint_sigma"=>3.0,"background_guard_nm"=>0.1,
        "min_background_pixels"=>200,"min_background_fraction"=>0.05,
        "acf_max_lag_px"=>12,"min_acf_pairs"=>100,"acf_threshold"=>exp(-1),
        "block_length_multiplier"=>2.0,"min_usable_blocks"=>4,
        "block_min_occupancy"=>0.8,"local_patch_radius_nm"=>0.4,"min_patch_pixels"=>20)
    for (k,v) in kwargs
        d[string(k)]=v
    end
    return AD.settings_from_dict(d)
end

# Independent generator with no wrapping/interpolation and an explicit sign.
function shifted_image(a,dx,dy;gain=1.0,offset=0.0)
    b=fill(NaN,size(a))
    for y in axes(a,1), x in axes(a,2)
        sy,sx=y-dy,x-dx
        if 1<=sy<=size(a,1) && 1<=sx<=size(a,2)
            b[y,x]=gain*a[sy,sx]+offset
        end
    end
    return b
end

@testset "native Julia 1.13 and explicit settings" begin
    @test v"1.13.0"<=VERSION<v"1.14.0"
    @test Threads.nthreads()==1
    s=synthetic_settings()
    @test s["channel"]=="Z"
    @test s["auxiliary_channel"]=="Current"
    for key in AD.REQUIRED_SETTINGS
        incomplete=copy(s); delete!(incomplete,key)
        @test_throws ErrorException AD.settings_from_dict(incomplete)
    end
    for (key,value) in (("lag_x_min_px",1),("lag_y_max_px",-1),
                        ("min_registration_pixels",true),("min_registration_pixels",2),
                        ("footprint_sigma",NaN),("footprint_sigma",0),
                        ("background_guard_nm",-1),("acf_threshold",1.0),
                        ("block_length_multiplier",0.5),("block_min_occupancy",1.1),
                        ("min_positive_correlation",0.0),("ambiguity_correlation_gap",-0.1),
                        ("min_background_fraction",0),("lag_bound_basis",""))
        invalid=copy(s); invalid[key]=value
        @test_throws ErrorException AD.settings_from_dict(invalid)
    end
    @test synthetic_settings(lag_bound_basis="Explicit finite synthetic window")["lag_bound_basis"]=="Explicit finite synthetic window"
    @test_throws ErrorException AcquisitionRunner.acquisition_options(["--truth","anything"])
    @test_throws ErrorException AcquisitionRunner.acquisition_options(String[])
end

@testset "signed affine metrics: inversion is a negative control" begin
    a=collect(-5.0:0.1:5.0)
    b=2.5 .* a .+ 7
    m=AD.pair_metrics(a,b)
    @test m.correlation≈1
    @test m.gain≈2.5
    @test m.offset≈7
    @test m.positive_affine_nrmse<1e-14
    @test m.direct_nrmse>0
    inv=AD.pair_metrics(a,-b)
    @test inv.correlation≈-1
    @test inv.gain==0
    @test inv.positive_affine_nrmse≈1
    @test AD.pair_metrics(ones(10),ones(10)).status=="insufficient_or_constant"
    @test AD.pair_metrics([1.0],[1.0]).n==1
    @test_throws DimensionMismatch AD.pair_metrics(a,[1.0])
    @test_throws ErrorException AD.pair_metrics([NaN,1,2],[1.0,2,3])
end

@testset "known 2D shift, constant support and nonfinite exclusion" begin
    rng=MersenneTwister(918)
    a=randn(rng,70,95)
    b=shifted_image(a,2,-1;gain=1.7,offset=3.2)
    a[30,40]=NaN; b[45,60]=NaN
    original_a,original_b=copy(a),copy(b)
    mask=falses(size(a)); mask[5:65,10:85].=true
    s=synthetic_settings()
    reg=AD.registration_scan(a,b,mask,s)
    @test reg.usable
    @test (reg.dx,reg.dy)==(2,-1)
    @test reg.best_correlation≈1
    @test length(reg.candidates)==45
    @test Set(r.n for r in reg.candidates)==Set([count(reg.support)])
    @test !reg.support[30,40]
    # One invalid backward sample excludes every forward location which can hit it.
    @test all(!reg.support[45-dy,60-dx] for dy in -2:2,dx in -4:4)
    @test isequal(a,original_a) && isequal(b,original_b)
    @test AD.pair_metrics(AD.sample_pairs(a,b,reg.support,reg.dx,reg.dy)...).positive_affine_nrmse<1e-13
    small=falses(size(a)); small[5:7,5:7].=true
    insufficient=AD.registration_scan(a,b,small,s)
    @test insufficient.status=="insufficient_support"
    @test isempty(insufficient.candidates)
    @test !insufficient.usable
    @test_throws DimensionMismatch AD.common_support(a,b,trues(2,2),-1:1,-1:1)
    @test_throws ErrorException AD.common_support(a,b,mask,Int[],0:0)
    allnan=AD.registration_scan(fill(NaN,20,20),zeros(20,20),trues(20,20),s)
    @test allnan.status=="insufficient_support"
end

@testset "boundary, ambiguity, inversion and independent null" begin
    rng=MersenneTwister(1909)
    a=randn(rng,80,90)
    mask=trues(size(a))
    boundary=AD.registration_scan(a,shifted_image(a,4,0),mask,synthetic_settings())
    @test boundary.dx==4
    @test boundary.boundary && !boundary.usable
    @test occursin("boundary",boundary.status)
    periodic=[isodd(x) ? 1.0 : -1.0 for y in 1:50,x in 1:60]
    ambiguous=AD.registration_scan(periodic,periodic,trues(size(periodic)),
        synthetic_settings(lag_y_min_px=0,lag_y_max_px=0))
    @test ambiguous.ambiguous && !ambiguous.usable
    @test ambiguous.near_maxima>=3
    zero=synthetic_settings(lag_x_min_px=0,lag_x_max_px=0,lag_y_min_px=0,lag_y_max_px=0)
    inverse=AD.registration_scan(a,-a,mask,zero)
    @test inverse.best_correlation≈-1
    @test occursin("nonpositive",inverse.status) && !inverse.usable
    null=AD.registration_scan(a,randn(rng,size(a)),mask,zero)
    @test abs(null.best_correlation)<0.05
    @test !null.usable
    constant=AD.registration_scan(zeros(40,40),ones(40,40),trues(40,40),zero)
    @test constant.status=="constant_signal"
end

@testset "background scale is off footprint, not fitted residual noise" begin
    rng=MersenneTwister(520)
    z=2 .* randn(rng,200,200)
    mask=trues(size(z)); mask[70:130,70:130].=false
    z[.!mask].+=1000
    s=synthetic_settings()
    bg=AD.background_diagnostic(z,mask,s)
    @test abs(bg.mad_sigma-2)<0.08
    @test abs(bg.standard_deviation-2)<0.08
    @test bg.n==count(mask)
    @test bg.length_x.pixels==1 && bg.length_y.pixels==1
    @test length(bg.acf)==2*(s["acf_max_lag_px"]+1)
    @test all(r->r.pairs>0,bg.acf)
    @test all(r->r.correlation≈1,filter(r->r.lag_px==0,bg.acf))
    bg2=AD.background_diagnostic(z .+ 100 .* .!mask,mask,s)
    @test bg2.mad_sigma==bg.mad_sigma
    insufficient=AD.background_diagnostic(z,falses(size(z)),s)
    @test insufficient.status=="insufficient_background"
    @test isnan(insufficient.mad_sigma) && isempty(insufficient.acf)
    frac=synthetic_settings(min_background_fraction=0.99)
    @test AD.background_diagnostic(z,mask,frac).status=="insufficient_background"
    allnan=AD.background_diagnostic(fill(NaN,30,30),trues(30,30),s)
    @test allnan.status=="insufficient_background"
end

@testset "anisotropic ACF lengths, pair counts, censoring and block feasibility" begin
    rng=MersenneTwister(1984)
    z=randn(rng,180,220)
    for y in axes(z,1), x in 2:size(z,2)
        z[y,x]=0.8*z[y,x-1]+sqrt(1-0.8^2)*z[y,x]
    end
    for y in 2:size(z,1), x in axes(z,2)
        z[y,x]=0.3*z[y-1,x]+sqrt(1-0.3^2)*z[y,x]
    end
    s=synthetic_settings()
    bg=AD.background_diagnostic(z,trues(size(z)),s)
    @test 3<=bg.length_x.pixels<=7
    @test 1<=bg.length_y.pixels<=2
    @test bg.length_x.pixels>bg.length_y.pixels
    rows=bg.acf
    @test only(filter(r->r.axis=="x"&&r.lag_px==3,rows)).pairs==180*217
    linear=[Float64(x+y) for y in 1:60,x in 1:70]
    censored=AD.background_diagnostic(linear,trues(size(linear)),s)
    @test censored.length_x.status=="right_censored_at_max_lag"
    @test isnan(censored.length_x.pixels)
    alternating=[isodd(x) ? 1.0 : -1.0 for y in 1:60,x in 1:70]
    anti=AD.background_diagnostic(alternating,trues(size(alternating)),s)
    @test anti.length_x.status=="right_censored_at_max_lag"
    resurgent=[(axis="x",lag_px=i,correlation=c) for (i,c) in enumerate([0.7,0.1,-0.8,-0.7])]
    @test AD.acf_length(resurgent,"x",s["acf_threshold"]).status=="resurgent_or_oscillatory_acf_tail"
    @test !AD.block_feasibility(trues(60,70),NaN,1.0,s).feasible
    blocks=AD.block_feasibility(trues(60,70),5.0,2.0,s)
    @test blocks.block_x_px==10 && blocks.block_y_px==4
    @test blocks.usable_blocks==105 && blocks.feasible
    thin=falses(60,70); thin[30,:].=true
    @test !AD.block_feasibility(thin,5.0,2.0,s).feasible
    very_sparse=falses(50,50); very_sparse[1:3,1:3].=true
    acf=AD.lag_acf(randn(rng,50,50),very_sparse,s)
    @test all(r->isnan(r.correlation),acf)
end

function write_synthetic_sxm(path,f,b; range_nm=(5.0,4.0), backward=true, current=nothing)
    ny,nx=size(f)
    header=":SCAN_PIXELS:\n$nx $ny\n:SCAN_RANGE:\n$(range_nm[1]*1e-9) $(range_nm[2]*1e-9)\n:SCAN_OFFSET:\n0 0\n:DATA_INFO:\nChannel Name Unit Direction Calibration Offset\n0 Z m $(backward ? "both" : "fwd") 1 0\n"
    current!==nothing && (header*="1 Current A both 1 0\n")
    header*=":SCANIT_END:\n"
    arrays=backward ? [(f,1e-9),(reverse(b;dims=2),1e-9)] : [(f,1e-9)]
    current!==nothing && append!(arrays,[(current[1],1e-12),(reverse(current[2];dims=2),1e-12)])
    open(path,"w") do io
        write(io,header)
        for (a,scale) in arrays
            for y in 1:ny,x in 1:nx
                bits=reinterpret(UInt32,Float32(a[y,x]*scale))
                write(io,hton(bits))
            end
        end
    end
    return path
end

function synthetic_artifacts(dir)
    rng=MersenneTwister(1201)
    nx,ny=96,80
    xs,ys=range(0,5;length=nx),range(0,4;length=ny)
    f=[0.7exp(-0.5*((x-2.0)^2/0.22^2+(y-2.0)^2/0.2^2)) +
       0.9exp(-0.5*((x-2.6)^2/0.22^2+(y-2.0)^2/0.2^2)) for y in ys,x in xs]
    f .+= 0.01randn(rng,ny,nx)
    b=f .+ 0.01randn(rng,ny,nx)
    f[15,20]=NaN
    raw=write_synthetic_sxm(joinpath(dir,"synthetic.sxm"),f,b)
    geo=joinpath(dir,"geometry.tsv")
    header=["file","N","lobe","amplitude","x_nm","y_nm","sigma_parallel_nm","sigma_perp_nm",
            "axis_x","axis_y","origin_x_nm","origin_y_nm","baseline","tilt_x","tilt_y"]
    open(geo,"w") do io
        println(io,join(header,'\t'))
        for (i,x,amp) in ((1,2.0,0.7),(2,2.6,0.9))
            println(io,join(["synthetic.sxm",2,i,amp,x,2.0,0.22,0.2,1.0,0.0,2.3,2.0,0.0,0.0,0.0],'\t'))
        end
    end
    physical=joinpath(dir,"physical.toml")
    cfg=Dict("model"=>Dict{String,Any}("spacing_min_nm"=>0.35,"spacing_max_nm"=>0.75,
        "fit_width_nm"=>0.45,"support_noise_k"=>2.5,"support_padding_nm"=>0.25,
        "support_min_length_nm"=>0.4,"support_baseline_quantile"=>0.1,"max_overlap"=>0.6,
        "global_maxtime"=>1.0,"global_maxiter"=>1,"sigma_parallel_min_nm"=>0.1,
        "sigma_parallel_max_nm"=>0.5,"selection_policy"=>"adaptive_support_rescue",
        "adaptive_rescue_support_noise_k"=>1.5,"adaptive_rescue_support_padding_nm"=>0.75),
        "preprocessing"=>Dict{String,Any}("stride"=>1,"flatten"=>"plane+rows","smooth_radius_px"=>1))
    open(physical,"w") do io
        TOML.print(io,cfg)
    end
    summary=joinpath(dir,"selected.tsv")
    write(summary,"filepath\tN_selected\trefined_policy\nsynthetic.sxm\t2\tadaptive_support_rescue_robust_guard\n")
    settings=joinpath(dir,"settings.toml")
    open(settings,"w") do io
        TOML.print(io,Dict("acquisition_noise"=>synthetic_settings()))
    end
    return (;raw,geo,physical,summary,settings,f,b)
end

@testset "native SXM and frozen geometry/support replay, no fitting" begin
    mktempdir() do dir
        a=synthetic_artifacts(dir)
        s=AD.load_settings(a.settings)
        img=STMSXMIO.read_sxm(a.raw)
        # The generator wrote acquisition-order backward; read_sxm flips exactly once.
        @test img.channels[2].data[35,26]*1e9≈a.b[35,26] rtol=1e-6
        @test_throws ErrorException AD.strict_channel(img,"missing","bwd")
        only_forward=write_synthetic_sxm(joinpath(dir,"forwardonly.sxm"),a.f,a.b;backward=false)
        @test_throws ErrorException AD.strict_channel(STMSXMIO.read_sxm(only_forward),"Z","bwd")
        case=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        @test case.selected==(n=2,use_rescue=true)
        @test case.effective_support_noise_k==1.5
        @test case.effective_support_padding_nm==0.75
        @test case.imputed_fwd==1 && case.imputed_bwd==0
        @test isnan(case.forward[15,20])
        @test !case.fused_valid[15,20]
        @test count(case.fit_support)>100
        @test count(case.fit_support)<=count(case.roi)
        before=copy(case.forward)
        d=AD.diagnose(case,s)
        @test isequal(case.forward,before)
        @test all(r->r.n==d.registration.n,d.registration.candidates)
        @test length(d.local_rows)==4
        @test !isempty(d.blocked_reasons)
        @test any(r->r.source=="fitted_residual_not_noise",d.noise_rows)
        @test any(r->r.source=="off_footprint_image",d.noise_rows)
        out=joinpath(dir,"output")
        inputs=Dict("raw"=>a.raw,"geometry"=>a.geo,"physical"=>a.physical,"summary"=>a.summary,"settings"=>a.settings)
        @test AD.write_diagnostic(out,case,d,s;inputs=inputs)==out
        @test sort(readdir(out))==sort(["registration_candidates.tsv","row_lags.tsv","local_views.tsv","noise_acf.tsv","summary.tsv","metadata.toml","auxiliary_views.tsv"])
        meta=TOML.parsefile(joinpath(out,"metadata.toml"))
        @test meta["selected_context"]["use_rescue"]
        @test meta["auxiliary_diagnostic"]["nuisance_projection"]=="not_performed"
        @test meta["auxiliary_diagnostic"]["forward"]["status"]=="missing_channel"
        @test meta["effective_support_padding_nm"]==0.75
        @test length(meta["inputs"]["raw"]["sha256"])==64
        @test !haskey(meta,"n_eff")
        @test length(meta["fit_support_sha256_column_major_bytes"])==64
        @test meta["fit_support_pixels"]==count(case.fit_support)
        out2=joinpath(dir,"cli_output")
        result=AcquisitionRunner.acquisition_main(["--file",a.raw,"--geometry",a.geo,"--config",a.physical,
            "--selected-summary",a.summary,"--settings",a.settings,"--outdir",out2])
        @test result.registration.status==d.registration.status
        for name in ("summary.tsv","registration_candidates.tsv","row_lags.tsv","local_views.tsv","noise_acf.tsv","auxiliary_views.tsv")
            @test read(joinpath(out,name))==read(joinpath(out2,name))
        end
        @test_throws ErrorException AD.write_diagnostic(out,case,d,s;inputs=inputs)
        link=joinpath(dir,"link"); symlink(out,link)
        @test_throws ErrorException AD.write_diagnostic(link,case,d,s;inputs=inputs)
        @test_throws ErrorException AcquisitionRunner.acquisition_options(["--file",a.raw,"--geometry",a.geo,"--config",a.physical,"--selected-summary",a.summary,"--settings",a.settings,"--outdir",out])
        # Ambiguous adaptive context must never infer selected support from N.
        write(a.summary,"filepath\tN_selected\nsynthetic.sxm\t2\n")
        @test_throws ErrorException AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        write(a.summary,"filepath\tN_selected\trefined_policy\nsynthetic.sxm\t2\tadaptive_support_rescue_keep\n")
        base=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        @test !base.selected.use_rescue
        @test base.effective_support_padding_nm==0.25
    end
end

@testset "saved geometry malformed/ambiguous inputs fail closed" begin
    mktempdir() do dir
        a=synthetic_artifacts(dir)
        original=read(a.geo,String)
        @test length(AD.read_geometry(a.geo,"synthetic.sxm",2))==2
        @test_throws ErrorException AD.read_geometry(a.geo,"synthetic.sxm",3)
        @test_throws ErrorException AD.read_geometry(a.geo,"absent.sxm",2)
        write(a.geo,original*split(original,'\n')[2]*"\n")
        @test_throws ErrorException AD.read_geometry(a.geo,"synthetic.sxm",2)
        for (from,to) in (("0.22","-0.22"),("0.22","NaN"),("axis_x","absent_axis"),("0.7","-0.7"))
            write(a.geo,replace(original,from=>to))
            @test_throws ErrorException AD.read_geometry(a.geo,"synthetic.sxm",2)
        end
        write(a.geo,replace(original,"tilt_y"=>"sequence"))
        @test_throws ErrorException AD.read_geometry(a.geo,"synthetic.sxm",2)
    end
end

@testset "Current direct coupling on exactly fixed Z support and lag" begin
    mktempdir() do dir
        a=synthetic_artifacts(dir)
        s=AD.load_settings(a.settings)
        f=copy(a.f); f[.!isfinite.(f)].=0.0
        cf=-2.5 .* f .+ 5.0
        cb=1.5 .* cf .+ 2.0
        write_synthetic_sxm(a.raw,f,f;current=(cf,cb))
        case=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        d=AD.diagnose(case,s)
        @test (d.registration.dx,d.registration.dy)==(0,0)
        @test case.auxiliary.forward.unit=="pA"
        @test case.auxiliary.backward.unit=="pA"
        @test case.auxiliary.forward.finite_pixels==length(f)
        @test case.auxiliary.forward.raw_std≈std(cf) rtol=1e-6
        @test length(d.auxiliary_rows)==7
        atlag=only(filter(r->r.comparison=="auxiliary_fwd_bwd" && r.condition=="at_fixed_Z_lag",d.auxiliary_rows))
        @test atlag.status=="ok"
        @test atlag.observed_pairs==atlag.z_support_pixels==d.registration.n
        @test (atlag.dx_px,atlag.dy_px)==(d.registration.dx,d.registration.dy)
        @test atlag.correlation≈1 atol=1e-10
        @test atlag.gain≈1.5 rtol=1e-5
        @test atlag.positive_affine_nrmse<1e-5
        coupling=filter(r->startswith(r.comparison,"Z_auxiliary"),d.auxiliary_rows)
        @test all(r->isapprox(r.correlation,-1;atol=1e-9),coupling)
        @test all(r->isapprox(r.linear_r2,1;atol=1e-9),coupling)
        @test all(r->isnan(r.positive_affine_nrmse)&&isnan(r.gain),coupling)
        # Fixed nonzero lag is inherited from Z, never re-estimated from Current.
        frozen=merge(d.registration,(dx=1,dy=-1))
        rows=AD.auxiliary_report(case,frozen,s)
        fixed=only(filter(r->r.comparison=="auxiliary_fwd_bwd"&&r.condition=="at_fixed_Z_lag",rows))
        @test (fixed.dx_px,fixed.dy_px)==(1,-1)
        @test fixed.z_support_pixels==d.registration.n
        # One missing Current sample must not trigger native median fill or
        # silently select an easier subset of the geometry for comparison.
        cf_bad=copy(cf); cf_bad[first(findall(d.registration.support))]=NaN
        write_synthetic_sxm(a.raw,f,f;current=(cf_bad,cb))
        bad=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        badrows=AD.auxiliary_report(bad,d.registration,s)
        invalid=only(filter(r->r.comparison=="auxiliary_fwd_bwd"&&r.condition=="at_fixed_Z_lag",badrows))
        @test invalid.status=="nonfinite_on_fixed_Z_support"
        @test invalid.observed_pairs==invalid.z_support_pixels-1
        @test isnan(invalid.correlation)
        # Zero channel carries no identifiable structure despite availability.
        write_synthetic_sxm(a.raw,f,f;current=(zeros(size(f)),zeros(size(f))))
        zero=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        zr=AD.auxiliary_report(zero,d.registration,s)
        @test zero.auxiliary.forward.status=="available"
        @test zero.auxiliary.forward.raw_std==0
        @test all(r->isnan(r.correlation),filter(r->r.comparison!="direction_scale",zr))
        @test all(r->r.status=="insufficient_or_constant",filter(r->r.comparison!="direction_scale",zr))
        # All-nonfinite or absent direction remains explicit, never a fallback.
        write_synthetic_sxm(a.raw,f,f;current=(fill(NaN,size(f)),cb))
        missing=AD.load_case(a.raw,a.geo,a.physical,a.summary,s)
        @test missing.auxiliary.forward.status=="all_nonfinite_channel"
        @test missing.auxiliary.forward.image===nothing
        mr=AD.auxiliary_report(missing,d.registration,s)
        @test any(r->r.status=="unavailable_direction",mr)
        img=STMSXMIO.read_sxm(a.raw)
        pcfg,_,_=AD.SavedSupport._configs(case.physical_config["model"],case.physical_config["preprocessing"],"";
            selected_context=case.selected)
        without_bwd=STMSXMIO.SXMImage(img.filepath,img.header,img.width,img.height,img.range_nm,img.offset_nm,
            filter(c->!(c.name=="Current"&&c.direction=="bwd"),img.channels))
        @test AD.load_auxiliary(without_bwd,pcfg,"Current","bwd").status=="missing_channel"
        duplicate=STMSXMIO.SXMImage(img.filepath,img.header,img.width,img.height,img.range_nm,img.offset_nm,
            vcat(img.channels,[img.channels[end]]))
        @test AD.load_auxiliary(duplicate,pcfg,"Current","bwd").status=="duplicate_direction_channels"
    end
end

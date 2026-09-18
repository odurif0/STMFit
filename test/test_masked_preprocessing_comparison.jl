#!/usr/bin/env julia
# Synthetic images/artifacts ONLY. No real input is read, fitted or diagnosed.
# julia --startup-file=no --threads=1 --project=. test/test_masked_preprocessing_comparison.jl
using Test, Random, Statistics, TOML, SHA, Serialization, LinearAlgebra
using STMSXMIO, GaussianFit2D
module ComparisonCLI
include(joinpath(@__DIR__, "diagnose_masked_preprocessing_real.jl"))
end
const MC = ComparisonCLI.MaskedPreprocessingComparison
const AD = MC.AD
const MRP = MC.MRP
const ROOT = normpath(joinpath(@__DIR__, ".."))
const MASK_PATH = joinpath(ROOT,"config","masked_robust_preprocessing.toml")
const ACQ_PATH = joinpath(ROOT,"config","label_free_exploration.toml")
const LOCK_PATHS = [MASK_PATH, ACQ_PATH, joinpath(ROOT,"Project.toml"),joinpath(ROOT,"Manifest.toml"),
    joinpath(ROOT,"config","chitosan.toml"),joinpath(ROOT,"test","lib","masked_robust_preprocessing.jl"),
    joinpath(ROOT,"test","lib","acquisition_noise_diagnostics.jl"),joinpath(ROOT,"test","extract_lobe_features.jl"),
    joinpath(ROOT,"packages","GaussianFit2D.jl","src","core.jl")]
const ORIGINAL_HASHES = Dict(p=>MC.filehash(p) for p in LOCK_PATHS)
BLAS.set_num_threads(1)
settings() = AD.load_settings(ACQ_PATH)
options() = MRP.read_options(TOML.parsefile(MASK_PATH))

function shifted(a, dx, dy)
    b = fill(NaN,size(a))
    for y in axes(a,1), x in axes(a,2)
        yy,xx = y-dy,x-dx
        1<=yy<=size(a,1) && 1<=xx<=size(a,2) && (b[y,x]=a[yy,xx])
    end
    return b
end

function rawrecord(z)
    return (z=copy(z), observed=BitMatrix(isfinite.(z)), fullimage_pixels=length(z),
        fullimage_observed=count(isfinite,z), fullimage_rows=size(z,1), fullimage_columns=size(z,2))
end

function native_preprocess(z)
    ny,nx = size(z)
    ch = SXMChannel("Z","m","fwd",z.*1e-9)
    img = SXMImage("synthetic://comparison",Dict{String,String}(),nx,ny,(5.0,4.0),(0.0,0.0),[ch])
    cfg = GaussianFit2D.PatternConfig(stride=1,flatten="plane+rows",smooth_radius_px=0)
    xs,ys,_,f,_,unit,_ = GaussianFit2D.preprocess_channel(img,ch,cfg)
    f[.!isfinite.(z)].=NaN
    return xs,ys,f,unit
end

function fixture(; badrows=false, missingness=true)
    rng=MersenneTwister(39421)
    rawf=0.04randn(rng,64,88)
    rawf .+= [0.1+0.03x+0.005y for y in 1:64,x in 1:88]
    rawb=shifted(rawf,2,-1) .+ 0.7
    if missingness
        rawf[5,:].=NaN
        rawf[35,40]=Inf
        rawb[40:41,55:57].=NaN
    end
    xs,ys,f,unit=native_preprocess(rawf)
    _,_,b,_=native_preprocess(rawb)
    fit=falses(size(f)); fit[10:54,16:70].=true
    roi=falses(size(f)); roi[20:43,32:54].=true
    lobes=[(lobe=i,x_nm=x,y_nm=2.0,amplitude=0.4,sigma_parallel_nm=0.16,sigma_perp_nm=0.12,
        axis_x=1.0,axis_y=0.0,origin_x_nm=2.5,origin_y_nm=2.0,baseline=0.0,tilt_x=0.0,tilt_y=0.0,skew_ratio=1.0)
        for (i,x) in enumerate((2.1,2.7))]
    case=(file="synthetic.sxm",xs=collect(xs),ys=collect(ys),forward=f,backward=b,
        fit_support=fit,roi=roi,lobes=lobes,unit=unit,selected=(n=2,use_rescue=false),
        support_meta=(support_method="synthetic_frozen",),effective_support_noise_k=2.5,
        effective_support_padding_nm=0.25)
    footprint=AD.geometry_footprint(case,settings())
    badrows && (footprint[24:25,:].=true)
    rawdata=(forward=rawrecord(rawf),backward=rawrecord(rawb))
    return (;case,rawdata,footprint)
end

function clone_views(case,f,b)
    rawdata=(forward=rawrecord(f),backward=rawrecord(b))
    views=Dict{String,Any}()
    for m in MC.METHODS
        # Deliberately injected corrected arrays test comparison semantics, not
        # estimator claims. Non-native slots keep a Boolean background field.
        views[m*"/fwd"]=merge(MC.native_view(f,rawdata.forward),
            (background_mask=copy(rawdata.forward.observed),))
        views[m*"/bwd"]=merge(MC.native_view(b,rawdata.backward),
            (background_mask=copy(rawdata.backward.observed),))
    end
    return views,rawdata
end

function write_synthetic_sxm(path,f,b;unit="m",range_nm=(5.0,4.0),backward=true)
    ny,nx=size(f)
    header=":SCAN_PIXELS:\n$nx $ny\n:SCAN_RANGE:\n$(range_nm[1]*1e-9) $(range_nm[2]*1e-9)\n:SCAN_OFFSET:\n0 0\n:DATA_INFO:\nChannel Name Unit Direction Calibration Offset\n0 Z $unit $(backward ? "both" : "fwd") 1 0\n:SCANIT_END:\n"
    arrays=backward ? [f,reverse(b;dims=2)] : [f]
    scale=unit=="m" ? 1e-9 : 1.0
    open(path,"w") do io
        write(io,header)
        for a in arrays, y in 1:ny,x in 1:nx
            write(io,hton(reinterpret(UInt32,Float32(a[y,x]*scale))))
        end
    end
    return path
end

function artifacts(dir; fake_raw=false)
    rng=MersenneTwister(1201)
    nx,ny=96,80
    xs,ys=range(0,5;length=nx),range(0,4;length=ny)
    f=[0.7exp(-0.5*((x-2.0)^2/0.22^2+(y-2.0)^2/0.2^2))+
       0.9exp(-0.5*((x-2.6)^2/0.22^2+(y-2.0)^2/0.2^2)) for y in ys,x in xs]
    f .+= 0.01randn(rng,ny,nx)
    b=f .+ 0.01randn(rng,ny,nx)
    f[15,20]=NaN
    raw=joinpath(dir,"synthetic.sxm")
    fake_raw ? write(raw,"THIS IS NOT SXM; metadata-only dry-run must not parse me") : write_synthetic_sxm(raw,f,b)
    geometry=joinpath(dir,"geometry.tsv")
    header=["file","N","lobe","amplitude","x_nm","y_nm","sigma_parallel_nm","sigma_perp_nm",
            "axis_x","axis_y","origin_x_nm","origin_y_nm","baseline","tilt_x","tilt_y"]
    open(geometry,"w") do io
        println(io,join(header,'\t'))
        for (i,x,amp) in ((1,2.0,0.7),(2,2.6,0.9))
            println(io,join(["synthetic.sxm",2,i,amp,x,2.0,0.22,0.2,1.0,0.0,2.3,2.0,0.0,0.0,0.0],'\t'))
        end
    end
    cfg=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    cfg["model"]["selection_policy"]="adaptive_support_rescue"
    cfg["model"]["adaptive_rescue_support_noise_k"]=1.5
    cfg["model"]["adaptive_rescue_support_padding_nm"]=0.75
    physical=joinpath(dir,"physical.toml")
    open(io->TOML.print(io,cfg),physical,"w")
    summary=joinpath(dir,"summary.tsv")
    write(summary,"filepath\tN_selected\trefined_policy\nsynthetic.sxm\t2\tadaptive_support_rescue_robust_guard\n")
    acquisition=joinpath(dir,"acquisition.toml"); cp(ACQ_PATH,acquisition)
    mask=joinpath(dir,"masked.toml"); cp(MASK_PATH,mask)
    out=joinpath(dir,"uncreated","out")
    args=["--file",raw,"--geometry",geometry,"--selected-summary",summary,"--config",physical,
        "--acquisition-settings",acquisition,"--settings",mask,"--outdir",out]
    return (;raw,geometry,summary,physical,acquisition,mask,out,args,f,b)
end

function all_plain(x)
    if x isa AbstractDict
        return all(k->k isa String,keys(x)) && all(all_plain,values(x))
    elseif x isa AbstractArray
        return all(all_plain,x)
    end
    return x isa Union{Bool,Integer,AbstractFloat,AbstractString}
end

@testset "masked real comparison (synthetic only)" begin
    @test v"1.13.0"<=VERSION<v"1.14.0"
    @test Threads.nthreads()==1
    @test BLAS.get_num_threads()==1
    @test MC.METHODS==("native_reference","finite_only_ols","guarded_ols","guarded_huber")

    @testset "metadata dry-run, strict CLI, original selected context" begin
        mktempdir() do dir
            a=artifacts(dir;fake_raw=true)
            before=Dict(p=>MC.filehash(p) for p in (a.raw,a.geometry,a.summary,a.physical,a.acquisition,a.mask))
            ctx=ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            @test ctx.selected==(n=2,use_rescue=true)
            @test ctx.effective_support_noise_k==1.5
            @test ctx.effective_support_padding_nm==0.75
            @test !ispath(a.out) && !ispath(dirname(a.out))
            @test read(a.raw,String)=="THIS IS NOT SXM; metadata-only dry-run must not parse me"
            @test all(MC.filehash(p)==h for (p,h) in before)
            @test ComparisonCLI.masked_comparison_options(a.args).dry_run==false
            @test_throws ErrorException ComparisonCLI.masked_comparison_options(String[])
            for args in ([a.args;"--file";a.raw],[a.args;"--unknown";a.raw],
                         [a.args;"--dry-run";"--dry-run"],[a.args;"--dry-run";"value"],
                         a.args[1:end-1],[a.args[1:end-1];"--dry-run"],[a.args[1:end-1];" "])
                @test_throws ErrorException ComparisonCLI.masked_comparison_options(args)
            end
            collision=[a.args[1:end-1];dir]
            @test_throws ErrorException ComparisonCLI.masked_comparison_options(collision)
            dangling=joinpath(dir,"dangling"); symlink(joinpath(dir,"missing"),dangling)
            @test_throws ErrorException MC.require_new_output(dangling)
            link=joinpath(dir,"linked"); symlink(dir,link)
            @test_throws ErrorException MC.require_new_output(joinpath(link,"new"))
            @test_throws ErrorException MC.require_new_output(joinpath(a.raw,"new"))
            io=IOBuffer(); ComparisonCLI.masked_comparison_help(io)
            help=String(take!(io))
            @test occursin("--acquisition-settings",help)
            @test occursin("--dry-run",help) && occursin("Julia 1.13",help)
            @test occursin("arrays.jls",help)
            original_summary=read(a.summary,String)
            write(a.summary,"filepath\tN_selected\nsynthetic.sxm\t2\n")
            @test_throws ErrorException ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            write(a.summary,replace(original_summary,"adaptive_support_rescue_robust_guard"=>"adaptive_support_rescue_keep"))
            keep=ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            @test !keep.selected.use_rescue
            @test keep.effective_support_padding_nm==0.25
            write(a.summary,replace(original_summary,"synthetic.sxm"=>"different.sxm"))
            @test_throws ErrorException ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            write(a.summary,original_summary)
            original_geo=read(a.geometry,String)
            write(a.geometry,original_geo*split(original_geo,'\n')[2]*"\n")
            @test_throws ErrorException ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            write(a.geometry,original_geo)
            cfg=TOML.parsefile(a.physical)
            for (key,value) in (("flatten","none"),("stride",true),("smooth_radius_px",-1))
                bad=deepcopy(cfg); bad["preprocessing"][key]=value
                open(io->TOML.print(io,bad),a.physical,"w")
                @test_throws ErrorException ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            end
            open(io->TOML.print(io,cfg),a.physical,"w")
            acq=TOML.parsefile(a.acquisition); acq["acquisition_noise"]["lag_x_max_px"]=9
            open(io->TOML.print(io,acq),a.acquisition,"w")
            @test_throws ErrorException ComparisonCLI.masked_comparison_main([a.args;"--dry-run"])
            @test !ispath(a.out)
        end
    end

    @testset "actual four estimators, coverage and fixed all/pair supports" begin
        t=fixture(badrows=true)
        original=deepcopy(t)
        views=MC.build_views(t.case,t.rawdata,t.footprint,options())
        @test length(views)==8
        d=MC.compare_views(t.case,t.rawdata,views,t.footprint,settings())
        @test isequal(t,original)
        @test length(d.view_rows)==8
        @test length(d.row_rows)==8*length(t.case.ys)
        @test length(d.summaries)==28
        @test length(d.candidates)==14*85
        @test length(d.local_rows)==14*2*length(t.case.lobes)
        @test length(d.background_rows)==14*2
        @test count(d.all_common)>100
        @test count(d.all_common)<count(d.own["native_reference"])
        for m in MC.METHODS, direction in MC.DIRECTIONS
            v=views[m*"/"*direction]
            @test !any(v.valid .& .!v.observed)
            @test v.valid==isfinite.(v.corrected)
            @test all(isnan,v.corrected[.!v.observed])
        end
        native=views["native_reference/fwd"]
        @test native.converged=="unknown" && isnan(native.iterations)
        @test isnan(native.coefficients.slope_x)
        @test all(isnan,native.row_offsets) && all(isnan,native.row_background_counts)
        for m in ("guarded_ols","guarded_huber"), direction in MC.DIRECTIONS
            v=views[m*"/"*direction]
            @test v.status==:PARTIAL
            @test all(==(:insufficient_background),v.row_status[24:25])
            @test !any(v.valid[24:25,:])
        end
        for scope in d.scopes
            @test all(!scope.support[I] || t.case.fit_support[I] for I in eachindex(scope.support))
            @test !any(scope.background_support .& t.footprint)
            if scope.kind=="all_common"
                @test scope.support==reduce((x,y)->x.&y,[d.own[m] for m in MC.METHODS])
            elseif scope.kind=="native_pair"
                @test scope.support==d.own[scope.methods[1]].&d.own[scope.methods[2]]
            end
            rows=filter(r->r.scope==scope.id,d.candidates)
            @test Set(r.n for r in rows)==Set([count(scope.support)])
            @test Set(r.support_sha256 for r in rows)==Set([MC.maskhash(scope.support)])
            @test Set(r.background_sha256 for r in rows)==Set([MC.maskhash(scope.background_support)])
            for m in scope.methods
                ownrows=filter(r->r.method==m,rows)
                @test length(ownrows)==85
                @test Set((r.dx_px,r.dy_px) for r in ownrows)==Set((dx,dy) for dx in -8:8 for dy in -2:2)
                @test length(unique(r.forward_level_nm for r in ownrows))==1
                @test length(unique(r.backward_level_nm for r in ownrows))==1
                f,b=views[m*"/fwd"].corrected,views[m*"/bwd"].corrected
                reg=AD.registration_scan(f,b,scope.support,settings())
                @test reg.support==scope.support
                for r in ownrows
                    metric=AD.pair_metrics(AD.sample_pairs(f,b,scope.support,r.dx_px,r.dy_px)...)
                    @test r.correlation==metric.correlation
                    @test r.difference_rms==metric.difference_rms
                end
            end
        end
        native_best=only(filter(r->r.scope=="own/native_reference"&&r.condition=="grid_best_diagnostic_only",d.summaries))
        @test (native_best.dx_px,native_best.dy_px)==(2,-1)
        @test all(!r.registration_accepted for r in d.summaries)
        for r in d.view_rows
            raw=r.direction=="fwd" ? t.rawdata.forward : t.rawdata.backward
            @test r.fullimage_missing_pixels==count(.!isfinite.(raw.z))
            @test r.raw_missing_pixels==r.fullimage_missing_pixels
        end
    end

    @testset "blocked Huber, missing rows and retained unevaluated grid" begin
        t=fixture(badrows=true)
        cfg=TOML.parsefile(MASK_PATH)
        cfg["masked_robust_preprocessing"]["maxiter"]=1
        forced=MRP.read_options(cfg)
        views=MC.build_views(t.case,t.rawdata,t.footprint,forced)
        @test views["guarded_huber/fwd"].status==:NONCONVERGED
        @test views["guarded_huber/bwd"].status==:NONCONVERGED
        d=MC.compare_views(t.case,t.rawdata,views,t.footprint,settings())
        @test count(d.all_common)==0
        @test count(d.own["native_reference"])>100
        @test count(d.own["guarded_ols"])>100
        @test length(d.view_rows)==8 && length(d.candidates)==1190
        blocked=filter(r->r.scope=="all_common",d.candidates)
        @test length(blocked)==340
        @test all(r->r.n==0 && isnan(r.correlation) && r.insufficient_support,blocked)
        @test all(r->r.status=="insufficient_support",blocked)
        best=filter(r->r.scope=="all_common"&&r.condition=="grid_best_diagnostic_only",d.local_rows)
        @test length(best)==4*length(t.case.lobes)
        @test all(r->isnan(r.dx_px)&&isnan(r.dy_px)&&r.status=="no_grid_maximum",best)
        pair=filter(r->r.scope=="native_pair/guarded_ols",d.candidates)
        @test length(pair)==170 && all(r->r.n>100&&isfinite(r.correlation),pair)
        @test all(r->r.all_common_pixels==0,pair)
        @test all(r->!r.registration_accepted,pair)
        all_excluded=trues(size(t.footprint))
        failed=MC.build_views(t.case,t.rawdata,all_excluded,options())
        @test failed["guarded_ols/fwd"].status==:UNAVAILABLE
        @test failed["finite_only_ols/fwd"].status==:PARTIAL
        df=MC.compare_views(t.case,t.rawdata,failed,all_excluded,settings())
        @test count(df.all_common)==0
        @test all(r->r.status=="insufficient_shared_observed_background",df.background_rows)
        @test all(r->isnan(r.level_aligned_rms_nm),df.candidates)
    end

    @testset "negative, boundary, ambiguous and constant signals remain visible" begin
        t=fixture(missingness=false)
        rng=MersenneTwister(763)
        z=randn(rng,size(t.case.forward))
        for (f,b,expected) in ((z,shifted(z,8,0),"boundary"),
                ([Float64(x+y) for y in 1:64,x in 1:88],[-Float64(x+y) for y in 1:64,x in 1:88],"negative"),
                ([isodd(x) ? 1.0 : -1.0 for y in 1:64,x in 1:88],[isodd(x) ? 1.0 : -1.0 for y in 1:64,x in 1:88],"ambiguous"),
                (ones(64,88),ones(64,88),"constant"))
            views,rawdata=clone_views(t.case,f,b)
            d=MC.compare_views(t.case,rawdata,views,t.footprint,settings())
            best=only(filter(r->r.scope=="all_common"&&r.method=="native_reference"&&r.condition=="grid_best_diagnostic_only",d.summaries))
            @test !best.registration_accepted
            @test !best.registration_usable_diagnostic
            if expected=="boundary"
                @test best.boundary && best.dx_px==8
                @test best.correlation≈1
            elseif expected=="negative"
                @test best.correlation≈-1
                @test occursin("nonpositive",best.registration_status)
                @test best.gain==0
                @test best.positive_affine_nrmse≈1
                @test any(r->r.correlation<0,d.local_rows)
            elseif expected=="ambiguous"
                @test best.ambiguous && best.near_maxima>1
            else
                @test best.status=="no_grid_maximum"
                @test isnan(best.dx_px) && isnan(best.dy_px)
                @test all(r->isnan(r.correlation),d.candidates)
                @test all(r->isnan(r.dx_px),filter(r->r.condition=="grid_best_diagnostic_only",d.local_rows))
            end
            @test length(d.candidates)==1190
        end
    end

    @testset "same observed background gauge, no truth or lag-specific level" begin
        t=fixture(missingness=false)
        rng=MersenneTwister(432)
        z=randn(rng,size(t.case.forward))
        views,rawdata=clone_views(t.case,z.+2,z.+7)
        d=MC.compare_views(t.case,rawdata,views,t.footprint,settings())
        zero=only(filter(r->r.scope=="all_common"&&r.method=="native_reference"&&r.condition=="zero_lag",d.summaries))
        @test zero.correlation≈1
        @test zero.difference_rms≈5
        @test zero.level_aligned_rms_nm<1e-14
        @test zero.backward_level_nm-zero.forward_level_nm≈5
        for r in d.candidates
            scope=only(filter(s->s.id==r.scope,d.scopes))
            f,b=views[r.method*"/fwd"].corrected,views[r.method*"/bwd"].corrected
            @test r.forward_level_nm==median(f[scope.background_support])
            @test r.backward_level_nm==median(b[scope.background_support])
        end
        # Changing only excluded foreground cannot affect the observable levels.
        altered=deepcopy(views)
        for m in MC.METHODS
            altered[m*"/fwd"].corrected[t.footprint].+=19
        end
        changed=MC.compare_views(t.case,rawdata,altered,t.footprint,settings())
        @test [r.forward_level_nm for r in d.candidates]==[r.forward_level_nm for r in changed.candidates]
        @test any(abs(d.candidates[i].difference_rms-changed.candidates[i].difference_rms)>1 for i in eachindex(d.candidates))
        for r in d.background_rows
            @test r.n>200 && r.fraction>=0.05
            @test r.change_level_aligned_rms_vs_native_nm==0
            @test occursin("descriptive_postfit",r.status)
        end
        # The fixed acquisition minimum background fraction/count gates levels,
        # not the already-valid correlation support. No relaxed replacement mask.
        small=trues(size(t.footprint)); small[1:2,1:20].=false
        sparse=MC.compare_views(t.case,rawdata,views,small,settings())
        @test all(r->r.anchor_pixels==40,sparse.candidates)
        @test all(r->r.anchor_status=="insufficient_shared_observed_background",sparse.candidates)
        @test all(r->isnan(r.level_aligned_rms_nm),sparse.candidates)
        @test any(r->isfinite(r.correlation),sparse.candidates)
        @test all(r->isnan(r.mean_nm)&&isnan(r.change_level_aligned_rms_vs_native_nm),sparse.background_rows)
    end

    @testset "native synthetic SXM, unchanged replay, plain Serialization and CLI" begin
        mktempdir() do dir
            a=artifacts(dir)
            ctx=MC.metadata_context(a.raw,a.geometry,a.physical,a.summary,a.acquisition,a.mask)
            ctx_before=deepcopy(ctx.s)
            loaded=MC.load_real_context(a.raw,a.geometry,a.physical,a.summary,ctx)
            case,rawdata,footprint=loaded.case,loaded.rawdata,loaded.footprint
            @test ctx.s==ctx_before
            @test !hasproperty(case,:auxiliary)
            native=AD.load_case(a.raw,a.geometry,a.physical,a.summary,ctx.s)
            @test isequal(case.forward,native.forward) && isequal(case.backward,native.backward)
            @test case.fit_support==native.fit_support && case.roi==native.roi
            @test case.selected==(n=2,use_rescue=true)
            @test footprint==AD.geometry_footprint(native,ctx.s)
            @test size(rawdata.forward.z)==size(a.f)
            @test rawdata.forward.observed==isfinite.(a.f)
            @test rawdata.backward.z≈a.b rtol=1e-6
            @test rawdata.backward.z[35,26]≈a.b[35,26] rtol=1e-6
            @test rawdata.forward.fullimage_observed==length(a.f)-1
            views=MC.build_views(case,rawdata,footprint,ctx.options)
            d=MC.compare_views(case,rawdata,views,footprint,ctx.s)
            reg=AD.registration_scan(native.forward,native.backward,native.fit_support,ctx.s)
            nr=filter(r->r.scope=="own/native_reference",d.candidates)
            @test length(nr)==length(reg.candidates)==85
            for (r,e) in zip(nr,reg.candidates)
                @test r.n==e.n && r.dx_px==e.dx_px && r.dy_px==e.dy_px
                @test isequal(r.correlation,e.correlation)
                @test isequal(r.positive_affine_nrmse,e.positive_affine_nrmse)
            end
            inputs=Dict("raw"=>a.raw,"geometry"=>a.geometry,"summary"=>a.summary,"physical"=>a.physical,
                "acquisition"=>a.acquisition,"masked"=>a.mask)
            @test MC.write_comparison(a.out,case,rawdata,views,footprint,d,ctx;inputs=inputs)==a.out
            expected=["view_status.tsv","row_status.tsv","comparison_summary.tsv","lag_candidates.tsv",
                "local_comparisons.tsv","background_summary.tsv","arrays.jls","metadata.toml"]
            @test sort(readdir(a.out))==sort(expected)
            snapshot=deserialize(joinpath(a.out,"arrays.jls"))
            @test all_plain(snapshot)
            saved_options=MRP.read_options(Dict("masked_robust_preprocessing"=>snapshot["masked_options"]))
            @test saved_options==ctx.options
            @test snapshot["masked_options"]["maxiter"] isa Integer
            @test snapshot["frozen"]["fit_support"]==case.fit_support
            @test isequal(snapshot["frozen"]["support_meta"],MC.plain(case.support_meta))
            @test snapshot["frozen"]["roi"]==case.roi
            @test snapshot["frozen"]["xs"]==collect(case.xs)
            @test snapshot["frozen"]["ys"]==collect(case.ys)
            @test snapshot["frozen"]["footprint"]==footprint
            @test isequal(snapshot["raw_scaled_analysis_nm"]["forward"]["z"],rawdata.forward.z)
            @test isequal(snapshot["views"]["guarded_huber/fwd"]["corrected"],views["guarded_huber/fwd"].corrected)
            @test snapshot["views"]["native_reference/fwd"]["converged"]=="unknown"
            @test !haskey(snapshot["frozen"],"auxiliary")
            # Reconstruct background/rows from serialized plain fields only.
            for method in MC.METHODS[2:end],direction in MC.DIRECTIONS
                v=snapshot["views"][method*"/"*direction]
                z=snapshot["raw_scaled_analysis_nm"][direction=="fwd" ? "forward" : "backward"]["z"]
                reconstructed=z.-v["plane"].-v["row_offsets"]
                @test all(isapprox.(reconstructed[v["valid"]],v["corrected"][v["valid"]];atol=1e-14,rtol=1e-14))
                @test all(isnan,v["corrected"][.!v["observed"]])
            end
            meta=TOML.parsefile(joinpath(a.out,"metadata.toml"))
            @test all_plain(meta)
            @test isequal(meta["support"],MC.plain(case.support_meta))
            @test any(c -> occursin("preprocesses auxiliary Current arrays, then discards them",c),meta["conventions"])
            @test meta["view_records"]==8 && meta["lag_records"]==1190
            @test meta["height_unit"]=="nm"
            @test meta["inputs"]["raw"]["sha256"]==MC.filehash(a.raw)
            @test length(meta["sources"])>7
            @test all(p["sha256"]==MC.filehash(p["path"]) for p in values(meta["sources"]))
            @test all(h==MC.filehash(joinpath(a.out,name)) for (name,h) in meta["artifacts_sha256"])
            @test_throws ErrorException MC.write_comparison(a.out,case,rawdata,views,footprint,d,ctx;inputs=inputs)
            @test_throws ErrorException MC.plain(ctx.options)
            result=ComparisonCLI.masked_comparison_main([a.args[1:end-1];joinpath(dir,"cli_output")])
            @test result.all_common==d.all_common
            for name in expected[1:6]
                @test read(joinpath(a.out,name))==read(joinpath(dir,"cli_output",name))
            end
            # Raw channel scale and unique f/b constraints fail structurally.
            write_synthetic_sxm(a.raw,a.f,a.b;unit="pA")
            @test_throws ErrorException MC.load_real_context(a.raw,a.geometry,a.physical,a.summary,ctx)
            write_synthetic_sxm(a.raw,a.f,a.b;backward=false)
            @test_throws ErrorException MC.load_real_context(a.raw,a.geometry,a.physical,a.summary,ctx)
        end
    end

    @testset "original configs/source locks unchanged" begin
        for (p,h) in ORIGINAL_HASHES
            @test MC.filehash(p)==h
        end
    end
end

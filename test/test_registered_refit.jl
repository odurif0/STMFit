using Test, TOML, Statistics, LinearAlgebra
using STMSXMIO: SXMImage, SXMChannel
include(joinpath(@__DIR__,"run_registered_refit_comparison.jl"))
const RR=RegisteredRefit
const G=RR.G
const PHYS=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
const ASSIGN=TOML.parsefile(joinpath(ROOT,"config","unit_assignment_patch_support.toml"))
const REFIT_SETTINGS=joinpath(ROOT,"config","registered_refit.toml")
const ORIGINAL_SETTINGS=joinpath(ROOT,"config","registered_refit_original_support.toml")
const ACQ_SETTINGS=joinpath(ROOT,"config","acquisition_registration.toml")
BLAS.set_num_threads(1)

function example_image(f,b)
    SXMImage("synthetic.sxm",Dict{String,String}(),size(f,2),size(f,1),(5.,2.),(0.,0.),
        [SXMChannel("Z","nm","fwd",f),SXMChannel("Z","nm","bwd",b)])
end

@testset "Zero-shift execution shares both real fit results exactly" begin
    for setting in (REFIT_SETTINGS,ORIGINAL_SETTINGS)
    mktempdir() do dir
        f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
        f[end,:].=NaN
        img=example_image(f,copy(f))
        # This is a synthetic image and a bounded synthetic optimizer test.
        rawcfg=deepcopy(PHYS)
        rawcfg["model"]["global_maxtime"]=.1; rawcfg["model"]["global_maxiter"]=10; rawcfg["model"]["max_iter"]=10
        rawcfg["preprocessing"]["flatten"]="none"; rawcfg["preprocessing"]["smooth_radius_px"]=0
        cfgpath=joinpath(dir,"count.toml"); open(io->TOML.print(io,rawcfg),cfgpath,"w")
        _,cfg,_=RR.F.Extractor._configs(rawcfg["model"],rawcfg["preprocessing"],"unused")
        pcfg,_,_=RR.F.Extractor._configs(rawcfg["model"],rawcfg["preprocessing"],"unused")
        original=RR.original_support(img,pcfg,cfg)
        axis=original.axis
        p=zeros(G._chain_nparams(3,cfg)); p[1]=.01
        r=G.ChainModelResult(n=3,params=p,amp_min=.03,amp_range=.07,gcv=.1)
        rows=RR.R.feature_rows("a.sxm",r,(axisctx=axis,),cfg,"ell")
        base=joinpath(dir,"base.tsv"); write_table(base,sort(collect(keys(first(rows)))),rows)
        shifts=joinpath(dir,"shifts.tsv"); write_table(shifts,["file","bwd_sample_dx_px"],[Dict("file"=>"a.sxm","bwd_sample_dx_px"=>"0")])
        touch(joinpath(dir,"a.sxm")); out=joinpath(dir,"fit.tsv")
        opts=RR.parse_cli(["--features",base,"--data-dir",dir,"--config",cfgpath,"--assignment-config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",setting,"--shifts",shifts,"--out",out])
        RR.execute(opts;reader=path->img)
        for profile in RR.PROFILES
            @test read(out*".control.$profile.tsv")==read(out*".registered.$profile.tsv")
        end
        @test read(joinpath(out*".fit_data","a.sxm.control.tsv"))==read(joinpath(out*".fit_data","a.sxm.registered.tsv"))
        _,audits=read_table(out*".fits.tsv")
        @test length(audits)==4
        @test all(r["reused_zero_shift"]=="true" && r["elapsed_s"]=="0" for r in audits if r["arm"]=="registered")
        @test all(r["valid"]=="true" for r in audits)
        @test isfile(out*".support.tsv")==RR.settings(setting).original_support
        if RR.settings(setting).original_support
            _,supports=read_table(out*".support.tsv")
            @test length(supports)==1
            @test all(r["support_tmin"]==supports[1]["support_tmin"] && r["support_tmax"]==supports[1]["support_tmax"] for r in audits)
            @test isfile(joinpath(out*".fit_data","a.sxm.original.tsv"))
        end
        @test_throws ErrorException RR.parse_cli(["--features",base,"--data-dir",dir,"--config",cfgpath,"--assignment-config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",REFIT_SETTINGS,"--shifts",shifts,"--out",out])
    end
    end
end

function seed_fixture(;profile="gaussian",circular=false)
    _,cfg,_=RR.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
    cfg.peak_profile=Symbol(profile); cfg.chain_circular_sigmas=circular
    axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
    xs=collect(range(-2.,2.;length=51)); ys=collect(range(-.5,.5;length=15))
    p=zeros(G._chain_nparams(3,cfg)); p[1]=.01
    z=reshape(G._chain_model_values(repeat(xs;inner=length(ys)),repeat(ys;outer=length(xs)),p,3,axis,cfg;
        amp_min=.03,amp_range=.07),length(ys),length(xs))
    return xs,ys,z,axis,cfg,p
end

@testset "Frozen native support stays original while objectives use only observed pixels" begin
    f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
    for dx in (-5,0,4), masked in (false,true)
        b=copy(f); fm=copy(f)
        masked && (fm[end,:].=NaN; fm[20,32]=NaN; b[22,35]=NaN)
        img=example_image(fm,b)
        pcfg=G.PatternConfig(flatten="none",smooth_radius_px=1,no_plot=true)
        _,cfg,_=RR.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
        v=RR.load_views(img,pcfg); original=RR.original_support(img,pcfg,cfg)
        before=deepcopy(original)
        bundle=RR.fused_data(img,pcfg,cfg,v,dx;original)
        @test isequal(original,before)
        @test bundle.data.axisctx===original.axis
        @test bundle.support===original.support
        @test bundle.mask==original.mask .& bundle.observed
        expected=[I for I in original.indices if original.fitmask[I] && bundle.observed[I]]
        @test bundle.indices==expected
        @test all(isfinite,bundle.data.z)
        @test all(bundle.data.z[j]==(v.zf[I[1],I[2]]+v.zb[I[1],I[2]+dx])/2-bundle.offset for (j,I) in enumerate(expected))
        @test length(bundle.data.z)<=count(original.fitmask)
        if dx==0 && !masked
            @test bundle.data.z==original.z[original.keep]
            @test bundle.data.zfull==original.z
        end
        row=Dict("axis_x"=>RR.@sprintf("%.8f",original.axis.axis[1]),"axis_y"=>RR.@sprintf("%.8f",original.axis.axis[2]),
            "origin_x_nm"=>RR.@sprintf("%.6f",original.axis.origin[1]),"origin_y_nm"=>RR.@sprintf("%.6f",original.axis.origin[2]))
        @test RR.check_original_frame(original,row)===nothing
        row["axis_x"]="0.00000000"
        @test_throws ErrorException RR.check_original_frame(original,row)
        @test_throws ErrorException RR.fused_data(img,pcfg,cfg,v,65;original)
    end
    @test RR.settings(ORIGINAL_SETTINGS).original_support
    @test !RR.settings(REFIT_SETTINGS).original_support
    mktempdir() do dir
        s=TOML.parsefile(ORIGINAL_SETTINGS); path=joinpath(dir,"invalid.toml")
        for (section,key,value) in (("preprocessing","geometry","extend_for_N"),("selection","expected_N",6),
                ("selection","count_policy","reselect"),("preprocessing","roi","native_finite_statistics_and_observed_mask"))
            bad=deepcopy(s); bad[section][key]=value
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException RR.settings(path)
        end
    end
    launcher=joinpath(ROOT,"hpc","compare_registered_refit_original_support.sbatch")
    @test success(`bash -n $launcher`)
    @test occursin("config/registered_refit_original_support.toml",read(launcher,String))
end

@testset "Observed opt-in preserves fully observed native ROI and seeds" begin
    f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
    img=example_image(f,copy(f))
    for radius in (0,1,2), stride in (1,2)
        pcfg=G.PatternConfig(flatten="none",stride=stride,smooth_radius_px=radius,no_plot=true)
        _,cfg,_=RR.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
        v=RR.load_views(img,pcfg)
        _,sf=RR.observed_shift(v.zf,v.cf,pcfg,0)
        @test G.molecule_roi_mask_fused(img,pcfg,sf)==G.molecule_roi_mask_fused(img,pcfg,sf;observed_only=true)
        b=RR.fused_data(img,pcfg,cfg,v,0)
        native=G._fused_roi_data(img,pcfg)
        @test b.data.zimg==native[3] && b.mask==native[4]
        @test b.data.noise==native[8]
        @test all(isfinite,b.data.z)
        @test all(b.data.z[j]==(b.f[I]+b.b[I])/2-b.offset for (j,I) in enumerate(b.indices))
    end
    for profile in ("gaussian","split"), circ in (false,true)
        xs,ys,z,axis,cfg,_=seed_fixture(;profile,circular=circ)
        before=copy(z)
        @test G._pack_chain_initial(xs,ys,z,3,axis,cfg)==G._pack_chain_initial(xs,ys,z,3,axis,cfg;observed_only=true)
        @test z==before
        z[1:2,:].=NaN; z[8,20]=NaN
        start,amin,arange=G._pack_chain_initial(xs,ys,z,3,axis,cfg;observed_only=true)
        @test all(isfinite,start) && isfinite(amin) && isfinite(arange)
        @test amin≈cfg.min_amplitude_fraction*maximum(filter(isfinite,vec(z)))
        feature=G.MolecularFeature(1.,xs[20],ys[8],.2,.2,1.)
        sample=only(G._nearest_values_on_grid(xs,ys,z,[feature];observed_only=true))
        indices=findall(isfinite,z)
        nearest=argmin(I->(xs[I[2]]-feature.x_nm)^2+(ys[I[1]]-feature.y_nm)^2,indices)
        @test sample==max(z[nearest],G.EPS)
        @test isnan(z[8,20]) # the fit data were not imputed
        @test_throws ErrorException G._pack_chain_initial(xs,ys,fill(NaN,size(z)),3,axis,cfg;observed_only=true)
    end
end

@testset "Known translation, raw holes and genuinely observed fit samples" begin
    f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
    for dx in (-5,0,4), radius in (0,1)
        b=copy(f)
        for x in 1:65
            b[:,x] .= 1<=x-dx<=65 ? f[:,x-dx] : f[:,1]
        end
        fmasked=copy(f); fmasked[19,32]=NaN; b[17,28]=NaN
        img=example_image(fmasked,b)
        pcfg=G.PatternConfig(flatten="none",stride=1,smooth_radius_px=radius,no_plot=true)
        _,cfg,_=RR.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
        views=RR.load_views(img,pcfg); rawcopy=copy(img.channels[1].data)
        bundle=RR.fused_data(img,pcfg,cfg,views,dx)
        @test isequal(rawcopy,img.channels[1].data)
        @test !bundle.observed[19,32]
        for (j,I) in enumerate(bundle.indices)
            y,x=Tuple(I)
            @test 1<=x+dx<=65 && isfinite(fmasked[y,x]) && isfinite(b[y,x+dx])
            @test bundle.data.z[j]==(views.zf[y,x]+views.zb[y,x+dx])/2-bundle.offset
        end
        @test all(isnan,bundle.data.zimg[:,dx<0 ? (1:-dx) : (66-dx:65)])
    end
    img=example_image(fill(NaN,10,10),fill(NaN,10,10))
    pcfg=G.PatternConfig(flatten="none",smooth_radius_px=0,no_plot=true)
    @test_throws Exception RR.load_views(img,pcfg) |> v->RR.fused_data(img,pcfg,G.ChainSweepConfig(),v,0)
end

@testset "Both profiles fit finite masked data at fixed N with full GCV" begin
    for profile in ("gaussian","split")
        xs,ys,z,axis,cfg,p=seed_fixture(;profile)
        z[1:2,:].=NaN; z[end,:].=NaN
        indices=findall(isfinite,z)
        x=[xs[I[2]] for I in indices]; y=[ys[I[1]] for I in indices]; zz=z[indices]
        data=(xs=xs,ys=ys,zimg=z,x=x,y=y,z=zz,zfull=zz,noise=.01,axisctx=axis)
        raw=deepcopy(PHYS); raw["model"]["global_maxtime"]=.1; raw["model"]["global_maxiter"]=10
        raw["model"]["max_iter"]=10
        fit=RR.fit_profile(3,(data=data,),raw,ASSIGN,profile,RR.settings(REFIT_SETTINGS))
        @test fit.result.n==3 && fit.result.success && fit.result.valid
        @test fit.cfg.peak_profile==Symbol(profile)
        @test fit.result.gcv≈length(zz)/(length(zz)-G._chain_nparams(3,fit.cfg))^2*fit.result.rss
        @test fit.result.gcv==minimum(r.gcv for r in values(fit.boot) if r.success && r.valid)
        @test length(RR.R.feature_rows("a.sxm",fit.result,data,fit.cfg,fit.source))==3
        @test all(isfinite,fit.result.params)
    end
end

@testset "Strict settings, complete matched arms, no labels or overwrite" begin
    @test_throws ErrorException RR.parse_cli(["--expected-N","6"])
    @test_throws ErrorException registered_refit_comparison(["--benchmark-manifest","truth"])
    mktempdir() do dir
        s=TOML.parsefile(REFIT_SETTINGS); path=joinpath(dir,"bad.toml")
        for section in keys(s), key in keys(s[section])
            bad=deepcopy(s); delete!(bad[section],key)
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException RR.settings(path)
        end
        bad=deepcopy(s); bad["selection"]["expected_N"]=6
        open(io->TOML.print(io,bad),path,"w")
        @test_throws ErrorException RR.settings(path)
        raw=joinpath(dir,"raw"); mkdir(raw); touch(joinpath(raw,"a.sxm"))
        xs,ys,z,axis,cfg,p=seed_fixture()
        d=(axisctx=axis,)
        result=G.ChainModelResult(n=3,params=p,amp_min=.03,amp_range=.07,gcv=.1)
        rows=RR.R.feature_rows("a.sxm",result,d,cfg,"ell"); header=sort(collect(keys(first(rows))))
        base=joinpath(dir,"base.tsv"); split=joinpath(dir,"split.tsv"); template=joinpath(dir,"templates.tsv")
        write_table(base,header,rows); cp(base,split); touch(template)
        out=joinpath(dir,"out")
        args=["--data-dir",raw,"--count-config",joinpath(ROOT,"config","chitosan.toml"),
            "--config",joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",ACQ_SETTINGS,"--refit-settings",REFIT_SETTINGS,
            "--features",base,"--split-features",split,"--templates",template,"--outdir",out]
        registered_refit_comparison(vcat(args,["--dry-run"]))
        @test !ispath(out)
        calls=[]
        runner=function(root,name,script,a;threads=4)
            push!(calls,(name,script,a)); op=Dict(a[i]=>a[i+1] for i in 1:2:length(a))
            mkdir(op["--outdir"]); cp(base,joinpath(op["--outdir"],"predictions.tsv"))
        end
        estimator=function(root,opts,rows;runner)
            write_table(joinpath(root,"shifts.tsv"),["file","bwd_sample_dx_px"],[Dict("file"=>"a.sxm","bwd_sample_dx_px"=>"-2")])
        end
        exporter=function(root,name,script,a,path;nfiles)
            @test nfiles==1 && script=="refit_registered_geometry.jl"
            write_table(path,header,rows)
            for arm in RR.ARMS, profile in RR.PROFILES; write_table(path*".$arm.$profile.tsv",header,rows); end
        end
        registered_refit_comparison(args;runner,estimator,exporter)
        @test first.(calls)==["reference","control","registered"]
        @test !("--acquisition-shifts" in calls[1][3])
        @test all("--acquisition-shifts" in c[3] for c in calls[2:end])
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        @test read(joinpath(out,"reference","features.tsv"))==read(base)
        @test_throws ErrorException registered_refit_comparison(args;runner,estimator,exporter)
        failed=replace.(args,out=>joinpath(dir,"failed"))
        @test_throws ErrorException registered_refit_comparison(failed;runner=(a...)->error("intentional"),estimator,exporter)
        @test isfile(joinpath(dir,"failed","failures.tsv"))
        touch(joinpath(raw,"extra.sxm"))
        @test_throws ErrorException registered_refit_comparison(vcat(replace.(args,out=>joinpath(dir,"extra")),["--dry-run"]))
    end
    launcher=joinpath(ROOT,"hpc","compare_registered_refit.sbatch")
    @test success(`bash -n $launcher`)
    shell=read(launcher,String)
    @test occursin("#SBATCH --time=02:00:00",shell) && occursin("SLURM_JOB_ID:?",shell)
    script=read(joinpath(@__DIR__,"refit_registered_geometry.jl"),String)
    @test !occursin("chain_gaussian_sweep(",script) && !occursin("_select_chain_model(",script)
end

using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"verify_image_registered_fit.jl"))
const V=VerifyImageRegisteredFit
const D=V.D
const H=D.H
const G=H.G
const ROOT=dirname(@__DIR__)
config(name)=joinpath(ROOT,"config",name*".toml")
const SETTINGS=D.settings(config("image_registered_fit"))
const AUDIT=H.C.settings(config("paired_convergence"))
const MODEL=H.P.settings(config("paired_acquisition"))
const SOLVER=merge(H.S.settings(config("paired_solver")),(slsqp_maxeval=1000,max_time_s=30.,evaluation_checkpoints=[1,1000]))
BLAS.set_num_threads(1)

function fixture(profile,circular;delta=(.375,-.625))
    cfg=G.ChainSweepConfig(peak_profile=profile,chain_circular_sigmas=circular,chain_tilted_baseline=true)
    n=2; axis=(origin=(0.,0.),axis=[.8,.6],perp=[-.6,.8],tmin=-.8,tmax=.8)
    xs=collect(range(-1.2,1.2;length=33)); ys=collect(range(-.8,.8;length=25))
    x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
    step=(xs[2]-xs[1],ys[2]-ys[1]); p=zeros(G._chain_nparams(n,cfg)); p[1:3]=[.01,.001,-.002]
    q=vcat(p,[.07,.0003,.0001,-.0002,delta...])
    pred=H.predictions(q,n,x,y,axis,cfg,.03,.07,step)
    z=pred.mean; data=(;x,y,z,zfull=z,noise=.01,axisctx=axis,xs,ys,zimg=reshape(z,length(ys),length(xs)))
    initial=G.ChainModelResult(n=n,params=p,success=true,amp_min=.03,amp_range=.07)
    return ((;initial,data,cfg,fwd=pred.fwd,bwd=pred.bwd),step,q)
end

@testset "Image-fixed and free predictors, independent arithmetic, all-parameter GCV" begin
    @test SETTINGS.start_seeds==[0,11,29,47]
    for profile in (:gaussian,:split),circular in (true,false)
        ctx,step,q=fixture(profile,circular); delta=Tuple(q[end-1:end]); limits=(2,3)
        zero=D.problem(ctx,"paired",MODEL,step,limits,delta)
        fixed=D.problem(ctx,"image_fixed",MODEL,step,limits,delta)
        free=D.problem(ctx,"shifted",MODEL,step,limits,delta)
        @test free.lo[end-1:end]==[-2,-3] && free.hi[end-1:end]==[2,3]
        @test fixed.target==free.target==zero.target
        @test free.predict(free.q0)==zero.predict(zero.q0)
        @test fixed.predict(q[1:end-2])==free.predict(q)
        @test fixed.lo==zero.lo && fixed.hi==zero.hi && fixed.q0==zero.q0
        @test H.start(free,11,.02)[end-1:end]==[0.,0.]
        for mode in ("image_fixed","shifted")
            p=mode=="image_fixed" ? fixed : free
            solved=H.S.solve(p,SOLVER,AUDIT); fit=D.finish(solved,p,ctx,mode)
            @test fit.valid && fit.rss<1e-8
            @test fit.np==length(q) && fit.gcv==fit.nd/(fit.nd-fit.np)^2*fit.rss
            independent=V.predict(fit.params,ctx,mode,step,delta)
            @test independent.fwd≈fit.predictions.fwd atol=1e-12
            @test independent.bwd≈fit.predictions.bwd atol=1e-12
            @test fit.result.rss≈sum(abs2,ctx.data.z.-independent.mean) atol=1e-12
            mode=="shifted" && (@test fit.params[end-1:end]≈q[end-1:end] atol=.002)
        end
        @test_throws ErrorException D.problem(ctx,"image_fixed",MODEL,step,limits,(3.,0.))
    end
    mktempdir() do dir
        for (section,key,v) in (("selection","expected_N",6),("model","labels","NKNNKN"),("model","start_seeds",[0,0]),("selection","free_image_agreement_px",0.))
            t=TOML.parsefile(config("image_registered_fit")); t[section][key]=v; path=joinpath(dir,"bad.toml")
            open(io->TOML.print(io,t),path,"w"); @test_throws ErrorException D.settings(path)
        end
    end
    @test_throws ErrorException D.parse_cli(["--labels","forbidden"])
end

@testset "Training-only start choice and independent image disagreement rejection" begin
    rows=Dict{String,String}[]
    for c in D.CASES,m in ("paired","image_fixed","shifted"),f in 0:2,seed in SETTINGS.start_seeds
        push!(rows,Dict("id"=>"$(c.profile).$(c.family).$m.$f.$seed","profile"=>c.profile,"family"=>c.family,
            "mode"=>m,"fold"=>string(f),"seed"=>string(seed),"rss"=>string(1+seed),"valid"=>"true",
            "gcv"=>m=="paired" ? "1.0" : "0.5","heldout_rss"=>m=="paired" ? "1.0" : "0.5",
            "shift_x_px"=>"0.125","shift_y_px"=>"0.0","image_dx_px"=>"0.125","image_dy_px"=>"0.0"))
    end
    summary=D.summaries(rows,SETTINGS,(11,11))
    @test length(summary.winners)==36 && all(r["seed"]=="0" for r in summary.winners)
    @test all(r["eligible"]=="true" for r in summary.pilots)
    changed=deepcopy(rows)
    for r in changed; r["mode"]=="shifted" && (r["shift_x_px"]="-1.0"); end
    other=D.summaries(changed,SETTINGS,(11,11))
    @test getindex.(other.winners,"id")==getindex.(summary.winners,"id")
    @test all(r["eligible"]==(r["mode"]=="image_fixed" ? "true" : "false") for r in other.pilots)
    @test_throws ErrorException D.summaries(rows[2:end],SETTINGS,(11,11))
end

@testset "Mocked saved-image binding, all shards and independent fit-output verification" begin
    mktempdir() do dir
        native=joinpath(dir,"native"); reference=joinpath(dir,"reference"); image=joinpath(dir,"image")
        mkpath(native); mkpath(reference); mkpath(joinpath(image,"fold0"))
        s=TOML.parsefile(config("image_registered_fit")); s["model"]["start_seeds"]=[0,11]; s["model"]["slsqp_maxeval"]=3
        fitsettings=joinpath(dir,"fit.toml"); open(io->TOML.print(io,s),fitsettings,"w")
        image_settings=TOML.parsefile(config("image_registration_2d"))
        image_settings["preprocessing"]["holdout_block_px"]=8; image_settings["preprocessing"]["holdout_buffer_px"]=0
        open(io->TOML.print(io,image_settings),joinpath(image,"settings.toml"),"w")
        cp(config("chitosan"),joinpath(image,"physical.toml"))
        raw=TOML.parsefile(config("chitosan")); file=s["selection"]["file"]; n=2
        _,ec,_=D.D.D.F.Extractor._configs(raw["model"],raw["preprocessing"],dir)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-.8,tmax=.8)
        xs=collect(range(-1.2,1.2;length=33)); ys=collect(range(-.8,.8;length=25)); ny,nx=length(ys),length(xs)
        x=repeat(xs;inner=ny); y=repeat(ys;outer=nx); step=(xs[2]-xs[1],ys[2]-ys[1])
        p=zeros(G._chain_nparams(n,ec)); p[1]=.01
        z=G._chain_model_values(x,y,p,n,axis,ec;amp_min=.03,amp_range=.07).+0.0001 .* sin.(eachindex(x))
        pixels=[Dict("row"=>string(mod1(j,ny)),"column"=>string(cld(j,ny)),"x_nm"=>string(x[j]),"y_nm"=>string(y[j]),
            "z_nm"=>string(z[j]),"fwd_nm"=>string(z[j]+.001),"bwd_nm"=>string(z[j]-.001)) for j in eachindex(x)]
        D.write_records(joinpath(native,"registered.pixels.tsv"),pixels)
        fits=Dict{String,String}[]; pars=Dict{String,String}[]; previous=Dict{String,String}[]
        for c in D.CASES
            cfg=deepcopy(ec); cfg.peak_profile=Symbol(c.profile); cfg.chain_circular_sigmas=c.family=="circ"
            pp=zeros(G._chain_nparams(n,cfg)); pp[1]=.01
            for arm in ("control","registered")
                push!(fits,Dict("file"=>file,"N"=>string(n),"arm"=>arm,"profile"=>c.profile,"family"=>c.family,"success"=>"true","dx_px"=>"0",
                    "axis_x"=>"1.0","axis_y"=>"0.0","origin_x_nm"=>"0.0","origin_y_nm"=>"0.0","support_tmin"=>"-.8","support_tmax"=>".8",
                    "amp_min"=>".03","amp_range"=>".07","noise"=>".01","offset"=>"0.0","n_pixels"=>string(length(z))))
                for j in eachindex(pp)
                    push!(pars,Dict("arm"=>arm,"profile"=>c.profile,"family"=>c.family,"stage"=>"final","start"=>"1","parameter"=>string(j),"value"=>string(pp[j])))
                end
            end
            for mode in ("fused","paired"); push!(previous,Dict("file"=>file,"profile"=>c.profile,"family"=>c.family,"mode"=>mode,"rss"=>".1")); end
        end
        D.write_records(joinpath(native,"fits.tsv"),fits); D.write_records(joinpath(native,"parameters.tsv"),pars)
        basehash=bytes2hex(sha256("mocked integer shift input"))
        hashes=[Dict("input"=>key,"sha256"=>bytes2hex(sha256(read(path)))) for (key,path) in (("--config",config("chitosan")),("--assignment-config",config("unit_assignment_patch_support")))]
        push!(hashes,Dict("input"=>"--shifts","sha256"=>basehash)); D.write_records(joinpath(native,"input_hashes.tsv"),hashes)
        D.write_records(joinpath(reference,"fits.tsv"),previous); D.write_records(joinpath(reference,"parameters.tsv"),pars)
        cp(config("paired_acquisition"),joinpath(reference,"paired_settings.toml"))
        D.write_records(joinpath(reference,"native_hashes.tsv"),[Dict("input"=>f,"sha256"=>bytes2hex(sha256(read(joinpath(native,f))))) for f in ("fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv")])
        # Calibration is deliberately mocked here: real image estimation/masks are
        # tested separately. This fixture checks the second runner's data binding.
        ir=[Dict("file"=>file,"fold"=>string(f),"accepted"=>"true","dx"=>"0.0","dy"=>"0.0","base_dx"=>"0","base_dy"=>"0",
            "step_x_nm"=>string(step[1]),"step_y_nm"=>string(step[2]),"limit_x_px"=>string(floor(Int,.16/step[1])),"limit_y_px"=>string(floor(Int,.16/step[2])),
            "nrows"=>string(ny),"ncols"=>string(nx)) for f in 0:2]
        D.write_records(joinpath(image,"summary.tsv"),ir); D.write_records(joinpath(image,"eligibility.tsv"),[Dict("eligible"=>"true","spread_px"=>"0.0")])
        D.write_records(joinpath(image,"input_hashes.tsv"),[Dict("input"=>"--base-shifts","sha256"=>basehash)])
        open(joinpath(image,"fold0","images.tsv"),"w") do io
            println(io,"row\tcolumn\tx_nm\ty_nm\tfwd_nm\tbwd_raw_nm")
            for j in eachindex(x); println(io,join((mod1(j,ny),cld(j,ny),x[j],y[j],z[j]+.001,z[j]-.001),'\t')); end
        end
        common=["--native-dir",native,"--reference-dir",reference,"--image-dir",image,"--physical-config",config("chitosan"),
            "--assignment-config",config("unit_assignment_patch_support"),"--model-settings",config("paired_acquisition"),
            "--settings",config("paired_convergence"),"--solver-settings",config("paired_solver"),"--fit-settings",fitsettings]
        dry=joinpath(dir,"dry"); D.execute(D.parse_cli(vcat(common,["--outdir",dry,"--dry-run"]))); @test !ispath(dry)
        out=joinpath(dir,"run"); mkpath(out)
        for i in 1:4
            args=vcat(common,["--chunk","$i/4","--outdir",joinpath(out,"chunk$i")]); D.execute(D.parse_cli(args))
            @test_throws ErrorException D.parse_cli(args)
        end
        D.merge_chunks(out); @test length(D.table(joinpath(out,"fits.tsv")))==72
        V.verify(native,reference,image,out); @test_throws Exception D.merge_chunks(out)
        # A nonidentified image is a hard stop, never permission for a free fit.
        ir[1]["accepted"]="false"
        open(joinpath(image,"summary.tsv"),"w") do io
            keys_=sort(collect(keys(first(ir)))); println(io,join(keys_,'\t'))
            for row in ir; println(io,join(getindex.(Ref(row),keys_),'\t')); end
        end
        opts=D.parse_cli(vcat(common,["--outdir",joinpath(dir,"rejected"),"--dry-run"]))
        @test_throws ErrorException D.inputs(opts)
    end
    @test success(`bash -n $(joinpath(ROOT,"hpc","diagnose_image_registered_fit.sbatch"))`)
end

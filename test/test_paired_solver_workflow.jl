using Test, TOML, SHA
include(joinpath(@__DIR__,"verify_paired_solver.jl"))
const V=VerifyPairedSolver
const D=V.D
const G=D.S.G
const ROOT=dirname(@__DIR__)

@testset "Synthetic saved-pixel solver workflow, merge and independent re-reading" begin
    mktempdir() do dir
        native=joinpath(dir,"native"); reference=joinpath(dir,"reference"); comparison=joinpath(dir,"comparison")
        mkpath(native); mkpath(reference); mkpath(comparison)
        physical=joinpath(ROOT,"config","chitosan.toml"); assignment=joinpath(ROOT,"config","unit_assignment_patch_support.toml")
        model=TOML.parsefile(joinpath(ROOT,"config","paired_acquisition.toml")); model["model"]["local_maxiter"]=3
        audit=TOML.parsefile(joinpath(ROOT,"config","paired_convergence.toml")); audit["model"]["control_maxiter"]=3
        audit["model"]["extended_maxiter"]=6; audit["model"]["checkpoints"]=[0,3,6]
        solver=TOML.parsefile(joinpath(ROOT,"config","paired_solver.toml")); solver["model"]["lm_maxiter"]=6
        solver["model"]["slsqp_maxeval"]=6; solver["model"]["evaluation_checkpoints"]=[1,3,6]
        modelpath=joinpath(dir,"model.toml"); auditpath=joinpath(dir,"audit.toml"); solverpath=joinpath(dir,"solver.toml")
        for (path,t) in ((modelpath,model),(auditpath,audit),(solverpath,solver)); open(io->TOML.print(io,t),path,"w"); end
        cp(modelpath,joinpath(reference,"paired_settings.toml")); cp(modelpath,joinpath(comparison,"model_settings.toml")); cp(auditpath,joinpath(comparison,"settings.toml"))
        file=solver["selection"]["file"]; n=3; raw=TOML.parsefile(physical)
        _,ec,_=D.D.F.Extractor._configs(raw["model"],raw["preprocessing"],dir)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=21)); ys=collect(range(-.5,.5;length=7)); x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(n,ec)); p[1]=.01
        z=G._chain_model_values(x,y,p,n,axis,ec;amp_min=.03,amp_range=.07).+0.0001 .* sin.(eachindex(x))
        pixels=[Dict("row"=>string(mod1(j,length(ys))),"column"=>string(cld(j,length(ys))),"x_nm"=>string(x[j]),"y_nm"=>string(y[j]),
            "z_nm"=>string(z[j]),"fwd_nm"=>string(z[j]+.001),"bwd_nm"=>string(z[j]-.001)) for j in eachindex(x)]
        D.write_records(joinpath(native,"registered.pixels.tsv"),pixels)
        fits=Dict{String,String}[]; pars=Dict{String,String}[]; prev=Dict{String,String}[]; prevpars=Dict{String,String}[]; oldfits=Dict{String,String}[]
        for profile in ("gaussian","split"),family in ("circ","ell")
            cfg=deepcopy(ec); cfg.peak_profile=Symbol(profile); cfg.chain_circular_sigmas=family=="circ"
            pp=zeros(G._chain_nparams(n,cfg)); pp[1]=.01
            for arm in ("control","registered")
                push!(fits,Dict("file"=>file,"N"=>string(n),"arm"=>arm,"profile"=>profile,"family"=>family,"success"=>"true",
                    "axis_x"=>"1.0","axis_y"=>"0.0","origin_x_nm"=>"0.0","origin_y_nm"=>"0.0","support_tmin"=>"-1.5","support_tmax"=>"1.5",
                    "amp_min"=>"0.03","amp_range"=>"0.07","noise"=>"0.01","offset"=>"0.0","n_pixels"=>string(length(z))))
                for j in eachindex(pp)
                    push!(pars,Dict("arm"=>arm,"profile"=>profile,"family"=>family,"stage"=>"final","start"=>"1","parameter"=>string(j),"value"=>string(pp[j])))
                end
            end
            for mode in ("fused","paired")
                base=Dict("file"=>file,"profile"=>profile,"family"=>family,"mode"=>mode,"rss"=>"0.1")
                push!(prev,base)
                for stage in ("control","extended"); push!(oldfits,merge(base,Dict("stage"=>stage))); end
                for j in 1:(length(pp)+(mode=="paired" ? 4 : 0))
                    push!(prevpars,Dict("profile"=>profile,"family"=>family,"mode"=>mode,"parameter"=>string(j),"value"=>string(j<=length(pp) ? pp[j] : 0.)))
                end
            end
        end
        D.write_records(joinpath(native,"fits.tsv"),fits); D.write_records(joinpath(native,"parameters.tsv"),pars)
        D.write_records(joinpath(native,"input_hashes.tsv"),[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(path)))) for (k,path) in (("--config",physical),("--assignment-config",assignment))])
        D.write_records(joinpath(reference,"fits.tsv"),prev); D.write_records(joinpath(reference,"parameters.tsv"),prevpars)
        D.write_records(joinpath(reference,"native_hashes.tsv"),[Dict("input"=>f,"sha256"=>bytes2hex(sha256(read(joinpath(native,f))))) for f in ("fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv")])
        D.write_records(joinpath(comparison,"fits.tsv"),oldfits)
        for (i,case) in enumerate(D.D.CASES)
            shard=joinpath(comparison,"chunk$(mod1(i,4))"); mkpath(shard)
            rows=filter(r->r["profile"]==case.profile && r["family"]==case.family && r["mode"]==case.mode,prevpars)
            D.write_records(joinpath(shard,"$(case.profile).$(case.family).$(case.mode).parameters.tsv"),
                [merge(r,Dict("stage"=>stage)) for r in rows for stage in ("control","extended")])
        end
        hashes=[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(path)))) for (k,path) in (("--physical-config",physical),("--assignment-config",assignment),("--model-settings",modelpath),("--settings",auditpath))]
        for (tag,path,files) in (("native",native,["fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv"]),
                ("reference",reference,["fits.tsv","parameters.tsv","native_hashes.tsv","paired_settings.toml"]))
            for f in files; push!(hashes,Dict("input"=>"$tag/$f","sha256"=>bytes2hex(sha256(read(joinpath(path,f)))))); end
        end
        D.write_records(joinpath(comparison,"input_hashes.tsv"),hashes)
        common=["--native-dir",native,"--reference-dir",reference,"--comparison-dir",comparison,"--physical-config",physical,
            "--assignment-config",assignment,"--model-settings",modelpath,"--settings",auditpath,"--solver-settings",solverpath]
        dry=vcat(common,["--outdir",joinpath(dir,"dry"),"--dry-run"]); D.execute(D.parse_cli(dry))
        @test !ispath(joinpath(dir,"dry"))
        out=joinpath(dir,"run"); mkpath(out)
        for i in 1:4
            args=vcat(common,["--outdir",joinpath(out,"chunk$i"),"--chunk","$i/4"])
            D.execute(D.parse_cli(args))
            @test_throws ErrorException D.parse_cli(args)
        end
        D.merge_chunks(out)
        @test length(D.table(joinpath(out,"fits.tsv")))==16
        @test_throws Exception D.merge_chunks(out)
        V.verify(native,reference,comparison,out)
        launcher=joinpath(ROOT,"hpc","diagnose_paired_solver.sbatch")
        @test success(`bash -n $launcher`)
        @test occursin("SLURM_JOB_ID:?",read(launcher,String))
    end
end

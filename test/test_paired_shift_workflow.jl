using Test, TOML, SHA
include(joinpath(@__DIR__,"verify_paired_shift.jl"))
const V=VerifyPairedShift
const D=V.D
const G=D.H.G
const ROOT=dirname(@__DIR__)

@testset "Synthetic saved-pixel pilot, all shards, merge and independent verification" begin
    mktempdir() do dir
        native=joinpath(dir,"native"); reference=joinpath(dir,"reference"); mkpath(native); mkpath(reference)
        config(k)=joinpath(ROOT,"config",k*".toml")
        physical=config("chitosan"); assignment=config("unit_assignment_patch_support")
        model=config("paired_acquisition"); audit=config("paired_convergence"); solver=config("paired_solver")
        settings=TOML.parsefile(config("paired_shift")); settings["model"]["start_seeds"]=[0,11]
        settings["model"]["slsqp_maxeval"]=3; settings["preprocessing"]["holdout_block_px"]=4
        settings["preprocessing"]["holdout_buffer_px"]=0
        shift=joinpath(dir,"shift.toml"); open(io->TOML.print(io,settings),shift,"w")
        raw=TOML.parsefile(physical); file=settings["selection"]["file"]; n=2
        _,ec,_=D.D.F.Extractor._configs(raw["model"],raw["preprocessing"],dir)
        axis=(origin=(0.,0.),axis=[1.,0.],perp=[0.,1.],tmin=-1.5,tmax=1.5)
        xs=collect(range(-2.,2.;length=21)); ys=collect(range(-.5,.5;length=13))
        x=repeat(xs;inner=length(ys)); y=repeat(ys;outer=length(xs))
        p=zeros(G._chain_nparams(n,ec)); p[1]=.01
        z=G._chain_model_values(x,y,p,n,axis,ec;amp_min=.03,amp_range=.07).+0.0001 .* sin.(eachindex(x))
        pixels=[Dict("row"=>string(mod1(j,length(ys))),"column"=>string(cld(j,length(ys))),"x_nm"=>string(x[j]),"y_nm"=>string(y[j]),
            "z_nm"=>string(z[j]),"fwd_nm"=>string(z[j]+.001),"bwd_nm"=>string(z[j]-.001)) for j in eachindex(x)]
        D.write_records(joinpath(native,"registered.pixels.tsv"),pixels)
        fits=Dict{String,String}[]; pars=Dict{String,String}[]; previous=Dict{String,String}[]
        for case in D.CASES
            cfg=deepcopy(ec); cfg.peak_profile=Symbol(case.profile); cfg.chain_circular_sigmas=case.family=="circ"
            pp=zeros(G._chain_nparams(n,cfg)); pp[1]=.01
            for arm in ("control","registered")
                push!(fits,Dict("file"=>file,"N"=>string(n),"arm"=>arm,"profile"=>case.profile,"family"=>case.family,"success"=>"true",
                    "axis_x"=>"1.0","axis_y"=>"0.0","origin_x_nm"=>"0.0","origin_y_nm"=>"0.0","support_tmin"=>"-1.5","support_tmax"=>"1.5",
                    "amp_min"=>"0.03","amp_range"=>"0.07","noise"=>"0.01","offset"=>"0.0","n_pixels"=>string(length(z))))
                for j in eachindex(pp)
                    push!(pars,Dict("arm"=>arm,"profile"=>case.profile,"family"=>case.family,"stage"=>"final","start"=>"1","parameter"=>string(j),"value"=>string(pp[j])))
                end
            end
            for mode in ("fused","paired")
                push!(previous,Dict("file"=>file,"profile"=>case.profile,"family"=>case.family,"mode"=>mode,"rss"=>"0.1"))
            end
        end
        D.write_records(joinpath(native,"fits.tsv"),fits); D.write_records(joinpath(native,"parameters.tsv"),pars)
        D.write_records(joinpath(native,"input_hashes.tsv"),[Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(path)))) for (k,path) in (("--config",physical),("--assignment-config",assignment))])
        D.write_records(joinpath(reference,"fits.tsv"),previous); D.write_records(joinpath(reference,"parameters.tsv"),pars)
        cp(model,joinpath(reference,"paired_settings.toml"))
        D.write_records(joinpath(reference,"native_hashes.tsv"),[Dict("input"=>f,"sha256"=>bytes2hex(sha256(read(joinpath(native,f))))) for f in ("fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv")])
        common=["--native-dir",native,"--reference-dir",reference,"--physical-config",physical,"--assignment-config",assignment,
            "--model-settings",model,"--settings",audit,"--solver-settings",solver,"--shift-settings",shift]
        dry=joinpath(dir,"dry"); D.execute(D.parse_cli(vcat(common,["--outdir",dry,"--dry-run"])))
        @test !ispath(dry)
        out=joinpath(dir,"run"); mkpath(out)
        for i in 1:4
            args=vcat(common,["--chunk","$i/4","--outdir",joinpath(out,"chunk$i")])
            D.execute(D.parse_cli(args))
            @test_throws ErrorException D.parse_cli(args)
        end
        D.merge_chunks(out)
        @test length(D.table(joinpath(out,"fits.tsv")))==56
        V.verify(native,reference,out)
        @test_throws Exception D.merge_chunks(out)
        # A missing result is a hard incomplete-cohort error, never a successful partial pilot.
        rows=D.table(joinpath(out,"fits.tsv"))
        @test_throws ErrorException D.summaries(rows[2:end],D.H.settings(shift))
    end
end

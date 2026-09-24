using Test, TOML, LinearAlgebra
include(joinpath(@__DIR__,"test_patch_preprocessing.jl"))
include(joinpath(@__DIR__,"run_observed_fit_comparison.jl"))
include(joinpath(@__DIR__,"verify_observed_fit.jl"))
const C=ObservedFitComparison
const O=C.ObservedGeometry
const G=O.G
const P=C.P
const PHYS=TOML.parsefile(COUNT_CONFIG)
const OBS_SETTINGS=joinpath(ROOT,"config","observed_fit.toml")
const ASSIGN=joinpath(ROOT,"config","unit_assignment_patch_support.toml")
BLAS.set_num_threads(1)

function image_fixture(;missing=false,swap=false)
    f=[.01+.12sum(exp(-((x-c)/5)^2-((y-21)/6)^2) for c in (22,32,42))+.0001sin(x+3y) for y in 1:41,x in 1:65]
    b=copy(f)
    missing && (f[end-2:end,:].=NaN; b[3,4]=Inf)
    swap && ((f,b)=(b,f))
    G.SXMImage("synthetic.sxm",Dict{String,String}(),size(f,2),size(f,1),(5.,2.),(0.,0.),
        [G.SXMChannel("Z","nm","fwd",f),G.SXMChannel("Z","nm","bwd",b)])
end

@testset "Observed preprocessing: unchanged finite path and explicit mask policy" begin
    s=O.settings(OBS_SETTINGS)
    for flatten in ("none","plane","rows","plane+rows"), radius in (0,1,2)
        img=image_fixture(); pcfg=G.PatternConfig(channel="Z",direction="fwd",stride=1,flatten=flatten,smooth_radius_px=radius)
        a=G.preprocess_channel(img,first(img.channels),pcfg)
        b=G.preprocess_observed_fit_channel(img,first(img.channels),pcfg;plane_rank_rtol=s.plane_rank_rtol)
        @test isequal(a,b)
        partial=image_fixture(missing=true)
        c=G.preprocess_observed_fit_channel(partial,first(partial.channels),pcfg;plane_rank_rtol=s.plane_rank_rtol)
        @test isfinite.(c[3])==isfinite.(c[4])==isfinite.(first(partial.channels).data)
        @test all(isfinite(c[4][I]) for I in findall(isfinite.(c[5])))
        @test count(isfinite,c[5])<=count(isfinite,c[4])
        vals=c[5][isfinite.(c[5])]
        @test c[7]==max(1.4826median(abs.(vals.-median(vals))),.1std(vals),G.EPS)
        cfg=deepcopy(pcfg); cfg.direction="bwd"
        @test_throws ErrorException G.preprocess_observed_fit_channel(img,first(img.channels),cfg;plane_rank_rtol=s.plane_rank_rtol)
    end
    mktempdir() do dir
        path=joinpath(dir,"config.toml")
        observed=C.candidate_config(PHYS,s)
        open(io->TOML.print(io,observed),path,"w")
        p=load_patch_preprocessing(path)
        @test p.missing_pixel_policy=="observed_only" && p.plane_rank_rtol==1e-12
        img=image_fixture(missing=true); cfg=G.PatternConfig(channel="Z",direction="fwd",stride=p.stride,flatten=p.flatten,smooth_radius_px=p.smooth_radius_px)
        @test isequal(preprocess_patch_channel(img,first(img.channels),cfg,p),
            G.preprocess_observed_fit_channel(img,first(img.channels),cfg;plane_rank_rtol=1e-12))
        delete!(observed["preprocessing"],"missing_pixel_policy"); delete!(observed["preprocessing"],"plane_rank_rtol")
        @test observed==PHYS
        for rank in (nothing,0.,-1.,NaN,1.,true)
            bad=C.candidate_config(PHYS,s)
            rank===nothing ? delete!(bad["preprocessing"],"plane_rank_rtol") : (bad["preprocessing"]["plane_rank_rtol"]=rank)
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ArgumentError load_patch_preprocessing(path)
        end
        bad=C.candidate_config(PHYS,s); bad["preprocessing"]["missing_pixel_policy"]="fill_then_mask"
        open(io->TOML.print(io,bad),path,"w")
        @test_throws ArgumentError load_patch_preprocessing(path)
    end
end

@testset "Fit inputs preserve finite control and forward/backward symmetry" begin
    pcfg,ell,_=O.F.Extractor._configs(PHYS["model"],PHYS["preprocessing"],"unused")
    s=O.settings(OBS_SETTINGS)
    for missing in (false,true)
        img=image_fixture(;missing); v=O.observed_views(img,pcfg,s)
        a=O.bundle(img,pcfg,ell,"control",v); b=O.bundle(img,pcfg,ell,"observed",v)
        if !missing
            @test O.same_data(a,b)
        else
            @test !O.same_data(a,b)
            @test all(v.observed[b.indices]) && all(isfinite,b.data.z)
        end
        swapped=image_fixture(;missing,swap=true); vs=O.observed_views(swapped,pcfg,s)
        bs=O.bundle(swapped,pcfg,ell,"observed",vs)
        @test isequal(v,vs) && O.same_data(b,bs)
    end
    img=image_fixture(); img.channels[1].data[:,:].=NaN
    @test_throws ErrorException O.observed_views(img,pcfg,s)
end

@testset "Both patch CLIs use the candidate background and preserve holes" begin
    for missing in (false,true)
        mktempdir() do dir
            sxm=synthetic_sxm(joinpath(dir,"synthetic.sxm");missing)
            features=joinpath(dir,"features.tsv")
            write(features,"file\tlobe\tt_nm\tu_nm\tamplitude\taxis_x\taxis_y\tx_nm\ty_nm\tsigma_parallel_nm\tsigma_perp_nm\tbaseline\ttilt_x\ttilt_y\tskew_ratio\n" *
                "synthetic.sxm\t1\t0\t0\t2\t1\t0\t1.2\t1.2\t0.2\t0.15\t0.4\t0.17\t0.08\t1\n")
            config=joinpath(dir,"observed.toml")
            open(io->TOML.print(io,C.candidate_config(PHYS,O.settings(OBS_SETTINGS))),config,"w")
            for (name,mod) in (("fwd",ForwardPatches),("bwd",BackwardPatches))
                output=joinpath(dir,name*".tsv")
                args=["--features",features,"--data-dir",dir,"--out",output,"--assignment-config",ASSIGN]
                extract_quietly(mod,vcat(args,["--config",COUNT_CONFIG])); legacy=read(output)
                extract_quietly(mod,vcat(args,["--config",config])); candidate=read(output)
                @test (candidate==legacy)==!missing
                _,rows=mod.ScriptUtils._read_tsv(output)
                @test length(rows)==1
                if missing
                    @test any(v=="NA" for (k,v) in only(rows) if occursin("res_p",k))
                end
            end
        end
    end
end

function fixture(dir)
    raw=joinpath(dir,"raw"); mkdir(raw); touch(joinpath(raw,"a.sxm"))
    summary=joinpath(dir,"summary.tsv")
    P.write_table(summary,["filepath","N_selected"],[Dict("filepath"=>"a.sxm","N_selected"=>"3")])
    t=joinpath(dir,"templates.tsv"); write(t,"fixture")
    cfg=deepcopy(PHYS)
    cfg["model"]["global_maxtime"]=.1; cfg["model"]["global_maxiter"]=10; cfg["model"]["max_iter"]=10
    cfg["preprocessing"]["flatten"]="none"; cfg["preprocessing"]["smooth_radius_px"]=0
    path=joinpath(dir,"count.toml"); open(io->TOML.print(io,cfg),path,"w")
    Dict("--data-dir"=>raw,"--count-config"=>path,"--config"=>ASSIGN,"--templates"=>t,
        "--selected-summary"=>summary,"--settings"=>OBS_SETTINGS,"--outdir"=>joinpath(dir,"out"))
end
refit_options(o,out)=Dict("--data-dir"=>o["--data-dir"],"--selected-summary"=>o["--selected-summary"],
    "--config"=>o["--count-config"],"--assignment-config"=>o["--config"],"--settings"=>o["--settings"],"--out"=>out)

@testset "Bounded synthetic fit: exact reuse, common-pixel RSS, no partial outputs" begin
    mktempdir() do dir
        o=fixture(dir); out=joinpath(dir,"fit.tsv"); ro=refit_options(o,out)
        O.execute(ro;reader=path->image_fixture())
        for profile in O.PROFILES
            @test read(out*".control.$profile.tsv")==read(out*".observed.$profile.tsv")
        end
        _,fits=P.read_table(out*".fits.tsv")
        @test length(fits)==4 && all(r["valid"]=="true" for r in fits)
        @test all(r["reused_control"]=="true" && r["elapsed_s"]=="0" for r in fits if r["arm"]=="observed")
        _,common=P.read_table(out*".common.tsv")
        @test length(common)==2 && all(r["control_rss"]==r["observed_rss"] for r in common)
        for r in common
            _,data=P.read_table(joinpath(out*".fit_data","a.sxm.$(r["profile"]).common.tsv"))
            rss=sum((parse(Float64,d["control_target_nm"])-parse(Float64,d["control_model_nm"]))^2 for d in data)
            @test rss≈parse(Float64,r["control_rss"]) rtol=1e-10
        end
        _,ds=P.read_table(out*".diagnostics.tsv")
        @test !isempty(ds) && all(r["arm"]=="control" for r in ds)
        @test all(r["lm_status"] in ("converged","not_converged","exception") for r in ds)
        failed=joinpath(dir,"failed.tsv")
        @test_throws ErrorException O.execute(refit_options(o,failed);reader=path->image_fixture(),fitter=(args...;kwargs...)->error("synthetic failure"))
        _,fs=P.read_table(failed*".failures.tsv")
        @test length(fs)==4 && !isfile(failed)
        @test all(!isfile(failed*".$a.$p.tsv") for a in O.ARMS for p in O.PROFILES)
        @test_throws ErrorException O.parse_cli(vcat([[k,v] for (k,v) in ro]...))
    end
end

@testset "Strict CLI and complete cohort before downstream assignment" begin
    @test C.options(["--help"])===nothing && O.parse_cli(["--help"])===nothing
    for key in ("--truth","--expected-N","--manifest","--features","--selection-policy","--sequence")
        @test_throws ErrorException C.options([key,"bad"])
        @test_throws ErrorException O.parse_cli([key,"bad"])
    end
    mktempdir() do dir
        o=fixture(dir)
        @test C.options(vcat([[k,v] for (k,v) in o]...))==o
        C.execute(merge(o,Dict("--dry-run"=>"true")))
        @test !ispath(o["--outdir"])
        calls=[]
        function exporter(out,stage,script,args,path;nfiles)
            @test script=="refit_observed_geometry.jl" && nfiles==1
            ro=refit_options(o,path)
            O.execute(ro;reader=p->image_fixture())
        end
        function runner(out,stage,script,args)
            push!(calls,(stage,args)); value(k)=args[findfirst(==(k),args)+1]
            @test value("--selected-summary")==o["--selected-summary"]
            @test value("--config")==o["--config"] && value("--templates")==o["--templates"]
            cfg=TOML.parsefile(value("--count-config"))
            base=TOML.parsefile(o["--count-config"])
            @test cfg==(stage=="control" ? base : C.candidate_config(base,O.settings(OBS_SETTINGS)))
            mkdir(value("--outdir"))
            for name in ("features_local","features_descriptor","features_predictor","patches_fwd17","patches_bwd17","patches_bwd9",
                    "training_support","score_fwd","score_bwd","fisher_cv","pred_gmm","pred_kmeans","predictions")
                cp(value("--features"),joinpath(value("--outdir"),name*".tsv"))
            end
        end
        C.execute(o;runner,exporter)
        @test first.(calls)==["control","observed"]
        v=ObservedFitVerification.verify(o["--outdir"];rawdir=o["--data-dir"])
        @test isempty(v.failures) && occursin("Both complete",ObservedFitVerification.report(v))
        @test_throws ErrorException C.options(vcat([[k,v] for (k,v) in o]...))
        empty!(calls)
        bad=merge(o,Dict("--outdir"=>joinpath(dir,"bad")))
        @test_throws ErrorException C.execute(bad;runner,exporter=(args...;kwargs...)->error("missing fits"))
        @test isempty(calls) && isfile(joinpath(bad["--outdir"],"failures.tsv"))
        partial=merge(o,Dict("--outdir"=>joinpath(dir,"partial")))
        function failed_exporter(out,stage,script,args,path;nfiles)
            O.execute(refit_options(o,path);reader=p->image_fixture(),
                fitter=(args...;kwargs...)->args[5]=="gaussian" ? error("synthetic Gaussian failure") : O.fit_profile(args...;kwargs...))
        end
        @test_throws ErrorException C.execute(partial;runner,exporter=failed_exporter)
        v=ObservedFitVerification.verify(partial["--outdir"];rawdir=o["--data-dir"])
        @test length(v.fits)==2 && length(v.failures)==2
        @test occursin("Incomplete cohort",ObservedFitVerification.report(v))
        forbidden=joinpath(dir,"forbidden.tsv")
        P.write_table(forbidden,["filepath","N_selected","expected_N"],
            [Dict("filepath"=>"a.sxm","N_selected"=>"3","expected_N"=>"3")])
        @test_throws ErrorException C.execute(merge(o,Dict("--outdir"=>joinpath(dir,"forbidden"),"--selected-summary"=>forbidden)))
    end
    sbatch=joinpath(ROOT,"hpc","compare_observed_fit.sbatch")
    @test success(`bash -n $sbatch`)
    source=read(sbatch,String)
    for key in ("--time=02:00:00","--cpus-per-task=4","--mem=16000MB","SLURM_JOB_ID","STMFIT_SELECTED_SUMMARY")
        @test occursin(key,source)
    end
end

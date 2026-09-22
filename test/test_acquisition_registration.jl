using Test, Random, Statistics, TOML
include(joinpath(@__DIR__,"lib","acquisition_registration.jl"))
include(joinpath(@__DIR__,"lib","patch_acquisition.jl"))
using .AcquisitionRegistration, .PatchAcquisition
using STMSXMIO: SXMChannel, _box_smooth
const SETTINGS=joinpath(@__DIR__,"..","config","acquisition_registration.toml")
const S=load_settings(SETTINGS)

function shifted_signal(a,lag;gain=1.,offset=0.)
    b=fill(NaN,size(a)); nx=size(a,2)
    for x in 1:nx
        1<=x+lag<=nx && (b[:,x+lag] .= gain.*a[:,x].+offset)
    end
    return b
end

@testset "Fixed support, signed direction, stability, ambiguity and rejection" begin
    rng=MersenneTwister(123)
    a=randn(rng,64,128)
    for lag in (-13,0,11)
        b=shifted_signal(a,lag;gain=2.3,offset=7.)
        a0,b0=copy(a),copy(b)
        r=scan(a,b,.04,S)
        @test r.accepted && r.applied_dx_px==lag
        @test r.peaks[1].correlation ≈ 1.
        @test all(p.best_dx_px==lag && p.usable for p in r.peaks)
        @test isequal(a,a0) && isequal(b,b0)
        @test all(isfinite(a[I]) && all(isfinite(b[I[1],I[2]+dx]) for dx in r.lags) for I in findall(r.support))
        again=scan(a,b,.04,S)
        @test isequal(r,again)
        swapped=scan(b,a,.04,S)
        @test swapped.accepted && swapped.applied_dx_px == -lag
        masked=copy(b); masked[10,65]=NaN
        rm=scan(a,masked,.04,S)
        @test rm.accepted && rm.applied_dx_px==lag
        @test all(!rm.support[10,x] for x in max(1,65-rm.limit_px):min(128,65+rm.limit_px))
        @test count(rm.support) < count(r.support)
    end
    opposite=scan(a,-a,.04,S)
    @test !opposite.accepted && opposite.applied_dx_px==0
    @test occursin("weak_or_negative",opposite.status)
    edge=scan(a,shifted_signal(a,32),.04,S)
    @test !edge.accepted && occursin("boundary_peak",edge.status)
    constant=scan(ones(64,128),ones(64,128),.04,S)
    @test !constant.accepted && occursin("constant_signal",constant.status)
    missing=scan(fill(NaN,64,128),a,.04,S)
    @test !missing.accepted && occursin("insufficient_support",missing.status)
    periodic=[sin(2pi*x/8) for y in 1:64,x in 1:128]
    periodic_result=scan(periodic,copy(periodic),.04,S)
    @test !periodic_result.accepted && occursin("ambiguous_peak",periodic_result.status)
    b=copy(a)
    for band in 1:4
        ys=(1+(band-1)*16):(band*16)
        b[ys,:] .= shifted_signal(a[ys,:],3band)
    end
    disagree=scan(a,b,.04,S)
    @test !disagree.accepted && occursin("disagreement",disagree.status)
    @test all(p.usable for p in disagree.peaks[2:end])
    # Row offsets do not count as spatial agreement.
    row_only=repeat(collect(1.:64.),1,128)
    @test !scan(row_only,copy(row_only),.04,S).accepted
    @test !scan(a,a,2.,S).accepted
    @test_throws ErrorException scan(a,a,0.,S)
    @test_throws DimensionMismatch scan(a,zeros(3,4),.04,S)
    @test_throws ErrorException common_support(a,a,-1,16)
    @test_throws DimensionMismatch common_support(a,zeros(2,2),1,16)
    m=common_support(a,a,5,16)
    for I in CartesianIndices(a)
        @test m[I] == (6<=I[2]<=123)
    end
end

@testset "Explicit settings, no truth-shaped options" begin
    mktempdir() do dir
        path=joinpath(dir,"config.toml"); cfg=TOML.parsefile(SETTINGS)
        for (section,d) in cfg, k in keys(d)
            bad=deepcopy(cfg); delete!(bad[section],k)
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException load_settings(path)
        end
        for (section,k,value) in (("model","reference_frame","backward"),
            ("model","method","absolute_correlation"),("selection","min_correlation",-1.),
            ("selection","min_band_pixels",true),("selection","bands",1),
            ("selection","min_band_rows",2.5),("selection","min_distant_peak_gap",NaN),
            ("selection","peak_neighborhood_nm",0.),("selection","unresolved_policy","drop"),
            ("preprocessing","missing_policy","impute"),("preprocessing","max_lag_width_fraction",.5))
            bad=deepcopy(cfg); bad[section][k]=value
            open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException load_settings(path)
        end
        bad=deepcopy(cfg); bad["selection"]["expected_N"]=6
        open(io->TOML.print(io,bad),path,"w")
        @test_throws ErrorException load_settings(path)
    end
end

@testset "Transform before Gaussian subtraction, no invented pixels" begin
    z=reshape(collect(1.:99.),9,11); ch=SXMChannel("Z","nm","bwd",copy(z))
    frozen=copy(z)
    for radius in (0,1,2), dx in (-2,0,3)
        pre=(stride=1,smooth_radius_px=radius)
        moved,smoothed=observed_shift(z,ch,pre,dx)
        for y in 1:9,x in 1:11
            @test isequal(moved[y,x],1<=x+dx<=11 ? z[y,x+dx] : NaN)
        end
        @test isequal(smoothed,_box_smooth(moved,radius))
        if dx==0
            @test moved==z && smoothed==_box_smooth(z,radius)
        end
        @test z==frozen && ch.data==frozen
    end
    hole=SXMChannel("Z","nm","bwd",copy(z)); hole.data[5,6]=NaN
    moved,smoothed=observed_shift(z,hole,(stride=1,smooth_radius_px=1),2)
    @test isnan(moved[5,4]) && all(isnan,smoothed[4:6,3:5])
    @test all(isnan,moved[:,10:11])
    @test_throws ErrorException observed_shift(z,ch,(stride=1,smooth_radius_px=0),11)
    @test_throws DimensionMismatch observed_shift(z,ch,(stride=2,smooth_radius_px=0),0)
    @test_throws ErrorException require_direction(ch,"fwd")
    @test require_direction(ch,"bwd")===nothing
    mktempdir() do dir
        p=joinpath(dir,"shifts.tsv")
        write(p,"file\tbwd_sample_dx_px\na.sxm\t-3\n")
        @test read_shifts(p,["a.sxm"])==Dict("a.sxm"=>-3)
        @test read_shifts(nothing,["a.sxm"])===nothing
        @test_throws ErrorException read_shifts(p,["a.sxm","b.sxm"])
        for content in ("file\texpected_N\na.sxm\t6\n","file\tbwd_sample_dx_px\na.sxm\t1.2\n",
            "file\tbwd_sample_dx_px\na.sxm\t1\na.sxm\t2\n","file\tbwd_sample_dx_px\n../a.sxm\t1\n")
            write(p,content); @test_throws ErrorException read_shifts(p,["a.sxm"])
        end
    end
end

module Comparison
include(joinpath(@__DIR__,"run_acquisition_registration_comparison.jl"))
end

@testset "Three full arms, no refit, dry-run and launcher boundaries" begin
    C=Comparison
    @test_throws ErrorException C.registration_comparison(["--benchmark-manifest","truth"])
    @test_throws ErrorException C.EstimateAcquisitionShifts.parse_cli(["--expected-N","6"])
    mktempdir() do dir
        raw=joinpath(dir,"raw"); mkdir(raw); touch(joinpath(raw,"a.sxm"))
        base=joinpath(dir,"base.tsv"); split=joinpath(dir,"split.tsv"); templates=joinpath(dir,"templates.tsv")
        write(base,"file\tlobe\tN\na.sxm\t1\t1\n"); cp(base,split); touch(templates)
        out=joinpath(dir,"run")
        args=["--data-dir",raw,"--count-config",joinpath(C.ROOT,"config","chitosan.toml"),
            "--config",joinpath(C.ROOT,"config","unit_assignment_patch_support.toml"),
            "--settings",SETTINGS,"--features",base,"--split-features",split,"--templates",templates,"--outdir",out]
        C.registration_comparison(vcat(args,["--dry-run"]))
        @test !ispath(out)
        calls=[]
        runner=function(root,name,script,a;threads=4)
            push!(calls,(name,script,a))
            op=Dict(a[i]=>a[i+1] for i in 1:2:length(a))
            @test read(op["--features"])==read(base) && read(op["--split-features"])==read(split)
            mkdir(op["--outdir"]); cp(base,joinpath(op["--outdir"],"predictions.tsv"))
        end
        estimator=function(root,opts,rows;runner)
            write(joinpath(root,"shifts.tsv"),"file\tbwd_sample_dx_px\na.sxm\t2\n")
        end
        C.registration_comparison(args;runner,estimator)
        @test first.(calls)==["reference","observed_control","registered"]
        @test all(c[2]=="run_reconstructed_chitosan.jl" for c in calls)
        @test !("--acquisition-shifts" in calls[1][3])
        @test "--acquisition-shifts" in calls[2][3] && "--acquisition-shifts" in calls[3][3]
        @test read(joinpath(out,"zero_shifts.tsv"),String)=="file\tbwd_sample_dx_px\na.sxm\t0\n"
        @test_throws ErrorException C.registration_comparison(args;runner,estimator)
        @test !isfile(joinpath(out,"failures.tsv"))
        for missing in ("--features","--split-features"), cached in (nothing,"--patches-bwd")
            op=Dict("--config"=>joinpath(C.ROOT,"config","unit_assignment_patch_support.toml"),
                "--features"=>base,"--split-features"=>split,"--acquisition-shifts"=>joinpath(out,"shifts.tsv"))
            delete!(op,missing)
            cached===nothing || (op[cached]="cache.tsv")
            @test_throws ErrorException C.execute_pipeline(op)
        end
        for mode in ("reference","observed_control","registered")
            @test read(joinpath(out,mode,"features.tsv"))==read(base)
            @test read(joinpath(out,mode,"features_split.tsv"))==read(split)
        end
        launcher=joinpath(C.ROOT,"hpc","compare_acquisition_registration.sbatch")
        @test success(`bash -n $launcher`)
        script=read(launcher,String)
        @test occursin("#SBATCH --time=02:00:00",script)
        @test occursin("SLURM_JOB_ID",script) && occursin("--dry-run",script)
        env=("STMFIT_PROJECT_DIR"=>C.ROOT,"STMFIT_CACHE_DIR"=>dir,"STMFIT_INPUT_DIR"=>dir,
            "STMFIT_OUTDIR"=>joinpath(dir,"another"),"JULIA_BIN"=>String(first(Base.julia_cmd().exec)))
        @test !success(pipeline(addenv(`bash $launcher`,env...,"SLURM_JOB_ID"=>"");stdout=devnull,stderr=devnull))
        # Multiple registration workers partition by file; all products merge
        # without filtering a file based on its acceptance/assignment.
        shardout=joinpath(dir,"shards"); mkpath(joinpath(shardout,"logs"))
        bpath=joinpath(dir,"five.tsv")
        C.write_table(bpath,["file","lobe","N"],[Dict("file"=>"f$i.sxm","lobe"=>"1","N"=>"1") for i in 1:5])
        _,bb=C.lobe_table(bpath); observed=String[]; nworkers=Ref(0)
        shardrunner=function(root,name,script,a;threads)
            @test threads==1; nworkers[]+=1
            op=Dict(a[i]=>a[i+1] for i in 1:2:length(a))
            _,part=C.lobe_table(op["--features"]); fs=sort(unique(first.(collect(keys(part)))))
            append!(observed,fs)
            for suffix in ("",".summary.tsv",".peaks.tsv",".scores.tsv")
                C.write_table(op["--out"]*suffix,["file","bwd_sample_dx_px"],
                    [Dict("file"=>f,"bwd_sample_dx_px"=>"0") for f in fs])
            end
        end
        opts=Dict("--features"=>bpath,"--data-dir"=>raw,"--count-config"=>args[4],"--settings"=>SETTINGS)
        C.estimate_shards(shardout,opts,bb;runner=shardrunner)
        @test nworkers[]==min(4,Threads.nthreads())
        @test sort(observed)==["f$i.sxm" for i in 1:5]
        @test length(C.EstimateAcquisitionShifts.PatchAcquisition.read_shifts(joinpath(shardout,"shifts.tsv"),observed))==5
        badargs=copy(args); badargs[end]=joinpath(dir,"failed")
        rejector=(root,opts,rows;runner)->error("synthetic acquisition failure")
        @test_throws ErrorException C.registration_comparison(badargs;runner,estimator=rejector)
        @test isfile(joinpath(dir,"failed","failures.tsv"))
        @test !isdir(joinpath(dir,"failed","reference"))
    end
end

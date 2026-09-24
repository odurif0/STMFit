# Synthetic fixtures only; no raw scan, benchmark or historical prediction.
using Test, Random, Statistics, LinearAlgebra, TOML, Printf
include(joinpath(@__DIR__,"verify_frozen_learning.jl"))
const D=FrozenLearningVerification.D
const CFG=joinpath(D.ROOT,"config","unit_assignment_patch_support.toml")
const SETTINGS=joinpath(D.ROOT,"config","frozen_learning.toml")

function feature_fixture(mod)
    rng=MersenneTwister(724)
    [mod.LobeRecord("synthetic_$f.sxm",j,10*f+2*(j>3)+.1randn(rng),
        Dict("x"=>(j>3 ? 2. : -1.)+.2randn(rng),"y"=>randn(rng),"constant"=>3.)) for f in 1:8 for j in 1:7]
end

@testset "Frozen statistics and heads: exact native reapplication, no target renaming" begin
    mktempdir() do dir
        for kind in ("gmm","kmeans")
            mod=kind=="gmm" ? D.Gaussian : D.KMeans
            rec=feature_fixture(mod); rec[end].features["x"]=NaN
            featurepath=joinpath(dir,kind*".tsv")
            D.RU.write_table(featurepath,["file","lobe","amplitude","x","y","constant"],
                [merge(Dict("file"=>r.file,"lobe"=>string(r.lobe),"amplitude"=>string(r.amplitude)),
                    Dict(f=>string(v) for (f,v) in r.features)) for r in rec])
            args=["--features",featurepath,"--out",joinpath(dir,"out.tsv"),"--seeds","2","--interactions"]
            kind=="gmm" && append!(args,["--config",CFG,"--selftrain","2"])
            opt=mod._parse_cli(args); features=["x","y","constant"]
            native=mod._view_probability(rec,features,opt)
            h=D.capture_head(kind,rec,features,opt)
            @test isequal(native,h.native)
            @test length(h.models)==2
            @test length(h.scales)==24
            @test all(s.scale==1. for s in h.scales if s.feature=="constant")
            X,valid=D.frozen_matrix(rec,features,h.scales;interactions=true)
            @test size(X)==(56,6) && !valid[end]
            score(r)=begin
                z,v=D.frozen_matrix(r,features,h.scales;interactions=true)
                kind=="gmm" ? D.frozen_gmm(z,v,h.models,opt) : D.frozen_kmeans(z,v,h.models)
            end
            @test isequal(score(reverse(rec)),reverse(native))
            target=deepcopy(rec)
            for r in target; r.amplitude=1e6-r.amplitude; end
            @test isequal(score(target),native) # No renaming from target amplitudes.
            target[1].features["x"]+=.7
            @test isequal(score(target)[2:end],native[2:end])
            @test score(target)[end] |> isnan
            @test_throws ErrorException D.frozen_matrix(rec,features,vcat(h.scales,h.scales[1:1]);interactions=true)
            @test_throws ErrorException D.frozen_matrix(rec,features,h.scales[2:end];interactions=true)
            @test all(isnan,kind=="gmm" ? D.frozen_gmm(X,valid,[],opt) : D.frozen_kmeans(X,valid,[]))
            if kind=="gmm"
                # Independent quadratic-form scores, including component weights.
                for m in h.models, i in findall(valid)
                    scores=[log(m.weights[c])-.5dot(X[i,:]-m.means[:,c],
                        (Symmetric(m.covariances[c])+1e-8I)\(X[i,:]-m.means[:,c])) for c in 1:2]
                    expected=Float64(argmax(scores)==m.high_cluster)
                    @test D.frozen_gmm(X[i:i,:],[true],[m],opt)[1]==expected
                end
            end
        end
    end
end

function patches_fixture()
    options=D.EF.load_fisher_config(CFG); grid=D.EF.fisher_grid(options); rng=MersenneTwister(935)
    keys=[("synthetic_$f.sxm",j) for f in 1:4 for j in 1:12]
    full=.05randn(rng,length(keys),grid.side^2)
    shape=[exp(-18(u*u+t*t))+.1sin(9t) for u in grid.coords for t in grid.coords]
    for (i,k) in enumerate(keys); full[i,:].+=(k[2]>7 ? 2. : -1.).*shape; end
    table=D.EF.PatchTable(keys,full[:,grid.disk_indices],full[:,grid.center_index],fill("",length(keys)),grid)
    options,table,full
end

@testset "Frozen Fisher: opposite-fold application and retained invalid keys" begin
    o,p,_=patches_fixture(); f=D.capture_fisher(p,o)
    @test length(f.models)==2
    @test [r.score for r in f.native]==[r.score for r in D.frozen_fisher(p,f.models)]
    before=deepcopy(f.models); x=copy(p.X); x[1,:].+=1.
    changed=D.EF.PatchTable(p.keys,x,p.amplitudes,p.invalid_reasons,p.grid)
    q=D.frozen_fisher(changed,f.models)
    @test [r.score for r in q[2:end]]==[r.score for r in f.native[2:end]]
    @test f.models==before
    for parity in 0:1
        m=f.models[parity]
        @test all(mod(p.keys[i][2],2)==parity for i in m.train)
        @test all(mod(p.keys[i][2],2)!=parity for i in m.held)
    end
    reasons=copy(p.invalid_reasons); reasons[1]="missing_pixel"
    bad=D.EF.PatchTable(p.keys,x,p.amplitudes,reasons,p.grid)
    q=D.frozen_fisher(bad,f.models)
    @test q[1].invalid_reason=="missing_pixel" && isnan(q[1].score)
    @test [(r.file,r.lobe) for r in q]==p.keys
    @test_throws ErrorException D.frozen_fisher(p,Dict())
end

function whole_fixture(root)
    mkpath(root); cfg=D.RU.load_config(CFG); cp(CFG,joinpath(root,"assignment.toml"))
    _,p,full=patches_fixture(); rng=MersenneTwister(553)
    rows=Dict{String,String}[]; splits=Dict{String,String}[]; scores=Dict{String,String}[]
    for (i,k) in enumerate(p.keys)
        r=Dict("file"=>k[1],"lobe"=>string(k[2]),"amplitude"=>string(10+k[2]+randn(rng)))
        for f in vcat(D.BASE,["patch_u_asym_reconstructed"]); r[f]=string(.2k[2]+randn(rng)); end
        push!(rows,r); push!(splits,Dict("file"=>k[1],"lobe"=>string(k[2]),"skew_ratio"=>string(.8+.4rand(rng))))
        push!(scores,Dict("file"=>k[1],"lobe"=>string(k[2]),"cost_margin"=>string(rand(rng))))
    end
    headers=vcat(["file","lobe","amplitude"],D.BASE,["patch_u_asym_reconstructed"])
    for arm in ("control","observed")
        dir=joinpath(root,arm); mkdir(dir); localrows=deepcopy(rows); localfull=copy(full)
        if arm=="observed"
            for i in 25:48
                localrows[i]["amp_rel"]=string(parse(Float64,localrows[i]["amp_rel"])+.13randn(rng))
                localfull[i,:].+=.1randn(rng,size(full,2))
            end
        end
        D.RU.write_table(joinpath(dir,"features_descriptor.tsv"),headers,localrows)
        D.RU.write_table(joinpath(dir,"features_split.tsv"),["file","lobe","skew_ratio"],splits)
        for stage in ("score_fwd","score_bwd"); D.RU.write_table(joinpath(dir,stage*".tsv"),["file","lobe","cost_margin"],scores); end
        for (stage,prefix) in (("patches_fwd17","res"),("patches_bwd17","bwd_res"))
            names=[@sprintf("%s_p%03d",prefix,j) for j in 1:size(full,2)]
            prs=[merge(Dict("file"=>k[1],"lobe"=>string(k[2])),Dict(names[j]=>string(localfull[i,j]) for j in eachindex(names))) for (i,k) in enumerate(p.keys)]
            D.RU.write_table(joinpath(dir,stage*".tsv"),vcat(["file","lobe"],names),prs)
        end
        D.EF.write_scores(joinpath(dir,"fisher_cv.tsv"),[D.EF.FisherScore(k[1],k[2],0.,"") for k in p.keys])
        D.joined(root,arm,joinpath(dir,"fisher_cv.tsv"),joinpath(dir,"features_predictor.tsv"),cfg)
        dummy=[Dict("file"=>k[1],"lobe"=>string(k[2]),"predicted"=>"0","probability_1"=>"0.1","confidence"=>"0.8","invalid_reason"=>"ok") for k in p.keys]
        for stage in ("pred_gmm","pred_kmeans","predictions")
            D.RU.write_table(joinpath(dir,stage*".tsv"),["file","lobe","predicted","probability_1","confidence","invalid_reason"],dummy)
        end
    end
    files=sort(unique(first.(p.keys)))
    D.RU.write_table(joinpath(root,"selected_summary.tsv"),["filepath","N_selected"],[Dict("filepath"=>f,"N_selected"=>"12") for f in files])
    audits=[Dict("file"=>f,"arm"=>a,"profile"=>pr,"valid"=>"true","reused_control"=>string(a=="observed" && f in files[1:2]),
        "observed_pixels"=>f in files[1:2] ? "100" : "90","total_pixels"=>"100") for f in files for a in ("control","observed") for pr in ("gaussian","split")]
    D.RU.write_table(joinpath(root,"refit_chunk1.tsv.fits.tsv"),["file","arm","profile","valid","reused_control","observed_pixels","total_pixels"],audits)
    root
end

@testset "Declared six-arm complete diagnostic, saved-state round trip and strict CLI" begin
    @test D.settings(SETTINGS)["model"]["repetitions"]==2
    @test length(D.ARMS)==6
    @test D.parse_cli(["--help"])===nothing
    mktempdir() do dir
        root=whole_fixture(joinpath(dir,"input")); out=joinpath(dir,"out")
        args=["--input",root,"--settings",SETTINGS,"--outdir",out]
        data=D.inputs(root)
        @test length(data.counts)==4 && sum(values(data.counts))==48
        @test count(==("fully_observed"),values(data.groups))==2
        @test D.main(vcat(args,["--dry-run"]))===nothing
        @test !ispath(out)
        launcher=joinpath(D.ROOT,"hpc","diagnose_frozen_learning.sbatch")
        env=["STMFIT_PROJECT_DIR"=>D.ROOT,"STMFIT_INPUT_DIR"=>root,"STMFIT_OUTDIR"=>joinpath(dir,"hpc"),
            "JULIA_BIN"=>first(Base.julia_cmd().exec)]
        @test success(addenv(`bash $launcher --dry-run`,env...))
        @test !ispath(joinpath(dir,"hpc"))
        @test !success(pipeline(ignorestatus(addenv(`bash $launcher`,env...,"SLURM_JOB_ID"=>""));stdout=devnull,stderr=devnull))
        @test occursin("#SBATCH --no-requeue",read(launcher,String))
        @test occursin("#SBATCH --time=01:00:00",read(launcher,String))
        for forbidden in ("--truth","--expected-n","--seed","--threshold","--manifest")
            @test_throws ErrorException D.parse_cli(vcat(args,[forbidden,"forbidden"]))
        end
        @test_throws ErrorException D.parse_cli(vcat(args,["--input",root]))
        D.execute(root,SETTINGS,out)
        @test all(isfile(joinpath(out,a.name,"predictions.tsv")) for a in D.ARMS)
        @test all(D.sha(p)==h for (p,h) in data.hashes)
        @test_throws ErrorException D.execute(root,SETTINGS,out)
        @test_throws ErrorException D.parse_cli(args)
        for a in ("control","observed"), stage in D.REFERENCES
            @test read(joinpath(out,a,stage*".tsv"))==read(joinpath(out,"native_"*a,stage*".tsv"))
        end
        for a in D.ARMS
            @test length(D.Pipeline.check_counts(joinpath(out,a.name,"predictions.tsv"),data.counts))==48
        end
        @test FrozenLearningVerification.verify(root,out)===nothing
        # Do not drop a failed or missing input row to make the attribution pass.
        open(joinpath(root,"control","score_fwd.tsv"),"a") do io
            println(io,"synthetic_1.sxm\t1\t0.2")
        end
        @test_throws ErrorException D.inputs(root)
    end
end

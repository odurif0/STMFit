# Also rerun the legacy finite/shared regression suite unchanged.
include(joinpath(@__DIR__,"test_mold_state_decode.jl"))
include(joinpath(@__DIR__,"run_mold_loo_comparison.jl"))
const LOO=MoldLOOComparison
const L=LOO.L
const LSET=joinpath(ROOT,"config/mold_leave_one_out.toml")
const LOPT=L.settings(LSET)

@testset "Leave-target-out state enumeration and no target influence" begin
    rng=MersenneTwister(92401)
    for n in 2:7, rep in 1:8
        _,rs,v=fixture(n)
        for c in values(v); c .= randn(rng,2); end
        original=deepcopy(v); decoded=L.decode(rs,v,LOPT)
        @test v==original
        for i in 1:n
            others=filter(j->j!=i,1:n)
            # Independent exhaustive chemical assignments for every geometric state.
            oracle=[minimum(sum(v[(rs[j].file,j,mod(j-1+p,2),m)][((mask>>(a-1))&1)+1]
                for (a,j) in enumerate(others)) for mask in 0:(2^(n-1)-1)) for p in 0:1 for m in 0:1]
            states=[(p,m) for p in 0:1 for m in 0:1]; best=states[argmin(oracle)]; b=decoded[i]
            @test (b.phase,b.mirror)==best
            @test b.objectives≈oracle atol=1e-12
            @test b.training_cost≈minimum(oracle) atol=1e-12
            @test b.training_lobes==n-1
            @test b.costs==v[(rs[i].file,i,mod(i-1+b.phase,2),b.mirror)]
            @test b.label==Int(b.costs[2]<b.costs[1]) && b.reason=="ok"
            changed=deepcopy(v)
            for (key,c) in changed
                key[2]==i && (c .= randn(rng,2)*1e180)
            end
            hold=L.decode(rs,changed,LOPT)[i]
            @test (hold.phase,hold.mirror,hold.objectives,hold.training_cost,hold.training_lobes)==
                  (b.phase,b.mirror,b.objectives,b.training_cost,b.training_lobes)
        end
    end
    _,rs,v=fixture(3)
    for (k,c) in v
        phase=mod(k[3]-k[2]+1,2)
        c .= [-0.1,0.1]
        k[2]==1 && phase==0 && k[4]==0 && (c .= [-10.,10.])
        k[2]!=1 && phase==1 && k[4]==1 && (c .= [-1.,1.])
    end
    full=M.decode(rs,[v,v],OPT;mode=:finite)[1]
    held=L.decode(rs,v,LOPT)
    @test (full.phase,full.mirror)==(0,0)
    @test (held[1].phase,held[1].mirror)==(1,1)
    @test all(b->(b.phase,b.mirror)==(0,0),held[2:3])
end

@testset "Missing evidence, ties, singleton and transparent serialization" begin
    _,rs,v=fixture(1); b=only(L.decode(rs,v,LOPT))
    @test b.reason=="no_remaining_observations" && b.training_lobes==0
    @test b.phase==b.mirror==b.label==-1 && b.costs==[Inf,Inf]
    _,rs,v=fixture(3); b=L.decode(rs,v,LOPT)
    @test all(x->x.phase==x.mirror==0 && x.training_lobes==2,b)
    for (k,c) in v; k[2]!=1 && (c .= Inf); end
    b=L.decode(rs,v,LOPT)
    @test b[1].reason=="no_remaining_observations" && b[1].costs==[Inf,Inf]
    @test all(x->x.reason=="insufficient_native_support" && x.training_lobes==1 && x.label==-1,b[2:3])
    io=IOBuffer(); L.write_decoded(io,rs,b); txt=String(take!(io))
    @test occursin("\t?\tInf\tInf\tNaN\t",txt)
    @test !any(startswith(c,"global_") || c=="file_cost" for c in L.SCORE_HEADER)
    v[(rs[1].file,2,0,0)]=[0.,1.]
    @test_throws ErrorException L.decode(rs,v,LOPT)
    v[(rs[1].file,2,0,0)]=[NaN,Inf]
    @test_throws ErrorException L.decode(rs,v,LOPT)
    @test_throws ErrorException L.decode(rs,v,merge(LOPT,(transition_penalty=0.1,)))
    @test_throws ErrorException L.decode(rs[2:end],v,LOPT)
    for field in ("empty_training","state_tie","missing_cost")
        mktempdir() do dir
            cfg=TOML.parsefile(LSET); cfg["selection"][field]="target_fallback"
            p=joinpath(dir,"bad.toml"); open(io->TOML.print(io,cfg),p,"w")
            @test_throws ErrorException L.settings(p)
        end
    end
end

@testset "Two-arm driver, independent checker, exact control and errors" begin
    mktempdir() do dir
        ref=joinpath(dir,"saved"); mkpath(joinpath(ref,"logs"))
        base,rs,v=fixture(3)
        for r in values(base)
            r["skew_ratio"]="1.0"; r["patch_u_asym_reconstructed"]="0.3"; r["score"]="0.2"
        end
        h=sort(collect(keys(first(values(base)))))
        for f in D.INPUTS; D.write_table(joinpath(ref,f),h,[base[k] for k in sort(collect(keys(base)))]); end
        ah=["file","lobe","parity","mirror","cost0","cost1","reason"]
        rows=[Dict("file"=>k[1],"lobe"=>string(k[2]),"parity"=>string(k[3]),"mirror"=>string(k[4]),
            "cost0"=>string(c[1]),"cost1"=>string(c[2]),"reason"=>"ok") for (k,c) in sort(collect(v))]
        for view in ("fwd","bwd"); D.write_table(joinpath(ref,"score_"*view*".tsv.audit.tsv"),ah,rows); end
        calls=[]
        runner=function(out,stage,script,args)
            push!(calls,(stage,script,args))
            if stage in ("gmm","kmeans")
                path=args[findfirst(==("--out"),args)+1]
                D.write_table(path,["file","lobe","predicted","probability_1"],
                    [Dict("file"=>r.file,"lobe"=>string(r.lobe),"predicted"=>"0","probability_1"=>"0.1") for r in rs])
            else
                @test stage=="validation"
            end
        end
        config=joinpath(ROOT,"config/unit_assignment_patch_support.toml")
        args=["--reference-dir",ref,"--config",config,"--settings",LSET,"--outdir",joinpath(dir,"out")]
        o=D.parse_cli(args); data=D.load_inputs(o;read_settings=L.settings)
        D.write_scores(ref,data,:finite)
        classifier=(dest,cfg,path)->D.classify(dest,cfg,path;runner)
        classifier(ref,data.cfg,config); empty!(calls)
        LOO.execute(merge(o,Dict("--dry-run"=>"true"));classifier)
        @test !ispath(o["--outdir"]) && isempty(calls)
        LOO.execute(o;classifier)
        @test [c[1] for c in calls]==repeat(["gmm","kmeans","validation"],2)
        for f in (D.INPUTS...,D.REPLAY_FILES...)
            @test read(joinpath(o["--outdir"],"reference",f))==read(joinpath(ref,f))
        end
        for f in D.INPUTS
            @test read(joinpath(o["--outdir"],"loo",f))==read(joinpath(ref,f))
        end
        _,audit=D.read_table(joinpath(o["--outdir"],"loo","state_selection.tsv"))
        @test length(audit)==24 && count(r->r["selected"]=="true",audit)==6
        checker=joinpath(ROOT,"test/verify_mold_loo_comparison.jl")
        @test success(`$(Base.julia_cmd()) --project=$ROOT $checker $ref $(o["--outdir"])`)
        @test_throws ErrorException LOO.execute(o;classifier)
        @test_throws ErrorException D.parse_cli(vcat(args,["--truth","forbidden"]))
        fail=merge(o,Dict("--outdir"=>joinpath(dir,"failed")))
        @test_throws ErrorException LOO.execute(fail;classifier=(a...)->error("Synthetic failure"))
        @test isfile(joinpath(fail["--outdir"],"failures.tsv"))
        @test !ispath(joinpath(fail["--outdir"],"loo"))
    end
    shell=joinpath(ROOT,"hpc/compare_mold_loo.sbatch")
    @test success(`bash -n $shell`)
    @test occursin("#SBATCH --time=02:00:00",read(shell,String))
    @test occursin("SLURM_JOB_ID:?",read(shell,String))
end

using Test, Random, TOML
include(joinpath(@__DIR__,"run_mold_state_comparison.jl"))
const D=MoldStateComparison
const M=D.M
const ROOT=dirname(@__DIR__)
const SETTINGS=joinpath(ROOT,"config/mold_state_comparison.toml")
const OPT=M.settings(SETTINGS)

function fixture(n)
    base=Dict(("synthetic.sxm",i)=>Dict("file"=>"synthetic.sxm","lobe"=>string(i),"N"=>string(n),"amplitude"=>"0.1") for i in 1:n)
    records=M.records_for(base,"synthetic.sxm")
    costs=Dict(("synthetic.sxm",i,p,m)=>[-0.1,0.1] for i in 1:n for p in 0:1 for m in 0:1)
    return base,records,costs
end
function scoretext(records,best)
    io=IOBuffer(); M.Connected.write_decoded(io,first(records).file,records,best); String(take!(io))
end

@testset "Missing lobes supply no evidence, not imputed predictions" begin
    base,rs,v=fixture(4)
    for (k,c) in v
        c .= k[2]<=2 ? [Inf,Inf] : k[4]==1 ? [-0.8,0.8] : [-0.1,0.1]
    end
    original=deepcopy(v)
    legacy=M.decode(rs,[v,v],OPT;mode=:reference)
    corrected=M.decode(rs,[v,v],OPT;mode=:finite)
    paired=M.decode(rs,[v,v],OPT;mode=:shared)
    @test all(b->b.total==Inf && b.mirror==0,legacy)
    @test all(b->b.total≈-1.6 && b.mirror==1,corrected)
    @test all(b->b.labels==[-1,-1,0,0],corrected)
    @test all(b->b.costs[1:2,:]==fill(Inf,2,2),corrected)
    @test isequal(v,original)
    @test scoretext(rs,corrected[1])==scoretext(rs,paired[1])
    @test occursin("\t-1\t0.1\t?\tInf\tInf\tNaN\t",scoretext(rs,corrected[1]))
    for c in values(v); c .= Inf; end
    for mode in (:finite,:shared)
        b=M.decode(rs,[v,v],OPT;mode)[1]
        @test b.total==0 && b.phase==0 && b.mirror==0
        @test b.labels==fill(-1,4)
        @test all(==(Inf),b.costs)
    end
    invalid=deepcopy(v); invalid[("synthetic.sxm",1,0,0)]=[NaN,Inf]
    @test_throws ErrorException M.decode(rs,[invalid,v],OPT;mode=:finite)
    invalid[("synthetic.sxm",1,0,0)]=[0.0,Inf]
    @test_throws ErrorException M.decode(rs,[invalid,v],OPT;mode=:shared)
    opt=M.Connected.Options("","",nothing,"","","ncc","contrast",0.1,1.0)
    @test_throws ErrorException M.Connected._decode_state(zeros(4,2),rs,0,0,0,nothing,opt;omit_unavailable=true)
end

@testset "Independent brute force, direction equivalence and view symmetry" begin
    rng=MersenneTwister(23923)
    for n in 1:7, repeat in 1:5
        _,rs,a=fixture(n); b=deepcopy(a)
        for v in (a,b),c in values(v); c .= randn(rng,2); end
        independent=M.decode(rs,[a,b],OPT;mode=:finite)
        legacy=M.decode(rs,[a,b],OPT;mode=:reference)
        @test all(scoretext(rs,x)==scoretext(rs,y) for (x,y) in zip(independent,legacy))
        best=M.decode(rs,[a,b],OPT;mode=:shared)
        swapped=M.decode(rs,[b,a],OPT;mode=:shared)
        @test scoretext(rs,best[1])==scoretext(rs,swapped[2])
        @test scoretext(rs,best[2])==scoretext(rs,swapped[1])
        totals=Dict()
        for p in 0:1,m in 0:1
            # Exhaust all chemical assignments, independently for each view.
            totals[(p,m)]=sum(minimum(sum(v[(r.file,r.lobe,mod(r.lobe-1+p,2),m)][((mask>>(r.lobe-1))&1)+1]
                for r in rs) for mask in 0:(2^n-1)) for v in (a,b))
        end
        p,m=best[1].phase,best[1].mirror
        @test (p,m)==(best[2].phase,best[2].mirror)
        @test sum(x.total for x in best)≈minimum(values(totals)) atol=1e-12
        @test totals[(p,m)]≈minimum(values(totals)) atol=1e-12
        for direction in 0:1,phase in 0:1,i in 1:n
            @test M.Connected._parity_for_lobe(i,n,direction,phase)==mod(i-1+phase+direction*(n-1),2)
        end
    end
    _,rs,a=fixture(3); b=deepcopy(a)
    for c in values(a); c .= [-0.9,0.9]; end
    for c in values(b); c .= [0.8,-0.8]; end
    best=M.decode(rs,[a,b],OPT;mode=:shared)
    @test best[1].labels==[0,0,0] && best[2].labels==[1,1,1]
    @test best[1].phase==0 && best[1].mirror==0
    # One absent view reduces exactly to the available view, including margins.
    for c in values(b); c .= Inf; end
    best=M.decode(rs,[a,b],OPT;mode=:shared)
    finite=M.decode(rs,[a,b],OPT;mode=:finite)
    @test scoretext(rs,best[1])==scoretext(rs,finite[1])
    @test best[2].labels==fill(-1,3)
    @test_throws ErrorException M.decode(rs,[a,b],OPT;mode=:unknown)
end

@testset "Strict complete cost-cache and CLI boundaries" begin
    mktempdir() do dir
        base,rs,v=fixture(3)
        rows=[Dict("file"=>k[1],"lobe"=>string(k[2]),"parity"=>string(k[3]),"mirror"=>string(k[4]),
            "cost0"=>string(c[1]),"cost1"=>string(c[2]),"reason"=>"ok") for (k,c) in sort(collect(v))]
        h=["file","lobe","parity","mirror","cost0","cost1","reason"]
        path=D.write_table(joinpath(dir,"score_fwd.tsv.audit.tsv"),h,rows)
        @test M.load_costs([path],base)==v
        @test M.audit_paths(dir,"fwd")==[path]
        @test_throws ErrorException M.audit_paths(dir,"bwd")
        @test_throws ErrorException M.load_costs([path,path],base)
        partial=D.write_table(joinpath(dir,"partial.tsv"),h,rows[2:end])
        @test_throws ErrorException M.load_costs([partial],base)
        missingrows=deepcopy(rows)
        for r in missingrows
            r["lobe"]=="1" || continue
            r["cost0"]="Inf";r["cost1"]="Inf";r["reason"]="insufficient_native_support"
        end
        missingpath=D.write_table(joinpath(dir,"missing.tsv"),h,missingrows)
        @test count(c->all(==(Inf),c),values(M.load_costs([missingpath],base)))==4
        missingrows[1]["cost0"]="NaN"
        invalid=D.write_table(joinpath(dir,"nan.tsv"),h,missingrows)
        @test_throws ErrorException M.load_costs([invalid],base)
        mixed=deepcopy(rows); mixed[1]["cost0"]="Inf"; mixed[1]["cost1"]="Inf"; mixed[1]["reason"]="insufficient_native_support"
        mixedpath=D.write_table(joinpath(dir,"mixed.tsv"),h,mixed)
        @test_throws ErrorException M.load_costs([mixedpath],base)
        forbidden=D.write_table(joinpath(dir,"forbidden.tsv"),vcat(h,["expected_N"]),rows)
        @test_throws ErrorException M.load_costs([forbidden],base)
        config=joinpath(ROOT,"config/unit_assignment_patch_support.toml")
        args=["--reference-dir",dir,"--config",config,"--settings",SETTINGS,"--outdir",joinpath(dir,"out")]
        @test D.parse_cli(args)["--settings"]==SETTINGS
        @test_throws ErrorException D.parse_cli(vcat(args,["--truth","forbidden"]))
        @test_throws ErrorException D.parse_cli(vcat(args,["--settings",SETTINGS]))
        @test_throws ErrorException D.parse_cli(replace.(args,joinpath(dir,"out")=>dir))
        for section in ("model","selection","preprocessing")
            c=TOML.parsefile(SETTINGS); c[section]["unexpected"]=1
            p=joinpath(dir,section*".toml"); open(io->TOML.print(io,c),p,"w")
            @test_throws ErrorException M.settings(p)
        end
    end
    shell=joinpath(ROOT,"hpc/compare_mold_states.sbatch")
    @test success(`bash -n $shell`)
    @test occursin("#SBATCH --time=02:00:00",read(shell,String))
    @test occursin("SLURM_JOB_ID:?",read(shell,String))
end

@testset "Three-arm driver and unchanged classifier commands, synthetic cache" begin
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
                @test stage=="validation" && script=="validate_unit_predictions.jl"
            end
        end
        config=joinpath(ROOT,"config/unit_assignment_patch_support.toml")
        args=["--reference-dir",ref,"--config",config,"--settings",SETTINGS,"--outdir",joinpath(dir,"run")]
        o=D.parse_cli(args); data=D.load_inputs(o)
        D.write_scores(ref,data,:reference)
        classifier=(dest,cfg,path)->D.classify(dest,cfg,path;runner)
        classifier(ref,data.cfg,config)
        @test [c[1] for c in calls]==["gmm","kmeans","validation"]
        @test calls[1][3][findfirst(==("--seeds"),calls[1][3])+1]=="10"
        @test calls[2][3][findfirst(==("--seeds"),calls[2][3])+1]=="20"
        @test "--interactions" in calls[1][3] && "--interactions" in calls[2][3]
        @test !("--training-support" in calls[1][3])
        empty!(calls)
        D.execute(merge(o,Dict("--dry-run"=>"true"));classifier)
        @test !ispath(o["--outdir"]) && isempty(calls)
        D.execute(o;classifier)
        @test [c[1] for c in calls]==repeat(["gmm","kmeans","validation"],3)
        for mode in ("reference","finite","shared"),f in (D.INPUTS...,D.REPLAY_FILES...)
            @test read(joinpath(o["--outdir"],mode,f))==read(joinpath(ref,f))
        end
        checker=joinpath(@__DIR__,"verify_mold_state_comparison.jl")
        @test success(`$(Base.julia_cmd()) --project=$ROOT $checker $ref $(o["--outdir"])`)
        @test_throws ErrorException D.parse_cli(args)
        fail=merge(o,Dict("--outdir"=>joinpath(dir,"failed")))
        @test_throws ErrorException D.execute(fail;classifier=(a...)->error("Synthetic failure"))
        @test isfile(joinpath(fail["--outdir"],"failures.tsv"))
        @test !ispath(joinpath(fail["--outdir"],"finite"))
    end
end

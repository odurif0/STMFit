#!/usr/bin/env julia
# Independent model arithmetic on saved outputs, without optimization or labels.
module VerifyPairedShift
using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_shift.jl"))
const D=PairedShiftDiagnostic
const H=D.H
const G=H.G
val(r,k)=parse(Float64,r[k])
eq(a,b)=isapprox(a,b;atol=1e-12,rtol=1e-9)

"Independent sum of decoded elliptical/split peaks, with the backward coordinate sign explicit."
function predict(q,ctx,mode,step)
    k=length(ctx.initial.params); p=q[1:k]; d=ctx.data; a=d.axisctx
    _,features,_,_,_,_=G._decode_chain(p,ctx.initial.n,a,ctx.cfg;amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
    dx,dy=mode=="shifted" ? (q[end-1]*step[1],q[end]*step[2]) : (0.,0.)
    g,d0,tx,ty=mode=="fused" ? (0.,0.,0.,0.) : Tuple(q[k+1:k+4])
    peak=function(x,y)
        result=zeros(length(x))
        for f in features,j in eachindex(x)
            t=(x[j]-f.x_nm)*a.axis[1]+(y[j]-f.y_nm)*a.axis[2]
            u=-(x[j]-f.x_nm)*a.axis[2]+(y[j]-f.y_nm)*a.axis[1]
            width=ctx.cfg.peak_profile==:gaussian ? f.sigma_x_nm : t<0 ? f.sigma_x_nm/sqrt(f.skew_ratio) : f.sigma_x_nm*sqrt(f.skew_ratio)
            result[j]+=f.amplitude*exp(-.5*((t/width)^2+(u/f.sigma_y_nm)^2))
        end
        result
    end
    f=(1+g).*peak(d.x,d.y).+(p[1]+d0).+(p[2]+tx).*d.x.+(p[3]+ty).*d.y
    b=(1-g).*peak(d.x.-dx,d.y.-dy).+(p[1]-d0).+(p[2]-tx).*d.x.+(p[3]-ty).*d.y
    return (mean=(f.+b)./2,fwd=f,bwd=b)
end

function verify(native,reference,out)
    opts=Dict("--native-dir"=>native,"--reference-dir"=>reference,"--outdir"=>out)
    for key in ("physical-config","assignment-config","model-settings","settings","solver-settings","shift-settings")
        opts["--"*key]=joinpath(out,key*".toml")
    end
    input=D.inputs(opts); rows=D.table(joinpath(out,"fits.tsv")); settings=input.settings
    @testset "Complete saved pilot and input identity" begin
        @test length(rows)==length(unique(getindex.(rows,"id")))==28length(settings.start_seeds)
        @test length(D.table(joinpath(out,"input_hashes.tsv")))==14
        for h in D.table(joinpath(out,"input_hashes.tsv"))
            k=h["input"]
            file=startswith(k,"--") ? opts[k] : joinpath(startswith(k,"native/") ? native : reference,split(k,'/')[2])
            @test bytes2hex(sha256(read(file)))==h["sha256"]
        end
        for i in 1:4
            shard=joinpath(out,"chunk$i")
            @test !isfile(joinpath(shard,"failures.tsv"))
            @test length(D.table(joinpath(shard,"fits.tsv")))==7length(settings.start_seeds)
            @test read(joinpath(shard,"input_hashes.tsv"))==read(joinpath(out,"input_hashes.tsv"))
        end
        summary=D.summaries(rows,settings)
        @test summary.winners==D.table(joinpath(out,"selected_starts.tsv"))
        @test summary.pilots==D.table(joinpath(out,"pilot.tsv"))
        for w in summary.winners
            group=filter(r->all(r[k]==w[k] for k in ("profile","family","mode","fold")),rows)
            @test val(w,"rss")==minimum(val.(group,"rss"))
        end
    end
    @testset "Independent predictions, all pixels, constraints, GCV and heldout RSS" begin
        for (i,case) in enumerate(D.CASES)
            ctx=D.D.context(input,case); d=ctx.data; k=length(ctx.initial.params)
            group=filter(r->r["profile"]==case.profile && r["family"]==case.family,rows)
            for r in group
                mode,fold,seed=r["mode"],parse(Int,r["fold"]),parse(Int,r["seed"])
                prefix=joinpath(out,"chunk$i",r["id"])
                @test only(D.table(prefix*".fit.tsv"))==r
                ps=sort(D.table(prefix*".parameters.tsv");by=s->parse(Int,s["parameter"]))
                q=val.(ps,"value"); q0=val.(ps,"initial"); lo=val.(ps,"lower"); hi=val.(ps,"upper")
                ix=fold==0 ? collect(eachindex(d.z)) : input.folds[fold].train
                test=fold==0 ? Int[] : input.folds[fold].test
                p=H.problem(ctx,mode,input.model_options,settings,input.step;indices=ix)
                @test lo==p.lo && hi==p.hi && all(lo.<=q.<=hi)
                @test q0==H.start(p,seed,settings.start_unit_box_radius)
                @test q0[1:3]==ctx.initial.params[1:3]
                @test mode=="fused" || all(iszero,q0[k+1:end])
                @test parse(Int,r["N"])==ctx.initial.n
                @test parse(Int,r["np"])==length(q)==k+(mode=="fused" ? 0 : mode=="paired" ? 4 : 6)
                pred=predict(q,ctx,mode,input.step)
                target=mode=="fused" ? d.z[ix] : vcat(ctx.fwd[ix],ctx.bwd[ix])
                objective=v->begin a=predict(v,ctx,mode,input.step); mode=="fused" ? a.mean[ix] : vcat(a.fwd[ix],a.bwd[ix]) end
                rss=sum(abs2,objective(q).-target); scale=max(sum(abs2,objective(p.q0).-target),input.options.objective_scale_floor)
                @test eq(scale,val(r,"objective_scale"))
                @test eq(sum(abs2,objective(q0).-target),val(r,"initial_rss"))
                @test eq(rss,val(r,"rss")) && parse(Int,r["nd"])==length(target)
                @test eq(length(target)/(length(target)-length(q))^2*rss,val(r,"gcv"))
                @test eq(sum(abs2,d.z.-pred.mean),val(r,"mean_rss"))
                @test eq(sum(abs2,ctx.fwd.-pred.fwd),val(r,"rss_fwd"))
                @test eq(sum(abs2,ctx.bwd.-pred.bwd),val(r,"rss_bwd"))
                peaks=[maximum(abs,d.z.-pred.mean),maximum(abs,ctx.fwd.-pred.fwd),maximum(abs,ctx.bwd.-pred.bwd)]./d.noise
                @test all(eq.(peaks,[val(r,k) for k in ("mean_peak_snr","peak_fwd","peak_bwd")]))
                @test val(r,"noise")==d.noise && val(r,"threshold")==ctx.cfg.residual_peak_snr_threshold
                # Geometry checks are shared; residual arithmetic and all acquisition terms are independent.
                nr=G.ChainModelResult(n=ctx.initial.n,params=q[1:k],success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
                G._chain_metrics!(nr,d.axisctx,ctx.cfg)
                expected=nr.overlap<=ctx.cfg.max_overlap && nr.endpoint_overrun_nm<=1e-6 && peaks[1]<=ctx.cfg.residual_peak_snr_threshold
                @test (r["mean_valid"]=="true")==expected
                @test (r["valid"]=="true")== (expected && (mode=="fused" || maximum(peaks[2:3])<=ctx.cfg.residual_peak_snr_threshold))
                if fold==0
                    @test isnan(val(r,"heldout_rss")) && r["heldout_n"]=="0"
                else
                    @test isempty(intersect(ix,test))
                    @test eq(sum(abs2,ctx.fwd[test].-pred.fwd[test])+sum(abs2,ctx.bwd[test].-pred.bwd[test]),val(r,"heldout_rss"))
                    @test parse(Int,r["heldout_n"])==2length(test)
                end
                saved=D.table(prefix*".residuals.tsv")
                @test length(saved)==length(input.pixels)
                @test getindex.(saved,"row")==getindex.(input.pixels,"row") && getindex.(saved,"column")==getindex.(input.pixels,"column")
                @test findall(==("true"),getindex.(saved,"training"))==ix
                @test findall(==("true"),getindex.(saved,"heldout"))==test
                for (name,values) in (("mean_prediction",pred.mean),("fwd_prediction",pred.fwd),("bwd_prediction",pred.bwd),
                    ("mean_residual",d.z.-pred.mean),("fwd_residual",ctx.fwd.-pred.fwd),("bwd_residual",ctx.bwd.-pred.bwd))
                    @test all(eq.(val.(saved,name),values))
                end
                trace=D.table(prefix*".trace.tsv")
                @test length(trace)==parse(Int,r["evaluations"])<=settings.slsqp_maxeval
                @test all(diff(val.(trace,"best_rss")).<=0)
                @test eq(val(first(trace),"rss"),val(r,"initial_rss"))
                @test parse(Int,r["gradient_evaluations"])==count(==("true"),getindex.(trace,"gradient_requested"))
                # Independent peak arithmetic also supplies the bound-gradient check.
                J1=H.C.finite_jacobian(objective,q,lo,hi,input.options.finite_difference_relative_step)
                J=H.C.finite_jacobian(objective,q,lo,hi,input.options.finite_difference_relative_step*.5)
                grad1=2 .* (hi.-lo).*(J1'*(objective(q).-target))./scale
                grad=2 .* (hi.-lo).*(J'*(objective(q).-target))./scale
                u=(q.-lo)./(hi.-lo); pg=u.-clamp.(u.-grad,0,1)
                @test isapprox(grad,val.(ps,"gradient");rtol=1e-5,atol=1e-7)
                @test isapprox(pg,val.(ps,"projected");rtol=1e-5,atol=1e-7)
                @test isapprox(norm(pg,Inf),val(r,"stationarity_half_step_projected_gradient");rtol=1e-5,atol=1e-7)
                tolerance=input.options.gradient_agreement_atol+input.options.gradient_agreement_rtol*max(norm(grad1,Inf),norm(grad,Inf))
                agreement=norm(grad1.-grad,Inf)<=tolerance
                pg1=u.-clamp.(u.-grad1,0,1)
                @test (r["stationarity_agreement"]=="true")==agreement
                @test (r["stationarity_passed"]=="true")== (agreement && max(norm(pg,Inf),norm(pg1,Inf))<=input.options.stationarity_tolerance)
            end
        end
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==3 || error("Usage: NATIVE_DIR REFERENCE_DIR RUN_DIR")
    VerifyPairedShift.verify(ARGS...)
end

#!/usr/bin/env julia
# Read-only replay: no optimizer, raw image, classifier or benchmark label.
module VerifyPairedConvergence
using Test, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_convergence.jl"))
module Previous
include(joinpath(@__DIR__,"verify_paired_acquisition_diagnostic.jl"))
end
const D=PairedConvergenceDiagnostic
const C=D.C
const ROOT=dirname(@__DIR__)
val(r,k)=parse(Float64,r[k])
eq(a,b)=isapprox(a,b;rtol=1e-9,atol=1e-12)
byparameter(rows)=sort(rows;by=r->parse(Int,r["parameter"]))

function verify(native,reference,out)
    opts=Dict("--native-dir"=>native,"--reference-dir"=>reference,
        "--physical-config"=>joinpath(ROOT,"config","chitosan.toml"),
        "--assignment-config"=>joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
        "--model-settings"=>joinpath(out,"model_settings.toml"),"--settings"=>joinpath(out,"settings.toml"))
    input=D.inputs(opts); options=input.options
    fits=D.table(joinpath(out,"fits.tsv")); hashes=D.table(joinpath(out,"input_hashes.tsv"))
    @testset "Saved convergence inputs, all cases and unchanged selection" begin
        @test length(fits)==16
        @test Set((r["profile"],r["family"],r["mode"],r["stage"]) for r in fits)==
            Set((c.profile,c.family,c.mode,s) for c in D.CASES for s in ("control","extended"))
        @test length(hashes)==12
        for r in hashes
            key=r["input"]
            path=startswith(key,"--") ? opts[key] :
                startswith(key,"native/") ? joinpath(native,split(key,'/')[2]) : joinpath(reference,split(key,'/')[2])
            @test bytes2hex(sha256(read(path)))==r["sha256"]
        end
        for i in 1:4
            dir=joinpath(out,"chunk$i")
            @test !isfile(joinpath(dir,"failures.tsv"))
            @test length(D.table(joinpath(dir,"fits.tsv")))==4
            for f in ("input_hashes.tsv","settings.toml","model_settings.toml")
                @test read(joinpath(dir,f))==read(joinpath(out,f))
            end
        end
        eligibility=D.table(joinpath(out,"eligibility.tsv"))
        @test length(eligibility)==8
        for r in eligibility
            valid=sort(filter(f->all(f[k]==r[k] for k in ("stage","mode","profile")) && f["valid"]=="true",fits);by=f->val(f,"gcv"))
            @test parse(Int,r["valid_families"])==length(valid)
            @test r["selected_family"]==(isempty(valid) ? "none" : first(valid)["family"])
            @test r["selected_lm_converged"]==(isempty(valid) ? "NA" : first(valid)["converged"])
            @test r["selected_stationary"]==(isempty(valid) ? "NA" : first(valid)["stationarity_passed"])
        end
    end
    @testset "Independent pixels, objective, checkpoint traces and final gradients" begin
        for (i,case) in enumerate(D.CASES)
            ctx=D.context(input,case); n=ctx.initial.n; k=length(ctx.initial.params)
            dir=joinpath(out,"chunk$(mod1(i,4))"); prefix="$(case.profile).$(case.family).$(case.mode)"
            pars=D.table(joinpath(dir,"$prefix.parameters.tsv")); grads=D.table(joinpath(dir,"$prefix.gradients.tsv"))
            traces=D.table(joinpath(dir,"$prefix.trace.tsv")); checkpoints=D.table(joinpath(dir,"$prefix.checkpoints.tsv"))
            cr=filter(r->r["stage"]=="control",traces); er=filter(r->r["stage"]=="extended",traces)
            @test length(er)>=length(cr)
            for (a,b) in zip(cr,er),key in ("iteration","rss","reported_gradient","trial_step_norm","lambda","accepted")
                @test a[key]==b[key]
            end
            predict=function(q)
                p=Previous.independently_predict(q[1:k],n,ctx.data.axisctx,ctx.cfg,
                    ctx.initial.amp_min,ctx.initial.amp_range,input.pixels)
                delta=case.mode=="paired" ? q[k+1].*p.lobes.+q[k+2].+q[k+3].*p.xs.+q[k+4].*p.ys : zeros(length(p.mean))
                return (;mean=p.mean,fwd=p.mean.+delta,bwd=p.mean.-delta)
            end
            objective=q->case.mode=="paired" ? vcat(predict(q).fwd,predict(q).bwd) : predict(q).mean
            target=case.mode=="paired" ? vcat(ctx.fwd,ctx.bwd) : ctx.data.z
            for stage in ("control","extended")
                r=only(f for f in fits if all(f[key]==getproperty(case,Symbol(key)) for key in ("profile","family","mode")) && f["stage"]==stage)
                @test r==only(f for f in D.table(joinpath(dir,"$prefix.fits.tsv")) if f["stage"]==stage)
                ps=byparameter(filter(p->p["stage"]==stage,pars)); gs=byparameter(filter(g->g["stage"]==stage,grads))
                q=val.(ps,"value"); initial=val.(ps,"initial"); lo=val.(ps,"lower"); hi=val.(ps,"upper")
                native_lo,native_hi=D.P.VP.native_raw_bounds(n,ctx.cfg)
                @test lo[1:k]==native_lo && hi[1:k]==native_hi
                @test initial[1:k]==ctx.initial.params
                if case.mode=="paired"
                    @test initial[k+1:end]==zeros(4)
                    @test hi[k+1:end]==[input.model_options.gain_delta_max,input.model_options.background_delta_max_nm,
                        input.model_options.tilt_delta_max,input.model_options.tilt_delta_max]
                    @test lo[k+1:end]==-hi[k+1:end]
                end
                @test all(lo.<=q.<=hi)
                @test parse(Int,r["N"])==n
                @test parse(Int,r["np"])==length(q)==k+(case.mode=="paired" ? 4 : 0)
                @test parse(Int,r["nd"])==length(target)
                pred=predict(q); residual=objective(q).-target; initial_rss=sum(abs2,objective(initial).-target)
                @test eq(initial_rss,val(r,"initial_rss"))
                @test eq(sum(abs2,residual),val(r,"rss"))
                @test eq(length(target)/(length(target)-length(q))^2*sum(abs2,residual),val(r,"gcv"))
                @test eq(sum(abs2,ctx.data.z.-pred.mean),val(r,"mean_rss"))
                @test eq(sum(abs2,ctx.fwd.-pred.fwd),val(r,"rss_fwd"))
                @test eq(sum(abs2,ctx.bwd.-pred.bwd),val(r,"rss_bwd"))
                @test eq(maximum(abs,ctx.data.z.-pred.mean)/ctx.data.noise,val(r,"mean_peak_snr"))
                @test eq(maximum(abs,ctx.fwd.-pred.fwd)/ctx.data.noise,val(r,"peak_fwd"))
                @test eq(maximum(abs,ctx.bwd.-pred.bwd)/ctx.data.noise,val(r,"peak_bwd"))
                @test val(r,"threshold")==ctx.cfg.residual_peak_snr_threshold==3.5
                native_result=D.G.ChainModelResult(n=n,params=q[1:k],success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
                D.P.VP.finalize_native!(native_result,ctx.data,ctx.cfg)
                @test native_result.valid==(r["native_valid"]=="true")
                @test (r["valid"]=="true")== (native_result.valid && (case.mode=="fused" || max(val(r,"peak_fwd"),val(r,"peak_bwd"))<=3.5))
                sr=D.table(joinpath(dir,"$prefix.$stage.residuals.tsv"))
                @test length(sr)==length(input.pixels)
                pixel_columns=(("mean_prediction",pred.mean),("fwd_prediction",pred.fwd),("bwd_prediction",pred.bwd),
                    ("mean_residual",ctx.data.z.-pred.mean),("fwd_residual",ctx.fwd.-pred.fwd),("bwd_residual",ctx.bwd.-pred.bwd))
                for j in eachindex(sr)
                    @test (sr[j]["row"],sr[j]["column"])==(input.pixels[j]["row"],input.pixels[j]["column"])
                    for (name,values) in pixel_columns
                        @test eq(val(sr[j],name),values[j])
                    end
                end
                tr=stage=="control" ? cr : er; iters=parse(Int,r["iterations"])
                @test parse.(Int,getindex.(tr,"iteration"))==collect(0:iters)
                @test all(diff(val.(tr,"rss")).<=0)
                @test eq(val(first(tr),"rss"),initial_rss) && eq(val(last(tr),"rss"),val(r,"rss"))
                @test 0<=iters<=(stage=="control" ? options.control_maxiter : options.extended_maxiter)
                @test (r["converged"]=="true")== (r["native_convergence_reason"]=="native_small_step_or_gradient")
                cp=filter(c->c["stage"]==stage,checkpoints)
                iterations=sort(unique(parse(Int,c["iteration"]) for c in cp))
                @test iterations==sort(unique(vcat(filter(i->i<=iters,options.checkpoints),[iters])))
                for iteration in iterations
                    point=val.(byparameter(filter(c->parse(Int,c["iteration"])==iteration,cp)),"value")
                    @test all(lo.<=point.<=hi)
                    @test eq(sum(abs2,objective(point).-target),val(tr[iteration+1],"rss"))
                    iteration==0 && (@test point==initial)
                    iteration==iters && (@test all(eq.(point,q)))
                end
                J1=C.finite_jacobian(objective,q,lo,hi,options.finite_difference_relative_step)
                J2=C.finite_jacobian(objective,q,lo,hi,options.finite_difference_relative_step*options.finite_difference_half_step_factor)
                scale=max(initial_rss,options.objective_scale_floor); width=hi.-lo; u=(q.-lo)./width
                g1=2 .* width .* (J1'*residual)./scale; g2=2 .* width .* (J2'*residual)./scale
                pg1=u.-clamp.(u.-g1,0,1); pg2=u.-clamp.(u.-g2,0,1)
                @test isapprox(g1,val.(gs,"gradient_step");rtol=1e-6,atol=1e-7)
                @test isapprox(g2,val.(gs,"gradient_half_step");rtol=1e-6,atol=1e-7)
                @test isapprox(pg2,val.(gs,"projected_half_step");rtol=1e-6,atol=1e-7)
                @test isapprox(norm(pg1,Inf),val(r,"stationarity_projected_gradient");rtol=1e-6,atol=1e-7)
                @test isapprox(norm(pg2,Inf),val(r,"stationarity_half_step_projected_gradient");rtol=1e-6,atol=1e-7)
                agreement=norm(g1.-g2,Inf)<=options.gradient_agreement_atol+options.gradient_agreement_rtol*max(norm(g1,Inf),norm(g2,Inf))
                @test agreement==(r["stationarity_agreement"]=="true")
                @test (agreement && max(norm(pg1,Inf),norm(pg2,Inf))<=options.stationarity_tolerance)==(r["stationarity_passed"]=="true")
                s=svdvals(J2.*reshape(width,1,:)); saved_s=val.(D.table(joinpath(dir,"$prefix.$stage.singular_values.tsv")),"value")
                @test isapprox(s,saved_s;rtol=1e-6,atol=1e-7)
                @test count(>(options.jacobian_rank_rtol*first(s)),s)==parse(Int,r["stationarity_rank"])
            end
        end
    end
end
function main(args=ARGS)
    length(args)==3 || error("Usage: verify_paired_convergence.jl NATIVE_DIR REFERENCE_DIR RUN_DIR")
    BLAS.set_num_threads(1)
    verify(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyPairedConvergence.main()

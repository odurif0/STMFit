#!/usr/bin/env julia
module VerifyPairedSolver
using Test, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_paired_solver.jl"))
module Independent
include(joinpath(@__DIR__,"verify_paired_acquisition_diagnostic.jl"))
end
const D=PairedSolverDiagnostic
const S=D.S
const C=S.C
const ROOT=dirname(@__DIR__)
val(r,k)=parse(Float64,r[k])
eq(a,b)=isapprox(a,b;rtol=1e-9,atol=1e-12)
sortpars(rows)=sort(rows;by=r->parse(Int,r["parameter"]))

function verify(native,reference,comparison,out)
    opts=Dict("--native-dir"=>native,"--reference-dir"=>reference,"--comparison-dir"=>comparison,
        "--physical-config"=>joinpath(ROOT,"config","chitosan.toml"),
        "--assignment-config"=>joinpath(ROOT,"config","unit_assignment_patch_support.toml"),
        "--model-settings"=>joinpath(out,"model_settings.toml"),"--settings"=>joinpath(out,"settings.toml"),
        "--solver-settings"=>joinpath(out,"solver_settings.toml"))
    input=D.inputs(opts); paths=D.input_paths(opts); fits=D.table(joinpath(out,"fits.tsv")); audit=input.options
    @testset "Solver comparison inputs, complete endpoints and unchanged selection" begin
        @test length(fits)==16
        @test Set((r["profile"],r["family"],r["mode"],r["solver"]) for r in fits)==
            Set((c.profile,c.family,c.mode,s) for c in D.D.CASES for s in ("lm","slsqp"))
        hashes=D.table(joinpath(out,"input_hashes.tsv"))
        @test length(hashes)==length(paths)==25
        for r in hashes; @test bytes2hex(sha256(read(paths[r["input"]])))==r["sha256"]; end
        for i in 1:4
            dir=joinpath(out,"chunk$i")
            @test !isfile(joinpath(dir,"failures.tsv"))
            @test length(D.table(joinpath(dir,"fits.tsv")))==4
            for f in ("input_hashes.tsv","settings.toml","solver_settings.toml","model_settings.toml")
                @test read(joinpath(dir,f))==read(joinpath(out,f))
            end
        end
        eligibility=D.table(joinpath(out,"eligibility.tsv")); @test length(eligibility)==8
        for r in eligibility
            candidates=sort(filter(f->all(f[k]==r[k] for k in ("solver","mode","profile")) && f["valid"]=="true",fits);by=f->val(f,"gcv"))
            @test parse(Int,r["valid_families"])==length(candidates)
            @test r["selected_family"]==(isempty(candidates) ? "none" : first(candidates)["family"])
            @test r["selected_stationary"]==(isempty(candidates) ? "NA" : first(candidates)["stationarity_passed"])
        end
    end
    @testset "Independent returned endpoints, residuals, gradients and trace checkpoints" begin
        for (i,case) in enumerate(D.D.CASES)
            ctx=D.D.context(input,case); k=length(ctx.initial.params); n=ctx.initial.n
            problem=S.problem(ctx.initial,ctx.data,ctx.cfg,ctx.fwd,ctx.bwd,case.mode,input.model_options)
            predict=function(q)
                p=Independent.independently_predict(q[1:k],n,ctx.data.axisctx,ctx.cfg,ctx.initial.amp_min,ctx.initial.amp_range,input.pixels)
                delta=case.mode=="paired" ? q[k+1].*p.lobes.+q[k+2].+q[k+3].*p.xs.+q[k+4].*p.ys : zeros(length(p.mean))
                return (;mean=p.mean,fwd=p.mean.+delta,bwd=p.mean.-delta)
            end
            objective=q->begin p=predict(q); case.mode=="paired" ? vcat(p.fwd,p.bwd) : p.mean end
            target=problem.target
            for solver in ("lm","slsqp")
                dir=joinpath(out,"chunk$(mod1(i,4))"); prefix="$(case.profile).$(case.family).$(case.mode).$solver"
                r=only(f for f in fits if all(f[key]==getproperty(case,Symbol(key)) for key in ("profile","family","mode")) && f["solver"]==solver)
                @test r==only(D.table(joinpath(dir,"$prefix.fit.tsv")))
                pars=sortpars(D.table(joinpath(dir,"$prefix.parameters.tsv"))); grads=sortpars(D.table(joinpath(dir,"$prefix.gradients.tsv")))
                q=val.(pars,"value"); lo=val.(pars,"lower"); hi=val.(pars,"upper"); initial=val.(pars,"initial")
                @test lo==problem.lo && hi==problem.hi && initial==problem.q0
                @test all(lo.<=q.<=hi)
                @test parse(Int,r["N"])==n && parse(Int,r["np"])==length(q) && parse(Int,r["nd"])==length(target)
                pred=predict(q); residual=objective(q).-target; rss=sum(abs2,residual); initial_rss=sum(abs2,objective(initial).-target)
                @test eq(rss,val(r,"rss")) && eq(initial_rss,val(r,"initial_rss"))
                @test eq(length(target)/(length(target)-length(q))^2*rss,val(r,"gcv"))
                @test eq(sum(abs2,ctx.data.z.-pred.mean),val(r,"mean_rss"))
                for (name,obs,pr) in (("fwd",ctx.fwd,pred.fwd),("bwd",ctx.bwd,pred.bwd))
                    @test eq(sum(abs2,obs.-pr),val(r,"rss_"*name))
                    @test eq(maximum(abs,obs.-pr)/ctx.data.noise,val(r,"peak_"*name))
                end
                @test eq(maximum(abs,ctx.data.z.-pred.mean)/ctx.data.noise,val(r,"mean_peak_snr"))
                @test val(r,"threshold")==ctx.cfg.residual_peak_snr_threshold==3.5
                native_result=S.G.ChainModelResult(n=n,params=q[1:k],success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
                S.P.VP.finalize_native!(native_result,ctx.data,ctx.cfg)
                @test native_result.valid==(r["native_valid"]=="true")
                @test (native_result.valid && (case.mode=="fused" || max(val(r,"peak_fwd"),val(r,"peak_bwd"))<=3.5))==(r["valid"]=="true")
                rows=D.table(joinpath(dir,"$prefix.residuals.tsv")); @test length(rows)==length(input.pixels)
                columns=(("mean_prediction",pred.mean),("fwd_prediction",pred.fwd),("bwd_prediction",pred.bwd),
                    ("mean_residual",ctx.data.z.-pred.mean),("fwd_residual",ctx.fwd.-pred.fwd),("bwd_residual",ctx.bwd.-pred.bwd))
                for j in eachindex(rows)
                    @test (rows[j]["row"],rows[j]["column"])==(input.pixels[j]["row"],input.pixels[j]["column"])
                    for (key,values) in columns; @test eq(val(rows[j],key),values[j]); end
                end
                J1=C.finite_jacobian(objective,q,lo,hi,audit.finite_difference_relative_step)
                J2=C.finite_jacobian(objective,q,lo,hi,audit.finite_difference_relative_step*audit.finite_difference_half_step_factor)
                scale=max(initial_rss,audit.objective_scale_floor); width=hi.-lo; u=(q.-lo)./width
                g1=2 .* width .* (J1'*residual)./scale; g2=2 .* width .* (J2'*residual)./scale
                pg1=u.-clamp.(u.-g1,0,1); pg2=u.-clamp.(u.-g2,0,1)
                @test isapprox(g1,val.(grads,"gradient_step");rtol=1e-6,atol=1e-7)
                @test isapprox(g2,val.(grads,"gradient_half_step");rtol=1e-6,atol=1e-7)
                @test isapprox(pg2,val.(grads,"projected_half_step");rtol=1e-6,atol=1e-7)
                @test isapprox(norm(pg1,Inf),val(r,"stationarity_projected_gradient");rtol=1e-6,atol=1e-7)
                @test isapprox(norm(pg2,Inf),val(r,"stationarity_half_step_projected_gradient");rtol=1e-6,atol=1e-7)
                agree=norm(g1.-g2,Inf)<=audit.gradient_agreement_atol+audit.gradient_agreement_rtol*max(norm(g1,Inf),norm(g2,Inf))
                @test agree==(r["stationarity_agreement"]=="true")
                @test (agree && max(norm(pg1,Inf),norm(pg2,Inf))<=audit.stationarity_tolerance)==(r["stationarity_passed"]=="true")
                s=svdvals(J2.*reshape(width,1,:))
                @test isapprox(s,val.(D.table(joinpath(dir,"$prefix.singular_values.tsv")),"value");rtol=1e-6,atol=1e-7)
                @test count(>(audit.jacobian_rank_rtol*first(s)),s)==parse(Int,r["stationarity_rank"])
                tr=D.table(joinpath(dir,"$prefix.trace.tsv")); cp=D.table(joinpath(dir,"$prefix.checkpoints.tsv"))
                key=solver=="lm" ? "iteration" : "evaluation"; steps=parse(Int,r[solver=="lm" ? "iterations" : "evaluations"])
                @test parse.(Int,getindex.(tr,key))==collect((solver=="lm" ? 0 : 1):steps)
                if solver=="lm"
                    @test all(diff(val.(tr,"rss")).<=0)
                    @test eq(val(last(tr),"rss"),rss)
                    @test steps<=input.solver_options.lm_maxiter
                    @test (r["lm_converged"]=="true")== (r["solver_return"]=="native_small_step_or_gradient")
                else
                    @test r["lm_converged"]==r["iterations"]=="NA"
                    @test parse(Int,r["gradient_evaluations"])==sum(t["gradient_requested"]=="true" for t in tr)
                    @test parse(Int,r["model_evaluations"])>=steps
                    @test all(eq.(val.(tr,"best_rss"),accumulate(min,val.(tr,"rss"))))
                end
                indices=sort(unique(parse(Int,c["index"]) for c in cp))
                declared=solver=="lm" ? audit.checkpoints : vcat(0,input.solver_options.evaluation_checkpoints)
                @test indices==sort(unique(vcat(filter(i->i<=steps,declared),[steps])))
                for index in indices
                    point=val.(sortpars(filter(c->parse(Int,c["index"])==index,cp)),"value")
                    @test all(lo.<=point.<=hi)
                    if index==0; @test point==initial; continue; end
                    saved=only(t for t in tr if parse(Int,t[key])==index)
                    @test eq(sum(abs2,objective(point).-target),val(saved,"rss"))
                end
                oldp=sortpars(filter(p->p["profile"]==case.profile && p["family"]==case.family && p["mode"]==case.mode && p["stage"]=="extended",input.oldparams))
                @test eq(maximum(abs.(q.-val.(oldp,"value"))),val(r,"previous_lm_parameter_difference"))
            end
        end
    end
end
function main(args=ARGS)
    length(args)==4 || error("Usage: verify_paired_solver.jl NATIVE_DIR REFERENCE_DIR PREVIOUS_COMPARISON_DIR RUN_DIR")
    BLAS.set_num_threads(1); verify(args...)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && VerifyPairedSolver.main()

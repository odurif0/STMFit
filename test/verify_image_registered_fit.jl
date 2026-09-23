#!/usr/bin/env julia
module VerifyImageRegisteredFit
using Test, SHA, LinearAlgebra
include(joinpath(@__DIR__,"diagnose_image_registered_fit.jl"))
include(joinpath(@__DIR__,"verify_paired_shift.jl"))
const D=ImageRegisteredFit
const H=D.H
const G=H.G
const V=VerifyPairedShift
val(r,k)=parse(Float64,r[k])
eq(a,b)=isapprox(a,b;atol=1e-12,rtol=1e-9)

function predict(q,ctx,mode,step,delta)
    mode=="image_fixed" ? V.predict(vcat(q,collect(delta)),ctx,"shifted",step) : V.predict(q,ctx,mode,step)
end

function verify(native,reference,image,out)
    opts=Dict("--native-dir"=>native,"--reference-dir"=>reference,"--image-dir"=>image,"--outdir"=>out)
    for key in ("physical-config","assignment-config","model-settings","settings","solver-settings","fit-settings")
        opts["--"*key]=joinpath(out,key*".toml")
    end
    input=D.inputs(opts); rows=D.table(joinpath(out,"fits.tsv")); s=input.settings
    @testset "Complete calibrated-fit schedule, starts and input binding" begin
        @test length(rows)==length(unique(getindex.(rows,"id")))==36length(s.start_seeds)
        hashes=D.table(joinpath(out,"input_hashes.tsv")); @test length(hashes)==20
        for h in hashes
            key=h["input"]
            if startswith(key,"--")
                path=opts[key]
            else
                tag,file=split(key,'/';limit=2)
                path=joinpath(Dict("native"=>native,"reference"=>reference,"image"=>image)[tag],file)
            end
            @test bytes2hex(sha256(read(path)))==h["sha256"]
        end
        summary=D.summaries(rows,s,input.limits)
        @test summary.winners==D.table(joinpath(out,"selected_starts.tsv"))
        @test summary.pilots==D.table(joinpath(out,"pilot.tsv"))
        for w in summary.winners
            group=filter(r->all(r[k]==w[k] for k in ("profile","family","mode","fold")),rows)
            @test val(w,"rss")==minimum(val.(group,"rss"))
        end
    end
    @testset "Independent translated peaks, mean/view guards, heldout scores and GCV" begin
        for (i,case) in enumerate(D.CASES)
            dir=joinpath(out,"chunk$i"); ctx=D.D.D.context(input,case); d=ctx.data; k=length(ctx.initial.params)
            @test !isfile(joinpath(dir,"failures.tsv"))
            @test read(joinpath(dir,"input_hashes.tsv"))==read(joinpath(out,"input_hashes.tsv"))
            group=filter(r->r["profile"]==case.profile && r["family"]==case.family,rows)
            @test length(group)==9length(s.start_seeds)
            for r in group
                mode,fold,seed=r["mode"],parse(Int,r["fold"]),parse(Int,r["seed"])
                prefix=joinpath(dir,r["id"]); @test only(D.table(prefix*".fit.tsv"))==r
                ps=sort(D.table(prefix*".parameters.tsv");by=p->parse(Int,p["parameter"]))
                q,q0,lo,hi=[val.(ps,key) for key in ("value","initial","lower","upper")]
                image_row=input.image_rows[fold+1]; delta=(val(image_row,"dx"),val(image_row,"dy"))
                ix=fold==0 ? collect(eachindex(d.z)) : input.folds[fold].train
                held=fold==0 ? Int[] : input.folds[fold].test
                p=D.problem(ctx,mode,input.model_options,input.step,input.limits,delta;indices=ix)
                @test lo==p.lo && hi==p.hi && all(lo.<=q.<=hi)
                @test q0==H.start(p,seed,s.start_unit_box_radius)
                @test q0[1:3]==ctx.initial.params[1:3] && all(iszero,q0[k+1:end])
                @test parse(Int,r["N"])==ctx.initial.n
                fitted=k+(mode=="shifted" ? 6 : 4); np=k+(mode=="paired" ? 4 : 6)
                @test length(q)==parse(Int,r["fitted_np"])==fitted && parse(Int,r["np"])==np
                @test (val(r,"image_dx_px"),val(r,"image_dy_px"))==delta
                expected_shift=mode=="paired" ? (0.,0.) : mode=="image_fixed" ? delta : Tuple(q[end-1:end])
                @test (val(r,"shift_x_px"),val(r,"shift_y_px"))==expected_shift
                pred=predict(q,ctx,mode,input.step,delta)
                objective=v->begin z=predict(v,ctx,mode,input.step,delta); vcat(z.fwd[ix],z.bwd[ix]) end
                target=vcat(ctx.fwd[ix],ctx.bwd[ix]); rss=sum(abs2,objective(q).-target)
                zero=D.problem(ctx,"paired",input.model_options,input.step,input.limits,delta;indices=ix)
                zpred=predict(zero.q0,ctx,"paired",input.step,delta)
                scale=max(sum(abs2,vcat(zpred.fwd[ix],zpred.bwd[ix]).-target),input.options.objective_scale_floor)
                @test eq(scale,val(r,"objective_scale")) && eq(sum(abs2,objective(q0).-target),val(r,"initial_rss"))
                @test eq(rss,val(r,"rss")) && parse(Int,r["nd"])==length(target)
                @test eq(length(target)/(length(target)-np)^2*rss,val(r,"gcv"))
                @test eq(sum(abs2,d.z.-pred.mean),val(r,"mean_rss"))
                @test eq(sum(abs2,ctx.fwd.-pred.fwd),val(r,"rss_fwd")) && eq(sum(abs2,ctx.bwd.-pred.bwd),val(r,"rss_bwd"))
                peaks=[maximum(abs,d.z.-pred.mean),maximum(abs,ctx.fwd.-pred.fwd),maximum(abs,ctx.bwd.-pred.bwd)]./d.noise
                @test all(eq.(peaks,[val(r,key) for key in ("mean_peak_snr","peak_fwd","peak_bwd")]))
                @test val(r,"noise")==d.noise && val(r,"threshold")==ctx.cfg.residual_peak_snr_threshold
                nr=G.ChainModelResult(n=ctx.initial.n,params=q[1:k],success=true,amp_min=ctx.initial.amp_min,amp_range=ctx.initial.amp_range)
                G._chain_metrics!(nr,d.axisctx,ctx.cfg)
                validmean=nr.overlap<=ctx.cfg.max_overlap && nr.endpoint_overrun_nm<=1e-6 && peaks[1]<=ctx.cfg.residual_peak_snr_threshold
                @test (r["mean_valid"]=="true")==validmean
                @test (r["valid"]=="true")== (validmean && maximum(peaks[2:3])<=ctx.cfg.residual_peak_snr_threshold)
                if fold==0
                    @test isnan(val(r,"heldout_rss")) && r["heldout_n"]=="0"
                else
                    @test isempty(intersect(ix,held))
                    @test eq(sum(abs2,ctx.fwd[held].-pred.fwd[held])+sum(abs2,ctx.bwd[held].-pred.bwd[held]),val(r,"heldout_rss"))
                    @test parse(Int,r["heldout_n"])==2length(held)
                end
                saved=D.table(prefix*".residuals.tsv")
                @test length(saved)==length(input.pixels)
                @test getindex.(saved,"row")==getindex.(input.pixels,"row") && getindex.(saved,"column")==getindex.(input.pixels,"column")
                @test findall(==("true"),getindex.(saved,"training"))==ix
                @test findall(==("true"),getindex.(saved,"heldout"))==held
                for (key,v) in (("mean_prediction",pred.mean),("fwd_prediction",pred.fwd),("bwd_prediction",pred.bwd),
                    ("mean_residual",d.z.-pred.mean),("fwd_residual",ctx.fwd.-pred.fwd),("bwd_residual",ctx.bwd.-pred.bwd))
                    @test all(eq.(val.(saved,key),v))
                end
                trace=D.table(prefix*".trace.tsv")
                @test length(trace)==parse(Int,r["evaluations"])<=s.slsqp_maxeval
                @test all(diff(val.(trace,"best_rss")).<=0) && eq(val(first(trace),"rss"),val(r,"initial_rss"))
                @test parse(Int,r["gradient_evaluations"])==count(==("true"),getindex.(trace,"gradient_requested"))
                J1=H.C.finite_jacobian(objective,q,lo,hi,input.options.finite_difference_relative_step)
                J=H.C.finite_jacobian(objective,q,lo,hi,input.options.finite_difference_relative_step*.5)
                grad1=2 .* (hi.-lo).*(J1'*(objective(q).-target))./scale
                grad=2 .* (hi.-lo).*(J'*(objective(q).-target))./scale
                u=(q.-lo)./(hi.-lo); pg=u.-clamp.(u.-grad,0,1); pg1=u.-clamp.(u.-grad1,0,1)
                @test isapprox(grad,val.(ps,"gradient");rtol=1e-5,atol=1e-7)
                @test isapprox(pg,val.(ps,"projected");rtol=1e-5,atol=1e-7)
                @test isapprox(norm(pg,Inf),val(r,"stationarity_half_step_projected_gradient");rtol=1e-5,atol=1e-7)
                tol=input.options.gradient_agreement_atol+input.options.gradient_agreement_rtol*max(norm(grad1,Inf),norm(grad,Inf))
                agreement=norm(grad1.-grad,Inf)<=tol
                @test (r["stationarity_agreement"]=="true")==agreement
                @test (r["stationarity_passed"]=="true")== (agreement && max(norm(pg,Inf),norm(pg1,Inf))<=input.options.stationarity_tolerance)
            end
        end
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==4 || error("Usage: NATIVE_DIR REFERENCE_DIR IMAGE_RUN_DIR FIT_RUN_DIR")
    VerifyImageRegisteredFit.verify(ARGS...)
end

#!/usr/bin/env julia
# Saved-output checks only: no optimizer, raw-image read, labels or grade.
# Usage: julia --project=. test/verify_paired_acquisition_diagnostic.jl NATIVE_DIR [PAIRED_DIR]
using Test, TOML, Statistics, LinearAlgebra, SHA
include(joinpath(@__DIR__,"diagnose_registered_fit_failure.jl"))
const D=RegisteredFailureDiagnostic
const G=D.G
const ROOT=dirname(@__DIR__)
val(r,k)=parse(Float64,r[k])
table(path)=last(D.read_table(path))
eq(a,b)=isapprox(a,b;rtol=1e-10,atol=1e-13)

function config(profile,family)
    raw=TOML.parsefile(joinpath(ROOT,"config","chitosan.toml"))
    assignment=TOML.parsefile(joinpath(ROOT,"config","unit_assignment_patch_support.toml"))
    raw["model"]["peak_profile"]=profile
    profile=="split" && (raw["model"]["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
    _,ec,cc=D.F.Extractor._configs(raw["model"],raw["preprocessing"],"unused")
    return family=="circ" ? cc : ec
end

function axis(row)
    ax,ay=val(row,"axis_x"),val(row,"axis_y")
    return (origin=(val(row,"origin_x_nm"),val(row,"origin_y_nm")),axis=[ax,ay],perp=[-ay,ax],
        tmin=val(row,"support_tmin"),tmax=val(row,"support_tmax"))
end

"Reconstruct pixels without the native prediction or paired prediction routines."
function independently_predict(p,n,a,cfg,amp_min,amp_range,pixels)
    _,features,_,_,_,_=G._decode_chain(p,n,a,cfg;amp_min,amp_range)
    xs=val.(pixels,"x_nm"); ys=val.(pixels,"y_nm")
    background=p[1].+p[2].*xs.+p[3].*ys
    lobes=zeros(length(pixels))
    for f in features, j in eachindex(lobes)
        dt=(xs[j]-f.x_nm)*a.axis[1]+(ys[j]-f.y_nm)*a.axis[2]
        du=-(xs[j]-f.x_nm)*a.axis[2]+(ys[j]-f.y_nm)*a.axis[1]
        width=cfg.peak_profile==:gaussian ? f.sigma_x_nm :
            dt<0 ? f.sigma_x_nm/sqrt(f.skew_ratio) : f.sigma_x_nm*sqrt(f.skew_ratio)
        lobes[j]+=f.amplitude*exp(-0.5*((dt/width)^2+(du/f.sigma_y_nm)^2))
    end
    return (;mean=background.+lobes,lobes,xs,ys)
end

function verify_native(dir)
    fits=table(joinpath(dir,"fits.tsv")); params=table(joinpath(dir,"parameters.tsv")); optimizers=table(joinpath(dir,"optimizers.tsv"))
    @testset "Native failure evidence independently reconstructed" begin
        @test length(fits)==length(optimizers)==8
        @test length(Set((r["arm"],r["profile"],r["family"]) for r in fits))==8
        @test length(unique(r["file"] for r in fits))==1
        for r in fits
            arm,profile,family=r["arm"],r["profile"],r["family"]
            cfg=config(profile,family); n=parse(Int,r["N"]); a=axis(r)
            sub=sort(filter(p->p["arm"]==arm && p["profile"]==profile && p["family"]==family && p["stage"]=="final",params);by=p->parse(Int,p["parameter"]))
            @test all(p["start"]=="1" for p in sub)
            p=val.(sub,"value"); lo,hi=D.VP.native_raw_bounds(n,cfg)
            @test length(p)==length(lo)==G._chain_nparams(n,cfg)
            @test all(lo.<=p.<=hi)
            pixels=table(joinpath(dir,"$arm.pixels.tsv"))
            @test length(pixels)==parse(Int,r["n_pixels"])
            pred=independently_predict(p,n,a,cfg,val(r,"amp_min"),val(r,"amp_range"),pixels)
            z=val.(pixels,"z_nm"); noise=val(r,"noise"); resid=z.-pred.mean
            @test all(eq.(z,(val.(pixels,"fwd_nm").+val.(pixels,"bwd_nm"))./2 .- val(r,"offset")))
            @test eq(sum(abs2,resid),val(r,"rss"))
            @test eq(length(z)/(length(z)-length(p))^2*sum(abs2,resid),val(r,"gcv"))
            @test eq(maximum(abs,resid)/noise,val(r,"residual_peak_snr"))
            @test eq(quantile(abs.(resid)./noise,.99),val(r,"p99_abs_residual_snr"))
            @test count(>(cfg.residual_peak_snr_threshold*noise),abs.(resid))==parse(Int,r["exceeding_pixels"])
            @test val(r,"threshold")==cfg.residual_peak_snr_threshold==3.5
            savedres=table(joinpath(dir,"$arm.$profile.$family.residuals.tsv"))
            @test length(savedres)==length(pixels)
            for j in eachindex(pixels)
                @test (savedres[j]["row"],savedres[j]["column"])==(pixels[j]["row"],pixels[j]["column"])
                @test eq(val(savedres[j],"prediction_nm"),pred.mean[j])
                @test eq(val(savedres[j],"residual_nm"),resid[j])
                @test eq(val(savedres[j],"residual_snr"),resid[j]/noise)
            end
            opt=only(q for q in optimizers if q["arm"]==arm && q["profile"]==profile && q["family"]==family)
            @test isempty(opt["global_error"]) && isempty(opt["lm_error"])
            @test 0<=parse(Int,opt["lm_iterations"])<=(family=="circ" ? cfg.max_iter : 50)
            @test (opt["lm_converged"]=="true")== (opt["lm_status"]=="converged")
        end
    end
    return fits
end

function verify_paired(dir,native)
    nf=table(joinpath(native,"fits.tsv")); fits=table(joinpath(dir,"fits.tsv")); params=table(joinpath(dir,"parameters.tsv"))
    pixels=table(joinpath(native,"registered.pixels.tsv"))
    @testset "Paired and matched fused evidence independently reconstructed" begin
        @test length(fits)==8
        @test !isfile(joinpath(dir,"failures.tsv"))
        for hash in table(joinpath(dir,"native_hashes.tsv"))
            @test bytes2hex(sha256(read(joinpath(native,hash["input"]))))==hash["sha256"]
        end
        for r in fits
            profile,family,mode=r["profile"],r["family"],r["mode"]
            cfg=config(profile,family); n=parse(Int,r["N"]); k=G._chain_nparams(n,cfg)
            nr=only(q for q in nf if q["arm"]=="registered" && q["profile"]==profile && q["family"]==family)
            sub=sort(filter(p->p["profile"]==profile && p["family"]==family && p["mode"]==mode,params);by=p->parse(Int,p["parameter"]))
            q=val.(sub,"value"); p=q[1:k]
            @test all(val.(sub,"lower").<=q.<=val.(sub,"upper"))
            pred=independently_predict(p,n,axis(nr),cfg,val(r,"amp_min"),val(r,"amp_range"),pixels)
            delta=mode=="paired" ? q[k+1].*pred.lobes.+q[k+2].+q[k+3].*pred.xs.+q[k+4].*pred.ys : zeros(length(pixels))
            fp,bp=pred.mean.+delta,pred.mean.-delta
            f=val.(pixels,"fwd_nm").-val(nr,"offset"); b=val.(pixels,"bwd_nm").-val(nr,"offset")
            mr=val.(pixels,"z_nm").-pred.mean; fr=f.-fp; br=b.-bp
            nd=mode=="paired" ? 2length(pixels) : length(pixels)
            rss=mode=="paired" ? sum(abs2,fr)+sum(abs2,br) : sum(abs2,mr)
            @test parse(Int,r["np"])==length(q)==k+(mode=="paired" ? 4 : 0)
            @test parse(Int,r["nd"])==nd
            @test eq(rss,val(r,"rss"))
            @test eq(nd/(nd-length(q))^2*rss,val(r,"gcv"))
            @test eq(sum(abs2,mr),val(r,"mean_rss"))
            @test eq(maximum(abs,mr)/val(r,"noise"),val(r,"mean_peak_snr"))
            @test eq(maximum(abs,fr)/val(r,"noise"),val(r,"peak_fwd"))
            @test eq(maximum(abs,br)/val(r,"noise"),val(r,"peak_bwd"))
            @test val(r,"threshold")==cfg.residual_peak_snr_threshold==3.5
            @test (r["valid"]=="true")==(r["native_valid"]=="true" && (mode=="fused" || max(val(r,"peak_fwd"),val(r,"peak_bwd"))<=3.5))
            sr=table(joinpath(dir,"$profile.$family.$mode.residuals.tsv"))
            @test length(sr)==length(pixels)
            for j in eachindex(pixels)
                @test (sr[j]["row"],sr[j]["column"])==(pixels[j]["row"],pixels[j]["column"])
                for (key,values) in (("mean_prediction",pred.mean),("fwd_prediction",fp),("bwd_prediction",bp),
                        ("mean_residual",mr),("fwd_residual",fr),("bwd_residual",br))
                    @test eq(val(sr[j],key),values[j])
                end
            end
        end
        for r in table(joinpath(dir,"eligibility.tsv"))
            selected=sort(filter(q->q["profile"]==r["profile"] && q["mode"]==r["mode"] && q["valid"]=="true",fits);by=q->val(q,"gcv"))
            @test length(selected)==parse(Int,r["valid_families"])
            @test r["selected_family"]==(isempty(selected) ? "none" : first(selected)["family"])
        end
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (1,2) || error("Usage: NATIVE_DIR [PAIRED_DIR]")
    verify_native(ARGS[1])
    length(ARGS)==2 && verify_paired(ARGS[2],ARGS[1])
end

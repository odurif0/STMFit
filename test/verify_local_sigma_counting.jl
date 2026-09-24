#!/usr/bin/env julia
# Saved-output physics and selection checks; no optimization or benchmark input.
module LocalSigmaVerification
include(joinpath(@__DIR__,"diagnose_local_sigma_counting.jl"))
using Test, TOML, SHA, Statistics
const C=LocalSigmaCounting
const P=C.P
const G=C.G
num(r,k)=parse(Float64,r[k])
int(r,k)=parse(Int,r[k])
rows(p)=last(P.read_table(p))
key(r)=(r["file"],int(r,"repetition"),r["arm"],int(r,"N"),r["family"])
close(a,b)=isapprox(a,b;atol=1e-12,rtol=1e-10)
snapshot(cfg)=Dict(string(k)=>(getfield(cfg,k) isa Symbol ? string(getfield(cfg,k)) : getfield(cfg,k)) for k in fieldnames(typeof(cfg)))

function verify(root;rawdir=nothing)
    dirs=isfile(joinpath(root,"cohort.tsv")) ? [root] :
        sort(filter(p->isdir(p)&&startswith(basename(p),"chunk"),readdir(root;join=true)))
    isempty(dirs) && error("No counting shards")
    seen=Set{String}();selections=Dict{String,String}[];allfits=Dict{String,String}[]
    @testset "Complete counting, actual widths, saved forward model, full GCV" begin
        for dir in dirs
            @test !isfile(joinpath(dir,"failures.tsv"))
            raw=TOML.parsefile(joinpath(dir,"count.toml"));s=C.settings(joinpath(dir,"settings.toml"))
            hashes=Dict(r["input"]=>r["sha256"] for r in rows(joinpath(dir,"inputs.tsv")))
            @test C.sha(joinpath(dir,"count.toml"))==hashes["--count-config"]
            @test C.sha(joinpath(dir,"settings.toml"))==hashes["--settings"]
            cfgs=Dict(a=>C.configs(raw,s,a) for a in C.ARMS)
            for a in C.ARMS,f in ("circ","ell")
                @test isequal(TOML.parsefile(joinpath(dir,"$a.$f.toml")),snapshot(getproperty(cfgs[a],Symbol(f))))
            end
            cohort=rows(joinpath(dir,"cohort.tsv"));files=Set(r["file"] for r in cohort)
            @test length(files)==length(cohort) && isempty(intersect(files,seen))
            union!(seen,files)
            if rawdir!==nothing
                paths=P.raw_index(rawdir)
                @test all(haskey(paths,r["file"]) && C.sha(paths[r["file"]])==r["raw_sha256"] for r in cohort)
            end
            contexts=Dict(r["file"]=>r for r in rows(joinpath(dir,"context.tsv")))
            @test Set(keys(contexts))==files
            @test all(contexts[r["file"]]["raw_sha256"]==r["raw_sha256"] for r in cohort)
            fits=rows(joinpath(dir,"candidates.tsv"));lobes=rows(joinpath(dir,"lobes.tsv"))
            selected=rows(joinpath(dir,"selected.tsv"));diag=rows(joinpath(dir,"diagnostics.tsv"))
            @test length(unique(key.(fits)))==length(fits)
            @test length(unique((key(r),int(r,"lobe")) for r in lobes))==length(lobes)
            @test Set(r["file"] for r in fits)==files
            @test Set(key(r) for r in lobes)==Set(key(r) for r in fits if r["success"]=="true")
            @test Set(key(r) for r in diag) ⊆ Set(key.(fits))
            @test length(selected)==length(files)*2length(C.ARMS)
            @test length(unique((r["file"],r["repetition"],r["arm"]) for r in selected))==length(selected)
            grouped=Dict(k=>[r for r in lobes if key(r)==k] for k in unique(key.(lobes)))
            for file in sort(collect(files))
                ctx=contexts[file];axis=(origin=(num(ctx,"origin_x_nm"),num(ctx,"origin_y_nm")),
                    axis=[num(ctx,"axis_x"),num(ctx,"axis_y")],perp=[-num(ctx,"axis_y"),num(ctx,"axis_x")],
                    tmin=num(ctx,"support_tmin_nm"),tmax=num(ctx,"support_tmax_nm"))
                @test close(sum(abs2,axis.axis),1.)
                @test 0<int(ctx,"observed_pixels")<=int(ctx,"total_pixels")
                data=rows(joinpath(dir,"fit_data",file*".tsv"))
                x=[num(r,"x_nm") for r in data];y=[num(r,"y_nm") for r in data];z=[num(r,"z_nm") for r in data]
                @test length(z)==int(ctx,"n_data") && all(isfinite,vcat(x,y,z))
                for a in C.ARMS,rep in 1:2
                    cfgbase=cfgs[a].circ
                    lower=a=="global_sigma_max" ? max(cfgbase.spacing_min_nm,sqrt(-2log(cfgbase.max_overlap))*cfgbase.sigma_parallel_max_nm) :
                        max(cfgbase.spacing_min_nm,sqrt(-2log(cfgbase.max_overlap))*cfgbase.sigma_parallel_min_nm)
                    upperN=min(14,max(1,floor(Int,(axis.tmax-axis.tmin)/lower)+1))
                    subset=[r for r in fits if r["file"]==file && r["arm"]==a && int(r,"repetition")==rep]
                    @test Set((int(r,"N"),r["family"]) for r in subset)==Set((n,f) for n in 2:upperN for f in ("circ","ell"))
                    for r in subset
                        n=int(r,"N");cfg=getproperty(cfgs[a],Symbol(r["family"]))
                        @test int(r,"n_data")==length(z)
                        pfull=3+(r["family"]=="circ" ? 4 : 5)*n # tilted baseline + free gaps/centers/widths
                        @test int(r,"p_full")==pfull
                        if r["success"]!="true"
                            @test r["valid"]=="false" && isempty(r["parameters"]) && !isempty(r["reason"])
                            continue
                        end
                        pars=parse.(Float64,split(r["parameters"],';'))
                        @test length(pars)==pfull && all(isfinite,pars)
                        lo,hi=C.VP.native_raw_bounds(n,cfg)
                        @test all(lo.-1e-12 .<= pars .<= hi.+1e-12)
                        ls=sort(grouped[key(r)];by=q->int(q,"lobe"))
                        @test [int(q,"lobe") for q in ls]==collect(1:n)
                        ts=[num(q,"t_nm") for q in ls];us=[num(q,"u_nm") for q in ls]
                        sp=[num(q,"sigma_parallel_nm") for q in ls];sq=[num(q,"sigma_perp_nm") for q in ls]
                        @test all(lower-1e-12 .<= diff(ts) .<= cfg.spacing_max_nm+1e-12)
                        @test all(cfg.sigma_parallel_min_nm-1e-12 .<= sp .<= cfg.sigma_parallel_max_nm+1e-12)
                        @test all(cfg.sigma_perp_min_nm-1e-12 .<= sq .<= cfg.sigma_perp_max_nm+1e-12)
                        @test first(ts)>=axis.tmin-1e-12 && last(ts)<=axis.tmax+1e-12
                        _,decoded,td,ud,pd,qd=G._decode_chain(pars,n,axis,cfg;amp_min=num(r,"amp_min"),amp_range=num(r,"amp_range"))
                        @test td==ts && ud==us && pd==sp && qd==sq
                        pred=num(r,"baseline").+num(r,"tilt_x").*x.+num(r,"tilt_y").*y
                        for i in 1:n
                            q=ls[i];xc=num(q,"x_nm");yc=num(q,"y_nm")
                            @test close(xc,axis.origin[1]+ts[i]*axis.axis[1]+us[i]*axis.perp[1])
                            @test close(yc,axis.origin[2]+ts[i]*axis.axis[2]+us[i]*axis.perp[2])
                            @test decoded[i].amplitude==num(q,"amplitude")
                            cap=min(i>1 ? ts[i]-ts[i-1] : Inf,i<n ? ts[i+1]-ts[i] : Inf)/sqrt(-2log(cfg.max_overlap))
                            @test close(cap,num(q,"local_sigma_cap_nm"))
                            a=="local_sigma_cap" && @test max(sp[i],sq[i])<=cap+1e-12
                            dt=(x.-xc).*axis.axis[1].+(y.-yc).*axis.axis[2]
                            du=(x.-xc).*axis.perp[1].+(y.-yc).*axis.perp[2]
                            pred .+= num(q,"amplitude") .* exp.(-0.5 .* ((dt ./ sp[i]).^2 .+ (du ./ sq[i]).^2))
                        end
                        rss=sum(abs2,z.-pred);gcv=length(z)>pfull ? rss*length(z)/(length(z)-pfull)^2 : Inf
                        @test close(rss,num(r,"rss")) && close(gcv,num(r,"gcv"))
                        ov=maximum(exp(-0.5*(hypot(ts[i]-ts[j],us[i]-us[j])/max(sp[i],sp[j],sq[i],sq[j]))^2) for i in 1:n-1 for j in i+1:n)
                        meanov=maximum(exp(-0.5*(hypot(ts[i]-ts[j],us[i]-us[j])/max(mean(sp),mean(sq)))^2) for i in 1:n-1 for j in i+1:n)
                        @test close(ov,num(r,"pair_overlap")) && close(meanov,num(r,"mean_overlap"))
                        nativeov=a=="local_sigma_cap" ? ov : meanov
                        @test close(nativeov,num(r,"overlap")) && ov<=cfg.max_overlap+1e-12
                        residual=maximum(abs,z.-pred)/num(ctx,"noise")
                        @test close(residual,num(r,"residual_peak_snr"))
                        @test (r["valid"]=="true")== (nativeov<=cfg.max_overlap && num(r,"endpoint_overrun_nm")<=1e-6 && residual<=cfg.residual_peak_snr_threshold && length(z)>pfull)
                        @test any(key(q)==key(r) for q in diag)
                    end
                    valid=sort(filter(r->r["success"]==r["valid"]=="true" && isfinite(num(r,"gcv")),subset);
                        by=r->(num(r,"gcv"),int(r,"N"),r["family"]=="ell" ? 0 : 1))
                    selection=only(filter(r->r["file"]==file && r["arm"]==a && int(r,"repetition")==rep,selected))
                    @test !isempty(valid) && selection["status"]=="ok"
                    isempty(valid) && continue
                    @test int(selection,"N_selected")==int(first(valid),"N") && selection["source"]==first(valid)["family"]
                    @test selection["gcv"]==first(valid)["gcv"]
                end
            end
            append!(selections,selected);append!(allfits,fits)
        end
        rawdir===nothing || @test seen==Set(keys(P.raw_index(rawdir)))
    end
    println("Verified $(length(seen)) scans, $(length(allfits)) candidate families and $(length(selections)) selections; no grade.")
    return (;files=seen,fits=allfits,selected=selections)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (1,2) || error("verify_local_sigma_counting.jl RUN [RAW_DIR]")
    LocalSigmaVerification.verify(ARGS[1];rawdir=length(ARGS)==2 ? ARGS[2] : nothing)
end

#!/usr/bin/env julia
# Two bounded linear profiles, used ONLY as subtraction models. No new geometry.
module CrossViewResidual
include(joinpath(@__DIR__, "profile_frozen_amplitudes.jl"))
include(joinpath(@__DIR__, "lib", "patch_acquisition.jl"))
using .FrozenAmplitudeProfile: lobe_table, write_table, number, fmt
using LinearAlgebra, Statistics, TOML, Printf
const F = FrozenAmplitudeProfile
const G = F.G

function read_settings(path)
    raw=TOML.parsefile(path)
    Set(keys(raw))==Set(["model","selection","preprocessing","cross_view_residual","counting_variable_projection"]) || error("Unexpected settings section")
    all(isempty(raw[k]) for k in ("model","selection","preprocessing")) || error("Use the unchanged count config")
    raw["cross_view_residual"]==Dict("validity"=>"native_fused_mean") || error("Native fused-mean validity required")
    return F.VP.read_options(raw["counting_variable_projection"])
end

function parse_cli(args)
    if "--help" in args
        println("profile_cross_view_residuals.jl --features TSV --data-dir DIR --config COUNT.toml --settings SETTINGS.toml --out NEW.tsv [--chunk I/N] [--dry-run]")
        println("Primary output is the fwd-trained subtraction model; .bwd.tsv is bwd-trained. Geometry and main classifier features never change. Julia 1.13; multi-file fits on Viper only.")
        return nothing
    end
    opts=F.parse_cli(args)
    for suffix in (".bwd.tsv",)
        p=opts["--out"]*suffix
        (ispath(p)||islink(p)) && error("Output exists: $p")
    end
    read_settings(opts["--settings"])
    return opts
end

"Native fused guards evaluated on literal cached geometry and mean coefficients."
function mean_diagnostics(rows,x,y,z,noise,axis,cfg,profiles)
    A=F.frozen_design(x,y,rows,cfg.chain_tilted_baseline)
    beta=(profiles.fwd.result.x+profiles.bwd.result.x)/2
    pred=A*beta; residual=z-pred
    sigma=max(mean(number(r,"sigma_parallel_nm") for r in rows),mean(number(r,"sigma_perp_nm") for r in rows),G.EPS)
    overlap=0.0
    for i in 1:length(rows)-1, j in i+1:length(rows)
        d2=(number(rows[i],"x_nm")-number(rows[j],"x_nm"))^2+(number(rows[i],"y_nm")-number(rows[j],"y_nm"))^2
        overlap=max(overlap,exp(-0.5*d2/sigma^2))
    end
    ts=[number(r,"t_nm") for r in rows]
    overrun=G.endpoint_overrun(ts,axis.tmin,axis.tmax)
    pfull=G._chain_nparams(length(rows),cfg)
    cv,_=G._chain_gcv_score(z,pred,noise,pfull,cfg.student_nu)
    snr=maximum(abs,residual)/max(noise,G.EPS)
    rss=sum(abs2,residual); gcv=length(z)/(length(z)-pfull)^2*rss
    valid=isfinite(cv) && overlap<=cfg.max_overlap && overrun<=1e-6 && snr<=cfg.residual_peak_snr_threshold
    return (;valid,overlap,overrun,cv,snr,rss,gcv,pfull)
end

function profile_pair(rows,x,y,zf,zb,lower,upper,cfg,options)
    # Conditional on the fixed design and bounds, neither solve sees its target
    # opposite view. Geometry, support and bounds themselves use the fused data.
    fwd=F.profile_rows(rows,x,y,zf,lower,upper,cfg,options)
    bwd=F.profile_rows(rows,x,y,zb,lower,upper,cfg,options)
    return (;fwd,bwd)
end

function execute(opts)
    BLAS.set_num_threads(1)
    raw=TOML.parsefile(opts["--config"]); options=read_settings(opts["--settings"])
    pcfg,cfg,_=F.Extractor._configs(raw["model"],raw["preprocessing"],dirname(opts["--out"]))
    G._chain_peak_profile(cfg)==:gaussian && cfg.chain_tilted_baseline || error("Tilted Gaussian base required")
    startswith(get(raw["model"],"selection_policy",""),"adaptive_support") && error("Adaptive support is outside this experiment")
    header,table=lobe_table(opts["--features"];required=F.REQUIRED)
    any(startswith(c,"model_axis") || c=="model_orientation" || c=="profile_view" for c in header) && error("Only original global-axis cache accepted")
    files=sort(unique(first.(collect(keys(table)))))
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(files) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty chunk")
    chains=Dict(f=>[table[k] for k in sort(collect(keys(table))) if first(k)==f] for f in files)
    for f in files
        F.validate_chain(chains[f]); isfile(joinpath(opts["--data-dir"],f)) || error("Missing raw $f")
    end
    haskey(opts,"--dry-run") && return println("Cross-view profiles: $(length(files)) files; metadata only, no pixels or output.")
    fitdir=opts["--out"]*".fit_data"; mkpath(fitdir)
    outputs=Dict(v=>Dict{String,String}[] for v in ("fwd","bwd"))
    audit=Dict{String,String}[]; coefficients=Dict{String,String}[]; failures=Dict{String,String}[]
    for (idx,f) in enumerate(files)
        try
            rs=chains[f]; n=length(rs); pcfg.filepath=joinpath(opts["--data-dir"],f)
            img=G.read_sxm(pcfg.filepath)
            cf=G.get_channel(img,"Z";direction="fwd"); cb=G.get_channel(img,"Z";direction="bwd")
            PatchAcquisition.require_direction(cf,"fwd"); PatchAcquisition.require_direction(cb,"bwd")
            xs,ys,zimg,mask,x,y,z,noise=G._fused_roi_data(img,pcfg)
            full=G._weighted_roi_axis(x,y,z)
            xf,yf,zm,axis,keep,_=G._chain_fit_data(x,y,z,full,cfg)
            for (key,value,isaxis) in (("axis_x",axis.axis[1],true),("axis_y",axis.axis[2],true),
                ("origin_x_nm",axis.origin[1],false),("origin_y_nm",axis.origin[2],false))
                printed=isaxis ? @sprintf("%.8f",value) : @sprintf("%.6f",value)
                printed==rs[1][key] || error("Cached $key differs from input/support")
            end
            x1,y1,_,zf,_,_,_=G.preprocess_channel(img,cf,pcfg)
            x2,y2,_,zb,_,_,_=G.preprocess_channel(img,cb,pcfg)
            xs==x1==x2 && ys==y1==y2 || error("Acquisition grids differ")
            offset=quantile(((zf+zb)/2)[mask],0.05)
            xx,yy,zff=G._flatten_roi(zf.-offset,mask,xs,ys)
            xxx,yyy,zbb=G._flatten_roi(zb.-offset,mask,xs,ys)
            xx==xxx==x && yy==yyy==y || error("Different observed fit pixels")
            zff=zff[keep]; zbb=zbb[keep]
            isapprox((zff+zbb)/2,zm;atol=options.mapping_atol,rtol=options.mapping_rtol) || error("View mean differs from native fused data")
            localcfg=deepcopy(cfg); localcfg.chain_circular_sigmas=rs[1]["source"]=="circ"
            lower,upper=F.coefficient_bounds(n,localcfg,axis,zimg)
            pair=profile_pair(rs,xf,yf,zff,zbb,lower,upper,localcfg,options)
            diag=mean_diagnostics(rs,xf,yf,zm,noise,axis,localcfg,pair)
            A=F.frozen_design(xf,yf,rs,true)
            open(joinpath(fitdir,f*".tsv"),"w") do io
                println(io,"x_nm\ty_nm\tz_fwd_nm\tz_bwd_nm\tz_fused_nm")
                for j in eachindex(zm); println(io,join(fmt.((xf[j],yf[j],zff[j],zbb[j],zm[j])),'\t')); end
            end
            for (view,obs,other,r) in (("fwd",zff,zbb,pair.fwd),("bwd",zbb,zff,pair.bwd))
                for row in r.rows; row["profile_view"]=view; end
                append!(outputs[view],r.rows)
                names=vcat(["baseline","tilt_x","tilt_y"],["amplitude_$j" for j in 1:n])
                for j in eachindex(names)
                    push!(coefficients,Dict("file"=>f,"view"=>view,"parameter"=>names[j],"initial"=>fmt(r.initial[j]),
                        "value"=>fmt(r.result.x[j]),"lower"=>fmt(lower[j]),"upper"=>fmt(upper[j]),"gradient"=>fmt(r.result.gradient[j])))
                end
                pred=A*r.result.x
                values=(file=f,view=view,N=n,n_data=length(zm),p_full=r.pfull,initial_rss=r.initial_rss,rss=r.result.rss,
                    gcv=r.gcv,cross_rss=sum(abs2,other-pred),own_peak_snr=maximum(abs,obs-pred)/noise,
                    cross_peak_snr=maximum(abs,other-pred)/noise,kkt_violation=r.result.kkt_violation,kkt_tolerance=r.result.kkt_tolerance,
                    iterations=r.result.iterations,status=r.result.status,initial_clamp=r.initial_clamp,
                    active_lower=join(r.result.active_lower,','),active_upper=join(r.result.active_upper,','),
                    fused_valid=diag.valid,fused_rss=diag.rss,fused_gcv=diag.gcv,fused_cv=diag.cv,
                    fused_peak_snr=diag.snr,overlap=diag.overlap,endpoint_overrun_nm=diag.overrun,
                    support_tmin=axis.tmin,support_tmax=axis.tmax,noise=noise,offset=offset,
                    amp_max_data=max(maximum(zimg),G.EPS))
                push!(audit,Dict(string(k)=>string(v) for (k,v) in pairs(values)))
            end
            diag.valid || error("Mean model fails unchanged native fused validity: SNR=$(diag.snr), overlap=$(diag.overlap), overrun=$(diag.overrun)")
            println("[$idx/$(length(files))] $f N=$n two KKT solves; mean SNR=$(diag.snr)"); flush(stdout)
        catch err
            push!(failures,Dict("file"=>f,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')))
            println(stderr,"FAILED $f: ",sprint(showerror,err))
        end
    end
    !isempty(audit) && write_table(opts["--out"]*".audit.tsv",sort(collect(keys(first(audit)))),audit)
    !isempty(coefficients) && write_table(opts["--out"]*".coefficients.tsv",["file","view","parameter","initial","value","lower","upper","gradient"],coefficients)
    if !isempty(failures)
        write_table(opts["--out"]*".failures.tsv",["file","reason"],failures)
        error("$(length(failures)) cross-view profiles failed; no incomplete model table emitted")
    end
    write_table(opts["--out"],vcat(header,["profile_view"]),outputs["fwd"])
    write_table(opts["--out"]*".bwd.tsv",vcat(header,["profile_view"]),outputs["bwd"])
end
function main(args=ARGS)
    opts=parse_cli(args); opts===nothing || execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && CrossViewResidual.main()

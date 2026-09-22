#!/usr/bin/env julia
# Native fixed-N base AND split refits. No labels, counts search or new optimizer.
module RegisteredRefit
include(joinpath(@__DIR__, "refine_fixed_n_geometry.jl"))
include(joinpath(@__DIR__, "lib", "patch_acquisition.jl"))
using .FixedNGeometryRefinement: G, VP
using .FixedNGeometryRefinement.FrozenAmplitudeProfile.ReconstructedUnitAssignment: lobe_table, read_table, write_table
using .PatchAcquisition
using TOML, Statistics, LinearAlgebra, Printf
const R = FixedNGeometryRefinement
const F = R.F
const OPTIONS = Set(["--features", "--data-dir", "--config", "--assignment-config", "--settings", "--shifts", "--out"])
const ARMS = ("control", "registered")
const PROFILES = ("gaussian", "split")

function settings(path)
    s = TOML.parsefile(path)
    expected = Dict(
        "model" => Dict("method"=>"native_circular_to_elliptical", "profiles"=>["gaussian","split"],
            "native_elliptical_maxiter"=>50, "one_dimensional_initialization"=>false),
        "selection" => Dict("count_policy"=>"reuse_saved_N", "family_policy"=>"valid_minimum_full_parameter_gcv",
            "zero_shift_policy"=>"reuse_control_fit_exactly", "failure_policy"=>"report_and_stop_no_partial_grade"),
        "preprocessing" => Dict("fusion"=>"mean_unsmoothed_observed_views", "roi"=>"native_finite_statistics_and_observed_mask",
            "seed_missing"=>"nearest_observed_initialization_only", "reference_frame"=>"forward"))
    frozen = deepcopy(expected)
    frozen["preprocessing"]["roi"] = "original_native_roi_intersect_observed_mask"
    frozen["preprocessing"]["geometry"] = "freeze_original_axis_tube_and_bounds"
    s in (expected, frozen) || error("Only the declared registered-refit experiments are supported")
    return (native_elliptical_maxiter=Int(s["model"]["native_elliptical_maxiter"]),
        original_support=s==frozen)
end

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(opts,k) && error("Repeated option $k")
        if k=="--dry-run"; opts[k]="true"; i+=1; continue; end
        k in union(OPTIONS,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden refit option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All registered-refit inputs required")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    isdir(opts["--data-dir"]) || error("Raw directory missing")
    for k in setdiff(OPTIONS,Set(["--data-dir","--out"])); isfile(opts[k]) || error("Missing input $k"); end
    for suffix in vcat(["", ".fits.tsv", ".parameters.tsv", ".bootstrap.tsv", ".support.tsv", ".failures.tsv", ".fit_data"],
            [".$a.$p.tsv" for a in ARMS for p in PROFILES])
        (ispath(opts["--out"]*suffix) || islink(opts["--out"]*suffix)) && error("Refit output exists")
    end
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] || error("Invalid chunk")
    settings(opts["--settings"])
    return opts
end

"Load each native preprocessing once; no second backward flip or gain fitting."
function load_views(img, pcfg)
    cf=G.get_channel(img,"Z";direction="fwd"); cb=G.get_channel(img,"Z";direction="bwd")
    require_direction(cf,"fwd"); require_direction(cb,"bwd")
    xs,ys,_,zf,_,_,nf=G.preprocess_channel(img,cf,pcfg)
    xb,yb,_,zb,_,_,nb=G.preprocess_channel(img,cb,pcfg)
    xs==xb && ys==yb || error("Acquisition grids differ")
    return (;xs,ys,zf,zb,cf,cb,noise=max(nf,nb,G.EPS))
end

"Replay native unregistered ROI/axis/tube/bounds, before restoring observation masks. No fit."
function original_support(img, pcfg, cfg)
    xs,ys,zimg,mask,x,y,z,noise=G._fused_roi_data(img,pcfg)
    full=G._weighted_roi_axis(x,y,z)
    _,_,_,axis,keep,support=G._chain_fit_data(x,y,z,full,cfg)
    indices=[CartesianIndex(iy,ix) for iy in eachindex(ys) for ix in eachindex(xs) if mask[iy,ix]]
    length(indices)==length(keep)==length(z) || error("Original support indexing mismatch")
    fitmask=falses(size(mask)); fitmask[indices[keep]].=true
    return (;xs,ys,zimg,mask,x,y,z,axis,keep,support,indices,fitmask,noise)
end

"Check the replayed frame at the original cache's serialization precision, not a fitted tolerance."
function check_original_frame(original, row)
    a=original.axis
    for (key,value,precision) in (("axis_x",a.axis[1],8),("axis_y",a.axis[2],8),
            ("origin_x_nm",a.origin[1],6),("origin_y_nm",a.origin[2],6))
        formatted=precision==8 ? @sprintf("%.8f",value) : @sprintf("%.6f",value)
        formatted==row[key] || error("Original support/cache frame mismatch: $key")
    end
end

"Observed unsmoothed mean for fitting; the native smoothed mean locates the ROI."
function fused_data(img, pcfg, cfg, v, dx; original=nothing)
    f,sf=observed_shift(v.zf,v.cf,pcfg,0)
    b,sb=observed_shift(v.zb,v.cb,pcfg,dx)
    z=(f.+b)./2; zs=(sf.+sb)./2
    observed=isfinite.(z); count(observed)>0 || error("Empty common observed view support")
    if original===nothing
        _,_,mask=G.molecule_roi_mask_fused(img,pcfg,zs;observed_only=true)
        mask .&= observed
        count(mask)==0 && (mask=observed .& isfinite.(zs))
    else
        original.xs==v.xs && original.ys==v.ys || error("Original support grid differs")
        mask=original.mask .& observed
    end
    count(mask)>0 || error("No observed ROI pixels")
    offset=quantile(z[mask],0.05)
    zimg=z.-offset
    x,y,zfull=G._flatten_roi(zimg,mask,v.xs,v.ys)
    roi_indices=[CartesianIndex(iy,ix) for iy in eachindex(v.ys) for ix in eachindex(v.xs) if mask[iy,ix]]
    if original===nothing
        axisfull=G._weighted_roi_axis(x,y,zfull)
        xf,yf,zf,axis,keep,support=G._chain_fit_data(x,y,zfull,axisfull,cfg)
    else
        keep=original.fitmask[roi_indices]
        count(keep)>0 || error("No observed pixels in original fit support")
        xf,yf,zf=x[keep],y[keep],zfull[keep]
        axis,support=original.axis,original.support
    end
    all(isfinite,zf) || error("Unobserved fit pixel")
    # _flatten_roi uses y outer / x inner order, unlike Julia's vec.
    indices=roi_indices[keep]
    data=(xs=v.xs,ys=v.ys,zimg=zimg,x=xf,y=yf,z=zf,zfull=zfull,noise=v.noise,axisctx=axis)
    return (;data,mask,observed,indices,f,b,offset,support)
end

function fit_profile(n, bundle, raw, assignment, profile, options)
    model=deepcopy(raw["model"])
    model["peak_profile"]=profile
    profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
    _,ell,circ=F.Extractor._configs(model,raw["preprocessing"],"unused")
    isempty(ell.init_centers_t) && isempty(circ.init_centers_t) || error("No 1D initializer allowed")
    start=time_ns()
    result,source,cfg,boot=R.initialize(n,bundle.data,circ,ell,options;observed_only=true)
    elapsed_s=(time_ns()-start)/1e9
    result.valid && result.success && result.n==n || error("No valid fixed-N result")
    return (;result,source,cfg,boot,elapsed_s)
end

function execute(opts; reader=G.read_sxm)
    BLAS.set_num_threads(1)
    options=settings(opts["--settings"])
    raw=TOML.parsefile(opts["--config"]); assignment=TOML.parsefile(opts["--assignment-config"])
    model=raw["model"]
    get(model,"selection_criterion","")==get(model,"cv_method","")=="gcv" || error("GCV required")
    startswith(get(model,"selection_policy",""),"adaptive_support") && error("Adaptive support is outside this experiment")
    for k in ("global_maxtime","global_maxiter","max_iter","multistart")
        haskey(model,k) && model[k]>0 || error("Explicit positive native budget required: $k")
    end
    header,base=lobe_table(opts["--features"];required=F.REQUIRED)
    allfiles=sort(unique(first.(collect(keys(base)))))
    shifts=read_shifts(opts["--shifts"],allfiles)
    chunk=parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(allfiles) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty refit shard")
    chains=Dict(f=>[base[k] for k in sort(collect(keys(base))) if first(k)==f] for f in files)
    for f in files
        F.validate_chain(chains[f]); isfile(joinpath(opts["--data-dir"],f)) || error("Missing raw $f")
    end
    haskey(opts,"--dry-run") && return println("Registered refit: $(length(files)) files, native Gaussian/split at saved N; original_support=$(options.original_support); metadata only.")
    out=opts["--out"]; mkpath(out*".fit_data")
    outputs=Dict((a,p)=>Dict{String,String}[] for a in ARMS for p in PROFILES)
    audits=Dict{String,String}[]; parameters=Dict{String,String}[]; bootstrap=Dict{String,String}[]; failures=Dict{String,String}[]
    supports=Dict{String,String}[]
    for (i,file) in enumerate(files)
        stage="load"
        try
            n=length(chains[file]); pcfg,ell,_=F.Extractor._configs(model,raw["preprocessing"],dirname(out))
            pcfg.filepath=joinpath(opts["--data-dir"],file)
            img=reader(pcfg.filepath); views=load_views(img,pcfg)
            original=nothing
            if options.original_support
                stage="original_support"
                original=original_support(img,pcfg,ell)
                check_original_frame(original,chains[file][1])
                a=original.axis
                push!(supports,Dict("file"=>file,"axis_x"=>F.fmt(a.axis[1]),"axis_y"=>F.fmt(a.axis[2]),
                    "origin_x_nm"=>F.fmt(a.origin[1]),"origin_y_nm"=>F.fmt(a.origin[2]),
                    "support_tmin"=>F.fmt(a.tmin),"support_tmax"=>F.fmt(a.tmax),
                    "roi_pixels"=>string(count(original.mask)),"fit_pixels"=>string(count(original.fitmask)),
                    "support_method"=>original.support.support_method))
                open(joinpath(out*".fit_data",file*".original.tsv"),"w") do io
                    println(io,"x_nm\ty_nm\tz_native_nm\trow\tcolumn\tin_fit")
                    for j in eachindex(original.z)
                        I=original.indices[j]
                        println(io,join((F.fmt(original.x[j]),F.fmt(original.y[j]),F.fmt(original.z[j]),I[1],I[2],original.keep[j]),'\t'))
                    end
                end
            end
            control_bundle=nothing; control_fits=Dict()
            for arm in ARMS
                stage=arm*"_support"
                reused=arm=="registered" && shifts[file]==0
                dx=arm=="control" ? 0 : shifts[file]
                bundle=reused ? control_bundle : fused_data(img,pcfg,ell,views,dx;original)
                arm=="control" && (control_bundle=bundle)
                d=bundle.data
                open(joinpath(out*".fit_data",file*".$arm.tsv"),"w") do io
                    println(io,"x_nm\ty_nm\tz_nm\tfwd_nm\tbwd_nm\trow\tcolumn")
                    for j in eachindex(d.z)
                        I=bundle.indices[j]
                        println(io,join((F.fmt(d.x[j]),F.fmt(d.y[j]),F.fmt(d.z[j]),F.fmt(bundle.f[I]),F.fmt(bundle.b[I]),I[1],I[2]),'\t'))
                    end
                end
                for profile in PROFILES
                    stage=arm*"_"*profile
                    fit=reused ? control_fits[profile] : fit_profile(n,bundle,raw,assignment,profile,options)
                    arm=="control" && (control_fits[profile]=fit)
                    r,cfg=fit.result,fit.cfg
                    append!(outputs[(arm,profile)],R.feature_rows(file,r,d,cfg,fit.source))
                    lo,hi=VP.native_raw_bounds(n,cfg)
                    for j in eachindex(r.params)
                        push!(parameters,Dict("file"=>file,"arm"=>arm,"profile"=>profile,"parameter"=>string(j),
                            "value"=>F.fmt(r.params[j]),"lower"=>F.fmt(lo[j]),"upper"=>F.fmt(hi[j])))
                    end
                    for family in ("circ","ell")
                        boot=getproperty(fit.boot,Symbol(family))
                        push!(bootstrap,Dict("file"=>file,"arm"=>arm,"profile"=>profile,"source"=>family,
                            "success"=>string(boot.success),"valid"=>string(boot.valid),"gcv"=>F.fmt(boot.gcv),"reason"=>boot.reason))
                    end
                    a=Dict("file"=>file,"arm"=>arm,"profile"=>profile,"N"=>string(n),"source"=>fit.source,
                        "dx_px"=>string(dx),"reused_zero_shift"=>string(reused),"elapsed_s"=>reused ? "0" : F.fmt(fit.elapsed_s),
                        "n_data"=>string(length(d.z)),"p_full"=>string(G._chain_nparams(n,cfg)),
                        "observed_pixels"=>string(count(bundle.observed)),"roi_pixels"=>string(count(bundle.mask)),
                        "offset"=>F.fmt(bundle.offset),"noise"=>F.fmt(d.noise),"amp_min"=>F.fmt(r.amp_min),"amp_range"=>F.fmt(r.amp_range),
                        "amplitude_max_data"=>F.fmt(max(maximum(filter(isfinite,vec(d.zimg))),G.EPS)),
                        "axis_x"=>F.fmt(d.axisctx.axis[1]),"axis_y"=>F.fmt(d.axisctx.axis[2]),
                        "origin_x_nm"=>F.fmt(d.axisctx.origin[1]),"origin_y_nm"=>F.fmt(d.axisctx.origin[2]),
                        "support_tmin"=>F.fmt(d.axisctx.tmin),"support_tmax"=>F.fmt(d.axisctx.tmax),
                        "optimizer_convergence"=>"unknown_native_not_exposed")
                    for key in (:rss,:gcv,:valid,:reason,:overlap,:endpoint_overrun_nm,:residual_peak_snr,:kappa_max_adj)
                        a[string(key)]=string(getproperty(r,key))
                    end
                    push!(audits,a)
                    println("[$i/$(length(files))] $file N=$n $arm $profile $(fit.source) RSS=$(r.rss) reused=$reused"); flush(stdout)
                end
            end
        catch err
            reason=replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
            push!(failures,Dict("file"=>file,"stage"=>stage,"reason"=>reason))
            println(stderr,"FAILED $file $stage: $reason"); flush(stderr)
        end
    end
    for (suffix,rows) in ((".fits.tsv",audits),(".parameters.tsv",parameters),(".bootstrap.tsv",bootstrap),(".support.tsv",supports),(".failures.tsv",failures))
        isempty(rows) || write_table(out*suffix,sort(collect(keys(first(rows)))),rows)
    end
    isempty(failures) || error("$(length(failures)) files failed; no incomplete refit arm emitted")
    for arm in ARMS, profile in PROFILES
        write_table(out*".$arm.$profile.tsv",header,outputs[(arm,profile)])
    end
    write_table(out,header,outputs[("registered","gaussian")])
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("refit_registered_geometry.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--chunk I/N] [--dry-run]")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && RegisteredRefit.main()

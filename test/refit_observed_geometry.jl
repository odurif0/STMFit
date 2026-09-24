#!/usr/bin/env julia
# Fixed saved N; compare native imputation with observed-only background + fitting.
module ObservedGeometry
include(joinpath(@__DIR__, "refine_fixed_n_geometry.jl"))
module Pipeline
include(joinpath(@__DIR__, "run_reconstructed_chitosan.jl"))
end
using .FixedNGeometryRefinement: G, VP
using .Pipeline: selected_counts, raw_index, read_table, write_table
using TOML, Statistics, LinearAlgebra
const R = FixedNGeometryRefinement
const F = R.F
const ARMS = ("control", "observed")
const PROFILES = ("gaussian", "split")
const REQUIRED = Set(["--data-dir", "--selected-summary", "--config", "--assignment-config", "--settings", "--out"])

function settings(path)
    s = TOML.parsefile(path)
    expected = Dict(
        "model" => Dict("method"=>"native_circular_to_elliptical", "profiles"=>collect(PROFILES),
            "native_elliptical_maxiter"=>50, "one_dimensional_initialization"=>false),
        "selection" => Dict("count_policy"=>"reuse_saved_N", "family_policy"=>"valid_minimum_full_parameter_gcv",
            "identical_data_policy"=>"reuse_control_fit_exactly", "failure_policy"=>"report_and_stop_no_partial_grade"),
        "preprocessing" => Dict("missing_pixel_policy"=>"observed_only", "plane_rank_rtol"=>1e-12,
            "fusion"=>"mean_unsmoothed_observed_views", "roi"=>"native_finite_statistics_and_observed_mask", "patch_background"=>"same_as_fit"))
    s == expected || error("Only the declared observed-pixel comparison is supported")
    return (native_elliptical_maxiter=Int(s["model"]["native_elliptical_maxiter"]),
        plane_rank_rtol=Float64(s["preprocessing"]["plane_rank_rtol"]))
end

function parse_cli(args)
    args == ["--help"] && return nothing
    o=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]; haskey(o,k) && error("Repeated option $k")
        if k=="--dry-run"; o[k]="true"; i+=1; continue; end
        k in union(REQUIRED,Set(["--chunk"])) && i<length(args) && !startswith(args[i+1],"--") || error("Missing/forbidden option $k")
        o[k]=args[i+1]; i+=2
    end
    all(haskey(o,k) for k in REQUIRED) || error("All observed-fit inputs required")
    isdir(o["--data-dir"]) || error("Raw directory missing")
    all(isfile(o[k]) for k in setdiff(REQUIRED,Set(["--data-dir","--out"]))) || error("Missing input")
    for suffix in vcat(["", ".fits.tsv", ".parameters.tsv", ".diagnostics.tsv", ".families.tsv", ".failures.tsv", ".common.tsv", ".fit_data"],
            [".$a.$p.tsv" for a in ARMS for p in PROFILES])
        (ispath(o["--out"]*suffix) || islink(o["--out"]*suffix)) && error("Output exists")
    end
    chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    length(chunk)==2 && 1<=chunk[1]<=chunk[2] || error("Invalid chunk")
    settings(o["--settings"])
    o
end

function inputs(o)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    counts=selected_counts(o["--selected-summary"]); paths=raw_index(o["--data-dir"])
    Set(keys(counts))==Set(keys(paths)) || error("Selected-N/raw cohorts differ")
    raw=TOML.parsefile(o["--config"]); model=raw["model"]
    model["selection_criterion"]==model["cv_method"]=="gcv" || error("Native GCV required")
    startswith(model["selection_policy"],"adaptive_support") && error("Adaptive support outside this experiment")
    get(raw["preprocessing"],"missing_pixel_policy","median_fill")=="median_fill" || error("Control must use native preprocessing")
    for key in ("global_maxtime","global_maxiter","max_iter","multistart")
        haskey(model,key) && model[key]>0 || error("Explicit native budget required: $key")
    end
    assignment=Pipeline.load_config(o["--assignment-config"])
    assignment["preprocessing"]["patch_residual_filter"]=="smooth_residual" || error("Matched residual filter required")
    return (;counts,paths,raw,assignment)
end

"Both actual views; no get_channel fallback, fill, registration or second flip."
function observed_views(img, pcfg, options)
    views=map(("fwd","bwd")) do direction
        cfg=deepcopy(pcfg); cfg.direction=direction
        ch=G.get_channel(img,cfg.channel;direction)
        G.preprocess_observed_fit_channel(img,ch,cfg;plane_rank_rtol=options.plane_rank_rtol)
    end
    f,b=views
    f[1]==b[1] && f[2]==b[2] && f[6]==b[6] || error("Acquisition grids/units differ")
    return (;xs=f[1],ys=f[2],raw=(f[3].+b[3])./2,z=(f[4].+b[4])./2,
        zs=(f[5].+b[5])./2,noise=max(f[7],b[7],G.EPS),observed=isfinite.(f[3]).&isfinite.(b[3]))
end

function bundle(img, pcfg, cfg, arm, views)
    arm in ARMS || error("Unknown arm")
    if arm=="control"
        xs,ys,zimg,mask,x,y,z,noise=G._fused_roi_data(img,pcfg)
    else
        xs,ys,noise=views.xs,views.ys,views.noise
        _,_,mask=G.molecule_roi_mask_fused(img,pcfg,views.zs;observed_only=true)
        mask .&= views.observed
        any(mask) || (mask=views.observed .& isfinite.(views.zs))
        any(mask) || error("No observed ROI pixels")
        zimg=views.z.-quantile(views.z[mask],0.05)
        x,y,z=G._flatten_roi(zimg,mask,xs,ys)
    end
    full=G._weighted_roi_axis(x,y,z)
    xf,yf,zf,axis,keep,support=G._chain_fit_data(x,y,z,full,cfg)
    indices=[CartesianIndex(iy,ix) for iy in eachindex(ys) for ix in eachindex(xs) if mask[iy,ix]][keep]
    all(isfinite,zf) || error("Nonfinite fit pixel")
    fitmask=falses(size(mask)); fitmask[indices].=true
    arm=="observed" && any(fitmask .& .!views.observed) && error("Invented fit observation")
    data=(;xs,ys,zimg,x=xf,y=yf,z=zf,zfull=z,noise,axisctx=axis)
    return (;data,mask,indices,fitmask,support)
end

same_data(a,b) = isequal(a.data,b.data) && a.indices==b.indices && a.mask==b.mask && isequal(a.support,b.support)

function fit_profile(n, d, raw, assignment, profile, options; observed_only, diagnostics)
    model=deepcopy(raw["model"]); model["peak_profile"]=profile
    profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
    _,ell,circ=F.Extractor._configs(model,raw["preprocessing"],"unused")
    isempty(ell.init_centers_t) && isempty(circ.init_centers_t) || error("No 1D initializer")
    started=time_ns()
    result,source,cfg,boot=R.initialize(n,d,circ,ell,options;observed_only,diagnostics)
    return (;result,source,cfg,boot,elapsed_s=(time_ns()-started)/1e9)
end

function common_pixels(control, observed, fits, raw)
    indices=findall(control.fitmask .& observed.fitmask .& isfinite.(raw))
    isempty(indices) && error("No common observed fit pixels")
    x=[control.data.xs[I[2]] for I in indices]; y=[control.data.ys[I[1]] for I in indices]
    pred=map(ARMS) do arm
        b=arm=="control" ? control : observed; fit=fits[arm]; r=fit.result
        G._chain_model_values(x,y,r.params,r.n,b.data.axisctx,fit.cfg;amp_min=r.amp_min,amp_range=r.amp_range)
    end
    target=(control.data.zimg[indices],observed.data.zimg[indices])
    # Raw prediction = fitted molecular model + each arm's removed background.
    # Equivalently residual = processed target - fitted model on these SAME pixels.
    rss=map(i->sum(abs2,target[i].-pred[i]),1:2)
    return (;indices,x,y,pred,target,rss)
end

function execute(o; reader=G.read_sxm, fitter=fit_profile)
    BLAS.set_num_threads(1)
    input=inputs(o); options=settings(o["--settings"])
    allfiles=sort(collect(keys(input.counts))); chunk=parse.(Int,split(get(o,"--chunk","1/1"),'/'))
    files=[f for (i,f) in enumerate(allfiles) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty shard")
    haskey(o,"--dry-run") && return println("Observed-fit comparison: $(length(files)) files, saved N, two native profiles; metadata only")
    out=o["--out"]; mkpath(out*".fit_data")
    outputs=Dict((a,p)=>Dict{String,String}[] for a in ARMS for p in PROFILES)
    audits=Dict{String,String}[]; parameters=Dict{String,String}[]; diagnostics=Dict{String,String}[]
    failures=Dict{String,String}[]; common=Dict{String,String}[]; families=Dict{String,String}[]
    fail(file,stage,e)=push!(failures,Dict("file"=>file,"stage"=>stage,"reason"=>replace(sprint(showerror,e),'\n'=>' ','\t'=>' ')))
    for (i,file) in enumerate(files)
        n=input.counts[file]; bundles=Dict(); fits=Dict(); views=nothing
        try
            pcfg,ell,_=F.Extractor._configs(input.raw["model"],input.raw["preprocessing"],dirname(out))
            pcfg.filepath=input.paths[file]; img=reader(pcfg.filepath)
            views=observed_views(img,pcfg,options)
            for arm in ARMS
                try
                    b=bundle(img,pcfg,ell,arm,views); bundles[arm]=b
                    if arm=="observed" && all(views.observed)
                        same_data(bundles["control"],b) || error("Fully acquired inputs differ between arms")
                    end
                    open(joinpath(out*".fit_data",file*".$arm.tsv"),"w") do io
                        println(io,"x_nm\ty_nm\tz_nm\trow\tcolumn")
                        for j in eachindex(b.indices)
                            I=b.indices[j]; println(io,join((F.fmt(b.data.x[j]),F.fmt(b.data.y[j]),F.fmt(b.data.z[j]),I[1],I[2]),'\t'))
                        end
                    end
                    for profile in PROFILES
                        try
                            reused=arm=="observed" && all(views.observed)
                            callback=(family,row)->push!(diagnostics,merge(Dict("file"=>file,"arm"=>arm,"profile"=>profile,"family"=>family),
                                Dict(string(k)=>replace(string(v),'\n'=>' ','\t'=>' ') for (k,v) in pairs(row))))
                            fit=reused ? fits[("control",profile)] : fitter(n,b.data,input.raw,input.assignment,profile,options;
                                observed_only=arm=="observed",diagnostics=callback)
                            r=fit.result; r.success && r.valid && r.n==n || error("Invalid fixed-N fit")
                            fits[(arm,profile)]=fit
                            append!(outputs[(arm,profile)],R.feature_rows(file,r,b.data,fit.cfg,fit.source))
                            for family in ("circ","ell")
                                boot=getproperty(fit.boot,Symbol(family))
                                push!(families,Dict("file"=>file,"arm"=>arm,"profile"=>profile,"family"=>family,
                                    "success"=>string(boot.success),"valid"=>string(boot.valid),"rss"=>F.fmt(boot.rss),"gcv"=>F.fmt(boot.gcv),"reason"=>boot.reason))
                            end
                            lo,hi=VP.native_raw_bounds(n,fit.cfg)
                            for j in eachindex(r.params)
                                push!(parameters,Dict("file"=>file,"arm"=>arm,"profile"=>profile,"parameter"=>string(j),
                                    "value"=>F.fmt(r.params[j]),"lower"=>F.fmt(lo[j]),"upper"=>F.fmt(hi[j])))
                            end
                            a=Dict("file"=>file,"arm"=>arm,"profile"=>profile,"N"=>string(n),"source"=>fit.source,
                                "reused_control"=>string(reused),"elapsed_s"=>reused ? "0" : F.fmt(fit.elapsed_s),
                                "n_data"=>string(length(b.data.z)),"observed_pixels"=>string(count(views.observed)),"total_pixels"=>string(length(views.observed)),
                                "noise"=>F.fmt(b.data.noise),"amp_min"=>F.fmt(r.amp_min),"amp_range"=>F.fmt(r.amp_range),
                                "axis_x"=>F.fmt(b.data.axisctx.axis[1]),"axis_y"=>F.fmt(b.data.axisctx.axis[2]),
                                "origin_x_nm"=>F.fmt(b.data.axisctx.origin[1]),"origin_y_nm"=>F.fmt(b.data.axisctx.origin[2]),
                                "p_full"=>string(G._chain_nparams(n,fit.cfg)),"support_tmin"=>F.fmt(b.data.axisctx.tmin),"support_tmax"=>F.fmt(b.data.axisctx.tmax))
                            for key in (:rss,:gcv,:valid,:reason,:overlap,:endpoint_overrun_nm,:residual_peak_snr,:kappa_max_adj)
                                a[string(key)]=string(getproperty(r,key))
                            end
                            push!(audits,a)
                            println("[$i/$(length(files))] $file N=$n $arm $profile $(fit.source) valid=$(r.valid) reused=$reused"); flush(stdout)
                        catch err
                            fail(file,arm*"_"*profile,err)
                        end
                    end
                catch err
                    fail(file,arm*"_support",err)
                end
            end
            for profile in PROFILES
                all(haskey(fits,(arm,profile)) for arm in ARMS) || continue
                c=common_pixels(bundles["control"],bundles["observed"],Dict(a=>fits[(a,profile)] for a in ARMS),views.raw)
                push!(common,Dict("file"=>file,"profile"=>profile,"n_common"=>string(length(c.indices)),
                    "control_rss"=>F.fmt(c.rss[1]),"observed_rss"=>F.fmt(c.rss[2])))
                open(joinpath(out*".fit_data",file*".$profile.common.tsv"),"w") do io
                    println(io,"x_nm\ty_nm\traw_mean_nm\tcontrol_target_nm\tobserved_target_nm\tcontrol_model_nm\tobserved_model_nm")
                    for j in eachindex(c.indices)
                        println(io,join(F.fmt.((c.x[j],c.y[j],views.raw[c.indices[j]],c.target[1][j],c.target[2][j],c.pred[1][j],c.pred[2][j])),'\t'))
                    end
                end
            end
        catch err
            fail(file,"load_or_comparison",err)
        end
    end
    for (suffix,rows) in ((".fits.tsv",audits),(".parameters.tsv",parameters),(".diagnostics.tsv",diagnostics),(".families.tsv",families),(".failures.tsv",failures),(".common.tsv",common))
        isempty(rows) || write_table(out*suffix,sort(collect(keys(first(rows)))),rows)
    end
    isempty(failures) || error("$(length(failures)) failures; no incomplete features or partial grade emitted")
    header=sort(collect(keys(first(outputs[("control","gaussian")]))))
    for arm in ARMS, profile in PROFILES
        write_table(out*".$arm.$profile.tsv",header,outputs[(arm,profile)])
    end
    write_table(out,header,outputs[("observed","gaussian")])
end

function main(args=ARGS)
    o=parse_cli(args)
    o===nothing ? println("refit_observed_geometry.jl ",join(sort(collect(REQUIRED))," VALUE ")," VALUE [--chunk I/N] [--dry-run]") : execute(o)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ObservedGeometry.main()

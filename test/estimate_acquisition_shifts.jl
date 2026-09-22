#!/usr/bin/env julia
module EstimateAcquisitionShifts
using GaussianFit2D: read_sxm, get_channel, preprocess_channel, PatternConfig
using Statistics
include(joinpath(@__DIR__,"lib","acquisition_registration.jl"))
include(joinpath(@__DIR__,"lib","patch_preprocessing.jl"))
include(joinpath(@__DIR__,"lib","patch_acquisition.jl"))
include(joinpath(@__DIR__,"lib","reconstructed_unit_assignment.jl"))
using .AcquisitionRegistration, .PatchPreprocessing, .PatchAcquisition, .ReconstructedUnitAssignment

const OPTIONS = Set(["--data-dir","--features","--config","--settings","--out"])
function parse_cli(args)
    "--help" in args && return nothing
    out=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]
        k in OPTIONS && i<length(args) || error("Missing/forbidden acquisition option")
        haskey(out,k) && error("Repeated acquisition option")
        out[k]=args[i+1]; i+=2
    end
    Set(keys(out)) == OPTIONS || error("All acquisition inputs are required")
    isdir(out["--data-dir"]) || error("Missing raw directory")
    for k in ("--features","--config","--settings")
        isfile(out[k]) || error("Missing input $k")
    end
    (ispath(out["--out"]) || islink(out["--out"])) && error("Refusing acquisition output overwrite")
    return out
end

function estimate_file(path,pre,s)
    img=read_sxm(path)
    cf=get_channel(img,"Z";direction="fwd"); cb=get_channel(img,"Z";direction="bwd")
    require_direction(cf,"fwd"); require_direction(cb,"bwd")
    cfg=PatternConfig(filepath=path,channel="Z",direction="fwd",stride=pre.stride,
        flatten=pre.flatten,smooth_radius_px=pre.smooth_radius_px,no_plot=true)
    xs,ys,_,zf,_,_,_=preprocess_channel(img,cf,cfg)
    xb,yb,_,zb,_,_,_=preprocess_channel(img,cb,cfg)
    xs==xb && ys==yb || error("Registration grids differ")
    # Use UNSMOOTHED processed Z. The raw masks, not preprocess_channel's already
    # imputed `raw` output, decide which samples are observations.
    zf[.!isfinite.(cf.data[1:pre.stride:end,1:pre.stride:end])] .= NaN
    zb[.!isfinite.(cb.data[1:pre.stride:end,1:pre.stride:end])] .= NaN
    pixel=median(diff(xs))
    r=scan(zf,zb,pixel,s)
    return r,pixel,count(!isfinite,cf.data),count(!isfinite,cb.data)
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("estimate_acquisition_shifts.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE")
    s=load_settings(opts["--settings"]); pre=load_patch_preprocessing(opts["--config"])
    _,features=lobe_table(opts["--features"])
    files=sort(unique(first.(collect(keys(features)))))
    summaries=Dict{String,String}[]; scores=Dict{String,String}[]; peaks=Dict{String,String}[]
    transforms=Dict{String,String}[]
    for (i,file) in enumerate(files)
        r,pixel,missing_f,missing_b=estimate_file(joinpath(opts["--data-dir"],file),pre,s)
        push!(transforms,Dict("file"=>file,"bwd_sample_dx_px"=>string(r.applied_dx_px)))
        push!(summaries,Dict("file"=>file,"accepted"=>string(r.accepted),"status"=>r.status,
            "best_dx_px"=>string(first(r.peaks).best_dx_px),"applied_dx_px"=>string(r.applied_dx_px),
            "pixel_nm"=>string(pixel),"applied_dx_nm"=>string(r.applied_dx_px*pixel),
            "limit_px"=>string(r.limit_px),"tolerance_px"=>string(r.tolerance_px),
            "agreement_px"=>string(r.agreement_px),"support_pixels"=>string(count(r.support)),
            "raw_missing_fwd"=>string(missing_f),"raw_missing_bwd"=>string(missing_b)))
        for (j,p) in enumerate(r.peaks)
            row=Dict(string(k)=>string(v) for (k,v) in pairs(p)); row["file"]=file; row["band"]=string(j-1)
            push!(peaks,row)
        end
        for band in 0:s["bands"], (j,lag) in enumerate(r.lags)
            push!(scores,Dict("file"=>file,"band"=>string(band),"dx_px"=>string(lag),
                "correlation"=>string(r.correlations[band+1,j])))
        end
        println("[$i/$(length(files))] $file dx=$(r.applied_dx_px), $(r.status)"); flush(stdout)
    end
    out=opts["--out"]
    write_table(out,["file","bwd_sample_dx_px"],transforms)
    write_table(out*".summary.tsv",["file","accepted","status","best_dx_px","applied_dx_px","pixel_nm",
        "applied_dx_nm","limit_px","tolerance_px","agreement_px","support_pixels","raw_missing_fwd","raw_missing_bwd"],summaries)
    write_table(out*".peaks.tsv",["file","band","best_dx_px","correlation","distant_peak_gap","boundary",
        "pixels","rows","usable","status"],peaks)
    write_table(out*".scores.tsv",["file","band","dx_px","correlation"],scores)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && EstimateAcquisitionShifts.main()

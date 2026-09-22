#!/usr/bin/env julia
# A focused feasibility check, not a benchmark and never a partial classifier.
module PairedAcquisitionDiagnostic
include(joinpath(@__DIR__,"diagnose_registered_fit_failure.jl"))
include(joinpath(@__DIR__,"lib","paired_acquisition_fit.jl"))
using .RegisteredFailureDiagnostic: RR, G, F, R, read_table, write_table, lobe_table
using .PairedAcquisitionFit
using TOML, SHA, LinearAlgebra
const D=RegisteredFailureDiagnostic
const P=PairedAcquisitionFit

function parse_cli(args)
    "--help" in args && return nothing
    extra=Dict{String,String}(); rest=String[]; i=1
    while i<=length(args)
        if args[i] in ("--native-diagnostic","--paired-settings")
            k=args[i]; i<length(args) && !haskey(extra,k) || error("Missing/repeated paired option")
            extra[k]=args[i+1]; i+=2
        else
            push!(rest,args[i]); i+=1
        end
    end
    Set(keys(extra))==Set(["--native-diagnostic","--paired-settings"]) || error("Paired settings and native diagnostic required")
    opts=D.parse_cli(rest); merge!(opts,extra)
    isfile(joinpath(opts["--native-diagnostic"],"fits.tsv")) || error("Completed native diagnostic required")
    P.settings(opts["--paired-settings"])
    return opts
end

function execute(opts;reader=G.read_sxm)
    BLAS.set_num_threads(1)
    settings=P.settings(opts["--paired-settings"])
    native=opts["--native-diagnostic"]; file=opts["--file"]
    _,saved=read_table(joinpath(native,"fits.tsv"))
    _,params=read_table(joinpath(native,"parameters.tsv"))
    _,hashes=read_table(joinpath(native,"input_hashes.tsv"))
    for row in hashes
        key=row["input"]
        haskey(opts,key) && bytes2hex(sha256(read(opts[key])))==row["sha256"] || error("Native diagnostic input changed: $key")
    end
    all(r["file"]==file for r in saved) || error("Native diagnostic contains another scan")
    _,base=lobe_table(opts["--features"];required=F.REQUIRED)
    chain=[base[k] for k in sort(collect(keys(base))) if first(k)==file]; F.validate_chain(chain); n=length(chain)
    shifts=RR.read_shifts(opts["--shifts"],sort(unique(first.(collect(keys(base))))))
    haskey(opts,"--dry-run") && return println("Paired feasibility: one scan $file, saved N=$n; same registered starts, matched fused/paired LM, no output.")
    out=opts["--outdir"]; mkpath(out)
    cp(opts["--paired-settings"],joinpath(out,"paired_settings.toml"))
    write_table(joinpath(out,"native_hashes.tsv"),["input","sha256"],
        [Dict("input"=>s,"sha256"=>bytes2hex(sha256(read(joinpath(native,s))))) for s in ("fits.tsv","parameters.tsv","registered.pixels.tsv","input_hashes.tsv")])
    raw=TOML.parsefile(opts["--config"]); assignment=TOML.parsefile(opts["--assignment-config"])
    pcfg,ell,_=F.Extractor._configs(raw["model"],raw["preprocessing"],out)
    img=reader(joinpath(opts["--data-dir"],file)); views=RR.load_views(img,pcfg)
    original=RR.original_support(img,pcfg,ell); RR.check_original_frame(original,first(chain))
    bundle=RR.fused_data(img,pcfg,ell,views,shifts[file];original); d=bundle.data
    _,pixels=read_table(joinpath(native,"registered.pixels.tsv"))
    length(pixels)==length(d.z) || error("Native objective support changed")
    for (j,I) in enumerate(bundle.indices)
        expected=Dict("row"=>string(I[1]),"column"=>string(I[2]),"x_nm"=>F.fmt(d.x[j]),"y_nm"=>F.fmt(d.y[j]),
            "z_nm"=>F.fmt(d.z[j]),"fwd_nm"=>F.fmt(bundle.f[I]),"bwd_nm"=>F.fmt(bundle.b[I]))
        pixels[j]==expected || error("Native objective changed at $j")
    end
    fwd=bundle.f[bundle.indices].-bundle.offset; bwd=bundle.b[bundle.indices].-bundle.offset
    rows=Dict{String,String}[]; allparameters=Dict{String,String}[]; failures=Dict{String,String}[]
    for profile in RR.PROFILES, family in ("circ","ell")
        model=deepcopy(raw["model"]); model["peak_profile"]=profile
        profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
        _,ec,cc=F.Extractor._configs(model,raw["preprocessing"],out); cfg=family=="circ" ? cc : ec
        row=only(r for r in saved if r["arm"]=="registered" && r["profile"]==profile && r["family"]==family)
        row["N"]==string(n) && row["success"]=="true" || error("No numerical native candidate")
        nativeparams=sort(filter(r->r["arm"]=="registered" && r["profile"]==profile && r["family"]==family && r["stage"]=="final",params);
            by=r->parse(Int,r["parameter"]))
        length(nativeparams)==G._chain_nparams(n,cfg) && all(r["start"]=="1" for r in nativeparams) || error("Unexpected native starts")
        initial=G.ChainModelResult(n=n,params=parse.(Float64,getindex.(nativeparams,"value")),success=true,
            amp_min=parse(Float64,row["amp_min"]),amp_range=parse(Float64,row["amp_range"]))
        for mode in ("fused","paired")
            prefix="$profile.$family.$mode"
            try
                fit=P.fit(initial,d,cfg,fwd,bwd,mode,settings)
                r=fit.result
                summary=Dict("file"=>file,"profile"=>profile,"family"=>family,"mode"=>mode,"N"=>string(n),
                    "mean_rss"=>F.fmt(r.rss),"mean_peak_snr"=>F.fmt(r.residual_peak_snr),
                    "amp_min"=>F.fmt(r.amp_min),"amp_range"=>F.fmt(r.amp_range),"noise"=>F.fmt(d.noise),
                    "threshold"=>F.fmt(cfg.residual_peak_snr_threshold),"native_valid"=>string(r.valid))
                for key in (:valid,:reason,:rss,:gcv,:nd,:np,:rss_fwd,:rss_bwd,:peak_fwd,:peak_bwd,:initial_rss,:elapsed_s,:converged,:iterations)
                    summary[string(key)]=string(getproperty(fit,key))
                end
                ps=[Dict("profile"=>profile,"family"=>family,"mode"=>mode,"parameter"=>string(j),
                    "value"=>F.fmt(fit.params[j]),"lower"=>F.fmt(fit.lower[j]),"upper"=>F.fmt(fit.upper[j])) for j in eachindex(fit.params)]
                push!(rows,summary); append!(allparameters,ps)
                write_table(joinpath(out,prefix*".fit.tsv"),sort(collect(keys(summary))),[summary])
                write_table(joinpath(out,prefix*".parameters.tsv"),sort(collect(keys(first(ps)))),ps)
                residuals=[Dict("row"=>string(I[1]),"column"=>string(I[2]),"mean_prediction"=>F.fmt(fit.predictions.mean[j]),
                    "fwd_prediction"=>F.fmt(fit.predictions.fwd[j]),"bwd_prediction"=>F.fmt(fit.predictions.bwd[j]),
                    "mean_residual"=>F.fmt(d.z[j]-fit.predictions.mean[j]),"fwd_residual"=>F.fmt(fwd[j]-fit.predictions.fwd[j]),
                    "bwd_residual"=>F.fmt(bwd[j]-fit.predictions.bwd[j])) for (j,I) in enumerate(bundle.indices)]
                write_table(joinpath(out,prefix*".residuals.tsv"),sort(collect(keys(first(residuals)))),residuals)
                println("$prefix valid=$(fit.valid) mean_peak=$(r.residual_peak_snr) views=$(fit.peak_fwd)/$(fit.peak_bwd) converged=$(fit.converged) iterations=$(fit.iterations)"); flush(stdout)
            catch err
                push!(failures,Dict("stage"=>prefix,"reason"=>replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')))
                println(stderr,"FAILED $prefix: ",sprint(showerror,err)); flush(stderr)
            end
        end
    end
    for (name,records) in (("fits",rows),("parameters",allparameters),("failures",failures))
        isempty(records) || write_table(joinpath(out,"$name.tsv"),sort(collect(keys(first(records)))),records)
    end
    eligibility=Dict{String,String}[]
    for mode in ("fused","paired"),profile in RR.PROFILES
        valid=filter(r->r["mode"]==mode && r["profile"]==profile && r["valid"]=="true",rows)
        sort!(valid;by=r->parse(Float64,r["gcv"]))
        push!(eligibility,Dict("mode"=>mode,"profile"=>profile,"valid_families"=>string(length(valid)),
            "selected_family"=>isempty(valid) ? "none" : first(valid)["family"]))
    end
    write_table(joinpath(out,"eligibility.tsv"),["mode","profile","valid_families","selected_family"],eligibility)
    println("Feasibility complete. No classifier or external grade. All four required selections available: ",
        isempty(failures) && all(r["selected_family"]!="none" for r in eligibility))
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_paired_acquisition.jl: native diagnostic options plus --native-diagnostic DIR --paired-settings TOML")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && PairedAcquisitionDiagnostic.main()

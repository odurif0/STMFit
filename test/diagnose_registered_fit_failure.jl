#!/usr/bin/env julia
# One input file selected by an existing numerical failure, never by a grade.
module RegisteredFailureDiagnostic
include(joinpath(@__DIR__, "refit_registered_geometry.jl"))
using .RegisteredRefit: G, VP, F, R, read_table, lobe_table, write_table
using TOML, LinearAlgebra, Statistics, SHA
const RR = RegisteredRefit
const OPTIONS = Set(["--file", "--data-dir", "--features", "--shifts", "--config",
    "--assignment-config", "--settings", "--outdir"])

function parse_cli(args)
    "--help" in args && return nothing
    opts=Dict{String,String}(); i=1
    while i<=length(args)
        k=args[i]
        if k=="--dry-run"; haskey(opts,k) && error("Repeated option"); opts[k]="true"; i+=1; continue; end
        k in OPTIONS && !haskey(opts,k) && i<length(args) || error("Missing/forbidden diagnostic option")
        opts[k]=args[i+1]; i+=2
    end
    all(haskey(opts,k) for k in OPTIONS) || error("All diagnostic inputs required")
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    basename(opts["--file"])==opts["--file"] || error("A single basename is required")
    for k in ("--features","--shifts","--config","--assignment-config","--settings")
        isfile(opts[k]) || error("Missing $k")
    end
    isfile(joinpath(opts["--data-dir"],opts["--file"])) || error("Missing raw scan")
    (ispath(opts["--outdir"]) || islink(opts["--outdir"])) && error("Output exists")
    RR.settings(opts["--settings"]).original_support || error("Original support required")
    return opts
end

function execute(opts;reader=G.read_sxm)
    BLAS.set_num_threads(1)
    file=opts["--file"]; _,base=lobe_table(opts["--features"];required=F.REQUIRED)
    chain=[base[k] for k in sort(collect(keys(base))) if first(k)==file]
    isempty(chain) && error("Missing cached scan"); F.validate_chain(chain); n=length(chain)
    allfiles=sort(unique(first.(collect(keys(base)))))
    shifts=RR.read_shifts(opts["--shifts"],allfiles)
    raw=TOML.parsefile(opts["--config"]); assignment=TOML.parsefile(opts["--assignment-config"])
    options=RR.settings(opts["--settings"])
    get(raw["model"],"cv_method","")==get(raw["model"],"selection_criterion","")=="gcv" || error("GCV required")
    haskey(opts,"--dry-run") && return println("One-file diagnostic: $file, saved N=$n; no fit/output/labels.")
    out=opts["--outdir"]; mkpath(out)
    write_table(joinpath(out,"input_hashes.tsv"),["input","sha256"],
        [Dict("input"=>k,"sha256"=>bytes2hex(sha256(read(opts[k])))) for k in sort(collect(OPTIONS)) if isfile(opts[k])])
    pcfg,ell,_=F.Extractor._configs(raw["model"],raw["preprocessing"],out)
    img=reader(joinpath(opts["--data-dir"],file))
    views=RR.load_views(img,pcfg); original=RR.original_support(img,pcfg,ell)
    RR.check_original_frame(original,first(chain))
    summaries=Dict{String,String}[]; parameters=Dict{String,String}[]; optimizers=Dict{String,String}[]
    for arm in RR.ARMS
        dx=arm=="control" ? 0 : shifts[file]
        bundle=RR.fused_data(img,pcfg,ell,views,dx;original); d=bundle.data
        write_table(joinpath(out,"$arm.pixels.tsv"),["row","column","x_nm","y_nm","z_nm","fwd_nm","bwd_nm"],
            [Dict("row"=>string(I[1]),"column"=>string(I[2]),"x_nm"=>F.fmt(d.x[j]),"y_nm"=>F.fmt(d.y[j]),
                "z_nm"=>F.fmt(d.z[j]),"fwd_nm"=>F.fmt(bundle.f[I]),"bwd_nm"=>F.fmt(bundle.b[I])) for (j,I) in enumerate(bundle.indices)])
        for profile in RR.PROFILES
            model=deepcopy(raw["model"]); model["peak_profile"]=profile
            profile=="split" && (model["skew_ratio_max"]=assignment["model"]["split_skew_ratio_max"])
            _,ec,cc=F.Extractor._configs(model,raw["preprocessing"],out)
            rc=nothing
            for family in ("circ","ell")
                cfg=family=="circ" ? cc : deepcopy(ec)
                observer=function(row)
                    r=Dict("arm"=>arm,"profile"=>profile,"family"=>family)
                    for k in (:start,:rss,:global_status,:global_error,:global_evaluations,:lm_status,:lm_error,:lm_converged,:lm_iterations)
                        r[string(k)]=replace(string(getproperty(row,k)),'\n'=>' ','\t'=>' ')
                    end
                    push!(optimizers,r)
                    for (stage,p) in (("initial",row.initial),("global",row.global_params),("final",row.params))
                        for j in eachindex(p)
                            push!(parameters,Dict("arm"=>arm,"profile"=>profile,"family"=>family,"stage"=>stage,
                                "start"=>string(row.start),"parameter"=>string(j),"value"=>F.fmt(p[j]),
                                "amp_min"=>F.fmt(row.amp_min),"amp_range"=>F.fmt(row.amp_range)))
                        end
                    end
                end
                started=time_ns()
                if family=="circ"
                    r=G._fit_chain_n(d.xs,d.ys,d.zimg,d.x,d.y,d.z,d.noise,n,d.axisctx,cfg;observed_only=true,diagnostics=observer)
                    rc=r
                elseif rc.success
                    cfg.skip_global=true; cfg.multistart=1; cfg.max_iter=options.native_elliptical_maxiter
                    pe=VP.circular_to_elliptical(rc.params,n,cfg)
                    r=G._fit_chain_n(d.xs,d.ys,d.zimg,d.x,d.y,d.z,d.noise,n,d.axisctx,cfg;
                        starts=1,warm_start=pe,observed_only=true,diagnostics=observer)
                else
                    r=G.ChainModelResult(n=n,reason="circular initialization failed")
                end
                VP.finalize_native!(r,d,cfg)
                a=d.axisctx
                summary=Dict("file"=>file,"arm"=>arm,"profile"=>profile,"family"=>family,"N"=>string(n),
                    "dx_px"=>string(dx),"n_pixels"=>string(length(d.z)),"noise"=>F.fmt(d.noise),"offset"=>F.fmt(bundle.offset),
                    "threshold"=>F.fmt(cfg.residual_peak_snr_threshold),"elapsed_s"=>F.fmt((time_ns()-started)/1e9),
                    "axis_x"=>F.fmt(a.axis[1]),"axis_y"=>F.fmt(a.axis[2]),"origin_x_nm"=>F.fmt(a.origin[1]),
                    "origin_y_nm"=>F.fmt(a.origin[2]),"support_tmin"=>F.fmt(a.tmin),"support_tmax"=>F.fmt(a.tmax))
                for k in (:success,:valid,:reason,:rss,:gcv,:residual_peak_snr,:overlap,:endpoint_overrun_nm,:amp_min,:amp_range)
                    summary[string(k)]=string(getproperty(r,k))
                end
                if r.success
                    pred=G._chain_model_values(d.x,d.y,r.params,n,a,cfg;amp_min=r.amp_min,amp_range=r.amp_range)
                    residual=d.z.-pred; jmax=argmax(abs.(residual))
                    summary["exceeding_pixels"]=string(count(>(cfg.residual_peak_snr_threshold*d.noise),abs.(residual)))
                    summary["peak_row"]=string(bundle.indices[jmax][1]); summary["peak_column"]=string(bundle.indices[jmax][2])
                    summary["p99_abs_residual_snr"]=F.fmt(quantile(abs.(residual)./d.noise,.99))
                    write_table(joinpath(out,"$arm.$profile.$family.residuals.tsv"),["row","column","prediction_nm","residual_nm","residual_snr"],
                        [Dict("row"=>string(I[1]),"column"=>string(I[2]),"prediction_nm"=>F.fmt(pred[j]),
                            "residual_nm"=>F.fmt(residual[j]),"residual_snr"=>F.fmt(residual[j]/d.noise)) for (j,I) in enumerate(bundle.indices)])
                else
                    for key in ("exceeding_pixels","peak_row","peak_column","p99_abs_residual_snr"); summary[key]="NA"; end
                end
                push!(summaries,summary)
                println("$arm $profile $family valid=$(r.valid) max_residual/noise=$(r.residual_peak_snr) reason=$(r.reason)"); flush(stdout)
                # Save every completed stage; a diagnostic failure cannot erase preceding evidence.
                for (name,rows) in (("fits",summaries),("parameters",parameters),("optimizers",optimizers))
                    part=filter(r->r["arm"]==arm && r["profile"]==profile && r["family"]==family,rows)
                    isempty(part) || write_table(joinpath(out,"$arm.$profile.$family.$name.tsv"),sort(collect(keys(first(part)))),part)
                end
            end
        end
    end
    for (name,rows) in (("fits",summaries),("parameters",parameters),("optimizers",optimizers))
        isempty(rows) || write_table(joinpath(out,"$name.tsv"),sort(collect(keys(first(rows)))),rows)
    end
    println("Diagnostic complete. No selection override, classifier, benchmark labels or grade.")
end

function main(args=ARGS)
    opts=parse_cli(args)
    opts===nothing && return println("diagnose_registered_fit_failure.jl ",join(sort(collect(OPTIONS))," VALUE ")," VALUE [--dry-run]")
    execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && RegisteredFailureDiagnostic.main()

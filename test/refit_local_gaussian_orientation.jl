#!/usr/bin/env julia
# One fresh global initialization; matched global/local elliptical continuations.
module LocalGaussianOrientation
include(joinpath(@__DIR__, "refine_fixed_n_geometry.jl"))
using .FixedNGeometryRefinement: G, VP, F, read_table, lobe_table, write_table
using LinearAlgebra, TOML, Printf
const R = FixedNGeometryRefinement
const AXIS_FIELDS = ["model_orientation", "model_axis_x", "model_axis_y"]

function settings(path)
    raw = TOML.parsefile(path)
    expected = Dict("model"=>Set(["chain_peak_orientation","chain_tangent_degree"]),
        "selection"=>Set(["criterion","initialization_elliptical_maxiter","continuation_maxiter","tie_breaker"]),
        "preprocessing"=>Set(["patch_orientation","residual_geometry"]))
    Set(keys(raw)) == Set(keys(expected)) || error("Unexpected orientation settings")
    for (s,fields) in expected
        raw[s] isa AbstractDict && Set(keys(raw[s])) == fields || error("Unexpected $s settings")
    end
    raw["model"]["chain_peak_orientation"] == "local_tangent" || error("Local-tangent candidate required")
    degree = raw["model"]["chain_tangent_degree"]
    degree isa Integer && !(degree isa Bool) && degree in (1,2) || error("Degree must be 1 or 2")
    for k in ("initialization_elliptical_maxiter","continuation_maxiter")
        v = raw["selection"][k]
        v isa Integer && !(v isa Bool) && v > 0 || error("Positive iteration limit required")
    end
    raw["selection"]["criterion"] == "valid_min_gcv" || error("Valid minimum GCV required")
    raw["selection"]["tie_breaker"] == "global" || error("Global tie precedence required")
    raw["preprocessing"]["patch_orientation"] == "global" || error("Patch rotation is outside this comparison")
    raw["preprocessing"]["residual_geometry"] == "exported_model_axes" || error("Fitted axes must reach subtraction")
    return (degree=Int(degree), initial_maxiter=Int(raw["selection"]["initialization_elliptical_maxiter"]),
        continuation_maxiter=Int(raw["selection"]["continuation_maxiter"]))
end

function parse_cli(args)
    "--help" in args && return println("refit_local_gaussian_orientation.jl --features TSV --data-dir DIR --config TOML --settings TOML --out NEW_PREFIX [--chunk I/N] [--dry-run]")
    opts = F.parse_cli(args)
    settings(opts["--settings"])
    for suffix in (".global.tsv",".gcv.tsv",".fits.tsv",".parameters.tsv",".selections.tsv")
        (ispath(opts["--out"]*suffix) || islink(opts["--out"]*suffix)) && error("Output exists")
    end
    return opts
end

function fit_endpoint(n, data, cfg; warm_start=nothing)
    observations = Any[]
    started = time_ns()
    r = G._fit_chain_n(data.xs,data.ys,data.zimg,data.x,data.y,data.z,data.noise,n,data.axisctx,cfg;
        warm_start, diagnostics=row->push!(observations,row))
    VP.finalize_native!(r,data,cfg)
    return (result=r, diagnostic=isempty(observations) ? nothing : only(observations),
        elapsed_s=(time_ns()-started)/1e9)
end

function best_valid(candidates)
    pool = filter(c->c.fit.result.success && c.fit.result.valid && isfinite(c.fit.result.gcv),candidates)
    isempty(pool) && error("No valid model at saved N; no partial cohort or invalid fallback")
    # Stable first minimum: all global candidates precede the local candidate.
    return pool[argmin([c.fit.result.gcv for c in pool])]
end

function fit_candidates(n,data,ell,circ,opt; fitter=fit_endpoint)
    all(c->c.chain_peak_orientation=="global" && G._chain_peak_profile(c)==:gaussian,(ell,circ)) ||
        error("Shared initialization must be Gaussian/global")
    candidates = Any[]
    fc = fitter(n,data,circ)
    push!(candidates,(name="initial_circular",cfg=circ,fit=fc))
    if fc.result.success
        initial_cfg = deepcopy(ell); initial_cfg.skip_global = true
        initial_cfg.multistart = 1; initial_cfg.max_iter = opt.initial_maxiter
        pe = VP.circular_to_elliptical(fc.result.params,n,initial_cfg)
        fe = fitter(n,data,initial_cfg;warm_start=copy(pe))
        push!(candidates,(name="initial_elliptical",cfg=initial_cfg,fit=fe))
        if fe.result.success
            # Same full-precision vector, data, bounds, solver and iteration cap.
            for mode in ("global","local_tangent")
                cfg = deepcopy(ell); cfg.skip_global = true; cfg.multistart = 1
                cfg.max_iter = opt.continuation_maxiter
                cfg.chain_peak_orientation = mode; cfg.chain_tangent_degree = opt.degree
                fitted = fitter(n,data,cfg;warm_start=copy(fe.result.params))
                push!(candidates,(name=mode*"_elliptical",cfg,fit=fitted))
            end
        end
    end
    control = best_valid(filter(c->c.cfg.chain_peak_orientation=="global",candidates))
    selected = best_valid(candidates)
    selected.fit.result.gcv <= control.fit.result.gcv || error("Expanded GCV pool worsened its minimum")
    return (; candidates, control, selected)
end

function export_rows(file, candidate, data)
    r,cfg = candidate.fit.result,candidate.cfg
    source = cfg.chain_circular_sigmas ? "circ" : "ell"
    rows = R.feature_rows(file,r,data,cfg,source)
    _,_,ts,us = G._decode_chain(r.params,r.n,data.axisctx,cfg;amp_min=r.amp_min,amp_range=r.amp_range)
    axes = G._chain_peak_axes(ts,us,data.axisctx,cfg)
    for (row,(a,b)) in zip(rows,axes)
        row["model_orientation"] = cfg.chain_peak_orientation
        # Preserve the old export's rounded global-axis norm. Centers/widths
        # already use the native export precision; axes are not re-estimated
        # from those separately rounded columns.
        ax,ay = F.number(row,"axis_x"),F.number(row,"axis_y")
        if cfg.chain_peak_orientation == "global"
            row["model_axis_x"] = row["axis_x"]; row["model_axis_y"] = row["axis_y"]
        else
            scale = hypot(ax,ay)/hypot(a,b)
            row["model_axis_x"] = F.fmt(a*scale); row["model_axis_y"] = F.fmt(b*scale)
        end
    end
    return rows
end

function execute(opts)
    BLAS.set_num_threads(1)
    opt = settings(opts["--settings"]); raw = TOML.parsefile(opts["--config"])
    model,pre = raw["model"],raw["preprocessing"]
    all(get(model,k,"")=="gcv" for k in ("selection_criterion","cv_method")) || error("GCV required")
    get(model,"multistart",1)==1 || error("One shared initialization required")
    get(model,"peak_profile","gaussian")=="gaussian" || error("Gaussian base only")
    startswith(get(model,"selection_policy",""),"adaptive_support") && error("Adaptive support is outside this comparison")
    header,base = lobe_table(opts["--features"];required=F.REQUIRED)
    isempty(intersect(Set(header),Set(AXIS_FIELDS))) || error("Original global cache required")
    files = sort(unique(first.(collect(keys(base)))))
    chunk = parse.(Int,split(get(opts,"--chunk","1/1"),'/'))
    files = [f for (i,f) in enumerate(files) if mod1(i,chunk[2])==chunk[1]]
    isempty(files) && error("Empty chunk")
    chains = Dict(f=>[base[k] for k in sort(collect(keys(base))) if first(k)==f] for f in files)
    for f in files
        F.validate_chain(chains[f]); isfile(joinpath(opts["--data-dir"],f)) || error("Missing raw $f")
    end
    haskey(opts,"--dry-run") && return println("Orientation refit: $(length(files)) files; saved N; metadata only, no images or output.")
    out = opts["--out"]; mkpath(out*".fit_data")
    outputs = Dict(mode=>Dict{String,String}[] for mode in ("global","gcv"))
    fits = Dict{String,String}[]; parameters = Dict{String,String}[]
    selections = Dict{String,String}[]; failures = Dict{String,String}[]
    for (idx,file) in enumerate(files)
        stage = "input"
        try
            n = length(chains[file])
            pcfg,ell,circ = F.Extractor._configs(model,pre,dirname(out))
            pcfg.filepath = joinpath(opts["--data-dir"],file)
            img = G.read_sxm(pcfg.filepath)
            xs,ys,zimg,mask,x,y,z,noise = G._fused_roi_data(img,pcfg)
            full = G._weighted_roi_axis(x,y,z)
            xf,yf,zf,axis,keep,support = G._chain_fit_data(x,y,z,full,ell)
            data = (xs=xs,ys=ys,zimg=zimg,x=xf,y=yf,z=zf,zfull=z,noise=noise,axisctx=axis)
            for (key,v,digits) in (("axis_x",axis.axis[1],8),("axis_y",axis.axis[2],8),
                ("origin_x_nm",axis.origin[1],6),("origin_y_nm",axis.origin[2],6))
                (digits==8 ? @sprintf("%.8f",v) : @sprintf("%.6f",v)) == chains[file][1][key] || error("Input/cache axis mismatch")
            end
            open(joinpath(out*".fit_data",file*".tsv"),"w") do io
                println(io,"x_nm\ty_nm\tz_nm")
                for j in eachindex(zf); println(io,join(F.fmt.((xf[j],yf[j],zf[j])),'\t')); end
            end
            stage = "fit"
            compared = fit_candidates(n,data,ell,circ,opt)
            for c in compared.candidates
                r,d,cfg = c.fit.result,c.fit.diagnostic,c.cfg
                audit = Dict("file"=>file,"N"=>string(n),"model"=>c.name,"orientation"=>cfg.chain_peak_orientation,
                    "circular"=>string(cfg.chain_circular_sigmas),"success"=>string(r.success),"valid"=>string(r.valid),
                    "reason"=>r.reason,"rss"=>F.fmt(r.rss),"gcv"=>F.fmt(r.gcv),"p_full"=>string(G._chain_nparams(n,cfg)),
                    "n_data"=>string(length(zf)),"noise"=>F.fmt(noise),"amp_min"=>F.fmt(r.amp_min),"amp_range"=>F.fmt(r.amp_range),
                    "axis_x"=>F.fmt(axis.axis[1]),"axis_y"=>F.fmt(axis.axis[2]),"origin_x_nm"=>F.fmt(axis.origin[1]),
                    "origin_y_nm"=>F.fmt(axis.origin[2]),"support_tmin"=>F.fmt(axis.tmin),"support_tmax"=>F.fmt(axis.tmax),
                    "overlap"=>F.fmt(r.overlap),"endpoint_overrun_nm"=>F.fmt(r.endpoint_overrun_nm),
                    "residual_peak_snr"=>F.fmt(r.residual_peak_snr),"kappa"=>F.fmt(r.kappa_max_adj),
                    "elapsed_s"=>F.fmt(c.fit.elapsed_s),"max_iter"=>string(cfg.max_iter),"skip_global"=>string(cfg.skip_global))
                for key in (:global_status,:global_error,:global_evaluations,:lm_status,:lm_error,:lm_converged,:lm_iterations)
                    audit[string(key)] = d===nothing ? "NA" : replace(string(getproperty(d,key)),'\n'=>' ','\t'=>' ')
                end
                push!(fits,audit)
                lo,hi = VP.native_raw_bounds(n,cfg)
                for j in eachindex(r.params)
                    push!(parameters,Dict("file"=>file,"model"=>c.name,"parameter"=>string(j),"value"=>F.fmt(r.params[j]),
                        "initial"=>d===nothing ? "NA" : F.fmt(d.initial[j]),
                        "lm_start"=>d===nothing ? "NA" : F.fmt(d.global_params[j]),"lower"=>F.fmt(lo[j]),"upper"=>F.fmt(hi[j])))
                end
            end
            for c in compared.candidates
                d = c.fit.diagnostic
                d === nothing && continue # Explicit support infeasibility, already recorded.
                isempty(d.global_error) && isempty(d.lm_error) || error("Optimizer exception in $(c.name); see fit audit")
            end
            for (mode,c) in (("global",compared.control),("gcv",compared.selected))
                append!(outputs[mode],export_rows(file,c,data))
            end
            push!(selections,Dict("file"=>file,"N"=>string(n),"control"=>compared.control.name,"selected"=>compared.selected.name,
                "control_gcv"=>F.fmt(compared.control.fit.result.gcv),"selected_gcv"=>F.fmt(compared.selected.fit.result.gcv)))
            println("[$idx/$(length(files))] $file N=$n global=$(compared.control.name) expanded=$(compared.selected.name)")
            flush(stdout)
        catch err
            reason = replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
            push!(failures,Dict("file"=>file,"stage"=>stage,"reason"=>reason))
            println(stderr,"FAILED $file $stage: $reason"); flush(stderr)
        end
    end
    for (suffix,rows) in ((".fits.tsv",fits),(".parameters.tsv",parameters),(".selections.tsv",selections),(".failures.tsv",failures))
        isempty(rows) || write_table(out*suffix,sort(collect(keys(first(rows)))),rows)
    end
    isempty(failures) || error("$(length(failures)) files failed; no incomplete feature arms emitted")
    for mode in ("global","gcv"); write_table(out*".$mode.tsv",vcat(header,AXIS_FIELDS),outputs[mode]); end
    write_table(out,vcat(header,AXIS_FIELDS),outputs["gcv"])
end

function main(args=ARGS)
    opts = parse_cli(args); opts===nothing || execute(opts)
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && LocalGaussianOrientation.main()

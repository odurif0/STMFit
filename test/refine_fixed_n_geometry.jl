#!/usr/bin/env julia
# Matched local optimization at saved N; no count sweep or benchmark inputs.
module FixedNGeometryRefinement
include(joinpath(@__DIR__, "profile_frozen_amplitudes.jl"))
using .FrozenAmplitudeProfile: G, VP
using .FrozenAmplitudeProfile.ReconstructedUnitAssignment: read_table, lobe_table, write_table
using LinearAlgebra, TOML, Printf
const F = FrozenAmplitudeProfile
const NL = VP.NL

function settings(path)
    raw = TOML.parsefile(path)
    all(haskey(raw, k) && isempty(raw[k]) for k in ("model", "selection", "preprocessing")) ||
        error("Diagnostic config must not override physical or assignment settings")
    r = raw["geometry_refinement"]
    Set(keys(r)) == Set(["algorithm", "objective"]) || error("Unexpected refinement settings")
    r["algorithm"] == "LN_BOBYQA" || error("Only the declared matched optimizer is supported")
    r["objective"] == "native_kappa_rss" || error("Native common objective required")
    return VP.read_options(raw["counting_variable_projection"]), Symbol(r["algorithm"])
end

function parse_cli(args)
    "--help" in args && return println("refine_fixed_n_geometry.jl --features TSV --data-dir DIR --config TOML --settings TOML --out NEW_PREFIX [--chunk I/N] [--dry-run]")
    # Reuse the strict Gaussian cache schema, forbidden-label checks and paths.
    opts = F.parse_cli(args)
    for suffix in (".joint.tsv", ".profiled.tsv", ".native.tsv", ".fits.tsv", ".parameters.tsv", ".trace.tsv", ".bootstrap.tsv")
        p = opts["--out"] * suffix
        (ispath(p) || islink(p)) && error("Output exists: $p")
    end
    settings(opts["--settings"])
    return opts
end

"The same native penalized objective and validity checks in both arms."
function evaluate(p, n, data, cfg, amp_min, amp_range)
    r = G.ChainModelResult(n=n, params=copy(p), success=true, amp_min=amp_min, amp_range=amp_range)
    VP.finalize_native!(r, data, cfg)
    factor = cfg.kappa_max > 0 && n > 1 ?
        1 + G.kappa_penalty(r.kappa_max_adj; kappa_max=cfg.kappa_max, weight=cfg.kappa_weight) : 1.0
    objective = r.rss * factor
    return (; result=r, objective, valid=r.valid && isfinite(objective))
end

"Shared fresh native circular→elliptical fit, selected by existing GCV at fixed N."
function initialize(n, data, circ, ell, options)
    rc = G._fit_chain_n(data.xs, data.ys, data.zimg, data.x, data.y, data.z, data.noise, n, data.axisctx, circ)
    VP.finalize_native!(rc, data, circ)
    re = G.ChainModelResult(n=n, success=false, reason="circular initialization failed")
    if rc.success
        localell = deepcopy(ell)
        localell.skip_global = true; localell.multistart = 1
        localell.max_iter = options.native_elliptical_maxiter
        pe = VP.circular_to_elliptical(rc.params, n, localell)
        re = G._fit_chain_n(data.xs, data.ys, data.zimg, data.x, data.y, data.z, data.noise, n,
            data.axisctx, localell; starts=1, warm_start=pe)
        VP.finalize_native!(re, data, localell)
    end
    best, source = F.Extractor._best_for_n(F.Extractor._best_by_n([re], "gcv"),
        F.Extractor._best_by_n([rc], "gcv"), n, "gcv")
    best === nothing && error("No valid native initialization at saved N=$n; circ=$(rc.reason), ell=$(re.reason)")
    return best, source, (source == "circ" ? circ : ell), (circ=rc, ell=re)
end

function refine(initial, data, cfg, mode, options, algorithm; trace=(row -> nothing))
    mode in ("joint", "profiled") || error("Unknown refinement mode")
    p0 = copy(initial.params); n = initial.n
    lo, hi = VP.native_raw_bounds(n, cfg)
    all(isfinite, p0) && all(lo .<= p0 .<= hi) || error("Initial raw vector outside native bounds")
    indices = mode == "joint" ? collect(eachindex(p0)) : collect(VP.linear_layout(n, cfg).geometry)
    initial_eval = evaluate(p0, n, data, cfg, initial.amp_min, initial.amp_range)
    initial_eval.valid || error("Invalid shared initial fit")
    best = initial_eval; best_profile = nothing
    evaluations = 0; rejected = 0; improvements = 0
    objective = function(q, gradient)
        isempty(gradient) || error("Derivative-free optimizer required")
        evaluations += 1
        p = copy(p0); p[indices] .= q
        prof = nothing
        if mode == "profiled"
            prof = VP.fixed_geometry_profile(p, n, data.x, data.y, data.z, data.axisctx, cfg;
                amp_min=initial.amp_min, amp_range=initial.amp_range, options)
            prof.success || error("Linear profile failure: $(prof.linear.status), mapping=$(prof.mapping.passed)")
            p = prof.params
        end
        candidate = evaluate(p, n, data, cfg, initial.amp_min, initial.amp_range)
        if candidate.valid && (candidate.objective < best.objective ||
                (mode == "profiled" && best_profile === nothing && candidate.objective == best.objective))
            best = candidate; best_profile = prof; improvements += 1
        end
        rejected += !candidate.valid
        trace((evaluation=evaluations, rss=candidate.result.rss, objective=candidate.objective,
            valid=candidate.valid, reason=candidate.result.reason, overlap=candidate.result.overlap,
            residual_peak_snr=candidate.result.residual_peak_snr, kappa=candidate.result.kappa_max_adj,
            linear_kkt=prof === nothing ? NaN : prof.linear.kkt_violation,
            linear_tolerance=prof === nothing ? NaN : prof.linear.kkt_tolerance))
        return candidate.valid ? candidate.objective : Inf
    end
    # This first call counts against the evaluation ceiling. It compiles each
    # callback before NLopt's wall-clock timer starts; setup time is reported
    # separately, not represented as part of that 30-second optimization cap.
    setup_start = time_ns()
    objective(p0[indices], Float64[])
    setup_s = (time_ns() - setup_start) / 1e9
    opt = NL.Opt(algorithm, length(indices))
    opt.lower_bounds = lo[indices]; opt.upper_bounds = hi[indices]
    opt.maxeval = options.outer_maxeval - 1
    opt.maxtime = options.outer_maxtime_s
    opt.xtol_rel = options.outer_xtol_rel; opt.ftol_rel = options.outer_ftol_rel
    opt.min_objective = objective
    started = time_ns()
    _, _, ret = NL.optimize(opt, p0[indices])
    elapsed_s = (time_ns() - started) / 1e9
    status = string(ret)
    status in ("SUCCESS", "FTOL_REACHED", "XTOL_REACHED", "STOPVAL_REACHED", "MAXEVAL_REACHED", "MAXTIME_REACHED", "ROUNDOFF_LIMITED") ||
        error("Unexpected optimizer termination: $status")
    evaluations <= options.outer_maxeval || error("Evaluation ceiling exceeded")
    best.valid && best.objective <= initial_eval.objective || error("No valid nonworsening objective")
    # A profiled solution must really have had its inner problem solved.
    mode == "profiled" && best_profile === nothing && error("No admissible profiled point; no silent native fallback")
    return (; best, profile=best_profile, status, evaluations, rejected, improvements, setup_s, elapsed_s,
        converged=status in ("SUCCESS", "FTOL_REACHED", "XTOL_REACHED", "STOPVAL_REACHED"),
        initial_objective=initial_eval.objective, initial_rss=initial_eval.result.rss)
end

"Use the unchanged extractor's serialization in all freshly computed tables."
function feature_rows(file, result, data, cfg, source)
    b0, feats, ts, us, sp, sq = G._decode_chain(result.params, result.n, data.axisctx, cfg;
        amp_min=result.amp_min, amp_range=result.amp_range)
    amax = maximum(f.amplitude for f in feats)
    ax, ay = data.axisctx.axis; ox, oy = data.axisctx.origin
    return [Dict("file"=>file, "N"=>string(result.n), "lobe"=>string(i), "source"=>source,
        "amplitude"=>@sprintf("%.8e", f.amplitude), "x_nm"=>@sprintf("%.6f", f.x_nm), "y_nm"=>@sprintf("%.6f", f.y_nm),
        "t_nm"=>@sprintf("%.6f", ts[i]), "u_nm"=>@sprintf("%.6f", us[i]),
        "sigma_parallel_nm"=>@sprintf("%.6f", sp[i]), "sigma_perp_nm"=>@sprintf("%.6f", sq[i]),
        "spacing_prev_nm"=>i == 1 ? "NA" : @sprintf("%.6f", ts[i]-ts[i-1]),
        "amp_rel"=>@sprintf("%.6f", f.amplitude / max(amax, 1e-30)), "skew_ratio"=>@sprintf("%.6f", f.skew_ratio),
        "axis_x"=>@sprintf("%.8f", ax), "axis_y"=>@sprintf("%.8f", ay),
        "origin_x_nm"=>@sprintf("%.6f", ox), "origin_y_nm"=>@sprintf("%.6f", oy),
        "baseline"=>@sprintf("%.8e", b0), "tilt_x"=>@sprintf("%.8e", cfg.chain_tilted_baseline ? result.params[2] : 0.),
        "tilt_y"=>@sprintf("%.8e", cfg.chain_tilted_baseline ? result.params[3] : 0.),
        "gcv"=>@sprintf("%.8e", result.gcv)) for (i, f) in enumerate(feats)]
end

function execute(opts)
    BLAS.set_num_threads(1)
    options, algorithm = settings(opts["--settings"])
    options.outer_maxeval >= 3 || error("At least three evaluations required")
    raw = TOML.parsefile(opts["--config"])
    get(raw["model"], "selection_criterion", "") == "gcv" || error("GCV required")
    get(raw["model"], "cv_method", "") == "gcv" || error("No hidden kfold refits allowed")
    startswith(get(raw["model"], "selection_policy", ""), "adaptive_support") && error("Adaptive support is outside this comparison")
    header, base = lobe_table(opts["--features"]; required=F.REQUIRED)
    files = sort(unique(first.(collect(keys(base)))))
    chunk = parse.(Int, split(get(opts, "--chunk", "1/1"), '/'))
    files = [f for (i,f) in enumerate(files) if mod1(i,chunk[2]) == chunk[1]]
    isempty(files) && error("Empty chunk")
    chains = Dict(f=>[base[k] for k in sort(collect(keys(base))) if first(k)==f] for f in files)
    for f in files
        F.validate_chain(chains[f]); isfile(joinpath(opts["--data-dir"], f)) || error("Missing raw $f")
    end
    haskey(opts, "--dry-run") && return println("Matched geometry refinement: $(length(files)) files; metadata only, no pixels or outputs.")
    out = opts["--out"]; mkpath(out * ".fit_data")
    outputs = Dict(m=>Dict{String,String}[] for m in ("native", "joint", "profiled"))
    fits = Dict{String,String}[]; parameters = Dict{String,String}[]; bootstrap = Dict{String,String}[]; failures = Dict{String,String}[]
    open(out * ".trace.tsv", "w") do traceio
        println(traceio, "file\tmethod\tevaluation\trss\tobjective\tvalid\treason\toverlap\tresidual_peak_snr\tkappa\tlinear_kkt\tlinear_tolerance")
        for (idx, file) in enumerate(files)
            stage = "data"
            try
                n = length(chains[file])
                pcfg, ell, circ = F.Extractor._configs(raw["model"], raw["preprocessing"], dirname(out))
                pcfg.filepath = joinpath(opts["--data-dir"], file)
                img = G.read_sxm(pcfg.filepath)
                xs, ys, zimg, mask, x, y, z, noise = G._fused_roi_data(img, pcfg)
                full = G._weighted_roi_axis(x, y, z)
                xf, yf, zf, axis, keep, support = G._chain_fit_data(x, y, z, full, ell)
                data = (xs=xs, ys=ys, zimg=zimg, x=xf, y=yf, z=zf, zfull=z, noise=noise, axisctx=axis)
                for (key,v,precision) in (("axis_x",axis.axis[1],8), ("axis_y",axis.axis[2],8),
                    ("origin_x_nm",axis.origin[1],6), ("origin_y_nm",axis.origin[2],6))
                    s = precision == 8 ? @sprintf("%.8f",v) : @sprintf("%.6f",v)
                    s == chains[file][1][key] || error("Cache/input axis or origin mismatch")
                end
                open(joinpath(out * ".fit_data", file * ".tsv"), "w") do io
                    println(io,"x_nm\ty_nm\tz_nm")
                    for j in eachindex(zf); println(io, join(F.fmt.((xf[j],yf[j],zf[j])), '\t')); end
                end
                stage = "initialization"
                initial, source, cfg, boot = initialize(n, data, circ, ell, options)
                for kind in ("circ", "ell")
                    r = getproperty(boot, Symbol(kind))
                    push!(bootstrap, Dict("file"=>file, "N"=>string(n), "source"=>kind,
                        "success"=>string(r.success), "valid"=>string(r.valid), "gcv"=>F.fmt(r.gcv), "reason"=>r.reason))
                end
                lo, hi = VP.native_raw_bounds(n, cfg)
                for mode in ("native", "joint", "profiled")
                    stage = mode
                    outcome = mode == "native" ? nothing : refine(initial, data, cfg, mode, options, algorithm;
                        trace=row->println(traceio, join((file, mode, values(row)...), '\t')))
                    result = mode == "native" ? initial : outcome.best.result
                    append!(outputs[mode], feature_rows(file, result, data, cfg, source))
                    audit = Dict("file"=>file, "N"=>string(n), "method"=>mode, "source"=>source,
                        "rss"=>F.fmt(result.rss), "gcv"=>F.fmt(result.gcv), "p_full"=>string(G._chain_nparams(n,cfg)),
                        "n_data"=>string(length(zf)), "amp_min"=>F.fmt(result.amp_min), "amp_range"=>F.fmt(result.amp_range),
                        "amp_max_data"=>F.fmt(max(maximum(zimg),G.EPS)), "noise"=>F.fmt(noise),
                        "support_tmin"=>F.fmt(axis.tmin), "support_tmax"=>F.fmt(axis.tmax),
                        "axis_x"=>F.fmt(axis.axis[1]), "axis_y"=>F.fmt(axis.axis[2]),
                        "origin_x_nm"=>F.fmt(axis.origin[1]), "origin_y_nm"=>F.fmt(axis.origin[2]),
                        "overlap"=>F.fmt(result.overlap), "endpoint_overrun_nm"=>F.fmt(result.endpoint_overrun_nm),
                        "residual_peak_snr"=>F.fmt(result.residual_peak_snr), "kappa"=>F.fmt(result.kappa_max_adj),
                        "valid"=>string(result.valid), "reason"=>result.reason,
                        "objective"=>F.fmt(evaluate(result.params,n,data,cfg,result.amp_min,result.amp_range).objective))
                    for key in (:status,:evaluations,:rejected,:improvements,:setup_s,:elapsed_s,:converged,:initial_objective,:initial_rss)
                        audit[string(key)] = outcome === nothing ? "NA" : string(getproperty(outcome,key))
                    end
                    prof = outcome === nothing ? nothing : outcome.profile
                    audit["linear_kkt"] = prof === nothing ? "NA" : F.fmt(prof.linear.kkt_violation)
                    audit["linear_tolerance"] = prof === nothing ? "NA" : F.fmt(prof.linear.kkt_tolerance)
                    push!(fits,audit)
                    for j in eachindex(result.params)
                        push!(parameters, Dict("file"=>file,"method"=>mode,"parameter"=>string(j),
                            "value"=>F.fmt(result.params[j]),"lower"=>F.fmt(lo[j]),"upper"=>F.fmt(hi[j])))
                    end
                    println("[$idx/$(length(files))] $file N=$n $mode rss=$(result.rss)",
                        outcome === nothing ? " shared initialization" : " $(outcome.status) evaluations=$(outcome.evaluations) rejected=$(outcome.rejected)")
                    flush(stdout); flush(traceio)
                end
            catch err
                reason = replace(sprint(showerror,err),'\n'=>' ','\t'=>' ')
                push!(failures,Dict("file"=>file,"stage"=>stage,"reason"=>reason))
                println(stderr,"FAILED $file $stage: $reason"); flush(stderr)
            end
        end
    end
    for (suffix,rows) in ((".fits.tsv",fits),(".parameters.tsv",parameters),(".bootstrap.tsv",bootstrap),(".failures.tsv",failures))
        isempty(rows) || write_table(out * suffix, sort(collect(keys(first(rows)))), rows)
    end
    isempty(failures) || error("$(length(failures)) files failed; no incomplete arm emitted")
    for mode in ("native", "joint", "profiled")
        write_table(out * ".$mode.tsv", header, outputs[mode])
    end
    # Main shard product is profiled; joint/native companions are merged separately.
    write_table(out, header, outputs["profiled"])
end

function main(args=ARGS)
    opts = parse_cli(args); opts === nothing || execute(opts)
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && FixedNGeometryRefinement.main()

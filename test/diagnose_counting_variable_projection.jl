#!/usr/bin/env julia
# Diagnostic only. One file and an explicit unlabeled candidate-N list.
# Never called by production; no selection, unit labels, or expected counts.
module CountingVariableProjectionCLI

using GaussianFit2D
using TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__, "lib", "counting_variable_projection.jl"))
using .CountingVariableProjection
const G = GaussianFit2D
module Extractor
include(joinpath(@__DIR__, "extract_lobe_features.jl"))
end

const HELP = """
Bounded opt-in Gaussian COUNTING variable-projection diagnostic (Julia 1.13).

julia --startup-file=no --threads=1 --project=. test/diagnose_counting_variable_projection.jl \\
  --file ONE.sxm --config MODEL.toml --diagnostic-config DIAGNOSTIC.toml \\
  --selected-summary ORIGINAL.tsv --candidate-ns N1,N2 --outdir NEW_DIRECTORY [--dry-run]

All six path/list arguments are required. Candidate Ns must be explicitly chosen
from unlabeled saved outputs; this command does not infer an N range or select N.
The original selected summary restores adaptive support, not an expected count.
DIAGNOSTIC.toml requires all fields of [counting_variable_projection] (see tests).
No output may exist, including an empty directory or dangling symlink.

For every N and circular/elliptical geometry, fits.tsv keeps native,
fixed_geometry_profile and profile_refinement rows, including failures.
parameters.tsv stores actual raw coordinates/bounds; lobes.tsv stores decoded
physical parameters; fit_data.tsv stores the exact common data/support.
metadata.toml records effective settings and input/source hashes.

This is profile refinement of a native result, NOT an equal-budget replacement.
The native circular/global budget comes from MODEL.toml. Elliptical native fitting
uses the existing circular→elliptical initialization, skip_global=true, and the
explicit native_elliptical_maxiter. Outer LN_BOBYQA profiles raw RSS only, matching
the native final local objective, not the global kappa-penalized objective.
Native optimizer convergence is not exposed by core and is reported unknown.
Budget-limited outer iterates remain present but are not called converged.
GCV keeps every original parameter. No N_selected output is generated.
--dry-run checks metadata/options/output boundaries only; it reads no SXM pixels.
"""

function parse_cli(args)
    "--help" in args && return nothing
    opts = Dict{String,String}()
    dry = false
    allowed = Set(["file", "config", "diagnostic-config", "selected-summary", "candidate-ns", "outdir"])
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--dry-run"
            dry && error("Duplicate --dry-run")
            dry = true; i += 1; continue
        end
        startswith(arg, "--") || error("Unexpected positional argument: $arg")
        pieces = split(arg[3:end], '='; limit=2)
        key = pieces[1]
        key in allowed || error("Unknown option --$key")
        haskey(opts, key) && error("Duplicate --$key")
        if length(pieces) == 2
            value = pieces[2]
        else
            i += 1
            i <= length(args) || error("Missing value for --$key")
            value = args[i]
        end
        isempty(value) && error("Empty --$key")
        opts[key] = value
        i += 1
    end
    Set(keys(opts)) == allowed || error("Required options: " * join(sort!(collect(allowed)), ", "))
    return (args=opts, dry_run=dry)
end

_sha(path) = bytes2hex(sha256(read(path)))
_toml_value(v) = v === nothing ? "nothing" : v isa Symbol ? String(v) : v
_snapshot(cfg) = Dict(string(k) => _toml_value(getfield(cfg, k)) for k in fieldnames(typeof(cfg)))

function prepare(args)
    parsed = parse_cli(args)
    parsed === nothing && return nothing
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 is required; got $VERSION")
    Threads.nthreads() == 1 || error("Use --threads=1 for this bounded diagnostic")
    a = parsed.args
    paths = Dict(k => abspath(a[k]) for k in ("file", "config", "diagnostic-config", "selected-summary", "outdir"))
    for k in ("file", "config", "diagnostic-config", "selected-summary")
        isfile(paths[k]) || error("Not a file (--$k): $(paths[k])")
    end
    ispath(paths["outdir"]) || islink(paths["outdir"]) ? error("Output already exists: $(paths["outdir"])") : nothing
    ns = [parse(Int, strip(s)) for s in split(a["candidate-ns"], ','; keepempty=true)]
    !isempty(ns) && all(>(0), ns) && length(unique(ns)) == length(ns) || error("Candidate Ns must be unique positive integers")
    cfg = TOML.parsefile(paths["config"])
    all(haskey(cfg, k) for k in ("model", "selection", "preprocessing")) || error("Config requires model, selection, preprocessing tables")
    diagnostic = TOML.parsefile(paths["diagnostic-config"])
    haskey(diagnostic, "counting_variable_projection") || error("Missing [counting_variable_projection]")
    options = read_options(diagnostic["counting_variable_projection"])
    contexts = Extractor._read_selected_context(paths["selected-summary"];
        selection_policy=get(cfg["model"], "selection_policy", ""))
    file = basename(paths["file"])
    haskey(contexts, file) || error("Original selected summary has no row for $file")
    context = contexts[file]
    pcfg, ell, circ = Extractor._configs(cfg["model"], cfg["preprocessing"], paths["outdir"];
                                        selected_context=context)
    pcfg.filepath = paths["file"]
    G._chain_peak_profile(ell) == :gaussian || error("Only Gaussian COUNTING is supported")
    ell.cv_method == "gcv" || error("cv_method must be gcv; no hidden kfold refits")
    for key in ("global_maxtime", "global_maxiter", "max_iter", "multistart")
        haskey(cfg["model"], key) || error("Explicit model.$key budget required")
        v = cfg["model"][key]
        v isa Real && !(v isa Bool) && isfinite(v) && v > 0 || error("Positive finite model.$key required")
    end
    # Same refinement schedule as the extractor, with its otherwise hard-coded
    # iteration budget made explicit in the separate diagnostic settings.
    ell.skip_global = true
    ell.multistart = 1
    ell.max_iter = options.native_elliptical_maxiter
    return (paths=paths, ns=ns, file=file, config=cfg, diagnostic=diagnostic,
            options=options, context=context, pcfg=pcfg, ell=ell, circ=circ,
            dry_run=parsed.dry_run)
end

"Load exactly the native fused/robust ROI and selected-support construction."
function load_data(prepared)
    img = G.read_sxm(prepared.paths["file"])
    pcfg, cfg = prepared.pcfg, prepared.circ
    has_bwd = any(c -> lowercase(c.name) == lowercase(pcfg.roi_channel) && lowercase(c.direction) == "bwd", img.channels)
    fused = cfg.fuse_z_bwd && has_bwd
    xs, ys, zimg, mask, x, y, z, noise = fused ? G._fused_roi_data(img, pcfg) : G._robust_roi_data(img, pcfg)
    full = G._weighted_roi_axis(x, y, z)
    xf, yf, zf, axis, keep, support = G._chain_fit_data(x, y, z, full, cfg)
    return (xs=xs, ys=ys, zimg=zimg, x=xf, y=yf, z=zf, zfull=z, noise=noise,
            axisctx=axis, axisctx_full=full, support=support, fit_keep=keep, fused=fused)
end

const FIT_COLUMNS = ["file", "N", "geometry", "method", "initialization", "status", "success", "valid", "reason",
    "rss", "gcv", "p_full", "p_linear", "p_outer", "n_data", "elapsed_s", "optimizer_status", "optimizer_converged",
    "outer_evaluations", "outer_failed_evaluations", "linear_status", "linear_converged", "linear_iterations",
    "linear_kkt_violation", "linear_kkt_tolerance", "linear_residual_norm", "linear_rank_free", "linear_active_lower", "linear_active_upper",
    "raw_active_lower", "raw_active_upper", "mapping_passed", "mapping_max_abs_error", "mapping_tolerance",
    "initial_mapping_passed", "initial_mapping_max_abs_error", "amp_min", "amp_range", "amplitude_lower", "amplitude_upper",
    "baseline_lower", "baseline_upper", "noise", "support_tmin_nm", "support_tmax_nm", "overlap", "endpoint_overrun_nm",
    "residual_peak_snr", "kappa_max_adj", "bic", "aicc", "error"]
const PARAM_COLUMNS = ["file", "N", "geometry", "method", "parameter", "value", "lower", "upper"]
const LOBE_COLUMNS = ["file", "N", "geometry", "method", "lobe", "amplitude", "x_nm", "y_nm", "t_nm", "u_nm",
    "sigma_parallel_nm", "sigma_perp_nm", "skew_ratio", "axis_x", "axis_y", "origin_x_nm", "origin_y_nm",
    "baseline", "tilt_x", "tilt_y", "gcv", "valid", "reason"]
_cell(x) = replace(string(x), '\t' => ' ', '\n' => ' ', '\r' => ' ')
_write(io, values) = (println(io, join(_cell.(values), '\t')); flush(io))

function write_parameters(io, prep, n, geometry, method, p, cfg)
    lo, hi = native_raw_bounds(n, cfg)
    for i in eachindex(p)
        _write(io, (prep.file, n, geometry, method, i, p[i], lo[i], hi[i]))
    end
end

function write_fit!(ios, prep, data, n, geometry, method, cfg;
                    result=nothing, elapsed_s=NaN, profile=nothing, outer=nothing, error="")
    row = Dict{String,Any}(k => "NA" for k in FIT_COLUMNS)
    row["file"] = prep.file; row["N"] = n; row["geometry"] = geometry; row["method"] = method
    row["initialization"] = method == "native" ? (geometry == "circ" ? "native_seed" : "native_circular_result") : "native_result"
    row["status"] = "failed"; row["success"] = false; row["valid"] = false; row["error"] = error
    row["reason"] = isempty(error) ? "no_result" : error
    row["elapsed_s"] = elapsed_s; row["p_full"] = G._chain_nparams(n, cfg)
    row["p_linear"] = length(linear_layout(n, cfg).names)
    row["p_outer"] = row["p_full"]-row["p_linear"]
    if data !== nothing
        row["n_data"] = length(data.z); row["noise"] = data.noise
        row["support_tmin_nm"] = data.axisctx.tmin; row["support_tmax_nm"] = data.axisctx.tmax
    end
    row["optimizer_status"] = method == "native" ? "unknown_native_not_exposed" : "not_run"
    row["optimizer_converged"] = method == "native" ? "unknown" : false
    if profile !== nothing
        lin, map = profile.linear, profile.mapping
        row["linear_status"] = lin.status; row["linear_converged"] = lin.converged
        row["linear_iterations"] = lin.iterations; row["linear_kkt_violation"] = lin.kkt_violation
        row["linear_kkt_tolerance"] = lin.kkt_tolerance; row["linear_residual_norm"] = lin.residual_norm
        row["linear_rank_free"] = lin.rank_free
        names = linear_layout(n, cfg).names
        row["linear_active_lower"] = join(names[lin.active_lower], ';')
        row["linear_active_upper"] = join(names[lin.active_upper], ';')
        row["mapping_passed"] = map.passed; row["mapping_max_abs_error"] = map.max_abs_error
        row["mapping_tolerance"] = map.tolerance
        row["initial_mapping_passed"] = profile.initial_mapping.passed
        row["initial_mapping_max_abs_error"] = profile.initial_mapping.max_abs_error
        row["amplitude_lower"] = profile.lower[end]; row["amplitude_upper"] = profile.upper[end]
        row["baseline_lower"] = profile.lower[1]; row["baseline_upper"] = profile.upper[1]
        if method == "fixed_geometry_profile"
            row["optimizer_status"] = lin.status; row["optimizer_converged"] = lin.converged
        end
    end
    if outer !== nothing
        row["optimizer_status"] = outer.status; row["optimizer_converged"] = outer.converged
        row["outer_evaluations"] = outer.evaluations; row["outer_failed_evaluations"] = outer.failed_evaluations
        row["error"] = outer.error
    end
    if result !== nothing
        r = result
        row["success"] = r.success; row["valid"] = r.valid; row["reason"] = r.reason
        row["status"] = !r.success ? "failed" : r.valid ? "ok" : "invalid"
        for key in ("rss", "gcv", "amp_min", "amp_range", "overlap", "endpoint_overrun_nm", "residual_peak_snr", "kappa_max_adj", "bic", "aicc")
            row[key] = getfield(r, Symbol(key))
        end
        if !isempty(r.params)
            lo, hi = native_raw_bounds(n, cfg)
            row["raw_active_lower"] = join(findall(abs.(r.params-lo) .<= prep.options.linear_bound_atol), ';')
            row["raw_active_upper"] = join(findall(abs.(hi-r.params) .<= prep.options.linear_bound_atol), ';')
            write_parameters(ios.params, prep, n, geometry, method, r.params, cfg)
            if data !== nothing
                kw = (amp_min=r.amp_min, amp_range=r.amp_range)
                if profile === nothing
                    A = design_matrix(data.x, data.y, r.params, n, data.axisctx, cfg; kw...)
                    coeff = linear_coefficients(r.params, n, data.axisctx, cfg; kw...)
                    map = mapping_check(A, coeff, data.x, data.y, r.params, n, data.axisctx, cfg; kw..., options=prep.options)
                    row["mapping_passed"] = map.passed; row["mapping_max_abs_error"] = map.max_abs_error; row["mapping_tolerance"] = map.tolerance
                    lower, upper = CountingVariableProjection._linear_bounds(r.params, n, data.axisctx, cfg; kw...)
                    row["amplitude_lower"] = lower[end]; row["amplitude_upper"] = upper[end]
                    row["baseline_lower"] = lower[1]; row["baseline_upper"] = upper[1]
                end
                b0, feats, ts, us, spars, sperps = G._decode_chain(r.params, n, data.axisctx, cfg; kw...)
                for i in eachindex(feats)
                    f = feats[i]
                    _write(ios.lobes, (prep.file, n, geometry, method, i, f.amplitude, f.x_nm, f.y_nm, ts[i], us[i], spars[i], sperps[i],
                        f.skew_ratio, data.axisctx.axis..., data.axisctx.origin..., b0,
                        cfg.chain_tilted_baseline ? r.params[2] : 0.0, cfg.chain_tilted_baseline ? r.params[3] : 0.0, r.gcv, r.valid, r.reason))
                end
            end
        end
    end
    _write(ios.fits, [row[k] for k in FIT_COLUMNS])
    return row
end

function process_geometry!(ios, prep, data, n, geometry, cfg; warm_start=nothing, unavailable="")
    r = nothing
    start = time_ns()
    try
        isempty(unavailable) || error(unavailable)
        r = G._fit_chain_n(data.xs, data.ys, data.zimg, data.x, data.y, data.z, data.noise, n, data.axisctx, cfg;
                          starts=cfg.multistart, warm_start=warm_start)
        finalize_native!(r, data, cfg)
        write_fit!(ios, prep, data, n, geometry, "native", cfg; result=r, elapsed_s=(time_ns()-start)/1e9)
    catch err
        write_fit!(ios, prep, data, n, geometry, "native", cfg; elapsed_s=(time_ns()-start)/1e9, error=sprint(showerror, err))
        r = nothing
    end
    if r === nothing || !r.success
        reason = r === nothing ? "native_result_unavailable" : "native_result_failed: " * r.reason
        for method in ("fixed_geometry_profile", "profile_refinement")
            write_fit!(ios, prep, data, n, geometry, method, cfg; error=reason)
        end
        return r, false
    end
    first = nothing
    kw = (amp_min=r.amp_min, amp_range=r.amp_range, options=prep.options)
    start = time_ns()
    try
        first = fixed_geometry_profile(r.params, n, data.x, data.y, data.z, data.axisctx, cfg; kw...)
        prof = G.ChainModelResult(n=n, params=first.params, success=first.success, amp_min=r.amp_min, amp_range=r.amp_range,
                                  rss=first.rss, gcv=full_gcv(first.rss, length(data.z), n, cfg),
                                  reason=first.success ? "" : "profile_linear_or_mapping_failure")
        finalize_native!(prof, data, cfg)
        write_fit!(ios, prep, data, n, geometry, "fixed_geometry_profile", cfg; result=prof, profile=first, elapsed_s=(time_ns()-start)/1e9)
    catch err
        write_fit!(ios, prep, data, n, geometry, "fixed_geometry_profile", cfg; elapsed_s=(time_ns()-start)/1e9, error=sprint(showerror, err))
    end
    start = time_ns()
    try
        first !== nothing && first.success || error("initial_profile_unavailable")
        outer = refine_profile(r.params, n, data.x, data.y, data.z, data.axisctx, cfg; kw..., initial=first)
        best = outer.best
        prof = G.ChainModelResult(n=n, params=best.params, success=best.success, amp_min=r.amp_min, amp_range=r.amp_range)
        finalize_native!(prof, data, cfg)
        write_fit!(ios, prep, data, n, geometry, "profile_refinement", cfg; result=prof, profile=best, outer=outer, elapsed_s=(time_ns()-start)/1e9)
        return r, best.success && outer.status != "exception"
    catch err
        write_fit!(ios, prep, data, n, geometry, "profile_refinement", cfg; elapsed_s=(time_ns()-start)/1e9, error=sprint(showerror, err))
        return r, false
    end
end

function run(prep)
    out = prep.paths["outdir"]
    # mkdir, not mkpath(out): reject output races/collisions a second time.
    mkpath(dirname(out)); mkdir(out)
    BLAS.set_num_threads(1)
    metadata = Dict{String,Any}("status" => "running", "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "comparison" => "profile refinement of native result; not equal-budget replacement",
        "outer_objective" => "raw RSS; same as native local final objective, not penalized global stage",
        "candidate_ns" => prep.ns, "candidate_source" => "explicit CLI list from unlabeled saved summary; no inference",
        "selected_context" => Dict("n" => prep.context.n, "use_rescue" => prep.context.use_rescue),
        "paths" => prep.paths, "input_sha256" => Dict(k => _sha(prep.paths[k]) for k in ("file", "config", "diagnostic-config", "selected-summary")),
        "source_sha256" => Dict("diagnostic_cli" => _sha(@__FILE__), "diagnostic_lib" => _sha(joinpath(@__DIR__, "lib", "counting_variable_projection.jl")),
                                "native_core" => _sha(joinpath(dirname(pathof(G)), "core.jl")), "extractor" => _sha(joinpath(@__DIR__, "extract_lobe_features.jl"))),
        "diagnostic_options" => _snapshot(prep.options), "effective_circular" => _snapshot(prep.circ),
        "effective_elliptical" => _snapshot(prep.ell), "effective_preprocessing" => _snapshot(prep.pcfg))
    save_metadata() = open(io -> TOML.print(io, metadata; sorted=true), joinpath(out, "metadata.toml"), "w")
    save_metadata()
    ios = (fits=open(joinpath(out, "fits.tsv"), "w"), params=open(joinpath(out, "parameters.tsv"), "w"), lobes=open(joinpath(out, "lobes.tsv"), "w"))
    _write(ios.fits, FIT_COLUMNS); _write(ios.params, PARAM_COLUMNS); _write(ios.lobes, LOBE_COLUMNS)
    complete = true
    try
        data = try
            load_data(prep)
        catch err
            reason = "data_preparation_failed: " * sprint(showerror, err)
            for n in prep.ns, (geometry, cfg) in (("circ", prep.circ), ("ell", prep.ell)), method in ("native", "fixed_geometry_profile", "profile_refinement")
                write_fit!(ios, prep, nothing, n, geometry, method, cfg; error=reason)
            end
            metadata["error"] = reason
            rethrow()
        end
        metadata["support"] = Dict(string(k) => v for (k,v) in pairs(data.support))
        metadata["fused_backward"] = data.fused
        metadata["axis"] = Dict("origin" => collect(data.axisctx.origin), "axis" => data.axisctx.axis,
                                "perp" => data.axisctx.perp, "tmin" => data.axisctx.tmin, "tmax" => data.axisctx.tmax)
        save_metadata()
        open(joinpath(out, "fit_data.tsv"), "w") do io
            _write(io, ("pixel", "x_nm", "y_nm", "z"))
            for i in eachindex(data.z); _write(io, (i, data.x[i], data.y[i], data.z[i])); end
        end
        for n in prep.ns
            println("Diagnostic N=$n (no selection)"); flush(stdout)
            circ, ok_c = process_geometry!(ios, prep, data, n, "circ", prep.circ)
            warm = circ !== nothing && circ.success ? circular_to_elliptical(circ.params, n, prep.ell) : nothing
            warm === nothing || write_parameters(ios.params, prep, n, "ell", "circular_warm_start", warm, prep.ell)
            _, ok_e = process_geometry!(ios, prep, data, n, "ell", prep.ell; warm_start=warm,
                unavailable=warm === nothing ? "native_circular_result_unavailable" : "")
            complete &= ok_c && ok_e
        end
        metadata["status"] = complete ? "complete" : "complete_with_failed_rows"
    catch err
        metadata["status"] = "failed"
        metadata["error"] = sprint(showerror, err)
        complete = false
    finally
        foreach(close, values(ios))
        save_metadata()
    end
    println("Saved diagnostic only: $out; status=$(metadata["status"])")
    return complete ? 0 : 1
end

function main(args=ARGS)
    prep = prepare(args)
    prep === nothing && (println(HELP); return 0)
    if prep.dry_run
        println("Metadata-only dry-run: file=$(prep.file), candidate_ns=$(prep.ns), selected_support_rescue=$(prep.context.use_rescue)")
        println("No SXM data read, fit, selection, or output creation.")
        return 0
    end
    return run(prep)
end
end # module

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(CountingVariableProjectionCLI.main())
end

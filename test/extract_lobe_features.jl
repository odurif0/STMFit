#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# extract_lobe_features.jl — Run the 2D chain Gaussian fit and extract per-lobe
# features for unit assignment diagnostics.
#
# This script re-runs the same fit as batch_full.jl (circular sweep →
# circ→ell refinement → GCV selection) for each SXM file and extracts the
# per-lobe Gaussian parameters of the selected model. It also saves the
# chain axis direction for downstream aligned feature and patch diagnostics.
#
# The output TSV has one row per lobe with columns:
#   file, N, lobe, amplitude, x_nm, y_nm, t_nm, u_nm,
#   sigma_parallel_nm, sigma_perp_nm, spacing_prev_nm, amp_rel,
#   axis_x, axis_y, origin_x_nm, origin_y_nm, gcv, source
#
# Usage:
#   STMFIT_DATA_DIR=/path/to/data julia -t 4 --project=. \
#       test/extract_lobe_features.jl \
#       --config config/chitosan.toml \
#       --out results/unit_separability/lobe_features.tsv
#
#   # Specific files only
#   STMFIT_DATA_DIR=/path/to/data julia --project=. \
#       test/extract_lobe_features.jl \
#       --config config/chitosan.toml \
#       --files 240817_002.sxm,240817_003.sxm \
#       --out results/unit_separability/lobe_features.tsv
#
#   # Chunk for HPC
#   STMFIT_DATA_DIR=/path/to/data julia --project=. \
#       test/extract_lobe_features.jl --config config/chitosan.toml \
#       --chunk 1/4 --out results/unit_separability/lobe_features_chunk01.tsv
# ──────────────────────────────────────────────────────────────────────────────

using GaussianFit2D
using GaussianFit2D: ChainSweepConfig, ChainModelResult
using Printf
using TOML
using Dates

const DEFAULT_DATA_DIR = get(ENV, "STMFIT_DATA_DIR", "")
const DEFAULT_CONFIG = "config/chitosan.toml"
const DEFAULT_OUT = "results/unit_separability/lobe_features.tsv"

function _parse_cli(args)
    config_file = DEFAULT_CONFIG
    data_dir = DEFAULT_DATA_DIR
    out_tsv = DEFAULT_OUT
    files::Union{Nothing,Vector{String}} = nothing
    chunk_idx = 1
    chunk_total = 1
    exclude_file = ""
    selected_summary = ""
    manifest_file = ""
    primary_only = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            config_file = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--data-dir"
            data_dir = args[i+1]; i += 2
        elseif startswith(arg, "--data-dir=")
            data_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            files = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--files=")
            files = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--chunk"
            parts = split(args[i+1], "/")
            chunk_idx = parse(Int, parts[1]); chunk_total = parse(Int, parts[2]); i += 2
        elseif startswith(arg, "--chunk=")
            parts = split(split(arg, "=", limit=2)[2], "/")
            chunk_idx = parse(Int, parts[1]); chunk_total = parse(Int, parts[2]); i += 1
        elseif arg == "--exclude-from"
            exclude_file = args[i+1]; i += 2
        elseif startswith(arg, "--exclude-from=")
            exclude_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--selected-summary"
            selected_summary = args[i+1]; i += 2
        elseif startswith(arg, "--selected-summary=")
            selected_summary = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--manifest"
            manifest_file = args[i+1]; i += 2
        elseif startswith(arg, "--manifest=")
            manifest_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--primary-only"
            primary_only = true; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: STMFIT_DATA_DIR=/data julia -t 4 --project=. test/extract_lobe_features.jl [options]

            Options:
              --config PATH      Calibration TOML [$(DEFAULT_CONFIG)]
              --data-dir PATH    SXM data directory [\$STMFIT_DATA_DIR]
              --out PATH         Output TSV [$(DEFAULT_OUT)]
              --files LIST       Comma-separated SXM filenames (default: all in data dir)
              --chunk I/N        Process chunk I of N (for HPC parallelism)
              --exclude-from F   Exclude files listed in F (one per line, # comments)
              --selected-summary S  Optional batch summary TSV. If provided,
                                    features are extracted at N_selected from
                                    the summary instead of raw best GCV.
              --manifest M      Optional benchmark manifest for quality filtering
              --primary-only    With --manifest, keep only clean/clean_target files
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    isempty(data_dir) && error("No data directory: set STMFIT_DATA_DIR or pass --data-dir")
    isdir(data_dir) || error("Data directory not found: $data_dir")
    isfile(config_file) || error("Config not found: $config_file")

    return config_file, data_dir, out_tsv, files, chunk_idx, chunk_total,
           exclude_file, selected_summary, manifest_file, primary_only
end

function _read_exclude_set(path::String)
    isempty(path) && return Set{String}()
    isfile(path) || return Set{String}()
    s = Set{String}()
    for line in readlines(path)
        strip(line) == "" && continue
        startswith(strip(line), '#') && continue
        push!(s, basename(strip(line)))
    end
    return s
end

function _read_selected_n(path::String)
    selected = Dict{String,Int}()
    isempty(path) && return selected
    isfile(path) || error("Selected summary not found: $path")
    lines = readlines(path)
    isempty(lines) && return selected
    header = split(lines[1], '\t'; keepempty=true)
    file_idx = findfirst(==("filepath"), header)
    n_idx = findfirst(==("N_selected"), header)
    file_idx === nothing && error("Summary missing filepath column: $path")
    n_idx === nothing && error("Summary missing N_selected column: $path")
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, '\t'; keepempty=true)
        length(vals) < max(file_idx, n_idx) && continue
        n = tryparse(Int, vals[n_idx])
        n === nothing && continue
        selected[basename(vals[file_idx])] = n
    end
    return selected
end

function _read_manifest_quality(path::String)
    quality = Dict{String,String}()
    isempty(path) && return quality
    isfile(path) || error("Manifest not found: $path")
    m = TOML.parsefile(path)
    files = get(m, "files", Dict{String,Any}())
    for (file, info_any) in files
        info = info_any isa Dict ? info_any : Dict{String,Any}()
        quality[file] = String(get(info, "quality", "clean"))
    end
    return quality
end

_is_primary_quality(q::String) = !(q in ("poor_quality", "excluded"))

_score(r, criterion) = lowercase(String(criterion)) == "gcv" ? r.gcv :
                       lowercase(String(criterion)) == "aicc" ? r.aicc :
                       lowercase(String(criterion)) == "cv" ? r.cv_nll_mean : r.bic
_valid(r, criterion) = r.success && r.valid && isfinite(_score(r, criterion))

function _best_by_n(results, criterion)
    by_n = Dict{Int,Any}()
    for r in results
        _valid(r, criterion) || continue
        if !haskey(by_n, r.n) || _score(r, criterion) < _score(by_n[r.n], criterion)
            by_n[r.n] = r
        end
    end
    return by_n
end

function _effective_best(by_ell, by_circ, criterion)
    best_n, best_r, best_source, best_s = 0, nothing, "NA", Inf
    for n in sort(unique(vcat(collect(keys(by_ell)), collect(keys(by_circ)))))
        r_ell = get(by_ell, n, nothing); r_circ = get(by_circ, n, nothing)
        s_ell = r_ell === nothing ? Inf : _score(r_ell, criterion)
        s_circ = r_circ === nothing ? Inf : _score(r_circ, criterion)
        if min(s_ell, s_circ) < best_s
            best_n, best_s = n, min(s_ell, s_circ)
            best_r, best_source = s_ell <= s_circ ? (r_ell, "ell") : (r_circ, "circ")
        end
    end
    return best_n, best_r, best_source
end

function _best_for_n(by_ell, by_circ, n::Int, criterion)
    r_ell = get(by_ell, n, nothing)
    r_circ = get(by_circ, n, nothing)
    r_ell === nothing && r_circ === nothing && return nothing, "NA"
    s_ell = r_ell === nothing ? Inf : _score(r_ell, criterion)
    s_circ = r_circ === nothing ? Inf : _score(r_circ, criterion)
    return s_ell <= s_circ ? (r_ell, "ell") : (r_circ, "circ")
end

# Refine circular results to elliptical (same as visual_inspect_chitosan_cases.jl)
function _refine_circ_to_ell(results_circ, img, pcfg, ccfg_ell, ctx_circ)
    refined = ChainModelResult[]
    xs, ys, zimg, _, x, y, z, noise = GaussianFit2D._fused_roi_data(img, pcfg)
    xfit, yfit, zfit, ac_fit, _, _ = GaussianFit2D._chain_fit_data(x, y, z, ctx_circ.axisctx_full, ccfg_ell)
    n_eff = max(10, length(zfit) ÷ 9)
    ccfg_refine = deepcopy(ccfg_ell)
    ccfg_refine.skip_global = true; ccfg_refine.max_iter = 50; ccfg_refine.multistart = 1
    for r_c in results_circ
        r_c.success || continue
        n = r_c.n
        n_prefix = 1 + (ccfg_refine.chain_tilted_baseline ? 2 : 0)
        split_idx = n_prefix + n + GaussianFit2D._chain_spacing_param_count(n, ccfg_refine) + n
        n_sigma_types = GaussianFit2D._chain_sigma_param_count(n, ccfg_refine)
        sigma_block = r_c.params[(split_idx+1):(split_idx+n_sigma_types)]
        tail_start = split_idx + n_sigma_types + 1
        tail = tail_start <= length(r_c.params) ? r_c.params[tail_start:end] : Float64[]
        p_init = vcat(r_c.params[1:split_idx], sigma_block, sigma_block, tail)
        try
            r = GaussianFit2D._fit_chain_n(xs, ys, zimg, xfit, yfit, zfit, noise, n, ac_fit, ccfg_refine; starts=1, warm_start=p_init)
            if r.success
                pred = GaussianFit2D._chain_model_values(xfit, yfit, r.params, n, ac_fit, ccfg_refine; amp_min=r.amp_min, amp_range=r.amp_range)
                GaussianFit2D._finalize_chain_result!(r, zfit, pred, noise, n, n_eff, z, xs, ys, zimg, xfit, yfit, ac_fit, ccfg_refine)
                push!(refined, r)
            end
        catch err
            @warn "elliptical refinement failed" n exception=(err, catch_backtrace())
        end
    end
    return refined
end

function _configs(model, preproc, output_dir)
    pcfg = PatternConfig(filepath="", channel="Z", direction="fwd",
        stride=get(preproc, "stride", 1),
        flatten=get(preproc, "flatten", "plane+rows"),
        smooth_radius_px=get(preproc, "smooth_radius_px", 1),
        output_dir=output_dir, no_plot=true)
    ccfg = ChainSweepConfig(n_min=Int(get(model, "n_min", 2)), n_max=Int(get(model, "n_max", 14)),
        spacing_min_nm=model["spacing_min_nm"], spacing_max_nm=model["spacing_max_nm"],
        fit_width_nm=model["fit_width_nm"],
        support_noise_k=model["support_noise_k"],
        support_padding_nm=model["support_padding_nm"],
        support_min_length_nm=get(model, "support_min_length_nm", 1.0),
        support_baseline_quantile=get(model, "support_baseline_quantile", 0.10),
        max_overlap=model["max_overlap"],
        global_maxtime=model["global_maxtime"], global_maxiter=model["global_maxiter"],
        max_iter=get(model, "max_iter", 300), multistart=get(model, "multistart", 1),
        cv_folds=get(model, "cv_folds", 5), cv_method=get(model, "cv_method", "gcv"),
        selection_criterion=get(model, "selection_criterion", "gcv"),
        sigma_parallel_min_nm=model["sigma_parallel_min_nm"],
        sigma_parallel_max_nm=model["sigma_parallel_max_nm"],
        sigma_perp_min_nm=get(model, "sigma_perp_min_nm", model["sigma_parallel_min_nm"]),
        sigma_perp_max_nm=get(model, "sigma_perp_max_nm", model["sigma_parallel_max_nm"]),
        kappa_max=get(model, "kappa_max", 10.0),
        kappa_weight=get(model, "kappa_weight", 1.0),
        min_amplitude_fraction=get(model, "min_amplitude_fraction", 0.3),
        shared_sigma_types=get(model, "shared_sigma_types", 0),
        chain_spacing_model=get(model, "chain_spacing_model", "free"),
        chain_tilted_baseline=get(model, "chain_tilted_baseline", true),
        peak_profile=Symbol(String(get(model, "peak_profile", "gaussian"))),
        skew_ratio_max=Float64(get(model, "skew_ratio_max", 2.0)),
        intelligent_sweep=true, fuse_z_bwd=true)
    ccfg_circ = deepcopy(ccfg); ccfg_circ.chain_circular_sigmas = true
    return pcfg, ccfg, ccfg_circ
end

function main()
    config_file, data_dir, out_tsv, files_opt, chunk_idx, chunk_total,
    exclude_file, selected_summary, manifest_file, primary_only = _parse_cli(ARGS)
    exclude_set = _read_exclude_set(exclude_file)
    selected_n = _read_selected_n(selected_summary)
    manifest_quality = _read_manifest_quality(manifest_file)

    # Determine file list
    if files_opt !== nothing
        all_files = String.(strip.(files_opt))
    elseif !isempty(selected_n)
        all_files = sort(collect(keys(selected_n)))
    else
        all_files = sort([f for f in readdir(data_dir) if endswith(lowercase(f), ".sxm")])
    end
    all_files = [f for f in all_files if !(basename(f) in exclude_set)]
    if primary_only
        all_files = [f for f in all_files if _is_primary_quality(get(manifest_quality, basename(f), "clean"))]
    end

    # Chunk selection
    if chunk_total > 1
        chunk_files = [f for (i, f) in enumerate(all_files) if mod1(i, chunk_total) == chunk_idx]
    else
        chunk_files = all_files
    end

    println("Extracting lobe features for ", length(chunk_files), " files (chunk ", chunk_idx, "/", chunk_total, ")")

    # Load config
    cfg = TOML.parsefile(config_file)
    model, preproc = cfg["model"], get(cfg, "preprocessing", Dict{String,Any}())
    criterion = get(model, "selection_criterion", "gcv")

    mkpath(dirname(out_tsv))

    open(out_tsv, "w") do io
        println(io, join([
            "file", "N", "lobe", "amplitude", "x_nm", "y_nm", "t_nm", "u_nm",
            "sigma_parallel_nm", "sigma_perp_nm", "spacing_prev_nm", "amp_rel",
            "skew_ratio",
            "axis_x", "axis_y", "origin_x_nm", "origin_y_nm",
            "baseline", "tilt_x", "tilt_y", "gcv", "source"
        ], '\t'))

        for (i_file, fn) in enumerate(chunk_files)
            fp = joinpath(data_dir, fn)
            isfile(fp) || (@warn "File not found: $fn"; continue)

            try
                img = read_sxm(fp)
                pcfg, ccfg, ccfg_circ = _configs(model, preproc, dirname(out_tsv))
                pcfg.filepath = fp
                if haskey(selected_n, fn)
                    # When a batch summary supplies N_selected, this extractor is
                    # a fixed-N refit/feature export. Avoid the full N sweep,
                    # especially for split-width diagnostics with extra params.
                    ccfg.n_min = selected_n[fn]
                    ccfg.n_max = selected_n[fn]
                    ccfg_circ.n_min = selected_n[fn]
                    ccfg_circ.n_max = selected_n[fn]
                end

                results_circ, _, ctx = chain_gaussian_sweep(img, pcfg, ccfg_circ)
                results_ell = _refine_circ_to_ell(results_circ, img, pcfg, ccfg, ctx)

                by_ell = _best_by_n(results_ell, criterion)
                by_circ = _best_by_n(results_circ, criterion)
                if haskey(selected_n, fn)
                    best_n = selected_n[fn]
                    best_r, best_source = _best_for_n(by_ell, by_circ, best_n, criterion)
                else
                    best_n, best_r, best_source = _effective_best(by_ell, by_circ, criterion)
                end

                if best_r === nothing
                    @warn "No valid fit for $fn"
                    continue
                end

                ccfg_sel = best_source == "circ" ? ccfg_circ : ccfg
                b0, feats, ts, us, spars, sperps = GaussianFit2D._decode_chain(
                    best_r.params, best_r.n, ctx.axisctx, ccfg_sel;
                    amp_min=best_r.amp_min, amp_range=best_r.amp_range)
                tilt_x = ccfg_sel.chain_tilted_baseline && length(best_r.params) >= 3 ? best_r.params[2] : 0.0
                tilt_y = ccfg_sel.chain_tilted_baseline && length(best_r.params) >= 3 ? best_r.params[3] : 0.0

                amax = maximum([f.amplitude for f in feats])
                ax, ay = ctx.axisctx.axis
                ox, oy = ctx.axisctx.origin

                for k in 1:length(feats)
                    spacing_prev = k == 1 ? NaN : ts[k] - ts[k-1]
                    println(io, join([
                        fn, string(best_r.n), string(k),
                        @sprintf("%.8e", feats[k].amplitude),
                        @sprintf("%.6f", feats[k].x_nm),
                        @sprintf("%.6f", feats[k].y_nm),
                        @sprintf("%.6f", ts[k]),
                        @sprintf("%.6f", us[k]),
                        @sprintf("%.6f", spars[k]),
                        @sprintf("%.6f", sperps[k]),
                        isnan(spacing_prev) ? "NA" : @sprintf("%.6f", spacing_prev),
                        @sprintf("%.6f", feats[k].amplitude / max(amax, 1e-30)),
                        @sprintf("%.6f", feats[k].skew_ratio),
                        @sprintf("%.8f", ax), @sprintf("%.8f", ay),
                        @sprintf("%.6f", ox), @sprintf("%.6f", oy),
                        @sprintf("%.8e", b0),
                        @sprintf("%.8e", tilt_x),
                        @sprintf("%.8e", tilt_y),
                        @sprintf("%.8e", best_r.gcv),
                        best_source
                    ], '\t'))
                end

                @printf("[%d/%d] %-24s  N=%d  source=%s  GCV=%.3e  lobes=%d\n",
                        i_file, length(chunk_files), fn, best_r.n, best_source, best_r.gcv, length(feats))
                flush(stdout)
                flush(io)
            catch e
                @warn "Failed" fn reason=sprint(showerror, e)
            end
        end
    end

    println("\nWrote: ", out_tsv)
    println("Done: ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
end

main()

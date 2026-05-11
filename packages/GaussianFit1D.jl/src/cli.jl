"""
Command-line interface for Multi-Gaussian fitting.

Usage:
    julia --project run.jl [options]
    julia --project run.jl data.txt --min-spacing 0.4
    julia --project run.jl --load cache.jld2 --plot-n 5
"""

using ArgParse

export parse_args, config_from_args, main_cli

function parse_command_line()
    s = ArgParseSettings(
        description="Multi-Gaussian fitting of 1D line profiles with physical constraints.",
        epilog="""
Examples:
  julia --project run.jl                                     # default config
  julia --project run.jl data.txt                            # custom file
  julia --project run.jl data.txt --min-spacing 0.4          # custom spacing
  julia --project run.jl --load cache.jld2 --plot-n 5        # replot from cache
""",
    )

    @add_arg_table! s begin
        "filepath"
            help = "Data file to fit (default: from CONFIG)"
            default = ""
            required = false

        "--min-spacing"
            help = "Min spacing between adjacent peaks (nm); default: 0.4"
            arg_type = Float64
            default = nothing

        "--max-spacing"
            help = "Max spacing between adjacent peaks (nm); default: 0.675"
            arg_type = Float64
            default = nothing

        "--fwhm-min"
            help = "Min FWHM of individual patterns (nm); default: 0.45"
            arg_type = Float64
            default = nothing

        "--fwhm-max"
            help = "Max FWHM of individual patterns (nm); default: 1.2"
            arg_type = Float64
            default = nothing

        "--max-overlap"
            help = "Max allowed overlap between adjacent widest peaks; derives an effective min spacing (default: 0.6)"
            arg_type = Float64
            default = nothing

        "--amplitude-min-fraction"
            help = "Min amplitude as fraction of max(y) (when amplitude-min not set; default: 0.7)"
            arg_type = Float64
            default = nothing

        "--amplitude-min"
            help = "Minimum peak amplitude (absolute); default: auto from fraction"
            arg_type = Float64
            default = nothing

        "--amplitude-max"
            help = "Maximum peak amplitude; default: no upper bound"
            arg_type = Float64
            default = nothing

        "--nlopt-algorithm"
            help = "NLopt algorithm (GN_ESCH, GN_ISRES, GN_DIRECT_L, GN_CRS2_LM)"
            arg_type = String
            default = nothing

        "--edge-sigma-min"
            help = "Min edge margin in units of sigma_max (default: 1.0)"
            arg_type = Float64
            default = nothing

        "--edge-sigma-max"
            help = "Max edge margin in units of sigma_max (default: 4.0)"
            arg_type = Float64
            default = nothing

        "--x-unit"
            help = "Distance unit for display (e.g. nm, Angstrom)"
            arg_type = String
            default = nothing

        "--global-maxiter"
            help = "Global optimization max evaluations"
            arg_type = Int
            default = nothing

        "--global-maxtime"
            help = "Global optimization max time (seconds)"
            arg_type = Float64
            default = nothing

        "--plot-n"
            help = "Plot a specific n_peaks model instead of the BIC-best"
            arg_type = Int
            default = nothing

        "--load"
            help = "Load cached results from JLD2 file (skip fitting)"
            arg_type = String
            default = nothing

        "--no-show"
            help = "Don't display plots (batch mode)"
            action = :store_true

        "--output-dir", "-o"
            help = "Output directory for plots and data files"
            arg_type = String
            default = nothing

        "--asymmetric-edges"
            help = "Use asymmetric (split) Gaussians for edge peaks"
            action = :store_true

        "--use-log-amplitude"
            help = "Fit amplitudes in log-space A = exp(p) (default: linear)"
            action = :store_true
    end

    return s
end

function config_from_args(args)
    """Build a FitConfig from parsed CLI arguments."""
    cfg = deepcopy(DEFAULT_CONFIG)

    # Filepath: only override if explicitly provided (non-empty)
    fp = args["filepath"]
    if !isempty(fp)
        cfg.filepath = fp
    end

    for (arg_name, cfg_key) in [
        ("min-spacing", "min_spacing"),
        ("max-spacing", "max_spacing"),
        ("fwhm-min", "fwhm_min"),
        ("fwhm-max", "fwhm_max"),
        ("max-overlap", "max_overlap"),
        ("amplitude-min-fraction", "amplitude_min_fraction"),
        ("amplitude-min", "amplitude_min"),
        ("amplitude-max", "amplitude_max"),
        ("edge-sigma-min", "edge_sigma_min"),
        ("edge-sigma-max", "edge_sigma_max"),
        ("x-unit", "x_unit"),
        ("global-maxiter", "global_maxiter"),
        ("global-maxtime", "global_maxtime"),
        ("nlopt-algorithm", "nlopt_algorithm"),
        ("output-dir", "output_dir"),
    ]
        val = args[arg_name]
        if val !== nothing
            sym = Symbol(cfg_key)
            v = sym == :nlopt_algorithm && val isa String ? Symbol(val) : val
            setfield!(cfg, sym, v)
        end
    end
    if args["no-show"]
        cfg.no_show = true
    end
    if args["asymmetric-edges"]
        cfg.asymmetric_edges = false
    end
    if args["use-log-amplitude"]
        cfg.use_log_amplitude = true
    end
    return cfg
end

function main_cli()
    s = parse_command_line()
    args = parse_args(s)

    loaded_from_cache = args["load"] !== nothing
    fr = nothing

    try
        if loaded_from_cache
            # ---- Load mode ----
            cache_file = args["load"]
            x, y, all_results, cfg = load_results(cache_file)

            # Apply CLI overrides on the loaded config
            for (arg_name, cfg_key) in [
                ("min-spacing", "min_spacing"),
                ("max-spacing", "max_spacing"),
                ("fwhm-min", "fwhm_min"),
                ("fwhm-max", "fwhm_max"),
                ("max-overlap", "max_overlap"),
                ("amplitude-min-fraction", "amplitude_min_fraction"),
                ("amplitude-min", "amplitude_min"),
                ("amplitude-max", "amplitude_max"),
                ("edge-sigma-min", "edge_sigma_min"),
                ("edge-sigma-max", "edge_sigma_max"),
                ("x-unit", "x_unit"),
                ("global-maxiter", "global_maxiter"),
                ("global-maxtime", "global_maxtime"),
                ("nlopt-algorithm", "nlopt_algorithm"),
                ("output-dir", "output_dir"),
            ]
                val = args[arg_name]
                if val !== nothing
                    sym = Symbol(cfg_key)
                    v = sym == :nlopt_algorithm && val isa String ? Symbol(val) : val
                    setfield!(cfg, sym, v)
                end
            end
            if args["no-show"]
                cfg.no_show = true
            end
            if args["asymmetric-edges"]
                cfg.asymmetric_edges = false
            end
            if args["use-log-amplitude"]
                cfg.use_log_amplitude = true
            end
            fp = args["filepath"]
            if !isempty(fp)
                cfg.filepath = fp
            end

            fr = FitRunResult(x, y, all_results, cfg, cache_file, nothing, nothing, nothing)
        else
            # ---- Run mode ----
            cfg = config_from_args(args)
            fr = run_fit(cfg; save_cache=true)
        end

    catch e
        println("ERROR: $e")
        return
    end

    best_r = best_result(fr)
    if best_r === nothing
        println("ERROR: No models could be fitted. Check your data and parameter settings.")
        return
    end

    # Select which model to plot
    plot_n = args["plot-n"]
    if plot_n !== nothing
        try
            plot_result = select_model_result(fr, plot_n)
        catch e
            println("ERROR: $e")
            return
        end
        println("\n  Plotting n_peaks=$plot_n (requested) instead of BIC-best ($(best_r.n_peaks))")
    else
        plot_result = best_r
    end

    print_summary(fr.all_results, plot_result, fr.cfg)

    export_results(fr.x, fr.y, fr.all_results, fr.cfg)
    plot_results(fr.x, fr.y, plot_result, fr.all_results, fr.cfg)
end

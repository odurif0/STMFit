function _print_info(img::SXMImage)
    println("File: $(img.filepath)")
    println("Pixels: $(img.width) × $(img.height)")
    println("Range: $(@sprintf("%.3f", img.range_nm[1])) nm × $(@sprintf("%.3f", img.range_nm[2])) nm")
    println("Channels:")
    for c in img.channels
        println("  - $(c.name) direction=$(c.direction) unit=$(c.unit) min=$(minimum(c.data)) max=$(maximum(c.data))")
    end
end

function _arg_settings()
    s = ArgParseSettings(description="STM molecular pattern recognition by preprocessing, blob detection and constrained Gaussian refinement.")
    @add_arg_table! s begin
        "filepath"
            help = "Fichier Nanonis .sxm"
            required = true
        "--info"
            action = :store_true
        "--channel"
            default = "Z"
        "--direction"
            default = "fwd"
        "--output-dir", "-o"
            default = "results"
        "--stride"
            arg_type = Int
            default = 2
        "--flatten"
            help = "none | plane | rows | plane+rows"
            default = "plane+rows"
        "--smooth-radius-px"
            arg_type = Int
            default = 1
        "--contrast"
            help = "auto | bright | dark"
            default = "auto"
        "--max-features"
            help = "Nombre maximal de candidats détectés/exportés"
            arg_type = Int
            default = 50
        "--max-fit-features"
            help = "Nombre maximal de candidats raffinés par fit non linéaire"
            arg_type = Int
            default = 8
        "--threshold-sigma"
            arg_type = Float64
            default = 2.5
        "--min-distance-px"
            arg_type = Int
            default = 10
        "--ignore-border-px"
            arg_type = Int
            default = 8
        "--chain-min-length"
            arg_type = Int
            default = 3
        "--chain-min-spacing-nm"
            arg_type = Float64
            default = 0.35
        "--chain-max-spacing-nm"
            arg_type = Float64
            default = 1.6
        "--chain-max-angle-deg"
            arg_type = Float64
            default = 45.0
        "--chain-max-spacing-cv"
            arg_type = Float64
            default = 0.45
        "--min-chain-score"
            arg_type = Float64
            default = 5.0
        "--max-chains"
            arg_type = Int
            default = 8
        "--max-path-branches"
            arg_type = Int
            default = 6
        "--repeat-spacing-nm"
            help = "Pas moléculaire attendu pour estimer le nombre d'unités (chitosane ≈ 0.52 nm)"
            arg_type = Float64
            default = 0.52
        "--end-extension-nm"
            help = "Correction d'enveloppe STM aux deux extrémités avant longueur/pas"
            arg_type = Float64
            default = 0.35
        "--no-fusion"
            help = "Ne pas fusionner les vues; détecter seulement le canal/direction choisi"
            action = :store_true
        "--fusion-channels"
            help = "Canaux utilisés pour la carte d'évidence fusionnée: Z, Current, Z,Current, ou all"
            default = "Z"
        "--min-sigma-nm"
            arg_type = Float64
            default = 0.05
        "--max-sigma-nm"
            arg_type = Float64
            default = 1.5
        "--initial-sigma-nm"
            arg_type = Float64
            default = nothing
        "--no-fit"
            action = :store_true
        "--no-plot"
            action = :store_true
        "--all-images"
            help = "Traiter tous les couples canal/direction du fichier SXM"
            action = :store_true
        "--inspect-images"
            help = "Créer seulement une planche raw/preprocessed pour toutes les images"
            action = :store_true
        "--threshold-sweep"
            help = "Balayage du seuil sans fit, format start:step:stop, ex. 2.0:0.25:4.0"
            default = ""
        "--selection-criterion"
            help = "Model selection criterion: bic | mdl | cv (default: bic)"
            default = "bic"
        "--chain-sweep"
            help = "Sélection N par chaîne 2D ordonnée"
            action = :store_true
        "--consensus"
            help = "Multi-channel consensus: fit on Z then on Current, compare N/spacing"
            action = :store_true
        "--batch"
            help = "Consensus batch processing on all .sxm files in a directory"
            arg_type = String
            default = ""
        "--chain-mode"
            help = "Chain fitting mode: sweep (default) | 1d-bootstrapped | direct"
            arg_type = String
            default = "sweep"
        "--chain-skip-global"
            help = "Skip NLopt global optimization when 1D bootstrap init is available (faster)"
            action = :store_true
        "--chain-init-centers-t"
            help = "Comma-separated list of t-axis centers (nm) from 1D fit, for 1d-bootstrapped mode"
            arg_type = String
            default = ""
        "--chain-init-amplitudes"
            help = "Comma-separated list of amplitudes from 1D fit, for 1d-bootstrapped mode"
            arg_type = String
            default = ""
        "--chain-init-sigma-parallel"
            help = "Mean sigma from 1D fit, for sigma_parallel initialization (nm)"
            arg_type = Float64
            default = NaN
        "--chain-init-sigma-perp"
            help = "Mean sigma from 1D fit, for sigma_perp initialization (nm)"
            arg_type = Float64
            default = NaN
        "--multistart"
            arg_type = Int
            default = 20
        "--cv-folds"
            arg_type = Int
            default = 5
        "--chain-n-min"
            arg_type = Int
            default = 2
        "--chain-n-max"
            arg_type = Int
            default = 14
        "--chain-spacing-min-nm"
            arg_type = Float64
            default = 0.35
        "--chain-spacing-max-nm"
            arg_type = Float64
            default = 0.75
        "--chain-lateral-max-nm"
            arg_type = Float64
            default = 0.35
        "--chain-fit-width-nm"
            help = "Demi-largeur du tube fitté autour de l'axe moléculaire; rend le fit comparable à une slice 1D"
            arg_type = Float64
            default = 0.45
        "--chain-support-baseline-quantile"
            help = "Quantile robuste de baseline du profil axial pour support auto"
            arg_type = Float64
            default = 0.10
        "--chain-support-threshold-fraction"
            help = "Seuil relatif du profil axial pour garder seulement le support longitudinal actif"
            arg_type = Float64
            default = 0.25
        "--chain-support-noise-k"
            help = "Seuil bruit du support auto: baseline + k*MAD; combiné au seuil contraste par max()"
            arg_type = Float64
            default = 2.5
        "--chain-support-padding-nm"
            arg_type = Float64
            default = 0.20
        "--chain-support-min-length-nm"
            help = "Longueur minimale d'un composant actif longitudinal"
            arg_type = Float64
            default = 1.0
        "--chain-t-min-nm"
            help = "Absolute minimum t coordinate along chain axis; use slide_metadata support_start_t_nm to match a 1D slice"
            arg_type = Float64
            default = NaN
        "--chain-t-max-nm"
            help = "Absolute maximum t coordinate along chain axis; use slide_metadata support_end_t_nm to match a 1D slice"
            arg_type = Float64
            default = NaN
        "--chain-sigma-parallel-min-nm"
            arg_type = Float64
            default = 0.12
        "--chain-sigma-parallel-max-nm"
            arg_type = Float64
            default = 0.45
        "--chain-sigma-perp-min-nm"
            arg_type = Float64
            default = 0.10
        "--chain-sigma-perp-max-nm"
            arg_type = Float64
            default = 0.55
        "--fuse-z-bwd"
            help = "Fuse Z fwd + Z bwd for SNR boost (default: true)"
            action = :store_true
            default = true
        "--no-fuse-z-bwd"
            help = "Disable Z fwd+bwd fusion, use Z fwd only"
            action = :store_true
            default = false
        "--global-maxtime"
            help = "NLopt global optimization max time per model (seconds)"
            arg_type = Float64
            default = 15.0
        "--global-maxiter"
            help = "NLopt global optimization max evaluations per model"
            arg_type = Int
            default = 5000
        "--global-tol"
            help = "NLopt relative tolerance"
            arg_type = Float64
            default = 1e-5
    end
    return s
end

function main_cli()
    args = parse_args(_arg_settings())
    img = read_sxm(args["filepath"])
    if args["info"]
        _print_info(img)
        return
    end
    cfg = PatternConfig(
        filepath=args["filepath"], channel=args["channel"], direction=args["direction"],
        output_dir=args["output-dir"], stride=args["stride"], flatten=args["flatten"],
        smooth_radius_px=args["smooth-radius-px"], contrast=args["contrast"],
        max_features=args["max-features"], max_fit_features=args["max-fit-features"],
        threshold_sigma=args["threshold-sigma"],
        min_distance_px=args["min-distance-px"], ignore_border_px=args["ignore-border-px"],
        fusion=!args["no-fusion"], chain_min_length=args["chain-min-length"],
        fusion_channels=args["fusion-channels"],
        chain_min_spacing_nm=args["chain-min-spacing-nm"],
        chain_max_spacing_nm=args["chain-max-spacing-nm"],
        chain_max_angle_deg=args["chain-max-angle-deg"],
        chain_max_spacing_cv=args["chain-max-spacing-cv"], min_chain_score=args["min-chain-score"],
        max_chains=args["max-chains"], max_path_branches=args["max-path-branches"],
        repeat_spacing_nm=args["repeat-spacing-nm"], end_extension_nm=args["end-extension-nm"],
        min_sigma_nm=args["min-sigma-nm"], max_sigma_nm=args["max-sigma-nm"],
        initial_sigma_nm=args["initial-sigma-nm"], no_fit=args["no-fit"], no_plot=args["no-plot"])
    cfg.all_images = args["all-images"]
    cfg.inspect_images = args["inspect-images"]

    if !isempty(args["threshold-sweep"])
        _run_threshold_sweep(img, cfg, args["threshold-sweep"])
        return
    end

    if args["chain-sweep"]
        chain_mode = lowercase(args["chain-mode"])
        # Parse optional 1D bootstrap parameters
        parse_csv(s) = isempty(s) ? Float64[] : parse.(Float64, split(s, ','))
        init_centers_t = parse_csv(args["chain-init-centers-t"])
        init_amps      = parse_csv(args["chain-init-amplitudes"])
        init_sigpar    = args["chain-init-sigma-parallel"]
        init_sigperp   = args["chain-init-sigma-perp"]

        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_baseline_quantile=args["chain-support-baseline-quantile"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                t_min_nm=args["chain-t-min-nm"],
                                t_max_nm=args["chain-t-max-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"],
                                init_centers_t=init_centers_t,
                                init_amplitudes=init_amps,
                                init_sigma_parallel=init_sigpar,
                                init_sigma_perp=init_sigperp,
                                skip_global=args["chain-skip-global"],
                                selection_criterion=args["selection-criterion"])
        if chain_mode in ("1d-bootstrapped", "direct")
            # Force n_min == n_max if bootstrapped with init data
            if !isempty(init_centers_t)
                n_boot = length(init_centers_t)
                ccfg.n_min = n_boot
                ccfg.n_max = n_boot
            end
            r, ctx = chain_direct_fit(img, cfg, ccfg)
            _write_chain_sweep([r], r, ctx, cfg, ccfg)
            cfg.no_plot || _plot_chain_sweep([r], r, ctx, cfg, ccfg)
            println("chain_direct_N: $(r.n)")
            println("chain_direct_valid: $(r.valid)")
            println("chain_direct_bic: $(r.bic)")
            println("chain_direct_bic: $(r.bic)")
            println("chain_direct_residual_peak_snr: $(r.residual_peak_snr)")
            println("chain_direct_train_nll: $(r.train_nll)")
            r.success && println("chain_direct_mean_spacing_nm: $(r.mean_spacing_nm)")
        else
            _run_chain_sweep(img, cfg, ccfg)
        end
        return
    end

    if args["consensus"]
        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_baseline_quantile=args["chain-support-baseline-quantile"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                t_min_nm=args["chain-t-min-nm"],
                                t_max_nm=args["chain-t-max-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"])
        consensus = fit_chain_consensus(img, cfg, ccfg)
        # Write Z results
        _write_chain_sweep(consensus.z.results, consensus.z.best, consensus.z.ctx, cfg, ccfg)
        cfg.no_plot || _plot_chain_sweep(consensus.z.results, consensus.z.best, consensus.z.ctx, cfg, ccfg)
        if consensus.current.best !== nothing
            println("consensus_z_N: $(consensus.z.best.n)")
            println("consensus_current_N: $(consensus.current.best.n)")
            println("consensus_agreement: $(consensus.agreement)")
            println("consensus_match: $(consensus.consensus)")
        end
        return
    end

    if !isempty(args["batch"])
        dir = args["batch"]
        files = filter(f -> endswith(lowercase(f), ".sxm"), readdir(dir; join=true))
        isempty(files) && error("No .sxm files found in $dir")
        ccfg = ChainSweepConfig(n_min=args["chain-n-min"], n_max=args["chain-n-max"],
                                multistart=args["multistart"], cv_folds=args["cv-folds"],
                                spacing_min_nm=args["chain-spacing-min-nm"],
                                spacing_max_nm=args["chain-spacing-max-nm"],
                                lateral_max_nm=args["chain-lateral-max-nm"],
                                fit_width_nm=args["chain-fit-width-nm"],
                                support_threshold_fraction=args["chain-support-threshold-fraction"],
                                support_noise_k=args["chain-support-noise-k"],
                                support_padding_nm=args["chain-support-padding-nm"],
                                support_min_length_nm=args["chain-support-min-length-nm"],
                                sigma_parallel_min_nm=args["chain-sigma-parallel-min-nm"],
                                sigma_parallel_max_nm=args["chain-sigma-parallel-max-nm"],
                                sigma_perp_min_nm=args["chain-sigma-perp-min-nm"],
                                sigma_perp_max_nm=args["chain-sigma-perp-max-nm"],
                                fuse_z_bwd=args["fuse-z-bwd"] && !args["no-fuse-z-bwd"],
                                global_maxtime=args["global-maxtime"],
                                global_maxiter=args["global-maxiter"],
                                global_tol=args["global-tol"])
        fit_chain_batch(files, cfg, ccfg; consensus=true)
        return
    end

    if cfg.inspect_images
        path = _plot_all_images(img, cfg)
        println("all_images_overview: $path")
        return
    end

    if cfg.all_images
        _run_all_images(img, cfg)
        return
    end

    result = fit_molecular_pattern(img, cfg)
    param_file = _write_parameters(result, cfg)
    plot_file = cfg.no_plot ? nothing : _plot_result(img, cfg, result)

    println("success: $(result.success)")
    println("estimated_repeats: $(result.estimated_repeats) (range $(result.estimated_repeat_range[1])-$(result.estimated_repeat_range[2]), ROI length $(@sprintf("%.3f", result.roi_length_nm)) nm)")
    println("chains: $(length(result.chains))")
    println("detected_features: $(length(result.raw_features))")
    println("fitted_features: $(length(result.features))")
    println("RSS: $(@sprintf("%.6g", result.rss))")
    println("R²: $(@sprintf("%.6f", result.r_squared))")
    println("BIC-like score: $(@sprintf("%.6g", result.bic))")
    println("parameters: $param_file")
    plot_file !== nothing && println("plot: $plot_file")
    for w in result.warnings
        println("warning: $w")
    end
end

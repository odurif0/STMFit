# Configuration Reference

## ChainSweepConfig (GaussianFit2D)

Full configuration for the 2D chain model sweep.

```julia
GaussianFit2D.ChainSweepConfig(
    # ── Sweep range ──
    n_min              = 2,       # Minimum N (safety bound)
    n_max              = 14,      # Maximum N (safety bound)
    intelligent_sweep  = true,    # Adaptive range from support length
    early_stop_patience = 3,      # Consecutive BIC increases before stop
    early_stop_dbic    = 100.0,   # BIC increase threshold for early stop

    # ── Physical constraints ──
    spacing_min_nm     = 0.35,    # Minimum inter-lobe spacing (nm)
    spacing_max_nm     = 0.75,    # Maximum inter-lobe spacing (nm)
    max_overlap        = 0.60,    # Maximum lobe overlap fraction
    fit_width_nm       = 0.15,    # Tube half-width around axis (nm)

    # ── Sigma bounds ──
    sigma_parallel_min_nm = 0.191, # Min axial sigma (FWHM 0.45 nm)
    sigma_parallel_max_nm = 0.509, # Max axial sigma (FWHM 1.20 nm)
    sigma_perp_min_nm   = 0.10,   # Min perpendicular sigma
    sigma_perp_max_nm   = 0.55,   # Max perpendicular sigma

    # ── Model variants ──
    chain_circular_sigmas = false, # σ∥=σ⟂ per lobe (simpler, more robust)
    chain_tilted_baseline = true,  # Linear tilt in 2D background

    # ── Optimization ──
    global_maxtime     = 10.0,    # NLopt timeout per N
    global_maxiter     = 10000,   # NLopt max iterations per N
    global_tol         = 1e-5,    # NLopt tolerance
    max_iter           = 300,     # LsqFit max iterations
    multistart         = 1,       # Number of random starts

    # ── Support detection ──
    support_threshold_fraction = 0.20,
    support_noise_k    = 2.5,
    support_padding_nm = 0.20,
    support_min_length_nm = 1.0,
    support_baseline_quantile = 0.10,

    # ── Penalties ──
    kappa_max          = 8.0,     # Condition number penalty threshold
    kappa_weight       = 1.0,     # Condition number penalty strength
    peak_profile       = :gaussian, # :gaussian (2D: only :gaussian supported)
    min_amplitude_fraction = 0.3, # Min lobe amplitude (fraction of max data)

    # ── Cross-validation ──
    cv_folds           = 5,       # Number of CV folds (kfold only)
    cv_method          = "gcv",   # "gcv" (analytical, free) | "kfold" (refit per fold)
    student_nu         = 4.0,     # Student-t degrees of freedom
    residual_peak_snr_threshold = 3.5,

    # ── Selection ──
    selection_criterion = "bic",  # "bic" | "aicc" | "cv"
)
```

## PatternConfig (GaussianFit2D)

Image preprocessing and blob detection configuration.

```julia
GaussianFit2D.PatternConfig(
    filepath   = "",          # SXM file path
    channel    = "Z",         # Channel name
    direction  = "fwd",       # Scan direction
    stride     = 1,           # Subsampling stride
    flatten    = "plane+rows",# Background flattening
    smooth_radius_px = 1,     # Preprocessing smoothing
    threshold_sigma = 2.5,    # Blob detection threshold
    min_distance_px = 10,     # Minimum blob separation
    fusion     = true,        # Fuse Z fwd+bwd channels
    fuse_z_bwd = true,        # Fuse Z forward+backward scans
)
```

## FitSlideConfig (STMMolecularFit)

1D slide profile fitting configuration.

```julia
STMMolecularFit.FitSlideConfig(
    min_spacing    = 0.35,    # Minimum peak spacing (nm)
    max_spacing    = 0.75,    # Maximum peak spacing (nm)
    fwhm_min       = 0.45,    # Minimum FWHM (nm)
    fwhm_max       = 1.20,    # Maximum FWHM (nm)
    max_overlap    = 0.60,    # Maximum peak overlap
    kappa_max      = 8.0,     # Condition number threshold
    peak_profile   = :gaussian,  # :gaussian | :lorentzian | :pseudo_voigt
    amplitude_min_fraction = 0.3,
    global_maxtime = 8.0,     # NLopt timeout (s)
    global_maxiter = 5000,    # NLopt max iterations
)
```

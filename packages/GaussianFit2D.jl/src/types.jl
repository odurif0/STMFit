struct SXMChannel
    name::String
    unit::String
    direction::String
    data::Matrix{Float64} # [y, x], in physical unit stored by Nanonis
end

struct SXMImage
    filepath::String
    header::Dict{String,String}
    width::Int
    height::Int
    range_nm::Tuple{Float64,Float64}
    offset_nm::Tuple{Float64,Float64}
    channels::Vector{SXMChannel}
end

struct MolecularFeature
    amplitude::Float64
    x_nm::Float64
    y_nm::Float64
    sigma_x_nm::Float64
    sigma_y_nm::Float64
    score::Float64
end

struct MolecularChain
    id::Int
    features::Vector{MolecularFeature}
    score::Float64
    mean_spacing_nm::Float64
    spacing_cv::Float64
    max_turn_angle_deg::Float64
end

Base.@kwdef mutable struct ChainSweepConfig
    n_min::Int = 2
    n_max::Int = 14
    multistart::Int = 20
    cv_folds::Int = 5
    rng_seed::Int = 4321
    student_nu::Float64 = 4.0
    spacing_min_nm::Float64 = 0.35
    spacing_max_nm::Float64 = 0.75
    lateral_max_nm::Float64 = 0.35
    fit_width_nm::Float64 = 0.45
    support_baseline_quantile::Float64 = 0.10
    support_threshold_fraction::Float64 = 0.25
    support_noise_k::Float64 = 2.5
    support_padding_nm::Float64 = 0.20
    support_min_length_nm::Float64 = 1.0
    t_min_nm::Float64 = NaN
    t_max_nm::Float64 = NaN
    sigma_parallel_min_nm::Float64 = 0.12
    sigma_parallel_max_nm::Float64 = 0.45
    sigma_perp_min_nm::Float64 = 0.10
    sigma_perp_max_nm::Float64 = 0.55
    # 1D bootstrap / external initialization
    init_centers_t::Vector{Float64} = Float64[]
    init_amplitudes::Vector{Float64} = Float64[]
    init_laterals::Vector{Float64} = Float64[]
    init_sigma_parallel::Float64 = NaN
    init_sigma_perp::Float64 = NaN
    skip_global::Bool = false
    boot_sweep_halfwidth::Int = 2
    intelligent_sweep::Bool = true
    early_stop_patience::Int = 3
    early_stop_dbic::Float64 = 100.0
    fuse_z_bwd::Bool = true
    residual_peak_snr_threshold::Float64 = 3.5
    max_overlap::Float64 = 0.60
    kappa_max::Float64 = 8.0      # condition-number penalty threshold (0 = disabled)
    kappa_weight::Float64 = 1.0   # penalty strength
    min_amplitude_fraction::Float64 = 0.3  # reject models with any peak amplitude < 30% of max (matches 1D)
    max_iter::Int = 300
    global_maxtime::Float64 = 15.0
    global_maxiter::Int = 5000
    global_tol::Float64 = 1e-5
    chain_circular_sigmas::Bool = false  # circular gaussians (spar=sperp per peak, fewer params)
    selection_criterion::String = "bic"  # model selection: "bic" | "aicc" | "cv"
end

Base.@kwdef mutable struct ChainModelResult
    n::Int = 0
    params::Vector{Float64} = Float64[]
    success::Bool = false
    train_nll::Float64 = Inf
    cv_nll_mean::Float64 = Inf
    cv_nll_std::Float64 = Inf
    bic::Float64 = Inf
    aicc::Float64 = Inf
    mdl::Float64 = Inf
    residual_peak_snr::Float64 = Inf
    mean_spacing_nm::Float64 = Inf
    spacing_cv::Float64 = Inf
    max_lateral_nm::Float64 = Inf
    sigma_parallel_nm::Float64 = Inf
    sigma_perp_nm::Float64 = Inf
    overlap::Float64 = Inf
    kappa_max_adj::Float64 = 1.0   # max adjacent condition number
    endpoint_overrun_nm::Float64 = Inf
    bound_like::Int = 0
    valid::Bool = false
    reason::String = ""
    rss::Float64 = Inf
    chi2_reduced::Float64 = Inf
    mad::Float64 = Inf
    amp_min::Float64 = NaN
    amp_range::Float64 = NaN
end

Base.@kwdef mutable struct PatternConfig
    filepath::String = ""
    channel::String = "Z"
    direction::String = "fwd"
    output_dir::String = "results"

    # Preprocessing
    stride::Int = 2
    flatten::String = "plane+rows" # none | plane | rows | plane+rows
    smooth_radius_px::Int = 1
    ignore_border_px::Int = 8

    # Detection
    contrast::String = "auto" # auto | bright | dark
    max_features::Int = 50
    max_fit_features::Int = 8
    min_features::Int = 1
    threshold_sigma::Float64 = 2.5
    min_distance_px::Int = 10

    # Chain validation. All image channels/directions are views of one measure;
    # fusion + chain validation is the main recognition mode.
    fusion::Bool = true
    fusion_channels::String = "Z" # comma-separated channel names, or "all"
    chain_min_length::Int = 3
    chain_min_spacing_nm::Float64 = 0.35
    chain_max_spacing_nm::Float64 = 1.6
    chain_max_angle_deg::Float64 = 45.0
    chain_max_spacing_cv::Float64 = 0.45
    min_chain_score::Float64 = 5.0
    max_chains::Int = 8
    max_path_branches::Int = 6

    # Molecule ROI / principal-axis mode. The views are measurements of the same
    # molecule, so first isolate the molecular object in topography, then search
    # features inside that object instead of chasing background chains.
    roi::Bool = true
    roi_channel::String = "Z"
    roi_threshold_fraction::Float64 = 0.35
    roi_dilate_px::Int = 10
    axis_profile::Bool = true
    axis_peak_sigma::Float64 = 1.5
    axis_min_peak_distance_nm::Float64 = 0.45
    repeat_spacing_nm::Float64 = 0.52
    end_extension_nm::Float64 = 0.35

    # Fit constraints
    min_sigma_nm::Float64 = 0.05
    max_sigma_nm::Float64 = 1.5
    initial_sigma_nm::Union{Nothing,Float64} = nothing
    max_iter::Int = 300
    no_fit::Bool = false
    no_plot::Bool = false
    all_images::Bool = false
    inspect_images::Bool = false
end

Base.@kwdef mutable struct PatternFitResult
    features::Vector{MolecularFeature} = MolecularFeature[]
    raw_features::Vector{MolecularFeature} = MolecularFeature[]
    chains::Vector{MolecularChain} = MolecularChain[]
    evidence_map::Matrix{Float64} = zeros(0, 0)
    roi_mask::BitMatrix = falses(0, 0)
    axis_peaks::Vector{MolecularFeature} = MolecularFeature[]
    evidence_unit::String = "SNR"
    roi_length_nm::Float64 = NaN
    estimated_repeats::Int = 0
    estimated_repeat_range::Tuple{Int,Int} = (0, 0)
    params_unconstrained::Vector{Float64} = Float64[]
    success::Bool = false
    rss::Float64 = Inf
    r_squared::Float64 = -Inf
    bic::Float64 = Inf
    warnings::Vector{String} = String[]
    parameter_file::Union{String,Nothing} = nothing
    plot_file::Union{String,Nothing} = nothing
end

"""
Type definitions for the multi-Gaussian fitting project.

Replaces the CONFIG `Dict{String,Any}` with typed structs.
"""

# ===========================================================================
# Constants
# ===========================================================================

const MGF_VERSION = "4.1.0"  # bump on breaking config/cache changes

# ===========================================================================
# FitConfig — immutable configuration
# ===========================================================================

Base.@kwdef mutable struct FitConfig
    # --- Input / Output ---
    filepath::String                = "240314_Cu100_039_Image.txt"
    output_dir::Union{String,Nothing} = nothing
    x_unit::String                  = "nm"
    y_unit::String                  = ""

    # --- Physical constraints ---
    max_spacing::Float64            = 0.675
    min_spacing::Float64            = 0.4
    fwhm_max::Float64               = 1.2
    fwhm_min::Float64               = 0.45
    max_overlap::Float64            = 0.6
    kappa_max::Float64              = 15.0   # condition-number penalty threshold (0 = disabled)
    kappa_weight::Float64           = 1.0    # penalty strength relative to RSS

    # --- Baseline ---
    offset_to_zero::Bool            = true

    # --- Amplitude bounds ---
    amplitude_min_fraction::Float64 = 0.3
    amplitude_min::Union{Float64,Nothing} = nothing
    amplitude_max::Union{Float64,Nothing} = nothing

    # --- Center bounds ---
    edge_sigma_min::Float64         = 1.0
    edge_sigma_max::Float64         = 4.0

    # --- Optimizer (NLopt global + LsqFit local) ---
    global_maxiter::Int             = 5000
    global_maxtime::Float64         = 15.0
    global_tol::Float64             = 1e-5
    nlopt_algorithm::Symbol         = :GN_DIRECT_L
    curve_fit_maxfev::Int           = 500

    # --- Model selection ---
    early_stop_patience::Int        = 3
    early_stop_dbic::Float64        = 100.0
    bic_competition_threshold::Float64 = 20.0
    student_nu::Float64             = 4.0
    use_student_bic::Bool           = true
    noise_estimate::Float64         = NaN

    # --- Amplitude parameterization ---
    use_log_amplitude::Bool         = false  # false=linear bounded [amp_min,amp_max], true=log-space A=exp(p)

    # --- Edge model ---
    asymmetric_edges::Bool          = false

    # --- Visualization ---
    fig_dpi::Int                    = 150
    fig_width::Int                  = 12
    fig_height::Int                 = 10
    fine_grid_points::Int           = 1000
    no_show::Bool                   = false
end

const DEFAULT_CONFIG = FitConfig()

# ===========================================================================
# FitResult — per-model result
# ===========================================================================

Base.@kwdef mutable struct FitResult
    n_peaks::Int                    = 0
    popt::Vector{Float64}           = Float64[]
    pcov::Union{Matrix{Float64},Nothing} = nothing
    perr::Vector{Float64}           = Float64[]
    y_fit::Vector{Float64}          = Float64[]
    success::Bool                   = false
    warnings::Vector{String}        = String[]
    plot_file::Union{String,Nothing} = nothing

    # Metrics
    bic::Float64                    = Inf
    student_bic::Float64            = Inf
    aic::Float64                    = Inf
    aicc::Float64                   = Inf
    r_squared::Float64              = -Inf
    chi2_red::Float64               = Inf
    dof::Int                        = 0
    rss::Float64                    = Inf
    n_params::Int                   = 0
    kappa_max_adj::Float64          = 1.0    # max adjacent condition number from fit

    # Model comparison
    competitive::Bool               = false
    delta_bic::Float64              = 0.0

    # Internal: params without y0 for reuse
    popt_inner::Vector{Float64}     = Float64[]
end

# ===========================================================================
# FitMetrics — lightweight metrics return
# ===========================================================================

struct FitMetrics
    bic::Float64
    student_bic::Float64
    aic::Float64
    aicc::Float64
    r_squared::Float64
    chi2_red::Float64
    dof::Int
    rss::Float64
    n_params::Int
end

# ===========================================================================
# FitRunResult — top-level result container
# ===========================================================================

mutable struct FitRunResult
    x::Vector{Float64}
    y::Vector{Float64}
    all_results::Vector{FitResult}
    cfg::FitConfig
    cache_file::Union{String,Nothing}
    export_file::Union{String,Nothing}
    benchmark_file::Union{String,Nothing}
    plot_files::Union{Vector{String},Nothing}
end

function FitRunResult(x, y, all_results, cfg)
    return FitRunResult(x, y, all_results, cfg, nothing, nothing, nothing, nothing)
end

# ===========================================================================
# Config helpers
# ===========================================================================

function build_config(overrides::AbstractDict{String}=Dict{String,Any}())
    """Return a FitConfig with non-nothing overrides applied."""
    cfg = FitConfig()  # fresh instance with all defaults
    for (k, v) in overrides
        v !== nothing || continue
        sym = Symbol(k)
        hasfield(FitConfig, sym) || continue
        if sym == :nlopt_algorithm && v isa String
            v = Symbol(v)
        end
        setfield!(cfg, sym, v)
    end
    return cfg
end



# ===========================================================================
# User data directory (cross-platform)
# ===========================================================================

function user_data_dir()
    # Linux:   ~/.multigaussianfit
    # macOS:   ~/Library/Application Support/GaussianFit1D
    # Windows: %LOCALAPPDATA%/GaussianFit1D
    if Sys.iswindows()
        return joinpath(get(ENV, "LOCALAPPDATA", joinpath(homedir(), "AppData", "Local")), "GaussianFit1D")
    elseif Sys.isapple()
        return joinpath(homedir(), "Library", "Application Support", "GaussianFit1D")
    else
        return joinpath(homedir(), ".multigaussianfit")
    end
end

_cache_dir()   = joinpath(user_data_dir(), "cache")
_results_dir() = joinpath(user_data_dir(), "results")

# ===========================================================================
# Path utility
# ===========================================================================

function output_path(cfg::FitConfig, suffix::String)
    base = splitext(basename(cfg.filepath))[1]
    out_dir = cfg.output_dir
    if out_dir === nothing || isempty(out_dir)
        out_dir = joinpath(_results_dir(), base)
    end
    mkpath(out_dir)
    return joinpath(out_dir, base * suffix)
end

# ===========================================================================
# API helpers
# ===========================================================================

function best_result(fr::FitRunResult)
    if isempty(fr.all_results)
        return nothing
    end
    use_sbic = fr.cfg.use_student_bic
    return use_sbic ? argmin(r -> r.student_bic, fr.all_results) : argmin(r -> r.bic, fr.all_results)
end

function result_for_n(fr::FitRunResult, n_peaks::Int)
    idx = findfirst(r -> r.n_peaks == n_peaks, fr.all_results)
    return idx !== nothing ? fr.all_results[idx] : nothing
end

function select_model_result(fr::FitRunResult, n_peaks::Union{Int,Nothing}=nothing)
    if n_peaks === nothing
        return best_result(fr)
    end
    r = result_for_n(fr, n_peaks)
    if r === nothing
        available = [res.n_peaks for res in fr.all_results]
        error("n_peaks=$n_peaks was not fitted. Available: $available")
    end
    return r
end

function update_model_rankings(all_results::Vector{FitResult}, cfg::FitConfig)
    if isempty(all_results)
        return
    end
    use_sbic = cfg.use_student_bic
    if use_sbic
        best_val = minimum(r.student_bic for r in all_results)
        for r in all_results
            r.competitive = (r.student_bic - best_val) <= cfg.bic_competition_threshold
            r.delta_bic = r.student_bic - best_val
        end
    else
        best_val = minimum(r.bic for r in all_results)
        for r in all_results
            r.competitive = (r.bic - best_val) <= cfg.bic_competition_threshold
            r.delta_bic = r.bic - best_val
        end
    end
end

# ===========================================================================
# Cache versioning
# ===========================================================================

function _config_hash(cfg::FitConfig)
    """Deterministic hash of the config fields that affect fit results.
    Used to detect stale caches."""
    # Only hash fields that affect the fit (not output/display settings)
    h = hash("MGF_VERSION=$(MGF_VERSION)")
    for fld in fieldnames(FitConfig)
        v = getfield(cfg, fld)
        if fld ∈ (:output_dir, :x_unit, :y_unit, :fig_dpi, :fig_width,
                   :fig_height, :fine_grid_points, :no_show)
            continue  # display/output only, don't affect fit
        end
        h = hash((fld, v), h)
    end
    return UInt64(h)
end

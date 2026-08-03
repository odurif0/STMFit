# Candidate-N paired-view geometry + residual-patch extraction
# Types and public helpers for Todo 2.

"""
    ViewData

A single preprocessed scan direction (fwd or bwd). `z` is the flattened image
on the (xs, ys) grid in nm, indexed `[y, x]`. `present=false` marks a missing
view (kept for audit, never decoded).
"""
struct ViewData
    label::String
    xs::Vector{Float64}
    ys::Vector{Float64}
    z::Matrix{Float64}
    noise::Float64
    present::Bool
end

"""
    ViewRecalibration

Per-view affine recalibration result: z_view ≈ a·model + b (OLS over the masked
finite pixels). `resid` is z_view - (a·model + b) on those pixels.
"""
struct ViewRecalibration
    label::String
    a::Float64
    b::Float64
    rss::Float64
    resid::Vector{Float64}
    n_pixels::Int
    valid::Bool
    reason::String
end

"""
    LobePatch

Decoded geometry for one lobe plus per-view aligned patches. `residual_patches`
are sampled from the recalibrated residual image (z - a·model - b); `raw_patches`
from z - model. Patch vectors are ordered to match `CandidateNView.patch_tu`.
"""
struct LobePatch
    index::Int
    x_nm::Float64
    y_nm::Float64
    t_nm::Float64
    u_nm::Float64
    amplitude::Float64
    sigma_parallel_nm::Float64
    sigma_perp_nm::Float64
    skew_ratio::Float64
    residual_patches::Dict{String,Vector{Float64}}
    raw_patches::Dict{String,Vector{Float64}}
end

"""
    CandidateNView

One retained candidate N with paired-view joint metrics and per-lobe patches.
`effective_n_views = view_count * effective_view_factor`.
"""
struct CandidateNView
    n::Int
    source_gcv::Float64
    success::Bool
    valid::Bool
    reason::String
    view_count::Int
    bwd_missing::Bool
    views::Vector{ViewRecalibration}
    joint_gcv::Float64
    joint_nrmse::Float64
    residual_corr::Float64
    effective_view_factor::Float64
    effective_n_views::Float64
    p_eff::Int
    nd_joint::Int
    lobes::Vector{LobePatch}
    patch_tu::Vector{NamedTuple{(:t, :u), Tuple{Float64, Float64}}}
end

struct SkippedCandidate
    n::Int
    reason::String
end

struct CandidateViewReport
    candidates::Vector{CandidateNView}
    skipped::Vector{SkippedCandidate}
    view_count::Int
    bwd_missing::Bool
    audit::Vector{String}
end

"""
    effective_view_factor(rho)

Effective independent-view multiplier: `1 / (1 + clamp(rho, 0, 0.95))`.
"""
function effective_view_factor(rho::Float64)
    isnan(rho) && return 1.0
    return 1.0 / (1.0 + clamp(rho, 0.0, 0.95))
end

"""
    build_view_data(img, channel, pcfg; direction) -> ViewData

Preprocess one scan direction from a real SXM image.
"""
function build_view_data(img::SXMImage, channel::String, pcfg::PatternConfig;
                         direction::String="fwd")::ViewData
    ch = get_channel(img, channel; direction=direction)
    xs, ys, _raw, z, _z_smooth, _unit, noise = preprocess_channel(img, ch, pcfg)
    return ViewData(direction, xs, ys, z, noise, true)
end

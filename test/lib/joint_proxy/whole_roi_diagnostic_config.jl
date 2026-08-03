# Whole-ROI common-template frozen-contrast diagnostic — registration config.
#
# Diagnostic only. Parses and validates
# config/joint_proxy_whole_roi_diagnostic.toml into a typed immutable record
# consumed by the chemistry-blind global registration (Todo 3) and the frozen
# type-contrast scorer (Todo 4). It never enters production selection,
# calibration, thresholds, abstention, or the production selector, and it never
# reads benchmark labels or expected composition. See
# .omo/plans/common-template-frozen-contrast-v2.md (Todo 1) and AGENTS.md.
#
# Design contract:
#   * The validator is the single source of truth for bound correctness; the
#     struct carries NO physical defaults. Every physical value comes from the
#     TOML; a missing value is an ArgumentError, never a silent default.
#   * Translation/rotation bounds are symmetric about zero (min == -max);
#     blur sigma bounds are non-negative (sigma has no sign).
#   * Every min/max/step is finite; every step is strictly positive.
#   * Every range width is an integer multiple of its step so the search grid
#     is deterministic and endpoint-closed.
#   * Each fine range span never exceeds its coarse span and the fine step
#     never exceeds the coarse step (coarse-to-fine refinement).
#   * Numerical tie tolerances are finite, non-negative, and not both zero.
module WholeRoiDiagnosticConfig

using TOML

export DiagnosticRange, DiagnosticParamBounds, DiagnosticConfig,
       parse_diagnostic_config, load_diagnostic_config, search_values,
       DEFAULT_DIAGNOSTIC_CONFIG_PATH

const DEFAULT_DIAGNOSTIC_CONFIG_PATH =
    normpath(joinpath(@__DIR__, "..", "..", "..", "config", "joint_proxy_whole_roi_diagnostic.toml"))

struct DiagnosticRange
    min::Float64
    max::Float64
    step::Float64
end

struct DiagnosticParamBounds
    coarse::DiagnosticRange
    fine::DiagnosticRange
end

struct DiagnosticConfig
    shift_t_nm::DiagnosticParamBounds
    shift_u_nm::DiagnosticParamBounds
    rotation_deg::DiagnosticParamBounds
    blur_sigma_nm::DiagnosticParamBounds
    tie_rtol::Float64
    tie_atol::Float64
end

# Validation kind per parameter. Translation/rotation are signed and constrained
# symmetric about zero; blur sigma is a magnitude and constrained non-negative.
const _PARAM_KINDS = Dict(
    "shift_t_nm"    => :symmetric,
    "shift_u_nm"    => :symmetric,
    "rotation_deg"  => :symmetric,
    "blur_sigma_nm" => :nonnegative,
)

const _PARAM_ORDER = ("shift_t_nm", "shift_u_nm", "rotation_deg", "blur_sigma_nm")

# Relative tolerance for "is this an integer step count" / span-equality checks.
# Wide enough to absorb Float64 rounding of decimal TOML literals, far too tight
# to mask a genuinely non-divisible or wider range.
const _FP_REL_TOL = 256 * eps(Float64)

function _require_table(v, label)
    v isa AbstractDict ||
        throw(ArgumentError("$label must be a table, got value of type $(typeof(v))"))
    return v
end

function _require_finite_number(v, label)
    v isa Real ||
        throw(ArgumentError("$label must be a number, got value of type $(typeof(v))"))
    fv = Float64(v)
    isfinite(fv) || throw(ArgumentError("$label must be finite, got $fv"))
    return fv
end

function _range_from_dict(d, param::AbstractString, level::AbstractString)
    d isa AbstractDict ||
        throw(ArgumentError("[registration.$level.$param] must be a table, got value of type $(typeof(d))"))
    haskey(d, "min")  || throw(ArgumentError("[registration.$level.$param].min missing"))
    haskey(d, "max")  || throw(ArgumentError("[registration.$level.$param].max missing"))
    haskey(d, "step") || throw(ArgumentError("[registration.$level.$param].step missing"))
    mn = _require_finite_number(d["min"],  "[registration.$level.$param].min")
    mx = _require_finite_number(d["max"],  "[registration.$level.$param].max")
    st = _require_finite_number(d["step"], "[registration.$level.$param].step")
    st > 0 || throw(ArgumentError("[registration.$level.$param].step must be positive, got $st"))
    mn <= mx || throw(ArgumentError("[registration.$level.$param]: min ($mn) must be <= max ($mx)"))
    kind = _PARAM_KINDS[param]
    if kind === :symmetric
        mn == -mx || throw(ArgumentError(
            "[registration.$level.$param] bounds must be symmetric (min == -max); got min=$mn max=$mx"))
    else # :nonnegative
        mn >= 0 || throw(ArgumentError("[registration.$level.$param] bounds must be non-negative; got min=$mn"))
        mx >= 0 || throw(ArgumentError("[registration.$level.$param] bounds must be non-negative; got max=$mx"))
    end
    # Determinism: the search grid must land on both endpoints, so the width
    # must be an integer multiple of the step (a degenerate min == max is allowed
    # and means the degree of freedom is frozen at that single value).
    if mx > mn
        nsteps = (mx - mn) / st
        rel = abs(nsteps - round(nsteps)) / max(1.0, abs(nsteps))
        rel <= _FP_REL_TOL || throw(ArgumentError(
            "[registration.$level.$param]: width ($(mx - mn)) must be an integer multiple of step ($st) for a deterministic grid"))
    end
    return DiagnosticRange(mn, mx, st)
end

function _check_fine_vs_coarse(param::AbstractString, fine::DiagnosticRange, coarse::DiagnosticRange)
    coarse_span = coarse.max - coarse.min
    fine_span   = fine.max   - fine.min
    fine_span <= coarse_span + _FP_REL_TOL * max(1.0, abs(coarse_span)) || throw(ArgumentError(
        "[registration.$param]: fine range span ($fine_span) is wider than coarse span ($coarse_span)"))
    fine.step <= coarse.step || throw(ArgumentError(
        "[registration.$param]: fine step ($(fine.step)) must be <= coarse step ($(coarse.step)) for coarse-to-fine refinement"))
    return nothing
end

"""
    parse_diagnostic_config(cfg::AbstractDict) -> DiagnosticConfig

Validate a parsed TOML dict and return a typed immutable `DiagnosticConfig`.
Throws `ArgumentError` naming the offending field for any missing section,
missing value, non-finite value, non-positive step, asymmetric or negative
range, non-deterministic grid, fine range wider/coarser than coarse, or invalid
numerical tie tolerance. Performs NO defaulting: every physical value must be
present in `cfg`.
"""
function parse_diagnostic_config(cfg)
    cfg isa AbstractDict ||
        throw(ArgumentError("diagnostic config must be a table, got value of type $(typeof(cfg))"))
    haskey(cfg, "registration") || throw(ArgumentError("config missing [registration] section"))
    reg = _require_table(cfg["registration"], "[registration]")
    haskey(reg, "coarse") || throw(ArgumentError("config missing [registration.coarse] section"))
    haskey(reg, "fine")   || throw(ArgumentError("config missing [registration.fine] section"))
    coarse_d = _require_table(reg["coarse"], "[registration.coarse]")
    fine_d   = _require_table(reg["fine"],   "[registration.fine]")

    bounds = Dict{String,DiagnosticParamBounds}()
    for param in _PARAM_ORDER
        haskey(coarse_d, param) ||
            throw(ArgumentError("[registration.coarse] missing parameter '$param'"))
        haskey(fine_d, param) ||
            throw(ArgumentError("[registration.fine] missing parameter '$param'"))
        cr = _range_from_dict(coarse_d[param], param, "coarse")
        fr = _range_from_dict(fine_d[param],   param, "fine")
        _check_fine_vs_coarse(param, fr, cr)
        bounds[param] = DiagnosticParamBounds(cr, fr)
    end

    haskey(cfg, "numerical_tie") || throw(ArgumentError("config missing [numerical_tie] section"))
    tie = _require_table(cfg["numerical_tie"], "[numerical_tie]")
    haskey(tie, "rtol") || throw(ArgumentError("[numerical_tie].rtol missing"))
    haskey(tie, "atol") || throw(ArgumentError("[numerical_tie].atol missing"))
    rtol = _require_finite_number(tie["rtol"], "numerical_tie.rtol")
    atol = _require_finite_number(tie["atol"], "numerical_tie.atol")
    rtol >= 0 || throw(ArgumentError("numerical_tie.rtol must be non-negative, got $rtol"))
    atol >= 0 || throw(ArgumentError("numerical_tie.atol must be non-negative, got $atol"))
    (rtol > 0 || atol > 0) || throw(ArgumentError(
        "numerical_tie: rtol and atol must not both be zero (no scale-aware tie tolerance)"))

    return DiagnosticConfig(
        bounds["shift_t_nm"], bounds["shift_u_nm"],
        bounds["rotation_deg"], bounds["blur_sigma_nm"],
        rtol, atol)
end

"""
    load_diagnostic_config(path::AbstractString) -> DiagnosticConfig

Read and validate a diagnostic registration TOML. Throws `ArgumentError` if the
file is missing or fails validation. A syntactically valid TOML that lacks the
required sections is rejected here rather than appearing to succeed.
"""
function load_diagnostic_config(path::AbstractString)
    isfile(path) || throw(ArgumentError("diagnostic config file not found: $path"))
    return parse_diagnostic_config(TOML.parsefile(path))
end

"""
    search_values(r::DiagnosticRange) -> Vector{Float64}

Deterministic, endpoint-closed search grid for one range. A degenerate range
(`min == max`) returns the single frozen value.
"""
function search_values(r::DiagnosticRange)
    r.min == r.max && return Float64[r.min]
    n = round(Int, (r.max - r.min) / r.step)
    n >= 1 || throw(ArgumentError(
        "range cannot be enumerated (width/step < 1): min=$(r.min) max=$(r.max) step=$(r.step)"))
    vals = Vector{Float64}(undef, n + 1)
    @inbounds for k in 0:n
        vals[k + 1] = k == n ? r.max : r.min + k * r.step
    end
    return vals
end

end # module

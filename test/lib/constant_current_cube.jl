module ConstantCurrentCube

export CubeGrid, FrameRef, CrossingResult, SurfaceResult,
       first_vacuum_crossing, constant_current_surface, isovalue_for_mean_height,
       DEFAULT_ISOVALUE_SCAN_INTERVALS

using Statistics
using LinearAlgebra

# --------------------------------------------------------------------------------------
# Validated grid metadata
# --------------------------------------------------------------------------------------

"`CubeGrid` holds a parsed Gaussian cube as origin, column-major axes, grid
sizes, and a flat value vector matching `prod(n)`. All fields are validated."
struct CubeGrid
    origin::Vector{Float64}
    axes::Matrix{Float64}
    n::NTuple{3,Int}
    values::Vector{Float64}

    function CubeGrid(origin::AbstractVector{<:Real}, axes::AbstractMatrix{<:Real},
                      n::NTuple{N,Integer}, values::AbstractVector{<:Real}) where {N}
        N == 3 ||
            throw(ArgumentError("n must have length 3, got $N"))
        length(origin) == 3 ||
            throw(ArgumentError("origin must have length 3, got $(length(origin))"))
        size(axes) == (3, 3) ||
            throw(ArgumentError("axes must be 3x3, got $(size(axes))"))
        all(>(0), n) ||
            throw(ArgumentError("n must be three positive integers, got $n"))
        expected = prod(n)
        length(values) == expected ||
            throw(ArgumentError("values length must equal prod(n)=$expected, got $(length(values))"))
        all(isfinite, origin) ||
            throw(ArgumentError("origin must be entirely finite"))
        all(isfinite, axes) ||
            throw(ArgumentError("axes must be entirely finite"))
        all(isfinite, values) ||
            throw(ArgumentError("cube values must be entirely finite; reject the cube rather than masking"))
        det(axes) != 0 ||
            throw(ArgumentError("cube axes must be non-singular (non-zero determinant)"))
        new(copy(Vector{Float64}(origin)), copy(Matrix{Float64}(axes)),
            Int.(n), copy(Vector{Float64}(values)))
    end
end

"`FrameRef` is a minimal, non-nullable local frame used for height calibration."
struct FrameRef
    origin::Vector{Float64}
    t_axis::Vector{Float64}
    u_axis::Vector{Float64}

    function FrameRef(origin::AbstractVector{<:Real}, t_axis::AbstractVector{<:Real},
                      u_axis::AbstractVector{<:Real})
        length(origin) == 3 || throw(ArgumentError("origin must have length 3"))
        length(t_axis) == 3 || throw(ArgumentError("t_axis must have length 3"))
        length(u_axis) == 3 || throw(ArgumentError("u_axis must have length 3"))
        all(isfinite, origin) || throw(ArgumentError("origin must be finite"))
        all(isfinite, t_axis) || throw(ArgumentError("t_axis must be finite"))
        all(isfinite, u_axis) || throw(ArgumentError("u_axis must be finite"))
        norm(t_axis) > 0 || throw(ArgumentError("t_axis must be nonzero"))
        # t and u must not be collinear (the Gram-Schmidt residual must survive).
        t_hat = t_axis ./ norm(t_axis)
        residual = u_axis .- dot(u_axis, t_hat) .* t_hat
        norm(residual) > 0 || throw(ArgumentError("u_axis is collinear with t_axis"))
        new(copy(Vector{Float64}(origin)), copy(Vector{Float64}(t_axis)),
            copy(Vector{Float64}(u_axis)))
    end
end

# --------------------------------------------------------------------------------------
# Result types
# --------------------------------------------------------------------------------------

"`CrossingResult` reports the first vacuum-side crossing of a 1D column."
struct CrossingResult
    z::Float64       # crossing height (NaN unless status === :found)
    bracket::Int     # upper (vacuum-side) index of the first bracket, 0 if none
    status::Symbol   # :found | :absent | :ambiguous | :nonfinite | :insufficient_z
    n_brackets::Int  # merged crossing count (diagnostic)
end

"`SurfaceResult` collects a full constant-current surface over the native grid."
struct SurfaceResult
    height::Matrix{Float64}     # crossing z (NaN unless valid)
    valid::Matrix{Bool}          # true iff status === :found
    bracket::Matrix{Int}         # first-bracket upper index per column
    status::Matrix{Symbol}       # per-column status
    counts::Dict{Symbol,Int}     # diagnostic status counts
end

# --------------------------------------------------------------------------------------
# first_vacuum_crossing
# --------------------------------------------------------------------------------------

"""
    first_vacuum_crossing(z_column, v_column, isovalue) -> CrossingResult

Search from vacuum (the last, highest-z element of `z_column`) toward the slab
(first element) and linearly interpolate the first bracketed crossing of
`isovalue`. Adjacent brackets that share an exact-touch grid point (a value
equal to `isovalue`) are merged into a single crossing. Multiple distinct
crossing zones mark the column as `:ambiguous` and the height is set to `NaN`.
"""
function first_vacuum_crossing(z_column::AbstractVector{<:Real},
                               v_column::AbstractVector{<:Real},
                               isovalue::Real)
    nz = length(z_column)
    length(v_column) == nz ||
        throw(ArgumentError("z_column and v_column must have equal length"))
    isfinite(Float64(isovalue)) ||
        throw(ArgumentError("isovalue must be finite"))
    nz >= 2 && return _first_vacuum_crossing_inner(z_column, v_column, Float64(isovalue), nz)
    return CrossingResult(NaN, 0, :insufficient_z, 0)
end

function _first_vacuum_crossing_inner(z_column, v_column, isovalue::Float64, nz::Int)
    for i in 2:nz
        z_column[i] > z_column[i - 1] ||
            throw(ArgumentError("z_column must be strictly ascending"))
    end
    for i in 1:nz
        isfinite(Float64(v_column[i])) ||
            return CrossingResult(NaN, 0, :nonfinite, 0)
    end
    raw = Int[]
    for i in nz:-1:2
        if (Float64(v_column[i]) - isovalue) * (Float64(v_column[i - 1]) - isovalue) <= 0.0
            push!(raw, i)
        end
    end
    brackets = Int[]
    for b in raw
        if !isempty(brackets) && brackets[end] == b + 1 && Float64(v_column[b]) == isovalue
            continue  # merge: shared exact-touch endpoint with the previous bracket
        end
        push!(brackets, b)
    end
    isempty(brackets) && return CrossingResult(NaN, 0, :absent, 0)
    first_i = brackets[1]
    n_found = length(brackets)
    v_hi = Float64(v_column[first_i])
    v_lo = Float64(v_column[first_i - 1])
    z_hi = Float64(z_column[first_i])
    z_lo = Float64(z_column[first_i - 1])
    denom = v_lo - v_hi
    if denom == 0.0
        return CrossingResult(NaN, first_i, :ambiguous, n_found)
    end
    frac = (isovalue - v_hi) / denom
    z_cross = z_hi + frac * (z_lo - z_hi)
    status = n_found == 1 ? :found : :ambiguous
    return CrossingResult(status === :found ? z_cross : NaN, first_i, status, n_found)
end

# --------------------------------------------------------------------------------------
# constant_current_surface
# --------------------------------------------------------------------------------------

"""
    constant_current_surface(cube, isovalue) -> SurfaceResult

For every native (ix, iy) column, find the first vacuum-side crossing of
`isovalue` along the cube z-axis. The cube's third axis must point toward
vacuum (`axes[3,3] > 0`), matching the convention that high z is vacuum and
low z is slab. Returns per-column height, validity, bracket index, status,
and aggregate diagnostic counts.
"""
function constant_current_surface(cube::CubeGrid, isovalue::Real)
    nx, ny, nz = cube.n
    z_step = cube.axes[3, 3]
    z_step > 0 ||
        throw(ArgumentError("cube z-axis must point toward vacuum (axes[3,3] > 0); got $z_step"))
    isfinite(Float64(isovalue)) ||
        throw(ArgumentError("isovalue must be finite"))
    iso = Float64(isovalue)
    height = Matrix{Float64}(undef, nx, ny)
    valid = Matrix{Bool}(undef, nx, ny)
    bracket = zeros(Int, nx, ny)
    status = Matrix{Symbol}(undef, nx, ny)
    counts = Dict{Symbol,Int}(:found => 0, :absent => 0, :ambiguous => 0,
                              :nonfinite => 0, :insufficient_z => 0)
    for ix in 1:nx, iy in 1:ny
        z_base = cube.origin[3] + (ix - 1) * cube.axes[3, 1] + (iy - 1) * cube.axes[3, 2]
        z_col = [z_base + (iz - 1) * z_step for iz in 1:nz]
        offset = ((ix - 1) * ny + (iy - 1)) * nz
        v_col = view(cube.values, offset + 1:offset + nz)
        result = _first_vacuum_crossing_inner(z_col, v_col, iso, nz)
        height[ix, iy] = result.z
        valid[ix, iy] = result.status === :found
        bracket[ix, iy] = result.bracket
        status[ix, iy] = result.status
        counts[result.status] = get(counts, result.status, 0) + 1
    end
    return SurfaceResult(height, valid, bracket, status, counts)
end

# --------------------------------------------------------------------------------------
# Frame helpers
# --------------------------------------------------------------------------------------

function _frame_unit_vectors(frame::FrameRef)
    t_hat = frame.t_axis ./ norm(frame.t_axis)
    u_proj = frame.u_axis .- dot(frame.u_axis, t_hat) .* t_hat
    u_hat = u_proj ./ norm(u_proj)
    n_hat = cross(t_hat, u_hat)
    return t_hat, u_hat, n_hat ./ norm(n_hat)
end

"`heights above the frame origin along the normal for every valid surface column."
function _valid_projected_heights(surf::SurfaceResult, cube::CubeGrid,
                                  frame::FrameRef, n_hat::Vector{Float64})
    nx, ny = cube.n[1], cube.n[2]
    hs = Float64[]
    for ix in 1:nx, iy in 1:ny
        surf.valid[ix, iy] || continue
        x = cube.origin[1] + (ix - 1) * cube.axes[1, 1] + (iy - 1) * cube.axes[1, 2]
        y = cube.origin[2] + (ix - 1) * cube.axes[2, 1] + (iy - 1) * cube.axes[2, 2]
        rx = x - frame.origin[1]
        ry = y - frame.origin[2]
        rz = surf.height[ix, iy] - frame.origin[3]
        push!(hs, n_hat[1] * rx + n_hat[2] * ry + n_hat[3] * rz)
    end
    return hs
end

# --------------------------------------------------------------------------------------
# isovalue_for_mean_height
# --------------------------------------------------------------------------------------

# Numerical discovery policy: inspect 1025 equally spaced isovalues, including
# both declared endpoints. This is a fixed ambiguity-detection resolution, not a
# physical branch preference; calibration rejects rather than choosing if the
# scan exposes more than one target root or a discontinuous target bracket.
const DEFAULT_ISOVALUE_SCAN_INTERVALS = 1024

function _validated_scan_intervals(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("scan_intervals must be a positive integer; got $(repr(value))"))
    value > 0 ||
        throw(ArgumentError("scan_intervals must be a positive integer; got $(repr(value))"))
    return Int(value)
end

function _validated_max_iter(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("max_iter must be a positive integer; got $(repr(value))"))
    value > 0 ||
        throw(ArgumentError("max_iter must be a positive integer; got $(repr(value))"))
    return Int(value)
end

function _validated_height_tol_factor(value)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("height_tol_factor must be a positive finite real; got $(repr(value))"))
    factor = Float64(value)
    isfinite(factor) && factor > 0 ||
        throw(ArgumentError("height_tol_factor must be a positive finite real; got $(repr(value))"))
    return factor
end

"""
    isovalue_for_mean_height(cube, frame, target_height_nm; iso_lo, iso_hi, ...) -> Float64

Find the isovalue whose mean valid crossing height (projected onto the frame
normal above the frame origin) equals `target_height_nm`. A bounded scan with
`scan_intervals + 1` points discovers interior target brackets before bisection;
`scan_intervals` defaults to 1024 and must be a positive integer. The scan never
selects among multiple roots, rejects exact roots at either scan-range endpoint,
requires finite immediate neighbours with unchanged support for an exact
interior root, and rejects a target bracket whose valid-column support changes.
Throws explicitly for ambiguous, unreachable, and entirely invalid responses. No
nA-to-cube conversion is claimed; the isovalue is a discrete ``|psi_n(r)|^2``
setpoint fixed by the 0.50 nm mean-height policy.
"""
function isovalue_for_mean_height(cube::CubeGrid, frame::FrameRef,
                                  target_height_nm::Real;
                                  iso_lo::Union{Nothing,Real}=nothing,
                                  iso_hi::Union{Nothing,Real}=nothing,
                                  max_iter=80, height_tol_factor=1e-6,
                                  scan_intervals=DEFAULT_ISOVALUE_SCAN_INTERVALS)
    n_max_iter = _validated_max_iter(max_iter)
    tol_factor = _validated_height_tol_factor(height_tol_factor)
    target = Float64(target_height_nm)
    isfinite(target) || throw(ArgumentError("target_height_nm must be finite"))
    n_scan_intervals = _validated_scan_intervals(scan_intervals)
    _, _, n_hat = _frame_unit_vectors(frame)
    finite_vals = cube.values
    v_min, v_max = extrema(finite_vals)
    v_range = v_max - v_min
    v_range > 0 || error("cube values are constant; cannot calibrate isovalue")
    lo = iso_lo === nothing ? v_min + 1e-3 * v_range : Float64(iso_lo)
    hi = iso_hi === nothing ? v_max - 1e-3 * v_range : Float64(iso_hi)
    lo < hi || error("invalid isovalue search range: [$lo, $hi]")

    z_step = cube.axes[3, 3]
    height_tol = abs(z_step) * tol_factor

    function mean_height(iso::Float64)
        surf = constant_current_surface(cube, iso)
        hs = _valid_projected_heights(surf, cube, frame, n_hat)
        return isempty(hs) ? NaN : mean(hs)
    end

    function response_with_support(iso::Float64)
        surf = constant_current_surface(cube, iso)
        hs = _valid_projected_heights(surf, cube, frame, n_hat)
        return isempty(hs) ? (NaN, BitVector(vec(surf.valid))) :
                             (mean(hs), BitVector(vec(surf.valid)))
    end

    scan_iso = collect(range(lo, hi; length=n_scan_intervals + 1))
    scan_h = Float64[mean_height(iso) for iso in scan_iso]
    finite_h = filter(isfinite, scan_h)
    isempty(finite_h) && error("no valid crossings for any of $(n_scan_intervals + 1) "
                               * "scanned isovalues in [$lo, $hi]")

    exact = findall(h -> isfinite(h) && abs(h - target) < height_tol, scan_h)
    brackets = Tuple{Int,Int}[]
    for i in 1:n_scan_intervals
        h1, h2 = scan_h[i], scan_h[i + 1]
        isfinite(h1) && isfinite(h2) || continue
        abs(h1 - target) < height_tol && continue
        abs(h2 - target) < height_tol && continue
        (h1 - target) * (h2 - target) < 0 && push!(brackets, (i, i + 1))
    end
    n_candidates = length(exact) + length(brackets)
    n_candidates > 1 &&
        error("ambiguous nonmonotone mean-height response: multiple target roots "
              * "($n_candidates) found in isovalue range [$lo, $hi]")
    if n_candidates == 0
        error("target height $target nm is outside the reachable mean-height range "
              * "[$(minimum(finite_h)), $(maximum(finite_h))] nm for isovalue in [$lo, $hi]")
    end

    if length(exact) == 1
        i = only(exact)
        if i == 1 || i == length(scan_h)
            error("ambiguous endpoint root at isovalue=$(scan_iso[i]): continuity "
                  * "and a crossing within the declared scan range cannot be established")
        end
        if !isfinite(scan_h[i - 1]) || !isfinite(scan_h[i + 1])
            error("ambiguous isolated-response continuity at exact interior root "
                  * "isovalue=$(scan_iso[i]): both immediate adjacent scan responses "
                  * "must be finite")
        end
        # A sampled touch with same-sided neighbours is a turning point, not a
        # unique monotone target crossing. A multi-sample tolerance plateau is
        # already rejected above as multiple candidates.
        if (scan_h[i - 1] - target) * (scan_h[i + 1] - target) > 0
            error("ambiguous nonmonotone mean-height response at isovalue=$(scan_iso[i])")
        end
        _, support_lo = response_with_support(scan_iso[i - 1])
        _, support = response_with_support(scan_iso[i])
        _, support_hi = response_with_support(scan_iso[i + 1])
        support_lo == support == support_hi ||
            error("ambiguous discontinuous mean-height response: valid-column support changes "
                  * "at target isovalue=$(scan_iso[i])")
        return scan_iso[i]
    end

    i_lo, i_hi = only(brackets)
    lo, hi = scan_iso[i_lo], scan_iso[i_hi]
    h_lo, support = response_with_support(lo)
    h_hi, support_hi = response_with_support(hi)
    support == support_hi ||
        error("ambiguous discontinuous mean-height response: valid-column support changes "
              * "across target bracket [$lo, $hi]")
    for _ in 1:n_max_iter
        mid = 0.5 * (lo + hi)
        h_mid, support_mid = response_with_support(mid)
        isfinite(h_mid) || error("ambiguous discontinuous mean-height response: "
                                 * "no valid crossings at isovalue=$mid during bisection")
        support_mid == support ||
            error("ambiguous discontinuous mean-height response: valid-column support changes "
                  * "at isovalue=$mid during bisection")
        abs(h_mid - target) < height_tol && return mid
        if (h_mid > target) == (h_lo > target)
            lo = mid
            h_lo = h_mid
        else
            hi = mid
            h_hi = h_mid
        end
    end
    return 0.5 * (lo + hi)
end

end # module

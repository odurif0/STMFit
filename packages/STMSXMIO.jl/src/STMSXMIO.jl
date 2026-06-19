module STMSXMIO

# Shared SXM (Nanonis) I/O layer for the STMFit monorepo.
#
# Owns the SXMImage/SXMChannel types, the .sxm reader, and the low-level
# preprocessing / ROI helpers that are identical between the GaussianFit2D and
# STMMolecularFit engines. Both engines `using STMSXMIO` instead of each
# defining their own copy.
#
# Two intentionally distinct row-flattening conventions are provided, because
# the two engines had silently diverged on this point (see Research Journal):
#   _row_median_flatten_global : preserve global median level (GaussianFit2D)
#   _row_median_flatten_zero   : zero each row's median (STMMolecularFit)
# Each engine imports the variant that matches its prior behaviour.

using Statistics

export SXMImage, SXMChannel, read_sxm, channel_names, get_channel

const EPS = 1e-12

# =============================================================================
# Types
# =============================================================================

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

# =============================================================================
# .sxm reader
# =============================================================================

function _parse_header(header_text::String)
    header = Dict{String,String}()
    current = nothing
    buf = IOBuffer()
    for line in split(header_text, '\n')
        if startswith(line, ":") && endswith(strip(line), ":")
            if current !== nothing
                header[current] = strip(String(take!(buf)))
            end
            current = strip(line, [':', ' ', '\r', '\n', '\t'])
        elseif current !== nothing
            println(buf, line)
        end
    end
    if current !== nothing
        header[current] = strip(String(take!(buf)))
    end
    return header
end

function _parse_pair(value::String, T=Float64)
    vals = split(strip(value))
    length(vals) >= 2 || error("Expected two values, got: $value")
    return parse(T, vals[1]), parse(T, vals[2])
end

function _parse_data_info(value::String)
    rows = split(value, '\n')
    infos = NamedTuple[]
    for row in rows
        s = strip(row)
        isempty(s) && continue
        startswith(s, "Channel") && continue
        parts = split(s)
        length(parts) >= 6 || continue
        push!(infos, (name=parts[2], unit=parts[3], direction=parts[4]))
    end
    return infos
end

function _read_be_float32(bytes::Vector{UInt8}, offset::Int, n::Int)
    io = IOBuffer(bytes[offset:(offset + 4n - 1)])
    vals = Vector{Float64}(undef, n)
    for i in 1:n
        u = read(io, UInt32)
        ENDIAN_BOM == 0x04030201 && (u = bswap(u))
        vals[i] = Float64(reinterpret(Float32, u))
    end
    return vals
end

function read_sxm(filepath::String)
    bytes = read(filepath)
    marker = Vector{UInt8}(":SCANIT_END:")
    idx = findfirst(marker, bytes)
    idx === nothing && error("SXM marker :SCANIT_END: not found in $filepath")
    header = _parse_header(String(bytes[1:(first(idx)-1)]))
    # Guard the mandatory Nanonis header keys: a truncated/malformed header used
    # to throw a cryptic KeyError instead of a readable message.
    for key in ("SCAN_PIXELS", "SCAN_RANGE", "DATA_INFO")
        haskey(header, key) || error("Malformed SXM header in $filepath: missing key '$key'")
    end

    width, height = _parse_pair(header["SCAN_PIXELS"], Int)
    rx, ry = _parse_pair(header["SCAN_RANGE"], Float64)
    ox, oy = _parse_pair(get(header, "SCAN_OFFSET", "0 0"), Float64)
    range_nm = (rx * 1e9, ry * 1e9)
    offset_nm = (ox * 1e9, oy * 1e9)

    expanded = Tuple{String,String,String}[]
    for info in _parse_data_info(header["DATA_INFO"])
        dirs = lowercase(info.direction) == "both" ? ["fwd", "bwd"] : [lowercase(info.direction)]
        for dir in dirs
            push!(expanded, (info.name, info.unit, dir))
        end
    end

    nvals = width * height * length(expanded)
    data_offset = length(bytes) - 4nvals + 1
    data_offset > 0 || error("Invalid SXM data size")
    vals = _read_be_float32(bytes, data_offset, nvals)

    channels = SXMChannel[]
    k = 1
    for (name, unit, dir) in expanded
        raw = vals[k:(k + width * height - 1)]
        k += width * height
        mat = permutedims(reshape(raw, width, height))
        # Nanonis stores backward scans in acquisition order. Flip x so fwd/bwd
        # share the same spatial coordinates. For the example file, Z fwd vs
        # bwd correlation changes from ≈ -0.35 to ≈ 0.96 after this flip.
        lowercase(dir) == "bwd" && (mat = reverse(mat; dims=2))
        push!(channels, SXMChannel(name, unit, dir, mat))
    end
    return SXMImage(filepath, header, width, height, range_nm, offset_nm, channels)
end

# =============================================================================
# Channel access
# =============================================================================

channel_names(img::SXMImage) = ["$(c.name) ($(c.direction), $(c.unit))" for c in img.channels]

function get_channel(img::SXMImage, name::String; direction::String="fwd")
    lname, ldir = lowercase(name), lowercase(direction)
    for c in img.channels
        lowercase(c.name) == lname && lowercase(c.direction) == ldir && return c
    end
    # Fallback: return any channel matching the name (ignore direction). This
    # preserves the historical STMMolecularFit behaviour, where a missing
    # direction is non-fatal as long as the channel exists.
    for c in img.channels
        lowercase(c.name) == lname && return c
    end
    error("Channel '$name' direction '$direction' not found. Available: $(join(channel_names(img), ", "))")
end

# =============================================================================
# Shared preprocessing helpers (used by both engines' preprocess_channel)
# =============================================================================

_coordinate_vectors(img::SXMImage; stride::Int=1) = (
    collect(range(0, img.range_nm[1], length=img.width)[1:stride:end]),
    collect(range(0, img.range_nm[2], length=img.height)[1:stride:end]),
)

function _value_scale(unit::String)
    u = lowercase(strip(unit))
    u == "m" && return 1e9, "nm"
    u == "a" && return 1e12, "pA"
    return 1.0, isempty(unit) ? "a.u." : unit
end

function _plane_fit(xs, ys, z::Matrix{Float64})
    xflat = repeat(xs, inner=length(ys))
    yflat = repeat(ys, outer=length(xs))
    zflat = vec(z) # z[y,x] aligned with xflat/yflat: y varies fastest
    coeff = hcat(ones(length(xflat)), xflat, yflat) \ zflat
    return [coeff[1] + coeff[2] * x + coeff[3] * y for y in ys, x in xs]
end

# Row flattening, variant A: preserve the global median level. This is the
# GaussianFit2D convention. Rows are equalized but the image-wide level is kept.
function _row_median_flatten_global(z::Matrix{Float64})
    out = copy(z)
    global_med = median(vec(out))
    for iy in axes(out, 1)
        out[iy, :] .-= median(view(out, iy, :)) - global_med
    end
    return out
end

# Row flattening, variant B: zero each row's own median. This is the
# STMMolecularFit convention. The global level shifts down to ≈0.
function _row_median_flatten_zero(z::Matrix{Float64})
    out = copy(z)
    for iy in axes(out, 1)
        out[iy, :] .-= median(out[iy, :])
    end
    return out
end

function _box_smooth(z::Matrix{Float64}, radius::Int)
    radius <= 0 && return copy(z)
    out = similar(z)
    ny, nx = size(z)
    for iy in 1:ny, ix in 1:nx
        ylo, yhi = max(1, iy-radius), min(ny, iy+radius)
        xlo, xhi = max(1, ix-radius), min(nx, ix+radius)
        out[iy, ix] = mean(@view z[ylo:yhi, xlo:xhi])
    end
    return out
end

# Otsu's automatic threshold: maximizes inter-class variance.
function _otsu_threshold(signal::AbstractMatrix{Float64})
    v = vec(signal)
    v = v[isfinite.(v) .& (v .> 0)]
    isempty(v) && return 0.0
    nbins = 128
    lo, hi = minimum(v), maximum(v)
    hi <= lo && return 0.0
    bin_edges = range(lo, hi, length=nbins+1)
    counts = zeros(Int, nbins)
    for x in v
        idx = min(nbins, max(1, Int(floor((x-lo)/(hi-lo)*nbins))+1))
        counts[idx] += 1
    end
    total = sum(counts)
    total == 0 && return 0.0
    sumB, wB, max_var, best = 0.0, 0.0, 0.0, 0
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2
    sumT = sum(bin_centers .* counts)
    for i in 1:nbins
        wB += counts[i]
        wB == 0 && continue
        wF = total - wB
        wF <= 0 && break
        sumB += bin_centers[i] * counts[i]
        mB = sumB / wB
        mF = (sumT - sumB) / wF
        var_between = wB * wF * (mB - mF)^2 / (total^2)
        if var_between > max_var
            max_var = var_between
            best = i
        end
    end
    return best > 0 ? bin_centers[best] : 0.0
end

function _largest_component(mask::BitMatrix)
    ny, nx = size(mask)
    seen = falses(ny, nx)
    best = Tuple{Int,Int}[]
    for iy in 1:ny, ix in 1:nx
        (!mask[iy, ix] || seen[iy, ix]) && continue
        comp = Tuple{Int,Int}[]
        stack = [(iy, ix)]
        seen[iy, ix] = true
        while !isempty(stack)
            y, x = pop!(stack)
            push!(comp, (y, x))
            for yy in max(1, y-1):min(ny, y+1), xx in max(1, x-1):min(nx, x+1)
                if mask[yy, xx] && !seen[yy, xx]
                    seen[yy, xx] = true
                    push!(stack, (yy, xx))
                end
            end
        end
        length(comp) > length(best) && (best = comp)
    end
    out = falses(ny, nx)
    for idx in best
        out[idx...] = true
    end
    return out
end

# Separable dilation by a square structuring element (dilate along rows then
# columns): equivalent to the (2r+1)^2 square but O(n*r) instead of O(n*r^2).
function _dilate_mask(mask::BitMatrix, radius::Int)
    radius <= 0 && return copy(mask)
    ny, nx = size(mask)
    row = falses(ny, nx)
    @inbounds for iy in 1:ny, ix in 1:nx
        if mask[iy, ix]
            row[iy, max(1, ix-radius):min(nx, ix+radius)] .= true
        end
    end
    out = falses(ny, nx)
    @inbounds for iy in 1:ny, ix in 1:nx
        if row[iy, ix]
            out[max(1, iy-radius):min(ny, iy+radius), ix] .= true
        end
    end
    return out
end

end # module

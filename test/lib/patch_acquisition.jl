module PatchAcquisition

# Opt-in acquisition transforms only. The saved Gaussian model stays in its
# original coordinate frame; shift the observation BEFORE model subtraction.
using STMSXMIO: _box_smooth
export read_shifts, observed_shift, require_direction

function read_shifts(path, files)
    path === nothing && return nothing
    lines = readlines(path)
    !isempty(lines) && first(lines) == "file\tbwd_sample_dx_px" || error("Invalid acquisition shift header")
    out = Dict{String,Int}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, '\t'; keepempty=true)
        length(parts) == 2 || error("Invalid acquisition shift row")
        file = String(parts[1])
        file == basename(file) && !isempty(file) || error("Invalid acquisition filename")
        haskey(out,file) && error("Duplicate acquisition filename")
        dx = tryparse(Int,parts[2]); dx === nothing && error("Shift must be an integer")
        out[file] = dx
    end
    Set(keys(out)) == Set(files) || error("Acquisition shifts and feature files differ")
    return out
end

function require_direction(ch, direction)
    lowercase(ch.name) == "z" && lowercase(ch.direction) == direction ||
        error("Actual Z $direction channel required; no direction fallback")
    return nothing
end

"""Return Z(y,x+dx), restoring actual missing pixels before native smoothing.
No wrap, edge padding, gain/offset fit, or new interpolation. A zero transform on
a fully observed image reproduces native arrays exactly. Flattening still comes
from native preprocessing (its earlier imputation influence is not undone).
"""
function observed_shift(z, ch, settings, dx::Int)
    observed = isfinite.(ch.data[1:settings.stride:end,1:settings.stride:end])
    size(z) == size(observed) || throw(DimensionMismatch("Raw mask differs from processed grid"))
    ny,nx = size(z)
    abs(dx) < nx || error("Acquisition shift exceeds image width")
    out = fill(NaN,ny,nx)
    for x in 1:nx
        source = x+dx
        1 <= source <= nx || continue
        for y in 1:ny
            observed[y,source] && (out[y,x] = z[y,source])
        end
    end
    return out, _box_smooth(out,settings.smooth_radius_px)
end

end

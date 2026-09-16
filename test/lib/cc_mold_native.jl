module CCMoldNative

using LinearAlgebra
using Printf
using Statistics
using TOML

export MoldSettings, Cube, Frame, load_settings, read_cube, read_frame,
       sample_volume, surface_grid, iso_for_mean_height_legacy,
       build_templates, build_molds

# The unit conversion and array conventions are those of cc_mold_builder.py.
const BOHR_NM = 0.05291772109

struct MoldSettings
    half_nm::Float64
    step_nm::Float64
    target_height_nm::Float64
    z_min_nm::Float64
    z_max_nm::Float64
    z_step_nm::Float64
    isovalue_min_log10::Float64
    isovalue_max_fraction::Float64
    isovalue_count::Int
    calibration::String
end

function _validate(s::MoldSettings)
    all(isfinite, (s.half_nm, s.step_nm, s.target_height_nm, s.z_min_nm,
                   s.z_max_nm, s.z_step_nm, s.isovalue_min_log10,
                   s.isovalue_max_fraction)) || throw(ArgumentError("nonfinite mold setting"))
    s.half_nm > 0 && s.step_nm > 0 || throw(ArgumentError("mold grid must have positive half-width and step"))
    q = 2 * s.half_nm / s.step_nm
    isfinite(q) && q >= 1 && q < typemax(Int) || throw(ArgumentError("invalid mold grid size"))
    abs(q - round(q)) <= 4 * eps(q) || throw(ArgumentError("mold step must divide the full grid width"))
    s.z_min_nm < s.z_max_nm && s.z_step_nm > 0 || throw(ArgumentError("invalid mold z range or step"))
    nz = (s.z_max_nm - s.z_min_nm) / s.z_step_nm
    isfinite(nz) && 1 < nz < typemax(Int) || throw(ArgumentError("mold z grid needs at least two samples"))
    s.isovalue_count >= 2 || throw(ArgumentError("mold isovalue scan needs at least two values"))
    0 < s.isovalue_max_fraction <= 1 || throw(ArgumentError("mold isovalue fraction must be in (0, 1]"))
    isfinite(10.0^s.isovalue_min_log10) && 10.0^s.isovalue_min_log10 > 0 ||
        throw(ArgumentError("invalid mold minimum isovalue"))
    s.calibration == "first_below_target" || throw(ArgumentError("only first_below_target calibration is supported"))
    return s
end

"""Read all mold settings from [model]. No scientific defaults are supplied."""
function load_settings(path::AbstractString)
    config = TOML.parsefile(path)
    model = get(config, "model", nothing)
    model isa AbstractDict || throw(ArgumentError("config requires [model]"))
    function number(key)
        value = get(model, key, nothing)
        value isa Real && !(value isa Bool) || throw(ArgumentError("missing or nonnumeric $key"))
        return Float64(value)
    end
    count = get(model, "mold_isovalue_count", nothing)
    count isa Integer && !(count isa Bool) || throw(ArgumentError("mold_isovalue_count must be an integer"))
    calibration = get(model, "mold_calibration", nothing)
    calibration isa AbstractString || throw(ArgumentError("missing mold_calibration"))
    return _validate(MoldSettings(
        number("mold_half_nm"), number("mold_step_nm"), number("mold_target_height_nm"),
        number("mold_z_min_nm"), number("mold_z_max_nm"), number("mold_z_step_nm"),
        number("mold_isovalue_min_log10"), number("mold_isovalue_max_fraction"),
        Int(count), String(calibration)))
end

struct Cube
    origin_nm::Vector{Float64}
    axes_nm::Matrix{Float64}  # columns are the three voxel step vectors
    dims::NTuple{3,Int}
    values::Vector{Float64}  # keep input token order; do not reshape a Gaussian cube
end

struct Frame
    origin_nm::Vector{Float64}
    t_axis::Vector{Float64}
    u_axis::Vector{Float64}
end

function _float(token::AbstractString)
    x = tryparse(Float64, replace(token, 'D' => 'E', 'd' => 'e'))
    x !== nothing && isfinite(x) || throw(ArgumentError("invalid finite number: $token"))
    return x
end

function _record(io::IO, description::AbstractString)
    eof(io) && throw(ArgumentError("truncated cube: missing $description"))
    return split(readline(io))
end

"""
Read a single scalar cube in bohr. Signed grid sizes are taken in absolute
value and still interpreted as bohr, as in the Python reference. Values remain
in token order. Interpolation deliberately uses i + n1*j + n1*n2*k (zero-based),
NOT Gaussian-standard last-axis-fast ordering. Malformed/extra data fail.
"""
function read_cube(path::AbstractString)
    return open(path, "r") do io
        _record(io, "comment 1"); _record(io, "comment 2")
        header = _record(io, "atom count and origin")
        length(header) in (4, 5) || throw(ArgumentError("invalid cube origin record"))
        nat = abs(parse(Int, header[1]))
        length(header) == 5 && parse(Int, header[5]) != 1 &&
            throw(ArgumentError("only single scalar cubes are supported"))
        origin = _float.(header[2:4]) .* BOHR_NM
        dims = Int[]
        axes = Matrix{Float64}(undef, 3, 3)
        for axis in 1:3
            row = _record(io, "axis $axis")
            length(row) == 4 || throw(ArgumentError("invalid cube axis record"))
            n = abs(parse(Int, row[1]))
            n > 0 || throw(ArgumentError("cube axis size must be positive"))
            push!(dims, n)
            axes[:, axis] = _float.(row[2:4]) .* BOHR_NM
        end
        inverse = try
            inv(axes)
        catch e
            e isa SingularException || rethrow()
            throw(ArgumentError("cube axes are singular"))
        end
        all(isfinite, inverse) || throw(ArgumentError("invalid cube axes"))
        for atom in 1:nat
            row = _record(io, "atom $atom")
            length(row) == 5 || throw(ArgumentError("invalid cube atom record"))
            parse(Int, row[1]); _float.(row[2:5])
        end
        nvalues = Base.checked_mul(Base.checked_mul(dims[1], dims[2]), dims[3])
        values = Vector{Float64}(undef, nvalues)
        n = 0
        for line in eachline(io), token in split(line)
            n += 1
            n <= nvalues || throw(ArgumentError("cube contains extra data; expected $nvalues scalar values"))
            values[n] = _float(token)
        end
        n == nvalues || throw(ArgumentError("truncated cube: expected $nvalues values, got $n"))
        return Cube(origin, axes, Tuple(dims), values)
    end
end

function _normal(frame::Frame)
    all(v -> length(v) == 3 && all(isfinite, v),
        (frame.origin_nm, frame.t_axis, frame.u_axis)) || throw(ArgumentError("invalid frame vectors"))
    nh = cross(frame.t_axis, frame.u_axis)
    scale = norm(nh)
    isfinite(scale) && scale > 0 || throw(ArgumentError("frame t/u axes must be nonzero and nonparallel"))
    # Do not normalize or orthogonalize t/u: the legacy builder does not.
    return nh ./ scale
end

"""Read origin_nm/t_axis/u_axis TSV records; other frame metadata is ignored."""
function read_frame(path::AbstractString)
    fields = Dict{String,Vector{Float64}}()
    for line in eachline(path)
        parts = split(strip(line), '\t'; keepempty=true)
        isempty(parts) && continue
        key = parts[1]
        key in ("origin_nm", "t_axis", "u_axis") || continue
        length(parts) == 2 || throw(ArgumentError("invalid frame record: $key"))
        haskey(fields, key) && throw(ArgumentError("duplicate frame record: $key"))
        values = split(parts[2], ',')
        length(values) == 3 || throw(ArgumentError("frame record must have three coordinates: $key"))
        fields[key] = _float.(values)
    end
    all(k -> haskey(fields, k), ("origin_nm", "t_axis", "u_axis")) ||
        throw(ArgumentError("frame requires origin_nm, t_axis and u_axis"))
    frame = Frame(fields["origin_nm"], fields["t_axis"], fields["u_axis"])
    _normal(frame)
    return frame
end

"""
Trilinear sampling of N×3 points. Partial boundary weights are NOT renormalized;
zero total supported weight gives NaN, matching the reference.
"""
function sample_volume(cube::Cube, points::AbstractMatrix{<:Real})
    size(points, 2) == 3 && all(isfinite, points) || throw(ArgumentError("points must be a finite N×3 matrix"))
    coords = (Matrix{Float64}(points) .- transpose(cube.origin_nm)) * transpose(inv(cube.axes_nm))
    n1, n2, n3 = cube.dims
    out = Vector{Float64}(undef, size(points, 1))
    for p in axes(coords, 1)
        x, y, z = coords[p, 1], coords[p, 2], coords[p, 3]
        i, j, k = floor(Int, x), floor(Int, y), floor(Int, z)
        fx, fy, fz = x - i, y - j, z - k
        total = 0.0
        wsum = 0.0
        for di in 0:1, dj in 0:1, dk in 0:1
            ii, jj, kk = i + di, j + dj, k + dk
            if 0 <= ii < n1 && 0 <= jj < n2 && 0 <= kk < n3
                w = (di == 0 ? 1 - fx : fx) * (dj == 0 ? 1 - fy : fy) * (dk == 0 ? 1 - fz : fz)
                idx = 1 + ii + n1 * jj + n1 * n2 * kk
                total += w * cube.values[idx]
                wsum += w
            end
        end
        out[p] = wsum == 0 ? NaN : total
    end
    return out
end

# NumPy arange uses the rounded (start + step) - start as its actual increment.
# Julia ranges use extra precision, which can change a discrete height decision.
function _numpy_grid(start::Float64, step::Float64, count::Int)
    delta = (start + step) - start
    isfinite(delta) && delta > 0 || throw(ArgumentError("grid step is not representable"))
    return [start + i * delta for i in 0:(count - 1)]
end

"""Sample columns in legacy meshgrid order: u outer, t inner; z ascending."""
function surface_grid(cube::Cube, frame::Frame, settings::MoldSettings)
    s = _validate(settings)
    side = round(Int, 2 * s.half_nm / s.step_nm) + 1
    xy = _numpy_grid(-s.half_nm, s.step_nm, side)
    # z_max is excluded (np.arange); the symmetric t/u grid includes both ends.
    zs = _numpy_grid(s.z_min_nm, s.z_step_nm,
                     ceil(Int, (s.z_max_nm - s.z_min_nm) / s.z_step_nm))
    nh = _normal(frame)
    points = Matrix{Float64}(undef, side * side * length(zs), 3)
    p = 0
    for u in xy, t in xy
        base = frame.origin_nm .+ t .* frame.t_axis .+ u .* frame.u_axis
        for z in zs
            p += 1
            for axis in 1:3
                points[p, axis] = base[axis] + z * nh[axis]
            end
        end
    end
    sampled = sample_volume(cube, points)
    vc = permutedims(reshape(sampled, length(zs), side * side))
    return zs, vc
end

"""
Scan the configured log-spaced isovalues, in order. Return the FIRST whose mean
highest occupied sample (strict LDOS > iso) is below target. No crossing
interpolation, support minimum, continuity check, or nearest-root replacement.
"""
function iso_for_mean_height_legacy(vc::AbstractMatrix{<:Real}, zs::AbstractVector{<:Real}, settings::MoldSettings)
    s = _validate(settings)
    size(vc, 2) == length(zs) && !isempty(zs) && size(vc, 1) > 0 ||
        throw(ArgumentError("incompatible or empty LDOS/height grid"))
    all(isfinite, zs) && all(>(0), diff(zs)) || throw(ArgumentError("z samples must be finite and strictly increasing"))
    any(isinf, vc) && throw(ArgumentError("infinite sampled LDOS"))
    finite_values = filter(isfinite, vec(vc))
    isempty(finite_values) && throw(ArgumentError("no finite LDOS samples"))
    maximum_ldos = maximum(finite_values)
    maximum_ldos > 0 || throw(ArgumentError("sampled LDOS has no positive maximum"))
    upper = log10(maximum_ldos * s.isovalue_max_fraction)
    isfinite(upper) || throw(ArgumentError("invalid upper isovalue"))
    step = (upper - s.isovalue_min_log10) / (s.isovalue_count - 1)
    for i in 0:(s.isovalue_count - 1)
        exponent = i == s.isovalue_count - 1 ? upper : s.isovalue_min_log10 + i * step
        iso = 10.0^exponent
        heights = fill(NaN, size(vc, 1))
        for row in axes(vc, 1)
            for k in reverse(eachindex(zs))
                if vc[row, k] > iso
                    heights[row] = zs[k]
                    break
                end
            end
        end
        good = filter(isfinite, heights)
        isempty(good) && continue
        mh = mean(good)
        if mh < s.target_height_nm
            return (iso=iso, mean_height=mh, nvalid=length(good), heights=heights)
        end
    end
    throw(ArgumentError("no scanned isovalue gives mean height below $(s.target_height_nm) nm"))
end

function _zscore(heights::AbstractVector{<:Real})
    good = filter(isfinite, heights)
    isempty(good) && throw(ArgumentError("height map has no finite values"))
    all(==(first(good)), good) && throw(ArgumentError("height map has zero population standard deviation"))
    mu = mean(good)
    sd = std(good; mean=mu, corrected=false)
    isfinite(sd) && sd > 0 || throw(ArgumentError("height map has invalid population standard deviation"))
    normalized = [isfinite(x) ? (x - mu) / sd : 0.0 for x in heights]
    all(isfinite, normalized) || throw(ArgumentError("nonfinite normalized mold"))
    return normalized
end

"""Build one type in memory, with four rows ordered (parity,mirror) = 00,01,10,11."""
function build_templates(cube::Cube, frame::Frame, settings::MoldSettings)
    zs, vc = surface_grid(cube, frame, settings)
    selection = iso_for_mean_height_legacy(vc, zs, settings)
    normalized = _zscore(selection.heights)
    side = round(Int, 2 * settings.half_nm / settings.step_nm) + 1
    variants = Matrix{Float64}(undef, 4, length(normalized))
    for parity in 0:1, mirror in 0:1, row in 1:side, col in 1:side
        src_row = parity == 0 ? row : side + 1 - row
        src_col = mirror == 0 ? col : side + 1 - col
        variants[1 + 2 * parity + mirror, (row - 1) * side + col] =
            normalized[(src_row - 1) * side + src_col]
    end
    return merge(selection, (normalized=normalized, variants=variants))
end

build_templates(cube_path::AbstractString, frame_path::AbstractString, settings::MoldSettings) =
    build_templates(read_cube(cube_path), read_frame(frame_path), settings)

function _check_new_output(path::AbstractString)
    isempty(path) && throw(ArgumentError("output path is empty"))
    (ispath(path) || islink(path)) && throw(ArgumentError("refusing to overwrite output: $path"))
end

"""
Build the pair and write a scorer-compatible TSV with %.7g pixels. All science
finishes before output creation. Exclusive creation refuses files, symlinks,
directories, and a destination created concurrently. Returns (type0, type1).
"""
function build_molds(; cube0::AbstractString, cube1::AbstractString,
                     frame0::AbstractString, frame1::AbstractString,
                     config::AbstractString, out::AbstractString)
    _check_new_output(out)
    settings = load_settings(config)
    type0 = build_templates(cube0, frame0, settings)
    type1 = build_templates(cube1, frame1, settings)
    buffer = IOBuffer()
    npixels = size(type0.variants, 2)
    println(buffer, join(vcat(["name", "type", "parity", "mirror"],
                             [@sprintf("p%03d", i) for i in 1:npixels]), '\t'))
    for (typ, label, result) in ((0, "GlcN", type0), (1, "GlcNAc", type1))
        for parity in 0:1, mirror in 0:1
            print(buffer, "$(label)_p$(parity)_m$(mirror)\t$typ\t$parity\t$mirror")
            for value in @view result.variants[1 + 2 * parity + mirror, :]
                @printf(buffer, "\t%.7g", value)
            end
            println(buffer)
        end
    end
    payload = take!(buffer)
    mkpath(dirname(abspath(out)))
    flags = Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL
    io = Base.Filesystem.open(out, flags, 0o644)
    try
        write(io, payload)
    finally
        close(io)
    end
    return (type0=type0, type1=type1)
end

end # module

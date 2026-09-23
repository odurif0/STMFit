module PatchFrames

using LinearAlgebra, Statistics, TOML

export load_frame_settings, local_frame_rows, write_frames, read_frames, patch_axis

const GEOMETRY = ["x_nm", "y_nm", "axis_x", "axis_y"]
const HEADER = vcat(["file", "lobe"], GEOMETRY,
    ["frame_axis_x", "frame_axis_y", "tangent_slope", "centerline_rms_nm"])

function load_frame_settings(path)
    cfg = TOML.parsefile(path)
    Set(keys(cfg)) == Set(["model", "selection", "preprocessing"]) || error("Unexpected frame settings")
    for (section, field) in (("model", "centerline_degree"), ("selection", "centerline_fit"),
                             ("preprocessing", "patch_orientation"))
        cfg[section] isa AbstractDict && Set(keys(cfg[section])) == Set([field]) ||
            error("Unexpected frame setting in $section")
    end
    degree = cfg["model"]["centerline_degree"]
    degree isa Integer && !(degree isa Bool) && degree in (1, 2) || error("Degree must be 1 or 2")
    cfg["selection"]["centerline_fit"] == "ordinary_least_squares" || error("Unknown centerline fit")
    cfg["preprocessing"]["patch_orientation"] == "local_tangent" || error("Unknown patch orientation")
    return (degree=Int(degree),)
end

number(row, key) = begin
    x = parse(Float64, row[key])
    isfinite(x) || error("Nonfinite frame geometry: $key")
    x
end
rowkey(row) = (basename(row["file"]), parse(Int, row["lobe"]))

function base_index(rows)
    index = Dict{Tuple{String,Int},Dict{String,String}}()
    for row in rows
        key = rowkey(row)
        !isempty(first(key)) && last(key) > 0 || error("Invalid frame key")
        haskey(index, key) && error("Duplicate frame key: $key")
        all(isfinite(number(row, k)) for k in GEOMETRY) || error("Invalid geometry")
        hypot(number(row, "axis_x"), number(row, "axis_y")) > 0 || error("Zero global axis")
        index[key] = row
    end
    isempty(index) && error("Empty frame geometry")
    return index
end

function local_frame_rows(rows, settings)
    index = base_index(rows)
    output = Dict{String,String}[]
    for file in sort(unique(first.(collect(keys(index)))))
        keys_file = sort([k for k in keys(index) if first(k) == file])
        last.(keys_file) == collect(1:length(keys_file)) || error("Nonconsecutive lobes: $file")
        chain = [index[k] for k in keys_file]
        ax, ay = number(chain[1], "axis_x"), number(chain[1], "axis_y")
        all(number(r, "axis_x") == ax && number(r, "axis_y") == ay for r in chain) ||
            error("Inconsistent global axes: $file")
        # Project the actual saved x/y centers, not separately rounded t/u exports.
        x0, y0 = number(chain[1], "x_nm"), number(chain[1], "y_nm")
        dx = [number(r, "x_nm") - x0 for r in chain]
        dy = [number(r, "y_nm") - y0 for r in chain]
        t, u = ax .* dx .+ ay .* dy, -ay .* dx .+ ax .* dy
        all(>(0), diff(t)) || error("Non-increasing centers in global frame: $file")
        degree = min(settings.degree, length(chain) - 1)
        slopes = zeros(length(chain)); rms = 0.0
        if degree > 0
            centered = t .- mean(t)
            scale = maximum(abs, centered)
            xi = centered ./ scale
            design = hcat([xi .^ j for j in 0:degree]...)
            beta = design \ u
            slopes = fill(beta[2] / scale, length(chain))
            degree == 2 && (slopes .+= (2beta[3] / scale) .* xi)
            rms = sqrt(mean(abs2, design * beta - u)) / hypot(ax, ay)
        end
        for (r, slope) in zip(chain, slopes)
            isfinite(slope) && isfinite(rms) || error("Nonfinite centerline fit: $file")
            # Proper rotation, preserving the rounded saved axis norm exactly in
            # principle. Never change the model axes used for Gaussian subtraction.
            a, b = (ax - slope * ay) / hypot(1, slope), (ay + slope * ax) / hypot(1, slope)
            row = Dict(k => r[k] for k in GEOMETRY)
            merge!(row, Dict("file"=>file, "lobe"=>r["lobe"], "frame_axis_x"=>string(a),
                "frame_axis_y"=>string(b), "tangent_slope"=>string(slope), "centerline_rms_nm"=>string(rms)))
            push!(output, row)
        end
    end
    return output
end

function write_frames(path, rows)
    open(path, "w") do io
        println(io, join(HEADER, '\t'))
        for row in rows
            println(io, join([row[k] for k in HEADER], '\t'))
        end
    end
    return path
end

read_frames(::Nothing, rows) = nothing
function read_frames(path, rows)
    index = base_index(rows)
    lines = readlines(path)
    !isempty(lines) && split(lines[1], '\t') == HEADER || error("Unexpected patch-frame header")
    frames = Dict{Tuple{String,Int},Tuple{Float64,Float64}}()
    for line in lines[2:end]
        vals = split(line, '\t'; keepempty=true)
        length(vals) == length(HEADER) || error("Malformed patch-frame row")
        row = Dict(zip(HEADER, vals)); key = rowkey(row)
        !haskey(frames, key) && haskey(index, key) || error("Duplicate or unexpected patch-frame key: $key")
        all(number(row, k) == number(index[key], k) for k in GEOMETRY) || error("Frame/base geometry mismatch: $key")
        ax, ay = number(row, "axis_x"), number(row, "axis_y")
        a, b, slope = number(row, "frame_axis_x"), number(row, "frame_axis_y"), number(row, "tangent_slope")
        number(row, "centerline_rms_nm") >= 0 || error("Invalid centerline residual")
        expected = ((ax - slope * ay) / hypot(1, slope), (ay + slope * ax) / hypot(1, slope))
        # Arithmetic tolerance only, not a fitted rotation/quality threshold.
        hypot(a - expected[1], b - expected[2]) <= 256eps(Float64) * hypot(ax, ay) ||
            error("Frame is not the declared proper rotation: $key")
        frames[key] = (a, b)
    end
    Set(keys(frames)) == Set(keys(index)) || error("Incomplete patch-frame coverage")
    return frames
end

patch_axis(::Nothing, row, ax, ay) = (ax, ay)
patch_axis(frames, row, ax, ay) = frames[rowkey(row)]

end

module ScriptUtils

using Statistics

export _ensure_parent, _parse_f, _parse_ints, _parse_vec3, _read_key_tsv,
       _read_tsv, _standardize

function _ensure_parent(path::AbstractString)
    dir = dirname(String(path))
    isempty(dir) || mkpath(dir)
    return nothing
end

function _read_tsv(path::AbstractString)
    lines = readlines(path)
    data = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '#'), lines)
    isempty(data) && return String[], Dict{String,String}[]
    header = split(data[1], '\t'; keepempty=true)
    rows = Dict{String,String}[]
    for line in data[2:end]
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (i, h) in enumerate(header)
            row[h] = i <= length(vals) ? vals[i] : ""
        end
        push!(rows, row)
    end
    return header, rows
end

function _read_key_tsv(path::AbstractString)
    isfile(path) || error("Missing key TSV: $path")
    d = Dict{String,String}()
    for line in readlines(path)
        t = strip(line)
        isempty(t) && continue
        startswith(t, '#') && continue
        parts = split(t, '\t'; limit=2)
        length(parts) == 2 || continue
        d[strip(parts[1])] = strip(parts[2])
    end
    return d
end

function _parse_f(s)
    t = strip(String(s))
    (isempty(t) || t in ("NA", "NaN", "nan")) && return NaN
    return parse(Float64, t)
end

function _parse_vec3(s::AbstractString)
    vals = [parse(Float64, strip(x)) for x in split(s, ',') if !isempty(strip(x))]
    length(vals) == 3 || error("expected comma-separated x,y,z vector, got: $s")
    return vals
end

function _parse_ints(s::AbstractString)
    vals = [parse(Int, strip(x)) for x in split(s, ',') if !isempty(strip(x))]
    isempty(vals) && error("empty index list")
    return vals
end

function _standardize(v::Vector{Float64})
    good = filter(isfinite, v)
    isempty(good) && return fill(NaN, length(v))
    μ = mean(good)
    σ = std(good)
    σ = σ > 0 ? σ : 1.0
    return [isfinite(x) ? (x - μ) / σ : NaN for x in v]
end

end

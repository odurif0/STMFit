# Minimal self-contained TSV/path helpers (the module does not depend on
# script_utils.jl so it can be included from any caller).

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

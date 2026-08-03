export ProxySource, ProxyTemplate, ProxyEntry, ProxyEnsemble,
       RegistryConfigError, sha256_file, canonical_payload_hash

struct RegistryConfigError <: Exception
    msg::String
end

Base.showerror(io::IO, e::RegistryConfigError) = print(io, "RegistryConfigError: ", e.msg)

struct ProxyTemplate
    type::Int
    parity::Int
    mirror::Int
    pixels::Vector{Float64}
end

struct ProxySource
    family::String
    name::String
    path::String
    sha256::String
    height_nm::Float64
    weight::Float64
    present::Bool
end

struct ProxyEntry
    source::ProxySource
    templates::Vector{ProxyTemplate}
end

struct ProxyEnsemble
    entries::Vector{ProxyEntry}
    grid_half_nm::Float64
    grid_step_nm::Float64
    grid_n::Int
    npix::Int
    payload_sha256::String
    audit::Vector{String}
end

function sha256_file(path::AbstractString)
    isfile(path) || throw(RegistryConfigError("cannot hash missing file: $path"))
    open(path, "r") do io
        return bytes2hex(SHA.sha256(io))
    end
end

function _canonical_payload_string(entries::Vector{ProxyEntry}, grid_half_nm::Float64,
                                   grid_step_nm::Float64, npix::Int)
    parts = String[]
    push!(parts, @sprintf("grid_half_nm=%.10g", grid_half_nm))
    push!(parts, @sprintf("grid_step_nm=%.10g", grid_step_nm))
    push!(parts, @sprintf("npix=%d", npix))
    for entry in sort(entries; by=e -> (e.source.family, e.source.name))
        s = entry.source
        push!(parts, @sprintf("source|%s|%s|%.10g|%.10g|%s",
                              s.family, s.name, s.height_nm, s.weight, s.sha256))
        for t in sort(entry.templates; by=t -> (t.type, t.parity, t.mirror))
            pix = join((@sprintf("%.10g", p) for p in t.pixels), ',')
            push!(parts, @sprintf("tmpl|%d|%d|%d|%s", t.type, t.parity, t.mirror, pix))
        end
    end
    return join(parts, '\n')
end

function canonical_payload_hash(entries::Vector{ProxyEntry}, grid_half_nm::Float64,
                               grid_step_nm::Float64, npix::Int)
    s = _canonical_payload_string(entries, grid_half_nm, grid_step_nm, npix)
    return bytes2hex(SHA.sha256(codeunits(s)))
end

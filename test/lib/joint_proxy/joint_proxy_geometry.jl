function _apply_flip(t::Float64, u::Float64, spec::String)
    spec == "none" && return t, u
    spec == "t" && return -t, u
    spec == "u" && return t, -u
    spec == "both" && return -t, -u
    throw(RegistryConfigError("invalid flip spec: $spec"))
end

function _normalize(v::Vector{Float64}, method::String)
    method == "none" && return v
    if method == "sum"
        s = sum(abs, v)
        return s > 0 ? v ./ s : v
    elseif method == "max"
        m = maximum(abs.(v))
        return m > 0 ? v ./ m : v
    elseif method == "zscore"
        σ = std(v)
        σ = σ > 0 ? σ : 1.0
        return (v .- mean(v)) ./ σ
    end
    throw(RegistryConfigError("invalid normalize method: $method"))
end

function _load_atoms(path::AbstractString)
    _, rows = _read_tsv(path)
    required = ("type", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm")
    atoms = Dict{Int,Vector{NamedTuple}}(0 => NamedTuple[], 1 => NamedTuple[])
    for row in rows
        for k in required
            haskey(row, k) || throw(RegistryConfigError("geometric sites TSV missing column '$k' in $path"))
        end
        typ = tryparse(Int, row["type"])
        typ === nothing && throw(RegistryConfigError("geometric sites TSV: invalid type '$(row["type"])' in $path"))
        typ in (0, 1) || throw(RegistryConfigError("geometric sites TSV: type must be 0 or 1, got $typ in $path"))
        w = _parse_f(row["weight"])
        sigt = _parse_f(row["sigma_t_nm"])
        sigu = _parse_f(row["sigma_u_nm"])
        (!isfinite(w) || w <= 0) && throw(RegistryConfigError("geometric sites TSV: non-positive/non-finite weight in $path"))
        (!isfinite(sigt) || sigt <= 0) && throw(RegistryConfigError("geometric sites TSV: non-positive/non-finite sigma_t_nm in $path"))
        (!isfinite(sigu) || sigu <= 0) && throw(RegistryConfigError("geometric sites TSV: non-positive/non-finite sigma_u_nm in $path"))
        push!(atoms[typ], (atom=get(row, "atom", "atom"), t=_parse_f(row["t_nm"]), u=_parse_f(row["u_nm"]),
                           weight=w, sigt=max(sigt, eps(Float64)), sigu=max(sigu, eps(Float64))))
    end
    isempty(atoms[0]) && throw(RegistryConfigError("geometric sites TSV: no rows for type=0 (GlcN) in $path"))
    isempty(atoms[1]) && throw(RegistryConfigError("geometric sites TSV: no rows for type=1 (GlcNAc) in $path"))
    return atoms
end

function _render_one(atoms, coords, parity::Int, mirror::Int,
                     parity_flip::String, mirror_flip::String, normalize::String)
    vals = Float64[]
    transformed = NamedTuple[]
    for a in atoms
        t, u = a.t, a.u
        parity == 1 && ((t, u) = _apply_flip(t, u, parity_flip))
        mirror == 1 && ((t, u) = _apply_flip(t, u, mirror_flip))
        push!(transformed, (t=t, u=u, weight=a.weight, sigt=a.sigt, sigu=a.sigu))
    end
    for u in coords, t in coords
        v = 0.0
        for a in transformed
            v += a.weight * exp(-0.5 * (((t - a.t) / a.sigt)^2 + ((u - a.u) / a.sigu)^2))
        end
        push!(vals, v)
    end
    return _normalize(vals, normalize)
end

function render_geometric_templates(sites_tsv::AbstractString;
                                    half_nm::Float64=0.32, step_nm::Float64=0.08,
                                    parity_flip::String="t", mirror_flip::String="u",
                                    normalize::String="zscore")
    isfile(sites_tsv) || throw(RegistryConfigError("required geometric sites TSV not found: $sites_tsv"))
    parity_flip in ("none", "t", "u", "both") || throw(RegistryConfigError("invalid parity_flip: $parity_flip"))
    mirror_flip in ("none", "t", "u", "both") || throw(RegistryConfigError("invalid mirror_flip: $mirror_flip"))
    normalize in ("none", "sum", "max", "zscore") || throw(RegistryConfigError("invalid normalize: $normalize"))
    step_nm > 0 || throw(RegistryConfigError("step_nm must be positive"))
    half_nm > 0 || throw(RegistryConfigError("half_nm must be positive"))
    atoms = _load_atoms(sites_tsv)
    coords = collect(-half_nm:step_nm:half_nm)
    tmpls = ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        pix = _render_one(atoms[typ], coords, parity, mirror, parity_flip, mirror_flip, normalize)
        push!(tmpls, ProxyTemplate(typ, parity, mirror, pix))
    end
    return tmpls
end

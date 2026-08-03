function load_stm_template_tsv(path::AbstractString; npix::Int=81)
    isfile(path) || throw(RegistryConfigError("STM template TSV not found: $path"))
    header, rows = _read_tsv(path)
    for k in ("name", "type", "parity", "mirror")
        k in header || throw(RegistryConfigError("STM template TSV missing column '$k' in $path"))
    end
    pix_cols = [c for c in header if occursin(r"^p\d+$", c)]
    length(pix_cols) == npix || throw(RegistryConfigError("STM template TSV has $(length(pix_cols)) pixel columns, expected $npix in $path"))
    sort!(pix_cols; by=c -> parse(Int, match(r"\d+", c).match))

    seen = Set{Tuple{Int,Int,Int}}()
    tmpls = ProxyTemplate[]
    for (i, row) in enumerate(rows)
        typ = tryparse(Int, row["type"])
        parity = tryparse(Int, row["parity"])
        mirror = tryparse(Int, row["mirror"])
        (typ === nothing || parity === nothing || mirror === nothing) && throw(RegistryConfigError("STM template TSV row $i: invalid type/parity/mirror in $path"))
        typ in (0, 1) || throw(RegistryConfigError("STM template TSV row $i: type must be 0 or 1, got $typ in $path"))
        parity in (0, 1) || throw(RegistryConfigError("STM template TSV row $i: parity must be 0 or 1, got $parity in $path"))
        mirror in (0, 1) || throw(RegistryConfigError("STM template TSV row $i: mirror must be 0 or 1, got $mirror in $path"))
        key = (typ, parity, mirror)
        key in seen && throw(RegistryConfigError("STM template TSV: duplicate (type=$typ,parity=$parity,mirror=$mirror) in $path"))
        push!(seen, key)
        pix = Float64[]
        for c in pix_cols
            v = _parse_f(row[c])
            isfinite(v) || throw(RegistryConfigError("STM template TSV row $i: non-finite pixel in column $c in $path"))
            push!(pix, v)
        end
        push!(tmpls, ProxyTemplate(typ, parity, mirror, pix))
    end
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        (typ, parity, mirror) in seen || throw(RegistryConfigError("STM template TSV: missing (type=$typ,parity=$parity,mirror=$mirror) in $path"))
    end
    return tmpls
end

function _ensure_stm_specs(cfg)::Vector{Tuple{Float64,String,Float64}}
    haskey(cfg, "stm_prelim") || throw(RegistryConfigError("config missing [stm_prelim] section"))
    stm_cfg = cfg["stm_prelim"]
    family_weight = Float64(get(stm_cfg, "family_weight", 0.5))
    (family_weight > 0 && family_weight <= 1.0) || throw(RegistryConfigError("stm_prelim.family_weight must be in (0, 1]"))
    sources_cfg = get(stm_cfg, "sources", [])
    specs = Tuple{Float64,String,Float64}[]
    for src in sources_cfg
        height_nm = Float64(get(src, "height_nm", 0.0))
        path = String(get(src, "path", ""))
        w = Float64(get(src, "weight", 1.0))
        isempty(path) && throw(RegistryConfigError("stm_prelim source: path is empty"))
        (!isfinite(w) || w <= 0) && throw(RegistryConfigError("stm_prelim source: non-positive/non-finite weight for $path"))
        push!(specs, (height_nm, path, w))
    end
    return specs
end

function _required_string(cfg, key::String, section::String)
    value = String(get(cfg, key, ""))
    isempty(value) && throw(RegistryConfigError("$section.$key is required"))
    return value
end

function _require_hash(value::String, label::String)
    occursin(r"^[0-9a-f]{64}$", value) || throw(RegistryConfigError("$label must be a lowercase SHA-256"))
    return value
end

function _sidecar_value(sidecar, key::String, path::String)
    haskey(sidecar, key) || throw(RegistryConfigError("production provenance missing '$key' in $path"))
    return sidecar[key]
end

function _load_production_entry(cfg, npix::Int, half_nm::Float64,
                                step_nm::Float64, audit::Vector{String})
    prod = cfg["stm_dft_v1"]
    prelim = get(cfg, "stm_prelim", nothing)
    prelim === nothing && throw(RegistryConfigError("production config must retain [stm_prelim] with enabled=false"))
    get(prelim, "enabled", true) === false ||
        throw(RegistryConfigError("stm_prelim.enabled must be false when stm_dft_v1 is configured"))

    path = _required_string(prod, "path", "stm_dft_v1")
    maps_path = _required_string(prod, "maps_path", "stm_dft_v1")
    provenance_path = _required_string(prod, "provenance_path", "stm_dft_v1")
    for required_path in (path, maps_path, provenance_path)
        isfile(required_path) || throw(RegistryConfigError("required production artifact not found: $required_path"))
    end

    expected_template_sha = _require_hash(_required_string(prod, "sha256", "stm_dft_v1"), "stm_dft_v1.sha256")
    expected_maps_sha = _require_hash(_required_string(prod, "maps_sha256", "stm_dft_v1"), "stm_dft_v1.maps_sha256")
    expected_provenance_sha = _require_hash(_required_string(prod, "provenance_sha256", "stm_dft_v1"), "stm_dft_v1.provenance_sha256")
    expected_glcn_sha = _require_hash(_required_string(prod, "glcn_cube_sha256", "stm_dft_v1"), "stm_dft_v1.glcn_cube_sha256")
    expected_glcnac_sha = _require_hash(_required_string(prod, "glcnac_cube_sha256", "stm_dft_v1"), "stm_dft_v1.glcnac_cube_sha256")
    actual_template_sha = sha256_file(path)
    actual_maps_sha = sha256_file(maps_path)
    actual_provenance_sha = sha256_file(provenance_path)
    actual_template_sha == expected_template_sha || throw(RegistryConfigError("stm_dft_v1 template SHA-256 mismatch"))
    actual_maps_sha == expected_maps_sha || throw(RegistryConfigError("stm_dft_v1 maps SHA-256 mismatch"))
    actual_provenance_sha == expected_provenance_sha || throw(RegistryConfigError("stm_dft_v1 provenance SHA-256 mismatch"))

    height_nm = Float64(get(prod, "height_nm", NaN))
    sample_bias_ev = Float64(get(prod, "sample_bias_ev", NaN))
    cube_units = _required_string(prod, "cube_units", "stm_dft_v1")
    family_weight = Float64(get(prod, "family_weight", 0.5))
    isfinite(height_nm) || throw(RegistryConfigError("stm_dft_v1.height_nm is required"))
    isfinite(sample_bias_ev) || throw(RegistryConfigError("stm_dft_v1.sample_bias_ev is required"))
    (family_weight > 0 && family_weight < 1) || throw(RegistryConfigError("stm_dft_v1.family_weight must be in (0, 1)"))

    sidecar = TOML.parsefile(provenance_path)
    _sidecar_value(sidecar, "schema", provenance_path) == "stmfit-qe-mold-provenance-v1" || throw(RegistryConfigError("unsupported production provenance schema"))
    _sidecar_value(sidecar, "provider", provenance_path) == "stm_dft_v1" || throw(RegistryConfigError("production provenance provider mismatch"))
    _sidecar_value(sidecar, "templates_sha256", provenance_path) == expected_template_sha || throw(RegistryConfigError("production provenance template hash mismatch"))
    _sidecar_value(sidecar, "maps_sha256", provenance_path) == expected_maps_sha || throw(RegistryConfigError("production provenance maps hash mismatch"))
    _sidecar_value(sidecar, "glcn_cube_sha256", provenance_path) == expected_glcn_sha || throw(RegistryConfigError("production provenance GlcN cube hash mismatch"))
    _sidecar_value(sidecar, "glcnac_cube_sha256", provenance_path) == expected_glcnac_sha || throw(RegistryConfigError("production provenance GlcNAc cube hash mismatch"))
    isapprox(Float64(_sidecar_value(sidecar, "height_nm", provenance_path)), height_nm; atol=1e-12) || throw(RegistryConfigError("production provenance height mismatch"))
    isapprox(Float64(_sidecar_value(sidecar, "sample_bias_ev", provenance_path)), sample_bias_ev; atol=1e-9) || throw(RegistryConfigError("production provenance sample bias mismatch"))
    isapprox(Float64(_sidecar_value(sidecar, "half_nm", provenance_path)), half_nm; atol=1e-12) || throw(RegistryConfigError("production provenance grid half-size mismatch"))
    isapprox(Float64(_sidecar_value(sidecar, "step_nm", provenance_path)), step_nm; atol=1e-12) || throw(RegistryConfigError("production provenance grid spacing mismatch"))
    String(_sidecar_value(sidecar, "cube_units", provenance_path)) == cube_units || throw(RegistryConfigError("production provenance cube units mismatch"))

    templates = load_stm_template_tsv(path; npix=npix)
    source = ProxySource("stm_dft_v1", "stm_dft_m030_h050_v1", path,
                         actual_template_sha, height_nm, family_weight, true)
    push!(audit, "[OK] production stm_dft_v1 source: $path")
    push!(audit, "     sha256=$actual_template_sha provenance=$actual_provenance_sha")
    push!(audit, "[INFO] stm_prelim configured but disabled after production promotion")
    return ProxyEntry(source, templates), family_weight
end

function load_registry(config_path::AbstractString)
    isfile(config_path) || throw(RegistryConfigError("config file not found: $config_path"))
    cfg = TOML.parsefile(config_path)

    haskey(cfg, "grid") || throw(RegistryConfigError("config missing [grid] section"))
    grid = cfg["grid"]
    half_nm = Float64(get(grid, "half_nm", 0.32))
    step_nm = Float64(get(grid, "step_nm", 0.08))
    half_nm > 0 || throw(RegistryConfigError("grid.half_nm must be positive"))
    step_nm > 0 || throw(RegistryConfigError("grid.step_nm must be positive"))
    coords = collect(-half_nm:step_nm:half_nm)
    grid_n = length(coords)
    npix = grid_n^2

    haskey(cfg, "geometric") || throw(RegistryConfigError("config missing [geometric] section (required)"))
    geom_cfg = cfg["geometric"]
    sites_tsv = get(geom_cfg, "sites_tsv", "")
    isempty(sites_tsv) && throw(RegistryConfigError("geometric.sites_tsv is empty (required)"))
    isfile(sites_tsv) || throw(RegistryConfigError("required geometric sites TSV not found: $sites_tsv"))
    parity_flip = String(get(geom_cfg, "parity_flip", "t"))
    mirror_flip = String(get(geom_cfg, "mirror_flip", "u"))
    normalize = String(get(geom_cfg, "normalize", "zscore"))

    geom_tmpls = render_geometric_templates(sites_tsv; half_nm=half_nm, step_nm=step_nm,
        parity_flip=parity_flip, mirror_flip=mirror_flip, normalize=normalize)
    geom_sha = sha256_file(sites_tsv)

    audit = String[]
    push!(audit, "[OK] geometric source: $sites_tsv")
    push!(audit, "     sha256=$geom_sha")
    push!(audit, "     grid=$(grid_n)x$(grid_n) ($(npix) pixels)")

    entries = ProxyEntry[]
    if haskey(cfg, "stm_dft_v1")
        production, family_weight = _load_production_entry(cfg, npix, half_nm, step_nm, audit)
        push!(entries, ProxyEntry(ProxySource("geometric", "geometric", sites_tsv,
            geom_sha, 0.0, 1.0 - family_weight, true), geom_tmpls))
        push!(entries, production)
        push!(audit, "[INFO] family weights: geometric=$(Printf.@sprintf("%.4f", 1.0 - family_weight)) stm_dft_v1=$(Printf.@sprintf("%.4f", family_weight))")
    else
        specs = _ensure_stm_specs(cfg)
        family_weight = Float64(get(cfg["stm_prelim"], "family_weight", 0.5))
        stm_entries = ProxyEntry[]
        present_weights = Float64[]
        for (height_nm, path, cw) in specs
            if isfile(path)
                tmpls = load_stm_template_tsv(path; npix=npix)
                sha = sha256_file(path)
                name = @sprintf("stm_h%03d", round(Int, height_nm * 100))
                push!(stm_entries, ProxyEntry(ProxySource("stm_prelim", name, path, sha, height_nm, 0.0, true), tmpls))
                push!(present_weights, cw)
                push!(audit, "[OK] stm_prelim source: $path (h=$(Printf.@sprintf("%.2f", height_nm)) nm)")
                push!(audit, "     sha256=$sha")
            else
                push!(audit, "[SKIP] missing optional stm_prelim source: $path (h=$(Printf.@sprintf("%.2f", height_nm)) nm) — audited and omitted")
            end
        end
        if isempty(stm_entries)
        push!(entries, ProxyEntry(ProxySource("geometric", "geometric", sites_tsv, geom_sha, 0.0, 1.0, true), geom_tmpls))
        push!(audit, "[INFO] no STM sources present; geometric weight = 1.0")
        else
            geom_w = 1.0 - family_weight
            push!(entries, ProxyEntry(ProxySource("geometric", "geometric", sites_tsv, geom_sha, 0.0, geom_w, true), geom_tmpls))
            total_present = sum(present_weights)
            for (entry, cw) in zip(stm_entries, present_weights)
                w = family_weight * (cw / total_present)
                src = ProxySource(entry.source.family, entry.source.name, entry.source.path,
                                  entry.source.sha256, entry.source.height_nm, w, true)
                push!(entries, ProxyEntry(src, entry.templates))
            end
            push!(audit, "[INFO] family weights: geometric=$(Printf.@sprintf("%.4f", geom_w)) stm_prelim=$(Printf.@sprintf("%.4f", family_weight)) (renormalized over $(length(stm_entries)) present sources)")
        end
    end

    bond_flag = get(get(cfg, "provenance", Dict()), "include_bond_templates", false)
    bond_flag === false || throw(RegistryConfigError("bond templates are excluded in v1; set include_bond_templates=false"))
    push!(audit, "[INFO] bond/pair templates excluded (v1: concatenation bonds duplicate unary evidence)")

    payload_hash = canonical_payload_hash(entries, half_nm, step_nm, npix)
    push!(audit, "[OK] payload sha256=$payload_hash")
    return ProxyEnsemble(entries, half_nm, step_nm, grid_n, npix, payload_hash, audit)
end

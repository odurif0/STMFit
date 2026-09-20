module ReconstructedUnitAssignment

using Printf, Statistics, TOML

export read_table, write_table, lobe_table, require_same_keys,
       transverse_asymmetry, transverse_descriptor, augment_descriptor, write_soft_vote, load_config

const TRANSVERSE_DESCRIPTORS = ("transverse_half_plane_asymmetry",
    "transverse_first_moment", "affine_residual_half_plane_asymmetry")

function read_table(path::AbstractString)
    lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), '#'), readlines(path))
    isempty(lines) && error("Empty table: $path")
    header = String.(split(first(lines), '\t'; keepempty=true))
    length(unique(header)) == length(header) || error("Duplicate columns: $path")
    for col in header
        key = replace(lowercase(col), r"[^a-z0-9]" => "")
        (key in ("sequence", "controlsequence", "expectedn", "targetn", "manifest", "grade", "report") ||
         occursin("truth", key) || occursin("benchmark", key) || endswith(key, "sequence")) &&
            error("Forbidden benchmark/control column: $col")
    end
    rows = Dict{String,String}[]
    for (i, line) in enumerate(lines[2:end])
        vals = String.(split(line, '\t'; keepempty=true))
        length(vals) == length(header) || error("Malformed row $(i+1): $path")
        push!(rows, Dict(zip(header, vals)))
    end
    isempty(rows) && error("No data rows: $path")
    return header, rows
end

function write_table(path::AbstractString, header, rows)
    (ispath(path) || islink(path)) && error("Output already exists: $path")
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, join(header, '\t'))
        for row in rows
            vals = [string(get(row, col, "NA")) for col in header]
            any(v -> occursin('\t', v) || occursin('\n', v) || occursin('\r', v), vals) &&
                error("Invalid TSV field")
            println(io, join(vals, '\t'))
        end
    end
    return path
end

function lobe_table(path::AbstractString; required=String[], contiguous=true)
    header, rows = read_table(path)
    for col in vcat(["file", "lobe"], collect(required))
        col in header || error("$path missing required column: $col")
    end
    table = Dict{Tuple{String,Int},Dict{String,String}}()
    by_file = Dict{String,Set{Int}}()
    for row in rows
        file = basename(strip(row["file"]))
        n = tryparse(Int, row["lobe"])
        !isempty(file) && n !== nothing && n > 0 || error("Invalid file/lobe key: $path")
        key = (file, n)
        haskey(table, key) && error("Duplicate file/lobe key $key: $path")
        row["file"] = file
        table[key] = row
        push!(get!(by_file, file, Set{Int}()), n)
    end
    if contiguous
        for (file, lobes) in by_file
            lobes == Set(1:maximum(lobes)) || error("Noncontiguous lobes: $file")
        end
    end
    return header, table
end

function require_same_keys(reference, other, label)
    missing = setdiff(Set(keys(reference)), Set(keys(other)))
    extra = setdiff(Set(keys(other)), Set(keys(reference)))
    isempty(missing) && isempty(extra) ||
        error("$label coverage mismatch: $(length(missing)) missing, $(length(extra)) extra")
    return nothing
end

function load_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    for section in ("model", "selection", "preprocessing")
        haskey(cfg, section) || error("Missing config section [$section]")
    end
    model, pre = cfg["model"], cfg["preprocessing"]
    model["descriptor"] in TRANSVERSE_DESCRIPTORS || error("Unsupported descriptor")
    model["descriptor_column"] == "patch_u_asym_reconstructed" || error("Unsupported descriptor column")
    model["descriptor_channel"] == "bwd_res" || error("Unsupported descriptor patch family")
    pre["pixel_order"] == "u_outer_t_inner" || error("Unsupported patch ordering")
    pre["descriptor_normalization"] == "median_sample_std" || error("Unsupported patch normalization")
    get(pre, "patch_residual_filter", nothing) in ("smooth_data_only", "smooth_residual") ||
        error("Explicit patch_residual_filter must be smooth_data_only or smooth_residual")
    model["descriptor_half_nm"] > 0 && model["descriptor_step_nm"] > 0 || error("Invalid descriptor grid")
    side = round(Int, 2model["descriptor_half_nm"] / model["descriptor_step_nm"]) + 1
    side == 9 || error("Reconstructed descriptor requires a 9x9 patch")
    isfinite(model["descriptor_zero_l1"]) && model["descriptor_zero_l1"] >= 0 || error("Invalid zero-signal tolerance")
    sel = cfg["selection"]
    for key in ("kmeans_seeds", "gmm_seeds")
        sel[key] isa Integer && !(sel[key] isa Bool) && sel[key] > 0 || error("Invalid $key")
    end
    sel["first_seed"] >= 0 || error("Invalid first_seed")
    0 <= sel["vote_threshold"] <= 1 || error("Invalid vote threshold")
    return cfg
end

_float(s) = something(tryparse(Float64, strip(String(s))), NaN)

"""Transverse half-plane asymmetry of an already-normalized patch.

Input ordering is the extractor's `for u in coords, t in coords`: t varies
fastest. The center row (u=0) contributes zero. No second normalization is
applied. Nonfinite samples or negligible L1 mass return (NaN, reason).
"""
function transverse_asymmetry(values::AbstractVector, coords::AbstractVector; zero_l1::Real)
    length(values) == length(coords)^2 || error("Patch/grid size mismatch")
    all(isfinite, coords) && issorted(coords) || error("Invalid patch coordinates")
    zero_l1 >= 0 && isfinite(zero_l1) || error("Invalid zero-signal tolerance")
    all(isfinite, values) || return (NaN, "nonfinite_patch")
    mass = sum(abs, values)
    mass > zero_l1 || return (NaN, "zero_patch_mass")
    result = 0.0
    i = 1
    for u in coords, _ in coords
        result += sign(u) * values[i]
        i += 1
    end
    return (result / mass, "ok")
end

"""Remove the least-squares affine plane on the complete square patch grid.

Only the descriptor uses this projection: raw pixels, fit geometry, CC/Fisher
patches and the k-means features remain untouched. The affine part of any true
molecular signal is removed too; this is not a calibrated background estimate.
"""
function _affine_residual_patch(values::AbstractVector, coords::AbstractVector)
    n = length(coords)
    length(values) == n^2 || error("Patch/grid size mismatch")
    all(isfinite, values) && all(isfinite, coords) && issorted(coords) || error("Invalid affine patch/grid")
    centered_coords = Float64.(coords .- mean(coords))
    u = repeat(centered_coords; inner=n)
    t = repeat(centered_coords; outer=n)
    coordinate_energy = n * sum(abs2, centered_coords)
    coordinate_energy > 0 || error("Affine patch requires a nonzero coordinate extent")
    slope_u = sum(u .* values) / coordinate_energy
    slope_t = sum(t .* values) / coordinate_energy
    return values .- mean(values) .- slope_u .* u .- slope_t .* t
end

"""Explicit descriptor choice on the same normalized u-outer/t-inner patch.

The first moment uses u/max(abs(u)); the affine variant applies the existing
half-plane formula to the least-squares plane residual, including its L1 norm.
None of these definitions recovers the unavailable historical producer.
"""
function transverse_descriptor(values::AbstractVector, coords::AbstractVector,
                               kind::AbstractString; zero_l1::Real)
    kind in TRANSVERSE_DESCRIPTORS || error("Unsupported descriptor: $kind")
    legacy = transverse_asymmetry(values, coords; zero_l1)
    kind == "transverse_half_plane_asymmetry" && return legacy
    last(legacy) == "ok" || return legacy
    if kind == "transverse_first_moment"
        extent = maximum(abs, coords)
        extent > 0 || error("First moment requires a nonzero coordinate extent")
        result = 0.0
        i = 1
        for u in coords, _ in coords
            result += (u / extent) * values[i]
            i += 1
        end
        return (result / sum(abs, values), "ok")
    end
    residual = _affine_residual_patch(values, coords)
    value, reason = transverse_asymmetry(residual, coords; zero_l1)
    return (value, reason == "zero_patch_mass" ? "zero_affine_residual_mass" : reason)
end

function augment_descriptor(features, patches, out, config)
    cfg = load_config(config)
    model = cfg["model"]
    header, base = lobe_table(features)
    col = model["descriptor_column"]
    col in header && error("Descriptor already present: $col")
    pixels = [@sprintf("bwd_res_p%03d", i) for i in 1:81]
    pheader, patch = lobe_table(patches; required=pixels)
    length(filter(c -> startswith(c, "bwd_res_p"), pheader)) == 81 || error("Expected exactly 81 residual pixels")
    require_same_keys(base, patch, "Descriptor patches")
    # Symmetric construction makes the central coordinate exactly zero.
    coords = collect(-4:4) .* model["descriptor_step_nm"]
    outrows = Dict{String,String}[]
    for key in sort(collect(keys(base)))
        row = copy(base[key])
        value, reason = transverse_descriptor([_float(patch[key][c]) for c in pixels], coords,
                                             model["descriptor"]; zero_l1=model["descriptor_zero_l1"])
        row[col] = isfinite(value) ? @sprintf("%.17g", value) : "NA"
        row["descriptor_reason"] = reason
        push!(outrows, row)
    end
    write_table(out, vcat(header, [col, "descriptor_reason"]), outrows)
end

function write_soft_vote(features, kmeans, gmm, out, config)
    cfg = load_config(config)
    _, base = lobe_table(features)
    components = [last(lobe_table(path; required=["predicted", "probability_1"], contiguous=false))
                  for path in (kmeans, gmm)]
    for rows in components
        isempty(setdiff(Set(keys(rows)), Set(keys(base)))) || error("Unexpected component prediction keys")
    end
    output = Dict{String,String}[]
    for key in sort(collect(keys(base)))
        ps = Float64[]
        reasons = String[]
        for (name, rows) in zip(("kmeans", "gmm"), components)
            row = get(rows, key, nothing)
            if row === nothing
                push!(reasons, "missing_" * name)
                continue
            end
            row["predicted"] in ("0", "1", "?") || error("Invalid component prediction")
            p = _float(row["probability_1"])
            if row["predicted"] == "?" || !isfinite(p)
                push!(reasons, "unavailable_" * name)
            else
                0 <= p <= 1 || error("Component probability outside [0,1]")
                push!(ps, p)
            end
        end
        row = Dict("file" => key[1], "lobe" => string(key[2]), "model" => cfg["model"]["name"])
        if isempty(reasons)
            p = sum(ps) / 2
            row["predicted"] = p >= cfg["selection"]["vote_threshold"] ? "1" : "0"
            row["confidence"] = @sprintf("%.8f", abs(p - 0.5) * 2)
            row["probability_1"] = @sprintf("%.8f", p)
            row["invalid_reason"] = "ok"
        else
            row["predicted"] = "?"
            row["confidence"] = "0.00000000"
            row["probability_1"] = "NA"
            row["invalid_reason"] = join(reasons, ",")
        end
        push!(output, row)
    end
    write_table(out, ["file", "lobe", "predicted", "confidence", "probability_1", "invalid_reason", "model"], output)
end

end

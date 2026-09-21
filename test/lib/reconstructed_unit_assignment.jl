module ReconstructedUnitAssignment

using Printf, Statistics, TOML

export read_table, write_table, lobe_table, require_same_keys,
       transverse_asymmetry, transverse_descriptor, augment_descriptor, mold_margin,
       write_soft_vote, load_config, load_training_policy, write_training_support,
       load_training_mask, validate_training_mask, load_gmm_normalization, load_gmm_weighting,
       load_gmm_seed_aggregation

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

"Explicit fitting eligibility, independent of scoring support and benchmark labels."
function load_training_policy(config::AbstractDict)
    policy = get(get(config, "selection", Dict()), "assignment_training_support", nothing)
    policy in ("all_admissible", "complete_patches") ||
        throw(ArgumentError("explicit assignment_training_support must be all_admissible or complete_patches"))
    return String(policy)
end

const TRAINING_PATCH_FAMILIES = (("fwd17", "res_p", 17), ("bwd17", "bwd_res_p", 17),
                               ("bwd9", "bwd_res_p", 9))

"Record observed pixel counts only; never impute pixels or select by a prediction."
function write_training_support(fwd17, bwd17, bwd9, output)
    tables = []
    for (path, (_, prefix, side)) in zip((fwd17, bwd17, bwd9), TRAINING_PATCH_FAMILIES)
        _, table = lobe_table(path; required=[prefix * lpad(string(i), 3, '0') for i in 1:side^2])
        isempty(tables) || require_same_keys(first(tables), table, "training patch support")
        push!(tables, table)
    end
    rows = Dict{String,String}[]
    for key in sort(collect(keys(first(tables))))
        row = Dict("file" => key[1], "lobe" => string(key[2]))
        for (table, (name, prefix, side)) in zip(tables, TRAINING_PATCH_FAMILIES)
            row[name * "_observed"] = string(count(i -> isfinite(_float(table[key][prefix * lpad(string(i), 3, '0')])), 1:side^2))
        end
        push!(rows, row)
    end
    return write_table(output, vcat(["file", "lobe"], [name * "_observed" for (name, _, _) in TRAINING_PATCH_FAMILIES]), rows)
end

function load_training_mask(policy, path, rowkeys, stage)
    stage in ("fisher", "gmm") || throw(ArgumentError("unknown training stage"))
    if policy == "all_admissible"
        isempty(path) || throw(ArgumentError("training support is only used with complete_patches"))
        return nothing
    end
    policy == "complete_patches" || throw(ArgumentError("unknown training policy"))
    isempty(path) && throw(ArgumentError("complete_patches requires --training-support"))
    columns = [name * "_observed" for (name, _, _) in TRAINING_PATCH_FAMILIES]
    _, table = lobe_table(path; required=columns, contiguous=false)
    length(unique(rowkeys)) == length(rowkeys) || error("Duplicate prediction keys")
    require_same_keys(Dict(k => nothing for k in rowkeys), table, "training support")
    mask = trues(length(rowkeys))
    for (i, key) in enumerate(rowkeys)
        complete = Bool[]
        for (name, _, side) in TRAINING_PATCH_FAMILIES
            observed = tryparse(Int, table[key][name * "_observed"])
            observed !== nothing && 0 <= observed <= side^2 || error("Invalid observed pixel count: $key, $name")
            push!(complete, observed == side^2)
        end
        mask[i] = stage == "fisher" ? first(complete) : all(complete)
    end
    return mask
end

function validate_training_mask(policy, mask, n)
    if policy == "all_admissible"
        mask === nothing || throw(ArgumentError("all_admissible must not filter training rows"))
        return trues(n)
    end
    policy == "complete_patches" || throw(ArgumentError("unknown training policy"))
    mask isa AbstractVector{Bool} && length(mask) == n ||
        throw(ArgumentError("complete_patches requires one Boolean training flag per row"))
    return mask
end

"Explicit per-file GMM feature scaling, independent of pixel normalization."
function load_gmm_normalization(config::AbstractDict)
    pre = get(config, "preprocessing", Dict())
    mode = get(pre, "gmm_feature_normalization", nothing)
    mode in ("mean_sample_std", "median_iqr") ||
        throw(ArgumentError("explicit gmm_feature_normalization must be mean_sample_std or median_iqr"))
    fallback = get(pre, "gmm_scale_fallback", nothing)
    fallback isa Real && !(fallback isa Bool) && isfinite(fallback) && fallback > 0 ||
        throw(ArgumentError("explicit gmm_scale_fallback must be positive and finite"))
    return (mode=String(mode), scale_fallback=Float64(fallback))
end

"Observation weights depend only on the number of usable training rows per scan."
function load_gmm_weighting(config::AbstractDict)
    mode = get(get(config, "selection", Dict()), "gmm_training_weighting", nothing)
    mode in ("equal_lobes", "equal_scans") ||
        throw(ArgumentError("explicit gmm_training_weighting must be equal_lobes or equal_scans"))
    mode == "equal_scans" && get(get(config, "model", Dict()), "gmm_final_covariance", nothing) != "ridge" &&
        throw(ArgumentError("equal_scans requires ridge; weighted shrinkage is not implemented"))
    return String(mode)
end

"Aggregation of existing per-seed scores; no change to fitting or group naming."
function load_gmm_seed_aggregation(config::AbstractDict)
    mode = get(get(config, "selection", Dict()), "gmm_seed_aggregation", nothing)
    mode in ("hard_vote", "mean_membership") ||
        throw(ArgumentError("explicit gmm_seed_aggregation must be hard_vote or mean_membership"))
    return String(mode)
end

function load_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    for section in ("model", "selection", "preprocessing")
        haskey(cfg, section) || error("Missing config section [$section]")
    end
    model, pre = cfg["model"], cfg["preprocessing"]
    load_gmm_normalization(cfg)
    load_gmm_weighting(cfg)
    load_gmm_seed_aggregation(cfg)
    model["descriptor"] in TRANSVERSE_DESCRIPTORS || error("Unsupported descriptor")
    model["descriptor_column"] == "patch_u_asym_reconstructed" || error("Unsupported descriptor column")
    model["descriptor_channel"] == "bwd_res" || error("Unsupported descriptor patch family")
    get(model, "fisher_score_center", nothing) in ("legacy_centered_mean", "training_mean") ||
        error("Explicit fisher_score_center must be legacy_centered_mean or training_mean")
    get(model, "gmm_final_covariance", nothing) in ("ridge", "ledoit_wolf") ||
        error("Explicit gmm_final_covariance must be ridge or ledoit_wolf")
    get(model, "gmm_final_score", nothing) in ("mahalanobis", "gaussian_density") ||
        error("Explicit gmm_final_score must be mahalanobis or gaussian_density")
    ridge = get(model, "gmm_covariance_ridge", nothing)
    ridge isa Real && !(ridge isa Bool) && isfinite(ridge) && ridge > 0 ||
        error("Explicit gmm_covariance_ridge must be positive and finite")
    pre["pixel_order"] == "u_outer_t_inner" || error("Unsupported patch ordering")
    pre["descriptor_normalization"] == "median_sample_std" || error("Unsupported patch normalization")
    get(model, "mold_margin_mode", nothing) in ("absolute_cost_margin", "signed_cost_difference") ||
        error("Explicit mold_margin_mode must be absolute_cost_margin or signed_cost_difference")
    get(pre, "fisher_patch_projection", nothing) in ("none", "affine_disk") ||
        error("Explicit fisher_patch_projection must be none or affine_disk")
    get(pre, "patch_residual_filter", nothing) in ("smooth_data_only", "smooth_residual") ||
        error("Explicit patch_residual_filter must be smooth_data_only or smooth_residual")
    get(pre, "assignment_patch_support", nothing) in ("full_square", "complete_disk_symmetric") ||
        error("Explicit assignment_patch_support must be full_square or complete_disk_symmetric")
    pre["assignment_patch_support"] == "complete_disk_symmetric" &&
        model["descriptor"] != "affine_residual_half_plane_asymmetry" &&
        error("Partial patch support is defined only for the affine residual descriptor")
    model["descriptor_half_nm"] > 0 && model["descriptor_step_nm"] > 0 || error("Invalid descriptor grid")
    side = round(Int, 2model["descriptor_half_nm"] / model["descriptor_step_nm"]) + 1
    side == 9 || error("Reconstructed descriptor requires a 9x9 patch")
    isfinite(model["descriptor_zero_l1"]) && model["descriptor_zero_l1"] >= 0 || error("Invalid zero-signal tolerance")
    sel = cfg["selection"]
    load_training_policy(cfg)
    model["gmm_final_covariance"] == "ledoit_wolf" && get(sel, "gmm_selftrain", 0) < 1 &&
        error("Final covariance shrinkage requires at least one hard self-training iteration")
    model["gmm_final_score"] == "gaussian_density" && get(sel, "gmm_selftrain", 0) < 1 &&
        error("Explicit final Gaussian score requires at least one hard self-training iteration")
    for key in ("kmeans_seeds", "gmm_seeds")
        sel[key] isa Integer && !(sel[key] isa Bool) && sel[key] > 0 || error("Invalid $key")
    end
    sel["first_seed"] >= 0 || error("Invalid first_seed")
    0 <= sel["vote_threshold"] <= 1 || error("Invalid vote threshold")
    return cfg
end

_float(s) = something(tryparse(Float64, strip(String(s))), NaN)

"""One fixed mold feature; experimental labels and decoded mold labels are unused.

Legacy mode reads the saved absolute margin without recomputing its rounding.
Signed mode is cost(GlcN)-cost(GlcNAc): positive favors the GlcNAc template under
the existing label-free geometric alignment. The sign is not calibrated truth.
"""
function mold_margin(row::AbstractDict, mode::AbstractString)
    mode == "absolute_cost_margin" && return _float(row["cost_margin"])
    mode == "signed_cost_difference" || error("Unsupported mold margin mode: $mode")
    c0, c1 = _float(row["cost_GlcN"]), _float(row["cost_GlcNAc"])
    return isfinite(c0) && isfinite(c1) ? c0 - c1 : NaN
end

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

"""Observed support closed under both grid reflections, with a complete disk.

The disk is defined on integer grid coordinates, not a fitted coverage fraction.
If any disk pixel is missing the patch remains unavailable. Outside the disk,
keep an observed pixel only if its entire (u,t) reflection orbit is observed.
No value is imputed. The retained design has full rank and no one-sided mask.
"""
function _symmetric_patch_support(values::AbstractVector, coords::AbstractVector)
    n = length(coords)
    n >= 3 && isodd(n) && length(values) == n^2 || error("Odd square patch required")
    half = n ÷ 2
    step = coords[half+2]
    step > 0 && coords == collect(-half:half) .* step || error("Symmetric uniform grid required")
    observed = isfinite.(values)
    for u in -half:half, t in -half:half
        u*u + t*t <= half^2 && !observed[(u+half)*n+t+half+1] &&
            return falses(n^2), "incomplete_descriptor_disk"
    end
    keep = copy(observed)
    for row in 1:n, col in 1:n
        keep[(row-1)*n+col] = observed[(row-1)*n+col] && observed[(n-row)*n+col] &&
            observed[(row-1)*n+n+1-col] && observed[(n-row)*n+n+1-col]
    end
    return keep, "ok_masked_symmetric"
end

function _masked_affine_descriptor(values, coords; zero_l1)
    keep, reason = _symmetric_patch_support(values, coords)
    reason == "ok_masked_symmetric" || return (NaN, reason)
    observed = values[keep]
    sum(abs, observed) > zero_l1 || return (NaN, "zero_patch_mass")
    n = length(coords)
    u, t = repeat(coords; inner=n)[keep], repeat(coords; outer=n)[keep]
    # Reflection closure makes [1,u,t] orthogonal on the observed support.
    residual = observed .- mean(observed) .-
        (sum(u .* observed) / sum(abs2, u)) .* u .-
        (sum(t .* observed) / sum(abs2, t)) .* t
    mass = sum(abs, residual)
    mass > zero_l1 || return (NaN, "zero_affine_residual_mass")
    return (sum(sign.(u) .* residual) / mass, reason)
end

"""Explicit descriptor choice on the same normalized u-outer/t-inner patch.

The first moment uses u/max(abs(u)); the affine variant applies the existing
half-plane formula to the least-squares plane residual, including its L1 norm.
None of these definitions recovers the unavailable historical producer.
"""
function transverse_descriptor(values::AbstractVector, coords::AbstractVector,
                               kind::AbstractString; zero_l1::Real, patch_support::AbstractString="full_square")
    kind in TRANSVERSE_DESCRIPTORS || error("Unsupported descriptor: $kind")
    patch_support in ("full_square", "complete_disk_symmetric") || error("Unsupported patch support")
    patch_support == "complete_disk_symmetric" && kind != "affine_residual_half_plane_asymmetry" &&
        error("Partial support requires the affine residual descriptor")
    legacy = transverse_asymmetry(values, coords; zero_l1)
    if patch_support == "complete_disk_symmetric" && last(legacy) == "nonfinite_patch"
        return _masked_affine_descriptor(values, coords; zero_l1)
    end
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
            model["descriptor"]; zero_l1=model["descriptor_zero_l1"],
            patch_support=cfg["preprocessing"]["assignment_patch_support"])
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

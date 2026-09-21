"""
Saved-patch representation arithmetic only. No fitting engine, classifier,
feature search, benchmark, or prediction reader is called. OLS here projects
serialized, already-normalized patches onto the fixed affine pixel basis;
it is not a new STM fit or a proposed replacement production feature.
"""
module ReconstructedRepresentationDiagnostics

using LinearAlgebra, Statistics, Printf, SHA
include(joinpath(@__DIR__, "reconstructed_unit_assignment.jl"))
include(joinpath(@__DIR__, "empirical_fisher_native.jl"))
using .ReconstructedUnitAssignment: load_config, lobe_table, require_same_keys,
    transverse_asymmetry, write_table
using .EmpiricalFisherNative: load_fisher_config, fisher_grid, flip_u_disk

export patch_audit, mass_parts, affine_components, response_decomposition,
       standardize_columns, negative_moments, run_audit, load_inputs

const BASE = ["amp_prominence", "amp_neighbor_ratio", "integrated_prominence", "amp_rel"]
const GMM_FEATURES = vcat(BASE, ["patch_u_asym_reconstructed", "mold_cc_fwd", "mold_cc_bwd", "emp_fisher"])
const ALL_FEATURES = vcat(GMM_FEATURES, ["split_log_skew", "bwd_neg_com_t", "bwd_neg_diag45"])
const PATCH_METRICS = ["mass_l1", "numerator", "d_raw", "mass_aligned", "mass_opposed",
    "mass_central", "fraction_aligned", "fraction_opposed", "fraction_central",
    "mass_identity_error", "d_identity_error", "one_minus_d_identity_error",
    "affine_constant", "gradient_t_per_nm", "gradient_u_per_nm",
    "centered_energy", "gradient_t_energy_fraction", "gradient_u_energy_fraction",
    "affine_residual_energy_fraction", "numerator_constant", "numerator_t",
    "numerator_u", "numerator_affine_residual", "numerator_attribution_error",
    "constant_over_raw_mass", "t_over_raw_mass", "u_over_raw_mass",
    "affine_residual_over_raw_mass", "u_share_of_numerator",
    "d_gradient_removed_diagnostic", "d_affine_removed_diagnostic",
    "gradient_removed_mass_l1", "affine_removed_mass_l1",
    "common_file_gradient_u_per_nm", "common_file_u_over_raw_mass",
    "common_file_u_share_of_numerator", "d_common_file_gradient_removed_diagnostic",
    "common_cohort_gradient_u_per_nm", "common_cohort_u_over_raw_mass",
    "common_cohort_u_share_of_numerator", "d_common_cohort_gradient_removed_diagnostic"]

number(x) = something(tryparse(Float64, strip(String(x))), NaN)
text(x::AbstractFloat) = isfinite(x) ? @sprintf("%.17g", x) : "NA"
text(x) = string(x)
ratio(x, y) = isfinite(y) && y != 0 ? x / y : NaN
function finite_summary(values)
    x = filter(isfinite, Float64.(collect(values)))
    isempty(x) && return (n=0, mean=NaN, std=NaN, min=NaN, q25=NaN, median=NaN, q75=NaN, max=NaN)
    return (n=length(x), mean=mean(x), std=length(x) > 1 ? std(x) : NaN,
            min=minimum(x), q25=quantile(x, 0.25), median=median(x), q75=quantile(x, 0.75), max=maximum(x))
end
function grid_vectors(coords)
    all(isfinite, coords) && issorted(coords) && length(coords) >= 3 || error("Invalid diagnostic grid")
    coords == -reverse(coords) || error("Diagnostic grid must be symmetric with exact central zero")
    t = [t for _ in coords for t in coords]
    u = [u for u in coords for _ in coords]
    return t, u
end

"Absolute mass identity, using the actual serialized u-outer/t-inner ordering."
function mass_parts(p, coords; zero_l1)
    t, u = grid_vectors(coords)
    length(p) == length(u) || error("Patch/grid size mismatch")
    d, reason = transverse_asymmetry(p, coords; zero_l1)
    reason == "ok" || return (; d, reason, mass=all(isfinite, p) ? sum(abs, p) : NaN,
                              numerator=NaN, aligned=NaN, opposed=NaN, central=NaN)
    signed = sign.(u) .* p
    return (; d, reason, mass=sum(abs, p), numerator=sum(signed),
            aligned=sum(abs(p[i]) for i in eachindex(p) if signed[i] > 0; init=0.0),
            opposed=sum(abs(p[i]) for i in eachindex(p) if signed[i] < 0; init=0.0),
            central=sum(abs(p[i]) for i in eachindex(p) if u[i] == 0; init=0.0))
end

"Orthogonal full-grid least-squares projection on [1,t,u], in normalized-pixel/nm units."
function affine_components(p, coords)
    t, u = grid_vectors(coords)
    length(p) == length(u) || error("Patch/grid size mismatch")
    all(isfinite, p) || error("Affine projection requires complete finite patch")
    c = mean(p)
    centered = p .- c
    bt = dot(t, centered) / dot(t, t)
    bu = dot(u, centered) / dot(u, u)
    constant = fill(c, length(p))
    along_t, along_u = bt .* t, bu .* u
    residual = p .- constant .- along_t .- along_u
    return (; c, bt, bu, constant, along_t, along_u, residual,
            centered_energy=sum(abs2, centered))
end

function patch_audit(p, coords; zero_l1::Real, identity_atol::Real)
    isfinite(identity_atol) && identity_atol >= 0 || error("Invalid identity tolerance")
    m = mass_parts(p, coords; zero_l1)
    r = Dict{String,Any}(k => NaN for k in PATCH_METRICS)
    r["reason"] = m.reason
    r["finite_pixels"] = count(isfinite, p)
    r["total_pixels"] = length(p)
    r["mass_l1"] = m.mass
    r["gradient_removed_reason"] = m.reason
    r["affine_removed_reason"] = m.reason
    r["common_file_removed_reason"] = m.reason
    r["common_cohort_removed_reason"] = m.reason
    m.reason == "ok" || return r
    r["numerator"], r["d_raw"] = m.numerator, m.d
    for (label, value) in (("aligned", m.aligned), ("opposed", m.opposed), ("central", m.central))
        r["mass_" * label] = value
        r["fraction_" * label] = value / m.mass
    end
    r["mass_identity_error"] = (m.aligned + m.opposed + m.central - m.mass) / m.mass
    r["d_identity_error"] = m.d - (m.aligned - m.opposed) / m.mass
    r["one_minus_d_identity_error"] = (1 - m.d) - (2m.opposed + m.central) / m.mass
    t, u = grid_vectors(coords)
    signs = sign.(u)
    a = affine_components(p, coords)
    r["affine_constant"], r["gradient_t_per_nm"], r["gradient_u_per_nm"] = a.c, a.bt, a.bu
    r["centered_energy"] = a.centered_energy
    r["gradient_t_energy_fraction"] = ratio(sum(abs2, a.along_t), a.centered_energy)
    r["gradient_u_energy_fraction"] = ratio(sum(abs2, a.along_u), a.centered_energy)
    r["affine_residual_energy_fraction"] = ratio(sum(abs2, a.residual), a.centered_energy)
    for (label, component) in (("constant", a.constant), ("t", a.along_t),
                               ("u", a.along_u), ("affine_residual", a.residual))
        numerator = dot(signs, component)
        r["numerator_" * label] = numerator
        r[label * "_over_raw_mass"] = numerator / m.mass
    end
    r["numerator_attribution_error"] =
        (sum(r["numerator_" * label] for label in ("constant", "t", "u", "affine_residual")) - m.numerator) / m.mass
    r["u_share_of_numerator"] = ratio(r["numerator_u"], m.numerator)
    for (label, component) in (("gradient", p .- a.along_u), ("affine", a.residual))
        value, reason = transverse_asymmetry(component, coords; zero_l1)
        r["d_" * label * "_removed_diagnostic"] = value
        r[label * "_removed_reason"] = reason
        r[label * "_removed_mass_l1"] = sum(abs, component)
    end
    for label in ("mass_identity_error", "d_identity_error", "one_minus_d_identity_error", "numerator_attribution_error")
        abs(r[label]) <= identity_atol || error("Arithmetic identity exceeds tolerance: $label=$(r[label])")
    end
    return r
end

"""
Exact max-reflection response identity, with the SAME centering offset on both
branches: max(w⋅(x-mid), w⋅(R*x-mid)) = w⋅(x+R*x)/2 - w⋅mid + |w⋅(x-R*x)/2|.
No re-centering of the reflected midpoint. For zero-padded legacy disk support,
R may not be an involution: these are half-sum/difference responses, not necessarily
parity eigenvectors. No model is trained and no real score is created here.
"""
function response_decomposition(x, w, mid, grid)
    length(x) == length(w) == length(mid) == length(grid.disk_indices) || error("Fisher length mismatch")
    all(isfinite, x) && all(isfinite, w) && all(isfinite, mid) || error("Nonfinite response input")
    reflected = flip_u_disk(x, grid)
    offset = -dot(mid, w)
    even_uncentered = dot((x + reflected) / 2, w)
    odd = dot((x - reflected) / 2, w)
    original, mirrored = dot(x - mid, w), dot(reflected - mid, w)
    return (; original, mirrored, shared_centering_offset=offset, even_uncentered,
            even_centered=even_uncentered + offset, odd, reflection_bonus=abs(odd),
            analytic_max=even_uncentered + offset + abs(odd), direct_max=max(original, mirrored))
end

"Exact existing predictor preprocessing, without fitting or predicting."
function standardize_columns(raw::AbstractMatrix, files; interactions::Bool)
    n, p = size(raw)
    length(files) == n || error("Wrong file count")
    z = fill(NaN, n, p)
    params = Dict{Tuple{String,Int},NamedTuple}()
    for file in sort(unique(files))
        idx = findall(==(file), files)
        for j in 1:p
            good = filter(isfinite, raw[idx, j])
            μ = isempty(good) ? NaN : mean(good)
            σ = isempty(good) ? NaN : std(good)
            scale = σ > 0 ? σ : 1.0
            params[(file, j)] = (n=length(good), mean=μ, sample_std=σ, scale=scale)
            for i in idx
                z[i, j] = isfinite(raw[i, j]) ? (raw[i, j] - μ) / scale : NaN
            end
        end
    end
    expanded = copy(z)
    if interactions && p > 1
        expanded = hcat(z, [z[:, a] .* z[:, b] for a in 1:(p-1) for b in (a+1):p]...)
    end
    valid = vec(all(isfinite, raw; dims=2)) .& vec(all(isfinite, expanded; dims=2))
    return (; z, expanded, valid, params)
end

"Existing backward negative-moment formulas, not a new descriptor."
function negative_moments(p)
    side = isqrt(length(p))
    side^2 == length(p) || error("Nonsquare patch")
    coords = side == 1 ? [0.0] : collect(range(-1.0, 1.0; length=side))
    weights = max.(-p, 0.0)
    total = sum(weights)
    total <= eps(Float64) && return (com_t=NaN, diag45=NaN)
    t, u = ([t for _ in coords for t in coords], [u for u in coords for _ in coords])
    # Preserve the source's loop accumulation (rather than BLAS dot reduction).
    com_t, diag45 = 0.0, 0.0
    for i in eachindex(p)
        com_t += weights[i] * t[i]
        diag45 += weights[i] * ((t[i] + u[i]) / sqrt(2.0))
    end
    return (com_t=com_t / total, diag45=diag45 / total)
end

function load_inputs(paths, production_config)
    cfg = load_config(production_config)
    cfg["selection"]["gmm_training_weighting"] == "equal_lobes" ||
        error("This historical representation diagnostic requires equal-lobe GMM weighting")
    cfg["preprocessing"]["gmm_feature_normalization"] == "mean_sample_std" &&
        cfg["preprocessing"]["gmm_scale_fallback"] == 1.0 ||
        error("This historical representation diagnostic requires mean_sample_std with scale fallback 1")
    inputs = Dict{String,Any}()
    for name in ("patches-bwd9", "patches-bwd17", "patches-fwd17", "descriptor", "predictor", "fisher")
        header, rows = lobe_table(paths[name])
        any(c -> lowercase(c) in ("predicted", "label", "labels", "probability_1", "confidence"), header) &&
            error("Not a saved-patch/feature/score input: $name")
        inputs[name] = (header=header, rows=rows)
    end
    base = inputs["descriptor"].rows
    for (name, data) in inputs
        require_same_keys(base, data.rows, name)
    end
    for (name, columns) in (("descriptor", ["patch_u_asym_reconstructed", "descriptor_reason"]),
                            ("predictor", vcat(GMM_FEATURES, ["split_log_skew"])),
                            ("fisher", ["score", "invalid_reason"]))
        all(in(inputs[name].header), columns) || error("Missing required $name columns")
    end
    for (name, prefix, side) in (("patches-bwd9", "bwd_res", 9), ("patches-bwd17", "bwd_res", 17),
                                  ("patches-fwd17", "res", 17))
        columns = [@sprintf("%s_p%03d", prefix, i) for i in 1:side^2]
        count(startswith(prefix * "_p"), inputs[name].header) == side^2 || error("Wrong $name grid")
        all(in(inputs[name].header), columns) || error("Missing $name pixel columns")
    end
    fg = fisher_grid(load_fisher_config(production_config))
    fg.side == 17 || error("Expected saved 17x17 Fisher patches")
    return inputs, cfg, fg
end

function add_common!(diagnostics, patches, coords, keys; zero_l1)
    good = [k for k in keys if diagnostics[k]["reason"] == "ok"]
    cohort_gradient = isempty(good) ? NaN : mean(diagnostics[k]["gradient_u_per_nm"] for k in good)
    _, u = grid_vectors(coords)
    signed_u_mass = sum(abs, u)
    for file in unique(first.(keys))
        in_file = [k for k in good if k[1] == file]
        file_gradient = isempty(in_file) ? NaN : mean(diagnostics[k]["gradient_u_per_nm"] for k in in_file)
        for key in filter(k -> k[1] == file, keys)
            row = diagnostics[key]
            row["reason"] == "ok" || continue
            for (scope, gradient) in (("file", file_gradient), ("cohort", cohort_gradient))
                numerator = gradient * signed_u_mass
                row["common_" * scope * "_gradient_u_per_nm"] = gradient
                row["common_" * scope * "_u_over_raw_mass"] = numerator / row["mass_l1"]
                row["common_" * scope * "_u_share_of_numerator"] = ratio(numerator, row["numerator"])
                value, reason = transverse_asymmetry(patches[key] .- gradient .* u, coords; zero_l1)
                row["d_common_" * scope * "_gradient_removed_diagnostic"] = value
                row["common_" * scope * "_removed_reason"] = reason
            end
        end
    end
end

function summary_row(scope, family, metric, values)
    s = finite_summary(values)
    row = Dict{String,Any}("file" => scope, "family" => family, "metric" => metric,
                           "rows" => length(values), "finite" => s.n, "invalid" => length(values)-s.n)
    for k in (:mean, :std, :min, :q25, :median, :q75, :max)
        row[string(k)] = getfield(s, k)
    end
    row["positive"] = count(x -> isfinite(x) && x > 0, values)
    row["negative"] = count(x -> isfinite(x) && x < 0, values)
    row["zero"] = count(iszero, values)
    return row
end

function write_rows(path, header, rows)
    write_table(path, header, [Dict(k => text(v) for (k, v) in row) for row in rows])
end

function run_audit(paths; production_config, outdir, identity_atol, descriptor_atol)
    VERSION.major == 1 && VERSION.minor == 13 || error("This diagnostic requires Julia 1.13")
    Threads.nthreads() == 1 || error("Use --threads=1 for this saved-table diagnostic")
    for tolerance in (identity_atol, descriptor_atol)
        isfinite(tolerance) && tolerance >= 0 || error("Tolerances must be finite and nonnegative")
    end
    (ispath(outdir) || islink(outdir)) && error("Output already exists: $outdir")
    source_paths = merge(Dict(paths), Dict("production-config" => production_config))
    hashes = Dict(k => bytes2hex(sha256(read(v))) for (k, v) in source_paths)
    inputs, cfg, fg = load_inputs(paths, production_config)
    keys = sort(collect(Base.keys(inputs["descriptor"].rows)))
    zero_l1 = cfg["model"]["descriptor_zero_l1"]
    families = (("bwd9_res", "patches-bwd9", "bwd_res", 9, cfg["model"]["descriptor_step_nm"]),
                ("bwd17_res", "patches-bwd17", "bwd_res", 17, cfg["model"]["mold_step_nm"]),
                ("fwd17_res", "patches-fwd17", "res", 17, cfg["model"]["mold_step_nm"]))
    all_diagnostics = Dict{String,Any}()
    all_patches = Dict{String,Any}()
    patch_rows = Dict{String,Any}[]
    summaries = Dict{String,Any}[]
    descriptor_matches = 0
    for (family, source, prefix, side, step) in families
        coords = collect(-(side÷2):(side÷2)) .* step
        # This is the extractor's actual Julia StepRangeLen grid, not Fisher's
        # historical NumPy-arange disk support. Assert that it agrees exactly.
        half = side == 9 ? cfg["model"]["descriptor_half_nm"] : cfg["model"]["mold_half_nm"]
        coords == collect(-half:step:half) || error("Configured extractor grid mismatch")
        columns = [@sprintf("%s_p%03d", prefix, i) for i in 1:side^2]
        patches = Dict(k => [number(inputs[source].rows[k][c]) for c in columns] for k in keys)
        diagnostics = Dict(k => patch_audit(patches[k], coords; zero_l1, identity_atol) for k in keys)
        add_common!(diagnostics, patches, coords, keys; zero_l1)
        all_diagnostics[family], all_patches[family] = diagnostics, patches
        for key in keys
            row = diagnostics[key]
            row["file"], row["lobe"], row["family"] = key[1], key[2], family
            row["saved_descriptor"], row["saved_descriptor_reason"] = "NA", "not_production_descriptor_grid"
            row["descriptor_difference"] = NaN
            if side == 9
                saved = inputs["descriptor"].rows[key]
                d, reason = number(saved["patch_u_asym_reconstructed"]), saved["descriptor_reason"]
                reason == row["reason"] || error("Descriptor reason mismatch: $key")
                if reason == "ok"
                    isfinite(d) && abs(d-row["d_raw"]) <= descriptor_atol || error("Descriptor value mismatch: $key")
                    row["descriptor_difference"] = row["d_raw"] - d
                else
                    !isfinite(d) || error("Finite invalid saved descriptor: $key")
                end
                descriptor_matches += 1
                row["saved_descriptor"], row["saved_descriptor_reason"] = saved["patch_u_asym_reconstructed"], reason
            end
            push!(patch_rows, row)
        end
        for file in vcat(["ALL"], sort(unique(first.(keys))))
            selected = file == "ALL" ? keys : filter(k -> k[1] == file, keys)
            for metric in ("d_raw", "fraction_aligned", "fraction_opposed", "fraction_central",
                           "affine_constant", "gradient_t_per_nm", "gradient_u_per_nm",
                           "gradient_u_energy_fraction", "affine_residual_energy_fraction",
                           "u_over_raw_mass", "u_share_of_numerator", "d_gradient_removed_diagnostic",
                           "d_affine_removed_diagnostic", "common_file_u_share_of_numerator",
                           "d_common_file_gradient_removed_diagnostic", "common_cohort_u_share_of_numerator",
                           "d_common_cohort_gradient_removed_diagnostic")
                push!(summaries, summary_row(file, family, metric, [diagnostics[k][metric] for k in selected]))
            end
        end
    end

    # Exact finite-column preprocessing from the two existing predictors.
    files = first.(keys)
    raw = fill(NaN, length(keys), length(ALL_FEATURES))
    raw_text = fill("NA", size(raw))
    for (i, key) in enumerate(keys)
        moments = negative_moments(all_patches["bwd17_res"][key])
        for (j, feature) in enumerate(ALL_FEATURES)
            value = feature == "bwd_neg_com_t" ? text(moments.com_t) :
                    feature == "bwd_neg_diag45" ? text(moments.diag45) : inputs["predictor"].rows[key][feature]
            raw[i, j], raw_text[i, j] = number(value), value
        end
        saved = number(inputs["descriptor"].rows[key]["patch_u_asym_reconstructed"])
        joined = raw[i, findfirst(==("patch_u_asym_reconstructed"), ALL_FEATURES)]
        isequal(saved, joined) || error("Predictor descriptor differs from saved descriptor: $key")
        fisher = number(inputs["fisher"].rows[key]["score"])
        empirical = raw[i, findfirst(==("emp_fisher"), ALL_FEATURES)]
        expected_empirical = isfinite(fisher) ? parse(Float64, @sprintf("%.8g", -fisher)) : NaN
        isequal(expected_empirical, empirical) || error("Predictor Fisher sign/join mismatch: $key")
    end
    standardized = standardize_columns(raw, files; interactions=false)
    standard_rows, scaling_rows = Dict{String,Any}[], Dict{String,Any}[]
    for (i, key) in enumerate(keys), (j, feature) in enumerate(ALL_FEATURES)
        p = standardized.params[(key[1], j)]
        push!(standard_rows, Dict("file" => key[1], "lobe" => key[2], "feature" => feature,
            "input_value" => raw_text[i, j], "within_file_z" => standardized.z[i, j],
            "reason" => isfinite(raw[i, j]) ? "ok" : "nonfinite_feature",
            "source" => startswith(feature, "bwd_neg_") ? "existing_bwd17_negative_moment" : "saved_predictor"))
    end
    for file in sort(unique(files)), (j, feature) in enumerate(ALL_FEATURES)
        p = standardized.params[(file, j)]
        idx = findall(==(file), files)
        push!(scaling_rows, Dict("file" => file, "feature" => feature, "rows" => length(idx),
            "finite" => p.n, "file_mean_removed" => p.mean, "sample_std" => p.sample_std,
            "scale_used" => p.scale, "z_mean" => finite_summary(standardized.z[idx, j]).mean,
            "z_sample_std" => finite_summary(standardized.z[idx, j]).std))
    end
    views = [("gmm_v_cc", GMM_FEATURES), ("kmeans_v_base", BASE),
             ("kmeans_v_split", vcat(BASE, ["split_log_skew"])),
             ("kmeans_v_comt", vcat(BASE, ["bwd_neg_com_t"])),
             ("kmeans_v_diag45", vcat(BASE, ["bwd_neg_diag45"]))]
    view_rows, interaction_rows = Dict{String,Any}[], Dict{String,Any}[]
    for (view, features) in views
        indices = [findfirst(==(name), ALL_FEATURES) for name in features]
        trace = standardize_columns(raw[:, indices], files; interactions=cfg["selection"]["interactions"])
        for (i, key) in enumerate(keys)
            push!(view_rows, Dict("file" => key[1], "lobe" => key[2], "view" => view,
                "valid_preprocessing" => trace.valid[i], "dimensions" => size(trace.expanded, 2),
                "invalid_features" => join([features[j] for j in eachindex(features) if !isfinite(trace.z[i,j])], ",")))
        end
        if cfg["selection"]["interactions"]
            for a in 1:length(features)-1, b in a+1:length(features)
                vals = trace.z[:, a] .* trace.z[:, b]
                row = summary_row("ALL", view, features[a] * "*" * features[b], vals)
                row["operation"] = "product_of_within_file_z_no_second_standardization"
                push!(interaction_rows, row)
            end
        end
    end
    # Trace ONLY the existing descriptor's additive numerator attribution through
    # its existing scale. These terms sum to its z; no residual d is substituted.
    attribution_rows = Dict{String,Any}[]
    descriptor_column = findfirst(==("patch_u_asym_reconstructed"), ALL_FEATURES)
    for file in sort(unique(files))
        idx = findall(==(file), files)
        p = standardized.params[(file, descriptor_column)]
        ug = [all_diagnostics["bwd9_res"][keys[i]]["u_over_raw_mass"] for i in idx]
        remainder = [raw[i, descriptor_column] - ug[j] for (j, i) in enumerate(idx)]
        mean_u, mean_r = finite_summary(ug).mean, finite_summary(remainder).mean
        for (j, i) in enumerate(idx)
            uz, rz = (ug[j] - mean_u) / p.scale, (remainder[j] - mean_r) / p.scale
            err = uz + rz - standardized.z[i, descriptor_column]
            isfinite(err) && abs(err) > identity_atol && error("Standardized attribution identity failed")
            push!(attribution_rows, Dict("file" => file, "lobe" => keys[i][2], "descriptor_z" => standardized.z[i,descriptor_column],
                "u_attribution_raw_denominator" => ug[j], "remainder_raw_denominator" => remainder[j],
                "u_attribution_centered_over_descriptor_scale" => uz,
                "remainder_centered_over_descriptor_scale" => rz, "sum_error" => err))
        end
        uz = [number(text(r["u_attribution_centered_over_descriptor_scale"])) for r in attribution_rows if r["file"] == file]
        rz = [number(text(r["remainder_centered_over_descriptor_scale"])) for r in attribution_rows if r["file"] == file]
        for (metric, vals) in (("descriptor_u_additive_z_term", uz), ("descriptor_remainder_additive_z_term", rz))
            push!(summaries, summary_row(file, "existing_predictor_attribution", metric, vals))
        end
    end

    fisher_rows = Dict{String,Any}[]
    for key in keys
        row = inputs["fisher"].rows[key]
        push!(fisher_rows, Dict("file" => key[1], "lobe" => key[2], "saved_score" => row["score"],
            "saved_invalid_reason" => row["invalid_reason"], "heldout_lobe_parity" => iseven(key[2]) ? "even" : "odd",
            "training_lobe_parity" => iseven(key[2]) ? "odd" : "even",
            "decomposition_status" => "unavailable_saved_weights_and_mid_not_exported"))
    end
    for file in vcat(["ALL"], sort(unique(files))), parity in ("all", "even", "odd")
        selected = filter(k -> (file == "ALL" || k[1] == file) &&
            (parity == "all" || (iseven(k[2]) ? "even" : "odd") == parity), keys)
        push!(summaries, summary_row(file, "saved_fisher_heldout_" * parity, "saved_score",
                                    [number(inputs["fisher"].rows[k]["score"]) for k in selected]))
    end
    # Verify source immutability before creating an exclusive output directory.
    all(bytes2hex(sha256(read(source_paths[k]))) == h for (k,h) in hashes) || error("Input changed during diagnostic")
    mkpath(dirname(abspath(outdir)))
    mkdir(outdir)
    summary_header = ["file", "family", "metric", "rows", "finite", "invalid", "positive", "negative", "zero",
                      "mean", "std", "min", "q25", "median", "q75", "max"]
    write_rows(joinpath(outdir, "patch_diagnostics.tsv"), vcat(["file", "lobe", "family", "reason", "finite_pixels", "total_pixels",
        "saved_descriptor", "saved_descriptor_reason", "descriptor_difference"], PATCH_METRICS,
        ["gradient_removed_reason", "affine_removed_reason", "common_file_removed_reason", "common_cohort_removed_reason"]), patch_rows)
    write_rows(joinpath(outdir, "per_file_and_cohort.tsv"), summary_header, summaries)
    write_rows(joinpath(outdir, "standardization_rows.tsv"), ["file", "lobe", "feature", "input_value", "within_file_z", "reason", "source"], standard_rows)
    write_rows(joinpath(outdir, "standardization_scales.tsv"), ["file", "feature", "rows", "finite", "file_mean_removed", "sample_std", "scale_used", "z_mean", "z_sample_std"], scaling_rows)
    write_rows(joinpath(outdir, "predictor_view_support.tsv"), ["file", "lobe", "view", "valid_preprocessing", "dimensions", "invalid_features"], view_rows)
    write_rows(joinpath(outdir, "existing_interaction_spread.tsv"), vcat(summary_header, ["operation"]), interaction_rows)
    write_rows(joinpath(outdir, "descriptor_standardized_attribution.tsv"), ["file", "lobe", "descriptor_z", "u_attribution_raw_denominator",
        "remainder_raw_denominator", "u_attribution_centered_over_descriptor_scale", "remainder_centered_over_descriptor_scale", "sum_error"], attribution_rows)
    write_rows(joinpath(outdir, "saved_fisher_audit.tsv"), ["file", "lobe", "saved_score", "saved_invalid_reason",
        "heldout_lobe_parity", "training_lobe_parity", "decomposition_status"], fisher_rows)
    write_rows(joinpath(outdir, "inputs.tsv"), ["input", "path", "sha256"],
        [Dict("input" => k, "path" => abspath(source_paths[k]), "sha256" => hashes[k]) for k in sort(collect(Base.keys(hashes)))])
    write_report(joinpath(outdir, "report.md"), keys, all_diagnostics, summaries, fg,
                 descriptor_matches, identity_atol, descriptor_atol, zero_l1, cfg)
    println("Saved representation audit: $(length(unique(files))) files / $(length(keys)) keys; descriptor matches $descriptor_matches/$(length(keys))")
    println("Output: ", outdir, "; no classifier, refit, prediction, grading, or replacement feature table created")
    return (keys=length(keys), files=length(unique(files)), descriptor_matches=descriptor_matches)
end

function write_report(path, keys, diagnostics, summaries, grid, matches, identity_atol, descriptor_atol, zero_l1, cfg)
    med(family, metric) = text(finite_summary([diagnostics[family][k][metric] for k in keys]).median)
    short(x) = isfinite(x) ? @sprintf("%.6g", x) : "NA"
    function cell(family, metric)
        s = finite_summary([diagnostics[family][k][metric] for k in keys])
        "$(short(s.median)) [$(short(s.q25)), $(short(s.q75))]"
    end
    original_map = collect(eachindex(grid.mirror_indices))
    reflected_twice = [j == 0 || grid.mirror_indices[j] == 0 ? 0 : grid.mirror_indices[j] for j in grid.mirror_indices]
    open(path, "w") do io
        println(io, "# Saved-patch representation diagnostic\n")
        println(io, "Julia $(VERSION); $(length(unique(first.(keys)))) files / $(length(keys)) keys in all six inputs. All keys retained. No labels, predictions, classifier, raw STM fitting, refitted Fisher weights, or external grades were read or produced.\n")
        println(io, "## Descriptor and affine attribution\n")
        println(io, "The production 9x9 backward descriptor and reason match at **$matches/$(length(keys)) keys**. Arithmetic tolerance: $identity_atol; descriptor comparison: $descriptor_atol; unchanged production zero-L1 cutoff: $zero_l1. Invalid and zero-mass rows are NA, never imputed.\n")
        println(io, "On the actual symmetric serialized u-outer/t-inner grid, A is absolute mass aligned with sign(u), O is opposed mass, and C is central-row mass. M=A+O+C; d=(A-O)/M; 1-d=(2O+C)/M. Near +1 means little opposed/central mass. A transverse ramp alone can yield +1 without chemistry.\n")
        println(io, "Unweighted affine projection is p=c+b_t*t+b_u*u+r on each complete normalized residual patch. Constant and t terms cancel in the numerator on this grid. Numerator terms divided by ORIGINAL M are additive; their fractions need not lie in [0,1]. Energy fractions use centered squared energy. A large numerator share is not explained variance.\n")
        println(io, "| Family | Complete | Raw d median [Q1,Q3] | u energy fraction | u numerator share | d after own-u removal | d after file-common-u removal |\n|---|---:|---:|---:|---:|---:|---:|")
        for family in ("bwd9_res", "bwd17_res", "fwd17_res")
            n = count(k -> diagnostics[family][k]["reason"] == "ok", keys)
            println(io, "| $family | $n/$(length(keys)) | $(cell(family,"d_raw")) | $(cell(family,"gradient_u_energy_fraction")) | $(cell(family,"u_share_of_numerator")) | $(cell(family,"d_gradient_removed_diagnostic")) | $(cell(family,"d_common_file_gradient_removed_diagnostic")) |")
        end
        println(io, "\n| Family | File-common u numerator share median [Q1,Q3] | Cohort-common u numerator share | Cohort-common removed d |\n|---|---:|---:|---:|")
        for family in ("bwd9_res", "bwd17_res", "fwd17_res")
            println(io, "| $family | $(cell(family,"common_file_u_share_of_numerator")) | $(cell(family,"common_cohort_u_share_of_numerator")) | $(cell(family,"d_common_cohort_gradient_removed_diagnostic")) |")
        end
        println(io, "\n`raw` means unmodified saved normalized residual, not a raw-height image. Gradient-removed d recomputes its denominator WITHOUT renormalizing or masking pixels. It is a diagnostic only, not a proposed feature. File/cohort common gradients are simple means of complete-patch b_u, not cross-validated estimates; they describe shared structure, not predictive ability. Ratios can exceed one or be unstable near zero numerator; they are not variance explained. All-file distributions and invalid counts are in per_file_and_cohort.tsv.\n")
        println(io, "A strong common transverse component can account for a large baseline of this statistic. The affine projection cannot identify its cause: acquisition/preprocessing, alignment, background subtraction, molecular geometry and chemistry remain confounded in these saved normalized patches. Missing raw normalization scales prevent recovering physical gradients. Correlated lobes/files are not independent replicates. Positivity proves neither chemistry nor a software bug.\n")
        println(io, "## What reaches existing predictor standardization\n")
        println(io, "Both existing predictors subtract each feature's within-file mean over finite values and divide by its within-file sample std (scale=1 for zero/undefined std). A file-wide constant offset therefore disappears; small within-file variation survives and can be expanded to unit std. Each column uses its own finite rows BEFORE joint view validity. Pairwise products are formed afterward with no second normalization (enabled=$(cfg["selection"]["interactions"])). No clustering was run.\n")
        println(io, "standardization_rows/scales.tsv show every actual existing feature; predictor_view_support.tsv retains all keys in each existing view. Existing backward negative moments are recomputed only to trace those established k-means views. existing_interaction_spread.tsv summarizes those existing products, not new features. descriptor_standardized_attribution.tsv decomposes the EXISTING z(d) into its centered u numerator / original M term and the centered remainder using the SAME descriptor scale. Large cancelling terms are possible; gradient attribution is not causal attribution and does not prove predictor dominance.\n")
        for (family, metric, label) in (("bwd9_res", "d_raw", "Within-file raw descriptor sample std"),
                ("existing_predictor_attribution", "descriptor_u_additive_z_term", "Within-file sample std of additive u term in existing z(d)"),
                ("existing_predictor_attribution", "descriptor_remainder_additive_z_term", "Within-file sample std of additive remainder in existing z(d)"))
            values = [row["std"] for row in summaries if row["file"] != "ALL" && row["family"] == family && row["metric"] == metric]
            s = finite_summary(values)
            println(io, "- $label: median $(short(s.median)) [Q1 $(short(s.q25)), Q3 $(short(s.q75))], $(s.n) files with defined std.")
        end
        println(io, "\n## Saved Fisher scores: identified and missing\n")
        println(io, "The native producer exports scores only, not the two opposite-parity fold weights or mid vectors. They are absent from the specified saved inputs. Hence real even/odd response, shared centering offset and max-reflection uplift are **not identifiable** here. No weights are refit, no missing model is reconstructed, and no new real score is produced. saved_fisher_audit.tsv copies original score/reason text with held-out/training parity; per_file_and_cohort.tsv reports score distributions. Positive scores alone do not establish chemical classes.\n")
        println(io, "For SAVED weights w and midpoint m, the exact identity is max(w⋅(x-m),w⋅(R*x-m)) = w⋅(x+R*x)/2 - w⋅m + |w⋅(x-R*x)/2|. The shared centering offset is -w⋅m, not a mirrored m and not the original training mean. Synthetic tests compare this directly with native maxmirror_score. The historical R reverses physical t, not u, and zero-fills omitted disk support. This grid has $(length(grid.disk_indices)) disk pixels, $(count(iszero, grid.mirror_indices)) zero-padded destinations, and R² differs from identity on $(count(!iszero, reflected_twice .- original_map)) pixels. With incomplete mirror support, half-sum/difference need not be true parity eigenvectors; the response identity still holds exactly. No historical convention is corrected.\n")
        for row in summaries
            row["file"] == "ALL" && row["family"] == "saved_fisher_heldout_all" || continue
            println(io, "Saved Fisher: $(row["finite"])/$(row["rows"]) finite, $(row["positive"]) positive, $(row["negative"]) negative, $(row["zero"]) zero; median $(short(row["median"])), range $(short(row["min"])) to $(short(row["max"])).")
        end
        println(io, "\n## Limits and deliverables\n")
        println(io, "This is a representation audit, not chemical validation, a count check, a feature campaign, a promotion decision or a replacement for raw acquisition controls. The three saved grids are diagnostic support checks, not a resolution sweep selected by outcome. inputs.tsv records exact input hashes; the originals were not changed. No production defaults or outputs are overwritten.")
    end
end

end # module

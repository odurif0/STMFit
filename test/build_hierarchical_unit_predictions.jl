#!/usr/bin/env julia

# Build label-free 0/1/? unit-assignment predictions using the composition-free
# hierarchical equal-prior model (T4 of the unit-assignment challenger plan).
#
# This CLI is the label-free counterpart to `build_labelfree_unit_predictions.jl`.
# It never reads benchmark truth, expected sequences, control motifs, lobe
# index/order, transitions, run lengths, per-chain class counts, top-k rules,
# occupancy priors, or composition counts. Cluster labels are mapped by the
# physical convention already used by the grader: the higher-amplitude cluster
# is GlcNAc (1). Use grade_unit_assignment.jl or report_unit_assignment_benchmark.jl
# only after this TSV has been written.
#
# Output columns: file, lobe, predicted, confidence, amplitude, probability_1,
# views_used, invalid_reason (the frozen existing schema), plus model,
# model_version, provenance_sha256 (the new challenger fields).

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
using .HierarchicalUnitAssignment
using Printf

const DEFAULT_FEATURES = "results/unit_separability/lobe_features_selectedN_primary_local.tsv"
const DEFAULT_SPLIT = ""
const DEFAULT_PATCHES = ""
const DEFAULT_OUT = "results/unit_assignment/hierarchical_unit_predictions.tsv"

struct Options
    features::String
    split_features::String
    patches::String
    out_tsv::String
    view_specs::Vector{Pair{String,Vector{String}}}
    first_seed::Int
    n_starts::Int
    cov_floor::Float64
    max_iter::Int
    tol::Float64
end

function _arg_value(args, i::Int, flag::String)
    i < length(args) || error("$flag requires a value")
    return args[i+1]
end

function _flag_name(arg::AbstractString)
    return split(String(arg), '='; limit=2)[1]
end

function _parse_view(s::AbstractString)
    parts = split(String(s), '='; limit=2)
    length(parts) == 2 || error("--view must be NAME=feature1,feature2,...")
    name = strip(parts[1])
    isempty(name) && error("--view name is empty")
    features = [strip(f) for f in split(parts[2], ',') if !isempty(strip(f))]
    isempty(features) && error("--view $name has no features")
    for fn in features
        is_forbidden_feature(fn) &&
            error("--view $name uses a forbidden feature (lobe-order/benchmark/truth leak): $fn")
    end
    return name => features
end

function _parse_cli(args)
    # Reject any forbidden benchmark/truth/lobe-order/sequence/transition flag
    # up front so it cannot influence the run.
    for arg in args
        if is_forbidden_flag(arg)
            error("Forbidden benchmark/label/lobe-order flag (label-free " *
                  "hierarchical model): $(_flag_name(arg))")
        end
    end

    features = DEFAULT_FEATURES
    split_features = DEFAULT_SPLIT
    patches = DEFAULT_PATCHES
    out_tsv = DEFAULT_OUT
    view_specs = Pair{String,Vector{String}}[]
    first_seed = 0
    n_starts = default_n_starts
    cov_floor = cov_floor_default
    max_iter = 200
    tol = 1e-6

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"
            features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--features=")
            features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--split-features"
            split_features = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--split-features=")
            split_features = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--patches"
            patches = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--patches=")
            patches = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = _arg_value(args, i, arg); i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "="; limit=2)[2]; i += 1
        elseif arg == "--view"
            push!(view_specs, _parse_view(_arg_value(args, i, arg))); i += 2
        elseif startswith(arg, "--view=")
            push!(view_specs, _parse_view(split(arg, "=", limit=2)[2])); i += 1
        elseif arg == "--first-seed"
            first_seed = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--first-seed=")
            first_seed = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--em-starts"
            n_starts = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--em-starts=")
            n_starts = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--cov-floor"
            cov_floor = parse(Float64, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--cov-floor=")
            cov_floor = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--em-max-iter"
            max_iter = parse(Int, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--em-max-iter=")
            max_iter = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--em-tol"
            tol = parse(Float64, _arg_value(args, i, arg)); i += 2
        elseif startswith(arg, "--em-tol=")
            tol = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            _show_help()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    n_starts > 0 || error("--em-starts must be positive")
    cov_floor > 0 || error("--cov-floor must be positive")
    max_iter > 0 || error("--em-max-iter must be positive")
    tol > 0 || error("--em-tol must be positive")
    isfile(features) || error("Feature TSV not found: $features")
    !isempty(split_features) && !isfile(split_features) &&
        error("Split feature TSV not found: $split_features")
    !isempty(patches) && !isfile(patches) &&
        error("Patch TSV not found: $patches")
    return Options(features, split_features, patches, out_tsv, view_specs,
                   first_seed, n_starts, cov_floor, max_iter, tol)
end

function _show_help()
    println("""
            Usage: julia --project=. test/build_hierarchical_unit_predictions.jl [options]

            Composition-free hierarchical equal-prior unit-assignment predictor.

            Label-free contract: this script reads no truth sequence, expected
            N, control motif, lobe index/order, transition structure, run
            length, per-chain class count, top-k rule, occupancy prior, or
            composition count. Equal class priors are fixed at exactly (0.5,
            0.5) — the equal prior policy is immutable and never updated. The
            higher-amplitude cluster is mapped to GlcNAc (1) by the physical
            convention. Benchmark grading is a separate post-hoc step.

            Options:
              --features PATH        Main per-lobe feature TSV [$(DEFAULT_FEATURES)]
              --split-features PATH  Optional split-width feature TSV; adds split_log_skew
              --patches PATH         Optional backward patch TSV; adds bwd_neg_* descriptors
              --out PATH             Output prediction TSV [$(DEFAULT_OUT)]
              --view NAME=LIST       Feature view to fit. Repeatable. If omitted,
                                      sensible label-free defaults are chosen from
                                      available columns.
              --first-seed INT       First EM multi-start index [0]
              --em-starts INT        Number of deterministic EM multi-starts [$(default_n_starts)]
              --cov-floor FLOAT      Diagonal variance floor [$(cov_floor_default)]
              --em-max-iter INT      EM max iterations [200]
              --em-tol FLOAT         EM convergence tolerance [1e-6]

            Output columns: file, lobe, predicted, confidence, amplitude,
            probability_1, views_used, invalid_reason, model, model_version,
            provenance_sha256.

            Forbidden (label-free firewall): --truth, --control-sequence,
            --manifest, --benchmark-manifest, --expected-N, --target-N,
            --transitions, --run-length, --top-k, --class-count, --occupancy,
            --lobe-index, --lobe-order, --lobe-position, --sequence, --full145,
            --control, and any feature name encoding sequence/position/order/
            transition/run-length/truth/benchmark information.
            """)
end

# --------------------------------------------------------------------------
# Optional merge helpers (split-width and backward patches) - mirror the
# existing labelfree builder but operate through the LobeRecord interface.
# --------------------------------------------------------------------------

function _parse_split_value(s)
    t = strip(String(s))
    (isempty(t) || t in ("NA", "NaN", "nan", "?")) && return NaN
    v = tryparse(Float64, t)
    return v === nothing ? NaN : v
end

function _merge_split_features!(records::Vector{LobeRecord}, path::AbstractString)
    isempty(path) && return nothing
    header, rows = HierarchicalUnitAssignment._read_tsv(path)
    "skew_ratio" in header || error("Split TSV $path missing skew_ratio column")
    by_key = Dict{Tuple{String,Int},Float64}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Split TSV $path has empty file value")
        haskey(row, "lobe") || error("Split TSV $path missing lobe column")
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe === nothing && error("Split TSV $path invalid lobe: $(row["lobe"])")
        lobe >= 1 || error("Split TSV $path lobe must be >= 1, got $lobe")
        key = (file, lobe)
        haskey(by_key, key) && error("Split TSV $path duplicate row for $file lobe $lobe")
        skew = _parse_split_value(row["skew_ratio"])
        by_key[key] = isfinite(skew) && skew > 0 ? log(skew) : NaN
    end
    for rec in records
        rec.features["split_log_skew"] = get(by_key, (rec.file, rec.lobe), NaN)
    end
    return nothing
end

function _patch_columns(header::Vector{String}, prefix::String)
    cols = [c for c in header if startswith(c, prefix)]
    isempty(cols) && return String[]
    side = round(Int, sqrt(length(cols)))
    side * side == length(cols) ||
        error("Patch columns for prefix $prefix are not a square grid: $(length(cols))")
    return cols
end

function _negative_moment(vals::Vector{Float64}, coords::Vector{Float64})
    side = length(coords)
    weights = max.(-vals, 0.0)
    total = sum(weights)
    total <= eps(Float64) && return (com_t=NaN, diag45=NaN, diag135=NaN)
    com_t = 0.0
    diag45 = 0.0
    diag135 = 0.0
    idx = 1
    for u in coords, t in coords
        w = weights[idx]
        com_t += w * t
        diag45 += w * ((t + u) / sqrt(2.0))
        diag135 += w * ((t - u) / sqrt(2.0))
        idx += 1
    end
    return (com_t=com_t / total, diag45=diag45 / total, diag135=diag135 / total)
end

function _merge_patches!(records::Vector{LobeRecord}, path::AbstractString)
    isempty(path) && return nothing
    header, rows = HierarchicalUnitAssignment._read_tsv(path)
    cols = _patch_columns(String.(header), "bwd_res_p")
    isempty(cols) && error("Patch TSV $path has no bwd_res_pNNN columns")
    side = round(Int, sqrt(length(cols)))
    coords = side == 1 ? [0.0] : collect(range(-1.0, 1.0; length=side))
    by_key = Dict{Tuple{String,Int},NamedTuple{(:com_t,:diag45,:diag135),Tuple{Float64,Float64,Float64}}}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Patch TSV $path has empty file value")
        haskey(row, "lobe") || error("Patch TSV $path missing lobe column")
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe === nothing && error("Patch TSV $path invalid lobe: $(row["lobe"])")
        lobe >= 1 || error("Patch TSV $path lobe must be >= 1, got $lobe")
        key = (file, lobe)
        haskey(by_key, key) && error("Patch TSV $path duplicate row for $file lobe $lobe")
        vals = Float64[_parse_split_value(row[c]) for c in cols]
        by_key[key] = _negative_moment(vals, coords)
    end
    for rec in records
        d = get(by_key, (rec.file, rec.lobe),
                (com_t=NaN, diag45=NaN, diag135=NaN))
        rec.features["bwd_neg_com_t"] = d.com_t
        rec.features["bwd_neg_diag45"] = d.diag45
        rec.features["bwd_neg_diag135"] = d.diag135
    end
    return nothing
end

# --------------------------------------------------------------------------
# Default views chosen from available columns. These match the production
# label-free feature lists frozen in config/unit_assignment_candidate.toml
# (features.base_local, base_gaussian, backward_descriptors, split_descriptor).
# --------------------------------------------------------------------------

function _available_features(records::Vector{LobeRecord})
    names = Set{String}()
    for rec in records
        union!(names, keys(rec.features))
    end
    return names
end

function _default_views(records::Vector{LobeRecord})
    available = _available_features(records)
    base_local = ["amp_prominence", "amp_neighbor_ratio",
                  "integrated_prominence", "amp_rel"]
    base_gaussian = ["amplitude", "sigma_parallel_nm",
                     "sigma_perp_nm", "integrated"]
    base = all(in(available), base_local) ? base_local :
           all(in(available), base_gaussian) ? base_gaussian :
           error("No default feature base is available; pass --view explicitly")

    views = Pair{String,Vector{String}}[]
    for (name, extra) in (("base_local", nothing),
                          ("base_local+bwd_neg_com_t", "bwd_neg_com_t"),
                          ("base_local+bwd_neg_diag45", "bwd_neg_diag45"),
                          ("base_local+split_log_skew", "split_log_skew"))
        if extra === nothing
            push!(views, name => copy(base))
        elseif extra in available
            push!(views, name => vcat(base, [extra]))
        end
    end
    isempty(views) && error("No usable default views; pass --view explicitly")
    return views
end

function main(args=ARGS)
    opt = _parse_cli(args)
    records = load_records(opt.features)
    _merge_split_features!(records, opt.split_features)
    _merge_patches!(records, opt.patches)
    provenance_inputs = HierarchicalUnitAssignment.ProvenanceInput[
        HierarchicalUnitAssignment.ProvenanceInput("primary_features", opt.features),
    ]
    if !isempty(opt.split_features)
        push!(provenance_inputs,
              HierarchicalUnitAssignment.ProvenanceInput("split_features", opt.split_features))
    end
    if !isempty(opt.patches)
        push!(provenance_inputs,
              HierarchicalUnitAssignment.ProvenanceInput("backward_patches", opt.patches))
    end
    views = isempty(opt.view_specs) ? _default_views(records) : opt.view_specs
    # Re-check firewall for every resolved view feature.
    for (name, feats) in views
        for fn in feats
            is_forbidden_feature(fn) &&
                error("View $name uses a forbidden feature: $fn")
            any(haskey(r.features, fn) for r in records) ||
                error("View $name uses an unavailable feature: $fn")
        end
    end
    run_pipeline(records, provenance_inputs, opt.out_tsv, views;
                 first_seed=opt.first_seed, n_starts=opt.n_starts,
                 cov_floor=opt.cov_floor, max_iter=opt.max_iter, tol=opt.tol)
    return nothing
end

main()

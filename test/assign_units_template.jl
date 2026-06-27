#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# assign_units_template.jl — Supervised empirical-template validation for
# GlcNAc/GlcN unit assignment.
#
# This script is intentionally SUPERVISED: it reads the ground-truth unit
# sequence, splits files into train/test, learns one feature centroid per unit
# type on train, and predicts test lobes by nearest centroid. It is used only to
# measure whether the STM features contain enough information to generalize;
# it must not be used in the fitting/model-selection path.
#
# Truth orientation must match the lobe order (t_nm increasing) unless
# --truth-orientation reverse is supplied. If orientation is uncertain, use this
# as a controlled validation tool after fixing a convention.
# ──────────────────────────────────────────────────────────────────────────────

using Printf
using Statistics
using Random
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_FEATURES = "results/unit_separability/lobe_features_selectedN.tsv"
const DEFAULT_TRUTH = "benchmarks/chitosan_240817_unit_sequences.tsv"
const DEFAULT_OUT = "results/unit_assignment/template_predictions.tsv"

struct Options
    features::String
    extra::Union{Nothing,String}
    truth::String
    out_tsv::String
    feature_names::Vector{String}
    train_fraction::Float64
    seed::Int
    orientation::String
end

function _parse_cli(args)
    features = DEFAULT_FEATURES
    extra::Union{Nothing,String} = nothing
    truth = DEFAULT_TRUTH
    out_tsv = DEFAULT_OUT
    feature_names = ["amplitude", "sigma_parallel_nm", "sigma_perp_nm", "integrated"]
    train_fraction = 0.7
    seed = 42
    orientation = "identity"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"; features = args[i+1]; i += 2
        elseif startswith(arg, "--features="); features = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--extra"; extra = args[i+1]; i += 2
        elseif startswith(arg, "--extra="); extra = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--truth"; truth = args[i+1]; i += 2
        elseif startswith(arg, "--truth="); truth = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--features-list"; feature_names = split(args[i+1], ","); i += 2
        elseif startswith(arg, "--features-list="); feature_names = split(split(arg, "=", limit=2)[2], ","); i += 1
        elseif arg == "--train-fraction"; train_fraction = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--train-fraction="); train_fraction = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--seed"; seed = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--seed="); seed = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--truth-orientation"; orientation = args[i+1]; i += 2
        elseif startswith(arg, "--truth-orientation="); orientation = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/assign_units_template.jl [options]

            Options:
              --features PATH          Lobe features TSV [$(DEFAULT_FEATURES)]
              --extra PATH             Optional residual features TSV
              --truth PATH             Truth sequence TSV [$(DEFAULT_TRUTH)]
              --out PATH               Output predictions TSV [$(DEFAULT_OUT)]
              --features-list LIST     Comma-separated feature names
              --train-fraction FLOAT   Fraction of files used for train [0.7]
              --seed INT               Split RNG seed [42]
              --truth-orientation STR  identity | reverse [identity]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return Options(features, extra, truth, out_tsv, String.(strip.(feature_names)), train_fraction, seed, orientation)
end

_basename_file(s) = basename(strip(s))
_parse_seq(s) = [parse(Int, c) for c in strip(s) if c in "01"]

function _load_truth(path, orientation)
    _, rows = _read_tsv(path)
    truth = Dict{String,Vector{Int}}()
    for row in rows
        seq = _parse_seq(get(row, "sequence", ""))
        isempty(seq) && continue
        orientation == "reverse" && (seq = reverse(seq))
        truth[_basename_file(row["file"])] = seq
    end
    return truth
end

function _load_extra(path)
    by = Dict{Tuple{String,Int},Dict{String,Float64}}()
    path === nothing && return by, String[]
    header, rows = _read_tsv(path)
    extra_cols = [c for c in header if !(c in ["file", "lobe", "amplitude", "t_nm", "u_nm", "sigma_parallel_nm", "sigma_perp_nm"])]
    for row in rows
        file = _basename_file(row["file"])
        lobe = parse(Int, row["lobe"])
        d = Dict{String,Float64}()
        for c in extra_cols
            d[c] = _parse_f(get(row, c, "NA"))
        end
        by[(file, lobe)] = d
    end
    return by, extra_cols
end

function main()
    opt = _parse_cli(ARGS)
    truth = _load_truth(opt.truth, opt.orientation)
    extra_by, extra_cols = _load_extra(opt.extra)
    _, rows = _read_tsv(opt.features)
    isempty(rows) && error("No feature rows: $(opt.features)")

    records = []
    for row in rows
        file = _basename_file(row["file"])
        lobe = parse(Int, row["lobe"])
        haskey(truth, file) || continue
        lobe <= length(truth[file]) || continue
        amp = _parse_f(row["amplitude"])
        spar = _parse_f(row["sigma_parallel_nm"])
        sperp = _parse_f(row["sigma_perp_nm"])
        feat = Dict{String,Float64}(
            "amplitude" => amp,
            "sigma_parallel_nm" => spar,
            "sigma_perp_nm" => sperp,
            "integrated" => amp * spar * sperp,
        )
        for (k, v) in row
            k in ("file", "N", "lobe", "source") && continue
            haskey(feat, k) && continue
            parsed = tryparse(Float64, strip(v))
            parsed === nothing || (feat[k] = parsed)
        end
        haskey(extra_by, (file, lobe)) && merge!(feat, extra_by[(file, lobe)])
        push!(records, (file=file, lobe=lobe, y=truth[file][lobe], feat=feat))
    end
    isempty(records) && error("No records with truth labels. Fill $(opt.truth) first.")

    feature_names = copy(opt.feature_names)
    for c in extra_cols
        c in feature_names || push!(feature_names, c)
    end

    files = sort(unique(r.file for r in records))
    rng = MersenneTwister(opt.seed)
    shuffled = shuffle(rng, files)
    n_train = max(1, min(length(files) - 1, round(Int, opt.train_fraction * length(files))))
    train_files = Set(shuffled[1:n_train])
    test_files = Set(shuffled[(n_train+1):end])

    # z-score by train statistics only
    train_records = [r for r in records if r.file in train_files]
    stats = Dict{String,Tuple{Float64,Float64}}()
    for f in feature_names
        vals = filter(isfinite, [get(r.feat, f, NaN) for r in train_records])
        isempty(vals) && error("Feature $f has no finite train values")
        σ = std(vals)
        stats[f] = (mean(vals), σ > 0 ? σ : 1.0)
    end
    function vecfeat(r)
        x = Float64[]
        for f in feature_names
            v = get(r.feat, f, NaN)
            isfinite(v) || return nothing
            μ, σ = stats[f]
            push!(x, (v - μ) / σ)
        end
        return x
    end

    centroids = Dict{Int,Vector{Float64}}()
    for y in (0, 1)
        xs = [vecfeat(r) for r in train_records if r.y == y]
        xs = [x for x in xs if x !== nothing]
        isempty(xs) && error("No train examples for class $y")
        centroids[y] = vec(mean(reduce(hcat, xs); dims=2))
    end

    mkpath(dirname(opt.out_tsv))
    correct = total = 0
    open(opt.out_tsv, "w") do io
        println(io, join(["file", "lobe", "truth", "predicted", "correct", "set"], '\t'))
        for r in records
            x = vecfeat(r)
            x === nothing && continue
            d0 = norm(x .- centroids[0])
            d1 = norm(x .- centroids[1])
            pred = d1 <= d0 ? 1 : 0
            setname = r.file in train_files ? "train" : "test"
            if setname == "test"
                total += 1
                correct += pred == r.y ? 1 : 0
            end
            println(io, join([r.file, r.lobe, r.y, pred, pred == r.y, setname], '\t'))
        end
    end

    println("Template assignment")
    println("  features: ", join(feature_names, ", "))
    println("  train files: ", length(train_files), "  test files: ", length(test_files))
    println("  test accuracy: ", correct, "/", total, " (", total > 0 ? @sprintf("%.1f%%", 100correct/total) : "NA", ")")
    println("  output: ", opt.out_tsv)
end

main()

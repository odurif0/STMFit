#!/usr/bin/env julia

# PCA/k-means separability diagnostics for aligned lobe patches. Label-free
# clustering is performed without truth. If --truth is provided, labels are used
# only for post-hoc external evaluation or supervised train/test diagnostics.

using Printf
using Statistics
using LinearAlgebra
using Random
using Clustering

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv

const DEFAULT_PATCHES = "results/unit_separability/lobe_patches_selectedN_primary.tsv"
const DEFAULT_OUT = "results/unit_separability/patch_analysis"

struct Options
    patches::String
    out_dir::String
    truth::Union{Nothing,String}
    prefix::String
    n_pcs::Int
    train_fraction::Float64
    seed::Int
end

function _parse_cli(args)
    patches = DEFAULT_PATCHES
    out_dir = DEFAULT_OUT
    truth::Union{Nothing,String} = nothing
    prefix = "res_p"
    n_pcs = 8
    train_fraction = 0.7
    seed = 42
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--patches"; patches = args[i+1]; i += 2
        elseif startswith(arg, "--patches="); patches = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_dir = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_dir = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--truth"; truth = args[i+1]; i += 2
        elseif startswith(arg, "--truth="); truth = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--prefix"; prefix = args[i+1]; i += 2
        elseif startswith(arg, "--prefix="); prefix = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--pcs"; n_pcs = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--pcs="); n_pcs = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--train-fraction"; train_fraction = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--train-fraction="); train_fraction = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--seed"; seed = parse(Int, args[i+1]); i += 2
        elseif startswith(arg, "--seed="); seed = parse(Int, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/analyze_lobe_patches.jl [options]

            Options:
              --patches PATH         Patch TSV from extract_lobe_patches.jl [$(DEFAULT_PATCHES)]
              --out PATH             Output directory [$(DEFAULT_OUT)]
              --truth PATH           Optional truth TSV, external evaluation only
              --prefix STR           Patch columns to use: raw_p or res_p [res_p]
              --pcs INT              Number of principal components [8]
              --train-fraction FLOAT Supervised diagnostic train file fraction [0.7]
              --seed INT             Split seed [42]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(patches) || error("Patch TSV not found: $patches")
    truth !== nothing && !isfile(truth) && error("Truth TSV not found: $truth")
    return Options(patches, out_dir, truth, prefix, n_pcs, train_fraction, seed)
end

_parse_seq(s) = [parse(Int, c) for c in strip(s) if c in "01"]

function _kmeans_bic(data::Matrix{Float64}, k::Int)
    nf, n = size(data)
    if k == 1
        center = mean(data, dims=2)
        sse = sum((data .- center).^2)
        params = nf
    else
        r = kmeans(data, k; maxiter=200, rng=MersenneTwister(42), display=:none)
        sse = r.totalcost
        params = k * nf + (k - 1)
    end
    return n * log(sse / n + 1e-30) + params * log(n), sse
end

function _load_truth(path)
    truth = Dict{String,Vector{Int}}()
    path === nothing && return truth
    _, rows = _read_tsv(path)
    for row in rows
        seq = _parse_seq(get(row, "sequence", ""))
        isempty(seq) || (truth[basename(row["file"])] = seq)
    end
    return truth
end

function _nearest_centroid(train_X, train_y, test_X)
    centroids = Dict{Int,Vector{Float64}}()
    for y in (0, 1)
        idx = findall(==(y), train_y)
        isempty(idx) && error("No train examples for class $y")
        centroids[y] = vec(mean(train_X[idx, :], dims=1))
    end
    pred = Int[]
    for i in 1:size(test_X, 1)
        x = vec(test_X[i, :])
        d0 = norm(x .- centroids[0])
        d1 = norm(x .- centroids[1])
        push!(pred, d1 <= d0 ? 1 : 0)
    end
    return pred
end

function main()
    opt = _parse_cli(ARGS)
    header, rows = _read_tsv(opt.patches)
    pix_cols = [c for c in header if startswith(c, opt.prefix)]
    isempty(pix_cols) && error("No columns with prefix $(opt.prefix)")
    n = length(rows); p = length(pix_cols)
    X = zeros(n, p)
    valid = trues(n)
    files = String[]; lobes = Int[]
    for (i, row) in enumerate(rows)
        push!(files, basename(row["file"]))
        push!(lobes, parse(Int, row["lobe"]))
        for (j, c) in enumerate(pix_cols)
            v = _parse_f(row[c])
            if !isfinite(v)
                valid[i] = false
                X[i, j] = 0.0
            else
                X[i, j] = v
            end
        end
    end
    X = X[valid, :]
    files = files[valid]
    lobes = lobes[valid]
    μ = vec(mean(X, dims=1))
    Xc = X .- μ'
    sv = svd(Xc; full=false)
    k = min(opt.n_pcs, size(sv.V, 2))
    scores = Xc * sv.V[:, 1:k]
    explained = sv.S.^2 ./ sum(sv.S.^2)

    bic1, sse1 = _kmeans_bic(permutedims(scores), 1)
    bic2, sse2 = _kmeans_bic(permutedims(scores), 2)
    km = kmeans(permutedims(scores), 2; maxiter=200, rng=MersenneTwister(42), display=:none)

    truth = _load_truth(opt.truth)
    y = Int[]; labeled_idx = Int[]
    for i in eachindex(files)
        seq = get(truth, files[i], Int[])
        if lobes[i] <= length(seq)
            push!(y, seq[lobes[i]])
            push!(labeled_idx, i)
        end
    end

    mkpath(opt.out_dir)
    open(joinpath(opt.out_dir, "patch_scores.tsv"), "w") do io
        println(io, join(vcat(["file", "lobe", "cluster"], ["pc$i" for i in 1:k]), '\t'))
        for i in eachindex(files)
            println(io, join(vcat([files[i], string(lobes[i]), string(km.assignments[i])],
                                  [@sprintf("%.7g", scores[i, j]) for j in 1:k]), '\t'))
        end
    end

    report = joinpath(opt.out_dir, "patch_report.txt")
    open(report, "w") do io
        println(io, "Patch Separability Report")
        println(io, "Patches TSV: ", opt.patches)
        println(io, "Prefix: ", opt.prefix)
        println(io, "Rows valid: ", size(X, 1), "  pixels: ", p, "  PCs: ", k)
        println(io, "Explained variance first PCs: ", join([@sprintf("%.3f", explained[i]) for i in 1:min(k, length(explained))], ", "))
        println(io)
        println(io, "Label-free k=1 vs k=2 on PC scores")
        println(io, "k=1 BIC: ", @sprintf("%.1f", bic1), "  SSE=", @sprintf("%.4f", sse1))
        println(io, "k=2 BIC: ", @sprintf("%.1f", bic2), "  SSE=", @sprintf("%.4f", sse2))
        println(io, "ΔBIC (k=1 - k=2): ", @sprintf("%.1f", bic1 - bic2))
        println(io, "Cluster sizes: [", count(==(1), km.assignments), ", ", count(==(2), km.assignments), "]")
        if opt.truth !== nothing && !isempty(y)
            rng = MersenneTwister(opt.seed)
            labeled_files = sort(unique(files[labeled_idx]))
            shuffled = shuffle(rng, labeled_files)
            n_train = max(1, min(length(shuffled) - 1, round(Int, opt.train_fraction * length(shuffled))))
            train_files = Set(shuffled[1:n_train])
            train = [i for i in labeled_idx if files[i] in train_files]
            test = [i for i in labeled_idx if !(files[i] in train_files)]
            y_by_i = Dict(labeled_idx[j] => y[j] for j in eachindex(labeled_idx))
            pred = _nearest_centroid(scores[train, :], [y_by_i[i] for i in train], scores[test, :])
            truth_test = [y_by_i[i] for i in test]
            correct = count(pred[i] == truth_test[i] for i in eachindex(pred))
            println(io)
            println(io, "External supervised diagnostic")
            println(io, "Train files: ", length(train_files), "  test files: ", length(setdiff(Set(labeled_files), train_files)))
            println(io, "Test accuracy: ", correct, "/", length(test), " (", @sprintf("%.1f%%", 100 * correct / max(1, length(test))), ")")
        end
    end
    println("Wrote report: ", report)
    println("ΔBIC (k=1 - k=2): ", @sprintf("%.1f", bic1 - bic2))
end

main()

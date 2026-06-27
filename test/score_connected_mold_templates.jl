#!/usr/bin/env julia

# Connected mold-template scoring for unit assignment.
#
# This script applies physically oriented GlcN/GlcNAc patch templates to the
# aligned lobe patches produced by extract_lobe_patches.jl. It tests a small set
# of global connectivity states (chain direction, pyranose parity phase, mirror)
# and decodes the whole chain with Viterbi. Optional bond templates add sliding
# pairwise costs for transitions 00/01/10/11 on every adjacent pair. No truth
# sequence and no composition count is used: the number of GlcNAc units emerges
# from template costs.
#
# Template TSV format (wide):
#   name  type  parity  mirror  p001  p002  ...
# where type is 0=GlcN or 1=GlcNAc, parity is 0/1, mirror is 0/1, and pNNN
# columns must match the chosen patch prefix/length after stripping the prefix.
# Optional bond TSV format:
#   name  left_type  right_type  parity  mirror  l_p001 ... r_p001 ...

using Printf
using Statistics
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv, _standardize

const DEFAULT_PATCHES = "results/unit_separability/lobe_patches_selectedN_primary.tsv"
const DEFAULT_TEMPLATES = "templates/chitosan_connected_molds.tsv"
const DEFAULT_OUT = "results/unit_assignment/connected_mold_predictions.tsv"

struct Options
    patches::String
    templates::String
    bond_templates::Union{Nothing,String}
    out_tsv::String
    prefix::String
    score::String
    template_mode::String
    transition_penalty::Float64
    bond_weight::Float64
end

function _parse_cli(args)
    patches = DEFAULT_PATCHES
    templates = DEFAULT_TEMPLATES
    bond_templates::Union{Nothing,String} = nothing
    out_tsv = DEFAULT_OUT
    prefix = "res_p"
    score = "ncc"
    template_mode = "full"
    transition_penalty = 0.0
    bond_weight = 1.0
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--patches"; patches = args[i+1]; i += 2
        elseif startswith(arg, "--patches="); patches = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--templates"; templates = args[i+1]; i += 2
        elseif startswith(arg, "--templates="); templates = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--bond-templates"; bond_templates = args[i+1]; i += 2
        elseif startswith(arg, "--bond-templates="); bond_templates = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"; out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out="); out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--prefix"; prefix = args[i+1]; i += 2
        elseif startswith(arg, "--prefix="); prefix = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--score"; score = args[i+1]; i += 2
        elseif startswith(arg, "--score="); score = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--template-mode"; template_mode = args[i+1]; i += 2
        elseif startswith(arg, "--template-mode="); template_mode = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--transition-penalty"; transition_penalty = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--transition-penalty="); transition_penalty = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--bond-weight"; bond_weight = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--bond-weight="); bond_weight = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/score_connected_mold_templates.jl [options]

            Options:
              --patches PATH       Patch TSV from extract_lobe_patches.jl [$(DEFAULT_PATCHES)]
              --templates PATH     Wide template TSV [$(DEFAULT_TEMPLATES)]
              --bond-templates PATH
                                   Optional sliding bond template TSV for
                                   transitions 00/01/10/11
              --out PATH           Output per-lobe predictions TSV [$(DEFAULT_OUT)]
              --prefix STR         Patch columns to score: raw_p or res_p [res_p]
              --score STR          ncc | sse [ncc]
              --template-mode STR  full | contrast. contrast subtracts the
                                   parity/mirror common mold before scoring [full]
              --transition-penalty FLOAT
                                   Optional pairwise type-switch penalty for
                                   Viterbi decoding within a global orientation [0]
              --bond-weight FLOAT  Weight for bond-template transition costs [1]

            Template TSV columns:
              name, type, parity, mirror, p001, p002, ...

            Optional bond TSV columns:
              name, left_type, right_type, parity, mirror, l_p001, ..., r_p001, ...

            Label-free constraints:
              Truth labels and composition counts are not read. The script tests
              only global connectivity states (direction, parity phase, mirror)
              and sliding adjacent-pair costs; it never tiles by dimers, so odd
              and even chains are both supported.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(patches) || error("Patch TSV not found: $patches")
    isfile(templates) || error("Template TSV not found: $templates")
    bond_templates !== nothing && !isfile(bond_templates) && error("Bond template TSV not found: $bond_templates")
    score = lowercase(strip(score))
    score in ("ncc", "sse") || error("--score must be ncc or sse")
    template_mode = lowercase(strip(template_mode))
    template_mode in ("full", "contrast") || error("--template-mode must be full or contrast")
    return Options(patches, templates, bond_templates, out_tsv, prefix, score, template_mode, transition_penalty, bond_weight)
end

function _score_patch(patch::Vector{Float64}, templ::Vector{Float64}, method::String)
    n_good = 0
    dot_pt = 0.0
    norm_p2 = 0.0
    norm_t2 = 0.0
    sse = 0.0
    for i in eachindex(patch)
        p = patch[i]
        t = templ[i]
        if isfinite(p) && isfinite(t)
            n_good += 1
            dot_pt += p * t
            norm_p2 += p * p
            norm_t2 += t * t
            d = p - t
            sse += d * d
        end
    end
    if n_good < max(5, cld(length(patch), 2))
        return Inf
    end
    if method == "sse"
        return sse / n_good
    end
    denom = sqrt(norm_p2 * norm_t2)
    denom <= 0 && return Inf
    # Cost: lower is better. Negative normalized correlation rewards alignment.
    return -dot_pt / denom
end

function _load_templates(path::String)
    header, rows = _read_tsv(path)
    pix_cols = [c for c in header if occursin(r"^p\d+$", c)]
    isempty(pix_cols) && error("Template TSV has no pNNN columns: $path")
    by_key = Dict{Tuple{Int,Int,Int},Vector{Float64}}()
    for row in rows
        typ = parse(Int, row["type"])
        parity = parse(Int, row["parity"])
        mirror = parse(Int, row["mirror"])
        typ in (0, 1) || error("Template type must be 0 or 1")
        parity in (0, 1) || error("Template parity must be 0 or 1")
        mirror in (0, 1) || error("Template mirror must be 0 or 1")
        vec = _standardize([_parse_f(row[c]) for c in pix_cols])
        by_key[(typ, parity, mirror)] = vec
    end
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        haskey(by_key, (typ, parity, mirror)) || error("Missing template type=$typ parity=$parity mirror=$mirror")
    end
    return by_key, pix_cols
end

function _load_bond_templates(path::Union{Nothing,String})
    path === nothing && return nothing, String[]
    header, rows = _read_tsv(path)
    left_cols = [c for c in header if occursin(r"^l_p\d+$", c)]
    right_cols = [c for c in header if occursin(r"^r_p\d+$", c)]
    isempty(left_cols) && error("Bond template TSV has no l_pNNN columns: $path")
    length(left_cols) == length(right_cols) || error("Bond template left/right pixel counts differ")
    by_key = Dict{Tuple{Int,Int,Int,Int},Vector{Float64}}()
    for row in rows
        lt = parse(Int, row["left_type"])
        rt = parse(Int, row["right_type"])
        parity = parse(Int, row["parity"])
        mirror = parse(Int, row["mirror"])
        lt in (0, 1) || error("left_type must be 0 or 1")
        rt in (0, 1) || error("right_type must be 0 or 1")
        parity in (0, 1) || error("parity must be 0 or 1")
        mirror in (0, 1) || error("mirror must be 0 or 1")
        vec = _standardize(vcat([_parse_f(row[c]) for c in left_cols],
                                [_parse_f(row[c]) for c in right_cols]))
        by_key[(lt, rt, parity, mirror)] = vec
    end
    for lt in (0, 1), rt in (0, 1), parity in (0, 1), mirror in (0, 1)
        haskey(by_key, (lt, rt, parity, mirror)) || error("Missing bond template left=$lt right=$rt parity=$parity mirror=$mirror")
    end
    return by_key, left_cols
end

function _contrast_unary_templates!(templates)
    for parity in (0, 1), mirror in (0, 1)
        k0 = (0, parity, mirror)
        k1 = (1, parity, mirror)
        common = 0.5 .* (templates[k0] .+ templates[k1])
        templates[k0] = _standardize(templates[k0] .- common)
        templates[k1] = _standardize(templates[k1] .- common)
    end
    return templates
end

function _contrast_bond_templates!(templates)
    templates === nothing && return nothing
    for parity in (0, 1), mirror in (0, 1)
        keys = [(lt, rt, parity, mirror) for lt in (0, 1) for rt in (0, 1)]
        common = zeros(length(templates[first(keys)]))
        for key in keys
            common .+= templates[key]
        end
        common ./= length(keys)
        for key in keys
            templates[key] = _standardize(templates[key] .- common)
        end
    end
    return templates
end

function _load_patches(path::String, prefix::String)
    header, rows = _read_tsv(path)
    pix_cols = [c for c in header if startswith(c, prefix)]
    isempty(pix_cols) && error("Patch TSV has no columns with prefix $prefix")
    by_file = Dict{String,Vector{NamedTuple}}()
    for row in rows
        file = basename(row["file"])
        lobe = parse(Int, row["lobe"])
        patch = _standardize([_parse_f(row[c]) for c in pix_cols])
        amplitude = haskey(row, "amplitude") ? _parse_f(row["amplitude"]) : NaN
        push!(get!(by_file, file, NamedTuple[]), (file=file, lobe=lobe, amplitude=amplitude, patch=patch))
    end
    for file in keys(by_file)
        sort!(by_file[file], by=r -> r.lobe)
    end
    return by_file, pix_cols
end

function _parity_for_lobe(lobe::Int, n::Int, direction::Int, phase::Int)
    idx = direction == 0 ? lobe : (n - lobe + 1)
    return mod(idx - 1 + phase, 2)
end

function _edge_template_key(prev_label::Int, label::Int, i::Int, records, n::Int, direction::Int, phase::Int, mirror::Int)
    # The DP iterates in observed lobe order. For reversed molecular direction,
    # the chemical left/right order of the edge is current -> previous.
    if direction == 0
        parity = _parity_for_lobe(records[i-1].lobe, n, direction, phase)
        return (prev_label, label, parity, mirror)
    else
        parity = _parity_for_lobe(records[i].lobe, n, direction, phase)
        return (label, prev_label, parity, mirror)
    end
end

function _bond_observation(prev_patch::Vector{Float64}, patch::Vector{Float64}, direction::Int)
    return direction == 0 ? vcat(prev_patch, patch) : vcat(patch, prev_patch)
end

function _viterbi(costs::Matrix{Float64}, records, direction::Int, phase::Int, mirror::Int,
                  bond_templates, score_method::String, switch_penalty::Float64, bond_weight::Float64)
    n = size(costs, 1)
    dp = fill(Inf, n, 2)
    prev = zeros(Int, n, 2)
    dp[1, :] .= costs[1, :]
    for i in 2:n, y in 1:2
        best_val = Inf
        best_prev = 1
        for yp in 1:2
            prev_label = yp - 1
            label = y - 1
            edge_cost = yp == y ? 0.0 : switch_penalty
            if bond_templates !== nothing
                key = _edge_template_key(prev_label, label, i, records, n, direction, phase, mirror)
                obs = _bond_observation(records[i-1].patch, records[i].patch, direction)
                edge_cost += bond_weight * _score_patch(obs, bond_templates[key], score_method)
            end
            val = dp[i-1, yp] + costs[i, y] + edge_cost
            if val < best_val
                best_val = val
                best_prev = yp
            end
        end
        dp[i, y] = best_val
        prev[i, y] = best_prev
    end
    labels = zeros(Int, n)
    labels[n] = dp[n, 2] < dp[n, 1] ? 1 : 0
    state = labels[n] + 1
    for i in (n-1):-1:1
        state = prev[i+1, state]
        labels[i] = state - 1
    end
    return labels, minimum(dp[n, :])
end

function _decode_file(records, templates, bond_templates, opt::Options)
    n = length(records)
    best = nothing
    for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
        costs = zeros(n, 2)
        for (i, rec) in enumerate(records)
            parity = _parity_for_lobe(rec.lobe, n, direction, phase)
            for typ in (0, 1)
                costs[i, typ+1] = _score_patch(rec.patch, templates[(typ, parity, mirror)], opt.score)
            end
        end
        labels, total = _viterbi(costs, records, direction, phase, mirror, bond_templates,
                                 opt.score, opt.transition_penalty, opt.bond_weight)
        margin = mean(abs.(costs[:, 1] .- costs[:, 2]))
        if best === nothing || total < best.total
            best = (total=total, labels=labels, costs=costs, direction=direction,
                    phase=phase, mirror=mirror, margin=margin)
        end
    end
    return best
end

function main()
    opt = _parse_cli(ARGS)
    templates, template_pix = _load_templates(opt.templates)
    bond_templates, bond_pix = _load_bond_templates(opt.bond_templates)
    if opt.template_mode == "contrast"
        _contrast_unary_templates!(templates)
        _contrast_bond_templates!(bond_templates)
    end
    by_file, patch_pix = _load_patches(opt.patches, opt.prefix)
    length(template_pix) == length(patch_pix) || error("Template pixel count $(length(template_pix)) != patch pixel count $(length(patch_pix))")
    bond_templates !== nothing && length(bond_pix) != length(patch_pix) && error("Bond template pixel count $(length(bond_pix)) != patch pixel count $(length(patch_pix))")
    mkpath(dirname(opt.out_tsv))

    open(opt.out_tsv, "w") do io
        println(io, join(["file", "lobe", "predicted", "amplitude", "physical_label", "cost_GlcN", "cost_GlcNAc",
                          "cost_margin", "global_direction", "global_phase", "global_mirror", "file_cost"], '\t'))
        for file in sort(collect(keys(by_file)))
            records = by_file[file]
            best = _decode_file(records, templates, bond_templates, opt)
            for (i, rec) in enumerate(records)
                c0, c1 = best.costs[i, 1], best.costs[i, 2]
                label = best.labels[i]
                println(io, join([file, rec.lobe, label, @sprintf("%.8g", rec.amplitude), label == 1 ? "GlcNAc" : "GlcN",
                                  @sprintf("%.8g", c0), @sprintf("%.8g", c1),
                                  @sprintf("%.8g", abs(c0 - c1)), best.direction,
                                  best.phase, best.mirror, @sprintf("%.8g", best.total)], '\t'))
            end
        end
    end

    println("Connected mold decoding")
    println("  patches:   ", opt.patches)
    println("  templates: ", opt.templates)
    opt.bond_templates !== nothing && println("  bonds:     ", opt.bond_templates, "  weight=", opt.bond_weight)
    println("  prefix:    ", opt.prefix)
    println("  score:     ", opt.score)
    println("  template:  ", opt.template_mode)
    println("  files:     ", length(by_file))
    println("  output:    ", opt.out_tsv)

    println("\nPredicted sequences:")
    _, pred_rows = _read_tsv(opt.out_tsv)
    by_pred = Dict{String,Vector{Tuple{Int,Int}}}()
    for row in pred_rows
        push!(get!(by_pred, row["file"], Tuple{Int,Int}[]), (parse(Int, row["lobe"]), parse(Int, row["predicted"])))
    end
    for file in sort(collect(keys(by_pred)))
        entries = sort(by_pred[file], by=x -> x[1])
        seq = join(string(e[2]) for e in entries)
        println("  ", file, "  seq=", seq, "  GlcNAc=", count(e[2] == 1 for e in entries))
    end
end

main()

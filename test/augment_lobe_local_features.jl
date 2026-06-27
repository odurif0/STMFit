#!/usr/bin/env julia

# Add label-free local/envelope-corrected per-lobe features to an existing
# extract_lobe_features.jl TSV. This does not use truth labels, target motifs,
# or any composition constraint.

using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv, _standardize

const DEFAULT_IN = "results/unit_separability/lobe_features_selectedN_primary.tsv"
const DEFAULT_OUT = "results/unit_separability/lobe_features_selectedN_primary_local.tsv"

function _parse_cli(args)
    input = DEFAULT_IN
    output = DEFAULT_OUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--features"
            input = args[i+1]; i += 2
        elseif startswith(arg, "--features=")
            input = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            output = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            output = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/augment_lobe_local_features.jl [options]

            Options:
              --features PATH   Input lobe features TSV [$(DEFAULT_IN)]
              --out PATH        Output augmented TSV [$(DEFAULT_OUT)]

            Adds only label-free chain-internal features. No truth sequence,
            expected motif, or number-of-units constraint is read or assumed.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(input) || error("Input TSV not found: $input")
    return input, output
end

_fmt(x) = isfinite(x) ? @sprintf("%.8g", x) : "NA"

function _poly_residuals(pos::Vector{Float64}, vals::Vector{Float64}, degree::Int)
    n = length(vals)
    n == 0 && return Float64[]
    valid = isfinite.(pos) .& isfinite.(vals)
    count(valid) < degree + 1 && return fill(NaN, n)
    X = zeros(count(valid), degree + 1)
    p = pos[valid]
    for d in 0:degree
        X[:, d+1] .= p .^ d
    end
    β = X \ vals[valid]
    residuals = fill(NaN, n)
    for i in 1:n
        if isfinite(pos[i]) && isfinite(vals[i])
            fit = sum(β[d+1] * pos[i]^d for d in 0:degree)
            residuals[i] = vals[i] - fit
        end
    end
    return residuals
end

function main()
    input, output = _parse_cli(ARGS)
    header, rows = _read_tsv(input)
    isempty(rows) && error("No rows in $input")

    by_file = Dict{String,Vector{Dict{String,String}}}()
    for row in rows
        push!(get!(by_file, basename(row["file"]), Dict{String,String}[]), row)
    end

    new_cols = [
        "amp_file_z", "amp_rel_file_z", "integrated", "integrated_file_z",
        "amp_prominence", "amp_neighbor_ratio", "amp_second_diff",
        "amp_env_resid_linear", "amp_env_resid_quadratic",
        "integrated_prominence", "sigma_parallel_file_z", "sigma_perp_file_z",
        "centered_pos", "edge_distance_norm"
    ]

    mkpath(dirname(output))
    open(output, "w") do io
        println(io, join(vcat(header, new_cols), '\t'))
        for file in sort(collect(keys(by_file)))
            rs = sort(by_file[file], by=r -> parse(Int, r["lobe"]))
            n = length(rs)
            amp = [_parse_f(r["amplitude"]) for r in rs]
            amp_rel = [_parse_f(r["amp_rel"]) for r in rs]
            spar = [_parse_f(r["sigma_parallel_nm"]) for r in rs]
            sperp = [_parse_f(r["sigma_perp_nm"]) for r in rs]
            integ = amp .* spar .* sperp
            pos = n == 1 ? [0.0] : [((i - 1) / (n - 1)) * 2 - 1 for i in 1:n]
            edge = [min(i - 1, n - i) / max(1, floor(Int, n / 2)) for i in 1:n]
            amp_z = _standardize(amp)
            amp_rel_z = _standardize(amp_rel)
            integ_z = _standardize(integ)
            spar_z = _standardize(spar)
            sperp_z = _standardize(sperp)
            lin_res = _poly_residuals(pos, amp_rel, 1)
            quad_res = _poly_residuals(pos, amp_rel, min(2, n - 1))

            for i in 1:n
                left = i > 1 ? amp_rel[i-1] : amp_rel[i]
                right = i < n ? amp_rel[i+1] : amp_rel[i]
                neigh = (left + right) / 2
                amp_prom = amp_rel[i] - neigh
                amp_ratio = amp_rel[i] / max(neigh, eps(Float64))
                amp_second = amp_rel[i] - neigh

                ileft = i > 1 ? integ[i-1] : integ[i]
                iright = i < n ? integ[i+1] : integ[i]
                integ_prom = integ[i] - (ileft + iright) / 2

                vals = [rs[i][h] for h in header]
                extras = [_fmt(x) for x in (
                    amp_z[i], amp_rel_z[i], integ[i], integ_z[i],
                    amp_prom, amp_ratio, amp_second,
                    lin_res[i], quad_res[i], integ_prom,
                    spar_z[i], sperp_z[i], pos[i], edge[i]
                )]
                println(io, join(vcat(vals, extras), '\t'))
            end
        end
    end
    println("Wrote: ", output)
end

main()

#!/usr/bin/env julia

using Printf
using Statistics

const DEFAULT_INPUT = joinpath(
    @__DIR__, "..", "templates",
    "chitosan_connected_molds_stm_dft_m030_h050_v1_half032.tsv",
)
const DEFAULT_OUTPUT = joinpath(
    @__DIR__, "..", "results", "joint_proxy_mold_visualization",
    "stm_dft_v1_pairwise_separability.tsv",
)

function load_molds(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("empty mold table: $path")
    header = split(lines[1], '\t')
    pixel_columns = findall(name -> occursin(r"^p\d{3}$", name), header)
    length(pixel_columns) == 81 || error("expected 81 pixels, found $(length(pixel_columns))")

    molds = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        push!(molds, (
            name=fields[1],
            type=parse(Int, fields[2]),
            parity=parse(Int, fields[3]),
            mirror=parse(Int, fields[4]),
            pixels=parse.(Float64, fields[pixel_columns]),
        ))
    end
    length(molds) == 8 || error("expected 8 molds, found $(length(molds))")
    return molds
end

rmse(a, b) = sqrt(mean((a .- b) .^ 2))

function pairwise_rows(molds)
    rows = NamedTuple[]
    for i in 1:(length(molds) - 1), j in (i + 1):length(molds)
        a, b = molds[i], molds[j]
        same_type = a.type == b.type
        same_state = a.parity == b.parity && a.mirror == b.mirror
        push!(rows, (
            mold_a=a.name,
            mold_b=b.name,
            relation=same_type ? "within_type" : "between_type",
            same_state=same_state,
            rmse=rmse(a.pixels, b.pixels),
            correlation=cor(a.pixels, b.pixels),
            max_abs_difference=maximum(abs, a.pixels .- b.pixels),
        ))
    end
    return rows
end

function write_report(rows, output_path::AbstractString)
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "mold_a\tmold_b\trelation\tsame_state\trmse\tcorrelation\tmax_abs_difference")
        for row in rows
            @printf(
                io, "%s\t%s\t%s\t%s\t%.9g\t%.9g\t%.9g\n",
                row.mold_a, row.mold_b, row.relation, row.same_state,
                row.rmse, row.correlation, row.max_abs_difference,
            )
        end
    end
end

function print_summary(molds, rows)
    within = filter(row -> row.relation == "within_type", rows)
    between = filter(row -> row.relation == "between_type", rows)
    matched = filter(row -> row.relation == "between_type" && row.same_state, rows)
    nearest_correct = count(molds) do mold
        others = filter(other -> other.name != mold.name, molds)
        nearest = argmin(other -> rmse(mold.pixels, other.pixels), others)
        nearest.type == mold.type
    end

    @printf("within-type RMSE range:  %.6f .. %.6f\n", extrema(getfield.(within, :rmse))...)
    @printf("between-type RMSE range: %.6f .. %.6f\n", extrema(getfield.(between, :rmse))...)
    @printf("matched-state RMSE mean: %.6f\n", mean(getfield.(matched, :rmse)))
    @printf("matched-state corr mean: %.6f\n", mean(getfield.(matched, :correlation)))
    println("nearest-template same-type: $nearest_correct/$(length(molds))")
end

function main(args)
    length(args) <= 2 || error("usage: julia --project=. test/analyze_joint_proxy_molds.jl [input.tsv] [output.tsv]")
    input_path = abspath(get(args, 1, DEFAULT_INPUT))
    output_path = abspath(get(args, 2, DEFAULT_OUTPUT))
    molds = load_molds(input_path)
    rows = pairwise_rows(molds)
    write_report(rows, output_path)
    print_summary(molds, rows)
    println("wrote $output_path")
end

main(ARGS)

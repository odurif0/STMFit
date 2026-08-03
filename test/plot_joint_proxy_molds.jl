#!/usr/bin/env julia

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Plots
using Printf

const DEFAULT_INPUT = joinpath(
    @__DIR__, "..", "templates",
    "chitosan_connected_molds_stm_dft_m030_h050_v1_half032.tsv",
)
const DEFAULT_OUTPUT = joinpath(
    @__DIR__, "..", "results", "joint_proxy_mold_visualization",
    "stm_dft_v1_glcn_glcnac_difference.png",
)

function load_molds(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("empty mold table: $path")
    header = split(lines[1], '\t')
    pixel_columns = findall(name -> occursin(r"^p\d{3}$", name), header)
    length(pixel_columns) == 81 || error("expected 81 pixels, found $(length(pixel_columns))")

    molds = Dict{Tuple{Int,Int,Int},Matrix{Float64}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        typ = parse(Int, fields[2])
        parity = parse(Int, fields[3])
        mirror = parse(Int, fields[4])
        values = parse.(Float64, fields[pixel_columns])
        molds[(typ, parity, mirror)] = reshape(values, 9, 9)'
    end
    length(molds) == 8 || error("expected 8 molds, found $(length(molds))")
    return molds
end

function render_molds(molds, output_path::AbstractString)
    coords = collect(range(-0.32, 0.32; length=9))
    raw_values = reduce(vcat, [vec(mold) for mold in values(molds)])
    raw_limits = extrema(raw_values)
    differences = [molds[(1, p, m)] - molds[(0, p, m)] for p in 0:1 for m in 0:1]
    difference_limit = maximum(abs, reduce(vcat, vec.(differences)))
    panels = Plots.Plot[]

    for parity in 0:1, mirror in 0:1
        glcn = molds[(0, parity, mirror)]
        glcnac = molds[(1, parity, mirror)]
        difference = glcnac - glcn
        row_name = @sprintf("parity=%d, mirror=%d", parity, mirror)

        push!(panels, heatmap(
            coords, coords, glcn;
            title="GlcN — $row_name", color=:viridis, clims=raw_limits,
            aspect_ratio=:equal, xlabel="t (nm)", ylabel="u (nm)",
        ))
        push!(panels, heatmap(
            coords, coords, glcnac;
            title="GlcNAc — $row_name", color=:viridis, clims=raw_limits,
            aspect_ratio=:equal, xlabel="t (nm)", ylabel="u (nm)",
        ))
        push!(panels, heatmap(
            coords, coords, difference;
            title="GlcNAc − GlcN — $row_name", color=:RdBu,
            clims=(-difference_limit, difference_limit),
            aspect_ratio=:equal, xlabel="t (nm)", ylabel="u (nm)",
        ))
    end

    mkpath(dirname(output_path))
    figure = plot(
        panels...;
        layout=(4, 3), size=(1500, 1700), dpi=180,
        plot_title="DFT-STM unary molds at −0.300 V and 0.50 nm",
        margin=4 * Plots.mm,
    )
    savefig(figure, output_path)
    return output_path
end

function main(args)
    length(args) <= 2 || error("usage: julia --project=. test/plot_joint_proxy_molds.jl [input.tsv] [output.png]")
    input_path = abspath(get(args, 1, DEFAULT_INPUT))
    output_path = abspath(get(args, 2, DEFAULT_OUTPUT))
    molds = load_molds(input_path)
    println("wrote ", render_molds(molds, output_path))
end

main(ARGS)

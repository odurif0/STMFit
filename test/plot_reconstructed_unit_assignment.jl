#!/usr/bin/env julia
# Native Julia 0/1/? maps. No benchmark or reference predictions are read.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using Plots

# At most eight full-size panels per image, in the original sorted file order.
overview_pages(n::Integer) = [i:min(i + 7, n) for i in 1:8:n]

function assignment_overview(panels, page, npages)
    1 <= length(panels) <= 8 || error("An overview page needs 1–8 panels")
    ncols = min(2, length(panels))
    nrows = cld(length(panels), ncols)
    height = 440nrows + 48
    return plot(panels...; layout=(nrows, ncols), size=(780ncols, height),
                plot_title="Reconstructed 0/1/? assignments — page $page/$npages",
                plot_titlefontsize=12, plot_titlevspan=48 / height,
                left_margin=10Plots.mm, widen=1.12)
end

function render_maps(features, predictions, outdir)
    _, geometry = lobe_table(features; required=["x_nm", "y_nm"])
    _, assigned = lobe_table(predictions; required=["predicted", "confidence"])
    require_same_keys(geometry, assigned, "Plot predictions")
    (ispath(outdir) || islink(outdir)) && error("Plot directory already exists: $outdir")
    mkpath(joinpath(outdir, "standalone"))
    files = sort(unique(first.(collect(keys(geometry)))))
    colors = Dict("0" => "#2166ac", "1" => "#b2182b", "?" => "#8c8c8c")
    labels = Dict("0" => "GlcN (0)", "1" => "GlcNAc (1)", "?" => "uncertain (?)")
    panels = Any[]
    for file in files
        keys_ = sort([key for key in keys(geometry) if first(key) == file]; by=last)
        xs = [parse(Float64, geometry[k]["x_nm"]) for k in keys_]
        ys = [parse(Float64, geometry[k]["y_nm"]) for k in keys_]
        all(isfinite, xs) && all(isfinite, ys) || error("Nonfinite plot coordinates: $file")
        pred = [assigned[k]["predicted"] for k in keys_]
        all(p -> haskey(colors, p), pred) || error("Invalid prediction in plot: $file")
        fig = plot(xs, ys; color=:gray, label="", aspect_ratio=:equal, yflip=true,
                   xlabel="x (nm)", ylabel="y (nm)", title="$file  N=$(length(keys_))",
                   legend=:outertopright, size=(780, 440), margin=5Plots.mm)
        for p in ("0", "1", "?")
            idx = findall(==(p), pred)
            scatter!(fig, xs[idx], ys[idx]; color=colors[p], label=labels[p],
                     markersize=9, markerstrokewidth=0)
        end
        for (i, p) in enumerate(pred)
            annotate!(fig, xs[i], ys[i], text(p, 9, p == "?" ? :black : :white))
        end
        filename = splitext(file)[1] * "_chain.png"
        savefig(fig, joinpath(outdir, "standalone", filename))
        push!(panels, fig)
    end
    # Keep the standalone panel footprint. A single very tall GR canvas squeezes
    # equal-aspect plot areas and a fractional title band grows with the cohort.
    # summary_grid.png is page 1; its title and the index make pagination explicit.
    pages = overview_pages(length(panels))
    open(joinpath(outdir, "summary_pages.tsv"), "w") do io
        println(io, "page\timage\tfile")
        for (page, indices) in enumerate(pages)
            grid = assignment_overview(panels[indices], page, length(pages))
            image = page == 1 ? "summary_grid.png" : "summary_grid_$(lpad(page, 3, '0')).png"
            savefig(grid, joinpath(outdir, image))
            for i in indices
                println(io, page, '\t', image, '\t', files[i])
            end
        end
    end
    return length(files)
end

function main(args=ARGS)
    if "--help" in args || "-h" in args
        println("Usage: julia --project=. test/plot_reconstructed_unit_assignment.jl --features PATH --predictions PATH --outdir DIR")
        return
    end
    iseven(length(args)) || error("Expected --option VALUE pairs")
    opts = Dict{String,String}()
    for i in 1:2:length(args)
        args[i] in ("--features", "--predictions", "--outdir") || error("Unknown or forbidden option: $(args[i])")
        haskey(opts, args[i]) && error("Repeated option: $(args[i])")
        opts[args[i]] = args[i+1]
    end
    all(k -> haskey(opts, k), ("--features", "--predictions", "--outdir")) || error("Missing required options")
    println("Mapped chains: ", render_maps(opts["--features"], opts["--predictions"], opts["--outdir"]))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

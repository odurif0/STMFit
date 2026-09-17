#!/usr/bin/env julia
# Focused plot-only regression: no fits, classifier, raw data or benchmark inputs.
ENV["GKSwstype"] = "100"
using Test
include(joinpath(@__DIR__, "plot_reconstructed_unit_assignment.jl"))

function plot_fixture(dir, n)
    files = ["synthetic_$(lpad(i, 3, '0')).sxm" for i in 1:n]
    features, predictions = joinpath(dir, "geometry.tsv"), joinpath(dir, "predictions.tsv")
    open(features, "w") do io
        println(io, "file\tlobe\tx_nm\ty_nm")
        for (i, file) in enumerate(files), j in 1:3
            # Include horizontal, vertical and oblique coordinate maps.
            x = i % 3 == 0 ? 1.0 : Float64(j)
            y = i % 3 == 1 ? 1.0 : Float64(j)
            println(io, "$file\t$j\t$x\t$y")
        end
    end
    open(predictions, "w") do io
        println(io, "file\tlobe\tpredicted\tconfidence")
        for file in files, (j, label) in enumerate(("0", "1", "?"))
            println(io, "$file\t$j\t$label\t0.0")
        end
    end
    return files, features, predictions
end

function check_panel_content(sp)
    @test sp[:aspect_ratio] == :equal
    @test sp[:yaxis][:flip]
    @test sp[:xaxis][:guide] == "x (nm)"
    @test sp[:yaxis][:guide] == "y (nm)"
    @test [a[3].str for a in sp[:annotations]] == ["0", "1", "?"]
    @test [s[:label] for s in sp.series_list[2:4]] == ["GlcN (0)", "GlcNAc (1)", "uncertain (?)"]
    @test [s[:markercolor] for s in sp.series_list[2:4]] ==
          [Plots.plot_color(c) for c in ("#2166ac", "#b2182b", "#8c8c8c")]
    @test all(s[:markersize] == 9 for s in sp.series_list[2:4])
    # These are rendered plot-area dimensions, not only the outer PNG size.
    box = Plots.plotarea(sp)
    width, height = Plots.width(box).value, Plots.height(box).value
    @test width > 90                 # mm: labels cannot consume the body
    @test height > 50
    xl, yl = xlims(sp), ylims(sp)
    @test width / (xl[2] - xl[1]) ≈ height / (yl[2] - yl[1])
    # Every colored endpoint has room for its disk, including horizontal maps.
    for s in sp.series_list[2:4], (x, y) in zip(s[:x], s[:y])
        xmargin = min(x - xl[1], xl[2] - x) / (xl[2] - xl[1]) * width
        ymargin = min(y - yl[1], yl[2] - y) / (yl[2] - yl[1]) * height
        @test min(xmargin, ymargin) * 100 / 25.4 >= 9
    end
end

@testset "Reconstructed assignment overview layout" begin
    @test length(overview_pages(146)) == 19
    @test reduce(vcat, overview_pages(146)) == collect(1:146)
    @test length.(overview_pages(146)) == vcat(fill(8, 18), 2)
    @test_throws ErrorException assignment_overview(Any[], 1, 1)
    @test_throws ErrorException assignment_overview(fill(nothing, 9), 1, 1)
    mktempdir() do tmp
        for n in (1, 3, 8, 146)
            dir = joinpath(tmp, "cohort_$n")
            mkpath(dir)
            files, features, predictions = plot_fixture(dir, n)
            before = (read(features), read(predictions))
            out = joinpath(dir, "plots")
            @test render_maps(features, predictions, out) == n
            @test (read(features), read(predictions)) == before
            @test sort(readdir(joinpath(out, "standalone"))) ==
                  [splitext(f)[1] * "_chain.png" for f in files]
            index = [split(line, '\t') for line in readlines(joinpath(out, "summary_pages.tsv"))[2:end]]
            @test [row[3] for row in index] == files
            @test [parse(Int, row[1]) for row in index] == [cld(i, 8) for i in 1:n]
            images = unique(row[2] for row in index)
            @test first(images) == "summary_grid.png"
            @test length(images) == cld(n, 8)
            @test all(filesize(joinpath(out, image)) > 5_000 for image in images)
            @test Set(filter(f -> endswith(f, ".png"), readdir(out))) == Set(images)
            grid = current()  # The actual final page created by render_maps.
            @test occursin("page $(cld(n, 8))/$(cld(n, 8))", grid[:plot_title])
            @test grid[:size][1] <= 1560 && grid[:size][2] <= 1808
            @test grid[:plot_titlevspan] * grid[:size][2] ≈ 48
            last_files = files[last(overview_pages(n))]
            for (sp, file) in zip(grid.subplots, last_files)
                @test sp[:title] == "$file  N=3"
                check_panel_content(sp)
            end
            @test_throws ErrorException render_maps(features, predictions, out)
        end
    end
end

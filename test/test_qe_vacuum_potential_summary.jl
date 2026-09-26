#!/usr/bin/env julia
# Lightweight saved-table checks; no cube loading or scientific recomputation.
using Test
include(joinpath(@__DIR__, "summarize_qe_vacuum_potential.jl"))
const M = SummarizeQEVacuumPotential
length(ARGS) == 2 || error("test_qe_vacuum_potential_summary.jl SAVED_RUN SAVED_SUMMARY_TSV")
run, saved = ARGS
rows = M.summaries(run)
@testset "Saved small-table summaries" begin
    @test length(rows) == 4
    for row in rows
        p = filter(r -> parse(Bool, r[row.domain]),
            M.R.V.tsv(joinpath(run, row.molecule, "analysis/planes.tsv")))
        b = filter(r -> parse(Bool, r[row.domain]),
            M.R.V.tsv(joinpath(run, row.molecule, "analysis/barriers.tsv")))
        @test length(p) == row.planes
        @test length(b) == row.band_planes
        @test sum(parse(Bool, r["all_sampled_barriers_positive"]) for r in b) == row.positive_band_planes
        for component in ("total", "electrostatic", "xc")
            av, spans, sig = Float64[], Float64[], Float64[]
            for r in p
                push!(av, parse(Float64, r["mean_"*component*"_ry"])*M.R.S.RY_EV)
                push!(spans, (parse(Float64, r["max_"*component*"_ry"]) -
                    parse(Float64, r["min_"*component*"_ry"]))*M.R.S.RY_EV)
                push!(sig, parse(Float64, r["std_"*component*"_ry"])*M.R.S.RY_EV)
            end
            for (kind, v) in (("mean", av), ("lateral_span", spans), ("lateral_sd", sig))
                @test getproperty(row, Symbol(component*"_"*kind*"_min_ev")) == minimum(v)
                @test getproperty(row, Symbol(component*"_"*kind*"_max_ev")) == maximum(v)
            end
        end
    end
    mktempdir() do dir
        file = joinpath(dir, "repeat.tsv")
        M.main([run, file])
        @test read(file) == read(saved)
        @test_throws ErrorException M.main([run, file])
    end
end

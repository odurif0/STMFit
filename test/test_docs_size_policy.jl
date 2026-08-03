#!/usr/bin/env julia

using Test

const MAKE_PATH = normpath(joinpath(@__DIR__, "..", "docs", "make.jl"))
const SOURCE = read(MAKE_PATH, String)
const HTML_CALL_PATTERN = r"Documenter\.HTML\s*\((?:[^()]|\([^()]*\))*\)"s

@testset "Documenter page-size policy" begin
    html_calls = collect(eachmatch(HTML_CALL_PATTERN, SOURCE))
    @test length(html_calls) == 1

    if length(html_calls) == 1
        html_call = only(html_calls).match
        @test occursin(
            r"size_threshold_ignore\s*=\s*\[\s*\"journal\.md\"\s*\]",
            html_call,
        )
        @test length(collect(eachmatch(r"\bsize_threshold_ignore\s*=", html_call))) == 1
    end

    @test length(collect(eachmatch(r"\bsize_threshold_ignore\s*=", SOURCE))) == 1
    @test !occursin(r"\bsize_threshold\s*=", SOURCE)
    @test !occursin(r"\bsize_threshold_warn\s*=", SOURCE)
end

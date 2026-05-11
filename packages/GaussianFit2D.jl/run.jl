#!/usr/bin/env julia

include(joinpath(@__DIR__, "src", "GaussianFit2D.jl"))
using .GaussianFit2D

GaussianFit2D.main_cli()

#!/usr/bin/env julia
#=
Multi-Gaussian fitting of STM/AFM line profiles.
Run: julia --project run.jl [data_file] [options]
=#

include(joinpath(@__DIR__, "src", "GaussianFit1D.jl"))
using .GaussianFit1D
GaussianFit1D.main_cli()

#!/usr/bin/env julia
#=
Multi-Gaussian fitting of STM/AFM line profiles.
Run: julia --project run.jl [data_file] [options]
=#

using GaussianFit1D
GaussianFit1D.main_cli()

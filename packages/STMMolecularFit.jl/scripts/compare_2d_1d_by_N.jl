#!/usr/bin/env julia
using STMMolecularFit
const F = length(ARGS) >= 1 ? ARGS[1] : error("Usage: julia --project=. scripts/compare_2d_1d_by_N.jl <filepath>")
STMMolecularFit.compare_2d_1d_by_N(F)

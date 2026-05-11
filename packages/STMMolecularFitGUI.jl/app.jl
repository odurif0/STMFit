"""
Multi-Gaussian Fitting GUI (Stipple.jl + PlotlyJS).

Usage:
    julia --project=. app.jl

Requires GaussianFit1D (core package) to be installed first:
    julia -e 'using Pkg; Pkg.develop(path="/home/durif/Git/GaussianFit1D.jl")'
"""

println("Starting Multi-Gaussian Fit GUI...")
println("  Loading packages (this may take a moment)...")

include(joinpath(@__DIR__, "src", "STMMolecularFitGUI.jl"))
using .STMMolecularFitGUI

println("  GUI module compiled.")

STMMolecularFitGUI.run()

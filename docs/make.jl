import Pkg
isempty(ARGS) || ARGS == ["--build-only"] || error("Usage: docs/make.jl [--build-only]")
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Documenter, GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore, STMSXMIO

makedocs(
    sitename = "STMFit",
    modules = [GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore, STMSXMIO],
    format = Documenter.HTML(
        prettyurls = false,
        edit_link = "main",
    ),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Pipeline" => "pipeline.md",
        "Model Selection" => "selection.md",
        "Unit Assignment" => "unit_assignment.md",
        "Calibration" => "calibration.md",
        "Configuration" => "config.md",
        "Runbook" => "chitosan_runbook.md",
        "DFT Molds" => "qe_stm_molds.md",
        "DFT Calculation Note" => "dft_calculation_note.md",
        "Mathematical Background" => "math.md",
        "Research Journal" => "journal.md",
        "API Reference" => "api.md",
    ],
)

if "--build-only" ∉ ARGS
    deploydocs(
        repo = "github.com/odurif0/STMFit.git",
        push_preview = true,
    )
end

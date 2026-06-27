import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Documenter, GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore, STMSXMIO

makedocs(
    sitename = "STMFit",
    modules = [GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore, STMSXMIO],
    format = Documenter.HTML(prettyurls = false, edit_link = "main"),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Pipeline" => "pipeline.md",
        "Calibration" => "calibration.md",
        "Chitosan Runbook" => "chitosan_runbook.md",
        "Unit Assignment" => "unit_assignment.md",
        "QE STM Molds" => "qe_stm_molds.md",
        "Model Selection" => "selection.md",
        "Mathematical Background" => "math.md",
        "Research Journal" => "journal.md",
        "Journal Archive" => "journal_archive.md",
        "Configuration" => "config.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/odurif0/STMFit.git",
    push_preview = true,
)

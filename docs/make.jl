using Documenter, GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore

makedocs(
    sitename = "STMFit",
    modules = [GaussianFit2D, GaussianFit1D, STMMolecularFit, STMFitCore],
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Pipeline" => "pipeline.md",
        "Model Selection" => "selection.md",
        "Mathematical Background" => "math.md",
        "Research Journal" => "journal.md",
        "Configuration" => "config.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(repo = "github.com/...")

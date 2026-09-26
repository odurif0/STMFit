#!/usr/bin/env julia
# Small saved tables only: all native PAW-free planes, no resampling or selection.
using Plots
include(joinpath(@__DIR__, "qe_density_cutoff.jl"))
const C = QEDensityCutoff

function main(args)
    args == ["--help"] && return println(
        "plot_qe_density_cutoff.jl SAVED_RUN glcn|glcnac NEW_PNG")
    VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
    length(args) == 3 || error("Use --help")
    run, molecule, out = args
    molecule in C.R.MOLECULES || error("Unknown molecule")
    ispath(out) && error("Output already exists")
    settings = C.settings(joinpath(run, "settings.toml"))
    dir = joinpath(run, molecule, "analysis")
    rows = C.R.V.tsv(joinpath(dir, "planes.tsv"))
    grids = C.R.V.tsv(joinpath(dir, "grids.tsv"))
    number(r, k) = parse(Float64, r[k])
    cutoffs = (settings["model"]["baseline_ecutrho_ry"],
        settings["model"]["candidate_ecutrho_ry"])
    midpoint = number(only(filter(r -> r["source"] == "baseline", grids)), "midpoint_nm")
    selected = Dict()
    for cutoff in cutoffs
        allrows = filter(r -> number(r, "ecutrho_ry") == cutoff, rows)
        grid = only(filter(r -> number(r, "ecutrho_ry") == cutoff, grids))
        number.(allrows, "k") == collect(0:Int(number(grid, "nz"))-1) ||
            error("Missing, duplicate or unordered native planes")
        selected[cutoff] = filter(r -> parse(Bool, r["paw_gap"]), allrows)
        isempty(selected[cutoff]) && error("Empty declared PAW-free domain")
    end
    fields = (("std_total_ry", "Total potential", C.R.S.RY_EV, "population spatial SD (eV)"),
        ("std_electrostatic_ry", "Electrostatic potential", C.R.S.RY_EV, "population spatial SD (eV)"),
        ("std_xc_ry", "XC = total - electrostatic", C.R.S.RY_EV, "population spatial SD (eV)"),
        ("negative_density_fraction", "Negative valence-density samples", 100., "fraction of native plane (%)"))
    panels = []
    for (field, title, scale, ylabel) in fields
        p = plot(; title, xlabel="accepted-cell z (nm)", ylabel,
            legend=:outertopright, left_margin=8Plots.mm, bottom_margin=7Plots.mm)
        for (cutoff, color, style) in zip(cutoffs, (:royalblue, :darkorange), (:solid, :dash))
            points = selected[cutoff]
            plot!(p, number.(points, "z_nm"), scale .* number.(points, field);
                color, linestyle=style, linewidth=1.4, marker=:circle, markersize=1.5,
                markerstrokewidth=0, label="$(Int(cutoff)) Ry ($(length(points)) planes)")
        end
        vline!(p, [midpoint]; color=:gray, linestyle=:dot, label="gap midpoint")
        field == "negative_density_fraction" && ylims!(p, (0, 100))
        push!(panels, p)
    end
    mkpath(dirname(out))
    savefig(plot(panels...; layout=(2, 2), size=(1600, 1100),
        plot_title="$molecule — cutoff sensitivity on native planes; not a convergence proof"), out)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)

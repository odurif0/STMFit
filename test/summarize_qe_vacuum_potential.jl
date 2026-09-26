#!/usr/bin/env julia
# Descriptive ranges only: no plane, barrier or chemical reference is selected.
module SummarizeQEVacuumPotential
include(joinpath(@__DIR__, "qe_vacuum_potential.jl"))
const R = QEVacuumPotential
number(r, k) = parse(Float64, r[k])
flag(r, k) = parse(Bool, r[k])

function summaries(run)
    R.settings(joinpath(run, "settings.toml"))
    records = NamedTuple[]
    for molecule in R.MOLECULES
        dir = joinpath(run, molecule)
        planes = R.V.tsv(joinpath(dir, "analysis/planes.tsv"))
        barriers = R.V.tsv(joinpath(dir, "analysis/barriers.tsv"))
        ef = R.TOML.parsefile(joinpath(dir, "metadata.toml"))["fermi_ev"]
        for domain in ("paw_gap", "lower_half")
            rows = filter(r -> flag(r, domain), planes)
            bands = filter(r -> flag(r, domain), barriers)
            isempty(rows) && error("Empty declared domain; no alternative is selected")
            positive = filter(r -> flag(r, "all_sampled_barriers_positive"), bands)
            isempty(positive) && error("No positive sampled barriers; inspect original table")
            nband = length(unique(number.(bands, "band")))
            length(bands) == nband * length(rows) || error("Incomplete band/plane table")
            stats = (; molecule, domain, planes=length(rows), selected_bands=nband,
                band_planes=length(bands), positive_band_planes=length(positive),
                z_min_nm=minimum(number.(rows, "z_nm")),
                z_max_nm=maximum(number.(rows, "z_nm")),
                mean_total_minus_fermi_min_ev=minimum(number.(rows, "mean_total_minus_fermi_ev")),
                mean_total_minus_fermi_max_ev=maximum(number.(rows, "mean_total_minus_fermi_ev")),
                global_total_minus_fermi_min_ev=minimum(number.(rows, "min_total_ry"))*R.S.RY_EV-ef,
                global_total_minus_fermi_max_ev=maximum(number.(rows, "max_total_ry"))*R.S.RY_EV-ef,
                barrier_global_min_ev=minimum(number.(bands, "barrier_min_ev")),
                barrier_global_max_ev=maximum(number.(bands, "barrier_max_ev")),
                positive_span_over_mean_barrier_min=minimum(number.(positive, "lateral_span_over_mean_barrier")),
                positive_span_over_mean_barrier_max=maximum(number.(positive, "lateral_span_over_mean_barrier")),
                positive_kappa_min_nm_inv=minimum(number.(positive, "kappa_min_nm_inv")),
                positive_kappa_max_nm_inv=maximum(number.(positive, "kappa_max_nm_inv")))
            for component in ("total", "electrostatic", "xc")
                values = (
                    "mean" => number.(rows, "mean_$(component)_ry") .* R.S.RY_EV,
                    "lateral_span" => (number.(rows, "max_$(component)_ry") .-
                        number.(rows, "min_$(component)_ry")) .* R.S.RY_EV,
                    "lateral_sd" => number.(rows, "std_$(component)_ry") .* R.S.RY_EV)
                for (name, v) in values
                    key = component * "_" * name
                    stats = merge(stats, (; Symbol(key*"_min_ev")=>minimum(v),
                        Symbol(key*"_max_ev")=>maximum(v)))
                end
            end
            push!(records, stats)
        end
    end
    records
end

function main(args)
    length(args) == 2 || error("summarize_qe_vacuum_potential.jl SAVED_RUN NEW_TSV")
    run, out = args
    ispath(out) && error("Output already exists")
    R.S.table(out, summaries(run))
    println("Four declared-domain summaries; positive-only descriptors explicitly conditional; no selection")
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && SummarizeQEVacuumPotential.main(ARGS)

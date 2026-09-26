#!/usr/bin/env julia
# Small-table report only; no volume loading, matching plane or potential change.
module SummarizeQEXCComponents
using TOML
include(joinpath(@__DIR__, "qe_xc_build_control.jl"))
const B = QEXCBuildControl
const DOMAINS = ("all_planes", "paw_gap", "lower_half")
const SUBSETS = ("all", "gga_active", "gga_inactive")
const COMPONENTS = ("xc_total", "lda", "gga_local", "gga_divergence")
number(r, k) = parse(Float64, r[k])
flag(r, k) = parse(Bool, r[k])
in_domain(r, d) = d == "all_planes" || flag(r, d)
extent(v) = isempty(v) ? (NaN, NaN) : extrema(v)

# For the saved native fields, XC = lda + local + divergence + error.
# Centering and the triangle inequality give this upper bound. It is not a
# fraction of explained variance, a physical tolerance, or a selection rule.
function residual_bound(sd_xc, sd_lda, sd_local, max_additivity_error)
    all(x -> isfinite(x) && x >= 0, (sd_xc, sd_lda, sd_local, max_additivity_error)) ||
        error("Invalid standard deviation or additivity bound")
    iszero(sd_xc) ? NaN : (sd_lda + sd_local + max_additivity_error) / sd_xc
end

function summarize_tables(molecule, rows, counts, max_additivity_error)
    components = NamedTuple[]; branches = NamedTuple[]; bounds = NamedTuple[]
    for domain in DOMAINS
        cc = filter(r -> in_domain(r, domain), counts)
        isempty(cc) && error("Empty declared domain; no alternative selected")
        keys = Set(r["k"] for r in cc)
        length(keys) == length(cc) || error("Duplicate count plane")
        location = (; molecule, domain, planes=length(cc),
            z_min_nm=minimum(number.(cc, "z_nm")), z_max_nm=maximum(number.(cc, "z_nm")))
        branch = location
        for field in ("gga_active", "negative_valence", "negative_xc_input", "below_lda", "below_gga_density")
            fraction = number.(cc, field) ./ number.(cc, "samples")
            all(x -> isfinite(x) && 0 <= x <= 1, fraction) || error("Invalid fraction")
            lo, hi = extrema(fraction)
            branch = merge(branch, (; Symbol(field*"_fraction_min")=>lo, Symbol(field*"_fraction_max")=>hi))
        end
        push!(branches, branch)
        for subset in SUBSETS
            lookup = Dict{String,Any}()
            for component in COMPONENTS
                rr = filter(r -> in_domain(r, domain) && r["subset"] == subset && r["component"] == component, rows)
                Set(r["k"] for r in rr) == keys && length(rr) == length(cc) || error("Incomplete component planes")
                nonempty = filter(r -> number(r, "samples") > 0, rr)
                empty = filter(r -> number(r, "samples") == 0, rr)
                all(r -> all(k -> isnan(number(r,k)), ("mean", "minimum", "maximum", "spatial_sd")), empty) || error("Invalid empty subset")
                record = merge(location, (; subset, component, nonempty_planes=length(nonempty), empty_planes=length(empty)))
                for stat in ("mean", "minimum", "maximum", "spatial_sd")
                    values = number.(nonempty, stat) .* B.S.RY_EV
                    all(isfinite, values) || error("Nonfinite nonempty statistic")
                    lo, hi = extent(values)
                    record = merge(record, (; Symbol(stat*"_min_ev")=>lo, Symbol(stat*"_max_ev")=>hi))
                end
                push!(components, record)
                lookup[component] = Dict(r["k"]=>r for r in rr)
            end
            ratios = Float64[]; empty_planes = 0; zero_xc_sd_planes = 0
            for k in keys
                rr = [lookup[c][k] for c in COMPONENTS]
                length(unique(number.(rr, "samples"))) == 1 || error("Different component subsets")
                if number(rr[1], "samples") == 0
                    empty_planes += 1
                else
                    value = residual_bound(number(rr[1], "spatial_sd"), number(rr[2], "spatial_sd"),
                        number(rr[3], "spatial_sd"), max_additivity_error)
                    isnan(value) ? (zero_xc_sd_planes += 1) : push!(ratios, value)
                end
            end
            lo, hi = extent(ratios)
            push!(bounds, merge(location, (; subset, defined_planes=length(ratios), empty_planes,
                zero_xc_sd_planes, centered_residual_relative_l2_upper_bound_min=lo,
                centered_residual_relative_l2_upper_bound_max=hi, max_additivity_error_ry=max_additivity_error)))
        end
    end
    (; components, branches, bounds)
end

function main(args)
    length(args) == 2 || error("summarize_qe_xc_components.jl SAVED_PAIRED_RUN NEW_DIRECTORY")
    run, out = args
    ispath(out) && error("Output exists")
    B.settings(joinpath(run, "settings.toml"))
    reports = map(B.R.MOLECULES) do molecule
        dir = joinpath(run, molecule)
        controls = TOML.parsefile(joinpath(dir, "controls/summary.toml"))
        analysis = TOML.parsefile(joinpath(dir, "analysis/summary.toml"))
        controls["same_build_controls_exact"] && analysis["exact_same_build_controls_and_repeats"] || error("Unqualified controls")
        !analysis["density_modified"] && !analysis["potential_adopted"] || error("Changed physical experiment")
        summarize_tables(molecule, B.V.tsv(joinpath(dir, "analysis/components.tsv")),
            B.V.tsv(joinpath(dir, "analysis/counts.tsv")), analysis["additivity"]["max_abs_error_ry"])
    end
    B.S.newdir(out)
    for key in (:components, :branches, :bounds)
        B.S.table(joinpath(out, string(key)*".tsv"), reduce(vcat, getproperty.(reports, key)))
    end
    println("All three domains and source subsets retained; triangle bounds are not variance attribution or selection")
end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && SummarizeQEXCComponents.main(ARGS)

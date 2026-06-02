#!/usr/bin/env julia
# ──────────────────────────────────────────────────────────────────────────────
# aggregate_synthetic_known_n.jl
#
# Aggregate synthetic known-N validation TSV results into a compact
# per-policy / per-group comparison table.
#
# Usage:
#   julia --project=. test/aggregate_synthetic_known_n.jl
#       [input.tsv] [--out output.tsv]
#
# Default input:  results/synthetic_known_n/summary.tsv
# Default output: results/synthetic_known_n/aggregate.tsv
# ──────────────────────────────────────────────────────────────────────────────

using DelimitedFiles, Printf, Statistics

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

const DEFAULT_IN  = "results/synthetic_known_n/summary.tsv"
const DEFAULT_OUT = "results/synthetic_known_n/aggregate.tsv"

function _parse_cli()
    inputs = String[]   # positional arguments = input TSV paths
    output = DEFAULT_OUT

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--out" && i < length(ARGS)
            output = ARGS[i+1]; i += 2
        elseif !startswith(arg, "-")
            push!(inputs, arg); i += 1
        else
            i += 1
        end
    end

    if isempty(inputs)
        push!(inputs, DEFAULT_IN)
    end

    return inputs, output
end

# ══════════════════════════════════════════════════════════════════════════════
# Reading
# ══════════════════════════════════════════════════════════════════════════════

struct SummaryRow
    case_id::String
    seed::Int
    true_N::Int
    artifact::String
    policy::String
    N_eff::Int
    N_selected::Union{Int,Nothing}  # nothing ↔ "NA"
    abs_error::Union{Float64,Nothing}
    status::String
    score_or_source::String
    noise_scale::Float64            # 1.0 for legacy 10-column TSVs
    mode::String                    # "circular" for legacy / default
end

function _read_tsv(path::String)
    data = readdlm(path, '\t'; skipstart=1)
    rows = SummaryRow[]
    for i in 1:size(data, 1)
        line = data[i, :]
        ncol = length(line)
        ncol >= 9 || continue

        case_id = string(line[1])
        seed    = parse(Int, string(line[2]))
        true_N  = parse(Int, string(line[3]))
        artifact = string(line[4])
        policy   = string(line[5])

        n_eff_str = string(line[6])
        n_eff = n_eff_str == "NA" ? 0 : parse(Int, n_eff_str)

        n_sel_str = string(line[7])
        n_selected = n_sel_str == "NA" ? nothing : parse(Int, n_sel_str)

        abs_str = string(line[8])
        abs_error = (abs_str == "NA" || abs_str == "NaN") ? nothing : parse(Float64, abs_str)

        status = string(line[9])
        score_or_source = ncol >= 10 ? string(line[10]) : ""

        # noise_scale: column 11 in new-format TSVs; default 1.0 for legacy 10-col
        noise_scale = 1.0
        if ncol >= 11
            ns_str = string(line[11])
            noise_scale = (ns_str == "NA" || ns_str == "NaN") ? 1.0 : parse(Float64, ns_str)
        end

        # mode: column 12 in newest-format TSVs; default "circular" for legacy
        mode = "circular"
        if ncol >= 12
            mode = string(line[12])
        end

        push!(rows, SummaryRow(case_id, seed, true_N, artifact, policy,
                               n_eff, n_selected, abs_error, status, score_or_source,
                               noise_scale, mode))
    end
    return rows
end

# ══════════════════════════════════════════════════════════════════════════════
# Metrics
# ══════════════════════════════════════════════════════════════════════════════

struct GroupMetrics
    group::String
    policy::String
    n_ok::Int
    n_total::Int
    exact::Int
    exact_rate::Float64
    mean_abs_error::Float64
    over_selected::Int
    under_selected::Int
    error_count::Int
end

function _compute_metrics(rows::Vector{SummaryRow})
    # Collect distinct groups
    true_Ns = sort(unique(r.true_N for r in rows))
    artifacts = sort(unique(r.artifact for r in rows))
    seeds = sort(unique(r.seed for r in rows))
    noise_scales = sort(unique(r.noise_scale for r in rows))
    modes = sort(unique(r.mode for r in rows))

    groups = String["all"]
    for tn in true_Ns; push!(groups, "true_N=$(tn)"); end
    for a in artifacts;  push!(groups, "artifact=$(a)"); end
    for s in seeds;      push!(groups, "seed=$(s)"); end
    for ns in noise_scales; push!(groups, @sprintf("noise_scale=%g", ns)); end
    for m in modes;       push!(groups, "mode=$(m)"); end

    results = GroupMetrics[]

    for g in groups
        # Filter rows belonging to this group
        subset = filter(rows) do r
            if g == "all"
                true
            elseif startswith(g, "true_N=")
                parse(Int, g[8:end]) == r.true_N
            elseif startswith(g, "artifact=")
                g[10:end] == r.artifact
            elseif startswith(g, "seed=")
                parse(Int, g[6:end]) == r.seed
            elseif startswith(g, "noise_scale=")
                abs(r.noise_scale - parse(Float64, g[13:end])) < 1e-9
            elseif startswith(g, "mode=")
                g[6:end] == r.mode
            else
                false
            end
        end

        for policy in sort(unique(r.policy for r in subset))
            pol_rows = filter(r -> r.policy == policy, subset)
            isempty(pol_rows) && continue

            n_total = length(pol_rows)
            ok_rows = filter(r -> r.status == "ok", pol_rows)
            n_ok = length(ok_rows)
            error_count = count(r -> startswith(r.status, "error"), pol_rows)

            exact = count(r -> r.status == "ok" && r.N_selected !== nothing &&
                              r.N_selected == r.true_N, pol_rows)
            exact_rate = n_ok > 0 ? exact / n_ok : NaN

            abs_vals = [r.abs_error for r in ok_rows if r.abs_error !== nothing]
            mean_ae = isempty(abs_vals) ? NaN : mean(abs_vals)

            over  = count(r -> r.status == "ok" && r.N_selected !== nothing &&
                              r.N_selected > r.true_N, pol_rows)
            under = count(r -> r.status == "ok" && r.N_selected !== nothing &&
                              r.N_selected < r.true_N, pol_rows)

            push!(results, GroupMetrics(g, policy, n_ok, n_total, exact,
                                        exact_rate, mean_ae, over, under, error_count))
        end
    end

    return results
end

# ══════════════════════════════════════════════════════════════════════════════
# Formatting & output
# ══════════════════════════════════════════════════════════════════════════════

const AGG_HEADER = [
    "group", "policy",
    "n_ok", "n_total", "exact", "exact_rate",
    "mean_abs_error", "over_selected", "under_selected", "errors",
]

function _fmt_val(v::Float64)
    isnan(v) && return "NA"
    return @sprintf("%.4f", v)
end

function _metrics_row(m::GroupMetrics)::Vector{String}
    return [
        m.group,
        m.policy,
        string(m.n_ok),
        string(m.n_total),
        string(m.exact),
        _fmt_val(m.exact_rate),
        _fmt_val(m.mean_abs_error),
        string(m.over_selected),
        string(m.under_selected),
        string(m.error_count),
    ]
end

function _print_table(metrics::Vector{GroupMetrics})
    # Header
    println(join(AGG_HEADER, "\t"))
    for m in metrics
        println(join(_metrics_row(m), "\t"))
    end
end

function _write_tsv(path::String, metrics::Vector{GroupMetrics})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(AGG_HEADER, "\t"))
        for m in metrics
            println(io, join(_metrics_row(m), "\t"))
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    input_paths, output_path = _parse_cli()

    all_rows = SummaryRow[]
    for p in input_paths
        isfile(p) || error("Input TSV not found: $p")
        @printf("Reading: %s\n", p)
        rows = _read_tsv(p)
        @printf("  %d rows loaded\n", length(rows))
        append!(all_rows, rows)
    end

    @printf("\nTotal rows: %d across %d file(s)\n", length(all_rows), length(input_paths))

    metrics = _compute_metrics(all_rows)
    @printf("  %d group×policy rows computed\n\n", length(metrics))

    _print_table(metrics)
    _write_tsv(output_path, metrics)
    @printf("\nWritten: %s\n", output_path)

    return metrics
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

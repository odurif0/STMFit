#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# audit_refined_selection.jl — External audit of conservative selection guards
#
# This script combines an existing primary selection TSV with an advisory TSV.
# It is NOT used by fitting.  The first supported policy is an overfit guard:
# keep the primary N unless an advisory criterion selects a LOWER N.  This tests
# whether a generic robust-IC diagnostic can correct over-splitting without
# allowing upward changes that would break clean controls such as 240817_026.
#
# No expected-N labels are used by this script.  Benchmark grading remains a
# separate external evaluation step.
# ──────────────────────────────────────────────────────────────────────────────

using Printf

const DEFAULT_PRIMARY = "results/best_plots/summary_overlap060_hard.tsv"
const DEFAULT_ADVISORY = "results/robust_rescore_audit/full_aicc_nu8.tsv"
const DEFAULT_OUT = "results/refined_selection/overfit_guard.tsv"

function _parse_cli(args)
    primary = DEFAULT_PRIMARY
    advisory = DEFAULT_ADVISORY
    out_tsv = DEFAULT_OUT
    primary_file_col = "filepath"
    primary_n_col = "N_eff"
    advisory_file_col = "file"
    advisory_n_col = "N"
    advisory_selected_col = "is_selected"
    policy = "down_only"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--primary"
            primary = args[i+1]; i += 2
        elseif startswith(arg, "--primary=")
            primary = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--advisory"
            advisory = args[i+1]; i += 2
        elseif startswith(arg, "--advisory=")
            advisory = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--primary-file-column"
            primary_file_col = args[i+1]; i += 2
        elseif arg == "--primary-n-column"
            primary_n_col = args[i+1]; i += 2
        elseif arg == "--advisory-file-column"
            advisory_file_col = args[i+1]; i += 2
        elseif arg == "--advisory-n-column"
            advisory_n_col = args[i+1]; i += 2
        elseif arg == "--advisory-selected-column"
            advisory_selected_col = args[i+1]; i += 2
        elseif arg == "--policy"
            policy = args[i+1]; i += 2
        elseif startswith(arg, "--policy=")
            policy = split(arg, "=", limit=2)[2]; i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/audit_refined_selection.jl [options]

            Options:
              --primary PATH                 Primary TSV [$(DEFAULT_PRIMARY)]
              --advisory PATH                Advisory selected-row TSV [$(DEFAULT_ADVISORY)]
              --out PATH                     Output TSV [$(DEFAULT_OUT)]
              --primary-file-column NAME     [filepath]
              --primary-n-column NAME        [N_eff]
              --advisory-file-column NAME    [file]
              --advisory-n-column NAME       [N]
              --advisory-selected-column NAME [is_selected]
              --policy down_only             Conservative overfit guard
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    policy == "down_only" || error("Only --policy down_only is implemented")
    return (; primary, advisory, out_tsv, primary_file_col, primary_n_col,
            advisory_file_col, advisory_n_col, advisory_selected_col, policy)
end

function _read_tsv(path)
    lines = readlines(path)
    isempty(lines) && error("Empty TSV: $path")
    header = split(lines[1], '\t'; keepempty=true)
    rows = Vector{Dict{String,String}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (i, h) in enumerate(header)
            row[h] = i <= length(vals) ? vals[i] : ""
        end
        push!(rows, row)
    end
    return rows
end

function _basename_file(s)
    b = basename(strip(s))
    m = match(r"(240817_\d+\.sxm)", b)
    return m === nothing ? b : m.captures[1]
end

_as_bool(s) = lowercase(strip(s)) in ("true", "t", "1", "yes", "y")
_parse_int(s) = parse(Int, strip(s))

function main(args=ARGS)
    opt = _parse_cli(args)
    primary_rows = _read_tsv(opt.primary)
    advisory_rows = _read_tsv(opt.advisory)

    advisory_n = Dict{String,Int}()
    advisory_source = Dict{String,String}()
    for row in advisory_rows
        haskey(row, opt.advisory_selected_col) && _as_bool(row[opt.advisory_selected_col]) || continue
        file = _basename_file(row[opt.advisory_file_col])
        advisory_n[file] = _parse_int(row[opt.advisory_n_col])
        advisory_source[file] = get(row, "source", "advisory")
    end

    mkpath(dirname(opt.out_tsv))
    n_changed = 0
    open(opt.out_tsv, "w") do io
        println(io, join(["file", "N_refined", "is_selected", "N_primary", "N_advisory",
                          "policy", "decision", "advisory_source"], '\t'))
        for row in primary_rows
            file = _basename_file(row[opt.primary_file_col])
            n_primary = _parse_int(row[opt.primary_n_col])
            n_adv = get(advisory_n, file, n_primary)
            n_refined = n_adv < n_primary ? n_adv : n_primary
            decision = n_refined == n_primary ? "keep_primary" : "downshift_advisory"
            decision == "downshift_advisory" && (n_changed += 1)
            println(io, join([file, n_refined, true, n_primary, n_adv, opt.policy,
                              decision, get(advisory_source, file, "")], '\t'))
        end
    end
    println("Refined selection audit")
    println("  primary:  ", opt.primary)
    println("  advisory: ", opt.advisory)
    println("  policy:   ", opt.policy)
    println("  changed:  ", n_changed)
    println("  output:   ", opt.out_tsv)
end

main()

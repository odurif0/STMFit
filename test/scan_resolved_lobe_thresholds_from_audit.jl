#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# scan_resolved_lobe_thresholds_from_audit.jl
#
# Fast external evaluation of resolved-lobe thresholds from one existing
# audit_resolved_lobes.jl TSV.  This avoids refitting for every threshold: the
# audit TSV already stores threshold-independent adjacent-pair metrics in
# `pair_details`, so this script recomputes N_resolved for many threshold
# combinations offline and grades each candidate against the benchmark manifest.
#
# This is methodological evaluation only.  Human validation labels are used here
# to compare generic post-processing rules, never inside fit-time selection.
# ──────────────────────────────────────────────────────────────────────────────

using Printf, TOML

const DEFAULT_AUDIT = "results/resolved_lobes_scan_sep/sep1_vsnr0_vfrac0/resolved_lobes.tsv"
const DEFAULT_MANIFEST = "benchmarks/chitosan_240817.toml"
const DEFAULT_OUT = "results/resolved_lobes_scan_offline/threshold_scores.tsv"

function _parse_csv_floats(s::AbstractString)
    vals = Float64[]
    for p in split(s, ',')
        v = tryparse(Float64, strip(p))
        v === nothing && error("Invalid float in grid: $p")
        push!(vals, v)
    end
    return vals
end

function _parse_cli(args)
    audit = DEFAULT_AUDIT
    manifest = DEFAULT_MANIFEST
    out_tsv = DEFAULT_OUT
    sep_grid = collect(1.0:0.05:2.0)
    vsnr_grid = [0.0, 0.5, 1.0, 2.0, 3.0]
    vfrac_grid = [0.0, 0.02, 0.05, 0.1, 0.2]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--audit"
            i < length(args) || error("--audit requires a path")
            audit = args[i+1]; i += 2
        elseif startswith(arg, "--audit=")
            audit = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--manifest"
            i < length(args) || error("--manifest requires a path")
            manifest = args[i+1]; i += 2
        elseif startswith(arg, "--manifest=")
            manifest = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out"
            i < length(args) || error("--out requires a path")
            out_tsv = args[i+1]; i += 2
        elseif startswith(arg, "--out=")
            out_tsv = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--sep-grid"
            i < length(args) || error("--sep-grid requires comma-separated floats")
            sep_grid = _parse_csv_floats(args[i+1]); i += 2
        elseif startswith(arg, "--sep-grid=")
            sep_grid = _parse_csv_floats(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--valley-snr-grid"
            i < length(args) || error("--valley-snr-grid requires comma-separated floats")
            vsnr_grid = _parse_csv_floats(args[i+1]); i += 2
        elseif startswith(arg, "--valley-snr-grid=")
            vsnr_grid = _parse_csv_floats(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--valley-frac-grid"
            i < length(args) || error("--valley-frac-grid requires comma-separated floats")
            vfrac_grid = _parse_csv_floats(args[i+1]); i += 2
        elseif startswith(arg, "--valley-frac-grid=")
            vfrac_grid = _parse_csv_floats(split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/scan_resolved_lobe_thresholds_from_audit.jl [options]

            Options:
              --audit PATH              Existing resolved_lobes.tsv [$(DEFAULT_AUDIT)]
              --manifest PATH           Benchmark manifest [$(DEFAULT_MANIFEST)]
              --out PATH                Aggregate output TSV [$(DEFAULT_OUT)]
              --sep-grid V1,V2          Separation thresholds [1.0:0.05:2.0]
              --valley-snr-grid V1,V2   Valley SNR thresholds [0,0.5,1,2,3]
              --valley-frac-grid V1,V2  Valley fraction thresholds [0,0.02,0.05,0.1,0.2]
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    return audit, manifest, out_tsv, sep_grid, vsnr_grid, vfrac_grid
end

function _read_tsv(path::AbstractString)
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

_as_bool(s::AbstractString) = lowercase(strip(s)) in ("true", "t", "1", "yes", "y")
_parse_int(s::AbstractString) = parse(Int, strip(s))
function _basename_file(s::AbstractString)
    b = basename(strip(s))
    m = match(r"(240817_\d+\.sxm)", b)
    return m === nothing ? b : m.captures[1]
end

function _quality(manifest, file)
    info = get(get(manifest, "files", Dict{String,Any}()), file, Dict{String,Any}())
    q = String(get(info, "quality", "clean"))
    excl = Bool(get(info, "exclude_primary", false))
    target = get(info, "target_N_for_human_validation", 6)
    return q, excl, Int(target)
end

function _parse_pairs(pair_details::AbstractString)
    pairs = Tuple{Int,Int,Float64,Float64,Float64}[]
    isempty(strip(pair_details)) && return pairs
    for item in split(pair_details, ';')
        isempty(strip(item)) && continue
        lhs_rhs = split(item, ':'; limit=2)
        length(lhs_rhs) == 2 || continue
        ij = split(lhs_rhs[1], '-')
        vals = split(lhs_rhs[2], '/')
        length(ij) == 2 && length(vals) >= 4 || continue
        i = parse(Int, ij[1])
        j = parse(Int, ij[2])
        sep_sigma = parse(Float64, vals[2])
        valley_snr = parse(Float64, vals[3])
        valley_frac = parse(Float64, vals[4])
        push!(pairs, (i, j, sep_sigma, valley_snr, valley_frac))
    end
    return pairs
end

function _n_resolved(n::Int, pairs, sep::Float64, vsnr::Float64, vfrac::Float64)
    parent = collect(1:n)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[rb] = ra)
    end
    for (i, j, sep_sigma, valley_snr, valley_frac) in pairs
        unresolved = (sep_sigma < sep) && (valley_snr < vsnr || valley_frac < vfrac)
        unresolved && union!(i, j)
    end
    return length(Set(find(i) for i in 1:n))
end

function _grade(selected_rows, manifest, sep, vsnr, vfrac)
    eval_cfg = get(manifest, "evaluation", Dict{String,Any}())
    primary_excl_q = Set(String.(get(eval_cfg, "primary_exclude_quality", ["poor_quality", "excluded"])))
    stress_q = Set(String.(get(eval_cfg, "stress_quality", ["poor_quality"])))
    primary_ok = primary_total = target_ok = target_total = stress_total = 0
    primary_fail = String[]
    target_fail = String[]
    observed_by_file = String[]
    for row in selected_rows
        file = _basename_file(row["file"])
        n = _parse_int(row["N"])
        observed = _n_resolved(n, _parse_pairs(get(row, "pair_details", "")), sep, vsnr, vfrac)
        q, excl, target = _quality(manifest, file)
        is_primary = !(excl || q in primary_excl_q)
        is_target = q == "clean_target"
        is_stress = q in stress_q
        passed = observed == target
        push!(observed_by_file, string(file, ":", observed))
        if is_primary
            primary_total += 1
            passed ? (primary_ok += 1) : push!(primary_fail, file)
        end
        if is_target
            target_total += 1
            passed ? (target_ok += 1) : push!(target_fail, file)
        end
        is_stress && (stress_total += 1)
    end
    return (; sep, vsnr, vfrac, primary_ok, primary_total, target_ok, target_total,
            stress_total, primary_fail=join(primary_fail, ";"),
            target_fail=join(target_fail, ";"), observed=join(observed_by_file, ";"))
end

function _pct(ok, total)
    total == 0 && return "NA"
    return @sprintf("%.1f", 100ok / total)
end

function main(args=ARGS)
    audit, manifest_path, out_tsv, sep_grid, vsnr_grid, vfrac_grid = _parse_cli(args)
    manifest = TOML.parsefile(manifest_path)
    rows = _read_tsv(audit)
    selected_rows = [r for r in rows if _as_bool(get(r, "is_selected", "false"))]
    isempty(selected_rows) && error("No selected rows in audit TSV: $audit")
    mkpath(dirname(out_tsv))

    grades = Any[]
    for sep in sep_grid, vsnr in vsnr_grid, vfrac in vfrac_grid
        push!(grades, _grade(selected_rows, manifest, sep, vsnr, vfrac))
    end
    sort!(grades; by = g -> (-g.primary_ok, -g.target_ok, g.sep, g.vsnr, g.vfrac))

    open(out_tsv, "w") do io
        println(io, join(["rank", "sep_threshold", "valley_snr_threshold", "valley_frac_threshold",
                          "primary_ok", "primary_total", "primary_pct", "target_ok",
                          "target_total", "target_pct", "stress_total", "primary_failures",
                          "target_failures", "observed_by_file"], '\t'))
        for (rank, g) in enumerate(grades)
            println(io, join([rank, @sprintf("%.4g", g.sep), @sprintf("%.4g", g.vsnr), @sprintf("%.4g", g.vfrac),
                              g.primary_ok, g.primary_total, _pct(g.primary_ok, g.primary_total),
                              g.target_ok, g.target_total, _pct(g.target_ok, g.target_total),
                              g.stress_total, g.primary_fail, g.target_fail, g.observed], '\t'))
        end
    end
    best = first(grades)
    println("Offline resolved-lobe threshold scan")
    println("  audit:   ", audit)
    println("  output:  ", out_tsv)
    println("  combos:  ", length(grades))
    println("  best:    primary ", best.primary_ok, "/", best.primary_total,
            " target ", best.target_ok, "/", best.target_total,
            " sep=", best.sep, " vsnr=", best.vsnr, " vfrac=", best.vfrac)
    !isempty(best.primary_fail) && println("  fails:   ", best.primary_fail)
end

main()

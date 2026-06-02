#!/usr/bin/env julia

# ──────────────────────────────────────────────────────────────────────────────
# scan_resolved_lobe_thresholds.jl — Generic resolved-lobe threshold scan
#
# This script evaluates how the N_resolved audit criterion behaves across a
# grid of (sep_threshold, valley_snr_threshold, valley_frac_threshold) values.
# For each threshold combination it:
#   1. runs audit_resolved_lobes.jl as a subprocess (fitting + pair metrics)
#   2. grades the selected rows against the chitosan benchmark manifest
#   3. writes an aggregate TSV across all threshold combinations
#
# IMPORTANT: This is external method evaluation, not fit-time selection.  The
# fitted chain models are the same for every threshold combination; only the
# resolvability post-processing changes.
#
# Usage:
#   julia --project=. test/scan_resolved_lobe_thresholds.jl \
#       --config config/chitosan.toml \
#       --sep-grid 1.0,1.5,2.0 \
#       --valley-snr-grid 1.0,3.0 \
#       --valley-frac-grid 0.1,0.2 \
#       --files 240817_019.sxm,240817_058.sxm \
#       --outdir results/resolved_lobes_scan
#
# For a dry run (prints what would be done without executing):
#   ... --dry-run
# ──────────────────────────────────────────────────────────────────────────────

using Printf, TOML

# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════

const SCRIPTS_DIR = @__DIR__
const PROJECT_DIR = dirname(SCRIPTS_DIR)
const AUDIT_SCRIPT = joinpath(SCRIPTS_DIR, "audit_resolved_lobes.jl")
const DEFAULT_MANIFEST = joinpath(PROJECT_DIR, "benchmarks", "chitosan_240817.toml")
const DEFAULT_CONFIG = "config/chitosan.toml"
const DEFAULT_OUTDIR = "results/resolved_lobes_scan"

const DEFAULT_TARGETS = [
    "240817_017.sxm", "240817_019.sxm", "240817_043.sxm", "240817_058.sxm"
]
const DEFAULT_CONTROLS = [
    "240817_002.sxm", "240817_003.sxm", "240817_004.sxm", "240817_005.sxm",
    "240817_018.sxm", "240817_021.sxm", "240817_039.sxm", "240817_060.sxm"
]
const DEFAULT_STRESS = [
    "240817_029.sxm", "240817_030.sxm", "240817_031.sxm", "240817_032.sxm",
    "240817_034.sxm", "240817_035.sxm", "240817_037.sxm", "240817_038.sxm",
    "240817_051.sxm"
]

# Small default grid (2 combos) to keep smoke tests fast
const DEFAULT_SEP_GRID = [1.0, 2.0]
const DEFAULT_VSNR_GRID = [3.0]
const DEFAULT_VFRAC_GRID = [0.2]

# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

function _parse_csv_floats(s::AbstractString)
    parts = split(s, ',')
    vals = Float64[]
    for p in parts
        v = tryparse(Float64, strip(p))
        v === nothing && error("Invalid float in grid: '$p'")
        push!(vals, v)
    end
    return vals
end

function _parse_cli(args)
    config_file = DEFAULT_CONFIG
    manifest_file = DEFAULT_MANIFEST
    outdir = DEFAULT_OUTDIR
    sep_grid = copy(DEFAULT_SEP_GRID)
    vsnr_grid = copy(DEFAULT_VSNR_GRID)
    vfrac_grid = copy(DEFAULT_VFRAC_GRID)
    files = vcat(DEFAULT_TARGETS, DEFAULT_CONTROLS, DEFAULT_STRESS)
    dry_run = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config"
            i < length(args) || error("--config requires a file path")
            config_file = args[i+1]; i += 2
        elseif startswith(arg, "--config=")
            config_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--manifest"
            i < length(args) || error("--manifest requires a file path")
            manifest_file = args[i+1]; i += 2
        elseif startswith(arg, "--manifest=")
            manifest_file = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--files"
            i < length(args) || error("--files requires comma-separated names")
            files = split(strip(args[i+1]), ','); i += 2
        elseif startswith(arg, "--files=")
            files = split(strip(split(arg, "=", limit=2)[2]), ','); i += 1
        elseif arg == "--outdir"
            i < length(args) || error("--outdir requires a path")
            outdir = args[i+1]; i += 2
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]; i += 1
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
        elseif arg == "--dry-run"
            dry_run = true; i += 1
        elseif arg in ("-h", "--help")
            _print_help()
            exit(0)
        else
            error("Unknown option: $arg  (use --help for usage)")
        end
    end

    return config_file, manifest_file, outdir, sep_grid, vsnr_grid, vfrac_grid, files, dry_run
end

function _print_help()
    println("""
    scan_resolved_lobe_thresholds.jl — Threshold scan for resolved-lobe audit

    Evaluates how different resolvability thresholds affect N_resolved against
    the chitosan benchmark manifest.

    Options:
      --config PATH             Fitting config TOML [$(DEFAULT_CONFIG)]
      --manifest PATH           Benchmark manifest [$(DEFAULT_MANIFEST)]
      --files F1,F2,...         Comma-separated .sxm files
                                Default: targets + controls + stress
      --outdir PATH             Output directory [$(DEFAULT_OUTDIR)]
      --sep-grid V1,V2,...      Comma-separated sep-threshold values
                                Default: $(join(string.(DEFAULT_SEP_GRID), ","))
      --valley-snr-grid V,...   Comma-separated valley-snr-threshold values
                                Default: $(join(string.(DEFAULT_VSNR_GRID), ","))
      --valley-frac-grid V,...  Comma-separated valley-frac-threshold values
                                Default: $(join(string.(DEFAULT_VFRAC_GRID), ","))
      --dry-run                 Print what would be done without executing
      -h, --help                This message

    Example:
      julia --project=. test/scan_resolved_lobe_thresholds.jl \\
          --sep-grid 1.0,1.5,2.0 \\
          --valley-snr-grid 1.0,3.0 \\
          --valley-frac-grid 0.1,0.2 \\
          --files 240817_019.sxm,240817_058.sxm
    """)
end

# ══════════════════════════════════════════════════════════════════════════════
# TSV reading helpers (mirrors grade_chitosan_benchmark.jl)
# ══════════════════════════════════════════════════════════════════════════════

function _read_tsv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("Empty TSV: $path")
    header = split(lines[1], '\t'; keepempty=true)
    rows = Vector{Dict{String,String}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, '\t'; keepempty=true)
        row = Dict{String,String}()
        for (j, name) in enumerate(header)
            row[name] = j <= length(vals) ? vals[j] : ""
        end
        push!(rows, row)
    end
    return header, rows
end

_as_bool(s::AbstractString) = lowercase(strip(s)) in ("true", "t", "1", "yes", "y")

function _parse_int_or_missing(s::AbstractString)
    t = strip(s)
    if isempty(t) || t in ("NA", "NaN", "missing")
        return missing
    end
    try
        return parse(Int, t)
    catch
        try
            return round(Int, parse(Float64, t))
        catch
            return missing
        end
    end
end

function _basename_file(s::AbstractString)
    b = basename(strip(s))
    if isempty(b)
        return b
    end
    m = match(r"(240817_\d+\.sxm)", b)
    return m === nothing ? b : m.captures[1]
end

# ══════════════════════════════════════════════════════════════════════════════
# Manifest helpers (mirrors grade_chitosan_benchmark.jl)
# ══════════════════════════════════════════════════════════════════════════════

function _load_quality_map(manifest::Dict{String,Any})
    """Build a map: basename -> (quality, exclude_primary, target_N_for_human_validation)"""
    files = get(manifest, "files", Dict{String,Any}())
    qmap = Dict{String,Tuple{String,Bool,Union{Int,Missing}}}()
    for (file_key, info) in files
        quality = String(get(info, "quality", "clean"))
        exclude_primary = Bool(get(info, "exclude_primary", false))
        target = get(info, "target_N_for_human_validation", missing)
        target = target === missing ? missing : Int(target)
        qmap[file_key] = (quality, exclude_primary, target)
    end
    return qmap
end

function _set_name(quality::String, exclude_primary::Bool, qmap, eval_cfg)
    stress_quality = Set(String.(get(eval_cfg, "stress_quality", ["poor_quality"])))
    primary_exclude_quality = Set(String.(get(eval_cfg, "primary_exclude_quality", ["poor_quality", "excluded"])))
    is_stress = quality in stress_quality
    is_primary = !(exclude_primary || quality in primary_exclude_quality)
    is_target = quality == "clean_target"
    return is_target ? "target" : (is_primary ? "primary" : (is_stress ? "stress" : "excluded"))
end

# ══════════════════════════════════════════════════════════════════════════════
# Grading a single audit TSV
# ══════════════════════════════════════════════════════════════════════════════

function _grade_audit_tsv(audit_tsv::String, manifest::Dict{String,Any};
                          default_target_n::Int=6)
    """
    Grade one audit result TSV against the benchmark manifest.
    Returns a named tuple with grade counts and failure lists.

    The audit TSV has one row per (file, N) pair with `is_selected` marking
    the effective-best row for each file.  We read `N_resolved` from selected
    rows and compare to the manifest's target_N_for_human_validation.
    """
    eval_cfg = get(manifest, "evaluation", Dict{String,Any}())
    primary_exclude_quality = Set(String.(get(eval_cfg, "primary_exclude_quality", ["poor_quality", "excluded"])))
    stress_quality = Set(String.(get(eval_cfg, "stress_quality", ["poor_quality"])))
    quality_map = _load_quality_map(manifest)

    # Read audit TSV
    !isfile(audit_tsv) && return nothing
    _, rows = _read_tsv(audit_tsv)
    isempty(rows) && return nothing

    # Build per-file selected row lookup
    selected_by_file = Dict{String,Dict{String,String}}()
    for row in rows
        haskey(row, "is_selected") || continue
        _as_bool(get(row, "is_selected", "false")) || continue
        file = _basename_file(get(row, "file", ""))
        isempty(file) && continue
        selected_by_file[file] = row
    end

    # Grade each file in the manifest that has a selected row
    primary_total = primary_ok = 0
    target_total = target_ok = 0
    stress_total = 0
    primary_failures = String[]
    target_failures = String[]
    missing_result = String[]

    for (file, row) in sort(collect(selected_by_file); by = p -> p.first)
        # Look up manifest info
        info = get(quality_map, file, ("clean", false, missing))
        quality, exclude_primary, manifest_target = info
        is_stress = quality in stress_quality
        is_primary = !(exclude_primary || quality in primary_exclude_quality)
        is_target = quality == "clean_target"
        target_n = manifest_target === missing ? default_target_n : manifest_target

        # Get observed N_resolved
        observed = _parse_int_or_missing(get(row, "N_resolved", ""))
        passed = observed !== missing && observed == target_n

        if observed === missing
            push!(missing_result, file)
        end
        if is_primary
            primary_total += 1
            if passed
                primary_ok += 1
            else
                push!(primary_failures, file)
            end
        end
        if is_target
            target_total += 1
            if passed
                target_ok += 1
            else
                push!(target_failures, file)
            end
        end
        if is_stress
            stress_total += 1
        end
    end

    return (primary_ok=primary_ok, primary_total=primary_total,
            target_ok=target_ok, target_total=target_total,
            stress_total=stress_total,
            primary_failures=join(primary_failures, ";"),
            target_failures=join(target_failures, ";"),
            missing_result=join(missing_result, ";"))
end

# ══════════════════════════════════════════════════════════════════════════════
# Subprocess invocation
# ══════════════════════════════════════════════════════════════════════════════

function _run_audit(config_file, files, out_tsv, sep, vsnr, vfrac)
    """Run audit_resolved_lobes.jl as a subprocess. Returns true on success."""
    files_str = join(files, ",")
    julia_exe = joinpath(Sys.BINDIR, "julia")
    # Use --project to pick up the project environment
    cmd = Cmd(`$julia_exe --project=$(PROJECT_DIR) $(AUDIT_SCRIPT)
              --config $(config_file)
              --files $(files_str)
              --out $(out_tsv)
              --sep-threshold $(sep)
              --valley-snr-threshold $(vsnr)
              --valley-frac-threshold $(vfrac)`)
    @printf("  Running: %s\n", cmd)
    try
        run(cmd)
        return true
    catch e
        @warn "Audit subprocess failed" sep vsnr vfrac exception=(e, catch_backtrace())
        return false
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Aggregate writer
# ══════════════════════════════════════════════════════════════════════════════

function _write_aggregate(outdir::String, grades::Vector)
    """Write aggregate TSV across all (completed) threshold combinations."""
    out_tsv = joinpath(outdir, "threshold_scores.tsv")
    header = [
        "sep_threshold", "valley_snr_threshold", "valley_frac_threshold",
        "primary_ok", "primary_total", "primary_pct",
        "target_ok", "target_total", "target_pct",
        "stress_total",
        "primary_failures", "target_failures", "missing_result",
        "audit_tsv"
    ]
    mkpath(outdir)
    open(out_tsv, "w") do io
        println(io, join(header, '\t'))
        for g in grades
            if g === nothing
                # Skipped/failed combo — write blank row with NA
                println(io, join([
                    "NA", "NA", "NA",
                    "NA", "NA", "NA",
                    "NA", "NA", "NA",
                    "NA",
                    "", "", "",
                    ""
                ], '\t'))
                continue
            end
            pct_prim = g.primary_total > 0 ?
                @sprintf("%.1f", 100g.primary_ok/g.primary_total) : "NA"
            pct_targ = g.target_total > 0 ?
                @sprintf("%.1f", 100g.target_ok/g.target_total) : "NA"
            println(io, join([
                @sprintf("%.4g", g.sep),
                @sprintf("%.4g", g.vsnr),
                @sprintf("%.4g", g.vfrac),
                g.primary_ok, g.primary_total, pct_prim,
                g.target_ok, g.target_total, pct_targ,
                g.stress_total,
                g.primary_failures,
                g.target_failures,
                g.missing_result,
                g.audit_tsv,
            ], '\t'))
        end
    end
    @printf("\nWrote aggregate: %s\n", out_tsv)
    return out_tsv
end

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

function main()
    args = _parse_cli(ARGS)
    config_file, manifest_file, outdir, sep_grid, vsnr_grid, vfrac_grid, files, dry_run = args

    # Load manifest
    manifest = TOML.parsefile(manifest_file)
    eval_cfg = get(manifest, "evaluation", Dict{String,Any}())
    default_target = 6  # chitosan calibration target

    total_combos = length(sep_grid) * length(vsnr_grid) * length(vfrac_grid)

    @printf("Resolved-lobe threshold scan\n")
    @printf("  config:     %s\n", config_file)
    @printf("  manifest:   %s\n", manifest_file)
    @printf("  files:      %d (targets=%d, controls=%d, stress=%d)\n",
            length(files),
            length(DEFAULT_TARGETS),
            length(DEFAULT_CONTROLS),
            length(DEFAULT_STRESS))
    @printf("  sep_grid:        %s\n", join([@sprintf("%.4g", v) for v in sep_grid], ", "))
    @printf("  valley_snr_grid: %s\n", join([@sprintf("%.4g", v) for v in vsnr_grid], ", "))
    @printf("  valley_frac_grid: %s\n", join([@sprintf("%.4g", v) for v in vfrac_grid], ", "))
    @printf("  total combos:    %d\n", total_combos)
    @printf("  outdir:     %s\n", outdir)
    @printf("  dry_run:    %s\n", dry_run)

    mkpath(outdir)

    # ── Build threshold combination grid ──
    combos = Tuple{Float64,Float64,Float64}[]
    for sep in sep_grid
        for vsnr in vsnr_grid
            for vfrac in vfrac_grid
                push!(combos, (sep, vsnr, vfrac))
            end
        end
    end

    # ── Iterate over threshold combinations ──
    grades = Vector{Any}(undef, total_combos)

    for (idx, (sep, vsnr, vfrac)) in enumerate(combos)
        @printf("\n[%d/%d] sep=%.4g  valley_snr=%.4g  valley_frac=%.4g\n",
                idx, total_combos, sep, vsnr, vfrac)

        # Build output path for this combo
        tag = @sprintf("sep%.4g_vsnr%.4g_vfrac%.4g", sep, vsnr, vfrac)
        audit_dir = joinpath(outdir, tag)
        audit_tsv = joinpath(audit_dir, "resolved_lobes.tsv")

        if dry_run
            @printf("  Would write: %s\n", audit_tsv)
            grades[idx] = nothing
            continue
        end

        # Run audit
        mkpath(audit_dir)
        success = _run_audit(config_file, files, audit_tsv, sep, vsnr, vfrac)
        if !success
            @warn "Audit failed for combo, recording blank grade" sep vsnr vfrac
            grades[idx] = nothing
            continue
        end

        # Grade
        grade = _grade_audit_tsv(audit_tsv, manifest; default_target_n=default_target)
        if grade === nothing
            @warn "Grading returned no data for" audit_tsv
            grades[idx] = nothing
            continue
        end

        # Attach threshold metadata
        grades[idx] = (sep=sep, vsnr=vsnr, vfrac=vfrac,
                       audit_tsv=audit_tsv,
                       primary_ok=grade.primary_ok,
                       primary_total=grade.primary_total,
                       target_ok=grade.target_ok,
                       target_total=grade.target_total,
                       stress_total=grade.stress_total,
                       primary_failures=grade.primary_failures,
                       target_failures=grade.target_failures,
                       missing_result=grade.missing_result)

        @printf("  primary: %d/%d  target: %d/%d  stress: %d\n",
                grade.primary_ok, grade.primary_total,
                grade.target_ok, grade.target_total,
                grade.stress_total)
        if !isempty(grade.primary_failures)
            @printf("  primary failures: %s\n", grade.primary_failures)
        end
        if !isempty(grade.target_failures)
            @printf("  target failures: %s\n", grade.target_failures)
        end
    end

    # ── Write aggregate ──
    if !dry_run
        _write_aggregate(outdir, grades)
    end

    if dry_run
        @printf("\nDry run complete.  %d combos would be processed (%d files each).\n",
                total_combos, length(files))
    end
end

main()

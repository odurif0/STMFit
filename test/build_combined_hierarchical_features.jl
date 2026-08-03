#!/usr/bin/env julia

# Build a combined hierarchical feature TSV that includes backward-view
# descriptors (bwd_neg_com_t, bwd_neg_diag45) and split-width skew
# (split_log_skew) alongside the existing base features.
#
# This exists because the T5 leave-date-out evaluator
# (evaluate_hierarchical_unit_assignment.jl) reads a single feature TSV and
# its _default_views only activates backward views when bwd_neg_* columns are
# present.  The original T5 run omitted them, so forward_backward_agreement
# was NaN for every row and no challenger was eligible at T7.
#
# The bwd_neg_* computation is identical to build_hierarchical_unit_predictions.jl
# (_negative_moment + _merge_patches!).  This script is label-free: it reads no
# truth, sequence, expected-N, or benchmark data.
#
# Usage:
#   julia --project=. test/build_combined_hierarchical_features.jl \
#     --features results/.../lobe_features_selectedN_primary_local.tsv \
#     --patches  results/.../lobe_patches_selectedN_primary_17x17_bwd.tsv \
#     [--split   results/.../lobe_features_selectedN_primary_split.tsv] \
#     --out      results/.../combined_features.tsv

include(joinpath(@__DIR__, "lib", "hierarchical_unit_assignment.jl"))
using .HierarchicalUnitAssignment: load_records, LobeRecord, _read_tsv,
      _parse_f, _parse_lobe, _key, is_forbidden_feature, NON_FEATURE_COLUMNS

using Printf

# --------------------------------------------------------------------------
# Lenient loader: drops forbidden columns from the INPUT instead of erroring.
# The combined OUTPUT is still firewall-checked in write_combined.
# --------------------------------------------------------------------------

function _load_records_lenient(path::AbstractString)
    isfile(path) || error("Feature TSV not found: $path")
    header, rows = _read_tsv(path)
    isempty(rows) && error("No rows in $path")
    for col in ("file", "lobe", "amplitude")
        col in header || error("Feature TSV $path missing required column: $col")
    end
    dropped = String[]
    kept_cols = filter(header) do col
        if col in NON_FEATURE_COLUMNS
            return true  # always keep structural columns
        elseif is_forbidden_feature(col)
            push!(dropped, col)
            return false
        else
            return true
        end
    end
    isempty(dropped) || @warn("dropped forbidden input columns from $path" *
        ":\n  $(join(dropped, ", "))")
    records = LobeRecord[]
    seen = Set{Tuple{String,Int}}()
    for row in rows
        # Preserve the ORIGINAL file path (with date-folder prefix like
        # 20240307_LHe_Cu100/240307_015.sxm).  The evaluator's
        # parse_leading_date needs the full path, not just the basename.
        file = strip(row["file"])
        isempty(file) && error("Feature TSV $path has empty file value")
        lobe = _parse_lobe(row, "Feature TSV $path")
        key = _key(file, lobe)  # _key uses basename internally for matching
        key in seen && error("Feature TSV $path has duplicate row for $(key[1]) lobe $(key[2])")
        push!(seen, key)
        amp = _parse_f(row["amplitude"])
        isfinite(amp) || error("Feature TSV $path has non-finite amplitude for $file lobe $lobe")
        feat = Dict{String,Float64}()
        for col in kept_cols
            col in NON_FEATURE_COLUMNS && continue
            parsed = tryparse(Float64, strip(row[col]))
            parsed === nothing && continue
            feat[col] = parsed
        end
        if !haskey(feat, "integrated") &&
           all(haskey(feat, k) for k in ("amplitude", "sigma_parallel_nm", "sigma_perp_nm"))
            feat["integrated"] = feat["amplitude"] * feat["sigma_parallel_nm"] * feat["sigma_perp_nm"]
        end
        push!(records, LobeRecord(file, lobe, amp, feat))
    end
    sort!(records, by = r -> (r.file, r.lobe))
    return records
end

# --------------------------------------------------------------------------
# Backward-patch negative moment (exact port of build_hierarchical_unit_predictions.jl)
# --------------------------------------------------------------------------

function _patch_columns(header::Vector{String}, prefix::String)
    cols = [c for c in header if startswith(c, prefix)]
    isempty(cols) && return String[]
    side = round(Int, sqrt(length(cols)))
    side * side == length(cols) ||
        error("Patch columns for prefix $prefix are not a square grid: $(length(cols))")
    return cols
end

function _parse_val(s)
    t = strip(String(s))
    (isempty(t) || t in ("NA", "NaN", "nan", "?")) && return NaN
    v = tryparse(Float64, t)
    return v === nothing ? NaN : v
end

function _negative_moment(vals::Vector{Float64}, coords::Vector{Float64})
    weights = max.(-vals, 0.0)
    total = sum(weights)
    total <= eps(Float64) && return (com_t=NaN, diag45=NaN, diag135=NaN)
    com_t = 0.0; diag45 = 0.0; diag135 = 0.0
    idx = 1
    for u in coords, t in coords
        w = weights[idx]
        com_t += w * t
        diag45 += w * ((t + u) / sqrt(2.0))
        diag135 += w * ((t - u) / sqrt(2.0))
        idx += 1
    end
    return (com_t=com_t / total, diag45=diag45 / total, diag135=diag135 / total)
end

function _merge_patches!(records::Vector{LobeRecord}, path::AbstractString)
    isempty(path) && return nothing
    isfile(path) || error("Patch TSV not found: $path")
    header, rows = _read_tsv(path)
    cols = _patch_columns(String.(header), "bwd_res_p")
    isempty(cols) && error("Patch TSV $path has no bwd_res_pNNN columns")
    side = round(Int, sqrt(length(cols)))
    coords = side == 1 ? [0.0] : collect(range(-1.0, 1.0; length=side))
    by_key = Dict{Tuple{String,Int},NamedTuple{(:com_t,:diag45,:diag135),Tuple{Float64,Float64,Float64}}}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Patch TSV $path has empty file value")
        haskey(row, "lobe") || error("Patch TSV $path missing lobe column")
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe === nothing && error("Patch TSV $path invalid lobe: $(row["lobe"])")
        lobe >= 1 || error("Patch TSV $path lobe must be >= 1, got $lobe")
        key = (file, lobe)
        haskey(by_key, key) && error("Patch TSV $path duplicate row for $file lobe $lobe")
        vals = Float64[_parse_val(row[c]) for c in cols]
        by_key[key] = _negative_moment(vals, coords)
    end
    n_finite = 0
    for rec in records
        d = get(by_key, _key(rec.file, rec.lobe),
                (com_t=NaN, diag45=NaN, diag135=NaN))
        rec.features["bwd_neg_com_t"] = d.com_t
        rec.features["bwd_neg_diag45"] = d.diag45
        rec.features["bwd_neg_diag135"] = d.diag135
        isfinite(d.com_t) && (n_finite += 1)
    end
    @printf("merged backward patches: %d/%d finite bwd_neg_com_t\n", n_finite, length(records))
    return nothing
end

function _merge_split!(records::Vector{LobeRecord}, path::AbstractString)
    isempty(path) && return nothing
    isfile(path) || error("Split TSV not found: $path")
    header, rows = _read_tsv(path)
    "skew_ratio" in header || error("Split TSV $path missing skew_ratio column")
    by_key = Dict{Tuple{String,Int},Float64}()
    for row in rows
        file = basename(strip(row["file"]))
        isempty(file) && error("Split TSV $path has empty file value")
        haskey(row, "lobe") || error("Split TSV $path missing lobe column")
        lobe = tryparse(Int, strip(row["lobe"]))
        lobe === nothing && error("Split TSV $path invalid lobe: $(row["lobe"])")
        lobe >= 1 || error("Split TSV $path lobe must be >= 1, got $lobe")
        key = (file, lobe)
        haskey(by_key, key) && error("Split TSV $path duplicate row for $file lobe $lobe")
        skew = _parse_val(row["skew_ratio"])
        by_key[key] = isfinite(skew) && skew > 0 ? log(skew) : NaN
    end
    n_finite = 0
    for rec in records
        ls = get(by_key, _key(rec.file, rec.lobe), NaN)
        rec.features["split_log_skew"] = ls
        isfinite(ls) && (n_finite += 1)
    end
    @printf("merged split skew: %d/%d finite split_log_skew\n", n_finite, length(records))
    return nothing
end

# --------------------------------------------------------------------------
# Write combined TSV
# --------------------------------------------------------------------------

function _fmt(x::Float64)
    isfinite(x) ? @sprintf("%.17g", x) : "NA"
end

function write_combined(records::Vector{LobeRecord}, out_path::AbstractString)
    # Collect all feature names (sorted for determinism), excluding forbidden
    # and non-feature columns.
    all_feats = Set{String}()
    for rec in records
        for k in keys(rec.features)
            push!(all_feats, k)
        end
    end
    feat_cols = sort!(collect(all_feats))
    # Final firewall check on the combined feature set.
    for c in feat_cols
        is_forbidden_feature(c) &&
            error("Combined feature column is forbidden: $c")
    end
    # Exclude "amplitude" from feature columns — it is already written as a
    # structural column (rec.amplitude) to avoid a duplicate-column TSV.
    filter!(!=("amplitude"), feat_cols)
    cols = ["file", "lobe", "amplitude", feat_cols...]
    parent = dirname(abspath(out_path))
    mkpath(parent)
    tmp, io = mktemp(parent; cleanup=false)
    try
        println(io, join(cols, '\t'))
        for rec in records
            vals = String[]
            push!(vals, rec.file)
            push!(vals, string(rec.lobe))
            push!(vals, _fmt(rec.amplitude))
            for fn in feat_cols
                push!(vals, _fmt(get(rec.features, fn, NaN)))
            end
            println(io, join(vals, '\t'))
        end
        flush(io); close(io)
        mv(tmp, abspath(out_path); force=true)
    catch
        close(io)
        rm(tmp; force=true)
        rethrow()
    end
    @printf("wrote %s: %d rows, %d cols (%d features)\n", out_path, length(records), length(cols), length(feat_cols))
    return out_path
end

# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

function _show_help()
    println("""
Usage: julia --project=. test/build_combined_hierarchical_features.jl [options]

Build a combined hierarchical feature TSV with backward-view descriptors.

Options:
  --features PATH   Base label-free feature TSV (required).
  --patches PATH    Optional backward patch TSV; adds bwd_neg_com_t/diag45/diag135.
  --split PATH      Optional split-width TSV; adds split_log_skew.
  --out PATH        Output combined feature TSV (required).

Label-free: reads no truth, sequence, expected-N, or benchmark data.
""")
end

function main(args = ARGS)
    vals = Dict("--features"=>"", "--patches"=>"", "--split"=>"", "--out"=>"")
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            _show_help(); return nothing
        elseif haskey(vals, arg)
            vals[arg] = args[i+1]; i += 2
        elseif startswith(arg, "--") && occursin('=', arg)
            key, value = split(arg, '='; limit=2)
            haskey(vals, key) || error("unknown argument: $arg")
            vals[key] = value; i += 1
        else
            error("unknown argument: $arg")
        end
    end
    isempty(vals["--features"]) && (_show_help(); error("--features is required"))
    isempty(vals["--out"]) && (_show_help(); error("--out is required"))

    records = _load_records_lenient(vals["--features"])
    @printf("loaded %d base records from %s\n", length(records), vals["--features"])
    _merge_patches!(records, vals["--patches"])
    _merge_split!(records, vals["--split"])
    write_combined(records, vals["--out"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

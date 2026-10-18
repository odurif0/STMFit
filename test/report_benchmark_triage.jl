#!/usr/bin/env julia
# Kept/discarded triage of raw 6-mer scans for the benchmark (grading side).
#
# Usage:
#   julia -t 4 --project=. test/report_benchmark_triage.jl --catalog CATALOG_TSV \
#       [--catalog CATALOG_TSV ...] --count-config config/chitosan.toml --sheets NEW_DIR
#
# Review sheets for the visual triage of raw scans. The catalogs
# (catalog_sxm_candidates.jl) give every scan once by sha256, the first
# catalog winning; the usable scans are drawn in an 8 x 8 nm window centred on
# the molecule (same scale on every tile, 2 nm bar, unscanned pixels blank),
# tagged [B] in the benchmark manifest, [R] reviewed before, [N] never
# reviewed. Decisions go to a triage table; test/report_data_sets.jl turns the
# tables into data/sets/ and data/INDEX.tsv. Nothing here feeds the method.
include(joinpath(@__DIR__, "report_benchmark_candidates.jl"))

const CROP_NM = 8.0

function triage_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); cats = String[]; i = 1
    while i <= length(args)
        k = args[i]
        k in ("--catalog", "--count-config", "--sheets") && i < length(args) ||
            error("Unknown option or missing value: $k")
        k == "--catalog" ? push!(cats, args[i+1]) : (o[k] = args[i+1]); i += 2
    end
    (isempty(cats) || !haskey(o, "--count-config") || !haskey(o, "--sheets")) && error("Required: --catalog --count-config --sheets")
    (ispath(o["--sheets"]) || islink(o["--sheets"])) && error("Output exists: $(o["--sheets"])")
    return o, cats
end

"Every scan once (by sha256) over the catalogs, the first catalog winning."
function scan_universe(cats)
    seen = Set{String}(); rows = Dict{String,String}[]
    for c in cats, r in grading_rows(c)
        (isempty(r["sha256"]) || r["sha256"] in seen) && continue
        push!(seen, r["sha256"]); push!(rows, r)
    end
    return rows
end

function review_tile(r, v, tag)
    z = copy(v.z); z[.!v.observed] .= NaN
    w = Float64.(v.mask .& v.observed); s = sum(w); h = CROP_NM / 2
    cx = s > 0 ? sum(w .* v.xs') / s : (first(v.xs) + last(v.xs)) / 2
    cy = s > 0 ? sum(w .* v.ys) / s : (first(v.ys) + last(v.ys)) / 2
    ix = findall(x -> abs(x - cx) <= h, v.xs); iy = findall(y -> abs(y - cy) <= h, v.ys)
    vals = filter(isfinite, vec(z[iy, ix]))
    lo, hi = isempty(vals) ? (0.0, 1.0) : (quantile(vals, 0.01), quantile(vals, 0.995))
    hi > lo || (hi = lo + 1e-6)
    title = "$tag $(r["file"])  $(r["range_nm"]) nm" * (isempty(r["flags"]) ? "" : "\n$(r["flags"])")
    fig = heatmap(v.xs[ix], v.ys[iy], z[iy, ix]; c=:inferno, clims=(lo, hi), aspect_ratio=:equal, colorbar=false,
                  axis=false, grid=false, xlims=(cx - h, cx + h), ylims=(cy - h, cy + h), title=title, titlefontsize=8)
    contour!(fig, v.xs[ix], v.ys[iy], w[iy, ix]; levels=[0.5], color=:cyan, linewidth=0.6, alpha=0.5, colorbar=false)
    plot!(fig, [cx - h + 0.4, cx - h + 2.4], [cy - h + 0.4, cy - h + 0.4]; color=:cyan, linewidth=2, label="")
    return fig
end

"Review sheets, 20 tiles per page, in the given order; `tag(r)` prefixes each title."
function review_sheets(rows, views, outdir, title; tag=r -> "")
    mkpath(outdir); per = 20; pages = cld(length(rows), per)
    for p in 1:pages
        part = rows[(p-1)*per+1:min(p * per, length(rows))]
        panels = [review_tile(r, views[r["sha256"]], tag(r)) for r in part]
        fig = plot(panels...; layout=(cld(length(panels), 5), 5), size=(1600, 340 * cld(length(panels), 5)),
                   plot_title="$title — page $p/$pages", plot_titlefontsize=12)
        savefig(fig, joinpath(outdir, @sprintf("page_%03d.png", p)))
    end
    return pages
end

function scan_views(rows, pcfg)
    out = Vector{Any}(undef, length(rows))
    Threads.@threads for k in eachindex(rows)
        out[k] = try last(inspect_scan(rows[k]["path"], pcfg, Inf)) catch; nothing end
    end
    return Dict(rows[k]["sha256"] => out[k] for k in eachindex(rows) if out[k] !== nothing)
end

function triage_main(args=ARGS)
    parsed = triage_options(args)
    parsed === nothing && return println("report_benchmark_triage.jl --catalog TSV [--catalog TSV] --count-config TOML --sheets NEW_DIR")
    o, cats = parsed
    pcfg = count_preprocessing(o["--count-config"])
    rows = sort(scan_universe(cats); by=r -> (r["session"], r["file"]))
    usable = [r for r in rows if r["status"] == "usable"]
    groups = Dict(r["sha256"] => r["group"] for r in candidate_groups(usable,
                  TOML.parsefile(DEFAULT_MANIFEST), grading_rows(DEFAULT_REVIEW)))
    tags = Dict("in_benchmark" => "[B]", "reviewed_before" => "[R]", "never_reviewed" => "[N]")
    pages = review_sheets(usable, scan_views(usable, pcfg), o["--sheets"], "Review"; tag=r -> tags[groups[r["sha256"]]])
    println("Review sheets: $(length(usable)) usable scans, $pages pages")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && triage_main()

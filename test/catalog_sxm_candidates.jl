#!/usr/bin/env julia
# Catalog and technical triage of raw STM scans (benchmark candidates, new data).
#
# Usage:
#   julia -t 4 --project=. test/catalog_sxm_candidates.jl --roots DIR[,DIR...] \
#       --consensus-config config/molecule_consensus.toml --count-config config/chitosan.toml \
#       --outdir NEW_DIR [--exclude TSV] [--copy-to CENTRAL_DIR]
#
# No fitting, no count and no label: only headers, pixels and the molecule ROI.
# Every .sxm under the roots gets one row in catalog.tsv:
#   status   usable | unreadable | no_z_fwd_bwd | overview (frame larger than
#            selection.collect_max_range_nm) | object_crosses_frame (the
#            observed ROI touches two opposite borders: a zoom on part of an
#            object or a substrate feature, no whole molecule) | excluded_by_review (listed in
#            --exclude, columns file and reason: a human visual decision) |
#            duplicate (same bytes, or same acquisition: date, time and
#            original SCAN_FILE, as an earlier file or a central scan)
#   flags    (usable scans, for visual review, never an exclusion):
#            partial_scan (scan stopped early: < 99% observed Z pixels),
#            roi_touches_unscanned (molecule may be cut by the end of the scan),
#            roi_touches_border, several_objects (second ROI component >= 30%
#            of the largest), fwd_bwd_mismatch (forward/backward correlation < 0.5)
# --copy-to copies every usable scan once into one flat central folder (basename,
# or session_basename if the basename is taken) and keeps its MANIFEST.tsv
# (provenance only). An existing central folder must hold that manifest: its
# scans count as earlier copies, so new roots only add scans it does not have.
# contact_sheets/ shows the usable scans (fused image, ROI contour), 24 per page.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using GaussianFit2D, STMSXMIO, SHA, Statistics, Printf, TOML, Plots, Dates
const G = GaussianFit2D

function catalog_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--roots", "--consensus-config", "--count-config", "--outdir", "--copy-to", "--exclude") && i < length(args) ||
            error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--roots", "--consensus-config", "--count-config", "--outdir")) ||
        error("Required: --roots --consensus-config --count-config --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    haskey(o, "--copy-to") && ispath(o["--copy-to"]) && !isfile(joinpath(o["--copy-to"], "MANIFEST.tsv")) &&
        error("Central folder exists without MANIFEST.tsv: $(o["--copy-to"])")
    return o
end

hv(h, k) = strip(get(h, k, ""))
nums(s) = [parse(Float64, t) for t in split(s) if tryparse(Float64, t) !== nothing]

"Connected-component areas (4-neighbour) of a boolean image, largest first."
function component_areas(m::AbstractMatrix{Bool})
    seen = falses(size(m)); areas = Int[]
    for start in CartesianIndices(m)
        (m[start] && !seen[start]) || continue
        stack = [start]; seen[start] = true; a = 0
        while !isempty(stack)
            c = pop!(stack); a += 1
            for d in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
                n = c + d
                checkbounds(Bool, m, n) && m[n] && !seen[n] && (seen[n] = true; push!(stack, n))
            end
        end
        push!(areas, a)
    end
    return sort(areas; rev=true)
end

function inspect_scan(path, pcfg, max_range)
    row = Dict{String,String}("path" => path, "file" => basename(path), "size_bytes" => string(filesize(path)),
                              "sha256" => bytes2hex(open(sha256, path)), "status" => "usable", "flags" => "")
    img = try read_sxm(path) catch e; row["status"] = "unreadable"; row["detail"] = sprint(showerror, e)[1:min(end, 120)]; return row, nothing end
    h = img.header
    scanfile = replace(hv(h, "SCAN_FILE"), '\\' => '/')
    row["rec_date"] = hv(h, "REC_DATE"); row["rec_time"] = hv(h, "REC_TIME")
    row["scan_file"] = scanfile; row["session"] = basename(dirname(path))
    row["range_nm"] = @sprintf("%.2fx%.2f", img.range_nm...); row["pixels"] = "$(img.width)x$(img.height)"
    row["angle_deg"] = hv(h, "SCAN_ANGLE"); row["bias_v"] = hv(h, "BIAS")
    z = [c for c in img.channels if lowercase(c.name) == "z"]
    dirs = Set(lowercase(c.direction) for c in z)
    ("fwd" in dirs && "bwd" in dirs) || (row["status"] = "no_z_fwd_bwd"; return row, nothing)
    finite = minimum(count(isfinite, c.data) / length(c.data) for c in z)
    row["finite_fraction"] = @sprintf("%.3f", finite)
    maximum(img.range_nm) > max_range && (row["status"] = "overview"; return row, nothing)
    flags = String[]
    finite < 0.99 && push!(flags, "partial_scan")
    xs, ys, z0, mask, x, y, zr, noise = G._fused_roi_data(img, pcfg)
    st = max(1, pcfg.stride)
    observed = reduce(.&, [isfinite.(c.data[1:st:end, 1:st:end]) for c in z])
    size(observed) == size(mask) || (observed = trues(size(mask)))
    if finite < 1
        near = copy(.!observed)
        for idx in findall(.!observed), d in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
            n = idx + d; checkbounds(Bool, near, n) && (near[n] = true)
        end
        any(mask .& near) && push!(flags, "roi_touches_unscanned")
    end
    any(mask[1, :]) || any(mask[end, :]) || any(mask[:, 1]) || any(mask[:, end]) ? push!(flags, "roi_touches_border") : nothing
    seen = mask .& observed
    crosses = (any(seen[1, :]) && any(seen[end, :])) || (any(seen[:, 1]) && any(seen[:, end]))
    crosses && (row["status"] = "object_crosses_frame")
    signal = z0 .- minimum(z0)
    thr = G._adaptive_roi_threshold(signal, noise, pcfg.roi_threshold_fraction, 2.5)
    areas = component_areas(signal .>= thr)
    length(areas) >= 2 && areas[2] >= 0.3 * areas[1] && push!(flags, "several_objects")
    ad = G.compute_image_artifact_diagnostics(img, pcfg)
    row["fwd_bwd_corr"] = @sprintf("%.3f", ad.fwd_bwd_corr)
    isfinite(ad.fwd_bwd_corr) && ad.fwd_bwd_corr < 0.5 && push!(flags, "fwd_bwd_mismatch")
    ac = G._weighted_roi_axis(x, y, zr)
    row["roi_length_nm"] = @sprintf("%.2f", ac.tmax - ac.tmin)
    row["flags"] = join(flags, ",")
    return row, (xs=xs, ys=ys, z=z0, mask=mask)
end

function contact_sheets(rows, views, outdir; title="Usable scans")
    usable = [r for r in rows if r["status"] == "usable"]
    sort!(usable; by=r -> (get(r, "session", ""), r["file"]))
    mkpath(outdir); per = 24; pages = cld(length(usable), per)
    for p in 1:pages
        part = usable[(p-1)*per+1:min(p * per, length(usable))]
        panels = map(part) do r
            v = views[r["path"]]
            fig = heatmap(v.xs, v.ys, v.z; c=:inferno, aspect_ratio=:equal, colorbar=false, axis=false,
                          title="$(r["file"])  $(r["range_nm"]) nm" * (isempty(r["flags"]) ? "" : "\n$(r["flags"])"),
                          titlefontsize=7)
            contour!(fig, v.xs, v.ys, Float64.(v.mask); levels=[0.5], color=:cyan, linewidth=1, colorbar=false)
            fig
        end
        fig = plot(panels...; layout=(cld(length(panels), 6), 6), size=(1800, 330 * cld(length(panels), 6)),
                   plot_title="$title — page $p/$pages", plot_titlefontsize=12)
        savefig(fig, joinpath(outdir, @sprintf("page_%03d.png", p)))
    end
    return pages
end

"Preprocessing of the counting config (as the fit sees the image)."
function count_preprocessing(count_config)
    pre = get(TOML.parsefile(count_config), "preprocessing", Dict{String,Any}())
    return G.PatternConfig(filepath="", channel="Z", direction="fwd", stride=get(pre, "stride", 1),
        flatten=get(pre, "flatten", "plane+rows"), smooth_radius_px=get(pre, "smooth_radius_px", 1), no_plot=true)
end

function main(args=ARGS)
    o = catalog_options(args)
    o === nothing && return println("catalog_sxm_candidates.jl --roots DIR[,DIR] --consensus-config TOML --count-config TOML --outdir NEW_DIR [--exclude TSV] [--copy-to CENTRAL_DIR]")
    max_range = Float64(TOML.parsefile(o["--consensus-config"])["selection"]["collect_max_range_nm"])
    pcfg = count_preprocessing(o["--count-config"])
    paths = String[]
    for root in split(o["--roots"], ","), (d, _, names) in walkdir(root), n in names
        endswith(lowercase(n), ".sxm") && push!(paths, joinpath(d, n))
    end
    sort!(paths)
    results = Vector{Any}(undef, length(paths))
    Threads.@threads for i in eachindex(paths)
        results[i] = try inspect_scan(paths[i], pcfg, max_range) catch e
            (Dict("path" => paths[i], "file" => basename(paths[i]), "status" => "unreadable",
                  "detail" => sprint(showerror, e)[1:min(end, 120)], "flags" => ""), nothing)
        end
    end
    rows = [first(r) for r in results]; views = Dict(paths[i] => last(results[i]) for i in eachindex(paths) if last(results[i]) !== nothing)
    if haskey(o, "--exclude")
        excluded = Dict(r["file"] => r["reason"] for r in last(read_table(o["--exclude"])))
        for r in rows
            haskey(excluded, r["file"]) && r["status"] == "usable" &&
                (r["status"] = "excluded_by_review"; r["detail"] = excluded[r["file"]])
        end
    end
    # duplicates: same bytes, or same acquisition (date, time, original file);
    # a scan already in the central folder comes first, then the first path
    acq(r) = join((r["rec_date"], r["rec_time"], basename(r["scan_file"])), "|")
    first_by = Dict{String,String}(); central = Dict{String,String}[]
    dest = get(o, "--copy-to", "")
    if isfile(joinpath(dest, "MANIFEST.tsv"))
        central = last(read_table(joinpath(dest, "MANIFEST.tsv")))
        for c in central
            first_by[c["sha256"]] = joinpath(dest, c["central_name"])
            isempty(c["scan_file"]) || (first_by[acq(c)] = joinpath(dest, c["central_name"]))
        end
    end
    for r in rows
        r["status"] in ("unreadable",) && continue
        keys_ = [r["sha256"]]
        haskey(r, "rec_date") && !isempty(get(r, "scan_file", "")) && push!(keys_, acq(r))
        dup = [first_by[k] for k in keys_ if haskey(first_by, k)]
        if !isempty(dup)
            r["status"] = "duplicate"; r["duplicate_of"] = first(dup)
        else
            for k in keys_; first_by[k] = r["path"]; end
        end
    end
    out = o["--outdir"]; mkpath(out)
    mcols = ["central_name", "sha256", "size_bytes", "session", "rec_date", "rec_time", "range_nm", "pixels",
             "flags", "scan_file", "source_path", "added"]
    if !isempty(dest)
        isempty(central) && !any(r -> r["status"] == "usable", rows) && error("No usable scan for a new central folder")
        mkpath(dest)
        taken = Set(c["central_name"] for c in central); added = Dict{String,String}[]
        names = Dict{String,Int}()
        for r in rows; r["status"] == "usable" && (names[r["file"]] = get(names, r["file"], 0) + 1); end
        for r in rows
            r["status"] == "usable" || continue
            name = names[r["file"]] == 1 && !(r["file"] in taken) ? r["file"] : string(get(r, "session", "unknown"), "_", r["file"])
            name in taken && error("Central name collision: $name")
            cp(r["path"], joinpath(dest, name)); push!(taken, name); r["central_name"] = name
            push!(added, merge(r, Dict("source_path" => r["path"], "added" => string(Dates.today()))))
        end
        manifest = vcat(central, [Dict(c => replace(get(r, c, ""), '\t' => ' ', '\n' => ' ') for c in mcols) for r in added])
        tmp = joinpath(dest, "MANIFEST.tsv.new")
        write_table(tmp, mcols, sort(manifest; by=c -> c["central_name"]))
        mv(tmp, joinpath(dest, "MANIFEST.tsv"); force=true)
        println("Central folder $dest: $(length(added)) scans added, $(length(manifest)) in total")
    end
    cols = ["file", "status", "flags", "duplicate_of", "session", "rec_date", "rec_time", "range_nm", "pixels",
            "angle_deg", "bias_v", "finite_fraction", "fwd_bwd_corr", "roi_length_nm", "central_name",
            "size_bytes", "sha256", "scan_file", "path", "detail"]
    write_table(joinpath(out, "catalog.tsv"), cols,
        [Dict(c => replace(get(r, c, ""), '\t' => ' ', '\n' => ' ') for c in cols) for r in rows])
    pages = contact_sheets(rows, views, joinpath(out, "contact_sheets"))
    tally = Dict{String,Int}(); for r in rows; tally[r["status"]] = get(tally, r["status"], 0) + 1; end
    fl = Dict{String,Int}(); for r in rows, f in split(r["flags"], ",", keepempty=false); r["status"] == "usable" && (fl[f] = get(fl, f, 0) + 1); end
    println("Scans: $(length(rows)); ", join(["$k $v" for (k, v) in sort(collect(tally))], ", "))
    println("Flags on usable scans: ", join(["$k $v" for (k, v) in sort(collect(fl))], ", "), "; contact sheets: $pages pages")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

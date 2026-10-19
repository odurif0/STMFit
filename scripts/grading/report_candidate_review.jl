#!/usr/bin/env julia
# Local review page for the 6-mer benchmark candidates (grading side, human decision).
#
# Usage:
#   julia -t 4 --project=. scripts/grading/report_candidate_review.jl --data data \
#       --triage benchmarks/chitosan_6mer_cu100_triage.tsv --outdir NEW_DIR [--suggestions TSV]
#
# Takes the scans of set 6mer_Cu100_candidates and groups repeat scans of one
# molecule: same session, consecutive acquisitions at most collect_max_gap_min
# apart, frame centres within collect_max_offset_nm (config/molecule_consensus.toml).
# Each scan is drawn twice (whole frame; 8 x 8 nm window on the molecule, 2 nm
# bar) in img/, and sheets/ holds contact sheets in group order. index.html is
# a local page: one card per scan with the suggestion (--suggestions: columns
# file, suggestion benchmark|ecarte, motif, note; default benchmark) and the
# first-pass note, buttons to decide, and an export of the decisions as TSV
# (file, decision, motif, note, suggestion). One benchmark scan per session is
# shown for reference. Nothing leaves the machine; decisions stay in the
# browser (localStorage) until exported.
include(joinpath(@__DIR__, "report_benchmark_triage.jl"))

const CANDIDATE_SET = "6mer_Cu100_candidates"
const BENCHMARK_SET = "6mer_Cu100_benchmark"
const MOTIFS = ["pas_chitosan", "pas_6mere_court", "pas_6mere_long", "plusieurs_objets", "molecule_coupee",
                "artefact_pointe", "pas_de_molecule", "incertain"]

function review_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--data", "--triage", "--outdir", "--suggestions") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--data", "--triage", "--outdir")) || error("Required: --data --triage --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    return o
end

acquired(img) = DateTime(strip(get(img.header, "REC_DATE", "01.01.2000")) * " " * strip(get(img.header, "REC_TIME", "00:00:00")),
                         dateformat"dd.mm.yyyy HH:MM:SS")
centre(img) = (img.offset_nm[1], img.offset_nm[2])

"Group ids: a new group when the session changes or the next scan is too late or too far."
function molecule_groups(scans, gap_min, dmax)
    ids = Int[]; g = 0; prev = nothing
    for s in scans
        if prev === nothing || s.session != prev.session ||
           Dates.value(s.time - prev.time) / 60000 > gap_min || hypot((s.centre .- prev.centre)...) > dmax
            g += 1
        end
        push!(ids, g); prev = s
    end
    return ids
end

function scan_figure(row, v, path)
    z = copy(v.z); z[.!v.observed] .= NaN
    vals = filter(isfinite, vec(z))
    lo, hi = isempty(vals) ? (0.0, 1.0) : (quantile(vals, 0.01), quantile(vals, 0.995))
    full = heatmap(v.xs, v.ys, z; c=:inferno, clims=(lo, max(hi, lo + 1e-6)), aspect_ratio=:equal, colorbar=false,
                   axis=false, grid=false, title="cadre complet $(row["range_nm"]) nm", titlefontsize=9)
    contour!(full, v.xs, v.ys, Float64.(v.mask .& v.observed); levels=[0.5], color=:cyan, linewidth=0.6, alpha=0.6, colorbar=false)
    crop = review_tile(merge(row, Dict("flags" => "", "range_nm" => "fenêtre 8 x 8")), v, "")
    savefig(plot(full, crop; layout=(1, 2), size=(980, 470), margin=2Plots.mm), path)
end

js(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ", "<" => "\\u003c") * "\""

function review_main(args=ARGS)
    o = review_options(args)
    o === nothing && return println("report_candidate_review.jl --data DIR --triage TSV --outdir NEW_DIR [--suggestions TSV]")
    root = dirname(dirname(@__DIR__))
    sel = TOML.parsefile(joinpath(root, "config", "molecule_consensus.toml"))["selection"]
    pcfg = count_preprocessing(joinpath(root, "config", "chitosan.toml"))
    _, tri = read_table(o["--triage"])
    sug = Dict{String,Dict{String,String}}()
    haskey(o, "--suggestions") && (sug = Dict(r["file"] => r for r in last(read_table(o["--suggestions"]))))
    raw = joinpath(o["--data"], "raw")
    cands = [t for t in tri if t["set"] == CANDIDATE_SET]
    refs = Dict{String,Dict{String,String}}()
    for t in sort([t for t in tri if t["set"] == BENCHMARK_SET]; by=t -> t["file"])
        haskey(refs, t["session"]) || (refs[t["session"]] = t)
    end
    todo = vcat(cands, collect(values(refs)))
    out = o["--outdir"]; mkpath(joinpath(out, "img"))
    meta = Vector{Any}(undef, length(todo))
    Threads.@threads for k in eachindex(todo)    # pixels and ROI in parallel; plotting below stays serial (GR)
        t = todo[k]; path = joinpath(raw, t["raw"])
        row, v = inspect_scan(path, pcfg, Inf)
        img = read_sxm(path)
        meta[k] = (file=t["file"], session=t["session"], time=acquired(img), centre=centre(img),
                   frame=row["range_nm"], flags=row["flags"], note=t["note"], row=row, view=v)
    end
    for m in meta
        scan_figure(merge(m.row, Dict("file" => m.file)), m.view, joinpath(out, "img", replace(m.file, r"\.sxm$"i => ".png")))
    end
    scans = sort([m for m in meta[1:length(cands)]]; by=m -> (m.session, m.time, m.file))
    groups = molecule_groups(scans, Float64(sel["collect_max_gap_min"]), Float64(sel["collect_max_offset_nm"]))
    review_sheets([merge(m.row, Dict("file" => m.file, "sha256" => m.file)) for m in scans],
                  Dict(m.file => m.view for m in scans), joinpath(out, "sheets"), "Candidats par molécule";
                  tag=r -> "g$(groups[findfirst(m -> m.file == r["file"], scans)])")
    items = map(enumerate(scans)) do (k, m)
        s = get(sug, m.file, Dict("suggestion" => "benchmark", "motif" => "", "note" => ""))
        "{f:$(js(m.file)),s:$(js(m.session)),g:$(groups[k]),i:$(js("img/" * replace(m.file, r"\.sxm$"i => ".png"))),fr:$(js(m.frame)),fl:$(js(m.flags)),n:$(js(m.note)),sg:$(js(s["suggestion"])),sm:$(js(s["motif"])),sn:$(js(s["note"]))}"
    end
    refitems = ["{f:$(js(m.file)),s:$(js(m.session)),i:$(js("img/" * replace(m.file, r"\.sxm$"i => ".png")))}"
                for m in sort([m for m in meta[length(cands)+1:end]]; by=m -> m.session)]
    html = read(joinpath(@__DIR__, "candidate_review_template.html"), String)
    html = replace(html, "/*SCANS*/" => join(items, ",\n"), "/*REFS*/" => join(refitems, ",\n"),
                   "/*MOTIFS*/" => join(js.(MOTIFS), ","))
    write(joinpath(out, "index.html"), html)
    println("Review page: $(length(scans)) candidates in $(maximum(groups)) molecule groups, $(length(refitems)) reference scans -> $(joinpath(out, "index.html"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && review_main()

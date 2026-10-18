#!/usr/bin/env julia
# Kept/discarded triage of raw 6-mer scans for the benchmark (grading side).
#
# Usage:
#   review sheets: julia -t 4 --project=. test/report_benchmark_triage.jl --catalog CATALOG_TSV \
#                      [--catalog CATALOG_TSV ...] --count-config config/chitosan.toml --sheets NEW_DIR
#   apply:         julia -t 4 --project=. test/report_benchmark_triage.jl --catalog CATALOG_TSV [...] \
#                      --count-config config/chitosan.toml --triage TRIAGE_TSV --dest NEW_DIR
#
# The catalogs (catalog_sxm_candidates.jl) give every scan once by sha256, the
# first catalog winning. --sheets draws the usable scans for review: an
# 8 x 8 nm window centred on the molecule (same scale on every tile, 2 nm bar,
# unscanned pixels blank), tagged [B] in the benchmark manifest, [R] reviewed
# before, [N] never reviewed. --triage (columns file, session, folder,
# decision garde|ecarte, motif, source, note, sha256, origin; one row per
# scan) is applied by copying every scan to DEST/<folder>/<file> (folder
# gardes/... for garde, ecartes/... for ecarte), with DEST/LISTE.tsv (one row
# per scan), DEST/LISEZMOI.md (counts and meaning of every folder) and review
# sheets in DEST/planches/<folder>. Benchmark membership is a human
# decision; nothing here feeds the method.
include(joinpath(@__DIR__, "report_benchmark_candidates.jl"))

const CROP_NM = 8.0
const FOLDER_TEXT = Dict(
    "gardes/benchmark" => "utilisés par les runs du benchmark (manifeste `benchmarks/chitosan_6mer_counting_confirmed.toml`, plus `240310_Cu100009` utilisé mais non noté)",
    "gardes/candidats" => "gardés au premier tri visuel, à valider avant d'entrer dans le benchmark",
    "ecartes/vue_large" => "vue d'ensemble, cadre > 20 nm",
    "ecartes/molecule_coupee" => "molécule coupée par le bord du cadre ou par la zone non scannée",
    "ecartes/incertain" => "cas limites du premier tri, à revoir en priorité",
    "ecartes/revue_anterieure_ambigu" => "jugés ambigus ou exclus lors de la revue antérieure",
    "ecartes/revue_anterieure_exclu" => "rejetés à l'œil lors de la revue antérieure (non chitosan, etc.)",
    "ecartes/pas_de_molecule" => "pas de molécule nette (bruit, marche, bord de scan)",
    "ecartes/artefact_pointe" => "artefact sur la molécule (saut de ligne, stries, image dédoublée, changement de pointe)",
    "ecartes/plusieurs_objets" => "plusieurs molécules, objets accolés ou agrégat",
    "ecartes/pas_6mere_long" => "chaîne nettement plus longue qu'un 6-mère",
    "ecartes/pas_6mere_court" => "objet nettement plus court qu'un 6-mère")

function write_readme(path, tri, triage_path)
    folders = sort(unique(t["folder"] for t in tri); by=f -> (startswith(f, "ecartes") , f))
    open(path, "w") do io
        println(io, "# Scans (NKN)2 sur Cu(100) : dossier de référence\n")
        println(io, "Chaque scan Cu(100) des sessions (NKN)2, du disque local et du serveur hive")
        println(io, "(`esi/DATA/STM/STM_DATA/(NKN)2`), est ici une seule fois (doublons retirés par sha256).")
        println(io, "Généré le $(Dates.today()) par `test/report_benchmark_triage.jl` (dépôt STMFit) depuis")
        println(io, "`$(basename(triage_path))`. Ne pas modifier à la main : changer le tableau, puis reconstruire.\n")
        println(io, "| Dossier | Scans | Contenu |\n|---|---:|---|")
        for f in folders
            println(io, "| `$f/` | $(count(t -> t["folder"] == f, tri)) | $(get(FOLDER_TEXT, f, "")) |")
        end
        ng = count(t -> t["decision"] == "garde", tri)
        println(io, "| **total** | **$(length(tri))** | $ng gardés, $(length(tri) - ng) écartés |\n")
        println(io, "- `LISTE.tsv` : une ligne par scan (fichier, dossier, décision, motif, qui a décidé, note,")
        println(io, "  session, date et heure d'acquisition, cadre, pixels, origine `<session>/<fichier>` dans")
        println(io, "  l'arbre (NKN)2, sha256).")
        println(io, "- `planches/<dossier>/` : vignettes de relecture, cadre 8 × 8 nm centré sur la molécule,")
        println(io, "  barre de 2 nm.")
        println(io, "- Qui a décidé : `auto` (règle technique), `benchmark_actuel`, `revue_anterieure`,")
        println(io, "  `premier_tri_assistant` (premier tri visuel, à vérifier).")
    end
end

function triage_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); cats = String[]; i = 1
    while i <= length(args)
        k = args[i]
        k in ("--catalog", "--count-config", "--sheets", "--triage", "--dest") && i < length(args) ||
            error("Unknown option or missing value: $k")
        k == "--catalog" ? push!(cats, args[i+1]) : (o[k] = args[i+1]); i += 2
    end
    (isempty(cats) || !haskey(o, "--count-config")) && error("Required: --catalog --count-config")
    xor(haskey(o, "--sheets"), haskey(o, "--triage") && haskey(o, "--dest")) || error("Use --sheets DIR or --triage TSV --dest DIR")
    for k in ("--sheets", "--dest")
        haskey(o, k) && (ispath(o[k]) || islink(o[k])) && error("Output exists: $(o[k])")
    end
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
    parsed === nothing && return println("report_benchmark_triage.jl --catalog TSV [--catalog TSV] --count-config TOML (--sheets NEW_DIR | --triage TSV --dest NEW_DIR)")
    o, cats = parsed
    pcfg = count_preprocessing(o["--count-config"])
    rows = sort(scan_universe(cats); by=r -> (r["session"], r["file"]))
    if haskey(o, "--sheets")
        usable = [r for r in rows if r["status"] == "usable"]
        groups = Dict(r["sha256"] => r["group"] for r in candidate_groups(usable,
                      TOML.parsefile(DEFAULT_MANIFEST), grading_rows(DEFAULT_REVIEW)))
        tags = Dict("in_benchmark" => "[B]", "reviewed_before" => "[R]", "never_reviewed" => "[N]")
        pages = review_sheets(usable, scan_views(usable, pcfg), o["--sheets"], "Review"; tag=r -> tags[groups[r["sha256"]]])
        return println("Review sheets: $(length(usable)) usable scans, $pages pages")
    end
    _, tri = read_table(o["--triage"])
    by_sha = Dict(r["sha256"] => r for r in rows)
    length(tri) == length(rows) && Set(t["sha256"] for t in tri) == Set(keys(by_sha)) ||
        error("The triage must list every catalog scan once ($(length(tri)) rows, $(length(rows)) scans)")
    length(unique(t["file"] for t in tri)) == length(tri) || error("Duplicate file names in the triage")
    all(t -> t["decision"] in ("garde", "ecarte"), tri) || error("decision must be garde or ecarte")
    for t in tri
        f = t["folder"]
        (occursin(r"^(gardes|ecartes)/[A-Za-z0-9_]+$", f) && startswith(f, t["decision"] * "s/")) ||
            error("Folder $f does not match decision $(t["decision"]) for $(t["file"])")
    end
    dest = o["--dest"]
    for t in tri
        mkpath(joinpath(dest, t["folder"]))
        cp(by_sha[t["sha256"]]["path"], joinpath(dest, t["folder"], t["file"]))
    end
    sorted = sort(tri; by=t -> (t["session"], t["file"]))
    cols = ["fichier", "dossier", "decision", "motif", "decide_par", "note", "session", "date", "heure",
            "cadre_nm", "pixels", "origine", "sha256"]
    write_table(joinpath(dest, "LISTE.tsv"), cols, [begin
        r = by_sha[t["sha256"]]
        Dict("fichier" => t["file"], "dossier" => t["folder"], "decision" => t["decision"], "motif" => t["motif"],
             "decide_par" => t["source"], "note" => t["note"], "session" => t["session"], "date" => get(r, "rec_date", ""),
             "heure" => get(r, "rec_time", ""), "cadre_nm" => get(r, "range_nm", ""), "pixels" => get(r, "pixels", ""),
             "origine" => t["origin"], "sha256" => t["sha256"])
    end for t in sorted])
    write_readme(joinpath(dest, "LISEZMOI.md"), tri, o["--triage"])
    views = scan_views(rows, pcfg)
    named(t) = merge(by_sha[t["sha256"]], Dict("file" => t["file"]))
    for f in sort(unique(t["folder"] for t in tri))
        part = [named(t) for t in sorted if t["folder"] == f && haskey(views, t["sha256"])]
        review_sheets(part, views, joinpath(dest, "planches", f), f)
    end
    n = Dict(d => count(t -> t["decision"] == d, tri) for d in ("garde", "ecarte"))
    println("Triage applied to $dest: $(n["garde"]) gardés, $(n["ecarte"]) écartés")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && triage_main()

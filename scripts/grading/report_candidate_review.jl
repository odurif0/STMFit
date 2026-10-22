#!/usr/bin/env julia
# Local review page for the 6-mer benchmark candidates (grading side, human decision).
#
# Usage:
#   julia -t 4 --project=. scripts/grading/report_candidate_review.jl --data data \
#       --triage benchmarks/chitosan_6mer_cu100_triage.tsv --outdir NEW_DIR [--suggestions TSV] \
#       [--profile 6mer|10_20mer]
#
# --profile picks the candidate and reference sets, the decision vocabulary,
# the motifs, the language and the window size (6mer, the default: French,
# sets 6mer_Cu100_candidates / 6mer_Cu100_benchmark, decisions
# benchmark|ecarte, 8 nm; 10_20mer: English for several reviewers, sets
# 10_20mer_candidates / 10_20mer_curated, decisions keep|discard, 16 nm, a
# reviewer name exported with every decision).
#
# Takes the scans of the candidate set and groups repeat scans of one
# molecule: same session, consecutive acquisitions at most collect_max_gap_min
# apart, frame centres within collect_max_offset_nm (config/molecule_consensus.toml).
# Each scan is drawn twice (whole frame; square window on the molecule, 2 nm
# bar) in img/, and sheets/ holds contact sheets in group order. index.html is
# a local page: one card per scan with the suggestion (--suggestions: columns
# file, suggestion KEEP|DISCARD, motif, note; default KEEP) and the
# first-pass note, buttons to decide, and an export of the decisions as TSV
# (file, decision, motif, note, suggestion). One scan of the reference set per
# session is shown for reference. Nothing leaves the machine; decisions stay in
# the browser (localStorage) until exported.
include(joinpath(@__DIR__, "report_benchmark_triage.jl"))

# French interface strings of the 6-mer page (rendered byte-identically to the
# original page); the 10-20mer page is in English for colleagues.
const FRENCH_UI = (lang="fr", discard="ecarte", discard_verb="écarter", discard_label="Écarter",
    default_motif="incertain", t_all="Tous", t_todo="Non décidés", t_disagree="Différents de la suggestion",
    t_fill="Suggestions → non décidés", t_export="Exporter (TSV)", t_reset="Tout effacer", t_molecule="Molécule",
    t_all_to="Tout →", t_suggestion="suggestion : ", t_note="note (facultatif)", t_decided="décidés", t_discarded="écartés",
    t_confirm="Effacer toutes les décisions de cette page ?", key_keep="b", key_discard="e",
    reviewer_input="", reviewer_js="", export_check="", export_head="", export_col="",
    frame_title="cadre complet", window_title="fenêtre", sheets_title="Candidats par molécule")

const PROFILES = Dict(
    "6mer" => merge(FRENCH_UI, (candidates="6mer_Cu100_candidates", reference="6mer_Cu100_benchmark", keep="benchmark",
        keep_label="Benchmark", window_nm=8.0, title="Candidats 6-mères", key="stmfit_6mer_candidates_v1",
        download="decisions_candidats_6mer.tsv", refs_title="Référence : un scan du benchmark actuel par session",
        motifs=["pas_chitosan", "pas_6mere_court", "pas_6mere_long", "plusieurs_objets", "molecule_coupee",
                "artefact_pointe", "pas_de_molecule", "incertain"],
        header=join(["  <h1>Candidats 6-mères (NKN)2 sur Cu(100) : entrée au benchmark</h1>",
            "  <p class=\"crit\">Benchmark = <b>une molécule (NKN)2 entière, isolée et entièrement imagée</b> (vérité NKNNKN, 6 unités).",
            "  Écarter les fragments, agrégats, molécules coupées, images dominées par un artefact, et ce qui n'est pas du chitosan.",
            "  On juge l'image, jamais les sorties de la méthode. Raccourcis : survoler une carte puis <b>b</b> (benchmark), <b>e</b> (écarter), <b>u</b> (annuler).</p>"], "\n"))),
    # English page for several reviewers: each enters a name, exported with every decision.
    "10_20mer" => (candidates="10_20mer_candidates", reference="10_20mer_curated", keep="keep", keep_label="Keep",
        window_nm=16.0, title="10–20mer candidates", key="stmfit_10_20mer_candidates_en_v1", download="triage_decisions_10_20mer.tsv",
        refs_title="Reference: one scan per session that the experimenters selected",
        motifs=["not_chitosan", "cut_chain", "folded_or_crossed", "several_objects", "short_fragment",
                "tip_artifact", "no_molecule", "uncertain"],
        header=join(["  <h1>Chitosan 10–20mer scans on Cu(100): triage for the application</h1>",
            "  <p class=\"crit\">Keep = <b>one whole chitosan chain, isolated, fully imaged and not folded</b> (10 to 20 units expected; the exact number is unknown).",
            "  Discard short fragments, aggregates, chains cut by the frame or by the end of the scan, folded or crossing chains, images dominated by an artifact, and anything that is not chitosan.",
            "  As a guide, the 6-mer control molecules measure about 3 nm on these images: a straight chain shorter than about 4 nm has the length of a 6–7-mer.",
            "  Judge the image only, never the outputs of the analysis. Each card shows the whole frame and a 16 nm window centred on the object (2 nm scale bar); cards are grouped by molecule (repeat scans of one molecule).",
            "  The suggestions come from a first pass by the assistant: please check them. Shortcuts: hover a card, then <b>k</b> (keep), <b>d</b> (discard), <b>u</b> (undo).",
            "  Enter your name, then export your decisions (TSV) and send the file back. Decisions are saved in this browser until exported; nothing leaves your computer.</p>"], "\n"),
        lang="en", discard="discard", discard_verb="discard", discard_label="Discard", default_motif="uncertain",
        t_all="All", t_todo="Undecided", t_disagree="Different from suggestion", t_fill="Suggestions → undecided",
        t_export="Export (TSV)", t_reset="Clear all", t_molecule="Molecule", t_all_to="All →", t_suggestion="suggestion: ",
        t_note="note (optional)", t_decided="decided", t_discarded="discarded", t_confirm="Clear all decisions on this page?",
        key_keep="k", key_discard="d",
        reviewer_input="\n    <label>Your name <input id=\"reviewer\" size=\"14\" placeholder=\"required for export\"></label>",
        reviewer_js=join(["", "const rv = document.getElementById(\"reviewer\");",
            "try { rv.value = localStorage.getItem(KEY + \"_reviewer\") || \"\"; } catch (e) {}",
            "rv.onchange = () => { try { localStorage.setItem(KEY + \"_reviewer\", rv.value); } catch (e) {} };",
            "const REVIEWER = () => rv.value.replace(/[\\t\\n]/g, \" \").trim();"], "\n"),
        export_check="\n  if (!REVIEWER()) { alert(\"Please enter your name in the top bar before exporting.\"); rv.focus(); return; }",
        export_head=" + \"\\treviewer\"", export_col=", REVIEWER()",
        frame_title="full frame", window_title="window", sheets_title="Candidates by molecule"))

function review_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        args[i] in ("--data", "--triage", "--outdir", "--suggestions", "--profile") && i < length(args) ||
            error("Unknown option or missing value: $(args[i])")
        o[args[i]] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--data", "--triage", "--outdir")) || error("Required: --data --triage --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists: $(o["--outdir"])")
    haskey(PROFILES, get(o, "--profile", "6mer")) || error("--profile must be one of $(join(sort(collect(keys(PROFILES))), ", "))")
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

function scan_figure(row, v, path; window_nm=CROP_NM, frame_title="cadre complet", window_title="fenêtre")
    z = copy(v.z); z[.!v.observed] .= NaN
    vals = filter(isfinite, vec(z))
    lo, hi = isempty(vals) ? (0.0, 1.0) : (quantile(vals, 0.01), quantile(vals, 0.995))
    full = heatmap(v.xs, v.ys, z; c=:inferno, clims=(lo, max(hi, lo + 1e-6)), aspect_ratio=:equal, colorbar=false,
                   axis=false, grid=false, title="$frame_title $(row["range_nm"]) nm", titlefontsize=9)
    contour!(full, v.xs, v.ys, Float64.(v.mask .& v.observed); levels=[0.5], color=:cyan, linewidth=0.6, alpha=0.6, colorbar=false)
    w = isinteger(window_nm) ? string(Int(window_nm)) : string(window_nm)
    crop = review_tile(merge(row, Dict("flags" => "", "range_nm" => "$window_title $w x $w")), v, ""; crop_nm=window_nm)
    savefig(plot(full, crop; layout=(1, 2), size=(980, 470), margin=2Plots.mm), path)
end

js(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ", "<" => "\\u003c") * "\""

function review_main(args=ARGS)
    o = review_options(args)
    o === nothing && return println("report_candidate_review.jl --data DIR --triage TSV --outdir NEW_DIR [--suggestions TSV] [--profile 6mer|10_20mer]")
    P = PROFILES[get(o, "--profile", "6mer")]
    root = dirname(dirname(@__DIR__))
    sel = TOML.parsefile(joinpath(root, "config", "molecule_consensus.toml"))["selection"]
    pcfg = count_preprocessing(joinpath(root, "config", "chitosan.toml"))
    _, tri = read_table(o["--triage"])
    sug = Dict{String,Dict{String,String}}()
    haskey(o, "--suggestions") && (sug = Dict(r["file"] => r for r in last(read_table(o["--suggestions"]))))
    raw = joinpath(o["--data"], "raw")
    cands = [t for t in tri if t["set"] == P.candidates]
    refs = Dict{String,Dict{String,String}}()
    for t in sort([t for t in tri if t["set"] == P.reference]; by=t -> t["file"])
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
        scan_figure(merge(m.row, Dict("file" => m.file)), m.view, joinpath(out, "img", replace(m.file, r"\.sxm$"i => ".png"));
                    window_nm=P.window_nm, frame_title=P.frame_title, window_title=P.window_title)
    end
    scans = sort([m for m in meta[1:length(cands)]]; by=m -> (m.session, m.time, m.file))
    groups = molecule_groups(scans, Float64(sel["collect_max_gap_min"]), Float64(sel["collect_max_offset_nm"]))
    review_sheets([merge(m.row, Dict("file" => m.file, "sha256" => m.file)) for m in scans],
                  Dict(m.file => m.view for m in scans), joinpath(out, "sheets"), P.sheets_title;
                  tag=r -> "g$(groups[findfirst(m -> m.file == r["file"], scans)])", crop_nm=P.window_nm)
    items = map(enumerate(scans)) do (k, m)
        s = get(sug, m.file, Dict("suggestion" => P.keep, "motif" => "", "note" => ""))
        "{f:$(js(m.file)),s:$(js(m.session)),g:$(groups[k]),i:$(js("img/" * replace(m.file, r"\.sxm$"i => ".png"))),fr:$(js(m.frame)),fl:$(js(m.flags)),n:$(js(m.note)),sg:$(js(s["suggestion"])),sm:$(js(s["motif"])),sn:$(js(s["note"]))}"
    end
    refitems = ["{f:$(js(m.file)),s:$(js(m.session)),i:$(js("img/" * replace(m.file, r"\.sxm$"i => ".png")))}"
                for m in sort([m for m in meta[length(cands)+1:end]]; by=m -> m.session)]
    html = read(joinpath(@__DIR__, "candidate_review_template.html"), String)
    html = replace(html, "/*SCANS*/" => join(items, ",\n"), "/*REFS*/" => join(refitems, ",\n"),
                   "/*MOTIFS*/" => join(js.(P.motifs), ","))
    html = render_profile(html, P)
    write(joinpath(out, "index.html"), html)
    println("Review page: $(length(scans)) candidates in $(maximum(groups)) molecule groups, $(length(refitems)) reference scans -> $(joinpath(out, "index.html"))")
end

"Fill the profile placeholders of the page template."
render_profile(html, P) = replace(html, "/*HEADER*/" => P.header, "/*TITLE*/" => P.title, "/*REFS_TITLE*/" => P.refs_title,
    "/*KEY*/" => P.key, "/*KEEP_LABEL*/" => P.keep_label, "/*KEEP*/" => P.keep, "/*DOWNLOAD*/" => P.download,
    "/*LANG*/" => P.lang, "/*DISCARD_VERB*/" => P.discard_verb, "/*DISCARD_LABEL*/" => P.discard_label, "/*DISCARD*/" => P.discard,
    "/*DEFAULT_MOTIF*/" => P.default_motif, "/*T_ALL_TO*/" => P.t_all_to, "/*T_ALL*/" => P.t_all, "/*T_TODO*/" => P.t_todo,
    "/*T_DISAGREE*/" => P.t_disagree, "/*T_FILL*/" => P.t_fill, "/*T_EXPORT*/" => P.t_export, "/*T_RESET*/" => P.t_reset,
    "/*T_MOLECULE*/" => P.t_molecule, "/*T_SUGGESTION*/" => P.t_suggestion, "/*T_NOTE*/" => P.t_note,
    "/*T_DECIDED*/" => P.t_decided, "/*T_DISCARDED*/" => P.t_discarded, "/*T_CONFIRM*/" => P.t_confirm,
    "/*KEY_KEEP*/" => P.key_keep, "/*KEY_DISCARD*/" => P.key_discard, "/*REVIEWER_INPUT*/" => P.reviewer_input,
    "/*REVIEWER_JS*/" => P.reviewer_js, "/*EXPORT_CHECK*/" => P.export_check, "/*EXPORT_HEAD*/" => P.export_head,
    "/*EXPORT_COL*/" => P.export_col)

abspath(PROGRAM_FILE) == abspath(@__FILE__) && review_main()

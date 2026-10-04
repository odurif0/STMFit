#!/usr/bin/env julia
# One plot per molecule scan of a run_molecule_consensus_chitosan.jl result.
#
# Usage:
#   julia --project=. test/plot_final_molecules.jl --run RUN_DIR --data-dir RAW_DIR \
#       --config COUNT_CONFIG --outdir NEW_DIR [--annotations TSV]
#
# Each plot shows the preprocessed STM image around the chain, the final lobes
# (colour = final call, label = lobe index and P(GlcNAc)), their 2-sigma position
# ellipses from positions.tsv, and a title with the final and per-scan counts,
# the consensus rule, the call sequence and the confidence. Label-free: plots
# go to all/ and to folders by consensus outcome (consensus_agrees,
# consensus_corrected, not_checked_by_consensus) and with_uncertain_calls.
#
# --annotations (written by an external grading step, never by the method):
# columns file, category, note, wrong_lobes (comma-separated lobe indices).
# Plots then go to all/ and to one folder per category; the note is added to
# the title and wrong lobes are circled.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
using .ReconstructedUnitAssignment
using GaussianFit2D, STMSXMIO, Plots, Printf, TOML
const G = GaussianFit2D

const COLORS = Dict("0" => "#2166ac", "1" => "#d6604d", "?" => "#bdbdbd")
const NAMES = Dict("0" => "GlcN (0)", "1" => "GlcNAc (1)", "?" => "uncertain (?)")
const CORRECTED = ("consensus_applied", "consensus_registered_roi")

function plot_options(args)
    ("--help" in args || "-h" in args || isempty(args)) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in ("--run", "--data-dir", "--config", "--outdir", "--annotations") && i < length(args) ||
            error("Unknown option or missing value: $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in ("--run", "--data-dir", "--config", "--outdir")) ||
        error("Required: --run --data-dir --config --outdir")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

const RULE_TEXT = Dict("agrees" => "agrees with the consensus", "consensus_applied" => "corrected by the consensus",
    "consensus_registered_roi" => "corrected by the consensus (registered ROI)",
    "track_too_short" => "fewer than 3 scans, not checked", "no_strict_majority" => "no strict majority, not checked",
    "footprint_outside_frame" => "molecule outside the frame, not checked",
    "consensus_fit_failed" => "consensus refit failed, scan count kept")

"Label-free folder of a scan from its consensus rule."
consensus_folder(rule) = rule == "agrees" ? "consensus_agrees" :
    rule in CORRECTED ? "consensus_corrected" : "not_checked_by_consensus"

function ellipse(x, y, ax, sa, sc; k=2.0, m=40)
    pe = (-ax[2], ax[1]); th = range(0, 2pi; length=m)
    ([x + k * (sa * cos(t) * ax[1] + sc * sin(t) * pe[1]) for t in th],
     [y + k * (sa * cos(t) * ax[2] + sc * sin(t) * pe[2]) for t in th])
end

function molecule_plot(img, pcfg, file, lobes, report, note, wrong)
    xs, ys, _, _, z, _, _ = G.preprocess_channel(img, get_channel(img, "Z"; direction="fwd"), pcfg)
    lx = [l.x for l in lobes]; ly = [l.y for l in lobes]; pad = 1.5
    ix = findall(x -> minimum(lx) - pad <= x <= maximum(lx) + pad, xs)
    iy = findall(y -> minimum(ly) - pad <= y <= maximum(ly) + pad, ys)
    isempty(ix) || isempty(iy) ? (ix, iy) = (eachindex(xs), eachindex(ys)) : nothing
    calls = join(l.call for l in lobes)
    nscan = report["N_scan"]; nfin = report["N_selected"]
    line1 = nscan == nfin ? "$file   N = $nfin" : "$file   N = $nfin (this scan alone: $nscan)"
    conf = something(tryparse(Float64, report["mean_confidence"]), NaN)
    rule = get(RULE_TEXT, report["count_rule"], report["count_rule"])
    line2 = "imaged $(report["track_size"])×: $rule   ·   calls $calls   ·   mean confidence " *
            (isnan(conf) ? "NA" : @sprintf("%.2f", conf))
    title = isempty(note) ? "$line1\n$line2" : "$line1\n$line2\n$note"
    fig = heatmap(xs[ix], ys[iy], z[iy, ix]; c=:inferno, aspect_ratio=:equal, colorbar=false,
                  xlabel="x (nm)", ylabel="y (nm)", title=title, titlefontsize=9,
                  xlims=(xs[first(ix)], xs[last(ix)]), ylims=(ys[first(iy)], ys[last(iy)]),
                  size=(820, 720), margin=4Plots.mm, legend=:outerbottom)
    plot!(fig, lx, ly; color=:white, alpha=0.6, linewidth=1, label="")
    ax = length(lobes) >= 2 ? (lx[end] - lx[1], ly[end] - ly[1]) : (1.0, 0.0)
    ax = ax ./ max(hypot(ax...), eps())
    for l in lobes
        isfinite(l.sda) && isfinite(l.sdc) && plot!(fig, ellipse(l.x, l.y, ax, l.sda, l.sdc)...; color=:white, alpha=0.8, linewidth=0.8, label="")
    end
    if !isempty(wrong)
        w = [l for l in lobes if l.lobe in wrong]
        scatter!(fig, [l.x for l in w], [l.y for l in w]; markersize=16, markershape=:circle,
                 markercolor=:transparent, markerstrokecolor=:cyan, markerstrokewidth=3, label="wrong call")
    end
    for c in ("0", "1", "?")
        idx = findall(l -> l.call == c, lobes); isempty(idx) && continue
        scatter!(fig, lx[idx], ly[idx]; color=COLORS[c], markersize=9, markerstrokecolor=:black, label=NAMES[c])
    end
    for l in lobes
        p = isnan(l.p1) ? "?" : @sprintf("%.2f", l.p1)
        annotate!(fig, l.x, l.y + 0.32, text("$(l.lobe): $p", 7, :white))
    end
    plot!(fig, [NaN], [NaN]; color=:white, label="2σ position ellipse")
    nlegend = length(unique(l.call for l in lobes)) + 1 + (isempty(wrong) ? 0 : 1)
    plot!(fig; legendcolumns=nlegend)
    return fig
end

function main(args=ARGS)
    o = plot_options(args)
    o === nothing && return println("plot_final_molecules.jl --run RUN_DIR --data-dir RAW_DIR --config COUNT_CONFIG --outdir NEW_DIR [--annotations TSV]")
    run = o["--run"]; pre = get(TOML.parsefile(o["--config"]), "preprocessing", Dict{String,Any}())
    pcfg = G.PatternConfig(filepath="", channel="Z", direction="fwd", stride=get(pre, "stride", 1),
        flatten=get(pre, "flatten", "plane+rows"), smooth_radius_px=get(pre, "smooth_radius_px", 1), no_plot=true)
    _, reports = read_table(joinpath(run, "chain_report.tsv"))
    _, pred = lobe_table(joinpath(run, "predictions.tsv"); required=["predicted", "probability_1"])
    _, geo = lobe_table(joinpath(run, "assignment", "features.tsv"); required=["x_nm", "y_nm"])
    pos = isfile(joinpath(run, "positions.tsv")) ? last(lobe_table(joinpath(run, "positions.tsv"))) : Dict()
    ann = Dict{String,Dict{String,String}}()
    if haskey(o, "--annotations")
        _, rows = read_table(o["--annotations"])
        ann = Dict(r["file"] => r for r in rows)
    end
    raw = Dict{String,String}()
    for (d, _, names) in walkdir(o["--data-dir"]), n in names
        endswith(lowercase(n), ".sxm") && (raw[n] = joinpath(d, n))
    end
    out = o["--outdir"]; mkpath(joinpath(out, "all"))
    index = Dict{String,String}[]
    f64(s) = something(tryparse(Float64, s), NaN)
    for r in sort(reports; by=r -> r["file"])
        file = r["file"]; n = parse(Int, r["N_selected"])
        lobes = [(lobe=i, x=f64(geo[(file, i)]["x_nm"]), y=f64(geo[(file, i)]["y_nm"]),
                  call=pred[(file, i)]["predicted"], p1=f64(pred[(file, i)]["probability_1"]),
                  sda=haskey(pos, (file, i)) ? f64(pos[(file, i)]["sd_along_nm"]) : NaN,
                  sdc=haskey(pos, (file, i)) ? f64(pos[(file, i)]["sd_across_nm"]) : NaN) for i in 1:n]
        a = get(ann, file, nothing)
        wrong = a === nothing || isempty(a["wrong_lobes"]) ? Int[] : parse.(Int, split(a["wrong_lobes"], ","))
        rep = Dict("N_scan" => parse(Int, r["N_scan"]), "N_selected" => n, "count_rule" => r["count_rule"],
                   "track_size" => r["track_size"], "mean_confidence" => r["mean_confidence"])
        fig = molecule_plot(read_sxm(raw[file]), pcfg, file, lobes, rep, a === nothing ? "" : a["note"], wrong)
        png = replace(file, r"\.sxm$"i => "") * ".png"
        savefig(fig, joinpath(out, "all", png))
        folders = a !== nothing ? [a["category"]] :
            vcat([consensus_folder(r["count_rule"])], any(l -> l.call == "?", lobes) ? ["with_uncertain_calls"] : String[])
        for f in folders
            mkpath(joinpath(out, f)); cp(joinpath(out, "all", png), joinpath(out, f, png))
        end
        push!(index, Dict("file" => file, "folders" => join(folders, ","), "N" => string(n),
                          "count_rule" => r["count_rule"], "calls" => join(l.call for l in lobes), "png" => "all/" * png))
    end
    write_table(joinpath(out, "index.tsv"), ["file", "folders", "N", "count_rule", "calls", "png"], index)
    counts = Dict{String,Int}()
    for r in index, f in split(r["folders"], ","); counts[f] = get(counts, f, 0) + 1; end
    println("Plotted $(length(index)) molecule scans: ", join(["$k $v" for (k, v) in sort(collect(counts))], ", "))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

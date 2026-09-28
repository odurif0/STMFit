#!/usr/bin/env julia
# Repeated-scan grouping and molecule-level count consensus (label-free).
#
# Usage:
#   julia --project=. test/build_molecule_consensus.jl --data-dir RAW_DIR \
#       --selected-summary counting_summary.tsv --features base_features.tsv \
#       --config config/molecule_consensus.toml --outdir NEW_DIR
#
# Inputs are raw SXM files, per-scan selected counts and the fixed-N geometry
# fitted at those counts. Outputs: pairs.tsv (consecutive registrations),
# consensus.tsv (per-scan provenance) and consensus_summary.tsv
# (filepath/status/N_selected for the downstream assignment runner).
# No benchmark label, expected count, sequence or class proportion is read.
include(joinpath(@__DIR__, "lib", "reconstructed_unit_assignment.jl"))
include(joinpath(@__DIR__, "lib", "molecule_consensus.jl"))
using .ReconstructedUnitAssignment, .MoleculeConsensus
using STMSXMIO, TOML, Printf, Statistics

const OPTIONS = Set(["--data-dir", "--selected-summary", "--features", "--config", "--outdir"])

function consensus_options(args)
    ("--help" in args || "-h" in args) && return nothing
    o = Dict{String,String}(); i = 1
    while i <= length(args)
        k = args[i]
        k in OPTIONS && i < length(args) || error("Unknown option or missing value: $k")
        haskey(o, k) && error("Repeated option $k")
        o[k] = args[i+1]; i += 2
    end
    all(haskey(o, k) for k in OPTIONS) || error("Required: $(join(sort(collect(OPTIONS)), ", "))")
    isdir(o["--data-dir"]) || error("Missing raw directory")
    all(isfile(o[k]) for k in ("--selected-summary", "--features", "--config")) || error("Missing input")
    (ispath(o["--outdir"]) || islink(o["--outdir"])) && error("Output exists")
    return o
end

function scan_counts(path)
    header, rows = read_table(path)
    all(c -> c in header, ("filepath", "N_selected")) || error("Summary requires filepath and N_selected")
    counts = Dict{String,Int}()
    for r in rows
        f = basename(r["filepath"]); haskey(counts, f) && error("Duplicate scan $f")
        get(r, "status", "ok") == "ok" || error("Counting did not succeed: $f")
        n = tryparse(Int, r["N_selected"]); n !== nothing && n > 0 || error("Invalid N_selected for $f")
        counts[f] = n
    end
    return counts
end

function raw_paths(dir)
    index = Dict{String,String}()
    for (d, _, names) in walkdir(dir), name in names
        endswith(lowercase(name), ".sxm") || continue
        haskey(index, name) && error("Duplicate raw basename $name")
        index[name] = joinpath(d, name)
    end
    return index
end

function build_consensus(o)
    st = load_consensus_settings(TOML.parsefile(o["--config"]))
    counts = scan_counts(o["--selected-summary"])
    _, lobes = lobe_table(o["--features"]; required=["x_nm", "y_nm"])
    geometry = Dict{String,Vector{NTuple{2,Float64}}}()
    for key in sort(collect(keys(lobes)))
        r = lobes[key]
        push!(get!(geometry, key[1], NTuple{2,Float64}[]), (parse(Float64, r["x_nm"]), parse(Float64, r["y_nm"])))
    end
    Set(keys(geometry)) == Set(keys(counts)) || error("Geometry and count cohorts differ")
    for (f, n) in counts
        length(geometry[f]) == n || error("Geometry rows differ from N_selected for $f")
    end
    raw = raw_paths(o["--data-dir"])
    all(f -> haskey(raw, f), keys(counts)) || error("Raw SXM missing for a counted scan")
    frames = Dict{String,ScanFrame}(); lobes_abs = Dict{String,Vector{NTuple{2,Float64}}}()
    centroids = Dict{String,NTuple{2,Float64}}(); cache = Dict{String,Matrix{Float64}}()
    for f in sort(collect(keys(counts)))
        img = read_sxm(raw[f]); fr = scan_frame(img); frames[f] = fr
        lobes_abs[f] = [absolute_xy(fr, p...) for p in geometry[f]]
        centroids[f] = (mean(first.(lobes_abs[f])), mean(last.(lobes_abs[f])))
    end
    images(f) = get!(() -> registration_image(read_sxm(raw[f]), st), cache, f)
    order, links = link_scans(frames, centroids, images, st)
    rows = consensus_counts(order, links, counts, lobes_abs, frames, st)
    out = abspath(o["--outdir"]); mkpath(out)
    fmt(x) = x isa AbstractFloat ? (isfinite(x) ? @sprintf("%.6f", x) : "NaN") : string(x)
    write_table(joinpath(out, "pairs.tsv"),
        ["a", "b", "center_distance_nm", "ncc", "tx_nm", "ty_nm", "linked", "reason"],
        [Dict(string(k) => fmt(v) for (k, v) in pairs(p)) for p in links])
    cols = ["file", "track", "track_size", "N_scan", "N_consensus", "agreement", "N_final", "reference_scan", "rule"]
    write_table(joinpath(out, "consensus.tsv"), cols, [Dict(k => fmt(r[k]) for k in cols) for r in rows])
    write_table(joinpath(out, "consensus_summary.tsv"),
        ["filepath", "status", "N_selected", "N_scan", "track", "track_size", "count_agreement", "count_rule"],
        [Dict("filepath" => r["file"], "status" => "ok", "N_selected" => string(r["N_final"]),
              "N_scan" => string(r["N_scan"]), "track" => string(r["track"]),
              "track_size" => string(r["track_size"]), "count_agreement" => fmt(r["agreement"]),
              "count_rule" => r["rule"]) for r in sort(rows; by=r -> r["file"])])
    changed = count(r -> r["N_final"] != r["N_scan"], rows)
    println("Molecule consensus: $(length(rows)) scans, $(length(unique(r["track"] for r in rows))) tracks, $changed counts changed")
    return joinpath(out, "consensus_summary.tsv")
end

function main(args=ARGS)
    o = consensus_options(args)
    o === nothing ? println("build_molecule_consensus.jl --data-dir DIR --selected-summary TSV --features TSV --config TOML --outdir NEW_DIR") :
        build_consensus(o)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

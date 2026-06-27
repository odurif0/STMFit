#!/usr/bin/env julia

# Label-free refinement of the manual geometric GlcN/GlcNAc mold.
#
# The script explores a small, explicitly parameterized family around
# templates/chitosan_geometric_sites.tsv. It never reads truth labels and never
# imposes a GlcNAc/GlcN composition. Candidate molds are scored only by the
# bimodality and margin of their template evidence on aligned patches.

using Printf
using Statistics
using LinearAlgebra
using Random
using Clustering

include(joinpath(@__DIR__, "lib", "script_utils.jl"))
using .ScriptUtils: _parse_f, _read_tsv, _standardize

const DEFAULT_PATCHES = "results/unit_separability/lobe_patches_selectedN_primary_half048.tsv"
const DEFAULT_SITES = "templates/chitosan_geometric_sites.tsv"
const DEFAULT_OUT_SITES = "templates/chitosan_geometric_sites_refined.tsv"
const DEFAULT_REPORT = "results/unit_assignment/geometric_mold_refinement.tsv"

struct Options
    patches::String
    sites::String
    out_sites::String
    report::String
    prefix::String
    step_nm::Float64
    acetyl_t_shifts::Vector{Float64}
    acetyl_u_shifts::Vector{Float64}
    acetyl_u_scales::Vector{Float64}
    acetyl_weight_scales::Vector{Float64}
    acetyl_sigma_scales::Vector{Float64}
end

function _parse_list(s::String)
    vals = Float64[]
    for part in split(s, ',')
        t = strip(part)
        isempty(t) || push!(vals, parse(Float64, t))
    end
    isempty(vals) && error("empty numeric list")
    return vals
end

function _parse_cli(args)
    patches = DEFAULT_PATCHES
    sites = DEFAULT_SITES
    out_sites = DEFAULT_OUT_SITES
    report = DEFAULT_REPORT
    prefix = "raw_p"
    step_nm = 0.08
    acetyl_t_shifts = [-0.08, -0.04, 0.0, 0.04, 0.08]
    acetyl_u_shifts = [-0.08, -0.04, 0.0, 0.04, 0.08]
    acetyl_u_scales = [0.85, 1.0, 1.15]
    acetyl_weight_scales = [0.5, 0.8, 1.1, 1.4, 1.8]
    acetyl_sigma_scales = [0.75, 1.0, 1.25]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--patches"; patches = args[i+1]; i += 2
        elseif startswith(arg, "--patches="); patches = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--sites"; sites = args[i+1]; i += 2
        elseif startswith(arg, "--sites="); sites = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--out-sites"; out_sites = args[i+1]; i += 2
        elseif startswith(arg, "--out-sites="); out_sites = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--report"; report = args[i+1]; i += 2
        elseif startswith(arg, "--report="); report = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--prefix"; prefix = args[i+1]; i += 2
        elseif startswith(arg, "--prefix="); prefix = split(arg, "=", limit=2)[2]; i += 1
        elseif arg == "--step-nm"; step_nm = parse(Float64, args[i+1]); i += 2
        elseif startswith(arg, "--step-nm="); step_nm = parse(Float64, split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--acetyl-t-shifts"; acetyl_t_shifts = _parse_list(args[i+1]); i += 2
        elseif startswith(arg, "--acetyl-t-shifts="); acetyl_t_shifts = _parse_list(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--acetyl-u-shifts"; acetyl_u_shifts = _parse_list(args[i+1]); i += 2
        elseif startswith(arg, "--acetyl-u-shifts="); acetyl_u_shifts = _parse_list(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--acetyl-u-scales"; acetyl_u_scales = _parse_list(args[i+1]); i += 2
        elseif startswith(arg, "--acetyl-u-scales="); acetyl_u_scales = _parse_list(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--acetyl-weight-scales"; acetyl_weight_scales = _parse_list(args[i+1]); i += 2
        elseif startswith(arg, "--acetyl-weight-scales="); acetyl_weight_scales = _parse_list(split(arg, "=", limit=2)[2]); i += 1
        elseif arg == "--acetyl-sigma-scales"; acetyl_sigma_scales = _parse_list(args[i+1]); i += 2
        elseif startswith(arg, "--acetyl-sigma-scales="); acetyl_sigma_scales = _parse_list(split(arg, "=", limit=2)[2]); i += 1
        elseif arg in ("-h", "--help")
            println("""
            Usage: julia --project=. test/refine_geometric_mold.jl [options]

            Options:
              --patches PATH                Wide patch TSV [$(DEFAULT_PATCHES)]
              --sites PATH                  Base geometric site TSV [$(DEFAULT_SITES)]
              --out-sites PATH              Refined site TSV [$(DEFAULT_OUT_SITES)]
              --report PATH                 Candidate report TSV [$(DEFAULT_REPORT)]
              --prefix STR                  Patch prefix raw_p or res_p [raw_p]
              --step-nm FLOAT               Patch/template grid spacing [0.08]
              --acetyl-t-shifts LIST        Comma-separated nm shifts for acetyl sites
              --acetyl-u-shifts LIST        Comma-separated nm shifts for acetyl sites
              --acetyl-u-scales LIST        Comma-separated transverse scale factors
              --acetyl-weight-scales LIST   Comma-separated acetyl weight scale factors
              --acetyl-sigma-scales LIST    Comma-separated acetyl sigma scale factors

            Label-free objective:
              For each candidate, the script scores patches against contrast
              templates, chooses only global direction/phase/mirror states, and
              ranks candidates by k=1 vs k=2 BIC of the resulting per-lobe
              template evidence margins. Truth labels and composition counts are
              never read.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
    isfile(patches) || error("Patch TSV not found: $patches")
    isfile(sites) || error("Site TSV not found: $sites")
    step_nm > 0 || error("--step-nm must be positive")
    return Options(patches, sites, out_sites, report, prefix, step_nm,
                   acetyl_t_shifts, acetyl_u_shifts, acetyl_u_scales,
                   acetyl_weight_scales, acetyl_sigma_scales)
end

function _load_sites(path::String)
    _, rows = _read_tsv(path)
    required = ("type", "atom", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm")
    sites = NamedTuple[]
    for row in rows
        all(haskey(row, k) for k in required) || error("Site TSV missing one of: $(join(required, ", "))")
        typ = parse(Int, row["type"])
        typ in (0, 1) || error("type must be 0 or 1")
        push!(sites, (
            typ=typ,
            atom=strip(row["atom"]),
            t=_parse_f(row["t_nm"]),
            u=_parse_f(row["u_nm"]),
            weight=_parse_f(row["weight"]),
            sigt=max(_parse_f(row["sigma_t_nm"]), eps(Float64)),
            sigu=max(_parse_f(row["sigma_u_nm"]), eps(Float64)),
        ))
    end
    any(s.typ == 0 for s in sites) || error("No type=0 sites")
    any(s.typ == 1 for s in sites) || error("No type=1 sites")
    return sites
end

function _is_acetyl_site(site)
    return site.typ == 1 && startswith(lowercase(site.atom), "acetyl")
end

function _transform_sites(sites, tshift, ushift, uscale, wscale, sigscale)
    out = NamedTuple[]
    for s in sites
        if _is_acetyl_site(s)
            push!(out, (typ=s.typ, atom=s.atom,
                        t=s.t + tshift,
                        u=s.u * uscale + ushift,
                        weight=s.weight * wscale,
                        sigt=s.sigt * sigscale,
                        sigu=s.sigu * sigscale))
        else
            push!(out, s)
        end
    end
    return out
end

function _apply_flip(t, u, spec::String)
    spec == "none" && return t, u
    spec == "t" && return -t, u
    spec == "u" && return t, -u
    return -t, -u
end

function _template_values(sites, coords, typ::Int, parity::Int, mirror::Int)
    vals = Float64[]
    selected = [s for s in sites if s.typ == typ]
    for u in coords, t in coords
        v = 0.0
        for s in selected
            st, su = s.t, s.u
            parity == 1 && ((st, su) = _apply_flip(st, su, "t"))
            mirror == 1 && ((st, su) = _apply_flip(st, su, "u"))
            v += s.weight * exp(-0.5 * (((t - st) / s.sigt)^2 + ((u - su) / s.sigu)^2))
        end
        push!(vals, v)
    end
    return _standardize(vals)
end

function _contrast_templates(sites, coords)
    templates = Dict{Tuple{Int,Int,Int},Vector{Float64}}()
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        templates[(typ, parity, mirror)] = _template_values(sites, coords, typ, parity, mirror)
    end
    for parity in (0, 1), mirror in (0, 1)
        k0 = (0, parity, mirror)
        k1 = (1, parity, mirror)
        common = 0.5 .* (templates[k0] .+ templates[k1])
        templates[k0] = _standardize(templates[k0] .- common)
        templates[k1] = _standardize(templates[k1] .- common)
    end
    return templates
end

function _load_patches(path::String, prefix::String)
    header, rows = _read_tsv(path)
    pix_cols = [c for c in header if startswith(c, prefix)]
    isempty(pix_cols) && error("Patch TSV has no columns with prefix $prefix")
    by_file = Dict{String,Vector{NamedTuple}}()
    for row in rows
        file = basename(row["file"])
        lobe = parse(Int, row["lobe"])
        patch = _standardize([_parse_f(row[c]) for c in pix_cols])
        push!(get!(by_file, file, NamedTuple[]), (file=file, lobe=lobe, patch=patch))
    end
    for file in keys(by_file)
        sort!(by_file[file], by=r -> r.lobe)
    end
    return by_file, pix_cols
end

function _parity_for_lobe(lobe::Int, n::Int, direction::Int, phase::Int)
    idx = direction == 0 ? lobe : (n - lobe + 1)
    return mod(idx - 1 + phase, 2)
end

function _score_patch(patch::Vector{Float64}, templ::Vector{Float64})
    n_good = 0
    dot_pt = 0.0
    norm_p2 = 0.0
    norm_t2 = 0.0
    for i in eachindex(patch)
        p = patch[i]
        t = templ[i]
        if isfinite(p) && isfinite(t)
            n_good += 1
            dot_pt += p * t
            norm_p2 += p * p
            norm_t2 += t * t
        end
    end
    if n_good < max(5, cld(length(patch), 2))
        return Inf
    end
    denom = sqrt(norm_p2 * norm_t2)
    denom <= 0 && return Inf
    return -dot_pt / denom
end

function _candidate_margins(by_file, templates)
    margins = Float64[]
    total_cost = 0.0
    for file in sort(collect(keys(by_file)))
        records = by_file[file]
        n = length(records)
        best_total = Inf
        best_margins = Float64[]
        for direction in (0, 1), phase in (0, 1), mirror in (0, 1)
            state_margins = Float64[]
            state_total = 0.0
            for rec in records
                parity = _parity_for_lobe(rec.lobe, n, direction, phase)
                c0 = _score_patch(rec.patch, templates[(0, parity, mirror)])
                c1 = _score_patch(rec.patch, templates[(1, parity, mirror)])
                state_total += min(c0, c1)
                push!(state_margins, c0 - c1)
            end
            if state_total < best_total
                best_total = state_total
                best_margins = state_margins
            end
        end
        total_cost += best_total
        append!(margins, best_margins)
    end
    return margins, total_cost
end

function _kmeans_bic_1d(vals::Vector{Float64})
    x = reshape(vals, 1, length(vals))
    n = length(vals)
    center = mean(vals)
    sse1 = sum((vals .- center).^2)
    bic1 = n * log(sse1 / n + 1e-30) + 1 * log(n)
    km = kmeans(x, 2; maxiter=200, rng=MersenneTwister(42), display=:none)
    sse2 = km.totalcost
    bic2 = n * log(sse2 / n + 1e-30) + 3 * log(n)
    sizes = [count(==(1), km.assignments), count(==(2), km.assignments)]
    return bic1 - bic2, sse1, sse2, minimum(sizes), maximum(sizes)
end

function _write_sites(path::String, sites)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Refined geometric proxy sites. Generated label-free by test/refine_geometric_mold.jl.")
        println(io, join(["type", "atom", "t_nm", "u_nm", "weight", "sigma_t_nm", "sigma_u_nm"], '\t'))
        for s in sites
            println(io, join([s.typ, s.atom, @sprintf("%.8g", s.t), @sprintf("%.8g", s.u),
                              @sprintf("%.8g", s.weight), @sprintf("%.8g", s.sigt),
                              @sprintf("%.8g", s.sigu)], '\t'))
        end
    end
end

function main()
    opt = _parse_cli(ARGS)
    base_sites = _load_sites(opt.sites)
    by_file, pix_cols = _load_patches(opt.patches, opt.prefix)
    grid_n = round(Int, sqrt(length(pix_cols)))
    grid_n^2 == length(pix_cols) || error("Patch pixel count $(length(pix_cols)) is not a square")
    coords = collect((-(grid_n - 1) / 2):(grid_n - 1) / 2) .* opt.step_nm

    candidates = NamedTuple[]
    total_candidates = length(opt.acetyl_t_shifts) * length(opt.acetyl_u_shifts) *
        length(opt.acetyl_u_scales) * length(opt.acetyl_weight_scales) * length(opt.acetyl_sigma_scales)
    idx = 0
    for dt in opt.acetyl_t_shifts, du in opt.acetyl_u_shifts,
        us in opt.acetyl_u_scales, ws in opt.acetyl_weight_scales,
        ss in opt.acetyl_sigma_scales
        idx += 1
        sites = _transform_sites(base_sites, dt, du, us, ws, ss)
        templates = _contrast_templates(sites, coords)
        margins, total_cost = _candidate_margins(by_file, templates)
        delta_bic, sse1, sse2, nminor, nmajor = _kmeans_bic_1d(margins)
        mean_abs_margin = mean(abs.(margins))
        push!(candidates, (delta_bic=delta_bic, mean_abs_margin=mean_abs_margin,
                           total_cost=total_cost, nminor=nminor, nmajor=nmajor,
                           acetyl_t_shift=dt, acetyl_u_shift=du,
                           acetyl_u_scale=us, acetyl_weight_scale=ws,
                           acetyl_sigma_scale=ss, sse1=sse1, sse2=sse2))
        if idx % 100 == 0 || idx == total_candidates
            @printf("[%d/%d] best ΔBIC=%.2f\n", idx, total_candidates, maximum(c.delta_bic for c in candidates))
        end
    end
    sort!(candidates, by=c -> (-c.delta_bic, -c.mean_abs_margin, c.total_cost))
    best = first(candidates)
    best_sites = _transform_sites(base_sites, best.acetyl_t_shift, best.acetyl_u_shift,
                                  best.acetyl_u_scale, best.acetyl_weight_scale,
                                  best.acetyl_sigma_scale)

    mkpath(dirname(opt.report))
    open(opt.report, "w") do io
        println(io, join(["rank", "delta_bic", "mean_abs_margin", "total_cost", "nminor", "nmajor",
                          "acetyl_t_shift", "acetyl_u_shift", "acetyl_u_scale",
                          "acetyl_weight_scale", "acetyl_sigma_scale", "sse1", "sse2"], '\t'))
        for (rank, c) in enumerate(candidates)
            println(io, join([rank, @sprintf("%.8g", c.delta_bic), @sprintf("%.8g", c.mean_abs_margin),
                              @sprintf("%.8g", c.total_cost), c.nminor, c.nmajor,
                              @sprintf("%.8g", c.acetyl_t_shift), @sprintf("%.8g", c.acetyl_u_shift),
                              @sprintf("%.8g", c.acetyl_u_scale), @sprintf("%.8g", c.acetyl_weight_scale),
                              @sprintf("%.8g", c.acetyl_sigma_scale), @sprintf("%.8g", c.sse1),
                              @sprintf("%.8g", c.sse2)], '\t'))
        end
    end
    _write_sites(opt.out_sites, best_sites)

    println("Geometric mold refinement")
    println("  patches:   ", opt.patches)
    println("  prefix:    ", opt.prefix)
    println("  candidates:", total_candidates)
    println("  report:    ", opt.report)
    println("  out sites: ", opt.out_sites)
    @printf("  best: ΔBIC=%.2f mean|margin|=%.4f minority=%d majority=%d dt=%.3f du=%.3f us=%.3f ws=%.3f ss=%.3f\n",
            best.delta_bic, best.mean_abs_margin, best.nminor, best.nmajor,
            best.acetyl_t_shift, best.acetyl_u_shift, best.acetyl_u_scale,
            best.acetyl_weight_scale, best.acetyl_sigma_scale)
end

main()

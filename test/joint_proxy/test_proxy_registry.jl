#!/usr/bin/env julia

# TDD test for the joint-proxy mold registry (Todo 1).
#
# Acceptance criteria (from .omo/plans/joint-proxy-mold-inference-v1.md):
#   - julia --project=. test/joint_proxy/test_proxy_registry.jl passes
#   - asserted ensemble has all 8 (type,parity,mirror) states
#   - finite 81-pixel templates
#   - weights summing to one
#   - deterministic hashes
#   - succeeds from geometric source alone when every optional STM file is absent
#
# QA scenarios:
#   Happy  — real config + local optional files: source list/weights/hashes.
#   Failure — duplicate state, invalid type, non-positive weight, missing
#             required geometric file throw typed/diagnostic errors; optional
#             missing files must NOT throw.
#
# Baseline characterization: the in-memory geometric rendering must match the
# output of test/generate_connected_mold_templates.jl run with the same grid
# params (half_nm=0.32, step_nm=0.08, parity_flip=t, mirror_flip=u, zscore).

using Test
using SHA
using TOML
using Printf

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "proxy_registry.jl"))
using .JointProxyRegistry:
    ProxyEnsemble, ProxySource, ProxyTemplate, ProxyEntry,
    load_registry, render_geometric_templates, load_stm_template_tsv,
    sha256_file, canonical_payload_hash, RegistryConfigError

include(joinpath(@__DIR__, "..", "lib", "joint_proxy", "whole_roi_end_to_end.jl"))

const REPO_ROOT = dirname(dirname(@__DIR__))
const CONFIG_PATH = joinpath(REPO_ROOT, "config", "joint_proxy_molds.toml")
const GEOMETRIC_SITES = joinpath(REPO_ROOT, "templates", "chitosan_geometric_sites.tsv")
const TEST_TMP = mktempdir(; prefix="stmfit-joint-proxy-")
atexit(() -> isdir(TEST_TMP) && rm(TEST_TMP; recursive=true, force=true))

# ─── Helper: write a minimal valid sites TSV ────────────────────────────────
function _write_minimal_sites(path)
    open(path, "w") do io
        println(io, "# minimal test sites")
        println(io, "type\tatom\tt_nm\tu_nm\tweight\tsigma_t_nm\tsigma_u_nm")
        println(io, "0\tring_center\t0.000\t0.000\t1.00\t0.150\t0.115")
        println(io, "0\ta\t0.100\t0.100\t0.50\t0.090\t0.090")
        println(io, "1\tring_center\t0.000\t0.000\t1.00\t0.150\t0.115")
        println(io, "1\tb\t0.100\t0.300\t0.50\t0.080\t0.075")
    end
end

# ─── Helper: write a minimal valid STM unary template TSV (81 px) ────────────
function _write_minimal_stm_template(path; npix=81)
    open(path, "w") do io
        cols = vcat(["name", "type", "parity", "mirror"],
                    [@sprintf("p%03d", i) for i in 1:npix])
        println(io, join(cols, '\t'))
        for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
            name = @sprintf("%s_p%d_m%d", typ == 0 ? "GlcN" : "GlcNAc", parity, mirror)
            vals = [(@sprintf("%.8g", (typ + parity + mirror + i) * 0.1)) for i in 1:npix]
            println(io, join(vcat([name, string(typ), string(parity), string(mirror)], vals), '\t'))
        end
    end
end

function _write_production_sidecar(path, maps_path, template_path;
                                   glcn_sha=repeat("1", 64), glcnac_sha=repeat("2", 64))
    payload = Dict(
        "schema" => "stmfit-qe-mold-provenance-v1",
        "provider" => "stm_dft_v1",
        "sample_bias_ev" => -0.3,
        "height_nm" => 0.5,
        "half_nm" => 0.32,
        "step_nm" => 0.08,
        "cube_units" => "bohr",
        "glcn_cube_sha256" => glcn_sha,
        "glcnac_cube_sha256" => glcnac_sha,
        "maps_sha256" => sha256_file(maps_path),
        "templates_sha256" => sha256_file(template_path),
    )
    open(path, "w") do io
        TOML.print(io, payload; sorted=true)
    end
end

const T2_CC_HALF = 0.32
const T2_CC_STEP = 0.08
const T2_CC_XS = collect(-1.0:0.08:1.0)
const T2_CC_YS = collect(-0.64:0.08:0.64)
const T2_CC_CENTERS = [(-0.35, 0.0), (0.0, 0.0), (0.35, 0.0)]
const T2_CC_SEQUENCES = [Int[0,1,0], Int[1,0,1], Int[0,0,1], Int[1,1,0]]

function _write_t2_cc_cube(path, typ; zstep=0.05)
    nx, ny, nz = 13, 13, 40
    xstep, ystep = 0.08, 0.08
    phase = typ == 0 ? 0.0 : pi
    open(path, "w") do io
        println(io, "synthetic constant-current cube (nm axes)")
        println(io, "T2 temporary provider fixture")
        println(io, "1 0.0 0.0 0.0")
        println(io, "$nx $xstep 0.0 0.0")
        println(io, "$ny 0.0 $ystep 0.0")
        println(io, "$nz 0.0 0.0 $zstep")
        println(io, "1 0.0 0.0 0.0 0.0 0.0 0.0")
        for ix in 1:nx, iy in 1:ny, iz in 1:nz
            x = (ix - 1) * xstep
            y = (iy - 1) * ystep
            z = (iz - 1) * zstep
            lateral = 1.0 + 0.30 * cos(2pi * x / 0.40 + phase) +
                0.10 * cos(2pi * y / 0.35 + phase)
            @printf(io, "%.10e\n", lateral * exp(-z / 0.30))
        end
    end
end

function _t2_cc_config(; shift=0.0, rotation=0.0, blur=0.0)
    signed(limit, step) = Main.WholeRoiDiagnosticConfig.DiagnosticRange(-limit, limit, step)
    positive(limit, step) = Main.WholeRoiDiagnosticConfig.DiagnosticRange(0.0, limit, step)
    shift_coarse = signed(shift, shift == 0 ? 0.01 : shift)
    shift_fine = signed(shift == 0 ? 0.0 : shift / 2, shift == 0 ? 0.01 : shift / 2)
    rot_coarse = signed(rotation, rotation == 0 ? 0.1 : rotation)
    rot_fine = signed(rotation == 0 ? 0.0 : rotation / 2,
                      rotation == 0 ? 0.1 : rotation / 2)
    blur_coarse = positive(blur, blur == 0 ? 0.01 : blur)
    blur_fine = positive(blur == 0 ? 0.0 : blur / 2,
                         blur == 0 ? 0.01 : blur / 2)
    bounds = Main.WholeRoiDiagnosticConfig.DiagnosticParamBounds
    return Main.WholeRoiDiagnosticConfig.DiagnosticConfig(
        bounds(shift_coarse, shift_fine), bounds(shift_coarse, shift_fine),
        bounds(rot_coarse, rot_fine), bounds(blur_coarse, blur_fine),
        1e-6, 1e-9)
end

function _t2_cc_backbone()
    return [sum(exp(-((x-cx)^2 / 0.22^2 + (y-cy)^2 / 0.15^2) / 2)
                for (cx,cy) in T2_CC_CENTERS) for y in T2_CC_YS, x in T2_CC_XS]
end

function _t2_cc_case(set, sequence; shift_t=0.0, shift_u=0.0,
                     rotation=0.0, blur=0.0, noise=0.0,
                     contrast_amp=0.35)
    return synthesize_whole_roi_observation(set, sequence, _t2_cc_backbone(),
        T2_CC_CENTERS, T2_CC_XS, T2_CC_YS;
        half_nm=T2_CC_HALF, step_nm=T2_CC_STEP,
        shift_t_nm=shift_t, shift_u_nm=shift_u,
        rotation_deg=rotation, blur_sigma_nm=blur,
        background=0.25, backbone_amp=0.7, common_amp=1.0,
        contrast_amp=contrast_amp, noise_sigma=noise, seed=20260728)
end

function _t2_cc_run(set, sequence, nuisance)
    if nuisance === :exact
        case, cfg = _t2_cc_case(set, sequence), _t2_cc_config()
    elseif nuisance === :noise
        case, cfg = _t2_cc_case(set, sequence; noise=0.01), _t2_cc_config()
    elseif nuisance === :drift
        case = _t2_cc_case(set, sequence; shift_t=0.02, shift_u=-0.02)
        cfg = _t2_cc_config(shift=0.04)
    elseif nuisance === :blur
        case, cfg = _t2_cc_case(set, sequence; blur=0.01), _t2_cc_config(blur=0.02)
    elseif nuisance === :frame
        case = _t2_cc_case(set, sequence; rotation=0.25)
        cfg = _t2_cc_config(rotation=0.5)
    else
        error("unknown T2 nuisance: $nuisance")
    end
    return run_whole_roi_common_contrast(case, set;
        config=cfg, half_nm=T2_CC_HALF, step_nm=T2_CC_STEP)
end

function _t2_cc_entry(template_path, height)
    templates = load_stm_template_tsv(template_path; npix=81)
    source = ProxySource("stm_dft_cc_diag", "cc_h$(round(Int, height * 100))",
        template_path, sha256_file(template_path), height, 1.0, true)
    return ProxyEntry(source, templates)
end

function _t2_swapped_entry(entry)
    templates = [ProxyTemplate(1 - t.type, t.parity, t.mirror, copy(t.pixels))
                 for t in entry.templates]
    return ProxyEntry(entry.source, templates)
end

function _t2_zero_contrast_entry(entry)
    base = Dict((t.parity, t.mirror) => t.pixels for t in entry.templates if t.type == 0)
    templates = ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        push!(templates, ProxyTemplate(typ, parity, mirror, copy(base[(parity, mirror)])))
    end
    return ProxyEntry(entry.source, templates)
end

@testset "proxy_registry: Todo 1" begin

    @testset "baseline characterization: in-memory render matches generator" begin
        # The registry must reuse the EXACT algorithm from
        # generate_connected_mold_templates.jl. We render in-memory and compare
        # against a fresh generator run with matching params.
        coords = collect(-0.32:0.08:0.32)
        @test length(coords) == 9
        tmpls = render_geometric_templates(GEOMETRIC_SITES;
                                           half_nm=0.32, step_nm=0.08,
                                           parity_flip="t", mirror_flip="u",
                                           normalize="zscore")
        @test length(tmpls) == 8
        @test all(length(t.pixels) == 81 for t in tmpls)
        states = Set((t.type, t.parity, t.mirror) for t in tmpls)
        @test length(states) == 8
        for typ in (0, 1), par in (0, 1), mir in (0, 1)
            @test (typ, par, mir) in states
        end
        @test all(all(isfinite, t.pixels) for t in tmpls)

        # Compare against the generator's own output (fresh run to temp file).
        # The generator is the canonical existing behavior we reuse.
        gen_out = joinpath(TEST_TMP, "baseline_gen_molds.tsv")
        run(`julia --project=. test/generate_connected_mold_templates.jl
             --atoms $(GEOMETRIC_SITES) --out $(gen_out)`)
        lines = readlines(gen_out)
        @test length(lines) == 9  # header + 8 rows
        header = split(lines[1], '\t')
        pix_cols = [c for c in header if occursin(r"^p\d+$", c)]
        @test length(pix_cols) == 81
        gen_map = Dict{Tuple{Int,Int,Int},Vector{Float64}}()
        for line in lines[2:end]
            vals = split(line, '\t')
            typ = parse(Int, vals[2])
            par = parse(Int, vals[3])
            mir = parse(Int, vals[4])
            pix = [parse(Float64, v) for v in vals[5:end]]
            gen_map[(typ, par, mir)] = pix
        end
        for t in tmpls
            gen_pix = gen_map[(t.type, t.parity, t.mirror)]
            @test all(isapprox.(t.pixels, gen_pix; atol=1e-7))
        end
    end

    @testset "happy: real config loads with all available sources" begin
        ens = load_registry(CONFIG_PATH)
        @test ens.grid_n == 9
        @test ens.npix == 81
        # Every entry has all 8 states with finite 81-pixel templates.
        for entry in ens.entries
            @test length(entry.templates) == 8
            states = Set((t.type, t.parity, t.mirror) for t in entry.templates)
            @test length(states) == 8
            for typ in (0, 1), par in (0, 1), mir in (0, 1)
                @test (typ, par, mir) in states
            end
            @test all(all(isfinite, t.pixels) for t in entry.templates)
            @test all(length(t.pixels) == 81 for t in entry.templates)
        end
        # Geometric is always present (required).
        geom = [e for e in ens.entries if e.source.family == "geometric"]
        @test length(geom) == 1
        @test geom[1].source.present
        # Weights sum to one.
        w = [e.source.weight for e in ens.entries]
        @test all(w .> 0)
        @test isapprox(sum(w), 1.0; atol=1e-12)
        # SHA-256 provenance present for every consumed source.
        for entry in ens.entries
            @test occursin(r"^[0-9a-f]{64}$", entry.source.sha256)
        end
        @test occursin(r"^[0-9a-f]{64}$", ens.payload_sha256)
        # Deterministic: reload and compare payload hash.
        ens2 = load_registry(CONFIG_PATH)
        @test ens.payload_sha256 == ens2.payload_sha256
    end

    @testset "geometric-only success when all STM files absent" begin
        # Build a temp config pointing geometric at the real sites and STM at
        # nonexistent paths. Must succeed with geometric-only weight = 1.0.
        tmp_sites = joinpath(TEST_TMP, "sites_geo_only.tsv")
        _write_minimal_sites(tmp_sites)
        tmp_stm = joinpath(TEST_TMP, "missing_h050.tsv")
        isfile(tmp_stm) && rm(tmp_stm)
        cfg = Dict(
            "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
            "geometric" => Dict("sites_tsv" => tmp_sites,
                                "parity_flip" => "t", "mirror_flip" => "u",
                                "normalize" => "zscore"),
            "stm_prelim" => Dict("family_weight" => 0.5,
                "sources" => [Dict("height_nm" => 0.50, "path" => tmp_stm,
                                   "weight" => 3.0)]),
            "provenance" => Dict("include_bond_templates" => false),
        )
        tmp_cfg = joinpath(TEST_TMP, "geo_only.toml")
        open(tmp_cfg, "w") do io
            TOML.print(io, cfg)
        end
        ens = load_registry(tmp_cfg)
        @test length(ens.entries) == 1
        @test ens.entries[1].source.family == "geometric"
        @test isapprox(ens.entries[1].source.weight, 1.0; atol=1e-12)
        @test length(ens.entries[1].templates) == 8
        # Audit must mention the missing optional file.
        @test any(occursin("missing", lowercase(a)) for a in ens.audit)
    end

    @testset "deterministic hashes across reloads" begin
        ens_a = load_registry(CONFIG_PATH)
        ens_b = load_registry(CONFIG_PATH)
        @test ens_a.payload_sha256 == ens_b.payload_sha256
        for (ea, eb) in zip(ens_a.entries, ens_b.entries)
            @test ea.source.sha256 == eb.source.sha256
            @test ea.source.weight == eb.source.weight
            for (ta, tb) in zip(ea.templates, eb.templates)
                @test ta.type == tb.type
                @test ta.parity == tb.parity
                @test ta.mirror == tb.mirror
                @test all(ta.pixels .== tb.pixels)
            end
        end
    end

    @testset "normalized priors: STM family renormalizes when some missing" begin
        # Two STM sources, one missing. Present one must absorb the STM family
        # mass (0.5). Geometric stays at 0.5.
        tmp_sites = joinpath(TEST_TMP, "sites_renorm.tsv")
        _write_minimal_sites(tmp_sites)
        tmp_stm_present = joinpath(TEST_TMP, "stm_present.tsv")
        _write_minimal_stm_template(tmp_stm_present)
        tmp_stm_missing = joinpath(TEST_TMP, "stm_missing.tsv")
        isfile(tmp_stm_missing) && rm(tmp_stm_missing)
        cfg = Dict(
            "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
            "geometric" => Dict("sites_tsv" => tmp_sites,
                                "parity_flip" => "t", "mirror_flip" => "u",
                                "normalize" => "zscore"),
            "stm_prelim" => Dict("family_weight" => 0.5,
                "sources" => [Dict("height_nm" => 0.45, "path" => tmp_stm_present,
                                   "weight" => 2.0),
                              Dict("height_nm" => 0.55, "path" => tmp_stm_missing,
                                   "weight" => 2.0)]),
            "provenance" => Dict("include_bond_templates" => false),
        )
        tmp_cfg = joinpath(TEST_TMP, "renorm.toml")
        open(tmp_cfg, "w") do io
            TOML.print(io, cfg)
        end
        ens = load_registry(tmp_cfg)
        # Only geometric + one STM entry present.
        @test length(ens.entries) == 2
        geom_w = [e.source.weight for e in ens.entries if e.source.family == "geometric"][1]
        stm_w = [e.source.weight for e in ens.entries if e.source.family == "stm_prelim"][1]
        @test isapprox(geom_w, 0.5; atol=1e-12)
        @test isapprox(stm_w, 0.5; atol=1e-12)  # renormalized within family
        @test isapprox(geom_w + stm_w, 1.0; atol=1e-12)
    end

    @testset "production DFT provider replaces preliminary evidence with pinned provenance" begin
        tmp = mktempdir()
        sites = joinpath(tmp, "sites.tsv")
        templates = joinpath(tmp, "dft_templates.tsv")
        maps = joinpath(tmp, "dft_maps.tsv")
        sidecar = joinpath(tmp, "dft_provenance.toml")
        prelim = joinpath(tmp, "prelim.tsv")
        _write_minimal_sites(sites)
        _write_minimal_stm_template(templates)
        _write_minimal_stm_template(prelim)
        write(maps, "fixture maps\n")
        _write_production_sidecar(sidecar, maps, templates)
        cfg = Dict(
            "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
            "geometric" => Dict("sites_tsv" => sites, "parity_flip" => "t",
                                "mirror_flip" => "u", "normalize" => "zscore"),
            "stm_dft_v1" => Dict(
                "family_weight" => 0.5, "path" => templates,
                "sha256" => sha256_file(templates), "maps_path" => maps,
                "maps_sha256" => sha256_file(maps), "provenance_path" => sidecar,
                "provenance_sha256" => sha256_file(sidecar), "height_nm" => 0.5,
                "sample_bias_ev" => -0.3, "cube_units" => "bohr",
                "glcn_cube_sha256" => repeat("1", 64),
                "glcnac_cube_sha256" => repeat("2", 64)),
            "stm_prelim" => Dict("enabled" => false, "family_weight" => 0.5,
                "sources" => [Dict("height_nm" => 0.5, "path" => prelim, "weight" => 1.0)]),
            "provenance" => Dict("include_bond_templates" => false),
        )
        config = joinpath(tmp, "production.toml")
        open(config, "w") do io; TOML.print(io, cfg); end
        ens = load_registry(config)
        @test Set(e.source.family for e in ens.entries) == Set(["geometric", "stm_dft_v1"])
        @test all(isapprox(e.source.weight, 0.5; atol=1e-12) for e in ens.entries)
        @test !any(e.source.path == prelim for e in ens.entries)
        @test any(occursin("disabled", lowercase(line)) for line in ens.audit)

        cfg["stm_dft_v1"]["sha256"] = repeat("0", 64)
        open(config, "w") do io; TOML.print(io, cfg); end
        @test_throws RegistryConfigError load_registry(config)
    end

    @testset "temporary constant-current provider remains diagnostic-only" begin
        config_sha_before = sha256_file(CONFIG_PATH)
        registry_before = load_registry(CONFIG_PATH)
        production_before = only(e for e in registry_before.entries
                                 if e.source.family == "stm_dft_v1")
        production_payload_before = canonical_payload_hash(
            [production_before], registry_before.grid_half_nm,
            registry_before.grid_step_nm, registry_before.npix)
        @test !any(e.source.family == "stm_dft_cc_diag" for e in registry_before.entries)

        mktempdir() do tmp
            cube0 = joinpath(tmp, "glcn.cube")
            cube1 = joinpath(tmp, "glcnac.cube")
            prefix = joinpath(tmp, "cc_diag")
            _write_t2_cc_cube(cube0, 0)
            _write_t2_cc_cube(cube1, 1)
            builder = joinpath(REPO_ROOT, "test", "build_constant_current_stm_maps.jl")

            frame0 = joinpath(tmp, "frame0.tsv")
            frame1 = joinpath(tmp, "frame1.tsv")
            write(frame0, "origin_nm\t0.48,0.48,0.0\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
            write(frame1, "origin_nm\t0.49,0.48,0.0\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
            mixed_prefix = joinpath(tmp, "mixed_frame")
            run(pipeline(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(builder)
                --cube0 $(cube0) --cube1 $(cube1) --cube-units nm
                --frame0 $(frame0) --frame1 $(frame1)
                --out-prefix $(mixed_prefix)`; stdout=devnull, stderr=devnull))
            mixed_provenance = TOML.parsefile(string(mixed_prefix, ".provenance.toml"))
            @test mixed_provenance["isovalue_scan_intervals"] == 1024
            @test !haskey(mixed_provenance, "reference_plane")
            mixed_frames = Dict(frame["type"] => frame for frame in mixed_provenance["type_frames"])
            @test Set(keys(mixed_frames)) == Set([0, 1])
            @test mixed_frames[0]["origin_nm"] == [0.48, 0.48, 0.0]
            @test mixed_frames[1]["origin_nm"] == [0.49, 0.48, 0.0]

            periodic_prefix = joinpath(tmp, "periodic_unbound")
            @test_throws ProcessFailedException run(pipeline(
                `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(builder)
                    --cube0 $(cube0) --cube1 $(cube1) --cube-units nm
                    --origin 0.48,0.48,0.0 --t-axis 1,0,0 --u-axis 0,1,0
                    --periodic-axes x --out-prefix $(periodic_prefix)`;
                stdout=devnull, stderr=devnull))

            invalid_scan_prefix = joinpath(tmp, "invalid_scan_intervals")
            @test_throws ProcessFailedException run(pipeline(
                `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(builder)
                    --cube0 $(cube0) --cube1 $(cube1) --cube-units nm
                    --origin 0.48,0.48,0.0 --t-axis 1,0,0 --u-axis 0,1,0
                    --isovalue-scan-intervals 0 --out-prefix $(invalid_scan_prefix)`;
                stdout=devnull, stderr=devnull))

            cube1_mismatched_z = joinpath(tmp, "glcnac-mismatched-z.cube")
            _write_t2_cc_cube(cube1_mismatched_z, 1; zstep=0.04)
            spacing_prefix = joinpath(tmp, "mismatched_z_spacing")
            @test_throws ProcessFailedException run(pipeline(
                `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(builder)
                    --cube0 $(cube0) --cube1 $(cube1_mismatched_z) --cube-units nm
                    --origin 0.48,0.48,0.0 --t-axis 1,0,0 --u-axis 0,1,0
                    --out-prefix $(spacing_prefix)`; stdout=devnull, stderr=devnull))

            victim = joinpath(tmp, "map-victim.tsv")
            nominal_map = string(prefix, "_h050.tsv")
            write(victim, "sentinel\n")
            symlink(victim, nominal_map)
            run(pipeline(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(builder)
                --cube0 $(cube0) --cube1 $(cube1) --cube-units nm
                --origin 0.48,0.48,0.0 --t-axis 1,0,0 --u-axis 0,1,0
                --isovalue-scan-intervals 64 --out-prefix $(prefix)`; stdout=devnull))
            @test read(victim, String) == "sentinel\n"
            @test !islink(nominal_map)

            provenance = TOML.parsefile(string(prefix, ".provenance.toml"))
            @test provenance["provider"] == "stm_dft_cc_diag"
            @test provenance["observable"] == "constant-current"
            @test provenance["isovalue_scan_intervals"] == 64
            @test length(provenance["type_frames"]) == 2
            type_isovalues = Dict(item["type"] => item for item in provenance["type_isovalues"])
            @test Set(keys(type_isovalues)) == Set([0, 1])
            @test all(isfinite(Float64(item["isovalue"])) for item in values(type_isovalues))
            @test !haskey(provenance, "reference_plane")
            @test !haskey(provenance, "isovalue")
            @test provenance["bracket_heights_nm"] ==
                  [0.40, 0.45, 0.50, 0.55, 0.60]
            @test provenance["glcn_cube_sha256"] == sha256_file(cube0)
            @test provenance["glcnac_cube_sha256"] == sha256_file(cube1)
            @test length(provenance["bracket_artifacts"]) == 4

            importer = joinpath(REPO_ROOT, "test", "import_stm_mold_maps.jl")
            provider_hashes = String[]
            for height in (0.40, 0.45, 0.50, 0.55, 0.60)
                tag = @sprintf("%03d", round(Int, 100height))
                map_path = string(prefix, "_h", tag, ".tsv")
                mask_path = string(prefix, "_h", tag, ".mask.tsv")
                template_path = string(prefix, "_h", tag, "_templates.tsv")
                run(pipeline(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(importer)
                    --maps $(map_path) --out $(template_path)
                    --half-nm 0.32 --step-nm 0.08 --normalize zscore`;
                    stdout=devnull))

                entry = _t2_cc_entry(template_path, height)
                push!(provider_hashes, entry.source.sha256)
                @test entry.source.family == "stm_dft_cc_diag"
                @test entry.source.family != "stm_dft_v1"
                @test entry.source.height_nm == height
                @test entry.source.sha256 == sha256_file(template_path)
                @test length(entry.templates) == 8
                @test all(length(t.pixels) == 81 for t in entry.templates)
                @test all(all(isfinite, t.pixels) for t in entry.templates)

                set = derive_common_contrast(entry)
                for sequence in T2_CC_SEQUENCES,
                    nuisance in (:exact, :noise, :drift, :blur, :frame)
                    result = _t2_cc_run(set, sequence, nuisance)
                    @test !result.scoring.abstain
                    @test result.scoring.best.sequence == sequence
                end

                common_case = _t2_cc_case(set, T2_CC_SEQUENCES[1]; contrast_amp=0.0)
                common_result = run_whole_roi_common_contrast(common_case, set;
                    config=_t2_cc_config(), half_nm=T2_CC_HALF, step_nm=T2_CC_STEP)
                @test common_result.scoring.abstain

                zero_entry = _t2_zero_contrast_entry(entry)
                zero_set = derive_common_contrast(zero_entry)
                zero_result = _t2_cc_run(zero_set, T2_CC_SEQUENCES[1], :exact)
                @test zero_result.scoring.abstain

                sequence = T2_CC_SEQUENCES[2]
                observed = _t2_cc_case(set, sequence)
                swapped_set = derive_common_contrast(_t2_swapped_entry(entry))
                swapped_result = run_whole_roi_common_contrast(observed, swapped_set;
                    config=_t2_cc_config(), half_nm=T2_CC_HALF, step_nm=T2_CC_STEP)
                @test !swapped_result.scoring.abstain
                @test swapped_result.scoring.best.sequence == 1 .- sequence

                @test isfile(map_path)
                @test isfile(mask_path)
            end
            @test length(unique(provider_hashes)) == 5

            nominal_mask = string(prefix, "_h050.mask.tsv")
            @test provenance["maps_sha256"] == sha256_file(nominal_map)
            @test provenance["invalid_mask_sha256"] == sha256_file(nominal_mask)
            nominal_template = string(prefix, "_h050_templates.tsv")
            nominal_binding = string(nominal_template, ".provenance.toml")
            run(pipeline(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(importer)
                --maps $(nominal_map) --out $(nominal_template)
                --half-nm 0.32 --step-nm 0.08 --normalize zscore
                --source-provenance $(string(prefix, ".provenance.toml"))
                --mold-provenance-out $(nominal_binding)`;
                stdout=devnull))
            mold_binding = TOML.parsefile(nominal_binding)
            @test mold_binding["schema"] == "stmfit-constant-current-mold-binding-v1"
            @test mold_binding["source_maps_sha256"] == provenance["maps_sha256"]
            @test mold_binding["molds_sha256"] == sha256_file(nominal_template)
            for artifact in provenance["bracket_artifacts"]
                tag = @sprintf("%03d", round(Int, 100artifact["height_nm"]))
                @test artifact["map_sha256"] == sha256_file(
                    string(prefix, "_h", tag, ".tsv"))
                @test artifact["map_path"] == string(prefix, "_h", tag, ".tsv")
                @test artifact["invalid_mask_sha256"] == sha256_file(
                    string(prefix, "_h", tag, ".mask.tsv"))
                @test artifact["invalid_mask_path"] == string(prefix, "_h", tag, ".mask.tsv")
            end
        end

        registry_after = load_registry(CONFIG_PATH)
        production_after = only(e for e in registry_after.entries
                                if e.source.family == "stm_dft_v1")
        production_payload_after = canonical_payload_hash(
            [production_after], registry_after.grid_half_nm,
            registry_after.grid_step_nm, registry_after.npix)
        @test sha256_file(CONFIG_PATH) == config_sha_before
        @test registry_after.payload_sha256 == registry_before.payload_sha256
        @test production_after.source.sha256 == production_before.source.sha256
        @test production_after.source.path == production_before.source.path
        @test production_after.source.height_nm == production_before.source.height_nm
        @test production_after.source.weight == production_before.source.weight
        @test production_payload_after == production_payload_before
        @test !any(e.source.family == "stm_dft_cc_diag" for e in registry_after.entries)
    end

    # ─── Failure scenarios (strict errors for malformed required input) ─────
    @testset "failure: missing required geometric sites file" begin
        tmp_cfg = joinpath(TEST_TMP, "missing_geom.toml")
        cfg = Dict(
            "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
            "geometric" => Dict("sites_tsv" => joinpath(TEST_TMP, "nonexistent.tsv"),
                                "parity_flip" => "t", "mirror_flip" => "u",
                                "normalize" => "zscore"),
            "stm_prelim" => Dict("family_weight" => 0.5, "sources" => []),
            "provenance" => Dict("include_bond_templates" => false),
        )
        open(tmp_cfg, "w") do io; TOML.print(io, cfg); end
        @test_throws RegistryConfigError load_registry(tmp_cfg)
    end

    @testset "failure: duplicate (type,parity,mirror) state in STM TSV" begin
        tmp_sites = joinpath(TEST_TMP, "sites_dup.tsv")
        _write_minimal_sites(tmp_sites)
        tmp_stm = joinpath(TEST_TMP, "stm_dup.tsv")
        open(tmp_stm, "w") do io
            cols = vcat(["name", "type", "parity", "mirror"],
                        [@sprintf("p%03d", i) for i in 1:81])
            println(io, join(cols, '\t'))
            # Duplicate (0,0,0) row.
            for _ in 1:2
                vals = [@sprintf("%.8g", i * 0.1) for i in 1:81]
                println(io, join(vcat(["GlcN_p0_m0", "0", "0", "0"], vals), '\t'))
            end
            for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
                (typ, parity, mirror) == (0, 0, 0) && continue
                vals = [@sprintf("%.8g", i * 0.1) for i in 1:81]
                name = @sprintf("%s_p%d_m%d", typ == 0 ? "GlcN" : "GlcNAc", parity, mirror)
                println(io, join(vcat([name, string(typ), string(parity), string(mirror)], vals), '\t'))
            end
        end
        @test_throws RegistryConfigError load_stm_template_tsv(tmp_stm; npix=81)
    end

    @testset "failure: invalid type in STM TSV" begin
        tmp_stm = joinpath(TEST_TMP, "stm_badtype.tsv")
        open(tmp_stm, "w") do io
            cols = vcat(["name", "type", "parity", "mirror"],
                        [@sprintf("p%03d", i) for i in 1:81])
            println(io, join(cols, '\t'))
            # type=2 is invalid.
            for typ in (0, 1, 2), parity in (0, 1), mirror in (0, 1)
                typ in (0, 1) && (typ, parity, mirror) != (0, 0, 0) && continue
                typ == 2 && (typ, parity, mirror) != (0, 0, 0) && continue
                vals = [@sprintf("%.8g", i * 0.1) for i in 1:81]
                name = @sprintf("T%d_p%d_m%d", typ, parity, mirror)
                println(io, join(vcat([name, string(typ), string(parity), string(mirror)], vals), '\t'))
            end
        end
        @test_throws RegistryConfigError load_stm_template_tsv(tmp_stm; npix=81)
    end

    @testset "failure: non-positive weight in config" begin
        tmp_sites = joinpath(TEST_TMP, "sites_negw.tsv")
        _write_minimal_sites(tmp_sites)
        tmp_stm = joinpath(TEST_TMP, "stm_negw.tsv")
        _write_minimal_stm_template(tmp_stm)
        cfg = Dict(
            "grid" => Dict("half_nm" => 0.32, "step_nm" => 0.08),
            "geometric" => Dict("sites_tsv" => tmp_sites,
                                "parity_flip" => "t", "mirror_flip" => "u",
                                "normalize" => "zscore"),
            "stm_prelim" => Dict("family_weight" => 0.5,
                "sources" => [Dict("height_nm" => 0.50, "path" => tmp_stm,
                                   "weight" => -1.0)]),  # negative
            "provenance" => Dict("include_bond_templates" => false),
        )
        tmp_cfg = joinpath(TEST_TMP, "negw.toml")
        open(tmp_cfg, "w") do io; TOML.print(io, cfg); end
        @test_throws RegistryConfigError load_registry(tmp_cfg)
    end

    @testset "failure: missing required states in STM TSV" begin
        tmp_stm = joinpath(TEST_TMP, "stm_partial.tsv")
        open(tmp_stm, "w") do io
            cols = vcat(["name", "type", "parity", "mirror"],
                        [@sprintf("p%03d", i) for i in 1:81])
            println(io, join(cols, '\t'))
            # Only write 3 of 8 states.
            for (typ, parity, mirror) in ((0,0,0), (0,1,0), (1,0,1))
                vals = [@sprintf("%.8g", i * 0.1) for i in 1:81]
                name = @sprintf("%s_p%d_m%d", typ == 0 ? "GlcN" : "GlcNAc", parity, mirror)
                println(io, join(vcat([name, string(typ), string(parity), string(mirror)], vals), '\t'))
            end
        end
        @test_throws RegistryConfigError load_stm_template_tsv(tmp_stm; npix=81)
    end

    @testset "failure: wrong pixel count in STM TSV" begin
        tmp_stm = joinpath(TEST_TMP, "stm_wrongpix.tsv")
        # 64 pixels instead of 81.
        _write_minimal_stm_template(tmp_stm; npix=64)
        @test_throws RegistryConfigError load_stm_template_tsv(tmp_stm; npix=81)
    end

    @testset "no bond templates loaded" begin
        ens = load_registry(CONFIG_PATH)
        # Every entry is a unary template set; no pair/bond entries.
        @test all(length(e.templates) == 8 for e in ens.entries)
        # Provenance flag recorded.
        @test occursin("bond", lowercase(join(ens.audit, " ")))
    end

    @testset "sha256_file helper" begin
        tmp = joinpath(TEST_TMP, "sha_test.txt")
        open(tmp, "w") do io; write(io, "hello\n"); end
        h = sha256_file(tmp)
        @test occursin(r"^[0-9a-f]{64}$", h)
        # Deterministic.
        @test sha256_file(tmp) == h
    end
end

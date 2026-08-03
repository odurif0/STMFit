#!/usr/bin/env julia

# Focused tests for the whole-ROI diagnostic registration config parser/validator.
# Part of the common-template-frozen-contrast-v2 diagnostic (Todo 1). Diagnostic
# only — does not touch production selection, calibration, or the mold registry.

using Test
using TOML

include(joinpath(dirname(@__DIR__), "lib", "joint_proxy", "whole_roi_diagnostic_config.jl"))
using .WholeRoiDiagnosticConfig: DiagnosticConfig, DiagnosticRange,
    DiagnosticParamBounds, parse_diagnostic_config, load_diagnostic_config,
    search_values, DEFAULT_DIAGNOSTIC_CONFIG_PATH

# ---- fixtures -------------------------------------------------------------

_good_range(mn, mx, st) = Dict("min" => mn, "max" => mx, "step" => st)

# A complete, well-formed config dict. Tests mutate a deepcopy to introduce
# exactly one defect, so every rejection is attributable to a single cause.
function _good_config()
    Dict(
        "registration" => Dict(
            "coarse" => Dict(
                "shift_t_nm"    => _good_range(-0.30, 0.30, 0.10),
                "shift_u_nm"    => _good_range(-0.20, 0.20, 0.10),
                "rotation_deg"  => _good_range(-5.0, 5.0, 1.0),
                "blur_sigma_nm" => _good_range(0.0, 0.10, 0.05),
            ),
            "fine" => Dict(
                "shift_t_nm"    => _good_range(-0.05, 0.05, 0.02),
                "shift_u_nm"    => _good_range(-0.04, 0.04, 0.02),
                "rotation_deg"  => _good_range(-1.0, 1.0, 0.25),
                "blur_sigma_nm" => _good_range(0.0, 0.03, 0.01),
            ),
        ),
        "numerical_tie" => Dict("rtol" => 1e-6, "atol" => 1e-9),
    )
end

# Deep-copy a Dict-of-Dict config (TOML.load gives JSON-like Dicts; copy is shallow
# at depth, so we walk one level for the mutation targets used below).
_copy(d::AbstractDict) = Dict(k => _copy(v) for (k, v) in d)
_copy(v) = v

"name the exception message produced by a failing parse, for evidence."
function _err_msg(f)
    try
        f()
        return "<no error thrown>"
    catch e
        return e isa Exception ? sprint(showerror, e) : string(e)
    end
end

@testset "whole-ROI diagnostic config: happy path" begin
    cfg = parse_diagnostic_config(_good_config())
    @test cfg isa DiagnosticConfig
    @test cfg.shift_t_nm.coarse isa DiagnosticRange
    @test cfg.shift_t_nm.coarse.min == -0.30
    @test cfg.shift_t_nm.coarse.max == 0.30
    @test cfg.shift_t_nm.coarse.step == 0.10
    @test cfg.shift_t_nm.fine.min == -0.05
    @test cfg.shift_u_nm.coarse.min == -0.20
    @test cfg.rotation_deg.fine.step == 0.25
    @test cfg.blur_sigma_nm.fine.min == 0.0
    @test cfg.blur_sigma_nm.fine.max == 0.03
    @test cfg.tie_rtol == 1e-6
    @test cfg.tie_atol == 1e-9
    # immutable: no field setters allowed
    @test_throws ErrorException (cfg.tie_rtol = 0.0)
end

@testset "parses the on-disk tracked config" begin
    @test isfile(DEFAULT_DIAGNOSTIC_CONFIG_PATH)
    on_disk = load_diagnostic_config(DEFAULT_DIAGNOSTIC_CONFIG_PATH)
    @test on_disk isa DiagnosticConfig
    # The tracked file must satisfy the same contract as the fixture.
    @test on_disk == parse_diagnostic_config(_good_config())
    # translation/rotation bounds symmetric; blur non-negative
    for b in (on_disk.shift_t_nm, on_disk.shift_u_nm, on_disk.rotation_deg)
        @test b.coarse.min == -b.coarse.max
        @test b.fine.min == -b.fine.max
    end
    @test on_disk.blur_sigma_nm.coarse.min >= 0
    @test on_disk.blur_sigma_nm.fine.min >= 0
end

@testset "search_values is deterministic and endpoint-closed" begin
    r = DiagnosticRange(-0.30, 0.30, 0.10)
    v1 = search_values(r)
    v2 = search_values(r)
    @test v1 == v2                         # determinism
    @test first(v1) == -0.30 && last(v1) == 0.30   # endpoints included
    @test length(v1) == 7                  # -0.30:0.10:0.30
    @test allunique(v1)
    # degenerate (frozen-DOF) range yields a single point, not an error
    @test search_values(DiagnosticRange(0.0, 0.0, 0.10)) == [0.0]
end

@testset "reject missing sections / parameters" begin
    for key in ("registration",)
        bad = _copy(_good_config()); pop!(bad, key)
        @test_throws ArgumentError parse_diagnostic_config(bad)
    end
    for key in ("coarse", "fine")
        bad = _copy(_good_config()); pop!(bad["registration"], key)
        @test_throws ArgumentError parse_diagnostic_config(bad)
    end
    bad = _copy(_good_config()); pop!(bad, "numerical_tie")
    @test_throws ArgumentError parse_diagnostic_config(bad)
    # missing parameter block within a stage
    bad = _copy(_good_config()); pop!(bad["registration"]["coarse"], "shift_t_nm")
    @test_throws ArgumentError parse_diagnostic_config(bad)
    # missing range key
    bad = _copy(_good_config()); pop!(bad["registration"]["coarse"]["shift_t_nm"], "step")
    @test_throws ArgumentError parse_diagnostic_config(bad)
    # missing tie tolerance
    bad = _copy(_good_config()); pop!(bad["numerical_tie"], "atol")
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "reject asymmetric / negative / inverted ranges" begin
    # asymmetric translation bounds
    bad = _copy(_good_config())
    bad["registration"]["coarse"]["shift_t_nm"]["max"] = 0.29
    @test_throws ArgumentError parse_diagnostic_config(bad)
    @test occursin("shift_t_nm", _err_msg(() -> parse_diagnostic_config(bad)))
    # negative blur sigma
    bad = _copy(_good_config())
    bad["registration"]["fine"]["blur_sigma_nm"]["min"] = -0.01
    @test_throws ArgumentError parse_diagnostic_config(bad)
    @test occursin("blur_sigma_nm", _err_msg(() -> parse_diagnostic_config(bad)))
    # inverted range (min > max)
    bad = _copy(_good_config())
    bad["registration"]["coarse"]["rotation_deg"]["min"] = 5.0
    bad["registration"]["coarse"]["rotation_deg"]["max"] = -5.0
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "reject nonpositive / non-finite steps and values" begin
    for bad_step in (0.0, -0.10)
        bad = _copy(_good_config())
        bad["registration"]["coarse"]["shift_t_nm"]["step"] = bad_step
        @test_throws ArgumentError parse_diagnostic_config(bad)
    end
    @test occursin("step", _err_msg(() -> parse_diagnostic_config(let d = _copy(_good_config())
        d["registration"]["coarse"]["shift_t_nm"]["step"] = 0.0; d end)))
    # non-finite value
    bad = _copy(_good_config())
    bad["registration"]["fine"]["shift_u_nm"]["max"] = Inf
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "reject non-deterministic grid (width not an integer multiple of step)" begin
    bad = _copy(_good_config())
    bad["registration"]["coarse"]["shift_t_nm"]["step"] = 0.07   # 0.60/0.07 not integer
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "reject fine range wider / coarser than coarse" begin
    # fine span > coarse span
    bad = _copy(_good_config())
    bad["registration"]["coarse"]["shift_t_nm"] = _good_range(-0.05, 0.05, 0.02)  # span 0.10
    bad["registration"]["fine"]["shift_t_nm"]   = _good_range(-0.10, 0.10, 0.05)  # span 0.20
    @test_throws ArgumentError parse_diagnostic_config(bad)
    @test occursin("shift_t_nm", _err_msg(() -> parse_diagnostic_config(bad)))
    # fine step > coarse step (spans kept legal)
    bad = _copy(_good_config())
    bad["registration"]["coarse"]["shift_t_nm"] = _good_range(-0.10, 0.10, 0.02)  # step 0.02
    bad["registration"]["fine"]["shift_t_nm"]   = _good_range(-0.05, 0.05, 0.10)  # step 0.10
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "reject invalid numerical tie tolerances" begin
    for kw in ("rtol", "atol")
        bad = _copy(_good_config()); bad["numerical_tie"][kw] = -1e-3
        @test_throws ArgumentError parse_diagnostic_config(bad)
        bad = _copy(_good_config()); bad["numerical_tie"][kw] = NaN
        @test_throws ArgumentError parse_diagnostic_config(bad)
    end
    # both zero: no scale-aware tie tolerance at all
    bad = _copy(_good_config()); bad["numerical_tie"]["rtol"] = 0.0; bad["numerical_tie"]["atol"] = 0.0
    @test_throws ArgumentError parse_diagnostic_config(bad)
end

@testset "no misleading success: missing data throws, never returns a default" begin
    @test_throws ArgumentError parse_diagnostic_config(Dict())            # empty
    @test_throws ArgumentError parse_diagnostic_config(Dict("numerical_tie" => Dict("rtol" => 1e-6, "atol" => 1e-9)))
    @test_throws ArgumentError parse_diagnostic_config(Dict("registration" => Dict()))  # stages absent
    # a syntactically-legal TOML lacking the required sections must not silently pass
    tmp = mkpath(joinpath(tempdir(), "wr_diag_cfg_test"))
    path = joinpath(tmp, "no_sections.toml")
    write(path, "[unrelated]\nkey = 1\n")
    @test_throws ArgumentError load_diagnostic_config(path)
    @test_throws ArgumentError load_diagnostic_config(joinpath(tmp, "does_not_exist.toml"))
    rm(tmp; force=true, recursive=true)
end

@testset "no stale state: parser holds no cache" begin
    d = _copy(_good_config())
    first = parse_diagnostic_config(d)
    @test first.tie_rtol == 1e-6
    d["numerical_tie"]["rtol"] = 7e-4            # mutate after first parse
    again = parse_diagnostic_config(d)
    @test again.tie_rtol == 7e-4                 # reflects latest input, no memoization
    @test first.tie_rtol == 1e-6                 # prior immutable result untouched
end

@testset "determinism: repeated parses are byte-equal" begin
    a = parse_diagnostic_config(_good_config())
    b = parse_diagnostic_config(_good_config())
    @test a == b
    @test search_values(a.rotation_deg.fine) == search_values(b.rotation_deg.fine)
end

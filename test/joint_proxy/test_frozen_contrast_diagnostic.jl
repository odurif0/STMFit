#!/usr/bin/env julia

using Test

include(joinpath(dirname(@__DIR__), "plot_frozen_contrast_diagnostic.jl"))

const RD_REG = Main.JointProxyRegistry
const RD_DC = Main.WholeRoiDiagnosticConfig

function rd_entry(; identical=false)
    M0 = zeros(3, 3); M0[2, 1] = 1.0; M0[1, 2] = 0.25
    M1 = identical ? copy(M0) : zeros(3, 3)
    identical || (M1[2, 3] = 1.1; M1[3, 2] = 0.2)
    templates = RD_REG.ProxyTemplate[]
    for typ in (0, 1), parity in (0, 1), mirror in (0, 1)
        matrix = typ == 0 ? M0 : M1
        push!(templates, RD_REG.ProxyTemplate(
            typ, parity, mirror, vec(permutedims(matrix))))
    end
    source = RD_REG.ProxySource("synthetic", "real_diagnostic_test", "", "", 0.0, 1.0, true)
    return RD_REG.ProxyEntry(source, templates)
end

function rd_config()
    symmetric(limit, step) = RD_DC.DiagnosticRange(-limit, limit, step)
    positive(limit, step) = RD_DC.DiagnosticRange(0.0, limit, step)
    return RD_DC.DiagnosticConfig(
        RD_DC.DiagnosticParamBounds(symmetric(0.02, 0.02), symmetric(0.01, 0.01)),
        RD_DC.DiagnosticParamBounds(symmetric(0.02, 0.02), symmetric(0.01, 0.01)),
        RD_DC.DiagnosticParamBounds(symmetric(0.5, 0.5), symmetric(0.25, 0.25)),
        RD_DC.DiagnosticParamBounds(positive(0.01, 0.01), positive(0.01, 0.01)),
        1e-6, 1e-9)
end

function rd_case(set; contrast_amp=0.3)
    xs = collect(-0.7:0.1:0.7)
    ys = collect(-0.4:0.1:0.4)
    centers = [(-0.25, 0.0), (0.25, 0.0)]
    backbone = [sum(exp(-((x-cx)^2 + y^2) / (2 * 0.25^2)) for (cx, _) in centers)
                for y in ys, x in xs]
    return synthesize_whole_roi_observation(set, [1, 0], backbone, centers, xs, ys;
        half_nm=0.2, step_nm=0.2, background=0.2, backbone_amp=0.6,
        common_amp=1.0, contrast_amp=contrast_amp)
end

@testset "frozen perturbations and shifted negative control" begin
    set = derive_common_contrast(rd_entry())
    bundle = run_frozen_real_controls(rd_case(set), set;
        config=rd_config(), half_nm=0.2, step_nm=0.2, control_shift_nm=0.2)
    @test length(bundle.perturbations) == 4
    @test Set(item.label for item in bundle.perturbations) ==
          Set(("shift_t_nm", "shift_u_nm", "rotation_deg", "blur_sigma_nm"))
    @test all(item -> item.registration_tsv == bundle.registration_tsv,
              bundle.controls)
    @test only(bundle.controls).label == "shifted_contrast"
    @test bundle.negative_control_beaten
    @test !bundle.final_abstain
    @test isempty(bundle.final_abstention_reasons)
end

@testset "registration search wall forces explicit abstention" begin
    set = derive_common_contrast(rd_entry())
    case = synthesize_whole_roi_observation(set, [1, 0], rd_case(set).backbone,
        rd_case(set).centers, rd_case(set).xs, rd_case(set).ys;
        half_nm=0.2, step_nm=0.2, shift_t_nm=0.04,
        background=0.2, backbone_amp=0.6, common_amp=1.0, contrast_amp=0.2)
    bundle = run_frozen_real_controls(case, set; config=rd_config(),
        half_nm=0.2, step_nm=0.2, control_shift_nm=0.2)
    @test :shift_t_nm in bundle.nominal.boundary_parameters
    @test bundle.final_abstain
    @test "registration_boundary" in bundle.final_abstention_reasons
end

@testset "null evidence remains an explicit abstention" begin
    set = derive_common_contrast(rd_entry(identical=true))
    bundle = run_frozen_real_controls(rd_case(set; contrast_amp=0.0), set;
        config=rd_config(), half_nm=0.2, step_nm=0.2, control_shift_nm=0.2)
    @test bundle.final_abstain
    @test "zero_contrast_gain" in bundle.final_abstention_reasons
    @test "negative_control_not_beaten" in bundle.final_abstention_reasons
end

@testset "required paths fail before output creation" begin
    mktempdir() do dir
        out = joinpath(dir, "summary.tsv")
        @test_throws ArgumentError validate_real_diagnostic_paths(
            joinpath(dir, "missing-summary.tsv"),
            joinpath(dir, "missing-molds.tsv"),
            joinpath(dir, "missing-config.toml"), dir)
        @test !ispath(out)
    end
end

@testset "diagnostic outputs replace symlinks without following their targets" begin
    mktempdir() do dir
        target_summary = joinpath(dir, "summary-target.tsv")
        target_controls = joinpath(dir, "controls-target.tsv")
        target_plot = joinpath(dir, "plot-target.png")
        write(target_summary, "protected summary\n")
        write(target_controls, "protected controls\n")
        write(target_plot, "protected plot\n")

        summary = joinpath(dir, "diagnostic_summary.tsv")
        controls = joinpath(dir, "diagnostic_controls.tsv")
        plot = joinpath(dir, "case_fused_fitmask.png")
        symlink(target_summary, summary)
        symlink(target_controls, controls)
        symlink(target_plot, plot)

        destinations = [summary, controls, plot]
        staged_seen = String[]
        _with_staged_outputs(destinations) do staged
            append!(staged_seen, values(staged))
            write(staged[summary], "new summary\n")
            write(staged[controls], "new controls\n")
            write(staged[plot], "new plot\n")
        end

        @test all(path -> dirname(path) == dir, staged_seen)
        @test length(unique(staged_seen)) == length(destinations)
        @test all(path -> !islink(path), destinations)
        @test read(summary, String) == "new summary\n"
        @test read(controls, String) == "new controls\n"
        @test read(plot, String) == "new plot\n"
        @test read(target_summary, String) == "protected summary\n"
        @test read(target_controls, String) == "protected controls\n"
        @test read(target_plot, String) == "protected plot\n"
    end
end

@testset "diagnostic failure before commit preserves every destination" begin
    mktempdir() do dir
        summary = joinpath(dir, "diagnostic_summary.tsv")
        controls = joinpath(dir, "diagnostic_controls.tsv")
        plot = joinpath(dir, "case_fused_fitmask.png")
        write(summary, "old summary\n")
        write(controls, "old controls\n")
        write(plot, "old plot\n")
        before = sort(readdir(dir))

        @test_throws ErrorException _with_staged_outputs([summary, controls, plot]) do staged
            write(staged[summary], "partial summary\n")
            write(staged[plot], "partial plot\n")
            error("injected failure before commit")
        end

        @test read(summary, String) == "old summary\n"
        @test read(controls, String) == "old controls\n"
        @test read(plot, String) == "old plot\n"
        @test sort(readdir(dir)) == before
    end
end

@testset "diagnostic install failure rolls back replaced files and symlinks" begin
    mktempdir() do dir
        summary = joinpath(dir, "diagnostic_summary.tsv")
        controls = joinpath(dir, "diagnostic_controls.tsv")
        plot = joinpath(dir, "case_fused_fitmask.png")
        victim = joinpath(dir, "victim.tsv")
        write(summary, "old summary\n")
        write(controls, "old controls\n")
        write(victim, "protected victim\n")
        symlink(victim, plot)

        @test_throws ErrorException _with_staged_outputs(
                [plot, controls, summary]; gate=summary,
                failpoint=:after_install_1) do staged
            write(staged[summary], "new summary\n")
            write(staged[controls], "new controls\n")
            write(staged[plot], "new plot\n")
        end

        @test read(summary, String) == "old summary\n"
        @test read(controls, String) == "old controls\n"
        @test islink(plot)
        @test readlink(plot) == victim
        @test read(victim, String) == "protected victim\n"
        @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(dir)))
    end
end

@testset "diagnostic interrupted transactions recover deterministically on rerun" begin
    for (failpoint, visible_generation) in (
            (:after_marker, :old),
            (:after_backup_1, :old),
            (:after_install_1, :old),
            (:after_gate, :old),
            (:after_commit, :new))
        mktempdir() do dir
            summary = joinpath(dir, "diagnostic_summary.tsv")
            controls = joinpath(dir, "diagnostic_controls.tsv")
            plot = joinpath(dir, "case_fused_fitmask.png")
            destinations = [plot, controls, summary]
            for (path, text) in zip(destinations, ("old plot\n", "old controls\n", "old summary\n"))
                write(path, text)
            end

            @test_throws SimulatedTransactionInterruption _with_staged_outputs(
                    destinations; gate=summary, interruptpoint=failpoint) do staged
                for (path, text) in zip(destinations, ("new plot\n", "new controls\n", "new summary\n"))
                    write(staged[path], text)
                end
            end
            if failpoint in (:after_backup_1, :after_install_1)
                @test !ispath(summary)
            end

            observed = Ref{Symbol}(:unset)
            _with_staged_outputs(destinations; gate=summary) do staged
                observed[] = startswith(read(summary, String), "old") ? :old : :new
                for (path, text) in zip(destinations, ("rerun plot\n", "rerun controls\n", "rerun summary\n"))
                    write(staged[path], text)
                end
            end
            @test observed[] == visible_generation
            @test read.(destinations, String) ==
                ["rerun plot\n", "rerun controls\n", "rerun summary\n"]
            @test isempty(filter(name -> occursin("stmfit-txn", name), readdir(dir)))
        end
    end
end

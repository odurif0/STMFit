using Test
VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
const FUSION_SCRIPT = joinpath(dirname(@__DIR__), "hpc", "compare_fusion.sbatch")

@testset "Paired fusion job regenerates both arms at identical cached N" begin
    @test success(`bash -n $FUSION_SCRIPT`)
    mktempdir() do dir
        control, symmetric = joinpath.(dir, ("control source", "symmetric source"))
        for project in (control, symmetric), file in ("Project.toml", "Manifest.toml",
                "config/chitosan.toml", "config/unit_assignment_reconstructed.toml")
            path = joinpath(project, file); mkpath(dirname(path)); write(path, "same\n")
        end
        input = joinpath(dir, "inputs"); mkpath(joinpath(input, "full146_raw"))
        write(joinpath(input, "templates_cc.tsv"), "synthetic\n")
        selected = joinpath(input, "selected.tsv")
        write(selected, "filepath\tN_selected\nsynthetic.sxm\t3\n")
        fake = joinpath(dir, "fake julia")
        write(fake, raw"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\0' "$PWD" "$OPENBLAS_NUM_THREADS" "$JULIA_NUM_THREADS" "$@" >> "${FUSION_TEST_CAPTURE:?}"
if [[ "${FUSION_TEST_FAIL:-}" == "$PWD" ]]; then exit 37; fi
""")
        chmod(fake, 0o755)
        function run_pair(name; args=["--dry-run"], extra=[])
            output = joinpath(dir, name)
            capture = joinpath(dir, name * ".calls")
            env = ["STMFIT_CONTROL_PROJECT"=>control, "STMFIT_SYMMETRIC_PROJECT"=>symmetric,
                "STMFIT_INPUT_DIR"=>input, "STMFIT_SELECTED_SUMMARY"=>selected,
                "STMFIT_OUTDIR"=>output, "JULIA_BIN"=>fake, "FUSION_TEST_CAPTURE"=>capture,
                "SLURM_JOB_ID"=>nothing, "SLURM_CPUS_PER_TASK"=>nothing]
            io = IOBuffer()
            p = run(pipeline(ignorestatus(addenv(`bash $FUSION_SCRIPT $args`, env..., extra...)); stdout=io, stderr=io))
            return p.exitcode, String(take!(io)), output,
                isfile(capture) ? split(read(capture, String), '\0'; keepempty=false) : String[]
        end
        code, log, output, argv = run_pair("dry")
        @test code == 0
        @test !ispath(output)
        expected = String[]
        for (arm, project) in (("control", control), ("symmetric", symmetric))
            append!(expected, [project, "1", "4", "--startup-file=no", "-t", "4", "--project=.",
                "test/run_reconstructed_chitosan.jl", "--data-dir", joinpath(input, "full146_raw"),
                "--count-config", "config/chitosan.toml", "--config", "config/unit_assignment_reconstructed.toml",
                "--selected-summary", selected, "--templates", joinpath(input, "templates_cc.tsv"),
                "--outdir", joinpath(output, arm), "--dry-run"])
        end
        @test argv == expected # No cached geometry, patches or benchmark arguments.
        code, log, output, argv = run_pair("no allocation"; args=String[])
        @test code != 0 && occursin("Slurm", log)
        @test isempty(argv) && !ispath(output)
        allocation = ["SLURM_JOB_ID"=>"synthetic", "SLURM_CPUS_PER_TASK"=>"4"]
        code, log, output, argv = run_pair("executed"; args=String[], extra=allocation)
        @test code == 0
        @test isfile(joinpath(output, "control.log")) && isfile(joinpath(output, "symmetric.log"))
        @test count(==("--selected-summary"), argv) == 2
        code, log, output, argv = run_pair("control failed"; args=String[],
            extra=vcat(allocation, ["FUSION_TEST_FAIL"=>control]))
        @test code != 0 && occursin("control failed (exit 37)", log)
        @test occursin("symmetric completed", log)
        @test isfile(joinpath(output, "control.log"))
        code, _, _, argv = run_pair("duplicate"; extra=["STMFIT_OUTDIR"=>output])
        @test code != 0 && isempty(argv)
        code, _, _, argv = run_pair("same source"; extra=["STMFIT_SYMMETRIC_PROJECT"=>control])
        @test code != 0 && isempty(argv)
        write(joinpath(symmetric,"config/chitosan.toml"), "different\n")
        code, _, output, argv = run_pair("mismatched config")
        @test code != 0 && isempty(argv) && !ispath(output)
    end
end

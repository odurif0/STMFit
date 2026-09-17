# julia --project=. test/test_reconstructed_full146_sbatch.jl
# Synthetic metadata and command capture only: no SXM parsing, fits or classifiers.
using Test
VERSION.major == 1 && VERSION.minor == 13 || error("Run this test with Julia 1.13")
const TEST_ROOT = dirname(@__DIR__)
const SCRIPT = joinpath(TEST_ROOT, "hpc", "reconstructed_full146.sbatch")
const JULIA = joinpath(Sys.BINDIR, "julia")

function fixture(dir)
    input = joinpath(dir, "input with spaces")
    raw = joinpath(input, "full146_raw"); mkpath(raw)
    summary = "filepath\tN_selected\tstatus\trefined_policy\n"
    base = "file\tN\tlobe\tskew_ratio\n"
    for i in 1:146
        file = "synthetic_$(lpad(i, 3, '0')).sxm"
        n = 2 + mod(i, 3)
        write(joinpath(raw, file), "Not an SXM; metadata checks must not read this.\n")
        summary *= "/synthetic/$file\t$n\tok\tsupport_midpoint_hybrid\n"
        base *= join("$file\t$n\t$lobe\t1.0\n" for lobe in 1:n)
    end
    paths = (summary=joinpath(input, "full146_selected_summary.tsv"),
             base=joinpath(input, "base_geometry_full146.tsv"),
             templates=joinpath(input, "templates_cc.tsv"))
    write(paths.summary, summary); write(paths.base, base)
    write(paths.templates, "Synthetic placeholder; dry-run checks paths only.\n")
    return (; input, raw, paths..., out=joinpath(dir, "not created", "full146 outputs"))
end
snapshot(input) = Dict(joinpath(dir, file) => read(joinpath(dir, file))
                       for (dir, _, files) in walkdir(input) for file in files)

function shell_check(f; args=["--dry-run"], out=f.out, julia=JULIA, env=[])
    cmd = Cmd(`bash $SCRIPT $args`; dir=dirname(f.input))
    cmd = addenv(cmd, "STMFIT_PROJECT_DIR" => TEST_ROOT, "JULIA_BIN" => julia,
        "STMFIT_INPUT_DIR" => f.input, "STMFIT_OUTDIR" => out,
        "STMFIT_SELECTED_SUMMARY" => nothing, "SLURM_CPUS_PER_TASK" => "4",
        "JULIA_NUM_THREADS" => "64", "OPENBLAS_NUM_THREADS" => "64",
        "OMP_NUM_THREADS" => "64", "GKSwstype" => "unset", env...)
    io = IOBuffer()
    process = run(pipeline(ignorestatus(cmd); stdout=io, stderr=io))
    return process.exitcode, String(take!(io))
end

@testset "Full146 single-job shell contract" begin
    @test success(`bash -n $SCRIPT`)
    source = read(SCRIPT, String)
    for directive in ("--nodes=1", "--ntasks=1", "--cpus-per-task=4", "--mem=16000MB", "--time=24:00:00")
        @test occursin("#SBATCH $directive\n", source)
    end
    @test !occursin(r"(?m)^#SBATCH\s+--(?:partition|array)", source)
    mktempdir() do dir
        f = fixture(dir); before = snapshot(f.input)
        fake = joinpath(dir, "fake julia"); capture = joinpath(dir, "calls.bin")
        write(fake, raw"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\0' "$PWD" "$JULIA_NUM_THREADS" "$OPENBLAS_NUM_THREADS" "$OMP_NUM_THREADS" "$GKSwstype" "$@" >> "${FULL146_TEST_CAPTURE:?}"
exit "${FULL146_TEST_EXIT:-0}"
""")
        chmod(fake, 0o755)
        capture_env = ["FULL146_TEST_CAPTURE" => capture, "FULL146_TEST_EXIT" => "0"]
        expected = ["-t", "4", "--project=.", "test/run_reconstructed_chitosan.jl",
            "--data-dir", f.raw, "--count-config", "config/chitosan.toml",
            "--config", "config/unit_assignment_reconstructed.toml",
            "--features", f.base, "--templates", f.templates, "--outdir", f.out]
        for (selected, summary) in ((false, nothing), (false, ""), (true, f.summary)), args in (String[], ["--dry-run"])
            rm(capture; force=true)
            env = vcat(capture_env, ["STMFIT_SELECTED_SUMMARY" => summary])
            code, _ = shell_check(f; julia=fake, args, env)
            @test code == 0
            fields = split(chomp(read(capture, String)), '\0'; keepempty=false)
            @test fields[1:5] == [TEST_ROOT, "4", "1", "1", "100"]
            selected_args = selected ? ["--selected-summary", f.summary] : String[]
            # Exact argv proves one native call, no count/base replay, no cached
            # all-one skew masquerading as split, and no comparison/grading CLI.
            @test fields[6:end] == vcat(expected, selected_args, args)
            @test !ispath(dirname(f.out))
        end
        # Unknown options never reach Julia, including attempts to inject a
        # split cache, labels, a different feature table, or another output path.
        for args in (["--split-features", f.base], ["--features", f.base],
                     ["--reference", "unused.tsv"], ["--expected-N", "6"],
                     ["--outdir", f.input], ["--dry-run", "--dry-run"], [""], ["--help"])
            rm(capture; force=true)
            code, log = shell_check(f; julia=fake, args, env=capture_env)
            @test code == 2
            @test occursin("Usage:", log)
            @test !ispath(capture)
        end
        for key in ("STMFIT_PROJECT_DIR", "JULIA_BIN", "STMFIT_INPUT_DIR", "STMFIT_OUTDIR")
            code, log = shell_check(f; julia=fake, env=vcat(capture_env, [key => nothing]))
            @test code != 0
            @test occursin(key, log)
            @test !ispath(capture)
        end
        code, _ = shell_check(f; julia=fake, env=["FULL146_TEST_CAPTURE" => capture, "FULL146_TEST_EXIT" => "37"])
        @test code == 37
        @test snapshot(f.input) == before
        @test !ispath(dirname(f.out))
    end
end

@testset "Actual native dry-run and input/output boundaries" begin
    mktempdir() do dir
        f = fixture(dir); before = snapshot(f.input)
        configs = [joinpath(TEST_ROOT, "config", file) for file in ("chitosan.toml", "unit_assignment_reconstructed.toml")]
        config_before = read.(configs)
        for env in ([], ["STMFIT_SELECTED_SUMMARY" => f.summary])
            code, log = shell_check(f; env)
            @test code == 0
            @test occursin("method=cc_soft_reconstructed_v1; files=146", log)
            @test occursin("No computation or output writing in dry-run. Output: " * f.out, log)
            @test !ispath(dirname(f.out))
        end
        # A relative STMFIT_OUTDIR is relative to the synced project, and still
        # names the direct output directory (there is no implicit child folder).
        code, log = shell_check(f; out=relpath(f.out, TEST_ROOT))
        @test code == 0
        @test occursin("Output: " * f.out, log)
        @test !ispath(dirname(f.out))
        for (path, reason, env) in ((f.raw, "Data directory not found", []),
                (f.base, "--features input does not exist", []),
                (f.templates, "--templates input does not exist", []),
                (f.summary, "--selected-summary input does not exist", ["STMFIT_SELECTED_SUMMARY" => f.summary]))
            held = path * ".held"; mv(path, held)
            try
                code, log = shell_check(f; env)
                @test code != 0
                @test occursin(reason, log)
                @test !ispath(dirname(f.out))
            finally
                mv(held, path)
            end
        end
        # An unrelated/incompatible saved summary is rejected, not allowed to
        # silently redefine the feature-derived N or trigger a count/base refit.
        original = read(f.summary, String)
        try
            write(f.summary, replace(original, "\t3\tok\t" => "\t4\tok\t"; count=1))
            code, log = shell_check(f)
            @test code == 0  # No auto-discovery of a summary beside the cache.
            @test occursin("method=cc_soft_reconstructed_v1; files=146", log)
            @test !ispath(dirname(f.out))
            code, log = shell_check(f; env=["STMFIT_SELECTED_SUMMARY" => f.summary])
            @test code != 0
            @test occursin("Feature/count coverage mismatch", log)
            @test !ispath(dirname(f.out))
            write(f.summary, replace(original, "refined_policy" => "truth"))
            code, log = shell_check(f; env=["STMFIT_SELECTED_SUMMARY" => f.summary])
            @test code != 0
            @test occursin("Forbidden benchmark/control column", log)
            @test !ispath(dirname(f.out))
        finally
            write(f.summary, original)
        end
        # Native parse_options rejects directories, files and dangling links in
        # both modes, before a fit or an output write can occur.
        dangling = joinpath(dir, "dangling output")
        symlink(joinpath(dir, "absent target"), dangling)
        for out in (f.input, f.base, dangling), args in (["--dry-run"], String[])
            code, log = shell_check(f; out, args)
            @test code != 0
            @test occursin("Output directory already exists", log)
            @test !ispath(dirname(f.out))
        end
        @test snapshot(f.input) == before
        @test read.(configs) == config_before
    end
end

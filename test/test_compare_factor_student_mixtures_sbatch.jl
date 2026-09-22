using Test
VERSION.major == 1 && VERSION.minor == 13 || error("Julia 1.13 required")
const SCRIPT = joinpath(dirname(@__DIR__), "hpc", "compare_factor_student_mixtures.sbatch")

@testset "Two isolated candidates and their unchanged support control" begin
    @test success(`bash -n $SCRIPT`)
    source = read(SCRIPT, String)
    @test occursin("#SBATCH --time=01:00:00", source)
    @test occursin("#SBATCH --cpus-per-task=4", source)
    @test occursin("#SBATCH --mem=16000MB", source)
    @test occursin("#SBATCH --no-requeue", source)
    mktempdir() do dir
        project, cache, input = joinpath.(dir, ("source checkout", "geometry cache", "inputs"))
        mkpath(project); mkpath(cache); mkpath(joinpath(input, "full146_raw"))
        for file in ("features.tsv", "features_split.tsv")
            write(joinpath(cache,file),"synthetic\n")
        end
        write(joinpath(input,"templates_cc.tsv"),"synthetic\n")
        fake = joinpath(dir,"fake julia")
        write(fake, raw"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\0' "$PWD" "$OPENBLAS_NUM_THREADS" "$JULIA_NUM_THREADS" "$@" >> "${SIGNAL_TEST_CAPTURE:?}"
for arg in "$@"; do
    if [[ "${SIGNAL_TEST_FAIL:-}" == "$arg" ]]; then exit 37; fi
done
""")
        chmod(fake,0o755)
        function run_comparison(name; args=["--dry-run"],extra=[])
            output,capture = joinpath(dir,name),joinpath(dir,name*".calls")
            env = ["STMFIT_PROJECT_DIR"=>project,"STMFIT_CACHE_DIR"=>cache,"STMFIT_INPUT_DIR"=>input,
                "STMFIT_OUTDIR"=>output,"JULIA_BIN"=>fake,"SIGNAL_TEST_CAPTURE"=>capture,
                "SLURM_JOB_ID"=>nothing,"SLURM_CPUS_PER_TASK"=>nothing]
            io=IOBuffer()
            p=run(pipeline(ignorestatus(addenv(`bash $SCRIPT $args`,env...,extra...));stdout=io,stderr=io))
            return p.exitcode,String(take!(io)),output,
                isfile(capture) ? split(read(capture,String),'\0';keepempty=false) : String[]
        end
        code,_,output,argv=run_comparison("dry")
        @test code==0
        @test !ispath(output)
        expected=String[]
        for (arm,cfg) in (("control","unit_assignment_patch_support.toml"),("factor_analyzer","unit_assignment_factor_analyzer.toml"),("student_t","unit_assignment_student_t.toml"))
            append!(expected,[project,"1","4","--startup-file=no","-t","4","--project=.",
                "test/run_reconstructed_chitosan.jl","--data-dir",joinpath(input,"full146_raw"),
                "--count-config","config/chitosan.toml","--config","config/"*cfg,
                "--templates",joinpath(input,"templates_cc.tsv"),"--features",joinpath(cache,"features.tsv"),
                "--split-features",joinpath(cache,"features_split.tsv"),"--outdir",joinpath(output,arm),"--dry-run"])
        end
        @test argv==expected # No patch-cache, label, benchmark or fitting-budget argument.
        code,log,output,argv=run_comparison("no allocation";args=String[])
        @test code!=0 && occursin("Slurm",log)
        @test isempty(argv) && !ispath(output)
        allocation=["SLURM_JOB_ID"=>"synthetic","SLURM_CPUS_PER_TASK"=>"4"]
        code,_,output,argv=run_comparison("run";args=String[],extra=allocation)
        @test code==0
        @test all(isfile(joinpath(output,arm*".log")) for arm in ("control","factor_analyzer","student_t"))
        @test count(==("--features"),argv)==3
        code,log,output,argv=run_comparison("failed";args=String[],extra=vcat(allocation,
            ["SIGNAL_TEST_FAIL"=>"config/unit_assignment_patch_support.toml"]))
        @test code!=0 && occursin("control failed (exit 37)",log)
        @test occursin("factor_analyzer completed",log) && occursin("student_t completed",log)
        @test isfile(joinpath(output,"control.log"))
        code,_,_,argv=run_comparison("exists";extra=["STMFIT_OUTDIR"=>output])
        @test code!=0 && isempty(argv)
        code,_,output,argv=run_comparison("relative";extra=["STMFIT_CACHE_DIR"=>"relative/cache"])
        @test code!=0 && isempty(argv) && !ispath(output)
        code,_,output,argv=run_comparison("small allocation";args=String[],
            extra=["SLURM_JOB_ID"=>"synthetic","SLURM_CPUS_PER_TASK"=>"2"])
        @test code!=0 && isempty(argv) && !ispath(output)
        rm(joinpath(cache,"features_split.tsv"))
        code,_,output,argv=run_comparison("missing geometry")
        @test code!=0 && isempty(argv) && !ispath(output)
    end
end
